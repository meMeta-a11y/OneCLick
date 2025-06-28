# streamlit_doi_extractor.py

import streamlit as st
import pandas as pd
import requests
from io import BytesIO
import base64
from rapidfuzz.fuzz import ratio as fuzz_ratio

st.set_page_config(page_title="DOI Extractor from Title", layout="centered")
st.title("🔍 Extract Missing DOIs using Title via Crossref")

uploaded_file = st.file_uploader("📄 Upload Excel File (must contain 'Title' and 'DI' columns)", type=["xlsx"])

CONFIDENCE_THRESHOLD = 90

# --- Function to Extract DOI ---
def get_best_doi(title, author=None, year=None):
    try:
        url = f"https://api.crossref.org/works?query.title={requests.utils.quote(title)}&rows=5"
        r = requests.get(url, timeout=10)
        if r.status_code == 200:
            data = r.json()
            items = data.get("message", {}).get("items", [])
            for item in items:
                api_title = item.get("title", [""])[0]
                similarity = fuzz_ratio(title.lower(), api_title.lower())
                if similarity >= CONFIDENCE_THRESHOLD:
                    # --- Optional: Check author match ---
                    if author:
                        api_authors = item.get("author", [])
                        api_author_names = [a.get("family", "").lower() for a in api_authors]
                        if not any(a in author.lower() for a in api_author_names):
                            continue
                    # --- Optional: Check year match ---
                    if year:
                        api_year = item.get("issued", {}).get("date-parts", [[None]])[0][0]
                        if api_year and abs(api_year - int(year)) > 1:
                            continue
                    return item.get("DOI", None)
    except Exception as e:
        print(f"Error fetching DOI for title: {title} — {e}")
    return None

# --- App Execution ---
if uploaded_file:
    df = pd.read_excel(uploaded_file, engine='openpyxl')

    if 'Title' not in df.columns:
        st.error("❌ The uploaded file must have a 'Title' column.")
    else:
        if 'DI' not in df.columns:
            df['DI'] = pd.NA

        if st.button("🚀 Extract DOIs"):
            with st.spinner("Looking up missing DOIs using Crossref with fuzzy matching and author/year cross-check..."):
                original_count = len(df)
                df['DI'] = df.apply(
                    lambda row: get_best_doi(
                        row['Title'],
                        author=row['Author'] if 'Author' in df.columns else None,
                        year=row['Year'] if 'Year' in df.columns else None
                    ) if pd.isna(row['DI']) and pd.notna(row['Title']) else row['DI'],
                    axis=1
                )

            st.success("✅ DOI extraction completed!")
            st.write(f"Rows processed: {original_count}")
            st.dataframe(df.head())

            output = BytesIO()
            df.to_excel(output, index=False, engine='openpyxl')
            b64 = base64.b64encode(output.getvalue()).decode()
            st.markdown(f'<a href="data:application/octet-stream;base64,{b64}" download="with_doi.xlsx">📥 Download Updated File</a>', unsafe_allow_html=True)
else:
    st.info("👆 Please upload an Excel file to begin.")
