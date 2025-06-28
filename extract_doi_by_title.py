import streamlit as st
import pandas as pd
import requests
from io import BytesIO
import base64
from rapidfuzz.fuzz import ratio as fuzz_ratio

st.set_page_config(page_title="DOI Extractor from Title", layout="centered")
st.title("🔍 Extract Missing DOIs using Title via Crossref")

uploaded_file = st.file_uploader("📄 Upload Excel File (must contain 'Author', 'Title', 'Year', 'DI', and 'Author Keywords')", type=["xlsx"])

CONFIDENCE_THRESHOLD = 95

# --- Function to Extract DOI with strict matching ---
def get_best_doi(title, author, year):
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
                    # Mandatory author match
                    api_authors = item.get("author", [])
                    api_author_names = [a.get("family", "").lower() for a in api_authors]
                    if not any(a in author.lower() for a in api_author_names):
                        continue

                    # Mandatory year match
                    api_year = item.get("issued", {}).get("date-parts", [[None]])[0][0]
                    if not api_year or abs(api_year - int(year)) > 1:
                        continue

                    return item.get("DOI", None)
    except Exception as e:
        print(f"Error fetching DOI for title: {title} — {e}")
    return None

# --- App Execution ---
if uploaded_file:
    df = pd.read_excel(uploaded_file, engine='openpyxl')

    required_columns = {'Author', 'Title', 'Year', 'DI', 'Author Keywords'}
    if not required_columns.issubset(df.columns):
        st.error("❌ The uploaded file must contain columns: 'Author', 'Title', 'Year', 'DI', and 'Author Keywords'.")
    else:
        if st.button("🚀 Extract DOIs"):
            with st.spinner("Looking up missing DOIs using Crossref with strict author/year matching..."):
                original_count = len(df)
                df['DI'] = df.apply(
                    lambda row: get_best_doi(
                        row['Title'],
                        author=row['Author'],
                        year=row['Year']
                    ) if pd.isna(row['DI']) and pd.notna(row['Title']) and pd.notna(row['Author']) and pd.notna(row['Year']) else row['DI'],
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
