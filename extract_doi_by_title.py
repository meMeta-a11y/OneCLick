import streamlit as st
import pandas as pd
import requests
from io import BytesIO
import base64
from rapidfuzz.fuzz import ratio as fuzz_ratio

# --- Page setup ---
st.set_page_config(page_title="DOI Extractor", layout="centered")
st.title("🔍 Extract Missing DOIs using Title + Author/Year")

# --- Upload file ---
uploaded_file = st.file_uploader("📄 Upload Excel file with columns: Author, Title, Year, DI, Author Keywords", type=["xlsx"])

CONFIDENCE_THRESHOLD = 95

# --- DOI Lookup Function ---
def get_best_doi(title, author, year):
    try:
        url = f"https://api.crossref.org/works?query.title={requests.utils.quote(title)}&rows=5"
        r = requests.get(url, timeout=10)
        if r.status_code == 200:
            for item in r.json().get("message", {}).get("items", []):
                api_title = item.get("title", [""])[0]
                similarity = fuzz_ratio(title.lower(), api_title.lower())
                if similarity >= CONFIDENCE_THRESHOLD:
                    api_authors = item.get("author", [])
                    api_author_names = [a.get("family", "").lower() for a in api_authors]
                    if not any(a in author.lower() for a in api_author_names):
                        continue
                    api_year = item.get("issued", {}).get("date-parts", [[None]])[0][0]
                    if not api_year or abs(api_year - int(year)) > 1:
                        continue
                    return item.get("DOI", None)
    except Exception as e:
        print(f"Error for title '{title}': {e}")
    return None

# --- Main App Logic ---
if uploaded_file:
    df = pd.read_excel(uploaded_file, engine='openpyxl')
    df.columns = df.columns.str.strip()  # Normalize column names

    required_cols = {"Author", "Title", "Year", "DI", "Author Keywords"}
    missing = required_cols - set(df.columns)

    if missing:
        st.error(f"❌ Missing required columns: {', '.join(missing)}")
    else:
        if st.button("🚀 Start DOI Extraction"):
            with st.spinner("Looking up missing DOIs using Crossref..."):
                df['DI'] = df.apply(
                    lambda row: get_best_doi(row['Title'], row['Author'], row['Year'])
                    if pd.isna(row['DI']) and pd.notna(row['Title']) and pd.notna(row['Author']) and pd.notna(row['Year'])
                    else row['DI'],
                    axis=1
                )

            st.success("✅ Done! Here's a preview:")
            st.dataframe(df.head())

            output = BytesIO()
            df.to_excel(output, index=False, engine='openpyxl')
            b64 = base64.b64encode(output.getvalue()).decode()
            st.markdown(f'<a href="data:application/octet-stream;base64,{b64}" download="with_doi.xlsx">📥 Download Updated Excel</a>', unsafe_allow_html=True)
else:
    st.info("👆 Upload a valid Excel file to begin.")
