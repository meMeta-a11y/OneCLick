import streamlit as st
import pandas as pd
import requests
from urllib.parse import quote
from io import BytesIO
import base64

st.set_page_config(page_title="Extract DOI via Crossref", layout="centered")
st.title("🔍 Extract DOI from Crossref using Title")

uploaded_file = st.file_uploader("📄 Upload Excel file (must contain 'TI' for Title, and 'DI' for DOI)", type=["xlsx"])

def get_doi_from_crossref(title):
    try:
        query_url = f"https://api.crossref.org/works?query.title={quote(title)}&rows=5"
        response = requests.get(query_url, timeout=10)
        response.raise_for_status()
        items = response.json().get("message", {}).get("items", [])
        for item in items:
            if item.get("title") and title.strip().lower() == item["title"][0].strip().lower():
                return item.get("DOI", None)
    except:
        return None
    return None

if uploaded_file:
    df = pd.read_excel(uploaded_file, engine='openpyxl')

    if 'TI' not in df.columns:
        st.error("❌ The file must contain a 'TI' column for Title.")
    else:
        if 'DI' not in df.columns:
            df['DI'] = pd.NA  # Add DI column if missing

        df['DI'] = df['DI'].astype(str).str.strip()
        df['DI'] = df['DI'].replace({'': pd.NA, 'nan': pd.NA})

        if st.button("🚀 Start DOI Extraction"):
            with st.spinner("Contacting Crossref..."):
                for idx, row in df.iterrows():
                    if pd.isna(row['DI']):
                        doi = get_doi_from_crossref(row['TI'])
                        df.at[idx, 'DI'] = doi if doi else pd.NA

            st.success("✅ DOI extraction complete.")
            st.dataframe(df.head())

            output = BytesIO()
            df.to_excel(output, index=False, engine='openpyxl')
            b64 = base64.b64encode(output.getvalue()).decode()
            st.markdown(
                f'<a href="data:application/octet-stream;base64,{b64}" download="updated_with_doi.xlsx">📥 Download Updated File</a>',
                unsafe_allow_html=True
            )
