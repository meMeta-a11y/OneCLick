# streamlit_doi_extractor.py

import streamlit as st
import pandas as pd
import requests
from io import BytesIO
import base64

st.set_page_config(page_title="DOI Extractor from Title", layout="centered")
st.title("🔍 Extract Missing DOIs using Title via Crossref")

uploaded_file = st.file_uploader("📄 Upload Excel File (must contain 'Title' and 'DI' columns)", type=["xlsx"])

def get_doi_from_title(title):
    try:
        url = f"https://api.crossref.org/works?query.title={requests.utils.quote(title)}&rows=1"
        r = requests.get(url, timeout=10)
        if r.status_code == 200:
            data = r.json()
            items = data.get("message", {}).get("items", [])
            if items:
                return items[0].get("DOI", None)
    except Exception as e:
        print(f"Error fetching DOI for title: {title} — {e}")
    return None

if uploaded_file:
    df = pd.read_excel(uploaded_file, engine='openpyxl')

    if 'Title' not in df.columns:
        st.error("❌ The uploaded file must have a 'Title' column.")
    else:
        if 'DI' not in df.columns:
            df['DI'] = pd.NA

        # Fill missing DI with DOIs from Crossref
        if st.button("🚀 Extract DOIs"):
            with st.spinner("Looking up missing DOIs using Crossref..."):
                original_count = len(df)
                df['DI'] = df.apply(
                    lambda row: get_doi_from_title(row['Title']) if pd.isna(row['DI']) and pd.notna(row['Title']) else row['DI'],
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
