# streamlit_app.py

import streamlit as st
import pandas as pd
import requests
import base64
from io import BytesIO, StringIO

st.set_page_config(page_title="OneCLick Keyword Retriever", layout="centered")
st.title("🚀 OneClick")

uploaded = st.file_uploader("📁 Upload Excel (.xlsx) with 'DI' (DOI)", type=["xlsx"])

# --- Fetch Keywords from OpenAlex ---
def fetch_openalex_keywords(doi):
    try:
        r = requests.get(f"https://api.openalex.org/works/doi:{doi}")
        if r.status_code != 200:
            return ""
        data = r.json()
        keywords = data.get("keywords", [])
        return '; '.join(sorted(set(k['display_name'].strip() for k in keywords)))
    except Exception:
        return ""

# --- File Download Helpers ---
def generate_excel_download(df):
    output = BytesIO()
    df.to_excel(output, index=False, engine='openpyxl')
    b64 = base64.b64encode(output.getvalue()).decode()
    href = f'<a href="data:application/octet-stream;base64,{b64}" download="openalex_keywords.xlsx">📥 Download .xlsx File</a>'
    return href

def generate_csv_download(df):
    output = StringIO()
    df.to_csv(output, index=False, encoding='utf-8-sig')
    b64 = base64.b64encode(output.getvalue().encode()).decode()
    href = f'<a href="data:file/csv;base64,{b64}" download="openalex_keywords.csv">📥 Download .csv File</a>'
    return href

# --- Main App Flow ---
if uploaded:
    df = pd.read_excel(uploaded, engine='openpyxl')

    if 'DI' not in df.columns:
        st.error("❌ 'DI' column (DOI) is required.")
    else:
        df['OpenAlex_KW'] = ""

        if st.button("🚀 Retrieve Keywords from OpenAlex"):
            with st.spinner("Fetching from OpenAlex..."):
                for idx, row in df.iterrows():
                    if pd.notna(row['DI']):
                        df.at[idx, 'OpenAlex_KW'] = fetch_openalex_keywords(row['DI'])

            st.success("✅ Keywords Retrieved!")
            st.dataframe(df[['DI', 'OpenAlex_KW']].head())

            st.markdown("### 📄 Download Enriched Dataset")
            st.markdown("The output retains your original structure with an additional `OpenAlex_KW` column for enriched keywords. Choose your preferred format:")
            st.markdown(generate_excel_download(df), unsafe_allow_html=True)
            st.markdown(generate_csv_download(df), unsafe_allow_html=True)
