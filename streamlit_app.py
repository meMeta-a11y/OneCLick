# streamlit_app.py

import streamlit as st
import pandas as pd
import requests
import base64
from io import BytesIO

st.set_page_config(page_title="One-Click Keyword Retriever", layout="centered")
st.title("🚀 One-Click Keyword Retriever")

uploaded_doi_file = st.file_uploader("📁 Step 1: Upload Excel (.xlsx) with 'DI' (DOI)", type=["xlsx"])
uploaded_existing_kw_file = st.file_uploader("📁 Step 2: Upload Dataset with Existing Keywords (Optional)", type=["xlsx"])

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

# --- Main App Flow ---
if uploaded_doi_file:
    df = pd.read_excel(uploaded_doi_file, engine='openpyxl')

    if 'DI' not in df.columns:
        st.error("❌ 'DI' column (DOI) is required in the first file.")
    else:
        df['Keywords'] = df.get('Keywords', "")

        if st.button("🚀 Retrieve and Combine Keywords"):
            with st.spinner("Fetching missing keywords from OpenAlex..."):
                for idx, row in df.iterrows():
                    if (not pd.notna(row['Keywords']) or not row['Keywords'].strip()) and pd.notna(row['DI']):
                        df.at[idx, 'Keywords'] = fetch_openalex_keywords(row['DI'])

            if uploaded_existing_kw_file:
                existing_df = pd.read_excel(uploaded_existing_kw_file, engine='openpyxl')

                if 'DI' not in existing_df.columns or 'Keywords' not in existing_df.columns:
                    st.error("❌ The second file must contain 'DI' and 'Keywords' columns.")
                else:
                    df = pd.merge(df, existing_df[['DI', 'Keywords']], on='DI', how='left', suffixes=('', '_existing'))
                    df['Keywords'] = df.apply(lambda row: row['Keywords_existing'] if pd.notna(row['Keywords_existing']) and row['Keywords_existing'].strip() else row['Keywords'], axis=1)
                    df.drop(columns=['Keywords_existing'], inplace=True)

            st.success("✅ Keywords Retrieved and Combined!")
            st.dataframe(df[['DI', 'Keywords']].head())

            output = BytesIO()
            df.to_excel(output, index=False, engine='openpyxl')
            b64 = base64.b64encode(output.getvalue()).decode()
            st.markdown(f'<a href="data:application/octet-stream;base64,{b64}" download="combined_keywords.xlsx"> Click to Download</a>', unsafe_allow_html=True)
