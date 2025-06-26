# streamlit_app.py

import streamlit as st
import pandas as pd
import requests
import base64
from io import BytesIO

st.set_page_config(page_title="One-Click Keyword Processor", layout="centered")
st.title("🧠 One-Click Keyword Processor")

step = st.radio("🔍 Select Operation Step:", [
    "🔹 Step 1: Fetch Keywords from OpenAlex",
    "🔹 Step 2: Merge Keywords into Existing Dataset"
])

# --- Helper Function ---
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

# --- Step 1: Fetch Missing Keywords ---
if step == "🔹 Step 1: Fetch Keywords from OpenAlex":
    st.subheader("Step 1️⃣: Fetch Missing Keywords")
    uploaded_doi_file = st.file_uploader("📄 Upload DOI-only Excel file (must contain 'DI')", type=["xlsx"], key="fetch")

    if uploaded_doi_file:
        df = pd.read_excel(uploaded_doi_file, engine='openpyxl')

        if 'DI' not in df.columns:
            st.error("❌ Uploaded file must contain a 'DI' column.")
        else:
            if st.button("🚀 Fetch Keywords from OpenAlex"):
                with st.spinner("Fetching keywords from OpenAlex..."):
                    df['OpenAlex_KW'] = df['DI'].apply(lambda x: fetch_openalex_keywords(x) if pd.notna(x) else "")

                st.success("✅ Keywords successfully fetched!")
                st.dataframe(df.head())

                output = BytesIO()
                df.to_excel(output, index=False, engine='openpyxl')
                b64 = base64.b64encode(output.getvalue()).decode()
                st.markdown(f'<a href="data:application/octet-stream;base64,{b64}" download="openalex_keywords.xlsx">📥 Download OpenAlex Keywords</a>', unsafe_allow_html=True)

# --- Step 2: Merge with Existing Dataset ---
elif step == "🔹 Step 2: Merge Keywords into Existing Dataset":
    st.subheader("Step 2️⃣: Merge OpenAlex Keywords into Existing Dataset")
    uploaded_keywords_file = st.file_uploader("📥 Upload 'openalex_keywords.xlsx' (must contain 'DI', 'OpenAlex_KW')", type=["xlsx"], key="merge_kw")
    uploaded_data_file = st.file_uploader("📥 Upload original dataset (must contain 'DI', 'Author Keywords')", type=["xlsx"], key="merge_data")

    if uploaded_keywords_file and uploaded_data_file:
        keywords_df = pd.read_excel(uploaded_keywords_file, engine='openpyxl')
        data_df = pd.read_excel(uploaded_data_file, engine='openpyxl')

        if 'DI' not in keywords_df.columns or 'OpenAlex_KW' not in keywords_df.columns:
            st.error("❌ 'openalex_keywords.xlsx' must contain 'DI' and 'OpenAlex_KW' columns.")
        elif 'DI' not in data_df.columns or 'Author Keywords' not in data_df.columns:
            st.error("❌ 'data.xlsx' must contain 'DI' and 'Author Keywords' columns.")
        else:
            merged_df = pd.merge(data_df, keywords_df[['DI', 'OpenAlex_KW']], on='DI', how='left')

            merged_df['Author Keywords'] = merged_df.apply(
                lambda row: row['OpenAlex_KW'] if (pd.isna(row['Author Keywords']) or not str(row['Author Keywords']).strip()) and pd.notna(row['OpenAlex_KW']) else row['Author Keywords'],
                axis=1
            )

            merged_df.drop(columns=['OpenAlex_KW'], inplace=True)

            st.success("✅ OpenAlex Keywords merged into 'Author Keywords' column!")
            st.dataframe(merged_df.head())

            output = BytesIO()
            merged_df.to_excel(output, index=False, engine='openpyxl')
            b64 = base64.b64encode(output.getvalue()).decode()
            st.markdown(f'<a href="data:application/octet-stream;base64,{b64}" download="data_with_keywords.xlsx">📥 Download Updated Dataset</a>', unsafe_allow_html=True)
    else:
        st.info("👆 Please upload both the keywords file and your dataset to continue.")
