import streamlit as st
import pandas as pd
import base64
from io import BytesIO

st.set_page_config(page_title="Step 2: Merge Keywords", layout="centered")
st.title("🔄 Step 2: Merge OpenAlex Keywords into Dataset")

uploaded_keywords_file = st.file_uploader("📥 Upload 'openalex_keywords.xlsx' (must contain 'DI' and 'OpenAlex_KW')", type=["xlsx"], key="kw")
uploaded_data_file = st.file_uploader("📥 Upload 'data.xlsx' (must contain 'DI' and 'Author Keywords')", type=["xlsx"], key="data")

if uploaded_keywords_file and uploaded_data_file:
    keywords_df = pd.read_excel(uploaded_keywords_file, engine='openpyxl')
    data_df = pd.read_excel(uploaded_data_file, engine='openpyxl')

    # Validate required columns
    if 'DI' not in keywords_df.columns or 'OpenAlex_KW' not in keywords_df.columns:
        st.error("❌ 'openalex_keywords.xlsx' must contain columns: 'DI' and 'OpenAlex_KW'")
    elif 'DI' not in data_df.columns or 'Author Keywords' not in data_df.columns:
        st.error("❌ 'data.xlsx' must contain columns: 'DI' and 'Author Keywords'")
    else:
        # Strip spaces and normalize for reliable matching
        keywords_df['DI'] = keywords_df['DI'].astype(str).str.strip()
        data_df['DI'] = data_df['DI'].astype(str).str.strip(
