# step2_merge_keywords.py

import streamlit as st
import pandas as pd
import base64
from io import BytesIO

st.set_page_config(page_title="Keyword Merger", layout="centered")
st.title("🔗 Step 2: Merge Keywords into Existing Dataset")

# Upload both files
uploaded_keywords_file = st.file_uploader("📥 Upload 'openalex_keywords.xlsx' (must contain 'DI' and 'OpenAlex_KW')", type=["xlsx"], key="kw")
uploaded_data_file = st.file_uploader("📥 Upload 'data.xlsx' (must contain 'DI' and 'Author Keywords')", type=["xlsx"], key="data")

if uploaded_keywords_file and uploaded_data_file:
    keywords_df = pd.read_excel(uploaded_keywords_file, engine='openpyxl')
    data_df = pd.read_excel(uploaded_data_file, engine='openpyxl')

    # Validate column existence
    if 'DI' not in keywords_df.columns or 'OpenAlex_KW' not in keywords_df.columns:
        st.error("❌ 'openalex_keywords.xlsx' must contain 'DI' and 'OpenAlex_KW' columns.")
    elif 'DI' not in data_df.columns or 'Author Keywords' not in data_df.columns:
        st.error("❌ 'data.xlsx' must contain 'DI' and 'Author Keywords' columns.")
    else:
        # Merge and fill missing keywords
        merged_df = pd.merge(data_df, keywords_df[['DI', 'OpenAlex_KW']], on='DI', how='left')

        merged_df['Author Keywords'] = merged_df.apply(
            lambda row: row['OpenAlex_KW'] if (pd.isna(row['Author Keywords']) or not str(row['Author Keywords']).strip()) and pd.notna(row['OpenAlex_KW']) else row['Author Keywords'],
            axis=1
        )

        merged_df.drop(columns=['OpenAlex_KW'], inplace=True)

        st.success("✅ Keywords successfully merged into 'Author Keywords' column!")
        st.dataframe(merged_df.head())

        # Export the result
        output = BytesIO()
        merged_df.to_excel(output, index=False, engine='openpyxl')
        b64 = base64.b64encode(output.getvalue()).decode()
        st.markdown(f'<a href="data:application/octet-stream;base64,{b64}" download="data_with_keywords.xlsx">📥 Download Updated File</a>', unsafe_allow_html=True)
else:
    st.info("👆 Please upload both required Excel files to continue.")
