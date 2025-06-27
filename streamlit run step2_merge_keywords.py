# step2_merge_keywords.py

import streamlit as st
import pandas as pd
import base64
from io import BytesIO

st.set_page_config(page_title="Step 2: Merge OpenAlex Keywords", layout="centered")
st.title("🔗 Step 2: Merge OpenAlex Keywords into Dataset")

st.markdown("Upload both files: the OpenAlex keywords file and the main dataset (with DOIs and optional Author Keywords).")

# Upload section
uploaded_keywords_file = st.file_uploader("📥 Upload OpenAlex Keywords File (must contain 'DI' and 'OpenAlex_KW')", type=["xlsx"])
uploaded_data_file = st.file_uploader("📥 Upload Main Dataset (must contain 'DI' and 'Author Keywords')", type=["xlsx"])

if uploaded_keywords_file and uploaded_data_file:
    keywords_df = pd.read_excel(uploaded_keywords_file, engine='openpyxl')
    data_df = pd.read_excel(uploaded_data_file, engine='openpyxl')

    # Validate required columns
    if 'DI' not in keywords_df.columns or 'OpenAlex_KW' not in keywords_df.columns:
        st.error("❌ OpenAlex keywords file must contain 'DI' and 'OpenAlex_KW'.")
    elif 'DI' not in data_df.columns or 'Author Keywords' not in data_df.columns:
        st.error("❌ Main dataset must contain 'DI' and 'Author Keywords'.")
    else:
        # Merge OpenAlex keywords into the main dataset without dropping rows
        merged_df = pd.merge(
            data_df,
            keywords_df[['DI', 'OpenAlex_KW']],
            on='DI',
            how='left'
        )

        # Update 'Author Keywords' only if empty AND DI exists AND OpenAlex_KW is available
        def fill_keywords(row):
            if pd.notna(row['DI']) and (pd.isna(row['Author Keywords']) or not str(row['Author Keywords']).strip()):
                return row['OpenAlex_KW'] if pd.notna(row['OpenAlex_KW']) else row['Author Keywords']
            return row['Author Keywords']

        merged_df['Author Keywords'] = merged_df.apply(fill_keywords, axis=1)

        # Drop OpenAlex_KW column after use
        merged_df.drop(columns=['OpenAlex_KW'], inplace=True)

        st.success("✅ Missing keywords filled where possible — no rows were deleted.")
        st.dataframe(merged_df.head())

        # Prepare download
        output = BytesIO()
        merged_df.to_excel(output, index=False, engine='openpyxl')
        b64 = base64.b64encode(output.getvalue()).decode()
        st.markdown(
            f'<a href="data:application/octet-stream;base64,{b64}" download="data_with_keywords.xlsx">📥 Download Updated Dataset</a>',
            unsafe_allow_html=True
        )
else:
    st.info("👆 Please upload both Excel files to continue.")
