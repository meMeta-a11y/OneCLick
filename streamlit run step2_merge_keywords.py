import streamlit as st
import pandas as pd
import base64
from io import BytesIO

st.set_page_config(page_title="Step 2: Merge Keywords", layout="centered")
st.title("🔗 Merge OpenAlex Keywords into Existing Dataset")

# Upload files
uploaded_keywords_file = st.file_uploader("📥 Upload 'openalex_keywords.xlsx' (must contain 'DI' & 'OpenAlex_KW')", type=["xlsx"])
uploaded_data_file = st.file_uploader("📥 Upload 'data.xlsx' (must contain 'DI' & 'Author Keywords')", type=["xlsx"])

if uploaded_keywords_file and uploaded_data_file:
    # Load files
    keywords_df = pd.read_excel(uploaded_keywords_file, engine='openpyxl')
    data_df = pd.read_excel(uploaded_data_file, engine='openpyxl')

    # Clean whitespace and standardize
    keywords_df['DI'] = keywords_df['DI'].astype(str).str.strip()
    data_df['DI'] = data_df['DI'].astype(str).str.strip()

    # Check for duplicates
    dupes_data = data_df['DI'][data_df['DI'].duplicated()].unique()
    dupes_keywords = keywords_df['DI'][keywords_df['DI'].duplicated()].unique()

    if len(dupes_data) > 0:
        st.warning("⚠️ Duplicate DOIs found in data.xlsx:")
        st.write(dupes_data.tolist())

    if len(dupes_keywords) > 0:
        st.warning("⚠️ Duplicate DOIs found in openalex_keywords.xlsx:")
        st.write(dupes_keywords.tolist())

    # Drop duplicates (keep first)
    data_df = data_df.drop_duplicates(subset='DI')
    keywords_df = keywords_df.drop_duplicates(subset='DI')

    # Ignore rows with empty or missing DI in data_df
    data_df = data_df[data_df['DI'].notna() & (data_df['DI'].str.strip() != '')]

    # Merge OpenAlex KW
    merged_df = pd.merge(data_df, keywords_df[['DI', 'OpenAlex_KW']], on='DI', how='left')

    # Fill Author Keywords only if empty
    merged_df['Author Keywords'] = merged_df.apply(
        lambda row: row['OpenAlex_KW'] if (pd.isna(row['Author Keywords']) or str(row['Author Keywords']).strip() == "") and pd.notna(row['OpenAlex_KW']) else row['Author Keywords'],
        axis=1
    )

    merged_df.drop(columns=['OpenAlex_KW'], inplace=True)

    # Display results
    st.success("✅ Merged successfully. OpenAlex keywords added where Author Keywords were missing.")
    st.dataframe(merged_df.head())

    # Download link
    output = BytesIO()
    merged_df.to_excel(output, index=False, engine='openpyxl')
    b64 = base64.b64encode(output.getvalue()).decode()
    st.markdown(
        f'<a href="data:application/octet-stream;base64,{b64}" download="data_with_keywords.xlsx">📥 Download Updated File</a>',
        unsafe_allow_html=True
    )
else:
    st.info("👆 Please upload both required files to continue.")
