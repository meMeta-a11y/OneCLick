import streamlit as st
import pandas as pd
from io import BytesIO
import base64

st.set_page_config(page_title="Step 2: Merge OpenAlex Keywords", layout="centered")
st.title("🔗 Step 2: Merge Keywords into Existing Dataset")

# Upload files
openalex_file = st.file_uploader("📄 Upload 'openalex_keywords.xlsx' (must contain 'DI' and 'OpenAlex_KW')", type=["xlsx"])
main_file = st.file_uploader("📄 Upload your dataset (must contain 'DI' and 'Author Keywords')", type=["xlsx"])

if openalex_file and main_file:
    # Load files
    openalex_df = pd.read_excel(openalex_file, engine='openpyxl')
    data_df = pd.read_excel(main_file, engine='openpyxl')

    # Normalize DOIs for matching
    openalex_df['DI'] = openalex_df['DI'].astype(str).str.strip().str.lower()
    data_df['DI'] = data_df['DI'].astype(str).str.strip().str.lower()

    # ✅ Check for duplicated DOIs in openalex file
    duplicated_openalex = openalex_df[openalex_df.duplicated(subset='DI', keep=False)]
    if not duplicated_openalex.empty:
        st.warning("⚠️ Duplicated DOIs found in 'openalex_keywords.xlsx':")
        st.dataframe(duplicated_openalex[['DI']].drop_duplicates())
        # Keep only first occurrence
        openalex_df = openalex_df.drop_duplicates(subset='DI', keep='first')

    # ✅ Check for duplicated DOIs in main data
    duplicated_main = data_df[data_df.duplicated(subset='DI', keep=False)]
    if not duplicated_main.empty:
        st.info("ℹ️ Note: Duplicated DOIs found in your dataset — they will be preserved, and keywords applied per row.")
        st.dataframe(duplicated_main[['DI']].drop_duplicates())

    # Create DOI-to-keywords mapping
    openalex_dict = dict(zip(openalex_df['DI'], openalex_df['OpenAlex_KW']))

    # Function to fill keywords where missing
    def fill_keywords(row):
        doi = row['DI']
        current_kw = str(row.get('Author Keywords', '')).strip()
        if doi and (not current_kw or current_kw.lower() == 'nan'):
            return openalex_dict.get(doi, current_kw)
        return current_kw

    # Apply keyword filling
    data_df['Author Keywords'] = data_df.apply(fill_keywords, axis=1)

    # Restore original column formatting (optional)
    data_df['DI'] = data_df['DI'].str.upper()

    st.success("✅ Keywords successfully merged!")
    st.dataframe(data_df.head())

    # Export result
    output = BytesIO()
    data_df.to_excel(output, index=False, engine='openpyxl')
    b64 = base64.b64encode(output.getvalue()).decode()
    st.markdown(
        f'<a href="data:application/octet-stream;base64,{b64}" download="data_with_keywords.xlsx">📥 Download Updated File</a>',
        unsafe_allow_html=True
    )
