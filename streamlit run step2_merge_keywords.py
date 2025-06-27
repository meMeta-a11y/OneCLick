import streamlit as st
import pandas as pd
import base64
from io import BytesIO

st.set_page_config(page_title="Step 2: Merge OpenAlex Keywords", layout="centered")
st.title("🔗 Step 2: Merge OpenAlex Keywords into Dataset")

# File uploads
openalex_file = st.file_uploader("📥 Upload 'openalex_keywords.xlsx' (must contain 'DI', 'OpenAlex_KW')", type=["xlsx"])
data_file = st.file_uploader("📥 Upload 'data.xlsx' (must contain 'DI' and 'Author Keywords')", type=["xlsx"])

if openalex_file and data_file:
    openalex_df = pd.read_excel(openalex_file, engine="openpyxl")
    data_df = pd.read_excel(data_file, engine="openpyxl")

    # Ensure consistent DOI formatting
    openalex_df['DI'] = openalex_df['DI'].astype(str).str.strip().str.lower()
    data_df['DI'] = data_df['DI'].astype(str).str.strip().str.lower()

    # Remove rows with empty DOI in OpenAlex (invalid source rows)
    openalex_df = openalex_df[openalex_df['DI'].notna() & (openalex_df['DI'] != '')]

    # Detect duplicate DOIs in either file
    duplicated_openalex = openalex_df[openalex_df.duplicated('DI', keep=False)]
    duplicated_data = data_df[data_df.duplicated('DI', keep=False)]

    if not duplicated_openalex.empty or not duplicated_data.empty:
        st.warning("⚠️ Duplicate DOIs found. These will be shown below and removed before merging.")

        if not duplicated_openalex.empty:
            st.subheader("🛑 Duplicate DOIs in openalex_keywords.xlsx")
            st.dataframe(duplicated_openalex[['DI']].drop_duplicates())
        if not duplicated_data.empty:
            st.subheader("🛑 Duplicate DOIs in data.xlsx")
            st.dataframe(duplicated_data[['DI']].drop_duplicates())

        # Remove duplicates to avoid merge conflicts
        openalex_df = openalex_df.drop_duplicates(subset='DI')
        data_df = data_df.drop_duplicates(subset='DI')

    # Identify rows in data_df with missing Author Keywords but valid DI
    data_df['Author Keywords'] = data_df['Author Keywords'].astype(str)
    data_df['DI'] = data_df['DI'].astype(str).str.strip().str.lower()

    merged_df = pd.merge(data_df, openalex_df[['DI', 'OpenAlex_KW']], on='DI', how='left')

    # Only fill if Author Keywords is empty and DI exists
    merged_df['Author Keywords'] = merged_df.apply(
        lambda row: row['OpenAlex_KW'] if row['DI'] and (not row['Author Keywords'].strip() or row['Author Keywords'].lower() == 'nan') and pd.notna(row['OpenAlex_KW']) else row['Author Keywords'],
        axis=1
    )
    merged_df.drop(columns=['OpenAlex_KW'], inplace=True)

    st.success("✅ Keywords merged successfully!")
    st.dataframe(merged_df.head())

    # Offer download of final merged file
    output = BytesIO()
    merged_df.to_excel(output, index=False, engine='openpyxl')
    b64 = base64.b64encode(output.getvalue()).decode()
    st.markdown(f'<a href="data:application/octet-stream;base64,{b64}" download="data_with_keywords.xlsx">📥 Download Updated Dataset</a>', unsafe_allow_html=True)

    # Optional: Show rows missing DI
    missing_doi_rows = data_df[data_df['DI'].isna() | (data_df['DI'].str.strip() == '')]
    if not missing_doi_rows.empty:
        st.warning(f"⚠️ {len(missing_doi_rows)} rows have missing DOI. These were ignored during the merge.")

        display_cols = [col for col in ['Author', 'Title', 'Year'] if col in missing_doi_rows.columns]
        st.dataframe(missing_doi_rows[display_cols] if display_cols else missing_doi_rows)

        output_missing = BytesIO()
        missing_doi_rows.to_excel(output_missing, index=False, engine='openpyxl')
        b64_missing = base64.b64encode(output_missing.getvalue()).decode()
        st.markdown(
            f'<a href="data:application/octet-stream;base64,{b64_missing}" download="missing_doi_rows.xlsx">📥 Download Missing DOI Rows</a>',
            unsafe_allow_html=True
        )
else:
    st.info("👆 Please upload both required Excel files to proceed.")
