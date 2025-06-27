### ✅ `step2_merge_keywords.py` (Clean & Final)

```python
import streamlit as st
import pandas as pd
import base64
from io import BytesIO

st.set_page_config(page_title="Step 2: Merge Keywords", layout="centered")
st.title("🧩 Step 2: Merge OpenAlex Keywords into Dataset")

# Upload input files
uploaded_keywords_file = st.file_uploader("📥 Upload 'openalex_keywords.xlsx' (with columns 'DI' and 'OpenAlex_KW')", type=["xlsx"])
uploaded_data_file = st.file_uploader("📥 Upload your main dataset (must contain 'DI' and 'Author Keywords')", type=["xlsx"])

if uploaded_keywords_file and uploaded_data_file:
    keywords_df = pd.read_excel(uploaded_keywords_file, engine='openpyxl')
    data_df = pd.read_excel(uploaded_data_file, engine='openpyxl')

    # Normalize column names
    keywords_df.columns = keywords_df.columns.str.strip()
    data_df.columns = data_df.columns.str.strip()

    # Validate necessary columns
    if 'DI' not in keywords_df.columns or 'OpenAlex_KW' not in keywords_df.columns:
        st.error("❌ 'openalex_keywords.xlsx' must contain 'DI' and 'OpenAlex_KW' columns.")
    elif 'DI' not in data_df.columns or 'Author Keywords' not in data_df.columns:
        st.error("❌ Main dataset must contain 'DI' and 'Author Keywords' columns.")
    else:
        # Strip spaces and normalize DOI values
        keywords_df['DI'] = keywords_df['DI'].astype(str).str.strip()
        data_df['DI'] = data_df['DI'].astype(str).str.strip()

        # Merge without duplicating rows with missing DI
        merged_df = pd.merge(
            data_df, 
            keywords_df[['DI', 'OpenAlex_KW']], 
            on='DI', 
            how='left', 
            validate='one_to_one'
        )

        # Fill missing 'Author Keywords' only where OpenAlex_KW exists
        merged_df['Author Keywords'] = merged_df.apply(
            lambda row: row['OpenAlex_KW'] if (
                (pd.isna(row['Author Keywords']) or not str(row['Author Keywords']).strip()) 
                and pd.notna(row['DI']) 
                and pd.notna(row['OpenAlex_KW'])
            ) else row['Author Keywords'],
            axis=1
        )

        # Drop helper column
        merged_df.drop(columns=['OpenAlex_KW'], inplace=True)

        st.success("✅ Keywords merged successfully!")
        st.dataframe(merged_df.head())

        # Download link
        output = BytesIO()
        merged_df.to_excel(output, index=False, engine='openpyxl')
        b64 = base64.b64encode(output.getvalue()).decode()
        st.markdown(
            f'<a href="data:application/octet-stream;base64,{b64}" download="merged_dataset.xlsx">📥 Download Merged Dataset</a>',
            unsafe_allow_html=True
        )
else:
    st.info("👆 Please upload both Excel files to proceed.")
```

