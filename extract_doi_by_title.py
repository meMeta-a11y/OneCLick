import streamlit as st
import pandas as pd
import requests
import time
from io import BytesIO
import base64

st.set_page_config(page_title="DOI Filler via Crossref", layout="centered")
st.title("🔍 Fill Missing DOIs using Crossref (Exact Title Match)")

uploaded_file = st.file_uploader("📄 Upload Excel file (must contain a 'Title' column)", type=["xlsx"])

# --- Crossref lookup function ---
def get_doi_from_crossref(title):
    try:
        url = f"https://api.crossref.org/works?query.bibliographic={requests.utils.quote(title)}&rows=1"
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            items = response.json().get("message", {}).get("items", [])
            if items:
                item = items[0]
                if 'title' in item and title.strip().lower() == item['title'][0].strip().lower():
                    return item.get("DOI", "")
        return ""
    except Exception:
        return ""

# --- Process Logic ---
if uploaded_file:
    df = pd.read_excel(uploaded_file, engine='openpyxl')

    if 'Title' not in df.columns:
        st.error("❌ The uploaded file must contain a column named 'Title'.")
    else:
        if 'DI' not in df.columns:
            df['DI'] = pd.NA  # create DOI column if missing

        missing_doi_mask = df['DI'].isna() | (df['DI'].astype(str).str.strip().str.lower() == 'nan') | (df['DI'].astype(str).str.strip() == '')
        missing_titles = df.loc[missing_doi_mask, 'Title'].dropna()

        if st.button("🚀 Fetch Missing DOIs"):
            updated = 0
            with st.spinner("Searching DOIs on Crossref..."):
                for idx in missing_titles.index:
                    title = df.at[idx, 'Title']
                    doi = get_doi_from_crossref(title)
                    if doi:
                        df.at[idx, 'DI'] = doi
                        updated += 1
                    time.sleep(1)  # prevent rate-limiting

            st.success(f"✅ Finished. {updated} DOIs updated.")
            st.dataframe(df)

            # Highlight rows still missing DOI
            missing_final = df['DI'].isna() | (df['DI'].astype(str).str.strip() == "")
            if missing_final.any():
                st.warning("⚠️ Some rows still missing DOI:")
                st.dataframe(df[missing_final][['Title']])

            # Download link
            output = BytesIO()
            df.to_excel(output, index=False, engine='openpyxl')
            b64 = base64.b64encode(output.getvalue()).decode()
            st.markdown(
                f'<a href="data:application/octet-stream;base64,{b64}" download="filled_doi.xlsx">📥 Download Result</a>',
                unsafe_allow_html=True
            )
