# streamlit_doi_extractor.py

import streamlit as st
import pandas as pd
import requests
import difflib
import unicodedata
import string
from io import BytesIO
import base64

st.set_page_config(page_title="DOI Filler via Crossref", layout="centered")
st.title("🔍 Fill Missing DOI using Crossref API")

# --- Helper Functions ---
def clean_title(title):
    if pd.isna(title):
        return ""
    title = str(title).strip().lower()
    title = unicodedata.normalize("NFKD", title)
    title = ''.join(char for char in title if char not in string.punctuation)
    title = ' '.join(title.split())  # Normalize whitespace
    return title

def query_crossref(title, max_results=5):
    try:
        url = f"https://api.crossref.org/works?query.title={requests.utils.quote(title)}&rows={max_results}"
        response = requests.get(url, timeout=10)
        if response.status_code != 200:
            return None
        results = response.json().get("message", {}).get("items", [])
        return results
    except Exception as e:
        return None

def get_best_doi(original_title, results):
    cleaned_original = clean_title(original_title)
    for item in results:
        item_title = item.get("title", [""])[0]
        if not item_title:
            continue
        cleaned_result = clean_title(item_title)
        if difflib.SequenceMatcher(None, cleaned_original, cleaned_result).ratio() > 0.9:
            return item.get("DOI")
    return None

# --- Streamlit Upload ---
uploaded_file = st.file_uploader("📄 Upload Excel File (must have 'Title' and 'DI' columns)", type=["xlsx"])

if uploaded_file:
    df = pd.read_excel(uploaded_file, engine='openpyxl')

    if 'Title' not in df.columns or 'DI' not in df.columns:
        st.error("❌ Excel file must contain both 'Title' and 'DI' columns.")
    else:
        missing_doi_df = df[df['DI'].isna()].copy()
        unmatched_titles = []

        if st.button("🚀 Start DOI Retrieval"):
            with st.spinner("Fetching DOIs from Crossref..."):
                for idx, row in missing_doi_df.iterrows():
                    title = row['Title']
                    results = query_crossref(title)
                    doi = get_best_doi(title, results) if results else None
                    if doi:
                        df.at[idx, 'DI'] = doi
                    else:
                        unmatched_titles.append(title)

            st.success("✅ DOI retrieval completed.")
            st.write(f"Filled {len(missing_doi_df) - len(unmatched_titles)} of {len(missing_doi_df)} missing DOIs")

            # Output results
            output = BytesIO()
            df.to_excel(output, index=False, engine='openpyxl')
            b64 = base64.b64encode(output.getvalue()).decode()
            st.markdown(f'<a href="data:application/octet-stream;base64,{b64}" download="updated_doi_dataset.xlsx">📥 Download Updated Excel</a>', unsafe_allow_html=True)

            if unmatched_titles:
                st.warning("Some titles did not return DOIs. See below:")
                st.dataframe(pd.DataFrame(unmatched_titles, columns=["Unmatched Titles"]))
