# streamlit_app.py

import streamlit as st
import pandas as pd
import requests
import base64
import time
import logging
from io import BytesIO, StringIO

# ---------------- Streamlit Config ----------------
st.set_page_config(page_title="OneCLick Keyword Retriever", layout="centered")
st.title("🚀 OneClick Keyword Retriever")

uploaded = st.file_uploader("📁 Upload Excel (.xlsx) with 'DI' (DOI)", type=["xlsx"])

# ---------------- Logging ----------------
logging.basicConfig(filename="error_log.txt", level=logging.INFO, format="%(asctime)s - %(message)s")

# ---------------- Fetch Keywords with Backoff ----------------
def fetch_openalex_keywords(doi, max_retries=3, backoff_factor=2):
    """Fetch OpenAlex keywords for a DOI with retry, backoff, and error logging."""
    url = f"https://api.openalex.org/works/doi:{doi}"
    retries = 0

    while retries < max_retries:
        try:
            r = requests.get(url, timeout=10)

            if r.status_code == 200:
                data = r.json()
                keywords = data.get("keywords", [])

                # Explicitly flag missing keywords
                if not keywords:
                    logging.warning(f"No keywords found for DOI {doi}")
                    return "[WARNING:NoKeywords]"

                return '; '.join(sorted(set(k['display_name'].strip() for k in keywords)))

            elif r.status_code == 429:  # Too many requests
                wait_time = int(r.headers.get("Retry-After", 2 ** retries))
                logging.warning(f"Rate limited for DOI {doi}. Retrying in {wait_time}s.")
                time.sleep(wait_time)

            else:
                logging.error(f"Failed DOI {doi} with status {r.status_code}")
                return f"[ERROR:{r.status_code}]"

        except requests.exceptions.Timeout:
            logging.error(f"Timeout error for DOI {doi}")
            return "[ERROR:Timeout]"

        except Exception as e:
            logging.error(f"Unexpected error for DOI {doi}: {str(e)}")
            return f"[ERROR:{str(e)}]"

        retries += 1

    return "[ERROR:MaxRetries]"

# ---------------- File Download Helpers ----------------
def generate_excel_download(df, filename="openalex_keywords.xlsx"):
    output = BytesIO()
    df.to_excel(output, index=False, engine='openpyxl')
    b64 = base64.b64encode(output.getvalue()).decode()
    return f'<a href="data:application/octet-stream;base64,{b64}" download="{filename}">📥 Download {filename}</a>'

def generate_csv_download(df, filename="openalex_keywords.csv"):
    output = StringIO()
    df.to_csv(output, index=False, encoding='utf-8-sig')
    b64 = base64.b64encode(output.getvalue().encode()).decode()
    return f'<a href="data:file/csv;base64,{b64}" download="{filename}">📥 Download {filename}</a>'

# ---------------- Main App Flow ----------------
if uploaded:
    df = pd.read_excel(uploaded, engine='openpyxl')

    if 'DI' not in df.columns:
        st.error("❌ Column 'DI' (DOI) is required.")
    else:
        df['OpenAlex_KW'] = ""

        if st.button("🚀 Retrieve Keywords from OpenAlex"):
            with st.spinner("Fetching from OpenAlex... please wait"):
                for idx, row in df.iterrows():
                    if pd.notna(row['DI']):
                        df.at[idx, 'OpenAlex_KW'] = fetch_openalex_keywords(str(row['DI']).strip())

            st.success("✅ Process completed!")

            # Preview clean dataset
            st.dataframe(df[['DI', 'OpenAlex_KW']].head())

            # --- Separate errors and warnings ---
            error_df = df[df['OpenAlex_KW'].str.contains("ERROR", na=False)].copy()
            warning_df = df[df['OpenAlex_KW'].str.contains("WARNING", na=False)].copy()

            # Clean enriched dataset (remove error/warning markers but keep traceable logs)
            enriched_df = df.copy()
            enriched_df.loc[enriched_df['OpenAlex_KW'].str.contains("ERROR|WARNING", na=False), 'OpenAlex_KW'] = ""

            # --- Downloads ---
            st.markdown("### 📄 Download Enriched Dataset")
            st.markdown("The enriched dataset includes a new `OpenAlex_KW` column. Errors and warnings are logged separately.")
            st.markdown(generate_excel_download(enriched_df), unsafe_allow_html=True)
            st.markdown(generate_csv_download(enriched_df), unsafe_allow_html=True)

            if not error_df.empty:
                st.error(f"❌ {len(error_df)} DOIs failed. Download error log below.")
                st.markdown(generate_csv_download(error_df, filename="Failed_DOIs.csv"), unsafe_allow_html=True)

            if not warning_df.empty:
                st.warning(f"⚠️ {len(warning_df)} DOIs had no keywords in OpenAlex. Download warning log below.")
                st.markdown(generate_csv_download(warning_df, filename="NoKeyword_DOIs.csv"), unsafe_allow_html=True)
