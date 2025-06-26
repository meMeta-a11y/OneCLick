# streamlit_app.py

import streamlit as st
import pandas as pd
import difflib
import requests
import base64
from io import BytesIO

st.set_page_config(page_title="One-Click Keyword Cleaner", layout="centered")
st.title("🚀 One-Click Keyword Processor")

uploaded = st.file_uploader("📁 Upload Excel (.xlsx) with 'DI' (DOI)", type=["xlsx"])

# --- Step 1: OpenAlex Retrieval ---
def fetch_openalex_keywords(doi):
    try:
        r = requests.get(f"https://api.openalex.org/works/doi:{doi}")
        if r.status_code != 200:
            return ""
        data = r.json()
        keywords = data.get("keywords", [])
        return '; '.join(sorted(set(k['display_name'].strip() for k in keywords)))
    except Exception:
        return ""

# --- Step 2: Clean Keyword Inconsistencies ---
def clean_keywords_column(df, colname='Final_KW', threshold=0.9):
    all_kws = set()
    for kwlist in df[colname].fillna(""):
        all_kws.update([k.strip() for k in kwlist.split('; ') if k.strip()])

    corr_map = {}
    kws = list(all_kws)
    for i in range(len(kws)):
        for j in range(i + 1, len(kws)):
            if difflib.SequenceMatcher(None, kws[i], kws[j]).ratio() >= threshold:
                corr_map[kws[j]] = kws[i]

    def correct(val):
        return '; '.join(sorted(set([corr_map.get(k.strip(), k.strip()) for k in val.split('; ') if k.strip()])))
    
    df['KW_CLEANED'] = df[colname].apply(correct)
    return df

# --- Main App Flow ---
if uploaded:
    df = pd.read_excel(uploaded, engine='openpyxl')

    if 'DI' not in df.columns:
        st.error("❌ 'DI' column (DOI) is required.")
    else:
        keyword_col = None
        for col in ['Keywords', 'DE']:
            if col in df.columns:
                keyword_col = col
                break

        if keyword_col:
            df['Final_KW'] = df[keyword_col].fillna("")
        else:
            df['Final_KW'] = ""

        if st.button("🚀 One-Click Process (Fetch + Clean + Export)"):
            with st.spinner("Processing..."):
                for idx, row in df.iterrows():
                    if not row['Final_KW'].strip() and pd.notna(row['DI']):
                        df.at[idx, 'Final_KW'] = fetch_openalex_keywords(row['DI'])
                df = clean_keywords_column(df, 'Final_KW')

            st.success("✅ All Done!")
            st.dataframe(df[['DI', 'Final_KW', 'KW_CLEANED']].head())

            output = BytesIO()
            df.to_excel(output, index=False, engine='openpyxl')
            b64 = base64.b64encode(output.getvalue()).decode()
            st.markdown(f'<a href="data:application/octet-stream;base64,{b64}" download="cleaned_keywords.xlsx">📅 Click to Download Cleaned File</a>', unsafe_allow_html=True)
