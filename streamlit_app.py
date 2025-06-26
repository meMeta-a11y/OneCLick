# streamlit_app.py

import streamlit as st
import pandas as pd
import difflib
import requests
import base64
import re
import unicodedata
from io import BytesIO

import spacy
from nltk.corpus import stopwords
from nltk import download as nltk_download

# Initialize NLP tools
nlp = spacy.load("en_core_web_sm")
nltk_download('stopwords')
stop_words = set(stopwords.words("english"))

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

# --- Step 2: Merge and Clean Keywords ---
def merge_and_clean_keywords(existing, fetched):
    merged = set()
    for kw in (existing + "; " + fetched).split(';'):
        cleaned = kw.strip()
        if cleaned:
            merged.add(cleaned)
    return '; '.join(sorted(merged))

# --- Step 3: Harmonize Keywords ---
def harmonize_keywords(df, colname='Final_KW', threshold=0.9):
    spelling_map = {
        "organisation": "organization",
        "behaviour": "behavior",
        "analyse": "analyze",
        "colour": "color",
        "modelling": "modeling",
    }

    synonym_map = {
        "ai": "artificial intelligence",
        "ml": "machine learning",
        "dl": "deep learning",
        "covid-19": "covid",
        "vr": "virtual reality",
    }

    def normalize_kw(kw):
        kw = re.sub(r"\\(([^)]+)\\)", "", kw)  # Remove content in brackets
        kw = unicodedata.normalize("NFKD", kw).encode("ascii", "ignore").decode("utf-8")
        kw = kw.lower().strip()
        for uk, us in spelling_map.items():
            kw = kw.replace(uk, us)
        if kw.endswith('ies'):
            kw = kw[:-3] + 'y'
        elif kw.endswith('s') and not kw.endswith('ss'):
            kw = kw[:-1]
        kw = ' '.join([w for w in kw.split() if w not in stop_words])
        kw = synonym_map.get(kw, kw)
        doc = nlp(kw)
        kw = ' '.join([token.lemma_ for token in doc])
        return kw

    all_kws = set()
    for kwlist in df[colname].fillna(""):
        all_kws.update([normalize_kw(k) for k in kwlist.split(';') if k.strip()])

    corr_map = {}
    kws = list(all_kws)
    for i in range(len(kws)):
        for j in range(i + 1, len(kws)):
            if difflib.SequenceMatcher(None, kws[i], kws[j]).ratio() >= threshold:
                corr_map[kws[j]] = kws[i]

    def correct(val):
        corrected = set()
        for k in val.split(';'):
            norm = normalize_kw(k)
            corrected.add(corr_map.get(norm, norm))
        return '; '.join(sorted(corrected))

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

        if st.button("🚀 One-Click Process (Fetch + Sync + Clean + Export)"):
            with st.spinner("Processing..."):
                for idx, row in df.iterrows():
                    if pd.notna(row['DI']):
                        fetched = fetch_openalex_keywords(row['DI'])
                        df.at[idx, 'Final_KW'] = merge_and_clean_keywords(row['Final_KW'], fetched)
                df = harmonize_keywords(df, 'Final_KW')

            st.success("✅ All Done!")
            st.dataframe(df[['DI', 'Final_KW', 'KW_CLEANED']].head())

            output = BytesIO()
            df.to_excel(output, index=False, engine='openpyxl')
            b64 = base64.b64encode(output.getvalue()).decode()
            st.markdown(f'<a href="data:application/octet-stream;base64,{b64}" download="cleaned_keywords.xlsx">📅 Click to Download Cleaned File</a>', unsafe_allow_html=True)
