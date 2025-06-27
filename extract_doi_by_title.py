import streamlit as st
import pandas as pd
import requests
from io import BytesIO
import base64
import time

st.set_page_config(page_title="DOI Extractor by Title", layout="centered")
st.title("📘 Extract Missing DOIs by Title (Exact Match via Crossref)")

uploaded_file = st.file_uploader("📤 Upload your Excel file (must contain 'Title' column)", type=["xlsx"])

# --- Function to query Crossref ---
def get_doi_from_crossref(title):
    try:
        url = f"https://api.crossref.org/works?query.bibliographic={requests.utils.quote(title)}&rows=1"
        r = requests.get(url, timeout=10)
        if r.status_code == 200:
            items = r.json().get("message", {}).get("items", [])
            if items:
                item = items[0]
                item_title = item.get('title', [""])[0].strip().lower()
                if item_title == title.strip().lower():
                    return item.get("DOI", "")
        return ""
    except:
        return ""

if uploaded_file:
    df = pd.read_excel(uploaded_file, engine="openpyxl")

    if "Title" not in df.columns:
        st.error("❌ Your file must contain a 'Title' column.")
    else:
        # Ensure 'DI' column exists
        if "DI" not in df.columns:
            df["DI"] = pd.NA

        missing_mask = df["DI"].isna() | (df["DI"].astype(str).str.strip() == "") | (df["DI"].astype(str).str.lower() == "nan")
        missing_titles = df[missing_mask]["Title"].dropna()

        if st.button("🔍 Start DOI Extraction"):
            updated_count = 0
            with st.spinner("Fetching missing DOIs..."):
                for idx in missing_titles.index:
                    title = df.at[idx, "Title"]
                    doi = get_doi_from_crossref(title)
                    if doi:
                        df.at[idx, "DI"] = doi
                        updated_count += 1
                    time.sleep(1)  # avoid Crossref rate limits

            st.success(f"✅ Completed. {updated_count} DOIs were filled.")
            st.dataframe(df)

            # Show rows still missing
            still_missing = df["DI"].isna() | (df["DI"].astype(str).str.strip() == "")
            if still_missing.any():
                st.warning("⚠️ These rows still have missing DOIs:")
                st.dataframe(df[still_missing][["Title"]])

            # Download updated Excel
            output = BytesIO()
            df.to_excel(output, index=False, engine="openpyxl")
            b64 = base64.b64encode(output.getvalue()).decode()
            st.markdown(
                f'<a href="data:application/octet-stream;base64,{b64}" download="doi_filled.xlsx">📥 Download Updated File</a>',
                unsafe_allow_html=True
            )
