import streamlit as st
import pandas as pd
import requests
from io import BytesIO
import base64

st.set_page_config(page_title="DOI Extractor by Title", layout="centered")
st.title("🔍 DOI Extractor from Crossref (by Title)")

uploaded_file = st.file_uploader("📄 Upload Excel File (must contain 'TI' column)", type=["xlsx"])

def get_doi_from_crossref(title):
    try:
        url = f"https://api.crossref.org/works?query.title={requests.utils.quote(title)}&rows=1"
        r = requests.get(url, timeout=10)
        if r.status_code == 200:
            items = r.json().get("message", {}).get("items", [])
            if items:
                found_title = items[0].get("title", [""])[0]
                if found_title.strip().lower() == title.strip().lower():
                    return items[0].get("DOI")
        return None
    except:
        return None

if uploaded_file:
    df = pd.read_excel(uploaded_file, engine="openpyxl")

    if "TI" not in df.columns:
        st.error("❌ Uploaded file must contain a 'TI' column (for Title).")
    else:
        if "DI" not in df.columns:
            df["DI"] = pd.NA

        if st.button("🚀 Extract Missing DOIs"):
            with st.spinner("Fetching DOIs from Crossref..."):
                for idx, row in df.iterrows():
                    if pd.isna(row["DI"]) or str(row["DI"]).strip() == "":
                        title = str(row["TI"])
                        doi = get_doi_from_crossref(title)
                        df.at[idx, "DI"] = doi if doi else pd.NA

            st.success("✅ DOI extraction completed!")
            st.dataframe(df.head())

            output = BytesIO()
            df.to_excel(output, index=False, engine="openpyxl")
            b64 = base64.b64encode(output.getvalue()).decode()
            st.markdown(
                f'<a href="data:application/octet-stream;base64,{b64}" download="with_doi.xlsx">📥 Download Updated File</a>',
                unsafe_allow_html=True
            )
