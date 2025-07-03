import streamlit as st
import pandas as pd
import requests
import time
import io

st.title("🔍 MeSH Term Annotator (Title + Abstract)")
st.markdown("This tool uses NLM's **Medical Text Indexer (MTI)** to extract MeSH terms based on your **Title and Abstract**. No PMID required!")

uploaded_file = st.file_uploader("📁 Upload Excel file (.xlsx)", type=["xlsx"])

API_URL = "https://lhncbc.nlm.nih.gov/ii/tools/MTI/Submitter.cgi"

def get_mesh_from_text(title, abstract):
    text = f"{title.strip()} {abstract.strip()}"
    try:
        response = requests.post(
            API_URL,
            data={
                "input": text,
                "tool": "MTI",
                "format": "json",
                "email": "anonymous@anonymous.com"  # Optional
            },
            timeout=30
        )
        if response.status_code == 200:
            data = response.json()
            terms = [entry["preferred"] for entry in data.get("MeshHeadings", [])]
            return "; ".join(sorted(set(terms)))
        else:
            return ""
    except Exception as e:
        return ""

def process_dataframe(df):
    df_result = df.copy()
    for idx, row in df.iterrows():
        title = str(row.get("Title", "")).strip()
        abstract = str(row.get("Abstract", "")).strip()

        if not title or not abstract:
            continue

        st.write(f"🔎 Processing row {idx+1}: {title[:50]}...")
        mesh_terms = get_mesh_from_text(title, abstract)
        df_result.at[idx, "Author Keywords"] = mesh_terms
        time.sleep(1)  # Be kind to the server

    return df_result

if uploaded_file:
    try:
        df = pd.read_excel(uploaded_file)

        required_columns = {"Title", "Abstract", "Author Keywords"}
        if not required_columns.issubset(df.columns):
            st.error(f"❌ Excel file must include columns: {required_columns}")
        else:
            st.info("Processing, please wait...")
            updated_df = process_dataframe(df)

            st.success("✅ Done!")
            st.dataframe(updated_df.head())

            output = io.BytesIO()
            with pd.ExcelWriter(output, engine="openpyxl") as writer:
                updated_df.to_excel(writer, index=False)

            st.download_button(
                label="📥 Download Updated File",
                data=output.getvalue(),
                file_name="mesh_terms_updated.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            )

    except Exception as e:
        st.error(f"⚠️ Error reading file: {e}")
