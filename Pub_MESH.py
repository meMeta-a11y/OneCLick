import streamlit as st
import pandas as pd
from Bio import Entrez
import io
import time

# Configure email for NCBI API access
st.title("🧠 MeSH Term Annotator (PubMed-based)")
st.markdown("Upload your Excel file with columns including **PMID**. This app will retrieve MeSH terms from PubMed and update the **Author Keywords** column.")

email_input = st.text_input("🔐 Enter your email (required for NCBI access):", type="default")

uploaded_file = st.file_uploader("📁 Upload Excel file (.xlsx)", type=["xlsx"])

# NCBI setup
if email_input:
    Entrez.email = email_input.strip()

# Retrieve MeSH terms from NCBI PubMed
def fetch_mesh_terms(pmid):
    try:
        handle = Entrez.efetch(db="pubmed", id=str(pmid), rettype="xml", retmode="text")
        records = Entrez.read(handle)
        handle.close()

        mesh_list = records['PubmedArticle'][0]['MedlineCitation'].get('MeshHeadingList', [])
        terms = [mesh['DescriptorName'] for mesh in mesh_list]
        return "; ".join(terms)
    except Exception as e:
        return None

# Process DataFrame
def process_dataframe(df):
    df_result = df.copy()
    for idx, row in df.iterrows():
        pmid = row.get("PMID", "")
        if pd.isna(pmid) or not str(pmid).isdigit():
            continue

        mesh_terms = fetch_mesh_terms(pmid)
        if mesh_terms:
            df_result.at[idx, "Author Keywords"] = mesh_terms

        time.sleep(0.34)  # Respect NCBI rate limit

    return df_result

# Main app logic
if uploaded_file and email_input:
    try:
        df = pd.read_excel(uploaded_file)

        required_columns = {"Title", "Abstract", "PMID", "Author Keywords"}
        if not required_columns.issubset(df.columns):
            st.error(f"❌ Missing one or more required columns: {required_columns}")
        else:
            st.info("🔄 Processing... Please wait.")
            updated_df = process_dataframe(df)

            # Preview
            st.subheader("✅ Preview of Updated File")
            st.dataframe(updated_df.head())

            # Download button
            output = io.BytesIO()
            with pd.ExcelWriter(output, engine="openpyxl") as writer:
                updated_df.to_excel(writer, index=False)

            st.download_button(
                label="📥 Download Updated Excel",
                data=output.getvalue(),
                file_name="mesh_updated.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            )
    except Exception as e:
        st.error(f"⚠️ Error processing file: {e}")
