import streamlit as st
from Bio import Entrez
import requests
from bs4 import BeautifulSoup
import pandas as pd
import io
import time

# --- UI ---
st.title("🔎 PubMed ID to MeSH Term Annotator")
st.markdown("Upload an Excel file containing a `PMID` column. The app will populate the `Author Keywords` column with MeSH terms.")

# --- User email input ---
email = st.text_input("Enter your email (required for NCBI access):", placeholder="you@example.com")

# --- File Upload ---
uploaded_file = st.file_uploader("Upload Excel file (.xlsx)", type=["xlsx"])

# --- Set Entrez email if provided ---
if email:
    Entrez.email = email

# --- Get Abstract from PubMed ---
def fetch_abstract_from_pubmed(pmid):
    try:
        handle = Entrez.efetch(db="pubmed", id=pmid, rettype="abstract", retmode="text")
        abstract = handle.read()
        handle.close()
        return abstract
    except Exception as e:
        st.warning(f"⚠️ Error fetching abstract for PMID {pmid}: {e}")
        return None

# --- Get MeSH terms from MeSH on Demand ---
def get_mesh_terms_from_text(text):
    try:
        url = "https://meshb.nlm.nih.gov/MeSHonDemand"
        response = requests.post(url, data={"search": text})
        soup = BeautifulSoup(response.text, "html.parser")

        mesh_terms = []
        table = soup.find("table", {"id": "meshResultTable"})
        if not table:
            return []

        rows = table.find_all("tr")[1:]  # skip header
        for row in rows:
            cols = row.find_all("td")
            if len(cols) >= 2:
                term = cols[1].get_text(strip=True)
                if term:
                    mesh_terms.append(term)
        return mesh_terms
    except Exception as e:
        st.warning(f"⚠️ Error retrieving MeSH terms: {e}")
        return []

# --- Process DataFrame ---
def process_dataframe(df):
    # Copy to avoid changing original
    df_result = df.copy()
    mesh_terms_list = []

    for i, row in df.iterrows():
        pmid = str(row.get("PMID")).strip() if pd.notnull(row.get("PMID")) else None

        if pmid and pmid.isdigit():
            abstract = fetch_abstract_from_pubmed(pmid)
            time.sleep(1)  # Respect NCBI rate limits

            if abstract and len(abstract.strip()) > 50:
                mesh_terms = get_mesh_terms_from_text(abstract)
                mesh_text = "; ".join(mesh_terms) if mesh_terms else ""
            else:
                mesh_text = ""
        else:
            mesh_text = ""  # leave blank if PMID missing or invalid

        mesh_terms_list.append(mesh_text)

    df_result["Author Keywords"] = mesh_terms_list
    return df_result

# --- Main logic ---
if uploaded_file and email:
    try:
        df = pd.read_excel(uploaded_file)

        if "PMID" not in df.columns:
            st.error("❌ The Excel file must contain a column named 'PMID'.")
        else:
            st.info("🔄 Processing... This may take a few seconds per PMID.")
            updated_df = process_dataframe(df)

            # Preview result
            st.success("✅ Completed MeSH annotation!")
            st.dataframe(updated_df.head())

            # Download link
            output = io.BytesIO()
            with pd.ExcelWriter(output, engine="openpyxl") as writer:
                updated_df.to_excel(writer, index=False)

            st.download_button(
                label="📥 Download Updated Excel File",
                data=output.getvalue(),
                file_name="mesh_annotated_output.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            )

    except Exception as e:
        st.error(f"❌ Error processing file: {e}")

elif uploaded_file and not email:
    st.warning("⚠️ Please enter your email to proceed (required by NCBI).")
