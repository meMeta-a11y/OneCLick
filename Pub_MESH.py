import streamlit as st
from Bio import Entrez
import requests
from bs4 import BeautifulSoup
import pandas as pd
import io
import time
import xml.etree.ElementTree as ET

# --- Streamlit UI ---
st.title("🔎 PubMed ID to MeSH Term Annotator")
st.markdown("Upload an Excel file with a `PMID` column. The app retrieves both title and abstract, and adds MeSH terms in the `Author Keywords` column.")

# --- Email Input ---
email = st.text_input("Enter your email (required for NCBI access):", placeholder="you@example.com")

# --- File Upload ---
uploaded_file = st.file_uploader("Upload Excel file (.xlsx)", type=["xlsx"])

# --- Set Entrez Email ---
if email:
    Entrez.email = email

# --- Fetch Title + Abstract using XML ---
def fetch_title_abstract(pmid):
    try:
        handle = Entrez.efetch(db="pubmed", id=pmid, rettype="medline", retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        if not records or 'PubmedArticle' not in records or not records['PubmedArticle']:
            return None

        article = records['PubmedArticle'][0]['MedlineCitation']['Article']
        title = article.get('ArticleTitle', "")
        abstract = ""
        if 'Abstract' in article and 'AbstractText' in article['Abstract']:
            abstract_parts = article['Abstract']['AbstractText']
            if isinstance(abstract_parts, list):
                abstract = " ".join(str(part) for part in abstract_parts)
            else:
                abstract = str(abstract_parts)

        return f"{title} {abstract}".strip()
    except Exception as e:
        st.warning(f"⚠️ Error retrieving article for PMID {pmid}: {e}")
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
    df_result = df.copy()
    mesh_terms_list = []

    for i, row in df.iterrows():
        pmid = str(row.get("PMID")).strip() if pd.notnull(row.get("PMID")) else None

        if pmid and pmid.isdigit():
            full_text = fetch_title_abstract(pmid)
            time.sleep(1)  # Respect NCBI rate limits

            if full_text and len(full_text) > 50:
                mesh_terms = get_mesh_terms_from_text(full_text)
                mesh_text = "; ".join(mesh_terms) if mesh_terms else ""
            else:
                mesh_text = ""
        else:
            mesh_text = ""

        mesh_terms_list.append(mesh_text)

    df_result["Author Keywords"] = mesh_terms_list
    return df_result

# --- Main App Logic ---
if uploaded_file and email:
    try:
        df = pd.read_excel(uploaded_file)

        if "PMID" not in df.columns:
            st.error("❌ Excel file must contain a column named 'PMID'.")
        else:
            st.info("🔄 Processing each PMID...")
            updated_df = process_dataframe(df)

            # Preview
            st.success("✅ MeSH term annotation completed.")
            st.dataframe(updated_df.head())

            # Output
            output = io.BytesIO()
            with pd.ExcelWriter(output, engine="openpyxl") as writer:
                updated_df.to_excel(writer, index=False)

            st.download_button(
                label="📥 Download MeSH-Annotated Excel",
                data=output.getvalue(),
                file_name="mesh_annotated_output.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            )
    except Exception as e:
        st.error(f"❌ Error processing file: {e}")

elif uploaded_file and not email:
    st.warning("⚠️ Please enter your email before processing.")
