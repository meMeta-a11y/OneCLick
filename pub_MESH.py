import streamlit as st
from Bio import Entrez
import requests
from bs4 import BeautifulSoup
import pandas as pd
import io
import time

# --- UI Header ---
st.title("🔍 Bulk PubMed ID to MeSH Term Annotator")
st.markdown("Upload an Excel file with a column `PMID`, enter your email, and get MeSH terms in the `Author Keywords` column.")

# --- User Email Input ---
email = st.text_input("Enter your email (required by NCBI):", placeholder="you@example.com")

# --- File Upload ---
uploaded_file = st.file_uploader("Upload Excel file (.xlsx)", type=["xlsx"])

# --- Set Entrez Email ---
if email:
    Entrez.email = email


# --- Functions ---

def fetch_abstract_from_pubmed(pmid):
    try:
        handle = Entrez.efetch(db="pubmed", id=pmid, rettype="abstract", retmode="text")
        abstract = handle.read()
        handle.close()
        return abstract
    except Exception as e:
        st.warning(f"Error fetching abstract for PMID {pmid}: {e}")
        return None

def get_mesh_terms_from_text(text):
    try:
        url = "https://meshb.nlm.nih.gov/MeSHonDemand"
        response = requests.post(url, data={"search": text})
        soup = BeautifulSoup(response.text, "html.parser")

        mesh_terms = []
        table = soup.find("table", {"id": "meshResultTable"})
        if table:
            for row in table.find_all("tr")[1:]:  # Skip header
                cols = row.find_all("td")
                if len(cols) >= 2:
                    mesh_terms.append(cols[1].get_text(strip=True))
        return mesh_terms
    except Exception as e:
        st.warning(f"Error retrieving MeSH terms: {e}")
        return []

def process_dataframe(df):
    mesh_list = []
    for i, row in df.iterrows():
        pmid = str(row.get("PMID")).strip()
        if pd.notnull(pmid):
            abstract = fetch_abstract_from_pubmed(pmid)
            time.sleep(1)  # be gentle with servers
            if abstract:
                mesh_terms = get_mesh_terms_from_text(abstract)
                mesh_string = "; ".join(mesh_terms)
                mesh_list.append(mesh_string)
            else:
                mesh_list.append("")
        else:
            mesh_list.append("")
    df["Author Keywords"] = mesh_list
    return df


# --- Main Logic ---
if uploaded_file and email:
    try:
        df = pd.read_excel(uploaded_file)

        if "PMID" not in df.columns:
            st.error("❌ Excel must contain a column named 'PMID'.")
        else:
            st.info("Processing... this may take a few minutes.")
            processed_df = process_dataframe(df)

            st.success("✅ Done! MeSH terms added to 'Author Keywords' column.")
            st.dataframe(processed_df.head())

            # Convert to Excel
            output = io.BytesIO()
            with pd.ExcelWriter(output, engine='openpyxl') as writer:
                processed_df.to_excel(writer, index=False)
            st.download_button(
                label="📥 Download Updated Excel",
                data=output.getvalue(),
                file_name="mesh_keywords_output.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            )
    except Exception as e:
        st.error(f"Something went wrong: {e}")
elif uploaded_file and not email:
    st.warning("Please enter your email to use NCBI services.")
