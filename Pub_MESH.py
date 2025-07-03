import streamlit as st
from Bio import Entrez
import requests
from bs4 import BeautifulSoup
import time

# --- Streamlit UI ---
st.title("🔍 PubMed ID to MeSH Term Finder")
st.markdown("Enter a PubMed ID (PMID) and your email to retrieve MeSH terms using NCBI + MeSH on Demand.")

# Email input
email = st.text_input("Enter your email (required by NCBI):", placeholder="you@example.com")

# PMID input
pmid = st.text_input("Enter a PubMed ID (PMID):", placeholder="e.g. 35533371")

# Set Entrez email
if email:
    Entrez.email = email


# Functions
def fetch_abstract_from_pubmed(pmid):
    try:
        handle = Entrez.efetch(db="pubmed", id=pmid, rettype="abstract", retmode="text")
        abstract = handle.read()
        handle.close()
        return abstract
    except Exception as e:
        st.error(f"❌ Error fetching abstract: {e}")
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
        st.error(f"❌ Error retrieving MeSH terms: {e}")
        return []


# Submit button
if st.button("🔍 Find MeSH Terms"):
    if not email or not pmid:
        st.warning("Please provide both your email and a valid PubMed ID.")
    else:
        with st.spinner("Fetching abstract and retrieving MeSH terms..."):
            abstract = fetch_abstract_from_pubmed(pmid)
            time.sleep(1)  # be nice to NCBI and MeSH servers
            if abstract:
                st.subheader("📄 Abstract")
                st.text(abstract.strip())
                mesh_terms = get_mesh_terms_from_text(abstract)
                if mesh_terms:
                    st.subheader("🔖 MeSH Terms")
                    for term in mesh_terms:
                        st.markdown(f"- {term}")
                else:
                    st.info("No MeSH terms were returned.")
