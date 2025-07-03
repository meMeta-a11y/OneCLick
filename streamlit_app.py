import streamlit as st
import pandas as pd
from Bio import Entrez
import io
import time

def get_pmid_from_doi(doi, email):
    """Query PubMed for a DOI and return the PMID."""
    Entrez.email = email  # Set dynamically from user input
    try:
        handle = Entrez.esearch(db="pubmed", term=doi, field="doi")
        record = Entrez.read(handle)
        handle.close()
        if record["IdList"]:
            return record["IdList"][0]
    except Exception as e:
        st.warning(f"Error for DOI {doi}: {e}")
    return None

def process_dataframe(df, email):
    """Adds a PMID column based on the DOI."""
    pmids = []
    for doi in df['DI']:
        if pd.notnull(doi):
            pmid = get_pmid_from_doi(str(doi).strip(), email)
            pmids.append(pmid)
            time.sleep(0.34)  # Respect NCBI rate limit
        else:
            pmids.append(None)
    df['PMID'] = pmids
    return df

# --- Streamlit UI ---
st.title("📘 DOI to PubMed ID Converter")

# Email input
email = st.text_input("Enter your email (required for PubMed access):", placeholder="name@example.com")

# File uploader
uploaded_file = st.file_uploader("Upload your Excel file (.xlsx)", type=["xlsx"])

if email and uploaded_file:
    try:
        df = pd.read_excel(uploaded_file)
        if "DI" not in df.columns:
            st.error("Excel file must contain a column named 'DI' for DOI.")
        else:
            st.info("Processing... please wait.")
            updated_df = process_dataframe(df, email)

            # Display preview
            st.dataframe(updated_df.head())

            # Convert to Excel in memory
            output = io.BytesIO()
            with pd.ExcelWriter(output, engine='openpyxl') as writer:
                updated_df.to_excel(writer, index=False)

            st.success("✅ Processing complete!")
            st.download_button(
                label="📥 Download Excel with PMIDs",
                data=output.getvalue(),
                file_name="output_with_pmids.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            )
    except Exception as e:
        st.error(f"Something went wrong: {e}")
elif uploaded_file and not email:
    st.warning("Please enter your email address before processing.")
