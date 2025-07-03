import streamlit as st
import pandas as pd
from Bio import Entrez
import io
import time

st.title("🧠 MeSH Term Annotator from Title & Abstract")
st.markdown("Upload an Excel file with `Title`, `Abstract`, and `Author Keywords`. The app will extract MeSH headings using PubMed XML and update the `Author Keywords` column.")

email = st.text_input("Enter your email (required for NCBI API):", placeholder="you@example.com")
uploaded_file = st.file_uploader("Upload Excel (.xlsx)", type=["xlsx"])

if email:
    Entrez.email = email.strip()

def fetch_mesh_terms_from_text(title, abstract):
    text = f"{title.strip()} {abstract.strip()}"
    # Use PubMed search for the title to get PMID
    try:
        search_handle = Entrez.esearch(db="pubmed", term=title, retmax=1)
        search_record = Entrez.read(search_handle)
        search_handle.close()
        pmid = search_record["IdList"][0] if search_record["IdList"] else None
        if not pmid:
            return ""
        
        fetch_handle = Entrez.efetch(db="pubmed", id=pmid, rettype="medline", retmode="text")
        medline_text = fetch_handle.read()
        fetch_handle.close()
        
        mesh_terms = []
        for line in medline_text.splitlines():
            if line.startswith("MH  - "):
                mesh_terms.append(line.replace("MH  - ", "").strip())
        return "; ".join(mesh_terms)
    except Exception as e:
        st.warning(f"⚠️ Error extracting MeSH for '{title[:50]}...': {e}")
        return ""

def process_dataframe(df):
    df_result = df.copy()
    for idx, row in df.iterrows():
        title = str(row.get("Title", "")).strip()
        abstract = str(row.get("Abstract", "")).strip()
        if title and abstract:
            st.write(f"🔍 Processing row {idx+1}: {title[:60]}...")
            mesh = fetch_mesh_terms_from_text(title, abstract)
            if mesh:
                df_result.at[idx, "Author Keywords"] = mesh
            time.sleep(0.34)
    return df_result

if uploaded_file and email:
    df = pd.read_excel(uploaded_file)
    required_cols = {"Title", "Abstract", "Author Keywords"}
    if not required_cols.issubset(df.columns):
        st.error(f"❌ Required columns: {required_cols}")
    else:
        updated = process_dataframe(df)
        st.success("✅ Annotation complete!")
        st.dataframe(updated.head())
        output = io.BytesIO()
        with pd.ExcelWriter(output, engine="openpyxl") as writer:
            updated.to_excel(writer, index=False)
        st.download_button("📥 Download Updated Excel", data=output.getvalue(),
                           file_name="mesh_terms_updated.xlsx",
                           mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
elif uploaded_file:
    st.warning("⚠️ Please enter your email to proceed.")
