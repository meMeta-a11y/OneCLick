import streamlit as st
import pandas as pd
import spacy
import scispacy
from scispacy.linking import EntityLinker
import io

st.set_page_config(page_title="ScispaCy MeSH Term Extractor")
st.title("🧠 Biomedical Keyword Extractor using ScispaCy")

# Load model and linker (cache to avoid reloading)
@st.cache_resource
def load_model():
    nlp = spacy.load("en_core_sci_sm")
    linker = EntityLinker(resolve_abbreviations=True, name="umls")
    nlp.add_pipe("scispacy_linker", config={"linker_name": "umls"})
    return nlp

nlp = load_model()

# Upload file
uploaded_file = st.file_uploader("Upload Excel File with Title and Abstract columns", type=["xlsx"])

if uploaded_file:
    df = pd.read_excel(uploaded_file)
    
    if not {"Title", "Abstract", "Author Keywords"}.issubset(df.columns):
        st.error("Excel must contain 'Title', 'Abstract', and 'Author Keywords' columns.")
    else:
        limit = st.number_input("Limit rows to process", min_value=1, max_value=len(df), value=min(10, len(df)))
        
        updated_keywords = []
        for i, row in df.head(limit).iterrows():
            title = str(row.get("Title", "")).strip()
            abstract = str(row.get("Abstract", "")).strip()

            if not abstract:
                updated_keywords.append("")
                continue

            doc = nlp(f"{title} {abstract}")
            terms = set()
            for ent in doc.ents:
                for umls_ent in ent._.kb_ents:
                    cui = umls_ent[0]
                    name = nlp.get_pipe("scispacy_linker").kb.cui_to_entity[cui].canonical_name
                    terms.add(name)
            keyword_str = "; ".join(sorted(terms))
            updated_keywords.append(keyword_str)

        df["Author Keywords"] = updated_keywords
        st.success("✅ Keywords updated successfully.")
        st.dataframe(df.head())

        output = io.BytesIO()
        with pd.ExcelWriter(output, engine='openpyxl') as writer:
            df.to_excel(writer, index=False)
        st.download_button("📥 Download Updated Excel", data=output.getvalue(),
                           file_name="keywords_scispacy.xlsx",
                           mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
