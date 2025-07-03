import streamlit as st
import pandas as pd
import io

st.set_page_config(page_title="Prepare MeSH Submission")
st.title("📄 Prepare Abstracts for MeSH on Demand")

st.markdown("""
Upload an Excel file with **Title** and **Abstract** columns.  
This app will generate text blocks you can manually paste into [MeSH on Demand](https://meshb.nlm.nih.gov/MeSHonDemand).
""")

uploaded_file = st.file_uploader("Upload Excel File", type=["xlsx"])

if uploaded_file:
    df = pd.read_excel(uploaded_file)

    if "Title" not in df.columns or "Abstract" not in df.columns:
        st.error("The Excel file must contain 'Title' and 'Abstract' columns.")
    else:
        st.success("File loaded! Here's a preview:")
        st.dataframe(df.head())

        blocks = []
        for idx, row in df.iterrows():
            title = str(row.get("Title", "")).strip()
            abstract = str(row.get("Abstract", "")).strip()
            if title and abstract:
                blocks.append(f"{title}\n\n{abstract}")

        st.markdown("### 📋 Text Blocks (Copy and paste to MeSH on Demand)")
        for i, block in enumerate(blocks):
            with st.expander(f"Entry {i+1}"):
                st.text_area("Text block:", value=block, height=200)

        st.markdown("""
        ---  
        🔗 [Click here to open MeSH on Demand](https://meshb.nlm.nih.gov/MeSHonDemand)  
        Paste the text block into the website to get MeSH terms.
        """)
