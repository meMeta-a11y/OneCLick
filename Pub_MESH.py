def process_dataframe(df):
    mesh_terms_list = []

    for i, row in df.iterrows():
        pmid = str(row.get("PMID")).strip()

        if pd.notnull(pmid) and pmid.isdigit():
            abstract = fetch_abstract_from_pubmed(pmid)
            time.sleep(1)

            if abstract and len(abstract.strip()) > 50:  # skip too-short abstracts
                mesh_terms = get_mesh_terms_from_text(abstract)
                mesh_text = "; ".join(mesh_terms) if mesh_terms else ""
            else:
                mesh_text = ""
        else:
            continue  # Skip rows without a valid PMID

        mesh_terms_list.append(mesh_text)

    # Align the new column with the original DataFrame (for matching index)
    df = df[df["PMID"].notnull() & df["PMID"].astype(str).str.isdigit()].copy()
    df["Author Keywords"] = mesh_terms_list
    return df
