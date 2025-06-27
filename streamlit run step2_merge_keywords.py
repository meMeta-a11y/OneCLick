# Ensure DI columns are strings
keywords_df['DI'] = keywords_df['DI'].astype(str).str.strip()
data_df['DI'] = data_df['DI'].astype(str).str.strip()

# Remove rows with empty DI from keywords_df (to prevent false matches)
keywords_df = keywords_df[keywords_df['DI'].notna() & (keywords_df['DI'] != "")]

# Merge (left join) — DI-only match
merged_df = pd.merge(data_df, keywords_df[['DI', 'OpenAlex_KW']], on='DI', how='left')

# Fill 'Author Keywords' if missing AND DI is present AND OpenAlex_KW exists
def fill_kw(row):
    if row['DI'] and pd.notna(row['DI']):
        if (pd.isna(row['Author Keywords']) or str(row['Author Keywords']).strip() == ''):
            if pd.notna(row['OpenAlex_KW']):
                return row['OpenAlex_KW']
    return row['Author Keywords']

merged_df['Author Keywords'] = merged_df.apply(fill_kw, axis=1)

# Drop temporary column
merged_df.drop(columns=['OpenAlex_KW'], inplace=True)
