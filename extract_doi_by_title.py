import pandas as pd
import requests
import time
import urllib.parse
from difflib import SequenceMatcher

def get_doi_from_title(title, threshold=0.95):
    encoded_title = urllib.parse.quote_plus(title)
    url = f"https://api.crossref.org/works?query.title={encoded_title}&rows=1"

    try:
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            items = response.json().get("message", {}).get("items", [])
            if items:
                returned_title = items[0].get("title", [""])[0]
                score = SequenceMatcher(None, title.lower(), returned_title.lower()).ratio()
                
                if score >= threshold:
                    return items[0].get("DOI", None), returned_title, score
                else:
                    return None, returned_title, score
    except Exception as e:
        print(f"Error retrieving DOI for title '{title}': {e}")
    
    return None, None, 0.0

def update_doi_in_excel(input_excel, output_excel):
    df = pd.read_excel(input_excel)

    if "MatchScore" not in df.columns:
        df["MatchScore"] = None
    if "ReturnedTitle" not in df.columns:
        df["ReturnedTitle"] = None

    for i, row in df[df['DI'].isna()].iterrows():
        title = row['TI']
        print(f"🔍 Searching DOI for: {title}")
        doi, returned_title, score = get_doi_from_title(title, threshold=0.95)

        df.at[i, "MatchScore"] = round(score, 3)
        df.at[i, "ReturnedTitle"] = returned_title

        if doi:
            df.at[i, 'DI'] = doi
            print(f"✅ High match ({score:.2f}) → DOI: {doi}")
        else:
            print(f"❌ Match too low ({score:.2f}) → Skipped")

        time.sleep(1)  # Respect Crossref API rate limit

    df.to_excel(output_excel, index=False)
    print(f"\n✅ Excel updated and saved as: {output_excel}")

# Example usage
if __name__ == "__main__":
    input_file = "your_dataset.xlsx"
    output_file = "doi_results_high_match.xlsx"
    update_doi_in_excel(input_file, output_file)
