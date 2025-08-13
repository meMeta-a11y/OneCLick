# OneClick: Streamlined Keyword Metadata Enrichment Using OpenAlex  
![GitHub release (latest by date)](https://img.shields.io/github/v/release/meMeta-a11y/OneCLick)  
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15876773.svg)](https://doi.org/10.5281/zenodo.15876773)  

## Overview  
OneClick is an interactive web app for keyword metadata enrichment using machine-inferred keywords from the OpenAlex '/works, keywords' field. It enables researchers to augment incomplete publication records—especially those missing author-assigned keywords—via a no-code interface. OneClick supports scalable, FAIR-compliant bibliometric workflows.  

## How It Works  
Input: Upload an '.xlsx' file containing a column of DOIs.  
Process: The system validates the file structure, then queries the OpenAlex API for each DOI to retrieve keywords generated through OpenAlex’s hierarchical topics system. These keywords are deduplicated, sorted, concatenated into a semicolon-separated string, and appended as a new column ('OpenAlex_KW') while preserving the original dataset structure.  
Output: Returns enriched metadata containing the assigned keywords. The keywords are provided without additional metadata such as confidence scores or hierarchical taxonomy levels, which are reserved for Future Work.  

[Watch usage demo (MP4)]
(https://github.com/meMeta-a11y/OneCLick/blob/main/OneCLick.mp4)  
(https://www.youtube.com/watch?v=NBFrPhuk67Y)

Sample input/output files:  
- ['sample_input.xlsx'](https://github.com/meMeta-a11y/OneCLick/blob/389148f5bce7c6bc934479c38bede7518182d932/input%20file.xlsx)  
- ['sample_output.xlsx'](https://github.com/meMeta-a11y/OneCLick/blob/389148f5bce7c6bc934479c38bede7518182d932/openalex_keywords.xlsx)  

---

## Reproducibility  
- Reproducible Capsule (Zenodo): [https://doi.org/10.5281/zenodo.15876773](https://doi.org/10.5281/zenodo.15876773)  
- Release Tag: ['meMeta-a11y/OneCLick: OneClick Bibliometrics'](https://github.com/meMeta-a11y/OneCLick/releases)  

---

## Tech Stack  
- Framework: [Streamlit](https://streamlit.io) v1.46.1  
- Language: Python  
- API: [OpenAlex](https://openalex.org)  
- License: Apache-2.0  
- Code Versioning: git  
- Repository: [GitHub – meMeta-a11y/OneCLick](https://github.com/meMeta-a11y/OneCLick)  

---

## Documentation  
The interface is designed for intuitive use. All functions are documented within the code, and dependencies are listed in 'requirements.txt'. Terminology and data field descriptions are fully aligned with the manuscript and reflect the current OpenAlex implementation.  

---

## Contact  
For questions or support, email: wynnelkw@gmail.com; lookw@utar.edu.my
