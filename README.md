# OneClick: Semantic Enrichment for Bibliographic Metadata

[![Streamlit App](https://img.shields.io/badge/launch-app-green)](https://memeta-a11y.github.io/OneCLick/)
[![License: Apache-2.0](https://img.shields.io/badge/license-Apache--2.0-blue.svg)](LICENSE)
[![GitHub Repo](https://img.shields.io/badge/GitHub-OneCLick-lightgrey?logo=github)](https://github.com/meMeta-a11y/OneCLick)

OneClick is an interactive web application for enriching bibliographic metadata using OpenAlex-inferred semantic concepts. Designed with accessibility and usability in mind, it allows researchers—regardless of technical background—to enhance large-scale publication datasets with consistent, machine-generated keywords for bibliometric analysis.

## Key Features

- Semantic Keyword Enrichment: Fills gaps in author-supplied metadata using OpenAlex concept embeddings.
- No-Code Web Interface: Built with Streamlit for instant browser-based use—no programming required.
- Compatible with Bibliometric Tools: Seamlessly integrates enriched metadata into VOSviewer, CiteSpace, and Bibliometrix workflows.
- Cross-Disciplinary Support: Handles heterogeneous corpora, including open-access and multilingual publications.
- FAIR-Aligned: Promotes Findable, Accessible, Interoperable, and Reusable research practices.

## Getting Started

### Requirements

- Python 3.8+
- Internet connection (to access OpenAlex API)

### Installation

```bash
# Clone the repository
git clone https://github.com/meMeta-a11y/OneCLick.git
cd OneCLick

# Install dependencies
pip install -r requirements.txt

# Launch the Streamlit app
streamlit run streamlit_app.py
How It Works
    1. Input: Upload a .csv or .xlsx file containing publication DOIs.
    2. Processing: OneClick queries OpenAlex to fetch semantic concepts associated with each DOI.
    3. Output: The enriched dataset includes standardized, field-normalized keywords for downstream analysis.
You can then export the data for use in co-word mapping, topic modeling, burst detection, and other bibliometric tasks.
Example Use Cases
    • Systematic Reviews: Improve keyword screening in low-metadata corpora.
    • Trend Analysis: Compare inferred vs. author-assigned terms over time.
    • Science of Science: Study conceptual drift and thematic emergence in research fields.
Documentation
The full usage guide, screenshots, and examples are available on our GitHub Pages site: https://github.com/meMeta-a11y/OneCLick

Project Structure
OneCLick/
├── .devcontainer/        # Optional VSCode container config
├── LICENSE               # Apache 2.0 License
├── requirements.txt      # Python dependencies
├── streamlit_app.py      # Main Streamlit application
├── data/                 # Sample input & output files
└── README.md             # This file
Reproducibility
A reproducible snapshot of this repository is available on Zenodo.
To cite OneClick, please use the DOI provided below.
License
This software is released under the Apache 2.0 License. You are free to use, modify, and distribute it with attribution.
Contact & Support
For questions, feature requests, or bug reports:
Email: wynnelkw@gmail.com; lookw@utar.edu.my
GitHub Issues: https://github.com/meMeta-a11y/OneCLick/issues
Citation
meMeta-a11y. (2025). meMeta-a11y/OneCLick: OneClick Bibliometrics (OneClick). Zenodo. https://doi.org/10.5281/zenodo.15876773
