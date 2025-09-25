# 🚀 OneCLick Keyword Retriever

**OneCLick** is a lightweight, open-source **Streamlit application** for enriching bibliometric datasets by retrieving **missing keywords** from the [OpenAlex API](https://openalex.org).  

It streamlines preprocessing for tools like **VOSviewer**, **CiteSpace**, and **Bibliometrix**, and outputs **reproducible, analysis-ready files** for downstream research.

---

## ✨ Features

- 📁 **Upload Input:**  
  - Accepts `.xlsx` files with a column named `DI` (DOI).

- 🔄 **DOI-Based Keyword Retrieval:**  
  - Sends DOI queries to OpenAlex and retrieves enriched metadata keywords.

- ⚠️ **Error & Warning Handling:**  
  - Invalid DOIs, network timeouts, or rate-limit errors → `Failed_DOIs.csv`  
  - Valid DOIs without keywords → flagged as `[WARNING:NoKeywords]` in `NoKeyword_DOIs.csv`

- 📊 **Latency Logging:**  
  - Per-DOI request latency is tracked and summarized (median & IQR)  
  - Full request log saved as `Latency_Log.csv`

- 📤 **Outputs:**  
  - Enriched dataset → `openalex_keywords.csv` / `.xlsx`  
  - Downloadable error & warning logs  
  - Latency logs for reproducibility and benchmarking

---

## 🛠️ Installation

```bash
# Clone the repository
git clone https://github.com/meMeta-a11y/oneclick.git
cd oneclick

# Install dependencies
pip install -r requirements.txt
```

Run the Streamlit application:

```bash
streamlit run streamlit_app.py
```

---

## 📊 Example Datasets

We provide three benchmark datasets for reproducibility and testing:

| Dataset       | DOIs | Purpose                           |
|--------------|------|-----------------------------------|
| Pilot        | 37   | Quick exploratory testing        |
| Intermediate | 100  | Lightweight benchmarking         |
| Validation   | 1,097| Large-scale robustness testing   |

All datasets are available in the [`example_data/`](./example_data/) directory and are permanently archived on [Zenodo](https://doi.org/10.5281/zenodo.15876773).

---

## ⚡ Performance Benchmarking

Validation experiments demonstrate **robustness and scalability** across dataset sizes:

| Dataset       | DOIs | Error DOIs | No-Keyword DOIs | Success Rate | Median Latency (s/DOI, IQR) |
|--------------|------|------------|------------------|---------------|----------------------------|
| Pilot        | 37   | 0          | 0                | 100.00%       | 0.30 (0.28 – 0.31)         |
| Intermediate | 100  | 0          | 26               | 74.00%        | 0.29 (0.28 – 0.31)         |
| Validation   | 1,097| 0          | 205              | 81.31%        | 0.29 (0.28 – 0.30)         |

> 🧪 Benchmarks were conducted on Streamlit Cloud (Aug 2025) using a Linux environment (4 GB RAM, 100 Mbps network). Performance may vary depending on network conditions and OpenAlex API rate limits.

---

## 🚀 Quick-Start: VOSviewer Workflow

Follow these steps to integrate **OneCLick** outputs with **VOSviewer**:

1. **Open VOSviewer** → *Create → Map based on network data*  
2. **Input:** Choose your enriched `openalex_keywords.csv`  
3. **Mapping:** Set `OpenAlex_KW` as the keyword field  
4. **Options:** Select co-occurrence counting method (e.g., *Full counting*)  
5. **Generate:** Visualize the keyword network

📚 Additional workflows are described in the manuscript (Section 2.12).

---

## 🧭 Limitations & Future Work

- ❌ **Topic hierarchies** and **confidence scores** are not yet exported — these are planned for future releases.  
- 🧠 Integration into bibliometric **clustering pipelines** (e.g., theme assignment) is under active development.

---

## 📁 Repository Structure

```
oneclick/
├── streamlit_app.py         # Main application
├── requirements.txt         # Dependencies
├── example_data/            # Example DOI datasets
│   ├── pilot_37.xlsx
│   ├── intermediate_100.xlsx
│   └── validation_1097.xlsx
├── outputs/                 # Example outputs & logs
│   ├── openalex_keywords.csv
│   ├── Failed_DOIs.csv
│   ├── NoKeyword_DOIs.csv
│   └── Latency_Log.csv
└── README.md               # Documentation
```

---

## 📚 Citation

If you use **OneCLick** in your research, please cite it using one of the following formats:

### 📘 APA  
Wei, L. K. (2025). *OneCLick: Streamlined metadata enrichment using machine-inferred keywords from OpenAlex*. *SoftwareX, 31*, 102353.

### 📙 Chicago  
Wei, Loo Keat. "OneCLick: Streamlined metadata enrichment using machine-inferred keywords from OpenAlex." *SoftwareX* 31 (2025): 102353.

### 📗 Harvard  
Wei, L.K., 2025. OneCLick: Streamlined metadata enrichment using machine-inferred keywords from OpenAlex. *SoftwareX*, 31, p.102353.

### 📕 Vancouver  
Wei LK. OneCLick: Streamlined metadata enrichment using machine-inferred keywords from OpenAlex. *SoftwareX*. 2025 Sep 1;31:102353.

---

## 🌐 About

- 🔗 **Live Demo:** [metakey.streamlit.app](https://metakey.streamlit.app/)  
- 📦 **Source Code:** [GitHub Repository](https://github.com/meMeta-a11y/oneclick)

---

## 📜 License

This project is licensed under the **Apache-2.0 License** - see the [LICENSE](./LICENSE) file for details.

---

💡 **OneCLick** is developed to accelerate bibliometric research workflows by automating metadata enrichment and preparing datasets for network analysis and visualization. Contributions, issues, and pull requests are always welcome!
