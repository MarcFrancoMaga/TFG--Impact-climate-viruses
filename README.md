# Impact-climate-viruses

## 📃 Table of Contents

1. [How It Works](#How-it-works)
2. [Technologies Used](#technologies-used)
3. [Images and Videos](#images-and-videos)
4. [Installation and Deployment Manual](#installation-and-deployment-manual)
5. [Roadmap](#roadmap)
6. [Contributions](#contributions)
7. [License](#license)

## ⚙️ How It Works

This project analyzes potential correlations between phage genomic data (viruses that infect bacteria) and meteorological conditions. Genomic data is sourced from **PhagesDB** and **NOAA**, while weather information is retrieved from **NCBI**. The full workflow is implemented in **Python**, leveraging data science and bioinformatics tools to extract, process, visualize, and correlate information.

### Objectives

- Retrieve genomic data from PhagesDB and NOAA.
- Extract meteorological data from NCBI.
- Identify patterns or correlations between climate variability and genomic changes.
- Generate visual and analytical outputs to support the findings.

## 🚀 Technologies Used

- **Python 3.10+**
- `pandas`, `numpy`: Data processing and analysis
- `biopython`: Genomic data parsing
- `matplotlib`, `seaborn`, `plotly`: Data visualization
- `scikit-learn`: Statistical analysis and modeling
- `requests`, `BeautifulSoup`: Web scraping for data extraction

## 📁 Installation and Deployment Manual

### Prerequisites
- Python 3.10+
- Git

### Steps

1. **Clone the repository**

```bash
git clone https://github.com/MarcFrancoMaga/TFG--Impact-climate-viruses
cd TFG--Impact-climate-viruses
```

2. **Run the scripts or notebooks**

```bash
# Jupyter Notebook (exploration and EDA)
jupyter notebook notebooks/

# Or run the pipeline scripts directly
python main.py
```

### Project Structure

```plaintext
TFG-IMPACT-CLIMATE-VIRUSES/
├── Data_analyze.py       # Analysis and correlation logic
├── Genome_Data.py        # Genomic data acquisition (PhagesDB & NCBI)
├── NOAA_data.py          # Meteorological data processing from NOAA
├── main.py               # Entry point script to run the analysis
├── README.md             # Project documentation
├── requirements.txt      # Python dependencies
├── notebooks/            # Jupyter notebooks for EDA
├── data/                 # Raw and processed datasets
├── results/              # Generated plots and analysis results
```

## 🔄 Roadmap

### Completed ✅
- [x] Extract data from PhagesDB and NOAA
- [x] Integrate meteorological data from NCBI
- [x] Data cleaning and preprocessing
- [x] Data visualization

### In Progress / To Do 🛠️
- [ ] Implement predictive models based on genomic-meteorological patterns
- [ ] Build an interactive dashboard
- [ ] Automate the entire pipeline

## ✍️ Contributions

We welcome contributions from the community! To contribute:
- Fork the repo
- Create a new branch (`git checkout -b feature-xyz`)
- Commit your changes
- Push to your fork
- Create a Pull Request

## 📄 License

This project is licensed under the [MIT License](LICENSE).

---

**Authors:** Your Name / Research Team  
**Contact:** [your.email@example.com](mailto:your.email@example.com)
