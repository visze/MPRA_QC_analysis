# 🧬 MPRA Analysis Pipeline

This repository contains a reproducible pipeline for analyzing **Massively Parallel Reporter Assay (MPRA)** data.  
It includes Jupyter notebooks and an environment file to make setup easy.  

---

## 🚀 Overview
The pipeline allows you to:
- Preprocess and quality-check MPRA data.  
- Perform statistical analysis of variant activity.  
- Generate plots and summary tables for downstream interpretation.  

It is designed for researchers who want a ready-to-use, reproducible framework for MPRA analysis.  

---

## 📂 Repository Structure
├── notebooks/
│ └── analysis_notebook.ipynb # Main Jupyter notebook
├── scripts/ # Helper scripts for preprocessing and analysis
├── data/ # Input data (not included in repo)
├── results/ # Example outputs
├── environment.yml # Conda environment file
└── README.md # Project documentation

---

## ⚙️ Installation
Clone the repository and set up the environment using **conda**:

```bash
Go to desired path
# Clone the repository
git clone https://github.com/GokhmanLabOrganization/MPRA_QC_analysis.git
cd MPRA_QC_analysis

# Create the conda environment
conda env create -f environment.yml

# Activate the environment
conda activate new_env
```
## ▶️ Usage

Start Jupyter and open the analysis notebook:
```bash
jupyter notebook notebooks/analysis_notebook.ipynb
```
## Contributing – Notebook Cleanliness

This repository uses [**pre-commit**](https://pre-commit.com/) with [**nbstripout**](https://github.com/kynan/nbstripout)  
to keep Jupyter notebooks clean in Git (no `execution_count` churn, no cell outputs).

### First-time setup after cloning
```bash
conda activate pipeline   # or your chosen env
pre-commit install
```
### Usage

On every git commit, notebook outputs and execution counts are stripped automatically.

If pre-commit reports that files were modified:

git add <file>

git commit -m "Clean notebooks"

Manual clean (optional)

Strip all notebooks at once:

pre-commit run --all-files

## 🙌 Acknowledgments

Special thanks to collaborators and prior work that inspired this pipeline.