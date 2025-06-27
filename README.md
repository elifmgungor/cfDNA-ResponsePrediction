# cfDNA-ResponsePrediction

## Project Overview

This repository contains a bioinformatics pipeline designed to predict treatment response in breast cancer patients using nucleosome positioning signals extracted from circulating cell-free DNA (cfDNA) sequencing data.

### Scientific Context

Circulating cell-free DNA (cfDNA) consists of fragmented DNA found in the bloodstream, primarily released through apoptosis or necrosis. While in healthy individuals cfDNA mainly originates from hematopoietic cells, in cancer patients a significant fraction comes from tumor cells.

Previous studies have demonstrated that nucleosome positioning signals in cfDNA can accurately determine the tissue of origin of these DNA fragments.

### Project Objective

In this project, our goal is not to infer the tissue-of-origin of cfDNA.  
Instead, we aim to evaluate whether nucleosome-derived signals can serve as predictive biomarkers of treatment response **before the initiation of therapy**, using a simple blood sample and without performing any detailed cell-of-origin analysis.

Specifically, we analyze cfDNA sequencing data to determine if signal patterns around selected genomic regions (e.g., Transcription Start Sites and Estrogen Response Elements) can distinguish between responder and non-responder patients.

---

## Pipeline Summary

The pipeline covers the following main steps:

1. **Data Preprocessing**
   - Conversion of raw FASTQ files to BAM format.
   - Quality filtering of low-confidence regions.

2. **Signal Extraction**
   - Calculation of nucleosome-associated signals: Fragment Depth (FD) and Windowed Protection Score (WPS).

3. **Normalization**
   - Evaluation and comparison of different normalization strategies.
   - Generation of normalized feature matrices for modeling.

4. **Supervised Learning**
   - Training and evaluation of multiple machine learning models (KNN, Decision Tree, SVM, Logistic Regression, Random Forest).
   - Comparison of model performances.

---

## Repository Structure


├── src/ # Python source code (preprocessing, signal extraction, normalization, modeling)

├── scripts/ # Bash scripts for pipeline execution

├── notebooks/ # Jupyter notebooks for exploratory analysis and result visualization

├── data/ # Raw and intermediate data (excluded from Git)

├── results/ # Model outputs and figures (excluded from Git)

├── docs/ # Additional documentation

├── requirements.txt # Python dependencies

├── LICENSE

└── README.md

---
## Installation

Clone this repository:

```bash
git clone https://github.com/elifmgungor/cfDNA-ResponsePrediction.git
cd cfDNA-ResponsePrediction
```

## Usage

Run each step of the pipeline sequentially using the provided Python and Bash scripts:

1. **Preprocessing:**  
Conversion of FASTQ files to BAM format and quality filtering of low-confidence regions.

2. **Signal Extraction:**  
Calculation of nucleosome-related signals such as Fragment Depth (FD) and Windowed Protection Score (WPS).

3. **Normalization:**  
Testing and application of different normalization strategies.

4. **Modeling:**  
Training and evaluation of classification models using the normalized features.

For detailed execution instructions, please refer to the specific scripts and notebooks within each directory.

---

## Data and Results Disclaimer

⚠️ **Important:**  
Raw sequencing data and generated results (e.g., BAM files, feature matrices, figures, trained models) are **not included** in this repository due to confidentiality and storage limitations.

---

## License

This project is licensed under the MIT License – see the [LICENSE](LICENSE) file for details.
