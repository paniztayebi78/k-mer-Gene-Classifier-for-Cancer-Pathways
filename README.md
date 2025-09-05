# k-mer-Gene-Classifier-for-Cancer-Pathways
A machine learning project that uses k-mer frequencies and a Random Forest classifier to distinguish between Cell Cycle Regulator and Receptor Tyrosine Kinase gene sequences, two groups critically implicated in cancer.

## 📖 Project Overview
This project implements a supervised machine learning pipeline to classify DNA sequences into two functional groups of cancer-related genes: **Cell Cycle Regulators (CCR)** and **Receptor Tyrosine Kinases (RTK)**. The model uses **k-mer frequency analysis** for feature extraction and a **Random Forest algorithm** for classification, achieving high accuracy in distinguishing between the two gene groups based solely on their sequence composition.

## 🧬 Gene Groups
*   **Cell Cycle Regulators (CCR - Label 1):** CCND1, RB1, CDC25A, CCNE1
*   **Receptor Tyrosine Kinases (RTK - Label 0):** PDGFRA, CD117, FGFR2, AXL

## 🛠️ Methodology
1.  **Data Acquisition:** Nucleotide sequences for each gene group were retrieved from the NCBI database using the `rentrez` R package, filtered for mammalian organisms and a specific length range (1500-5000 bp).
2.  **Data Preprocessing:** Sequences were deduplicated and filtered to remove low-quality entries (those containing ambiguous 'N' bases).
3.  **Feature Engineering:** K-mer frequency analysis (for k=4) was performed on each sequence using the `Biostrings` package to convert raw DNA data into a numerical feature matrix.
4.  **Model Training:** A Random Forest classifier was trained on 80% of the processed data.
5.  **Evaluation:** The model's performance was evaluated on a held-out test set (20% of data) using accuracy, a confusion matrix, precision, recall, and F1-score. Feature importance was plotted to identify the most discriminatory k-mers.

## 📊 Key Results
The Random Forest model successfully classified the gene sequences with high accuracy. The most important k-mers for classification were identified and visualized, providing potential insights into sequence motifs characteristic of each functional group.

## 📁 Files in this Repository
*   `Project_Report.pdf`: The complete project report, including the full R code in the appendix, introduction, methodology, results, and discussion.
*   `kmer_rf_classifier.R`: The original R script containing the code for the entire project pipeline.

## How to Use

1.  Ensure you have R and RStudio installed.
2.  Install the required R packages listed in the script.
3.  Run the `kmer_random_forest_classifier.R` script. The script is designed to first fetch data from NCBI (commented out by default to avoid redundant calls) and then use the local `CCR.fasta` and `RTK.fasta` files for analysis.


## 📚 Dependencies (R Packages)
The project utilized the following R packages:
*   `tidyverse`/`dplyr` (data manipulation)
*   `rentrez` (NCBI data retrieval)
*   `Biostrings` (biological sequence manipulation)
*   `randomForest` (machine learning model)
*   `ggplot2` (visualization)

## 👩‍💻 Author
**Paniz Tayebi** \
*University of Guelph, Bioinformatics*
