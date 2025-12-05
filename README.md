# 🧬 Comparative Survey of Eukaryotic Gene Prediction Tools

**Authors:** Mohammed Abraar Khan & Arul Sathya Rajasrinivasan  
**Course:** Bioinformatics  
**Institution:** University of Florida  
**Emails:** mohammed.abraar@ufl.edu, arulsath.rajasri@ufl.edu

---

## 📋 Abstract

This project presents a comparative study of four widely-used ab initio eukaryotic gene prediction tools—**Genscan**, **GlimmerHMM**, **SNAP**, and **AUGUSTUS**—on human genomic DNA segments with known gene annotations.

## 🎯 Objectives

1. Construct a dataset of ~50 human genomic regions with diverse gene structures
2. Develop wrapper scripts to execute multiple gene finders and normalize outputs
3. Define evaluation metrics at exon, gene, and nucleotide levels
4. Compare tools under identical conditions and analyze strengths/weaknesses

## 🛠️ Tools Compared

| Tool | Type | Reference |
|------|------|-----------|
| **AUGUSTUS** | Generalized HMM | Stanke & Waack, 2003 |
| **SNAP** | HMM-based | Korf, 2004 |
| **GlimmerHMM** | Interpolated Markov Model | Salzberg et al., 1999 |
| **Genscan** | Classic HMM | Burge & Karlin, 1997 |

## 📊 Dataset

- **Source:** Human GRCh38 reference genome
- **Annotations:** GENCODE/Ensembl
- **Regions:** ~50 genomic segments
- **Complexity:**
  - Simple (1-2 exons): 10 genes
  - Moderate (3-10 exons): 25 genes
  - Complex (11+ exons): 15 genes

## 📈 Evaluation Metrics

### Exon-Level
- Sensitivity, Precision, F1 Score
- Exact matching and IoU-based matching (≥0.5)

### Gene-Level
- Perfect prediction rate
- Partial match rate

### Nucleotide-Level
- Coding sensitivity
- Non-coding specificity
- Matthews Correlation Coefficient (MCC)

## 🚀 Installation & Usage

### Prerequisites
- Python 3.8+
- No external dependencies required (uses only standard library)

### Quick Start
```bash
# Clone the repository
git clone https://github.com/YOUR_USERNAME/gene-prediction-comparison.git
cd gene-prediction-comparison

# Create virtual environment
python3 -m venv venv
source venv/bin/activate

# Run the pipeline
python3 gene_prediction_project.py

# View results
open visualizations/dashboard.html
```

## 📁 Project Structure
```
gene-prediction-comparison/
├── gene_prediction_project.py    # Main pipeline script
├── README.md                     # This file
├── data/
│   ├── sequences/                # Generated FASTA files
│   ├── annotations/              # Gene annotations
│   └── metadata.json             # Dataset metadata
├── results/
│   ├── augustus/                 # AUGUSTUS predictions
│   ├── snap/                     # SNAP predictions
│   ├── glimmerhmm/               # GlimmerHMM predictions
│   ├── genscan/                  # Genscan predictions
│   └── evaluation_results.json   # Aggregated metrics
└── visualizations/
    └── dashboard.html            # Interactive results dashboard
```

## 📊 Key Findings

1. **AUGUSTUS and SNAP** achieve highest exon-level F1 scores, particularly on multi-exon genes
2. **Performance gap widens** as gene complexity increases
3. **Common error patterns:** missed short exons, boundary shifts, merged/split genes
4. **All tools feasible** for large-scale analysis; trade-off between speed and accuracy

## 📚 References

1. Burge & Karlin (1997). "Prediction of complete gene structures in human genomic DNA." J. Mol. Biol.
2. Salzberg et al. (1999). "Interpolated Markov models for eukaryotic gene finding." Genomics.
3. Korf (2004). "Gene finding in novel genomes." BMC Bioinformatics.
4. Stanke & Waack (2003). "Gene prediction with a hidden Markov model." Bioinformatics.

## 📄 License

This project is for educational purposes as part of the Bioinformatics course at University of Florida.

## 👥 Contributors

- **Mohammed Abraar Khan** - Data acquisition, tool setup, wrapper scripts
- **Arul Sathya Rajasrinivasan** - Evaluation framework, metrics, report
