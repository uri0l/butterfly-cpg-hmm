# Distribution of CpG Islands on a Butterfly Genome

> A computational genomics project developed as part of the *Algorithms in Biology* course at Universitat Politècnica de Barcelona, May 2025

[![Python](https://img.shields.io/badge/Language-Python-blue.svg)](https://www.python.org/)
[![Bioinformatics](https://img.shields.io/badge/Field-Bioinformatics-green.svg)](https://en.wikipedia.org/wiki/Bioinformatics)
[![HMM](https://img.shields.io/badge/Method-HMM%20%7C%20Viterbi-orange.svg)](https://en.wikipedia.org/wiki/Hidden_Markov_model)
[![Genomics](https://img.shields.io/badge/Data-NCBI%20Genomes-lightgrey.svg)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_933228805.2/)

---

## Overview

This project implements a Hidden Markov Model (HMM)-based tool to detect and characterize CpG islands across chromosomes of *Hipparchia semele* (the grayling butterfly), using the chromosome-level genome assembly [GCA_933228805.2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_933228805.2/) from NCBI.

CpG islands are short genomic regions with a higher-than-expected frequency of CpG dinucleotides. They are commonly found near gene promoters and play a key role in epigenetic regulation of gene expression through DNA methylation. The project focuses on comparing CpG island distributions between the sex chromosomes (Z and W) and the autosome chromosome 28, combining a sliding window detection approach with HMM-based probabilistic modelling.

---

## Goals

### General Goal
Implement an HMM-based tool to determine the location of CpG islands in *Hipparchia semele* and compare their distribution between sex chromosomes and an autosome.

### Specific Goals
- Extract chromosome-level FASTA sequences for chromosomes Z, W, and 28 from the full genome assembly.
- Detect CpG islands using a sliding window approach with two detection thresholds (10% and 15%).
- Map and visualize the genomic location and distribution of CpG islands across chromosomes.
- Model CpG island detection as an HMM with two hidden states (CpG / non-CpG) and compute transition and emission probability matrices.
- Perform a comparative analysis of CpG island density between sex chromosomes and the autosome.

---

## Dataset

- **Source:** *Hipparchia semele* genome assembly [GCA_933228805.2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_933228805.2/) (NCBI).
- **Chromosomes analyzed:** Z (16,045,015 bp, 36.5% GC), W (8,523,649 bp, 38.5% GC), and 28 (6,818,840 bp, 38% GC).
- **Format:** FASTA (`.fna`).

The genome assembly file (`GCA_933228805.2_ilHipSeme1.2_genomic.fna`) and the extracted chromosome FASTA files are not included in this repository due to file size. They can be downloaded directly from NCBI at the link above. Once downloaded, run `read_first.py` to extract the chromosomes of interest.

---

## Methods & Pipeline

### 1. Chromosome Extraction
`read_first.py` parses the full multi-FASTA genome assembly and extracts the sequences for chromosomes Z, W, and 28 into separate FASTA files (`chrZ.fasta`, `chrW.fasta`, `chr28.fasta`) based on chromosome descriptors in the sequence headers.

### 2. CpG Island Detection - Sliding Window
A non-overlapping sliding window approach (window size = 300–500 bp) scans each chromosome and classifies each window as CpG or non-CpG based on the fraction of CG dinucleotides. Two thresholds were tested:
- **10%** - higher sensitivity, detects weaker CpG islands.
- **15%** - higher specificity, focuses on well-defined CpG islands.

### 3. CpG Island Mapping & Visualization
Individual nucleotide positions belonging to CpG windows are extracted and plotted as event plots, showing the distribution of CpG islands along each chromosome. Positions are normalized to chromosome length (%) to allow direct comparison across chromosomes of different sizes.

### 4. Hidden Markov Model
Each nucleotide is assigned a hidden state: **C** (CpG island region) or **N** (non-CpG region), derived from the sliding window results. Two probability matrices are then computed:
- **Transition matrix:** probability of moving from one state to the next (C→C, C→N, N→C, N→N).
- **Emission matrix:** probability of observing each nucleotide (A, C, G, T) given each hidden state.

A full Viterbi-based HMM implementation (`cpg_hmm_complete.py`) is also provided for comparison with the sliding window method.

### 5. Comparative Analysis
CpG island counts and densities (islands per kb) are compared across the three chromosomes. A bar plot on a logarithmic scale visualizes the stark difference between CpG and non-CpG window counts. Results are discussed in the context of the biological differences between sex chromosomes and autosomes.

---

## Key Results

- CpG islands are very sparse across all three chromosomes, regardless of the detection threshold.
- Chromosome 28 had the highest absolute number of CpG islands, followed by Z and then W.
- When normalized for chromosome length, chromosome W showed the highest CpG island density despite having the lowest absolute count.
- The transition and emission matrices confirmed clear compositional differences between CpG and non-CpG regions across all chromosomes.
- The absence of a genome annotation prevented mapping CpG islands to specific genes or promoters, which remains a key limitation and direction for future work.

---

## Repository Structure

```
.
├── scripts/
│   ├── read_first.py                      # step 1 — chromosome extraction
│   ├── comparison_cpg_vs_non__10__.py     # step 2 — CpG vs non-CpG, 10% threshold
│   ├── comparison_cpg_vs_non__15__.py     # step 2 — CpG vs non-CpG, 15% threshold
│   ├── locations_10.py                    # step 3 — island positions + HMM matrices, 10%
│   ├── locations_normalized_10.py         # step 3 — normalized location plot, 10%
│   ├── locations_normalized_15.py         # step 3 — normalized location plot, 15%
│   └── cpg_hmm_complete.py               # step 4 — full HMM class with Viterbi decoder
├── plots/
│   ├── cpg_island_comparison_10.png
│   ├── cpg_island_comparison_15.png
│   ├── location_zw28_10.png
│   └── ...
└── README.md
```

> All figures generated by the pipeline are saved to the `plots/` folder and committed to the repository. If you encounter any issues running the code locally, you can browse the pre-generated outputs there directly.

---

## Dependencies

```bash
pip install biopython matplotlib numpy scipy pandas
```

---

## How to Run

1. Clone this repository.
2. Download the genome assembly from [NCBI](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_933228805.2/) and place `GCA_933228805.2_ilHipSeme1.2_genomic.fna` in the project root.
3. Run `read_first.py` to extract the chromosome FASTA files:
```bash
python read_first.py
```
4. Run any of the analysis scripts:
```bash
python ZW28_locations_remastered_10__.py
python ZW28_comparison_cpg_vs_non__15__.py
```

---

## Limitations

- The sliding window approach uses a fixed window size, which cannot guarantee that detected windows correspond to continuous CpG islands.
- HMM states are derived from sliding window labels rather than from a ground-truth annotation, making the model semi-supervised.
- No genome annotation is available for *H. semele*, preventing direct mapping of CpG islands to genes or promoters.
- Binary classification (CpG / non-CpG) does not account for repeated or mutated sequences.

---

## Reference

NCBI Genome Assembly: *Hipparchia semele* - GCA_933228805.2. https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_933228805.2/
