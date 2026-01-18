# Non-Local Gravitational Leakage: A Spectral Propagator Ansatz from the Dark Dimension

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![arXiv](https://img.shields.io/badge/arXiv-240X.XXXXX-b31b1b.svg)](https://arxiv.org/)
[![Status](https://img.shields.io/badge/Status-Submitted-green.svg)]()

**Author:** Ibrahim Gul  
**Correspondence:** [gulpoetika@gmail.com](mailto:gulpoetika@gmail.com)

## 📌 Overview

This repository contains the numerical analysis pipeline, data subsets, and plotting scripts for the paper **"Non-Local Gravitational Leakage: A Spectral Propagator Ansatz from the Dark Dimension"**.

We investigate a UV-motivated modified gravity model embedded within a 5D Einstein-Dilaton framework. This code reproduces the key results of the paper, including:
1.  **Galactic Dynamics:** Bayesian inference on SPARC galaxies favoring the leakage model ($Z \approx 3.81\sigma$).
2.  **Bullet Cluster:** WKB tunneling calculations explaining the ~200 kpc offset.
3.  **Cosmology:** Modified growth structure predicting $S_8 \approx 0.78$.

## 📂 Repository Structure

* `analyze.py`: Main analysis script performing Bayesian inference (Nested Sampling) on galaxy data.
* `generate_final_figures.py`: Script to reproduce Figures 1-5 as seen in the manuscript.
* `data/`: Contains the processed SPARC subset (`sparc_subset_118.csv`).
* `chains/`: Contains posterior chain files from `dynesty`.
* `class_mu_patch.diff`: Diff file to patch the Boltzmann code `CLASS` for non-local gravity (optional for reproduction of cosmological plots).

## 🚀 Installation & Usage

### Prerequisites
The analysis requires Python 3.8+ and the following libraries:

```bash
pip install -r requirements.txt
