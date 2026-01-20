# Non-Local Gravitational Leakage: A Spectral Propagator Ansatz from the Dark Dimension

**Author:** Ibrahim Gul  
**Correspondence:** gulpoetika@gmail.com  
**Status:** Under Review (**EPJ C**)

## 📌 Overview
This repository contains the source code, data, and numerical analysis for the paper **"Non-Local Gravitational Leakage: A Spectral Propagator Ansatz from the Dark Dimension"**. 

The model introduces a massive graviton spectral propagator that facilitates gravitational leakage into a dark dimension at low acceleration scales. This code reproduces the key results, demonstrating that the model resolves galactic and cosmological tensions without particle dark matter.

## 🏆 Key Results
1.  **Galactic Dynamics:** Bayesian inference on the SPARC dataset favors the leakage model ($n > 1$) over standard MOND ($n=1$) with high statistical significance.
2.  **Extended Sample Analysis:** While the primary paper focuses on the core **118-galaxy** "Golden Sample," this repository includes an updated analysis of **125 high-quality galaxies**. The results ($n \approx 2.15$) remain consistent and provide even stronger evidence for gravitational leakage.
3.  **Cosmology ($S_8$):** The provided `class_mu_patch.diff` implements the model in **CLASS**, yielding $S_8 \approx 0.78$, alleviating the $S_8$ tension.

## ⚠️ Data Selection (118 vs 125 Galaxies)
The primary manuscript is based on a strict 118-galaxy sample. In this repository, the sample has been expanded to **125 galaxies** that satisfy the quality criteria (Inclination $i > 30^\circ$ and Quality Flag $Q=1$). This extended analysis serves as a robustness check, confirming that the physical leakage signature is not sensitive to slight variations in the sample size and that the model's predictive power remains superior over larger datasets.

## 📂 Repository Structure
* **`analyze.py`**: The main Bayesian analysis script using the Dynesty sampler.
* **`sparc_raw_full_3389_points.csv`**: Raw rotation curve data from the SPARC database.
* **`sparc_subset_118_list.txt`**: List of Galaxy IDs used for filtering (contains 125 unique entries).
* **`chain_results.txt` / `analysis_summary.txt`**: Output data from the latest 125-galaxy analysis run.
* **`class_mu_patch.diff`**: Patch for the CLASS Boltzmann code to implement the non-local leakage propagator.
* **`Figure_1` - `Figure_5`**: Official visualizations presented in the manuscript.

## 🚀 Usage
1.  Install requirements:  
    `pip install -r requirements.txt`
2.  Run analysis:  
    `python analyze.py`

## 📜 Citation
```bibtex
@article{Gul2026Leakage,
  title={Non-Local Gravitational Leakage: A Spectral Propagator Ansatz from the Dark Dimension},
  author={Gul, Ibrahim},
  journal={arXiv preprint / Submitted to EPJ C},
  year={2026}
}
