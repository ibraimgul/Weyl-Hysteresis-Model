"""
Non-Local Gravity / Dark Dimension Leakage Analysis
===================================================
Author: Ibrahim Gul
Description:
  Bayesian analysis of SPARC rotation curves using Dark Dimension leakage model.
"""

import numpy as np
import pandas as pd
import sys
import os

try:
    import dynesty
    from dynesty import utils as dyfunc
except ImportError:
    print("Error: 'dynesty' library not found. Please install it.")
    sys.exit(1)

# --- CONFIGURATION & CONSTANTS ---
KPC_TO_M = 3.086e19
KM_TO_M = 1000.0
CONV_FACTOR = (KM_TO_M**2) / KPC_TO_M

UPSILON_DISK = 0.5
UPSILON_BULGE = 0.7

# --- DOSYA YOLLARI GÜNCELLENDİ ---
DATA_FILE = "sparc_raw_full_3389_points.csv"
LIST_FILE = "sparc_subset_125_list.txt"  # 125 olarak güncellendi
OUTPUT_CHAIN = "chain_results.txt"
OUTPUT_SUMMARY = "analysis_summary.txt"

def load_data():
    if not os.path.exists(DATA_FILE):
        print(f"Error: Data file {DATA_FILE} not found.")
        sys.exit(1)
    
    df = pd.read_csv(DATA_FILE)
    
    if os.path.exists(LIST_FILE):
        with open(LIST_FILE, 'r') as f:
            valid_galaxies = sorted(list(set([line.strip() for line in f if line.strip()])))
        
        df = df[df['Galaxy'].isin(valid_galaxies)].copy()
        print(f"--> Analysis set: {len(valid_galaxies)} galaxies loaded from list.")
    else:
        print(f"Error: {LIST_FILE} not found!")
        sys.exit(1)
    
    df = df[df['Radius'] > 0].copy()
    R_kpc = df['Radius'].values
    df['g_obs'] = (df['Vobs'].values**2 / R_kpc) * CONV_FACTOR
    df['term_gas'] = (df['Vgas'].values**2 / R_kpc) * CONV_FACTOR
    df['term_disk'] = (df['Vdisk'].values**2 / R_kpc) * CONV_FACTOR
    df['term_bul'] = (df['Vbul'].values**2 / R_kpc) * CONV_FACTOR
    df = df.replace([np.inf, -np.inf], np.nan).dropna()
    return df

def theory_leakage(g_bar, n, log_a5):
    a5 = 10**log_a5
    g_bar = np.maximum(g_bar, 1e-20)
    y = g_bar / a5
    exponent = -np.power(np.sqrt(y), 2.0/n)
    denominator = 1.0 - np.exp(exponent)
    denominator = np.maximum(denominator, 1e-10)
    return g_bar / denominator

def log_likelihood(theta, data):
    n, log_a5, sigma_int, f_out = theta
    g_bar = data['term_gas'] + (UPSILON_DISK * data['term_disk']) + (UPSILON_BULGE * data['term_bul'])
    g_pred = theory_leakage(g_bar, n, log_a5)
    
    mask = (data['g_obs'] > 1e-20) & (g_pred > 1e-20)
    if np.sum(mask) < 5: return -1e10
    
    res = np.log10(data['g_obs'][mask]) - np.log10(g_pred[mask])
    sigma_bad = 0.5
    ln_p_good = -0.5 * (res**2 / sigma_int**2 + np.log(2 * np.pi * sigma_int**2))
    ln_p_bad  = -0.5 * (res**2 / sigma_bad**2 + np.log(2 * np.pi * sigma_bad**2))
    log_L = np.logaddexp(np.log(1 - f_out) + ln_p_good, np.log(f_out) + ln_p_bad)
    return np.sum(log_L)

def prior_transform(u):
    return (1.0 + u[0]*2.0, -13.0 + u[1]*5.0, 0.05 + u[2]*0.2, u[3]*0.1)

if __name__ == "__main__":
    print("--> Loading data...")
    df = load_data()
    print(f"--> Starting Dynesty run on {len(df)} data points...")
    dsampler = dynesty.DynamicNestedSampler(
        lambda t: log_likelihood(t, df),
        prior_transform,
        ndim=4,
        bound='multi',
        sample='rwalk',
        nlive=500
    )
    dsampler.run_nested(dlogz_init=0.05, print_progress=True)
    
    print("--> Saving results...")
    res = dsampler.results
    data_out = np.column_stack((res.samples, res.logl))
    np.savetxt(OUTPUT_CHAIN, data_out, header="n log_a5 sigma f_out logL", comments="")
    
    mean, cov = dyfunc.mean_and_cov(res.samples, np.exp(res.logwt - res.logz[-1]))
    std = np.sqrt(np.diag(cov))
    
    with open(OUTPUT_SUMMARY, 'w') as f:
        f.write(f"# Analysis Summary\n")
        f.write(f"Galaxies: {df['Galaxy'].nunique()}\n")
        f.write(f"Data Points: {len(df)}\n")
        f.write(f"Log Evidence: {res.logz[-1]:.2f}\n")
        f.write(f"n: {mean[0]:.4f} +/- {std[0]:.4f}\n")
        f.write(f"log_a5: {mean[1]:.4f} +/- {std[1]:.4f}\n")
        f.write(f"sigma: {mean[2]:.4f} +/- {std[2]:.4f}\n")
    print("--> Done.")
