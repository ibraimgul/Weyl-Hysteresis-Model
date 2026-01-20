"""
Non-Local Gravity / Dark Dimension Leakage Analysis
===================================================
Author: Ibrahim Gul
Status: Final Production Version
Description:
  This script performs a Bayesian analysis on galaxy rotation curves to constrain
  the leakage parameters (n, a5) using the SPARC 'Golden Sample'.
  It implements the Spectral Propagator Ansatz derived from 5D Einstein-Dilaton theory.

Outputs:
  - chain_results.txt: Raw posterior chains for plotting.
  - analysis_summary.txt: Statistical summary (Mean, Error, Bayes Factor).
"""

import numpy as np
import pandas as pd
import sys
import os
import time

# Kütüphane Kontrolü
try:
    import dynesty
    from dynesty import utils as dyfunc
except ImportError:
    print("CRITICAL ERROR: 'dynesty' library not found.")
    print("Please install it via: pip install dynesty")
    sys.exit(1)

# --- 1. FİZİKSEL SABİTLER VE AYARLAR ---
# Birim Dönüşümü: (km/s)^2 / kpc -> m/s^2
KPC_TO_M = 3.086e19
KM_TO_M = 1000.0
CONV_FACTOR = (KM_TO_M**2) / KPC_TO_M

# Kütle-Işık Oranları (Stellar Mass-to-Light Ratios)
# SPARC "Standard" değerleri (PopIII/Kroupa IMF ile tutarlı)
UPSILON_DISK = 0.5
UPSILON_BULGE = 0.7

# Dosya Yolları
DATA_FILE = "sparc_raw_full_3389_points.csv"
LIST_FILE = "sparc_subset_118_list.txt"
OUTPUT_CHAIN = "chain_results.txt"
OUTPUT_SUMMARY = "analysis_summary.txt"

# --- 2. VERİ YÜKLEME VE İŞLEME ---
def load_data():
    print(f"--> Loading raw data from {DATA_FILE}...")
    
    if not os.path.exists(DATA_FILE):
        print(f"ERROR: {DATA_FILE} not found!")
        sys.exit(1)
    
    df_full = pd.read_csv(DATA_FILE)
    
    # Filtreleme Listesi Kontrolü
    if os.path.exists(LIST_FILE):
        with open(LIST_FILE, 'r') as f:
            # Tekrarlayan isimleri (duplicates) temizle ve sırala
            raw_list = [line.strip() for line in f if line.strip()]
            valid_galaxies = sorted(list(set(raw_list)))
        
        # Sadece listedeki galaksileri seç
        df = df_full[df_full['Galaxy'].isin(valid_galaxies)].copy()
        
        print(f"--> Filter Applied: {len(valid_galaxies)} unique IDs found in list.")
        print(f"--> Final Selection: Using {df['Galaxy'].nunique()} galaxies found in dataset.")
    else:
        print("--> WARNING: Filter list not found. Using FULL dataset (Caution!).")
        df = df_full.copy()

    # İvme Hesaplamaları (Fiziksel Birim: m/s^2)
    # R <= 0 olan hatalı noktaları temizle
    df = df[df['Radius'] > 0].copy()
    
    R_kpc = df['Radius'].values
    
    # Gözlenen İvme: g_obs = V_obs^2 / R
    df['g_obs'] = (df['Vobs'].values**2 / R_kpc) * CONV_FACTOR
    
    # Baryonik Bileşenler: g_x = V_x^2 / R
    # (Bileşenleri ayrı ayrı saklıyoruz, likelihood'da birleştireceğiz)
    df['term_gas'] = (df['Vgas'].values**2 / R_kpc) * CONV_FACTOR
    df['term_disk'] = (df['Vdisk'].values**2 / R_kpc) * CONV_FACTOR
    df['term_bul'] = (df['Vbul'].values**2 / R_kpc) * CONV_FACTOR
    
    # NaN veya Sonsuz değerleri temizle
    df = df.replace([np.inf, -np.inf], np.nan).dropna()
    
    print(f"--> Data Ready: {len(df)} total radial points loaded.")
    return df

# --- 3. TEORİK MODEL (DARK DIMENSION LEAKAGE) ---
def theory_leakage(g_bar, n, log_a5):
    """
    Hesaplar: g_eff = g_bar * nu(g_bar/a5)
    Model: Spectral Propagator Ansatz (Makale Eq. 2)
    Formül: nu(y) = 1 / (1 - exp(-sqrt(y)^(2/n)))
    """
    a5 = 10**log_a5
    
    # Sayısal kararlılık (0'a bölmeyi önle)
    g_bar = np.maximum(g_bar, 1e-20)
    y = g_bar / a5
    
    # Bessel Modlarından Türetilen Geçiş Fonksiyonu
    # n = 2*nu (Order of Bessel function related)
    exponent = -np.power(np.sqrt(y), 2.0/n)
    
    # Overflow koruması
    denominator = 1.0 - np.exp(exponent)
    denominator = np.maximum(denominator, 1e-10) # Asla 0 olmasın
    
    nu = 1.0 / denominator
    return g_bar * nu

# --- 4. BAYESIAN LIKELIHOOD & PRIOR ---
def log_likelihood(theta, data):
    # theta: [n, log_a5, sigma_int, f_out]
    n, log_a5, sigma_int, f_out = theta
    
    # 1. Baryonik İvme (g_bar)
    g_bar = data['term_gas'] + (UPSILON_DISK * data['term_disk']) + (UPSILON_BULGE * data['term_bul'])
    
    # 2. Model Tahmini (g_pred)
    g_pred = theory_leakage(g_bar, n, log_a5)
    
    # 3. Residuals (Logaritmik Uzayda)
    mask = (data['g_obs'] > 1e-20) & (g_pred > 1e-20)
    if np.sum(mask) < 5: return -1e10 # Veri çok kötüyse öldür
    
    obs = data['g_obs'][mask]
    pred = g_pred[mask]
    
    res = np.log10(obs) - np.log10(pred)
    
    # 4. Robust Likelihood (Mixture Model)
    # Makaledeki %10 outlier toleransını uygular.
    sigma_bad = 0.5 
    
    ln_p_good = -0.5 * (res**2 / sigma_int**2 + np.log(2 * np.pi * sigma_int**2))
    ln_p_bad  = -0.5 * (res**2 / sigma_bad**2 + np.log(2 * np.pi * sigma_bad**2))
    
    log_L = np.logaddexp(np.log(1 - f_out) + ln_p_good, np.log(f_out) + ln_p_bad)
    
    return np.sum(log_L)

def prior_transform(u):
    """
    Unit Cube [0,1] -> Physical Parameters
    Aralıklar makale sonuçlarını (n~1.46, log_a5~-10.85) kapsayacak şekilde ayarlandı.
    """
    n = 1.0 + u[0] * 2.0        # n: [1.0, 3.0]
    log_a5 = -13.0 + u[1]*5.0   # log_a5: [-13, -8]
    sigma = 0.05 + u[2]*0.2     # sigma: [0.05, 0.25]
    f_out = u[3] * 0.1          # f_out: [0.0, 0.1]
    return n, log_a5, sigma, f_out

# --- 5. ANA ÇALIŞMA BLOĞU ---
if __name__ == "__main__":
    start_time = time.time()
    
    # 1. Veri Hazırlığı
    df = load_data()
    
    # 2. Dynesty Kurulumu
    print("\n--> Initializing Dynamic Nested Sampling (Dynesty)...")
    dsampler = dynesty.DynamicNestedSampler(
        lambda t: log_likelihood(t, df),
        prior_transform,
        ndim=4,
        bound='multi',
        sample='rwalk',
        nlive=500  # Makale kalitesinde örnekleme
    )
    
    # 3. Analizi Başlat
    print("--> Running Sampling (This may take several minutes)...")
    dsampler.run_nested(dlogz_init=0.05, print_progress=True)
    
    # 4. Sonuçları Kaydet
    res = dsampler.results
    samples = res.samples
    
    # Dosyaya yaz (Grafik çizimi için ham zincir)
    save_data = np.column_stack((samples, res.logl))
    np.savetxt(OUTPUT_CHAIN, save_data, header="n log_a5 sigma f_out logL", comments="")
    
    print(f"\n--> ANALYSIS COMPLETE. Results saved to {OUTPUT_CHAIN}")
