"""
Non-Local Gravity / Dark Dimension Leakage Analysis
===================================================
Author: Ibrahim Gul
Description:
  This script performs a Bayesian analysis on galaxy rotation curves to constrain
  the leakage parameters (n, a5). It reads raw data points from the SPARC database
  subset (Golden Sample) and compares observed accelerations with theoretical predictions.

Data Source:
  - sparc_raw_full_3389_points.csv (Full raw rotation curves)
  - golden_sample_list.txt (List of 118 high-quality galaxies)

Method:
  - Dynamic Nested Sampling (Dynesty)
  - Robust Likelihood with Outlier Modeling
"""

import numpy as np
import pandas as pd
import sys
import os

# Dynesty Kontrolü
try:
    import dynesty
    from dynesty import utils as dyfunc
except ImportError:
    print("Error: 'dynesty' library not found. Install via: pip install dynesty")
    sys.exit(1)

# --- FİZİKSEL SABİTLER ---
# Birim Dönüşümü: (km/s)^2 / kpc -> m/s^2
KPC_TO_M = 3.086e19
KM_TO_M = 1000.0
CONV_FACTOR = (KM_TO_M**2) / KPC_TO_M

# Varsayılan Kütle-Işık Oranları (PopIII/Kroupa IMF referanslı)
UPSILON_DISK = 0.5
UPSILON_BULGE = 0.7

# --- 1. VERİ YÜKLEME VE FİLTRELEME ---
def load_data(csv_path="sparc_raw_full_3389_points.csv", list_path="sparc_subset_118_list.txt"):
    print(f"Loading raw data from {csv_path}...")
    
    # 1. Ana Ham Veriyi Oku
    if not os.path.exists(csv_path):
        print(f"CRITICAL ERROR: {csv_path} not found!")
        sys.exit(1)
    
    df_full = pd.read_csv(csv_path)
    
    # 2. 118 Galaksilik 'Altın Liste'yi Uygula
    # Eğer liste dosyası varsa onu kullan, yoksa hepsini kullan (ama uyar)
    if os.path.exists(list_path):
        with open(list_path, 'r') as f:
            valid_galaxies = [line.strip() for line in f if line.strip()]
        
        df = df_full[df_full['Galaxy'].isin(valid_galaxies)].copy()
        print(f"Applied Filter: Selected {len(valid_galaxies)} galaxies (Golden Sample).")
        print(f"Total Data Points used for Analysis: {len(df)}")
    else:
        print("WARNING: 'sparc_subset_118_list.txt' not found. Using ALL 175 galaxies (Results may be noisy).")
        df = df_full.copy()

    # 3. İvme Hesaplamaları (Ham Hızlardan)
    # g = V^2 / R formülü uygulanıyor
    R_kpc = df['Radius'].values
    valid_mask = R_kpc > 0  # Yarıçapı 0 olanları at
    df = df[valid_mask]
    R_kpc = df['Radius'].values

    # Gözlenen İvme
    g_obs = (df['Vobs'].values**2 / R_kpc) * CONV_FACTOR
    
    # Baryonik İvme Bileşenleri (Disk, Gaz, Bulge)
    # V_bar^2 = V_gas^2 + V_disk^2 * Y_disk + V_bul^2 * Y_bul
    V_gas_sq = df['Vgas'].values**2
    V_disk_sq = df['Vdisk'].values**2
    V_bul_sq = df['Vbul'].values**2
    
    # Bileşenleri dataframe'e ekle (Likelihood içinde M/L ile çarpacağız)
    df['g_obs'] = g_obs
    df['term_gas'] = (V_gas_sq / R_kpc) * CONV_FACTOR
    df['term_disk'] = (V_disk_sq / R_kpc) * CONV_FACTOR
    df['term_bul'] = (V_bul_sq / R_kpc) * CONV_FACTOR
    
    return df

# --- 2. TEORİK MODEL (Sızıntı) ---
def theory_leakage(g_bar, n, log_a5):
    """
    Hesaplar: g_eff = g_bar * nu(g_bar/a5, n)
    """
    a5 = 10**log_a5
    # Sızıntı Fonksiyonu (Spectral Propagator Ansatz)
    # nu(y) = 1 / (1 - exp(-sqrt(y)^(2/n)))
    # y = g_bar / a5
    
    # Sayısal kararlılık için çok küçük g_bar değerlerini koru
    g_bar = np.maximum(g_bar, 1e-18)
    y = g_bar / a5
    
    # Üs hesabı
    exponent = -np.power(np.sqrt(y), 2.0/n)
    denominator = 1.0 - np.exp(exponent)
    
    nu = 1.0 / denominator
    return g_bar * nu

# --- 3. LOG-LIKELIHOOD ---
def log_likelihood(theta, data):
    # theta: [n, log_a5, intrinsic_scatter, f_outlier]
    n, log_a5, sigma_int, f_out = theta
    
    # Baryonik İvme Hesabı (Sabit M/L varsayımı ile)
    # İstersen M/L'yi de parametre yapabilirsin ama şimdilik standart tutuyoruz.
    g_bar = data['term_gas'] + (UPSILON_DISK * data['term_disk']) + (UPSILON_BULGE * data['term_bul'])
    
    # Teorik Tahmin
    g_pred = theory_leakage(g_bar, n, log_a5)
    
    # Logaritmik Farklar (Residuals)
    g_obs = data['g_obs']
    
    # Sadece pozitif ivmelerde log alınabilir
    mask = (g_obs > 1e-15) & (g_pred > 1e-15)
    res = np.log10(g_obs[mask]) - np.log10(g_pred[mask])
    
    # Mixture Model (Gürültü yönetimi)
    # Model 1: İyi veri (Sigma_total^2 = Sigma_obs^2 + Sigma_int^2)
    # Burada basitleştirme için global sigma kullanıyoruz (makalendeki gibi)
    sigma_total = sigma_int 
    
    # Model 2: Outlier dağılımı (Daha geniş sigma)
    sigma_bad = 0.5 # 0.5 dex saçılma (çok kötü veri)
    
    # Log-Likelihood Calculation
    ln_p_good = -0.5 * (res**2 / sigma_total**2 + np.log(2 * np.pi * sigma_total**2))
    ln_p_bad  = -0.5 * (res**2 / sigma_bad**2   + np.log(2 * np.pi * sigma_bad**2))
    
    # Toplam Olasılık: (1-f)*Good + f*Bad
    log_L = np.logaddexp(np.log(1 - f_out) + ln_p_good, np.log(f_out) + ln_p_bad)
    
    return np.sum(log_L)

# --- 4. PRIOR DÖNÜŞÜMÜ ---
def prior_transform(u):
    """
    Unit cube [0,1] -> Physical Parameters
    """
    n = 1.0 + u[0] * 2.0       # n: [1.0, 3.0]
    log_a5 = -13.0 + u[1]*5.0  # log_a5: [-13, -8]
    sigma = 0.05 + u[2]*0.2    # sigma: [0.05, 0.25]
    f_out = u[3] * 0.1         # f_out: [0.0, 0.1] (Max %10 outlier kabulü)
    return n, log_a5, sigma, f_out

# --- MAIN ÇALIŞTIRMA ---
if __name__ == "__main__":
    # Veriyi Yükle
    df = load_data()
    
    # Sampler Ayarları
    print("\nInitializing Dynamic Nested Sampling...")
    dsampler = dynesty.DynamicNestedSampler(
        lambda t: log_likelihood(t, df),
        prior_transform,
        ndim=4,
        bound='multi',
        sample='rwalk',
        nlive=500
    )
    
    # Çalıştır
    print("Running Analysis (This may take time)...")
    dsampler.run_nested(dlogz_init=0.05)
    
    # Sonuçları Kaydet
    res = dsampler.results
    
    # 1. chain_results.txt (Ham zincir - Grafikler için)
    # Format: n, log_a5, log_likelihood (Ağırlıklı örneklerden)
    print("Saving chain results...")
    samples = res.samples  # param samples
    logl = res.logl        # log likelihoods
    
    # Dosyaya yaz (n, a5, logL sütunları)
    # generate_final_figures.py bu formatı bekliyor
    output_data = np.column_stack((samples[:, 0], samples[:, 1], logl))
    np.savetxt("chain_results.txt", output_data, header="n_index log10_a5 log_likelihood", comments="")
    
    print("\n--- ANALYSIS COMPLETE ---")
    print(f"Results saved to chain_results.txt")
