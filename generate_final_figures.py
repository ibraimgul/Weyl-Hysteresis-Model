import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import corner
import os
import sys

# --- AYARLAR ---
# Makale kalitesinde grafik ayarları
plt.rcParams.update({
    'font.size': 14,
    'axes.linewidth': 1.5,
    'xtick.major.width': 1.5,
    'ytick.major.width': 1.5,
    'font.family': 'serif'
})

KPC_TO_M = 3.086e19
KM_TO_M = 1000.0
CONV_FACTOR = (KM_TO_M**2) / KPC_TO_M

# --- VERİ YÜKLEME ---
def load_golden_sample():
    print("Loading data for plotting...")
    # 1. Ham Veriyi Oku
    if not os.path.exists("sparc_raw_full_3389_points.csv"):
        print("Error: sparc_raw_full_3389_points.csv not found.")
        sys.exit(1)
    df_full = pd.read_csv("sparc_raw_full_3389_points.csv")
    
    # 2. Listeyi Oku ve Filtrele
    if not os.path.exists("sparc_subset_118_list.txt"):
        print("Error: sparc_subset_118_list.txt not found.")
        sys.exit(1)
        
    with open("sparc_subset_118_list.txt", 'r') as f:
        valid_galaxies = [line.strip() for line in f if line.strip()]
    
    df = df_full[df_full['Galaxy'].isin(valid_galaxies)].copy()
    
    # 3. İvme Hesapla
    # Sadece R > 0 olanları al
    df = df[df['Radius'] > 0].copy()
    
    # g_obs
    df['g_obs'] = (df['Vobs']**2 / df['Radius']) * CONV_FACTOR
    
    # g_bar (Sabit M/L varsayımı: Disk=0.5, Bulge=0.7)
    V_bar_sq = df['Vgas']**2 + (0.5 * df['Vdisk']**2) + (0.7 * df['Vbul']**2)
    df['g_bar'] = (V_bar_sq / df['Radius']) * CONV_FACTOR
    
    return df

# --- TEORİK FONKSİYON ---
def theory_curve(g_bar, n, log_a5):
    a5 = 10**log_a5
    # nu(y) = 1 / (1 - exp(-sqrt(y)^(2/n)))
    y = g_bar / a5
    nu = 1.0 / (1.0 - np.exp(-np.sqrt(y)**(2.0/n)))
    return g_bar * nu

# --- FİGÜR 1: RAR + RESIDUALS ---
def plot_rar(df, n_best, log_a5_best):
    print("Generating Figure 1: RAR Relation...")
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 10), gridspec_kw={'height_ratios': [3, 1]}, sharex=True)
    
    # Veri Noktaları (Log Scale)
    x = np.log10(df['g_bar'])
    y = np.log10(df['g_obs'])
    
    # Hata barları (basit görselleştirme için şeffaf)
    ax1.scatter(x, y, alpha=0.15, c='gray', s=10, label='SPARC Data (118 Galaxies)')
    
    # Teorik Eğri
    g_bar_grid = np.logspace(-12, -8, 100)
    g_pred = theory_curve(g_bar_grid, n_best, log_a5_best)
    
    ax1.plot(np.log10(g_bar_grid), np.log10(g_pred), 'r-', lw=2.5, label=f'Leakage Model (n={n_best:.2f})')
    
    # 1:1 Çizgisi (Newton)
    ax1.plot([-12, -8], [-12, -8], 'k--', lw=1.5, label='Newtonian (1:1)')
    
    ax1.set_ylabel(r'$\log(g_{obs})$ [m/s$^2$]')
    ax1.legend(loc='upper left')
    ax1.set_title("Radial Acceleration Relation (Golden Sample)")
    ax1.grid(True, which='both', linestyle='--', alpha=0.3)

    # Residuals (Observed - Theory)
    # Veri noktaları için teori tahmini
    g_theory_points = theory_curve(df['g_bar'].values, n_best, log_a5_best)
    res = np.log10(df['g_obs'].values) - np.log10(g_theory_points)
    
    ax2.scatter(x, res, alpha=0.15, c='gray', s=10)
    ax2.axhline(0, color='r', linestyle='-', lw=2)
    ax2.set_ylabel('Residuals (dex)')
    ax2.set_xlabel(r'$\log(g_{bar})$ [m/s$^2$]')
    ax2.set_ylim(-0.5, 0.5)
    ax2.grid(True, which='both', linestyle='--', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig("Figure_1_RAR_Residuals.png", dpi=300)
    print("Saved: Figure_1_RAR_Residuals.png")

# --- FİGÜR 2: CORNER PLOT ---
def plot_corner():
    print("Generating Figure 2: Corner Plot...")
    
    if not os.path.exists("chain_results.txt"):
        print("Warning: chain_results.txt not found. Skipping Corner Plot.")
        return None, None

    # Dosyayı oku (n, log_a5, logL)
    try:
        data = np.loadtxt("chain_results.txt", skiprows=1)
        samples = data[:, :2] # Sadece n ve log_a5 al
        
        # En iyi değerleri (median) bul
        n_med = np.median(samples[:, 0])
        a5_med = np.median(samples[:, 1])
        
        # Corner Plot
        fig = corner.corner(
            samples, 
            labels=[r"$n$ (Index)", r"$\log a_5$ (Scale)"],
            quantiles=[0.16, 0.5, 0.84],
            show_titles=True,
            title_kwargs={"fontsize": 12},
            color="#1f77b4"
        )
        fig.savefig("Figure_2_Corner.png", dpi=300)
        print("Saved: Figure_2_Corner.png")
        return n_med, a5_med
        
    except Exception as e:
        print(f"Corner plot error: {e}")
        return 1.46, -10.85 # Hata olursa varsayılan dön

# --- MAIN ---
if __name__ == "__main__":
    # 1. Corner Plot çiz ve en iyi n değerini al
    n_best, log_a5_best = plot_corner()
    
    if n_best is None:
        # Eğer chain dosyası yoksa varsayılan değerleri kullan
        n_best = 1.46
        log_a5_best = -10.85
    
    # 2. Veriyi yükle ve RAR grafiğini çiz
    df = load_golden_sample()
    plot_rar(df, n_best, log_a5_best)
    
    print("\nAll figures generated successfully!")
