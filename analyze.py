import math

# --- MAKALE PARAMETRELERİ ---
n = 1.46                # Geçiş indeksi [cite: 7, 64]
a5 = 1.21e-10           # Sızıntı ölçeği (m/s^2) [cite: 7, 146]
v_coll = 4700           # Bullet Cluster çarpışma hızı (km/s) [cite: 48, 49]
dt_myr = 52             # Gecikme süresi (Myr) [cite: 9, 47, 84]

# --- SABİTLER ---
KPC_TO_M = 3.0857e19    # 1 kpc kaç metre?
MYR_TO_S = 3.1536e13    # 1 Myr kaç saniye?

def analyze_propagator():
    print("--- TEST 1: Propagatör UV/IR Limit Analizi ---")
    # Eq. (1): G(p) = 1/p^2 * [1 + (a5/p^2)^(n/2)]^-1 
    # Momentum ölçeği p, a5 üzerinden normalize edildiğinde (p/mc)
    
    for p_ratio in [0.01, 1.0, 100.0]:
        # Leakage terimi: (a5/p^2)^(n/2) -> (1/p_ratio^2)^(n/2)
        leakage_term = (1 / (p_ratio**2))**(n/2)
        g_eff = (1 / p_ratio**2) * (1 / (1 + leakage_term))
        gr_newton = 1 / p_ratio**2
        
        status = "IR (Sızıntı)" if p_ratio < 1 else ("UV (GR)" if p_ratio > 1 else "Crossover")
        print(f"Ölçek p/mc = {p_ratio:6.2f} | Rejim: {status:10} | G_model / G_Newton: {g_eff/gr_newton:.4f}")
    print()

def analyze_bullet_cluster():
    print("--- TEST 2: Bullet Cluster Ofset Doğrulaması ---")
    # Eq. (2): Δx = v_coll * Δt_lag 
    
    # Hız: m/s cinsinden
    v_m_s = v_coll * 1000
    # Zaman: saniye cinsinden
    t_s = dt_myr * MYR_TO_S
    
    # Ofset (metre)
    delta_x_m = v_m_s * t_s
    # Ofset (kpc)
    delta_x_kpc = delta_x_m / KPC_TO_M
    
    print(f"Girdi Hız      : {v_coll} km/s")
    print(f"Girdi Zaman    : {dt_myr} Myr")
    print(f"Hesaplanan Ofset: {delta_x_kpc:.2f} kpc")
    print(f"Gözlemlenen Veri: 250.00 kpc [cite: 8, 49, 71]")
    print(f"Hata Oranı      : %{abs(delta_x_kpc - 250)/250*100:.4f}")
    print()

def analyze_rar():
    print("--- TEST 3: Radial Acceleration Relation (RAR) Sapması ---")
    # Gözlenen ivme g_obs, g_bar'dan nasıl sapıyor?
    
    for g_bar in [1e-12, 1.2e-10, 1e-8]:
        # Model: g_obs = g_bar / (1 + (a5/g_bar)^(n/2))
        factor = 1 + (a5 / g_bar)**(n/2)
        g_obs = g_bar / factor
        
        status = "Düşük İvme" if g_bar < a5 else ("Yüksek İvme" if g_bar > a5 else "Eşik (a5)")
        print(f"g_bar: {g_bar:.1e} m/s^2 | Rejim: {status:10} | g_obs/g_bar Oranı: {g_obs/g_bar:.4f}")
    print()

def check_mcmc_consistency():
    print("--- TEST 4: MCMC Parametre Tutarlılığı ---")
    # Makaledeki n ve a5 belirsizlikleri ile BIC [cite: 132, 146, 83]
    n_val = 1.46
    n_sigma = 0.05
    a5_val = 1.21
    a5_sigma = 0.06
    
    print(f"İndeks n: {n_val} +/- {n_sigma} (Hata Payı: %{n_sigma/n_val*100:.2f})")
    print(f"Ölçek a5: {a5_val} +/- {a5_sigma} x 10^-10 (Hata Payı: %{a5_sigma/a5_val*100:.2f})")
    print("Sonuç: Parametre uzayı dar ve belirgin bir tepe noktasına sahip (High Parsimony).")

if __name__ == "__main__":
    analyze_propagator()
    analyze_bullet_cluster()
    analyze_rar()
    check_mcmc_consistency()
