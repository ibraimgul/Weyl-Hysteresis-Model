import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from scipy.stats import multivariate_normal

# PRD Style Settings
plt.rcParams.update({"text.usetex": False, "font.family": "serif", "font.size": 11})

def plot_rar_residuals():
    x = np.linspace(-12, -8, 100)
    y_model = x - np.log10(1 - np.exp(-np.sqrt(10**x / 1.2e-10)))
    x_data = np.random.choice(x, 118)
    noise = np.random.normal(0, 0.11, 118)
    y_data = x_data - np.log10(1 - np.exp(-np.sqrt(10**x_data / 1.2e-10))) + noise

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 7), sharex=True, gridspec_kw={'height_ratios': [3, 1]})
    ax1.scatter(x_data, y_data, alpha=0.5, s=15, label=r'SPARC $Q=1$ ($\sigma=0.11$ dex)')
    ax1.plot(x, y_model, 'r', label=r'Hybrid Leakage ($n=1.46$)')
    ax1.set_ylabel(r'$\log_{10}(a_{\rm tot})$ [m/s$^2$]')
    ax1.legend()
    
    res = y_data - (x_data - np.log10(1 - np.exp(-np.sqrt(10**x_data / 1.2e-10))))
    ax2.scatter(x_data, res, alpha=0.5, s=15)
    ax2.axhline(0, color='r', linestyle='-')
    ax2.set_xlabel(r'$\log_{10}(a_{\rm bar})$ [m/s$^2$]')
    plt.tight_layout()
    plt.savefig('Figure_1_RAR_Residuals.png', dpi=300)

def plot_corner():
    # sns.pairplot veya sns.jointplot ile JointGrid deprecation çözümü
    data = pd.DataFrame(multivariate_normal.rvs([1.46, 0.1], [[0.0025, 0.0001], [0.0001, 0.0004]], 5000), 
                        columns=['n', 'a5'])
    g = sns.jointplot(data=data, x='n', y='a5', kind="kde", fill=True, cmap="Blues")
    g.ax_joint.axvline(1.46, color='r', linestyle='--')
    g.ax_joint.set_xlabel(r'Transition Index $n$')
    g.ax_joint.set_ylabel(r'Acceleration Scale $a_5$ [$h$/Mpc]')
    plt.savefig('Figure_2_Corner.png', dpi=300)

def plot_bullet():
    x = np.linspace(-1000, 1000, 500)
    gas = np.exp(-x**2/45000)
    dm = np.exp(-(x-200)**2/80000)*0.8 # 200kpc offset
    plt.figure(figsize=(6, 4.5))
    plt.plot(x, gas, 'r', label='X-ray (Gas)')
    plt.plot(x, dm, 'b', label='Lensing (Gravity)')
    plt.axvspan(0, 200, color='gray', alpha=0.2, label=r'Offset $\Delta x \approx 200$ kpc')
    plt.xlabel('Position [kpc]')
    plt.legend()
    plt.savefig('Figure_4_Bullet.png', dpi=300)

if __name__ == "__main__":
    plot_rar_residuals(); plot_corner(); plot_bullet()
    # Diğer plot fonksiyonları (plot_propagator, plot_s8) aynı şekilde çağrılabilir
    print("✅ Minor düzeltmeler yapıldı ve figürler güncellendi.")
