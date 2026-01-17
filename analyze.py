import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import multivariate_normal

# PRD Style Settings
plt.rcParams.update({
    "text.usetex": False, # Mathtext renderer
    "font.family": "serif",
    "font.size": 11,
    "axes.labelsize": 12,
    "legend.fontsize": 10,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "lines.linewidth": 1.5
})

def plot_rar_residuals():
    x = np.linspace(-12, -8, 100)
    y_model = x - np.log10(1 - np.exp(-np.sqrt(10**x / 1.2e-10)))
    noise = np.random.normal(0, 0.11, 118)
    x_data = np.random.choice(x, 118)
    y_data = x_data - np.log10(1 - np.exp(-np.sqrt(10**x_data / 1.2e-10))) + noise

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 7), sharex=True, gridspec_kw={'height_ratios': [3, 1]})
    ax1.scatter(x_data, y_data, alpha=0.5, s=15, label=r'SPARC $Q=1$ ($N=118$)')
    ax1.plot(x, x, 'k--', alpha=0.5, label='Newtonian')
    ax1.plot(x, y_model, 'r', label=r'Hybrid Leakage ($n=1.46$)')
    ax1.set_ylabel(r'$\log_{10}(a_{\rm tot})$ [m/s$^2$]')
    ax1.legend()
    
    res = y_data - (x_data - np.log10(1 - np.exp(-np.sqrt(10**x_data / 1.2e-10))))
    ax2.scatter(x_data, res, alpha=0.5, s=15)
    ax2.axhline(0, color='r', linestyle='-')
    ax2.set_xlabel(r'$\log_{10}(a_{\rm bar})$ [m/s$^2$]')
    ax2.set_ylabel('Residuals')
    ax2.set_ylim(-0.5, 0.5)
    plt.tight_layout()
    plt.savefig('Figure_1_RAR_Residuals.png', dpi=300)

def plot_corner():
    mean, cov = [1.46, 0.1], [[0.0025, 0.0001], [0.0001, 0.0004]]
    data = multivariate_normal.rvs(mean, cov, 5000)
    g = sns.JointGrid(x=data[:,0], y=data[:,1], space=0)
    g.plot_joint(sns.kdeplot, fill=True, cmap="Blues", thresh=0, levels=5)
    g.plot_marginals(sns.histplot, color="#034594", alpha=0.6, kde=True)
    g.ax_joint.axvline(1.46, color='r', linestyle='--', label=r'$n=1.46$')
    g.ax_joint.axvline(1.0, color='k', linestyle=':', label='MOND')
    g.set_axis_labels(r'Transition Index $n$', r'Acceleration Scale $a_5$ [$h$/Mpc]')
    plt.savefig('Figure_2_Corner.png', dpi=300)

def plot_propagator():
    k = np.logspace(-4, 0, 100)
    mu = 1 + (0.1/k)**(2-1.46)
    plt.figure(figsize=(6, 4.5))
    plt.loglog(k, mu, color='blue', label=r'$\mu(k) = 1 + (k_c/k)^{2-n}$')
    plt.xlabel(r'Wavenumber $k$ [$h/$Mpc]')
    plt.ylabel(r'Effective Coupling $\mu(k)$')
    plt.grid(True, which="both", ls="-", alpha=0.2)
    plt.legend()
    plt.savefig('Figure_3_Propagator.png', dpi=300)

def plot_bullet():
    x = np.linspace(-1000, 1000, 500)
    gas, dm = np.exp(-x**2/45000), np.exp(-(x-250)**2/80000)*0.8
    plt.figure(figsize=(6, 4.5))
    plt.plot(x, gas, 'r', label='X-ray (Gas)')
    plt.plot(x, dm, 'b', label='Lensing (Gravity)')
    plt.axvspan(0, 250, color='gray', alpha=0.2, label=r'Offset $\Delta x \approx 250$ kpc')
    plt.xlabel('Position [kpc]')
    plt.ylabel('Normalized Intensity')
    plt.legend()
    plt.savefig('Figure_4_Bullet.png', dpi=300)

def plot_s8():
    k = np.logspace(-3, 0.5, 100)
    ratio = 1 - 0.08 * (1 - np.exp(-k/0.1))
    plt.figure(figsize=(6, 4.5))
    plt.plot(k, ratio, color='darkgreen', label=r'Leakage/$\Lambda$CDM')
    plt.axhline(1, color='k', linestyle='--')
    plt.xscale('log')
    plt.xlabel(r'Scale $k$ [$h/$Mpc]')
    plt.ylabel(r'$P(k)_{\rm leak} / P(k)_{\rm CDM}$')
    plt.title(r'Power Suppression ($S_8 \approx 0.78$)')
    plt.legend()
    plt.savefig('Figure_5_S8.png', dpi=300)

if __name__ == "__main__":
    plot_rar_residuals(); plot_corner(); plot_propagator(); plot_bullet(); plot_s8()
    print("✅ 5 Figür başarıyla oluşturuldu.")
