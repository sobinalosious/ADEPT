# tg_postproc_with_plots.py
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

# ------------------ Inputs ------------------
if len(sys.argv) < 3:
    raise SystemExit("Usage: python tg_postproc_with_plots.py <PID> <rough_Tg[K]>")

PID = sys.argv[1]
rough_Tg = float(sys.argv[2])  # rough Tg from command line

density_file = f"{PID}_tg_density.dat"
result_file  = f"{PID}_Tg_result.dat"

# Match your original knobs
n_fit_points_list = [5, 6, 7, 8]   # try multiple low/high fit sizes
skip_low  = 0
skip_high = 0

# Sampling layout (must match how you logged in LAMMPS)
n_eq   = 5000
n_prod = 5000
nstep  = n_eq + n_prod

# ------------------ Load & prep ------------------
data = np.loadtxt(density_file, comments="#")
temps_all = data[:, 1]
dens_all  = data[:, 2]

n_total   = len(temps_all)
n_Tpoints = n_total // nstep
if n_Tpoints < 2:
    raise ValueError("Not enough temperature points found for Tg fit.")

# Average production window at each temperature
avg_T, avg_rho = [], []
for i in range(n_Tpoints):
    a = i * nstep + n_eq
    b = a + n_prod
    if b > n_total:
        break
    avg_T.append(np.mean(temps_all[a:b]))
    avg_rho.append(np.mean(dens_all[a:b]))

avg_T   = np.array(avg_T)
avg_rho = np.array(avg_rho)

# Sort by temperature (ascending)
idx = np.argsort(avg_T)
avg_T = avg_T[idx]
avg_rho = avg_rho[idx]

# Optionally skip extremes
if skip_low + skip_high >= len(avg_T) - 2:
    raise ValueError("Too many points skipped for Tg fit.")
T_use   = avg_T[skip_low:len(avg_T)-skip_high]
rho_use = avg_rho[skip_low:len(avg_rho)-skip_high]

# ------------------ Tg estimation (same logic) ------------------
tg_vals = []
max_fit = len(T_use) // 2
fit_list = [n for n in n_fit_points_list if 1 <= n <= max_fit]
if not fit_list:
    raise ValueError("n_fit_points_list leaves no valid fit window sizes.")

for n in fit_list:
    low_idx  = np.arange(n)
    high_idx = np.arange(len(T_use) - n, len(T_use))

    sl_low, ic_low, *_ = stats.linregress(T_use[low_idx],  rho_use[low_idx])
    sl_hi,  ic_hi,  *_ = stats.linregress(T_use[high_idx], rho_use[high_idx])

    Tg = (ic_hi - ic_low) / (sl_low - sl_hi)
    tg_vals.append(Tg)
    rho_Tg = sl_low * Tg + ic_low

    plt.figure(figsize=(7, 5))
    plt.plot(T_use, rho_use, 'o-', label='Averaged Density')
    plt.plot(T_use[low_idx],  rho_use[low_idx],  's', label='Low-T Fit Region')
    plt.plot(T_use[high_idx], rho_use[high_idx], 's', label='High-T Fit Region')
    plt.plot(T_use, sl_low*T_use + ic_low,  '--', label='Low-T Fit')
    plt.plot(T_use, sl_hi*T_use  + ic_hi,   '--', label='High-T Fit')
    plt.plot(Tg, rho_Tg, 'ro', markersize=8, label=f'Tg = {Tg:.1f} K')
    plt.xlabel('Temperature (K)')
    plt.ylabel('Density (g/cm³)')
    plt.title(f'{PID}: Tg stepwise (n_fit_points={n})')
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'{PID}_tg_fit_{n}.png', dpi=300)
    plt.close()

tg_vals = np.array(tg_vals)
tg_mean = float(np.mean(tg_vals))
tg_std  = float(np.std(tg_vals, ddof=0))

# ------------------ Output only if within range ------------------
lower_bound = rough_Tg - 150
upper_bound = rough_Tg + 150

# Condition: all tg_vals must be inside range AND tg_mean inside range
all_within_range = np.all((tg_vals >= lower_bound) & (tg_vals <= upper_bound))
mean_within_range = lower_bound <= tg_mean <= upper_bound

if all_within_range and mean_within_range:
    with open(result_file, "w") as f:
        f.write(f"{tg_mean:.2f}\n")
        f.write(f"{tg_std:.2f}\n")
    print(f"Result written to {result_file} (Tg_mean={tg_mean:.2f})")
else:
    print(f"Skipped writing result: Tg values={tg_vals}, Tg_mean={tg_mean:.2f} not all within [{lower_bound},{upper_bound}]")
