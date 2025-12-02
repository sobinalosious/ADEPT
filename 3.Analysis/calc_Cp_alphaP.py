#!/usr/bin/env python3
import os
import sys
import numpy as np


def compute_cp_alpha_emd(filename, pid, n_windows=8):
    """
    EMD estimates from an NPT trajectory printed as:
        T[K], E_tot[kcal/mol], (pressure ignored), V[Å^3], rho[g/cm^3]

    Definitions (per window w):
      H(t)  = E(t) + P*V(t)                    [J], with P = 1 bar = 1.0e5 Pa
      m_w   = <rho*V>_w                        [kg]
      Cp_w  = Var(H)_w / (kB * <T>_w^2 * m_w)  [J/(kg·K)]
      αP_w  = Cov(V, H)_w / (kB * <T>_w^2 * <V>_w)   [1/K]
      αPL_w = αP_w / 3                         [1/K]

    We return only the mean values over windows:
      Cp_mean, alphaP_mean, alphaPL_mean
    """

    # -------- 1) Load data --------
    data = np.loadtxt(filename)
    if data.ndim == 1:
        data = data.reshape(-1, 5)

    T_all         = data[:, 0]   # K
    E_kcalmol_all = data[:, 1]   # kcal/mol
    V_ang3_all    = data[:, 3]   # Å^3
    rho_gcc_all   = data[:, 4]   # g/cm^3

    n_frames = len(data)
    if n_frames < n_windows:
        raise ValueError(f"Not enough frames ({n_frames}) to split into {n_windows} windows.")

    # -------- 2) Constants --------
    ANG3_TO_M3 = 1e-30
    KCAL_TO_J  = 4184.0
    NA         = 6.02214076e23
    KB         = 1.380649e-23     # J/K
    P_CONST_PA = 1.0e5            # 1 bar in Pa

    # -------- 3) Per-frame conversions --------
    V_m3_all     = V_ang3_all * ANG3_TO_M3
    rho_kgm3_all = rho_gcc_all * 1000.0
    m_kg_all     = rho_kgm3_all * V_m3_all            # kg
    E_J_all      = E_kcalmol_all * KCAL_TO_J / NA     # J
    H_J_all      = E_J_all + P_CONST_PA * V_m3_all    # J

    # -------- 4) Split into equal windows --------
    window_size = n_frames // n_windows
    n_used = window_size * n_windows
    if n_used != n_frames:
        print(f"[WARN] Using first {n_used} of {n_frames} frames to form {n_windows} equal windows.")

    T_all = T_all[:n_used]
    V_all = V_m3_all[:n_used]
    H_all = H_J_all[:n_used]
    m_all = m_kg_all[:n_used]

    # -------- 5) Window-wise stats --------
    Cp_list = []
    alphaP_list = []

    for w in range(n_windows):
        i0 = w * window_size
        i1 = (w + 1) * window_size

        Tw = float(np.mean(T_all[i0:i1]))
        Vw = V_all[i0:i1]
        Hw = H_all[i0:i1]
        mw = float(np.mean(m_all[i0:i1]))

        # Unbiased variance and covariance
        var_Hw = float(np.var(Hw, ddof=1))
        cov_VH = float(np.cov(Vw, Hw, ddof=1)[0, 1])

        Cp_w    = var_Hw / (KB * (Tw ** 2) * mw)             # J/(kg·K)
        alphaP_w = cov_VH / (KB * (Tw ** 2) * float(np.mean(Vw)))  # 1/K

        Cp_list.append(Cp_w)
        alphaP_list.append(alphaP_w)

    Cp_arr    = np.array(Cp_list)
    alphaP_arr = np.array(alphaP_list)
    alphaPL_arr = alphaP_arr / 3.0

    Cp_mean      = float(np.mean(Cp_arr))
    alphaP_mean  = float(np.mean(alphaP_arr))
    alphaPL_mean = float(np.mean(alphaPL_arr))

    return Cp_mean, alphaP_mean, alphaPL_mean


def main():
    """
    Usage:
        python emd_cp_alpha.py <pid> [n_windows]

    Expects: <pid>_npt_properties.dat in cwd.
    Outputs (single-line files with the mean values):
        <pid>_Cp_EMD_result.dat       (Cp_mean, J/(kg·K))
        <pid>_alphaP_EMD_result.dat   (alphaP_mean, 1/K)
        <pid>_alphaPL_EMD_result.dat  (alphaPL_mean, 1/K)
    """
    if len(sys.argv) < 2:
        print("Usage: python emd_cp_alpha.py <pid> [n_windows]")
        sys.exit(1)

    pid = sys.argv[1]
    n_windows = 8

    # Optional integer for n_windows
    for arg in sys.argv[2:]:
        if arg.strip().isdigit():
            n_windows = int(arg)

    infile = os.path.join(os.getcwd(), f"{pid}_npt_properties.dat")
    if not os.path.isfile(infile):
        print(f"Error: input file not found: {infile}")
        sys.exit(1)

    try:
        Cp_mean, alphaP_mean, alphaPL_mean = compute_cp_alpha_emd(
            infile, pid, n_windows=n_windows
        )

        # Print only final mean values
        print(f"PID          : {pid}")
        print(f"n_windows    : {n_windows}")
        print(f"Cp_mean      : {Cp_mean:.8f} J/(kg·K)")
        print(f"alphaP_mean  : {alphaP_mean:.10e} 1/K")
        print(f"alphaPL_mean : {alphaPL_mean:.10e} 1/K")

        # Save only final mean values for pipelines
        np.savetxt(f"{pid}_Cp_EMD_result.dat",      [Cp_mean],      fmt="%.8f")
        np.savetxt(f"{pid}_alphaP_EMD_result.dat",  [alphaP_mean],  fmt="%.10e")
        #np.savetxt(f"{pid}_alphaPL_EMD_result.dat", [alphaPL_mean], fmt="%.10e")

        print("\n[DONE] Saved Cp, alphaP, alphaPL mean values.")

    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
