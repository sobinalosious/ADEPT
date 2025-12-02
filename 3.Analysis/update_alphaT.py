#!/usr/bin/env python3
# Usage: python update_alpha_diffusivity.py PID
#
# Computes thermal diffusivity alpha = k / (rho * Cp)
# and upserts into: ../../RESULTS/ALPHAT_MD.csv
#
# Inputs (../../RESULTS/):
#   - CP_MD.csv     (PID, SMILES, CP)       [J/(kg·K)]
#   - TC_MD.csv     (PID, SMILES, TC)       [W/(m·K)]
#   - RHO_MD.csv    (PID, SMILES, RHO)      [kg/m^3]
#
# Output:
#   - ALPHAT_MD.csv  (PID, SMILES, ALPHAT)  [m^2/s]

import sys, csv
from pathlib import Path
from typing import Set, Optional

RESULTS_DIR = Path("../../RESULTS")
SMILES_CSV  = Path("../../..") / "SMILES.csv"

CP_CSV    = RESULTS_DIR / "CP_MD.csv"
TC_CSV    = RESULTS_DIR / "TC_MD.csv"
RHO_CSV   = RESULTS_DIR / "RHO_MD.csv"
OUT_CSV   = RESULTS_DIR / "ALPHAT_MD.csv"

PID_COL, SMILES_COL = "PID", "SMILES"
CP_COL, TC_COL, RHO_COL, ALPHAT_COL = "CP", "TC", "RHO", "ALPHAT"

# Aliases (lowercase) to be tolerant of varied headers
CP_ALIASES   = {"cp", "cp_md", "cp_value", "cp_j_per_kgk"}
TC_ALIASES   = {"tc", "tc_md", "k", "thermal_conductivity"}
RHO_ALIASES  = {"rho", "rho_md", "density", "rho_kg_m3"}

def usage():
    print("Usage: python update_alpha_diffusivity.py PID")
    sys.exit(1)

def _open_csv(path: Path):
    try:
        return path.open(newline="", encoding="utf-8")
    except UnicodeError:
        return path.open(newline="", encoding="utf-8-sig")

def _lower_map(fields):
    return {f.lower(): f for f in (fields or [])}

def _read_value(csv_path: Path, pid: str, aliases: Set[str]) -> Optional[float]:
    """Find numeric value for PID in CSV with possible alias columns."""
    if not csv_path.exists():
        return None
    with _open_csv(csv_path) as fh:
        rdr = csv.DictReader(fh)
        if not rdr.fieldnames:
            return None
        low = _lower_map(rdr.fieldnames)
        pid_key = low.get(PID_COL.lower())
        if not pid_key:
            return None
        val_key = None
        for a in aliases:
            if a in low:
                val_key = low[a]
                break
        if not val_key:
            return None
        for row in rdr:
            if (row.get(pid_key) or "").strip() == pid:
                raw = (row.get(val_key) or "").strip()
                try:
                    return float(raw)
                except Exception:
                    return None
    return None

def _load_smiles(pid: str) -> str:
    if not SMILES_CSV.exists():
        return ""
    with _open_csv(SMILES_CSV) as fh:
        rdr = csv.DictReader(fh)
        if not rdr.fieldnames:
            return ""
        low = _lower_map(rdr.fieldnames)
        pid_k, smi_k = low.get("pid"), low.get("smiles")
        if not pid_k or not smi_k:
            return ""
        for row in rdr:
            if (row.get(pid_k) or "").strip() == pid:
                return (row.get(smi_k) or "").strip()
    return ""

def _upsert(csv_path: Path, pid: str, smiles: str, value: float):
    data = {}
    if csv_path.exists():
        with _open_csv(csv_path) as fh:
            rdr = csv.DictReader(fh)
            if rdr.fieldnames:
                for row in rdr:
                    k = (row.get(PID_COL) or "").strip()
                    if k:
                        data[k] = {
                            PID_COL: k,
                            SMILES_COL: row.get(SMILES_COL, ""),
                            ALPHAT_COL: row.get(ALPHAT_COL, ""),
                        }
    data[pid] = {
        PID_COL: pid,
        SMILES_COL: smiles,
        ALPHAT_COL: f"{value:.6e}",
    }

    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=[PID_COL, SMILES_COL, ALPHAT_COL])
        w.writeheader()
        for k in sorted(data.keys()):
            w.writerow(data[k])

def main():
    if len(sys.argv) != 2:
        usage()
    pid = sys.argv[1].strip()
    if not pid:
        usage()

    cp  = _read_value(CP_CSV,  pid, CP_ALIASES | {CP_COL.lower()})
    tc  = _read_value(TC_CSV,  pid, TC_ALIASES | {TC_COL.lower()})
    rho = _read_value(RHO_CSV, pid, RHO_ALIASES | {RHO_COL.lower()})

    if cp is None:
        print(f"[SKIP] No Cp for {pid}")
        return
    if tc is None:
        print(f"[SKIP] No TC for {pid}")
        return
    if rho is None:
        print(f"[SKIP] No RHO for {pid}")
        return

    try:
        alphat = tc / (rho * cp)   # Cp is strictly J/(kg·K)
    except Exception as e:
        print(f"[SKIP] Failed calc for {pid}: {e}")
        return

    smiles = _load_smiles(pid)
    _upsert(OUT_CSV, pid, smiles, alphat)
    print(
        f"[OK] PID={pid}: Cp={cp:.3f} J/(kg·K), TC={tc:.3f} W/(m·K), "
        f"RHO={rho:.3f} kg/m^3 → ALPHAT={alphat:.6e} m^2/s"
    )

if __name__ == "__main__":
    main()
