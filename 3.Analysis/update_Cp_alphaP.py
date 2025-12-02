#!/usr/bin/env python3
# Usage: python update_properties_csv.py PID
#
# Updates:
#   ../../RESULTS/CP_MD.csv       (PID, SMILES, CP)
#   ../../RESULTS/ALPHAP_MD.csv   (PID, SMILES, ALPHAP)   [1/K]
#   ../../RESULTS/ALPHAL_MD.csv   (PID, SMILES, ALPHAL)   [1/K] = ALPHAP/3
#
# Reads, if present:
#   <PID>_Cp_result.dat
#   <PID>_alphaP_EMD_result.dat

import sys, csv, re
from pathlib import Path

# --- Config ---
SMILES_CSV = Path("../../..") / "SMILES.csv"
OUT_DIR    = Path("../../RESULTS")
CP_CSV     = OUT_DIR / "CP_MD.csv"
ALPHAP_CSV = OUT_DIR / "ALPHAP_MD.csv"
ALPHAL_CSV = OUT_DIR / "ALPHAL_MD.csv"

PID_COL, SMILES_COL = "PID", "SMILES"
CP_COL, ALPHAP_COL, ALPHAL_COL = "CP", "ALPHAP", "ALPHAL"

FLOAT_RE = re.compile(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?")

def usage():
    print("Usage: python update_properties_csv.py PID")
    sys.exit(1)

def first_number(text):
    m = FLOAT_RE.search(text)
    return m.group(0) if m else None

def load_smiles(pid: str) -> str:
    if not SMILES_CSV.exists():
        return ""
    with SMILES_CSV.open(newline="", encoding="utf-8") as fh:
        rdr = csv.DictReader(fh)
        if not rdr.fieldnames:
            return ""
        lower = {k.lower(): k for k in rdr.fieldnames}
        pid_k, smi_k = lower.get("pid"), lower.get("smiles")
        if not pid_k or not smi_k:
            return ""
        for row in rdr:
            if (row.get(pid_k) or "").strip() == pid:
                return (row.get(smi_k) or "").strip()
    return ""

def read_scalar(path: Path):
    if not path.exists():
        return None
    try:
        txt = path.read_text(encoding="utf-8", errors="ignore")
    except Exception:
        return None
    s = first_number(txt)
    if s is None:
        return None
    try:
        return float(s)
    except Exception:
        return None

def upsert(csv_path: Path, pid: str, smiles: str, colname: str, value_str: str):
    data = {}
    if csv_path.exists():
        with csv_path.open(newline="", encoding="utf-8") as fh:
            rdr = csv.DictReader(fh)
            if rdr.fieldnames and {PID_COL, SMILES_COL, colname}.issubset(rdr.fieldnames):
                for row in rdr:
                    k = (row.get(PID_COL) or "").strip()
                    if k:
                        data[k] = {
                            PID_COL: k,
                            SMILES_COL: row.get(SMILES_COL, ""),
                            colname: row.get(colname, "")
                        }
    data[pid] = {PID_COL: pid, SMILES_COL: smiles, colname: value_str}
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=[PID_COL, SMILES_COL, colname])
        w.writeheader()
        for v in data.values():
            w.writerow(v)

def main():
    if len(sys.argv) != 2:
        usage()
    pid = sys.argv[1].strip()
    if not pid:
        usage()

    cp_file     = Path(f"{pid}_Cp_EMD_result.dat")
    alphap_file = Path(f"{pid}_alphaP_EMD_result.dat")

    cp_val     = read_scalar(cp_file)
    alphap_val = read_scalar(alphap_file)
    alphal_val = (alphap_val / 3.0) if alphap_val is not None else None

    smiles = load_smiles(pid)

    any_written = False
    if cp_val is not None:
        upsert(CP_CSV, pid, smiles, CP_COL, f"{cp_val:.8f}")
        print(f"[OK] {cp_file.name}     → {CP_CSV.name}:    {CP_COL}    = {cp_val:.8f}")
        any_written = True
    else:
        print(f"[SKIP] {cp_file.name}: not found or unreadable")

    if alphap_val is not None:
        upsert(ALPHAP_CSV, pid, smiles, ALPHAP_COL, f"{alphap_val:.10e}")
        print(f"[OK] {alphap_file.name} → {ALPHAP_CSV.name}: {ALPHAP_COL}= {alphap_val:.10e}")
        any_written = True
    else:
        print(f"[SKIP] {alphap_file.name}: not found or unreadable")

    if alphal_val is not None:
        upsert(ALPHAL_CSV, pid, smiles, ALPHAL_COL, f"{alphal_val:.10e}")
        print(f"[OK] {alphap_file.name} → {ALPHAL_CSV.name}: {ALPHAL_COL}= {alphal_val:.10e} (ALPHAP/3)")
        any_written = True
    else:
        print(f"[SKIP] ALPHAL: needs {alphap_file.name}; no ALPHAP available")

    if not any_written:
        print("[DONE] Nothing updated.")
    else:
        print("[DONE] Updates complete.")

if __name__ == "__main__":
    main()
