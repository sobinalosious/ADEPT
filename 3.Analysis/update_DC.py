#!/usr/bin/env python3
# Usage: python update_props_csv.py PID
# Example: python update_props_csv.py P010003

import sys, csv, re
from pathlib import Path

# --- Column names (fixed) ---
PID_COL, SMILES_COL = "PID", "SMILES"

# --- Locations (relative to POLYMER_DATA/DIELECTRIC_CONSTANT/<PID>/) ---
SMILES_CSV = Path("../../..") / "SMILES.csv"   # up to the project root
OUT_DIR    = Path("../../RESULTS")             # results directory

# --- Regex to grab the first float in a text file ---
FLOAT_RE = re.compile(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?")

# --- Properties to process (filename tail, output CSV, value column) ---
PROPS = [
    ("refractive_index_polymer_result.dat",               "RI_MD.csv",  "RI"),
    ("dielectric_constant_dipole_result.dat",             "DCD_MD.csv", "DCD"),
    ("dielectric_constant_electronic_polymer_result.dat", "DCE_MD.csv", "DCE"),
    ("dielectric_constant_total_result.dat",              "DC_MD.csv",  "DC"),
    ("permittivity_total_result.dat",                     "PE_MD.csv",  "PE"),
]

def usage():
    print("Usage: python update_props_csv.py PID")
    sys.exit(1)

def first_number(text: str):
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

def upsert(csv_path: Path, pid: str, smiles: str, value_col: str, value: str):
    data = {}
    # read existing rows if the header matches
    if csv_path.exists():
        with csv_path.open(newline="", encoding="utf-8") as fh:
            rdr = csv.DictReader(fh)
            needed = {PID_COL, SMILES_COL, value_col}
            if rdr.fieldnames and needed.issubset(rdr.fieldnames):
                for row in rdr:
                    k = (row.get(PID_COL) or "").strip()
                    if k:
                        data[k] = {
                            PID_COL:    k,
                            SMILES_COL: row.get(SMILES_COL, ""),
                            value_col:  row.get(value_col, "")
                        }
    # upsert this PID
    data[pid] = {PID_COL: pid, SMILES_COL: smiles, value_col: value}

    # (re)write
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=[PID_COL, SMILES_COL, value_col])
        w.writeheader()
        for v in data.values():
            w.writerow(v)

def process_prop(pid: str, smiles: str, src_tail: str, out_csv_name: str, value_col: str):
    src = Path(f"{pid}_{src_tail}")
    if not src.exists():
        print(f"[SKIP] {src.name}: not found"); return
    try:
        val = first_number(src.read_text(encoding="utf-8", errors="ignore"))
    except Exception:
        print(f"[SKIP] {src.name}: cannot read"); return
    if val is None:
        print(f"[SKIP] {src.name}: no numeric value found"); return
    out_csv = OUT_DIR / out_csv_name
    upsert(out_csv, pid, smiles, value_col, val)
    print(f"[OK] {src.name} → {out_csv.name}: PID={pid} {value_col}={val}")

def main():
    if len(sys.argv) != 2:
        usage()
    pid = sys.argv[1].strip()
    if not pid:
        usage()

    smiles = load_smiles(pid)
    for src_tail, out_csv_name, value_col in PROPS:
        process_prop(pid, smiles, src_tail, out_csv_name, value_col)

if __name__ == "__main__":
    main()
