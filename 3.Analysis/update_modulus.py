#!/usr/bin/env python3
# update_modulus_csv.py
# Usage: python update_modulus_csv.py PID
#
# Reads PID-prefixed modulus result files in the *current directory*:
#   {PID}_K_result.dat, {PID}_G1_result.dat, {PID}_G2_result.dat,
#   {PID}_G_result.dat, {PID}_nu_result.dat, {PID}_E_result.dat
#
# Writes/updates six separate CSVs in ../../RESULTS/:
#   K_MD.csv, G1_MD.csv, G2_MD.csv, G_MD.csv, NU_MD.csv, E_MD.csv
#
# Each CSV has columns: PID, SMILES, <property>
# SMILES is looked up from ../../../SMILES.csv (expects headers: PID,SMILES)

import sys, csv, re
from pathlib import Path

# --- Config ---
OUT_DIR        = Path("../../RESULTS")
SMILES_CSV     = Path("../../..") / "SMILES.csv"  # run from POLYMER_DATA/MODULUS/<PID>/
PID_COL        = "PID"
SMILES_COL     = "SMILES"

# Input .dat files (PID-prefixed)
PROP_FILES_BASE = {
    "K(GPa)" : "K_result.dat",
    "G1(GPa)": "G1_result.dat",
    "G2(GPa)": "G2_result.dat",
    "G(GPa)" : "G_result.dat",   # combined shear modulus (Shear Modulus = ...)
    "nu"     : "nu_result.dat",
    "E(GPa)" : "E_result.dat",
}

# Output CSV filenames per property
PROP_OUTFILE = {
    "K(GPa)" : "K_MD.csv",
    "G1(GPa)": "G1_MD.csv",
    "G2(GPa)": "G2_MD.csv",
    "G(GPa)" : "G_MD.csv",       # new CSV for combined shear modulus
    "nu"     : "NU_MD.csv",
    "E(GPa)" : "E_MD.csv",
}

FLOAT_RE = re.compile(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?")

def usage():
    print("Usage: python update_modulus_csv.py PID")
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

def read_existing_prop(csv_path: Path, prop_col: str):
    """Read an existing per-property CSV into a dict keyed by PID."""
    rows = {}
    if csv_path.exists():
        with csv_path.open(newline="", encoding="utf-8") as fh:
            rdr = csv.DictReader(fh)
            for row in rdr:
                k = (row.get(PID_COL) or "").strip()
                if not k:
                    continue
                rows[k] = {
                    PID_COL: k,
                    SMILES_COL: (row.get(SMILES_COL) or "").strip(),
                    prop_col: (row.get(prop_col) or "").strip(),
                }
    return rows

def write_prop(csv_path: Path, rows: dict, prop_col: str):
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=[PID_COL, SMILES_COL, prop_col])
        w.writeheader()
        # Write in stable PID order
        for pid_key in sorted(rows.keys()):
            w.writerow(rows[pid_key])

def main():
    if len(sys.argv) != 2:
        usage()
    pid = sys.argv[1].strip()
    if not pid:
        usage()

    # Gather values from PID-prefixed files in current directory
    updates = {}
    for prop_col, base in PROP_FILES_BASE.items():
        fname = f"{pid}_{base}"
        p = Path(fname)
        if not p.exists():
            print(f"[SKIP] {fname}: not found")
            continue
        try:
            txt = p.read_text(encoding="utf-8", errors="ignore")
        except Exception:
            print(f"[SKIP] {fname}: cannot read")
            continue
        val = first_number(txt)
        if val is None:
            print(f"[SKIP] {fname}: no numeric value found")
            continue
        updates[prop_col] = val
        print(f"[OK] {fname}: {prop_col}={val}")

    if not updates:
        print("[INFO] No modulus values found; nothing to write.")
        return

    smiles = load_smiles(pid)

    # Upsert each property into its own CSV
    for prop_col, val in updates.items():
        out_csv = OUT_DIR / PROP_OUTFILE[prop_col]
        rows = read_existing_prop(out_csv, prop_col)

        row = rows.get(pid, {PID_COL: pid, SMILES_COL: "", prop_col: ""})
        row[PID_COL]    = pid
        # only overwrite SMILES if we have a non-empty lookup, else keep existing
        row[SMILES_COL] = smiles or row.get(SMILES_COL, "")
        row[prop_col]   = val
        rows[pid] = row

        write_prop(out_csv, rows, prop_col)
        print(f"[DONE] Updated {out_csv.name} for PID={pid}")

if __name__ == "__main__":
    main()
