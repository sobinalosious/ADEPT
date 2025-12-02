#!/usr/bin/env python3
# update_md_csvs.py
# Usage: python update_md_csvs.py PID
#
# Reads MD result files in the *current directory*:
#   {PID}_density_result.dat  -> RHO_MD.csv  (column: RHO)
#   {PID}_Rg_result.dat       -> RG_MD.csv   (column: RG)
#
# Looks up SMILES for PID from ../../../SMILES.csv (expects headers: PID,SMILES)
# Appends/updates CSVs under ../../RESULTS/
# Schema: PID,SMILES,<VALUE_COL>  where VALUE_COL is RHO or RG

import sys, csv, re
from pathlib import Path

PID_COL, SMILES_COL = "PID", "SMILES"
FLOAT_RE = re.compile(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?")

# Fixed locations (relative to POLYMER_DATA/OPTIMIZATION/<PID>/)
SMILES_CSV = Path("../../..") / "SMILES.csv"
OUT_DIR    = Path("../../RESULTS")

# Source → target mapping
JOBS = [
    # (source_file_pattern, out_csv_name, value_col_name)
    ("{pid}_density_result.dat", "RHO_MD.csv", "RHO"),
    ("{pid}_Rg_result.dat",      "RG_MD.csv",  "RG"),
]

def usage():
    print("Usage: python update_md_csvs.py PID")
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
        lower_map = {k.lower(): k for k in rdr.fieldnames}
        pid_key = lower_map.get("pid"); smi_key = lower_map.get("smiles")
        if not pid_key or not smi_key:
            return ""
        for row in rdr:
            if (row.get(pid_key) or "").strip() == pid:
                return (row.get(smi_key) or "").strip()
    return ""

def upsert(csv_path: Path, pid: str, smiles: str, value: str, value_col: str):
    rows = []
    found = False
    # Read existing (if any) and normalize to desired value_col
    if csv_path.exists():
        with csv_path.open(newline="", encoding="utf-8") as fh:
            rdr = csv.DictReader(fh)
            existing_val_col = value_col
            if rdr.fieldnames:
                for c in rdr.fieldnames:
                    if c not in (PID_COL, SMILES_COL):
                        existing_val_col = c
                        break
            for row in rdr:
                v = row.get(existing_val_col, "")
                norm = {
                    PID_COL: (row.get(PID_COL) or "").strip(),
                    SMILES_COL: row.get(SMILES_COL) or "",
                    value_col: v,
                }
                if norm[PID_COL] == pid:
                    norm[SMILES_COL] = smiles
                    norm[value_col] = value
                    found = True
                rows.append(norm)
    if not found:
        rows.append({PID_COL: pid, SMILES_COL: smiles, value_col: value})

    tmp = csv_path.with_suffix(csv_path.suffix + ".tmp")
    with tmp.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=[PID_COL, SMILES_COL, value_col])
        w.writeheader()
        for row in rows:
            w.writerow({
                PID_COL: row.get(PID_COL, ""),
                SMILES_COL: row.get(SMILES_COL, ""),
                value_col: row.get(value_col, ""),
            })
    tmp.replace(csv_path)

def main():
    if len(sys.argv) != 2:
        usage()
    pid = sys.argv[1].strip()
    if not pid:
        usage()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    smiles = load_smiles(pid)

    for pattern, out_csv_name, value_col in JOBS:
        src = Path(pattern.format(pid=pid))
        if not src.exists():
            print(f"[SKIP] {src.name}: not found")
            continue

        try:
            text = src.read_text(encoding="utf-8", errors="ignore")
        except Exception:
            print(f"[SKIP] {src.name}: cannot read")
            continue

        val = first_number(text)
        if val is None:
            print(f"[SKIP] {src.name}: no numeric value found")
            continue

        out_csv = OUT_DIR / out_csv_name
        upsert(out_csv, pid, smiles, val, value_col)
        print(f"[OK] {src.name} → {out_csv_name}: PID={pid} {value_col}={val}")

if __name__ == "__main__":
    main()
