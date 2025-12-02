#!/usr/bin/env python3
# update_visc_diff_csv.py
# Usage: python update_visc_diff_csv.py PID
#
# Reads result files in the *current directory*:
#   {PID}_viscosity_result.dat   -> ../../RESULTS/VISC_MD.csv  (columns: TEMP, VISC)
#   {PID}_diffusion_result.dat   -> ../../RESULTS/D_MD.csv     (columns: TEMP, D)
#
# Looks up SMILES for PID from ../../../SMILES.csv (expects headers: PID,SMILES)

import sys, csv, re
from pathlib import Path
from typing import Optional, Tuple, List

PID_COL, SMILES_COL, TEMP_COL = "PID", "SMILES", "TEMP"
FLOAT_RE = re.compile(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?")

# Locations (relative to POLYMER_DATA/VISCOSITY/<PID>/)
SMILES_CSV = Path("../../..") / "SMILES.csv"
OUT_DIR    = Path("../../RESULTS")

# Source → target mapping
JOBS = [
    # (source_file_pattern, out_csv_name, value_col_name)
    ("{pid}_viscosity_result.dat", "VISC_MD.csv", "VISC"),
    ("{pid}_diffusion_result.dat", "D_MD.csv",    "D"),
]

def usage():
    print("Usage: python update_visc_diff_csv.py PID")
    sys.exit(1)

def first_two_numbers(text: str) -> Optional[Tuple[str, str]]:
    nums: List[str] = FLOAT_RE.findall(text)
    if len(nums) < 2:
        return None
    return nums[0], nums[1]  # (TEMP, VALUE)

def load_smiles(pid: str) -> str:
    if not SMILES_CSV.exists():
        return ""
    with SMILES_CSV.open(newline="", encoding="utf-8") as fh:
        rdr = csv.DictReader(fh)
        if not rdr.fieldnames:
            return ""
        lower = {k.lower(): k for k in rdr.fieldnames}
        pid_key = lower.get("pid"); smi_key = lower.get("smiles")
        if not pid_key or not smi_key:
            return ""
        for row in rdr:
            if (row.get(pid_key) or "").strip() == pid:
                return (row.get(smi_key) or "").strip()
    return ""

def upsert(csv_path: Path, pid: str, smiles: str, temp: str, value: str, value_col: str):
    """
    Upsert by PID into csv_path with standardized columns:
      [PID, SMILES, TEMP, <value_col>]
    If an existing file has a different value column name, normalize it to `value_col`.
    If an existing file lacks TEMP, create/normalize it.
    """
    rows = []
    found = False
    existing_val_col = value_col
    existing_has_temp = False

    if csv_path.exists():
        with csv_path.open(newline="", encoding="utf-8") as fh:
            rdr = csv.DictReader(fh)
            if rdr.fieldnames:
                # Determine existing value column (first non PID/SMILES/TEMP)
                for c in rdr.fieldnames:
                    if c == TEMP_COL:
                        existing_has_temp = True
                for c in rdr.fieldnames:
                    if c not in (PID_COL, SMILES_COL, TEMP_COL):
                        existing_val_col = c
                        break

            for row in rdr:
                norm = {
                    PID_COL: (row.get(PID_COL) or "").strip(),
                    SMILES_COL: row.get(SMILES_COL) or "",
                    TEMP_COL: (row.get(TEMP_COL) or "") if existing_has_temp else "",
                    value_col: row.get(existing_val_col, ""),
                }
                if norm[PID_COL] == pid:
                    # Update this PID
                    norm[SMILES_COL] = smiles
                    norm[TEMP_COL] = temp
                    norm[value_col] = value
                    found = True
                rows.append(norm)

    if not found:
        rows.append({PID_COL: pid, SMILES_COL: smiles, TEMP_COL: temp, value_col: value})

    # Write back with standardized header
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = csv_path.with_suffix(csv_path.suffix + ".tmp")
    with tmp.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=[PID_COL, SMILES_COL, TEMP_COL, value_col])
        w.writeheader()
        for row in rows:
            w.writerow({
                PID_COL: row.get(PID_COL, ""),
                SMILES_COL: row.get(SMILES_COL, ""),
                TEMP_COL: row.get(TEMP_COL, ""),
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

        pair = first_two_numbers(text)
        if pair is None:
            print(f"[SKIP] {src.name}: expected two numbers (TEMP and {value_col})")
            continue
        temp, val = pair

        out_csv = OUT_DIR / out_csv_name
        upsert(out_csv, pid, smiles, temp, val, value_col)
        print(f"[OK] {src.name} → {out_csv_name}: PID={pid} TEMP={temp} {value_col}={val}")

if __name__ == "__main__":
    main()
