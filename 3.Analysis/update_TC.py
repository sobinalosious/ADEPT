#!/usr/bin/env python3
# update_tc_csv.py
# Usage: python update_tc_csv.py PID
#
# Reads thermal conductivity result in the *current directory*:
#   {PID}_TC_result.dat  -> ../../RESULTS/TC_MD.csv  (columns: PID,SMILES,TC)
#
# Looks up SMILES for PID from ../../../SMILES.csv (expects headers: PID,SMILES)

import sys, csv, re
from pathlib import Path

PID_COL, SMILES_COL, VALUE_COL = "PID", "SMILES", "TC"
SMILES_CSV = Path("../../..") / "SMILES.csv"  # from POLYMER_DATA/THERMAL_CONDUCTIVITY/<PID>/
OUT_DIR    = Path("../../RESULTS")
OUT_CSV    = OUT_DIR / "TC_MD.csv"
FLOAT_RE   = re.compile(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?")

def usage():
    print("Usage: python update_tc_csv.py PID")
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
        pid_key = lower.get("pid"); smi_key = lower.get("smiles")
        if not pid_key or not smi_key:
            return ""
        for row in rdr:
            if (row.get(pid_key) or "").strip() == pid:
                return (row.get(smi_key) or "").strip()
    return ""

def upsert(csv_path: Path, pid: str, smiles: str, value: str):
    # Simple upsert keyed by PID; expects columns PID,SMILES,TC
    data = {}
    if csv_path.exists():
        with csv_path.open(newline="", encoding="utf-8") as fh:
            rdr = csv.DictReader(fh)
            if rdr.fieldnames and {PID_COL, SMILES_COL, VALUE_COL}.issubset(rdr.fieldnames):
                for row in rdr:
                    k = (row.get(PID_COL) or "").strip()
                    if k:
                        data[k] = {
                            PID_COL: k,
                            SMILES_COL: row.get(SMILES_COL, ""),
                            VALUE_COL: row.get(VALUE_COL, ""),
                        }

    data[pid] = {PID_COL: pid, SMILES_COL: smiles, VALUE_COL: value}

    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=[PID_COL, SMILES_COL, VALUE_COL])
        w.writeheader()
        for row in data.values():
            w.writerow(row)

def main():
    if len(sys.argv) != 2:
        usage()
    pid = sys.argv[1].strip()
    if not pid:
        usage()

    src = Path(f"{pid}_TC_result.dat")
    if not src.exists():
        print(f"[SKIP] {src.name}: not found")
        sys.exit(0)

    try:
        text = src.read_text(encoding="utf-8", errors="ignore")
    except Exception:
        print(f"[SKIP] {src.name}: cannot read")
        sys.exit(0)

    val = first_number(text)
    if val is None:
        print(f"[SKIP] {src.name}: no numeric value found")
        sys.exit(0)

    smiles = load_smiles(pid)
    upsert(OUT_CSV, pid, smiles, val)
    print(f"[OK] {src.name} → {OUT_CSV.name}: PID={pid} {VALUE_COL}={val}")

if __name__ == "__main__":
    main()
