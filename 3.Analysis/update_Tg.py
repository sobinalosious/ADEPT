#!/usr/bin/env python3
# Usage: python update_tg_csv.py PID

import sys, csv, re
from pathlib import Path

PID_COL, SMILES_COL, VALUE_COL = "PID", "SMILES", "TG"
SMILES_CSV = Path("../../..") / "SMILES.csv"   # from POLYMER_DATA/GLASS_TRANSITION/<PID>/
OUT_DIR    = Path("../../RESULTS")
OUT_CSV    = OUT_DIR / "TG_MD.csv"
FLOAT_RE   = re.compile(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?")

def usage():
    print("Usage: python update_tg_csv.py PID"); sys.exit(1)

def first_number(text): 
    m = FLOAT_RE.search(text); return m.group(0) if m else None

def load_smiles(pid: str) -> str:
    if not SMILES_CSV.exists(): return ""
    with SMILES_CSV.open(newline="", encoding="utf-8") as fh:
        rdr = csv.DictReader(fh)
        if not rdr.fieldnames: return ""
        lower = {k.lower(): k for k in rdr.fieldnames}
        pid_k, smi_k = lower.get("pid"), lower.get("smiles")
        if not pid_k or not smi_k: return ""
        for row in rdr:
            if (row.get(pid_k) or "").strip() == pid:
                return (row.get(smi_k) or "").strip()
    return ""

def upsert(csv_path: Path, pid: str, smiles: str, value: str):
    data = {}
    if csv_path.exists():
        with csv_path.open(newline="", encoding="utf-8") as fh:
            rdr = csv.DictReader(fh)
            if rdr.fieldnames and {PID_COL, SMILES_COL, VALUE_COL}.issubset(rdr.fieldnames):
                for row in rdr:
                    k = (row.get(PID_COL) or "").strip()
                    if k:
                        data[k] = {PID_COL:k, SMILES_COL:row.get(SMILES_COL,""), VALUE_COL:row.get(VALUE_COL,"")}
    data[pid] = {PID_COL: pid, SMILES_COL: smiles, VALUE_COL: value}
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=[PID_COL, SMILES_COL, VALUE_COL])
        w.writeheader(); [w.writerow(v) for v in data.values()]

def main():
    if len(sys.argv)!=2: usage()
    pid = sys.argv[1].strip()
    if not pid: usage()
    src = Path(f"{pid}_Tg_result.dat")
    if not src.exists(): print(f"[SKIP] {src.name}: not found"); return
    try: val = first_number(src.read_text(encoding="utf-8", errors="ignore"))
    except Exception: print(f"[SKIP] {src.name}: cannot read"); return
    if val is None: print(f"[SKIP] {src.name}: no numeric value found"); return
    smiles = load_smiles(pid)
    upsert(OUT_CSV, pid, smiles, val)
    print(f"[OK] {src.name} → {OUT_CSV.name}: PID={pid} {VALUE_COL}={val}")

if __name__ == "__main__": main()
