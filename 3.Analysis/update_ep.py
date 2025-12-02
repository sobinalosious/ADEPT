#!/usr/bin/env python3
# update_electronic_csvs.py
# Usage: python update_electronic_csvs.py PID
#
# Reads electronic property files in the *current directory* named like:
#   {PID}_{key}.dat  (or fallback: {key}.dat)
# Looks up SMILES for PID from ../../../SMILES.csv (expects headers: PID,SMILES)
# Appends/updates one CSV per property under ../../RESULTS/
# Schema for each CSV: PID,SMILES,<METRIC>  where <METRIC> = CSV name without _DFT
#
# Short CSV names (exactly one "_" before DFT):
#   avg_polarizability                     -> ALPHA_DFT.csv   (column: ALPHA)
#   refractive_index_monomer               -> ERIMONO_DFT.csv  (column: RIMONO)
#   electronic_dielectric_constant_monomer -> EDCMONO_DFT.csv (column: EDCMONO)
#   total_energy                           -> ETOTAL_DFT.csv  (column: ETOTAL)
#   homo                                   -> HOMO_DFT.csv    (column: HOMO)
#   lumo                                   -> LUMO_DFT.csv    (column: LUMO)
#   dipole_moment                          -> MU_DFT.csv      (column: MU)
# Derived:
#   permittivity (ε = ε0 * dielectric constant) -> EPERMONO_DFT.csv (column: EPERMONO)
#   bandgap (Eg = LUMO - HOMO)                   -> BANDGAP_DFT.csv  (column: BANDGAP)

import sys, csv, re
from pathlib import Path

# ---------- config ----------
PROPS = [
    ("avg_polarizability", "ALPHA_DFT.csv"),
    ("refractive_index_monomer", "ERIMONO_DFT.csv"),
    ("electronic_dielectric_constant_monomer", "EDCMONO_DFT.csv"),
    ("total_energy", "ETOTAL_DFT.csv"),
    ("homo", "HOMO_DFT.csv"),
    ("lumo", "LUMO_DFT.csv"),
    ("dipole_moment", "MU_DFT.csv"),
]

# Derived CSV names
CSV_PERMITTIVITY = "EPERMONO_DFT.csv"   # ε = ε0 * ε_r
CSV_BANDGAP      = "BANDGAP_DFT.csv"    # Eg = LUMO - HOMO

EPS0 = 8.854187817e-12  # F/m

SMILES_CSV = Path("../../..") / "SMILES.csv"   # from POLYMER_DATA/MONOMER_ELECTRONIC/<PID>/
OUT_DIR    = Path("../../RESULTS")             # write aggregate CSVs here

PID_COL, SMILES_COL = "PID", "SMILES"
FLOAT_RE = re.compile(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?")

def usage():
    print("Usage: python update_electronic_csvs.py PID")
    sys.exit(1)

def first_number(text: str):
    m = FLOAT_RE.search(text)
    return m.group(0) if m else None

def value_col_from_csv_name(csv_name: str) -> str:
    """
    Turn 'ALPHA_DFT.csv' -> 'ALPHA', 'ETOTAL_DFT.csv' -> 'ETOTAL', 'EDCMONO_DFT.csv' -> 'EDCMONO'
    """
    stem = csv_name.rsplit(".", 1)[0]
    if stem.upper().endswith("_DFT"):
        return stem[:-4]  # drop '_DFT'
    return stem

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
    """
    Upsert PID row with schema [PID, SMILES, <value_col>].
    If file exists with a different value column, normalize to value_col.
    """
    rows = []
    found = False
    existing_val_col = value_col
    if csv_path.exists():
        with csv_path.open(newline="", encoding="utf-8") as fh:
            rdr = csv.DictReader(fh)
            if rdr.fieldnames:
                for c in rdr.fieldnames:
                    if c not in (PID_COL, SMILES_COL):
                        existing_val_col = c
                        break
            for row in rdr:
                val_here = row.get(existing_val_col, "")
                norm_row = {
                    PID_COL: (row.get(PID_COL) or "").strip(),
                    SMILES_COL: row.get(SMILES_COL) or "",
                    value_col: val_here,
                }
                if norm_row[PID_COL] == pid:
                    norm_row[SMILES_COL] = smiles
                    norm_row[value_col] = value
                    found = True
                rows.append(norm_row)
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

    cwd = Path.cwd()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    smiles = load_smiles(pid)

    # Read base properties and cache numeric values
    cached = {}  # key -> float
    for key, csv_name in PROPS:
        f1 = cwd / f"{pid}_{key}.dat"
        f2 = cwd / f"{key}.dat"
        src = f1 if f1.exists() else f2 if f2.exists() else None
        if src is None:
            print(f"[SKIP] {key}: file not found ({f1.name} / {f2.name})")
            continue

        try:
            text = src.read_text(encoding="utf-8", errors="ignore")
        except Exception:
            print(f"[SKIP] {key}: could not read {src.name}")
            continue

        val_str = first_number(text)
        if val_str is None:
            print(f"[SKIP] {key}: no numeric value in {src.name}")
            continue

        try:
            cached[key] = float(val_str)
        except ValueError:
            pass

        out_csv = OUT_DIR / csv_name
        value_col = value_col_from_csv_name(csv_name)
        upsert(out_csv, pid, smiles, val_str, value_col)
        print(f"[OK] {key} → {csv_name}: PID={pid} {value_col}={val_str}")

    # Derived: PERMITTIVITY from electronic dielectric constant (ε = ε0 * ε_r)
    eps_r = cached.get("electronic_dielectric_constant_monomer")
    if eps_r is not None:
        eps_abs = EPS0 * eps_r
        csv_name = CSV_PERMITTIVITY                      # EPERMONO_DFT.csv
        out_csv = OUT_DIR / csv_name
        value_col = value_col_from_csv_name(csv_name)    # "EPERMONO"
        upsert(out_csv, pid, smiles, f"{eps_abs}", value_col)
        print(f"[OK] permittivity → {csv_name}: PID={pid} {value_col}={eps_abs}")
    else:
        print(f"[SKIP] permittivity: missing electronic_dielectric_constant_monomer")

    # Derived: BANDGAP from HOMO/LUMO (Eg = LUMO - HOMO)
    homo = cached.get("homo")
    lumo = cached.get("lumo")
    if (homo is not None) and (lumo is not None):
        bandgap = lumo - homo
        csv_name = CSV_BANDGAP
        out_csv = OUT_DIR / csv_name
        value_col = value_col_from_csv_name(csv_name)    # "BANDGAP"
        upsert(out_csv, pid, smiles, f"{bandgap}", value_col)
        print(f"[OK] bandgap → {csv_name}: PID={pid} {value_col}={bandgap}")
    else:
        print(f"[SKIP] bandgap: missing homo and/or lumo")

if __name__ == "__main__":
    main()
