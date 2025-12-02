#!/usr/bin/env python3
"""
Predict Tg using a pre-trained MLP + MACCS model for a given PID.

Requirements:
    - Trained model : model_mlp_tg.pkl
    - SMILES.csv file with columns: PID, SMILES

Usage:
    python predict_tg.py P000001
"""

import sys
import os
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import MACCSkeys, DataStructs
import joblib

# ----------------- Config -----------------
# Path to the MLP model you trained in train_tg_mlp_maccs.py
MODEL_PATH = "Files/model_mlp_tg.pkl"
SMILES_CSV_PATH = "SMILES.csv"
# ------------------------------------------


def load_model(model_path: str):
    """Load the pre-trained MLP model from disk."""
    if not os.path.exists(model_path):
        raise FileNotFoundError(f"Model file not found: {model_path}")
    return joblib.load(model_path)


def load_smiles_for_pid(pid: str, smiles_csv_path: str) -> str:
    """Return the SMILES string corresponding to the given PID."""
    if not os.path.exists(smiles_csv_path):
        raise FileNotFoundError(f"SMILES CSV not found: {smiles_csv_path}")

    df = pd.read_csv(smiles_csv_path)

    if "PID" not in df.columns or "SMILES" not in df.columns:
        raise ValueError("SMILES.csv must contain 'PID' and 'SMILES' columns.")

    # Compare as strings to be robust
    df["PID"] = df["PID"].astype(str)
    pid_str = str(pid)

    row = df.loc[df["PID"] == pid_str]
    if row.empty:
        raise ValueError(f"PID '{pid_str}' not found in {smiles_csv_path}.")

    smiles = row["SMILES"].iloc[0]
    if not isinstance(smiles, str) or smiles.strip() == "":
        raise ValueError(f"Empty or invalid SMILES for PID '{pid_str}'.")

    return smiles


def smiles_to_maccs(smiles: str) -> np.ndarray:
    """Convert a single SMILES string to a MACCS fingerprint row (shape: 1 x n_bits)."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string: {smiles}")

    fp = MACCSkeys.GenMACCSKeys(mol)  # 167 bits in RDKit
    arr = np.zeros((fp.GetNumBits(),), dtype=np.int8)
    DataStructs.ConvertToNumpyArray(fp, arr)

    # Shape (1, n_features) for sklearn
    return np.array([arr], dtype=np.float32)


def main():
    if len(sys.argv) != 2:
        print("Usage: python predict_tg_mlp_maccs.py <PID>")
        sys.exit(1)

    pid = sys.argv[1]

    # Load pre-trained MLP model
    model = load_model(MODEL_PATH)

    # Get SMILES for this PID
    smiles = load_smiles_for_pid(pid, SMILES_CSV_PATH)

    # Convert SMILES to MACCS fingerprint
    X = smiles_to_maccs(smiles)

    # Predict Tg
    tg_pred = model.predict(X)[0]

    # Print only the numeric value (easy to capture in bash)
    print(f"{tg_pred:.6f}")


if __name__ == "__main__":
    main()
