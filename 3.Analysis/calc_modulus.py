#!/usr/bin/env python3
import re
import sys
import os

# -------- regexes (robust to spaces & scientific notation) --------
FLOAT = r'([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)'

# Start-of-block: Bulk Modulus line
RE_K      = re.compile(rf'^\s*Bulk\s*Modulus\s*\d*\s*=\s*{FLOAT}', re.I)

# Individual shear moduli
RE_G1     = re.compile(rf'^\s*Shear\s*Modulus\s*1\s*=\s*{FLOAT}', re.I)
RE_G2     = re.compile(rf'^\s*Shear\s*Modulus\s*2\s*=\s*{FLOAT}', re.I)

# Combined shear modulus "Shear Modulus = ..."
RE_G      = re.compile(rf'^\s*Shear\s*Modulus\s*=\s*{FLOAT}', re.I)

RE_NU     = re.compile(rf'^\s*Poisson\s*Ratio\s*=\s*{FLOAT}', re.I)
# Allow both "Youngs" and "Young's"
RE_E      = re.compile(rf"^\s*Young'?s\s*Modulus\s*=\s*{FLOAT}", re.I)

def main():
    if len(sys.argv) != 2:
        print(f"Usage: {os.path.basename(sys.argv[0])} PID", file=sys.stderr)
        sys.exit(1)
    pid = sys.argv[1]

    logfile = "log.elastic_T"
    if not os.path.exists(logfile):
        print(f"Error: {logfile} not found in current directory.", file=sys.stderr)
        sys.exit(1)

    # Store per-block values
    K_vals, G1_vals, G2_vals, G_vals, NU_vals, E_vals = [], [], [], [], [], []

    current = None  # holds the current block's values

    with open(logfile, "r") as f:
        for line in f:
            # --- Start of a new block: Bulk Modulus line ---
            m = RE_K.search(line)
            if m:
                # If previous block is complete, save it
                if current and all(v is not None for v in current.values()):
                    K_vals.append(current["K"])
                    G1_vals.append(current["G1"])
                    G2_vals.append(current["G2"])
                    G_vals.append(current["G"])
                    NU_vals.append(current["NU"])
                    E_vals.append(current["E"])

                # Start a new block
                current = {"K": None, "G1": None, "G2": None, "G": None, "NU": None, "E": None}
                current["K"] = float(m.group(1))
                continue

            # If we haven't started a block yet, skip
            if current is None:
                continue

            # --- Within a block: try to match each property line ---

            # Shear Modulus 1
            m = RE_G1.search(line)
            if m:
                current["G1"] = float(m.group(1))
                continue

            # Shear Modulus 2
            m = RE_G2.search(line)
            if m:
                current["G2"] = float(m.group(1))
                continue

            # Combined Shear Modulus
            m = RE_G.search(line)
            if m:
                # This line is "Shear Modulus = ..."
                current["G"] = float(m.group(1))
                continue

            # Poisson Ratio
            m = RE_NU.search(line)
            if m:
                current["NU"] = float(m.group(1))
                continue

            # Youngs Modulus
            m = RE_E.search(line)
            if m:
                current["E"] = float(m.group(1))
                continue

    # Flush the last block if it is complete
    if current and all(v is not None for v in current.values()):
        K_vals.append(current["K"])
        G1_vals.append(current["G1"])
        G2_vals.append(current["G2"])
        G_vals.append(current["G"])
        NU_vals.append(current["NU"])
        E_vals.append(current["E"])

    if not K_vals:
        print("Error: no complete Bulk/Shear/Poisson/Youngs blocks found.", file=sys.stderr)
        sys.exit(2)

    n_blocks = len(K_vals)

    # Compute averages
    K_avg  = sum(K_vals)  / n_blocks
    G1_avg = sum(G1_vals) / n_blocks
    G2_avg = sum(G2_vals) / n_blocks
    G_avg  = sum(G_vals)  / n_blocks
    NU_avg = sum(NU_vals) / n_blocks
    E_avg  = sum(E_vals)  / n_blocks

    def write_scalar(filename, value):
        with open(filename, "w") as fh:
            fh.write(f"{value:.10g}\n")  # single average value

    # Write one averaged value to each *_result.dat file
    write_scalar(f"{pid}_K_result.dat",   K_avg)
    write_scalar(f"{pid}_G1_result.dat",  G1_avg)
    write_scalar(f"{pid}_G2_result.dat",  G2_avg)
    write_scalar(f"{pid}_G_result.dat",   G_avg)   # combined shear modulus
    write_scalar(f"{pid}_nu_result.dat",  NU_avg)
    write_scalar(f"{pid}_E_result.dat",   E_avg)

    # Optional: quick confirmation
    print(
        f"Parsed {n_blocks} complete blocks. "
        f"Wrote averaged values to: "
        f"{pid}_K_result.dat, {pid}_G1_result.dat, {pid}_G2_result.dat, "
        f"{pid}_G_result.dat, {pid}_nu_result.dat, {pid}_E_result.dat"
    )

if __name__ == "__main__":
    main()
