import os
import sys
import numpy as np

def process_density(pid_file, min_points=8000, fraction=0.25):
    # Read density values, skipping header or comments
    values = []
    with open(pid_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):  # skip empty and comment lines
                continue
            try:
                values.append(float(line))
            except ValueError:
                print(f"Skipping non-numeric line in {pid_file}: {line}")

    if len(values) < min_points:
        print(f"{pid_file} has only {len(values)} points (< {min_points}). Skipping.")
        return None

    # Take the last fraction of values (e.g., 25%)
    n_tail = int(len(values) * fraction)
    values_tail = values[-n_tail:]

    mean_density = np.mean(values_tail)
    return mean_density

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python densCal.py <PID>")
        sys.exit(1)

    pid = sys.argv[1]
    dens_file = f"{pid}_density.dat"

    if not os.path.exists(dens_file):
        print(f"Density file {dens_file} not found.")
        sys.exit(1)

    mean_density = process_density(dens_file, min_points=8000, fraction=0.25)

    if mean_density is not None:
        output_file = f"{pid}_density_result.dat"
        with open(output_file, 'w') as f:
            f.write(f"{mean_density:.6f}\n")
        print(f"Density result saved to {output_file}")
