import os
import sys

def extract_last_rg_values(pid_file, start_line=5598, end_line=5603):
    with open(pid_file, 'r') as file:
        lines = file.readlines()

    # Extract the specified lines
    data_lines = lines[start_line - 1:end_line]

    if not data_lines or len(data_lines) != (end_line - start_line + 1):
        print("Error: Specified lines are not valid or incomplete.")
        return None

    try:
        rg_values = [float(line.strip().split()[1]) for line in data_lines]
    except (ValueError, IndexError) as e:
        print(f"Error parsing Rg values: {e}")
        return None

    rg_mean = sum(rg_values) / len(rg_values)
    return rg_mean


def process_rg_results(pid):
    pid_file = f"{pid}_Rg.data"
    if not os.path.exists(pid_file):
        print(f"Missing: {pid_file}")
        return

    rg_mean = extract_last_rg_values(pid_file, 5598, 5603)
    if rg_mean is None:
        return

    # Save to new file
    output_file = f"{pid}_Rg_result.dat"
    with open(output_file, 'w') as f:
        f.write(f"{rg_mean:.6f}\n")

    print(f"Rg mean value for {pid} saved to {output_file}")


# Entry point
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python update_rg.py <PID>")
    else:
        process_rg_results(sys.argv[1])
