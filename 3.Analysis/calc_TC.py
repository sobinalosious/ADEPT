import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy import stats

# Step 1: Split thermo.dat into hf*.txt
def split_hf_file(input_path, output_path, pid, num_points=2000):
    thermo_file = os.path.join(input_path, f"{pid}_thermo.dat")
    pid_output_path = os.path.join(output_path, pid + '_HT')

    if not os.path.exists(thermo_file):
        print(f"File {thermo_file} not found.")
        return False

    if not os.path.exists(pid_output_path):
        os.makedirs(pid_output_path)

    # Read thermo file directly
    data = pd.read_csv(thermo_file, sep=r"\s+", header=None,
                       names=['Step','Temp','E_pair','E_mol','TotEng','Press',
                              'Density','Volume','Lx','Ly','Lz','f_4','f_5'])

    # Need at least 2000 data points
    if len(data) < num_points:
        print(f"Not enough data points in {thermo_file}. Expected {num_points} data points.")
        return False

    # Select last `num_points` rows (instead of searching headers)
    data = data.tail(num_points)

    # Split into 10 equal parts
    batch = int(len(data) / 10)
    for i in range(10):
        file_name = os.path.join(pid_output_path, f"hf{i+1}.txt")
        data.iloc[i*batch:(i+1)*batch].to_csv(file_name, sep=" ", index=False, header=False)

    return True

# Step 2: Split temperature profile
def generate_temp_file(input_path, output_path, pid):
    pid_output_path = os.path.join(output_path, pid + '_HT')

    profile_file = os.path.join(input_path, f"{pid}_temp.profile")
    if not os.path.exists(profile_file):
        print(f"File {profile_file} not found.")
        return False

    with open(profile_file, 'r') as f:
        contents = f.readlines()

    contents_for_use = contents[3:]  # skip header lines
    for i in range(10):
        batch = int(len(contents_for_use) / 10)
        file_name = os.path.join(pid_output_path, f"temp{i+1}.txt")
        position = []
        with open(file_name, 'w') as file:
            for j in range(batch):
                line = contents_for_use[i*batch:(i+1)*batch][j]
                if len(line.split()) != 4:
                    continue
                temperature = float(line.split()[-1])
                if temperature >= 250:
                    position.append(j)
            if len(position) < 2:
                continue
            head = position[1]
            tail = position[-1]
            for line in contents_for_use[i*batch:(i+1)*batch][head:tail]:
                file.write(line)

    return True

# Step 3: Heat flux regression
def heat_flux(file_name):
    data = pd.read_csv(file_name, header=None, sep=r"\s+",
                       names=['Step','Temp','E_pair','E_mol','TotEng','Press',
                              'Density','Volume','Lx','Ly','Lz','f_4','f_5'])

    x_length = data['Lx'].iloc[0] * 1e-10
    y_length = data['Ly'].iloc[0] * 1e-10
    z_length = data['Lz'].iloc[0] * 1e-10

    time = data['Step'] * 0.25 * 1e-15
    heat_source = data['f_4']
    heat_sink = data['f_5']
    energy = ((heat_sink - heat_source) / (2 * 23.06)) * 1.6022e-19

    slope, intercept, r_value, p_value, std_err = stats.linregress(time, energy)
    return slope, x_length, y_length, z_length, time, heat_source, heat_sink

# Step 4: Temperature gradient regression
def temp_gradient(file_name):
    with open(file_name, 'r') as f:
        data = f.readlines()

    data = [x.strip().split() for x in data]
    distance = [float(row[1]) for row in data]
    temperature = [float(row[3]) for row in data]

    distance = np.array(distance)
    temperature = np.array(temperature)

    slope_, intercept, r_value, p_value, std_err = stats.linregress(distance, temperature)
    return slope_, temperature, distance

# Step 5: Thermal conductivity
def calculate_TC(slope, slope_, x_length, y_length, z_length):
    return slope / (y_length * z_length * abs(slope_ / x_length))

# Plotting
def plot(time, heat_source, heat_sink, distance, temperature):
    fig, (ax1, ax2) = plt.subplots(figsize=(10, 4), ncols=2)
    ax1.plot(time * 1e10, -heat_source, label='heat source')
    ax1.plot(time * 1e10, heat_sink, label='heat sink')
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Energy (Kcal/mole)')
    ax1.legend(loc='best')

    ax2.scatter(distance, temperature, label='temperature')
    ax2.set_ylabel('Temperature (K)')
    ax2.set_xlabel('Position')
    ax2.legend(loc='best')

    plt.show()

# Main
if __name__ == '__main__':
    num_sets = 8

    if len(sys.argv) != 2:
        print("Usage: python script.py <pid>")
        sys.exit(1)

    pid = sys.argv[1]
    input_path = os.getcwd()
    output_path = os.getcwd()

    hf_success = split_hf_file(input_path, output_path, pid)

    if not hf_success:
        mean_TC, std_TC = 0.00, 0.00
        output_filename = os.path.join(output_path, f"{pid}_TC_results.dat")
        np.savetxt(output_filename, [mean_TC, std_TC], delimiter=' ', fmt='%.6f')
        print(f"Thermal conductivity values saved to {output_filename}")
        sys.exit(0)

    generate_temp_file(input_path, output_path, pid)

    pid_output_path = os.path.join(output_path, f'{pid}_HT')
    print(f"Processing files in {pid_output_path}")
    TCs = []

    start_index = 11 - num_sets
    for i in range(start_index, 11):
        hf_file = os.path.join(pid_output_path, f'hf{i}.txt')
        temp_file = os.path.join(pid_output_path, f'temp{i}.txt')

        slope, x_length, y_length, z_length, time, heat_source, heat_sink = heat_flux(hf_file)
        slope_, temperature, distance = temp_gradient(temp_file)
        TC = calculate_TC(slope, slope_, x_length, y_length, z_length)
        print(f"TC for hf{i}.txt: {TC}")
        TCs.append(TC)

    mean_TC = np.mean(TCs)
    std_TC = np.std(TCs)
    TCs_array = np.array([mean_TC, std_TC])

    print('Thermal conductivity array with mean and std:', TCs_array)

    output_filename = os.path.join(output_path, f"{pid}_TC_result.dat")
    np.savetxt(output_filename, TCs_array, delimiter=' ', fmt='%.6f')
    print(f"Thermal conductivity values saved to {output_filename}")
