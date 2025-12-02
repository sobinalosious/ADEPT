import sys
import os
import numpy as np

# Constants
epsilon_0 = 8.854187817e-12  # Permittivity of vacuum in C^2/(N·m^2) (F/m)
kB = 1.380649e-23  # Boltzmann constant in J/K
T = 300.0  # Temperature in Kelvin

# Conversion factors for LAMMPS "real" units to SI units
conversion_factor_dipole = 1.602176634e-29  # e·Å to C·m
conversion_factor_volume = 1e-30  # Å³ to m³

def electronic_dielectric_constant(number_of_monomers, V_avg_angstrom3, alpha, replica=6):

    number_density_avg=number_of_monomers*replica/V_avg_angstrom3
    print(number_density_avg)
    print(alpha)
    RHS = (4 * np.pi * number_density_avg * alpha) / 3
    print(RHS)
    if (1 - RHS) == 0:
        raise ValueError("Division by zero encountered in refractive index calculation.")
    n_squared = (1 + 2 * RHS) / (1 - RHS)
    if n_squared < 0:
        raise ValueError("Invalid value encountered in refractive index calculation.")

    refractive_index = np.sqrt(n_squared)
    dielectric_constant = n_squared

    return refractive_index, dielectric_constant

def calculate_dielectric_constant_component(dipole_moments, V_avg_angstrom3):
    """
    Calculate the dielectric constant contribution from dipole fluctuations.

    Args:
        dipole_moments (ndarray): Array of dipole moments.
        V_avg_angstrom3 (float): Average volume in Å³.

    Returns:
        float: Dielectric constant.
    """
    M_squared_avg = np.mean(np.sum(dipole_moments**2, axis=1))
    M_avg_squared = np.sum(np.mean(dipole_moments, axis=0)**2)

    V_avg = V_avg_angstrom3 * conversion_factor_volume

    dielectric_constant = 1 + ((M_squared_avg - M_avg_squared) * (conversion_factor_dipole**2)) / (epsilon_0 * kB * T * V_avg * 3)

    return dielectric_constant

def read_dipole_and_volume_data(filename):
    """
    Read dipole moment and volume data from a file.

    Args:
        filename (str): Path to the data file.

    Returns:
        tuple: (dipole_moments, V_avg_angstrom3)
    """
    try:
        data = np.loadtxt(filename, usecols=(0, 1, 2, 4))
        dipole_moments = data[:, :3]
        volume_angstrom3 = data[:, 3]
        V_avg_angstrom3 = np.mean(volume_angstrom3)

        return dipole_moments, V_avg_angstrom3
    except FileNotFoundError:
        print(f"Error: File {filename} not found.")
        return None, None

def save_results(filename, value):
    """
    Save results to a file.

    Args:
        filename (str): Path to the output file.
        value (float): Value to save.
    """
    np.savetxt(filename, [value], delimiter=' ', fmt='%.6f')

def save_if_missing(variable, output_file, default_value=0.00):
    """
    Save a default value if a variable is missing.

    Args:
        variable: The variable to check.
        output_file (str): Path to the output file.
        default_value (float): Value to save if variable is None.
    """
    if variable is None:
        save_results(output_file, default_value)
        return True
    return False

def main():
    if len(sys.argv) != 2:
        print("Usage: python script.py <pid>")
        sys.exit(1)

    pid = sys.argv[1]
    input_path_alpha = os.path.join(os.getcwd(), "../../MONOMER_ELECTRONIC", pid)
    input_path_dipole = os.getcwd()
    output_path = os.getcwd()
    input_path_model = os.path.join(os.getcwd(), "../../MODEL", pid)
    input_dipole_file = os.path.join(input_path_dipole, f"{pid}_dipole_fluctuation.dat")
    input_number_of_monomers = os.path.join(input_path_model, f"{pid}_number_of_monomers.dat")
    input_alpha_file = os.path.join(input_path_alpha, f"{pid}_avg_polarizability.dat")

    dipole_moments, V_avg_angstrom3 = read_dipole_and_volume_data(input_dipole_file)

    if save_if_missing(dipole_moments, os.path.join(output_path, f"{pid}_dielectric_constant_dipole.dat")) or \
       save_if_missing(V_avg_angstrom3, os.path.join(output_path, f"{pid}_dielectric_constant_dipole.dat")):
        sys.exit(0)

    dielectric_constant_dipole = calculate_dielectric_constant_component(dipole_moments, V_avg_angstrom3)

    try:
        number_of_monomers = np.loadtxt(input_number_of_monomers)
    except FileNotFoundError:
        print(f"Error: File {input_number_of_monomers} not found.")
        sys.exit(1)

    try:
        alpha = np.loadtxt(input_alpha_file)
    except FileNotFoundError:
        print(f"Error: File {input_alpha_file} not found.")
        sys.exit(1)
        
    refractive_index, dielectric_constant_electronic = electronic_dielectric_constant(number_of_monomers,V_avg_angstrom3, alpha)

    dielectric_constant_total = dielectric_constant_electronic + dielectric_constant_dipole - 1

    # Calculate permittivity
    permittivity_total = dielectric_constant_total * epsilon_0 * 1e12 

    print(f"Refractive Index: {refractive_index:.6f}")
    print(f"Dielectric Constant (dipole): {dielectric_constant_dipole:.6f}")
    print(f"Dielectric Constant (electronic): {dielectric_constant_electronic:.6f}")
    print(f"Dielectric Constant (total): {dielectric_constant_total:.6f}")
    print(f"Permittivity (total): {permittivity_total:.6e} pF/m")

    save_results(os.path.join(output_path, f"{pid}_refractive_index_polymer_result.dat"), refractive_index)
    save_results(os.path.join(output_path, f"{pid}_dielectric_constant_dipole_result.dat"), dielectric_constant_dipole)
    save_results(os.path.join(output_path, f"{pid}_dielectric_constant_electronic_polymer_result.dat"), dielectric_constant_electronic)
    save_results(os.path.join(output_path, f"{pid}_dielectric_constant_total_result.dat"), dielectric_constant_total)
    save_results(os.path.join(output_path, f"{pid}_permittivity_total_result.dat"), permittivity_total)

if __name__ == "__main__":
    main()
