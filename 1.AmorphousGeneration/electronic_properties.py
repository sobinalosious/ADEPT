import os
import csv
import sys
import psi4
import time
from rdkit import Chem
from rdkit.Chem import AllChem


def generate_molecule(smiles):
    """
    Generates a 3D molecule from a SMILES string using the ETKDG method.
    Replaces attachment points (*) with tritium (H[3]).
    """
    smiles = smiles.replace("*", "[3H]")  # Replace attachment points
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string: {smiles}")
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv2()
    AllChem.EmbedMolecule(mol, params)
    return mol


def generate_xyz_string(mol):
    """
    Generates the XYZ format string for an RDKit molecule suitable for Psi4.
    
    Args:
        mol: An RDKit molecule object with 3D coordinates.
    
    Returns:
        xyz_content (str): The molecule in XYZ file format as a string.
    """
    conf = mol.GetConformer()
    atoms = mol.GetAtoms()
    lines = [f"{len(atoms)}", ""]  # First line: atom count, second line: blank

    for i, atom in enumerate(atoms):
        pos = conf.GetAtomPosition(i)
        lines.append(f"{atom.GetSymbol()} {pos.x:.6f} {pos.y:.6f} {pos.z:.6f}")

    xyz_content = "\n".join(lines)
    return xyz_content


def geometry_optimization(xyz_file, num_processors=1):
    """
    Optimize molecular geometry using Psi4.

    Args:
        xyz_file (str): Path to the input XYZ file.
        num_processors (int): Number of processors to use (default: 1).

    Returns:
        optimized_geometry (str): The optimized geometry in XYZ format.
    """
    import psi4
    import time

    # Load molecule from XYZ file
    with open(xyz_file, "r") as f:
        mol = psi4.geometry(f.read())

    # Set memory and number of threads
    psi4.set_memory("4 GB")
    psi4.set_num_threads(num_processors)

    print("Performing geometry optimization...")
    start_optimization = time.time()

    # Step 1: Initial Semi-Empirical Optimization
    print("Performing initial semi-empirical optimization...")
    try:
        psi4.set_options({'geom_maxiter': 30, 'maxiter': 50})
        psi4.optimize('hf/sto-3g', molecule=mol)
        print("Initial semi-empirical optimization converged successfully.")
    except psi4.OptimizationConvergenceError:
        print("WARNING: Initial semi-empirical optimization did not converge. Using the latest geometry.")

    # Step 2: Hartree-Fock Optimization
    print("Performing Hartree-Fock optimization...")
    try:
        psi4.set_options({'geom_maxiter': 50, 'maxiter': 100})
        psi4.optimize('hf/6-31G', molecule=mol)
        print("Hartree-Fock optimization converged successfully.")
    except psi4.OptimizationConvergenceError:
        print("WARNING: Hartree-Fock optimization did not converge. Using the latest geometry.")

    # Step 3: Final DFT Refinement
    print("Performing final DFT refinement...")
    try:
        psi4.set_options({'geom_maxiter': 100, 'maxiter': 150})
        psi4.optimize('wb97m-d3bj/6-311+G(2d,p)', molecule=mol)
        print("Final DFT refinement converged successfully.")
    except psi4.OptimizationConvergenceError:
        print("WARNING: Final DFT refinement did not converge. Using the latest geometry.")

    optimized_geometry = mol.save_string_xyz()  # Get optimized geometry as XYZ formatted string
    print(f"Geometry optimization completed in {time.time() - start_optimization:.2f} seconds.")
    return mol, optimized_geometry


def dipole_polarizability(mol_input, delta=1e-4, nproc=1):
    import numpy as np
    import psi4
    import datetime
    from packaging.version import parse as version_parse

    # Unit conversion and constants
    atomic_conv = 1.648777e-41
    conversion_factor = (atomic_conv * (1e10)**3) / (4 * np.pi * 8.854187817e-12)
    au_to_debye = 2.541746

    # Initialize arrays and error flag
    tensor_diff = np.zeros((3, 3))
    dipole_results = np.zeros((2, 3, 3))
    calculation_error = False

    # Create Psi4 molecule from XYZ string
    #molecule = psi4.geometry(geom_xyz_str)
    molecule = mol_input

    # Determine appropriate basis set based on atoms present
    atoms_in_mol = {molecule.symbol(i) for i in range(molecule.natom())}
    if "I" in atoms_in_mol:
        selected_basis = "LanL2DZ"
    elif "Br" in atoms_in_mol:
        selected_basis = "6-311G(d,p)"
    else:
        selected_basis = "6-311+G(2d,p)"

    # Configure Psi4 settings
    psi4.set_memory("4 GB")
    psi4.set_num_threads(nproc)
    psi4.set_options({
        "basis": selected_basis,
        "reference": "rhf",
        "scf_type": "df",
        "e_convergence": 1e-6,
        "d_convergence": 1e-6,
    })

    print(f"Using basis set: {selected_basis}")
    print("Commencing finite field polarizability calculation using Psi4...")
    start_time = datetime.datetime.now()

    def compute_dipole(perturb_val, axis_label, mol_instance):
        try:
            field_vector = [0.0, 0.0, 0.0]
            axis_index = "xyz".index(axis_label.lower())
            field_vector[axis_index] = perturb_val
            psi4.set_options({
                "perturb_h": True,
                "perturb_with": "dipole",
                "perturb_dipole": field_vector,
            })
            energy_val, wavefunction = psi4.energy('wb97m-d3bj/' + selected_basis, molecule=mol_instance, return_wfn=True)
            psi4.oeprop(wavefunction, "DIPOLE")
            if version_parse(psi4.__version__) < version_parse('1.3.100'):
                dipole_values = np.array([
                    psi4.variable("SCF DIPOLE X") / au_to_debye,
                    psi4.variable("SCF DIPOLE Y") / au_to_debye,
                    psi4.variable("SCF DIPOLE Z") / au_to_debye,
                ])
            else:
                dipole_values = np.array(psi4.variable("SCF DIPOLE"))
            return dipole_values, False
        except psi4.SCFConvergenceError:
            field_vector[axis_index] *= 0.5
            try:
                psi4.set_options({
                    "perturb_h": True,
                    "perturb_with": "dipole",
                    "perturb_dipole": field_vector,
                })
                energy_val, wavefunction = psi4.energy('wb97m-d3bj/' + selected_basis, molecule=mol_instance, return_wfn=True)
                psi4.oeprop(wavefunction, "DIPOLE")
                if version_parse(psi4.__version__) < version_parse('1.3.100'):
                    dipole_values = np.array([
                        psi4.variable("SCF DIPOLE X") / au_to_debye,
                        psi4.variable("SCF DIPOLE Y") / au_to_debye,
                        psi4.variable("SCF DIPOLE Z") / au_to_debye,
                    ])
                else:
                    dipole_values = np.array(psi4.variable("SCF DIPOLE"))
                return dipole_values, False
            except Exception as exc:
                print(f"Failed calculation on axis {axis_label} with error: {exc}")
                return np.array([np.nan, np.nan, np.nan]), True
        except Exception as err:
            print(f"Unexpected error for axis {axis_label}: {err}")
            return np.array([np.nan, np.nan, np.nan]), True

    # Prepare perturbation arguments for each axis direction and sign
    perturb_args = [(delta, ax, molecule.clone()) for ax in "xyz"] + \
                   [(-delta, ax, molecule.clone()) for ax in "xyz"]

    # Perform finite field calculations sequentially
    computed_results = [compute_dipole(val, ax, mol_inst) for val, ax, mol_inst in perturb_args]

    # Aggregate results and check for errors
    for idx, (dip_outcome, err_flag) in enumerate(computed_results):
        group = idx // 3
        component = idx % 3
        dipole_results[group, component] = dip_outcome
        if err_flag:
            calculation_error = True

    # Compute the polarizability tensor and its average
    tensor_diff = -(dipole_results[0] - dipole_results[1]) / (2 * delta) * conversion_factor
    avg_polarizability = np.mean(np.diag(tensor_diff))

    end_time = datetime.datetime.now()
    if calculation_error:
        print("Some errors occurred during the polarizability calculations.")
    else:
        print(f"Calculations completed normally in {end_time - start_time}.")

    return avg_polarizability, tensor_diff

def electronic_properties(mol_input, nproc=1):
    """
    Calculates electronic properties using Psi4 from an optimized geometry string.
    
    The total energy, HOMO, LUMO, and dipole moment are calculated with a single-point
    calculation using the ωB97M-D3BJ functional. The 6–311G(d,p) basis set is used for 
    H, C, N, O, F, P, S, Cl, and Br atoms, while the LanL2DZ basis set is used for I atoms.
    
    Args:
        optimized_geom_str (str): Optimized molecular geometry in XYZ format.
        nproc (int): Number of processors to use.
    
    Returns:
        dict: A dictionary containing total energy (kcal/mol), HOMO and LUMO (eV), 
              and dipole moment (Debye).
    """
    import psi4

    psi4.core.clean()  # Clean any leftover Psi4 resources
    psi4.set_num_threads(nproc)  # Set number of processors dynamically
    psi4.set_memory('4 GB')

    # Create Psi4 molecule from the optimized geometry string
    #psi4_mol = psi4.geometry(optimized_geom_str)
    psi4_mol = mol_input

    # Set options for single-point energy calculation using ωB97M-D3BJ functional 
    # and the 6-311G(d,p) basis set by default.
    psi4.set_options({
        'basis': '6-311G(d,p)',  # Basis set for H, C, N, O, F, P, S, Cl, Br atoms
        'scf_type': 'df',        # Density fitting for faster calculations
        'reference': 'rhf',      # Restricted Hartree-Fock
        'd_convergence': 1e-6,   # Tight SCF convergence
        'e_convergence': 1e-6,   # Tight energy convergence
    })

    # Adjust basis set for iodine if present
    if "I" in [psi4_mol.symbol(i) for i in range(psi4_mol.natom())]:
        psi4.set_basis("LanL2DZ", "I")

    # Perform single-point energy calculation with ωB97M-D3BJ functional
    energy, wavefunction = psi4.energy('wb97m-d3bj/6-311G(d,p)', molecule=psi4_mol, return_wfn=True)

    # Convert energy to kcal/mol
    total_energy_kcal = energy * 627.509

    # Extract HOMO and LUMO from molecular orbital energies
    mo_energies = wavefunction.epsilon_a().np  # Alpha molecular orbital energies
    homo = mo_energies[wavefunction.nalpha() - 1]  # HOMO is the last occupied orbital
    lumo = mo_energies[wavefunction.nalpha()]      # LUMO is the first unoccupied orbital

    # Convert HOMO and LUMO to eV
    homo_eV = homo * 27.2114
    lumo_eV = lumo * 27.2114

    # Calculate dipole moment
    dipole_components = wavefunction.variable('SCF DIPOLE')  # Returns [X, Y, Z] components
    dipole_moment_au = sum(d**2 for d in dipole_components)**0.5  # Magnitude in atomic units

    # Convert dipole moment to Debye
    dipole_moment_debye = dipole_moment_au * 2.541746

    return {
        'Total Energy (kcal/mol)': total_energy_kcal,
        'HOMO (eV)': homo_eV,
        'LUMO (eV)': lumo_eV,
        'Dipole Moment (Debye)': dipole_moment_debye,
    }

def refractive_index_with_density(psi4_mol, density, alpha):
    """
    Estimate refractive index and dielectric constant using given density, 
    average polarizability, and molecular mass from psi4_mol.
    
    Args:
        psi4_mol: A Psi4 molecule object.
        density (float): Density in g/cm³.
        alpha (float): Average polarizability in Å³.
    
    Returns:
        tuple: (refractive_index (unitless), dielectric_constant)
    """
    import numpy as np

    # Expanded atomic weights for common elements
    atomic_weights = {
        'H': 1.008,  'C': 12.011, 'N': 14.007, 'O': 15.999,
        'F': 18.998, 'P': 30.974, 'S': 32.06,  'Cl': 35.45,
        'Br': 79.904, 'I': 126.90, 'Si': 28.085, 'B': 10.81,
        'Na': 22.990, 'Mg': 24.305, 'Al': 26.982, 'K': 39.098,
        'Ca': 40.078, 'Fe': 55.845, 'Cu': 63.546, 'Zn': 65.38
        # Add more elements as needed
    }
    
    natoms = psi4_mol.natom()
    missing_atoms = set()
    mass_u = 0.0
    for i in range(natoms):
        symbol = psi4_mol.symbol(i)
        weight = atomic_weights.get(symbol)
        if weight is None:
            missing_atoms.add(symbol)
            weight = 0.0
        mass_u += weight

    if missing_atoms:
        print(f"Warning: Atomic weights for {', '.join(missing_atoms)} are not defined.")
    
    # Molecular mass in g/mol
    molecular_mass = mass_u

    # Avogadro's number
    NA = 6.022e23

    # Calculate number density N (molecules/cm³) using given density
    N = (density * NA) / molecular_mass

    # Convert alpha from Å³ to cm³
    alpha_cm3 = alpha * 1e-24

    # Apply Lorentz-Lorenz equation
    RHS = (4 * np.pi * N * alpha_cm3) / 3
    if (1 - RHS) == 0:
        raise ValueError("Division by zero encountered in refractive index calculation.")
    n_squared = (1 + 2*RHS) / (1 - RHS)
    if n_squared < 0:
        raise ValueError("Invalid value encountered in refractive index calculation.")

    refractive_index = np.sqrt(n_squared)
    dielectric_constant = n_squared  # dielectric constant equals n²

    return refractive_index, dielectric_constant



def refractive_index(psi4_mol, alpha):
    """
    Estimate refractive index and dielectric constant using a simple bounding box 
    volume estimation and the Lorentz-Lorenz equation.
    
    Args:
        psi4_mol: A Psi4 molecule object.
        alpha (float): Average polarizability in Å³.
    
    Returns:
        tuple: (refractive_index (unitless), dielectric_constant)
    """
    import numpy as np

    natoms = psi4_mol.natom()
    coords_list = []
    for i in range(natoms):
        vec = psi4_mol.xyz(i)
        coords_list.append([float(vec[0]), float(vec[1]), float(vec[2])])
    
    coords = np.array(coords_list)

    # Calculate bounding box volume in Å³
    min_coords = coords.min(axis=0)
    max_coords = coords.max(axis=0)
    volume_angstrom_cubed = np.prod(max_coords - min_coords)
    
    # Apply Lorentz-Lorenz equation simplified:
    # Since N = 1/volume and alpha in Å³, RHS simplifies to 4*pi*alpha/(3*volume)
    RHS = (4 * np.pi * alpha) / (3 * volume_angstrom_cubed)

    if (1 - RHS) == 0:
        raise ValueError("Division by zero encountered in refractive index calculation.")
    n_squared = (1 + 2*RHS) / (1 - RHS)
    if n_squared < 0:
        raise ValueError("Invalid value encountered in refractive index calculation.")

    refractive_index = np.sqrt(n_squared)
    dielectric_constant = n_squared  # Dielectric constant is n²

    return refractive_index, dielectric_constant




def process_smiles(pid, smiles, output_dir, nproc=1):
    """
    Processes a single SMILES string: generates the molecule, saves the initial 
    geometry as an XYZ file, optimizes the geometry, calculates dipole polarizability 
    and electronic properties, and saves results to files.
    """
    try:
        # Generate molecule and save initial geometry
        mol = generate_molecule(smiles)
        initial_xyz_str = generate_xyz_string(mol)
        initial_xyz_file = os.path.join(output_dir, f"{pid}_initial.xyz")
        with open(initial_xyz_file, "w") as f:
            f.write(initial_xyz_str)
        print(f"Saved initial geometry to {initial_xyz_file}")

        # Optimize geometry using the modified function
        optimized_mol, optimized_geom_str = geometry_optimization(initial_xyz_file, num_processors=nproc)
        print(f"Optimized geometry obtained for PID {pid}.")

        # Save optimized geometry to file
        optimized_xyz_file = os.path.join(output_dir, f"{pid}_optimized.xyz")
        with open(optimized_xyz_file, "w") as f:
            f.write(optimized_geom_str)
        print(f"Optimized geometry saved to {optimized_xyz_file}")

        # Calculate dipole polarizability
        alpha, polarizability_tensor = dipole_polarizability(optimized_mol, nproc=nproc)
        print(f"PID {pid} - Average polarizability: {alpha:.3f} angstrom^3")
        print(f"PID {pid} - Polarizability tensor:\n{polarizability_tensor}")

        # Save polarizability results
        with open(os.path.join(output_dir, f"{pid}_avg_polarizability.dat"), "w") as f:
            f.write(f"{alpha:.3f}\n")
        with open(os.path.join(output_dir, f"{pid}_polarizability_tensor.dat"), "w") as f:
            f.write(f"{polarizability_tensor}\n")

        # Calculate refractive index and dielectric constant
        refractive_ind, dielectric_constant = refractive_index(optimized_mol, alpha)
        print(f" Refractive Index: {refractive_ind:.3f}, Dielectric Constant: {dielectric_constant:.3f}")

        # Save refractive index and dielectric constant
        with open(os.path.join(output_dir, f"{pid}_refractive_index_monomer.dat"), "w") as f:
            f.write(f"{refractive_ind}\n")
        with open(os.path.join(output_dir, f"{pid}_electronic_dielectric_constant_monomer.dat"), "w") as f:
            f.write(f"{dielectric_constant}\n")

        # Calculate electronic properties
        properties = electronic_properties(optimized_mol, nproc=nproc)
        print(f"Electronic properties for PID {pid}:")
        for key, value in properties.items():
            print(f"  {key}: {value}")

        # Save electronic properties
        property_files = {
            "total_energy.dat": "Total Energy (kcal/mol)",
            "homo.dat": "HOMO (eV)",
            "lumo.dat": "LUMO (eV)",
            "dipole_moment.dat": "Dipole Moment (Debye)"
        }
        for file_name, property_key in property_files.items():
            with open(os.path.join(output_dir, f"{pid}_{file_name}"), "w") as f:
                f.write(f"{properties[property_key]}\n")

    except Exception as e:
        print(f"Error processing PID {pid}: {e}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <PID> <NPROC>")
        sys.exit(1)

    pid = sys.argv[1]
    nproc = int(sys.argv[2])
    smiles_csv = "../../../SMILES.csv"

    output_dir = os.getcwd() 
    #os.makedirs(output_dir, exist_ok=True)

    start_time = time.time()

    # Load SMILES from CSV
    smiles = None
    try:
        with open(smiles_csv, "r") as f:
            reader = csv.reader(f)
            for row in reader:
                if row[0] == pid:
                    smiles = row[1]
                    break
        if smiles is None:
            raise ValueError(f"PID {pid} not found in {smiles_csv}")

        print(f"Processing PID: {pid} with SMILES: {smiles}")
        process_smiles(pid, smiles, output_dir, nproc=nproc)

    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Total script execution time: {elapsed_time:.2f} seconds")

