#!/usr/bin/env python
# coding: utf-8

import csv
import os
import sys
from rdkit import rdBase, Chem
from rdkit.Chem import AllChem
from pysimm import system, lmps, forcefield
from pysimm.apps.random_walk import random_walk, copolymer
import time
import numpy as np

mol_file_gen = True
gen_pid = sys.argv[1]  
nproc = int(sys.argv[2])


all_gen = False
natoms = 600
replica = 6
density = 0.01
charge_type = "resp"

ch3_charges = np.array([-0.54422096, 0.18133882, 0.18154373, 0.18133841])
ch3_positions = [
        (-2.5721942e-05, 1.7851831e-05, 4.4830327e-05),
        (-0.16312222344, -1.022431951556, 0.307698571486),
        (0.89613005499, 0.51962506855, 0.305933966401),
        (-0.732701564789, 0.502594324309, -0.614166324907),
    ]


smiles_csv = "../SMILES.csv"
mol_file_dir = "./test/"
script_dir = os.path.dirname(os.path.abspath(__file__))
ch3_mol_path = os.path.join(script_dir, 'CH3.mol')
ter = system.read_mol(ch3_mol_path)




def GenMolFile(file_name, smiles):
    try:
        e = ''
        smiles = smiles.replace('*', '[3H]')
        mol = Chem.MolFromSmiles(smiles)
        mol, e = ETKDG(mol, version=2)
        if e:
            print('Polymer ID = '+str(gen_pid)+'\n'+str(e))
        else:
            Chem.MolToMolFile(mol, file_name, kekulize=False)
    except:
        pass


def ETKDG(mol, version=1):
    mh = Chem.AddHs(mol)
    if version == 1:
        p = AllChem.ETKDG()
    elif version == 2:
        p = AllChem.ETKDGv2()
    else:
        print('invalid input')

    try:
        AllChem.EmbedMolecule(mh, p)
    except Exception as e:
        return [mh, e]
        
    return [mh, '']

def GetHeadTailAtoms(mol_file):
    mol = Chem.MolFromMolFile(mol_file)
    idx_list = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "H" and atom.GetIsotope() == 3:
            idx = atom.GetNeighbors()[0].GetIdx()
            idx+=1
            idx_list.append(idx)

    return idx_list
    
def assign_resp_charges(
    pysimm_obj,
    output_xyz="output.xyz",
    nproc=1,
    opt=True,  # Unified parameter for geometry optimization and PySIMM update
    opt_method='wb97m-d3bj',
    opt_basis='6-31G(d,p)',
    geom_iter=200,
    charge_method='HF',
    charge_basis='6-31G(d)',
    total_charge=0,
    total_multiplicity=1,
    log_name='charge'
):
    """
    Assign RESP charges and optionally optimized geometry to a PySIMM object using Psi4 and the RESP plugin.
    Writes an XYZ file, calculates RESP charges, and assigns charges directly to particles.
    If `opt` is True, geometry optimization is performed and the PySIMM object is updated with optimized coordinates.
    """
    import psi4
    import resp

    psi4.set_num_threads(nproc)
    psi4.set_memory('4 GB')

    try:
        # Step 1: Write XYZ file from PySIMM object
        pysimm_obj.write_xyz(output_xyz)
        print(f"XYZ file written to: {output_xyz}")

        # Step 2: Read XYZ file for Psi4
        with open(output_xyz, 'r') as f:
            mol_str = f.read()
        psi4_mol = psi4.geometry(mol_str)
        psi4_mol.set_molecular_charge(total_charge)
        psi4_mol.set_multiplicity(total_multiplicity)

        # Step 3: Set Psi4 options
        #psi4.set_options({
            #'reference': 'uhf' if total_multiplicity > 1 else 'rhf',  # Use UHF for radicals
            #'basis': charge_basis,
            #'scf_type': 'df',
       # })
        psi4.set_options({
            'reference': 'uhf' if total_multiplicity > 1 else 'rhf',
            'basis': charge_basis,
            'scf_type': 'df',
            'maxiter': 200,
            'diis': True,
            'damping_percentage': 30,
            'guess': 'sad',
        })

        # Step 4: Geometry Optimization and Update (if opt=True)
        if opt:
            psi4.set_options({'geom_maxiter': geom_iter})
            psi4.optimize(f"{opt_method}/{opt_basis}", molecule=psi4_mol)

            # Update PySIMM object with optimized geometry
            optimized_geometry = psi4_mol.save_string_xyz()
            with open(f"optimized_monomer.xyz", "w") as f:
                f.write(optimized_geometry)
            positions = []
            for line in optimized_geometry.splitlines()[1:]:  # Skip first two lines of XYZ format
                parts = line.split()
                if len(parts) == 4:
                    x, y, z = map(float, parts[1:])
                    positions.append((x, y, z))

            # Validate position count
            if len(positions) != len(pysimm_obj.particles):
                raise ValueError("Number of positions in the Psi4 molecule does not match the number of particles.")

            # Assign the new positions to the particles
            for particle, pos in zip(pysimm_obj.particles, positions):
                particle.x, particle.y, particle.z = pos

            
            print(f"Optimized geometry updated in PySIMM object.")
        else:
            print(f"Skipping geometry optimization and PySIMM geometry update.")

        # Step 5: RESP Fitting
        options = {
            'RESP_A': 0.0005,
            'RESP_B': 0.1,
            'VDW_SCALE_FACTORS': [1.4, 1.6, 1.8, 2.0],
            'VDW_POINT_DENSITY': 1.0,
        }
        charges = resp.resp([psi4_mol], options)

        # Step 6: Verify and Assign RESP charges to PySIMM particles
        if len(charges[1]) != len(pysimm_obj.particles):
            raise ValueError(
                f"RESP charges ({len(charges[1])}) do not match PySIMM particles ({len(pysimm_obj.particles)})."
            )

        for i, particle in enumerate(pysimm_obj.particles):
            particle.charge = charges[1][i]

        print(f"RESP charges successfully assigned to PySIMM particles.")
        return pysimm_obj  # Return updated PySIMM object with RESP charges and geometry (if optimized)

    except Exception as e:
        print(f"Error in RESP charge assignment: {e}")
        return None


def PolymerGen(mol_file, head, tail, length=20, natoms=None, density=0.1, nproc=1, debug=False):

    os.environ['OMP_NUM_THREADS'] = str(nproc)

    s = system.read_mol(mol_file)
    if natoms:
        na = s.particles.count
        length = int(natoms/(na - 2) + 0.5)
        print("Polymer length = %d, Num. of atoms = %d" % (length, length*(na - 2)))

    p_head = s.particles[head]
    p_tail = s.particles[tail]
    p_head.linker = 'head'
    p_tail.linker = 'tail'

    for b in p_head.bonds:
        if b.a.elem == 'H' or b.b.elem == 'H':
            pb = b.a if b.b is p_head else b.b 
            s.particles.remove(pb.tag, update=False)
            break

    for b in p_tail.bonds:
        if b.a.elem == 'H' or b.b.elem == 'H':
            pb = b.a if b.b is p_tail else b.b 
            s.particles.remove(pb.tag, update=False)
            break

    s.remove_spare_bonding()

    f = forcefield.Gaff2()
    s.apply_forcefield(f)
    if debug:
        for p in s.particles:
            print(p.type.name)
    s.add_particle_bonding()
    
    if charge_type == "resp":
        # Use RESP charge assignment
        s = assign_resp_charges(pysimm_obj=s,output_xyz="monomer.xyz",nproc=nproc,opt=True,total_multiplicity=1)
    
    else:
        # Use Gasteiger charge assignment
        s.apply_charges(f, charges='gasteiger')
    
    
    
    lmps.quick_min(s, min_style='sd')
    lmps.quick_min(s, min_style='fire')

    print('Building polymer chain by random walk')
    polymer = random_walk(s, nmon=length, forcefield=f, reassign=False)

    # Fix bond order
    for b in polymer.bonds:
        if b.order == None:
            b.order = 1

    # Termination of a polymer chain by -CH3
    # Getting the terminal atoms in the generated polymer
    p_count = polymer.particles.count
    flag = False
    for p in polymer.particles:
        if p.linker == 'tail' and not flag:  # Updated: 'is' replaced with '==' 
            p_tail = p
            flag = True
        elif p.linker == 'head':  # Updated: 'is' replaced with '==' 
            p_head = p
        p.linker = None

    p_head.linker = 'head'
    p_tail.linker = 'tail'
    print("Polymer head = %s\ttail = %s\tcount = %s" % (p_head.tag, p_tail.tag, p_count))

    ter = system.read_mol(ch3_mol_path)
    ter.particles[1].linker = 'tail'
    ter.apply_forcefield(f)
    
    if charge_type == "resp":
        # Use RESP charge assignment
        #ter = assign_resp_charges(pysimm_obj=ter,output_xyz="ter.xyz",nproc=nproc,opt=True,total_multiplicity=2)
        ch3_particles = ter.particles[-4:]
        for i, particle in enumerate(ch3_particles):
            particle.charge = ch3_charges[i]
            particle.x, particle.y, particle.z = ch3_positions[i]
    
    else:
        # Use Gasteiger charge assignment
        ter.apply_charges(f, charges="gasteiger")
    

    # Termination process 1
    print('Terminating polymer chain, process 1')
    c1_polymer = copolymer([polymer, ter], nmon=1, forcefield=f, traj=False)

    # Getting the terminal atoms in the generated polymer
    p_count = c1_polymer.particles.count
    flag = False
    for p in c1_polymer.particles:
        if p.linker == 'tail' and not flag:  # Updated: 'is' replaced with '==' 
            p_tail = p
            flag = True
        elif p.linker == 'head':  # Updated: 'is' replaced with '==' 
            p_head = p
        p.linker = None

    p_head.linker = 'tail'  # Replacement of head to tail
    p_tail.linker = 'head'
    print("Polymer head = %s\ttail = %s\tcount = %s" % (p_head.tag, p_tail.tag, p_count))

    # Termination process 2
    print('Terminating polymer chain, process 2')
    c2_polymer = copolymer([c1_polymer, ter], nmon=1, forcefield=f, traj=False)
    if debug:
        for p in c2_polymer.particles:
            print(p.type.name)

    # Fix linker and bond order
    for p in c2_polymer.particles:
        p.linker = None
    for b in c2_polymer.bonds:
        if b.order == None:
            b.order = 1
    # Re-assignment of forcefield and charge for a terminated polymer
    print("Re-assignment of forcefield and charge")
    c2_polymer.apply_forcefield(f)
    if debug:
        for p in c2_polymer.particles:
            print(p.type.name)
    
    if charge_type == "resp":
        # Use RESP charge assignment
        #c2_polymer = assign_resp_charges(pysimm_obj=c2_polymer,output_xyz="c2_polymer.xyz",nproc=nproc,opt=True,total_multiplicity=1)
        print("Skipping RESP charge assignment")
    
    else:
        # Use Gasteiger charge assignment
        c2_polymer.apply_charges(f, charges='gasteiger')
       
    c2_polymer.set_mm_dist()
    lmps.quick_min(c2_polymer, min_style='sd')
    lmps.quick_min(c2_polymer, min_style='fire')
    c2_polymer.write_lammps('polymer.lmps')
    c2_polymer.write_xyz('polymer.xyz')

    return c2_polymer
def calculate_polymer_properties(mol_file, pid, natoms=600,replica=6):
 
    mol = Chem.MolFromMolFile(mol_file)
    if mol is None:
        raise ValueError("Invalid molecular file or format.")

    # Add explicit hydrogens
    mol = Chem.AddHs(mol)

    # Reset isotopic masses to defaults
    for atom in mol.GetAtoms():
        atom.SetIsotope(0)  # Reset all isotopes

    # Identify head and tail atoms connected to isotopic hydrogens (if present)
    head_tail_atoms = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "H" and atom.GetIsotope() == 3:
            head_tail_atoms.append(atom.GetNeighbors()[0].GetIdx())
    print(f"Head and Tail Atoms: {head_tail_atoms}")

    # Calculate monomer molar mass
    monomer_weight = 0.0
    num_atoms_in_monomer = mol.GetNumAtoms()
    for atom in mol.GetAtoms():
        monomer_weight += atom.GetMass()

    # Adjust for head and tail hydrogens (2 hydrogens removed)
    removed_h_weight = 2 * 1.008  # Hydrogen atomic mass
    adjusted_monomer_weight = monomer_weight - removed_h_weight

    # Calculate the number of monomers in the chain
    n_monomers = int(natoms / (num_atoms_in_monomer - 2) + 0.5)

    # Add termination group weight (CH3 on both ends)
    termination_group_weight = 2 * (12.01 + 3 * 1.008)  # CH3 group on both sides
    total_chain_weight = n_monomers * adjusted_monomer_weight + termination_group_weight
    ampolymer_weight=total_chain_weight*replica
    # Print results
    print(f"Monomer Molar Mass (g/mol): {adjusted_monomer_weight:.2f}")
    print(f"Number of Monomers in Chain: {n_monomers}")
    print(f"Total Chain Molar Mass (g/mol): {total_chain_weight:.2f}")
    print(f"Amorphous Polymer Molar Mass (g/mol): {ampolymer_weight:.2f}")
    # Save each value to separate files
    with open(f"{pid}_monomer_weight.dat", "w") as f:
        f.write(f"{adjusted_monomer_weight:.2f}\n")

    with open(f"{pid}_number_of_monomers.dat", "w") as f:
        f.write(f"{n_monomers}\n")

    with open(f"{pid}_total_chain_weight.dat", "w") as f:
        f.write(f"{total_chain_weight:.2f}\n")

    with open(f"{pid}_ampolymer_weight.dat", "w") as f:
        f.write(f"{ampolymer_weight:.2f}\n")

    print(f"Values saved to separate files with prefix '{pid}'")


def AmorphousGen(polymer, pid, replica=10, density=0.1, nproc=1, debug=False):
    os.environ['OMP_NUM_THREADS'] = str(nproc)

    print('Building amorphous cell of polymer')
    amo_polymer = system.replicate(polymer, replica, density=density, rand=True)
    print('amo gen done')
    amo_polymer.set_mm_dist()
    amo_polymer.write_lammps('amorphous_polymer_{}.lmps'.format(pid))
    amo_polymer.write_xyz('amorphous_polymer_{}.xyz'.format(pid))
    
    return amo_polymer

print(gen_pid)

if mol_file_gen:
    with open(smiles_csv, encoding='ISO-8859-1') as f:  # Updated: Added encoding to prevent UnicodeDecodeError
        reader = csv.reader(f)
        data_found = False
        for row in reader:
            if row[0] == gen_pid:
                file_name = mol_file_dir + row[0] + '.mol'
                smiles = row[1]
                GenMolFile(file_name, smiles)
                data_found = True
                break  # Stop after processing the matching PID
                
        if not data_found:
            print(f"PID {gen_pid} not found in SMILES.csv")

# Proceed with the rest of the script for generating polymers and amorphous structures...

# Determining atom indexes of head and tail atom
head_tail = {}
pid_list = []
file_list = sorted(os.listdir(mol_file_dir))
for file_name in file_list:
    if os.path.isfile(mol_file_dir+file_name) and file_name.endswith(".mol"):
        try:
            tail, head = GetHeadTailAtoms(mol_file_dir+file_name)
            pid = file_name.rstrip(".mol")
            pid_list.append(pid)
            head_tail[pid] = {"head": head, "tail": tail}
            print("Polymer ID = "+str(pid)+"\tHead = "+str(head_tail[pid]["head"])+"\tTail = "+str(head_tail[pid]["tail"])+"\n")
        except:
            print("Polymer ID = "+str(pid)+" Failed!")

if all_gen:
    pid_list = pid_list
else:
    pid_list = [gen_pid]  # Just the single PID you passed as an argument

for pid in pid_list:
    if os.path.isdir(mol_file_dir + pid) and all_gen:
        continue

    print(pid)
    if not os.path.exists(mol_file_dir + pid):
        os.makedirs(mol_file_dir + pid)

    polymer_data_dir = '../POLYMER_DATA/MODEL/' + pid
    if not os.path.exists(polymer_data_dir):
        os.makedirs(polymer_data_dir)

    cwd = os.getcwd()
    os.chdir(polymer_data_dir)

    try:
        start_time = time.time()
        mol_file_path = os.path.join(script_dir, 'test', f'{pid}.mol')
        polymer = PolymerGen(mol_file_path, head_tail[pid]["head"], head_tail[pid]["tail"],
                             natoms=natoms, density=density, nproc=nproc, debug=True)
        calculate_polymer_properties(mol_file_path, pid)
        AmorphousGen(polymer, pid, replica=replica, density=density, nproc=nproc, debug=True)
        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"Total script execution time: {elapsed_time:.2f} seconds")
    except Exception as e:
        print(e)

    os.chdir(cwd)
