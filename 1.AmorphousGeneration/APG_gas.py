#!/usr/bin/env python
# coding: utf-8

import csv
import os
import sys
from rdkit import rdBase, Chem
from rdkit.Chem import AllChem
from pysimm import system, lmps, forcefield
from pysimm.apps.random_walk import random_walk, copolymer

nproc = 6
smiles_csv = "../SMILES.csv"
mol_file_dir = "./test/"
mol_file_gen = True
# Use this variable to take a single gen_pid from the command line
gen_pid = sys.argv[1]  # Command-line argument to capture single gen_pid

script_dir = os.path.dirname(os.path.abspath(__file__))
ch3_mol_path = os.path.join(script_dir, 'CH3.mol')
ter = system.read_mol(ch3_mol_path)

all_gen = False
natoms = 600
replica = 6
density = 0.01
def rewrite_atom_block_clean(filepath):
    with open(filepath, 'r') as f:
        lines = f.readlines()

    counts_line = lines[3]
    natoms = int(counts_line[:3])
    nbonds = int(counts_line[3:6])
    print(f"Cleaning .mol file: Atoms = {natoms}, Bonds = {nbonds}")

    atom_lines = lines[4:4+natoms]
    cleaned = []

    for line in atom_lines:
        fields = line.strip().split()
        if len(fields) < 4:
            raise ValueError(f"Malformed atom line: {line}")
        x, y, z, element = fields[:4]
        extras = fields[4:] + ['0'] * (12 - len(fields[4:]))
        cleaned_line = f"{float(x):>10.4f}{float(y):>10.4f}{float(z):>10.4f} {element:<3}" + ''.join(f"{int(val):3d}" for val in extras[:12]) + '\n'
        cleaned.append(cleaned_line)

    lines[4:4+natoms] = cleaned
    with open(filepath, 'w') as f:
        f.writelines(lines)

    print(f"Finished cleaning {filepath}")

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
    ter.apply_charges(f, charges='gasteiger')

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
    if os.path.isfile(mol_file_dir + file_name) and file_name.endswith(".mol"):
        pid = file_name.rstrip(".mol")  # ✅ define pid early
        try:
            mol_path = mol_file_dir + file_name

            # Clean the .mol file before parsing
            rewrite_atom_block_clean(mol_path)

            # Get head/tail from cleaned file
            tail, head = GetHeadTailAtoms(mol_path)

            pid_list.append(pid)
            head_tail[pid] = {"head": head, "tail": tail}
            print(f"Polymer ID = {pid}\tHead = {head}\tTail = {tail}\n")
        except Exception as e:
            print(f"Polymer ID = {pid} Failed! Reason: {e}")
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
        mol_file_path = os.path.join(script_dir, 'test', f'{pid}.mol')
        polymer = PolymerGen(mol_file_path, head_tail[pid]["head"], head_tail[pid]["tail"],
                             natoms=natoms, density=density, nproc=nproc, debug=True)
        calculate_polymer_properties(mol_file_path, pid)
        AmorphousGen(polymer, pid, replica=replica, density=density, nproc=nproc, debug=True)
    except Exception as e:
        print(e)

    os.chdir(cwd)
