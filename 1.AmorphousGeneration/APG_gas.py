#!/usr/bin/env python
# coding: utf-8

import csv
import os
import sys
from rdkit import rdBase, Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdPartialCharges
from pysimm import system, lmps, forcefield
from pysimm.apps.random_walk import random_walk, copolymer

nproc = 6
smiles_csv = "../SMILES.csv"
mol_file_dir = "./test/"
mol_file_gen = True
charge_method = os.environ.get('APG_CHARGE_METHOD', 'rdkit_gasteiger').strip().lower()
charge_fallback = os.environ.get('APG_CHARGE_FALLBACK', 'fail').strip().lower()
# Use this variable to take a single gen_pid from the command line
gen_pid = sys.argv[1]  # Command-line argument to capture single gen_pid

script_dir = os.path.dirname(os.path.abspath(__file__))
ch3_mol_path = os.path.join(script_dir, 'CH3.mol')
ter = system.read_mol(ch3_mol_path)

all_gen = False
natoms = 600
replica = 6
density = 0.01


def expected_output_files(pid):
    return [
        'polymer.lmps',
        'polymer.xyz',
        f'amorphous_polymer_{pid}.lmps',
        f'amorphous_polymer_{pid}.xyz',
    ]


def write_status_marker(status_dir, pid, success, message=''):
    status_dir = os.path.abspath(status_dir)
    os.makedirs(status_dir, exist_ok=True)
    success_path = os.path.join(status_dir, f'{pid}_SUCCESS')
    fail_path = os.path.join(status_dir, f'{pid}_FAIL')
    target_path = success_path if success else fail_path
    stale_path = fail_path if success else success_path

    if os.path.exists(stale_path):
        os.remove(stale_path)

    with open(target_path, 'w') as f:
        f.write('SUCCESS\n' if success else 'FAIL\n')
        if message:
            f.write(f'{message}\n')

    print(f'Status marker written: {target_path}')


def validate_run_outputs(status_dir, pid):
    status_dir = os.path.abspath(status_dir)
    missing = [
        name for name in expected_output_files(pid)
        if not os.path.exists(os.path.join(status_dir, name))
    ]
    return missing


def log_charge_event(stage, message, logfile='charge_assignment.log'):
    line = f'[{stage}] {message}'
    print(line)
    with open(os.path.abspath(logfile), 'a') as f:
        f.write(line + '\n')


def log_forcefield_event(stage, message, logfile='forcefield_assignment.log'):
    line = f'[{stage}] {message}'
    print(line)
    with open(os.path.abspath(logfile), 'a') as f:
        f.write(line + '\n')


def assign_zero_charges(s, stage):
    for p in s.particles:
        p.charge = 0.0
    log_charge_event(stage, f'Assigned zero charges to {s.particles.count} particles.')


def infer_bond_order_from_r0(elem_a, elem_b, r0):
    if r0 is None:
        return Chem.rdchem.BondType.SINGLE

    a_sym, b_sym = sorted((elem_a, elem_b))

    if a_sym in {'H', 'F', 'Cl', 'Br', 'I'} or b_sym in {'H', 'F', 'Cl', 'Br', 'I'}:
        return Chem.rdchem.BondType.SINGLE
    if (a_sym, b_sym) == ('C', 'C'):
        if r0 >= 1.41:
            return Chem.rdchem.BondType.SINGLE
        if 1.39 <= r0 < 1.41:
            return Chem.rdchem.BondType.AROMATIC
        if 1.33 <= r0 < 1.39:
            return Chem.rdchem.BondType.DOUBLE
        if 1.22 <= r0 < 1.33:
            return Chem.rdchem.BondType.SINGLE
        return Chem.rdchem.BondType.TRIPLE
    if (a_sym, b_sym) == ('C', 'N'):
        if r0 >= 1.34:
            return Chem.rdchem.BondType.SINGLE
        if 1.339 <= r0 < 1.340:
            return Chem.rdchem.BondType.AROMATIC
        if 1.29 <= r0 < 1.339:
            return Chem.rdchem.BondType.SINGLE
        if 1.16 <= r0 < 1.29:
            return Chem.rdchem.BondType.DOUBLE
        return Chem.rdchem.BondType.TRIPLE
    if (a_sym, b_sym) == ('C', 'O'):
        return Chem.rdchem.BondType.SINGLE if r0 >= 1.22 else Chem.rdchem.BondType.DOUBLE
    if (a_sym, b_sym) == ('N', 'O'):
        return Chem.rdchem.BondType.SINGLE if r0 >= 1.25 else Chem.rdchem.BondType.DOUBLE
    if (a_sym, b_sym) == ('N', 'N'):
        if r0 >= 1.30:
            return Chem.rdchem.BondType.SINGLE
        if 1.16 < r0 < 1.30:
            return Chem.rdchem.BondType.DOUBLE
        return Chem.rdchem.BondType.TRIPLE
    if (a_sym, b_sym) == ('O', 'O'):
        return Chem.rdchem.BondType.SINGLE if r0 >= 1.44 else Chem.rdchem.BondType.DOUBLE
    if (a_sym, b_sym) == ('C', 'P'):
        return Chem.rdchem.BondType.SINGLE if r0 >= 1.70 else Chem.rdchem.BondType.DOUBLE
    if (a_sym, b_sym) == ('N', 'P'):
        return Chem.rdchem.BondType.SINGLE if r0 >= 1.65 else Chem.rdchem.BondType.DOUBLE
    if (a_sym, b_sym) == ('O', 'P'):
        return Chem.rdchem.BondType.SINGLE if r0 >= 1.53 else Chem.rdchem.BondType.DOUBLE
    if (a_sym, b_sym) == ('P', 'P'):
        return Chem.rdchem.BondType.SINGLE if r0 >= 1.80 else Chem.rdchem.BondType.DOUBLE
    if (a_sym, b_sym) == ('C', 'S'):
        return Chem.rdchem.BondType.SINGLE if r0 >= 1.64 else Chem.rdchem.BondType.DOUBLE
    if (a_sym, b_sym) == ('N', 'S'):
        return Chem.rdchem.BondType.SINGLE if r0 >= 1.58 else Chem.rdchem.BondType.DOUBLE
    if (a_sym, b_sym) == ('O', 'S'):
        return Chem.rdchem.BondType.SINGLE if r0 >= 1.56 else Chem.rdchem.BondType.DOUBLE
    if (a_sym, b_sym) == ('S', 'S'):
        return Chem.rdchem.BondType.SINGLE
    return Chem.rdchem.BondType.SINGLE


def particle_type_name(p):
    return getattr(getattr(p, 'type', None), 'name', '') or ''


def oxygen_neighbor_count(p):
    return sum(1 for nb in getattr(p, 'bonded_to', []) if getattr(nb, 'elem', '') == 'O')


def is_nitro_like_nitrogen(p):
    if getattr(p, 'elem', '') != 'N':
        return False
    oxygen_neighbors = [nb for nb in getattr(p, 'bonded_to', []) if getattr(nb, 'elem', '') == 'O']
    return len(oxygen_neighbors) >= 2


def aromatic_particle(p):
    type_name = particle_type_name(p)
    if is_nitro_like_nitrogen(p):
        return False
    aromatic_types = {
        'ca', 'cp', 'cq', 'cc', 'cd', 'ce', 'cf', 'cg', 'ch',
        'cu', 'cv', 'nb', 'nc', 'nd', 'ne', 'nf', 'nh', 'no'
    }
    return type_name in aromatic_types


def guess_bond_type(a, b, bond):
    order = getattr(bond, 'order', None)
    if order in (1, 1.0):
        return Chem.rdchem.BondType.SINGLE
    if order in (2, 2.0):
        return Chem.rdchem.BondType.DOUBLE
    if order in (3, 3.0):
        return Chem.rdchem.BondType.TRIPLE
    if order in ('A', 'a', 1.5):
        return Chem.rdchem.BondType.AROMATIC

    if aromatic_particle(a) and aromatic_particle(b):
        return Chem.rdchem.BondType.AROMATIC

    if is_nitro_like_nitrogen(a) and getattr(b, 'elem', '') == 'O':
        return Chem.rdchem.BondType.DOUBLE if oxygen_neighbor_count(b) == 1 and particle_type_name(b) == 'o' else Chem.rdchem.BondType.SINGLE
    if is_nitro_like_nitrogen(b) and getattr(a, 'elem', '') == 'O':
        return Chem.rdchem.BondType.DOUBLE if oxygen_neighbor_count(a) == 1 and particle_type_name(a) == 'o' else Chem.rdchem.BondType.SINGLE

    r0 = getattr(getattr(bond, 'type', None), 'r0', None)
    if r0 is not None:
        r0 = float(r0)
    return infer_bond_order_from_r0(a.elem, b.elem, r0)


def particle_formal_charge(p):
    type_name = getattr(getattr(p, 'type', None), 'name', '') or ''
    elem = getattr(p, 'elem', '')
    bonded_to = list(getattr(p, 'bonded_to', []))

    if elem == 'N' and (type_name in {'n4', 'nx', 'ny', 'nz', 'n+'} or len(bonded_to) >= 4):
        return 1
    if is_nitro_like_nitrogen(p):
        return 1
    if elem == 'O':
        for nb in bonded_to:
            if is_nitro_like_nitrogen(nb):
                return -1 if particle_type_name(p) != 'o' else 0
    return 0


def system_to_rdkit_mol(s):
    rwmol = Chem.RWMol()
    particles = sorted(list(s.particles), key=lambda p: p.tag)
    tag_to_idx = {}

    for p in particles:
        atom = Chem.Atom(p.elem)
        atom.SetNoImplicit(True)
        atom.SetFormalCharge(particle_formal_charge(p))
        idx = rwmol.AddAtom(atom)
        tag_to_idx[p.tag] = idx

    for b in s.bonds:
        if b.a.tag not in tag_to_idx or b.b.tag not in tag_to_idx:
            continue

        bond_type = guess_bond_type(b.a, b.b, b)

        rwmol.AddBond(tag_to_idx[b.a.tag], tag_to_idx[b.b.tag], bond_type)
        if bond_type == Chem.rdchem.BondType.AROMATIC:
            rwmol.GetAtomWithIdx(tag_to_idx[b.a.tag]).SetIsAromatic(True)
            rwmol.GetAtomWithIdx(tag_to_idx[b.b.tag]).SetIsAromatic(True)
            rwmol.GetBondBetweenAtoms(tag_to_idx[b.a.tag], tag_to_idx[b.b.tag]).SetIsAromatic(True)

    mol = rwmol.GetMol()
    sanitize_result = Chem.SanitizeMol(mol, catchErrors=True)
    if sanitize_result == Chem.SanitizeFlags.SANITIZE_NONE:
        return mol, particles

    if sanitize_result == Chem.SanitizeFlags.SANITIZE_KEKULIZE:
        ops = Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE
        retry_result = Chem.SanitizeMol(mol, sanitizeOps=ops, catchErrors=True)
        if retry_result == Chem.SanitizeFlags.SANITIZE_NONE:
            return mol, particles

    if sanitize_result == Chem.SanitizeFlags.SANITIZE_PROPERTIES:
        mol.UpdatePropertyCache(strict=False)
        ops = Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE
        retry_result = Chem.SanitizeMol(mol, sanitizeOps=ops, catchErrors=True)
        if retry_result == Chem.SanitizeFlags.SANITIZE_NONE:
            return mol, particles

    raise ValueError(f'RDKit could not sanitize topology for charge assignment: {sanitize_result}')


def assign_rdkit_gasteiger_charges(s, stage):
    mol, particles = system_to_rdkit_mol(s)
    rdPartialCharges.ComputeGasteigerCharges(mol)

    missing = []
    for atom, particle in zip(mol.GetAtoms(), particles):
        charge = atom.GetProp('_GasteigerCharge') if atom.HasProp('_GasteigerCharge') else None
        if charge is None or charge.lower() in {'nan', '-nan', 'inf', '-inf'}:
            missing.append(particle.tag)
            continue
        particle.charge = float(charge)

    if missing:
        raise ValueError(f'RDKit Gasteiger charge assignment failed for particles {missing[:10]}')

    log_charge_event(stage, f'Assigned RDKit Gasteiger charges to {len(particles)} particles.')


def apply_charges_safe(s, ff, stage):
    method = charge_method or 'gasteiger'
    if method == 'zero':
        assign_zero_charges(s, stage)
        return 'zero'

    try:
        log_charge_event(stage, f"Applying {method} charges.")
        if method == 'rdkit_gasteiger':
            assign_rdkit_gasteiger_charges(s, stage)
        else:
            s.apply_charges(ff, charges=method)
            missing = [p.tag for p in s.particles if getattr(p, 'charge', None) is None]
            if missing:
                raise ValueError(f'Charge assignment left None charges on particles {missing[:10]}')
        return method
    except Exception as e:
        log_charge_event(stage, f'Charge assignment failed: {e}')
        if charge_fallback == 'zero':
            log_charge_event(stage, 'Falling back to zero charges because Gasteiger parameters are missing.')
            assign_zero_charges(s, stage)
            return 'zero'
        raise


GAFF2_SUPPORTED_ELEMENTS = {'H', 'C', 'N', 'O', 'F', 'P', 'S', 'Cl', 'Br', 'I'}


def system_elements(s):
    return sorted({getattr(p, 'elem', None) for p in s.particles if getattr(p, 'elem', None)})


def validate_typed_particles(s, ff_name):
    missing = [p.tag for p in s.particles if getattr(p, 'type', None) is None]
    if missing:
        raise ValueError(f'{ff_name} left untyped particles {missing[:10]}')


def describe_interaction(obj):
    labels = []
    for attr in ('a', 'b', 'c', 'd'):
        atom = getattr(obj, attr, None)
        if atom is None:
            continue
        labels.append(getattr(getattr(atom, 'type', None), 'name', getattr(atom, 'elem', '?')))
    if labels:
        return ','.join(labels)
    return repr(obj)


def validate_forcefield_terms(s, ff_name):
    missing = []
    for collection_name in ('bonds', 'angles', 'dihedrals', 'impropers'):
        collection = getattr(s, collection_name, None)
        if not collection:
            continue
        for obj in collection:
            if getattr(obj, 'type', None) is None:
                missing.append(f'{collection_name[:-1]}:{describe_interaction(obj)}')
                if len(missing) >= 10:
                    break
        if len(missing) >= 10:
            break
    if missing:
        raise ValueError(f'{ff_name} left untyped bonded terms {missing}')


def forcefield_name(ff):
    name = getattr(ff, 'name', None)
    if name:
        return str(name).lower()
    return ff.__class__.__name__.lower()


def instantiate_forcefield(name):
    if name == 'dreiding':
        return forcefield.Dreiding()
    if name == 'gaff2':
        return forcefield.Gaff2()
    raise ValueError(f'Unknown force field {name}')


def apply_existing_forcefield(s, ff, stage):
    ff_name = forcefield_name(ff)
    log_forcefield_event(stage, f'Applying {ff_name} force field.')
    s.apply_forcefield(ff)
    validate_typed_particles(s, ff_name)
    validate_forcefield_terms(s, ff_name)


def choose_forcefield(s, stage):
    elems = system_elements(s)
    unsupported = sorted([elem for elem in elems if elem not in GAFF2_SUPPORTED_ELEMENTS])
    candidates = []
    if unsupported:
        log_forcefield_event(stage, 'Using dreiding because GAFF2 does not cover elements: ' + ', '.join(unsupported))
        candidates = ['dreiding']
    else:
        candidates = ['gaff2', 'dreiding']

    last_error = None
    for name in candidates:
        ff = instantiate_forcefield(name)
        try:
            apply_existing_forcefield(s, ff, stage)
            return ff
        except Exception as exc:
            last_error = exc
            log_forcefield_event(stage, f'{name} force field assignment failed: {exc}')

    raise RuntimeError(f'No usable force field found for {stage}: {last_error}')


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
    mol = Chem.MolFromMolFile(mol_file, removeHs=False)
    if mol is None:
        raise ValueError(f'Could not read MOL file {mol_file}')
    idx_list = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "H" and atom.GetIsotope() == 3:
            idx = atom.GetNeighbors()[0].GetIdx()
            idx+=1
            idx_list.append(idx)
    if len(idx_list) != 2:
        raise ValueError(f'Expected 2 labeled head/tail atoms, found {len(idx_list)}')
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

    f = choose_forcefield(s, 'monomer')
    if debug:
        for p in s.particles:
            print(getattr(getattr(p, 'type', None), 'name', None))
    s.add_particle_bonding()
    apply_charges_safe(s, f, 'monomer')
    lmps.quick_min(s, min_style='sd')

    lmps.quick_min(s, min_style='fire')

    print('Building polymer chain by random walk')
    rw_reassign = forcefield_name(f) != 'gaff2'
    try:
        polymer = random_walk(s, nmon=length, forcefield=f, reassign=rw_reassign)
    except Exception as exc:
        if not rw_reassign:
            log_forcefield_event('random_walk', f'Initial random_walk failed: {exc}. Retrying with reassign=True.')
            polymer = random_walk(s, nmon=length, forcefield=f, reassign=True)
        else:
            raise

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
    apply_existing_forcefield(ter, f, 'terminator')
    apply_charges_safe(ter, f, 'terminator')

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
    apply_existing_forcefield(c2_polymer, f, 'terminated_polymer')
    if debug:
        for p in c2_polymer.particles:
            print(getattr(getattr(p, 'type', None), 'name', None))
    apply_charges_safe(c2_polymer, f, 'terminated_polymer')

    c2_polymer.set_mm_dist()
    lmps.quick_min(c2_polymer, min_style='sd')
    lmps.quick_min(c2_polymer, min_style='fire')
    c2_polymer.write_lammps('polymer.lmps')
    c2_polymer.write_xyz('polymer.xyz')

    return c2_polymer


def calculate_polymer_properties(mol_file, pid, natoms=600,replica=6):
 
    mol = Chem.MolFromMolFile(mol_file, removeHs=False)
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
if all_gen:
    file_list = sorted(os.listdir(mol_file_dir))
else:
    file_list = [f'{gen_pid}.mol']
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

    polymer_data_dir = os.path.abspath('../POLYMER_DATA/MODEL/' + pid)
    if not os.path.exists(polymer_data_dir):
        os.makedirs(polymer_data_dir)

    cwd = os.getcwd()
    os.chdir(polymer_data_dir)

    try:
        mol_file_path = os.path.join(script_dir, 'test', f'{pid}.mol')
        if pid not in head_tail:
            raise KeyError(f'No head/tail labels found for {pid}')
        polymer = PolymerGen(mol_file_path, head_tail[pid]["head"], head_tail[pid]["tail"],
                             natoms=natoms, density=density, nproc=nproc, debug=True)
        calculate_polymer_properties(mol_file_path, pid)
        AmorphousGen(polymer, pid, replica=replica, density=density, nproc=nproc, debug=True)
        missing_outputs = validate_run_outputs(polymer_data_dir, pid)
        if missing_outputs:
            write_status_marker(
                polymer_data_dir,
                pid,
                False,
                'Missing expected output files: ' + ', '.join(missing_outputs)
            )
        else:
            write_status_marker(polymer_data_dir, pid, True)
    except Exception as e:
        print(e)
        write_status_marker(polymer_data_dir, pid, False, str(e))
    finally:
        os.chdir(cwd)
