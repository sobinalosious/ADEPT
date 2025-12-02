#ADEPT: Automated molecular Dynamics Engine for Polymer simulaTions
# Polymer Workflow: MD + DFT Property Pipeline

This workflow automates the generation of polymer structures from SMILES, performs molecular dynamics (MD) simulations using LAMMPS, and computes a wide range of polymer properties including density, Tg, Cp, dielectric constant, thermal conductivity, viscosity, and mechanical moduli. It also computes DFT-level monomer electronic properties using Psi4.

The workflow is designed for an HPC cluster with:
- Sun Grid Engine (SGE) or compatible scheduler
- LAMMPS installed as a module
- Conda for Python environment management


---

## 1. Create the Conda Environment

Load the conda module on your HPC:

```bash
module load conda
```

Create the environment:

```bash
conda env create -f polymer-workflow.yml
```

Activate it:

```bash
conda activate polymer-workflow
```
---

## 2. Install PySIMM


conda activate polymer-workflow
git clone https://github.com/polysimtools/pysimm.git
cd pysimm
pip install .
Add the following line to your ~/.bashrc file:
export LAMMPS_EXEC="export LAMMPS_EXEC=/usr/bin/lmp"   (change as per your path )
source ~/.bashrc

OR

Inside the active conda environment:

git clone https://github.com/polysimtools/pysimm
python pysimm/complete_install.py --pysimm $PWD

PySIMM handles polymer chain construction, random-walk polymerization, GAFF2 parameter assignment, RESP/Gasteiger charges, and LAMMPS data file generation.


---

## 3. Load LAMMPS on HPC

LAMMPS is loaded using your cluster's module system. For example:

```bash
module load lammps
```

LAMMPS does not need to be installed inside the conda environment.

---


## 5. Using the Submission Script

To run the first 50 polymers:

```bash
qsub -t 1-50 submit.sh
```

Logs are stored in:

```text
Log_Files/
```

---

## 6. Property Toggles

Inside `submit.sh`, you will find:

```bash
DO_APG=
DO_OPT=
DO_MONO_ELECTRONIC=
DO_DC=
DO_TC=
DO_TG=
DO_VISC=
DO_EMD=
```

Set each option to:
- `1` = enable
- `0` = disable

### Meaning of Each Toggle

**DO_APG**

Generates amorphous polymer structures using PySIMM. Builds chains, assigns GAFF2 parameters, RESP/Gasteiger charges, and writes LAMMPS `.lmps` data files.

**DO_OPT**

Runs MD equilibration, calculates density (`RHO_MD.csv`), radius of gyration (`RG_MD.csv`).

**DO_MONO_ELECTRONIC**

Runs Psi4 DFT to compute monomer properties: HOMO, LUMO, bandgap, dipole moment, polarizability, monomer refractive index, monomer dielectric constant, total energy.

**DO_DC**

Computes polymer dielectric properties from MD:
- dipole component
- electronic component
- total dielectric constant
- refractive index
- permittivity

**DO_TC**

Computes thermal conductivity using MD (Green–Kubo or NEMD depending on your LAMMPS inputs).

**DO_TG**

Reads Tg from `TG_EXP.csv` or from your MD-based Tg method and generates `TG_MD.csv`.

**DO_VISC**

Computes viscosity, diffusion coefficient, and MSD.

**DO_EMD**

Computes thermodynamic and mechanical properties:
- specific heat capacity Cp
- thermal diffusivity
- volume expansion coefficient
- linear expansion coefficient
- bulk modulus
- shear modulus
- Young’s modulus
- Poisson’s ratio

---

## 7. Directory Structure (Auto-created)

The workflow will create folders similar to:

```text
POLYMER_DATA/RESULTS/
    RHO_MD.csv
    CP_MD.csv
    DC_MD.csv
    TC_MD.csv
    TG_MD.csv
    VISC_MD.csv

```

Each PID has its own subfolder under `POLYMER_DATA`.

---

