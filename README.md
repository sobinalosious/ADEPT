# ADEPT

**A**utomated molecular **D**ynamics **E**ngine for **P**olymer simula**T**ions

ADEPT is an end-to-end, high-throughput workflow that turns a list of polymer SMILES into a rich set of computed properties. Starting from SMILES strings, it builds amorphous polymer structures, runs molecular dynamics (MD) simulations with LAMMPS, computes DFT-level monomer electronic properties with Psi4, and post-processes the trajectories into tabulated polymer properties.

It is designed to run as array jobs on an HPC cluster, processing many polymers in parallel with a single submission.

## Table of Contents

- [Features](#features)
- [Computed Properties](#computed-properties)
- [Repository Structure](#repository-structure)
- [Requirements](#requirements)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Configuration](#configuration)
- [Workflow Overview](#workflow-overview)
- [Output Files](#output-files)
- [Citation](#citation)
- [License](#license)

## Features

- High-throughput generation of amorphous polymer structures directly from SMILES using PySIMM (chain construction, random-walk polymerization, GAFF2 parameters, RESP or Gasteiger charges, and LAMMPS data files).
- Molecular dynamics equilibration and production runs with LAMMPS.
- DFT-level monomer electronic properties using Psi4.
- A single submission script with per-property toggles, so you can compute only the properties you need.
- ML-assisted glass transition temperature (Tg) prediction used as an input guess to MD.
- Automatic organization of inputs and results into a clean directory tree, with per-property CSV summaries.

## Computed Properties

| Category | Properties |
| --- | --- |
| Structural | Density, radius of gyration |
| Monomer electronic (DFT) | HOMO, LUMO, band gap, dipole moment, dipole polarizability, total energy, monomer dielectric constant and refractive index |
| Polymer dielectric (MD) | Dielectric constant (dipole, electronic, and total components), refractive index, permittivity |
| Thermal | Thermal conductivity, glass transition temperature (Tg), specific heat capacity (Cp), thermal diffusivity, volume/linear expansion coefficients |
| Transport | Viscosity, diffusion coefficient, mean square displacement (MSD) |
| Mechanical | Bulk modulus, shear modulus, Young's modulus, Poisson's ratio |

## Repository Structure

```text
ADEPT/
├── 1.AmorphousGeneration/      # Structure generation and DFT
│   ├── APG_resp.py             # Amorphous polymer generation (RESP charges)
│   ├── APG_gas.py              # Amorphous polymer generation (Gasteiger charges)
│   └── electronic_properties.py# Monomer DFT properties (Psi4)
├── 2.Simulations/              # LAMMPS input files
│   ├── lammps_eq1.in           # Equilibration stage 1
│   ├── lammps_eq2.in           # Equilibration stage 2
│   ├── lammps_dc.in            # Dielectric constant
│   ├── lammps_TC.in            # Thermal conductivity (NEMD)
│   ├── lammps_Tg.in            # Glass transition temperature
│   ├── lammps_visc.in          # Viscosity
│   ├── lammps_emd.in           # Thermodynamic/mechanical properties
│   └── modulus_files/          # Elastic modulus LAMMPS templates
├── 3.Analysis/                 # Post-processing and property calculation
│   ├── calc_*.py               # Compute properties from trajectories
│   ├── update_*.py             # Append results to RESULTS/*.csv
│   └── predict_tg.py           # ML Tg predictor
├── Files/
│   ├── files.sh                # Per-PID workflow orchestrator
│   └── model_mlp_tg.pkl        # Pre-trained Tg model
├── SMILES.csv                  # Input: PID, SMILES
├── polymer-workflow.yml        # Conda environment definition
├── submit.sh                   # HPC submission script (toggles + scheduler)
└── README.md
```

## Requirements

- An HPC cluster with the **Sun Grid Engine (SGE)** scheduler (or a compatible scheduler).
- **LAMMPS** installed as an environment module.
- **Conda** for Python environment management.
- **PySIMM** (installed separately, see below).

The Conda environment (`polymer-workflow.yml`) provides Python 3.10, NumPy, pandas, RDKit, Psi4, scikit-learn, joblib, and the RESP charge tooling.

## Installation

### 1. Create the Conda environment

```bash
module load conda
conda env create -f polymer-workflow.yml
conda activate polymer-workflow
```

### 2. Install PySIMM

PySIMM is installed separately inside the same Conda environment. Choose one option.

**Option 1 — pip install**

```bash
conda activate polymer-workflow
git clone https://github.com/polysimtools/pysimm.git
cd pysimm
pip install .
```

Then add the LAMMPS executable path to your `~/.bashrc` (edit the path as needed) and reload it:

```bash
export LAMMPS_EXEC="/usr/bin/lmp"
source ~/.bashrc
```

**Option 2 — complete_install.py**

```bash
conda activate polymer-workflow
git clone https://github.com/polysimtools/pysimm.git
python pysimm/complete_install.py --pysimm "$PWD"
```

### 3. Load LAMMPS

LAMMPS is provided by your cluster's module system and does not need to be installed inside the Conda environment:

```bash
module load lammps
```

## Quick Start

Prepare `SMILES.csv` with a header row followed by one polymer per line:

```text
PID,SMILES
PID00001,*C*
PID00002,*CC(*)C
```

Submit the workflow as an SGE array job. Each task index maps to a row in `SMILES.csv`. For example, to process the first 50 polymers:

```bash
qsub -t 1-50 submit.sh
```

Job logs are written to `Log_Files/`, and results accumulate under `POLYMER_DATA/`.

## Configuration

All run settings live in `submit.sh`. Edit the property toggles to choose which calculations to run (`1` = enable, `0` = disable):

| Toggle | Step | Description |
| --- | --- | --- |
| `DO_APG` | Amorphous generation | Build the amorphous polymer structure (`.lmps`) from SMILES. |
| `DO_OPT` | Equilibration | MD equilibration; computes density and radius of gyration. |
| `DO_MONO_ELECTRONIC` | Monomer DFT | Psi4 monomer electronic properties (HOMO, LUMO, band gap, dipole, etc.). |
| `DO_DC` | Dielectric | Polymer dielectric constant, refractive index, permittivity from MD. |
| `DO_TC` | Thermal conductivity | Thermal conductivity via NEMD. |
| `DO_TG` | Glass transition | Glass transition temperature (uses ML-predicted Tg as input to MD). |
| `DO_VISC` | Viscosity | Viscosity, diffusion coefficient, and MSD. |
| `DO_EMD` | Thermo/mechanical | Cp, thermal/volume/linear expansion, and mechanical moduli. |

Additional settings in `submit.sh`:

- `CHARGE_METHOD` — `"RESP"` or `"GASTEIGER"`. RESP is more accurate (recommended when computing the dielectric constant) but slower than Gasteiger.
- `SMILES_FILE` — input CSV (default `SMILES.csv`).
- `CSV_START_LINE` — `2` to skip the header (default), or `1` if the file has no header.

The scheduler directives at the top of `submit.sh` (queue, parallel environment, module versions) may need to be adjusted to match your cluster.

## Workflow Overview

For each polymer (PID), `Files/files.sh` runs the enabled steps in order:

1. **Amorphous Polymer Generation (APG)** — build the structure and assign force-field parameters and charges.
2. **Optimization / Equilibration** — equilibrate the structure; compute density and radius of gyration.
3. **Monomer Electronic Properties** — run Psi4 DFT on the monomer.
4. **Dielectric Constant** — MD-based dielectric properties.
5. **Thermal Conductivity** — NEMD thermal conductivity.
6. **Glass Transition Temperature** — MD Tg, seeded with an ML-predicted Tg.
7. **Viscosity** — viscosity and diffusion (prefers MD Tg, falls back to ML-predicted Tg).
8. **EMD** — specific heat, expansion coefficients, and mechanical moduli.

## Output Files

Each property step writes per-polymer working files under a dedicated subfolder in `POLYMER_DATA/` (e.g. `POLYMER_DATA/OPTIMIZATION/<PID>/`) and appends a summary row to a CSV in `POLYMER_DATA/RESULTS/`:

```text
POLYMER_DATA/
├── OPTIMIZATION/<PID>/
├── MONOMER_ELECTRONIC/<PID>/
├── DIELECTRIC_CONSTANT/<PID>/
├── THERMAL_CONDUCTIVITY/<PID>/
├── GLASS_TRANSITION/<PID>/
├── VISCOSITY/<PID>/
├── EMD/<PID>/
└── RESULTS/
    ├── RHO_MD.csv      # Density
    ├── RG_MD.csv       # Radius of gyration
    ├── HOMO_DFT.csv    # Monomer electronic properties (HOMO, LUMO, ...)
    ├── DC_MD.csv       # Dielectric constant (+ DCD, DCE, RI, PE)
    ├── TC_MD.csv       # Thermal conductivity
    ├── TG_MD.csv       # Glass transition temperature
    ├── VISC_MD.csv     # Viscosity (+ diffusion coefficient)
    └── CP_MD.csv       # Cp, expansion coefficients, moduli, ...
```

The `RESULTS/` CSVs are keyed by `PID`, so results from many array tasks collect into a single table per property.

## Citation


## License

