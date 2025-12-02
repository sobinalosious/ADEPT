#!/bin/bash
#===============================
# Submission Script
#==============================================================
#$ -S /bin/bash
#$ -pe mpi-48 48
#$ -q hpc
#$ -o Log_Files/$JOB_NAME.o$JOB_ID.$TASK_ID
#$ -e Log_Files/$JOB_NAME.e$JOB_ID.$TASK_ID
#$ -cwd
#==============================================================


# --- Environment Setup ---

module purge
module load conda/24.7.1
conda activate polymer-workflow

module load lammps/29Aug24/intel/24.2
export MKL_INTERFACE_LAYER=""
export OMP_NUM_THREADS=1
#==============================================================
export NSLOTS="$NSLOTS"
initial_dir=$(pwd)

# SMILES CSV (for array jobs)
SMILES_FILE="SMILES.csv"  

#==============================================================
# --- USER CONFIG  ---
#==============================================================
#
# === USAGE EXAMPLE ===
#
# Run from SMILES.csv (array jobs), e.g. first 10 rows (after header):
#    qsub -t 1-10 submit.sh
#
#===============================


#==============================================================
# --- PROPERTY TOGGLES (1=run, 0=skip; 
#==============================================================

DO_APG=1    # Amorphous generation:
#   ✅ Build amorphous polymer structure (.lmps)
	# --- Charge method for APG (choose one) ---
	CHARGE_METHOD="RESP"     # Options: "RESP" or "GASTEIGER" 
								  #(Use RESP if Dielectric Constant is calculating)
								  # Note: RESP is time consuming compared to GASTEIGER

DO_OPT=1    # Optimization / Equilibration:					(Output Files: POLYMER_DATA/RESULTS)
#   ✅ Equilibration runs 									(POLYMER_DATA/OPTIMIZATION/PID/PID_eq1.data, PID_eq2.data)
#   ✅ Density calculation 									(RHO_MD.csv)
#   ✅ Radius of gyration 									(RG_MD.csv)

# Property calculations 
DO_MONO_ELECTRONIC=1  # Monomer electronic properties:(DFT)
#   ✅ HOMO (Highest Occupied Molecular Orbital) 			(HOMO_DFT.csv)
#   ✅ LUMO (Lowest Unoccupied Molecular Orbital) 			(LUMO_DFT.csv)
#   ✅ Band gap energy 										(BANDGAP_DFT.csv)
#   ✅ Dipole moment 										(MU_DFT.csv)
#   ✅ Dipole polarizability 								(ALPHA_DFT.csv)
#   ✅ Total energy 										(ETOTAL_DFT.csv)
#   ✅ Dielectric constant-electronic component (monomer)   (EDCMONO_DFT.csv) (Rough estimate only, use MD data)
#   ✅ Refractive index- (monomer) 							(ERIMONO_DFT.csv) (Rough estimate only, use MD data)

DO_DC=1  # Polymer electronic properties (MD)
#   ✅ Dielectric constant- dipole component (polymer)		(DCD_MD.csv)
#   ✅ Dielectric constant- electronic component (polymer)	(DCE_MD.csv)
#   ✅ Dielectric constant- total (Polymer)					(DC_MD.csv)
#   ✅ Refractive index- (Polymer)							(RI_MD.csv)
#   ✅ Permittivity - (Polymer)								(PE_MD.csv)

DO_TC=1     # Thermal conductivity							(TC_MD.csv)

DO_TG=1     # Glass transition temperature:
#   ✅ Passes Tg from TG_EXP.csv (PID,Tg); 					(TG_MD.csv)

DO_VISC=1   # Viscosity: Passes Tg from TG_EXP.csv (PID,Tg); 
#   ✅ Viscosity (per your LAMMPS input)					(VISC_MD.csv)
#   ✅ Diffusion coefficient								(D_MD.csv)
#   ✅ Mean square displacement (MSD)

DO_EMD=1    
#	✅ Specific heat (Cp)									(CP_MD.csv)
#	✅ Thermal diffusivity (alphaT)							(ALPHAT_MD.csv)
#	✅ Volume expansion coefficient (alphaP)				(ALPHAP_MD.csv)
#	✅ Linear expansion coefficient (alphaL)				(ALPHAL_MD.csv)
#   ✅ Bulk Modulus											(K_MD.csv)
#   ✅ Shear Modulus 										(G_MD.csv)
#   ✅ Youngs Modulus										(E_MD.csv)
#   ✅ Poisson ratio										(NU_MD.csv)
#==============================================================

CSV_START_LINE=2   # 2 = skip header; set to 1 if no header
# Export config so child jobs (manual list mode) inherit these values
export SMILES_FILE CSV_START_LINE
export DO_APG DO_MONO_ELECTRONIC DO_OPT DO_TC DO_TG DO_VISC DO_DC DO_EMD
export CHARGE_METHOD
export initial_dir

#==============================================================
# --- Run (files.sh handles dispatch vs per-PID workflow) ---
#==============================================================
source Files/files.sh
