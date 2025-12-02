#!/bin/bash

#==============================================================
# Per-PID workflow (sourced by submit.sh)
#==============================================================

# Use initial_dir from environment (set in submit.sh), or fall back to current dir
initial_dir="${initial_dir:-$PWD}"

# ---------------- Helpers ----------------
is_true() {
  case "$(echo "${1:-0}" | tr '[:upper:]' '[:lower:]')" in
    1|y|yes|true|on) return 0 ;; *) return 1 ;;
  esac
}

is_number() {
  [[ "${1:-}" =~ ^-?[0-9]+([.][0-9]+)?$ ]]
}

# Path to Python Tg predictor (RF + MACCS), inside 3.Analysis
PREDICT_TG_SCRIPT="${PREDICT_TG_SCRIPT:-$initial_dir/3.Analysis/predict_tg.py}"

# Return Tg (K) for a PID using RF+MACCS predictor (Python script)
# Expects predictor to print a single numeric Tg on stdout.
get_tg_for_pid() {
  local pid="$1"
  local script="$PREDICT_TG_SCRIPT"

  # Check predictor script exists
  if [[ ! -f "$script" ]]; then
    echo "[WARN] Tg predictor script not found: $script" >&2
    return 0
  fi

  # Run Python predictor; expect a single numeric Tg on stdout
  local val
  val=$(python "$script" "$pid" 2>/dev/null | head -n1 | tr -d '[:space:]')

  # Validate and print Tg (K) with 2 decimals, or return empty on failure
  if [[ -n "$val" ]] && is_number "$val"; then
    awk -v t="$val" 'BEGIN { printf "%.2f\n", t + 0 }'
  else
    echo "[WARN] Tg predictor returned no valid value for PID=$pid" >&2
    return 0
  fi
}

# Return MD Tg (K) for a PID from RESULTS/TG_MD.csv (PID,Tg_MD or PID,Tg)
# Can override path via TG_MD_FILE env; default RESULTS/TG_MD.csv.
get_tg_md_for_pid() {
  local pid="$1"
  local file="${TG_MD_FILE:-RESULTS/TG_MD.csv}"

  # If file doesn't exist, just return (no Tg)
  [[ -f "$file" ]] || { return 0; }

  local val
  val=$(awk -F',' -v pid="$pid" '
    BEGIN { pid_col=-1; tg_col=-1 }
    NR==1 {
      for (i=1;i<=NF;i++){
        h=$i; gsub(/\r/,"",h); tl=tolower(h);
        if (tl=="pid")   pid_col=i;
        if (tl=="tg_md" || tl=="tg") tg_col=i;
      }
      next
    }
    {
      for (i=1;i<=NF;i++){
        gsub(/\r/,"",$i);
        gsub(/^ +| +$/,"",$i);
      }
      if (pid_col>0 && tg_col>0 && $pid_col==pid) {
        print $tg_col;
        exit
      }
    }
  ' "$file")

  # Trim spaces
  val="$(echo -n "$val" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')"

  # If numeric, print as Tg in K with 2 decimals
  if [[ -n "$val" ]] && is_number "$val"; then
    awk -v t="$val" 'BEGIN { printf "%.2f\n", t + 0 }'
  else
    return 0
  fi
}

#==============================================================
# --- Job submission / PID resolution (no nested qsub) ---
#==============================================================

# If PID is already provided (single-job mode), run workflow
if [[ -n "${PID:-}" ]]; then
  echo "[INFO] PID already set (${PID}); running workflow."
else
  # CSV mode: require SGE_TASK_ID and map it to SMILES.csv
  if [[ -z "${SGE_TASK_ID:-}" ]]; then
    echo "[ERROR] No PID and no SGE_TASK_ID found."
    echo "       Submit as an array job, e.g.: qsub -t 1-N submit.sh"
    exit 1
  fi

  line=$((CSV_START_LINE + SGE_TASK_ID - 1))
  PID=$(sed -n "${line}p" "$SMILES_FILE" | cut -d',' -f1 | tr -d '\r')
  echo "[INFO] CSV mode: ${SMILES_FILE} line ${line} → PID=${PID}"
fi

# Safety check
if [[ -z "$PID" ]]; then
  echo "[ERROR] No PID resolved! Exiting."
  exit 1
fi

#==============================================================
# --- Workflow (per PID) with property toggles ---
#==============================================================

# Step 1: Amorphous Polymer Generation (APG)
if is_true "$DO_APG"; then
  method_norm="$(echo "$CHARGE_METHOD" | tr '[:lower:]' '[:upper:]')"
  case "$method_norm" in
    RESP)      apg_script="APG_resp.py" ;;
    GASTEIGER) apg_script="APG_gas.py"  ;;
    *)
      echo "[ERROR] CHARGE_METHOD must be RESP or GASTEIGER (got: '$CHARGE_METHOD')"
      exit 1
      ;;
  esac

  echo "[RUN ] APG for $PID with $method_norm charges ($apg_script)"
  cd 1.AmorphousGeneration
  python "$apg_script" "$PID" "$NSLOTS"
  cd "$initial_dir"
else
  echo "[SKIP] APG disabled (DO_APG=$DO_APG)"
fi

# Step 2: Optimization (always run if DO_OPT=true)
if is_true "$DO_OPT"; then
  echo "[RUN ] Optimization for $PID"
  mkdir -p "POLYMER_DATA/OPTIMIZATION/$PID"
  cd "POLYMER_DATA/OPTIMIZATION/$PID"
  mpirun -n "$NSLOTS" lmp -in "../../../2.Simulations/lammps_eq1.in" -var name "$PID"
  mpirun -n "$NSLOTS" lmp -in "../../../2.Simulations/lammps_eq2.in" -var name "$PID"
  python "../../../3.Analysis/calc_Rg.py" "$PID"
  python "../../../3.Analysis/calc_density.py" "$PID"
  python "../../../3.Analysis/update_Rg_dens.py" "$PID"
  cd "$initial_dir"
else
  echo "[SKIP] Optimization disabled (DO_OPT=$DO_OPT)"
fi


# Step 3: Monomer Electronic Properties
if is_true "$DO_MONO_ELECTRONIC"; then
  echo "[RUN ] Monomer electronic for $PID"
  mkdir -p "POLYMER_DATA/MONOMER_ELECTRONIC/$PID"
  cd "POLYMER_DATA/MONOMER_ELECTRONIC/$PID"
  python "../../../1.AmorphousGeneration/electronic_properties.py" "$PID" "$NSLOTS"
  python "../../../3.Analysis/update_ep.py" "$PID"
  cd "$initial_dir"
else
  echo "[SKIP] Monomer electronic disabled (DO_MONO_ELECTRONIC=$DO_MONO_ELECTRONIC)"
fi

# Step 4: Dielectric Constant
if is_true "$DO_DC"; then
  echo "[RUN ] Dielectric constant for $PID"
  mkdir -p "POLYMER_DATA/DIELECTRIC_CONSTANT/$PID"
  cd "POLYMER_DATA/DIELECTRIC_CONSTANT/$PID"
  mpirun -n "$NSLOTS" lmp -in "../../../2.Simulations/lammps_dc.in" -var name "$PID"
  python "../../../3.Analysis/calc_DC.py" "$PID"
  python "../../../3.Analysis/update_DC.py" "$PID"
  cd "$initial_dir"
else
  echo "[SKIP] Dielectric constant disabled (DO_DC=$DO_DC)"
fi

# Step 5: Thermal Conductivity
if is_true "$DO_TC"; then
  echo "[RUN ] TC for $PID"
  mkdir -p "POLYMER_DATA/THERMAL_CONDUCTIVITY/$PID"
  cd "POLYMER_DATA/THERMAL_CONDUCTIVITY/$PID"
  mpirun -n "$NSLOTS" lmp -in "../../../2.Simulations/lammps_TC.in" -var name "$PID"
  python "../../../3.Analysis/calc_TC.py" "$PID"
  python "../../../3.Analysis/update_TC.py" "$PID"
  cd "$initial_dir"
else
  echo "[SKIP] TC disabled (DO_TC=$DO_TC)"
fi

# Step 6: Glass Transition Temperature (Tg, uses ML-predicted Tg as input guess to MD)
if is_true "$DO_TG"; then
  TG_VAL="$(get_tg_for_pid "$PID")"
  if [[ -z "$TG_VAL" ]]; then
    echo "[SKIP] Tg: no Tg value predicted for PID=$PID → skipping Tg LAMMPS run."
  else
    echo "[RUN ] Tg for $PID using ML Tg=${TG_VAL} K (as input to MD)"
    mkdir -p "POLYMER_DATA/GLASS_TRANSITION/$PID"
    cd "POLYMER_DATA/GLASS_TRANSITION/$PID"
    mpirun -n "$NSLOTS" lmp -in "../../../2.Simulations/lammps_Tg.in" \
          -var name "$PID" -var Tg_input "$TG_VAL"
    python "../../../3.Analysis/calc_Tg.py" "$PID" "$TG_VAL"
    python "../../../3.Analysis/update_Tg.py" "$PID"
    cd "$initial_dir"
  fi
else
  echo "[SKIP] Tg disabled (DO_TG=$DO_TG)"
fi

# Step 7: Viscosity (prefers MD Tg, falls back to ML-predicted Tg)
if is_true "$DO_VISC"; then
  # 1) Try MD Tg from RESULTS/TG_MD.csv
  TG_VAL_MD="$(get_tg_md_for_pid "$PID")"
  TG_VAL=""

  if [[ -n "$TG_VAL_MD" ]]; then
    TG_VAL="$TG_VAL_MD"
    echo "[INFO] Viscosity: using MD Tg from RESULTS/TG_MD.csv for PID=$PID → Tg=${TG_VAL} K"
  else
    # 2) Fallback: ML-predicted Tg via RF+MACCS
    TG_VAL="$(get_tg_for_pid "$PID")"
    if [[ -n "$TG_VAL" ]]; then
      echo "[INFO] Viscosity: MD Tg not found, using ML-predicted Tg for PID=$PID → Tg=${TG_VAL} K"
    fi
  fi

  if [[ -z "$TG_VAL" ]]; then
    echo "[SKIP] Viscosity: no Tg (MD or ML) available for PID=$PID → skipping viscosity LAMMPS run."
  else
    echo "[RUN ] Viscosity for $PID using Tg=${TG_VAL} K"
    mkdir -p "POLYMER_DATA/VISCOSITY/$PID"
    cd "POLYMER_DATA/VISCOSITY/$PID"
    mpirun -n "$NSLOTS" lmp -in "../../../2.Simulations/lammps_visc.in" \
          -var name "$PID" -var Tg "$TG_VAL"
    python "../../../3.Analysis/update_visc.py" "$PID"
    cd "$initial_dir"
  fi
else
  echo "[SKIP] Viscosity disabled (DO_VISC=$DO_VISC)"
fi

# Step 8: EMD (bundle: Cp, alpha, modulus, etc. from an EMD run)
if is_true "$DO_EMD"; then
  echo "[RUN ] EMD for $PID"
  mkdir -p "POLYMER_DATA/EMD/$PID"
  cd "POLYMER_DATA/EMD/$PID"
  mpirun -n "$NSLOTS" lmp -in "../../../2.Simulations/lammps_emd.in" -var name "$PID"

  # Post-process: modulus, Cp, alpha from EMD trajectory
  python "../../../3.Analysis/calc_modulus.py" "$PID"
  python "../../../3.Analysis/update_modulus.py" "$PID"
  python "../../../3.Analysis/calc_Cp_alphaP.py" "$PID"
  python "../../../3.Analysis/update_Cp_alphaP.py" "$PID"
  python "../../../3.Analysis/update_alphaT.py" "$PID"

  cd "$initial_dir"
else
  echo "[SKIP] EMD disabled (DO_EMD=$DO_EMD)"
fi

echo "Workflow complete for PID=$PID"
