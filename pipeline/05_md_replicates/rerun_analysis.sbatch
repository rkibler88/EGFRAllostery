#!/bin/bash
#SBATCH --job-name=EGFR_rep_dccm
#SBATCH --output=logs/dccm_rep%a_%A.out
#SBATCH --error=logs/dccm_rep%a_%A.err
#SBATCH --array=2-5
#SBATCH --partition=hpg-turin
#SBATCH --account=raquel.dias
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16gb
#SBATCH --time=02:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=rkibler@ufl.edu
# ---------------------------------------------------------------------------
# Recovery for array 50632117: the 100 ns MD of replicates 2-5 COMPLETED, and
# every task then died at the last line of md_replicates.sbatch with
#   FileNotFoundError: covar.ascii not found.
# so only the analysis tail is redone here. No GPU is requested and mdrun is
# never called: md_production.tpr / md_clean.xtc are reused as they stand.
#
# Submit from the directory that holds md_replicates.sbatch and
# dccm_from_covar.py, and that already contains logs/ -- slurm opens the
# --output path before this script runs:
#     mkdir -p logs && sbatch rerun_analysis.sbatch
# One replicate only:  sbatch --array=3 rerun_analysis.sbatch
# ---------------------------------------------------------------------------
set -euo pipefail

PROJECT="/blue/raquel.dias/rkibler/EGFR"
WORK_ROOT="${PROJECT}/replicates"
REP=${SLURM_ARRAY_TASK_ID}
RUN="${WORK_ROOT}/rep${REP}"

# --- GROMACS (same conda route as the production array) ---------------------
CONDA_ENV="gmx2025"
set +u
module load conda
conda activate "${CONDA_ENV}"
set -u
GMX="$(command -v gmx || command -v gmx_mpi || true)"
[ -n "${GMX}" ] || { echo "ERROR: no gmx on PATH in env ${CONDA_ENV}" >&2; exit 1; }
"${GMX}" --version >/dev/null 2>&1 || { echo "ERROR: ${GMX} will not execute" >&2
                                        ldd "${GMX}" | grep -i 'not found' >&2; exit 1; }
echo "gmx: ${GMX}   conda: ${CONDA_PREFIX:-<unset>}"

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export GMX_MAXBACKUP=-1

[ -d "${RUN}" ] || { echo "ERROR: ${RUN} does not exist" >&2; exit 1; }
cd "${RUN}"
cp -f "${SLURM_SUBMIT_DIR}/dccm_from_covar.py" .

# --- 0. the trajectory must already be there --------------------------------
[ -s md_production.tpr ] || { echo "ERROR: md_production.tpr missing in ${RUN}" >&2; exit 1; }
if [ ! -s md_clean.xtc ]; then
    echo "md_clean.xtc absent -- rebuilding it from the production trajectory"
    if [ ! -s md_full.xtc ]; then
        if ls md_production.part*.xtc >/dev/null 2>&1; then
            "${GMX}" trjcat -f md_production.part*.xtc -o md_full.xtc -cat
        elif [ -s md_production.xtc ]; then
            cp -f md_production.xtc md_full.xtc
        else
            echo "ERROR: no production .xtc in ${RUN} -- the MD output is gone" >&2
            ls -lt | head -20 >&2
            exit 1
        fi
    fi
    printf 'Protein\nSystem\n' | "${GMX}" trjconv -s md_production.tpr -f md_full.xtc \
                                            -o md_clean.xtc -pbc mol -center
fi

# --- 1. find the covariance matrix, whatever gmx called it ------------------
# `-ascii` is a .dat-typed output: GROMACS appends its own extension when the
# name it is handed carries an unrecognised one, so the 100 ns run may have
# left covar.ascii.dat (or the default covar.dat) rather than covar.ascii.
find_ascii() {
    for f in covar.ascii covar.ascii.dat covar.dat; do
        [ -s "$f" ] && { echo "$f"; return 0; }
    done
    return 1
}
if ASCII=$(find_ascii); then
    echo "found existing covariance matrix: ${ASCII} ($(stat -c %s "${ASCII}") bytes)"
else
    echo "no covariance matrix on disk -- re-running gmx covar"
    printf 'C-alpha\nC-alpha\n' | "${GMX}" covar -s md_production.tpr -f md_clean.xtc \
                                            -o eigenval.xvg -ascii covar.ascii -nomwa
    ASCII=$(find_ascii) || {
        echo "ERROR: gmx covar exited 0 but wrote no ascii matrix. Directory:" >&2
        ls -lt | head -20 >&2
        exit 1
    }
    echo "gmx covar wrote: ${ASCII}"
fi

# --- 2. per-residue fluctuations (skipped if the first run got this far) ----
if [ ! -s "rmsf_rep${REP}.xvg" ]; then
    printf 'C-alpha\n' | "${GMX}" rmsf -s md_production.tpr -f md_clean.xtc \
                                   -o "rmsf_rep${REP}.xvg" -res
fi

# --- 3. DCCM ----------------------------------------------------------------
# 259 C-alpha in the construct (712-970); scored subset 714-966 = 253 residues,
# identical to the published matrix.
python3 dccm_from_covar.py --covar "${ASCII}" --n-ca 259 \
        --traj-first 712 --score-first 714 --score-last 966 \
        --out "${WORK_ROOT}/egfr_dccm_matrix_rep${REP}.csv"

cp -f "rmsf_rep${REP}.xvg" "${WORK_ROOT}/"
[ -s "provenance_rep${REP}.txt" ] && cp -f "provenance_rep${REP}.txt" "${WORK_ROOT}/"
echo "replicate ${REP} analysis complete"
