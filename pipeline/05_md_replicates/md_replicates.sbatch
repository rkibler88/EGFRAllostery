#!/bin/bash
#SBATCH --job-name=EGFR_100ns_rep
#SBATCH --output=logs/md_rep%a_%A.out
#SBATCH --error=logs/md_rep%a_%A.err
#SBATCH --array=2-5                 # replicate IDs; replicate 1 = published run
#SBATCH --partition=hpg-turin
#SBATCH --account=raquel.dias
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=30gb
#SBATCH --time=48:00:00
#SBATCH --mail-type=END,FAIL        # BEGIN dropped: 4 array tasks would mail 4x
#SBATCH --mail-user=rkibler@ufl.edu
# ---------------------------------------------------------------------------
# EGFR kinase domain -- CORRECTIONS.md Section 9.1, replicate production MD.
#
# Scheduler settings, modules, paths and mdrun flags are copied verbatim from
# the original run_md.sh / ca_extract.sh on HiPerGator. The ONE change is that
# each array task generates its own initial velocities, so the replicates are
# statistically independent. Production parameters (md_rep.mdp) are
# byte-identical to the published md.mdp.
#
# System: 48,489 atoms, AMBER99SB-ILDN / TIP3P, 14,736 SOL, 45 NA, 47 CL,
#         triclinic box 8.881 x 8.881 x 6.280 nm, 310 K, 1 bar.
# Construct: residues 712-970 (259 C-alpha); DCCM is subset to the 253 scored
#            residues 714-966 (0-based C-alpha indices 2:255).
# ---------------------------------------------------------------------------
set -euo pipefail

PROJECT="/blue/raquel.dias/rkibler/EGFR"   # from ca_extract.sh
INPUT_DIR="${PROJECT}"                      # npt.gro npt.cpt topol.top posre.itp index.ndx
WORK_ROOT="${PROJECT}/replicates"

REP=${SLURM_ARRAY_TASK_ID}
SEED=$(( 20250000 + REP ))          # deterministic, distinct, recorded below
RUN="${WORK_ROOT}/rep${REP}"

# --- GROMACS source ---------------------------------------------------------
# gromacs/2025.2 is UNUSABLE: its openmpi/5.0.7 dependency has hard DT_NEEDED
# references to HPC-X libraries that are not installed on hpg-turin --
# libhcoll.so.1, libocoms.so.0 -- so gmx exits 127 before doing anything
# (confirmed on a compute node, job 50623271).
#
# Set GROMACS_SOURCE to whichever route is known to work:
#   "module" -> a site build; set GROMACS_MODULES to a version that executes
#   "conda"  -> conda-forge thread-MPI + CUDA build in CONDA_ENV
GROMACS_SOURCE="conda"
GROMACS_MODULES="cuda/12.8.1 gcc/14.2.0 openmpi/5.0.7 gromacs/2026.0"
CONDA_ENV="gmx2025"   # conda-forge gromacs 2025.2 nompi_cuda + numpy

if [ "${GROMACS_SOURCE}" = "module" ]; then
    # No `module purge`: run_md.sh does not purge, and purging can drop the
    # base environment the module hierarchy relies on.
    module load ${GROMACS_MODULES}
else
    # conda ships its own CUDA runtime, so the cuda module is deliberately NOT
    # loaded -- a system cuda module on LD_LIBRARY_PATH can shadow it.
    # libcuda.so.1 comes from the compute node's NVIDIA driver either way.
    # `set -u` is lifted across activation: conda's shell function dereferences
    # unbound variables (PS1 among them) and aborts under nounset.
    set +u
    module load conda
    conda activate "${CONDA_ENV}"
    set -u
    echo "conda env: ${CONDA_PREFIX:-<unset>}"
fi

# Resolve the GROMACS binary explicitly and fail loudly if it is absent.
# An unresolved `gmx` is exit 127, which is how the first submission died.
GMX="$(command -v gmx || command -v gmx_mpi || true)"
if [ -z "${GMX}" ]; then
    echo "ERROR: neither gmx nor gmx_mpi is on PATH (source=${GROMACS_SOURCE})" >&2
    module list 2>&1 >&2 || true
    echo "--- PATH ---" >&2; echo "${PATH}" | tr ':' '\n' >&2
    exit 1
fi
# Locating the binary is not enough: a missing shared library also exits 127,
# which is how submission 50293979 died. Probe that it actually executes.
if ! GMX_PROBE="$("${GMX}" --version 2>&1)"; then
    echo "ERROR: ${GMX} is on PATH but will not execute:" >&2
    echo "${GMX_PROBE}" >&2
    echo "--- unresolved shared libraries ---" >&2
    ldd "${GMX}" 2>&1 | grep -i 'not found' >&2 || echo "(none reported by ldd)" >&2
    exit 1
fi
# Confirm this build can actually use the GPU; the mdrun flags below require it.
if ! echo "${GMX_PROBE}" | grep -qi 'GPU support:.*CUDA'; then
    echo "ERROR: ${GMX} has no CUDA support -- mdrun uses -nb/-pme/-bonded/-update gpu" >&2
    echo "${GMX_PROBE}" | grep -i -E 'GROMACS version|GPU support|SIMD' >&2
    echo "reinstall with: conda install -n egfr_env -c conda-forge 'gromacs=*=nompi_cuda*'" >&2
    exit 1
fi
echo "using GROMACS binary: ${GMX}"
echo "${GMX_PROBE}" | grep -i -E 'GROMACS version|GPU support' || true

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export GMX_MAXBACKUP=-1             # never silently shuffle files on restart

mkdir -p "${RUN}" "${SLURM_SUBMIT_DIR}/logs"
cd "${RUN}"

# --- stage inputs (copy, never symlink: the run must be self-contained) -----
for f in npt.gro npt.cpt topol.top posre.itp index.ndx; do
    cp -f "${INPUT_DIR}/${f}" .
done
cp -f "${SLURM_SUBMIT_DIR}/eq_seed.mdp" "${SLURM_SUBMIT_DIR}/md_rep.mdp" .
cp -f "${SLURM_SUBMIT_DIR}/dccm_from_covar.py" .
sed "s/__SEED__/${SEED}/g" eq_seed.mdp > eq_rep${REP}.mdp

# --- provenance -------------------------------------------------------------
{
  echo "replicate      ${REP}"
  echo "gen_seed       ${SEED}"
  echo "slurm_job      ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
  echo "node           $(hostname)"
  echo "date_utc       $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "gmx_binary     ${GMX}"
  echo "conda_prefix   ${CONDA_PREFIX:-<unset>}"
  "${GMX}" --version | sed -n '1,12p'
  echo "--- input sha256 ---"
  sha256sum npt.gro npt.cpt topol.top eq_rep${REP}.mdp md_rep.mdp
} > provenance_rep${REP}.txt

# --- stage A: 1 ns re-equilibration, fresh velocities -----------------------
# The published md.mdp has gen_vel = no and restarts from npt.cpt, so rerunning
# it reproduces the same trajectory. This stage is what makes a replicate.
if [ ! -f eq.gro ]; then
    "${GMX}" grompp -f eq_rep${REP}.mdp -c npt.gro -p topol.top -n index.ndx \
               -o eq.tpr -po eq_mdout.mdp -maxwarn 1
    "${GMX}" mdrun -v -deffnm eq \
        -nb gpu \
        -pme gpu \
        -bonded gpu \
        -update gpu \
        -cpt 15
fi

# --- stage B: 100 ns production, flags identical to run_md.sh ---------------
if [ ! -f md_production.tpr ]; then
    "${GMX}" grompp -f md_rep.mdp -c eq.gro -t eq.cpt -p topol.top -n index.ndx \
               -o md_production.tpr -po prod_mdout.mdp
fi

CPI=""
[ -f md_production.cpt ] && CPI="-cpi md_production.cpt"
# -maxh/-cpi added to run_md.sh's flags: leaves a clean checkpoint before the
# 48 h wall so a resubmit of the same array index resumes instead of restarting.
"${GMX}" mdrun -v -deffnm md_production ${CPI} \
    -nb gpu \
    -pme gpu \
    -bonded gpu \
    -update gpu \
    -cpt 15 \
    -maxh 47.5 \
    -noappend

# concatenate parts if the run was resumed (-noappend writes .partNNNN files)
if ls md_production.part*.xtc >/dev/null 2>&1; then
    "${GMX}" trjcat -f md_production.part*.xtc -o md_full.xtc -cat
else
    cp -f md_production.xtc md_full.xtc
fi

# --- analysis ---------------------------------------------------------------
# Group names on stdin, as in ca_extract.sh -- no dependence on index.ndx order.
printf 'Protein\nSystem\n' | "${GMX}" trjconv -s md_production.tpr -f md_full.xtc \
                                         -o md_clean.xtc -pbc mol -center
printf 'C-alpha\nC-alpha\n' | "${GMX}" covar -s md_production.tpr -f md_clean.xtc \
                                        -o eigenval.xvg -ascii covar.ascii -nomwa
printf 'C-alpha\n'          | "${GMX}" rmsf  -s md_production.tpr -f md_clean.xtc \
                                        -o rmsf_rep${REP}.xvg -res
"${GMX}" check -f md_clean.xtc

python3 dccm_from_covar.py --covar covar.ascii --n-ca 259 \
        --traj-first 712 --score-first 714 --score-last 966 \
        --out "${WORK_ROOT}/egfr_dccm_matrix_rep${REP}.csv"
cp -f rmsf_rep${REP}.xvg provenance_rep${REP}.txt "${WORK_ROOT}/"

echo "replicate ${REP} complete"
