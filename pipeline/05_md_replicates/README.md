# Replicate production MD — CORRECTIONS.md §9.1

Four independent 100 ns replicates of the EGFR kinase-domain monomer, to put
§2.7 (DCCM) and §2.8 (RMSF hierarchy) on more than one trajectory.

## Cluster configuration — copied from the original scripts

Taken verbatim from `run_md.sh` and `ca_extract.sh` (UF HiPerGator):

| Setting | Value | Source |
|---|---|---|
| partition | `hpg-turin` | `run_md.sh` |
| account | `raquel.dias` | `run_md.sh` |
| resources | 1 node, 1 task, 8 CPUs, 1 GPU, 30 gb, 48 h | `run_md.sh` |
| mail | `rkibler@ufl.edu` | `run_md.sh` |
| modules | `conda` only — see GROMACS note below | — |
| project dir | `/blue/raquel.dias/rkibler/EGFR` | `ca_extract.sh` |
| mdrun flags | `-nb gpu -pme gpu -bonded gpu -update gpu -cpt 15` | `run_md.sh` |

Replicates run in `/blue/raquel.dias/rkibler/EGFR/replicates/repN`, reading
inputs from the project directory and writing nothing back into it.

Two deliberate departures from `run_md.sh`:

- `--mail-type=END,FAIL` rather than `BEGIN,END,FAIL`, since a 4-task array
  would otherwise send four BEGIN mails.
- `-maxh 47.5 -noappend` added alongside `-cpi`, so a task that reaches the
  48 h wall leaves a clean checkpoint and a resubmit of the same array index
  resumes rather than restarting. `gmx trjcat` stitches the parts.

## GROMACS source — a departure from the original, and it is methods-relevant

The `gromacs/2025.2` module **cannot be used**. Its `openmpi/5.0.7` dependency
carries hard `DT_NEEDED` references to HPC-X libraries that are not installed:

    libhcoll.so.1 => not found
    libocoms.so.0 => not found

`gmx` therefore fails to load and exits 127. Confirmed on an `hpg-turin`
compute node (job 50623271), not only on the login node, and neither library
appears in `ldconfig -p`, under `/opt/mellanox`, `/opt/hpcx*`, or `/apps`. This
is site breakage, not a script fault, and it means the published trajectory can
no longer be reproduced with the module stack recorded in `run_md.sh`. Worth a
ticket to UF Research Computing independently of this analysis.

Replicates therefore use GROMACS from a dedicated conda environment:

    CONDA_OVERRIDE_CUDA=12.8 conda create -n gmx2025 \
      -c conda-forge --override-channels \
      'gromacs=2025.2=nompi_cuda*' numpy

`CONDA_OVERRIDE_CUDA` is required because login nodes have no NVIDIA driver,
so conda's `__cuda` virtual package is absent and CUDA builds appear
uninstallable. `--override-channels` drops `bioconda` and `defaults`, which
otherwise stall the solve. Installing into `egfr_env` instead was tried and is
not recommended: the solve hangs, and an unpinned attempt silently resolved to
gromacs 2024.5 with cudatoolkit 11.8 while downgrading `icu` in the
environment the alignment and network analysis ran in.

Verified on an `hpg-turin` node (job 50630758):

    GROMACS version:     2025.2-conda_forge
    GPU support:         CUDA
    SIMD instructions:   AVX2_256
    numpy 2.5.4

The `nompi_cuda` build is thread-MPI with CUDA. It links no `libmpi`, so it
cannot hit the HPC-X fault, and the job is single-node/single-GPU
(`-ntmpi 1 -ntomp 8`) so MPI was never needed. The cuda module is deliberately
not loaded — conda ships its own CUDA runtime, and `libcuda.so.1` comes from
the node driver.

**State this in the methods.** The GROMACS *release* matches the published
replicate 1 exactly (2025.2), so only the build provenance differs: conda-forge
rather than the site module, hence different compiler and SIMD flags
(`AVX2_256`). The `.mdp` parameters are identical. One sentence covers it --
e.g. "replicates 2-5 used GROMACS 2025.2 (conda-forge build) after the site
module became unusable" -- and it is far easier to defend than a cross-release
comparison would have been. The script records the resolved binary,
`CONDA_PREFIX` and the full version banner in `provenance_repN.txt` for every
replicate.

The site modules were exhausted first: `gromacs/2026.0` also requires the
broken `openmpi/5.0.7`, and the only builds on a working MPI are
`gromacs/2024.4` (`cuda/12.4.1 gcc/12.2.0 openmpi/4.1.6`) and `gromacs/2023.2`
(`cuda/12.2.2 gcc/12.2.0 openmpi/4.1.5`) -- both a release behind. Either can
be used by setting `GROMACS_SOURCE="module"` and `GROMACS_MODULES`.

The script fails fast with a readable error if `gmx` is missing, if it is
present but will not execute, or if the build reports no CUDA support.

## What this reproduces, and the one thing it changes

GROMACS 2025.2, AMBER99SB-ILDN, TIP3P, PME, Parrinello-Rahman, 310 K, 1 bar;
48,489 atoms; 14,736 SOL, 45 NA, 47 CL; triclinic 8.881 x 8.881 x 6.280 nm.

`md_rep.mdp` is byte-identical to the published `md.mdp`. The published run
inherited its velocities from `npt.cpt` (`gen_vel = no`, `continuation = yes`),
so simply rerunning it reproduces the same trajectory. Each array task
therefore inserts a **1 ns re-equilibration (`eq_seed.mdp`) with
`gen_vel = yes` and a distinct `gen_seed`** before production. That stage is
the only difference between replicates. Seeds are `20250000 + REP`, recorded
with the GROMACS version and input SHA256s in `provenance_repN.txt`.

`eq_seed.mdp` uses C-rescale rather than Parrinello-Rahman, because PR is not
robust immediately after velocity generation; production restores PR. Because
production starts from an equilibrated per-seed state, the full 100 ns is
analysed with no leading frames discarded, keeping each replicate directly
comparable to the published one.

## Inputs read from the project directory

    npt.gro  npt.cpt  topol.top  posre.itp  index.ndx

The force field is not staged — `topol.top` includes `amber99sb-ildn.ff/...`,
which GROMACS resolves from its own share directory, and no local `.ff`
directory exists in the working folder.

`gmx trjconv`, `gmx covar` and `gmx rmsf` select groups by **name** on stdin
(`Protein`, `System`, `C-alpha`), following `ca_extract.sh`, so nothing depends
on the group ordering inside `index.ndx`. `index.ndx` is still passed to
`grompp`, as in `run_md.sh`.

## Running

    sbatch md_replicates.sbatch

Array IDs 2–5 produce replicates 2–5; replicate 1 is the published run.

## Outputs, per replicate

    egfr_dccm_matrix_repN.csv   253 x 253, unit diagonal, residues 714-966
    rmsf_repN.xvg               per-residue, C-alpha
    provenance_repN.txt         seed, GROMACS version, node, input hashes

`dccm_from_covar.py` replaces `calculate_dccm_fast.py`, which the committed
`pipeline/04_run_md.sh` calls but which was never added to the repository. It
computes `DCCM_ij = Tr(Cov_ij) / sqrt(Tr(Cov_ii) Tr(Cov_jj))` from the
3N x 3N `gmx covar -ascii` output and subsets the 259-residue construct
(712-970) to the 253 scored residues (714-966; 0-based C-alpha indices 2:255).
Verified against a direct block-trace reference to 5e-7.

## After the array finishes

    python3 aggregate_replicates.py --published egfr_dccm_matrix.csv

Writes `replicate_dccm_correlation.csv` (pairwise Pearson r over the strict
upper triangle) and `replicate_dccm_mean.csv`. This is the between-replicate
DCCM correlation §9.1 asks for.

Still to be run separately, per §9.1 and §6: the distance-matched enrichment
null applied to each replicate, so per-replicate enrichment can be reported
alongside the published 1.02-fold value.

## Provenance discrepancy worth fixing in the repo

The committed `pipeline/04_run_md.sh` is a reconstruction and does not match
the script that was actually submitted. Against the real `run_md.sh` it:

- omits `-n index.ndx` from `grompp`;
- writes `md.tpr` but runs `mdrun -deffnm md_production`, which reads
  `md_production.tpr` — the real script writes `md_production.tpr` directly;
- omits the module versions (`cuda/12.8.1`, `gcc/14.2.0`, `openmpi/5.0.7`),
  giving only `module load gromacs`;
- omits `-update gpu` and `-cpt 15`;
- folds in `trjconv`/`covar`/`rmsf` steps that in practice lived in
  `ca_extract.sh`, and shows them without the group selections they require.

Anyone reproducing from the repository alone would not be able to run it. The
real scheduler script, module versions and group selections belong in the
committed pipeline.
