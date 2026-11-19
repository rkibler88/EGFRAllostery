#!/bin/bash
# Post-mortem for the EGFR replicate MD array (CORRECTIONS.md 9.1) on HiPerGator.
# Companion to monitor_replicates.sh: that one watches a healthy run, this one
# explains a dead one. Everything happens over ssh; nothing is written remotely.
#
#   ./diagnose_replicates.sh                  newest array named EGFR_100ns_rep
#   ./diagnose_replicates.sh 50712345         a specific array job id
#   HOST=user@host ./diagnose_replicates.sh
#   ./diagnose_replicates.sh > crash_report.txt 2>&1
#
# Reads, in order: slurm accounting -> stdout/stderr -> GROMACS fatal errors ->
# per-replicate file inventory -> staging inputs -> gmx environment -> quota.
set -uo pipefail
HOST="${HOST:-rkibler@hpg.rc.ufl.edu}"
JOBID="${1:-}"

# `bash -ls` is a LOGIN shell: without it /opt/slurm/bin is off PATH and sacct
# is not found. `-s --` passes the job id as $1 to the script read on stdin.
# The quoted heredoc delimiter stops the LOCAL shell touching $P, $r or globs.
ssh "${HOST}" "bash -ls -s -- '${JOBID}'" <<'REMOTE'
JOBID="${1:-}"
P=/blue/raquel.dias/rkibler/EGFR
W="$P/replicates"
REPS="2 3 4 5"
hr() { printf '\n===== %s =====\n' "$1"; }

# --- 0. identify the array job ---------------------------------------------
if [ -z "$JOBID" ]; then
    JOBID=$(sacct -u "$USER" --name=EGFR_100ns_rep -S now-30days -X -n -P -o JobID 2>/dev/null \
            | sed 's/_.*$//; s/\..*$//' | sort -un | tail -1)
fi
if [ -z "$JOBID" ]; then
    echo "No job named EGFR_100ns_rep in the last 30 days for $USER."
    echo "If the name differs, pick the id from this list and pass it as an argument:"
    sacct -u "$USER" -S now-30days -X -o JobID%14,JobName%24,State%18,End%20 2>/dev/null | tail -25
else
    echo "array job: $JOBID   (user $USER, $(date -u +%Y-%m-%dT%H:%M:%SZ))"
fi

if [ -n "$JOBID" ]; then
# --- 1. what slurm recorded -------------------------------------------------
hr "slurm accounting"
sacct -j "$JOBID" --units=M \
  -o JobID%18,State%20,ExitCode%8,Elapsed%10,MaxRSS%9,ReqMem%8,NodeList%14,Start%20,End%20 2>/dev/null

hr "verdict per task"
# Read the allocation rows only (-X); the .batch/.extern steps repeat the state.
sacct -j "$JOBID" -X -n -P -o JobID,State,ExitCode,Elapsed 2>/dev/null |
while IFS='|' read -r id st ec el; do
    case "$st" in
      TIMEOUT*)      v="hit the 48 h wall. Resubmit that index: -cpi resumes from md_production.cpt";;
      OUT_OF_ME*)    v="killed by cgroup OOM: raise --mem above 30gb";;
      NODE_FAIL*)    v="the node died, not your script. Resubmit";;
      CANCELLED*)    v="cancelled: scancel, a QOS/account limit, or preemption";;
      PENDING*)      v="still queued, never ran";;
      RUNNING*)      v="still running";;
      COMPLETED*)    v="completed";;
      FAILED*)
        case "${ec%%:*}" in
          127) v="exit 127: gmx not on PATH, or it is but a shared library is unresolved";;
          1)   v="exit 1: a guard in the sbatch tripped (see stdout), or GROMACS called a fatal error";;
          2)   v="exit 2: shell/syntax or a missing file under 'set -e'";;
          *)   v="exit ${ec}";;
        esac;;
      *) v="$st";;
    esac
    # A nonzero signal field (after the colon) means it was killed, not exited.
    case "$ec" in *:0) ;; *:*) v="$v  [killed by signal ${ec##*:}]";; esac
    printf '  %-16s %-14s exit=%-8s %-10s %s\n' "$id" "$st" "$ec" "$el" "$v"
done

# --- 2. stdout / stderr -----------------------------------------------------
hr "log files"
LOGDIR=""
for cand in "$P/logs" "$HOME/logs" "$W/logs"; do
    [ -n "$(ls "$cand"/md_rep*_*.out 2>/dev/null)" ] && { LOGDIR="$cand"; break; }
done
if [ -z "$LOGDIR" ]; then
    echo "No md_rep*.out anywhere obvious. Submit directory according to slurm:"
    sacct -j "$JOBID" -n -P -o WorkDir 2>/dev/null | sort -u | sed 's/^/  /'
    find "$P" "$HOME" -maxdepth 4 -name 'md_rep*.out' -mtime -30 2>/dev/null | sed 's/^/  found: /'
    echo
    echo "  If there is NO output file at all, that is itself the answer:"
    echo "  --output=logs/md_rep%a_%A.out is relative to the submit directory, and"
    echo "  slurm opens it BEFORE the script's mkdir -p runs. Submit from a"
    echo "  directory that already contains logs/, or the task dies writing nothing."
else
    echo "log directory: $LOGDIR"
fi

for r in $REPS; do
    o=$(ls -t "$LOGDIR"/md_rep${r}_${JOBID}.out 2>/dev/null | head -1)
    [ -z "$o" ] && o=$(ls -t "$LOGDIR"/md_rep${r}_*.out 2>/dev/null | head -1)
    hr "rep$r stdout/stderr"
    if [ -z "$o" ]; then echo "  (no stdout file)"; continue; fi
    echo "--- $o (last 40 lines)"
    tr '\r' '\n' < "$o" | grep -v '^[[:space:]]*$' | tail -40 | sed 's/^/  /'
    e="${o%.out}.err"
    if [ -s "$e" ]; then
        echo "--- $e (last 40 lines)"
        tr '\r' '\n' < "$e" | grep -v '^[[:space:]]*$' | tail -40 | sed 's/^/  /'
    else
        echo "--- ${e##*/}: empty or absent"
    fi
done

# --- 3. targeted scan for known failure modes -------------------------------
hr "fatal-error scan (logs + GROMACS .log files)"
PAT='error while loading shared libraries|command not found|No such file or directory|Disk quota exceeded|No space left on device|Cannot allocate memory|oom-kill|Out Of Memory|slurmstepd:|CondaError|EnvironmentLocationNotFound|EnvironmentNameNotFound|CUDA error|cudaError|cudaErrorNoDevice|No compatible GPUs|Segmentation fault|core dumped|DUE TO TIME LIMIT|CANCELLED AT'
for f in "$LOGDIR"/md_rep*_*.out "$LOGDIR"/md_rep*_*.err "$W"/rep*/*.log; do
    [ -e "$f" ] || continue
    hits=$(tr '\r' '\n' < "$f" | grep -n -E "$PAT" | sort -u -t: -k2 | head -6)
    [ -n "$hits" ] && { echo "--- ${f#$P/}"; echo "$hits" | sed 's/^/  /'; }
done
# GROMACS prints the real message on the lines AFTER "Fatal error:".
for f in "$LOGDIR"/md_rep*_*.out "$LOGDIR"/md_rep*_*.err "$W"/rep*/*.log; do
    [ -e "$f" ] || continue
    blk=$(tr '\r' '\n' < "$f" | grep -A6 -E '^Fatal error' | head -20)
    [ -n "$blk" ] && { echo "--- GROMACS fatal error in ${f#$P/}"; echo "$blk" | sed 's/^/  /'; }
done
echo "(nothing printed above this line means no known pattern matched)"

# --- 4. how far each replicate got ------------------------------------------
hr "per-replicate file inventory"
for r in $REPS; do
    d="$W/rep$r"
    printf 'rep%-2s ' "$r"
    if [ ! -d "$d" ]; then
        echo "no run directory -- died before staging (module/conda/gmx probe, or mkdir/cp)"
        continue
    fi
    got=""
    for f in npt.gro npt.cpt topol.top index.ndx eq.tpr eq.gro md_production.tpr \
             md_production.log md_full.xtc md_clean.xtc covar.ascii; do
        [ -e "$d/$f" ] && got="$got $f"
    done
    [ -e "$W/egfr_dccm_matrix_rep$r.csv" ] && got="$got DCCM.csv"
    echo "have:${got:- (empty directory)}"
    log=$(ls -t "$d"/md_production*.log "$d"/eq*.log 2>/dev/null | head -1)
    if [ -n "$log" ]; then
        st=$(grep -A1 -E '^ *Step +Time' "$log" | tail -1 | awk 'NF==2 {printf "step %s (%.1f ps)", $1, $2}')
        printf '      last log %s  %s  mtime %s\n' "${log##*/}" "${st:-no step line}" \
               "$(date -r "$log" -u +%Y-%m-%dT%H:%MZ 2>/dev/null)"
        tail -5 "$log" | sed 's/^/      | /'
    fi
done
fi  # end JOBID block

# --- 5. inputs the sbatch stages (a missing one kills it instantly) ---------
hr "staging inputs in $P"
for f in npt.gro npt.cpt topol.top posre.itp index.ndx \
         eq_seed.mdp md_rep.mdp dccm_from_covar.py md_replicates.sbatch; do
    if [ -e "$P/$f" ]; then ls -lh "$P/$f" | awk '{printf "  ok      %-26s %6s  %s %s %s\n",$9,$5,$6,$7,$8}'
    else echo "  MISSING $f"; fi
done
echo "  (the .mdp/.py/.sbatch three are read from the SUBMIT directory, not \$P --"
echo "   if you sbatch'd from elsewhere, that is where they must be)"

# --- 6. can gmx even start here? -------------------------------------------
hr "gromacs environment (login node; GPU checks can only run on a compute node)"
set +u
module load conda >/dev/null 2>&1
if conda env list 2>/dev/null | awk '{print $1}' | grep -qx gmx2025; then
    PRE=$(conda env list 2>/dev/null | awk '$1=="gmx2025"{print $NF}')
    echo "  conda env gmx2025: $PRE"
    if [ -x "$PRE/bin/gmx" ]; then
        echo "  binary: $PRE/bin/gmx"
        miss=$(ldd "$PRE/bin/gmx" 2>/dev/null | grep -i 'not found')
        [ -n "$miss" ] && { echo "  UNRESOLVED LIBRARIES (this is exit 127):"; echo "$miss" | sed 's/^/    /'; }
        "$PRE/bin/gmx" --version 2>&1 | grep -i -E 'GROMACS version|GPU support|SIMD' | sed 's/^/    /' \
            || echo "    gmx --version failed to run at all"
    else
        echo "  NO gmx in $PRE/bin -- the env exists but GROMACS is not installed in it"
    fi
else
    echo "  conda env gmx2025 NOT FOUND (GROMACS_SOURCE=conda in the sbatch expects it)"
    conda env list 2>/dev/null | sed 's/^/    /'
fi
set -u

# --- 7. storage ------------------------------------------------------------
hr "quota and space"
command -v blue_quota >/dev/null && blue_quota 2>/dev/null | sed 's/^/  /'
command -v home_quota >/dev/null && home_quota 2>/dev/null | sed 's/^/  /'
df -h "$P" 2>/dev/null | sed 's/^/  /'
echo
echo "A full /blue allocation fails every task at the same point and looks like a"
echo "code bug; check the numbers above before changing anything in the sbatch."
REMOTE
