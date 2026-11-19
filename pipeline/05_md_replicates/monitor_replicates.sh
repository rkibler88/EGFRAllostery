#!/bin/bash
# Monitor the EGFR replicate MD array (CORRECTIONS.md 9.1) on HiPerGator.
# Run from a local machine; all work happens over ssh on the cluster.
#
#   ./monitor_replicates.sh        one snapshot
#   ./monitor_replicates.sh 120    refresh every 120 s until Ctrl-C
#
# Override the login target with:  HOST=user@host ./monitor_replicates.sh
set -uo pipefail
HOST="${HOST:-rkibler@hpg.rc.ufl.edu}"
INTERVAL="${1:-0}"

snapshot() {
# `bash -ls` is a LOGIN shell: without it /opt/slurm/bin is off PATH and squeue
# is not found. The quoted heredoc delimiter stops the LOCAL shell expanding
# $P, $r and the globs -- zsh aborts outright on an unmatched glob.
ssh "${HOST}" 'bash -ls' <<'REMOTE'
P=/blue/raquel.dias/rkibler/EGFR
echo "===== queue ====="
squeue -u rkibler -o "%.12i %.8T %.10M %.4C %b %R" 2>/dev/null || echo "squeue unavailable"

echo "===== progress ====="
# Progress is read from the GROMACS .log, not from mdrun -v on stderr: the
# step counter is only reliably emitted to a terminal, so the .err file can
# stay empty for a long time while the run is perfectly healthy.
for r in 2 3 4 5; do
    d="$P/replicates/rep$r"
    printf "rep%-3s " "$r"
    if [ -f "$P/replicates/egfr_dccm_matrix_rep$r.csv" ]; then
        echo "DONE - DCCM written"
        continue
    fi
    log=$(ls -t "$d"/md_production*.log 2>/dev/null | head -1); stage="production  "
    if [ -z "$log" ]; then
        log=$(ls -t "$d"/eq*.log 2>/dev/null | head -1); stage="equilibration"
    fi
    if [ -z "$log" ]; then
        if [ -d "$d" ]; then echo "staged, no GROMACS log yet"; else echo "not started"; fi
        continue
    fi
    st=$(grep -A1 -E '^ *Step +Time' "$log" | tail -1 \
         | awk 'NF==2 {printf "step %s, %.1f ps", $1, $2}')
    perf=$(grep -E 'Performance:' "$log" 2>/dev/null | tail -1 | awk '{print $2" ns/day"}')
    echo "$stage  ${st:-starting}${perf:+   [$perf]}"
done

echo "===== startup banner (rep2) ====="
grep -h -E 'conda env|GROMACS binary|GROMACS version|GPU support' \
     "$P"/logs/md_rep2_*.out 2>/dev/null | head -6 || echo "no stdout yet"

echo "===== any stderr ====="
for r in 2 3 4 5; do
    f=$(ls -t "$P"/logs/md_rep"$r"_*.err 2>/dev/null | head -1)
    [ -n "$f" ] && [ -s "$f" ] && { printf "rep%-3s " "$r"; tail -c 300 "$f" | tr '\r' '\n' | grep -v '^$' | tail -1; }
done
true
REMOTE
}

case "${INTERVAL}" in
    ''|*[!0-9]*|0) snapshot ;;
    *)             while true; do clear; date; snapshot; sleep "${INTERVAL}"; done ;;
esac
