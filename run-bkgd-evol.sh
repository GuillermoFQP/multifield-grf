#!/usr/bin/env bash
set -e
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"

indices=({1..1024})       # Or, for example: indices=(3 17 42)
workers=8

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1
#mkdir -p gradhess-new

# Give each background job its own process group.
set -m
stop_jobs() {
    local pid
    for pid in $(jobs -pr); do
        kill -KILL -- "-$pid" 2>/dev/null || true
    done
    wait 2>/dev/null || true
}
trap stop_jobs EXIT
trap 'exit 130' INT
trap 'exit 143' TERM

for i in "${indices[@]}"; do
    while :; do
        active=($(jobs -pr))
        (( ${#active[@]} < workers )) && break
        sleep 0.1
    done
    ./mass_projs "data/smpout-new/smpout-${i}.fits" > "plots/figure-5/gradhess-new/gradhess-${i}.txt" &
done
wait
