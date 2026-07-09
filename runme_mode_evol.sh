#!/bin/bash
N=${1:?Usage: $0 N}
mkdir -p modes

run_one() {
    local i="$1"
    local N="$2"
    echo "[$i/$N] Computing mode evolution..."
    if ! gtimeout 100 ./ps_mode "smpout_tour_0/smpout_${i}.fits" > "modes/ps_full_grf_${i}.txt"; then
        echo "[$i/$N] Exceeded 60s or failed, skipping."
        rm -f "modes/ps_full_grf_${i}.txt"
    else
        echo "[$i/$N] Done."
    fi
}
export -f run_one

seq 1 "$N" | xargs -P 8 -I{} bash -c 'run_one "$0" "$1"' {} "$N"

echo "Done. Computed $N mode evolutions."
