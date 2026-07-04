#!/bin/bash

N=${1:?Usage: $0 N}

for i in $(seq 15 $N); do
    echo ""
    echo "Gaussian random field potential realization $i/$N"
    echo "Computing trajectory..."
    if ! gtimeout 60 bash -c "./bkgd_onetraj smpout_${i}.fits > traj_${i}.txt"; then
        echo "Realization $i exceeded 60s or failed, skipping."
        continue
    fi
    echo "Computing power spectra..."
    ./ps_full smpout_${i}.fits > ps_full_grf_${i}.txt
    echo ""
done

echo "Done. Computed $N power spectra."
