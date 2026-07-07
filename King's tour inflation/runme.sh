#!/bin/bash

N=${1:?Usage: $0 N}

seq 1 $N | xargs -P 8 -I {} bash -c 'echo "Computing realization {} of '"$N"'..."; python make-grf.py smpout_{}.fits {}'

echo "Done. Generated $N realizations."
