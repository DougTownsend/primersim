#!/bin/bash
#BSUB -q shared_memory
#BSUB -n 16
#BSUB -R "span[hosts=1]"
#BSUB -W 360
#BSUB -J "prime_sweep[1-8]"
#BSUB -o /share/tuck/dktownse/prime-anneal/lsf.%J.%I.out
#BSUB -e /share/tuck/dktownse/prime-anneal/lsf.%J.%I.err

# Single ./test_eq invocation per array element. n_trials=0 means
# "loop forever inside the binary"; LSF SIGTERMs the process at -W
# (60 min). Trials in flight at SIGTERM (≤ num_cpu = 32) just don't
# produce a CSV; completed trials are intact on disk.
# Each array element writes to its own sweep_job_NN/ subdir.

set -euo pipefail

cd /share/tuck/dktownse/prime-anneal

JOB=$(printf "%02d" "${LSB_JOBINDEX:-0}")
out_dir="sweep_job_${JOB}"
mkdir -p "${out_dir}"

echo "host: $(hostname)  job array index: ${LSB_JOBINDEX}  out: ${out_dir}"

./test_eq "${LSB_DJOB_NUMPROC:-32}" sweep 0 "$out_dir"
