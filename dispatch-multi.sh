#!/bin/bash
#
# SGE (Sun Grid Engine) submission script for an array of runs of the MC code.
#
# To submit this job to be run on the cluster, type:
#
# qsub -t 1-<array_size> ./dispatch-multi.sh [number_of_trials] [output_file_prefix]
#
# at the command line. This will queue an array of size <array_size> jobs, each of
# which will comprises number_of_trials trials (default=1,000,000). Output will be
# written to a file in the data/ directory that begins with output_file_prefix and
# ends with job<n>.dat, where <n> = 1, 2, ..., array_size. The default prefix is pi-1e6
#
#####################################################################
#
# Lines beginning with '#$' below provide arguments the qsub command,
# and save you typing them out each time.
#
# (1) Run the job in the same directory that qsub is issued
#
#$ -cwd
#
# (2) Give the job a meaningful name (this will appear in qstat listings)
#
#$ -N hmc
#
# (3) Send a notification email when your job begins (b) and ends (e)
# being run
#
#$ -m be
#
# (4) Set the queue that the job should run on. As undergraduates, you have
# access to the sopa.1.day queue.
#
#$ -q sopa.1.day
#
# (5) Give an estimate of how long the job should run for, expressed as
# hours:minutes:seconds. This can be up to 24 hours on the sopa.1.day queue.
# NB: your job will be halted when this time expires, so it is best to slightly
# over-estimate.
#
#$ -l h_rt=10:00:00

mkdir -p data # This makes sure the data directory can be written to
./mcpi.py -t ${1:-1000000} -o data/${2:-pi-1e6}-job${SGE_TASK_ID}.dat

# The ${n:-default} syntax looks for the n'th argument supplied to the
# script, and substitutes default if this is empty.
