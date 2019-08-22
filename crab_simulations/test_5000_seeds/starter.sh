#!/bin/bash
#SBATCH --partition=deeplearning
#SBATCH --job-name=test
SEED=$1 python exec.py
