#!/bin/bash
#SBATCH -M es
#SBATCH -A bil-fire3
#SBATCH -J cmaq2hemco
#SBATCH -o logs/%x-%j.out
#SBATCH -t 00:20:00
#SBATCH -p dtn_f5_f6
#SBATCH -N 1

d=$1
pushd /gpfs/f6/bil-fire3/world-shared/Barry.Baker/cmaq2hemco

/gpfs/f6/bil-fire3/world-shared/python/envs/cmaq2hemco/bin/python examples/aws_mp2022.py --date $d
