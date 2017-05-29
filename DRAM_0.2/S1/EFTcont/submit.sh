#!/bin/bash
#$ -q cli.q
#$ -pe openmp 16
#$ -cwd
#$ -S /bin/bash
#$ -o output.log
#$ -e error.log
#$ -l h_rt=06:00:00

./DRAM
