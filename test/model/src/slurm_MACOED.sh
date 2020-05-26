#!/bin/bash
#SBATCH --job-name=MACOED_ME
#SBATCH --output=../bin/MACOED_ME.txt
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mail-type=ALL
#SBATCH --mem-per-cpu=20000
#SBATCH --mail-user=david.blumenthal@wzw.tum.de

export OMP_NUM_THREADS=10

mkdir -p ../res/MACOED/ME/
#srun --ntasks 1 --cpus-per-task $OMP_NUM_THREADS ../bin/compare_models --input-directory ../../../data/MACOED/ME/ --input-format CSV_SNPS_AS_COLUMNS_FIRST --pheno-type CATEGORICAL --num-threads 10 --output-directory ../res/MACOED/ME/ --disease-snps 0 1
../bin/compare_models --input-directory ../../../data/MACOED/ME/ --input-format CSV_SNPS_AS_COLUMNS_FIRST --pheno-type CATEGORICAL --num-threads 10 --output-directory ../res/MACOED/ME/ --disease-snps 0 1
mkdir -p ../res/MACOED/NME/
#srun --ntasks 1 --cpus-per-task $OMP_NUM_THREADS ../bin/compare_models --input-directory ../../../data/MACOED/NME/ --input-format CSV_SNPS_AS_COLUMNS_FIRST --pheno-type CATEGORICAL --num-threads 10 --output-directory ../res/MACOED/NME/ --disease-snps 0 1
../bin/compare_models --input-directory ../../../data/MACOED/NME/ --input-format CSV_SNPS_AS_COLUMNS_FIRST --pheno-type CATEGORICAL --num-threads 10 --output-directory ../res/MACOED/NME/ --disease-snps 0 1