#!/bin/bash
#SBATCH --job-name=BOOST_ME
#SBATCH --output=../bin/BOOST_ME.txt
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mail-type=ALL
#SBATCH --mem-per-cpu=20000
#SBATCH --mail-user=david.blumenthal@wzw.tum.de

export OMP_NUM_THREADS=10

for SIZE in 400 800 1200 1600
do
	for ID in {1..12}
	do
		mkdir -p ../res/BOOST/ME/$SIZE/$ID
#		srun --ntasks 1 --cpus-per-task $OMP_NUM_THREADS ../bin/compare_models --input-directory ../../../data/BOOST/ME/$SIZE/$ID/ --input-format CSV_SNPS_AS_COLUMNS_FIRST --pheno-type CATEGORICAL --num-threads 10 --output-directory ../res/BOOST/ME/$SIZE/$ID/ --disease-snps 0 999
		../bin/compare_models --input-directory ../../../data/BOOST/ME/$SIZE/$ID/ --input-format CSV_SNPS_AS_COLUMNS_FIRST --pheno-type CATEGORICAL --num-threads 10 --output-directory ../res/BOOST/ME/$SIZE/$ID/ --disease-snps 0 999
	done
done