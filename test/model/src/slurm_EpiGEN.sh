##!/bin/bash
#SBATCH --job-name=EpiGEN
#SBATCH --output=../bin/EpiGEN.txt
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mail-type=ALL
#SBATCH --mem-per-cpu=20000
#SBATCH --mail-user=david.blumenthal@wzw.tum.de

export OMP_NUM_THREADS=10

#for PHENOTYPE in "quantitative" "dichotomous"
#do
#	for SIZE in 2 3 4
#	do
#		for MODEL in "exponential" "multiplicative" "joint-dominant" "joint-recessive"
#		do
#			mkdir -p ../res/EpiGEN/${PHENOTYPE}/${SIZE}_disease_snps/${MODEL}
			srun --ntasks 1 --cpus-per-task $OMP_NUM_THREADS ../bin/compare_models --input-directory ../../../data/EpiGEN/${PHENOTYPE}/${SIZE}_disease_snps/${MODEL}/ --input-format JSON_EPIGEN --pheno-type ${PHENOTYPE^^} --num-threads 10 --output-directory ../res/EpiGEN/${PHENOTYPE}/${SIZE}_disease_snps/${MODEL}/
#			../bin/compare_models --input-directory ../../../data/EpiGEN/${PHENOTYPE}/${SIZE}_disease_snps/${MODEL}/ --input-format JSON_EPIGEN --pheno-type ${PHENOTYPE^^} --num-threads 10 --output-directory ../res/EpiGEN/${PHENOTYPE}/${SIZE}_disease_snps/${MODEL}/
#		done
#	done
#done

for PHENOTYPE in "categorical"
do
	for SIZE in 2 3 4
	do
		for MODEL in "model_1" "model_2"
		do
			mkdir -p ../res/EpiGEN/${PHENOTYPE}/${SIZE}_disease_snps/${MODEL}
#			srun --ntasks 1 --cpus-per-task $OMP_NUM_THREADS ../bin/compare_models --input-directory ../../../data/EpiGEN/${PHENOTYPE}/${SIZE}_disease_snps/${MODEL}/ --input-format JSON_EPIGEN --pheno-type ${PHENOTYPE^^} --num-threads 10 --output-directory ../res/EpiGEN/${PHENOTYPE}/${SIZE}_disease_snps/${MODEL}/
			../bin/compare_models --input-directory ../../../data/EpiGEN/${PHENOTYPE}/${SIZE}_disease_snps/${MODEL}/ --input-format JSON_EPIGEN --pheno-type ${PHENOTYPE^^} --num-categories 3 --num-threads 10 --output-directory ../res/EpiGEN/${PHENOTYPE}/${SIZE}_disease_snps/${MODEL}/
		done
	done
done