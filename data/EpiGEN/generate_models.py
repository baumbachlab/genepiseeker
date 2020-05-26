import os.path
import argparse
import random

parser = argparse.ArgumentParser(description="Generate categorical EpiGEN models.")
parser.add_argument("--num-disease-snps", help="Number of disease SNPs.", required=True, type=int)
parser.add_argument("--num-models", help="Number of models.", default=2, type=int)
parser.add_argument("--num-categories", help="Number of categories.", default=3, type=int)
args = parser.parse_args()

models = ["model_{}".format(model_id) for model_id in range(1, args.num_models + 1)]
filenames = [os.path.join(os.path.dirname(os.path.abspath(__file__)), "categorical", "{}_disease_snps".format(args.num_disease_snps), model, "model.ini") for model in models]
genotypes = [""]
for snp in range(args.num_disease_snps):
    genotypes = ["{}{},".format(genotype, genotype_at_snp) for genotype in genotypes for genotype_at_snp in range(3)]
genotypes = [genotype[:-1] for genotype in genotypes]

for filename in filenames:
    with open(filename, "w") as f:
        f.write("[Model Type]\n")
        f.write("size = {}\n".format(args.num_disease_snps))
        f.write("phenotype = {}\n".format(args.num_categories))
        f.write("\n")
        f.write("[Model Definition]\n")
        for genotype in genotypes:
            distr =[random.random() for c in range(args.num_categories)]
            sum_distr = sum(distr)
            distr = [prob / sum_distr for prob in distr]
            line = genotype + " = " + str(distr[0])
            for pos in range(1,len(distr)):
                line += "," + str(distr[pos])
            f.write("{}\n".format(line))
