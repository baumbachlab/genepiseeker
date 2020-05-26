import os.path
import subprocess
import argparse

parser = argparse.ArgumentParser(description="Generate EpiGEN data.")
parser.add_argument("--epigen-root", help="Path to EpiGEN's root directory.", required=True)
args = parser.parse_args()

model_types = ["dichotomous", "quantitative"]
nums_disease_snps = ["2_disease_snps", "3_disease_snps", "4_disease_snps"]
interaction_types = ["exponential", "joint-dominant", "joint-recessive", "multiplicative"]

subdirs = [os.path.join(os.path.dirname(os.path.abspath(__file__)), m_type, snps, i_type) for m_type in model_types for snps in nums_disease_snps for i_type in interaction_types]

result_files = os.path.join(args.epigen_root, "sim", "*_1_ASW.json")
#for subdir in subdirs:
#    epigen_options = "--corpus-id 1 --pop ASW --snps 100 --inds 1600 --num-sims 100 --disease-maf-range 0.2 0.4 --model {}".format(os.path.join(subdir, "model.xml"))
#    subprocess.call("cd {}; python simulate_data.py {}".format(args.epigen_root, epigen_options), shell=True)
#    subprocess.call("mv {} {}".format(result_files, subdir), shell=True)

models = ["model_1", "model_2"]
categorical_subdirs = [os.path.join(os.path.dirname(os.path.abspath(__file__)), "categorical", snps, model) for snps in nums_disease_snps for model in models]
for subdir in categorical_subdirs:
    epigen_options = "--corpus-id 1 --pop ASW --snps 100 --inds 1600 --num-sims 100 --disease-maf-range 0.2 0.4 --model {}".format(os.path.join(subdir, "model.ini"))
    subprocess.call("cd {}; python simulate_data.py {}".format(args.epigen_root, epigen_options), shell=True)
    subprocess.call("mv {} {}".format(result_files, subdir), shell=True)