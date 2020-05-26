/*******************************************************************************
 *                                                                             *
 *   Copyright (C) 2020 by David B. Blumenthal                                 *
 *                                                                             *
 *   This file is part of GenEpiSeeker.                                        *
 *                                                                             *
 *   GenEpiSeeker is free software: you can redistribute it and/or modify it   *
 *   under the terms of the GNU General Public License as published by         *
 *   the Free Software Foundation, either version 3 of the License, or         *
 *   (at your option) any later version.                                       *
 *                                                                             *
 *   GenEpiSeeker is distributed in the hope that it will be useful,           *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the              *
 *   GNU General Public License for more details.                              *
 *                                                                             *
 *   You should have received a copy of the GNU General Public License         *
 *   along with GenEpiSeeker. If not, see <http://www.gnu.org/licenses/>.      *
 *                                                                             *
 ******************************************************************************/


/*!
 * @file  bayesian_model.ipp
 * @brief Definition of epi::BayesianModel.
 */

#ifndef SRC_MODEL_BAYESIAN_MODEL_IPP_
#define SRC_MODEL_BAYESIAN_MODEL_IPP_

namespace epi {

template<class PhenoType>
BayesianModel<PhenoType>::
BayesianModel(Instance<PhenoType> * instance):
EpistasisModel<PhenoType>(instance),
mu_{0.0},
sigma_{0.0},
num_bins_{100} {}

template<class PhenoType>
BayesianModel<PhenoType>::
~BayesianModel() {}

template<>
void
BayesianModel<QuantitativePhenoType>::
initialize_() {
	mu_ = 0.0;
	for (Ind ind{0}; ind < this->instance_->num_inds(); ind++) {
		mu_ += this->instance_->phenotype(ind);
	}
	mu_ /= static_cast<double>(this->instance_->num_inds());
	sigma_ = 0.0;
	for (Ind ind{0}; ind < this->instance_->num_inds(); ind++) {
		sigma_ += (this->instance_->phenotype(ind) - mu_) * (this->instance_->phenotype(ind) - mu_);
	}
	sigma_ /= static_cast<double>(this->instance_->num_inds());
}

template<>
void
BayesianModel<CategoricalPhenoType>::
initialize_() {
	num_bins_= instance_->num_categories();
}

template<>
void
BayesianModel<QuantitativePhenoType>::
initialize_binned_phenotypes_(std::vector<CategoricalPhenoType> & binned_phenotypes) const {
	CategoricalPhenoType bin;
	for (Ind ind{0}; ind < this->instance_->num_inds(); ind++) {
		bin = static_cast<CategoricalPhenoType>(std::floor(misc::normal_cdf(this->instance_->phenotype(ind), mu_, sigma_) * static_cast<double>(num_bins_)));
		if (bin == num_bins_) {
			bin--;
		}
		binned_phenotypes.emplace_back(bin);
	}
}

template<>
void
BayesianModel<CategoricalPhenoType>::
initialize_binned_phenotypes_(std::vector<CategoricalPhenoType> & binned_phenotypes) const {
	for (Ind ind{0}; ind < this->instance_->num_inds(); ind++) {
		binned_phenotypes.emplace_back(this->instance_->phenotype(ind));
	}
}

template<>
bool
BayesianModel<CategoricalPhenoType>::
parse_option_(const std::string & option, const std::string & arg) {
	if (option == "num-bins") {
		return true;
	}
	return false;
}

template<>
bool
BayesianModel<QuantitativePhenoType>::
parse_option_(const std::string & option, const std::string & arg) {
	if (option == "num-bins") {
		try {
			num_bins_ = std::stoul(arg);
		}
		catch (...) {
			throw Error(std::string("Invalid argument \"") + arg  + "\" for option \"--" + option + "\". Usage: options = \"[--" + option + " <convertible to int greater equal 2>] [...]\"");
		}
		if (num_bins_ < 2) {
			throw Error(std::string("Invalid argument \"") + arg  + "\" for option \"--" + option + "\". Usage: options = \"[--" + option + " <convertible to int greater equal 2>] [...]\"");
		}
		return true;
	}
	return false;
}

template<>
void
BayesianModel<QuantitativePhenoType>::
set_default_options_() {
	num_bins_ = 100;
}

template<>
void
BayesianModel<CategoricalPhenoType>::
set_default_options_() {
	num_bins_ = this->instance_->num_categories();
}

template<class PhenoType>
std::string
BayesianModel<PhenoType>::
valid_options_() const {
	return "[--num-bins <convertible to int greater equal 2>]";
}

template<class PhenoType>
options::ModelSense
BayesianModel<PhenoType>::
model_sense_() const {
	return options::ModelSense::MINIMIZE;
}

template<class PhenoType>
bool
BayesianModel<PhenoType>::
is_predictive_() const {
	return false;
}

template<class PhenoType>
double
BayesianModel<PhenoType>::
evaluate_(const std::vector<SNP> & snp_set, bool prepare_prediction) {

	// Initialize binned phenotypes.
	std::vector<CategoricalPhenoType> binned_phenotypes;
	initialize_binned_phenotypes_(binned_phenotypes);

	// Construct the penetrance table.
	std::vector<std::vector<PhenoType>> penetrance_table(misc::size_penetrance_table(snp_set.size()), std::vector<PhenoType>());
	for (Ind ind{0}; ind < this->instance_->num_inds(); ind++) {
		penetrance_table.at(this->instance_->genotype_at_snp_set(snp_set, ind)).emplace_back(binned_phenotypes.at(ind));
	}

	// Compute log-K2 score.
	double score{0.0};
	std::vector<std::size_t> num_cell_inds_in_categories(num_bins_, 0);
	for (const auto & cell : penetrance_table) {
		for (auto & count : num_cell_inds_in_categories) {
			count = 0;
		}
		for (const auto & phenotype : cell) {
			num_cell_inds_in_categories.at(phenotype)++;
		}
		score += log_factorial_(cell.size() + 1);
		for (CategoricalPhenoType phenotype{0}; phenotype < num_bins_; phenotype++) {
			score -= log_factorial_(num_cell_inds_in_categories.at(phenotype));
		}
	}

	// Return the log-K2 score.
	return score;
}

template<class PhenoType>
double
BayesianModel<PhenoType>::
log_factorial_(std::size_t n) {
	if (n == 0) {
		return 0.0;
	}
	double result{0.0};
	for (std::size_t i{1}; i <= n; i++) {
		result += std::log(static_cast<double>(i));
	}
	return result;
}

}



#endif /* SRC_MODEL_BAYESIAN_MODEL_IPP_ */
