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
 * @file regression_model.ipp
 * @brief Definition of epi::RegressionModel.
 */

#ifndef SRC_MODEL_REGRESSION_MODEL_IPP_
#define SRC_MODEL_REGRESSION_MODEL_IPP_


namespace epi {

template<class PhenoType>
RegressionModel<PhenoType>::
RegressionModel(Instance<PhenoType> * instance) :
EpistasisModel<PhenoType>(instance),
score_("NLL"),
learning_rate_{0.1},
epsilon_{0.01},
max_itrs_{500},
snp_set_(),
interaction_model_() {}

template<class PhenoType>
RegressionModel<PhenoType>::
~RegressionModel() {}

template<class PhenoType>
options::ModelSense
RegressionModel<PhenoType>::
model_sense_() const {
	if (gain_score_() or score_ == "LLH") {
		return options::ModelSense::MAXIMIZE;
	}
	return options::ModelSense::MINIMIZE;
}

template<class PhenoType>
bool
RegressionModel<PhenoType>::
is_predictive_() const {
	return true;
}

template<>
QuantitativePhenoType
RegressionModel<QuantitativePhenoType>::
predict_(Ind ind) const {
	Eigen::MatrixXd interaction_features;
	construct_features_(snp_set_, ind, true, interaction_features);
	return (interaction_features * interaction_model_)(0, 0);
}

template<>
CategoricalPhenoType
RegressionModel<CategoricalPhenoType>::
predict_(Ind ind) const {
	Eigen::MatrixXd interaction_features;
	construct_features_(snp_set_, ind, true, interaction_features);
	if (this->instance_->num_categories() == 2) {
		if (1.0 / (1.0 + std::exp(-(interaction_features * interaction_model_)(0, 0))) < 0.5) {
			return 0;
		}
		return 1;
	}
	Eigen::MatrixXd interaction_hypothesis(1, this->instance_->num_categories());
	interaction_hypothesis = (interaction_features * interaction_model_).array().exp().matrix();
	interaction_hypothesis /= interaction_hypothesis.sum();
	CategoricalPhenoType prediction{0};
	interaction_hypothesis.row(0).maxCoeff(&prediction);
	return prediction;
}

template<class PhenoType>
bool
RegressionModel<PhenoType>::
parse_option_(const std::string & option, const std::string & arg) {
	if (option == "score") {
		score_ = arg;
		if (score_ != "LLH" and score_ != "NLL" and score_ != "AIC" and score_ != "BIC" and not gain_score_()) {
			throw Error(std::string("Invalid argument \"") + arg  + "\" for option \"--" + option + "\". Usage: options = \"[--" + option + " LLLH|LLH-GAIN|NLL|NLL-GAIN|AIC|AIC-GAIN|BIC|BIC-GAIN] [...]\"");
		}
		return true;
	}
	if (option == "max-itrs") {
		try {
			max_itrs_ = std::stoul(arg);
		}
		catch (...) {
			throw Error(std::string("Invalid argument \"") + arg  + "\" for option \"--" + option + "\". Usage: options = \"[--" + option + " <convertible to int greater equal 0>] [...]\"");
		}
		if (max_itrs_ == 0) {
			max_itrs_ = std::numeric_limits<std::size_t>::max();
		}
		return true;
	}
	if (option == "learning-rate") {
		try {
			learning_rate_ = std::stod(arg);
		}
		catch (...) {
			throw Error(std::string("Invalid argument \"") + arg  + "\" for option \"--" + option + "\". Usage: options = \"[--" + option + " <convertible to double greater 0>] [...]\"");
		}
		if (learning_rate_ <= 0) {
			throw Error(std::string("Invalid argument \"") + arg  + "\" for option \"--" + option + "\". Usage: options = \"[--" + option + " <convertible to double greater 0>] [...]\"");
		}
		return true;
	}
	if (option == "epsilon") {
		try {
			epsilon_ = std::stod(arg);
		}
		catch (...) {
			throw Error(std::string("Invalid argument \"") + arg  + "\" for option \"--" + option + "\". Usage: options = \"[--" + option + " <convertible to double greater 0>] [...]\"");
		}
		if (epsilon_ <= 0) {
			throw Error(std::string("Invalid argument \"") + arg  + "\" for option \"--" + option + "\". Usage: options = \"[--" + option + " <convertible to double greater 0>] [...]\"");
		}
		return true;
	}
	return false;
}

template<class PhenoType>
void
RegressionModel<PhenoType>::
set_default_options_() {
	score_ = "NLL";
	learning_rate_ = 0.1;
	epsilon_ = 0.01;
	max_itrs_ = 500;
}


template<class PhenoType>
std::string
RegressionModel<PhenoType>::
valid_options_() const {
	return "[--score LLH|LLH-GAIN|NLL|NLL-GAIN|AIC|AIC-GAIN|BIC|BIC-GAIN] [--max-itrs <convertible to int greater equal 0>] [--learning-rate <convertible to double greater 0>] [--epsilon <convertible to double greater 0>]";
}

template<class PhenoType>
void
RegressionModel<PhenoType>::
save_(const std::string & filename) const {
	std::ofstream file(filename.c_str());
	file << "[meta]\n";
	file << "model-type = RegressionModel\n";
	file << meta_info_();
	file << "snp-set = ";
	std::size_t pos_next_item{1};
	for (SNP snp : snp_set_) {
		file << snp;
		if (pos_next_item++ < snp_set_.size()) {
			file << ",";
		}
	}
	file << "\n";
	file << "[parameters]\n";
	file << "num-features = " << interaction_model_.rows() << "\n";
	file << "num-parameters-per-feature = " << interaction_model_.cols() << "\n";
	for (std::size_t feature_id{0}; feature_id < static_cast<std::size_t>(interaction_model_.rows()); feature_id++) {
		file << feature_id << " = ";
		for (std::size_t param_id{0}; param_id < static_cast<std::size_t>(interaction_model_.cols()); param_id++) {
			file << interaction_model_(feature_id, param_id);
			if (param_id < static_cast<std::size_t>(interaction_model_.cols()) - 1) {
				file << ",";
			}
		}
		file << "\n";
	}
	file.close();
}

template<class PhenoType>
void
RegressionModel<PhenoType>::
load_(const std::string & filename) {
	boost::property_tree::ptree pt;
	boost::property_tree::ini_parser::read_ini(filename, pt);
	if (pt.get<std::string>("meta.model-type") != "RegressionModel") {
		throw Error(std::string("Unexpected model type ") + pt.get<std::string>("meta.model-type") + ".Expected: RegressionModel.");
	}
	std::vector<std::string> tokens;
	misc::tokenize(pt.get<std::string>("meta.snp-set"), ',', '\'', true, tokens);
	snp_set_.clear();
	for (const auto & snp_as_string : tokens) {
		snp_set_.emplace_back(std::stoul(snp_as_string));
	}
	interaction_model_.resize(pt.get<std::size_t>("parameters.num-features"), pt.get<std::size_t>("parameters.num-parameters-per-feature"));
	for (std::size_t feature_id{0}; feature_id < static_cast<std::size_t>(interaction_model_.rows()); feature_id++) {
		misc::tokenize(pt.get<std::string>(std::string("parameters.") + std::to_string(feature_id)), ',', '\'', true, tokens);
		for (std::size_t param_id{0}; param_id < static_cast<std::size_t>(interaction_model_.cols()); param_id++) {
			interaction_model_(feature_id, param_id) = std::stod(tokens.at(param_id));
		}
	}
	check_loaded_model_();
}

template<class PhenoType>
double
RegressionModel<PhenoType>::
evaluate_(const std::vector<SNP> & snp_set, bool prepare_prediction) {

	// Compute interaction model.
	Eigen::MatrixXd interaction_features;
	Eigen::MatrixXd interaction_model;
	Eigen::MatrixXd interaction_hypothesis;
	double interaction_sigma{0};
	construct_feature_matrix_(snp_set, true, interaction_features);
	fit_model_(interaction_features, interaction_model, interaction_hypothesis, interaction_sigma);

	// Compute additive model if required.
	Eigen::MatrixXd additive_features;
	Eigen::MatrixXd additive_model;
	Eigen::MatrixXd additive_hypothesis;
	double additive_sigma{0};
	if (gain_score_()) {
		construct_feature_matrix_(snp_set, false, additive_features);
		fit_model_(additive_features, additive_model, additive_hypothesis, additive_sigma);
	}

	// Prepare prediction if required.
	if (prepare_prediction) {
		snp_set_ = snp_set;
		interaction_model_ = interaction_model;
	}

	// If the likelihood (LLH) is the selected score, compute and return it.
	if (score_ == "LLH" or score_ == "LLH-GAIN") {
		double interaction_score{1.0};
		double additive_score{1.0};
		for (Ind ind{0}; ind < this->instance_->num_inds(); ind++) {
			interaction_score *= likelihood_(ind, interaction_hypothesis, interaction_sigma);
			if (gain_score_()) {
				additive_score *= likelihood_(ind, additive_hypothesis, additive_sigma);
			}
		}
		if (gain_score_()) {
			return interaction_score - additive_score;
		}
		return interaction_score;
	}

	// Compute the negative log-likelihood.
	double interaction_score{0.0};
	double additive_score{0.0};
	for (Ind ind{0}; ind < this->instance_->num_inds(); ind++) {
		interaction_score -= std::log(likelihood_(ind, interaction_hypothesis, interaction_sigma));
		if (gain_score_()) {
			additive_score -= std::log(likelihood_(ind, additive_hypothesis, additive_sigma));
		}
	}

	// If the negative log-likelihood (NLL) or the NLL gain is the selected score, return it.
	if (score_ == "NLL") {
		return interaction_score;
	}
	if (score_ == "NLL-GAIN") {
		return additive_score - interaction_score;
	}

	// Compute the degrees of freedom of the model, i.e., the number of estimated parameters.
	double interaction_degrees_of_freedom{static_cast<double>(interaction_features.cols() * degrees_of_freedom_per_feature_())};
	double additive_degrees_of_freedom{static_cast<double>(additive_features.cols() * degrees_of_freedom_per_feature_())};
	if (this->instance_->quantitative_phenotypes()) {
		interaction_degrees_of_freedom += 1.0;
		additive_degrees_of_freedom += 1.0;
	}

	// If the Aikake information criterion (AIC) or the AIC gain is the selected score, return it.
	if (score_ == "AIC") {
		return 2.0 * interaction_score + 2.0 * interaction_degrees_of_freedom;
	}
	if (score_ == "AIC-GAIN") {
		return 2.0 * additive_score + 2.0 * additive_degrees_of_freedom - 2.0 * interaction_score + 2.0 * interaction_degrees_of_freedom;
	}

	// If the Bayesian information criterion (BIC) or the BIC gain is the selected score, return it.
	double log_num_inds{std::log(static_cast<double>(this->instance_->num_inds()))};
	if (score_ == "BIC") {
		return 2.0 * interaction_score + log_num_inds * interaction_degrees_of_freedom;
	}
	return 2.0 * additive_score + log_num_inds * additive_degrees_of_freedom - 2.0 * interaction_score + log_num_inds * interaction_degrees_of_freedom;
}

template<class PhenoType>
void
RegressionModel<PhenoType>::
construct_feature_matrix_(const std::vector<SNP> & snp_set, bool construct_interaction_features, Eigen::MatrixXd & feature_matrix) const {
	std::size_t num_features{snp_set.size() + 1};
	if (construct_interaction_features) {
		num_features += (snp_set.size() * (snp_set.size() - 1)) / 2;
	}
	feature_matrix.resize(this->instance_->num_inds(), num_features);
	feature_matrix.col(0).setOnes();
	for (std::size_t snp_pos{0}; snp_pos < snp_set.size(); snp_pos++) {
		for (Ind ind{0}; ind < this->instance_->num_inds(); ind++) {
			feature_matrix(ind, snp_pos + 1) = this->instance_->genotype_at_snp(snp_set.at(snp_pos), ind);
		}
	}
	if (not construct_interaction_features) {
		return;
	}
	std::size_t index{snp_set.size() + 1};
	Eigen::ArrayXd col_1(this->instance_->num_inds());
	Eigen::ArrayXd col_2(this->instance_->num_inds());
	for (std::size_t snp_pos_1{0}; snp_pos_1 < snp_set.size() - 1; snp_pos_1++) {
		col_1 = feature_matrix.col(snp_pos_1 + 1).array();
		for (std::size_t snp_pos_2{snp_pos_1 + 1}; snp_pos_2 < snp_set.size(); snp_pos_2++) {
			col_2 = feature_matrix.col(snp_pos_2 + 1).array();
			feature_matrix.col(index++) = (col_1 * col_2).matrix();
		}
	}
}

template<class PhenoType>
void
RegressionModel<PhenoType>::
construct_features_(const std::vector<SNP> & snp_set, Ind ind, bool construct_interaction_features, Eigen::MatrixXd & features) const {
	std::size_t num_features{snp_set.size() + 1};
	if (construct_interaction_features) {
		num_features += (snp_set.size() * (snp_set.size() - 1)) / 2;
	}
	features.resize(1, num_features);
	features(0,0) = 1;
	for (std::size_t snp_pos{0}; snp_pos < snp_set.size(); snp_pos++) {
		features(0, snp_pos + 1) = this->instance_->genotype_at_snp(snp_set.at(snp_pos), ind);
	}
	if (not construct_interaction_features) {
		return;
	}
	std::size_t index{snp_set.size() + 1};
	double feature_1{0.0};
	double feature_2{0.0};
	for (std::size_t snp_pos_1{0}; snp_pos_1 < snp_set.size() - 1; snp_pos_1++) {
		feature_1 = features(0, snp_pos_1 + 1);
		for (std::size_t snp_pos_2{snp_pos_1 + 1}; snp_pos_2 < snp_set.size(); snp_pos_2++) {
			feature_2 = features(0, snp_pos_2 + 1);
			features(0, index++) = feature_1 * feature_2;
		}
	}
}

template<>
void
RegressionModel<QuantitativePhenoType>::
compute_linear_model_sigma_(const Eigen::MatrixXd & hypothesis, double & sigma) const {
	double error{0.0};
	sigma = 0;
	for (Ind ind{0}; ind < this->instance_->num_inds(); ind++) {
		error = hypothesis(ind, 0) - this->instance_->phenotype(ind);
		sigma += error * error;
	}
	sigma /= static_cast<double>(this->instance_->num_inds());
}

template<>
void
RegressionModel<QuantitativePhenoType>::
fit_linear_model_(const Eigen::MatrixXd & feature_matrix, Eigen::MatrixXd & model, Eigen::MatrixXd & hypothesis) const {
	std::size_t num_features{static_cast<std::size_t>(feature_matrix.cols())};
	std::size_t num_inds{this->instance_->num_inds()};
	model.resize(num_features, 1);
	hypothesis.resize(num_inds, 1);
	double loss{std::numeric_limits<double>::max()};
	double loss_old{loss};
	Eigen::MatrixXd gradient(num_features, 1);
	model.setOnes();
	double error{0};
	for (std::size_t itr{0}; itr < max_itrs_; itr++) {
		loss_old = loss;
		hypothesis = feature_matrix * model;
		loss = 0;
		for (Ind ind{0}; ind < num_inds; ind++) {
			error = hypothesis(ind, 0) - this->instance_->phenotype(ind);
			loss += error * error;
		}
		if (loss_old - loss < epsilon_) {
			break;
		}
		gradient = Eigen::VectorXd::Zero(num_features);
		for (Ind ind{0}; ind < num_inds; ind++) {
			gradient += (hypothesis(ind, 0) - this->instance_->phenotype(ind)) * feature_matrix.row(ind).transpose();
		}
		gradient /= static_cast<double>(num_inds);
		model -= learning_rate_ * gradient;
	}
}

template<>
void
RegressionModel<CategoricalPhenoType>::
fit_logistic_model_(const Eigen::MatrixXd & feature_matrix, Eigen::MatrixXd & model, Eigen::MatrixXd & hypothesis) const {
	std::size_t num_features{static_cast<std::size_t>(feature_matrix.cols())};
	std::size_t num_inds{this->instance_->num_inds()};
	model.resize(num_features, 1);
	hypothesis.resize(num_inds, 1);
	double loss{std::numeric_limits<double>::max()};
	double loss_old{loss};
	Eigen::MatrixXd gradient(num_features, 1);
	model.setOnes();
	for (std::size_t itr{0}; itr < max_itrs_ and loss > epsilon_; itr++) {
		loss_old = loss;
		hypothesis = ((-(feature_matrix * model)).array().exp() + 1.0).inverse().matrix();
		loss = 0;
		for (Ind ind{0}; ind < num_inds; ind++) {
			if (this->instance_->phenotype(ind) == 0) {
				loss -= std::log(1.0 - hypothesis(ind, 0));
			}
			else {
				loss -= std::log(hypothesis(ind, 0));
			}
		}
		if (loss_old - loss < epsilon_) {
			break;
		}
		gradient.setZero();
		for (Ind ind{0}; ind < num_inds; ind++) {
			if (this->instance_->phenotype(ind) == 0) {
				gradient += hypothesis(ind, 0) * feature_matrix.row(ind).transpose();
			}
			else {
				gradient += (hypothesis(ind, 0) - 1.0) * feature_matrix.row(ind).transpose();
			}
		}
		gradient /= static_cast<double>(num_inds);
		model -= learning_rate_ * gradient;
	}
	Eigen::MatrixXd temp(hypothesis);
	hypothesis.resize(num_inds, 2);
	hypothesis.col(0) = Eigen::MatrixXd::Ones(num_inds, 1) - temp.col(0);
	hypothesis.col(1) = temp.col(0);
}

template<>
void
RegressionModel<QuantitativePhenoType>::
fit_model_(const Eigen::MatrixXd & feature_matrix, Eigen::MatrixXd & model, Eigen::MatrixXd & hypothesis, double & sigma) const {
	fit_linear_model_(feature_matrix, model, hypothesis);
	compute_linear_model_sigma_(hypothesis, sigma);
}

template<>
void
RegressionModel<CategoricalPhenoType>::
fit_multinomial_model_(const Eigen::MatrixXd & feature_matrix, Eigen::MatrixXd & model, Eigen::MatrixXd & hypothesis) const {
	std::size_t num_features{static_cast<std::size_t>(feature_matrix.cols())};
	std::size_t num_inds{this->instance_->num_inds()};
	std::size_t num_categories{this->instance_->num_categories()};
	model.resize(num_features, num_categories);
	hypothesis.resize(num_inds, num_categories);
	double loss{std::numeric_limits<double>::max()};
	double loss_old{loss};
	Eigen::MatrixXd gradient(num_features, num_categories);
	model.setOnes();
	for (std::size_t itr{0}; itr < max_itrs_ and loss > epsilon_; itr++) {
		loss_old = loss;
		hypothesis = (feature_matrix * model).array().exp().matrix();
		for (Ind ind{0}; ind < num_inds; ind++) {
			hypothesis.row(ind) /= hypothesis.row(ind).sum();
		}
		loss = 0;
		for (Ind ind{0}; ind < num_inds; ind++) {
			loss -= std::log(hypothesis(ind, this->instance_->phenotype(ind)));
		}
		if (loss_old - loss < epsilon_) {
			break;
		}
		gradient.setZero();
		for (Ind ind{0}; ind < num_inds; ind++) {
			for (std::size_t k{0}; k < num_categories; k++) {
				if (k == static_cast<std::size_t>(this->instance_->phenotype(ind))) {
					gradient.col(k) += (hypothesis(ind, k) - 1.0) * feature_matrix.row(ind).transpose();
				}
				else {
					gradient.col(k) += hypothesis(ind, k) * feature_matrix.row(ind).transpose();
				}
			}
		}
		gradient /= static_cast<double>(num_inds);
		model -= learning_rate_ * gradient;
	}
}

template<>
void
RegressionModel<CategoricalPhenoType>::
fit_model_(const Eigen::MatrixXd & feature_matrix, Eigen::MatrixXd & model, Eigen::MatrixXd & hypothesis, double & sigma) const {
	if (this->instance_->num_categories() == 2) {
		fit_logistic_model_(feature_matrix, model, hypothesis);
	}
	else {
		fit_multinomial_model_(feature_matrix, model, hypothesis);
	}
}

template<>
double
RegressionModel<QuantitativePhenoType>::
likelihood_(Ind ind, const Eigen::MatrixXd & hypothesis, double sigma) const {
	return misc::normal_pdf(this->instance_->phenotype(ind), hypothesis(ind, 0), sigma);
}

template<>
double
RegressionModel<CategoricalPhenoType>::
likelihood_(Ind ind, const Eigen::MatrixXd & hypothesis, double sigma) const {
	return misc::ensure_valid_continuous_probability(hypothesis(ind, this->instance_->phenotype(ind)));
}

template<class PhenoType>
bool
RegressionModel<PhenoType>::
gain_score_() const {
	return score_ == "LLH-GAIN" or score_ == "NLL-GAIN" or score_ == "AIC-GAIN" or score_ == "BIC-GAIN";
}

template<>
std::size_t
RegressionModel<QuantitativePhenoType>::
degrees_of_freedom_per_feature_() const {
	return 1;
}

template<>
std::size_t
RegressionModel<CategoricalPhenoType>::
degrees_of_freedom_per_feature_() const {
	return this->instance_->num_categories() - 1;
}

template<>
std::string
RegressionModel<QuantitativePhenoType>::
meta_info_() const {
	return "phenotype = quantitative\n";
}

template<>
std::string
RegressionModel<CategoricalPhenoType>::
meta_info_() const {
	return std::string("phenotype = categorical\ncategories = ") + std::to_string(this->instance_->num_categories()) + "\n";
}

template<>
void
RegressionModel<QuantitativePhenoType>::
check_loaded_model_() const {
	std::size_t num_features{snp_set_.size() + 1 + (snp_set_.size() * (snp_set_.size() - 1)) / 2};
	if (static_cast<std::size_t>(interaction_model_.rows()) != snp_set_.size() + 1 + (snp_set_.size() * (snp_set_.size() - 1)) / 2) {
		throw Error(std::string("Wrong number of features in loaded model for SNP set of size ") + std::to_string(snp_set_.size()) + ". Expected: " + std::to_string(num_features) + ". Actual: " + std::to_string(interaction_model_.rows()) + ".");
	}
	if (interaction_model_.cols() != 1) {
		throw Error(std::string("Wrong number of parameters per feature. Expected: 1. Actual: ") + std::to_string(interaction_model_.cols()) + ".");
	}
}

template<>
void
RegressionModel<CategoricalPhenoType>::
check_loaded_model_() const {
	std::size_t num_features{snp_set_.size() + 1 + (snp_set_.size() * (snp_set_.size() - 1)) / 2};
	if (static_cast<std::size_t>(interaction_model_.rows()) != snp_set_.size() + 1 + (snp_set_.size() * (snp_set_.size() - 1)) / 2) {
		throw Error(std::string("Wrong number of features in loaded model for SNP set of size ") + std::to_string(snp_set_.size()) + ". Expected: " + std::to_string(num_features) + ". Actual: " + std::to_string(interaction_model_.rows()) + ".");
	}
	if (this->instance_->num_categories() == 2) {
		if (interaction_model_.cols() != 1) {
			throw Error(std::string("Wrong number of parameters per feature. Expected: 1. Actual: ") + std::to_string(interaction_model_.cols()) + ".");
		}
	}
	else if (static_cast<std::size_t>(interaction_model_.cols()) != this->instance_->num_categories()) {
		throw Error(std::string("Wrong number of parameters per feature. Expected: ") + std::to_string(this->instance_->num_categories()) + ". Actual: " + std::to_string(interaction_model_.cols()) + ".");
	}
}



}


#endif /* SRC_MODEL_REGRESSION_MODEL_IPP_ */
