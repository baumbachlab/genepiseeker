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
 * @file  instance.ipp
 * @brief Definition of epi::Instance.
 */

#ifndef SRC_MODEL_INSTANCE_IPP_
#define SRC_MODEL_INSTANCE_IPP_

namespace epi {

template<class PhenoType>
Instance<PhenoType>::
const_snp_iterator::
const_snp_iterator(std::vector<GenoType>::const_iterator itr, std::size_t offset, std::size_t pos):
itr_(itr),
offset_{offset},
pos_{pos} {}

template<class PhenoType>
typename Instance<PhenoType>::const_snp_iterator
Instance<PhenoType>::
const_snp_iterator::
operator++() {
	itr_ += offset_;
	pos_ += offset_;
	return *this;
}

template<class PhenoType>
typename  Instance<PhenoType>::const_snp_iterator
Instance<PhenoType>::
const_snp_iterator::
operator++(int) {
	const_snp_iterator temp(*this);
	itr_ += offset_;
	pos_ += offset_;
	return temp;
}

template<class PhenoType>
const typename Instance<PhenoType>::const_snp_iterator::value_type &
Instance<PhenoType>::
const_snp_iterator::
operator*() {
	return *itr_;
}

template<class PhenoType>
bool
Instance<PhenoType>::
const_snp_iterator::
operator==(const const_snp_iterator & rhs) {
	return pos_ == rhs.pos_;
}

template<class PhenoType>
bool
Instance<PhenoType>::
const_snp_iterator::
operator!=(const const_snp_iterator & rhs) {
	return pos_ != rhs.pos_;
}

template<class PhenoType>
Instance<PhenoType>::
Instance(std::size_t num_categories):
num_categories_{num_categories},
num_snps_{undefined_uint()},
num_inds_{undefined_uint()},
genotypes_(),
phenotypes_(),
original_phenotypes_(),
rs_ids_(),
disease_snps_(),
urng_() {}

template<class PhenoType>
void
Instance<PhenoType>::
load(options::InputFormat input_format, const std::string & filename, std::size_t num_folds, std::size_t fold_id, options::DataPurpose cv_data) {
	genotypes_.clear();
	phenotypes_.clear();
	original_phenotypes_.clear();
	rs_ids_.clear();
	switch (input_format) {
	case options::InputFormat::CSV_SNPS_AS_ROWS_FIRST:
		load_csv_(filename, true, true, num_folds, fold_id, cv_data);
		break;
	case options::InputFormat::CSV_SNPS_AS_ROWS_LAST:
		load_csv_(filename, true, false, num_folds, fold_id, cv_data);
		break;
	case options::InputFormat::CSV_SNPS_AS_COLUMNS_FIRST:
		load_csv_(filename, false, true, num_folds, fold_id, cv_data);
		break;
	case options::InputFormat::CSV_SNPS_AS_COLUMNS_LAST:
		load_csv_(filename, false, false, num_folds, fold_id, cv_data);
		break;
	case options::InputFormat::JSON_EPIGEN:
		load_json_(filename, num_folds, fold_id, cv_data);
		break;
	default:
		throw Error("Unsupported input format.");
	}
	original_phenotypes_ = phenotypes_;
}

template<class PhenoType>
std::size_t
Instance<PhenoType>::
num_snps() const {
	return num_snps_;
}

template<class PhenoType>
std::size_t
Instance<PhenoType>::
num_inds() const {
	return num_inds_;
}

template<>
std::size_t
Instance<CategoricalPhenoType>::
num_categories() const {
	return num_categories_;
}

template<class PhenoType>
GenoType
Instance<PhenoType>::
genotype_at_snp(SNP snp, Ind ind) const {
	return genotypes_.at(snp * num_inds_ + ind);
}

template<class PhenoType>
void
Instance<PhenoType>::
genotype_at_snp_set(const std::vector<SNP> & snp_set, Ind ind, std::vector<GenoType> & genotype) const {
	genotype.clear();
	for (SNP snp : snp_set) {
		genotype.emplace_back(genotype_at_snp(snp, ind));
	}
}

template<class PhenoType>
std::size_t
Instance<PhenoType>::
genotype_at_snp_set(const std::vector<SNP> & snp_set, Ind ind) const {
	std::size_t genotype_id{0};
	std::size_t exponent{snp_set.size() - 1};
	for (SNP snp : snp_set) {
		genotype_id += genotype_at_snp(snp, ind) * misc::uint_pow(3, exponent--);
	}
	return genotype_id;
}

template<class PhenoType>
void
Instance<PhenoType>::
inds_with_genotype_at_snp_set(const std::vector<SNP> & snp_set, const std::vector<GenoType> & genotype, std::vector<Ind> & inds) const {
	inds.clear();
	for (Ind ind{0}; ind < num_inds_; ind++) {
		bool match{true};
		for (std::size_t pos{0}; pos < snp_set.size(); pos++) {
			if (genotype_at_snp(snp_set.at(pos), ind) != genotype.at(pos)) {
				match = false;
				break;
			}
		}
		if (match) {
			inds.emplace_back(ind);
		}
	}
}

template<class PhenoType>
void
Instance<PhenoType>::
inds_with_genotype_at_snp_set(const std::vector<SNP> & snp_set, std::size_t genotype_id, std::vector<Ind> & inds) const {
	std::vector<GenoType> genotype;
	misc::id_to_genotype(genotype_id, snp_set.size(), genotype);
	return inds_with_genotype_at_snp_set(snp_set, genotype, inds);
}

template<class PhenoType>
std::size_t
Instance<PhenoType>::
num_inds_with_genotype_at_snp_set(const std::vector<SNP> & snp_set, const std::vector<GenoType> & genotype) const {
	std::size_t num_matches{0};
	for (Ind ind{0}; ind < num_inds_; ind++) {
		bool match{true};
		for (std::size_t pos{0}; pos < snp_set.size(); pos++) {
			if (genotype_at_snp(snp_set.at(pos), ind) != genotype.at(pos)) {
				match = false;
				break;
			}
		}
		if (match) {
			num_matches++;
		}
	}
	return num_matches;
}

template<class PhenoType>
std::size_t
Instance<PhenoType>::
num_inds_with_genotype_at_snp_set(const std::vector<SNP> & snp_set, std::size_t genotype_id) const {
	std::vector<GenoType> genotype;
	misc::id_to_genotype(genotype_id, snp_set.size(), genotype);
	return num_inds_with_genotype_at_snp_set(snp_set, genotype);
}

template<class PhenoType>
PhenoType
Instance<PhenoType>::
phenotype(Ind ind) const {
	return phenotypes_.at(ind);
}

template<class PhenoType>
typename Instance<PhenoType>::const_ind_iterator
Instance<PhenoType>::
genotypes_of_all_inds_begin(SNP snp) const {
	return (genotypes_.cbegin() + (snp * num_inds_));
}

template<class PhenoType>
typename Instance<PhenoType>::const_ind_iterator
Instance<PhenoType>::
genotypes_of_all_inds_end(SNP snp) const {
	return (genotypes_.cbegin() + ((snp + 1) * num_inds_));
}

template<class PhenoType>
typename Instance<PhenoType>::const_snp_iterator
Instance<PhenoType>::
genotypes_at_all_snps_begin(Ind ind) const {
	return const_snp_iterator(genotypes_.cbegin() + ind, num_inds_, ind);
}

template<class PhenoType>
typename Instance<PhenoType>::const_snp_iterator
Instance<PhenoType>::
genotypes_at_all_snps_end(Ind ind) const {
	return const_snp_iterator(genotypes_.cbegin() + ind, num_inds_, ind + num_inds_ * num_snps_);
}

template<class PhenoType>
void
Instance<PhenoType>::
shuffle_phenotypes() {
	std::shuffle(phenotypes_.begin(), phenotypes_.end(), urng_);
}

template<class PhenoType>
void
Instance<PhenoType>::
restore_phenotypes() {
	phenotypes_ = original_phenotypes_;
}

template<class PhenoType>
void
Instance<PhenoType>::
set_seed(std::size_t seed) {
	urng_.seed(seed);
}

template<>
bool
Instance<QuantitativePhenoType>::
quantitative_phenotypes() const {
	return true;
}

template<>
bool
Instance<CategoricalPhenoType>::
quantitative_phenotypes() const {
	return false;
}

template<class PhenoType>
bool
Instance<PhenoType>::
categorical_phenotypes() const {
	return not quantitative_phenotypes();
}

template<class PhenoType>
const std::vector<SNP> &
Instance<PhenoType>::
disease_snps() const {
	if (disease_snps_.empty()) {
		throw Error("No disease SNP set available.");
	}
	return disease_snps_;
}

template<class PhenoType>
void
Instance<PhenoType>::
set_disease_snps(const std::vector<SNP> & disease_snps) {
	disease_snps_ = disease_snps;
	for (SNP snp : disease_snps_) {
		if (snp >= num_snps_) {
			throw Error("The selected disease SNP " + std::to_string(snp) + " does not exist.");
		}
	}
	if (disease_snps_.empty()) {
		throw Error("The selected disease SNP set is empty.");
	}
	std::set<double> test_unique(disease_snps_.begin(), disease_snps_.end());
	if (test_unique.size() != disease_snps_.size()) {
		throw Error("The selected disease SNP set contains duplicates.");
	}
}

template<class PhenoType>
std::string
Instance<PhenoType>::
snp_descriptor(SNP snp) const {
	if (snp >= num_snps_) {
		throw Error("The selected SNP " + std::to_string(snp) + " does not exist.");
	}
	return rs_ids_.at(snp);
}

template<class PhenoType>
void
Instance<PhenoType>::
load_csv_(const std::string & filename, bool snps_as_rows, bool snp_info_comes_first, std::size_t num_folds, std::size_t fold_id, options::DataPurpose cv_data) {
	CSVParser csv_parser;
	csv_parser.parse(filename);
	num_snps_ = snps_as_rows ? csv_parser.num_rows() - 1 : csv_parser.num_columns() - 1;
	num_inds_ = snps_as_rows ? csv_parser.num_columns() - 1 : csv_parser.num_rows() - 1;
	std::string field;
	std::string index;

	// Determine which individuals should be skipped.
	std::vector<bool> skip;
	std::size_t num_skipped_inds{construct_folds_(num_folds, fold_id, cv_data, skip)};

	// Determine indices of first and of last individual as well as of SNP information.
	Ind first_ind{0};
	Ind last_ind{num_inds_ - 1};
	std::size_t index_snp_info{num_inds_};
	std::size_t skip_ind_shift{0};
	if (snp_info_comes_first) {
		skip_ind_shift = 1;
		first_ind = 1;
		last_ind = num_inds_;
		index_snp_info = 0;
	}

	// Load genotypes.
	for (SNP snp{0}; snp < num_snps_; snp++) {
		for (Ind ind{first_ind}; ind <= last_ind + num_skipped_inds; ind++) {
			if (not skip.at(ind - skip_ind_shift)) {
				field = snps_as_rows ? csv_parser.cell(snp, ind) : csv_parser.cell(ind, snp);
				try {
					genotypes_.emplace_back(static_cast<GenoType>(std::stoul(field)));
				}
				catch (...) {
					throw Error("The input file contains invalid genotype for individual " + std::to_string(ind) + " at SNP " + std::to_string(snp) + ". Expected: 0, 1, or 2.");
				}
				if (genotypes_.back() > 2) {
					throw Error("The input file contains invalid genotype for individual " + std::to_string(ind) + " at SNP " + std::to_string(snp) + ". Expected: 0, 1, or 2.");
				}
			}
		}
	}

	// Load phenotypes.
	for (Ind ind{first_ind}; ind <= last_ind + num_skipped_inds; ind++) {
		if (not skip.at(ind - skip_ind_shift)) {
			field = snps_as_rows ? csv_parser.cell(num_snps_, ind) : csv_parser.cell(ind, num_snps_);
			phenotypes_.emplace_back(parse_phenotype_(field, ind));
		}
	}

	// Load RS IDs.
	for (SNP snp{0}; snp < num_snps_; snp++) {
		field = snps_as_rows ? csv_parser.cell(snp, index_snp_info) : csv_parser.cell(index_snp_info, snp);
		rs_ids_.emplace_back(field);
	}
}

template<class PhenoType>
void
Instance<PhenoType>::
load_json_(const std::string & filename, std::size_t num_folds, std::size_t fold_id, options::DataPurpose cv_data) {
	boost::property_tree::ptree root;
	try {
		boost::property_tree::read_json(filename, root);
	}
	catch (...) {
		throw Error("The file " + filename + " cannot be opened.");
	}
	try {
		num_snps_ = root.get<std::size_t>("num_snps");
	}
	catch (...) {
		throw Error("The input file must contain a field \"num_snps\" whose value must be convertible to an int greater 0.");
	}
	if (num_snps_ <= 0) {
		throw Error("The input file must contain a field \"num_snps\" whose value must be convertible to an int greater 0.");
	}
	try {
		num_inds_ = root.get<std::size_t>("num_inds");
	}
	catch (...) {
		throw Error("The input file must contain a field \"num_inds\" whose value must be convertible to an int greater 0.");
	}
	if (num_inds_ <= 0) {
		throw Error("The input file must contain a field \"num_inds\" whose value must be convertible to an int greater 0.");
	}

	// Determine which individuals should be skipped.
	std::vector<bool> skip;
	std::size_t num_skipped_inds{construct_folds_(num_folds, fold_id, cv_data, skip)};

	// Load genotypes.
	if (root.find("genotype") == root.not_found()) {
		throw Error("The input file must contain a field \"genotype\".");
	}
	SNP snp{0};
	Ind ind{0};
	for (const auto & row : root.get_child("genotype")) {
		ind = 0;
		for (const auto & cell : row.second) {
			if (not skip.at(ind)) {

				try {
					genotypes_.emplace_back(cell.second.get_value<GenoType>());
				}
				catch (...) {
					throw Error("The input file contains invalid genotype for individual " + std::to_string(ind) + " at SNP " + std::to_string(snp) + ". Expected: 0, 1, or 2.");
				}
				if (genotypes_.back() > 2) {
					throw Error("The input file contains invalid genotype for individual " + std::to_string(ind) + " at SNP " + std::to_string(snp) + ". Expected: 0, 1, or 2.");
				}
			}
			ind++;
		}
		if (ind != num_inds_ + num_skipped_inds) {
			throw Error("The actual number of individuals " + std::to_string(ind) + " for SNP " + std::to_string(snp) + " in the \"genotype\" field does not match the number of individuals " + std::to_string(num_inds_ + num_skipped_inds) + " specified in the \"num_inds\" field.");
		}
		snp++;
	}
	if (snp != num_snps_) {
		throw Error("The actual number of SNPs in the \"genotype\" field  does not match the specified number of SNPs specified in the \"num_snps\" field.");
	}

	// Load phenotypes.
	if (root.find("phenotype") == root.not_found()) {
		throw Error("The input file must contain a field \"phenotype\".");
	}
	ind = 0;
	for (const auto & phenotype : root.get_child("phenotype")) {
		if (not skip.at(ind)) {
			phenotypes_.emplace_back(parse_phenotype_(phenotype.second.get_value<std::string>(), ind));
		}
		ind++;
	}
	if (ind != num_inds_ + num_skipped_inds) {
		throw Error("The actual number of individuals in the \"phenotype\" field  does not match the specified number of individuals specified in the \"num_inds\" field.");
	}

	// Load RS IDs.
	if (root.find("snps") == root.not_found()) {
		throw Error("The input file must contain a field \"snps\".");
	}
	snp = 0;
	for (const auto & snp_info : root.get_child("snps")) {
		for (const auto & cell : snp_info.second) {
			rs_ids_.emplace_back(cell.second.get_value<std::string>());
			break;
		}
		snp++;
	}
	if (snp != num_snps_) {
		throw Error("The actual number of SNPs in the \"snps\" field  does not match the specified number of SNPs specified in the \"num_snps\" field.");
	}

	// Load disease SNPs (if available).
	if (root.find("disease_snps") == root.not_found()) {
		return;
	}
	disease_snps_.clear();
	for (const auto & disease_snp : root.get_child("disease_snps")) {
		disease_snps_.emplace_back(disease_snp.second.get_value<SNP>());
		if (disease_snps_.back() >= num_snps_) {
			throw Error("The selected disease SNP " + std::to_string(disease_snps_.back()) + " does not exist");
		}
	}
	std::set<double> test_unique(disease_snps_.begin(), disease_snps_.end());
	if (test_unique.size() != disease_snps_.size()) {
		throw Error("The selected disease SNP set contains duplicates.");
	}
}

template<>
CategoricalPhenoType
Instance<CategoricalPhenoType>::
parse_phenotype_(const std::string & phenotype_str, Ind ind) const {
	CategoricalPhenoType phenotype;
	try {
		phenotype = static_cast<CategoricalPhenoType>(std::stoul(phenotype_str));
	}
	catch (...) {
		throw Error("The input file contains invalid phenotype for individual " + std::to_string(ind) + ". Expected: convertible to int between 0 and " + std::to_string(num_categories_ - 1) + ".");
	}
	if (phenotype >= num_categories_) {
		throw Error("The input file contains invalid phenotype for individual " + std::to_string(ind) + ". Expected: convertible to int between 0 and " + std::to_string(num_categories_ - 1) + ".");
	}
	return phenotype;
}

template<>
QuantitativePhenoType
Instance<QuantitativePhenoType>::
parse_phenotype_(const std::string & phenotype_str, Ind ind) const {
	QuantitativePhenoType phenotype;
	try {
		phenotype = static_cast<QuantitativePhenoType>(std::stod(phenotype_str));
	}
	catch (...) {
		throw Error("The input file contains invalid phenotype for individual " + std::to_string(ind) + ". Expected: convertible to double.");
	}
	return phenotype;
}

template<class PhenoType>
std::size_t
Instance<PhenoType>::
construct_folds_(std::size_t num_folds, std::size_t fold_id, options::DataPurpose cv_data, std::vector<bool> & skip) {
	if (num_folds > num_inds_) {
		throw Error("The specified number of folds exceeds the number of individuals in the input data.");
	}
	if (fold_id >= num_folds) {
		throw Error("The specified fold ID has to be smaller than the specified number of folds.");
	}
	bool skip_most_inds{cv_data == options::DataPurpose::TRAINING or num_folds == 1 ? false : true};
	skip = std::vector<bool>(num_inds_, skip_most_inds);
	if (num_folds == 1) {
		return 0;
	}
	std::size_t size_fold{num_inds_ / num_folds};
	std::size_t remainder{num_inds_ % num_folds};
	Ind ind_in_fold{fold_id * size_fold + std::min(fold_id, remainder)};
	Ind last_ind_in_fold{(fold_id + 1) * size_fold - 1 + std::min(fold_id + 1, remainder)};
	std::size_t num_inds_in_fold{0};
	while (ind_in_fold <= last_ind_in_fold) {
		skip[ind_in_fold++] = not skip_most_inds;
		num_inds_in_fold++;
	}
	std::size_t num_skipped_inds{skip_most_inds ? num_inds_ - num_inds_in_fold : num_inds_in_fold};
	if (skip_most_inds) {
		num_inds_ = num_inds_in_fold;
	}
	else {
		num_inds_ -= num_inds_in_fold;
	}
	return num_skipped_inds;
}


}



#endif /* SRC_MODEL_INSTANCE_IPP_ */
