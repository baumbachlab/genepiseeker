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
 * @file  util.hpp
 * @brief Declares helper structures used by compare_models.cpp.
 */

#ifndef TEST_MODEL_SRC_UTIL_HPP_
#define TEST_MODEL_SRC_UTIL_HPP_

#include "../../../src/util/types.hpp"

namespace epi {

/*!
 * @brief Command line options for model tests.
 */
struct ModelTestCLIOptions {

	std::string input_directory;

	std::string output_directory;

	options::InputFormat input_format;

	std::vector<std::string> input_files;

	options::PhenoType pheno_type;

	std::size_t num_categories;

	std::vector<SNP> disease_snps;

	std::size_t num_random_solutions;

	std::vector<std::size_t> size_range_random_solutions;

	std::size_t seed;

	std::vector<epi::options::EpistasisModel> excluded_models;

	std::size_t num_threads;

};

/*!
 * @brief Specifications of and results for epistasis models compared in model tests.
 */
struct ModelTestComparedModel {

	options::EpistasisModel model;

	std::string options;

	bool predict;

	std::vector<double> runtimes;

	double mean_runtime;

	double score_disease_snps;

	std::vector<double> scores_random_solutions;

	double p_value_from_score;

	double mean_prediction_time;

	double prediction_score;

	double mean_correlation;
};

}

std::ostream & operator<<(std::ostream & os, const epi::ModelTestComparedModel & model) {
	os << model.model << "[" << model.options << "],";
	os << model.mean_runtime << ",";
	os << -std::log(model.p_value_from_score) << ",";
	os << model.mean_correlation << ",";
	if (model.predict) {
		os << model.mean_prediction_time << ",";
		os << model.prediction_score << "\n";
	}
	else {
		os << "NA,NA\n";
	}
	return os;
}



#endif /* TEST_MODEL_SRC_UTIL_HPP_ */
