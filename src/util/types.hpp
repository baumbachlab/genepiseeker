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
 * @file types.hpp
 * @brief Declarations and inclusions of types used by various classes.
 */

#ifndef SRC_UTIL_TYPES_HPP_
#define SRC_UTIL_TYPES_HPP_

// Include standard libraries.
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <cctype>
#include <map>
#include <set>
#include <stdexcept>
#include <algorithm>
#include <type_traits>
#include <cmath>
#include <random>
#include <functional>
#include <dirent.h>

// Include Boost libraries.
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/fisher_f.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/students_t.hpp>
#include <boost/filesystem.hpp>

// Include OpenMP.
#ifdef _OPENMP
#include <omp.h>
#endif

// Include Eigen.
#include <Eigen/Dense>

/*!
 * @brief Global namespace for GenEpiSeeker.
 */
namespace epi {

// Type definitions.
/*!
 * @brief Type of SNPs.
 */
typedef std::size_t SNP;

/*!
 * @brief Type of individuals.
 */
typedef std::size_t Ind;

/*!
 * @brief Type of entries in the genotype matrix.
 */
typedef std::uint_fast8_t GenoType;

/*!
 * @brief Use this type if you want to use GenEpiSeeker for quantitative phenotypes.
 */
typedef double QuantitativePhenoType;

/*!
 * @brief Use this type if you want to use GenEpiSeeker for categorical phenotypes.
 */
typedef std::uint_fast8_t CategoricalPhenoType;

/*!
 * @brief Type for reporting runtimes in seconds.
 */
typedef std::chrono::duration<double> Seconds;

// Constant expressions.
/*!
 * @brief Used to denote undefined unsigned integers.
 * @return An unsigned integer used to denote undefined integers.
 */
constexpr std::size_t undefined_uint() {return std::numeric_limits<std::size_t>::max();}

/*!
 * @brief Used to denote undefined doubles.
 * @return A double used to denote undefined doubles.
 */
constexpr double undefined_double() {return std::numeric_limits<double>::infinity();}

/*!
 * @brief The constant @f$\pi@f$.
 * @return The constant @f$\pi@f$.
 */
constexpr double pi() {return 3.14159265358979323846;}

// Error class.
/*!
 * @brief Exceptions thrown by GenEpiSeeker are of this type.
 */
struct Error : public std::runtime_error {
	/*!
	 * @brief Constructor.
	 * @param[in] message Error message
	 */
	Error(const std::string & message) : std::runtime_error(message) {}
};

/*!
 * @brief Contains enum classes used to specify options.
 */
namespace options {

/*!
 * @brief Specifies the SNP type.
 */
enum class SNPType {
	NON_CODING,       //!< Non-coding SNPs.
	CODING_SYNONYMOUS,//!< Coding, synonymous SNPs.
	CODING_MISSENSE,  //!< Coding, missense SNPs.
	CODING_NONSENSE   //!< Coding, nonsense SNPs.
};

/*!
 * @brief Specifies the input format.
 */
enum class InputFormat {
	JSON_EPIGEN,              //!< Input is JSON file in the format produced by EpiGEN.
	CSV_SNPS_AS_ROWS_FIRST,   //!< Input is CSV file where the SNPs are represented by the rows and the first column contains SNP information.
	CSV_SNPS_AS_ROWS_LAST,    //!< Input is CSV file where the SNPs are represented by the rows and the last column contains SNP information.
	CSV_SNPS_AS_COLUMNS_FIRST,//!< Input is CSV file where the SNPs are represented by the columns and the first row contains SNP information.
	CSV_SNPS_AS_COLUMNS_LAST  //!< Input is CSV file where the SNPs are represented by the columns and the last row contains SNP information.
};

/*!
 * @brief Selects the epistasis model.
 */
enum class EpistasisModel {
	BAYESIAN_MODEL,  //!< Selects epi::BayesianModel.
	PENETRANCE_MODEL,//!< Selects epi::PenetranceModel.
	REGRESSION_MODEL,//!< Selects epi::RegressionModel.
	VARIANCE_MODEL   //!< Selects epi::VarianceModel.
};

/*!
 * @brief Specifies the sense of epistasis models a.k.a. objective functions.
 */
enum class ModelSense {
	MINIMIZE,//!< The objective should be minimized.
	MAXIMIZE //!< The objective should be maximized.
};

/*!
 * @brief Specified whether the phenotypes are quantitative or categorical.
 */
enum class PhenoType {
	QUANTITATIVE,//!< Quantitative phenotypes.
	CATEGORICAL  //!< Categorical phenotypes.
};

/*!
 * @brief Specifies whether the data loaded into the instance should be used for training or validation.
 */
enum class DataPurpose {
	TRAINING, //!< The data should be used for training.
	VALIDATION//!< The data should be used for validation.
};

/*!
 * @brief Specifies direction of statistical test.
 */
enum class TestDirection {
	TWO_TAILED,  //!< Two tailed test, i.e. p-value = 2 * Pr(D > |test_statistic|), with D distributed according null hypothesis.
	LOWER_TAILED,//!< Lower tailed test, i.e. p-value = Pr(D < test_statistic), with D distributed according null hypothesis
	UPPER_TAILED //!< Upper tailed test, i.e. p-value = Pr(D > test_statistic), with D distributed according null hypothesis
};

}

}

/*!
 * @brief Streams epi::Options::EpistasisModel.
 * @param[in] os Output stream.
 * @param[in] model Epistasis model selector.
 * @return Output stream.
 */
std::ostream & operator<<(std::ostream & os, const epi::options::EpistasisModel & model) {
	switch (model) {
	case epi::options::EpistasisModel::BAYESIAN_MODEL:
		os << "BAYESIAN_MODEL";
		break;
	case epi::options::EpistasisModel::PENETRANCE_MODEL:
		os << "PENETRANCE_MODEL";
		break;
	case epi::options::EpistasisModel::REGRESSION_MODEL:
		os << "REGRESSION_MODEL";
		break;
	case epi::options::EpistasisModel::VARIANCE_MODEL:
		os << "VARIANCE_MODEL";
		break;
	}
	return os;
}


#endif /* SRC_UTIL_TYPES_HPP_ */
