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
 * @file  snp_analyzer.hpp
 * @brief 
 */

#ifndef SRC_UTIL_SNP_ANALYZER_HPP_
#define SRC_UTIL_SNP_ANALYZER_HPP_

#include "types.hpp"

namespace epi {

class SNPAnalyzer {

public:

	/*!
	 * @brief Constructor.
	 */
	SNPAnalyzer();

	// One getter for each attribute of SNPAnalyzer::SNPInfo_.

	/*!
	 * @brief Returns SNP type.
	 * @param[in] rsid The SNP's rsID.
	 * @return The type of the SNP.
	 */
	options::SNPType snp_type(const std::string & rsid);

	/*!
	 * @brief Returns name of associated gene or pseudo-gene.
	 * @param[in] rsid The SNP's rsID.
	 * @return Name of containing gene for coding SNPs and ? for non-coding genes.
	 */
	std::string associated_gene(const std::string & rsid);

	/*!
	 * @brief Returns locus of associated gene or pseudo-gene.
	 * @param[in] rsid The SNP's rsID.
	 * @return Locus of the (pseudo-)gene returned by epi::SNPAnalyzer::associated_gene().
	 */
	std::string associated_gene_locus(const std::string & rsid);

	/*!
	 * @brief Checks if the SNP is associated to a pseudo-gene.
	 * @param[in] rsid The SNP's rsID.
	 * @return Boolean @p true if the SNP is associated to a pseudo-gene and @p false otherwise.
	 */
	bool is_associated_to_pseudo_gene(const std::string & rsid);

private:

	// The structure where the SNP information is stored.

	struct SNPInfo_ {

		options::SNPType snp_type;

		bool is_associated_to_pseudo_gene;

		std::string associated_gene;

		std::string associated_gene_locus;

		// maybe even more information, e.g. chromosome and location on chromosome
	};

	// The map where the data is stored.

	std::map<std::string, SNPInfo_> snp_data_;

	// Helper function which returns reference to SNPInfo_ object and if necessary calls retrieve_snp_info_.

	const SNPInfo_ & snp_info_(const std::string & rsid);

	// Retrieves SNP information from NCBI Variation Services.

	void retrieve_snp_info_(const std::string & rsid, SNPInfo_ & snp_info) const;

};

}

#include "snp_analyzer.ipp"

#endif /* SRC_UTIL_SNP_ANALYZER_HPP_ */
