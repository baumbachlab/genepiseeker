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
 * @file  snp_analyzer.ipp
 * @brief 
 */

#ifndef SRC_UTIL_SNP_ANALYZER_IPP_
#define SRC_UTIL_SNP_ANALYZER_IPP_

namespace epi {

SNPAnalyzer::
SNPAnalyzer():
snp_data_() {}

options::SNPType
SNPAnalyzer::
snp_type(const std::string & rsid) {
	return snp_info_(rsid).snp_type;
}

std::string
SNPAnalyzer::
associated_gene(const std::string & rsid) {
	return snp_info_(rsid).associated_gene;
}

std::string
SNPAnalyzer::
associated_gene_locus(const std::string & rsid) {
	return snp_info_(rsid).associated_gene_locus;
}

bool
SNPAnalyzer::
is_associated_to_pseudo_gene(const std::string & rsid) {
	return snp_info_(rsid).is_associated_to_pseudo_gene;
}

const SNPAnalyzer::SNPInfo_ &
SNPAnalyzer::
snp_info_(const std::string & rsid) {
	std::string rsid_digits(rsid);
	if (rsid.substr(0, 2) == "rs") {
		rsid_digits = rsid.substr(2);
	}
	if (rsid_digits.empty() or not std::all_of(rsid_digits.begin(), rsid_digits.end(), ::isdigit)) {
		throw Error("Invalid rsID " + rsid + ".");
	}
	auto itr = snp_data_.find(rsid_digits);
	if (itr == snp_data_.end()) {
		SNPInfo_ snp_info;
		retrieve_snp_info_(rsid_digits, snp_info);
		itr = snp_data_.emplace(rsid_digits, snp_info).first;
	}
	return itr->second;
}

void
SNPAnalyzer::
retrieve_snp_info_(const std::string & rsid, SNPInfo_ & snp_info) const {
	// todo: retrieve SNP information from NCBI and save it in snp_info
}

}



#endif /* SRC_UTIL_SNP_ANALYZER_IPP_ */
