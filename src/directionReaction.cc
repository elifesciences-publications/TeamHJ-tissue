///
/// Filename     : directionReaction.cc
/// Description  : Classes describing some reaction updates related to directions
/// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
/// Created      : May 2008
/// Revision     : $Id:$
///

#include"directionReaction.h"
#include"tissue.h"
#include"baseReaction.h"
#include<cmath>

ContinousMTDirection::ContinousMTDirection(std::vector<double> &paraValue,
																					 std::vector< std::vector<size_t> > &indValue)
{
	if (paraValue.size() != 1) {
		std::cerr << "ContinousMTDirection::ContinousMTDirection() " 
							<< "Uses one parameter: k_rate" << std::endl;
		exit(EXIT_FAILURE);
	}
	
	if (indValue.size() != 2 || indValue[0].size() != 1 ||  indValue[1].size() != 1) {
		std::cerr << "ContinousMTDirection::ContinousMTDirection() " << std::endl
							<< "First level gives target direction index (input)." << std::endl
							<< "Second level gives real direction index." << std::endl;
		exit(EXIT_FAILURE);
	}
	
	setId("ContinousMTDirection");
	setParameter(paraValue);
	setVariableIndex(indValue);
	
	std::vector<std::string> tmp(numParameter());
	tmp[0] = "k_rate";
	
	setParameterId(tmp);
}

void ContinousMTDirection::derivs(Tissue &T,
																	std::vector< std::vector<double> > &cellData,
																	std::vector< std::vector<double> > &wallData,
																	std::vector< std::vector<double> > &vertexData,
																	std::vector< std::vector<double> > &cellDerivs,
																	std::vector< std::vector<double> > &wallDerivs,
																	std::vector< std::vector<double> > &vertexDerivs)
{
	size_t dimension=vertexData[0].size();
	if (dimension!=2) {
		std::cerr << "ContinuosMTDirection::derivs() Only implemented for two dimensions." << std::endl;
		exit(EXIT_FAILURE);
	}
	size_t target = variableIndex(0, 0);
	size_t real = variableIndex(1, 0);
	
	double k_rate = parameter(0);
	
	for (size_t n = 0; n < T.numCell(); ++n) {
		Cell cell = T.cell(n);
		size_t index = cell.index();
		
		double x = cellData[index][real + 0];
		double y = cellData[index][real + 1];
		
		double sigma = std::atan2(y, x);
		
		while (sigma > 0.5 * M_PI || sigma <= -0.5 * M_PI) {
			if (sigma > 0.5 * M_PI) {
				sigma -= M_PI;
			}
			if (sigma <= -0.5 * M_PI) {
				sigma += M_PI;
			}
		}
		
		double dx = cellData[index][target + 0];
		double dy = cellData[index][target + 1];
		
		double dsigma = std::atan2(dy, dx);
		
		while (dsigma > 0.5 * M_PI || dsigma <= -0.5 * M_PI) {
			if (dsigma > 0.5 * M_PI) {
				dsigma -= M_PI;
			}
			if (dsigma <= -0.5 * M_PI) {
				dsigma += M_PI;
			}
		}
		
		double angle = dsigma - sigma;
		
		while (angle > M_PI || angle <= -M_PI) {
			if (angle > M_PI) {
				angle -= 2.0 * M_PI;
			}
			if (angle <= -M_PI) {
				angle += 2.0 * M_PI;
			}
		}
		
		double speed = k_rate * std::abs(angle) / (0.25 * M_PI + std::abs(angle));
		
		speed *= (angle >= 0) ? +1 : -1;
		
		cellDerivs[index][real + 0] += -y * speed;
		cellDerivs[index][real + 1] += x * speed;
	}
}

UpdateMTDirection::UpdateMTDirection(std::vector<double> &paraValue,
																		 std::vector< std::vector<size_t> > &indValue)
{
	if (paraValue.size() != 1) {
		std::cerr << "UpdateMTDirection::UpdateMTDirection() " 
							<< "Uses one parameter: k_rate" << std::endl;
		exit(EXIT_FAILURE);
	}
	
	if (indValue.size() != 2 || indValue[0].size() != 1 ||  indValue[1].size() != 1) {
		std::cerr << "UpdateMTDirection::UpdateMTDirection() " << std::endl
							<< "First level gives target direction index (input)." << std::endl
							<< "Second level gives real direction index." << std::endl;
		exit(EXIT_FAILURE);
	}
	
	setId("UpdateMTDirection");
	setParameter(paraValue);
	setVariableIndex(indValue);
	
	std::vector<std::string> tmp(numParameter());
	tmp[0] = "k_rate";
	
	setParameterId(tmp);
}

void UpdateMTDirection::initiate(Tissue &T,
																 std::vector< std::vector<double> > &cellData,
																 std::vector< std::vector<double> > &wallData,
																 std::vector< std::vector<double> > &vertexData) 
{
	size_t numCell=cellData.size();
	size_t dimension=vertexData[0].size();
	size_t inIndex=variableIndex(0,0);
	size_t outIndex=variableIndex(1,0);
	for (size_t i=0; i<numCell; ++i)
		for (size_t d=0; d<dimension; ++d)
			cellData[i][outIndex+d] = cellData[i][inIndex+d];
}

void UpdateMTDirection::derivs(Tissue &T,
															 std::vector< std::vector<double> > &cellData,
															 std::vector< std::vector<double> > &wallData,
															 std::vector< std::vector<double> > &vertexData,
															 std::vector< std::vector<double> > &cellDerivs,
															 std::vector< std::vector<double> > &wallDerivs,
															 std::vector< std::vector<double> > &vertexDerivs ) {}

void UpdateMTDirection::update(Tissue &T,
															 std::vector< std::vector<double> > &cellData,
															 std::vector< std::vector<double> > &wallData,
															 std::vector< std::vector<double> > &vertexData,
															 double h) 
{
	size_t numCell=cellData.size();
	size_t dimension=vertexData[0].size();
	size_t inIndex=variableIndex(0,0);
	size_t outIndex=variableIndex(1,0);
	for (size_t i=0; i<numCell; ++i) {
		for (size_t d=0; d<dimension; ++d)
			cellData[i][outIndex+d] += parameter(0)*h*(cellData[i][inIndex+d]-cellData[i][outIndex+d]);
		// Normalize
		double norm=0.0;
		for (size_t d=0; d<dimension; ++d)
			norm += cellData[i][outIndex+d]*cellData[i][outIndex+d];
		norm = 1.0/std::sqrt(norm);
		for (size_t d=0; d<dimension; ++d)
			cellData[i][outIndex+d] *= norm;
	}
}
