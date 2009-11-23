#include "cellTime.h"
#include <cstdlib>

CellTimeDerivative::CellTimeDerivative(std::vector<double> &paraValue, 
	std::vector< std::vector<size_t> > 
	&indValue)
{
	if (paraValue.size() != 0) {
		std::cerr << "CellTimeDerivative::CellTimeDerivative() uses no parameters.\n";
		std::exit(EXIT_FAILURE);
	}

	if (indValue.size() != 1 || indValue[0].size() != 1) {
		std::cerr << "CellTimeDerivative::CellTimeDerivative(): "
		<< "First index of first level sets cell time index.\n";
		std::exit(EXIT_FAILURE);
	}

	setId("CellTimeDerivative");
	setParameter(paraValue);  
	setVariableIndex(indValue);
  
	std::vector<std::string> tmp(numParameter());
	setParameterId( tmp );
}

void CellTimeDerivative::derivs(Tissue &T,
	std::vector< std::vector<double> > &cellData,
	std::vector< std::vector<double> > &wallData,
	std::vector< std::vector<double> > &vertexData,
	std::vector< std::vector<double> > &cellDerivs,
	std::vector< std::vector<double> > &wallDerivs,
	std::vector< std::vector<double> > &vertexDerivs)
{
	for (size_t index = 0; index < T.numCell(); ++index)
	{
		Cell &cell = T.cell(index);
		
		cellDerivs[cell.index()][variableIndex(0, 0)] = 1.0;
	}
}
