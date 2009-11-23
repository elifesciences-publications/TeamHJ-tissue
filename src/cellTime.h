#ifndef CELLTIME_H
#define CELLTIME_H

#include "baseReaction.h"
#include "tissue.h"

class CellTimeDerivative : public BaseReaction
{
public:
	CellTimeDerivative(std::vector<double> &paraValue, std::vector< std::vector<size_t> > &indValue);
  
	void derivs(Tissue &T,
		std::vector< std::vector<double> > &cellData,
		std::vector< std::vector<double> > &wallData,
		std::vector< std::vector<double> > &vertexData,
		std::vector< std::vector<double> > &cellDerivs,
		std::vector< std::vector<double> > &wallDerivs,
		std::vector< std::vector<double> > &vertexDerivs);
};

#endif /* CELLTIME_H */
