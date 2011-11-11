#ifndef CELLTIME_H
#define CELLTIME_H

#include "baseReaction.h"
#include "tissue.h"

class CellTimeDerivative : public BaseReaction
{
public:
	CellTimeDerivative(std::vector<double> &paraValue, std::vector< std::vector<size_t> > &indValue);
  
	void derivs(Tissue &T,
		    DataMatrix &cellData,
		    DataMatrix &wallData,
		    DataMatrix &vertexData,
		    DataMatrix &cellDerivs,
		    DataMatrix &wallDerivs,
		    DataMatrix &vertexDerivs);
};

#endif /* CELLTIME_H */
