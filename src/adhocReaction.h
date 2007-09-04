/**
 * Filename     : adhocReaction.h
 * Description  : Classes describing some ad hoc updates
 * Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
 * Created      : September 2007
 * Revision     : $Id:$
 */
#ifndef ADHOCREACTION_H
#define ADHOCREACTION_H

#include"tissue.h"
#include"baseReaction.h"
#include<cmath>

//!Sets positional derivatives to zero for vertices in specified region  
class VertexNoUpdateFromPosition : public BaseReaction {
  
 public:
  
  VertexNoUpdateFromPosition(std::vector<double> &paraValue, 
			     std::vector< std::vector<size_t> > 
			     &indValue );
  
  void derivs(Tissue &T,
	      std::vector< std::vector<double> > &cellData,
	      std::vector< std::vector<double> > &wallData,
	      std::vector< std::vector<double> > &vertexData,
	      std::vector< std::vector<double> > &cellDerivs,
	      std::vector< std::vector<double> > &wallDerivs,
	      std::vector< std::vector<double> > &vertexDerivs );
};

class VertexTranslateToMax : public BaseReaction {
  
public:
  
  VertexTranslateToMax(std::vector<double> &paraValue, 
											 std::vector< std::vector<size_t> > 
											 &indValue );
  
  void derivs(Tissue &T,
							std::vector< std::vector<double> > &cellData,
							std::vector< std::vector<double> > &wallData,
							std::vector< std::vector<double> > &vertexData,
							std::vector< std::vector<double> > &cellDerivs,
							std::vector< std::vector<double> > &wallDerivs,
							std::vector< std::vector<double> > &vertexDerivs );
	
	void update(Tissue &T,
							std::vector< std::vector<double> > &cellData,
							std::vector< std::vector<double> > &wallData,
							std::vector< std::vector<double> > &vertexData,
							double h);
};

#endif //ADHOCREACTION_H
