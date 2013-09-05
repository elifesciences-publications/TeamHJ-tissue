//
// Filename     : bending.cc
// Description  : Classes describing reactions related to bending moments
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : September 2013
// Revision     : $Id:$
//
#include<cmath>
#include"baseReaction.h"
#include"bending.h"
#include"tissue.h"

namespace Bending {
  
  NeighborCenter::
  NeighborCenter(std::vector<double> &paraValue, 
		 std::vector< std::vector<size_t> > 
		 &indValue )
  {
    //Do some checks on the parameters and variable indices
    //
    if( paraValue.size()!=1 ) {
      std::cerr << "Bending::NeighborCenter::"
		<< "NeighborCenter() "
		<< "One parameter, k_bend, should be provided." << std::endl;
      exit(EXIT_FAILURE);
    }
    if( indValue.size() != 1 || indValue[0].size() != 1 ) {
      std::cerr << "VertexFromCellPressure::"
		<< "VertexFromCellPressure() "
		<< "One index level with one index (wall length) given" 
		<< std::endl;
      exit(EXIT_FAILURE);
    }
    //Set the variable values
    //
    setId("Bending::NeighborCenter");
    setParameter(paraValue);  
    setVariableIndex(indValue);
  }
  
  void NeighborCenter::
  derivs(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs )
  {
    size_t numCells = T.numCell();
    size_t dimension = T.vertex(0).numPosition();
    size_t Li = variableIndex(0,0);
    for (size_t i=0; i<numCells; ++i) {
      size_t numWalls = T.cell(i).numWall();
      for (size_t k=0; k<numWalls; ++k) {
	size_t kp = k<numWalls-1 ? k+1 : 0;
	size_t km = k>0 ? k-1 : numWalls-1;

	// Get global vertex indices
	size_t j = T.cell(i).vertex(k)->index();
	size_t jp = T.cell(i).vertex(kp)->index();
	size_t jm = T.cell(i).vertex(km)->index();
	
	// Get edges
	size_t ep = T.cell(i).wall(k)->index();
	size_t em = T.cell(i).wall(km)->index();
	
	// Update for each dimension
	for (size_t d=0; d<dimension; ++d) {
	  double derivs = -parameter(0)*(vertexData[j][d] - 
					 (vertexData[jp][d]*wallData[em][Li] + 
					  vertexData[jm][d]*wallData[ep][Li])/
					 (wallData[em][Li]+wallData[ep][Li]) );
	  vertexDerivs[k][d] += derivs; 
	}
      }      
    }
  }
}
