//
// Filename     : degradation.cc
// Description  : Classes describing molecular degradation updates
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : January 2011
// Revision     : $Id:$
//
#include"tissue.h"
#include"baseReaction.h"
#include"degradation.h"
#include<cmath>

DegradationOne::
DegradationOne(std::vector<double> &paraValue, 
	     std::vector< std::vector<size_t> > 
	     &indValue ) 
{  
  // Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=1 ) {
    std::cerr << "DegradationOne::"
	      << "DegradationOne() "
	      << "Uses one parameter k_d (constant degradation rate)." << std::endl;
    exit(0);
  }
  if( indValue.size() != 1 || indValue[0].size() != 1 ) {
    std::cerr << "DegradationOne::"
	      << "DegradationOne() "
	      << "Index for variable to be updated (degraded) given." << std::endl;
    exit(0);
  }
  //Set the variable values
  //
  setId("DegradationOne");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "k_d";
  setParameterId( tmp );
}

void DegradationOne::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) {
  
  //Do the update for each cell
  size_t numCells = T.numCell();
  
  size_t cIndex = variableIndex(0,0);
  double k_d = parameter(0);
  //For each cell
  for (size_t cellI = 0; cellI < numCells; ++cellI) {      
    cellDerivs[cellI][cIndex] -= k_d * cellData[cellI][cIndex];
  }
}

DegradationTwo::
DegradationTwo(std::vector<double> &paraValue, 
	       std::vector< std::vector<size_t> > 
	       &indValue ) 
{  
  // Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=1 ) {
    std::cerr << "DegradationTwo::"
	      << "DegradationTwo() "
	      << "Uses one parameter k_d (degradation rate)." << std::endl;
    exit(0);
  }
  if( indValue.size() != 2 || indValue[0].size() != 1 || indValue[1].size() != 1 ) {
    std::cerr << "DegradationTwo::"
	      << "DegradationTwo() "
	      << "Index for variable to be updated given in first row and "
	      << "index for degradation-dependent variable in 2nd." << std::endl;
    exit(0);
  }
  //Set the variable values
  //
  setId("DegradationTwo");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "k_d";
  setParameterId( tmp );
}

void DegradationTwo::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) {
  
  //Do the update for each cell
  size_t numCells = T.numCell();

  size_t cIndex = variableIndex(0,0);
  size_t xIndex = variableIndex(1,0);
  double k_d = parameter(0);
  //For each cell
  for (size_t cellI = 0; cellI < numCells; ++cellI) {      
    cellDerivs[cellI][cIndex] -= k_d * cellData[cellI][xIndex] * cellData[cellI][cIndex]; 
  }
}

