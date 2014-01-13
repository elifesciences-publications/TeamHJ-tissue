//
// Filename     : creation.cc
// Description  : Classes describing molecular production/creation updates
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : January 2011
// Revision     : $Id:$
//
#include"tissue.h"
#include"baseReaction.h"
#include"creation.h"
#include<cmath>

CreationZero::
CreationZero(std::vector<double> &paraValue, 
	     std::vector< std::vector<size_t> > 
	     &indValue ) 
{  
  // Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=1 ) {
    std::cerr << "CreationZero::"
	      << "CreationZero() "
	      << "Uses one parameter k_c (constant production rate)." << std::endl;
    exit(0);
  }
  if( indValue.size() != 1 || indValue[0].size() != 1 ) {
    std::cerr << "CreationZero::"
	      << "CreationZero() "
	      << "Index for variable to be updated given." << std::endl;
    exit(0);
  }
  //Set the variable values
  //
  setId("CreationZero");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "k_c";
  setParameterId( tmp );
}

void CreationZero::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{  
  //Do the update for each cell
  size_t numCells = T.numCell();

  size_t cIndex = variableIndex(0,0);
  double k_c = parameter(0);
  //For each cell
  for (size_t cellI = 0; cellI < numCells; ++cellI) {      
    cellDerivs[cellI][cIndex] += k_c;
  }
}

CreationOne::
CreationOne(std::vector<double> &paraValue, 
	     std::vector< std::vector<size_t> > 
	     &indValue ) 
{  
  // Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=1 ) {
    std::cerr << "CreationOne::"
	      << "CreationOne() "
	      << "Uses one parameter k_c (linear production rate)." << std::endl;
    exit(0);
  }
  if( indValue.size() != 2 || indValue[0].size() != 1 || indValue[1].size() != 1 ) {
    std::cerr << "CreationOne::"
	      << "CreationOne() "
	      << "Index for variable to be updated given in first row and "
	      << "index for production-dependent variable in 2nd." << std::endl;
    exit(0);
  }
  //Set the variable values
  //
  setId("CreationOne");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "k_c";
  setParameterId( tmp );
}

void CreationOne::
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
  double k_c = parameter(0);
  //For each cell
  for (size_t cellI = 0; cellI < numCells; ++cellI) {      
    cellDerivs[cellI][cIndex] += k_c * cellData[cellI][xIndex];
  }
}





