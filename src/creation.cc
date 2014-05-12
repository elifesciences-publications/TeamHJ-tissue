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

CreationSpatialSphere::
CreationSpatialSphere(std::vector<double> &paraValue, 
	     std::vector< std::vector<size_t> > 
	     &indValue ) 
{  

  // Do some checks on the parameters and variable indeces
  if( paraValue.size()!=4 ) {
    std::cerr << "CreationSpatialSphere::CreationSpatialSphere() "
	      << "Uses four parameters V_max R(K_Hill) n_Hill and R_sign\n";
    exit(0);
  }
  if( indValue.size() != 1 || indValue[0].size() != 1 ) {
    std::cerr << "CreationSpatialSphere::"
	      << "CreationSpatialSphere() "
	      << "Index for variable to be updated given." << std::endl;
    exit(0);
  }
  // Sign should be -/+1
  if( paraValue[3] != -1 && paraValue[3] != 1 ) {
    std::cerr << "CreationSpatialSphere::CreationSpatialSphere() "
	      << "R_sign should be +/-1 to set production inside/outside R"
	      << std::endl;
    exit(0);
  }
	
  // Set the variable values
  setId("creationSpatialSphere");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  // Set the parameter identities
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "V_max";
  tmp[1] = "R (K_Hill)";
  tmp[2] = "n_Hill";
  tmp[3] = "R_sign";  
  setParameterId( tmp );
}

void CreationSpatialSphere::
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
  double powK_ = std::pow(parameter(1),parameter(2));
  //For each cell
  for (size_t cellI = 0; cellI < numCells; ++cellI) {

    //Calculate cell center from vertices positions
    std::vector<double> cellCenter;
    cellCenter = T.cell(cellI).positionFromVertex(vertexData);
    assert( cellCenter.size() == vertexData[0].size() );
    double r=0.0;
    for( size_t d=0 ; d<cellCenter.size() ; ++d )
      r += cellCenter[d]*cellCenter[d];
    r = std::sqrt(r);

    double powR = std::pow(r,parameter(2));
	
    if (parameter(3)>0.0)
      cellDerivs[cellI][cIndex] += parameter(0)*powR/(powK_+powR);
    else
      cellDerivs[cellI][cIndex] += parameter(0)*powK_/(powK_+powR);


  }
}
