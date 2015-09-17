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


void DegradationOne::
derivsWithAbs(Tissue &T,
        DataMatrix &cellData,
        DataMatrix &wallData,
        DataMatrix &vertexData,
        DataMatrix &cellDerivs,
        DataMatrix &wallDerivs,
        DataMatrix &vertexDerivs,
        DataMatrix &sdydtCell,
        DataMatrix &sdydtWall,
        DataMatrix &sdydtVertex ) 
{
  //Do the update for each cell
  size_t numCells = T.numCell();
  
  size_t cIndex = variableIndex(0,0);
  size_t xIndex = variableIndex(1,0);
  double k_c = parameter(0);
  //For each cell
  for (size_t cellI = 0; cellI < numCells; ++cellI) {  
      double value = k_c * cellData[cellI][xIndex];
    
    cellDerivs[cellI][cIndex]  -= value;
    sdydtCell[cellI][cIndex]  += value;
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

DegradationHill::DegradationHill(std::vector<double> &paraValue, 
				 std::vector< std::vector<size_t> > 
				 &indValue ) 
{  
  // Do some checks on the parameters and variable indeces
  if( paraValue.size()!=3 ) {
    std::cerr << "DegradationHill::DegradationHill() "
	      << "Uses three parameters k_deg K_hill n_hill\n";
    exit(0);
  }
  if( indValue.size() != 2 || indValue[0].size() != 1 || indValue[1].size() != 1 ) {
    std::cerr << "DegradationHill::DegradationHill() "
	      << "One variable index in first layer for variable to be updated "
	      << "and one in the second that activate degradation must be given." 
	      << std::endl;
    exit(0);
  }
  
  // Set the variable values
  setId("DegradationHill");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  // Set the parameter identities
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "k_deg";
  tmp[1] = "K_hill";
  tmp[2] = "n_hill";
  setParameterId( tmp );
}

void DegradationHill::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs )
{  
  size_t degradedIndex = variableIndex(0,0);
  size_t degraderIndex = variableIndex(1,0);
  double Kpow = std::pow(parameter(1),parameter(2));

  for (size_t cellI=0; cellI<T.numCell(); ++cellI) {
    double fac = parameter(0)*cellData[cellI][degradedIndex];
    double yPow = std::pow(cellData[cellI][degraderIndex],parameter(2));
    cellDerivs[cellI][degradedIndex] -= fac*yPow/(Kpow+yPow);
  }
}



DegradationTwoGeometric::
DegradationTwoGeometric(std::vector<double> &paraValue, 
	       std::vector< std::vector<size_t> > 
	       &indValue ) 
{  
  // Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=1 ) {
    std::cerr << "DegradationTwoGeometric::"
	      << "DegradationTwoGeometric() "
	      << "Uses one parameter k_d (degradation rate)." << std::endl;
    exit(0);
  }
  if( indValue.size() != 2 || indValue[0].size() != 1 || indValue[1].size() != 1 ) {
    std::cerr << "DegradationTwoGeometric::"
	      << "DegradationTwoGeometric() "
	      << "Index for variable to be updated given in first row and "
	      << "index for degradation-dependent variable in 2nd." << std::endl;
    exit(0);
  }
  //Set the variable values
  //
  setId("DegradationTwoGeometric");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "k_d";
  setParameterId( tmp );
}

void DegradationTwoGeometric::
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
   double cellVolume = T.cell(cellI).calculateVolume(vertexData);
    cellDerivs[cellI][cIndex] -= cellVolume* k_d * cellData[cellI][xIndex] * cellData[cellI][cIndex]; 
  }
}
