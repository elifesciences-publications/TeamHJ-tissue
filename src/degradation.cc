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
  double k_c = parameter(0);
  //For each cell
  for (size_t cellI = 0; cellI < numCells; ++cellI) {  
      double value = k_c * cellData[cellI][cIndex];
    
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

DegradationN::
DegradationN(std::vector<double> &paraValue, 
	       std::vector< std::vector<size_t> > 
	       &indValue ) 
{  
  // Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=1 ) {
    std::cerr << "DegradationN::"
	      << "DegradationN() "
	      << "Uses one parameter k_d (degradation rate)." << std::endl;
    exit(0);
  }
  if( indValue.size() != 2 || indValue[0].size() != 1) {
    std::cerr << "DegradationN::"
	      << "DegradationN() "
	      << "Index for variable to be updated given in first row and "
	      << "indices for degradation-dependent variables in 2nd." << std::endl;
    exit(0);
  }
  //Set the variable values
  //
  setId("DegradationN");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "k_d";
  setParameterId( tmp );
}

void DegradationN::
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

    double product = k_d;
    for( size_t i=0 ; i<numVariableIndex(1) ; i++ ) { // for each rate-limiting variable
      product *= cellData[cellI][variableIndex(1,i)]; // multiply rate with its concentration
    }

    cellDerivs[cellI][cIndex] -= product*cellData[cellI][cIndex]; 
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

DegradationHillN::DegradationHillN(std::vector<double> &paraValue, 
				 std::vector< std::vector<size_t> > 
				 &indValue ) 
{  
  // Do some checks on the parameters and variable indeces


  if( indValue.size() != 3 || indValue[0].size() != 1 ) {
    std::cerr << "DegradationHillN::DegradationHillN() "
	      << "Three levels of variable indices are used, "
	      << "the first for the updated molecule, 2nd for activators and 3rd for repressors"
	      << std::endl;
    exit(EXIT_FAILURE);
  }
  if( 1 + 2*(indValue[1].size()+indValue[2].size()) != paraValue.size() ) {
    std::cerr << "Hill::Hill() "
	      << "Number of parameters does not agree with number of "
	      << "activators/repressors.\n"
	      << indValue[1].size() << " activators, " 
	      << indValue[2].size() << " repressors, and "
	      << paraValue.size() << " parameters\n"
 	      << "One plus pairs of parameters (V_max,K_half1,n_Hill1,"
 	      << "K2,n2,...) must be given.\n";
    exit(EXIT_FAILURE);
  }
  
  // Set the variable values
  setId("DegradationHill");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  // Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "d_max";
  size_t count=1;
  for (size_t i=0; i<indValue[1].size(); ++i) {
    tmp[count++] = "K_A";
    tmp[count++] = "n_A";
  }
  for (size_t i=0; i<indValue[2].size(); ++i) {
    tmp[count++] = "K_R";
    tmp[count++] = "n_R";
  }
  setParameterId( tmp );
}

void DegradationHillN::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs )
{  

  // Do the update for each cell
  size_t numCells = T.numCell();
  size_t cIndex = variableIndex(0,0);
  //For each cell
  for (size_t cellI = 0; cellI < numCells; ++cellI) {      
  //
    double contribution=parameter(0);
    size_t parameterIndex=1;
    // Activator contributions
    for( size_t i=0 ; i<numVariableIndex(1) ; i++ ) {
      double c = std::pow(cellData[cellI][variableIndex(1,i)], 
      	   parameter(parameterIndex+1));
      contribution *= c
      / ( std::pow(parameter(parameterIndex),parameter(parameterIndex+1)) + c );
      parameterIndex+=2;
    }
    // Repressor contributions
   for( size_t i=0 ; i<numVariableIndex(2) ; i++ ) {
     double c = std::pow(parameter(parameterIndex),parameter(parameterIndex+1));
     contribution *= c /
       ( c + std::pow(cellData[cellI][variableIndex(2,i)],
		      parameter(parameterIndex+1)) );
     parameterIndex+=2;
   }   
   cellDerivs[cellI][cIndex] -= contribution*cellData[cellI][cIndex]; 
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


DegradationOneWall::
DegradationOneWall(std::vector<double> &paraValue, 
	     std::vector< std::vector<size_t> > 
	     &indValue ) 
{  
  // Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=1 ) {
    std::cerr << "DegradationOneWall::"
	      << "DegradationOneWall() "
	      << "Uses one parameter k_cw (constant degradation rate)." << std::endl;
    exit(0);
  }
  if( indValue.size() != 1 || indValue[0].size() != 1 ) {
    std::cerr << "DegradationOneWall::"
	      << "DegradationOneWall() "
	      << "Index for wall variable to be updated (degraded) given." << std::endl;
    exit(0);
  }
  //Set the variable values
  //
  setId("DegradationOneWall");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "k_cw";
  setParameterId( tmp );
}

void DegradationOneWall::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) {


  size_t numWalls = T.numWall();
  size_t wIndex = variableIndex(0,0);
  double k_cw = parameter(0);
  for (size_t k=0; k<numWalls; ++k) {
    wallDerivs[k][wIndex] -= k_cw * wallData[k][wIndex];
    //std::cerr << k << "\t" << wIndex << "\t" << wallDerivs[k][wIndex] << "\n";
  }



}


void DegradationOneWall::
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


  size_t numWalls = T.numWall();
  size_t wIndex = variableIndex(0,0);
  double k_cw = parameter(0);
  for (size_t k=0; k<numWalls; ++k) {
    double value = k_cw * wallData[k][wIndex];
    wallDerivs[k][wIndex] -= value;
    sdydtWall[k][wIndex] += value;
  }

}
