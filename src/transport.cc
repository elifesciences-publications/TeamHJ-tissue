//
// Filename     : transport.cc
// Description  : Classes describing transport reactions
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : September 2013
// Revision     : $Id:$
//

#include"transport.h"

MembraneDiffusionSimple::
MembraneDiffusionSimple(std::vector<double> &paraValue, 
		  std::vector< std::vector<size_t> > 
		  &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=1 ) {
    std::cerr << "MembraneDiffusionSimple::"
	      << "MembraneDiffusionSimple() "
	      << "One parameters (Diffusion constant) used." << std::endl;
    exit(EXIT_FAILURE);
  }
  if( indValue.size() != 1 || indValue[0].size() != 1 ) {
    std::cerr << "MembraneDiffusionSimple::"
	      << "MembraneDiffusionSimple() "
	      << "One variable (pair) used for diffusing wall variable." << std::endl;
    exit(EXIT_FAILURE);
  }
  //Set the variable values
  //
  setId("MembraneDiffusionSimple");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "D_mem";
  
  setParameterId( tmp );
}

void MembraneDiffusionSimple::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{  
  
  size_t numCells = T.numCell();
  size_t pwI = variableIndex(0,0);//diffusing molecule (membrane/wall)
  
  assert(pwI<wallData[0].size());

  for (size_t i=0; i<numCells; ++i) {
    size_t numWalls = T.cell(i).numWall();
    for (size_t k=0; k<numWalls; ++k) {
      size_t j = T.cell(i).wall(k)->index();
      
      size_t pwIadd=1;//Two variables per wall, setting which to use
      if( T.cell(i).wall(k)->cell1()->index() == i ) {
	pwIadd=0;
      }
      double fac = parameter(3)*wallData[j][pwI+pwIadd];
      wallDerivs[j][pwI+pwIadd] -= fac;

      size_t kNext = (k+1)%numWalls;     
      size_t jNext = T.cell(i).wall(kNext)->index();   
      if( T.cell(i).wall(kNext)->cell1()->index() == i) {
	wallDerivs[jNext][pwI] +=fac;
      }
      else {
	wallDerivs[jNext][pwI+1] +=fac;
      }	
      size_t kBefore = k > 0 ? k-1 : numWalls-1; 
      size_t jBefore = T.cell(i).wall(kBefore)->index(); 	
      if( T.cell(i).wall(kBefore)->cell1()->index() == i) {
	wallDerivs[jBefore][pwI] +=fac;
      }
      else {
	wallDerivs[jBefore][pwI + 1] +=fac;
      }	
    }
  }
}
