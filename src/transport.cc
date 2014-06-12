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
      double fac = parameter(0)*wallData[j][pwI+pwIadd];
      wallDerivs[j][pwI+pwIadd] -= 2.*fac;

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



MembraneDiffusionSimple2::
MembraneDiffusionSimple2(std::vector<double> &paraValue, 
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

void MembraneDiffusionSimple2::
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

  for (size_t i=5; i<6; ++i) {
    size_t numWalls = T.cell(i).numWall();
    for (size_t k=0; k<numWalls; ++k) {
      size_t j = T.cell(i).wall(k)->index();
      
      size_t pwIadd=1;//Two variables per wall, setting which to use
      if( T.cell(i).wall(k)->cell1()->index() == i ) {
	pwIadd=0;
      }
      double fac = parameter(0)*wallData[j][pwI+pwIadd];
      wallDerivs[j][pwI+pwIadd] -= 2.*fac;

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

for (size_t i=9; i<10; ++i) {
    size_t numWalls = T.cell(i).numWall();
    for (size_t k=0; k<numWalls; ++k) {
      size_t j = T.cell(i).wall(k)->index();
      
      size_t pwIadd=1;//Two variables per wall, setting which to use
      if( T.cell(i).wall(k)->cell1()->index() == i ) {
	pwIadd=0;
      }
      double fac = parameter(0)*wallData[j][pwI+pwIadd];
      wallDerivs[j][pwI+pwIadd] -= 2.*fac;

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




 ActiveTransportCellEfflux::
 ActiveTransportCellEfflux(std::vector<double> &paraValue, 
 	      std::vector< std::vector<size_t> > 
 	      &indValue ) {
  
   //Do some checks on the parameters and variable indeces
   //
   if( paraValue.size()!=1 ) {
     std::cerr << "ActiveTransportCellEfflux::"
 	      << "ActiveTransportCellEfflux() "
 	      << "1 parameters used (see transport.h)\n";
     exit(0);
   }
   if( indValue.size() != 2 || indValue[0].size() != 1 || indValue[1].size() != 1 ) {
     std::cerr << "ActiveTransportCellEfflux::"
 	      << "ActiveTransportCellEfflux() "
 	      << "One cell variable indices (auxin) and one wall variable"
 	      << " indices are used (PIN)." << std::endl;
    exit(0);
  }
  //Set the variable values
  //
  setId("ActiveTransportCellEfflux");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "T";
   setParameterId( tmp );
}

void ActiveTransportCellEfflux::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{  
  size_t numCells = T.numCell();
  size_t aI = variableIndex(0,0);//auxin
  size_t pwI = variableIndex(1,0);//pin (membrane/wall)


  assert( aI<cellData[0].size() &&
	  pwI<wallData[0].size() );

 for (size_t i=0; i<numCells; ++i) {
    
    //Auxin transport and protein cycling
    size_t numWalls = T.cell(i).numWall();
    for (size_t k=0; k<numWalls; ++k) {
      size_t j = T.cell(i).wall(k)->index();
      if( T.cell(i).wall(k)->cell1()->index() == i && T.cell(i).wall(k)->cell2() != T.background() ) {
	// cell-cell transport
	size_t iNeighbor = T.cell(i).wall(k)->cell2()->index();
	double fac = parameter(0)*cellData[i][aI]*wallData[j][pwI]; //p_2 A_i + p_4 A_i P_ij
	
	cellDerivs[i][aI] -= fac;
	cellDerivs[iNeighbor][aI] += fac;
      }
      else if( T.cell(i).wall(k)->cell2()->index() == i && T.cell(i).wall(k)->cell1() != T.background() ) {
	// cell-cell transport
	size_t iNeighbor = T.cell(i).wall(k)->cell1()->index();
	double fac =  parameter(0)*cellData[i][aI]*wallData[j][pwI+1];
	
	cellDerivs[i][aI] -= fac;
	cellDerivs[iNeighbor][aI] += fac;
      }
    }
  }
}
