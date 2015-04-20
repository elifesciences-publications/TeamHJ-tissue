//
// Filename     : network.cc
// Description  : Classes describing complete models
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : November 2006
// Revision     : $Id:$
//
#include "network.h"
#include "baseReaction.h"

AuxinModelSimple1::
AuxinModelSimple1(std::vector<double> &paraValue, 
		  std::vector< std::vector<size_t> > 
		  &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=12 ) {
    std::cerr << "AuxinModelSimple1::"
	      << "AuxinModelSimple1() "
	      << "Twelve parameters used.\n\n";
    std::cerr << "dA_i/dt = p0*M_i + p1 - p2*A_i +p5*Sum_{neigh} (A_n-A_i) +\n" 
	      << "p4*Sum_{neigh} (P_ni*A_n-P_in*A_i)\n\n" 
	      << "dP_i/dt = p6 - p7*P_i\n\n"
	      << "dX_i/dt = p8*A_i - p9*X_i\n\n"
	      << "dM_i/dt = p10*Theta_L1 - p11*M_i\n\n"
	      << "P_in = P_i*X_n/(p_3+Sum_{k,neigh}X_k)\n";
    exit(EXIT_FAILURE);
  }
  if( (indValue.size() != 1 && indValue.size() !=2) || indValue[0].size() != 4 ||
      (indValue.size()==2 &&indValue[1].size() != 1) ) {
    std::cerr << "AuxinModelSimple1::"
	      << "AuxinModelSimple1() "
	      << "Four variable indices are used (auxin,pin,X,M) in first level."
	      << std::endl << "Optionally a wall index is given for saving (pair of) PIN"
	      << " at membranes." << std::endl;
    exit(EXIT_FAILURE);
  }
  //Set the variable values
  //
  setId("AuxinModelSimple1");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "p_auxin(M)";
  tmp[1] = "p_auxin";
  tmp[2] = "d_auxin";
  tmp[3] = "p_pol";
  tmp[4] = "T_auxin";
  tmp[5] = "D_auxin";
  tmp[6] = "p_pin";
  tmp[7] = "d_pin";
  tmp[8] = "p_X";
  tmp[9] = "d_X";
  tmp[10] = "p_M";
  tmp[11] = "d_M";
  
  setParameterId( tmp );
}

void AuxinModelSimple1::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{  
  size_t numCells = T.numCell();
  size_t aI = variableIndex(0,0);
  size_t pI = variableIndex(0,1);
  size_t xI = variableIndex(0,2);
  size_t mI = variableIndex(0,3);
  assert( aI<cellData[0].size() &&
	  pI<cellData[0].size() &&
	  xI<cellData[0].size() &&
	  mI<cellData[0].size() );
  
  if (numVariableIndexLevel()==2) {
    assert(variableIndex(1,0)<wallData[0].size());
  }

  for( size_t i=0 ; i<numCells ; ++i ) {
    
    //Production and degradation
    cellDerivs[i][aI] += parameter(0)*cellData[i][mI] + parameter(1) - 
      parameter(2)*cellData[i][aI];
    
    cellDerivs[i][pI] += parameter(6) - parameter(7)*cellData[i][pI];
    
    cellDerivs[i][xI] += parameter(8)*cellData[i][aI] 
      - parameter(9)*cellData[i][xI];
    
    cellDerivs[i][mI] -= parameter(11)*cellData[i][mI];
    if( T.cell(i).isNeighbor(T.background()) )
      cellDerivs[i][mI] += parameter(10);
    
    //Transport
    size_t numWalls=T.cell(i).numWall();
    //Polarization coefficient normalization constant
    double sum=0.0;
    size_t numActualWalls=0;
    //pin[i].resize( numWalls+1 );
    for( size_t n=0 ; n<numWalls ; ++n ) {
      if( T.cell(i).wall(n)->cell1() != T.background() &&
	  T.cell(i).wall(n)->cell2() != T.background() ) { 
	numActualWalls++;
	if( T.cell(i).wall(n)->cell1()->index()==i )
	  sum += cellData[ T.cell(i).wall(n)->cell2()->index() ][ xI ];
	else
	  sum += cellData[ T.cell(i).wall(n)->cell1()->index() ][ xI ];
      }
    }
    //sum /= numActualWalls;//For adjusting for different num neigh
    sum += parameter(3);
    
    for( size_t n=0 ; n<numWalls ; ++n ) {
      //if( !T.cell(i).isNeighbor(T.background()) ) { 
      if( T.cell(i).wall(n)->cell1() != T.background() &&
	  T.cell(i).wall(n)->cell2() != T.background() ) { 
	size_t neighIndex;
	size_t pinIndexAdd=0; 
	if( T.cell(i).wall(n)->cell1()->index()==i )
	  neighIndex = T.cell(i).wall(n)->cell2()->index();				
	else {
	  neighIndex = T.cell(i).wall(n)->cell1()->index();				
	  pinIndexAdd=1;
	}
	double polRate=0.0;
	
	if( sum != 0.0 ) {
	  polRate = cellData[i][pI] * cellData[neighIndex][xI] / sum;
	}
	else {
	  polRate = 1.;
	}
	if (numVariableIndexLevel()==2) {//store PIN value for membrane (in wall pair data)
	  wallData[T.cell(i).wall(n)->index()][variableIndex(1,0)+pinIndexAdd] = polRate;
	}
	//pin[i][n+1] = polRate;
	cellDerivs[i][aI] -= (parameter(4)*polRate+parameter(5))*cellData[i][aI];
	cellDerivs[neighIndex][aI] += (parameter(4)*polRate+parameter(5))*cellData[i][aI];
      }
    }
  }
}

AuxinModelStress::
AuxinModelStress(std::vector<double> &paraValue, 
		       std::vector< std::vector<size_t> > 
		       &indValue ) 
{  
  //
  // Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=9 && paraValue.size()!=11 ) {
    std::cerr << "AuxinModelStress::"
	      << "AuxinModelStress() "
	      << "Seven plus three optional parameters used.\n\n";
    std::cerr << "dA_i/dt = p0 - p1*A_i +p2*Sum_{neigh} (A_n-A_i) +\n" 
	      << "p3*Sum_{neigh} (P_ni*A_n-P_in*A_i)\n\n" 
	      << "dP_i/dt = p4 - p5*P_i\n\n"
	      << "P_in = P_i*X_in/(p6+Sum_k X_ik)\n\n"
	      << "S_in = k_in*F/(k_in+k_ni)\n\n"
	      << "X_in = S_in^p8/(S_in^p8+p7^p8)\n\n"
	      << "k_in = p9 + p10/(p11+A_i^p12)  (optional if updated)\n";
    exit(0);
  }
  if( indValue.size() != 2 || indValue[0].size() != 2 || indValue[1].size() != 3 ) {
    std::cerr << "AuxinModelStress::"
	      << "AuxinModelStress() "
	      << "Two cell variable indices are used in first level (auxin,pin),\n"
	      << "and three wall indices are used in second level (F,k_1,k_2),\n";
    exit(0);
  }
  //Set the variable values
  //
  setId("AuxinModelStress");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "p_auxin";
  tmp[1] = "d_auxin";
  tmp[2] = "D_auxin";
  tmp[3] = "T_auxin";
  tmp[4] = "p_pin";
  tmp[5] = "d_pin";
  tmp[6] = "K_pin";
  if( numParameter()==10 ) {
    tmp[7] = "k0_k";
    tmp[8] = "k1_k";
    tmp[9] = "K_k";
    tmp[10] = "n_k";
    tmp[11] = "K_stress";
    tmp[12] = "n_stress";
  }
  setParameterId( tmp );
}

void AuxinModelStress::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{  
  size_t numCells = T.numCell();
  size_t aI = variableIndex(0,0);
  size_t pI = variableIndex(0,1);
  size_t FI = variableIndex(1,0);
  std::vector<size_t> kI(2);
  kI[0] = variableIndex(1,1);
  kI[1] = variableIndex(1,2);
  assert( aI<cellData[0].size() &&
	  pI<cellData[0].size() &&
	  FI<wallData[0].size() &&
	  kI[0]<wallData[0].size() &&
	  kI[1]<cellData[0].size() );
  size_t dimension = T.vertex(0).numPosition();
  std::vector<double> pos(dimension);
  
  for( size_t i=0 ; i<numCells ; ++i ) {
    
    //Production and degradation
    pos = T.cell(i).positionFromVertex(vertexData);
    double r=0;
    for(size_t d=0; d<dimension; ++d)
      r += pos[d]*pos[d];
    r = std::sqrt(r);
    if (r>3.0)
      cellDerivs[i][aI] += parameter(0);
    
    cellDerivs[i][aI] -= parameter(1)*cellData[i][aI];    
    cellDerivs[i][pI] += parameter(4) - parameter(5)*cellData[i][pI];
    
    //Transport
    size_t numWalls=T.cell(i).numWall();
    //Polarization coefficient normalization constant
    double sum=0.0;
    std::vector<double> pin(numWalls);
    for( size_t n=0 ; n<numWalls ; ++n ) {
      size_t kII=1;
      size_t k = T.cell(i).wall(n)->index();
      if( T.cell(i).wall(n)->cell1()->index() == i ) {
	kII=0;
      }
      double Sn = wallData[k][FI]*wallData[k][kI[kII]]/(wallData[k][kI[0]]+wallData[k][kI[1]]);
      if (Sn>=0.0) {
	//sum += pin[n] = ( Sn >= 0.0 ? Sn : 0.0 );
	Sn = std::pow(Sn,parameter(12));
	Sn = Sn/(Sn+std::pow(parameter(11),parameter(12)));
	sum += pin[n] = Sn;
      }
      else
	sum += pin[n] = 0.0;
    }
    //sum /= numWalls;//For adjusting for different num neigh
    sum += parameter(6);
    // Actual passive and active transport
    for( size_t n=0 ; n<numWalls ; ++n ) {
      if( T.cell(i).wall(n)->cell1() != T.background() &&
	  T.cell(i).wall(n)->cell2() != T.background() ) { 
	size_t neighIndex; 
	if( T.cell(i).wall(n)->cell1()->index()==i )
	  neighIndex = T.cell(i).wall(n)->cell2()->index();				
	else
	  neighIndex = T.cell(i).wall(n)->cell1()->index();				
	
	double polRate=0.0;
	if( sum != 0.0 )
	  polRate = cellData[i][pI] * pin[n] / sum;
	
	cellDerivs[i][aI] -= (parameter(3)*polRate+parameter(2))*cellData[i][aI];
	cellDerivs[neighIndex][aI] += (parameter(3)*polRate+parameter(2))*cellData[i][aI];
      }
    }
    // Calculate the new k values from the auxin concentration
    if( numParameter() == 10 ) {
      for( size_t n=0 ; n<numWalls ; ++n ) {
	size_t kII=1;
	size_t k = T.cell(i).wall(n)->index();
	if( T.cell(i).wall(n)->cell1()->index() == i ) {
	  kII=0;
				}
	//wallData[k][kI[kII]] = parameter(7) + parameter(8)/(parameter(9)+cellData[i][aI]);
	wallData[k][kI[kII]] = parameter(7) + parameter(8)/(parameter(9)+std::pow(cellData[i][aI],parameter(10)));
      }    
    }
  }
}

AuxinModelSimpleStress::
AuxinModelSimpleStress(std::vector<double> &paraValue, 
		       std::vector< std::vector<size_t> > 
		       &indValue ) 
{  
  //
  // Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=7 && paraValue.size()!=10 ) {
    std::cerr << "AuxinModelSimpleStress::"
	      << "AuxinModelSimpleStress() "
	      << "Seven plus three optional parameters used.\n\n";
    std::cerr << "dA_i/dt = p0 - p1*A_i +p2*Sum_{neigh} (A_n-A_i) +\n" 
	      << "p3*Sum_{neigh} (P_ni*A_n-P_in*A_i)\n\n" 
	      << "dP_i/dt = p4 - p5*P_i\n\n"
	      << "P_in = P_i*X_in/(p6+Sum_k X_ik)\n\n"
	      << "X_in = k_in*F/(k_in+k_ni)\n\n"
	      << "k_in = p7 + p8/(p9+A_i)  (optional if updated)\n";
    exit(0);
  }
  if( (indValue.size() != 2 && indValue.size() != 3) || 
      indValue[0].size() != 2 || indValue[1].size() != 3 || 
      (indValue.size() == 3 && indValue[2].size() != 1) ) {
    std::cerr << "AuxinModelSimpleStress::"
	      << "AuxinModelSimpleStress() "
	      << "Two cell variable indices are used in first level (auxin,pin),\n"
	      << " three wall indices are used in second level (F,k_1,k_2),\n"
	      << " and one index may be used in third layer for wall pin storage\n";
    exit(0);
  }
  //Set the variable values
  //
  setId("AuxinModelSimpleStress");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "p_auxin";
  tmp[1] = "d_auxin";
  tmp[2] = "D_auxin";
  tmp[3] = "T_auxin";
  tmp[4] = "p_pin";
  tmp[5] = "d_pin";
  tmp[6] = "K_pin";
  if( numParameter()==10 ) {
    tmp[7] = "k0_k";
    tmp[8] = "k1_k";
    tmp[9] = "K_k";
  }
  
  setParameterId( tmp );
}

void AuxinModelSimpleStress::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{  
  size_t numCells = T.numCell();
  size_t aI = variableIndex(0,0);
  size_t pI = variableIndex(0,1);
  size_t FI = variableIndex(1,0);
  std::vector<size_t> kI(2);
  kI[0] = variableIndex(1,1);
  kI[1] = variableIndex(1,2);
  assert( aI<cellData[0].size() &&
	  pI<cellData[0].size() &&
	  FI<wallData[0].size() &&
	  kI[0]<wallData[0].size() &&
	  kI[1]<cellData[0].size() );
  
  for( size_t i=0 ; i<numCells ; ++i ) {
    
    //Production and degradation
		
    cellDerivs[i][aI] += parameter(0) - parameter(1)*cellData[i][aI];
    
    cellDerivs[i][pI] += parameter(4) - parameter(5)*cellData[i][pI];
    
    //Transport
    size_t numWalls=T.cell(i).numWall();
    //Polarization coefficient normalization constant
    double sum=0.0;
    std::vector<double> pin(numWalls);
    for( size_t n=0 ; n<numWalls ; ++n ) {
      size_t kII=1;
      size_t k = T.cell(i).wall(n)->index();
      if( T.cell(i).wall(n)->cell1()->index() == i ) {
	kII=0;
      }
      double Sn = wallData[k][FI]*wallData[k][kI[kII]]/(wallData[k][kI[0]]+wallData[k][kI[1]]);
      //sum += pin[n] = ( Sn >= 0.0 ? Sn : 0.0 );
      sum += pin[n] = ( Sn>=0.0 ? std::pow(Sn,3) : 0.0 );
    }
    //sum /= numWalls;//For adjusting for different num neigh
    sum += parameter(6);
    // Actual passive and active transport
    for( size_t n=0 ; n<numWalls ; ++n ) {
      if( T.cell(i).wall(n)->cell1() != T.background() &&
	  T.cell(i).wall(n)->cell2() != T.background() ) { 
	size_t neighIndex; 
	if( T.cell(i).wall(n)->cell1()->index()==i )
	  neighIndex = T.cell(i).wall(n)->cell2()->index();				
	else
	  neighIndex = T.cell(i).wall(n)->cell1()->index();				
	
	double polRate=0.0;
	if( sum != 0.0 )
	  polRate = cellData[i][pI] * pin[n] / sum;

	// Store PIN1
	if (numVariableIndexLevel()==3) {
	  size_t pinIndex = variableIndex(2,0);
	  if (neighIndex<i) 
	    pinIndex++;
	  wallData[T.cell(i).wall(n)->index()][pinIndex]=polRate;
	}
	
	cellDerivs[i][aI] -= (parameter(3)*polRate+parameter(2))*cellData[i][aI];
	cellDerivs[neighIndex][aI] += (parameter(3)*polRate+parameter(2))*cellData[i][aI];
      }
    }
    // Calculate the new k values from the auxin concentration
    if( numParameter() == 10 ) {
      for( size_t n=0 ; n<numWalls ; ++n ) {
	size_t kII=1;
	size_t k = T.cell(i).wall(n)->index();
	if( T.cell(i).wall(n)->cell1()->index() == i ) {
	  kII=0;
	}
	//wallData[k][kI[kII]] = parameter(7) + parameter(8)/(parameter(9)+cellData[i][aI]);
	wallData[k][kI[kII]] = parameter(7) + parameter(8)/(parameter(9)+std::pow(cellData[i][aI],3));
      }    
    }
  }
}

AuxinModelSimple1Wall::
AuxinModelSimple1Wall(std::vector<double> &paraValue, 
		      std::vector< std::vector<size_t> > 
		      &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=7 ) {
    std::cerr << "AuxinModelSimple1::"
	      << "AuxinModelSimple1() "
	      << "Twelve parameters used.\n\n";
    std::cerr << "dA_i/dt = p0 - p1*A_i +p4*Sum_{neigh} (A_n-A_i) +\n" 
	      << "p3*Sum_{neigh} (P_ni*A_n-P_in*A_i)\n\n" 
	      << "dP_i/dt = p5 - p6*P_i\n\n"
	      << "P_in = P_i*X_in/(p_2+Sum_{k,neigh}X_ik)\n";
    exit(0);
  }
  if( indValue.size() != 1 || indValue[0].size() != 3 ) {
    std::cerr << "AuxinModelSimple1::"
	      << "AuxinModelSimple1() "
	      << "Three variable indices are used (auxin,pin,X in wall).\n";
    exit(0);
  }
  //Set the variable values
  //
  setId("AuxinModelSimple1Wall");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "p_auxin";
  tmp[1] = "d_auxin";
  tmp[2] = "p_pol";
  tmp[3] = "T_auxin";
  tmp[4] = "D_auxin";
  tmp[5] = "p_pin";
  tmp[6] = "d_pin";
  setParameterId( tmp );
}

void AuxinModelSimple1Wall::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{  
  size_t numCells = T.numCell();
	size_t aI = variableIndex(0,0);
	size_t pI = variableIndex(0,1);
	size_t xI = variableIndex(0,2);
	std::vector<double> pin;
  assert( aI<cellData[0].size() &&
					pI<cellData[0].size() &&
					xI<cellData[0].size() );
  
  for( size_t i=0 ; i<numCells ; ++i ) {
		
		//Production and degradation
    cellDerivs[i][aI] += parameter(0) - parameter(1)*cellData[i][aI];
		cellDerivs[i][pI] += parameter(5) - parameter(6)*cellData[i][pI];
		
		//Transport
		size_t numWalls=T.cell(i).numWall();
		//Polarization coefficient normalization constant
		double sum=0.0;
		pin.resize( numWalls );
		double minPin=0.0;
		for( size_t n=0 ; n<numWalls ; ++n ) {
			sum += pin[n] = wallData[ T.cell(i).wall(n)->index() ][ xI ];
			if (pin[n]<minPin) {
				minPin = pin[n];
			}
			if (minPin<0.0) {
				sum += numWalls*minPin;
				for( size_t n=0 ; n<numWalls ; ++n ) {
					pin[n] += minPin;
				}
			}
		}
		sum += parameter(2);
		
		for( size_t n=0 ; n<numWalls ; ++n ) {
			if( T.cell(i).wall(n)->cell1() != T.background() &&
					T.cell(i).wall(n)->cell2() != T.background() ) { 
				size_t neighIndex; 
				if( T.cell(i).wall(n)->cell1()->index()==i )
					neighIndex = T.cell(i).wall(n)->cell2()->index();				
				else
					neighIndex = T.cell(i).wall(n)->cell1()->index();				
				double polRate=0.0;			
				if( sum > 0.0 )
					polRate = cellData[i][pI] * pin[n] / sum;
				//pin[i][n+1] = polRate;			
				cellDerivs[i][aI] -= (parameter(3)*polRate+parameter(4))*cellData[i][aI];
				cellDerivs[neighIndex][aI] += (parameter(3)*polRate+parameter(4))*cellData[i][aI];
			}
		}
	}
}

AuxinModelSimple2::
AuxinModelSimple2(std::vector<double> &paraValue, 
									std::vector< std::vector<size_t> > 
									&indValue ) 
{ 
  //Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=32 ) {
    std::cerr << "AuxinModelSimple2::"
	      << "AuxinModelSimple2() "
	      << "32 parameters used (see network.h)\n";
    exit(0);
  }
  if( indValue.size() != 1 || indValue[0].size() != 7 ) {
    std::cerr << "AuxinModelSimple2::"
	      << "AuxinModelSimple2() "
	      << "Seven variable indices are used "
	      << "(auxin,pin,aux1,pid,X,L1,M).\n";
    exit(0);
  }
  //Set the variable values
  //
  setId("AuxinModelSimple2");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "p_auxin(M)";
  tmp[1] = "p_auxin";
  tmp[2] = "d_auxin";
  tmp[3] = "p_pol";
  tmp[4] = "T_auxin";
  tmp[5] = "D_auxin";
  tmp[6] = "p_pin";
  tmp[7] = "d_pin";
  tmp[8] = "p_X";
  tmp[9] = "d_X";
  tmp[10] = "p_M";
  tmp[11] = "d_M";
  tmp[12] = "d_M";
  tmp[13] = "d_M";
  tmp[14] = "d_M";
  tmp[15] = "d_M";
  tmp[16] = "d_M";
  tmp[17] = "d_M";
  tmp[18] = "d_M";
  tmp[19] = "d_M";
  tmp[20] = "d_M";
  tmp[21] = "d_M";
  tmp[22] = "d_M";
  tmp[23] = "d_M";
  tmp[24] = "d_M";
  tmp[25] = "d_M";
  tmp[26] = "d_M";
  tmp[27] = "d_M";
  tmp[28] = "d_M";
  tmp[29] = "d_M";
  tmp[30] = "D_A";
  tmp[31] = "D_A";
  
  setParameterId( tmp );
}

void AuxinModelSimple2::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{  
  size_t numCells = T.numCell();
  //Setting the indices
  size_t auxinI = variableIndex(0,0);
  size_t pinI = variableIndex(0,1);
  size_t auxI = variableIndex(0,2);
  size_t pidI = variableIndex(0,3);
  size_t xI = variableIndex(0,4);
  size_t l1I = variableIndex(0,5);
  size_t mI;
  mI = variableIndex(0,6);
  assert( auxinI<cellData[0].size() &&
	  pinI<cellData[0].size() &&
	  auxI<cellData[0].size() &&
	  pidI<cellData[0].size() &&
	  xI<cellData[0].size() &&
	  l1I<cellData[0].size() &&
	  mI<cellData[0].size() );
  
  for( size_t i=0 ; i<numCells ; ++i ) {
    
    //Production and degradation
    cellDerivs[i][auxinI] += parameter(0)* 
      ( (1.0-parameter(1)) + parameter(1)*cellData[i][l1I] ) - 
      parameter(2)*cellData[i][auxinI];
    
    cellDerivs[i][pinI] += parameter(3)* 
      ( (1-parameter(4) ) + parameter(4)*cellData[i][l1I] ) * 
      ( (1-parameter(5)) + parameter(5)*cellData[i][auxinI] / 
	(parameter(6)+cellData[i][auxinI]) ) - 
      parameter(7)*cellData[i][pinI];
    
    cellDerivs[i][auxI] += parameter(8)* 
      ( (1-parameter(9) ) + parameter(9)*cellData[i][l1I] ) * 
      ( (1-parameter(10)) + parameter(10)*cellData[i][auxinI] / 
	(parameter(11)+cellData[i][auxinI]) ) - 
      parameter(12)*cellData[i][auxI];
    
    static double KpowN = std::pow(parameter(15),parameter(16)); 
    cellDerivs[i][pidI] += parameter(13)*( (1.0-parameter(14)) + 
					   parameter(14)*KpowN/
					   (KpowN+std::pow(cellData[i][auxinI],parameter(16))))- 
      parameter(17)*cellData[i][pidI];
    
    cellDerivs[i][xI] += parameter(18)*((1.0-parameter(19))+parameter(19)*cellData[i][l1I])*
      cellData[i][auxinI] - parameter(20)*cellData[i][xI];
    
    cellDerivs[i][l1I] -= parameter(22)*cellData[i][l1I];
    if( T.cell(i).isNeighbor(T.background()) )
      cellDerivs[i][l1I] += parameter(21);
    
    //cellDerivs[i][mI] +=...;
    
    //Transport
    //
    size_t numWalls=T.cell(i).numWall();
    //PID factor
    double tmpPow = std::pow(cellData[i][pidI],parameter(30));
    double Ci = tmpPow/(tmpPow+std::pow(parameter(29),parameter(30)));
    
    //Polarization coefficient normalization constant
    double sum=0.0;
    std::vector<double> Pij(numWalls);
    for( size_t n=0 ; n<numWalls ; ++n ) {
      if( T.cell(i).wall(n)->cell1() != T.background() &&
	  T.cell(i).wall(n)->cell2() != T.background() ) { 
	size_t neighI;
	if( T.cell(i).wall(n)->cell1()->index()==i )
	  neighI = T.cell(i).wall(n)->cell2()->index();
	else
	  neighI = T.cell(i).wall(n)->cell1()->index();
	double powX = std::pow(cellData[ neighI ][ xI ],parameter(28));
	double Cij = powX/(std::pow(parameter(27),parameter(28))+powX);
	sum += Pij[n] = (1.0-parameter(25)) + 
	  parameter(25)*(Ci*Cij+(1.0-Ci)*(1.0-Cij));
      }
      else 
	sum += Pij[n] = (1.0-parameter(25));
    }
    //sum /= numWalls;//For adjusting for different num neigh
    sum += parameter(26);
    
    for( size_t n=0 ; n<numWalls ; ++n ) {
      //if( !T.cell(i).isNeighbor(T.background()) ) { 
      if( T.cell(i).wall(n)->cell1() != T.background() &&
	  T.cell(i).wall(n)->cell2() != T.background() ) { 
	size_t neighI; 
	if( T.cell(i).wall(n)->cell1()->index()==i )
	  neighI = T.cell(i).wall(n)->cell2()->index();				
	else
	  neighI = T.cell(i).wall(n)->cell1()->index();				
	double pol=0.0;
	if( sum != 0.0 )
	  pol = cellData[i][pinI] * Pij[n] / sum;
	double transportRate = parameter(23)*cellData[neighI][auxI]*
	  pol*cellData[i][auxinI] /
	  ( (parameter(24)+cellData[i][auxinI])*
	    (cellData[i][auxI]+cellData[neighI][auxI]) );
	cellDerivs[i][auxinI] -= transportRate + 
	  parameter(31)*cellData[i][auxinI];
	cellDerivs[neighI][auxinI] += transportRate +
	  parameter(31)*cellData[i][auxinI];
      }
    }
  }
}

AuxinModelSimple3::
AuxinModelSimple3(std::vector<double> &paraValue, 
		  std::vector< std::vector<size_t> > 
		  &indValue ) 
{ 
  //Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=32 ) {
    std::cerr << "AuxinModelSimple3::"
	      << "AuxinModelSimple3() "
	      << "32 parameters used (see network.h)\n";
    exit(0);
  }
  if( indValue.size() != 1 || indValue[0].size() != 7 ) {
    std::cerr << "AuxinModelSimple3::"
	      << "AuxinModelSimple3() "
	      << "Seven variable indices are used "
	      << "(auxin,pin,aux1,pid,X,L1,M).\n";
    exit(0);
  }
  //Set the variable values
  //
  setId("AuxinModelSimple3");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "p_auxin(M)";
  tmp[1] = "p_auxin";
  tmp[2] = "d_auxin";
  tmp[3] = "p_pol";
  tmp[4] = "T_auxin";
  tmp[5] = "D_auxin";
  tmp[6] = "p_pin";
  tmp[7] = "d_pin";
  tmp[8] = "p_X";
  tmp[9] = "d_X";
  tmp[10] = "p_M";
  tmp[11] = "d_M";
  tmp[12] = "d_M";
  tmp[13] = "d_M";
  tmp[14] = "d_M";
  tmp[15] = "d_M";
  tmp[16] = "d_M";
  tmp[17] = "d_M";
  tmp[18] = "d_M";
  tmp[19] = "d_M";
  tmp[20] = "d_M";
  tmp[21] = "d_M";
  tmp[22] = "d_M";
  tmp[23] = "d_M";
  tmp[24] = "d_M";
  tmp[25] = "d_M";
  tmp[26] = "d_M";
  tmp[27] = "d_M";
  tmp[28] = "d_M";
  tmp[29] = "d_M";
  tmp[30] = "D_A";
  tmp[31] = "D_A";
  
  setParameterId( tmp );
}

void AuxinModelSimple3::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{  
  size_t numCells = T.numCell();
  //Setting the indices
  size_t auxinI = variableIndex(0,0);
  size_t pinI = variableIndex(0,1);
  size_t auxI = variableIndex(0,2);
  size_t pidI = variableIndex(0,3);
  size_t xI = variableIndex(0,4);
  size_t l1I = variableIndex(0,5);
  size_t mI;
  mI = variableIndex(0,6);

  assert( auxinI<cellData[0].size() &&
	  pinI<cellData[0].size() &&
	  auxI<cellData[0].size() &&
	  pidI<cellData[0].size() &&
	  xI<cellData[0].size() &&
	  l1I<cellData[0].size() &&
	  mI<cellData[0].size() );
  
  for( size_t i=0 ; i<numCells ; ++i ) {
    
    //Production and degradation
    cellDerivs[i][auxinI] += parameter(0)* 
      ( (1.0-parameter(1)) + parameter(1)*cellData[i][l1I] ) - 
      parameter(2)*cellData[i][auxinI];
    
    cellDerivs[i][pinI] += parameter(3)* 
      ( (1-parameter(4) ) + parameter(4)*cellData[i][l1I] ) * 
      ( (1-parameter(5)) + parameter(5)*cellData[i][auxinI] / 
	(parameter(6)+cellData[i][auxinI]) ) - 
      parameter(7)*cellData[i][pinI];
    
    cellDerivs[i][auxI] += parameter(8)* 
      ( (1-parameter(9) ) + parameter(9)*cellData[i][l1I] ) * 
      ( (1-parameter(10)) + parameter(10)*cellData[i][auxinI] / 
	(parameter(11)+cellData[i][auxinI]) ) - 
      parameter(12)*cellData[i][auxI];
    
    static double KpowN = std::pow(parameter(15),parameter(16)); 
    cellDerivs[i][pidI] += parameter(13)*( (1.0-parameter(14)) + 
					   parameter(14)*KpowN/
					   (KpowN+std::pow(cellData[i][auxinI],parameter(16))))- 
      parameter(17)*cellData[i][pidI];
    
    cellDerivs[i][xI] += parameter(18)*((1.0-parameter(19))+parameter(19)*cellData[i][l1I])*
      cellData[i][auxinI] - parameter(20)*cellData[i][xI];
    
    cellDerivs[i][l1I] -= parameter(22)*cellData[i][l1I];
    if( T.cell(i).isNeighbor(T.background()) )
      cellDerivs[i][l1I] += parameter(21);
    
    //cellDerivs[i][mI] +=...;
    
    //Transport
    //
    size_t numWalls=T.cell(i).numWall();
    //PID factor
    double tmpPow = std::pow(cellData[i][pidI],parameter(30));
    double Ci = tmpPow/(tmpPow+std::pow(parameter(29),parameter(30)));
    
    //Polarization coefficient normalization constant
    double sum=0.0;
    std::vector<double> Pij(numWalls);
    for( size_t n=0 ; n<numWalls ; ++n ) {
      if( T.cell(i).wall(n)->cell1() != T.background() &&
	  T.cell(i).wall(n)->cell2() != T.background() ) { 
	size_t neighI;
	if( T.cell(i).wall(n)->cell1()->index()==i )
	  neighI = T.cell(i).wall(n)->cell2()->index();
	else
	  neighI = T.cell(i).wall(n)->cell1()->index();
	//double powX = std::pow(cellData[ neighI ][ xI ],parameter(28));
	//double Cij = powX/(std::pow(parameter(27),parameter(28))+powX);
	double Cij = cellData[ neighI ][ xI ];
	sum += Pij[n] = (1.0-parameter(25)) + 
	  parameter(25)*(Ci*Cij+(1.0-Ci)*(1.0-Cij));
	//sum += Pij[n] = (1.0-parameter(25)) + 
	//parameter(25)*cellData[ neighI ][xI];
      }
      else 
	sum += Pij[n] = (1.0-parameter(25));
    }
    //sum /= numWalls;//For adjusting for different num neigh
    sum += parameter(26);
    
    for( size_t n=0 ; n<numWalls ; ++n ) {
      //if( !T.cell(i).isNeighbor(T.background()) ) { 
      if( T.cell(i).wall(n)->cell1() != T.background() &&
	  T.cell(i).wall(n)->cell2() != T.background() ) { 
	size_t neighI; 
	if( T.cell(i).wall(n)->cell1()->index()==i )
	  neighI = T.cell(i).wall(n)->cell2()->index();				
	else
	  neighI = T.cell(i).wall(n)->cell1()->index();				
	double pol=0.0;
	if( sum != 0.0 )
	  pol = cellData[i][pinI] * Pij[n] / sum;
	// 				double transportRate = parameter(23)*cellData[neighI][auxI]*
	// 					pol*cellData[i][auxinI] /
	// 					( (parameter(24)+cellData[i][auxinI])*
	// 						(cellData[i][auxI]+cellData[neighI][auxI]) );
	double transportRate = parameter(23)*pol*cellData[i][auxinI]*
	  cellData[neighI][auxI] / 
	  (cellData[i][auxI]+cellData[neighI][auxI]);
	
	cellDerivs[i][auxinI] -= transportRate + 
	  parameter(31)*cellData[i][auxinI];
	cellDerivs[neighI][auxinI] += transportRate +
	  parameter(31)*cellData[i][auxinI];
      }
    }
  }
}


AuxinModel4::
AuxinModel4(std::vector<double> &paraValue, 
	    std::vector< std::vector<size_t> > 
	    &indValue ) 
{ 
  //Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=37 ) {
    std::cerr << "AuxinModel4::"
	      << "AuxinModel4() "
	      << "37 parameters used (see network.h)\n";
    exit(0);
  }
  if( indValue.size() != 1 || indValue[0].size() != 7 ) {
    std::cerr << "AuxinModel4::"
	      << "AuxinModel4() "
	      << "Seven variable indices are used "
	      << "(auxin,pin,aux1,pid,X,L1,M).\n";
    exit(0);
  }
  //Set the variable values
  //
  setId("AuxinModel4");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "p_auxin(M)";
  tmp[1] = "p_auxin";
  tmp[2] = "d_auxin";
  tmp[3] = "p_pol";
  tmp[4] = "T_auxin";
  tmp[5] = "D_auxin";
  tmp[6] = "p_pin";
  tmp[7] = "d_pin";
  tmp[8] = "p_X";
  tmp[9] = "d_X";
  tmp[10] = "p_M";
  tmp[11] = "d_M";
  tmp[12] = "d_M";
  tmp[13] = "d_M";
  tmp[14] = "d_M";
  tmp[15] = "d_M";
  tmp[16] = "d_M";
  tmp[17] = "d_M";
  tmp[18] = "d_M";
  tmp[19] = "d_M";
  tmp[20] = "d_M";
  tmp[21] = "d_M";
  tmp[22] = "d_M";
  tmp[23] = "d_M";
  tmp[24] = "d_M";
  tmp[25] = "d_M";
  tmp[26] = "d_M";
  tmp[27] = "d_M";
  tmp[28] = "d_M";
  tmp[29] = "d_M";
  tmp[30] = "D_A";
  tmp[31] = "D_A";
  tmp[32] = "D_A";
  tmp[33] = "D_A";
  tmp[34] = "D_A";
  tmp[35] = "D_A";
  tmp[36] = "D_A";
}

void AuxinModel4::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{  
  size_t numCells = T.numCell();
  //Setting the indices
  size_t auxinI = variableIndex(0,0);
  size_t pinI = variableIndex(0,1);
  size_t auxI = variableIndex(0,2);
  size_t pidI = variableIndex(0,3);
  size_t xI = variableIndex(0,4);
  size_t l1I = variableIndex(0,5);
  size_t mI;
  mI = variableIndex(0,6);
  assert( auxinI<cellData[0].size() &&
	  pinI<cellData[0].size() &&
	  auxI<cellData[0].size() &&
	  pidI<cellData[0].size() &&
	  xI<cellData[0].size() &&
	  l1I<cellData[0].size() &&
	  mI<cellData[0].size() );
  
  for( size_t i=0 ; i<numCells ; ++i ) {
    
    //Production and degradation
    double tmpPow=std::pow(cellData[i][auxinI],parameter(4));
    cellDerivs[i][auxinI] += parameter(0)* 
      ( (1.0-parameter(1)) + parameter(1)*cellData[i][l1I] )*
      ((1-parameter(2)) + parameter(2)*tmpPow/
       (tmpPow+std::pow(parameter(3),parameter(4)))) - 
      parameter(5)*cellData[i][auxinI];
    
    cellDerivs[i][pinI] += parameter(6)* 
      ( (1-parameter(7) ) + parameter(7)*cellData[i][l1I] ) * 
      ( (1-parameter(8)) + parameter(8)*cellData[i][auxinI] / 
	(parameter(9)+cellData[i][auxinI]) ) - 
      parameter(10)*cellData[i][pinI];
    
    cellDerivs[i][auxI] += parameter(11)* 
      ( (1-parameter(12) ) + parameter(12)*cellData[i][l1I] ) * 
      ( (1-parameter(13)) + parameter(13)*cellData[i][auxinI] / 
	(parameter(14)+cellData[i][auxinI]) ) - 
      parameter(15)*cellData[i][auxI];
    
    static double KpowN = std::pow(parameter(19),parameter(20)); 
    tmpPow = std::pow(cellData[i][auxinI],parameter(20));
    cellDerivs[i][pidI] += parameter(16)*( (1-parameter(17)) +
					   parameter(17)*
					   cellData[i][l1I] )*
      ( (1.0-parameter(18)) + 
	parameter(18)*
	KpowN/(KpowN+tmpPow) ) - 
      parameter(21)*cellData[i][pidI];
    
    cellDerivs[i][xI] += parameter(22)*((1.0-parameter(23))+parameter(23)*cellData[i][l1I])*
      cellData[i][auxinI] - parameter(24)*cellData[i][xI];
    
    cellDerivs[i][l1I] += parameter(25)-parameter(27)*cellData[i][l1I];
    if( T.cell(i).isNeighbor(T.background()) )
      cellDerivs[i][l1I] += parameter(26);
    
    //cellDerivs[i][mI] +=...;
    
    //Transport
    //
    size_t numWalls=T.cell(i).numWall();
    //PID factor
    tmpPow = std::pow(cellData[i][pidI],parameter(35));
    double Ci = tmpPow/(tmpPow+std::pow(parameter(34),parameter(35)));
    //Find max Cij (tmp solution for liear Cij)
    double maxCij=0.0;
    for( size_t n=0 ; n<numWalls ; ++n ) {
      if( T.cell(i).wall(n)->cell1() != T.background() &&
	  T.cell(i).wall(n)->cell2() != T.background() ) { 
	size_t neighI;
	if( T.cell(i).wall(n)->cell1()->index()==i )
	  neighI = T.cell(i).wall(n)->cell2()->index();
	else
	  neighI = T.cell(i).wall(n)->cell1()->index();
	double Cij = cellData[ neighI ][ xI ];
	if( Cij>maxCij ) maxCij=Cij;
      }
    }

    //Polarization coefficient normalization constant
    double sum=0.0;
    std::vector<double> Pij(numWalls);
    for( size_t n=0 ; n<numWalls ; ++n ) {
      if( T.cell(i).wall(n)->cell1() != T.background() &&
	  T.cell(i).wall(n)->cell2() != T.background() ) { 
	size_t neighI;
	if( T.cell(i).wall(n)->cell1()->index()==i )
	  neighI = T.cell(i).wall(n)->cell2()->index();
	else
	  neighI = T.cell(i).wall(n)->cell1()->index();
	//double powX = std::pow(cellData[ neighI ][ xI ],parameter(33));
	//double Cij = powX/(std::pow(parameter(32),parameter(33))+powX);
	double Cij = cellData[ neighI ][ xI ];
	sum += Pij[n] = (1.0-parameter(30)) + 
	  parameter(30)*(Ci*Cij+(1.0-Ci)*(maxCij-Cij));
	//sum += Pij[n] = (1.0-parameter(30)) + 
	//parameter(30)*cellData[ neighI ][xI];
      }
      else 
	sum += Pij[n] = (1.0-parameter(30));
    }
    //sum /= numWalls;//For adjusting for different num neigh
    sum += parameter(31);
    
    for( size_t n=0 ; n<numWalls ; ++n ) {
      //if( !T.cell(i).isNeighbor(T.background()) ) { 
      if( T.cell(i).wall(n)->cell1() != T.background() &&
	  T.cell(i).wall(n)->cell2() != T.background() ) { 
	size_t neighI; 
	if( T.cell(i).wall(n)->cell1()->index()==i )
	  neighI = T.cell(i).wall(n)->cell2()->index();				
	else
	  neighI = T.cell(i).wall(n)->cell1()->index();				
	double pol=0.0;
	if( sum != 0.0 )
	  pol = cellData[i][pinI] * Pij[n] / sum;
	// 				double transportRate = parameter(28)*cellData[neighI][auxI]*
	// 					pol*cellData[i][auxinI] /
	// 					( (parameter(29)+cellData[i][auxinI])*
	// 						(cellData[i][auxI]+cellData[neighI][auxI]) );
	double transportRate = parameter(28)*pol*cellData[i][auxinI]*
	  cellData[neighI][auxI] / 
	  (cellData[i][auxI]+cellData[neighI][auxI]);
	
	cellDerivs[i][auxinI] -= transportRate + 
	  parameter(36)*cellData[i][auxinI];
	cellDerivs[neighI][auxinI] += transportRate +
	  parameter(36)*cellData[i][auxinI];
      }
    }
  }
}

AuxinModel5::
AuxinModel5(std::vector<double> &paraValue, 
	    std::vector< std::vector<size_t> > 
	    &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=16 ) {
    std::cerr << "AuxinModel5::"
	      << "AuxinModel5() "
	      << "Sixteen parameters used (see network.h)\n";
    exit(0);
  }
  if( indValue.size() != 1 || indValue[0].size() != 4 ) {
    std::cerr << "AuxinModel5::"
	      << "AuxinModel5() "
	      << "Four variable indices are used (auxin,pin,X,M).\n";
    exit(0);
  }
  //Set the variable values
  //
  setId("AuxinModel5");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "p_auxin(M)";
  tmp[1] = "p_auxin";
  tmp[2] = "d_auxin";
  tmp[3] = "p_pol";
  tmp[4] = "T_auxin";
  tmp[5] = "D_auxin";
  tmp[6] = "p_pin";
  tmp[7] = "d_pin";
  tmp[8] = "p_X";
  tmp[9] = "d_X";
  tmp[10] = "p_M";
  tmp[11] = "K_M";
  tmp[12] = "n_M";	
  tmp[13] = "K_M2";
  tmp[14] = "n_M2";	
  tmp[15] = "d_M";
  
  setParameterId( tmp );
}

void AuxinModel5::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{  
  size_t numCells = T.numCell();
  size_t aI = variableIndex(0,0);
  size_t pI = variableIndex(0,1);
  size_t xI = variableIndex(0,2);
  size_t mI = variableIndex(0,3);
  assert( aI<cellData[0].size() &&
	  pI<cellData[0].size() &&
	  xI<cellData[0].size() &&
	  mI<cellData[0].size() );
  size_t dimension = vertexData[0].size();
  double powK = std::pow(parameter(11),parameter(12));
  double powK2 = std::pow(parameter(13),parameter(14));
  
  for( size_t i=0 ; i<numCells ; ++i ) {
    
    //Production and degradation
    cellDerivs[i][aI] += parameter(0)*cellData[i][mI] + parameter(1) - 
      parameter(2)*cellData[i][aI];
    
    cellDerivs[i][pI] += parameter(6) - parameter(7)*cellData[i][pI];
    
    //cellDerivs[i][xI] += parameter(8)*cellData[i][aI]*cellData[i][aI]/
    //(2.0+cellData[i][aI]*cellData[i][aI]) 
    //- parameter(9)*cellData[i][xI];
    cellDerivs[i][xI] += parameter(8)*cellData[i][aI]
      - parameter(9)*cellData[i][xI];

    
    std::vector<double> cellCenter = T.cell(i).positionFromVertex(vertexData);
    double R=0.0;
    size_t dimIndexMax = dimension;
    if (dimension==3)
      dimIndexMax -= 1;
    for (size_t d=0; d<dimIndexMax; ++d)
      R += cellCenter[d]*cellCenter[d];
    R = std::sqrt(R);
    double powR = std::pow(R,parameter(12));
    double powR2 = std::pow(R,parameter(14));
    cellDerivs[i][mI] += parameter(10)*(powR*powK2) / ((powK+powR)*(powK2+powR2)) 
      - parameter(15)*cellData[i][mI];
    
    //Transport
    size_t numWalls=T.cell(i).numWall();
    //Polarization coefficient normalization constant
    double sum=0.0;
    size_t numActualWalls=0;
    //pin[i].resize( numWalls+1 );
    for( size_t n=0 ; n<numWalls ; ++n ) {
      if( T.cell(i).wall(n)->cell1() != T.background() &&
	  T.cell(i).wall(n)->cell2() != T.background() ) { 
	numActualWalls++;
	if( T.cell(i).wall(n)->cell1()->index()==i )
	  sum += cellData[ T.cell(i).wall(n)->cell2()->index() ][ xI ];
	else
	  sum += cellData[ T.cell(i).wall(n)->cell1()->index() ][ xI ];
      }
      else 
	sum += 0.5;
    }
    //sum /= numActualWalls;//For adjusting for different num neigh
    sum += parameter(3);
    
    for( size_t n=0 ; n<numWalls ; ++n ) {
      //if( !T.cell(i).isNeighbor(T.background()) ) { 
      if( T.cell(i).wall(n)->cell1() != T.background() &&
	  T.cell(i).wall(n)->cell2() != T.background() ) { 
	size_t neighIndex; 
	if( T.cell(i).wall(n)->cell1()->index()==i )
	  neighIndex = T.cell(i).wall(n)->cell2()->index();				
	else
	  neighIndex = T.cell(i).wall(n)->cell1()->index();				
	double polRate=0.0;
	
	if( sum != 0.0 )
	  polRate = cellData[i][pI] * cellData[neighIndex][xI] / sum;
	//pin[i][n+1] = polRate;
	cellDerivs[i][aI] -= (parameter(4)*polRate+parameter(5))*cellData[i][aI];
	cellDerivs[neighIndex][aI] += (parameter(4)*polRate+parameter(5))*cellData[i][aI];
      }
    }
  }
}

AuxinModel6::
AuxinModel6(std::vector<double> &paraValue, 
						std::vector< std::vector<size_t> > 
						&indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=16 ) {
    std::cerr << "AuxinModel6::"
							<< "AuxinModel6() "
							<< "16 parameters used (see network.h)\n";
    exit(0);
  }
  if( indValue.size() != 1 || indValue[0].size() != 5 ) {
    std::cerr << "AuxinModel6::"
							<< "AuxinModel6() "
							<< "Five variable indices are used (auxin,pin,aux,X,M).\n";
    exit(0);
  }
  //Set the variable values
  //
  setId("AuxinModel6");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "p_auxin(M)";
  tmp[1] = "p0_auxin";
  tmp[2] = "d_auxin";
  tmp[3] = "p0_pin";
  tmp[4] = "d_pin";
  tmp[5] = "p0_aux";
  tmp[6] = "d_aux";
  tmp[7] = "p0_X";
  tmp[8] = "p_X(auxin)";
  tmp[9] = "d_X";
  tmp[10] = "p_M";
  tmp[11] = "K_M";
  tmp[12] = "n_M";	
  tmp[13] = "K_M2";
  tmp[14] = "n_M2";	
  tmp[15] = "d_M";
	
  setParameterId( tmp );
}

void AuxinModel6::
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
  size_t pI = variableIndex(0,1);//pin
  size_t AI = variableIndex(0,2);//aux
  size_t xI = variableIndex(0,3);//X
  size_t mI = variableIndex(0,4);//M
  assert( aI<cellData[0].size() &&
	  pI<cellData[0].size() &&
	  AI<cellData[0].size() &&
	  xI<cellData[0].size() &&
	  mI<cellData[0].size() );
  size_t dimension = vertexData[0].size();
  double powK = std::pow(parameter(11),parameter(12));
  double powK2 = std::pow(parameter(13),parameter(14));
  
  for( size_t i=0 ; i<numCells ; ++i ) {
    
    //Production and degradation
    cellDerivs[i][aI] += parameter(0)*cellData[i][mI] + parameter(1) - 
      parameter(2)*cellData[i][aI];
    
    cellDerivs[i][pI] += parameter(3) - parameter(4)*cellData[i][pI];
    
    cellDerivs[i][AI] += parameter(5) - parameter(6)*cellData[i][pI];
    
    cellDerivs[i][xI] += parameter(7) + parameter(8)*cellData[i][aI]*
      cellData[i][aI]/(2.0+cellData[i][aI]*cellData[i][aI]) 
      - parameter(9)*cellData[i][xI];
    
    std::vector<double> cellCenter = T.cell(i).positionFromVertex(vertexData);
    double R=0.0;
    for (size_t d=0; d<dimension; ++d)
      R += cellCenter[d]*cellCenter[d];
    R = std::sqrt(R);
    double powR = std::pow(R,parameter(12));
    double powR2 = std::pow(R,parameter(14));
    cellDerivs[i][mI] += parameter(10)*(powR*powK2) / ((powK+powR)*(powK2+powR2)) 
      - parameter(15)*cellData[i][mI];
  }
}

AuxinModel7::
AuxinModel7(std::vector<double> &paraValue, 
	    std::vector< std::vector<size_t> > 
	    &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=22 ) {
    std::cerr << "AuxinModel7::"
	      << "AuxinModel7() "
	      << "22 parameters used (see network.h)\n";
    exit(0);
  }
  if( indValue.size() != 1 || indValue[0].size() != 5 ) {
    std::cerr << "AuxinModel7::"
	      << "AuxinModel7() "
	      << "Five variable indices are used (auxin,pin,aux,X,M).\n";
    exit(0);
  }
  //Set the variable values
  //
  setId("AuxinModel7");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "p_auxin(M)";
  tmp[1] = "p0_auxin";
  tmp[2] = "d_auxin";
  tmp[3] = "p0_pin";
  tmp[4] = "f_Pauxin";
  tmp[5] = "K_P";
  tmp[6] = "n_P";
  tmp[7] = "d_pin";
  tmp[8] = "p0_aux";
  tmp[9] = "f_Aauxin";
  tmp[10] = "K_A";
  tmp[11] = "n_A";
  tmp[12] = "d_aux";
  tmp[13] = "p0_X";
  tmp[14] = "p_X(auxin)";
  tmp[15] = "d_X";
  tmp[16] = "p_M";
  tmp[17] = "K_M";
  tmp[18] = "n_M";	
  tmp[19] = "K_M2";
  tmp[20] = "n_M2";	
  tmp[21] = "d_M";
	
  setParameterId( tmp );
}

void AuxinModel7::
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
	size_t pI = variableIndex(0,1);//pin
	size_t AI = variableIndex(0,2);//aux
	size_t xI = variableIndex(0,3);//X
	size_t mI = variableIndex(0,4);//M
  assert( aI<cellData[0].size() &&
					pI<cellData[0].size() &&
					AI<cellData[0].size() &&
					xI<cellData[0].size() &&
					mI<cellData[0].size() );
  size_t dimension = vertexData[0].size();
	double powK = std::pow(parameter(17),parameter(18));
	double powK2 = std::pow(parameter(19),parameter(20));
	double powPK = std::pow(parameter(5),parameter(6));
	double powAK = std::pow(parameter(10),parameter(11));
	
  for( size_t i=0 ; i<numCells ; ++i ) {
		
		//Production and degradation
    cellDerivs[i][aI] += parameter(0)*cellData[i][mI] + parameter(1) - 
			parameter(2)*cellData[i][aI];
		
		double powPA = std::pow(cellData[i][aI],parameter(6));
    cellDerivs[i][pI] += parameter(3)*((1-parameter(4))+ parameter(4)*powPA/(powPK+powPA))
			- parameter(7)*cellData[i][pI];
		
		double powAA = std::pow(cellData[i][aI],parameter(11));
    cellDerivs[i][AI] += parameter(8)*((1-parameter(9))+parameter(9)*powAA/(powAK+powAA)) 
			- parameter(12)*cellData[i][pI];
		
    cellDerivs[i][xI] += parameter(13) + parameter(14)*cellData[i][aI]*
			cellData[i][aI]/(1.0+cellData[i][aI]*cellData[i][aI]) 
			- parameter(15)*cellData[i][xI];
		
		std::vector<double> cellCenter = T.cell(i).positionFromVertex(vertexData);
		double R=0.0;
		for (size_t d=0; d<dimension; ++d)
			R += cellCenter[d]*cellCenter[d];
		R = std::sqrt(R);
		double powR = std::pow(R,parameter(18));
		double powR2 = std::pow(R,parameter(20));
		cellDerivs[i][mI] += parameter(16)*(powR*powK2) / ((powK+powR)*(powK2+powR2)) 
			- parameter(21)*cellData[i][mI];
	}
}

AuxinTransportCellCellNoGeometry::
AuxinTransportCellCellNoGeometry(std::vector<double> &paraValue, 
				 std::vector< std::vector<size_t> > &indValue )
{
  if (paraValue.size() != 7) {
    std::cerr << "AuxinTransportCellCellNoGeometry::AuxinTransportCellCellNoGeometry"
	      << "Uses seven parameters: " << std::endl
	      << "d, "             // 0
	      << "T, "             // 1
	      << "k1, "            // 2
	      << "k2, "            // 3
	      << "K_H, "           // 4
	      << "n and "          // 5
	      << "K_M."            // 6
	      << std::endl;
    exit(EXIT_FAILURE);
  }
  
  if (indValue.size() != 1 || indValue[0].size() != 4) {
    std::cerr << "AuxinTransportCellCellNoGeometry::AuxinTransportCellCellNoGeometry "
	      << "Four variable indices are needed: " << std::endl
	      << "Auxin, "      // 0
	      << "PIN1, "       // 1
	      << "AUX1, and "   // 2
	      << "X(feedback)." // 3
	      << std::endl;
    exit(EXIT_FAILURE);
  }
  
  setId("AuxinTransportCellCellNoGeometry");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  std::vector<std::string> tmp(numParameter());
  tmp.resize(numParameter());
  
  tmp[0] = "d";
  tmp[1] = "T";
  tmp[2] = "k1";
  tmp[3] = "k2";
  tmp[4] = "K_H";
  tmp[5] = "n";
  tmp[6] = "K_M";
  
  setParameterId( tmp );
}

void AuxinTransportCellCellNoGeometry::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{  
  size_t aI = variableIndex(0,0);
  size_t PI = variableIndex(0,1);
  size_t AI = variableIndex(0,2);
  size_t xI = variableIndex(0,3);
  
  for (size_t i=0; i<cellData.size(); ++i) {
    size_t numActualWalls=0;
    double sum = 1.0;
    std::vector<double> Pin(T.cell(i).numWall());
    for (size_t n = 0; n < T.cell(i).numWall(); ++n) {
      if( T.cell(i).wall(n)->cell1() != T.background() &&
	  T.cell(i).wall(n)->cell2() != T.background() ) { 
	numActualWalls++;
	if( T.cell(i).wall(n)->cell1()->index()==i ) {
	  double tmp = std::pow(cellData[ T.cell(i).wall(n)->cell2()->index() ][ xI ],parameter(5));
	  tmp = parameter(2) + 
	    parameter(3)*tmp/(std::pow(parameter(4),parameter(5))+tmp);
	  sum += Pin[n] = tmp;
	}
	else {
	  double tmp = std::pow(cellData[ T.cell(i).wall(n)->cell1()->index() ][ xI ],parameter(5));
	  tmp = parameter(2) + 
	    parameter(3)*tmp/(std::pow(parameter(4),parameter(5))+tmp);
	  sum += Pin[n] = tmp;
	}
      }
      else 
	sum += Pin[n] = parameter(2);
    }
    
    for (size_t n = 0; n < T.cell(i).numWall(); ++n) {
      if( T.cell(i).wall(n)->cell1() != T.background() &&
	  T.cell(i).wall(n)->cell2() != T.background() ) { 
	size_t j=T.cell(i).wall(n)->cell1()->index();
	if( T.cell(i).wall(n)->cell1()->index()==i )
	  j = T.cell(i).wall(n)->cell2()->index();
	
	double passive = parameter(0) * cellData[i][aI];
	
	double Pij = cellData[i][PI] * Pin[n] / sum;		
	
	double active = Pij*parameter(1)*cellData[j][AI]*cellData[i][aI]/ 
	  ((parameter(6) + cellData[i][aI])*(cellData[i][AI]+cellData[j][AI]));
	
	cellDerivs[i][aI] -= (passive + active);
	cellDerivs[j][aI] += (passive + active);
      }
    }
  }
}	

AuxinWallModel::
AuxinWallModel(std::vector<double> &paraValue, 
	      std::vector< std::vector<size_t> > 
	      &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=14 ) {
    std::cerr << "AuxinWallModel::"
	      << "AuxinWallModel() "
	      << "14 parameters used (see documentation or network.h)\n";
    exit(EXIT_FAILURE);
  }
  if( indValue.size() != 2 || indValue[0].size() != 4 || indValue[1].size() != 2 ) {
    std::cerr << "AuxinWallModel::"
	      << "AuxinWallModel() "
	      << "Four cell variable indices (auxin, PIN, AUX, X(feedback); first row) and two wall variable"
	      << " indices (paired) are used (auxin, PIN)." << std::endl;
    exit(EXIT_FAILURE);
  }
  //Set the variable values
  //
  setId("AuxinWallModel");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "c_IAA";
  tmp[1] = "d_IAA";
  tmp[2] = "p_IAAH(in)";
  tmp[3] = "p_IAA-(in,AUX)";
  tmp[4] = "p_IAAH(out)";
  tmp[5] = "p_IAA-(out,PIN)";
  tmp[6] = "D_IAA";
  tmp[7] = "c_PIN";
  tmp[8] = "d_PIN";
  tmp[9] = "endo_PIN";
  tmp[10] = "exo_PIN_const";
  tmp[11] = "exo_PIN_auxin";
  tmp[12] = "c_X";
  tmp[13] = "d_X";
  setParameterId( tmp );
}

void AuxinWallModel::
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
  size_t pI = variableIndex(0,1);//PIN
  size_t auxI = variableIndex(0,2);//AUX
  size_t xI = variableIndex(0,3);//X
  size_t awI = variableIndex(1,0);//auxin (wall)
  size_t pwI = variableIndex(1,1);//PIN (membrane/wall)

  assert( aI<cellData[0].size() &&
	  pI<cellData[0].size() &&
	  auxI<cellData[0].size() &&
	  xI<cellData[0].size() &&
	  awI<wallData[0].size() &&
	  pwI<wallData[0].size() );

  for (size_t i=0; i<numCells; ++i) {
    
    //Production and degradation
    cellDerivs[i][aI] += parameter(0) - parameter(1)*cellData[i][aI];
    cellDerivs[i][pI] += parameter(7) - parameter(8)*cellData[i][pI];
    cellDerivs[i][xI] += parameter(11)*cellData[i][aI] - parameter(12)*cellData[i][xI];
    
    //Auxin transport and protein cycling
    size_t numWalls = T.cell(i).numWall();
    for (size_t k=0; k<numWalls; ++k) {
      size_t j = T.cell(i).wall(k)->index();
      if( T.cell(i).wall(k)->cell1()->index() == i && T.cell(i).wall(k)->cell2() != T.background() ) {
	size_t cellNeigh = T.cell(i).wall(k)->cell2()->index();
	// cell-wall transport
	double fac = (parameter(4)+parameter(5)*wallData[j][pwI]) * cellData[i][aI] 
	  - (parameter(2)+parameter(3)*cellData[i][auxI]) * wallData[j][awI];
	wallDerivs[j][awI] += fac;
	cellDerivs[i][aI] -= fac;
	// wall-wall diffusion
	fac = parameter(6)*wallData[j][awI];
	wallDerivs[j][awI] -= fac;
	wallDerivs[j][awI+1] += fac;
	
	//PIN cycling
	fac = parameter(9)*wallData[j][pwI] - 
	  (parameter(10)+parameter(11)*cellData[cellNeigh][xI])*cellData[i][pI];
	wallDerivs[j][pwI] -= fac;
	cellDerivs[i][pI] += fac;
      }
      else if( T.cell(i).wall(k)->cell2()->index() == i && T.cell(i).wall(k)->cell1() != T.background() ) {
	size_t cellNeigh = T.cell(i).wall(k)->cell1()->index();
	// cell-wall transport
	double fac = (parameter(4)+parameter(5)*wallData[j][pwI+1]) * cellData[i][aI] 
	  - (parameter(2) + parameter(3)*cellData[i][auxI]) * wallData[j][awI+1];
		
	wallDerivs[j][awI+1] += fac;
	cellDerivs[i][aI] -= fac;
	// wall-wall diffusion
	wallDerivs[j][awI+1] -= parameter(6)*wallData[j][awI+1];
	wallDerivs[j][awI] += parameter(6)*wallData[j][awI+1];

	//PIN cycling
	fac = parameter(9)*wallData[j][pwI+1] - 
	  (parameter(10)+parameter(11)*cellData[cellNeigh][xI])*cellData[i][pI];
	wallDerivs[j][pwI+1] -= fac;
	cellDerivs[i][pI] += fac;
      }
      //else {
      //std::cerr << "AuxinWallModel::derivs() Cell-wall neighborhood wrong." 
      //	  << std::endl;
      //exit(-1);
      //}
    }
  }
}

AuxinROPModel::
AuxinROPModel(std::vector<double> &paraValue, 
	      std::vector< std::vector<size_t> > 
	      &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=16 ) {
    std::cerr << "AuxinROPModel::"
	      << "AuxinROPModel() "
	      << "16 parameters used (see network.h)\n";
    exit(0);
  }
  if( indValue.size() != 2 || indValue[0].size() != 3 || indValue[1].size() != 3 ) {
    std::cerr << "AuxinROPModel::"
	      << "AuxinROPModel() "
	      << "Three cell variable indices (first row) and three wall variable"
	      << " indices are used (auxin,PIN,ROP)." << std::endl;
    exit(0);
  }
  //Set the variable values
  //
  setId("AuxinROPModel");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "c_IAA";
  tmp[1] = "d_IAA";
  tmp[2] = "p_IAAH(in)";
  tmp[3] = "p_IAAH(out)";
  tmp[4] = "p_IAA-";
  tmp[5] = "D_IAA";
  tmp[6] = "c_PIN";
  tmp[7] = "d_PIN";
  tmp[8] = "endo_PIN";
  tmp[9] = "exo_PIN";
  tmp[10] = "c_ROP";
  tmp[11] = "d_ROP";
  tmp[12] = "endo_ROP";
  tmp[13] = "exo_ROP";
  tmp[14] = "K_hill";
  tmp[15] = "n_hill";
	
  setParameterId( tmp );
}

void AuxinROPModel::
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
  size_t pI = variableIndex(0,1);//pin
  size_t rI = variableIndex(0,2);//rop
  size_t awI = variableIndex(1,0);//auxin (wall)
  size_t pwI = variableIndex(1,1);//pin (membrane/wall)
  size_t rwI = variableIndex(1,2);//rop (membrane/wall)

  assert( aI<cellData[0].size() &&
	  pI<cellData[0].size() &&
	  rI<cellData[0].size() &&
	  awI<wallData[0].size() &&
	  pwI<wallData[0].size() &&
	  rwI<wallData[0].size() );

  for (size_t i=0; i<numCells; ++i) {
	  
    //Production and degradation
    cellDerivs[i][aI] += parameter(0) - parameter(1)*cellData[i][aI];
    cellDerivs[i][pI] += parameter(6) - parameter(7)*cellData[i][pI];
    cellDerivs[i][rI] += parameter(10) - parameter(11)*cellData[i][rI];
    
    //Auxin transport and protein cycling
    size_t numWalls = T.cell(i).numWall();
    for (size_t k=0; k<numWalls; ++k) {
      size_t j = T.cell(i).wall(k)->index();
      if( T.cell(i).wall(k)->cell1()->index() == i && T.cell(i).wall(k)->cell2() != T.background() ) {
	// cell-wall transport
	double fac = parameter(3)*cellData[i][aI] - parameter(2)*wallData[j][awI] +
	  parameter(4)*cellData[i][aI]*wallData[j][pwI];
	
	wallDerivs[j][awI] += fac;
	cellDerivs[i][aI] -= fac;
	// wall-wall diffusion
	wallDerivs[j][awI] -= parameter(5)*wallData[j][awI];
	wallDerivs[j][awI+1] += parameter(5)*wallData[j][awI];
	
	//PIN cycling
	fac = parameter(8)*wallData[j][pwI] - parameter(9)*cellData[i][pI]*wallData[j][rwI];
	wallDerivs[j][pwI] -= fac;
	cellDerivs[i][pI] += fac;

	//ROP cycling
	//fac = parameter(12)*wallData[j][rwI] - parameter(13)*cellData[i][rI]*wallData[j][awI];
	fac = parameter(12)*wallData[j][rwI]
	  *std::pow(wallData[j][rwI+1],parameter(15))/
	  ( std::pow(parameter(14),parameter(15)) + std::pow(wallData[j][rwI+1],parameter(15)) ) -
	  parameter(13)*cellData[i][rI]*wallData[j][awI];
	wallDerivs[j][rwI] -= fac;
	cellDerivs[i][rI] += fac;

      }
      else if( T.cell(i).wall(k)->cell2()->index() == i && T.cell(i).wall(k)->cell1() != T.background() ) {
	// cell-wall transport
	double fac = parameter(3)*cellData[i][aI] - parameter(2)*wallData[j][awI+1] +
	  parameter(4)*cellData[i][aI]*wallData[j][pwI+1];
	
	wallDerivs[j][awI+1] += fac;
	cellDerivs[i][aI] -= fac;
	// wall-wall diffusion
	wallDerivs[j][awI+1] -= parameter(5)*wallData[j][awI+1];
	wallDerivs[j][awI] += parameter(5)*wallData[j][awI+1];

	//PIN cycling
	fac = parameter(8)*wallData[j][pwI+1] - parameter(9)*cellData[i][pI]*wallData[j][rwI+1];
	wallDerivs[j][pwI+1] -= fac;
	cellDerivs[i][pI] += fac;

	//ROP cycling
	//fac = parameter(12)*wallData[j][rwI+1] - parameter(13)*cellData[i][rI]*wallData[j][awI+1];
	fac = parameter(12)*wallData[j][rwI+1]
	  *std::pow(wallData[j][rwI],parameter(15))/
	  ( std::pow(parameter(14),parameter(15)) + std::pow(wallData[j][rwI],parameter(15)) ) -
	  parameter(13)*cellData[i][rI]*wallData[j][awI+1];
	wallDerivs[j][rwI+1] -= fac;
	cellDerivs[i][rI] += fac;

      }
      //else {
      //std::cerr << "AuxinROPModel::derivs() Cell-wall neighborhood wrong." 
      //	  << std::endl;
      //exit(-1);
      //}
    }
  }
}

AuxinROPModel2::
AuxinROPModel2(std::vector<double> &paraValue,
               std::vector< std::vector<size_t> >
               &indValue ) {

  //Do some checks on the parameters and variable indeces
  //       
  if( paraValue.size()!=19 ) {
    std::cerr << "AuxinROPModel2::"
              << "AuxinROPModel2() "
              << "19 parameters used (see network.h)\n";
    exit(0);
  }
  if( indValue.size() != 2 || indValue[0].size() != 3 || indValue[1].size() != 3 ) {
    std::cerr << "AuxinROPModel2::"
              << "AuxinROPModel2() "
              << "Three cell variable indices (first row) and three wall variable"
              << " indices are used (auxin,PIN,ROP)." << std::endl;
    exit(0);
  }
  //Set the variable values
  //
  setId("AuxinROPModel2");
  setParameter(paraValue);
  setVariableIndex(indValue);

  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "c_IAA";
  tmp[1] = "d_IAA";
  tmp[2] = "p_IAAH(in)";
  tmp[3] = "p_IAAH(out)";
  tmp[4] = "p_IAA-";
  tmp[5] = "K_A";
  tmp[6] = "D_IAA";
  tmp[7] = "c_PIN";
  tmp[8] = "d_PIN";
  tmp[9] = "endo_PIN";
  tmp[10] = "exo_PIN";
  tmp[11] = "K_RP";
  tmp[12] = "m";
  tmp[13] = "c_ROP";
  tmp[14] = "d_ROP";
  tmp[15] = "endo_ROP";
  tmp[16] = "K_RR";
  tmp[17] = "n";
  tmp[18] = "exo_ROP";
  setParameterId( tmp );
}

void AuxinROPModel2::
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
  size_t pI = variableIndex(0,1);//pin                                                                                               
  size_t rI = variableIndex(0,2);//rop                                                                                               
  size_t awI = variableIndex(1,0);//auxin (wall)                                                                                     
  size_t pwI = variableIndex(1,1);//pin (membrane/wall)                                                                              
  size_t rwI = variableIndex(1,2);//rop (membrane/wall)                                                                              

  assert( aI<cellData[0].size() &&
          pI<cellData[0].size() &&
          rI<cellData[0].size() &&
          awI<wallData[0].size() &&
          pwI<wallData[0].size() &&
          rwI<wallData[0].size() );

  for (size_t i=0; i<numCells; ++i) {

    //Production and degradation                                                               
    cellDerivs[i][aI] += parameter(0) - parameter(1)*cellData[i][aI];
    cellDerivs[i][pI] += parameter(7) - parameter(8)*cellData[i][pI];
    cellDerivs[i][rI] += parameter(13) - parameter(14)*cellData[i][rI];

    //Auxin transport and protein cycling                                                      
    size_t numWalls = T.cell(i).numWall();
    for (size_t k=0; k<numWalls; ++k) {
      size_t j = T.cell(i).wall(k)->index();
      if( T.cell(i).wall(k)->cell1()->index() == i && T.cell(i).wall(k)->cell2() !=   
	  T.background() ) {
        // cell-wall transport                                                                 
        double fac = parameter(3)*cellData[i][aI] - parameter(2)*wallData[j][awI] +
          parameter(4)*wallData[j][pwI]*cellData[i][aI]/(parameter(5) + cellData[i][aI]);

        wallDerivs[j][awI] += fac;
        cellDerivs[i][aI] -= fac;
	// wall-wall diffusion
        wallDerivs[j][awI] -= parameter(6)*wallData[j][awI];
        wallDerivs[j][awI+1] += parameter(6)*wallData[j][awI];

        //PIN cycling                                                                          
        fac = parameter(9)*wallData[j][pwI]-
	  parameter(10)*cellData[i][pI]*std::pow(wallData[j][rwI],parameter(12))/
	  ( std::pow(parameter(11),parameter(12)) + std::pow(wallData[j][rwI],parameter(12)) );
        wallDerivs[j][pwI] -= fac;
        cellDerivs[i][pI] += fac;
	
        //ROP cycling                                                      
        fac = parameter(15)*wallData[j][rwI] * std::pow(wallData[j][rwI+1],parameter(17))/   
	  ( std::pow(parameter(16),parameter(17)) + std::pow(wallData[j][rwI+1],parameter(17)) ) - 
	  parameter(18)*cellData[i][rI]*wallData[j][awI];
	wallDerivs[j][rwI] -= fac;
	cellDerivs[i][rI] += fac;	
      }
      else if( T.cell(i).wall(k)->cell2()->index() == i && T.cell(i).wall(k)->cell1()
	       != T.background() ) {
	// cell-wall transport                                                                 
	double fac = parameter(3)*cellData[i][aI] - parameter(2)*wallData[j][awI+1] +
          parameter(4)*wallData[j][pwI+1]*cellData[i][aI]/
	  ( parameter(5) + cellData[i][aI] );
	
	wallDerivs[j][awI+1] += fac;
	cellDerivs[i][aI] -= fac;
	// wall-wall diffusion
	wallDerivs[j][awI+1] -= parameter(6)*wallData[j][awI+1];
	wallDerivs[j][awI] += parameter(6)*wallData[j][awI+1];
	
	//PIN cycling                                                                          
	fac = parameter(9)*wallData[j][pwI+1]-
	  parameter(10)*cellData[i][pI]*std::pow(wallData[j][rwI+1],parameter(12))/
	  ( std::pow(parameter(11),parameter(12)) + std::pow(wallData[j][rwI+1],parameter(12)) );
	wallDerivs[j][pwI+1] -= fac;
	cellDerivs[i][pI] += fac;
	
	//ROP cycling                                                                          
	fac = parameter(15)*wallData[j][rwI+1]
	  *std::pow(wallData[j][rwI],parameter(17))/
	  ( std::pow(parameter(16),parameter(17)) + std::pow(wallData[j][rwI],parameter(17)) ) -
	  parameter(18)*cellData[i][rI]*wallData[j][awI+1];
	wallDerivs[j][rwI+1] -= fac;
	cellDerivs[i][rI] += fac;
      }
    }
  }
}

AuxinROPModel3::
AuxinROPModel3(std::vector<double> &paraValue, 
	      std::vector< std::vector<size_t> > 
	      &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=17 ) {
    std::cerr << "AuxinROPModel3::"
	      << "AuxinROPModel3() "
	      << "17 parameters used (see documentation or network.h)" << std::endl;
    exit(0);
  }
  if( indValue.size() != 2 || indValue[0].size() != 4 || indValue[1].size() != 3 ) {
    std::cerr << "AuxinROPModel3::"
	      << "AuxinROPModel3() "
	      << "Four cell variable indices (auxin, PIN, ROP, AUX, first row) and three wall variable"
	      << " indices (paired) are used (auxin,PIN,ROP)." << std::endl;
    exit(0);
  }
  //Set the variable values
  //
  setId("AuxinROPModel3");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "c_IAA";
  tmp[1] = "d_IAA";
  tmp[2] = "p_IAAH(in)";
  tmp[3] = "p_IAAH(out)";
  tmp[4] = "p_IAA-";
  tmp[5] = "D_IAA";
  tmp[6] = "c_PIN";
  tmp[7] = "d_PIN";
  tmp[8] = "endo_PIN";
  tmp[9] = "exo_PIN";
  tmp[10] = "c_ROP";
  tmp[11] = "d_ROP";
  tmp[12] = "endo_ROP";
  tmp[13] = "exo_ROP";
  tmp[14] = "K_hill";
  tmp[15] = "n_hill";
  tmp[16] = "d_IAA^AUX";
	
  setParameterId( tmp );
}

void AuxinROPModel3::
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
  size_t pI = variableIndex(0,1);//pin
  size_t rI = variableIndex(0,2);//rop
  size_t auxI = variableIndex(0,3);//aux
  size_t awI = variableIndex(1,0);//auxin (wall)
  size_t pwI = variableIndex(1,1);//pin (membrane/wall)
  size_t rwI = variableIndex(1,2);//rop (membrane/wall)

  assert( aI<cellData[0].size() &&
	  pI<cellData[0].size() &&
	  rI<cellData[0].size() &&
	  awI<wallData[0].size() &&
	  pwI<wallData[0].size() &&
	  rwI<wallData[0].size() );

  for (size_t i=0; i<numCells; ++i) {
	  
    //Production and degradation
    cellDerivs[i][aI] += parameter(0) - parameter(1)*cellData[i][aI];
    cellDerivs[i][pI] += parameter(6) - parameter(7)*cellData[i][pI];
    cellDerivs[i][rI] += parameter(10) - parameter(11)*cellData[i][rI];
    
    //Auxin transport and protein cycling
    size_t numWalls = T.cell(i).numWall();
    for (size_t k=0; k<numWalls; ++k) {
      size_t j = T.cell(i).wall(k)->index();
      if( T.cell(i).wall(k)->cell1()->index() == i && T.cell(i).wall(k)->cell2() != T.background() ) {
	// cell-wall transport
	double fac = (parameter(3)+parameter(4)*wallData[j][pwI])*cellData[i][aI] 
	  - (parameter(2)+parameter(16)*cellData[i][auxI])*wallData[j][awI];
	  	
	wallDerivs[j][awI] += fac;
	cellDerivs[i][aI] -= fac;
	// wall-wall diffusion
	wallDerivs[j][awI] -= parameter(5)*wallData[j][awI];
	wallDerivs[j][awI+1] += parameter(5)*wallData[j][awI];
	
	//PIN cycling
	fac = parameter(8)*wallData[j][pwI] - parameter(9)*cellData[i][pI]*wallData[j][rwI];
	wallDerivs[j][pwI] -= fac;
	cellDerivs[i][pI] += fac;

	//ROP cycling
	//fac = parameter(12)*wallData[j][rwI] - parameter(13)*cellData[i][rI]*wallData[j][awI];
	fac = parameter(12)*wallData[j][rwI]
	  *std::pow(wallData[j][rwI+1],parameter(15))/
	  ( std::pow(parameter(14),parameter(15)) + std::pow(wallData[j][rwI+1],parameter(15)) ) -
	  parameter(13)*cellData[i][rI]*wallData[j][awI];
	wallDerivs[j][rwI] -= fac;
	cellDerivs[i][rI] += fac;

      }
      else if( T.cell(i).wall(k)->cell2()->index() == i && T.cell(i).wall(k)->cell1() != T.background() ) {
	// cell-wall transport
	double fac = parameter(3)*cellData[i][aI] - parameter(2)*wallData[j][awI+1] +
	  parameter(4)*cellData[i][aI]*wallData[j][pwI+1];
	
	wallDerivs[j][awI+1] += fac;
	cellDerivs[i][aI] -= fac;
	// wall-wall diffusion
	wallDerivs[j][awI+1] -= parameter(5)*wallData[j][awI+1];
	wallDerivs[j][awI] += parameter(5)*wallData[j][awI+1];

	//PIN cycling
	fac = parameter(8)*wallData[j][pwI+1] - parameter(9)*cellData[i][pI]*wallData[j][rwI+1];
	wallDerivs[j][pwI+1] -= fac;
	cellDerivs[i][pI] += fac;

	//ROP cycling
	//fac = parameter(12)*wallData[j][rwI+1] - parameter(13)*cellData[i][rI]*wallData[j][awI+1];
	fac = parameter(12)*wallData[j][rwI+1]
	  *std::pow(wallData[j][rwI],parameter(15))/
	  ( std::pow(parameter(14),parameter(15)) + std::pow(wallData[j][rwI],parameter(15)) ) -
	  parameter(13)*cellData[i][rI]*wallData[j][awI+1];
	wallDerivs[j][rwI+1] -= fac;
	cellDerivs[i][rI] += fac;

      }
      //else {
      //std::cerr << "AuxinROPModel::derivs() Cell-wall neighborhood wrong." 
      //	  << std::endl;
      //exit(-1);
      //}
    }
  }
}

AuxinPINBistabilityModel::
AuxinPINBistabilityModel(std::vector<double> &paraValue, 
	      std::vector< std::vector<size_t> > 
	      &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=7 ) {
    std::cerr << "AuxinPINBistabilityModel::"
	      << "AuxinPINBistabilityModel() "
	      << "7 parameters used (see network.h)\n";
    exit(0);
  }
  if( indValue.size() != 2 || indValue[0].size() != 2 || indValue[1].size() != 2 ) {
    std::cerr << "AuxinPINBistabilityModel::"
	      << "AuxinPINBistabilityModel() "
	      << "Two cell variable indices (first row) and two wall variable"
	      << " indices are used (auxin,PIN)." << std::endl;
    exit(0);
  }
  //Set the variable values
  //
  setId("AuxinPINBistabilityModel");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "c_IAA";
  tmp[1] = "d_IAA";
  tmp[2] = "D_w";
  tmp[3] = "D_c";
  tmp[4] = "phi";
  tmp[5] = "c_PIN";
  tmp[6] = "d_PIN";
	
  setParameterId( tmp );
}

void AuxinPINBistabilityModel::
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
  size_t pI = variableIndex(0,1);//pin
  size_t awI = variableIndex(1,0);//auxin (wall)
  size_t pwI = variableIndex(1,1);//pin (membrane/wall)

  assert( aI<cellData[0].size() &&
	  pI<cellData[0].size() &&
	  awI<wallData[0].size() &&
	  pwI<wallData[0].size() );


  for (size_t i=0; i<numCells; ++i) {
	  
    //Production and degradation
    cellDerivs[i][aI] += parameter(0) - parameter(1)*cellData[i][aI]; //p_0 - p_1 A_i
    cellDerivs[i][pI] += parameter(5) - parameter(6)*cellData[i][pI]; //p_5 - p_6 P_i
    
    //Auxin transport and protein cycling
    size_t numWalls = T.cell(i).numWall();
    for (size_t k=0; k<numWalls; ++k) {
      size_t j = T.cell(i).wall(k)->index();
      // Checks if cell i is first or second neighbor to wall (wall variables are stored as pairs, and membrane
      // parameters are stored in walls, e.g. wallData[j][pwI] and wallData[j][pwI+1]) 
      if( T.cell(i).wall(k)->cell1()->index() == i && T.cell(i).wall(k)->cell2() != T.background() ) {
	// cell-wall transport
	double fac = parameter(3)*cellData[i][aI] - parameter(2)*wallData[j][awI] +
	  parameter(4)*cellData[i][aI]*wallData[j][pwI]; //p_3 A_i - p_2 A_ij + p_4 A_i P_ij
	
	wallDerivs[j][awI] += fac;
	wallDerivs[j][awI+1] += fac;
	cellDerivs[i][aI] -= fac;
	
	//PIN cycling
	fac = 0.5*(2. - wallData[j][pwI]*wallData[j][pwI+1] + wallData[j][pwI+1]*wallData[j][pwI+1])*wallData[j][pwI] - 
	  cellData[i][pI]*wallData[j][awI];//(1 - P_ij P_ji + P_ji P_ji) P_ij - P_i A_ij
	wallDerivs[j][pwI] -= fac;
	cellDerivs[i][pI] += fac;
      }
      else if( T.cell(i).wall(k)->cell2()->index() == i && T.cell(i).wall(k)->cell1() != T.background() ) {
	// cell-wall transport
	double fac = parameter(3)*cellData[i][aI] - parameter(2)*wallData[j][awI+1] +
	  parameter(4)*cellData[i][aI]*wallData[j][pwI+1];
	
	wallDerivs[j][awI] += fac;
	wallDerivs[j][awI+1] += fac;
	cellDerivs[i][aI] -= fac;
	// wall-wall diffusion
	//wallDerivs[j][awI+1] -= parameter(5)*wallData[j][awI+1];
	//wallDerivs[j][awI] += parameter(5)*wallData[j][awI+1];

	//PIN cycling
	fac = 0.5*(2. - wallData[j][pwI]*wallData[j][pwI+1] + wallData[j][pwI]*wallData[j][pwI])*wallData[j][pwI+1] - 
	  cellData[i][pI]*wallData[j][awI+1];
	wallDerivs[j][pwI+1] -= fac;
	cellDerivs[i][pI] += fac;
      }
      //else {
      //std::cerr << "AuxinPINBistabilityModel::derivs() Cell-wall neighborhood wrong." 
      //	  << std::endl;
      //exit(-1);
      //}
    }
  }
}

AuxinPINBistabilityModelCell::
AuxinPINBistabilityModelCell(std::vector<double> &paraValue, 
			     std::vector< std::vector<size_t> > 
			     &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=8 ) {
    std::cerr << "AuxinPINBistabilityModelCell::"
	      << "AuxinPINBistabilityModelCell() "
	      << "8 parameters used (see network.h)\n";
    exit(0);
  }
  if( indValue.size() != 2 || indValue[0].size() != 2 || indValue[1].size() != 1 ) {
    std::cerr << "AuxinPINBistabilityModelCell::"
	      << "AuxinPINBistabilityModelCell() "
	      << "Two cell variable indices (first row) and one wall variable"
	      << " index are used (PIN)." << std::endl;
    exit(0);
  }
  //Set the variable values
  //
  setId("AuxinPINBistabilityModelCell");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "c_IAA";
  tmp[1] = "d_IAA";
  tmp[2] = "D";
  tmp[3] = "alpha";
  tmp[4] = "phi";
  tmp[5] = "c_PIN";
  tmp[6] = "c_PIN(*Auxin)";
  tmp[7] = "d_PIN (delta)";
	
  setParameterId( tmp );
}

void AuxinPINBistabilityModelCell::
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
  size_t pI = variableIndex(0,1);//pin
  // size_t awI = variableIndex(1,0);//auxin (wall)
  size_t pwI = variableIndex(1,1);//pin (membrane/wall)

  assert( aI<cellData[0].size() &&
	  pI<cellData[0].size() &&
// 	  awI<wallData[0].size() &&
	  pwI<wallData[0].size() );


  for (size_t i=0; i<numCells; ++i) {
	  
    //Production and degradation
    cellDerivs[i][aI] += parameter(0) - parameter(1)*cellData[i][aI]; //p_0 - p_1 A_i
    cellDerivs[i][pI] += parameter(5) + parameter(6)*cellData[i][aI] - parameter(7)*cellData[i][pI]; //p_5 + p_6 A_i - p_7 P_i
    
    //Auxin transport and protein cycling
    size_t numWalls = T.cell(i).numWall();
    for (size_t k=0; k<numWalls; ++k) {
      size_t j = T.cell(i).wall(k)->index();
      // Checks if cell i is first or second neighbor to wall (wall variables are stored as pairs, and membrane
      // parameters are stored in walls, e.g. wallData[j][pwI] and wallData[j][pwI+1]) 
      if( T.cell(i).wall(k)->cell1()->index() == i && T.cell(i).wall(k)->cell2() != T.background() ) {
	// cell-cell transport
	size_t iNeighbor = T.cell(i).wall(k)->cell2()->index();
	double fac = parameter(2)*cellData[i][aI] +
	  parameter(4)*cellData[i][aI]*wallData[j][pwI]; //p_2 A_i + p_4 A_i P_ij
	
	cellDerivs[i][aI] -= fac;
	cellDerivs[iNeighbor][aI] += fac;
	
	//PIN cycling
	fac = 0.5*(2. - wallData[j][pwI]*wallData[j][pwI+1] + wallData[j][pwI+1]*wallData[j][pwI+1]);
	if (fac<=0.) {//to avoid negative concentrations
	  fac = 0.;
	}
	fac = fac*wallData[j][pwI] - parameter(3)*cellData[i][pI];//0.5*(2 - P_ij P_ji + P_ji P_ji) P_ij - alpha P_i
	wallDerivs[j][pwI] -= fac;
	cellDerivs[i][pI] += fac;
      }
      else if( T.cell(i).wall(k)->cell2()->index() == i && T.cell(i).wall(k)->cell1() != T.background() ) {
	// cell-cell transport
	size_t iNeighbor = T.cell(i).wall(k)->cell1()->index();
	double fac = parameter(2)*cellData[i][aI] +
	  parameter(4)*cellData[i][aI]*wallData[j][pwI+1];
	
	cellDerivs[i][aI] -= fac;
	cellDerivs[iNeighbor][aI] += fac;

	//PIN cycling
	fac = 0.5*(2. - wallData[j][pwI]*wallData[j][pwI+1] + wallData[j][pwI]*wallData[j][pwI]);
	if (fac<=0.) {//to avoid negative concentrations
	  fac = 0.;
	}
	fac = fac*wallData[j][pwI+1] - parameter(3)*cellData[i][pI];
	wallDerivs[j][pwI+1] -= fac;
	cellDerivs[i][pI] += fac;
      }
      //else {
      //std::cerr << "AuxinPINBistabilityModelCell::derivs() Cell-wall neighborhood wrong." 
      //	  << std::endl;
      //exit(-1);
      //}
    }
  }
}


AuxinExoBistability::
AuxinExoBistability(std::vector<double> &paraValue, 
			     std::vector< std::vector<size_t> > 
			     &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=9 ) {
    std::cerr << "AuxinExoBistability::"
	      << "AuxinPINBistability() "
	      << "9 parameters used (see network.h)\n";
    exit(0);
  }
  if( indValue.size() != 2 || indValue[0].size() != 2 || indValue[1].size() != 2 ) {
    std::cerr << "AuxinExoBistability::"
	      << "AuxinExoBistability() "
	      << "Two cell variable indices (first row) and one wall variable"
	      << " index are used (PIN)." << std::endl;
    exit(0);
  }
  //Set the variable values
  //
  setId("AuxinExoBistability");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "D_a (auxin diffusion)";
  tmp[1] = "phi (auxin driven PIN transport)";
  tmp[2] = "gamma (PIN cycling rate)";
  tmp[3] = "D_p (PIN diffusion on membrane)";
  tmp[4] = "mu_a (auxin cytosolic regulation)";
  tmp[5] = "mu_p (PIN cytosolic regulation)";
  tmp[6] = "rho_max";
  tmp[7] = "a_X";
  tmp[8] = "rho_B";
  setParameterId( tmp );
}


void AuxinExoBistability::
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
  size_t pI = variableIndex(0,1);//pin
  size_t pwI = variableIndex(1,0);//pin (membrane/wall)
  size_t rwI = variableIndex(1,1);//pin endo rate (membrane/wall)

  assert( aI<cellData[0].size() &&
	  pI<cellData[0].size() &&
	  pwI<wallData[0].size() &&
	  rwI<wallData[0].size() );

  size_t i, j, k, numWalls, iNeighbor, kNext, kBefore, jNext, jBefore;
  double lengthWall, lengthWallNeighbor, g_ij, g_ji, a_ij, cellVolume, fac;	

  for (i=0; i<numCells; ++i) {
    //cytosolic auxin concentration
    cellDerivs[i][aI] += parameter(4)* (1 - cellData[i][aI]); 
    //cytosolic PIN regulation
    cellDerivs[i][pI] += parameter(5)* (cellData[i][aI] - cellData[i][pI]); 
    
    numWalls = T.cell(i).numWall();
    cellVolume = T.cell(i).calculateVolume(vertexData);	
    for (k=0; k<numWalls; ++k) {
      j = T.cell(i).wall(k)->index();
      lengthWall = T.cell(i).wall(k)->length();
      g_ij = lengthWall/cellVolume;
      // Checks if cell i is first or second neighbor to wall (wall variables are stored as pairs)
      if( T.cell(i).wall(k)->cell1()->index() == i && T.cell(i).wall(k)->cell2() != T.background() ) {
	iNeighbor = T.cell(i).wall(k)->cell2()->index();
	g_ji = lengthWall/T.cell(iNeighbor).calculateVolume(vertexData);
	// cell-cell auxin transport
	fac = parameter(0)*cellData[i][aI] +
	  parameter(1)*cellData[i][aI]*wallData[j][pwI]; 
	cellDerivs[i][aI] -= g_ij * fac;
	cellDerivs[iNeighbor][aI] += g_ji * fac;
	//PIN cycling
	fac = parameter(2)*(wallData[j][pwI] - wallData[j][rwI]*cellData[i][pI]);
	wallDerivs[j][pwI] -= fac;
	cellDerivs[i][pI] += g_ij * fac;
        //PIN endocytosis rate
        a_ij = (cellData[i][aI] + cellData[iNeighbor][aI])/2; 
	wallDerivs[j][rwI] += parameter(6)*a_ij /(parameter(7) + a_ij) - wallData[j][rwI]  + 0.5*(-wallData[j][rwI+1] + wallData[j][rwI])*wallData[j][rwI]*wallData[j][rwI + 1]/parameter(8);
       
      }
      else if( T.cell(i).wall(k)->cell2()->index() == i && T.cell(i).wall(k)->cell1() != T.background() ){ 
	iNeighbor = T.cell(i).wall(k)->cell1()->index();
	g_ji = lengthWall/T.cell(iNeighbor).calculateVolume(vertexData);
        // cell-cell auxin transport
	fac = parameter(0)*cellData[i][aI] +
	  parameter(1)*cellData[i][aI]*wallData[j][pwI+1];
	cellDerivs[i][aI] -= g_ij * fac;
	cellDerivs[iNeighbor][aI] +=  g_ji * fac;
	//PIN cycling
	fac = parameter(2)*(wallData[j][pwI+1] - wallData[j][rwI+1] * cellData[i][pI]);
	wallDerivs[j][pwI+1] -= fac;
	cellDerivs[i][pI] += g_ij * fac;
	//PIN endocytosis rate
        a_ij = (cellData[i][aI] + cellData[iNeighbor][aI])/2; 
	wallDerivs[j][rwI+1] += parameter(6)*a_ij/(parameter(7) + a_ij) - wallData[j][rwI+1] + 0.5*(-wallData[j][rwI] + wallData[j][rwI+1])*wallData[j][rwI+1]*wallData[j][rwI]/parameter(8);

      
      }      
      //
      //PIN diffusion on the membrane
      //
      size_t pwIadd=1;//Two variables per wall, setting which to use
      if( T.cell(i).wall(k)->cell1()->index() == i ) {
	pwIadd=0;
      }
      kNext = (k +1 )%numWalls;     
      jNext = T.cell(i).wall(kNext)->index();   
      lengthWallNeighbor = T.cell(i).wall(kNext)->length();
      fac = parameter(3)/(lengthWall* lengthWallNeighbor) * wallData[j][pwI+pwIadd];
      wallDerivs[j][pwI+pwIadd] -= fac;
      if( T.cell(i).wall(kNext)->cell1()->index() == i) {
	wallDerivs[jNext][pwI] +=fac;
      }
      else {
	wallDerivs[jNext][pwI + 1] +=fac;
      }	
      kBefore = k > 0 ? k-1 : numWalls-1; 
      jBefore = T.cell(i).wall(kBefore)->index(); 	
      lengthWallNeighbor = T.cell(i).wall(kBefore)->length();
      fac = parameter(3)/(lengthWall* lengthWallNeighbor) * wallData[j][pwI+pwIadd];
      wallDerivs[j][pwI+pwIadd] -= fac;
      if( T.cell(i).wall(kBefore)->cell1()->index() == i) {
	wallDerivs[jBefore][pwI] +=fac;
      }
      else {
	wallDerivs[jBefore][pwI + 1] +=fac;
      }	
    }
  }
}


SimpleROPModel::
SimpleROPModel(std::vector<double> &paraValue, 
	      std::vector< std::vector<size_t> > 
	      &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=13 ) {
    std::cerr << "SimpleROPModel::"
	      << "SimpleROPModel() "
	      << "13 parameters used (see network.h)\n";
    exit(0);
  }
  if( indValue.size() != 2 || indValue[0].size() != 2 || indValue[1].size() != 2 ) {
    std::cerr << "SimpleROPModel::"
	      << "SimpleROPModel() "
	      << "Two cell variable indices (first row) and two wall variable"
	      << " indices are used (auxin,PIN,ROP)." << std::endl;
    exit(0);
  }
  //Set the variable values
  //
  setId("SimpleROPModel");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "c_IAA";
  tmp[1] = "d_IAA";
  tmp[2] = "p_IAAH(in)";
  tmp[3] = "p_IAAH(out)";
  tmp[4] = "p_IAA-";
  tmp[5] = "D_IAA";
  tmp[6] = "c_PIN";
  tmp[7] = "d_PIN";
  tmp[8] = "endo_PIN";
  tmp[9] = "exo_PIN";
  tmp[10] = "K_hill";
  tmp[11] = "n_hill";
  tmp[12] = "endo_PIN_back";	
  setParameterId( tmp );
}

void SimpleROPModel::
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
  size_t pI = variableIndex(0,1);//pin
  size_t awI = variableIndex(1,0);//auxin (wall)
  size_t pwI = variableIndex(1,1);//pin (membrane/wall)


  assert( aI<cellData[0].size() &&
	  pI<cellData[0].size() &&
	  awI<wallData[0].size() &&
	  pwI<wallData[0].size() );

  for (size_t i=0; i<numCells; ++i) {
	  
    //Production and degradation
    cellDerivs[i][aI] += parameter(0) - parameter(1)*cellData[i][aI];
    cellDerivs[i][pI] += parameter(6) - parameter(7)*cellData[i][pI];

    
    //Auxin transport and protein cycling
    size_t numWalls = T.cell(i).numWall();
    for (size_t k=0; k<numWalls; ++k) {
      size_t j = T.cell(i).wall(k)->index();
      if( T.cell(i).wall(k)->cell1()->index() == i && T.cell(i).wall(k)->cell2() != T.background() ) {
	// cell-wall transport
	double fac = parameter(3)*cellData[i][aI] - parameter(2)*wallData[j][awI] +
	  parameter(4)*cellData[i][aI]*wallData[j][pwI];
	
	wallDerivs[j][awI] += fac;
	cellDerivs[i][aI] -= fac;
	// wall-wall diffusion
	wallDerivs[j][awI] -= parameter(5)*wallData[j][awI];
	wallDerivs[j][awI+1] += parameter(5)*wallData[j][awI];
	
	//PIN cycling
	fac = parameter(12)*wallData[j][pwI]+ parameter(8)*wallData[j][pwI]
	 *std::pow(wallData[j][pwI+1],parameter(11))/
	( std::pow(parameter(10),parameter(11)) + std::pow(wallData[j][pwI+1],parameter(11)) ) -
	parameter(9)*cellData[i][pI]*wallData[j][awI];
	wallDerivs[j][pwI] -= fac;
	cellDerivs[i][pI] += fac;


      }
      else if( T.cell(i).wall(k)->cell2()->index() == i && T.cell(i).wall(k)->cell1() != T.background() ) {
	// cell-wall transport
	double fac = parameter(3)*cellData[i][aI] - parameter(2)*wallData[j][awI+1] +
	  parameter(4)*cellData[i][aI]*wallData[j][pwI+1];
	
	wallDerivs[j][awI+1] += fac;
	cellDerivs[i][aI] -= fac;
	// wall-wall diffusion
	wallDerivs[j][awI+1] -= parameter(5)*wallData[j][awI+1];
	wallDerivs[j][awI] += parameter(5)*wallData[j][awI+1];

	//PIN cycling

	fac = parameter(12)*wallData[j][pwI+1]+ parameter(8)*wallData[j][pwI+1]
	 *std::pow(wallData[j][pwI],parameter(11))/
	( std::pow(parameter(10),parameter(11)) + std::pow(wallData[j][pwI],parameter(11)) ) -
	parameter(9)*cellData[i][pI]*wallData[j][awI+1];

	wallDerivs[j][pwI+1] -= fac;
	cellDerivs[i][pI] += fac;
      }
      //else {
      //std::cerr << "SimpleROPModel::derivs() Cell-wall neighborhood wrong." 
      //	  << std::endl;
      //exit(-1);
      //}
    }
  }
}



SimpleROPModel2::
SimpleROPModel2(std::vector<double> &paraValue, 
	      std::vector< std::vector<size_t> > 
	      &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=13 ) {
    std::cerr << "SimpleROPModel2::"
	      << "SimpleROPModel2() "
	      << "13 parameters used (see network.h)\n";
    exit(0);
  }
  if( indValue.size() != 2 || indValue[0].size() != 2 || indValue[1].size() != 2 ) {
    std::cerr << "SimpleROPModel2::"
	      << "SimpleROPModel2() "
	      << "Two cell variable indices (first row) and two wall variable"
	      << " indices are used (auxin,PIN)." << std::endl;
    exit(0);
  }
  //Set the variable values
  //
  setId("SimpleROPModel2");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "c_IAA";
  tmp[1] = "d_IAA";
  tmp[2] = "p_IAAH(in)";
  tmp[3] = "p_IAAH(out)";
  tmp[4] = "p_IAA-";
  tmp[5] = "D_IAA";
  tmp[6] = "c_PIN";
  tmp[7] = "d_PIN";
  tmp[8] = "endo_PIN";
  tmp[9] = "exo_PIN";
  tmp[10] = "K_hill";
  tmp[11] = "n_hill";
  tmp[12] = "endo_PIN_back";
	
  setParameterId( tmp );
}

void SimpleROPModel2::
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
  size_t pI = variableIndex(0,1);//pin
  size_t awI = variableIndex(1,0);//auxin (wall)
  size_t pwI = variableIndex(1,1);//pin (membrane/wall)


  assert( aI<cellData[0].size() &&
	  pI<cellData[0].size() &&
	  awI<wallData[0].size() &&
	  pwI<wallData[0].size() );

  for (size_t i=0; i<numCells; ++i) {
	  
    //Production and degradation
    cellDerivs[i][aI] += parameter(0) - parameter(1)*cellData[i][aI];
    cellDerivs[i][pI] += parameter(6) - parameter(7)*cellData[i][pI];

    
    //Auxin transport and protein cycling
    size_t numWalls = T.cell(i).numWall();
    for (size_t k=0; k<numWalls; ++k) {
      size_t j = T.cell(i).wall(k)->index();
      if( T.cell(i).wall(k)->cell1()->index() == i && T.cell(i).wall(k)->cell2() != T.background() ) {
	// cell-wall transport
	double fac = parameter(3)*cellData[i][aI] - parameter(2)*wallData[j][awI] +
	  parameter(4)*cellData[i][aI]*wallData[j][pwI];
	
	wallDerivs[j][awI] += fac;
	cellDerivs[i][aI] -= fac;
	// wall-wall diffusion
	wallDerivs[j][awI] -= parameter(5)*wallData[j][awI];
	wallDerivs[j][awI+1] += parameter(5)*wallData[j][awI];
	
	//PIN cycling
	fac =   parameter(12)*wallData[j][pwI]+ parameter(8)*wallData[j][pwI]
	 *std::pow(wallData[j][pwI+1],parameter(11))/
	( std::pow(parameter(10),parameter(11)) + std::pow(wallData[j][pwI+1],parameter(11)) ) -
	parameter(9)*cellData[i][pI];
	wallDerivs[j][pwI] -= fac;
	cellDerivs[i][pI] += fac;
      }
      else if( T.cell(i).wall(k)->cell2()->index() == i && T.cell(i).wall(k)->cell1() != T.background() ) {
	// cell-wall transport
	double fac = parameter(3)*cellData[i][aI] - parameter(2)*wallData[j][awI+1] +
	  parameter(4)*cellData[i][aI]*wallData[j][pwI+1];
	
	wallDerivs[j][awI+1] += fac;
	cellDerivs[i][aI] -= fac;
	// wall-wall diffusion
	wallDerivs[j][awI+1] -= parameter(5)*wallData[j][awI+1];
	wallDerivs[j][awI] += parameter(5)*wallData[j][awI+1];

	//PIN cycling

	fac = parameter(12)*wallData[j][pwI+1]+ parameter(8)*wallData[j][pwI+1]
	 *std::pow(wallData[j][pwI],parameter(11))/
	( std::pow(parameter(10),parameter(11)) + std::pow(wallData[j][pwI],parameter(11)) ) -
	parameter(9)*cellData[i][pI];
	wallDerivs[j][pwI+1] -= fac;
	cellDerivs[i][pI] += fac;
      }
      //else {
      //std::cerr << "SimpleROPModel2::derivs() Cell-wall neighborhood wrong." 
      //	  << std::endl;
      //exit(-1);
      //}
    }
  }
}


SimpleROPModel3::
SimpleROPModel3(std::vector<double> &paraValue, 
	      std::vector< std::vector<size_t> > 
	      &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=13 ) {
    std::cerr << "SimpleROPModel3::"
	      << "SimpleROPModel3() "
	      << "13 parameters used (see network.h)\n";
    exit(0);
  }
  if( indValue.size() != 2 || indValue[0].size() != 2 || indValue[1].size() != 2 ) {
    std::cerr << "SimpleROPModel3::"
	      << "SimpleROPModel3() "
	      << "Two cell variable indices (first row) and two wall variable"
	      << " indices are used (auxin,PIN)." << std::endl;
    exit(0);
  }
  //Set the variable values
  //
  setId("SimpleROPModel3");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "c_IAA";
  tmp[1] = "d_IAA";
  tmp[2] = "p_IAAH(in)";
  tmp[3] = "p_IAAH(out)";
  tmp[4] = "p_IAA-";
  tmp[5] = "D_IAA";
  tmp[6] = "c_PIN";
  tmp[7] = "d_PIN";
  tmp[8] = "endo_PIN";
  tmp[9] = "exo_PIN";
  tmp[10] = "K_hill";
  tmp[11] = "n_hill";
	 tmp[12] = "endo_back";
  setParameterId( tmp );
}

void SimpleROPModel3::
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
  size_t pI = variableIndex(0,1);//pin
  size_t awI = variableIndex(1,0);//auxin (wall)
  size_t pwI = variableIndex(1,1);//pin (membrane/wall)


  assert( aI<cellData[0].size() &&
	  pI<cellData[0].size() &&
	  awI<wallData[0].size() &&
	  pwI<wallData[0].size() );

 for (size_t i=0; i<numCells; ++i) {
	  
    //Production and degradation
    cellDerivs[i][aI] += parameter(0) - parameter(1)*cellData[i][aI];
    cellDerivs[i][pI] += parameter(6) - parameter(7)*cellData[i][pI];

    
    //Auxin transport and protein cycling
    size_t numWalls = T.cell(i).numWall();
    for (size_t k=0; k<numWalls; ++k) {
      size_t j = T.cell(i).wall(k)->index();
      if( T.cell(i).wall(k)->cell1()->index() == i && T.cell(i).wall(k)->cell2() != T.background() ) {
	// cell-wall transport
	double fac = parameter(3)*cellData[i][aI] - parameter(2)*wallData[j][awI] +
	  parameter(4)*cellData[i][aI]*wallData[j][pwI];
	
	wallDerivs[j][awI] += fac;
	cellDerivs[i][aI] -= fac;
	// wall-wall diffusion
	wallDerivs[j][awI] -= parameter(5)*wallData[j][awI];
	wallDerivs[j][awI+1] += parameter(5)*wallData[j][awI];
	
	//PIN cycling
	fac = parameter(8)*wallData[j][pwI]*
	   std::pow(wallData[j][pwI+1],parameter(11))  /
	  	( std::pow(parameter(10),parameter(11)) + std::pow(wallData[j][pwI+1],parameter(11)) ) -
	  parameter(9)*cellData[i][pI]+parameter(12)*wallData[j][pwI];
	wallDerivs[j][pwI] -= fac;
	cellDerivs[i][pI] += fac;


      }
      else if( T.cell(i).wall(k)->cell2()->index() == i && T.cell(i).wall(k)->cell1() != T.background() ) {
	// cell-wall transport
	double fac = parameter(3)*cellData[i][aI] - parameter(2)*wallData[j][awI+1] +
	  parameter(4)*cellData[i][aI]*wallData[j][pwI+1];
	
	wallDerivs[j][awI+1] += fac;
	cellDerivs[i][aI] -= fac;
	// wall-wall diffusion
	wallDerivs[j][awI+1] -= parameter(5)*wallData[j][awI+1];
	wallDerivs[j][awI] += parameter(5)*wallData[j][awI+1];

	//PIN cycling

	fac = parameter(8)*std::pow(wallData[j][pwI],parameter(11))  /
	  ( std::pow(parameter(10),parameter(11)) + std::pow(wallData[j][pwI],parameter(11)) )*  wallData[j][pwI+1]
	 -
	parameter(9)*
cellData[i][pI]+parameter(12)*wallData[j][pwI+1];

	wallDerivs[j][pwI+1] -= fac;
	cellDerivs[i][pI] += fac;
      }
      //else {
      //std::cerr << "SimpleROPModel::derivs() Cell-wall neighborhood wrong." 
      //	  << std::endl;
      //exit(-1);
      //}
    }
  }

}




SimpleROPModel4::
SimpleROPModel4(std::vector<double> &paraValue, 
	      std::vector< std::vector<size_t> > 
	      &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=13 ) {
    std::cerr << "SimpleROPModel4::"
	      << "SimpleROPModel4() "
	      << "13 parameters used (see network.h)\n";
    exit(0);
  }
  if( indValue.size() != 2 || indValue[0].size() != 2 || indValue[1].size() != 2) {
    std::cerr << "SimpleROPModel4::"
	      << "SimpleROPModel4() "
	      << "Two cell variable indices (first row) and two wall variable"
	      << " indices are used (auxin,PIN)." << std::endl;
    exit(0);
  }
  //Set the variable values
  //
  setId("SimpleROPModel4");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "c_IAA";
  tmp[1] = "d_IAA";
  tmp[2] = "p_IAAH(in)";
  tmp[3] = "p_IAAH(out)";
  tmp[4] = "p_IAA-";
  tmp[5] = "D_IAA";
  tmp[6] = "c_PIN";
  tmp[7] = "d_PIN";
  tmp[8] = "endo_PIN";
  tmp[9] = "exo_PIN";
  tmp[10] = "K_hill";
  tmp[11] = "n_hill";
  tmp[12] = "endo_back";
  setParameterId( tmp );
}

void SimpleROPModel4::
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
  size_t pI = variableIndex(0,1);//pin
  size_t awI = variableIndex(1,0);//auxin (wall)
  size_t pwI = variableIndex(1,1);//pin (membrane/wall)


  assert( aI<cellData[0].size() &&
	  pI<cellData[0].size() &&
	  awI<wallData[0].size() &&
	  pwI<wallData[0].size() );

 for (size_t i=0; i<numCells; ++i) {
	  
    //Production and degradation
    cellDerivs[i][aI] += parameter(0) - parameter(1)*cellData[i][aI];
    cellDerivs[i][pI] += parameter(6)*cellData[i][aI]/(parameter(4)+cellData[i][aI]) - parameter(7)*cellData[i][pI];

    
    //Auxin transport and protein cycling
    size_t numWalls = T.cell(i).numWall();
    for (size_t k=0; k<numWalls; ++k) {
      size_t j = T.cell(i).wall(k)->index();
      if( T.cell(i).wall(k)->cell1()->index() == i && T.cell(i).wall(k)->cell2() != T.background() ) {
	// cell-cell transport
	size_t iNeighbor = T.cell(i).wall(k)->cell2()->index();
	double fac = parameter(2)*cellData[i][aI] +
	  parameter(3)*cellData[i][aI]*wallData[j][pwI]; //p_2 A_i + p_4 A_i P_ij
	
	cellDerivs[i][aI] -= fac;
	cellDerivs[iNeighbor][aI] += fac;
	
	
	//PIN cycling
	fac = parameter(8)*wallData[j][pwI]*std::pow(wallData[j][pwI+1],parameter(11))/( std::pow(parameter(10),parameter(11)) + std::pow(wallData[j][pwI+1],parameter(11)) ) 
              -	  parameter(9)*cellData[i][pI]+parameter(12)*wallData[j][pwI];
	wallDerivs[j][pwI] -= fac;
	cellDerivs[i][pI] += fac;


      }
      else if( T.cell(i).wall(k)->cell2()->index() == i && T.cell(i).wall(k)->cell1() != T.background()) {
	// cell-cell transport
	size_t iNeighbor = T.cell(i).wall(k)->cell1()->index();
	double fac = parameter(2)*cellData[i][aI] +
	  parameter(3)*cellData[i][aI]*wallData[j][pwI+1];
	
	cellDerivs[i][aI] -= fac;
	cellDerivs[iNeighbor][aI] += fac;


	//PIN cycling

	fac = parameter(8)*wallData[j][pwI+1]*std::pow(wallData[j][pwI],parameter(11))/( std::pow(parameter(10),parameter(11)) + std::pow(wallData[j][pwI],parameter(11)) ) 
	 -parameter(9)*cellData[i][pI]+parameter(12)*wallData[j][pwI+1];

	wallDerivs[j][pwI+1] -= fac;
	cellDerivs[i][pI] += fac;
      }
    }
  }
}



SimpleROPModel5::
SimpleROPModel5(std::vector<double> &paraValue, 
	      std::vector< std::vector<size_t> > 
	      &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=16 ) {
    std::cerr << "SimpleROPModel4::"
	      << "SimpleROPModel4() "
	      << "16 parameters used (see network.h)\n";
    exit(0);
  }
  if( indValue.size() != 2 || indValue[0].size() != 2 || indValue[1].size() != 2 ) {
    std::cerr << "SimpleROPModel4::"
	      << "SimpleROPModel4() "
	      << "Two cell variable indices (first row) and two wall variable"
	      << " indices are used (auxin,PIN)." << std::endl;
    exit(0);
  }
  //Set the variable values
  //
  setId("SimpleROPModel5");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "c_IAA";
  tmp[1] = "d_IAA";
  tmp[2] = "p_IAAH(in)";
  tmp[3] = "p_IAAH(out)";
  tmp[4] = "p_IAA-";
  tmp[5] = "D_IAA";
  tmp[6] = "c_PIN";
  tmp[7] = "d_PIN";
  tmp[8] = "endo_PIN";
  tmp[9] = "exo_PIN";
  tmp[10] = "K_hill";
  tmp[11] = "n_hill";
	 tmp[12] = "endo_back";
  tmp[13] = "K_hill";
  tmp[14] = "n_hill";
	 tmp[15] = "endo_back";
  setParameterId( tmp );
}

void SimpleROPModel5::
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
  size_t pI = variableIndex(0,1);//pin
  size_t awI = variableIndex(1,0);//auxin (wall)
  size_t pwI = variableIndex(1,1);//pin (membrane/wall)


  assert( aI<cellData[0].size() &&
	  pI<cellData[0].size() &&
	  awI<wallData[0].size() &&
	  pwI<wallData[0].size() );

 for (size_t i=0; i<numCells; ++i) {
	  
    //Production and degradation
    cellDerivs[i][aI] += parameter(0) - parameter(1)*cellData[i][aI];
    cellDerivs[i][pI] += parameter(6)*cellData[i][aI]/(parameter(4)+cellData[i][aI]) - parameter(7)*cellData[i][pI];

    
    //Auxin transport and protein cycling
    size_t numWalls = T.cell(i).numWall();
    for (size_t k=0; k<numWalls; ++k) {
      size_t j = T.cell(i).wall(k)->index();
      if( T.cell(i).wall(k)->cell1()->index() == i && T.cell(i).wall(k)->cell2() != T.background() ) {
	// cell-cell transport
	size_t iNeighbor = T.cell(i).wall(k)->cell2()->index();
	double fac = parameter(2)*cellData[i][aI] +
	  parameter(3)*cellData[i][aI]*wallData[j][pwI]; //p_2 A_i + p_4 A_i P_ij
	
	cellDerivs[i][aI] -= fac;
	cellDerivs[iNeighbor][aI] += fac;
	
	//PIN cycling
	fac = parameter(8)*wallData[j][pwI]*
	   std::pow(wallData[j][pwI+1],parameter(11))  /
	  	( std::pow(parameter(10),parameter(11)) + std::pow(wallData[j][pwI+1],parameter(11)) ) -  parameter(9)*cellData[i][pI]+parameter(12)*wallData[j][pwI]-parameter(13)*cellData[i][pI]*std::pow(wallData[j][pwI],parameter(14) ) /
	  ( std::pow(parameter(15),parameter(14)) + std::pow(wallData[j][pwI],parameter(14)) );;
	wallDerivs[j][pwI] -= fac;
	cellDerivs[i][pI] += fac;


      }
      else if( T.cell(i).wall(k)->cell2()->index() == i && T.cell(i).wall(k)->cell1() != T.background() ) {
	// cell-cell transport
	size_t iNeighbor = T.cell(i).wall(k)->cell1()->index();
	double fac = parameter(2)*cellData[i][aI] +
	  parameter(3)*cellData[i][aI]*wallData[j][pwI+1];
	
	cellDerivs[i][aI] -= fac;
	cellDerivs[iNeighbor][aI] += fac;


	//PIN cycling

	fac = parameter(8)*std::pow(wallData[j][pwI],parameter(11))  /
	  ( std::pow(parameter(10),parameter(11)) + std::pow(wallData[j][pwI],parameter(11)) )*  wallData[j][pwI+1] -parameter(9)*cellData[i][pI]+parameter(12)*wallData[j][pwI+1]-
parameter(13)*cellData[i][pI]*std::pow(wallData[j][pwI+1],parameter(14) ) /
	  ( std::pow(parameter(15),parameter(14)) + std::pow(wallData[j][pwI+1],parameter(14)) );

	wallDerivs[j][pwI+1] -= fac;
	cellDerivs[i][pI] += fac;
      }
    }
  }
}






SimpleROPModel6::
SimpleROPModel6(std::vector<double> &paraValue, 
	      std::vector< std::vector<size_t> > 
	      &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=13 ) {
    std::cerr << "SimpleROPModel4::"
	      << "SimpleROPModel4() "
	      << "13 parameters used (see network.h)\n";
    exit(0);
  }
  if( indValue.size() != 2 || indValue[0].size() != 2 || indValue[1].size() != 2 ) {
    std::cerr << "SimpleROPModel4::"
	      << "SimpleROPModel4() "
	      << "Two cell variable indices (first row) and two wall variable"
	      << " indices are used (auxin,PIN)." << std::endl;
    exit(0);
  }
  //Set the variable values
  //
  setId("SimpleROPModel6");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "c_IAA";
  tmp[1] = "d_IAA";
  tmp[2] = "p_IAAH(in)";
  tmp[3] = "p_IAAH(out)";
  tmp[4] = "p_IAA-";
  tmp[5] = "D_IAA";
  tmp[6] = "c_PIN";
  tmp[7] = "d_PIN";
  tmp[8] = "endo_PIN";
  tmp[9] = "exo_PIN";
  tmp[10] = "K_hill";
  tmp[11] = "n_hill";
  tmp[12] = "endo_back";
 
  setParameterId( tmp );
}

void SimpleROPModel6::
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
  size_t pI = variableIndex(0,1);//pin
  size_t awI = variableIndex(1,0);//auxin (wall)
  size_t pwI = variableIndex(1,1);//pin (membrane/wall)


  assert( aI<cellData[0].size() &&
	  pI<cellData[0].size() &&
	  awI<wallData[0].size() &&
	  pwI<wallData[0].size() );

 for (size_t i=0; i<numCells; ++i) {
	  
    //Production and degradation
    cellDerivs[i][aI] += parameter(0) - parameter(1)*cellData[i][aI];
    cellDerivs[i][pI] += parameter(6)*cellData[i][aI]/(parameter(4)+cellData[i][aI]) - parameter(7)*cellData[i][pI];

    
    //Auxin transport and protein cycling
    size_t numWalls = T.cell(i).numWall();
    for (size_t k=0; k<numWalls; ++k) {


      size_t j = T.cell(i).wall(k)->index();
      if( T.cell(i).wall(k)->cell1()->index() == i && T.cell(i).wall(k)->cell2() != T.background() ) {
	// cell-cell transport
	size_t iNeighbor = T.cell(i).wall(k)->cell2()->index();
	double fac = parameter(2)*cellData[i][aI] +
	  parameter(3)*cellData[i][aI]*wallData[j][pwI]; //p_2 A_i + p_4 A_i P_ij
	
	cellDerivs[i][aI] -= fac;
	cellDerivs[iNeighbor][aI] += fac;
	
	//PIN cycling
	fac = parameter(8)*wallData[j][pwI]*
	   std::pow(wallData[j][pwI+1],parameter(11))  /
	  ( std::pow(parameter(10),parameter(11)) + std::pow(wallData[j][pwI+1],parameter(11)) ) -  parameter(9)*cellData[i][pI]+parameter(12)*wallData[j][pwI];
	wallDerivs[j][pwI] -= fac;
	cellDerivs[i][pI] += fac;


      }
      else if( T.cell(i).wall(k)->cell2()->index() == i && T.cell(i).wall(k)->cell1() != T.background() ) {
	// cell-cell transport
	size_t iNeighbor = T.cell(i).wall(k)->cell1()->index();
	double fac = parameter(2)*cellData[i][aI] +
	  parameter(3)*cellData[i][aI]*wallData[j][pwI+1];
	
	cellDerivs[i][aI] -= fac;
	cellDerivs[iNeighbor][aI] += fac;


	//PIN cycling

	fac = parameter(8)*std::pow(wallData[j][pwI],parameter(11))  /
	  ( std::pow(parameter(10),parameter(11)) + std::pow(wallData[j][pwI],parameter(11)) )*  wallData[j][pwI+1] -parameter(9)*cellData[i][pI]+parameter(12)*wallData[j][pwI+1];
	wallDerivs[j][pwI+1] -= fac;
	cellDerivs[i][pI] += fac;
      }
    else{
	//PIN cycling

	double fac =  -parameter(9)*cellData[i][pI]+parameter(12)*wallData[j][pwI];
	wallDerivs[j][pwI] -= fac;
	cellDerivs[i][pI] += fac;
      }
    }
  }

}


SimpleROPModel7::
SimpleROPModel7(std::vector<double> &paraValue, 
	      std::vector< std::vector<size_t> > 
	      &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=15 ) {
    std::cerr << "SimpleROPModel7::"
	      << "SimpleROPModel7() "
	      << "13 parameters used (see network.h)\n";
    exit(0);
  }
  if( indValue.size() != 2 || indValue[0].size() != 2 || indValue[1].size() != 3 ) {
    std::cerr << "SimpleROPModel7::"
	      << "SimpleROPModel7() "
	      << "Two cell variable indices (first row) and two wall variable"
	      << " indices are used (auxin,PIN)." << std::endl;
    exit(0);
  }
  //Set the variable values
  //
  setId("SimpleROPModel7");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "c_IAA";
  tmp[1] = "d_IAA";
  tmp[2] = "p_IAAH(in)";
  tmp[3] = "p_IAAH(out)";
  tmp[4] = "p_IAA-";
  tmp[5] = "D_IAA";
  tmp[6] = "c_PIN";
  tmp[7] = "d_PIN";
  tmp[8] = "endo_PIN";
  tmp[9] = "exo_PIN";
  tmp[10] = "K_hill";
  tmp[11] = "n_hill";
  tmp[12] = "endo_back";
 tmp[13] = "n_hill";
  tmp[14] = "endo_back";
  setParameterId( tmp );
}

void SimpleROPModel7::
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
  size_t pI = variableIndex(0,1);//pin
  size_t awI = variableIndex(1,0);//auxin (wall)
  size_t pwI = variableIndex(1,1);//pin (membrane/wall)
 size_t mwI = variableIndex(1,2);//pin (membrane/wall)

  assert( aI<cellData[0].size() &&
	  pI<cellData[0].size() &&
	  awI<wallData[0].size() &&
	  pwI<wallData[0].size()  &&
	  mwI<wallData[0].size() );

 for (size_t i=0; i<numCells; ++i) {
	  
    //Production and degradation
    cellDerivs[i][aI] += parameter(0) - parameter(1)*cellData[i][aI];
    cellDerivs[i][pI] += parameter(6)*cellData[i][aI]/(parameter(4)+cellData[i][aI]) - parameter(7)*cellData[i][pI];

    
    //Auxin transport and protein cycling
    size_t numWalls = T.cell(i).numWall();
    for (size_t k=0; k<numWalls; ++k) {
      size_t j = T.cell(i).wall(k)->index();
      if( T.cell(i).wall(k)->cell1()->index() == i && T.cell(i).wall(k)->cell2() != T.background()  && wallData[j][mwI]==1) {
	// cell-cell transport
	size_t iNeighbor = T.cell(i).wall(k)->cell2()->index();
	double fac = parameter(2)*cellData[i][aI] +
	  parameter(3)*cellData[i][aI]*wallData[j][pwI]; //p_2 A_i + p_4 A_i P_ij
	
	cellDerivs[i][aI] -= fac;
	cellDerivs[iNeighbor][aI] += fac;
	
	
	//PIN cycling
	fac = parameter(8)*wallData[j][pwI]*std::pow(wallData[j][pwI+1],parameter(11))/( std::pow(parameter(10),parameter(11)) + std::pow(wallData[j][pwI+1],parameter(11)) ) 
              -	  parameter(9)*cellData[i][pI]+parameter(12)*wallData[j][pwI];
	wallDerivs[j][pwI] -= fac;
	cellDerivs[i][pI] += fac;


      }
      else if( T.cell(i).wall(k)->cell2()->index() == i && T.cell(i).wall(k)->cell1() != T.background()  && wallData[j][mwI]==1) {
	// cell-cell transport
	size_t iNeighbor = T.cell(i).wall(k)->cell1()->index();
	double fac = parameter(2)*cellData[i][aI] +
	  parameter(3)*cellData[i][aI]*wallData[j][pwI+1];
	
	cellDerivs[i][aI] -= fac;
	cellDerivs[iNeighbor][aI] += fac;


	//PIN cycling

	fac = parameter(8)*wallData[j][pwI+1]*std::pow(wallData[j][pwI],parameter(11))/( std::pow(parameter(10),parameter(11)) + std::pow(wallData[j][pwI],parameter(11)) ) 
	 -parameter(9)*cellData[i][pI]+parameter(12)*wallData[j][pwI+1];

	wallDerivs[j][pwI+1] -= fac;
	cellDerivs[i][pI] += fac;
      }

  else if( T.cell(i).wall(k)->cell2()->index() == i && T.cell(i).wall(k)->cell1() != T.background()  && wallData[j][mwI]==0) {
	// cell-cell transport
	size_t iNeighbor = T.cell(i).wall(k)->cell1()->index();
	double fac = parameter(13)*cellData[i][aI];
		double fac2 = parameter(14)*cellData[i][pI];
	cellDerivs[i][aI] -= fac;
        cellDerivs[i][pI] -= fac2;
	cellDerivs[iNeighbor][aI] += fac;
	cellDerivs[iNeighbor][pI] += fac2;
      }

  else if( T.cell(i).wall(k)->cell1()->index() == i && T.cell(i).wall(k)->cell2() != T.background()  && wallData[j][mwI]==0) {
	// cell-cell transport
	size_t iNeighbor = T.cell(i).wall(k)->cell1()->index();
	double fac = parameter(13)*cellData[i][aI];
	double fac2 = parameter(14)*cellData[i][pI];
	
	cellDerivs[i][aI] -= fac;
         cellDerivs[i][pI] -= fac2;
	cellDerivs[iNeighbor][aI] += fac;
cellDerivs[iNeighbor][pI] += fac2;

      }
    }
  }
}


AuxinModelSimple4::
AuxinModelSimple4(std::vector<double> &paraValue, 
		  std::vector< std::vector<size_t> > 
		  &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=12 ) {
    std::cerr << "AuxinModelSimple4::"
	      << "AuxinModelSimple4() "
	      << "Twelve parameters used.\n\n";
    std::cerr << "dA_i/dt = p0*M_i + p1 - p2*A_i +p5*Sum_{neigh} (A_n-A_i) +\n" 
	      << "p4*Sum_{neigh} (P_ni*A_n-P_in*A_i)\n\n" 
	      << "dP_i/dt = p6 - p7*P_i\n\n"
	      << "dX_i/dt = p8*A_i - p9*X_i\n\n"
	      << "dM_i/dt = p10*Theta_L1 - p11*M_i\n\n"
	      << "P_in = P_i*X_n/(p_3+Sum_{k,neigh}X_k)\n";
    exit(EXIT_FAILURE);
  }
  if( (indValue.size() != 1 && indValue.size() !=2) || indValue[0].size() != 4 ||
      (indValue.size()==2 &&indValue[1].size() != 1) ) {
    std::cerr << "AuxinModelSimple4::"
	      << "AuxinModelSimple4() "
	      << "Four variable indices are used (auxin,pin,X,M) in first level."
	      << std::endl << "Optionally a wall index is given for saving PIN"
	      << " at membranes." << std::endl;
    exit(EXIT_FAILURE);
  }
  //Set the variable values
  //
  setId("AuxinModelSimple4");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "p_auxin(M)";
  tmp[1] = "p_auxin";
  tmp[2] = "d_auxin";
  tmp[3] = "p_pol";
  tmp[4] = "T_auxin";
  tmp[5] = "D_auxin";
  tmp[6] = "p_pin";
  tmp[7] = "d_pin";
  tmp[8] = "p_X";
  tmp[9] = "d_X";
  tmp[10] = "p_M";
  tmp[11] = "d_M";
  
  setParameterId( tmp );
}

void AuxinModelSimple4::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{  
  size_t numCells = T.numCell();
  size_t aI = variableIndex(0,0);
  size_t pI = variableIndex(0,1);
  size_t xI = variableIndex(0,2);
  size_t mI = variableIndex(0,3);
  assert( aI<cellData[0].size() &&
	  pI<cellData[0].size() &&
	  xI<cellData[0].size() &&
	  mI<cellData[0].size() );
  
  if (numVariableIndexLevel()==2) {
    assert(variableIndex(1,0)<wallData[0].size());
  }

  for( size_t i=0 ; i<numCells ; ++i ) {
    
    //Production and degradation
    cellDerivs[i][aI] += parameter(0)*cellData[i][mI] + parameter(1) - 
      parameter(2)*cellData[i][aI];
    
    cellDerivs[i][pI] += parameter(6) - parameter(7)*cellData[i][pI];
    
    cellDerivs[i][xI] += parameter(8)*cellData[i][aI] 
      - parameter(9)*cellData[i][xI];
    
    cellDerivs[i][mI] -= parameter(11)*cellData[i][mI];
    if( T.cell(i).isNeighbor(T.background()) )
      cellDerivs[i][mI] += parameter(10);
    
    //Transport
    size_t numWalls=T.cell(i).numWall();
    //Polarization coefficient normalization constant
    double sum=0.0;
    size_t numActualWalls=0;
    //pin[i].resize( numWalls+1 );
    for( size_t n=0 ; n<numWalls ; ++n ) {
      if( T.cell(i).wall(n)->cell1() != T.background() &&
	  T.cell(i).wall(n)->cell2() != T.background() ) { 
	numActualWalls++;
	if( T.cell(i).wall(n)->cell1()->index()==i )
	  sum += cellData[ T.cell(i).wall(n)->cell2()->index() ][ xI ];
	else
	  sum += cellData[ T.cell(i).wall(n)->cell1()->index() ][ xI ];
      }
    }
    //sum /= numActualWalls;//For adjusting for different num neigh
    sum += parameter(3);
    
    for( size_t n=0 ; n<numWalls ; ++n ) {

      // transport from cell i, wall n to neighbor

      //if( !T.cell(i).isNeighbor(T.background()) ) { 
      if( T.cell(i).wall(n)->cell1() != T.background() &&
	  T.cell(i).wall(n)->cell2() != T.background() ) { 
	size_t neighIndex;
	size_t pinIndexAdd=0; 

	if( T.cell(i).wall(n)->cell1()->index()==i )
	  // if cell1 of the wall is the current cell, cell2 is the neighbor
	  neighIndex = T.cell(i).wall(n)->cell2()->index();				
	else {
	  // if cell1 of the wall is not the current cell, cell1 must be the neighbor
	  pinIndexAdd=1;
	  neighIndex = T.cell(i).wall(n)->cell1()->index();				
	}
	double polRate=0.0;
	
	if( sum != 0.0 ) {
	  polRate = cellData[i][pI] * cellData[neighIndex][xI] / sum;
	}
	else {
	  polRate = 1.;
	}

	// TODO: save PIN (each wall compartment only has one PIN concentration...)

	if (numVariableIndexLevel()==2) {//store PIN value for membrane (in wall data)
	  wallData[T.cell(i).wall(n)->index()][variableIndex(1,0) + pinIndexAdd] = polRate;
	  //wallData[T.cell(i).wall(n)->index()][variableIndex(1,0)+pinIndexAdd] = 1.0;
	}
	//pin[i][n+1] = polRate;
	cellDerivs[i][aI] -= (parameter(4)*polRate+parameter(5))*cellData[i][aI];
	cellDerivs[neighIndex][aI] += (parameter(4)*polRate+parameter(5))*cellData[i][aI];
      }
    }
  }
}

AuxinModelSimple5::
AuxinModelSimple5(std::vector<double> &paraValue, 
		  std::vector< std::vector<size_t> > 
		  &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=10 ) {
    std::cerr << "AuxinModelSimple5::"
	      << "AuxinModelSimple5() "
	      << "Ten parameters used.\n\n";
    std::cerr << "dA_i/dt = p1 - p2*A_i +p5*Sum_{neigh} (A_n-A_i) +\n" 
	      << "p4*Sum_{neigh} (P_ni*A_n-P_in*A_i)\n\n" 
	      << "dP_i/dt = p6 - p7*P_i\n\n"
	      << "dX_i/dt = p8*A_i - p9*X_i\n\n"
	      << "P_in = P_i*(X_n + p0) /(p_3+Sum_{k,neigh}( X_k + p_0))\n";
    exit(EXIT_FAILURE);
  }
  if( (indValue.size() != 1 && indValue.size() !=2) || indValue[0].size() != 3 ||
      (indValue.size()==2 &&indValue[1].size() != 1) ) {
    std::cerr << "AuxinModelSimple1::"
	      << "AuxinModelSimple1() "
	      << "Three variable indices are used (auxin,pin,X) in first level."
	      << std::endl << "Optionally a wall index is given for saving (pair of) PIN"
	      << " at membranes." << std::endl;
    exit(EXIT_FAILURE);
  }
  //Set the variable values
  //
  setId("AuxinModelSimple5");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "c_pinexo";
  tmp[1] = "p_auxin";
  tmp[2] = "d_auxin";
  tmp[3] = "p_pol";
  tmp[4] = "T_auxin";
  tmp[5] = "D_auxin";
  tmp[6] = "p_pin";
  tmp[7] = "d_pin";
  tmp[8] = "p_X";
  tmp[9] = "d_X";
  
  setParameterId( tmp );
}

void AuxinModelSimple5::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{  
  size_t numCells = T.numCell();
  size_t aI = variableIndex(0,0);
  size_t pI = variableIndex(0,1);
  size_t xI = variableIndex(0,2);
  assert( aI<cellData[0].size() &&
	  pI<cellData[0].size() &&
	  xI<cellData[0].size() &&);
  
  if (numVariableIndexLevel()==2) {
    assert(variableIndex(1,0)<wallData[0].size());
  }

  for( size_t i=0 ; i<numCells ; ++i ) {
    
    //Production and degradation
    cellDerivs[i][aI] += parameter(1) - parameter(2)*cellData[i][aI];
    
    cellDerivs[i][pI] += parameter(6) - parameter(7)*cellData[i][pI];
    
    cellDerivs[i][xI] += parameter(8)*cellData[i][aI] 
      - parameter(9)*cellData[i][xI];
    
    //Transport
    size_t numWalls=T.cell(i).numWall();
    //Polarization coefficient normalization constant
    double sum=0.0;
    size_t numActualWalls=0;
    //pin[i].resize( numWalls+1 );
    for( size_t n=0 ; n<numWalls ; ++n ) {
      if( T.cell(i).wall(n)->cell1() != T.background() &&
	  T.cell(i).wall(n)->cell2() != T.background() ) { 
	numActualWalls++;
	if( T.cell(i).wall(n)->cell1()->index()==i )
	  sum += cellData[ T.cell(i).wall(n)->cell2()->index() ][ xI ] + parameter(0); // added constant exo
	else
	  sum += cellData[ T.cell(i).wall(n)->cell1()->index() ][ xI ]+ parameter(0);   // added constant exo
      }
    }
    //sum /= numActualWalls;//For adjusting for different num neigh
    sum += parameter(3);
    
    for( size_t n=0 ; n<numWalls ; ++n ) {
      //if( !T.cell(i).isNeighbor(T.background()) ) { 
      if( T.cell(i).wall(n)->cell1() != T.background() &&
	  T.cell(i).wall(n)->cell2() != T.background() ) { 
	size_t neighIndex;
	size_t pinIndexAdd=0; 
	if( T.cell(i).wall(n)->cell1()->index()==i )
	  neighIndex = T.cell(i).wall(n)->cell2()->index();				
	else {
	  neighIndex = T.cell(i).wall(n)->cell1()->index();				
	  pinIndexAdd=1;
	}
	double polRate=0.0;
	
	if( sum != 0.0 ) {
	  polRate = cellData[i][pI] * (cellData[neighIndex][xI] + + parameter(0)) / sum;  // added constant exo
	}
	else {
	  polRate = 1.;
	}
	if (numVariableIndexLevel()==2) {//store PIN value for membrane (in wall pair data)
	  wallData[T.cell(i).wall(n)->index()][variableIndex(1,0)+pinIndexAdd] = polRate;
	}
	//pin[i][n+1] = polRate;
	cellDerivs[i][aI] -= (parameter(4)*polRate+parameter(5))*cellData[i][aI];
	cellDerivs[neighIndex][aI] += (parameter(4)*polRate+parameter(5))*cellData[i][aI];
      }
    }
  }
}




UpInternalGradientModel::
UpInternalGradientModel(std::vector<double> &paraValue, 
	      std::vector< std::vector<size_t> > 
	      &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=15 ) {
    std::cerr << "UpInternalGradientModel::"
	      << "UpInternalGradientModel() "
	      << "13 parameters used (see network.h)\n";
    exit(0);
  }
  if( indValue.size() != 2 || indValue[0].size() != 2 || indValue[1].size() != 3 ) {
    std::cerr << "UpInternalGradientModel::"
	      << "UpInternalGradientModel() "
	      << "Two cell variable indices (first row) and two wall variable"
	      << " indices are used (auxin,PIN)." << std::endl;
    exit(0);
  }
  //Set the variable values
  //
  setId("UpInternalGradientModel");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "c_IAA";
  tmp[1] = "d_IAA";
  tmp[2] = "p_IAAH(in)";
  tmp[3] = "p_IAAH(out)";
  tmp[4] = "p_IAA-";
  tmp[5] = "D_IAA";
  tmp[6] = "c_PIN";
  tmp[7] = "d_PIN";
  tmp[8] = "endo_PIN";
  tmp[9] = "exo_PIN";
  tmp[10] = "K_hill";
  tmp[11] = "n_hill";
  tmp[12] = "endo_back";
  tmp[13] = "n_hill";
  tmp[14] = "endo_back";
  setParameterId( tmp );
}

void UpInternalGradientModel::
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
  size_t pI = variableIndex(0,1);//pin
  size_t awI = variableIndex(1,0);//auxin (wall)
  size_t pwI = variableIndex(1,1);//pin (membrane/wall)
  size_t mwI = variableIndex(1,2);//(membrane marker)

  assert( aI<cellData[0].size() &&
	  pI<cellData[0].size() &&
	  awI<wallData[0].size() &&
	  pwI<wallData[0].size()  &&
	  mwI<wallData[0].size() );

 for (size_t i=0; i<numCells; ++i) {
	  
    //Production and degradation
    cellDerivs[i][aI] += parameter(0) - parameter(1)*cellData[i][aI];
    cellDerivs[i][pI] += parameter(6)*cellData[i][aI]/(parameter(4)+cellData[i][aI]) - parameter(7)*cellData[i][pI];

    
    //Auxin transport and protein cycling
    size_t numWalls = T.cell(i).numWall();
    for (size_t k=0; k<numWalls; ++k) {
      size_t j = T.cell(i).wall(k)->index();
      if( T.cell(i).wall(k)->cell1()->index() == i && T.cell(i).wall(k)->cell2() != T.background()  && wallData[j][mwI]==1) {
	// cell-cell transport
	size_t iNeighbor = T.cell(i).wall(k)->cell2()->index();
	double fac = parameter(2)*cellData[i][aI] +
	  parameter(3)*cellData[i][aI]*wallData[j][pwI]; //p_2 A_i + p_4 A_i P_ij
	
	cellDerivs[i][aI] -= fac;
	cellDerivs[iNeighbor][aI] += fac;
	
	
	//PIN cycling
	fac = parameter(8)*wallData[j][pwI]
              -	  parameter(9)*cellData[i][pI]*std::pow(cellData[i][aI],parameter(11))/( std::pow(parameter(10),parameter(11)) + std::pow(cellData[i][aI],parameter(11)) );
	wallDerivs[j][pwI] -= fac;
	cellDerivs[i][pI] += fac;


      }
      else if( T.cell(i).wall(k)->cell2()->index() == i && T.cell(i).wall(k)->cell1() != T.background()  && wallData[j][mwI]==1) {
	// cell-cell transport
	size_t iNeighbor = T.cell(i).wall(k)->cell1()->index();
	double fac = parameter(2)*cellData[i][aI] +
	  parameter(3)*cellData[i][aI]*wallData[j][pwI+1];
	
	cellDerivs[i][aI] -= fac;
	cellDerivs[iNeighbor][aI] += fac;


	//PIN cycling

	fac = parameter(8)*wallData[j][pwI+1]
	 -parameter(9)*cellData[i][pI]*std::pow(cellData[i][aI],parameter(11))/( std::pow(parameter(10),parameter(11)) + std::pow(cellData[i][aI],parameter(11)) ) ;

	wallDerivs[j][pwI+1] -= fac;
	cellDerivs[i][pI] += fac;
      }

  else if( T.cell(i).wall(k)->cell2()->index() == i && T.cell(i).wall(k)->cell1() != T.background()  && wallData[j][mwI]==0) {
	// cell-cell transport
	size_t iNeighbor = T.cell(i).wall(k)->cell1()->index();
	double fac = parameter(13)*cellData[i][aI];
	double fac2 = parameter(14)*cellData[i][pI];
	cellDerivs[i][aI] -= fac;
        cellDerivs[i][pI] -= fac2;
	cellDerivs[iNeighbor][aI] += fac;
	cellDerivs[iNeighbor][pI] += fac2;
      }

  else if( T.cell(i).wall(k)->cell1()->index() == i && T.cell(i).wall(k)->cell2() != T.background()  && wallData[j][mwI]==0) {
	// cell-cell transport
	size_t iNeighbor = T.cell(i).wall(k)->cell2()->index();
	double fac = parameter(13)*cellData[i][aI];
	double fac2 = parameter(14)*cellData[i][pI];
	
	cellDerivs[i][aI] -= fac;
        cellDerivs[i][pI] -= fac2;
	cellDerivs[iNeighbor][aI] += fac;
	cellDerivs[iNeighbor][pI] += fac2;

      }
    }
  }
}





DownInternalGradientModel::
DownInternalGradientModel(std::vector<double> &paraValue, 
	      std::vector< std::vector<size_t> > 
	      &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=15 ) {
    std::cerr << "DownInternalGradientModel::"
	      << "DownInternalGradientModel() "
	      << "15 parameters used (see network.h)\n";
    exit(0);
  }
  if( indValue.size() != 2 || indValue[0].size() != 2 || indValue[1].size() != 3 ) {
    std::cerr << "DownInternalGradientModel::"
	      << "DownInternalGradientModel() "
	      << "Two cell variable indices (first row) and three wall variable"
	      << " indices are used (auxin,PIN, membrane marker)." << std::endl;
    exit(0);
  }
  //Set the variable values
  //
  setId("DownInternalGradientModel");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "c_IAA";
  tmp[1] = "d_IAA";
  tmp[2] = "p_IAAH(in)";
  tmp[3] = "p_IAAH(out)";
  tmp[4] = "p_IAA-";
  tmp[5] = "D_IAA";
  tmp[6] = "c_PIN";
  tmp[7] = "d_PIN";
  tmp[8] = "endo_PIN";
  tmp[9] = "exo_PIN";
  tmp[10] = "K_hill";
  tmp[11] = "n_hill";
  tmp[12] = "endo_back";
  tmp[13] = "D_c_auxin";
  tmp[14] = "D_c_pin";
  setParameterId( tmp );
}

void DownInternalGradientModel::
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
  size_t pI = variableIndex(0,1);//pin
  size_t awI = variableIndex(1,0);//auxin (wall)
  size_t pwI = variableIndex(1,1);//pin (membrane/wall)
  size_t mwI = variableIndex(1,2);//(membrane marker)

  assert( aI<cellData[0].size() &&
	  pI<cellData[0].size() &&
	  awI<wallData[0].size() &&
	  pwI<wallData[0].size()  &&
	  mwI<wallData[0].size() );

 for (size_t i=0; i<numCells; ++i) {
	  
    //Production and degradation
    cellDerivs[i][aI] += parameter(0) - parameter(1)*cellData[i][aI];
    cellDerivs[i][pI] += parameter(6)*cellData[i][aI]/(parameter(4)+cellData[i][aI]) - parameter(7)*cellData[i][pI];

    
    //Auxin transport and protein cycling
    size_t numWalls = T.cell(i).numWall();
    for (size_t k=0; k<numWalls; ++k) {
      size_t j = T.cell(i).wall(k)->index();
      if( T.cell(i).wall(k)->cell1()->index() == i && T.cell(i).wall(k)->cell2() != T.background()  && wallData[j][mwI]==1) {
	// cell-cell transport
	size_t iNeighbor = T.cell(i).wall(k)->cell2()->index();
	double fac = parameter(2)*cellData[i][aI] +
	  parameter(3)*cellData[i][aI]*wallData[j][pwI]; //p_2 A_i + p_4 A_i P_ij
	
	cellDerivs[i][aI] -= fac;
	cellDerivs[iNeighbor][aI] += fac;
	
	
	//PIN cycling
	fac = parameter(8)*wallData[j][pwI]
	  -	  parameter(9)*cellData[i][pI]*std::pow(parameter(10),parameter(11))/( std::pow(parameter(10),parameter(11)) + std::pow(cellData[i][aI],parameter(11)) );
	wallDerivs[j][pwI] -= fac;
	cellDerivs[i][pI] += fac;


      }
      else if( T.cell(i).wall(k)->cell2()->index() == i && T.cell(i).wall(k)->cell1() != T.background()  && wallData[j][mwI]==1) {
	// cell-cell transport
	size_t iNeighbor = T.cell(i).wall(k)->cell1()->index();
	double fac = parameter(2)*cellData[i][aI] +
	  parameter(3)*cellData[i][aI]*wallData[j][pwI+1];
	
	cellDerivs[i][aI] -= fac;
	cellDerivs[iNeighbor][aI] += fac;


	//PIN cycling

	fac = parameter(8)*wallData[j][pwI+1]
	  -parameter(9)*cellData[i][pI]*std::pow(parameter(10),parameter(11))/( std::pow(parameter(10),parameter(11)) + std::pow(cellData[i][aI],parameter(11)) ) ;

	wallDerivs[j][pwI+1] -= fac;
	cellDerivs[i][pI] += fac;
      }

  else if( T.cell(i).wall(k)->cell2()->index() == i && T.cell(i).wall(k)->cell1() != T.background()  && wallData[j][mwI]==0) {
	// cell-cell transport
	size_t iNeighbor = T.cell(i).wall(k)->cell1()->index();
	double fac = parameter(13)*cellData[i][aI];
	double fac2 = parameter(14)*cellData[i][pI];
	cellDerivs[i][aI] -= fac;
        cellDerivs[i][pI] -= fac2;
	cellDerivs[iNeighbor][aI] += fac;
	cellDerivs[iNeighbor][pI] += fac2;
      }

  else if( T.cell(i).wall(k)->cell1()->index() == i && T.cell(i).wall(k)->cell2() != T.background()  && wallData[j][mwI]==0) {
	// cell-cell transport
	size_t iNeighbor = T.cell(i).wall(k)->cell2()->index();
	double fac = parameter(13)*cellData[i][aI];
	double fac2 = parameter(14)*cellData[i][pI];
	
	cellDerivs[i][aI] -= fac;
        cellDerivs[i][pI] -= fac2;
	cellDerivs[iNeighbor][aI] += fac;
	cellDerivs[iNeighbor][pI] += fac2;

      }
    }
  }
}


UpExternalGradientModel::
UpExternalGradientModel(std::vector<double> &paraValue,
                          std::vector< std::vector<size_t> >
                          &indValue ) {
    
    //Do some checks on the parameters and variable indeces
    //
    if( paraValue.size()!=15 ) {
        std::cerr << "UpExternalGradientModel::"
        << "UpExternalGradientModel() "
        << "15 parameters used (see network.h)\n";
        exit(0);
    }
    if( indValue.size() != 2 || indValue[0].size() != 2 || indValue[1].size() != 3 ) {
        std::cerr << "UpExternalGradientModel::"
        << "UpExternalGradientModel() "
        << "Two cell variable indices (first row) and three wall variable"
        << " indices are used (auxin,PIN, membrane marker)." << std::endl;
        exit(0);
    }
    //Set the variable values
    //
    setId("UpExternalGradientModel");
    setParameter(paraValue);
    setVariableIndex(indValue);
    
    //Set the parameter identities
    //
    std::vector<std::string> tmp( numParameter() );
    tmp.resize( numParameter() );
    tmp[0] = "c_IAA";
    tmp[1] = "d_IAA";
    tmp[2] = "p_IAAH(in)";
    tmp[3] = "p_IAAH(out)";
    tmp[4] = "p_IAA-";
    tmp[5] = "D_IAA";
    tmp[6] = "c_PIN";
    tmp[7] = "d_PIN";
    tmp[8] = "endo_PIN";
    tmp[9] = "exo_PIN";
    tmp[10] = "K_hill";
    tmp[11] = "n_hill";
    tmp[12] = "endo_back";
    tmp[13] = "D_c_auxin";
    tmp[14] = "D_c_pin";
    setParameterId( tmp );
}

void UpExternalGradientModel::
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
    size_t pI = variableIndex(0,1);//pin
    size_t awI = variableIndex(1,0);//auxin (wall)
    size_t pwI = variableIndex(1,1);//pin (membrane/wall)
    size_t mwI = variableIndex(1,2);//(membrane marker)
    
    assert( aI<cellData[0].size() &&
           pI<cellData[0].size() &&
           awI<wallData[0].size() &&
           pwI<wallData[0].size()  &&
           mwI<wallData[0].size() );
    
    for (size_t i=0; i<numCells; ++i) {
        
        //Production and degradation
        cellDerivs[i][aI] += parameter(0) - parameter(1)*cellData[i][aI];
        cellDerivs[i][pI] += parameter(6)*cellData[i][aI]/(parameter(4)+cellData[i][aI]) - parameter(7)*cellData[i][pI];
        
        
        //Auxin transport and protein cycling
        size_t numWalls = T.cell(i).numWall();
        for (size_t k=0; k<numWalls; ++k) {
            size_t j = T.cell(i).wall(k)->index();
            if( T.cell(i).wall(k)->cell1()->index() == i && T.cell(i).wall(k)->cell2() != T.background()  && wallData[j][mwI]==1) {
                // cell-cell transport
                size_t iNeighbor = T.cell(i).wall(k)->cell2()->index();
                double fac = parameter(2)*cellData[i][aI] +
                parameter(3)*cellData[i][aI]*wallData[j][pwI]; //p_2 A_i + p_4 A_i P_ij
                
                cellDerivs[i][aI] -= fac;
                cellDerivs[iNeighbor][aI] += fac;
                
                
                //PIN cycling
                fac = parameter(8)*wallData[j][pwI]
                -	  parameter(9)*cellData[i][pI]*std::pow(cellData[iNeighbor][aI],parameter(11))/( std::pow(parameter(10),parameter(11)) + std::pow(cellData[iNeighbor][aI],parameter(11)) );
                wallDerivs[j][pwI] -= fac;
                cellDerivs[i][pI] += fac;
                
                
            }
            else if( T.cell(i).wall(k)->cell2()->index() == i && T.cell(i).wall(k)->cell1() != T.background()  && wallData[j][mwI]==1) {
                // cell-cell transport
                size_t iNeighbor = T.cell(i).wall(k)->cell1()->index();
                double fac = parameter(2)*cellData[i][aI] +
                parameter(3)*cellData[i][aI]*wallData[j][pwI+1];
                
                cellDerivs[i][aI] -= fac;
                cellDerivs[iNeighbor][aI] += fac;
                
                
                //PIN cycling
                
                fac = parameter(8)*wallData[j][pwI+1]
                -parameter(9)*cellData[i][pI]*std::pow(cellData[iNeighbor][aI],parameter(11))/( std::pow(parameter(10),parameter(11)) + std::pow(cellData[iNeighbor][aI],parameter(11)) ) ;
                
                wallDerivs[j][pwI+1] -= fac;
                cellDerivs[i][pI] += fac;
            }
            
            else if( T.cell(i).wall(k)->cell2()->index() == i && T.cell(i).wall(k)->cell1() != T.background()  && wallData[j][mwI]==0) {
                // cell-cell transport
                size_t iNeighbor = T.cell(i).wall(k)->cell1()->index();
                double fac = parameter(13)*cellData[i][aI];
                double fac2 = parameter(14)*cellData[i][pI];
                cellDerivs[i][aI] -= fac;
                cellDerivs[i][pI] -= fac2;
                cellDerivs[iNeighbor][aI] += fac;
                cellDerivs[iNeighbor][pI] += fac2;
            }
            
            else if( T.cell(i).wall(k)->cell1()->index() == i && T.cell(i).wall(k)->cell2() != T.background()  && wallData[j][mwI]==0) {
                // cell-cell transport
                size_t iNeighbor = T.cell(i).wall(k)->cell2()->index();
                double fac = parameter(13)*cellData[i][aI];
                double fac2 = parameter(14)*cellData[i][pI];
                
                cellDerivs[i][aI] -= fac;
                cellDerivs[i][pI] -= fac2;
                cellDerivs[iNeighbor][aI] += fac;
                cellDerivs[iNeighbor][pI] += fac2;
                
            }
        }
    }
}




DownInternalGradientModelSingleCell::
DownInternalGradientModelSingleCell(std::vector<double> &paraValue, 
	      std::vector< std::vector<size_t> > 
	      &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=15 ) {
    std::cerr << "DownInternalGradientModelSingleCell::"
	      << "DownInternalGradientModelSingleCell() "
	      << "15 parameters used (see network.h)\n";
    exit(0);
  }
  if( indValue.size() != 2 || indValue[0].size() != 2 || indValue[1].size() != 3 ) {
    std::cerr << "DownInternalGradientModelSingleCell::"
	      << "DownInternalGradientModelSingleCell() "
	      << "Two cell variable indices (first row) and three wall variable"
	      << " indices are used (auxin,PIN, membrane marker)." << std::endl;
    exit(0);
  }
  //Set the variable values
  //
  setId("DownInternalGradientModelSingleCell");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "c_IAA";
  tmp[1] = "d_IAA";
  tmp[2] = "p_IAAH(in)";
  tmp[3] = "p_IAAH(out)";
  tmp[4] = "p_IAA-";
  tmp[5] = "D_IAA";
  tmp[6] = "c_PIN";
  tmp[7] = "d_PIN";
  tmp[8] = "endo_PIN";
  tmp[9] = "exo_PIN";
  tmp[10] = "K_hill";
  tmp[11] = "n_hill";
  tmp[12] = "endo_back";
  tmp[13] = "D_c_auxin";
  tmp[14] = "D_c_pin";
  setParameterId( tmp );
}

void DownInternalGradientModelSingleCell::
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
  size_t pI = variableIndex(0,1);//pin
  size_t awI = variableIndex(1,0);//auxin (wall)
  size_t pwI = variableIndex(1,1);//pin (membrane/wall)
  size_t mwI = variableIndex(1,2);//(membrane marker)

  assert( aI<cellData[0].size() &&
	  pI<cellData[0].size() &&
	  awI<wallData[0].size() &&
	  pwI<wallData[0].size()  &&
	  mwI<wallData[0].size() );

 for (size_t i=0; i<numCells; ++i) {
	  
    //Production and degradation
    cellDerivs[i][aI] += parameter(0) - parameter(1)*cellData[i][aI];
    cellDerivs[i][pI] += parameter(6)*cellData[i][aI]/(parameter(4)+cellData[i][aI]) - parameter(7)*cellData[i][pI];

    
    //Auxin transport and protein cycling
    size_t numWalls = T.cell(i).numWall();
    for (size_t k=0; k<numWalls; ++k) {
      size_t j = T.cell(i).wall(k)->index();
      if( T.cell(i).wall(k)->cell1()->index() == i  && wallData[j][mwI]==1) {
	// cell-cell transport
	size_t iNeighbor = T.cell(i).wall(k)->cell2()->index();
	
	
	//PIN cycling
	double 	fac = parameter(8)*wallData[j][pwI]-parameter(9)*cellData[i][pI]*std::pow(parameter(10),parameter(11))/( std::pow(parameter(10),parameter(11)) + std::pow(cellData[i][aI],parameter(11)) );
	wallDerivs[j][pwI] -= fac;
	cellDerivs[i][pI] += fac;


      }

if( T.cell(i).wall(k)->cell2()->index() == i  && wallData[j][mwI]==1) {
	// cell-cell transport
	size_t iNeighbor = T.cell(i).wall(k)->cell2()->index();
	
double fac = parameter(8)*wallData[j][pwI+1]
	  -parameter(9)*cellData[i][pI]*std::pow(parameter(10),parameter(11))/( std::pow(parameter(10),parameter(11)) + std::pow(cellData[i][aI],parameter(11)) ) ;

	wallDerivs[j][pwI+1] -= fac;
	cellDerivs[i][pI] += fac;


      }
     
  
  else if( T.cell(i).wall(k)->cell1()->index() == i && T.cell(i).wall(k)->cell2() != T.background()  && wallData[j][mwI]==0) {
	// cell-cell transport
	size_t iNeighbor = T.cell(i).wall(k)->cell2()->index();
	double fac = parameter(13)*cellData[i][aI];
	double fac2 = parameter(14)*cellData[i][pI];
	
	cellDerivs[i][aI] -= fac;
        cellDerivs[i][pI] -= fac2;
	cellDerivs[iNeighbor][aI] += fac;
	cellDerivs[iNeighbor][pI] += fac2;

      }
    }
  }
}




AuxinFluxModel::
AuxinFluxModel(std::vector<double> &paraValue, 
	      std::vector< std::vector<size_t> > 
	      &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=11 ) {
    std::cerr << "AuxinFluxModel::"
	      << "AuxinFluxModel() "
	      << "11 parameters used (see network.h)\n";
    exit(0);
  }
  if( indValue.size() != 2 || indValue[0].size() != 2 || indValue[1].size() != 2 ) {
    std::cerr << "AuxinFluxModel::"
	      << "AuxinFluxModel() "
	      << "Two cell variable indices (first row) and two wall variable"
	      << " indices are used (auxin,PIN)." << std::endl;
    exit(0);
  }
  //Set the variable values
  //
  setId("AuxinFluxModel");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "c_IAA";
  tmp[1] = "d_IAA";
  tmp[2] = "D";
  tmp[3] = "T";
  tmp[4] = "K";
  tmp[5] = "n";
  tmp[6] = "c_PIN";
  tmp[7] = "d_PIN";
  tmp[8] = "endo_PIN";
  tmp[9] = "exo_PIN";
  tmp[10] = "endo_back";
  setParameterId( tmp );
}

void AuxinFluxModel::
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
  size_t pI = variableIndex(0,1);//pin
  size_t awI = variableIndex(1,0);//auxin (wall)
  size_t pwI = variableIndex(1,1);//pin (membrane/wall)


  assert( aI<cellData[0].size() &&
	  pI<cellData[0].size() &&
	  awI<wallData[0].size() &&
	  pwI<wallData[0].size() );

 for (size_t i=0; i<numCells; ++i) {
	  
    //Production and degradation
    cellDerivs[i][aI] += parameter(0) - parameter(1)*cellData[i][aI];
    cellDerivs[i][pI] += parameter(6)*cellData[i][aI]/(parameter(4)+cellData[i][aI]) - parameter(7)*cellData[i][pI];

    
    //Auxin transport and protein cycling
    size_t numWalls = T.cell(i).numWall();
    for (size_t k=0; k<numWalls; ++k) {
      size_t j = T.cell(i).wall(k)->index();
      if( T.cell(i).wall(k)->cell1()->index() == i && T.cell(i).wall(k)->cell2() != T.background() ) {
	// cell-cell transport
	size_t iNeighbor = T.cell(i).wall(k)->cell2()->index();
	double fac = parameter(2)*cellData[i][aI] + parameter(3)*cellData[i][aI]*wallData[j][pwI]; //p_2 A_i + p_4 A_i P_ij
	
	cellDerivs[i][aI] -= fac;
	cellDerivs[iNeighbor][aI] += fac;
	if (parameter(2)*cellData[i][aI]-parameter(2)*cellData[iNeighbor][aI]+parameter(3)*cellData[i][aI]*wallData[j][pwI]-parameter(3)*cellData[iNeighbor][aI]*wallData[j][pwI+1]>=0) {
	//PIN cycling
	fac = parameter(8)*wallData[j][pwI] -parameter(10)*cellData[i][pI]
	 - parameter(9)*cellData[i][pI]*std::pow((parameter(2)*cellData[i][aI]-parameter(2)*cellData[iNeighbor][aI])+parameter(3)*cellData[i][aI]*wallData[j][pwI]-parameter(3)*cellData[iNeighbor][aI]*wallData[j][pwI+1],2);}

	else if ((parameter(3)*cellData[i][aI]*wallData[j][pwI]-parameter(3)*cellData[iNeighbor][aI]*wallData[j][pwI+1])<0){
	//PIN cycling
	fac = parameter(8)*wallData[j][pwI] -parameter(10)*cellData[i][pI];}

	wallDerivs[j][pwI] -= fac;
	cellDerivs[i][pI] += fac;


      }
      else if( T.cell(i).wall(k)->cell2()->index() == i && T.cell(i).wall(k)->cell1() != T.background() ) {
	// cell-cell transport
	size_t iNeighbor = T.cell(i).wall(k)->cell1()->index();
	double fac = parameter(2)*cellData[i][aI] +
	  parameter(3)*cellData[i][aI]*wallData[j][pwI+1];
	
	cellDerivs[i][aI] -= fac;
	cellDerivs[iNeighbor][aI] += fac;


	//PIN cycling
	if (parameter(2)*cellData[i][aI]-parameter(2)*cellData[iNeighbor][aI]-parameter(3)*cellData[iNeighbor][aI]*wallData[j][pwI]+parameter(3)*cellData[i][aI]*wallData[j][pwI+1]>=0) {

	  fac = parameter(8)* wallData[j][pwI+1]-parameter(10)*cellData[i][pI]
      -parameter(9)*cellData[i][pI]*std::pow((parameter(2)*cellData[i][aI]-parameter(2)*cellData[iNeighbor][aI])+(-parameter(3)*cellData[iNeighbor][aI]*wallData[j][pwI]+parameter(3)*cellData[i][aI]*wallData[j][pwI+1]),2);}
 else if (-parameter(3)*cellData[iNeighbor][aI]*wallData[j][pwI]+parameter(3)*cellData[i][aI]*wallData[j][pwI+1]<0) {

	fac = parameter(8)*wallData[j][pwI+1] -parameter(10)*cellData[i][pI];}


	wallDerivs[j][pwI+1] -= fac;
	cellDerivs[i][pI] += fac;
      }
    }
  }
}




IntracellularPartitioning::
IntracellularPartitioning(std::vector<double> &paraValue, 
	      std::vector< std::vector<size_t> > 
	      &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=5 ) {
    std::cerr << "IntracellularPartitioning::"
	      << "IntracellularPartitioning() "
	      << "5 parameters used (see network.h)\n";
    exit(0);
  }
  if( indValue.size() != 2 || indValue[0].size() != 2 || indValue[1].size() != 2 ) {
    std::cerr << "IntracellularPartitioning::"
	      << "IntracellularPartitioning() "
	      << "Two cell variable indices (first row) and two wall variable"
	      << " indices are used (auxin,PIN)." << std::endl;
    exit(0);
  }
  //Set the variable values
  //
  setId("IntracellularPartitioning");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "phi";
  tmp[1] = "eta";
  tmp[2] = "mu";
  tmp[3] = "alpha";
  tmp[4] = "D";
  
 
  setParameterId( tmp );
}

void IntracellularPartitioning::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{  
  size_t numCells = T.numCell();
  size_t aI = variableIndex(0,0);
  size_t bI = variableIndex(0,1);
  size_t astarwI = variableIndex(1,0);
  size_t bstarwI = variableIndex(1,1);


  assert( aI<cellData[0].size() &&
	  bI<cellData[0].size() &&
	  astarwI<wallData[0].size() &&
	  bstarwI<wallData[0].size() );

 for (size_t i=0; i<numCells; ++i) {
	  
  
    
    //Auxin transport and protein cycling
    size_t numWalls = T.cell(i).numWall();
    for (size_t k=0; k<numWalls; ++k) {
      size_t j = T.cell(i).wall(k)->index();
      if( T.cell(i).wall(k)->cell1()->index() == i && T.cell(i).wall(k)->cell2() != T.background() ) {


	//PIN cycling
	double fac = (parameter(0)+parameter(1)*wallData[j][astarwI])*cellData[i][aI]-(parameter(2)+parameter(3)*wallData[j][bstarwI])*wallData[j][astarwI];

	double facB = (parameter(0)+parameter(1)*wallData[j][bstarwI])*cellData[i][bI]-(parameter(2)+parameter(3)*wallData[j][astarwI])*wallData[j][bstarwI];



	wallDerivs[j][astarwI] += fac;
	cellDerivs[i][aI] -= fac;
	wallDerivs[j][bstarwI] += facB;
	cellDerivs[i][bI] -= facB;

      }
      else if( T.cell(i).wall(k)->cell2()->index() == i && T.cell(i).wall(k)->cell1() != T.background() ) {

	double fac = (parameter(0)+parameter(1)*wallData[j][astarwI+1])*cellData[i][aI]-(parameter(2)+parameter(3)*wallData[j][bstarwI+1])*wallData[j][astarwI+1];

	double facB = (parameter(0)+parameter(1)*wallData[j][bstarwI+1])*cellData[i][bI]-(parameter(2)+parameter(3)*wallData[j][astarwI+1])*wallData[j][bstarwI+1];



	wallDerivs[j][astarwI+1] += fac;
	cellDerivs[i][aI] -= fac;
	wallDerivs[j][bstarwI+1] += facB;
	cellDerivs[i][bI] -= facB;

      }
    }
  }
}




IntracellularCoupling::
IntracellularCoupling(std::vector<double> &paraValue, 
	      std::vector< std::vector<size_t> > 
	      &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=6 ) {
    std::cerr << "IntracellularCoupling::"
	      << "IntracellularCoupling() "
	      << "6 parameters used (see network.h)\n";
    exit(0);
  }
  if( indValue.size() != 2 || indValue[0].size() != 2 || indValue[1].size() != 2 ) {
    std::cerr << "IntracellularCoupling::"
	      << "IntracellularCoupling() "
	      << "Two cell variable indices (first row) and two wall variable"
	      << " indices are used (auxin,PIN)." << std::endl;
    exit(0);
  }
  //Set the variable values
  //
  setId("IntracellularCoupling");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "phi";
  tmp[1] = "eta";
  tmp[2] = "mu";
  tmp[3] = "alpha";
  tmp[4] = "D";
  tmp[5] = "D";
  
 
  setParameterId( tmp );
}

void IntracellularCoupling::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{  
  size_t numCells = T.numCell();
  size_t aI = variableIndex(0,0);
  size_t bI = variableIndex(0,1);
  size_t astarwI = variableIndex(1,0);
  size_t bstarwI = variableIndex(1,1);


  assert( aI<cellData[0].size() &&
	  bI<cellData[0].size() &&
	  astarwI<wallData[0].size() &&
	  bstarwI<wallData[0].size() );

 for (size_t i=0; i<numCells; ++i) {
	  
  
    
    //Auxin transport and protein cycling
    size_t numWalls = T.cell(i).numWall();
    for (size_t k=0; k<numWalls; ++k) {
      size_t j = T.cell(i).wall(k)->index();
      if( T.cell(i).wall(k)->cell1()->index() == i && T.cell(i).wall(k)->cell2() != T.background() ) {


	//PIN cycling
	double fac = (parameter(0)+parameter(1)*wallData[j][astarwI])*cellData[i][aI]-(parameter(2)+parameter(3)*wallData[j][bstarwI]+parameter(5)*wallData[j][bstarwI]*wallData[j][astarwI+1])*wallData[j][astarwI];

      double facB = (parameter(0)+parameter(1)*wallData[j][bstarwI])*cellData[i][bI]-(parameter(2)+parameter(3)*wallData[j][astarwI])*wallData[j][bstarwI];



	wallDerivs[j][astarwI] += fac;
	cellDerivs[i][aI] -= fac;
	wallDerivs[j][bstarwI] += facB;
	cellDerivs[i][bI] -= facB;

      }
      else if( T.cell(i).wall(k)->cell2()->index() == i && T.cell(i).wall(k)->cell1() != T.background() ) {

	double fac = (parameter(0)+parameter(1)*wallData[j][astarwI+1])*cellData[i][aI]-(parameter(2)+parameter(3)*wallData[j][bstarwI+1] +parameter(5)*wallData[j][bstarwI+1]*wallData[j][astarwI])*wallData[j][astarwI+1];

	double facB = (parameter(0)+parameter(1)*wallData[j][bstarwI+1])*cellData[i][bI]-(parameter(2)+parameter(3)*wallData[j][astarwI+1])*wallData[j][bstarwI+1];



	wallDerivs[j][astarwI+1] += fac;
	cellDerivs[i][aI] -= fac;
	wallDerivs[j][bstarwI+1] += facB;
	cellDerivs[i][bI] -= facB;

      }
    }
  }
}





IntracellularIndirectCoupling::
IntracellularIndirectCoupling(std::vector<double> &paraValue, 
	      std::vector< std::vector<size_t> > 
	      &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=12 ) {
    std::cerr << "IntracellularIndirectCoupling::"
	      << "IntracellularIndirectCoupling() "
	      << "11 parameters used (see network.h)\n";
    exit(0);
  }
  if( indValue.size() != 2 || indValue[0].size() != 3 || indValue[1].size() != 3 ) {
    std::cerr << "IntracellularIndirectCoupling::"
	      << "IntracellularIndirectCoupling() "
	      << "Two cell variable indices (first row) and two wall variable"
	      << " indices are used (auxin,PIN)." << std::endl;
    exit(0);
  }
  //Set the variable values
  //
  setId("IntracellularIndirectCoupling");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "phi";
  tmp[1] = "eta";
  tmp[2] = "mu";
  tmp[3] = "alpha";
  tmp[4] = "D";
  tmp[5] = "D";
  tmp[6] = "alpha";
  tmp[7] = "D";
  tmp[8] = "D";
  tmp[9] = "D";
  tmp[10] = "D";
  tmp[11] = "D";
 
  setParameterId( tmp );
}

void IntracellularIndirectCoupling::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{  
  size_t numCells = T.numCell();
  size_t aI = variableIndex(0,0);
  size_t bI = variableIndex(0,1);
  size_t mI = variableIndex(0,2);
  size_t astarwI = variableIndex(1,0);
  size_t bstarwI = variableIndex(1,1);
  size_t mwI = variableIndex(1,2);


  assert( aI<cellData[0].size() &&
	  bI<cellData[0].size() &&
	  mI<cellData[0].size() &&
	  astarwI<wallData[0].size() &&
	  bstarwI<wallData[0].size()  &&
          mwI<wallData[0].size()) ;


 for (size_t i=0; i<numCells; ++i) {
	  
  
 //Production and degradation
    cellDerivs[i][mI] += parameter(9) - parameter(10)*cellData[i][mI];
 
    
    //Auxin transport and protein cycling
    size_t numWalls = T.cell(i).numWall();
    for (size_t k=0; k<numWalls; ++k) {
      size_t j = T.cell(i).wall(k)->index();
      if( T.cell(i).wall(k)->cell1()->index() == i && T.cell(i).wall(k)->cell2() != T.background() ) {
	double fac = parameter(11)*wallData[j][mwI];
	wallDerivs[j][mwI] -= fac;
	wallDerivs[j][mwI+1] += fac;
	//PIN cycling
	 fac = (parameter(0)+parameter(1)*wallData[j][astarwI])*cellData[i][aI]-(parameter(2)+parameter(3)*wallData[j][bstarwI]+parameter(5)*wallData[j][mwI])*wallData[j][astarwI];

      double facB = (parameter(0)+parameter(1)*wallData[j][bstarwI])*cellData[i][bI]-(parameter(2)+parameter(3)*wallData[j][astarwI])*wallData[j][bstarwI];
      double facC = -parameter(6)*wallData[j][mwI]+parameter(7)*cellData[i][mI]+parameter(8)*wallData[j][bstarwI]*cellData[i][mI];


	wallDerivs[j][astarwI] += fac;
	cellDerivs[i][aI] -= fac;
	wallDerivs[j][bstarwI] += facB;
	cellDerivs[i][bI] -= facB;
	wallDerivs[j][mwI] += facC;
	cellDerivs[i][mI] -= facC;

      }
      else if( T.cell(i).wall(k)->cell2()->index() == i && T.cell(i).wall(k)->cell1() != T.background() ) {
	// wall-wall diffusion
	wallDerivs[j][mwI+1] -= parameter(11)*wallData[j][mwI+1];
	wallDerivs[j][mwI] += parameter(11)*wallData[j][mwI+1];
	double fac = (parameter(0)+parameter(1)*wallData[j][astarwI+1])*cellData[i][aI]-(parameter(2)+parameter(3)*wallData[j][bstarwI+1] +parameter(5)*wallData[j][mwI+1])*wallData[j][astarwI+1];

	double facB = (parameter(0)+parameter(1)*wallData[j][bstarwI+1])*cellData[i][bI]-(parameter(2)+parameter(3)*wallData[j][astarwI+1])*wallData[j][bstarwI+1];

    double facC = -parameter(6)*wallData[j][mwI+1]+parameter(7)*cellData[i][mI]+parameter(8)*wallData[j][bstarwI+1]*cellData[i][mI];

	wallDerivs[j][astarwI+1] += fac;
	cellDerivs[i][aI] -= fac;
	wallDerivs[j][bstarwI+1] += facB;
	cellDerivs[i][bI] -= facB;
	wallDerivs[j][mwI+1] += facC;
	cellDerivs[i][mI] -= facC;



      }
    }
  }
}























DownInternalGradientModelGeometric::
DownInternalGradientModelGeometric(std::vector<double> &paraValue, 
	      std::vector< std::vector<size_t> > 
	      &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=15 ) {
    std::cerr << "DownInternalGradientModelGeometric::"
	      << "DownInternalGradientModelGeometric() "
	      << "15 parameters used (see network.h)\n";
    exit(0);
  }
  if( indValue.size() != 2 || indValue[0].size() != 2 || indValue[1].size() != 3 ) {
    std::cerr << "DownInternalGradientModelGeometric::"
	      << "DownInternalGradientModelGeometric() "
	      << "Two cell variable indices (first row) and three wall variable"
	      << " indices are used (auxin,PIN, membrane marker)." << std::endl;
    exit(0);
  }
  //Set the variable values
  //
  setId("DownInternalGradientModelGeometric");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "c_IAA";
  tmp[1] = "d_IAA";
  tmp[2] = "p_IAAH(in)";
  tmp[3] = "p_IAAH(out)";
  tmp[4] = "p_IAA-";
  tmp[5] = "D_IAA";
  tmp[6] = "c_PIN";
  tmp[7] = "d_PIN";
  tmp[8] = "endo_PIN";
  tmp[9] = "exo_PIN";
  tmp[10] = "K_hill";
  tmp[11] = "n_hill";
  tmp[12] = "endo_back";
  tmp[13] = "D_c_auxin";
  tmp[14] = "D_c_pin";
  setParameterId( tmp );
}

void DownInternalGradientModelGeometric::
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
  size_t pI = variableIndex(0,1);//pin
  size_t awI = variableIndex(1,0);//auxin (wall)
  size_t pwI = variableIndex(1,1);//pin (membrane/wall)
  size_t mwI = variableIndex(1,2);//(membrane marker)

  assert( aI<cellData[0].size() &&
	  pI<cellData[0].size() &&
	  awI<wallData[0].size() &&
	  pwI<wallData[0].size()  &&
	  mwI<wallData[0].size() );

 
 

 for (size_t i=0; i<numCells; ++i) {
	  
    //Production and degradation
    cellDerivs[i][aI] += parameter(0) - parameter(1)*cellData[i][aI];
    cellDerivs[i][pI] += parameter(6)*cellData[i][aI]/(parameter(4)+cellData[i][aI]) - parameter(7)*cellData[i][pI];

    
    //Auxin transport and protein cycling
    size_t numWalls = T.cell(i).numWall();
    double cellVolume = T.cell(i).calculateVolume(vertexData);	
    for (size_t k=0; k<numWalls; ++k) {
      size_t j = T.cell(i).wall(k)->index();
      double lengthWall = T.cell(i).wall(k)->length();
      double g_ij = lengthWall/cellVolume;



      // Checks if cell i is first or second neighbor to wall (wall variables are stored as pairs)
      if( T.cell(i).wall(k)->cell1()->index() == i && T.cell(i).wall(k)->cell2() != T.background()  && wallData[j][mwI]==1) {

	
	// cell-cell transport
	size_t iNeighbor = T.cell(i).wall(k)->cell2()->index();
  	double g_ji = lengthWall/T.cell(iNeighbor).calculateVolume(vertexData);
	double fac = parameter(2)*cellData[i][aI] +parameter(3)*cellData[i][aI]*wallData[j][pwI]; //p_2 A_i + p_4 A_i P_ij
	
	cellDerivs[i][aI] -=  g_ij*fac;
	cellDerivs[iNeighbor][aI] += g_ji*fac;
	
	
	//PIN cycling
	fac = parameter(8)*wallData[j][pwI]
	  -	  parameter(9)*cellData[i][pI]*std::pow(parameter(10),parameter(11))/( std::pow(parameter(10),parameter(11)) + std::pow(cellData[i][aI],parameter(11)) );
	wallDerivs[j][pwI] -= fac;
	cellDerivs[i][pI] += g_ij*fac;


      }
      else if( T.cell(i).wall(k)->cell2()->index() == i && T.cell(i).wall(k)->cell1() != T.background()  && wallData[j][mwI]==1) {
	// cell-cell transport
	size_t iNeighbor = T.cell(i).wall(k)->cell1()->index();
        double g_ji = lengthWall/T.cell(iNeighbor).calculateVolume(vertexData);
	double fac = parameter(2)*cellData[i][aI] +
	  parameter(3)*cellData[i][aI]*wallData[j][pwI+1];
	
	cellDerivs[i][aI] -=  g_ij*fac;
	cellDerivs[iNeighbor][aI] +=  g_ji*fac;


	//PIN cycling

	fac = parameter(8)*wallData[j][pwI+1]
	  -parameter(9)*cellData[i][pI]*std::pow(parameter(10),parameter(11))/( std::pow(parameter(10),parameter(11)) + std::pow(cellData[i][aI],parameter(11)) ) ;

	wallDerivs[j][pwI+1] -= fac;
	cellDerivs[i][pI] += g_ij*fac;
      }

  else if( T.cell(i).wall(k)->cell2()->index() == i && T.cell(i).wall(k)->cell1() != T.background()  && wallData[j][mwI]==0) {
	// cell-cell transport
	size_t iNeighbor = T.cell(i).wall(k)->cell1()->index();
        double g_ji = lengthWall/T.cell(iNeighbor).calculateVolume(vertexData);

	double fac = parameter(13)*cellData[i][aI];
	double fac2 = parameter(14)*cellData[i][pI];
	cellDerivs[i][aI] -= g_ij* fac;
        cellDerivs[i][pI] -= g_ij* fac2;
	cellDerivs[iNeighbor][aI] += g_ji*fac;
	cellDerivs[iNeighbor][pI] += g_ji*fac2;
      }

  else if( T.cell(i).wall(k)->cell1()->index() == i && T.cell(i).wall(k)->cell2() != T.background()  && wallData[j][mwI]==0) {
	// cell-cell transport
	size_t iNeighbor = T.cell(i).wall(k)->cell2()->index();
        double  g_ji = lengthWall/T.cell(iNeighbor).calculateVolume(vertexData);
	double fac = parameter(13)*cellData[i][aI];
	double fac2 = parameter(14)*cellData[i][pI];
	
	cellDerivs[i][aI] -=g_ij* fac;
        cellDerivs[i][pI] -=g_ij* fac2;
	cellDerivs[iNeighbor][aI] +=g_ji* fac;
	cellDerivs[iNeighbor][pI] +=g_ji* fac2;

      }
    }
  }
}







