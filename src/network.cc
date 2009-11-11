///
/// Filename     : network.cc
/// Description  : Classes describing complete models
/// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
/// Created      : November 2006
/// Revision     : $Id:$
///
#include "network.h"
#include "baseReaction.h"

AuxinModelSimple1::
AuxinModelSimple1(std::vector<double> &paraValue, 
		  std::vector< std::vector<size_t> > 
		  &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //////////////////////////////////////////////////////////////////////
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
    exit(0);
  }
  if( indValue.size() != 1 || indValue[0].size() != 4 ) {
    std::cerr << "AuxinModelSimple1::"
	      << "AuxinModelSimple1() "
	      << "Four variable indices are used (auxin,pin,X,M).\n";
    exit(0);
  }
  //Set the variable values
  //////////////////////////////////////////////////////////////////////
  setId("AuxinModelSimple1");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //////////////////////////////////////////////////////////////////////
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
       std::vector< std::vector<double> > &cellData,
       std::vector< std::vector<double> > &wallData,
       std::vector< std::vector<double> > &vertexData,
       std::vector< std::vector<double> > &cellDerivs,
       std::vector< std::vector<double> > &wallDerivs,
       std::vector< std::vector<double> > &vertexDerivs ) 
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

AuxinModelStress::
AuxinModelStress(std::vector<double> &paraValue, 
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
  if( indValue.size() != 2 || indValue[0].size() != 2 || indValue[1].size() != 3 ) {
    std::cerr << "AuxinModelSimpleStress::"
							<< "AuxinModelSimpleStress() "
							<< "Two cell variable indices are used in first level (auxin,pin),\n"
							<< "and three wall indices are used in second level (F,k_1,k_2),\n";
    exit(0);
  }
  //Set the variable values
  //////////////////////////////////////////////////////////////////////
  setId("AuxinModelStress");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //////////////////////////////////////////////////////////////////////
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
       std::vector< std::vector<double> > &cellData,
       std::vector< std::vector<double> > &wallData,
       std::vector< std::vector<double> > &vertexData,
       std::vector< std::vector<double> > &cellDerivs,
       std::vector< std::vector<double> > &wallDerivs,
       std::vector< std::vector<double> > &vertexDerivs ) 
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
  if( indValue.size() != 2 || indValue[0].size() != 2 || indValue[1].size() != 3 ) {
    std::cerr << "AuxinModelSimpleStress::"
							<< "AuxinModelSimpleStress() "
							<< "Two cell variable indices are used in first level (auxin,pin),\n"
							<< "and three wall indices are used in second level (F,k_1,k_2),\n";
    exit(0);
  }
  //Set the variable values
  //////////////////////////////////////////////////////////////////////
  setId("AuxinModelSimpleStress");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //////////////////////////////////////////////////////////////////////
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
       std::vector< std::vector<double> > &cellData,
       std::vector< std::vector<double> > &wallData,
       std::vector< std::vector<double> > &vertexData,
       std::vector< std::vector<double> > &cellDerivs,
       std::vector< std::vector<double> > &wallDerivs,
       std::vector< std::vector<double> > &vertexDerivs ) 
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
  //////////////////////////////////////////////////////////////////////
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
  //////////////////////////////////////////////////////////////////////
  setId("AuxinModelSimple1Wall");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //////////////////////////////////////////////////////////////////////
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
       std::vector< std::vector<double> > &cellData,
       std::vector< std::vector<double> > &wallData,
       std::vector< std::vector<double> > &vertexData,
       std::vector< std::vector<double> > &cellDerivs,
       std::vector< std::vector<double> > &wallDerivs,
       std::vector< std::vector<double> > &vertexDerivs ) 
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

//!Constructor
AuxinModelSimple2::
AuxinModelSimple2(std::vector<double> &paraValue, 
									std::vector< std::vector<size_t> > 
									&indValue ) 
{ 
  //Do some checks on the parameters and variable indeces
  //////////////////////////////////////////////////////////////////////
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
  //////////////////////////////////////////////////////////////////////
  setId("AuxinModelSimple2");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //////////////////////////////////////////////////////////////////////
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

//! Derivative contribution for the growth
/*! Deriving the time derivative contribution for the growth for all
  walls in the tissue.
*/
void AuxinModelSimple2::
derivs(Tissue &T,
       std::vector< std::vector<double> > &cellData,
       std::vector< std::vector<double> > &wallData,
       std::vector< std::vector<double> > &vertexData,
       std::vector< std::vector<double> > &cellDerivs,
       std::vector< std::vector<double> > &wallDerivs,
       std::vector< std::vector<double> > &vertexDerivs ) 
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
		//////////////////////////////////////////////////////////////////////
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

//!Constructor
AuxinModelSimple3::
AuxinModelSimple3(std::vector<double> &paraValue, 
									std::vector< std::vector<size_t> > 
									&indValue ) 
{ 
  //Do some checks on the parameters and variable indeces
  //////////////////////////////////////////////////////////////////////
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
  //////////////////////////////////////////////////////////////////////
  setId("AuxinModelSimple3");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //////////////////////////////////////////////////////////////////////
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

//! Derivative contribution for the growth
/*! Deriving the time derivative contribution for the growth for all
  walls in the tissue.
*/
void AuxinModelSimple3::
derivs(Tissue &T,
       std::vector< std::vector<double> > &cellData,
       std::vector< std::vector<double> > &wallData,
       std::vector< std::vector<double> > &vertexData,
       std::vector< std::vector<double> > &cellDerivs,
       std::vector< std::vector<double> > &wallDerivs,
       std::vector< std::vector<double> > &vertexDerivs ) 
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
		//////////////////////////////////////////////////////////////////////
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

//!Constructor
AuxinModel4::
AuxinModel4(std::vector<double> &paraValue, 
	    std::vector< std::vector<size_t> > 
	    &indValue ) 
{ 
  //Do some checks on the parameters and variable indeces
  //////////////////////////////////////////////////////////////////////
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
  //////////////////////////////////////////////////////////////////////
  setId("AuxinModel4");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //////////////////////////////////////////////////////////////////////
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

//! Derivative contribution for the growth
/*! Deriving the time derivative contribution for the growth for all
  walls in the tissue.
*/
void AuxinModel4::
derivs(Tissue &T,
       std::vector< std::vector<double> > &cellData,
       std::vector< std::vector<double> > &wallData,
       std::vector< std::vector<double> > &vertexData,
       std::vector< std::vector<double> > &cellDerivs,
       std::vector< std::vector<double> > &wallDerivs,
       std::vector< std::vector<double> > &vertexDerivs ) 
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
    //////////////////////////////////////////////////////////////////////
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
  //////////////////////////////////////////////////////////////////////
  if( paraValue.size()!=16 ) {
    std::cerr << "AuxinModel5::"
							<< "AuxinModel5() "
							<< "Ten parameters used (see network.h)\n";
    exit(0);
  }
  if( indValue.size() != 1 || indValue[0].size() != 4 ) {
    std::cerr << "AuxinModel5::"
							<< "AuxinModel5() "
							<< "Four variable indices are used (auxin,pin,X,M).\n";
    exit(0);
  }
  //Set the variable values
  //////////////////////////////////////////////////////////////////////
  setId("AuxinModel5");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //////////////////////////////////////////////////////////////////////
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
       std::vector< std::vector<double> > &cellData,
       std::vector< std::vector<double> > &wallData,
       std::vector< std::vector<double> > &vertexData,
       std::vector< std::vector<double> > &cellDerivs,
       std::vector< std::vector<double> > &wallDerivs,
       std::vector< std::vector<double> > &vertexDerivs ) 
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
		
    cellDerivs[i][xI] += parameter(8)*cellData[i][aI]*cellData[i][aI]/(2.0+cellData[i][aI]*cellData[i][aI]) 
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
  //////////////////////////////////////////////////////////////////////
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
  //////////////////////////////////////////////////////////////////////
  setId("AuxinModel6");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //////////////////////////////////////////////////////////////////////
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
       std::vector< std::vector<double> > &cellData,
       std::vector< std::vector<double> > &wallData,
       std::vector< std::vector<double> > &vertexData,
       std::vector< std::vector<double> > &cellDerivs,
       std::vector< std::vector<double> > &wallDerivs,
       std::vector< std::vector<double> > &vertexDerivs ) 
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
  //////////////////////////////////////////////////////////////////////
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
  //////////////////////////////////////////////////////////////////////
  setId("AuxinModel7");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //////////////////////////////////////////////////////////////////////
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
       std::vector< std::vector<double> > &cellData,
       std::vector< std::vector<double> > &wallData,
       std::vector< std::vector<double> > &vertexData,
       std::vector< std::vector<double> > &cellDerivs,
       std::vector< std::vector<double> > &wallDerivs,
       std::vector< std::vector<double> > &vertexDerivs ) 
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
       std::vector< std::vector<double> > &cellData,
       std::vector< std::vector<double> > &wallData,
       std::vector< std::vector<double> > &vertexData,
       std::vector< std::vector<double> > &cellDerivs,
       std::vector< std::vector<double> > &wallDerivs,
       std::vector< std::vector<double> > &vertexDerivs ) 
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

