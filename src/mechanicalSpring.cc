//
// Filename     : mechanicalSpring.cc
// Description  : Classes describing updates due to mechanical spring interactions
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : September 2007
// Revision     : $Id:$
//
#include <utility>
#include <vector>
#include "baseReaction.h"
#include "mechanicalSpring.h"
#include "tissue.h"

VertexFromWallSpring::
VertexFromWallSpring(std::vector<double> &paraValue, 
		     std::vector< std::vector<size_t> > 
		     &indValue ) 
{  
  // Do some checks on the parameters and variable indeces
  if( paraValue.size()!=2 && paraValue.size()!=3 ) {
    std::cerr << "VertexFromWallSpring::"
	      << "VertexFromWallSpring() "
	      << "Uses two parameters K_force frac_adhesion.\n"
	      << "or \n"
	      << "Uses three parameters K_force1,frac_adhesion, K_force2.\n";
    exit(0);
  }
  
  if( (indValue.size() !=1  &&  indValue.size() !=2  && 
       (indValue.size()!=3 ||  indValue[2].size() != 1 ||( indValue[1].size() != 0  &&  indValue[1].size() != 1 )) 
       ) || indValue[0].size() != 1)  {
    std::cerr << "VertexFromWallSpring::"
     	      << "VertexFromWallSpring() "
     	      << "Wall length index given in first level. "
     	      << "If two levels given, force save index given in second level. "
	      << "If three levels given optionalforce save index given in second level "
	      << "and wall_type_index given in third level. \n";
    exit(0);
  }
  
  // Set the variable values
  setId("VertexFromWallSpring");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  // Set the parameter identities
  std::vector<std::string> tmp( numParameter() );
  if(numParameter()==2 ){
    tmp[0] = "K_force";
    tmp[1] = "frac_adh";
  }
  if(numParameter()==3 ){
    tmp[0] = "K_force1";
    tmp[1] = "frac_adh";
    tmp[2] = "K_force2";
  }

  setParameterId( tmp );
}

void VertexFromWallSpring::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) {
  
  //Do the update for each wall
  size_t numWalls = T.numWall();
  size_t wallLengthIndex = variableIndex(0,0);
  
  // internal wall indices aorta templates
  // std::cerr<<".............begin................"<<std::endl;
  // for( size_t i=0 ; i<numWalls ; ++i ) {
  //   size_t v1 = T.wall(i).vertex1()->index();
  //   size_t v2 = T.wall(i).vertex2()->index();
  //   size_t dimension = vertexData[v1].size();
  //   assert( vertexData[v2].size()==dimension );
  //   double distance1=0.0;
  //   double distance2=0.0;
  //   for( size_t d=0 ; d<dimension ; d++ ){
  //     distance1 += vertexData[v1][d]*vertexData[v1][d];
  //     distance2 += vertexData[v2][d]*vertexData[v2][d];
  //   }
  //   if (distance1<.98 || distance2<.98) std::cerr<< i <<std::endl;
  // }

  // std::cerr<<".............end ................."<<std::endl;


  for( size_t i=0 ; i<numWalls ; ++i ) {
    size_t v1 = T.wall(i).vertex1()->index();
    size_t v2 = T.wall(i).vertex2()->index();
    size_t dimension = vertexData[v1].size();
    assert( vertexData[v2].size()==dimension );
    //Calculate shared factors
    double distance=0.0;
    for( size_t d=0 ; d<dimension ; d++ )
      distance += (vertexData[v1][d]-vertexData[v2][d])*
	(vertexData[v1][d]-vertexData[v2][d]);
    distance = std::sqrt(distance);
    double wallLength=wallData[i][wallLengthIndex];
    double coeff = parameter(0)*((1.0/wallLength)-(1.0/distance));
    
    // Use different spring elasticity if wall type is provided in wall vector
    if(numParameter()==3 && wallData[i][variableIndex(2,0)] ==1 ){
      coeff = parameter(2)*((1.0/wallLength)-(1.0/distance));
    }
      
    if( distance <= 0.0 && wallLength <=0.0 ) {
      //std::cerr << i << " - " << wallLength << " " << distance << std::endl;
      coeff = 0.0;
    }
    if( distance>wallLength )
      coeff *=parameter(1);
    
    //Save force in wall variable if appropriate
    if( numVariableIndexLevel()==2 ||(numParameter()==3 && numVariableIndex(1)==1) ) 
        wallData[i][variableIndex(1,0)] = coeff*distance;
    
    //Update both vertices for each dimension
    for(size_t d=0 ; d<dimension ; d++ ) {
      double div = (vertexData[v1][d]-vertexData[v2][d])*coeff;
      vertexDerivs[v1][d] -= div;
      vertexDerivs[v2][d] += div;
    }
  }
}

namespace CenterTriangulation {
  EdgeSpring::
  EdgeSpring(std::vector<double> &paraValue, 
	     std::vector< std::vector<size_t> > 
	     &indValue ) 
  {  
    // Do some checks on the parameters and variable indeces
    if( paraValue.size()!=2 ) {
      std::cerr << "CenterTriangulation::EdgeSpring"
		<< "EdgeSpring() "
		<< "Uses two parameters K_force frac_adhesion.\n";
      exit(EXIT_FAILURE);
    }
    if( indValue.size() != 1 || indValue[0].size() != 1 ) { 
      std::cerr << "CenterTriangulation::EdgeSpring"
		<< "EdgeSpring() "
		<< "Start of additional Cell variable indices (center(x,y,z) "
		<< "L_1,...,L_n, n=num vertex) is given in first level." 
		<< std::endl;
      exit(EXIT_FAILURE);
    }
    
    // Set the variable values
    setId("CenterTriangulation::EdgeSpring");
    setParameter(paraValue);  
    setVariableIndex(indValue);
    
    // Set the parameter identities
    std::vector<std::string> tmp( numParameter() );
    tmp[0] = "K_force";
    tmp[1] = "frac_adh";
    setParameterId( tmp );
  }
  
  void EdgeSpring::
  derivs(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs ) {
    
    //Do the update for each internal edge for eache cell
    size_t numCells = T.numCell();
    size_t posIndex = variableIndex(0,0);
    size_t dimension = vertexData[0].size();
    assert( 3==dimension );//assuming 3D
    size_t lengthIndex = posIndex+dimension;
    
    for (size_t i=0; i<numCells; ++i) {
      for (size_t k=0; k<T.cell(i).numVertex(); ++k) {
	size_t v = T.cell(i).vertex(k)->index();
	
	//Calculate shared factors
	double distance=0.0;
	for( size_t d=0 ; d<dimension ; d++ ) {
	  distance += (vertexData[v][d]-cellData[i][posIndex+d])*
	    (vertexData[v][d]-cellData[i][posIndex+d]);
	}
	distance = std::sqrt(distance);
	double edgeLength=cellData[i][lengthIndex+k];
	double coeff = parameter(0)*((1.0/edgeLength)-(1.0/distance));
	if( distance <= 0.0 && edgeLength <=0.0 ) {
	  //std::cerr << i << " " << k << " - " << edgeLength << " " 
	  //<< distance << std::endl;
	  coeff = 0.0;
	}
	if( distance>edgeLength )
	  coeff *=parameter(1);
	
	// Save force in wall variable if appropriate
	//if( numVariableIndexLevel()>1 )
	//wallData[i][variableIndex(1,0)] = coeff*distance;
	
	//Update both vertices for each dimension
	for(size_t d=0 ; d<dimension ; d++ ) {
	  double div = (vertexData[v][d]-cellData[i][posIndex+d])*coeff;
	  vertexDerivs[v][d] -= div;
	  cellDerivs[i][posIndex+d] += div;
	}
      }
    }
  }
}

VertexFromDoubleWallSpring::
VertexFromDoubleWallSpring(std::vector<double> &paraValue, 
			   std::vector< std::vector<size_t> > 
			   &indValue ) 
{  
  // Do some checks on the parameters and variable indeces
  if( paraValue.size()!=2 ) {
    std::cerr << "VertexFromDoubleWallSpring::"
	      << "VertexFromDoubleWallSpring() "
	      << "Uses two parameters K_force frac_adhesion.\n";
    exit(0);
  }
  if( indValue.size() < 2 || indValue.size() > 3 
      || indValue[0].size() != 1
      || indValue[1].size() != 2
      || (indValue.size()==3 && indValue[2].size() != 1) ) {
    std::cerr << "VertexFromDoubleWallSpring::"
	      << "VertexFromDoubleWallSpring() "
	      << "Wall length index given in first level,"
	      << " the two wall k variable indices in second,"
	      << " and optionally wall variable save index in third.\n";
    exit(0);
  }
  
  // Set the variable values
  setId("VertexFromDoubleWallSpring");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  // Set the parameter identities
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "K_force";
  tmp[1] = "frac_adh";
  setParameterId( tmp );
}

void VertexFromDoubleWallSpring::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{  
  //Do the update for each wall
  size_t numWalls = T.numWall();
  size_t wallLengthIndex = variableIndex(0,0);
  size_t wall1kIndex = variableIndex(1,0);
  size_t wall2kIndex = variableIndex(1,1);
  for( size_t i=0 ; i<numWalls ; ++i ) {
    size_t v1 = T.wall(i).vertex1()->index();
    size_t v2 = T.wall(i).vertex2()->index();
    size_t dimension = vertexData[v1].size();
    assert( vertexData[v2].size()==dimension );

    //Calculate shared factors
    double distance=0.0;
    for( size_t d=0 ; d<dimension ; d++ )
      distance += (vertexData[v1][d]-vertexData[v2][d])*
	(vertexData[v1][d]-vertexData[v2][d]);
    distance = std::sqrt(distance);
    double wallLength=wallData[i][wallLengthIndex];
    double coeff = (wallData[i][wall1kIndex]+wallData[i][wall2kIndex])*parameter(0)*((1.0/wallLength)-(1.0/distance));
    if( distance <= 0.0 || wallLength <=0.0 ) {
      std::cerr << i << " - " << wallLength << " " << distance << std::endl;
      coeff = 0.0;
    }
    if( distance>wallLength )
      coeff *=parameter(1);
    
    //Save force in wall variable if appropriate
    if( numVariableIndexLevel()>2 )
      wallData[i][variableIndex(2,0)] = coeff*distance;
    
    //Update both vertices for each dimension
    for(size_t d=0 ; d<dimension ; d++ ) {
      double div = (vertexData[v1][d]-vertexData[v2][d])*coeff;
      vertexDerivs[v1][d] -= div;
      vertexDerivs[v2][d] += div;
    }
  }
}

VertexFromWallSpringSpatial::
VertexFromWallSpringSpatial(std::vector<double> &paraValue, 
														std::vector< std::vector<size_t> > 
														&indValue ) 
{  
  // Do some checks on the parameters and variable indeces
  if( paraValue.size()!=5 ) {
    std::cerr << "VertexFromWallSpringSpatial::"
							<< "VertexFromWallSpringSpatial() "
							<< "Uses five parameters k_min k_max K_spatial n_spatial frac_adhesion.\n";
    exit(0);
  }
  if( indValue.size() < 1 || indValue.size() > 2 
			|| indValue[0].size() != 2 
			|| (indValue.size()==2 && indValue[1].size() != 1) ) {
    std::cerr << "VertexFromWallSpringSpatial::"
							<< "VertexFromWallSpringSpatial() "
							<< "Wall length index and spatial coordinate given in first level,"
							<< " and optionally wall variable save index in second.\n";
    exit(0);
  }
	
  // Set the variable values
  setId("VertexFromWallSpringSpatial");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  // Set the parameter identities
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "k_min";
  tmp[1] = "k_max";
  tmp[2] = "K_spatial";
  tmp[3] = "n_spatial";
  tmp[4] = "frac_adh";
  setParameterId( tmp );
	Kpow_=std::pow(paraValue[2],paraValue[3]);
}

void VertexFromWallSpringSpatial::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) {
  
  // Do the update for each wall
  size_t numWalls = T.numWall();
  size_t wallLengthIndex = variableIndex(0,0);
	size_t dimension = vertexData[0].size();
  
	// Initiate positional factor
	size_t sI = variableIndex(0,1);
	size_t numVertices = T.numVertex();
	assert (sI<vertexData[0].size());
	double max = vertexData[0][sI];
	size_t maxI = 0;
	for (size_t i=1; i<numVertices; ++i)
		if (vertexData[i][sI]>max) {
			max=vertexData[i][sI];
			maxI = i;
		}
	std::vector<double> maxPos(dimension);
	for (size_t d=0; d<dimension; ++d)
		maxPos[d] = vertexData[maxI][d];
	
	// Calculate update from each wall
  for( size_t i=0 ; i<numWalls ; ++i ) {
    size_t v1 = T.wall(i).vertex1()->index();
    size_t v2 = T.wall(i).vertex2()->index();
    assert( vertexData[v1].size()==dimension && vertexData[v2].size()==dimension);
    //Calculate shared factors
    double distance=0.0;
    for( size_t d=0 ; d<dimension ; d++ )
      distance += (vertexData[v1][d]-vertexData[v2][d])*
				(vertexData[v1][d]-vertexData[v2][d]);
    distance = std::sqrt(distance);
    double wallLength=wallData[i][wallLengthIndex];
    double coeff = ((1.0/wallLength)-(1.0/distance));
    if( distance <= 0.0 && wallLength <=0.0 ) {
      //std::cerr << i << " - " << wallLength << " " << distance << std::endl;
      coeff = 0.0;
    }
		//Calculate positional factor
		double maxDistance = 0.0;
		for (size_t d=0; d<dimension; ++d) {
			double pos = 0.5*(vertexData[v1][sI]+vertexData[v2][sI]);
			maxDistance += (maxPos[d]-pos)*(maxPos[d]-pos);
		}
		maxDistance = std::sqrt(maxDistance);
		double posFactor = std::pow(maxDistance,parameter(3));
		posFactor = parameter(0) + parameter(1)*posFactor/(Kpow_+posFactor); 
		coeff *= posFactor;

    if( distance>wallLength )
      coeff *=parameter(4);
		
		//Save force in wall variable if appropriate
		if( numVariableIndexLevel()>1 )
			wallData[i][variableIndex(1,0)] = coeff*distance;
    
    //Update both vertices for each dimension
    for(size_t d=0 ; d<dimension ; d++ ) {
      double div = (vertexData[v1][d]-vertexData[v2][d])*coeff;
      vertexDerivs[v1][d] -= div;
      vertexDerivs[v2][d] += div;
    }
  }
}

VertexFromWallSpringMT::
VertexFromWallSpringMT(std::vector<double> &paraValue, 
													std::vector< std::vector<size_t> > 
													&indValue ) 
{  
  // Do some checks on the parameters and variable indeces
  if( paraValue.size()!=3 ) {
    std::cerr << "VertexFromWallSpringMT::"
	      << "VertexFromWallSpringMT() "
	      << "Uses three parameters K_force^min K_force^max "
	      << "frac_adhesion.\n";
    exit(0);
  }
  if( indValue.size() < 1 || indValue.size() > 2 
      || indValue[0].size() != 2 
      || (indValue.size()==2 && indValue[1].size() != 1) ) {
    std::cerr << "VertexFromWallSpringMT::"
	      << "VertexFromWallSpringMT() "
	      << "Wall length index and cell MT direction start index"
	      << "given at first level,"
	      << " and optionally wall variable save index in second.\n";
    exit(0);
  }
  //Set the variable values
  setId("VertexFromWallSpringMT");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "K_force^min";
  tmp[1] = "K_force^max";
  tmp[2] = "frac_adh";
  setParameterId( tmp );	
}

void VertexFromWallSpringMT::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{  
  //Do the update for each wall
  size_t numWalls = T.numWall();
  size_t wallLengthIndex = variableIndex(0,0);
	size_t directionIndex = variableIndex(0,1);
	size_t dimension = vertexData[0].size();
  
  for( size_t i=0 ; i<numWalls ; ++i ) {
    size_t v1 = T.wall(i).vertex1()->index();
    size_t v2 = T.wall(i).vertex2()->index();
    assert( vertexData[v2].size()==dimension );
    //Calculate shared factors
    double distance=0.0,c1Norm=0.0,c2Norm=0.0;
		std::vector<double> n_w(dimension),n_c1(dimension),n_c2(dimension);
    for( size_t d=0 ; d<dimension ; d++ ) {
			n_w[d] = vertexData[v2][d]-vertexData[v1][d];
			distance += n_w[d]*n_w[d];
			if( T.wall(i).cell1() != T.background() && 
					cellData[T.wall(i).cell1()->index()][directionIndex+dimension]>0.5 ) {
				n_c1[d] = cellData[T.wall(i).cell1()->index()][directionIndex+d];
				c1Norm += n_c1[d]*n_c1[d];
			}
			if( T.wall(i).cell2() != T.background() &&
					cellData[T.wall(i).cell2()->index()][directionIndex+dimension]>0.5 ) {
				n_c2[d] = cellData[T.wall(i).cell2()->index()][directionIndex+d];
				c2Norm += n_c2[d]*n_c2[d];			
			}
		}
    distance = std::sqrt( distance );
		c1Norm = std::sqrt( c1Norm );
		c2Norm = std::sqrt( c2Norm );
		double c1Fac=0.0,c2Fac=0.0;
		if( T.wall(i).cell1() != T.background() &&
				cellData[T.wall(i).cell1()->index()][directionIndex+dimension]>0.5 ) {
			for( size_t d=0 ; d<dimension ; d++ )		
				c1Fac += n_c1[d]*n_w[d];
			//c1Fac = std::fabs(c1Fac)/(c1Norm*distance);
			c1Fac /= (c1Norm*distance);
			c1Fac = c1Fac*c1Fac;
		}
		else
			c1Fac = 1.0;//0.5
		if( T.wall(i).cell2() != T.background() &&
				cellData[T.wall(i).cell2()->index()][directionIndex+dimension]>0.5 ) {
			for( size_t d=0 ; d<dimension ; d++ )		
				c2Fac += n_c2[d]*n_w[d];
			//c2Fac = std::fabs(c2Fac)/(c2Norm*distance);
			c2Fac /= (c2Norm*distance);
			c2Fac = c2Fac*c2Fac;
		}
		else
			c2Fac = 1.0;//0.5
		
    double wallLength=wallData[i][wallLengthIndex];
    double coeff = (parameter(0)+parameter(1)*(2.0-c1Fac-c2Fac))*
			((1.0/wallLength)-(1.0/distance));
    if( distance <= 0.0 && wallLength <=0.0 ) {
      //std::cerr << i << " - " << wallLength << " " << distance << std::endl;
      coeff = 0.0;
    }
    if( distance>wallLength )
      coeff *=parameter(2);
    
		//Save force in wall variable if appropriate
		if( numVariableIndexLevel()>1 )
			wallData[i][variableIndex(1,0)] = coeff*distance;
		
    //Update both vertices for each dimension
    for(size_t d=0 ; d<dimension ; d++ ) {
      double div = (vertexData[v1][d]-vertexData[v2][d])*coeff;
      vertexDerivs[v1][d] -= div;
      vertexDerivs[v2][d] += div;
    }
  }
}

VertexFromWallSpringMTSpatial::
VertexFromWallSpringMTSpatial(std::vector<double> &paraValue, 
			      std::vector< std::vector<size_t> > 
			      &indValue ) 
{  
  // Do some checks on the parameters and variable indeces
  if( paraValue.size()!=6 ) {
    std::cerr << "VertexFromWallSpringMTSpatial::"
							<< "VertexFromWallSpringMTSpatial() "
							<< "Uses six parameters k_0 k_minMT k_maxMT K_spatial n_spatial frac_adhesion."
							<< std::endl;
    exit(0);
  }
  if( indValue.size() < 1 || indValue.size() > 2 
			|| indValue[0].size() != 3 
			|| (indValue.size()==2 && indValue[1].size() != 1) ) {
    std::cerr << "VertexFromWallSpringMTSpatial::"
							<< "VertexFromWallSpringMTSpatial() "
							<< "Wall length index, MT direction start index and spatial max index given in"
							<< " first level, and optionally wall variable save index in second."
							<< std::endl;
    exit(0);
  }
	
  // Set the variable values
  setId("VertexFromWallSpringMTSpatial");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  // Set the parameter identities
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "k_0";
  tmp[1] = "k_minMT";
  tmp[2] = "k_maxMT";
  tmp[3] = "K_spatial";
  tmp[4] = "n_spatial";
  tmp[5] = "frac_adh";
  setParameterId( tmp );
	Kpow_=std::pow(parameter(3),parameter(4));
}

void VertexFromWallSpringMTSpatial::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) {
  
  // Do the update for each wall
  size_t numWalls = T.numWall();
  size_t wallLengthIndex = variableIndex(0,0);
	size_t directionIndex = variableIndex(0,1);
	size_t dimension = vertexData[0].size();
  
	// Initiate positional factor
	size_t numVertices = T.numVertex();
	size_t sI=variableIndex(0,2);
	assert (sI<vertexData[0].size());
	double max = vertexData[0][sI];
	size_t maxI=0;
	for (size_t i=1; i<numVertices; ++i)
		if (vertexData[i][sI]>max) {
			max=vertexData[i][sI];
			maxI=i;
		}
	std::vector<double> maxPos(dimension);
	for (size_t d=0; d<dimension; ++d)
		maxPos[d] = vertexData[maxI][d];

	
	// Calculate update from each wall
  for( size_t i=0 ; i<numWalls ; ++i ) {
    size_t v1 = T.wall(i).vertex1()->index();
    size_t v2 = T.wall(i).vertex2()->index();
    assert( vertexData[v2].size()==dimension );
		
    //Calculate shared factors
    double distance=0.0,c1Norm=0.0,c2Norm=0.0;
		std::vector<double> n_w(dimension),n_c1(dimension),n_c2(dimension);
    for( size_t d=0 ; d<dimension ; d++ ) {
			n_w[d] = vertexData[v2][d]-vertexData[v1][d];
			distance += n_w[d]*n_w[d];
			if( T.wall(i).cell1() != T.background() && 
					cellData[T.wall(i).cell1()->index()][directionIndex+dimension]>0.5 ) {
				n_c1[d] = cellData[T.wall(i).cell1()->index()][directionIndex+d];
				c1Norm += n_c1[d]*n_c1[d];
			}
			if( T.wall(i).cell2() != T.background() &&
					cellData[T.wall(i).cell2()->index()][directionIndex+dimension]>0.5 ) {
				n_c2[d] = cellData[T.wall(i).cell2()->index()][directionIndex+d];
				c2Norm += n_c2[d]*n_c2[d];			
			}
		}
    distance = std::sqrt( distance );
		c1Norm = std::sqrt( c1Norm );
		c2Norm = std::sqrt( c2Norm );

		// Calculate MT factors
		double c1Fac=0.0,c2Fac=0.0;
		if( T.wall(i).cell1() != T.background() &&
				cellData[T.wall(i).cell1()->index()][directionIndex+dimension]>0.5 ) {
			for( size_t d=0 ; d<dimension ; d++ )		
				c1Fac += n_c1[d]*n_w[d];
			c1Fac = std::fabs(c1Fac)/(c1Norm*distance);
		}
		else
			c1Fac = 0.5;//1.0;
		if( T.wall(i).cell2() != T.background() &&
				cellData[T.wall(i).cell2()->index()][directionIndex+dimension]>0.5 ) {
			for( size_t d=0 ; d<dimension ; d++ )		
				c2Fac += n_c2[d]*n_w[d];
			c2Fac = std::fabs(c2Fac)/(c2Norm*distance);
		}
		else
			c2Fac = 0.5;//1.0;
		double mtFactor = (parameter(1)+parameter(2)*(2.0-c1Fac-c2Fac));
		
		//Calculate positional factor
		double maxDistance=0.0;
		for (size_t d=0; d<dimension; ++d) {
			double pos = 0.5*(vertexData[v1][d]+vertexData[v2][d]);
			maxDistance += (maxPos[d]-pos)*(maxPos[d]-pos);
		}
		maxDistance = std::sqrt(distance);
		double posFactor = std::pow(maxDistance,parameter(4));
		
		double wallLength=wallData[i][wallLengthIndex];
		double coeff = ((1.0/wallLength)-(1.0/distance));
    if( distance <= 0.0 && wallLength <=0.0 ) {
      //std::cerr << i << " - " << wallLength << " " << distance << std::endl;
      coeff = 0.0;
    }
		// Multiply with positional and MT factor
		coeff *= parameter(0) + mtFactor*posFactor/(Kpow_+posFactor); 
		
		// Multiply with compression factor
    if( distance>wallLength )
      coeff *=parameter(5);
		
		//Save force in wall variable if appropriate
		if( numVariableIndexLevel()>1 )
			wallData[i][variableIndex(1,0)] = coeff*distance;
    
    //Update both vertices for each dimension
    for(size_t d=0 ; d<dimension ; d++ ) {
      double div = (vertexData[v1][d]-vertexData[v2][d])*coeff;
      vertexDerivs[v1][d] -= div;
      vertexDerivs[v2][d] += div;
    }
  }
}

VertexFromWallSpringMTHistory::
VertexFromWallSpringMTHistory(std::vector<double> &paraValue, 
															std::vector< std::vector<size_t> > 
															&indValue ) 
{  
  // Do some checks on the parameters and variable indeces
  if( paraValue.size()!=4 ) {
    std::cerr << "VertexFromWallSpringMTHistory::"
							<< "VertexFromWallSpringMTHistory() "
							<< "Uses four parameters K_force^min K_force^max "
							<< "frac_adhesion and k_rate.\n";
    exit(0);
  }
  if( indValue.size() < 2 || indValue.size()>3 || indValue[0].size() != 2 
			|| indValue[1].size() != 1 ||
			(indValue.size()==3 && indValue[2].size() != 1) ) {
		std::cerr << "VertexFromWallSpringMTHistory::"
							<< "VertexFromWallSpringMTHistory() "
							<< "Wall length index and cell MT direction start index"
							<< "given at first level,"
							<< " wall spring constant index at second, and optinally "
							<< " force saving index at level three." << std::endl;
    exit(0);
  }
  //Set the variable values
  setId("VertexFromWallSpringMTHistory");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "K_force^min";
  tmp[1] = "K_force^max";
  tmp[2] = "frac_adh";
	tmp[3] = "k_rate";
  setParameterId( tmp );	
}

void VertexFromWallSpringMTHistory::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{  
  //Do the update for each wall
  size_t numWalls = T.numWall();
  size_t wallLengthIndex = variableIndex(0,0);
	size_t directionIndex = variableIndex(0,1);
	size_t dimension = vertexData[0].size();
  
  for( size_t i=0 ; i<numWalls ; ++i ) {
    size_t v1 = T.wall(i).vertex1()->index();
    size_t v2 = T.wall(i).vertex2()->index();
    assert( vertexData[v2].size()==dimension );
    //Calculate shared factors
    double distance=0.0,c1Norm=0.0,c2Norm=0.0;
		std::vector<double> n_w(dimension),n_c1(dimension),n_c2(dimension);
    for( size_t d=0 ; d<dimension ; d++ ) {
			n_w[d] = vertexData[v2][d]-vertexData[v1][d];
			distance += n_w[d]*n_w[d];
			if( T.wall(i).cell1() != T.background() && 
					cellData[T.wall(i).cell1()->index()][directionIndex+dimension]>0.5 ) {
				n_c1[d] = cellData[T.wall(i).cell1()->index()][directionIndex+d];
				c1Norm += n_c1[d]*n_c1[d];
			}
			if( T.wall(i).cell2() != T.background() &&
					cellData[T.wall(i).cell2()->index()][directionIndex+dimension]>0.5 ) {
				n_c2[d] = cellData[T.wall(i).cell2()->index()][directionIndex+d];
				c2Norm += n_c2[d]*n_c2[d];			
			}
		}
    distance = std::sqrt( distance );
		c1Norm = std::sqrt( c1Norm );
		c2Norm = std::sqrt( c2Norm );
		double c1Fac=0.0,c2Fac=0.0;
		if( T.wall(i).cell1() != T.background() &&
				cellData[T.wall(i).cell1()->index()][directionIndex+dimension]>0.5 ) {
			for( size_t d=0 ; d<dimension ; d++ )		
				c1Fac += n_c1[d]*n_w[d];
			c1Fac = std::fabs(c1Fac)/(c1Norm*distance);
		}
		else
			c1Fac = 0.5;//1.0;
		if( T.wall(i).cell2() != T.background() &&
				cellData[T.wall(i).cell2()->index()][directionIndex+dimension]>0.5 ) {
			for( size_t d=0 ; d<dimension ; d++ )		
				c2Fac += n_c2[d]*n_w[d];
			c2Fac = std::fabs(c2Fac)/(c2Norm*distance);
		}
		else
			c2Fac = 0.5;//1.0;
		
    double wallLength=wallData[i][wallLengthIndex];
		double springConstant = (parameter(0)+parameter(1)*(2.0-c1Fac-c2Fac));
    double coeff = wallData[i][variableIndex(1,0)]*
			((1.0/wallLength)-(1.0/distance));
    if( distance <= 0.0 && wallLength <=0.0 ) {
      //std::cerr << i << " - " << wallLength << " " << distance << std::endl;
      coeff = 0.0;
    }
    if( distance>wallLength )
      coeff *=parameter(2);
    
		//Save force in wall variable if appropriate
		if( numVariableIndexLevel()>2 )
			wallData[i][variableIndex(2,0)] = coeff*distance;
		
    // Update both vertices for each dimension
    for(size_t d=0 ; d<dimension ; d++ ) {
      double div = (vertexData[v1][d]-vertexData[v2][d])*coeff;
      vertexDerivs[v1][d] -= div;
      vertexDerivs[v2][d] += div;
    }
		// Update the spring constant
		wallDerivs[i][variableIndex(1,0)] += parameter(3)*
			(springConstant-wallData[i][variableIndex(1,0)]); 
  }
}

void VertexFromWallSpringMTHistory::
initiate(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs)
{
  //Do the initiation for each wall
  size_t numWalls = T.numWall();
  size_t directionIndex = variableIndex(0,1);
  size_t dimension = T.vertex(0).numPosition();
  
  for( size_t i=0 ; i<numWalls ; ++i ) {
    Vertex* v1 = T.wall(i).vertex1();
    Vertex* v2 = T.wall(i).vertex2();
    
    //Calculate shared factors
    double distance=0.0,c1Norm=0.0,c2Norm=0.0;
    std::vector<double> n_w(dimension),n_c1(dimension),n_c2(dimension);
    for( size_t d=0 ; d<dimension ; d++ ) {
      n_w[d] = v2->position(d)-v1->position(d);
      distance += n_w[d]*n_w[d];
      if( T.wall(i).cell1() != T.background() && 
	  T.wall(i).cell1()->variable(directionIndex+dimension)>0.5 ) {
	n_c1[d] = T.wall(i).cell1()->variable(directionIndex+d);
	c1Norm += n_c1[d]*n_c1[d];
      }
      if( T.wall(i).cell2() != T.background() &&
	  T.wall(i).cell2()->variable(directionIndex+dimension)>0.5 ) {
	n_c2[d] = T.wall(i).cell2()->variable(directionIndex+d);
	c2Norm += n_c2[d]*n_c2[d];			
      }
    }
    distance = std::sqrt( distance );
    c1Norm = std::sqrt( c1Norm );
    c2Norm = std::sqrt( c2Norm );
    double c1Fac=0.0,c2Fac=0.0;
    if( T.wall(i).cell1() != T.background() &&
	T.wall(i).cell1()->variable(directionIndex+dimension)>0.5 ) {
      for( size_t d=0 ; d<dimension ; d++ )		
	c1Fac += n_c1[d]*n_w[d];
      c1Fac = std::fabs(c1Fac)/(c1Norm*distance);
    }
    else
      c1Fac = 0.5;//1.0;
    if( T.wall(i).cell2() != T.background() &&
	T.wall(i).cell2()->variable(directionIndex+dimension)>0.5 ) {
      for( size_t d=0 ; d<dimension ; d++ )		
	c2Fac += n_c2[d]*n_w[d];
      c2Fac = std::fabs(c2Fac)/(c2Norm*distance);
    }
    else
      c2Fac = 0.5;//1.0;
    
    wallData[i][variableIndex(1,0)] = parameter(0)+parameter(1) *
      (2.0-c1Fac-c2Fac);
  }
}

VertexFromEpidermalWallSpring::
VertexFromEpidermalWallSpring(std::vector<double> &paraValue, 
					std::vector< std::vector<size_t> > 
					&indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //////////////////////////////////////////////////////////////////////
  if( paraValue.size()!=2 ) {
    std::cerr << "VertexFromEpidermalWallSpring::"
	      << "VertexFromEpidermalWallSpring() "
	      << "Uses two parameters K_force frac_adhesion.\n";
    exit(0);
  }

  if( indValue.size() < 1 || indValue.size() > 2 
			|| indValue[0].size() != 1 
			|| (indValue.size()==2 && indValue[1].size() != 1) ) {
    std::cerr << "VertexFromEpidermalWallSpring::"
							<< "VertexFromEpidermalWallSpring() "
							<< "Wall length index given in first level,"
							<< " and optionally wall variable save index in second.\n";
    exit(0);
  }
  //Set the variable values
  //////////////////////////////////////////////////////////////////////
  setId("VertexFromEpidermalWallSpring");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //////////////////////////////////////////////////////////////////////
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "K_force";
  tmp[1] = "frac_adh";
  setParameterId( tmp );
}

//! Derivative contribution for asymmetric wall springs on vertices
/*! 
*/
void VertexFromEpidermalWallSpring::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) {
  
  //Do the update for each wall
  size_t numWalls = T.numWall();
  size_t wallLengthIndex = variableIndex(0,0);
  
  for( size_t i=0 ; i<numWalls ; ++i ) {
    if( !( T.wall(i).cell1() != T.background() &&
					 T.wall(i).cell2() != T.background() ) ) {
      size_t v1 = T.wall(i).vertex1()->index();
      size_t v2 = T.wall(i).vertex2()->index();
      size_t dimension = vertexData[v1].size();
      assert( vertexData[v2].size()==dimension );
      //Calculate shared factors
      double distance=0.0;
      for( size_t d=0 ; d<dimension ; d++ )
				distance += (vertexData[v1][d]-vertexData[v2][d])*
					(vertexData[v1][d]-vertexData[v2][d]);
      distance = std::sqrt(distance);
      double wallLength=wallData[i][wallLengthIndex];
      double coeff = parameter(0)*((1.0/wallLength)-(1.0/distance));
      if( distance <= 0.0 && wallLength <=0.0 ) {
				//std::cerr << i << " - " << wallLength << " " << distance << std::endl;
				coeff = 0.0;
      }
      if( distance>wallLength )
				coeff *=parameter(1);

			//Save force in wall variable if appropriate
			if( numVariableIndexLevel()>1 )
				wallData[i][variableIndex(1,0)] = coeff*distance;
      
      //Update both vertices for each dimension
      for(size_t d=0 ; d<dimension ; d++ ) {
				double div = (vertexData[v1][d]-vertexData[v2][d])*coeff;
				vertexDerivs[v1][d] -= div;
				vertexDerivs[v2][d] += div;
      }
    }
		else if( numVariableIndexLevel()>1 )
			wallData[i][variableIndex(1,0)] = 0.0;
  }
}

VertexFromEpidermalCellWallSpring::
VertexFromEpidermalCellWallSpring(std::vector<double> &paraValue, 
																						std::vector< std::vector<size_t> > 
																						&indValue ) 
{  
  //Do some checks on the parameters and variable indeces
  //////////////////////////////////////////////////////////////////////
  if( paraValue.size()!=2 ) {
    std::cerr << "VertexFromEpidermalCellWallSpring::"
							<< "VertexFromEpidermalCellWallSpring() "
							<< "Uses two parameters K_force frac_adhesion.\n";
    exit(0);
  }
  if( indValue.size() < 1 || indValue.size() > 2 
			|| indValue[0].size() != 1 
			|| (indValue.size()==2 && indValue[1].size() != 1) ) {
    std::cerr << "VertexFromEpidermalCellWallSpring::"
							<< "VertexFromEpidermalCellWallSpring() "
							<< "Wall length index given in first level,"
							<< " and optionally wall variable save index in second.\n";
    exit(0);
  }
  //Set the variable values
  //////////////////////////////////////////////////////////////////////
  setId("VertexFromEpidermalCellWallSpring");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //////////////////////////////////////////////////////////////////////
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "K_force";
  tmp[1] = "frac_adh";
  setParameterId( tmp );
}

void VertexFromEpidermalCellWallSpring::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{
  
  //Do the update for each wall
  size_t numWalls = T.numWall();
  size_t wallLengthIndex = variableIndex(0,0);
  
  for( size_t i=0 ; i<numWalls ; ++i ) {
    if( T.wall(i).cell1() == T.background() ||
				T.wall(i).cell1()->isNeighbor(T.background()) ||
				T.wall(i).cell2() == T.background() || 
				T.wall(i).cell2()->isNeighbor(T.background()) ) {
      size_t v1 = T.wall(i).vertex1()->index();
      size_t v2 = T.wall(i).vertex2()->index();
      size_t dimension = vertexData[v1].size();
      assert( vertexData[v2].size()==dimension );
      //Calculate shared factors
      double distance=0.0;
      for( size_t d=0 ; d<dimension ; d++ )
				distance += (vertexData[v1][d]-vertexData[v2][d])*
					(vertexData[v1][d]-vertexData[v2][d]);
      distance = std::sqrt(distance);
      double wallLength=wallData[i][wallLengthIndex];
      double coeff = parameter(0)*((1.0/wallLength)-(1.0/distance));
      if( distance <= 0.0 && wallLength <=0.0 ) {
				//std::cerr << i << " - " << wallLength << " " << distance << std::endl;
				coeff = 0.0;
      }
      if( distance>wallLength )
				coeff *=parameter(1);
      
			//Save force in wall variable if appropriate
			if( numVariableIndexLevel()>1 )
				wallData[i][variableIndex(1,0)] = coeff*distance;
			
      //Update both vertices for each dimension
      for(size_t d=0 ; d<dimension ; d++ ) {
				double div = (vertexData[v1][d]-vertexData[v2][d])*coeff;
				vertexDerivs[v1][d] -= div;
				vertexDerivs[v2][d] += div;
      }
    }
		else if( numVariableIndexLevel()>1 )
			wallData[i][variableIndex(1,0)] = 0.0;
  }
}

VertexFromWallSpringExperimental::
VertexFromWallSpringExperimental(std::vector<double> &paraValue,
																 std::vector< std::vector<size_t> > &indValue)
{
	if (paraValue.size() != 1) {
		std::cerr << "VertexFromWallSpringExperimental::VertexFromWallSpringExperimental() "
				<< "Uses one parameter: k" << std::endl;
		exit(EXIT_FAILURE);
	}

     if (indValue.size() == 0 || indValue.size() > 2 || indValue[0].size() != 1) {
		std::cerr << "VertexFromWallSpringExperimental::VertexFromWallSpringExperimental() "
                    << "Wall length index given.\n";
          exit(EXIT_FAILURE);
     }

     if (indValue.size() == 1 && indValue[0].size() != 1) {
		std::cerr << "VertexFromWallSpringExperimental::VertexFromWallSpringExperimental() -"
                    << "Second level of indices gives index for storage of force.\n";
          exit(EXIT_FAILURE);
     }

	setId("VertexFromWallSpringExperimental");
	setParameter(paraValue);  
	setVariableIndex(indValue);
	
	std::vector<std::string> tmp(numParameter());
	tmp[0] = "k";
	setParameterId(tmp);
}

void VertexFromWallSpringExperimental::
derivs(Tissue &T,
			 DataMatrix &cellData,
			 DataMatrix &wallData,
			 DataMatrix &vertexData,
			 DataMatrix &cellDerivs,
			 DataMatrix &wallDerivs,
			 DataMatrix &vertexDerivs)
{
	for (size_t i = 0; i < T.numWall(); ++i) {
		size_t vertex1Index = T.wall(i).vertex1()->index();
		size_t vertex2Index = T.wall(i).vertex2()->index();
		size_t dimensions = vertexData[vertex1Index].size();
		
		double distance = 0.0;
		for (size_t d = 0; d < dimensions; ++d) {
			distance += (vertexData[vertex1Index][d] - vertexData[vertex2Index][d])
				* (vertexData[vertex1Index][d] - vertexData[vertex2Index][d]);
		}
		distance = std::sqrt(distance);

 		for (size_t d = 0; d < dimensions; ++d) {
			double dx1dt = parameter(0) * (vertexData[vertex2Index][d] - vertexData[vertex1Index][d])
				* (1.0/wallData[i][variableIndex(0, 0)] - 1.0/distance);

			vertexDerivs[vertex1Index][d] += dx1dt;
			vertexDerivs[vertex2Index][d] -= dx1dt;
		}
		if (numVariableIndexLevel() == 2) {
			wallData[T.wall(i).index()][variableIndex(1, 0)] = (parameter(0) / wallData[T.wall(i).index()][variableIndex(0, 0)]);
			wallData[T.wall(i).index()][variableIndex(1, 0)] *= (distance - wallData[T.wall(i).index()][variableIndex(0, 0)]);
		}
 	}
}

VertexFromWallSpringConcentrationHill::
VertexFromWallSpringConcentrationHill(std::vector<double> &paraValue, 
															 std::vector< std::vector<size_t> > 
															 &indValue ) 
{  
  // Do some checks on the parameters and variable indeces
  if( paraValue.size()!=5 ) {
    std::cerr << "VertexFromWallSpringConcentrationHill::"
	      << "VertexFromWallSpringConcentrationHill() "
	      << "Uses five parameters K_min, K_max, K_Hill, n_Hill and frac_adhesion.\n";
    exit(0);
  }
  if (indValue.size() < 1 || indValue.size() > 2 
			|| indValue[0].size() != 2 
			|| (indValue.size()==2 && indValue[1].size() != 1) ) {
    std::cerr << "VertexFromWallSpringConcentrationHill::"
							<< "VertexFromWallSpringConcentrationHill() "
							<< "Wall length index and cell concentration index given in first level,"
							<< " and optionally wall force save index in second.\n";
    exit(0);
  }
	
  // Set the variable values
  setId("VertexFromWallSpringConcentrationHill");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  // Set the parameter identities
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "K_min";
  tmp[1] = "K_max";
  tmp[2] = "K_Hill";
  tmp[3] = "n_Hill";
  tmp[4] = "frac_adh";
  setParameterId( tmp );
}

//! Derivative contribution for asymmetric wall springs on vertices
/*! 
*/
void VertexFromWallSpringConcentrationHill::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) {
  
  //Do the update for each wall
  size_t numWalls = T.numWall();
  size_t wallLengthIndex = variableIndex(0,0);
  size_t concentrationIndex = variableIndex(0,1);
  size_t dimension = vertexData[0].size();
  
  for( size_t i=0 ; i<numWalls ; ++i ) {
    size_t v1 = T.wall(i).vertex1()->index();
    size_t v2 = T.wall(i).vertex2()->index();
    assert( vertexData[v2].size()==dimension );
    //Calculate shared factors
    double distance=0.0;
    std::vector<double> n_w(dimension),n_c1(dimension),n_c2(dimension);
    for( size_t d=0 ; d<dimension ; d++ ) {
      n_w[d] = vertexData[v2][d]-vertexData[v1][d];
      distance += n_w[d]*n_w[d];
    }
    distance = std::sqrt( distance );
    double c1Fac=0.0,c2Fac=0.0,KPow = std::pow(parameter(2),parameter(3));
    if( T.wall(i).cell1() != T.background() ) {
      double conc = cellData[T.wall(i).cell1()->index()][concentrationIndex];
      c1Fac = KPow/(KPow+std::pow(conc,parameter(3)));
    }
    if( T.wall(i).cell2() != T.background() ) {
      double conc = cellData[T.wall(i).cell2()->index()][concentrationIndex];
      c2Fac = KPow/(KPow+std::pow(conc,parameter(3)));
    }
    
    double wallLength=wallData[i][wallLengthIndex];
    double coeff = (parameter(0)+parameter(1)*(c1Fac+c2Fac))*
      ((1.0/wallLength)-(1.0/distance));
    if( distance <= 0.0 && wallLength <=0.0 ) {
      //std::cerr << i << " - " << wallLength << " " << distance << std::endl;
      coeff = 0.0;
    }
    if( distance>wallLength )
      coeff *=parameter(4);
    
    //Save force in wall variable if appropriate
    if( numVariableIndexLevel()>1 )
      wallData[i][variableIndex(1,0)] = coeff*distance;
    
    //Update both vertices for each dimension
    for(size_t d=0 ; d<dimension ; d++ ) {
      double div = (vertexData[v1][d]-vertexData[v2][d])*coeff;
      vertexDerivs[v1][d] -= div;
      vertexDerivs[v2][d] += div;
    }
  }
}

VertexFromWallSpringMTConcentrationHill::
VertexFromWallSpringMTConcentrationHill(std::vector<double> &paraValue, 
					std::vector< std::vector<size_t> > 
					&indValue ) 
{  
  // Do some checks on the parameters and variable indeces
  if( paraValue.size()!=6 ) {
    std::cerr << "VertexFromWallSpringMTConcentrationHill::"
	      << "VertexFromWallSpringMTConcentrationHill() "
	      << "Uses six parameters K_0 frac_MT frac_conc K_Hill n_Hill "
	      << "frac_adhesion.\n";
    exit(0);
  }
  if( indValue.size() < 2 || indValue.size() > 3 
      || indValue[0].size() != 2
      || indValue[1].size() != 1
      || (indValue.size()==3 && indValue[2].size() != 1) ) {
    std::cerr << "VertexFromWallSpringMTConcentrationHill::"
	      << "VertexFromWallSpringMTConcentrationHill() "
	      << "Wall length index and cell MT direction start index"
	      << "given at first level, conc index at second,"
	      << " and optionally wall variable save index in third.\n";
    exit(0);
  }
  //Set the variable values
  setId("VertexFromWallSpringMTConcentrationHill");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "K_0";
  tmp[1] = "frac_MT";
  tmp[2] = "frac_conc";
  tmp[3] = "K_Hill";
  tmp[4] = "n_Hill";
  tmp[5] = "frac_adh";
  setParameterId( tmp );	
}

void VertexFromWallSpringMTConcentrationHill::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{  
  //Do the update for each wall
  size_t numWalls = T.numWall();
  size_t wallLengthIndex = variableIndex(0,0);
  size_t directionIndex = variableIndex(0,1);
  size_t dimension = vertexData[0].size();
  double KPow = std::pow(parameter(3),parameter(4));
  
  for( size_t i=0 ; i<numWalls ; ++i ) {
    size_t v1 = T.wall(i).vertex1()->index();
    size_t v2 = T.wall(i).vertex2()->index();
    assert( vertexData[v1].size()==dimension );
    assert( vertexData[v2].size()==dimension );
    //Calculate shared factors
    double distance=0.0,c1Norm=0.0,c2Norm=0.0;
    std::vector<double> n_w(dimension),n_c1(dimension),n_c2(dimension);
    for( size_t d=0 ; d<dimension ; d++ ) {
      n_w[d] = vertexData[v2][d]-vertexData[v1][d];
      distance += n_w[d]*n_w[d];
      if( T.wall(i).cell1() != T.background() && 
	  cellData[T.wall(i).cell1()->index()][directionIndex+dimension]>0.5 ) {
	n_c1[d] = cellData[T.wall(i).cell1()->index()][directionIndex+d];
	c1Norm += n_c1[d]*n_c1[d];
      }
      if( T.wall(i).cell2() != T.background() &&
	  cellData[T.wall(i).cell2()->index()][directionIndex+dimension]>0.5 ) {
	n_c2[d] = cellData[T.wall(i).cell2()->index()][directionIndex+d];
	c2Norm += n_c2[d]*n_c2[d];			
      }
    }
    distance = std::sqrt( distance );
    c1Norm = std::sqrt( c1Norm );
    c2Norm = std::sqrt( c2Norm );
    double c1Fac=0.0,c2Fac=0.0,c1FacConc=0.0,c2FacConc=0.0;
    if( T.wall(i).cell1() != T.background() &&
	cellData[T.wall(i).cell1()->index()][directionIndex+dimension]>0.5 ) {
      for( size_t d=0 ; d<dimension ; d++ )		
	c1Fac += n_c1[d]*n_w[d];
      c1Fac /= (c1Norm*distance);
      c1Fac = c1Fac*c1Fac;
      double conc=cellData[T.wall(i).cell1()->index()][variableIndex(1,0)];
      c1FacConc = KPow/(KPow+std::pow(conc,parameter(4)));
    }
    else
      c1Fac = 0.5;//1.0;
    if( T.wall(i).cell2() != T.background() &&
	cellData[T.wall(i).cell2()->index()][directionIndex+dimension]>0.5 ) {
      for( size_t d=0 ; d<dimension ; d++ )		
	c2Fac += n_c2[d]*n_w[d];
      c2Fac /= (c2Norm*distance);
      c2Fac = c2Fac*c2Fac;
      double conc=cellData[T.wall(i).cell2()->index()][variableIndex(1,0)];
      c2FacConc = KPow/(KPow+std::pow(conc,parameter(4)));
    }
    else
      c2Fac = 0.5;//1.0;
    
    double wallLength=wallData[i][wallLengthIndex];
    double coeff = parameter(0)*((1.0-parameter(1))+parameter(1)*0.5*(2.0-c1Fac-c2Fac))*
      ((1.0-parameter(2))+parameter(2)*0.5*(c1FacConc+c2FacConc))*
      ((1.0/wallLength)-(1.0/distance));
    if( distance <= 0.0 && wallLength <=0.0 ) {
      //std::cerr << i << " - " << wallLength << " " << distance << std::endl;
      coeff = 0.0;
    }
    if( distance>wallLength )
      coeff *=parameter(5);
    
    //Save force in wall variable if appropriate
    if( numVariableIndexLevel()>2 )
      wallData[i][variableIndex(2,0)] = coeff*distance;
    
    //Update both vertices for each dimension
    for(size_t d=0 ; d<dimension ; d++ ) {
      double div = (vertexData[v1][d]-vertexData[v2][d])*coeff;
      vertexDerivs[v1][d] -= div;
      vertexDerivs[v2][d] += div;
    }
  }
}

VertexFromDoubleWallSpringMTConcentrationHill::
VertexFromDoubleWallSpringMTConcentrationHill(std::vector<double> &paraValue,											 std::vector< std::vector<size_t> > &indValue ) 
{  
  // Do some checks on the parameters and variable indeces
  if( paraValue.size()!=6 ) {
    std::cerr << "VertexFromDoubleWallSpringMTConcentrationHill::"
	      << "VertexFromDoubleWallSpringMTConcentrationHill() "
	      << "Uses six parameters K_0 frac_MT frac_conc K_Hill n_Hill "
	      << "frac_adhesion.\n";
    exit(0);
  }
  if( indValue.size() < 2 || indValue.size() > 3 
      || indValue[0].size() != 2
      || indValue[1].size() != 1
      || (indValue.size()==3 && indValue[2].size() != 1 && indValue[2].size() != 3) ) {
    std::cerr << "VertexFromDoubleWallSpringMTConcentrationHill::"
	      << "VertexFromDoubleWallSpringMTConcentrationHill() "
	      << "Wall length index and cell MT direction start index"
	      << "given at first level, conc index at second,"
	      << " and optionally wall variable save index (F) or "
	      << "indices (F,k1,k2) in third.\n";
    exit(0);
  }
  //Set the variable values
  setId("VertexFromDoubleWallSpringMTConcentrationHill");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "K_0";
  tmp[1] = "frac_MT";
  tmp[2] = "frac_conc";
  tmp[3] = "K_Hill";
  tmp[4] = "n_Hill";
  tmp[5] = "frac_adh";
  setParameterId( tmp );	
}

void VertexFromDoubleWallSpringMTConcentrationHill::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{  
  //Do the update for each wall
  size_t numWalls = T.numWall();
  size_t wallLengthIndex = variableIndex(0,0);
	size_t directionIndex = variableIndex(0,1);
	size_t dimension = vertexData[0].size();
  double KPow = std::pow(parameter(3),parameter(4));
	
  for( size_t i=0 ; i<numWalls ; ++i ) {
    size_t v1 = T.wall(i).vertex1()->index();
    size_t v2 = T.wall(i).vertex2()->index();
    assert( vertexData[v1].size()==dimension );
    assert( vertexData[v2].size()==dimension );
    //Calculate shared factors
    double distance=0.0,c1Norm=0.0,c2Norm=0.0;
    std::vector<double> n_w(dimension),n_c1(dimension),n_c2(dimension);
    for( size_t d=0 ; d<dimension ; d++ ) {
      n_w[d] = vertexData[v2][d]-vertexData[v1][d];
      distance += n_w[d]*n_w[d];
      if( T.wall(i).cell1() != T.background() && 
	  cellData[T.wall(i).cell1()->index()][directionIndex+dimension]>0.5 ) {
	n_c1[d] = cellData[T.wall(i).cell1()->index()][directionIndex+d];
	c1Norm += n_c1[d]*n_c1[d];
      }
      if( T.wall(i).cell2() != T.background() &&
	  cellData[T.wall(i).cell2()->index()][directionIndex+dimension]>0.5 ) {
	n_c2[d] = cellData[T.wall(i).cell2()->index()][directionIndex+d];
	c2Norm += n_c2[d]*n_c2[d];			
      }
    }
    distance = std::sqrt( distance );
    c1Norm = std::sqrt( c1Norm );
    c2Norm = std::sqrt( c2Norm );
    double c1Fac=0.0,c2Fac=0.0,c1FacConc=0.0,c2FacConc=0.0;
    if( T.wall(i).cell1() != T.background() &&
	cellData[T.wall(i).cell1()->index()][directionIndex+dimension]>0.5 ) {
      for( size_t d=0 ; d<dimension ; d++ )		
	c1Fac += n_c1[d]*n_w[d];
      c1Fac /= (c1Norm*distance);
      c1Fac = c1Fac*c1Fac;
      double conc=cellData[T.wall(i).cell1()->index()][variableIndex(1,0)];
      c1FacConc = KPow/(KPow+std::pow(conc,parameter(4)));
    }
    else
      c1Fac = 0.5;//1.0;
    if( T.wall(i).cell2() != T.background() &&
	cellData[T.wall(i).cell2()->index()][directionIndex+dimension]>0.5 ) {
      for( size_t d=0 ; d<dimension ; d++ )		
	c2Fac += n_c2[d]*n_w[d];
      c2Fac /= (c2Norm*distance);
      c2Fac = c2Fac*c2Fac;
      double conc=cellData[T.wall(i).cell2()->index()][variableIndex(1,0)];
      c2FacConc = KPow/(KPow+std::pow(conc,parameter(4)));
    }
    else
      c2Fac = 0.5;//1.0;
    
    double wallLength=wallData[i][wallLengthIndex];
    double k1 = 0.5*parameter(0)*((1.0-parameter(1))+parameter(1)*(1.0-c1Fac))*
      ((1.0-parameter(2))+parameter(2)*c1FacConc);
    double k2 = 0.5*parameter(0)*((1.0-parameter(1))+parameter(1)*(1.0-c2Fac))*
      ((1.0-parameter(2))+parameter(2)*c2FacConc);
    
    double coeff = (k1+k2)*((1.0/wallLength)-(1.0/distance));
    if( distance <= 0.0 && wallLength <=0.0 ) {
      //std::cerr << i << " - " << wallLength << " " << distance << std::endl;
      coeff = 0.0;
    }
    if( distance>wallLength && parameter(5) != 1.0 ) {
      coeff *= parameter(5);
      k1 *= parameter(5);
      k2 *= parameter(5);
    }
    //Save force as well as k parameters in wall variable if appropriate
    if( numVariableIndexLevel()>2 ) {
      wallData[i][variableIndex(2,0)] = coeff*distance;
      if( numVariableIndex(2) == 3 ) {
	wallData[i][variableIndex(2,1)] = k1;
	wallData[i][variableIndex(2,2)] = k2;
      }
    }
    //Update both vertices for each dimension
    for(size_t d=0 ; d<dimension ; d++ ) {
      double div = (vertexData[v1][d]-vertexData[v2][d])*coeff;
      vertexDerivs[v1][d] -= div;
      vertexDerivs[v2][d] += div;
    }
  }
}



VertexFromExternalSpring::
VertexFromExternalSpring(std::vector<double> &paraValue, 
			 std::vector< std::vector<size_t> > &indValue ) 
{  
  // Do some checks on the parameters and variable indeces
  if( paraValue.size()!=4 ) {
    std::cerr << "VertexFromExternalSpring::"
	      << "VertexFromExternalSpring() "
	      << "Uses four parameters spring constant K, frac_adhesion, Lmaxfactor and growth_rate.\n";
    exit(0);
  }
  if( indValue.size() != 3  || indValue[1].size() !=indValue[2].size() || indValue[0].size()!=1 ) {
    std::cerr << "VertexFromExternalSpring::"
	      << "VertexFromExternalSpring() "
	      << "vertex indices come as pairs each from one of the two sets with " 
              << "the same number of elements in the second and third levels. \n";
    exit(0);
  }
  //Set the variable values
  setId("VertexFromExternalSpring");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "K_sp";
  tmp[1] = "frac_adh";
  tmp[2] = "Lmaxfactor";
  tmp[3] = "growth_rate";
  
  setParameterId( tmp );

  Npairs=indValue[1].size();
}



void VertexFromExternalSpring::initiate(Tissue &T,
					DataMatrix &cellData,
					DataMatrix &wallData,
					DataMatrix &vertexData,
					DataMatrix &cellDerivs,
					DataMatrix &wallDerivs,
					DataMatrix &vertexDerivs) { 

  //size_t Npairs=numVariableIndex(1);  
  restinglength.resize(Npairs);
  Kspring.resize(Npairs);
  
  size_t dimension=vertexData[0].size();
  
  
  //double distance=std::sqrt(distance);
  
  for (size_t i=0 ; i< Npairs ; i++){
    Kspring[i]=parameter(0);
    size_t v1=variableIndex(1,i);
    size_t v2=variableIndex(2,i);
    restinglength[i]=0;
    for( size_t d=0 ; d<dimension ; d++ ) {
      restinglength[i] += (vertexData[v2][d]-vertexData[v1][d])* (vertexData[v2][d]-vertexData[v1][d]);
    }
    
    restinglength[i]=std::sqrt(restinglength[i]);
  }
   
}




void VertexFromExternalSpring::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{  
  double fad=parameter(1);
 
  //size_t Npairs=indValue[1].size();
  size_t dimension=vertexData[0].size();

  for (size_t i=0 ; i< Npairs ; i++){
    
    size_t v1=variableIndex(1,i);
    size_t v2=variableIndex(2,i);
    
    double distance=0;
    
    for( size_t d=0 ; d<dimension ; d++ ) 
      distance += (vertexData[v2][d]-vertexData[v1][d])*(vertexData[v2][d]-vertexData[v1][d]);
    
    distance=std::sqrt(distance);

    double coeff=0;
    if ( restinglength[i]>0.0)
      coeff=Kspring[i]*((1.0/restinglength[i])-(1.0/distance));
    
    if( distance <= 0.0 && restinglength[i] <=0.0 ) 
      coeff = 0.0;
    
    if( distance>restinglength[i])
      coeff *=fad;
    //Update both vertices for each dimension
    for(size_t d=0 ; d<dimension ; d++ ) {
      double div = coeff*(vertexData[v2][d]-vertexData[v1][d]);
      if ( restinglength[i]>0){
	vertexDerivs[v1][d] += div;
	vertexDerivs[v2][d] -= div;
      }
      if ( distance<.1*restinglength[i]){
	double force=vertexDerivs[v1][d];
	vertexDerivs[v1][d] +=vertexDerivs[v2][d];
	vertexDerivs[v2][d] +=force;
      }
      
    }
  }  
}

void VertexFromExternalSpring::update(Tissue &T,
				      DataMatrix &cellData,
				      DataMatrix &wallData,
				      DataMatrix &vertexData, 
				      double h) 
{
  double Lmaxfactor=parameter(2);
  double Kgrowth=parameter(3);
  
  //size_t Npairs=indValue[1].size();
  size_t dimension=vertexData[0].size();
  
  for (size_t i=0 ; i< Npairs ; i++){
    size_t v1=variableIndex(1,i);
    size_t v2=variableIndex(2,i);
    double distance=0;
    for( size_t d=0 ; d<dimension ; d++ ) 
      distance += (vertexData[v2][d]-vertexData[v1][d])*(vertexData[v2][d]-vertexData[v1][d]);
    
    distance=std::sqrt(distance);
    
    if (distance>Lmaxfactor*restinglength[i])  Kspring[i]=0;  
    if (distance<Lmaxfactor*restinglength[i])  Kspring[i]=parameter(0);  
  
    if (Kspring[i]!=0 && restinglength[i]>0){ 
      if(variableIndex(0,0)==1) restinglength[i]+=h*Kgrowth ;
      if(variableIndex(0,0)==2) restinglength[i]+=h*Kgrowth*restinglength[i] ;
      if(variableIndex(0,0)==3) restinglength[i]+=h*Kgrowth*(distance-restinglength[i]) ;
      if(variableIndex(0,0)==4 && distance>restinglength[i]) restinglength[i]+=h*Kgrowth ;
      if(variableIndex(0,0)==5 && distance>restinglength[i]) restinglength[i]+=h*Kgrowth*restinglength[i] ;
      if(variableIndex(0,0)==6 && distance>restinglength[i]) restinglength[i]+=h*Kgrowth*(distance-restinglength[i]) ;
      
      
    }

    
  }
  //std::cerr<<"resting   "<<restinglength[14]<<std::endl; 
}


VertexFromExternalSpringFromPerpVertex::
VertexFromExternalSpringFromPerpVertex(std::vector<double> &paraValue, 
			 std::vector< std::vector<size_t> > &indValue ) 
{  
  // Do some checks on the parameters and variable indeces
  if( paraValue.size()!=8 ) {
    std::cerr << "VertexFromExternalSpringFromPerpVertex::"
	      << "VertexFromExternalSpringFromPerpVertex() "
	      << "Uses eight parameters spring constant K, frac_adhesion,"
	      << "Lmaxfactor, growth_rate, intraction-angle, corner_angle, "
	      << "growthDecay, growth_rate_stress. "<< std::endl;

    exit(EXIT_FAILURE);
  }
  if( indValue.size() != 1 || indValue[0].size()!=4 ) {
    std::cerr << " VertexFromExternalSpringFromPerpVertex::"
	      << " VertexFromExternalSpringFromPerpVertex()"
	      << " one level four indices for gowth flag" 
	      << " connection_flag (constraint on first vertex(1)"
	      << " or both vertices(2)), exclude_corner_flag,"
	      << " and initiate flag (0: for all vertices, 1: only"
	      << " the vertices in the sisterVertex list)" 
	      << " constraint on connections." << std::endl;
    exit(EXIT_FAILURE);
  }
  //Set the variable values
  setId("VertexFromExternalSpringFromPerpVertex");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "K_sp";
  tmp[1] = "frac_adh";
  tmp[2] = "Lmaxfactor";
  tmp[3] = "growth_rate";
  tmp[4] = "intraction_angle";
  tmp[5] = "corner_angle";
  tmp[6] = "growth_rate_decay_rate";
  tmp[7] = "growth_rate_stress";

  setParameterId( tmp ); 
}


void VertexFromExternalSpringFromPerpVertex::
initiate(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs) 
{   
  size_t numCells = T.numCell();
  connections.resize(numCells);
  for (size_t cellIndex = 0; cellIndex <numCells; cellIndex++){
    size_t numCellVertices= T.cell(cellIndex).numVertex();
    connections[cellIndex].resize(numCellVertices);
    for (size_t vertexCellIndex = 0; vertexCellIndex < numCellVertices; vertexCellIndex++)
      {
	//size_t vertex = T.cell(cellIndex).vertex(verIndex)->index();
	connections[cellIndex][vertexCellIndex].resize(numCellVertices,0);
      }
  }
  
  std:: vector<int> hasSister;
  hasSister.resize(T.numVertex());
  size_t NsisterPairs = T.numSisterVertex();
  for(size_t is=0 ; is<NsisterPairs ; is++) {
    hasSister[T.sisterVertex(is,0)]=1; 
    hasSister[T.sisterVertex(is,1)]=1; 
  }
  
  vertexVec.resize(T.numVertex());
  
  for(size_t i=0 ; i<T.numVertex() ; i++) 
    vertexVec[i].resize(3,0);

  // calculating normals to the membrane
  for (size_t cellIndex=0 ; cellIndex< numCells ; cellIndex++){
    size_t numCellVertices= T.cell(cellIndex).numVertex();
    for (size_t vertexCellIndex=0 ; vertexCellIndex< numCellVertices ; vertexCellIndex++){
      
      size_t vertexIndex = T.cell(cellIndex).vertex(vertexCellIndex)->index();
      size_t vertexPlusIndex;
      size_t vertexMinusIndex;
      if (vertexCellIndex!=0 )
	vertexMinusIndex = T.cell(cellIndex).vertex(vertexCellIndex-1)->index();
      else
	vertexMinusIndex= T.cell(cellIndex).vertex(numCellVertices-1)->index();
      if ( vertexCellIndex!=numCellVertices-1 )
	vertexPlusIndex = T.cell(cellIndex).vertex(vertexCellIndex+1)->index();
      else
	vertexPlusIndex = T.cell(cellIndex).vertex(0)->index();
      
      DataMatrix position(3,vertexData[vertexMinusIndex]);
      position[1] = vertexData[vertexIndex];
      position[2] = vertexData[vertexPlusIndex];
      // position[0][2] z for vertexMinus
      double right[3]={position[2][0]-position[1][0] ,
		       position[2][1]-position[1][1] ,
		       position[2][2]-position[1][2] };
      double left[3]={position[0][0]-position[1][0] ,
		      position[0][1]-position[1][1] ,
		      position[0][2]-position[1][2] };
      double tmp=std::sqrt(right[0]*right[0]+right[1]*right[1]+right[2]*right[2]);
      if (tmp!=0){
	right[0]/=tmp;
	right[1]/=tmp;
	right[2]/=tmp;
      }
      tmp=std::sqrt(left[0]*left[0]+left[1]*left[1]+left[2]*left[2]);
      if (tmp!=0){
	left[0]/=tmp;
	left[1]/=tmp;
	left[2]/=tmp;
      }
      vertexVec[vertexIndex][0]=-(right[1]-left[1]);
      vertexVec[vertexIndex][1]=right[0]-left[0];
      tmp=std::sqrt(
		    vertexVec[vertexIndex][0]*vertexVec[vertexIndex][0]+
		    vertexVec[vertexIndex][1]*vertexVec[vertexIndex][1]);
      if (tmp!=0){ 
	vertexVec[vertexIndex][0]/=tmp;
	vertexVec[vertexIndex][1]/=tmp;
      }
      else{ // if two consecutive edges overlay the right one is sellected!!
	vertexVec[vertexIndex][0]=right[0];
	vertexVec[vertexIndex][1]=right[1];
      }
    }
  }
  
  //setting internal actins
  for (size_t cellIndex=0 ; cellIndex< numCells ; cellIndex++){
    size_t numCellVertices= T.cell(cellIndex).numVertex();
    for (size_t vertexCellIndex1=0 ; vertexCellIndex1< numCellVertices ; vertexCellIndex1++){
      
      size_t vertexIndex1 = T.cell(cellIndex).vertex(vertexCellIndex1)->index();
      DataMatrix position(2,vertexData[vertexIndex1]);
      
      for (size_t vertexCellIndex2=0 ; vertexCellIndex2< numCellVertices ; vertexCellIndex2++)
  	if (vertexCellIndex2!=vertexCellIndex1){
  	  size_t vertexIndex2= T.cell(cellIndex).vertex(vertexCellIndex2)->index();
	  position[1] = vertexData[vertexIndex2];
	  
  	  // double  N1N2=
  	  //   vertexVec[vertexIndex1][0]*vertexVec[vertexIndex2][0]+
  	  //   vertexVec[vertexIndex1][1]*vertexVec[vertexIndex2][1];
          
	  double v1v2[2]={vertexData[vertexIndex2][0]-vertexData[vertexIndex1][0],
			  vertexData[vertexIndex2][1]-vertexData[vertexIndex1][1]};
	  
	  double tmp=std::sqrt(v1v2[0]*v1v2[0]+v1v2[1]*v1v2[1]);
	  if (tmp!=0){
	    v1v2[0]/=tmp;
	    v1v2[1]/=tmp;
	  }
	  double teta1=std::acos(
				 vertexVec[vertexIndex1][0]*v1v2[0]+
				 vertexVec[vertexIndex1][1]*v1v2[1]
				 );
	  double teta2=std::acos(
				 -vertexVec[vertexIndex2][0]*v1v2[0]+
				 -vertexVec[vertexIndex2][1]*v1v2[1]
				 );
	  
	  //if (tmp==0) 
	  //  std::cerr<<"cell "<<cellIndex<<" N  " << vertexIndex1 <<" N " << vertexIndex2 <<" "<< N1N2<<" teta "<< teta<<std::endl;
	  if(variableIndex(0,3)==0 || (variableIndex(0,3)==1 && hasSister[vertexIndex1]==1 )){
	    if(variableIndex(0,2)==0){       // no constrain on corners
	      if (variableIndex(0,1)==1 && teta1<(parameter(4)*3.1416/180) ) 
		connections[cellIndex][vertexCellIndex1][vertexCellIndex2]
		  = std::sqrt((position[0][0]-position[1][0])*(position[0][0]-position[1][0])+
			      (position[0][1]-position[1][1])*(position[0][1]-position[1][1])); 
	      
	      if (variableIndex(0,1)==2 && teta1<(parameter(4)*3.1416/180) && teta2<(parameter(4)*3.1416/180)) 
		connections[cellIndex][vertexCellIndex1][vertexCellIndex2]
		  = std::sqrt((position[0][0]-position[1][0])*(position[0][0]-position[1][0])+
			      (position[0][1]-position[1][1])*(position[0][1]-position[1][1])); 
	    }
	    
	    if(variableIndex(0,2)==1) {      // exclude_corner
	      
	      // size_t vertexIndex1 = T.cell(cellIndex).vertex(vertex1)->index();
	      size_t vertexPlusIndex;
	      size_t vertexMinusIndex;
	      
	      if (vertexCellIndex1!=0 )
		vertexMinusIndex = T.cell(cellIndex).vertex(vertexCellIndex1-1)->index();
	      else
		vertexMinusIndex= T.cell(cellIndex).vertex(numCellVertices-1)->index();
	      if ( vertexCellIndex1!=numCellVertices-1 )
		vertexPlusIndex = T.cell(cellIndex).vertex(vertexCellIndex1+1)->index();
	      else
		vertexPlusIndex = T.cell(cellIndex).vertex(0)->index();
	      
	      double cornerAngle=vertexVec[vertexMinusIndex][0]*vertexVec[vertexPlusIndex][0]+
		vertexVec[vertexMinusIndex][1]*vertexVec[vertexPlusIndex][1];
	      if (cornerAngle > std::cos(3.1416*(180-parameter(5))/180)){
		if (variableIndex(0,1)==1 && teta1<(parameter(4)*3.1416/180) ) 
		  connections[cellIndex][vertexCellIndex1][vertexCellIndex2]=
		    std::sqrt((position[0][0]-position[1][0])*(position[0][0]-position[1][0])+
			      (position[0][1]-position[1][1])*(position[0][1]-position[1][1])); 
		
		if (variableIndex(0,1)==2 && teta1<(parameter(4)*3.1416/180) && teta2<(parameter(4)*3.1416/180)) 
		  connections[cellIndex][vertexCellIndex1][vertexCellIndex2]= 
		    std::sqrt((position[0][0]-position[1][0])*(position[0][0]-position[1][0])+
			      (position[0][1]-position[1][1])*(position[0][1]-position[1][1])); 
	      }
	    }
	  }
	  connections[cellIndex][vertexCellIndex2][vertexCellIndex1]=
	    connections[cellIndex][vertexCellIndex1][vertexCellIndex2];
	  
  	}
      
    }
  }
  
  
  // for(size_t zx=0; zx<numCells;zx++){
  //   std::cerr<<"cell  "<<zx<<std::endl;
  //   size_t numV= T.cell(zx).numVertex();
  //   for (int a=1;a<numV;a++)
  //     for (int b=1;b<numV;b++)
  // 	if (connections[zx][a][b]!=0)
  // 	  {
  // 	  std::cerr<<"a,b  "<<a<<" , "<<b<<" connections["<<zx<<"][a][b] is "<<connections[zx][a][b]<<std::endl;
	 
  // 	  }
  // }
  
  //size_t Npairs=variableIndex(1).size();    
  // Print the edges gnuplot style to check connectivity
  size_t printFlag=0;
  if (printFlag) {
    for (size_t i=0; i<numCells; ++i) {
      size_t numCellVertices= T.cell(i).numVertex();
      for (size_t j=0; j<numCellVertices; ++j) {
	for (size_t k=j+1; k<numCellVertices; ++k) {
	  if (connections[i][j][k]) {
	    std::cerr << "In cell " << i << " vertices " << j << " and " << k << " connected." << std::endl;  
	    size_t numDimension = vertexData[0].size();
	    size_t v1 = T.cell(i).vertex(j)->index();
	    size_t v2 = T.cell(i).vertex(k)->index();
	    for (size_t d=0; d<numDimension; ++d) {
	      std::cout << vertexData[v1][d] << " ";
	    }
	    std::cout << i << std::endl;
	    for (size_t d=0; d<numDimension; ++d) {
	      std::cout << vertexData[v2][d] << " ";
	    }
	    std::cout << i << std::endl;
	    std::cout << std::endl;
	  }
	}
      }
    }
    std::cerr << "End of printing." << std::endl;
    // END printing
  }
}


void VertexFromExternalSpringFromPerpVertex::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{  
  double fad=parameter(1);
 
  //size_t Npairs=indValue[1].size();
  size_t dimension=vertexData[0].size();

  size_t numCells = T.numCell();

  for (size_t cellIndex=0 ; cellIndex< numCells ; cellIndex++){
    size_t numVertices= T.cell(cellIndex).numVertex();
    for (size_t verIndex1=0 ; verIndex1< numVertices-1 ; verIndex1++){
      
      size_t vertex1 = T.cell(cellIndex).vertex(verIndex1)->index();
     
      
      for (size_t verIndex2=verIndex1+1 ; verIndex2< numVertices ; verIndex2++)
  	if (verIndex2!=verIndex1 && connections[cellIndex][verIndex1][verIndex2]!=0){
  	  size_t vertex2= T.cell(cellIndex).vertex(verIndex2)->index();
  	 
	  
	  double distance=0;
	  
	  for( size_t d=0 ; d<dimension ; d++ ) 
	    distance += (vertexData[vertex2][d]-vertexData[vertex1][d])*(vertexData[vertex2][d]-vertexData[vertex1][d]);
	  
	  distance=std::sqrt(distance);
  	double coeff=0;
	if ( connections[cellIndex][verIndex1][verIndex2]>0.0)
	  coeff=parameter(0)*((1.0/connections[cellIndex][verIndex1][verIndex2])-(1.0/distance));
	
	if( distance <= 0.0 && connections[cellIndex][verIndex1][verIndex2] <=0.0 ) 
	  coeff = 0.0;
	
	if( distance>connections[cellIndex][verIndex1][verIndex2])
	  coeff *=fad;
	//Update both vertices for each dimension
	for(size_t d=0 ; d<dimension ; d++ ) {
	  double div = coeff*(vertexData[vertex2][d]-vertexData[vertex1][d]);
	  if ( connections[cellIndex][verIndex1][verIndex2]>0){
	    vertexDerivs[vertex1][d] += div;
	    vertexDerivs[vertex2][d] -= div;
	  }
	  // if ( distance<.1*restinglength[i]){
	  //   double force=vertexDerivs[vertex1][d];
	  //   vertexDerivs[vertex1][d] +=vertexDerivs[v2][d];
	  //   vertexDerivs[vertex2][d] +=force;
	  // }  
	}
	}
      
    }
    
    
  }  
}

void VertexFromExternalSpringFromPerpVertex::update(Tissue &T,
				      DataMatrix &cellData,
				      DataMatrix &wallData,
				      DataMatrix &vertexData, 
				      double h) 
{
  // double Lmaxfactor=parameter(2);
  double Kgrowth=parameter(3);
  double KgrowthStress=parameter(7);

  //size_t Npairs=indValue[1].size();
  size_t dimension=vertexData[0].size();

  std:: vector<int> hasSister; 
  hasSister.resize(T.numVertex());
  size_t NsisterPairs = T.numSisterVertex();
  for(size_t is=0 ; is<NsisterPairs ; is++) {
    hasSister[T.sisterVertex(is,0)]=1; 
    hasSister[T.sisterVertex(is,1)]=1; 
  }
  
  size_t numCells = T.numCell();
  
  for (size_t cellIndex=0 ; cellIndex< numCells ; cellIndex++){
    size_t numVertices= T.cell(cellIndex).numVertex();
    for (size_t verIndex1=0 ; verIndex1< numVertices-1 ; verIndex1++){
      size_t vertex1 = T.cell(cellIndex).vertex(verIndex1)->index();
      for (size_t verIndex2=verIndex1+1 ; verIndex2< numVertices ; verIndex2++)
	if (verIndex2!=verIndex1 && connections[cellIndex][verIndex1][verIndex2]!=0){
	  double restinglength=connections[cellIndex][verIndex1][verIndex2];
	  size_t vertex2= T.cell(cellIndex).vertex(verIndex2)->index();
	  
	  double distance=0;
	  for( size_t d=0 ; d<dimension ; d++ ) 
	    distance += (vertexData[vertex2][d]-vertexData[vertex1][d])*
	      (vertexData[vertex2][d]-vertexData[vertex1][d]);
	  
	  distance=std::sqrt(distance);
	  
	  if (restinglength>0){ 
	    if(variableIndex(0,0)==1) 
	      restinglength+=h*Kgrowth;
	    else if(variableIndex(0,0)==2) 
	      restinglength+=h*Kgrowth*restinglength ;
	    else if(variableIndex(0,0)==3) 
	      restinglength+=h*Kgrowth*(distance-restinglength) ;
	    else if(variableIndex(0,0)==4 && distance>restinglength) 
	      restinglength+=h*Kgrowth ;
	    else if(variableIndex(0,0)==5 && distance>restinglength) 
	      restinglength+=h*Kgrowth*restinglength ;
	    else if(variableIndex(0,0)==6 && distance>restinglength) 
	      restinglength+=h*Kgrowth*(distance-restinglength) ;	  
	    else if(variableIndex(0,0)==7) {
	      restinglength += h*KgrowthStress*(distance-restinglength)
		+ h*Kgrowth*restinglength;
	    }
	    else if(variableIndex(0,0)==8) {
	      if(hasSister[vertex1]==1 || hasSister[vertex2]==1)
		restinglength+=h*Kgrowth*restinglength ;
	      else{
		restinglength+=h*KgrowthStress*restinglength;
		//std::cerr<<"broken sisters"<<std::endl;
		
	      }
	    }
	    
	    else {
	      std::cerr << "VertexFromExternalSpringFromPerpVertex::update()"
			<< std::endl << "Wrong growth rule index given!"
			<< std::endl;
	      std::exit(EXIT_FAILURE);
	    }
	    connections[cellIndex][verIndex1][verIndex2]=connections[cellIndex][verIndex2][verIndex1]=restinglength;
	  }					
	}
    }
  }

 
  setParameter(3,parameter(3)*parameter(6));
  //std::cerr << parameter(3) << std::endl;
}

void VertexFromExternalSpringFromPerpVertex::
printState(Tissue *T,
	   DataMatrix &cellData,
	   DataMatrix &wallData,
	   DataMatrix &vertexData, 
	   std::ostream &os)
{
  static size_t index=0;
  size_t counter=0; 
  size_t dimension = vertexData[0].size();
  size_t numCells = T->numCell();
  for (size_t cellIndex=0 ; cellIndex< numCells ; cellIndex++){
    size_t numVertices= T->cell(cellIndex).numVertex();
    for (size_t verIndex1=0 ; verIndex1< numVertices-1 ; verIndex1++){
      size_t vertex1 = T->cell(cellIndex).vertex(verIndex1)->index();
      for (size_t verIndex2=verIndex1+1 ; verIndex2< numVertices ; verIndex2++)
	if (verIndex2!=verIndex1 && connections[cellIndex][verIndex1][verIndex2]!=0){
	  double restingLength=connections[cellIndex][verIndex1][verIndex2];
	  size_t vertex2= T->cell(cellIndex).vertex(verIndex2)->index();
	  double distance=0;
	  for( size_t d=0 ; d<dimension ; d++ )
	    distance += (vertexData[vertex2][d]-vertexData[vertex1][d])*
	      (vertexData[vertex2][d]-vertexData[vertex1][d]);	  
	  distance=std::sqrt(distance);

	  //vertex1
	  os << index << " " << counter << " ";
	  for( size_t d=0 ; d<dimension ; d++ )
	    os << vertexData[vertex1][d] << " ";
	  os << restingLength << " " << distance << " " 
	     << parameter(3)*restingLength+parameter(7)*(distance-restingLength) << std::endl;
	  //vertex2
	  os << index << " " << counter << " ";
	  for( size_t d=0 ; d<dimension ; d++ )
	    os << vertexData[vertex2][d] << " ";
	  os << restingLength << " " << distance
	     << parameter(3)*restingLength+parameter(7)*(distance-restingLength) << std::endl;
	  os << std::endl;
	  counter++;
	}
    }
  }
  index++;
  //std::cerr<<" index "<<index<<"  counter "<< counter<< std::endl;
}


////////////////////////////////////////////////////////////////

VertexFromExternalSpringFromPerpVertexDynamic::
VertexFromExternalSpringFromPerpVertexDynamic(std::vector<double> &paraValue, 
			 std::vector< std::vector<size_t> > &indValue ) 
{  
  // Do some checks on the parameters and variable indeces
  if( paraValue.size()!=8 ) {
    std::cerr << "VertexFromExternalSpringFromPerpVertexDynamic::"
	      << "VertexFromExternalSpringFromPerpVertexDynamic() "
	      << "Uses eight parameters spring constant K, frac_adhesion,"
	      << "Lmaxfactor, growth_rate, intraction-angle, corner_angle, "
	      << "growthDecay, growth_rate_stress. "<< std::endl;

    exit(EXIT_FAILURE);
  }
  if( indValue.size() != 1 || indValue[0].size()!=4 ) {
    std::cerr << " VertexFromExternalSpringFromPerpVertexDynamic::"
	      << " VertexFromExternalSpringFromPerpVertexDynamic()"
	      << " one level four indices for gowth flag" 
	      << " connection_flag (constraint on first vertex(1)"
	      << " or both vertices(2)), exclude_corner_flag,"
	      << " and initiate flag (0: for all vertices, 1: only"
	      << " the vertices in the sisterVertex list)" 
	      << " constraint on connections." << std::endl;
    exit(EXIT_FAILURE);
  }
  //Set the variable values
  setId("VertexFromExternalSpringFromPerpVertexDynamic");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "K_sp";
  tmp[1] = "frac_adh";
  tmp[2] = "Lmaxfactor";
  tmp[3] = "growth_rate";
  tmp[4] = "intraction_angle";
  tmp[5] = "corner_angle";
  tmp[6] = "growth_rate_decay_rate";
  tmp[7] = "growth_rate_stress";

  setParameterId( tmp ); 
}

void VertexFromExternalSpringFromPerpVertexDynamic::
initiate(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs) 
{   
  size_t numCells = T.numCell();
  connections.resize(numCells);
  for (size_t cellIndex = 0; cellIndex <numCells; cellIndex++){
    size_t numCellVertices= T.cell(cellIndex).numVertex();
    connections[cellIndex].resize(numCellVertices);
    for (size_t vertexCellIndex = 0; vertexCellIndex < numCellVertices; vertexCellIndex++)
      {
	//size_t vertex = T.cell(cellIndex).vertex(verIndex)->index();
	connections[cellIndex][vertexCellIndex].resize(numCellVertices,0);
      }
  }
  
  std:: vector<int> hasSister;
  hasSister.resize(T.numVertex());
  size_t NsisterPairs = T.numSisterVertex();
  for(size_t is=0 ; is<NsisterPairs ; is++) {
    hasSister[T.sisterVertex(is,0)]=1; 
    hasSister[T.sisterVertex(is,1)]=1; 
  }
  
  vertexVec.resize(T.numVertex());
  
  for(size_t i=0 ; i<T.numVertex() ; i++) 
    vertexVec[i].resize(3,0);

  // calculating normals to the membrane
  for (size_t cellIndex=0 ; cellIndex< numCells ; cellIndex++){
    size_t numCellVertices= T.cell(cellIndex).numVertex();
    for (size_t vertexCellIndex=0 ; vertexCellIndex< numCellVertices ; vertexCellIndex++){
      
      size_t vertexIndex = T.cell(cellIndex).vertex(vertexCellIndex)->index();
      size_t vertexPlusIndex;
      size_t vertexMinusIndex;
      if (vertexCellIndex!=0 )
	vertexMinusIndex = T.cell(cellIndex).vertex(vertexCellIndex-1)->index();
      else
	vertexMinusIndex= T.cell(cellIndex).vertex(numCellVertices-1)->index();
      if ( vertexCellIndex!=numCellVertices-1 )
	vertexPlusIndex = T.cell(cellIndex).vertex(vertexCellIndex+1)->index();
      else
	vertexPlusIndex = T.cell(cellIndex).vertex(0)->index();
      
      DataMatrix position(3,vertexData[vertexMinusIndex]);
      position[1] = vertexData[vertexIndex];
      position[2] = vertexData[vertexPlusIndex];
      // position[0][2] z for vertexMinus
      double right[3]={position[2][0]-position[1][0] ,
		       position[2][1]-position[1][1] ,
		       position[2][2]-position[1][2] };
      double left[3]={position[0][0]-position[1][0] ,
		      position[0][1]-position[1][1] ,
		      position[0][2]-position[1][2] };
      double tmp=std::sqrt(right[0]*right[0]+right[1]*right[1]+right[2]*right[2]);
      if (tmp!=0){
	right[0]/=tmp;
	right[1]/=tmp;
	right[2]/=tmp;
      }
      tmp=std::sqrt(left[0]*left[0]+left[1]*left[1]+left[2]*left[2]);
      if (tmp!=0){
	left[0]/=tmp;
	left[1]/=tmp;
	left[2]/=tmp;
      }
      vertexVec[vertexIndex][0]=-(right[1]-left[1]);
      vertexVec[vertexIndex][1]=right[0]-left[0];
      tmp=std::sqrt(
		    vertexVec[vertexIndex][0]*vertexVec[vertexIndex][0]+
		    vertexVec[vertexIndex][1]*vertexVec[vertexIndex][1]);
      if (tmp!=0){ 
	vertexVec[vertexIndex][0]/=tmp;
	vertexVec[vertexIndex][1]/=tmp;
      }
      else{ // if two consecutive edges overlay the right one is sellected!!
	vertexVec[vertexIndex][0]=right[0];
	vertexVec[vertexIndex][1]=right[1];
      }
    }
  }
  
  //setting internal actins
  for (size_t cellIndex=0 ; cellIndex< numCells ; cellIndex++){
    size_t numCellVertices= T.cell(cellIndex).numVertex();
    for (size_t vertexCellIndex1=0 ; vertexCellIndex1< numCellVertices ; vertexCellIndex1++){
      
      size_t vertexIndex1 = T.cell(cellIndex).vertex(vertexCellIndex1)->index();
      DataMatrix position(2,vertexData[vertexIndex1]);
      
      for (size_t vertexCellIndex2=0 ; vertexCellIndex2< numCellVertices ; vertexCellIndex2++)
  	if (vertexCellIndex2!=vertexCellIndex1){
  	  size_t vertexIndex2= T.cell(cellIndex).vertex(vertexCellIndex2)->index();
	  position[1] = vertexData[vertexIndex2];
	  
  	  // double  N1N2=
  	  //   vertexVec[vertexIndex1][0]*vertexVec[vertexIndex2][0]+
  	  //   vertexVec[vertexIndex1][1]*vertexVec[vertexIndex2][1];
          
	  double v1v2[2]={vertexData[vertexIndex2][0]-vertexData[vertexIndex1][0],
			  vertexData[vertexIndex2][1]-vertexData[vertexIndex1][1]};
	  
	  double tmp=std::sqrt(v1v2[0]*v1v2[0]+v1v2[1]*v1v2[1]);
	  if (tmp!=0){
	    v1v2[0]/=tmp;
	    v1v2[1]/=tmp;
	  }
	  double teta1=std::acos(
				 vertexVec[vertexIndex1][0]*v1v2[0]+
				 vertexVec[vertexIndex1][1]*v1v2[1]
				 );
	  double teta2=std::acos(
				 -vertexVec[vertexIndex2][0]*v1v2[0]+
				 -vertexVec[vertexIndex2][1]*v1v2[1]
				 );
	  
	  //if (tmp==0) 
	  //  std::cerr<<"cell "<<cellIndex<<" N  " << vertexIndex1 <<" N " << vertexIndex2 <<" "<< N1N2<<" teta "<< teta<<std::endl;
	  if(variableIndex(0,3)==0 || (variableIndex(0,3)==1 && hasSister[vertexIndex1]==1 )){
	    if(variableIndex(0,2)==0){       // no constrain on corners
	      if (variableIndex(0,1)==1 && teta1<(parameter(4)*3.1416/180) ) 
		connections[cellIndex][vertexCellIndex1][vertexCellIndex2]
		  = std::sqrt((position[0][0]-position[1][0])*(position[0][0]-position[1][0])+
			      (position[0][1]-position[1][1])*(position[0][1]-position[1][1])); 
	      
	      if (variableIndex(0,1)==2 && teta1<(parameter(4)*3.1416/180) && teta2<(parameter(4)*3.1416/180)) 
		connections[cellIndex][vertexCellIndex1][vertexCellIndex2]
		  = std::sqrt((position[0][0]-position[1][0])*(position[0][0]-position[1][0])+
			      (position[0][1]-position[1][1])*(position[0][1]-position[1][1])); 
	    }
	    
	    if(variableIndex(0,2)==1) {      // exclude_corner
	      
	      // size_t vertexIndex1 = T.cell(cellIndex).vertex(vertex1)->index();
	      size_t vertexPlusIndex;
	      size_t vertexMinusIndex;
	      
	      if (vertexCellIndex1!=0 )
		vertexMinusIndex = T.cell(cellIndex).vertex(vertexCellIndex1-1)->index();
	      else
		vertexMinusIndex= T.cell(cellIndex).vertex(numCellVertices-1)->index();
	      if ( vertexCellIndex1!=numCellVertices-1 )
		vertexPlusIndex = T.cell(cellIndex).vertex(vertexCellIndex1+1)->index();
	      else
		vertexPlusIndex = T.cell(cellIndex).vertex(0)->index();
	      
	      double cornerAngle=vertexVec[vertexMinusIndex][0]*vertexVec[vertexPlusIndex][0]+
		vertexVec[vertexMinusIndex][1]*vertexVec[vertexPlusIndex][1];
	      if (cornerAngle > std::cos(3.1416*(180-parameter(5))/180)){
		if (variableIndex(0,1)==1 && teta1<(parameter(4)*3.1416/180) ) 
		  connections[cellIndex][vertexCellIndex1][vertexCellIndex2]=
		    std::sqrt((position[0][0]-position[1][0])*(position[0][0]-position[1][0])+
			      (position[0][1]-position[1][1])*(position[0][1]-position[1][1])); 
		
		if (variableIndex(0,1)==2 && teta1<(parameter(4)*3.1416/180) && teta2<(parameter(4)*3.1416/180)) 
		  connections[cellIndex][vertexCellIndex1][vertexCellIndex2]= 
		    std::sqrt((position[0][0]-position[1][0])*(position[0][0]-position[1][0])+
			      (position[0][1]-position[1][1])*(position[0][1]-position[1][1])); 
	      }
	    }
	  }
	  connections[cellIndex][vertexCellIndex2][vertexCellIndex1]=
	    connections[cellIndex][vertexCellIndex1][vertexCellIndex2];
	  
  	}
      
    }
  }
  
  
  // for(size_t zx=0; zx<numCells;zx++){
  //   std::cerr<<"cell  "<<zx<<std::endl;
  //   size_t numV= T.cell(zx).numVertex();
  //   for (int a=1;a<numV;a++)
  //     for (int b=1;b<numV;b++)
  // 	if (connections[zx][a][b]!=0)
  // 	  {
  // 	  std::cerr<<"a,b  "<<a<<" , "<<b<<" connections["<<zx<<"][a][b] is "<<connections[zx][a][b]<<std::endl;
	 
  // 	  }
  // }
  
  //size_t Npairs=variableIndex(1).size();    
  // Print the edges gnuplot style to check connectivity
  size_t printFlag=0;
  if (printFlag) {
    for (size_t i=0; i<numCells; ++i) {
      size_t numCellVertices= T.cell(i).numVertex();
      for (size_t j=0; j<numCellVertices; ++j) {
	for (size_t k=j+1; k<numCellVertices; ++k) {
	  if (connections[i][j][k]) {
	    std::cerr << "In cell " << i << " vertices " << j << " and " << k << " connected." << std::endl;  
	    size_t numDimension = vertexData[0].size();
	    size_t v1 = T.cell(i).vertex(j)->index();
	    size_t v2 = T.cell(i).vertex(k)->index();
	    for (size_t d=0; d<numDimension; ++d) {
	      std::cout << vertexData[v1][d] << " ";
	    }
	    std::cout << i << std::endl;
	    for (size_t d=0; d<numDimension; ++d) {
	      std::cout << vertexData[v2][d] << " ";
	    }
	    std::cout << i << std::endl;
	    std::cout << std::endl;
	  }
	}
      }
    }
    std::cerr << "End of printing." << std::endl;
    // END printing
  }
}


void VertexFromExternalSpringFromPerpVertexDynamic::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{  
  double fad=parameter(1);
  //size_t Npairs=indValue[1].size();
  size_t dimension=vertexData[0].size();
  
  size_t numCells = T.numCell();
  for (size_t cellIndex=0 ; cellIndex< numCells ; cellIndex++){
    size_t numCellVertices= T.cell(cellIndex).numVertex();
    for (size_t vertexCellIndex1=0 ; vertexCellIndex1< numCellVertices-1 ; vertexCellIndex1++){
      
      size_t vertexIndex1 = T.cell(cellIndex).vertex(vertexCellIndex1)->index();
      
      for (size_t vertexCellIndex2=vertexCellIndex1+1 ; vertexCellIndex2< numCellVertices ; vertexCellIndex2++)
  	if (vertexCellIndex2!=vertexCellIndex1 && connections[cellIndex][vertexCellIndex1][vertexCellIndex2]!=0){
  	  size_t vertexIndex2= T.cell(cellIndex).vertex(vertexCellIndex2)->index();
	  
	  double distance=0;
	  
	  for( size_t d=0 ; d<dimension ; d++ ) 
	    distance += (vertexData[vertexIndex2][d]-vertexData[vertexIndex1][d])
	      *(vertexData[vertexIndex2][d]-vertexData[vertexIndex1][d]);
	  
	  distance=std::sqrt(distance);
	  double coeff=0;
	  if ( connections[cellIndex][vertexCellIndex1][vertexCellIndex2]>0.0)
	    coeff=parameter(0)*((1.0/connections[cellIndex][vertexCellIndex1][vertexCellIndex2])-(1.0/distance));
	  
	  if( distance <= 0.0 && connections[cellIndex][vertexCellIndex1][vertexCellIndex2] <=0.0 ) 
	    coeff = 0.0;
	  
	  if( distance>connections[cellIndex][vertexCellIndex1][vertexCellIndex2])
	    coeff *=fad;
	  //Update both vertices for each dimension
	  for(size_t d=0 ; d<dimension ; d++ ) {
	    double div = coeff*(vertexData[vertexIndex2][d]-vertexData[vertexIndex1][d]);
	    if ( connections[cellIndex][vertexCellIndex1][vertexCellIndex2]>0){
	      vertexDerivs[vertexIndex1][d] += div;
	      vertexDerivs[vertexIndex2][d] -= div;
	    }
	    // if ( distance<.1*restinglength[i]){
	    //   double force=vertexDerivs[vertex1][d];
	    //   vertexDerivs[vertex1][d] +=vertexDerivs[v2][d];
	    //   vertexDerivs[vertex2][d] +=force;
	    // }  
	  }
	}
    }
  }  
}

void VertexFromExternalSpringFromPerpVertexDynamic::update(Tissue &T,
				      DataMatrix &cellData,
				      DataMatrix &wallData,
				      DataMatrix &vertexData, 
				      double h) 
{ // double Lmaxfactor=parameter(2);
  double Kgrowth=parameter(3);
  double KgrowthStress=parameter(7);
  
  //size_t Npairs=indValue[1].size();
  size_t dimension=vertexData[0].size();
  
  size_t numCells = T.numCell(); 
  //size_t numTotalVertices = T.numVertex();
  
  // setting flag for sistership
  std:: vector<int> hasSister; 
  hasSister.resize(T.numVertex());
  size_t NsisterPairs = T.numSisterVertex();
  for(size_t is=0 ; is<NsisterPairs ; is++) {
    hasSister[T.sisterVertex(is,0)]=1; 
    hasSister[T.sisterVertex(is,1)]=1; 
  }
  
  // internal actin and membrane normal update
  // updating normals to the membrane Counterclockwise sorting of vertices
  for (size_t cellIndex=0 ; cellIndex< numCells ; cellIndex++){
    size_t numCellVertices= T.cell(cellIndex).numVertex();
    for (size_t vertexCellIndex=0 ; vertexCellIndex< numCellVertices ; vertexCellIndex++){
      
      size_t vertexIndex = T.cell(cellIndex).vertex(vertexCellIndex)->index();
      size_t vertexPlusIndex;
      size_t vertexMinusIndex;
      if (vertexCellIndex!=0 )
  	vertexMinusIndex = T.cell(cellIndex).vertex(vertexCellIndex-1)->index();
      else
  	vertexMinusIndex= T.cell(cellIndex).vertex(numCellVertices-1)->index();
      if ( vertexCellIndex!=numCellVertices-1 )
  	vertexPlusIndex = T.cell(cellIndex).vertex(vertexCellIndex+1)->index();
      else
  	vertexPlusIndex = T.cell(cellIndex).vertex(0)->index();
      DataMatrix position(3,vertexData[vertexMinusIndex]);
      position[1] = vertexData[vertexIndex];
      position[2] = vertexData[vertexPlusIndex];
      // position[0][2] z for vertexMinus
      double right[3]={position[2][0]-position[1][0] ,
  		       position[2][1]-position[1][1] ,
  		       position[2][2]-position[1][2] };
      double left[3]={position[0][0]-position[1][0] ,
  		      position[0][1]-position[1][1] ,
  		      position[0][2]-position[1][2] };
      double tmp=std::sqrt(right[0]*right[0]+right[1]*right[1]+right[2]*right[2]);
      if (tmp!=0){
  	right[0]/=tmp;
  	right[1]/=tmp;
  	right[2]/=tmp;
      }
      tmp=std::sqrt(left[0]*left[0]+left[1]*left[1]+left[2]*left[2]);
      if (tmp!=0){
  	left[0]/=tmp;
  	left[1]/=tmp;
  	left[2]/=tmp;
      }
      vertexVec[vertexIndex][0]=-(right[1]-left[1]);
      vertexVec[vertexIndex][1]=right[0]-left[0];
      tmp=std::sqrt(
  		    vertexVec[vertexIndex][0]*vertexVec[vertexIndex][0]+
  		    vertexVec[vertexIndex][1]*vertexVec[vertexIndex][1]);
      if (tmp!=0){ 
  	vertexVec[vertexIndex][0]/=tmp;
  	vertexVec[vertexIndex][1]/=tmp;
      }
      else{ // if two consecutive edges overlay the right one is sellected!!
  	vertexVec[vertexIndex][0]=right[0];
  	vertexVec[vertexIndex][1]=right[1];
      }
    }
  }
  
  //updating internal actins
  std:: vector<std::vector<std::vector<double> > >  connectionsNew;
  connectionsNew.resize(numCells);
  for (size_t cellIndex = 0; cellIndex <numCells; cellIndex++){
    size_t numCellVertices= T.cell(cellIndex).numVertex();
    connectionsNew[cellIndex].resize(numCellVertices);
    for (size_t vertexCellIndex = 0; vertexCellIndex < numCellVertices; vertexCellIndex++)
      {
	//size_t vertex = T.cell(cellIndex).vertex(verIndex)->index();
	connectionsNew[cellIndex][vertexCellIndex].resize(numCellVertices,0);
      }
  }
  
  for (size_t cellIndex=0 ; cellIndex< numCells ; cellIndex++){
    size_t numCellVertices= T.cell(cellIndex).numVertex();
    for (size_t vertexCellIndex1=0 ; vertexCellIndex1< numCellVertices ; vertexCellIndex1++){
      
      size_t vertexIndex1 = T.cell(cellIndex).vertex(vertexCellIndex1)->index();
      DataMatrix position(2,vertexData[vertexIndex1]);
      
      for (size_t vertexCellIndex2=0 ; vertexCellIndex2< numCellVertices ; vertexCellIndex2++)
  	if (vertexCellIndex2!=vertexCellIndex1){
  	  size_t vertexIndex2= T.cell(cellIndex).vertex(vertexCellIndex2)->index();
  	  position[1] = vertexData[vertexIndex2];
	  
  	  // double  N1N2=
  	  //   vertexVec[vertexIndex1][0]*vertexVec[vertexIndex2][0]+
  	  //   vertexVec[vertexIndex1][1]*vertexVec[vertexIndex2][1];
          
  	  double v1v2[2]={vertexData[vertexIndex2][0]-vertexData[vertexIndex1][0],
  			  vertexData[vertexIndex2][1]-vertexData[vertexIndex1][1]};
	  
  	  double tmp=std::sqrt(v1v2[0]*v1v2[0]+v1v2[1]*v1v2[1]);
  	  if (tmp!=0){
  	    v1v2[0]/=tmp;
  	    v1v2[1]/=tmp;
  	  }
  	  double teta1=std::acos(
  				 vertexVec[vertexIndex1][0]*v1v2[0]+
  				 vertexVec[vertexIndex1][1]*v1v2[1]
  				 );
  	  double teta2=std::acos(
  				 -vertexVec[vertexIndex2][0]*v1v2[0]+
  				 -vertexVec[vertexIndex2][1]*v1v2[1]
  				 );
	  
  	  //if (tmp==0) 
  	  //  std::cerr<<"cell "<<cellIndex<<" N  " << vertexIndex1 <<" N " << vertexIndex2 <<" "<< N1N2<<" teta "<< teta<<std::endl;
	  
  	  if(variableIndex(0,3)==0 || (variableIndex(0,3)==1 && hasSister[vertexIndex1]==1 )){
	    
  	    if(variableIndex(0,2)==0){       // no constrain on corners
	      if (variableIndex(0,1)==1 && teta1<(parameter(4)*3.1416/180) ) 
  		connectionsNew[cellIndex][vertexCellIndex1][vertexCellIndex2]= 
		  std::sqrt( (position[0][0]-position[1][0])*(position[0][0]-position[1][0])+
			     (position[0][1]-position[1][1])*(position[0][1]-position[1][1])); 
	      
  	      if (variableIndex(0,1)==2 && teta1<(parameter(4)*3.1416/180) && teta2<(parameter(4)*3.1416/180)) 
  		connectionsNew[cellIndex][vertexCellIndex1][vertexCellIndex2]= 
	  	 std::sqrt((position[0][0]-position[1][0])*(position[0][0]-position[1][0])+
			   (position[0][1]-position[1][1])*(position[0][1]-position[1][1])); 
	      
  	    }
	    
  	    if(variableIndex(0,2)==1) {      // exclude_corner
	      
  	      // size_t vertexIndex1 = T.cell(cellIndex).vertex(vertex1)->index();
  	      size_t vertexPlusIndex;
  	      size_t vertexMinusIndex;
	      
  	      if (vertexCellIndex1!=0 )
  		vertexMinusIndex = T.cell(cellIndex).vertex(vertexCellIndex1-1)->index();
  	      else
  		vertexMinusIndex= T.cell(cellIndex).vertex(numCellVertices-1)->index();
  	      if ( vertexCellIndex1!=numCellVertices-1 )
  		vertexPlusIndex = T.cell(cellIndex).vertex(vertexCellIndex1+1)->index();
  	      else
  		vertexPlusIndex = T.cell(cellIndex).vertex(0)->index();
	       
  	      double cornerAngle=vertexVec[vertexMinusIndex][0]*vertexVec[vertexPlusIndex][0]+
  		vertexVec[vertexMinusIndex][1]*vertexVec[vertexPlusIndex][1];


  	      if (cornerAngle > std::cos(3.1416*(180-parameter(5))/180)){	 
  		if (variableIndex(0,1)==1 && teta1<(parameter(4)*3.1416/180) ) 
  		  connectionsNew[cellIndex][vertexCellIndex1][vertexCellIndex2]= 
		    std::sqrt((position[0][0]-position[1][0])*(position[0][0]-position[1][0])+
			      (position[0][1]-position[1][1])*(position[0][1]-position[1][1])); 
	
  		if (variableIndex(0,1)==2 && teta1<(parameter(4)*3.1416/180) && teta2<(parameter(4)*3.1416/180)) 
  		  connectionsNew[cellIndex][vertexCellIndex1][vertexCellIndex2]= 
                    std::sqrt((position[0][0]-position[1][0])*(position[0][0]-position[1][0])+
			      (position[0][1]-position[1][1])*(position[0][1]-position[1][1])); 
  	      }
  	    }
  	  }
	  //std::cerr<< "up to here"<<std::endl;
  	  connectionsNew[cellIndex][vertexCellIndex2][vertexCellIndex1]=
	    connectionsNew[cellIndex][vertexCellIndex1][vertexCellIndex2]; 
  	}
    }
  }
  
  // setting new connections if there were no connection already and keeping the old ones
  for(size_t iNew=0; iNew< numCells ;iNew++){
    size_t numCellVertices = T.cell(iNew).numVertex();
    for(size_t jNew=0; jNew< numCellVertices ;jNew++)
      for(size_t kNew=0; kNew< numCellVertices ;kNew++){
  	if(connections[iNew][jNew][kNew]==0)
  	  connections[iNew][jNew][kNew]=connectionsNew[iNew][jNew][kNew];
  	else if(connectionsNew[iNew][jNew][kNew]==0)
  	  connections[iNew][jNew][kNew]=0;
      }
  }
  
  // End of internal actin and membrane normal update
 
  // growth rules  
  for (size_t cellIndex=0 ; cellIndex< numCells ; cellIndex++){
    size_t numCellVertices= T.cell(cellIndex).numVertex();
    for (size_t vertexCellIndex1=0 ; vertexCellIndex1< numCellVertices ; vertexCellIndex1++){
      size_t vertexIndex1 = T.cell(cellIndex).vertex(vertexCellIndex1)->index();
      for (size_t vertexCellIndex2=0 ; vertexCellIndex2< numCellVertices ; vertexCellIndex2++)
	if (vertexCellIndex2!=vertexCellIndex1 && connections[cellIndex][vertexCellIndex1][vertexCellIndex2]!=0){
	  double restinglength=connections[cellIndex][vertexCellIndex1][vertexCellIndex2];
	  size_t vertexIndex2= T.cell(cellIndex).vertex(vertexCellIndex2)->index();
	  double distance=0;
	  for( size_t d=0 ; d<dimension ; d++ ) 
	    distance += (vertexData[vertexIndex2][d]-vertexData[vertexIndex1][d])*
	      (vertexData[vertexIndex2][d]-vertexData[vertexIndex1][d]);
	  
	  distance=std::sqrt(distance);
	  
	  if (restinglength>0){ 
	    if(variableIndex(0,0)==1) 
	      restinglength+=h*Kgrowth;
	    else if(variableIndex(0,0)==2) 
	      restinglength+=h*Kgrowth*restinglength ;
	    else if(variableIndex(0,0)==3){ 
	      restinglength+=h*Kgrowth*(distance-restinglength) ;
	    }
	    else if(variableIndex(0,0)==4 && distance>restinglength) 
	      restinglength+=h*Kgrowth ;
	    else if(variableIndex(0,0)==5 && distance>restinglength) 
	      restinglength+=h*Kgrowth*restinglength ;
	    else if(variableIndex(0,0)==6 && distance>restinglength) 
	      restinglength+=h*Kgrowth*(distance-restinglength) ;	  
	    else if(variableIndex(0,0)==7) {
	      restinglength += h*KgrowthStress*(distance-restinglength)
		+ h*Kgrowth*restinglength;
	    }
	    else if(variableIndex(0,0)==8) {
	      if(hasSister[vertexIndex1]==1 || hasSister[vertexIndex2]==1)
		restinglength+=h*Kgrowth*restinglength ;
	      else{
		restinglength+=h*KgrowthStress*restinglength;
		//std::cerr<<"broken sisters"<<std::endl;
	      }
	    }
	    else {
	      std::cerr << "VertexFromExternalSpringFromPerpVertexDynamic::update()"
			<< std::endl << "Wrong growth rule index given!"
			<< std::endl;
	      std::exit(EXIT_FAILURE);
	    }
	    connections[cellIndex][vertexCellIndex1][vertexCellIndex2]=restinglength;
	  }					
	}
    }
  } // end of growth rules
  // growth decay
  setParameter(3,parameter(3)*parameter(6));
  //std::cerr << parameter(3) << std::endl;
}

void VertexFromExternalSpringFromPerpVertexDynamic::
printState(Tissue *T,
	   DataMatrix &cellData,
	   DataMatrix &wallData,
	   DataMatrix &vertexData, 
	   std::ostream &os)
{
  static size_t index=0;
  size_t counter=0; 
  size_t dimension = vertexData[0].size();
  size_t numCells = T->numCell();
  for (size_t cellIndex=0 ; cellIndex< numCells ; cellIndex++){
    size_t numCellVertices= T->cell(cellIndex).numVertex();
    for (size_t vertexCellIndex1=0 ; vertexCellIndex1< numCellVertices ; vertexCellIndex1++){
      size_t vertexIndex1 = T->cell(cellIndex).vertex(vertexCellIndex1)->index();
      for (size_t vertexCellIndex2=0 ; vertexCellIndex2< numCellVertices ; vertexCellIndex2++)
	if (vertexCellIndex2!=vertexCellIndex1 && connections[cellIndex][vertexCellIndex1][vertexCellIndex2]!=0){
	  double restingLength=connections[cellIndex][vertexCellIndex1][vertexCellIndex2];
	  size_t vertexIndex2= T->cell(cellIndex).vertex(vertexCellIndex2)->index();
	  double distance=0;
	  for( size_t d=0 ; d<dimension ; d++ )
	    distance += (vertexData[vertexIndex2][d]-vertexData[vertexIndex1][d])*
	      (vertexData[vertexIndex2][d]-vertexData[vertexIndex1][d]);	  
	  distance=std::sqrt(distance);

	  //vertex1
	  os << index << " " << counter << " ";
	  for( size_t d=0 ; d<dimension ; d++ )
	    os << vertexData[vertexIndex1][d] << " ";
	  os << restingLength << " " << distance << " " 
	     << parameter(3)*restingLength+parameter(7)*(distance-restingLength) << std::endl;
	  //vertex2
	  os << index << " " << counter << " ";
	  for( size_t d=0 ; d<dimension ; d++ )
	    os << vertexData[vertexIndex2][d] << " ";
	  os << restingLength << " " << distance
	     << parameter(3)*restingLength+parameter(7)*(distance-restingLength) << std::endl;
	  os << std::endl;
	  counter++;
	}
    }
  }
  index++;
}









