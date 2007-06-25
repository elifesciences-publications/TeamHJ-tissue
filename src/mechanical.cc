//
// Filename     : mechanical.cc
// Description  : Classes describing updates due to mechanical interaction
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : April 2006
// Revision     : $Id:$
//
#include <utility>
#include "baseReaction.h"
#include "mechanical.h"
#include "tissue.h"

//!Constructor for the SpringAsymmetric class
VertexFromWallSpringAsymmetric::
VertexFromWallSpringAsymmetric(std::vector<double> &paraValue, 
															 std::vector< std::vector<size_t> > 
															 &indValue ) 
{  
  // Do some checks on the parameters and variable indeces
  if( paraValue.size()!=2 ) {
    std::cerr << "VertexFromWallSpringAsymmetric::"
							<< "VertexFromWallSpringAsymmetric() "
							<< "Uses two parameters K_force frac_adhesion.\n";
    exit(0);
  }
  if( indValue.size() < 1 || indValue.size() > 2 
			|| indValue[0].size() != 1 
			|| (indValue.size()==2 && indValue[1].size() != 1) ) {
    std::cerr << "VertexFromWallSpringAsymmetric::"
							<< "VertexFromWallSpringAsymmetric() "
							<< "Wall length index given in first level,"
							<< " and optionally wall variable save index in second.\n";
    exit(0);
  }

  // Set the variable values
  setId("VertexFromWallSpringAsymmetric");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  // Set the parameter identities
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "K_force";
  tmp[1] = "frac_adh";
  setParameterId( tmp );
}

//! Derivative contribution for asymmetric wall springs on vertices
/*! 
*/
void VertexFromWallSpringAsymmetric::
derivs(Tissue &T,
       std::vector< std::vector<double> > &cellData,
       std::vector< std::vector<double> > &wallData,
       std::vector< std::vector<double> > &vertexData,
       std::vector< std::vector<double> > &cellDerivs,
       std::vector< std::vector<double> > &wallDerivs,
       std::vector< std::vector<double> > &vertexDerivs ) {
  
  //Do the update for each wall
  size_t numWalls = T.numWall();
  size_t wallLengthIndex = variableIndex(0,0);
  
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
       std::vector< std::vector<double> > &cellData,
       std::vector< std::vector<double> > &wallData,
       std::vector< std::vector<double> > &vertexData,
       std::vector< std::vector<double> > &cellDerivs,
       std::vector< std::vector<double> > &wallDerivs,
       std::vector< std::vector<double> > &vertexDerivs ) 
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
				c1Fac += n_c1[d]*n_w[d]/(c1Norm*distance);
		}
		else
			c1Fac = 0.5;
		if( T.wall(i).cell2() != T.background() &&
				cellData[T.wall(i).cell2()->index()][directionIndex+dimension]>0.5 ) {
			for( size_t d=0 ; d<dimension ; d++ )		
				c2Fac += n_c2[d]*n_w[d]/(c2Norm*distance);
		}
		else
			c2Fac = 0.5;
		
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

void VertexFromWallSpringMT::initiate(Tissue &T)
{
}

VertexFromEpidermalWallSpringAsymmetric::
VertexFromEpidermalWallSpringAsymmetric(std::vector<double> &paraValue, 
					std::vector< std::vector<size_t> > 
					&indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //////////////////////////////////////////////////////////////////////
  if( paraValue.size()!=2 ) {
    std::cerr << "VertexFromEpidermalWallSpringAsymmetric::"
	      << "VertexFromEpidermalWallSpringAsymmetric() "
	      << "Uses two parameters K_force frac_adhesion.\n";
    exit(0);
  }

  if( indValue.size() < 1 || indValue.size() > 2 
			|| indValue[0].size() != 1 
			|| (indValue.size()==2 && indValue[1].size() != 1) ) {
    std::cerr << "VertexFromEpidermalWallSpringAsymmetric::"
							<< "VertexFromEpidermalWallSpringAsymmetric() "
							<< "Wall length index given in first level,"
							<< " and optionally wall variable save index in second.\n";
    exit(0);
  }
  //Set the variable values
  //////////////////////////////////////////////////////////////////////
  setId("VertexFromEpidermalWallSpringAsymmetric");
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
void VertexFromEpidermalWallSpringAsymmetric::
derivs(Tissue &T,
       std::vector< std::vector<double> > &cellData,
       std::vector< std::vector<double> > &wallData,
       std::vector< std::vector<double> > &vertexData,
       std::vector< std::vector<double> > &cellDerivs,
       std::vector< std::vector<double> > &wallDerivs,
       std::vector< std::vector<double> > &vertexDerivs ) {
  
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
  }
}

VertexFromEpidermalCellWallSpringAsymmetric::
VertexFromEpidermalCellWallSpringAsymmetric(std::vector<double> &paraValue, 
																						std::vector< std::vector<size_t> > 
																						&indValue ) 
{  
  //Do some checks on the parameters and variable indeces
  //////////////////////////////////////////////////////////////////////
  if( paraValue.size()!=2 ) {
    std::cerr << "VertexFromEpidermalCellWallSpringAsymmetric::"
							<< "VertexFromEpidermalCellWallSpringAsymmetric() "
							<< "Uses two parameters K_force frac_adhesion.\n";
    exit(0);
  }
  if( indValue.size() < 1 || indValue.size() > 2 
			|| indValue[0].size() != 1 
			|| (indValue.size()==2 && indValue[1].size() != 1) ) {
    std::cerr << "VertexFromEpidermalCellWallSpringAsymmetric::"
							<< "VertexFromEpidermalCellWallSpringAsymmetric() "
							<< "Wall length index given in first level,"
							<< " and optionally wall variable save index in second.\n";
    exit(0);
  }
  //Set the variable values
  //////////////////////////////////////////////////////////////////////
  setId("VertexFromEpidermalCellWallSpringAsymmetric");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //////////////////////////////////////////////////////////////////////
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "K_force";
  tmp[1] = "frac_adh";
  setParameterId( tmp );
}

void VertexFromEpidermalCellWallSpringAsymmetric::
derivs(Tissue &T,
       std::vector< std::vector<double> > &cellData,
       std::vector< std::vector<double> > &wallData,
       std::vector< std::vector<double> > &vertexData,
       std::vector< std::vector<double> > &cellDerivs,
       std::vector< std::vector<double> > &wallDerivs,
       std::vector< std::vector<double> > &vertexDerivs ) 
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
  }
}

//!Constructor for the SpringAsymmetric class
VertexFromCellPowerdiagram::
VertexFromCellPowerdiagram(std::vector<double> &paraValue, 
			   std::vector< std::vector<size_t> > 
			   &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //////////////////////////////////////////////////////////////////////
  if( paraValue.size()!=1 ) {
    std::cerr << "VertexFromCellPowerdiagram::"
	      << "VertexFromCellPowerdiagram() "
	      << "Uses one parameter K_force.\n";
    exit(0);
  }
  if( indValue.size() != 1 || indValue[0].size() != 1 ) {
    std::cerr << "VertexFromCellPowerdiagram::"
	      << "VertexFromCellPowerdiagram() "
	      << "Cell radius index given.\n";
    exit(0);
  }
  //Set the variable values
  //////////////////////////////////////////////////////////////////////
  setId("VertexFromCellPowerdiagram");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //////////////////////////////////////////////////////////////////////
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "K_force";
  setParameterId( tmp );
}

//! Derivative contribution for asymmetric wall springs on vertices
/*! 
*/
void VertexFromCellPowerdiagram::
derivs(Tissue &T,
       std::vector< std::vector<double> > &cellData,
       std::vector< std::vector<double> > &wallData,
       std::vector< std::vector<double> > &vertexData,
       std::vector< std::vector<double> > &cellDerivs,
       std::vector< std::vector<double> > &wallDerivs,
       std::vector< std::vector<double> > &vertexDerivs ) {
  
  //Do the update for each vertex
  //size_t numCells = T.numCell();
  size_t numVertices = T.numVertex();
  size_t cellRIndex = variableIndex(0,0);
  size_t dimension = T.vertex(0).numPosition(); 
  
  for( size_t i=0 ; i<numVertices ; ++i ) {
    //Only calculate if vertex connected to three cells
    size_t numCellForVertex = T.vertex(i).numCell();
    if(  numCellForVertex == 3 ) {
      //Calculate position for all three cells
      std::vector< std::vector<double> > cellPos(numCellForVertex);
      std::vector<double> cellR(numCellForVertex);
      for( size_t cellI=0 ; cellI<numCellForVertex ; ++cellI ) {
	cellR[cellI] = cellData[ T.vertex(i).cell(cellI)->index() ][ cellRIndex ];
	//cellR[cellI] = T.cell( T.vertex(i).cell(cellI)->index() ).variable(0);
	cellPos[cellI] = T.cell( T.vertex(i).cell(cellI)->index() ).
	  positionFromVertex(vertexData);
      }
      //std::cerr << "* " << cellPos[0][0] << " " << cellPos[0][1] << "  "
      //	<< cellPos[1][0] << " " << cellPos[1][1] << "  "
      //	<< cellPos[2][0] << " " << cellPos[2][1] << std::endl;

      //Calculate optimal position from cell positions and raddii
      //according to power diagram formulation
      if( dimension != 2 ) {
	std::cerr << "VertexFromCellPowerdiagram::derivs() Only defined "
		  << "for two dimensions so far..." << std::endl;
	exit(-1);
      }
      double Kji = cellPos[1][0]*cellPos[1][0]+cellPos[1][1]*cellPos[1][1]-
	cellPos[0][0]*cellPos[0][0]-cellPos[0][1]*cellPos[0][1]-
	(cellR[1]*cellR[1]-cellR[0]*cellR[0]);
      double Kki = cellPos[2][0]*cellPos[2][0]+cellPos[2][1]*cellPos[2][1]-
	cellPos[0][0]*cellPos[0][0]-cellPos[0][1]*cellPos[0][1]-
	(cellR[2]*cellR[2]-cellR[0]*cellR[0]);
      std::vector<double> powPos(dimension);
      powPos[1] = 0.5*( Kji*(cellPos[2][0]-cellPos[0][0]) -
			Kki*(cellPos[1][0]-cellPos[0][0]) ) /
	( (cellPos[1][1]-cellPos[0][1])*(cellPos[2][0]-cellPos[0][0]) -
	  (cellPos[2][1]-cellPos[0][1])*(cellPos[1][0]-cellPos[0][0]) );
      powPos[0] = 0.5*Kji/(cellPos[1][0]-cellPos[0][0]) -
	powPos[1]*(cellPos[1][1]-cellPos[0][1])/(cellPos[1][0]-cellPos[0][0]);
      //std::cerr << i << " " << vertexData[i][0] << " " << vertexData[i][1] << "\t"
      //	<< powPos[0] << " " << powPos[1] << "\t"
      //	<< cellPos[0][0] << " " << cellPos[0][1] << "  "
      //	<< cellPos[1][0] << " " << cellPos[1][1] << "  "
      //	<< cellPos[2][0] << " " << cellPos[2][1] << "\t"
      //	<< cellR[0] << " " << cellR[1] << " " << cellR[2] << std::endl;

      //Update vertex for each dimension
      for(size_t d=0 ; d<dimension ; d++ )
	vertexDerivs[i][d] -= parameter(0)*(vertexData[i][d]-powPos[d]);
    }
  }
}

//!Constructor
VertexFromCellPressure::
VertexFromCellPressure(std::vector<double> &paraValue, 
			   std::vector< std::vector<size_t> > 
			   &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //////////////////////////////////////////////////////////////////////
  if( paraValue.size()!=1 ) {
    std::cerr << "VertexFromCellPressure::"
	      << "VertexFromCellPressure() "
	      << "Uses one parameter K_force.\n";
    exit(0);
  }
  if( indValue.size() != 0 ) {
    std::cerr << "VertexFromCellPressure::"
	      << "VertexFromCellPressure() "
	      << "No index given.\n";
    exit(0);
  }
  //Set the variable values
  //////////////////////////////////////////////////////////////////////
  setId("VertexFromCellPressure");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //////////////////////////////////////////////////////////////////////
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "K_force";
  setParameterId( tmp );
}

//! Derivative contribution for asymmetric wall springs on vertices
/*! 
*/
void VertexFromCellPressure::
derivs(Tissue &T,
       std::vector< std::vector<double> > &cellData,
       std::vector< std::vector<double> > &wallData,
       std::vector< std::vector<double> > &vertexData,
       std::vector< std::vector<double> > &cellDerivs,
       std::vector< std::vector<double> > &wallDerivs,
       std::vector< std::vector<double> > &vertexDerivs ) {
  
  //Do the update for each vertex via each wall in each cell
  size_t numCells = T.numCell();
  size_t dimension = T.vertex(0).numPosition(); 
	
	//Assumming vertices and walls are sorted
	//////////////////////////////////////////////////////////////////////
	assert( dimension==2 );
  //For each cell
  for( size_t cellI=0 ; cellI<numCells ; ++cellI ) {
		
    Cell &tmpCell = T.cell(cellI);
    
    //calculate volume (area)
    double cellVolume=0.0;
    for( size_t k=0 ; k<tmpCell.numVertex() ; ++k ) {
      size_t v1I = tmpCell.vertex(k)->index();
      size_t v2I = tmpCell.vertex((k+1)%(tmpCell.numVertex()))->index();
      cellVolume += vertexData[v1I][0]*vertexData[v2I][1]-
	vertexData[v1I][1]*vertexData[v2I][0];
    }
    cellVolume = 0.5*cellVolume;
    double factor=0.5*parameter(0);
    if( cellVolume<0.0 ) 
      factor =-factor;
    
    for( size_t k=0 ; k<tmpCell.numVertex() ; ++k ) {
      size_t v1I = tmpCell.vertex(k)->index();
      size_t v1PlusI = tmpCell.vertex((k+1)%(tmpCell.numVertex()))->index();
      size_t v1MinusK = k>0 ? k-1 : tmpCell.numVertex()-1;
      size_t v1MinusI = tmpCell.vertex(v1MinusK)->index();
      
      vertexDerivs[v1I][0] += factor*(vertexData[v1PlusI][1]-
				      vertexData[v1MinusI][1]);
      vertexDerivs[v1I][1] += factor*(vertexData[v1MinusI][0]-
				      vertexData[v1PlusI][0]);
    }
  }
  
	//Assuming vertices and walls not sorted (old version)
	//////////////////////////////////////////////////////////////////////
//   //For each cell
//   for( size_t cellI=0 ; cellI<numCells ; ++cellI ) {
    
//     Cell &tmpCell = T.cell(cellI);
//     //Calculate cell position from vertices
//     std::vector<double> xCenter = tmpCell.positionFromVertex(vertexData);
//     assert( xCenter.size()==dimension );
    
//     //Calculate derivative contributions to vertices from each wall
//     for( size_t k=0 ; k<tmpCell.numWall() ; ++k ) {
//       Wall &tmpWall = tmpCell.wallRef(k);
//       size_t v1I = tmpWall.vertex1()->index();
//       size_t v2I = tmpWall.vertex2()->index();
//       std::vector<double> n(dimension),dx(dimension), x0(dimension);
//       double b=0;
//       for( size_t d=0 ; d<dimension ; ++d ) {
// 				n[d] = vertexData[v2I][d]-vertexData[v1I][d];
// 				b += n[d]*n[d];
// 				x0[d] = 0.5*(vertexData[v1I][d] + vertexData[v2I][d]);
// 				dx[d] = xCenter[d]-x0[d];
//       }
//       assert( b>0.0 );
//       b = std::sqrt(b);
//       for( size_t d=0 ; d<dimension ; ++d )
// 				n[d] /= b;
//       double bInv = 1.0/b;
//       double h = dx[0]*dx[0] + dx[1]*dx[1]
// 				-(n[0]*dx[0]+n[1]*dx[1])*(n[0]*dx[0]+n[1]*dx[1]);
//       assert( h>0.0 );
//       h = std::sqrt(h);
//       double hInv = 1.0/h;
//       double fac = parameter(0)*0.5;
      
//       vertexDerivs[v1I][0] += fac*( 0.5*b*hInv*
// 																		( -dx[0] - 2.0*(n[0]*dx[0]+n[1]*dx[1])
// 																			*(-n[1]*n[1]*bInv*dx[0]
// 																				-0.5*n[0]
// 																				+n[0]*n[1]*bInv*dx[1]) ) 
// 																		- h*n[0] ); 
//       vertexDerivs[v2I][0] += fac*( 0.5*b*hInv*
// 																		( -dx[0] - 2.0*(n[0]*dx[0]+n[1]*dx[1])
// 																			*(n[1]*n[1]*bInv*dx[0]
// 																				-0.5*n[0]
// 																				-n[0]*n[1]*bInv*dx[1]) )
// 																		+ h*n[0] );
//       vertexDerivs[v1I][1] += fac*( 0.5*b*hInv*
// 																		( -dx[1] - 2.0*(n[0]*dx[0]+n[1]*dx[1])
// 																			*(-n[0]*n[0]*bInv*dx[1]
// 																				-0.5*n[1]
// 																				+n[0]*n[1]*bInv*dx[0]) ) 
// 																		- h*n[1] ); 
//       vertexDerivs[v2I][1] += fac*( 0.5*b*hInv*
// 																		( -dx[1] - 2.0*(n[0]*dx[0]+n[1]*dx[1])
// 																			*(n[0]*n[0]*bInv*dx[1]
// 																				-0.5*n[1]
// 																				-n[0]*n[1]*bInv*dx[0]) )
// 																		+ h*n[1] );
      
//       //vertexDerivs[v1I][0] += fac*( b*hInv*(bInv*(n[0]*dx[0]+n[1]*dx[1])*
//       //				    (n[1]*n[1]*dx[0]+b*n[0]-n[0]*n[1]*dx[1]) + dx[0] ) - h*n[0] );
//       //vertexDerivs[v1I][1] += fac*( b*hInv*(bInv*(n[1]*dx[1]+n[0]*dx[0])*
//       //				    (-n[0]*n[0]*dx[1]+b*n[1]+n[1]*n[0]*dx[0]) + dx[1] ) - h*n[1] );
//       //vertexDerivs[v2I][0] += fac*( n[0]*h + n[1]*hInv*(n[0]*dx[0]+n[1]*dx[1])*(n[0]*dx[1]-n[1]*dx[0]) );
//       //vertexDerivs[v2I][1] += fac*( n[1]*h + n[0]*hInv*(n[1]*dx[1]+n[0]*dx[0])*(n[1]*dx[0]-n[0]*dx[1]) );      
//     }
//   }
}

//!Constructor
VertexFromCellPressureVolumeNormalized::
VertexFromCellPressureVolumeNormalized(std::vector<double> &paraValue, 
			   std::vector< std::vector<size_t> > 
			   &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //////////////////////////////////////////////////////////////////////
  if( paraValue.size()!=1 ) {
    std::cerr << "VertexFromCellPressureVolumeNormalized::"
	      << "VertexFromCellPressureVolumeNormalized() "
	      << "Uses one parameter K_force.\n";
    exit(0);
  }
  if( indValue.size() != 0 ) {
    std::cerr << "VertexFromCellPressureVolumeNormalized::"
	      << "VertexFromCellPressureVolumeNormalized() "
	      << "No index given.\n";
    exit(0);
  }
  //Set the variable values
  //////////////////////////////////////////////////////////////////////
  setId("VertexFromCellPressureVolumeNormalized");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //////////////////////////////////////////////////////////////////////
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "K_force";
  setParameterId( tmp );
}

//! Derivative contribution for asymmetric wall springs on vertices
/*! 
 */
void VertexFromCellPressureVolumeNormalized::
derivs(Tissue &T,
       std::vector< std::vector<double> > &cellData,
       std::vector< std::vector<double> > &wallData,
       std::vector< std::vector<double> > &vertexData,
       std::vector< std::vector<double> > &cellDerivs,
       std::vector< std::vector<double> > &wallDerivs,
       std::vector< std::vector<double> > &vertexDerivs ) {
  
  //Do the update for each vertex via each wall in each cell
  size_t numCells = T.numCell();
  size_t dimension = T.vertex(0).numPosition(); 
	
  //For each cell
  for( size_t cellI=0 ; cellI<numCells ; ++cellI ) {
    
    Cell &tmpCell = T.cell(cellI);
    //Calculate cell position from vertices
    std::vector<double> xCenter = tmpCell.positionFromVertex(vertexData);
    assert( xCenter.size()==dimension );
    double cellVolume = tmpCell.calculateVolume(vertexData);
		
    //Calculate derivative contributions to vertices from each wall
    for( size_t k=0 ; k<tmpCell.numWall() ; ++k ) {
      Wall &tmpWall = tmpCell.wallRef(k);
      size_t v1I = tmpWall.vertex1()->index();
      size_t v2I = tmpWall.vertex2()->index();
      std::vector<double> n(dimension),dx(dimension), x0(dimension);
      double b=0;
      for( size_t d=0 ; d<dimension ; ++d ) {
				n[d] = vertexData[v2I][d]-vertexData[v1I][d];
				b += n[d]*n[d];
				x0[d] = 0.5*(vertexData[v1I][d] + vertexData[v2I][d]);
				dx[d] = xCenter[d]-x0[d];
      }
      assert( b>0.0 );
      b = std::sqrt(b);
      for( size_t d=0 ; d<dimension ; ++d )
				n[d] /= b;
      double bInv = 1.0/b;
      double h = dx[0]*dx[0] + dx[1]*dx[1]
				-(n[0]*dx[0]+n[1]*dx[1])*(n[0]*dx[0]+n[1]*dx[1]);
      assert( h>0.0 );
      h = std::sqrt(h);
      double hInv = 1.0/h;
      double fac = parameter(0)*0.5/cellVolume;
      
      vertexDerivs[v1I][0] += fac*( 0.5*b*hInv*
																		( -dx[0] - 2.0*(n[0]*dx[0]+n[1]*dx[1])
																			*(-n[1]*n[1]*bInv*dx[0]
																				-0.5*n[0]
																				+n[0]*n[1]*bInv*dx[1]) ) 
																		- h*n[0] ); 
      vertexDerivs[v2I][0] += fac*( 0.5*b*hInv*
																		( -dx[0] - 2.0*(n[0]*dx[0]+n[1]*dx[1])
																			*(n[1]*n[1]*bInv*dx[0]
																				-0.5*n[0]
																				-n[0]*n[1]*bInv*dx[1]) )
																		+ h*n[0] );
      vertexDerivs[v1I][1] += fac*( 0.5*b*hInv*
																		( -dx[1] - 2.0*(n[0]*dx[0]+n[1]*dx[1])
																			*(-n[0]*n[0]*bInv*dx[1]
																				-0.5*n[1]
																				+n[0]*n[1]*bInv*dx[0]) ) 
																		- h*n[1] ); 
      vertexDerivs[v2I][1] += fac*( 0.5*b*hInv*
																		( -dx[1] - 2.0*(n[0]*dx[0]+n[1]*dx[1])
																			*(n[0]*n[0]*bInv*dx[1]
																				-0.5*n[1]
																				-n[0]*n[1]*bInv*dx[0]) )
																		+ h*n[1] );      
    }
  }
}

//!Constructor
VertexFromCellPressureThresholdFromMaxPos::
VertexFromCellPressureThresholdFromMaxPos(std::vector<double> &paraValue, 
			   std::vector< std::vector<size_t> > 
			   &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //////////////////////////////////////////////////////////////////////
  if( paraValue.size()!=2 ) {
    std::cerr << "VertexFromCellPressureThresholdFromMaxPos::"
	      << "VertexFromCellPressureThresholdFromMaxPos() "
	      << "Uses two parameters K_force and X_th.\n";
    exit(0);
  }
  if( indValue.size() != 1 || indValue[0].size() != 1 ) {
    std::cerr << "VertexFromCellPressureThresholdFromMaxPos::"
	      << "VertexFromCellPressureThresholdFromMaxPos() "
	      << "One index given (direction).\n";
    exit(0);
  }
  //Set the variable values
  //////////////////////////////////////////////////////////////////////
  setId("VertexFromCellPressureThresholdFromMaxPos");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //////////////////////////////////////////////////////////////////////
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "K_force";
  tmp[1] = "X_th";
  setParameterId( tmp );
}

//! Derivative contribution for asymmetric wall springs on vertices
/*! 
*/
void VertexFromCellPressureThresholdFromMaxPos::
derivs(Tissue &T,
       std::vector< std::vector<double> > &cellData,
       std::vector< std::vector<double> > &wallData,
       std::vector< std::vector<double> > &vertexData,
       std::vector< std::vector<double> > &cellDerivs,
       std::vector< std::vector<double> > &wallDerivs,
       std::vector< std::vector<double> > &vertexDerivs ) {
  
  //Do the update for each vertex via each wall in each cell
  size_t numCells = T.numCell();
  size_t dimension = T.vertex(0).numPosition(); 

  //For each cell
  for( size_t cellI=0 ; cellI<numCells ; ++cellI ) {
    
    Cell &tmpCell = T.cell(cellI);
    //Calculate cell position from vertices
    std::vector<double> xCenter = tmpCell.positionFromVertex(vertexData);
    assert( xCenter.size()==dimension );
    assert( variableIndex(0,0)<dimension );
		//Find max pos in given direction;
		double max=vertexData[0][variableIndex(0,0)];
		for (size_t i=1; i<vertexData.size(); ++i ) {
			if ( vertexData[i][variableIndex(0,0)]>max )
				max = vertexData[i][variableIndex(0,0)];
		}
		//Only if close to apex
    if ( max-xCenter[variableIndex(0,0)]<parameter(1) ) {
			
			//Calculate derivative contributions to vertices from each wall
			for( size_t k=0 ; k<tmpCell.numWall() ; ++k ) {
				Wall &tmpWall = tmpCell.wallRef(k);
				size_t v1I = tmpWall.vertex1()->index();
				size_t v2I = tmpWall.vertex2()->index();
				std::vector<double> n(dimension),dx(dimension), x0(dimension);
				double b=0;
				for( size_t d=0 ; d<dimension ; ++d ) {
					n[d] = vertexData[v2I][d]-vertexData[v1I][d];
					b += n[d]*n[d];
					x0[d] = 0.5*(vertexData[v1I][d] + vertexData[v2I][d]);
					dx[d] = xCenter[d]-x0[d];
				}
				assert( b>0.0 );
				b = std::sqrt(b);
				for( size_t d=0 ; d<dimension ; ++d )
					n[d] /= b;
				double bInv = 1.0/b;
				double h = dx[0]*dx[0] + dx[1]*dx[1]
					-(n[0]*dx[0]+n[1]*dx[1])*(n[0]*dx[0]+n[1]*dx[1]);
				assert( h>0.0 );
				h = std::sqrt(h);
				double hInv = 1.0/h;
				double fac = parameter(0)*0.5;
				
				vertexDerivs[v1I][0] += fac*( 0.5*b*hInv*
																			( -dx[0] - 2.0*(n[0]*dx[0]+n[1]*dx[1])
																				*(-n[1]*n[1]*bInv*dx[0]
																					-0.5*n[0]
																					+n[0]*n[1]*bInv*dx[1]) ) 
																			- h*n[0] ); 
				vertexDerivs[v2I][0] += fac*( 0.5*b*hInv*
																			( -dx[0] - 2.0*(n[0]*dx[0]+n[1]*dx[1])
																				*(n[1]*n[1]*bInv*dx[0]
																					-0.5*n[0]
																					-n[0]*n[1]*bInv*dx[1]) )
																			+ h*n[0] );
				vertexDerivs[v1I][1] += fac*( 0.5*b*hInv*
																			( -dx[1] - 2.0*(n[0]*dx[0]+n[1]*dx[1])
																				*(-n[0]*n[0]*bInv*dx[1]
																					-0.5*n[1]
																					+n[0]*n[1]*bInv*dx[0]) ) 
																			- h*n[1] ); 
				vertexDerivs[v2I][1] += fac*( 0.5*b*hInv*
																			( -dx[1] - 2.0*(n[0]*dx[0]+n[1]*dx[1])
																				*(n[0]*n[0]*bInv*dx[1]
																					-0.5*n[1]
																					-n[0]*n[1]*bInv*dx[0]) )
																			+ h*n[1] );
				
				//vertexDerivs[v1I][0] += fac*( b*hInv*(bInv*(n[0]*dx[0]+n[1]*dx[1])*
				//				    (n[1]*n[1]*dx[0]+b*n[0]-n[0]*n[1]*dx[1]) + dx[0] ) - h*n[0] );
				//vertexDerivs[v1I][1] += fac*( b*hInv*(bInv*(n[1]*dx[1]+n[0]*dx[0])*
				//				    (-n[0]*n[0]*dx[1]+b*n[1]+n[1]*n[0]*dx[0]) + dx[1] ) - h*n[1] );
				//vertexDerivs[v2I][0] += fac*( n[0]*h + n[1]*hInv*(n[0]*dx[0]+n[1]*dx[1])*(n[0]*dx[1]-n[1]*dx[0]) );
				//vertexDerivs[v2I][1] += fac*( n[1]*h + n[0]*hInv*(n[1]*dx[1]+n[0]*dx[0])*(n[1]*dx[0]-n[0]*dx[1]) );      
			}
    }
  }
}

//!Constructor
VertexFromCellInternalPressure::
VertexFromCellInternalPressure(std::vector<double> &paraValue, 
			   std::vector< std::vector<size_t> > 
			   &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //////////////////////////////////////////////////////////////////////
  if( paraValue.size()!=1 ) {
    std::cerr << "VertexFromCellInternalPressure::"
	      << "VertexFromCellInternalPressure() "
	      << "Uses one parameter K_force.\n";
    exit(0);
  }
  if( indValue.size() != 0 ) {
    std::cerr << "VertexFromCellInternalPressure::"
	      << "VertexFromCellInternalPressure() "
	      << "No index given.\n";
    exit(0);
  }
  //Set the variable values
  //////////////////////////////////////////////////////////////////////
  setId("VertexFromCellInternalPressure");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //////////////////////////////////////////////////////////////////////
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "K_force";
  setParameterId( tmp );
}

//! Derivative contribution for asymmetric wall springs on vertices
/*! 
*/
void VertexFromCellInternalPressure::
derivs(Tissue &T,
       std::vector< std::vector<double> > &cellData,
       std::vector< std::vector<double> > &wallData,
       std::vector< std::vector<double> > &vertexData,
       std::vector< std::vector<double> > &cellDerivs,
       std::vector< std::vector<double> > &wallDerivs,
       std::vector< std::vector<double> > &vertexDerivs ) {
  
  //Do the update for each vertex via each wall in each cell
  size_t numCells = T.numCell();
  size_t dimension = T.vertex(0).numPosition(); 

  //For each cell
  for( size_t cellI=0 ; cellI<numCells ; ++cellI ) {
    Cell &tmpCell = T.cell(cellI);
    if( !(tmpCell.isNeighbor(T.background())) ) {
      //Calculate cell position from vertices
      std::vector<double> xCenter = tmpCell.positionFromVertex(vertexData);
      assert( xCenter.size()==dimension );
      
      //Calculate derivative contributions to vertices from each wall
      for( size_t k=0 ; k<tmpCell.numWall() ; ++k ) {
	Wall &tmpWall = tmpCell.wallRef(k);
	size_t v1I = tmpWall.vertex1()->index();
	size_t v2I = tmpWall.vertex2()->index();
	std::vector<double> n(dimension),dx(dimension), x0(dimension);
	double b=0;
	for( size_t d=0 ; d<dimension ; ++d ) {
	  n[d] = vertexData[v2I][d]-vertexData[v1I][d];
	  b += n[d]*n[d];
	  x0[d] = 0.5*(vertexData[v1I][d] + vertexData[v2I][d]);
	  dx[d] = xCenter[d]-x0[d];
	}
	assert( b>0.0 );
	b = std::sqrt(b);
	for( size_t d=0 ; d<dimension ; ++d )
	  n[d] /= b;
	double bInv = 1.0/b;
	double h = dx[0]*dx[0] + dx[1]*dx[1]
	  -(n[0]*dx[0]+n[1]*dx[1])*(n[0]*dx[0]+n[1]*dx[1]);
	assert( h>0.0 );
	h = std::sqrt(h);
	double hInv = 1.0/h;
	double fac = parameter(0)*0.5;
	
	vertexDerivs[v1I][0] += fac*( 0.5*b*hInv*
				      ( -dx[0] - 2.0*(n[0]*dx[0]+n[1]*dx[1])
					*(-n[1]*n[1]*bInv*dx[0]
					  -0.5*n[0]
					  +n[0]*n[1]*bInv*dx[1]) ) 
				      - h*n[0] ); 
	vertexDerivs[v2I][0] += fac*( 0.5*b*hInv*
				      ( -dx[0] - 2.0*(n[0]*dx[0]+n[1]*dx[1])
					*(n[1]*n[1]*bInv*dx[0]
					  -0.5*n[0]
					  -n[0]*n[1]*bInv*dx[1]) )
				      + h*n[0] );
	vertexDerivs[v1I][1] += fac*( 0.5*b*hInv*
				      ( -dx[1] - 2.0*(n[0]*dx[0]+n[1]*dx[1])
					*(-n[0]*n[0]*bInv*dx[1]
					  -0.5*n[1]
					  +n[0]*n[1]*bInv*dx[0]) ) 
				      - h*n[1] ); 
	vertexDerivs[v2I][1] += fac*( 0.5*b*hInv*
				      ( -dx[1] - 2.0*(n[0]*dx[0]+n[1]*dx[1])
					*(n[0]*n[0]*bInv*dx[1]
					  -0.5*n[1]
					  -n[0]*n[1]*bInv*dx[0]) )
				      + h*n[1] );
	
	//vertexDerivs[v1I][0] += fac*( b*hInv*(bInv*(n[0]*dx[0]+n[1]*dx[1])*
	//				    (n[1]*n[1]*dx[0]+b*n[0]-n[0]*n[1]*dx[1]) + dx[0] ) - h*n[0] );
	//vertexDerivs[v1I][1] += fac*( b*hInv*(bInv*(n[1]*dx[1]+n[0]*dx[0])*
	//				    (-n[0]*n[0]*dx[1]+b*n[1]+n[1]*n[0]*dx[0]) + dx[1] ) - h*n[1] );
	//vertexDerivs[v2I][0] += fac*( n[0]*h + n[1]*hInv*(n[0]*dx[0]+n[1]*dx[1])*(n[0]*dx[1]-n[1]*dx[0]) );
	//vertexDerivs[v2I][1] += fac*( n[1]*h + n[0]*hInv*(n[1]*dx[1]+n[0]*dx[0])*(n[1]*dx[0]-n[0]*dx[1]) );      
      }
    }
  }
}

//!Constructor
VertexForceOrigoFromIndex::
VertexForceOrigoFromIndex(std::vector<double> &paraValue, 
			  std::vector< std::vector<size_t> > 
			  &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //////////////////////////////////////////////////////////////////////
  if( paraValue.size()!=2 ) {
    std::cerr << "VertexForceOrigoFromIndex::"
	      << "VertexForceOrigoFromIndex() "
	      << "Uses two parameters K_force direction(-1 -> inwards)" 
	      << std::endl;
    exit(0);
  }
  if( paraValue[1] != 1.0 && paraValue[1] != -1.0 ) {
    std::cerr << "VertexForceOrigoFromIndex::"
	      << "VertexForceOrigoFromIndex() "
	      << "direction (second parameter) need to be 1 (outward) "
	      << "or -1 (inwards)." << std::endl;
    exit(0);
  }
  if( indValue.size() != 1 || indValue[0].size()<1 ) {
    std::cerr << "VertexForceOrigoFromIndex::"
	      << "VertexForceOrigoFromIndex() "
	      << "Vertex indices in first level." << std::endl;
    exit(0);
  }
  //Set the variable values
  //////////////////////////////////////////////////////////////////////
  setId("VertexForceOrigoFromIndex");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //////////////////////////////////////////////////////////////////////
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "K_force";
  tmp[1] = "direction";
  setParameterId( tmp );
}

//!Derivative contribution for force towards or from origo
/*! 
 */
void VertexForceOrigoFromIndex::
derivs(Tissue &T,
       std::vector< std::vector<double> > &cellData,
       std::vector< std::vector<double> > &wallData,
       std::vector< std::vector<double> > &vertexData,
       std::vector< std::vector<double> > &cellDerivs,
       std::vector< std::vector<double> > &wallDerivs,
       std::vector< std::vector<double> > &vertexDerivs ) {
  
  size_t dimension = vertexData[variableIndex(0,0)].size();
  double coeff = parameter(0)*parameter(1);
  //For each vertex in list
  for( size_t k=0 ; k<numVariableIndex(0) ; ++k ) {
    size_t i=variableIndex(0,k);
    for( size_t d=0 ; d<dimension ; d++ )
      vertexDerivs[i][d] -= coeff*vertexData[i][d];
  }
}

//!Constructor
CellForceOrigoFromIndex::
CellForceOrigoFromIndex(std::vector<double> &paraValue, 
			  std::vector< std::vector<size_t> > 
			  &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //////////////////////////////////////////////////////////////////////
  if( paraValue.size()!=2 ) {
    std::cerr << "CellForceOrigoFromIndex::"
	      << "CellForceOrigoFromIndex() "
	      << "Uses two parameters K_force direction(-1 -> inwards)" 
	      << std::endl;
    exit(0);
  }
  if( paraValue[1] != 1.0 && paraValue[1] != -1.0 ) {
    std::cerr << "CellForceOrigoFromIndex::"
	      << "CellForceOrigoFromIndex() "
	      << "direction (second parameter) needs to be 1 (outward) "
	      << "or -1 (inwards)." << std::endl;
    exit(0);
  }
  if( indValue.size() != 1 || indValue[0].size()<1 ) {
    std::cerr << "CellForceOrigoFromIndex::"
	      << "CellForceOrigoFromIndex() "
	      << "Cell indices in first level." << std::endl;
    exit(0);
  }
  //Set the variable values
  //////////////////////////////////////////////////////////////////////
  setId("CellForceOrigoFromIndex");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //////////////////////////////////////////////////////////////////////
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "K_force";
  tmp[1] = "direction";
  setParameterId( tmp );
}

//!Derivative contribution for force towards or from origo
/*! 
 */
void CellForceOrigoFromIndex::
derivs(Tissue &T,
       std::vector< std::vector<double> > &cellData,
       std::vector< std::vector<double> > &wallData,
       std::vector< std::vector<double> > &vertexData,
       std::vector< std::vector<double> > &cellDerivs,
       std::vector< std::vector<double> > &wallDerivs,
       std::vector< std::vector<double> > &vertexDerivs ) {
  
  double coeff = parameter(0)*parameter(1);
  size_t dimension = vertexData[0].size();
  //For each vertex in cell list
  for( size_t k=0 ; k<numVariableIndex(0) ; ++k ) {
    size_t cellI=variableIndex(0,k);
    for( size_t l=0 ; l<T.cell(cellI).numVertex() ; ++l ) {
      size_t i = T.cell(cellI).vertex(l)->index();      
      for( size_t d=0 ; d<dimension ; d++ )
	vertexDerivs[i][d] += coeff*vertexData[i][d];
    }
  }
}

//!Constructor
CylinderForce::
CylinderForce(std::vector<double> &paraValue, 
							std::vector< std::vector<size_t> > 
							&indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //////////////////////////////////////////////////////////////////////
  if( paraValue.size()!=2 ) {
    std::cerr << "CylinderForce::"
							<< "CylinderForce() "
							<< "Uses two parameters K_force direction(-1 -> inwards)" 
							<< std::endl;
    exit(0);
  }
  if( paraValue[1] != 1.0 && paraValue[1] != -1.0 ) {
    std::cerr << "CylinderForce::"
							<< "CylinderForce() "
							<< "direction (second parameter) needs to be 1 (outward) "
							<< "or -1 (inwards)." << std::endl;
    exit(0);
  }
  if( indValue.size() != 0 ) {
    std::cerr << "CylinderForce::"
							<< "CylinderForce() "
							<< "No indices used." << std::endl;
    exit(0);
  }
  //Set the variable values
  //////////////////////////////////////////////////////////////////////
  setId("CylinderForce");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //////////////////////////////////////////////////////////////////////
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "K_force";
  tmp[1] = "direction";
  setParameterId( tmp );
}

//!Derivative contribution for force towards or from origo
/*! 
 */
void CylinderForce::
derivs(Tissue &T,
       std::vector< std::vector<double> > &cellData,
       std::vector< std::vector<double> > &wallData,
       std::vector< std::vector<double> > &vertexData,
       std::vector< std::vector<double> > &cellDerivs,
       std::vector< std::vector<double> > &wallDerivs,
       std::vector< std::vector<double> > &vertexDerivs ) {
  
  double coeff = parameter(0)*parameter(1);
  size_t dimension = vertexData[0].size();
  size_t lastPosIndex = dimension-1;
  //For each vertex 
  for( size_t i=0 ; i<T.numVertex() ; ++i ) {
		//On cylinder
		double norm = 0.0;
		for( size_t d=0 ; d<lastPosIndex ; d++ )
			norm += vertexData[i][d]*vertexData[i][d];
		if( norm>0.0 )
			norm = 1.0/std::sqrt(norm);
		else norm=0.0;
		for( size_t d=0 ; d<lastPosIndex ; d++ )
			vertexDerivs[i][d] += coeff*norm*vertexData[i][d];
	}
}

//!Constructor
SphereCylinderForce::
SphereCylinderForce(std::vector<double> &paraValue, 
			  std::vector< std::vector<size_t> > 
			  &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //////////////////////////////////////////////////////////////////////
  if( paraValue.size()!=2 ) {
    std::cerr << "SphereCylinderForce::"
	      << "SphereCylinderForce() "
	      << "Uses two parameters K_force direction(-1 -> inwards)" 
	      << std::endl;
    exit(0);
  }
  if( paraValue[1] != 1.0 && paraValue[1] != -1.0 ) {
    std::cerr << "SphereCylinderForce::"
	      << "SphereCylinderForce() "
	      << "direction (second parameter) needs to be 1 (outward) "
	      << "or -1 (inwards)." << std::endl;
    exit(0);
  }
  if( indValue.size() != 0 ) {
    std::cerr << "SphereCylinderForce::"
	      << "SphereCylinderForce() "
	      << "Cell indices in first level." << std::endl;
    exit(0);
  }
  //Set the variable values
  //////////////////////////////////////////////////////////////////////
  setId("SphereCylinderForce");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //////////////////////////////////////////////////////////////////////
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "K_force";
  tmp[1] = "direction";
  setParameterId( tmp );
}

//!Derivative contribution for force towards or from origo
/*! 
 */
void SphereCylinderForce::
derivs(Tissue &T,
       std::vector< std::vector<double> > &cellData,
       std::vector< std::vector<double> > &wallData,
       std::vector< std::vector<double> > &vertexData,
       std::vector< std::vector<double> > &cellDerivs,
       std::vector< std::vector<double> > &wallDerivs,
       std::vector< std::vector<double> > &vertexDerivs ) {
  
  double coeff = parameter(0)*parameter(1);
  size_t dimension = vertexData[0].size();
  size_t lastPosIndex = dimension-1;
  //For each vertex 
  for( size_t i=0 ; i<T.numVertex() ; ++i ) {
    if( vertexData[i][lastPosIndex]>0.0 ) {
      //On sphere
      double norm = 0.0;
      for( size_t d=0 ; d<dimension ; d++ )
	norm += vertexData[i][d]*vertexData[i][d];
      if( norm>0.0 )
	norm = 1.0/std::sqrt(norm);
      else norm=0.0;
      for( size_t d=0 ; d<dimension ; d++ )
	vertexDerivs[i][d] += coeff*norm*vertexData[i][d];
    }
    else {
      //On cylinder
      double norm = 0.0;
      for( size_t d=0 ; d<dimension-1 ; d++ )
	norm += vertexData[i][d]*vertexData[i][d];
      if( norm>0.0 )
	norm = 1.0/std::sqrt(norm);
      else norm=0.0;
      for( size_t d=0 ; d<dimension-1 ; d++ )
	vertexDerivs[i][d] += coeff*norm*vertexData[i][d];
    }
  }
}

//!Constructor
SphereCylinderForceFromRadius::
SphereCylinderForceFromRadius(std::vector<double> &paraValue, 
			      std::vector< std::vector<size_t> > 
			      &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //////////////////////////////////////////////////////////////////////
  if( paraValue.size()!=3 ) {
    std::cerr << "SphereCylinderForceFromRadius::"
	      << "SphereCylinderForceFromRadius() "
	      << "Uses three parameters F_out, F_in and radius." 
	      << std::endl;
    exit(0);
  }
  if( paraValue[0] < 0.0 || paraValue[1] < 0.0 || paraValue[2] < 0.0 ) {
    std::cerr << "SphereCylinderForceFromRadius::"
	      << "SphereCylinderForceFromRadius() "
	      << "Parameters need to be positive." << std::endl;
    exit(0);
  }
  if( indValue.size() != 0 ) {
    std::cerr << "SphereCylinderForceFromRadius::"
	      << "SphereCylinderForceFromRadius() "
	      << "No indices used." << std::endl;
    exit(0);
  }
  //Set the variable values
  //////////////////////////////////////////////////////////////////////
  setId("SphereCylinderForceFromRadius");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //////////////////////////////////////////////////////////////////////
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "F_out";
  tmp[1] = "F_in";
  tmp[2] = "R";
  setParameterId( tmp );
}

//!Derivative contribution for force towards or from origo
/*! 
 */
void SphereCylinderForceFromRadius::
derivs(Tissue &T,
       std::vector< std::vector<double> > &cellData,
       std::vector< std::vector<double> > &wallData,
       std::vector< std::vector<double> > &vertexData,
       std::vector< std::vector<double> > &cellDerivs,
       std::vector< std::vector<double> > &wallDerivs,
       std::vector< std::vector<double> > &vertexDerivs ) {
  
  size_t dimension = vertexData[0].size();
  size_t lastPosIndex = dimension-1;
  double F_out = parameter(0);
  double F_in = parameter(1);
  double R = parameter(2);
  //For each vertex 
  for( size_t i=0 ; i<T.numVertex() ; ++i ) {
    if( vertexData[i][lastPosIndex]>0.0 ) {
      //On sphere
      double norm = 0.0;
      for( size_t d=0 ; d<dimension ; d++ )
				norm += vertexData[i][d]*vertexData[i][d];
      norm = std::sqrt(norm);
      if( norm<R ) {
				double coeff = F_out*(R-norm)/norm;
				for( size_t d=0 ; d<dimension ; d++ )
					vertexDerivs[i][d] += coeff*vertexData[i][d];
      }
      else {
				double coeff = F_in*(R-norm)/norm;
				for( size_t d=0 ; d<dimension ; d++ )
					vertexDerivs[i][d] += coeff*vertexData[i][d];
      }
    }
    else {
      //On cylinder
      double norm = 0.0;
      for( size_t d=0 ; d<dimension-1 ; d++ )
				norm += vertexData[i][d]*vertexData[i][d];
      norm = std::sqrt(norm);
      if( norm<R ) {
				double coeff = F_out*(R-norm)/norm;
				for( size_t d=0 ; d<dimension-1 ; d++ )
					vertexDerivs[i][d] += coeff*vertexData[i][d];
      }
      else {
				double coeff = F_in*(R-norm)/norm;
				for( size_t d=0 ; d<dimension-1 ; d++ )
					vertexDerivs[i][d] += coeff*vertexData[i][d];
      }      
    }
  }
}

//!Constructor
VertexNoUpdateFromPosition::
VertexNoUpdateFromPosition(std::vector<double> &paraValue, 
			   std::vector< std::vector<size_t> > 
			   &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //////////////////////////////////////////////////////////////////////
  if( paraValue.size()!=2 ) {
    std::cerr << "VertexNoUpdateFromPosition::"
	      << "VertexNoUpdateFromPosition() "
	      << "Uses two parameters threshold and direction "
	      << "(-1 -> less than)" 
	      << std::endl;
    exit(0);
  }
  if( paraValue[1] != 1.0 && paraValue[1] != -1.0 ) {
    std::cerr << "VertexNoUpdateFromPosition::"
	      << "VertexNoUpdateFromPosition() "
	      << "direction (second parameter) need to be 1 (greater than) "
	      << "or -1 (less than)." << std::endl;
    exit(0);
  }
  if( indValue.size() != 1 || indValue[0].size() != 1 ) {
    std::cerr << "VertexNoUpdateFromPosition::"
	      << "VertexNoUpdateFromPosition() "
	      << "Pos index in first level." << std::endl;
    exit(0);
  }
  //Set the variable values
  //////////////////////////////////////////////////////////////////////
  setId("VertexNoUpdateFromPosition");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //////////////////////////////////////////////////////////////////////
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "threshold";
  tmp[1] = "direction";
  setParameterId( tmp );
}

//!Derivative contribution for force towards or from origo
/*! 
 */
void VertexNoUpdateFromPosition::
derivs(Tissue &T,
       std::vector< std::vector<double> > &cellData,
       std::vector< std::vector<double> > &wallData,
       std::vector< std::vector<double> > &vertexData,
       std::vector< std::vector<double> > &cellDerivs,
       std::vector< std::vector<double> > &wallDerivs,
       std::vector< std::vector<double> > &vertexDerivs ) {
  
  //Check the cancelation for every vertex
  size_t numVertices = T.numVertex();
  size_t posIndex = variableIndex(0,0);
  size_t dimension = vertexData[0].size();
  assert( posIndex < vertexData[0].size() );
  
  for( size_t i=0 ; i<numVertices ; ++i )
    if( (parameter(1)>0 && vertexData[i][posIndex]>parameter(0)) || 
	(parameter(1)<0 && vertexData[i][posIndex]<parameter(0)) )
      for( size_t d=0 ; d<dimension ; ++d )
	vertexDerivs[i][d] = 0.0;
}

//!Constructor
InfiniteWallForce::
InfiniteWallForce(std::vector<double> &paraValue, 
		  std::vector< std::vector<size_t> > 
		  &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //////////////////////////////////////////////////////////////////////
  if( paraValue.size()!=3 ) {
    std::cerr << "InfiniteWallForce::"
	      << "InfiniteWallForce() "
	      << "Uses three parameters k_spring threshold and direction "
	      << "(-1 -> less than)" 
	      << std::endl;
    exit(0);
  }
  if( paraValue[2] != 1.0 && paraValue[2] != -1.0 ) {
    std::cerr << "InfiniteWallForce::"
	      << "InfiniteWallForce() "
	      << "direction (second parameter) need to be 1 (greater than) "
	      << "or -1 (less than)." << std::endl;
    exit(0);
  }
  if( indValue.size() != 1 || indValue[0].size() != 1 ) {
    std::cerr << "InfiniteWallForce::"
	      << "InfiniteWallForce() "
	      << "Pos index in first level." << std::endl;
    exit(0);
  }
  //Set the variable values
  //////////////////////////////////////////////////////////////////////
  setId("InfiniteWallForce");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //////////////////////////////////////////////////////////////////////
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "k_spring";
  tmp[1] = "threshold";
  tmp[2] = "direction";
  setParameterId( tmp );
}

//!Derivative contribution for force towards or from origo
/*! A spring force in a perpendicular direction is applied. Note, the
  wall can only be defined along coordinate axes.
*/
void InfiniteWallForce::
derivs(Tissue &T,
       std::vector< std::vector<double> > &cellData,
       std::vector< std::vector<double> > &wallData,
       std::vector< std::vector<double> > &vertexData,
       std::vector< std::vector<double> > &cellDerivs,
       std::vector< std::vector<double> > &wallDerivs,
       std::vector< std::vector<double> > &vertexDerivs ) {
  
  //Check the cancelation for every vertex
  size_t numVertices = T.numVertex();
  size_t posIndex = variableIndex(0,0);
  //size_t dimension = vertexData[0].size();
  assert( posIndex < vertexData[0].size() );
  
  for( size_t i=0 ; i<numVertices ; ++i ) {
    if( parameter(2)>0 && vertexData[i][posIndex]>parameter(1) ) {
      vertexDerivs[i][posIndex] -= parameter(0)*
	(vertexData[i][posIndex]-parameter(1));
    }
    else if( parameter(2)<0 && vertexData[i][posIndex]<parameter(1) ) {
	    vertexDerivs[i][posIndex] -= parameter(0)*
		    (vertexData[i][posIndex]-parameter(1));
    }
  }
}

//!Constructor
EpidermalVertexForce::
EpidermalVertexForce(std::vector<double> &paraValue, 
		  std::vector< std::vector<size_t> > 
		  &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //////////////////////////////////////////////////////////////////////
  if( paraValue.size()!=2 ) {
    std::cerr << "EpidermalVertexForce::"
	      << "EpidermalVertexForce() "
	      << "Uses two parameters k_strength direction "
	      << "(-1 -> down)" 
	      << std::endl;
    exit(0);
  }
  if( paraValue[1] != 1.0 && paraValue[1] != -1.0 ) {
    std::cerr << "EpidermalVertexForce::"
	      << "EpidermalVertexForce() "
	      << "direction (second parameter) need to be 1 (pos) "
	      << "or -1 (neg direction)." << std::endl;
    exit(0);
  }
  if( indValue.size() != 1 || indValue[0].size() != 1 ) {
    std::cerr << "EpidermalVertexForce::"
	      << "EpidermalVertexForce() "
	      << "Pos index in first level." << std::endl;
    exit(0);
  }
  //Set the variable values
  //////////////////////////////////////////////////////////////////////
  setId("EpidermalVertexForce");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //////////////////////////////////////////////////////////////////////
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "k_spring";
  tmp[1] = "direction";
  setParameterId( tmp );
}

//!Derivative contribution for force towards or from origo
/*! A spring force in a perpendicular direction is applied. Note, the
  wall can only be defined along coordinate axes.
*/
void EpidermalVertexForce::
derivs(Tissue &T,
       std::vector< std::vector<double> > &cellData,
       std::vector< std::vector<double> > &wallData,
       std::vector< std::vector<double> > &vertexData,
       std::vector< std::vector<double> > &cellDerivs,
       std::vector< std::vector<double> > &wallDerivs,
       std::vector< std::vector<double> > &vertexDerivs ) {
  
  //Check the cancelation for every vertex
  size_t numVertices = T.numVertex();
  size_t posIndex = variableIndex(0,0);
  //size_t dimension = vertexData[0].size();
  assert( posIndex < vertexData[0].size() );
  
  for( size_t i=0 ; i<numVertices ; ++i ) {
		int epidermisFlag=0;
		for( size_t k=0 ; k<T.vertex(i).numWall() ; ++k ) {
			if( T.vertex(i).wall(k)->cell1()==T.background() ||
					T.vertex(i).wall(k)->cell1()==T.background() ) {
				epidermisFlag++;
			}
		}		
    if( epidermisFlag ) {
      vertexDerivs[i][posIndex] += parameter(0)*parameter(1);
    }
  }
}

VertexFromPressureExperimental::VertexFromPressureExperimental(std::vector<double> &paraValue, 
												   std::vector< std::vector<size_t> > &indValue)
{
	if (paraValue.size() != 2) {
		std::cerr << "VertexFromPressureExperimental::VertexFromPressureExperimental() "
							<< "Uses two parameter: k no_contraction_flag" << std::endl;
		exit(EXIT_FAILURE);
	}
	
	if (indValue.size() != 1 || indValue[0].size() != 1) {
		std::cerr << "VertexFromPressureExperimental::VertexFromPressureExperimental() "
							<< "Water volume index given.\n";
		exit(EXIT_FAILURE);
	}
	
	setId("VertexFromPressureExperimental");
	setParameter(paraValue);  
	setVariableIndex(indValue);
	
	std::vector<std::string> tmp(numParameter());
	tmp[0] = "k";
	tmp[1] = "no_contraction_flag";
	setParameterId(tmp);
}

void VertexFromPressureExperimental::
derivs(Tissue &T,
			 std::vector< std::vector<double> > &cellData,
			 std::vector< std::vector<double> > &wallData,
			 std::vector< std::vector<double> > &vertexData,
			 std::vector< std::vector<double> > &cellDerivs,
			 std::vector< std::vector<double> > &wallDerivs,
			 std::vector< std::vector<double> > &vertexDerivs)
{
	// ATTENTION! This function assume two dimensions (x, y) and that
	// the vertices are sorted in a clock-wise order.
	
	assert(T.vertex(0).numPosition());
	
	double area;
	for (size_t n = 0; n < T.numCell(); ++n) {
		Cell cell = T.cell(n);
		std::vector< std::pair<double, double> > vertices;
		for (size_t i = 0; i < cell.numVertex(); ++i) {
			std::pair<double, double> vertex;
			vertex.first = vertexData[cell.vertex(i)->index()][0];
			vertex.second = vertexData[cell.vertex(i)->index()][1];
			vertices.push_back(vertex);
		}
		area = polygonArea(vertices);
		double sa = area < 0 ? -1 : +1;
		area = std::fabs(area);
		
		if (parameter(1) != 0.0 && cellData[n][variableIndex(0, 0)] - area < 0) {
			//Std::cerr << "Lower water volume than volume: " << n << " " 
			//				<< area << " " << cellData[n][variableIndex(0,0)]
			//				<< std::endl;
			continue;
		}
		for (size_t i = 1; i < (cell.numVertex() + 1); ++i) {
			Vertex *vertex = cell.vertex(i % cell.numVertex()); // Current vertex in polygon.
			Vertex *pvertex = cell.vertex((i - 1) % cell.numVertex()); // Previous vertex in polygon.
			Vertex *nvertex = cell.vertex((i + 1) % cell.numVertex()); // Next vertex in polygon.
			
			double px = vertexData[pvertex->index()][0];
			double py = vertexData[pvertex->index()][1];
			double nx = vertexData[nvertex->index()][0];
			double ny = vertexData[nvertex->index()][1];
			
			double dAdx = sa * 0.5 * (-py + ny);
			double dAdy = sa * 0.5 * (px - nx);
			
			vertexDerivs[vertex->index()][0] += parameter(0) * (1 - area / cellData[n][variableIndex(0, 0)]) * dAdx;
			vertexDerivs[vertex->index()][1] += parameter(0) * (1 - area / cellData[n][variableIndex(0, 0)]) * dAdy;
		}
	}
}

double VertexFromPressureExperimental::polygonArea(std::vector< std::pair<double, double> > vertices)
{
	double area = 0.0;
	size_t N = vertices.size();
	for (size_t n = 0; n < N; ++n) {
		area += vertices[n].first * vertices[(n + 1) % N].second;
		area -= vertices[(n + 1) % N].first * vertices[n].second;
	}
	area *= 0.5;
	return area;
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
			 std::vector< std::vector<double> > &cellData,
			 std::vector< std::vector<double> > &wallData,
			 std::vector< std::vector<double> > &vertexData,
			 std::vector< std::vector<double> > &cellDerivs,
			 std::vector< std::vector<double> > &wallDerivs,
			 std::vector< std::vector<double> > &vertexDerivs)
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

CellVolumeExperimental::
CellVolumeExperimental(std::vector<double> &paraValue,
											 std::vector< std::vector<size_t> > &indValue)
{
	if (paraValue.size() != 3) {
		std::cerr << "CellVolumeExperimental::CellVolumeExperimental() "
							<< "Uses three parameters: k_p, P_max and k_pp" << std::endl;
		exit(EXIT_FAILURE);
	}
	
	if (indValue.size() < 2 || indValue.size() > 3
			|| indValue[0].size() != 2 
			|| (indValue.size()==3 && indValue[2].size() != 1 ) ) {
		std::cerr << "CellVolumeExperimental::CellVolumeExperimental() "
							<< "Wall length index and cell volume index must be given in "
							<< "first level.\n"
							<< "Force indices must ge given in second level. "
							<< "Optionally index for saving the pressure can be"
							<< " given at third level." << std::endl; 		
		exit(EXIT_FAILURE);
	}
	
	setId("VertexFromWallSpringExperimental");
	setParameter(paraValue);
	setVariableIndex(indValue);
	
	std::vector<std::string> tmp(numParameter());
	tmp[0] = "k_p";
	tmp[1] = "P_max";
	tmp[2] = "k_pp";
	
	setParameterId(tmp);
}

void CellVolumeExperimental::
derivs(Tissue &T,
			 std::vector< std::vector<double> > &cellData,
			 std::vector< std::vector<double> > &wallData,
			 std::vector< std::vector<double> > &vertexData,
			 std::vector< std::vector<double> > &cellDerivs,
			 std::vector< std::vector<double> > &wallDerivs,
			 std::vector< std::vector<double> > &vertexDerivs)
{
	for (size_t n = 0; n < T.numCell(); ++n) {
		Cell cell = T.cell(n);
		
		double P = 0.0;
		double sum = 0.0;
		for (size_t i = 0; i < cell.numWall(); ++i) {
			size_t vertex1Index = cell.wall(i)->vertex1()->index();
			size_t vertex2Index = cell.wall(i)->vertex2()->index();
			size_t dimensions = vertexData[vertex1Index].size();
			
			double distance = 0.0;
			for (size_t d = 0; d < dimensions; ++d) {
				distance += (vertexData[vertex1Index][d] - vertexData[vertex2Index][d])
					* (vertexData[vertex1Index][d] - vertexData[vertex2Index][d]);
			}
			distance = std::sqrt(distance);
			
			for (size_t j = 0; j < numVariableIndex(1); ++j)
				P += wallData[cell.wall(i)->index()][variableIndex(1, j)];
			sum += distance;
		}
		P *= parameter(2);

		if (numVariableIndexLevel()==3)
			cellData[n][variableIndex(2,0)]=P;
		
		if( parameter(1)-P>0.0 )
			cellDerivs[cell.index()][variableIndex(0, 1)] += 
				parameter(0) * (parameter(1) - P) * sum;
	}
}


