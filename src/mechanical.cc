//
// Filename     : mechanical.cc
// Description  : Classes describing updates due to mechanical interaction
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : April 2006
// Revision     : $Id:$
//
#include <utility>
#include <vector>
#include "baseReaction.h"
#include "mechanical.h"
#include "tissue.h"

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
  size_t dimension;
  dimension = T.vertex(0).numPosition(); 
	
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
	      << "Uses two parameters K_force direction (-1 -> inwards)" 
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

CellVolumeExperimental::
CellVolumeExperimental(std::vector<double> &paraValue,
											 std::vector< std::vector<size_t> > &indValue)
{
	if (paraValue.size() != 4) {
		std::cerr << "CellVolumeExperimental::CellVolumeExperimental() "
							<< "Uses four parameters: k_p, P_max, k_pp and "
							<< "allowShrink_flag." << std::endl;
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
	tmp[3] = "allowShrink_flag";

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
				P += wallData[cell.wall(i)->index()][variableIndex(1, j)]/distance;
			sum += distance;
		}
		P *= parameter(2);

		if (numVariableIndexLevel()==3)
			cellData[n][variableIndex(2,0)]=P;
		
		if( parameter(3) || parameter(1)-P>0.0 )
			cellDerivs[cell.index()][variableIndex(0, 1)] += 
				parameter(0) * (parameter(1) - P) * sum;
	}
}

EpidermalRadialForce::EpidermalRadialForce(std::vector<double> &paraValue,
								   std::vector< std::vector<size_t> > &indValue)
{
	if (paraValue.size() != 1) {
		std::cerr << "EpidermalRadialForce::EpidermalRadialForce() " 
							<< "Uses one parameter: force" << std::endl;
		exit(EXIT_FAILURE);
	}
	
	if (indValue.size() != 0) {
		std::cerr << "EpidermalRadialForce::EpidermalRadialForce() "
                    << "No indices are given." << std::endl;
		exit(EXIT_FAILURE);
	}
	
	setId("EpidermalRadialForce");
	setParameter(paraValue);
	setVariableIndex(indValue);
	
	std::vector<std::string> tmp(numParameter());
	tmp[0] = "force";

	setParameterId(tmp);
}

void EpidermalRadialForce::derivs(Tissue &T,
						    std::vector< std::vector<double> > &cellData,
						    std::vector< std::vector<double> > &wallData,
						    std::vector< std::vector<double> > &vertexData,
						    std::vector< std::vector<double> > &cellDerivs,
						    std::vector< std::vector<double> > &wallDerivs,
						    std::vector< std::vector<double> > &vertexDerivs)
{
	for (size_t n = 0; n < T.numVertex(); ++n) {
		Vertex vertex = T.vertex(n);
		
		// If vertex is not in the epidermal layer then skip.
		bool isEpidermalVertex = false;
		for (size_t i = 0; i < vertex.numWall(); ++i) {
			Wall *wall = vertex.wall(i);
			if ((wall->cell1()->index() == (size_t) -1) || (wall->cell2()->index() == (size_t) -1)) {
				isEpidermalVertex = true;
				break;
			}
		}
		if (isEpidermalVertex == false) {
			continue;

		}

		double x = vertex.position(0);
		double y = vertex.position(1);
		double A = std::sqrt(x * x + y * y);
		if (A == 0)
			continue;
		x /= A;
		y /= A;

// 		std::cerr << "Vertex " << vertex.index() << std::endl;
// 		std::cerr << " x = " << x << std::endl;
// 		std::cerr << " y = " << y << std::endl;

		vertexDerivs[vertex.index()][0] += -parameter(0) * x;
		vertexDerivs[vertex.index()][1] += -parameter(0) * y;
	}
}

PerpendicularWallPressure::PerpendicularWallPressure(std::vector<double> &paraValue,
										   std::vector< std::vector<size_t> > &indValue)
{
	if (paraValue.size() != 1) {
		std::cerr << "PerpendicularWallPressure::PerpendicularWallPressure() " 
				<< "Uses one parameter: k_force" << std::endl;
		exit(EXIT_FAILURE);
	}
	
	if (indValue.size() != 1 || indValue[0].size() != 1) {
		std::cerr << "PerpendicularWallPressure::PerpendicularWallPressure() "
                    << "First level gives pressure index in cell." << std::endl;
		exit(EXIT_FAILURE);
	}
	
	setId("PerpendicularWallPressure");
	setParameter(paraValue);
	setVariableIndex(indValue);
	
	std::vector<std::string> tmp(numParameter());
	tmp[0] = "k_force";

	setParameterId(tmp);
}

void PerpendicularWallPressure::derivs(Tissue &T,
							    std::vector< std::vector<double> > &cellData,
							    std::vector< std::vector<double> > &wallData,
							    std::vector< std::vector<double> > &vertexData,
							    std::vector< std::vector<double> > &cellDerivs,
							    std::vector< std::vector<double> > &wallDerivs,
							    std::vector< std::vector<double> > &vertexDerivs)
{
	size_t dimension = vertexData[0].size();
	if (dimension==2) {
		for (size_t n = 0; n < T.numWall(); ++n) {
			Wall wall = T.wall(n);
			std::vector<Cell *> cells(2);
			cells[0] = wall.cell1();
			cells[1] = wall.cell2();
			std::vector<Vertex *> vertices(2);
			vertices[0] = wall.vertex1();
			vertices[1] = wall.vertex2();
			
			double nx = vertexData[vertices[1]->index()][0] - vertexData[vertices[0]->index()][0];
			double ny = vertexData[vertices[1]->index()][1] - vertexData[vertices[0]->index()][1];
			double A = std::sqrt(nx * nx + ny * ny);
			nx = nx / A;
			ny = ny / A;
			
			for (size_t i = 0; i < cells.size(); ++i) {
				Cell *cell = cells[i];
				
				if (cell == T.background()) {
					continue;
				}
				
				double force = parameter(0) * 
					cellData[cell->index()][variableIndex(0, 0)] *  
					wall.lengthFromVertexPosition(vertexData);
				
				double x = 0.0;
				double y = 0.0;
				for (size_t j = 0; j < cell->numVertex(); ++j) {
					x += vertexData[cell->vertex(j)->index()][0];
					y += vertexData[cell->vertex(j)->index()][1];
				}
				x = (x / cell->numVertex()) - vertexData[vertices[0]->index()][0];
				y = (y / cell->numVertex()) - vertexData[vertices[0]->index()][1];
				
				int sign = (nx * y - ny * x) >= 0 ? 1 : -1;
				
				for (size_t k = 0; k < vertices.size(); ++k) {
					Vertex *vertex = vertices[k];
					vertexDerivs[vertex->index()][0] += ny * sign * force;
					vertexDerivs[vertex->index()][1] += - nx * sign * force;
				}
			}
		}
	}
	else if (dimension==3) {

		size_t numCell = T.numCell();
		for (size_t i=0; i<numCell; ++i) {
			Cell* cell1 = &(T.cell(i));
			size_t c1I = cell1->index();
			std::vector<double> n1 = cell1->getNormalToPCAPlane();
			int n1Sign = cell1->vectorSignFromSort(n1,vertexData);
			std::vector<double> n = n1;
			std::vector<double> x1 = cell1->positionFromVertex(vertexData);

			size_t numNeigh = cell1->numWall();
			for (size_t k=0; k<numNeigh; ++k) {
				Cell* cell2 = cell1->cellNeighbor(k);
				size_t c2I = cell2->index();
				// Get normal to second cell (if not background) and calculate combined normal (n)
				if (cell2!=T.background() && c2I<c1I) {
					continue;
				}
				else if (cell2!=T.background()) {
					std::vector<double> n2(dimension);
					n2 = cell2->getNormalToPCAPlane();
					int n2Sign = cell2->vectorSignFromSort(n2,vertexData);
					for (size_t d=0; d<dimension; ++d) 
						n[d] = n1Sign*n1[d]+n2Sign*n2[d];
				}
				// Get wall direction
				size_t v1I = cell1->wall(k)->vertex1()->index();
				size_t v2I = cell1->wall(k)->vertex2()->index();
				std::vector<double> nw(dimension),xw(dimension);
				for (size_t d=0; d<dimension; ++d) {
					nw[d] = vertexData[v2I][d]-vertexData[v1I][d];
					xw[d] = 0.5*(vertexData[v2I][d]+vertexData[v1I][d]);					
				}
				// Normalize n and nw
				double norm=0.0,normW=0.0;
				for (size_t d=0; d<dimension; ++d) { 
					norm += n[d]*n[d];
					normW += nw[d]*nw[d];
				}
				assert(norm>0.0 && normW>0.0);
				norm = std::sqrt(norm);
				normW = std::sqrt(normW);
				double normFac=1.0/norm,normWFac=1.0/normW;
				for (size_t d=0; d<dimension; ++d) {
					n[d] *= normFac;
					nw[d] *= normWFac;
				}
				
				// Calculate force direction perpendicular to n and nw, plus sign
				std::vector<double> nF(dimension);
				nF[0] = n[1]*nw[2]-n[2]*nw[1];
				nF[1] = n[2]*nw[0]-n[0]*nw[2];
				nF[2] = n[0]*nw[1]-n[1]*nw[0];
				norm=0.0;
				for (size_t d=0; d<dimension; ++d) { 
					norm += nF[d]*nF[d];
				}
				assert(norm>0.0);
				norm = std::sqrt(norm);
				normFac=1.0/norm;
				for (size_t d=0; d<dimension; ++d) {
					nF[d] *= normFac;
				}

				// sign from scalar product?
				double scalar=0.0;
				for (size_t d=0; d<dimension; ++d)
					scalar += (xw[d]-x1[d])*nF[d];
				int sign=1;
				if (scalar<0.0)
					sign=-1;
				double forceFactor = parameter(0)*cellData[c1I][variableIndex(0,0)]*normW;
				if (cell2!=T.background())
					forceFactor -= parameter(0)*cellData[c2I][variableIndex(0,0)]*normW;
				
				// Print result (wall and nF) for debugging
				//double fac=0.5*sign;
				//std::cerr << c1I << " " << vertexData[v1I][0] << " " 
				//				<< vertexData[v1I][1] << " " 
				//				<< vertexData[v1I][2] << " " 
				//				<< xw[0] << " " << xw[1] << " " << xw[2] << std::endl;
				//std::cerr << c1I << " " << vertexData[v2I][0] << " " 
				//				<< vertexData[v2I][1] << " " 
				//				<< vertexData[v2I][2] << " " 
				//				<< xw[0]+fac*nF[0] << " " << xw[1]+fac*nF[1] << " " << xw[2]+fac*nF[2] 
				//				<< std::endl << std::endl << std::endl;
				
				for (size_t d=0; d<dimension; ++d) {
					vertexDerivs[v1I][d] += forceFactor*sign*nF[d];
					vertexDerivs[v2I][d] += forceFactor*sign*nF[d];
				}				
			}
		}
	}
	else {
		std::cerr << "PerpendicularWallPressure::derivs() Only applicable for two or three"
							<< " dimensions." << std::endl;
		exit(-1);
	}
}

ContinousMTDirection::ContinousMTDirection(std::vector<double> &paraValue,
																					 std::vector< std::vector<size_t> > &indValue)
{
	if (paraValue.size() != 1) {
		std::cerr << "ContinousMTDirection::ContinousMTDirection() " 
							<< "Uses one parameter: k_rate" << std::endl;
		exit(EXIT_FAILURE);
	}
	
	if (indValue.size() != 2 || indValue[0].size() != 1 ||  indValue[1].size() != 1) {
		std::cerr << "ContinousMTDirection::ContinousMTDirection() " << std::endl
							<< "First level gives target direction index." << std::endl
			  << "Second level gives real direction index." << std::endl;
		exit(EXIT_FAILURE);
	}
	
	setId("ContinousMTDirection");
	setParameter(paraValue);
	setVariableIndex(indValue);
	
	std::vector<std::string> tmp(numParameter());
	tmp[0] = "k_rate";
	
	setParameterId(tmp);
}

void ContinousMTDirection::derivs(Tissue &T,
				  std::vector< std::vector<double> > &cellData,
				  std::vector< std::vector<double> > &wallData,
				  std::vector< std::vector<double> > &vertexData,
				  std::vector< std::vector<double> > &cellDerivs,
				  std::vector< std::vector<double> > &wallDerivs,
				  std::vector< std::vector<double> > &vertexDerivs)
{
	size_t target = variableIndex(0, 0);
	size_t real = variableIndex(1, 0);

	double k_rate = parameter(0);

	for (size_t n = 0; n < T.numCell(); ++n) {
		Cell cell = T.cell(n);
		size_t index = cell.index();

		double x = cellData[index][real + 0];
		double y = cellData[index][real + 1];

		double sigma = std::atan2(y, x);

		while (sigma > 0.5 * M_PI || sigma <= -0.5 * M_PI) {
			if (sigma > 0.5 * M_PI) {
				sigma -= M_PI;
			}
			if (sigma <= -0.5 * M_PI) {
				sigma += M_PI;
			}
		}

		double dx = cellData[index][target + 0];
		double dy = cellData[index][target + 1];

		double dsigma = std::atan2(dy, dx);

		while (dsigma > 0.5 * M_PI || dsigma <= -0.5 * M_PI) {
			if (dsigma > 0.5 * M_PI) {
				dsigma -= M_PI;
			}
			if (dsigma <= -0.5 * M_PI) {
				dsigma += M_PI;
			}
		}

		double angle = dsigma - sigma;

		while (angle > M_PI || angle <= -M_PI) {
			if (angle > M_PI) {
				angle -= 2.0 * M_PI;
			}
			if (angle <= -M_PI) {
				angle += 2.0 * M_PI;
			}
		}

		double speed = k_rate * std::abs(angle) / (0.25 * M_PI + std::abs(angle));

		speed *= (angle >= 0) ? +1 : -1;

		cellDerivs[index][real + 0] += -y * speed;
		cellDerivs[index][real + 1] += x * speed;
	}
}

VertexFromCellPlane::
VertexFromCellPlane(std::vector<double> &paraValue,
																	std::vector< std::vector<size_t> > &indValue)
{
	if (paraValue.size() != 1) {
		std::cerr << "VertexFromCellPlane::VertexFromCellPlane() " 
							<< "Uses one parameter: k_force" << std::endl;
		exit(EXIT_FAILURE);
	}
	
	if (indValue.size() != 0) {
		std::cerr << "VertexFromCellPlane::VertexFromCellPlane() " 
							<< std::endl
							<< "No variable index used." << std::endl;
		exit(EXIT_FAILURE);
	}
	
	setId("VertexFromCellPlane");
	setParameter(paraValue);
	setVariableIndex(indValue);
	
	std::vector<std::string> tmp(numParameter());
	tmp[0] = "k_force";
	
	setParameterId(tmp);
}

void VertexFromCellPlane::
derivs(Tissue &T,
			 std::vector< std::vector<double> > &cellData,
			 std::vector< std::vector<double> > &wallData,
			 std::vector< std::vector<double> > &vertexData,
			 std::vector< std::vector<double> > &cellDerivs,
			 std::vector< std::vector<double> > &wallDerivs,
			 std::vector< std::vector<double> > &vertexDerivs)
{
	size_t dimension = vertexData[0].size();
	if (dimension!=3) {
		std::cerr << "VertexFromCellPlane::VertexFromCellPlane() " 
							<< "Only implemented for three dimensions." << std::endl;
		exit(EXIT_FAILURE);
	}
	unsigned int numFlipNormal=0;
	for (size_t n = 0; n < T.numCell(); ++n) {
		Cell cell = T.cell(n);
		// This calculation should now be done in reaction CalculatePCAPlane
		//cell.calculatePCAPlane(vertexData);
		unsigned int flipFlag=0;
		
		std::vector<double> normal = cell.getNormalToPCAPlane();
		double norm=0.0;
		for (size_t d=0; d<dimension; ++d)
			norm += normal[d]*normal[d];
		if (norm != 1.0) {
			norm = std::sqrt(norm);
			assert(norm>0.0);
			double normFac = 1.0/norm; 
			for (size_t d=0; d<dimension; ++d)
				normal[d] *= normFac;
		}
		
		std::vector<int> scalarProdSign(cell.numVertex());
		std::vector<double> scalarProdVal(cell.numVertex());
		double scalarProdSum=0.0;
		for (size_t k=0; k<cell.numVertex(); ++k) {
			size_t k2=(k+1)%cell.numVertex();
			size_t k3=(k+2)%cell.numVertex();
			//Make sure ((v2-v1)x(v3-v2))n has same sign for all cells
			assert(cell.numVertex()>2);
			std::vector<double> nw1(dimension),nw2(dimension);
			for (size_t d=0; d<dimension; ++d) {
				nw1[d] = vertexData[cell.vertex(k2)->index()][d]-vertexData[cell.vertex(k)->index()][d];
				nw2[d] = vertexData[cell.vertex(k3)->index()][d]-vertexData[cell.vertex(k2)->index()][d];
			}
			//cross product
			double scalarProd=0.0;
			for (size_t d1=0; d1<dimension; ++d1) {
				size_t d2=(d1+1)%dimension;
				size_t d3=(d1+2)%dimension;
				scalarProd += (nw1[d1]*nw2[d2]-nw1[d2]*nw2[d1])*normal[d3];
			}
			scalarProdVal[k] = scalarProd;
			scalarProdSum += scalarProd;
			if (scalarProd>0.0)
				scalarProdSign[k]=1;
			else
				scalarProdSign[k]=-1;
		}
//  		for (size_t k=1; k<cell.numVertex(); ++k)
//  			if (scalarProdSign[k]!=scalarProdSign[0]) {
//  				std::cerr << "Cell " << n << " has diverging signs on scalar product." << std::endl;
//  				break;
//  			}
		int scalarProdSignSum=0;
		for (size_t k=0; k<scalarProdSign.size(); ++k)
			scalarProdSignSum += scalarProdSign[k];
		
		if (scalarProdSignSum<0) {
			numFlipNormal++;
			flipFlag=1;
			for (size_t d=0; d<dimension; ++d)
				normal[d] = -normal[d];
		}
		else if (scalarProdSignSum==0) {
			//std::cerr << "Cell " << n << " has no majority sign in right hand rule expression." 
			//				<< std::endl;
			if (std::fabs(scalarProdSum)>0.01) {
				if (scalarProdSum<0.0) {
					numFlipNormal++;
					flipFlag=1;
					for (size_t d=0; d<dimension; ++d)
						normal[d] = -normal[d];
				}
			}
			else {
				std::vector<double> center = cell.positionFromVertex(vertexData);
				// Print all walls
				for (size_t k=0; k<cell.numWall(); ++k) {
					std::cerr << "0 "; 
					for (size_t d=0; d<dimension; ++d)
						std::cerr << vertexData[cell.wall(k)->vertex1()->index()][d] << " ";
					std::cerr << std::endl << "0 "; 
					for (size_t d=0; d<dimension; ++d)
						std::cerr << vertexData[cell.wall(k)->vertex2()->index()][d] << " ";
					std::cerr << std::endl << std::endl;
				}		 
				std::cerr << "1 "; 
				for (size_t d=0; d<dimension; ++d)
					std::cerr << vertexData[cell.wall(0)->vertex1()->index()][d] << " ";
				std::cerr << std::endl << std::endl;
				std::cerr << "2 "; 
				for (size_t d=0; d<dimension; ++d)
					std::cerr << center[d] << " ";
				std::cerr << std::endl << "2 ";
				for (size_t d=0; d<dimension; ++d)
					std::cerr << center[d]+normal[d] << " ";
				std::cerr << std::endl << std::endl;
				std::cerr << "3 "; 
				for (size_t d=0; d<dimension; ++d)
					std::cerr << normal[d] << " ";
				std::cerr << std::endl << std::endl;
				std::cerr << "4 "; 
				for (size_t i=0; i<scalarProdVal.size(); ++i)
					std::cerr << scalarProdVal[i] << " ";
				std::cerr << std::endl << std::endl;
				std::cerr << "5 "; 
				for (size_t i=0; i<scalarProdSign.size(); ++i)
					std::cerr << scalarProdSign[i] << " ";
				std::cerr << std::endl << std::endl;
				
				exit(-1);
			}
		}
		// Get the cell size
		//double A = cell.calculateVolume(vertexData);
		double A=1.0;
		//update the vertex derivatives
		for (size_t d=0; d<dimension; ++d) {
			//escellData[n][d]=normal[d];
			for (size_t k=0; k<cell.numVertex(); ++k)
				vertexDerivs[cell.vertex(k)->index()][d] += parameter(0) * A * normal[d];
		}
	}	
	//std::cerr << numFlipNormal << " cells out of " << T.numCell() << " has flipped normal."
	//				<< std::endl;
}

VertexFromCellPlaneNormalized::
VertexFromCellPlaneNormalized(std::vector<double> &paraValue,
															std::vector< std::vector<size_t> > &indValue)
{
	if (paraValue.size() != 1) {
		std::cerr << "VertexFromCellPlaneNormalized::VertexFromCellPlaneNormalized() " 
							<< "Uses one parameter: k_force" << std::endl;
		exit(EXIT_FAILURE);
	}
	
	if (indValue.size() != 0) {
		std::cerr << "VertexFromCellPlaneNormalized::VertexFromCellPlaneNormalized() " 
							<< std::endl
							<< "No variable index used." << std::endl;
		exit(EXIT_FAILURE);
	}
	
	setId("VertexFromCellPlaneNormalized");
	setParameter(paraValue);
	setVariableIndex(indValue);
	
	std::vector<std::string> tmp(numParameter());
	tmp[0] = "k_force";
	
	setParameterId(tmp);
}

void VertexFromCellPlaneNormalized::
derivs(Tissue &T,
			 std::vector< std::vector<double> > &cellData,
			 std::vector< std::vector<double> > &wallData,
			 std::vector< std::vector<double> > &vertexData,
			 std::vector< std::vector<double> > &cellDerivs,
			 std::vector< std::vector<double> > &wallDerivs,
			 std::vector< std::vector<double> > &vertexDerivs)
{
	size_t dimension = vertexData[0].size();
	if (dimension!=3) {
		std::cerr << "VertexFromCellPlaneNormalized::VertexFromCellPlaneNormalized() " 
							<< "Only implemented for three dimensions." << std::endl;
		exit(EXIT_FAILURE);
	}
	std::vector<double> tmpD(dimension);
	std::vector< std::vector<double> > vertexDerivsTmp(vertexDerivs.size(),tmpD);
	unsigned int numFlipNormal=0;
	for (size_t n = 0; n < T.numCell(); ++n) {
		Cell cell = T.cell(n);
		// This calculation should now be done in reaction CalculatePCAPlane
		//cell.calculatePCAPlane(vertexData);
		unsigned int flipFlag=0;
		
		std::vector<double> normal = cell.getNormalToPCAPlane();
		double norm=0.0;
		for (size_t d=0; d<dimension; ++d)
			norm += normal[d]*normal[d];
		if (norm != 1.0) {
			norm = std::sqrt(norm);
			assert(norm>0.0);
			double normFac = 1.0/norm; 
			for (size_t d=0; d<dimension; ++d)
				normal[d] *= normFac;
		}
		
		std::vector<int> scalarProdSign(cell.numVertex());
		std::vector<double> scalarProdVal(cell.numVertex());
		double scalarProdSum=0.0;
		for (size_t k=0; k<cell.numVertex(); ++k) {
			size_t k2=(k+1)%cell.numVertex();
			size_t k3=(k+2)%cell.numVertex();
			//Make sure ((v2-v1)x(v3-v2))n has same sign for all cells
			assert(cell.numVertex()>2);
			std::vector<double> nw1(dimension),nw2(dimension);
			for (size_t d=0; d<dimension; ++d) {
				nw1[d] = vertexData[cell.vertex(k2)->index()][d]-vertexData[cell.vertex(k)->index()][d];
				nw2[d] = vertexData[cell.vertex(k3)->index()][d]-vertexData[cell.vertex(k2)->index()][d];
			}
			//cross product
			double scalarProd=0.0;
			for (size_t d1=0; d1<dimension; ++d1) {
				size_t d2=(d1+1)%dimension;
				size_t d3=(d1+2)%dimension;
				scalarProd += (nw1[d1]*nw2[d2]-nw1[d2]*nw2[d1])*normal[d3];
			}
			scalarProdVal[k] = scalarProd;
			scalarProdSum += scalarProd;
			if (scalarProd>0.0)
				scalarProdSign[k]=1;
			else
				scalarProdSign[k]=-1;
		}
//  		for (size_t k=1; k<cell.numVertex(); ++k)
//  			if (scalarProdSign[k]!=scalarProdSign[0]) {
//  				std::cerr << "Cell " << n << " has diverging signs on scalar product." << std::endl;
//  				break;
//  			}
		int scalarProdSignSum=0;
		for (size_t k=0; k<scalarProdSign.size(); ++k)
			scalarProdSignSum += scalarProdSign[k];
		
		if (scalarProdSignSum<0) {
			numFlipNormal++;
			flipFlag=1;
			for (size_t d=0; d<dimension; ++d)
				normal[d] = -normal[d];
		}
		else if (scalarProdSignSum==0) {
			//std::cerr << "Cell " << n << " has no majority sign in right hand rule expression." 
			//				<< std::endl;
			if (std::fabs(scalarProdSum)>0.01) {
				if (scalarProdSum<0.0) {
					numFlipNormal++;
					flipFlag=1;
					for (size_t d=0; d<dimension; ++d)
						normal[d] = -normal[d];
				}
			}
			else {
				std::vector<double> center = cell.positionFromVertex(vertexData);
				// Print all walls
				for (size_t k=0; k<cell.numWall(); ++k) {
					std::cerr << "0 "; 
					for (size_t d=0; d<dimension; ++d)
						std::cerr << vertexData[cell.wall(k)->vertex1()->index()][d] << " ";
					std::cerr << std::endl << "0 "; 
					for (size_t d=0; d<dimension; ++d)
						std::cerr << vertexData[cell.wall(k)->vertex2()->index()][d] << " ";
					std::cerr << std::endl << std::endl;
				}		 
				std::cerr << "1 "; 
				for (size_t d=0; d<dimension; ++d)
					std::cerr << vertexData[cell.wall(0)->vertex1()->index()][d] << " ";
				std::cerr << std::endl << std::endl;
				std::cerr << "2 "; 
				for (size_t d=0; d<dimension; ++d)
					std::cerr << center[d] << " ";
				std::cerr << std::endl << "2 ";
				for (size_t d=0; d<dimension; ++d)
					std::cerr << center[d]+normal[d] << " ";
				std::cerr << std::endl << std::endl;
				std::cerr << "3 "; 
				for (size_t d=0; d<dimension; ++d)
					std::cerr << normal[d] << " ";
				std::cerr << std::endl << std::endl;
				std::cerr << "4 "; 
				for (size_t i=0; i<scalarProdVal.size(); ++i)
					std::cerr << scalarProdVal[i] << " ";
				std::cerr << std::endl << std::endl;
				std::cerr << "5 "; 
				for (size_t i=0; i<scalarProdSign.size(); ++i)
					std::cerr << scalarProdSign[i] << " ";
				std::cerr << std::endl << std::endl;
				
				exit(-1);
			}
		}
		// Get the cell size
		//double A = cell.calculateVolume(vertexData);
		double A=1.0;
		//update the vertex derivatives
		for (size_t d=0; d<dimension; ++d) {
			//escellData[n][d]=normal[d];
			for (size_t k=0; k<cell.numVertex(); ++k)
				vertexDerivsTmp[cell.vertex(k)->index()][d] += A * normal[d];
		}
	}	
	// Normalize and update all vertices
	size_t nVertex = vertexDerivs.size();
	for (size_t i=0; i<nVertex; ++i) {
		double norm=0.0;
		for (size_t d=0; d<dimension; ++d)
			norm += vertexDerivsTmp[i][d]*vertexDerivsTmp[i][d];
		double normFac = 1.0/std::sqrt(norm);
		for (size_t d=0; d<dimension; ++d) {
			vertexDerivsTmp[i][d] *= normFac;
			vertexDerivs[i][d] += parameter(0) * vertexDerivsTmp[i][d];
		}
	}
	//std::cerr << numFlipNormal << " cells out of " << T.numCell() << " has flipped normal."
	//				<< std::endl;
}

VertexFromCellPlaneSphereCylinder::
VertexFromCellPlaneSphereCylinder(std::vector<double> &paraValue,
																	std::vector< std::vector<size_t> > &indValue)
{
	if (paraValue.size() != 1) {
		std::cerr << "VertexFromCellPlaneSphereCylinder::VertexFromCellPlaneSphereCylinder() " 
							<< "Uses one parameter: k_force" << std::endl;
		exit(EXIT_FAILURE);
	}
	
	if (indValue.size() != 0) {
		std::cerr << "VertexFromCellPlaneSphereCylinder::VertexFromCellPlaneSphereCylinder() " 
							<< std::endl
							<< "No variable index used." << std::endl;
		exit(EXIT_FAILURE);
	}
	
	setId("VertexFromCellPlaneSphereCylinder");
	setParameter(paraValue);
	setVariableIndex(indValue);
	
	std::vector<std::string> tmp(numParameter());
	tmp[0] = "k_force";
	
	setParameterId(tmp);
}

void VertexFromCellPlaneSphereCylinder::
derivs(Tissue &T,
			 std::vector< std::vector<double> > &cellData,
			 std::vector< std::vector<double> > &wallData,
			 std::vector< std::vector<double> > &vertexData,
			 std::vector< std::vector<double> > &cellDerivs,
			 std::vector< std::vector<double> > &wallDerivs,
			 std::vector< std::vector<double> > &vertexDerivs)
{
	size_t dimension = vertexData[0].size();
	if (dimension!=3) {
		std::cerr << "VertexFromCellPlaneSphereCylinder::VertexFromCellPlaneSphereCylinder() " 
							<< "Only implemented for three dimensions." << std::endl;
		exit(EXIT_FAILURE);
	}
	for (size_t n = 0; n < T.numCell(); ++n) {
		Cell cell = T.cell(n);
		// This calculation should now be done in CalculatePCAPlane reaction.
		//cell.calculatePCAPlane(vertexData);
		std::vector<double> normal = cell.getNormalToPCAPlane();
		double norm=0.0;
		for (size_t d=0; d<dimension; ++d)
			norm += normal[d]*normal[d];
		if (norm != 1.0) {
			norm = std::sqrt(norm);
			assert(norm>0.0);
			for (size_t d=0; d<dimension; ++d)
				normal[d] /= norm;
		}
		std::vector<double> cellPos = cell.positionFromVertex(vertexData);

		std::vector<double> scNorm(dimension);
		for (size_t d=0; d<dimension; ++d)
			scNorm[d] = cellPos[d];
		if( cellPos[2]<0.0 )
			scNorm[2]=0.0;
		double scalarProd=0.0;
		for (size_t d=0; d<dimension; ++d)
			scalarProd += scNorm[d]*normal[d];
		if (scalarProd<0.0) {
			for (size_t d=0; d<dimension; ++d)
				normal[d] = -normal[d];
		}
		// Introduce cell size
		//double A = cell.calculateVolume(vertexData); 
		double A = 1.0;
		for (size_t k=0; k<cell.numVertex(); ++k)
			for (size_t d=0; d<dimension; ++d)
				vertexDerivs[cell.vertex(k)->index()][d] += parameter(0) * A * normal[d];
	}	
}

VertexFromCellPlaneSphereCylinderConcentrationHill::
VertexFromCellPlaneSphereCylinderConcentrationHill(std::vector<double> &paraValue,
																	std::vector< std::vector<size_t> > &indValue)
{
	if (paraValue.size() != 4) {
		std::cerr << "VertexFromCellPlaneSphereCylinderConcentrationHill::"
							<< "VertexFromCellPlaneSphereCylinderConcentrationHill() " 
							<< "Uses four parameters: k_forceConst, k_forceHill, K_Hill, n_Hill" << std::endl;
		exit(EXIT_FAILURE);
	}
	
	if (indValue.size() != 1 || indValue[0].size() != 1) {
		std::cerr << "VertexFromCellPlaneSphereCylinderConcentrationHill::"
							<< "VertexFromCellPlaneSphereCylinderConcentrationHill() " 
							<< std::endl
							<< "One variable index for concentration index is used." << std::endl;
		exit(EXIT_FAILURE);
	}
	
	setId("VertexFromCellPlaneSphereCylinderConcentrationHill");
	setParameter(paraValue);
	setVariableIndex(indValue);
	
	std::vector<std::string> tmp(numParameter());
	tmp[0] = "k_forceConst";
	tmp[1] = "k_forceHill";
	tmp[2] = "K_Hill";
	tmp[3] = "n_Hill";
	
	setParameterId(tmp);
}

void VertexFromCellPlaneSphereCylinderConcentrationHill::
derivs(Tissue &T,
			 std::vector< std::vector<double> > &cellData,
			 std::vector< std::vector<double> > &wallData,
			 std::vector< std::vector<double> > &vertexData,
			 std::vector< std::vector<double> > &cellDerivs,
			 std::vector< std::vector<double> > &wallDerivs,
			 std::vector< std::vector<double> > &vertexDerivs)
{
	size_t dimension = vertexData[0].size();
	if (dimension!=3) {
		std::cerr << "VertexFromCellPlaneSphereCylinderConcentrationHill::"
							<< "VertexFromCellPlaneSphereCylinderConcentrationHill() " 
							<< "Only implemented for three dimensions." << std::endl;
		exit(EXIT_FAILURE);
	}
	double Kpow = std::pow(parameter(2),parameter(3));
	for (size_t n = 0; n < T.numCell(); ++n) {
		Cell cell = T.cell(n);
		// This calculation should now be done in reaction CalculatePCAPlane
		//cell.calculatePCAPlane(vertexData);
		double concpow = std::pow(cellData[cell.index()][variableIndex(0,0)],parameter(3));
		std::vector<double> normal = cell.getNormalToPCAPlane();
		double norm=0.0;
		for (size_t d=0; d<dimension; ++d)
			norm += normal[d]*normal[d];
		if (norm != 1.0) {
			norm = std::sqrt(norm);
			assert(norm>0.0);
			for (size_t d=0; d<dimension; ++d)
				normal[d] /= norm;
		}
		std::vector<double> cellPos = cell.positionFromVertex(vertexData);
		
		std::vector<double> scNorm(dimension);
		for (size_t d=0; d<dimension; ++d)
			scNorm[d] = cellPos[d];
		if( cellPos[2]<0.0 )
			scNorm[2]=0.0;
		double scalarProd=0.0;
		for (size_t d=0; d<dimension; ++d)
			scalarProd += scNorm[d]*normal[d];
		if (scalarProd<0.0) {
			for (size_t d=0; d<dimension; ++d)
				normal[d] = -normal[d];
		}
		double coeff = parameter(0) + parameter(1)*concpow/(concpow+Kpow);
		for (size_t k=0; k<cell.numVertex(); ++k)
			for (size_t d=0; d<dimension; ++d)
				vertexDerivs[cell.vertex(k)->index()][d] += coeff * normal[d];
	}	
}

DebugReaction::DebugReaction(std::vector<double> &paraValue,
			     std::vector< std::vector<size_t> > &indValue)
{

}

void DebugReaction::derivs(Tissue &T,
			   std::vector< std::vector<double> > &cellData,
			   std::vector< std::vector<double> > &wallData,
			   std::vector< std::vector<double> > &vertexData,
			   std::vector< std::vector<double> > &cellDerivs,
			   std::vector< std::vector<double> > &wallDerivs,
			   std::vector< std::vector<double> > &vertexDerivs)
{
	for (size_t i = 0; i < T.numCell(); ++i) {
		Cell cell = T.cell(i);
		
// 		std::cerr << "Cell " << cell.index() << std::endl;

// 		for (size_t i = 0; i < cell.numVertex(); ++i) {
// 			Vertex *vertex = cell.vertex(i);
// 			std::cerr << "   Vertex: " << i << std::endl;
// 			std::cerr << "      x = " << vertexData[vertex->index()][0] << std::endl;
// 			std::cerr << "      y = " << vertexData[vertex->index()][1] << std::endl;
// 			std::cerr << "      z = " << vertexData[vertex->index()][2] << std::endl;
// 		}

		// This calculation should now be done in reaction CalculatePCAPlane
		//cell.calculatePCAPlane(vertexData);
		
		std::vector<double> N = cell.getNormalToPCAPlane();

		std::cerr << "x = " << N[0] << std::endl;
		std::cerr << "y = " << N[1] << std::endl;
		std::cerr << "z = " << N[2] << std::endl;
	}
}
