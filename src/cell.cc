/**
 * Filename     : cell.cc
 * Description  : A class describing a two-dimensional cell
 * Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
 * Created      : April 2006
 * Revision     : $Id:$
 */

#include <assert.h>
#include "cell.h"
#include "tissue.h"
#include "vertex.h"


//!The empty cel constructor
Cell::Cell() {
  
  mitosisFlag_=0;
}

//!The copy constructor
Cell::Cell( const Cell & cellCopy ) {

  index_ = cellCopy.index();
  id_ = cellCopy.id();
  volume_ = cellCopy.volume();
  mitosisFlag_ = cellCopy.mitosisFlag();
  wall_ = cellCopy.wall();
  vertex_ = cellCopy.vertex();
	variable_ = cellCopy.variable();
}

Cell::Cell(size_t indexVal,std::string idVal) {
  setIndex(indexVal);
  setId(idVal);
}

Cell::~Cell() {}

void Cell::sortWallAndVertex(Tissue &T) {
	
	assert( numWall()==numVertex() );
	
	//std::cerr << "Cell " << index() << std::endl;
	//for( size_t i=0 ; i<numVertex() ; ++i )
	//std::cerr << vertex(i)->index() << std::endl; 
	//for( size_t i=0 ; i<numWall() ; ++i )
	//std::cerr << wall(i)->index() << "\t" << wall(i)->vertex1()->index() << " "
	//					<< wall(i)->vertex2()->index() << std::endl; 
	
	size_t directionalWallFlag=0;
	Wall* tmpDirectionalWall;
	if( T.numDirectionalWall() && T.directionalWall(index())<numWall() ) {
		directionalWallFlag=1;
		tmpDirectionalWall = wall(T.directionalWall(index()));
	}
	std::vector<Wall*> tmpWall( numWall() );
	std::vector<Vertex*> tmpVertex( numVertex() );
	size_t wallIndex=0;
	size_t vertexIndex=0;
	tmpWall[wallIndex] = wall(wallIndex);
	//std::cerr << "Adding wall " << wall(wallIndex)->index() 
	//				<< " at position " << wallIndex << std::endl;
	tmpVertex[vertexIndex] = wall(wallIndex)->vertex1();
	//std::cerr << "Adding vertex " << wall(wallIndex)->vertex1()->index() 
	//				<< " at position " << vertexIndex << std::endl;
	++vertexIndex;
	tmpVertex[vertexIndex] = wall(wallIndex)->vertex2();
	//std::cerr << "Adding vertex " << wall(wallIndex)->vertex2()->index() 
	//				<< " at position " << vertexIndex << std::endl;
	
	while( wallIndex<numWall()-1 ) {
		for( size_t wI=1 ; wI<numWall() ; ++wI ) {
			if( wall(wI) != tmpWall[wallIndex] &&
					wall(wI)->hasVertex( tmpVertex[vertexIndex] ) ) {
				wallIndex++;
				tmpWall[wallIndex]=wall(wI);				
				//std::cerr << "Adding wall " << wall(wI)->index() 
				//				<< " at position " << wallIndex << std::endl;
				if( wallIndex<numWall()-1 ) {
					if( tmpVertex[vertexIndex]==wall(wI)->vertex1() )  {
						++vertexIndex;
						tmpVertex[vertexIndex]=wall(wI)->vertex2();
						//std::cerr << "Adding vertex " << wall(wI)->vertex2()->index() 
						//				<< " at position " << vertexIndex << std::endl;
						
					}
					else if( tmpVertex[vertexIndex]==wall(wI)->vertex2() ) {
						++vertexIndex;
						tmpVertex[vertexIndex]=wall(wI)->vertex1();
						//std::cerr << "Adding vertex " << wall(wI)->vertex1()->index() 
						//				<< " at position " << vertexIndex << std::endl;
					}
					else {
						std::cerr << "Cell::sortWallAndVertex() "
											<< "Wrong vertex index in wall." << std::endl;
						exit(-1);
					}
				}
				break;
			}
		}
	}
	assert( wallIndex==vertexIndex );
	
	setWall( tmpWall );
	setVertex( tmpVertex );
	
	if( directionalWallFlag ) {
		//Find the wall in the new list
		size_t foundWall=0;
		for( size_t k=0; k<numWall(); ++k )
			if( wall(k)==tmpDirectionalWall ) {
				foundWall++;
				T.setDirectionalWall(index(),k);
			}
		if( foundWall != 1 ) {
			std::cerr << "Cell::sortWallAndVertex() Wrong number of "
								<< "directional walls found (." << foundWall 
								<< ")." << std::endl;
			exit(-1);
		}
		assert( numWall()>T.directionalWall(index()) );
	}
	//std::cerr << "Cell " << index() << std::endl;
	//for( size_t i=0 ; i<numVertex() ; ++i )
	//std::cerr << vertex(i)->index() << std::endl; 
	//for( size_t i=0 ; i<numWall() ; ++i )
	//std::cerr << wall(i)->index() << "\t" << wall(i)->vertex1()->index() << " "
	//					<< wall(i)->vertex2()->index() << std::endl; 
	
}

//!Calculates the volume from vertex positions
double Cell::calculateVolume( size_t signFlag ) 
{	
  assert( numVertex() );	
  size_t dimension = vertex(0)->numPosition();
	if( dimension == 2 ) {
		//Assuming vertices are sorted
		double tmpVolume=0.0;
		for( size_t k=0 ; k<numVertex() ; ++k ) {
			size_t kk = (k+1)%(numVertex());
			tmpVolume += vertex(k)->position(0)*vertex(kk)->position(1) -
				vertex(k)->position(1)*vertex(kk)->position(0);
		}
		tmpVolume *= 0.5;
		volume_ = std::fabs(tmpVolume);
		if( signFlag ) 
			return tmpVolume;
		else
			return volume_;
	}
	else {
		std::cerr << "Cell::calculateVolume(signFlag) "
							<< " only implemented for dimension=2." << std::endl;
		exit(-1);
	}
}

//!Calculates the volume from vertex positions
double Cell::calculateVolume( std::vector< std::vector<double> > 
															&vertexData, size_t signFlag ) {
	
  assert( numVertex() );
  size_t dimension = vertex(0)->numPosition();
  
	if( dimension == 2 ) {
		//Assuming vertices are sorted
		double tmpVolume=0.0;
		for( size_t k=0 ; k<numVertex() ; ++k ) {
			size_t v1I = vertex(k)->index();
			size_t v2I = vertex((k+1)%(numVertex()))->index();
			tmpVolume += vertexData[v1I][0]*vertexData[v2I][1]-
				vertexData[v1I][1]*vertexData[v2I][0];
		}
		tmpVolume *= 0.5;
		volume_ = std::fabs(tmpVolume);
		if( signFlag ) 
			return tmpVolume;
		else
			return volume_;
	}
	else {
		//Caveat:Old version to be changed to projected version of 2D variant
		//Calculate cell position from vertices
		std::vector<double> xCenter = positionFromVertex(vertexData);
		assert( xCenter.size()==dimension );
		
		//Calculate volume from vertex positions for each wall
		volume_=0.0;
		for( size_t k=0 ; k<numWall() ; ++k ) {
			Wall *tmpWall = wall(k);
			size_t v1I = tmpWall->vertex1()->index();
			size_t v2I = tmpWall->vertex2()->index();
			std::vector<double> n(dimension);
			double b=0.0;
			for( size_t d=0 ; d<dimension ; ++d ) {
				n[d] = vertexData[v2I][d]-vertexData[v1I][d];
				b += n[d]*n[d];
			}
			assert( b>0.0 );
			b = std::sqrt(b);
			for( size_t d=0 ; d<dimension ; ++d )
				n[d] /= b;
			double h = (xCenter[0]-vertexData[v1I][0])*(xCenter[0]-vertexData[v1I][0])
				+ (xCenter[1]-vertexData[v1I][1])*(xCenter[1]-vertexData[v1I][1])
				- ( n[0]*(xCenter[0]-vertexData[v1I][0])+n[1]*(xCenter[1]-vertexData[v1I][1]) )*
				( n[0]*(xCenter[0]-vertexData[v1I][0])+n[1]*(xCenter[1]-vertexData[v1I][1]) );
			assert( h>0.0 );
			h = std::sqrt(h);
			volume_ += b*h;
		}
		volume_ *= 0.5;
		return volume_;
	}
}

//!Calculates the cell position from average vertex position
std::vector<double> Cell::positionFromVertex() 
{
	assert( numVertex() );
  size_t dimension=vertex(0)->numPosition();
	assert(dimension==2);
  std::vector<double> pos( dimension );
	double area = calculateVolume(1);
	
  for (size_t i=0; i<numVertex(); ++i) {
		size_t ii = (i+1)%(numVertex());
		double factor = vertex(i)->position(0)*vertex(ii)->position(1)-
			vertex(ii)->position(0)*vertex(i)->position(1);
    for (size_t d=0; d<dimension; ++d) {
      pos[d] += factor*(vertex(i)->position(d)+vertex(ii)->position(d));
    }
  }
  for (size_t d=0; d<dimension; ++d)
    pos[d] /= 6*area;
  
  return pos;
}

//!Calculates the cell position from average vertex position
std::vector<double> Cell::
positionFromVertex( std::vector< 
										std::vector<double> > &vertexData ) 
{  
	assert( numVertex() );
  size_t dimension=vertexData[0].size();
	assert(dimension==2);
  std::vector<double> pos( dimension );
	double area = calculateVolume(vertexData,1);
	
  for (size_t i=0; i<numVertex(); ++i) {
		size_t vI = vertex(i)->index();
		size_t vIPlus = vertex( (i+1)%(numVertex()) )->index();
		double factor = vertexData[vI][0]*vertexData[vIPlus][1]-
			vertexData[vIPlus][0]*vertexData[vI][1];
    for (size_t d=0; d<dimension; ++d) {
      pos[d] += factor*(vertexData[vI][d]+vertexData[vIPlus][d]);
    }
  }
  for (size_t d=0; d<dimension; ++d)
    pos[d] /= 6*area;
  
  return pos;
}

