//
// Filename     : cell.cc
// Description  : A class describing a two-dimensional cell
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : April 2006
// Revision     : $Id:$
//

#include <assert.h>
#include "cell.h"
#include <limits>
#include "tissue.h"
#include "vertex.h"
#include "myMath.h"
#include "myRandom.h"

Cell::Cell() {
  
  mitosisFlag_=0;
}

Cell::Cell( const Cell & cellCopy ) {

  index_ = cellCopy.index();
  id_ = cellCopy.id();
  volume_ = cellCopy.volume();
  mitosisFlag_ = cellCopy.mitosisFlag();
  wall_ = cellCopy.wall();
  vertex_ = cellCopy.vertex();
  variable_ = cellCopy.variable();
  E_ = cellCopy.getPCAPlane();
}

Cell::Cell(size_t indexVal,std::string idVal) {
  setIndex(indexVal);
  setId(idVal);
}

Cell::~Cell() {}

void Cell::setVariable( std::vector<double> &val ) 
{
  if (numVariable() != val.size()) {
    std::cerr << "Cell::setVariable(vector) Not the same number of variables in the cell as in the given vector."
	      << std::endl;
    exit(EXIT_FAILURE);
  }
  for (size_t i=0; i<val.size(); ++i)
    variable_[i]=val[i];
}

void Cell::sortWallAndVertexOld(Tissue &T) {
	
	assert( numWall()==numVertex() );
	
	//std::cerr << "Cell " << index() << std::endl;
	//for( size_t i=0 ; i<numVertex() ; ++i )
	//std::cerr << vertex(i)->index() << std::endl; 
	//for( size_t i=0 ; i<numWall() ; ++i )
	//std::cerr << wall(i)->index() << "\t" << wall(i)->vertex1()->index() << " "
	//					<< wall(i)->vertex2()->index() << std::endl; 
	
	size_t directionalWallFlag=0;
	Wall* tmpDirectionalWall=NULL;
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
		assert(tmpDirectionalWall);
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

void Cell::sortWallAndVertex(Tissue &T) {
	
	assert( numWall()==numVertex() );
	
// 	std::cerr << "Cell " << index() << std::endl;
// 	for( size_t i=0 ; i<numVertex() ; ++i )
// 		std::cerr << vertex(i)->index() << " "; 
// 	std::cerr << std::endl;
// 	for( size_t i=0 ; i<numWall() ; ++i )
// 		std::cerr << wall(i)->index() << "\t" << wall(i)->vertex1()->index() << " "
// 							<< wall(i)->vertex2()->index() << std::endl; 
	
	size_t directionalWallFlag=0;
	Wall* tmpDirectionalWall=NULL;
	if( T.numDirectionalWall() && T.directionalWall(index())<numWall() ) {
		directionalWallFlag=1;
		tmpDirectionalWall = wall(T.directionalWall(index()));
	}
	std::vector<Wall*> tmpWall( numWall() );
	std::vector<Vertex*> tmpVertex( numVertex() );
	//Find initial wall and initial vertices
	size_t wallIndex=0,wallIndexStart=0,vertexIndex=0;
	size_t foundCellSortFlag=0;
	while (wallIndex<numWall() && !wall(wallIndex)->cellSort1())
	  ++wallIndex;
	if (wallIndex<numWall()) {
	  foundCellSortFlag=1;
	  wallIndexStart=vertexIndex=wallIndex;
	}
	else
	  wallIndex=wallIndexStart=0;
	
	// 	if (foundCellSortFlag)
	// 		std::cerr << "Found cell sort flag." << std::endl;
// 	std::cerr << "WallIndexStart: " << wallIndexStart << std::endl;
	unsigned int numWallSorted=0,numVertexSorted=0;
	// Set sorting order for first wall
	if (!foundCellSortFlag) {//No previous information about sorting
		tmpWall[wallIndex] = wall(wallIndex);
		++numWallSorted;
// 		std::cerr << "Adding wall " << wall(wallIndex)->index() 
// 							<< " at position " << wallIndex << std::endl;
		tmpVertex[vertexIndex] = wall(wallIndex)->vertex1();
		++numVertexSorted;
// 		std::cerr << "Adding vertex " << wall(wallIndex)->vertex1()->index() 
// 							<< " at position " << vertexIndex << std::endl;
		++vertexIndex;
		tmpVertex[vertexIndex] = wall(wallIndex)->vertex2();
		++numVertexSorted;
// 		std::cerr << "Adding vertex " << wall(wallIndex)->vertex2()->index() 
// 							<< " at position " << vertexIndex << std::endl;
		if (wall(wallIndex)->cell1()==this) {
			wall(wallIndex)->setCellSort1(1);
			if (wall(wallIndex)->cell2()!=T.background())
				wall(wallIndex)->setCellSort2(-1);
		}
		else if (wall(wallIndex)->cell2()==this) {
			wall(wallIndex)->setCellSort2(1);
			if (wall(wallIndex)->cell1()!=T.background())
				wall(wallIndex)->setCellSort1(-1);
		}
		else {
			std::cerr << "cell.sortWallAndVertex() Wall to be sorted does not belong to the cell"
								<< std::endl;
			exit(EXIT_FAILURE);
		}
	}
	else {//First wall chosen from earlier sorting information
		tmpWall[wallIndex] = wall(wallIndex);
		++numWallSorted;
// 		std::cerr << "Adding wall " << wall(wallIndex)->index() 
// 							<< " at position " << wallIndex << std::endl;
		if (wall(wallIndex)->cell1()==this) {
			if (wall(wallIndex)->cellSort1()==1) {
				tmpVertex[vertexIndex] = wall(wallIndex)->vertex1();
				++numVertexSorted;
// 				std::cerr << "Adding vertex " << wall(wallIndex)->vertex1()->index() 
// 									<< " at position " << vertexIndex << std::endl;
				++vertexIndex;
				vertexIndex = vertexIndex%numVertex();
				tmpVertex[vertexIndex] = wall(wallIndex)->vertex2();
				++numVertexSorted;
// 				std::cerr << "Adding vertex " << wall(wallIndex)->vertex2()->index() 
// 									<< " at position " << vertexIndex << std::endl;
			}
			else if (wall(wallIndex)->cellSort1()==-1) {
				tmpVertex[vertexIndex] = wall(wallIndex)->vertex2();
				++numVertexSorted;
// 				std::cerr << "Adding vertex " << wall(wallIndex)->vertex2()->index() 
// 									<< " at position " << vertexIndex << std::endl;				
				++vertexIndex;
				vertexIndex = vertexIndex%numVertex();
				tmpVertex[vertexIndex] = wall(wallIndex)->vertex1();
				++numVertexSorted;
// 				std::cerr << "Adding vertex " << wall(wallIndex)->vertex1()->index() 
// 									<< " at position " << vertexIndex << std::endl;
			}
			else {
				std::cerr << "cell.sortWallAndVertex() Wall to be sorted marked as but is not sorted."
									<< std::endl;
				exit(EXIT_FAILURE);
			}
		}
		else if (wall(wallIndex)->cell2()==this) {
			if (wall(wallIndex)->cellSort2()==1) {
				tmpVertex[vertexIndex] = wall(wallIndex)->vertex1();
				++numVertexSorted;
// 				std::cerr << "Adding vertex " << wall(wallIndex)->vertex1()->index() 
// 									<< " at position " << vertexIndex << std::endl;
				++vertexIndex;
				vertexIndex = vertexIndex%numVertex();
				tmpVertex[vertexIndex] = wall(wallIndex)->vertex2();
				++numVertexSorted;
// 				std::cerr << "Adding vertex " << wall(wallIndex)->vertex2()->index() 
// 									<< " at position " << vertexIndex << std::endl;
			}
			else if (wall(wallIndex)->cellSort2()==-1) {
				tmpVertex[vertexIndex] = wall(wallIndex)->vertex2();
				++numVertexSorted;
// 				std::cerr << "Adding vertex " << wall(wallIndex)->vertex2()->index() 
// 									<< " at position " << vertexIndex << std::endl;				
				++vertexIndex;
				vertexIndex = vertexIndex%numVertex();
				tmpVertex[vertexIndex] = wall(wallIndex)->vertex1();
				++numVertexSorted;
// 				std::cerr << "Adding vertex " << wall(wallIndex)->vertex1()->index() 
// 									<< " at position " << vertexIndex << std::endl;
			}
			else {
				std::cerr << "cell.sortWallAndVertex() Wall to be sorted marked as but is not sorted."
									<< std::endl;
				exit(EXIT_FAILURE);
			}			
		}
		else {
			std::cerr << "cell.sortWallAndVertex() Wall to be sorted does not belong to the cell"
								<< std::endl;
			exit(EXIT_FAILURE);
		}
	}
	
	//size_t wallIndexMod = wallIndex%numWall();
	// Sort the remaining walls
	while( wallIndex<wallIndexStart+numWall()-1 ) {
		for( size_t wI=0 ; wI<numWall() ; ++wI ) {
			if( wall(wI) != tmpWall[wallIndex%numWall()] &&
					wall(wI)->hasVertex( tmpVertex[vertexIndex] ) ) {
				++wallIndex;
				tmpWall[wallIndex%numWall()]=wall(wI);				
				++numWallSorted;
// 				std::cerr << "Adding wall " << wall(wI)->index() 
// 									<< " at position " << wallIndex%numWall() << std::endl;
				if( wallIndex<wallIndexStart+numWall()-1 ) {
					if( tmpVertex[vertexIndex]==wall(wI)->vertex1() )  {
						++vertexIndex;
						vertexIndex=vertexIndex%numVertex();
						tmpVertex[vertexIndex]=wall(wI)->vertex2();
						++numVertexSorted;
// 						std::cerr << "Adding vertex " << wall(wI)->vertex2()->index() 
// 											<< " at position " << vertexIndex << std::endl;						
					}
					else if( tmpVertex[vertexIndex]==wall(wI)->vertex2() ) {
						++vertexIndex;
						vertexIndex=vertexIndex%numVertex();
						tmpVertex[vertexIndex]=wall(wI)->vertex1();
						++numVertexSorted;
// 						std::cerr << "Adding vertex " << wall(wI)->vertex1()->index() 
// 											<< " at position " << vertexIndex << std::endl;
					}
					else {
						std::cerr << "Cell::sortWallAndVertex() "
							  << "Wrong vertex index in wall." << std::endl;
						exit(-1);
					}
				}
				else {//make sure vertexIndex points at initial vertex (which is not added again)
					++vertexIndex;
					vertexIndex=vertexIndex%numVertex();
				}
				// Check if this wall has marked sorting and either check if it complies with
				// current sorting or add sorting variables
				if (wall(wI)->cell1()==this) {
				  if (wall(wI)->cellSort1()==1 ) {//sorted along wall vertices
				    if (wall(wI)->vertex2() != tmpVertex[vertexIndex] ) {
				      std::cerr << "Cell::sortWallAndVertex() "
						<< "Current sorting direction for wall does not comply with previous."
						<< std::endl;
				      std::cerr << "Cell sort: " << wall(wI)->cellSort1() << " wall vertices: "
						<< wall(wI)->vertex1()->index() << "," 
						<< wall(wI)->vertex2()->index() << std::endl
						<< "vertexIndex: " << vertexIndex << std::endl;
				      std::cerr << "tmpVertex ";
				      for (size_t k=0; k<tmpVertex.size(); ++k)
					if (tmpVertex[k]) 
					  std::cerr << tmpVertex[k]->index() << " (" << k << ") ";
				      std::cerr << std::endl;
				      std::cerr << "tmpWall ";
				      for (size_t k=0; k<tmpWall.size(); ++k)
					if (tmpWall[k]) 
					  std::cerr << tmpWall[k]->index() << " (" << k << ") ";
				      std::cerr << std::endl;
				      exit(-1);
				    }
				  }
				  else if (wall(wI)->cellSort1()==-1 ) {//sorted 'against' wall vertices
				    if (wall(wI)->vertex1() != tmpVertex[vertexIndex] ) {
				      std::cerr << "Cell::sortWallAndVertex() "
						<< "Current sorting direction for wall does not comply with previous."
						<< std::endl;
				      std::cerr << "Cell sort: " << wall(wI)->cellSort1() << " wall vertices: "
						<< wall(wI)->vertex1()->index() << "," 
						<< wall(wI)->vertex2()->index() << std::endl
						<< "vertexIndex: " << vertexIndex << std::endl;
				      std::cerr << "tmpVertex ";
				      for (size_t k=0; k<tmpVertex.size(); ++k)
					if (tmpVertex[k]) 
					  std::cerr << tmpVertex[k]->index() << " (" << k << ") ";
				      std::cerr << std::endl;
				      std::cerr << "tmpWall ";
				      for (size_t k=0; k<tmpWall.size(); ++k)
					if (tmpWall[k]) 
					  std::cerr << tmpWall[k]->index() << " (" << k << ") ";
				      std::cerr << std::endl;
				      exit(-1);
				    }
				  }
				  else {//not sorted before
				    if (wall(wI)->vertex2()==tmpVertex[vertexIndex] ) {
				      wall(wI)->setCellSort1(1);
				      if (wall(wI)->cell2()!=T.background())
					wall(wI)->setCellSort2(-1);
				    }
				    else if (wall(wI)->vertex1()==tmpVertex[vertexIndex] ) {
				      wall(wI)->setCellSort1(-1);
				      if (wall(wI)->cell2()!=T.background())
					wall(wI)->setCellSort2(1);
				    }
				    else {
				      std::cerr << "Cell::sortWallAndVertex() "
						<< "tmpVertex not in wall sorted!" << std::endl;
				      exit(-1);
				    }
				  }					
				}
				else if (wall(wI)->cell2()==this) {
				  if (wall(wI)->cellSort2()==1 ) {//sorted along wall vertices
				    if (wall(wI)->vertex2() != tmpVertex[vertexIndex] ) {
				      std::cerr << "Cell::sortWallAndVertex() "
						<< "Current sorting direction for wall does not comply with previous."
						<< std::endl;
				      std::cerr << "Cell sort: " << wall(wI)->cellSort2() << " wall vertices: "
						<< wall(wI)->vertex1()->index() << "," 
						<< wall(wI)->vertex2()->index() << std::endl
						<< "vertexIndex: " << vertexIndex << std::endl;
				      std::cerr << "tmpVertex ";
				      for (size_t k=0; k<tmpVertex.size(); ++k)
					if (tmpVertex[k]) 
					  std::cerr << tmpVertex[k]->index() << " (" << k << ") ";
				      std::cerr << std::endl;
				      std::cerr << "tmpWall ";
				      for (size_t k=0; k<tmpWall.size(); ++k)
					if (tmpWall[k]) 
					  std::cerr << tmpWall[k]->index() << " (" << k << ") ";
				      std::cerr << std::endl;
				      exit(-1);
				    }
				  }
				  else if (wall(wI)->cellSort2()==-1 ) {//sorted 'against' wall vertices
						if (wall(wI)->vertex1() != tmpVertex[vertexIndex] ) {
							std::cerr << "Cell::sortWallAndVertex() "
												<< "Current sorting direction for wall does not comply with previous."
												<< std::endl;
							std::cerr << "Cell sort: " << wall(wI)->cellSort2() << " wall vertices: "
												<< wall(wI)->vertex1()->index() << "," 
												<< wall(wI)->vertex2()->index() << std::endl
												<< "vertexIndex: " << vertexIndex << std::endl;
							std::cerr << "tmpVertex ";
							for (size_t k=0; k<tmpVertex.size(); ++k)
								if (tmpVertex[k]) 
									std::cerr << tmpVertex[k]->index() << " (" << k << ") ";
							std::cerr << std::endl;
							std::cerr << "tmpWall ";
							for (size_t k=0; k<tmpWall.size(); ++k)
								if (tmpWall[k]) 
									std::cerr << tmpWall[k]->index() << " (" << k << ") ";
							std::cerr << std::endl;
							exit(-1);
						}
					}
					else {//not sorted before
						if (wall(wI)->vertex2()==tmpVertex[vertexIndex] ) {
							wall(wI)->setCellSort2(1);
							if (wall(wI)->cell1()!=T.background())
								wall(wI)->setCellSort1(-1);
						}
						else if (wall(wI)->vertex1()==tmpVertex[vertexIndex] ) {
							wall(wI)->setCellSort2(-1);
							if (wall(wI)->cell1()!=T.background())
								wall(wI)->setCellSort1(1);
						}
						else {
							std::cerr << "Cell::sortWallAndVertex() "
												<< "tmpVertex not in wall sorted!" << std::endl;
							exit(-1);
						}
					}					
				}					
				else {
					std::cerr << "cell.sortWallAndVertex() Wall to be sorted does not belong to the cell"
										<< std::endl;
					exit(EXIT_FAILURE);
				}
				break;
			}
		}
	}
// 	std::cerr << "numWallSorted: " << numWallSorted << " (" << numWall() << ") numVertexSorted: "
// 						<< numVertexSorted << " (" << numVertex() << ")" << std::endl;
	if (numWallSorted!=numWall() || numVertexSorted!=numVertex()) {
		std::cerr << "Cell " << index() << " has sorted " << numWallSorted 
							<< " walls out of " << numWall()
							<<" and " << numVertexSorted << " vertices out of " 
							<< numVertex() << std::endl;
		for (size_t k=0; k<tmpWall.size(); ++k)
			std::cerr << tmpWall[k]->index() << " ";
		std::cerr << std::endl;
		for (size_t k=0; k<tmpVertex.size(); ++k)
			std::cerr << tmpVertex[k]->index() << " ";
		std::cerr << std::endl;
		exit(-1);
	}
	assert( numWallSorted==numVertexSorted );
	
	setWall( tmpWall );
	setVertex( tmpVertex );
	
	if( directionalWallFlag ) {
		//Find the wall in the new list
		assert(tmpDirectionalWall);
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
	else if (dimension==3) {
		//Caveat:Old version to be changed to projected version of 2D variant
		//Calculate cell position from vertices
		std::vector<double> xCenter = positionFromVertex();
		assert( xCenter.size()==dimension );
		
		//Calculate volume from vertex positions for each wall
		volume_=0.0;
		for( size_t k=0 ; k<numWall() ; ++k ) {
			std::vector<double> r1(dimension),r2(dimension);
			for( size_t d=0 ; d<dimension ; ++d ) {
				r1[d] = wall(k)->vertex1()->position(d) - xCenter[d];
				r2[d] = wall(k)->vertex2()->position(d) - xCenter[d];
			}
			double triArea=0.0;
			for( size_t d=0 ; d<dimension ; ++d ) {
				size_t d1 = (d+1)%dimension;
				size_t d2 = (d+2)%dimension;
				double r1crossr2 = r1[d1]*r2[d2]-r1[d2]*r2[d1];
				triArea += r1crossr2*r1crossr2;
			}
			triArea = std::sqrt(triArea);
			volume_ += 0.5*triArea;
		}
		return volume_;
	}
	else {
		std::cerr << "Cell::calculateVolume() Only applicable for two or three"
							<< " dimensions." << std::endl;
		exit(-1);
	}
}

double Cell::calculateVolume( DataMatrix 				                                  &vertexData, size_t signFlag ) 
{	
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
  else if (dimension==3) {
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
      std::vector<double> r1(dimension),r2(dimension);
      for( size_t d=0 ; d<dimension ; ++d ) {
	r1[d] = vertexData[v1I][d]-xCenter[d];
	r2[d] = vertexData[v2I][d]-xCenter[d];
      }
      double triArea=0.0;
      for( size_t d=0 ; d<dimension ; ++d ) {
	size_t d1 = (d+1)%dimension;
	size_t d2 = (d+2)%dimension;
	double r1crossr2 = r1[d1]*r2[d2]-r1[d2]*r2[d1];
	triArea += r1crossr2*r1crossr2;
      }
      triArea = std::sqrt(triArea);
      volume_ += 0.5*triArea;
    }
    return volume_;
  }
  else {
    std::cerr << "Cell::calculateVolume(vertexData) Only applicable for two or three"
	      << " dimensions." << std::endl;
    exit(-1);
  }
}

double Cell::
calculateVolumeCenterTriangulation( DataMatrix 
				    &vertexData,
				    DataMatrix &cellData,
				    size_t centerIndex)
{
  double area = 0.0;
  size_t dimension=vertexData[0].size();
  std::vector<double> p1(dimension);
  for (size_t d=0; d<dimension; ++d)
    p1[d] = cellData[index()][centerIndex+d];
  for (size_t i=0; i<numVertex()-1; ++i) {    
    area += myMath::areaTriangle(p1,vertexData[vertex(i)->index()],vertexData[vertex(i+1)->index()]);
  }
  return area;
}

std::vector<double> Cell::positionFromVertex() 
{
	assert( numVertex() );
	size_t dimension=vertex(0)->numPosition();
	std::vector<double> pos(dimension, 0);
	if (dimension==2) {
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
	}
	else if (dimension==3) {
		double sumWallLength=0.0;
		for (size_t w=0; w<numWall(); ++w) {
			double wallLength=0.0;
			for (size_t d=0; d<dimension; ++d)
				wallLength += (wall(w)->vertex1()->position(d) -
											 wall(w)->vertex2()->position(d) ) *
					(wall(w)->vertex1()->position(d) -
					 wall(w)->vertex2()->position(d) );
			wallLength = std::sqrt(wallLength);
			sumWallLength += wallLength;
			for (size_t d=0; d<dimension; ++d)
				pos[d] += 0.5 * (wall(w)->vertex1()->position(d) +
												 wall(w)->vertex2()->position(d) ) * wallLength;
		}
		for (size_t d=0; d<dimension; ++d)
			pos[d] /= sumWallLength;
	}
	else {
		std::cerr << "Cell::positionFromVertex() only implemented for two and"
							<< " three dimensions" << std::endl;
		exit(-1);
	}
	return pos;
}

std::vector<double> Cell::
positionFromVertex(DataMatrix &vertexData) 
{  
	assert( numVertex() );
	size_t dimension=vertexData[0].size();
	std::vector<double> pos (dimension, 0);
	if (dimension==2) {
		double area = calculateVolume(vertexData,1);
		
		for (size_t i=0; i<numVertex(); ++i) {
			
			size_t vI = vertex(i)->index();
			// 		std::cerr << " *** Vertex " << i << " = (" << vertexData[vI][0] << ", " << vertexData[vI][1] << ")" << std::endl;
			
			size_t vIPlus = vertex( (i+1)%(numVertex()) )->index();
			double factor = vertexData[vI][0]*vertexData[vIPlus][1]-
				vertexData[vIPlus][0]*vertexData[vI][1];
			for (size_t d=0; d<dimension; ++d) {
				pos[d] += factor*(vertexData[vI][d]+vertexData[vIPlus][d]);
			}
		}
		for (size_t d=0; d<dimension; ++d)
			pos[d] /= 6*area;
	}
	else if (dimension==3) {
		double sumWallLength=0.0;
		for (size_t w=0; w<numWall(); ++w) {
			double wallLength=0.0;
			for (size_t d=0; d<dimension; ++d)
				wallLength += (vertexData[wall(w)->vertex1()->index()][d] -
											 vertexData[wall(w)->vertex2()->index()][d] ) *
					(vertexData[wall(w)->vertex1()->index()][d] -
					 vertexData[wall(w)->vertex2()->index()][d] );
			wallLength = std::sqrt(wallLength);
			sumWallLength += wallLength;
			for (size_t d=0; d<dimension; ++d)
				pos[d] += 0.5*(vertexData[wall(w)->vertex1()->index()][d] +
											 vertexData[wall(w)->vertex2()->index()][d] ) *
					wallLength;
		}
		for (size_t d=0; d<dimension; ++d)
			pos[d] /= sumWallLength;
	}
	else {
		std::cerr << "Cell::positionFromVertex() only implemented for two and"
							<< " three dimensions" << std::endl;
		exit(-1);
	}
	return pos;
}

std::vector<double> Cell::randomPositionInCell(const DataMatrix &vertexData, const int numberOfTries) 
{
	typedef std::vector<double> Vector;

	const size_t dimensions = vertexData[0].size();

	if (dimensions != 2)
	{
		std::cerr << "Cell::randomPositionInCell only supports two dimensions.\n";
		std::exit(EXIT_FAILURE);
	}

	double xmin = std::numeric_limits<double>::max();
	double xmax = std::numeric_limits<double>::min();
	double ymin = std::numeric_limits<double>::max();
	double ymax = std::numeric_limits<double>::min();

	for (size_t vertexIndex = 0; vertexIndex < numVertex(); ++vertexIndex)
	{
		const Vertex &v = *vertex(vertexIndex);

		const double &x = vertexData[v.index()][0];

		if (x < xmin)
		{
			xmin = x;
		}
		if (x > xmax)
		{
			xmax = x;
		}

		const double &y = vertexData[v.index()][1];

		if (y < ymin)
		{
			ymin = y;
		}
		if (y > ymax)
		{
			ymax = y;
		}
	}

	for (int tryCounter = 0; tryCounter < numberOfTries; ++tryCounter)
	{
		double rx = xmin + (xmax - xmin) * myRandom::Rnd();
		double ry = ymin + (ymax - ymin) * myRandom::Rnd();

		const size_t numberOfVertices = numVertex();

		signed int sign = 0;

		bool success = true;

		for (size_t vertexIndex = 0; vertexIndex < numberOfVertices; ++vertexIndex)
		{
			const Vertex &v1 = *vertex(vertexIndex);
			const Vertex &v2 = *vertex((vertexIndex + 1) % numberOfVertices);

			const double vx = vertexData[v2.index()][0] - vertexData[v1.index()][0];
			const double vy = vertexData[v2.index()][1] - vertexData[v1.index()][1];
				
			const double dx = rx - vertexData[v1.index()][0];
			const double dy = ry - vertexData[v1.index()][1];

			const signed int s = myMath::sign(vx * dy - vy * dx);

			if (!sign)
			{
				sign = s;
			}
			else
			{
				if (sign == s)
				{
					continue;
				}
				else
				{
					success = false;
					break;
				}
			}
		}

		if (success)
		{
			std::vector<double> result(2);

			result[0] = rx;
			result[1] = ry;

			return result;
		}
	}

	throw FailedToFindRandomPositionInCellException();
}

void Cell::calculatePCAPlane(DataMatrix &vertexData)
{
	size_t dimensions = vertexData[0].size();
	size_t numberOfVertices = vertex_.size();

	// Copy vertex data to temporary container and calculate mean values.
 
	DataMatrix vertices(numberOfVertices);
	std::vector<double> mean(dimensions, 0.0);

	for (size_t i = 0; i < numberOfVertices; ++i) {
		Vertex *v = vertex(i);
		vertices[i].resize(dimensions);
		for (size_t j = 0; j < dimensions; ++j) {
			vertices[i][j] = vertexData[v->index()][j];
			mean[j] += vertexData[v->index()][j];
		}
	}

	for (size_t i = 0; i < dimensions; ++i) {
		mean[i] /= numberOfVertices;
	}

	// Subtract mean from data to get an expectation value equal to zero.

	for (size_t i = 0; i < numberOfVertices; ++i) {
		for (size_t j = 0; j < dimensions; ++j) {
			vertices[i][j] -= mean[j];
		}
	}

	// Calculate the correlation matrix.
	
	DataMatrix R(dimensions);

	for (size_t i = 0; i < dimensions; ++i) {
	  R[i].resize(dimensions);
	  for (size_t j = 0; j < dimensions; ++j) {
	    R[i][j] = 0.0;
	  }
	}

	for (size_t k = 0; k < dimensions; ++k) {
		for (size_t l = 0; l < dimensions; ++l) {
			for (size_t i = 0; i < numberOfVertices; ++i) {
				R[k][l] += vertices[i][k] * vertices[i][l];
			}
			R[k][l] /= numberOfVertices;
		}
	}

	// Find the eigenvectors with the two greatests corresponding eigenvalues.

	DataMatrix candidates;

	DataMatrix V;
	std::vector<double> d;

	myMath::jacobiTransformation(R , V, d);
       
	double max = 0.0;
	size_t max1 = d.size();
	size_t max2 = d.size();

	max = 0.0;
	for (size_t i = 0; i < d.size(); ++i) {
		if (std::abs(d[i]) >= max) {
			max1 = i;
			max = std::abs(d[i]);
		}
	}

	max = 0.0;
	for (size_t i = 0; i < d.size(); ++i) {
		if (std::abs(d[i]) >= max && i != max1) {
			max2 = i;
			max = std::abs(d[i]);
		}
	}

	if (max1 == d.size() || max2 == d.size()) {
		std::cerr << "Cell::calculatePCAPlane(): Unexpected behaviour." << std::endl;
		exit(EXIT_FAILURE);
	}

 	// Find orthonormal basis.

	E_.resize(2);
	E_[0].resize(dimensions);
	E_[1].resize(dimensions);

	for (size_t i = 0; i < dimensions; ++i) {
		//		E_[0][i] = V[max1][i];
		E_[0][i] = V[i][max1];
	}
	
	double s = 0.0;
	
	double numerator = 0.0;
	double denominator = 0.0;
	for (size_t i = 0; i < dimensions; ++i) {
		// 		numerator += V[max1][i] * V[max2][i];
		// 		denominator += V[max1][i] * V[max1][i];
		numerator += V[i][max1] * V[i][max2];
		denominator += V[i][max1] * V[i][max1];
	}
	assert(denominator>0.0);
	s = -numerator / denominator;
	
	for (size_t i = 0; i < dimensions; ++i) {
		// 		E_[1][i] = s * E_[0][i] + V[max2][i];
		E_[1][i] = s * E_[0][i] + V[i][max2];
	}

	for (size_t i = 0; i < E_.size(); ++i) {
		double sum = 0.0;
		for (size_t j = 0; j < dimensions; ++j) {
			sum += E_[i][j] * E_[i][j];
		}
		for (size_t j = 0; j < dimensions; ++j) {
			E_[i][j] /= std::sqrt(sum);
		}
	}	
}

DataMatrix Cell::getPCAPlane(void) const
{
	return E_;
}

std::vector< std::pair<double, double> > Cell::
projectVerticesOnPCAPlane(DataMatrix &vertexData)
{
	if (E_.size()!=2 || E_[0].size() != 3) {
		std::cerr << "Cell::projectVerticesOnPCAPlane(): "
							<< "PCAPlane not calculated or dimension not equal to three." << std::endl;
		exit(EXIT_FAILURE);
	}
	size_t dimensions = vertexData[0].size();
	size_t numberOfVertices = vertex_.size();	

	// Copy vertex data to temporary container and calculate mean values.
 
	DataMatrix vertices(numberOfVertices);
	std::vector<double> mean(dimensions, 0.0);
	
	for (size_t i = 0; i < numberOfVertices; ++i) {
		Vertex *v = vertex(i);
		vertices[i].resize(dimensions);
		for (size_t j = 0; j < dimensions; ++j) {
			vertices[i][j] = vertexData[v->index()][j];
			mean[j] += vertexData[v->index()][j];
		}
	}
	
	for (size_t i = 0; i < dimensions; ++i) {
		mean[i] /= numberOfVertices;
	}
	
	// Subtract mean from data to get an expectation value equal to zero.
	
	for (size_t i = 0; i < numberOfVertices; ++i) {
		for (size_t j = 0; j < dimensions; ++j) {
			vertices[i][j] -= mean[j];
		}
	}
	
	std::vector< std::pair<double, double> > coordinates;
	
	for (size_t i = 0; i < numberOfVertices; ++i) {
		std::pair<double, double> tmp(0.0, 0.0);
		for (size_t j = 0; j < dimensions; ++j) {
			tmp.first += E_[0][j] * vertices[i][j];
			tmp.second += E_[1][j] * vertices[i][j];
		}
		coordinates.push_back(tmp);
	}
	
	return coordinates;
}

std::vector<double> Cell::getNormalToPCAPlane(void)
{
  if (E_.size()!=2 || E_[0].size() != 3) {
    std::cerr << "Cell::getNormalToPCAPlane(): "
	      << "PCAPlane not calculated or dimension not equal to three." << std::endl;
    exit(EXIT_FAILURE);
  }
  size_t dimension = E_[0].size();
  
  std::vector<double> N(dimension);
  
  N[0] = E_[0][1] * E_[1][2] - E_[0][2] * E_[1][1];
  N[1] = E_[0][2] * E_[1][0] - E_[0][0] * E_[1][2];
  N[2] = E_[0][0] * E_[1][1] - E_[0][1] * E_[1][0];
  
  return N;
}

std::vector<double> Cell::getNormalTriangular(DataMatrix &vertexData)
{
  size_t dimension = vertexData[0].size();
  if (numVertex()!=3 || dimension!=3) {
    std::cerr << "Cell::getNormalTriangular(): "
	      << "Only works for triangular cells in three dimensions." << std::endl;
    exit(EXIT_FAILURE);
  }
  
  std::vector<double> N(dimension);
  DataMatrix E(2);
  E[0].resize(dimension);
  E[1].resize(dimension);
  size_t v_0 = vertex(0)->index();
  size_t v_1 = vertex(1)->index();
  size_t v_2 = vertex(2)->index();
  for (size_t d=0; d<dimension; ++d) {
    E[0][d] = vertexData[v_1][d]-vertexData[v_0][d]; 
    E[1][d] = vertexData[v_2][d]-vertexData[v_1][d]; 
  }  
  N[0] = E[0][1] * E[1][2] - E[0][2] * E[1][1];
  N[1] = E[0][2] * E[1][0] - E[0][0] * E[1][2];
  N[2] = E[0][0] * E[1][1] - E[0][1] * E[1][0];
  
  return N;
}

int Cell::
vectorSignFromSort(std::vector<double> &n,
		   DataMatrix &vertexData) 
{
  size_t numV=numVertex();
  size_t dimension=vertexData[0].size();
  assert(numV>2);
  int sign=1;
  std::vector<int> scalarProdSign(numV);
  std::vector<double> scalarProdVal(numV);
  double scalarProdSum=0.0;
  for (size_t k=0; k<numV; ++k) {
    size_t k2=(k+1)%numV;
    size_t k3=(k+2)%numV;
    //Make sure ((v2-v1)x(v3-v2))n has same sign for all cells
    std::vector<double> nw1(dimension),nw2(dimension);
    for (size_t d=0; d<dimension; ++d) {
      nw1[d] = vertexData[vertex(k2)->index()][d]-vertexData[vertex(k)->index()][d];
      nw2[d] = vertexData[vertex(k3)->index()][d]-vertexData[vertex(k2)->index()][d];
    }
    //cross product
    double scalarProd=0.0;
    for (size_t d1=0; d1<dimension; ++d1) {
      size_t d2=(d1+1)%dimension;
      size_t d3=(d1+2)%dimension;
      scalarProd += (nw1[d1]*nw2[d2]-nw1[d2]*nw2[d1])*n[d3];
    }
    scalarProdVal[k] = scalarProd;
    scalarProdSum += scalarProd;
    if (scalarProd>0.0)
      scalarProdSign[k]=1;
    else
      scalarProdSign[k]=-1;
  }
  int scalarProdSignSum=0;
  for (size_t k=0; k<scalarProdSign.size(); ++k)
    scalarProdSignSum += scalarProdSign[k];
  
  if (scalarProdSignSum<0) {
    sign=-1;
  }
  else if (scalarProdSignSum==0) {
    //std::cerr << "Cell " << n << " has no majority sign in right hand rule expression." 
    //				<< std::endl;
    if (std::fabs(scalarProdSum)>0.01) {
      if (scalarProdSum<0.0) {
	sign=-1;
      }
    }
    else {
      std::cerr << "Cell::vectorSignFromSort() Cell " 
		<< index() << " has no majority sign in right hand rule expression." 
		<< std::endl;
      exit(-1);
    }
  }
  return sign;
}


bool Cell::isConcave(DataMatrix &vertexData, const double tolerance)
{
  signed int s = 0;
  
  for (size_t k = 0; k < numVertex(); ++k)
    {
      Vertex *vertex1 = vertex((k + 0) % numVertex());
      Vertex *vertex2 = vertex((k + 1) % numVertex());
      Vertex *vertex3 = vertex((k + 2) % numVertex());
      
      const size_t vertexIndex1 = vertex1->index();
      const size_t vertexIndex2 = vertex2->index();
      const size_t vertexIndex3 = vertex3->index();
      
      const double ux = vertexData[vertexIndex2][0] - vertexData[vertexIndex1][0];
      const double uy = vertexData[vertexIndex2][1] - vertexData[vertexIndex1][1];
      
      const double vx = vertexData[vertexIndex3][0] - vertexData[vertexIndex2][0];
      const double vy = vertexData[vertexIndex3][1] - vertexData[vertexIndex2][1];
      
      const double u = std::sqrt(ux * ux + uy * uy);
      const double v = std::sqrt(vx * vx + vy * vy);
      
      const double uv = ux *vx + uy * vy;
      
      const double costheta = uv / (u * v);
      
      if (std::abs(costheta - 1.0) < tolerance)
	{
	  continue;
	}
      
      const double tmp = myMath::sign(ux * vy - uy * vx);
      
      if (s == 0)
	{
	  s = tmp;
	}
      else if (s != tmp)
	{
	  return true;
	}
    }
  
  return false;
}

bool Cell::isFolded(DataMatrix &vertexData)
{
  for (size_t k = 0; k < numWall(); ++k)
    {
      const Wall *wall1 = wall(k);
      
      const size_t a_index = wall1->vertex1()->index();
      const size_t b_index = wall1->vertex2()->index();
      
      const double vx = vertexData[b_index][0] - vertexData[a_index][0];
      const double vy = vertexData[b_index][1] - vertexData[a_index][1];
      
      for (size_t l = k + 2; l < numWall(); ++l)
	{
	  if (k == 0 && l == numWall() - 1)
	    {
	      continue;
	    }
	  
	  const Wall *wall2 = wall(l);
	  
	  const size_t c_index = wall2->vertex1()->index();
	  const size_t d_index = wall2->vertex2()->index();
	  
	  const double ux = vertexData[d_index][0] - vertexData[c_index][0];
	  const double uy = vertexData[d_index][1] - vertexData[c_index][1];
	  
	  const double wx = vertexData[c_index][0] - vertexData[a_index][0];
	  const double wy = vertexData[c_index][1] - vertexData[a_index][1];
	  
	  const double detA = vx * uy - ux * vy;
	  
	  if (detA == 0.0)
	    {
	      continue;
	    }
	  
	  const double p = (uy * wx - ux * wy) / detA;
	  const double q = (vy * wx - vx * wy) / detA;
	  
	  if (p > 0.0 && p < 1.0 && q > 0.0 && q < 1.0)
	    {
	      return true;
	    }
	}
      
    }
  
  return false;
}
