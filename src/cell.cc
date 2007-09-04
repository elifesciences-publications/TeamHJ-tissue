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
#include "math.h"

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

void Cell::sortWallAndVertexNew(Tissue &T) {
	
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
	size_t wallIndex=0,wallIndexStart=0;
	size_t vertexIndex=0;
	//Find initial wall and initial vertices
	size_t foundCellSortFlag=0;
	while (wallIndex<numWall() && !wall(wallIndex)->cellSort1())
		++wallIndex;
	if (wallIndex<numWall()) {
		foundCellSortFlag=1;
		wallIndexStart=vertexIndex=wallIndex;
	}
	else
		wallIndex=wallIndexStart=0;
	
	// Set sorting order for first wall
	if (!foundCellSortFlag) {
		tmpWall[wallIndex] = wall(wallIndex);
		tmpVertex[vertexIndex] = wall(wallIndex)->vertex1();
		++vertexIndex;
		tmpVertex[vertexIndex] = wall(wallIndex)->vertex2();
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
	else {
		tmpWall[wallIndex] = wall(wallIndex);
		if (wall(wallIndex)->cell1()==this) {
			if (wall(wallIndex)->cellSort1()==1) {
				tmpVertex[vertexIndex] = wall(wallIndex)->vertex1();
				++vertexIndex;
				vertexIndex = vertexIndex%numVertex();
				tmpVertex[vertexIndex] = wall(wallIndex)->vertex2();
			}
			else if (wall(wallIndex)->cellSort1()==-1) {
				tmpVertex[vertexIndex] = wall(wallIndex)->vertex2();
				++vertexIndex;
				vertexIndex = vertexIndex%numVertex();
				tmpVertex[vertexIndex] = wall(wallIndex)->vertex1();
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
				++vertexIndex;
				vertexIndex = vertexIndex%numVertex();
				tmpVertex[vertexIndex] = wall(wallIndex)->vertex2();
			}
			else if (wall(wallIndex)->cellSort2()==-1) {
				tmpVertex[vertexIndex] = wall(wallIndex)->vertex2();
				++vertexIndex;
				vertexIndex = vertexIndex%numVertex();
				tmpVertex[vertexIndex] = wall(wallIndex)->vertex1();
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
	
	//CONTINUE HERE!!!
	//size_t wallIndexMod = wallIndex%numWall();
	while( wallIndex<wallIndexStart+numWall()-1 ) {
		for( size_t wI=0 ; wI<numWall() ; ++wI ) {
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

//!Calculates the cell position from average vertex position
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
positionFromVertex(std::vector< std::vector<double> > &vertexData) 
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


void Cell::calculatePCAPlane(std::vector< std::vector<double> > &vertexData)
{
	size_t dimensions = vertexData[0].size();
	size_t numberOfVertices = vertex_.size();

	// Copy vertex data to temporary container and calculate mean values.
 
	std::vector< std::vector<double> > vertices(numberOfVertices, dimensions);
	std::vector<double> mean(dimensions, 0.0);

	for (size_t i = 0; i < numberOfVertices; ++i) {
		Vertex *v = vertex(i);
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
	
	std::vector< std::vector<double> > R(dimensions, dimensions);
	for (size_t k = 0; k < dimensions; ++k) {
		for (size_t l = 0; l < dimensions; ++l) {
			for (size_t i = 0; i < numberOfVertices; ++i) {
				R[k][l] += vertices[i][k] * vertices[i][l];
			}
			R[k][l] /= numberOfVertices;
		}
	}

	// Find the eigenvectors with the two greatests corresponding eigenvalues.

	std::vector< std::vector<double> > candidates;

	std::vector< std::vector<double> > V;
	std::vector<double> d;

	jacobiTransformation(R , V, d);
       
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

	E.resize(2);
	E[0].resize(dimensions);
	E[1].resize(dimensions);

	for (size_t i = 0; i < dimensions; ++i) {
		//		E[0][i] = V[max1][i];
		E[0][i] = V[i][max1];
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
		// 		E[1][i] = s * E[0][i] + V[max2][i];
		E[1][i] = s * E[0][i] + V[i][max2];
	}

	for (size_t i = 0; i < E.size(); ++i) {
		double sum = 0.0;
		for (size_t j = 0; j < dimensions; ++j) {
			sum += E[i][j] * E[i][j];
		}
		for (size_t j = 0; j < dimensions; ++j) {
			E[i][j] /= std::sqrt(sum);
		}
	}	
}

std::vector< std::vector<double> > Cell::getPCAPlane(void)
{
	return E;
}

std::vector< std::pair<double, double> > Cell::projectVerticesOnPCAPlane(std::vector< std::vector<double> > &vertexData)
{
	size_t dimensions = vertexData[0].size();
	size_t numberOfVertices = vertex_.size();	

	// Copy vertex data to temporary container and calculate mean values.
 
	std::vector< std::vector<double> > vertices(numberOfVertices, dimensions);
	std::vector<double> mean(dimensions, 0.0);
	
	for (size_t i = 0; i < numberOfVertices; ++i) {
		Vertex *v = vertex(i);
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
			tmp.first += E[0][j] * vertices[i][j];
			tmp.second += E[1][j] * vertices[i][j];
		}
		coordinates.push_back(tmp);
	}
	
	return coordinates;
}

std::vector<double> Cell::getNormalToPCAPlane(void)
{
	size_t dimension = E[0].size();
	if (dimension != 3) {
		std::cerr << "Cell::getNormalToPCAPlane(): Dimension is not equal to three." << std::endl;
		exit(EXIT_FAILURE);
	}
	
	std::vector<double> N(dimension);

	N[0] = E[0][1] * E[1][2] - E[0][2] * E[1][1];
	N[1] = E[0][2] * E[1][0] - E[0][0] * E[1][2];
	N[2] = E[0][0] * E[1][1] - E[0][1] * E[1][0];

	return N;
}
