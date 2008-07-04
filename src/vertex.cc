/**
 * Filename     : vertex.cc
 * Description  : A class describing a vertex
 * Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
 * Created      : April 2006
 * Revision     : $Id:$
 */

#include "vertex.h"
#include "wall.h"
#include "math.h"
#include <cmath>
#include <vector>

Vertex::Vertex() 
{
}

Vertex::Vertex( const Vertex & vertexCopy ) 
{
  index_ = vertexCopy.index();
  id_ = vertexCopy.id();
  position_ = vertexCopy.position();
  cell_ = vertexCopy.cell();
  wall_ = vertexCopy.wall();
  stressDirection_ = vertexCopy.getStressDirection();
}
  
Vertex::~Vertex() 
{
}

int Vertex::removeCell( Cell* val ) 
{
	for (size_t k=0; k<cell_.size(); ++k)
		if (cell_[k]==val) {
			cell_[k]=cell_[cell_.size()-1];
			cell_.pop_back();
			return 1;
		}
	return 0;
}

int Vertex::removeWall( Wall* val ) 
{
	for (size_t k=0; k<wall_.size(); ++k)
		if (wall_[k]==val) {
			wall_[k]=wall_[wall_.size()-1];
			wall_.pop_back();
			return 1;
		}
	return 0;
}

int Vertex::isBoundary(Cell *background) const
{
	for (size_t wI=0; wI<numWall(); ++wI)
		if (wall_[wI]->hasCell(background))
			return 1;
	return 0;
}

void Vertex::calculateStressDirection(std::vector< std::vector<double> > &vertexData,
	std::vector< std::vector<double> > &wallData, std::vector<size_t> wallForceIndexes)
{
	size_t dimensions = vertexData[0].size();
	size_t numberOfWalls = wall_.size();
	
	// Copy wall force data to temporary container and calculate mean values.
 
	std::vector< std::vector<double> > walls(2 * numberOfWalls, dimensions);

	for (size_t i = 0; i < wall_.size(); ++i) {
		Wall *w = wall_[i];
		Vertex *v1 = w->vertex1();
		Vertex *v2 = w->vertex2();

		double force = 0.0;
		for (size_t j = 0; j < wallForceIndexes.size(); ++j) {
			force += wallData[w->index()][wallForceIndexes[j]];
		}

		std::vector<double> n(dimensions);
		double A = 0.0;

		for (size_t j = 0; j < dimensions; ++j) {
			if (v1 == this) {
				n[j] = vertexData[v2->index()][j] - vertexData[v1->index()][j];
			} else {
				n[j] = vertexData[v1->index()][j] - vertexData[v2->index()][j];
			}
			A += n[j] * n[j];
		}

		A = std::sqrt(A);

		for (size_t j = 0; j < dimensions; ++j) {
			n[j] /= A;
		}
		
		for (size_t j = 0; j < dimensions; ++j) {
			walls[2 * i + 0][j] = force * n[j];
			walls[2 * i + 1][j] = -force * n[j];
		}
	}

	//
	// Mean is always equal to zero so no need to subtract it.
	//

 	// Calculate the correlation matrix.
	
 	std::vector< std::vector<double> > R(dimensions, dimensions);

	for (size_t i = 0; i < dimensions; ++i) {
		for (size_t j = 0; j < dimensions; ++j) {
			R[i][j] = 0.0;
		}
	}

 	for (size_t k = 0; k < dimensions; ++k) {
 		for (size_t l = 0; l < dimensions; ++l) {
 			for (size_t i = 0; i < numberOfWalls; ++i) {
 				R[k][l] += walls[i][k] * walls[i][l];
 			}
 			R[k][l] /= numberOfWalls;
 		}
 	}

 	// Find the eigenvectors with the two greatests corresponding eigenvalues.

 	std::vector< std::vector<double> > candidates;

 	std::vector< std::vector<double> > V;
 	std::vector<double> d;

 	jacobiTransformation(R , V, d);
       
 	double max = 0.0;
 	size_t max1 = d.size();

 	max = 0.0;
 	for (size_t i = 0; i < d.size(); ++i) {
 		if (std::abs(d[i]) >= max) {
 			max1 = i;
 			max = std::abs(d[i]);
 		}
 	}

 	if (max1 == d.size()) {
 		std::cerr << "Vertex::calculateStressDirection(): Unexpected behaviour." << std::endl;
 		exit(EXIT_FAILURE);
 	}

  	// Find orthonormal basis.
 
	std::vector< std::vector<double> > E_(1);
	E_[0].resize(dimensions);

 	for (size_t i = 0; i < dimensions; ++i) {
  		//		E_[0][i] = V[max1][i];
 		E_[0][i] = V[i][max1];
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

	stressDirection_.resize(dimensions);
	for (size_t i = 0; i < dimensions; ++i) {
		stressDirection_[i] = d[max1] * E_[0][i];
	}
}

std::vector<double> Vertex::getStressDirection(void) const
{
	return stressDirection_;
}
