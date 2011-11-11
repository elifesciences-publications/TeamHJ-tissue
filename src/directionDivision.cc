//
// Filename     : directionDivision.cc
// Description  : Classes describing direction updates
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : June 2007
// Revision     : $Id:$
//
#include"directionDivision.h"
#include"baseDirectionDivision.h"
#include "myRandom.h"

//!Constructor
ParallellDirection::
ParallellDirection(std::vector<double> &paraValue, 
								std::vector< std::vector<size_t> > 
								&indValue ) 
{  
	//
  // Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=0 ) {
    std::cerr << "ParallellDirection::"
							<< "ParallellDirection() "
							<< "No parameters used.\n";
    exit(0);
  }
  if( indValue.size() != 1 || indValue[0].size() != 1 ) {
    std::cerr << "ParallellDirection::"
	      << "ParallellDirection() "
	      << "One variable index is used (start of cell direction).\n";
    exit(0);
  }
  //Set the variable values
  //////////////////////////////////////////////////////////////////////
  setId("ParallellDirection");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //////////////////////////////////////////////////////////////////////
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  //tmp[0] = "k_growth";
  setParameterId( tmp );
}

void ParallellDirection::
update(Tissue &T,size_t cellI,
			 DataMatrix &cellData,
			 DataMatrix &wallData,
			 DataMatrix &vertexData,
			 DataMatrix &cellDerivs,
			 DataMatrix &wallDerivs,
			 DataMatrix &vertexDerivs ) {
  
	// Nothing to do with the direction, but if direction connected to wall
	// the directionalWall vector needs to be extended and updated
	if( T.numDirectionalWall() == cellData.size()-1 ) {
		size_t cellNI=T.numDirectionalWall();
		size_t dimension = vertexData[0].size();
		std::vector<double> tmpN(dimension);
		T.setNumDirectionalWall(cellData.size());
		T.setDirectionalWall(cellNI,T.directionalWall(cellI));
		if( T.directionalWall(cellI)<T.cell(cellI).numWall() ) {
			// Find wall with direction closest to direction stored in cellData
			// For both cellI and cellNI
			std::vector<size_t> cellList(2);
			cellList[0]=cellI;
			cellList[1]=cellNI;
			for( size_t ii=0; ii<cellList.size(); ++ii ) {
				size_t i=cellList[ii];
				if ( cellData[i][variableIndex(0,0)+dimension] > 0.0 ) {
					double normW = 0.0;
					size_t v1I = T.cell(i).wall(0)->vertex1()->index();
					size_t v2I = T.cell(i).wall(0)->vertex2()->index();
					for (size_t dim=0; dim<dimension; ++dim) {
						tmpN[dim] = vertexData[v1I][dim] - vertexData[v2I][dim];
						normW += tmpN[dim]*tmpN[dim];
					}
					normW = std::sqrt( normW );
					if (normW<=0.0) {
						std::cerr << "ParallellDirection::update() Normalization=0!"
											<< std::endl;
						exit(-1);
					}
					normW = 1.0/normW;
					double prod=0.0;
					for (size_t dim=0; dim<dimension; ++dim) {
						tmpN[dim] *= normW;
						prod += tmpN[dim]*cellData[i][dim+variableIndex(0,0)];
					}
					size_t maxK=0;
					double maxProd = std::fabs(prod);
				
					for (size_t k=1; k<T.cell(i).numWall(); ++k) {
						normW = 0.0;
						v1I = T.cell(i).wall(k)->vertex1()->index();
						v2I = T.cell(i).wall(k)->vertex2()->index();
						for (size_t dim=0; dim<dimension; ++dim) {
							tmpN[dim] = vertexData[v1I][dim] - vertexData[v2I][dim];
							normW += tmpN[dim]*tmpN[dim];
						}
						normW = std::sqrt( normW );
						if (normW<=0.0) {
							std::cerr << "ParallellDirection::update() Normalization=0!"
												<< std::endl;
							exit(-1);
						}
						normW = 1.0/normW;
						prod=0.0;
						for (size_t dim=0; dim<dimension; ++dim) {
							tmpN[dim] *= normW;
							prod += tmpN[dim]*cellData[i][dim+variableIndex(0,0)];
						}
						prod = std::fabs(prod);
						if( prod>maxProd ) {
							maxProd=prod;
							maxK=k;
						}
					}
					T.setDirectionalWall(i,maxK);
				}
			}
		}
	}
	else if( T.numDirectionalWall() ) {
		std::cerr << "ParallellDirection::update() Strange number of directions."
							<< std::endl;
		exit(-1);
	}
}

//!Constructor
PerpendicularDirection::
PerpendicularDirection(std::vector<double> &paraValue, 
											 std::vector< std::vector<size_t> > 
											 &indValue ) 
{  
	//
  // Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=0 ) {
    std::cerr << "PerpendicularDirection::"
							<< "PerpendicularDirection() "
							<< "No parameters used.\n";
    exit(0);
  }
  if( indValue.size() != 1 || indValue[0].size() != 1 ) {
    std::cerr << "PerpendicularDirection::"
							<< "PerpendicularDirection() "
							<< "One variable index is used (start of cell direction).\n";
    exit(0);
  }
  //Set the variable values
  //////////////////////////////////////////////////////////////////////
  setId("PerpendicularDirection");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //////////////////////////////////////////////////////////////////////
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  //tmp[0] = "k_growth";
  setParameterId( tmp );
}

void PerpendicularDirection::
update(Tissue &T,size_t cellI,
			 DataMatrix &cellData,
			 DataMatrix &wallData,
			 DataMatrix &vertexData,
			 DataMatrix &cellDerivs,
			 DataMatrix &wallDerivs,
			 DataMatrix &vertexDerivs ) {

	// Extract the perpendicular direction
	size_t cellNI=cellData.size()-1;
	size_t dimension = vertexData[0].size();
	assert( dimension==2 );
	std::vector<double> tmpDirection(dimension);
	tmpDirection[0] = cellData[cellI][variableIndex(0,0)+1];
	tmpDirection[1] = -cellData[cellI][variableIndex(0,0)];
	cellData[cellI][variableIndex(0,0)] = cellData[cellNI][variableIndex(0,0)]
		= tmpDirection[0];
	cellData[cellI][variableIndex(0,0)+1] = cellData[cellNI][variableIndex(0,0)+1]
		= tmpDirection[1];

	// If direction connected to wall
	// the directionalWall vector needs to be extended and updated
	if( T.numDirectionalWall() == cellData.size()-1 ) {
		std::vector<double> tmpN(dimension);
		T.setNumDirectionalWall(cellData.size());
		T.setDirectionalWall(cellNI,T.directionalWall(cellI));
		if( T.directionalWall(cellI)<T.cell(cellI).numWall() ) {

			// Find wall with direction closest to direction stored in cellData
			// For both cellI and cellNI
			std::vector<size_t> cellList(2);
			cellList[0]=cellI;
			cellList[1]=cellNI;
			for( size_t ii=0; ii<cellList.size(); ++ii ) {
				size_t i=cellList[ii];
				if ( cellData[i][variableIndex(0,0)+dimension] > 0.0 ) {
					double normW = 0.0;
					size_t v1I = T.cell(i).wall(0)->vertex1()->index();
					size_t v2I = T.cell(i).wall(0)->vertex2()->index();
					for (size_t dim=0; dim<dimension; ++dim) {
						tmpN[dim] = vertexData[v1I][dim] - vertexData[v2I][dim];
						normW += tmpN[dim]*tmpN[dim];
					}
					normW = std::sqrt( normW );
					if (normW<=0.0) {
						std::cerr << "PerpendicularDirection::update() Normalization=0!"
											<< std::endl;
						exit(-1);
					}
					normW = 1.0/normW;
					double prod=0.0;
					for (size_t dim=0; dim<dimension; ++dim) {
						tmpN[dim] *= normW;
						prod += tmpN[dim]*cellData[i][dim+variableIndex(0,0)];
					}
					size_t maxK=0;
					double maxProd = std::fabs(prod);
					
					for (size_t k=1; k<T.cell(i).numWall(); ++k) {
						normW = 0.0;
						v1I = T.cell(i).wall(k)->vertex1()->index();
						v2I = T.cell(i).wall(k)->vertex2()->index();
						for (size_t dim=0; dim<dimension; ++dim) {
							tmpN[dim] = vertexData[v1I][dim] - vertexData[v2I][dim];
							normW += tmpN[dim]*tmpN[dim];
						}
						normW = std::sqrt( normW );
						if (normW<=0.0) {
							std::cerr << "PerpendicularDirection::update() Normalization=0!"
												<< std::endl;
							exit(-1);
						}
						normW = 1.0/normW;
						prod=0.0;
						for (size_t dim=0; dim<dimension; ++dim) {
							tmpN[dim] *= normW;
							prod += tmpN[dim]*cellData[i][dim+variableIndex(0,0)];
						}
						prod = std::fabs(prod);
						if( prod>maxProd ) {
							maxProd=prod;
							maxK=k;
						}
					}
					T.setDirectionalWall(i,maxK);
				}
			}
		}
	}
	else if( T.numDirectionalWall() ) {
		std::cerr << "PerpendicularDirection::update() Strange number of directions."
							<< std::endl;
		exit(-1);
	}
}

RandomDirection::RandomDirection(std::vector<double> &paraValue, std::vector< std::vector<size_t> > &indValue) 
{  
	//
	// Do some checks on the parameters and variable indeces
	//
	if (paraValue.size() != 0) {
		std::cerr << "RandomDirection::RandomDirection() No parameters used.\n";
		std::exit(EXIT_FAILURE);
	}

	if (indValue.size() != 1 || indValue[0].size() != 1) {
		std::cerr << "RandomDirection::RandomDirection() "
		<< "One variable index is used (start of cell direction).\n";
		std::exit(EXIT_FAILURE);
	}

	//Set the variable values
	//////////////////////////////////////////////////////////////////////
	setId("RandomDirection");
	setParameter(paraValue);  
	setVariableIndex(indValue);
	
	//Set the parameter identities
	//////////////////////////////////////////////////////////////////////
	std::vector<std::string> tmp(numParameter());
	tmp.resize( numParameter() );
	setParameterId(tmp);
}

void RandomDirection::update(Tissue &T, size_t cellIndex,
	DataMatrix &cellData,
	DataMatrix &wallData,
	DataMatrix &vertexData,
	DataMatrix &cellDerivs,
	DataMatrix &wallDerivs,
	DataMatrix &vertexDerivs)
{
	size_t dimension = vertexData[0].size();

	if (dimension != 2) {
		std::cerr << "RandomDirection only support two dimensions.\n";
		std::exit(EXIT_FAILURE);
	}

	const size_t xIndex = variableIndex(0, 0) + 0;
	const size_t yIndex = variableIndex(0, 0) + 1;

	const double angle = 2.0 * M_PI * myRandom::Rnd();

	const double x = std::cos(angle);
	const double y = std::sin(angle);

	cellData[cellIndex][xIndex] = x;
	cellData[cellIndex][yIndex] = y;
}



