//
// Filename     : directionUpdate.cc
// Description  : Classes describing direction updates
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : June 2007
// Revision     : $Id:$
//
#include"directionUpdate.h"
#include"baseDirectionUpdate.h"

//!Constructor
StaticDirection::
StaticDirection(std::vector<double> &paraValue, 
								std::vector< std::vector<size_t> > 
								&indValue ) 
{  
  // Do some checks on the parameters and variable indeces
  if( paraValue.size()!=0 ) {
    std::cerr << "StaticDirection::"
							<< "StaticDirection() "
							<< "No parameters used.\n";
    exit(0);
  }
  if( indValue.size() != 0 ) {
    std::cerr << "StaticDirection::"
							<< "StaticDirection() "
							<< "No variable index is used." << std::endl;
    exit(0);
  }

  // Set the variable values
	setId("StaticDirection");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  // Set the parameter identities
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  //tmp[0] = "k_growth";
  setParameterId( tmp );
}

void StaticDirection::
initiate(Tissue &T,
				 std::vector< std::vector<double> > &cellData,
				 std::vector< std::vector<double> > &wallData,
				 std::vector< std::vector<double> > &vertexData,
				 std::vector< std::vector<double> > &cellDerivs,
				 std::vector< std::vector<double> > &wallDerivs,
				 std::vector< std::vector<double> > &vertexDerivs ) {
	
}

void StaticDirection::
update(Tissue &T, double h,
			 std::vector< std::vector<double> > &cellData,
			 std::vector< std::vector<double> > &wallData,
			 std::vector< std::vector<double> > &vertexData,
			 std::vector< std::vector<double> > &cellDerivs,
			 std::vector< std::vector<double> > &wallDerivs,
			 std::vector< std::vector<double> > &vertexDerivs ) {
  
}

//!Constructor
WallDirection::
WallDirection(std::vector<double> &paraValue, 
							std::vector< std::vector<size_t> > 
							&indValue ) 
{  
  // Do some checks on the parameters and variable indeces
  if( paraValue.size()!=0 ) {
    std::cerr << "WallDirection::"
							<< "WallDirection() "
							<< "No parameters used.\n";
    exit(0);
  }
  if( indValue.size() != 1 || indValue[0].size() != 1 ) {
    std::cerr << "WallDirection::"
							<< "WallDirection() "
							<< "One variable index is used (start of cell direction).\n";
    exit(0);
  }

  // Set the variable values
	setId("WallDirection");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  // Set the parameter identities
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  //tmp[0] = "k_growth";
  setParameterId( tmp );
}

void WallDirection::
initiate(Tissue &T,
				 std::vector< std::vector<double> > &cellData,
				 std::vector< std::vector<double> > &wallData,
				 std::vector< std::vector<double> > &vertexData,
				 std::vector< std::vector<double> > &cellDerivs,
				 std::vector< std::vector<double> > &wallDerivs,
				 std::vector< std::vector<double> > &vertexDerivs ) {

	// Find walls with direction closest to given direction
	size_t dimension = T.numDimension();
	T.setNumDirectionalWall(T.numCell());
	std::vector<double> tmpN(dimension);
	
	for (size_t i=0; i<T.numCell(); ++i) {
		if ( T.cell(i).variable(variableIndex(0,0)+dimension) > 0.0 ) {
			double normW = 0.0;
			for (size_t dim=0; dim<dimension; ++dim) {
				tmpN[dim] = T.cell(i).wall(0)->vertex1()->position(dim) -
					T.cell(i).wall(0)->vertex2()->position(dim);
				normW += tmpN[dim]*tmpN[dim];
			}
			normW = std::sqrt( normW );
			if (normW<=0.0) {
				std::cerr << "VertexFromWallSpringMT::initiate Normalization=0!"
									<< std::endl;
				exit(-1);
			}
			normW = 1.0/normW;
			double prod=0.0;
			for (size_t dim=0; dim<dimension; ++dim) {
				tmpN[dim] *= normW;
				prod += tmpN[dim]*T.cell(i).variable(dim+variableIndex(0,0));
			}
			size_t maxK=0;
			double maxProd = std::fabs(prod);
			
			for (size_t k=1; k<T.cell(i).numWall(); ++k) {
				normW = 0.0;
				for (size_t dim=0; dim<dimension; ++dim) {
					tmpN[dim] = T.cell(i).wall(k)->vertex1()->position(dim) -
						T.cell(i).wall(k)->vertex2()->position(dim);
					normW += tmpN[dim]*tmpN[dim];
				}
				normW = std::sqrt( normW );
				if (normW<=0.0) {
					std::cerr << "VertexFromWallSpringMT::initiate Normalization=0!"
										<< std::endl;
					exit(-1);
				}
				normW = 1.0/normW;
				prod=0.0;
				for (size_t dim=0; dim<dimension; ++dim) {
					tmpN[dim] *= normW;
					prod += tmpN[dim]*T.cell(i).variable(dim+variableIndex(0,0));
				}
				prod = std::fabs(prod);
				if( prod>maxProd ) {
					maxProd=prod;
					maxK=k;
				}
			}
			T.setDirectionalWall(i,maxK);
			assert( T.cell(i).numWall()>T.directionalWall(i) );
		}
		else
			T.setDirectionalWall(i,static_cast<size_t>(-1));
	}	  
}

void WallDirection::
update(Tissue &T, double h,
			 std::vector< std::vector<double> > &cellData,
			 std::vector< std::vector<double> > &wallData,
			 std::vector< std::vector<double> > &vertexData,
			 std::vector< std::vector<double> > &cellDerivs,
			 std::vector< std::vector<double> > &wallDerivs,
			 std::vector< std::vector<double> > &vertexDerivs ) {
  
	size_t dimension=vertexData[0].size(); 
	for (size_t i=0; i<T.numDirectionalWall(); ++i) {
		if (T.directionalWall(i)<T.cell(i).numWall()) {
			std::vector<double> tmpN(dimension);
			size_t v1I = T.cell(i).wall(T.directionalWall(i))->vertex1()->index();
			size_t v2I = T.cell(i).wall(T.directionalWall(i))->vertex2()->index();
			double normW=0.0;
			for (size_t dim=0; dim<dimension; ++dim) {
				tmpN[dim] = vertexData[v2I][dim] - vertexData[v1I][dim];
				normW += tmpN[dim]*tmpN[dim];
			}
			normW = std::sqrt(normW);
			if( normW<=0.0 ) {
				std::cerr << "WallDirection::update() Wrong norm factor"
									<< std::endl;
				exit(-1);
			}
			normW = 1.0/normW;
			for (size_t dim=0; dim<dimension; ++dim) {
				tmpN[dim] *= normW;
				cellData[i][dim+variableIndex(0,0)] = tmpN[dim];
			}
		}
	}
}

//!Constructor
StrainDirection::
StrainDirection(std::vector<double> &paraValue, 
							std::vector< std::vector<size_t> > 
							&indValue ) 
{  
	//
  // Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=0 ) {
    std::cerr << "StrainDirection::"
							<< "StrainDirection() "
							<< "No parameters used.\n";
    exit(0);
  }
  if( indValue.size() != 1 || indValue[0].size() != 1 ) {
    std::cerr << "StrainDirection::"
							<< "StrainDirection() "
							<< "One variable index is used (start of cell direction).\n";
    exit(0);
  }
  //Set the variable values
  //////////////////////////////////////////////////////////////////////
  setId("StrainDirection");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //////////////////////////////////////////////////////////////////////
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  //tmp[0] = "k_growth";
  setParameterId( tmp );
}

void StrainDirection::
initiate(Tissue &T,
				 std::vector< std::vector<double> > &cellData,
				 std::vector< std::vector<double> > &wallData,
				 std::vector< std::vector<double> > &vertexData,
				 std::vector< std::vector<double> > &cellDerivs,
				 std::vector< std::vector<double> > &wallDerivs,
				 std::vector< std::vector<double> > &vertexDerivs ) {
  
}

void StrainDirection::
update(Tissue &T, double h,
			 std::vector< std::vector<double> > &cellData,
			 std::vector< std::vector<double> > &wallData,
			 std::vector< std::vector<double> > &vertexData,
			 std::vector< std::vector<double> > &cellDerivs,
			 std::vector< std::vector<double> > &wallDerivs,
			 std::vector< std::vector<double> > &vertexDerivs ) {
  
} 

//!Constructor
GradientDirection::
GradientDirection(std::vector<double> &paraValue, 
							std::vector< std::vector<size_t> > 
							&indValue ) 
{  
	//
  // Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=0 ) {
    std::cerr << "GradientDirection::"
							<< "GradientDirection() "
							<< "No parameters used.\n";
    exit(0);
  }
  if( indValue.size() != 1 || indValue[0].size() != 2 ) {
    std::cerr << "GradientDirection::"
							<< "GradientDirection() "
							<< "Two variable indices are used (start of cell direction,"
							<< " gradient variable).\n";
    exit(0);
  }
  //Set the variable values
  //////////////////////////////////////////////////////////////////////
  setId("GradientDirection");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //////////////////////////////////////////////////////////////////////
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  //tmp[0] = "k_growth";
  setParameterId( tmp );
}

void GradientDirection::
initiate(Tissue &T,
				 std::vector< std::vector<double> > &cellData,
				 std::vector< std::vector<double> > &wallData,
				 std::vector< std::vector<double> > &vertexData,
				 std::vector< std::vector<double> > &cellDerivs,
				 std::vector< std::vector<double> > &wallDerivs,
				 std::vector< std::vector<double> > &vertexDerivs ) {
  
}

void GradientDirection::
update(Tissue &T, double h,
			 std::vector< std::vector<double> > &cellData,
			 std::vector< std::vector<double> > &wallData,
			 std::vector< std::vector<double> > &vertexData,
			 std::vector< std::vector<double> > &cellDerivs,
			 std::vector< std::vector<double> > &wallDerivs,
			 std::vector< std::vector<double> > &vertexDerivs ) {
  
}

ForceDirection::ForceDirection(std::vector<double> &paraValue, std::vector< std::vector<size_t> > &indValue)
{
	if (paraValue.size() != 0) {
		std::cerr << "ForceDirection::ForceDirection() " 
				<< "No parameters are used.\n";
		exit(EXIT_FAILURE);
	}

	if (indValue.size() != 1 || indValue[0].size() != 2) {
		std::cerr << "ForceDirection::ForceDirection() "
				<< "Wall force index and start of cell direction index are used.\n";
		exit(EXIT_FAILURE);
	}

	setId("ForceDirection");
	setParameter(paraValue);  
	setVariableIndex(indValue);
	
	std::vector<std::string> tmp(numParameter());
	tmp.resize(numParameter());
	setParameterId(tmp);

}
  
void ForceDirection::initiate(Tissue &T,
						std::vector< std::vector<double> > &cellData,
						std::vector< std::vector<double> > &wallData,
						std::vector< std::vector<double> > &vertexData,
						std::vector< std::vector<double> > &cellDerivs,
						std::vector< std::vector<double> > &wallDerivs,
						std::vector< std::vector<double> > &vertexDerivs)
{
	// No initialization
}

void ForceDirection::update(Tissue &T, double h,
					   std::vector< std::vector<double> > &cellData,
					   std::vector< std::vector<double> > &wallData,
					   std::vector< std::vector<double> > &vertexData,
					   std::vector< std::vector<double> > &cellDerivs,
					   std::vector< std::vector<double> > &wallDerivs,
					   std::vector< std::vector<double> > &vertexDerivs)
{
	for (size_t n = 0; n < T.numCell(); ++n) {
		Cell cell = T.cell(n);
		double x = 0.0;
		double y = 0.0;

		if (cellData[cell.index()][variableIndex(0, 1) + 2] == 0)
			continue;

		for (size_t i = 0; i < cell.numWall(); ++i) {
			Wall *wall = cell.wall(i);
			double wx = wall->vertex1()->position(0) - wall->vertex2()->position(0);
			double wy = wall->vertex1()->position(1) - wall->vertex2()->position(1);
			double Aw = std::sqrt(wx * wx  + wy * wy);
			if (wx > 0) {
				wx = wx / Aw;
				wy = wy / Aw;
			} else {
				wx = -1 * wx / Aw;
				wy = -1 * wy / Aw;
			}
			double force = wallData[wall->index()][variableIndex(0, 1)];

			x += wx * force;
			y += wy * force;
		}
		
		double A = std::sqrt(x * x + y * y);
		cellData[cell.index()][variableIndex(0, 1)] = x / A;
		cellData[cell.index()][variableIndex(0, 1) + 1] = y / A;
	}
}



