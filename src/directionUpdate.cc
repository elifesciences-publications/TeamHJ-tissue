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
  if( paraValue.size()!=1 ) {
    std::cerr << "StrainDirection::"
							<< "StrainDirection() "
							<< "One parameter used flag_perpendicular."
							<< std::endl;
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
  tmp[0] = "flag_perpendicular";
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
			 std::vector< std::vector<double> > &vertexDerivs ) 
{
  size_t dimension = vertexData[0].size();
	assert( dimension==2 );
	
	//
	//Calculate strain directions and print walls and strain vectors
	//by using x,x+dt*dx/dt as two points
	//
	T.derivs(cellData,wallData,vertexData,
					 cellDerivs,wallDerivs,vertexDerivs);
	
	//
	// Update all cells
	//
	for (size_t cellI=0; cellI<T.numCell(); ++cellI) {
		//Create temporary x,y,dx positions
		size_t numV = T.cell(cellI).numVertex(); 
		std::vector< std::vector<double> > x(numV),y(numV),dx(numV),
			xM(numV),yM(numV),dxM(numV);
	
		double dt=1.0;
// 		std::vector<double> xMean(dimension),yMean(dimension),
// 			dxMean(dimension);
		std::vector<double> xMean(dimension),yMean(dimension);
	
		for( size_t i=0 ; i<numV ; ++i ) {
			size_t vI = T.cell(cellI).vertex(i)->index();
			x[i] = vertexData[vI];
			dx[i] = vertexDerivs[vI];
			std::vector<double> tmp(dimension);
			tmp[0] = x[i][0]+dt*dx[i][0];
			tmp[1] = x[i][1]+dt*dx[i][1];
			y[i] = tmp;
			xMean[0] += x[i][0];
			xMean[1] += x[i][1];
			yMean[0] += y[i][0];
			yMean[1] += y[i][1];
// 			dxMean[0] += dx[i][0];			
// 			dxMean[1] += dx[i][1];			
		}
		xMean[0] /=numV;
		xMean[1] /=numV;
		yMean[0] /=numV;
		yMean[1] /=numV;
// 		dxMean[0] /=numV;
// 		dxMean[1] /=numV;
		for( size_t i=0 ; i<numV ; ++i ) {
			xM[i].resize(dimension);
			xM[i][0] =x[i][0]-xMean[0];
			xM[i][1] =x[i][1]-xMean[1];
			yM[i].resize(dimension);
			yM[i][0] =y[i][0]-yMean[0];
			yM[i][1] =y[i][1]-yMean[1];
			dxM[i].resize(dimension);
// 			dxM[i][0] =dx[i][0]-dxMean[0];
// 			dxM[i][1] =dx[i][1]-dxMean[1];
			dxM[i][0] =dx[i][0];
			dxM[i][1] =dx[i][1];

		}
		
		//Calculate A = (x^t x)^{-1} (x^t y)
		std::vector<std::vector<double> > xTx(dimension),xTy(dimension),
			xTxM(dimension),A(dimension);
		for( size_t i=0 ; i<dimension ; ++i ) {
			xTx[i].resize(dimension);
			xTy[i].resize(dimension);
			xTxM[i].resize(dimension);
			A[i].resize(dimension);
		}
		
		for( size_t i=0 ; i<dimension ; ++i ) {
			for( size_t j=0 ; j<dimension ; ++j ) {
				for( size_t v=0 ; v<numV ; ++v ) {
					xTx[i][j] += xM[v][i]*xM[v][j];
					xTy[i][j] += xM[v][i]*yM[v][j];
				}
			}
		}
		double detM = xTx[0][0]*xTx[1][1]-xTx[0][1]*xTx[1][0];
		detM = 1.0/detM;
		xTxM[0][0] = detM*xTx[1][1];
		xTxM[1][1] = detM*xTx[0][0];
// 		xTxM[0][1] = -detM*xTx[1][0];
// 		xTxM[1][0] = -detM*xTx[0][1];
		xTxM[0][1] = -detM*xTx[0][1];
		xTxM[1][0] = -detM*xTx[1][0];
	
		//Calculate A
		A[0][0] = xTxM[0][0]*xTy[0][0] + xTxM[0][1]*xTy[1][0];
		A[0][1] = xTxM[0][0]*xTy[0][1] + xTxM[0][1]*xTy[1][1];
		A[1][0] = xTxM[1][0]*xTy[0][0] + xTxM[1][1]*xTy[1][0];
		A[1][1] = xTxM[1][0]*xTy[0][1] + xTxM[1][1]*xTy[1][1];

		//Apply SVD to A
		//
		
		//Make sure determinant is non-zero
		double detA = A[0][0]*A[1][1] - A[0][1]*A[1][0];
		if( detA==0 ) {
			std::cerr << "StrainDirection::update() Determinant zero\n";
			exit(-1);
		}
		double tau = std::atan2( A[0][0]-A[1][1],A[0][1]+A[1][0] );
		double omega = std::atan2( A[0][0]+A[1][1],A[0][1]-A[1][0] );
		double theta = 0.5*(tau-omega);
		
		//Create direction for update
		std::vector<double> n(dimension);
		double v = theta;
		if( parameter(0)==1.0 )		
			v = theta - 0.5 * M_PI;

// 		double a = A[0][0];
// 		double b = A[0][1];
// 		double c = A[1][0];
// 		double d = A[1][1];

// 		double t = std::sqrt((a + d) * (a + d) + (b - c) * (b - c));
// 		double w = std::sqrt((a - d) * (a - d) + (b + c) * (b + c));

// 		double p = 0.5 * (t + w);
// 		double q = 0.5 * (t - w);

// 		assert(std::abs(p) >= std::abs(q));
		n[0]=std::cos(v);
		n[1]=std::sin(v);
		for (size_t dim=0; dim<dimension; ++dim) 
			cellData[cellI][variableIndex(0,0)+dim] = n[dim];
	}
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

ForceDirection::
ForceDirection(std::vector<double> &paraValue, std::vector< std::vector<size_t> > &indValue)
{
	if (paraValue.size() != 1) {
		std::cerr << "ForceDirection::ForceDirection() " 
							<< "One parameter is used orientation_flag (0 for direction parallel with "
							<< "force, 1 for direction perpendicular to force)" << std::endl;
		exit(EXIT_FAILURE);
	}

	if (indValue.size() != 2 || indValue[0].size() != 1) {
		std::cerr << "ForceDirection::ForceDirection() \n"
							<< "First level: Start of cell direction index are used.\n"
							<< "Second level: Wall force indices\n";
		exit(EXIT_FAILURE);
	}

	setId("ForceDirection");
	setParameter(paraValue);  
	setVariableIndex(indValue);
	
	std::vector<std::string> tmp(numParameter());
	tmp.resize(numParameter());
	tmp[0] = "orientation_flag";
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
		
		if (cellData[cell.index()][variableIndex(0, 0) + 2] == 0) {
			continue;
		}
		double enumerator = 0.0;
		double denominator = 0.0;

		for (size_t i = 0; i < cell.numWall(); ++i) {
			Wall *wall = cell.wall(i);

			double wx = wall->vertex1()->position(0) - wall->vertex2()->position(0);
			double wy = wall->vertex1()->position(1) - wall->vertex2()->position(1);

			// Dodgy error check. Might have to improve it.
			if (wx < 0) {
				wx *= -1.0;
				wy *= -1.0;
			}
			double sigma = std::atan2(wy, wx);

			double c = std::cos(2.0 * sigma);
			double s = std::sin(2.0 * sigma);

			double force = 0.0;
			for (size_t j = 0; j < numVariableIndex(1); ++j) {
				force += wallData[wall->index()][variableIndex(1, j)];
			}
			enumerator += force * s;
			denominator += force * c;
		}
		
		double angle = std::atan2(enumerator, denominator);
		
		double x = std::cos(0.5 * angle);
		double y = std::sin(0.5 * angle);

		if (parameter(0) == 0) {
			cellData[cell.index()][variableIndex(0, 0) + 0] = x;
			cellData[cell.index()][variableIndex(0, 0) + 1] = y;
		} else {
			cellData[cell.index()][variableIndex(0, 0) + 0] = - y;
			cellData[cell.index()][variableIndex(0, 0) + 1] = x;
		}
	}
}

StretchDirection::
StretchDirection(std::vector<double> &paraValue, std::vector< std::vector<size_t> > &indValue)
{
	if (paraValue.size() != 1) {
		std::cerr << "StretchDirection::StretchDirection() " 
							<< "One parameter is used orientation_flag (0 for direction parallel with "
							<< "stretch, 1 for direction perpendicular to stretch)" << std::endl;
		exit(EXIT_FAILURE);
	}
	
	if (indValue.size() != 2 || indValue[0].size() != 1 || 
			indValue[1].size() != 1) {
		std::cerr << "StretchDirection::StretchDirection() \n"
							<< "First level: Start of cell direction index is used.\n"
							<< "Second level: Wall length index." << std::endl;
		exit(EXIT_FAILURE);
	}
	
	setId("StretchDirection");
	setParameter(paraValue);  
	setVariableIndex(indValue);
	
	std::vector<std::string> tmp(numParameter());
	tmp.resize(numParameter());
	tmp[0] = "orientation_flag";
	setParameterId(tmp);
}

void StretchDirection::
initiate(Tissue &T,
				 std::vector< std::vector<double> > &cellData,
				 std::vector< std::vector<double> > &wallData,
				 std::vector< std::vector<double> > &vertexData,
				 std::vector< std::vector<double> > &cellDerivs,
				 std::vector< std::vector<double> > &wallDerivs,
				 std::vector< std::vector<double> > &vertexDerivs)
{
	// No initialization
}

void StretchDirection::
update(Tissue &T, double h,
			 std::vector< std::vector<double> > &cellData,
			 std::vector< std::vector<double> > &wallData,
			 std::vector< std::vector<double> > &vertexData,
			 std::vector< std::vector<double> > &cellDerivs,
			 std::vector< std::vector<double> > &wallDerivs,
			 std::vector< std::vector<double> > &vertexDerivs)
{
	size_t dimension = vertexData[0].size();
	if (dimension==2) { 
		for (size_t n = 0; n < T.numCell(); ++n) {
			Cell cell = T.cell(n);

			if (cellData[cell.index()][variableIndex(0, 0) + 2] == 0) {
				continue;
			}

			double enumerator = 0.0;
			double denominator = 0.0;
			
			for (size_t i = 0; i < cell.numWall(); ++i) {
				Wall *wall = cell.wall(i);

				double wx = wall->vertex1()->position(0) - wall->vertex2()->position(0);
				double wy = wall->vertex1()->position(1) - wall->vertex2()->position(1);


				// Dodgy error check. Might have to improve it.
				if (wx < 0) {
					wx *= -1.0;
					wy *= -1.0;
				}
				double sigma = std::atan2(wy, wx);
				
				double c = std::cos(2.0 * sigma);
				double s = std::sin(2.0 * sigma);
				
				enumerator += s;
				denominator += c; 
			}
			
			double angle = std::atan2(enumerator, denominator);
			
			double x = std::cos(0.5 * angle);
			double y = std::sin(0.5 * angle);
		
			if (parameter(0) == 0) {
				cellData[cell.index()][variableIndex(0, 0) + 0] = x;
				cellData[cell.index()][variableIndex(0, 0) + 1] = y;
			} else {
				cellData[cell.index()][variableIndex(0, 0) + 0] = - y;
				cellData[cell.index()][variableIndex(0, 0) + 1] = x;
			}
		}
	}
	else if (dimension==3) {
		if (parameter(0) != 0) {
			std::cerr << "StretchDirection::update() Not yet implemented for three dimensions and"
								<< " perpendicular direction." << std::endl;
			exit(-1);
		}
		for (size_t n = 0; n < T.numCell(); ++n) {
			Cell cell = T.cell(n);
			double x = 0.0;
			double y = 0.0;
			double z = 0.0;

			if (cellData[cell.index()][variableIndex(0, 0) + dimension] == 0)
				continue;
			
			for (size_t i = 0; i < cell.numWall(); ++i) {
				Wall *wall = cell.wall(i);
				double wx = wall->vertex1()->position(0) - wall->vertex2()->position(0);
				double wy = wall->vertex1()->position(1) - wall->vertex2()->position(1);
				double wz = wall->vertex1()->position(2) - wall->vertex2()->position(2);
				double Aw = std::sqrt(wx * wx  + wy * wy + wz * wz);
				if (wx > 0) {
					wx /= Aw;
					wy /= Aw;
					wz /= Aw;
				} else {
					wx /= -Aw;
					wy /= -Aw;
					wz /= -Aw;
				}
				double restLength =wallData[wall->index()][variableIndex(1,0)];
				double stretch = (Aw-restLength)/restLength;
				
				x += wx * stretch;
				y += wy * stretch;
				z += wz * stretch;
			}
			double A = std::sqrt(x * x + y * y + z * z);
			
			// If the stretch is of zero magnitude, leave it as it is.
			if (A == 0.0)
				continue; 
			
			cellData[cell.index()][variableIndex(0, 0)] = x / A;
			cellData[cell.index()][variableIndex(0, 0) + 1] = y / A;
			cellData[cell.index()][variableIndex(0, 0) + 2] = z / A;
		}		
	}
}
