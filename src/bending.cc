//
// Filename     : bending.cc
// Description  : Classes describing reactions related to bending moments
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : September 2013
// Revision     : $Id:$
//
#include<cmath>
#include"baseReaction.h"
#include"bending.h"
#include"tissue.h"

namespace Bending {
  
  NeighborCenter::
  NeighborCenter(std::vector<double> &paraValue, 
		 std::vector< std::vector<size_t> > 
		 &indValue )
  {
    //Do some checks on the parameters and variable indices
    //
    if( paraValue.size()!=1 ) {
      std::cerr << "Bending::NeighborCenter::"
		<< "NeighborCenter() "
		<< "One parameter, k_bend, should be provided." << std::endl;
      exit(EXIT_FAILURE);
    }
    if( indValue.size() != 1 || indValue[0].size() != 1 ) {
      std::cerr << "Bending::NeighborCenter::"
		<< "NeighborCenter() "
		<< "One index level with one index (wall length) given" 
		<< std::endl;
      exit(EXIT_FAILURE);
    }
    //Set the variable values
    //
    setId("Bending::NeighborCenter");
    setParameter(paraValue);  
    setVariableIndex(indValue);
  }
  
  void NeighborCenter::
  derivs(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs )
  {
    size_t numCells = T.numCell();
    size_t dimension = T.vertex(0).numPosition();
    size_t Li = variableIndex(0,0);
    for (size_t i=0; i<numCells; ++i) {
      size_t numWalls = T.cell(i).numWall();
      for (size_t k=0; k<numWalls; ++k) {
	size_t kp = k<numWalls-1 ? k+1 : 0;
	size_t km = k>0 ? k-1 : numWalls-1;

	// Get global vertex indices
	size_t j = T.cell(i).vertex(k)->index();
	size_t jp = T.cell(i).vertex(kp)->index();
	size_t jm = T.cell(i).vertex(km)->index();
	
	// Get edges
	size_t ep = T.cell(i).wall(k)->index();
	size_t em = T.cell(i).wall(km)->index();
	
	// Update for each dimension
	for (size_t d=0; d<dimension; ++d) {
	  double derivs = -parameter(0)*(vertexData[j][d] - 
					 (vertexData[jp][d]*wallData[em][Li] + 
					  vertexData[jm][d]*wallData[ep][Li])/
					 (wallData[em][Li]+wallData[ep][Li]) );
	  vertexDerivs[k][d] += derivs; 
	}
      }      
    }
  }

  Angle::
  Angle(std::vector<double> &paraValue, 
	std::vector< std::vector<size_t> > 
	&indValue )
  {
    //Do some checks on the parameters and variable indices
    //
    if( paraValue.size()!=1 ) {
      std::cerr << "Bending::Angle::Angle() "
		<< "One parameter, k_bend, should be provided." << std::endl;
      exit(EXIT_FAILURE);
    }
    if( indValue.size() != 1 || indValue[0].size() != 1 ) {
      std::cerr << "Bending::Angle::Angle() "
		<< "One index level with one index (angle) given." 
		<< std::endl;
      exit(EXIT_FAILURE);
    }
    //Set the variable values
    //
    setId("Bending::Angle");
    setParameter(paraValue);  
    setVariableIndex(indValue);
  }
  
  void Angle::
  derivs(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs )
  {
    size_t numCells = T.numCell();
    size_t dimension = T.vertex(0).numPosition();
    size_t Ti = variableIndex(0,0);
    for (size_t i=0; i<numCells; ++i) {
      size_t numWalls = T.cell(i).numWall();
      for (size_t k=0; k<numWalls; ++k) {
	size_t kp = k<numWalls-1 ? k+1 : 0;
	size_t km = k>0 ? k-1 : numWalls-1;

	// Get global vertex indices
	size_t j = T.cell(i).vertex(k)->index();
	size_t jp = T.cell(i).vertex(kp)->index();
	size_t jm = T.cell(i).vertex(km)->index();
	
	// Get edges
	//size_t ep = T.cell(i).wall(k)->index();
	size_t em = T.cell(i).wall(km)->index();
	
	// Calculate angle \Theta = acos(n1 \dot n2) - \pi
	double f = 0.; //scalar product
	double Lp = 0.;
	double Lm = 0.;
	for (size_t d=0; d<dimension; ++d) {
	  double dm = vertexData[j][d]-vertexData[jm][d]; 
	  double dp = vertexData[jp][d]-vertexData[j][d]; 
	  f += dm*dp;
	  Lp += dp*dp; 
	  Lm += dm*dm;
	}
	Lp = std::sqrt(Lp);
	Lm = std::sqrt(Lm);
	double gDenom = 1.0/(Lp*Lm);// 1./g 
	double F = f*gDenom;
	// Ad hoc way to avoid inf for derivative...
	if (F>0.999) F=0.999;
	else if (F<-0.999) F=-0.999;

	double theta = std::acos(F) - 3.14159; 

	double f0 = parameter(0)*(theta-wallData[em][Ti])*gDenom/(std::sqrt(1-F*F));
	double f1 = f*gDenom*Lm/Lp;
	double f2 = f*gDenom*Lp/Lm;
	//std::cerr << i << " " << j << " " << em << "\t" << theta << " " << wallData[em][Ti] << "\t"
	//<< F << " " << f0 << " " << f1 << " " << f2 << std::endl;  

	// Update for each dimension the contribution to vertices jm, j, jp
	for (size_t d=0; d<dimension; ++d) {
	  double dm = vertexData[j][d]-vertexData[jm][d]; 
	  double dp = vertexData[jp][d]-vertexData[j][d]; 
	  vertexDerivs[jm][d] += f0*(dm*f2-dp); 
	  vertexDerivs[jm][d] += f0*(dp*(1.+f1)-dm*(1+f2)); 
	  vertexDerivs[jp][d] += f0*(dm-dp*f1); 
	}
      }      
    }
  }

  AngleInitiate::
  AngleInitiate(std::vector<double> &paraValue, 
		std::vector< std::vector<size_t> > 
		&indValue )
  {
    //Do some checks on the parameters and variable indices
    //
    if( paraValue.size()!=0 ) {
      std::cerr << "Bending::AngleInitiate::AngleInitiate() "
		<< "No parameters are used." << std::endl;
      exit(EXIT_FAILURE);
    }
    if( indValue.size() != 1 || indValue[0].size() != 1 ) {
      std::cerr << "Bending::AngleInitiate::AngleInitiate() "
		<< "One index level with angle index given (angle stored as wall variable)." 
		<< std::endl;
      exit(EXIT_FAILURE);
    }
    //Set the variable values
    //
    setId("Bending::AngleInitiate");
    setParameter(paraValue);  
    setVariableIndex(indValue);
  }
  
  void AngleInitiate::
  derivs(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs )
  {
  }

  void AngleInitiate::
  initiate(Tissue &T,
	   DataMatrix &cellData,
	   DataMatrix &wallData,
	   DataMatrix &vertexData,
	   DataMatrix &cellDerivs,
	   DataMatrix &wallDerivs,
	   DataMatrix &vertexDerivs )
  {
    size_t numCells = T.numCell();
    size_t dimension = T.vertex(0).numPosition();
    size_t Ti = variableIndex(0,0);
    for (size_t i=0; i<numCells; ++i) {
      size_t numWalls = T.cell(i).numWall();
      for (size_t k=0; k<numWalls; ++k) {
	size_t kp = k<numWalls-1 ? k+1 : 0;
	size_t km = k>0 ? k-1 : numWalls-1;
	
	// Get global vertex indices
	size_t j = T.cell(i).vertex(k)->index();
	size_t jp = T.cell(i).vertex(kp)->index();
	size_t jm = T.cell(i).vertex(km)->index();
	
	// Get edges
	//size_t ep = T.cell(i).wall(k)->index();
	size_t em = T.cell(i).wall(km)->index();
	
	// Calculate angle \Theta = acos(n1 \dot n2) - \pi
	double scalarProd = 0.;
	double Lp = 0.;
	double Lm = 0.;
	for (size_t d=0; d<dimension; ++d) {
	  double dm = vertexData[j][d]-vertexData[jm][d]; 
	  double dp = vertexData[jp][d]-vertexData[j][d]; 
	  scalarProd += dm*dp;
	  Lp += dp*dp; 
	  Lm += dm*dm;
	}
	//std::cerr << em << " " << scalarProd << " " << Lp << " " << Lm << "\t";
	Lp = std::sqrt(Lp);
	Lm = std::sqrt(Lm);
	scalarProd /= Lp*Lm;
	if (scalarProd>1.0) scalarProd=1.0;//To make sure it is not giving NaN
	//std::cerr << em << " " << scalarProd << " " << Lp << " " << Lm << "\t";
	double theta = std::acos(scalarProd) - 3.14159; 
	// Setting 'stored' angle to current angle
	wallData[em][Ti] = theta;
	//std::cerr << em << " " << wallData[em][Ti] << std::endl;
      }
    }      
  }
  
  AngleRelax::
  AngleRelax(std::vector<double> &paraValue, 
	     std::vector< std::vector<size_t> > 
	     &indValue )
  {
    //Do some checks on the parameters and variable indices
    //
    if( paraValue.size()!=1 ) {
      std::cerr << "Bending::AngleRelax::AngleRelax() "
		<< "One parameter, k_bend, used." << std::endl;
      exit(EXIT_FAILURE);
    }
    if( indValue.size() != 1 || indValue[0].size() != 1 ) {
      std::cerr << "Bending::AngleRelax::AngleRelax() "
		<< "One index level with angle index given (angle stored as wall variable)." 
		<< std::endl;
      exit(EXIT_FAILURE);
    }
    //Set the variable values
    //
    setId("Bending::AngleRelax");
    setParameter(paraValue);  
    setVariableIndex(indValue);
  }
  
  void AngleRelax::
  derivs(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs )
  {
    size_t numCells = T.numCell();
    size_t dimension = T.vertex(0).numPosition();
    size_t Ti = variableIndex(0,0);
    for (size_t i=0; i<numCells; ++i) {
      size_t numWalls = T.cell(i).numWall();
      for (size_t k=0; k<numWalls; ++k) {
	size_t kp = k<numWalls-1 ? k+1 : 0;
	size_t km = k>0 ? k-1 : numWalls-1;
	
	// Get global vertex indices
	size_t j = T.cell(i).vertex(k)->index();
	size_t jp = T.cell(i).vertex(kp)->index();
	size_t jm = T.cell(i).vertex(km)->index();
	
	// Get edges
	//size_t ep = T.cell(i).wall(k)->index();
	size_t em = T.cell(i).wall(km)->index();
	
	// Calculate angle \Theta = acos(n1 \dot n2) - \pi
	double scalarProd = 0.;
	double Lp = 0.;
	double Lm = 0.;
	for (size_t d=0; d<dimension; ++d) {
	  double dm = vertexData[j][d]-vertexData[jm][d]; 
	  double dp = vertexData[jp][d]-vertexData[j][d]; 
	  scalarProd += dm*dp;
	  Lp += dp*dp; 
	  Lm += dm*dm;
	}
	Lp = std::sqrt(Lp);
	Lm = std::sqrt(Lm);
	scalarProd /= Lp*Lm;
	if (scalarProd>1.0) scalarProd=1.0;
	double theta = std::acos(scalarProd) - 3.14159; 
	// Update stored angle towards current angle
	wallDerivs[em][Ti] -= parameter(0)*(wallData[em][Ti]-theta); 
      }
    }      
  }

}


  
