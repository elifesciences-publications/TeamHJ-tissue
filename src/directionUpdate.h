//
// Filename     : directionUpdate.h
// Description  : Classes describing growth updates
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : April 2006
// Revision     : $Id:$
//
#ifndef DIRECTIONUPDATE_H
#define DIRECTIONUPDATE_H

#include<cmath>

#include"tissue.h"
#include"baseDirectionUpdate.h"

///
/// @brief No update of the direction is applied, but a direction is created
///
/// This can be used to create a direction but when it is not updated.
///
class StaticDirection : public BaseDirectionUpdate {
  
 public:
  
  StaticDirection(std::vector<double> &paraValue, 
									std::vector< std::vector<size_t> > 
									&indValue );
  
  void initiate(Tissue &T,
								std::vector< std::vector<double> > &cellData,
								std::vector< std::vector<double> > &wallData,
								std::vector< std::vector<double> > &vertexData,
								std::vector< std::vector<double> > &cellDerivs,
								std::vector< std::vector<double> > &wallDerivs,
								std::vector< std::vector<double> > &vertexDerivs );
  void update(Tissue &T, double h,
							std::vector< std::vector<double> > &cellData,
							std::vector< std::vector<double> > &wallData,
							std::vector< std::vector<double> > &vertexData,
							std::vector< std::vector<double> > &cellDerivs,
							std::vector< std::vector<double> > &wallDerivs,
							std::vector< std::vector<double> > &vertexDerivs );
};

///
/// @brief The direction follows the direction of one of the cell walls
///
class WallDirection : public BaseDirectionUpdate {
  
 public:
  
  WallDirection(std::vector<double> &paraValue, 
								std::vector< std::vector<size_t> > 
								&indValue );
  
  void initiate(Tissue &T,
								std::vector< std::vector<double> > &cellData,
								std::vector< std::vector<double> > &wallData,
								std::vector< std::vector<double> > &vertexData,
								std::vector< std::vector<double> > &cellDerivs,
								std::vector< std::vector<double> > &wallDerivs,
								std::vector< std::vector<double> > &vertexDerivs );
  void update(Tissue &T, double h,
								std::vector< std::vector<double> > &cellData,
								std::vector< std::vector<double> > &wallData,
								std::vector< std::vector<double> > &vertexData,
								std::vector< std::vector<double> > &cellDerivs,
								std::vector< std::vector<double> > &wallDerivs,
								std::vector< std::vector<double> > &vertexDerivs );
};

///
/// @brief The direction follows the direction of the cell strain
///
/// It is based on the Goodall and Green (1986) calculation of the 
/// maximal strain direction.
///
class StrainDirection : public BaseDirectionUpdate {
  
 public:
  
  StrainDirection(std::vector<double> &paraValue, 
									std::vector< std::vector<size_t> > 
									&indValue );
  
  void initiate(Tissue &T,
								std::vector< std::vector<double> > &cellData,
								std::vector< std::vector<double> > &wallData,
								std::vector< std::vector<double> > &vertexData,
								std::vector< std::vector<double> > &cellDerivs,
								std::vector< std::vector<double> > &wallDerivs,
								std::vector< std::vector<double> > &vertexDerivs );
  void update(Tissue &T, double h,
								std::vector< std::vector<double> > &cellData,
								std::vector< std::vector<double> > &wallData,
								std::vector< std::vector<double> > &vertexData,
								std::vector< std::vector<double> > &cellDerivs,
								std::vector< std::vector<double> > &wallDerivs,
								std::vector< std::vector<double> > &vertexDerivs );
};

///
/// @brief Direction is updated from strain averaged over the wall strains
///
/// The strain is calculated for each wall and is used to calculate an
/// averaged direction for the cell.
/// 
class StrainDirectionWall : public BaseDirectionUpdate
{
 public:
  StrainDirectionWall(std::vector<double> &paraValue, std::vector< std::vector<size_t> > &indValue );
  
  void initiate(Tissue &T,
								std::vector< std::vector<double> > &cellData,
								std::vector< std::vector<double> > &wallData,
								std::vector< std::vector<double> > &vertexData,
								std::vector< std::vector<double> > &cellDerivs,
								std::vector< std::vector<double> > &wallDerivs,
								std::vector< std::vector<double> > &vertexDerivs);
	
  void update(Tissue &T, double h,
							std::vector< std::vector<double> > &cellData,
							std::vector< std::vector<double> > &wallData,
							std::vector< std::vector<double> > &vertexData,
							std::vector< std::vector<double> > &cellDerivs,
							std::vector< std::vector<double> > &wallDerivs,
							std::vector< std::vector<double> > &vertexDerivs);
};

///
/// @brief The direction follows the direction of one molecule
///
class GradientDirection : public BaseDirectionUpdate {
  
 public:
  
  GradientDirection(std::vector<double> &paraValue, 
										std::vector< std::vector<size_t> > 
										&indValue );
  
  void initiate(Tissue &T,
								std::vector< std::vector<double> > &cellData,
								std::vector< std::vector<double> > &wallData,
								std::vector< std::vector<double> > &vertexData,
								std::vector< std::vector<double> > &cellDerivs,
								std::vector< std::vector<double> > &wallDerivs,
								std::vector< std::vector<double> > &vertexDerivs );
  void update(Tissue &T, double h,
							std::vector< std::vector<double> > &cellData,
							std::vector< std::vector<double> > &wallData,
							std::vector< std::vector<double> > &vertexData,
							std::vector< std::vector<double> > &cellDerivs,
							std::vector< std::vector<double> > &wallDerivs,
							std::vector< std::vector<double> > &vertexDerivs );
};

class WallStressDirection : public BaseDirectionUpdate
{
public:
  WallStressDirection(std::vector<double> &paraValue, std::vector< std::vector<size_t> > &indValue );
  
  void initiate(Tissue &T,
								std::vector< std::vector<double> > &cellData,
								std::vector< std::vector<double> > &wallData,
								std::vector< std::vector<double> > &vertexData,
								std::vector< std::vector<double> > &cellDerivs,
								std::vector< std::vector<double> > &wallDerivs,
								std::vector< std::vector<double> > &vertexDerivs);

  void update(Tissue &T, double h,
							std::vector< std::vector<double> > &cellData,
							std::vector< std::vector<double> > &wallData,
							std::vector< std::vector<double> > &vertexData,
							std::vector< std::vector<double> > &cellDerivs,
							std::vector< std::vector<double> > &wallDerivs,
							std::vector< std::vector<double> > &vertexDerivs);
};

class DoubleWallStressDirection : public BaseDirectionUpdate
{
 public:
  DoubleWallStressDirection(std::vector<double> &paraValue, std::vector< std::vector<size_t> > &indValue );
  
  void initiate(Tissue &T,
								std::vector< std::vector<double> > &cellData,
								std::vector< std::vector<double> > &wallData,
								std::vector< std::vector<double> > &vertexData,
								std::vector< std::vector<double> > &cellDerivs,
								std::vector< std::vector<double> > &wallDerivs,
								std::vector< std::vector<double> > &vertexDerivs);

  void update(Tissue &T, double h,
							std::vector< std::vector<double> > &cellData,
							std::vector< std::vector<double> > &wallData,
							std::vector< std::vector<double> > &vertexData,
							std::vector< std::vector<double> > &cellDerivs,
							std::vector< std::vector<double> > &wallDerivs,
							std::vector< std::vector<double> > &vertexDerivs);
};

class StretchDirection : public BaseDirectionUpdate
{
 public:
  StretchDirection(std::vector<double> &paraValue, std::vector< std::vector<size_t> > &indValue );
  
  void initiate(Tissue &T,
			 std::vector< std::vector<double> > &cellData,
			 std::vector< std::vector<double> > &wallData,
			 std::vector< std::vector<double> > &vertexData,
			 std::vector< std::vector<double> > &cellDerivs,
			 std::vector< std::vector<double> > &wallDerivs,
			 std::vector< std::vector<double> > &vertexDerivs);

  void update(Tissue &T, double h,
		    std::vector< std::vector<double> > &cellData,
		    std::vector< std::vector<double> > &wallData,
		    std::vector< std::vector<double> > &vertexData,
		    std::vector< std::vector<double> > &cellDerivs,
		    std::vector< std::vector<double> > &wallDerivs,
		    std::vector< std::vector<double> > &vertexDerivs);
};

class PCAPlaneDirection : public BaseDirectionUpdate
{
public:
  PCAPlaneDirection(std::vector<double> &paraValue, std::vector< std::vector<size_t> > &indValue );
  
	void initiate(Tissue &T,
								std::vector< std::vector<double> > &cellData,
								std::vector< std::vector<double> > &wallData,
								std::vector< std::vector<double> > &vertexData,
								std::vector< std::vector<double> > &cellDerivs,
								std::vector< std::vector<double> > &wallDerivs,
								std::vector< std::vector<double> > &vertexDerivs);
	
	void update(Tissue &T, double h,
							std::vector< std::vector<double> > &cellData,
							std::vector< std::vector<double> > &wallData,
							std::vector< std::vector<double> > &vertexData,
							std::vector< std::vector<double> > &cellDerivs,
							std::vector< std::vector<double> > &wallDerivs,
							std::vector< std::vector<double> > &vertexDerivs);
};

class VertexStressDirection : public BaseDirectionUpdate
{
public:
	VertexStressDirection(std::vector<double> &paraValue, std::vector< std::vector<size_t> > &indValue);
	
	void initiate(Tissue &T,
		std::vector< std::vector<double> > &cellData,
		std::vector< std::vector<double> > &wallData,
		std::vector< std::vector<double> > &vertexData,
		std::vector< std::vector<double> > &cellDerivs,
		std::vector< std::vector<double> > &wallDerivs,
		std::vector< std::vector<double> > &vertexDerivs);
	
	void update(Tissue &T, double h,
		std::vector< std::vector<double> > &cellData,
		std::vector< std::vector<double> > &wallData,
		std::vector< std::vector<double> > &vertexData,
		std::vector< std::vector<double> > &cellDerivs,
		std::vector< std::vector<double> > &wallDerivs,
		std::vector< std::vector<double> > &vertexDerivs);
};

#endif
