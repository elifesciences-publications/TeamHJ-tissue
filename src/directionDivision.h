//
// Filename     : directionDivision.h
// Description  : Classes describing growth updates
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : June 2007
// Revision     : $Id:$
//
#ifndef DIRECTIONDIVISION_H
#define DIRECTIONDIVISION_H

#include<cmath>

#include"tissue.h"
#include"baseDirectionDivision.h"

///
/// @brief The direction stays the same
///
class ParallellDirection : public BaseDirectionDivision {
  
 public:
  
  ParallellDirection(std::vector<double> &paraValue, 
		     std::vector< std::vector<size_t> > 
		     &indValue );
  
  void update(Tissue &T,size_t cellI,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      DataMatrix &cellDerivs,
	      DataMatrix &wallDerivs,
	      DataMatrix &vertexDerivs );
};

///
/// @brief The direction becomes perpendicular
///
class PerpendicularDirection : public BaseDirectionDivision {
  
 public:
  
  PerpendicularDirection(std::vector<double> &paraValue, 
			 std::vector< std::vector<size_t> > 
			 &indValue );
  
  void update(Tissue &T,size_t cellI,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      DataMatrix &cellDerivs,
	      DataMatrix &wallDerivs,
	      DataMatrix &vertexDerivs );
};

class RandomDirection : public BaseDirectionDivision
{
 public:
  RandomDirection(std::vector<double> &paraValue, std::vector< std::vector<size_t> > &indValue);
  
	void update(Tissue &T, size_t cellI,
		    DataMatrix &cellData,
		    DataMatrix &wallData,
		    DataMatrix &vertexData,
		    DataMatrix &cellDerivs,
		    DataMatrix &wallDerivs,
		    DataMatrix &vertexDerivs);
};

#endif
