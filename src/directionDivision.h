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
class StaticDirection : public BaseDirectionDivision {
  
 public:
  
  StaticDirection(std::vector<double> &paraValue, 
									std::vector< std::vector<size_t> > 
									&indValue );
  
  void update(Tissue &T,size_t cellI,
							std::vector< std::vector<double> > &cellData,
							std::vector< std::vector<double> > &wallData,
							std::vector< std::vector<double> > &vertexData,
							std::vector< std::vector<double> > &cellDerivs,
							std::vector< std::vector<double> > &wallDerivs,
							std::vector< std::vector<double> > &vertexDerivs );
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
							std::vector< std::vector<double> > &cellData,
							std::vector< std::vector<double> > &wallData,
							std::vector< std::vector<double> > &vertexData,
							std::vector< std::vector<double> > &cellDerivs,
							std::vector< std::vector<double> > &wallDerivs,
							std::vector< std::vector<double> > &vertexDerivs );
};

#endif
