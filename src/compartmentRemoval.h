/**
 * Filename     : compartmentRemoval.h
 * Description  : Classes describing compartment removal updates
 * Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
 * Created      : November 2006
 * Revision     : $Id:$
 */
#ifndef COMPARTMENTREMOVAL_H
#define COMPARTMENTREMOVAL_H

#include <cmath>

#include "tissue.h"
#include "baseCompartmentChange.h"

//!Removes a cell when position outside a radius from origo
class RemovalOutsideRadius : public BaseCompartmentChange {
  
 public:
  
  RemovalOutsideRadius(std::vector<double> &paraValue, 
											 std::vector< std::vector<size_t> > 
											 &indValue );
  
  int flag(Tissue *T,size_t i,
					 std::vector< std::vector<double> > &cellData,
					 std::vector< std::vector<double> > &wallData,
					 std::vector< std::vector<double> > &vertexData,
					 std::vector< std::vector<double> > &cellDerivs,
					 std::vector< std::vector<double> > &wallDerivs,
					 std::vector< std::vector<double> > &vertexDerivs );
  void update(Tissue* T,size_t i,
							std::vector< std::vector<double> > &cellData,
							std::vector< std::vector<double> > &wallData,
							std::vector< std::vector<double> > &vertexData,
							std::vector< std::vector<double> > &cellDerivs,
							std::vector< std::vector<double> > &wallDerivs,
							std::vector< std::vector<double> > &vertexDerivs );  
};

//!Removes a cell when position outside a radius from origo
class RemovalOutsideRadiusEpidermis : public BaseCompartmentChange {
  
 public:
  
  RemovalOutsideRadiusEpidermis(std::vector<double> &paraValue, 
																std::vector< std::vector<size_t> > 
																&indValue );
  
  int flag(Tissue *T,size_t i,
					 std::vector< std::vector<double> > &cellData,
					 std::vector< std::vector<double> > &wallData,
					 std::vector< std::vector<double> > &vertexData,
					 std::vector< std::vector<double> > &cellDerivs,
					 std::vector< std::vector<double> > &wallDerivs,
					 std::vector< std::vector<double> > &vertexDerivs );
  void update(Tissue* T,size_t i,
							std::vector< std::vector<double> > &cellData,
							std::vector< std::vector<double> > &wallData,
							std::vector< std::vector<double> > &vertexData,
							std::vector< std::vector<double> > &cellDerivs,
							std::vector< std::vector<double> > &wallDerivs,
							std::vector< std::vector<double> > &vertexDerivs );  
};

//!Removes a cell when position outside a distance from max
class RemovalOutsideMaxDistanceEpidermis : public BaseCompartmentChange {
  
 public:
  
  RemovalOutsideMaxDistanceEpidermis(std::vector<double> &paraValue, 
																		 std::vector< std::vector<size_t> > 
																		 &indValue );
  
  int flag(Tissue *T,size_t i,
					 std::vector< std::vector<double> > &cellData,
					 std::vector< std::vector<double> > &wallData,
					 std::vector< std::vector<double> > &vertexData,
					 std::vector< std::vector<double> > &cellDerivs,
					 std::vector< std::vector<double> > &wallDerivs,
					 std::vector< std::vector<double> > &vertexDerivs );
  void update(Tissue* T,size_t i,
							std::vector< std::vector<double> > &cellData,
							std::vector< std::vector<double> > &wallData,
							std::vector< std::vector<double> > &vertexData,
							std::vector< std::vector<double> > &cellDerivs,
							std::vector< std::vector<double> > &wallDerivs,
							std::vector< std::vector<double> > &vertexDerivs );  

 private:
	double max_;
};

//!Removes a cell when position outside a threshold value
class RemovalOutsidePosition : public BaseCompartmentChange {
  
public:
  
  RemovalOutsidePosition(std::vector<double> &paraValue, 
												 std::vector< std::vector<size_t> > 
												 &indValue );
  
  int flag(Tissue *T,size_t i,
					 std::vector< std::vector<double> > &cellData,
					 std::vector< std::vector<double> > &wallData,
					 std::vector< std::vector<double> > &vertexData,
					 std::vector< std::vector<double> > &cellDerivs,
					 std::vector< std::vector<double> > &wallDerivs,
					 std::vector< std::vector<double> > &vertexDerivs );
  void update(Tissue* T,size_t i,
							std::vector< std::vector<double> > &cellData,
							std::vector< std::vector<double> > &wallData,
							std::vector< std::vector<double> > &vertexData,
							std::vector< std::vector<double> > &cellDerivs,
							std::vector< std::vector<double> > &wallDerivs,
							std::vector< std::vector<double> > &vertexDerivs );  
};

#endif
