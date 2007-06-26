/**
 * Filename     : compartmentDivision.h
 * Description  : Classes describing compartmentDivision updates
 * Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
 * Created      : July 2006
 * Revision     : $Id:$
 */
#ifndef COMPARTMENTDIVISION_H
#define COMPARTMENTDIVISION_H

#include <cmath>

#include "tissue.h"
#include "baseCompartmentChange.h"

//!Divides a cell when volume above a threshold
/*!Divides a cell when volume above a threshold. New wall is created
  prependicular to the longest cell wall.
 */
class DivisionVolumeViaLongestWall : public BaseCompartmentChange {
  
 public:
  
  DivisionVolumeViaLongestWall(std::vector<double> &paraValue, 
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

//!Divides a cell when volume above a threshold in 3D
/*!Divides a cell when volume above a threshold. Same as
  DivisionVolumeViaLongestWall but used for surfaces in 3D.
*/
class DivisionVolumeViaLongestWall3D : public BaseCompartmentChange {
  
 public:
  
  DivisionVolumeViaLongestWall3D(std::vector<double> &paraValue, 
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

//!Divides a cell when volume above a threshold
/*!Divides a cell when volume above a threshold. New wall is created
  prependicular to maximal strain rate.
 */
class DivisionVolumeViaStrain : public BaseCompartmentChange {
  
 public:
  
  DivisionVolumeViaStrain(std::vector<double> &paraValue, 
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

//!Divides a cell when volume above a threshold
/*!Divides a cell when volume above a threshold. New wall is created
  prependicular to direction given as cell variable.
 */
class DivisionVolumeViaDirection : public BaseCompartmentChange {
  
 public:
  
  DivisionVolumeViaDirection(std::vector<double> &paraValue, 
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

//!Divides a cell when volume above a threshold
/*!Divides a cell when volume above a threshold. New wall is created
  in a random direction through center of mass.
 */
class DivisionVolumeRandomDirection : public BaseCompartmentChange {
  
 public:
  
  DivisionVolumeRandomDirection(std::vector<double> &paraValue, 
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


class DivisionForceDirection : public BaseCompartmentChange
{
 public:
	DivisionForceDirection(std::vector<double> &paraValue, 
					   std::vector< std::vector<size_t> > &indValue);
	
	int flag(Tissue *T, size_t i,
		    std::vector< std::vector<double> > &cellData,
		    std::vector< std::vector<double> > &wallData,
		    std::vector< std::vector<double> > &vertexData,
		    std::vector< std::vector<double> > &cellDerivs,
		    std::vector< std::vector<double> > &wallDerivs,
		    std::vector< std::vector<double> > &vertexDerivs);
	void update(Tissue* T, size_t i,
			  std::vector< std::vector<double> > &cellData,
			  std::vector< std::vector<double> > &wallData,
			  std::vector< std::vector<double> > &vertexData,
			  std::vector< std::vector<double> > &cellDerivs,
			  std::vector< std::vector<double> > &wallDerivs,
			  std::vector< std::vector<double> > &vertexDerivs);  
};

#endif
