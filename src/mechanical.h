/**
 * Filename     : mechanical.h
 * Description  : Classes describing mechanical updates
 * Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
 * Created      : April 2006
 * Revision     : $Id:$
 */
#ifndef MECHANICAL_H
#define MECHANICAL_H

#include"tissue.h"
#include"baseReaction.h"
#include<cmath>

//!Updates vertices from a cell pressure potential
class VertexFromCellPressure : public BaseReaction {
  
 public:
  
  VertexFromCellPressure(std::vector<double> &paraValue, 
			 std::vector< std::vector<size_t> > &indValue );
  
  void derivs(Tissue &T,
	      std::vector< std::vector<double> > &cellData,
	      std::vector< std::vector<double> > &wallData,
	      std::vector< std::vector<double> > &vertexData,
	      std::vector< std::vector<double> > &cellDerivs,
	      std::vector< std::vector<double> > &wallDerivs,
	      std::vector< std::vector<double> > &vertexDerivs );
};

//!Updates vertices from a cell pressure potential
class VertexFromCellPressureVolumeNormalized : public BaseReaction {
  
 public:
  
  VertexFromCellPressureVolumeNormalized(std::vector<double> &paraValue, 
																				 std::vector< std::vector<size_t> > &indValue );
  
  void derivs(Tissue &T,
							std::vector< std::vector<double> > &cellData,
							std::vector< std::vector<double> > &wallData,
							std::vector< std::vector<double> > &vertexData,
							std::vector< std::vector<double> > &cellDerivs,
							std::vector< std::vector<double> > &wallDerivs,
							std::vector< std::vector<double> > &vertexDerivs );
};

//!Updates vertices from a cell pressure potential
class VertexFromCellPressureThresholdFromMaxPos : public BaseReaction {
  
 public:
  
  VertexFromCellPressureThresholdFromMaxPos(std::vector<double> &paraValue, 
																						std::vector< std::vector<size_t> > &indValue );
  
  void derivs(Tissue &T,
							std::vector< std::vector<double> > &cellData,
							std::vector< std::vector<double> > &wallData,
							std::vector< std::vector<double> > &vertexData,
							std::vector< std::vector<double> > &cellDerivs,
							std::vector< std::vector<double> > &wallDerivs,
							std::vector< std::vector<double> > &vertexDerivs );
};

//!Updates vertices from a cell 'pressure' potential for internal cells
class VertexFromCellInternalPressure : public BaseReaction {
  
 public:
  
  VertexFromCellInternalPressure(std::vector<double> &paraValue, 
				 std::vector< std::vector<size_t> > &indValue );
  
  void derivs(Tissue &T,
	      std::vector< std::vector<double> > &cellData,
	      std::vector< std::vector<double> > &wallData,
	      std::vector< std::vector<double> > &vertexData,
	      std::vector< std::vector<double> > &cellDerivs,
	      std::vector< std::vector<double> > &wallDerivs,
	      std::vector< std::vector<double> > &vertexDerivs );
};

//!Updates vertices from cells via a power diagram potential
class VertexFromCellPowerdiagram : public BaseReaction {
  
 public:
  
  VertexFromCellPowerdiagram(std::vector<double> &paraValue, 
			     std::vector< std::vector<size_t> > 
			     &indValue );
  
  void derivs(Tissue &T,
	      std::vector< std::vector<double> > &cellData,
	      std::vector< std::vector<double> > &wallData,
	      std::vector< std::vector<double> > &vertexData,
	      std::vector< std::vector<double> > &cellDerivs,
	      std::vector< std::vector<double> > &wallDerivs,
	      std::vector< std::vector<double> > &vertexDerivs );
};

//!Applies a force towards or from origo on vertices specified by indices
class VertexForceOrigoFromIndex : public BaseReaction {
  
 public:
  
  VertexForceOrigoFromIndex(std::vector<double> &paraValue, 
			  std::vector< std::vector<size_t> > 
			  &indValue );
  
  void derivs(Tissue &T,
	      std::vector< std::vector<double> > &cellData,
	      std::vector< std::vector<double> > &wallData,
	      std::vector< std::vector<double> > &vertexData,
	      std::vector< std::vector<double> > &cellDerivs,
	      std::vector< std::vector<double> > &wallDerivs,
	      std::vector< std::vector<double> > &vertexDerivs );
};

//!Applies a force towards or from origo on vertices of cells
class CellForceOrigoFromIndex : public BaseReaction {
  
 public:
  
  CellForceOrigoFromIndex(std::vector<double> &paraValue, 
			  std::vector< std::vector<size_t> > 
			  &indValue );
  
  void derivs(Tissue &T,
	      std::vector< std::vector<double> > &cellData,
	      std::vector< std::vector<double> > &wallData,
	      std::vector< std::vector<double> > &vertexData,
	      std::vector< std::vector<double> > &cellDerivs,
	      std::vector< std::vector<double> > &wallDerivs,
	      std::vector< std::vector<double> > &vertexDerivs );
};

//!Applies a force towards or from a Cylinder surface
class CylinderForce : public BaseReaction {
  
 public:
  
  CylinderForce(std::vector<double> &paraValue, 
								std::vector< std::vector<size_t> > 
								&indValue );
  
  void derivs(Tissue &T,
							std::vector< std::vector<double> > &cellData,
							std::vector< std::vector<double> > &wallData,
							std::vector< std::vector<double> > &vertexData,
							std::vector< std::vector<double> > &cellDerivs,
							std::vector< std::vector<double> > &wallDerivs,
							std::vector< std::vector<double> > &vertexDerivs );
};

//!Applies a force towards or from a SphereCylinder surface
class SphereCylinderForce : public BaseReaction {
  
 public:
  
  SphereCylinderForce(std::vector<double> &paraValue, 
		      std::vector< std::vector<size_t> > 
		      &indValue );
  
  void derivs(Tissue &T,
	      std::vector< std::vector<double> > &cellData,
	      std::vector< std::vector<double> > &wallData,
	      std::vector< std::vector<double> > &vertexData,
	      std::vector< std::vector<double> > &cellDerivs,
	      std::vector< std::vector<double> > &wallDerivs,
	      std::vector< std::vector<double> > &vertexDerivs );
};

//!Applies a force towards a spherecylinder surface with defined radius
class SphereCylinderForceFromRadius : public BaseReaction {
  
 public:
  
  SphereCylinderForceFromRadius(std::vector<double> &paraValue, 
				std::vector< std::vector<size_t> > 
				&indValue );
  
  void derivs(Tissue &T,
	      std::vector< std::vector<double> > &cellData,
	      std::vector< std::vector<double> > &wallData,
	      std::vector< std::vector<double> > &vertexData,
	      std::vector< std::vector<double> > &cellDerivs,
	      std::vector< std::vector<double> > &wallDerivs,
	      std::vector< std::vector<double> > &vertexDerivs );
};

//!Applies a force perpendicular to a defined wall of infinite size
class InfiniteWallForce : public BaseReaction {
  
 public:
  
  InfiniteWallForce(std::vector<double> &paraValue, 
		    std::vector< std::vector<size_t> > 
		    &indValue );
  
  void derivs(Tissue &T,
	      std::vector< std::vector<double> > &cellData,
	      std::vector< std::vector<double> > &wallData,
	      std::vector< std::vector<double> > &vertexData,
	      std::vector< std::vector<double> > &cellDerivs,
	      std::vector< std::vector<double> > &wallDerivs,
	      std::vector< std::vector<double> > &vertexDerivs );
};

//!Applies a force on epidermal vertices
class EpidermalVertexForce : public BaseReaction {
  
 public:
  
  EpidermalVertexForce(std::vector<double> &paraValue, 
		    std::vector< std::vector<size_t> > 
		    &indValue );
  
  void derivs(Tissue &T,
	      std::vector< std::vector<double> > &cellData,
	      std::vector< std::vector<double> > &wallData,
	      std::vector< std::vector<double> > &vertexData,
	      std::vector< std::vector<double> > &cellDerivs,
	      std::vector< std::vector<double> > &wallDerivs,
	      std::vector< std::vector<double> > &vertexDerivs );
};

class VertexFromPressureExperimental : public BaseReaction
{  
 public:
	VertexFromPressureExperimental(std::vector<double> &paraValue, 
							 std::vector< std::vector<size_t> > &indValue);
	
	void derivs(Tissue &T,
			  std::vector< std::vector<double> > &cellData,
			  std::vector< std::vector<double> > &wallData,
			  std::vector< std::vector<double> > &vertexData,
			  std::vector< std::vector<double> > &cellDerivs,
			  std::vector< std::vector<double> > &wallDerivs,
			  std::vector< std::vector<double> > &vertexDerivs);
	double polygonArea(std::vector< std::pair<double, double> > vertices);
};

class CellVolumeExperimental : public BaseReaction
{
 public:
	CellVolumeExperimental(std::vector<double> &paraValue,
					   std::vector< std::vector<size_t> > &indValue);
	
	void derivs(Tissue &T,
			  std::vector< std::vector<double> > &cellData,
			  std::vector< std::vector<double> > &wallData,
			  std::vector< std::vector<double> > &vertexData,
			  std::vector< std::vector<double> > &cellDerivs,
			  std::vector< std::vector<double> > &wallDerivs,
			  std::vector< std::vector<double> > &vertexDerivs);
};

class EpidermalRadialForce : public BaseReaction
{
 public:
	EpidermalRadialForce(std::vector<double> &paraValue,
					 std::vector< std::vector<size_t> > &indValue);
	
	void derivs(Tissue &T,
			  std::vector< std::vector<double> > &cellData,
			  std::vector< std::vector<double> > &wallData,
			  std::vector< std::vector<double> > &vertexData,
			  std::vector< std::vector<double> > &cellDerivs,
			  std::vector< std::vector<double> > &wallDerivs,
			  std::vector< std::vector<double> > &vertexDerivs);
};

class PerpendicularWallPressure : public BaseReaction
{
 public:
	PerpendicularWallPressure(std::vector<double> &paraValue,
						 std::vector< std::vector<size_t> > &indValue);
	
	void derivs(Tissue &T,
			  std::vector< std::vector<double> > &cellData,
			  std::vector< std::vector<double> > &wallData,
			  std::vector< std::vector<double> > &vertexData,
			  std::vector< std::vector<double> > &cellDerivs,
			  std::vector< std::vector<double> > &wallDerivs,
			  std::vector< std::vector<double> > &vertexDerivs);
};

class ContinousMTDirection : public BaseReaction
{
 public:
	ContinousMTDirection(std::vector<double> &paraValue,
			     std::vector< std::vector<size_t> > &indValue);
	
	void derivs(Tissue &T,
		    std::vector< std::vector<double> > &cellData,
		    std::vector< std::vector<double> > &wallData,
		    std::vector< std::vector<double> > &vertexData,
		    std::vector< std::vector<double> > &cellDerivs,
		    std::vector< std::vector<double> > &wallDerivs,
		    std::vector< std::vector<double> > &vertexDerivs);
};

class VertexFromCellPlaneSphereCylinder : public BaseReaction
{
 public:
	VertexFromCellPlaneSphereCylinder(std::vector<double> &paraValue,
																		std::vector< std::vector<size_t> > &indValue);
	
	void derivs(Tissue &T,
		    std::vector< std::vector<double> > &cellData,
		    std::vector< std::vector<double> > &wallData,
		    std::vector< std::vector<double> > &vertexData,
		    std::vector< std::vector<double> > &cellDerivs,
		    std::vector< std::vector<double> > &wallDerivs,
		    std::vector< std::vector<double> > &vertexDerivs);
};

class VertexFromCellPlaneSphereCylinderConcentrationHill : public BaseReaction
{
 public:
	VertexFromCellPlaneSphereCylinderConcentrationHill(std::vector<double> &paraValue,
																										 std::vector< std::vector<size_t> > &indValue);
	
	void derivs(Tissue &T,
							std::vector< std::vector<double> > &cellData,
							std::vector< std::vector<double> > &wallData,
							std::vector< std::vector<double> > &vertexData,
							std::vector< std::vector<double> > &cellDerivs,
							std::vector< std::vector<double> > &wallDerivs,
							std::vector< std::vector<double> > &vertexDerivs);
};

// Do not use this reaction. Restricted area (unless you are a developer).
class DebugReaction : public BaseReaction
{
public:
	DebugReaction(std::vector<double> &paraValue,
								std::vector< std::vector<size_t> > &indValue);
	
	void derivs(Tissue &T,
							std::vector< std::vector<double> > &cellData,
							std::vector< std::vector<double> > &wallData,
							std::vector< std::vector<double> > &vertexData,
							std::vector< std::vector<double> > &cellDerivs,
							std::vector< std::vector<double> > &wallDerivs,
							std::vector< std::vector<double> > &vertexDerivs);
};

#endif
