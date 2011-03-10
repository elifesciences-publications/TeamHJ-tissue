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
  
  ///
  /// @brief Main constructor
  ///
  /// This is the main constructor which sets the parameters and variable
  /// indices that defines the reaction.
  ///
  /// @param paraValue vector with parameters
  ///
  /// @param indValue vector of vectors with variable indices
  ///
  /// @see BaseReaction::createReaction(std::vector<double> &paraValue,...)
  ///
  VertexFromCellPressure(std::vector<double> &paraValue, 
			 std::vector< std::vector<size_t> > &indValue );
  
  ///
  /// @brief Derivative function for this reaction class
  ///
  /// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
  ///
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
/*! A spring force in a perpendicular direction is applied. Note, the
  wall can only be defined along coordinate axes.
*/
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
/*! A spring force in a perpendicular direction is applied. Note, the
  wall can only be defined along coordinate axes.
*/
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

///
/// @brief Updates vertices from a 'pressure' term defined to act in the cell normal
/// direction.
///
/// This function calculates the area of a cell and then distribute a force 'outwards'
/// among the cell vertices. It relies on that the PCA cell planes have been calculated.
/// A cell contributes to a vertex update with
///
/// @f[ \frac{dx_{i}}{dt} = p_{0} A n_{i} / N_{vertex} @f]
///
/// where @f$p_{0}@f$ is a 'pressure' parameter, A is the cell area @f$n_{i}@f$ is the 
/// cell normal component and @f$N_{vertex}@f$ is the number of vertices for the cell.
/// An additional parameter @f$p_{2}@f$ can be used to not include the area factor if
/// set to zero (normally it should be set to 1).
///
/// In a model file the reaction is defined as
///
/// @verbatim
/// VertexFromCellPlane 2 0
/// P A_flag
/// @endverbatim
///
/// @see CalculatePCAPlane
///
class VertexFromCellPlane : public BaseReaction
{
 public:
  VertexFromCellPlane(std::vector<double> &paraValue,
		      std::vector< std::vector<size_t> > &indValue);
	
  void derivs(Tissue &T,
	      std::vector< std::vector<double> > &cellData,
	      std::vector< std::vector<double> > &wallData,
	      std::vector< std::vector<double> > &vertexData,
	      std::vector< std::vector<double> > &cellDerivs,
	      std::vector< std::vector<double> > &wallDerivs,
	      std::vector< std::vector<double> > &vertexDerivs);
};

class VertexFromCellPlaneSpatial : public BaseReaction
{
 private:
  
  double Kpow_;
  
 public:
  VertexFromCellPlaneSpatial(std::vector<double> &paraValue,
			     std::vector< std::vector<size_t> > &indValue);
  
  void derivs(Tissue &T,
	      std::vector< std::vector<double> > &cellData,
	      std::vector< std::vector<double> > &wallData,
	      std::vector< std::vector<double> > &vertexData,
	      std::vector< std::vector<double> > &cellDerivs,
	      std::vector< std::vector<double> > &wallDerivs,
	      std::vector< std::vector<double> > &vertexDerivs);
};

///
/// @brief Same as VertexFromCellPlane but with the strength dependent
/// on a molecular concentration
///
/// This class is the same 'pressure' from inside update as
/// VertexFromCellPlane with the difference that the strength of the
/// resulting force depends on the cellular concentration of a
/// molecule described with a Hill formalism.
///
/// The force is described by:
///
/// \f[ 
/// F = \frac{A_{cell}}{N_{vertex}}(p_{0} + 
/// p_{1} \frac{C^{p_{3}}}{p_{2}^{p_{3}}+C^{p_{3}}}
/// \f]
///
/// where the area factor \f$A_{cell}\f$ is present if \f$p_{4}=1\f$,
/// and C is the molecular concentration.
///
class VertexFromCellPlaneConcentrationHill : public BaseReaction
{
 private:
  
  double Kpow_;
  
 public:
  VertexFromCellPlaneConcentrationHill(std::vector<double> &paraValue,
				       std::vector< std::vector<size_t> > &indValue);
  
  void derivs(Tissue &T,
	      std::vector< std::vector<double> > &cellData,
	      std::vector< std::vector<double> > &wallData,
	      std::vector< std::vector<double> > &vertexData,
	      std::vector< std::vector<double> > &cellDerivs,
	      std::vector< std::vector<double> > &wallDerivs,
	      std::vector< std::vector<double> > &vertexDerivs);
};

class VertexFromCellPlaneNormalized : public BaseReaction
{
 public:
  VertexFromCellPlaneNormalized(std::vector<double> &paraValue,
				std::vector< std::vector<size_t> > &indValue);
  
  void derivs(Tissue &T,
	      std::vector< std::vector<double> > &cellData,
	      std::vector< std::vector<double> > &wallData,
	      std::vector< std::vector<double> > &vertexData,
	      std::vector< std::vector<double> > &cellDerivs,
	      std::vector< std::vector<double> > &wallDerivs,
	      std::vector< std::vector<double> > &vertexDerivs);
};

class VertexFromCellPlaneNormalizedSpatial : public BaseReaction
{
 private:
  
  double Kpow_;
  
 public:
  VertexFromCellPlaneNormalizedSpatial(std::vector<double> &paraValue,
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
