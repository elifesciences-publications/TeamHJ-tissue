//
// Filename     : mechanical.h
// Description  : Classes describing mechanical updates
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : April 2006
// Revision     : $Id:$
//
#ifndef MECHANICAL_H
#define MECHANICAL_H

#include"tissue.h"
#include"baseReaction.h"
#include<cmath>

///
/// @brief Updates vertices from a cell pressure potential
///
/// 
///
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
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      DataMatrix &cellDerivs,
	      DataMatrix &wallDerivs,
	      DataMatrix &vertexDerivs );
};

///
/// @brief Updates vertices from a cell pressure potential
///
/// This function determines the direction of the pressure force term
/// is from the position of the central mesh cell vertex to the center of the wall. 
/// Applies a force proportional to the pressure (parameter(0)) and the 
/// size of the wall.
///
/// @note Maybe it should rather be normal to the wall in the plane of the triangle?
///
class VertexFromCellPressurecenterTriangulation : public BaseReaction {
  
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
  VertexFromCellPressurecenterTriangulation(std::vector<double> &paraValue, 
					    std::vector< std::vector<size_t> > &indValue );
  
  ///
  /// @brief Derivative function for this reaction class
  ///
  /// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
  ///
  void derivs(Tissue &T,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      DataMatrix &cellDerivs,
	      DataMatrix &wallDerivs,
	      DataMatrix &vertexDerivs );
};

//!Updates vertices from a cell pressure potential
class VertexFromCellPressureVolumeNormalized : public BaseReaction {
  
 public:
  
  VertexFromCellPressureVolumeNormalized(std::vector<double> &paraValue, 
					 std::vector< std::vector<size_t> > &indValue );
  
  void derivs(Tissue &T,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      DataMatrix &cellDerivs,
	      DataMatrix &wallDerivs,
	      DataMatrix &vertexDerivs );
};

//!Updates vertices from a cell pressure potential
class VertexFromCellPressureThresholdFromMaxPos : public BaseReaction {
  
 public:
  
  VertexFromCellPressureThresholdFromMaxPos(std::vector<double> &paraValue, 
					    std::vector< std::vector<size_t> > &indValue );
  
  void derivs(Tissue &T,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      DataMatrix &cellDerivs,
	      DataMatrix &wallDerivs,
	      DataMatrix &vertexDerivs );
};

//!Updates vertices from a cell 'pressure' potential for internal cells
class VertexFromCellInternalPressure : public BaseReaction {
  
 public:
  
  VertexFromCellInternalPressure(std::vector<double> &paraValue, 
				 std::vector< std::vector<size_t> > &indValue );
  
  void derivs(Tissue &T,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      DataMatrix &cellDerivs,
	      DataMatrix &wallDerivs,
	      DataMatrix &vertexDerivs );
};

//!Updates vertices from cells via a power diagram potential
class VertexFromCellPowerdiagram : public BaseReaction {
  
 public:
  
  VertexFromCellPowerdiagram(std::vector<double> &paraValue, 
			     std::vector< std::vector<size_t> > 
			     &indValue );
  
  void derivs(Tissue &T,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      DataMatrix &cellDerivs,
	      DataMatrix &wallDerivs,
	      DataMatrix &vertexDerivs );
};

//!Applies a force towards or from origo on vertices specified by indices
class VertexForceOrigoFromIndex : public BaseReaction {
  
 public:
  
  VertexForceOrigoFromIndex(std::vector<double> &paraValue, 
			    std::vector< std::vector<size_t> > 
			    &indValue );
  
  void derivs(Tissue &T,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      DataMatrix &cellDerivs,
	      DataMatrix &wallDerivs,
	      DataMatrix &vertexDerivs );
};

//!Applies a force towards or from origo on vertices of cells
class CellForceOrigoFromIndex : public BaseReaction {
  
 public:
  
  CellForceOrigoFromIndex(std::vector<double> &paraValue, 
			  std::vector< std::vector<size_t> > 
			  &indValue );
  
  void derivs(Tissue &T,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      DataMatrix &cellDerivs,
	      DataMatrix &wallDerivs,
	      DataMatrix &vertexDerivs );
};

//!Applies a force towards or from a Cylinder surface
class CylinderForce : public BaseReaction {
  
 public:
  
  CylinderForce(std::vector<double> &paraValue, 
		std::vector< std::vector<size_t> > 
		&indValue );
  
  void derivs(Tissue &T,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      DataMatrix &cellDerivs,
	      DataMatrix &wallDerivs,
	      DataMatrix &vertexDerivs );
};

///
/// @brief Applies a force towards or from a SphereCylinder surface
///
class SphereCylinderForce : public BaseReaction {
  
 public:
  
  SphereCylinderForce(std::vector<double> &paraValue, 
		      std::vector< std::vector<size_t> > 
		      &indValue );
  
  void derivs(Tissue &T,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      DataMatrix &cellDerivs,
	      DataMatrix &wallDerivs,
	      DataMatrix &vertexDerivs );
};

//!Applies a force towards a spherecylinder surface with defined radius
class SphereCylinderForceFromRadius : public BaseReaction {
  
 public:
  
  SphereCylinderForceFromRadius(std::vector<double> &paraValue, 
				std::vector< std::vector<size_t> > 
				&indValue );
  
  void derivs(Tissue &T,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      DataMatrix &cellDerivs,
	      DataMatrix &wallDerivs,
	      DataMatrix &vertexDerivs );
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
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      DataMatrix &cellDerivs,
	      DataMatrix &wallDerivs,
	      DataMatrix &vertexDerivs );
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
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      DataMatrix &cellDerivs,
	      DataMatrix &wallDerivs,
	      DataMatrix &vertexDerivs );
};

class VertexFromPressureExperimental : public BaseReaction
{  
 public:
  VertexFromPressureExperimental(std::vector<double> &paraValue, 
				 std::vector< std::vector<size_t> > &indValue);
  
  void derivs(Tissue &T,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      DataMatrix &cellDerivs,
	      DataMatrix &wallDerivs,
	      DataMatrix &vertexDerivs);
  double polygonArea(std::vector< std::pair<double, double> > vertices);
};

class CellVolumeExperimental : public BaseReaction
{
 public:
  CellVolumeExperimental(std::vector<double> &paraValue,
			 std::vector< std::vector<size_t> > &indValue);
  
  void derivs(Tissue &T,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      DataMatrix &cellDerivs,
	      DataMatrix &wallDerivs,
	      DataMatrix &vertexDerivs);
};

class EpidermalRadialForce : public BaseReaction
{
 public:
  EpidermalRadialForce(std::vector<double> &paraValue,
		       std::vector< std::vector<size_t> > &indValue);
  
  void derivs(Tissue &T,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      DataMatrix &cellDerivs,
	      DataMatrix &wallDerivs,
	      DataMatrix &vertexDerivs);
};

class PerpendicularWallPressure : public BaseReaction
{
 public:
  PerpendicularWallPressure(std::vector<double> &paraValue,
			    std::vector< std::vector<size_t> > &indValue);
  
  void derivs(Tissue &T,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      DataMatrix &cellDerivs,
	      DataMatrix &wallDerivs,
	      DataMatrix &vertexDerivs);
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
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      DataMatrix &cellDerivs,
	      DataMatrix &wallDerivs,
	      DataMatrix &vertexDerivs);
};

class VertexFromCellPlaneSpatial : public BaseReaction
{
 private:
  
  double Kpow_;
  
 public:
  VertexFromCellPlaneSpatial(std::vector<double> &paraValue,
			     std::vector< std::vector<size_t> > &indValue);
  
  void derivs(Tissue &T,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      DataMatrix &cellDerivs,
	      DataMatrix &wallDerivs,
	      DataMatrix &vertexDerivs);
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
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      DataMatrix &cellDerivs,
	      DataMatrix &wallDerivs,
	      DataMatrix &vertexDerivs);
};

class VertexFromCellPlaneNormalized : public BaseReaction
{
 public:
  VertexFromCellPlaneNormalized(std::vector<double> &paraValue,
				std::vector< std::vector<size_t> > &indValue);
  
  void derivs(Tissue &T,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      DataMatrix &cellDerivs,
	      DataMatrix &wallDerivs,
	      DataMatrix &vertexDerivs);
};

class VertexFromCellPlaneNormalizedSpatial : public BaseReaction
{
 private:
  
  double Kpow_;
  
 public:
  VertexFromCellPlaneNormalizedSpatial(std::vector<double> &paraValue,
				       std::vector< std::vector<size_t> > &indValue);
  
  void derivs(Tissue &T,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      DataMatrix &cellDerivs,
	      DataMatrix &wallDerivs,
	      DataMatrix &vertexDerivs);
};

class VertexFromCellPlaneSphereCylinder : public BaseReaction
{
 public:
  VertexFromCellPlaneSphereCylinder(std::vector<double> &paraValue,
				    std::vector< std::vector<size_t> > &indValue);
  
  void derivs(Tissue &T,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      DataMatrix &cellDerivs,
	      DataMatrix &wallDerivs,
	      DataMatrix &vertexDerivs);
};

class VertexFromCellPlaneSphereCylinderConcentrationHill : public BaseReaction
{
 public:
  VertexFromCellPlaneSphereCylinderConcentrationHill(std::vector<double> &paraValue,
						     std::vector< std::vector<size_t> > &indValue);
  
  void derivs(Tissue &T,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      DataMatrix &cellDerivs,
	      DataMatrix &wallDerivs,
	      DataMatrix &vertexDerivs);
};

///
/// @brief Updates vertices from a 'pressure' term defined to act in the cell normal
/// direction for triangular cells only.
///
/// This function calculates the area of a cell and then distribute a force 'outwards'
/// among the cell vertices. It relies on that the cells are triangular.
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
/// VertexFromCellPlaneTriangular 2 0
/// P A_flag
/// @endverbatim
///
///
class VertexFromCellPlaneTriangular : public BaseReaction
{
 public:
  VertexFromCellPlaneTriangular(std::vector<double> &paraValue,
				std::vector< std::vector<size_t> > &indValue);
  
  void derivs(Tissue &T,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      DataMatrix &cellDerivs,
	      DataMatrix &wallDerivs,
	      DataMatrix &vertexDerivs);
};

///
/// @brief Updates list of vertices with a given force applied
///
/// In a model file the reaction is defined as
///
/// @verbatim
/// VertexFromForce 1/2/3(dimension) 1 (no of vertices)
/// Force component(s)
/// 1st vertex index
/// 2nd vertex index
/// ...
/// @endverbatim
/// 
///
class VertexFromForce : public BaseReaction {
  
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
  VertexFromForce(std::vector<double> &paraValue, 
			 std::vector< std::vector<size_t> > &indValue );
  
  ///
  /// @brief Derivative function for this reaction class
  ///
  /// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
  ///
  void derivs(Tissue &T,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      DataMatrix &cellDerivs,
	      DataMatrix &wallDerivs,
	      DataMatrix &vertexDerivs );
};

///
/// @brief Updates position of vertices assuming that they are constrained with a ball from above
///
/// In a model file the reaction is defined as
///
/// @verbatim
/// VertexFromBall 4 0
/// Radius Xc Yc Zc
/// @endverbatim
/// 
/// or
///
/// @verbatim
/// VertexFromBall 7 0
/// Radius Xc Yc Zc dXc dYc dZc
/// @endverbatim
///
/// where radius is the size of the 'ball' pushing at the tissue, Xc,Yc,Zc is the center 
/// of the ball, and the optional dXc,dYc,dZc are the rates for moving the ball along the different
/// directions (the movement is defined in the update function).
/// 
class VertexFromBall : public BaseReaction {
  
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
  VertexFromBall(std::vector<double> &paraValue, 
			 std::vector< std::vector<size_t> > &indValue );
  
  ///
  /// @brief Derivative function for this reaction class
  ///
  /// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
  ///
  void derivs(Tissue &T,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      DataMatrix &cellDerivs,
	      DataMatrix &wallDerivs,
	      DataMatrix &vertexDerivs );

  ///
  /// @brief Update function for this reaction class
  ///
  /// @see BaseReaction::update(Tissue &T,...)
  ///
  void update(Tissue &T,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      double h);
  
};

//-------------------------------------------------


// Do not use this reaction. Restricted area (unless you are a developer).
class DebugReaction : public BaseReaction
{
 public:
  DebugReaction(std::vector<double> &paraValue,
		std::vector< std::vector<size_t> > &indValue);
  
  void derivs(Tissue &T,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      DataMatrix &cellDerivs,
	      DataMatrix &wallDerivs,
	      DataMatrix &vertexDerivs);
};

#endif
