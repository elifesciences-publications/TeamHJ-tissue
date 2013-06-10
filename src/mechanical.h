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
/// @brief Updates vertices from a cell pressure potential, i.e. forces normal to edges
///
/// A area rule in two dimensions is used to calculate the forces on vertex @$v@$ from a cell is
///
/// @f[\frac{dx_v}{dt} = 0.5*p_0 (y_{v_r} - y_{v_l}) @f]
/// @f[\frac{dy_v}{dt} = 0.5*p_0 (x_{v_l} - x_{v_r}) @f]
///
/// where @$v_r,v_l@$ are right and left vertices in the sorted order. @$p_0$ represents the pressure,
/// and if @$p_1=1@$, the pressure will be divided by the cell volume.
///
/// In a model file, the reaction is given by
///
/// @verbatim
/// VertexFromCellPressure 2 0
/// P V_normflag(=0/1)
/// @endverbatim
///
/// @note Requires two dimensions with vertices sorted.
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

namespace CenterTriangulation {
  ///
  /// @brief Updates vertices from a cell pressure potential
  ///
  /// This function determines the direction of the pressure force term
  /// from the position of the central mesh cell vertex to the center of the wall. 
  /// Applies a force proportional to the pressure (parameter(0)) and the 
  /// size of the wall. parameter(1) equal to 1 normalizes the force with cell volume. 
  /// (0 otherwise).
  ///
  /// In a model file, the reaction is given by
  ///
  /// @verbatim
  /// CenterTriangulation:VertexFromCellPressure 2 1 1
  /// P V_normflag(=0/1)
  /// startIndex 
  /// @endverbatim
  ///
  /// where the startindex is marking the start of internal edge varibales (x,y,z,L_1,...).
  ///
  /// @note Maybe it should rather be normal to the wall in the plane of the triangle?
  /// @note Assumes three dimensions as all CenterTriangulation functions.
  /// @see VertexFromCellPressure
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
  /// @brief Updates vertices from a cell pressure potential linearly increasing in a given time span
  ///
  /// This function determines the direction of the pressure force term
  /// from the position of the central mesh cell vertex to the center of the wall. 
  /// Applies a force proportional to the pressure (parameter(0)) and the 
  /// size of the wall.
  ///
  /// In a model file the reaction is defined as
  ///
  /// @verbatim
  /// CenterTriangulation:VertexFromCellPressureLinear 3 1 1
  /// P A_flag deltaT
  /// @endverbatim
  /// @note Maybe it should rather be normal to the wall in the plane of the triangle?
  ///
  class VertexFromCellPressureLinear : public BaseReaction {
  private:
    
    double timeFactor_; 
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
    VertexFromCellPressureLinear(std::vector<double> &paraValue, 
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
} // end namespace CenterTriangulation


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


///
/// @brief 
///
/// In a model file the reaction is defined as
///
/// @verbatim
/// VertexFromPressureExperimental
///
///  
/// 
/// @endverbatim
///
///
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


///
/// @brief 
///
/// In a model file the reaction is defined as
///
/// @verbatim
/// CellVolumeExperimental 4 2 2 n
///
/// k_p  P_max  k_pp allowShrink_flag
///
/// Wall_length_index  cell_volume_index
/// Force indices
///
/// or
///
/// CellVolumeExperimental 4 3 2 n 1
///
/// k_p  P_max  k_pp allowShrink_flag
///
/// Wall_length_index  cell_volume_index
/// Force indices
/// Optionally_index_for_saving_the_pressure

/// @endverbatim
///
///
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


///
/// @brief Updates vertices from a 'pressure' term defined to act in the cell normal
/// direction. The pressure is applied increasingly (... linear in a given time span).
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
/// VertexFromCellPlaneLinear 3 0
/// P A_flag deltaT
/// @endverbatim
///
/// @see CalculatePCAPlane
///
class VertexFromCellPlaneLinear : public BaseReaction
{
private:
  
  double timeFactor_;
  
public:
  VertexFromCellPlaneLinear(std::vector<double> &paraValue,
                            std::vector< std::vector<size_t> > &indValue);
  
  void derivs(Tissue &T,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      DataMatrix &cellDerivs,
	      DataMatrix &wallDerivs,
	      DataMatrix &vertexDerivs);
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


///
/// @brief Updates vertices from a 'pressure' term defined to act in the normal 
/// direction  to  triangular  elements  of  the cell force for each triangular 
/// element is calculated according to  the  element's  current area and in the 
/// direction of normal  to each  triangular  element. The force is distributed 
/// equally on  nodes (including the centeral node). The  pressure  is  applied 
/// increasingly (... linear in a given time span).
/// It does  not rely  on  PCA  plane  in contrast with VertexFromCellPlane and 
/// VertexFromCellPlaneLinear .
///
/// A cell contributes to a vertex update with
///
/// @f[ \frac{dx_{i}}{dt} = p_{0} A n_{i} / 3 @f]
///
/// where @f$p_{0}@f$ is a 'pressure' parameter, A is the triangular element area @f$n_{i}@f$ is the 
/// triangular element normal vector .
/// An additional parameter @f$p_{2}@f$ can be used to not include the area factor if
/// set to zero (normally it should be set to 1).
///
/// In a model file the reaction is defined as
///
/// @verbatim
///
/// VertexFromCellPlaneLinearCenterTriangulation 3 1 1
///
/// P 
/// Area_flag
/// deltaT
///
/// InternalVarStartIndex
///
/// @endverbatim
///
/// @see CalculatePCAPlane
///
class VertexFromCellPlaneLinearCenterTriangulation : public BaseReaction
{
private:
  
  double timeFactor_;
  
public:
  VertexFromCellPlaneLinearCenterTriangulation(std::vector<double> &paraValue,
                            std::vector< std::vector<size_t> > &indValue);
  
  void derivs(Tissue &T,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      DataMatrix &cellDerivs,
	      DataMatrix &wallDerivs,
	      DataMatrix &vertexDerivs);
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





///
/// @brief Updates vertices from a 'pressure' term defined to act in the cell normal
/// direction.
///
/// This function calculates the area of a cell and then distribute a force 'outwards'
/// among the cell vertices. It relies on that the PCA cell planes have been calculated.
/// A cell contributes to a  ...............


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
/// as an indication for equilibrium state in case of elastic deformations the function calculates the volume between template and Z=Z0 plane.this volume is not used in calculations so Z0 value can be choosen arbitrarily. 
///
/// In a model file the reaction is defined as
///
/// @verbatim
/// VertexFromCellPlaneTriangular 3 0
/// P A_flag Z0
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
/// @brief Updates list of vertices with a given force applied where the force is 
/// linearly increased from zero across a given time span (deltaT)
///
/// In a model file the reaction is defined as
///
/// @verbatim
/// VertexFromForceLinear 1/2/3(dimension+1) 1 (no of vertices)
/// Force component(s) deltaT
/// 1st vertex index
/// 2nd vertex index
/// ...
/// @endverbatim
/// 
///
class VertexFromForceLinear : public BaseReaction {
  
 private:

  double timeFactor_;

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
  VertexFromForceLinear(std::vector<double> &paraValue, 
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

///
/// @brief Updates position of vertices assuming that a ball is moving with a given velocity vector into the ball
///The force applied outward respect to ball proportional to (overlap)^(3/2)
///
/// In a model file the reaction is defined as
///
/// @verbatim
/// VertexFromBall 5 0
/// Radius Xc Yc Zc Kforce
/// @endverbatim
/// 
/// or
///
/// @verbatim
/// VertexFromBall 8 0
/// Radius Xc Yc Zc Kforce dXc dYc dZc
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


///
/// @brief Updates position of vertices assuming that a parabolid is moving with a given velocity(z) into the template
///The force applied outward respect to the parabolid
///
/// In a model file the reaction is defined as
///
/// @verbatim
/// VertexFromParabolid 5 0
/// a Xc Yc b Kforce
/// 
/// or
///
/// VertexFromParabolid 8 0
/// a Xc Yc b Kforce VelocityZ
/// @endverbatim
///
/// where the parabolid is defined by z=a((x-xc)2 +(y-yc)2)+b radius is the size of the 'ball' pushing at the tissue, Xc,Yc,Zc is the center 
/// of the ball, and the optional dXc,dYc,dZc are the rates for moving the ball along the different
/// directions (the movement is defined in the update function).
/// 
class VertexFromParabolid : public BaseReaction {
  
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
  VertexFromParabolid(std::vector<double> &paraValue, 
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


///
/// @brief Updates position of vertices assuming that an external wall is moving with a given velocity vector 
/// toward the meristem
/// The force applied outward respect to wall proportional to (overlap)^(3/2)
///
/// In a model file the reaction is defined as

/// @verbatim
/// VertexFromExternalWall 12 0
/// X0 Y0 Z0 nx ny nz Zmin Zmax dXc dYc dZc Kforce
/// @endverbatim
///
/// where n is the normal vector to the 'wall' pushing at the tissue, X0,Y0,Z0 is a point on the wall
/// and dXc,dYc,dZc are the rates for moving the wall along the different
/// directions (the movement is defined in the update function).
/// 
class VertexFromExternalWall : public BaseReaction {
  
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
  VertexFromExternalWall(std::vector<double> &paraValue, 
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



///
/// @brief Calculates change in template volume and its time derivative 
/// and total Derivative and stores them in the given indices in cellData vector
///
/// In a model file the reaction is defined as
///
/// @verbatim
/// TemplateVolumeChange 0 1 6
/// cell-index-VolumeChange       component-index-VolumeChange
/// cell-index-deltaVolumeChange  component-index-deltaVolumeChange
/// cell-index-totalDerivative    component-index-totalDerivative
///
/// @endverbatim
/// 
class TemplateVolumeChange : public BaseReaction 
{

 private:

  DataMatrix vertexDataRest;
  double VolumeChange;
  double deltaVolumeChange;
  double totalDerivative;
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
  TemplateVolumeChange(std::vector<double> &paraValue, 
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
  /// @brief Reaction initiation applied before simulation starts
  ///
  /// @see BaseReaction::initiate(Tissue &T,...)
  ///
  void initiate(Tissue &T,
                DataMatrix &cellData,
                DataMatrix &wallData,
                DataMatrix &vertexData,
                DataMatrix &cellDerivs,
                DataMatrix &wallDerivs,
                DataMatrix &vertexDerivs);
                 
  ///
  /// @brief Update function for this reaction class
  ///
  /// @see BaseReaction::update(Tissue &T,...)
  ///
  void update(Tissue &T,
              DataMatrix &cellData,
	      DataMatrix &vertexData,
              DataMatrix &wallData,
              DataMatrix &vertexDerivs,
              double h);  
};





///
/// @brief Calculates abs(cos(...)) of angle between two 3d vectors
/// (starting from given indices) in cellData vector and stores it in the given 
/// index in cellData vector, uses no parameter 
///
/// In a model file the reaction is defined as
///
/// @verbatim
/// CalculateAngleVectors 0 2 2 1
/// start-index(1st vector)   start-index(2nd vector) 
/// store-index(angle-deg) 
///
/// @endverbatim
/// 
class CalculateAngleVectors : public BaseReaction 
{

 private:

  DataMatrix vertexDataRest;
  double VolumeChange;
  double deltaVolumeChange;
  double totalDerivative;
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
  CalculateAngleVectors(std::vector<double> &paraValue, 
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
/// @brief Calculates the angle between a 3d vector
/// (starting from given indices) in cellData vector and a given axes(x,y,z) 
/// and stores it in the given index in cellData vector, uses one parameter for specifying the axes 
///
/// In a model file the reaction is defined as
///
/// @verbatim
/// AngleVector 1 2 1 1
/// axes_flag (0:X, 1:Y, 2:Z)
/// start-index(the vector)   
/// store-index(angle-deg) 
///
/// @endverbatim
/// 
class AngleVector : public BaseReaction 
{

 private:

  DataMatrix vertexDataRest;
  double VolumeChange;
  double deltaVolumeChange;
  double totalDerivative;
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
  AngleVector(std::vector<double> &paraValue, 
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
