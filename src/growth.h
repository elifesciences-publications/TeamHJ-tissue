//
// Filename     : growth.h
// Description  : Classes describing growth updates
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : April 2006
// Revision     : $Id:$
//
#ifndef GROWTH_H
#define GROWTH_H

#include<cmath>

#include"tissue.h"
#include"baseReaction.h"

///
/// @brief Exponential wall growth which truncates at threshold length
///
/// Exponential growth (proportional to length) truncated at threshold length. 
/// The wall lengths, L, are updated according to
///
/// @f[ \frac{dL}{dt} = p_0 L (1-\frac{L}{p_1}) @f]
///
/// @f$ p_0 @f$ is the growth rate (@f$ k_{growth} @f$), @f$ p_1 @f$ is the 
/// threshold for growth truncation (@f$ L_{trunc} @f$).
/// The column index for the wall length should be given at the first level.
///
/// In a model file the reaction is defined as
///
/// @verbatim
/// WallGrowthExponentialTruncated 2 1 1
/// k_growth L_trunc
/// L
/// @endverbatim
///
class WallGrowthExponentialTruncated : public BaseReaction {
  
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
  WallGrowthExponentialTruncated(std::vector<double> &paraValue, 
				 std::vector< std::vector<size_t> > 
				 &indValue );
  ///
  /// @brief Derivative function for this reaction class
  ///
  /// @see BaseReaction::derivs(Tissue &T,...)
  ///  
  void derivs(Tissue &T,
	      std::vector< std::vector<double> > &cellData,
	      std::vector< std::vector<double> > &wallData,
	      std::vector< std::vector<double> > &vertexData,
	      std::vector< std::vector<double> > &cellDerivs,
	      std::vector< std::vector<double> > &wallDerivs,
	      std::vector< std::vector<double> > &vertexDerivs );
};

///
/// @brief Exponential strech-driven wall growth which truncates at threshold length
///
/// Exponential growth driven by a streched wall and truncated at
/// threshold length. The wall lengths, L, are updated only if the
/// length is shorter than the distance between the vertices of the
/// wall and then according to
///
/// @f[ \frac{dL}{dt} = p_0 (d_v-L) (1 - \frac{L}{p_1}) @f]
///
/// @f$ p_0 @f$ is the growth rate (@f$ k_{growth} @f$), @f$ p_1 @f$ is the 
/// threshold for growth truncation (@f$ L_{trunc} @f$).
/// @f$ d_v @f$ is the distance between the two wall vertices.
///  
/// The column index for the wall length should be given at the first level.
///
/// In a model file the reaction is defined as
///
/// @verbatim
/// WallGrowthExponentialStressTruncated 2 1 1
/// k_growth L_trunc
/// L
/// @endverbatim
///
/// @note Should be used with concern since it uses the stretch/strain (not stress 
/// if spring constant can vary) and always exponential growth (prop to length).
/// Should probably use the more general WallGrowthStress instead, but it is still 
/// here since it has the truncation feature.
///
/// @see WallGrowthStress
///
class WallGrowthExponentialStressTruncated : public BaseReaction {
  
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
  WallGrowthExponentialStressTruncated(std::vector<double> &paraValue, 
				       std::vector< std::vector<size_t> > 
				       &indValue );
  ///
  /// @brief Derivative function for this reaction class
  ///
  /// @see BaseReaction::derivs(Tissue &T,...)
  ///  
  void derivs(Tissue &T,
	      std::vector< std::vector<double> > &cellData,
	      std::vector< std::vector<double> > &wallData,
	      std::vector< std::vector<double> > &vertexData,
	      std::vector< std::vector<double> > &cellDerivs,
	      std::vector< std::vector<double> > &wallDerivs,
	      std::vector< std::vector<double> > &vertexDerivs );
};

///
/// @brief Constant stress/strain-driven wall growth dependent on a threshold
///
/// Constant (constant mode) or exponential (linear mode) growth driven by a 
/// streched/stressed wall. The wall 
/// lengths are updated only if the length plus a threshold value is shorter 
/// than the distance between the vertices of the wall (in strain mode), and if 
/// the total stress is above the threshold (in stress mode). Hence the update follows:
///
///  @f[ \frac{dL}{dt} = p_0 L (S-p_1) @f] if @f$ S > p_1 @f$ and 0 otherwise.
///
/// where S is the stress/strain and L is the wall length.
/// @f$ p_0 @f$ is the growth rate (@f$ k_{growth} @f$). 
/// @f$ p_1 @f$ is a threshold (@f$ s_{threshold} @f$) in stress or strain depending 
/// on @f$ p_2 @f$. 
/// @f$ p_2 @f$ is a flag (@f$ strain_{flag} @f$) for using stretch/strain (@f$p_2=1 @f$) 
/// or stress (@f$p_2=0 @f$). Strain is calculated by @f$ S=(d-L)/L @f$, where d is the distance
/// between the vertices. Stress is read from the second layer of variable indices.  
/// @f$ p_3 @f$ is a flag (@f$ linear_{flag} @f$) for using growth proportional to 
/// wall length (@f$ p_3=1 @f$, as in equation above. If @f$ p_3=0 @f$ the above equation
/// will not include the length L on the right-hand side.
///
///  The column index for the wall length should be given in the first layer, and if stress
///  is to be used, a second layer with calculated wall stress variables has to be given.
///
/// In a model file the reaction is defined as
///
/// @verbatim
/// WallGrowthStress 4 1/2 1 [N]
/// k_growth s_threshold stretch_flag linear_flag 
/// L
/// [stress1 ... stressN]
/// @endverbatim
///
/// If stress is used (streth_flag=0) a second level of wall stresses has to be read.
///
class WallGrowthStress : public BaseReaction {
  
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
  WallGrowthStress(std::vector<double> &paraValue, 
		   std::vector< std::vector<size_t> > 
		   &indValue );
  ///
  /// @brief Derivative function for this reaction class
  ///
  /// @see BaseReaction::derivs(Tissue &T,...)
  ///  
  void derivs(Tissue &T,
	      std::vector< std::vector<double> > &cellData,
	      std::vector< std::vector<double> > &wallData,
	      std::vector< std::vector<double> > &vertexData,
	      std::vector< std::vector<double> > &cellDerivs,
	      std::vector< std::vector<double> > &wallDerivs,
	      std::vector< std::vector<double> > &vertexDerivs );
};

///
/// @brief Constant stress/strain-driven wall growth dependent on a
/// distance to an maximal coordinate (e.g. tip)
///
/// Constant growth driven by a streched wall. The wall lengths, L, are
/// updated only if the length is shorter than the distance between the
///  vertices of the wall and then according to
///
///  dL/dt = p_0*(d_v-L-p_1)*p_2^p_3/(p_2^p_3+d^p_3) iff (d_v-L) > p_1
///
/// p_0 is the growth rate. 
/// p_1 is a stress/strain threshold 
/// p_2 is the K_Hill of the spatial factor
/// p_3 is the n_Hill of the spatial factor
/// p_4 is a flag for using stretch/strain instead of stress
/// p_5 is a flag for using growth proportional to wall length (not constant)
/// d_v is the distance between the two wall vertices.
/// d is the distance between the max value and wall
///
///  In addition, the column index for the wall length, the distance
///  coordinate should be given at first level and stress index in second.
///
class WallGrowthStressSpatial : public BaseReaction {
  
 private:
  
  double Kpow_;
  
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
  WallGrowthStressSpatial(std::vector<double> &paraValue, 
			  std::vector< std::vector<size_t> > 
			  &indValue );
  ///
  /// @brief Derivative function for this reaction class
  ///
  /// @see BaseReaction::derivs(Tissue &T,...)
  ///  
  void derivs(Tissue &T,
	      std::vector< std::vector<double> > &cellData,
	      std::vector< std::vector<double> > &wallData,
	      std::vector< std::vector<double> > &vertexData,
	      std::vector< std::vector<double> > &cellDerivs,
	      std::vector< std::vector<double> > &wallDerivs,
	      std::vector< std::vector<double> > &vertexDerivs );
};

///
/// @brief Constant stress/strain-driven wall growth dependent on a
/// distance to an maximal coordinate (e.g. tip)
///
/// Constant growth driven by a streched wall. The wall lengths, L, are
/// updated only if the length is shorter than the distance between the
///  vertices of the wall and then according to
///
///  dL/dt = p_0*(d_v-L-p_1)*p_2^p_3/(p_2^p_3+d^p_3) iff (d_v-L) > p_1
///
/// p_0 is the growth rate. 
/// p_1 is a stress/strain threshold 
/// p_2 is the K_Hill of the spatial factor
/// p_3 is the n_Hill of the spatial factor
/// p_4 is a flag for using stretch/strain instead of stress
/// p_5 is a flag for using growth proportional to wall length (not constant)
/// d_v is the distance between the two wall vertices.
/// d is the distance between the max value and wall
///
///  In addition, the column index for the wall length, the distance
///  coordinate should be given at first level and stress index in second.
///
class WallGrowthStressSpatialSingle : public BaseReaction {
  
 private:
  
  double Kpow_;
  
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
  WallGrowthStressSpatialSingle(std::vector<double> &paraValue, 
				std::vector< std::vector<size_t> > 
				&indValue );
  ///
  /// @brief Derivative function for this reaction class
  ///
  /// @see BaseReaction::derivs(Tissue &T,...)
  ///  
  void derivs(Tissue &T,
	      std::vector< std::vector<double> > &cellData,
	      std::vector< std::vector<double> > &wallData,
	      std::vector< std::vector<double> > &vertexData,
	      std::vector< std::vector<double> > &cellDerivs,
	      std::vector< std::vector<double> > &wallDerivs,
	      std::vector< std::vector<double> > &vertexDerivs );
};

///
/// @brief Constant stress/strain-driven wall growth dependent on a
/// concentration level in the cell
///
/// Constant growth driven by a streched wall. The wall lengths, L, are
/// updated only if the length is shorter than the distance between the
///  vertices of the wall and then according to
///
///  dL/dt = (p_0+p_1*f(c,p_2,p_3)) * (d_v-L) if (d_v-L) > p_4
///
/// p_0 is the constant growth rate. 
/// p_1 is the maximal added growth rate (V_max in the Hill function)
/// p_2 is the K_Hill
/// p_3 is the n_Hill
/// p_4 is a threshold 
/// p_5 is a flag for using stretch/strain instead of stress
/// p_6 is a flag for using growth proportional to wall length (not constant)
/// f is the hill function (increasing)
/// c is the concentration
/// d_v is the distance between the two wall vertices.
/// 
///
///  In addition, the column index for the wall length and the cell
///  concentration should be given.
///
class WallGrowthStressConcentrationHill : public BaseReaction {
  
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
  WallGrowthStressConcentrationHill(std::vector<double> &paraValue, 
				    std::vector< std::vector<size_t> > 
				    &indValue );
  ///
  /// @brief Derivative function for this reaction class
  ///
  /// @see BaseReaction::derivs(Tissue &T,...)
  ///  
  void derivs(Tissue &T,
	      std::vector< std::vector<double> > &cellData,
	      std::vector< std::vector<double> > &wallData,
	      std::vector< std::vector<double> > &vertexData,
	      std::vector< std::vector<double> > &cellDerivs,
	      std::vector< std::vector<double> > &wallDerivs,
	      std::vector< std::vector<double> > &vertexDerivs );
};

//!Constant strech-driven wall growth 
/*!Constant growth driven by a streched wall. The wall lengths, L, are
  updated only if the length is shorter than the distance between the
  vertices of the wall and then according to
  
  dL/dt = p_0*(d_v-L)*f_e
  
  p_0 is the growth rate. 
  f_e fraction for epidermal walls
  d_v is the distance between the two wall vertices.
  
  In addition, the column index for the wall length should be given.
*/
class WallGrowthConstantStressEpidermalAsymmetric : public BaseReaction {
  
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
  WallGrowthConstantStressEpidermalAsymmetric(std::vector<double> &paraValue, 
					      std::vector< std::vector<size_t> > 
					      &indValue );
  ///
  /// @brief Derivative function for this reaction class
  ///
  /// @see BaseReaction::derivs(Tissue &T,...)
  ///  
  void derivs(Tissue &T,
	      std::vector< std::vector<double> > &cellData,
	      std::vector< std::vector<double> > &wallData,
	      std::vector< std::vector<double> > &vertexData,
	      std::vector< std::vector<double> > &cellDerivs,
	      std::vector< std::vector<double> > &wallDerivs,
	      std::vector< std::vector<double> > &vertexDerivs );
};

///
/// @brief Growth via vertex movement radially outwards
///
/// The tissue grows from vertex movement radially outwards. The update can be
///
/// @f[ \frac{dr}{dt} = p_{0} @f] (if @f$ p_1=0 @f$) or 
/// @f[ \frac{dr}{dt} = p_{0}*r @f] (if @f$ p_{1}=1 @f$) 
///
/// @f$ p_{0} @f$ is the rate (@f$ k_{growth} @f$), 
/// @f$ p_{1} @f$ {0,1} is a flag determining which function to be used (@f$ r_{pow} @f$). 
///
/// In a model file the reaction is defined as
///
/// @verbatim
/// MoveVertexRadially 2 0
/// p_0 p_1
/// @endverbatim
///
class MoveVertexRadially : public BaseReaction {
  
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
  MoveVertexRadially(std::vector<double> &paraValue, 
		     std::vector< std::vector<size_t> > 
		     &indValue );
  ///
  /// @brief Derivative function for this reaction class
  ///
  /// @see BaseReaction::derivs(Tissue &T,...)
  ///  
  void derivs(Tissue &T,
	      std::vector< std::vector<double> > &cellData,
	      std::vector< std::vector<double> > &wallData,
	      std::vector< std::vector<double> > &vertexData,
	      std::vector< std::vector<double> > &cellDerivs,
	      std::vector< std::vector<double> > &wallDerivs,
	      std::vector< std::vector<double> > &vertexDerivs );
};

///
/// @brief Growth via vertex movement along sphereCylinder
///
/// The tissue grows from vertex movement outwards (from apex) along sphereCylinder. 
/// The update can be described by the angular (v) movement
///
/// dv/dt = p_0 (p_1=0) or
/// dv/dt = p_0*r (p_1=1) 
///
/// where 
///
/// v is the angle from the apex
/// r is the sphere radius
/// p_0 is the rate, 
/// p_1 is a flag determining function. 
///
/// On the cylinder the vertex is moved downwards (in -z direction).
/// 
class MoveVertexSphereCylinder : public BaseReaction {
  
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
  MoveVertexSphereCylinder(std::vector<double> &paraValue, 
			   std::vector< std::vector<size_t> > 
			   &indValue );
  ///
  /// @brief Derivative function for this reaction class
  ///
  /// @see BaseReaction::derivs(Tissue &T,...)
  ///  
  void derivs(Tissue &T,
	      std::vector< std::vector<double> > &cellData,
	      std::vector< std::vector<double> > &wallData,
	      std::vector< std::vector<double> > &vertexData,
	      std::vector< std::vector<double> > &cellDerivs,
	      std::vector< std::vector<double> > &wallDerivs,
	      std::vector< std::vector<double> > &vertexDerivs );
};

class WallLengthGrowExperimental : public BaseReaction
{
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
  WallLengthGrowExperimental(std::vector<double> &paraValue,
			     std::vector< std::vector<size_t> > &indValue);
  ///
  /// @brief Derivative function for this reaction class
  ///
  /// @see BaseReaction::derivs(Tissue &T,...)
  ///  
  void derivs(Tissue &T,
	      std::vector< std::vector<double> > &cellData,
	      std::vector< std::vector<double> > &wallData,
	      std::vector< std::vector<double> > &vertexData,
	      std::vector< std::vector<double> > &cellDerivs,
	      std::vector< std::vector<double> > &wallDerivs,
	      std::vector< std::vector<double> > &vertexDerivs);
};

///
/// @brief Updates the water volume variable given osmotic and turgor
/// potentials
///
/// This function uses a constant osmotic potential and calculates the turgor
/// potential to calculate water intake into the cell according to
///
/// @f[ \frac{V_w}{dt} = p_0 A (p_1-p_2T) @f]
///
/// where V_w is the water volume, T is the turgor,p_0 is the rate, p_1 is the
/// osmotic potential and p_2 is an scaling factor. Also p_3=denyShrink_flag
/// and p_4=allowNegTurgor_flag can be set to restrict the behavior.
///
/// The turgor, T, is calculated as T=V_w-V.
///
class WaterVolumeFromTurgor : public BaseReaction
{
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
  WaterVolumeFromTurgor(std::vector<double> &paraValue,
			std::vector< std::vector<size_t> > &indValue);
  ///
  /// @brief Derivative function for this reaction class
  ///
  /// @see BaseReaction::derivs(Tissue &T,...)
  ///  
  void derivs(Tissue &T,
	      std::vector< std::vector<double> > &cellData,
	      std::vector< std::vector<double> > &wallData,
	      std::vector< std::vector<double> > &vertexData,
	      std::vector< std::vector<double> > &cellDerivs,
	      std::vector< std::vector<double> > &wallDerivs,
	      std::vector< std::vector<double> > &vertexDerivs);
};

///
/// @brief Updates 'concentration' variables according to volume changes from
/// derivatives of the vertex positions.
///
/// The dilution of 'concentration' variables, C_i are calculated according to
///
/// @f$ \frac{C_i}{dt} = - \frac{C_i}{V}\frac{dV}{dt} @f$
///
/// where V is the volume, and dV/dt is calculated from vertex position
/// derivatives.
///
/// @note Since this function uses the derivatives of the vertex positions it
/// needs to be applied after all other derivatives applied to the vertices.
///
class DilutionFromVertexDerivs : public BaseReaction
{
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
  DilutionFromVertexDerivs(std::vector<double> &paraValue,
			   std::vector< std::vector<size_t> > &indValue);
  ///
  /// @brief Derivative function for this reaction class
  ///
  /// @see BaseReaction::derivs(Tissue &T,...)
  ///  
  void derivs(Tissue &T,
	      std::vector< std::vector<double> > &cellData,
	      std::vector< std::vector<double> > &wallData,
	      std::vector< std::vector<double> > &vertexData,
	      std::vector< std::vector<double> > &cellDerivs,
	      std::vector< std::vector<double> > &wallDerivs,
	      std::vector< std::vector<double> > &vertexDerivs);
};

#endif
