//
// Filename     : mechanicalTRBS.h
// Description  : Classes describing mechanical updates for wall springs
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : September 2007
// Revision     : $Id:$
//
#ifndef MECHANICALTRBS_H
#define MECHANICALTRBS_H

#include"tissue.h"
#include"baseReaction.h"
#include<cmath>

///
/// @brief Triangular spring model for plates (2D walls) assuming
/// triangular walls/cells.
///
/// The update (in all dimensions) are given by
///
/// @f[ \frac{dx_i}{dt} = ... @f]
///
/// ...
///
/// The theory of the mechanical model comes from H. Delingette,
/// Triangular springs for modelling non-linear membranes, IEEE Trans
/// Vis Comput Graph 14, 329-41 (2008)
///
/// In a model file the reaction is defined as
///
/// @verbatim
/// VertexFromTRBS 2 1 1
/// Y_modulus P_coeff
/// L_ij-index
/// @endverbatim
///
///
class VertexFromTRBS : public BaseReaction {
  
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
  VertexFromTRBS(std::vector<double> &paraValue, 
		       std::vector< std::vector<size_t> > 
		       &indValue );  
  ///
  /// @brief Derivative function for this reaction class
  ///
  /// @see BaseReaction::derivs(Tissue &T,...)
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
/// @brief Triangular spring model for plates (2D walls) assuming
/// triangulation with a central point on the 2D wall/cell.
///
/// The update (in all dimensions) are given by
///
/// @f[ \frac{dx_i}{dt} = ... @f]
///
/// ...
///
/// The theory of the mechanical model comes from H. Delingette,
/// Triangular springs for modelling non-linear membranes, IEEE Trans
/// Vis Comput Graph 14, 329-41 (2008)
///
/// In a model file the reaction is defined as
///
/// @verbatim
/// VertexFromTRBScenterTriangulation 2 2 1 1
/// Y_modulus P_coeff
/// L_ij-index
/// InternalVarStartIndex
/// or
/// VertexFromTRBScenterTriangulation 2 4 1 1 1/0 1/0
/// Y_modulus P_coeff
/// L_ij-index
/// InternalVarStartIndex
/// Optional index for storing strain
/// Optional index for storing stress
/// @endverbatim
///
class VertexFromTRBScenterTriangulation : public BaseReaction {
  
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
  VertexFromTRBScenterTriangulation(std::vector<double> &paraValue, 
				    std::vector< std::vector<size_t> > 
				    &indValue );  
  ///
  /// @brief Derivative function for this reaction class
  ///
  /// @see BaseReaction::derivs(Tissue &T,...)
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
		DataMatrix &vertexDerivs );  
};

///
/// @brief Triangular spring model for plates (2D walls) assuming
/// triangulation with a central point on the 2D wall/cell
/// concentration(auxin)-dependent Young's modulus.
///
/// The update (in all dimensions) are given by
///
/// @f[ \frac{dx_i}{dt} = ... @f]
///
/// ...
///
/// The theory of the mechanical model comes from H. Delingette,
/// Triangular springs for modelling non-linear membranes, IEEE Trans
/// Vis Comput Graph 14, 329-41 (2008)
///
/// In a model file the reaction is defined as
///
/// @verbatim
/// VertexFromTRBScenterTriangulationConcentrationHill 5 2 2 1
/// Y_modulus_min Y_modulus_max P_coeff K_hill n_hill
/// L_ij-index  concentration-index
/// InternalVarStartIndex
/// @endverbatim
///
class VertexFromTRBScenterTriangulationConcentrationHill : public BaseReaction {
  
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
  VertexFromTRBScenterTriangulationConcentrationHill(std::vector<double> &paraValue, 
						     std::vector< std::vector<size_t> > 
						     &indValue );  
  ///
  /// @brief Derivative function for this reaction class
  ///
  /// @see BaseReaction::derivs(Tissue &T,...)
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
		DataMatrix &vertexDerivs );  
};

///
/// @brief Triangular spring model with anisotropy for plates (2D walls) assuming
/// triangular walls/cells.
///
/// The update (in all dimensions) are given by
///
/// @f[ \frac{dx_i}{dt} = ... @f]
///
/// ...
///
/// The theory of the mechanical model comes from H. Delingette,
/// Triangular springs for modelling non-linear membranes, IEEE Trans
/// Vis Comput Graph 14, 329-41 (2008), wich is developed into an
/// anisotropic model.
///
/// In a model file the reaction is defined as
///
/// @verbatim
///
/// VertexFromTRBSMT 10 1 9
/// 
/// Y_matrix 
/// Y_fiber 
/// Poisson_Long
/// Poisson_Trans
/// MF_flag(0/1) 
/// neighborWeight 
/// unusedparameter 
/// plane-strain/stress-flag 
/// MT-angle MT-feedback-flag 
/// 
/// L_ij-index 
/// MT_cellIndex 
/// strainAnisotropy-Index
/// stressAnisotropy-Index 
/// areaRatioIndex 
/// isoEnergyIndex
/// anisoEnergyIndex 
/// YoungL-index
/// MTstress
///
/// or
///
/// VertexFromTRBSMT 10 3 8 0/1/2/3 0/1/2
///
/// Y_matrix
/// Y_fiber
/// Poisson_Long
/// Poisson_Trans
/// MF_flag(0/1)
/// neighborWeight
/// unusedparameter
/// plane-strain/stress-flag
/// MT-angle
/// MT-feedback-flag
/// 
/// L_ij-index
/// MT_cellIndex
/// strainAnisotropy-Index
/// stressAnisotropy-Index
/// areaRatioIndex 
/// isoEnergyIndex
/// anisoEnergyIndex
/// YoungL-index
/// MTstress
///
/// optional index for storing strain(0: no strain,
///                                   1: strain, 
///                                   2: strain/perpendicular strain, 
///                                   3: strain/perpendicular strain/2nd strain)
///
/// optional index for storing stress(0: no stress, 
///                                   1: stress, 
///                                   2: stress/2nd stress)
///
///  @endverbatim
/// In case of storing strain/stress direction/value, in 3(2) dimensions, 
/// strain/stress values will be stored after  (3) components of vectors.  
/// The value for perpendicular strain is maximal strain value.
class VertexFromTRBSMT : public BaseReaction {
  
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
  VertexFromTRBSMT(std::vector<double> &paraValue, 
		   std::vector< std::vector<size_t> > 
		   &indValue );  
  ///
  /// @brief Derivative function for this reaction class
  ///
  /// @see BaseReaction::derivs(Tissue &T,...)
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
/// @brief Triangular spring model for plates (2D walls) assuming
/// triangulation with a central point on the 2D wall/cell.
///
/// The update (in all dimensions) are given by
///
/// @f[ \frac{dx_i}{dt} = ... @f]
///
/// ...
///
/// The theory of the mechanical model comes from H. Delingette,
/// Triangular springs for modelling non-linear membranes, IEEE Trans
/// Vis Comput Graph 14, 329-41 (2008)
///
/// In a model file the reaction is defined as
///
/// @verbatim
/// VertexFromTRBScenterTriangulationMT 11 2 11 1 
/// 
/// Y_matrix 
/// Y_fiber 
/// Poisson_Long
/// Poisson_Trans
/// MF_flag(0/1) 
/// neighborWeight 
/// max_stress(if 0 absolute stress anisotropy is calculated)
/// plane-strain/stress-flag 
/// MT-angle 
/// MT-feedback-flag 
/// unused parameter 
///
/// L_ij-index 
/// MT_cellIndex 
/// strainAnisotropy-Index 
/// stressAnisotropy-Index
/// areaRatioIndex 
/// isoEnergyIndex 
/// anisoEnergyIndex 
/// youngL-index 
/// MTstressIndex 
/// stressTensorIndex 
/// normalVectorIndex
///
/// InternalVarStartIndex
/// 
/// or
/// 
/// VertexFromTRBScenterTriangulationMT 11 4 11 1 0/1/2/3 0/1/2
///
/// Y_matrix 
/// Y_fiber 
/// Poisson_Long  
/// Poisson_Trans 
/// MF_flag(0/1) 
/// neighborWeight 
/// unusedparameter 
/// plane-strain/stress-flag 
/// MT-angle 
/// MT-feedback-flag
/// unused parameter
/// 
/// L_ij-index 
/// MT_cellIndex 
/// strainAnisotropy-Index 
/// stressAnisotropy-Index
/// areaRatioIndex 
/// isoEnergyIndex 
/// anisoEnergyIndex 
/// youngL-index 
/// MTstressIndex 
/// stressTensorIndex 
/// normalVectorIndex
///
/// InternalVarStartIndex
///
/// optional indices for storing strain(0: no strain, 
///                                     1: strain, 
///                                     2: strain/perpendicular strain, 
///                                     3: strain/perpendicular strain/2nd strain)
/// optional indices for storing stress(0: no stress, 
///                                     1: stress, 
///                                     2: stress/2nd stress)
/// @endverbatim
/// In case of storing strain/stress direction/value, in 3(2) dimensions, 
/// strain/stress values will be stored after (3) components of vectors.
/// The value for perpendicular strain is maximal strain value.  

class VertexFromTRBScenterTriangulationMT : public BaseReaction {
  
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
  VertexFromTRBScenterTriangulationMT(std::vector<double> &paraValue, 
				    std::vector< std::vector<size_t> > 
				    &indValue );  
  ///
  /// @brief Derivative function for this reaction class
  ///
  /// @see BaseReaction::derivs(Tissue &T,...)
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
  
// void initiate(Tissue &T,
// 		DataMatrix &cellData,
// 		DataMatrix &wallData,
// 		DataMatrix &vertexData,
// 		DataMatrix &cellDerivs,
// 		DataMatrix &wallDerivs,
// 		DataMatrix &vertexDerivs );  


};



///
/// @brief Triangular spring model for plates (2D walls) assuming
/// triangulation with a central point on the 2D wall/cell.
/// concentration(auxin)-dependent Young's modulus.
///
/// The update (in all dimensions) are given by
///
/// @f[ \frac{dx_i}{dt} = ... @f]
///
/// ...
///
/// The theory of the mechanical model comes from H. Delingette,
/// Triangular springs for modelling non-linear membranes, IEEE Trans
/// Vis Comput Graph 14, 329-41 (2008)
///
/// In a model file the reaction is defined as
///
/// @verbatim
/// VertexFromTRBScenterTriangulationConcentrationHillMT 8 2 3 1
/// Y_modulus_Longitudinal_min Y_modulus_Longitudinal_max P_coeff_Longitudinal 
/// Y_modulus_Transverse_min Y_modulus_Transverse_max P_coeff_Transverse  
/// K_hill n_hill
/// L_ij-index  concentration-index MT_cellIndex
/// InternalVarStartIndex
///or
/// VertexFromTRBScenterTriangulationConcentrationHillMT 8 6 3 1 1/0 1/0 1/0 1/0
/// Y_modulus_Longitudinal_min Y_modulus_Longitudinal_max P_coeff_Longitudinal 
/// Y_modulus_Transverse_min Y_modulus_Transverse_max P_coeff_Transverse  
/// K_hill n_hill
/// L_ij-index  concentration-index MT_cellIndex
/// InternalVarStartIndex
/// optional index for storing strain
/// optional index for storing 2nd strain
/// optional index for storing stress
/// optional index for storing 2nd stress
/// @endverbatim
/// In case of storing strain/stress direction/value, in 3(2) dimensions, 
/// strain/stress values will be stored after  2(3) components of vectors.  

class VertexFromTRBScenterTriangulationConcentrationHillMT : public BaseReaction {
  
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
  VertexFromTRBScenterTriangulationConcentrationHillMT(std::vector<double> &paraValue, 
						     std::vector< std::vector<size_t> > 
						     &indValue );  
  ///
  /// @brief Derivative function for this reaction class
  ///
  /// @see BaseReaction::derivs(Tissue &T,...)
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
		DataMatrix &vertexDerivs );  
};





///
/// @brief 
/// Updates Young modulus of cells within "update" based on linear or nonlinear Fiber_model
/// It uses anisotropy(stress or strain) that should be calculated by VertexFromTRBS... functions.
/// and velocity that should be calculated by "UpdateMTdirectionEquilibrium" function. It cooperates 
/// with TRBS functions via cell vector component wich is introduced for storing longitudinal Young modulus.
/// 
/// In a model file the reaction is defined as
///
/// @verbatim
///
/// FiberModel 8 3 1 1 1
/// 
///  k_rate
///  velocity threshold
///  liniear-hill-flag 0/1
///  k_hill
///  n_hill
///  Y_matrix
///  Y_fiber
///  initiate_flag (0:no initiation ,1: initiate with isotropic ,2: initiate with anisotropy from aniso_index)
///
///  anisotropy index.
///
///  Young_Longitudinal index
///
///  store index for velocity from "UpdateMTdirectionEquilibrium"
///
/// @endverbatim
///


class FiberModel : public BaseReaction {
  
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
  FiberModel(std::vector<double> &paraValue, 
             std::vector< std::vector<size_t> > 
             &indValue );  
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
		DataMatrix &vertexDerivs );  
  ///
  /// @brief Derivative function for this reaction class
  ///
  /// @see BaseReaction::derivs(Tissue &T,...)
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


#endif
