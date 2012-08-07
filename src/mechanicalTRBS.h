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
/// VertexFromTRBSMT 5 1 3
/// Y_modulus_Matrix Y_modulus_Fibre Poisson_Longit. P_coeff_Trans. Matrix-Fiber-flag 
/// L_ij-index MT_cellIndex anisotropyIndex
/// or
/// VertexFromTRBSMT 5 3 3 0/1/2/3 0/1/2
/// Y_modulus_Matrix Y_modulus_Fibre Poisson_Longit. P_coeff_Trans. Matrix-Fiber-flag 
/// L_ij-index MT_cellIndex anisotropyIndex
/// optional index(indices) for storing strain(strain(1), perpendicular to strain(2) and 2nd strain(3))
/// optional index(indices) for storing stress(stress(1) and 2nd stress(2))
/// @endverbatim
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
/// VertexFromTRBScenterTriangulationMT 5 2 3 1  
/// Y_matrix Y_fiber Poisson_Long  Poisson_Trans MF_flag 
/// L_ij-index MT_cellIndex Anisotropy-Index
/// InternalVarStartIndex
/// or
/// VertexFromTRBScenterTriangulationMT 5 4 3 1 0/1/2/3 0/1/2
/// Y_matrix Y_fiber Poisson_Long  Poisson_Trans MF_flag  
/// L_ij-index MT_cellIndex Anisotropy-Index
/// InternalVarStartIndex
/// optional index for storing strain(0: no strain, 1: strain, 2: strain/perpendicular strain, 3: strain/perpendicular strain/2nd strain)
/// optional index for storing stress(0: no stress, 1: stress, 2: stress/2nd stress)
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

#endif
