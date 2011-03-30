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
/// Alternatively if no force save index is supplied the first line
/// can be replaced by 'VertexFromWallTRBS 2 1 1'.
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
	      std::vector< std::vector<double> > &cellData,
	      std::vector< std::vector<double> > &wallData,
	      std::vector< std::vector<double> > &vertexData,
	      std::vector< std::vector<double> > &cellDerivs,
	      std::vector< std::vector<double> > &wallDerivs,
	      std::vector< std::vector<double> > &vertexDerivs );
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
	      std::vector< std::vector<double> > &cellData,
	      std::vector< std::vector<double> > &wallData,
	      std::vector< std::vector<double> > &vertexData,
	      std::vector< std::vector<double> > &cellDerivs,
	      std::vector< std::vector<double> > &wallDerivs,
	      std::vector< std::vector<double> > &vertexDerivs );

  ///
  /// @brief Reaction initiation applied before simulation starts
  ///
  /// @see BaseReaction::initiate(Tissue &T,...)
  ///
  void initiate(Tissue &T,
		std::vector< std::vector<double> > &cellData,
		std::vector< std::vector<double> > &wallData,
		std::vector< std::vector<double> > &vertexData,
		std::vector< std::vector<double> > &cellDerivs,
		std::vector< std::vector<double> > &wallDerivs,
		std::vector< std::vector<double> > &vertexDerivs );  
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
/// VertexFromTRBScenterTriangulationConcentrationHill 5 2 1 1
/// Y_modulus_min Y_modulus_max P_coeff K_hill n_hill
/// L_ij-index
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
	      std::vector< std::vector<double> > &cellData,
	      std::vector< std::vector<double> > &wallData,
	      std::vector< std::vector<double> > &vertexData,
	      std::vector< std::vector<double> > &cellDerivs,
	      std::vector< std::vector<double> > &wallDerivs,
	      std::vector< std::vector<double> > &vertexDerivs );

  ///
  /// @brief Reaction initiation applied before simulation starts
  ///
  /// @see BaseReaction::initiate(Tissue &T,...)
  ///
  void initiate(Tissue &T,
		std::vector< std::vector<double> > &cellData,
		std::vector< std::vector<double> > &wallData,
		std::vector< std::vector<double> > &vertexData,
		std::vector< std::vector<double> > &cellDerivs,
		std::vector< std::vector<double> > &wallDerivs,
		std::vector< std::vector<double> > &vertexDerivs );  
};

#endif
