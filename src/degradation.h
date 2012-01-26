//
// Filename     : degradation.h
// Description  : Classes describing molecular production/degradation updates
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : January 2011
// Revision     : $Id:$
//
#ifndef DEGRADATION_H
#define DEGRADATION_H

#include"tissue.h"
#include"baseReaction.h"
#include<cmath>

///
/// @brief In each cell a molecule is degraded with a constant rate (dependent on its own conc)
///
/// The variable update is for each cell given by 
///
/// @f[ \frac{dc}{dt} = - k_c c@f]
///
/// where @f$ k_c @f$ is a constant parameter and @f$ c @f$ is the (mulecular) variable/concentration
/// to be updated.
///
/// In a model file the reaction is defined as
///
/// @verbatim
/// degradationOne 1 1 1
/// k_c
/// c_index
/// @endverbatim
///
class DegradationOne : public BaseReaction {
  
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
  DegradationOne(std::vector<double> &paraValue, 
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
/// @brief In each cell a molecule is degraded with a rate dependent on another molecule.
///
/// The variable update is for each cell given by 
///
/// @f[ \frac{dc}{dt} = - k_c c X @f]
///
/// where @f$ k_c @f$ is a constant parameter, @f$ c @f$ is the variable to be updated,
/// and @f$ X @f$ is the concentration of the 
/// production-dependent molecule.
///
/// In a model file the reaction is defined as
///
/// @verbatim
/// degradationTwo 1 2 1 1
/// k_c
/// c_index
/// X_index
/// @endverbatim
///
class DegradationTwo : public BaseReaction {
  
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
  DegradationTwo(std::vector<double> &paraValue, 
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

#endif
