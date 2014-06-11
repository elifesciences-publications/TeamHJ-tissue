//
// Filename     : creation.h
// Description  : Classes describing molecular production/creation updates
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : January 2011
// Revision     : $Id:$
//
#ifndef CREATION_H
#define CREATION_H

#include"tissue.h"
#include"baseReaction.h"
#include<cmath>

///
/// @brief In each cell a molecule is produced/created with a constant rate
///
/// The variable update is for each cell given by 
///
/// @f[ \frac{dc}{dt} = k_c @f]
///
/// where @f$ k_c @f$ is a constant parameter, and @f$ c @f$ is the variable to be updated.
///
/// In a model file the reaction is defined as
///
/// @verbatim
/// creationZero 1 1 1
/// k_c
/// c_index
/// @endverbatim
///
class CreationZero : public BaseReaction {
  
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
  CreationZero(std::vector<double> &paraValue, 
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
/// @brief In each cell a molecule is produced/created with a rate dependent on another molecule.
///
/// The variable update is for each cell given by 
///
/// @f[ \frac{dc}{dt} = - k_c X @f]
///
/// where @f$ k_c @f$ is a constant parameter, @f$ c @f$ is the variable to be updated,
/// and @f$ X @f$ is the concentration of the production-dependent molecule.
///
/// In a model file the reaction is defined as
///
/// @verbatim
/// creationZero 1 2 1 1
/// k_c
/// c_index
/// X_index
/// @endverbatim
///
class CreationOne : public BaseReaction {
  
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
  CreationOne(std::vector<double> &paraValue, 
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

class CreationTwo : public BaseReaction {
  
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
  CreationTwo(std::vector<double> &paraValue, 
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
/// @brief In each cell a molecule is produced/created with rate dependent on the distance of the cell from the center
///
/// The variable update is for each cell given by ( SIGN= -1, production inside the sphere)
///
/// @f[ \frac{dc}{dt} = V \frac{r^n + R^n}{R^n} @f]
///
/// or (SIGN = +1, production outside the sphere),
///
/// @f[ \frac{dc}{dt} = V \frac{r^n + R^n}{R^n} @f]
///
/// where @f$ V, R, n, SIGN@f$ are constant parameters, @f$ c @f$ is the variable to be updated and @f$ r @f$ the distance of the cell to the center of the template.
///
/// In a model file the reaction is defined as
///
/// @verbatim
/// creationSpatialSphere 4 1 1
/// V R n SIGN
/// c_index
/// @endverbatim
///
class CreationSpatialSphere : public BaseReaction {
  
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
  CreationSpatialSphere(std::vector<double> &paraValue, 
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


