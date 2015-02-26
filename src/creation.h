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
  ///
  /// @brief Derivative function for this reaction class calculating the absolute value for noise solvers
  ///
  /// @see BaseReaction::derivsWithAbs(Compartment &compartment,size_t species,...)
  ///
  void derivsWithAbs(Tissue &T,
		     DataMatrix &cellData,
		     DataMatrix &wallData,
		     DataMatrix &vertexData,
		     DataMatrix &cellDerivs,
		     DataMatrix &wallDerivs,
		     DataMatrix &vertexDerivs,
		     DataMatrix &sdydtCell,
		     DataMatrix &sdydtWall,
		     DataMatrix &sdydtVertex );
  
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
/// creationOne 1 2 1 1
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
  ///
  /// @brief Derivative function for this reaction class calculating the absolute value for noise solvers
  ///
  /// @see BaseReaction::derivsWithAbs(Compartment &compartment,size_t species,...)
  ///
  void derivsWithAbs(Tissue &T,
		     DataMatrix &cellData,
		     DataMatrix &wallData,
		     DataMatrix &vertexData,
		     DataMatrix &cellDerivs,
		     DataMatrix &wallDerivs,
		     DataMatrix &vertexDerivs,
		     DataMatrix &sdydtCell,
		     DataMatrix &sdydtWall,
		     DataMatrix &sdydtVertex );
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


///
/// @brief In each cell a molecule is produced/created with rate dependent on the distance of the cell from a ring
///
/// The variable update is for each cell given by ( SIGN= -1, production inside the ring)
///
/// @f[ \frac{dc}{dt} = V \frac{K^n}{K^n + d^n} @f]
///
/// or (SIGN = +1, production away from the ring),
///
/// @f[ \frac{dc}{dt} = V \frac{d^n}{K^n + d^n} @f]
///
/// where @f$ V, K, n, SIGN@f$ are constant parameters, @f$ c @f$ is the variable to be updated and @f$ d @f$ the distance of the cell to ring:
///
/// @f[ d = |r_{cell} - r_{ring}| @f]
///
/// In the above definition, @f$ r_{cell} @f$ is the distance of the cell to the center of the template, and @f$ r_{ring} @f$ is a constant parameter
///
/// In a model file the reaction is defined as
///
/// @verbatim
/// creationSpatialRing 5 1 1
/// V R r_ring n SIGN
/// c_index
/// @endverbatim
///
class CreationSpatialRing : public BaseReaction {
  
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
  CreationSpatialRing(std::vector<double> &paraValue, 
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
/// @brief In each cell a molecule is produced/created with rate dependent on one of
/// the coordinates of each cell
///
/// The variable update is for each cell given by ( SIGN= -1, higher production for
/// lower values of the coordinate)
///
/// @f[ \frac{dc}{dt} = V \frac{X^n}{X^n + x^n} @f]
///
/// or (SIGN = +1, higher production at larger coordinate values),
///
/// @f[ \frac{dc}{dt} = V \frac{x^n}{x^n + X^n} @f]
///
/// where @f$ V, X, n, SIGN@f$ are constant parameters, @f$ c @f$ is the variable to
/// be updated and @f$ x @f$ is the cell coordinate used to set the rate-dependent
/// production.
///
/// In a model file the reaction is defined as
///
/// @verbatim
/// creationSpatialSphere 4 2 1 1
/// V X n SIGN
/// c_index
/// x_index
/// @endverbatim
///
class CreationSpatialCoordinate: public BaseReaction {
  
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
  CreationSpatialCoordinate(std::vector<double> &paraValue, 
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
/// @brief In the cells with given indices a molecule is produced/created with constant rate.
///
/// The variable update is for each cell given by ( SIGN= -1, higher production for
/// lower values of the coordinate)
///
/// @f[ \frac{dc}{dt} = k_c @f] if cell index is in the given list
///
///
/// In a model file the reaction is defined as
///
/// @verbatim
/// creationFromList 1 2 1 n
/// k_c
/// c_index
/// a list of indices with n members
/// @endverbatim
///
class CreationFromList: public BaseReaction {
  
  
private: 
  size_t proCells;

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
  CreationFromList(std::vector<double> &paraValue, 
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


