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
/// where @f$ k_c @f$ is a constant parameter and @f$ c @f$ is the (mulecular) cell variable/concentration
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
/// @brief In each cell a molecule is degraded with a rate dependent on another molecule.
///
/// The variable update is for each cell is given by 
///
/// @f[ \frac{dc}{dt} = - k_c c X @f]
///
/// where @f$ k_c @f$ is a constant parameter, @f$ c @f$ is the cell variable to be updated,
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

///
/// @brief In each cell a molecule is degraded with a rate dependent on its own conc and N other variables
///
/// The variable update is for each cell given by 
///
/// @f[ \frac{dc}{dt} = - k_c X_0 ... X_n c@f]
///
/// where @f$ k_c @f$ is a constant parameter, @f$ c @f$ is the (molecular) cell variable/concentration
/// to be updated and the @f$ X_i @f$ are @f$ N=(n+1) @f$ variables on which the degradation rate also depends on
///
/// In a model file the reaction is defined as
///
/// @verbatim
/// degradationN 1 2 1 N
/// k_c
/// c_index
/// X_0 .. X_n
/// @endverbatim
///
class DegradationN : public BaseReaction {
  
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
  DegradationN(std::vector<double> &paraValue, 
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
/// @brief Implements a degradation proportional to its own concentration and
/// a Hill function of another alternatively a Hill function of its own concentration.
///
/// The update to the variable is
///
/// @f[ \frac{dc}{dt} -= k_{d} c \frac{C_{1}^{n}}{(C_{1}^{n}+K^{n})} @f]
///
/// or
///
/// @f[ \frac{dc}{dt} -= k_{d} \frac{c^{n}}{(c^{n}+K^{n})} @f]
///
/// where @f$ k_d @f$ is the degradation rate, and @f$ C_1 @f$ is the
/// user supplied variable. @f$c@f$ is the concentration of the degraded variable.
///
/// In a model file the reaction is defined as
///
/// @verbatim
/// degradationHill 3 2 1 1
/// k_d K n
/// index1
/// index2
/// @endverbatim
///
/// where index1 is the degraded molecule and index2 is the other (e.g miRNA) molecule in the Hill.
///
class DegradationHill : public BaseReaction {
  
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
  DegradationHill(std::vector<double> &paraValue, 
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
/// @brief Implements a degradation proportional to its own concentration and
/// a Hill function that depends on the concentration of N other molecules.
///
/// @details This reaction is a Michaelis-Menten inspired formalism for degradation
/// where the equation is given by
///
/// @f[ \frac{dy_{ij}}{dt} = - p0 \frac{y_{ik}^{p_2}}{p_1^{p_2}+y_{ik}^{p_2}}...
/// \frac{p_{1'}^{p_{2'}}}{p_{1'}^{p_{2'}}+y_{ik'}^{p_{2'}}} y_{ij} @f] 
///
/// where p_0 is the degradation rate (@f$d@f$), and @f$p_1@f$ is the Hill
/// constant (@f$K_{half}@f$), and @f$p_2@f$ is the Hill coefficient (n). The
/// k index is given in the first level of varIndex and corresponds to the
/// activators (increasing the degradation rate), and the k' in the second level
/// which corresponds to the repressors (decreasing the degradation rate). In both
/// layers each index corresponds to a pair of parameters
/// (K,n) that are preceded by the @f$d@f$ parameter.
/// Note that the degradation rate is also once multiplied with the concentration of
/// the affected variable (see equation above).
///
/// In a model file the reaction is, with @f$N@f$ being the number of activators,
/// @f$M@f$ being the number of repressors, and @f$n_p = 1+2(N+M)@f$ is the number of parameters
///
/// @verbatim
/// degradationHillN n_p 3 1 N M
/// d                         # degradation rate
/// K_A0 n_A0 .. K_AN n_AN    # K_half and n parameters for the activators
/// K_R0 n_R0 .. K_RM n_RM    # K_half and n parameters for the repressors
/// c                         # index of degraded molecule
/// k_0 .. k_N                # indices of activators
/// k'_0 .. k'_M              # indices of repressors
/// @endverbatim
///
class DegradationHillN : public BaseReaction {
  
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
  DegradationHillN(std::vector<double> &paraValue, 
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
/// @brief In each cell a molecule is degraded with a rate dependent on another molecule and cell volume
///
/// The variable update is for each cell given by 
///
/// @f[ \frac{dc}{dt} = -volume* k_c c X @f]
///
/// where @f$ k_c @f$ is a constant parameter, @f$ c @f$ is the variable to be updated,
/// and @f$ X @f$ is the concentration of the 
/// production-dependent molecule.
///
/// In a model file the reaction is defined as
///
/// @verbatim
/// degradationTwoGeometric 1 2 1 1
/// k_c
/// c_index
/// X_index
/// @endverbatim
///
class DegradationTwoGeometric : public BaseReaction {
  
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
  DegradationTwoGeometric(std::vector<double> &paraValue, 
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
/// @brief In each cell a wall molecule is degraded with a constant rate (dependent on its own conc)
///
/// The variable update is for each wall variable given by 
///
/// @f[ \frac{dc_w}{dt} = - k_{cw} c_w@f]
///
/// where @f$ k_{cw} @f$ is a constant parameter and @f$ c_w @f$ is the (mulecular) variable/concentration
/// to be updated.
///
/// In a model file the reaction is defined as
///
/// @verbatim
/// degradationOneWall 1 1 1
/// k_cw
/// cw_index
/// @endverbatim
///
class DegradationOneWall : public BaseReaction {
  
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
  DegradationOneWall(std::vector<double> &paraValue, 
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
/// @brief In each boundary cell a molecule is degraded with a constant rate (dependent on its own conc)
///
/// The variable update is for each boundary cell given by 
///
/// @f[ \frac{dc}{dt} = - k_c c@f]
///
/// where @f$ k_c @f$ is a constant parameter and @f$ c @f$ is the (mulecular) cell variable/concentration
/// to be updated. Boundary cells are cells lies at the boundary of the template, i.e., has the background
/// as a neighbor
///
/// In a model file the reaction is defined as
///
/// @verbatim
/// degradationOneBoundary 1 1 1
/// k_c
/// c_index
/// @endverbatim
///
class DegradationOneBoundary : public BaseReaction {
  
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
  DegradationOneBoundary(std::vector<double> &paraValue, 
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
/// @brief In the cells with given indices a molecule is degraded with a constant rate (dependent on its own conc)
///
/// The variable update is for each cell given by 
///
/// @f[ \frac{dc}{dt} = - k_c c@f]
///
/// where @f$ k_c @f$ is a constant parameter and @f$ c @f$ is the (mulecular) cell variable/concentration
/// to be updated.
///
/// In a model file the reaction is defined as
///
/// @verbatim
/// degradationOneFromList 1 2 1 N
/// k_c
/// c_index
/// cell_index_0 .. cell_index_N
/// @endverbatim
///
class DegradationOneFromList : public BaseReaction {

 private:
  size_t numCellI;
  
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
  DegradationOneFromList(std::vector<double> &paraValue, 
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

#endif
