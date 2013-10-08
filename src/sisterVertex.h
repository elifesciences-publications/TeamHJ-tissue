//
// Filename     : sisterVertex.h
// Description  : Classes describing reactions related to sisterVertices
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : July 2013
// Revision     : $Id:$
//
#ifndef SISTERVERTEX_H
#define SISTERVERTEX_H

#include<cmath>

#include"tissue.h"
#include"baseReaction.h"

///
/// @brief SisterVertex is used to connect vertices not connected by edges
///
/// SisterVertex can be seen as an additional way to connect certices except by normal wall edges. 
/// This can be used to connect cell walls that are allowed to break etc. Also,
/// the SisterVertex can be used to connect two 3D cells.
///
namespace SisterVertex {
  
  /// 
  /// @brief Initiates sister vertices by reading a list of pair indices from a file
  ///
  /// This reaction initiates the sister vertices by reading pair indices from a file 
  /// named 'sister' in the current directory. It uses no parameters or variableIndices
  /// and is called from a model file as
  ///
  /// @verbatim
  /// SisterVertex::InitiateFromFile 0 0
  /// @endverbatim
  ///
  /// It expects a file with a first line with number of sister pairs followed by
  /// pairs of vertex indices that are to become sisters:
  /// @verbatim
  /// N_pairs
  /// vertexIndex1 vertexIndex2
  /// ...
  /// @endverbatim
  ///
  /// After the initiation it
  /// does not update the list.
  ///
  /// @note The list of vertices must be stored in a file called sister
  ///
  class InitiateFromFile : public BaseReaction {
    
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
    InitiateFromFile(std::vector<double> &paraValue, 
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
		  DataMatrix &vertexDerivs);
    
    ///
    /// @brief Derivative function for this reaction class
    ///
    /// For this reaction nothing is added to the derivative
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
  /// @brief Initiates sister vertices by adding vertex pairs that are close in space
  ///
  /// This reaction will go through all vertices and add all closer than a distance, @f$d_{max}@f$
  /// , from each other in the sisterVertex list. It will only initiate the list and not
  /// do any further updates. In a model file it is defined as
  /// @verbatim
  /// SisterVertex::InitiateFromDistance 1 0
  /// d_{max}
  /// @endverbatim
  ///
  class InitiateFromDistance : public BaseReaction {
    
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
    InitiateFromDistance(std::vector<double> &paraValue, 
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
		  DataMatrix &vertexDerivs);    
    ///
    /// @brief Derivative function for this reaction class
    ///
    /// For this reaction nothing is added to the derivative
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
  /// @brief A mechanical spring between all defined sister vertices
  ///
  /// A simple spring force is added between all sister vertices, and their spatial
  /// directions are updated according to
  ///
  /// @f[ \frac{dx_i}{dt} = - k_{spring} (x_{i}-x_{j}) @f]
  /// 
  /// where @f$i,j@f$ are the two sister vertices and it is done in all 
  /// spatial directions (assuming resting length 0).
  /// In a model file the reaction is given as
  /// 
  /// @verbatim
  /// SisterVertex::Spring 1 0
  /// k_{spring}
  /// @endverbatim
  ///
  class Spring : public BaseReaction {
    
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
    Spring(std::vector<double> &paraValue, 
	   std::vector< std::vector<size_t> > 
	   &indValue );
    
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
  /// @brief Combines (adds) the derivatives for two sister vertices
  ///
  /// This reaction adds up all derivative contribution for one sister node to the other
  /// and vice versa. This leads to that they always move in concert.
  ///
  /// In a model file the reaction is given by
  ///
  /// @verbatim
  /// SisterVertex::CombineDerivatives 0 0
  /// @endverbatim
  ///
  /// @note Since the derivative values are added directly, this reaction must be placed
  /// after all reactions updating vertex positions.
  ///
  class CombineDerivatives : public BaseReaction {
    
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
    CombineDerivatives(std::vector<double> &paraValue, 
		       std::vector< std::vector<size_t> > 
		       &indValue );
    
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
}
#endif
