//
// Filename     : bending.h
// Description  : Classes describing reactions related to bending moments
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : September 2013
// Revision     : $Id:$
//
#ifndef BENDING_H
#define BENDING_H

#include"tissue.h"
#include"baseReaction.h"

///
/// @brief Bending describes reactions that generates mechanical updates of bending moments
///
namespace Bending {
  
  /// 
  /// @brief Creates a bending resistance by a force towards the line connecting the neighboring vertices 
  ///
  /// This reaction applies a force acting on a vertex towards the weighted (by edge length) middle
  /// of the line connecting the left and right vertex neighbors within the cell. The force is a spring force
  /// acting as soon as the vertex is outside the line (i.e. the edges connected to the vertex have an angle
  /// in between). The update is given by
  ///
  /// @f[ \frac{dx_i}{dt} = - k_{spring} (x_{i}-\frac{1}{L_{-}+L_{+}}( L_{-}x_{+} + L_{+}x_{-})) @f]
  ///
  /// where
  ///
  /// @f[ \frac{1}{L_{-}+L_{+}}( L_{-}x_{+} + L_{+}x_{-}) @f]
  ///
  /// is the weighted middle on the line between left(-) and right(+) vertex neighbors, @$L_{-},L_{+}@$ 
  /// are the edge lengths towards the left and right vertex neighbors and @$x_{-},x_{+}@$ ar the positions.
  /// The update is done in all dimensions.
  ///
  /// The reaction expects a strength of the spring force, and the edge variable index for the length (=0). 
  /// In a model file it is defined as:
  /// @verbatim
  /// Bending::NeighborCenter 1 1 1
  /// k_spring
  /// L_index
  /// @endverbatim
  ///
  /// @Note It might be good to add a maximal angle to allow for a bent structure.
  ///
  class NeighborCenter : public BaseReaction {
    
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
    NeighborCenter(std::vector<double> &paraValue, 
		   std::vector< std::vector<size_t> > 
		   &indValue );
    
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
}
#endif
