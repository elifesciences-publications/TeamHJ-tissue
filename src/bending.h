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
/// @brief Bending describes reactions that generates mechanical updates of bending moments.
///
namespace Bending {
  
  /// 
  /// @brief Creates a bending resistance by a force towards the line connecting the neighboring vertices 
  ///
  /// @details This reaction applies a force acting on a vertex towards the weighted (by edge length) middle
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
  /// @note It might be good to add a maximal angle to allow for a bent structure.
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

  /// 
  /// @brief Creates a bending resistance by forces to keep an angle at a specific value
  ///
  /// @details Each vertex has a prefered angle, and a bendng moment is introduced to
  /// move vertices to keep this angle. The update for the angle at a specific vertex comes from the potential
  /// 
  /// @f[ V(\Theta,\Theta^{t}) = \frac{1}{2}k_{\Theta} (\Theta-\Theta^{t})^2 @f]
  ///
  /// where vertex positions are updated according to
  ///
  /// @f[ \frac{dx_i}{dt} = -\frac{dV}{dx} = \frac{dV}{d\Theta}\frac{d\Theta}{dF}\frac{dF}{dx_i}  @f]
  ///
  /// where the angle is calculated by
  ///
  /// @f[ \Theta = acos(F) \pm \Pi @f]
  ///
  /// with
  ///
  /// @f[ F = \frac{x_1 \dot x_2}{|x_1||x_2|} = \frac{f}{g} @f]
  ///
  /// where x1, x2, are the two edges joined at the vertex, -/+ comes from if the structure is convex/concave. 
  /// For each angle three edges are involved and contribute (in all dimensions) to the update (assuming x-,x,x+
  /// are positions for the vertices:
  ///
  /// @f[ \frac{dx_{-}}{dt} = k_{\Theta} (\Theta-\Theta^{t}) \frac{1}{(1 - F^2)^{0.5}} 
  /// \frac{1}{g} [(x-x_+)+(x-x_-) \frac{f x_2}{g x_1}] @f]
  /// @f[ \frac{dx}{dt} = k_{\Theta} (\Theta-\Theta^{t}) \frac{1}{(1 - F^2)^{0.5}}
  /// \frac{1}{g} [(x_+ - x)(1+\frac{f x_1}{g x_2}) + (x_- - x) (1+\frac{f x_2}{g x_1})] @f]
  /// @f[ \frac{dx_{+}}{dt} = k_{\Theta} (\Theta-\Theta^{t}) \frac{1}{(1 - F^2)^{0.5}}
  /// \frac{1}{g} [(x-x_-)+(x-x_+) \frac{f x_1}{g x_2}] @f]
  ///
  /// The reaction expects a strength of the bending force, and the edge variable index for the angle. 
  /// The prefered angle is stored as a second wall variable. In a model file it is defined as:
  /// @verbatim
  /// Bending::Angle 1 1 1
  /// k_theta
  /// Theta^t_index
  /// @endverbatim
  ///
  /// @see Bending::AngleRelax for how to update the prefered angle towards the current.
  /// @see Bending::AngleInitiate will initiate the stored angle to the current measured angle.
  ///
  class Angle : public BaseReaction {
    
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
    Angle(std::vector<double> &paraValue, 
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

  /// 
  /// @brief Initiates an angle variable to the current value of the angle at vertices
  ///
  /// @details Sets the prefered angle for a vertex to the initial calculated value.
  /// The angle is calculated by
  ///
  /// @f[ \Theta = acos(F) \pm \Pi @f]
  ///
  /// with
  ///
  /// @f[ F = \frac{x_1 \dot x_2}{|x_1||x_2|} = \frac{f}{g} @f]
  ///
  /// where x1, x2, are the two edges joined at the vertex, -/+ comes from if the structure is convex/concave. 
  ///
  /// In a model file the reaction is defined by:
  /// @verbatim
  /// Bending::AngleInitiate 0 1 1
  /// Theta^t_index
  /// @endverbatim
  /// where the index is where the angle is stored (as wall variable).
  ///
  /// @see Bending::Angle for how to update the vertices towards this angle 
  /// @see Bending::AngleRelax for how to update the prefered angle towards the current.
  ///
  class AngleInitiate : public BaseReaction {
    
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
    AngleInitiate(std::vector<double> &paraValue, 
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
  /// @brief Initiates an angle variable to the current value of the angle at vertices
  ///
  /// @details Sets the prefered angle for a vertex to the initial calculated value.
  /// The angle is calculated by
  ///
  /// @f[ \Theta = acos(F) \pm \Pi @f]
  ///
  /// with
  ///
  /// @f[ F = \frac{x_1 \dot x_2}{|x_1||x_2|} = \frac{f}{g} @f]
  ///
  /// where x1, x2, are the two edges joined at the vertex, -/+ comes from if the structure is convex/concave. 
  ///
  /// In a model file the reaction is defined by:
  /// @verbatim
  /// Bending::AngleRelax 0 1 1
  /// Theta^t_index
  /// @endverbatim
  /// where the index is where the angle is stored (in wall variable).
  ///
  /// @see Bending::Angle for how to update the vertices towards this angle 
  /// @see Bending::AngleInitiate for how to initiate the angle variable to the current at start.
  ///
  class AngleRelax : public BaseReaction {
    
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
    AngleRelax(std::vector<double> &paraValue, 
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
