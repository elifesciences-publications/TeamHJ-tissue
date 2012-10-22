//
// Filename     : centerTriangulation.h
// Description  : Classes describing updates for tissues with cells storing a central point (and internal edges)
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : October 2012
// Revision     : $Id:$
//
#ifndef CENTERTRIANGULATION_H
#define CENTERTRIANGULATION_H

#include"tissue.h"
#include"baseReaction.h"
#include<cmath>

namespace CenterTriangulation
{

  ///
  /// @brief Initiates a center triangulation, either from cells in tissue or from scratch
  ///
  /// This reaction does not update the tissue. It only initiates a central triangulation of cells
  /// and add variables to cellData at initiation. This is mainly used for triangular biquadratic
  /// spring models. One parameter (flag) can be provided and is setr to one if the central point and
  /// internal edges should be initiated from scratch even if they have been provided to the tissue
  /// when reading the init. One variable index is provided for compability with an old version, and
  /// should represents the end of the cellData vector (cell(i).numVariable()). 
  ///
  /// In a model file the reaction is given by
  ///
  /// @verbatim
  /// CeterTriangulation::Initiate 0/1 1 1
  /// [overrideFlag]
  /// InternalVarStartIndex
  /// @endverbatim
  ///
  /// @note The order of reactions may be of importance when using this reaction, e.g. if vertex are added on the walls.
  /// @see Tissue::readInitCenterTri()
  /// @see Cell::centerPosition()
  ///
  class Initiate : public BaseReaction {
  
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
    Initiate(std::vector<double> &paraValue, 
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
} // end namespace CenterTriangulation
  
#endif //CENTERTRIANGULATION_H
