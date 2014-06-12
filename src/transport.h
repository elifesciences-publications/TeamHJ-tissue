//
// Filename     : transport.h
// Description  : Classes describing transport reactions
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : September 2013
// Revision     : $Id:$
//
#ifndef TRANSPORT_H
#define TRANSPORT_H

#include<cmath>

#include"tissue.h"
#include"baseReaction.h"

///
/// @brief A membrane diffusion reaction
///
/// A reaction for passive diffusion of molecules localized in the membrane. The
/// transport is between neighboring membrane compartments within the same cell
/// described by:
///  
/// @f[ \frac{dP_{ij}}{dt} = - p_0 ( 2 P_{ij} - P_{ij} - P_{ij}) @f] 
///  
/// where p_0 is the diffusion rate, i is the cell, j a membrane section, anf j+/- neighboring membrane sections.
///  
/// In a model file the reaction is defined as
///
/// @verbatim
/// MembraneDiffusionSimple 1 1 1
/// p_0
/// P_{wallindex}
/// @endverbatim
///
/// where the reaction assumes that each wall keeps two variables per membrane molecule.
///
/// @note The Simple in the name reflects the fact that no geometric factors are included.
///
class MembraneDiffusionSimple : public BaseReaction {
  
 public:
  
  MembraneDiffusionSimple(std::vector<double> &paraValue, 
			  std::vector< std::vector<size_t> > 
			  &indValue );
  
  void derivs(Tissue &T,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      DataMatrix &cellDerivs,
	      DataMatrix &wallDerivs,
	      DataMatrix &vertexDerivs );
};


///
/// @brief A membrane diffusion reaction
///
/// A reaction for passive diffusion of molecules localized in the membrane. The
/// transport is between neighboring membrane compartments within the same cell
/// described by:
///  
/// @f[ \frac{dP_{ij}}{dt} = - p_0 ( 2 P_{ij} - P_{ij} - P_{ij}) @f] 
///  
/// where p_0 is the diffusion rate, i is the cell, j a membrane section, anf j+/- neighboring membrane sections.
///  
/// In a model file the reaction is defined as
///
/// @verbatim
/// MembraneDiffusionSimple 1 1 1
/// p_0
/// P_{wallindex}
/// @endverbatim
///
/// where the reaction assumes that each wall keeps two variables per membrane molecule.
///
/// @note The Simple in the name reflects the fact that no geometric factors are included.
///

class MembraneDiffusionSimple2 : public BaseReaction {
  
 public:
  
  MembraneDiffusionSimple2(std::vector<double> &paraValue, 
			  std::vector< std::vector<size_t> > 
			  &indValue );
  
  void derivs(Tissue &T,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      DataMatrix &cellDerivs,
	      DataMatrix &wallDerivs,
	      DataMatrix &vertexDerivs );
};






///
/// @brief A cell to cell transport reaction
///
/// A reaction for transport molecules from cell to cell dependent on a membrane localised efflux carrier. The
/// transport is between neighboring cells
/// is described by:
///  
///  @f[ \frac{dA_i}{dt} = @f[ p_0 \sum_{neigh} (P_{ni} A_n - P_{in} A_i) @f] 
///
///  
/// where p_0 is the tranpsort rate, i is the cell and n neighboring membrane sections.
///  
/// In a model file the reaction is defined as
///
/// @verbatim
/// ActiveTransportCellEfflux 1 2 1 1
/// p_0
/// A_{cellIndex}
/// P_{wallindex}
/// @endverbatim
///
/// where the reaction assumes that each wall keeps two variables per membrane molecule.
///
///
///
 

class ActiveTransportCellEfflux  : public BaseReaction {
  
 public:
  
  ActiveTransportCellEfflux(std::vector<double> &paraValue, 
			  std::vector< std::vector<size_t> > 
			  &indValue );
  
  void derivs(Tissue &T,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      DataMatrix &cellDerivs,
	      DataMatrix &wallDerivs,
	      DataMatrix &vertexDerivs );
};





#endif


