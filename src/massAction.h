// Filename     : MassAction.h
// Description  : Classes describing mass action reactions
// Author(s)    : Laura Brown (laura.brown@slcu.cam.ac.uk))
// Created      : August 2014
// Revision     : $Id:$
//
#ifndef MASSACTION_H
#define  MASSACTION_H

#include<cmath>

#include"tissue.h"
#include"baseReaction.h"



///
/// @brief A collect of mass action reactions
///
namespace MassAction {





///
/// @brief A one way mass action reaction assume one reactants and two products.
///
///
///  
///
/// @f[ \frac{dP}{dt} = k_1 R_1 R_2 @f] 
///  
/// @f[ \frac{dR_i}{dt} = - k_1 R_1 R_2 @f]
///
///
///  
/// In the model file the reaction is given by:
/// @verbatim
/// MassAction::OneToTwo 1 2 3 0
/// p_0
/// r_cell p1_cell P2_cell 
/// @endverbatim
///
class OneToTwo : public BaseReaction {
  
 public:
  
 OneToTwo(std::vector<double> &paraValue, 
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
/// @brief A one way mass action reaction assume two reactants and a single product
///
///
///  
///
/// @f[ \frac{dP}{dt} = k_1 R_1 R_2 @f] 
///  
/// @f[ \frac{dR_i}{dt} = - k_1 R_1 R_2 @f]
///
///
///  
/// In the model file the reaction is given by:
/// @verbatim
/// MassAction::TwoToOne 1 2 3 0
/// p_0
/// r1_cell r2_cell P_cell 
/// @endverbatim
///
class TwoToOne : public BaseReaction {
  
 public:
  
 TwoToOne(std::vector<double> &paraValue, 
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





}


#endif

