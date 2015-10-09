
// Filename     : MembraneCyclingAll.h
// Description  : Classes describing cycling to/from mebrane
// Author(s)    : Laura Brown (laura.brown@slcu.cam.ac.uk))
// Created      : August 2014
// Revision     : $Id:$

// This file is a modification of MembraneCycling.h file by Pau
#ifndef MEMBRANECYCLINGAll_H
#define MEMBRANECYCLINGAll_H

#include<cmath>

#include"tissue.h"
#include"baseReaction.h"



///
/// @brief Membrane Cycling describes reactions that give the cycling of a protein to and from the membrane/Wall.
///
namespace MembraneCyclingAll {




///
/// @brief A function describing the constant exocytosis and endocytosis of PIN (or another protein) from the cytosol to the cell membrane at a constant rate. 
///
/// It uses two compartments for each wall and a single for the cells. p0 gives exocytosis rate, p1 endocytosis rate.
/// PIN  molecules are updated according to:
///
/// @f[ \frac{dP_i}{dt} = -p_0 P_i +\sum_{j} p_1 P_{ij} @f] 
///  
/// @f[ \frac{dP_{ij}}{dt} = p_0 P_i -p_1 P_{ij} @f]
///  
/// In the model file the reaction is given by:
/// @verbatim
/// MembraneCyclingAll::Constant 2 2 1 1
/// p_0 p_1
/// ci_PIN 
/// wi_PIN 
/// @endverbatim
///
class Constant : public BaseReaction {
  
 public:
  
 Constant(std::vector<double> &paraValue, 
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
/// @brief A function describing the exocytosis and endocytosis of PIN (or another protein) from the cell membrane to the cytosol at a  rate
/// dependent on the amount of auxin in the wall compartment.
///
/// It uses two compartments for each wall and a single for the cells. p0 gives exocytosis rate, p1 endocytosis rate.
/// PIN  molecules are updated according to:
/// @f[ \frac{dP_i}{dt} = \sum_{j}- p_0 P_i \frac{X_{ij}^{p_3}}{X_{ij}^{p_3}+{p_2}^{p_3}}+ p_1 P_{ij} \frac{X_{ij}^{p_3}}{X_{ij}^{p_3}+{p_2}^{p_3}} @f] 
///  
/// @f[ \frac{dP_{ij}}{dt} =  p_0 P_i \frac{X_{ij}^{p_3}}{X_{ij}^{p_3}+{p_2}^{p_3}}- p_1 P_{ij} \frac{X_{ij}^{p_3}}{X_{ij}^{p_3}+{p_2}^{p_3}} @f]  @f]
///
/// In the model file the reaction is given by:
/// @verbatim
/// MembraneCyclingAll::LocalWallFeedbackNonLinear 4 2 1 2
/// p_0 ..p_2
/// ci_PIN 
/// wi_X  Wi_PIN 
/// @endverbatim
///
class LocalWallFeedbackNonLinear : public BaseReaction {
  
 public:
  
 LocalWallFeedbackNonLinear(std::vector<double> &paraValue, 
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



