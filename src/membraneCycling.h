// Filename     : MembraneCycling.h
// Description  : Classes describing cycling to/from mebrane
// Author(s)    : Laura Brown (laura.brown@slcu.cam.ac.uk))
// Created      : August 2014
// Revision     : $Id:$
//
#ifndef MEMBRANECYCLING_H
#define MEMBRANECYCLING_H

#include<cmath>

#include"tissue.h"
#include"baseReaction.h"



///
/// @brief Membrane Cycling describes reactions that give the cycling of a protein to anf from the membrane/Wall.
///
namespace MembraneCycling {




///
/// @brief A function describing the constant exocytosis and endocytosis of PIN (or another protein) from the cytosol to the cell membrane at a constant rate. 
///
/// It uses two compartments for each wall and a single for the cells. p0 gives exocytosis, p1 endocytosis.
/// PIN  molecules are updated according to:
///  
///
/// @f[ \frac{dP_i}{dt} = -p_0 P_i +p_1 P_ij @f] 
///  
/// @f[ \frac{dP_{ij}}{dt} = p_0 P_i -p_0 P_ij @f]
///
///
///  
/// In the model file the reaction is given by:
/// @verbatim
/// MembraneCycling::Constant 2 2 1 1
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
/// @brief A function describing the endocytosis of PIN (or another protein) from the cell membrane to the cytosol at a  rate
/// Dependent on the amount of PIN on the adjacent membrane of a neighbouring cell, as a non-linear function.
///
/// It uses two compartments for each wall and a single for the cells. p0 gives exocytosis, p1 endocytosis.
/// PIN  molecules are updated according to:
///  
///
/// @f[ \frac{dP_i}{dt} = -p_0 P_i \frac{Pji^{p_3}}{Pji^{p_3}+p_2^{p_3}}+ p_1 P_ij \frac{Pji^{p_3}}{Pji^{p_3}+{p_2}^{p_3}} @f] 
///  
/// @f[ \frac{dP_{ij}}{dt} = p_0 P_i \frac{Pji^{p_3}}{Pji^{p_3}+p_2^{p_3}}- p_1 P_ij \frac{Pji^{p_3}}{Pji^{p_3}+{p_2}^{p_3}} @f]
///
///
///  
/// In the model file the reaction is given by:
/// @verbatim
/// MembraneCycling::CrossMembraneNonLinear 4 2 1 1
/// p_0 ..p_3
/// ci_PIN 
/// wi_PIN 
/// @endverbatim
///
class CrossMembraneNonLinear : public BaseReaction {
  
 public:
  
 CrossMembraneNonLinear(std::vector<double> &paraValue, 
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
/// @brief A function describing the endocytosis of PIN (or another protein) from the cell membrane to the cytosol at a  rate
/// Dependent on the amount of PIN on the adjacent membrane of a neighbouring cell, as a linear function.
///
/// It uses two compartments for each wall and a single for the cells. p0 gives exocytosis, p1 endocytosis.
/// PIN  molecules are updated according to:
///  
///
/// @f[ \frac{dP_i}{dt} =- p_0 P_i Pji+ p_1 P_ij Pji @f] 
///  
/// @f[ \frac{dP_{ij}}{dt} = p_0 P_i Pji- p_1 P_i Pji @f]
///
///
///  
/// In the model file the reaction is given by:
/// @verbatim
/// MembraneCycling::CrossMembraneLinear 2 2 1 1
/// p_0 p_1
/// ci_PIN 
///  wi_PIN 
/// @endverbatim
///
class CrossMembraneLinear : public BaseReaction {
  
 public:
  
 CrossMembraneLinear(std::vector<double> &paraValue, 
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
/// @brief A function describing the exocytosis of PIN (or another protein) from the cell membrane to the cytosol at a  rate
/// Dependent on the amount of PIN/auxin in the wall compartment.
///
/// It uses two compartments for each wall and a single for the cells. p0 gives exocytosis, p1 endocytosis.
/// PIN  molecules are updated according to:
///  
///
/// @f[ \frac{dP_i}{dt} = p_0 P_ij \frac{Xij^{p_2}}{Xij^{p_2}+{p_1}^{p_2}} @f] 
///  
/// @f[ \frac{dP_{ij}}{dt} = -p_0 P_ij \frac{Xij^{p_2}}{Xij^{p_2}+{p_1}^{p_2}} @f]
///
///
///  
/// In the model file the reaction is given by:
/// @verbatim
/// MembraneCycling::LocalWallFeedbackNonLinear 4 2 1 2
/// p_0 ..p_2
/// ci_PIN 
/// wi_PIN  Wi_X 
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


///
/// @brief A function describing the exocytosis of PIN (or another protein) from the cell membrane to the cytosol at a  rate
/// Dependent on the amount of PIN/auxin in the wall compartment.
///
/// It uses two compartments for each wall and a single for the cells. p0 gives exocytosis, p1 endocytosis.
/// PIN molecules are updated according to:
///  
///
/// @f[ \frac{dP_i}{dt} = p_0 P_ij Xij @f] 
///  
/// @f[ \frac{dP_{ij}}{dt} = -p_0 P_ij Xij @f]
///
///
///  
/// In the model file the reaction is given by:
/// @verbatim
/// MembraneCycling::LocalWallFeedbackLinear 2 2 1 2
/// p_0 p_1
/// ci_PIN 
///  wi_PIN wi_X 
/// @endverbatim
///
class LocalWallFeedbackLinear : public BaseReaction {
  
 public:
  
 LocalWallFeedbackLinear(std::vector<double> &paraValue, 
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
/// @brief A function describing the exocytosis of PIN (or another protein) from the cell membrane to the cytosol at a  rate
/// Dependent on the amount of auxin in neighbouring cell.
///
/// It uses two compartments for each wall and a single for the cells.
/// PIN  molecules are updated according to:
///  
///
/// @f[ \frac{dP_i}{dt} = p_0 P_ij \frac{Pji^{p_2}}{Pji^{p_2}+{p_1}^{p_2}} @f] 
///  
/// @f[ \frac{dP_{ij}}{dt} = -p_0 P_ij \frac{Pji^{p_2}}{Pji^{p_2}+{p_1}^{p_2}} @f]
///
///
///  
/// In the model file the reaction is given by:
/// @verbatim
/// membranCycling::CellUpTheGradientNonLinear 4 2 2 1
/// p_0 ..p_2
/// ci_Auxin ci_PIN 
/// wi_PIN 
/// @endverbatim
///
class CellUpTheGradientNonLinear : public BaseReaction {
  
 public:
  
 CellUpTheGradientNonLinear(std::vector<double> &paraValue, 
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
/// @brief A function describing the exocytosis of PIN (or another protein) from the cell membrane to the cytosol at a  rate
/// Dependent on the amount of auxin in neighbouring cell.
///
/// It uses two compartments for each wall and a single for the cells.
/// PIN  molecules are updated according to:
///  
///
/// @f[ \frac{dP_i}{dt} = p_0 P_ij Aji @f] 
///  
/// @f[ \frac{dP_{ij}}{dt} = -p_0 P_ij Aji @f]
///
///
///  
/// In the model file the reaction is given by:
/// @verbatim
/// membraneCycling::CellUpTheGradientNonLinear 4 2 2 1
/// p_0 ..p_2
/// ci_Auxin ci_PIN 
/// wi_PIN 
/// @endverbatim
///
class CellUpTheGradientLinear : public BaseReaction {
  
 public:
  
 CellUpTheGradientLinear(std::vector<double> &paraValue, 
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
/// @brief A function discribing the exocytosis of PIN (or another protein) from the cell membrane to the cytosol at a  rate
/// Dependent on the amount of PIN on the adjacent membrane of a neighbouring cell.
///
/// It uses two compartments for each wall and a single for the cells.
/// PIN  molecules are updated according to:
///  
///
/// @f[ \frac{dP_i}{dt} = p_0 P_ij \frac{Pji^{p_2}}{Pji^{p_2}+{p_1}^{p_2}} @f] 
///  
/// @f[ \frac{dP_{ij}}{dt} = -p_0 P_ij \frac{Pji^{p_2}}{Pji^{p_2}+{p_1}^{p_2}} @f]
///
///
///  
/// In the model file the reaction is given by:
/// @verbatim
/// CrossMembraneEndocytosis 3 2 1 1
/// p_0 ..p_2
/// ci_PIN 
///  wi_PIN 
/// @endverbatim
///
class CellFluxExocytosis : public BaseReaction {
  
 public:
  
 CellFluxExocytosis(std::vector<double> &paraValue, 
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
/// @brief A function describing the exocytosis of PIN (or another protein) from the cytosol to the membrane at a rate
/// dependent on the concentration of another molecule in the cytosol.
///
/// It uses two compartments for each wall and a single for the cells. p0 gives exocytosis, p1 endocytosis.
/// PIN molecules are updated according to:
///  
///
/// @f[ \frac{dP_i}{dt} = p_0 P_ij - p_1 P_i X_i@f] 
///  
/// @f[ \frac{dP_{ij}}{dt} = -p_0 P_ij + p_1 P_i X_i @f]
///
///
///  
/// In the model file the reaction is given by:
/// @verbatim
/// MembraneCycling::LocalWallFeedbackLinear 2 2 2 1
/// p_0 p_1
/// ci_PIN ci_X
///  wi_PIN
/// @endverbatim
///
class LocalCellWallFeedbackLinear : public BaseReaction {
  
 public:
  
 LocalCellWallFeedbackLinear(std::vector<double> &paraValue, 
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


