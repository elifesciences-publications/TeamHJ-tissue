
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
/// @brief A function discribing the constant exocytosis and endocytosis of PIN (or another protein) from the cytosol to the cell membrane at a constant rate. p0 gives exocytosis, p1 endocytosis.
///
/// It uses two compartments for each wall and a single for the cells.
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
/// @brief A function discribing the endocytosis of PIN (or another protein) from the cell membrane to the cytosol at a  rate
/// Dependent on the amount of PIN on the adjacent membrane of a neighbouring cell.
///
/// It uses two compartments for each wall and a single for the cells.
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
/// @brief A function discribing the endocytosis of PIN (or another protein) from the cell membrane to the cytosol at a  rate
/// Dependent on the amount of PIN on the adjacent membrane of a neighbouring cell.
///
/// It uses two compartments for each wall and a single for the cells.
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
/// @brief A function discribing the exocytosis of PIN (or another protein) from the cell membrane to the cytosol at a  rate
/// Dependent on the amount of auxin in the neighbouring cell.
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
/// CrossMembraneExocytosis 3 2 1 1
/// p_0 ..p_2
/// ci_PIN 
///  wi_PIN 
/// @endverbatim
///
class CellUpTheGradientExocytosis : public BaseReaction {
  
 public:
  
 CellUpTheGradientExocytosis(std::vector<double> &paraValue, 
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


}


#endif


