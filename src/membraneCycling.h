
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
/// @brief Membrane Cycling describes reactions that give the cycling of a protein to and from the membrane/Wall.
///
namespace MembraneCycling {




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
/// @brief A function describing the exocytosis and endocytosis of PIN (or another protein) from the cytosol to the cell
/// membrane at a rate dependent on the amount of PIN on the adjacent membrane of a neighbouring cell, as a non-linear 
/// function.
///
/// It uses two compartments for each wall and a single for the cells. p0 gives maximal exocytosis rate, p1 the maximal endocytosis rate.
/// PIN  molecules are updated according to:
///  
///
/// @f[ \frac{dP_i}{dt} =\sum_{j} -p_0 P_i \frac{P_{ji}^{p_3}}{P_{ji}^{p_3}+p_2^{p_3}}+\sum_{j} p_1 P_{ij} \frac{P_{ji}^{p_3}}{P_{ji}^{p_3}+{p_2}^{p_3}} @f] 
///  
/// @f[ \frac{dP_{ij}}{dt} = p_0 P_i \frac{P_{ji}^{p_3}}{P_{ji}^{p_3}+p_2^{p_3}}- p_1 P_{ij} \frac{P_{ji}^{p_3}}{P_{ji}^{p_3}+{p_2}^{p_3}} @f]
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



/// @brief A function describing the exocytosis and endocytosis of PIN (or another protein) from the cytosol to the cell
/// membrane at a rate dependent on the amount of PIN on the adjacent membrane of a neighbouring cell, as a linear 
/// function.
///
/// It uses two compartments for each wall and a single for the cells. p0 gives maximal exocytosis rate, p1 the maximal endocytosis rate.
/// PIN  molecules are updated according to:
///
/// @f[ \frac{dP_i}{dt} =\sum_{j} -p_0 P_i P_{ji}+ p_1 P_ij P_{ji} @f] 
///  
/// @f[ \frac{dP_{ij}}{dt} = p_0 P_i P_{ji}- p_1 P_i P_{ji} @f]
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
/// MembraneCycling::LocalWallFeedbackNonLinear 4 2 1 2
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


///
/// @brief A function describing the exocytosis and endocytosis of PIN (or another protein) from the cell membrane to the cytosol at a  rate
/// dependent on the amount of auxin in the wall compartment.
///
/// It uses two compartments for each wall and a single for the cells. p0 gives exocytosis, p1 endocytosis.
/// PIN  molecules are updated according to:
///
/// @f[ \frac{dP_i}{dt} =\sum_{j}- p_0 P_i X_{ij}+ p_1 P_{ij} X_{ij} @f] 
///  
/// @f[ \frac{dP_{ij}}{dt} = p_0 P_i X_{ij}- p_1 P_{ij} X_{ij} @f]
///
/// In the model file the reaction is given by:
/// @verbatim
/// MembraneCycling::LocalWallFeedbackLinear 2 2 1 2
/// p_0 p_1
/// ci_PIN 
/// wi_X  wi_PIN
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
/// @brief A function describing the exocytosis and endocytosis of PIN (or another protein) from the cell membrane to the cytosol at a  rate
/// dependent on the amount of auxin in neighbouring cell.
///
/// It uses two compartments for each wall and a single for the cells.
/// PIN  molecules are updated according to:
///  
///
/// @f[ \frac{dP_i}{dt} =\sum_{j} - p_0 P_{i} \frac{A_{ji}^{p_3}}{A_{ji}^{p_3}+{p_2}^{p_3}}+ p_1 P_{ij} \frac{A_{ji}^{p_3}}{A_{ji}^{p_3}+{p_2}^{p_3}} @f] 
///  
/// @f[ \frac{dP_{ij}}{dt} = p_0 P_{i} \frac{A_{ji}^{p_3}}{A_{ji}^{p_3}+{p_2}^{p_3}}- p_1 P_{ij} \frac{A_{ji}^{p_3}}{A_{ji}^{p_3}+{p_2}^{p_3}} @f] 
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
/// @f[ \frac{dP_i}{dt} =\sum_{j} -p_0 P_{i} A_{ji}+p_1 P_{ij} A_{ji} @f] 
///  
/// @f[ \frac{dP_{ij}}{dt} =  p_0 P_{i} A_{ji}-p_1 P_{ij} A_{ji} @f]
///
/// In the model file the reaction is given by:
/// @verbatim
/// membraneCycling::CellUpTheGradientLinear 2 2 2 1
/// p_0 p_1
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
/// @brief A function describing the exocytosis of PIN (or another protein) from the cell membrane to the cytosol at a  rate
/// dependent on the amount of auxin in the cell.
///
/// It uses two compartments for each wall and a single for the cells.
/// PIN  molecules are updated according to:
///
/// @f[ \frac{dP_i}{dt} =- p_0 P_i \frac{A_{ij}^{p_3}}{A_{ij}^{p_3}+{p_2}^{p_3}}+ sum_{j} p_1 P_{ij} \frac{A_{ij}^{p_3}}{A_{ij}^{p_3}+{p_2}^{p_3}} @f] 
///  
/// @f[ \frac{dP_{ij}}{dt} = p_0 P_i \frac{A_{ij}^{p_3}}{A_{ij}^{p_3}+{p_2}^{p_3}}- p_1 P_{ij} \frac{A_{ij}^{p_3}}{A_{ij}^{p_3}+{p_2}^{p_3}} @f]
///
/// In the model file the reaction is given by:
/// @verbatim
/// membraneCycling::InternalCellNonLinear 4 2 2 1
/// p_0 ..p_2
/// ci_Auxin ci_PIN 
/// wi_PIN 
/// @endverbatim
///
class InternalCellNonLinear : public BaseReaction {
  
 public:
  
 InternalCellNonLinear(std::vector<double> &paraValue, 
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
/// dependent on the amount of auxin in the  cell.
///
/// It uses two compartments for each wall and a single for the cells.
/// PIN  molecules are updated according to:
///
/// @f[ \frac{dP_i}{dt} = - p_0 P_i A_{ij}+\sum_{j}p_1 P_{ij} A_{ij} @f] 
///  
/// @f[ \frac{dP_{ij}}{dt} =  p_0 P_i A_{ij}-p_1 P_{ij} A_{ij} @f]
///
/// In the model file the reaction is given by:
/// @verbatim
/// membraneCycling::InternalCellLinear 2 2 2 1
/// p_0 p_1
/// ci_Auxin ci_PIN 
/// wi_PIN 
/// @endverbatim
///
class InternalCellLinear : public BaseReaction {
  
 public:
  
 InternalCellLinear(std::vector<double> &paraValue, 
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
/// dependent on the the flux of auxin.
///
/// It uses two compartments for each wall and a single for the cells.
/// PIN  molecules are updated according to:
///  To do
///
/// In the model file the reaction is given by:
/// @verbatim
/// CellFluxExocytosis 3 2 1 1
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
/// @brief A function describing the exocytosis and endocytosis of PIN (or another protein) from the cell membrane to the cytosol at a  rate
/// dependent on the amount of PIN in the wall compartment.
///
/// It uses two compartments for each wall and a single for the cells. p0 gives exocytosis, p1 endocytosis.
/// PIN  molecules are updated according to:
///  

/// @f[ \frac{dP_i}{dt} =\sum_{j} -p_0 P_{ij} \frac{P_{ij}^{p_3}}{P_{ij}^{p_3}+{p_2}^{p_3}} +p_1 P_{ij} \frac{P_{ij}^{p_3}}{P_{ij}^{p_3}+{p_2}^{p_3}} @f] 
///  
/// @f[ \frac{dP_{ij}}{dt} =  p_0 P_{ij} \frac{P_{ij}^{p_3}}{P_{ij}^{p_3}+{p_2}^{p_3}}-p_1 P_{ij} \frac{P_{ij}^{p_3}}{P_{ij}^{p_3}+{p_2}^{p_3}} @f] 
///
///
///
///  
/// In the model file the reaction is given by:
/// @verbatim
/// MembraneCycling:: PINFeedbackNonLinear 4 2 1 1
/// p_0 ..p_3
/// ci_PIN 
/// wi_PIN  
/// @endverbatim
///
class PINFeedbackNonLinear : public BaseReaction {
  
 public:
  
 PINFeedbackNonLinear(std::vector<double> &paraValue, 
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
/// dependent on the amount of PIN in the wall compartment.
///
/// It uses two compartments for each wall and a single for the cells. p0 gives exocytosis, p1 endocytosis.
/// PIN  molecules are updated according to:
///  
///
/// @f[ \frac{dP_i}{dt} =\sum_{j}- p_0 P_i X_{ij}+ p_1 P_i X_{ij} @f] 
///  
/// @f[ \frac{dP_{ij}}{dt} = p_0 P_i X_{ij}-p_1 P_i X_{ij} @f]
///  
/// In the model file the reaction is given by:
/// @verbatim
/// MembraneCycling::PINFeedbackLinear 2 2 1 1
/// p_0 p_1
/// ci_PIN 
/// wi_PIN
/// @endverbatim
///
class PINFeedbackLinear : public BaseReaction {
  
 public:
  
 PINFeedbackLinear(std::vector<double> &paraValue, 
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



