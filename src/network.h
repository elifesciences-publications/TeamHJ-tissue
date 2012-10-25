//
// Filename     : network.h
// Description  : Classes describing complete models updates
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : November 2006
// Revision     : $Id:$
//
#ifndef NETWORK_H
#define NETWORK_H

#include<cmath>

#include"tissue.h"
#include"baseReaction.h"

/// @brief A stress-based PIN1 and MT polarization model
///
///
class AuxinModelStress : public BaseReaction {
  
 public:
  
  AuxinModelStress(std::vector<double> &paraValue, 
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
/// @brief A linear polarization cell-based auxin transport model
///
/// A complete pattern generating auxin model based on only cellular
/// compartments. The four molecules A(uxin), P(IN), X(auxin induced
/// molecule), and M(epidermally expressed molecule) are updated
/// according to:
///  
/// @f[ \frac{dA_i}{dt} = p_0 M_i + p_1 - p_2 A_i + p_4 \sum_{n}^{neigh} (P_{ni} A_n - P_{in} A_i) + p_5 \sum_{n}^{neigh} (A_n-A_i) @f] 
///  
/// @f[ \frac{dP_i}{dt} = p_6 - p_7 P_i @f] 
///  
/// @f[ \frac{dX_i}{dt} = p_8 A_i - p_9 X_i @f]
///  
/// @f[ \frac{dM_i}{dt} = p_{10} \Theta_{L1} - p_{11} M_i @f]
///  
/// @f[ P_{in} = \frac{P_i X_n}{(p_3 + \sum_{k}^{neigh} X_k)} @f]
///  
/// In a model file the reaction is defined as
///
/// @verbatim
/// AuxinModelSimple1 12 1 4
/// p_0 ... p_11
/// A_index P_index X_index M_index
/// @endverbatim
///
class AuxinModelSimple1 : public BaseReaction {
  
 public:
  
  AuxinModelSimple1(std::vector<double> &paraValue, 
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

//!A wall-based auxin transport model
/*!A complete auxin transport model based on PIN polarization from a linear
  feedback from a wall variable. PIN and auxin are updated
  according to:
  
  dA_i/dt = p0 - p1*A_i +p4*\Sum_{neigh} (A_n-A_i) + 
  p3*\Sum_{neigh} (P_ni*A_n-P_in*A_i) 
  
  dP_i/dt = p5 - p6*P_i 
  
  P_in = P_i*X_in/(p_2+\Sum_{k,neigh}X_ik)
  
  where X_in is the variable in the wall.
  In addition, the column index for auxin, PIN, X (in wall) should be given.
*/
class AuxinModelSimple1Wall : public BaseReaction {
  
 public:
  
  AuxinModelSimple1Wall(std::vector<double> &paraValue, 
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
/// @brief A cell-based auxin transport model where the PIN polarization comes from wall stresses
///
/// An auxin transport model for cellular auxin and where PIN polarization is based on stresses
/// in a double wall compartment, i.e. the wall compartment between a cell pair is divided into
/// two compartments, and stress is measured in the wall compartment neighbor to the specific cell.
/// Auxin concentration is affecting the spring constants in wall pieces neighboring the cell.
///
/// PIN and auxin are updated according to:
///
///	dA_i/dt = p0 - p1*A_i +p2*\Sum_{neigh} (A_n-A_i) + 
///	p3*\Sum_{neigh} (P_ni*A_n-P_in*A_i) 
///
///	dP_i/dt = p4 - p5*P_i 
///
///	P_in = P_i*X_in/(p_6+\Sum_{k,neigh}X_ik)
///
///     k_in = p_7 + p_8/(p_9+A_i)
///
/// where X_in is the stress in the wall (X_in = k_in*F/(k_in+k_ni).
/// The column indices for auxin and PIN in cells are at first level, 
/// and F k_1 k_2 in wall are at the second level of indices.
/// k_1 and k_2 are set via the cell auxin concentrations and 1,2 are the order of cell neighbors.

class AuxinModelSimpleStress : public BaseReaction {
  
 public:
  
  AuxinModelSimpleStress(std::vector<double> &paraValue, 
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

//!A cell-based auxin transport model including AUX1 and PID
/*!A complete pattern generating auxin model based on only cellular
  compartments. The four molecules are updated according to:
  
  dA_i/dt = p0*M_i p1 - p2*A_i +p5*\Sum_{neigh} (A_n-A_i) + 
  p4*\Sum_{neigh} (P_ni*A_n-P_in*A_i) 
  
  dP_i/dt = p6 - p7*P_i 
  
  dX_i/dt = p8*A_i - p9*X_i
  
  dM_i/dt = p10*\Theta_L1 - p11*M_i
  
  In addition, the column index for auxin, PIN, AUX1, PID, X, 
  and M should be given.
*/
class AuxinModelSimple2 : public BaseReaction {
  
 public:
  
  AuxinModelSimple2(std::vector<double> &paraValue, 
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

//!A cell-based auxin transport model including AUX1 and PID
/*!A complete pattern generating auxin model based on only cellular
  compartments. The four molecules are updated according to:
  
  dA_i/dt = p0*M_i p1 - p2*A_i +p5*\Sum_{neigh} (A_n-A_i) + 
  p4*\Sum_{neigh} (P_ni*A_n-P_in*A_i) 
  
  dP_i/dt = p6 - p7*P_i 
  
  dX_i/dt = p8*A_i - p9*X_i
  
  dM_i/dt = p10*\Theta_L1 - p11*M_i
  
  In addition, the column index for auxin, PIN, AUX1, PID, X, 
  and M should be given.
*/
class AuxinModelSimple3 : public BaseReaction {
  
 public:
  
  AuxinModelSimple3(std::vector<double> &paraValue, 
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

class AuxinModel4 : public BaseReaction {
  
 public:
  
  AuxinModel4(std::vector<double> &paraValue, 
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

class AuxinModel5 : public BaseReaction {
  
 public:
  
  AuxinModel5(std::vector<double> &paraValue, 
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

class AuxinModel6 : public BaseReaction {
  
 public:
  
  AuxinModel6(std::vector<double> &paraValue, 
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

class AuxinModel7 : public BaseReaction {
  
 public:
  
  AuxinModel7(std::vector<double> &paraValue, 
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

class AuxinTransportCellCellNoGeometry : public BaseReaction {
  
 public:
  
  AuxinTransportCellCellNoGeometry(std::vector<double> &paraValue, 
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
/// @brief A cell-wall based auxin transport model including PINs and ROPs
///
/// A complete pattern generating auxin model based on cell and wall
/// compartments. It uses two compartments for each wall and a single for the cells.
/// Auxin and PIN molecules are updated according to:
///  
/// @f[ \frac{A_i}{dt} = p_0 - p_1 A_i + \sum_{j} ( (p_2+p_3 [AUX]_i) A_{ij} ) 
/// - \sum_{j} (p_4+ p_5 P_{ij}) A_i @f]
///  
/// @f[ \frac{dA_{ij}}{dt} = (from above) + p_6 (A_{ji}-A_{ij}) @f]
///
/// @f[ \frac{dP_i}{dt} = p_7 - p_8 P_i + \sum_j (p_9 P_{ij} - (p_{10} + p_{11} A_{j}) P_i @f] 
///  
/// @f[ \frac{dP_{ij}}{dt} = (from above) @f]
///
/// The (here static) symmetric AUX contribution is applied via AUX in the cells.
///   
/// In the model file the reaction is given by:
/// @verbatim
/// AuxinWallModel 12 2 3 2
/// p_0 ... p_11
/// ci_auxin ci_PIN ci_AUX
/// wi_auxin wi_PIN
/// @endverbatim
///
class AuxinWallModel : public BaseReaction {
  
 public:
  
  AuxinWallModel(std::vector<double> &paraValue, 
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
/// @brief A cell-wall based auxin transport model including PINs and ROPs
///
/// A complete (hopefully) pattern generating auxin model based on cell and wall
/// compartments. It uses two compartments for each wall and a single for the cells.
/// Auxin PIN and ROP molecules are updated according to:
///  
/// @f[ \frac{A_i}{dt} = p_0 - p_1 A_i + p_2 \sum_{j} (A_{ij}) - p_3 \sum_{j} (A_i) - 
///  p_4 \sum_{j} (P_{ij} A_i) @f]
///  
/// @f[ \frac{dA_{ij}}{dt} = (from above) + p_5 (A_{ji}-A_{ij}) @f]
///
/// @f[ \frac{dP_i}{dt} = p_6 - p_7 P_i + \sum_j (p_8 P_{ij} - p_9 P_i R_{ij} @f] 
///  
/// @f[ \frac{dP_{ij}}{dt} = (from above) @f]
///
/// @f[ \frac{dR_{i}}{dt} = p_{10} – p_{11} R_i + \sum_j (p_{12} R_{ij} 
/// ((R_{ji})^(p_{15}))/((p_{14})^(p_{15}) + (R_{ji})^(p_{15}))) – p_{13} R_i A_{ij} @f]
///  
/// @f[ \frac{dR_{ij}}{dt} = (from above) @f]
///  
/// In the model file the reaction is given by:
/// @verbatim
/// AuxinROPModel 16 2 3 3
/// p_0 ... p_15
/// ci_auxin ci_PIN ci_ROP
/// wi_auxin wi_PIN wi_ROP
/// @endverbatim
///
class AuxinROPModel : public BaseReaction {
  
 public:
  
  AuxinROPModel(std::vector<double> &paraValue, 
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
/// A cell-wall based auxin transport model including PINs and ROPs.
///
/// Auxin model based on cell and wall compartments. It uses two compartments for  
/// each wall and a single for the cells. It contains several updates to the  
/// AuxinROPModel above; most notably using MM and Hill terms. Auxin, PIN and ROP 
/// molecules are updated according to:
///
/// @f[ \frac{A_i}{dt} = p_0 - p_1 A_i + p_2 \sum_{j} (A_{ij}) – p_3 \sum_{j}  
/// (A_i) – p_4 \sum_{j} (P_{ij} (A_i / (p_5 + A_i))) @f]
///
/// @f[ \frac{dA_{ij}}{dt} = (from above) + p_6 (A_{ji} – A_{ij}) @f]
///
/// @f[ \frac{dP_i}{dt} = p_7 – p_8 P_i + \sum_j (p_9 P_{ij} – p_{10} P_i    
/// (R_{ij})^{p_{12}}/(p_{11}^{p_{12}} + R_{ij}^{p_{12}}) @f]
///
/// @f[ \frac{dP_{ij}}{dt} = (from above) @f]
///
/// @f[ \frac{dR_{i}}{dt} = p_{13} – p_{14} R_i + \sum_j (p_{15} R_{ij} 
/// ((R_{ji})^(p_{17}))/((p_{16})^(p_{17}) + (R_{ji})^(p_{17}))) – p_{18} R_i A_{ij} @f]
///
/// @f[ \frac{dR_{ij}}{dt} = (from above) @f]
///
/// In the model file the reaction is given by:
/// @verbatim
/// AuxinROPModel2 19 2 3 3
/// p_0 ... p_18
/// ci_auxin ci_PIN ci_ROP
/// wi_auxin wi_PIN wi_ROP
/// @endverbatim
///
class AuxinROPModel2 : public BaseReaction {

  public:

  AuxinROPModel2(std::vector<double> &paraValue,
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
/// A cell-wall based auxin transport model including PINs and ROPs
///
/// A complete (hopefully) pattern generating auxin model based on cell and wall
/// compartments. It uses two compartments for each wall and a single for the cells.
/// Auxin PIN and ROP molecules are updated according to:
///  
/// @f[ \frac{A_i}{dt} = p_0 - p_1 A_i + p_2 \sum_{j} (A_{ij}) - p_3 \sum_{j} (A_i) - 
/// p_4 \sum_{j} (P_{ij} A_i) + p_{16} \sum_{j} (AUX_i A_{ij})@f]
///  
/// @f[ \frac{dA_{ij}}{dt} = (from above) + p_5 (A_{ji}-A_{ij}) @f]
///
/// @f[ \frac{dP_i}{dt} = p_6 - p_7 P_i + \sum_j (p_8 P_{ij} - p_9 P_i R_{ij} @f] 
///
/// @f[ \frac{dP_{ij}}{dt} = (from above) @f]
///
/// @f[ \frac{dR_i}{dt} = p_{10} - p_{11} R_i + \sum_j ( p_{12} R_{ij} - p_{13} R_i A_{ij}) @f]
///  
/// @f[ \frac{dR_{ij}}{dt} = (from above) @f]
///
/// The (here static) and symmetric AUX dependence is the only difference from AuxinROPModel, and is
/// via parameter p_16.  
///
/// In the model file the reaction is given by:
/// @verbatim
/// AuxinROPModel3 17 2 4 3
/// p_0 ... p_16
/// ci_auxin ci_PIN ci_ROP ci_AUX
/// wi_auxin wi_PIN wi_ROP
/// @endverbatim
///
class AuxinROPModel3 : public BaseReaction {
  
 public:
  
  AuxinROPModel3(std::vector<double> &paraValue, 
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
/// @brief A cell-wall based auxin transport model including PINs and an mutual repression
/// between PINs at neighboring membranes in adjacent cells
///
/// A complete (hopefully) pattern generating auxin model based on cell and wall
/// compartments. It uses two compartments for each wall and a single for the cells.
/// Auxin and PIN molecules are updated according to:
///  
/// @f[ \frac{A_i}{dt} = p_0 - p_1 A_i + p_2 \sum_{j} (A_{ij}) - p_3 \sum_{j} (A_i) - 
///  p_4 \sum_{j} (P_{ij} A_i) @f]
///  
/// @f[ \frac{dA_{ij}}{dt} = (from above) @f]
///
/// @f[ \frac{dP_i}{dt} = p_5 - p_6 P_i + \sum_j (1 - P_{ij}P_{ji}/2 + P_{ji}^2/2)P_{ij} - \sum_j P_i A_{ij} @f] 
///  
/// @f[ \frac{dP_{ij}}{dt} = (from above) @f]
///
/// In the model file the reaction is given by:
/// @verbatim
/// AuxinPINBistabilityModel 7 2 2 2
/// p_0 ... p_6
/// ci_auxin ci_PIN
/// wi_auxin wi_PIN
/// @endverbatim
///
class AuxinPINBistabilityModel : public BaseReaction {
  
 public:
  
  AuxinPINBistabilityModel(std::vector<double> &paraValue, 
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
