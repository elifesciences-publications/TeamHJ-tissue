
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

///
/// @brief A stress-based PIN1 and MT polarization model
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
/// @details A complete pattern generating auxin model based on only cellular
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
/// @f[ \frac{dM_i}{dt} = p_{10} \theta_{L1} - p_{11} M_i @f]
///  
/// @f[ P_{in} = \frac{P_i X_n}{(p_3 + \sum_{k}^{neigh} X_k)} @f]
///  
/// In a model file the reaction is defined as:
/// @verbatim
/// AuxinModelSimple1 12 1 4
/// p_0 ... p_11
/// A_index P_index X_index M_index
/// @endverbatim
/// or alternatively
/// @verbatim
/// AuxinModelSimple1 12 2 4 1
/// p_0 ... p_11
/// A_index P_index X_index M_index
/// P_wall (save index pair)
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

///
/// @brief A linear polarization from wall signal'-based auxin transport model
///
/// @details A complete auxin transport model based on PIN polarization from a linear
/// feedback from a wall variable (X). PIN and auxin are updated
/// according to:
///
/// @f[ \frac{dA_i}{dt} = p_0 - p_1 A_i + p_4 \sum_{neigh} (A_n-A_i) + @f] 
/// @f[ p_3 \sum_{neigh} (P_{ni} A_n - P_{in} A_i) @f] 
///
/// @f[ \frac{dP_i}{dt} = p_5 - p_6 P_i @f] 
///
/// @f[ P_{in} = P_i \frac{X_{in}}{(p_2 + \sum_{k,neigh} X_{ik})} @f]
///
/// where @f$X_{in}@f$ is the polarizing variable in the wall.
/// The column index for auxin, PIN, and X (in wall) should be given.
///
/// In a model file, the reaction is given as:
/// @verbatim
/// to come...
/// @endverbatim
///
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
/// @details An auxin transport model for cellular auxin and where PIN polarization is based on stresses
/// in a double wall compartment, i.e. the wall compartment between a cell pair is divided into
/// two compartments, and stress is measured in the wall compartment neighbor to the specific cell.
/// Auxin concentration is affecting the spring constants in wall pieces neighboring the cell.
///
/// PIN and auxin are updated according to:
///
/// @f[ \frac{dA_i}{dt} = p_0 - p_1 A_i +p_2 \sum_{neigh} (A_n-A_i) + @f] 
/// @f[ p_3 \sum_{neigh} (P_{ni} A_n - P_{in} A_i) @f] 
///
/// and
///
/// @f[ \frac{dP_i}{dt} = p_4 - p_5 P_i @f] 
///
/// where
///
/// @f[P_{in} = P_i \frac{X_{in}}{(p_6 + \sum_{k,neigh} X_{ik})} @f]
///
/// @f[ k_{in} = p_7 + \frac{p_8}{(p_9+A_i)} @f]
///
/// where @f$X_{in}@f$ is the stress in the wall (@f$X_{in} = k_{in} F/(k_{in}+k_{ni})@f$).
/// The column indices for auxin and PIN in cells are at first level, 
/// and @f$ F, k_1, k_2 @f$ in wall are at the second level of indices.
/// @f$ k_1 @f$ and @f$ k_2 @f$ are set via the cell auxin concentrations and 1,2 
/// are the order of cell neighbors.
///
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

///
/// @brief A cell-based auxin transport model including AUX1 and PID
///
/// @details A complete pattern generating auxin model based on only cellular
/// compartments. The four molecules are updated according to:
///  
/// @f[ \frac{dA_i}{dt} = p_0 M_i p_1 - p_2 A_i + p_5 \sum_{neigh} (A_n-A_i) + @f] 
/// @f[ p_4 \sum_{neigh} (P_{ni} A_n - P_{in} A_i) @f] 
///
/// @f[ \frac{dP_i}{dt} = p_6 - p_7 P_i @f] 
///  
/// @f[ \frac{dX_i}{dt} = p_8 A_i - p_9 X_i @f]
///
/// @f[ \frac{dM_i}{dt} = p_{10} \theta_{L1} - p_{11} M_i @f]
///
/// In addition, the column index for auxin, PIN, AUX1, PID, X, 
/// and M should be given in a model file.
///
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

///
/// @brief A cell-based auxin transport model including AUX1 and PID
///
/// @details A complete pattern generating auxin model based on only cellular
/// compartments. The four molecules are updated according to:
///
/// @f[ \frac{dA_i}{dt} = p_0 M_i p_1 - p_2 A_i + p_5 \sum_{neigh} (A_n-A_i) + @f] 
/// @f[ p_4 \sum_{neigh} (P_{ni} A_n - P_{in} A_i) @f] 
///
/// @f[ \frac{dP_i}{dt} = p_6 - p_7 P_i @f] 
///
/// @f[ \frac{dX_i}{dt} = p_8 A_i - p_9 X_i @f]
///
/// @f[ \frac{dM_i}{dt} = p_{10} \theta_{L1} - p_{11} M_i @f]
///
/// In addition to the parameter values, the column index for auxin, PIN, AUX1, PID, X, 
/// and M should be given in the model file.
///
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
/// @brief A cell-wall based auxin transport model including PINs and AUXs
///
/// A complete pattern generating auxin model based on cell and wall
/// compartments. It uses two compartments for each wall and a single for the cells.
/// Auxin and PIN molecules are updated according to:
///  
/// @f[ \frac{dA_i}{dt} = p_0 - p_1 A_i + \sum_{j} ( (p_2+p_3 [AUX]_i) A_{ij} ) 
/// - \sum_{j} (p_4+ p_5 P_{ij}) A_i @f]
///  
/// @f[ \frac{dA_{ij}}{dt} = (from above) + p_6 (A_{ji}-A_{ij}) @f]
///
/// @f[ \frac{dP_i}{dt} = p_7 - p_8 P_i + \sum_j (p_9 P_{ij} - (p_{10} + p_{11} X_{j}) P_i) @f] 
///  
/// @f[ \frac{dP_{ij}}{dt} = (from above) @f]
///
/// @f[ \frac{dX_i}{dt} = p_{12} A_{i} - p_{13} X_{i} @f]
///
/// The (here static) symmetric AUX contribution is applied via AUX in the cells.
/// Polarization feedback comes from the molecule X that is activated by auxin.
///   
/// In the model file the reaction is given by:
/// @verbatim
/// AuxinWallModel 14 2 4 2
/// p_0 ... p_11
/// ci_auxin ci_PIN ci_AUX ci_X
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
/// @f[ \frac{dA_i}{dt} = p_0 - p_1 A_i + p_2 \sum_{j} (A_{ij}) - p_3 \sum_{j} (A_i) - 
///  p_4 \sum_{j} (P_{ij} A_i) @f]
///  
/// @f[ \frac{dA_{ij}}{dt} = (from above) + p_5 (A_{ji}-A_{ij}) @f]
///
/// @f[ \frac{dP_i}{dt} = p_6 - p_7 P_i + \sum_j (p_8 P_{ij} - p_9 P_i R_{ij} @f] 
///  
/// @f[ \frac{dP_{ij}}{dt} = (from above) @f]
///
/// @f[ \frac{dR_{i}}{dt} = p_{10} – p_{11} R_i + @f]
/// @f[ \sum_j (p_{12} R_{ij} \frac{R_{ji}^{p_{15}}}{p_{14}^{p_{15}} + R_{ji}^{p_{15}}} @f]
/// @f[ - p_{13} R_i A_{ij} ) @f]
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
/// @brief A cell-wall based auxin transport model including PINs and ROPs.
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
/// @brief A cell-wall based auxin transport model including PINs and ROPs
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

///
/// @brief A cell-wall based auxin transport model including PINs and an mutual repression
/// between PINs at neighboring membranes in adjacent cells
///
/// A complete (hopefully) pattern generating auxin model based on cell and wall
/// compartments. It uses two compartments for each wall and a single for the cells.
/// Auxin and PIN molecules are updated according to:
///  
/// @f[ \frac{dA_i}{dt} = p_0 - p_1 A_i - p_2 \sum_{j} (A_i-A_j) - 
///  p_4 \sum_{j} (P_{ij} A_i - P_{ji} A_j) @f]
///  
/// @f[ \frac{dP_i}{dt} = p_5 + p_6 A_i - p_7 P_i + \sum_j (1 - P_{ij}P_{ji}/2 + P_{ji}^2/2)P_{ij}^+ - p_3 \sum_j P_i @f] 
///  
/// @f[ \frac{dP_{ij}}{dt} = (from above) @f]
///
/// In the model file the reaction is given by:
/// @verbatim
/// AuxinPINBistabilityModelCell 8 2 2 1
/// p_0 ... p_7
/// ci_auxin ci_PIN
/// wi_PIN
/// @endverbatim
///
class AuxinPINBistabilityModelCell : public BaseReaction {
  
 public:
  
  AuxinPINBistabilityModelCell(std::vector<double> &paraValue, 
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
/// @brief A cell-wall based auxin transport model including PINs on walls, an mutual repression
/// between exocytosis rate of PINs at neighboring membranes in adjacent cells and PIN membrane diffusion
///
/// This model was previously named auxinPINBistabilityModelCellNew (this name should still run the reaction)
/// A complete (hopefully) pattern generating auxin model based on cell and wall
/// compartments. It uses two compartments for each wall and a single for the cells.
/// Auxin and PIN molecules and PIN exocytosis rates are updated according to:
///  
/// @f[ \frac{dA_i}{dt} =  p_0 \sum_{j} (A_i-A_j) - 
///  p_1 \sum_{j} (P_{ij} A_i - P_{ji} A_j)+p_4 (1-A_i) @f]
///  
/// @f[ \frac{dP_i}{dt} =\sum_j (p_2 (P_{ij}-R_{ij}P_i))+p_5 (A_i-P_i) @f] 
///  
/// @f[ \frac{dP_{ij}}{dt} = (from above) - \sum_k p_3 (p_{ik}-p_{ij}) @f]
///
/// @f[ \frac{dR_{ij}}{dt} = R_X(0.5*(A_i+A_j))-R_{ij}+\frac{R_{ij}-R_{ji}}{2}\frac{R_{ij}R_{ji}}{p_8} @f]
///
/// @f[ R_X(x) = \frac{p_6 x}{p_7+ x} @f]
///
/// In the model file the reaction is given by:
/// @verbatim
/// AuxinExoBistability 9 2 2 1
/// p_0 ... p_8
/// ci_auxin ci_PIN
/// wi_PIN wi_Exo
/// @endverbatim
///
/// @note PIN is allowed to diffuse in the membrane.
///
class AuxinExoBistability : public BaseReaction {
  
 public:
  
  AuxinExoBistability(std::vector<double> &paraValue, 
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
/// between PINs at neighboring membranes in adjacent cells auxin production and degradation occurs in source or sink cells, respectively
///
/// A complete (hopefully) pattern generating auxin model based on cell and wall
/// compartments. It uses two compartments for each wall and a single for the cells.
/// Auxin and PIN molecules are updated according to:
///  
/// @f[ \frac{dA_i}{dt} = p_0 - p_1 A_i - p_2 \sum_{j} (A_i-A_j) - 
///  p_4 \sum_{j} (P_{ij} A_i - P_{ji} A_j) @f]
///  
/// @f[ \frac{dP_i}{dt} = p_5 + p_6 A_i - p_7 P_i + \sum_j (1 - P_{ij}P_{ji}/2 + P_{ji}^2/2)P_{ij}^+ - p_3 \sum_j P_i @f] 
///  
/// @f[ \frac{dP_{ij}}{dt} = (from above) @f]
///
/// In the model file the reaction is given by:
/// @verbatim
/// AuxinPINBistabilityModelCellSourceSink 8 2 4 1
/// p_0 ... p_7
/// ci_auxin ci_PIN ci_Source ci_Sink
/// wi_PIN
/// @endverbatim
///
class AuxinPINBistabilityModelCellSourceSink : public BaseReaction {
  
 public:
  
  AuxinPINBistabilityModelCellSourceSink(std::vector<double> &paraValue, 
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
/// @brief A cell-wall based auxin transport model including PINs and AUX1 (L) and an mutual repression
/// between PINs at neighboring membranes in adjacent cells auxin production and degradation occurs in source or sink cells, respectively
///
/// A complete (hopefully) pattern generating auxin model based on cell and wall
/// compartments. It uses two compartments for each wall and a single for the cells.
/// Auxin and PIN molecules are updated according to:
///  
/// @f[ \frac{dA_i}{dt} = p_0 - p_1 A_i - p_2 \sum_{j} (A_i-A_j) - 
///  p_4 \sum_{j} (P_{ij} A_i L_j/(L_i+L_j) - P_{ji} A_j  L_i/(L_i+L_j)) @f]
///  
/// @f[ \frac{dP_i}{dt} = p_5 + p_6 A_i - p_7 P_i + \sum_j (1 - P_{ij}P_{ji}/2 + P_{ji}^2/2)P_{ij}^+ - p_3 \sum_j P_i @f] 
///  
/// @f[ \frac{dP_{ij}}{dt} = (from above) @f]
///
/// In the model file the reaction is given by:
/// @verbatim
/// AuxinPINBistabilityInfluxModel 8 2 5 1
/// p_0 ... p_7
/// ci_auxin ci_PIN ci_Source ci_Sink cI_AUX
/// wi_PIN
/// @endverbatim
///
class AuxinPINBistabilityInfluxModel : public BaseReaction {
  
 public:
  
  AuxinPINBistabilityInfluxModel(std::vector<double> &paraValue, 
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
/// @brief A cell-wall based auxin transport model including PINs with cross membrane interaction (inspired by ROPs). Here PIN exocytosis also depends on auxin in neighbouring wall compartment.
///
/// A complete (hopefully) pattern generating auxin model based on cell and wall
/// compartments. It uses two compartments for each wall and a single for the cells.
/// Auxin and PIN  molecules are updated according to:
///  
///  
/// @f[ \frac{dA_i}{dt} = p_0 - p_1 A_i +\sum_{j}(p_2A_{ij}- p_3 A_i) - 
///  p_4 \sum_{j} (P_{ij} A_i) @f]
///  
/// @f[ \frac{dA_{ij}}{dt} = (from above) + p_5 (A_{ji}-A_{ij}) @f]
///
/// @f[ \frac{dP_i}{dt} = p_6 - p_7 P_i + \sum_j (p_8 P_{ij}\frac{P_{ji}^{p_9}}{p_{10}^{p_9}+P_{ji}^{p_9}}+p_{11}P_{ij} - p_{12} P_i A_{ij} @f] 
///  
/// @f[ \frac{dP_{ij}}{dt} = (from above) @f]
///
///In this model PIN in neighbouring cell acts by enhancing PIN endocytosis.
///  
/// In the model file the reaction is given by:
/// @verbatim
/// SimpleROPModel 13 2 2 2
/// p_0 ... p_12
/// ci_auxin ci_PIN 
/// wi_auxin wi_PIN 
/// @endverbatim
///
class SimpleROPModel : public BaseReaction {
  
 public:
  
  SimpleROPModel(std::vector<double> &paraValue, 
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
/// @brief A cell-wall based auxin transport model including PINs with cross membrane interaction (inspired by ROPs). Here PIN exocytosis does not depend on auxin in neighbouring wall compartment.
///
/// A complete (hopefully) pattern generating auxin model based on cell and wall
/// compartments. It uses two compartments for each wall and a single for the cells.
/// Auxin and PIN  molecules are updated according to:
///  
///  
/// @f[ \frac{dA_i}{dt} = p_0 - p_1 A_i +\sum_{j}(p_2A_{ij}- p_3 A_i) - 
///  p_4 \sum_{j} (P_{ij} A_i) @f]
///  
/// @f[ \frac{dA_{ij}}{dt} = (from above) + p_5 (A_{ji}-A_{ij}) @f]
///
/// @f[ \frac{dP_i}{dt} = p_6 - p_7 P_i + \sum_j (p_8 P_{ij}\frac{P_{ji}^{p_9}}{p_{10}^{p_9}+P_{ji}^{p_9}}+p_{11}P_{ij} - p_{12} P_i @f] 
///  
/// @f[ \frac{dP_{ij}}{dt} = (from above) @f]
///
///In this model PIN in neighbouring cell acts by enhancing PIN endocytosis.
///  
/// In the model file the reaction is given by:
/// @verbatim
/// SimpleROPModel2 13 2 2 2
/// p_0 ... p_12
/// ci_auxin ci_PIN 
/// wi_auxin wi_PIN 
/// @endverbatim
///
class SimpleROPModel2 : public BaseReaction {
  
 public:
  
  SimpleROPModel2(std::vector<double> &paraValue, 
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
/// @brief A cell-wall based auxin transport model including PINs with cross membrane interaction (inspired by ROPs). Here PIN exocytosis does not depend on auxin in neighbouring wall compartment.
///
/// A complete (hopefully) pattern generating auxin model based on cell and wall
/// compartments. It uses two compartments for each wall and a single for the cells.
/// Auxin and PIN  molecules are updated according to:
///  
///  
/// @f[ \frac{dA_i}{dt} = p_0 - p_1 A_i +\sum_{j}(p_2A_{ij}- p_3 A_i) - 
///  p_4 \sum_{j} (P_{ij} A_i) @f]
///  
/// @f[ \frac{dA_{ij}}{dt} = (from above) + p_5 (A_{ji}-A_{ij}) @f]
///
/// @f[ \frac{dP_i}{dt} = p_6 - p_7 P_i + \sum_j (p_8 P_{ij} - p_{11} P_i\frac{p_{10}^{p_9}}{p_{10}^{p_9}+P_{ji}^{p_9}}) @f] 
///  
/// @f[ \frac{dP_{ij}}{dt} = (from above) @f]
///
///In this model PIN in neighbouring cell acts by repressing PIN exocytosis.
///  
/// In the model file the reaction is given by:
/// @verbatim
/// SimpleROPModel3 12 2 2 2
/// p_0 ... p_11
/// ci_auxin ci_PIN 
/// wi_auxin wi_PIN 
/// @endverbatim
///
class SimpleROPModel3 : public BaseReaction {
  
 public:
  
  SimpleROPModel3(std::vector<double> &paraValue, 
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
/// @brief A cell-wall based auxin transport model including PINs with cross membrane interaction (inspired by ROPs). Here PIN exocytosis does not depend on auxin in neighbouring wall compartment.
///
/// A complete (hopefully) pattern generating auxin model based on cell and wall
/// compartments. It uses two compartments for each wall and a single for the cells.
/// Auxin and PIN  molecules are updated according to:
///  
///  
/// @f[ \frac{dA_i}{dt} = p_0 - p_1 A_i +\sum_{j}(p_2A_{ij}- p_3 A_i) - 
///  p_4 \sum_{j} (P_{ij} A_i) @f]
///  
/// @f[ \frac{dA_{ij}}{dt} = (from above) + p_5 (A_{ji}-A_{ij}) @f]
///
/// @f[ \frac{dP_i}{dt} = p_6 - p_7 P_i + \sum_j (p_8 P_{ij} - p_{11} P_i\frac{p_{10}^{p_9}}{p_{10}^{p_9}+P_{ji}^{p_9}}) @f] 
///  
/// @f[ \frac{dP_{ij}}{dt} = (from above) @f]
///
///In this model PIN in neighbouring cell acts by repressing PIN exocytosis.
///  
/// In the model file the reaction is given by:
/// @verbatim
/// SimpleROPModel3 12 2 2 2
/// p_0 ... p_11
/// ci_auxin ci_PIN 
/// wi_auxin wi_PIN 
/// @endverbatim
///
class SimpleROPModel4 : public BaseReaction {
  
 public:
  
  SimpleROPModel4(std::vector<double> &paraValue, 
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
/// @brief A cell-wall based auxin transport model including PINs with cross membrane interaction (inspired by ROPs). Here PIN exocytosis does not depend on auxin in neighbouring wall compartment.
///


/// 2D grid with reactions in cell 4 only.
/// A complete (hopefully) pattern generating auxin model based on cell and wall
/// compartments. It uses two compartments for each wall and a single for the cells.
/// Auxin and PIN  molecules are updated according to:
///  
///  
/// @f[ \frac{dA_i}{dt} = p_0 - p_1 A_i +\sum_{j}(p_2A_{ij}- p_3 A_i) - 
///  p_4 \sum_{j} (P_{ij} A_i) @f]
///  
/// @f[ \frac{dA_{ij}}{dt} = (from above) + p_5 (A_{ji}-A_{ij}) @f]
///
/// @f[ \frac{dP_i}{dt} = p_6 - p_7 P_i + \sum_j (p_8 P_{ij} - p_{11} P_i\frac{p_{10}^{p_9}}{p_{10}^{p_9}+P_{ji}^{p_9}}) @f] 
///  
/// @f[ \frac{dP_{ij}}{dt} = (from above) @f]
///
///In this model PIN in neighbouring cell acts by repressing PIN exocytosis.
///  
/// In the model file the reaction is given by:
/// @verbatim
/// SimpleROPModel5 16 2 2 2
/// p_0 ... p_11
/// ci_auxin ci_PIN 
/// wi_auxin wi_PIN 
/// @endverbatim
///
class SimpleROPModel5 : public BaseReaction {
  
 public:
  
  SimpleROPModel5(std::vector<double> &paraValue, 
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
/// @brief A cell-wall based auxin transport model including PINs with cross membrane interaction (inspired by ROPs). Here PIN exocytosis does not depend on auxin in neighbouring wall compartment.
///


/// 2D grid with reactions in cell 4 only.
/// A complete (hopefully) pattern generating auxin model based on cell and wall
/// compartments. It uses two compartments for each wall and a single for the cells.
/// Auxin and PIN  molecules are updated according to:
///  
///  
/// @f[ \frac{dA_i}{dt} = p_0 - p_1 A_i +\sum_{j}(p_2A_{ij}- p_3 A_i) - 
///  p_4 \sum_{j} (P_{ij} A_i) @f]
///  
/// @f[ \frac{dA_{ij}}{dt} = (from above) + p_5 (A_{ji}-A_{ij}) @f]
///
/// @f[ \frac{dP_i}{dt} = p_6 - p_7 P_i + \sum_j (p_8 P_{ij} - p_{11} P_i\frac{p_{10}^{p_9}}{p_{10}^{p_9}+P_{ji}^{p_9}}) @f] 
///  
/// @f[ \frac{dP_{ij}}{dt} = (from above) @f]
///
///In this model PIN in neighbouring cell acts by repressing PIN exocytosis.
///  
/// In the model file the reaction is given by:
/// @verbatim
/// SimpleROPModel5 16 2 2 2
/// p_0 ... p_11
/// ci_auxin ci_PIN 
/// wi_auxin wi_PIN 
/// @endverbatim
///
class SimpleROPModel6 : public BaseReaction {
  
 public:
  
  SimpleROPModel6(std::vector<double> &paraValue, 
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

//
/// @brief A cell-wall based auxin transport model including PINs with cross membrane interaction (inspired by ROPs). Here PIN exocytosis does not depend on auxin in neighbouring wall compartment.
///


/// 2D grid with reactions in cell 4 only.
/// A complete (hopefully) pattern generating auxin model based on cell and wall
/// compartments. It uses two compartments for each wall and a single for the cells.
/// Auxin and PIN  molecules are updated according to:
///  
///  
/// @f[ \frac{dA_i}{dt} = p_0 - p_1 A_i +\sum_{j}(p_2A_{ij}- p_3 A_i) - 
///  p_4 \sum_{j} (P_{ij} A_i) @f]
///  
/// @f[ \frac{dA_{ij}}{dt} = (from above) + p_5 (A_{ji}-A_{ij}) @f]
///
/// @f[ \frac{dP_i}{dt} = p_6 - p_7 P_i + \sum_j (p_8 P_{ij} - p_{11} P_i\frac{p_{10}^{p_9}}{p_{10}^{p_9}+P_{ji}^{p_9}}) @f] 
///  
/// @f[ \frac{dP_{ij}}{dt} = (from above) @f]
///
///In this model PIN in neighbouring cell acts by repressing PIN exocytosis.
///  
/// In the model file the reaction is given by:
/// @verbatim
/// SimpleROPModel5 16 2 2 2
/// p_0 ... p_11
/// ci_auxin ci_PIN 
/// wi_auxin wi_PIN 
/// @endverbatim
///
class SimpleROPModel7 : public BaseReaction {
  
 public:
  
  SimpleROPModel7(std::vector<double> &paraValue, 
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
/// @details A complete pattern generating auxin model based on only cellular
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
/// @f[ \frac{dM_i}{dt} = p_{10} \theta_{L1} - p_{11} M_i @f]
///  
/// @f[ P_{in} = \frac{P_i X_n}{(p_3 + \sum_{k}^{neigh} X_k)} @f]
///  
/// In a model file the reaction is defined as:
/// @verbatim
/// AuxinModelSimple4 12 1 4
/// p_0 ... p_11
/// A_index P_index X_index M_index
/// @endverbatim
/// or alternatively
/// @verbatim
/// AuxinModelSimple4 12 2 4 1
/// p_0 ... p_11
/// A_index P_index X_index M_index
/// P_wall (save wall index)
/// @endverbatim
///
/// Note that this reaction is the same as AuxinModelSimple1, except that this only saves one PIN concentration to each wall compartment
///
///
class AuxinModelSimple4 : public BaseReaction {
  
 public:
  
  AuxinModelSimple4(std::vector<double> &paraValue, 
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
/// @details A complete pattern generating auxin model based on only cellular
/// compartments. The four molecules A(uxin), P(IN), X(auxin induced
/// molecule) are updated according to:
///  
/// @f[ \frac{dA_i}{dt} = p_1 - p_2 A_i + p_4 \sum_{n}^{neigh} (P_{ni} A_n - P_{in} A_i) + p_5 \sum_{n}^{neigh} (A_n-A_i) @f] 
///  
/// @f[ \frac{dP_i}{dt} = p_6 - p_7 P_i @f] 
///  
/// @f[ \frac{dX_i}{dt} = p_8 A_i - p_9 X_i @f]
///  
/// @f[ P_{in} = \frac{P_i f( X_n )}{(p_3 + \sum_{k}^{neigh} f( X_k))} @f]
///  
/// In a model file the reaction is defined as:
/// @verbatim
/// AuxinModelSimple1 10 1 3
/// p_0 ... p_9
/// A_index P_index X_index
/// @endverbatim
/// or alternatively
/// @verbatim
/// AuxinModelSimple1 10 2 3 1
/// p_0 ... p_9
/// A_index P_index X_index
/// P_wall (save index pair)
/// @endverbatim
///
class AuxinModelSimple5 : public BaseReaction {
  
 public:
  
  AuxinModelSimple5(std::vector<double> &paraValue, 
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
/// @brief A cell-wall based auxin transport model including PINs with down the inernal gradient. Here PIN exocytosis does depend on auxin in neighbouring wall compartment.
/// Documantation to follow

class UpInternalGradientModel : public BaseReaction {
  
 public:
  
 UpInternalGradientModel(std::vector<double> &paraValue, 
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
/// @brief A cell-wall based auxin transport model including PINs with down the inernal gradient. Here PIN exocytosis does depend on auxin in neighbouring wall compartment.
/// Documantation to follow

class DownInternalGradientModel : public BaseReaction {
  
 public:
  
 DownInternalGradientModel(std::vector<double> &paraValue, 
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
/// @brief A cell-wall based auxin transport model including PINs with down the inernal gradient. Here PIN exocytosis does depend on auxin in neighbouring wall compartment.
/// Documantation to follow

class UpExternalGradientModel : public BaseReaction {
  
 public:
  
 UpExternalGradientModel(std::vector<double> &paraValue, 
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
/// @brief A cell-wall based auxin transport model including PINs with down the inernal gradient. Here PIN exocytosis does depend on auxin in neighbouring wall compartment.
/// Documantation to follow

class DownInternalGradientModelSingleCell : public BaseReaction {
  
 public:
  
 DownInternalGradientModelSingleCell(std::vector<double> &paraValue, 
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
/// @brief A cell-wall based auxin transport model including PINs with cross membrane interaction (inspired by ROPs). Here PIN exocytosis does depend on auxin in neighbouring wall compartment.
///
/// A complete (hopefully) pattern generating auxin model based on cell and wall
/// compartments. It uses two compartments for each wall and a single for the cells.
/// Auxin and PIN  molecules are updated according to:
///  
///  
/// @f[ \frac{dA_i}{dt} = p_0 - p_1 A_i +\sum_{j}(p_2A_{ij}- p_3 A_i) - 
///  p_4 \sum_{j} (P_{ij} A_i) @f]
///  
/// @f[ \frac{dA_{ij}}{dt} = (from above) + p_5 (A_{ji}-A_{ij}) @f]
///
/// @f[ \frac{dP_i}{dt} = p_6 - p_7 P_i*A_ij + \sum_j (p_8 P_{ij} - p_{11} P_i(Ai Pij-Aj Pji)) @f] 
///  
/// @f[ \frac{dP_{ij}}{dt} = (from above) @f]
///
/// In this model PIN in neighbouring cell acts by repressing PIN exocytosis.
///  
/// In the model file the reaction is given by:
/// @verbatim
/// AuxinFluxModel 12 2 2 2
/// p_0 ... p_11
/// ci_auxin ci_PIN 
/// wi_auxin wi_PIN 
/// @endverbatim
///
class AuxinFluxModel : public BaseReaction {
  
 public:
  
  AuxinFluxModel(std::vector<double> &paraValue, 
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
/// @brief A polarity model based on Abley et all Development 2013.
///
/// A and B  molecules are updated according to(needs editing):
///  
///  
/// @f[ \frac{dA_i}{dt} = p_0 - p_1 A_i +\sum_{j}(p_2A_{ij}- p_3 A_i) - 
///  p_4 \sum_{j} (P_{ij} A_i) @f]
///  
/// @f[ \frac{dA_{ij}}{dt} = (from above) + p_5 (A_{ji}-A_{ij}) @f]
///
/// @f[ \frac{dB_i}{dt} = p_6 - p_7 P_i*A_ij + \sum_j (p_8 P_{ij} - p_{11} P_i*(Ai*Pij-Aj*Pji)) @f] 
///  
/// @f[ \frac{dB_{ij}}{dt} = (from above) @f]
///
///In this model PIN in neighbouring cell acts by repressing PIN exocytosis.
///  
/// In the model file the reaction is given by:
/// @verbatim
/// IntracellularPartitioning 5 2 2 2
/// p_0 ... p_11
/// ci_auxin ci_PIN 
/// wi_auxin wi_PIN 
/// @endverbatim
///
class IntracellularPartitioning : public BaseReaction {
  
 public:
  
  IntracellularPartitioning(std::vector<double> &paraValue, 
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
/// @brief A polarity model based on Abley et all Development 2013.
///
/// A and B  molecules are updated according to(needs editing):
///  
/// @f[ \frac{dA_i}{dt} = p_0 - p_1 A_i +\sum_{j}(p_2A_{ij}- p_3 A_i) - 
///  p_4 \sum_{j} (P_{ij} A_i) @f]
///  
/// @f[ \frac{dA_{ij}}{dt} = (from above) + p_5 (A_{ji}-A_{ij}) @f]
///
/// @f[ \frac{dP_i}{dt} = p_6 - p_7 P_i*A_ij + \sum_j (p_8 P_{ij} - p_{11} P_i\frac{p_{10}^{p_9}}{p_{10}^{p_9}+P_{ji}^{p_9}}) @f] 
///  
/// @f[ \frac{dP_{ij}}{dt} = (from above) @f]
///
///In this model PIN in neighbouring cell acts by repressing PIN exocytosis.
///  
/// In the model file the reaction is given by:
/// @verbatim
///IntracellularCoupling 6 2 2 2
/// p_0 ... p_11
/// ci_auxin ci_PIN 
/// wi_auxin wi_PIN 
/// @endverbatim
///
class IntracellularCoupling : public BaseReaction {
  
 public:
  
  IntracellularCoupling(std::vector<double> &paraValue, 
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
/// @brief A polarity model based on Abley et all Development 2013.
///
/// A and B  molecules are updated according to(needs editing):
///  
///  
/// @f[ \frac{dA_i}{dt} = p_0 - p_1 A_i +\sum_{j}(p_2A_{ij}- p_3 A_i) - 
///  p_4 \sum_{j} (P_{ij} A_i) @f]
///  
/// @f[ \frac{dA_{ij}}{dt} = (from above) + p_5 (A_{ji}-A_{ij}) @f]
///
/// @f[ \frac{dP_i}{dt} = p_6 - p_7 P_i*A_ij + \sum_j (p_8 P_{ij} - p_{11} P_i\frac{p_{10}^{p_9}}{p_{10}^{p_9}+P_{ji}^{p_9}}) @f] 
///  
/// @f[ \frac{dP_{ij}}{dt} = (from above) @f]
///
///In this model PIN in neighbouring cell acts by repressing PIN exocytosis.
///  
/// In the model file the reaction is given by:
/// @verbatim
///  IntracellularIndirectCoupling 12 2 3 3
/// p_0 ... p_11
/// ci_auxin ci_PIN 
/// wi_auxin wi_PIN 
/// @endverbatim
///
class IntracellularIndirectCoupling : public BaseReaction {
  
 public:
  
  IntracellularIndirectCoupling(std::vector<double> &paraValue, 
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
/// @brief A cell-wall based auxin transport model including PINs with down the internal gradient. Here PIN exocytosis does depend on auxin in neighbouring wall compartment. includes gemoetric considerations
/// Documantation to follow

class DownInternalGradientModelGeometric : public BaseReaction {
  
 public:
  
 DownInternalGradientModelGeometric(std::vector<double> &paraValue, 
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






