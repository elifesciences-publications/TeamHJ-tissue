/**
 * Filename     : network.h
 * Description  : Classes describing complete models updates
 * Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
 * Created      : November 2006
 * Revision     : $Id:$
 */
#ifndef NETWORK_H
#define NETWORK_H

#include<cmath>

#include"tissue.h"
#include"baseReaction.h"

//!A stress-based PIN1 and MT polarization model
class AuxinModelStress : public BaseReaction {
  
 public:
  
  AuxinModelStress(std::vector<double> &paraValue, 
		   std::vector< std::vector<size_t> > 
		   &indValue );
  
  void derivs(Tissue &T,
	      std::vector< std::vector<double> > &cellData,
	      std::vector< std::vector<double> > &wallData,
	      std::vector< std::vector<double> > &vertexData,
	      std::vector< std::vector<double> > &cellDerivs,
	      std::vector< std::vector<double> > &wallDerivs,
	      std::vector< std::vector<double> > &vertexDerivs );
};

//!A linear polarization cell-based auxin transport model
/*!A complete pattern generating auxin model based on only cellular
  compartments. The four molecules are updated according to:
  
  dA_i/dt = p0*M_i + p1 - p2*A_i +p5*\Sum_{neigh} (A_n-A_i) + 
  p4*\Sum_{neigh} (P_ni*A_n-P_in*A_i) 
  
  dP_i/dt = p6 - p7*P_i 
  
  dX_i/dt = p8*A_i - p9*X_i
  
  dM_i/dt = p10*\Theta_L1 - p11*M_i
  
  P_in = P_i*X_n/(p_3+\Sum_{k,neigh}X_k)
  
  In addition, the column index for auxin, PIN, X, and M should be given.
*/
class AuxinModelSimple1 : public BaseReaction {
  
 public:
  
  AuxinModelSimple1(std::vector<double> &paraValue, 
		    std::vector< std::vector<size_t> > 
		    &indValue );
  
  void derivs(Tissue &T,
	      std::vector< std::vector<double> > &cellData,
	      std::vector< std::vector<double> > &wallData,
	      std::vector< std::vector<double> > &vertexData,
	      std::vector< std::vector<double> > &cellDerivs,
	      std::vector< std::vector<double> > &wallDerivs,
	      std::vector< std::vector<double> > &vertexDerivs );
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
	      std::vector< std::vector<double> > &cellData,
	      std::vector< std::vector<double> > &wallData,
	      std::vector< std::vector<double> > &vertexData,
	      std::vector< std::vector<double> > &cellDerivs,
	      std::vector< std::vector<double> > &wallDerivs,
	      std::vector< std::vector<double> > &vertexDerivs );
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
	      std::vector< std::vector<double> > &cellData,
	      std::vector< std::vector<double> > &wallData,
	      std::vector< std::vector<double> > &vertexData,
	      std::vector< std::vector<double> > &cellDerivs,
	      std::vector< std::vector<double> > &wallDerivs,
	      std::vector< std::vector<double> > &vertexDerivs );
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
	      std::vector< std::vector<double> > &cellData,
	      std::vector< std::vector<double> > &wallData,
	      std::vector< std::vector<double> > &vertexData,
	      std::vector< std::vector<double> > &cellDerivs,
	      std::vector< std::vector<double> > &wallDerivs,
	      std::vector< std::vector<double> > &vertexDerivs );
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
	      std::vector< std::vector<double> > &cellData,
	      std::vector< std::vector<double> > &wallData,
	      std::vector< std::vector<double> > &vertexData,
	      std::vector< std::vector<double> > &cellDerivs,
	      std::vector< std::vector<double> > &wallDerivs,
	      std::vector< std::vector<double> > &vertexDerivs );
};

class AuxinModel4 : public BaseReaction {
  
 public:
  
  AuxinModel4(std::vector<double> &paraValue, 
	      std::vector< std::vector<size_t> > 
	      &indValue );
  
  void derivs(Tissue &T,
	      std::vector< std::vector<double> > &cellData,
	      std::vector< std::vector<double> > &wallData,
	      std::vector< std::vector<double> > &vertexData,
	      std::vector< std::vector<double> > &cellDerivs,
	      std::vector< std::vector<double> > &wallDerivs,
	      std::vector< std::vector<double> > &vertexDerivs );
};

class AuxinModel5 : public BaseReaction {
  
 public:
  
  AuxinModel5(std::vector<double> &paraValue, 
	      std::vector< std::vector<size_t> > 
	      &indValue );
  
  void derivs(Tissue &T,
	      std::vector< std::vector<double> > &cellData,
	      std::vector< std::vector<double> > &wallData,
	      std::vector< std::vector<double> > &vertexData,
	      std::vector< std::vector<double> > &cellDerivs,
	      std::vector< std::vector<double> > &wallDerivs,
	      std::vector< std::vector<double> > &vertexDerivs );
};

class AuxinModel6 : public BaseReaction {
  
 public:
  
  AuxinModel6(std::vector<double> &paraValue, 
	      std::vector< std::vector<size_t> > 
	      &indValue );
  
  void derivs(Tissue &T,
	      std::vector< std::vector<double> > &cellData,
	      std::vector< std::vector<double> > &wallData,
	      std::vector< std::vector<double> > &vertexData,
	      std::vector< std::vector<double> > &cellDerivs,
	      std::vector< std::vector<double> > &wallDerivs,
	      std::vector< std::vector<double> > &vertexDerivs );
};

class AuxinModel7 : public BaseReaction {
  
 public:
  
  AuxinModel7(std::vector<double> &paraValue, 
	      std::vector< std::vector<size_t> > 
	      &indValue );
  
  void derivs(Tissue &T,
	      std::vector< std::vector<double> > &cellData,
	      std::vector< std::vector<double> > &wallData,
	      std::vector< std::vector<double> > &vertexData,
	      std::vector< std::vector<double> > &cellDerivs,
	      std::vector< std::vector<double> > &wallDerivs,
	      std::vector< std::vector<double> > &vertexDerivs );
};

class AuxinTransportCellCellNoGeometry : public BaseReaction {
  
 public:
  
  AuxinTransportCellCellNoGeometry(std::vector<double> &paraValue, 
				   std::vector< std::vector<size_t> > 
				   &indValue );
  
  void derivs(Tissue &T,
	      std::vector< std::vector<double> > &cellData,
	      std::vector< std::vector<double> > &wallData,
	      std::vector< std::vector<double> > &vertexData,
	      std::vector< std::vector<double> > &cellDerivs,
	      std::vector< std::vector<double> > &wallDerivs,
	      std::vector< std::vector<double> > &vertexDerivs );
};

#endif
