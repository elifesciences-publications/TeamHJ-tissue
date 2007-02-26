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

//!A cell-based auxin transport model
/*!A complete pattern generating auxin model based on only cellular
  compartments. The four molecules are updated according to:

	dA_i/dt = p0*M_i p1 - p2*A_i +p3*\Sum_{neigh} (A_n-A_i) + 
	p4*\Sum_{neigh} (P_ni*A_n-P_in*A_i) 

	dP_i/dt = p6 - p7*P_i 

	dX_i/dt = p8*A_i - p9*X_i

	dM_i/dt = p10*\Theta_L1 - p11*M_i

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

#endif
