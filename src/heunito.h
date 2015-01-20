#ifndef HEUNITO_H
#define HEUNITO_H
//
// Filename     : heunito.h
// Description  : Heun numerical solver in the Ito interpretation (Carrillo et al. 2003 PRE)
// Author(s)    : Henrik Jonsson (henrik.jonsson@slcu.cam.ac.uk)
// Created      : January 2015
// Revision     : $Id:$
//

#include "baseSolver.h"

///
/// @brief Heun numerical solver in the Ito interpretation
///
/// This class implements the solver described in (Carrillo et al. 2003 PRE)
///
class HeunIto : public BaseSolver {

 private:
  
  double h_;
  double vol_;
  
 public:
  ///
  /// @brief Main constructor
  ///
  HeunIto(Tissue *T,std::ifstream &IN);
  
  ///
  /// @brief Reads the parameters used by the HeunIto solver
  ///
  /// This function is responsible for reading parameters used by the
  /// HeunIto algorithm. The parameter file sent to the simulator
  /// binary looks like:
  ///
  /// <pre> 
  /// HeunIto
  /// T_start T_end 
  /// printFlag printNum 
  /// h vol
  /// </pre> 
  ///
  /// where HeunIto is the identity string used by BaseSolver::getSolver
  /// to identify that the HeunIto algorithm should be used. T_start
  /// (T_end) is the start (end) time for the simulation, printFlag is an
  /// integer which sets the output format (read by BaseSolver::print()),
  /// printNum is the number of equally spread time points to be printed. h
  /// is the step size for each HeunIto step and vol is the effective volume.
  ///
  /// Comments can be included in the parameter file by starting the line with
  /// an #. Caveat: No check on the validity of the read data is applied.
  ///
  /// @see BaseSolver::getSolver()
  /// @see BaseSolver::print()
  ///
  void readParameterFile(std::ifstream &IN);
	
  ///
  /// @brief Runs a simulation of an organism model
  ///
  void simulate(size_t verbose=0);
	
  ///
  /// @brief An HeunIto step
  ///
  void heunito(DataMatrix &sdydtCell,
	       DataMatrix &sdydtWall,
	       DataMatrix &sdydtVertex,
	       DataMatrix &stCell,
	       DataMatrix &stWall,
	       DataMatrix &stVertex,
	       DataMatrix &y1Cell,
	       DataMatrix &y1Wall,
	       DataMatrix &y1Vertex,
	       DataMatrix &dydt2Cell,
	       DataMatrix &dydt2Wall,
	       DataMatrix &dydt2Vertex);
};

#endif /* HEUNITO_H */

