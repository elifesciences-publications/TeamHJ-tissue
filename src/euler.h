#ifndef EULER_H
#define EULER_H

#include "baseSolver.h"

///
/// @brief An explicit Euler solver
///
/// Simple explicit Euler method for numrical integration of the system.
/// It updates the state (y={cellData,WallData,vertexData}) from the derivatives
/// given by the defined reactions (dydt={cellDerivs,wallDerivs,vertexDerivs})
/// according to:
///
/// @f[ y_{n+1} = y_n + dydt_n h @f]
///
/// where h is the time step provided in the parameter file Euler::readParameterFile().
///
class Euler : public BaseSolver {

private:

	double h_;

public:
	///
	/// @brief Main constructor
	///
	Euler(Tissue *T,std::ifstream &IN);
	
	///
	/// @brief Reads the parameters used by the Euler algorithm
	///
	/// This function is responsible for reading parameters used by the fifth
	/// order Runge-Kutta algorithm. The parameter file sent to the simulator
	/// binary looks like:
	///
	/// <pre> 
	/// Euler
	/// T_start T_end 
	/// printFlag printNum 
	/// h1
	/// </pre> 
	///
	/// where Euler is the identity string used by BaseSolver::getSolver
	/// to identify that the Euler algorithm should be used. T_start
	/// (T_end) is the start (end) time for the simulation, printFlag is an
	/// integer which sets the output format (read by BaseSolver::print()),
	/// printNum is the number of equally spread time points to be printed. h1
	/// is the step size for each Euler step, and eps sets
	///
	/// Comments can be included in the parameter file by starting the line with
	/// an #. Caveat: No check on the validity of the read data is applied.
	///
	/// @see BaseSolver::getSolver()
	/// @see BaseSolver::print()
	///
	void readParameterFile(std::ifstream &IN);
	
	void simulate(size_t verbose=0);
			
	///
	/// @brief A single Euler step
	///
	void eulerStep();
};


#endif //EULER_H

