#ifndef RUNGEKUTTA_H
#define RUNGEKUTTA_H

#include "baseSolver.h"

///
/// @brief A fifth order Runge-Kutta solver for ODEs
///
class RK5Adaptive : public BaseSolver {

private:

	double eps_;
	double h1_;

public:
	///
	/// @brief Main constructor
	///
	RK5Adaptive(Tissue *T,std::ifstream &IN);
	
	///
	/// @brief Reads the parameters used by the RK5Adaptive algorithm
	///
	/// This function is responsible for reading parameters used by the fifth
	/// order Runge-Kutta algorithm. The parameter file sent to the simulator
	/// binary looks like:
	///
	/// <pre> 
	/// RK5Adaptive
	/// T_start T_end 
	/// printFlag printNum 
	/// h1 eps 
	/// </pre> 
	///
	/// where RK5Adaptive is the identity string used by BaseSolver::getSolver
	/// to identify that the RK5Adaptive algorithm should be used. T_start
	/// (T_end) is the start (end) time for the simulation, printFlag is an
	/// integer which sets the output format (read by BaseSolver::print()),
	/// printNum is the number of equally spread time points to be printed. h1
	/// is the maximal (and initial) step size for each RK5 step, and eps sets
	/// the error threshold (should be <<1.0).
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
	/// @brief Fifth order Runge-Kutta adaptive stepper.
	///
	void rkqs(double hTry, double &hDid, double &hNext,
						std::vector< std::vector<double> > &yScalC,
						std::vector< std::vector<double> > &yScalW,
						std::vector< std::vector<double> > &yScalV,
						std::vector< std::vector<double> > &yTempC,
						std::vector< std::vector<double> > &yTempW,
						std::vector< std::vector<double> > &yTempV,
						std::vector< std::vector<double> > &yErrC,
						std::vector< std::vector<double> > &yErrW,
						std::vector< std::vector<double> > &yErrV,
						std::vector< std::vector<double> > &ak2C,
						std::vector< std::vector<double> > &ak2W,
						std::vector< std::vector<double> > &ak2V,
						std::vector< std::vector<double> > &ak3C,
						std::vector< std::vector<double> > &ak3W,
						std::vector< std::vector<double> > &ak3V,
						std::vector< std::vector<double> > &ak4C,
						std::vector< std::vector<double> > &ak4W,
						std::vector< std::vector<double> > &ak4V,
						std::vector< std::vector<double> > &ak5C,
						std::vector< std::vector<double> > &ak5W,
						std::vector< std::vector<double> > &ak5V,
						std::vector< std::vector<double> > &ak6C,
						std::vector< std::vector<double> > &ak6W,
						std::vector< std::vector<double> > &ak6V,
						std::vector< std::vector<double> > &yTempRkckC,
						std::vector< std::vector<double> > &yTempRkckW,
						std::vector< std::vector<double> > &yTempRkckV);

	void rkck(double h,
						std::vector< std::vector<double> > &yOutC,
						std::vector< std::vector<double> > &yOutW,
						std::vector< std::vector<double> > &yOutV,
						std::vector< std::vector<double> > &yErrC,
						std::vector< std::vector<double> > &yErrW,
						std::vector< std::vector<double> > &yErrV,
						std::vector< std::vector<double> > &ak2C,
						std::vector< std::vector<double> > &ak2W,
						std::vector< std::vector<double> > &ak2V,
						std::vector< std::vector<double> > &ak3C,
						std::vector< std::vector<double> > &ak3W,
						std::vector< std::vector<double> > &ak3V,
						std::vector< std::vector<double> > &ak4C,
						std::vector< std::vector<double> > &ak4W,
						std::vector< std::vector<double> > &ak4V,
						std::vector< std::vector<double> > &ak5C,
						std::vector< std::vector<double> > &ak5W,
						std::vector< std::vector<double> > &ak5V,
						std::vector< std::vector<double> > &ak6C,
						std::vector< std::vector<double> > &ak6W,
						std::vector< std::vector<double> > &ak6V,
						std::vector< std::vector<double> > &yTempRkckC, 
						std::vector< std::vector<double> > &yTempRkckW, 
						std::vector< std::vector<double> > &yTempRkckV ); 
	
	double maxDerivative();
};

///
/// @brief This class solves the ODE using a fourth order Runge-Kutta solver
///
class RK4 : public BaseSolver {

private:

	double h_;
	unsigned int numStep_;

public:
	
	///
	/// @brief Main constructor.
	///
	RK4(Tissue *T,std::ifstream &IN);
	
	///
	/// @brief Reads the parameters used by the RK4 algorithm
	///
	/// This function is responsible for reading parameters used by the fourth
	/// order Runge-Kutta algorithm. The parameter file sent to the simulator
	/// binary looks like:
	///
	/// <pre> 
	/// RK4 
	/// T_start T_end 
	/// printFlag printNum 
	/// h 
	/// </pre> 
	///
	/// where RK4 is the identity string used by BaseSolver::getSolver to
	/// identify that the RK4 algorithm should be used. T_start (T_end) is the
	/// start (end) time for the simulation, printFlag is an integer which sets
	/// the output format (read by BaseSolver::print()), printNum is the number
	/// of equally spread time points to be printed, and h is the step size for
	/// each RK4 step.
	///
	/// Comments can be included in the parameter file by starting the line with
	/// an #. Caveat: No check on the validity of the read data is applied.
	///
	/// @see BaseSolver::getSolver()
	/// @see BaseSolver::print()
	///
	void readParameterFile(std::ifstream &IN);
	
	///
	/// @brief Runs a simulation of an tissue model
	///
	void simulate(size_t verbose=0);
	
	///
	/// @brief Fourth order Runge-Kutta stepper
	///
	void rk4(std::vector< std::vector<double> > &ytCell,
					 std::vector< std::vector<double> > &ytWall,
					 std::vector< std::vector<double> > &ytVertex,
					 std::vector< std::vector<double> > &dytCell,
					 std::vector< std::vector<double> > &dytWall,
					 std::vector< std::vector<double> > &dytVertex,
					 std::vector< std::vector<double> > &dymCell,
					 std::vector< std::vector<double> > &dymWall,
					 std::vector< std::vector<double> > &dymVertex );
	
	double maxDerivative();
};

#endif /* RUNGEKUTTA_H */

