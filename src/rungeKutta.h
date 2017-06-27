#ifndef RUNGEKUTTA_H
#define RUNGEKUTTA_H

#include "baseSolver.h"

///
/// @brief A fifth order Runge-Kutta solver for ODEs
///
/// This method follows the Numerical Recipes implementation with added cell division functionality.
/// Between each RK step, updates can be made and cell divisions can be implemented. The starting step
/// is also used as a maximal step.
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
	/// @verbatim 
	/// RK5Adaptive
	/// T_start T_end 
	/// printFlag printNum 
	/// h1 eps 
	/// @verbatim
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
		  DataMatrix &yScalC,
		  DataMatrix &yScalW,
		  DataMatrix &yScalV,
		  DataMatrix &yTempC,
		  DataMatrix &yTempW,
		  DataMatrix &yTempV,
		  DataMatrix &yErrC,
		  DataMatrix &yErrW,
		  DataMatrix &yErrV,
		  DataMatrix &ak2C,
		  DataMatrix &ak2W,
		  DataMatrix &ak2V,
		  DataMatrix &ak3C,
		  DataMatrix &ak3W,
		  DataMatrix &ak3V,
		  DataMatrix &ak4C,
		  DataMatrix &ak4W,
		  DataMatrix &ak4V,
		  DataMatrix &ak5C,
		  DataMatrix &ak5W,
		  DataMatrix &ak5V,
		  DataMatrix &ak6C,
		  DataMatrix &ak6W,
		  DataMatrix &ak6V,
		  DataMatrix &yTempRkckC,
		  DataMatrix &yTempRkckW,
		  DataMatrix &yTempRkckV);
	
	void rkck(double h,
		  DataMatrix &yOutC,
		  DataMatrix &yOutW,
		  DataMatrix &yOutV,
		  DataMatrix &yErrC,
		  DataMatrix &yErrW,
		  DataMatrix &yErrV,
		  DataMatrix &ak2C,
		  DataMatrix &ak2W,
		  DataMatrix &ak2V,
		  DataMatrix &ak3C,
		  DataMatrix &ak3W,
		  DataMatrix &ak3V,
		  DataMatrix &ak4C,
		  DataMatrix &ak4W,
		  DataMatrix &ak4V,
		  DataMatrix &ak5C,
		  DataMatrix &ak5W,
		  DataMatrix &ak5V,
		  DataMatrix &ak6C,
		  DataMatrix &ak6W,
		  DataMatrix &ak6V,
		  DataMatrix &yTempRkckC, 
		  DataMatrix &yTempRkckW, 
		  DataMatrix &yTempRkckV ); 
	
	double maxDerivative();
};

///
/// @brief This class solves the ODE using a fourth order Runge-Kutta solver
///
class RK4 : public BaseSolver {

 private:

  double h_;

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
  /// @verbatim 
  /// RK4 
  /// T_start T_end 
  /// printFlag printNum 
  /// h 
  /// @verbatim
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
  void rk4(DataMatrix &ytCell,
	   DataMatrix &ytWall,
	   DataMatrix &ytVertex,
	   DataMatrix &dytCell,
	   DataMatrix &dytWall,
	   DataMatrix &dytVertex,
	   DataMatrix &dymCell,
	   DataMatrix &dymWall,
	   DataMatrix &dymVertex );
  
  double maxDerivative();
};

#endif /* RUNGEKUTTA_H */
