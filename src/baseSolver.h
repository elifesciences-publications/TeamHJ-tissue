//
// Filename     : baseSolver.h
// Description  : Base class for solvers 
// Author(s)    : Patrik Sahlin (sahlin@thep.lu.se)
//              : Henrik Jonsson (henrik@thep.lu.se)
// Created      : June 2007
// Revision     : $Id:$
//
#ifndef SOLVER_H
#define SOLVER_H

#include "tissue.h"

///
/// @brief A factory class for classes describing different numerical solvers
/// for the ordinary differential equations
///
/// Detailed description to come...
///
class BaseSolver {

protected:
	Tissue *T_;
	std::vector< std::vector<double> > cellData_;
	std::vector< std::vector<double> > wallData_;
	std::vector< std::vector<double> > vertexData_;
	std::vector< std::vector<double> > cellDerivs_; 
	std::vector< std::vector<double> > wallDerivs_; 
	std::vector< std::vector<double> > vertexDerivs_; 
	std::vector< std::vector<std::vector<double> > > cellDataCopy_;
	double t_;
	double startTime_;
	double endTime_;
	int printFlag_;
	int numPrint_;
	unsigned int numOk_, numBad_;
	bool debugFlag_;
	//size_t numSimulation_;
	
public:
	BaseSolver();
	BaseSolver(Tissue *T,std::ifstream &IN);
	virtual ~BaseSolver();
	
	///
	/// @brief This function implements the factory method for initiating a
	/// numerical solver.
	///
	/// This function reads a solver id in the provided file and initiates a
	/// numerical solver of the id type, which is returned. The Tissue (model)
	/// pointer is used to introduce the model to the solver. The initiation is
	/// continued by the individual solver classes since they use different
	/// parameters (see links below). The parameter file sent to the simulator
	/// binary looks like:
	///
	/// @verbatim 
	/// solverId 
	/// parameter_i ...
	/// ...
	/// @endverbatim
	///
	/// where the solverId is the name of the numerical method (class) used. The
	/// parameters used can be found in the links below which lists the
	/// currently available methods/classes.
	///
	/// @see RK5Adaptive::readParameterFile()
	/// @see RK4::readParameterFile()
	///
	static BaseSolver *getSolver(Tissue *T, const std::string &file);
	
	size_t debugCount() const;
	void getInit();
	
	///
	/// @brief General printing function
	///
	/// Caveat: Not yet general, but will be...
	///
	void print(std::ostream &os=std::cout);
	void printInit(std::ostream &os) const;
	void printDebug(std::ostream &os) const;
	virtual void readParameterFile(std::ifstream &IN);
	virtual void simulate(void);
	
	inline double startTime() const;
	inline double endTime() const;
	bool debugFlag() const;
	inline void readInit(const std::string &initFile);
	inline Tissue *getTissue();
	inline void setTissue(Tissue *T);
};

inline double BaseSolver::startTime() const
{
	return startTime_;
}

inline double BaseSolver::endTime() const 
{
	return endTime_;
}

inline bool BaseSolver::debugFlag() const
{
	return debugFlag_;
}

inline void BaseSolver::readInit(const std::string &initFile)
{
	T_->readInit(initFile);
}

inline Tissue *BaseSolver::getTissue()
{
	return T_;
}

inline void BaseSolver::setTissue(Tissue *T)
{
	T_= T;
}

#endif /* SOLVER_H */

