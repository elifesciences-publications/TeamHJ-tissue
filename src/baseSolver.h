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
/// Each solver is implemented inheriting the BaseSolver, and main differences will 
/// be the constructor that is typically from reading infomration from a file and
/// simulate() that implements the solver. Note that the different solvers may have
/// different numbers of copies of the data structures cellData, wallDeata, vertexData
/// as needed by the numerical algorithms.
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
  static BaseSolver* getSolver(Tissue *T, const std::string &file);
  
  size_t debugCount() const;
  
  ///
  /// @brief Sets internal variables from values in the tissue.
  ///
  void getInit();
  
  ///
  /// @brief Updates the tissue variables from the current state of the internal variables
  ///
  void setTissueVariables();
  
  ///
  /// @brief General printing function
  ///
  /// This is the main print function for output data during a simulation. It has a couple
  /// standard modes given as a flag in the input file (together with the number of time
  /// points to print):
  /// @verbatim
  /// 0) Standard output for openGL developed plotting of cell and wall variables 
  /// 1) Standard output for openGL developed plotting of cell variables 
  /// 2) Standard output for openGL developed plotting of wall variables
  /// 3) Output that can be used for plotting in gnuplot.
  ///
  /// In addition there are several methods for plotting also membrane data (e.g. PIN1),
  /// as well as specific methods.
  /// @endverbatime 
  /// Caveat: Not yet general, but will be...?
  ///
  void print(std::ostream &os=std::cout);
  ///
  /// @brief Prints standard tissue init
  ///
  /// Prints the current state in init format using the data matrices.
  /// It also checks that sizes in data and tissue are equal.
  ///
  void printInit(std::ostream &os) const;
  /// 
  /// @brief Prints init in Pawels FEM format
  ///
  /// Prints the current state in Pawels FEM init format using vertexData and
  /// tissue connections.
  ///
  void printInitFem(std::ostream &os) const;
  void printDebug(std::ostream &os) const;
  
  virtual void readParameterFile(std::ifstream &IN);
  ///
  /// @brief The main function doing the actual numerical integration from startTime to endTime
  ///
  /// Each solver has this function implemented according to its algorithm.
  ///
  virtual void simulate(size_t verbose=0);
  
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

