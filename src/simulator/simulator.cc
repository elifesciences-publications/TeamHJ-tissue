//
// Filename     : simulator.cc
// Description  : Simulates a tissue using a 4th or 5th order Runge-Kutta
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : June 2007
// Revision     : $Id:$
//
#include <fstream>

#include "../baseSolver.h"
#include "../cell.h"
#include "../myConfig.h"
#include "../mySignal.h"
#include "../myTimes.h"
#include "../tissue.h"
#include "../vertex.h"
#include "../wall.h"
#include "../pvd_file.h"

int main(int argc,char *argv[]) {
  
  //Command line handling
  myConfig::registerOption("init_output", 1);
  myConfig::registerOption("init_output_format", 1);
  //myConfig::registerOption("rk2", 0);
  myConfig::registerOption("help", 0);
  myConfig::registerOption("merry_init", 0);
  //myConfig::registerOption("wallOutput", 0);
  myConfig::registerOption("verbose", 1);
  myConfig::registerOption("debug_output", 1);
  
  int verboseFlag=1;
  std::string verboseString;
  verboseString = myConfig::getValue("verbose", 0);
  if( !verboseString.empty() ) {
    verboseFlag = atoi( verboseString.c_str() );
    if( verboseFlag != 0 || verboseFlag !=1 || verboseFlag !=2) {
      verboseFlag=0;
      std::cerr << "Flag given to -verbose not recognized (0, 1, 2 allowed)."
		<< " Setting it to zero (silent)." << std::endl;
    }
  }
	
  // Get current time (at start of program)
  myTimes::getTime();
	
  std::string configFile(getenv("HOME"));
  configFile.append("/.tissue");
  myConfig::initConfig(argc, argv, configFile);
  
  if (myConfig::getBooleanValue("help")) {
    std::cerr << std::endl 
	      << "Usage: " << argv[0] << " modelFile initFile "
	      << "simulatorParaFile." << std::endl 
	      << std::endl;
    std::cerr << "Possible additional flags are:" << std::endl;
    std::cerr << "-merry_init - Init file format is set to the "
	      << "output generated from merryproj." << std::endl;
    std::cerr << "-init_output file - Set filename for output of"
	      << " final state in init file format." << std::endl;
    std::cerr << "-init_output_format format - Sets format for output of"
	      << " final state in specified init file format." << std::endl
	      << "Available formats are tissue (default), fem, and triTissue." << std::endl;
    std::cerr << "-verbose flag - Set flag for verbose (flag=1) or "
	      << "silent (0) output mode to stderr." << std::endl; 
    std::cerr << "-debug_output file - Saves the last ten variable"
	      << " states before exiting." << std::endl;
    std::cerr << "-help - Shows this message." << std::endl;
    exit(EXIT_FAILURE);
  } else if (myConfig::argc() != 4 ) {
    std::cerr << "Type '" << argv[0] << " -help' for usage." << std::endl;
    exit(EXIT_FAILURE);
  }
  
  // Create the tissue and read init and model files
  std::string modelFile = myConfig::argv(1);
  std::string initFile = myConfig::argv(2);
  std::string simPara = myConfig::argv(3);
  
  Tissue T;
  if (verboseFlag)
    std::cerr << "Reading model file " << modelFile << std::endl;
  T.readModel(modelFile.c_str(),verboseFlag);
  if (verboseFlag)
    std::cerr << "Reading init file " << initFile << std::endl;	
  if (!myConfig::getBooleanValue("merry_init")) 
    T.readInit(initFile.c_str(),verboseFlag);
  else {
    std::cerr << "Using merryproj init file format" << std::endl;
    T.readMerryInit(initFile.c_str(),verboseFlag);
  }
  
  // Create solver and initiate values
  if (verboseFlag)
    std::cerr << "Generating solver from file " << simPara << std::endl;
  BaseSolver *S = BaseSolver::getSolver(&T, simPara);
  
  // Add solver to signal handler.
  mySignal::addSolver(S);
  
  // Simulate with updates of the neighborhood
  if (verboseFlag)
    std::cerr << "Initiating solver from tissue." << std::endl; 
  S->getInit();
  std::cerr << "Start simulation." << std::endl;
  S->simulate();
  
  // Print init in specified format if applicable
  std::string fileName;
  fileName = myConfig::getValue("init_output", 0);
  if (!fileName.empty()) {
    std::ofstream OUT(fileName.c_str());
    if (!OUT) {
      std::cerr << "Warning: main() -"
		<< "Cannot open file for init output ("
		<< fileName << ")" << std::endl;
    } 
    else {
      std::cerr << "Setting tissue variables from simulator data." << std::endl;
      S->setTissueVariables();
      std::string initFormat;
      initFormat = myConfig::getValue("init_output_format",0);
      if (initFormat.empty() || initFormat.compare("tissue")==0) {
	std::cerr << "Printing init in file " << fileName << " using tissue format." << std::endl;
	S->printInit(OUT);
	OUT.close();
      }
      else if (initFormat.compare("fem")==0) {
	std::cerr << "Printing init in file " << fileName << " using fem format." << std::endl;
	S->printInitFem(OUT);
	OUT.close();
      }
      else if (initFormat.compare("triTissue")==0) {
	std::cerr << "Printing init in file " << fileName << " using triangulated tissue format." << std::endl;
	S->printInitTri(OUT);
	OUT.close();
      }
      else {
	std::cerr << "Warning: main() - Format " << initFormat << " not recognized. "
		  << "No init file written." << std::endl;
      }
    }
  }
}
