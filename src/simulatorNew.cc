//
// Filename     : simulator.cc
// Description  : Simulates a tissue using a 4th or 5th order Runge-Kutta
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : June 2007
// Revision     : $Id:$
//
#include <fstream>

#include "baseSolver.h"
#include "cell.h"
#include "tissue.h"
#include "myConfig.h"
#include "mySignal.h"
#include "myTimes.h"
#include "vertex.h"
#include "wall.h"

int main(int argc,char *argv[]) {
	
  //Command line handling
	myConfig::registerOption("init_output", 1);
	//myConfig::registerOption("rk2", 0);
	myConfig::registerOption("help", 0);
	myConfig::registerOption("merry_init", 0);
	//myConfig::registerOption("wallOutput", 0);
	myConfig::registerOption("verbose", 1);
	myConfig::registerOption("debug_output", 1);

	int verboseFlag=0;
	std::string verboseString;
	verboseString = myConfig::getValue("verbose", 0);
	if( !verboseString.empty() ) {
		verboseFlag = atoi( verboseString.c_str() );
		if( verboseFlag != 0 || verboseFlag !=1 ) {
			verboseFlag=0;
			std::cerr << "Flag given to -verbose not recognized (0, 1 allowed)."
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
		std::cerr << "-init_output file - Set filename for output of"
							<< " final state in init file format." << std::endl;
		std::cerr << "-merry_init - Init file format is set to the "
							<< "output generated from merryproj." << std::endl;
		//std::cerr << "-rk2 - Overrides default numerical solver (rk4)"
		//				<< " with a 2nd order Runge-Kutta solver." << std::endl;
		//std::cerr << "-wallOutput - Overrides default cell output with "
		//				<< "wall output." << std::endl;
		std::cerr << "-verbose flag - Set flag for verbose (flag=1) or "
							<< "silent (0) output mode to stderr." << std::endl; 
		std::cerr << "-debug_output file - Saves the last ten variable"
							<< " states before exiting." << std::endl;
    exit(EXIT_FAILURE);
	} else if (myConfig::argc() != 4 ) {
		std::cout << "Type '" << argv[0] << " -help' for usage." << std::endl;
		exit(EXIT_FAILURE);
  }
	
  // Create the tissue and read init and model files
	std::string modelFile = myConfig::argv(1);
  std::string initFile = myConfig::argv(2);
  std::string simPara = myConfig::argv(3);

  Tissue T;
  T.readModel(modelFile.c_str(),verboseFlag);
	if (!myConfig::getBooleanValue("merry_init")) 
		T.readInit(initFile.c_str(),verboseFlag);
	else {
		std::cerr << "Using merryproj init file format" << std::endl;
		T.readMerryInit(initFile.c_str(),verboseFlag);
	}
	
	// Create solver and initiate values
	BaseSolver *S = BaseSolver::getSolver(&T, simPara);

  // Add solver to signal handler.
  mySignal::addSolver(S);
	
  // Simulate with updates of the neighborhood
  std::cerr << "Start simulation.\n";
  S->getInit();
	S->simulate();
	
  // Print init if applicable
	std::string fileName;
	fileName = myConfig::getValue("init_output", 0);
	if(!fileName.empty()) {
		std::ofstream OUT(fileName.c_str());
		if (!OUT) {
			std::cerr << "Warning: main() -"
								<< "Cannot open file for init output ("
								<< fileName << ")" << std::endl;
		} 
		else {
			S->printInit(OUT);
			OUT.close();
		}
	}
}
