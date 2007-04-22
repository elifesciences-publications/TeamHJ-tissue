//
// Filename     : simulator.cc
// Description  : Simulates a tissue using a 4th or 2nd order Runge-Kutta
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : April 2007
// Revision     : $Id:$
//
#include <fstream>

#include "cell.h"
#include "tissue.h"
#include "myConfig.h"
#include "vertex.h"
#include "wall.h"

int main(int argc,char *argv[]) {
	
  //Command line handling
	myConfig::registerOption("init_output", 1);
	myConfig::registerOption("rk2", 0);
	myConfig::registerOption("help", 0);
	myConfig::registerOption("merry_init", 0);
	
	size_t verbose=1;
	std::string configFile(getenv("HOME"));
	configFile.append("/.tissue");
	myConfig::initConfig(argc, argv, configFile);

	if (myConfig::getBooleanValue("help")) {
		std::cerr << std::endl 
							<< "Usage: " << argv[0] << " modelFile initFile "
							<< "[simulatorParaFile]." << std::endl 
							<< std::endl;
		std::cerr << "Possible additional flags are:" << std::endl;
		std::cerr << "-init_output file - Set filename for output of"
							<< " final state in init file format." << std::endl;
		std::cerr << "-merry_init - Init file format is set to the "
							<< "output generated from merryproj." << std::endl;
		std::cerr << "-rk2 - Overrides default numerical solver (rk4)"
							<< " with a 2nd order Runge-Kutta solver." << std::endl;
    exit(EXIT_FAILURE);
	} else if (myConfig::argc() < 3 || myConfig::argc() > 4) {
		std::cout << "Type '" << argv[0] << " -help' for usage." << std::endl;
		exit(EXIT_FAILURE);
  }

  double startTime=0.0,endTime=100.0,step=0.01;
  size_t printNum=100;
  if (argc>3) {
		std::string simPara = myConfig::argv(3);
		std::ifstream IN(simPara.c_str());
    if (!IN) {
      std::cerr << "Cannot open file " << argv[2] << std::endl; exit(-1);}
    IN >> startTime;
    IN >> endTime;
    IN >> step;
    IN >> printNum;
    IN.close();
  }
  //Create the tissue and read init and model files
	std::string modelFile = myConfig::argv(1);
  std::string initFile = myConfig::argv(2);
  Tissue T;
  T.readModel(modelFile.c_str());
	if (!myConfig::getBooleanValue("merry_init")) 
		T.readInit(initFile.c_str(),verbose);
	else {
		std::cerr << "Using merryproj init file format" << std::endl;
		T.readMerryInit(initFile.c_str(),verbose);
	}
	if (!myConfig::getBooleanValue("rk2")) 
		T.simulateRk4(startTime,endTime,step,printNum);
	else
		T.simulateRk2(startTime,endTime,step,printNum);
	
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
			T.printInit(OUT);
			OUT.close();
		}
	}
}
