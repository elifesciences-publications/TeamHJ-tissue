//
// Filename     : simulator.cc
// Description  : Simulates a tissue using a 4th or 5th order Runge-Kutta
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : June 2007
// Revision     : $Id:$
//
#include <fstream>

//#include "baseSolver.h"
#include "cell.h"
#include "tissue.h"
#include "myConfig.h"
#include "mySignal.h"
#include "myTimes.h"
#include "vertex.h"
#include "wall.h"

int addCellVariable(Tissue &T,std::vector<double> &p,std::string &type);
int addWallVariable(Tissue &T,std::vector<double> &p,std::string &type);

int main(int argc,char *argv[]) {
	
  //Command line handling
	myConfig::registerOption("init_output_format", 1);
	//myConfig::registerOption("rk2", 0);
	myConfig::registerOption("help", 0);
	myConfig::registerOption("merry_init", 0);
	myConfig::registerOption("verbose", 1);
	
	int verboseFlag=1;
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
							<< "Usage: " << argv[0] << " initFile " << std::endl
							<< std::endl;
		std::cerr << "Possible additional flags are:" << std::endl;
		std::cerr << "-merry_init - Init file format is set to the "
							<< "output generated from merryproj." << std::endl;
		std::cerr << "-init_output_format format - Sets format for output of"
							<< " final state in specified init file format." << std::endl
							<< "Available formats are tissue (default), and fem." << std::endl;
		
		std::cerr << "-verbose flag - Set flag for verbose (flag=1) or "
							<< "silent (0) output mode to stderr." << std::endl; 
		std::cerr << "-help - Shows this message." << std::endl;
    exit(EXIT_FAILURE);
	} else if (myConfig::argc() != 2 ) {
		std::cerr << "Type '" << argv[0] << " -help' for usage." << std::endl;
		exit(EXIT_FAILURE);
  }
	
  // Create the tissue and read init and model files
  std::string initFile = myConfig::argv(1);
	
  Tissue T;
	if (verboseFlag)
		std::cerr << "Reading init file " << initFile << std::endl;	
	if (!myConfig::getBooleanValue("merry_init")) 
		T.readInit(initFile.c_str(),verboseFlag);
	else {
		std::cerr << "Using merryproj init file format" << std::endl;
		T.readMerryInit(initFile.c_str(),verboseFlag);
	}
	
	//
	// Do the manipulations
	//
	std::vector<double> p;
	std::string type("uniform");
	
	p.resize(1);
	p[0]=0.0;
	//type = new std::string("uniform");
	addWallVariable(T,p,type);
	addWallVariable(T,p,type);
	
	p[0] = 1.0;
	addCellVariable(T,p,type);
	addCellVariable(T,p,type);
	
	


	//
  // Print init in specified format
	//
	std::string initFormat;
	initFormat = myConfig::getValue("init_output_format",0);
	if (initFormat.empty() || initFormat.compare("tissue")==0) {
		std::cerr << "Printing init to standard out using tissue format." << std::endl;
		T.printInit(std::cout);
	}
	else if (initFormat.compare("fem")==0) {
		std::cerr << "Printing init to standard out using fem format." << std::endl;
		std::cerr << "NOT YET!" << std::endl;
		//T.printInitFem(std::cout);
	}
	else {
		std::cerr << "Warning: main() - Format " << initFormat << " not recognized. "
							<< "No init file written." << std::endl;
	}
}

int addCellVariable(Tissue &T,std::vector<double> &p,std::string &type)
{
	if (type.compare("uniform")==0) {
		assert( p.size()==1 );
		size_t numC=T.numCell();
		for (size_t i=0; i<numC; ++i)
			T.cell(i).addVariable(p[0]);
	}
	else {
		std::cerr << "addCellVariable does not accept type " << type << std::endl;
		exit(-1);
	}
	return 0;
}

int addWallVariable(Tissue &T,std::vector<double> &p,std::string &type)
{
	if (type.compare("uniform")==0) {
		assert( p.size()==1 );
		size_t numW=T.numWall();
		for (size_t i=0; i<numW; ++i)
			T.wall(i).addVariable(p[0]);
	}
	else {
		std::cerr << "addWallVariable does not accept type " << type << std::endl;
		exit(-1);
	}	
	return 0;
}
