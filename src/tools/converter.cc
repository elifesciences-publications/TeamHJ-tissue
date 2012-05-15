//
// Filename     : converter.cc
// Description  : Reads and prints states in different formats
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : May 2012
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
  myConfig::registerOption("input_format", 1);
  myConfig::registerOption("output_format", 1);
  //myConfig::registerOption("rk2", 0);
  myConfig::registerOption("help", 0);
  myConfig::registerOption("merry_init", 0);
  myConfig::registerOption("verbose", 1);
  
  //myConfig::initConfig(argc, argv, configFile);
  
  if (myConfig::getBooleanValue("help")) {
    std::cerr << std::endl 
	      << "Usage: " << argv[0] << " initFile > outputFile" << std::endl;
    std::cerr << "Additional flags are:" << std::endl;
    std::cerr << "-input_format format - Sets format in input file." << std::endl
	      << "Available input formats are tissue (default), "
	      << "organism (organism file assuming spheres and only position+radii), " 
	      << "voronoi (voronoi format from qhull output), "
	      << "MGXTriVtu (MGX exported mesh in vtu format, before making cells), " 
	      << "MGXCellVtu (MGX exported mesh in vtu format, after making cells [TO COME]), " 
	      << "merryProj (Montpellier (openAlea) format [now obselete?])." 
	      << std::endl;
    std::cerr << "-output_format format - Sets format for output of"
	      << " state in specified file format. " << std::endl 
	      << "Available output formats are tissue (default), "
	      << "triTissue (tissue with central triangulation), "
	      << "fem (Pawel's FEM simulation format), "
	      << "organism (organism init file including neighborhood), " 
	      << "vtu1 (vtk format with single wall compartment variables), " 
	      << "vtu2 (vtk format with two wall compartment variables)."
	      << std::endl;
    std::cerr << "-verbose flag - Set flag for verbose (flag=1, default) more verbose (2) or "
	      << "silent (0) output mode to stderr." << std::endl; 
    std::cerr << "-help - Shows this message." << std::endl;
    exit(EXIT_FAILURE);
  } else if (myConfig::argc() != 3 ) {
    std::cerr << "Type '" << argv[0] << " -help' for usage." << std::endl;
    exit(EXIT_FAILURE);
  }
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
  
  // Create the tissue and read init file
  Tissue T;
  std::string initFile = myConfig::argv(1);
  std::string inputFormat;
  inputFormat = myConfig::getValue("input_format",0);
  if (inputFormat.empty() || inputFormat.compare("tissue")==0) {
    if (verboseFlag) {
      std::cerr << "Reading init from file " << initFile << " assuming tissue format." << std::endl;
    }
    T.readInit(initFile.c_str(),verboseFlag);
  }
  else if (inputFormat.compare("organism")==0) {
    if (verboseFlag) {
      std::cerr << "Reading init from file " << initFile << " assuming organism (sphere) format." 
		<< std::endl;
    }
    T.readSphereInit(initFile.c_str(),verboseFlag);
  }
  else if (inputFormat.compare("voronoi")==0) {
    if (verboseFlag) {
      std::cerr << "Reading init from file " << initFile << " assuming voronoi (qhull) format." 
		<< std::endl;
    }
    T.readVoronoiInit(initFile.c_str(),verboseFlag);
  }
  else if (inputFormat.compare("MGXTriVtu")==0) {
    if (verboseFlag) {
      std::cerr << "Reading init from file " << initFile << " assuming MGXVtuTri format." 
		<< std::endl;
    }
    T.readMGXTriVtuInit(initFile.c_str(),verboseFlag);
  }
  else if (inputFormat.compare("MGXCellVtu")==0) {
    if (verboseFlag) {
      std::cerr << "Reading init from file " << initFile << " assuming MGXCellVtu format." 
		<< std::endl;
    }
    std::cerr << inputFormat << " not yet implemented for reading! Contact HJ." << std::endl;
    exit(EXIT_FAILURE);
    //T.readMGXCellVtuInit(initFile.c_str(),verboseFlag);
  }
  else if (inputFormat.compare("merryProj")==0) {
    if (verboseFlag) {
      std::cerr << "Reading init from file " << initFile << " assuming merryProj format." 
		<< std::endl;
    }
    T.readMerryInit(initFile.c_str(),verboseFlag);
  }
  else {
    std::cerr << "Input format " << inputFormat << " not recognized. Use '-help' for allowed formats."
	      << std::endl;
    exit(EXIT_FAILURE);
  }
    
  // Print in specified output format
  std::string outputFormat;
  outputFormat = myConfig::getValue("output_format",0);
  if (outputFormat.empty() || outputFormat.compare("tissue")==0) {
    if (verboseFlag) {
      std::cerr << "Printing output using tissue format." << std::endl;
    }
    T.printInit(std::cout);
  }
  else if (inputFormat.compare("triTissue")==0) {
    if (verboseFlag) {
      std::cerr << "Printing output using triangulated tissue format." 
		<< std::endl;
    }
    T.printInitTri(std::cout);
  }
  else if (inputFormat.compare("fem")==0) {
    if (verboseFlag) {
      std::cerr << "Printing output using fem format." << std::endl;
    }
    T.printInitFem(std::cout);
  }
  else {
    std::cerr << "Warning: main() - Format " << outputFormat << " not recognized. "
	      << "No outputwritten." << std::endl;
      }
}
