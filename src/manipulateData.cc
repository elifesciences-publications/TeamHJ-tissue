//
// Filename     : manipulateData.cc
// Description  : Reads a data file, makes some manipulations, and prints a new init file.
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : August 2008
// Revision     : $Id:$
//
#include <iostream>
#include <fstream>
#include<vector>

#include "cell.h"
#include "tissue.h"
#include "myConfig.h"
#include "myFiles.h"
//#include "mySignal.h"
#include "myTimes.h"
#include "vertex.h"
#include "wall.h"

///
/// @brief Extracts a timePoint state and adds it into a Tissue
///
void putTimePointIntoTissue(size_t timePoint);
///
/// @brief Reads a tissue data file (with cells and walls)
///
void readData( std::string dataFile,
							 std::vector< std::vector< std::vector<double> > > &cellData,
							 std::vector< std::vector< std::vector<double> > > &wallData,
							 std::vector< std::vector< std::vector<double> > > &vertexPos,
							 std::vector< std::vector< std::vector<size_t> > > &cellVertex,
							 std::vector< std::vector< std::vector<size_t> > > &wallVertex);

int main(int argc,char *argv[]) {
	
  //Command line handling
	myConfig::registerOption("init_output_format", 1);
	//myConfig::registerOption("rk2", 0);
	myConfig::registerOption("help", 0);
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
							<< "Usage: " << argv[0] << " dataFile " << std::endl
							<< std::endl;
		std::cerr << "Possible additional flags are:" << std::endl;
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
  std::string dataFile = myConfig::argv(1);
	std::vector< std::vector< std::vector<double> > > cellData,wallData,vertexPos;
	std::vector< std::vector< std::vector<size_t> > > cellVertex,wallVertex;
	readData(dataFile,cellData,wallData,vertexPos,cellVertex,wallVertex);


	size_t t=10;
	Tissue T(cellData[t],wallData[t],vertexPos[t],cellVertex[t],wallVertex[t]);
	

	// Print time point in init format
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
 	std::cerr << "Data manipulation done." << std::endl;
}


void readData( std::string dataFile,
							 std::vector< std::vector< std::vector<double> > > &cellData,
							 std::vector< std::vector< std::vector<double> > > &wallData,
							 std::vector< std::vector< std::vector<double> > > &vertexPos,
							 std::vector< std::vector< std::vector<size_t> > > &cellVertex,
							 std::vector< std::vector< std::vector<size_t> > > &wallVertex)
{
  std::ifstream IN;
  IN.open( dataFile.c_str() );
  if( !IN ) {
    std::cerr << "readData() Cannot open file " << dataFile << "\n";
    exit(-1);
  }
  
  size_t numTimeStep,numVertex,d,dimension;
	std::vector<double> timePoint;
  
  IN >> numTimeStep;
  timePoint.resize(numTimeStep);
  cellData.resize(numTimeStep);
  wallData.resize(numTimeStep);
  cellVertex.resize(numTimeStep);
  wallVertex.resize(numTimeStep);
  vertexPos.resize(numTimeStep);
	
  for( size_t t=0 ; t<numTimeStep ; ++t ) {
    //Read vertex positions
		IN >> numVertex;
		IN >> d;
		if (t==0) {
			dimension=d;
		}
		if (d!=dimension) {
			std::cerr << "readData() Dimension at time point " << t << " wrong." 
								<< std::endl;
			exit(-1);
		}
    IN >> timePoint[t];
    vertexPos[t].resize(numVertex);
    for( size_t i=0 ; i<numVertex ; ++i ) {
      vertexPos[t][i].resize( dimension );
      for( size_t j=0 ; j<d ; ++j )
				IN >> vertexPos[t][i][j];
    }
    size_t numCell,numVar;
    //Read cell data
    IN >> numCell;
    IN >> numVar;
    cellVertex[t].resize(numCell);
    cellData[t].resize(numCell);
    for( size_t i=0 ; i<numCell ; ++i ) {
      size_t numCellVertex;
      IN >> numCellVertex;
      cellVertex[t][i].resize(numCellVertex);
      for( size_t k=0 ; k<numCellVertex ; ++k )
				IN >> cellVertex[t][i][k];
      cellData[t][i].resize(numVar);
      for( size_t k=0 ; k<numVar ; ++k )
				IN >> cellData[t][i][k];
    }
    size_t numWall,numVarWall;
    //Read wall data
    IN >> numWall;
    IN >> numVarWall;
    wallVertex[t].resize(numWall);
    wallData[t].resize(numWall);
    for( size_t i=0 ; i<numWall ; ++i ) {
      size_t numWallVertex=2;
      //IN >> numWallVertex;
      wallVertex[t][i].resize(numWallVertex);
      for( size_t k=0 ; k<numWallVertex ; ++k )
				IN >> wallVertex[t][i][k];
      wallData[t][i].resize(numVarWall);
      for( size_t k=0 ; k<numVarWall ; ++k )
				IN >> wallData[t][i][k];
    }		
  }
}

void putTimePointIntoTissue(size_t timePoint) {

}
		
