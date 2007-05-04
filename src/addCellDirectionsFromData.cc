//
// Filename     : addCellDirectionsFromData.cc
// Description  : Reads a tissue init and a file containing directions vectors
//              : and prints a new init with directions saved as cell data
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
	myConfig::registerOption("help", 0);
	myConfig::registerOption("merry_init", 0);
	myConfig::registerOption("init_output", 1);
	
	size_t verbose=1;
	std::string configFile(getenv("HOME"));
	configFile.append("/.tissue");
	myConfig::initConfig(argc, argv, configFile);

	if (myConfig::getBooleanValue("help")) {
		std::cerr << std::endl 
							<< "Usage: " << argv[0] << " initFile directionFile."
							<< std::endl 
							<< std::endl;
		std::cerr << "Possible additional flags are:" << std::endl;
		std::cerr << "-init_output file - Set filename for output of"
							<< " final state in init file format (stdout default)." 
							<< std::endl;
		std::cerr << "-merry_init - Init file format is set to the "
							<< "output generated from merryproj." << std::endl;
    exit(EXIT_FAILURE);
	} 
	else if (myConfig::argc() < 3 || myConfig::argc() > 4) {
		std::cout << "Type '" << argv[0] << " -help' for usage." << std::endl;
		exit(EXIT_FAILURE);
  }
	
  //Create the tissue and read init file
  std::string initFile = myConfig::argv(1);
	std::string directionFile = myConfig::argv(2);
  Tissue T;
  //T.readModel(modelFile.c_str());
	if (!myConfig::getBooleanValue("merry_init")) 
		T.readInit(initFile.c_str(),verbose);
	else {
		std::cerr << "Using merryproj init file format" << std::endl;
		T.readMerryInit(initFile.c_str(),verbose);
	}
	
	//Read directions, find corresponding cell and put data in cell variables
	size_t numCell = T.numCell();
	size_t numDirection;
	size_t dimension = T.numDimension();
	std::vector< std::vector<double> > cellPosition(numCell),direction,directionPos;
	std::vector<size_t> cellFlag(numCell),directionCell;
	for (size_t i=0; i<numCell; ++i) {
		cellPosition[i] = T.cell(i).positionFromVertex();
	}
	std::ifstream IN(directionFile.c_str());
	if (!IN) {
		std::cerr << "Warning: main() -"
							<< "Cannot open file for reading directions ("
							<< directionFile << ")" << std::endl;
		exit(EXIT_FAILURE);
	} 
	else {
		size_t tmpI;
		double tmpD;
		std::vector<double> tmpPos1(dimension),tmpPos2(dimension);
		IN >> numDirection;
		direction.resize(numDirection);
		directionPos.resize(numDirection);
		directionCell.resize(numDirection);
		for (size_t d=0; d<numDirection; ++d) {
			direction[d].resize(dimension);
			directionPos[d].resize(dimension);
			IN >> tmpI;//index
			assert(d==tmpI);
			for (size_t dim=0; dim<dimension; ++dim)//first pos
				IN >> tmpPos1[dim];
			for (size_t dim=dimension; dim<3; ++dim)//additional zeros
				IN >> tmpD;
			for (size_t dim=0; dim<dimension; ++dim)//second pos
				IN >> tmpPos2[dim];
			for (size_t dim=dimension; dim<3; ++dim)//additional zeros
				IN >> tmpD;
			
			double norm=0.0;
			for (size_t dim=0; dim<dimension; ++dim) {
				directionPos[d][dim] = 0.5*(tmpPos1[dim]+tmpPos2[dim]);
				direction[d][dim] = tmpPos2[dim]-tmpPos1[dim];
				norm += direction[d][dim]*direction[d][dim];
			}
			norm = std::sqrt(norm);
			assert(norm>0.0);
			norm = 1.0/norm;
			for (size_t dim=0; dim<dimension; ++dim)
				direction[d][dim] *= norm;
		}
	}
	IN.close();

	//Find cells connected to directions
	for (size_t d=0; d<numDirection; ++d) {
		double distance=0.0;
		for (size_t dim=0; dim<dimension; ++dim)
			distance += (directionPos[d][dim]-cellPosition[0][dim])*
				(directionPos[d][dim]-cellPosition[0][dim]);
		double minDistance=distance;
		size_t cellIndex=0;
		for (size_t i=1; i<numCell; ++i) {
			distance=0.0;
			for (size_t dim=0; dim<dimension; ++dim)
				distance += (directionPos[d][dim]-cellPosition[i][dim])*
					(directionPos[d][dim]-cellPosition[i][dim]);
			if (distance<minDistance) {
				minDistance = distance;
				cellIndex = i;
			}
		}
		cellFlag[cellIndex]++;
		directionCell[d] = cellIndex;
	} 
	for (size_t i=0; i<numCell; ++i)
		assert( cellFlag[i]==0 || cellFlag[i]==1 );
	
	//Add directions to cell variables (initiate with zeros)
	size_t numOldVar=T.cell(0).numVariable();
	for (size_t i=0; i<numCell; ++i) {
		for (size_t dim=0; dim<dimension; ++dim)
			T.cell(i).addVariable(0.0);
		T.cell(i).addVariable(0.0);
	}
	for (size_t d=0; d<numDirection; ++d) {
		size_t tmpCellI=directionCell[d];
		for (size_t dim=0; dim<dimension; ++dim)
			T.cell(tmpCellI).setVariable(dim+numOldVar,direction[d][dim]);
		T.cell(tmpCellI).setVariable(numOldVar+dimension,1.0);
	}
	
	// Print init to stdout or file given
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
	else {
		T.printInit(std::cout);
	}
}
