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

///
/// @brief Adds a variable to each cell.
///
/// Adds a variable to each cell, where the values are defined by type:
///
/// type=uniform: A constant value, given in p[0]
///
int addCellVariable(Tissue &T,std::vector<double> &p,std::string &type);
///
/// @brief Adds a variable to each wall.
///
/// @see addCellVariable()
///
int addWallVariable(Tissue &T,std::vector<double> &p,std::string &type);
///
/// @brief Sets a variable value in each wall.
///
/// The constant value given in p[0] is set in wall variables with index p[1].
///
int setWallVariable(Tissue &T,std::vector<double> &p,std::string &type);
///
/// @brief Adjust vertex positions such that all walls have a length above a threshold value.
///
/// p[0] sets a threshold, and all vetex positions are adjusted such that all cells are at least 
/// p[0] in length.
///
int minimalWallLength(Tissue &T,std::vector<double> &p);
///
/// @brief Adjust wall resting lengths such that all walls have a length above a threshold value.
///
/// p[0] sets a threshold, and all wall resting lengths are adjusted such that all walls are at 
/// most p[0] in resting length.
///
int maximalWallLength(Tissue &T,std::vector<double> &p);
///
/// @brief Sets all wall lengths to be equal to a factor times the distance between the vertices.
///
int wallLengthFromDistance(Tissue &T,std::vector<double> &p);
///
/// @brief Removes all 2-vertices and adjust the tissue accordingly
///
/// A 2-vertex is only connected to two walls (cells) which could be replaced with a single wall.
///
int removeTwoVertices(Tissue &T,std::vector<double> &p);


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
	
	//p.resize(2);
	//p[0]=0.0;
	//p[1]=0;
	//type = new std::string("uniform");
	//setWallVariable(T,p,type);

	p.resize(1);
	addWallVariable(T,p,type);
	addWallVariable(T,p,type);
	
	p[0] = 0.0;
	addCellVariable(T,p,type);
 	addCellVariable(T,p,type);
 	p[0] = 1.0;
 	addCellVariable(T,p,type);
 	addCellVariable(T,p,type);
 	addCellVariable(T,p,type);
	
	p[0]=1.0;
	minimalWallLength(T,p);
	
	p[0]=1.0;
	wallLengthFromDistance(T,p);

	p[0]=10.0;
	maximalWallLength(T,p);

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
	std::cerr << "Init manipulation done." << std::endl;
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

int setWallVariable(Tissue &T,std::vector<double> &p,std::string &type)
{
	if (type.compare("uniform")==0) {
		assert( p.size()==2 );
		size_t index = size_t(p[1]);
		size_t numW=T.numWall();
		for (size_t i=0; i<numW; ++i)
			T.wall(i).setVariable(index,p[0]);
	}
	else {
		std::cerr << "addWallVariable does not accept type " << type << std::endl;
		exit(-1);
	}	
	return 0;
}

int minimalWallLength(Tissue &T,std::vector<double> &p)
{
	assert( p.size()==1 );
	int count=0;
	size_t numW=T.numWall();
	for (size_t i=0; i<numW; ++i) {
		double distance=0.0;
		Vertex *v1=T.wall(i).vertex1();
		Vertex *v2=T.wall(i).vertex2();
		size_t dimension=v1->numPosition();
		std::vector<double> n(dimension);
		for (size_t d=0; d<dimension; ++d) {
			distance += ((v2->position(d) - v1->position(d)) *
									 (v2->position(d) - v1->position(d)));
			n[d] = v2->position(d) - v1->position(d);
		}
		distance = std::sqrt(distance);
		if (distance<p[0]) {
			double normFactor = 1.0/distance; 
			for (size_t d=0; d<dimension; ++d)
				n[d] *= normFactor;
			double addDistance=0.5*(p[0]-distance);
			std::vector<double> v1New(dimension),v2New(dimension);
			for (size_t d=0; d<dimension; ++d) {
				v1New[d] = v1->position(d)-addDistance*n[d];
				v2New[d] = v2->position(d)+addDistance*n[d];
			}
			v1->setPosition( v1New );
			v2->setPosition( v2New );
			++count;
		}
	}
	std::cerr << "minimalWallLength(T,p) adjusted (lengthened) " << count 
						<< " walls to length " << p[0] << " by moving vertices." << std::endl;	
	return 0;
}

int maximalWallLength(Tissue &T,std::vector<double> &p)
{
	assert( p.size()==1 );
	int count=0;
	size_t numW=T.numWall();
	for (size_t i=0; i<numW; ++i) {
		double distance=0.0;
		Vertex *v1=T.wall(i).vertex1();
		Vertex *v2=T.wall(i).vertex2();
		size_t dimension=v1->numPosition();
		std::vector<double> n(dimension);
		for (size_t d=0; d<dimension; ++d) {
			distance += ((v2->position(d) - v1->position(d)) *
									 (v2->position(d) - v1->position(d)));
			n[d] = v2->position(d) - v1->position(d);
		}
		distance = std::sqrt(distance);
		if (distance>p[0]) {
			T.wall(i).setLength( p[0] );
			++count;
		}
	}
	std::cerr << "maximalWallLength(T,p) adjusted (shortened) " << count 
						<< " walls to resting length " << p[0] << std::endl;
	return 0;
}

int wallLengthFromDistance(Tissue &T,std::vector<double> &p)
{
	assert( p.size()==1 );
	size_t numW=T.numWall();
	for (size_t i=0; i<numW; ++i) {
		double distance=0.0;
		Vertex *v1=T.wall(i).vertex1();
		Vertex *v2=T.wall(i).vertex2();
		for (size_t d=0; d<v1->numPosition(); ++d)
			distance += ((v2->position(d) - v1->position(d)) *
									 (v2->position(d) - v1->position(d)));
		distance = std::sqrt(distance);
		T.wall(i).setLength(p[0]*distance);
	}
	return 0;
}

int removeTwoVertices(Tissue &T,std::vector<double> &p) 
{
	assert( p.size()==0 );
	size_t numV=T.numVertex();
	for (size_t i=0; i<numV; ++i) {
		if (T.vertex(i).numWall()==2)
			T.removeTwoVertex(i);
	}
}
