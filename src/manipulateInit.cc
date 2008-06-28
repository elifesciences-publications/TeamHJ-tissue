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
///
/// @brief Scales space (vertex positions and wall resting lengths) such that
/// the maximal area is about p[0]
///
/// Scales vertex positions and wall resting lengths with a factor
/// sqrt(p[0])/sqrt(A_{max}) where A_{max} is the maximal area of a cell.
///
int scaleSpaceToMaxArea(Tissue &T,std::vector<double> &p);
///
/// @brief Scales space (vertex positions and wall resting lengths) such that
/// the maximal wall length is p[0]
///
/// Scales vertex positions and wall resting lengths with a factor
/// p[0]/l_{max} where l_{max} is the maximal length of a wall (defined from
/// its two vertex positions).
///
int scaleSpaceToMaxWallLength(Tissue &T,std::vector<double> &p);
///
/// @brief Flips a variable value around a given value
///
/// Takes a variable and flips the value for all cells/walls/vertices around a
/// given value.
///
/// p[0] is a cell (0) wall (1) vertex (2) flag, p[1] is the variable index,
/// and p[2] is the center value.
///
int flipVariable(Tissue &T,std::vector<double> &p);
///
/// @brief Translates a variable value with given value
///
/// Takes a variable and translates the value for all cells/walls/vertices
/// such that the maximal (minimal) value is a given value.
///
/// p[0] is a cell (0) wall (1) vertex (2) flag, p[1] is the variable index,
/// p[2] is the border value, and p[3] is the max (0) or min (1) flag.
///
int translateVariableToBorder(Tissue &T,std::vector<double> &p);

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

	// For synthetic templates ///////////////////////
	//p.resize(0);
	//removeTwoVertices(T,p);
	//////////////////////////////////////////////////

	p.resize(1);
	// Wall variables 0 0
	p[0]=0.0;
	addWallVariable(T,p,type);
	addWallVariable(T,p,type);
	// Direction 1 (already in .original files)
 	//p[0] = 1.0;
 	//addCellVariable(T,p,type);
	//p[0] = 0.0;
	//addCellVariable(T,p,type);
 	//addCellVariable(T,p,type);
 	//p[0] = 1.0;
 	//addCellVariable(T,p,type);
	// Direction 2
	p[0] = 1.0;
 	addCellVariable(T,p,type);
	p[0] = 0.0;
	addCellVariable(T,p,type);
 	addCellVariable(T,p,type);
 	p[0] = 1.0;
 	addCellVariable(T,p,type);
	// Additional variables 1 1 0 1 1
 	addCellVariable(T,p,type);
 	addCellVariable(T,p,type);
	p[0] = 0.0;
 	addCellVariable(T,p,type);
	p[0] = 1.0;
 	addCellVariable(T,p,type);
 	addCellVariable(T,p,type);
	
	// For experimental template /////////////////////
	p.resize(3);
	p[0] = 2;// vertex
	p[1] = 2;// z
	p[2] = 0.0;// flip around 0
	flipVariable(T,p);
	p.resize(4);
	p[0] = 2;// vertex
	p[1] = 2;// z
	p[2] = 0.0;// move to 0
	p[3] = 1;// move min
	translateVariableToBorder(T,p);
	//////////////////////////////////////////////////

	p.resize(1);
	p[0]=1.0;
	scaleSpaceToMaxArea(T,p);

	p[0]=0.1;
	minimalWallLength(T,p);
	
	p.resize(1);
	p[0]=1.0;
	wallLengthFromDistance(T,p);

	//p[0]=10.0;
	//maximalWallLength(T,p);

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
		std::cerr << "addCellVariable() Uniform cell variable with value " << p[0] << " added."
							<< std::endl;
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
		std::cerr << "addWallVariable() Uniform wall variable with value " << p[0] << " added."
							<< std::endl;
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
	int count=0;
	for (size_t i=0; i<T.numVertex(); ++i) {
		if (T.vertex(i).numWall()==2) {
			T.removeTwoVertex(i--);
			++count;
		}
	}
	std::cerr << "removeTwoVertices() " << count << " two-vertices removed." << std::endl;
	return 0;
}

int scaleSpaceToMaxArea(Tissue &T,std::vector<double> &p)
{
	assert( p.size()==1 );
	// Find maximal area
	double maxA=0.0;
	size_t numC = T.numCell();
	for (size_t i=0; i<numC; ++i) {
		double tmpA = T.cell(i).calculateVolume();
		if (tmpA>maxA)
			maxA=tmpA;
	}
	assert(maxA>0.0);
	// Rescale space with a factor sqrt(p)/sqrt(maxA)
	// space is in vector positions and wall lengths
	double scaleFactor = std::sqrt(p[0]/maxA);
	size_t numV = T.numVertex();
	size_t numW = T.numWall();
	size_t dimension = T.vertex(0).numPosition();
	std::vector<double> tmpPosition(dimension);
	for (size_t i=0; i<numV; ++i) {
		for (size_t d=0; d<dimension; ++d)
			tmpPosition[d] = scaleFactor*T.vertex(i).position(d);
		T.vertex(i).setPosition(tmpPosition);
	}
	for (size_t i=0; i<numW; ++i)
		T.wall(i).setLength(scaleFactor*T.wall(i).length());

	std::cerr << "Scaled maximal area from " << maxA << " to ~" << p[0] << std::endl;
	return 0;
}

int scaleSpaceToMaxWallLength(Tissue &T,std::vector<double> &p)
{
	assert( p.size()==1 );
	// Find maximal wall length (from vertex positions)
	double maxL=0.0;
	size_t numW = T.numWall();
	for (size_t i=0; i<numW; ++i) {
		double tmpL = T.wall(i).lengthFromVertexPosition();
		if (tmpL>maxL)
			maxL=tmpL;
	}
	assert(maxL>0.0);
	// Rescale space with a factor p/maxL
	// space is in vector positions and wall lengths
	double scaleFactor = p[0]/maxL;
	size_t numV = T.numVertex();
	size_t dimension = T.vertex(0).numPosition();
	std::vector<double> tmpPosition(dimension);
	for (size_t i=0; i<numV; ++i) {
		for (size_t d=0; d<dimension; ++d)
			tmpPosition[d] = scaleFactor*T.vertex(i).position(d);
		T.vertex(i).setPosition(tmpPosition);
	}
	for (size_t i=0; i<numW; ++i)
		T.wall(i).setLength(scaleFactor*T.wall(i).length());
	
	std::cerr << "Scaled maximal wallDistance from " << maxL << " to ~" << p[0] << std::endl;
	return 0;
}

int flipVariable(Tissue &T,std::vector<double> &p)
{
	assert( p.size()==3 );

	size_t flag = size_t(p[0]);
	size_t index = size_t(p[1]);
	double centerValue = p[2];

	if (flag==0) {//cell variable
		std::cerr << "flipVariable() not implemented for cells yet." << std::endl;
		exit(-1);
	}
	else if (flag==1) {//wall variable
		std::cerr << "flipVariable() not implemented for walls yet." << std::endl;
		exit(-1);
	}
	else if (flag==2) {//vertex variable (position)
		size_t dimension = T.vertex(0).numPosition();
		assert(dimension>index);
		size_t numV = T.numVertex();
		for (size_t i=0; i<numV; ++i) {
			double delta = T.vertex(i).position(index)-centerValue;
			T.vertex(i).setPosition(index,T.vertex(i).position(index) - 2*delta);
		}
		std::cerr << "flipVariable() Position " << index << " flipped around " << centerValue
							<< "." << std::endl;
		return 0;
	}
	else {
		std::cerr << "flipVariable() wrong flag given as p[0]." << std::endl;
		exit(-1);
	}
}

int translateVariableToBorder(Tissue &T,std::vector<double> &p)
{
	assert( p.size()==4 );

	size_t flag = size_t(p[0]);
	size_t index = size_t(p[1]);
	double border = p[2];
	size_t maxMinFlag = size_t(p[3]);

	if (flag==0) {//cell variable
		std::cerr << "translateVariableToBorder() not implemented for cells yet." << std::endl;
		exit(-1);
	}
	else if (flag==1) {//wall variable
		std::cerr << "translateVariableToBorder() not implemented for walls yet." << std::endl;
		exit(-1);
	}
	else if (flag==2) {//vertex variable (position)
		size_t dimension = T.vertex(0).numPosition();
		assert(dimension>index);
		size_t numV = T.numVertex();

		if (maxMinFlag==0) {//translate such that max is at given border
			double max = T.vertex(0).position(index);
			for (size_t i=1; i<numV; ++i) {
				if (T.vertex(i).position(index)>max)
					max=T.vertex(i).position(index);
			}
			double delta = max-border;
			for (size_t i=0; i<numV; ++i) {
				T.vertex(i).setPosition(index,T.vertex(i).position(index) - delta);
			}
		}
		else if (maxMinFlag==1) {//translate such that min is at given border
			double min = T.vertex(0).position(index);
			for (size_t i=1; i<numV; ++i) {
				if (T.vertex(i).position(index)<min)
					min=T.vertex(i).position(index);
			}
			double delta = min-border;
			for (size_t i=0; i<numV; ++i) {
				T.vertex(i).setPosition(index,T.vertex(i).position(index) - delta);
			}
		}
		else {
			std::cerr << "translateVariableToBorder() wrong minmaxflag given as p[3]." << std::endl;
			exit(-1);
		}
		std::cerr << "translateVariableToBorder() Vertex positions in dimension " << index 
							<< " translated such that ";
		if (maxMinFlag==0)
			std::cerr << "max";
		else
			std::cerr << "min";
		std::cerr << " translated to " << border << "." << std::endl;
		return 0;
	}
	else {
		std::cerr << "translateVariableToBorder() wrong flag given as p[0]." << std::endl;
		exit(-1);
	}
}

		
