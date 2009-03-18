/**
 * Filename     : createInitFromVoronoiFile.cc
 * Description  : Converts a voronoi output to a tissue init file
 * Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
 * Created      : April 2006
 * Revision     : $Id:$
 */

#include<assert.h>
#include<fstream>
#include"tissue.h"
#include"cell.h"
#include"wall.h"
#include"vertex.h"

int main(int argc,char *argv[]) {
	
  //Input handling
  if( argc<2 || argc>2 ) {
    std::cerr << "Usage: " << argv[0] << " inFile" 
							<< std::endl;
    exit(0);
  }
  //Read the voronoi output data
  std::ifstream IN(argv[1]);
  if( !IN ) {
    std::cerr << "Cannot open file " << argv[1] << std::endl; exit(-1);}
  size_t dimension, numVertex,numCell;
	int tmpI=0;
	double tmpD=0.0;
	//Read header information
	IN >> dimension;
	assert(dimension==2 || dimension==3);
  IN >> numVertex;
	--numVertex;
  IN >> numCell;
  IN >> tmpI;
	for( size_t j=0 ; j<dimension ; ++j )
		IN >> tmpD;
	
	//Read vertexPositions
	std::vector< std::vector<double> > vertexPos(numVertex);
	for( size_t i=0 ; i<vertexPos.size() ; ++i ) {
		vertexPos[i].resize(dimension);
		for( size_t j=0 ; j<dimension ; ++j )
			IN >> vertexPos[i][j];
	}
	
	//Read cellVertex connections
	std::vector< std::vector<size_t> > cellVertex(numCell);
	size_t numCellVertex=0;
	for( size_t i=0 ; i<cellVertex.size() ; ++i ) {
		IN >> numCellVertex;
		cellVertex[i].resize(numCellVertex);
		for( size_t j=0 ; j<numCellVertex ; ++j )
			IN >> cellVertex[i][j];
	}
	
	
  //Create the tissue from the data and print the init
	//////////////////////////////////////////////////////////////////////
  int verbose=1;
  Tissue T;
  T.createTissueFromVoronoi(vertexPos,cellVertex,verbose);  
  T.printInit(std::cout);
}
