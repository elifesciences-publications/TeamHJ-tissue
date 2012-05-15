/**
 * Filename     : createInitFromSphereFile.cc
 * Description  : Converts a file with sphere cells to a tissue init file
 * Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
 * Created      : April 2006
 * Revision     : $Id:$
 */

#include<fstream>
#include"tissue.h"
#include"cell.h"
#include"wall.h"
#include"vertex.h"

int main(int argc,char *argv[]) {

  //Input handling
  if( argc<2 || argc>4 ) {
    std::cerr << "Usage: " << argv[0] << " inFile [rFac=1] [dimension]" 
	      << std::endl;
    exit(0);
  }
  double rFac=1.0;
  if( argc>2 )
    rFac = std::atof(argv[2]); 
  size_t dimension=0;
  size_t dimensionFlag=0;
  if( argc>3 ) {
    dimension = std::atoi(argv[3]); 
    dimensionFlag=1;
  }
  
  //Read the sphere data file
  std::ifstream IN(argv[1]);
  if( !IN ) {
    std::cerr << "Cannot open file " << argv[1] << std::endl; exit(-1);}
  size_t numCell,numCol;
  IN >> numCell;
  IN >> numCol;
  //If not given asssume only positions and radii in file
  if( !dimensionFlag )
    dimension = numCol-1;
  
  if( dimension<2 || dimension>3 ) {
    std::cerr << "Dimension needs to be 2 or 3." << std::endl;
    exit(0);
  }
  std::vector< std::vector<double> > data(numCell);
  double tmp;
  for( size_t i=0 ; i<numCell ; ++i ) {
    data[i].resize( dimension+1 );
    for( size_t j=0 ; j<data[i].size() ; ++j )
      IN >> data[i][j];
    for( size_t j=data[i].size() ; j<numCol ; ++j )
      IN >> tmp;		
  }
  
  //Create the tissue from the sphere data and print the init
  int verbose=1;
  Tissue T;
  T.createTissueFromSpheres(data,rFac,verbose);  
  T.printInit(std::cout);
}
