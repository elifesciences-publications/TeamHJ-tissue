//
// Filename     : findCellPeaks.cc
// Description  : Finds peaks within the cells in the tissue 
// Author(s)    : Henrik Jonsson (henrik at thep.lu.se)
// Created      : February 2007
// Revision     : $Id:$
//

#include<fstream>
#include"tissue.h"
#include"cell.h"
#include"wall.h"
#include"vertex.h"

int main(int argc,char *argv[]) {
  
  //Input handling
  if( argc<3 || argc>3 ) {
    std::cerr << "Usage: " << argv[0] << " initFile col" 
							<< std::endl;
    exit(0);
  }
  
  //Create the tissue and read init file
  Tissue T;
  T.readInit(argv[1],2);
  size_t col = atoi( argv[2] );
  
  std::vector< std::vector<double> > cellData( T.numCell() );
  for( size_t i=0 ; i<T.numCell() ; ++i ) {
    cellData[i].resize( T.cell(i).numVariable() );
    for( size_t j=0 ; j<cellData[i].size() ; ++j )
      cellData[i][j]=T.cell(i).variable(j);
  } 
  std::vector<size_t> cellMax;
  std::vector<size_t> flag; 
  T.findPeaksGradientAscent( cellData, col, cellMax, flag );
	
	std::cerr << cellMax.size() << " cell maxima found in column "
						<< col << std::endl;
	for( size_t i=0; i<cellMax.size(); ++i )
		std::cout << cellMax[i] << std::endl; 
	
}
