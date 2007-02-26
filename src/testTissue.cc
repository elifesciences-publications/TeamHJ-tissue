/**
 * Filename     : testTissue.cc
 * Description  : A class for testing a tissue
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

  //Read from init file and print
  //////////////////////////////////////////////////////////////////////
  //Tissue T("square_2_2.init",1);
  //T.printVertexAndCell(std::cout);
  //T.printWall(std::cout);

  //Create sphere data and crete tissue from these data
  //////////////////////////////////////////////////////////////////////
  //size_t N=6,dim=2;
  //std::vector< std::vector<double> > y(N);
  //for( size_t i=0 ; i<N ; i++ )
  //y[i].resize(dim+1);
  //double r=0.8;
  //y[0][0]=0.0;y[0][1]=0.0;y[0][2]=r;
  //y[1][0]=0.0;y[1][1]=1.0;y[1][2]=r;
  //y[2][0]=1.0;y[2][1]=0.0;y[2][2]=r;
  //y[3][0]=1.0;y[3][1]=1.0;y[3][2]=r;
  //y[4][0]=2.0;y[4][1]=0.0;y[4][2]=r;
  //y[5][0]=2.0;y[5][1]=1.0;y[5][2]=r;
  //Tissue T;
  //T.createTissueFromSpheres(y);

  //Read sphere data from file and create tissue
  //////////////////////////////////////////////////////////////////////
  //if( argc != 2 ) {
  //std::cerr << "Usage: " << argv[0] << " inFile" << std::endl;
  //exit(0);
  //}
  //std::ifstream IN(argv[1]);
  //if( !IN ) {
  //std::cerr << "Cannot open file " << argv[1] << std::endl; exit(-1);}
  //size_t numCell,numCol,dimension;
  //IN >> numCell;
  //IN >> numCol;
  //dimension = numCol-1;
  //if( dimension<2 || dimension>3 ) {
  //std::cerr << "Dimension needs to be 2 or 3." << std::endl;
  //exit(0);
  //}
  //std::vector< std::vector<double> > data(numCell);
  //for( size_t i=0 ; i<numCell ; ++i ) {
  //data[i].resize( numCol );
  //for( size_t j=0 ; j<numCol ; ++j )
  //  IN >> data[i][j];
  //}
  //Tissue T;
  //T.createTissueFromSpheres(data,1.2,1);  
  //T.printVertexAndCell(std::cout);
  //T.printWall(std::cout);

  //Read and print init file
  //////////////////////////////////////////////////////////////////////
  if( argc != 2 ) {
    std::cerr << "Usage: " << argv[0] << " inFile" << std::endl;
    exit(0);
  }
  //Create the tissue and read init and model files
  Tissue T;
  T.readInit(argv[1],2);
  T.printWall(std::cout);
  //T.printInit();


}
