/**
 * Filename     : simulatorRk2.cc
 * Description  : Simulates a tissue using a 2nd order Runge-Kutta
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
  if( argc<3 || argc>4 ) {
    std::cerr << "Usage: " << argv[0] << " modelFile initFile [paraFile]" 
	      << std::endl;
    exit(0);
  }
  double startTime=0.0,endTime=100.0,step=0.01;
  size_t printNum=100;
  if( argc>3 ) {
    std::ifstream IN(argv[3]);
    if( !IN ) {
      std::cerr << "Cannot open file " << argv[2] << std::endl; exit(-1);}
    IN >> startTime;
    IN >> endTime;
    IN >> step;
    IN >> printNum;
    IN.close();
  }
  //Create the tissue and read init and model files
  Tissue T;
  T.readInit(argv[2],2);
  T.readModel(argv[1]);

  T.simulateRk2(startTime,endTime,step,printNum);
}
