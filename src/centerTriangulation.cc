//
// Filename     : centerTriangulation.cc
// Description  : Classes describing updates for tissues with cells storing a central point (and internal edges)
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : October 2012
// Revision     : $Id:$
//
#include <vector>
#include "baseReaction.h"
#include "centerTriangulation.h"
#include "tissue.h"

#include <ctime>
using namespace std;

namespace CenterTriangulation {

  Initiate::
  Initiate(std::vector<double> &paraValue, 
	   std::vector< std::vector<size_t> > &indValue)
  {
    // Do some checks on the parameters and variable indeces
    if( paraValue.size()!=0 && paraValue.size()!=1 ) {
      std::cerr << "CenterTriangulation::Initiate::Initiate() "
		<< "Uses zero or one parameter. If one provided it is a flag equal to one for overriding "
		<< "a central triangulation stored in tissue (generated from scratch instead)." << std::endl;
      exit(EXIT_FAILURE);
    }
    if( indValue.size()!=1 || indValue[0].size()!=1 ) { 
      std::cerr << "CenterTriangulation::Initiate::Initiate() "
		<< "Start of internal cell variables index given in first level (=cell.numVariable)." << std::endl;
      exit(EXIT_FAILURE);
    }
    
    // Set the variable values
    setId("CenterTriangulation::Initiate");
    setParameter(paraValue);  
    setVariableIndex(indValue);
    
    // Set the parameter identities
    std::vector<std::string> tmp( numParameter() );
    if (paraValue.size())
      tmp[0] = "overRideFlag";
    setParameterId( tmp );    
  }
  
  void Initiate::
  initiate(Tissue &T,
	   DataMatrix &cellData,
	   DataMatrix &wallData,
	   DataMatrix &vertexData,
	   DataMatrix &cellDerivs,
	   DataMatrix &wallDerivs,
	   DataMatrix &vertexDerivs)
  {
    size_t dimension=3; //Only implemented for 3D models
    assert (dimension==vertexData[0].size());
    size_t numVariable = T.cell(0).numVariable();
    assert (numVariable==cellData[0].size());

    // Create the new variables
    if (variableIndex(0,0) != numVariable) {
      std::cerr << "CenterTriangulation::Initiate::initiate() "
		<< "Wrong index given as start index for additional variables."
		<< std::endl;
      exit(EXIT_FAILURE);
    }
    size_t numCell = cellData.size();
    std::vector<double> com(dimension);    
    assert (numCell==T.numCell());
    if (!T.cell(0).numCenterPosition() || (numParameter() && parameter(0)==1) ) {
      //No centers defined in cells in tissue or overridden (create central triangulation from scratch)
      for (size_t cellIndex=0; cellIndex<numCell; ++cellIndex) {	
        size_t numInternalWall = T.cell(cellIndex).numVertex();
	cellData[cellIndex].resize(numVariable+dimension+numInternalWall);
	cellDerivs[cellIndex].resize(numVariable+dimension+numInternalWall);
	com = T.cell(cellIndex).positionFromVertex(vertexData);
	// Set center position to com of the cell
	for (size_t d=0; d<dimension; ++d)
	  cellData[cellIndex][numVariable+d] = com[d];    
	// Set internal wall lengths to the distance btw com and the vertex
	for (size_t k=0; k<numInternalWall; ++k) {
	  Vertex *tmpVertex = T.cell(cellIndex).vertex(k); 
	  size_t vertexIndex = tmpVertex->index();
	  double distance = std::sqrt( (com[0]-vertexData[vertexIndex][0])*
				       (com[0]-vertexData[vertexIndex][0])+
				       (com[1]-vertexData[vertexIndex][1])*
				       (com[1]-vertexData[vertexIndex][1])+
				       (com[2]-vertexData[vertexIndex][2])*
				       (com[2]-vertexData[vertexIndex][2]) );   
	  cellData[cellIndex][numVariable+dimension+k] = distance;
	}

        // sleep(1); //initiating MT dirrections randomly in xy plane
        // srand((unsigned)time(NULL));
        
        // cellData[cellIndex][0]=((double)rand() / (RAND_MAX+1)) ;
        // cellData[cellIndex][1]=((double)rand() / (RAND_MAX+1)) ;
        // double tmp=-std::sqrt(cellData[cellIndex][0]*cellData[cellIndex][0]+cellData[cellIndex][1]*cellData[cellIndex][1]); 
        
        // cellData[cellIndex][0]/=tmp;
        // cellData[cellIndex][1]/=tmp;
        // std::cerr<< "MT vector"<< cellData[cellIndex][0]<<" "<< cellData[cellIndex][1]<< std::endl;
      }
    }
    else { 
      // Copy central position and edge length data from cells in Tissue
      for (size_t cellIndex=0; cellIndex<numCell; ++cellIndex) {
	size_t numInternalWall = T.cell(cellIndex).numVertex();
	cellData[cellIndex].resize(numVariable+dimension+numInternalWall);
	cellDerivs[cellIndex].resize(numVariable+dimension+numInternalWall);
	com = T.cell(cellIndex).centerPosition();
	// Set center position to com of the cell
	for (size_t d=0; d<dimension; ++d)
	  cellData[cellIndex][numVariable+d] = com[d];    
	// Set internal wall lengths to the distance btw com and the vertex
	for (size_t k=0; k<numInternalWall; ++k) {
	  cellData[cellIndex][numVariable+dimension+k] = T.cell(cellIndex).edgeLength(k);
	}
      }
    }
  }
  
  void Initiate::
  derivs(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs )
  {
    // Does nothing
  }
  
  
}
