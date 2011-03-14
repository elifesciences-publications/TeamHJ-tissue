//
// Filename     : mechanicalTRBS.cc
// Description  : Classes describing updates due to mechanical triangular biquadratic springs
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : February 2011
// Revision     : $Id:$
//
#include <utility>
#include <vector>
#include "baseReaction.h"
#include "mechanicalTRBS.h"
#include "tissue.h"
#include <cmath>

VertexFromTRBS::
VertexFromTRBS(std::vector<double> &paraValue, 
	       std::vector< std::vector<size_t> > 
	       &indValue ) 
{  
  // Do some checks on the parameters and variable indeces
  if( paraValue.size()!=2 ) {
    std::cerr << "VertexFromTRBS::"
	      << "VertexFromTRBS() "
	      << "Uses two parameters young modulus and poisson coefficient.\n";
    exit(0);
  }
  if( indValue.size()!=1 || indValue[0].size()!=1 ) { 
    std::cerr << "VertexFromTRBS::"
	      << "VertexFromTRBS() "
	      << "Only wall length index given in first level." << std::endl;
    exit(0);
  }
  
  // Set the variable values
  setId("VertexFromTRBS");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  // Set the parameter identities
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "Y_mod";
  tmp[1] = "P_ratio";
  setParameterId( tmp );
}

void VertexFromTRBS::
derivs(Tissue &T,
       std::vector< std::vector<double> > &cellData,
       std::vector< std::vector<double> > &wallData,
       std::vector< std::vector<double> > &vertexData,
       std::vector< std::vector<double> > &cellDerivs,
       std::vector< std::vector<double> > &wallDerivs,
       std::vector< std::vector<double> > &vertexDerivs ) {
  
  //Do the update for each cell
  size_t numCells = T.numCell();
  size_t wallLengthIndex = variableIndex(0,0);
  size_t numWalls = 3; // defined only for triangles at the moment
  
  for( size_t i=0 ; i<numCells ; ++i ) {
    if( T.cell(i).numWall() != numWalls ) {
      std::cerr << "VertexFromTRBS::derivs() only defined for triangular cells."
		<< std::endl;
      exit(-1);
    }
    
    double young = parameter(0);
    double poisson =parameter(1);
    size_t v1 = T.cell(i).vertex(0)->index();
    size_t v2 = T.cell(i).vertex(1)->index();
    size_t v3 = T.cell(i).vertex(2)->index();
    size_t w1 = T.cell(i).wall(0)->index();
    size_t w2 = T.cell(i).wall(1)->index();
    size_t w3 = T.cell(i).wall(2)->index();
    std::vector<double> restingLength(numWalls);
    restingLength[0] = wallData[w1][wallLengthIndex];
    restingLength[1] = wallData[w2][wallLengthIndex];
    restingLength[2] = wallData[w3][wallLengthIndex];

    std::vector< std::vector<double> > position(3,vertexData[v1]);
    position[1] = vertexData[v2];
    position[2] = vertexData[v3];
    //position[0][2] z for vertex 1 (of the cell)
   
    std::vector<double> length(numWalls);
    length[0] = T.wall(w1).lengthFromVertexPosition(vertexData);
    length[1] = T.wall(w2).lengthFromVertexPosition(vertexData);
    length[2] = T.wall(w3).lengthFromVertexPosition(vertexData);
    
    // Lame coefficients (can be defined out of loop)
    double lambda=young*poisson/(1-poisson*poisson);
    double mio=young/(1+poisson);

    //Area of the element (Heron's formula)                                    
    double Area=std::sqrt( ( restingLength[0]+restingLength[1]+restingLength[2])*
                           (-restingLength[0]+restingLength[1]+restingLength[2])*
                           ( restingLength[0]-restingLength[1]+restingLength[2])*
                           ( restingLength[0]+restingLength[1]-restingLength[2])  )*0.25;
    
    //Angles of the element ( assuming the order: 0,L0,1,L1,2,L2 )
    std::vector<double> Angle(3);
    Angle[0]=std::asin(2*Area/(restingLength[0]*restingLength[2]));
    // can be ommited by cotan(A)=.25*sqrt(4*b*b*c*c/K-1)
    Angle[1]=std::asin(2*Area/(restingLength[0]*restingLength[1]));        
    Angle[2]=std::asin(2*Area/(restingLength[1]*restingLength[2]));
    
    //Tensile Stiffness
    double tensileStiffness[3];
    double const temp = 1.0/(Area*16);                                      
    double cotan[3] = {1.0/std::tan(Angle[0]),1.0/std::tan(Angle[1]),1.0/std::tan(Angle[2])};    
    tensileStiffness[0]=(2*cotan[2]*cotan[2]*(lambda+mio)+mio)*temp;
    tensileStiffness[1]=(2*cotan[0]*cotan[0]*(lambda+mio)+mio)*temp;
    tensileStiffness[2]=(2*cotan[1]*cotan[1]*(lambda+mio)+mio)*temp;
    
    //Angular Stiffness
    double angularStiffness[3];
    angularStiffness[0]=(2*cotan[1]*cotan[2]*(lambda+mio)-mio)*temp;
    angularStiffness[1]=(2*cotan[0]*cotan[2]*(lambda+mio)-mio)*temp;
    angularStiffness[2]=(2*cotan[0]*cotan[1]*(lambda+mio)-mio)*temp;
    
    //Calculate biquadratic strains  
    std::vector<double> Delta(3);
    Delta[0]=(length[0])*(length[0])-(restingLength[0])*(restingLength[0]);
    Delta[1]=(length[1])*(length[1])-(restingLength[1])*(restingLength[1]);
    Delta[2]=(length[2])*(length[2])-(restingLength[2])*(restingLength[2]);
    //Forces of vertices
    double Force[3][3];                                           
          
    Force[0][0]= (tensileStiffness[0]*Delta[0]+angularStiffness[1]*Delta[1]+angularStiffness[0]*Delta[2])*(position[1][0]-position[0][0])
                +(tensileStiffness[2]*Delta[2]+angularStiffness[2]*Delta[1]+angularStiffness[0]*Delta[0])*(position[2][0]-position[0][0]); 
    Force[0][1]= (tensileStiffness[0]*Delta[0]+angularStiffness[1]*Delta[1]+angularStiffness[0]*Delta[2])*(position[1][1]-position[0][1])
                +(tensileStiffness[2]*Delta[2]+angularStiffness[2]*Delta[1]+angularStiffness[0]*Delta[0])*(position[2][1]-position[0][1]); 
    Force[0][2]= (tensileStiffness[0]*Delta[0]+angularStiffness[1]*Delta[1]+angularStiffness[0]*Delta[2])*(position[1][2]-position[0][2])
                +(tensileStiffness[2]*Delta[2]+angularStiffness[2]*Delta[1]+angularStiffness[0]*Delta[0])*(position[2][2]-position[0][2]); 

    Force[1][0]= (tensileStiffness[0]*Delta[0]+angularStiffness[0]*Delta[2]+angularStiffness[1]*Delta[1])*(position[0][0]-position[1][0])
                +(tensileStiffness[1]*Delta[1]+angularStiffness[2]*Delta[2]+angularStiffness[1]*Delta[0])*(position[2][0]-position[1][0]); 
    Force[1][1]= (tensileStiffness[0]*Delta[0]+angularStiffness[0]*Delta[2]+angularStiffness[1]*Delta[1])*(position[0][1]-position[1][1])
                +(tensileStiffness[1]*Delta[1]+angularStiffness[2]*Delta[2]+angularStiffness[1]*Delta[0])*(position[2][1]-position[1][1]); 
    Force[1][2]= (tensileStiffness[0]*Delta[0]+angularStiffness[0]*Delta[2]+angularStiffness[1]*Delta[1])*(position[0][2]-position[1][2])
                +(tensileStiffness[1]*Delta[1]+angularStiffness[2]*Delta[2]+angularStiffness[1]*Delta[0])*(position[2][2]-position[1][2]); 

    Force[2][0]= (tensileStiffness[2]*Delta[2]+angularStiffness[0]*Delta[0]+angularStiffness[2]*Delta[1])*(position[0][0]-position[2][0])
                +(tensileStiffness[1]*Delta[1]+angularStiffness[1]*Delta[0]+angularStiffness[2]*Delta[2])*(position[1][0]-position[2][0]); 
    Force[2][1]= (tensileStiffness[2]*Delta[2]+angularStiffness[0]*Delta[0]+angularStiffness[2]*Delta[1])*(position[0][1]-position[2][1])
                +(tensileStiffness[1]*Delta[1]+angularStiffness[1]*Delta[0]+angularStiffness[2]*Delta[2])*(position[1][1]-position[2][1]); 
    Force[2][2]= (tensileStiffness[2]*Delta[2]+angularStiffness[0]*Delta[0]+angularStiffness[2]*Delta[1])*(position[0][2]-position[2][2])
                +(tensileStiffness[1]*Delta[1]+angularStiffness[1]*Delta[0]+angularStiffness[2]*Delta[2])*(position[1][2]-position[2][2]); 

    // adding TRBS forces to the total vertexDerivs
    
    vertexDerivs[v1][0]+= Force[0][0];
    vertexDerivs[v1][1]+= Force[0][1];
    vertexDerivs[v1][2]+= Force[0][2];
    
    vertexDerivs[v2][0]+= Force[1][0];
    vertexDerivs[v2][1]+= Force[1][1];
    vertexDerivs[v2][2]+= Force[1][2];
    
    vertexDerivs[v3][0]+= Force[2][0];
    vertexDerivs[v3][1]+= Force[2][1];
    vertexDerivs[v3][2]+= Force[2][2];
  }
}
