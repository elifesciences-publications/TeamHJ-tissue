//
// Filename     : mechanicalTRBS.cc
// Description  : Classes describing updates due to mechanical triangular biquadratic springs
// Author(s)    : Behruz Bozorg, Henrik Jonsson (henrik@thep.lu.se)
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
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) {
  
  //Do the update for each cell
  size_t numCells = T.numCell();
  size_t wallLengthIndex = variableIndex(0,0);
  size_t numWalls = 3; // defined only for triangles at the moment
  
  for( size_t i=0 ; i<numCells ; ++i ) {
    if( T.cell(i).numWall() != numWalls ) {
      std::cerr << "VertexFromTRBS::derivs() only defined for triangular cells."
		<< " Not for cells with " << T.cell(i).numWall() << " walls!"
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

    DataMatrix position(3,vertexData[v1]);
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
    
    // Area of the element (using Heron's formula)                                      
    double Area=std::sqrt( ( restingLength[0]+restingLength[1]+restingLength[2])*
                           (-restingLength[0]+restingLength[1]+restingLength[2])*
                           ( restingLength[0]-restingLength[1]+restingLength[2])*
                           ( restingLength[0]+restingLength[1]-restingLength[2])  )*0.25;
    
    //Angles of the element ( assuming the order: 0,L0,1,L1,2,L2 )
    std::vector<double> Angle(3);
    Angle[0]=std::acos(  (restingLength[0]*restingLength[0]+restingLength[2]*restingLength[2]-restingLength[1]*restingLength[1])/
                         (restingLength[0]*restingLength[2]*2)    );
    Angle[1]=std::acos(  (restingLength[0]*restingLength[0]+restingLength[1]*restingLength[1]-restingLength[2]*restingLength[2])/
                         (restingLength[0]*restingLength[1]*2)    );
    Angle[2]=std::acos(  (restingLength[1]*restingLength[1]+restingLength[2]*restingLength[2]-restingLength[0]*restingLength[0])/
                         (restingLength[1]*restingLength[2]*2)    );
    // can be ommited by cotan(A)=.25*sqrt(4*b*b*c*c/K-1)
    
    //Tensile Stiffness
    double tensileStiffness[3];
    double const temp = 1.0/(Area*16);                                      
    double cotan[3] = {1.0/std::tan(Angle[0]),1.0/std::tan(Angle[1]),1.0/std::tan(Angle[2])};    
    tensileStiffness[0]=(2*cotan[2]*cotan[2]*(lambda+mio)+mio)*temp;
    tensileStiffness[1]=(2*cotan[0]*cotan[0]*(lambda+mio)+mio)*temp;
    tensileStiffness[2]=(2*cotan[1]*cotan[1]*(lambda+mio)+mio)*temp;
    
    //std::cerr <<"young, parameter0, poisson, parameter1, lambda and mio :  " 
    //          << young <<" , " << parameter(0) <<" "<< poisson <<" , " << parameter(1)<< " "<< lambda <<" , "<< mio << std::endl;
    
    //std::cerr <<" angles, element " << i << std::endl;                                   
    //std::cerr <<"[0] "<< Angle[0] << std::endl
    //          <<"[1] "<< Angle[1] << std::endl
    //          <<"[2] "<< Angle[2] << std::endl << std::endl;
    
    
    //std::cerr <<"cotan of angles, element " << i << std::endl;     
    //std::cerr <<"[0] "<< cotan[0] << std::endl
    //          <<"[1] "<< cotan[1] << std::endl
    //          <<"[2] "<< cotan[2] << std::endl << std::endl;

    //std::cerr <<"tan of angles, element " << i << std::endl;                                    
    //std::cerr <<"[0] "<< std::tan(Angle[0]) << std::endl
    //          <<"[1] "<< std::tan(Angle[1]) << std::endl
    //          <<"[2] "<< std::tan(Angle[2]) << std::endl << std::endl;


    //std::cerr <<"Tensile stiffness, element " << i << std::endl;                                   
    //std::cerr <<"[0] "<< tensileStiffness[0] << std::endl
    //	        <<"[1] "<< tensileStiffness[1] << std::endl
    // 	        <<"[2] "<< tensileStiffness[2] << std::endl << std::endl;

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

    //std::cerr <<"Forces on Vertices , Element  "<< i << std::endl;                                                  //<<<<<<<<<<
    //std::cerr <<"vertex [0] "<< Force[0][0] <<" "<<Force[0][1] <<" "<<Force[0][2] << std::endl
    //          <<"vertex [1] "<< Force[1][0] <<" "<<Force[1][1] <<" "<<Force[1][2] << std::endl
    //          <<"vertex [2] "<< Force[2][0] <<" "<<Force[2][1] <<" "<<Force[2][2] << std::endl<< std::endl;
   

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


VertexFromTRBScenterTriangulation::
VertexFromTRBScenterTriangulation(std::vector<double> &paraValue, 
	       std::vector< std::vector<size_t> > 
	       &indValue ) 
{  
  // Do some checks on the parameters and variable indeces
  if( paraValue.size()!=2 ) {
    std::cerr << "VertexFromTRBScenterTriangulation::"
	      << "VertexFromTRBScenterTriangulation() "
	      << "Uses two parameters young modulus and poisson coefficient.\n";
    exit(0);
  }
  if( (indValue.size()!=2 && indValue.size()!=4) || 
      indValue[0].size()!=1 || indValue[1].size()!=1 ||
      (indValue.size()==4 && (indValue[2].size()!=0 && indValue[2].size()!=1)) ||
      (indValue.size()==4 && (indValue[3].size()!=0 && indValue[3].size()!=1)) 
      ){
    std::cerr << "VertexFromTRBScenterTriangulation::"
	      << "VertexFromTRBScenterTriangulation() "
	      << "Wall length index is given in first level." 
	      << "Start of additional Cell variable indices (center(x,y,z) "
	      << "L_1,...,L_n, n=num vertex) is given in second level (typically at end)." 
              << "Optionally two additional levels can be given where the strain and stress "
	      << "directions can be stored at given indices. If index given at third level, "
	      << "strain direction will be stored starting at this (cell) variable index, "
	      << "and for fourth level stress will be stored."
	      << std::endl;
    exit(0);
  }
  
  // Set the variable values
  setId("VertexFromTRBScenterTriangulation");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  // Set the parameter identities
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "Y_mod";
  tmp[1] = "P_ratio";
  setParameterId( tmp );
}


void VertexFromTRBScenterTriangulation::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) {
  
  //Do the update for each cell
  size_t dimension = 3;
  assert (dimension==vertexData[0].size());
  size_t numCells = T.numCell();
  size_t wallLengthIndex = variableIndex(0,0);
  size_t comIndex = variableIndex(1,0);
  size_t lengthInternalIndex = comIndex+dimension;
  
  
  for (size_t cellIndex=0 ; cellIndex<numCells ; ++cellIndex) {
    size_t numWalls = T.cell(cellIndex).numWall(); 
    
    if(  T.cell(cellIndex).numVertex()!= numWalls ) {
      std::cerr << "VertexFromTRBScenterTriangulation::derivs() same number of vertices and walls."
		<< " Not for cells with " << T.cell(cellIndex).numWall() << " walls and "
		<< T.cell(cellIndex).numVertex() << " vertices!"	
		<< std::endl;
      exit(-1);
    }
    
    double young = parameter(0);
    double poisson =parameter(1);
    
    double StrainCellGlobal[3][3]={{0,0,0},{0,0,0},{0,0,0}};
    double StressCellGlobal[3][3]={{0,0,0},{0,0,0},{0,0,0}};
    double TotalCellRestingArea=0;

    // One triangle per 'vertex' in cyclic order
    for (size_t k=0; k<numWalls; ++k) { 
      size_t kPlusOneMod = (k+1)%numWalls;
      //size_t v1 = com;
      size_t v2 = T.cell(cellIndex).vertex(k)->index();
      size_t v3 = T.cell(cellIndex).vertex(kPlusOneMod)->index();
      //size_t w1 = internal k
      size_t w2 = T.cell(cellIndex).wall(k)->index();
      //size_t w3 = internal k+1

      // Position matrix holds in rows positions for com, vertex(k), vertex(k+1)
      DataMatrix position(3,vertexData[v2]);
      for (size_t d=0; d<dimension; ++d)
	position[0][d] = cellData[cellIndex][comIndex+d]; // com position
      //position[1] = vertexData[v2]; // given by initiation
      position[2] = vertexData[v3];
      
      // Resting lengths are from com-vertex(k), vertex(k)-vertex(k+1) (wall(k)), com-vertex(k+1)
      std::vector<double> restingLength(numWalls);
      restingLength[0] = cellData[cellIndex][lengthInternalIndex + k];
      restingLength[1] = wallData[w2][wallLengthIndex];
      restingLength[2] = cellData[cellIndex][lengthInternalIndex + kPlusOneMod];
      
      // Lengths are from com-vertex(k), vertex(k)-vertex(k+1) (wall(k)), com-vertex(k+1)
      std::vector<double> length(numWalls);
      length[0] = std::sqrt( (position[0][0]-position[1][0])*(position[0][0]-position[1][0]) +
			     (position[0][1]-position[1][1])*(position[0][1]-position[1][1]) +
			     (position[0][2]-position[1][2])*(position[0][2]-position[1][2]) );
      
      length[1] = T.wall(w2).lengthFromVertexPosition(vertexData);
        
      length[2] = std::sqrt( (position[0][0]-position[2][0])*(position[0][0]-position[2][0]) +
			     (position[0][1]-position[2][1])*(position[0][1]-position[2][1]) +
			     (position[0][2]-position[2][2])*(position[0][2]-position[2][2]) );

            
      // Lame coefficients (can be defined out of loop)
      double lambda=young*poisson/(1-poisson*poisson);
      double mio=young/(1+poisson);
      
      // resting Area of the element (using Heron's formula)                                      
      double restingArea=std::sqrt( ( restingLength[0]+restingLength[1]+restingLength[2])*
                                    (-restingLength[0]+restingLength[1]+restingLength[2])*
                                    ( restingLength[0]-restingLength[1]+restingLength[2])*
                                    ( restingLength[0]+restingLength[1]-restingLength[2])  )*0.25;
            

      //Angles of the element ( assuming the order: 0,L0,1,L1,2,L2 )
      std::vector<double> Angle(3);
       // can be ommited by cotan(A)=.25*sqrt(4*b*b*c*c/K-1)
      Angle[0]=std::acos(  (restingLength[0]*restingLength[0]+restingLength[2]*restingLength[2]-restingLength[1]*restingLength[1])/
                           (restingLength[0]*restingLength[2]*2)    );
      Angle[1]=std::acos(  (restingLength[0]*restingLength[0]+restingLength[1]*restingLength[1]-restingLength[2]*restingLength[2])/
                           (restingLength[0]*restingLength[1]*2)    );
      Angle[2]=std::acos(  (restingLength[1]*restingLength[1]+restingLength[2]*restingLength[2]-restingLength[0]*restingLength[0])/
                           (restingLength[1]*restingLength[2]*2)    );
      
      //Tensile Stiffness
      double tensileStiffness[3];
      double temp = 1.0/(restingArea*16);                                      
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
  
      //Area of the element (using Heron's formula)                                      
      double Area=std::sqrt( ( length[0]+length[1]+length[2])*
                             (-length[0]+length[1]+length[2])*
                             ( length[0]-length[1]+length[2])*
                             ( length[0]+length[1]-length[2])  )*0.25;
        

      //Current shape local coordinate of the element  (counterclockwise ordering of nodes/edges)
      double CurrentAngle1=std::acos(  (length[0]*length[0]+length[1]*length[1]-length[2]*length[2])/
                                       (length[0]*length[1]*2)    );

      double Qa=std::cos(CurrentAngle1)*length[0];
      double Qc=std::sin(CurrentAngle1)*length[0];
      double Qb=length[1];
      // shape vector matrix = inverse of coordinate matrix ( only first two elements i.e. ShapeVector[3][2] )      
      // double ShapeVectorCurrent[3][3]={ {  0   ,       1/Qc      , 0 }, 
      //                                   {-1/Qb , (Qa-Qb)/(Qb*Qc) , 1 },       
      //                                   { 1/Qb ,     -Qa/(Qb*Qc) , 0 }  };
            
            
      //Shape vectors in current (clockwise ordering of nodes/edges)
      // double CurrentAngle2=std::acos(  (length[1]*length[1]+length[2]*length[2]-length[0]*length[0])/
      //                                  (length[1]*length[2]*2)    );

      // double Qa=std::cos(CurrentAngle2)*length[2];
      // double Qb=length[1];
      // double Qc=std::sin(CurrentAngle2)*length[2];
      
      // double ShapeVectorCurrent[3][2]={ {  0   ,       1/Qc      }, 
      //                                   { 1/Qb ,     -Qa/(Qb*Qc) },       
      //                                   {-1/Qb , (Qa-Qb)/(Qb*Qc) }  };

      // Local coordinates of the resting shape ( counterclockwise )
      double RestingAngle1=std::acos(  (restingLength[0]*restingLength[0]+restingLength[1]*restingLength[1]-restingLength[2]*restingLength[2])/
                                       (restingLength[0]*restingLength[1]*2)    );

      double Pa=std::cos(RestingAngle1)*restingLength[0];
      double Pc=std::sin(RestingAngle1)*restingLength[0];
      double Pb=restingLength[1];

      // shape vector matrix in resting shape in local coordinate system  = inverse of coordinate matrix ( only first two elements i.e. ShapeVectorResting[3][2] )      
      double ShapeVectorResting[3][3]={ {  0   ,       1/Pc      , 0 }, 
                                        {-1/Pb , (Pa-Pb)/(Pb*Pc) , 1 },       
                                        { 1/Pb ,     -Pa/(Pb*Pc) , 0 }  };

      // Local coordinates of the resting shape (clockwise )
      //....
      //....

      //square of radius of circumstancing circle in resting shape
      //double Rcirc2=(0.25*restingLength[0]*restingLength[1]*restingLength[2]/Area)*(0.25*restingLength[0]*restingLength[1]*restingLength[2]/Area);  
      

 // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> STRAIN and STRESS TENSOR (BEGIN) <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      
      // deformation gradiant tensor F =Sigma i=1,2,3 Qi x Di
      // strain tensor in resting shape E=0.5(FtF-I)
      // trE
      // B=FFt
      double trE=( Delta[1]*cotan[0]+ Delta[2]*cotan[1]+Delta[0]*cotan[2])/(4*restingArea);
            
      double positionLocal[3][2]={ {Qa , Qc}, 
                                   {0  , 0 },  
                                   {Qb , 0 }  };
      
      double DeformGrad[2][2]={{0,0},{0,0}}; // F= Qi x Di
      for ( int i=0 ; i<3 ; ++i ) {
        DeformGrad[0][0]=DeformGrad[0][0]+positionLocal[i][0]*ShapeVectorResting[i][0];
        DeformGrad[1][0]=DeformGrad[1][0]+positionLocal[i][1]*ShapeVectorResting[i][0];
        DeformGrad[0][1]=DeformGrad[0][1]+positionLocal[i][0]*ShapeVectorResting[i][1];
        DeformGrad[1][1]=DeformGrad[1][1]+positionLocal[i][1]*ShapeVectorResting[i][1];
      } 

      double LeftCauchy[2][2]; // B=FFt
      LeftCauchy[0][0]=DeformGrad[0][0]*DeformGrad[0][0]+DeformGrad[0][1]*DeformGrad[0][1];
      LeftCauchy[1][0]=DeformGrad[1][0]*DeformGrad[0][0]+DeformGrad[1][1]*DeformGrad[0][1];
      LeftCauchy[0][1]=DeformGrad[0][0]*DeformGrad[1][0]+DeformGrad[0][1]*DeformGrad[1][1];
      LeftCauchy[1][1]=DeformGrad[1][0]*DeformGrad[1][0]+DeformGrad[1][1]*DeformGrad[1][1];


      double Egreen[2][2];//E=0.5(C-I)
      Egreen[0][0]=0.5*(DeformGrad[0][0]*DeformGrad[0][0]+DeformGrad[1][0]*DeformGrad[1][0]-1);
      Egreen[1][0]=0.5*(DeformGrad[0][1]*DeformGrad[0][0]+DeformGrad[1][1]*DeformGrad[1][0]);
      Egreen[0][1]=0.5*(DeformGrad[0][0]*DeformGrad[0][1]+DeformGrad[1][0]*DeformGrad[1][1]);
      Egreen[1][1]=0.5*(DeformGrad[0][1]*DeformGrad[0][1]+DeformGrad[1][1]*DeformGrad[1][1]-1);

      
      double StrainAlmansi[2][2]; // e=0.5(1-B^-1)  True strain tensor
      temp=LeftCauchy[0][0]*LeftCauchy[1][1]-LeftCauchy[1][0]*LeftCauchy[0][1]; // det(B)
      StrainAlmansi[0][0]=0.5*(1-(LeftCauchy[1][1]/temp));
      StrainAlmansi[1][0]=0.5*LeftCauchy[1][0]/temp;
      StrainAlmansi[0][1]=0.5*LeftCauchy[0][1]/temp;  
      StrainAlmansi[1][1]=0.5*(1-(LeftCauchy[0][0]/temp));
      
      
      
      double B2[2][2];// LeftCauchy^2
      B2[0][0]=LeftCauchy[0][0]*LeftCauchy[0][0]+LeftCauchy[0][1]*LeftCauchy[1][0];
      B2[1][0]=LeftCauchy[1][0]*LeftCauchy[0][0]+LeftCauchy[1][1]*LeftCauchy[1][0];
      B2[0][1]=LeftCauchy[0][0]*LeftCauchy[0][1]+LeftCauchy[0][1]*LeftCauchy[1][1];
      B2[1][1]=LeftCauchy[1][0]*LeftCauchy[0][1]+LeftCauchy[1][1]*LeftCauchy[1][1];

      double StressTensor[3][3]; // true stress tensor (isotropic term) based on lambdaT and mioT
      StressTensor[0][0]=(Area/restingArea)*((lambda*trE-mio/2)*LeftCauchy[0][0]+(mio/2)*B2[0][0]);
      StressTensor[1][0]=(Area/restingArea)*((lambda*trE-mio/2)*LeftCauchy[1][0]+(mio/2)*B2[1][0]);
      StressTensor[0][1]=(Area/restingArea)*((lambda*trE-mio/2)*LeftCauchy[0][1]+(mio/2)*B2[0][1]);
      StressTensor[1][1]=(Area/restingArea)*((lambda*trE-mio/2)*LeftCauchy[1][1]+(mio/2)*B2[1][1]);


      // std::cerr <<"stress tensor " << std::endl;
      // std::cerr <<" Sxx  "<< StressTensor[0][0] <<" Sxy  "<< StressTensor[0][1] <<" Sxz  "<< StressTensor[0][2] << std::endl
      //           <<" Syx  "<< StressTensor[1][0] <<" Syy  "<< StressTensor[1][1] <<" Syz  "<< StressTensor[1][2] << std::endl
      //           <<" Szx  "<< StressTensor[2][0] <<" Szy  "<< StressTensor[2][1] <<" Szz  "<< StressTensor[2][2] << std::endl <<std::endl;


      //Shape vectors in Current shape (counterclockwise ordering of nodes/edges)     ShapeVectorCurrent[3][3]  calculated above   
    
          

      //square of radius of circumstancing circle in resting shape
      //double Rcirc2Resting=(0.25*restingLength[0]*restingLength[1]*restingLength[2]/Area)*(0.25*restingLength[0]*restingLength[1]*restingLength[2]/Area);  
      
      double StrainTensor[3][3]; // there are other alternatives than StrainAlmansi for strain tensor
      StrainTensor[0][0]=StrainAlmansi[0][0];
      StrainTensor[1][0]=StrainAlmansi[1][0];
      StrainTensor[0][1]=StrainAlmansi[0][1];
      StrainTensor[1][1]=StrainAlmansi[1][1];


      StrainTensor[0][2]=0;  // adding 3rd dimension which is zero, the tensor is still in element plane
      StrainTensor[1][2]=0;
      StrainTensor[2][2]=0;
      StrainTensor[2][0]=0;
      StrainTensor[2][1]=0;
      
      StressTensor[0][2]=0;  // adding 3rd dimension which is zero, the tensor is still in element plane
      StressTensor[1][2]=0;
      StressTensor[2][2]=0;
      StressTensor[2][0]=0;
      StressTensor[2][1]=0;

      //rotation matrix for going from local to global coordinate system based on counterclockwise ordering
      double rotation[3][3];  

      double tempA=std::sqrt((position[2][0]-position[1][0])*(position[2][0]-position[1][0])+
                             (position[2][1]-position[1][1])*(position[2][1]-position[1][1])+
                             (position[2][2]-position[1][2])*(position[2][2]-position[1][2])  );

      double tempB=std::sqrt((position[0][0]-position[1][0])*(position[0][0]-position[1][0])+
                             (position[0][1]-position[1][1])*(position[0][1]-position[1][1])+
                             (position[0][2]-position[1][2])*(position[0][2]-position[1][2])  );

      double Xcurrent[3];      
      Xcurrent[0]= (position[2][0]-position[1][0])/tempA;
      Xcurrent[1]= (position[2][1]-position[1][1])/tempA;
      Xcurrent[2]= (position[2][2]-position[1][2])/tempA;
      
      double Bcurrent[3];      
      Bcurrent[0]= (position[0][0]-position[1][0])/tempB;
      Bcurrent[1]= (position[0][1]-position[1][1])/tempB;
      Bcurrent[2]= (position[0][2]-position[1][2])/tempB;

      double Zcurrent[3];      
      Zcurrent[0]= Xcurrent[1]*Bcurrent[2]-Xcurrent[2]*Bcurrent[1];
      Zcurrent[1]= Xcurrent[2]*Bcurrent[0]-Xcurrent[0]*Bcurrent[2];
      Zcurrent[2]= Xcurrent[0]*Bcurrent[1]-Xcurrent[1]*Bcurrent[0];
      
      tempA=std:: sqrt(Zcurrent[0]*Zcurrent[0]+Zcurrent[1]*Zcurrent[1]+Zcurrent[2]*Zcurrent[2]);
      Zcurrent[0]=Zcurrent[0]/tempA;
      Zcurrent[1]=Zcurrent[1]/tempA;
      Zcurrent[2]=Zcurrent[2]/tempA;

      double Ycurrent[3];      
      Ycurrent[0]= Zcurrent[1]*Xcurrent[2]-Zcurrent[2]*Xcurrent[1];
      Ycurrent[1]= Zcurrent[2]*Xcurrent[0]-Zcurrent[0]*Xcurrent[2];
      Ycurrent[2]= Zcurrent[0]*Xcurrent[1]-Zcurrent[1]*Xcurrent[0];


      rotation[0][0]=Xcurrent[0];
      rotation[1][0]=Xcurrent[1];
      rotation[2][0]=Xcurrent[2];

      rotation[0][1]=Ycurrent[0];
      rotation[1][1]=Ycurrent[1];
      rotation[2][1]=Ycurrent[2];

      rotation[0][2]=Zcurrent[0];
      rotation[1][2]=Zcurrent[1];
      rotation[2][2]=Zcurrent[2];      

      // rotariong strain tensor to the global coordinate system
      double tempR[3][3]={{0,0,0},{0,0,0},{0,0,0}};
      for (int r=0 ; r<3 ; r++) {
        for (int s=0 ; s<3 ; s++) {
          for(int w=0 ; w<3 ; w++) {
            tempR[r][s]=tempR[r][s]+rotation[r][w]*StrainTensor[w][s];
          }
        }
      }
      for (int r=0 ; r<3 ; r++) {
        for (int s=0 ; s<3 ; s++) {
          StrainTensor[r][s]=0;
          for(int w=0 ; w<3 ; w++) {
            StrainTensor[r][s]=StrainTensor[r][s]+tempR[r][w]*rotation[s][w];
          }
        }
      }

      // rotating stress tensor to the global coordinate system
      for (int r=0 ; r<3 ; r++) 
        for (int s=0 ; s<3 ; s++) 
          tempR[r][s]=0;
        
      for (int r=0 ; r<3 ; r++) 
        for (int s=0 ; s<3 ; s++) 
          for(int w=0 ; w<3 ; w++) 
            tempR[r][s]=tempR[r][s]+rotation[r][w]*StressTensor[w][s];
          
      for (int r=0 ; r<3 ; r++) 
        for (int s=0 ; s<3 ; s++) {
          StressTensor[r][s]=0;
          for(int w=0 ; w<3 ; w++) 
            StressTensor[r][s]=StressTensor[r][s]+tempR[r][w]*rotation[s][w]; 
        }


      // std::cerr <<" Sxx  "<< StrainTensor[0][0] <<" Sxy  "<< StrainTensor[0][1] << std::endl
      //           <<" Syx  "<< StrainTensor[1][0] <<" Syy  "<< StrainTensor[1][1] << std::endl<<std::endl;
     
      for (int r=0 ; r<3 ; r++) {
        for (int s=0 ; s<3 ; s++) {   
          StrainCellGlobal[r][s]= StrainCellGlobal[r][s]+restingArea*StrainTensor[r][s];
        }
      }
      
      for (int r=0 ; r<3 ; r++) 
        for (int s=0 ; s<3 ; s++)    
          StressCellGlobal[r][s]= StressCellGlobal[r][s]+restingArea*StressTensor[r][s];
      
      TotalCellRestingArea=TotalCellRestingArea+restingArea;
    
      // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> STRAIN and STRESS TENSORS (END) <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

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
      
      // adding forces to the total vertexDerivs
      
      cellDerivs[cellIndex][comIndex  ] += Force[0][0];
      cellDerivs[cellIndex][comIndex+1] += Force[0][1];
      cellDerivs[cellIndex][comIndex+2] += Force[0][2];
      
      vertexDerivs[v2][0] += Force[1][0];
      vertexDerivs[v2][1] += Force[1][1];
      vertexDerivs[v2][2] += Force[1][2];
      
      vertexDerivs[v3][0] += Force[2][0];
      vertexDerivs[v3][1] += Force[2][1];
      vertexDerivs[v3][2] += Force[2][2];
    }

    for (int r=0 ; r<3 ; r++) 
      for (int s=0 ; s<3 ; s++)    
        StrainCellGlobal[r][s]= StrainCellGlobal[r][s]/TotalCellRestingArea;

    for (int r=0 ; r<3 ; r++) 
      for (int s=0 ; s<3 ; s++)    
        StressCellGlobal[r][s]= StressCellGlobal[r][s]/TotalCellRestingArea; 
   
    
    // eigenvalue/eigenvectors of averaged STRAIN and STRESS tensors in global coordinate system. (Jacobi method)

    // STRAIN:

    double eigenVectorStrain[3][3]={{1,0,0},{0,1,0},{0,0,1}};
    double pivot=1;
    double pi=3.1415;
    int I,J;
    double RotAngle,Si,Co;
    while (pivot>0.00001) {
      pivot=std::abs(StrainCellGlobal[1][0]);
      I=1;
      J=0;
      if (std::abs(StrainCellGlobal[2][0])>pivot) {
        pivot=std::abs(StrainCellGlobal[2][0]);
        I=2;
        J=0;
      }
      if (std::abs(StrainCellGlobal[2][1])>pivot) {
        pivot=std::abs(StrainCellGlobal[2][1]);
        I=2;
        J=1;
      }
      if (std::abs(StrainCellGlobal[I][I]-StrainCellGlobal[J][J])<0.00001) {
          RotAngle=pi/4;
      }            
      else {
        RotAngle=0.5*std::atan((2*StrainCellGlobal[I][J])/(StrainCellGlobal[J][J]-StrainCellGlobal[I][I]));
      }
        Si=std::sin(RotAngle);
        Co=std::cos(RotAngle);
        double tempRot[3][3]={{1,0,0},{0,1,0},{0,0,1}};
        tempRot[I][I]=Co;
        tempRot[J][J]=Co;
        tempRot[I][J]=Si;
        tempRot[J][I]=-Si;
        double tempStrain[3][3]={{0,0,0},{0,0,0},{0,0,0}};
        for (int r=0 ; r<3 ; r++) 
          for (int s=0 ; s<3 ; s++) 
            for(int w=0 ; w<3 ; w++) 
              tempStrain[r][s]=tempStrain[r][s]+StrainCellGlobal[r][w]*tempRot[w][s];
        
        for (int r=0 ; r<3 ; r++) 
          for (int s=0 ; s<3 ; s++) 
            StrainCellGlobal[r][s]=0;
         
        for (int r=0 ; r<3 ; r++) 
          for (int s=0 ; s<3 ; s++) 
            for(int w=0 ; w<3 ; w++) 
              StrainCellGlobal[r][s]=StrainCellGlobal[r][s]+tempRot[w][r]*tempStrain[w][s];
                
        for (int r=0 ; r<3 ; r++) 
          for (int s=0 ; s<3 ; s++) 
            tempStrain[r][s]=eigenVectorStrain[r][s];
        
        for (int r=0 ; r<3 ; r++) 
          for (int s=0 ; s<3 ; s++) 
            eigenVectorStrain[r][s]=0;
        
        for (int r=0 ; r<3 ; r++) 
          for (int s=0 ; s<3 ; s++) 
            for(int w=0 ; w<3 ; w++) 
              eigenVectorStrain[r][s]=eigenVectorStrain[r][s]+tempStrain[r][w]*tempRot[w][s];
            
    }
       
      
      // maximal strain direction
      double maximalStrainValue=StrainCellGlobal[0][0];
      int Istrain=0;
      if (StrainCellGlobal[1][1]>maximalStrainValue) 
        {
          maximalStrainValue=StrainCellGlobal[1][1];
          Istrain=1;
        }
      if (StrainCellGlobal[2][2]>maximalStrainValue) 
        {
          maximalStrainValue=StrainCellGlobal[2][2];
          Istrain=2;
        }
      // std::cerr<<"maximal Strain direction "<< eigenVectorStrain[0][Istrain] <<" "<< eigenVectorStrain[1][Istrain] <<" "<< eigenVectorStrain[2][Istrain] <<std::endl;  
      // std::cerr<<"maximal Strain value "<< maximalStrainValue <<std::endl;  
      
      // STRESS:
      
      double eigenVectorStress[3][3]={{1,0,0},{0,1,0},{0,0,1}};
      pivot=1;
      //double RotAngle,Si,Co;
      while (pivot>0.00001) {
        pivot=std::abs(StressCellGlobal[1][0]);
        I=1;
        J=0;
        if (std::abs(StressCellGlobal[2][0])>pivot) {
          pivot=std::abs(StressCellGlobal[2][0]);
          I=2;
          J=0;
          }
        if (std::abs(StressCellGlobal[2][1])>pivot) {
          pivot=std::abs(StressCellGlobal[2][1]);
          I=2;
          J=1;
        }
        if (std::abs(StressCellGlobal[I][I]-StressCellGlobal[J][J])<0.00001) {
          RotAngle=pi/4;
        }            
        else {
          RotAngle=0.5*std::atan((2*StressCellGlobal[I][J])/(StressCellGlobal[J][J]-StressCellGlobal[I][I]));
        }
        Si=std::sin(RotAngle);
        Co=std::cos(RotAngle);
        double tempRot[3][3]={{1,0,0},{0,1,0},{0,0,1}};
        tempRot[I][I]=Co;
        tempRot[J][J]=Co;
        tempRot[I][J]=Si;
        tempRot[J][I]=-Si;

        double tempStress[3][3]={{0,0,0},{0,0,0},{0,0,0}};
        for (int r=0 ; r<3 ; r++) 
          for (int s=0 ; s<3 ; s++) 
            for(int w=0 ; w<3 ; w++) 
              tempStress[r][s]=tempStress[r][s]+StressCellGlobal[r][w]*tempRot[w][s];
                        
        for (int r=0 ; r<3 ; r++) 
          for (int s=0 ; s<3 ; s++) 
            StressCellGlobal[r][s]=0;
          
        for (int r=0 ; r<3 ; r++) 
          for (int s=0 ; s<3 ; s++) 
            for(int w=0 ; w<3 ; w++) 
              StressCellGlobal[r][s]=StressCellGlobal[r][s]+tempRot[w][r]*tempStress[w][s];
        
        for (int r=0 ; r<3 ; r++) 
          for (int s=0 ; s<3 ; s++) 
            tempStress[r][s]=eigenVectorStress[r][s];
          
	for (size_t ii=0; ii<3; ++ii) {
	  for (size_t jj=0; jj<3; ++jj) {
	    eigenVectorStress[ii][jj] = 0.0;
	  }
	}
        
        for (int r=0 ; r<3 ; r++) 
          for (int s=0 ; s<3 ; s++) 
            for(int w=0 ; w<3 ; w++) 
              eigenVectorStress[r][s]=eigenVectorStress[r][s]+tempStress[r][w]*tempRot[w][s];
      }
      
      
      // maximal stress direction
      double maximalStressValue=StressCellGlobal[0][0];
      int Istress=0;
      if (StressCellGlobal[1][1]>maximalStressValue) 
        {
          maximalStressValue=StressCellGlobal[1][1];
          Istress=1;
        }
      if (StressCellGlobal[2][2]>maximalStressValue) 
        {
          maximalStressValue=StressCellGlobal[2][2];
          Istress=2;
        }
      // std::cerr<<"maximal Stress direction "<< eigenVectorStress[0][Istress] <<" "<< eigenVectorStress[1][Istress] <<" "<< eigenVectorStress[2][Istress] <<std::endl;  
      // std::cerr<<"maximal Stress value "<< maximalStressValue <<std::endl;  
      
    


      // storing normal dirrection to  strain in cellData  
     
      // normal to the cell plane in global direction is Zcurrent[], vector product gives the perpendicular strain direction
     
      // double PerpStrain[3];
      // PerpStrain[0]=NormalCurrent[1]*eigenVectorStrain[2][Istrain]-NormalCurrent[2]*eigenVectorStrain[1][Istrain];
      // PerpStrain[1]=NormalCurrent[2]*eigenVectorStrain[0][Istrain]-NormalCurrent[0]*eigenVectorStrain[2][Istrain];
      // PerpStrain[2]=NormalCurrent[0]*eigenVectorStrain[1][Istrain]-NormalCurrent[1]*eigenVectorStrain[0][Istrain];
      // temp=std::sqrt(PerpStrain[0]*PerpStrain[0]+PerpStrain[1]*PerpStrain[1]+PerpStrain[2]*PerpStrain[2]);     

      // if(std::abs(temp)>0.001)
      //   { // if aniso vector is not perpendicular to the cell plane
      //   if (dimension==2)
      //     { cellData[cellIndex][0]=PerpStrain[0];        
      //       cellData[cellIndex][1]=PerpStrain[1];
      //     }
      //   if (dimension==3)
      //     { cellData[cellIndex][0]=PerpStrain[0];
      //       cellData[cellIndex][1]=PerpStrain[1];
      //       cellData[cellIndex][2]=PerpStrain[2];
      //       // cellData[cellIndex][3]=10*maximalStrainValue;  //NOTE maximal Strain and Stress Values can be used this is an option
      //     }
      // }
      

      if (numVariableIndexLevel()==4 && numVariableIndex(2) ) {// storing maximal strain
        if (dimension==2)
          {
            cellData[cellIndex][variableIndex(2,0)]  =eigenVectorStrain[0][Istrain];
            cellData[cellIndex][variableIndex(2,0)+1]=eigenVectorStrain[1][Istrain];
            cellData[cellIndex][variableIndex(2,0)+2]=maximalStrainValue;  //maximal Strain Value is stored after its eigenvector
          }
        if (dimension==3)
          {
            cellData[cellIndex][variableIndex(2,0)]  =eigenVectorStrain[0][Istrain];
            cellData[cellIndex][variableIndex(2,0)+1]=eigenVectorStrain[1][Istrain];
            cellData[cellIndex][variableIndex(2,0)+2]=eigenVectorStrain[2][Istrain];
            cellData[cellIndex][variableIndex(2,0)+3]=maximalStrainValue;  //maximal Strain Value is stored after its eigenvector
          }
      }


      if (numVariableIndexLevel()==4 && numVariableIndex(3) ) { // storing maximal stress
	if (dimension==2)
	  {
	    cellData[cellIndex][variableIndex(3,0)]  =eigenVectorStress[0][Istress];
	    cellData[cellIndex][variableIndex(3,0)+1]=eigenVectorStress[1][Istress];
            cellData[cellIndex][variableIndex(3,0)+2]=maximalStressValue;  //maximal Stress Value is stored after its eigenvector
	  }
	if (dimension==3)
	  {
	    cellData[cellIndex][variableIndex(3,0)]  =eigenVectorStress[0][Istress];
	    cellData[cellIndex][variableIndex(3,0)+1]=eigenVectorStress[1][Istress];
	    cellData[cellIndex][variableIndex(3,0)+2]=eigenVectorStress[2][Istress];
            cellData[cellIndex][variableIndex(3,0)+3]=maximalStressValue;  //maximal Stress Value is stored after its eigenvector
	  }
      }
  }      
}     




void VertexFromTRBScenterTriangulation::
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
  if (variableIndex(1,0) != numVariable) {
    std::cerr << "VertexFromTRBScenterTriangulation::initiate() "
	      << "Wrong index given as start index for additional variables."
	      << std::endl;
    exit(-1);
  }
  size_t numCell = cellData.size();
  assert (numCell==T.numCell());
  std::vector<double> com(dimension);
  
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
  }
}

VertexFromTRBScenterTriangulationConcentrationHill::
VertexFromTRBScenterTriangulationConcentrationHill(std::vector<double> &paraValue, 
						   std::vector< std::vector<size_t> > 
						   &indValue ) 
{  
  // Do some checks on the parameters and variable indeces
  if( paraValue.size()!=5 ) {
    std::cerr << "VertexFromTRBScenterTriangulationConcentrationHill::"
	      << "VertexFromTRBScenterTriangulationConcentrationHill() "
	      << "Uses five parameters Young_modulus_min, Young_modulus_max, "
	      << " poisson coefficient, K_hill, and n_hill.\n";
    exit(0);
  }
  if( indValue.size()!=2 || indValue[0].size()!=2 || indValue[1].size()!=1 ) { 
    std::cerr << "VertexFromTRBScenterTriangulationConcentrationHill::"
	      << "VertexFromTRBScenterTriangulationConcentrationHill() "
	      << "Wall length and concentration indices given in first level." << std::endl
	      << "Start of additional Cell variable indices (center(x,y,z) "
	      << "L_1,...,L_n, n=num vertex) is given in second level." 
	      << std::endl;
    exit(0);
  }
  
  // Set the variable values
  setId("VertexFromTRBScenterTriangulationConcentrationHill");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  // Set the parameter identities
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "Y_mod_min";
  tmp[1] = "Y_mod_max";
  tmp[2] = "P_ratio";
  tmp[3] = "K_hill";
  tmp[4] = "n_hill";
  setParameterId( tmp );
}

void VertexFromTRBScenterTriangulationConcentrationHill::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) {
  
  //Do the update for each cell
  size_t dimension = 3;
  assert (dimension==vertexData[0].size());
  size_t numCells = T.numCell();
  size_t wallLengthIndex = variableIndex(0,0);
  size_t concIndex = variableIndex(0,1);
  size_t comIndex = variableIndex(1,0);
  size_t lengthInternalIndex = comIndex+dimension;
  double Kpow = std::pow(parameter(3),parameter(4));
  
  for (size_t cellIndex=0 ; cellIndex<numCells ; ++cellIndex) {
    size_t numWalls = T.cell(cellIndex).numWall(); 
    
    if(  T.cell(cellIndex).numVertex()!= numWalls ) {
      std::cerr << "VertexFromTRBScenterTriangulationConcentrationHill::derivs() "
		<< "same number of vertices and walls."
		<< " Not for cells with " << T.cell(cellIndex).numWall() << " walls and "
		<< T.cell(cellIndex).numVertex() << " vertices!"	
		<< std::endl;
      exit(-1);
    }
    
    double young = parameter(0) + 
      parameter(1)*Kpow/( Kpow+std::pow(cellData[cellIndex][concIndex],parameter(4)) );
    double poisson =parameter(1);
    
    // One triangle per 'vertex' in cyclic order
    for (size_t k=0; k<numWalls; ++k) { 
      size_t kPlusOneMod = (k+1)%numWalls;
      //size_t v1 = com;
      size_t v2 = T.cell(cellIndex).vertex(k)->index();
      size_t v3 = T.cell(cellIndex).vertex(kPlusOneMod)->index();
      //size_t w1 = internal k
      size_t w2 = T.cell(cellIndex).wall(k)->index();
      //size_t w3 = internal k+1
      
      // Position matrix holds in rows positions for com, vertex(k), vertex(k+1)
      DataMatrix position(3,vertexData[v2]);
      for (size_t d=0; d<dimension; ++d)
	position[0][d] = cellData[cellIndex][comIndex+d]; // com position
      //position[1] = vertexData[v2]; // given by initiation
      position[2] = vertexData[v3];
      
      // Resting lengths are from com-vertex(k), vertex(k)-vertex(k+1) (wall(k)), com-vertex(k+1)
      std::vector<double> restingLength(numWalls);
      restingLength[0] = cellData[cellIndex][lengthInternalIndex + k];
      restingLength[1] = wallData[w2][wallLengthIndex];
      restingLength[2] = cellData[cellIndex][lengthInternalIndex + kPlusOneMod];
      
      // Lengths are from com-vertex(k), vertex(k)-vertex(k+1) (wall(k)), com-vertex(k+1)
      std::vector<double> length(numWalls);
      length[0] = std::sqrt( (position[0][0]-position[1][0])*
			     (position[0][0]-position[1][0]) +
			     (position[0][1]-position[1][1])*
			     (position[0][1]-position[1][1]) +
			     (position[0][2]-position[1][2])*
			     (position[0][2]-position[1][2]) );
      length[1] = T.wall(w2).lengthFromVertexPosition(vertexData);
      length[2] = std::sqrt( (position[0][0]-position[2][0])*
			     (position[0][0]-position[2][0]) +
			     (position[0][1]-position[2][1])*
			     (position[0][1]-position[2][1]) +
			     (position[0][2]-position[2][2])*
			     (position[0][2]-position[2][2]) );
      
      // Lame coefficients (can be defined out of loop)
      double lambda=young*poisson/(1-poisson*poisson);
      double mio=young/(1+poisson);
      
      // Area of the element (using Heron's formula)                                      
      double Area=std::sqrt( ( restingLength[0]+restingLength[1]+restingLength[2])*
			     (-restingLength[0]+restingLength[1]+restingLength[2])*
			     ( restingLength[0]-restingLength[1]+restingLength[2])*
			     ( restingLength[0]+restingLength[1]-restingLength[2])  )*0.25;
      
      //Angles of the element ( assuming the order: 0,L0,1,L1,2,L2 )
      std::vector<double> Angle(3);
      // can be ommited by cotan(A)=.25*sqrt(4*b*b*c*c/K-1)
      Angle[0]=std::acos(  (restingLength[0]*restingLength[0]+restingLength[2]*restingLength[2]-restingLength[1]*restingLength[1])/
                           (restingLength[0]*restingLength[2]*2)    );
      Angle[1]=std::acos(  (restingLength[0]*restingLength[0]+restingLength[1]*restingLength[1]-restingLength[2]*restingLength[2])/
                           (restingLength[0]*restingLength[1]*2)    );
      Angle[2]=std::acos(  (restingLength[1]*restingLength[1]+restingLength[2]*restingLength[2]-restingLength[0]*restingLength[0])/
                           (restingLength[1]*restingLength[2]*2)    );
      
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
      
      // adding forces to the total vertexDerivs
      
      cellDerivs[cellIndex][comIndex]+= Force[0][0];
      cellDerivs[cellIndex][comIndex+1]+= Force[0][1];
      cellDerivs[cellIndex][comIndex+2]+= Force[0][2];
      
      vertexDerivs[v2][0]+= Force[1][0];
      vertexDerivs[v2][1]+= Force[1][1];
      vertexDerivs[v2][2]+= Force[1][2];
      
      vertexDerivs[v3][0]+= Force[2][0];
      vertexDerivs[v3][1]+= Force[2][1];
      vertexDerivs[v3][2]+= Force[2][2];
    }
  }
}

void VertexFromTRBScenterTriangulationConcentrationHill::
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
  if (variableIndex(1,0) != numVariable) {
    std::cerr << "VertexFromTRBScenterTriangulationConcentrationHill::initiate() "
	      << "Wrong index given as start index for additional variables."
	      << std::endl;
    exit(-1);
  }
  size_t numCell = cellData.size();
  assert (numCell==T.numCell());
  std::vector<double> com(dimension);
  
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
  }
}



VertexFromTRBSMT::
VertexFromTRBSMT(std::vector<double> &paraValue, 
	       std::vector< std::vector<size_t> > 
	       &indValue ) 
{  
  // Do some checks on the parameters and variable indeces
  if( paraValue.size()!=10 ) {
    std::cerr << "VertexFromTRBSMT::"
	      << "VertexFromTRBSMT() "
              << "Uses ten parameters: "
              << "0,1 : young modulus(matrix and fibre) " 
              << "2,3 : poisson ratio (longitudinal (MT) and transverse directions)"
	      << "4 : MF flag(0 constant material anisotropy ,1: material anisotropy via FiberModel "
              << "5 : K_Hill, 6: n_Hill  " 
              << "7 : flag(0: plane strain, 1: plane stress) " 
              << "8 : MT direction angle"
              << "9 : flag(0: for no feedback or direct feedback by indices,"
              << "         1: for MT direction from 6th parameter TETA, 2: force to Stress,  "
              << "         3: force to Strain , 4: force to perp-strain "<< std::endl;
    exit(0);
  }

  if( (indValue.size()!=1 && indValue.size()!=3) ||
      indValue[0].size()!=10 ||
      (indValue.size()==3 && (indValue[1].size()!=0 && indValue[1].size()!=1 && indValue[1].size()!=2 && indValue[1].size()!=3 )) ||
      (indValue.size()==3 && (indValue[2].size()!=0 && indValue[2].size()!=1 && indValue[2].size()!=2)) 
      ) { 
    std::cerr << "VertexFromTRBSMT::"
	      << "VertexFromTRBSMT() "
              << "resting length index and MT direction initial index and strain and stress anisotropy indices and " 
              << "indices for storing area ratio, isotropic energy and anisotropic energy and young_fiber and MTstress and stress tensor(6 components) given in first level." 
              << "Optionally two additional levels can be given where the strain and perpendicular direction to strain in cell plane and 2nd strain" 
              << "directions/values(dx dy dz value) can be stored at given indices at second level."
              << "If no index given at second level, strain will not be stored, if one index given strain will be stored"
              << "and if two indices are given maximal and perpendicular to strain will be stored at first and second"
              << "indices at second level respectively"
              << "and if 3 indices are given maximal and perpendicular to strain and 2nd strain will be stored at 1st, 2nd and 3rd "
              << "indices at second level respectively"
              << "If no index given at 3rd level, stress will not be stored, if one index given stress will be stored"
              << "and if two indices are given maximal and 2nd stress will be stored at first and second"
              << "indices at 3rd level respectively"
              << std::endl;
    exit(0);
  }
  
  // Set the variable values
  setId("VertexFromTRBSMT");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  

  // Set the parameter identities
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "Y_mod_M";   // Matrix Young modulus
  tmp[1] = "Y_mod_F";   // Fiber Young modulus
  tmp[2] = "P_ratio_L"; // Longitudinal Poisson ratio
  tmp[3] = "P_ratio_T"; // Transverse Poisson ratio
  tmp[4] = "MF flag";
  tmp[5] = "neighborweight";
  tmp[6] = "parameter6";
  tmp[7] = "Strain-Stress flag";
  tmp[8] = "TETA anisotropy";
  tmp[9] = "MT update flag";

  setParameterId( tmp );
 
  if( parameter(2)<0 || parameter(2)>=0.5 || parameter(3)<0 || parameter(3)>=0.5 ) {
    std::cerr << "VertexFromTRBSMT::"
 	      << "VertexFromTRBSMT() "
 	      << "poisson ratios must be 0 <= p < 0.5 " << std::endl;
    exit(0);
  }

  if( parameter(4)!=0 && parameter(4)!=1 ) {
    std::cerr << " VertexFromTRBSMT::"
	      << " VertexFromTRBSMT() "
	      << " 5th parameter must be 0 or 1 " 
              << " 0: constant material anisotropy, 1: material anisotropy via FiberModel " << std::endl;
    exit(0);
    exit(0);
  }

  // if( parameter(5)<0 || parameter(5)>=1 ) {
  //   std::cerr << " VertexFromTRBSMT::"
  //             << " VertexFromTRBSMT() "
  //             << " K parameter in Hill function must be 0 < K < 1 " << std::endl;
  //   exit(0);
  // }

  if( parameter(7)!=0 && parameter(7)!=1 ) {
    std::cerr << " VertexFromTRBSMT::"
	      << " VertexFromTRBSMT() "
	      << " 6th parameter must be 0 or 1(0:plane strain, 1:plane stress) " << std::endl;
    exit(0);
  }

  if( parameter(9)!=0 && parameter(9)!=1 && parameter(9)!=2 && parameter(9)!=3 && parameter(9)!=4) {
    std::cerr << " VertexFromTRBSMT::"
              << " VertexFromTRBSMT() "
              << " 8th parameter must be 0/1/2/3/4"
              << " 0: for no feedback or direct feedback by indices "
              << " 1: for MT direction from 7th parameter TETA "
              << " 2: force to Stress "
              << " 3: force to S1train "
              << " 4: force to perp-strain " << std::endl;
    exit(0);
  }
}

void VertexFromTRBSMT::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) {
  
  //Do the update for each cell
  size_t dimension = 3;
  assert (dimension==vertexData[0].size());
  size_t numCells = T.numCell();
  size_t numWalls = 3; // defined only for triangles 
  
  // double TotalVolume=0;
  // double deltaVolume=0;

  size_t wallLengthIndex   =variableIndex(0,0);
  size_t MTindex           =variableIndex(0,1);	 
  size_t strainAnIndex     =variableIndex(0,2);	
  size_t stressAnIndex     =variableIndex(0,3);	
  size_t areaRatioIndex    =variableIndex(0,4);	
  size_t isoEnergyIndex    =variableIndex(0,5);	
  size_t anisoEnergyIndex  =variableIndex(0,6);	
  size_t youngLIndex       =variableIndex(0,7);	
  size_t MTstressIndex     =variableIndex(0,8);	
  size_t stressTensorIndex =variableIndex(0,9);	
 


  for( size_t cellIndex=0 ; cellIndex<numCells ; ++cellIndex ) {
    if( T.cell(cellIndex).numWall() != numWalls ) {
      std::cerr << "VertexFromTRBSMT::derivs() only defined for triangular cells."
		<< " Not for cells with " << T.cell(cellIndex).numWall() << " walls!"
		<< std::endl;
      exit(-1);
    }
    
    double youngMatrix= parameter(0);    
    double youngFiber = parameter(1); 
    double poissonL   = parameter(2);    
    double poissonT   = parameter(3);
    double TETA= parameter(8); 


    size_t v1 = T.cell(cellIndex).vertex(0)->index();
    size_t v2 = T.cell(cellIndex).vertex(1)->index();
    size_t v3 = T.cell(cellIndex).vertex(2)->index();
    size_t w1 = T.cell(cellIndex).wall(0)->index();
    size_t w2 = T.cell(cellIndex).wall(1)->index();
    size_t w3 = T.cell(cellIndex).wall(2)->index();




    //std::cerr<< "cell "<< cellIndex<< " vertices  "<< v1<<" "<< v2 << " "<< v3 << " walls  "<< w1 <<" "<< w2 << " "<< w3<< std::endl;
    
   

    double youngL=1;
    double youngT=1;    

    if( parameter(4)==0 ){ // constant anisotropic material
      youngL = youngMatrix+youngFiber;
      youngT = youngMatrix; 
    }
    if( parameter(4)==1){  // material anisotropy via FiberModel
      youngL = cellData[cellIndex][youngLIndex]; 
      youngT = 2*youngMatrix+youngFiber-youngL; 
    }
    
   
    double lambdaL, mioL, lambdaT, mioT;
    
    if (parameter(7)==0){      
      // Lame coefficients based on plane strain (for 3D 0<poisson<0.5)
      lambdaL=youngL*poissonL/((1+poissonL)*(1-2*poissonL));
      mioL=youngL/(2*(1+poissonL));
      lambdaT=youngT*poissonT/((1+poissonT)*(1-2*poissonT));
      mioT=youngT/(2*(1+poissonT));
    } 
    else{      
      // Lame coefficients based on plane stress (for 3D 0<poisson<0.5)
      lambdaL=youngL*poissonL/(1-poissonL*poissonL);
      mioL=youngL/(2*(1+poissonL));
      lambdaT=youngT*poissonT/(1-poissonT*poissonT);
      mioT=youngT/(2*(1+poissonT));
    }
    
    // Lame coefficients based on delin. paper (for 2D 0<poisson<1)
    // double lambdaL=youngL*poissonL/(1-poissonL*poissonL);
    // double mioL=youngL/(1+poissonL);
    // double lambdaT=youngT*poissonT/(1-poissonT*poissonT);
    // double mioT=youngT/(1+poissonT);
     
 
    
    double EnergyIso=0;                      
    double EnergyAniso=0;
    double strainZ=0;


    double TETA1=parameter(5);
    if ( parameter(9)==1 ) {  // aniso direction from TETA 
      if ( cellIndex==0 || cellIndex==2) {
	cellData[cellIndex][MTindex ]=std::cos(TETA1);  
	cellData[cellIndex][MTindex +1]=std::sin(TETA1);
	cellData[cellIndex][MTindex +2]=0;
      }
      if ( cellIndex==1 || cellIndex==3) {
	cellData[cellIndex][MTindex ]=std::cos(TETA);  
	cellData[cellIndex][MTindex +1]=std::sin(TETA);
	cellData[cellIndex][MTindex +2]=0;
      }
    }
    // if ( parameter(9)==1 ) {  // aniso direction from TETA 
    //   cellData[cellIndex][variableIndex(0,1)]=std::cos(TETA);  
    //   cellData[cellIndex][variableIndex(0,1)+1]=std::sin(TETA);
    //   cellData[cellIndex][variableIndex(0,1)+2]=0;
    // }
    
    // Aniso vector in current shape in global coordinate system
    double AnisoCurrGlob[3];
    AnisoCurrGlob[0] = cellData[cellIndex][MTindex ];
    AnisoCurrGlob[1] = cellData[cellIndex][MTindex +1];
    AnisoCurrGlob[2] = cellData[cellIndex][MTindex +2];
    

    std::vector<double> restingLength(numWalls);
    restingLength[0] = wallData[w1][wallLengthIndex];
    restingLength[1] = wallData[w2][wallLengthIndex];
    restingLength[2] = wallData[w3][wallLengthIndex];

    DataMatrix position(3,vertexData[v1]);
    position[1] = vertexData[v2];
    position[2] = vertexData[v3];
    //position[0][2] z for vertex 1 (of the cell)
   
    std::vector<double> length(numWalls);
    length[0] = T.wall(w1).lengthFromVertexPosition(vertexData);
    length[1] = T.wall(w2).lengthFromVertexPosition(vertexData);
    length[2] = T.wall(w3).lengthFromVertexPosition(vertexData);
    
    //Anisotropic Correction is based on difference between Lam Coefficients of Longitudinal and Transverse dirrections:
    double deltaLam=lambdaL-lambdaT;
    double deltaMio=mioL-mioT;
    // double deltaMio1=(youngL-youngT)/2;
    // double deltaLam=AnisoMeasure*(lambdaL-lambdaT);
    // double deltaMio=AnisoMeasure*(mioL-mioT);
    
    //Resting area of the element (using Heron's formula)                                      
    double restingArea=std::sqrt( ( restingLength[0]+restingLength[1]+restingLength[2])*
                                  (-restingLength[0]+restingLength[1]+restingLength[2])*
                                  ( restingLength[0]-restingLength[1]+restingLength[2])*
                                  ( restingLength[0]+restingLength[1]-restingLength[2])  )*0.25;
    
    //Angles of the element ( assuming the order: 0,L0,1,L1,2,L2 clockwise )
    std::vector<double> Angle(3);
    // can be ommited by cotan(A)=.25*sqrt(4*b*b*c*c/K-1)
    Angle[0]=std::acos(  (restingLength[0]*restingLength[0]+restingLength[2]*restingLength[2]-restingLength[1]*restingLength[1])/
                         (restingLength[0]*restingLength[2]*2)    );
    Angle[1]=std::acos(  (restingLength[0]*restingLength[0]+restingLength[1]*restingLength[1]-restingLength[2]*restingLength[2])/
                         (restingLength[0]*restingLength[1]*2)    );
    Angle[2]=std::acos(  (restingLength[1]*restingLength[1]+restingLength[2]*restingLength[2]-restingLength[0]*restingLength[0])/
                         (restingLength[1]*restingLength[2]*2)    );
    
    //Tensile Stiffness
    double tensileStiffness[3];
    double temp = 1.0/(restingArea*16);                                      
    std::vector<double> cotan(3);
    cotan[0] = 1.0/std::tan(Angle[0]);
    cotan[1] = 1.0/std::tan(Angle[1]);
    cotan[2] = 1.0/std::tan(Angle[2]);    
    //the force is calculated based on Transverse coefficients
    //Longitudinal coefficients are considered in deltaF
    
    tensileStiffness[0]=(2*cotan[2]*cotan[2]*(lambdaT+2*mioT)+2*mioT)*temp;
    tensileStiffness[1]=(2*cotan[0]*cotan[0]*(lambdaT+2*mioT)+2*mioT)*temp;
    tensileStiffness[2]=(2*cotan[1]*cotan[1]*(lambdaT+2*mioT)+2*mioT)*temp;
    
    //Angular Stiffness
    std::vector<double> angularStiffness(3);
    angularStiffness[0]=(2*cotan[1]*cotan[2]*(lambdaT+2*mioT)-2*mioT)*temp;                          
    angularStiffness[1]=(2*cotan[0]*cotan[2]*(lambdaT+2*mioT)-2*mioT)*temp;
    angularStiffness[2]=(2*cotan[0]*cotan[1]*(lambdaT+2*mioT)-2*mioT)*temp;
    
    //Calculate biquadratic strains  
    std::vector<double> Delta(3);
    Delta[0]=(length[0])*(length[0])-(restingLength[0])*(restingLength[0]);
    Delta[1]=(length[1])*(length[1])-(restingLength[1])*(restingLength[1]);
    Delta[2]=(length[2])*(length[2])-(restingLength[2])*(restingLength[2]);

    //Area of the element (using Heron's formula)                                      
    double Area=std::sqrt( ( length[0]+length[1]+length[2])*
                           (-length[0]+length[1]+length[2])*
                           ( length[0]-length[1]+length[2])*
                           ( length[0]+length[1]-length[2])  )*0.25;


    // calculating the angles between shape vectors and anisotropy direction in resting shape when anisotropy vector is provided in current shape
    
    //Current shape local coordinate of the element  (counterclockwise ordering of nodes/edges)
      double CurrentAngle1=std::acos(  (length[0]*length[0]+length[1]*length[1]-length[2]*length[2])/
                                       (length[0]*length[1]*2)    );

      double Qa=std::cos(CurrentAngle1)*length[0];
      double Qc=std::sin(CurrentAngle1)*length[0];
      double Qb=length[1];
      // shape vector matrix = inverse of coordinate matrix ( only first two elements i.e. ShapeVector[3][2] )      
      // double ShapeVectorCurrent[3][3]={ {  0   ,       1/Qc      , 0 }, 
      //                                   {-1/Qb , (Qa-Qb)/(Qb*Qc) , 1 },       
      //                                   { 1/Qb ,     -Qa/(Qb*Qc) , 0 }  };
            
    
      // Local coordinates of the resting shape ( counterclockwise )
      double RestingAngle1=std::acos(  (restingLength[0]*restingLength[0]+restingLength[1]*restingLength[1]-restingLength[2]*restingLength[2])/
                                       (restingLength[0]*restingLength[1]*2)    );

      double Pa=std::cos(RestingAngle1)*restingLength[0];
      double Pc=std::sin(RestingAngle1)*restingLength[0];
      double Pb=restingLength[1];

      // shape vector matrix in resting shape in local coordinate system  = inverse of coordinate matrix ( only first two elements i.e. ShapeVectorResting[3][2] )      
      double ShapeVectorResting[3][3]={ {  0   ,       1/Pc      , 0 }, 
                                        {-1/Pb , (Pa-Pb)/(Pb*Pc) , 1 },       
                                        { 1/Pb ,     -Pa/(Pb*Pc) , 0 }  };
     
      
      // //Strain tensor  (clockwise ordering of nodes/edges)
      // double CurrentAngle2=std::acos(  (length[1]*length[1]+length[2]*length[2]-length[0]*length[0])/
      //                                  (length[1]*length[2]*2)    );

      // double Qa=std::cos(CurrentAngle2)*length[2];
      // double Qb=length[1];
      // double Qc=std::sin(CurrentAngle2)*length[2];
      
      // double ShapeVectorCurrent[3][2]={ {  0   ,       1/Qc      }, 
      //                                   { 1/Qb ,     -Qa/(Qb*Qc) },       
      //                                   {-1/Qb , (Qa-Qb)/(Qb*Qc) }  };

     
    

      // Rotation Matrix for changing coordinate systems for both Local to Global( Strain Tensor) and Global to Local( Aniso Vector in the current shape)
      double rotation[3][3];  

      double tempA=std::sqrt((position[2][0]-position[1][0])*(position[2][0]-position[1][0])+
                             (position[2][1]-position[1][1])*(position[2][1]-position[1][1])+
                             (position[2][2]-position[1][2])*(position[2][2]-position[1][2])  );

      double tempB=std::sqrt((position[0][0]-position[1][0])*(position[0][0]-position[1][0])+
                             (position[0][1]-position[1][1])*(position[0][1]-position[1][1])+
                             (position[0][2]-position[1][2])*(position[0][2]-position[1][2])  );
      
      double Xcurrent[3];      
      Xcurrent[0]= (position[2][0]-position[1][0])/tempA;
      Xcurrent[1]= (position[2][1]-position[1][1])/tempA;
      Xcurrent[2]= (position[2][2]-position[1][2])/tempA;
      
      double Bcurrent[3];      
      Bcurrent[0]= (position[0][0]-position[1][0])/tempB;
      Bcurrent[1]= (position[0][1]-position[1][1])/tempB;
      Bcurrent[2]= (position[0][2]-position[1][2])/tempB;

      double Zcurrent[3];      
      Zcurrent[0]= Xcurrent[1]*Bcurrent[2]-Xcurrent[2]*Bcurrent[1];
      Zcurrent[1]= Xcurrent[2]*Bcurrent[0]-Xcurrent[0]*Bcurrent[2];
      Zcurrent[2]= Xcurrent[0]*Bcurrent[1]-Xcurrent[1]*Bcurrent[0];
      
      tempA=std:: sqrt(Zcurrent[0]*Zcurrent[0]+Zcurrent[1]*Zcurrent[1]+Zcurrent[2]*Zcurrent[2]);
      Zcurrent[0]=Zcurrent[0]/tempA;
      Zcurrent[1]=Zcurrent[1]/tempA;
      Zcurrent[2]=Zcurrent[2]/tempA;

      double Ycurrent[3];      
      Ycurrent[0]= Zcurrent[1]*Xcurrent[2]-Zcurrent[2]*Xcurrent[1];
      Ycurrent[1]= Zcurrent[2]*Xcurrent[0]-Zcurrent[0]*Xcurrent[2];
      Ycurrent[2]= Zcurrent[0]*Xcurrent[1]-Zcurrent[1]*Xcurrent[0];


      rotation[0][0]=Xcurrent[0];
      rotation[1][0]=Xcurrent[1];
      rotation[2][0]=Xcurrent[2];

      rotation[0][1]=Ycurrent[0];
      rotation[1][1]=Ycurrent[1];
      rotation[2][1]=Ycurrent[2];

      rotation[0][2]=Zcurrent[0];
      rotation[1][2]=Zcurrent[1];
      rotation[2][2]=Zcurrent[2];      
       
      // rotating the anisotropy vector from global coordinate system to the local one in the current shape
      double AnisoCurrLocal[3];
      AnisoCurrLocal[0]=rotation[0][0]*AnisoCurrGlob[0]+rotation[1][0]*AnisoCurrGlob[1]+rotation[2][0]*AnisoCurrGlob[2];
      AnisoCurrLocal[1]=rotation[0][1]*AnisoCurrGlob[0]+rotation[1][1]*AnisoCurrGlob[1]+rotation[2][1]*AnisoCurrGlob[2];
      AnisoCurrLocal[2]=rotation[0][2]*AnisoCurrGlob[0]+rotation[1][2]*AnisoCurrGlob[1]+rotation[2][2]*AnisoCurrGlob[2];
     


      double positionLocal[3][2]={ {Qa , Qc}, 
                                   {0  , 0 },  
                                   {Qb , 0 }  };
      
      double DeformGrad[2][2]={{0,0},{0,0}}; // F= Qi x Di
      for ( int i=0 ; i<3 ; ++i ) {
        DeformGrad[0][0]=DeformGrad[0][0]+positionLocal[i][0]*ShapeVectorResting[i][0];
        DeformGrad[1][0]=DeformGrad[1][0]+positionLocal[i][1]*ShapeVectorResting[i][0];
        DeformGrad[0][1]=DeformGrad[0][1]+positionLocal[i][0]*ShapeVectorResting[i][1];
        DeformGrad[1][1]=DeformGrad[1][1]+positionLocal[i][1]*ShapeVectorResting[i][1];
      }      

      double AnisoRestLocal[3]={0,0,0};

      AnisoRestLocal[0]=DeformGrad[0][0]*AnisoCurrLocal[0]+DeformGrad[1][0]*AnisoCurrLocal[1];
      AnisoRestLocal[1]=DeformGrad[0][1]*AnisoCurrLocal[0]+DeformGrad[1][1]*AnisoCurrLocal[1];
     
     
      double tempAn=std::sqrt(AnisoRestLocal[0]*AnisoRestLocal[0]+AnisoRestLocal[1]*AnisoRestLocal[1]+AnisoRestLocal[2]*AnisoRestLocal[2]);
      AnisoRestLocal[0]/=tempAn;
      AnisoRestLocal[1]/=tempAn;
      


    



      double AnisoMeasure=std::sqrt(AnisoRestLocal[0]*AnisoRestLocal[0]+AnisoRestLocal[1]*AnisoRestLocal[1]);
      // std::cerr<< "cell "<< cellIndex<<" AnisoMeasure "<< AnisoMeasure<<std::endl;
      
      // choosing a random normalized dirrection for anisotropy if AnisoVector is close to perpendicular to the cell plane
      if ( AnisoMeasure<0.0001) {
        double randomAngle=((rand()%360)*2*3.14159265)/360;
        
        AnisoRestLocal[0]=std::cos(randomAngle);
        AnisoRestLocal[1]=std::sin(randomAngle);
      }  
      else {// normalizing the anisoVector if it is not random
        AnisoRestLocal[0]=AnisoRestLocal[0]/AnisoMeasure;
        AnisoRestLocal[1]=AnisoRestLocal[1]/AnisoMeasure;        
      }
      
      
      
      //Angles between anisotropy vector and shape vectors for calculating the terms like a.Di , teta(k) = acos((dot(Anisorest,Dk))/(norm(Anisorest)*norm(Dk))),
      std::vector<double> teta(3);
      teta[0] = std::acos(  (ShapeVectorResting[0][0]*AnisoRestLocal[0]+ShapeVectorResting[0][1]*AnisoRestLocal[1])/
                            std::sqrt(ShapeVectorResting[0][0]*ShapeVectorResting[0][0]+ShapeVectorResting[0][1]*ShapeVectorResting[0][1]+0.0000001) );
      
      teta[1] = std::acos(  (ShapeVectorResting[1][0]*AnisoRestLocal[0]+ShapeVectorResting[1][1]*AnisoRestLocal[1])/
                            std::sqrt(ShapeVectorResting[1][0]*ShapeVectorResting[1][0]+ShapeVectorResting[1][1]*ShapeVectorResting[1][1]+0.0000001) );
      
      teta[2] = std::acos(  (ShapeVectorResting[2][0]*AnisoRestLocal[0]+ShapeVectorResting[2][1]*AnisoRestLocal[1])/
                            std::sqrt(ShapeVectorResting[2][0]*ShapeVectorResting[2][0]+ShapeVectorResting[2][1]*ShapeVectorResting[2][1]+0.0000001) );

      //  std::cerr<< "cell "<< cellIndex<<"  numerator  " << (ShapeVectorResting[2][0]*AnisoRestLocal[0]+ShapeVectorResting[2][1]*AnisoRestLocal[1])/
      //                      std::sqrt(ShapeVectorResting[2][0]*ShapeVectorResting[2][0]+ShapeVectorResting[2][1]*ShapeVectorResting[2][1])  << std::endl;    
    



      
      // >>>>>>>>>>>>>>>>>>>>>> STRAIN and STRESS TENSOR (BEGIN) <<<<<<<<<<<<<<<<<<<<
      
      // deformation gradiant tensor F =Sigma i=1,2,3 Qi x Di
      // strain tensor in resting shape E=0.5(FtF-I)
      // trE
      // B=FFt
      // axa (direct product of aniso vector in resting shape)
      // atEa
      // E(axa) and (axa)E
      double trE=( Delta[1]*cotan[0]+ Delta[2]*cotan[1]+Delta[0]*cotan[2])/(4*restingArea);
      
      double directAniso[2][2]={{AnisoRestLocal[0]*AnisoRestLocal[0],AnisoRestLocal[0]*AnisoRestLocal[1]},
                                {AnisoRestLocal[1]*AnisoRestLocal[0],AnisoRestLocal[1]*AnisoRestLocal[1]}};
      

      double LeftCauchy[2][2]; // B=FFt
      LeftCauchy[0][0]=DeformGrad[0][0]*DeformGrad[0][0]+DeformGrad[0][1]*DeformGrad[0][1];
      LeftCauchy[1][0]=DeformGrad[1][0]*DeformGrad[0][0]+DeformGrad[1][1]*DeformGrad[0][1];
      LeftCauchy[0][1]=DeformGrad[0][0]*DeformGrad[1][0]+DeformGrad[0][1]*DeformGrad[1][1];
      LeftCauchy[1][1]=DeformGrad[1][0]*DeformGrad[1][0]+DeformGrad[1][1]*DeformGrad[1][1];


      double Egreen[2][2];//E=0.5(C-I)
      Egreen[0][0]=0.5*(DeformGrad[0][0]*DeformGrad[0][0]+DeformGrad[1][0]*DeformGrad[1][0]-1);
      Egreen[1][0]=0.5*(DeformGrad[0][1]*DeformGrad[0][0]+DeformGrad[1][1]*DeformGrad[1][0]);
      Egreen[0][1]=0.5*(DeformGrad[0][0]*DeformGrad[0][1]+DeformGrad[1][0]*DeformGrad[1][1]);
      Egreen[1][1]=0.5*(DeformGrad[0][1]*DeformGrad[0][1]+DeformGrad[1][1]*DeformGrad[1][1]-1);


      double E2[2][2]; // used for energy calculation only
      E2[0][0]=Egreen[0][0]*Egreen[0][0]+Egreen[0][1]*Egreen[1][0];
      E2[1][0]=Egreen[1][0]*Egreen[0][0]+Egreen[1][1]*Egreen[1][0];
      E2[0][1]=Egreen[0][0]*Egreen[0][1]+Egreen[0][1]*Egreen[1][1];
      E2[1][1]=Egreen[1][0]*Egreen[0][1]+Egreen[1][1]*Egreen[1][1];
      
      double I2=E2[0][0]+E2[1][1]; //trE2 used for energy calculation only
      
       //atE2a used for energy calculation only
      double I5= AnisoRestLocal[0]*AnisoRestLocal[0]*E2[0][0]   +  
        AnisoRestLocal[0]*AnisoRestLocal[1]*(E2[0][1]+E2[1][0]) +
        AnisoRestLocal[1]*AnisoRestLocal[1]*E2[1][1];
      
      
      double StrainAlmansi[2][2]; // e=0.5(1-B^-1)  True strain tensor
      temp=LeftCauchy[0][0]*LeftCauchy[1][1]-LeftCauchy[1][0]*LeftCauchy[0][1]; // det(B)
      StrainAlmansi[0][0]=0.5*(1-(LeftCauchy[1][1]/temp));
      StrainAlmansi[1][0]=0.5*LeftCauchy[1][0]/temp;
      StrainAlmansi[0][1]=0.5*LeftCauchy[0][1]/temp;  
      StrainAlmansi[1][1]=0.5*(1-(LeftCauchy[0][0]/temp));
      
      double atEa=AnisoRestLocal[0]*AnisoRestLocal[0]*Egreen[0][0]
                 +AnisoRestLocal[0]*AnisoRestLocal[1]*(Egreen[0][1]+Egreen[1][0])
                 +AnisoRestLocal[1]*AnisoRestLocal[1]*Egreen[1][1];

      double I4=atEa;      

      double Eaa[2][2];
      Eaa[0][0]= Egreen[0][0]*directAniso[0][0]+Egreen[0][1]*directAniso[1][0];        
      Eaa[1][0]= Egreen[1][0]*directAniso[0][0]+Egreen[1][1]*directAniso[1][0];        
      Eaa[0][1]= Egreen[0][0]*directAniso[0][1]+Egreen[0][1]*directAniso[1][1];        
      Eaa[1][1]= Egreen[1][0]*directAniso[0][1]+Egreen[1][1]*directAniso[1][1];        

      double aaE[2][2];
      aaE[0][0]= directAniso[0][0]*Egreen[0][0]+directAniso[0][1]*Egreen[1][0];        
      aaE[1][0]= directAniso[1][0]*Egreen[0][0]+directAniso[1][1]*Egreen[1][0];        
      aaE[0][1]= directAniso[0][0]*Egreen[0][1]+directAniso[0][1]*Egreen[1][1];        
      aaE[1][1]= directAniso[1][0]*Egreen[0][1]+directAniso[1][1]*Egreen[1][1];        
      
      double B2[2][2];// LeftCauchy^2
      B2[0][0]=LeftCauchy[0][0]*LeftCauchy[0][0]+LeftCauchy[0][1]*LeftCauchy[1][0];
      B2[1][0]=LeftCauchy[1][0]*LeftCauchy[0][0]+LeftCauchy[1][1]*LeftCauchy[1][0];
      B2[0][1]=LeftCauchy[0][0]*LeftCauchy[0][1]+LeftCauchy[0][1]*LeftCauchy[1][1];
      B2[1][1]=LeftCauchy[1][0]*LeftCauchy[0][1]+LeftCauchy[1][1]*LeftCauchy[1][1];


      double areaFactor=restingArea/Area; // 1/detF
      //double areaFactor=restingArea/Area; // detF

      double Sigma[2][2]; // true stress tensor (isotropic term) based on lambdaT and mioT
      Sigma[0][0]=areaFactor*((lambdaT*trE-mioT)*LeftCauchy[0][0]+(mioT)*B2[0][0]);
      Sigma[1][0]=areaFactor*((lambdaT*trE-mioT)*LeftCauchy[1][0]+(mioT)*B2[1][0]);
      Sigma[0][1]=areaFactor*((lambdaT*trE-mioT)*LeftCauchy[0][1]+(mioT)*B2[0][1]);
      Sigma[1][1]=areaFactor*((lambdaT*trE-mioT)*LeftCauchy[1][1]+(mioT)*B2[1][1]);
     
      // double deltaS[2][2]; // based on Delin. paper
      // deltaS[0][0]=deltaLam*(trE*directAniso[0][0]+atEa)+(2*deltaMio)*(Eaa[0][0]+aaE[0][0])-(deltaLam+2*deltaMio)*atEa*directAniso[0][0];
      // deltaS[1][0]=deltaLam*(trE*directAniso[1][0]     )+(2*deltaMio)*(Eaa[1][0]+aaE[1][0])-(deltaLam+2*deltaMio)*atEa*directAniso[1][0];
      // deltaS[0][1]=deltaLam*(trE*directAniso[0][1]     )+(2*deltaMio)*(Eaa[0][1]+aaE[0][1])-(deltaLam+2*deltaMio)*atEa*directAniso[0][1];
      // deltaS[1][1]=deltaLam*(trE*directAniso[1][1]+atEa)+(2*deltaMio)*(Eaa[1][1]+aaE[1][1])-(deltaLam+2*deltaMio)*atEa*directAniso[1][1];

     
      double deltaS[2][2]; // based on  equipartition energy 
      deltaS[0][0]=(deltaLam/2)*(trE*directAniso[0][0]+atEa)+(deltaMio)*(Eaa[0][0]+aaE[0][0]);
      deltaS[1][0]=(deltaLam/2)*(trE*directAniso[1][0]     )+(deltaMio)*(Eaa[1][0]+aaE[1][0]);
      deltaS[0][1]=(deltaLam/2)*(trE*directAniso[0][1]     )+(deltaMio)*(Eaa[0][1]+aaE[0][1]);
      deltaS[1][1]=(deltaLam/2)*(trE*directAniso[1][1]+atEa)+(deltaMio)*(Eaa[1][1]+aaE[1][1]);


      // double deltaS[2][2]; // based on  ... 
      // deltaS[0][0]=(deltaLam/2)*(trE*directAniso[0][0]+atEa)+2*(deltaMio)* atEa * directAniso[0][0];
      // deltaS[1][0]=(deltaLam/2)*(trE*directAniso[1][0]     )+2*(deltaMio)* atEa * directAniso[1][0];
      // deltaS[0][1]=(deltaLam/2)*(trE*directAniso[0][1]     )+2*(deltaMio)* atEa * directAniso[0][1];
      // deltaS[1][1]=(deltaLam/2)*(trE*directAniso[1][1]+atEa)+2*(deltaMio)* atEa * directAniso[1][1];




      // double deltaS[2][2]; // based on  ... 
      // deltaS[0][0]=2*(deltaMio1)* atEa * directAniso[0][0];
      // deltaS[1][0]=2*(deltaMio1)* atEa * directAniso[1][0];
      // deltaS[0][1]=2*(deltaMio1)* atEa * directAniso[0][1];
      // deltaS[1][1]=2*(deltaMio1)* atEa * directAniso[1][1];


        
      strainZ =1-poissonT*((2*lambdaT*trE+2*mioT*trE)+deltaS[0][0]+deltaS[1][1])/youngT;

      //<<<<<<<<<<<<<<<<<<<isotropic force from stress tensor <<<<<<<<<<<<<<<<<<<<<<<<<<
      // double ss[2][2];//lambda(trE)I+2mioE
      // ss[0][0]=lambdaT*trE+2*mioT*Egreen[0][0];
      // ss[0][1]=            2*mioT*Egreen[0][1];
      // ss[1][0]=            2*mioT*Egreen[1][0];
      // ss[1][1]=lambdaT*trE+2*mioT*Egreen[1][1];

       

      // double TPK[2][2];// 2nd Piola Kirchhoff stress tensor 
      // TPK[0][0]=restingArea*(DeformGrad[0][0]*ss[0][0]+DeformGrad[0][1]*ss[1][0]);
      // TPK[1][0]=restingArea*(DeformGrad[1][0]*ss[0][0]+DeformGrad[1][1]*ss[1][0]);
      // TPK[0][1]=restingArea*(DeformGrad[0][0]*ss[0][1]+DeformGrad[0][1]*ss[1][1]);
      // TPK[1][1]=restingArea*(DeformGrad[1][0]*ss[0][1]+DeformGrad[1][1]*ss[1][1]);

      // //deltaFTPKlocal[i][0]= TPK[0][0]*ShapeVectorResting[i][0]+TPK[0][1]*ShapeVectorResting[i][1];
      // //deltaFTPKlocal[i][1]= TPK[1][0]*ShapeVectorResting[i][0]+TPK[1][1]*ShapeVectorResting[i][1];
     
      // double deltaFTPKlocal[2][2];
      // deltaFTPKlocal[0][0]= TPK[0][0]*ShapeVectorResting[0][0]+TPK[0][1]*ShapeVectorResting[0][1];
      // deltaFTPKlocal[0][1]= TPK[1][0]*ShapeVectorResting[0][0]+TPK[1][1]*ShapeVectorResting[0][1];
     
      // double deltaFTPK[2][2]; 
      // deltaFTPK[0][0]= rotation[0][0]*deltaFTPKlocal[0][0]+rotation[0][1]*deltaFTPKlocal[0][1];
      // deltaFTPK[0][1]= rotation[1][0]*deltaFTPKlocal[0][0]+rotation[1][1]*deltaFTPKlocal[0][1];
      // deltaFTPK[0][2]= rotation[2][0]*deltaFTPKlocal[0][0]+rotation[2][1]*deltaFTPKlocal[0][1];
      //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      // //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      double TPK[2][2];// 2nd Piola Kirchhoff stress tensor 
      TPK[0][0]=restingArea*(DeformGrad[0][0]*deltaS[0][0]+DeformGrad[0][1]*deltaS[1][0]);
      TPK[1][0]=restingArea*(DeformGrad[1][0]*deltaS[0][0]+DeformGrad[1][1]*deltaS[1][0]);
      TPK[0][1]=restingArea*(DeformGrad[0][0]*deltaS[0][1]+DeformGrad[0][1]*deltaS[1][1]);
      TPK[1][1]=restingArea*(DeformGrad[1][0]*deltaS[0][1]+DeformGrad[1][1]*deltaS[1][1]);

      //deltaFTPKlocal[i][0]= TPK[0][0]*ShapeVectorResting[i][0]+TPK[0][1]*ShapeVectorResting[i][1];
      //deltaFTPKlocal[i][1]= TPK[1][0]*ShapeVectorResting[i][0]+TPK[1][1]*ShapeVectorResting[i][1];
     
      double deltaFTPKlocal[3][2];
      deltaFTPKlocal[0][0]= TPK[0][0]*ShapeVectorResting[0][0]+TPK[0][1]*ShapeVectorResting[0][1];
      deltaFTPKlocal[0][1]= TPK[1][0]*ShapeVectorResting[0][0]+TPK[1][1]*ShapeVectorResting[0][1];
     
      deltaFTPKlocal[1][0]= TPK[0][0]*ShapeVectorResting[1][0]+TPK[0][1]*ShapeVectorResting[1][1];
      deltaFTPKlocal[1][1]= TPK[1][0]*ShapeVectorResting[1][0]+TPK[1][1]*ShapeVectorResting[1][1];

      deltaFTPKlocal[2][0]= TPK[0][0]*ShapeVectorResting[2][0]+TPK[0][1]*ShapeVectorResting[2][1];
      deltaFTPKlocal[2][1]= TPK[1][0]*ShapeVectorResting[2][0]+TPK[1][1]*ShapeVectorResting[2][1];


      double deltaFTPK[3][3]; 
      deltaFTPK[0][0]= rotation[0][0]*deltaFTPKlocal[0][0]+rotation[0][1]*deltaFTPKlocal[0][1];
      deltaFTPK[0][1]= rotation[1][0]*deltaFTPKlocal[0][0]+rotation[1][1]*deltaFTPKlocal[0][1];
      deltaFTPK[0][2]= rotation[2][0]*deltaFTPKlocal[0][0]+rotation[2][1]*deltaFTPKlocal[0][1];

      deltaFTPK[1][0]= rotation[0][0]*deltaFTPKlocal[1][0]+rotation[0][1]*deltaFTPKlocal[1][1];
      deltaFTPK[1][1]= rotation[1][0]*deltaFTPKlocal[1][0]+rotation[1][1]*deltaFTPKlocal[1][1];
      deltaFTPK[1][2]= rotation[2][0]*deltaFTPKlocal[1][0]+rotation[2][1]*deltaFTPKlocal[1][1];
      
      deltaFTPK[2][0]= rotation[0][0]*deltaFTPKlocal[2][0]+rotation[0][1]*deltaFTPKlocal[2][1];
      deltaFTPK[2][1]= rotation[1][0]*deltaFTPKlocal[2][0]+rotation[1][1]*deltaFTPKlocal[2][1];
      deltaFTPK[2][2]= rotation[2][0]*deltaFTPKlocal[2][0]+rotation[2][1]*deltaFTPKlocal[2][1];


      // //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


      double deltaSFt[2][2];
      deltaSFt[0][0]=deltaS[0][0]*DeformGrad[0][0]+deltaS[0][1]*DeformGrad[0][1];
      deltaSFt[1][0]=deltaS[1][0]*DeformGrad[0][0]+deltaS[1][1]*DeformGrad[0][1];
      deltaSFt[0][1]=deltaS[0][0]*DeformGrad[1][0]+deltaS[0][1]*DeformGrad[1][1];
      deltaSFt[1][1]=deltaS[1][0]*DeformGrad[1][0]+deltaS[1][1]*DeformGrad[1][1];
      
      double deltaSigma[2][2];// true stress tensor (anisotropic correction term)deltaLambda and deltaMio (Longitudinal-Transverse)
      deltaSigma[0][0]=areaFactor*(DeformGrad[0][0]*deltaSFt[0][0]+DeformGrad[0][1]*deltaSFt[1][0]);
      deltaSigma[1][0]=areaFactor*(DeformGrad[1][0]*deltaSFt[0][0]+DeformGrad[1][1]*deltaSFt[1][0]);
      deltaSigma[0][1]=areaFactor*(DeformGrad[0][0]*deltaSFt[0][1]+DeformGrad[0][1]*deltaSFt[1][1]);
      deltaSigma[1][1]=areaFactor*(DeformGrad[1][0]*deltaSFt[0][1]+DeformGrad[1][1]*deltaSFt[1][1]);
    

      //double anisoimpact=(2*( youngL - youngMatrix)/youngFiber)-1;
      //std::cerr <<"aniso impact"<< anisoimpact << std::endl;

      // double sigratio=0.2; 
      //double sigmafactor=1;
      // if( anisoimpact <=sigratio)   sigmafactor=(sigratio-anisoimpact)/(1-anisoimpact);
      // // cellData[cellIndex][stressAnIndex ] :  total stress aniso
      // else  sigmafactor=0;
     
      double StressTensor[3][3];
      StressTensor[0][0]=Sigma[0][0]+deltaSigma[0][0];
      StressTensor[1][0]=Sigma[1][0]+deltaSigma[1][0];
      StressTensor[0][1]=Sigma[0][1]+deltaSigma[0][1];
      StressTensor[1][1]=Sigma[1][1]+deltaSigma[1][1];
     
 
      //Shape vectors in Current shape (counterclockwise ordering of nodes/edges)     ShapeVectorCurrent[3][3]  calculated above   
      
          
 
      //square of radius of circumstancing circle in resting shape
      //double Rcirc2Resting=(0.25*restingLength[0]*restingLength[1]*restingLength[2]/Area)*(0.25*restingLength[0]*restingLength[1]*restingLength[2]/Area);  
      
      double StrainTensor[3][3];
      //1
      StrainTensor[0][0]=StrainAlmansi[0][0];
      StrainTensor[1][0]=StrainAlmansi[1][0];
      StrainTensor[0][1]=StrainAlmansi[0][1];
      StrainTensor[1][1]=StrainAlmansi[1][1];


      StrainTensor[0][2]=0;  // adding 3rd dimension which is zero, the tensor is still in element plane
      StrainTensor[1][2]=0;
      StrainTensor[2][2]=0;
      StrainTensor[2][0]=0;
      StrainTensor[2][1]=0;
      
      StressTensor[0][2]=0;  // adding 3rd dimension which is zero, the tensor is still in element plane
      StressTensor[1][2]=0;
      StressTensor[2][2]=0;
      StressTensor[2][0]=0;
      StressTensor[2][1]=0;
 

      

      //rotation matrix to go to global coordinate system based on counterclockwise ordering;   rotation[3][3] calculated above  
      //double testRatio1=Sigma[0][0]*
      // double testRatio1=StressTensor[0][0]*(1-2*StrainTensor[0][0])*(1-2*StrainTensor[0][0])/((mioT*StrainTensor[0][0]+lambdaT*trE*(1-2*StrainTensor[0][0]))*(Area/restingArea));
      // std::cerr<< "test1 "<<testRatio1<<std::endl;
      // double testRatio2=StressTensor[1][1]*(1-2*StrainTensor[1][1])*(1-2*StrainTensor[1][1])/((mioT*StrainTensor[1][1]+lambdaT*trE*(1-2*StrainTensor[1][1]))*(Area/restingArea));
      // std::cerr<< "test2 "<<testRatio2<<std::endl;
      // std::cerr<< "trE1 "<<trE<<std::endl;
      // std::cerr<< "trE2 "<<Egreen[0][0]+Egreen[1][1]<<std::endl;
      // std::cerr<< "mioT "<<mioT<<std::endl;
      // std::cerr<< "Area "<<Area<<std::endl;
      //std::cerr<< "restingArea "<<restingArea<<std::endl;

      //   std::cerr<< "alfa0           "<<StrainAlmansi[0][0]<<std::endl;

      // std::cerr<< "0.5(1-(1/beta0)) "<<0.5*(1-(1/LeftCauchy[0][0]))<<std::endl;
      // std::cerr<< "0.5(1-(1/beta1)) "<<0.5*(1-(1/LeftCauchy[1][1]))<<std::endl;
      // std::cerr<< "almansi 01           "<<StrainAlmansi[0][1]<<std::endl;
      // std::cerr<< "almansi 10           "<<StrainAlmansi[1][0]<<std::endl;
      
      // std::cerr <<"rotation tensor " << std::endl;
      // std::cerr <<" Rxx  "<< rotation[0][0] <<" Rxy  "<< rotation[0][1] <<" Rxz  "<< rotation[0][2] << std::endl
      //           <<" Ryx  "<< rotation[1][0] <<" Ryy  "<< rotation[1][1] <<" Ryz  "<< rotation[1][2] << std::endl
      //           <<" Rzx  "<< rotation[2][0] <<" Rzy  "<< rotation[2][1] <<" Rzz  "<< rotation[2][2] << std::endl <<std::endl;
      //   std::cerr<< "deltaSigma[0][0] and [1][1]  "<<deltaSigma[0][0]<<"  "<<deltaSigma[1][1]<<std::endl;
      // std::cerr<< "Sigma[0][0] and [1][1]  "<<Sigma[0][0]<<"  "<<Sigma[1][1]<<std::endl;
      
      // double test1=((mioT*StrainTensor[0][0]+lambdaT*trE*(1-2*StrainTensor[0][0]))*(Area/restingArea))/((1-2*StrainTensor[0][0])*(1-2*StrainTensor[0][0]));
      // double test2=((mioT*StrainTensor[1][1]+lambdaT*trE*(1-2*StrainTensor[1][1]))*(Area/restingArea))/((1-2*StrainTensor[1][1])*(1-2*StrainTensor[1][1]));
      // std::cerr<< "test1  test2  "<<test1<<" "<<test2<<std::endl;
      
 
         // std::cerr<< "Egreen[0][0] and Egreen[1][1]  "<<Egreen[0][0]<<"  "<<Egreen[1][1]<<std::endl;
         // std::cerr<< "cauchyStress[0][0] and cauchyStressEgreen[1][1]  "<<lambdaT*trE+mioT*Egreen[0][0]<<"  "<<lambdaT*trE+mioT*Egreen[1][1]<<std::endl;
      
 
      // rotating strain tensor to the global coordinate system
      double tempR[3][3]={{0,0,0},{0,0,0},{0,0,0}};
      for (int r=0 ; r<3 ; r++) {
        for (int s=0 ; s<3 ; s++) {
          for(int w=0 ; w<3 ; w++) {
            tempR[r][s]=tempR[r][s]+rotation[r][w]*StrainTensor[w][s];
          }
        }
      }
      for (int r=0 ; r<3 ; r++) {
        for (int s=0 ; s<3 ; s++) {
          StrainTensor[r][s]=0;
          for(int w=0 ; w<3 ; w++) {
            StrainTensor[r][s]=StrainTensor[r][s]+tempR[r][w]*rotation[s][w];
          }
        }
      }
      
  
      // rotating stress tensor to the global coordinate system
  
      for (int r=0 ; r<3 ; r++) {
        for (int s=0 ; s<3 ; s++) {
          tempR[r][s]=0;
        }
      }
 
      for (int r=0 ; r<3 ; r++) {
        for (int s=0 ; s<3 ; s++) {
          for(int w=0 ; w<3 ; w++) {
            tempR[r][s]=tempR[r][s]+rotation[r][w]*StressTensor[w][s];
          }
        }
      }
      for (int r=0 ; r<3 ; r++) {
        for (int s=0 ; s<3 ; s++) {
          StressTensor[r][s]=0;
          for(int w=0 ; w<3 ; w++) {
            StressTensor[r][s]=StressTensor[r][s]+tempR[r][w]*rotation[s][w];
          }
        }
      }

      

      //Cell *  cell1=&(T.cell(cellIndex));
      // std::cerr<<" cell   "<<cellIndex<<"   wall   "<<w1<<  "  cell neighbohr  "<<(cell1->cellNeighbor(0))->index() <<std::endl;
      // std::cerr<<" cell   "<<cellIndex<<"   wall   "<<w2<<  "  cell neighbohr  "<<(cell1->cellNeighbor(1))->index() <<std::endl;
      // std::cerr<<" cell   "<<cellIndex<<"   wall   "<<w3<<  "  cell neighbohr  "<<(cell1->cellNeighbor(2))->index() <<std::endl;
      cellData[cellIndex][stressTensorIndex]=StressTensor[0][0];
      cellData[cellIndex][stressTensorIndex+1]=StressTensor[1][1];
      cellData[cellIndex][stressTensorIndex+2]=StressTensor[2][2];
      cellData[cellIndex][stressTensorIndex+3]=StressTensor[0][1];
      cellData[cellIndex][stressTensorIndex+4]=StressTensor[0][2];
      cellData[cellIndex][stressTensorIndex+5]=StressTensor[1][2];
      
      
      
      
      // stress component along MT direction 
      temp = std::sqrt ( cellData[cellIndex][MTindex   ] * cellData[cellIndex][MTindex   ] +
			 cellData[cellIndex][MTindex +1] * cellData[cellIndex][MTindex +1] +
			 cellData[cellIndex][MTindex +2] * cellData[cellIndex][MTindex +2]  );
      cellData[cellIndex][MTindex   ] /=temp;
      cellData[cellIndex][MTindex +1] /=temp;
      cellData[cellIndex][MTindex +2] /=temp;
      
      cellData[cellIndex][MTstressIndex] =
	cellData[cellIndex][MTindex   ] *cellData[cellIndex][MTindex   ] *StressTensor[0][0]  +
	cellData[cellIndex][MTindex   ] *cellData[cellIndex][MTindex +1] *StressTensor[0][1]  +
	cellData[cellIndex][MTindex   ] *cellData[cellIndex][MTindex +2] *StressTensor[0][2]  +
	cellData[cellIndex][MTindex +1] *cellData[cellIndex][MTindex   ] *StressTensor[1][0]  +
	cellData[cellIndex][MTindex +1] *cellData[cellIndex][MTindex +1] *StressTensor[1][1]  +
	cellData[cellIndex][MTindex +1] *cellData[cellIndex][MTindex +2] *StressTensor[1][2]  +
	cellData[cellIndex][MTindex +2] *cellData[cellIndex][MTindex   ] *StressTensor[2][0]  +
	cellData[cellIndex][MTindex +2] *cellData[cellIndex][MTindex +1] *StressTensor[2][1]  +
	cellData[cellIndex][MTindex +2] *cellData[cellIndex][MTindex +2] *StressTensor[2][2]       ;
  
     
      // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> STRAIN and STRESS TENSORS (END) <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      
      // eigenvalues/eigenvectors of strain tensor in global coordinate system. (Jacobi method)
      
      double strainEpcilon =0.000001;
      

      double eigenVectorStrain[3][3]={{1,0,0},{0,1,0},{0,0,1}};
      double pivot=1;
      double pi=3.14159265;
      int I,J;
      double RotAngle,Si,Co;
      
      pivot=std::abs(StrainTensor[1][0]);
      I=1;
      J=0;
      if (std::abs(StrainTensor[2][0])>pivot) {
        pivot=std::abs(StrainTensor[2][0]);
        I=2;
        J=0;
      }
      if (std::abs(StrainTensor[2][1])>pivot) {
        pivot=std::abs(StrainTensor[2][1]);
        I=2;
        J=1;
      }

      while (pivot>strainEpcilon) {
      
        if (std::abs(StrainTensor[I][I]-StrainTensor[J][J])<strainEpcilon) {
          RotAngle=pi/4;
        }            
        else {
          RotAngle=0.5*std::atan((2*StrainTensor[I][J])/(StrainTensor[J][J]-StrainTensor[I][I]));
        }
        Si=std::sin(RotAngle);
        Co=std::cos(RotAngle);
        double tempRot[3][3]={{1,0,0},{0,1,0},{0,0,1}};
        tempRot[I][I]=Co;
        tempRot[J][J]=Co;
        tempRot[I][J]=Si;
        tempRot[J][I]=-Si;
        double tempStrain[3][3]={{0,0,0},{0,0,0},{0,0,0}};
        for (int r=0 ; r<3 ; r++) {
          for (int s=0 ; s<3 ; s++) {
            for(int w=0 ; w<3 ; w++) {
              tempStrain[r][s]=tempStrain[r][s]+StrainTensor[r][w]*tempRot[w][s];
            }
          }
        }
        
        for (int r=0 ; r<3 ; r++) {
          for (int s=0 ; s<3 ; s++) {
            StrainTensor[r][s]=0;
          }
        }
        
        for (int r=0 ; r<3 ; r++) {
          for (int s=0 ; s<3 ; s++) {
            for(int w=0 ; w<3 ; w++) {
              StrainTensor[r][s]=StrainTensor[r][s]+tempRot[w][r]*tempStrain[w][s];
            }
          }
        }
        
        for (int r=0 ; r<3 ; r++) 
          {
            for (int s=0 ; s<3 ; s++) 
              {
                tempStrain[r][s]=eigenVectorStrain[r][s];
              }
          }

	for (size_t ii=0; ii<3; ++ii) {
	  for (size_t jj=0; jj<3; ++jj) {
	    eigenVectorStrain[ii][jj] = 0.0; 
	  }
	}	
        for (int r=0 ; r<3 ; r++) {
          for (int s=0 ; s<3 ; s++) {
            for(int w=0 ; w<3 ; w++) {
              eigenVectorStrain[r][s]=eigenVectorStrain[r][s]+tempStrain[r][w]*tempRot[w][s];
            }
          }
        }
        pivot=std::abs(StrainTensor[1][0]);
        I=1;
        J=0;
        if (std::abs(StrainTensor[2][0])>pivot) {
          pivot=std::abs(StrainTensor[2][0]);
          I=2;
          J=0;
        }
        if (std::abs(StrainTensor[2][1])>pivot) {
          pivot=std::abs(StrainTensor[2][1]);
          I=2;
          J=1;
        }
      }
      
      
      // maximal strain direction
      double maximalStrainValue=StrainTensor[0][0];
      int Istrain=0;
      if (std::abs(StrainTensor[1][1])>std::abs(maximalStrainValue)) 
        {
          maximalStrainValue=StrainTensor[1][1];
          Istrain=1;
        }
      if (std::abs(StrainTensor[2][2])>std::abs(maximalStrainValue)) 
        {
          maximalStrainValue=StrainTensor[2][2];
          Istrain=2;
        }
      //std::cerr<<"maximal Strain direction "<< eigenVectorStrain[0][Istrain] <<" "<< eigenVectorStrain[1][Istrain] <<" "<< eigenVectorStrain[2][Istrain] <<std::endl;  
      //std::cerr<<"maximal Strain value "<< maximalStrainValue <<std::endl;  
      
      // 2nd maximalstrain direction/value
      double maximalStrainValue2;
      int Istrain2,Istrain3;
      if (Istrain==0) {
        Istrain2=1;
        Istrain3=2;
      }
      if (Istrain==1) {
        Istrain2=0;
        Istrain3=2;
      }
      if (Istrain==2) {
        Istrain2=0;
        Istrain3=1;
      }
      if(std::abs(StrainTensor[Istrain3][Istrain3])>std::abs(StrainTensor[Istrain2][Istrain2])) {
        Istrain2=Istrain3;
      }
      maximalStrainValue2=StrainTensor[Istrain2][Istrain2];
      
      
      

   
      //perpendicular direction to strain in cellData
      
      // normal to the cell plane in global direction is Zcurrent[], vector product gives the perpendicular strain direction
      double PerpStrain[3];
      PerpStrain[0]=Zcurrent[1]*eigenVectorStrain[2][Istrain]-Zcurrent[2]*eigenVectorStrain[1][Istrain];
      PerpStrain[1]=Zcurrent[2]*eigenVectorStrain[0][Istrain]-Zcurrent[0]*eigenVectorStrain[2][Istrain];
      PerpStrain[2]=Zcurrent[0]*eigenVectorStrain[1][Istrain]-Zcurrent[1]*eigenVectorStrain[0][Istrain];
      temp=std::sqrt(PerpStrain[0]*PerpStrain[0]+PerpStrain[1]*PerpStrain[1]+PerpStrain[2]*PerpStrain[2]);     
      
       if(std::abs(temp)<0.0001){ // if maximal strain is normal to the cell plane storing strain direction instead as it should not be used
         PerpStrain[0]=eigenVectorStrain[0][Istrain];
         PerpStrain[1]=eigenVectorStrain[1][Istrain];
         PerpStrain[2]=eigenVectorStrain[2][Istrain];
       }
       
  

       // storing a measure for strain anisotropy in cell vector
       if (std::abs(maximalStrainValue)<0.0000001) cellData[cellIndex][strainAnIndex]=0;
       if (std::abs(maximalStrainValue)>= 0.0000001) cellData[cellIndex][strainAnIndex]=1-std::abs(maximalStrainValue2/maximalStrainValue);
       
       // storing   strain direction/value in cellData
       if (numVariableIndexLevel()==3 && (numVariableIndex(1)==1 || numVariableIndex(1)==2 || numVariableIndex(1)==3) ) {// storing maximal strain
         if (dimension==2)
           {
             cellData[cellIndex][variableIndex(1,0)]  =eigenVectorStrain[0][Istrain];
             cellData[cellIndex][variableIndex(1,0)+1]=eigenVectorStrain[1][Istrain];
             cellData[cellIndex][variableIndex(1,0)+3]=maximalStrainValue;  //maximal Strain Value is stored after its eigenvector
           }
         if (dimension==3)
           {
             cellData[cellIndex][variableIndex(1,0)]  =eigenVectorStrain[0][Istrain];
             cellData[cellIndex][variableIndex(1,0)+1]=eigenVectorStrain[1][Istrain];
             cellData[cellIndex][variableIndex(1,0)+2]=eigenVectorStrain[2][Istrain];
             cellData[cellIndex][variableIndex(1,0)+3]=maximalStrainValue;  //maximal Strain Value is stored after its eigenvector
           }
       }
       
       if (numVariableIndexLevel()==3 &&  numVariableIndex(1)==3  ) {//storing 2nd maximal strain
         if (dimension==2)
           {
             cellData[cellIndex][variableIndex(1,2)]  =eigenVectorStrain[0][Istrain2];
             cellData[cellIndex][variableIndex(1,2)+1]=eigenVectorStrain[1][Istrain2];
             cellData[cellIndex][variableIndex(1,2)+3]=maximalStrainValue2;  //2nd maximal Strain Value is stored after its eigenvector
           }
         if (dimension==3)
           {
             cellData[cellIndex][variableIndex(1,2)]  =eigenVectorStrain[0][Istrain2];
             cellData[cellIndex][variableIndex(1,2)+1]=eigenVectorStrain[1][Istrain2];
             cellData[cellIndex][variableIndex(1,2)+2]=eigenVectorStrain[2][Istrain2];
             cellData[cellIndex][variableIndex(1,2)+3]=maximalStrainValue2;  //2nd maximal Strain Value is stored after its eigenvector
           }
       }
       if (numVariableIndexLevel()==3 && ( numVariableIndex(1)==2 || numVariableIndex(1)==3 ) ) {//storing perpendicular to maximal strain
         if (dimension==2)
           {
             cellData[cellIndex][variableIndex(1,1)]  =PerpStrain[0];
             cellData[cellIndex][variableIndex(1,1)+1]=PerpStrain[1];
             cellData[cellIndex][variableIndex(1,1)+3]=maximalStrainValue;//maximal Strain Value is stored after its eigenvector
             //cellData[cellIndex][variableIndex(1,1)+3]=Area/restingArea;
           }
         if (dimension==3)
           {
             cellData[cellIndex][variableIndex(1,1)]  =PerpStrain[0];
             cellData[cellIndex][variableIndex(1,1)+1]=PerpStrain[1];
             cellData[cellIndex][variableIndex(1,1)+2]=PerpStrain[2];
             cellData[cellIndex][variableIndex(1,1)+3]=maximalStrainValue;//maximal Strain Value is stored after its eigenvector
             //cellData[cellIndex][variableIndex(1,1)+3]=Area/restingArea;
           }
       }


       

   
    
       
       //---- Anisotropic Correction Force-------------------------------
       double deltaF[3][3];
       
        
         for ( int i=0 ; i<3 ; ++i )  // from stress tensor(equipartitioning energy)
          for ( int j=0 ; j<3 ; ++j )
            deltaF[i][j]=(-deltaFTPK[i][j]);
      
        double  I1=trE;

        // energy
        EnergyIso =( (lambdaT/2)*I1*I1 + mioT*I2 );//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        EnergyAniso =( (deltaLam/2)*I4*I1 + deltaMio*I5 ); //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      
        //Forces of vertices   
        double Force[3][3];                                           
    
        Force[0][0]= (tensileStiffness[0]*Delta[0]+angularStiffness[1]*Delta[1]+angularStiffness[0]*Delta[2])*(position[1][0]-position[0][0])
          +(tensileStiffness[2]*Delta[2]+angularStiffness[2]*Delta[1]+angularStiffness[0]*Delta[0])*(position[2][0]-position[0][0])
          + deltaF[0][0]; 
        Force[0][1]= (tensileStiffness[0]*Delta[0]+angularStiffness[1]*Delta[1]+angularStiffness[0]*Delta[2])*(position[1][1]-position[0][1])
          +(tensileStiffness[2]*Delta[2]+angularStiffness[2]*Delta[1]+angularStiffness[0]*Delta[0])*(position[2][1]-position[0][1])
          + deltaF[0][1];  
        Force[0][2]= (tensileStiffness[0]*Delta[0]+angularStiffness[1]*Delta[1]+angularStiffness[0]*Delta[2])*(position[1][2]-position[0][2])
          +(tensileStiffness[2]*Delta[2]+angularStiffness[2]*Delta[1]+angularStiffness[0]*Delta[0])*(position[2][2]-position[0][2])
          + deltaF[0][2]; 
        Force[1][0]= (tensileStiffness[0]*Delta[0]+angularStiffness[0]*Delta[2]+angularStiffness[1]*Delta[1])*(position[0][0]-position[1][0])
          +(tensileStiffness[1]*Delta[1]+angularStiffness[2]*Delta[2]+angularStiffness[1]*Delta[0])*(position[2][0]-position[1][0])
          + deltaF[1][0];  
        Force[1][1]= (tensileStiffness[0]*Delta[0]+angularStiffness[0]*Delta[2]+angularStiffness[1]*Delta[1])*(position[0][1]-position[1][1])
          +(tensileStiffness[1]*Delta[1]+angularStiffness[2]*Delta[2]+angularStiffness[1]*Delta[0])*(position[2][1]-position[1][1])
          + deltaF[1][1];  
        Force[1][2]= (tensileStiffness[0]*Delta[0]+angularStiffness[0]*Delta[2]+angularStiffness[1]*Delta[1])*(position[0][2]-position[1][2])
          +(tensileStiffness[1]*Delta[1]+angularStiffness[2]*Delta[2]+angularStiffness[1]*Delta[0])*(position[2][2]-position[1][2])
          + deltaF[1][2];  
        
        Force[2][0]= (tensileStiffness[2]*Delta[2]+angularStiffness[0]*Delta[0]+angularStiffness[2]*Delta[1])*(position[0][0]-position[2][0])
          +(tensileStiffness[1]*Delta[1]+angularStiffness[1]*Delta[0]+angularStiffness[2]*Delta[2])*(position[1][0]-position[2][0])
          + deltaF[2][0];  
        Force[2][1]= (tensileStiffness[2]*Delta[2]+angularStiffness[0]*Delta[0]+angularStiffness[2]*Delta[1])*(position[0][1]-position[2][1])
          +(tensileStiffness[1]*Delta[1]+angularStiffness[1]*Delta[0]+angularStiffness[2]*Delta[2])*(position[1][1]-position[2][1])
          + deltaF[2][1];  
        Force[2][2]= (tensileStiffness[2]*Delta[2]+angularStiffness[0]*Delta[0]+angularStiffness[2]*Delta[1])*(position[0][2]-position[2][2])
          +(tensileStiffness[1]*Delta[1]+angularStiffness[1]*Delta[0]+angularStiffness[2]*Delta[2])*(position[1][2]-position[2][2])
          + deltaF[2][2];  
        
      
        
        // adding TRBSMT forces to the total vertexDerivs
        
        vertexDerivs[v1][0]+= Force[0][0];
        vertexDerivs[v1][1]+= Force[0][1];
        vertexDerivs[v1][2]+= Force[0][2];
        
        vertexDerivs[v2][0]+= Force[1][0];
        vertexDerivs[v2][1]+= Force[1][1];
        vertexDerivs[v2][2]+= Force[1][2];
        
        vertexDerivs[v3][0]+= Force[2][0];
        vertexDerivs[v3][1]+= Force[2][1];
        vertexDerivs[v3][2]+= Force[2][2];

 

        cellData[cellIndex][areaRatioIndex]  =Area/restingArea; //area ratio
        cellData[cellIndex][isoEnergyIndex]  =EnergyIso;        //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        cellData[cellIndex][anisoEnergyIndex]=EnergyAniso;  
        
  }
   

  double totalEnergyIso=0;
  double totalEnergyAniso=0;

  for( size_t cellIndex=0 ; cellIndex<numCells ; ++cellIndex ) {
   
    totalEnergyIso +=cellData[cellIndex][variableIndex(0,5)];
    totalEnergyAniso +=cellData[cellIndex][variableIndex(0,6)];

    Cell *  cell1=&(T.cell(cellIndex));
	
    std::vector<size_t>  neighbor(3);  
    neighbor[0]=(cell1->cellNeighbor(0))->index();
    neighbor[1]=(cell1->cellNeighbor(1))->index();
    neighbor[2]=(cell1->cellNeighbor(2))->index();

    double neighborweight=parameter(5);

    double StressTensor[3][3]={{0,0,0},{0,0,0},{0,0,0}};
    
    int counter=0;
    for (int nn=0 ; nn<3 ; nn++){
      if (neighbor[nn]<numCells && neighbor[nn]>-1){
	StressTensor[0][0]+=neighborweight*cellData[neighbor[nn]][stressTensorIndex];
	StressTensor[1][1]+=neighborweight*cellData[neighbor[nn]][stressTensorIndex+1];
	StressTensor[2][2]+=neighborweight*cellData[neighbor[nn]][stressTensorIndex+2];
	StressTensor[0][1]+=neighborweight*cellData[neighbor[nn]][stressTensorIndex+3];
	StressTensor[2][0]+=neighborweight*cellData[neighbor[nn]][stressTensorIndex+4];
	StressTensor[1][2]+=neighborweight*cellData[neighbor[nn]][stressTensorIndex+5];
	counter+=1;
      }
    }
    if(counter !=0){
      StressTensor[0][0]/=counter;
      StressTensor[1][1]/=counter;
      StressTensor[2][2]/=counter;
      StressTensor[0][1]/=counter;
      StressTensor[2][0]/=counter;
      StressTensor[1][2]/=counter;
    }

    StressTensor[0][0]+=(1-neighborweight)*cellData[cellIndex][stressTensorIndex];
    StressTensor[1][1]+=(1-neighborweight)*cellData[cellIndex][stressTensorIndex+1];
    StressTensor[2][2]+=(1-neighborweight)*cellData[cellIndex][stressTensorIndex+2];
    StressTensor[0][1]+=(1-neighborweight)*cellData[cellIndex][stressTensorIndex+3];
    StressTensor[2][0]+=(1-neighborweight)*cellData[cellIndex][stressTensorIndex+4];
    StressTensor[1][2]+=(1-neighborweight)*cellData[cellIndex][stressTensorIndex+5];
        

    //std:: cerr<<counter<<std::endl;

    StressTensor[0][2]=StressTensor[2][0];
    StressTensor[1][0]=StressTensor[0][1];
    StressTensor[2][1]=StressTensor[1][2];
       
    // eigenvalue/eigenvectors of  stress tensor in global coordinate system. (Jacobi method)
      
    double pi=3.14159265;
    double stressEpcilon=0.000001;    

    double RotAngle,Si,Co;
    double eigenVectorStress[3][3]={{1,0,0},{0,1,0},{0,0,1}};
    int I,J;    
    double pivot=1;
    
    while (pivot>stressEpcilon) {
      pivot=std::abs(StressTensor[1][0]);
      I=1;
      J=0;
      if (std::abs(StressTensor[2][0])>pivot) {
	pivot=std::abs(StressTensor[2][0]);
	I=2;
	J=0;
      }
      if (std::abs(StressTensor[2][1])>pivot) {
	pivot=std::abs(StressTensor[2][1]);
	I=2;
	J=1;
      }
      if (std::abs(StressTensor[I][I]-StressTensor[J][J])<stressEpcilon) {
	RotAngle=pi/4;
      }            
      else {
	RotAngle=0.5*std::atan((2*StressTensor[I][J])/(StressTensor[J][J]-StressTensor[I][I]));
      }
      Si=std::sin(RotAngle);
      Co=std::cos(RotAngle);
      double tempRot[3][3]={{1,0,0},{0,1,0},{0,0,1}};
      tempRot[I][I]=Co;
      tempRot[J][J]=Co;
      tempRot[I][J]=Si;
      tempRot[J][I]=-Si;
      double tempStress[3][3]={{0,0,0},{0,0,0},{0,0,0}};
      for (int r=0 ; r<3 ; r++) 
	for (int s=0 ; s<3 ; s++) 
	  for(int w=0 ; w<3 ; w++) 
	    tempStress[r][s]=tempStress[r][s]+StressTensor[r][w]*tempRot[w][s];
      
      for (int r=0 ; r<3 ; r++) 
	for (int s=0 ; s<3 ; s++) 
	  StressTensor[r][s]=0;
      
      for (int r=0 ; r<3 ; r++) 
	for (int s=0 ; s<3 ; s++) 
	  for(int w=0 ; w<3 ; w++) 
	    StressTensor[r][s]=StressTensor[r][s]+tempRot[w][r]*tempStress[w][s];
      
      for (int r=0 ; r<3 ; r++) 
	for (int s=0 ; s<3 ; s++) 
	  tempStress[r][s]=eigenVectorStress[r][s];
      
      for (size_t ii=0; ii<3; ++ii) 
	for (size_t jj=0; jj<3; ++jj) 
	  eigenVectorStress[ii][jj] = 0.0;
      
      for (int r=0 ; r<3 ; r++) 
	for (int s=0 ; s<3 ; s++) 
	  for(int w=0 ; w<3 ; w++) 
	    eigenVectorStress[r][s]=eigenVectorStress[r][s]+tempStress[r][w]*tempRot[w][s];
      
    }
    
    
    // maximal stress direction
    double maximalStressValue=StressTensor[0][0];
    int Istress=0;
    if (std::abs(StressTensor[1][1])>std::abs(maximalStressValue)) 
      {
	maximalStressValue=StressTensor[1][1];
	Istress=1;
      }
    if (std::abs(StressTensor[2][2])>std::abs(maximalStressValue)) 
      {
	maximalStressValue=StressTensor[2][2];
	Istress=2;
      }
    
       // 2nd maximal stress direction/value
      double maximalStressValue2;
      int Istress2,Istress3;
      if (Istress==0) {
        Istress2=1;
        Istress3=2;
      }
      if (Istress==1) {
        Istress2=0;
        Istress3=2;
      }
      if (Istress==2) {
        Istress2=0;
        Istress3=1;
      }
      if(std::abs(StressTensor[Istress3][Istress3])>std::abs(StressTensor[Istress2][Istress2])) {
        Istress2=Istress3;
      }
      maximalStressValue2=StressTensor[Istress2][Istress2];
      
      
    
    
    // storing a measure for stress anisotropy in cell vector
    if (std::abs(maximalStressValue)<0.000001) cellData[cellIndex][stressAnIndex ]=0;
    if (std::abs(maximalStressValue)>= 0.000001) cellData[cellIndex][stressAnIndex ]=1-std::abs(maximalStressValue2/maximalStressValue);
    
    
        
   //  // ---------------angles stress with circumferental direction- begin
   // double  temp=std::sqrt(eigenVectorStress[0][Istress]*eigenVectorStress[0][Istress]+                                // 14 --> Stress
   // 		   eigenVectorStress[1][Istress]*eigenVectorStress[1][Istress] );
   //  if(temp<0.000001){
   //    cellData[cellIndex][14]=pi/2;
   //  }
   //  else{
   //    cellData[cellIndex][14]=std::atan(eigenVectorStress[2][Istress]/temp);
   //  } 
   //  // ---------------angles stress with circumferental direction- end
    
    
    
    

    if (numVariableIndexLevel()==3 && (numVariableIndex(2)==1 || numVariableIndex(2)==2)) { // storing maximal stress
      if (dimension==2)
	{
	  cellData[cellIndex][variableIndex(2,0)]  =eigenVectorStress[0][Istress];
	  cellData[cellIndex][variableIndex(2,0)+1]=eigenVectorStress[1][Istress];
	  cellData[cellIndex][variableIndex(2,0)+3]=maximalStressValue;  //maximal Stress Value is stored after its eigenvector
	}
      if (dimension==3)
	{
	  cellData[cellIndex][variableIndex(2,0)]  =eigenVectorStress[0][Istress];
	  cellData[cellIndex][variableIndex(2,0)+1]=eigenVectorStress[1][Istress];
	  cellData[cellIndex][variableIndex(2,0)+2]=eigenVectorStress[2][Istress];
	  cellData[cellIndex][variableIndex(2,0)+3]=maximalStressValue;  //maximal Stress Value is stored after its eigenvector
	}
    }
    
    if (numVariableIndexLevel()==3 && numVariableIndex(2)==2 ) { // storing 2nd maximal stress
      if (dimension==2)
	{
	  cellData[cellIndex][variableIndex(2,1)]  =eigenVectorStress[0][Istress2];
	  cellData[cellIndex][variableIndex(2,1)+1]=eigenVectorStress[1][Istress2];
	  cellData[cellIndex][variableIndex(2,1)+3]=maximalStressValue2;  //2nd maximal Stress Value is stored after its eigenvector
	}
      if (dimension==3)
	{
	  cellData[cellIndex][variableIndex(2,1)]  =eigenVectorStress[0][Istress2];
	  cellData[cellIndex][variableIndex(2,1)+1]=eigenVectorStress[1][Istress2];
	  cellData[cellIndex][variableIndex(2,1)+2]=eigenVectorStress[2][Istress2];
	  cellData[cellIndex][variableIndex(2,1)+3]=maximalStressValue2;  //2nd maximal Stress Value is stored after its eigenvector
	}
    }
    
    
  }

    cellData[0][variableIndex(0,5)]=totalEnergyIso ;
    cellData[0][variableIndex(0,6)]=totalEnergyAniso ;

}
 
    

VertexFromTRBScenterTriangulationMT::
VertexFromTRBScenterTriangulationMT(std::vector<double> &paraValue, 
                                    std::vector< std::vector<size_t> > 
                                    &indValue ) 
{  
  // Do some checks on the parameters and variable indeces
  if( paraValue.size()!=11 ) { 
    std::cerr << "VertexFromTRBScenterTriangulationMT::"
              << "VertexFromTRBScenterTriangulationMT() "
              << "Uses 11 parameters: "
	      << "0,1: young modulus(matrix and fibre) " 
              << "2,3 : poisson ratio (longitudinal (MT) and transverse directions)"
	      << "4 : MF flag(0 constant material anisotropy ,1: material anisotropy via FiberModel " 
              << "5 : neighbor weight for correlation " 
	      << "6 : max stress/strain  " 
              << "7 : 2nd flag(0: plane strain, 1: plane stress) " 
              << "8 : MT direction angle"
              << "9 (MT update flag): 0:for no feedback or direct feedback by indices,"
              << "                    1:for MT direction from 7th parameter TETA, 2:force to Stress,  "
              << "                    3: force to Strain ,4:force to perp-strain "
	      << "10 : unused parameter for temporary use "<< std::endl;
    
    exit(0);
  }
    
  if( (indValue.size()!=2 && indValue.size()!=4) || 
      indValue[0].size()!=11 || indValue[1].size()!=1 ||
      (indValue.size()==4 && (indValue[2].size()!=0 && indValue[2].size()!=1 && indValue[2].size()!=2 && indValue[2].size()!=3)) ||
      (indValue.size()==4 && (indValue[3].size()!=0 && indValue[3].size()!=1 && indValue[3].size()!=2 ))
      ) { 
    std::cerr << "VertexFromTRBScenterTriangulationMT::"
	      << "VertexFromTRBScenterTriangulationMT() "
	      << "11 indices including (0)Wall length index and (1)MT direction initial index and (2)strain and (3)stress  "
	      << "anisotropy indices and indices for storing (4)area ratio, (5)isotropic energy and (6)anisotropic energy  " 
	      << "and (7)young_fiber and (8)MTstress  and strart indices for storing (9)stress tensor(6 elements) and      "
	      << "(10)normal vector to cell plane (3 elements) given in first level." 
	      << "Start of additional Cell variable indices (center(x,y,z) "
	      << "L_1,...,L_n, n=num vertex) is given in second level (typically at the end)." 
              << "Optionally two additional levels can be given where the strain, perpendicular to strain and  2nd strain "
	      << "can be stored in 3rd, stress and 2nd stress can be stored in 4th level "
	      << "directions/values(dx dy dz value) can be stored at given indices."
              << "If no index given at 3rd level, strain will not be stored, if one index given strain will be stored and  " 
	      << "if two indices are given maximal and perpendicular strain will be stored and if 3 indices are given 2nd  "
	      << "strain direction values will be stored at 3rd index"
              << "If no index given at 4th level, stress will not be stored, if one index given stress will be stored  "
              << "and if two indices are given maximal and 2nd stress will be stored at 1st and 2nd index respectively "
              << std::endl;
    exit(0);
  }
  


  // Set the variable values
  setId("VertexFromTRBScenterTriangulationMT");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  // Set the parameter identities
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "Y_mod_M";   // Matrix Young modulus
  tmp[1] = "Y_mod_F";   // Fiber Young modulus
  tmp[2] = "P_ratio_L"; // Longitudinal Poisson ratio
  tmp[3] = "P_ratio_T"; // Transverse Poisson ratio
  tmp[4] = "MF flag";
  tmp[5] = "neigborweight";
  tmp[6] = "stressmax";
  tmp[7] = "Strain-Stress flag";
  tmp[8] = "TETA anisotropy";
  tmp[9] = "MT update flag";
  tmp[10] = "unused";
 
  setParameterId( tmp );
  
  if( parameter(2)<0 || parameter(2)>=0.5 || parameter(3)<0 || parameter(3)>=0.5 ) {
    std::cerr << " VertexFromTRBScenterTriangulationMT::"
              << " VertexFromTRBScenterTriangulationMT() "
              << " poisson ratio must be 0 <= p < 0.5 " << std::endl;
    exit(0);
  }
  
  // if( parameter(4)!=0 && parameter(4)!=1 && parameter(4)!=2 ) {
  //   std::cerr << " VertexFromTRBScenterTriangulationMT::"
  // 	      << " VertexFromTRBScenterTriangulationMT() "
  // 	      << " 5th parameter must be 0 or 1  "
  //             << " 0: constant material anisotropy, 1: material anisotropy via FiberModel  2: for energy landscape(see the code, do not use until you are sure)" << std::endl;
  //   exit(0);
  // }

  
  if( parameter(7)!=0 && parameter(7)!=1 ) {
    std::cerr << " VertexFromTRBScenterTriangulationMT::"
	      << " VertexFromTRBScenterTriangulationMT() "
	      << " 6th parameter must be 0 or 1(0:plane strain, 1:plane stress) " << std::endl;
    exit(0);
  }

  if( parameter(9)!=0 && parameter(9)!=1 && parameter(9)!=2 && parameter(9)!=3 && parameter(9)!=4) {
    std::cerr << " VertexFromTRBScenterTriangulationMT::"
              << " VertexFromTRBScenterTriangulationMT() "
              << " 8th parameter must be 0/1/2/3/4"
              << " 0: for no feedback or direct feedback by indices "
              << " 1: for MT direction from 7th parameter TETA "
              << " 2: force to Stress "
              << " 3: force to Strain "
              << " 4: force to perp-strain " << std::endl;
    exit(0);
  }
}


void VertexFromTRBScenterTriangulationMT::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) {
 

  
  //Do the update for each cell
  size_t dimension = 3;
  assert (dimension==vertexData[0].size());
  size_t numCells = T.numCell();
  size_t numVertices = T.numVertex();
  size_t wallLengthIndex = variableIndex(0,0);
  size_t comIndex = variableIndex(1,0);
  size_t lengthInternalIndex = comIndex+dimension;
  
  
  double TotalVolume=0;
  double deltaVolume=0;
  for(size_t vertexIndex=0; vertexIndex<numVertices; ++vertexIndex){ // stimating volume for equilibrium 
    TotalVolume +=std::sqrt(vertexData[vertexIndex][0]*vertexData[vertexIndex][0] +
                            vertexData[vertexIndex][1]*vertexData[vertexIndex][1] +
                            vertexData[vertexIndex][2]*vertexData[vertexIndex][2] );
  }
  deltaVolume=TotalVolume-cellData[0][25];
  // cellData[0][25]=TotalVolume;
  // cellData[0][24]=deltaVolume;
 
  size_t MTindex           =variableIndex(0,1);	 
  size_t strainAnIndex     =variableIndex(0,2);	
  size_t stressAnIndex     =variableIndex(0,3);	
  size_t areaRatioIndex    =variableIndex(0,4);	
  size_t isoEnergyIndex    =variableIndex(0,5);	
  size_t anisoEnergyIndex  =variableIndex(0,6);	
  size_t youngLIndex       =variableIndex(0,7);	
  size_t MTstressIndex     =variableIndex(0,8);	
  size_t stressTensorIndex =variableIndex(0,9);	
  size_t normalVectorIndex =variableIndex(0,10);



  //std::cerr<<"from here"<<std::endl;
  
  for (size_t cellIndex=0 ; cellIndex<numCells ; ++cellIndex) {
    size_t numWalls = T.cell(cellIndex).numWall();
    

    if(  T.cell(cellIndex).numVertex()!= numWalls ) {
     
      std::cerr << "VertexFromTRBScenterTriangulationMT::derivs() same number of vertices and walls."
		<< " Not for cells with " << T.cell(cellIndex).numWall() << " walls and "
		<< T.cell(cellIndex).numVertex() << " vertices!"	
		<< std::endl;
      exit(-1);
    }
    double youngMatrix= parameter(0);    
    double youngFiber = parameter(1); 
    double poissonL   = parameter(2);    
    double poissonT   = parameter(3);
    double TETA       = parameter(8);  

  
 

    // if ( std::sqrt(cellData[cellIndex][comIndex  ]*cellData[cellIndex][comIndex  ]
    //                +cellData[cellIndex][comIndex+1]*cellData[cellIndex][comIndex+1])<30)
    //   cellData[cellIndex][youngLIndex]=20;
    // else
    //   cellData[cellIndex][youngLIndex]=200;

 
    double youngL=1;
    double youngT=1;
    if( parameter(4)==0 ){ // constant anisotropic material
      youngL = youngMatrix+youngFiber;
      youngT = youngMatrix; 
    }
    if( parameter(4)==1){  // material anisotropy via FiberModel
      youngL = cellData[cellIndex][youngLIndex]; 
      youngT = 2*youngMatrix+youngFiber-youngL; 
    }
    if( parameter(4)==2){  // for varrying material anisotropy with constant overall stiffness for energy landscape
      youngL =youngFiber; 
      youngT =youngMatrix-youngL;  // here youngMatrix is total stiffness
    }
    if( parameter(4)==3){  // for varrying material anisotropy with constant overall stiffness for energy landscape
      double totalElast=140;
      double slope=1;
      //youngL =youngMatrix+youngFiber; 
      // youngT =youngMatrix+std::sqrt(
      //                               ((youngFiber-totalElast)*4*totalElast)
      //                               /
      //                               ((2-3.1415)*3.1415)
      //                               );  // here youngMatrix is total stiffness
      //youngT =youngMatrix+totalElast*std::sqrt(
      //                                         -(std::log((youngFiber/totalElast)-0.55))
      //                                         /6
      //                                         );
      //youngT =youngMatrix+2*(totalElast-youngFiber)/(3.1415-2);

      youngL=youngMatrix+youngFiber; 
      youngT =youngMatrix+(youngFiber-slope*totalElast)/(1-2*slope); 
    }

    // if( parameter(4)<0){  // for heterogeneous stiffness(adhoc)
    //   double hFactor=0;
    //   double Hthreshold=0.01;
    //   if (std::abs(cellData[cellIndex][youngLIndex])>Hthreshold){
    //     hFactor=Hthreshold;
    //   }
    //   else{
    //     hFactor=std::abs(cellData[cellIndex][youngLIndex]); // take the heterogeneous info from this index
    //   }
    //   hFactor *=std::abs(parameter(4)); // factor for heterogeneity      
      
    //   youngL = youngMatrix+youngFiber-hFactor*(youngMatrix+youngFiber); 
    //   youngT = youngMatrix-hFactor*(youngMatrix); 
    // }


    if( parameter(4)<0){  // for spatial elasticity (very adhoc)
      double x0=90;
      double y0=120;

      double xx=cellData[cellIndex][comIndex]-x0;
      double yy=cellData[cellIndex][comIndex+1]-y0;

      youngL = youngMatrix+youngFiber*(1-std::exp(parameter(4)*(xx*xx+yy*yy)));
      youngT = youngMatrix; 
      cellData[cellIndex][13]=youngL;
    }
   


    double lambdaL, mioL, lambdaT, mioT; // ,lambdaTmatrix, mioTmatrix;
    
    if (parameter(7)==0){      
      // Lame coefficients based on plane strain (for 3D 0<poisson<0.5)
      lambdaL=youngL*poissonL/((1+poissonL)*(1-2*poissonL));
      mioL=youngL/(2*(1+poissonL));
      lambdaT=youngT*poissonT/((1+poissonT)*(1-2*poissonT));
      mioT=youngT/(2*(1+poissonT));

   
    } 
    else{      
      // Lame coefficients based on plane stress (for 3D 0<poisson<0.5)
      lambdaL=youngL*poissonL/(1-poissonL*poissonL);
      mioL=youngL/(2*(1+poissonL));
      lambdaT=youngT*poissonT/(1-poissonT*poissonT);
      mioT=youngT/(2*(1+poissonT));
     
      // lambdaTmatrix=youngMatrixA*poissonT/(1-poissonT*poissonT);
      // mioTmatrix=youngMatrixA/(2*(1+poissonT));
    }
    
    // Lame coefficients based on delin. paper (for 2D 0<poisson<1)
    // double lambdaL=youngL*poissonL/(1-poissonL*poissonL);
    // double mioL=youngL/(1+poissonL);
    // double lambdaT=youngT*poissonT/(1-poissonT*poissonT);
    // double mioT=youngT/(1+poissonT);
    

    double StrainCellGlobal[3][3]={{0,0,0},{0,0,0},{0,0,0}};
    double StressCellGlobal[3][3]={{0,0,0},{0,0,0},{0,0,0}};
    double normalGlob[3]={0,0,0};
    double TotalCellRestingArea=0;
    double TotalCellArea=0;
    double EnergyIso=0;     
    //double EnergyIsoFiber=0;            
    double EnergyAniso=0;
    double strainZ=0;

    if ( parameter(9)==1 ) {  // aniso direction from TETA 
      
      cellData[cellIndex][variableIndex(0,1)]=std::cos(TETA);    
      cellData[cellIndex][variableIndex(0,1)+1]=std::sin(TETA);
      cellData[cellIndex][variableIndex(0,1)+2]=0;
      if(parameter(10)!=100){
	if(cellIndex==0 || cellIndex==2 ){
	  cellData[cellIndex][variableIndex(0,1)]=std::cos(TETA);    
	  cellData[cellIndex][variableIndex(0,1)+1]=std::sin(TETA);
	  cellData[cellIndex][variableIndex(0,1)+2]=0;
	}
	if(cellIndex==1 || cellIndex==3 ){
	  cellData[cellIndex][variableIndex(0,1)]=std::cos(parameter(10));    
	  cellData[cellIndex][variableIndex(0,1)+1]=std::sin(parameter(10));
	  cellData[cellIndex][variableIndex(0,1)+2]=0;
	}
      }
    }
    
   

    // Aniso vector in current shape in global coordinate system
      double AnisoCurrGlob[3];
      AnisoCurrGlob[0] = cellData[cellIndex][variableIndex(0,1)];  
      AnisoCurrGlob[1] = cellData[cellIndex][variableIndex(0,1)+1];
      AnisoCurrGlob[2] = cellData[cellIndex][variableIndex(0,1)+2];

      
     
    // One triangle per 'vertex' in cyclic order
    for (size_t wallindex=0; wallindex<numWalls; ++wallindex) { 
      size_t kPlusOneMod = (wallindex+1)%numWalls;
      //size_t v1 = com;
      size_t v2 = T.cell(cellIndex).vertex(wallindex)->index();
      size_t v3 = T.cell(cellIndex).vertex(kPlusOneMod)->index();
      //size_t w1 = internal wallindex
      size_t w2 = T.cell(cellIndex).wall(wallindex)->index();
      //size_t w3 = internal wallindex+1

      // Position matrix holds in rows positions for com, vertex(wallindex), vertex(wallindex+1)
      DataMatrix position(3,vertexData[v2]);
      for (size_t d=0; d<dimension; ++d)
	position[0][d] = cellData[cellIndex][comIndex+d]; // com position
      //position[1] = vertexData[v2]; // given by initiation
      position[2] = vertexData[v3];
      //position[0][2] z for vertex 1 of the current element
      
      

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////


      // Resting lengths are from com-vertex(wallindex), vertex(wallindex)-vertex(wallindex+1) (wall(wallindex)), com-vertex(wallindex+1)
      std::vector<double> restingLength(3);
      restingLength[0] = cellData[cellIndex][lengthInternalIndex + wallindex];
      restingLength[1] = wallData[w2][wallLengthIndex];
      restingLength[2] = cellData[cellIndex][lengthInternalIndex + kPlusOneMod];
      
      // Lengths are from com-vertex(wallindex), vertex(wallindex)-vertex(wallindex+1) (wall(wallindex)), com-vertex(wallindex+1)
      std::vector<double> length(3);
      length[0] = std::sqrt( (position[0][0]-position[1][0])*(position[0][0]-position[1][0]) +
			     (position[0][1]-position[1][1])*(position[0][1]-position[1][1]) +
			     (position[0][2]-position[1][2])*(position[0][2]-position[1][2]) );
      
      length[1] = T.wall(w2).lengthFromVertexPosition(vertexData);
        
      length[2] = std::sqrt( (position[0][0]-position[2][0])*(position[0][0]-position[2][0]) +
			     (position[0][1]-position[2][1])*(position[0][1]-position[2][1]) +
			     (position[0][2]-position[2][2])*(position[0][2]-position[2][2]) );

      
      // Area of the element (using Heron's formula)                                      
      double restingArea=std::sqrt( ( restingLength[0]+restingLength[1]+restingLength[2])*
                                    (-restingLength[0]+restingLength[1]+restingLength[2])*
                                    ( restingLength[0]-restingLength[1]+restingLength[2])*
                                    ( restingLength[0]+restingLength[1]-restingLength[2])  )*0.25;
      
      //double currentArea=std::sqrt( ( length[0]+length[1]+length[2])*
      //                              (-length[0]+length[1]+length[2])*
      //                              ( length[0]-length[1]+length[2])*
      //                              ( length[0]+length[1]-length[2])  )*0.25;
      
      
      //Angles of the element ( assuming the order: 0,L0,1,L1,2,L2 )
      std::vector<double> Angle(3);
      // can be changed by cotan(A)=.25*sqrt(4*b*b*c*c/K-1)
      Angle[0]=std::acos(  (restingLength[0]*restingLength[0]+restingLength[2]*restingLength[2]-restingLength[1]*restingLength[1])/
                           (restingLength[0]*restingLength[2]*2)    );
      Angle[1]=std::acos(  (restingLength[0]*restingLength[0]+restingLength[1]*restingLength[1]-restingLength[2]*restingLength[2])/
                           (restingLength[0]*restingLength[1]*2)    );
      Angle[2]=std::acos(  (restingLength[1]*restingLength[1]+restingLength[2]*restingLength[2]-restingLength[0]*restingLength[0])/
                           (restingLength[1]*restingLength[2]*2)    );
      
      //Tensile Stiffness
      std::vector<double>  tensileStiffness(3);
      //double tensileStiffness[3];
      double temp = 1.0/(restingArea*16);                                      
      std::vector<double> cotan(3);
      cotan[0] = 1.0/std::tan(Angle[0]);
    cotan[1] = 1.0/std::tan(Angle[1]);
    cotan[2] = 1.0/std::tan(Angle[2]);

    
    // if(Angle[0]<0.024 ||Angle[1]<0.024 ||Angle[2]<0.024 ){
    //   std::cerr<<"in mechanicalTRBS cell "<<cellIndex<<" wall  "<<wallindex<<" angles 0,1,2  "<<Angle[0]<<"  "<<Angle[1]<<"  "<<Angle[2]<<std::endl;
    //   std::cerr<<"in mechanicalTRBS cotan  0,1,2  "<<cotan[0]<<"  "<<cotan[1]<<"  "<<cotan[2]<<std::endl;
    // } 
    //the force is calculated based on Transverse coefficients
    //Longitudinal coefficients are considered in deltaF
    tensileStiffness[0]=(2*cotan[2]*cotan[2]*(lambdaT+2*mioT)+2*mioT)*temp;
    tensileStiffness[1]=(2*cotan[0]*cotan[0]*(lambdaT+2*mioT)+2*mioT)*temp;
    tensileStiffness[2]=(2*cotan[1]*cotan[1]*(lambdaT+2*mioT)+2*mioT)*temp;
    
    //Angular Stiffness
    std::vector<double> angularStiffness(3);
    angularStiffness[0]=(2*cotan[1]*cotan[2]*(lambdaT+2*mioT)-2*mioT)*temp;                          
    angularStiffness[1]=(2*cotan[0]*cotan[2]*(lambdaT+2*mioT)-2*mioT)*temp;
    angularStiffness[2]=(2*cotan[0]*cotan[1]*(lambdaT+2*mioT)-2*mioT)*temp;
    
    //Calculate biquadratic strains  
    std::vector<double> Delta(3);
    Delta[0]=(length[0])*(length[0])-(restingLength[0])*(restingLength[0]);
    Delta[1]=(length[1])*(length[1])-(restingLength[1])*(restingLength[1]);
    Delta[2]=(length[2])*(length[2])-(restingLength[2])*(restingLength[2]);

    //Area of the element (using Heron's formula)                                      
    double Area=std::sqrt( ( length[0]+length[1]+length[2])*
                           (-length[0]+length[1]+length[2])*
                           ( length[0]-length[1]+length[2])*
                           ( length[0]+length[1]-length[2])  )*0.25;
    
    // calculating the angles between shape vectors and anisotropy direction in resting shape when anisotropy vector is provided in current shape
    
    //Current shape local coordinate of the element  (counterclockwise ordering of nodes/edges)
      double CurrentAngle1=std::acos(  (length[0]*length[0]+length[1]*length[1]-length[2]*length[2])/
                                       (length[0]*length[1]*2)    );

      double Qa=std::cos(CurrentAngle1)*length[0];
      double Qc=std::sin(CurrentAngle1)*length[0];
      double Qb=length[1];
      // shape vector matrix = inverse of coordinate matrix ( only first two elements i.e. ShapeVector[3][2] )      
      double ShapeVectorCurrent[3][3]={ {  0   ,       1/Qc      , 0 }, 
                                        {-1/Qb , (Qa-Qb)/(Qb*Qc) , 1 },       
                                        { 1/Qb ,     -Qa/(Qb*Qc) , 0 }  };
            
     
      
      // Local coordinates of the resting shape ( counterclockwise )
      double RestingAngle1=std::acos(  (restingLength[0]*restingLength[0]+restingLength[1]*restingLength[1]-restingLength[2]*restingLength[2])/
                                       (restingLength[0]*restingLength[1]*2)    );

      double Pa=std::cos(RestingAngle1)*restingLength[0];
      double Pc=std::sin(RestingAngle1)*restingLength[0];
      double Pb=restingLength[1];

      // shape vector matrix in resting shape in local coordinate system  = inverse of coordinate matrix ( only first two elements i.e. ShapeVectorResting[3][2] )      
      double ShapeVectorResting[3][3]={ {  0   ,       1/Pc      , 0 }, 
                                        {-1/Pb , (Pa-Pb)/(Pb*Pc) , 1 },       
                                        { 1/Pb ,     -Pa/(Pb*Pc) , 0 }  };
    
      // //Strain tensor  (clockwise ordering of nodes/edges)
      // double CurrentAngle2=std::acos(  (length[1]*length[1]+length[2]*length[2]-length[0]*length[0])/
      //                                  (length[1]*length[2]*2)    );

      // double Qa=std::cos(CurrentAngle2)*length[2];
      // double Qb=length[1];
      // double Qc=std::sin(CurrentAngle2)*length[2];
      
      // double ShapeVectorCurrent[3][2]={ {  0   ,       1/Qc      }, 
      //                                   { 1/Qb ,     -Qa/(Qb*Qc) },       
      //                                   {-1/Qb , (Qa-Qb)/(Qb*Qc) }  };

      

      
     

      // Rotation Matrix for changing coordinate systems for both Local to Global( Strain Tensor) and Global to Local( Aniso Vector in the current shape)
      double rotation[3][3];  

      double tempA=std::sqrt((position[2][0]-position[1][0])*(position[2][0]-position[1][0])+
                             (position[2][1]-position[1][1])*(position[2][1]-position[1][1])+
                             (position[2][2]-position[1][2])*(position[2][2]-position[1][2])  );

      double tempB=std::sqrt((position[0][0]-position[1][0])*(position[0][0]-position[1][0])+
                             (position[0][1]-position[1][1])*(position[0][1]-position[1][1])+
                             (position[0][2]-position[1][2])*(position[0][2]-position[1][2])  );


      double Xcurrent[3];      
      Xcurrent[0]= (position[2][0]-position[1][0])/tempA;
      Xcurrent[1]= (position[2][1]-position[1][1])/tempA;
      Xcurrent[2]= (position[2][2]-position[1][2])/tempA;
      
      double Bcurrent[3];      
      Bcurrent[0]= (position[0][0]-position[1][0])/tempB;
      Bcurrent[1]= (position[0][1]-position[1][1])/tempB;
      Bcurrent[2]= (position[0][2]-position[1][2])/tempB;

      double Zcurrent[3];      
      Zcurrent[0]= Xcurrent[1]*Bcurrent[2]-Xcurrent[2]*Bcurrent[1];
      Zcurrent[1]= Xcurrent[2]*Bcurrent[0]-Xcurrent[0]*Bcurrent[2];
      Zcurrent[2]= Xcurrent[0]*Bcurrent[1]-Xcurrent[1]*Bcurrent[0];
      
      tempA=std:: sqrt(Zcurrent[0]*Zcurrent[0]+Zcurrent[1]*Zcurrent[1]+Zcurrent[2]*Zcurrent[2]);
      Zcurrent[0]=Zcurrent[0]/tempA;
      Zcurrent[1]=Zcurrent[1]/tempA;
      Zcurrent[2]=Zcurrent[2]/tempA;

      double Ycurrent[3];      
      Ycurrent[0]= Zcurrent[1]*Xcurrent[2]-Zcurrent[2]*Xcurrent[1];
      Ycurrent[1]= Zcurrent[2]*Xcurrent[0]-Zcurrent[0]*Xcurrent[2];
      Ycurrent[2]= Zcurrent[0]*Xcurrent[1]-Zcurrent[1]*Xcurrent[0];


      rotation[0][0]=Xcurrent[0];
      rotation[1][0]=Xcurrent[1];
      rotation[2][0]=Xcurrent[2];

      rotation[0][1]=Ycurrent[0];
      rotation[1][1]=Ycurrent[1];
      rotation[2][1]=Ycurrent[2];

      rotation[0][2]=Zcurrent[0];
      rotation[1][2]=Zcurrent[1];
      rotation[2][2]=Zcurrent[2];  

      // AnisoCurrGlobPlus[2]     AnisoCurrGlobMinus[3];    
       
      // rotating the anisotropy vector from global coordinate system to the local one in the current shape
      double AnisoCurrLocal[3];
      AnisoCurrLocal[0]=rotation[0][0]*AnisoCurrGlob[0]+rotation[1][0]*AnisoCurrGlob[1]+rotation[2][0]*AnisoCurrGlob[2];
      AnisoCurrLocal[1]=rotation[0][1]*AnisoCurrGlob[0]+rotation[1][1]*AnisoCurrGlob[1]+rotation[2][1]*AnisoCurrGlob[2];
      AnisoCurrLocal[2]=rotation[0][2]*AnisoCurrGlob[0]+rotation[1][2]*AnisoCurrGlob[1]+rotation[2][2]*AnisoCurrGlob[2];
      
     
     
      //std::cerr<< "cell "<< cellIndex<< " anisoVector current local "<<AnisoCurrLocal[0]<<"  "<<AnisoCurrLocal[1]<<"  "<<AnisoCurrLocal[2]<< std::endl;
      // Center of Mass for current shape in local coordinate
      double CMCurrentLocal[2]={(Qa+Qb)/3, Qc/3};
      
      // Tip of the ansotropy vector drown from the center of mass in the current shape
      double  ACurrLocal[2] = {CMCurrentLocal[0]+AnisoCurrLocal[0],CMCurrentLocal[1]+AnisoCurrLocal[1]};
     
    
      // Baricentric Coordinates of tip of anisotropy vector in the current shape wihich is equivalent to the baricentric coordinate of the corresponding point in the resting shape
      double ABari[3];
      ABari[0]=ShapeVectorCurrent[0][0]*ACurrLocal[0]+ShapeVectorCurrent[0][1]*ACurrLocal[1]+ShapeVectorCurrent[0][2];
      ABari[1]=ShapeVectorCurrent[1][0]*ACurrLocal[0]+ShapeVectorCurrent[1][1]*ACurrLocal[1]+ShapeVectorCurrent[1][2];
      ABari[2]=ShapeVectorCurrent[2][0]*ACurrLocal[0]+ShapeVectorCurrent[2][1]*ACurrLocal[1]+ShapeVectorCurrent[2][2];
      //std::cerr<< "cell "<< cellIndex<< " baricentric coor of tip of anisoVector in local current "<<ABari[0]<<"  "<<ABari[1]<<"  "<<ABari[2]<< std::endl;
      
      
      
      // Local coordinates of tip of AnisoVector in the resting shape from multyplying ABari by position matrix in the local resting shape coordinates           
      double ARestLocal[2];
      ARestLocal[0]=Pa*ABari[0]+Pb*ABari[2];
      ARestLocal[1]=Pc*ABari[0];
      //std::cerr<< "cell "<< cellIndex<< " tip of anisoVector in local rest "<<ARestLocal[0]<<"  "<<ARestLocal[1]<< std::endl;     
      
     
      // Center of Mass for resting shape in local coordinate
      double CMRestingLocal[2]={(Pa+Pb)/3, Pc/3};
            
      // Aniso Vector in the resting shape in local coordinate system
      double AnisoRestLocal[2]={ARestLocal[0]-CMRestingLocal[0],ARestLocal[1]-CMRestingLocal[1]};
      //std::cerr<< "cell "<< cellIndex<< " anisoVector Rest "<<AnisoRestLocal[0]<<"  "<<AnisoRestLocal[1]<<"  "<<AnisoRestLocal[2]<< std::endl;    
      
      
      // Anisotropy measure or magnitude of the projection of anisotropy vector on the cell plane in the resting shape 
      //this measure is considered as a factor controling the anisotropy properties of the cell
      double AnisoMeasure=std::sqrt(AnisoRestLocal[0]*AnisoRestLocal[0]+AnisoRestLocal[1]*AnisoRestLocal[1]);
      // std::cerr<< "cell "<< cellIndex<<" AnisoMeasure "<< AnisoMeasure<<std::endl;
      
      
        
      // choosing a random normalized dirrection for anisotropy if AnisoVector is close to perpendicular to the cell plane
      if ( AnisoMeasure<0.0001) {
        double randomAngle=((rand()%360)*2*3.14159265)/360;
        
        AnisoRestLocal[0]=std::cos(randomAngle);
        AnisoRestLocal[1]=std::sin(randomAngle);
      }  
      else {// normalizing the anisoVector if it is not random
        AnisoRestLocal[0]=AnisoRestLocal[0]/AnisoMeasure;
        AnisoRestLocal[1]=AnisoRestLocal[1]/AnisoMeasure;        
      }
      

      //Anisotropic Correction is based on difference between Lam Coefficients of Longitudinal and Transverse dirrections:
      // double deltaLam=AnisoMeasure*(lambdaL-lambdaT);
      // double deltaMio=AnisoMeasure*(mioL-mioT);
      double deltaLam=lambdaL-lambdaT;
      double deltaMio=mioL-mioT;  

      // double deltaLamIsoFiber=lambdaT-lambdaTmatrix;
      // double deltaMioIsoFiber=mioT-mioTmatrix;  
 

      //Angles between anisotropy vector and shape vectors for calculating the terms like a.Di , 
      //teta(k) = acos((dot(Anisorest,Dk))/(norm(Anisorest)*norm(Dk))),
      std::vector<double> teta(3);
      teta[0] = std::acos(  (ShapeVectorResting[0][0]*AnisoRestLocal[0]
                             +ShapeVectorResting[0][1]*AnisoRestLocal[1])/
                            std::sqrt(ShapeVectorResting[0][0]*ShapeVectorResting[0][0]
                                      +ShapeVectorResting[0][1]*ShapeVectorResting[0][1]+0.00000001) );
      
      teta[1] = std::acos(  (ShapeVectorResting[1][0]*AnisoRestLocal[0]
                             +ShapeVectorResting[1][1]*AnisoRestLocal[1])/
                            std::sqrt(ShapeVectorResting[1][0]*ShapeVectorResting[1][0]
                                      +ShapeVectorResting[1][1]*ShapeVectorResting[1][1]+0.00000001) );
      
      teta[2] = std::acos(  (ShapeVectorResting[2][0]*AnisoRestLocal[0]
                             +ShapeVectorResting[2][1]*AnisoRestLocal[1])/
                            std::sqrt(ShapeVectorResting[2][0]*ShapeVectorResting[2][0]
                                      +ShapeVectorResting[2][1]*ShapeVectorResting[2][1]+0.00000001) );

     
    
       // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> STRAIN and STRESS TENSOR (BEGIN) <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      
      // deformation gradiant tensor F =Sigma i=1,2,3 Qi x Di
      // strain tensor in resting shape E=0.5(FtF-I)
      // trE
      // B=FFt
      // axa (direct product of aniso vector in resting shape)
      // atEa
      // E(axa) and (axa)E
      double trE=( Delta[1]*cotan[0]+ Delta[2]*cotan[1]+Delta[0]*cotan[2])/(4*restingArea);
      
     

      double directAniso[2][2]={{AnisoRestLocal[0]*AnisoRestLocal[0],AnisoRestLocal[0]*AnisoRestLocal[1]},
                                {AnisoRestLocal[1]*AnisoRestLocal[0],AnisoRestLocal[1]*AnisoRestLocal[1]}};
      
      double positionLocal[3][2]={ {Qa , Qc}, 
                                   {0  , 0 },  
                                   {Qb , 0 }  };
      
      double DeformGrad[2][2]={{0,0},{0,0}}; // F=Sigma i Qi x Di   <-------------------------------------------------------------------
      for ( int i=0 ; i<3 ; ++i ) {
        DeformGrad[0][0]=DeformGrad[0][0]+positionLocal[i][0]*ShapeVectorResting[i][0];
        DeformGrad[1][0]=DeformGrad[1][0]+positionLocal[i][1]*ShapeVectorResting[i][0];
        DeformGrad[0][1]=DeformGrad[0][1]+positionLocal[i][0]*ShapeVectorResting[i][1];
        DeformGrad[1][1]=DeformGrad[1][1]+positionLocal[i][1]*ShapeVectorResting[i][1];
      } 

      double LeftCauchy[2][2]; // B=FFt
      LeftCauchy[0][0]=DeformGrad[0][0]*DeformGrad[0][0]+DeformGrad[0][1]*DeformGrad[0][1];
      LeftCauchy[1][0]=DeformGrad[1][0]*DeformGrad[0][0]+DeformGrad[1][1]*DeformGrad[0][1];
      LeftCauchy[0][1]=DeformGrad[0][0]*DeformGrad[1][0]+DeformGrad[0][1]*DeformGrad[1][1];
      LeftCauchy[1][1]=DeformGrad[1][0]*DeformGrad[1][0]+DeformGrad[1][1]*DeformGrad[1][1];


      double Egreen[2][2];//E=0.5(C-I)
      Egreen[0][0]=0.5*(DeformGrad[0][0]*DeformGrad[0][0]+DeformGrad[1][0]*DeformGrad[1][0]-1);
      Egreen[1][0]=0.5*(DeformGrad[0][1]*DeformGrad[0][0]+DeformGrad[1][1]*DeformGrad[1][0]);
      Egreen[0][1]=0.5*(DeformGrad[0][0]*DeformGrad[0][1]+DeformGrad[1][0]*DeformGrad[1][1]);
      Egreen[1][1]=0.5*(DeformGrad[0][1]*DeformGrad[0][1]+DeformGrad[1][1]*DeformGrad[1][1]-1);
      
      double E2[2][2]; // used for energy calculation only
      E2[0][0]=Egreen[0][0]*Egreen[0][0]+Egreen[0][1]*Egreen[1][0];
      E2[1][0]=Egreen[1][0]*Egreen[0][0]+Egreen[1][1]*Egreen[1][0];
      E2[0][1]=Egreen[0][0]*Egreen[0][1]+Egreen[0][1]*Egreen[1][1];
      E2[1][1]=Egreen[1][0]*Egreen[0][1]+Egreen[1][1]*Egreen[1][1];

      double I2=E2[0][0]+E2[1][1]; //trE2 used for energy calculation only

      
      double I5=AnisoRestLocal[0]*AnisoRestLocal[0]*E2[0][0]   //atE2a used for energy calculation only
        +AnisoRestLocal[0]*AnisoRestLocal[1]*(E2[0][1]+E2[1][0])
        +AnisoRestLocal[1]*AnisoRestLocal[1]*E2[1][1];
      
    

      double StrainAlmansi[2][2]; // e=0.5(1-B^-1)  True strain tensor
      temp=LeftCauchy[0][0]*LeftCauchy[1][1]-LeftCauchy[1][0]*LeftCauchy[0][1]; // det(B)
      StrainAlmansi[0][0]=0.5*(1-(LeftCauchy[1][1]/temp));
      StrainAlmansi[1][0]=0.5*LeftCauchy[1][0]/temp;
      StrainAlmansi[0][1]=0.5*LeftCauchy[0][1]/temp;  
      StrainAlmansi[1][1]=0.5*(1-(LeftCauchy[0][0]/temp));
      




      double atEa=AnisoRestLocal[0]*AnisoRestLocal[0]*Egreen[0][0]
        +AnisoRestLocal[0]*AnisoRestLocal[1]*(Egreen[0][1]+Egreen[1][0])
        +AnisoRestLocal[1]*AnisoRestLocal[1]*Egreen[1][1];

      double I4=atEa;


    

  
      double Eaa[2][2];
      Eaa[0][0]= Egreen[0][0]*directAniso[0][0]+Egreen[0][1]*directAniso[1][0];        
      Eaa[1][0]= Egreen[1][0]*directAniso[0][0]+Egreen[1][1]*directAniso[1][0];        
      Eaa[0][1]= Egreen[0][0]*directAniso[0][1]+Egreen[0][1]*directAniso[1][1];        
      Eaa[1][1]= Egreen[1][0]*directAniso[0][1]+Egreen[1][1]*directAniso[1][1];        

      double aaE[2][2];
      aaE[0][0]= directAniso[0][0]*Egreen[0][0]+directAniso[0][1]*Egreen[1][0];        
      aaE[1][0]= directAniso[1][0]*Egreen[0][0]+directAniso[1][1]*Egreen[1][0];        
      aaE[0][1]= directAniso[0][0]*Egreen[0][1]+directAniso[0][1]*Egreen[1][1];        
      aaE[1][1]= directAniso[1][0]*Egreen[0][1]+directAniso[1][1]*Egreen[1][1];        
      
      double B2[2][2];// LeftCauchy^2
      B2[0][0]=LeftCauchy[0][0]*LeftCauchy[0][0]+LeftCauchy[0][1]*LeftCauchy[1][0];
      B2[1][0]=LeftCauchy[1][0]*LeftCauchy[0][0]+LeftCauchy[1][1]*LeftCauchy[1][0];
      B2[0][1]=LeftCauchy[0][0]*LeftCauchy[0][1]+LeftCauchy[0][1]*LeftCauchy[1][1];
      B2[1][1]=LeftCauchy[1][0]*LeftCauchy[0][1]+LeftCauchy[1][1]*LeftCauchy[1][1];

      double areaFactor=restingArea/Area; // 1/detF
      //double areaFactor=Area/restingArea; // detF
      
      double Sigma[2][2]; // true stress tensor (isotropic term) 
      Sigma[0][0]=areaFactor*((lambdaT*trE-mioT)*LeftCauchy[0][0]+(mioT)*B2[0][0]);
      Sigma[1][0]=areaFactor*((lambdaT*trE-mioT)*LeftCauchy[1][0]+(mioT)*B2[1][0]);
      Sigma[0][1]=areaFactor*((lambdaT*trE-mioT)*LeftCauchy[0][1]+(mioT)*B2[0][1]);
      Sigma[1][1]=areaFactor*((lambdaT*trE-mioT)*LeftCauchy[1][1]+(mioT)*B2[1][1]);
       
      

      // ENERGY MODELS  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


      // double deltaS[2][2]; // based on Delin. paper
      // deltaS[0][0]=deltaLam*(trE*directAniso[0][0]+atEa)+(2*deltaMio)*(Eaa[0][0]+aaE[0][0])-(deltaLam+2*deltaMio)*atEa*directAniso[0][0];
      // deltaS[1][0]=deltaLam*(trE*directAniso[1][0]     )+(2*deltaMio)*(Eaa[1][0]+aaE[1][0])-(deltaLam+2*deltaMio)*atEa*directAniso[1][0];
      // deltaS[0][1]=deltaLam*(trE*directAniso[0][1]     )+(2*deltaMio)*(Eaa[0][1]+aaE[0][1])-(deltaLam+2*deltaMio)*atEa*directAniso[0][1];
      // deltaS[1][1]=deltaLam*(trE*directAniso[1][1]+atEa)+(2*deltaMio)*(Eaa[1][1]+aaE[1][1])-(deltaLam+2*deltaMio)*atEa*directAniso[1][1];

      double deltaS[2][2]; // based on  equipartitioning
      deltaS[0][0]=(deltaLam/2)*(trE*directAniso[0][0]+atEa)+(deltaMio)*(Eaa[0][0]+aaE[0][0]);
      deltaS[1][0]=(deltaLam/2)*(trE*directAniso[1][0]     )+(deltaMio)*(Eaa[1][0]+aaE[1][0]);
      deltaS[0][1]=(deltaLam/2)*(trE*directAniso[0][1]     )+(deltaMio)*(Eaa[0][1]+aaE[0][1]);
      deltaS[1][1]=(deltaLam/2)*(trE*directAniso[1][1]+atEa)+(deltaMio)*(Eaa[1][1]+aaE[1][1]);


      // double deltaMio1=(youngL-youngT)/2;
      // double deltaS[2][2]; // based on  ... 
      // deltaS[0][0]=2*(deltaMio1)* atEa * directAniso[0][0];
      // deltaS[1][0]=2*(deltaMio1)* atEa * directAniso[1][0];
      // deltaS[0][1]=2*(deltaMio1)* atEa * directAniso[0][1];
      // deltaS[1][1]=2*(deltaMio1)* atEa * directAniso[1][1];




      strainZ +=restingArea*(1-poissonT*((2*lambdaT*trE+2*mioT*trE)+deltaS[0][0]+deltaS[1][1])/youngT);
      //std::cerr<< "cell "<< cellIndex<< " thickness :  " << strainZ << std::endl;
      
      //<<<<<<<<<<<<<<<<<<<isotropic force from stress tensor <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      // double ss[2][2];//lambda(trE)I+2mioE
      // ss[0][0]=lambdaT*trE+2*mioT*Egreen[0][0];
      // ss[0][1]=            2*mioT*Egreen[0][1];
      // ss[1][0]=            2*mioT*Egreen[1][0];
      // ss[1][1]=lambdaT*trE+2*mioT*Egreen[1][1];

       

      // double TPK[2][2];// 2nd Piola Kirchhoff stress tensor 
      // TPK[0][0]=restingArea*(DeformGrad[0][0]*ss[0][0]+DeformGrad[0][1]*ss[1][0]);
      // TPK[1][0]=restingArea*(DeformGrad[1][0]*ss[0][0]+DeformGrad[1][1]*ss[1][0]);
      // TPK[0][1]=restingArea*(DeformGrad[0][0]*ss[0][1]+DeformGrad[0][1]*ss[1][1]);
      // TPK[1][1]=restingArea*(DeformGrad[1][0]*ss[0][1]+DeformGrad[1][1]*ss[1][1]);

      // //deltaFTPKlocal[i][0]= TPK[0][0]*ShapeVectorResting[i][0]+TPK[0][1]*ShapeVectorResting[i][1];
      // //deltaFTPKlocal[i][1]= TPK[1][0]*ShapeVectorResting[i][0]+TPK[1][1]*ShapeVectorResting[i][1];
     
      // double deltaFTPKlocal[2][2];
      // deltaFTPKlocal[0][0]= TPK[0][0]*ShapeVectorResting[0][0]+TPK[0][1]*ShapeVectorResting[0][1];
      // deltaFTPKlocal[0][1]= TPK[1][0]*ShapeVectorResting[0][0]+TPK[1][1]*ShapeVectorResting[0][1];
     
      // double deltaFTPK[2][2]; 
      // deltaFTPK[0][0]= rotation[0][0]*deltaFTPKlocal[0][0]+rotation[0][1]*deltaFTPKlocal[0][1];
      // deltaFTPK[0][1]= rotation[1][0]*deltaFTPKlocal[0][0]+rotation[1][1]*deltaFTPKlocal[0][1];
      // deltaFTPK[0][2]= rotation[2][0]*deltaFTPKlocal[0][0]+rotation[2][1]*deltaFTPKlocal[0][1];
      //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      double TPK[2][2];// 2nd Piola Kirchhoff stress tensor 
      TPK[0][0]=restingArea*(DeformGrad[0][0]*deltaS[0][0]+DeformGrad[0][1]*deltaS[1][0]);
      TPK[1][0]=restingArea*(DeformGrad[1][0]*deltaS[0][0]+DeformGrad[1][1]*deltaS[1][0]);
      TPK[0][1]=restingArea*(DeformGrad[0][0]*deltaS[0][1]+DeformGrad[0][1]*deltaS[1][1]);
      TPK[1][1]=restingArea*(DeformGrad[1][0]*deltaS[0][1]+DeformGrad[1][1]*deltaS[1][1]);

      //deltaFTPKlocal[i][0]= TPK[0][0]*ShapeVectorResting[i][0]+TPK[0][1]*ShapeVectorResting[i][1];
      //deltaFTPKlocal[i][1]= TPK[1][0]*ShapeVectorResting[i][0]+TPK[1][1]*ShapeVectorResting[i][1];
     
      double deltaFTPKlocal[3][2];
      deltaFTPKlocal[0][0]= TPK[0][0]*ShapeVectorResting[0][0]+TPK[0][1]*ShapeVectorResting[0][1];
      deltaFTPKlocal[0][1]= TPK[1][0]*ShapeVectorResting[0][0]+TPK[1][1]*ShapeVectorResting[0][1];
     
      deltaFTPKlocal[1][0]= TPK[0][0]*ShapeVectorResting[1][0]+TPK[0][1]*ShapeVectorResting[1][1];
      deltaFTPKlocal[1][1]= TPK[1][0]*ShapeVectorResting[1][0]+TPK[1][1]*ShapeVectorResting[1][1];

      deltaFTPKlocal[2][0]= TPK[0][0]*ShapeVectorResting[2][0]+TPK[0][1]*ShapeVectorResting[2][1];
      deltaFTPKlocal[2][1]= TPK[1][0]*ShapeVectorResting[2][0]+TPK[1][1]*ShapeVectorResting[2][1];


      double deltaFTPK[3][3]; 
      deltaFTPK[0][0]= rotation[0][0]*deltaFTPKlocal[0][0]+rotation[0][1]*deltaFTPKlocal[0][1];
      deltaFTPK[0][1]= rotation[1][0]*deltaFTPKlocal[0][0]+rotation[1][1]*deltaFTPKlocal[0][1];
      deltaFTPK[0][2]= rotation[2][0]*deltaFTPKlocal[0][0]+rotation[2][1]*deltaFTPKlocal[0][1];

      deltaFTPK[1][0]= rotation[0][0]*deltaFTPKlocal[1][0]+rotation[0][1]*deltaFTPKlocal[1][1];
      deltaFTPK[1][1]= rotation[1][0]*deltaFTPKlocal[1][0]+rotation[1][1]*deltaFTPKlocal[1][1];
      deltaFTPK[1][2]= rotation[2][0]*deltaFTPKlocal[1][0]+rotation[2][1]*deltaFTPKlocal[1][1];
      
      deltaFTPK[2][0]= rotation[0][0]*deltaFTPKlocal[2][0]+rotation[0][1]*deltaFTPKlocal[2][1];
      deltaFTPK[2][1]= rotation[1][0]*deltaFTPKlocal[2][0]+rotation[1][1]*deltaFTPKlocal[2][1];
      deltaFTPK[2][2]= rotation[2][0]*deltaFTPKlocal[2][0]+rotation[2][1]*deltaFTPKlocal[2][1];


      // //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



      double deltaSFt[2][2];
      deltaSFt[0][0]=deltaS[0][0]*DeformGrad[0][0]+deltaS[0][1]*DeformGrad[0][1];
      deltaSFt[1][0]=deltaS[1][0]*DeformGrad[0][0]+deltaS[1][1]*DeformGrad[0][1];
      deltaSFt[0][1]=deltaS[0][0]*DeformGrad[1][0]+deltaS[0][1]*DeformGrad[1][1];
      deltaSFt[1][1]=deltaS[1][0]*DeformGrad[1][0]+deltaS[1][1]*DeformGrad[1][1];
      
      double deltaSigma[2][2];// true stress tensor (anisotropic correction term)deltaLambda and deltaMio (Longitudinal-Transverse)
      deltaSigma[0][0]=areaFactor*(DeformGrad[0][0]*deltaSFt[0][0]+DeformGrad[0][1]*deltaSFt[1][0]);
      deltaSigma[1][0]=areaFactor*(DeformGrad[1][0]*deltaSFt[0][0]+DeformGrad[1][1]*deltaSFt[1][0]);
      deltaSigma[0][1]=areaFactor*(DeformGrad[0][0]*deltaSFt[0][1]+DeformGrad[0][1]*deltaSFt[1][1]);
      deltaSigma[1][1]=areaFactor*(DeformGrad[1][0]*deltaSFt[0][1]+DeformGrad[1][1]*deltaSFt[1][1]);

      double sigmafactor=1;
      double deltasigmafactor=1;

      double StressTensor[3][3];
      StressTensor[0][0]=sigmafactor*Sigma[0][0]+deltasigmafactor*deltaSigma[0][0];
      StressTensor[1][0]=sigmafactor*Sigma[1][0]+deltasigmafactor*deltaSigma[1][0];
      StressTensor[0][1]=sigmafactor*Sigma[0][1]+deltasigmafactor*deltaSigma[0][1];
      StressTensor[1][1]=sigmafactor*Sigma[1][1]+deltasigmafactor*deltaSigma[1][1];


      


      //Shape vectors in Current shape (counterclockwise ordering of nodes/edges)     ShapeVectorCurrent[3][3]  calculated above   
      //.............................. ( or clockwise ordering of nodes/edges)
          

      //square of radius of circumstancing circle in resting shape
      //double Rcirc2Resting=(0.25*restingLength[0]*restingLength[1]*restingLength[2]/Area)*(0.25*restingLength[0]*restingLength[1]*restingLength[2]/Area);  
      
      double StrainTensor[3][3]; // there are other alternatives than StrainAlmansi for strain tensor
      StrainTensor[0][0]=StrainAlmansi[0][0];
      StrainTensor[1][0]=StrainAlmansi[1][0];
      StrainTensor[0][1]=StrainAlmansi[0][1];
      StrainTensor[1][1]=StrainAlmansi[1][1];


      StrainTensor[0][2]=0;  // adding 3rd dimension which is zero, the tensor is still in element plane
      StrainTensor[1][2]=0;
      StrainTensor[2][2]=0;
      StrainTensor[2][0]=0;
      StrainTensor[2][1]=0;
      
      StressTensor[0][2]=0;  // adding 3rd dimension which is zero, the tensor is still in element plane
      StressTensor[1][2]=0;
      StressTensor[2][2]=0;
      StressTensor[2][0]=0;
      StressTensor[2][1]=0;

      //rotation matrix to go to global coordinate system based on counterclockwise ordering;   rotation[3][3] calculated above  

      // rotating strain tensor to the global coordinate system
      double tempR[3][3]={{0,0,0},{0,0,0},{0,0,0}};
      for (int r=0 ; r<3 ; r++) 
        for (int s=0 ; s<3 ; s++) 
          for(int w=0 ; w<3 ; w++) 
            tempR[r][s] += rotation[r][w]*StrainTensor[w][s];
          
      for (int r=0 ; r<3 ; r++) 
        for (int s=0 ; s<3 ; s++) {
          StrainTensor[r][s]=0;
          for(int w=0 ; w<3 ; w++) 
            StrainTensor[r][s] += tempR[r][w]*rotation[s][w]; 
        }

      // rotating stress tensor to the global coordinate system
      for (int r=0 ; r<3 ; r++) 
        for (int s=0 ; s<3 ; s++) 
          tempR[r][s]=0;
        
      for (int r=0 ; r<3 ; r++) 
        for (int s=0 ; s<3 ; s++) 
          for(int w=0 ; w<3 ; w++) 
            tempR[r][s] += rotation[r][w]*StressTensor[w][s];
          
      for (int r=0 ; r<3 ; r++) 
        for (int s=0 ; s<3 ; s++) {
          StressTensor[r][s]=0;
          for(int w=0 ; w<3 ; w++) 
            StressTensor[r][s] += tempR[r][w]*rotation[s][w]; 
        }
      
 


           
     
      // accumulating strain and stress tensors and normal to cell plane vector to be averaged later
      for (int r=0 ; r<3 ; r++) 
        for (int s=0 ; s<3 ; s++)    
          StrainCellGlobal[r][s] += Area*StrainTensor[r][s];

      for (int r=0 ; r<3 ; r++) 
        for (int s=0 ; s<3 ; s++)    
          StressCellGlobal[r][s] += Area*StressTensor[r][s];
      
      for (int r=0 ; r<3 ; r++) 
        normalGlob[r] += Area*Zcurrent[r];

      TotalCellRestingArea=TotalCellRestingArea+restingArea;
      TotalCellArea=TotalCellArea+Area;
      // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> STRAIN and STRESS TENSORS (END) <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

   
    //---- Anisotropic Correction Force-------------------------------
       double deltaF[3][3];




        for ( int i=0 ; i<3 ; ++i )  // from stress tensor(equipartitioning energy)
          for ( int j=0 ; j<3 ; ++j )
            deltaF[i][j]=(-deltaFTPK[i][j]);
      
        double  I1=trE;

        // energy
        EnergyIso +=( (lambdaT/2)*I1*I1 + mioT*I2 )*restingArea; //<<<<<<<<<<<<<<<<<<<<<
      
        EnergyAniso +=( (deltaLam/2)*I4*I1 + deltaMio*I5 )*restingArea; //<<<<<<<<<<<<<<<<<<<<<<<<<<<
        
        
        
        //Forces of vertices   
        double Force[3][3]={{0,0,0},{0,0,0},{0,0,0}};                                          
        

        Force[0][0]= (tensileStiffness[0]*Delta[0]+angularStiffness[1]*Delta[1]+angularStiffness[0]*Delta[2])*(position[1][0]-position[0][0])
          +(tensileStiffness[2]*Delta[2]+angularStiffness[2]*Delta[1]+angularStiffness[0]*Delta[0])*(position[2][0]-position[0][0])
          + deltaF[0][0]; 
        Force[0][1]= (tensileStiffness[0]*Delta[0]+angularStiffness[1]*Delta[1]+angularStiffness[0]*Delta[2])*(position[1][1]-position[0][1])
          +(tensileStiffness[2]*Delta[2]+angularStiffness[2]*Delta[1]+angularStiffness[0]*Delta[0])*(position[2][1]-position[0][1])
          + deltaF[0][1];  
        Force[0][2]= (tensileStiffness[0]*Delta[0]+angularStiffness[1]*Delta[1]+angularStiffness[0]*Delta[2])*(position[1][2]-position[0][2])
          +(tensileStiffness[2]*Delta[2]+angularStiffness[2]*Delta[1]+angularStiffness[0]*Delta[0])*(position[2][2]-position[0][2])
          + deltaF[0][2]; 

        Force[1][0]= (tensileStiffness[0]*Delta[0]+angularStiffness[0]*Delta[2]+angularStiffness[1]*Delta[1])*(position[0][0]-position[1][0])
          +(tensileStiffness[1]*Delta[1]+angularStiffness[2]*Delta[2]+angularStiffness[1]*Delta[0])*(position[2][0]-position[1][0])
          + deltaF[1][0];  
        Force[1][1]= (tensileStiffness[0]*Delta[0]+angularStiffness[0]*Delta[2]+angularStiffness[1]*Delta[1])*(position[0][1]-position[1][1])
          +(tensileStiffness[1]*Delta[1]+angularStiffness[2]*Delta[2]+angularStiffness[1]*Delta[0])*(position[2][1]-position[1][1])
          + deltaF[1][1];  
        Force[1][2]= (tensileStiffness[0]*Delta[0]+angularStiffness[0]*Delta[2]+angularStiffness[1]*Delta[1])*(position[0][2]-position[1][2])
          +(tensileStiffness[1]*Delta[1]+angularStiffness[2]*Delta[2]+angularStiffness[1]*Delta[0])*(position[2][2]-position[1][2])
          + deltaF[1][2];  
        
        Force[2][0]= (tensileStiffness[2]*Delta[2]+angularStiffness[0]*Delta[0]+angularStiffness[2]*Delta[1])*(position[0][0]-position[2][0])
          +(tensileStiffness[1]*Delta[1]+angularStiffness[1]*Delta[0]+angularStiffness[2]*Delta[2])*(position[1][0]-position[2][0])
          + deltaF[2][0];  
        Force[2][1]= (tensileStiffness[2]*Delta[2]+angularStiffness[0]*Delta[0]+angularStiffness[2]*Delta[1])*(position[0][1]-position[2][1])
          +(tensileStiffness[1]*Delta[1]+angularStiffness[1]*Delta[0]+angularStiffness[2]*Delta[2])*(position[1][1]-position[2][1])
          + deltaF[2][1];  
        Force[2][2]= (tensileStiffness[2]*Delta[2]+angularStiffness[0]*Delta[0]+angularStiffness[2]*Delta[1])*(position[0][2]-position[2][2])
          +(tensileStiffness[1]*Delta[1]+angularStiffness[1]*Delta[0]+angularStiffness[2]*Delta[2])*(position[1][2]-position[2][2])
          + deltaF[2][2];
        
        
        bool isSliver=false;
        if (Angle[0]<0.1 || Angle[0]<0.1 ||   Angle[0]<0.1) {
          double tmp1=cellDerivs[cellIndex][comIndex  ]+vertexDerivs[v2][0]+vertexDerivs[v3][0],
            tmp2=cellDerivs[cellIndex][comIndex+1]+vertexDerivs[v2][1]+vertexDerivs[v3][1],
            tmp3=cellDerivs[cellIndex][comIndex+2]+vertexDerivs[v2][2]+vertexDerivs[v3][2];

          cellDerivs[cellIndex][comIndex  ] +=tmp1;
          cellDerivs[cellIndex][comIndex+1] +=tmp3;
          cellDerivs[cellIndex][comIndex+2] +=tmp3;
          
          vertexDerivs[v2][0] += tmp1;
          vertexDerivs[v2][1] += tmp2;
          vertexDerivs[v2][2] += tmp3;
          
          vertexDerivs[v3][0] += tmp1;
          vertexDerivs[v3][1] += tmp2;
          vertexDerivs[v3][2] += tmp3;

          isSliver=true;
        }  
        



        // if (cellIndex==88){
        //   std::cerr<<std::endl;                
        //   std::cerr<<std::endl;                
        //   std::cerr<<std::endl;      
        //   std::cerr<<"  position 0 "<<position[0][0]<<"  "<<position[0][1]<<"  "<<position[0][2]<<std::endl;                
        //   std::cerr<<"  position 1 "<<position[1][0]<<"  "<<position[1][1]<<"  "<<position[1][2]<<std::endl;                
        //   std::cerr<<"  position 2 "<<position[2][0]<<"  "<<position[2][1]<<"  "<<position[2][2]<<std::endl;                
        //   std::cerr<<std::endl;                
        //   std::cerr<<"  tensileStiffness "<<tensileStiffness[0]<<"  "<<tensileStiffness[1]<<"  "<<tensileStiffness[2]<<std::endl;                
        //   std::cerr<<"  angularStiffness "<<angularStiffness[0]<<"  "<<angularStiffness[1]<<"  "<<angularStiffness[2]<<std::endl;                
        //   std::cerr<<"  Delta2           "<<Delta[0]<<"  "<<Delta[1]<<"  "<<Delta[2]<<std::endl;                
        //   std::cerr<<std::endl;                
        //   std::cerr<<" force 0 "<<Force[0][0]<<"  "<<Force[0][1]<<"  "<<Force[0][2]<<std::endl;                
        //   std::cerr<<" force 1 "<<Force[1][0]<<"  "<<Force[1][1]<<"  "<<Force[1][2]<<std::endl;                
        //   std::cerr<<" force 2 "<<Force[2][0]<<"  "<<Force[2][1]<<"  "<<Force[2][2]<<std::endl;              
        //   std::cerr<<std::endl;                
        //   std::cerr<<"  test   "<<(tensileStiffness[2]*0.4+angularStiffness[0]*0.9+angularStiffness[2]*0.1)*(position[0][2]-position[2][2])
        //     +(tensileStiffness[1]*0.1+angularStiffness[1]*0.9+angularStiffness[2]*0.4)*(position[1][2]-position[2][2])
        //            <<std::endl;                
          
        // }
        
        // adding TRBSMT forces to the total vertexDerives
        if (! isSliver) {
          cellDerivs[cellIndex][comIndex  ] += Force[0][0];
          cellDerivs[cellIndex][comIndex+1] += Force[0][1];
          cellDerivs[cellIndex][comIndex+2] += Force[0][2];
          
          vertexDerivs[v2][0] += Force[1][0];
          vertexDerivs[v2][1] += Force[1][1];
          vertexDerivs[v2][2] += Force[1][2];
          
          vertexDerivs[v3][0] += Force[2][0];
          vertexDerivs[v3][1] += Force[2][1];
          vertexDerivs[v3][2] += Force[2][2];
        }
        
    }
    

    for (int r=0 ; r<3 ; r++) 
      for (int s=0 ; s<3 ; s++)    
        StrainCellGlobal[r][s]= StrainCellGlobal[r][s]/TotalCellArea; 
    
    
    
    for (int r=0 ; r<3 ; r++) 
      for (int s=0 ; s<3 ; s++)    
        StressCellGlobal[r][s]= StressCellGlobal[r][s]/TotalCellArea; 

    for (int r=0 ; r<3 ; r++)   
      normalGlob[r]/= TotalCellArea;
    
    double temp=std::sqrt(normalGlob[0]*normalGlob[0]+normalGlob[1]*normalGlob[1]+normalGlob[2]*normalGlob[2]);
    
    if(temp>0)
      for (int r=0 ; r<3 ; r++)   
        normalGlob[r]/=temp;
   
    double areaRatio=TotalCellArea/ TotalCellRestingArea; 
    
    strainZ=strainZ/TotalCellRestingArea; 
  
   
    // eigenvalue/eigenvectors of averaged STRAIN and STRESS tensors in global coordinate system. (Jacobi method)

    cellData[cellIndex][stressTensorIndex  ]=StressCellGlobal[0][0];
    cellData[cellIndex][stressTensorIndex+1]=StressCellGlobal[1][1];
    cellData[cellIndex][stressTensorIndex+2]=StressCellGlobal[2][2];
    cellData[cellIndex][stressTensorIndex+3]=StressCellGlobal[0][1];
    cellData[cellIndex][stressTensorIndex+4]=StressCellGlobal[0][2];
    cellData[cellIndex][stressTensorIndex+5]=StressCellGlobal[1][2];
    
    cellData[cellIndex][normalVectorIndex  ]=normalGlob[0];
    cellData[cellIndex][normalVectorIndex+1]=normalGlob[1];
    cellData[cellIndex][normalVectorIndex+2]=normalGlob[2];
   
   
  
    // STRAIN:
    double strainEpcilon =0.000001;

    double eigenVectorStrain[3][3]={{1,0,0},{0,1,0},{0,0,1}};
    double pivot=1;
    double pi=3.1415;
    int I,J;
    double RotAngle,Si,Co;
    


    pivot=std::abs(StrainCellGlobal[1][0]);
    I=1;
    J=0;
    if (std::abs(StrainCellGlobal[2][0])>pivot) {
      pivot=std::abs(StrainCellGlobal[2][0]);
      I=2;
      J=0;
    }
    if (std::abs(StrainCellGlobal[2][1])>pivot) {
      pivot=std::abs(StrainCellGlobal[2][1]);
      I=2;
      J=1;
    }


    while (pivot> strainEpcilon) {
      
      if (std::abs(StrainCellGlobal[I][I]-StrainCellGlobal[J][J])< strainEpcilon ) {
          RotAngle=pi/4;
      }            
      else {
        RotAngle=0.5*std::atan((2*StrainCellGlobal[I][J])/(StrainCellGlobal[J][J]-StrainCellGlobal[I][I]));
      }
        Si=std::sin(RotAngle);
        Co=std::cos(RotAngle);
        double tempRot[3][3]={{1,0,0},{0,1,0},{0,0,1}};
        tempRot[I][I]=Co;
        tempRot[J][J]=Co;
        tempRot[I][J]=Si;
        tempRot[J][I]=-Si;
        double tempStrain[3][3]={{0,0,0},{0,0,0},{0,0,0}};
        for (int r=0 ; r<3 ; r++) 
          for (int s=0 ; s<3 ; s++) 
            for(int w=0 ; w<3 ; w++) 
              tempStrain[r][s]=tempStrain[r][s]+StrainCellGlobal[r][w]*tempRot[w][s];
        
        for (int r=0 ; r<3 ; r++) 
          for (int s=0 ; s<3 ; s++) 
            StrainCellGlobal[r][s]=0;
         
        for (int r=0 ; r<3 ; r++) 
          for (int s=0 ; s<3 ; s++) 
            for(int w=0 ; w<3 ; w++) 
              StrainCellGlobal[r][s]=StrainCellGlobal[r][s]+tempRot[w][r]*tempStrain[w][s];
                
        for (int r=0 ; r<3 ; r++) 
          for (int s=0 ; s<3 ; s++) 
            tempStrain[r][s]=eigenVectorStrain[r][s];
        
        for (int r=0 ; r<3 ; r++) 
          for (int s=0 ; s<3 ; s++) 
            eigenVectorStrain[r][s]=0;
        
        for (int r=0 ; r<3 ; r++) 
          for (int s=0 ; s<3 ; s++) 
            for(int w=0 ; w<3 ; w++) 
              eigenVectorStrain[r][s]=eigenVectorStrain[r][s]+tempStrain[r][w]*tempRot[w][s];
            
        
        pivot=std::abs(StrainCellGlobal[1][0]);
        I=1;
        J=0;
        if (std::abs(StrainCellGlobal[2][0])>pivot) {
          pivot=std::abs(StrainCellGlobal[2][0]);
          I=2;
          J=0;
        }
        if (std::abs(StrainCellGlobal[2][1])>pivot) {
          pivot=std::abs(StrainCellGlobal[2][1]);
          I=2;
          J=1;
        }
    }
    
    // normalizing eigenvectors (remove if not necessary)  
    temp=std::sqrt(eigenVectorStrain[0][0]*eigenVectorStrain[0][0] +
                   eigenVectorStrain[1][0]*eigenVectorStrain[1][0] +
                   eigenVectorStrain[2][0]*eigenVectorStrain[2][0] );
    if(temp>0){
      eigenVectorStrain[0][0]/=temp;
      eigenVectorStrain[1][0]/=temp;
      eigenVectorStrain[2][0]/=temp;
    }
    temp=std::sqrt(eigenVectorStrain[0][1]*eigenVectorStrain[0][1] +
                   eigenVectorStrain[1][1]*eigenVectorStrain[1][1] +
                   eigenVectorStrain[2][1]*eigenVectorStrain[2][1] );
    if(temp>0){
      eigenVectorStrain[0][1]/=temp;
      eigenVectorStrain[1][1]/=temp;
      eigenVectorStrain[2][1]/=temp;
    }
    temp=std::sqrt(eigenVectorStrain[0][2]*eigenVectorStrain[0][2] +
                   eigenVectorStrain[1][2]*eigenVectorStrain[1][2] +
                   eigenVectorStrain[2][2]*eigenVectorStrain[2][2] );
    if(temp>0){
      eigenVectorStrain[0][2]/=temp;
      eigenVectorStrain[1][2]/=temp;
      eigenVectorStrain[2][2]/=temp;
    }
 
      // maximal strain direction
      double maximalStrainValue=StrainCellGlobal[0][0];
      int Istrain=0;
      if (StrainCellGlobal[1][1]>maximalStrainValue) 
        {
          maximalStrainValue=StrainCellGlobal[1][1];
          Istrain=1;
        }
      if (StrainCellGlobal[2][2]>maximalStrainValue) 
        {
          maximalStrainValue=StrainCellGlobal[2][2];
          Istrain=2;
        }
   

      // 2nd maximalstrain direction/value
      double maximalStrainValue2,maximalStrainValue3;
      int Istrain2,Istrain3;
      if (Istrain==0) {
        Istrain2=1;
        Istrain3=2;
      }
      if (Istrain==1) {
        Istrain2=0;
        Istrain3=2;
      }
      if (Istrain==2) {
        Istrain2=0;
        Istrain3=1;
      }
      if(StrainCellGlobal[Istrain3][Istrain3] > StrainCellGlobal[Istrain2][Istrain2] ) {
        temp=Istrain2;
        Istrain2=Istrain3;
        Istrain3=temp;
      }
      maximalStrainValue2=StrainCellGlobal[Istrain2][Istrain2];
      maximalStrainValue3=StrainCellGlobal[Istrain3][Istrain3];
      
      
      if(std::abs(normalGlob[0]*eigenVectorStrain[0][Istrain]+   //<<<<<<<<<<<<<<<<<<
                  normalGlob[1]*eigenVectorStrain[1][Istrain]+
                  normalGlob[2]*eigenVectorStrain[2][Istrain]) > .2) {

	std::cerr << "max strain vector is out of cell plane in the cell " <<cellIndex <<std::endl;
        // Istrain=Istrain2; 
        // Istrain2=Istrain3;
        // maximalStrainValue=maximalStrainValue2;
        // maximalStrainValue2=maximalStrainValue3;
        
        // for (size_t wallindex=0; wallindex<numWalls; ++wallindex){
          
        //   size_t kPlusOneMod = (wallindex+1)%numWalls;

        //   //size_t v1 = com;
        //   size_t v2 = T.cell(cellIndex).vertex(wallindex)->index();
        //   size_t v3 = T.cell(cellIndex).vertex(kPlusOneMod)->index();
        //   //size_t w1 = internal wallindex
        //   size_t w2 = T.cell(cellIndex).wall(wallindex)->index();
        //   //size_t w3 = internal wallindex+1
          
        //   // Position matrix holds in rows positions for com, vertex(wallindex), vertex(wallindex+1)
        //   DataMatrix position(3,vertexData[v2]);
        //   for (size_t d=0; d<dimension; ++d)
        //     position[0][d] = cellData[cellIndex][comIndex+d]; // com position
        //   //position[1] = vertexData[v2]; // given by initiation
        //   position[2] = vertexData[v3];
        //   //position[0][2] z for vertex 1 of the current element
          
        //   std::vector<double> restingLength(3);
        //   restingLength[0] = cellData[cellIndex][lengthInternalIndex + wallindex];
        //   restingLength[1] = wallData[w2][wallLengthIndex];
        //   restingLength[2] = cellData[cellIndex][lengthInternalIndex + kPlusOneMod];
          
          
        //   std::vector<double> length(3);
        //   length[0] = std::sqrt( (position[0][0]-position[1][0])*(position[0][0]-position[1][0]) +
        //                          (position[0][1]-position[1][1])*(position[0][1]-position[1][1]) +
        //                          (position[0][2]-position[1][2])*(position[0][2]-position[1][2]) );
          
        //   length[1] = T.wall(w2).lengthFromVertexPosition(vertexData);
          
        //   length[2] = std::sqrt( (position[0][0]-position[2][0])*(position[0][0]-position[2][0]) +
        //                          (position[0][1]-position[2][1])*(position[0][1]-position[2][1]) +
        //                          (position[0][2]-position[2][2])*(position[0][2]-position[2][2]) );
        //   std::cerr<<" resting length 1,2,3 :   "
        //            <<restingLength[0]<<"  "<<restingLength[1]<<"  "<<restingLength[2]<<std::endl;
        //   std::cerr<<"         length 1,2,3 :   "
        //            <<length[0]<<"  "<<length[1]<<"  "<<length[2]<<std::endl;
        //   cellData[cellIndex][23]=1000;
        // }
        // std::cerr<<" strain 0,1,2 :   " 
       //            <<maximalStrainValue<<"  "<<maximalStrainValue2<<"  "<<maximalStrainValue3<<std::endl;
      }
      
      //if(maximalStrainValue != maximalStrainValue){
      //if(std::abs(cellData[cellIndex][comIndex+1]) >100 ||std::abs(cellData[cellIndex][comIndex]) >100){
      // if(cellIndex==88){
      //   for (size_t wallindex=0; wallindex<numWalls; ++wallindex){
          
      //     size_t kPlusOneMod = (wallindex+1)%numWalls;
          
      //     //size_t v1 = com;
      //     size_t v2 = T.cell(cellIndex).vertex(wallindex)->index();
      //     size_t v3 = T.cell(cellIndex).vertex(kPlusOneMod)->index();
      //     //size_t w1 = internal wallindex
      //     size_t w2 = T.cell(cellIndex).wall(wallindex)->index();
      //     //size_t w3 = internal wallindex+1
          
      //     // Position matrix holds in rows positions for com, vertex(wallindex), vertex(wallindex+1)
      //     DataMatrix position(3,vertexData[v2]);
      //     for (size_t d=0; d<dimension; ++d)
      //       position[0][d] = cellData[cellIndex][comIndex+d]; // com position
      //     //position[1] = vertexData[v2]; // given by initiation
      //     position[2] = vertexData[v3];
      //     //position[0][2] z for vertex 1 of the current element
          
      //     std::vector<double> restingLength(3);
      //     restingLength[0] = cellData[cellIndex][lengthInternalIndex + wallindex];
      //     restingLength[1] = wallData[w2][wallLengthIndex];
      //     restingLength[2] = cellData[cellIndex][lengthInternalIndex + kPlusOneMod];
          
          
      //     std::vector<double> length(3);
      //     length[0] = std::sqrt( (position[0][0]-position[1][0])*(position[0][0]-position[1][0]) +
      //                            (position[0][1]-position[1][1])*(position[0][1]-position[1][1]) +
      //                            (position[0][2]-position[1][2])*(position[0][2]-position[1][2]) );
          
      //     length[1] = T.wall(w2).lengthFromVertexPosition(vertexData);
          
      //     length[2] = std::sqrt( (position[0][0]-position[2][0])*(position[0][0]-position[2][0]) +
      //                            (position[0][1]-position[2][1])*(position[0][1]-position[2][1]) +
      //                            (position[0][2]-position[2][2])*(position[0][2]-position[2][2]) );
          
      //     std::cerr<<" resting length 1,2,3 :   "
      //              <<restingLength[0]<<"  "<<restingLength[1]<<"  "<<restingLength[2]<<std::endl;
      //     std::cerr<<"         length 1,2,3 :   "
      //              <<length[0]<<"  "<<length[1]<<"  "<<length[2]<<std::endl;
      //     cellData[cellIndex][23]=1000;
      //   }
      //   std::cerr<<" in the cell :   "<<cellIndex<<std::endl;
      //   std::cerr<<" strain 0,1,2 :   " 
      //             <<maximalStrainValue<<"  "<<maximalStrainValue2<<"  "<<maximalStrainValue3<<std::endl;
      // }
      
      
      if(maximalStrainValue != maximalStrainValue){
        
        std::cerr<<"mechanicalTRBS maximal strain "<<maximalStrainValue<<std::endl;
        std::cerr<<"in the cell "<<cellIndex<<std::endl;
        exit(-1);
      }


      // /////// ad-hoc begin
      // for (size_t wallindex=0; wallindex<numWalls; ++wallindex) { 
      //   size_t kPlusOneMod = (wallindex+1)%numWalls;
      //   //size_t v1 = com;
      //   size_t v2 = T.cell(cellIndex).vertex(wallindex)->index();
      //   size_t v3 = T.cell(cellIndex).vertex(kPlusOneMod)->index();
      //   //size_t w1 = internal wallindex
      //   size_t w2 = T.cell(cellIndex).wall(wallindex)->index();
      //   //size_t w3 = internal wallindex+1
        
      //   // Position matrix holds in rows positions for com, vertex(wallindex), vertex(wallindex+1)
      //   DataMatrix position(3,vertexData[v2]);
      //   for (size_t d=0; d<dimension; ++d)
      //     position[0][d] = cellData[cellIndex][comIndex+d]; // com position
      //   //position[1] = vertexData[v2]; // given by initiation
      //   position[2] = vertexData[v3];
      //   //position[0][2] z for vertex 1 of the current element
        
        
        
      //   std::vector<double> restingLength(numWalls);
      //   restingLength[0] = cellData[cellIndex][lengthInternalIndex + wallindex];
      //   restingLength[1] = wallData[w2][wallLengthIndex];
      //   restingLength[2] = cellData[cellIndex][lengthInternalIndex + kPlusOneMod];
        
      
      //   std::vector<double> length(numWalls);
      //   length[0] = std::sqrt( (position[0][0]-position[1][0])*(position[0][0]-position[1][0]) +
      //                          (position[0][1]-position[1][1])*(position[0][1]-position[1][1]) +
      //                          (position[0][2]-position[1][2])*(position[0][2]-position[1][2]) );
        
      //   length[1] = T.wall(w2).lengthFromVertexPosition(vertexData);
        
      //   length[2] = std::sqrt( (position[0][0]-position[2][0])*(position[0][0]-position[2][0]) +
      //                          (position[0][1]-position[2][1])*(position[0][1]-position[2][1]) +
      //                          (position[0][2]-position[2][2])*(position[0][2]-position[2][2]) );
      //   std:: vector<std::vector<double> > edges(3);
      //   for (size_t s=0; s<3 ; s++)
      //     edges[s].resize(3);
      //   for (size_t d=0; d<dimension; ++d){
      //     edges[0][d]=position[1][d]-position[0][d];
      //     edges[1][d]=position[2][d]-position[1][d];
      //     edges[2][d]=position[0][d]-position[2][d];
      //   }
      //   //go 3 dimensional
      //   std::vector<double> restingFromStrain(3);
      //   for (size_t d=0; d<dimension; ++d)
      //     restingFromStrain[d]=std::sqrt(
      //                                    (edges[d][0]*eigenVectorStrain[0][Istrain]+
      //                                     edges[d][1]*eigenVectorStrain[1][Istrain]+
      //                                     edges[d][2]*eigenVectorStrain[2][Istrain])*
      //                                    (edges[d][0]*eigenVectorStrain[0][Istrain]+
      //                                     edges[d][1]*eigenVectorStrain[1][Istrain]+
      //                                     edges[d][2]*eigenVectorStrain[2][Istrain])*
      //                                    (1-maximalStrainValue)*(1-maximalStrainValue)+
                                         
      //                                    (edges[d][0]*eigenVectorStrain[0][Istrain2]+
      //                                     edges[d][1]*eigenVectorStrain[1][Istrain2]+
      //                                     edges[d][2]*eigenVectorStrain[2][Istrain2])*
      //                                    (edges[d][0]*eigenVectorStrain[0][Istrain2]+
      //                                     edges[d][1]*eigenVectorStrain[1][Istrain2]+
      //                                     edges[d][2]*eigenVectorStrain[2][Istrain2])*
      //                                    (1-maximalStrainValue2)*(1-maximalStrainValue2)+
                                         
      //                                    (edges[d][0]*eigenVectorStrain[0][Istrain3]+
      //                                     edges[d][1]*eigenVectorStrain[1][Istrain3]+
      //                                     edges[d][2]*eigenVectorStrain[2][Istrain3])*
      //                                    (edges[d][0]*eigenVectorStrain[0][Istrain3]+
      //                                     edges[d][1]*eigenVectorStrain[1][Istrain3]+
      //                                     edges[d][2]*eigenVectorStrain[2][Istrain3])*
      //                                    (1-maximalStrainValue3)*(1-maximalStrainValue3)
      //                                    );
      //   for (size_t d=0; d<dimension; ++d)
      //     std::cerr<<restingLength[d]<<"  "<<restingFromStrain[d]<<"  "<<length[d]<<std::endl;
      // }
      // /////// ad-hoc end
      





 
      double growthStrain=maximalStrainValue2;
      if ( cellData[cellIndex][MTindex  ]*eigenVectorStrain[0][Istrain] +
	   cellData[cellIndex][MTindex+1]*eigenVectorStrain[1][Istrain] +
	   cellData[cellIndex][MTindex+2]*eigenVectorStrain[2][Istrain] < 0.01  ){
      growthStrain=maximalStrainValue;
      }

      // normal to the cell plane in global coordinate is averaged Zcurrent[], vector product gives the perpendicular strain direction
      double PerpStrain[3];
      PerpStrain[0]=normalGlob[1]*eigenVectorStrain[2][Istrain]-normalGlob[2]*eigenVectorStrain[1][Istrain];
      PerpStrain[1]=normalGlob[2]*eigenVectorStrain[0][Istrain]-normalGlob[0]*eigenVectorStrain[2][Istrain];
      PerpStrain[2]=normalGlob[0]*eigenVectorStrain[1][Istrain]-normalGlob[1]*eigenVectorStrain[0][Istrain];
      temp=std::sqrt(PerpStrain[0]*PerpStrain[0]+PerpStrain[1]*PerpStrain[1]+PerpStrain[2]*PerpStrain[2]);     
      
       if(std::abs(temp)<0.0001){ // if maximal strain is normal to the cell plane storing strain direction instead as it should not be used
         PerpStrain[0]=eigenVectorStrain[0][Istrain];
         PerpStrain[1]=eigenVectorStrain[1][Istrain];
         PerpStrain[2]=eigenVectorStrain[2][Istrain];
       }

       
    // storing a measure for strain anisotropy in cell vector
    if (std::abs(maximalStrainValue) <  0.000001) cellData[cellIndex][strainAnIndex]=0;
    if (std::abs(maximalStrainValue) >= 0.000001) cellData[cellIndex][strainAnIndex]=1-std::abs(maximalStrainValue2/maximalStrainValue);// relative
    //if (std::abs(maximalStrainValue) >= 0.000001) cellData[cellIndex][strainAnIndex]=(maximalStrainValue-maximalStrainValue2)/0.08;// absolute

    // storing strain/stress direction/value in cellData
    if (numVariableIndexLevel()==4 && (numVariableIndex(2)==1 || numVariableIndex(2)==2 || numVariableIndex(2)==3) ) {// storing maximal strain
      if (dimension==2)
        {
          cellData[cellIndex][variableIndex(2,0)]  =eigenVectorStrain[0][Istrain];
          cellData[cellIndex][variableIndex(2,0)+1]=eigenVectorStrain[1][Istrain];
          cellData[cellIndex][variableIndex(2,0)+3]=maximalStrainValue;  //maximal Strain Value is stored after its eigenvector
        }
      if (dimension==3)
        {
          cellData[cellIndex][variableIndex(2,0)]  =eigenVectorStrain[0][Istrain];
          cellData[cellIndex][variableIndex(2,0)+1]=eigenVectorStrain[1][Istrain];
          cellData[cellIndex][variableIndex(2,0)+2]=eigenVectorStrain[2][Istrain];
          cellData[cellIndex][variableIndex(2,0)+3]=maximalStrainValue;  //maximal Strain Value is stored after its eigenvector
        }
    }
    //std::cerr<< maximalStrainValue<< std::endl;
    if (numVariableIndexLevel()==4 && numVariableIndex(2)==3) {//storing 2nd maximal strain
      if (dimension==2)
        {
          cellData[cellIndex][variableIndex(2,2)]  =eigenVectorStrain[0][Istrain2];
          cellData[cellIndex][variableIndex(2,2)+1]=eigenVectorStrain[1][Istrain2];
          cellData[cellIndex][variableIndex(2,2)+3]=maximalStrainValue2;  //2nd maximal Strain Value is stored after its eigenvector
        }
      if (dimension==3)
        {
          cellData[cellIndex][variableIndex(2,2)]  =eigenVectorStrain[0][Istrain2];
          cellData[cellIndex][variableIndex(2,2)+1]=eigenVectorStrain[1][Istrain2];
          cellData[cellIndex][variableIndex(2,2)+2]=eigenVectorStrain[2][Istrain2];
          cellData[cellIndex][variableIndex(2,2)+3]=maximalStrainValue2;  //2nd maximal Strain Value is stored after its eigenvector
        }
    }
    
    //cellData[cellIndex][20]=0;
    
    if (numVariableIndexLevel()==4 && ( numVariableIndex(2)==2 || numVariableIndex(2)==3) ) {//storing perpendicular to maximal strain
      if (dimension==2)
        {
          cellData[cellIndex][variableIndex(2,1)]  =PerpStrain[0];
          cellData[cellIndex][variableIndex(2,1)+1]=PerpStrain[1];
          cellData[cellIndex][variableIndex(2,1)+3]=maximalStrainValue;  //maximal Strain Value is stored after its eigenvector
          
        }
      if (dimension==3)
        { 
          cellData[cellIndex][variableIndex(2,1)]  =PerpStrain[0];
          cellData[cellIndex][variableIndex(2,1)+1]=PerpStrain[1];
          cellData[cellIndex][variableIndex(2,1)+2]=PerpStrain[2];
          cellData[cellIndex][variableIndex(2,1)+3]=growthStrain;  //growth Strain Value is stored after its eigenvector
          
        }
    }
    
    
    
    // normalizing MT vector
    temp=std::sqrt(cellData[cellIndex][MTindex  ]*cellData[cellIndex][MTindex  ] +
                   cellData[cellIndex][MTindex+1]*cellData[cellIndex][MTindex+1] +
                   cellData[cellIndex][MTindex+2]*cellData[cellIndex][MTindex+2]   );
    if(temp>0){      
      cellData[cellIndex][MTindex  ]=cellData[cellIndex][MTindex  ]/temp;
      cellData[cellIndex][MTindex+1]=cellData[cellIndex][MTindex+1]/temp;
      cellData[cellIndex][MTindex+2]=cellData[cellIndex][MTindex+2]/temp;
    }

 
    //<<<<<<<<<<<<<<<<<<<<<<<< angles between  vectors and circumferential direction <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   
    //cellData[cellIndex][15]= cellData[cellIndex][comIndex+2]; // z coordinate of central vertex of the cell // 15 --> Z coordinate
        
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  

    // cellData[cellIndex][20]= strainZ;
    // cellData[cellIndex][21]= TetaStress;
    // cellData[cellIndex][22]= TetaStrain; 
    // cellData[cellIndex][23]= TetaPerp;
    // cellData[cellIndex][19]= TETA;
    // cellData[cellIndex][20]=EnergyIso+EnergyAniso;    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
    //cellData[cellIndex][areaRatioIndex  ]= areaRatio;
    cellData[cellIndex][areaRatioIndex  ]= youngL/youngT;
    cellData[cellIndex][isoEnergyIndex  ]= EnergyIso;    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    cellData[cellIndex][anisoEnergyIndex]= EnergyAniso;
    
  }
  //std::cerr<<"up to here"<<std::endl; 
  
  double totalEnergyIso=0;
  double totalEnergyAniso=0;
  
  for( size_t cellIndex=0 ; cellIndex<numCells ; ++cellIndex ) {
    
    const size_t numWalls = T.cell(cellIndex).numWall();
    //if(cellData[cellIndex][39+2]>-120) {
    totalEnergyIso   +=cellData[cellIndex][isoEnergyIndex];
    totalEnergyAniso +=cellData[cellIndex][anisoEnergyIndex];
    // }

    
    Cell *  cell1=&(T.cell(cellIndex));
    	
    //std::vector<int> neighbor(numWalls);
    std::vector<size_t> neighbor(numWalls);
    for   ( size_t wallIndex=0 ; wallIndex<numWalls ; ++wallIndex ){
      neighbor[wallIndex]=(cell1->cellNeighbor(wallIndex))->index();
    }
    
    // std::cerr<<" cell   "<<cellIndex << "cell neighbors    ";
    // for   ( size_t wallIndex=0 ; wallIndex<numWalls ; ++wallIndex ){
    //   std::cerr<< "  "<<neighbor[wallIndex]<<" " ;
    // }
    

    double neighborweight=parameter(5);
    double normalGlob[3]={0,0,0};
    
    normalGlob[0]= cellData[cellIndex][normalVectorIndex  ];
    normalGlob[1]= cellData[cellIndex][normalVectorIndex+1];
    normalGlob[2]= cellData[cellIndex][normalVectorIndex+2]; 
     
  
    
    double StressTensor[3][3]={{0,0,0},{0,0,0},{0,0,0}};
    
    int counter=0;
    for (size_t nn=0 ; nn<numWalls ; nn++){
      if (neighbor[nn]<numCells && neighbor[nn]>-1){

  	StressTensor[0][0]+=neighborweight*cellData[neighbor[nn]][stressTensorIndex  ];
  	StressTensor[1][1]+=neighborweight*cellData[neighbor[nn]][stressTensorIndex+1];
  	StressTensor[2][2]+=neighborweight*cellData[neighbor[nn]][stressTensorIndex+2];
  	StressTensor[0][1]+=neighborweight*cellData[neighbor[nn]][stressTensorIndex+3];
  	StressTensor[2][0]+=neighborweight*cellData[neighbor[nn]][stressTensorIndex+4];
  	StressTensor[1][2]+=neighborweight*cellData[neighbor[nn]][stressTensorIndex+5];
	counter+=1;
      }
    }
    if(counter !=0){
      StressTensor[0][0]/=counter;
      StressTensor[1][1]/=counter;
      StressTensor[2][2]/=counter;
      StressTensor[0][1]/=counter;
      StressTensor[2][0]/=counter;
      StressTensor[1][2]/=counter;
    }
   
    StressTensor[0][0]+=(1-neighborweight)*cellData[cellIndex][stressTensorIndex  ];
    StressTensor[1][1]+=(1-neighborweight)*cellData[cellIndex][stressTensorIndex+1];
    StressTensor[2][2]+=(1-neighborweight)*cellData[cellIndex][stressTensorIndex+2];
    StressTensor[0][1]+=(1-neighborweight)*cellData[cellIndex][stressTensorIndex+3];
    StressTensor[2][0]+=(1-neighborweight)*cellData[cellIndex][stressTensorIndex+4];
    StressTensor[1][2]+=(1-neighborweight)*cellData[cellIndex][stressTensorIndex+5];
   
 
    
    StressTensor[0][2]=StressTensor[2][0];
    StressTensor[1][0]=StressTensor[0][1];
    StressTensor[2][1]=StressTensor[1][2];
    
    
    // stress component along MT direction 
    cellData[cellIndex][MTstressIndex ] =
      cellData[cellIndex][MTindex  ] *cellData[cellIndex][MTindex  ] *StressTensor[0][0]  +
      cellData[cellIndex][MTindex  ] *cellData[cellIndex][MTindex+1] *StressTensor[0][1]  +
      cellData[cellIndex][MTindex  ] *cellData[cellIndex][MTindex+2] *StressTensor[0][2]  +
      cellData[cellIndex][MTindex+1] *cellData[cellIndex][MTindex  ] *StressTensor[1][0]  +
      cellData[cellIndex][MTindex+1] *cellData[cellIndex][MTindex+1] *StressTensor[1][1]  +
      cellData[cellIndex][MTindex+1] *cellData[cellIndex][MTindex+2] *StressTensor[1][2]  +
      cellData[cellIndex][MTindex+2] *cellData[cellIndex][MTindex  ] *StressTensor[2][0]  +
      cellData[cellIndex][MTindex+2] *cellData[cellIndex][MTindex+1] *StressTensor[2][1]  +
      cellData[cellIndex][MTindex+2] *cellData[cellIndex][MTindex+2] *StressTensor[2][2]       ;
    
  
    
    // eigenvalue/eigenvectors of  stress tensor in global coordinate system. (Jacobi method)
    
    double pi=3.14159265;
    double stressEpcilon=0.000001;    
    int I,J;    
    double pivot=1;
    double RotAngle,Si,Co;
    double eigenVectorStress[3][3]={{1,0,0},{0,1,0},{0,0,1}};
    pivot=1;
    
    while (pivot>stressEpcilon) {
      pivot=std::abs(StressTensor[1][0]);
      I=1;
      J=0;
      if (std::abs(StressTensor[2][0])>pivot) {
  	pivot=std::abs(StressTensor[2][0]);
  	I=2;
  	J=0;
      }
      if (std::abs(StressTensor[2][1])>pivot) {
  	pivot=std::abs(StressTensor[2][1]);
  	I=2;
  	J=1;
      }
      if (std::abs(StressTensor[I][I]-StressTensor[J][J])<stressEpcilon) {
  	RotAngle=pi/4;
      }            
      else {
  	RotAngle=0.5*std::atan((2*StressTensor[I][J])/(StressTensor[J][J]-StressTensor[I][I]));
      }
      Si=std::sin(RotAngle);
      Co=std::cos(RotAngle);
      double tempRot[3][3]={{1,0,0},{0,1,0},{0,0,1}};
      tempRot[I][I]=Co;
      tempRot[J][J]=Co;
      tempRot[I][J]=Si;
      tempRot[J][I]=-Si;
      
      double tempStress[3][3]={{0,0,0},{0,0,0},{0,0,0}};
      for (int r=0 ; r<3 ; r++) 
  	for (int s=0 ; s<3 ; s++) 
  	  for(int w=0 ; w<3 ; w++) 
  	    tempStress[r][s]=tempStress[r][s]+StressTensor[r][w]*tempRot[w][s];
      
      for (int r=0 ; r<3 ; r++) 
  	for (int s=0 ; s<3 ; s++) 
  	  StressTensor[r][s]=0;
      
      for (int r=0 ; r<3 ; r++) 
  	for (int s=0 ; s<3 ; s++) 
  	  for(int w=0 ; w<3 ; w++) 
  	    StressTensor[r][s]=StressTensor[r][s]+tempRot[w][r]*tempStress[w][s];
      
      for (int r=0 ; r<3 ; r++) 
  	for (int s=0 ; s<3 ; s++) 
  	  tempStress[r][s]=eigenVectorStress[r][s];
      
      for (size_t ii=0; ii<3; ++ii) {
  	for (size_t jj=0; jj<3; ++jj) {
  	  eigenVectorStress[ii][jj] = 0.0;
  	}
      }
      
      for (int r=0 ; r<3 ; r++) 
  	for (int s=0 ; s<3 ; s++) 
  	  for(int w=0 ; w<3 ; w++) 
  	    eigenVectorStress[r][s]=eigenVectorStress[r][s]+tempStress[r][w]*tempRot[w][s];
    }
  
    // normalizing eigenvectors (remove if not necessary)  
    double temp=std::sqrt(eigenVectorStress[0][0]*eigenVectorStress[0][0] +
  		   eigenVectorStress[1][0]*eigenVectorStress[1][0] +
  		   eigenVectorStress[2][0]*eigenVectorStress[2][0] );
    if(temp>0){
      eigenVectorStress[0][0]/=temp;
      eigenVectorStress[1][0]/=temp;
      eigenVectorStress[2][0]/=temp;
    }
    temp=std::sqrt(eigenVectorStress[0][1]*eigenVectorStress[0][1] +
  		   eigenVectorStress[1][1]*eigenVectorStress[1][1] +
  		   eigenVectorStress[2][1]*eigenVectorStress[2][1] );
    if(temp>0){
      eigenVectorStress[0][1]/=temp;
      eigenVectorStress[1][1]/=temp;
      eigenVectorStress[2][1]/=temp;
    }
    temp=std::sqrt(eigenVectorStress[0][2]*eigenVectorStress[0][2] +
  		   eigenVectorStress[1][2]*eigenVectorStress[1][2] +
  		   eigenVectorStress[2][2]*eigenVectorStress[2][2] );
    if(temp>0){
      eigenVectorStress[0][2]/=temp;
      eigenVectorStress[1][2]/=temp;
      eigenVectorStress[2][2]/=temp;
    } 
    
    // maximal stress direction
    double maximalStressValue=StressTensor[0][0];
    int Istress=0;
    if ( StressTensor[1][1] > maximalStressValue ) 
      {
  	maximalStressValue=StressTensor[1][1];
  	Istress=1;
      }
    if ( StressTensor[2][2] > maximalStressValue ) 
      {
  	maximalStressValue=StressTensor[2][2];
  	Istress=2;
      }
    
    
    // 2nd maximalstress direction/value
    double maximalStressValue2, maximalStressValue3;
    int Istress2,Istress3;
    if (Istress==0) {
      Istress2=1;
      Istress3=2;
    }
    if (Istress==1) {
      Istress2=2;
      Istress3=0;
    }
    if (Istress==2) {
      Istress2=0;
      Istress3=1;
    }
    if(StressTensor[Istress3][Istress3] > StressTensor[Istress2][Istress2] ) {
      temp=Istress2;
      Istress2=Istress3;
      Istress3=temp;
    }
    maximalStressValue2=StressTensor[Istress2][Istress2];
    maximalStressValue3=StressTensor[Istress3][Istress3];
    
    
    if(std::abs(normalGlob[0]*eigenVectorStress[0][Istress]+
  		normalGlob[1]*eigenVectorStress[1][Istress]+
  		normalGlob[2]*eigenVectorStress[2][Istress]) > .7) {
      // std::cerr << "max stress normal to the cell plane "<<cellIndex<<std::endl;
      Istress=Istress2; 
      Istress2=Istress3;
      maximalStressValue=maximalStressValue2;
      maximalStressValue2=maximalStressValue3; 
    }
    
    


    // storing a measure for stress anisotropy in cell vector
    if (std::abs(maximalStressValue)<  0.000001) cellData[cellIndex][stressAnIndex]=0;

    if (std::abs(maximalStressValue)>= 0.000001) {
      if(parameter(6)==0)
	cellData[cellIndex][stressAnIndex]=1-std::abs(maximalStressValue2/maximalStressValue);
      else
	cellData[cellIndex][stressAnIndex]=(1-std::abs(maximalStressValue2/maximalStressValue))*(maximalStressValue/parameter(6));
    }    

 
    
    
    
    if (numVariableIndexLevel()==4 &&(numVariableIndex(3)==1 || numVariableIndex(3)==2 ) ) { // storing maximal stress
      if (dimension==2)
  	{
  	  cellData[cellIndex][variableIndex(3,0)]  =eigenVectorStress[0][Istress];
  	  cellData[cellIndex][variableIndex(3,0)+1]=eigenVectorStress[1][Istress];
  	  cellData[cellIndex][variableIndex(3,0)+3]=maximalStressValue;  //maximal Stress Value is stored after its eigenvector
  	}
      if (dimension==3)
  	{ 
  	  cellData[cellIndex][variableIndex(3,0)]  =eigenVectorStress[0][Istress];
  	  cellData[cellIndex][variableIndex(3,0)+1]=eigenVectorStress[1][Istress];
  	  cellData[cellIndex][variableIndex(3,0)+2]=eigenVectorStress[2][Istress];
	  
  	  cellData[cellIndex][variableIndex(3,0)+3]=maximalStressValue;  //maximal Stress Value is stored after its eigenvector
  	}
    }
    
    
    if (numVariableIndexLevel()==4 && numVariableIndex(3)==2 ) { // storing 2nd maximal stress
      if (dimension==2)
  	{
  	  cellData[cellIndex][variableIndex(3,1)]  =eigenVectorStress[0][Istress2];
  	  cellData[cellIndex][variableIndex(3,1)+1]=eigenVectorStress[1][Istress2];
  	  cellData[cellIndex][variableIndex(3,1)+3]=maximalStressValue2;  //2nd maximal Stress Value is stored after its eigenvector
  	}
      if (dimension==3)
  	{
  	  cellData[cellIndex][variableIndex(3,1)]  =eigenVectorStress[0][Istress2];
  	  cellData[cellIndex][variableIndex(3,1)+1]=eigenVectorStress[1][Istress2];
  	  cellData[cellIndex][variableIndex(3,1)+2]=eigenVectorStress[2][Istress2];
  	  cellData[cellIndex][variableIndex(3,1)+3]=maximalStressValue2;  //2nd maximal Stress Value is stored after its eigenvector
  	}
    }
    
    
    
    
  }


 
     
  
  //cellData[0][14]=totalEnergyIso ;
  //cellData[1][14]=totalEnergyAniso ;   
    
  //std::cerr<<"here.................................."<<std::endl;
  //std::cerr<<"here.................................."<<std::endl; 
    
}     



//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

// void VertexFromTRBScenterTriangulationMT::
// initiate(Tissue &T,
// 	 DataMatrix &cellData,
// 	 DataMatrix &wallData,
// 	 DataMatrix &vertexData,
// 	 DataMatrix &cellDerivs,
// 	 DataMatrix &wallDerivs,
// 	 DataMatrix &vertexDerivs)
// {
//   size_t dimension=3; //Only implemented for 3D models
//   assert (dimension==vertexData[0].size());
//   size_t numVariable = T.cell(0).numVariable();
//   assert (numVariable==cellData[0].size());
//   // Create the new variables
//   if (variableIndex(1,0) != numVariable) {
//     std::cerr << "VertexFromTRBScenterTriangulation::initiate() "
// 	      << "Wrong index given as start index for additional variables."
// 	      << std::endl;
//     exit(-1);
//   }
//   size_t numCell = cellData.size();
//   assert (numCell==T.numCell());
//   std::vector<double> com(dimension);
  
//   for (size_t cellIndex=0; cellIndex<numCell; ++cellIndex) {
//     size_t numInternalWall = T.cell(cellIndex).numVertex();
//     cellData[cellIndex].resize(numVariable+dimension+numInternalWall);
//     cellDerivs[cellIndex].resize(numVariable+dimension+numInternalWall);
//     com = T.cell(cellIndex).positionFromVertex(vertexData);
//     // Set center position to com of the cell
//     for (size_t d=0; d<dimension; ++d)
//       cellData[cellIndex][numVariable+d] = com[d];    
//     // Set internal wall lengths to the distance btw com and the vertex
//     for (size_t wallindex=0; wallindex<numInternalWall; ++wallindex) {
//       Vertex *tmpVertex = T.cell(cellIndex).vertex(wallindex); 
//       size_t vertexIndex = tmpVertex->index();
//       double distance = std::sqrt( (com[0]-vertexData[vertexIndex][0])*
// 				   (com[0]-vertexData[vertexIndex][0])+
// 				   (com[1]-vertexData[vertexIndex][1])*
// 				   (com[1]-vertexData[vertexIndex][1])+
// 				   (com[2]-vertexData[vertexIndex][2])*
// 				   (com[2]-vertexData[vertexIndex][2]) );   
//       cellData[cellIndex][numVariable+dimension+wallindex] = distance;
//     }
//   }
// }



VertexFromTRBScenterTriangulationConcentrationHillMT::
VertexFromTRBScenterTriangulationConcentrationHillMT(std::vector<double> &paraValue, 
                                                     std::vector< std::vector<size_t> > 
                                                     &indValue ) 
{  
  // Do some checks on the parameters and variable indeces
  if( paraValue.size()!=8 ) {
    std::cerr << "VertexFromTRBScenterTriangulationConcentrationHillMT::"
	      << "VertexFromTRBScenterTriangulationConcentrationHillMT() "
	      << "Uses eight parameters Young_modulus_Long_min, Young_modulus_Long_max, "                           
	      << " poisson coefficient_long,Young_modulus_Trans_min, Young_modulus_Trans_max, "
	      << " poisson coefficient_Trans, K_hill, and n_hill.\n";
    exit(0);
  }
  
  if( (indValue.size()!=2 && indValue.size()!=6) || 
      indValue[0].size()!=3 || indValue[1].size()!=1 ||
      (indValue.size()==6 && (indValue[2].size()!=0 && indValue[2].size()!=1)) ||
      (indValue.size()==6 && (indValue[3].size()!=0 && indValue[3].size()!=1)) ||
      (indValue.size()==6 && (indValue[4].size()!=0 && indValue[4].size()!=1)) ||
      (indValue.size()==6 && (indValue[5].size()!=0 && indValue[5].size()!=1)) 
      ) { 
    std::cerr << "VertexFromTRBScenterTriangulationConcentrationHillMT::"
	      << "VertexFromTRBScenterTriangulationConcentrationHillMT() "
	      << "Wall length and concentration indices and MT direction initial index given in first level." << std::endl
	      << "Start of additional Cell variable indices (center(x,y,z) "
	      << "L_1,...,L_n, n=num vertex) is given in second level (typically at end)." 
              << "Optionally four additional levels can be given where the strain, 2nd strain, stress and 2nd stress "
	      << "directions/values(dx dy dz value) can be stored at given indices. If index given at third level, "
              << "strain direction will be stored starting at this (cell) variable index, "
	      << "for fourth 2nd strain , for fifth level stress and for sixth level 2nd stress will be stored."
	      << std::endl;
    exit(0);
  }
  
  // Set the variable values
  setId("VertexFromTRBScenterTriangulationConcentrationHillMT");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  // Set the parameter identities
  std::vector<std::string> tmp( numParameter() );   
  tmp[0] = "Y_mod_L_min";
  tmp[1] = "Y_mod_L_max";
  tmp[2] = "P_ratio_L";
  tmp[3] = "Y_mod_T_min";
  tmp[4] = "Y_mod_T_max";
  tmp[5] = "P_ratio_T";
  tmp[6] = "K_hill";
  tmp[7] = "n_hill";
  setParameterId( tmp );
}

void VertexFromTRBScenterTriangulationConcentrationHillMT::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) {
  
  //Do the update for each cell
  size_t dimension = 3;
  assert (dimension==vertexData[0].size());
  size_t numCells = T.numCell();
  size_t wallLengthIndex = variableIndex(0,0);                       
  size_t concIndex = variableIndex(0,1);
  //variableIndex(0,2) index for MT direction
  size_t comIndex = variableIndex(1,0);
  size_t lengthInternalIndex = comIndex+dimension;
  double Kpow = std::pow(parameter(6),parameter(7));                    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  for (size_t cellIndex=0 ; cellIndex<numCells ; ++cellIndex) {
    size_t numWalls = T.cell(cellIndex).numWall(); 
    
    if(  T.cell(cellIndex).numVertex()!= numWalls ) {
      std::cerr << "VertexFromTRBScenterTriangulationConcentrationHillMT::derivs() "
		<< "same number of vertices and walls."
		<< " Not for cells with " << T.cell(cellIndex).numWall() << " walls and "
		<< T.cell(cellIndex).numVertex() << " vertices!"	
		<< std::endl;
      exit(-1);
    }
    
    double youngL = parameter(0) + 
      parameter(1)*Kpow/( Kpow+std::pow(cellData[cellIndex][concIndex],parameter(7)) );
    double poissonL =parameter(2);
    double youngT = parameter(3) + 
      parameter(4)*Kpow/( Kpow+std::pow(cellData[cellIndex][concIndex],parameter(7)) );
    double poissonT =parameter(5);
    
    double StrainCellGlobal[3][3]={{0,0,0},{0,0,0},{0,0,0}};
    double StressCellGlobal[3][3]={{0,0,0},{0,0,0},{0,0,0}};
    double TotalCellRestingArea=0;
    
    
    // One triangle per 'vertex' in cyclic order
    for (size_t k=0; k<numWalls; ++k) { 
      size_t kPlusOneMod = (k+1)%numWalls;
      //size_t v1 = com;
      size_t v2 = T.cell(cellIndex).vertex(k)->index();
      size_t v3 = T.cell(cellIndex).vertex(kPlusOneMod)->index();
      //size_t w1 = internal k
      size_t w2 = T.cell(cellIndex).wall(k)->index();
      //size_t w3 = internal k+1
      
      // Position matrix holds in rows positions for com, vertex(k), vertex(k+1)
      DataMatrix position(3,vertexData[v2]);
      for (size_t d=0; d<dimension; ++d)
	position[0][d] = cellData[cellIndex][comIndex+d]; // com position
      //position[1] = vertexData[v2]; // given by initiation
      position[2] = vertexData[v3];
      
      // Resting lengths are from com-vertex(k), vertex(k)-vertex(k+1) (wall(k)), com-vertex(k+1)
      std::vector<double> restingLength(3);
      restingLength[0] = cellData[cellIndex][lengthInternalIndex + k];
      restingLength[1] = wallData[w2][wallLengthIndex];
      restingLength[2] = cellData[cellIndex][lengthInternalIndex + kPlusOneMod];
      
      // Lengths are from com-vertex(k), vertex(k)-vertex(k+1) (wall(k)), com-vertex(k+1)
      std::vector<double> length(3);
      length[0] = std::sqrt( (position[0][0]-position[1][0])*(position[0][0]-position[1][0]) +
			     (position[0][1]-position[1][1])*(position[0][1]-position[1][1]) +
			     (position[0][2]-position[1][2])*(position[0][2]-position[1][2]) );
      
      length[1] = T.wall(w2).lengthFromVertexPosition(vertexData);

      length[2] = std::sqrt( (position[0][0]-position[2][0])*(position[0][0]-position[2][0]) +
			     (position[0][1]-position[2][1])*(position[0][1]-position[2][1]) +
			     (position[0][2]-position[2][2])*(position[0][2]-position[2][2]) );
      
      // Lame coefficients (can be defined out of loop)
      double lambdaL=youngL*poissonL/(1-poissonL*poissonL);
      double mioL=youngL/(1+poissonL);
      double lambdaT=youngT*poissonT/(1-poissonT*poissonT);
      double mioT=youngT/(1+poissonT);
      
      // Area of the element (using Heron's formula)                                      
      double restingArea=std::sqrt( ( restingLength[0]+restingLength[1]+restingLength[2])*
			     (-restingLength[0]+restingLength[1]+restingLength[2])*
			     ( restingLength[0]-restingLength[1]+restingLength[2])*
			     ( restingLength[0]+restingLength[1]-restingLength[2])  )*0.25;
     
      //double currentArea=std::sqrt( ( length[0]+length[1]+length[2])*
      //                              (-length[0]+length[1]+length[2])*
      //                              ( length[0]-length[1]+length[2])*
      //                              ( length[0]+length[1]-length[2])  )*0.25;
      
      
      //Angles of the element ( assuming the order: 0,L0,1,L1,2,L2 )
      std::vector<double> Angle(3);
      // can be ommited by cotan(A)=.25*sqrt(4*b*b*c*c/K-1)
      Angle[0]=std::acos(  (restingLength[0]*restingLength[0]+restingLength[2]*restingLength[2]-restingLength[1]*restingLength[1])/
                           (restingLength[0]*restingLength[2]*2)    );
      Angle[1]=std::acos(  (restingLength[0]*restingLength[0]+restingLength[1]*restingLength[1]-restingLength[2]*restingLength[2])/
                           (restingLength[0]*restingLength[1]*2)    );
      Angle[2]=std::acos(  (restingLength[1]*restingLength[1]+restingLength[2]*restingLength[2]-restingLength[0]*restingLength[0])/
                           (restingLength[1]*restingLength[2]*2)    );
      
      //Tensile Stiffness
      double tensileStiffness[3];
      double temp = 1.0/(restingArea*16);                                      
      double cotan[3] = {1.0/std::tan(Angle[0]),1.0/std::tan(Angle[1]),1.0/std::tan(Angle[2])};    
      tensileStiffness[0]=(2*cotan[2]*cotan[2]*(lambdaT+mioT)+mioT)*temp;
      tensileStiffness[1]=(2*cotan[0]*cotan[0]*(lambdaT+mioT)+mioT)*temp;
      tensileStiffness[2]=(2*cotan[1]*cotan[1]*(lambdaT+mioT)+mioT)*temp;
      
      //Angular Stiffness
      double angularStiffness[3];
      angularStiffness[0]=(2*cotan[1]*cotan[2]*(lambdaT+mioT)-mioT)*temp;
      angularStiffness[1]=(2*cotan[0]*cotan[2]*(lambdaT+mioT)-mioT)*temp;
      angularStiffness[2]=(2*cotan[0]*cotan[1]*(lambdaT+mioT)-mioT)*temp;
      
      //Calculate biquadratic strains  
      std::vector<double> Delta(3);
      Delta[0]=(length[0])*(length[0])-(restingLength[0])*(restingLength[0]);
      Delta[1]=(length[1])*(length[1])-(restingLength[1])*(restingLength[1]);
      Delta[2]=(length[2])*(length[2])-(restingLength[2])*(restingLength[2]);

    //Area of the element (using Heron's formula)                                      
    double Area=std::sqrt( ( length[0]+length[1]+length[2])*
                           (-length[0]+length[1]+length[2])*
                           ( length[0]-length[1]+length[2])*
                           ( length[0]+length[1]-length[2])  )*0.25;
    
    // calculating the angles between shape vectors and anisotropy direction in resting shape when anisotropy vector is provided in current shape
    
    //Current shape local coordinate of the element  (counterclockwise ordering of nodes/edges)
      double CurrentAngle1=std::acos(  (length[0]*length[0]+length[1]*length[1]-length[2]*length[2])/
                                       (length[0]*length[1]*2)    );

      double Qa=std::cos(CurrentAngle1)*length[0];
      double Qc=std::sin(CurrentAngle1)*length[0];
      double Qb=length[1];
      // shape vector matrix = inverse of coordinate matrix ( only first two elements i.e. ShapeVector[3][2] )      
      double ShapeVectorCurrent[3][3]={ {  0   ,       1/Qc      , 0 }, 
                                        {-1/Qb , (Qa-Qb)/(Qb*Qc) , 1 },       
                                        { 1/Qb ,     -Qa/(Qb*Qc) , 0 }  };
            
      // std::cerr<< "cell "<< cellIndex<< " shape vactor 0 : "<< ShapeVectorCurrent[0][0]<<"  "<< ShapeVectorCurrent[0][1]<<"  "<< ShapeVectorCurrent[0][2]<< std::endl;
      // std::cerr<< "cell "<< cellIndex<< " shape vactor 1 : "<< ShapeVectorCurrent[1][0]<<"  "<< ShapeVectorCurrent[1][1]<<"  "<< ShapeVectorCurrent[1][2]<< std::endl;
      // std::cerr<< "cell "<< cellIndex<< " shape vactor 2 : "<< ShapeVectorCurrent[2][0]<<"  "<< ShapeVectorCurrent[2][1]<<"  "<< ShapeVectorCurrent[2][2]<< std::endl;
      
      // Local coordinates of the resting shape ( counterclockwise )
      double RestingAngle1=std::acos(  (restingLength[0]*restingLength[0]+restingLength[1]*restingLength[1]-restingLength[2]*restingLength[2])/
                                       (restingLength[0]*restingLength[1]*2)    );

      double Pa=std::cos(RestingAngle1)*restingLength[0];
      double Pc=std::sin(RestingAngle1)*restingLength[0];
      double Pb=restingLength[1];

      // shape vector matrix in resting shape in local coordinate system  = inverse of coordinate matrix ( only first two elements i.e. ShapeVectorResting[3][2] )      
      double ShapeVectorResting[3][3]={ {  0   ,       1/Pc      , 0 }, 
                                        {-1/Pb , (Pa-Pb)/(Pb*Pc) , 1 },       
                                        { 1/Pb ,     -Pa/(Pb*Pc) , 0 }  };
      // std::cerr<< "cell "<< cellIndex<< " shape vactor 0 : "<< ShapeVectorResting[0][0]<<"  "<< ShapeVectorResting[0][1]<<"  "<< ShapeVectorResting[0][2]<< std::endl;
      // std::cerr<< "cell "<< cellIndex<< " shape vactor 1 : "<< ShapeVectorResting[1][0]<<"  "<< ShapeVectorResting[1][1]<<"  "<< ShapeVectorResting[1][2]<< std::endl;
      // std::cerr<< "cell "<< cellIndex<< " shape vactor 2 : "<< ShapeVectorResting[2][0]<<"  "<< ShapeVectorResting[2][1]<<"  "<< ShapeVectorResting[2][2]<< std::endl;
      
      // //Strain tensor  (clockwise nodes/edges)
      // double CurrentAngle2=std::acos(  (length[1]*length[1]+length[2]*length[2]-length[0]*length[0])/
      //                                  (length[1]*length[2]*2)    );

      // double Qa=std::cos(CurrentAngle2)*length[2];
      // double Qb=length[1];
      // double Qc=std::sin(CurrentAngle2)*length[2];
      
      // double ShapeVectorCurrent[3][2]={ {  0   ,       1/Qc      }, 
      //                                   { 1/Qb ,     -Qa/(Qb*Qc) },       
      //                                   {-1/Qb , (Qa-Qb)/(Qb*Qc) }  };

      // Aniso vector in current shape in global coordinate system
      double AnisoCurrGlob[3];
      AnisoCurrGlob[0] = cellData[cellIndex][variableIndex(0,2)];  // this was written by Henrik, Behruz changed it
      AnisoCurrGlob[1] = cellData[cellIndex][variableIndex(0,2)+1];
      AnisoCurrGlob[2] = cellData[cellIndex][variableIndex(0,2)+2];
     
     
     
      // Rotation Matrix for changing coordinate systems for both Local to Global( Strain Tensor) and Global to Local( Aniso Vector in the current shape)
      double rotation[3][3];  

      double tempA=std::sqrt((position[2][0]-position[1][0])*(position[2][0]-position[1][0])+
                             (position[2][1]-position[1][1])*(position[2][1]-position[1][1])+
                             (position[2][2]-position[1][2])*(position[2][2]-position[1][2])  );

      double tempB=std::sqrt((position[0][0]-position[1][0])*(position[0][0]-position[1][0])+
                             (position[0][1]-position[1][1])*(position[0][1]-position[1][1])+
                             (position[0][2]-position[1][2])*(position[0][2]-position[1][2])  );

      double Xcurrent[3];      
      Xcurrent[0]= (position[2][0]-position[1][0])/tempA;
      Xcurrent[1]= (position[2][1]-position[1][1])/tempA;
      Xcurrent[2]= (position[2][2]-position[1][2])/tempA;
      
      double Bcurrent[3];      
      Bcurrent[0]= (position[0][0]-position[1][0])/tempB;
      Bcurrent[1]= (position[0][1]-position[1][1])/tempB;
      Bcurrent[2]= (position[0][2]-position[1][2])/tempB;

      double Zcurrent[3];      
      Zcurrent[0]= Xcurrent[1]*Bcurrent[2]-Xcurrent[2]*Bcurrent[1];
      Zcurrent[1]= Xcurrent[2]*Bcurrent[0]-Xcurrent[0]*Bcurrent[2];
      Zcurrent[2]= Xcurrent[0]*Bcurrent[1]-Xcurrent[1]*Bcurrent[0];
      
      tempA=std:: sqrt(Zcurrent[0]*Zcurrent[0]+Zcurrent[1]*Zcurrent[1]+Zcurrent[2]*Zcurrent[2]);
      Zcurrent[0]=Zcurrent[0]/tempA;
      Zcurrent[1]=Zcurrent[1]/tempA;
      Zcurrent[2]=Zcurrent[2]/tempA;

      double Ycurrent[3];      
      Ycurrent[0]= Zcurrent[1]*Xcurrent[2]-Zcurrent[2]*Xcurrent[1];
      Ycurrent[1]= Zcurrent[2]*Xcurrent[0]-Zcurrent[0]*Xcurrent[2];
      Ycurrent[2]= Zcurrent[0]*Xcurrent[1]-Zcurrent[1]*Xcurrent[0];


      rotation[0][0]=Xcurrent[0];
      rotation[1][0]=Xcurrent[1];
      rotation[2][0]=Xcurrent[2];

      rotation[0][1]=Ycurrent[0];
      rotation[1][1]=Ycurrent[1];
      rotation[2][1]=Ycurrent[2];

      rotation[0][2]=Zcurrent[0];
      rotation[1][2]=Zcurrent[1];
      rotation[2][2]=Zcurrent[2];      
       
      // rotating the anisotropy vector from global coordinate system to the local one in the current shape
      double AnisoCurrLocal[3];
      AnisoCurrLocal[0]=rotation[0][0]*AnisoCurrGlob[0]+rotation[1][0]*AnisoCurrGlob[1]+rotation[2][0]*AnisoCurrGlob[2];
      AnisoCurrLocal[1]=rotation[0][1]*AnisoCurrGlob[0]+rotation[1][1]*AnisoCurrGlob[1]+rotation[2][1]*AnisoCurrGlob[2];
      AnisoCurrLocal[2]=rotation[0][2]*AnisoCurrGlob[0]+rotation[1][2]*AnisoCurrGlob[1]+rotation[2][2]*AnisoCurrGlob[2];
     
      //std::cerr<< "cell "<< cellIndex<< " anisoVector current local "<<AnisoCurrLocal[0]<<"  "<<AnisoCurrLocal[1]<<"  "<<AnisoCurrLocal[2]<< std::endl;
      // Center of Mass for current shape in local coordinate
      double CMCurrentLocal[2]={(Qa+Qb)/3, Qc/3};
      
      // Tip of the ansotropy vector drown from the center of mass in the current shape
      double  ACurrLocal[2] = {CMCurrentLocal[0]+AnisoCurrLocal[0],CMCurrentLocal[1]+AnisoCurrLocal[1]};

      //std::cerr<< "cell "<< cellIndex<< " tip of anisoVector in local current "<<ACurrLocal[0]<<"  "<<ACurrLocal[1]<< std::endl;

      // Baricentric Coordinates of tip of anisotropy vector in the current shape wihich is equivalent to the baricentric coordinate of the corresponding point in the resting shape
      double ABari[3];
      ABari[0]=ShapeVectorCurrent[0][0]*ACurrLocal[0]+ShapeVectorCurrent[0][1]*ACurrLocal[1]+ShapeVectorCurrent[0][2];
      ABari[1]=ShapeVectorCurrent[1][0]*ACurrLocal[0]+ShapeVectorCurrent[1][1]*ACurrLocal[1]+ShapeVectorCurrent[1][2];
      ABari[2]=ShapeVectorCurrent[2][0]*ACurrLocal[0]+ShapeVectorCurrent[2][1]*ACurrLocal[1]+ShapeVectorCurrent[2][2];
      //std::cerr<< "cell "<< cellIndex<< " baricentric coor of tip of anisoVector in local current "<<ABari[0]<<"  "<<ABari[1]<<"  "<<ABari[2]<< std::endl;

      
      // Local coordinates of tip of AnisoVector in the resting shape from multyplying ABari by position matrix in the local resting shape coordinates      
      double ARestLocal[2];
      ARestLocal[0]=Pa*ABari[0]+Pb*ABari[2];
      ARestLocal[1]=Pc*ABari[0];
      //std::cerr<< "cell "<< cellIndex<< " tip of anisoVector in local rest "<<ARestLocal[0]<<"  "<<ARestLocal[1]<< std::endl;     


      // Center of Mass for resting shape in local coordinate
      double CMRestingLocal[2]={(Pa+Pb)/3, Pc/3};
            
      // Aniso Vector in the resting shape in local coordinate system
      double AnisoRestLocal[2]={ARestLocal[0]-CMRestingLocal[0],ARestLocal[1]-CMRestingLocal[1]};
      //std::cerr<< "cell "<< cellIndex<< " anisoVector Rest "<<AnisoRestLocal[0]<<"  "<<AnisoRestLocal[1]<<"  "<<AnisoRestLocal[2]<< std::endl;    

      // Anisotropy measure or magnitude of the projection of anisotropy vector on the cell plane in the resting shape 
      //this measure is considered as a factor controling the anisotropy properties of the cell
      double AnisoMeasure=std::sqrt(AnisoRestLocal[0]*AnisoRestLocal[0]+AnisoRestLocal[1]*AnisoRestLocal[1]);
      // std::cerr<< "cell "<< cellIndex<<" AnisoMeasure "<< AnisoMeasure<<std::endl;
      
      // choosing a random normalized dirrection for anisotropy if AnisoVector is close to perpendicular to the cell plane
      if ( AnisoMeasure<0.001) {
        double randomAngle=((rand()%360)*2*3.14159265)/360;
        
        AnisoRestLocal[0]=std::cos(randomAngle);
        AnisoRestLocal[1]=std::sin(randomAngle);
      }  
      else {// normalizing the anisoVector if it is not random
        AnisoRestLocal[0]=AnisoRestLocal[0]/AnisoMeasure;
        AnisoRestLocal[1]=AnisoRestLocal[1]/AnisoMeasure;        
      }
      
      //Anisotropic Correction is based on difference between Lam Coefficients of Longitudinal and Transverse dirrections:
      double deltaLam=AnisoMeasure*(lambdaL-lambdaT);
      double deltaMio=AnisoMeasure*(mioL-mioT);
      
      //Angles between anisotropy vector and shape vectors for calculating the terms like a.Di , teta(k) = acos((dot(Anisorest,Dk))/(norm(Anisorest)*norm(Dk))),
      std::vector<double> teta(3);
      teta[0] = std::acos(  (ShapeVectorResting[0][0]*AnisoRestLocal[0]+ShapeVectorResting[0][1]*AnisoRestLocal[1])/
                            std::sqrt(ShapeVectorResting[0][0]*ShapeVectorResting[0][0]+ShapeVectorResting[0][1]*ShapeVectorResting[0][1]+0.0000001) );
      
      teta[1] = std::acos(  (ShapeVectorResting[1][0]*AnisoRestLocal[0]+ShapeVectorResting[1][1]*AnisoRestLocal[1])/
                            std::sqrt(ShapeVectorResting[1][0]*ShapeVectorResting[1][0]+ShapeVectorResting[1][1]*ShapeVectorResting[1][1]+0.0000001) );
      
      teta[2] = std::acos(  (ShapeVectorResting[2][0]*AnisoRestLocal[0]+ShapeVectorResting[2][1]*AnisoRestLocal[1])/
                            std::sqrt(ShapeVectorResting[2][0]*ShapeVectorResting[2][0]+ShapeVectorResting[2][1]*ShapeVectorResting[2][1]+0.0000001) );

      //  std::cerr<< "cell "<< cellIndex<<"  numerator  " << (ShapeVectorResting[2][0]*AnisoRestLocal[0]+ShapeVectorResting[2][1]*AnisoRestLocal[1])/
      //                      std::sqrt(ShapeVectorResting[2][0]*ShapeVectorResting[2][0]+ShapeVectorResting[2][1]*ShapeVectorResting[2][1])  << std::endl;    
      //  std::cerr<< "cell "<< cellIndex<< " Q 0, 1 , 2:  " << Qa<<" , "<<Qb << " , " << Qc  << std::endl;
      //  std::cerr<< "cell "<< cellIndex<< " P 0, 1 , 2:  " << Pa<<" , "<<Pb << " , " << Pc  << std::endl;
      //  std::cerr<<"cell "<<cellIndex<< "   tet 0, 1 , 2:  " << teta[0]<<" , "<<teta[1] << " , " << teta[2]  << std::endl;
         
       // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> STRAIN and STRESS TENSOR (BEGIN) <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      
      // deformation gradiant tensor F =Sigma i=1,2,3 Qi x Di
      // strain tensor in resting shape E=0.5(FtF-I)
      // trE
      // B=FFt
      // axa (direct product of aniso vector in resting shape)
      // atEa
      // E(axa) and (axa)E
      double trE=( Delta[1]*cotan[0]+ Delta[2]*cotan[1]+Delta[0]*cotan[2])/(4*restingArea);
      
      double directAniso[2][2]={{AnisoRestLocal[0]*AnisoRestLocal[0],AnisoRestLocal[0]*AnisoRestLocal[1]},
                                {AnisoRestLocal[1]*AnisoRestLocal[0],AnisoRestLocal[1]*AnisoRestLocal[1]}};
      
      double positionLocal[3][2]={ {Qa , Qc}, 
                                   {0  , 0 },  
                                   {Qb , 0 }  };
      
      double DeformGrad[2][2]={{0,0},{0,0}}; // F= Qi x Di
      for ( int i=0 ; i<3 ; ++i ) {
        DeformGrad[0][0]=DeformGrad[0][0]+positionLocal[i][0]*ShapeVectorResting[i][0];
        DeformGrad[1][0]=DeformGrad[1][0]+positionLocal[i][1]*ShapeVectorResting[i][0];
        DeformGrad[0][1]=DeformGrad[0][1]+positionLocal[i][0]*ShapeVectorResting[i][1];
        DeformGrad[1][1]=DeformGrad[1][1]+positionLocal[i][1]*ShapeVectorResting[i][1];
      } 

      double LeftCauchy[2][2]; // B=FFt
      LeftCauchy[0][0]=DeformGrad[0][0]*DeformGrad[0][0]+DeformGrad[0][1]*DeformGrad[0][1];
      LeftCauchy[1][0]=DeformGrad[1][0]*DeformGrad[0][0]+DeformGrad[1][1]*DeformGrad[0][1];
      LeftCauchy[0][1]=DeformGrad[0][0]*DeformGrad[1][0]+DeformGrad[0][1]*DeformGrad[1][1];
      LeftCauchy[1][1]=DeformGrad[1][0]*DeformGrad[1][0]+DeformGrad[1][1]*DeformGrad[1][1];


      double Egreen[2][2];//E=0.5(C-I)
      Egreen[0][0]=0.5*(DeformGrad[0][0]*DeformGrad[0][0]+DeformGrad[1][0]*DeformGrad[1][0]-1);
      Egreen[1][0]=0.5*(DeformGrad[0][1]*DeformGrad[0][0]+DeformGrad[1][1]*DeformGrad[1][0]);
      Egreen[0][1]=0.5*(DeformGrad[0][0]*DeformGrad[0][1]+DeformGrad[1][0]*DeformGrad[1][1]);
      Egreen[1][1]=0.5*(DeformGrad[0][1]*DeformGrad[0][1]+DeformGrad[1][1]*DeformGrad[1][1]-1);

      
      double StrainAlmansi[2][2]; // e=0.5(1-B^-1)  True strain tensor
      temp=LeftCauchy[0][0]*LeftCauchy[1][1]-LeftCauchy[1][0]*LeftCauchy[0][1]; // det(B)
      StrainAlmansi[0][0]=0.5*(1-(LeftCauchy[1][1]/temp));
      StrainAlmansi[1][0]=0.5*LeftCauchy[1][0]/temp;
      StrainAlmansi[0][1]=0.5*LeftCauchy[0][1]/temp;  
      StrainAlmansi[1][1]=0.5*(1-(LeftCauchy[0][0]/temp));
      
      double atEa=AnisoRestLocal[0]*AnisoRestLocal[0]*Egreen[0][0]
                 +AnisoRestLocal[0]*AnisoRestLocal[1]*(Egreen[0][1]+Egreen[1][0])
                 +AnisoRestLocal[1]*AnisoRestLocal[1]*Egreen[1][1];
      
      double Eaa[2][2];
      Eaa[0][0]= Egreen[0][0]*directAniso[0][0]+Egreen[0][1]*directAniso[1][0];        
      Eaa[1][0]= Egreen[1][0]*directAniso[0][0]+Egreen[1][1]*directAniso[1][0];        
      Eaa[0][1]= Egreen[0][0]*directAniso[0][1]+Egreen[0][1]*directAniso[1][1];        
      Eaa[1][1]= Egreen[1][0]*directAniso[0][1]+Egreen[1][1]*directAniso[1][1];        

      double aaE[2][2];
      aaE[0][0]= directAniso[0][0]*Egreen[0][0]+directAniso[0][1]*Egreen[1][0];        
      aaE[1][0]= directAniso[1][0]*Egreen[0][0]+directAniso[1][1]*Egreen[1][0];        
      aaE[0][1]= directAniso[0][0]*Egreen[0][1]+directAniso[0][1]*Egreen[1][1];        
      aaE[1][1]= directAniso[1][0]*Egreen[0][1]+directAniso[1][1]*Egreen[1][1];        
      
      double B2[2][2];// LeftCauchy^2
      B2[0][0]=LeftCauchy[0][0]*LeftCauchy[0][0]+LeftCauchy[0][1]*LeftCauchy[1][0];
      B2[1][0]=LeftCauchy[1][0]*LeftCauchy[0][0]+LeftCauchy[1][1]*LeftCauchy[1][0];
      B2[0][1]=LeftCauchy[0][0]*LeftCauchy[0][1]+LeftCauchy[0][1]*LeftCauchy[1][1];
      B2[1][1]=LeftCauchy[1][0]*LeftCauchy[0][1]+LeftCauchy[1][1]*LeftCauchy[1][1];

      double Sigma[2][2]; // true stress tensor (isotropic term) based on lambdaT and mioT
      Sigma[0][0]=(Area/restingArea)*((lambdaT*trE-mioT/2)*LeftCauchy[0][0]+(mioT/2)*B2[0][0]);
      Sigma[1][0]=(Area/restingArea)*((lambdaT*trE-mioT/2)*LeftCauchy[1][0]+(mioT/2)*B2[1][0]);
      Sigma[0][1]=(Area/restingArea)*((lambdaT*trE-mioT/2)*LeftCauchy[0][1]+(mioT/2)*B2[0][1]);
      Sigma[1][1]=(Area/restingArea)*((lambdaT*trE-mioT/2)*LeftCauchy[1][1]+(mioT/2)*B2[1][1]);

      double deltaS[2][2];
      deltaS[0][0]=deltaLam*(trE*directAniso[0][0]+atEa)+(deltaMio/2)*(Eaa[0][0]+aaE[0][0])-(deltaLam+deltaMio)*atEa*directAniso[0][0];
      deltaS[1][0]=deltaLam*(trE*directAniso[1][0]     )+(deltaMio/2)*(Eaa[1][0]+aaE[1][0])-(deltaLam+deltaMio)*atEa*directAniso[1][0];
      deltaS[0][1]=deltaLam*(trE*directAniso[0][1]     )+(deltaMio/2)*(Eaa[0][1]+aaE[0][1])-(deltaLam+deltaMio)*atEa*directAniso[0][1];
      deltaS[1][1]=deltaLam*(trE*directAniso[1][1]+atEa)+(deltaMio/2)*(Eaa[1][1]+aaE[1][1])-(deltaLam+deltaMio)*atEa*directAniso[1][1];

      double deltaSFt[2][2];
      deltaSFt[0][0]=deltaS[0][0]*DeformGrad[0][0]+deltaS[0][1]*DeformGrad[0][1];
      deltaSFt[1][0]=deltaS[1][0]*DeformGrad[0][0]+deltaS[1][1]*DeformGrad[0][1];
      deltaSFt[0][1]=deltaS[0][0]*DeformGrad[1][0]+deltaS[0][1]*DeformGrad[1][1];
      deltaSFt[1][1]=deltaS[1][0]*DeformGrad[1][0]+deltaS[1][1]*DeformGrad[1][1];
      
      double deltaSigma[2][2];// true stress tensor (anisotropic correction term)deltaLambda and deltaMio (Longitudinal-Transverse)
      deltaSigma[0][0]=(Area/restingArea)*(DeformGrad[0][0]*deltaSFt[0][0]+DeformGrad[0][1]*deltaSFt[1][0]);
      deltaSigma[1][0]=(Area/restingArea)*(DeformGrad[1][0]*deltaSFt[0][0]+DeformGrad[1][1]*deltaSFt[1][0]);
      deltaSigma[0][1]=(Area/restingArea)*(DeformGrad[0][0]*deltaSFt[0][1]+DeformGrad[0][1]*deltaSFt[1][1]);
      deltaSigma[1][1]=(Area/restingArea)*(DeformGrad[1][0]*deltaSFt[0][1]+DeformGrad[1][1]*deltaSFt[1][1]);

      double StressTensor[3][3];
      StressTensor[0][0]=Sigma[0][0]+deltaSigma[0][0];
      StressTensor[1][0]=Sigma[1][0]+deltaSigma[1][0];
      StressTensor[0][1]=Sigma[0][1]+deltaSigma[0][1];
      StressTensor[1][1]=Sigma[1][1]+deltaSigma[1][1];


      // std::cerr <<"stress tensor " << std::endl;
      // std::cerr <<" Sxx  "<< StressTensor[0][0] <<" Sxy  "<< StressTensor[0][1] <<" Sxz  "<< StressTensor[0][2] << std::endl
      //           <<" Syx  "<< StressTensor[1][0] <<" Syy  "<< StressTensor[1][1] <<" Syz  "<< StressTensor[1][2] << std::endl
      //           <<" Szx  "<< StressTensor[2][0] <<" Szy  "<< StressTensor[2][1] <<" Szz  "<< StressTensor[2][2] << std::endl <<std::endl;


      //Shape vectors in Current shape (counterclockwise ordering of nodes/edges)     ShapeVectorCurrent[3][3]  calculated above   
      //.............................. ( or clockwise ordering of nodes/edges)
          

      //square of radius of circumstancing circle in resting shape
      //double Rcirc2Resting=(0.25*restingLength[0]*restingLength[1]*restingLength[2]/Area)*(0.25*restingLength[0]*restingLength[1]*restingLength[2]/Area);  
      
      double StrainTensor[3][3]; // there are other alternatives than StrainAlmansi for strain tensor
      StrainTensor[0][0]=StrainAlmansi[0][0];
      StrainTensor[1][0]=StrainAlmansi[1][0];
      StrainTensor[0][1]=StrainAlmansi[0][1];
      StrainTensor[1][1]=StrainAlmansi[1][1];


      StrainTensor[0][2]=0;  // adding 3rd dimension which is zero, the tensor is still in element plane
      StrainTensor[1][2]=0;
      StrainTensor[2][2]=0;
      StrainTensor[2][0]=0;
      StrainTensor[2][1]=0;
      
      StressTensor[0][2]=0;  // adding 3rd dimension which is zero, the tensor is still in element plane
      StressTensor[1][2]=0;
      StressTensor[2][2]=0;
      StressTensor[2][0]=0;
      StressTensor[2][1]=0;

      //rotation matrix to go to global coordinate system based on counterclockwise ordering;   rotation[3][3] calculated above  

      // rotating strain tensor to the global coordinate system
      double tempR[3][3]={{0,0,0},{0,0,0},{0,0,0}};
      for (int r=0 ; r<3 ; r++) 
        for (int s=0 ; s<3 ; s++) 
          for(int w=0 ; w<3 ; w++) 
            tempR[r][s]=tempR[r][s]+rotation[r][w]*StrainTensor[w][s];
          
      for (int r=0 ; r<3 ; r++) 
        for (int s=0 ; s<3 ; s++) {
          StrainTensor[r][s]=0;
          for(int w=0 ; w<3 ; w++) 
            StrainTensor[r][s]=StrainTensor[r][s]+tempR[r][w]*rotation[s][w]; 
        }

      // rotating stress tensor to the global coordinate system
      for (int r=0 ; r<3 ; r++) 
        for (int s=0 ; s<3 ; s++) 
          tempR[r][s]=0;
        
      for (int r=0 ; r<3 ; r++) 
        for (int s=0 ; s<3 ; s++) 
          for(int w=0 ; w<3 ; w++) 
            tempR[r][s]=tempR[r][s]+rotation[r][w]*StressTensor[w][s];
          
      for (int r=0 ; r<3 ; r++) 
        for (int s=0 ; s<3 ; s++) {
          StressTensor[r][s]=0;
          for(int w=0 ; w<3 ; w++) 
            StressTensor[r][s]=StressTensor[r][s]+tempR[r][w]*rotation[s][w]; 
        }

      // std::cerr <<"strain tensor in global coordinate system" << std::endl;
      // std::cerr <<" Sxx  "<< StrainTensor[0][0] <<" Sxy  "<< StrainTensor[0][1] <<" Sxz  "<< StrainTensor[0][2] << std::endl
      //           <<" Syx  "<< StrainTensor[1][0] <<" Syy  "<< StrainTensor[1][1] <<" Syz  "<< StrainTensor[1][2] << std::endl
      //           <<" Szx  "<< StrainTensor[2][0] <<" Szy  "<< StrainTensor[2][1] <<" Szz  "<< StrainTensor[2][2] << std::endl <<std::endl;

      //  std::cerr <<"stress tensor in global coordinate system" << std::endl;
      // std::cerr <<" Sxx  "<< StressTensor[0][0] <<" Sxy  "<< StressTensor[0][1] <<" Sxz  "<< StressTensor[0][2] << std::endl
      //           <<" Syx  "<< StressTensor[1][0] <<" Syy  "<< StressTensor[1][1] <<" Syz  "<< StressTensor[1][2] << std::endl
      //           <<" Szx  "<< StressTensor[2][0] <<" Szy  "<< StressTensor[2][1] <<" Szz  "<< StressTensor[2][2] << std::endl <<std::endl;


      // accumulating strain and stress tensors to be averaged later
      for (int r=0 ; r<3 ; r++) 
        for (int s=0 ; s<3 ; s++)    
          StrainCellGlobal[r][s]= StrainCellGlobal[r][s]+restingArea*StrainTensor[r][s];

      for (int r=0 ; r<3 ; r++) 
        for (int s=0 ; s<3 ; s++)    
          StressCellGlobal[r][s]= StressCellGlobal[r][s]+restingArea*StressTensor[r][s];
      
      TotalCellRestingArea=TotalCellRestingArea+restingArea;
    
      // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> STRAIN and STRESS TENSORS (END) <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

   
    //---- Anisotropic Correction Force-------------------------------
      double deltaF[3][3];
      double  Rcirc2=(0.25*length[0]*length[1]*length[2]/Area)*(0.25*length[0]*length[1]*length[2]/Area);  // square of radius of circumscribed circle in current shape
      
      double derIprim1[3][3];         // Invariants and their derivatives
      double derIprim4[3][3];
      double derIprim5[3][3];
      
      double DiDm;                    // inner products between shape vectors
      double DnDr;
      double DsDp;
      
      double QiQj;                    // inner products between position vectors of vertices
      double QrQs;
      
      double aDi;                     // inner products between shape vectors and anisotropy vector(direction)
      double aDj;
      double aDm;
      double aDp;
      double aDr;
      double aDs;
      double aDn;
      
      temp=restingLength[0];
      restingLength[0]=restingLength[1];
      restingLength[1]=restingLength[2];
      restingLength[2]=temp;
      
      temp=length[0];
      length[0]=length[1];
      length[1]=length[2];
      length[2]=temp;
      
      
      int kPerm;
      for ( int m=0 ; m<3 ; ++m )
        {
          for ( int coor=0 ; coor<3 ; ++coor ) 
            derIprim1[m][coor]=0;
          for ( int i=0 ; i<3 ; ++i )
            {            
              if ((i==0 && m==1)||(i==1 && m==0)) kPerm=2;
              if ((i==0 && m==2)||(i==2 && m==0)) kPerm=1;
              if ((i==1 && m==2)||(i==2 && m==1)) kPerm=0; 
              //else {
              //std::cerr << "mechanicalTRBS::derivs() kPerm not given a value..." << std::endl;
              //exit(-1);
              //}
              if (i!=m) DiDm=-0.5*cotan[kPerm]/restingArea;
              if (i==m) DiDm=0.25*restingLength[i]*restingLength[i] / (restingArea*restingArea);
              for ( int coor=0 ; coor<3 ; ++coor ) 
                derIprim1[m][coor]=derIprim1[m][coor]+2*DiDm*position[m][coor];
            }
        }
      
      double Iprim4=0;
      for ( int i=0 ; i<3 ; ++i )
        { 
          for ( int j=0 ; j<3 ; ++j )
            {
              if ((i==0 && j==1)||(i==1 && j==0)) kPerm=2; 
              if ((i==0 && j==2)||(i==2 && j==0)) kPerm=1; 
              if ((i==1 && j==2)||(i==2 && j==1)) kPerm=0;  
              if (i!=j) QiQj=Rcirc2-(length[kPerm]*length[kPerm])*0.5; 
              if (i==j) QiQj=Rcirc2;              
              aDi=0.5*cos(teta[i])*restingLength[i]/restingArea;
              aDj=0.5*cos(teta[j])*restingLength[j]/restingArea;
              Iprim4=Iprim4+ QiQj*aDi*aDj;
            }
        }
      
        for ( int p=0 ; p<3 ; ++p )
          {
            for ( int coor=0 ; coor<3 ; ++coor ) derIprim4[p][coor]=0;
            
            for ( int m=0 ; m<3 ; ++m )
              {
                aDm=0.5*cos(teta[m])*restingLength[m]/restingArea;
                for ( int coor=0 ; coor<3 ; ++coor ) derIprim4[p][coor]=derIprim4[p][coor]+aDm*position[m][coor];
              }
            aDp=0.5*cos(teta[p])*restingLength[p]/restingArea;
            for ( int coor=0 ; coor<3 ; ++coor ) derIprim4[p][coor]=2*aDp*derIprim4[p][coor];
          }
        
        
        for ( int p=0 ; p<3 ; ++p )                                                                             
          { 
            for ( int coor=0 ; coor<3 ; ++coor ) derIprim5[p][coor]=0;
            for ( int n=0 ; n<3 ; ++n )
              { 
                for ( int r=0 ; r<3 ; ++r )
                  { for ( int s=0 ; s<3 ; ++s )
                      {     
                        if ((r==0 && s==1)||(r==1 && s==0)) kPerm=2;
                        if ((r==0 && s==2)||(r==2 && s==0)) kPerm=1;
                        if ((r==1 && s==2)||(r==2 && s==1)) kPerm=0; 
                        //if ( s!=r ) QrQs=Rcirc2-(length[kPerm]*length[kPerm])*0.5;
                        //if ( s==r ) QrQs=Rcirc2;
                        QrQs=position[r][0]*position[s][0]+position[r][1]*position[s][1]+position[r][2]*position[s][2];
                        
                        if ((n==0 && r==1)||(n==1 && r==0)) kPerm=2; 
                        if ((n==0 && r==2)||(n==2 && r==0)) kPerm=1; 
                        if ((n==1 && r==2)||(n==2 && r==1)) kPerm=0;    
                        if ( n!=r )  DnDr=-0.5*cotan[kPerm]/restingArea;
                        if ( n==r )  DnDr=0.25*restingLength[n]*restingLength[n] / (restingArea*restingArea);
                        
                        if ((s==0 && p==1)||(s==1 && p==0)) kPerm=2; 
                        if ((s==0 && p==2)||(s==2 && p==0)) kPerm=1; 
                        if ((s==1 && p==2)||(s==2 && p==1)) kPerm=0;   
                        if ( s!=p ) DsDp=-0.5*cotan[kPerm]/restingArea;
                        if ( s==p ) DsDp=0.25*restingLength[s]*restingLength[s] / (restingArea*restingArea);
                        
                        aDs=0.5*cos(teta[s])*restingLength[s]/restingArea;
                        aDp=0.5*cos(teta[p])*restingLength[p]/restingArea;
                        aDr=0.5*cos(teta[r])*restingLength[r]/restingArea;
                        aDn=0.5*cos(teta[n])*restingLength[n]/restingArea;
                        
                        for ( int coor=0 ; coor<3 ; ++coor ) 
                          derIprim5[p][coor] = derIprim5[p][coor] +
                            2*(DnDr*aDs*aDp+DsDp*aDr*aDn)*QrQs*position[n][coor];
                      }
                  }
              }
          }		      
        
        double derI1[3][3];             // Invariants and their derivatives
        double derI4[3][3];
        double derI5[3][3];
        
        double I1=( Delta[1]*cotan[0]+ Delta[2]*cotan[1]+Delta[0]*cotan[2])/(4*restingArea);
        double I4=0.5*Iprim4-0.5;
        
        for ( int i=0 ; i<3 ; ++i ) 
          for ( int j=0 ; j<3 ; ++j ) {
            derI1[i][j]=0.5*derIprim1[i][j];
            derI4[i][j]=0.5*derIprim4[i][j];
            derI5[i][j]=0.25*derIprim5[i][j]-0.5*derIprim4[i][j];
          }   
        for ( int i=0 ; i<3 ; ++i ) 
          for ( int j=0 ; j<3 ; ++j )
            deltaF[i][j]=(-deltaLam*(I4*derI1[i][j]+I1*derI4[i][j])-deltaMio*derI5[i][j]+(deltaMio+deltaLam)*I4*derI4[i][j])*restingArea;
   
        //Forces of vertices   
        double Force[3][3];                                           
    
        Force[0][0]= (tensileStiffness[0]*Delta[0]+angularStiffness[1]*Delta[1]+angularStiffness[0]*Delta[2])*(position[1][0]-position[0][0])
          +(tensileStiffness[2]*Delta[2]+angularStiffness[2]*Delta[1]+angularStiffness[0]*Delta[0])*(position[2][0]-position[0][0])
          + deltaF[0][0]; 
        Force[0][1]= (tensileStiffness[0]*Delta[0]+angularStiffness[1]*Delta[1]+angularStiffness[0]*Delta[2])*(position[1][1]-position[0][1])
          +(tensileStiffness[2]*Delta[2]+angularStiffness[2]*Delta[1]+angularStiffness[0]*Delta[0])*(position[2][1]-position[0][1])
          + deltaF[0][1];  
        Force[0][2]= (tensileStiffness[0]*Delta[0]+angularStiffness[1]*Delta[1]+angularStiffness[0]*Delta[2])*(position[1][2]-position[0][2])
          +(tensileStiffness[2]*Delta[2]+angularStiffness[2]*Delta[1]+angularStiffness[0]*Delta[0])*(position[2][2]-position[0][2])
          + deltaF[0][2]; 
        Force[1][0]= (tensileStiffness[0]*Delta[0]+angularStiffness[0]*Delta[2]+angularStiffness[1]*Delta[1])*(position[0][0]-position[1][0])
          +(tensileStiffness[1]*Delta[1]+angularStiffness[2]*Delta[2]+angularStiffness[1]*Delta[0])*(position[2][0]-position[1][0])
          + deltaF[1][0];  
        Force[1][1]= (tensileStiffness[0]*Delta[0]+angularStiffness[0]*Delta[2]+angularStiffness[1]*Delta[1])*(position[0][1]-position[1][1])
          +(tensileStiffness[1]*Delta[1]+angularStiffness[2]*Delta[2]+angularStiffness[1]*Delta[0])*(position[2][1]-position[1][1])
          + deltaF[1][1];  
        Force[1][2]= (tensileStiffness[0]*Delta[0]+angularStiffness[0]*Delta[2]+angularStiffness[1]*Delta[1])*(position[0][2]-position[1][2])
          +(tensileStiffness[1]*Delta[1]+angularStiffness[2]*Delta[2]+angularStiffness[1]*Delta[0])*(position[2][2]-position[1][2])
          + deltaF[1][2];  
        
        Force[2][0]= (tensileStiffness[2]*Delta[2]+angularStiffness[0]*Delta[0]+angularStiffness[2]*Delta[1])*(position[0][0]-position[2][0])
          +(tensileStiffness[1]*Delta[1]+angularStiffness[1]*Delta[0]+angularStiffness[2]*Delta[2])*(position[1][0]-position[2][0])
          + deltaF[2][0];  
        Force[2][1]= (tensileStiffness[2]*Delta[2]+angularStiffness[0]*Delta[0]+angularStiffness[2]*Delta[1])*(position[0][1]-position[2][1])
          +(tensileStiffness[1]*Delta[1]+angularStiffness[1]*Delta[0]+angularStiffness[2]*Delta[2])*(position[1][1]-position[2][1])
          + deltaF[2][1];  
        Force[2][2]= (tensileStiffness[2]*Delta[2]+angularStiffness[0]*Delta[0]+angularStiffness[2]*Delta[1])*(position[0][2]-position[2][2])
          +(tensileStiffness[1]*Delta[1]+angularStiffness[1]*Delta[0]+angularStiffness[2]*Delta[2])*(position[1][2]-position[2][2])
          + deltaF[2][2];  
        
        // std::cerr << "Forces (cell " << cellIndex << "):" << std::endl 
        // 	      << Force[0][0] << " " << Force[0][1] << " " << Force[0][2] << std::endl
        // 	      << Force[1][0] << " " << Force[1][1] << " " << Force[1][2] << std::endl
        // 	      << Force[2][0] << " " << Force[2][1] << " " << Force[2][2] << std::endl;
        
        // adding TRBSMT forces to the total vertexDerives
        
        cellDerivs[cellIndex][comIndex  ] += Force[0][0];
        cellDerivs[cellIndex][comIndex+1] += Force[0][1];
        cellDerivs[cellIndex][comIndex+2] += Force[0][2];
        
        vertexDerivs[v2][0] += Force[1][0];
        vertexDerivs[v2][1] += Force[1][1];
        vertexDerivs[v2][2] += Force[1][2];
        
        vertexDerivs[v3][0] += Force[2][0];
        vertexDerivs[v3][1] += Force[2][1];
        vertexDerivs[v3][2] += Force[2][2];
        
    }
    for (int r=0 ; r<3 ; r++) 
      for (int s=0 ; s<3 ; s++)    
        StrainCellGlobal[r][s]= StrainCellGlobal[r][s]/TotalCellRestingArea; 
    
    
    
    for (int r=0 ; r<3 ; r++) 
      for (int s=0 ; s<3 ; s++)    
        StressCellGlobal[r][s]= StressCellGlobal[r][s]/TotalCellRestingArea; 
      
    
    // std::cerr <<" Sxx  "<< StrainCellGlobal[0][0] <<" Sxy  "<< StrainCellGlobal[0][1] <<" Sxz  "<< StrainCellGlobal[0][2] << std::endl
    //           <<" Syx  "<< StrainCellGlobal[1][0] <<" Syy  "<< StrainCellGlobal[1][1] <<" Syz  "<< StrainCellGlobal[1][2] << std::endl
    //           <<" Szx  "<< StrainCellGlobal[2][0] <<" Szy  "<< StrainCellGlobal[2][1] <<" Szz  "<< StrainCellGlobal[2][2] << std::endl <<std::endl;
    
    // std::cerr <<" Sxx  "<< StressCellGlobal[0][0] <<" Sxy  "<< StressCellGlobal[0][1] <<" Sxz  "<< StressCellGlobal[0][2] << std::endl
    //           <<" Syx  "<< StressCellGlobal[1][0] <<" Syy  "<< StressCellGlobal[1][1] <<" Syz  "<< StressCellGlobal[1][2] << std::endl
    //           <<" Szx  "<< StressCellGlobal[2][0] <<" Szy  "<< StressCellGlobal[2][1] <<" Szz  "<< StressCellGlobal[2][2] << std::endl <<std::endl;
    
   
    // eigenvalue/eigenvectors of averaged STRAIN and STRESS tensors in global coordinate system. (Jacobi method)

    // STRAIN:

    double eigenVectorStrain[3][3]={{1,0,0},{0,1,0},{0,0,1}};
    double pivot=1;
    double pi=3.1415;
    int I,J;
    double RotAngle,Si,Co;
    while (pivot>0.00001) {
      pivot=std::abs(StrainCellGlobal[1][0]);
      I=1;
      J=0;
      if (std::abs(StrainCellGlobal[2][0])>pivot) {
        pivot=std::abs(StrainCellGlobal[2][0]);
        I=2;
        J=0;
      }
      if (std::abs(StrainCellGlobal[2][1])>pivot) {
        pivot=std::abs(StrainCellGlobal[2][1]);
        I=2;
        J=1;
      }
      if (std::abs(StrainCellGlobal[I][I]-StrainCellGlobal[J][J])<0.00001) {
          RotAngle=pi/4;
      }            
      else {
        RotAngle=0.5*std::atan((2*StrainCellGlobal[I][J])/(StrainCellGlobal[J][J]-StrainCellGlobal[I][I]));
      }
        Si=std::sin(RotAngle);
        Co=std::cos(RotAngle);
        double tempRot[3][3]={{1,0,0},{0,1,0},{0,0,1}};
        tempRot[I][I]=Co;
        tempRot[J][J]=Co;
        tempRot[I][J]=Si;
        tempRot[J][I]=-Si;
        double tempStrain[3][3]={{0,0,0},{0,0,0},{0,0,0}};
        for (int r=0 ; r<3 ; r++) 
          for (int s=0 ; s<3 ; s++) 
            for(int w=0 ; w<3 ; w++) 
              tempStrain[r][s]=tempStrain[r][s]+StrainCellGlobal[r][w]*tempRot[w][s];
        
        for (int r=0 ; r<3 ; r++) 
          for (int s=0 ; s<3 ; s++) 
            StrainCellGlobal[r][s]=0;
         
        for (int r=0 ; r<3 ; r++) 
          for (int s=0 ; s<3 ; s++) 
            for(int w=0 ; w<3 ; w++) 
              StrainCellGlobal[r][s]=StrainCellGlobal[r][s]+tempRot[w][r]*tempStrain[w][s];
                
        for (int r=0 ; r<3 ; r++) 
          for (int s=0 ; s<3 ; s++) 
            tempStrain[r][s]=eigenVectorStrain[r][s];
        
        for (int r=0 ; r<3 ; r++) 
          for (int s=0 ; s<3 ; s++) 
            eigenVectorStrain[r][s]=0;
        
        for (int r=0 ; r<3 ; r++) 
          for (int s=0 ; s<3 ; s++) 
            for(int w=0 ; w<3 ; w++) 
              eigenVectorStrain[r][s]=eigenVectorStrain[r][s]+tempStrain[r][w]*tempRot[w][s];
            
    }
       
      
      // maximal strain direction
      double maximalStrainValue=StrainCellGlobal[0][0];
      int Istrain=0;
      if (StrainCellGlobal[1][1]>maximalStrainValue) 
        {
          maximalStrainValue=StrainCellGlobal[1][1];
          Istrain=1;
        }
      if (StrainCellGlobal[2][2]>maximalStrainValue) 
        {
          maximalStrainValue=StrainCellGlobal[2][2];
          Istrain=2;
        }
      // std::cerr<<"maximal Strain direction "<< eigenVectorStrain[0][Istrain] <<" "<< eigenVectorStrain[1][Istrain] <<" "<< eigenVectorStrain[2][Istrain] <<std::endl;  
      // std::cerr<<"maximal Strain value "<< maximalStrainValue <<std::endl;  
      
      // 2nd maximalstrain direction/value
      double maximalStrainValue2;
      int Istrain2,Istrain3;
      if (Istrain==0) {
        Istrain2=1;
        Istrain3=2;
      }
      if (Istrain==1) {
        Istrain2=0;
        Istrain3=2;
      }
      if (Istrain==2) {
        Istrain2=0;
        Istrain3=1;
      }
      if(StrainCellGlobal[Istrain3][Istrain3]>StrainCellGlobal[Istrain2][Istrain2]) {
        Istrain2=Istrain3;
      }
      maximalStrainValue2=StrainCellGlobal[Istrain2][Istrain2];


      // STRESS:
      
      double eigenVectorStress[3][3]={{1,0,0},{0,1,0},{0,0,1}};
      pivot=1;
      //double RotAngle,Si,Co;
      while (pivot>0.00001) {
        pivot=std::abs(StressCellGlobal[1][0]);
        I=1;
        J=0;
        if (std::abs(StressCellGlobal[2][0])>pivot) {
          pivot=std::abs(StressCellGlobal[2][0]);
          I=2;
          J=0;
        }
        if (std::abs(StressCellGlobal[2][1])>pivot) {
          pivot=std::abs(StressCellGlobal[2][1]);
          I=2;
          J=1;
        }
        if (std::abs(StressCellGlobal[I][I]-StressCellGlobal[J][J])<0.00001) {
          RotAngle=pi/4;
        }            
        else {
          RotAngle=0.5*std::atan((2*StressCellGlobal[I][J])/(StressCellGlobal[J][J]-StressCellGlobal[I][I]));
        }
        Si=std::sin(RotAngle);
        Co=std::cos(RotAngle);
        double tempRot[3][3]={{1,0,0},{0,1,0},{0,0,1}};
        tempRot[I][I]=Co;
        tempRot[J][J]=Co;
        tempRot[I][J]=Si;
        tempRot[J][I]=-Si;

        double tempStress[3][3]={{0,0,0},{0,0,0},{0,0,0}};
        for (int r=0 ; r<3 ; r++) 
          for (int s=0 ; s<3 ; s++) 
            for(int w=0 ; w<3 ; w++) 
              tempStress[r][s]=tempStress[r][s]+StressCellGlobal[r][w]*tempRot[w][s];
                        
        for (int r=0 ; r<3 ; r++) 
          for (int s=0 ; s<3 ; s++) 
            StressCellGlobal[r][s]=0;
          
        for (int r=0 ; r<3 ; r++) 
          for (int s=0 ; s<3 ; s++) 
            for(int w=0 ; w<3 ; w++) 
              StressCellGlobal[r][s]=StressCellGlobal[r][s]+tempRot[w][r]*tempStress[w][s];
        
        for (int r=0 ; r<3 ; r++) 
          for (int s=0 ; s<3 ; s++) 
            tempStress[r][s]=eigenVectorStress[r][s];
          
	for (size_t ii=0; ii<3; ++ii) {
	  for (size_t jj=0; jj<3; ++jj) {
	    eigenVectorStress[ii][jj] = 0.0;
	  }
	}
        
        for (int r=0 ; r<3 ; r++) 
          for (int s=0 ; s<3 ; s++) 
            for(int w=0 ; w<3 ; w++) 
              eigenVectorStress[r][s]=eigenVectorStress[r][s]+tempStress[r][w]*tempRot[w][s];
      }
      
      
      // maximal stress direction
      double maximalStressValue=StressCellGlobal[0][0];
      int Istress=0;
      if (StressCellGlobal[1][1]>maximalStressValue) 
        {
          maximalStressValue=StressCellGlobal[1][1];
          Istress=1;
        }
      if (StressCellGlobal[2][2]>maximalStressValue) 
        {
          maximalStressValue=StressCellGlobal[2][2];
          Istress=2;
        }
      // std::cerr<<"maximal Stress direction "<< eigenVectorStress[0][Istress] <<" "<< eigenVectorStress[1][Istress] <<" "<< eigenVectorStress[2][Istress] <<std::endl;  
      // std::cerr<<"maximal Stress value "<< maximalStressValue <<std::endl;  
      
      // 2nd maximalstress direction/value
      double maximalStressValue2;
      int Istress2,Istress3;
      if (Istress==0) {
        Istress2=1;
        Istress3=2;
      }
      if (Istress==1) {
        Istress2=0;
        Istress3=2;
      }
      if (Istress==2) {
        Istress2=0;
        Istress3=1;
      }
      if(StressCellGlobal[Istress3][Istress3]>StressCellGlobal[Istress2][Istress2]) {
        Istress2=Istress3;
      }
      maximalStressValue2=StressCellGlobal[Istress2][Istress2];

      
      // storing normal dirrection to  strain in cellData  ????????? not ready YET
     
      // normal to the cell plane in global direction is Zcurrent[], vector product gives the perpendicular strain direction
      // double NormalCurrent[3]; //= normal to the PCA plane in the current cell shape?????????????????????????????????????????????? not ready YET
      // double PerpStrain[3];
      // PerpStrain[0]=NormalCurrent[1]*eigenVectorStrain[2][Istrain]-NormalCurrent[2]*eigenVectorStrain[1][Istrain];
      // PerpStrain[1]=NormalCurrent[2]*eigenVectorStrain[0][Istrain]-NormalCurrent[0]*eigenVectorStrain[2][Istrain];
      // PerpStrain[2]=NormalCurrent[0]*eigenVectorStrain[1][Istrain]-NormalCurrent[1]*eigenVectorStrain[0][Istrain];
      // temp=std::sqrt(PerpStrain[0]*PerpStrain[0]+PerpStrain[1]*PerpStrain[1]+PerpStrain[2]*PerpStrain[2]);     

      // if(std::abs(temp)>0.001)
      //   { // if aniso vector is not perpendicular to the cell plane
      //   if (dimension==2)
      //     { cellData[cellIndex][0]=PerpStrain[0];        
      //       cellData[cellIndex][1]=PerpStrain[1];
      //     }
      //   if (dimension==3)
      //     { cellData[cellIndex][0]=PerpStrain[0];
      //       cellData[cellIndex][1]=PerpStrain[1];
      //       cellData[cellIndex][2]=PerpStrain[2];
      //       // cellData[cellIndex][3]=10*maximalStrainValue;  //NOTE maximal Strain and Stress Values can be used this is an option
      //     }
      // }
      

      if (numVariableIndexLevel()==6 && numVariableIndex(2) ) {// storing maximal strain
        if (dimension==2)
          {
            cellData[cellIndex][variableIndex(2,0)]  =eigenVectorStrain[0][Istrain];
            cellData[cellIndex][variableIndex(2,0)+1]=eigenVectorStrain[1][Istrain];
            cellData[cellIndex][variableIndex(2,0)+3]=maximalStrainValue;  //maximal Strain Value is stored after its eigenvector
          }
        if (dimension==3)
          {
            cellData[cellIndex][variableIndex(2,0)]  =eigenVectorStrain[0][Istrain];
            cellData[cellIndex][variableIndex(2,0)+1]=eigenVectorStrain[1][Istrain];
            cellData[cellIndex][variableIndex(2,0)+2]=eigenVectorStrain[2][Istrain];
            cellData[cellIndex][variableIndex(2,0)+3]=maximalStrainValue;  //maximal Strain Value is stored after its eigenvector
          }
      }

      
     
      if (numVariableIndexLevel()==6 && numVariableIndex(3) ) {//storing 2nd maximal strain
        if (dimension==2)
          {
            cellData[cellIndex][variableIndex(3,0)]  =eigenVectorStrain[0][Istrain2];
            cellData[cellIndex][variableIndex(3,0)+1]=eigenVectorStrain[1][Istrain2];
            cellData[cellIndex][variableIndex(3,0)+3]=maximalStrainValue2;  //2nd maximal Strain Value is stored after its eigenvector
          }
        if (dimension==3)
          {
            cellData[cellIndex][variableIndex(3,0)]  =eigenVectorStrain[0][Istrain2];
            cellData[cellIndex][variableIndex(3,0)+1]=eigenVectorStrain[1][Istrain2];
            cellData[cellIndex][variableIndex(3,0)+2]=eigenVectorStrain[2][Istrain2];
            cellData[cellIndex][variableIndex(3,0)+3]=maximalStrainValue2;  //2nd maximal Strain Value is stored after its eigenvector
          }
      }

      if (numVariableIndexLevel()==6 && numVariableIndex(4) ) { // storing maximal stress
	if (dimension==2)
	  {
	    cellData[cellIndex][variableIndex(4,0)]  =eigenVectorStress[0][Istress];
	    cellData[cellIndex][variableIndex(4,0)+1]=eigenVectorStress[1][Istress];
            cellData[cellIndex][variableIndex(4,0)+3]=maximalStressValue;  //maximal Stress Value is stored after its eigenvector
	  }
	if (dimension==3)
	  {
	    cellData[cellIndex][variableIndex(4,0)]  =eigenVectorStress[0][Istress];
	    cellData[cellIndex][variableIndex(4,0)+1]=eigenVectorStress[1][Istress];
	    cellData[cellIndex][variableIndex(4,0)+2]=eigenVectorStress[2][Istress];
            cellData[cellIndex][variableIndex(4,0)+3]=maximalStressValue;  //maximal Stress Value is stored after its eigenvector
	  }
      }

      if (numVariableIndexLevel()==6 && numVariableIndex(5) ) { // storing 2nd maximal stress
	if (dimension==2)
	  {
	    cellData[cellIndex][variableIndex(5,0)]  =eigenVectorStress[0][Istress2];
	    cellData[cellIndex][variableIndex(5,0)+1]=eigenVectorStress[1][Istress2];
            cellData[cellIndex][variableIndex(5,0)+3]=maximalStressValue2;  //2nd maximal Stress Value is stored after its eigenvector
	  }
	if (dimension==3)
	  {
	    cellData[cellIndex][variableIndex(5,0)]  =eigenVectorStress[0][Istress2];
	    cellData[cellIndex][variableIndex(5,0)+1]=eigenVectorStress[1][Istress2];
	    cellData[cellIndex][variableIndex(5,0)+2]=eigenVectorStress[2][Istress2];
            cellData[cellIndex][variableIndex(5,0)+3]=maximalStressValue2;  //2nd maximal Stress Value is stored after its eigenvector
	  }
      }

  }      
}     


void VertexFromTRBScenterTriangulationConcentrationHillMT::
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
  if (variableIndex(1,0) != numVariable) {
    std::cerr << "VertexFromTRBScenterTriangulationConcentrationHillMT::initiate() "
	      << "Wrong index given as start index for additional variables."
	      << std::endl;
    exit(-1);
  }
  size_t numCell = cellData.size();
  assert (numCell==T.numCell());
  std::vector<double> com(dimension);
  
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
  }
}


FiberModel::FiberModel(std::vector<double> &paraValue,
				     std::vector< std::vector<size_t> > &indValue)
{
  if (paraValue.size() != 8 && paraValue.size() != 9) {
    std::cerr << "FiberModel::FiberModel() " 
	      << "Uses eight or nine parameters: k_rate, equilibrium threshold , linear-hill flag, K_hill, n_hill, young_matrix, young_fiber, initialization flag and Poisson ratio. " << std::endl;
    exit(EXIT_FAILURE);
  }
  
  if (indValue.size() != 3 || indValue[0].size() != 1 ||  indValue[1].size() != 1 || indValue[2].size() != 1) {
    std::cerr << "FiberModel::FiberModel() " << std::endl
	      << "First level gives stress/strain anisotropy index." << std::endl
              << "Second level gives Young_Longitudinal index." << std::endl
	      << "Third level gives store index for velocity." << std::endl;
    exit(EXIT_FAILURE);
  }
  
  setId("FiberModel");
  setParameter(paraValue);
  setVariableIndex(indValue);
  
  std::vector<std::string> tmp(numParameter());
  tmp[0] = "k_rate";
  tmp[1] = "velocitythreshold";
  tmp[2] = "liniear-hill-flag";
  tmp[3] = "k_hill";
  tmp[4] = "n_hill";
  tmp[5] = "Y_matrix";
  tmp[6] = "Y_fiber";
  tmp[7] = "init_flag";
  if (parameter(7)==3)
    tmp[8] = "Poisson";

  setParameterId(tmp);
}

void FiberModel::initiate(Tissue &T,
                          DataMatrix &cellData,
                          DataMatrix &wallData,
                          DataMatrix &vertexData,
                          DataMatrix &cellDerivs,
                          DataMatrix &wallDerivs,
                          DataMatrix &vertexDerivs) {
  size_t numCell=cellData.size();
  size_t AnisoIndex=variableIndex(0,0);
  size_t YoungLIndex=variableIndex(1,0);
  double Kh=parameter(3);
  double Nh=parameter(4);
  double youngMatrix=parameter(5);
  double youngFiber=parameter(6);

  if (parameter(7)==1){
    for (size_t cellIndex=0; cellIndex<numCell; ++cellIndex) { // initiating with 0 anisotropy and isotropic material
      cellData[cellIndex][AnisoIndex] = 0;
      cellData[cellIndex][YoungLIndex] = youngMatrix+0.5*youngFiber; // youngL = youngMatrix + 0.5*youngFiber;
      //std::cerr<< cellData[cellIndex][variableIndex(1,0)] << std::endl;
    }
  }
  if (parameter(7)==2){
    for (size_t cellIndex=0; cellIndex<numCell; ++cellIndex) { // initiating with 0 anisotropy and isotropic material
      double anisotropy=cellData[cellIndex][AnisoIndex];
      if ( parameter(2)==0 ) // linear
	cellData[cellIndex][YoungLIndex]=youngMatrix+0.5*(1+anisotropy)* youngFiber;
      
      if ( parameter(2)==1 ) // Hill
	cellData[cellIndex][YoungLIndex] =  youngMatrix+ 
	  0.5*(1+(std::pow(anisotropy,Nh) /(std::pow((1-anisotropy),Nh)*std::pow(Kh,Nh)
					    +std::pow(anisotropy,Nh))))* youngFiber;

      if ( parameter(2)==2 )  // minimum energy
	cellData[cellIndex][YoungLIndex] =  youngMatrix+ 
	  0.5*(1+(std::pow(anisotropy,Nh) /(std::pow((1-anisotropy),Nh)*std::pow(Kh,Nh)
					    +std::pow(anisotropy,Nh))))* youngFiber; 
      //std::cerr<< cellData[cellIndex][variableIndex(1,0)] << std::endl;
    }
  }

 
  
}

void FiberModel::derivs(Tissue &T,
                        DataMatrix &cellData,
                        DataMatrix &wallData,
                        DataMatrix &vertexData,
                        DataMatrix &cellDerivs,
                        DataMatrix &wallDerivs,
                        DataMatrix &vertexDerivs ) {}


void FiberModel::update(Tissue &T,
                        DataMatrix &cellData,
                        DataMatrix &wallData,
                        DataMatrix &vertexData, 
                        double h) 
{
  size_t numCell=cellData.size();
  size_t AnisoIndex=variableIndex(0,0);
  size_t YoungLIndex=variableIndex(1,0);
  size_t velocityIndex=variableIndex(2,0);
  double Kh=parameter(3);
  double Nh=parameter(4);
  double youngMatrix=parameter(5);
  double youngFiber=parameter(6);
  double poisson=parameter(7);

  if (parameter(0)==0.0)
    return;
  for (size_t cellIndex=0; cellIndex<numCell; ++cellIndex) {
    double anisotropy=cellData[cellIndex][AnisoIndex];
    if ( parameter(2)==0 
	 && cellData[cellIndex][velocityIndex] < parameter(1) 
	 && cellData[cellIndex][YoungLIndex] < youngMatrix+0.5*(1+anisotropy)* youngFiber
	 && cellData[cellIndex][YoungLIndex] < youngMatrix+youngFiber // not to exceed maximum when using absolute stress anisotropy 
	 ){ // linear Fiber model
       cellData[cellIndex][YoungLIndex] += parameter(0)*h*0.5*anisotropy * youngFiber;
      //cellData[cellIndex][YoungLIndex] = youngMatrix+0.5*(1+anisotropy)* youngFiber;
    }
    if ( parameter(2)==1 
	 && cellData[cellIndex][velocityIndex] < parameter(1) 
	 && cellData[cellIndex][YoungLIndex] < 
	 youngMatrix+
	 0.5*(1+(std::pow(anisotropy,Nh)
		 /(std::pow((1-anisotropy),Nh)*std::pow(Kh,Nh)+std::pow(anisotropy,Nh))))* youngFiber
	 && cellData[cellIndex][YoungLIndex] < youngMatrix+youngFiber // not to exceed maximum when using absolute stress anisotropy 
	 ){ // Hill-like Fiber model

      // cellData[cellIndex][] += parameter(0)*h*(cellData[cellIndex][]-cellData[cellIndex][]);
      // cellData[cellIndex][YoungLIndex] = youngMatrix+
      //   0.5*(1+(std::pow(anisotropy,Nh)/(std::pow((1-anisotropy),Nh)*std::pow(Kh,Nh)+std::pow(anisotropy,Nh))))* youngFiber;
      cellData[cellIndex][YoungLIndex] 
	+= parameter(0)*h*
	(std::pow(anisotropy,Nh)/(std::pow(((1-anisotropy)*Kh),Nh)+std::pow(anisotropy,Nh)))* youngFiber*0.5;
      
      //std::cerr<< h << std::endl;
      // std::cerr<< cellData[cellIndex][YoungLIndex] << std::endl;
    }
  }
}









// //functions including bending energy

// VertexFromTRBSMTbending::
// VertexFromTRBSMTbending(std::vector<double> &paraValue, 
// 	       std::vector< std::vector<size_t> > 
// 	       &indValue ) 
// {  
//   // Do some checks on the parameters and variable indeces
//   if( paraValue.size()!=10 ) {
//     std::cerr << "VertexFromTRBSMTbending::"
// 	      << "VertexFromTRBSMTbending() "
//               << "Uses ten parameters: "
//               << "0,1 : young modulus(matrix and fibre) " 
//               << "2,3 : poisson ratio (longitudinal (MT) and transverse directions)"
// 	      << "4 : MF flag(0 constant material anisotropy ,1: material anisotropy via FiberModel "
//               << "5 : K_Hill, 6: n_Hill  " 
//               << "7 : flag(0: plane strain, 1: plane stress) " 
//               << "8 : MT direction angle"
//               << "9 : flag(0: for no feedback or direct feedback by indices,"
//               << "         1: for MT direction from 6th parameter TETA, 2: force to Stress,  "
//               << "         3: force to Strain , 4: force to perp-strain "<< std::endl;
//     exit(0);
//   }

//   if( (indValue.size()!=1 && indValue.size()!=3) ||
//       indValue[0].size()!=9 ||
//       (indValue.size()==3 && (indValue[1].size()!=0 && indValue[1].size()!=1 && indValue[1].size()!=2 && indValue[1].size()!=3 )) ||
//       (indValue.size()==3 && (indValue[2].size()!=0 && indValue[2].size()!=1 && indValue[2].size()!=2)) 
//       ) { 
//     std::cerr << "VertexFromTRBSMTbending::"
// 	      << "VertexFromTRBSMTbending() "
//               << "resting length index and MT direction initial index and strain and stress anisotropy indices and " 
//               << "indices for storing area ratio, isotropic energy and anisotropic energy and young_fiber and MTstress given in first level." 
//               << "Optionally two additional levels can be given where the strain and perpendicular direction to strain in cell plane and 2nd strain" 
//               << "directions/values(dx dy dz value) can be stored at given indices at second level."
//               << "If no index given at second level, strain will not be stored, if one index given strain will be stored"
//               << "and if two indices are given maximal and perpendicular to strain will be stored at first and second"
//               << "indices at second level respectively"
//               << "and if 3 indices are given maximal and perpendicular to strain and 2nd strain will be stored at 1st, 2nd and 3rd "
//               << "indices at second level respectively"
//               << "If no index given at 3rd level, stress will not be stored, if one index given stress will be stored"
//               << "and if two indices are given maximal and 2nd stress will be stored at first and second"
//               << "indices at 3rd level respectively"
//               << std::endl;
//     exit(0);
//   }
  
//   // Set the variable values
//   setId("VertexFromTRBSMTbending");
//   setParameter(paraValue);  
//   setVariableIndex(indValue);
  
  

//   // Set the parameter identities
//   std::vector<std::string> tmp( numParameter() );
//   tmp[0] = "Y_mod_M";   // Matrix Young modulus
//   tmp[1] = "Y_mod_F";   // Fiber Young modulus
//   tmp[2] = "P_ratio_L"; // Longitudinal Poisson ratio
//   tmp[3] = "P_ratio_T"; // Transverse Poisson ratio
//   tmp[4] = "MF flag";
//   tmp[5] = "neighborweight";
//   tmp[6] = "parameter6";
//   tmp[7] = "Strain-Stress flag";
//   tmp[8] = "TETA anisotropy";
//   tmp[9] = "MT update flag";

//   setParameterId( tmp );
 
//   if( parameter(2)<0 || parameter(2)>=0.5 || parameter(3)<0 || parameter(3)>=0.5 ) {
//     std::cerr << "VertexFromTRBSMTbending::"
//  	      << "VertexFromTRBSMTbending() "
//  	      << "poisson ratios must be 0 <= p < 0.5 " << std::endl;
//     exit(0);
//   }

//   if( parameter(4)!=0 && parameter(4)!=1 ) {
//     std::cerr << " VertexFromTRBSMTbending::"
// 	      << " VertexFromTRBSMTbending() "
// 	      << " 5th parameter must be 0 or 1 " 
//               << " 0: constant material anisotropy, 1: material anisotropy via FiberModel " << std::endl;
//     exit(0);
//     exit(0);
//   }

//   // if( parameter(5)<0 || parameter(5)>=1 ) {
//   //   std::cerr << " VertexFromTRBSMTbending::"
//   //             << " VertexFromTRBSMTbending() "
//   //             << " K parameter in Hill function must be 0 < K < 1 " << std::endl;
//   //   exit(0);
//   // }

//   if( parameter(7)!=0 && parameter(7)!=1 ) {
//     std::cerr << " VertexFromTRBSMTbending::"
// 	      << " VertexFromTRBSMTbending() "
// 	      << " 6th parameter must be 0 or 1(0:plane strain, 1:plane stress) " << std::endl;
//     exit(0);
//   }

//   if( parameter(9)!=0 && parameter(9)!=1 && parameter(9)!=2 && parameter(9)!=3 && parameter(9)!=4) {
//     std::cerr << " VertexFromTRBSMTbending::"
//               << " VertexFromTRBSMTbending() "
//               << " 8th parameter must be 0/1/2/3/4"
//               << " 0: for no feedback or direct feedback by indices "
//               << " 1: for MT direction from 7th parameter TETA "
//               << " 2: force to Stress "
//               << " 3: force to S1train "
//               << " 4: force to perp-strain " << std::endl;
//     exit(0);
//   }
// }

// void VertexFromTRBSMTbending::
// derivs(Tissue &T,
//        DataMatrix &cellData,
//        DataMatrix &wallData,
//        DataMatrix &vertexData,
//        DataMatrix &cellDerivs,
//        DataMatrix &wallDerivs,
//        DataMatrix &vertexDerivs ) {
  
//   //Do the update for each cell
//   size_t dimension = 3;
//   assert (dimension==vertexData[0].size());
//   size_t numCells = T.numCell();
//   size_t wallLengthIndex = variableIndex(0,0);
//   size_t numWalls = 3; // defined only for triangles 
  
//   double TotalVolume=0;
//   double deltaVolume=0;


//   for( size_t cellIndex=0 ; cellIndex<numCells ; ++cellIndex ) {
//     if( T.cell(cellIndex).numWall() != numWalls ) {
//       std::cerr << "VertexFromTRBSMTbending::derivs() only defined for triangular cells."
// 		<< " Not for cells with " << T.cell(cellIndex).numWall() << " walls!"
// 		<< std::endl;
//       exit(-1);
//     }
    
//     double youngMatrix= parameter(0);    
//     double youngFiber = parameter(1); 
//     double poissonL   = parameter(2);    
//     double poissonT   = parameter(3);
//     double TETA= parameter(8); 

//     size_t v1 = T.cell(cellIndex).vertex(0)->index();
//     size_t v2 = T.cell(cellIndex).vertex(1)->index();
//     size_t v3 = T.cell(cellIndex).vertex(2)->index();
//     size_t w1 = T.cell(cellIndex).wall(0)->index();
//     size_t w2 = T.cell(cellIndex).wall(1)->index();
//     size_t w3 = T.cell(cellIndex).wall(2)->index();

//     //std::cerr<< "cell "<< cellIndex<< " vertices  "<< v1<<" "<< v2 << " "<< v3 << " walls  "<< w1 <<" "<< w2 << " "<< w3<< std::endl;
    
   

//     double youngL;
//     double youngT;    

//     if( parameter(4)==0 ){ // constant anisotropic material
//       youngL = youngMatrix+youngFiber;
//       youngT = youngMatrix; 
//     }
//     if( parameter(4)==1){  // material anisotropy via FiberModel
//       youngL = cellData[cellIndex][variableIndex(0,7)]; 
//       youngT = 2*youngMatrix+youngFiber-youngL; 
//     }
    
   
//     double lambdaL, mioL, lambdaT, mioT;
    
//     if (parameter(7)==0){      
//       // Lame coefficients based on plane strain (for 3D 0<poisson<0.5)
//       lambdaL=youngL*poissonL/((1+poissonL)*(1-2*poissonL));
//       mioL=youngL/(2*(1+poissonL));
//       lambdaT=youngT*poissonT/((1+poissonT)*(1-2*poissonT));
//       mioT=youngT/(2*(1+poissonT));
//     } 
//     else{      
//       // Lame coefficients based on plane stress (for 3D 0<poisson<0.5)
//       lambdaL=youngL*poissonL/(1-poissonL*poissonL);
//       mioL=youngL/(2*(1+poissonL));
//       lambdaT=youngT*poissonT/(1-poissonT*poissonT);
//       mioT=youngT/(2*(1+poissonT));
//     }
    
//     // Lame coefficients based on delin. paper (for 2D 0<poisson<1)
//     // double lambdaL=youngL*poissonL/(1-poissonL*poissonL);
//     // double mioL=youngL/(1+poissonL);
//     // double lambdaT=youngT*poissonT/(1-poissonT*poissonT);
//     // double mioT=youngT/(1+poissonT);
     
 
    
//     double EnergyIso=0;                      
//     double EnergyAniso=0;
//     double strainZ=0;


//     double TETA1=parameter(5);
//     if ( parameter(9)==1 ) {  // aniso direction from TETA 
//       if ( cellIndex==0 || cellIndex==2) {
// 	cellData[cellIndex][variableIndex(0,1)]=std::cos(TETA1);  
// 	cellData[cellIndex][variableIndex(0,1)+1]=std::sin(TETA1);
// 	cellData[cellIndex][variableIndex(0,1)+2]=0;
//       }
//       if ( cellIndex==1 || cellIndex==3) {
// 	cellData[cellIndex][variableIndex(0,1)]=std::cos(TETA);  
// 	cellData[cellIndex][variableIndex(0,1)+1]=std::sin(TETA);
// 	cellData[cellIndex][variableIndex(0,1)+2]=0;
//       }
//     }
//     // if ( parameter(9)==1 ) {  // aniso direction from TETA 
//     //   cellData[cellIndex][variableIndex(0,1)]=std::cos(TETA);  
//     //   cellData[cellIndex][variableIndex(0,1)+1]=std::sin(TETA);
//     //   cellData[cellIndex][variableIndex(0,1)+2]=0;
//     // }
    
//     // Aniso vector in current shape in global coordinate system
//     double AnisoCurrGlob[3];
//     AnisoCurrGlob[0] = cellData[cellIndex][variableIndex(0,1)];
//     AnisoCurrGlob[1] = cellData[cellIndex][variableIndex(0,1)+1];
//     AnisoCurrGlob[2] = cellData[cellIndex][variableIndex(0,1)+2];
    

//     std::vector<double> restingLength(numWalls);
//     restingLength[0] = wallData[w1][wallLengthIndex];
//     restingLength[1] = wallData[w2][wallLengthIndex];
//     restingLength[2] = wallData[w3][wallLengthIndex];

//     DataMatrix position(3,vertexData[v1]);
//     position[1] = vertexData[v2];
//     position[2] = vertexData[v3];
//     //position[0][2] z for vertex 1 (of the cell)
   
//     std::vector<double> length(numWalls);
//     length[0] = T.wall(w1).lengthFromVertexPosition(vertexData);
//     length[1] = T.wall(w2).lengthFromVertexPosition(vertexData);
//     length[2] = T.wall(w3).lengthFromVertexPosition(vertexData);
    
//     //Anisotropic Correction is based on difference between Lam Coefficients of Longitudinal and Transverse dirrections:
//     double deltaLam=lambdaL-lambdaT;
//     double deltaMio=mioL-mioT;
//     double deltaMio1=(youngL-youngT)/2;
//     // double deltaLam=AnisoMeasure*(lambdaL-lambdaT);
//     // double deltaMio=AnisoMeasure*(mioL-mioT);
    
//     //Resting area of the element (using Heron's formula)                                      
//     double restingArea=std::sqrt( ( restingLength[0]+restingLength[1]+restingLength[2])*
//                                   (-restingLength[0]+restingLength[1]+restingLength[2])*
//                                   ( restingLength[0]-restingLength[1]+restingLength[2])*
//                                   ( restingLength[0]+restingLength[1]-restingLength[2])  )*0.25;
    
//     //Angles of the element ( assuming the order: 0,L0,1,L1,2,L2 clockwise )
//     std::vector<double> Angle(3);
//     // can be ommited by cotan(A)=.25*sqrt(4*b*b*c*c/K-1)
//     Angle[0]=std::acos(  (restingLength[0]*restingLength[0]+restingLength[2]*restingLength[2]-restingLength[1]*restingLength[1])/
//                          (restingLength[0]*restingLength[2]*2)    );
//     Angle[1]=std::acos(  (restingLength[0]*restingLength[0]+restingLength[1]*restingLength[1]-restingLength[2]*restingLength[2])/
//                          (restingLength[0]*restingLength[1]*2)    );
//     Angle[2]=std::acos(  (restingLength[1]*restingLength[1]+restingLength[2]*restingLength[2]-restingLength[0]*restingLength[0])/
//                          (restingLength[1]*restingLength[2]*2)    );
    
//     //Tensile Stiffness
//     double tensileStiffness[3];
//     double temp = 1.0/(restingArea*16);                                      
//     std::vector<double> cotan(3);
//     cotan[0] = 1.0/std::tan(Angle[0]);
//     cotan[1] = 1.0/std::tan(Angle[1]);
//     cotan[2] = 1.0/std::tan(Angle[2]);    
//     //the force is calculated based on Transverse coefficients
//     //Longitudinal coefficients are considered in deltaF
    
//     tensileStiffness[0]=(2*cotan[2]*cotan[2]*(lambdaT+2*mioT)+2*mioT)*temp;
//     tensileStiffness[1]=(2*cotan[0]*cotan[0]*(lambdaT+2*mioT)+2*mioT)*temp;
//     tensileStiffness[2]=(2*cotan[1]*cotan[1]*(lambdaT+2*mioT)+2*mioT)*temp;
    
//     //Angular Stiffness
//     std::vector<double> angularStiffness(3);
//     angularStiffness[0]=(2*cotan[1]*cotan[2]*(lambdaT+2*mioT)-2*mioT)*temp;                          
//     angularStiffness[1]=(2*cotan[0]*cotan[2]*(lambdaT+2*mioT)-2*mioT)*temp;
//     angularStiffness[2]=(2*cotan[0]*cotan[1]*(lambdaT+2*mioT)-2*mioT)*temp;
    
//     //Calculate biquadratic strains  
//     std::vector<double> Delta(3);
//     Delta[0]=(length[0])*(length[0])-(restingLength[0])*(restingLength[0]);
//     Delta[1]=(length[1])*(length[1])-(restingLength[1])*(restingLength[1]);
//     Delta[2]=(length[2])*(length[2])-(restingLength[2])*(restingLength[2]);

//     //Area of the element (using Heron's formula)                                      
//     double Area=std::sqrt( ( length[0]+length[1]+length[2])*
//                            (-length[0]+length[1]+length[2])*
//                            ( length[0]-length[1]+length[2])*
//                            ( length[0]+length[1]-length[2])  )*0.25;


//     // calculating the angles between shape vectors and anisotropy direction in resting shape when anisotropy vector is provided in current shape
    
//     //Current shape local coordinate of the element  (counterclockwise ordering of nodes/edges)
//       double CurrentAngle1=std::acos(  (length[0]*length[0]+length[1]*length[1]-length[2]*length[2])/
//                                        (length[0]*length[1]*2)    );

//       double Qa=std::cos(CurrentAngle1)*length[0];
//       double Qc=std::sin(CurrentAngle1)*length[0];
//       double Qb=length[1];
//       // shape vector matrix = inverse of coordinate matrix ( only first two elements i.e. ShapeVector[3][2] )      
//       double ShapeVectorCurrent[3][3]={ {  0   ,       1/Qc      , 0 }, 
//                                         {-1/Qb , (Qa-Qb)/(Qb*Qc) , 1 },       
//                                         { 1/Qb ,     -Qa/(Qb*Qc) , 0 }  };
            
   
//       // Local coordinates of the resting shape ( counterclockwise )
//       double RestingAngle1=std::acos(  (restingLength[0]*restingLength[0]+restingLength[1]*restingLength[1]-restingLength[2]*restingLength[2])/
//                                        (restingLength[0]*restingLength[1]*2)    );

//       double Pa=std::cos(RestingAngle1)*restingLength[0];
//       double Pc=std::sin(RestingAngle1)*restingLength[0];
//       double Pb=restingLength[1];

//       // shape vector matrix in resting shape in local coordinate system  = inverse of coordinate matrix ( only first two elements i.e. ShapeVectorResting[3][2] )      
//       double ShapeVectorResting[3][3]={ {  0   ,       1/Pc      , 0 }, 
//                                         {-1/Pb , (Pa-Pb)/(Pb*Pc) , 1 },       
//                                         { 1/Pb ,     -Pa/(Pb*Pc) , 0 }  };
   
//       // Rotation Matrix for changing coordinate systems for both Local to Global( Strain Tensor) and Global to Local( Aniso Vector in the current shape)
//       double rotation[3][3];  

//       double tempA=std::sqrt((position[2][0]-position[1][0])*(position[2][0]-position[1][0])+
//                              (position[2][1]-position[1][1])*(position[2][1]-position[1][1])+
//                              (position[2][2]-position[1][2])*(position[2][2]-position[1][2])  );

//       double tempB=std::sqrt((position[0][0]-position[1][0])*(position[0][0]-position[1][0])+
//                              (position[0][1]-position[1][1])*(position[0][1]-position[1][1])+
//                              (position[0][2]-position[1][2])*(position[0][2]-position[1][2])  );
      
//       double Xcurrent[3];      
//       Xcurrent[0]= (position[2][0]-position[1][0])/tempA;
//       Xcurrent[1]= (position[2][1]-position[1][1])/tempA;
//       Xcurrent[2]= (position[2][2]-position[1][2])/tempA;
      
//       double Bcurrent[3];      
//       Bcurrent[0]= (position[0][0]-position[1][0])/tempB;
//       Bcurrent[1]= (position[0][1]-position[1][1])/tempB;
//       Bcurrent[2]= (position[0][2]-position[1][2])/tempB;

//       double Zcurrent[3];      
//       Zcurrent[0]= Xcurrent[1]*Bcurrent[2]-Xcurrent[2]*Bcurrent[1];
//       Zcurrent[1]= Xcurrent[2]*Bcurrent[0]-Xcurrent[0]*Bcurrent[2];
//       Zcurrent[2]= Xcurrent[0]*Bcurrent[1]-Xcurrent[1]*Bcurrent[0];
      
//       tempA=std:: sqrt(Zcurrent[0]*Zcurrent[0]+Zcurrent[1]*Zcurrent[1]+Zcurrent[2]*Zcurrent[2]);
//       Zcurrent[0]=Zcurrent[0]/tempA;
//       Zcurrent[1]=Zcurrent[1]/tempA;
//       Zcurrent[2]=Zcurrent[2]/tempA;

//       double Ycurrent[3];      
//       Ycurrent[0]= Zcurrent[1]*Xcurrent[2]-Zcurrent[2]*Xcurrent[1];
//       Ycurrent[1]= Zcurrent[2]*Xcurrent[0]-Zcurrent[0]*Xcurrent[2];
//       Ycurrent[2]= Zcurrent[0]*Xcurrent[1]-Zcurrent[1]*Xcurrent[0];


//       rotation[0][0]=Xcurrent[0];
//       rotation[1][0]=Xcurrent[1];
//       rotation[2][0]=Xcurrent[2];

//       rotation[0][1]=Ycurrent[0];
//       rotation[1][1]=Ycurrent[1];
//       rotation[2][1]=Ycurrent[2];

//       rotation[0][2]=Zcurrent[0];
//       rotation[1][2]=Zcurrent[1];
//       rotation[2][2]=Zcurrent[2];      
       
//       // rotating the anisotropy vector from global coordinate system to the local one in the current shape
//       double AnisoCurrLocal[3];
//       AnisoCurrLocal[0]=rotation[0][0]*AnisoCurrGlob[0]+rotation[1][0]*AnisoCurrGlob[1]+rotation[2][0]*AnisoCurrGlob[2];
//       AnisoCurrLocal[1]=rotation[0][1]*AnisoCurrGlob[0]+rotation[1][1]*AnisoCurrGlob[1]+rotation[2][1]*AnisoCurrGlob[2];
//       AnisoCurrLocal[2]=rotation[0][2]*AnisoCurrGlob[0]+rotation[1][2]*AnisoCurrGlob[1]+rotation[2][2]*AnisoCurrGlob[2];
     


//       double positionLocal[3][2]={ {Qa , Qc}, 
//                                    {0  , 0 },  
//                                    {Qb , 0 }  };
      
//       double DeformGrad[2][2]={{0,0},{0,0}}; // F= Qi x Di
//       for ( int i=0 ; i<3 ; ++i ) {
//         DeformGrad[0][0]=DeformGrad[0][0]+positionLocal[i][0]*ShapeVectorResting[i][0];
//         DeformGrad[1][0]=DeformGrad[1][0]+positionLocal[i][1]*ShapeVectorResting[i][0];
//         DeformGrad[0][1]=DeformGrad[0][1]+positionLocal[i][0]*ShapeVectorResting[i][1];
//         DeformGrad[1][1]=DeformGrad[1][1]+positionLocal[i][1]*ShapeVectorResting[i][1];
//       }      

//       double AnisoRestLocal[3]={0,0,0};

//       AnisoRestLocal[0]=DeformGrad[0][0]*AnisoCurrLocal[0]+DeformGrad[1][0]*AnisoCurrLocal[1];
//       AnisoRestLocal[1]=DeformGrad[0][1]*AnisoCurrLocal[0]+DeformGrad[1][1]*AnisoCurrLocal[1];
     
     
//       double tempAn=std::sqrt(AnisoRestLocal[0]*AnisoRestLocal[0]+AnisoRestLocal[1]*AnisoRestLocal[1]+AnisoRestLocal[2]*AnisoRestLocal[2]);
//       AnisoRestLocal[0]/=tempAn;
//       AnisoRestLocal[1]/=tempAn;
      




//       double AnisoMeasure=std::sqrt(AnisoRestLocal[0]*AnisoRestLocal[0]+AnisoRestLocal[1]*AnisoRestLocal[1]);
//       // std::cerr<< "cell "<< cellIndex<<" AnisoMeasure "<< AnisoMeasure<<std::endl;
      
//       // choosing a random normalized dirrection for anisotropy if AnisoVector is close to perpendicular to the cell plane
//       if ( AnisoMeasure<0.0001) {
//         double randomAngle=((rand()%360)*2*3.14159265)/360;
        
//         AnisoRestLocal[0]=std::cos(randomAngle);
//         AnisoRestLocal[1]=std::sin(randomAngle);
//       }  
//       else {// normalizing the anisoVector if it is not random
//         AnisoRestLocal[0]=AnisoRestLocal[0]/AnisoMeasure;
//         AnisoRestLocal[1]=AnisoRestLocal[1]/AnisoMeasure;        
//       }
      
      
      
//       //Angles between anisotropy vector and shape vectors for calculating the terms like a.Di , teta(k) = acos((dot(Anisorest,Dk))/(norm(Anisorest)*norm(Dk))),
//       std::vector<double> teta(3);
//       teta[0] = std::acos(  (ShapeVectorResting[0][0]*AnisoRestLocal[0]+ShapeVectorResting[0][1]*AnisoRestLocal[1])/
//                             std::sqrt(ShapeVectorResting[0][0]*ShapeVectorResting[0][0]+ShapeVectorResting[0][1]*ShapeVectorResting[0][1]+0.0000001) );
      
//       teta[1] = std::acos(  (ShapeVectorResting[1][0]*AnisoRestLocal[0]+ShapeVectorResting[1][1]*AnisoRestLocal[1])/
//                             std::sqrt(ShapeVectorResting[1][0]*ShapeVectorResting[1][0]+ShapeVectorResting[1][1]*ShapeVectorResting[1][1]+0.0000001) );
      
//       teta[2] = std::acos(  (ShapeVectorResting[2][0]*AnisoRestLocal[0]+ShapeVectorResting[2][1]*AnisoRestLocal[1])/
//                             std::sqrt(ShapeVectorResting[2][0]*ShapeVectorResting[2][0]+ShapeVectorResting[2][1]*ShapeVectorResting[2][1]+0.0000001) );


      
//       // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> STRAIN and STRESS TENSOR (BEGIN) <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      
//       // deformation gradiant tensor F =Sigma i=1,2,3 Qi x Di
//       // strain tensor in resting shape E=0.5(FtF-I)
//       // trE
//       // B=FFt
//       // axa (direct product of aniso vector in resting shape)
//       // atEa
//       // E(axa) and (axa)E
//       double trE=( Delta[1]*cotan[0]+ Delta[2]*cotan[1]+Delta[0]*cotan[2])/(4*restingArea);
      
//       double directAniso[2][2]={{AnisoRestLocal[0]*AnisoRestLocal[0],AnisoRestLocal[0]*AnisoRestLocal[1]},
//                                 {AnisoRestLocal[1]*AnisoRestLocal[0],AnisoRestLocal[1]*AnisoRestLocal[1]}};
      

//       double LeftCauchy[2][2]; // B=FFt
//       LeftCauchy[0][0]=DeformGrad[0][0]*DeformGrad[0][0]+DeformGrad[0][1]*DeformGrad[0][1];
//       LeftCauchy[1][0]=DeformGrad[1][0]*DeformGrad[0][0]+DeformGrad[1][1]*DeformGrad[0][1];
//       LeftCauchy[0][1]=DeformGrad[0][0]*DeformGrad[1][0]+DeformGrad[0][1]*DeformGrad[1][1];
//       LeftCauchy[1][1]=DeformGrad[1][0]*DeformGrad[1][0]+DeformGrad[1][1]*DeformGrad[1][1];


//       double Egreen[2][2];//E=0.5(C-I)
//       Egreen[0][0]=0.5*(DeformGrad[0][0]*DeformGrad[0][0]+DeformGrad[1][0]*DeformGrad[1][0]-1);
//       Egreen[1][0]=0.5*(DeformGrad[0][1]*DeformGrad[0][0]+DeformGrad[1][1]*DeformGrad[1][0]);
//       Egreen[0][1]=0.5*(DeformGrad[0][0]*DeformGrad[0][1]+DeformGrad[1][0]*DeformGrad[1][1]);
//       Egreen[1][1]=0.5*(DeformGrad[0][1]*DeformGrad[0][1]+DeformGrad[1][1]*DeformGrad[1][1]-1);


//       double E2[2][2]; // used for energy calculation only
//       E2[0][0]=Egreen[0][0]*Egreen[0][0]+Egreen[0][1]*Egreen[1][0];
//       E2[1][0]=Egreen[1][0]*Egreen[0][0]+Egreen[1][1]*Egreen[1][0];
//       E2[0][1]=Egreen[0][0]*Egreen[0][1]+Egreen[0][1]*Egreen[1][1];
//       E2[1][1]=Egreen[1][0]*Egreen[0][1]+Egreen[1][1]*Egreen[1][1];
      
//       double I2=E2[0][0]+E2[1][1]; //trE2 used for energy calculation only
      
//        //atE2a used for energy calculation only
//       double I5= AnisoRestLocal[0]*AnisoRestLocal[0]*E2[0][0]   +  
//         AnisoRestLocal[0]*AnisoRestLocal[1]*(E2[0][1]+E2[1][0]) +
//         AnisoRestLocal[1]*AnisoRestLocal[1]*E2[1][1];
      
      
//       double StrainAlmansi[2][2]; // e=0.5(1-B^-1)  True strain tensor
//       temp=LeftCauchy[0][0]*LeftCauchy[1][1]-LeftCauchy[1][0]*LeftCauchy[0][1]; // det(B)
//       StrainAlmansi[0][0]=0.5*(1-(LeftCauchy[1][1]/temp));
//       StrainAlmansi[1][0]=0.5*LeftCauchy[1][0]/temp;
//       StrainAlmansi[0][1]=0.5*LeftCauchy[0][1]/temp;  
//       StrainAlmansi[1][1]=0.5*(1-(LeftCauchy[0][0]/temp));
      
//       double atEa=AnisoRestLocal[0]*AnisoRestLocal[0]*Egreen[0][0]
//                  +AnisoRestLocal[0]*AnisoRestLocal[1]*(Egreen[0][1]+Egreen[1][0])
//                  +AnisoRestLocal[1]*AnisoRestLocal[1]*Egreen[1][1];

//       double I4=atEa;      

//       double Eaa[2][2];
//       Eaa[0][0]= Egreen[0][0]*directAniso[0][0]+Egreen[0][1]*directAniso[1][0];        
//       Eaa[1][0]= Egreen[1][0]*directAniso[0][0]+Egreen[1][1]*directAniso[1][0];        
//       Eaa[0][1]= Egreen[0][0]*directAniso[0][1]+Egreen[0][1]*directAniso[1][1];        
//       Eaa[1][1]= Egreen[1][0]*directAniso[0][1]+Egreen[1][1]*directAniso[1][1];        

//       double aaE[2][2];
//       aaE[0][0]= directAniso[0][0]*Egreen[0][0]+directAniso[0][1]*Egreen[1][0];        
//       aaE[1][0]= directAniso[1][0]*Egreen[0][0]+directAniso[1][1]*Egreen[1][0];        
//       aaE[0][1]= directAniso[0][0]*Egreen[0][1]+directAniso[0][1]*Egreen[1][1];        
//       aaE[1][1]= directAniso[1][0]*Egreen[0][1]+directAniso[1][1]*Egreen[1][1];        
      
//       double B2[2][2];// LeftCauchy^2
//       B2[0][0]=LeftCauchy[0][0]*LeftCauchy[0][0]+LeftCauchy[0][1]*LeftCauchy[1][0];
//       B2[1][0]=LeftCauchy[1][0]*LeftCauchy[0][0]+LeftCauchy[1][1]*LeftCauchy[1][0];
//       B2[0][1]=LeftCauchy[0][0]*LeftCauchy[0][1]+LeftCauchy[0][1]*LeftCauchy[1][1];
//       B2[1][1]=LeftCauchy[1][0]*LeftCauchy[0][1]+LeftCauchy[1][1]*LeftCauchy[1][1];


//       double areaFactor=restingArea/Area; // 1/detF
//       //double areaFactor=restingArea/Area; // detF

//       double Sigma[2][2]; // true stress tensor (isotropic term) based on lambdaT and mioT
//       Sigma[0][0]=areaFactor*((lambdaT*trE-mioT)*LeftCauchy[0][0]+(mioT)*B2[0][0]);
//       Sigma[1][0]=areaFactor*((lambdaT*trE-mioT)*LeftCauchy[1][0]+(mioT)*B2[1][0]);
//       Sigma[0][1]=areaFactor*((lambdaT*trE-mioT)*LeftCauchy[0][1]+(mioT)*B2[0][1]);
//       Sigma[1][1]=areaFactor*((lambdaT*trE-mioT)*LeftCauchy[1][1]+(mioT)*B2[1][1]);
     
//       // double deltaS[2][2]; // based on Delin. paper
//       // deltaS[0][0]=deltaLam*(trE*directAniso[0][0]+atEa)+(2*deltaMio)*(Eaa[0][0]+aaE[0][0])-(deltaLam+2*deltaMio)*atEa*directAniso[0][0];
//       // deltaS[1][0]=deltaLam*(trE*directAniso[1][0]     )+(2*deltaMio)*(Eaa[1][0]+aaE[1][0])-(deltaLam+2*deltaMio)*atEa*directAniso[1][0];
//       // deltaS[0][1]=deltaLam*(trE*directAniso[0][1]     )+(2*deltaMio)*(Eaa[0][1]+aaE[0][1])-(deltaLam+2*deltaMio)*atEa*directAniso[0][1];
//       // deltaS[1][1]=deltaLam*(trE*directAniso[1][1]+atEa)+(2*deltaMio)*(Eaa[1][1]+aaE[1][1])-(deltaLam+2*deltaMio)*atEa*directAniso[1][1];

     
//       double deltaS[2][2]; // based on  equipartition energy 
//       deltaS[0][0]=(deltaLam/2)*(trE*directAniso[0][0]+atEa)+(deltaMio)*(Eaa[0][0]+aaE[0][0]);
//       deltaS[1][0]=(deltaLam/2)*(trE*directAniso[1][0]     )+(deltaMio)*(Eaa[1][0]+aaE[1][0]);
//       deltaS[0][1]=(deltaLam/2)*(trE*directAniso[0][1]     )+(deltaMio)*(Eaa[0][1]+aaE[0][1]);
//       deltaS[1][1]=(deltaLam/2)*(trE*directAniso[1][1]+atEa)+(deltaMio)*(Eaa[1][1]+aaE[1][1]);


//       // double deltaS[2][2]; // based on  ... 
//       // deltaS[0][0]=(deltaLam/2)*(trE*directAniso[0][0]+atEa)+2*(deltaMio)* atEa * directAniso[0][0];
//       // deltaS[1][0]=(deltaLam/2)*(trE*directAniso[1][0]     )+2*(deltaMio)* atEa * directAniso[1][0];
//       // deltaS[0][1]=(deltaLam/2)*(trE*directAniso[0][1]     )+2*(deltaMio)* atEa * directAniso[0][1];
//       // deltaS[1][1]=(deltaLam/2)*(trE*directAniso[1][1]+atEa)+2*(deltaMio)* atEa * directAniso[1][1];




//       // double deltaS[2][2]; // based on  ... 
//       // deltaS[0][0]=2*(deltaMio1)* atEa * directAniso[0][0];
//       // deltaS[1][0]=2*(deltaMio1)* atEa * directAniso[1][0];
//       // deltaS[0][1]=2*(deltaMio1)* atEa * directAniso[0][1];
//       // deltaS[1][1]=2*(deltaMio1)* atEa * directAniso[1][1];


        
//       strainZ =1-poissonT*((2*lambdaT*trE+2*mioT*trE)+deltaS[0][0]+deltaS[1][1])/youngT;

//       //<<<<<<<<<<<<<<<<<<<isotropic force from stress tensor <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

//       // double ss[2][2];//lambda(trE)I+2mioE
//       // ss[0][0]=lambdaT*trE+2*mioT*Egreen[0][0];
//       // ss[0][1]=            2*mioT*Egreen[0][1];
//       // ss[1][0]=            2*mioT*Egreen[1][0];
//       // ss[1][1]=lambdaT*trE+2*mioT*Egreen[1][1];

       

//       // double TPK[2][2];// 2nd Piola Kirchhoff stress tensor 
//       // TPK[0][0]=restingArea*(DeformGrad[0][0]*ss[0][0]+DeformGrad[0][1]*ss[1][0]);
//       // TPK[1][0]=restingArea*(DeformGrad[1][0]*ss[0][0]+DeformGrad[1][1]*ss[1][0]);
//       // TPK[0][1]=restingArea*(DeformGrad[0][0]*ss[0][1]+DeformGrad[0][1]*ss[1][1]);
//       // TPK[1][1]=restingArea*(DeformGrad[1][0]*ss[0][1]+DeformGrad[1][1]*ss[1][1]);

//       // //deltaFTPKlocal[i][0]= TPK[0][0]*ShapeVectorResting[i][0]+TPK[0][1]*ShapeVectorResting[i][1];
//       // //deltaFTPKlocal[i][1]= TPK[1][0]*ShapeVectorResting[i][0]+TPK[1][1]*ShapeVectorResting[i][1];
     
//       // double deltaFTPKlocal[2][2];
//       // deltaFTPKlocal[0][0]= TPK[0][0]*ShapeVectorResting[0][0]+TPK[0][1]*ShapeVectorResting[0][1];
//       // deltaFTPKlocal[0][1]= TPK[1][0]*ShapeVectorResting[0][0]+TPK[1][1]*ShapeVectorResting[0][1];
     
//       // double deltaFTPK[2][2]; 
//       // deltaFTPK[0][0]= rotation[0][0]*deltaFTPKlocal[0][0]+rotation[0][1]*deltaFTPKlocal[0][1];
//       // deltaFTPK[0][1]= rotation[1][0]*deltaFTPKlocal[0][0]+rotation[1][1]*deltaFTPKlocal[0][1];
//       // deltaFTPK[0][2]= rotation[2][0]*deltaFTPKlocal[0][0]+rotation[2][1]*deltaFTPKlocal[0][1];
//       //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

//       // //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//       double TPK[2][2];// 2nd Piola Kirchhoff stress tensor 
//       TPK[0][0]=restingArea*(DeformGrad[0][0]*deltaS[0][0]+DeformGrad[0][1]*deltaS[1][0]);
//       TPK[1][0]=restingArea*(DeformGrad[1][0]*deltaS[0][0]+DeformGrad[1][1]*deltaS[1][0]);
//       TPK[0][1]=restingArea*(DeformGrad[0][0]*deltaS[0][1]+DeformGrad[0][1]*deltaS[1][1]);
//       TPK[1][1]=restingArea*(DeformGrad[1][0]*deltaS[0][1]+DeformGrad[1][1]*deltaS[1][1]);

//       //deltaFTPKlocal[i][0]= TPK[0][0]*ShapeVectorResting[i][0]+TPK[0][1]*ShapeVectorResting[i][1];
//       //deltaFTPKlocal[i][1]= TPK[1][0]*ShapeVectorResting[i][0]+TPK[1][1]*ShapeVectorResting[i][1];
     
//       double deltaFTPKlocal[3][2];
//       deltaFTPKlocal[0][0]= TPK[0][0]*ShapeVectorResting[0][0]+TPK[0][1]*ShapeVectorResting[0][1];
//       deltaFTPKlocal[0][1]= TPK[1][0]*ShapeVectorResting[0][0]+TPK[1][1]*ShapeVectorResting[0][1];
     
//       deltaFTPKlocal[1][0]= TPK[0][0]*ShapeVectorResting[1][0]+TPK[0][1]*ShapeVectorResting[1][1];
//       deltaFTPKlocal[1][1]= TPK[1][0]*ShapeVectorResting[1][0]+TPK[1][1]*ShapeVectorResting[1][1];

//       deltaFTPKlocal[2][0]= TPK[0][0]*ShapeVectorResting[2][0]+TPK[0][1]*ShapeVectorResting[2][1];
//       deltaFTPKlocal[2][1]= TPK[1][0]*ShapeVectorResting[2][0]+TPK[1][1]*ShapeVectorResting[2][1];


//       double deltaFTPK[3][3]; 
//       deltaFTPK[0][0]= rotation[0][0]*deltaFTPKlocal[0][0]+rotation[0][1]*deltaFTPKlocal[0][1];
//       deltaFTPK[0][1]= rotation[1][0]*deltaFTPKlocal[0][0]+rotation[1][1]*deltaFTPKlocal[0][1];
//       deltaFTPK[0][2]= rotation[2][0]*deltaFTPKlocal[0][0]+rotation[2][1]*deltaFTPKlocal[0][1];

//       deltaFTPK[1][0]= rotation[0][0]*deltaFTPKlocal[1][0]+rotation[0][1]*deltaFTPKlocal[1][1];
//       deltaFTPK[1][1]= rotation[1][0]*deltaFTPKlocal[1][0]+rotation[1][1]*deltaFTPKlocal[1][1];
//       deltaFTPK[1][2]= rotation[2][0]*deltaFTPKlocal[1][0]+rotation[2][1]*deltaFTPKlocal[1][1];
      
//       deltaFTPK[2][0]= rotation[0][0]*deltaFTPKlocal[2][0]+rotation[0][1]*deltaFTPKlocal[2][1];
//       deltaFTPK[2][1]= rotation[1][0]*deltaFTPKlocal[2][0]+rotation[1][1]*deltaFTPKlocal[2][1];
//       deltaFTPK[2][2]= rotation[2][0]*deltaFTPKlocal[2][0]+rotation[2][1]*deltaFTPKlocal[2][1];


//       // //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


//       double deltaSFt[2][2];
//       deltaSFt[0][0]=deltaS[0][0]*DeformGrad[0][0]+deltaS[0][1]*DeformGrad[0][1];
//       deltaSFt[1][0]=deltaS[1][0]*DeformGrad[0][0]+deltaS[1][1]*DeformGrad[0][1];
//       deltaSFt[0][1]=deltaS[0][0]*DeformGrad[1][0]+deltaS[0][1]*DeformGrad[1][1];
//       deltaSFt[1][1]=deltaS[1][0]*DeformGrad[1][0]+deltaS[1][1]*DeformGrad[1][1];
      
//       double deltaSigma[2][2];// true stress tensor (anisotropic correction term)deltaLambda and deltaMio (Longitudinal-Transverse)
//       deltaSigma[0][0]=areaFactor*(DeformGrad[0][0]*deltaSFt[0][0]+DeformGrad[0][1]*deltaSFt[1][0]);
//       deltaSigma[1][0]=areaFactor*(DeformGrad[1][0]*deltaSFt[0][0]+DeformGrad[1][1]*deltaSFt[1][0]);
//       deltaSigma[0][1]=areaFactor*(DeformGrad[0][0]*deltaSFt[0][1]+DeformGrad[0][1]*deltaSFt[1][1]);
//       deltaSigma[1][1]=areaFactor*(DeformGrad[1][0]*deltaSFt[0][1]+DeformGrad[1][1]*deltaSFt[1][1]);
    

  
     
//       double StressTensor[3][3];
//       StressTensor[0][0]=Sigma[0][0]+deltaSigma[0][0];
//       StressTensor[1][0]=Sigma[1][0]+deltaSigma[1][0];
//       StressTensor[0][1]=Sigma[0][1]+deltaSigma[0][1];
//       StressTensor[1][1]=Sigma[1][1]+deltaSigma[1][1];
     
 
//       //Shape vectors in Current shape (counterclockwise ordering of nodes/edges)     ShapeVectorCurrent[3][3]  calculated above   
//       //.............................. ( or clockwise ordering of nodes/edges)
          
 
//       //square of radius of circumstancing circle in resting shape
//       //double Rcirc2Resting=(0.25*restingLength[0]*restingLength[1]*restingLength[2]/Area)*(0.25*restingLength[0]*restingLength[1]*restingLength[2]/Area);  
      
//       double StrainTensor[3][3];
//       //1
//       StrainTensor[0][0]=StrainAlmansi[0][0];
//       StrainTensor[1][0]=StrainAlmansi[1][0];
//       StrainTensor[0][1]=StrainAlmansi[0][1];
//       StrainTensor[1][1]=StrainAlmansi[1][1];


//       StrainTensor[0][2]=0;  // adding 3rd dimension which is zero, the tensor is still in element plane
//       StrainTensor[1][2]=0;
//       StrainTensor[2][2]=0;
//       StrainTensor[2][0]=0;
//       StrainTensor[2][1]=0;
      
//       StressTensor[0][2]=0;  // adding 3rd dimension which is zero, the tensor is still in element plane
//       StressTensor[1][2]=0;
//       StressTensor[2][2]=0;
//       StressTensor[2][0]=0;
//       StressTensor[2][1]=0;
 

      

//       //rotation matrix to go to global coordinate system based on counterclockwise ordering;   rotation[3][3] calculated above  
//       //double testRatio1=Sigma[0][0]*
//       // double testRatio1=StressTensor[0][0]*(1-2*StrainTensor[0][0])*(1-2*StrainTensor[0][0])/((mioT*StrainTensor[0][0]+lambdaT*trE*(1-2*StrainTensor[0][0]))*(Area/restingArea));
//       // std::cerr<< "test1 "<<testRatio1<<std::endl;
//       // double testRatio2=StressTensor[1][1]*(1-2*StrainTensor[1][1])*(1-2*StrainTensor[1][1])/((mioT*StrainTensor[1][1]+lambdaT*trE*(1-2*StrainTensor[1][1]))*(Area/restingArea));
//       // std::cerr<< "test2 "<<testRatio2<<std::endl;
//       // std::cerr<< "trE1 "<<trE<<std::endl;
//       // std::cerr<< "trE2 "<<Egreen[0][0]+Egreen[1][1]<<std::endl;
//       // std::cerr<< "mioT "<<mioT<<std::endl;
//       // std::cerr<< "Area "<<Area<<std::endl;
//       //std::cerr<< "restingArea "<<restingArea<<std::endl;

//       //   std::cerr<< "alfa0           "<<StrainAlmansi[0][0]<<std::endl;

//       // std::cerr<< "0.5(1-(1/beta0)) "<<0.5*(1-(1/LeftCauchy[0][0]))<<std::endl;
//       // std::cerr<< "0.5(1-(1/beta1)) "<<0.5*(1-(1/LeftCauchy[1][1]))<<std::endl;
//       // std::cerr<< "almansi 01           "<<StrainAlmansi[0][1]<<std::endl;
//       // std::cerr<< "almansi 10           "<<StrainAlmansi[1][0]<<std::endl;
      
//       // std::cerr <<"rotation tensor " << std::endl;
//       // std::cerr <<" Rxx  "<< rotation[0][0] <<" Rxy  "<< rotation[0][1] <<" Rxz  "<< rotation[0][2] << std::endl
//       //           <<" Ryx  "<< rotation[1][0] <<" Ryy  "<< rotation[1][1] <<" Ryz  "<< rotation[1][2] << std::endl
//       //           <<" Rzx  "<< rotation[2][0] <<" Rzy  "<< rotation[2][1] <<" Rzz  "<< rotation[2][2] << std::endl <<std::endl;
//       //   std::cerr<< "deltaSigma[0][0] and [1][1]  "<<deltaSigma[0][0]<<"  "<<deltaSigma[1][1]<<std::endl;
//       // std::cerr<< "Sigma[0][0] and [1][1]  "<<Sigma[0][0]<<"  "<<Sigma[1][1]<<std::endl;
      
//       // double test1=((mioT*StrainTensor[0][0]+lambdaT*trE*(1-2*StrainTensor[0][0]))*(Area/restingArea))/((1-2*StrainTensor[0][0])*(1-2*StrainTensor[0][0]));
//       // double test2=((mioT*StrainTensor[1][1]+lambdaT*trE*(1-2*StrainTensor[1][1]))*(Area/restingArea))/((1-2*StrainTensor[1][1])*(1-2*StrainTensor[1][1]));
//       // std::cerr<< "test1  test2  "<<test1<<" "<<test2<<std::endl;
      
 
//          // std::cerr<< "Egreen[0][0] and Egreen[1][1]  "<<Egreen[0][0]<<"  "<<Egreen[1][1]<<std::endl;
//          // std::cerr<< "cauchyStress[0][0] and cauchyStressEgreen[1][1]  "<<lambdaT*trE+mioT*Egreen[0][0]<<"  "<<lambdaT*trE+mioT*Egreen[1][1]<<std::endl;
      
 
//       // rotating strain tensor to the global coordinate system
//       double tempR[3][3]={{0,0,0},{0,0,0},{0,0,0}};
//       for (int r=0 ; r<3 ; r++) {
//         for (int s=0 ; s<3 ; s++) {
//           for(int w=0 ; w<3 ; w++) {
//             tempR[r][s]=tempR[r][s]+rotation[r][w]*StrainTensor[w][s];
//           }
//         }
//       }
//       for (int r=0 ; r<3 ; r++) {
//         for (int s=0 ; s<3 ; s++) {
//           StrainTensor[r][s]=0;
//           for(int w=0 ; w<3 ; w++) {
//             StrainTensor[r][s]=StrainTensor[r][s]+tempR[r][w]*rotation[s][w];
//           }
//         }
//       }
      
  
//       // rotating stress tensor to the global coordinate system
  
//       for (int r=0 ; r<3 ; r++) {
//         for (int s=0 ; s<3 ; s++) {
//           tempR[r][s]=0;
//         }
//       }
 
//       for (int r=0 ; r<3 ; r++) {
//         for (int s=0 ; s<3 ; s++) {
//           for(int w=0 ; w<3 ; w++) {
//             tempR[r][s]=tempR[r][s]+rotation[r][w]*StressTensor[w][s];
//           }
//         }
//       }
//       for (int r=0 ; r<3 ; r++) {
//         for (int s=0 ; s<3 ; s++) {
//           StressTensor[r][s]=0;
//           for(int w=0 ; w<3 ; w++) {
//             StressTensor[r][s]=StressTensor[r][s]+tempR[r][w]*rotation[s][w];
//           }
//         }
//       }

      

//       //Cell *  cell1=&(T.cell(cellIndex));
//       // std::cerr<<" cell   "<<cellIndex<<"   wall   "<<w1<<  "  cell neighbohr  "<<(cell1->cellNeighbor(0))->index() <<std::endl;
//       // std::cerr<<" cell   "<<cellIndex<<"   wall   "<<w2<<  "  cell neighbohr  "<<(cell1->cellNeighbor(1))->index() <<std::endl;
//       // std::cerr<<" cell   "<<cellIndex<<"   wall   "<<w3<<  "  cell neighbohr  "<<(cell1->cellNeighbor(2))->index() <<std::endl;
//       cellData[cellIndex][26]=StressTensor[0][0];
//       cellData[cellIndex][27]=StressTensor[1][1];
//       cellData[cellIndex][28]=StressTensor[2][2];
//       cellData[cellIndex][29]=StressTensor[0][1];
//       cellData[cellIndex][30]=StressTensor[0][2];
//       cellData[cellIndex][31]=StressTensor[1][2];
      
      
      
      
//       // stress component along MT direction 
//       temp = std::sqrt ( cellData[cellIndex][variableIndex(0,1)  ] * cellData[cellIndex][variableIndex(0,1)  ] +
// 			 cellData[cellIndex][variableIndex(0,1)+1] * cellData[cellIndex][variableIndex(0,1)+1] +
// 			 cellData[cellIndex][variableIndex(0,1)+2] * cellData[cellIndex][variableIndex(0,1)+2]  );
//       cellData[cellIndex][variableIndex(0,1)  ] /=temp;
//       cellData[cellIndex][variableIndex(0,1)+1] /=temp;
//       cellData[cellIndex][variableIndex(0,1)+2] /=temp;
      
//       cellData[cellIndex][variableIndex(0,8)  ] =
// 	cellData[cellIndex][variableIndex(0,1)     ] *cellData[cellIndex][variableIndex(0,1)     ] *StressTensor[0][0]  +
// 	cellData[cellIndex][variableIndex(0,1)     ] *cellData[cellIndex][variableIndex(0,1)+1] *StressTensor[0][1]  +
// 	cellData[cellIndex][variableIndex(0,1)     ] *cellData[cellIndex][variableIndex(0,1)+2] *StressTensor[0][2]  +
// 	cellData[cellIndex][variableIndex(0,1)+1] *cellData[cellIndex][variableIndex(0,1)     ] *StressTensor[1][0]  +
// 	cellData[cellIndex][variableIndex(0,1)+1] *cellData[cellIndex][variableIndex(0,1)+1] *StressTensor[1][1]  +
// 	cellData[cellIndex][variableIndex(0,1)+1] *cellData[cellIndex][variableIndex(0,1)+2] *StressTensor[1][2]  +
// 	cellData[cellIndex][variableIndex(0,1)+2] *cellData[cellIndex][variableIndex(0,1)     ] *StressTensor[2][0]  +
// 	cellData[cellIndex][variableIndex(0,1)+2] *cellData[cellIndex][variableIndex(0,1)+1] *StressTensor[2][1]  +
// 	cellData[cellIndex][variableIndex(0,1)+2] *cellData[cellIndex][variableIndex(0,1)+2] *StressTensor[2][2]       ;
  
    
//       // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> STRAIN and STRESS TENSORS (END) <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      
//       // eigenvalues/eigenvectors of strain tensor in global coordinate system. (Jacobi method)
      
//       double strainEpcilon =0.000001;
      

//       double eigenVectorStrain[3][3]={{1,0,0},{0,1,0},{0,0,1}};
//       double pivot=1;
//       double pi=3.14159265;
//       int I,J;
//       double RotAngle,Si,Co;
      
//       pivot=std::abs(StrainTensor[1][0]);
//       I=1;
//       J=0;
//       if (std::abs(StrainTensor[2][0])>pivot) {
//         pivot=std::abs(StrainTensor[2][0]);
//         I=2;
//         J=0;
//       }
//       if (std::abs(StrainTensor[2][1])>pivot) {
//         pivot=std::abs(StrainTensor[2][1]);
//         I=2;
//         J=1;
//       }

//       while (pivot>strainEpcilon) {
      
//         if (std::abs(StrainTensor[I][I]-StrainTensor[J][J])<strainEpcilon) {
//           RotAngle=pi/4;
//         }            
//         else {
//           RotAngle=0.5*std::atan((2*StrainTensor[I][J])/(StrainTensor[J][J]-StrainTensor[I][I]));
//         }
//         Si=std::sin(RotAngle);
//         Co=std::cos(RotAngle);
//         double tempRot[3][3]={{1,0,0},{0,1,0},{0,0,1}};
//         tempRot[I][I]=Co;
//         tempRot[J][J]=Co;
//         tempRot[I][J]=Si;
//         tempRot[J][I]=-Si;
//         double tempStrain[3][3]={{0,0,0},{0,0,0},{0,0,0}};
//         for (int r=0 ; r<3 ; r++) {
//           for (int s=0 ; s<3 ; s++) {
//             for(int w=0 ; w<3 ; w++) {
//               tempStrain[r][s]=tempStrain[r][s]+StrainTensor[r][w]*tempRot[w][s];
//             }
//           }
//         }
        
//         for (int r=0 ; r<3 ; r++) {
//           for (int s=0 ; s<3 ; s++) {
//             StrainTensor[r][s]=0;
//           }
//         }
        
//         for (int r=0 ; r<3 ; r++) {
//           for (int s=0 ; s<3 ; s++) {
//             for(int w=0 ; w<3 ; w++) {
//               StrainTensor[r][s]=StrainTensor[r][s]+tempRot[w][r]*tempStrain[w][s];
//             }
//           }
//         }
        
//         for (int r=0 ; r<3 ; r++) 
//           {
//             for (int s=0 ; s<3 ; s++) 
//               {
//                 tempStrain[r][s]=eigenVectorStrain[r][s];
//               }
//           }

// 	for (size_t ii=0; ii<3; ++ii) {
// 	  for (size_t jj=0; jj<3; ++jj) {
// 	    eigenVectorStrain[ii][jj] = 0.0; 
// 	  }
// 	}	
//         for (int r=0 ; r<3 ; r++) {
//           for (int s=0 ; s<3 ; s++) {
//             for(int w=0 ; w<3 ; w++) {
//               eigenVectorStrain[r][s]=eigenVectorStrain[r][s]+tempStrain[r][w]*tempRot[w][s];
//             }
//           }
//         }
//         pivot=std::abs(StrainTensor[1][0]);
//         I=1;
//         J=0;
//         if (std::abs(StrainTensor[2][0])>pivot) {
//           pivot=std::abs(StrainTensor[2][0]);
//           I=2;
//           J=0;
//         }
//         if (std::abs(StrainTensor[2][1])>pivot) {
//           pivot=std::abs(StrainTensor[2][1]);
//           I=2;
//           J=1;
//         }
//       }
      
      
//       // maximal strain direction
//       double maximalStrainValue=StrainTensor[0][0];
//       int Istrain=0;
//       if (std::abs(StrainTensor[1][1])>std::abs(maximalStrainValue)) 
//         {
//           maximalStrainValue=StrainTensor[1][1];
//           Istrain=1;
//         }
//       if (std::abs(StrainTensor[2][2])>std::abs(maximalStrainValue)) 
//         {
//           maximalStrainValue=StrainTensor[2][2];
//           Istrain=2;
//         }
//       //std::cerr<<"maximal Strain direction "<< eigenVectorStrain[0][Istrain] <<" "<< eigenVectorStrain[1][Istrain] <<" "<< eigenVectorStrain[2][Istrain] <<std::endl;  
//       //std::cerr<<"maximal Strain value "<< maximalStrainValue <<std::endl;  
      
//       // 2nd maximalstrain direction/value
//       double maximalStrainValue2;
//       int Istrain2,Istrain3;
//       if (Istrain==0) {
//         Istrain2=1;
//         Istrain3=2;
//       }
//       if (Istrain==1) {
//         Istrain2=0;
//         Istrain3=2;
//       }
//       if (Istrain==2) {
//         Istrain2=0;
//         Istrain3=1;
//       }
//       if(std::abs(StrainTensor[Istrain3][Istrain3])>std::abs(StrainTensor[Istrain2][Istrain2])) {
//         Istrain2=Istrain3;
//       }
//       maximalStrainValue2=StrainTensor[Istrain2][Istrain2];
      
      
      

   
//       //perpendicular direction to strain in cellData
      
//       // normal to the cell plane in global direction is Zcurrent[], vector product gives the perpendicular strain direction
//       double PerpStrain[3];
//       PerpStrain[0]=Zcurrent[1]*eigenVectorStrain[2][Istrain]-Zcurrent[2]*eigenVectorStrain[1][Istrain];
//       PerpStrain[1]=Zcurrent[2]*eigenVectorStrain[0][Istrain]-Zcurrent[0]*eigenVectorStrain[2][Istrain];
//       PerpStrain[2]=Zcurrent[0]*eigenVectorStrain[1][Istrain]-Zcurrent[1]*eigenVectorStrain[0][Istrain];
//       temp=std::sqrt(PerpStrain[0]*PerpStrain[0]+PerpStrain[1]*PerpStrain[1]+PerpStrain[2]*PerpStrain[2]);     
      
//        if(std::abs(temp)<0.0001){ // if maximal strain is normal to the cell plane storing strain direction instead as it should not be used
//          PerpStrain[0]=eigenVectorStrain[0][Istrain];
//          PerpStrain[1]=eigenVectorStrain[1][Istrain];
//          PerpStrain[2]=eigenVectorStrain[2][Istrain];
//        }
       
//        //   if (dimension==2)
//        //     {
//        //       cellData[cellIndex][0]=PerpStrain[0];        
//        //       cellData[cellIndex][1]=PerpStrain[1];
//        //     }
//        //   if (dimension==3)
//        //     {
//        //       cellData[cellIndex][0]=PerpStrain[0];
//        //       cellData[cellIndex][1]=PerpStrain[1];
//        //       cellData[cellIndex][2]=PerpStrain[2];
//        //       // cellData[cellIndex][3]=10*maximalStrainValue;  //NOTE maximal Strain and Stress Values should be used somewhere this is an option
//        //     }
       

//        // storing a measure for strain anisotropy in cell vector
//        if (std::abs(maximalStrainValue)<0.0000001) cellData[cellIndex][variableIndex(0,2)]=0;
//        if (std::abs(maximalStrainValue)>= 0.0000001) cellData[cellIndex][variableIndex(0,2)]=1-std::abs(maximalStrainValue2/maximalStrainValue);
       
//        // storing   strain direction/value in cellData
//        if (numVariableIndexLevel()==3 && (numVariableIndex(1)==1 || numVariableIndex(1)==2 || numVariableIndex(1)==3) ) {// storing maximal strain
//          if (dimension==2)
//            {
//              cellData[cellIndex][variableIndex(1,0)]  =eigenVectorStrain[0][Istrain];
//              cellData[cellIndex][variableIndex(1,0)+1]=eigenVectorStrain[1][Istrain];
//              cellData[cellIndex][variableIndex(1,0)+3]=maximalStrainValue;  //maximal Strain Value is stored after its eigenvector
//            }
//          if (dimension==3)
//            {
//              cellData[cellIndex][variableIndex(1,0)]  =eigenVectorStrain[0][Istrain];
//              cellData[cellIndex][variableIndex(1,0)+1]=eigenVectorStrain[1][Istrain];
//              cellData[cellIndex][variableIndex(1,0)+2]=eigenVectorStrain[2][Istrain];
//              cellData[cellIndex][variableIndex(1,0)+3]=maximalStrainValue;  //maximal Strain Value is stored after its eigenvector
//            }
//        }
       
//        if (numVariableIndexLevel()==3 &&  numVariableIndex(1)==3  ) {//storing 2nd maximal strain
//          if (dimension==2)
//            {
//              cellData[cellIndex][variableIndex(1,2)]  =eigenVectorStrain[0][Istrain2];
//              cellData[cellIndex][variableIndex(1,2)+1]=eigenVectorStrain[1][Istrain2];
//              cellData[cellIndex][variableIndex(1,2)+3]=maximalStrainValue2;  //2nd maximal Strain Value is stored after its eigenvector
//            }
//          if (dimension==3)
//            {
//              cellData[cellIndex][variableIndex(1,2)]  =eigenVectorStrain[0][Istrain2];
//              cellData[cellIndex][variableIndex(1,2)+1]=eigenVectorStrain[1][Istrain2];
//              cellData[cellIndex][variableIndex(1,2)+2]=eigenVectorStrain[2][Istrain2];
//              cellData[cellIndex][variableIndex(1,2)+3]=maximalStrainValue2;  //2nd maximal Strain Value is stored after its eigenvector
//            }
//        }
//        if (numVariableIndexLevel()==3 && ( numVariableIndex(1)==2 || numVariableIndex(1)==3 ) ) {//storing perpendicular to maximal strain
//          if (dimension==2)
//            {
//              cellData[cellIndex][variableIndex(1,1)]  =PerpStrain[0];
//              cellData[cellIndex][variableIndex(1,1)+1]=PerpStrain[1];
//              cellData[cellIndex][variableIndex(1,1)+3]=maximalStrainValue;//maximal Strain Value is stored after its eigenvector
//              //cellData[cellIndex][variableIndex(1,1)+3]=Area/restingArea;
//            }
//          if (dimension==3)
//            {
//              cellData[cellIndex][variableIndex(1,1)]  =PerpStrain[0];
//              cellData[cellIndex][variableIndex(1,1)+1]=PerpStrain[1];
//              cellData[cellIndex][variableIndex(1,1)+2]=PerpStrain[2];
//              cellData[cellIndex][variableIndex(1,1)+3]=maximalStrainValue;//maximal Strain Value is stored after its eigenvector
//              //cellData[cellIndex][variableIndex(1,1)+3]=Area/restingArea;
//            }
//        }


       

       
//        //<<<<<<<<<<<<<<<<<<<<<<<< angles between  vectors and circumferential direction <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      
//        temp=std::sqrt(cellData[cellIndex][variableIndex(0,1)]*cellData[cellIndex][variableIndex(0,1)]+     // 12 --> MT
//                       cellData[cellIndex][variableIndex(0,1)+1]*cellData[cellIndex][variableIndex(0,1)+1]);
//        if(temp<0.000001){
//          cellData[cellIndex][12]=pi/2;
//        }
//        else{
//          cellData[cellIndex][12]=std::atan(cellData[cellIndex][variableIndex(0,1)+2]/temp);
//        }
       
//        temp=std::sqrt(eigenVectorStrain[0][Istrain]*eigenVectorStrain[0][Istrain]+                                   // 13 --> Strain
//                       eigenVectorStrain[1][Istrain]*eigenVectorStrain[1][Istrain] );
//        if(temp<0.000000001){
//          cellData[cellIndex][13]=pi/2;
//        }
//        else{
//          cellData[cellIndex][13]=std::atan(eigenVectorStrain[2][Istrain]/temp);
//        }

      
//        cellData[cellIndex][15]=position[0][2]; // z coordinate of 0th vertex of the cell                                  // 15 --> Z coordinate
      
       
       
//        //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
       
       
//        //---- Anisotropic Correction Force-------------------------------
//        double deltaF[3][3];
       
        
//          for ( int i=0 ; i<3 ; ++i )  // from stress tensor(equipartitioning energy)
//           for ( int j=0 ; j<3 ; ++j )
//             deltaF[i][j]=(-deltaFTPK[i][j]);
      
//         double  I1=trE;

//         // energy
//         EnergyIso =( (lambdaT/2)*I1*I1 + mioT*I2 );//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//         EnergyAniso =( (deltaLam/2)*I4*I1 + deltaMio*I5 ); //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      
//         //Forces of vertices   
//         double Force[3][3];                                           
    
//         Force[0][0]= (tensileStiffness[0]*Delta[0]+angularStiffness[1]*Delta[1]+angularStiffness[0]*Delta[2])*(position[1][0]-position[0][0])
//           +(tensileStiffness[2]*Delta[2]+angularStiffness[2]*Delta[1]+angularStiffness[0]*Delta[0])*(position[2][0]-position[0][0])
//           + deltaF[0][0]; 
//         Force[0][1]= (tensileStiffness[0]*Delta[0]+angularStiffness[1]*Delta[1]+angularStiffness[0]*Delta[2])*(position[1][1]-position[0][1])
//           +(tensileStiffness[2]*Delta[2]+angularStiffness[2]*Delta[1]+angularStiffness[0]*Delta[0])*(position[2][1]-position[0][1])
//           + deltaF[0][1];  
//         Force[0][2]= (tensileStiffness[0]*Delta[0]+angularStiffness[1]*Delta[1]+angularStiffness[0]*Delta[2])*(position[1][2]-position[0][2])
//           +(tensileStiffness[2]*Delta[2]+angularStiffness[2]*Delta[1]+angularStiffness[0]*Delta[0])*(position[2][2]-position[0][2])
//           + deltaF[0][2]; 
//         Force[1][0]= (tensileStiffness[0]*Delta[0]+angularStiffness[0]*Delta[2]+angularStiffness[1]*Delta[1])*(position[0][0]-position[1][0])
//           +(tensileStiffness[1]*Delta[1]+angularStiffness[2]*Delta[2]+angularStiffness[1]*Delta[0])*(position[2][0]-position[1][0])
//           + deltaF[1][0];  
//         Force[1][1]= (tensileStiffness[0]*Delta[0]+angularStiffness[0]*Delta[2]+angularStiffness[1]*Delta[1])*(position[0][1]-position[1][1])
//           +(tensileStiffness[1]*Delta[1]+angularStiffness[2]*Delta[2]+angularStiffness[1]*Delta[0])*(position[2][1]-position[1][1])
//           + deltaF[1][1];  
//         Force[1][2]= (tensileStiffness[0]*Delta[0]+angularStiffness[0]*Delta[2]+angularStiffness[1]*Delta[1])*(position[0][2]-position[1][2])
//           +(tensileStiffness[1]*Delta[1]+angularStiffness[2]*Delta[2]+angularStiffness[1]*Delta[0])*(position[2][2]-position[1][2])
//           + deltaF[1][2];  
        
//         Force[2][0]= (tensileStiffness[2]*Delta[2]+angularStiffness[0]*Delta[0]+angularStiffness[2]*Delta[1])*(position[0][0]-position[2][0])
//           +(tensileStiffness[1]*Delta[1]+angularStiffness[1]*Delta[0]+angularStiffness[2]*Delta[2])*(position[1][0]-position[2][0])
//           + deltaF[2][0];  
//         Force[2][1]= (tensileStiffness[2]*Delta[2]+angularStiffness[0]*Delta[0]+angularStiffness[2]*Delta[1])*(position[0][1]-position[2][1])
//           +(tensileStiffness[1]*Delta[1]+angularStiffness[1]*Delta[0]+angularStiffness[2]*Delta[2])*(position[1][1]-position[2][1])
//           + deltaF[2][1];  
//         Force[2][2]= (tensileStiffness[2]*Delta[2]+angularStiffness[0]*Delta[0]+angularStiffness[2]*Delta[1])*(position[0][2]-position[2][2])
//           +(tensileStiffness[1]*Delta[1]+angularStiffness[1]*Delta[0]+angularStiffness[2]*Delta[2])*(position[1][2]-position[2][2])
//           + deltaF[2][2];  
        
      
        
//         // adding TRBSMT forces to the total vertexDerivs
        
//         vertexDerivs[v1][0]+= Force[0][0];
//         vertexDerivs[v1][1]+= Force[0][1];
//         vertexDerivs[v1][2]+= Force[0][2];
        
//         vertexDerivs[v2][0]+= Force[1][0];
//         vertexDerivs[v2][1]+= Force[1][1];
//         vertexDerivs[v2][2]+= Force[1][2];
        
//         vertexDerivs[v3][0]+= Force[2][0];
//         vertexDerivs[v3][1]+= Force[2][1];
//         vertexDerivs[v3][2]+= Force[2][2];

 


//         //cellData[cellIndex][19]=TetaPerpStress;
//         //cellData[cellIndex][20]= strainZ;
//         //cellData[cellIndex][21]= areaRatio;
//         // cellData[cellIndex][21]= TetaStress;
//         // cellData[cellIndex][22]= TetaStrain; 
//         // cellData[cellIndex][23]= TetaPerp;
//         // cellData[cellIndex][19]= TETA;
//         // cellData[cellIndex][20]=EnergyIso+EnergyAniso;    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//         cellData[cellIndex][variableIndex(0,4)]= Area/restingArea; //area ratio
//         cellData[cellIndex][variableIndex(0,5)]=EnergyIso;    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//         cellData[cellIndex][variableIndex(0,6)]=EnergyAniso;  
        
//   }
   

//   double totalEnergyIso=0;
//   double totalEnergyAniso=0;

//   for( size_t cellIndex=0 ; cellIndex<numCells ; ++cellIndex ) {
   
//     totalEnergyIso +=cellData[cellIndex][variableIndex(0,5)];
//     totalEnergyAniso +=cellData[cellIndex][variableIndex(0,6)];

//     Cell *  cell1=&(T.cell(cellIndex));
 	
//     int neighbor[3];  
//     neighbor[0]=(cell1->cellNeighbor(0))->index();
//     neighbor[1]=(cell1->cellNeighbor(1))->index();
//     neighbor[2]=(cell1->cellNeighbor(2))->index();

//     double neighborweight=parameter(5);

//     double StressTensor[3][3]={{0,0,0},{0,0,0},{0,0,0}};
    
//     int counter=0;
//     for (int nn=0 ; nn<3 ; nn++){
//       if (neighbor[nn]<numCells && neighbor[nn]>-1){
// 	StressTensor[0][0]+=neighborweight*cellData[neighbor[nn]][26];
// 	counter+=1;
//       }
//     }
//     StressTensor[0][0]/=counter;
//     StressTensor[0][0]+=(1-neighborweight)*cellData[cellIndex][26];

//     counter=0;
//     for (int nn=0 ; nn<3 ; nn++){
//       if (neighbor[nn]<numCells && neighbor[nn]>-1){
// 	StressTensor[1][1]+=neighborweight*cellData[neighbor[nn]][27];
// 	counter+=1;
//       }
//     }
//     StressTensor[1][1]/=counter;
//     StressTensor[1][1]+=(1-neighborweight)*cellData[cellIndex][27];

//     counter=0;
//     for (int nn=0 ; nn<3 ; nn++){
//       if (neighbor[nn]<numCells && neighbor[nn]>-1){
// 	StressTensor[2][2]+=neighborweight*cellData[neighbor[nn]][28];
// 	counter+=1;
//       }
//     }
//     StressTensor[2][2]/=counter;
//     StressTensor[2][2]+=(1-neighborweight)*cellData[cellIndex][28];

//     counter=0;
//     for (int nn=0 ; nn<3 ; nn++){
//       if (neighbor[nn]<numCells && neighbor[nn]>-1){
// 	StressTensor[0][1]+=neighborweight*cellData[neighbor[nn]][29];
// 	counter+=1;
//       }
//     }
//     StressTensor[0][1]/=counter;
//     StressTensor[0][1]+=(1-neighborweight)*cellData[cellIndex][29];

//     counter=0;
//     for (int nn=0 ; nn<3 ; nn++){
//       if (neighbor[nn]<numCells && neighbor[nn]>-1){
// 	StressTensor[2][0]+=neighborweight*cellData[neighbor[nn]][30];
// 	counter+=1;
//       }
//     }
//     StressTensor[2][0]/=counter;
//     StressTensor[2][0]+=(1-neighborweight)*cellData[cellIndex][30];

//     counter=0;
//     for (int nn=0 ; nn<3 ; nn++){
//       if (neighbor[nn]<numCells && neighbor[nn]>-1){
// 	StressTensor[1][2]+=neighborweight*cellData[neighbor[nn]][31];
// 	counter+=1;
//       }
//     }
//     StressTensor[1][2]/=counter;
//     StressTensor[1][2]+=(1-neighborweight)*cellData[cellIndex][31];
    
//     StressTensor[0][2]=StressTensor[2][0];
//     StressTensor[1][0]=StressTensor[0][1];
//     StressTensor[2][1]=StressTensor[1][2];
       
//     // eigenvalue/eigenvectors of  stress tensor in global coordinate system. (Jacobi method)
      
//     double pi=3.14159265;
//     double stressEpcilon=0.000001;    

//     double RotAngle,Si,Co;
//     double eigenVectorStress[3][3]={{1,0,0},{0,1,0},{0,0,1}};
//     int I,J;    
//     double pivot=1;
    
//     while (pivot>stressEpcilon) {
//       pivot=std::abs(StressTensor[1][0]);
//       I=1;
//       J=0;
//       if (std::abs(StressTensor[2][0])>pivot) {
// 	pivot=std::abs(StressTensor[2][0]);
// 	I=2;
// 	J=0;
//       }
//       if (std::abs(StressTensor[2][1])>pivot) {
// 	pivot=std::abs(StressTensor[2][1]);
// 	I=2;
// 	J=1;
//       }
//       if (std::abs(StressTensor[I][I]-StressTensor[J][J])<stressEpcilon) {
// 	RotAngle=pi/4;
//       }            
//       else {
// 	RotAngle=0.5*std::atan((2*StressTensor[I][J])/(StressTensor[J][J]-StressTensor[I][I]));
//       }
//       Si=std::sin(RotAngle);
//       Co=std::cos(RotAngle);
//       double tempRot[3][3]={{1,0,0},{0,1,0},{0,0,1}};
//       tempRot[I][I]=Co;
//       tempRot[J][J]=Co;
//       tempRot[I][J]=Si;
//       tempRot[J][I]=-Si;
//       double tempStress[3][3]={{0,0,0},{0,0,0},{0,0,0}};
//       for (int r=0 ; r<3 ; r++) 
// 	for (int s=0 ; s<3 ; s++) 
// 	  for(int w=0 ; w<3 ; w++) 
// 	    tempStress[r][s]=tempStress[r][s]+StressTensor[r][w]*tempRot[w][s];
      
//       for (int r=0 ; r<3 ; r++) 
// 	for (int s=0 ; s<3 ; s++) 
// 	  StressTensor[r][s]=0;
      
//       for (int r=0 ; r<3 ; r++) 
// 	for (int s=0 ; s<3 ; s++) 
// 	  for(int w=0 ; w<3 ; w++) 
// 	    StressTensor[r][s]=StressTensor[r][s]+tempRot[w][r]*tempStress[w][s];
      
//       for (int r=0 ; r<3 ; r++) 
// 	for (int s=0 ; s<3 ; s++) 
// 	  tempStress[r][s]=eigenVectorStress[r][s];
      
//       for (size_t ii=0; ii<3; ++ii) 
// 	for (size_t jj=0; jj<3; ++jj) 
// 	  eigenVectorStress[ii][jj] = 0.0;
      
//       for (int r=0 ; r<3 ; r++) 
// 	for (int s=0 ; s<3 ; s++) 
// 	  for(int w=0 ; w<3 ; w++) 
// 	    eigenVectorStress[r][s]=eigenVectorStress[r][s]+tempStress[r][w]*tempRot[w][s];
      
//     }
    
    
//     // maximal stress direction
//     double maximalStressValue=StressTensor[0][0];
//     int Istress=0;
//     if (std::abs(StressTensor[1][1])>std::abs(maximalStressValue)) 
//       {
// 	maximalStressValue=StressTensor[1][1];
// 	Istress=1;
//       }
//     if (std::abs(StressTensor[2][2])>std::abs(maximalStressValue)) 
//       {
// 	maximalStressValue=StressTensor[2][2];
// 	Istress=2;
//       }
    
//        // 2nd maximal stress direction/value
//       double maximalStressValue2;
//       int Istress2,Istress3;
//       if (Istress==0) {
//         Istress2=1;
//         Istress3=2;
//       }
//       if (Istress==1) {
//         Istress2=0;
//         Istress3=2;
//       }
//       if (Istress==2) {
//         Istress2=0;
//         Istress3=1;
//       }
//       if(std::abs(StressTensor[Istress3][Istress3])>std::abs(StressTensor[Istress2][Istress2])) {
//         Istress2=Istress3;
//       }
//       maximalStressValue2=StressTensor[Istress2][Istress2];
      
      
    
    
//     // storing a measure for stress anisotropy in cell vector
//     if (std::abs(maximalStressValue)<0.000001) cellData[cellIndex][variableIndex(0,3)]=0;
//     if (std::abs(maximalStressValue)>= 0.000001) cellData[cellIndex][variableIndex(0,3)]=1-std::abs(maximalStressValue2/maximalStressValue);
    
    
        
//     // ---------------angles stress with circumferental direction- begin
//    double  temp=std::sqrt(eigenVectorStress[0][Istress]*eigenVectorStress[0][Istress]+                                // 14 --> Stress
// 		   eigenVectorStress[1][Istress]*eigenVectorStress[1][Istress] );
//     if(temp<0.000001){
//       cellData[cellIndex][14]=pi/2;
//     }
//     else{
//       cellData[cellIndex][14]=std::atan(eigenVectorStress[2][Istress]/temp);
//     } 
//     // ---------------angles stress with circumferental direction- end
    
    
    


//     if (numVariableIndexLevel()==3 && (numVariableIndex(2)==1 || numVariableIndex(2)==2)) { // storing maximal stress
//       if (dimension==2)
// 	{
// 	  cellData[cellIndex][variableIndex(2,0)]  =eigenVectorStress[0][Istress];
// 	  cellData[cellIndex][variableIndex(2,0)+1]=eigenVectorStress[1][Istress];
// 	  cellData[cellIndex][variableIndex(2,0)+3]=maximalStressValue;  //maximal Stress Value is stored after its eigenvector
// 	}
//       if (dimension==3)
// 	{
// 	  cellData[cellIndex][variableIndex(2,0)]  =eigenVectorStress[0][Istress];
// 	  cellData[cellIndex][variableIndex(2,0)+1]=eigenVectorStress[1][Istress];
// 	  cellData[cellIndex][variableIndex(2,0)+2]=eigenVectorStress[2][Istress];
// 	  cellData[cellIndex][variableIndex(2,0)+3]=maximalStressValue;  //maximal Stress Value is stored after its eigenvector
// 	}
//     }
    
//     if (numVariableIndexLevel()==3 && numVariableIndex(2)==2 ) { // storing 2nd maximal stress
//       if (dimension==2)
// 	{
// 	  cellData[cellIndex][variableIndex(2,1)]  =eigenVectorStress[0][Istress2];
// 	  cellData[cellIndex][variableIndex(2,1)+1]=eigenVectorStress[1][Istress2];
// 	  cellData[cellIndex][variableIndex(2,1)+3]=maximalStressValue2;  //2nd maximal Stress Value is stored after its eigenvector
// 	}
//       if (dimension==3)
// 	{
// 	  cellData[cellIndex][variableIndex(2,1)]  =eigenVectorStress[0][Istress2];
// 	  cellData[cellIndex][variableIndex(2,1)+1]=eigenVectorStress[1][Istress2];
// 	  cellData[cellIndex][variableIndex(2,1)+2]=eigenVectorStress[2][Istress2];
// 	  cellData[cellIndex][variableIndex(2,1)+3]=maximalStressValue2;  //2nd maximal Stress Value is stored after its eigenvector
// 	}
//     }
    
    
//   }//-----------------------------

//     cellData[0][variableIndex(0,5)]=totalEnergyIso ;
//     cellData[0][variableIndex(0,6)]=totalEnergyAniso ;
  // size_t numVertices = T.numVertex();
  // for(size_t vertexIndex=0; vertexIndex<numVertices; ++vertexIndex){ // stimating volume for equilibrium 
  //   TotalVolume +=std::sqrt(vertexData[vertexIndex][0]*vertexData[vertexIndex][0] +
  //                           vertexData[vertexIndex][1]*vertexData[vertexIndex][1] +
  //                           vertexData[vertexIndex][2]*vertexData[vertexIndex][2] );
  // }
  // deltaVolume=TotalVolume-cellData[0][25];
  // cellData[0][25]=TotalVolume;
  // cellData[0][24]=deltaVolume;
// }

// void VertexFromTRBSMTbending::initiate(Tissue &T,
// 				       DataMatrix &cellData,
// 				       DataMatrix &wallData,
// 				       DataMatrix &vertexData,
// 				       DataMatrix &cellDerivs,
// 				       DataMatrix &wallDerivs,
// 				       DataMatrix &vertexDerivs) {
//   size_t numWall=wallData.size();

//   double pi=3.14159265359;
  
//   for (size_t wallIndex=0; wallIndex<numWall; ++wallIndex) 
//     { // initiating with resting angle for each wall 
   

//       if( T.wall(wallIndex).cell1() != T.background() && T.wall(wallIndex).cell2() != T.background())
// 	{
	  
	  
// 	  size_t v1 = T.wall(wallIndex).vertex1()->index();
// 	  size_t v2 = T.wall(wallIndex).vertex2()->index();
	  
// 	  cell1index=T.wall(walIndex).cell1()->index();
// 	  cell2index=T.wall(walIndex).cell2()->index();
	  
// 	  size_t v1C1 = T.cell(cell1Index).vertex(0)->index();
// 	  size_t v2C1 = T.cell(cell1Index).vertex(1)->index();
// 	  size_t v3C1 = T.cell(cell1Index).vertex(2)->index();
	  
// 	  size_t v1C2 = T.cell(cell2Index).vertex(0)->index();
// 	  size_t v2C2 = T.cell(cell2Index).vertex(1)->index();
// 	  size_t v3C2 = T.cell(cell2Index).vertex(2)->index();
// 	  if (v1==v1C1)
// 	    { 
// 	      if (v2==v2C1){V1=v1; V2=v2}      
// 	      if (v2==v3C1){V1=v2; V2=v1}
// 	    }
// 	  if (v1==v2C1)
// 	    { 
// 	      if (v2==v3C1){V1=v1; V2=v2}      
// 	      if (v2==v1C1){V1=v2; V2=v1}
// 	    }
// 	  if (v1==v3C1)
// 	    { 
// 	      if (v2==v1C1){V1=v1; V2=v2}      
// 	      if (v2==v2C1){V1=v2; V2=v1}
// 	    }
// 	  // V1->V2 wall vector counterclockwise for cell 1 
// 	  double wallVector[3]={0,0,0};
// 	  wallVector[0]=vertexData[V2][0]-vertexData[V1][0];
// 	  wallVector[1]=vertexData[V2][1]-vertexData[V1][1];
// 	  wallVector[2]=vertexData[V2][2]-vertexData[V1][2];

 
// 	  double normal_1[3]={0,0,0};
// 	  double normal_2[3]={0,0,0};
// 	  normal_1[0]= (vertexData[v2C1][1]-vertexData[v1C1][1])*(vertexData[v3C1][2]-vertexData[v1C1][2])
// 	              -(vertexData[v2C1][2]-vertexData[v1C1][2])*(vertexData[v3C1][1]-vertexData[v1C1][1]);
// 	  normal_1[1]= (vertexData[v2C1][2]-vertexData[v1C1][2])*(vertexData[v3C1][0]-vertexData[v1C1][0])
// 	              -(vertexData[v2C1][0]-vertexData[v1C1][0])*(vertexData[v3C1][2]-vertexData[v1C1][2]);
// 	  normal_1[2]= (vertexData[v2C1][0]-vertexData[v1C1][0])*(vertexData[v3C1][1]-vertexData[v1C1][1])
// 	              -(vertexData[v2C1][1]-vertexData[v1C1][1])*(vertexData[v3C1][0]-vertexData[v1C1][0]);
// 	  double tmp=std::sqrt(normal_1[0]*normal_1[0]+normal_1[1]*normal_1[1]+normal_1[2]*normal_1[2]);
//           normal_1[0]/=tmp;          
// 	  normal_1[1]/=tmp;          
// 	  normal_1[2]/=tmp;          

	  
// 	  normal_2[0]= (vertexData[v2C2][1]-vertexData[v1C2][1])*(vertexData[v3C2][2]-vertexData[v1C2][2])
// 	              -(vertexData[v2C2][2]-vertexData[v1C2][2])*(vertexData[v3C2][1]-vertexData[v1C2][1]);
// 	  normal_2[1]= (vertexData[v2C2][2]-vertexData[v1C2][2])*(vertexData[v3C2][0]-vertexData[v1C2][0])
// 	              -(vertexData[v2C2][0]-vertexData[v1C2][0])*(vertexData[v3C2][2]-vertexData[v1C2][2]);
// 	  normal_2[2]= (vertexData[v2C2][0]-vertexData[v1C2][0])*(vertexData[v3C2][1]-vertexData[v1C2][1])
// 	              -(vertexData[v2C2][1]-vertexData[v1C2][1])*(vertexData[v3C2][0]-vertexData[v1C2][0]);
// 	  tmp=std::sqrt(normal_2[0]*normal_2[0]+normal_2[1]*normal_2[1]+normal_2[2]*normal_2[2]);
// 	  normal_2[0]/=tmp;          
// 	  normal_2[1]/=tmp;          
// 	  normal_2[2]/=tmp;          
          

// 	  tmp= wallVector[0]*(normal_1[1]*normal_2[2]-normal_1[2]*normal_2[1])
//               +wallVector[1]*(normal_1[2]*normal_2[0]-normal_1[0]*normal_2[2])
// 	      +wallVector[2]*(normal_1[0]*normal_2[1]-normal_1[1]*normal_2[0]);

// 	  double teta=std::acos(normal_1[0]*normal_2[0]+normal_1[1]*normal_2[1]+normal_1[2]*normal_2[2]);

// 	  if (tmp>0)  wallData[wallIndex][2]=pi+teta;
// 	  if (tmp<0)  wallData[wallIndex][2]=pi-teta;
// 	  // (normal1 x normal2).wallvector > 0  => teta > 180 :: teta= pi+acos(normal_1.normal_2)
// 	  //                                   else teta < 180 :: teta= pi-acos(normal_1.normal_2) 
// 	}
      

//       else{
// 	wallData[wallIndex][2]=-10;// for walls at the boundary   
//       }
      
      
//     }
  
  
  
// }
