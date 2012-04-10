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
  if( indValue.size()!=2 || indValue[0].size()!=1 || indValue[1].size()!=1 ) { 
    std::cerr << "VertexFromTRBScenterTriangulation::"
	      << "VertexFromTRBScenterTriangulation() "
	      << "Wall length index is given in first level." << std::endl
	      << "Start of additional Cell variable indices (center(x,y,z) "
	      << "L_1,...,L_n, n=num vertex) is given in second level." 
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
  
  
  for (size_t i=0 ; i<numCells ; ++i) {
    size_t numWalls = T.cell(i).numWall(); 
    
    if(  T.cell(i).numVertex()!= numWalls ) {
      std::cerr << "VertexFromTRBScenterTriangulation::derivs() same number of vertices and walls."
		<< " Not for cells with " << T.cell(i).numWall() << " walls and "
		<< T.cell(i).numVertex() << " vertices!"	
		<< std::endl;
      exit(-1);
    }
    
    double young = parameter(0);
    double poisson =parameter(1);
    
    double StrainCellGlobal[3][3]={{0,0,0},{0,0,0},{0,0,0}};
    // One triangle per 'vertex' in cyclic order
    for (size_t k=0; k<numWalls; ++k) { 
      size_t kPlusOneMod = (k+1)%numWalls;
      //size_t v1 = com;
      size_t v2 = T.cell(i).vertex(k)->index();
      size_t v3 = T.cell(i).vertex(kPlusOneMod)->index();
      //size_t w1 = internal k
      size_t w2 = T.cell(i).wall(k)->index();
      //size_t w3 = internal k+1

      // Position matrix holds in rows positions for com, vertex(k), vertex(k+1)
      DataMatrix position(3,vertexData[v2]);
      for (size_t d=0; d<dimension; ++d)
	position[0][d] = cellData[i][comIndex+d]; // com position
      //position[1] = vertexData[v2]; // given by initiation
      position[2] = vertexData[v3];
      
      // Resting lengths are from com-vertex(k), vertex(k)-vertex(k+1) (wall(k)), com-vertex(k+1)
      std::vector<double> restingLength(numWalls);
      restingLength[0] = cellData[i][lengthInternalIndex + k];
      restingLength[1] = wallData[w2][wallLengthIndex];
      restingLength[2] = cellData[i][lengthInternalIndex + kPlusOneMod];
      
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
      
      // Area of the element (using Heron's formula)                                      
      double Area=std::sqrt( ( restingLength[0]+restingLength[1]+restingLength[2])*
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
      
      
      //Strain tensor- reverse calculation- sign of eigen values should be changed as the roles of current and resting shape are switched  
      //Shape vectors in   (counterclockwise ordering of nodes/edges)
      double CurrentAngle1=std::acos(  (length[0]*length[0]+length[1]*length[1]-length[2]*length[2])/
                           (length[0]*length[1]*2)    );

      double Pa=std::cos(CurrentAngle1)*length[0];
      double Pc=std::sin(CurrentAngle1)*length[0];
      double Pb=length[1];
      // shape vector matrix = inverse of coordinate matrix      
      double ShapeVector[3][2]={ {  0   ,       1/Pc      }, 
                                 {-1/Pb , (Pa-Pb)/(Pb*Pc) },       
                                 { 1/Pb ,     -Pa/(Pb*Pc) }  };
            
      
      // //Strain tensor  (clockwise ordering of nodes/edges)
      // double CurrentAngle2=std::acos(  (length[1]*length[1]+length[2]*length[2]-length[0]*length[0])/
      //                                  (length[1]*length[2]*2)    );

      // double Pa=std::cos(CurrentAngle2)*length[2];
      // double Pb=length[1];
      // double Pc=std::sin(CurrentAngle2)*length[2];
      
      // double ShapeVector[3][2]={ {  0   ,       1/Pc      }, 
      //                            { 1/Pb ,     -Pa/(Pb*Pc) },       
      //                            {-1/Pb , (Pa-Pb)/(Pb*Pc) }  };


      //square of radius of circumstancing circle in resting shape
      double Rcirc2=(0.25*restingLength[0]*restingLength[1]*restingLength[2]/Area)*(0.25*restingLength[0]*restingLength[1]*restingLength[2]/Area);  
      
      double StrainTensor[3][3];
      int kk;
      double QrQs;
      StrainTensor[0][0]=0;
      StrainTensor[0][1]=0;
      StrainTensor[1][0]=0;
      StrainTensor[1][1]=0;

      for ( int r=0 ; r<3 ; ++r )   // QrQs is  in resting shape i.e. PrPs and shape vectors are in the current shape instead 
        { for ( int s=0 ; s<3 ; ++s )
            {     
              if ((r==0 && s==1)||(r==1 && s==0)) kk=0;
              if ((r==0 && s==2)||(r==2 && s==0)) kk=2;
              if ((r==1 && s==2)||(r==2 && s==1)) kk=1; 
              if ( s!=r ) QrQs=Rcirc2-(restingLength[kk]*restingLength[kk])*0.5;
              if ( s==r ) QrQs=Rcirc2;
              
              StrainTensor[0][0]=StrainTensor[0][0]+QrQs*ShapeVector[r][0]*ShapeVector[s][0];  
              StrainTensor[0][1]=StrainTensor[0][1]+QrQs*ShapeVector[r][0]*ShapeVector[s][1];
              StrainTensor[1][0]=StrainTensor[1][0]+QrQs*ShapeVector[r][1]*ShapeVector[s][0];
              StrainTensor[1][1]=StrainTensor[1][1]+QrQs*ShapeVector[r][1]*ShapeVector[s][1];
            }
        }

      StrainTensor[0][0]=StrainTensor[0][0]-1; // E=(C-I)/2 
      StrainTensor[1][1]=StrainTensor[1][1]-1;
      
      StrainTensor[0][0]=0.5*StrainTensor[0][0];
      StrainTensor[0][1]=0.5*StrainTensor[0][1];
      StrainTensor[1][0]=0.5*StrainTensor[1][0];
      StrainTensor[1][1]=0.5*StrainTensor[1][1];

      StrainTensor[0][2]=0;  // adding the 3rd dimension which is zero as the tensor is still in the element plane
      StrainTensor[1][2]=0;
      StrainTensor[2][2]=0;
      StrainTensor[2][0]=0;
      StrainTensor[2][1]=0;
      

      //rotation matrix to go to global coordinate system based on counterclockwise ordering
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
      
      // std::cerr <<" Sxx  "<< StrainTensor[0][0] <<" Sxy  "<< StrainTensor[0][1] << std::endl
      //           <<" Syx  "<< StrainTensor[1][0] <<" Syy  "<< StrainTensor[1][1] << std::endl<<std::endl;
     
      for (int r=0 ; r<3 ; r++) {
        for (int s=0 ; s<3 ; s++) {   
          StrainCellGlobal[r][s]= StrainCellGlobal[r][s]+StrainTensor[r][s];
        }
      }
      
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
      
      cellDerivs[i][comIndex  ] += Force[0][0];
      cellDerivs[i][comIndex+1] += Force[0][1];
      cellDerivs[i][comIndex+2] += Force[0][2];
      
      vertexDerivs[v2][0] += Force[1][0];
      vertexDerivs[v2][1] += Force[1][1];
      vertexDerivs[v2][2] += Force[1][2];
      
      vertexDerivs[v3][0] += Force[2][0];
      vertexDerivs[v3][1] += Force[2][1];
      vertexDerivs[v3][2] += Force[2][2];
    }
    for (int r=0 ; r<3 ; r++) {
      for (int s=0 ; s<3 ; s++) {   
        StrainCellGlobal[r][s]= StrainCellGlobal[r][s]/numWalls;
      }
    }
    // std::cerr <<" Sxx  "<< StrainCellGlobal[0][0] <<" Sxy  "<< StrainCellGlobal[0][1] <<" Sxz  "<< StrainCellGlobal[0][2] << std::endl
    //           <<" Syx  "<< StrainCellGlobal[1][0] <<" Syy  "<< StrainCellGlobal[1][1] <<" Syz  "<< StrainCellGlobal[1][2] << std::endl
    //           <<" Szx  "<< StrainCellGlobal[2][0] <<" Szy  "<< StrainCellGlobal[2][1] <<" Szz  "<< StrainCellGlobal[2][2] << std::endl <<std::endl;
    
    // // eigenvalue/eigenvectors of averaged strain tensor in global coordinate system. (Jacobi method)
    
    double eigenVector[3][3]={{1,0,0},{0,1,0},{0,0,1}};
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
        for (int r=0 ; r<3 ; r++) {
          for (int s=0 ; s<3 ; s++) {
            for(int w=0 ; w<3 ; w++) {
              tempStrain[r][s]=tempStrain[r][s]+StrainCellGlobal[r][w]*tempRot[w][s];
            }
          }
        }

        for (int r=0 ; r<3 ; r++) {
          for (int s=0 ; s<3 ; s++) {
            StrainCellGlobal[r][s]=0;
          }
        }

        for (int r=0 ; r<3 ; r++) {
          for (int s=0 ; s<3 ; s++) {
            for(int w=0 ; w<3 ; w++) {
              StrainCellGlobal[r][s]=StrainCellGlobal[r][s]+tempRot[w][r]*tempStrain[w][s];
            }
          }
        }
        
        for (int r=0 ; r<3 ; r++) {
          for (int s=0 ; s<3 ; s++) {
            tempStrain[r][s]=eigenVector[r][s];
          }
        }
        eigenVector= {{0,0,0},{0,0,0},{0,0,0}}; 
        for (int r=0 ; r<3 ; r++) {
          for (int s=0 ; s<3 ; s++) {
            for(int w=0 ; w<3 ; w++) {
              eigenVector[r][s]=eigenVector[r][s]+tempStrain[r][w]*tempRot[w][s];
            }
          }
        }
    }
    // maximal strain direction
    double maximalStrainValue=-StrainCellGlobal[0][0];
    I=0;
    if (-StrainCellGlobal[1][1]>maximalStrainValue) {
      maximalStrainValue=-StrainCellGlobal[1][1];
      I=1;
    }
    if (-StrainCellGlobal[2][2]>maximalStrainValue) {
      maximalStrainValue=-StrainCellGlobal[2][2];
      I=2;
    }
    // std::cerr<<"maximal direction "<< eigenVector[0][I] <<" "<< eigenVector[1][I] <<" "<< eigenVector[2][I] <<std::endl;  
    // std::cerr<<"maximal strain value "<< maximalStrainValue <<std::endl;  
    
    if (dimension==2){
      cellData[i][0]=eigenVector[0][I];
      cellData[i][1]=eigenVector[1][I];
    }
    if (dimension==3){
      cellData[i][0]=eigenVector[0][I];
      cellData[i][1]=eigenVector[1][I];
      cellData[i][2]=eigenVector[2][I];
      //cellData[i][3]=10*maximalStrainValue;
    }
    
    // maximal Strain value shoul be used
    
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
  
  for (size_t i=0; i<numCell; ++i) {
    size_t numInternalWall = T.cell(i).numVertex();
    cellData[i].resize(numVariable+dimension+numInternalWall);
    cellDerivs[i].resize(numVariable+dimension+numInternalWall);
    com = T.cell(i).positionFromVertex(vertexData);
    // Set center position to com of the cell
    for (size_t d=0; d<dimension; ++d)
      cellData[i][numVariable+d] = com[d];    
    // Set internal wall lengths to the distance btw com and the vertex
    for (size_t k=0; k<numInternalWall; ++k) {
      Vertex *tmpVertex = T.cell(i).vertex(k); 
      size_t vertexIndex = tmpVertex->index();
      double distance = std::sqrt( (com[0]-vertexData[vertexIndex][0])*
				   (com[0]-vertexData[vertexIndex][0])+
				   (com[1]-vertexData[vertexIndex][1])*
				   (com[1]-vertexData[vertexIndex][1])+
				   (com[2]-vertexData[vertexIndex][2])*
				   (com[2]-vertexData[vertexIndex][2]) );   
      cellData[i][numVariable+dimension+k] = distance;
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
	      << "Uses two parameters Young_modulus_min, Young_modulus_max, "
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
  
  for (size_t i=0 ; i<numCells ; ++i) {
    size_t numWalls = T.cell(i).numWall(); 
    
    if(  T.cell(i).numVertex()!= numWalls ) {
      std::cerr << "VertexFromTRBScenterTriangulationConcentrationHill::derivs() "
		<< "same number of vertices and walls."
		<< " Not for cells with " << T.cell(i).numWall() << " walls and "
		<< T.cell(i).numVertex() << " vertices!"	
		<< std::endl;
      exit(-1);
    }
    
    double young = parameter(0) + 
      parameter(1)*Kpow/( Kpow+std::pow(cellData[i][concIndex],parameter(4)) );
    double poisson =parameter(1);
    
    // One triangle per 'vertex' in cyclic order
    for (size_t k=0; k<numWalls; ++k) { 
      size_t kPlusOneMod = (k+1)%numWalls;
      //size_t v1 = com;
      size_t v2 = T.cell(i).vertex(k)->index();
      size_t v3 = T.cell(i).vertex(kPlusOneMod)->index();
      //size_t w1 = internal k
      size_t w2 = T.cell(i).wall(k)->index();
      //size_t w3 = internal k+1
      
      // Position matrix holds in rows positions for com, vertex(k), vertex(k+1)
      DataMatrix position(3,vertexData[v2]);
      for (size_t d=0; d<dimension; ++d)
	position[0][d] = cellData[i][comIndex+d]; // com position
      //position[1] = vertexData[v2]; // given by initiation
      position[2] = vertexData[v3];
      
      // Resting lengths are from com-vertex(k), vertex(k)-vertex(k+1) (wall(k)), com-vertex(k+1)
      std::vector<double> restingLength(numWalls);
      restingLength[0] = cellData[i][lengthInternalIndex + k];
      restingLength[1] = wallData[w2][wallLengthIndex];
      restingLength[2] = cellData[i][lengthInternalIndex + kPlusOneMod];
      
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
      
      cellDerivs[i][comIndex]+= Force[0][0];
      cellDerivs[i][comIndex+1]+= Force[0][1];
      cellDerivs[i][comIndex+2]+= Force[0][2];
      
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
  
  for (size_t i=0; i<numCell; ++i) {
    size_t numInternalWall = T.cell(i).numVertex();
    cellData[i].resize(numVariable+dimension+numInternalWall);
    cellDerivs[i].resize(numVariable+dimension+numInternalWall);
    com = T.cell(i).positionFromVertex(vertexData);
    // Set center position to com of the cell
    for (size_t d=0; d<dimension; ++d)
      cellData[i][numVariable+d] = com[d];    
    // Set internal wall lengths to the distance btw com and the vertex
    for (size_t k=0; k<numInternalWall; ++k) {
      Vertex *tmpVertex = T.cell(i).vertex(k); 
      size_t vertexIndex = tmpVertex->index();
      double distance = std::sqrt( (com[0]-vertexData[vertexIndex][0])*
				   (com[0]-vertexData[vertexIndex][0])+
				   (com[1]-vertexData[vertexIndex][1])*
				   (com[1]-vertexData[vertexIndex][1])+
				   (com[2]-vertexData[vertexIndex][2])*
				   (com[2]-vertexData[vertexIndex][2]) );   
      cellData[i][numVariable+dimension+k] = distance;
    }
  }
}

VertexFromTRBSMT::
VertexFromTRBSMT(std::vector<double> &paraValue, 
	       std::vector< std::vector<size_t> > 
	       &indValue ) 
{  
  // Do some checks on the parameters and variable indeces
  if( paraValue.size()!=4 ) {
    std::cerr << "VertexFromTRBSMT::"
	      << "VertexFromTRBSMT() "
	      << "Uses four parameters young modulus and poisson coefficients in "
	      << "longitudinal and transverse directions." << std::endl;
    exit(0);
  }
  if( indValue.size()!=1 || indValue[0].size()!=2 ) { 
    std::cerr << "VertexFromTRBSMT::"
	      << "VertexFromTRBSMT() "
	      << "Wall length index and cell MT direction start index"
	      << "given at first level." << std::endl;
    exit(0);
  }
  
  // Set the variable values
  setId("VertexFromTRBSMT");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  // Set the parameter identities
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "Y_mod_L";// Longitudinal components of parameters
  tmp[1] = "P_ratio_L";
  tmp[2] = "Y_mod_T";// Transverse components of parameters
  tmp[3] = "P_ratio_T";

  setParameterId( tmp );
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
  size_t wallLengthIndex = variableIndex(0,0);
  size_t numWalls = 3; // defined only for triangles at the moment
  
  for( size_t cellIndex=0 ; cellIndex<numCells ; ++cellIndex ) {
    if( T.cell(cellIndex).numWall() != numWalls ) {
      std::cerr << "VertexFromTRBSMT::derivs() only defined for triangular cells."
		<< " Not for cells with " << T.cell(cellIndex).numWall() << " walls!"
		<< std::endl;
      exit(-1);
    }
    
    double youngL   = parameter(0);//two more parameters here 
    double poissonL = parameter(1);
    double youngT   = parameter(2);
    double poissonT = parameter(3);
    
    size_t v1 = T.cell(cellIndex).vertex(0)->index();
    size_t v2 = T.cell(cellIndex).vertex(1)->index();
    size_t v3 = T.cell(cellIndex).vertex(2)->index();
    size_t w1 = T.cell(cellIndex).wall(0)->index();
    size_t w2 = T.cell(cellIndex).wall(1)->index();
    size_t w3 = T.cell(cellIndex).wall(2)->index();

    // std::cerr<< "cell "<< cellIndex<< " vertices  "<< v1<<" "<< v2 << " "<< v3 << " walls  "<< w1 <<" "<< w2 << " "<< w3<< std::endl;
    

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
    double lambdaL=youngL*poissonL/(1-poissonL*poissonL);
    double mioL=youngL/(1+poissonL);
    double lambdaT=youngT*poissonT/(1-poissonT*poissonT);
    double mioT=youngT/(1+poissonT);
    
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
    tensileStiffness[0]=(2*cotan[2]*cotan[2]*(lambdaT+mioT)+mioT)*temp;
    tensileStiffness[1]=(2*cotan[0]*cotan[0]*(lambdaT+mioT)+mioT)*temp;
    tensileStiffness[2]=(2*cotan[1]*cotan[1]*(lambdaT+mioT)+mioT)*temp;
    
    //Angular Stiffness
    std::vector<double> angularStiffness(3);
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
      
      // //Strain tensor  (clockwise ordering of nodes/edges)
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
      // AnisoCurrGlob[0] = cellData[cellIndex][variableIndex(0,1)];
      // AnisoCurrGlob[1] = cellData[cellIndex][variableIndex(0,1)+1];
      // AnisoCurrGlob[2] = cellData[cellIndex][variableIndex(0,1)+2];
     
      // static aniso direction
      // AnisoCurrGlob[0] =1;
      // AnisoCurrGlob[1] =0;
      // AnisoCurrGlob[2] =0;

      // dynamic aniso direction      
      AnisoCurrGlob[0] =cellData[cellIndex][0];
      AnisoCurrGlob[1] =cellData[cellIndex][1];
      AnisoCurrGlob[2] =cellData[cellIndex][2];

      //std::cerr<< "cell "<< cellIndex<< std::endl;
      //std::cerr<<position[0][0]<<"  "<<position[0][1]<<"  "<<position[0][2]<< std::endl;
      //std::cerr<<position[1][0]<<"  "<<position[1][1]<<"  "<<position[1][2]<< std::endl;
      //std::cerr<<position[2][0]<<"  "<<position[2][1]<<"  "<<position[2][2]<< std::endl;      
      //std::cerr<< "cell "<< cellIndex<< " anisoVector current global "<<AnisoCurrGlob[0]<<"  "<<AnisoCurrGlob[1]<<"  "<<AnisoCurrGlob[2]<< std::endl;

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
      
      //---- Anisotropic Correction-------------------------------
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
         
        // a sample collection of fixed angles 

        // if (cellIndex==1){
        //   teta[0]=0;
        //   teta[1]=(3*3.14159265)/4;
        //   teta[2]=3.14159265/2;
        // }
        // if (cellIndex==0){
        //   teta[0]=3.14159265;
        //   teta[1]=(3.14159265)/4;
        //   teta[2]=3.14159265/2;
        // }




      // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> STRESS TENSOR (BEGIN) <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      
      // deformation gradiant tensor F =Sigma i=1,2,3 Di x Qi
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
      
      double DeformGrad[2][2]={{0,0},{0,0}}; // F
      for ( int i=0 ; i<3 ; ++i ) {
        DeformGrad[0][0]=DeformGrad[0][0]+ShapeVectorResting[i][0]*positionLocal[i][0];
        DeformGrad[1][0]=DeformGrad[1][0]+ShapeVectorResting[i][1]*positionLocal[i][0];
        DeformGrad[0][1]=DeformGrad[0][1]+ShapeVectorResting[i][0]*positionLocal[i][1];
        DeformGrad[1][1]=DeformGrad[1][1]+ShapeVectorResting[i][1]*positionLocal[i][1];
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



      // double B[2][2]={{0,0},{0,0}};// another way for B
      // double DrDs;
      
      // for ( int r=0 ; r<3 ; ++r )  
      //   { for ( int s=0 ; s<3 ; ++s )
      //       {     
      //         if ((r==0 && s==1)||(r==1 && s==0)) DrDs=-0.5*cotan[2]/restingArea;
      //         if ((r==0 && s==2)||(r==2 && s==0)) DrDs=-0.5*cotan[1]/restingArea;
      //         if ((r==1 && s==2)||(r==2 && s==1)) DrDs=-0.5*cotan[0]/restingArea; 
      //         if (r==0 && s==0)  DrDs=0.25*restingLength[1]*restingLength[1] / (restingArea*restingArea);
      //         if (r==1 && s==1)  DrDs=0.25*restingLength[2]*restingLength[2] / (restingArea*restingArea);
      //         if (r==2 && s==2)  DrDs=0.25*restingLength[0]*restingLength[0] / (restingArea*restingArea);
                  
      //         B[0][0]=B[0][0]+DrDs*positionLocal[r][0]*positionLocal[s][0];  
      //         B[0][1]=B[0][1]+DrDs*positionLocal[r][0]*positionLocal[s][1];
      //         B[1][0]=B[1][0]+DrDs*positionLocal[r][1]*positionLocal[s][0];
      //         B[1][1]=B[1][1]+DrDs*positionLocal[r][1]*positionLocal[s][1];
      //       }
      //   }



      // std::cerr <<"stress tensor " << std::endl;
      // std::cerr <<" Sxx  "<< StressTensor[0][0] <<" Sxy  "<< StressTensor[0][1] <<" Sxz  "<< StressTensor[0][2] << std::endl
      //           <<" Syx  "<< StressTensor[1][0] <<" Syy  "<< StressTensor[1][1] <<" Syz  "<< StressTensor[1][2] << std::endl
      //           <<" Szx  "<< StressTensor[2][0] <<" Szy  "<< StressTensor[2][1] <<" Szz  "<< StressTensor[2][2] << std::endl <<std::endl;




      // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> STRESS TENSOR (END) <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      //Shape vectors in Current shape (counterclockwise ordering of nodes/edges)     ShapeVectorCurrent[3][3]  calculated above   
      //.............................. ( or clockwise ordering of nodes/edges)
          

      //square of radius of circumstancing circle in resting shape
      //double Rcirc2Resting=(0.25*restingLength[0]*restingLength[1]*restingLength[2]/Area)*(0.25*restingLength[0]*restingLength[1]*restingLength[2]/Area);  
      
      double StrainTensor[3][3];
      //1


      //     FOR TEST
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








      //  std::cerr <<"strain tensor " << std::endl;
      // std::cerr <<" Sxx  "<< StrainTensor[0][0] <<" Sxy  "<< StrainTensor[0][1] <<" Sxz  "<< StrainTensor[0][2] << std::endl
      //           <<" Syx  "<< StrainTensor[1][0] <<" Syy  "<< StrainTensor[1][1] <<" Syz  "<< StrainTensor[1][2] << std::endl
      //           <<" Szx  "<< StrainTensor[2][0] <<" Szy  "<< StrainTensor[2][1] <<" Szz  "<< StrainTensor[2][2] << std::endl <<std::endl;
    
      // // eigenvalue/eigenvectors of averaged strain tensor in global coordinate system. (Jacobi method)
      
      double eigenVectorStrain[3][3]={{1,0,0},{0,1,0},{0,0,1}};
      double pivot=1;
      double pi=3.1415;
      int I,J;
      double RotAngle,Si,Co;
      while (pivot>0.00001) {
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
        if (std::abs(StrainTensor[I][I]-StrainTensor[J][J])<0.00001) {
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

        eigenVectorStrain= {{0,0,0},{0,0,0},{0,0,0}}; 
	//HJ: commented due to compilation problems...all values set below anyway?
        for (int r=0 ; r<3 ; r++) {
          for (int s=0 ; s<3 ; s++) {
            for(int w=0 ; w<3 ; w++) {
              eigenVectorStrain[r][s]=eigenVectorStrain[r][s]+tempStrain[r][w]*tempRot[w][s];
            }
          }
        }
      }





      // maximal strain direction for REVERSE calculation of strain
      // double maximalStrainValue=-StressTensor[0][0];
      // I=0;
      // if (-StrainTensor[1][1]>maximalStrainValue) {
      //   maximalStrainValue=-StrainTensor[1][1];
      //   I=1;
      // }
      // if (-StrainTensor[2][2]>maximalStrainValue) {
      //   maximalStrainValue=-StrainTensor[2][2];
      //   I=2;
      // }
      
      
      // maximal strain direction
      double maximalStrainValue=StrainTensor[0][0];
      int Istrain=0;
      if (StrainTensor[1][1]>maximalStrainValue) 
        {
          maximalStrainValue=StrainTensor[1][1];
          Istrain=1;
        }
      if (StrainTensor[2][2]>maximalStrainValue) 
        {
          maximalStrainValue=StrainTensor[2][2];
          Istrain=2;
        }
      // std::cerr<<"maximal Strain direction "<< eigenVectorStrain[0][Istrain] <<" "<< eigenVectorStrain[1][Istrain] <<" "<< eigenVectorStrain[2][Istrain] <<std::endl;  
      // std::cerr<<"maximal Strain value "<< maximalStrainValue <<std::endl;  
      
      


      // eigenvalue/eigenvectors of averaged stress tensor in global coordinate system. (Jacobi method)
      
      double eigenVectorStress[3][3]={{1,0,0},{0,1,0},{0,0,1}};
      pivot=1;
      //double RotAngle,Si,Co;
      while (pivot>0.00001) {
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
        if (std::abs(StressTensor[I][I]-StressTensor[J][J])<0.00001) {
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
        for (int r=0 ; r<3 ; r++) {
          for (int s=0 ; s<3 ; s++) {
            for(int w=0 ; w<3 ; w++) {
              tempStress[r][s]=tempStress[r][s]+StressTensor[r][w]*tempRot[w][s];
            }
          }
        }
        
        for (int r=0 ; r<3 ; r++) {
          for (int s=0 ; s<3 ; s++) {
            StressTensor[r][s]=0;
          }
        }
        
        for (int r=0 ; r<3 ; r++) {
          for (int s=0 ; s<3 ; s++) {
            for(int w=0 ; w<3 ; w++) {
              StressTensor[r][s]=StressTensor[r][s]+tempRot[w][r]*tempStress[w][s];
            }
          }
        }
        
        for (int r=0 ; r<3 ; r++) 
          {
            for (int s=0 ; s<3 ; s++) 
              {
                tempStress[r][s]=eigenVectorStress[r][s];
              }
          }
        eigenVectorStress= {{0,0,0},{0,0,0},{0,0,0}}; 
        //HJ: commented due to compilation problems...all values set below anyway?
        for (int r=0 ; r<3 ; r++) {
          for (int s=0 ; s<3 ; s++) {
            for(int w=0 ; w<3 ; w++) {
              eigenVectorStress[r][s]=eigenVectorStress[r][s]+tempStress[r][w]*tempRot[w][s];
            }
          }
        }
      }
      
      
      // maximal stress direction
      double maximalStressValue=StressTensor[0][0];
      int Istress=0;
      if (StressTensor[1][1]>maximalStressValue) 
        {
          maximalStressValue=StressTensor[1][1];
          Istress=1;
        }
      if (StressTensor[2][2]>maximalStressValue) 
        {
          maximalStressValue=StressTensor[2][2];
          Istress=2;
        }
      // std::cerr<<"maximal Stress direction "<< eigenVectorStress[0][Istress] <<" "<< eigenVectorStress[1][Istress] <<" "<< eigenVectorStress[2][Istress] <<std::endl;  
      // std::cerr<<"maximal Stress value "<< maximalStressValue <<std::endl;  
      
     




      // storing normal dirrection to  strain in cellData
     
      // normal to the cell plane in global direction is Zcurrent[], vector product gives the perpendicular strain direction
      double PerpStrain[3];
      PerpStrain[0]=Zcurrent[1]*eigenVectorStrain[2][Istrain]-Zcurrent[2]*eigenVectorStrain[1][Istrain];
      PerpStrain[1]=Zcurrent[2]*eigenVectorStrain[0][Istrain]-Zcurrent[0]*eigenVectorStrain[2][Istrain];
      PerpStrain[2]=Zcurrent[0]*eigenVectorStrain[1][Istrain]-Zcurrent[1]*eigenVectorStrain[0][Istrain];
      temp=std::sqrt(PerpStrain[0]*PerpStrain[0]+PerpStrain[1]*PerpStrain[1]+PerpStrain[2]*PerpStrain[2]);     

      if(std::abs(temp)>0.001){ // if aniso vector is not perpendicular to the cell plane
        if (dimension==2)
          {
            cellData[cellIndex][0]=PerpStrain[0];        
            cellData[cellIndex][1]=PerpStrain[1];
          }
        if (dimension==3)
          {
            cellData[cellIndex][0]=PerpStrain[0];
            cellData[cellIndex][1]=PerpStrain[1];
            cellData[cellIndex][2]=PerpStrain[2];
            // cellData[cellIndex][3]=10*maximalStrainValue;  //NOTE maximal Strain and Stress Values should be used somewhere this is an option
          }
      }
      
      // storing maximal strain in cellData
      // if (dimension==2)
      //   {
      //     cellData[cellIndex][0]=eigenVectorStrain[0][Istrain];
      //     cellData[cellIndex][1]=eigenVectorStrain[1][Istrain];
      //   }
      // if (dimension==3)
      //   {
      //     cellData[cellIndex][0]=eigenVectorStrain[0][Istrain];
      //     cellData[cellIndex][1]=eigenVectorStrain[1][Istrain];
      //     cellData[cellIndex][2]=eigenVectorStrain[2][Istrain];
      //     // cellData[cellIndex][3]=10*maximalStrainValue;  //NOTE maximal Strain and Stress Values should be used somewhere this is an option
      //   }




      // storing maximal stress in cellData
      // if (dimension==2)
      //   {
      //     cellData[cellIndex][0]=eigenVectorStress[0][Istress];
      //     cellData[cellIndex][1]=eigenVectorStress[1][Istress];
      //   }
      // if (dimension==3)
      //   {
      //     cellData[cellIndex][0]=eigenVectorStress[0][Istress];
      //     cellData[cellIndex][1]=eigenVectorStress[1][Istress];
      //     cellData[cellIndex][2]=eigenVectorStress[2][Istress];
      //     // cellData[cellIndex][3]=10*maximalStressValue;  //NOTE maximal Strain and Stress Values should be used somewhere this is an option
      //   }
      
      
      
      
      
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
      
      // double temp=0.5/(restingArea*restingArea);
      // double temp1=-1/restingArea;

      // derIprim1[0][0]=(temp*restingLength[1]*restingLength[1])*position[0][0]
      //                +(temp1*cotan[2])*position[1][0]                                
      //                +(temp1*cotan[1])*position[2][0];
      // derIprim1[1][0]=(temp1*cotan[2])*position[0][0]                                
      //                +(temp*restingLength[2]*restingLength[2])*position[1][0]
      //                +(temp1*cotan[0])*position[2][0];
      // derIprim1[2][0]=(temp1*cotan[1])*position[0][0]                                
      //                +(temp1*cotan[0])*position[1][0]                                
      //                +(temp*restingLength[0]*restingLength[0])*position[2][0];
      
      // derIprim1[0][0]=(temp*restingLength[1]*restingLength[1])*position[0][1]
      //                +(temp1*cotan[2])*position[1][1]                                
      //                +(temp1*cotan[1])*position[2][1];
      // derIprim1[1][0]=(temp1*cotan[2])*position[0][1]                                
      //                +(temp*restingLength[2]*restingLength[2])*position[1][1]
      //                +(temp1*cotan[0])*position[2][1];
      // derIprim1[2][0]=(temp1*cotan[1])*position[0][1]                                
      //                +(temp1*cotan[0])*position[1][1]                                
      //                +(temp*restingLength[0]*restingLength[0])*position[2][1];
      
      // derIprim1[0][0]=(temp*restingLength[1]*restingLength[1])*position[0][2]
      //                +(temp1*cotan[2])*position[1][2]                                
      //                +(temp1*cotan[1])*position[2][2];
      // derIprim1[1][0]=(temp1*cotan[2])*position[0][2]                                
      //                +(temp*restingLength[2]*restingLength[2])*position[1][2]
      //                +(temp1*cotan[0])*position[2][2];
      // derIprim1[2][0]=(temp1*cotan[1])*position[0][2]                                
      //                +(temp1*cotan[0])*position[1][2]                                
      //                +(temp*restingLength[0]*restingLength[0])*position[2][2];
      
      // derIprim1[0][1]=(2*D1D1)*position[0][1]+(2*D1D2)*position[1][1]+(2*D0D1)*position[2][1];
      // derIprim1[1][1]=(2*D1D2)*position[0][1]+(2*D2D2)*position[1][1]+(2*D0D2)*position[2][1];
      // derIprim1[2][1]=(2*D1D0)*position[0][1]+(2*D2D0)*position[1][1]+(2*D0D0)*position[2][1];
      
      // derIprim1[0][2]=(2*D1D1)*position[0][2]+(2*D1D2)*position[1][2]+(2*D0D1)*position[2][2];
      // derIprim1[1][2]=(2*D1D2)*position[0][2]+(2*D2D2)*position[1][2]+(2*D0D2)*position[2][2];
      // derIprim1[2][2]=(2*D1D0)*position[0][2]+(2*D2D0)*position[1][2]+(2*D0D0)*position[2][2];
      

      
      int k;
      for ( int m=0 ; m<3 ; ++m )
        {
          for ( int coor=0 ; coor<3 ; ++coor ) 
            derIprim1[m][coor]=0;
          for ( int i=0 ; i<3 ; ++i )
            {            
              if ((i==0 && m==1)||(i==1 && m==0)) k=2;
              if ((i==0 && m==2)||(i==2 && m==0)) k=1;
              if ((i==1 && m==2)||(i==2 && m==1)) k=0; 
              //else {
              //std::cerr << "mechanicalTRBS::derivs() k not given a value..." << std::endl;
              //exit(-1);
              //}
              if (i!=m) DiDm=-0.5*cotan[k]/restingArea;
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
              if ((i==0 && j==1)||(i==1 && j==0)) k=2; 
              if ((i==0 && j==2)||(i==2 && j==0)) k=1; 
              if ((i==1 && j==2)||(i==2 && j==1)) k=0;  
              if (i!=j) QiQj=Rcirc2-(length[k]*length[k])*0.5; 
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
                        if ((r==0 && s==1)||(r==1 && s==0)) k=2;
                        if ((r==0 && s==2)||(r==2 && s==0)) k=1;
                        if ((r==1 && s==2)||(r==2 && s==1)) k=0; 
                        //if ( s!=r ) QrQs=Rcirc2-(length[k]*length[k])*0.5;
                        //if ( s==r ) QrQs=Rcirc2;
                        QrQs=position[r][0]*position[s][0]+position[r][1]*position[s][1]+position[r][2]*position[s][2];
                        
                        if ((n==0 && r==1)||(n==1 && r==0)) k=2; 
                        if ((n==0 && r==2)||(n==2 && r==0)) k=1; 
                        if ((n==1 && r==2)||(n==2 && r==1)) k=0;    
                        if ( n!=r )  DnDr=-0.5*cotan[k]/restingArea;
                        if ( n==r )  DnDr=0.25*restingLength[n]*restingLength[n] / (restingArea*restingArea);
                        
                        if ((s==0 && p==1)||(s==1 && p==0)) k=2; 
                        if ((s==0 && p==2)||(s==2 && p==0)) k=1; 
                        if ((s==1 && p==2)||(s==2 && p==1)) k=0;   
                        if ( s!=p ) DsDp=-0.5*cotan[k]/restingArea;
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
    //   }
    // else {// if anisotropy vector is perpendicular to the elements we put anisotropic correction 0
    //   for ( int i=0 ; i<3 ; ++i ) 
    //     for ( int j=0 ; j<3 ; ++j )
    //       deltaF[i][j]=0;
    // }  
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
  }
}
// 1
  // int kk;
      // double PrPs;
      // StrainTensor[0][0]=0;
      // StrainTensor[0][1]=0;
      // StrainTensor[1][0]=0;
      // StrainTensor[1][1]=0;

      // for ( int r=0 ; r<3 ; ++r )   // PrPs is in resting shape and shape vectors are in the current shape instead 
      //   { for ( int s=0 ; s<3 ; ++s )
      //       {     
      //         if ((r==0 && s==1)||(r==1 && s==0)) kk=0;
      //         if ((r==0 && s==2)||(r==2 && s==0)) kk=2;
      //         if ((r==1 && s==2)||(r==2 && s==1)) kk=1; 
      //         if ( s!=r ) PrPs=Rcirc2Resting-(restingLength[kk]*restingLength[kk])*0.5;
      //         if ( s==r ) PrPs=Rcirc2Resting;
              
      //         StrainTensor[0][0]=StrainTensor[0][0]+PrPs*ShapeVectorCurrent[r][0]*ShapeVectorCurrent[s][0];  
      //         StrainTensor[0][1]=StrainTensor[0][1]+PrPs*ShapeVectorCurrent[r][0]*ShapeVectorCurrent[s][1];
      //         StrainTensor[1][0]=StrainTensor[1][0]+PrPs*ShapeVectorCurrent[r][1]*ShapeVectorCurrent[s][0];
      //         StrainTensor[1][1]=StrainTensor[1][1]+PrPs*ShapeVectorCurrent[r][1]*ShapeVectorCurrent[s][1];
      //       }
      //   }
      
      // std::cerr <<"strain tensor " << std::endl;
      // std::cerr <<" Sxx  "<< StrainTensor[0][0] <<" Sxy  "<< StrainTensor[0][1]  << std::endl
      //           <<" Syx  "<< StrainTensor[1][0] <<" Syy  "<< StrainTensor[1][1] << std::endl;
    

VertexFromTRBScenterTriangulationMT::
VertexFromTRBScenterTriangulationMT(std::vector<double> &paraValue, 
	       std::vector< std::vector<size_t> > 
	       &indValue ) 
{  
  // Do some checks on the parameters and variable indeces
  if( paraValue.size()!=4 ) {
    std::cerr << "VertexFromTRBScenterTriangulationMT::"
	      << "VertexFromTRBScenterTriangulationMT() "
              << "Uses four parameters young modulus and poisson coefficients in "
	      << "longitudinal (MT) and transverse directions." << std::endl;
	      
    exit(0);
  }
  if( (indValue.size()!=2 && indValue.size()!=4) || 
      indValue[0].size()!=2 || indValue[1].size()!=1 ||
      (indValue.size()==4 && (indValue[2].size()!=0 && indValue[2].size()!=1)) ||
      (indValue.size()==4 && (indValue[3].size()!=0 && indValue[3].size()!=1)) 
      ) { 
    std::cerr << "VertexFromTRBScenterTriangulationMT::"
	      << "VertexFromTRBScenterTriangulationMT() "
	      << "Wall length index and MT direction initial index are given in first level." 
	      << std::endl
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
  setId("VertexFromTRBScenterTriangulationMT");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  // Set the parameter identities
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "Y_mod_L";// Longitudinal components of parameters
  tmp[1] = "P_ratio_L";
  tmp[2] = "Y_mod_T";// Transverse components of parameters
  tmp[3] = "P_ratio_T";
  setParameterId( tmp );
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

    double youngL   = parameter(0);//two more parameters here 
    double poissonL = parameter(1);
    double youngT   = parameter(2);
    double poissonT = parameter(3);
    
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
      //position[0][2] z for vertex 1 of the current element

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

            
      // Lame coefficients 
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
    std::vector<double> cotan(3);
    cotan[0] = 1.0/std::tan(Angle[0]);
    cotan[1] = 1.0/std::tan(Angle[1]);
    cotan[2] = 1.0/std::tan(Angle[2]);    
    //the force is calculated based on Transverse coefficients
    //Longitudinal coefficients are considered in deltaF
    tensileStiffness[0]=(2*cotan[2]*cotan[2]*(lambdaT+mioT)+mioT)*temp;
    tensileStiffness[1]=(2*cotan[0]*cotan[0]*(lambdaT+mioT)+mioT)*temp;
    tensileStiffness[2]=(2*cotan[1]*cotan[1]*(lambdaT+mioT)+mioT)*temp;
    
    //Angular Stiffness
    std::vector<double> angularStiffness(3);
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
      
      // //Strain tensor  (clockwise ordering of nodes/edges)
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
      // AnisoCurrGlob[0] = cellData[cellIndex][variableIndex(0,1)];  // this was written by Henrik, Behruz changed it
      // AnisoCurrGlob[1] = cellData[cellIndex][variableIndex(0,1)+1];
      // AnisoCurrGlob[2] = cellData[cellIndex][variableIndex(0,1)+2];
     
      // static aniso direction
      AnisoCurrGlob[0] =1;
      AnisoCurrGlob[1] =0;
      AnisoCurrGlob[2] =0;

      // dynamic aniso direction      
      // AnisoCurrGlob[0] =cellData[cellIndex][0];
      // AnisoCurrGlob[1] =cellData[cellIndex][1];
      // AnisoCurrGlob[2] =cellData[cellIndex][2];

      //std::cerr<< "cell "<< cellIndex<< std::endl;
      //std::cerr<<position[0][0]<<"  "<<position[0][1]<<"  "<<position[0][2]<< std::endl;
      //std::cerr<<position[1][0]<<"  "<<position[1][1]<<"  "<<position[1][2]<< std::endl;
      //std::cerr<<position[2][0]<<"  "<<position[2][1]<<"  "<<position[2][2]<< std::endl;      
      //std::cerr<< "cell "<< cellIndex<< " anisoVector current global "<<AnisoCurrGlob[0]<<"  "<<AnisoCurrGlob[1]<<"  "<<AnisoCurrGlob[2]<< std::endl;

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
      
      // deformation gradiant tensor F =Sigma i=1,2,3 Di x Qi
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
      
      double DeformGrad[2][2]={{0,0},{0,0}}; // F
      for ( int i=0 ; i<3 ; ++i ) {
        DeformGrad[0][0]=DeformGrad[0][0]+ShapeVectorResting[i][0]*positionLocal[i][0];
        DeformGrad[1][0]=DeformGrad[1][0]+ShapeVectorResting[i][1]*positionLocal[i][0];
        DeformGrad[0][1]=DeformGrad[0][1]+ShapeVectorResting[i][0]*positionLocal[i][1];
        DeformGrad[1][1]=DeformGrad[1][1]+ShapeVectorResting[i][1]*positionLocal[i][1];
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

      //  std::cerr <<"strain tensor in global coordinate system" << std::endl;
      // std::cerr <<" Sxx  "<< StrainTensor[0][0] <<" Sxy  "<< StrainTensor[0][1] <<" Sxz  "<< StrainTensor[0][2] << std::endl
      //           <<" Syx  "<< StrainTensor[1][0] <<" Syy  "<< StrainTensor[1][1] <<" Syz  "<< StrainTensor[1][2] << std::endl
      //           <<" Szx  "<< StrainTensor[2][0] <<" Szy  "<< StrainTensor[2][1] <<" Szz  "<< StrainTensor[2][2] << std::endl <<std::endl;

      //  std::cerr <<"stress tensor in global coordinate system" << std::endl;
      // std::cerr <<" Sxx  "<< StressTensor[0][0] <<" Sxy  "<< StressTensor[0][1] <<" Sxz  "<< StressTensor[0][2] << std::endl
      //           <<" Syx  "<< StressTensor[1][0] <<" Syy  "<< StressTensor[1][1] <<" Syz  "<< StressTensor[1][2] << std::endl
      //           <<" Szx  "<< StressTensor[2][0] <<" Szy  "<< StressTensor[2][1] <<" Szz  "<< StressTensor[2][2] << std::endl <<std::endl;


      // accumulating strain and stress tensors to be averaged later
      for (int r=0 ; r<3 ; r++) {
        for (int s=0 ; s<3 ; s++) {   
          StrainCellGlobal[r][s]= StrainCellGlobal[r][s]+restingArea*StrainTensor[r][s];
        }
      }

      for (int r=0 ; r<3 ; r++) {
        for (int s=0 ; s<3 ; s++) {   
          StressCellGlobal[r][s]= StressCellGlobal[r][s]+restingArea*StressTensor[r][s];
        }
      }

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
      
      // double temp=0.5/(restingArea*restingArea);
      // double temp1=-1/restingArea;

      // derIprim1[0][0]=(temp*restingLength[1]*restingLength[1])*position[0][0]
      //                +(temp1*cotan[2])*position[1][0]                                
      //                +(temp1*cotan[1])*position[2][0];
      // derIprim1[1][0]=(temp1*cotan[2])*position[0][0]                                
      //                +(temp*restingLength[2]*restingLength[2])*position[1][0]
      //                +(temp1*cotan[0])*position[2][0];
      // derIprim1[2][0]=(temp1*cotan[1])*position[0][0]                                
      //                +(temp1*cotan[0])*position[1][0]                                
      //                +(temp*restingLength[0]*restingLength[0])*position[2][0];
      
      // derIprim1[0][0]=(temp*restingLength[1]*restingLength[1])*position[0][1]
      //                +(temp1*cotan[2])*position[1][1]                                
      //                +(temp1*cotan[1])*position[2][1];
      // derIprim1[1][0]=(temp1*cotan[2])*position[0][1]                                
      //                +(temp*restingLength[2]*restingLength[2])*position[1][1]
      //                +(temp1*cotan[0])*position[2][1];
      // derIprim1[2][0]=(temp1*cotan[1])*position[0][1]                                
      //                +(temp1*cotan[0])*position[1][1]                                
      //                +(temp*restingLength[0]*restingLength[0])*position[2][1];
      
      // derIprim1[0][0]=(temp*restingLength[1]*restingLength[1])*position[0][2]
      //                +(temp1*cotan[2])*position[1][2]                                
      //                +(temp1*cotan[1])*position[2][2];
      // derIprim1[1][0]=(temp1*cotan[2])*position[0][2]                                
      //                +(temp*restingLength[2]*restingLength[2])*position[1][2]
      //                +(temp1*cotan[0])*position[2][2];
      // derIprim1[2][0]=(temp1*cotan[1])*position[0][2]                                
      //                +(temp1*cotan[0])*position[1][2]                                
      //                +(temp*restingLength[0]*restingLength[0])*position[2][2];
      
      // derIprim1[0][1]=(2*D1D1)*position[0][1]+(2*D1D2)*position[1][1]+(2*D0D1)*position[2][1];
      // derIprim1[1][1]=(2*D1D2)*position[0][1]+(2*D2D2)*position[1][1]+(2*D0D2)*position[2][1];
      // derIprim1[2][1]=(2*D1D0)*position[0][1]+(2*D2D0)*position[1][1]+(2*D0D0)*position[2][1];
      
      // derIprim1[0][2]=(2*D1D1)*position[0][2]+(2*D1D2)*position[1][2]+(2*D0D1)*position[2][2];
      // derIprim1[1][2]=(2*D1D2)*position[0][2]+(2*D2D2)*position[1][2]+(2*D0D2)*position[2][2];
      // derIprim1[2][2]=(2*D1D0)*position[0][2]+(2*D2D0)*position[1][2]+(2*D0D0)*position[2][2];
      

      
      int k;
      for ( int m=0 ; m<3 ; ++m )
        {
          for ( int coor=0 ; coor<3 ; ++coor ) 
            derIprim1[m][coor]=0;
          for ( int i=0 ; i<3 ; ++i )
            {            
              if ((i==0 && m==1)||(i==1 && m==0)) k=2;
              if ((i==0 && m==2)||(i==2 && m==0)) k=1;
              if ((i==1 && m==2)||(i==2 && m==1)) k=0; 
              //else {
              //std::cerr << "mechanicalTRBS::derivs() k not given a value..." << std::endl;
              //exit(-1);
              //}
              if (i!=m) DiDm=-0.5*cotan[k]/restingArea;
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
              if ((i==0 && j==1)||(i==1 && j==0)) k=2; 
              if ((i==0 && j==2)||(i==2 && j==0)) k=1; 
              if ((i==1 && j==2)||(i==2 && j==1)) k=0;  
              if (i!=j) QiQj=Rcirc2-(length[k]*length[k])*0.5; 
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
                        if ((r==0 && s==1)||(r==1 && s==0)) k=2;
                        if ((r==0 && s==2)||(r==2 && s==0)) k=1;
                        if ((r==1 && s==2)||(r==2 && s==1)) k=0; 
                        //if ( s!=r ) QrQs=Rcirc2-(length[k]*length[k])*0.5;
                        //if ( s==r ) QrQs=Rcirc2;
                        QrQs=position[r][0]*position[s][0]+position[r][1]*position[s][1]+position[r][2]*position[s][2];
                        
                        if ((n==0 && r==1)||(n==1 && r==0)) k=2; 
                        if ((n==0 && r==2)||(n==2 && r==0)) k=1; 
                        if ((n==1 && r==2)||(n==2 && r==1)) k=0;    
                        if ( n!=r )  DnDr=-0.5*cotan[k]/restingArea;
                        if ( n==r )  DnDr=0.25*restingLength[n]*restingLength[n] / (restingArea*restingArea);
                        
                        if ((s==0 && p==1)||(s==1 && p==0)) k=2; 
                        if ((s==0 && p==2)||(s==2 && p==0)) k=1; 
                        if ((s==1 && p==2)||(s==2 && p==1)) k=0;   
                        if ( s!=p ) DsDp=-0.5*cotan[k]/restingArea;
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
    //   }
    // else {// if anisotropy vector is perpendicular to the elements we put anisotropic correction 0
    //   for ( int i=0 ; i<3 ; ++i ) 
    //     for ( int j=0 ; j<3 ; ++j )
    //       deltaF[i][j]=0;
    // }  
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
    
    // adding TRBSMT forces to the total vertexDerivs
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
          
        eigenVectorStress= {{0,0,0},{0,0,0},{0,0,0}}; 
        //HJ: commented due to compilation problems...all values set below anyway?
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
      

      // storing maximal strain in cellData
      // if (dimension==2)
      //   {
      //     cellData[cellIndex][0]=eigenVectorStrain[0][Istrain];
      //     cellData[cellIndex][1]=eigenVectorStrain[1][Istrain];
      //   }
      // if (dimension==3)
      //   {
      //     cellData[cellIndex][0]=eigenVectorStrain[0][Istrain];
      //     cellData[cellIndex][1]=eigenVectorStrain[1][Istrain];
      //     cellData[cellIndex][2]=eigenVectorStrain[2][Istrain];
          // cellData[cellIndex][3]=10*maximalStrainValue;  //NOTE maximal Strain and Stress Values can be used - this is an option
      //  }



      if (numVariableIndexLevel()==4 && numVariableIndex(3) ) {
	// storing maximal stress in cellData
	if (dimension==2)
	  {
	    cellData[cellIndex][variableIndex(3,0)]=eigenVectorStress[0][Istress];
	    cellData[cellIndex][variableIndex(3,0)+1]=eigenVectorStress[1][Istress];
	  }
	if (dimension==3)
	  {
	    cellData[cellIndex][variableIndex(3,0)]=eigenVectorStress[0][Istress];
	    cellData[cellIndex][variableIndex(3,0)+1]=eigenVectorStress[1][Istress];
	    cellData[cellIndex][variableIndex(3,0)+2]=eigenVectorStress[2][Istress];
	    // cellData[cellIndex][3]=10*maximalStressValue;   //NOTE maximal Strain and Stress Values can be used this is an option
	  }
      }
  }      
}     


void VertexFromTRBScenterTriangulationMT::
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

