//
// Filename     : mechanicalTRBSAniso.cc
// Description  : Classes describing updates due to mechanical triangular biquadratic springs with anisotropic property corresponding to microtubals
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : July 2011
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
  tmp[0] = "Y_mod_L";                                            // Longitudinal components of parameters
  tmp[1] = "P_ratio_L";
  tmp[2] = "Y_mod_T";                                            // Transverse components of parameters
  tmp[3] = "P_ratio_T";

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
		<< " Not for cells with " << T.cell(i).numWall() << " walls!"
		<< std::endl;
      exit(-1);
    }
    
    double youngL   = parameter(0);                                                             //two more parameters here 
    double poissonL = parameter(1);
    double youngT   = parameter(2);                                                                      
    double poissonT = parameter(3);

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
    double lambdaL=youngL*poissonL/(1-poissonL*poissonL);
    double mioL=youngL/(1+poissonL);
    double lambdaT=youngT*poissonT/(1-poissonT*poissonT);                                    // Longitudinal and Transverse Poisson coefficients 
    double mioT=youngT/(1+poissonT);
    
    //---- Anisotropic Correction-------------------------------
    double const PI=3.14159265358979;
    double deltaLam=lambdaL-lambdaT;
    double deltaMio=mioL-mioT;
    double teta[3];
    



    teta[0]=2*PI/3;          // angle between anisotropy vector and D1(dirrection of D is perpendicular to triangle edge
    teta[1]=2*PI/3;          // these angles are the only data we need for considering anisotropy direction.These amounts are
    teta[2]=0;               // just for test.
    //------------------------------------------------------------------------

    //Resting area of the element (using Heron's formula)                                      
    double restingArea=std::sqrt( ( restingLength[0]+restingLength[1]+restingLength[2])*
                                  (-restingLength[0]+restingLength[1]+restingLength[2])*
                                  ( restingLength[0]-restingLength[1]+restingLength[2])*
                                  ( restingLength[0]+restingLength[1]-restingLength[2])  )*0.25;
    
    //Angles of the element ( assuming the order: 0,L0,1,L1,2,L2 )
    std::vector<double> Angle(3);
    Angle[0]=std::asin(2*restingArea/(restingLength[0]*restingLength[2]));
    // can be ommited by cotan(A)=.25*sqrt(4*b*b*c*c/K-1)
    Angle[1]=std::asin(2*restingArea/(restingLength[0]*restingLength[1]));        
    Angle[2]=std::asin(2*restingArea/(restingLength[1]*restingLength[2]));
    
    //Tensile Stiffness
    double tensileStiffness[3];
    double const temp = 1.0/(restingArea*16);                                      
    double cotan[3] = {1.0/std::tan(Angle[0]),1.0/std::tan(Angle[1]),1.0/std::tan(Angle[2])};    
    tensileStiffness[0]=(2*cotan[2]*cotan[2]*(lambdaT+mioT)+mioT)*temp;                     //the force is calculated based on Transverse coefficients
    tensileStiffness[1]=(2*cotan[0]*cotan[0]*(lambdaT+mioT)+mioT)*temp;                     //Longitudinal coefficients are considered in deltaF
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
    
//---- Anisotropic Correction-------------------------------
    
    double  Rcirc2=(0.25*length[0]*length[1]*length[2]/Area)*(0.25*length[0]*length[1]*length[2]/Area);  // square of radius of circumstancing circle
   
    
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

    int k;
    for ( int m=0 ; m<3 ; ++m )
      {
        for ( int coor=0 ; coor<3 ; ++coor ) derIprim1[m][coor]=0;
	for ( int i=0 ; i<3 ; ++i )
	  {            
	    if ((i==0 && m==1)||(i==1 && m==0)) k=2;
	    if ((i==0 && m==2)||(i==2 && m==0)) k=1;
	    if ((i==1 && m==2)||(i==2 && m==1)) k=0; 
	    if (i!=m) DiDm=-0.5*cotan[k]/restingArea;
	    if (i==m) DiDm=0.25*restingLength[i]*restingLength[i] / (restingArea*restingArea);
	    for ( int coor=0 ; coor<3 ; ++coor ) derIprim1[m][coor]=derIprim1[m][coor]+2*DiDm*position[m][coor];
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
		    if ( s!=r ) QrQs=Rcirc2-(length[k]*length[k])*0.5;
		    if ( s==r ) QrQs=Rcirc2;
		    
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
		    
		    for ( int coor=0 ; coor<3 ; ++coor ) derIprim5[p][coor]=derIprim5[p][coor]+2*(DnDr*aDs*aDp+DsDp*aDr*aDn)*QrQs*position[n][coor];
		  }
	      }
	  }
      }		       
    
    double derI1[3][3];             // Invariants and their derivatives
    double derI4[3][3];
    double derI5[3][3];
    
    
    double I1=( Delta[0]*cotan[0]+ Delta[1]*cotan[1]+Delta[2]*cotan[2])/(2*restingArea);
    double I4=0.5*Iprim4-0.5;
    
    for ( int i=0 ; i<3 ; ++i ) for ( int j=0 ; j<3 ; ++j ) derI1[i][j]=0.5*derIprim1[i][j];
    for ( int i=0 ; i<3 ; ++i ) for ( int j=0 ; j<3 ; ++j ) derI4[i][j]=0.5*derIprim4[i][j];
    for ( int i=0 ; i<3 ; ++i ) for ( int j=0 ; j<3 ; ++j ) derI5[i][j]=0.25*derIprim5[i][j]-derI4[i][j];
  
   
    double deltaF[3][3]; 
    
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
	    
	    
