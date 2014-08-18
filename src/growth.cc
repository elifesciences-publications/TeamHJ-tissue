//
// Filename     : growth.cc
// Description  : Classes describing growth updates
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : April 2006
// Revision     : $Id:$
//
#include <vector>
#include <cmath>
#include"growth.h"
#include"baseReaction.h"

namespace WallGrowth {
  Constant::
  Constant(std::vector<double> &paraValue, 
	   std::vector< std::vector<size_t> > 
	   &indValue ) {
    
    //Do some checks on the parameters and variable indeces
    //
    if( paraValue.size()!=2 && paraValue.size() !=3) {
      std::cerr << "WallGrowth::Constant::"
		<< "Constant() "
		<< "Two or three parameters used  k_growth, linearFlag, [L_trunc]"
		<< std::endl;
      exit(EXIT_FAILURE);
    }
    if( indValue.size() != 1 || indValue[0].size() != 1 ) {
      std::cerr << "WallGrowth::Const::"
		<< "Const() "
		<< "One variable index is used for specifying resting length."
		<< std::endl;
      exit(EXIT_FAILURE);
    }
    //Set the variable values
    //
    setId("WallGrowth::Constant");
    setParameter(paraValue);  
    setVariableIndex(indValue);
    
    //Set the parameter identities
    //
    std::vector<std::string> tmp( numParameter() );
    tmp.resize( numParameter() );
    tmp[0] = "k_growth";
    tmp[1] = "linearFlag";
    if (numParameter()>2) {
      tmp[1] = "L_trunc";
    }
    setParameterId( tmp );
  }
  
  void Constant::
  derivs(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs ) {
    
    size_t numWalls = T.numWall();
    size_t lengthIndex = variableIndex(0,0);
    
    for( size_t i=0 ; i<numWalls ; ++i ) {
      double arg = parameter(0);
      if (parameter(1)==1) {//linearFlag
	arg *= wallData[i][lengthIndex];
      }
      if (numParameter()>2) {//truncated at maximal length
	arg *= (1 - wallData[i][lengthIndex]/parameter(2));
      }
      wallDerivs[i][lengthIndex] += arg;
    }
  }

  Stress::
  Stress(std::vector<double> &paraValue, 
	 std::vector< std::vector<size_t> > 
	 &indValue ) {
    
    // Do some checks on the parameters and variable indeces
    //
    if( paraValue.size()!=4 && paraValue.size()!=5 ) {
      std::cerr << "WallGrowth::Stress::"
		<< "Stress() "
		<< "Uses four or five parameters k_growth, stress_threshold "
		<< "stretch_flag (0 for stress (read from wall variable),1 for strain),"
		<< " linear_flag (0 const, 1 prop to wall length), and [L_threshold]." 
		<< std::endl;
      exit(EXIT_FAILURE);
    }
    if( paraValue[2] != 0.0 && paraValue[2] != 1.0 ) {
      std::cerr << "WallGrowth::Stress::"
		<< "Stress() "
		<< "stretch_flag parameter must be 0 (stress used) or " 
		<< "1 (stretch used)." << std::endl;
      exit(EXIT_FAILURE);
    }
    if( paraValue[3] != 0.0 && paraValue[3] != 1.0 ) {
      std::cerr << "WallGrowth::Stress::"
		<< "Stress() "
		<< "linear_flag parameter must be 0 (constant growth) or " 
		<< "1 (length dependent growth)." << std::endl;
      exit(EXIT_FAILURE);
    }
    
    if( (indValue.size()!=1 && indValue.size()!=2) || indValue[0].size() != 1 
	|| (paraValue[2]==0 && (indValue.size()!=2 || !indValue[1].size())) ) {
      std::cerr << "WallGrowth::Stress::"
		<< "Stress() "
		<< "One variable index is used (wall length index) at first "
		<< "level, and stress variable index/indices at second (if strain_flag not set)."
		<< std::endl;
      exit(EXIT_FAILURE);
    }
    //Set the variable values
    //
    setId("WallGrowth::Stress");
    setParameter(paraValue);  
    setVariableIndex(indValue);
    
    //Set the parameter identities
    //
    std::vector<std::string> tmp( numParameter() );
    tmp.resize( numParameter() );
    tmp[0] = "k_growth";
    tmp[1] = "s_threshold";
    tmp[2] = "strain_flag";
    tmp[3] = "linear_flag";
    if (numParameter()>4) {
      tmp[4] = "L_threshold";
    }
    setParameterId( tmp );
  }
  
  void Stress::
  derivs(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs ) {
    
    size_t numWalls = T.numWall();
    size_t lengthIndex = variableIndex(0,0);
    
    for( size_t i=0 ; i<numWalls ; ++i ) {
      size_t v1 = T.wall(i).vertex1()->index();
      size_t v2 = T.wall(i).vertex2()->index();
      double stress=0.0;
      if (!parameter(2)) {//Stress used, read from saved data in the wall
	for (size_t k=0; k<numVariableIndex(1); ++k)
	  stress += wallData[i][variableIndex(1,k)];
      }
      else { //Strain/stretch used
	double distance=0.0;
	for( size_t d=0 ; d<vertexData[v1].size() ; d++ )
	  distance += (vertexData[v1][d]-vertexData[v2][d])*
	    (vertexData[v1][d]-vertexData[v2][d]);
	distance = std::sqrt(distance);
	stress = (distance-wallData[i][lengthIndex]) /
	  wallData[i][lengthIndex];
      }
      if (parameter(1)==0.0 || stress > parameter(1)) {
	double growthRate = parameter(0)*(stress - parameter(1));
	if (parameter(3))
	  growthRate *= wallData[i][lengthIndex];
	if (numParameter()>4) {
	  growthRate *= (1.0 - wallData[i][lengthIndex]/parameter(4));
	}
	wallDerivs[i][lengthIndex] += growthRate;
      }
    }
  }
  
  namespace CenterTriangulation {
    Constant::
    Constant(std::vector<double> &paraValue, 
	     std::vector< std::vector<size_t> > 
	     &indValue ) {
      
      //Do some checks on the parameters and variable indeces
      //
      if( paraValue.size()!=2 && paraValue.size() !=3) {
	std::cerr << "WallGrowth::CenterTriangulation::Constant::"
		  << "Constant() "
		  << "Two or three parameters used  k_growth, linearFlag, [L_trunc]"
		  << std::endl;
	exit(EXIT_FAILURE);
      }
      if( indValue.size() != 1 || indValue[0].size() != 1 ) {
	std::cerr << "WallGrowth::CenterTriangulation::Const::"
		  << "Const() "
		  << "Start of additional Cell variable indices (center(x,y,z) "
		  << "L_1,...,L_n, n=num vertex) is given in first level." 
		  << std::endl;
	exit(EXIT_FAILURE);
      }
      //Set the variable values
      //
      setId("WallGrowth::CenterTriangulation::Constant");
      setParameter(paraValue);  
      setVariableIndex(indValue);
      
      //Set the parameter identities
      //
      std::vector<std::string> tmp( numParameter() );
      tmp.resize( numParameter() );
      tmp[0] = "k_growth";
      tmp[1] = "linearFlag";
      if (numParameter()>2) {
	tmp[1] = "L_trunc";
      }
      setParameterId( tmp );
    }
    
    void Constant::
    derivs(Tissue &T,
	   DataMatrix &cellData,
	   DataMatrix &wallData,
	   DataMatrix &vertexData,
	   DataMatrix &cellDerivs,
	   DataMatrix &wallDerivs,
	   DataMatrix &vertexDerivs ) {
      
      size_t numCells = T.numCell();
      size_t lengthIndex = variableIndex(0,0);
      size_t lengthStartIndex = lengthIndex+3;//assuming 3D
      
      for (size_t i=0; i<T.numCell(); ++i) {
	for (size_t k=0; k<T.cell(i).numVertex(); ++k) {
	  double arg = parameter(0);
	  if (parameter(1)==1) {//linearFlag (prop to length)
	    arg *= cellData[i][k+lengthStartIndex];
	  }
	  if (numParameter()>2) {//truncated at maximal length
	    arg *= (1 - cellData[i][k+lengthStartIndex]/parameter(2));
	  }
	  cellDerivs[i][k+lengthStartIndex] += arg;
	}
      }
    }
    
    Stress::
    Stress(std::vector<double> &paraValue, 
	   std::vector< std::vector<size_t> > 
	   &indValue ) {
      
      //Do some checks on the parameters and variable indeces
      //
      if( paraValue.size()!=4 ) {
	std::cerr << "WallGrowthStresscenterTriangulation::"
		  << "WallGrowthStresscenterTriangulation() "
		  << "Uses four parameters k_growth, stress_threshold "
		  << "stretch_flag and linear_flag (0 const, 1 prop to wall length)" 
		  << std::endl;
	exit(0);
      }
      if( paraValue[2] != 0.0 && paraValue[2] != 1.0 ) {
	std::cerr << "WallGrowthStresscenterTriangulation::"
		  << "WallGrowthStresscenterTriangulation() "
		  << "stretch_flag parameter must be 0 (stress used) or " 
		  << "1 (stretch used)." << std::endl;
	exit(0);
      }
      if( paraValue[2] == 0.0 ) {
	std::cerr << "WallGrowthStresscenterTriangulation::"
		  << "WallGrowthStresscenterTriangulation() "
		  << "stretch_flag parameter must be 1 (stretch used) (not implemented for" 
		  << " stress yet..." 
		  << std::endl;
	exit(0);
      }
      if( paraValue[3] != 0.0 && paraValue[3] != 1.0 ) {
	std::cerr << "WallGrowthStresscenterTriangulation::"
		  << "WallGrowthStresscenterTriangulation() "
		  << "linear_flag parameter must be 0 (constant growth) or " 
		  << "1 (length dependent growth)." << std::endl;
	exit(0);
      }
      
      if( (indValue.size()!=1 && indValue.size()!=2) || indValue[0].size() != 1 
	  || (paraValue[2]==0 && (indValue.size()!=2 || !indValue[1].size())) ) {
	std::cerr << "WallGrowthStresscenterTriangulation::"
		  << "WallGrowthStresscenterTriangulation() "
		  << "Start of additional Cell variable indices (center(x,y,z) "
		  << "L_1,...,L_n, n=num vertex) is given in first level, " 
		  << "and stress variable indices at second (if strain_flag not set)."
		  << std::endl;
	exit(0);
      }
      //Set the variable values
      //
      setId("WallGrowthStresscenterTriangulation");
      setParameter(paraValue);  
      setVariableIndex(indValue);
      
      //Set the parameter identities
      //
      std::vector<std::string> tmp( numParameter() );
      tmp.resize( numParameter() );
      tmp[0] = "k_growth";
      tmp[1] = "s_threshold";
      tmp[2] = "strain_flag";
      tmp[3] = "linear_flag";
      setParameterId( tmp );
    }
    
    void Stress::
    derivs(Tissue &T,
	   DataMatrix &cellData,
	   DataMatrix &wallData,
	   DataMatrix &vertexData,
	   DataMatrix &cellDerivs,
	   DataMatrix &wallDerivs,
	   DataMatrix &vertexDerivs ) {
      
      size_t numCells = T.numCell();
      size_t posStartIndex = variableIndex(0,0);
      size_t lengthStartIndex = posStartIndex+3;
      
      for (size_t i=0; i<numCells; ++i) {
	for (size_t k=0; k<T.cell(i).numVertex(); ++k) {
	  size_t v = T.cell(i).vertex(k)->index();
	  double stress=0.0;
	  if (!parameter(2)) {//Stress used, read from saved data in the wall
	    std::cerr << "WallGrowthStresscenterTriangulation::derivs() " << std::endl
		      << "Strain (and not stress) is the only implemented version sofar."
		      << std::endl;
	    //for (size_t k=0; k<numVariableIndex(1); ++k)
	    //stress += wallData[i][variableIndex(1,k)];
	  }
	  else { //Strain/stretch used
	    double distance=0.0;
	    for( size_t d=0 ; d<vertexData[v].size() ; d++ )
	      distance += (vertexData[v][d]-cellData[i][d+posStartIndex])*
		(vertexData[v][d]-cellData[i][d+posStartIndex]);
	    distance = std::sqrt(distance);
	    stress = (distance-cellData[i][k+lengthStartIndex]) /
	      cellData[i][k+lengthStartIndex];
	  }
	  if (stress > parameter(1)) {
	    double growthRate = parameter(0)*(stress - parameter(1));
	    if (parameter(3))
	      growthRate *= cellData[i][k+lengthStartIndex];
	    cellDerivs[i][k+lengthStartIndex] += growthRate;
	  }
	}
      }
    }


    ///////////////////////////////////strainTRBS begin

    StrainTRBS::
    StrainTRBS(std::vector<double> &paraValue, 
	   std::vector< std::vector<size_t> > 
	   &indValue ) {
      
      //Do some checks on the parameters and variable indeces
      //
      if( paraValue.size()!=2 ) {
	std::cerr << "WallGrowthStrainTRBScenterTriangulation::"
		  << "WallGrowthStrainTRBScenterTriangulation() "
		  << "Uses two parameters k_growth, strain_threshold "
                  << std::endl;
	exit(0);
      }
        
      if( indValue.size() !=3 || 
          indValue[0].size() !=1 || 
          indValue[1].size() !=1 || 
          indValue[2].size() !=3  ) {
	std::cerr << "WallGrowthStrainTRBScenterTriangulation::"
		  << "WallGrowthStrainTRBScenterTriangulation() "
                  << "wall length index is given in first level,"
		  << "Start of additional Cell variable indices (center(x,y,z) "
		  << "L_1,...,L_n, n=num vertex) is given in second level, " 
		  << "and strain1_value_strore_index and strain2_value_strore_index and "
                  <<" strain_vector_start_store_index(3 components) at second level."
		  << std::endl;
	exit(0);
      }
      //Set the variable values
      //
      setId("WallGrowthStrainTRBScenterTriangulation");
      setParameter(paraValue);  
      setVariableIndex(indValue);
      
      //Set the parameter identities
      //
      std::vector<std::string> tmp( numParameter() );
      tmp.resize( numParameter() );
      tmp[0] = "k_growth";
      tmp[1] = "s_threshold";
      setParameterId( tmp );
    }
    
    void StrainTRBS::
    derivs(Tissue &T,
	   DataMatrix &cellData,
	   DataMatrix &wallData,
	   DataMatrix &vertexData,
	   DataMatrix &cellDerivs,
	   DataMatrix &wallDerivs,
	   DataMatrix &vertexDerivs ) {}

    void StrainTRBS::
    update(Tissue &T,
	   DataMatrix &cellData,
	   DataMatrix &wallData,
	   DataMatrix &vertexData,
           double h ) {

      size_t dimension = 3;
      size_t numCells = T.numCell();
      size_t numWalls = T.numWall();
      size_t wallLengthIndex= variableIndex(0,0);
      size_t comIndex = variableIndex(1,0);
      size_t lengthInternalIndex = comIndex+dimension;
      size_t strainValIndex1 = variableIndex(2,0);
      size_t strainValIndex2 = variableIndex(2,1);
      size_t strainVecIndex = variableIndex(2,2);
      double strainThreshold=parameter(1);

      std::vector<std::vector<double> > mainWalls(numWalls);
      std::vector<std::vector<std::vector<double> > > internalWalls(numCells);
      
      for (size_t wallIndex=0 ; wallIndex<numWalls ; ++wallIndex)
        mainWalls[wallIndex].resize(2);
      for (size_t cellIndex=0 ; cellIndex<numCells ; ++cellIndex){
        size_t numCellWalls = T.cell(cellIndex).numWall();
        internalWalls[cellIndex].resize(numCellWalls);
        for (size_t cellWallIndex=0 ; cellWallIndex<numCellWalls ; ++cellWallIndex)
          internalWalls[cellIndex][cellWallIndex].resize(2);
      }
      
      for (size_t cellIndex=0 ; cellIndex<numCells ; ++cellIndex) {
        size_t numCellWalls = T.cell(cellIndex).numWall(); 
        // internalWallLenth[cellIndex].resize(numCellWalls);
        // externalWallLength[cellIndex].resize(numCellWalls);
      

        for (size_t wallIndex=0; wallIndex<numCellWalls; ++wallIndex) { 
          size_t wallIndexPlusOneMod = (wallIndex+1)%numCellWalls;
          //size_t v1 = com;
          size_t v2 = T.cell(cellIndex).vertex(wallIndex)->index();
          size_t v3 = T.cell(cellIndex).vertex(wallIndexPlusOneMod)->index();
          //size_t w1 = internal wallIndex
          size_t w2 = T.cell(cellIndex).wall(wallIndex)->index();
          //size_t w3 = internal wallIndex+1

          DataMatrix position(3,vertexData[v2]);
          for (size_t d=0; d<dimension; ++d)
            position[0][d] = cellData[cellIndex][comIndex+d]; // com position
          //position[1] = vertexData[v2]; // given by initiation
          position[2] = vertexData[v3];
          
          std::vector<double> restingLength(3);
          restingLength[0] = cellData[cellIndex][lengthInternalIndex + wallIndex];
          restingLength[1] = wallData[w2][wallLengthIndex];
          restingLength[2] = cellData[cellIndex][lengthInternalIndex + wallIndexPlusOneMod];
         
          double restingArea=std::sqrt( ( restingLength[0]+restingLength[1]+restingLength[2])*
                                        (-restingLength[0]+restingLength[1]+restingLength[2])*
                                        ( restingLength[0]-restingLength[1]+restingLength[2])*
                                        ( restingLength[0]+restingLength[1]-restingLength[2])  )*0.25;
          
          std::vector<double> length(3);
          length[0] = std::sqrt( (position[0][0]-position[1][0])*(position[0][0]-position[1][0]) +
                                 (position[0][1]-position[1][1])*(position[0][1]-position[1][1]) +
                                 (position[0][2]-position[1][2])*(position[0][2]-position[1][2]) );
          
          length[1] = T.wall(w2).lengthFromVertexPosition(vertexData);
          
          length[2] = std::sqrt( (position[0][0]-position[2][0])*(position[0][0]-position[2][0]) +
                                 (position[0][1]-position[2][1])*(position[0][1]-position[2][1]) +
                                 (position[0][2]-position[2][2])*(position[0][2]-position[2][2]) );
          
          std::vector<double> strainVec(3);
          //double strainVecL;

          for (size_t d=0; d< dimension; ++d)
            strainVec[d]=cellData[cellIndex][strainVecIndex+d];
          // strainVecL=std::sqrt(strainVec[0]*strainVec[0]+
          //                      strainVec[1]*strainVec[1]+
          //                      strainVec[2]*strainVec[2]);
          
          //Angles of the element ( assuming the order: 0,L0,1,L1,2,L2 clockwise )
          // std::vector<double> Angle(3);
          // // can be ommited by cotan(A)=.25*sqrt(4*b*b*c*c/K-1)
          // Angle[0]=std::acos(  (restingLength[0]*restingLength[0]+
          //                       restingLength[2]*restingLength[2]-
          //                       restingLength[1]*restingLength[1])/
          //                      (restingLength[0]*restingLength[2]*2)    );
          // Angle[1]=std::acos(  (restingLength[0]*restingLength[0]+
          //                       restingLength[1]*restingLength[1]-
          //                       restingLength[2]*restingLength[2])/
          //                      (restingLength[0]*restingLength[1]*2)    );
          // Angle[2]=std::acos(  (restingLength[1]*restingLength[1]+
          //                       restingLength[2]*restingLength[2]-
          //                       restingLength[0]*restingLength[0])/
          //                      (restingLength[1]*restingLength[2]*2)    );
    
          //Current shape local coordinate of the element  (counterclockwise ordering of nodes/edges)
          double CurrentAngle1=std::acos(  (length[0]*length[0]+
                                            length[1]*length[1]-
                                            length[2]*length[2])/
                                           (length[0]*length[1]*2)    );
          
          double Qa=std::cos(CurrentAngle1)*length[0];
          double Qc=std::sin(CurrentAngle1)*length[0];
          double Qb=length[1];
          

          double RestingAngle1=std::acos(  (restingLength[0]*restingLength[0]+
                                            restingLength[1]*restingLength[1]-
                                            restingLength[2]*restingLength[2])/
                                           (restingLength[0]*restingLength[1]*2)    );
          
          double Pa=std::cos(RestingAngle1)*restingLength[0];
          double Pc=std::sin(RestingAngle1)*restingLength[0];
          double Pb=restingLength[1];
          
          // shape vector matrix in resting shape in local coordinate system  = 
          // inverse of coordinate matrix ( only first two elements i.e. ShapeVectorResting[3][2] )      
          double ShapeVectorResting[3][3]={ {  0   ,       1/Pc      , 0 }, 
                                            {-1/Pb , (Pa-Pb)/(Pb*Pc) , 1 },       
                                            { 1/Pb ,     -Pa/(Pb*Pc) , 0 }  };
          
          double positionLocal[3][2]={ {Qa , Qc}, 
                                       {0  , 0 },  
                                       {Qb , 0 }  };
          
          double DeformGrad[2][2]={{0,0},{0,0}}; // F= Qi x Di
          for ( int ii=0 ; ii<3 ; ++ii ) {
            DeformGrad[0][0]=DeformGrad[0][0]+positionLocal[ii][0]*ShapeVectorResting[ii][0];
            DeformGrad[1][0]=DeformGrad[1][0]+positionLocal[ii][1]*ShapeVectorResting[ii][0];
            DeformGrad[0][1]=DeformGrad[0][1]+positionLocal[ii][0]*ShapeVectorResting[ii][1];
            DeformGrad[1][1]=DeformGrad[1][1]+positionLocal[ii][1]*ShapeVectorResting[ii][1];
          }

          double Egreen[2][2];//E=0.5(C-I)
          Egreen[0][0]=0.5*(DeformGrad[0][0]*DeformGrad[0][0]+DeformGrad[1][0]*DeformGrad[1][0]-1);
          Egreen[1][0]=0.5*(DeformGrad[0][1]*DeformGrad[0][0]+DeformGrad[1][1]*DeformGrad[1][0]);
          Egreen[0][1]=0.5*(DeformGrad[0][0]*DeformGrad[0][1]+DeformGrad[1][0]*DeformGrad[1][1]);
          Egreen[1][1]=0.5*(DeformGrad[0][1]*DeformGrad[0][1]+DeformGrad[1][1]*DeformGrad[1][1]-1);
          
          
          double det=Egreen[0][0]*Egreen[1][1]-Egreen[0][1]*Egreen[1][0];
          double tr=Egreen[0][0]+Egreen[1][1];

          double strainValue1=tr/2+std::sqrt(((tr*tr)/4)-det);
          double strainValue2=tr/2-std::sqrt(((tr*tr)/4)-det);

          std::cerr<<"strain values"<<strainValue1<<"  "<<strainValue2<<std::endl;
          double strainRestLocal1[2]={0,0};
          double strainRestLocal2[2]={0,0};
          if (Egreen[1][0]!=0){
            strainRestLocal1[0]=strainValue1-Egreen[1][1];
            strainRestLocal1[1]=Egreen[1][0];
            strainRestLocal2[0]=strainValue2-Egreen[1][1];
            strainRestLocal2[1]=Egreen[1][0];
          }
          else if (Egreen[0][1]!=0){
            strainRestLocal1[0]=Egreen[0][1];
            strainRestLocal1[1]=strainValue1-Egreen[0][0];
            strainRestLocal2[0]=Egreen[0][1];
            strainRestLocal2[1]=strainValue2-Egreen[0][0];
          }
          else {
            strainRestLocal1[0]=1;
            strainRestLocal1[1]=0;
            strainRestLocal2[0]=0;
            strainRestLocal2[1]=1;
          }
 
          double tempAn=std::sqrt(strainRestLocal1[0]*strainRestLocal1[0]+
                                  strainRestLocal1[1]*strainRestLocal1[1]);
          if (tempAn !=0 ){
            strainRestLocal1[0]/=tempAn;
            strainRestLocal1[1]/=tempAn;
          }
          tempAn=std::sqrt(strainRestLocal2[0]*strainRestLocal2[0]+
                           strainRestLocal2[1]*strainRestLocal2[1]);
          if (tempAn !=0 ){
            strainRestLocal2[0]/=tempAn;
            strainRestLocal2[1]/=tempAn;
          }
          
          if (strainValue2>strainValue1){
            double temp=strainValue1;
            strainValue1=strainValue2;
            strainValue2=temp;
            temp=strainRestLocal1[0];
            strainRestLocal1[0]=strainRestLocal2[0];
            strainRestLocal2[0]=temp;        
            temp=strainRestLocal1[1];
            strainRestLocal1[1]=strainRestLocal2[1];
            strainRestLocal2[1]=temp;          
          }

          std::vector<std::vector<double> > edgeRestLocal(3);
          
          for (size_t d=0; d< 3; ++d)
            edgeRestLocal[d].resize(3);
          
          
          edgeRestLocal[0][0]= -Pa;   //positionLocal[][0]-positionLocal[][0];
          edgeRestLocal[0][1]= -Pc;   //positionLocal[][1]-positionLocal[][1];
          //edgeRestLocal[0][2]= positionLocal[][2]-positionLocal[][2];
          
          edgeRestLocal[1][0]= Pb;    //positionLocal[][0]-positionLocal[][0];
          edgeRestLocal[1][1]= 0;     //positionLocal[][1]-positionLocal[][1];
          //edgeRestLocal[1][2]= positionLocal[][2]-positionLocal[][2];
          
          edgeRestLocal[2][0]= Pa-Pb; //positionLocal[][0]-positionLocal[][0];
          edgeRestLocal[2][1]= Pc;    //positionLocal[][1]-positionLocal[][1];
          //edgeRestLocal[2][2]= positionLocal[][2]-positionLocal[][2];
          

          std:: vector<double> cosTet(3);
          std:: vector<double> sinTet(3);
          for (size_t j=0; j< 3; ++j){
            cosTet[j]=std::abs((strainRestLocal1[0]*edgeRestLocal[j][0]+
                                strainRestLocal1[1]*edgeRestLocal[j][1])/
                               restingLength[j]);
            sinTet[j]=std::sqrt(1-cosTet[j]*cosTet[j]);
    
      }

          //std::cerr<<"  "<< cosTet[0]<<"  "<<cosTet[1]<<"  "<<cosTet[2]<<std::endl;
          
          // double strainValue1=cellData[cellIndex][strainValIndex1];
          // double strainValue2=cellData[cellIndex][strainValIndex2];
      
          std::vector<std::vector<double> > restingComp(3);
          for (size_t j=0; j< 3; ++j)
            restingComp[j].resize(2);
          
          for (size_t j=0; j< 3; ++j){
            restingComp[j][0]=restingLength[j]*cosTet[j];
            restingComp[j][1]=restingLength[j]*sinTet[j];
          }
          std::cerr<<" 1 "<<std::endl;
          if (strainValue1>strainThreshold && strainValue2<strainThreshold){
            std::cerr<<" 2 "<<std::endl;            
            for (size_t j=0; j< 3; ++j)
              restingComp[j][0]+=restingComp[j][0]*parameter(0)*(strainValue1-strainThreshold);
          }
          
          if (strainValue1>strainThreshold && strainValue2>strainThreshold){
             std::cerr<<" 3 "<<std::endl;
            for (size_t j=0; j< 3; ++j){
              restingComp[j][0]+=restingComp[j][0]*parameter(0)*(strainValue1-strainThreshold);
              restingComp[j][1]+=restingComp[j][1]*parameter(0)*(strainValue2-strainThreshold);
            }
            
          }
          
    
          double internalTemp=std::sqrt(restingComp[0][0]*restingComp[0][0]+
                                        restingComp[0][1]*restingComp[0][1]);
          
          double externalTemp=std::sqrt(restingComp[1][0]*restingComp[1][0]+
                                        restingComp[1][1]*restingComp[1][1]);
          
          double internalTempPlusOne=std::sqrt(restingComp[2][0]*restingComp[2][0]+
                                               restingComp[2][1]*restingComp[2][1]);

          std::cerr<<" edges "<<internalTemp<<"  "<<externalTemp<<"  "<<internalTempPlusOne<<std::endl;
          
          // WITH AREA AVERAGING

          size_t wallGlobalInd= T.cell(cellIndex).wall(wallIndex) ->index();
          std::cerr<<" normalization area factor before "<<mainWalls[wallGlobalInd][0]<<std::endl;
          if (mainWalls[wallGlobalInd][0]==0){
            mainWalls[wallGlobalInd][0]=restingArea;
            mainWalls[wallGlobalInd][1]=restingArea*externalTemp;
            std::cerr<<" normalization area factor middle "<<mainWalls[wallGlobalInd][0]<<std::endl;
          }
          else if (mainWalls[wallGlobalInd][0]!=0){
            mainWalls[wallGlobalInd][0]+=restingArea;
            mainWalls[wallGlobalInd][1]=
              (mainWalls[wallGlobalInd][1]+restingArea*externalTemp);
          }
          std::cerr<<" normalization area factor after "<<mainWalls[wallGlobalInd][0]<<std::endl;
          std::cerr<<" main wall "<<mainWalls[wallGlobalInd][1]<<std::endl;
          //wallIndexPlusOneMod
            
          if (internalWalls[cellIndex][wallIndex][0]==0){
            internalWalls[cellIndex][wallIndex][0]=restingArea;
            internalWalls[cellIndex][wallIndex][1]=restingArea*internalTemp;
          }
          else if (internalWalls[cellIndex][wallIndex][0]!=0){
            internalWalls[cellIndex][wallIndex][0]+=restingArea;
            internalWalls[cellIndex][wallIndex][1]=
              (internalWalls[cellIndex][wallIndex][1]+restingArea*internalTemp);
          }         
          
          if (internalWalls[cellIndex][wallIndexPlusOneMod][0]==0){
            internalWalls[cellIndex][wallIndexPlusOneMod][0]=restingArea;
            internalWalls[cellIndex][wallIndexPlusOneMod][1]=restingArea*internalTempPlusOne;
          }
          else if (internalWalls[cellIndex][wallIndexPlusOneMod][0]!=0){
            internalWalls[cellIndex][wallIndexPlusOneMod][0]+=restingArea;
            internalWalls[cellIndex][wallIndexPlusOneMod][1]=
              (internalWalls[cellIndex][wallIndexPlusOneMod][1]+restingArea*internalTempPlusOne);
              
          }         
          

        } // walls
        
        // updating wall length
                
      } // cells
     
      for (size_t cellIndex=0; cellIndex< T.numCell(); cellIndex++)
        for (size_t wallIndex=0; wallIndex< T.cell(cellIndex).numWall(); wallIndex++){
          cellData[cellIndex][lengthInternalIndex + wallIndex]= 
            internalWalls[cellIndex][wallIndex][1]
            /internalWalls[cellIndex][wallIndex][0]; 
          size_t wallGlobalInd=T.cell(cellIndex).wall(wallIndex)->index();
          wallData[wallGlobalInd][wallLengthIndex]=
            mainWalls[wallGlobalInd][1]/mainWalls[wallGlobalInd][0];
        }
      
    }

    // void StrainTRBS::
    // update(Tissue &T,
    //        DataMatrix &cellData,
    //        DataMatrix &wallData,
    //        DataMatrix &vertexData,
    //        double h ) {

    //   size_t dimension = 3;
    //   size_t numCells = T.numCell();
    //   size_t wallLengthIndex= variableIndex(0,0);
    //   size_t comIndex = variableIndex(1,0);
    //   size_t lengthInternalIndex = comIndex+dimension;
    //   size_t strainValIndex1 = variableIndex(2,0);
    //   size_t strainValIndex2 = variableIndex(2,1);
    //   size_t strainVecIndex = variableIndex(2,2);
    //   double strainThreshold=parameter(1);

    //   std::vector<std::vector<double> > internalWallLenth(numCells),externalWallLength(numCells);

    //   for (size_t cellIndex=0 ; cellIndex<numCells ; ++cellIndex) {
    //     size_t numWalls = T.cell(cellIndex).numWall(); 
    //     internalWallLenth[cellIndex].resize(numWalls);
    //     externalWallLength[cellIndex].resize(numWalls);
    //     for (size_t wallIndex=0; wallIndex<numWalls; ++wallIndex) { 
    //       size_t wallIndexPlusOneMod = (wallIndex+1)%numWalls;
    //       //size_t v1 = com;
    //       size_t v2 = T.cell(cellIndex).vertex(wallIndex)->index();
    //       size_t v3 = T.cell(cellIndex).vertex(wallIndexPlusOneMod)->index();
    //       //size_t w1 = internal wallIndex
    //       size_t w2 = T.cell(cellIndex).wall(wallIndex)->index();
    //       //size_t w3 = internal wallIndex+1

    //       DataMatrix position(3,vertexData[v2]);
    //       for (size_t d=0; d<dimension; ++d)
    //         position[0][d] = cellData[cellIndex][comIndex+d]; // com position
    //       //position[1] = vertexData[v2]; // given by initiation
    //       position[2] = vertexData[v3];
          
    //       std::vector<double> restingLength(3);
    //       restingLength[0] = cellData[cellIndex][lengthInternalIndex + wallIndex];
    //       restingLength[1] = wallData[w2][wallLengthIndex];
    //       restingLength[2] = cellData[cellIndex][lengthInternalIndex + wallIndexPlusOneMod];
         
    //       std::vector<double> length(3);
    //       length[0] = std::sqrt( (position[0][0]-position[1][0])*(position[0][0]-position[1][0]) +
    //                              (position[0][1]-position[1][1])*(position[0][1]-position[1][1]) +
    //                              (position[0][2]-position[1][2])*(position[0][2]-position[1][2]) );
          
    //       length[1] = T.wall(w2).lengthFromVertexPosition(vertexData);
          
    //       length[2] = std::sqrt( (position[0][0]-position[2][0])*(position[0][0]-position[2][0]) +
    //                              (position[0][1]-position[2][1])*(position[0][1]-position[2][1]) +
    //                              (position[0][2]-position[2][2])*(position[0][2]-position[2][2]) );
          
    //       std::vector<double> strainVec(3);
    //       //double strainVecL;

    //       for (size_t d=0; d< dimension; ++d)
    //         strainVec[d]=cellData[cellIndex][strainVecIndex+d];
    //       // strainVecL=std::sqrt(strainVec[0]*strainVec[0]+
    //       //                      strainVec[1]*strainVec[1]+
    //       //                      strainVec[2]*strainVec[2]);
          


    //       //Angles of the element ( assuming the order: 0,L0,1,L1,2,L2 clockwise )
    //       // std::vector<double> Angle(3);
    //       // // can be ommited by cotan(A)=.25*sqrt(4*b*b*c*c/K-1)
    //       // Angle[0]=std::acos(  (restingLength[0]*restingLength[0]+
    //       //                       restingLength[2]*restingLength[2]-
    //       //                       restingLength[1]*restingLength[1])/
    //       //                      (restingLength[0]*restingLength[2]*2)    );
    //       // Angle[1]=std::acos(  (restingLength[0]*restingLength[0]+
    //       //                       restingLength[1]*restingLength[1]-
    //       //                       restingLength[2]*restingLength[2])/
    //       //                      (restingLength[0]*restingLength[1]*2)    );
    //       // Angle[2]=std::acos(  (restingLength[1]*restingLength[1]+
    //       //                       restingLength[2]*restingLength[2]-
    //       //                       restingLength[0]*restingLength[0])/
    //       //                      (restingLength[1]*restingLength[2]*2)    );
    
    //       //Current shape local coordinate of the element  (counterclockwise ordering of nodes/edges)
    //       double CurrentAngle1=std::acos(  (length[0]*length[0]+
    //                                         length[1]*length[1]-
    //                                         length[2]*length[2])/
    //                                        (length[0]*length[1]*2)    );
          
    //       double Qa=std::cos(CurrentAngle1)*length[0];
    //       double Qc=std::sin(CurrentAngle1)*length[0];
    //       double Qb=length[1];
          

    //       double RestingAngle1=std::acos(  (restingLength[0]*restingLength[0]+
    //                                         restingLength[1]*restingLength[1]-
    //                                         restingLength[2]*restingLength[2])/
    //                                        (restingLength[0]*restingLength[1]*2)    );
          
          

    //       double Pa=std::cos(RestingAngle1)*restingLength[0];
    //       double Pc=std::sin(RestingAngle1)*restingLength[0];
    //       double Pb=restingLength[1];
          
    //       // shape vector matrix in resting shape in local coordinate system  = 
    //       // inverse of coordinate matrix ( only first two elements i.e. ShapeVectorResting[3][2] )      
    //       double ShapeVectorResting[3][3]={ {  0   ,       1/Pc      , 0 }, 
    //                                         {-1/Pb , (Pa-Pb)/(Pb*Pc) , 1 },       
    //                                         { 1/Pb ,     -Pa/(Pb*Pc) , 0 }  };
          
          

    //       // Rotation Matrix for changing coordinate systems for both Local to 
    //       // Global( Strain Tensor) and Global to Local( Aniso Vector in the current shape)
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
          
    //       // rotating the growth vector from global coordinate system to the local in the current shape
    //       double strainVecLocal[3];
    //       strainVecLocal[0]=
    //         rotation[0][0]*strainVec[0]+
    //         rotation[1][0]*strainVec[1]+
    //         rotation[2][0]*strainVec[2];
    //       strainVecLocal[1]=
    //         rotation[0][1]*strainVec[0]+
    //         rotation[1][1]*strainVec[1]+
    //         rotation[2][1]*strainVec[2];
    //       strainVecLocal[2]=
    //         rotation[0][2]*strainVec[0]+
    //         rotation[1][2]*strainVec[1]+
    //         rotation[2][2]*strainVec[2];
          
    //       double positionLocal[3][2]={ {Qa , Qc}, 
    //                                    {0  , 0 },  
    //                                    {Qb , 0 }  };
          
    //       double DeformGrad[2][2]={{0,0},{0,0}}; // F= Qi x Di
    //       for ( int ii=0 ; ii<3 ; ++ii ) {
    //         DeformGrad[0][0]=DeformGrad[0][0]+positionLocal[ii][0]*ShapeVectorResting[ii][0];
    //         DeformGrad[1][0]=DeformGrad[1][0]+positionLocal[ii][1]*ShapeVectorResting[ii][0];
    //         DeformGrad[0][1]=DeformGrad[0][1]+positionLocal[ii][0]*ShapeVectorResting[ii][1];
    //         DeformGrad[1][1]=DeformGrad[1][1]+positionLocal[ii][1]*ShapeVectorResting[ii][1];
    //       }

          
          
    //       double strainRestLocal[3]={0,0,0};
          
    //       strainRestLocal[0]=DeformGrad[0][0]*strainVecLocal[0]+DeformGrad[1][0]*strainVecLocal[1];
    //       strainRestLocal[1]=DeformGrad[0][1]*strainVecLocal[0]+DeformGrad[1][1]*strainVecLocal[1];
    //       strainRestLocal[2]=strainVecLocal[2];
    //       double tempAn=std::sqrt(strainRestLocal[0]*strainRestLocal[0]+
    //                               strainRestLocal[1]*strainRestLocal[1]+
    //                               strainRestLocal[2]*strainRestLocal[2]);
    //       strainRestLocal[0]/=tempAn;
    //       strainRestLocal[1]/=tempAn;
    //       strainRestLocal[2]/=tempAn;
          


    //       std::vector<std::vector<double> > edgeRestLocal(3);
          
    //       for (size_t d=0; d< 3; ++d)
    //         edgeRestLocal[d].resize(3);
          
          
    //       edgeRestLocal[0][0]= -Pa;   //positionLocal[][0]-positionLocal[][0];
    //       edgeRestLocal[0][1]= -Pc;   //positionLocal[][1]-positionLocal[][1];
    //       //edgeRestLocal[0][2]= positionLocal[][2]-positionLocal[][2];
          
    //       edgeRestLocal[1][0]= Pb;    //positionLocal[][0]-positionLocal[][0];
    //       edgeRestLocal[1][1]= 0;     //positionLocal[][1]-positionLocal[][1];
    //       //edgeRestLocal[1][2]= positionLocal[][2]-positionLocal[][2];
          
    //       edgeRestLocal[2][0]= Pa-Pb; //positionLocal[][0]-positionLocal[][0];
    //       edgeRestLocal[2][1]= Pc;    //positionLocal[][1]-positionLocal[][1];
    //       //edgeRestLocal[2][2]= positionLocal[][2]-positionLocal[][2];
          

    //       std:: vector<double> cosTet(3);
    //       std:: vector<double> sinTet(3);
    //       for (size_t j=0; j< 3; ++j){
    //         cosTet[j]=std::abs((strainRestLocal[0]*edgeRestLocal[j][0]+
    //                             strainRestLocal[1]*edgeRestLocal[j][1])/
    //                            restingLength[j]);
    //         sinTet[j]=std::sqrt(1-cosTet[j]*cosTet[j]);
    
    //   }

    //       //std::cerr<<"  "<< cosTet[0]<<"  "<<cosTet[1]<<"  "<<cosTet[2]<<"  "<<
          
          
    //       double strainValue1=cellData[cellIndex][strainValIndex1];
    //       double strainValue2=cellData[cellIndex][strainValIndex2];
      
    //       std::vector<std::vector<double> > restingComp(3);
    //       for (size_t j=0; j< 3; ++j)
    //         restingComp[j].resize(2);
          
    //       for (size_t j=0; j< 3; ++j){
    //         restingComp[j][0]=restingLength[j]*cosTet[j];
    //         restingComp[j][1]=restingLength[j]*sinTet[j];
    //       }
          

         

    //       //  if (strainValue1>strainThreshold ){
    //       //   for (size_t j=0; j< 3; ++j)
    
    //       //     restingComp[j][0]+=restingComp[j][0]*parameter(0)*(strainValue1-strainThreshold);
            
            
    //       // }


    //       if (strainValue1>strainThreshold && strainValue2<strainThreshold){
    //         for (size_t j=0; j< 3; ++j)
    
    //           restingComp[j][0]+=restingComp[j][0]*parameter(0)*(strainValue1-strainThreshold);
            
            
    //       }
    //       if (strainValue1>strainThreshold && strainValue2>strainThreshold){
    //         for (size_t j=0; j< 3; ++j){
    //           restingComp[j][0]+=restingComp[j][0]*parameter(0)*(strainValue1-strainThreshold);
    //           restingComp[j][1]+=restingComp[j][1]*parameter(0)*(strainValue2-strainThreshold);
    //         }
            
    //       }
          
    
          
    //       internalWallLenth[cellIndex][wallIndex]=std::sqrt(restingComp[0][0]*restingComp[0][0]+
    //                                              restingComp[0][1]*restingComp[0][1]);
    //       externalWallLength[cellIndex][wallIndex]=std::sqrt(restingComp[1][0]*restingComp[1][0]+
    //                                               restingComp[1][1]*restingComp[1][1]);
          
    
    //     } // walls

    //     // updating wall length
                
    //   } // cells
     
    //   for (size_t cellIndex=0; cellIndex< T.numCell(); ++cellIndex)
    //     for (size_t wallIndex=0; wallIndex< T.cell(cellIndex).numWall(); ++wallIndex){
    //       cellData[cellIndex][lengthInternalIndex + wallIndex]= internalWallLenth[cellIndex][wallIndex]; 
    //       wallData[T.cell(cellIndex).wall(wallIndex)->index()][wallLengthIndex]=externalWallLength[cellIndex][wallIndex];             
    //     }
      
    // }
    ///////////////////////////////////strainTRBS end
    
  }
  
  StressSpatial::
  StressSpatial(std::vector<double> &paraValue, 
		std::vector< std::vector<size_t> > 
		&indValue ) {
    
    // Do some checks on the parameters and variable indeces
    //
    if( paraValue.size()!=6 ) {
      std::cerr << "WallGrowth::StressSpatial::"
		<< "StressSpatial() "
		<< "Uses six parameters k_growth, stress(stretch)_threshold "
		<< "K_hill n_Hill "
		<< "stretch_flag and linear_flag" << std::endl;
      exit(0);
    }
    if( paraValue[4] != 0.0 && paraValue[4] != 1.0 ) {
      std::cerr << "WallGrowth::StressSpatial::"
		<< "StressSpatial() "
		<< "stretch_flag parameter must be 0 (stress used) or " 
		<< "1 (stretch used)." << std::endl;
      exit(0);
    }
    if( paraValue[5] != 0.0 && paraValue[5] != 1.0 ) {
      std::cerr << "WallGrowth::StressSpatial::"
		<< "StressSpatial() "
		<< "linear_flag parameter must be 0 (constant growth) or " 
		<< "1 (length dependent growth)." << std::endl;
      exit(0);
    }
    
    if( indValue.size() != 2 || indValue[0].size() != 2 ) {
      std::cerr << "WallGrowth::StressSpatial::"
		<< "StressSpatial() "
		<< "Two variable index is used (wall length,spatial coordinate) at first "
		<< "level, and force variable index at second."
		<< std::endl;
      exit(0);
    }
    // Set the variable values
    //
    setId("WallGrowth::StressSpatial");
    setParameter(paraValue);  
    setVariableIndex(indValue);
    Kpow_=std::pow(paraValue[2],paraValue[3]);
    
    // Set the parameter identities
    //
    std::vector<std::string> tmp( numParameter() );
    tmp.resize( numParameter() );
    tmp[0] = "k_growth";
    tmp[1] = "stress_threshold";
    tmp[2] = "K_Hill";
    tmp[3] = "n_Hill";
    tmp[4] = "stretch_flag";
    tmp[5] = "linear_flag";
    setParameterId( tmp );
  }
  
  void StressSpatial::
  derivs(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs ) {
    
    size_t numWalls = T.numWall();
    size_t lengthIndex = variableIndex(0,0);
    size_t dimension = vertexData[0].size();
    
    // Prepare spatial factor
    size_t sI=variableIndex(0,1);
    assert (sI<vertexData[0].size());
    size_t numVertices = vertexData.size();
    double sMax= vertexData[0][sI];
    size_t maxI=0;
    for (size_t i=1; i<numVertices; ++i)
      if (vertexData[i][sI]>sMax) {
	sMax=vertexData[i][sI];
	maxI=i;
      }
    std::vector<double> maxPos(dimension);
    for (size_t d=0; d<dimension; ++d)
      maxPos[d] = vertexData[maxI][d];
    
    for( size_t i=0 ; i<numWalls ; ++i ) {
      size_t v1 = T.wall(i).vertex1()->index();
      size_t v2 = T.wall(i).vertex2()->index();
      double stress=0.0;
      if (!parameter(4)) {
	for (size_t k=0; k<numVariableIndex(1); ++k)
	  stress += wallData[i][variableIndex(1,k)];
      }
      else {
	double distance=0.0;
	for( size_t d=0 ; d<dimension ; ++d )
	  distance += (vertexData[v1][d]-vertexData[v2][d])*
	    (vertexData[v1][d]-vertexData[v2][d]);
	distance = std::sqrt(distance);
	stress = (distance-wallData[i][lengthIndex]) /
	  wallData[i][lengthIndex];
      }
      if (stress > parameter(1)) {
	// Calculate spatial factor
	double maxDistance = 0.0;
	for (size_t d=0; d<dimension; ++d) {
	  double pos = 0.5*(vertexData[v1][d]+vertexData[v2][d]);
	  maxDistance += (maxPos[d]-pos)*(maxPos[d]-pos);
	}
	maxDistance = std::sqrt(maxDistance);
	double spatialFactor = Kpow_/(Kpow_+std::pow(maxDistance,parameter(3)));
	
	double growthRate = parameter(0)*(stress - parameter(1))*spatialFactor;
		
	if (parameter(5))
	  growthRate *= wallData[i][lengthIndex];
	wallDerivs[i][lengthIndex] += growthRate;
      }
    }
  }
  
  StressSpatialSingle::
  StressSpatialSingle(std::vector<double> &paraValue, 
		      std::vector< std::vector<size_t> > 
		      &indValue ) {
    
    // Do some checks on the parameters and variable indeces
    //
    if( paraValue.size()!=6 ) {
      std::cerr << "WallGrowth::StressSpatialSingle::"
		<< "StressSpatialSingle() "
		<< "Uses six parameters k_growth, stress(stretch)_threshold "
		<< "K_hill n_Hill "
		<< "stretch_flag and linear_flag" << std::endl;
      exit(0);
    }
    if( paraValue[4] != 0.0 && paraValue[4] != 1.0 ) {
      std::cerr << "WallGrowth::StressSpatialSingle::"
		<< "StressSpatialSingle() "
		<< "stretch_flag parameter must be 0 (stress used) or " 
		<< "1 (stretch used)." << std::endl;
      exit(0);
    }
    if( paraValue[5] != 0.0 && paraValue[5] != 1.0 ) {
      std::cerr << "WallGrowth::StressSpatialSingle::"
		<< "StressSpatialSingle() "
		<< "linear_flag parameter must be 0 (constant growth) or " 
		<< "1 (length dependent growth)." << std::endl;
      exit(0);
    }
    
    if( indValue.size() != 2 || indValue[0].size() != 2 ) {
      std::cerr << "WallGrowth::StressSpatialSingle::"
		<< "StressSpatialSingle() "
		<< "Two variable index is used (wall length,spatial coordinate) at first "
		<< "level, and force variable index at second."
		<< std::endl;
      exit(0);
    }
    // Set the variable values
    //
    setId("WallGrowth::StressSpatialSingle");
    setParameter(paraValue);  
    setVariableIndex(indValue);
    Kpow_=std::pow(paraValue[2],paraValue[3]);
    
    // Set the parameter identities
    //
    std::vector<std::string> tmp( numParameter() );
    tmp.resize( numParameter() );
    tmp[0] = "k_growth";
    tmp[1] = "stress_threshold";
    tmp[2] = "K_Hill";
    tmp[3] = "n_Hill";
    tmp[4] = "stretch_flag";
    tmp[5] = "linear_flag";
    setParameterId( tmp );
  }
  
  void StressSpatialSingle::
  derivs(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs ) {
    
    size_t numWalls = T.numWall();
    size_t lengthIndex = variableIndex(0,0);
    size_t dimension = vertexData[0].size();
    
    // Prepare spatial factor
    size_t sI=variableIndex(0,1);
    assert (sI<vertexData[0].size());
    size_t numVertices = vertexData.size();
    double sMax= vertexData[0][sI];
    size_t maxI=0;
    for (size_t i=1; i<numVertices; ++i)
      if (vertexData[i][sI]>sMax) {
	sMax=vertexData[i][sI];
	maxI=i;
      }
    
    for( size_t i=0 ; i<numWalls ; ++i ) {
      size_t v1 = T.wall(i).vertex1()->index();
      size_t v2 = T.wall(i).vertex2()->index();
      double stress=0.0;
      if (!parameter(4)) {
	for (size_t k=0; k<numVariableIndex(1); ++k)
	  stress += wallData[i][variableIndex(1,k)];
      }
      else {
	double distance=0.0;
	for( size_t d=0 ; d<dimension ; ++d )
	  distance += (vertexData[v1][d]-vertexData[v2][d])*
	    (vertexData[v1][d]-vertexData[v2][d]);
	distance = std::sqrt(distance);
	stress = (distance-wallData[i][lengthIndex]) /
	  wallData[i][lengthIndex];
      }
      if (stress > parameter(1)) {
	// Calculate spatial factor
	double maxDistance = sMax - 0.5*(vertexData[v1][sI]+vertexData[v2][sI]);;
	double spatialFactor = Kpow_/(Kpow_+std::pow(maxDistance,parameter(3)));
	
	double growthRate = parameter(0)*(stress - parameter(1))*spatialFactor;
	
	if (parameter(5))
	  growthRate *= wallData[i][lengthIndex];
	wallDerivs[i][lengthIndex] += growthRate;
      }
    }
  }
  
  StressConcentrationHill::
  StressConcentrationHill(std::vector<double> &paraValue, 
			  std::vector< std::vector<size_t> > 
			  &indValue ) {
    
    // Do some checks on the parameters and variable indeces
    //
    if( paraValue.size()!=7 ) {
      std::cerr << "WallGrowth::StressConcentrationHill::"
		<< "StressConcentrationHill() "
		<< "Uses seven parameters k_growthConst, k_growthHill, K_Hill, n_Hill,"
		<< " stretch_threshold stretch_flag and linear_flag" << std::endl;
      exit(0);
    }
    if( paraValue[5] != 0.0 && paraValue[5] != 1.0 ) {
      std::cerr << "WallGrowth::StressConcentrationHill::"
		<< "StressConcentrationHill() "
		<< "stretch_flag parameter must be 0 (stress used) or " 
		<< "1 (stretch used)." << std::endl;
      exit(0);
    }
    if( paraValue[6] != 0.0 && paraValue[6] != 1.0 ) {
      std::cerr << "WallGrowth::StressConcentrationHill::"
		<< "StressConcentrationHill() "
		<< "linear_flag parameter must be 0 (constant growth) or " 
		<< "1 (length dependent growth)." << std::endl;
      exit(0);
    }
    
    if( indValue.size() != 2 || indValue[0].size() != 2 ) {
      std::cerr << "WallGrowth::StressConcentrationHill::"
		<< "StressConcentrationHill() "
		<< "wall length index and concentration index at first "
		<< "level, and spring constant variable indices at second"
		<< std::endl;
      exit(0);
    }
    // Set the variable values
    //
    setId("WallGrowth::StressConcentrationHill");
    setParameter(paraValue);  
    setVariableIndex(indValue);
    
    // Set the parameter identities
    //
    std::vector<std::string> tmp( numParameter() );
    tmp.resize( numParameter() );
    tmp[0] = "k_growth";
    tmp[1] = "stress_threshold";
    tmp[2] = "stretch_flag";
    tmp[3] = "linear_flag";
    setParameterId( tmp );
  }

  void StressConcentrationHill::
  derivs(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs ) {
    
    size_t numWalls = T.numWall();
    size_t lengthIndex = variableIndex(0,0);
    size_t concIndex = variableIndex(0,1);
    double Kpow = std::pow(parameter(2),parameter(3));
    
    for( size_t i=0 ; i<numWalls ; ++i ) {
      size_t v1 = T.wall(i).vertex1()->index();
      size_t v2 = T.wall(i).vertex2()->index();
      double stress=0.0;
      if (!parameter(5)) {
	for (size_t k=0; k<numVariableIndex(1); ++k)
	  stress += wallData[i][variableIndex(1,k)];
      }
      else {
	double distance=0.0;
	for( size_t d=0 ; d<vertexData[v1].size() ; d++ )
	  distance += (vertexData[v1][d]-vertexData[v2][d])*
	    (vertexData[v1][d]-vertexData[v2][d]);
	distance = std::sqrt(distance);
	stress = (distance-wallData[i][lengthIndex]) /
	  wallData[i][lengthIndex];
      }
      if (stress > parameter(4)) {
	// Get the Hill factor from the two cells
	double hillFactor=0.0;
	if (T.wall(i).cell1() != T.background()) {
	  double concpow = std::pow(cellData[T.wall(i).cell1()->index()][concIndex],parameter(3));
	  hillFactor += concpow/(Kpow+concpow); 
	}
	if (T.wall(i).cell2() != T.background()) {
	  double concpow = std::pow(cellData[T.wall(i).cell2()->index()][concIndex],parameter(3));
	  hillFactor += concpow/(Kpow+concpow); 
	}
	double growthRate = (parameter(0)+hillFactor*parameter(1))*(stress - parameter(4));
	if (parameter(6))
	  growthRate *= wallData[i][lengthIndex];
	wallDerivs[i][lengthIndex] += growthRate;
      }
    }
  }
  
  ConstantStressEpidermalAsymmetric::
  ConstantStressEpidermalAsymmetric(std::vector<double> &paraValue, 
				    std::vector< std::vector<size_t> > 
				    &indValue ) {
    
    // Do some checks on the parameters and variable indeces
    //
    if( paraValue.size()!=2 ) {
      std::cerr << "WallGrowth::ConstantStressEpidermalAsymmetric::"
		<< "ConstantStressEpidermalAsymmetric() "
		<< "Uses two parameters k_growth and frac_epi\n";
      exit(0);
    }  
    if( indValue.size() != 1 || indValue[0].size() != 1 ) {
      std::cerr << "WallGrowth::ConstantStressEpidermalAsymmetric::"
		<< "ConstantStressEpidermalAsymmetric() "
		<< "One variable index is used.\n";
      exit(0);
    }
    //Set the variable values
    //
    setId("WallGrowth::ConstantStressEpidermalAsymmetric");
    setParameter(paraValue);  
    setVariableIndex(indValue);
    
    //Set the parameter identities
    //
    std::vector<std::string> tmp( numParameter() );
    tmp.resize( numParameter() );
    tmp[0] = "k_growth";
    tmp[0] = "frac_epi";
    setParameterId( tmp );
  }
  
  void ConstantStressEpidermalAsymmetric::
  derivs(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs ) {
    
    size_t numWalls = T.numWall();
    size_t lengthIndex = variableIndex(0,0);
    
    for( size_t i=0 ; i<numWalls ; ++i ) {
      size_t v1 = T.wall(i).vertex1()->index();
      size_t v2 = T.wall(i).vertex2()->index();
      double kGrowth = parameter(0);
      if( T.wall(i).cell1() == T.background() || 
	  T.wall(i).cell2() == T.background() )
	kGrowth *= parameter(1);
      double distance=0.0;
      for( size_t d=0 ; d<vertexData[v1].size() ; d++ )
	distance += (vertexData[v1][d]-vertexData[v2][d])*
	  (vertexData[v1][d]-vertexData[v2][d]);
      distance = std::sqrt(distance);
      if( distance>wallData[i][lengthIndex] )
	wallDerivs[i][lengthIndex] += kGrowth*
	  (distance-wallData[i][lengthIndex]);
    }
  }
  
  Force::Force(std::vector<double> &paraValue,
	       std::vector< std::vector<size_t> > &indValue)
  {
    if (paraValue.size() != 2) {
      std::cerr << "WallGrowth::Force::Force() "
		<< "Uses two parameters: k_growth and Force_threshold" << std::endl;
      exit(EXIT_FAILURE);
    }
    
    if (indValue.size() != 2 || indValue[0].size() != 1) {
      std::cerr << "WallGrowth::Force::Force() "
		<< "Wall length index must be given in first level.\n"
		<< "Wall force index/indices must be given in second level.\n";
      exit(EXIT_FAILURE);
    }
    
    setId("VertexFromWallSpringExperimental");
    setParameter(paraValue);
    setVariableIndex(indValue);
    
    std::vector<std::string> tmp(numParameter());
    tmp[0] = "k_L";
    tmp[1] = "phi";
    
    setParameterId(tmp);
  }

  void Force::
  derivs(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs)
  {
    for (size_t i = 0; i < T.numWall(); ++i) {
      size_t vertex1Index = T.wall(i).vertex1()->index();
      size_t vertex2Index = T.wall(i).vertex2()->index();
      size_t dimensions = vertexData[vertex1Index].size();
      
      double distance = 0.0;
      for (size_t d = 0; d < dimensions; ++d) {
	distance += (vertexData[vertex1Index][d] - vertexData[vertex2Index][d])
	  * (vertexData[vertex1Index][d] - vertexData[vertex2Index][d]);
      }
      distance = std::sqrt(distance);
      
      double F = 0.0;
      for (size_t j = 0; j < numVariableIndex(1); ++j)
	F += wallData[T.wall(i).index()][variableIndex(1, j)];
      
      double arg = F - parameter(1);
      if (arg > 0)
	wallDerivs[T.wall(i).index()][variableIndex(0, 0)] 
	  += parameter(0) * arg; 
      //*wallData[T.wall(i).index()][variableIndex(0, 0)];
    }
  }
}

MoveVertexRadially::
MoveVertexRadially(std::vector<double> &paraValue, 
		   std::vector< std::vector<size_t> > 
		   &indValue ) {
  
  // Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=2 || ( paraValue[1]!=0 && paraValue[1]!=1) ) {
    std::cerr << "MoveVertexRadially::"
	      << "MoveVertexRadially() "
	      << "Uses two parameters k_growth and r_pow (0,1)\n";
    exit(0);
  }  
  if( indValue.size() != 0 ) {
    std::cerr << "MoveVertexRadially::"
	      << "MoveVertexRadially() "
	      << "No variable index is used.\n";
    exit(0);
  }
  // Set the variable values
  //
  setId("MoveVertexRadially");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  // Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "k_growth";
  tmp[0] = "r_pow";
  setParameterId( tmp );
}

void MoveVertexRadially::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) {
  
  size_t numVertices = T.numVertex();
  size_t dimension=vertexData[0].size();
  
  for( size_t i=0 ; i<numVertices ; ++i ) {
    double fac=parameter(0);
    if( parameter(1)==0.0 ) {
      double r=0.0;
      for( size_t d=0 ; d<dimension ; ++d )
	r += vertexData[i][d]*vertexData[i][d];
      if( r>0.0 )
	r = std::sqrt(r);
      if( r>0.0 )
	fac /= r;
      else
	fac=0.0;
    }
    for( size_t d=0 ; d<dimension ; ++d )
      vertexDerivs[i][d] += fac*vertexData[i][d];
  }
}

MoveVerteX::
MoveVerteX(std::vector<double> &paraValue, 
		   std::vector< std::vector<size_t> > 
		   &indValue ) {
  
  // Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=2 || ( paraValue[1]!=0 && paraValue[1]!=1) ) {
    std::cerr << "MoveVertexX::"
	      << "MoveVertexX() "
	      << "Uses two parameters k_growth and growth_mode (0,1)\n";
    exit(0);
  }  
  if( indValue.size() != 0 ) {
    std::cerr << "MoveVerteX::"
	      << "MoveVerteX() "
	      << "No variable index is used.\n";
    exit(0);
  }
  // Set the variable values
  //
  setId("MoveVertexRadially");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  // Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "k_growth";
  tmp[1] = "growth_mode";
  setParameterId( tmp );
}

void MoveVerteX::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) {
  
  size_t numVertices = T.numVertex();
  size_t s_i = 1; // spatial index
  size_t dimension=vertexData[s_i].size();
  double fac=parameter(0);
  size_t growth_mode = parameter(1);

  
  for( size_t i=0 ; i<numVertices ; ++i ) {
    double x= std::sqrt(vertexData[i][s_i]*vertexData[i][s_i]);
    if( growth_mode == 1 ) {
      fac *= vertexData[i][s_i];
    }

    vertexDerivs[i][s_i] += fac;
  }
}


MoveVertexRadiallycenterTriangulation::
MoveVertexRadiallycenterTriangulation(std::vector<double> &paraValue, 
				      std::vector< std::vector<size_t> > 
				      &indValue ) {
  
  // Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=2 || ( paraValue[1]!=0 && paraValue[1]!=1) ) {
    std::cerr << "MoveVertexRadiallycenterTriangulation::"
	      << "MoveVertexRadiallycenterTriangulation() "
	      << "Uses two parameters k_growth and r_pow (0,1)\n";
    exit(0);
  }  
  if( indValue.size() != 1 || indValue[0].size() != 1 ) {
    std::cerr << "MoveVertexRadiallycenterTriangulation::"
	      << "MoveVertexRadiallycenterTriangulation() " << std::endl
	      << "Start of additional Cell variable indices (center(x,y,z) "
	      << "L_1,...,L_n, n=num vertex) is given in first level." 
	      << std::endl;
    exit(0);
  }
  // Set the variable values
  //
  setId("MoveVertexRadiallycenterTriangulation");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  // Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "k_growth";
  tmp[0] = "r_pow";
  setParameterId( tmp );
}

void MoveVertexRadiallycenterTriangulation::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) {
  
  size_t numVertices = T.numVertex();
  size_t numCells = T.numCell();
  size_t dimension=vertexData[0].size();
  
  // Move vertices
  for( size_t i=0 ; i<numVertices ; ++i ) {
    double fac=parameter(0);
    if( parameter(1)==0.0 ) {
      double r=0.0;
      for( size_t d=0 ; d<dimension ; ++d )
	r += vertexData[i][d]*vertexData[i][d];
      if( r>0.0 )
	r = std::sqrt(r);
      if( r>0.0 )
	fac /= r;
      else
	fac=0.0;
    }
    for( size_t d=0 ; d<dimension ; ++d )
      vertexDerivs[i][d] += fac*vertexData[i][d];
  }
  // Move vertices defined in cell centers
  for( size_t i=0 ; i<numCells ; ++i ) {
    double fac=parameter(0);
    if( parameter(1)==0.0 ) {
      double r=0.0;
      for( size_t d=variableIndex(0,0) ; d<variableIndex(0,0)+dimension ; ++d )
	r += cellData[i][d]*cellData[i][d];
      if( r>0.0 )
	r = std::sqrt(r);
      if( r>0.0 )
	fac /= r;
      else
	fac=0.0;
    }
    for( size_t d=variableIndex(0,0) ; d<variableIndex(0,0)+dimension ; ++d )
      cellDerivs[i][d] += fac*cellData[i][d];
  }
}

MoveVertexSphereCylinder::
MoveVertexSphereCylinder(std::vector<double> &paraValue, 
			 std::vector< std::vector<size_t> > 
			 &indValue ) {
  
  // Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=2 || ( paraValue[1]!=0 && paraValue[1]!=1) ) {
    std::cerr << "MoveVertexSphereCylinder::"
	      << "MoveVertexSphereCylinder() "
	      << "Uses two parameters k_growth and r_pow (0,1)\n";
    exit(0);
  }  
  if( indValue.size() != 0 ) {
    std::cerr << "MoveVertexSphereCylinder::"
	      << "MoveVertexSphereCylinder() "
	      << "No variable index is used.\n";
    exit(0);
  }
  // Set the variable values
  //
  setId("MoveVertexSphereCylinder");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  // Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "k_growth";
  tmp[0] = "r_pow";
  setParameterId( tmp );
}

void MoveVertexSphereCylinder::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) {
  
  size_t numVertices = T.numVertex();
  if (vertexData[0].size()!=3) {
    std::cerr << "MoveVertexSphereCylinder:: Only works for 3 dimensions." << std::endl;
    exit(-1);
  }
  size_t xI=0;
  size_t yI=1;
  size_t zI=2;
 
  for( size_t i=0 ; i<numVertices ; ++i ) {
    if (vertexData[i][zI]<0.0) { // on cylinder
      if( parameter(1)==0.0 ) {
	vertexDerivs[i][zI] -= parameter(0);
      }
      else {
	double r = std::sqrt(vertexData[i][xI]*vertexData[i][xI]+
			     vertexData[i][yI]*vertexData[i][yI]);
	vertexDerivs[i][zI] -= parameter(0)*(3.14159265*0.5*r-vertexData[i][zI]);
      }
    }
    else { // on half sphere
      double r = std::sqrt(vertexData[i][xI]*vertexData[i][xI]+
			   vertexData[i][yI]*vertexData[i][yI]+
			   vertexData[i][zI]*vertexData[i][zI]);
      double rPrime = std::sqrt(vertexData[i][xI]*vertexData[i][xI]+
				vertexData[i][yI]*vertexData[i][yI]);
      double theta = std::asin(rPrime/r);

      double fac=parameter(0)*theta;
      if (parameter(0)==1) {
	fac *= r;
      }
      vertexDerivs[i][xI] += fac*vertexData[i][xI]*vertexData[i][zI]/rPrime;
      vertexDerivs[i][yI] += fac*vertexData[i][yI]*vertexData[i][zI]/rPrime;
      vertexDerivs[i][zI] -= fac*rPrime;
    }
  }
}

WaterVolumeFromTurgor::
WaterVolumeFromTurgor(std::vector<double> &paraValue,
		      std::vector< std::vector<size_t> > &indValue)
{
  if (paraValue.size() != 5) {
    std::cerr << "WaterVolumeFromTurgor::WaterVolumeFromTurgor() "
	      << "Uses five parameters: k_p, P_max, k_pp and "
	      << "denyShrink_flag allowNegTurgor_flag." << std::endl;
    exit(EXIT_FAILURE);
  }
  
  if (indValue.size() < 1 || indValue.size() > 2
      || indValue[0].size() != 1 
      || (indValue.size()==2 && indValue[1].size() != 1 ) ) {
    std::cerr << "WaterVolumeFromTurgor::WaterVolumeFromTurgor() "
	      << "Water volume index must be given in "
	      << "first level.\n"
	      << "Optionally index for saving the turgor pressure can be"
	      << " given at second level." << std::endl; 		
    exit(EXIT_FAILURE);
  }
  
  setId("WaterVolumeFromTurgor");
  setParameter(paraValue);
  setVariableIndex(indValue);
  
  std::vector<std::string> tmp(numParameter());
  tmp[0] = "k_p";
  tmp[1] = "P_max";
  tmp[2] = "k_pp";
  tmp[3] = "denyShrink_flag";
  tmp[4] = "allowNegTurgor_flag";
  setParameterId(tmp);
}

void WaterVolumeFromTurgor::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs)
{
  for (size_t n = 0; n < T.numCell(); ++n) {
    Cell cell = T.cell(n);
    
    double P = 0.0;
    double totalLength = 0.0;
    for (size_t i = 0; i < cell.numWall(); ++i) {
      size_t vertex1Index = cell.wall(i)->vertex1()->index();
      size_t vertex2Index = cell.wall(i)->vertex2()->index();
      size_t dimensions = vertexData[vertex1Index].size();
      
      double distance = 0.0;
      for (size_t d = 0; d < dimensions; ++d) {
	distance += (vertexData[vertex1Index][d] - vertexData[vertex2Index][d])
	  * (vertexData[vertex1Index][d] - vertexData[vertex2Index][d]);
      }
      distance = std::sqrt(distance);
      totalLength += distance;
      
      // Old turgor measure from wall extensions
      //for (size_t j = 0; j < numVariableIndex(1); ++j)
      //P += wallData[cell.wall(i)->index()][variableIndex(1, j)]/distance;
    }
    
    // Calculate turgor measure from volume and 'water volume'
    // P ~ p_2(V_w-V)/V
    double cellVolume =cell.calculateVolume(vertexData);
    P = (cellData[cell.index()][variableIndex(0,0)]-cellVolume) / cellVolume;
    if (P<0.0 && !parameter(4))
      P=0.0;
    
    P *= parameter(2);
    
    if (numVariableIndexLevel()==2)
      cellData[n][variableIndex(1,0)]=P;
    
    if( !parameter(3) || parameter(1)-P>0.0 )
      cellDerivs[cell.index()][variableIndex(0,0)] += 
	parameter(0) * (parameter(1) - P) * totalLength;
  }
}

DilutionFromVertexDerivs::
DilutionFromVertexDerivs(std::vector<double> &paraValue,
			 std::vector< std::vector<size_t> > &indValue)
{
  if (paraValue.size()) {
    std::cerr << "DilutionFromVertexDerivs::DilutionFromVertexDerivs() "
	      << "Uses no parameters." << std::endl;
    exit(EXIT_FAILURE);
  }
  
  if (indValue.size() != 1 || indValue[0].size() < 1) {
    std::cerr << "DilutionFromVertexDerivs::DilutionFromVertexDerivs() "
	      << "List of concentration variable index must be given in "
	      << "first level." << std::endl;
    exit(EXIT_FAILURE);
  }
  
  setId("DilutionFromVertexDerivs");
  setParameter(paraValue);
  setVariableIndex(indValue);
  
  std::vector<std::string> tmp(numParameter());
  setParameterId(tmp);
}

void DilutionFromVertexDerivs::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs)
{
  size_t dimension;
  dimension = vertexData[0].size();
  assert(dimension==2);
  
  for (size_t n = 0; n < T.numCell(); ++n) {
    Cell cell = T.cell(n);
    double area = cell.calculateVolume(vertexData,1);
    
    double areaDerivs=0.0;
    for( size_t k=0 ; k<cell.numVertex() ; ++k ) {
      size_t vI = cell.vertex(k)->index();
      size_t vIPlus = cell.vertex((k+1)%(cell.numVertex()))->index();
      areaDerivs += vertexData[vIPlus][1]*vertexDerivs[vI][0] - 
	vertexData[vI][1]*vertexDerivs[vIPlus][0] -
	vertexData[vIPlus][0]*vertexDerivs[vI][1] +
	vertexData[vI][0]*vertexDerivs[vIPlus][1];
    }
    
    double fac = areaDerivs/area;
    for (size_t k=0; k<numVariableIndex(0); ++k)
      cellDerivs[n][variableIndex(0,k)] -= cellData[n][variableIndex(0,k)]*
	fac;
  }
}
