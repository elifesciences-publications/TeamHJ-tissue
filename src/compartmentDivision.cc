//
// Filename     : compartmentDivision.cc
// Description  : Classes describing compartmentDivision updates
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : July 2006
// Revision     : $Id:$
//
#include<limits>

#include "baseCompartmentChange.h"
#include "compartmentDivision.h"
#include "myRandom.h"
#include "myMath.h"

namespace Division {

  VolumeViaLongestWall::
  VolumeViaLongestWall(std::vector<double> &paraValue, 
		       std::vector< std::vector<size_t> > 
		       &indValue ) {
    //
    // Do some checks on the parameters and variable indeces
    //
    if( paraValue.size()!=3 ) {
      std::cerr << "DivisionVolumeViaLongestWall::"
		<< "DivisionVolumeViaLongestWall() "
		<< "Three parameters used V_threshold, LWall_frac, and "
		<< "Lwall_threshold." << std::endl;
      exit(EXIT_FAILURE);
  }
    if( indValue.size() != 1 ) {
      std::cerr << "DivisionVolumeViaLongestWall::"
		<< "DivisionVolumeViaLongestWall() "
		<< "Variable indices for volume dependent cell "
		<< "variables is used.\n";
      exit(EXIT_FAILURE);
    }
    //
    // Set the variable values
    //
    setId("Division::VolumeViaLongestWall");
    setNumChange(1);
    setParameter(paraValue);  
    setVariableIndex(indValue);
    //
    // Set the parameter identities
    //
    std::vector<std::string> tmp( numParameter() );
    tmp.resize( numParameter() );
    tmp[0] = "V_threshold";
    tmp[1] = "LWall_frac";
    tmp[2] = "LWall_threshold";
    setParameterId( tmp );
  }

  int VolumeViaLongestWall::
  flag(Tissue *T,size_t i,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) {
    
    if( T->cell(i).calculateVolume(vertexData) > parameter(0) ) {
      std::cerr << "Cell " << i << " marked for division with volume " 
		<< T->cell(i).volume() << std::endl;
      return 1;
    } 
    return 0;
  }
  
  void VolumeViaLongestWall::
  update(Tissue *T,size_t i,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
       DataMatrix &vertexData,
	 DataMatrix &cellDeriv,
	 DataMatrix &wallDeriv,
	 DataMatrix &vertexDeriv ) {
    
    Cell *divCell = &(T->cell(i));
    size_t dimension = vertexData[0].size();
    assert( divCell->numWall() > 1 );
    assert( dimension==2 || dimension==3 ); 
    //
    // Find longest wall
    // 
    size_t wI=0,w3I=divCell->numWall();
    double maxLength = divCell->wall(0)->setLengthFromVertexPosition(vertexData);
    for( size_t k=1 ; k<divCell->numWall() ; ++k ) {
      double tmpLength = divCell->wall(k)->setLengthFromVertexPosition(vertexData);
      if( tmpLength > maxLength ) {
	wI=k;
	maxLength = tmpLength;
      }
    }   
    
    //
    // Find position for first new vertex
    //
    std::vector<double> nW(dimension),nW2(dimension),v1Pos(dimension),
      v2Pos(dimension);
    size_t v1wI = divCell->wall(wI)->vertex1()->index();
    size_t v2wI = divCell->wall(wI)->vertex2()->index();
    for( size_t d=0 ; d<dimension ; ++d ) {
      nW[d] = (vertexData[v1wI][d]-vertexData[v2wI][d])/maxLength;
      v1Pos[d] = 0.5*(vertexData[v1wI][d]+vertexData[v2wI][d]);
    }
    //
    // Find intersection with another wall via vector perpendicular to first wall
    //
    if (dimension==2) {
      nW2[1] = nW[0];
      nW2[0] = -nW[1];
    }
    else if (dimension==3) {
      nW2[0]=nW[0];
      nW2[1]=nW[1];
      nW2[2]=nW[2];
    }
    if (findSecondDivisionWall(vertexData, divCell, wI, w3I, v1Pos, nW2, v2Pos)) {
      std::cerr << "Division::VolumeViaLongestWall::update "
		<< "failed to find the second wall for division!" << std::endl;
      exit(EXIT_FAILURE);
    }
    //
    // Do the division (add one cell, three walls, and two vertices)
    //
    size_t numWallTmp=wallData.size();
    assert( numWallTmp==T->numWall() );
    //Divide
    T->divideCell(divCell,wI,w3I,v1Pos,v2Pos,cellData,wallData,vertexData,
		  cellDeriv,wallDeriv,vertexDeriv,variableIndex(0),
		  parameter(2));
    assert( numWallTmp+3 == T->numWall() );
    
    //Change length of new wall between the divided daugther cells 
    wallData[numWallTmp][0] *= parameter(1);
    
    //Check that the division did not mess up the data structure
    //T->checkConnectivity(1);	
  }
  
  VolumeViaLongestWallCenterTriangulation::
  VolumeViaLongestWallCenterTriangulation(std::vector<double> &paraValue, 
					  std::vector< std::vector<size_t> > 
					  &indValue ) {
    //
    // Do some checks on the parameters and variable indeces
    //
    if( paraValue.size()!=3 && paraValue.size()!=4) {
      std::cerr << "DivisionVolumeViaLongestWallCenterTriangulation::"
		<< "DivisionVolumeViaLongestWallCenterTriangulation() "
		<< "Three parameters used V_threshold, LWall_frac, and "
		<< "Lwall_threshold and optionally forth parameter(=1) "
		<< "if double resting length option used" 
		<< std::endl;
      exit(0);
    }
    if( indValue.size() != 2 || indValue[1].size() !=2 ) {
      std::cerr << "DivisionVolumeViaLongestWallCenterTriangulation::"
		<< "DivisionVolumeViaLongestWallCenterTriangulation() "
		<< "Variable indices for volume dependent cell "
		<< "at first level, Start of additional Cell variable indices (center(x,y,z) "
		<< "and Wall length index and at second level" 
		<< std::endl;
      exit(0);
    }
    //
    // Set the variable values
    //
    setId("Division::VolumeViaLongestWallCenterTriangulation");
    setNumChange(1);
    setParameter(paraValue);  
    setVariableIndex(indValue);
    //
    // Set the parameter identities
    //
    std::vector<std::string> tmp( numParameter() );
    tmp.resize( numParameter() );
    tmp[0] = "V_threshold";
    tmp[1] = "LWall_frac";
    tmp[2] = "LWall_threshold";
    if (numParameter()==4)
      tmp[3]="doubleFlag";
    setParameterId( tmp );
  }
  
  int VolumeViaLongestWallCenterTriangulation::
  flag(Tissue *T,size_t i,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) {
    
    if( T->cell(i).calculateVolume(vertexData) > parameter(0) ) {
      std::cerr << "Cell " << i << " marked for division with volume " 
		<< T->cell(i).volume() << std::endl;
      return 1;
    } 
    return 0;
  }
  
  void VolumeViaLongestWallCenterTriangulation::
  update(Tissue *T,size_t i,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDeriv,
	 DataMatrix &wallDeriv,
	 DataMatrix &vertexDeriv ) {
    
    Cell *divCell = &(T->cell(i));
    size_t dimension = vertexData[0].size();
    assert( divCell->numWall() > 1 );
    assert( dimension==2 || dimension==3 ); 
    //
    // Find longest wall
    // 
    size_t wI=0,w3I=divCell->numWall();
    double maxLength = divCell->wall(0)->setLengthFromVertexPosition(vertexData);
    for( size_t k=1 ; k<divCell->numWall() ; ++k ) {
      double tmpLength = divCell->wall(k)->setLengthFromVertexPosition(vertexData);
      if( tmpLength > maxLength ) {
	wI=k;
	maxLength = tmpLength;
      }
    }   
    
    //
    // Find position for first new vertex
    //
    std::vector<double> nW(dimension),nW2(dimension),v1Pos(dimension),
      v2Pos(dimension);
    size_t v1wI = divCell->wall(wI)->vertex1()->index();
    size_t v2wI = divCell->wall(wI)->vertex2()->index();
    for( size_t d=0 ; d<dimension ; ++d ) {
      nW[d] = (vertexData[v1wI][d]-vertexData[v2wI][d])/maxLength;
      v1Pos[d] = 0.5*(vertexData[v1wI][d]+vertexData[v2wI][d]);
    }
    //
    // Find intersection with another wall via vector perpendicular to first wall
    //
    if (dimension==2) {
      nW2[1] = nW[0];
      nW2[0] = -nW[1];
    }
    else if (dimension==3) {
      nW2[0]=nW[0];
      nW2[1]=nW[1];
      nW2[2]=nW[2];
    }
    if (findSecondDivisionWall(vertexData, divCell, wI, w3I, v1Pos, nW2, v2Pos)) {
      std::cerr << "DivisionVolumeViaLongestWallCenterTriangulation::update "
		<< "failed to find the second wall for division!" << std::endl;
      exit(-1);
    }
    //
    // Do the division (add one cell, three walls, and two vertices)
    //
    size_t numWallTmp=wallData.size();
    assert( numWallTmp==T->numWall() );
    //Divide
    if (numParameter()==3)
      T->divideCellCenterTriangulation(divCell,wI,w3I,
				       variableIndex(1,0),variableIndex(1,1),
				       v1Pos,v2Pos,cellData,wallData,vertexData,
				       cellDeriv,wallDeriv,vertexDeriv,variableIndex(0),
				       parameter(2),
				       0);
    if (numParameter()==4 && parameter(3)==1)
      T->divideCellCenterTriangulation(divCell,wI,w3I,
				       variableIndex(1,0),variableIndex(1,1),
				       v1Pos,v2Pos,cellData,wallData,vertexData,
				       cellDeriv,wallDeriv,vertexDeriv,variableIndex(0),
				       parameter(2),
				       1);
    assert( numWallTmp+3 == T->numWall() );
  
    //Change length of new wall between the divided daugther cells 
    wallData[numWallTmp][0] *= parameter(1);
    
    //Check that the division did not mess up the data structure
    //T->checkConnectivity(1);	
  }
  
  VolumeViaLongestWall3DCenterTriangulation::
  VolumeViaLongestWall3DCenterTriangulation(std::vector<double> &paraValue, 
					    std::vector< std::vector<size_t> > 
					    &indValue ) {
    
    // Do some checks on the parameters and variable indeces
    //
    if( paraValue.size()!=3 ) {
      std::cerr << "DivisionVolumeViaLongestWall3DCenterTriangulation::"
		<< "DivisionVolumeViaLongestWall3DCenterTriangulation() "
		<< "Three parameters used V_threshold, LWall_frac LWall_threshold\n";
      exit(0);
    }
    if( indValue.size() != 2 || indValue[1].size() !=2 ) {
      std::cerr << "DivisionVolumeViaLongestWall3DCenterTriangulation::"
		<< "DivisionVolumeViaLongestWall3DCenterTriangulation() "
		<< "Variable indices for volume dependent cell "
		<< "at first level, Start of additional Cell variable indices (center(x,y,z) "
		<< "and Wall length index and at second level" 
		<< std::endl;
      exit(0);
    }
    //Set the variable values
    //
    setId("DivisionVolumeViaLongestWall3DCenterTriangulation");
    setNumChange(1);
    setParameter(paraValue);  
    setVariableIndex(indValue);
    //
    //Set the parameter identities
    //
    std::vector<std::string> tmp( numParameter() );
    tmp.resize( numParameter() );
    tmp[0] = "V_threshold";
    tmp[1] = "LWall_frac";
    tmp[2] = "LWall_threshold";
    setParameterId( tmp );
  }
  
  int VolumeViaLongestWall3DCenterTriangulation::
  flag(Tissue *T,size_t i,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) {
    
    if( T->cell(i).calculateVolume(vertexData) > parameter(0) ) {
      std::cerr << "Cell " << i << " marked for division with volume " 
		<< T->cell(i).volume() << std::endl;
      return 1;
    } 
    return 0;
  }
  
  void VolumeViaLongestWall3DCenterTriangulation::
  update(Tissue *T,size_t i,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDeriv,
	 DataMatrix &wallDeriv,
	 DataMatrix &vertexDeriv ) {
    
    Cell *divCell = &(T->cell(i));
    size_t dimension = vertexData[0].size();
    assert( divCell->numWall() > 1 );
    assert( dimension==3 );
    
    //Find longest wall
    //
    size_t wI=0, w3I=divCell->numWall();
    double maxLength = divCell->wall(0)->setLengthFromVertexPosition(vertexData);
    for( size_t k=1 ; k<divCell->numWall() ; ++k ) {
      double tmpLength = divCell->wall(k)->setLengthFromVertexPosition(vertexData);
      if( tmpLength > maxLength ) {
	wI=k;
	maxLength = tmpLength;
      }
    }   
    
    std::vector<double> nW(dimension),nW2(dimension),v1Pos(dimension),
      v2Pos(dimension);
    size_t v1wI = divCell->wall(wI)->vertex1()->index();
    size_t v2wI = divCell->wall(wI)->vertex2()->index();
    
    for( size_t d=0 ; d<dimension ; ++d ) {
      nW[d] = (vertexData[v1wI][d]-vertexData[v2wI][d])/maxLength;
      v1Pos[d] = 0.5*(vertexData[v1wI][d]+vertexData[v2wI][d]);
    }
    
    //Find intersection with another wall by looking at perpendicular plane
    //
    
    //double minDist,w3s;
    std::vector<size_t> w3Tmp;
    std::vector<double> w3tTmp;
    int flag=0;
    for( size_t k=0 ; k<divCell->numWall() ; ++k ) {
      if( k!=wI ) {
	size_t v1w3Itmp = divCell->wall(k)->vertex1()->index();
	size_t v2w3Itmp = divCell->wall(k)->vertex2()->index();
	std::vector<double> w3(dimension),w0(dimension);
	double fac1=0.0,fac2=0.0;
	for( size_t d=0 ; d<dimension ; ++d ) {
	  w3[d] = vertexData[v2w3Itmp][d]-vertexData[v1w3Itmp][d];
	  fac1 += nW[d]*(v1Pos[d]-vertexData[v1w3Itmp][d]);
	  fac2 += nW[d]*w3[d]; 
	}
	if( fac2 != 0.0 ) {//else parallell and not applicable
	  double t = fac1/fac2;
	  if( t>=0.0 && t<1.0 ) {//within wall
	    std::cerr << "wall " << k << " (t=" << t << " " << fac1 << "/"
		      << fac2 
		      << ") chosen as second wall"
		      << std::endl;					
	    flag++;
	    w3I = k;
	    w3Tmp.push_back(k);
	    w3tTmp.push_back(t);
	  }
	}
      }
    }
    assert( w3I != wI );
    assert( w3I != divCell->numWall() );
    if( flag != 1 ) {
      std::cerr << "divideVolumeViaLongestWall3DCenterTriangulation::update() Warning"
		<< " more than one wall possible as connection "
		<< "for cell " 
		<< i << std::endl; 
      for( size_t k=0 ; k<divCell->numWall() ; ++k ) {
	std::cerr << "0 " 
		  << vertexData[divCell->wall(k)->vertex1()->index()][0]
		  << " " 
		  << vertexData[divCell->wall(k)->vertex1()->index()][1]
		  << " " 
		  << vertexData[divCell->wall(k)->vertex1()->index()][2]
		  << "\n0 " 
		  << vertexData[divCell->wall(k)->vertex2()->index()][0]
		  << " " 
		  << vertexData[divCell->wall(k)->vertex2()->index()][1]
		  << " " 
		  << vertexData[divCell->wall(k)->vertex2()->index()][2]
		  << "\n\n\n";
      }
      for( size_t kk=0 ; kk<w3Tmp.size() ; ++kk ) {
	size_t k = w3Tmp[kk];
	std::cerr << "1 " 
		  << vertexData[divCell->wall(k)->vertex1()->index()][0]
		  << " " 
		  << vertexData[divCell->wall(k)->vertex1()->index()][1]
		  << " " 
		  << vertexData[divCell->wall(k)->vertex1()->index()][2]
		  << "\n1 " 
		  << vertexData[divCell->wall(k)->vertex2()->index()][0]
		  << " " 
		  << vertexData[divCell->wall(k)->vertex2()->index()][1]
		  << " " 
		  << vertexData[divCell->wall(k)->vertex2()->index()][2]
		  << "\n\n\n";
      }
      std::cerr << "2 " 
		<< vertexData[divCell->wall(wI)->vertex1()->index()][0]
		<< " " 
		<< vertexData[divCell->wall(wI)->vertex1()->index()][1]
		<< " " 
		<< vertexData[divCell->wall(wI)->vertex1()->index()][2]
		<< "\n2 " 
		<< vertexData[divCell->wall(wI)->vertex2()->index()][0]
		<< " " 
		<< vertexData[divCell->wall(wI)->vertex2()->index()][1]
		<< " " 
		<< vertexData[divCell->wall(wI)->vertex2()->index()][2]
		<< "\n\n\n";
      std::cerr << "3 " 
		<< vertexData[divCell->wall(w3I)->vertex1()->index()][0]
		<< " " 
		<< vertexData[divCell->wall(w3I)->vertex1()->index()][1]
		<< " " 
		<< vertexData[divCell->wall(w3I)->vertex1()->index()][2]
		<< "\n3 " 
		<< vertexData[divCell->wall(w3I)->vertex2()->index()][0]
		<< " " 
		<< vertexData[divCell->wall(w3I)->vertex2()->index()][1]
		<< " " 
		<< vertexData[divCell->wall(w3I)->vertex2()->index()][2]
		<< "\n\n\n";
      std::cerr << "4 " 
		<< 0.5*(vertexData[divCell->wall(wI)->vertex1()->index()][0]+
			vertexData[divCell->wall(wI)->vertex2()->index()][0])
		<< " " 
		<< 0.5*(vertexData[divCell->wall(wI)->vertex1()->index()][1]+
			vertexData[divCell->wall(wI)->vertex2()->index()][1])
		<< " " 
		<< 0.5*(vertexData[divCell->wall(wI)->vertex1()->index()][2]+
			vertexData[divCell->wall(wI)->vertex2()->index()][2])
		<< "\n4 "
		<< 0.5*(vertexData[divCell->wall(w3I)->vertex1()->index()][0]+
			vertexData[divCell->wall(w3I)->vertex2()->index()][0])
		<< " " 
		<< 0.5*(vertexData[divCell->wall(w3I)->vertex1()->index()][1]+
			vertexData[divCell->wall(w3I)->vertex2()->index()][1])
		<< " " 
		<< 0.5*(vertexData[divCell->wall(w3I)->vertex1()->index()][2]+
			vertexData[divCell->wall(w3I)->vertex2()->index()][2])
		<< "\n\n\n";
      exit(-1);
    }
    
    size_t v1w3I = divCell->wall(w3I)->vertex1()->index();
    size_t v2w3I = divCell->wall(w3I)->vertex2()->index();
    for( size_t d=0 ; d<dimension ; ++d )
      v2Pos[d] = vertexData[v1w3I][d] + 
	w3tTmp[w3tTmp.size()-1]*(vertexData[v2w3I][d]-vertexData[v1w3I][d]);
    //v2Pos[d] = 0.5*(vertexData[v1w3I][d]+vertexData[v2w3I][d]);
    
    //Add one cell, three walls, and two vertices
    //
    //Save number of walls
    size_t numWallTmp=wallData.size();
    assert( numWallTmp==T->numWall() );
    //Divide
    T->divideCellCenterTriangulation(divCell,wI,w3I,
				     variableIndex(1,0),variableIndex(1,1),
				     v1Pos,v2Pos,cellData,wallData,vertexData,
				     cellDeriv,wallDeriv,vertexDeriv,variableIndex(0),
				     parameter(2),
                                   0);
    
    assert( numWallTmp+3 == T->numWall() );
    
    //Change length of new wall between the divided daugther cells 
    wallData[numWallTmp][0] *= parameter(1);
    
    //Check that the division did not messed up the data structure
    //T->checkConnectivity(1);	
  }
  
  Branching::
  Branching(std::vector<double> &paraValue, 
	    std::vector< std::vector<size_t> > 
	    &indValue ) {
    //
    // Do some checks on the parameters and variable indeces
    //
    if( paraValue.size()!=4 ) {
      std::cerr << "Branching::"
		<< "Branching() "
		<< "four parameters used firstConc_threshold, secondConc_threshold, LWall_frac, and "
		<< "Lwall_threshold " << std::endl;
      // if( paraValue.size()!=5 ) {
      //   std::cerr << "Branching::"
      //             << "Branching() "
      //             << "five parameters used firstConc_threshold, secondConc_threshold, LWall_frac, and "
      //             << "Lwall_threshold and branching_wall_bending_stiffness." << std::endl;
      
      exit(0);
    }
    if( indValue.size() != 1 ) {
      std::cerr << "Branching::"
		<< "Branching() "
		<< "five variable indices for cell branching flag and first and second concentrations  "
		<< "and mt-vector initial index and wallbranching flag are used.\n";
      exit(0);
    }
    //
    // Set the variable values
    //
    setId("Branching");
    setNumChange(1);
    setParameter(paraValue);  
    setVariableIndex(indValue);
    //
    // Set the parameter identities
    //
    std::vector<std::string> tmp( numParameter() );
    tmp.resize( numParameter() );
    tmp[0] = "firstConc_threshold";
    tmp[1] = "secondConc_threshold";
    tmp[2] = "LWall_frac";
    tmp[3] = "LWall_threshold";
    //tmp[4] = "branching_wall_bending_stiffness";
    setParameterId( tmp );
  }
  
  int Branching::
  flag(Tissue *T,size_t i,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) {
    
    //  if( T->cell(i).calculateVolume(vertexData) > parameter(0) && cellData[i][11]==0 ) {
    if( cellData[i][variableIndex(0,0)]==0 && cellData[i][variableIndex(0,1)] > parameter(0)
	&& cellData[i][variableIndex(0,2)] < parameter(1) ) {
      std::cerr << "Cell " << i << " marked for branching " 
		<< T->cell(i).volume() << std::endl;
      cellData[i][variableIndex(0,0)]=1;
      return 1;
    }
    
    return 0;
  }
  
  void Branching::
  update(Tissue *T,size_t i,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDeriv,
	 DataMatrix &wallDeriv,
	 DataMatrix &vertexDeriv ) {
    
    //size_t wallFlag=variableIndex(0,4);
    
    Cell *brCell = &(T->cell(i));
    size_t dimension = vertexData[0].size();
    if( dimension > 2 ) {
      std::cerr << "Branching::"
		<< "Branching() "
		<< "only for two-dimensional tissue "
		<< std::endl;
      exit(0);
    }
    // assert( divCell->numWall() > 1 );
    // assert( dimension==2 || dimension==3 ); 
    
    //i=1;
    // the cell has a branch        
    //
    // Find longest wall
    // 
    size_t wI=0,w3I=brCell->numWall();
    double maxLength = brCell->wall(0)->setLengthFromVertexPosition(vertexData);
    for( size_t k=1 ; k<brCell->numWall() ; ++k ) {
      double tmpLength = brCell->wall(k)->setLengthFromVertexPosition(vertexData);
      if( tmpLength > maxLength ) {
	wI=k;
	maxLength = tmpLength;
      }
    }   
  
  
    std::vector<size_t> wallsBack; 
    // finding the walls that are facing the background
    for( size_t k=0 ; k<brCell->numWall() ; ++k ) 
      if( (brCell ->wall(k)-> cell1() == T-> background()) || (brCell ->wall(k)-> cell2() == T -> background()) )
	wallsBack.push_back(k);
  
    // randomly choose one of them if more than one wall is facing the background 
    size_t s=wallsBack.size();
    if (wallsBack.size()==1)
      wI=wallsBack[0];
    else{ 
      size_t s=wallsBack.size();
      wI=wallsBack[rand() % s];
    }
  
    maxLength = brCell->wall(wI)->setLengthFromVertexPosition(vertexData);
  
    double minLength=brCell->wall((wI+1)%(brCell->numWall()))->setLengthFromVertexPosition(vertexData);  
  
  
    // Find position for first two new vertices on the wall
    std::vector<double> nW(dimension),nW2(dimension),v1Pos(dimension),
      v2Pos(dimension),v3Pos(dimension),v4Pos(dimension);
    size_t v1wI = brCell->wall(wI)->vertex1()->index();
    size_t v2wI = brCell->wall(wI)->vertex2()->index();
    for( size_t d=0 ; d<dimension ; ++d ) {
      nW[d] = (vertexData[v1wI][d]-vertexData[v2wI][d])/maxLength;
      v1Pos[d] = 0.3*vertexData[v1wI][d]+0.7*vertexData[v2wI][d];
      v2Pos[d] = (0.3+minLength/maxLength)*vertexData[v1wI][d]+(0.7-minLength/maxLength)*vertexData[v2wI][d];
    
    }

    // Finding the positions of other two vertices
    //next wall
    size_t wIplusOne=(wI+1)%(brCell ->numWall());
    size_t v1wIplusOne = brCell->wall(wIplusOne)->vertex1()->index();
    size_t v2wIplusOne = brCell->wall(wIplusOne)->vertex2()->index();

    //outward normal to the wall wI
    std::vector<double>  wallvecNorm(dimension);
    std::vector<double>  wallvec(dimension);
    if ( v1wIplusOne==v1wI ){
      wallvec[0]=vertexData[v2wI][0]-vertexData[v1wI][0];
      wallvec[1]=vertexData[v2wI][1]-vertexData[v1wI][1];
    
      wallvecNorm[0]=-vertexData[v2wIplusOne][0]+vertexData[v1wIplusOne][0];
      wallvecNorm[1]=-vertexData[v2wIplusOne][1]+vertexData[v1wIplusOne][1];
    }


    else if( v1wIplusOne==v2wI ){
      wallvec[0]=vertexData[v1wI][0]-vertexData[v2wI][0];
      wallvec[1]=vertexData[v1wI][1]-vertexData[v2wI][1];
    
      wallvecNorm[0]=-vertexData[v2wIplusOne][0]+vertexData[v1wIplusOne][0];
      wallvecNorm[1]=-vertexData[v2wIplusOne][1]+vertexData[v1wIplusOne][1];
    } 
    else if( v2wIplusOne==v1wI ){
      wallvec[0]=vertexData[v2wI][0]-vertexData[v1wI][0];
      wallvec[1]=vertexData[v2wI][1]-vertexData[v1wI][1];
    
      wallvecNorm[0]=-vertexData[v1wIplusOne][0]+vertexData[v2wIplusOne][0];
      wallvecNorm[1]=-vertexData[v1wIplusOne][1]+vertexData[v2wIplusOne][1];
    } 
    else if( v2wIplusOne==v2wI ){
      wallvec[0]=vertexData[v1wI][0]-vertexData[v2wI][0];
      wallvec[1]=vertexData[v1wI][1]-vertexData[v2wI][1];
    
      wallvecNorm[0]=-vertexData[v1wIplusOne][0]+vertexData[v2wIplusOne][0];
      wallvecNorm[1]=-vertexData[v1wIplusOne][1]+vertexData[v2wIplusOne][1];
    } 





    double temp=std::sqrt(wallvec[0]*wallvec[0]+wallvec[1]*wallvec[1]);
    if (temp!=0){
      wallvec[0]/=temp;
      wallvec[1]/=temp;
    }


    temp=wallvec[0]*wallvecNorm[0]+wallvec[1]*wallvecNorm[1];
  
    wallvecNorm[0]-=temp*wallvec[0];
    wallvecNorm[1]-=temp*wallvec[1];

    temp=std::sqrt(wallvecNorm[0]*wallvecNorm[0]+wallvecNorm[1]*wallvecNorm[1]);
    if (temp!=0){
      wallvecNorm[0]/=temp;
      wallvecNorm[1]/=temp;
    }
  
    v4Pos[0] = v1Pos[0]+.2*wallvecNorm[0];
    v4Pos[1] = v1Pos[1]+.2*wallvecNorm[1];

    v3Pos[0] = v2Pos[0]+.2*wallvecNorm[0];
    v3Pos[1] = v2Pos[1]+.2*wallvecNorm[1];


    //std::cerr<<"wallvec   "<<wallvec[0]<<"  " <<wallvec[1] <<std::endl;
    //std::cerr<<"wallvec norm   "<<wallvecNorm[0]<<"  " <<wallvecNorm[1] <<std::endl;
    //
    // Do the branching (add one cell, five walls, and four vertices)
    //
    size_t numWallTmp=wallData.size();
    size_t numCellTmp=cellData.size();
    assert( numWallTmp==T->numWall() );
 
   
    //branching
    T->branchCell(brCell,wI,v1Pos,v2Pos,v3Pos,v4Pos,cellData,wallData,vertexData,
		  cellDeriv,wallDeriv,vertexDeriv, parameter(2));
    assert( numWallTmp+5 == T->numWall() );
  
    for(size_t kk=0; kk< cellData[0].size(); ++kk)
      cellData[numCellTmp][kk]=0;


    // double tmp=wallvecNorm[0];
    // double angle=((double) rand()/(RAND_MAX));
    // wallvecNorm[0]-=wallvecNorm[1]*20*(angle-0.5);
    // wallvecNorm[1]+=tmp*20*(angle-0.5);
    // tmp=std::sqrt(wallvecNorm[0]*wallvecNorm[0]+wallvecNorm[1]*wallvecNorm[1]);
    // wallvecNorm[0]/=tmp;
    // wallvecNorm[1]/=tmp;        
  

    cellData[numCellTmp][variableIndex(0,3)]  =wallvecNorm[0];
    cellData[numCellTmp][variableIndex(0,3)+1]=wallvecNorm[1];
    cellData[numCellTmp][variableIndex(0,3)+2]=1;
    cellData[numCellTmp][9]=1;
    cellData[numCellTmp][10]=1;



    std::cerr<<numCellTmp<<"  " <<cellData[numCellTmp][11] <<std::endl;

    //Change length of new wall between the divided daugther cells 
    //llData[numWallTmp][0] *= parameter(1);
  
 
    //if (i==2){
    //  std::cerr<< "cell "<< i << "auxin " << cellData[i][5]<<std::endl;
    //  std::cerr<< "cell "<< i+1 << "auxin " << cellData[i+1][5]<<" "<<wallData[numWallTmp][1]<<std::endl;
    //}
    //Check that the division did not mess up the data structure
    //T->checkConnectivity(1);	
  }

  VolumeViaLongestWallSpatial::
  VolumeViaLongestWallSpatial(std::vector<double> &paraValue, 
			      std::vector< std::vector<size_t> > 
			      &indValue ) {
    //
    // Do some checks on the parameters and variable indeces
    //
    if( paraValue.size()!=4 ) {
      std::cerr << "DivisionVolumeViaLongestWallSpatial::"
		<< "DivisionVolumeViaLongestWallSpatial() "
		<< "Four parameters used V_threshold, LWall_frac, "
		<< "Lwall_threshold, and spatial threshold." << std::endl;
      exit(0);
    }
    if( indValue.size() != 2 || indValue[0].size() != 1) {
      std::cerr << "DivisionVolumeViaLongestWallSpatial::"
		<< "DivisionVolumeViaLongestWallSpatial() "
		<< "Spatial index in first level and variable indices for volume dependent cell "
		<< "variables in second need to be provided." << std::endl;
      exit(0);
    }
    //
    // Set the variable values
    //
    setId("DivisionVolumeViaLongestWallSpatial");
    setNumChange(1);
    setParameter(paraValue);  
    setVariableIndex(indValue);
    //
    // Set the parameter identities
    //
    std::vector<std::string> tmp( numParameter() );
    tmp.resize( numParameter() );
    tmp[0] = "V_threshold";
    tmp[1] = "LWall_frac";
    tmp[2] = "LWall_threshold";
    tmp[3] = "spatial_threshold";
    setParameterId( tmp );
  }

  int VolumeViaLongestWallSpatial::
  flag(Tissue *T,size_t i,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) {
  
    size_t sI=variableIndex(0,0);
    assert( sI<vertexData[0].size());
    if (i==0) {//Calculate max position
      sMax_=vertexData[0][sI];
      size_t numV=vertexData.size();
      for (size_t i=1; i<numV; ++i)
	if (vertexData[i][sI]>sMax_)
	  sMax_=vertexData[i][sI];
    }
  
    std::vector<double> position=T->cell(i).positionFromVertex(vertexData);
    double sDistance = sMax_-position[sI];
    if( T->cell(i).calculateVolume(vertexData) > parameter(0) &&
	sDistance<parameter(3) ) {
      std::cerr << "Cell " << i << " marked for division with volume " 
		<< T->cell(i).volume() << std::endl;
      return 1;
    } 
    return 0;
  }

  void VolumeViaLongestWallSpatial::
  update(Tissue *T,size_t i,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDeriv,
	 DataMatrix &wallDeriv,
	 DataMatrix &vertexDeriv ) {
  
    Cell *divCell = &(T->cell(i));
    size_t dimension = vertexData[0].size();
    assert( divCell->numWall() > 1 );
    assert( dimension==2 || dimension==3 ); 
    //
    // Find longest wall
    // 
    size_t wI=0,w3I=divCell->numWall();
    double maxLength = divCell->wall(0)->setLengthFromVertexPosition(vertexData);
    for( size_t k=1 ; k<divCell->numWall() ; ++k ) {
      double tmpLength = divCell->wall(k)->setLengthFromVertexPosition(vertexData);
      if( tmpLength > maxLength ) {
	wI=k;
	maxLength = tmpLength;
      }
    }   
  
    //
    // Find position for first new vertex
    //
    std::vector<double> nW(dimension),nW2(dimension),v1Pos(dimension),
      v2Pos(dimension);
    size_t v1wI = divCell->wall(wI)->vertex1()->index();
    size_t v2wI = divCell->wall(wI)->vertex2()->index();
    for( size_t d=0 ; d<dimension ; ++d ) {
      nW[d] = (vertexData[v1wI][d]-vertexData[v2wI][d])/maxLength;
      v1Pos[d] = 0.5*(vertexData[v1wI][d]+vertexData[v2wI][d]);
    }
    //
    // Find intersection with another wall via vector perpendicular to first wall
    //
    if (dimension==2) {
      nW2[1] = nW[0];
      nW2[0] = -nW[1];
    }
    else if (dimension==3) {
      nW2[0]=nW[0];
      nW2[1]=nW[1];
      nW2[2]=nW[2];
    }
    if (findSecondDivisionWall(vertexData, divCell, wI, w3I, v1Pos, nW2, v2Pos)) {
      std::cerr << "DivisionVolumeViaLongestWallSpatial::update "
		<< "failed to find the second wall for division!" << std::endl;
      exit(-1);
    }
    //
    // Do the division (add one cell, three walls, and two vertices)
    //
    size_t numWallTmp=wallData.size();
    assert( numWallTmp==T->numWall() );
    //Divide
    T->divideCell(divCell,wI,w3I,v1Pos,v2Pos,cellData,wallData,vertexData,
		  cellDeriv,wallDeriv,vertexDeriv,variableIndex(1),
		  parameter(2));
    assert( numWallTmp+3 == T->numWall() );
  
    //Change length of new wall between the divided daugther cells 
    wallData[numWallTmp][0] *= parameter(1);
  
    //Check that the division did not mess up the data structure
    //T->checkConnectivity(1);	
  }

  VolumeViaLongestWall3D::
  VolumeViaLongestWall3D(std::vector<double> &paraValue, 
			 std::vector< std::vector<size_t> > 
			 &indValue ) {
  
    //Do some checks on the parameters and variable indeces
    //////////////////////////////////////////////////////////////////////
    if( paraValue.size()!=3 ) {
      std::cerr << "DivisionVolumeViaLongestWall3D::"
		<< "DivisionVolumeViaLongestWall3D() "
		<< "Two parameters used V_threshold, LWall_frac LWall_threshold\n";
      exit(0);
    }
    if( indValue.size() != 1 ) {
      std::cerr << "DivisionVolumeViaLongestWall3D::"
		<< "DivisionVolumeViaLongestWall3D() "
		<< "Variable indices for volume dependent cell "
		<< "variables is used." << std::endl;
      exit(0);
    }
    //Set the variable values
    //////////////////////////////////////////////////////////////////////
    setId("DivisionVolumeViaLongestWall3D");
    setNumChange(1);
    setParameter(paraValue);  
    setVariableIndex(indValue);
  
    //Set the parameter identities
    //////////////////////////////////////////////////////////////////////
    std::vector<std::string> tmp( numParameter() );
    tmp.resize( numParameter() );
    tmp[0] = "V_threshold";
    tmp[1] = "LWall_frac";
    tmp[2] = "LWall_threshold";
    setParameterId( tmp );
  }

  //! Flags a cell for division if the volume above threshold
  /*! 
   */
  int VolumeViaLongestWall3D::
  flag(Tissue *T,size_t i,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) {
  
    if( T->cell(i).calculateVolume(vertexData) > parameter(0) ) {
      std::cerr << "Cell " << i << " marked for division with volume " 
		<< T->cell(i).volume() << std::endl;
      return 1;
    } 
    return 0;
  }

  //! Updates the dividing cell by adding a prependicular wall from the longest
  /*! 
   */
  void VolumeViaLongestWall3D::
  update(Tissue *T,size_t i,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDeriv,
	 DataMatrix &wallDeriv,
	 DataMatrix &vertexDeriv ) {
  
    Cell *divCell = &(T->cell(i));
    size_t dimension = vertexData[0].size();
    assert( divCell->numWall() > 1 );
    assert( dimension==3 );
  
    //Find longest wall
    //////////////////////////////////////////////////////////////////////
    size_t wI=0;
    double maxLength = divCell->wall(0)->setLengthFromVertexPosition(vertexData);
    for( size_t k=1 ; k<divCell->numWall() ; ++k ) {
      double tmpLength = divCell->wall(k)->setLengthFromVertexPosition(vertexData);
      if( tmpLength > maxLength ) {
	wI=k;
	maxLength = tmpLength;
      }
    }   
  
    std::vector<double> nW(dimension),nW2(dimension),v1Pos(dimension),
      v2Pos(dimension);
    size_t v1wI = divCell->wall(wI)->vertex1()->index();
    size_t v2wI = divCell->wall(wI)->vertex2()->index();
  
    for( size_t d=0 ; d<dimension ; ++d ) {
      nW[d] = (vertexData[v1wI][d]-vertexData[v2wI][d])/maxLength;
      v1Pos[d] = 0.5*(vertexData[v1wI][d]+vertexData[v2wI][d]);
    }
  
    //Find intersection with another wall by looking at perpendicular plane
    //////////////////////////////////////////////////////////////////////
    size_t w3I=divCell->numWall();
    //double minDist,w3s;
    std::vector<size_t> w3Tmp;
    std::vector<double> w3tTmp;
    int flag=0;
    for( size_t k=0 ; k<divCell->numWall() ; ++k ) {
      if( k!=wI ) {
	size_t v1w3Itmp = divCell->wall(k)->vertex1()->index();
	size_t v2w3Itmp = divCell->wall(k)->vertex2()->index();
	std::vector<double> w3(dimension),w0(dimension);
	double fac1=0.0,fac2=0.0;
	for( size_t d=0 ; d<dimension ; ++d ) {
	  w3[d] = vertexData[v2w3Itmp][d]-vertexData[v1w3Itmp][d];
	  fac1 += nW[d]*(v1Pos[d]-vertexData[v1w3Itmp][d]);
	  fac2 += nW[d]*w3[d]; 
	}
	if( fac2 != 0.0 ) {//else parallell and not applicable
	  double t = fac1/fac2;
	  if( t>=0.0 && t<1.0 ) {//within wall
	    std::cerr << "wall " << k << " (t=" << t << " " << fac1 << "/"
		      << fac2 
		      << ") chosen as second wall"
		      << std::endl;					
	    flag++;
	    w3I = k;
	    w3Tmp.push_back(k);
	    w3tTmp.push_back(t);
	  }
	}
      }
    }
    assert( w3I != wI );
    assert( w3I != divCell->numWall() );
    if( flag != 1 ) {
      std::cerr << "divideVolumeViaLongestWall3D::update() Warning"
		<< " more than one wall possible as connection "
		<< "for cell " 
		<< i << std::endl; 
      for( size_t k=0 ; k<divCell->numWall() ; ++k ) {
	std::cerr << "0 " 
		  << vertexData[divCell->wall(k)->vertex1()->index()][0]
		  << " " 
		  << vertexData[divCell->wall(k)->vertex1()->index()][1]
		  << " " 
		  << vertexData[divCell->wall(k)->vertex1()->index()][2]
		  << "\n0 " 
		  << vertexData[divCell->wall(k)->vertex2()->index()][0]
		  << " " 
		  << vertexData[divCell->wall(k)->vertex2()->index()][1]
		  << " " 
		  << vertexData[divCell->wall(k)->vertex2()->index()][2]
		  << "\n\n\n";
      }
      for( size_t kk=0 ; kk<w3Tmp.size() ; ++kk ) {
	size_t k = w3Tmp[kk];
	std::cerr << "1 " 
		  << vertexData[divCell->wall(k)->vertex1()->index()][0]
		  << " " 
		  << vertexData[divCell->wall(k)->vertex1()->index()][1]
		  << " " 
		  << vertexData[divCell->wall(k)->vertex1()->index()][2]
		  << "\n1 " 
		  << vertexData[divCell->wall(k)->vertex2()->index()][0]
		  << " " 
		  << vertexData[divCell->wall(k)->vertex2()->index()][1]
		  << " " 
		  << vertexData[divCell->wall(k)->vertex2()->index()][2]
		  << "\n\n\n";
      }
      std::cerr << "2 " 
		<< vertexData[divCell->wall(wI)->vertex1()->index()][0]
		<< " " 
		<< vertexData[divCell->wall(wI)->vertex1()->index()][1]
		<< " " 
		<< vertexData[divCell->wall(wI)->vertex1()->index()][2]
		<< "\n2 " 
		<< vertexData[divCell->wall(wI)->vertex2()->index()][0]
		<< " " 
		<< vertexData[divCell->wall(wI)->vertex2()->index()][1]
		<< " " 
		<< vertexData[divCell->wall(wI)->vertex2()->index()][2]
		<< "\n\n\n";
      std::cerr << "3 " 
		<< vertexData[divCell->wall(w3I)->vertex1()->index()][0]
		<< " " 
		<< vertexData[divCell->wall(w3I)->vertex1()->index()][1]
		<< " " 
		<< vertexData[divCell->wall(w3I)->vertex1()->index()][2]
		<< "\n3 " 
		<< vertexData[divCell->wall(w3I)->vertex2()->index()][0]
		<< " " 
		<< vertexData[divCell->wall(w3I)->vertex2()->index()][1]
		<< " " 
		<< vertexData[divCell->wall(w3I)->vertex2()->index()][2]
		<< "\n\n\n";
      std::cerr << "4 " 
		<< 0.5*(vertexData[divCell->wall(wI)->vertex1()->index()][0]+
			vertexData[divCell->wall(wI)->vertex2()->index()][0])
		<< " " 
		<< 0.5*(vertexData[divCell->wall(wI)->vertex1()->index()][1]+
			vertexData[divCell->wall(wI)->vertex2()->index()][1])
		<< " " 
		<< 0.5*(vertexData[divCell->wall(wI)->vertex1()->index()][2]+
			vertexData[divCell->wall(wI)->vertex2()->index()][2])
		<< "\n4 "
		<< 0.5*(vertexData[divCell->wall(w3I)->vertex1()->index()][0]+
			vertexData[divCell->wall(w3I)->vertex2()->index()][0])
		<< " " 
		<< 0.5*(vertexData[divCell->wall(w3I)->vertex1()->index()][1]+
			vertexData[divCell->wall(w3I)->vertex2()->index()][1])
		<< " " 
		<< 0.5*(vertexData[divCell->wall(w3I)->vertex1()->index()][2]+
			vertexData[divCell->wall(w3I)->vertex2()->index()][2])
		<< "\n\n\n";
      exit(-1);
    }
  
    size_t v1w3I = divCell->wall(w3I)->vertex1()->index();
    size_t v2w3I = divCell->wall(w3I)->vertex2()->index();
    for( size_t d=0 ; d<dimension ; ++d )
      v2Pos[d] = vertexData[v1w3I][d] + 
	w3tTmp[w3tTmp.size()-1]*(vertexData[v2w3I][d]-vertexData[v1w3I][d]);
    //v2Pos[d] = 0.5*(vertexData[v1w3I][d]+vertexData[v2w3I][d]);
  
    //Add one cell, three walls, and two vertices
    //////////////////////////////////////////////////////////////////////
    //Save number of walls
    size_t numWallTmp=wallData.size();
    assert( numWallTmp==T->numWall() );
    //Divide
    T->divideCell(divCell,wI,w3I,v1Pos,v2Pos,cellData,wallData,vertexData,
		  cellDeriv,wallDeriv,vertexDeriv,variableIndex(0),
		  parameter(2));
    assert( numWallTmp+3 == T->numWall() );
  
    //Change length of new wall between the divided daugther cells 
    wallData[numWallTmp][0] *= parameter(1);
  
    //Check that the division did not messed up the data structure
    //T->checkConnectivity(1);	
  }

  VolumeViaLongestWall3DSpatial::
  VolumeViaLongestWall3DSpatial(std::vector<double> &paraValue, 
				std::vector< std::vector<size_t> > 
				&indValue ) {
  
    //Do some checks on the parameters and variable indeces
    //////////////////////////////////////////////////////////////////////
    if( paraValue.size()!=4 ) {
      std::cerr << "DivisionVolumeViaLongestWall3DSpatial::"
		<< "DivisionVolumeViaLongestWall3DSpatial() "
		<< "Four parameters used V_threshold, LWall_frac LWall_threshold "
		<< "and spatial_threshold." << std::endl;
      exit(0);
    }
    if( indValue.size() != 2 || indValue[0].size() != 1 ) {
      std::cerr << "DivisionVolumeViaLongestWall3DSpatial::"
		<< "DivisionVolumeViaLongestWall3DSpatial() "
		<< "Variable indices for spatial direction in first and volume dependent cell "
		<< "variables in second are needed." << std::endl;
      exit(0);
    }
    //Set the variable values
    //////////////////////////////////////////////////////////////////////
    setId("DivisionVolumeViaLongestWall3DSpatial");
    setNumChange(1);
    setParameter(paraValue);  
    setVariableIndex(indValue);
  
    //Set the parameter identities
    //////////////////////////////////////////////////////////////////////
    std::vector<std::string> tmp( numParameter() );
    tmp.resize( numParameter() );
    tmp[0] = "V_threshold";
    tmp[1] = "LWall_frac";
    tmp[2] = "LWall_threshold";
    tmp[3] = "spatial_threshold";
    setParameterId( tmp );
  }

  int VolumeViaLongestWall3DSpatial::
  flag(Tissue *T,size_t i,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) {
  
    size_t sI=variableIndex(0,0);
    assert( sI<vertexData[0].size());
    if (i==0) {//Calculate max position
      sMax_=vertexData[0][sI];
      size_t numV=vertexData.size();
      for (size_t i=1; i<numV; ++i)
	if (vertexData[i][sI]>sMax_)
	  sMax_=vertexData[i][sI];
    }
  
    std::vector<double> position=T->cell(i).positionFromVertex(vertexData);
    double sDistance = sMax_-position[sI];
    if( T->cell(i).calculateVolume(vertexData) > parameter(0) &&
	sDistance<parameter(3) ) {
      std::cerr << "Cell " << i << " marked for division with volume " 
		<< T->cell(i).volume() << std::endl;
      return 1;
    } 
    return 0;
  }

  void VolumeViaLongestWall3DSpatial::
  update(Tissue *T,size_t i,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDeriv,
	 DataMatrix &wallDeriv,
	 DataMatrix &vertexDeriv ) {
  
    Cell *divCell = &(T->cell(i));
    size_t dimension = vertexData[0].size();
    assert( divCell->numWall() > 1 );
    assert( dimension==3 );
  
    //Find longest wall
    //////////////////////////////////////////////////////////////////////
    size_t wI=0;
    double maxLength = divCell->wall(0)->setLengthFromVertexPosition(vertexData);
    for( size_t k=1 ; k<divCell->numWall() ; ++k ) {
      double tmpLength = divCell->wall(k)->setLengthFromVertexPosition(vertexData);
      if( tmpLength > maxLength ) {
	wI=k;
	maxLength = tmpLength;
      }
    }   
  
    std::vector<double> nW(dimension),nW2(dimension),v1Pos(dimension),
      v2Pos(dimension);
    size_t v1wI = divCell->wall(wI)->vertex1()->index();
    size_t v2wI = divCell->wall(wI)->vertex2()->index();
  
    for( size_t d=0 ; d<dimension ; ++d ) {
      nW[d] = (vertexData[v1wI][d]-vertexData[v2wI][d])/maxLength;
      v1Pos[d] = 0.5*(vertexData[v1wI][d]+vertexData[v2wI][d]);
    }
  
    //Find intersection with another wall by looking at perpendicular plane
    //////////////////////////////////////////////////////////////////////
    size_t w3I=divCell->numWall();
    //double minDist,w3s;
    std::vector<size_t> w3Tmp;
    std::vector<double> w3tTmp;
    int flag=0;
    for( size_t k=0 ; k<divCell->numWall() ; ++k ) {
      if( k!=wI ) {
	size_t v1w3Itmp = divCell->wall(k)->vertex1()->index();
	size_t v2w3Itmp = divCell->wall(k)->vertex2()->index();
	std::vector<double> w3(dimension),w0(dimension);
	double fac1=0.0,fac2=0.0;
	for( size_t d=0 ; d<dimension ; ++d ) {
	  w3[d] = vertexData[v2w3Itmp][d]-vertexData[v1w3Itmp][d];
	  fac1 += nW[d]*(v1Pos[d]-vertexData[v1w3Itmp][d]);
	  fac2 += nW[d]*w3[d]; 
	}
	if( fac2 != 0.0 ) {//else parallell and not applicable
	  double t = fac1/fac2;
	  if( t>=0.0 && t<1.0 ) {//within wall
	    std::cerr << "wall " << k << " (t=" << t << " " << fac1 << "/"
		      << fac2 
		      << ") chosen as second wall"
		      << std::endl;					
	    flag++;
	    w3I = k;
	    w3Tmp.push_back(k);
	    w3tTmp.push_back(t);
	  }
	}
      }
    }
    assert( w3I != wI );
    assert( w3I != divCell->numWall() );
    if( flag != 1 ) {
      std::cerr << "divideVolumeViaLongestWall3DSpatial::update() Warning"
		<< " more than one wall possible as connection "
		<< "for cell " 
		<< i << std::endl; 
      for( size_t k=0 ; k<divCell->numWall() ; ++k ) {
	std::cerr << "0 " 
		  << vertexData[divCell->wall(k)->vertex1()->index()][0]
		  << " " 
		  << vertexData[divCell->wall(k)->vertex1()->index()][1]
		  << " " 
		  << vertexData[divCell->wall(k)->vertex1()->index()][2]
		  << "\n0 " 
		  << vertexData[divCell->wall(k)->vertex2()->index()][0]
		  << " " 
		  << vertexData[divCell->wall(k)->vertex2()->index()][1]
		  << " " 
		  << vertexData[divCell->wall(k)->vertex2()->index()][2]
		  << "\n\n\n";
      }
      for( size_t kk=0 ; kk<w3Tmp.size() ; ++kk ) {
	size_t k = w3Tmp[kk];
	std::cerr << "1 " 
		  << vertexData[divCell->wall(k)->vertex1()->index()][0]
		  << " " 
		  << vertexData[divCell->wall(k)->vertex1()->index()][1]
		  << " " 
		  << vertexData[divCell->wall(k)->vertex1()->index()][2]
		  << "\n1 " 
		  << vertexData[divCell->wall(k)->vertex2()->index()][0]
		  << " " 
		  << vertexData[divCell->wall(k)->vertex2()->index()][1]
		  << " " 
		  << vertexData[divCell->wall(k)->vertex2()->index()][2]
		  << "\n\n\n";
      }
      std::cerr << "2 " 
		<< vertexData[divCell->wall(wI)->vertex1()->index()][0]
		<< " " 
		<< vertexData[divCell->wall(wI)->vertex1()->index()][1]
		<< " " 
		<< vertexData[divCell->wall(wI)->vertex1()->index()][2]
		<< "\n2 " 
		<< vertexData[divCell->wall(wI)->vertex2()->index()][0]
		<< " " 
		<< vertexData[divCell->wall(wI)->vertex2()->index()][1]
		<< " " 
		<< vertexData[divCell->wall(wI)->vertex2()->index()][2]
		<< "\n\n\n";
      std::cerr << "3 " 
		<< vertexData[divCell->wall(w3I)->vertex1()->index()][0]
		<< " " 
		<< vertexData[divCell->wall(w3I)->vertex1()->index()][1]
		<< " " 
		<< vertexData[divCell->wall(w3I)->vertex1()->index()][2]
		<< "\n3 " 
		<< vertexData[divCell->wall(w3I)->vertex2()->index()][0]
		<< " " 
		<< vertexData[divCell->wall(w3I)->vertex2()->index()][1]
		<< " " 
		<< vertexData[divCell->wall(w3I)->vertex2()->index()][2]
		<< "\n\n\n";
      std::cerr << "4 " 
		<< 0.5*(vertexData[divCell->wall(wI)->vertex1()->index()][0]+
			vertexData[divCell->wall(wI)->vertex2()->index()][0])
		<< " " 
		<< 0.5*(vertexData[divCell->wall(wI)->vertex1()->index()][1]+
			vertexData[divCell->wall(wI)->vertex2()->index()][1])
		<< " " 
		<< 0.5*(vertexData[divCell->wall(wI)->vertex1()->index()][2]+
			vertexData[divCell->wall(wI)->vertex2()->index()][2])
		<< "\n4 "
		<< 0.5*(vertexData[divCell->wall(w3I)->vertex1()->index()][0]+
			vertexData[divCell->wall(w3I)->vertex2()->index()][0])
		<< " " 
		<< 0.5*(vertexData[divCell->wall(w3I)->vertex1()->index()][1]+
			vertexData[divCell->wall(w3I)->vertex2()->index()][1])
		<< " " 
		<< 0.5*(vertexData[divCell->wall(w3I)->vertex1()->index()][2]+
			vertexData[divCell->wall(w3I)->vertex2()->index()][2])
		<< "\n\n\n";
      exit(-1);
    }
  
    size_t v1w3I = divCell->wall(w3I)->vertex1()->index();
    size_t v2w3I = divCell->wall(w3I)->vertex2()->index();
    for( size_t d=0 ; d<dimension ; ++d )
      v2Pos[d] = vertexData[v1w3I][d] + 
	w3tTmp[w3tTmp.size()-1]*(vertexData[v2w3I][d]-vertexData[v1w3I][d]);
    //v2Pos[d] = 0.5*(vertexData[v1w3I][d]+vertexData[v2w3I][d]);
  
    //Add one cell, three walls, and two vertices
    //////////////////////////////////////////////////////////////////////
    //Save number of walls
    size_t numWallTmp=wallData.size();
    assert( numWallTmp==T->numWall() );
    //Divide
    T->divideCell(divCell,wI,w3I,v1Pos,v2Pos,cellData,wallData,vertexData,
		  cellDeriv,wallDeriv,vertexDeriv,variableIndex(1),
		  parameter(2));
    assert( numWallTmp+3 == T->numWall() );
  
    //Change length of new wall between the divided daugther cells 
    wallData[numWallTmp][0] *= parameter(1);
  
    //Check that the division did not messed up the data structure
    //T->checkConnectivity(1);	
  }

  //!Constructor
  VolumeViaStrain::
  VolumeViaStrain(std::vector<double> &paraValue, 
		  std::vector< std::vector<size_t> > 
		  &indValue ) {
  
    //Do some checks on the parameters and variable indeces
    //////////////////////////////////////////////////////////////////////
    if( paraValue.size()!=4 ) {
      std::cerr << "DivisionVolumeViaStrain::"
		<< "DivisionVolumeViaStrain() "
		<< "Four parameters used V_threshold, LWall_frac, "
		<< "Lwall_threshold, and Parallell_flag" << std::endl;
      exit(0);
    }
    if( indValue.size() != 1 ) {
      std::cerr << "DivisionVolumeViaStrain::"
		<< "DivisionVolumeViaStrain() "
		<< "Variable indices for volume dependent cell "
		<< "variables is used.\n";
      exit(0);
    }
    //Set the variable values
    //////////////////////////////////////////////////////////////////////
    setId("DivisionVolumeViaStrain");
    setNumChange(1);
    setParameter(paraValue);  
    setVariableIndex(indValue);
  
    //Set the parameter identities
    //////////////////////////////////////////////////////////////////////
    std::vector<std::string> tmp( numParameter() );
    tmp.resize( numParameter() );
    tmp[0] = "V_threshold";
    tmp[1] = "LWall_frac";
    tmp[2] = "LWall_threshold";
    tmp[3] = "Parallell_flag";
    setParameterId( tmp );
  }

  //! Flags a cell for division if the volume above threshold
  /*! 
   */
  int VolumeViaStrain::
  flag(Tissue *T,size_t i,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) {
  
    if( T->cell(i).calculateVolume(vertexData) > parameter(0) ) {
      std::cerr << "Cell " << i << " marked for division with volume " 
		<< T->cell(i).volume() << std::endl;
      return 1;
    } 
    return 0;
  }

  //! Updates the dividing cell by adding a prependicular wall from the longest
  /*! 
   */
  void VolumeViaStrain::
  update(Tissue *T,size_t cellI,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDeriv,
	 DataMatrix &wallDeriv,
	 DataMatrix &vertexDeriv ) {
  
    Cell *divCell = &(T->cell(cellI));
    size_t dimension = vertexData[0].size();
    assert( divCell->numWall() > 2 );
    assert( dimension==2 );
  
    //////////////////////////////////////////////////////////////////////
    //Calculate strain directions and print walls and strain vectors
    //by using x,x+dt*dx/dt as two points
    //////////////////////////////////////////////////////////////////////
  
    T->derivs(cellData,wallData,vertexData,cellDeriv,wallDeriv,vertexDeriv);
  
    //Create temporary x,y,dx positions
    size_t numV = divCell->numVertex(); 
    DataMatrix x(numV),y(numV),dx(numV),
      xM(numV),yM(numV),dxM(numV);
  
    double dt=1.0;
    std::vector<double> xMean(dimension),yMean(dimension),dxMean(dimension);
  
    for( size_t i=0 ; i<numV ; ++i ) {
      size_t vI = divCell->vertex(i)->index();
      x[i] = vertexData[vI];
      dx[i] = vertexDeriv[vI];
      std::vector<double> tmp(dimension);
      tmp[0] = x[i][0]+dt*dx[i][0];
      tmp[1] = x[i][1]+dt*dx[i][1];
      y[i] = tmp;
      xMean[0] += x[i][0];
      xMean[1] += x[i][1];
      yMean[0] += y[i][0];
      yMean[1] += y[i][1];
      dxMean[0] += dx[i][0];			
      dxMean[1] += dx[i][1];			
    }
    xMean[0] /=numV;
    xMean[1] /=numV;
    yMean[0] /=numV;
    yMean[1] /=numV;
    dxMean[0] /=numV;
    dxMean[1] /=numV;
    for( size_t i=0 ; i<numV ; ++i ) {
      xM[i].resize(dimension);
      xM[i][0] =x[i][0]-xMean[0];
      xM[i][1] =x[i][1]-xMean[1];
      yM[i].resize(dimension);
      yM[i][0] =y[i][0]-yMean[0];
      yM[i][1] =y[i][1]-yMean[1];
      dxM[i].resize(dimension);
      dxM[i][0] =dx[i][0]-dxMean[0];
      dxM[i][1] =dx[i][1]-dxMean[1];
    }
  
    //Calculate A = (x^t x)^{-1} (x^t y)
    std::vector<std::vector<double> > xTx(dimension),xTy(dimension),
      xTxM(dimension),A(dimension);
    for( size_t i=0 ; i<dimension ; ++i ) {
      xTx[i].resize(dimension);
      xTy[i].resize(dimension);
      xTxM[i].resize(dimension);
      A[i].resize(dimension);
    }
  
    for( size_t i=0 ; i<dimension ; ++i ) {
      for( size_t j=0 ; j<dimension ; ++j ) {
	for( size_t v=0 ; v<numV ; ++v ) {
	  xTx[i][j] += xM[v][i]*xM[v][j];
	  xTy[i][j] += xM[v][i]*yM[v][j];
	}
      }
    }
    double detM = xTx[0][0]*xTx[1][1]-xTx[0][1]*xTx[1][0];
    detM = 1.0/detM;
    xTxM[0][0] = detM*xTx[1][1];
    xTxM[1][1] = detM*xTx[0][0];
    xTxM[0][1] = -detM*xTx[1][0];
    xTxM[1][0] = -detM*xTx[0][1];
  
    //Calculate A
    A[0][0] = xTxM[0][0]*xTy[0][0] + xTxM[0][1]*xTy[1][0];
    A[0][1] = xTxM[0][0]*xTy[0][1] + xTxM[0][1]*xTy[1][1];
    A[1][0] = xTxM[1][0]*xTy[0][0] + xTxM[1][1]*xTy[1][0];
    A[1][1] = xTxM[1][0]*xTy[0][1] + xTxM[1][1]*xTy[1][1];
  
    //Apply SVD to A
    //////////////////////////////////////////////////////////////////////
  
    //Make sure determinant is non-zero
    double detA = A[0][0]*A[1][1] - A[0][1]*A[1][0];
    if( detA==0 ) {
      std::cerr << "DivisionVolumeViaStrain::update() Determinant zero\n";
      exit(-1);
    }
    //double t = std::sqrt( (A[0][0]+A[1][1])*(A[0][0]+A[1][1]) +
    //(A[0][1]-A[1][0])*(A[0][1]-A[1][0]) );
    //double w = std::sqrt( (A[0][0]-A[1][1])*(A[0][0]-A[1][1]) +
    //(A[0][1]+A[1][0])*(A[0][1]+A[1][0]) );
    double tau = std::atan2( A[0][0]-A[1][1],A[0][1]+A[1][0] );
    double omega = std::atan2( A[0][0]+A[1][1],A[0][1]-A[1][0] );
  
    //double p = 0.5*(t+w);
    //double q = 0.5*(t-w);
    double theta = 0.5*(tau-omega);
    //double phi = 0.5*(tau+omega);
  
    //std::cerr << "Cell: " << divCell->index() << "\n";
    //std::cerr << "A:\n";
    //std::cerr << A[0][0] << " " << A[0][1] << "\n";
    //std::cerr << A[1][0] << " " << A[1][1] << "\n";
    //std::cerr << p << " " << q << "  " << t << " " << w << "\n";
    //std::cout << divCell->index() << " " << p << " " << q << " " << theta
    //					<< " " << phi << "\n";
  
    // Find walls and vertex positions needed for the division
    //////////////////////////////////////////////////////////////////////
  
    //Create direction for new wall (xMean+t*n)
    std::vector<double> n(dimension);
    //rotate angle v=theta-90
    double v = theta;
    if( parameter(3) != 1.0 )
      v = theta - 0.5*3.14159;
    n[0]=std::cos(v);
    n[1]=std::sin(v);
  
    //Find two (and two only) intersecting walls
    //////////////////////////////////////////////////////////////////////
    std::vector<size_t> wI(2);
    std::vector<double> s(2);
    wI[0]=0;
    wI[1]=divCell->numWall();
    s[0]=s[1]=-1.0;
    //double minDist,w3s;
    std::vector<size_t> w3Tmp;
    std::vector<double> w3tTmp;
    int flag=0;
    for( size_t k=0 ; k<divCell->numWall() ; ++k ) {
      size_t v1Tmp = divCell->wall(k)->vertex1()->index();
      size_t v2Tmp = divCell->wall(k)->vertex2()->index();
      std::vector<double> w3(dimension),w0(dimension);
      for( size_t d=0 ; d<dimension ; ++d ) {
	w3[d] = vertexData[v2Tmp][d]-vertexData[v1Tmp][d];
	w0[d] = xMean[d]-vertexData[v1Tmp][d];
      }
      double a=0.0,b=0.0,c=0.0,d=0.0,e=0.0;//a=1.0
      for( size_t dim=0 ; dim<dimension ; ++dim ) {
	a += n[dim]*n[dim];
	b += n[dim]*w3[dim];
	c += w3[dim]*w3[dim];
	d += n[dim]*w0[dim];
	e += w3[dim]*w0[dim];
      }
      double fac=a*c-b*b;//a*c-b*b
      if( fac>0.0 ) {//else parallell and not applicable
	fac = 1.0/fac;
	//double s = fac*(b*e-c*d);
	double t = fac*(a*e-b*d);//fac*(a*e-b*d)
	if( t>=0.0 && t<1.0 ) {//within wall
	  //double dx0 = w0[0] +fac*((b*e-c*d)*nW2[0]+()*w3[0]); 					
	  w3Tmp.push_back(k);
	  w3tTmp.push_back(t);
	  std::cerr << "Dividing cell " << divCell->index() << " via wall "
		    << k << " at t=" << t << std::endl;
	  if( flag<2 ) {
	    s[flag] = t;
	    wI[flag] = k;
	  }				
	  flag++;
	}
      }
    }
    assert( wI[1] != divCell->numWall() && wI[0] != wI[1] );
    if( flag != 2 ) {
      std::cerr << "divideVolumeVisStrain::update Warning"
		<< " not two walls possible as connection "
		<< "for cell " 
		<< cellI << std::endl; 
      for( size_t k=0 ; k<divCell->numWall() ; ++k ) {
	std::cerr << "0 " 
		  << vertexData[divCell->wall(k)->vertex1()->index()][0]
		  << " " 
		  << vertexData[divCell->wall(k)->vertex1()->index()][1]
		  << "\n0 " 
		  << vertexData[divCell->wall(k)->vertex2()->index()][0]
		  << " " 
		  << vertexData[divCell->wall(k)->vertex2()->index()][1]
		  << "\n\n\n";
      }
      for( size_t kk=0 ; kk<w3Tmp.size() ; ++kk ) {
	size_t k = w3Tmp[kk];
	std::cerr << "1 " 
		  << vertexData[divCell->wall(k)->vertex1()->index()][0]
		  << " " 
		  << vertexData[divCell->wall(k)->vertex1()->index()][1]
		  << "\n1 " 
		  << vertexData[divCell->wall(k)->vertex2()->index()][0]
		  << " " 
		  << vertexData[divCell->wall(k)->vertex2()->index()][1]
		  << "\n\n\n";
      }
      std::cerr << "2 " 
		<< vertexData[divCell->wall(wI[0])->vertex1()->index()][0]
		<< " " 
		<< vertexData[divCell->wall(wI[0])->vertex1()->index()][1]
		<< "\n2 " 
		<< vertexData[divCell->wall(wI[0])->vertex2()->index()][0]
		<< " " 
		<< vertexData[divCell->wall(wI[0])->vertex2()->index()][1]
		<< "\n\n\n";
      std::cerr << "3 " 
		<< vertexData[divCell->wall(wI[1])->vertex1()->index()][0]
		<< " " 
		<< vertexData[divCell->wall(wI[1])->vertex1()->index()][1]
		<< "\n3 " 
		<< vertexData[divCell->wall(wI[1])->vertex2()->index()][0]
		<< " " 
		<< vertexData[divCell->wall(wI[1])->vertex2()->index()][1]
		<< "\n\n\n";
      std::cerr << "4 " 
		<< 0.5*(vertexData[divCell->wall(wI[0])->vertex1()->index()][0]+
			vertexData[divCell->wall(wI[0])->vertex2()->index()][0])
		<< " " 
		<< 0.5*(vertexData[divCell->wall(wI[0])->vertex1()->index()][1]+
			vertexData[divCell->wall(wI[0])->vertex2()->index()][1])
		<< "\n4 "
		<< 0.5*(vertexData[divCell->wall(wI[1])->vertex1()->index()][0]+
			vertexData[divCell->wall(wI[1])->vertex2()->index()][0])
		<< " " 
		<< 0.5*(vertexData[divCell->wall(wI[1])->vertex1()->index()][1]+
			vertexData[divCell->wall(wI[1])->vertex2()->index()][1])
		<< "\n\n\n";
      exit(-1);
    }	
    //Addition of new vertices at walls at position 's' 
    std::vector<double> v1Pos(dimension),v2Pos(dimension);
    size_t v1I = divCell->wall(wI[0])->vertex1()->index();
    size_t v2I = divCell->wall(wI[0])->vertex2()->index();
    for( size_t d=0 ; d<dimension ; ++d )
      v1Pos[d] = vertexData[v1I][d]+ s[0]*(vertexData[v2I][d]-vertexData[v1I][d]);
    v1I = divCell->wall(wI[1])->vertex1()->index();
    v2I = divCell->wall(wI[1])->vertex2()->index();
    for( size_t d=0 ; d<dimension ; ++d )
      v2Pos[d] = vertexData[v1I][d]+s[1]*(vertexData[v2I][d]-vertexData[v1I][d]);
  
    //Add one cell, three walls, and two vertices
    //////////////////////////////////////////////////////////////////////
    //Save number of walls
    size_t numWallTmp=wallData.size();
    assert( numWallTmp==T->numWall() );
    //Divide
    T->divideCell(divCell,wI[0],wI[1],v1Pos,v2Pos,cellData,wallData,vertexData,
		  cellDeriv,wallDeriv,vertexDeriv,variableIndex(0),
		  parameter(2));
    assert( numWallTmp+3 == T->numWall() );
  
    //Change length of new wall between the divided daugther cells 
    wallData[numWallTmp][0] *= parameter(1);
  
    //Check that the division did not mess up the data structure
    //T->checkConnectivity(1);		
  }

  //!Constructor
  VolumeViaDirection::
  VolumeViaDirection(std::vector<double> &paraValue, 
		     std::vector< std::vector<size_t> > 
		     &indValue ) {
  
    //Do some checks on the parameters and variable indeces
    //////////////////////////////////////////////////////////////////////
    if (paraValue.size() != 5) {
      std::cerr << "DivisionVolumeViaDirection::"
		<< "DivisionVolumeViaDirection() "
		<< "Four parameters used V_threshold, LWall_frac, "
		<< "Lwall_threshold, Parallell_flag, and COM (1 = COM, 0 = Random)." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if( indValue.size() != 2 || indValue[0].size() != 1 ) {
      std::cerr << "DivisionVolumeViaDirection::"
		<< "DivisionVolumeViaDirection() "
		<< "Variable index for cell direction start index in first"
		<< "level and indices for volume dependent cell "
		<< "variables is at second level.\n";
      exit(0);
    }
    //Set the variable values
    //////////////////////////////////////////////////////////////////////
    setId("DivisionVolumeViaDirection");
    setNumChange(1);
    setParameter(paraValue);  
    setVariableIndex(indValue);
  
    //Set the parameter identities
    //////////////////////////////////////////////////////////////////////
    std::vector<std::string> tmp( numParameter() );
    tmp.resize( numParameter() );
    tmp[0] = "V_threshold";
    tmp[1] = "LWall_frac";
    tmp[2] = "LWall_threshold";
    tmp[3] = "Parallell_flag";
    tmp[4] = "COM_flag";
    setParameterId( tmp );
  }

  //! Flags a cell for division if the volume above threshold
  /*! 
   */
  int VolumeViaDirection::
  flag(Tissue *T,size_t i,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) {
  
    if( T->cell(i).calculateVolume(vertexData) > parameter(0) ) {
      std::cerr << "Cell " << i << " marked for division with volume " 
		<< T->cell(i).volume() << std::endl;
      return 1;
    } 
    return 0;
  }

  //! Updates the dividing cell by adding a prependicular wall from the defined direction
  /*! 
   */
  void VolumeViaDirection::
  update(Tissue *T,size_t cellI,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDeriv,
	 DataMatrix &wallDeriv,
	 DataMatrix &vertexDeriv ) {
  
    Cell *divCell = &(T->cell(cellI));
    size_t dimension = vertexData[0].size();
    assert( divCell->numWall() > 2 );
    assert( dimension==2 || dimension==3);
  
    std::vector<double> com(dimension);
  
    if (parameter(4) == 1)
      {
	com = divCell->positionFromVertex(vertexData);
      }
    else
      {
	try
	  {
	    com = divCell->randomPositionInCell(vertexData);
	  }
	catch (Cell::FailedToFindRandomPositionInCellException)
	  {
	    return;
	  }
      }
  
    // Get direction from cell variable (or random if not present
    std::vector<double> n(dimension);
    if( cellData[cellI][variableIndex(0,0)+dimension] ) {
      if (dimension==2) {
	if( parameter(3) != 1.0 ) {//Perpendicular to given direction
	  n[0]=cellData[cellI][variableIndex(0,0)+1];
	  n[1]=-cellData[cellI][variableIndex(0,0)];
	}
	else { //Parallell to given direction
	  n[0]=cellData[cellI][variableIndex(0,0)];
	  n[1]=cellData[cellI][variableIndex(0,0)+1];		
	}
      }
      else {
	//dimension=3
	if( parameter(3) != 1.0 ) {//Perpendicular to given direction, ie plane vector parallell
	  n[0]=cellData[cellI][variableIndex(0,0)];
	  n[1]=cellData[cellI][variableIndex(0,0)+1];
	  n[2]=cellData[cellI][variableIndex(0,0)+2];
	}
	else {
	  // Parallell to given direction is undefined in 3D
	  // @Todo: Maybe should introduce the cell plane 
	  // and create the direction within the plane.
	  std::cerr << "DivisionVolumeViaDirection::update "
		    << "Parallell direction (Parallell_flag=1) and 3D "
		    << "not yet implemented!" << std::endl;
	  exit(-1);
	}
      }
    }
    else {
      //Random
      if (dimension==2) {
	double phi=2*3.14*myRandom::Rnd();
	n[0] = std::sin(phi);
	n[1] = std::cos(phi);
      }
      else {
	//dimension=3
	// @Todo: Should this be random within the cell plane
	double phi=2*3.14*myRandom::Rnd();
	double theta=2*3.14*myRandom::Rnd();
	n[0] = std::sin(phi)*std::sin(theta);
	n[1] = std::cos(phi)*std::sin(theta);
	n[2] = std::cos(theta);
      }
    }
  
    // Find two (and two only) intersecting walls
    //
    std::vector<size_t> wI(2);
    std::vector<double> v1Pos(dimension),v2Pos(dimension);
  
    if (findTwoDivisionWalls(vertexData,divCell,wI,com,n,v1Pos,v2Pos)) {
      return;
      // std::cerr << "DivisionVolumeViaDirection::update "
      // << "failed to find two walls for division!" << std::endl;
      // exit(-1);
    }
  
    // Do the division (add one cell, three walls, and two vertices)
    //
    size_t numWallTmp=wallData.size();
    assert( numWallTmp==T->numWall() );
    //Divide
    T->divideCell(divCell,wI[0],wI[1],v1Pos,v2Pos,cellData,wallData,vertexData,
		  cellDeriv,wallDeriv,vertexDeriv,variableIndex(1),
		  parameter(2));
    assert( numWallTmp+3 == T->numWall() );
  
    //Change length of new wall between the divided daugther cells 
    wallData[numWallTmp][0] *= parameter(1);
  
    //Check that the division did not mess up the data structure
    //T->checkConnectivity(1);		
  }

  VolumeRandomDirection::
  VolumeRandomDirection(std::vector<double> &paraValue, 
			std::vector< std::vector<size_t> > 
			&indValue ) {
  
    //Do some checks on the parameters and variable indeces
    //
    if ( paraValue.size() != 4) {
      std::cerr << "DivisionVolumeRandomDirection::"
		<< "DivisionVolumeRandomDirection() "
		<< "Four parameters used V_threshold, LWall_frac, "
		<< "Lwall_threshold, and COM (1 = COM, 0 = Random).\n";
      std::exit(EXIT_FAILURE);
    }
    if( indValue.size() != 1 ) {
      std::cerr << "DivisionVolumeRandomDirection::"
		<< "DivisionVolumeRandomDirection() "
		<< "Indices for volume dependent cell "
		<< "variables are at first level.\n";
      exit(0);
    }
    //Set the variable values
    //
    setId("DivisionVolumeRandomDirection");
    setNumChange(1);
    setParameter(paraValue);  
    setVariableIndex(indValue);
  
    //Set the parameter identities
    //
    std::vector<std::string> tmp( numParameter() );
    tmp.resize( numParameter() );
    tmp[0] = "V_threshold";
    tmp[1] = "LWall_frac";
    tmp[2] = "LWall_threshold";
    tmp[3] = "COM";
    setParameterId( tmp );
  }

  int VolumeRandomDirection::
  flag(Tissue *T,size_t i,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) {
  
    if( T->cell(i).calculateVolume(vertexData) > parameter(0) ) {
      std::cerr << "Cell " << i << " marked for division with volume " 
		<< T->cell(i).volume() << std::endl;
      return 1;
    } 
    return 0;
  }

  void VolumeRandomDirection::
  update(Tissue *T,size_t cellI,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDeriv,
	 DataMatrix &wallDeriv,
	 DataMatrix &vertexDeriv ) {
  
    Cell *divCell = &(T->cell(cellI));
    size_t dimension = vertexData[0].size();
    //size_t numV = divCell->numVertex();
    assert( divCell->numWall() > 2 );
    assert( dimension==2 );
  
    std::vector<double> com;
  
    if (parameter(3) == 1)
      {
	com = divCell->positionFromVertex(vertexData);
      }
    else
      {
	try
	  {
	    com = divCell->randomPositionInCell(vertexData);
	  }
	catch (Cell::FailedToFindRandomPositionInCellException)
	  {
	    return;
	  }
      }
  
    std::vector<double> n(dimension);
    double phi=2*3.14*myRandom::Rnd();
    n[0] = std::sin(phi);
    n[1] = std::cos(phi);
  
    //Find two (and two only) intersecting walls
    //
    std::vector<size_t> wI(2);
    std::vector<double> s(2);
    wI[0]=0;
    wI[1]=divCell->numWall();
    s[0]=s[1]=-1.0;
    //double minDist,w3s;
    std::vector<size_t> w3Tmp;
    std::vector<double> w3tTmp;
    int flag=0;
    for( size_t k=0 ; k<divCell->numWall() ; ++k ) {
      size_t v1Tmp = divCell->wall(k)->vertex1()->index();
      size_t v2Tmp = divCell->wall(k)->vertex2()->index();
      std::vector<double> w3(dimension),w0(dimension);
      for( size_t d=0 ; d<dimension ; ++d ) {
	w3[d] = vertexData[v2Tmp][d]-vertexData[v1Tmp][d];
	w0[d] = com[d]-vertexData[v1Tmp][d];
      }
      double a=0.0,b=0.0,c=0.0,d=0.0,e=0.0;//a=1.0
      for( size_t dim=0 ; dim<dimension ; ++dim ) {
	a += n[dim]*n[dim];
	b += n[dim]*w3[dim];
	c += w3[dim]*w3[dim];
	d += n[dim]*w0[dim];
	e += w3[dim]*w0[dim];
      }
      double fac=a*c-b*b;//a*c-b*b
      if( fac>0.0 ) {//else parallell and not applicable
	fac = 1.0/fac;
	//double s = fac*(b*e-c*d);
	double t = fac*(a*e-b*d);//fac*(a*e-b*d)
	if( t>=0.0 && t<1.0 ) {//within wall
	  //double dx0 = w0[0] +fac*((b*e-c*d)*nW2[0]+()*w3[0]); 					
	  w3Tmp.push_back(k);
	  w3tTmp.push_back(t);
	  std::cerr << "Dividing cell " << divCell->index() << " via wall "
		    << k << " at t=" << t << std::endl;
	  if( flag<2 ) {
	    s[flag] = t;
	    wI[flag] = k;
	  }				
	  flag++;
	}
      }
    }
    assert( wI[1] != divCell->numWall() && wI[0] != wI[1] );
    if( flag != 2 ) {
      return;
      // std::cerr << "divideVolumeVisStrain::update Warning"
      // 	      << " not two walls possible as connection "
      // 	      << "for cell " 
      // 	      << cellI << std::endl; 
      // for( size_t k=0 ; k<divCell->numWall() ; ++k ) {
      //   std::cerr << "0 " 
      // 		<< vertexData[divCell->wall(k)->vertex1()->index()][0]
      // 		<< " " 
      // 		<< vertexData[divCell->wall(k)->vertex1()->index()][1]
      // 		<< "\n0 " 
      // 		<< vertexData[divCell->wall(k)->vertex2()->index()][0]
      // 		<< " " 
      // 		<< vertexData[divCell->wall(k)->vertex2()->index()][1]
      // 		<< "\n\n\n";
      // }
      // for( size_t kk=0 ; kk<w3Tmp.size() ; ++kk ) {
      //   size_t k = w3Tmp[kk];
      //   std::cerr << "1 " 
      // 		<< vertexData[divCell->wall(k)->vertex1()->index()][0]
      // 		<< " " 
      // 		<< vertexData[divCell->wall(k)->vertex1()->index()][1]
      // 		<< "\n1 " 
      // 		<< vertexData[divCell->wall(k)->vertex2()->index()][0]
      // 		<< " " 
      // 		<< vertexData[divCell->wall(k)->vertex2()->index()][1]
      // 		<< "\n\n\n";
      // }
      // std::cerr << "2 " 
      // 	      << vertexData[divCell->wall(wI[0])->vertex1()->index()][0]
      // 	      << " " 
      // 	      << vertexData[divCell->wall(wI[0])->vertex1()->index()][1]
      // 	      << "\n2 " 
      // 	      << vertexData[divCell->wall(wI[0])->vertex2()->index()][0]
      // 	      << " " 
      // 	      << vertexData[divCell->wall(wI[0])->vertex2()->index()][1]
      // 	      << "\n\n\n";
      // std::cerr << "3 " 
      // 	      << vertexData[divCell->wall(wI[1])->vertex1()->index()][0]
      // 	      << " " 
      // 	      << vertexData[divCell->wall(wI[1])->vertex1()->index()][1]
      // 	      << "\n3 " 
      // 	      << vertexData[divCell->wall(wI[1])->vertex2()->index()][0]
      // 	      << " " 
      // 	      << vertexData[divCell->wall(wI[1])->vertex2()->index()][1]
      // 	      << "\n\n\n";
      // std::cerr << "4 " 
      // 	      << 0.5*(vertexData[divCell->wall(wI[0])->vertex1()->index()][0]+
      // 		      vertexData[divCell->wall(wI[0])->vertex2()->index()][0])
      // 	      << " " 
      // 	      << 0.5*(vertexData[divCell->wall(wI[0])->vertex1()->index()][1]+
      // 		      vertexData[divCell->wall(wI[0])->vertex2()->index()][1])
      // 	      << "\n4 "
      // 	      << 0.5*(vertexData[divCell->wall(wI[1])->vertex1()->index()][0]+
      // 		      vertexData[divCell->wall(wI[1])->vertex2()->index()][0])
      // 	      << " " 
      // 	      << 0.5*(vertexData[divCell->wall(wI[1])->vertex1()->index()][1]+
      // 		      vertexData[divCell->wall(wI[1])->vertex2()->index()][1])
      // 	      << "\n\n\n";
      // exit(-1);
    }	
    //Addition of new vertices at walls at position 's' 
    std::vector<double> v1Pos(dimension),v2Pos(dimension);
    size_t v1I = divCell->wall(wI[0])->vertex1()->index();
    size_t v2I = divCell->wall(wI[0])->vertex2()->index();
    for( size_t d=0 ; d<dimension ; ++d )
      v1Pos[d] = vertexData[v1I][d]+ s[0]*(vertexData[v2I][d]-vertexData[v1I][d]);
    v1I = divCell->wall(wI[1])->vertex1()->index();
    v2I = divCell->wall(wI[1])->vertex2()->index();
    for( size_t d=0 ; d<dimension ; ++d )
      v2Pos[d] = vertexData[v1I][d]+s[1]*(vertexData[v2I][d]-vertexData[v1I][d]);
  
    //Add one cell, three walls, and two vertices
    //////////////////////////////////////////////////////////////////////
    //Save number of walls
    size_t numWallTmp=wallData.size();
    assert( numWallTmp==T->numWall() );
    //Divide
    T->divideCell(divCell,wI[0],wI[1],v1Pos,v2Pos,cellData,wallData,vertexData,
		  cellDeriv,wallDeriv,vertexDeriv,variableIndex(0),
		  parameter(2));
    assert( numWallTmp+3 == T->numWall() );
  
    //Change length of new wall between the divided daugther cells 
    wallData[numWallTmp][0] *= parameter(1);
  
    //Check that the division did not mess up the data structure
    //T->checkConnectivity(1);		
  }



  VolumeRandomDirectionCenterTriangulation::
  VolumeRandomDirectionCenterTriangulation(std::vector<double> &paraValue, 
					   std::vector< std::vector<size_t> > 
					   &indValue ) {
  
    //Do some checks on the parameters and variable indeces
    //
    std::cerr << "DivisionVolumeRandomDirectionCenterTriangulation:: "
	      << " is not working adaptation is needed. "
	      <<std::endl;

    exit(0);

    if ( paraValue.size() != 1) {
      std::cerr << "DivisionVolumeRandomDirectionCenterTriangulation::"
		<< "DivisionVolumeRandomDirectionCenterTriangulation() "
		<< "One parameter used V_threshold, " << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if( indValue.size() != 2 || indValue[0].size() != 1) {
      std::cerr << "DivisionVolumeRandomDirectionCenterTriangulation::"
		<< "DivisionVolumeRandomDirectionCenterTriangulation() "
		<< "First level should store the index where the central vertex "
		<< "is located (followed by resting lengths of internal edges)."
		<< "Indices for volume dependent cell "
		<< "variables are at second level." << std::endl;
      exit(0);
    }
    //Set the variable values
    //
    setId("DivisionVolumeRandomDirectionCenterTriangulation");
    setNumChange(1);
    setParameter(paraValue);  
    setVariableIndex(indValue);
  
    //Set the parameter identities
    //
    std::vector<std::string> tmp( numParameter() );
    tmp.resize( numParameter() );
    tmp[0] = "V_threshold";
    setParameterId( tmp );
  }

  int VolumeRandomDirectionCenterTriangulation::
  flag(Tissue *T,size_t i,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) {
  
    if( T->cell(i).calculateVolumeCenterTriangulation(vertexData,cellData,variableIndex(0,0)) > 
	parameter(0) ) {
      std::cerr << "Cell " << i << " marked for division with volume " 
		<< T->cell(i).volume() << std::endl;
      return 1;
    } 
    return 0;
  }

  void VolumeRandomDirectionCenterTriangulation::
  update(Tissue *T,size_t cellI,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDeriv,
	 DataMatrix &wallDeriv,
	 DataMatrix &vertexDeriv ) {
  
    Cell *divCell = &(T->cell(cellI));
    size_t cellIndex = divCell->index();
    size_t dimension = vertexData[0].size();
    size_t numV = divCell->numVertex();
    assert( divCell->numWall() > 2 );
    assert( dimension==3 );
  
    // Find first vertex (random)
    size_t b;
    size_t counter = 0;
    do {
      b = (size_t) myRandom::Rnd()*numV; 
      counter++;
    } while (divCell->vertex(b)->numWall()>2 && counter<1000);
    if (counter>999) {
      std::cerr << "DivisionVolumeRandomDirectionCenterTriangulation::"
		<< "DivisionVolumeRandomDirectionCenterTriangulation() " << std::endl
		<< " Could not find vertex with two walls." << std::endl;
      exit(EXIT_FAILURE);
    }

    // Select second vertex creating as equally sized daughters as possible.
    double halfArea = 0.5*divCell->volume();
    double tmpArea = 0.0;

    std::vector<double> center(dimension);
    for (size_t i=0; i<dimension; i++)
      center[i] = cellData[cellIndex][variableIndex(0,0)+i];
    size_t best=0;
    double minArea = divCell->volume();
    for (size_t i=1; i<numV; i++)
      { 
	std::vector<double> v1pos=vertexData[ divCell->vertex((b+i-1)%numV)->index() ];
	std::vector<double> v2pos=vertexData[divCell->vertex((b+i)%numV)->index()];
	tmpArea += myMath::areaTriangle(center,v1pos,v2pos);
	if ( std::abs(tmpArea-halfArea)< minArea && 
	     divCell->vertex((b+i)%numV)->numWall()==2 ) {  
	  best=i;
	  minArea = std::abs( tmpArea-halfArea );
	}
      }
    size_t e = (b+best)%(numV); // e is the index of the end vertex of division wall 
  
    // Divide the cell
    //size_t numWallTmp = T->numWall();
    //size_t numVertexTmp = T->numVertex();
    //size_t numCellTmp = T->numCell();
 
    // T->divideCellCenterTriangulation(divCell,b,e,variableIndex(0,0),
    //       			   cellData,wallData,vertexData,
    //       			   cellDeriv,wallDeriv,vertexDeriv,
    //       			   variableIndex(1));
  
    //assert( numWallTmp+6 == T->numWall() );
    //assert( numVertexTmp+5 == T->numVertex() );
    //assert( numCellTmp+1 == T->numCell() );
  
    //Check that the division did not mess up the data structure
    //T->checkConnectivity(1);		
  }


  ForceDirection::ForceDirection(std::vector<double> &paraValue,
				 std::vector< std::vector<size_t> > &indValue)
  {
    if (paraValue.size() != 4) {
      std::cerr << "DivisionForceDirection::ForceDirection() "
		<< "Four parameters are used V_threshold, Lwall_fraction, Lwall_threshold and orientation_flag." << std::endl;
      exit(EXIT_FAILURE);
    }
	
    if (indValue.size() != 2 || indValue[0].size() != 1) {
      std::cerr << "DivisionForceDirection::ForceDirection() "
		<< "First level: Variable indices for volume dependent cell "
		<< "variables are used.\n"
		<< "Second level: Variable indices for forces.\n";
      exit(EXIT_FAILURE);
    }
	
    setId("DivisionForceDirection");
    setNumChange(1);
    setParameter(paraValue);  
    setVariableIndex(indValue);
  
    std::vector<std::string> tmp(numParameter());
    tmp.resize (numParameter());
    tmp[0] = "V_threshold";
    tmp[1] = "Lwall_threshold";
    tmp[2] = "orientation_flag";
    setParameterId(tmp);
  }

  int ForceDirection::flag(Tissue *T, size_t i,
			   DataMatrix &cellData,
			   DataMatrix &wallData,
			   DataMatrix &vertexData,
			   DataMatrix &cellDerivs,
			   DataMatrix &wallDerivs,
			   DataMatrix &vertexDerivs)
  {
    if (T->cell(i).calculateVolume(vertexData) > parameter(0)) {
      std::cerr << "Cell " << i << " marked for division with volume " 
		<< T->cell(i).volume() << std::endl;
      return 1;
    } 
    return 0;
  }

  void ForceDirection::update(Tissue* T, size_t i,
			      DataMatrix &cellData,
			      DataMatrix &wallData,
			      DataMatrix &vertexData,
			      DataMatrix &cellDerivs,
			      DataMatrix &wallDerivs,
			      DataMatrix &vertexDerivs)
  {
    Cell cell = T->cell(i);
    assert(cell.numWall() > 1);
    assert(vertexData[0].size() == 2); // Make sure dimension == 2
  
  
    // Calculate force direction (nx, ny)
    double nx = 0.0;
    double ny = 0.0;
  
    for (size_t i = 0; i < cell.numWall(); ++i) {
      Wall *wall = cell.wall(i);
      double wx = wall->vertex1()->position(0) - wall->vertex2()->position(0);
      double wy = wall->vertex1()->position(1) - wall->vertex2()->position(1);
      double Aw = std::sqrt(wx * wx  + wy * wy);
      if (wx > 0) {
	wx = wx / Aw;
	wy = wy / Aw;
      } else {
	wx = -1 * wx / Aw;
	wy = -1 * wy / Aw;
      }
    
      double force = 0.0;
      for (size_t j = 0; j < numVariableIndex(1); ++j)
	force += wallData[wall->index()][variableIndex(1, j)];
    
      nx += wx * force;
      ny += wy * force;
    }
  
    double A = std::sqrt(nx * nx + ny * ny);
  
    double nxp;
    double nyp;
    if (parameter(3) == 0) {
      nxp = nx / A;
      nyp = ny / A;
    } else {
      nxp = -ny / A;
      nyp = nx / A;
    }
    nx = nxp;
    ny = nyp;
  
    // Calculate mean vertex position (mx, my)
    // by center of mass).
    double mx = 0.0;
    double my = 0.0;
  
    for (size_t i = 0; i < cell.numVertex(); ++i) {
      Vertex *vertex = cell.vertex(i);
      mx += vertex->position(0);
      my += vertex->position(1);
    }
    mx /= cell.numVertex();
    my /= cell.numVertex();
  
    // Find candidate walls for division (Patrik: See RR014)
    std::vector<size_t> candidateWalls;
    DataMatrix verticesPosition;
  
    for (size_t i = 0; i < cell.numWall(); ++i) {
      Wall *wall = cell.wall(i);
      Vertex *v1 = wall->vertex1();
      Vertex *v2 = wall->vertex2();
    
      double x1x = v1->position(0);
      double x1y = v1->position(1);
      double x2x = v2->position(0);
      double x2y = v2->position(1);
    
      double detA = nx * (x1y - x2y) - ny * (x1x - x2x);
    
      if (detA == 0)
	continue;
    
      // double s = ((x1y - x2y) * (x1x - mx) - (x1x - x2x) * (x1y - my)) / detA;
      double t = (-ny * (x1x - mx) + nx * (x1y - my)) / detA;
    
      if (t <= 0.0 || t >= 1.0)
	continue;
      else {
	candidateWalls.push_back(i);
	std::vector<double> position(2);
	position[0] = x1x + t * (x2x - x1x);
	position[1] = x1y + t * (x2y - x1y);
	verticesPosition.push_back(position);
      }
    }
  
    if (candidateWalls.size() != 2) {
      std::cerr << "DivisionForceDirection::update() "
		<< "More than two or less than one candidate walls for division." 
		<< std::endl;
      exit(EXIT_FAILURE);
    }
  
    T->divideCell(&cell, candidateWalls[0], candidateWalls[1],
		  verticesPosition[0], verticesPosition[1],
		  cellData, wallData, vertexData, cellDerivs, wallDerivs, vertexDerivs,
		  variableIndex(0), parameter(2));
  
    //Change length of new wall between the divided daugther cells 
    wallData[T->numWall()-1][0] *= parameter(1);
  }

  //!Constructor
  VolumeViaShortestPath::
  VolumeViaShortestPath(std::vector<double> &paraValue, 
			std::vector< std::vector<size_t> > 
			&indValue ) 
  {  
    // Do some checks on the parameters and variable indeces
    if( paraValue.size()!=3 ) {
      std::cerr << "DivisionVolumeViaShortestPath::"
		<< "DivisionVolumeViaShortestPath() "
		<< "Three parameters used V_threshold, LWall_frac, and"
		<< "Lwall_threshold." << std::endl;
      exit(0);
    }
    if( indValue.size() != 1 ) {
      std::cerr << "DivisionVolumeViaShortestPath::"
		<< "DivisionVolumeViaShortestPath() "
		<< "Variable indices for volume dependent cell "
		<< "variables is used.\n";
      exit(0);
    }
    // Set the variable values
    setId("DivisionVolumeViaShortestPath");
    setNumChange(1);
    setParameter(paraValue);  
    setVariableIndex(indValue);
  
    // Set the parameter identities
    std::vector<std::string> tmp( numParameter() );
    tmp.resize( numParameter() );
    tmp[0] = "V_threshold";
    tmp[1] = "LWall_frac";
    tmp[2] = "LWall_threshold";
    setParameterId( tmp );
  }

  int VolumeViaShortestPath::flag(Tissue *T,size_t i,
				  DataMatrix &cellData,
				  DataMatrix &wallData,
				  DataMatrix &vertexData,
				  DataMatrix &cellDerivs,
				  DataMatrix &wallDerivs,
				  DataMatrix &vertexDerivs ) 
  {	
    if( T->cell(i).calculateVolume(vertexData) > parameter(0) ) {
      std::cerr << "Cell " << i << " marked for division with volume " 
		<< T->cell(i).volume() << std::endl;
      return 1;
    } 
    return 0;
  }
  
  void VolumeViaShortestPath::
  update(Tissue *T,size_t i,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDeriv,
	 DataMatrix &wallDeriv,
	 DataMatrix &vertexDeriv ) 
  {  
    Cell *divCell = &(T->cell(i));
    size_t dimension = vertexData[0].size();
    assert( divCell->numWall() > 1 );
    assert( dimension==2 );
    
    // Find shortest path for each wall pair
    size_t wI=0;
    //  size_t wI2=1;
    for (size_t k=0; k<divCell->numWall()-1; ++k) {
      for (size_t k2=k+1; k<divCell->numWall(); ++k) {
	std::vector<double> x0(dimension),x1(dimension),x2(dimension),n1(dimension),n2(dimension);
	x1[0] = divCell->wall(k)->vertex1()->position(0);
	x1[1] = divCell->wall(k)->vertex1()->position(1);
	x2[0] = divCell->wall(k2)->vertex1()->position(0);
	x2[1] = divCell->wall(k2)->vertex1()->position(1);
	
	
	n1[0] = divCell->wall(k)->vertex2()->position(0)-x1[0];
	n1[1] = divCell->wall(k)->vertex2()->position(1)-x1[1];
	n2[0] = divCell->wall(k2)->vertex2()->position(0)-x2[0];
	n2[1] = divCell->wall(k2)->vertex2()->position(1)-x2[1];
	double wL1 = std::sqrt(n1[0]*n1[0]+n1[1]*n1[1]);
	double wL2 = std::sqrt(n2[0]*n2[0]+n2[1]*n2[1]);
	
	//std::cerr << "x1[0]: " << x1[0] << "\n";
	//std::cerr << "x1[1]: " << x1[1] << "\n";
	//std::cerr << "pos0: " << divCell->wall(k2)->vertex2()->position(0) << "\n";
	//std::cerr << "pos1: " << divCell->wall(k2)->vertex2()->position(1) << "\n";
	//std::cerr << "wL1: " << wL1 << "\n";
	
	double t1=0.0, t2=0.0;
	double denominator = n2[0]*n2[1]-n2[0]*n1[1];
	if (denominator != 0.0) {
	  t1 = (n2[0]*(x2[1]-x1[1])-n2[1]*(x2[0]-x1[0]))/denominator;
	  Vertex *v1=divCell->wall(k)->vertex1();
	  Vertex *v1Other=divCell->wall(k)->vertex2();
	  if (t1>0.0) {
	    t1 = wL1-t1;
	    x1[0] = divCell->wall(k)->vertex2()->position(0);
	    x1[1] = divCell->wall(k)->vertex2()->position(1);
	    n1[0] = -n1[0];
	    n1[1] = -n1[1];
	    v1=divCell->wall(k)->vertex2();
	    v1Other=divCell->wall(k)->vertex1();
	  }
	  t2 = (n1[0]*(x2[1]-x1[1])-n1[1]*(x2[0]-x1[0]))/denominator;
	  Vertex *v2=divCell->wall(k2)->vertex1();
	  Vertex *v2Other=divCell->wall(k2)->vertex2();
	  if (t2>0.0) {
	    t2 = wL2-t2;
	    x2[0] = divCell->wall(k2)->vertex2()->position(0);
	    x2[1] = divCell->wall(k2)->vertex2()->position(1);
	    n2[0] = -n2[0];
	    n2[1] = -n2[1];
	    v2=divCell->wall(k2)->vertex2();
	    v2Other=divCell->wall(k2)->vertex1();
	  }
	  if( v1==v2 ) 
	    continue;
	  double ACell = divCell->calculateVolume(vertexData);
	  double Ahalf = 0.5*ACell;
	  //
	  // Calculate A2
	  //
	  size_t vk=0;
	  size_t oppositeVolumeFlag=0;
	  // Find start index for volume calculation
	  while (divCell->vertex(vk) != v1 && divCell->vertex(vk) != v2)
	    ++vk;
	  std::vector<Vertex*> vertices;
	  vertices.push_back(divCell->vertex(vk++));
	  while (divCell->vertex(vk) != v1 && divCell->vertex(vk) != v2) {
	    if ( divCell->vertex(vk) == v1Other || 
		 divCell->vertex(vk) == v2Other )
	      oppositeVolumeFlag++;
	    vertices.push_back(divCell->vertex(vk++));
	  }
	  vertices.push_back(divCell->vertex(vk));
	  assert( vertices[0]==v1 || vertices[0]==v2 );
	  assert( vertices[vertices.size()-1]==v1 || 
		  vertices[vertices.size()-1]==v2 );
	  // Calculate area from extracted vertices
	  double A2=0.0;
	  for (vk=0; vk<vertices.size(); ++vk) {
	    size_t vkPlus=(vk+1)%vertices.size();
	    A2 += vertices[vk]->position(0)*vertices[vkPlus]->position(1) -
	      vertices[vkPlus]->position(0)*vertices[vk]->position(1);
	  }
	  A2 = 0.5*std::fabs(A2);
	  if (oppositeVolumeFlag)
	    A2 = ACell-A2;
	  if (A2>Ahalf)
	    std::cerr << "Cell " << divCell->index() << " walls " << k << "," 
		      << k2 << " not applicable" << std::endl;
	  x0[0] = x1[0] + n1[0]*t1;
	  x0[1] = x1[1] + n1[1]*t1;
	  double n1n2 = n1[0]*n2[0]+n1[1]*n2[1];
	  double alpha = std::acos(n1n2/(wL1*wL2));
	  double A0 = 0.5*((x1[0]-x0[0])*(x2[1]-x0[1])-
			   (x1[1]-x0[1])*(x2[0]-x0[0]));

	  double root = std::sqrt(std::fabs((A0+Ahalf-A2)/(std::cos(alpha)*std::sin(alpha)))); // not needed?
	  //added fabs of argument \andre
	  double t1A = std::sqrt((x1[0]-x0[0])*(x1[0]-x0[0])+(x1[1]-x0[1])*(x1[1]-x0[1]));
	  double t2A = std::sqrt((x2[0]-x0[0])*(x2[0]-x0[0])+(x2[1]-x0[1])*(x2[1]-x0[1]));
	  // double length = 2*(2*t1A+root)*std::cos(alpha);

	  //std::cerr << "A0: " << A0 << " Ahalf: " << Ahalf << " A2: " << A2 << "\n";
	  //std::cerr << "root**2: " << (A0+Ahalf-A2)/(std::cos(alpha)*std::sin(alpha)) << "\n";
	  //std::cerr << "root: " << root << "\n";
	  //std::cerr << "t1A: " << t1A << "\n";
	  //std::cerr << "root: " << root << "\n";

	  std::cerr << divCell->index() << " " << k << " " << k2 
		    << "\t" << wL1 << " " << t1A+root << " (" << t1A-root << ")  "
		    << wL2 << " " << t2A+root << " (" << t2A-root << ")" << std::endl;
	}
	else {
	}	      
	//double tmpLength = divCell->wall(k)->setLengthFromVertexPosition(vertexData);
	//if( tmpLength > maxLength ) {
	//wI=k;
	//maxLength = tmpLength;
      }
    }
    //exit(0);
    double maxLength=0.0;
    
    std::vector<double> nW(dimension),nW2(dimension),v1Pos(dimension),
      v2Pos(dimension);
    size_t v1wI = divCell->wall(wI)->vertex1()->index();
    size_t v2wI = divCell->wall(wI)->vertex2()->index();
  
    for( size_t d=0 ; d<dimension ; ++d ) {
      nW[d] = (vertexData[v1wI][d]-vertexData[v2wI][d])/maxLength;
      v1Pos[d] = 0.5*(vertexData[v1wI][d]+vertexData[v2wI][d]);
    }
    nW2[1] = nW[0];
    nW2[0] = -nW[1];
  
    // Find intersection with another wall
    //
    size_t w3I=divCell->numWall();
    //double minDist,w3s;
    std::vector<size_t> w3Tmp;
    std::vector<double> w3tTmp;
    int flag=0,vertexFlag=0;
    for( size_t k=0 ; k<divCell->numWall() ; ++k ) {
      if( k!=wI ) {
	size_t v1w3Itmp = divCell->wall(k)->vertex1()->index();
	size_t v2w3Itmp = divCell->wall(k)->vertex2()->index();
	std::vector<double> w3(dimension),w0(dimension);
	for( size_t d=0 ; d<dimension ; ++d ) {
	  w3[d] = vertexData[v2w3Itmp][d]-vertexData[v1w3Itmp][d];
	  w0[d] = v1Pos[d]-vertexData[v1w3Itmp][d];
	}
	double a=0.0,b=0.0,c=0.0,d=0.0,e=0.0;//a=1.0
	for( size_t dim=0 ; dim<dimension ; ++dim ) {
	  a += nW2[dim]*nW2[dim];
	  b += nW2[dim]*w3[dim];
	  c += w3[dim]*w3[dim];
	  d += nW2[dim]*w0[dim];
	  e += w3[dim]*w0[dim];
	}
	double fac=a*c-b*b;//a*c-b*b
	if( fac>1e-10 ) {//else parallell and not applicable
	  fac = 1.0/fac;
	  //double s = fac*(b*e-c*d);
	  double t = fac*(a*e-b*d);//fac*(a*e-b*d)
	  if( t>0.0 && t<=1.0 ) {//within wall
	    //double dx0 = w0[0] +fac*((b*e-c*d)*nW2[0]+()*w3[0]); 					
	    flag++;
	    if( t==1.0 )
	      vertexFlag++;
	    w3I = k;
	    w3Tmp.push_back(k);
	    w3tTmp.push_back(t);
	  }
	}
      }
    }
    assert( w3I != divCell->numWall() && w3I != wI );
    if( flag != 1 && !(flag==2 && vertexFlag)) {
      std::cerr << "divideVolumeViaShortestPath::update Warning"
		<< " more than one wall possible as connection "
		<< "for cell " 
		<< i << std::endl; 
      std::cerr << flag << " " << vertexFlag << std::endl; 
      for( size_t k=0 ; k<divCell->numWall() ; ++k ) {
	std::cerr << "0 " 
		  << vertexData[divCell->wall(k)->vertex1()->index()][0]
		  << " " 
		  << vertexData[divCell->wall(k)->vertex1()->index()][1]
		  << "\n0 " 
		  << vertexData[divCell->wall(k)->vertex2()->index()][0]
		  << " " 
		  << vertexData[divCell->wall(k)->vertex2()->index()][1]
		  << "\n\n\n";
      }
      for( size_t kk=0 ; kk<w3Tmp.size() ; ++kk ) {
	size_t k = w3Tmp[kk];
	std::cerr << "1 " 
		  << vertexData[divCell->wall(k)->vertex1()->index()][0]
		  << " " 
		  << vertexData[divCell->wall(k)->vertex1()->index()][1]
		  << "\n1 " 
		  << vertexData[divCell->wall(k)->vertex2()->index()][0]
		  << " " 
		  << vertexData[divCell->wall(k)->vertex2()->index()][1]
		  << "\n\n\n";
      }
      std::cerr << "2 " 
		<< vertexData[divCell->wall(wI)->vertex1()->index()][0]
		<< " " 
		<< vertexData[divCell->wall(wI)->vertex1()->index()][1]
		<< "\n2 " 
		<< vertexData[divCell->wall(wI)->vertex2()->index()][0]
		<< " " 
		<< vertexData[divCell->wall(wI)->vertex2()->index()][1]
		<< "\n\n\n";
      std::cerr << "3 " 
		<< vertexData[divCell->wall(w3I)->vertex1()->index()][0]
		<< " " 
		<< vertexData[divCell->wall(w3I)->vertex1()->index()][1]
		<< "\n3 " 
		<< vertexData[divCell->wall(w3I)->vertex2()->index()][0]
		<< " " 
		<< vertexData[divCell->wall(w3I)->vertex2()->index()][1]
		<< "\n\n\n";
      std::cerr << "4 " 
		<< 0.5*(vertexData[divCell->wall(wI)->vertex1()->index()][0]+
			vertexData[divCell->wall(wI)->vertex2()->index()][0])
		<< " " 
		<< 0.5*(vertexData[divCell->wall(wI)->vertex1()->index()][1]+
			vertexData[divCell->wall(wI)->vertex2()->index()][1])
		<< "\n4 "
		<< 0.5*(vertexData[divCell->wall(w3I)->vertex1()->index()][0]+
			vertexData[divCell->wall(w3I)->vertex2()->index()][0])
		<< " " 
		<< 0.5*(vertexData[divCell->wall(w3I)->vertex1()->index()][1]+
			vertexData[divCell->wall(w3I)->vertex2()->index()][1])
		<< "\n\n\n";
      exit(-1);
    }	
    size_t v1w3I = divCell->wall(w3I)->vertex1()->index();
    size_t v2w3I = divCell->wall(w3I)->vertex2()->index();
    //Set the vertex using the collected t
    for( size_t d=0 ; d<dimension ; ++d )
      v2Pos[d] = vertexData[v1w3I][d] + 
	w3tTmp[w3tTmp.size()-1]*(vertexData[v2w3I][d]-vertexData[v1w3I][d]);
    //for( size_t d=0 ; d<dimension ; ++d )
    //v2Pos[d] = 0.5*(vertexData[v2w3I][d]+vertexData[v1w3I][d]);
  
  
    // Add one cell, three walls, and two vertices
    //
    //Save number of walls
    size_t numWallTmp=wallData.size();
    assert( numWallTmp==T->numWall() );
    //Divide
    T->divideCell(divCell,wI,w3I,v1Pos,v2Pos,cellData,wallData,vertexData,
		  cellDeriv,wallDeriv,vertexDeriv,variableIndex(0),
		  parameter(2));
    assert( numWallTmp+3 == T->numWall() );
  
    //Change length of new wall between the divided daugther cells 
    wallData[numWallTmp][0] *= parameter(1);
  
    //Check that the division did not mess up the data structure
    //T->checkConnectivity(1);	
  }

  ShortestPath::ShortestPath(std::vector<double> &paraValue, 
			     std::vector< std::vector<size_t> > &indValue)
  {
    if ( paraValue.size() != 4 &&  paraValue.size() != 6 ) {
      std::cerr << "DivisionShortestPath::DivisionShortestPath() "
		<< "Four or six parameters are used V_threshold, Lwall_fraction, "
		<< "Lwall_threshold, and COM (1 = COM, 0 = Random) "
		<< "If six parameters are used, two additional parameters are for " 
		<< "centerTriangulation(1) and double resting length (1: double, 0:single). "
		<<std::endl;
      std::exit(EXIT_FAILURE);
    }
  
    if ((indValue.size() == 2 && indValue[1].size() != 1) || 
	(indValue.size() == 3 && indValue[2].size() != 2)) {
      std::cerr << "DivisionShortestPath::DivisionShortestPath() "
		<< "First level: Variable indices for volume dependent cell variables are used.\n"
		<< "Second level (optional): Cell time index."
		<< "Third level (if centerTriangulated): first:Com index, second: wall length index"
		<<std::endl;
      exit(EXIT_FAILURE);
    }
  
    setId("DivisionShortestPath");
    setNumChange(1);
    setParameter(paraValue);  
    setVariableIndex(indValue);
  
    std::vector<std::string> tmp(numParameter());
    tmp.resize (numParameter());
    tmp[0] = "V_threshold";
    tmp[1] = "Lwall_fraction";
    tmp[2] = "Lwall_threshold";
    tmp[3] = "COM";
    if(numParameter()==6){
      tmp[4] = "centerTriangulationFlag";
      tmp[5] = "doubleLengthFlag";
    }
    setParameterId(tmp);
  }

  int ShortestPath::flag(Tissue *T, size_t i,
			 DataMatrix &cellData,
			 DataMatrix &wallData,
			 DataMatrix &vertexData,
			 DataMatrix &cellDerivs,
			 DataMatrix &wallDerivs,
			 DataMatrix &vertexDerivs)
  {
    if (T->cell(i).calculateVolume(vertexData) > parameter(0))
      {
	return 1;
      } 
    else
      {
	return 0;
      }
  }

  void ShortestPath::update(Tissue* T, size_t i,
			    DataMatrix &cellData,
			    DataMatrix &wallData,
			    DataMatrix &vertexData,
			    DataMatrix &cellDerivs,
			    DataMatrix &wallDerivs,
			    DataMatrix &vertexDerivs)
  {
    Cell &cell = T->cell(i);
  
    // if (vertexData[0].size() != 2)
    //   {
    //     std::cerr << "Division::ShortestPath only supports two dimensions.\n";
    //     std::exit(EXIT_FAILURE);
    //   }
  

    // projecting the cell to the perpendicular plane to the averaged normal
    double rot[3][3]={{1,0,0},{0,1,0},{0,0,1}};
    std::vector<std::vector<double> > verticesPosition;
    verticesPosition.resize(cell.numVertex());
    std::vector<double> centerTmp(3), COMTmp(3); 
  
    if (vertexData[0].size() == 3){
      // storing the vertex positions
      for(size_t k=0; k< cell.numVertex(); ++k) {
	verticesPosition[k].resize(3);
	size_t Vind=cell.vertex(k) -> index();
	verticesPosition[k][0]=vertexData[Vind][0];
	verticesPosition[k][1]=vertexData[Vind][1];
	verticesPosition[k][2]=vertexData[Vind][2];
      }
      // storing COM position    
      if(numParameter()==6 && parameter(4)==1) {// centerTriangulation
	COMTmp[0]=cellData[i][variableIndex(2,0)  ];
	COMTmp[1]=cellData[i][variableIndex(2,0)+1];
	COMTmp[2]=cellData[i][variableIndex(2,0)+2];
      }    
      // Finding the average normal vector to the cell plane
      std::vector<double> normal(3);
      normal[0]=0;
      normal[1]=0;
      normal[2]=0;
    
      //  double area=0;
      for (size_t k=1; k< cell.numWall()-1; ++k) {
	std::vector<double> x0(3),x1(3),x2(3);
	size_t ind1=cell.vertex(0) -> index();
	size_t ind2=cell.vertex(k) -> index();
	size_t ind3=cell.vertex(k+1) -> index();
	//std::cerr<<ind1<<" "<<ind2<<" "<<ind3<<" "<<std::endl;
	//normal to the element
	std::vector<double> temp(3);
	temp[0]=(vertexData[ind2][1]-vertexData[ind1][1])*(vertexData[ind3][2]-vertexData[ind1][2])
	  -(vertexData[ind2][2]-vertexData[ind1][2])*(vertexData[ind3][1]-vertexData[ind1][1]);
	temp[1]=(vertexData[ind2][2]-vertexData[ind1][2])*(vertexData[ind3][0]-vertexData[ind1][0])
	  -(vertexData[ind2][0]-vertexData[ind1][0])*(vertexData[ind3][2]-vertexData[ind1][2]);
	temp[2]=(vertexData[ind2][0]-vertexData[ind1][0])*(vertexData[ind3][1]-vertexData[ind1][1])
	  -(vertexData[ind2][1]-vertexData[ind1][1])*(vertexData[ind3][0]-vertexData[ind1][0]);
    
	normal[0]+=temp[0];
	normal[1]+=temp[1];
	normal[2]+=temp[2];
	//std::cerr<<normal[0]<<" "<<normal[1]<<" "<<normal[2]<<" "<<area<<std::endl;  
      
	//total area of the cell
	//area +=std::sqrt(temp[0]*temp[0]+temp[1]*temp[1]+temp[2]*temp[2])/2;
      }
      double norm=std::sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
      // normalizing the normal vector
      normal[0]/=norm;
      normal[1]/=norm;
      normal[2]/=norm;
      //std::cerr<<normal[0]<<" "<<normal[1]<<" "<<normal[2]<<" "<<area<<std::endl;  
    
    
      //std::cerr<<deltaTet<<" "<<nplane<<std::endl;  
      // rotation matrix for going to the perpendicular plane to thr averaged normal
  
      if(normal[2]!=1){ // if there rotation matrix is not 1
	rot[0][0]=normal[2]+((normal[1]*normal[1])/(normal[2]+1));
	rot[1][1]=normal[2]+((normal[0]*normal[0])/(normal[2]+1));
	rot[0][1]= - (normal[0]*normal[1])/(normal[2]+1) ;
	rot[1][0]= - (normal[1]*normal[0])/(normal[2]+1) ;
      
	rot[0][2]=-normal[0];
	rot[2][0]= normal[0];
	rot[1][2]=-normal[1];
	rot[2][1]= normal[1];
	rot[2][2]= normal[2];
      
      }
      //rotating the cell vertices
      //std::vector<double> tmpVertex(3);
    
      for(size_t k=0; k< cell.numVertex(); ++k){
      
	size_t Vind=cell.vertex(k) -> index();
	for(size_t ii=0; ii< 3; ++ii) {
	  vertexData[Vind][ii]=0;
	  for(size_t jj=0; jj< 3; ++jj){ 
	    vertexData[Vind][ii]+=rot[ii][jj]*verticesPosition[k][jj];
	  }
	}
      }

      // rotating the center is cell is centertriangulated
      if(numParameter()==6 && parameter(4)==1) {// centerTriangulation
	cellData[i][variableIndex(2,0)  ]=rot[0][0]*COMTmp[0]+rot[0][1]*COMTmp[1]+rot[0][2]*COMTmp[2];
	cellData[i][variableIndex(2,0)+1]=rot[1][0]*COMTmp[0]+rot[1][1]*COMTmp[1]+rot[1][2]*COMTmp[2];
	cellData[i][variableIndex(2,0)+2]=rot[2][0]*COMTmp[0]+rot[2][1]*COMTmp[1]+rot[2][2]*COMTmp[2];
      }

      centerTmp=cell.positionFromVertex(vertexData);


      // std::cerr<<"normal "<<normal[0]<<" "<<normal[1]<<" "<<normal[2]<<std::endl;
      //  for(size_t k=0; k< cell.numVertex(); ++k){
      //    size_t Vind=cell.vertex(k) -> index();
      //    std::cerr<<"vertex "<<k<<" "
      //             <<vertexData[Vind][0]<<" "
      //             <<vertexData[Vind][1]<<" "
      //             <<vertexData[Vind][2]<<std::endl;
      //  }
      //  std::cerr<<"cell center "
      //           <<cellData[i][variableIndex(2,0)  ]<<" "
      //           <<cellData[i][variableIndex(2,0)+1]<<" "
      //           <<cellData[i][variableIndex(2,0)+2]<<std::endl;
      //  std::cerr<<"cell centroid "
      //           <<centerTmp[0]<<" "
      //           <<centerTmp[1]<<" "
      //           <<centerTmp[2]<<std::endl;

    
  
  
    }
  
  

  
    //  // Find intersection with another wall
    //  //
    //  size_t w3I=divCell->numWall();
    //  //double minDist,w3s;
    // std::vector<size_t> w3Tmp;
    //  std::vector<double> w3tTmp;
    //  int flag=0,vertexFlag=0;
    //  for( size_t k=0 ; k<divCell->numWall() ; ++k ) {
    //    if( k!=wI ) {
    //      size_t v1w3Itmp = divCell->wall(k)->vertex1()->index();
    //      size_t v2w3Itmp = divCell->wall(k)->vertex2()->index();
    //      std::vector<double> w3(dimension),w0(dimension);
    //      for( size_t d=0 ; d<dimension ; ++d ) {
    //        w3[d] = vertexData[v2w3Itmp][d]-vertexData[v1w3Itmp][d];
    //        w0[d] = v1Pos[d]-vertexData[v1w3Itmp][d];

    std::vector<Candidate> candidates = 
      getCandidates(T, i, cellData, wallData, vertexData, cellDerivs, wallDerivs, vertexDerivs);
  
    if (candidates.size() == 0)
      {
	return;
      }
    
    Candidate winner = { std::numeric_limits<double>::max(), 0, 0, 0, 0, 0, 0 };
    for (size_t i = 0; i < candidates.size(); ++i) {
      if (candidates[i].distance < winner.distance) {
	winner = candidates[i];
      }
    }
  
    // std::cerr << "Winner: " << std::endl
    //           << " distance = " << winner.distance << std::endl
    //           << " p = (" << winner.px << ", " << winner.py << ")" << std::endl
    //           << " q = (" << winner.qx << ", " << winner.qy << ")" << std::endl;
  
    size_t numWallTmp = wallData.size();
    assert(numWallTmp == T->numWall());
  
  
    std::vector<double> p(3);
    p[0] = winner.px;
    p[1] = winner.py;
    std::vector<double> q(3);
    q[0] = winner.qx;
    q[1] = winner.qy;
  
    if (vertexData[0].size() == 3){

      double px=p[0],py=p[1],qx=q[0], qy=q[1];
      // p.resize(3);
      // q.resize(3);




      // size_t V1W1=(cell.wall(winner.wall1) -> vertex1()) -> index();
      // size_t V2W1=(cell.wall(winner.wall1) -> vertex2()) -> index();
      // size_t V1W2=(cell.wall(winner.wall2) -> vertex1()) -> index();
      // size_t V2W2=(cell.wall(winner.wall2) -> vertex2()) -> index();
    
      // std::cerr<<std::endl;
      // std::cerr<<" v1w1 "<<vertexData[V1W1][0]<<"  "<<vertexData[V1W1][1]<<"  "<<vertexData[V1W1][2]<<std::endl;
      // std::cerr<<" v2w1 "<<vertexData[V2W1][0]<<"  "<<vertexData[V2W1][1]<<"  "<<vertexData[V2W1][2]<<std::endl;
      // std::cerr<<" v1w2 "<<vertexData[V1W2][0]<<"  "<<vertexData[V1W2][1]<<"  "<<vertexData[V1W2][2]<<std::endl;
      // std::cerr<<" v2w2 "<<vertexData[V2W2][0]<<"  "<<vertexData[V2W2][1]<<"  "<<vertexData[V2W2][2]<<std::endl;
      // std::cerr<<" p "<<p[0]<<"  "<<p[1]<<"  "<<p[2]<<std::endl;
      // std::cerr<<" q "<<q[0]<<"  "<<q[1]<<"  "<<q[2]<<std::endl;
      // std::cerr<<" com "<<COMTmp[0]<<"  "<<COMTmp[1]<<"  "<<COMTmp[2]<<std::endl;

      // //testing the rotation matrix    
      // std::vector<std::vector<double> > rotTemp(3);
      // for(size_t iii=0; iii< 3; ++iii)
      //   rotTemp[iii].resize(3);
      // for(size_t iii=0; iii< 3; ++iii)
      //   for(size_t jjj=0; jjj<3; ++jjj)
      //     for(size_t kkk=0; kkk<3; ++kkk)
      //       rotTemp[iii][jjj]+=rot[iii][kkk]*rot[jjj][kkk];
      // for(size_t iii=0; iii< 3; ++iii)
      //   std::cerr<<"rot   "<<rotTemp[iii][0]<<" "<<rotTemp[iii][1]<<" "<<rotTemp[iii][2]<<std::endl;
    
      //rotating back the found positions for division points
    
      p[0]= rot[0][0]*px + rot[1][0]*py + rot[2][0]*centerTmp[2];
      p[1]= rot[0][1]*px + rot[1][1]*py + rot[2][1]*centerTmp[2];
      p[2]= rot[0][2]*px + rot[1][2]*py + rot[2][2]*centerTmp[2];
    
      q[0]= rot[0][0]*qx + rot[1][0]*qy + rot[2][0]*centerTmp[2];
      q[1]= rot[0][1]*qx + rot[1][1]*qy + rot[2][1]*centerTmp[2];
      q[2]= rot[0][2]*qx + rot[1][2]*qy + rot[2][2]*centerTmp[2];
      // these positions might be out of walls but close, depending on the flatness of the cell plane

      for(size_t k=0; k< cell.numVertex(); ++k){// copying back the original positions
	size_t Vind=cell.vertex(k) -> index();
	for(size_t ii=0; ii< 3; ++ii) 
	  vertexData[Vind][ii]=verticesPosition[k][ii];
      }
      // copynig back the centerCOM
      if(numParameter()==6 && parameter(4)==1) {// centerTriangulation
	cellData[i][variableIndex(2,0)  ]=COMTmp[0];
	cellData[i][variableIndex(2,0)+1]=COMTmp[1];
	cellData[i][variableIndex(2,0)+2]=COMTmp[2];
      }
    
   

      size_t V1W1=(cell.wall(winner.wall1) -> vertex1()) -> index();
      size_t V2W1=(cell.wall(winner.wall1) -> vertex2()) -> index();
      size_t V1W2=(cell.wall(winner.wall2) -> vertex1()) -> index();
      size_t V2W2=(cell.wall(winner.wall2) -> vertex2()) -> index();
    
      // //closest point on wall1 to p 
      std::vector<double> tmpVec(3), wallVec(3);
      for(size_t k=0; k< 3; ++k){
	tmpVec[k]=p[k]-vertexData[V1W1][k];
	wallVec[k]=vertexData[V2W1][k]-vertexData[V1W1][k];
      }


      double normW=std::sqrt(wallVec[0]*wallVec[0]+wallVec[1]*wallVec[1]+wallVec[2]*wallVec[2]);
      for(size_t k=0; k< 3; ++k)
	if(normW!=0)
	  wallVec[k]/=normW;
	else{
	  std::cerr<<" in DivisionShortestPath first wall vector is wrong"<<std::endl;
	  exit(0);
	}
    

      // projecting the vector between first wall vertex and p on wall vector
      double norm1=tmpVec[0]*wallVec[0]+tmpVec[1]*wallVec[1]+tmpVec[2]*wallVec[2];

      // normW=std::sqrt(tmpVec[0]*tmpVec[0]+tmpVec[1]*tmpVec[1]+tmpVec[2]*tmpVec[2]);
      // std::cerr<<std::endl<<" angle 1  "<<norm1/normW<<std::endl;
      if (normW>0)
	for(size_t k=0; k< 3; ++k)      
	  p[k]=vertexData[V1W1][k]+norm1*wallVec[k];    
      else{
	std::cerr<<" in DivisionShortestPath p is wrong"<<std::endl;
	exit(0);    
      }
    
      //closest point on wall2 to q 
      for(size_t k=0; k< 3; ++k){
	tmpVec[k]=q[k]-vertexData[V1W2][k];
	wallVec[k]=vertexData[V2W2][k]-vertexData[V1W2][k];
      }
    
      normW=std::sqrt(wallVec[0]*wallVec[0]+wallVec[1]*wallVec[1]+wallVec[2]*wallVec[2]);
      for(size_t k=0; k< 3; ++k)
	if(normW!=0)
	  wallVec[k]/=normW;
	else{
	  std::cerr<<" in DivisionShortestPath second wall vector is wrong"<<std::endl;
	  exit(0);
	}


      // projecting the vector between second wall vertex and q on wall vector
      norm1=tmpVec[0]*wallVec[0]+tmpVec[1]*wallVec[1]+tmpVec[2]*wallVec[2];
    
      // normW=std::sqrt(tmpVec[0]*tmpVec[0]+tmpVec[1]*tmpVec[1]+tmpVec[2]*tmpVec[2]);
      // std::cerr<<std::endl<<" angle 2  "<<norm1/normW<<std::endl;
      if (normW>0)
	for(size_t k=0; k< 3; ++k)      
	  q[k]=vertexData[V1W2][k]+norm1*wallVec[k];    
      else{
	std::cerr<<" in DivisionShortestPath q is wrong"<<std::endl;
	exit(0);    
      }
    

    
    }

    // std::cerr<<" before:: "<<std::endl;
    // for(size_t k=0; k< cellData[cell.index()].size(); ++k)
    //   std::cerr<<cellData[cell.index()][k]<<"  ";
    // std::cerr<<std::endl;

    if (numVariableIndexLevel() == 2)
      {
	const size_t timeIndex = variableIndex(1, 0);
      
	const double age = cellData[cell.index()][timeIndex];
      
	std::cerr << "Cell age at division is " << age << "\n";
      
	cellData[cell.index()][timeIndex] = 0.0;
      }


  
    if(numParameter()==6 && parameter(4)==1) {// centerTriangulation
      if(parameter(5)==0 || parameter(5)==1)
	T->divideCellCenterTriangulation(&cell, winner.wall1, winner.wall2,
					 variableIndex(2,0),variableIndex(2,1),
					 p, q, cellData, wallData, vertexData,
					 cellDerivs,wallDerivs,vertexDerivs,variableIndex(0),
					 parameter(2),
					 parameter(5));
      else{
	std::cerr<<"parameter(5) should be 0 or 1"
		 <<std::endl;
	exit(-1);
      }
    }
    else
      T->divideCell(&cell, winner.wall1, winner.wall2, 
		    p, q, cellData, wallData, vertexData,
		    cellDerivs, wallDerivs, vertexDerivs, variableIndex(0), 
		    parameter(2));
  
    assert (numWallTmp + 3 == T->numWall());
  

    //  std::cerr<<"  after::  "<<std::endl;
    //  for(size_t k=0; k< cellData[cell.index()].size(); ++k)
    //    std::cerr<<cellData[cell.index()][k]<<"  ";
    //  std::cerr<<std::endl;
    // std::cerr<<std::endl;
    //  for(size_t k=0; k< cellData[(T -> numCell())-1].size(); ++k)
    //    std::cerr<<cellData[(T -> numCell())-1][k]<<"  ";
    //  std::cerr<<std::endl;



    //Check that the division did not mess up the data structure
    //T->checkConnectivity(1);
  }

  std::vector<ShortestPath::Candidate> 
  ShortestPath::getCandidates(Tissue* T, size_t i,
			      DataMatrix &cellData,
			      DataMatrix &wallData,
			      DataMatrix &vertexData,
			      DataMatrix &cellDerivs,
			      DataMatrix &wallDerivs,
			      DataMatrix &vertexDerivs)
  {
    Cell cell = T->cell(i);
  
    assert(cell.numWall() > 1);
  
    std::vector<double> o;
  
    if (parameter(3) == 1)
      {
	o = cell.positionFromVertex(vertexData);
      }
    else
      {
	try
	  {
	    o = cell.randomPositionInCell(vertexData);
	  }
	catch (Cell::FailedToFindRandomPositionInCellException)
	  {
	    return std::vector<Candidate>();
	  }
      }
  
    double ox = o[0];
    double oy = o[1];
  

    std::vector<Candidate> candidates;
  
    for (size_t i = 0; i < cell.numWall() - 1; ++i) {
      for (size_t j = i + 1; j < cell.numWall(); ++j) {
	Wall *wall1 = cell.wall(i);
	Wall *wall2 = cell.wall(j);
	size_t wall1Index = i;
	size_t wall2Index = j;
      
	// 			std::cerr << "i = " << wall1->index() << " : j = " << wall2->index() << std::endl;
	//  			std::cerr << "o = (" << ox << ", " << oy << ")" << std::endl;
      
	double x1x;
	double x1y;
	double x2x;
	double x2y;
      
	double vx;
	double vy;
      
	double x1px;
	double x1py;
	double x2px;
	double x2py;
      
	double ux;
	double uy;
      
	bool flippedVectors;
      
	do {
	  flippedVectors = false;
        
	  x1x = vertexData[wall1->vertex1()->index()][0];
	  x1y = vertexData[wall1->vertex1()->index()][1];
	  x2x = vertexData[wall1->vertex2()->index()][0];
	  x2y = vertexData[wall1->vertex2()->index()][1];
	
	  vx = x2x - x1x;
	  vy = x2y - x1y;
	
	  if (vx * (oy - x1y) - vy * (ox - x1x) > 0) {
	    //std::cerr << "Change v" << std::endl;
	    double tmpx = x1x;
	    double tmpy = x1y;
	    x1x = x2x;
	    x1y = x2y;
	    x2x = tmpx;
	    x2y = tmpy;
	    vx = -vx;
	    vy = -vy;
	  }
        
        
	  x1px = vertexData[wall2->vertex1()->index()][0];
	  x1py = vertexData[wall2->vertex1()->index()][1];
	  x2px = vertexData[wall2->vertex2()->index()][0];
	  x2py = vertexData[wall2->vertex2()->index()][1];
	
	  ux = x2px - x1px;
	  uy = x2py - x1py;
	
	  if (ux * (oy - x1py) - uy * (ox - x1px) < 0) {
	    //  					std::cerr << "Change u" << std::endl;
	    double tmpx = x1px;
	    double tmpy = x1py;
	    x1px = x2px;
	    x1py = x2py;
	    x2px = tmpx;
	    x2py = tmpy;
	    ux = -ux;
	    uy = -uy;
	  }
        
	  if (vx * uy - vy * ux > 0) {
	    //  					std::cerr << "Flipped walls" << std::endl;
	    Wall *tmp = wall1;
	    wall1 = wall2;
	    wall2 = tmp;
	    size_t tmpIndex = wall1Index;
	    wall1Index = wall2Index;
	    wall2Index = tmpIndex;
	    flippedVectors = true;
	  }
	} while (flippedVectors == true);
      
	double wx = ox - x1x;
	double wy = oy - x1y;
	double wpx = ox - x1px;
	double wpy = oy - x1py;
      
	double dvx = wx - ((vx * wx + vy * wy) / (vx * vx + vy * vy)) * vx;
	double dvy = wy - ((vx * wx + vy * wy) / (vx * vx + vy * vy)) * vy;
	double dux = wpx - ((ux * wpx + uy * wpy) / (ux * ux + uy * uy)) * ux;
	double duy = wpy - ((ux * wpx + uy * wpy) / (ux * ux + uy * uy)) * uy;
      
      
	// 			std::cerr << " x1 = (" << x1x << ", " << x1y << ")" << std::endl;
	// 			std::cerr << " x2 = (" << x2x << ", " << x2y << ")" << std::endl;
	// 			std::cerr << " x1p = (" << x1px << ", " << x1py << ")" << std::endl;
	// 			std::cerr << " x2p = (" << x2px << ", " << x2py << ")" << std::endl;
      
	// 			std::cerr << " v = (" << vx << ", " << vy << ")" << std::endl;
	// 			std::cerr << " u = (" << ux << ", " << uy << ")" << std::endl;
	// 			std::cerr << " w = (" << wx << ", " << wy << ")" << std::endl;
	// 			std::cerr << " wp = (" << wpx << ", " << wpy << ")" << std::endl;
	// 			std::cerr << " dv = (" << dvx << ", " << dvy << ")" << std::endl;
	// 			std::cerr << " du = (" << dux << ", " << duy << ")" << std::endl;
      
      
	double A = std::sqrt(dvx * dvx + dvy * dvy);
	double B = std::sqrt(dux * dux + duy * duy);
      
	double sigma = std::acos((vx * ux + vy * uy) / 
				 (std::sqrt(vx * vx + vy * vy) * std::sqrt(ux * ux + uy * uy)));
      
	double alpha = astar(sigma, A, B);
	double beta = myMath::pi() + sigma - alpha;
      
	double t = (vx * wx + vy * wy) / (vx * vx + vy * vy);
	double tp = t + (1.0 / std::sqrt(vx * vx + vy * vy)) * A * 
	  std::sin(alpha - 0.50 * myMath::pi()) / std::sin(alpha);
      
	double s = (ux * wpx + uy * wpy) / (ux * ux + uy * uy);
	double sp = s + (1.0 / std::sqrt(ux * ux + uy * uy)) * B * 
	  std::sin(beta - 0.50 * myMath::pi()) / std::sin(beta);
      
	double px = x1x + tp * vx;
	double py = x1y + tp * vy;
      
	double qx = x1px + sp * ux;
	double qy = x1py + sp * uy;
      
	// 			std::cerr << " sigma = " << sigma << std::endl
	// 				  << " alpha = " << alpha << std::endl
	// 				  << " beta = " << beta << std::endl
	// 				  << " A = " << A << std::endl
	// 				  << " B = " << B << std::endl
	// 				  << " px = " << px << std::endl
	// 				  << " py = " << py << std::endl
	// 				  << " qx = " << qx << std::endl
	// 				  << " qy = " << qy << std::endl
	// 				  << " tp = " << tp << std::endl
	// 				  << " sp = " << sp << std::endl;
      
      
	double distance = std::sqrt((qx - px) * (qx - px) + (qy - py) * (qy - py));
      
	//   			std::cerr << " distance = " << distance << std::endl;
      
	if (tp <= 0.0 || tp >= 1.0 || sp <= 0.0 || sp >= 1.0) {
	  //   				std::cerr << "Discard" << std::endl;
	  continue;
	} else {
	  //  				std::cerr << "Keep" << std::endl;
	  Candidate candidate;
	  candidate.distance = distance;
	  candidate.px = px;
	  candidate.py = py;
	  candidate.qx = qx;
	  candidate.qy = qy;
	  candidate.wall1 = wall1Index;
	  candidate.wall2 = wall2Index;
        
	  candidates.push_back(candidate);
	}
      }
    }	
  
    return candidates;
  }

  double ShortestPath::astar(double sigma, double A, double B)
  {
    double a = 0;
    double b = myMath::pi();
    double e = b - a;
    double u = f(a, sigma, A, B);
    double v = f(b, sigma, A, B);
    double c;
  
    if (myMath::sign(u) == myMath::sign(v)) {
      return 0;
    }
  
    for (size_t k = 0; k < 10; ++k) {
      e = 0.5 * e;
      c = a + e;
      double w = f(c, sigma, A, B);
    
      if (myMath::sign(w) != myMath::sign(u)) {
	b = c;
	v = w;
      } else {
	a = c;
	u = w;
      }
    }
    return c;
  }

  double ShortestPath::f(double a, double sigma, double A, double B)
  {
    double tmp = - A * std::cos(a) / (std::sin(a) * std::sin(a));
    tmp += B * std::cos(myMath::pi() + sigma - a) / (std::sin(sigma - a) * std::sin(sigma - a));
    return tmp;
  }

  STAViaShortestPath::STAViaShortestPath(std::vector<double> &paraValue, 
			     std::vector< std::vector<size_t> > &indValue)
  {
    if ( paraValue.size() != 4 &&  paraValue.size() != 6 ) {
      std::cerr << "DivisionSTAViaShortestPath::DivisionSTAViaShortestPath() "
		<< "Four or six parameters are used V_threshold, Lwall_fraction, "
		<< "Lwall_threshold, and COM (1 = COM, 0 = Random) "
		<< "If six parameters are used, two additional parameters are for " 
		<< "centerTriangulation(1) and double resting length (1: double, 0:single). "
		<<std::endl;
      std::exit(EXIT_FAILURE);
    }
  
    if ((indValue.size() == 2 && indValue[1].size() != 1) || 
	(indValue.size() == 3 && indValue[2].size() != 2)) {
      std::cerr << "DivisionSTAViaShortestPath::DivisionSTAViaShortestPath() "
		<< "First level: Variable indices for volume dependent cell variables are used.\n"
		<< "Second level (optional): Cell time index."
		<< "Third level (if centerTriangulated): first:Com index, second: wall length index"
		<<std::endl;
      exit(EXIT_FAILURE);
    }
  
    setId("DivisionSTAViaShortestPath");
    setNumChange(1);
    setParameter(paraValue);  
    setVariableIndex(indValue);
  
    std::vector<std::string> tmp(numParameter());
    tmp.resize (numParameter());
    tmp[0] = "V_threshold";
    tmp[1] = "Lwall_fraction";
    tmp[2] = "Lwall_threshold";
    tmp[3] = "COM";
    if(numParameter()==6){
      tmp[4] = "centerTriangulationFlag";
      tmp[5] = "doubleLengthFlag";
    }
    setParameterId(tmp);
  }

  int STAViaShortestPath::flag(Tissue *T, size_t i,
			 DataMatrix &cellData,
			 DataMatrix &wallData,
			 DataMatrix &vertexData,
			 DataMatrix &cellDerivs,
			 DataMatrix &wallDerivs,
			 DataMatrix &vertexDerivs)
  {
    if (T->cell(i).calculateVolume(vertexData) > parameter(0))
      {
	return 1;
      } 
    else
      {
	return 0;
      }
  }

  void STAViaShortestPath::update(Tissue* T, size_t i,
			    DataMatrix &cellData,
			    DataMatrix &wallData,
			    DataMatrix &vertexData,
			    DataMatrix &cellDerivs,
			    DataMatrix &wallDerivs,
			    DataMatrix &vertexDerivs)
  {
    Cell &cell = T->cell(i);
  
    // if (vertexData[0].size() != 2)
    //   {
    //     std::cerr << "Division::STAViaShortestPath only supports two dimensions.\n";
    //     std::exit(EXIT_FAILURE);
    //   }
  

    // projecting the cell to the perpendicular plane to the averaged normal
    double rot[3][3]={{1,0,0},{0,1,0},{0,0,1}};
    std::vector<std::vector<double> > verticesPosition;
    verticesPosition.resize(cell.numVertex());
    std::vector<double> centerTmp(3), COMTmp(3); 
  
    if (vertexData[0].size() == 3){
      // storing the vertex positions
      for(size_t k=0; k< cell.numVertex(); ++k) {
	verticesPosition[k].resize(3);
	size_t Vind=cell.vertex(k) -> index();
	verticesPosition[k][0]=vertexData[Vind][0];
	verticesPosition[k][1]=vertexData[Vind][1];
	verticesPosition[k][2]=vertexData[Vind][2];
      }
      // storing COM position    
      if(numParameter()==6 && parameter(4)==1) {// centerTriangulation
	COMTmp[0]=cellData[i][variableIndex(2,0)  ];
	COMTmp[1]=cellData[i][variableIndex(2,0)+1];
	COMTmp[2]=cellData[i][variableIndex(2,0)+2];
      }    
      // Finding the average normal vector to the cell plane
      std::vector<double> normal(3);
      normal[0]=0;
      normal[1]=0;
      normal[2]=0;
    
      //  double area=0;
      for (size_t k=1; k< cell.numWall()-1; ++k) {
	std::vector<double> x0(3),x1(3),x2(3);
	size_t ind1=cell.vertex(0) -> index();
	size_t ind2=cell.vertex(k) -> index();
	size_t ind3=cell.vertex(k+1) -> index();
	//std::cerr<<ind1<<" "<<ind2<<" "<<ind3<<" "<<std::endl;
	//normal to the element
	std::vector<double> temp(3);
	temp[0]=(vertexData[ind2][1]-vertexData[ind1][1])*(vertexData[ind3][2]-vertexData[ind1][2])
	  -(vertexData[ind2][2]-vertexData[ind1][2])*(vertexData[ind3][1]-vertexData[ind1][1]);
	temp[1]=(vertexData[ind2][2]-vertexData[ind1][2])*(vertexData[ind3][0]-vertexData[ind1][0])
	  -(vertexData[ind2][0]-vertexData[ind1][0])*(vertexData[ind3][2]-vertexData[ind1][2]);
	temp[2]=(vertexData[ind2][0]-vertexData[ind1][0])*(vertexData[ind3][1]-vertexData[ind1][1])
	  -(vertexData[ind2][1]-vertexData[ind1][1])*(vertexData[ind3][0]-vertexData[ind1][0]);
    
	normal[0]+=temp[0];
	normal[1]+=temp[1];
	normal[2]+=temp[2];
	//std::cerr<<normal[0]<<" "<<normal[1]<<" "<<normal[2]<<" "<<area<<std::endl;  
      
	//total area of the cell
	//area +=std::sqrt(temp[0]*temp[0]+temp[1]*temp[1]+temp[2]*temp[2])/2;
      }
      double norm=std::sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
      // normalizing the normal vector
      normal[0]/=norm;
      normal[1]/=norm;
      normal[2]/=norm;
      //std::cerr<<normal[0]<<" "<<normal[1]<<" "<<normal[2]<<" "<<area<<std::endl;  
    
    
      //std::cerr<<deltaTet<<" "<<nplane<<std::endl;  
      // rotation matrix for going to the perpendicular plane to thr averaged normal
  
      if(normal[2]!=1){ // if there rotation matrix is not 1
	rot[0][0]=normal[2]+((normal[1]*normal[1])/(normal[2]+1));
	rot[1][1]=normal[2]+((normal[0]*normal[0])/(normal[2]+1));
	rot[0][1]= - (normal[0]*normal[1])/(normal[2]+1) ;
	rot[1][0]= - (normal[1]*normal[0])/(normal[2]+1) ;
      
	rot[0][2]=-normal[0];
	rot[2][0]= normal[0];
	rot[1][2]=-normal[1];
	rot[2][1]= normal[1];
	rot[2][2]= normal[2];
      
      }
      //rotating the cell vertices
      //std::vector<double> tmpVertex(3);
    
      for(size_t k=0; k< cell.numVertex(); ++k){
      
	size_t Vind=cell.vertex(k) -> index();
	for(size_t ii=0; ii< 3; ++ii) {
	  vertexData[Vind][ii]=0;
	  for(size_t jj=0; jj< 3; ++jj){ 
	    vertexData[Vind][ii]+=rot[ii][jj]*verticesPosition[k][jj];
	  }
	}
      }

      // rotating the center is cell is centertriangulated
      if(numParameter()==6 && parameter(4)==1) {// centerTriangulation
	cellData[i][variableIndex(2,0)  ]=rot[0][0]*COMTmp[0]+rot[0][1]*COMTmp[1]+rot[0][2]*COMTmp[2];
	cellData[i][variableIndex(2,0)+1]=rot[1][0]*COMTmp[0]+rot[1][1]*COMTmp[1]+rot[1][2]*COMTmp[2];
	cellData[i][variableIndex(2,0)+2]=rot[2][0]*COMTmp[0]+rot[2][1]*COMTmp[1]+rot[2][2]*COMTmp[2];
      }

      centerTmp=cell.positionFromVertex(vertexData);


      // std::cerr<<"normal "<<normal[0]<<" "<<normal[1]<<" "<<normal[2]<<std::endl;
      //  for(size_t k=0; k< cell.numVertex(); ++k){
      //    size_t Vind=cell.vertex(k) -> index();
      //    std::cerr<<"vertex "<<k<<" "
      //             <<vertexData[Vind][0]<<" "
      //             <<vertexData[Vind][1]<<" "
      //             <<vertexData[Vind][2]<<std::endl;
      //  }
      //  std::cerr<<"cell center "
      //           <<cellData[i][variableIndex(2,0)  ]<<" "
      //           <<cellData[i][variableIndex(2,0)+1]<<" "
      //           <<cellData[i][variableIndex(2,0)+2]<<std::endl;
      //  std::cerr<<"cell centroid "
      //           <<centerTmp[0]<<" "
      //           <<centerTmp[1]<<" "
      //           <<centerTmp[2]<<std::endl;  
    }  
    std::vector<Candidate> candidates = 
      getCandidates(T, i, cellData, wallData, vertexData, cellDerivs, wallDerivs, vertexDerivs);
  
    if (candidates.size() == 0)
      {
	return;
      }
  
    Candidate winner = { std::numeric_limits<double>::max(), 0, 0, 0, 0, 0, 0 };
    for (size_t i = 0; i < candidates.size(); ++i) {
      if (candidates[i].distance < winner.distance) {
	winner = candidates[i];
      }
    }
  
    // std::cerr << "Winner: " << std::endl
    //           << " distance = " << winner.distance << std::endl
    //           << " p = (" << winner.px << ", " << winner.py << ")" << std::endl
    //           << " q = (" << winner.qx << ", " << winner.qy << ")" << std::endl;
  
    size_t numWallTmp = wallData.size();
    assert(numWallTmp == T->numWall());
  
  
    std::vector<double> p(3);
    p[0] = winner.px;
    p[1] = winner.py;
    std::vector<double> q(3);
    q[0] = winner.qx;
    q[1] = winner.qy;
  
    if (vertexData[0].size() == 3){

      double px=p[0],py=p[1],qx=q[0], qy=q[1];
      // p.resize(3);
      // q.resize(3);




      // size_t V1W1=(cell.wall(winner.wall1) -> vertex1()) -> index();
      // size_t V2W1=(cell.wall(winner.wall1) -> vertex2()) -> index();
      // size_t V1W2=(cell.wall(winner.wall2) -> vertex1()) -> index();
      // size_t V2W2=(cell.wall(winner.wall2) -> vertex2()) -> index();
    
      // std::cerr<<std::endl;
      // std::cerr<<" v1w1 "<<vertexData[V1W1][0]<<"  "<<vertexData[V1W1][1]<<"  "<<vertexData[V1W1][2]<<std::endl;
      // std::cerr<<" v2w1 "<<vertexData[V2W1][0]<<"  "<<vertexData[V2W1][1]<<"  "<<vertexData[V2W1][2]<<std::endl;
      // std::cerr<<" v1w2 "<<vertexData[V1W2][0]<<"  "<<vertexData[V1W2][1]<<"  "<<vertexData[V1W2][2]<<std::endl;
      // std::cerr<<" v2w2 "<<vertexData[V2W2][0]<<"  "<<vertexData[V2W2][1]<<"  "<<vertexData[V2W2][2]<<std::endl;
      // std::cerr<<" p "<<p[0]<<"  "<<p[1]<<"  "<<p[2]<<std::endl;
      // std::cerr<<" q "<<q[0]<<"  "<<q[1]<<"  "<<q[2]<<std::endl;
      // std::cerr<<" com "<<COMTmp[0]<<"  "<<COMTmp[1]<<"  "<<COMTmp[2]<<std::endl;

      // //testing the rotation matrix    
      // std::vector<std::vector<double> > rotTemp(3);
      // for(size_t iii=0; iii< 3; ++iii)
      //   rotTemp[iii].resize(3);
      // for(size_t iii=0; iii< 3; ++iii)
      //   for(size_t jjj=0; jjj<3; ++jjj)
      //     for(size_t kkk=0; kkk<3; ++kkk)
      //       rotTemp[iii][jjj]+=rot[iii][kkk]*rot[jjj][kkk];
      // for(size_t iii=0; iii< 3; ++iii)
      //   std::cerr<<"rot   "<<rotTemp[iii][0]<<" "<<rotTemp[iii][1]<<" "<<rotTemp[iii][2]<<std::endl;
    
      //rotating back the found positions for division points
    
      p[0]= rot[0][0]*px + rot[1][0]*py + rot[2][0]*centerTmp[2];
      p[1]= rot[0][1]*px + rot[1][1]*py + rot[2][1]*centerTmp[2];
      p[2]= rot[0][2]*px + rot[1][2]*py + rot[2][2]*centerTmp[2];
    
      q[0]= rot[0][0]*qx + rot[1][0]*qy + rot[2][0]*centerTmp[2];
      q[1]= rot[0][1]*qx + rot[1][1]*qy + rot[2][1]*centerTmp[2];
      q[2]= rot[0][2]*qx + rot[1][2]*qy + rot[2][2]*centerTmp[2];
      // these positions might be out of walls but close, depending on the flatness of the cell plane

      for(size_t k=0; k< cell.numVertex(); ++k){// copying back the original positions
	size_t Vind=cell.vertex(k) -> index();
	for(size_t ii=0; ii< 3; ++ii) 
	  vertexData[Vind][ii]=verticesPosition[k][ii];
      }
      // copynig back the centerCOM
      if(numParameter()==6 && parameter(4)==1) {// centerTriangulation
	cellData[i][variableIndex(2,0)  ]=COMTmp[0];
	cellData[i][variableIndex(2,0)+1]=COMTmp[1];
	cellData[i][variableIndex(2,0)+2]=COMTmp[2];
      }
    
   

      size_t V1W1=(cell.wall(winner.wall1) -> vertex1()) -> index();
      size_t V2W1=(cell.wall(winner.wall1) -> vertex2()) -> index();
      size_t V1W2=(cell.wall(winner.wall2) -> vertex1()) -> index();
      size_t V2W2=(cell.wall(winner.wall2) -> vertex2()) -> index();
    
      // //closest point on wall1 to p 
      std::vector<double> tmpVec(3), wallVec(3);
      for(size_t k=0; k< 3; ++k){
	tmpVec[k]=p[k]-vertexData[V1W1][k];
	wallVec[k]=vertexData[V2W1][k]-vertexData[V1W1][k];
      }


      double normW=std::sqrt(wallVec[0]*wallVec[0]+wallVec[1]*wallVec[1]+wallVec[2]*wallVec[2]);
      for(size_t k=0; k< 3; ++k)
	if(normW!=0)
	  wallVec[k]/=normW;
	else{
	  std::cerr<<" in DivisionSTAViaShortestPath first wall vector is wrong"<<std::endl;
	  exit(0);
	}
    

      // projecting the vector between first wall vertex and p on wall vector
      double norm1=tmpVec[0]*wallVec[0]+tmpVec[1]*wallVec[1]+tmpVec[2]*wallVec[2];

      // normW=std::sqrt(tmpVec[0]*tmpVec[0]+tmpVec[1]*tmpVec[1]+tmpVec[2]*tmpVec[2]);
      // std::cerr<<std::endl<<" angle 1  "<<norm1/normW<<std::endl;
      if (normW>0)
	for(size_t k=0; k< 3; ++k)      
	  p[k]=vertexData[V1W1][k]+norm1*wallVec[k];    
      else{
	std::cerr<<" in DivisionSTAViaShortestPath p is wrong"<<std::endl;
	exit(0);    
      }
    
      //closest point on wall2 to q 
      for(size_t k=0; k< 3; ++k){
	tmpVec[k]=q[k]-vertexData[V1W2][k];
	wallVec[k]=vertexData[V2W2][k]-vertexData[V1W2][k];
      }
    
      normW=std::sqrt(wallVec[0]*wallVec[0]+wallVec[1]*wallVec[1]+wallVec[2]*wallVec[2]);
      for(size_t k=0; k< 3; ++k)
	if(normW!=0)
	  wallVec[k]/=normW;
	else{
	  std::cerr<<" in DivisionSTAViaShortestPath second wall vector is wrong"<<std::endl;
	  exit(0);
	}


      // projecting the vector between second wall vertex and q on wall vector
      norm1=tmpVec[0]*wallVec[0]+tmpVec[1]*wallVec[1]+tmpVec[2]*wallVec[2];
    
      // normW=std::sqrt(tmpVec[0]*tmpVec[0]+tmpVec[1]*tmpVec[1]+tmpVec[2]*tmpVec[2]);
      // std::cerr<<std::endl<<" angle 2  "<<norm1/normW<<std::endl;
      if (normW>0)
	for(size_t k=0; k< 3; ++k)      
	  q[k]=vertexData[V1W2][k]+norm1*wallVec[k];    
      else{
	std::cerr<<" in DivisionSTAViaShortestPath q is wrong"<<std::endl;
	exit(0);    
      }
    

    
    }

    // std::cerr<<" before:: "<<std::endl;
    // for(size_t k=0; k< cellData[cell.index()].size(); ++k)
    //   std::cerr<<cellData[cell.index()][k]<<"  ";
    // std::cerr<<std::endl;

    if (numVariableIndexLevel() == 2)
      {
	const size_t timeIndex = variableIndex(1, 0);
      
	const double age = cellData[cell.index()][timeIndex];
      
	std::cerr << "Cell age at division is " << age << "\n";
      
	cellData[cell.index()][timeIndex] = 0.0;
      }


  
    if(numParameter()==6 && parameter(4)==1) {// centerTriangulation
      if(parameter(5)==0 || parameter(5)==1)
	T->divideCellCenterTriangulation(&cell, winner.wall1, winner.wall2,
					 variableIndex(2,0),variableIndex(2,1),
					 p, q, cellData, wallData, vertexData,
					 cellDerivs,wallDerivs,vertexDerivs,variableIndex(0),
					 parameter(2),
					 parameter(5));
      else{
	std::cerr<<"parameter(5) should be 0 or 1"
		 <<std::endl;
	exit(-1);
      }
    }
    else
      T->divideCell(&cell, winner.wall1, winner.wall2, 
		    p, q, cellData, wallData, vertexData,
		    cellDerivs, wallDerivs, vertexDerivs, variableIndex(0), 
		    parameter(2));
  
    assert (numWallTmp + 3 == T->numWall());
  
    //Change length of new wall between the divided daugther cells
    wallData[numWallTmp][0] *= parameter(1);
  

    //  std::cerr<<"  after::  "<<std::endl;
    //  for(size_t k=0; k< cellData[cell.index()].size(); ++k)
    //    std::cerr<<cellData[cell.index()][k]<<"  ";
    //  std::cerr<<std::endl;
    // std::cerr<<std::endl;
    //  for(size_t k=0; k< cellData[(T -> numCell())-1].size(); ++k)
    //    std::cerr<<cellData[(T -> numCell())-1][k]<<"  ";
    //  std::cerr<<std::endl;



    //Check that the division did not mess up the data structure
    //T->checkConnectivity(1);
  }

  std::vector<STAViaShortestPath::Candidate> 
  STAViaShortestPath::getCandidates(Tissue* T, size_t i,
			      DataMatrix &cellData,
			      DataMatrix &wallData,
			      DataMatrix &vertexData,
			      DataMatrix &cellDerivs,
			      DataMatrix &wallDerivs,
			      DataMatrix &vertexDerivs)
  {
    Cell cell = T->cell(i);
  
    assert(cell.numWall() > 1);
  
    std::vector<double> o;
  
    if (parameter(3) == 1)
      {
	o = cell.positionFromVertex(vertexData);
      }
    else
      {
	try
	  {
	    o = cell.randomPositionInCell(vertexData);
	  }
	catch (Cell::FailedToFindRandomPositionInCellException)
	  {
	    return std::vector<Candidate>();
	  }
      }
  
    double ox = o[0];
    double oy = o[1];
  

    std::vector<Candidate> candidates;
  
    for (size_t i = 0; i < cell.numWall() - 1; ++i) {
      for (size_t j = i + 1; j < cell.numWall(); ++j) {
	Wall *wall1 = cell.wall(i);
	Wall *wall2 = cell.wall(j);
	size_t wall1Index = i;
	size_t wall2Index = j;
      
	// 			std::cerr << "i = " << wall1->index() << " : j = " << wall2->index() << std::endl;
	//  			std::cerr << "o = (" << ox << ", " << oy << ")" << std::endl;
      
	double x1x;
	double x1y;
	double x2x;
	double x2y;
      
	double vx;
	double vy;
      
	double x1px;
	double x1py;
	double x2px;
	double x2py;
      
	double ux;
	double uy;
      
	bool flippedVectors;
      
	do {
	  flippedVectors = false;
        
	  x1x = vertexData[wall1->vertex1()->index()][0];
	  x1y = vertexData[wall1->vertex1()->index()][1];
	  x2x = vertexData[wall1->vertex2()->index()][0];
	  x2y = vertexData[wall1->vertex2()->index()][1];
	
	  vx = x2x - x1x;
	  vy = x2y - x1y;
	
	  if (vx * (oy - x1y) - vy * (ox - x1x) > 0) {
	    //std::cerr << "Change v" << std::endl;
	    double tmpx = x1x;
	    double tmpy = x1y;
	    x1x = x2x;
	    x1y = x2y;
	    x2x = tmpx;
	    x2y = tmpy;
	    vx = -vx;
	    vy = -vy;
	  }
        
        
	  x1px = vertexData[wall2->vertex1()->index()][0];
	  x1py = vertexData[wall2->vertex1()->index()][1];
	  x2px = vertexData[wall2->vertex2()->index()][0];
	  x2py = vertexData[wall2->vertex2()->index()][1];
	
	  ux = x2px - x1px;
	  uy = x2py - x1py;
	
	  if (ux * (oy - x1py) - uy * (ox - x1px) < 0) {
	    //  					std::cerr << "Change u" << std::endl;
	    double tmpx = x1px;
	    double tmpy = x1py;
	    x1px = x2px;
	    x1py = x2py;
	    x2px = tmpx;
	    x2py = tmpy;
	    ux = -ux;
	    uy = -uy;
	  }
        
	  if (vx * uy - vy * ux > 0) {
	    //  					std::cerr << "Flipped walls" << std::endl;
	    Wall *tmp = wall1;
	    wall1 = wall2;
	    wall2 = tmp;
	    size_t tmpIndex = wall1Index;
	    wall1Index = wall2Index;
	    wall2Index = tmpIndex;
	    flippedVectors = true;
	  }
	} while (flippedVectors == true);
      
	double wx = ox - x1x;
	double wy = oy - x1y;
	double wpx = ox - x1px;
	double wpy = oy - x1py;
      
	double dvx = wx - ((vx * wx + vy * wy) / (vx * vx + vy * vy)) * vx;
	double dvy = wy - ((vx * wx + vy * wy) / (vx * vx + vy * vy)) * vy;
	double dux = wpx - ((ux * wpx + uy * wpy) / (ux * ux + uy * uy)) * ux;
	double duy = wpy - ((ux * wpx + uy * wpy) / (ux * ux + uy * uy)) * uy;
      
      
	// 			std::cerr << " x1 = (" << x1x << ", " << x1y << ")" << std::endl;
	// 			std::cerr << " x2 = (" << x2x << ", " << x2y << ")" << std::endl;
	// 			std::cerr << " x1p = (" << x1px << ", " << x1py << ")" << std::endl;
	// 			std::cerr << " x2p = (" << x2px << ", " << x2py << ")" << std::endl;
      
	// 			std::cerr << " v = (" << vx << ", " << vy << ")" << std::endl;
	// 			std::cerr << " u = (" << ux << ", " << uy << ")" << std::endl;
	// 			std::cerr << " w = (" << wx << ", " << wy << ")" << std::endl;
	// 			std::cerr << " wp = (" << wpx << ", " << wpy << ")" << std::endl;
	// 			std::cerr << " dv = (" << dvx << ", " << dvy << ")" << std::endl;
	// 			std::cerr << " du = (" << dux << ", " << duy << ")" << std::endl;
      
      
	double A = std::sqrt(dvx * dvx + dvy * dvy);
	double B = std::sqrt(dux * dux + duy * duy);
      
	double sigma = std::acos((vx * ux + vy * uy) / 
				 (std::sqrt(vx * vx + vy * vy) * std::sqrt(ux * ux + uy * uy)));
      
	double alpha = astar(sigma, A, B);
	double beta = myMath::pi() + sigma - alpha;
      
	double t = (vx * wx + vy * wy) / (vx * vx + vy * vy);
	double tp = t + (1.0 / std::sqrt(vx * vx + vy * vy)) * A * 
	  std::sin(alpha - 0.50 * myMath::pi()) / std::sin(alpha);
      
	double s = (ux * wpx + uy * wpy) / (ux * ux + uy * uy);
	double sp = s + (1.0 / std::sqrt(ux * ux + uy * uy)) * B * 
	  std::sin(beta - 0.50 * myMath::pi()) / std::sin(beta);
      
	double px = x1x + tp * vx;
	double py = x1y + tp * vy;
      
	double qx = x1px + sp * ux;
	double qy = x1py + sp * uy;
      
	// 			std::cerr << " sigma = " << sigma << std::endl
	// 				  << " alpha = " << alpha << std::endl
	// 				  << " beta = " << beta << std::endl
	// 				  << " A = " << A << std::endl
	// 				  << " B = " << B << std::endl
	// 				  << " px = " << px << std::endl
	// 				  << " py = " << py << std::endl
	// 				  << " qx = " << qx << std::endl
	// 				  << " qy = " << qy << std::endl
	// 				  << " tp = " << tp << std::endl
	// 				  << " sp = " << sp << std::endl;
      
      
	double distance = std::sqrt((qx - px) * (qx - px) + (qy - py) * (qy - py));
      
	//   			std::cerr << " distance = " << distance << std::endl;
      
	if (tp <= 0.0 || tp >= 1.0 || sp <= 0.0 || sp >= 1.0) {
	  //   				std::cerr << "Discard" << std::endl;
	  continue;
	} else {
	  //  				std::cerr << "Keep" << std::endl;
	  Candidate candidate;
	  candidate.distance = distance;
	  candidate.px = px;
	  candidate.py = py;
	  candidate.qx = qx;
	  candidate.qy = qy;
	  candidate.wall1 = wall1Index;
	  candidate.wall2 = wall2Index;
        
	  candidates.push_back(candidate);
	}
      }
    }	
  
    return candidates;
  }

  double STAViaShortestPath::astar(double sigma, double A, double B)
  {
    double a = 0;
    double b = myMath::pi();
    double e = b - a;
    double u = f(a, sigma, A, B);
    double v = f(b, sigma, A, B);
    double c;
  
    if (myMath::sign(u) == myMath::sign(v)) {
      return 0;
    }
  
    for (size_t k = 0; k < 10; ++k) {
      e = 0.5 * e;
      c = a + e;
      double w = f(c, sigma, A, B);
    
      if (myMath::sign(w) != myMath::sign(u)) {
	b = c;
	v = w;
      } else {
	a = c;
	u = w;
      }
    }
    return c;
  }

  double STAViaShortestPath::f(double a, double sigma, double A, double B)
  {
    double tmp = - A * std::cos(a) / (std::sin(a) * std::sin(a));
    tmp += B * std::cos(myMath::pi() + sigma - a) / (std::sin(sigma - a) * std::sin(sigma - a));
    return tmp;
  }

  FlagResetShortestPath::FlagResetShortestPath(std::vector<double> &paraValue, 
           std::vector< std::vector<size_t> > &indValue)
  {
    if ( paraValue.size() != 4 &&  paraValue.size() != 6 ) {
      std::cerr << "DivisionFlagResetShortestPath::DivisionFlagResetShortestPath() "
    << "Four or six parameters are used V_threshold, Lwall_fraction, "
    << "Lwall_threshold, and COM (1 = COM, 0 = Random) "
    << "If six parameters are used, two additional parameters are for " 
    << "centerTriangulation(1) and double resting length (1: double, 0:single). "
    <<std::endl;
      std::exit(EXIT_FAILURE);
    }
  
    if ((indValue.size() == 2 && indValue[1].size() != 1) || 
  (indValue.size() == 3 && indValue[2].size() != 2)) {
      std::cerr << "DivisionFlagResetShortestPath::DivisionFlagResetShortestPath() "
    << "First level: Variable indices for volume dependent cell variables are used.\n"
    << "Second level (optional): Cell time index."
    << "Third level (if centerTriangulated): first:Com index, second: wall length index"
    <<std::endl;
      exit(EXIT_FAILURE);
    }
  
    setId("DivisionFlagResetShortestPath");
    setNumChange(1);
    setParameter(paraValue);  
    setVariableIndex(indValue);
  
    std::vector<std::string> tmp(numParameter());
    tmp.resize (numParameter());
    tmp[0] = "V_threshold";
    tmp[1] = "Lwall_fraction";
    tmp[2] = "Lwall_threshold";
    tmp[3] = "COM";
    if(numParameter()==6){
      tmp[4] = "centerTriangulationFlag";
      tmp[5] = "doubleLengthFlag";
    }
    setParameterId(tmp);
  }

  int FlagResetShortestPath::flag(Tissue *T, size_t i,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs)
  { if (cellData[i][variableIndex(3,0)]==1)
      {
  return 1;
      } 
    else
      {
  return 0;
      }
  }

  void FlagResetShortestPath::update(Tissue* T, size_t i,
          DataMatrix &cellData,
          DataMatrix &wallData,
          DataMatrix &vertexData,
          DataMatrix &cellDerivs,
          DataMatrix &wallDerivs,
          DataMatrix &vertexDerivs)
  {
    Cell &cell = T->cell(i);
  
    // if (vertexData[0].size() != 2)
    //   {
    //     std::cerr << "Division::FlagResetShortestPath only supports two dimensions.\n";
    //     std::exit(EXIT_FAILURE);
    //   }
  

    // projecting the cell to the perpendicular plane to the averaged normal
    double rot[3][3]={{1,0,0},{0,1,0},{0,0,1}};
    std::vector<std::vector<double> > verticesPosition;
    verticesPosition.resize(cell.numVertex());
    std::vector<double> centerTmp(3), COMTmp(3); 
  
    if (vertexData[0].size() == 3){
      // storing the vertex positions
      for(size_t k=0; k< cell.numVertex(); ++k) {
  verticesPosition[k].resize(3);
  size_t Vind=cell.vertex(k) -> index();
  verticesPosition[k][0]=vertexData[Vind][0];
  verticesPosition[k][1]=vertexData[Vind][1];
  verticesPosition[k][2]=vertexData[Vind][2];
      }
      // storing COM position    
      if(numParameter()==6 && parameter(4)==1) {// centerTriangulation
  COMTmp[0]=cellData[i][variableIndex(2,0)  ];
  COMTmp[1]=cellData[i][variableIndex(2,0)+1];
  COMTmp[2]=cellData[i][variableIndex(2,0)+2];
      }    
      // Finding the average normal vector to the cell plane
      std::vector<double> normal(3);
      normal[0]=0;
      normal[1]=0;
      normal[2]=0;
    
      //  double area=0;
      for (size_t k=1; k< cell.numWall()-1; ++k) {
  std::vector<double> x0(3),x1(3),x2(3);
  size_t ind1=cell.vertex(0) -> index();
  size_t ind2=cell.vertex(k) -> index();
  size_t ind3=cell.vertex(k+1) -> index();
  //std::cerr<<ind1<<" "<<ind2<<" "<<ind3<<" "<<std::endl;
  //normal to the element
  std::vector<double> temp(3);
  temp[0]=(vertexData[ind2][1]-vertexData[ind1][1])*(vertexData[ind3][2]-vertexData[ind1][2])
    -(vertexData[ind2][2]-vertexData[ind1][2])*(vertexData[ind3][1]-vertexData[ind1][1]);
  temp[1]=(vertexData[ind2][2]-vertexData[ind1][2])*(vertexData[ind3][0]-vertexData[ind1][0])
    -(vertexData[ind2][0]-vertexData[ind1][0])*(vertexData[ind3][2]-vertexData[ind1][2]);
  temp[2]=(vertexData[ind2][0]-vertexData[ind1][0])*(vertexData[ind3][1]-vertexData[ind1][1])
    -(vertexData[ind2][1]-vertexData[ind1][1])*(vertexData[ind3][0]-vertexData[ind1][0]);
    
  normal[0]+=temp[0];
  normal[1]+=temp[1];
  normal[2]+=temp[2];
  //std::cerr<<normal[0]<<" "<<normal[1]<<" "<<normal[2]<<" "<<area<<std::endl;  
      
  //total area of the cell
  //area +=std::sqrt(temp[0]*temp[0]+temp[1]*temp[1]+temp[2]*temp[2])/2;
      }
      double norm=std::sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
      // normalizing the normal vector
      normal[0]/=norm;
      normal[1]/=norm;
      normal[2]/=norm;
      //std::cerr<<normal[0]<<" "<<normal[1]<<" "<<normal[2]<<" "<<area<<std::endl;  
    
    
      //std::cerr<<deltaTet<<" "<<nplane<<std::endl;  
      // rotation matrix for going to the perpendicular plane to thr averaged normal
  
      if(normal[2]!=1){ // if there rotation matrix is not 1
  rot[0][0]=normal[2]+((normal[1]*normal[1])/(normal[2]+1));
  rot[1][1]=normal[2]+((normal[0]*normal[0])/(normal[2]+1));
  rot[0][1]= - (normal[0]*normal[1])/(normal[2]+1) ;
  rot[1][0]= - (normal[1]*normal[0])/(normal[2]+1) ;
      
  rot[0][2]=-normal[0];
  rot[2][0]= normal[0];
  rot[1][2]=-normal[1];
  rot[2][1]= normal[1];
  rot[2][2]= normal[2];
      
      }
      //rotating the cell vertices
      //std::vector<double> tmpVertex(3);
    
      for(size_t k=0; k< cell.numVertex(); ++k){
      
  size_t Vind=cell.vertex(k) -> index();
  for(size_t ii=0; ii< 3; ++ii) {
    vertexData[Vind][ii]=0;
    for(size_t jj=0; jj< 3; ++jj){ 
      vertexData[Vind][ii]+=rot[ii][jj]*verticesPosition[k][jj];
    }
  }
      }

      // rotating the center is cell is centertriangulated
      if(numParameter()==6 && parameter(4)==1) {// centerTriangulation
  cellData[i][variableIndex(2,0)  ]=rot[0][0]*COMTmp[0]+rot[0][1]*COMTmp[1]+rot[0][2]*COMTmp[2];
  cellData[i][variableIndex(2,0)+1]=rot[1][0]*COMTmp[0]+rot[1][1]*COMTmp[1]+rot[1][2]*COMTmp[2];
  cellData[i][variableIndex(2,0)+2]=rot[2][0]*COMTmp[0]+rot[2][1]*COMTmp[1]+rot[2][2]*COMTmp[2];
      }

      centerTmp=cell.positionFromVertex(vertexData);


      // std::cerr<<"normal "<<normal[0]<<" "<<normal[1]<<" "<<normal[2]<<std::endl;
      //  for(size_t k=0; k< cell.numVertex(); ++k){
      //    size_t Vind=cell.vertex(k) -> index();
      //    std::cerr<<"vertex "<<k<<" "
      //             <<vertexData[Vind][0]<<" "
      //             <<vertexData[Vind][1]<<" "
      //             <<vertexData[Vind][2]<<std::endl;
      //  }
      //  std::cerr<<"cell center "
      //           <<cellData[i][variableIndex(2,0)  ]<<" "
      //           <<cellData[i][variableIndex(2,0)+1]<<" "
      //           <<cellData[i][variableIndex(2,0)+2]<<std::endl;
      //  std::cerr<<"cell centroid "
      //           <<centerTmp[0]<<" "
      //           <<centerTmp[1]<<" "
      //           <<centerTmp[2]<<std::endl;

    
  
  
    }  
    std::vector<Candidate> candidates = 
      getCandidates(T, i, cellData, wallData, vertexData, cellDerivs, wallDerivs, vertexDerivs);
    
    if (candidates.size() == 0)
      {
	return;
      }
    
    Candidate winner = { std::numeric_limits<double>::max(), 0, 0, 0, 0, 0, 0 };
    for (size_t i = 0; i < candidates.size(); ++i) {
      if (candidates[i].distance < winner.distance) {
	winner = candidates[i];
      }
    }
    
    // std::cerr << "Winner: " << std::endl
    //           << " distance = " << winner.distance << std::endl
    //           << " p = (" << winner.px << ", " << winner.py << ")" << std::endl
    //           << " q = (" << winner.qx << ", " << winner.qy << ")" << std::endl;
  
    size_t numWallTmp = wallData.size();
    assert(numWallTmp == T->numWall());
  
  
    std::vector<double> p(3);
    p[0] = winner.px;
    p[1] = winner.py;
    std::vector<double> q(3);
    q[0] = winner.qx;
    q[1] = winner.qy;
  
    if (vertexData[0].size() == 3){

      double px=p[0],py=p[1],qx=q[0], qy=q[1];
      // p.resize(3);
      // q.resize(3);




      // size_t V1W1=(cell.wall(winner.wall1) -> vertex1()) -> index();
      // size_t V2W1=(cell.wall(winner.wall1) -> vertex2()) -> index();
      // size_t V1W2=(cell.wall(winner.wall2) -> vertex1()) -> index();
      // size_t V2W2=(cell.wall(winner.wall2) -> vertex2()) -> index();
    
      // std::cerr<<std::endl;
      // std::cerr<<" v1w1 "<<vertexData[V1W1][0]<<"  "<<vertexData[V1W1][1]<<"  "<<vertexData[V1W1][2]<<std::endl;
      // std::cerr<<" v2w1 "<<vertexData[V2W1][0]<<"  "<<vertexData[V2W1][1]<<"  "<<vertexData[V2W1][2]<<std::endl;
      // std::cerr<<" v1w2 "<<vertexData[V1W2][0]<<"  "<<vertexData[V1W2][1]<<"  "<<vertexData[V1W2][2]<<std::endl;
      // std::cerr<<" v2w2 "<<vertexData[V2W2][0]<<"  "<<vertexData[V2W2][1]<<"  "<<vertexData[V2W2][2]<<std::endl;
      // std::cerr<<" p "<<p[0]<<"  "<<p[1]<<"  "<<p[2]<<std::endl;
      // std::cerr<<" q "<<q[0]<<"  "<<q[1]<<"  "<<q[2]<<std::endl;
      // std::cerr<<" com "<<COMTmp[0]<<"  "<<COMTmp[1]<<"  "<<COMTmp[2]<<std::endl;

      // //testing the rotation matrix    
      // std::vector<std::vector<double> > rotTemp(3);
      // for(size_t iii=0; iii< 3; ++iii)
      //   rotTemp[iii].resize(3);
      // for(size_t iii=0; iii< 3; ++iii)
      //   for(size_t jjj=0; jjj<3; ++jjj)
      //     for(size_t kkk=0; kkk<3; ++kkk)
      //       rotTemp[iii][jjj]+=rot[iii][kkk]*rot[jjj][kkk];
      // for(size_t iii=0; iii< 3; ++iii)
      //   std::cerr<<"rot   "<<rotTemp[iii][0]<<" "<<rotTemp[iii][1]<<" "<<rotTemp[iii][2]<<std::endl;
    
      //rotating back the found positions for division points
    
      p[0]= rot[0][0]*px + rot[1][0]*py + rot[2][0]*centerTmp[2];
      p[1]= rot[0][1]*px + rot[1][1]*py + rot[2][1]*centerTmp[2];
      p[2]= rot[0][2]*px + rot[1][2]*py + rot[2][2]*centerTmp[2];
    
      q[0]= rot[0][0]*qx + rot[1][0]*qy + rot[2][0]*centerTmp[2];
      q[1]= rot[0][1]*qx + rot[1][1]*qy + rot[2][1]*centerTmp[2];
      q[2]= rot[0][2]*qx + rot[1][2]*qy + rot[2][2]*centerTmp[2];
      // these positions might be out of walls but close, depending on the flatness of the cell plane

      for(size_t k=0; k< cell.numVertex(); ++k){// copying back the original positions
  size_t Vind=cell.vertex(k) -> index();
  for(size_t ii=0; ii< 3; ++ii) 
    vertexData[Vind][ii]=verticesPosition[k][ii];
      }
      // copynig back the centerCOM
      if(numParameter()==6 && parameter(4)==1) {// centerTriangulation
  cellData[i][variableIndex(2,0)  ]=COMTmp[0];
  cellData[i][variableIndex(2,0)+1]=COMTmp[1];
  cellData[i][variableIndex(2,0)+2]=COMTmp[2];
      }
    
   

      size_t V1W1=(cell.wall(winner.wall1) -> vertex1()) -> index();
      size_t V2W1=(cell.wall(winner.wall1) -> vertex2()) -> index();
      size_t V1W2=(cell.wall(winner.wall2) -> vertex1()) -> index();
      size_t V2W2=(cell.wall(winner.wall2) -> vertex2()) -> index();
    
      // //closest point on wall1 to p 
      std::vector<double> tmpVec(3), wallVec(3);
      for(size_t k=0; k< 3; ++k){
  tmpVec[k]=p[k]-vertexData[V1W1][k];
  wallVec[k]=vertexData[V2W1][k]-vertexData[V1W1][k];
      }


      double normW=std::sqrt(wallVec[0]*wallVec[0]+wallVec[1]*wallVec[1]+wallVec[2]*wallVec[2]);
      for(size_t k=0; k< 3; ++k)
  if(normW!=0)
    wallVec[k]/=normW;
  else{
    std::cerr<<" in DivisionFlagResetShortestPath first wall vector is wrong"<<std::endl;
    exit(0);
  }
    

      // projecting the vector between first wall vertex and p on wall vector
      double norm1=tmpVec[0]*wallVec[0]+tmpVec[1]*wallVec[1]+tmpVec[2]*wallVec[2];

      // normW=std::sqrt(tmpVec[0]*tmpVec[0]+tmpVec[1]*tmpVec[1]+tmpVec[2]*tmpVec[2]);
      // std::cerr<<std::endl<<" angle 1  "<<norm1/normW<<std::endl;
      if (normW>0)
  for(size_t k=0; k< 3; ++k)      
    p[k]=vertexData[V1W1][k]+norm1*wallVec[k];    
      else{
  std::cerr<<" in DivisionShortestPath p is wrong"<<std::endl;
  exit(0);    
      }
    
      //closest point on wall2 to q 
      for(size_t k=0; k< 3; ++k){
  tmpVec[k]=q[k]-vertexData[V1W2][k];
  wallVec[k]=vertexData[V2W2][k]-vertexData[V1W2][k];
      }
    
      normW=std::sqrt(wallVec[0]*wallVec[0]+wallVec[1]*wallVec[1]+wallVec[2]*wallVec[2]);
      for(size_t k=0; k< 3; ++k)
  if(normW!=0)
    wallVec[k]/=normW;
  else{
    std::cerr<<" in DivisionFlagResetShortestPath second wall vector is wrong"<<std::endl;
    exit(0);
  }


      // projecting the vector between second wall vertex and q on wall vector
      norm1=tmpVec[0]*wallVec[0]+tmpVec[1]*wallVec[1]+tmpVec[2]*wallVec[2];
    
      // normW=std::sqrt(tmpVec[0]*tmpVec[0]+tmpVec[1]*tmpVec[1]+tmpVec[2]*tmpVec[2]);
      // std::cerr<<std::endl<<" angle 2  "<<norm1/normW<<std::endl;
      if (normW>0)
  for(size_t k=0; k< 3; ++k)      
    q[k]=vertexData[V1W2][k]+norm1*wallVec[k];    
      else{
  std::cerr<<" in DivisionFlagResetShortestPath q is wrong"<<std::endl;
  exit(0);    
      }
    

    
    }

    // std::cerr<<" before:: "<<std::endl;
    // for(size_t k=0; k< cellData[cell.index()].size(); ++k)
    //   std::cerr<<cellData[cell.index()][k]<<"  ";
    // std::cerr<<std::endl;

    if (numVariableIndexLevel() == 2)
      {
  const size_t timeIndex = variableIndex(1, 0);
      
  const double age = cellData[cell.index()][timeIndex];
      
  std::cerr << "Cell age at division is " << age << "\n";
      
  cellData[cell.index()][timeIndex] = 0.0;
      }


  
    if(numParameter()==6 && parameter(4)==1) {// centerTriangulation
      if(parameter(5)==0 || parameter(5)==1)
  T->divideCellCenterTriangulation(&cell, winner.wall1, winner.wall2,
           variableIndex(2,0),variableIndex(2,1),
           p, q, cellData, wallData, vertexData,
           cellDerivs,wallDerivs,vertexDerivs,variableIndex(0),
           parameter(2),
           parameter(5));
      else{
  std::cerr<<"parameter(5) should be 0 or 1"
     <<std::endl;
  exit(-1);
      }
    }
    else
      { //if numParameter()==6 && parameter(4)==1 is not fulfilled
      double model0=1.0;
      if (model0==1.0)
       {

       cellData[i][8]=0;        //resetting variable to 0 when dividing
       cellData[i][9]=0.0;      //resetting variable to 0 when dividing
       cellData[i][10]=0.0;     //resetting variable to 0 when dividing

       }
      else
       {
       cellData[i][14]=0.0; //resetting variable to 0 when dividing
       }

      //std::cerr<<"dividing "<<i<<"\n"<<std::endl;
      //std::cerr<<cellData[i][5]<<"\n"<<std::endl;
      //std::cerr<<"Former size"<<cellData.size()<<i<<"\n"<<std::endl;
      //std::cerr<<cellData[i].size()<<"\n"<<std::endl;
      //int numcelli=(T -> numCell());
      //std::cerr<<"numcells "<<numcelli<<"\n"<<std::endl;

      T->divideCell(&cell, winner.wall1, winner.wall2, 
        p, q, cellData, wallData, vertexData,
        cellDerivs, wallDerivs, vertexDerivs, variableIndex(0), 
        parameter(2));

      //std::cerr<<"divided "<<i<<"\n"<<std::endl;
      //std::cerr<<cellData[i][5]<<"\n"<<std::endl;
      //std::cerr<<"Latter size"<<cellData.size()<<i<<"\n"<<std::endl;

      //std::cerr<<"divided conc "<<cellData[(T -> numCell())-1][5]<<"\n"<<std::endl;
     // std::cerr<<"X divided conc "<<cellData[(T -> numCell())-2][5]<<"\n"<<std::endl;
      int numcell=(T -> numCell());
      //std::cerr<<"numcells "<<numcell<<"\n"<<std::endl;
      
      int maxcellnum=numcell;

      for(size_t kk=0; kk< numcell; ++kk)
        {
          if (maxcellnum<cellData[kk][17])
            {maxcellnum=cellData[kk][17];}
          if (maxcellnum<cellData[kk][18])
            {maxcellnum=cellData[kk][18];}
      }

      //inheriting the mother cell index from the lineage indexing before it is updated
      cellData[i][19]=cellData[i][18];
      cellData[(T -> numCell())-1][19]=cellData[i][18];
      // lineage indexing
      cellData[i][18]=maxcellnum+1;
      cellData[(T -> numCell())-1][18]=maxcellnum+2;

    }
  
    assert (numWallTmp + 3 == T->numWall());
  
    //Change length of new wall between the divided daugther cells
    wallData[numWallTmp][0] *= parameter(1);
  

    //  std::cerr<<"  after::  "<<std::endl;
    //  for(size_t k=0; k< cellData[cell.index()].size(); ++k)
    //    std::cerr<<cellData[cell.index()][k]<<"  ";
    //  std::cerr<<std::endl;
    // std::cerr<<std::endl;
    //  for(size_t k=0; k< cellData[(T -> numCell())-1].size(); ++k)
    //    std::cerr<<cellData[(T -> numCell())-1][k]<<"  ";
    //  std::cerr<<std::endl;



    //Check that the division did not mess up the data structure
    //T->checkConnectivity(1);
  }

  std::vector<FlagResetShortestPath::Candidate> 
  FlagResetShortestPath::getCandidates(Tissue* T, size_t i,
            DataMatrix &cellData,
            DataMatrix &wallData,
            DataMatrix &vertexData,
            DataMatrix &cellDerivs,
            DataMatrix &wallDerivs,
            DataMatrix &vertexDerivs)
  {
    Cell cell = T->cell(i);
  
    assert(cell.numWall() > 1);
  
    std::vector<double> o;
  
    if (parameter(3) == 1)
      {
  o = cell.positionFromVertex(vertexData);
      }
    else
      {
  try
    {
      o = cell.randomPositionInCell(vertexData);
    }
  catch (Cell::FailedToFindRandomPositionInCellException)
    {
      return std::vector<Candidate>();
    }
      }
  
    double ox = o[0];
    double oy = o[1];
  

    std::vector<Candidate> candidates;
  
    for (size_t i = 0; i < cell.numWall() - 1; ++i) {
      for (size_t j = i + 1; j < cell.numWall(); ++j) {
  Wall *wall1 = cell.wall(i);
  Wall *wall2 = cell.wall(j);
  size_t wall1Index = i;
  size_t wall2Index = j;
      
  //      std::cerr << "i = " << wall1->index() << " : j = " << wall2->index() << std::endl;
  //        std::cerr << "o = (" << ox << ", " << oy << ")" << std::endl;
      
  double x1x;
  double x1y;
  double x2x;
  double x2y;
      
  double vx;
  double vy;
      
  double x1px;
  double x1py;
  double x2px;
  double x2py;
      
  double ux;
  double uy;
      
  bool flippedVectors;
      
  do {
    flippedVectors = false;
        
    x1x = vertexData[wall1->vertex1()->index()][0];
    x1y = vertexData[wall1->vertex1()->index()][1];
    x2x = vertexData[wall1->vertex2()->index()][0];
    x2y = vertexData[wall1->vertex2()->index()][1];
  
    vx = x2x - x1x;
    vy = x2y - x1y;
  
    if (vx * (oy - x1y) - vy * (ox - x1x) > 0) {
      //std::cerr << "Change v" << std::endl;
      double tmpx = x1x;
      double tmpy = x1y;
      x1x = x2x;
      x1y = x2y;
      x2x = tmpx;
      x2y = tmpy;
      vx = -vx;
      vy = -vy;
    }
        
        
    x1px = vertexData[wall2->vertex1()->index()][0];
    x1py = vertexData[wall2->vertex1()->index()][1];
    x2px = vertexData[wall2->vertex2()->index()][0];
    x2py = vertexData[wall2->vertex2()->index()][1];
  
    ux = x2px - x1px;
    uy = x2py - x1py;
  
    if (ux * (oy - x1py) - uy * (ox - x1px) < 0) {
      //            std::cerr << "Change u" << std::endl;
      double tmpx = x1px;
      double tmpy = x1py;
      x1px = x2px;
      x1py = x2py;
      x2px = tmpx;
      x2py = tmpy;
      ux = -ux;
      uy = -uy;
    }
        
    if (vx * uy - vy * ux > 0) {
      //            std::cerr << "Flipped walls" << std::endl;
      Wall *tmp = wall1;
      wall1 = wall2;
      wall2 = tmp;
      size_t tmpIndex = wall1Index;
      wall1Index = wall2Index;
      wall2Index = tmpIndex;
      flippedVectors = true;
    }
  } while (flippedVectors == true);
      
  double wx = ox - x1x;
  double wy = oy - x1y;
  double wpx = ox - x1px;
  double wpy = oy - x1py;
      
  double dvx = wx - ((vx * wx + vy * wy) / (vx * vx + vy * vy)) * vx;
  double dvy = wy - ((vx * wx + vy * wy) / (vx * vx + vy * vy)) * vy;
  double dux = wpx - ((ux * wpx + uy * wpy) / (ux * ux + uy * uy)) * ux;
  double duy = wpy - ((ux * wpx + uy * wpy) / (ux * ux + uy * uy)) * uy;
      
      
  //      std::cerr << " x1 = (" << x1x << ", " << x1y << ")" << std::endl;
  //      std::cerr << " x2 = (" << x2x << ", " << x2y << ")" << std::endl;
  //      std::cerr << " x1p = (" << x1px << ", " << x1py << ")" << std::endl;
  //      std::cerr << " x2p = (" << x2px << ", " << x2py << ")" << std::endl;
      
  //      std::cerr << " v = (" << vx << ", " << vy << ")" << std::endl;
  //      std::cerr << " u = (" << ux << ", " << uy << ")" << std::endl;
  //      std::cerr << " w = (" << wx << ", " << wy << ")" << std::endl;
  //      std::cerr << " wp = (" << wpx << ", " << wpy << ")" << std::endl;
  //      std::cerr << " dv = (" << dvx << ", " << dvy << ")" << std::endl;
  //      std::cerr << " du = (" << dux << ", " << duy << ")" << std::endl;
      
      
  double A = std::sqrt(dvx * dvx + dvy * dvy);
  double B = std::sqrt(dux * dux + duy * duy);
      
  double sigma = std::acos((vx * ux + vy * uy) / 
         (std::sqrt(vx * vx + vy * vy) * std::sqrt(ux * ux + uy * uy)));
      
  double alpha = astar(sigma, A, B);
  double beta = myMath::pi() + sigma - alpha;
      
  double t = (vx * wx + vy * wy) / (vx * vx + vy * vy);
  double tp = t + (1.0 / std::sqrt(vx * vx + vy * vy)) * A * 
    std::sin(alpha - 0.50 * myMath::pi()) / std::sin(alpha);
      
  double s = (ux * wpx + uy * wpy) / (ux * ux + uy * uy);
  double sp = s + (1.0 / std::sqrt(ux * ux + uy * uy)) * B * 
    std::sin(beta - 0.50 * myMath::pi()) / std::sin(beta);
      
  double px = x1x + tp * vx;
  double py = x1y + tp * vy;
      
  double qx = x1px + sp * ux;
  double qy = x1py + sp * uy;
      
  //      std::cerr << " sigma = " << sigma << std::endl
  //          << " alpha = " << alpha << std::endl
  //          << " beta = " << beta << std::endl
  //          << " A = " << A << std::endl
  //          << " B = " << B << std::endl
  //          << " px = " << px << std::endl
  //          << " py = " << py << std::endl
  //          << " qx = " << qx << std::endl
  //          << " qy = " << qy << std::endl
  //          << " tp = " << tp << std::endl
  //          << " sp = " << sp << std::endl;
      
      
  double distance = std::sqrt((qx - px) * (qx - px) + (qy - py) * (qy - py));
      
  //        std::cerr << " distance = " << distance << std::endl;
      
  if (tp <= 0.0 || tp >= 1.0 || sp <= 0.0 || sp >= 1.0) {
    //          std::cerr << "Discard" << std::endl;
    continue;
  } else {
    //          std::cerr << "Keep" << std::endl;
    Candidate candidate;
    candidate.distance = distance;
    candidate.px = px;
    candidate.py = py;
    candidate.qx = qx;
    candidate.qy = qy;
    candidate.wall1 = wall1Index;
    candidate.wall2 = wall2Index;
        
    candidates.push_back(candidate);
  }
      }
    } 
  
    return candidates;
  }

  double FlagResetShortestPath::astar(double sigma, double A, double B)
  {
    double a = 0;
    double b = myMath::pi();
    double e = b - a;
    double u = f(a, sigma, A, B);
    double v = f(b, sigma, A, B);
    double c;
  
    if (myMath::sign(u) == myMath::sign(v)) {
      return 0;
    }
    for (size_t k = 0; k < 10; ++k) {
      e = 0.5 * e;
      c = a + e;
      double w = f(c, sigma, A, B);   
      if (myMath::sign(w) != myMath::sign(u)) {
	b = c;
	v = w;
      } else {
	a = c;
	u = w;
      }
      return c;
    }
  }

  double FlagResetShortestPath::f(double a, double sigma, double A, double B)
  {
    double tmp = - A * std::cos(a) / (std::sin(a) * std::sin(a));
    tmp += B * std::cos(myMath::pi() + sigma - a) / (std::sin(sigma - a) * std::sin(sigma - a));
    return tmp;
  }


 Random::Random(std::vector<double> &paraValue, std::vector< std::vector<size_t> > &indValue)
 {
   if (paraValue.size() != 3) {
     std::cerr << "DivisionRandom::DivisionRandom() "
	       << "Three parameters are used V_threshold, Lwall_fraction, Lwall_threshold." << std::endl;
     exit(EXIT_FAILURE);
   }
   
   if (indValue.size() != 1) {
     std::cerr << "Division::Random::Random() "
	       << "First level: Variable indices for volume dependent cell "
	       << "variables are used.\n";
     exit(EXIT_FAILURE);
   }
   
   setId("Division::Random");
   setNumChange(1);
   setParameter(paraValue);  
   setVariableIndex(indValue);
   
   std::vector<std::string> tmp(numParameter());
   tmp.resize (numParameter());
   tmp[0] = "V_threshold";
   tmp[1] = "Lwall_threshold";
   setParameterId(tmp);
 }

 int Random::flag(Tissue *T, size_t i,
		  DataMatrix &cellData,
		  DataMatrix &wallData,
		  DataMatrix &vertexData,
		  DataMatrix &cellDerivs,
		  DataMatrix &wallDerivs,
		  DataMatrix &vertexDerivs)
 {
   if (T->cell(i).calculateVolume(vertexData) > parameter(0)) {
     std::cerr << "Cell " << i << " marked for division with volume " 
	       << T->cell(i).volume() << std::endl;
     return 1;
   } 
   return 0;
 }
 
 void Random::update(Tissue* T, size_t i,
		     DataMatrix &cellData,
		     DataMatrix &wallData,
		     DataMatrix &vertexData,
		     DataMatrix &cellDerivs,
		     DataMatrix &wallDerivs,
		     DataMatrix &vertexDerivs)
 {
   Cell cell = T->cell(i);
   
   assert(vertexData[0].size() == 2); // Make sure dimension == 2
   
   size_t wall1Index = random(cell.numWall());
   size_t wall2Index;
   
   while (true) {
     wall2Index = random(cell.numWall());
     
     if (wall1Index != wall2Index) {
       break;
     }
   }
   
   std::vector<double> p(2);
   
   // Find first vertex.
   {
     Wall *wall = cell.wall(wall1Index);
     
     Vertex *vertex1 = wall->vertex1();
     Vertex *vertex2 = wall->vertex2();
     
     double r = 0.0;
     
     while (r < 0.01 || r > 0.99) {
       r = myRandom::Rnd();
     }
     
     p[0] = vertexData[vertex1->index()][0] + r * (vertexData[vertex2->index()][0] - vertexData[vertex1->index()][0]);
     p[1] = vertexData[vertex1->index()][1] + r * (vertexData[vertex2->index()][1] - vertexData[vertex1->index()][1]);
   }
   
   // Find second vertex.
   std::vector<double> q(2);
   
   {
     Wall *wall = cell.wall(wall2Index);
     
     Vertex *vertex1 = wall->vertex1();
     Vertex *vertex2 = wall->vertex2();
     
     double r = 0.0;
     
     while (r < 0.01 || r > 0.99) {
       r = myRandom::Rnd();
     }
     
     q[0] = vertexData[vertex1->index()][0] + r * (vertexData[vertex2->index()][0] - vertexData[vertex1->index()][0]);
     q[1] = vertexData[vertex1->index()][1] + r * (vertexData[vertex2->index()][1] - vertexData[vertex1->index()][1]);
   }
   
   T->divideCell(&cell, wall1Index, wall2Index, p, q, cellData, wallData, vertexData, cellDerivs, wallDerivs, vertexDerivs, variableIndex(0), parameter(2));
   
   // Change length of new wall between the divided daugther cells.
   wallData[T->numWall()-1][0] *= parameter(1);
 }
 
 int Random::random(int n)
 {
   double r = myRandom::Rnd();
   
   int result = (int) floor(n * r);
   
   if (result == n) {
     return result - 1;
   } else {
     return result;
   }
 }
 
 MainAxis::MainAxis(std::vector<double> &paraValue, 
		    std::vector< std::vector<size_t> > &indValue)
 {
   //Do some checks on the parameters and variable indeces
   //////////////////////////////////////////////////////////////////////
   if (paraValue.size() != 4) {
     std::cerr << "DivisionMainAxis::DivisionMainAxis() "
		<< "Four parameters are used V_threshold, Lwall_fraction, Lwall_threshold, and direction flag (0 = perpendicular to main axis, 1 = parallel to main axis)\n";
      exit(EXIT_FAILURE);
    }
	
    if (indValue.size() != 1) {
      std::cerr << "DivisionMainAxis::DivisionMainAxis() "
		<< "First level: Variable indices for volume dependent cell variables are used.\n";
      exit(EXIT_FAILURE);
    }

    //Set the variable values
    //////////////////////////////////////////////////////////////////////
    setId("MainAxis");
    setNumChange(1);
    setParameter(paraValue);  
    setVariableIndex(indValue);
	
    //Set the parameter identities
    //////////////////////////////////////////////////////////////////////
    std::vector<std::string> tmp( numParameter() );
    tmp.resize( numParameter() );
    tmp[0] = "V_threshold";
    tmp[1] = "LWall_frac";
    tmp[2] = "LWall_threshold";
    tmp[3] = "direction flag";
    setParameterId( tmp );
  }

  int MainAxis::flag(Tissue *T, size_t i,
		     DataMatrix &cellData,
		     DataMatrix &wallData,
		     DataMatrix &vertexData,
		     DataMatrix &cellDerivs,
		     DataMatrix &wallDerivs,
		     DataMatrix &vertexDerivs)
  {
    if (T->cell(i).calculateVolume(vertexData) > parameter(0)) {
      std::cerr << "Cell " << i << " marked for division with volume " << T->cell(i).volume() << std::endl;
      return 1;
    } else { 
      return 0;
    }
  }

  void MainAxis::update(Tissue *T, size_t cellI,
			DataMatrix &cellData,
			DataMatrix &wallData,
			DataMatrix &vertexData,
			DataMatrix &cellDeriv,
			DataMatrix &wallDeriv,
			DataMatrix &vertexDeriv)
  {
	
    Cell &cell = T->cell(cellI);
    size_t dimension = vertexData[0].size();
	
    if (dimension != 2) {
      std::cerr << "DivisionMainAxis only supports two dimensions.\n";
      exit(EXIT_FAILURE);
    }
	
    std::vector<double> com(2);
	
    for (size_t i = 0, e = cell.numVertex(); i < e; ++i) {
      size_t vertexIndex = cell.vertex(i)->index();
		
      com[0] += vertexData[vertexIndex][0];
      com[1] += vertexData[vertexIndex][1];
    }
	
    com[0] /= cell.numVertex();
    com[1] /= cell.numVertex();
	
    std::vector<double> n = getMainAxis(cell, vertexData);

    if (parameter(3) == 0)
      {
	double tmp = n[0];
	n[0] = -n[1];
	n[1] = tmp;
      }

    std::vector<Candidate> candidates;

    for (size_t i = 0, e = cell.numWall(); i < e; ++i) {
      Wall *wall = cell.wall(i);
		
      double ax = vertexData[wall->vertex1()->index()][0];
      double ay = vertexData[wall->vertex1()->index()][1];
      double bx = vertexData[wall->vertex2()->index()][0];
      double by = vertexData[wall->vertex2()->index()][1];
		
      double s = ((ax - com[0]) * n[1] - (ay - com[1]) * n[0]) / ((ax - bx) * n[1] - (ay - by) * n[0]);

      if (s != s) {
	continue;
      }

      if (s >= 0.0 && s <= 1.0) {
	Candidate candidate;
			
	candidate.s = s;
	candidate.index  = i;
	candidate.p.resize(2);
	candidate.p[0] = ax + s * (bx - ax);
	candidate.p[1] = ay + s * (by - ay);
			
	candidates.push_back(candidate);
      }
    }

    if (candidates.size() < 2) {
      std::cerr << "DivisionMainAxis: Unable to find enough candidates.\n";
      exit(EXIT_FAILURE);
    }
	
    if (candidates.size() > 2) {
      std::cerr << "DivisionMainAxis: Warning, more than two candidates could be found.\n";
    }
	
    std::sort(candidates.begin(), candidates.end(), CompareCandidate());
	
    Candidate &c1 = candidates[0];
    Candidate &c2 = candidates[1];
	
    size_t numWallTmp=wallData.size();
	
    T->divideCell(&cell, c1.index, c2.index, c1.p, c2.p, cellData, wallData, vertexData, cellDeriv, wallDeriv, vertexDeriv, variableIndex(0), parameter(2));
	
    //Change length of new wall between the divided daugther cells 
    wallData[numWallTmp][0] *= parameter(1);
  }

  std::vector<double> MainAxis::getMainAxis(Cell &cell, DataMatrix &vertexData)
  {
    size_t dimensions = vertexData[0].size();
    size_t numberOfVertices = cell.numVertex();
	
    // Copy vertex data to temporary container and calculate mean values.
	
    DataMatrix vertices(numberOfVertices);
    for (size_t i = 0; i < numberOfVertices; ++i) {
      vertices[i].resize(dimensions);
    }
    std::vector<double> mean(dimensions, 0.0);
	
    for (size_t i = 0; i < numberOfVertices; ++i) {
      Vertex *v = cell.vertex(i);
      for (size_t j = 0; j < dimensions; ++j) {
	vertices[i][j] = vertexData[v->index()][j];
	mean[j] += vertexData[v->index()][j];
      }
    }
	
    for (size_t i = 0; i < dimensions; ++i) {
      mean[i] /= numberOfVertices;
    }
	
    // Subtract mean from data to get an expectation value equal to zero.
	
    for (size_t i = 0; i < numberOfVertices; ++i) {
      for (size_t j = 0; j < dimensions; ++j) {
	vertices[i][j] -= mean[j];
      }
    }
	
    // Calculate the correlation matrix.
	
    DataMatrix R(dimensions);
	
    for (size_t i = 0; i < dimensions; ++i) {
      R[i].resize(dimensions);
      for (size_t j = 0; j < dimensions; ++j) {
	R[i][j] = 0.0;
      }
    }
	
    for (size_t k = 0; k < dimensions; ++k) {
      for (size_t l = 0; l < dimensions; ++l) {
	for (size_t i = 0; i < numberOfVertices; ++i) {
	  R[k][l] += vertices[i][k] * vertices[i][l];
	}
	R[k][l] /= numberOfVertices;
      }
    }
	
    // Find the eigenvectors with the two greatests corresponding eigenvalues.
	
    DataMatrix candidates;
	
    DataMatrix V;
    std::vector<double> d;
	
    myMath::jacobiTransformation(R , V, d);
	
    size_t max = d.size();
	
    for (size_t i = 0; i < d.size(); ++i) {
      if (std::abs(d[i]) >= std::abs(d[max])) {
	max = i;
      }
    }
	
    return V[max];
  }

  VolumeRandomDirectionGiantCells::VolumeRandomDirectionGiantCells(std::vector<double> &paraValue, 
								   std::vector< std::vector<size_t> > 
								   &indValue)
  {
    //Do some checks on the parameters and variable indeces
    //////////////////////////////////////////////////////////////////////
    if (paraValue.size() != 5) {
      std::cerr << "DivisionVolumeRandomDirectionGiantCells::DivisionVolumeRandomDirectionGiantCells() "
		<< "Four parameters used V_threshold, LWall_fraction, Lwall_threshold, COM (1 = COM, 0 = Random), and giant cell factor to V_threshold.\n";
      std::exit(EXIT_FAILURE);
    }
	
    if (indValue.size() != 2 || (indValue.size() == 2 && indValue[1].size() !=1 )) {
      std::cerr << "DivisionVolumeRandomDirectionGiantCells::DivisionVolumeRandomDirectionGiantCells() "
 		<< "First level: Variable indices for volume dependent cell variables are used.\n"
 		<< "Second level: Varible index for giant cell flag.\n";
      std::exit(EXIT_FAILURE);
    }
	
    //Set the variable values
    //////////////////////////////////////////////////////////////////////
    setId("DivisionVolumeRandomDirectionGiantCells");
    setNumChange(1);
    setParameter(paraValue);  
    setVariableIndex(indValue);
	
    //Set the parameter identities
    //////////////////////////////////////////////////////////////////////
    std::vector<std::string> tmp(numParameter());
    tmp.resize(numParameter());
    tmp[0] = "V_threshold";
    tmp[1] = "LWall_frac";
    tmp[2] = "LWall_threshold";
    tmp[3] = "COM";
    tmp[4] = "giant_frac";
    setParameterId(tmp);
  }

  int VolumeRandomDirectionGiantCells::flag(Tissue *T, size_t i,
					    DataMatrix &cellData,
					    DataMatrix &wallData,
					    DataMatrix &vertexData,
					    DataMatrix &cellDerivs,
					    DataMatrix &wallDerivs,
					    DataMatrix &vertexDerivs)
  {
    if (cellData[i][variableIndex(1, 0)] && T->cell(i).calculateVolume(vertexData) > parameter(0) * parameter(4))
      {
	std::cerr << "Giant Cell " << i << " marked for division with volume " << T->cell(i).volume() << std::endl;
	return 1;
		
      }
    else if (!cellData[i][variableIndex(1, 0)] && T->cell(i).calculateVolume(vertexData) > parameter(0))
      {
	std::cerr << "Cell " << i << " marked for division with volume " << T->cell(i).volume() << std::endl;
	return 1;
      } 

    return 0;
  }

  void VolumeRandomDirectionGiantCells::
  update(Tissue *T,size_t cellI,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDeriv,
	 DataMatrix &wallDeriv,
	 DataMatrix &vertexDeriv ) {
  
    Cell *divCell = &(T->cell(cellI));
    size_t dimension = vertexData[0].size();
    //size_t numV = divCell->numVertex();
    assert( divCell->numWall() > 2 );
    assert( dimension==2 );
  
    std::vector<double> com;

    if (parameter(3) == 1)
      {
	com = divCell->positionFromVertex(vertexData);
      }
    else
      {
	try
	  {
	    com = divCell->randomPositionInCell(vertexData);
	  }
	catch (Cell::FailedToFindRandomPositionInCellException)
	  {
	    return;
	  }
      }
  
    std::vector<double> n(dimension);
    double phi=2*3.14*myRandom::Rnd();
    n[0] = std::sin(phi);
    n[1] = std::cos(phi);
  
    //Find two (and two only) intersecting walls
    //////////////////////////////////////////////////////////////////////
    std::vector<size_t> wI(2);
    std::vector<double> s(2);
    wI[0]=0;
    wI[1]=divCell->numWall();
    s[0]=s[1]=-1.0;
    //double minDist,w3s;
    std::vector<size_t> w3Tmp;
    std::vector<double> w3tTmp;
    int flag=0;
    for( size_t k=0 ; k<divCell->numWall() ; ++k ) {
      size_t v1Tmp = divCell->wall(k)->vertex1()->index();
      size_t v2Tmp = divCell->wall(k)->vertex2()->index();
      std::vector<double> w3(dimension),w0(dimension);
      for( size_t d=0 ; d<dimension ; ++d ) {
	w3[d] = vertexData[v2Tmp][d]-vertexData[v1Tmp][d];
	w0[d] = com[d]-vertexData[v1Tmp][d];
      }
      double a=0.0,b=0.0,c=0.0,d=0.0,e=0.0;//a=1.0
      for( size_t dim=0 ; dim<dimension ; ++dim ) {
	a += n[dim]*n[dim];
	b += n[dim]*w3[dim];
	c += w3[dim]*w3[dim];
	d += n[dim]*w0[dim];
	e += w3[dim]*w0[dim];
      }
      double fac=a*c-b*b;//a*c-b*b
      if( fac>0.0 ) {//else parallell and not applicable
	fac = 1.0/fac;
	//double s = fac*(b*e-c*d);
	double t = fac*(a*e-b*d);//fac*(a*e-b*d)
	if( t>=0.0 && t<1.0 ) {//within wall
	  //double dx0 = w0[0] +fac*((b*e-c*d)*nW2[0]+()*w3[0]); 					
	  w3Tmp.push_back(k);
	  w3tTmp.push_back(t);
	  std::cerr << "Dividing cell " << divCell->index() << " via wall "
		    << k << " at t=" << t << std::endl;
	  if( flag<2 ) {
	    s[flag] = t;
	    wI[flag] = k;
	  }				
	  flag++;
	}
      }
    }
    assert( wI[1] != divCell->numWall() && wI[0] != wI[1] );
    if( flag != 2 ) {
      return;
      // std::cerr << "divideVolumeVisStrain::update Warning"
      // 	      << " not two walls possible as connection "
      // 	      << "for cell " 
      // 	      << cellI << std::endl; 
      // for( size_t k=0 ; k<divCell->numWall() ; ++k ) {
      //   std::cerr << "0 " 
      // 		<< vertexData[divCell->wall(k)->vertex1()->index()][0]
      // 		<< " " 
      // 		<< vertexData[divCell->wall(k)->vertex1()->index()][1]
      // 		<< "\n0 " 
      // 		<< vertexData[divCell->wall(k)->vertex2()->index()][0]
      // 		<< " " 
      // 		<< vertexData[divCell->wall(k)->vertex2()->index()][1]
      // 		<< "\n\n\n";
      // }
      // for( size_t kk=0 ; kk<w3Tmp.size() ; ++kk ) {
      //   size_t k = w3Tmp[kk];
      //   std::cerr << "1 " 
      // 		<< vertexData[divCell->wall(k)->vertex1()->index()][0]
      // 		<< " " 
      // 		<< vertexData[divCell->wall(k)->vertex1()->index()][1]
      // 		<< "\n1 " 
      // 		<< vertexData[divCell->wall(k)->vertex2()->index()][0]
      // 		<< " " 
      // 		<< vertexData[divCell->wall(k)->vertex2()->index()][1]
      // 		<< "\n\n\n";
      // }
      // std::cerr << "2 " 
      // 	      << vertexData[divCell->wall(wI[0])->vertex1()->index()][0]
      // 	      << " " 
      // 	      << vertexData[divCell->wall(wI[0])->vertex1()->index()][1]
      // 	      << "\n2 " 
      // 	      << vertexData[divCell->wall(wI[0])->vertex2()->index()][0]
      // 	      << " " 
      // 	      << vertexData[divCell->wall(wI[0])->vertex2()->index()][1]
      // 	      << "\n\n\n";
      // std::cerr << "3 " 
      // 	      << vertexData[divCell->wall(wI[1])->vertex1()->index()][0]
      // 	      << " " 
      // 	      << vertexData[divCell->wall(wI[1])->vertex1()->index()][1]
      // 	      << "\n3 " 
      // 	      << vertexData[divCell->wall(wI[1])->vertex2()->index()][0]
      // 	      << " " 
      // 	      << vertexData[divCell->wall(wI[1])->vertex2()->index()][1]
      // 	      << "\n\n\n";
      // std::cerr << "4 " 
      // 	      << 0.5*(vertexData[divCell->wall(wI[0])->vertex1()->index()][0]+
      // 		      vertexData[divCell->wall(wI[0])->vertex2()->index()][0])
      // 	      << " " 
      // 	      << 0.5*(vertexData[divCell->wall(wI[0])->vertex1()->index()][1]+
      // 		      vertexData[divCell->wall(wI[0])->vertex2()->index()][1])
      // 	      << "\n4 "
      // 	      << 0.5*(vertexData[divCell->wall(wI[1])->vertex1()->index()][0]+
      // 		      vertexData[divCell->wall(wI[1])->vertex2()->index()][0])
      // 	      << " " 
      // 	      << 0.5*(vertexData[divCell->wall(wI[1])->vertex1()->index()][1]+
      // 		      vertexData[divCell->wall(wI[1])->vertex2()->index()][1])
      // 	      << "\n\n\n";
      // exit(-1);
    }	
    //Addition of new vertices at walls at position 's' 
    std::vector<double> v1Pos(dimension),v2Pos(dimension);
    size_t v1I = divCell->wall(wI[0])->vertex1()->index();
    size_t v2I = divCell->wall(wI[0])->vertex2()->index();
    for( size_t d=0 ; d<dimension ; ++d )
      v1Pos[d] = vertexData[v1I][d]+ s[0]*(vertexData[v2I][d]-vertexData[v1I][d]);
    v1I = divCell->wall(wI[1])->vertex1()->index();
    v2I = divCell->wall(wI[1])->vertex2()->index();
    for( size_t d=0 ; d<dimension ; ++d )
      v2Pos[d] = vertexData[v1I][d]+s[1]*(vertexData[v2I][d]-vertexData[v1I][d]);
  
    //Add one cell, three walls, and two vertices
    //////////////////////////////////////////////////////////////////////
    //Save number of walls
    size_t numWallTmp=wallData.size();
    assert( numWallTmp==T->numWall() );
    //Divide
    T->divideCell(divCell,wI[0],wI[1],v1Pos,v2Pos,cellData,wallData,vertexData,
		  cellDeriv,wallDeriv,vertexDeriv,variableIndex(0),
		  parameter(2));

    const size_t daughterIndex = T->numCell() - 1;
	
    if (myRandom::Rnd() < 0.5)
      {
	cellData[daughterIndex][variableIndex(1, 0)] = 1;
      } 
    else
      {
	cellData[daughterIndex][variableIndex(1, 0)] = 0;
      }

    assert( numWallTmp+3 == T->numWall() );
  
    //Change length of new wall between the divided daugther cells 
    wallData[numWallTmp][0] *= parameter(1);
  
    //Check that the division did not mess up the data structure
    //T->checkConnectivity(1);		
  }



  ShortestPathGiantCells::ShortestPathGiantCells(std::vector<double> &paraValue, 
						 std::vector< std::vector<size_t> > &indValue)
  {
    if (paraValue.size() != 5) {
      std::cerr << "DivisionShortestPathGiantCells::DivisionShortestPathGiantCells() "
		<< "Five parameters are used V_threshold, Lwall_fraction, Lwall_threshold, COM (1 = COM, 0 = Random), and giant cell factor to V_threshold.\n";
      std::exit(EXIT_FAILURE);
    }
	
    if (indValue.size() != 2 || (indValue.size() == 2 && indValue[1].size() != 1)) {
      std::cerr << "DivisionShortestPathGiantCells::DivisionShortestPathGiantCells() "
		<< "First level: Variable indices for volume dependent cell variables are used.\n"
		<< "Second level: Varible index for giant cell flag.\n";
      std::exit(EXIT_FAILURE);
    }

    setId("DivisionShortestPathGiantCells");
    setNumChange(1);
    setParameter(paraValue);  
    setVariableIndex(indValue);
  
    std::vector<std::string> tmp(numParameter());
    tmp.resize (numParameter());
    tmp[0] = "V_threshold";
    tmp[1] = "Lwall_fraction";
    tmp[2] = "Lwall_threshold";
    tmp[3] = "COM";
    tmp[4] = "giantcell_factor";
    setParameterId(tmp);
  }

  int ShortestPathGiantCells::flag(Tissue *T, size_t i,
				   DataMatrix &cellData,
				   DataMatrix &wallData,
				   DataMatrix &vertexData,
				   DataMatrix &cellDerivs,
				   DataMatrix &wallDerivs,
				   DataMatrix &vertexDerivs)
  {
    if (cellData[i][variableIndex(1, 0)] && T->cell(i).calculateVolume(vertexData) > parameter(0) * parameter(4))
      {
	std::cerr << "Giant Cell " << i << " marked for division with volume " << T->cell(i).volume() << std::endl;
	return 1;
      }
    else if (!cellData[i][variableIndex(1, 0)] && T->cell(i).calculateVolume(vertexData) > parameter(0))
      {
	std::cerr << "Cell " << i << " marked for division with volume " << T->cell(i).volume() << std::endl;
	return 1;
      }

    return 0;
  }

  void ShortestPathGiantCells::update(Tissue* T, size_t i,
				      DataMatrix &cellData,
				      DataMatrix &wallData,
				      DataMatrix &vertexData,
				      DataMatrix &cellDerivs,
				      DataMatrix &wallDerivs,
				      DataMatrix &vertexDerivs)
  {
    Cell &cell = T->cell(i);
	
    if (vertexData[0].size() != 2)
      {
	std::cerr << "DivisionShortestPathGiantCells only supports two dimensions.\n";
	std::exit(EXIT_FAILURE);
      }

    std::vector<Candidate> candidates = getCandidates(T, i, cellData, wallData, vertexData, cellDerivs, wallDerivs, vertexDerivs);
	
    if (candidates.size() == 0)
      {
	return;
      }
  
    Candidate winner = { std::numeric_limits<double>::max(), 0, 0, 0, 0, 0, 0 };
    for (size_t i = 0; i < candidates.size(); ++i) {
      if (candidates[i].distance < winner.distance) {
	winner = candidates[i];
      }
    }
	
    // 	std::cerr << "Winner: " << std::endl
    // 		  << " distance = " << winner.distance << std::endl
    // 		  << " p = (" << winner.px << ", " << winner.py << ")" << std::endl
    // 		  << " q = (" << winner.qx << ", " << winner.qy << ")" << std::endl;
	
    size_t numWallTmp = wallData.size();
    assert(numWallTmp == T->numWall());
	
    std::vector<double> p(2);
    p[0] = winner.px;
    p[1] = winner.py;
    std::vector<double> q(2);
    q[0] = winner.qx;
    q[1] = winner.qy;
	
    T->divideCell(&cell, winner.wall1, winner.wall2, p, q, cellData, wallData, vertexData,
		  cellDerivs, wallDerivs, vertexDerivs, variableIndex(0), parameter(2));

    const size_t daughterIndex = T->numCell() - 1;
	
    if (myRandom::Rnd() < 0.5) 
      {
	cellData[daughterIndex][variableIndex(1, 0)] = 1;
      }
    else 
      {
	cellData[daughterIndex][variableIndex(1, 0)] = 0;
      }
	
    assert (numWallTmp + 3 == T->numWall());
	
    //Change length of new wall between the divided daugther cells
    wallData[numWallTmp][0] *= parameter(1);
	
    //Check that the division did not mess up the data structure
    //T->checkConnectivity(1);
  }

  std::vector<ShortestPathGiantCells::Candidate> 
  ShortestPathGiantCells::getCandidates(Tissue* T, size_t i,
					DataMatrix &cellData,
					DataMatrix &wallData,
					DataMatrix &vertexData,
					DataMatrix &cellDerivs,
					DataMatrix &wallDerivs,
					DataMatrix &vertexDerivs)
  {
    Cell cell = T->cell(i);
  
    assert(cell.numWall() > 1);
  
    std::vector<double> o;
  
    if (parameter(3) == 1)
      {
	o = cell.positionFromVertex(vertexData);
      }
    else
      {
	try
	  {
	    o = cell.randomPositionInCell(vertexData);
	  }
	catch (Cell::FailedToFindRandomPositionInCellException)
	  {
	    return std::vector<Candidate>();
	  }
      }
  
    double ox = o[0];
    double oy = o[1];
  
    std::vector<Candidate> candidates;
  
    for (size_t i = 0; i < cell.numWall() - 1; ++i) {
      for (size_t j = i + 1; j < cell.numWall(); ++j) {
	Wall *wall1 = cell.wall(i);
	Wall *wall2 = cell.wall(j);
	size_t wall1Index = i;
	size_t wall2Index = j;
      
	// 			std::cerr << "i = " << wall1->index() << " : j = " << wall2->index() << std::endl;
	//  			std::cerr << "o = (" << ox << ", " << oy << ")" << std::endl;
      
	double x1x;
	double x1y;
	double x2x;
	double x2y;
      
	double vx;
	double vy;
      
	double x1px;
	double x1py;
	double x2px;
	double x2py;
      
	double ux;
	double uy;
      
	bool flippedVectors;
      
	do {
	  flippedVectors = false;
	
	  x1x = vertexData[wall1->vertex1()->index()][0];
	  x1y = vertexData[wall1->vertex1()->index()][1];
	  x2x = vertexData[wall1->vertex2()->index()][0];
	  x2y = vertexData[wall1->vertex2()->index()][1];
	
	  vx = x2x - x1x;
	  vy = x2y - x1y;
	
	  if (vx * (oy - x1y) - vy * (ox - x1x) > 0) {
	    //  					std::cerr << "Change v" << std::endl;
	    double tmpx = x1x;
	    double tmpy = x1y;
	    x1x = x2x;
	    x1y = x2y;
	    x2x = tmpx;
	    x2y = tmpy;
	    vx = -vx;
	    vy = -vy;
	  }
	
	
	  x1px = vertexData[wall2->vertex1()->index()][0];
	  x1py = vertexData[wall2->vertex1()->index()][1];
	  x2px = vertexData[wall2->vertex2()->index()][0];
	  x2py = vertexData[wall2->vertex2()->index()][1];
	
	  ux = x2px - x1px;
	  uy = x2py - x1py;
	
	  if (ux * (oy - x1py) - uy * (ox - x1px) < 0) {
	    //  					std::cerr << "Change u" << std::endl;
	    double tmpx = x1px;
	    double tmpy = x1py;
	    x1px = x2px;
	    x1py = x2py;
	    x2px = tmpx;
	    x2py = tmpy;
	    ux = -ux;
	    uy = -uy;
	  }
	
	  if (vx * uy - vy * ux > 0) {
	    //  					std::cerr << "Flipped walls" << std::endl;
	    Wall *tmp = wall1;
	    wall1 = wall2;
	    wall2 = tmp;
	    size_t tmpIndex = wall1Index;
	    wall1Index = wall2Index;
	    wall2Index = tmpIndex;
	    flippedVectors = true;
	  }
	} while (flippedVectors == true);
      
	double wx = ox - x1x;
	double wy = oy - x1y;
	double wpx = ox - x1px;
	double wpy = oy - x1py;
      
	double dvx = wx - ((vx * wx + vy * wy) / (vx * vx + vy * vy)) * vx;
	double dvy = wy - ((vx * wx + vy * wy) / (vx * vx + vy * vy)) * vy;
	double dux = wpx - ((ux * wpx + uy * wpy) / (ux * ux + uy * uy)) * ux;
	double duy = wpy - ((ux * wpx + uy * wpy) / (ux * ux + uy * uy)) * uy;
      
      
	// 			std::cerr << " x1 = (" << x1x << ", " << x1y << ")" << std::endl;
	// 			std::cerr << " x2 = (" << x2x << ", " << x2y << ")" << std::endl;
	// 			std::cerr << " x1p = (" << x1px << ", " << x1py << ")" << std::endl;
	// 			std::cerr << " x2p = (" << x2px << ", " << x2py << ")" << std::endl;
      
	// 			std::cerr << " v = (" << vx << ", " << vy << ")" << std::endl;
	// 			std::cerr << " u = (" << ux << ", " << uy << ")" << std::endl;
	// 			std::cerr << " w = (" << wx << ", " << wy << ")" << std::endl;
	// 			std::cerr << " wp = (" << wpx << ", " << wpy << ")" << std::endl;
	// 			std::cerr << " dv = (" << dvx << ", " << dvy << ")" << std::endl;
	// 			std::cerr << " du = (" << dux << ", " << duy << ")" << std::endl;
      
      
	double A = std::sqrt(dvx * dvx + dvy * dvy);
	double B = std::sqrt(dux * dux + duy * duy);
      
	double sigma = std::acos((vx * ux + vy * uy) / (std::sqrt(vx * vx + vy * vy) * std::sqrt(ux * ux + uy * uy)));
      
	double alpha = astar(sigma, A, B);
	double beta = myMath::pi() + sigma - alpha;
      
	double t = (vx * wx + vy * wy) / (vx * vx + vy * vy);
	double tp = t + (1.0 / std::sqrt(vx * vx + vy * vy)) * A * std::sin(alpha - 0.50 * myMath::pi()) / std::sin(alpha);
      
	double s = (ux * wpx + uy * wpy) / (ux * ux + uy * uy);
	double sp = s + (1.0 / std::sqrt(ux * ux + uy * uy)) * B * std::sin(beta - 0.50 * myMath::pi()) / std::sin(beta);
      
	double px = x1x + tp * vx;
	double py = x1y + tp * vy;
      
	double qx = x1px + sp * ux;
	double qy = x1py + sp * uy;
      
	// 			std::cerr << " sigma = " << sigma << std::endl
	// 				  << " alpha = " << alpha << std::endl
	// 				  << " beta = " << beta << std::endl
	// 				  << " A = " << A << std::endl
	// 				  << " B = " << B << std::endl
	// 				  << " px = " << px << std::endl
	// 				  << " py = " << py << std::endl
	// 				  << " qx = " << qx << std::endl
	// 				  << " qy = " << qy << std::endl
	// 				  << " tp = " << tp << std::endl
	// 				  << " sp = " << sp << std::endl;
      
      
	double distance = std::sqrt((qx - px) * (qx - px) + (qy - py) * (qy - py));
      
	//   			std::cerr << " distance = " << distance << std::endl;
      
	if (tp <= 0.0 || tp >= 1.0 || sp <= 0.0 || sp >= 1.0) {
	  //   				std::cerr << "Discard" << std::endl;
	  continue;
	} else {
	  //  				std::cerr << "Keep" << std::endl;
	  Candidate candidate;
	  candidate.distance = distance;
	  candidate.px = px;
	  candidate.py = py;
	  candidate.qx = qx;
	  candidate.qy = qy;
	  candidate.wall1 = wall1Index;
	  candidate.wall2 = wall2Index;
	
	  candidates.push_back(candidate);
	}
      }
    }	
  
    return candidates;
  }

  double ShortestPathGiantCells::astar(double sigma, double A, double B)
  {
    double a = 0;
    double b = myMath::pi();
    double e = b - a;
    double u = f(a, sigma, A, B);
    double v = f(b, sigma, A, B);
    double c;

    if (myMath::sign(u) == myMath::sign(v)) {
      return 0;
    }

    for (size_t k = 0; k < 10; ++k) {
      e = 0.5 * e;
      c = a + e;
      double w = f(c, sigma, A, B);

      if (myMath::sign(w) != myMath::sign(u)) {
	b = c;
	v = w;
      } else {
	a = c;
	u = w;
      }
    }
    return c;
  }

  double ShortestPathGiantCells::f(double a, double sigma, double A, double B)
  {
    double tmp = - A * std::cos(a) / (std::sin(a) * std::sin(a));
    tmp += B * std::cos(myMath::pi() + sigma - a) / (std::sin(sigma - a) * std::sin(sigma - a));
    return tmp;
  }

  FlagResetViaLongestWall::
  FlagResetViaLongestWall(std::vector<double> &paraValue, 
           std::vector< std::vector<size_t> > 
           &indValue ) {
    //
    // Do some checks on the parameters and variable indeces
    //
    if( paraValue.size()!=3 ) {
      std::cerr << "DivisionFlagResetViaLongestWall::"
    << "DivisionFlagResetViaLongestWall() "
    << "Three parameters used V_threshold, LWall_frac, and "
    << "Lwall_threshold." << std::endl;
      exit(EXIT_FAILURE);
  }
    if( indValue.size() != 1 ) {
      std::cerr << "DivisionFlagResetViaLongestWall::"
    << "DivisionFlagResetViaLongestWall() "
    << "Variable indices for volume dependent cell "
    << "variables is used.\n";
      exit(EXIT_FAILURE);
    }
    //
    // Set the variable values
    //
    setId("Division::FlagResetViaLongestWall");
    setNumChange(1);
    setParameter(paraValue);  
    setVariableIndex(indValue);
    //
    // Set the parameter identities
    //
    std::vector<std::string> tmp( numParameter() );
    tmp.resize( numParameter() );
    tmp[0] = "V_threshold";
    tmp[1] = "LWall_frac";
    tmp[2] = "LWall_threshold";
    setParameterId( tmp );
  }

  int FlagResetViaLongestWall::
  flag(Tissue *T,size_t i,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) {
    
    //if( cellData[i][11]==1 && cellData[i][10]>5 && cellData[i][7]==0) {
   //if(cellData[i][10]>5 && cellData[i][7]==0) {
   //if(cellData[i][11]==1 && cellData[i][7]==0) {
  if(cellData[i][variableIndex(0,0)]==1) {
      std::cerr << "Cell " << i << " marked for division" << std::endl;
      return 1;
    } 
      return 0;
  }
  
  void FlagResetViaLongestWall::
  update(Tissue *T,size_t i,
   DataMatrix &cellData,
   DataMatrix &wallData,
       DataMatrix &vertexData,
   DataMatrix &cellDeriv,
   DataMatrix &wallDeriv,
   DataMatrix &vertexDeriv ) {
    
    Cell *divCell = &(T->cell(i));
    size_t dimension = vertexData[0].size();
    assert( divCell->numWall() > 1 );
    assert( dimension==2 || dimension==3 ); 

    cellData[i][8]=0;     //resetting variable to 0 when dividing
    cellData[i][9]=0.0;   //resetting variable to 0 when dividing
    cellData[i][10]=0.0;  //resetting variable to 0 when dividing

    
    //
    // Find longest wall
    // 
    size_t wI=0,w3I=divCell->numWall();
    double maxLength = divCell->wall(0)->setLengthFromVertexPosition(vertexData);
    for( size_t k=1 ; k<divCell->numWall() ; ++k ) {
      double tmpLength = divCell->wall(k)->setLengthFromVertexPosition(vertexData);
      if( tmpLength > maxLength ) {
  wI=k;
  maxLength = tmpLength;
      }
    }   
    
    //
    // Find position for first new vertex
    //
    std::vector<double> nW(dimension),nW2(dimension),v1Pos(dimension),
      v2Pos(dimension);
    size_t v1wI = divCell->wall(wI)->vertex1()->index();
    size_t v2wI = divCell->wall(wI)->vertex2()->index();
    for( size_t d=0 ; d<dimension ; ++d ) {
      nW[d] = (vertexData[v1wI][d]-vertexData[v2wI][d])/maxLength;
      v1Pos[d] = 0.5*(vertexData[v1wI][d]+vertexData[v2wI][d]);
    }
    //
    // Find intersection with another wall via vector perpendicular to first wall
    //
    if (dimension==2) {
      nW2[1] = nW[0];
      nW2[0] = -nW[1];
    }
    else if (dimension==3) {
      nW2[0]=nW[0];
      nW2[1]=nW[1];
      nW2[2]=nW[2];
    }
    if (findSecondDivisionWall(vertexData, divCell, wI, w3I, v1Pos, nW2, v2Pos)) {
      std::cerr << "Division::VolumeViaLongestWall::update "
    << "failed to find the second wall for division!" << std::endl;
      exit(EXIT_FAILURE);
    }
    //
    // Do the division (add one cell, three walls, and two vertices)
    //
    size_t numWallTmp=wallData.size();
    assert( numWallTmp==T->numWall() );
    //Divide
    T->divideCell(divCell,wI,w3I,v1Pos,v2Pos,cellData,wallData,vertexData,
      cellDeriv,wallDeriv,vertexDeriv,variableIndex(0),
      parameter(2));
    assert( numWallTmp+3 == T->numWall() );
    
    //Change length of new wall between the divided daugther cells 
    wallData[numWallTmp][0] *= parameter(1);
    
    //Check that the division did not mess up the data structure
    //T->checkConnectivity(1);  
  }

} // end namespace Division
