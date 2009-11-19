/**
 * Filename     : compartmentDivision.cc
 * Description  : Classes describing compartmentDivision updates
 * Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
 * Created      : July 2006
 * Revision     : $Id:$
 */
#include<limits>

#include "baseCompartmentChange.h"
#include "compartmentDivision.h"
#include "myRandom.h"
#include "myMath.h"

DivisionVolumeViaLongestWall::
DivisionVolumeViaLongestWall(std::vector<double> &paraValue, 
														 std::vector< std::vector<size_t> > 
														 &indValue ) {
  //
  // Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=3 ) {
    std::cerr << "DivisionVolumeViaLongestWall::"
							<< "DivisionVolumeViaLongestWall() "
							<< "Three parameters used V_threshold, LWall_frac, and"
							<< "Lwall_threshold." << std::endl;
    exit(0);
  }
  if( indValue.size() != 1 ) {
    std::cerr << "DivisionVolumeViaLongestWall::"
							<< "DivisionVolumeViaLongestWall() "
							<< "Variable indices for volume dependent cell "
							<< "variables is used.\n";
    exit(0);
  }
  //
	// Set the variable values
  //
  setId("DivisionVolumeViaLongestWall");
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

int DivisionVolumeViaLongestWall::
flag(Tissue *T,size_t i,
     std::vector< std::vector<double> > &cellData,
     std::vector< std::vector<double> > &wallData,
     std::vector< std::vector<double> > &vertexData,
     std::vector< std::vector<double> > &cellDerivs,
     std::vector< std::vector<double> > &wallDerivs,
     std::vector< std::vector<double> > &vertexDerivs ) {
	
  if( T->cell(i).calculateVolume(vertexData) > parameter(0) ) {
    std::cerr << "Cell " << i << " marked for division with volume " 
							<< T->cell(i).volume() << std::endl;
    return 1;
  } 
  return 0;
}

void DivisionVolumeViaLongestWall::
update(Tissue *T,size_t i,
       std::vector< std::vector<double> > &cellData,
       std::vector< std::vector<double> > &wallData,
       std::vector< std::vector<double> > &vertexData,
       std::vector< std::vector<double> > &cellDeriv,
       std::vector< std::vector<double> > &wallDeriv,
       std::vector< std::vector<double> > &vertexDeriv ) {
  
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
		std::cerr << "DivisionVolumeViaLongestWall::update "
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
								cellDeriv,wallDeriv,vertexDeriv,variableIndex(0),
								parameter(2));
	assert( numWallTmp+3 == T->numWall() );
	
	//Change length of new wall between the divided daugther cells 
	wallData[numWallTmp][0] *= parameter(1);
	
	//Check that the division did not mess up the data structure
	//T->checkConnectivity(1);	
}

DivisionVolumeViaLongestWallSpatial::
DivisionVolumeViaLongestWallSpatial(std::vector<double> &paraValue, 
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

int DivisionVolumeViaLongestWallSpatial::
flag(Tissue *T,size_t i,
     std::vector< std::vector<double> > &cellData,
     std::vector< std::vector<double> > &wallData,
     std::vector< std::vector<double> > &vertexData,
     std::vector< std::vector<double> > &cellDerivs,
     std::vector< std::vector<double> > &wallDerivs,
     std::vector< std::vector<double> > &vertexDerivs ) {
	
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

void DivisionVolumeViaLongestWallSpatial::
update(Tissue *T,size_t i,
       std::vector< std::vector<double> > &cellData,
       std::vector< std::vector<double> > &wallData,
       std::vector< std::vector<double> > &vertexData,
       std::vector< std::vector<double> > &cellDeriv,
       std::vector< std::vector<double> > &wallDeriv,
       std::vector< std::vector<double> > &vertexDeriv ) {
  
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

DivisionVolumeViaLongestWall3D::
DivisionVolumeViaLongestWall3D(std::vector<double> &paraValue, 
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
int DivisionVolumeViaLongestWall3D::
flag(Tissue *T,size_t i,
     std::vector< std::vector<double> > &cellData,
     std::vector< std::vector<double> > &wallData,
     std::vector< std::vector<double> > &vertexData,
     std::vector< std::vector<double> > &cellDerivs,
     std::vector< std::vector<double> > &wallDerivs,
     std::vector< std::vector<double> > &vertexDerivs ) {
	
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
void DivisionVolumeViaLongestWall3D::
update(Tissue *T,size_t i,
       std::vector< std::vector<double> > &cellData,
       std::vector< std::vector<double> > &wallData,
       std::vector< std::vector<double> > &vertexData,
       std::vector< std::vector<double> > &cellDeriv,
       std::vector< std::vector<double> > &wallDeriv,
       std::vector< std::vector<double> > &vertexDeriv ) {
  
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

DivisionVolumeViaLongestWall3DSpatial::
DivisionVolumeViaLongestWall3DSpatial(std::vector<double> &paraValue, 
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

int DivisionVolumeViaLongestWall3DSpatial::
flag(Tissue *T,size_t i,
     std::vector< std::vector<double> > &cellData,
     std::vector< std::vector<double> > &wallData,
     std::vector< std::vector<double> > &vertexData,
     std::vector< std::vector<double> > &cellDerivs,
     std::vector< std::vector<double> > &wallDerivs,
     std::vector< std::vector<double> > &vertexDerivs ) {
	
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

void DivisionVolumeViaLongestWall3DSpatial::
update(Tissue *T,size_t i,
       std::vector< std::vector<double> > &cellData,
       std::vector< std::vector<double> > &wallData,
       std::vector< std::vector<double> > &vertexData,
       std::vector< std::vector<double> > &cellDeriv,
       std::vector< std::vector<double> > &wallDeriv,
       std::vector< std::vector<double> > &vertexDeriv ) {
  
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
DivisionVolumeViaStrain::
DivisionVolumeViaStrain(std::vector<double> &paraValue, 
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
int DivisionVolumeViaStrain::
flag(Tissue *T,size_t i,
     std::vector< std::vector<double> > &cellData,
     std::vector< std::vector<double> > &wallData,
     std::vector< std::vector<double> > &vertexData,
     std::vector< std::vector<double> > &cellDerivs,
     std::vector< std::vector<double> > &wallDerivs,
     std::vector< std::vector<double> > &vertexDerivs ) {
	
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
void DivisionVolumeViaStrain::
update(Tissue *T,size_t cellI,
       std::vector< std::vector<double> > &cellData,
       std::vector< std::vector<double> > &wallData,
       std::vector< std::vector<double> > &vertexData,
       std::vector< std::vector<double> > &cellDeriv,
       std::vector< std::vector<double> > &wallDeriv,
       std::vector< std::vector<double> > &vertexDeriv ) {
  
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
	std::vector< std::vector<double> > x(numV),y(numV),dx(numV),
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
DivisionVolumeViaDirection::
DivisionVolumeViaDirection(std::vector<double> &paraValue, 
													 std::vector< std::vector<size_t> > 
													 &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //////////////////////////////////////////////////////////////////////
  if( paraValue.size()!=4 ) {
    std::cerr << "DivisionVolumeViaDirection::"
							<< "DivisionVolumeViaDirection() "
							<< "Four parameters used V_threshold, LWall_frac, "
							<< "Lwall_threshold and Parallell_flag" << std::endl;
    exit(0);
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
  setParameterId( tmp );
}

//! Flags a cell for division if the volume above threshold
/*! 
 */
int DivisionVolumeViaDirection::
flag(Tissue *T,size_t i,
     std::vector< std::vector<double> > &cellData,
     std::vector< std::vector<double> > &wallData,
     std::vector< std::vector<double> > &vertexData,
     std::vector< std::vector<double> > &cellDerivs,
     std::vector< std::vector<double> > &wallDerivs,
     std::vector< std::vector<double> > &vertexDerivs ) {
	
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
void DivisionVolumeViaDirection::
update(Tissue *T,size_t cellI,
       std::vector< std::vector<double> > &cellData,
       std::vector< std::vector<double> > &wallData,
       std::vector< std::vector<double> > &vertexData,
       std::vector< std::vector<double> > &cellDeriv,
       std::vector< std::vector<double> > &wallDeriv,
       std::vector< std::vector<double> > &vertexDeriv ) {
  
  Cell *divCell = &(T->cell(cellI));
  size_t dimension = vertexData[0].size();
  assert( divCell->numWall() > 2 );
	assert( dimension==2 || dimension==3);
	
	// Get center of mass of the dividing cell
	std::vector<double> com(dimension);
	com = divCell->positionFromVertex(vertexData);
	
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
		std::cerr << "DivisionVolumeViaDirection::update "
							<< "failed to find two walls for division!" << std::endl;
		exit(-1);
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

//!Constructor
DivisionVolumeRandomDirection::
DivisionVolumeRandomDirection(std::vector<double> &paraValue, 
			     std::vector< std::vector<size_t> > 
			     &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //////////////////////////////////////////////////////////////////////
	if ( paraValue.size() != 4) {
		std::cerr << "DivisionVolumeRandomDirection::"
		<< "DivisionVolumeRandomDirection() "
		<< "Four parameters used V_threshold, LWall_frac, Lwall_threshold, and COM (1 = COM, 0 = Random).\n";
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
  //////////////////////////////////////////////////////////////////////
  setId("DivisionVolumeRandomDirection");
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
  tmp[3] = "COM";
  setParameterId( tmp );
}

//! Flags a cell for division if the volume above threshold
/*! 
 */
int DivisionVolumeRandomDirection::
flag(Tissue *T,size_t i,
     std::vector< std::vector<double> > &cellData,
     std::vector< std::vector<double> > &wallData,
     std::vector< std::vector<double> > &vertexData,
     std::vector< std::vector<double> > &cellDerivs,
     std::vector< std::vector<double> > &wallDerivs,
     std::vector< std::vector<double> > &vertexDerivs ) {
	
  if( T->cell(i).calculateVolume(vertexData) > parameter(0) ) {
    std::cerr << "Cell " << i << " marked for division with volume " 
							<< T->cell(i).volume() << std::endl;
    return 1;
  } 
  return 0;
}

void DivisionVolumeRandomDirection::
update(Tissue *T,size_t cellI,
       std::vector< std::vector<double> > &cellData,
       std::vector< std::vector<double> > &wallData,
       std::vector< std::vector<double> > &vertexData,
       std::vector< std::vector<double> > &cellDeriv,
       std::vector< std::vector<double> > &wallDeriv,
       std::vector< std::vector<double> > &vertexDeriv ) {
  
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
	  com = divCell->randomPositionInCell(vertexData);
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

DivisionForceDirection::DivisionForceDirection(std::vector<double> &paraValue,
									  std::vector< std::vector<size_t> > &indValue)
{
	if (paraValue.size() != 4) {
		std::cerr << "DivisionForceDirection::DivisionForceDirection() "
				<< "Four parameters are used V_threshold, Lwall_fraction, Lwall_threshold and orientation_flag." << std::endl;
		exit(EXIT_FAILURE);
	}
	
	if (indValue.size() != 2 || indValue[0].size() != 1) {
		std::cerr << "DivisionForceDirection::DivisionForceDirection() "
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

int DivisionForceDirection::flag(Tissue *T, size_t i,
						   std::vector< std::vector<double> > &cellData,
						   std::vector< std::vector<double> > &wallData,
						   std::vector< std::vector<double> > &vertexData,
						   std::vector< std::vector<double> > &cellDerivs,
						   std::vector< std::vector<double> > &wallDerivs,
						   std::vector< std::vector<double> > &vertexDerivs)
{
	if (T->cell(i).calculateVolume(vertexData) > parameter(0)) {
		std::cerr << "Cell " << i << " marked for division with volume " 
				<< T->cell(i).volume() << std::endl;
		return 1;
	} 
	return 0;
}

void DivisionForceDirection::update(Tissue* T, size_t i,
							 std::vector< std::vector<double> > &cellData,
							 std::vector< std::vector<double> > &wallData,
							 std::vector< std::vector<double> > &vertexData,
							 std::vector< std::vector<double> > &cellDerivs,
							 std::vector< std::vector<double> > &wallDerivs,
							 std::vector< std::vector<double> > &vertexDerivs)
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
	std::vector< std::vector<double> > verticesPosition;

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
DivisionVolumeViaShortestPath::
DivisionVolumeViaShortestPath(std::vector<double> &paraValue, 
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

//! Flags a cell for division if the volume above threshold
/*! 
 */
int DivisionVolumeViaShortestPath::
flag(Tissue *T,size_t i,
     std::vector< std::vector<double> > &cellData,
     std::vector< std::vector<double> > &wallData,
     std::vector< std::vector<double> > &vertexData,
     std::vector< std::vector<double> > &cellDerivs,
     std::vector< std::vector<double> > &wallDerivs,
     std::vector< std::vector<double> > &vertexDerivs ) 
{	
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
void DivisionVolumeViaShortestPath::
update(Tissue *T,size_t i,
       std::vector< std::vector<double> > &cellData,
       std::vector< std::vector<double> > &wallData,
       std::vector< std::vector<double> > &vertexData,
       std::vector< std::vector<double> > &cellDeriv,
       std::vector< std::vector<double> > &wallDeriv,
       std::vector< std::vector<double> > &vertexDeriv ) 
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
				double root = std::sqrt((A0+Ahalf-A2)/(std::cos(alpha)*std::sin(alpha)));
				double t1A = std::sqrt((x1[0]-x0[0])*(x1[0]-x0[0])+(x1[1]-x0[1])*(x1[1]-x0[1]));
				double t2A = std::sqrt((x2[0]-x0[0])*(x2[0]-x0[0])+(x2[1]-x0[1])*(x2[1]-x0[1]));
				//				double length = 2*(2*t1A+root)*std::cos(alpha);
				std::cerr << divCell->index() << " " << k << " " << k2 << "\t" << wL1 << " " << t1A+root << " (" << t1A-root << ")  "
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
	exit(0);
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
	
  //Find intersection with another wall
	//////////////////////////////////////////////////////////////////////
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
	if( flag != 1 && !(flag==2 && vertexFlag) ) {
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
	
	//Check that the division did not mess up the data structure
	//T->checkConnectivity(1);	
}

DivisionShortestPath::DivisionShortestPath(std::vector<double> &paraValue, 
								   std::vector< std::vector<size_t> > &indValue)
{
	if ( paraValue.size() != 4) {
		std::cerr << "DivisionShortestPath::DivisionShortestPath() "
		<< "Four parameters are used V_threshold, Lwall_fraction, Lwall_threshold, and COM (1 = COM, 0 = Random).\n";
		std::exit(EXIT_FAILURE);
	}
	
	if (indValue.size() != 1) {
		std::cerr << "DivisionShortestPath::DivisionShortestPath() "
				<< "First level: Variable indices for volume dependent cell "
				<< "variables are used.\n";
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
	setParameterId(tmp);
}

int DivisionShortestPath::flag(Tissue *T, size_t i,
						 std::vector< std::vector<double> > &cellData,
						 std::vector< std::vector<double> > &wallData,
						 std::vector< std::vector<double> > &vertexData,
						 std::vector< std::vector<double> > &cellDerivs,
						 std::vector< std::vector<double> > &wallDerivs,
						 std::vector< std::vector<double> > &vertexDerivs)
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

void DivisionShortestPath::update(Tissue* T, size_t i,
						    std::vector< std::vector<double> > &cellData,
						    std::vector< std::vector<double> > &wallData,
						    std::vector< std::vector<double> > &vertexData,
						    std::vector< std::vector<double> > &cellDerivs,
						    std::vector< std::vector<double> > &wallDerivs,
						    std::vector< std::vector<double> > &vertexDerivs)
{
	Cell &cell = T->cell(i);
	
	if (vertexData[0].size() != 2)
	{
		std::cerr << "DivisionShortestPath only supports two dimensions.\n";
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
	
	assert (numWallTmp + 3 == T->numWall());
	
	//Change length of new wall between the divided daugther cells
	wallData[numWallTmp][0] *= parameter(1);
	
	//Check that the division did not mess up the data structure
	//T->checkConnectivity(1);
}

std::vector<DivisionShortestPath::Candidate> DivisionShortestPath::getCandidates(Tissue* T, size_t i,
	std::vector< std::vector<double> > &cellData,
	std::vector< std::vector<double> > &wallData,
	std::vector< std::vector<double> > &vertexData,
	std::vector< std::vector<double> > &cellDerivs,
	std::vector< std::vector<double> > &wallDerivs,
	std::vector< std::vector<double> > &vertexDerivs)
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
		o = cell.randomPositionInCell(vertexData);
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
			double beta = M_PI + sigma - alpha;

			double t = (vx * wx + vy * wy) / (vx * vx + vy * vy);
			double tp = t + (1.0 / std::sqrt(vx * vx + vy * vy)) * A * std::sin(alpha - 0.50 * M_PI) / std::sin(alpha);

			double s = (ux * wpx + uy * wpy) / (ux * ux + uy * uy);
			double sp = s + (1.0 / std::sqrt(ux * ux + uy * uy)) * B * std::sin(beta - 0.50 * M_PI) / std::sin(beta);

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

double DivisionShortestPath::astar(double sigma, double A, double B)
{
     double a = 0;
     double b = M_PI;
     double e = b - a;
     double u = f(a, sigma, A, B);
     double v = f(b, sigma, A, B);
     double c;

     if (sign(u) == sign(v)) {
	     return 0;
     }

     for (size_t k = 0; k < 10; ++k) {
          e = 0.5 * e;
          c = a + e;
          double w = f(c, sigma, A, B);

          if (sign(w) != sign(u)) {
               b = c;
               v = w;
          } else {
               a = c;
               u = w;
          }
     }
     return c;
}

double DivisionShortestPath::f(double a, double sigma, double A, double B)
{
     double tmp = - A * std::cos(a) / (std::sin(a) * std::sin(a));
     tmp += B * std::cos(M_PI + sigma - a) / (std::sin(sigma - a) * std::sin(sigma - a));
     return tmp;
}

int DivisionShortestPath::sign(double a)
{
     return (a >= 0) ? +1 : -1;
}



DivisionRandom::DivisionRandom(std::vector<double> &paraValue, std::vector< std::vector<size_t> > &indValue)
{
	if (paraValue.size() != 3) {
		std::cerr << "DivisionRandom::DivisionRandom() "
		<< "Three parameters are used V_threshold, Lwall_fraction, Lwall_threshold." << std::endl;
		exit(EXIT_FAILURE);
	}
	
	if (indValue.size() != 1) {
		std::cerr << "DivisionRandom::DivisionRandom() "
		<< "First level: Variable indices for volume dependent cell "
		<< "variables are used.\n";
		exit(EXIT_FAILURE);
	}
	
	setId("DivisionRandom");
	setNumChange(1);
	setParameter(paraValue);  
	setVariableIndex(indValue);
	
	std::vector<std::string> tmp(numParameter());
	tmp.resize (numParameter());
	tmp[0] = "V_threshold";
	tmp[1] = "Lwall_threshold";
	setParameterId(tmp);
}

int DivisionRandom::flag(Tissue *T, size_t i,
	std::vector< std::vector<double> > &cellData,
	std::vector< std::vector<double> > &wallData,
	std::vector< std::vector<double> > &vertexData,
	std::vector< std::vector<double> > &cellDerivs,
	std::vector< std::vector<double> > &wallDerivs,
	std::vector< std::vector<double> > &vertexDerivs)
{
	if (T->cell(i).calculateVolume(vertexData) > parameter(0)) {
		std::cerr << "Cell " << i << " marked for division with volume " 
		<< T->cell(i).volume() << std::endl;
		return 1;
	} 
	return 0;
}

void DivisionRandom::update(Tissue* T, size_t i,
	std::vector< std::vector<double> > &cellData,
	std::vector< std::vector<double> > &wallData,
	std::vector< std::vector<double> > &vertexData,
	std::vector< std::vector<double> > &cellDerivs,
	std::vector< std::vector<double> > &wallDerivs,
	std::vector< std::vector<double> > &vertexDerivs)
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

int DivisionRandom::random(int n)
{
	double r = myRandom::Rnd();

	int result = (int) floor(n * r);

	if (result == n) {
		return result - 1;
	} else {
		return result;
	}
}

DivisionMainAxis::DivisionMainAxis(std::vector<double> &paraValue, 
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
	setId("DivisionMainAxis");
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

int DivisionMainAxis::flag(Tissue *T, size_t i,
	std::vector< std::vector<double> > &cellData,
	std::vector< std::vector<double> > &wallData,
	std::vector< std::vector<double> > &vertexData,
	std::vector< std::vector<double> > &cellDerivs,
	std::vector< std::vector<double> > &wallDerivs,
	std::vector< std::vector<double> > &vertexDerivs)
{
	if (T->cell(i).calculateVolume(vertexData) > parameter(0)) {
		std::cerr << "Cell " << i << " marked for division with volume " << T->cell(i).volume() << std::endl;
		return 1;
	} else { 
		return 0;
	}
}

void DivisionMainAxis::update(Tissue *T, size_t cellI,
	std::vector< std::vector<double> > &cellData,
	std::vector< std::vector<double> > &wallData,
	std::vector< std::vector<double> > &vertexData,
	std::vector< std::vector<double> > &cellDeriv,
	std::vector< std::vector<double> > &wallDeriv,
	std::vector< std::vector<double> > &vertexDeriv)
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

std::vector<double> DivisionMainAxis::getMainAxis(Cell &cell, std::vector< std::vector<double> > &vertexData)
{
	size_t dimensions = vertexData[0].size();
	size_t numberOfVertices = cell.numVertex();
	
	// Copy vertex data to temporary container and calculate mean values.
	
	std::vector< std::vector<double> > vertices(numberOfVertices);
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
	
	std::vector< std::vector<double> > R(dimensions);
	
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
	
	std::vector< std::vector<double> > candidates;
	
	std::vector< std::vector<double> > V;
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

DivisionVolumeRandomDirectionGiantCells::DivisionVolumeRandomDirectionGiantCells(std::vector<double> &paraValue, 
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

int DivisionVolumeRandomDirectionGiantCells::flag(Tissue *T, size_t i,
	std::vector< std::vector<double> > &cellData,
	std::vector< std::vector<double> > &wallData,
	std::vector< std::vector<double> > &vertexData,
	std::vector< std::vector<double> > &cellDerivs,
	std::vector< std::vector<double> > &wallDerivs,
	std::vector< std::vector<double> > &vertexDerivs)
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

void DivisionVolumeRandomDirectionGiantCells::
update(Tissue *T,size_t cellI,
       std::vector< std::vector<double> > &cellData,
       std::vector< std::vector<double> > &wallData,
       std::vector< std::vector<double> > &vertexData,
       std::vector< std::vector<double> > &cellDeriv,
       std::vector< std::vector<double> > &wallDeriv,
       std::vector< std::vector<double> > &vertexDeriv ) {
  
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
	  com = divCell->randomPositionInCell(vertexData);
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



DivisionShortestPathGiantCells::DivisionShortestPathGiantCells(std::vector<double> &paraValue, 
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

int DivisionShortestPathGiantCells::flag(Tissue *T, size_t i,
						 std::vector< std::vector<double> > &cellData,
						 std::vector< std::vector<double> > &wallData,
						 std::vector< std::vector<double> > &vertexData,
						 std::vector< std::vector<double> > &cellDerivs,
						 std::vector< std::vector<double> > &wallDerivs,
						 std::vector< std::vector<double> > &vertexDerivs)
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

void DivisionShortestPathGiantCells::update(Tissue* T, size_t i,
						    std::vector< std::vector<double> > &cellData,
						    std::vector< std::vector<double> > &wallData,
						    std::vector< std::vector<double> > &vertexData,
						    std::vector< std::vector<double> > &cellDerivs,
						    std::vector< std::vector<double> > &wallDerivs,
						    std::vector< std::vector<double> > &vertexDerivs)
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

std::vector<DivisionShortestPathGiantCells::Candidate> DivisionShortestPathGiantCells::getCandidates(Tissue* T, size_t i,
	std::vector< std::vector<double> > &cellData,
	std::vector< std::vector<double> > &wallData,
	std::vector< std::vector<double> > &vertexData,
	std::vector< std::vector<double> > &cellDerivs,
	std::vector< std::vector<double> > &wallDerivs,
	std::vector< std::vector<double> > &vertexDerivs)
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
		o = cell.randomPositionInCell(vertexData);
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
			double beta = M_PI + sigma - alpha;

			double t = (vx * wx + vy * wy) / (vx * vx + vy * vy);
			double tp = t + (1.0 / std::sqrt(vx * vx + vy * vy)) * A * std::sin(alpha - 0.50 * M_PI) / std::sin(alpha);

			double s = (ux * wpx + uy * wpy) / (ux * ux + uy * uy);
			double sp = s + (1.0 / std::sqrt(ux * ux + uy * uy)) * B * std::sin(beta - 0.50 * M_PI) / std::sin(beta);

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

double DivisionShortestPathGiantCells::astar(double sigma, double A, double B)
{
     double a = 0;
     double b = M_PI;
     double e = b - a;
     double u = f(a, sigma, A, B);
     double v = f(b, sigma, A, B);
     double c;

     if (sign(u) == sign(v)) {
	     return 0;
     }

     for (size_t k = 0; k < 10; ++k) {
          e = 0.5 * e;
          c = a + e;
          double w = f(c, sigma, A, B);

          if (sign(w) != sign(u)) {
               b = c;
               v = w;
          } else {
               a = c;
               u = w;
          }
     }
     return c;
}

double DivisionShortestPathGiantCells::f(double a, double sigma, double A, double B)
{
     double tmp = - A * std::cos(a) / (std::sin(a) * std::sin(a));
     tmp += B * std::cos(M_PI + sigma - a) / (std::sin(sigma - a) * std::sin(sigma - a));
     return tmp;
}

int DivisionShortestPathGiantCells::sign(double a)
{
     return (a >= 0) ? +1 : -1;
}
