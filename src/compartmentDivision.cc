/**
 * Filename     : compartmentDivision.cc
 * Description  : Classes describing compartmentDivision updates
 * Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
 * Created      : July 2006
 * Revision     : $Id:$
 */
#include "baseCompartmentChange.h"
#include "compartmentDivision.h"
#include "myRandom.h"

//!Constructor
DivisionVolumeViaLongestWall::
DivisionVolumeViaLongestWall(std::vector<double> &paraValue, 
			     std::vector< std::vector<size_t> > 
			     &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //////////////////////////////////////////////////////////////////////
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
  //Set the variable values
  //////////////////////////////////////////////////////////////////////
  setId("DivisionVolumeViaLongestWall");
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

//! Updates the dividing cell by adding a prependicular wall from the longest
/*! 
 */
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
	assert( dimension==2 );
	
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
		std::cerr << "divideVolumeViaLongestWall::update Warning"
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
	T->checkConnectivity(1);	
}

//!Constructor
DivisionVolumeViaLongestWall3D::
DivisionVolumeViaLongestWall3D(std::vector<double> &paraValue, 
															 std::vector< std::vector<size_t> > 
															 &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //////////////////////////////////////////////////////////////////////
  if( paraValue.size()!=2 ) {
    std::cerr << "DivisionVolumeViaLongestWall3D::"
							<< "DivisionVolumeViaLongestWall3D() "
							<< "Two parameters used V_threshold, LWall_frac\n";
    exit(0);
  }
  if( indValue.size() != 1 ) {
    std::cerr << "DivisionVolumeViaLongestWall3D::"
							<< "DivisionVolumeViaLongestWall3D() "
							<< "Variable indices for volume dependent cell "
							<< "variables is used.\n";
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
    v2Pos[d] = 0.5*(vertexData[v1w3I][d]+vertexData[v2w3I][d]);
	
  //Add one cell, three walls, and two vertices
  //////////////////////////////////////////////////////////////////////
	//Save number of walls
	size_t numWallTmp=wallData.size();
	assert( numWallTmp==T->numWall() );
	//Divide
	T->divideCell(divCell,wI,w3I,v1Pos,v2Pos,cellData,wallData,vertexData,
									 cellDeriv,wallDeriv,vertexDeriv,variableIndex(0));
	assert( numWallTmp+3 == T->numWall() );
	
	//Change length of new wall between the divided daugther cells 
	wallData[numWallTmp][0] *= parameter(1);
	
	//Check that the division did not messed up the data structure
	T->checkConnectivity(1);	
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
	T->checkConnectivity(1);		
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
	size_t numV = divCell->numVertex();
  assert( divCell->numWall() > 2 );
	assert( dimension==2 );
	
	std::vector<double> xMean(dimension),yMean(dimension);
	
	for( size_t i=0 ; i<numV ; ++i ) {
		size_t vI = divCell->vertex(i)->index();
		xMean[0] += vertexData[vI][0];
		xMean[1] += vertexData[vI][1];
	}
	xMean[0] /= numV;
	xMean[1] /= numV;
	
	std::vector<double> n(dimension);
	if( cellData[cellI][variableIndex(0,0)+dimension] ) {
		//Perpendicular to given direction
		n[0]=cellData[cellI][variableIndex(0,0)];
		n[1]=cellData[cellI][variableIndex(0,0)+1];		
		if( parameter(3) != 1.0 ) {
			n[0]=cellData[cellI][variableIndex(0,0)+1];
			n[1]=-cellData[cellI][variableIndex(0,0)];
		}
	}
	else {
		//Random
		double phi=2*3.14*myRandom::Rnd();
		n[0] = std::sin(phi);
		n[1] = std::cos(phi);
	}
	
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
	std::cerr << "Dividing cell has " << divCell->numWall() << " walls."
						<< std::endl;
  for( size_t k=0 ; k<divCell->numWall() ; ++k ) {
		size_t v1Tmp = divCell->wall(k)->vertex1()->index();
		size_t v2Tmp = divCell->wall(k)->vertex2()->index();
		std::vector<double> w3(dimension),w0(dimension);
		for( size_t dim=0 ; dim<dimension ; ++dim ) {
			w3[dim] = vertexData[v2Tmp][dim]-vertexData[v1Tmp][dim];
			w0[dim] = xMean[dim]-vertexData[v1Tmp][dim];
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
		if( fac>1e-10 ) {//else parallell and not applicable
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
			//else {
			//std::cerr << "Dividing cell " << divCell->index() << " not via wall "
			//					<< k << " at t=" << t << std::endl;
			//}
		}		
		//else {
		//std::cerr << "Dividing cell " << divCell->index() << " not via wall "
		//					<< k << " since parallell." << std::endl;
		//}
	}
  assert( wI[1] != divCell->numWall() && wI[0] != wI[1] );
	if( flag != 2 ) {
		std::cerr << "divideVolumeViaDirection::update() Warning:"
							<< " not two, but " << flag << " walls chosen as "
							<< "connection for cell " 
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
								cellDeriv,wallDeriv,vertexDeriv,variableIndex(1),
								parameter(2));
	assert( numWallTmp+3 == T->numWall() );
	
	//Change length of new wall between the divided daugther cells 
	wallData[numWallTmp][0] *= parameter(1);
	
	//Check that the division did not mess up the data structure
	T->checkConnectivity(1);		
}

//!Constructor
DivisionVolumeRandomDirection::
DivisionVolumeRandomDirection(std::vector<double> &paraValue, 
			     std::vector< std::vector<size_t> > 
			     &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //////////////////////////////////////////////////////////////////////
  if( paraValue.size()!=3 ) {
    std::cerr << "DivisionVolumeRandomDirection::"
							<< "DivisionVolumeRandomDirection() "
							<< "Two parameters used V_threshold, LWall_frac, and "
							<< "Lwall_threshold" << std::endl;
    exit(0);
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

//! Updates the dividing cell by adding a prependicular wall from the longest
/*! 
 */
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
	size_t numV = divCell->numVertex();
  assert( divCell->numWall() > 2 );
	assert( dimension==2 );
	
	std::vector<double> xMean(dimension),yMean(dimension);
	
	for( size_t i=0 ; i<numV ; ++i ) {
		size_t vI = divCell->vertex(i)->index();
		xMean[0] += vertexData[vI][0];
		xMean[1] += vertexData[vI][1];
	}
	xMean[0] /= numV;
	xMean[1] /= numV;
	
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
	T->checkConnectivity(1);		
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
  size_t wI=0, wI2=1;
  for (size_t k=0; k<divCell->numWall()-1; ++k) {
		for (size_t k2=k+1; k<divCell->numWall(); ++k) {
			std::vector<double> x0(dimension),x1(dimension),x2(dimension),n1(dimension),n2(dimension);
			x1[0] = divCell->wall(k)->vertex1().position(0);
			x1[1] = divCell->wall(k)->vertex1().position(1);
			x2[0] = divCell->wall(k2)->vertex1().position(0);
			x2[1] = divCell->wall(k2)->vertex1().position(1);
			n1[0] = divCell->wall(k)->vertex2().position(0)-x1[0];
			n1[1] = divCell->wall(k)->vertex2().position(1)-x1[1];
			n2[0] = divCell->wall(k2)->vertex2().position(0)-x2[0];
			n2[1] = divCell->wall(k2)->vertex2().position(1)-x2[1];
			double wL1 = std::sqrt(n1[0]*n1[0]+n1[1]*n1[1]);
			double wL2 = std::sqrt(n2[0]*n2[0]+n2[1]*n2[1]);
			
			double t1=0.0, t2=0.0;
			double denominator = n2[0]*n2[1]-n2[0]*n1[1];
			if (denominator != 0.0) {
				t1 = (n2[0]*(x2[1]-x1[1])-n2[1]*(x2[0]-x1[0]))/denominator;
				if (t1>0.0) {
					t1 = wL1-t1;
					x1[0] = divCell->wall(k)->vertex2().position(0);
					x1[1] = divCell->wall(k)->vertex2().position(1);
					n1[0] = -n1[0];
					n1[1] = -n1[1];
				}
				t2 = (n1[0]*(x2[1]-x1[1])-n1[1]*(x2[0]-x1[0]))/denominator;
				if (t2>0.0) {
					t2 = wL2-t2;
					x2[0] = divCell->wall(k2)->vertex2().position(0);
					x2[1] = divCell->wall(k2)->vertex2().position(1);
					n2[0] = -n2[0];
					n2[1] = -n2[1];
				}
				x0[0] = x1[0] + n1[0]*t1;
				x0[1] = x1[1] + n1[1]*t1;
				double n1n2 = n1[0]*n2[0]+n1[1]*n2[1];
				double alpha = std::acos(n1n2/(wL1*wL2));
				double A0 = 0.5*((x1[0]-x0[0])*(x2[1]-x0[1])-
												 (x1[1]-x0[1])*(x2[0]-x0[0]));
				double A2;
				double Ahalf = divCell->calculateVolume(vertexData);
				double root = std::sqrt((A0+Ahalf-A2)/(std::cos(alpha)*std::sin(alpha)));
				double t1A = std::sqrt((x1[0]-x0[0])*(x1[0]-x0[0])+(x1[1]-x0[1])*(x1[1]-x0[1]));
				double t2A = std::sqrt((x2[0]-x0[0])*(x2[0]-x0[0])+(x2[1]-x0[1])*(x2[1]-x0[1]));
				double length = 2*(2*t1A+root)*std::cos(alpha);
			}
			else {
			}


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
	T->checkConnectivity(1);	
}

