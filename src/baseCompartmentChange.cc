/**
 * Filename     : baseCompartmentChange.cc
 * Description  : A base class describing variable updates
 * Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
 * Created      : October 2003
 * Revision     : $Id: baseCompartmentChange.cc,v 1.25 2006/03/18 00:05:14 henrik Exp $
 */
#include<vector>

#include"baseCompartmentChange.h"
#include"compartmentDivision.h"
#include"compartmentRemoval.h"

BaseCompartmentChange::~BaseCompartmentChange(){}

//!Factory creator, all creation should be mapped onto this one 
/*! Given the idValue a compartmentChange of the defined type is returned
  (using new Class).*/
BaseCompartmentChange *
BaseCompartmentChange::
createCompartmentChange(std::vector<double> &paraValue,
												std::vector< std::vector<size_t> > &indValue, 
												std::string idValue ) {
  
  //Cell divisions
  //compartmentDivision.h,compartmentDivision.cc
  if(idValue=="DivisionVolumeViaLongestWall")
    return new DivisionVolumeViaLongestWall(paraValue,indValue);
  else if(idValue=="DivisionVolumeViaLongestWall3D")
    return new DivisionVolumeViaLongestWall3D(paraValue,indValue);
	else if(idValue=="DivisionVolumeViaStrain")
    return new DivisionVolumeViaStrain(paraValue,indValue);
	else if(idValue=="DivisionVolumeViaDirection")
    return new DivisionVolumeViaDirection(paraValue,indValue);
	else if(idValue=="DivisionVolumeRandomDirection")
    return new DivisionVolumeRandomDirection(paraValue,indValue);
	else if(idValue=="DivisionVolumeViaShortestPath")
    return new DivisionVolumeViaShortestPath(paraValue,indValue);
	//compartmentRemoval.h,compartmentRemoval.cc
  else if(idValue=="RemovalOutsideRadius")
    return new RemovalOutsideRadius(paraValue,indValue);
  else if(idValue=="RemovalOutsideRadiusEpidermis")
    return new RemovalOutsideRadiusEpidermis(paraValue,indValue);
  else if(idValue=="RemovalOutsideMaxDistanceEpidermis")
    return new RemovalOutsideMaxDistanceEpidermis(paraValue,indValue);
  else if (idValue == "DivisionForceDirection")
	  return new DivisionForceDirection(paraValue, indValue);
  else if (idValue == "DivisionShortestPath")
	  return new DivisionShortestPath(paraValue, indValue);

  //Default, if nothing found
  else {
    std::cerr << "\nBaseCompartmentChange::createCompartmentChange()"
							<< " WARNING: CompartmentChangetype " 
							<< idValue << " not known, no compartmentChange created.\n\7";
    exit(-1);
  }
}

//!This creator reads from an open file and then calls for the main creator
BaseCompartmentChange* 
BaseCompartmentChange::createCompartmentChange(std::istream &IN ) {
  
  std::string idVal;
  size_t pNum,levelNum;
  IN >> idVal;
  IN >> pNum;
  IN >> levelNum;
  std::vector<size_t> varIndexNum( levelNum );
  for( size_t i=0 ; i<levelNum ; i++ )
    IN >> varIndexNum[i];
  
  std::vector<double> pVal( pNum );
  for( size_t i=0 ; i<pNum ; i++ )
    IN >> pVal[i];
  
  std::vector< std::vector<size_t> > varIndexVal( levelNum );
  for( size_t i=0 ; i<levelNum ; i++ )
    varIndexVal[i].resize( varIndexNum[i] );
  
  for( size_t i=0 ; i<levelNum ; i++ )
    for( size_t j=0 ; j<varIndexNum[i] ; j++ )
      IN >> varIndexVal[i][j];
  
  return createCompartmentChange(pVal,varIndexVal,idVal);
}

int BaseCompartmentChange::
flag(Tissue *T,size_t i,std::vector< std::vector<double> > &cellData,
     std::vector< std::vector<double> > &walldata,
     std::vector< std::vector<double> > &vertexData,
     std::vector< std::vector<double> > &cellderivs, 
     std::vector< std::vector<double> > &wallderivs,
     std::vector< std::vector<double> > &vertexDerivs ) {
  std::cerr << "BaseCompartmentChange::flag() should not be used. "
						<< "Should always be mapped onto one of the real types.\n";
  exit(0);
}  

void BaseCompartmentChange::
update(Tissue *T,size_t i,std::vector< std::vector<double> > &cellData,
       std::vector< std::vector<double> > &walldata,
       std::vector< std::vector<double> > &vertexData,
       std::vector< std::vector<double> > &cellderivs, 
       std::vector< std::vector<double> > &wallderivs,
       std::vector< std::vector<double> > &vertexDerivs ) {
  std::cerr << "BaseCompartmentChange::update() should not be used. "
						<< "Should always be mapped onto one of the real types.\n";
  exit(0);
}  

void BaseCompartmentChange::
printCellWallError(std::vector< std::vector<double> > &vertexData,
									 Cell *divCell, 
									 std::vector<size_t> &w3Tmp, 
									 size_t &wI, 
									 size_t &w3I,
									 std::ostream &os) {
	
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
}

void BaseCompartmentChange::
findSecondDivisionWall(std::vector< std::vector<double> > &vertexData, 
											 Cell *divCell, size_t &wI, size_t &w3I, 
											 size_t &flag, size_t &vertexFlag,
											 std::vector<double> &v1Pos, 
											 std::vector<double> &nW2, 
											 std::vector<size_t> &w3Tmp, 
											 std::vector<double> &w3tTmp) {
	
	size_t dimension=vertexData[0].size();
	if (dimension==2) {
		for (size_t k=0; k<divCell->numWall(); ++k) {
			if (k!=wI) {
				size_t v1w3Itmp = divCell->wall(k)->vertex1()->index();
				size_t v2w3Itmp = divCell->wall(k)->vertex2()->index();
				std::vector<double> w3(dimension),w0(dimension);
				for (size_t d=0; d<dimension; ++d) {
					w3[d] = vertexData[v2w3Itmp][d]-vertexData[v1w3Itmp][d];
					w0[d] = v1Pos[d]-vertexData[v1w3Itmp][d];
				}
				double a=0.0,b=0.0,c=0.0,d=0.0,e=0.0;//a=1.0
				for (size_t dim=0; dim<dimension; ++dim) {
					a += nW2[dim]*nW2[dim];
					b += nW2[dim]*w3[dim];
					c += w3[dim]*w3[dim];
					d += nW2[dim]*w0[dim];
					e += w3[dim]*w0[dim];
				}
				double fac=a*c-b*b;//a*c-b*b
				if (fac>1e-10) {//else parallell and not applicable
					fac = 1.0/fac;
					//double s = fac*(b*e-c*d);
					double t = fac*(a*e-b*d);//fac*(a*e-b*d)
					if (t>0.0 && t<=1.0) {//within wall
						//double dx0 = w0[0] +fac*((b*e-c*d)*nW2[0]+()*w3[0]); 					
						++flag;
						if (t==1.0)
							vertexFlag++;
						w3I = k;
						w3Tmp.push_back(k);
						w3tTmp.push_back(t);
					}
				}
			}
		}
	}
}
