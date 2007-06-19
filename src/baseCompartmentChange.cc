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
	//compartmentRemoval.h,compartmentRemoval.cc
  else if(idValue=="RemovalOutsideRadius")
    return new RemovalOutsideRadius(paraValue,indValue);
  else if(idValue=="RemovalOutsideRadiusEpidermis")
    return new RemovalOutsideRadiusEpidermis(paraValue,indValue);
  else if(idValue=="RemovalOutsideMaxDistanceEpidermis")
    return new RemovalOutsideMaxDistanceEpidermis(paraValue,indValue);
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
