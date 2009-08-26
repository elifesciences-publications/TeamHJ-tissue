//
// Filename     : baseDirectionUpdate.cc
// Description  : A base class describing directional variable updates
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : October 2003
// Revision     : $Id: baseReaction.cc,v 1.25 2006/03/18 00:05:14 henrik Exp $
//
#include<vector>

#include"baseDirectionUpdate.h"
#include"directionUpdate.h"

BaseDirectionUpdate::~BaseDirectionUpdate(){}

//!Factory creator, all creation should be mapped onto this one 
/*! Given the idValue a directionUpdate of the defined type is returned
(using new Class).*/
BaseDirectionUpdate *
BaseDirectionUpdate::createDirectionUpdate(std::vector<double> &paraValue,
																					 std::vector< std::vector<size_t> > &indValue, 
																					 std::string idValue ) 
{  
  // All directionUpdate classes are in directionUpdate.cc .h
  if(idValue=="StaticDirection")
    return new StaticDirection(paraValue,indValue);
  else if(idValue=="WallDirection")
    return new WallDirection(paraValue,indValue);
  else if(idValue=="StrainDirection")
    return new StrainDirection(paraValue,indValue);
  else if(idValue=="StrainDirectionWall")
    return new StrainDirectionWall(paraValue,indValue);
  else if(idValue=="GradientDirection")
    return new GradientDirection(paraValue,indValue);
  else if (idValue == "WallStressDirection")
	  return new WallStressDirection(paraValue, indValue);
  else if (idValue == "DoubleWallStressDirection")
	  return new DoubleWallStressDirection(paraValue, indValue);
  else if (idValue == "StretchDirection")
	  return new StretchDirection(paraValue, indValue);
  else if (idValue == "PCAPlaneDirection")
	  return new PCAPlaneDirection(paraValue, indValue);
  else if (idValue == "VertexStressDirection")
	  return new VertexStressDirection(paraValue, indValue);
  else if (idValue == "ForceDirection") {
		std::cerr << "ForceDirection renamed into WallStressDirection." << std::endl;
		exit(-1);
	}
  //Default, if nothing found
  else {
	  std::cerr << "\nBaseDirectionUpdate::createDirectionUpdate() WARNING: DirectionUpdatetype "
							<< idValue << " not known, no directionUpdate created.\n\7";
	  exit(-1);
  }
}

//!This creator reads from an open file and then calls for the main creator
BaseDirectionUpdate* 
BaseDirectionUpdate::createDirectionUpdate(std::istream &IN) {
  
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
  
  return createDirectionUpdate(pVal,varIndexVal,idVal);
}

void BaseDirectionUpdate::
initiate(Tissue &T,
				 std::vector< std::vector<double> > &cellData,
				 std::vector< std::vector<double> > &walldata,
				 std::vector< std::vector<double> > &vertexData,
				 std::vector< std::vector<double> > &cellderivs, 
				 std::vector< std::vector<double> > &wallderivs,
				 std::vector< std::vector<double> > &vertexDerivs ) {
  std::cerr << "BaseDirectionUpdate::derivs() should not be used. "
						<< "Should always be mapped onto one of the real types.\n";
  exit(0);
}  

void BaseDirectionUpdate::
update(Tissue &T, double h,
			 std::vector< std::vector<double> > &cellData,
			 std::vector< std::vector<double> > &walldata,
			 std::vector< std::vector<double> > &vertexData,
			 std::vector< std::vector<double> > &cellderivs, 
			 std::vector< std::vector<double> > &wallderivs,
			 std::vector< std::vector<double> > &vertexDerivs ) {
  std::cerr << "BaseDirectionUpdate::derivs() should not be used. "
						<< "Should always be mapped onto one of the real types.\n";
  exit(0);
}  

void BaseDirectionUpdate::print( std::ofstream &os ) {
  std::cerr << "BaseDirectionUpdate::print(ofstream) should not be used. "
						<< "Should always be mapped onto one of the real types.\n";
  exit(0);
}
