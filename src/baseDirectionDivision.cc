//
// Filename     : baseDirectionDivision.cc
// Description  : A base class describing directional variable divisions
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : June 2007
// Revision     : $Id:$
//
#include<vector>

#include "baseDirectionDivision.h"
#include "directionDivision.h"

BaseDirectionDivision::~BaseDirectionDivision(){}

//!Factory creator, all creation should be mapped onto this one 
/*! Given the idValue a directionDivision of the defined type is returned
(using new Class).*/
BaseDirectionDivision *
BaseDirectionDivision::createDirectionDivision(std::vector<double> &paraValue,
																					 std::vector< std::vector<size_t> > &indValue, 
																					 std::string idValue ) 
{  
  // All directionDivision classes are in directionDivision.cc .h
  if(idValue=="StaticDirection")
    return new StaticDirection(paraValue,indValue);
  else if(idValue=="PerpendicularDirection")
    return new PerpendicularDirection(paraValue,indValue);
  //Default, if nothing found
  else {
    std::cerr << "\nBaseDirectionDivision::createDirectionDivision() WARNING: DirectionDivisiontype " 
							<< idValue << " not known, no directionDivision created.\n\7";
    exit(-1);
  }
}

//!This creator reads from an open file and then calls for the main creator
BaseDirectionDivision* 
BaseDirectionDivision::createDirectionDivision(std::istream &IN ) {
  
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
  
  return createDirectionDivision(pVal,varIndexVal,idVal);
}

void BaseDirectionDivision::
update(Tissue &T,size_t cellI,
			 std::vector< std::vector<double> > &cellData,
			 std::vector< std::vector<double> > &walldata,
			 std::vector< std::vector<double> > &vertexData,
			 std::vector< std::vector<double> > &cellderivs, 
			 std::vector< std::vector<double> > &wallderivs,
			 std::vector< std::vector<double> > &vertexDerivs ) {
  std::cerr << "BaseDirectionDivision::derivs() should not be used. "
						<< "Should always be mapped onto one of the real types.\n";
  exit(0);
}  

void BaseDirectionDivision::print( std::ofstream &os ) {
  std::cerr << "BaseDirectionDivision::print(ofstream) should not be used. "
						<< "Should always be mapped onto one of the real types.\n";
  exit(0);
}
