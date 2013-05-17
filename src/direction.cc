//
// Filename     : direction.cc
// Description  : A class describing a direction
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : June 2007
// Revision     : $Id:$
//
#include <cstdlib>
#include <fstream>
#include <iostream>

#include "baseDirectionUpdate.h"
#include "baseDirectionDivision.h"
#include "direction.h"

Direction::Direction() 
{
    directionUpdate_ = NULL;
    directionDivision_ = NULL;
  //std::cerr << "Direction::Direction().\n";
}

Direction::Direction( Direction & directionCopy ) 
{ 
  //std::cerr << "Direction::Direction(Direction).\n";
  setId( directionCopy.id() );
  setDirectionUpdate( directionCopy.directionUpdate() );
  setDirectionDivision( directionCopy.directionDivision() );
}

Direction::Direction( char *inFile ) 
{  
  std::ifstream IN( inFile );
  if( !IN ) {
    std::cerr << "Direction::Direction() - "
	      << "Cannot open file " << inFile << "\n\n\7";exit(-1);}
  readDirection(IN);
}

Direction::Direction( const std::string &inFile ) 
{  
  const char *tmp = inFile.c_str();
  std::ifstream IN( tmp );
  if( !IN ) {
    std::cerr << "Direction::Direction() - "
							<< "Cannot open file " << inFile << "\n\n\7";exit(-1);}
  readDirection(IN);
}

Direction::Direction( std::ifstream &IN ) 
{  
  readDirection(IN);
}

Direction::~Direction() 
{
	delete directionUpdate_;
	delete directionDivision_;
}

//!Read a direction from an open filestream
int Direction::readDirection( std::ifstream &IN ) 
{  
	if( !IN ) return -1;
	if( addUpdate(IN) || addDivision(IN) )
		return -1;
  std::string idVal=directionUpdate()->id()+"_"+directionDivision()->id();
	return 0;
}

//!Adds an update rule
int Direction::addUpdate( std::istream &IN ) {
  if( !IN )
    return -1;
  directionUpdate_ = BaseDirectionUpdate::createDirectionUpdate(IN);
  return 0;
}

//!Adds an division rule
int Direction::addDivision( std::istream &IN ) {
  if( !IN )
    return -1;
  directionDivision_ = BaseDirectionDivision::createDirectionDivision(IN);
  return 0;
}

//!Calls directionUpdate for initiation of the direction
void Direction::initiate(Tissue &T, 
			 DataMatrix &cellData,
			 DataMatrix &wallData,
			 DataMatrix &vertexData,
			 DataMatrix &cellDerivs,
			 DataMatrix &wallDerivs,
			 DataMatrix &vertexDerivs )
{
  directionUpdate()->initiate(T,cellData,wallData,vertexData,cellDerivs,
			      wallDerivs,vertexDerivs);
}

//!Calls directionUpdate for update of the direction (during simulation)
void Direction::update(Tissue &T,double step,
		       DataMatrix &cellData,
		       DataMatrix &wallData,
		       DataMatrix &vertexData,
		       DataMatrix &cellDerivs,
		       DataMatrix &wallDerivs,
		       DataMatrix &vertexDerivs )
{
  directionUpdate()->update(T,step,cellData,wallData,vertexData,cellDerivs,
			    wallDerivs,vertexDerivs);
}

//!Calls directionDivision for direction division rule
void Direction::divide(Tissue &T,size_t cellI,
		       DataMatrix &cellData,
		       DataMatrix &wallData,
		       DataMatrix &vertexData,
		       DataMatrix &cellDerivs,
		       DataMatrix &wallDerivs,
											 DataMatrix &vertexDerivs ) 
{
	directionDivision()->update(T,cellI,cellData,wallData,vertexData,
				    cellDerivs,wallDerivs,vertexDerivs);
}
