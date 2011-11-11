//
// Filename     : direction.h
// Description  : A class describing cell directions
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : June 2007
// Revision     : $Id:$
//
#ifndef DIRECTION_H
#define DIRECTION_H

#include<vector>
#include<string>
#include<iostream>
#include<fstream>
#include"baseDirectionUpdate.h"
#include"baseDirectionDivision.h"

class Tissue;

///
/// @brief This class defines all rules for specific directions in the cells.
///
/// Directions can be defined for cells in a tissue to hold vectorial information, e.g. MT/MF directions.
/// This class is a base class holding two main update base classes. One for dynamic updates
/// during a simulation, and on for rules relating the direction to cell divisions.
///
/// @see BaseDirectionUpdate
/// @see BaseDirectionDivision
/// @see Tissue
/// @see Cell
/// 
class Direction {
  
 private:
  
  std::string id_;
  BaseDirectionUpdate* directionUpdate_;           
  BaseDirectionDivision* directionDivision_;           
  
 public:
  
  Direction();
  Direction( Direction &directionCopy );
  Direction( char *inFile );
  Direction( const std::string &inFile );
  Direction( std::ifstream &IN );
  
  ~Direction();
  
  // Get values
  inline std::string id() const;
  
  inline BaseDirectionUpdate* directionUpdate();
  inline BaseDirectionDivision* directionDivision();
  
  // Set values
  inline void setId(const std::string &value);
  inline void setDirectionUpdate(BaseDirectionUpdate* value);
  inline void setDirectionDivision(BaseDirectionDivision* value);
  
  //Other functions
  int readDirection( std::ifstream &IN);
  int addUpdate( std::istream &IN );
  int addDivision( std::istream &IN );
  void initiate(Tissue &T,
		DataMatrix &cellData,
		DataMatrix &wallData,
		DataMatrix &vertexData,
		DataMatrix &cellDerivs,
		DataMatrix &wallDerivs,
		DataMatrix &vertexDerivs );
  
  void update(Tissue &T,double step,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      DataMatrix &cellDerivs,
	      DataMatrix &wallDerivs,
	      DataMatrix &vertexDerivs );
  
  void divide(Tissue &T,size_t cellI,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      DataMatrix &cellDerivs,
	      DataMatrix &wallDerivs,
	      DataMatrix &vertexDerivs );
};

//!Returns the id string
inline std::string Direction::id() const {
  return id_;
}

//!Returns a directionUpdate
inline BaseDirectionUpdate* Direction::directionUpdate()
{
  return directionUpdate_;
}

//!Returns a directionDivision
inline BaseDirectionDivision* Direction::directionDivision()
{
  return directionDivision_;
}

//!Sets the id string
inline void Direction::setId(const std::string &value) {
  id_=value;
}

//!Sets the directionUpdate
inline void Direction::setDirectionUpdate(BaseDirectionUpdate* value) {
  directionUpdate_=value;
}

//!Sets the directionDivision
inline void Direction::setDirectionDivision(BaseDirectionDivision* value) {
  directionDivision_=value;
}

#endif
