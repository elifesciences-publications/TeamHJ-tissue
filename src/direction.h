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
								std::vector< std::vector<double> > &cellData,
								std::vector< std::vector<double> > &wallData,
								std::vector< std::vector<double> > &vertexData,
								std::vector< std::vector<double> > &cellDerivs,
								std::vector< std::vector<double> > &wallDerivs,
								std::vector< std::vector<double> > &vertexDerivs );

  void update(Tissue &T,double step,
							std::vector< std::vector<double> > &cellData,
							std::vector< std::vector<double> > &wallData,
							std::vector< std::vector<double> > &vertexData,
							std::vector< std::vector<double> > &cellDerivs,
							std::vector< std::vector<double> > &wallDerivs,
							std::vector< std::vector<double> > &vertexDerivs );

  void divide(Tissue &T,size_t cellI,
							std::vector< std::vector<double> > &cellData,
							std::vector< std::vector<double> > &wallData,
							std::vector< std::vector<double> > &vertexData,
							std::vector< std::vector<double> > &cellDerivs,
							std::vector< std::vector<double> > &wallDerivs,
							std::vector< std::vector<double> > &vertexDerivs );
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