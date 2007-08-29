/**
 * Filename     : cell.h
 * Description  : A class describing a two-dimensional cell
 * Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
 * Created      : April 2006
 * Revision     : $Id:$
 */
#ifndef CELL_H
#define CELL_H

#include<assert.h>
#include<cmath>
#include<fstream>
#include<iostream>
#include<list>
#include<string>
#include<vector>

#include"wall.h"

//class Wall;
class Vertex;
class Tissue;

//!Describes the properties of a two-dimensional cell
/*!*/ 
class Cell {
  
 private:
  
  size_t index_;
  std::string id_;           
  
  std::vector<Wall*> wall_;
  std::vector<Vertex*> vertex_;
  std::vector<double> variable_;
  double volume_;
  size_t mitosisFlag_;
  Wall* directionWall_;

  // The vectors obtained from PCA.
  std::vector< std::vector<double> > E;
  
 public:
  
  Cell();
  Cell( const Cell & cellCopy );
  Cell(size_t indexVal,std::string idVal);
  ~Cell();
  
  inline size_t index() const;
  inline std::string id() const;
  inline double volume() const;
  inline size_t mitosisFlag() const;
	
  inline size_t numWall() const;
  inline size_t numVertex() const;
  inline size_t numVariable() const;
  
  inline const std::vector<Wall*> & wall() const;
  inline Wall* wall( size_t k );
  inline Wall& wallRef( size_t k );
	inline int hasWall( Wall *val );
  inline const std::vector<Vertex*> & vertex() const;
  inline Vertex* vertex( size_t v );
  inline int hasVertex( Vertex *val );
  inline const std::vector<double> & variable() const;
  inline double variable(size_t i);
  inline void setId(std::string value);
  inline void setIndex(size_t value);
  inline void setWall( std::vector<Wall*> &val );
	inline void setWall( size_t index,Wall* val);
  inline void addWall( Wall* val );
  inline void setVertex( std::vector<Vertex*> &val );
	inline void setVertex( size_t index,Vertex* val);
  inline void addVertex( Vertex* val );
  inline void setVariable( size_t index,double val );
  inline void addVariable( double val );
  inline int isNeighbor( Cell *neighbor );
	inline Cell* cellNeighbor(size_t k);
  
	void sortWallAndVertex(Tissue &T);

  double calculateVolume( size_t signFlag=0 );
  double calculateVolume( std::vector< std::vector<double> > 
													&vertexData, size_t signFlag=0 );
  std::vector<double> positionFromVertex();
  std::vector<double> 
	positionFromVertex( std::vector< std::vector<double> > &vertexData );	


  // These functions handles projection to a PCAPlane. Important: Make
  // sure you call calculatePCAPlane() before any of the other
  // functions.
  void calculatePCAPlane(std::vector< std::vector<double> > &vertexData);
  std::vector< std::vector<double> > getPCAPlane(void);
  std::vector< std::pair<double, double> > projectVerticesOnPCAPlane(std::vector< std::vector<double> > &vertexData);
  std::vector<double> getNormalToPCAPlane(void);
};

//!Returns the cell index
inline size_t Cell::index() const 
{ 
	return index_; 
}

//!The cell id
inline std::string Cell::id() const 
{ 
	return id_; 
}

//!The cell volume (if updated)
inline double Cell::volume() const 
{ 
	return volume_; 
}

//!A flag that is set if the cell can divide 
inline size_t Cell::mitosisFlag() const 
{ 
	return mitosisFlag_; 
}

//!Number of walls for the cell
inline size_t Cell::numWall() const 
{ 
	return wall_.size(); 
}

//!Number of vertices for the cell
inline size_t Cell::numVertex() const 
{ 
	return vertex_.size(); 
}

//!Number of variables in the cell
inline size_t Cell::numVariable() const 
{ 
	return variable_.size(); 
}

//!Returns a reference to the wall vector
inline const std::vector<Wall*> & Cell::wall() const 
{ 
	return wall_; 
}

//!Returns the wall pointer from position k
inline Wall* Cell::wall( size_t k ) 
{ 
	return wall_[k]; 
}

//!Returns a reference to a wall from position k
inline Wall& Cell::wallRef( size_t k ) 
{ 
	return *wall_[k]; 
}

//!Returns 1 if the cell has the wall
inline int Cell::hasWall( Wall *val ) 
{
  std::vector<Wall*>::iterator it = find(wall_.begin(),wall_.end(),val);
  if( it != wall_.end() ) return 1;
  return 0;
}

//!Returns a reference to the vertex vector
inline const std::vector<Vertex*> & Cell::vertex() const 
{ 
	return vertex_; 
}

//!Returns the vertex pointer from position v
inline Vertex* Cell::vertex( size_t v ) 
{ 
	return vertex_[v]; 
}

//!Returns true if the cell has the vertex
inline int Cell::hasVertex( Vertex *val ) 
{
  std::vector<Vertex*>::iterator it = find(vertex_.begin(),vertex_.end(),val);
  if( it != vertex_.end() ) return 1;
  return 0;
}

//!Returns a reference to the variable vector
inline const std::vector<double> & Cell::variable() const 
{ 
	return variable_; 
}

//!Returns variable i
inline double Cell::variable(size_t i) 
{ 
	return variable_[i];
}

//!Sets the id string
inline void Cell::setId(std::string value) 
{
  id_ = value;
}

//!Sets the index 
inline void Cell::setIndex(size_t value) 
{
  index_ = value;
}

//!Sets pointers to walls
inline void Cell::setWall( std::vector<Wall*> &val ) 
{
  wall_=val;
}

//!Sets a pointer to a wall at index 
inline void Cell::setWall( size_t index,Wall* val) 
{
	assert( index<wall_.size() );
	wall_[index]=val;
}

//!Adds a wall to the cell
inline void Cell::addWall( Wall* val ) 
{
  wall_.push_back(val);
}

//!Sets pointers to vertices
inline void Cell::setVertex( std::vector<Vertex*> &val ) 
{
  vertex_=val;
}

//!Sets a pointer to a vertex at index 
inline void Cell::setVertex( size_t index,Vertex* val) 
{
	assert( index<vertex_.size() );
	vertex_[index]=val;
}

//!Adds a vertex to the cell
inline void Cell::addVertex( Vertex* val ) 
{
  vertex_.push_back(val);
}

//!Adds a variable to the variable vector
inline void Cell::setVariable( size_t index,double val ) 
{
  variable_[index]=val;
}

//!Adds a variable to the variable vector
inline void Cell::addVariable( double val ) 
{
  variable_.push_back(val);
}

//!Checks if any of the walls connects to the given cell
inline int Cell::isNeighbor( Cell *neighbor ) 
{
  for (size_t i=0; i<numWall(); ++i )
    if (wall(i)->cell1()==neighbor || wall(i)->cell2()==neighbor) 
      return 1;
  return 0;
}

///
/// @brief returns the cell neighbor on the other side of wall k
///
inline Cell* Cell::cellNeighbor(size_t k)
{
	return wall(k)->cell1()==this ? wall(k)->cell2() : wall(k)->cell1(); 
}

#endif
