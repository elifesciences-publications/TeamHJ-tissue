/**
 * Filename     : wall.h
 * Description  : A class describing a one-dimensional wall (for 2d tissue)
 * Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
 * Created      : April 2006
 * Revision     : $Id:$
 */
#ifndef WALL_H
#define WALL_H

#include<utility>
#include<vector>
#include<string>
#include<iostream>
#include<fstream>
#include"vertex.h"

class Cell;

//!Describes the properties of a one-dimensional wall
/*!*/ 
class Wall {
  
 private:
  
  size_t index_;
  std::string id_;           
  
  std::pair<Cell*,Cell*> cell_;
  std::pair<Vertex*,Vertex*> vertex_;
  double length_;
  std::vector<double> variable_;

 public:
  
  Wall();
  Wall( const Wall& wallCopy );
  
  ~Wall();
  
  inline size_t index() const;
  inline std::string id() const;
  inline Cell* cell1() const;
  inline Cell* cell2() const;
  inline Vertex* vertex1() const;
  inline Vertex* vertex2() const;
  inline double length() const;
  inline size_t numVariable() const;
  inline double variable(size_t i) const;

  inline void setVariable(std::vector<double> variable);

  inline void setIndex( size_t value );
  inline void setCell(Cell *v1,Cell *v2);
  inline void setCell1(Cell *v1);
  inline void setCell2(Cell *v2);
  inline void setVertex(Vertex *v1,Vertex *v2);
  inline void setVertex1(Vertex *v1);
  inline void setVertex2(Vertex *v2);
	inline int hasCell( Cell *val );
	inline int hasVertex( Vertex *val );
  inline void setLength(double val);
  double setLengthFromVertexPosition();
  double setLengthFromVertexPosition( std::vector< std::vector<double> > 
				      &vertexData);

};

//!Returns the wall index
inline size_t Wall::index() const { return index_; }

//!The wall id
inline std::string Wall::id() const { return id_; }

//!The first Cell
inline Cell* Wall::cell1() const {
  return cell_.first;
}

//!The second Cell
inline Cell* Wall::cell2() const {
  return cell_.second;
}

//!The first Vertex
inline Vertex* Wall::vertex1() const {
  return vertex_.first;
}

//!The second Vertex
inline Vertex* Wall::vertex2() const {
  return vertex_.second;
}

//!The wall length (preferred but not always actual)
inline double Wall::length() const { return length_; }

//!Returns the number of variables
inline size_t Wall::numVariable() const { return variable_.size(); }

//!Returns variable i
inline double Wall::variable(size_t i) const { return variable_[i];}

inline void Wall::setVariable(std::vector<double> variable)
{
	variable_ = variable;
}


//!Sets the index variable
inline void Wall::setIndex( size_t value ) { index_ = value; }

//!Sets the cells
inline void Wall::setCell(Cell *c1,Cell *c2) {
  cell_.first = c1;
  cell_.second =c2;
}

//!Sets cell1
inline void Wall::setCell1(Cell *c1) {
  cell_.first = c1;
}

//!Sets cell2
inline void Wall::setCell2(Cell *c2) {
  cell_.second = c2;
}

//!Sets the vertecis
inline void Wall::setVertex(Vertex *v1,Vertex *v2) {
  vertex_.first = v1;
  vertex_.second =v2;
}

//!Sets vertex1
inline void Wall::setVertex1(Vertex *v1) {
  vertex_.first = v1;
}

//!Sets vertex2
inline void Wall::setVertex2(Vertex *v2) {
  vertex_.second = v2;
}

inline int Wall::hasCell( Cell *val ) 
{
	if( cell1()==val || cell2()==val ) 
		return 1;
	return 0;
}

inline int Wall::hasVertex( Vertex *val ) 
{
	if( vertex1()==val || vertex2()==val ) 
		return 1;
	return 0;
}

//!Sets the wall length
inline void Wall::setLength(double val) { length_=val; }


#endif
