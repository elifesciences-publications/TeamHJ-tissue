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
	std::pair<int,int> cellSort_;
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
	inline int cellSort1() const;
	inline int cellSort2() const;
  inline Vertex* vertex1() const;
  inline Vertex* vertex2() const;
  inline double length() const;
	///
	/// @brief Returns the wall length calculated from the vertex positions stored in T.vertex()
	///
  double lengthFromVertexPosition();
	///
	/// @brief Returns the wall length calculated from the vertex positions in vertexData
	///
  double lengthFromVertexPosition( std::vector< std::vector<double> > 
																	 &vertexData);
  inline size_t numVariable() const;
  inline double variable(size_t i) const;
	///
	/// @brief Sets all wall variables to values provided in the vector
	///
  inline void setVariable(std::vector<double> variable);
	///
	/// @brief Sets a wall variable with index iVal to value provided in val
	///
  inline void setVariable(size_t iVal,double val);	
	///
	/// @brief Adds a new variable to the variable vector
	///
	inline void addVariable(double val);

  inline void setIndex( size_t value );
  inline void setCell(Cell *v1,Cell *v2);
  inline void setCell1(Cell *v1);
  inline void setCell2(Cell *v2);
	inline void setCellSort(int value1,int value2);
	inline void setCellSort1(int value);
	inline void setCellSort2(int value);
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

/// 
/// @brief The sort direction for the first cell
///
inline int Wall::cellSort1() const {
	return cellSort_.first;
}

/// 
/// @brief The sort direction for the second cell
///
inline int Wall::cellSort2() const {
	return cellSort_.second;
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

inline void Wall::setVariable(size_t iVal,double val)
{
	variable_[iVal] = val;
}

inline void Wall::setVariable(std::vector<double> variable)
{
	variable_ = variable;
}

inline void Wall::addVariable(double val)
{
	variable_.push_back(val);
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

///
/// @brief Sets the cell sorting directions
///
inline void Wall::setCellSort(int value1,int value2) {
	cellSort_.first = value1;
	cellSort_.second = value2;
}

///
/// @brief Sets the cell sorting direction for cell 1
///
inline void Wall::setCellSort1(int value) {
	cellSort_.first = value;
}

///
/// @brief Sets the cell sorting direction for cell 2
///
inline void Wall::setCellSort2(int value) {
	cellSort_.second = value;
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
