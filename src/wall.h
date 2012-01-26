//
// Filename     : wall.h
// Description  : A class describing a one-dimensional wall (for 2d tissue)
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : April 2006
// Revision     : $Id$
//
#ifndef WALL_H
#define WALL_H

#include<utility>
#include<vector>
#include<string>
#include<iostream>
#include<fstream>
#include"vertex.h"

class Cell;

///
/// @brief Defines the properties of a one-dimensional wall element
///
/// Walls are one-dimensional elements connected to two (2D) cells and two (1D) vertices.
/// Walls at the boundary of the tissue are also connected to two cells, but one of these 
/// is defined as a background via an index=-1. The walls hold a specific variable for the length
/// which is the optimal (resting) length, while its actual lenth is defined via the two
/// vertex positions. A wall can also store a vector of variables representing e.g. molecular
/// concentrations.
///
/// @see Tissue::background()
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
  
  ///
  /// @brief Empty constructor
  ///
  /// Only sets the two cellSort variables to 0.
  ///
  /// @see cellSort1()
  ///
  Wall();
  ///
  /// @brief Copy constructor
  ///
  /// Copies the infomration from the provided wall.
  ///
  Wall( const Wall& wallCopy );
  
  ///
  /// @brief Empty destructor.
  ///
  /// Is in the current implementation not doing anything.
  ///
  ~Wall();
  
  ///
  /// @brief Returns the wall index
  ///
  /// The index can be used for finding connection information.
  ///
  inline size_t index() const;
  ///
  /// @brief The wall id
  ///
  /// This is the id/name of the wall. It is optional and the information is not used 
  /// in any important functions.
  ///
  inline std::string id() const;
  ///
  /// @brief The first Cell
  ///
  /// Each wall is connected to two cells (one on each side of the wall).
  ///
  /// @see Cell
  /// @see Tissue::background()
  ///
  inline Cell* cell1() const;
  ///
  /// @brief The second Cell
  ///
  /// @see cell1()
  ///
  inline Cell* cell2() const;
  /// 
  /// @brief The sort direction for the first cell
  ///
  /// The sort direction...
  ///
  inline int cellSort1() const;
  /// 
  /// @brief The sort direction for the second cell
  ///
  /// @see cellSort1()
  ///
  inline int cellSort2() const;
  ///
  /// The first Vertex
  ///
  /// A wall is defined by two (1D) vertices from where spatial information such as actual length
  /// and position of the wall can be extracted.
  ///
  /// @see Vertex
  ///
  inline Vertex* vertex1() const;
  ///
  /// @brief The second Vertex
  ///
  /// @see vertex1()
  ///
  inline Vertex* vertex2() const;
  ///
  /// @brief Returns the wall length variable, ('preferred' but not always actual).
  ///
  /// This variable is used to store length information, e.g. resting length if the wall
  /// is used as a spring. This variable does not have to be the same as the actual length
  /// calculated from the distance between the two vertices defining the wall.
  ///
  inline double length() const;
  ///
  /// @brief Returns the wall length calculated from the vertex positions stored in T.vertex()
  ///
  /// Note that if the vertex positions are updated this value might be obsolete.
  ///
  /// @see lengthFromVertexPosition( DataMatrix&)
  /// 
  double lengthFromVertexPosition();
  ///
  /// @brief Returns the wall length calculated from the vertex positions in vertexData
  ///
  double lengthFromVertexPosition( DataMatrix 
				   &vertexData);
  ///
  /// @brief Returns the number of variables stored by the Wall.
  ///
  inline size_t numVariable() const;
  ///
  /// @brief Sets (resizes) the number of variables stored by the Wall.
  ///
  inline size_t setNumVariable(size_t n);
  ///
  /// @brief Returns variable stored at i in the variable vector in the Wall.
  ///
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
  
  //!Sets the index variable
  inline void setIndex( size_t value );
  //!Sets the cells
  inline void setCell(Cell *v1,Cell *v2);
  //!Sets cell1
  inline void setCell1(Cell *v1);
  //!Sets cell2
  inline void setCell2(Cell *v2);
  ///
  /// @brief Sets the cell sorting directions
  ///
  inline void setCellSort(int value1,int value2);
  ///
  /// @brief Sets the cell sorting direction for cell 1
  ///
  inline void setCellSort1(int value);
  ///
  /// @brief Sets the cell sorting direction for cell 2
  ///
  inline void setCellSort2(int value);
  ///
  /// @brief Sets the vertecis into the pair of vertices defining the wall
  ///
  /// The vertex are provided as pointers to Vertex elements.
  ///
  /// @see vertex1()
  ///
  inline void setVertex(Vertex *v1,Vertex *v2);
  ///
  /// @brief Sets vertex1
  ///
  /// @see setVertex(Vertex*,Vertex*)
  ///
  inline void setVertex1(Vertex *v1);
  ///
  /// @brief Sets vertex2
  ///
  /// @see setVertex(Vertex*,Vertex*)
  ///
  inline void setVertex2(Vertex *v2);
  ///
  /// @brief Returns 1 if the cell is one of the two cells connected to the wall
  ///
  /// It checks whether the cell is one of the two cells connected to the wall. It returns 1
  /// if this is the case, and 0 otherwise.
  ///
  inline int hasCell( Cell *val );
  ///
  /// @brief Returns 1 if the vertex is one of the two vertices connected to the wall
  ///
  /// It checks whether the vertex is one of the two vertices defining the wall. It returns 1
  /// if this is the case, and 0 otherwise.
  ///
  inline int hasVertex( Vertex *val );
  ///
  /// @brief Sets the wall length variable
  ///
  /// @see length()
  ///
  inline void setLength(double val);
  ///
  /// @brief Sets the length to the distance between the two vertices defining the wall
  ///
  /// The distance is calculated from vertex position data stored in the Tissue. Note that
  /// this data may be obsolete if the positions have been updated/moved.
  ///
  /// @see vertex1()
  /// @see Tissue
  /// @see setLengthFromVertexPosition( DataMatrix&) 
  ///
  double setLengthFromVertexPosition();
  ///
  /// @brief Sets the length to the distance between the two vertices defining the wall
  ///
  /// The distance is calculated from vertex position data provided in vertexData.
  ///
  /// @see vertex1()
  /// @see Tissue
  ///
  double setLengthFromVertexPosition( DataMatrix 
				      &vertexData);
  
};

inline size_t Wall::index() const { return index_; }

inline std::string Wall::id() const { return id_; }

inline Cell* Wall::cell1() const {
  return cell_.first;
}

inline Cell* Wall::cell2() const {
  return cell_.second;
}

inline int Wall::cellSort1() const {
	return cellSort_.first;
}

inline int Wall::cellSort2() const {
	return cellSort_.second;
}

inline Vertex* Wall::vertex1() const {
  return vertex_.first;
}

inline Vertex* Wall::vertex2() const {
  return vertex_.second;
}

inline double Wall::length() const { return length_; }

inline size_t Wall::numVariable() const { return variable_.size(); }

inline size_t Wall::setNumVariable(size_t n) { 
  variable_.resize(n); return variable_.size(); }

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

inline void Wall::setIndex( size_t value ) { index_ = value; }

inline void Wall::setCell(Cell *c1,Cell *c2) {
  cell_.first = c1;
  cell_.second =c2;
}

inline void Wall::setCell1(Cell *c1) {
  cell_.first = c1;
}

inline void Wall::setCell2(Cell *c2) {
  cell_.second = c2;
}

inline void Wall::setCellSort(int value1,int value2) {
	cellSort_.first = value1;
	cellSort_.second = value2;
}

inline void Wall::setCellSort1(int value) {
	cellSort_.first = value;
}

inline void Wall::setCellSort2(int value) {
	cellSort_.second = value;
}

inline void Wall::setVertex(Vertex *v1,Vertex *v2) {
  vertex_.first = v1;
  vertex_.second =v2;
}

inline void Wall::setVertex1(Vertex *v1) {
  vertex_.first = v1;
}

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

inline void Wall::setLength(double val) { length_=val; }


#endif
