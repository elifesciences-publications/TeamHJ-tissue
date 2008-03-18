/**
 * Filename     : vertex.h
 * Description  : A class describing a vertex
 * Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
 * Created      : April 2006
 * Revision     : $Id:$
 */
#ifndef VERTEX_H
#define VERTEX_H

#include<vector>
#include<list>
#include<string>
#include<iostream>
#include<fstream>

class Cell;
class Wall;

//!Describes the properties of a vertex including updatable position
/*!*/ 
class Vertex {
  
 private:
  
  size_t index_;
  std::string id_;           
  
  std::vector<Cell*> cell_;
  std::vector<Wall*> wall_;
  
  std::vector<double> position_;
  
  
 public:
  
  Vertex();
  Vertex( const Vertex & cellCopy );
  
  ~Vertex();
  
  //!Returns the cell index
  inline size_t index() const;
  //!The cell id
  inline std::string id() const;

  //!Number of cells connected to the vertex
  inline size_t numCell() const;
  //!Number of walls connected to the vertex
  inline size_t numWall() const;
  //!Number of dimensions for the position vector
  inline size_t numPosition() const;

//!Returns a reference to the cell vector
  inline const std::vector<Cell*> & cell() const;
//!Returns the cell pointer from position k
  inline Cell* cell( size_t k );
//!Returns a reference to the wall vector
  inline const std::vector<Wall*> & wall() const;
//!Returns the wall pointer from position k
  inline Wall* wall( size_t k ) const;
  
//!The position vector
  inline const std::vector<double> & position() const;
//!The position in dimension d
  inline double position(size_t d);

//!Sets the index variable
  inline void setIndex( size_t value );
//!Sets a cell pointer
  inline void setCell( size_t index,Cell* val );
//!Sets pointers to cells
  inline void setCell( std::vector<Cell*> &val );
//!Adds a cell to the vertex
  inline void addCell( Cell* val );
//!Removes a cell from the vertex
  int removeCell( Cell* val );
//!Sets a wall pointer
  inline void setWall( size_t index,Wall* val );
//!Sets pointers to walls
  inline void setWall( std::vector<Wall*> &val );
//!Adds a wall to the vertex
  inline void addWall( Wall* val );
//!Removes a wall from the vertex
  int removeWall( Wall* val );
//!Sets the position from a vector
  inline void setPosition(std::vector<double> &pos);
	///
	/// @brief Check if the vertex is at the boundary of the tissue.
	///
	/// A boundary vertex is defined from whether any of the connected walls is
	/// connected to the background.
	///
	/// @para Cell* background is a pointer to the tissue background.
	///
	/// @return 1 if boundary and 0 otherwise.
	int isBoundary(Cell *background) const;
};

inline size_t Vertex::index() const { return index_; }
inline std::string Vertex::id() const { return id_; }
inline size_t Vertex::numCell() const { return cell_.size(); }
inline size_t Vertex::numWall() const { return wall_.size(); }
inline size_t Vertex::numPosition() const { return position_.size(); }
inline const std::vector<Cell*> & Vertex::cell() const { return cell_; }
inline Cell* Vertex::cell( size_t k ) { return cell_[k]; }
inline const std::vector<Wall*> & Vertex::wall() const { return wall_; }
inline Wall* Vertex::wall( size_t k ) const { return wall_[k]; }
inline const std::vector<double> & Vertex::position() const { return position_; }
inline double Vertex::position(size_t d) { return position_[d]; }
inline void Vertex::setIndex( size_t value ) { index_ = value; }
inline void Vertex::setCell( size_t index,Cell* val ) { cell_[index]=val; }
inline void Vertex::setCell( std::vector<Cell*> &val ) { cell_=val; }
inline void Vertex::addCell( Cell* val ) { cell_.push_back(val); }
inline void Vertex::setWall( size_t index,Wall* val ) { wall_[index]=val; }
inline void Vertex::setWall( std::vector<Wall*> &val ) { wall_=val; }
inline void Vertex::addWall( Wall* val ) { wall_.push_back(val); }
inline void Vertex::setPosition(std::vector<double> &pos) { position_=pos; }

#endif
