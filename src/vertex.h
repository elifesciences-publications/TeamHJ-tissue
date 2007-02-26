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
  
  inline size_t index() const;
  inline std::string id() const;

  inline size_t numCell() const;
  inline size_t numWall() const;
  inline size_t numPosition() const;

  inline const std::vector<Cell*> & cell() const;
  inline Cell* cell( size_t k );
  inline const std::vector<Wall*> & wall() const;
  inline Wall* wall( size_t k );
  
  inline const std::vector<double> & position() const;
  inline double position(size_t d);

  inline void setIndex( size_t value );
  inline void setCell( size_t index,Cell* val );
  inline void setCell( std::vector<Cell*> &val );
  inline void addCell( Cell* val );
  inline int removeCell( Cell* val );
  inline void setWall( size_t index,Wall* val );
  inline void setWall( std::vector<Wall*> &val );
  inline void addWall( Wall* val );
  inline int removeWall( Wall* val );
  inline void setPosition(std::vector<double> &pos);

};

//!Returns the cell index
inline size_t Vertex::index() const { return index_; }

//!The cell id
inline std::string Vertex::id() const { return id_; }

//!Number of cells connected to the vertex
inline size_t Vertex::numCell() const { return cell_.size(); }

//!Number of walls connected to the vertex
inline size_t Vertex::numWall() const { return wall_.size(); }

//!Number of dimensions for the position vector
inline size_t Vertex::numPosition() const { return position_.size(); }

//!Returns a reference to the cell vector
inline const std::vector<Cell*> & Vertex::cell() const { return cell_; }

//!Returns the cell pointer from position k
inline Cell* Vertex::cell( size_t k ) { return cell_[k]; }

//!Returns a reference to the wall vector
inline const std::vector<Wall*> & Vertex::wall() const { return wall_; }

//!Returns the wall pointer from position k
inline Wall* Vertex::wall( size_t k ) { return wall_[k]; }

//!The position vector
inline const std::vector<double> & Vertex::position() const { return position_; }

//!The position in dimension d
inline double Vertex::position(size_t d) { return position_[d]; }

//!Sets the index variable
inline void Vertex::setIndex( size_t value ) { index_ = value; }

//!Sets a cell pointer
inline void Vertex::setCell( size_t index,Cell* val ) {
  cell_[index]=val;
}

//!Sets pointers to cells
inline void Vertex::setCell( std::vector<Cell*> &val ) {
  cell_=val;
}

//!Adds a cell to the vertex
inline void Vertex::addCell( Cell* val ) {
  cell_.push_back(val);
}

//!Removes a cell from the vertex
inline int Vertex::removeCell( Cell* val ) {
	for( size_t k=0 ; k<cell_.size() ; ++k )
		if( cell_[k] == val ) {
			cell_[k]=cell_[cell_.size()-1];
			cell_.pop_back();
			return 1;
		}
	return 0;
}

//!Sets a wall pointer
inline void Vertex::setWall( size_t index,Wall* val ) {
  wall_[index]=val;
}

//!Sets pointers to walls
inline void Vertex::setWall( std::vector<Wall*> &val ) {
  wall_=val;
}

//!Adds a wall to the vertex
inline void Vertex::addWall( Wall* val ) {
  wall_.push_back(val);
}

//!Removes a wall from the vertex
inline int Vertex::removeWall( Wall* val ) {
	for( size_t k=0 ; k<wall_.size() ; ++k )
		if( wall_[k] == val ) {
			wall_[k]=wall_[wall_.size()-1];
			wall_.pop_back();
			return 1;
		}
	return 0;
}

//!Sets the position from a vector
inline void Vertex::setPosition(std::vector<double> &pos) {
  position_=pos;
}
#endif
