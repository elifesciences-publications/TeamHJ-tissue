//
// Filename     : vertex.h
// Description  : A class describing a vertex
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : April 2006
// Revision     : $Id$
//
#ifndef VERTEX_H
#define VERTEX_H

#include<vector>
#include<list>
#include<string>
#include<iostream>
#include<fstream>

class Cell;
class Wall;

///
/// @brief Defines the properties of a (1D) vertex element, including updatable spatial position
///
class Vertex {
  
 private:
  
  size_t index_;
  std::string id_;           
  
  std::vector<Cell*> cell_;
  std::vector<Wall*> wall_;
  
  std::vector<double> position_;
  
  std::vector<double> stressDirection_;
  
 public:
  
  Vertex();
  Vertex( const Vertex & cellCopy );
  
  ~Vertex();
  
  ///
  /// @brief Returns the vertex index
  ///
  /// The index may be used when checking connectivity. In normal situations, the indices
  /// are automatically generated when creating a Tissue.
  ///
  inline size_t index() const;
  ///
  /// @brief The vertex id
  ///
  /// The id/name is optional and is not used in any update functions.
  ///
  inline std::string id() const;
  ///
  /// @brief Number of cells connected to the vertex
  ///
  inline size_t numCell() const;
  ///
  /// @ brief Number of walls connected to the vertex
  ///
  inline size_t numWall() const;
  ///
  /// @brief Number of dimensions for the position vector
  ///
  /// For most usage, the number of dimensions should be 2 or 3. Tissue uses this to
  /// find out its dimension.
  ///
  /// @see Tissue::dimension()
  /// @see position()
  ///
  inline size_t numPosition() const;
  ///
  /// @brief Returns a reference to the cell vector
  ///
  /// A vertex is connected to a number of Cells (at least one) stored in a vector. 
  /// There are no formal restrictions on the size of this vector.
  ///
  /// @see Cell
  ///
  inline const std::vector<Cell*> & cell() const;
  ///
  /// @brief Returns the cell pointer from position k in the Cell vector for the vertex.
  ///
  /// @see cell()
  ///
  inline Cell* cell( size_t k );
  ///
  /// @brief Returns a reference to the wall vector
  ///
  /// A vertex is connected to a number of Walls (at least two) stored in a vector. 
  /// There are no formal restrictions on the size of this vector.
  /// 
  /// @see Wall
  ///
  inline const std::vector<Wall*> & wall() const;
  ///
  /// @brief Returns the wall pointer from position k in the Wall vector for the Vertex.
  ///
  inline Wall* wall( size_t k ) const;
  ///  
  /// @brief Returns the position vector
  ///
  /// The position vector is the main variables of the Vertex, and is avector of dimension 2 or 3.
  ///
  inline const std::vector<double> & position() const;
  ///
  /// @brief The position in dimension d
  ///
  /// @see position()
  ///
  inline double position(size_t d);
  
  ///
  /// @brief Sets the index variable
  ///
  /// @see index()
  ///
  inline void setIndex( size_t value );
  ///
  /// @brief Sets a cell pointer
  ///
  /// @see cell()
  ///
  inline void setCell( size_t index,Cell* val );
  ///
  /// @brief Sets pointers to cells
  ///
  /// Sets a new Cell vector for the Vertex (and removes any old information).
  ///
  /// @see cell()
  ///
  inline void setCell( std::vector<Cell*> &val );
  ///
  /// @brief Adds a Cell to the end of the vector for the Vertex.
  ///
  /// @see cell()
  ///
  inline void addCell( Cell* val );
  ///
  /// @brief Removes a Cell from the vector for the Vertex.
  ///
  /// It looks for the provided Cell and removes it from the cell vector if it is present.
  /// If the cell is removed, it returns 1, otherwise 0. Note that it does not update the 
  /// connection from the Cell to the Vertex.
  ///
  /// @see cell()
  ///
  int removeCell( Cell* val );
  ///
  /// @brief Sets a wall pointer at position index of the wall vector for the Vertex.
  ///
  /// @see wall()
  ///
  inline void setWall( size_t index,Wall* val );
  ///
  /// @brief Sets pointers to walls
  ///
  /// Sets a new wall vector for the vertex, and remove any previous walls.
  ///
  /// @see wall()
  ///
  inline void setWall( std::vector<Wall*> &val );
  ///
  /// @brief Adds a Wall to the vector for the Vertex.
  ///
  /// @see wall()
  ///
  inline void addWall( Wall* val );
  ///
  /// @brief Removes a Wall from the vector for the Vertex
  ///
  /// Removes the Wall provided if it is in the vector for the Vertex. It returns 1 if the wall
  /// is removed, and 0 otherwise. Note that it does not update the connection from the Wall to
  /// the vertex.
  ///
  /// @see wall()
  ///
  int removeWall( Wall* val );
  ///
  /// @brief Sets the position from a vector
  ///
  /// @see position()
  ///
  inline void setPosition(std::vector<double> &pos);
  ///
  /// @brief Sets a vertex position in dimension d.
  ///
  inline void setPosition(size_t d, double pos);
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
  
  ///
  /// @brief Calculate stress direction for the vertex from stresses in the connected walls.
  ///
  /// A stress direction is calculated for the vertex by combining stresses in the walls
  /// connected to the vertex. The calculation is based on Dumains and Kwiatkowska,
  /// Plan Journal (2002) !HJ: have to be checked if this ref is correct!
  ///
  void calculateStressDirection(std::vector< std::vector<double> > &vertexData, 
				std::vector< std::vector<double> > &wallData, 
				std::vector<size_t> &wallForceIndexes);  
  ///
  /// @brief Returns a stress direction vector (stored) for the vertex.
  ///
  /// Note: the value in the stored vector may be obselete if the tissue has been updated 
  /// since the last calculation.
  ///
  /// @see calculateStressDirection(std::vector< std::vector<double> >&,std::vector< std::vector<double> >&, std::vector<size_t>&); 
  std::vector<double> stressDirection(void) const;
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
inline void Vertex::setPosition(size_t d,double pos) { position_[d]=pos; }

#endif
