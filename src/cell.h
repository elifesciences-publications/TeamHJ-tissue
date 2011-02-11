/**
 * Filename     : cell.h
 * Description  : A class describing a two-dimensional cell
 * Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
 * Created      : April 2006
 * Revision     : $Id:$
 */
#ifndef CELL_H
#define CELL_H

#include<algorithm>
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

///
/// @brief Describes the properties of a two-dimensional cell within a tissue
///
/// A cell includes all information for the faces (2D cells) of the tissues,
/// including its connectivity to walls and vertices. It can store a number
/// of variables representing e.g. molecular concentrations, and have functions
/// for calculating properties such as volume (area).
///
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
  std::vector< std::vector<double> > E_;
  
 public:
  
	///
	/// @brief The empty cell constructor.
	///
	/// sets only the mitosisFlag to 0.
	///
  Cell();
	///
	/// @brief Copy constructor
	///
	/// Copies all relevant variables from the provided cell.
	///
  Cell( const Cell & cellCopy );
	///
	/// @brief Constructor setting index and id/name for the cell.
	///
  Cell(size_t indexVal,std::string idVal);
	///
	/// @brief Destructor (empty in the current implementation).
	///
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
  
	void sortWallAndVertexOld(Tissue &T);
	void sortWallAndVertex(Tissue &T);

	bool isConcave(std::vector< std::vector<double> > &vertexData, const double tolerance = 0.01);
	bool isFolded(std::vector< std::vector<double> > &vertexData);
	///
	/// @brief Calculates the cell volume(area) from the vertex positions.
	///
	/// Uses a cross-product rule of the wall directions to calculate the cell
	/// area. Works only if the cell vertices are sorted/cyclic. The vertex
	/// positions used are taken from the vertices directly.
	/// 
  double calculateVolume( size_t signFlag=0 );
	///
	/// @brief Calculates the cell volume(area) from the vertex positions.
	///
	/// Uses a cross-product rule of the wall directions to calculate the cell
	/// area. Works only if the cell vertices are sorted/cyclic. The vertex
	/// positions used are taken from the provided matrix.
	/// 
  double calculateVolume( std::vector< std::vector<double> > 
													&vertexData, size_t signFlag=0 );
	///
	/// @brief Calculates the cell center-of-mass position.
	///
	/// Uses a cross-product rule to define the center of mass of the cell. Works
	/// only if the cell vertices are sorted/cyclic. The vertex
	/// positions used are taken from the vertices directly.
	///
  std::vector<double> positionFromVertex();
  std::vector<double> 
	///
	/// @brief Calculates the cell center-of-mass position.
	///
	/// Uses a cross-product rule to define the center of mass of the cell. Works
	/// only if the cell vertices are sorted/cyclic. The vertex
	/// positions used are taken from the provided matrix.
	///
	positionFromVertex( std::vector< std::vector<double> > &vertexData );	

	class FailedToFindRandomPositionInCellException
	{
		
	};

	std::vector<double> randomPositionInCell(const std::vector< std::vector<double> > &vertexData, const int numberOfTries = 10000);



  // These functions handles projection to a PCAPlane. Important: Make
  // sure you call calculatePCAPlane() before any of the other
  // functions.
	///
	/// @brief Calculates the PCA plane (two principal directions) from cell vertex positions
	///
	/// This functions creates a PCA plane from the cell vertex positions, hence
	/// creating a 2D representation of the cell used for 3D simulations. The function
	/// stores the two 'largest' principal directions in the cell E_ variable.
	/// The PCA plane is used to apply some 2D rules in 3D simulations. The usage
	/// of the PCA plane assumes that the vertices of a cell are close to a plane,
	/// although the calculations do not require it (but otherwise estimates of)
	/// e.g. cell volume(area) would have errors).
 	/// 
  void calculatePCAPlane(std::vector< std::vector<double> > &vertexData);
	///
	/// @brief Returns the (two) vectors spanning the PCA plane defined by the cell vertex positions
	///
	/// Returns E that stores the PCA plane (two first principal directions) 
	/// calculated from the cell vertex poitions (in 3D). It relies upon that
	/// calculatePCAPlane has been called since last movement of the vertices.
	///
	/// @see calculatePCAPlane()
	///
  std::vector< std::vector<double> > getPCAPlane(void) const;
	///
	/// @brief Returns the vertex positions on the PCA plane defined by the cell vertex positions
	///
	/// Returns the vertex positions (2D) as projected onto the PCA plane 
	/// calculated from the cell vertex poitions (in 3D). It relies upon that
	/// calculatePCAPlane has been called since last movement of the vertices.
	///
	/// @see calculatePCAPlane()
	///
  std::vector< std::pair<double, double> > projectVerticesOnPCAPlane(std::vector< std::vector<double> > &vertexData);
	///
	/// @brief Returns the normal to the PCA plane defined by the cell vertex positions
	///
	/// Returns the normal vector of the PCA plane 
	/// calculated from the cell vertex poitions (in 3D). It relies upon that
	/// calculatePCAPlane has been called since last movement of the vertices.
	///
	/// @see calculatePCAPlane()
	///
  std::vector<double> getNormalToPCAPlane(void);
	int vectorSignFromSort(std::vector<double> &n,
												 std::vector< std::vector<double> > &vertexData);
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
