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

#include"myTypedefs.h"
#include"wall.h"

//class Wall;
class Vertex;
class Tissue;

///
/// @brief Describes the properties of a two-dimensional cell within a tissue
///
/// A cell includes all information for faces (2D cells) of the tissues,
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
  DataMatrix E_;
  
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
  ///
  /// @brief Returns the cell index
  ///
  /// The index is ofimportance since it is used in connection queries. Normally, the index
  /// is set automatically when a tissue is created.
  ///
  inline size_t index() const;
  ///
  /// @brief returns the cell id (name)
  ///
  /// A name is optional and is not used at any important instances.
  ///
  inline std::string id() const;
  ///
  /// @brief The cell volume
  ///
  /// Returns the volume variable from the cell that should hold the cell volume(area).
  /// Note that if the tissue has changed (i.e. the vertex positions have moved) this value
  /// can be obsolete. If not sure use calculateVolume(DataMatrix&,size_t)
  /// to calculate it directly.
  ///
  /// @see calculateVolume(DataMatrix&,size_t)
  ///
  inline double volume() const;
  ///
  /// @brief A flag that is set if the cell can divide 
  ///
  /// Its value should be 1 or 0.
  ///
  inline size_t mitosisFlag() const;
  ///
  /// @brief Number of walls for the cell	
  ///
  inline size_t numWall() const;
  ///
  /// @brief Number of vertices for the cell
  ///
  inline size_t numVertex() const;
  ///
  /// @brief Number of variables in the cell
  ///
  inline size_t numVariable() const;
  ///
  /// @brief Returns a reference to the wall vector
  ///
  inline const std::vector<Wall*> & wall() const;
  ///
  /// @brief Returns a wall pointer from position k in the vector
  ///
  inline Wall* wall( size_t k );
  ///
  /// @brief Returns a reference to the wall at position k in the vector
  ///
  inline Wall& wallRef( size_t k );
  ///
  /// @brief Returns 1 if the cell has the wall
  ///
  inline int hasWall( Wall *val );
  ///
  /// @brief Returns a reference to the vertex vector
  ///
  inline const std::vector<Vertex*> & vertex() const;
  ///
  /// @brief Returns a vertex pointer from position v in the vector
  ///
  inline Vertex* vertex( size_t v );
  ///
  /// @brief Returns true if the cell has the vertex
  ///
  inline int hasVertex( Vertex *val );
  ///
  /// @brief Returns a reference to the variable vector
  ///
  inline const std::vector<double> & variable() const;
  ///
  /// @brief Returns variable i from the variable vector
  ///
  inline double variable(size_t i);
  ///
  /// @brief Sets the id (name) string
  ///
  inline void setId(std::string value);
  ///
  /// @brief Sets the index of a cell
  ///  
  inline void setIndex(size_t value);
  ///
  /// @brief Sets the wall vector of a cell
  ///  
  inline void setWall( std::vector<Wall*> &val );
  ///
  /// @brief Sets a wall at index in the wall vector
  ///  
  inline void setWall( size_t index,Wall* val);
  ///
  /// @brief Adds a wall to the end of the vector in a cell
  ///  
  inline void addWall( Wall* val );
  ///
  /// @brief Sets the wall vector of a cell
  ///  
  inline void setVertex( std::vector<Vertex*> &val );
  ///
  /// @brief Sets a vertex at index in the vertex vector
  ///  
  inline void setVertex( size_t index,Vertex* val);
  ///
  /// @brief Adds a vertex to the end of the vector in a cell
  ///  
  inline void addVertex( Vertex* val );
  ///
  /// @brief Sets a variable at index in the variable vector
  ///  
  inline void setVariable( size_t index,double val );
  ///
  /// @brief Sets all variables in the variable vector to the values provided
  ///  
  void setVariable( std::vector<double> &val );
  ///
  /// @brief Adds a variable to the end of the vector in a cell
  ///  
  inline void addVariable( double val );
  ///
  /// @brief Checks if any of the walls connects to the given cell
  ///
  /// Returns 1 if the Cell 'neighbor' is a neighbor to the Cell, and 0 otherwise.
  /// A cell neighbor is defined by having the same wall.
  ///
  inline int isNeighbor( Cell *neighbor );
  ///
  /// @brief Returns a pointer to the cell neighbor on the other side of wall k
  ///
  inline Cell* cellNeighbor(size_t k);
  
  ///
  /// @brief Sorts the vertices and walls connected to the cell such that they are cyclic
  ///
  /// Obsolete version that sort cells independently of each other and hence might lead to 
  /// cell normal vectors in opposite direction.
  ///
  /// @see sortWallAndVertex(Tissue&);
  ///
  void sortWallAndVertexOld(Tissue &T);
  ///
  /// @brief Sorts the vertices and walls connected to the cell such that they are cyclic
  ///
  /// Uses also information (if available) such that cells in the tissue are sorted 'together'
  /// to generate a consistent cell normal direction (the normals point in 'the same' direction
  /// for all cells).
  ///
  void sortWallAndVertex(Tissue &T);
  
  ///
  /// @brief Returns TRUE if the cell has any 'concave' wall pairs
  ///
  bool isConcave(DataMatrix &vertexData, const double tolerance = 0.01);
  ///
  /// @brief Returns TRUE if any of the cell walls cross each other
  ///
  bool isFolded(DataMatrix &vertexData);
  ///
  /// @brief Calculates the cell volume(area) from the vertex positions.
  ///
  /// Uses a cross-product rule of the wall directions to calculate the cell
  /// area. Works only if the cell vertices are sorted/cyclic. The vertex
  /// positions used are taken from the vertices directly (and hence might be obselete
  /// if the vertex positions have been moved.
  ///
  /// @see calculateVolume(DataMatrix&,size_t)
  /// 
  double calculateVolume( size_t signFlag=0 );
  ///
  /// @brief Calculates the cell volume(area) from the vertex positions.
  ///
  /// Uses a cross-product rule of the wall directions to calculate the cell
  /// area. In two dimensions the calculation is:
  ///
  /// @f[ A = \frac{1}{2} | \sum_i^{vertex} (x_i y_{i+1} - y_i x_{i+1}) | @f]
  /// 
  /// In three dimensions, the rule is:
  ///
  /// @f[ A = \frac{1}{2} \sqrt{\sum_i^{vertex} (r_i \times r_{i+1})^2} @f]
  ///
  /// where @f$r_i @f$ is the vector from the cell center to vertex i.
  ///  
  /// Works only if the cell vertices are sorted/cyclic. The vertex
  /// positions used are taken from the provided matrix. The volume is returned
  /// and also stored in a volume variable in the Cell.
  ///
  /// Note: The 3D version adds the area of triangles meshed with the
  /// center and maybe should be replaced with the cross-product rule
  /// on the PCA-extracted plane. (the 2D version).
  /// 
  /// @see volume()
  ///
  double calculateVolume( DataMatrix 
			  &vertexData, size_t signFlag=0 );
  ///
  /// @brief Calculates the cell volume(area) from the vertex positions using triangles.
  ///
  /// Assumes center 'vertex' defined and uses triangles to calculate the area (2D volume).
  ///
  /// @see volume()
  ///
  double calculateVolumeCenterTriangulation( DataMatrix 
					     &vertexData,
					     DataMatrix &cellData,
					     size_t centerIndex);
  ///
  /// @brief Calculates the cell center-of-mass position.
  ///
  /// Uses a cross-product rule to define the center of mass of the cell. Works
  /// only if the cell vertices are sorted/cyclic. The vertex
  /// positions used are taken from the vertices directly. The volume is returned
  /// and also stored in a volume variable in the Cell.
  ///
  /// @see volume()
  ///
  std::vector<double> positionFromVertex();
  ///
  /// @brief Calculates the cell center-of-mass position.
  ///
  /// Uses a cross-product rule to define the center of mass of the cell. Works
  /// only if the cell vertices are sorted/cyclic. The vertex
  /// positions used are taken from the provided matrix.
  ///
  std::vector<double> positionFromVertex( DataMatrix &vertexData );	
  
  class FailedToFindRandomPositionInCellException
  {  
  };
  
  std::vector<double> randomPositionInCell(const DataMatrix &vertexData, const int numberOfTries = 10000);
  
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
  void calculatePCAPlane(DataMatrix &vertexData);
  ///
  /// @brief Returns the (two) vectors spanning the PCA plane defined by the cell vertex positions
  ///
  /// Returns E that stores the PCA plane (two first principal directions) 
  /// calculated from the cell vertex poitions (in 3D). It relies upon that
  /// calculatePCAPlane has been called since last movement of the vertices.
  ///
  /// @see calculatePCAPlane()
  ///
  DataMatrix getPCAPlane(void) const;
  ///
  /// @brief Returns the vertex positions on the PCA plane defined by the cell vertex positions
  ///
  /// Returns the vertex positions (2D) as projected onto the PCA plane 
  /// calculated from the cell vertex poitions (in 3D). It relies upon that
  /// calculatePCAPlane has been called since last movement of the vertices.
  ///
  /// @see calculatePCAPlane()
  ///
  std::vector< std::pair<double, double> > projectVerticesOnPCAPlane(DataMatrix &vertexData);
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
			 DataMatrix &vertexData);
  ///
  /// @brief Returns the normal to the cell plane for triangular cells
  ///
  /// Returns the normal vector of the cell plane 
  /// calculated from the cell vertex poitions (in 3D). It relies on that
  /// the vertices are sorted (cyclic) for the cell.
  ///
  std::vector<double> getNormalTriangular(DataMatrix &vertexData);
};

inline size_t Cell::index() const 
{ 
  return index_; 
}

inline std::string Cell::id() const 
{ 
  return id_; 
}

inline double Cell::volume() const 
{ 
  return volume_; 
}

inline size_t Cell::mitosisFlag() const 
{ 
  return mitosisFlag_; 
}

inline size_t Cell::numWall() const 
{ 
  return wall_.size(); 
}

inline size_t Cell::numVertex() const 
{ 
  return vertex_.size(); 
}

inline size_t Cell::numVariable() const 
{ 
  return variable_.size(); 
}

inline const std::vector<Wall*> & Cell::wall() const 
{ 
  return wall_; 
}

inline Wall* Cell::wall( size_t k ) 
{ 
  return wall_[k]; 
}

inline Wall& Cell::wallRef( size_t k ) 
{ 
  return *wall_[k]; 
}

inline int Cell::hasWall( Wall *val ) 
{
  std::vector<Wall*>::iterator it = find(wall_.begin(),wall_.end(),val);
  if( it != wall_.end() ) return 1;
  return 0;
}

inline const std::vector<Vertex*> & Cell::vertex() const 
{ 
  return vertex_; 
}

inline Vertex* Cell::vertex( size_t v ) 
{ 
  return vertex_[v]; 
}

inline int Cell::hasVertex( Vertex *val ) 
{
  std::vector<Vertex*>::iterator it = find(vertex_.begin(),vertex_.end(),val);
  if( it != vertex_.end() ) return 1;
  return 0;
}

inline const std::vector<double> & Cell::variable() const 
{ 
  return variable_; 
}

inline double Cell::variable(size_t i) 
{ 
  return variable_[i];
}

inline void Cell::setId(std::string value) 
{
  id_ = value;
}

inline void Cell::setIndex(size_t value) 
{
  index_ = value;
}

inline void Cell::setWall( std::vector<Wall*> &val ) 
{
  wall_=val;
}

inline void Cell::setWall( size_t index,Wall* val) 
{
  assert( index<wall_.size() );
  wall_[index]=val;
}

inline void Cell::addWall( Wall* val ) 
{
  wall_.push_back(val);
}

inline void Cell::setVertex( std::vector<Vertex*> &val ) 
{
  vertex_=val;
}

inline void Cell::setVertex( size_t index,Vertex* val) 
{
  assert( index<vertex_.size() );
  vertex_[index]=val;
}

inline void Cell::addVertex( Vertex* val ) 
{
  vertex_.push_back(val);
}

inline void Cell::setVariable( size_t index,double val ) 
{
  variable_[index]=val;
}

inline void Cell::addVariable( double val ) 
{
  variable_.push_back(val);
}

inline int Cell::isNeighbor( Cell *neighbor ) 
{
  for (size_t i=0; i<numWall(); ++i )
    if (wall(i)->cell1()==neighbor || wall(i)->cell2()==neighbor) 
      return 1;
  return 0;
}

inline Cell* Cell::cellNeighbor(size_t k)
{
  return wall(k)->cell1()==this ? wall(k)->cell2() : wall(k)->cell1(); 
}

#endif
