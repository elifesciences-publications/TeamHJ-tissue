/**
 * Filename     : tissue.h
 * Description  : A class describing a two-dimensional tissue of cells
 * Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
 * Created      : April 2006
 * Revision     : $Id:$
 */
#ifndef TISSUE_H
#define TISSUE_H

#include <assert.h>
#include <fstream>
#include <iostream>
#include <list>
#include <string>
#include <vector>
#include "baseReaction.h"
#include "baseCompartmentChange.h"
#include "cell.h"
#include "direction.h"
#include "vertex.h"
#include "wall.h"

//!Describes the properties of a two-dimensional cell tissue
/*!*/ 
class Tissue {
  
 private:
  
  std::string id_;           
  
  std::vector<Cell> cell_;
  std::vector<Wall> wall_;
  std::vector<Vertex> vertex_;
  Cell background_;
	Direction direction_;
  std::vector<BaseReaction*> reaction_;
  std::vector<BaseCompartmentChange*> compartmentChange_;
	//std::vector< std::vector<double> > tmpCellData_;

	std::vector<size_t> directionalWall_;
	
 public:
  
  Tissue();
  Tissue( const Tissue & tissueCopy );
  Tissue( const std::vector<Cell> &cellVal,
	  const std::vector<Wall> &wallVal,
	  const std::vector<Vertex> &vertexVal );
  Tissue( const char *initFile, int verbose=0 );
  Tissue( std::string initFile, int verbose=0 );
  
  ~Tissue();

  void readInit(const char *initFile,int verbose=0);
  void readInit(std::string initFile,int verbose=0);
  void readInit(std::ifstream &IN,int verbose=0);
  void readMerryInit(const char *initFile,int verbose=0);

  void readModel(const char *modelFile,int verbose=0);
  void readModel(std::string modelFile,int verbose=0);
  void readModel(std::ifstream &IN,int verbose=0);
  
  void createTissueFromSpheres(std::vector< std::vector<double> > &y,
			       double rFac=1.0, int verbose=0);
  void createTissueFromVoronoi(std::vector< std::vector<double> > &vertexPos,
			       std::vector< std::vector<size_t> > &cellVertex,
			       int verbose=0);
  
  
  inline std::string id() const;
  //inline double volume() const;
  //inline size_t mitosisFlag() const;
  
  inline size_t numCell() const;
  inline size_t numWall() const;
  inline size_t numVertex() const;
  inline size_t numReaction() const;
  inline size_t numCompartmentChange() const;
	inline size_t numDirectionalWall() const;

	inline size_t directionalWall(size_t i) const;
  inline const std::vector<Cell> & cell() const;
  inline const Cell & cell(size_t i) const;
  inline Cell & cell(size_t i);
	inline Cell * cellP(size_t i);
  inline void addCell( Cell val );
  inline void removeCell( size_t index );
  void removeEpidermalCells(std::vector< std::vector<double> > &cellData,
			    std::vector< std::vector<double> > &wallData,
			    std::vector< std::vector<double> > &vertexData,
			    std::vector< std::vector<double> > &cellDeriv,
			    std::vector< std::vector<double> > &wallDeriv,
			    std::vector< std::vector<double> > &vertexDeriv,
			    double radialThreshold=0.0);
  void removeEpidermalCellsAtDistance(std::vector< std::vector<double> > &cellData,
				      std::vector< std::vector<double> > &wallData,
				      std::vector< std::vector<double> > &vertexData,
				      std::vector< std::vector<double> > &cellDeriv,
				      std::vector< std::vector<double> > &wallDeriv,
				      std::vector< std::vector<double> > &vertexDeriv,
				      double radialThreshold,double max,
				      size_t direction);
  inline Cell* background();
	inline Direction* direction();
  inline const std::vector<Wall> & wall() const;
  inline const Wall & wall(size_t i) const;
  inline Wall & wall(size_t i);
  inline void addWall( Wall val );
  inline void removeWall( size_t index );
  inline const std::vector<Vertex> & vertex() const;
  inline const Vertex & vertex(size_t i) const;
  inline Vertex & vertex(size_t i);
	inline size_t numDimension();
  inline void addVertex( Vertex val );
  inline void removeVertex( size_t index );
  inline BaseReaction* reaction(size_t i) const;
  inline BaseCompartmentChange* compartmentChange(size_t i) const;
  
  inline void setId(std::string value);
  inline void setNumCell(size_t val);
  inline void setNumWall(size_t val);
  inline void setNumVertex(size_t val);
  inline void setNumDirectionalWall(size_t val);

  //inline const std::vector< std::vector<double> > & tmpCellData() const;
  //inline void setTmpCellData(std::vector< std::vector<double> > &val);
  
  inline void setWallLengthFromVertexPosition();
	inline void setDirectionalWall(size_t i,size_t val);
  int addReaction( std::istream &IN );
  int addCompartmentChange( std::istream &IN );
  
  //Functions related to model simulations
  void derivs( std::vector< std::vector<double> > &cellData,
							 std::vector< std::vector<double> > &wallData,
							 std::vector< std::vector<double> > &vertexData,
							 std::vector< std::vector<double> > &cellDeriv,
							 std::vector< std::vector<double> > &wallDeriv,
							 std::vector< std::vector<double> > &vertexDeriv );
	void initiateReactions(std::vector< std::vector<double> > &cellData,
												 std::vector< std::vector<double> > &wallData,
												 std::vector< std::vector<double> > &vertexData);
	void updateReactions(double step);
	void initiateDirection(std::vector< std::vector<double> > &cellData,
												 std::vector< std::vector<double> > &wallData,
												 std::vector< std::vector<double> > &vertexData,
												 std::vector< std::vector<double> > &cellDerivs,
												 std::vector< std::vector<double> > &wallDerivs,
												 std::vector< std::vector<double> > &vertexDerivs );
	void updateDirection(double step,							
											 std::vector< std::vector<double> > &cellData,
											 std::vector< std::vector<double> > &wallData,
											 std::vector< std::vector<double> > &vertexData,
											 std::vector< std::vector<double> > &cellDerivs,
											 std::vector< std::vector<double> > &wallDerivs,
											 std::vector< std::vector<double> > &vertexDerivs );
	void updateDirectionDivision(size_t cellI,
															 std::vector< std::vector<double> > &cellData,
															 std::vector< std::vector<double> > &wallData,
															 std::vector< std::vector<double> > &vertexData,
															 std::vector< std::vector<double> > &cellDerivs,
															 std::vector< std::vector<double> > &wallDerivs,
															 std::vector< std::vector<double> > &vertexDerivs);
  void checkCompartmentChange(std::vector< std::vector<double> > &cellData,
															std::vector< std::vector<double> > &wallData,
															std::vector< std::vector<double> > &vertexData,
															std::vector< std::vector<double> > &cellDeriv,
															std::vector< std::vector<double> > &wallDeriv,
															std::vector< std::vector<double> > &vertexDeriv );
  //!Updates topology and variables for a cell removal
  void removeCell(size_t cellIndex,
									std::vector< std::vector<double> > &cellData,
									std::vector< std::vector<double> > &wallData,
									std::vector< std::vector<double> > &vertexData,
									std::vector< std::vector<double> > &cellDeriv,
									std::vector< std::vector<double> > &wallDeriv,
									std::vector< std::vector<double> > &vertexDeriv );			
  //! Updates topology and variables for cell division
  void divideCell( Cell *divCell, size_t w1, size_t w2, 
									 std::vector<double> &v1Pos,
									 std::vector<double> &v2Pos,
									 std::vector< std::vector<double> > &cellData,
									 std::vector< std::vector<double> > &wallData,
									 std::vector< std::vector<double> > &vertexData,
									 std::vector< std::vector<double> > &cellDeriv,
									 std::vector< std::vector<double> > &wallDeriv,
									 std::vector< std::vector<double> > &vertexDeriv,
									 std::vector<size_t> &volumeChangeList,
									 double threshold=0.0);
  
  //!Sorts cell.wall and cell.vertex vectors to be cyclic 
  void sortCellWallAndCellVertex(Cell* cell=NULL);
	void sortCellRecursive( Cell* cell, std::vector<size_t> &sortedFlag, size_t &numSorted);
  //!Checks all connectivities for errors
  void checkConnectivity(size_t verbose=0);
  
  ///Finds maxima in a variable column for the cells 
  unsigned int findPeaksGradientAscent( std::vector< std::vector<double> > &cellData, 
																				size_t col, std::vector<size_t> &cellMax,
																				std::vector<size_t> &flag );
  
  //Print functions
  void printInit(std::ostream &os=std::cout);
  void printInit(std::vector< std::vector<double> > &cellData,
								 std::vector< std::vector<double> > &wallData,
								 std::vector< std::vector<double> > &vertexData,
								 std::ostream &os);
  void printVertex(std::ostream &os=std::cout);
  void printWall(std::ostream &os=std::cout);
  void printVertexAndCell(std::ostream &os=std::cout);
  void printVertexAndCell(std::vector< std::vector<double> > &cellData,
													std::vector< std::vector<double> > &vertexData,
													std::ostream &os=std::cout);  
  void printVertexAndWall(std::vector< std::vector<double> > &wallData,
													std::vector< std::vector<double> > &vertexData,
													std::ostream &os=std::cout);  
};

//!The tissue id
inline std::string Tissue::id() const { return id_; }

//!Number of cells for the tissue
inline size_t Tissue::numCell() const { return cell_.size(); }

//!Number of walls for the tissue
inline size_t Tissue::numWall() const { return wall_.size(); }

//!Number of vertices for the tissue
inline size_t Tissue::numVertex() const { return vertex_.size(); }

//!Number of reactions in the tissue model
inline size_t Tissue::numReaction() const { return reaction_.size(); }

//!Number of compartment changes defined in the tissue model
inline size_t Tissue::numCompartmentChange() const 
{ return compartmentChange_.size(); }

//!Number of directional walls
inline size_t Tissue::numDirectionalWall() const 
{ 
	return directionalWall_.size(); 
}

//!Returns the directional wall for cell i
inline size_t Tissue::directionalWall(size_t i) const 
{ 
	return directionalWall_[i];
}

//!Returns a reference to the cell vector
inline const std::vector<Cell> & Tissue::cell() const { return cell_; }

//!Returns a const reference to cell i 
inline const Cell & Tissue::cell(size_t i) const { return cell_[i]; }

//!Returns a reference to cell i 
inline Cell & Tissue::cell(size_t i) { return cell_[i]; }

///
/// @brief Return a pointer to cell i
//
inline Cell* Tissue::cellP(size_t i) {return &cell_[i];}

//!Adds a cell to the vector
inline void Tissue::addCell( Cell val ) { cell_.push_back(val);}

//!Removes cell at index and moves last cell into that position 
inline void Tissue::removeCell( size_t index ) {
	assert(index<numCell());
	Cell *cp = &cell_[numCell()-1];
	Cell *cpNew = &cell_[index];
	if( cp != cpNew ) {
		cell_[index] = cell_[numCell()-1];
		cell_[index].setIndex(index);
		//Reconnect wall neighbors
		for( size_t k=0 ; k<cell_[index].numWall() ; ++k )
			if( cell_[index].wall(k)->cell1() == cp )
				cell_[index].wall(k)->setCell1(cpNew);
			else if( cell_[index].wall(k)->cell2() == cp )
				cell_[index].wall(k)->setCell2(cpNew);
		//Reconnect vertex neighbors
		for( size_t k=0 ; k<cell_[index].numVertex() ; ++k )
			for( size_t l=0 ; l<cell_[index].vertex(k)->numCell() ; ++l )
				if( cell_[index].vertex(k)->cell(l) == cp )
					cell_[index].vertex(k)->setCell(l,cpNew);
	}
	cell_.pop_back();
}

//!Returns a pointer to the background
inline Cell* Tissue::background() { return &background_;}

//!Returns a pointer to the direction
inline Direction* Tissue::direction() { return &direction_;}

//!Returns a reference to the wall vector
inline const std::vector<Wall> & Tissue::wall() const { return wall_; }

//!Returns a const reference to wall i 
inline const Wall & Tissue::wall(size_t i) const { return wall_[i]; }

//!Returns a reference to wall i 
inline Wall & Tissue::wall(size_t i) { return wall_[i]; }

//!Adds a wall to the vector
inline void Tissue::addWall( Wall val ) { wall_.push_back(val);}

//!Removes wall at index and moves last wall into that position 
inline void Tissue::removeWall( size_t index ) {
	assert(index<numWall());
	Wall *wp = &wall_[numWall()-1];
	Wall *wpNew = &wall_[index];
	if( wp != wpNew ) {
		wall_[index] = wall_[numWall()-1];
		wall_[index].setIndex(index);
		//Reconnect cell neighbors
		for( size_t l=0 ; l<wall_[index].cell1()->numWall() ; ++l )
			if( wall_[index].cell1()->wall(l) == wp )
				wall_[index].cell1()->setWall(l,wpNew);
		for( size_t l=0 ; l<wall_[index].cell2()->numWall() ; ++l )
			if( wall_[index].cell2()->wall(l) == wp )
				wall_[index].cell2()->setWall(l,wpNew);
		//Reconnect vertex neighbors
		for( size_t l=0 ; l<wall_[index].vertex1()->numWall() ; ++l )
			if( wall_[index].vertex1()->wall(l) == wp )
				wall_[index].vertex1()->setWall(l,wpNew);
		for( size_t l=0 ; l<wall_[index].vertex2()->numWall() ; ++l )
			if( wall_[index].vertex2()->wall(l) == wp )
				wall_[index].vertex2()->setWall(l,wpNew);	
	}
	wall_.pop_back();
}

//!Returns a reference to the vertex vector
inline const std::vector<Vertex> & Tissue::vertex() const { return vertex_; }

//!Returns a const reference to vertex i 
inline const Vertex & Tissue::vertex(size_t i) const { return vertex_[i]; }

//!Returns a reference to vertex i 
inline Vertex & Tissue::vertex(size_t i) { return vertex_[i]; }

inline size_t Tissue::numDimension() 
{
	return vertex(0).numPosition();
}

//!Adds a vertex to the vector
inline void Tissue::addVertex( Vertex val ) { vertex_.push_back(val);}

//!Removes vertex at index and moves last vertex into that position 
inline void Tissue::removeVertex( size_t index ) {
	assert(index<numVertex());
	Vertex *vp = &vertex_[numVertex()-1];
	Vertex *vpNew = &vertex_[index];
	std::cerr << "Vertex " << vpNew->index() << " removed and "
						<< vp->index() << " moved" << std::endl;
	if( vp != vpNew ) {
		vertex_[index] = vertex_[numVertex()-1];
		vertex_[index].setIndex(index);
		std::cerr << "In cell (" << vertex_[index].numCell() << ") ";
		//Reconnect cell neighbors
		for( size_t k=0 ; k<vertex_[index].numCell() ; ++k ) {
			size_t moved=0;
			for( size_t l=0 ; l<vertex_[index].cell(k)->numVertex() ; ++l )
				if( vertex_[index].cell(k)->vertex(l) == vp ) {
					vertex_[index].cell(k)->setVertex(l,vpNew);
					std::cerr << vertex_[index].cell(k)->index() << " ";
					++moved;
				}
			if( !moved )
				std::cerr << vertex_[index].cell(k)->index() << "(not) ";
		}
		std::cerr << std::endl;
		//Reconnect wall neighbors
		for( size_t k=0 ; k<vertex_[index].numWall() ; ++k )
			if( vertex_[index].wall(k)->vertex1() == vp )
				vertex_[index].wall(k)->setVertex1(vpNew);
			else if( vertex_[index].wall(k)->vertex2() == vp )
				vertex_[index].wall(k)->setVertex2(vpNew);
	}
	vertex_.pop_back();
}

//!Returns a pointer to reaction i
inline BaseReaction* Tissue::reaction(size_t i) const { 
  assert( i<reaction_.size() );
  return reaction_[i];
}

//!Returns a pointer to compartmentChange i
inline BaseCompartmentChange* Tissue::compartmentChange(size_t i) const { 
  assert( i<compartmentChange_.size() );
  return compartmentChange_[i];
}

//!Sets the id string
inline void Tissue::setId(std::string value) {
  id_ = value;
}

//!Resizes the cell_ vector
inline void Tissue::setNumCell(size_t val) {
  cell_.resize(val);
}

//!Resizes the wall_ vector
inline void Tissue::setNumWall(size_t val) {
  wall_.resize(val);
}

//!Resizes the vertex_ vector
inline void Tissue::setNumVertex(size_t val) {
  vertex_.resize(val);
}

//!Resizes the directionalWall_ vector
inline void Tissue::setNumDirectionalWall(size_t val) 
{
  directionalWall_.resize(val);
}

//!Sets a directional wall for cell i
inline void Tissue::setDirectionalWall(size_t i,size_t val) 
{
  directionalWall_[i]=val;
}

//inline const std::vector< std::vector<double> > & Tissue::
//tmpCellData() const 
//{
//	return tmpCellData_;
//}

//inline void Tissue::
//setTmpCellData(std::vector< std::vector<double> > &val)
//{
//	tmpCellData_ = val;
//}


#endif
