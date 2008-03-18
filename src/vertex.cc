/**
 * Filename     : vertex.cc
 * Description  : A class describing a vertex
 * Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
 * Created      : April 2006
 * Revision     : $Id:$
 */

#include"vertex.h"
#include"wall.h"
  
Vertex::Vertex() 
{
}

Vertex::Vertex( const Vertex & vertexCopy ) 
{
  index_ = vertexCopy.index();
  id_ = vertexCopy.id();
  position_ = vertexCopy.position();
  cell_ = vertexCopy.cell();
  wall_ = vertexCopy.wall();
}
  
Vertex::~Vertex() 
{
}

int Vertex::removeCell( Cell* val ) 
{
	for (size_t k=0; k<cell_.size(); ++k)
		if (cell_[k]==val) {
			cell_[k]=cell_[cell_.size()-1];
			cell_.pop_back();
			return 1;
		}
	return 0;
}

int Vertex::removeWall( Wall* val ) 
{
	for (size_t k=0; k<wall_.size(); ++k)
		if (wall_[k]==val) {
			wall_[k]=wall_[wall_.size()-1];
			wall_.pop_back();
			return 1;
		}
	return 0;
}

int Vertex::isBoundary(Cell *background) const
{
	for (size_t wI=0; wI<numWall(); ++wI)
		if (wall_[wI]->hasCell(background))
			return 1;
	return 0;
}

