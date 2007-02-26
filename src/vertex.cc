/**
 * Filename     : vertex.cc
 * Description  : A class describing a vertex
 * Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
 * Created      : April 2006
 * Revision     : $Id:$
 */

#include"vertex.h"
  
Vertex::Vertex() {

}

Vertex::Vertex( const Vertex & vertexCopy ) {

  index_ = vertexCopy.index();
  id_ = vertexCopy.id();
  position_ = vertexCopy.position();
  cell_ = vertexCopy.cell();
  wall_ = vertexCopy.wall();

}
  
Vertex::~Vertex() {

}
