//
// Filename     : wall.cc
// Description  : A class describing a one-dimensional wall (for 2d tissue)
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : April 2006
// Revision     : $Id$
//

#include"wall.h"
#include<cmath>  


Wall::Wall() {

  cellSort_.first = 0;
  cellSort_.second = 0;
}

Wall::Wall( const Wall& wallCopy ) {

  index_ = wallCopy.index();
  id_ = wallCopy.id();
  length_ = wallCopy.length();
  cell_.first = wallCopy.cell1();
  cell_.second = wallCopy.cell2();
  cellSort_.first = wallCopy.cellSort1();
  cellSort_.second = wallCopy.cellSort2();
  vertex_.first = wallCopy.vertex1();
  vertex_.second = wallCopy.vertex2();
  variable_.resize(wallCopy.numVariable());
  for( size_t i=0 ; i<variable_.size() ; ++i )
    variable_[i] = wallCopy.variable(i);
}
  
Wall::~Wall() {
}

double Wall::setLengthFromVertexPosition() {
  size_t dimension = vertex1()->numPosition();
  double distance=0.0;
  for( size_t d=0 ; d<dimension ; ++d )
    distance += ( vertex1()->position(d) - vertex2()->position(d) ) *
      ( vertex1()->position(d) - vertex2()->position(d) );
  setLength( std::sqrt(distance) );
  return length();
}

double Wall::
setLengthFromVertexPosition( DataMatrix &vertexData) 
{
  size_t dimension = vertex1()->numPosition();
  size_t v1I=vertex1()->index();
  size_t v2I=vertex2()->index();
  double distance=0.0;
  for( size_t d=0 ; d<dimension ; ++d )
    distance += ( vertexData[v1I][d]-vertexData[v2I][d] ) *
      ( vertexData[v1I][d]-vertexData[v2I][d] );
  setLength( std::sqrt(distance) );
  return length();
}

double Wall::
lengthFromVertexPosition()
{
  size_t dimension = vertex1()->numPosition();
  double distance=0.0;
  for( size_t d=0 ; d<dimension ; ++d )
    distance += ( vertex1()->position(d)-vertex2()->position(d) ) *
      ( vertex1()->position(d)-vertex2()->position(d) );
  return std::sqrt(distance);
}

double Wall::
lengthFromVertexPosition( DataMatrix &vertexData)
{
  size_t dimension = vertex1()->numPosition();
  size_t v1I=vertex1()->index();
  size_t v2I=vertex2()->index();
  double distance=0.0;
  for( size_t d=0 ; d<dimension ; ++d )
    distance += ( vertexData[v1I][d]-vertexData[v2I][d] ) *
      ( vertexData[v1I][d]-vertexData[v2I][d] );
  return std::sqrt(distance);
}
