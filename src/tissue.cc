//
// Filename     : tissue.cc
// Description  : A class describing a two-dimensional cell tissue
// Author(s)    : Henrik Jonsson (henrik at thep.lu.se)
// Created      : April 2006
// Revision     : $Id:$
//
#include <algorithm>
#include <assert.h>
#include <cmath>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>
#include "tissue.h"
#include "wall.h"
#include "myFiles.h"
#include "myMath.h"
//#include "ply_reader.h"

Tissue::Tissue() {  
  cell_.reserve(200000);
  wall_.reserve(200000);
  vertex_.reserve(200000);
  Cell tmpCell(static_cast<size_t>(-1),static_cast<std::string>("Background"));
  background_ = tmpCell;
}

Tissue::Tissue( const Tissue & tissueCopy ) {
  cell_.reserve(200000);
  wall_.reserve(200000);
  vertex_.reserve(200000);
  Cell tmpCell(static_cast<size_t>(-1),static_cast<std::string>("Background"));
  background_ = tmpCell;
}

Tissue::Tissue( const std::vector<Cell> &cellVal,
		const std::vector<Wall> &wallVal,
		const std::vector<Vertex> &vertexVal ) {
  
  cell_.reserve(200000);
  wall_.reserve(200000);
  vertex_.reserve(200000);

  Cell tmpCell(static_cast<size_t>(-1),static_cast<std::string>("Background"));
  background_ = tmpCell;
  cell_ = cellVal;
  wall_ = wallVal;
  vertex_ = vertexVal;
}

Tissue::Tissue( const char *initFile, int verbose ) {
  cell_.reserve(200000);
  wall_.reserve(200000);
  vertex_.reserve(200000);

  Cell tmpCell(static_cast<size_t>(-1),static_cast<std::string>("Background"));
  background_ = tmpCell;
  readInit(initFile,verbose);
}


Tissue::Tissue( std::string initFile, int verbose ) 
{
  cell_.reserve(200000);
  wall_.reserve(200000);
  vertex_.reserve(200000);
	
  Cell tmpCell(static_cast<size_t>(-1),static_cast<std::string>("Background"));
  background_ = tmpCell;
  readInit(initFile,verbose);
}

Tissue::Tissue( DataMatrix &cellData,
								DataMatrix &wallData,
								DataMatrix &vertexData,
								std::vector< std::vector<size_t> > &cellVertex,
								std::vector< std::vector<size_t> > &wallVertex,
								int verbose)
{
  cell_.reserve(200000);
  wall_.reserve(200000);
  vertex_.reserve(200000);
	
  Cell tmpCell(static_cast<size_t>(-1),static_cast<std::string>("Background"));
  background_ = tmpCell;
	
	size_t numCell = cellData.size();
	size_t numWall = wallData.size();
	size_t numVertex = vertexData.size();
  setNumCell( numCell );
  setNumWall( numWall );
  setNumVertex( numVertex );
	assert( numCell );
  assert( numWall );
  assert( numVertex );
  
  //Set all indices to the placement in the vectors
  for( size_t i=0 ; i<numCell ; ++i )
    cell(i).setIndex(i);
  for( size_t i=0 ; i<numWall ; ++i )
    wall(i).setIndex(i);
  for( size_t i=0 ; i<numVertex ; ++i )
    vertex(i).setIndex(i);
  //
  // set cell variables
  //
  size_t numCellVar = cellData[0].size();
  if( numCellVar ) {
    for( size_t i=0 ; i<numCell ; ++i ) {
      if (cellData[i].size() != numCellVar) {
        std::cerr << "Tissue::Tissue(cellData,wallData,vertexData,cellVertex,wallVertex) "
                  << "Wrong number of variables in cell " << i << std::endl;
        exit(-1);
      }
      for( size_t j=0 ; j<numCellVar ; ++j ) {
        cell(i).addVariable(cellData[i][j]);
      }
    }
  }
  //
  // Set wall data
  //
  if (wallData[0].size() < 1) {
    std::cerr << "Tissue::Tissue(cellData,wallData,vertexData,cellVertex,wallVertex) "
              << " At least wall length must be given in wall variables." << std::endl; 
    exit(-1);
  }
  size_t numWallVar=wallData.size();
  for (size_t i = 0; i < numWall; ++i) {
    if (wallData[i].size() != numWallVar) {
      std::cerr << "Tissue::Tissue(cellData,wallData,vertexData,cellVertex,wallVertex) "
                << "Wrong number of variables in wall " << i << std::endl;
      exit(-1);
    }
    wall(i).setLength(wallData[i][0]);
    for (size_t j=1; j<numWallVar; ++j) {
      wall(i).addVariable(wallData[i][j]);
    }
  }
  //
  //Set vertex positions
  //
  size_t dimension = vertexData[0].size();
  assert( dimension==2 || dimension==3 );
  for( size_t i=0 ; i<numVertex ; ++i ) {
    if (vertexData[i].size() != dimension) {
      std::cerr << "Tissue::Tissue(cellData,wallData,vertexData,cellVertex,wallVertex) "
                << "Wrong dimension in vertex " << i << std::endl;
      exit(-1);
    }
    vertex(i).setPosition(vertexData[i]);
  }
  //
  // Set connectivity
  //
  // Cell-Vertex
  //
  if (numCell!=cellVertex.size()) {
    std::cerr << "Tissue::Tissue(cellData,wallData,vertexData,cellVertex,wallVertex) "
              << "Cell number in cell variables not same as in cell vertex." << std::endl;
    exit(-1);
  }
  for (size_t i=0; i<numCell; ++i) {
    size_t numCellVertex=cellVertex.size();
    for (size_t k=0; k<numCellVertex; ++k) {
      size_t j=cellVertex[i][k];
      cell(i).addVertex( vertexP(j) );
      vertex(j).addCell( cellP(i) );
    }
  }
  //
  // Wall-Vertex
  //
  if (numWall!=wallVertex.size()) {
    std::cerr << "Tissue::Tissue(cellData,wallData,vertexData,cellVertex,wallVertex) "
              << "Wall number in wall variables not same as in wall vertex." << std::endl;
    exit(-1);
  }
  for (size_t i=0; i<numWall; ++i) {
    assert (wallVertex.size()==2);
    size_t j1=wallVertex[i][0];
    size_t j2=wallVertex[i][1];
    vertex(j1).addWall( wallP(i) );
    vertex(j2).addWall( wallP(i) );
    wall(i).setVertex( vertexP(j1),vertexP(j2) );
  }
  //
  // Cell-Wall
  //
  // For each wall, find one or two cells that have the same vertex pair. 
  // If two cells connect wall with these, and if one connect wall with this and background.
  //
  for (size_t wI=0; wI<numWall; ++wI) {
    size_t cellCount=0;
    std::vector<size_t> cI(2,numCell);
    size_t vI1=wallVertex[wI][0];
    size_t vI2=wallVertex[wI][1];
    for (size_t i=0; i<numCell; ++i) {
      size_t cellVertexCount=0;
      for (size_t k=0; k<cellVertex[i].size(); ++k) {
        if (cellVertex[i][k]==vI1 || cellVertex[i][k]==vI2) {
          ++cellVertexCount;
        }
      }
      if (cellVertexCount==2) {//cell and wall connected
        if (cellCount>1) {
          std::cerr << "Tissue::Tissue(cellData,wallData,vertexData,cellVertex,wallVertex) "
                    << " More than two cells found to wall " << wI << std::endl;
          exit(-1);
        }
        cell(i).addWall(wallP(wI));
        cI[cellCount] = i;
        ++cellCount;
      }
    }
    if (cellCount==1) {
      wall(wI).setCell( cellP(cI[0]),background() );
    }
    else if (cellCount==2) {
      wall(wI).setCell( cellP(cI[0]),cellP(cI[1]) );
    }
    else {
      std::cerr << "Tissue::Tissue(cellData,wallData,vertexData,cellVertex,wallVertex) "
                << " Wrong number of cells found to wall " << wI << std::endl;
      exit(-1);
    }		
  }
  
  //Sort all cellWalls and cellVertices to comply with area calculations
  //and plotting
  sortCellWallAndCellVertex();
  checkConnectivity(verbose);
}


Tissue::~Tissue() {
}

void Tissue::setWallLengthFromVertexPosition() {
  for( size_t i=0 ; i<numWall() ; ++i )
    wall(i).setLengthFromVertexPosition();
}

int Tissue::addReaction( std::istream &IN ) {
  if( !IN ) return -1;
  reaction_.push_back( BaseReaction::createReaction(IN) );
  return 0;
}

int Tissue::addCompartmentChange( std::istream &IN ) {
  if( !IN ) return -1;
  compartmentChange_.push_back( BaseCompartmentChange::createCompartmentChange(IN) );
  return 0;
}

void Tissue::addSisterVertex(size_t val1,size_t val2) 
{
  std::vector<size_t> tmpIndex(2);
  tmpIndex[0]=val1;
  tmpIndex[1]=val2;
  sisterVertexIndex_.push_back(tmpIndex);
}

void Tissue::readInit(std::istream &IN,int verbose) {
  
  unsigned int numCellVal,numWallVal,numVertexVal;
  //std::string idVal;
  
  //IN >> idVal;
  IN >> numCellVal;
  IN >> numWallVal;
  IN >> numVertexVal;

  //setId( idVal );
  setNumCell( numCellVal );
  setNumWall( numWallVal );
  setNumVertex( numVertexVal );

  assert( numCellVal );
  assert( numWallVal );
  assert( numVertexVal );
  
  //Set all indeces to the placement in the vectors
  for( size_t i=0 ; i<numCellVal ; ++i )
    cell(i).setIndex(i);
  for( size_t i=0 ; i<numWallVal ; ++i )
    wall(i).setIndex(i);
  for( size_t i=0 ; i<numVertexVal ; ++i )
    vertex(i).setIndex(i);
  
  //Read connective topology
  //
  if( verbose )
    std::cerr << "Tissue::readInit(IN) - reading connectivity topology" << std::endl;
  int wI,c1I,c2I,v1I,v2I;
  size_t w,c1,c2,v1,v2;
  for( size_t nW=0 ; nW<numWall() ; ++nW ) {
    //Read the connections
    IN >> wI;
    IN >> c1I;
    IN >> c2I;
    IN >> v1I;
    IN >> v2I;
    w = static_cast<size_t>(wI);
    c1 = static_cast<size_t>(c1I);
    c2 = static_cast<size_t>(c2I);
    v1 = static_cast<size_t>(v1I);
    v2 = static_cast<size_t>(v2I);
    //Assert all data is ok
    if( verbose>1 ) {
      std::cerr << wI << " " << c1I << " " << c2I << " " << v1I 
		<< " " << v2I << std::endl;    
      std::cerr << w << " " << c1 << " " << c2 << " " << v1 << " " << v2 
		<< " " << static_cast<size_t>(-1) << std::endl 
		<< std::endl;
    }    
    assert( w==nW );
    assert( c1==static_cast<size_t>(-1) || c1<numCell() );
    assert( c2==static_cast<size_t>(-1) || c2<numCell() );
    assert( v1<numVertex() );
    assert( v2<numVertex() );
    //Make the connections by first extracting the pointers
    Wall *wp = &(wall(w));
    Vertex *v1p = &(vertex(v1)),*v2p=&(vertex(v2));
    Cell *c1p,*c2p;
    if( c1 != static_cast<size_t>(-1) )
      c1p = &(cell(c1));
    else
      c1p = &(background_);
    if( c2 != static_cast<size_t>(-1) )
      c2p = &(cell(c2));
    else
      c2p = &(background_);
    //vertex-wall
    wall(w).setVertex( v1p,v2p );
    vertex(v1).addWall( wp );
    vertex(v2).addWall( wp );
    //cell-wall
    wall(w).setCell(c1p,c2p);
    if( c1 != static_cast<size_t>(-1) )
      cell(c1).addWall( wp );
    if( c2 != static_cast<size_t>(-1) )
      cell(c2).addWall( wp );
    //cell-vertex
    if( c1 != static_cast<size_t>(-1) ) {
      if( !cell(c1).hasVertex(v1p) ) {
				cell(c1).addVertex( v1p );
				vertex(v1).addCell( c1p );
      }
      if( !cell(c1).hasVertex(v2p) ) {
				cell(c1).addVertex( v2p );
				vertex(v2).addCell( c1p );
      }
    }
    if( c2 != static_cast<size_t>(-1) ) {
      if( !cell(c2).hasVertex(v1p) ) {
				cell(c2).addVertex( v1p );
				vertex(v1).addCell( c2p );
      }
      if( !cell(c2).hasVertex(v2p) ) {
				cell(c2).addVertex( v2p );
				vertex(v2).addCell( c2p );
      }
    }
  }
  //Read vertex positions
  //
  if( verbose )
    std::cerr << "Tissue::readInit(IN) - reading vertex positions" 
	      << std::endl;
  size_t numVertexTmp,dimension;
  IN >> numVertexTmp;
  IN >> dimension;
  assert( numVertexTmp==numVertex() );
  assert( dimension==2 || dimension==3 );
  if (verbose) {
    std::cerr << numVertexTmp << "(" << numVertex() << ") vertices in " << dimension << " dimensions." << std::endl;
  }

  std::vector<double> pos(dimension);
  for( size_t i=0 ; i<numVertex() ; ++i ) {
    for( size_t j=0 ; j<dimension ; ++j )
      IN >> pos[j];
    vertex(i).setPosition(pos);
  }
	
  //Read wall data
  //
  if( verbose )
    std::cerr << "Tissue::readInit(IN) - reading wall data" << std::endl;
	
  size_t numWallTmp,numLength,numVar;
  IN >> numWallTmp;
  IN >> numLength;
  IN >> numVar;
  assert( numWallTmp==numWall() );
  assert( numLength==1 );
  if (verbose) {
    std::cerr << numWallTmp << "(" << numWall() << ") walls." << std::endl;
  }
  double length;
  double value;
  for (size_t i = 0; i < numWall(); ++i) {
	  IN >> length;
	  wall(i).setLength(length);
	  std::vector<double> variable;
	  for (size_t j = 0; j < numVar; ++j) {
	    IN >> value;
	    variable.push_back(value);
	  }
	  wall(i).setVariable(variable);
  }

  // Read cell variables
  //
  if( verbose )
    std::cerr << "Tissue::readInit(IN) - reading cell variables" << std::endl;
  size_t numCellTmp,numCellVar;
  IN >> numCellTmp;
  IN >> numCellVar;
  assert( numCellTmp==numCell() );
  double var=0.0;
  if( numCellVar )
    for( size_t i=0 ; i<numCell() ; ++i ) {
      for( size_t j=0 ; j<numCellVar ; ++j ) {
	IN >> var;
	cell(i).addVariable(var);
      }
    }
  //Sort all cellWalls and cellVertices to comply with area calculations
  //and plotting
  sortCellWallAndCellVertex();
  checkConnectivity(verbose);
}

void Tissue::readInit( const char *initFile, int verbose ) {

  std::istream *IN = myFiles::openFile(initFile);
  //std::ifstream IN(initFile);
  if( !IN ) {
    std::cerr << "Tissue::readInit(char*) - "
	      << "Cannot open file " << initFile << std::endl; exit(-1);}
  if( verbose )
    std::cerr << "Tissue::readInit(char*) - calling readInit(IN)" << std::endl;
  readInit((std::istream &) *IN,verbose);
}

void Tissue::readInit( std::string initFile, int verbose ) {

  const char* iFile = initFile.c_str();
  std::istream *IN = myFiles::openFile(iFile);
  //std::ifstream IN(iFile);
  if( !IN ) {
    std::cerr << "Tissue::readInit(string) - "
	      << "Cannot open file " << initFile << std::endl; exit(-1);}
  if( verbose )
    std::cerr << "Tissue::readInit(string) - calling readInit(IN)" << std::endl;
  readInit((std::istream &) *IN,verbose);
}

void Tissue::readInitCenterTri(const char *initFile,int verbose) {
  
  std::ifstream IN(initFile);
  if( !IN ) {
    std::cerr << "Tissue::readInitCenterTri(char*,int) - "
	      << "Cannot open file " << initFile 
	      << std::endl; 
    exit(EXIT_FAILURE);
  }

  unsigned int numCellVal,numWallVal,numVertexVal;
  //std::string idVal;
  
  //IN >> idVal;
  IN >> numCellVal;
  IN >> numWallVal;
  IN >> numVertexVal;

  //setId( idVal );
  setNumCell( numCellVal );
  setNumWall( numWallVal );
  setNumVertex( numVertexVal );

  assert( numCellVal );
  assert( numWallVal );
  assert( numVertexVal );
  
  //Set all indeces to the placement in the vectors
  for( size_t i=0 ; i<numCellVal ; ++i )
    cell(i).setIndex(i);
  for( size_t i=0 ; i<numWallVal ; ++i )
    wall(i).setIndex(i);
  for( size_t i=0 ; i<numVertexVal ; ++i )
    vertex(i).setIndex(i);
  
  //Read connective topology
  //
  if( verbose )
    std::cerr << "Tissue::readInit(IN) - reading connectivity topology" << std::endl;
  int wI,c1I,c2I,v1I,v2I;
  size_t w,c1,c2,v1,v2;
  for( size_t nW=0 ; nW<numWall() ; ++nW ) {
    //Read the connections
    IN >> wI;
    IN >> c1I;
    IN >> c2I;
    IN >> v1I;
    IN >> v2I;
    w = static_cast<size_t>(wI);
    c1 = static_cast<size_t>(c1I);
    c2 = static_cast<size_t>(c2I);
    v1 = static_cast<size_t>(v1I);
    v2 = static_cast<size_t>(v2I);
    //Assert all data is ok
    if( verbose>1 ) {
      std::cerr << wI << " " << c1I << " " << c2I << " " << v1I 
		<< " " << v2I << std::endl;    
      std::cerr << w << " " << c1 << " " << c2 << " " << v1 << " " << v2 
		<< " " << static_cast<size_t>(-1) << std::endl 
		<< std::endl;
    }    
    assert( w==nW );
    assert( c1==static_cast<size_t>(-1) || c1<numCell() );
    assert( c2==static_cast<size_t>(-1) || c2<numCell() );
    assert( v1<numVertex() );
    assert( v2<numVertex() );
    //Make the connections by first extracting the pointers
    Wall *wp = &(wall(w));
    Vertex *v1p = &(vertex(v1)),*v2p=&(vertex(v2));
    Cell *c1p,*c2p;
    if( c1 != static_cast<size_t>(-1) )
      c1p = &(cell(c1));
    else
      c1p = &(background_);
    if( c2 != static_cast<size_t>(-1) )
      c2p = &(cell(c2));
    else
      c2p = &(background_);
    //vertex-wall
    wall(w).setVertex( v1p,v2p );
    vertex(v1).addWall( wp );
    vertex(v2).addWall( wp );
    //cell-wall
    wall(w).setCell(c1p,c2p);
    if( c1 != static_cast<size_t>(-1) )
      cell(c1).addWall( wp );
    if( c2 != static_cast<size_t>(-1) )
      cell(c2).addWall( wp );
    //cell-vertex
    if( c1 != static_cast<size_t>(-1) ) {
      if( !cell(c1).hasVertex(v1p) ) {
        cell(c1).addVertex( v1p );
        vertex(v1).addCell( c1p );
      }
      if( !cell(c1).hasVertex(v2p) ) {
        cell(c1).addVertex( v2p );
        vertex(v2).addCell( c1p );
      }
    }
    if( c2 != static_cast<size_t>(-1) ) {
      if( !cell(c2).hasVertex(v1p) ) {
        cell(c2).addVertex( v1p );
        vertex(v1).addCell( c2p );
      }
      if( !cell(c2).hasVertex(v2p) ) {
        cell(c2).addVertex( v2p );
        vertex(v2).addCell( c2p );
      }
    }
  }
  //Read vertex positions
  //
  if( verbose )
    std::cerr << "Tissue::readInit(IN) - reading vertex positions" 
	      << std::endl;
  size_t numVertexTmp,dimension;
  IN >> numVertexTmp;
  IN >> dimension;
  assert( numVertexTmp==numVertex() );
  assert( dimension==2 || dimension==3 );
  if (verbose) {
    std::cerr << numVertexTmp << "(" << numVertex() << ") vertices in " << dimension << " dimensions." << std::endl;
  }

  std::vector<double> pos(dimension);
  for( size_t i=0 ; i<numVertex() ; ++i ) {
    for( size_t j=0 ; j<dimension ; ++j )
      IN >> pos[j];
    vertex(i).setPosition(pos);
  }
  
  //Read wall data
  //
  if( verbose )
    std::cerr << "Tissue::readInit(IN) - reading wall data" << std::endl;
  
  size_t numWallTmp,numLength,numVar;
  IN >> numWallTmp;
  IN >> numLength;
  IN >> numVar;
  assert( numWallTmp==numWall() );
  assert( numLength==1 );
  if (verbose) {
    std::cerr << numWallTmp << "(" << numWall() << ") walls." << std::endl;
  }
  double length;
  double value;
  for (size_t i = 0; i < numWall(); ++i) {
    IN >> length;
    wall(i).setLength(length);
    std::vector<double> variable;
    for (size_t j = 0; j < numVar; ++j) {
      IN >> value;
      variable.push_back(value);
    }
    wall(i).setVariable(variable);
  }
  
  // Read cell variables
  //
  if( verbose )
    std::cerr << "Tissue::readInit(IN) - reading cell variables" << std::endl;
  size_t numCellTmp,numCellVar;
  IN >> numCellTmp;
  IN >> numCellVar;
  assert( numCellTmp==numCell() );
  double var=0.0;
  if (verbose) {
    std::cerr << numCellTmp << "(" << numCell() << ") cells with " << numCellVar 
	      << "variables plus center position and internal edges." << std::endl;
  }
  if( numCellVar )
    for( size_t i=0 ; i<numCell() ; ++i ) {
      // read cell variables
      for( size_t j=0 ; j<numCellVar ; ++j ) {
	IN >> var;
	cell(i).addVariable(var);
      }
      // read center positions for the cell
      if (dimension!=3) {
	std::cerr << "Tissue::readInitCenterTri() Center triangulation implemented for three"
		  << " dimensions only." << std::endl;
	exit(EXIT_FAILURE);
      }
      for (size_t d=0; d<dimension; ++d) {
	IN >> var;
	cell(i).addCenterPosition(var);
      }
      size_t numInternalEdge = cell(i).numVertex();
      for (size_t j=0; j<numInternalEdge; ++j) {
	IN >> var;
	cell(i).addEdgeLength(var);
      }
    }
  //Sort all cellWalls and cellVertices to comply with area calculations
  //and plotting
  sortCellWallAndCellVertex();
  checkConnectivity(verbose);
}

void Tissue::readInitMerryProj( const char *initFile, int verbose ) 
{
  std::ifstream IN(initFile);
  if( !IN ) {
    std::cerr << "Tissue::readInitMerryProj(char*) - "
	      << "Cannot open file " << initFile 
	      << std::endl; 
    exit(EXIT_FAILURE);
  }
  unsigned int numVertexVal,dimension;
  IN >> numVertexVal;
  IN >> dimension;
  setNumVertex( numVertexVal );
  for( size_t i=0 ; i<numVertexVal ; ++i )
    vertex(i).setIndex(i);
  
  std::vector<double> pos(dimension);
  std::vector<size_t> cellName,vertexName;
  // Read information about vertices
  for( size_t i=0 ; i<numVertexVal ; ++i ) {
    size_t tmp,numVertexCell;
    IN >> tmp;
    vertexName.push_back(tmp);
    for( size_t dim=0 ; dim<dimension ; ++dim )
      IN >> pos[dim];
    vertex(i).setPosition(pos);
    IN >> numVertexCell;
    for( size_t j=0 ; j<numVertexCell ; ++j ) {
      size_t tmpCellIndex,cellIndex=numCell(),newCellFlag=1;
      IN >> tmpCellIndex;
      for( size_t c=0 ; c<cellName.size() ; ++c ) {
	if( tmpCellIndex==cellName[c] ) {
	  cellIndex=c;
	  newCellFlag=0;
	  break;
	}
      }
      if( newCellFlag ) {
	Cell tmpCell(cellIndex,"");
	cellName.push_back(tmpCellIndex);
	addCell(tmpCell);
      }
      vertex(i).addCell(&(cell(cellIndex)));
      cell(cellIndex).addVertex(&(vertex(i)));
    }
  }
  // Reading the walls as well
  size_t numWallVal;
  IN >> numWallVal;
  // Store vertexName,vertexIndex pairs in a hash table
  std::map<size_t,size_t> vertexNameToIndex;
  for (size_t i=0; i<numVertex(); ++i)
    vertexNameToIndex[vertexName[i]] = i;
  
  for (size_t i=0; i<numWallVal; ++i) {
    Wall tmpWall;
    tmpWall.setIndex(i);
    size_t v1Name,v2Name;
    IN >> v1Name;
    IN >> v2Name;
    if (vertexNameToIndex.find(v1Name) == vertexNameToIndex.end() ||
	vertexNameToIndex.find(v2Name) == vertexNameToIndex.end() ) {
      std::cerr << "Tissue::readInitMerryProj() Vertex read from file (in wall list) not found."
		<< std::endl;
      exit(-1);
    }
    Vertex *tmpVertex1(vertexP(vertexNameToIndex[v1Name]));
    Vertex *tmpVertex2(vertexP(vertexNameToIndex[v2Name]));
    tmpWall.setVertex1(tmpVertex1);
    tmpWall.setVertex2(tmpVertex2);
    // Find neighboring cells as well
    size_t numFoundCell=0;
    for (size_t c1I=0; c1I<tmpVertex1->numCell(); ++c1I)
      for (size_t c2I=0; c2I<tmpVertex2->numCell(); ++c2I)
	if (tmpVertex1->cell(c1I)==tmpVertex2->cell(c2I)) {
	  if (numFoundCell==0) {
	    tmpWall.setCell1(tmpVertex1->cell(c1I));
	  }
	  else if (numFoundCell==1) {
	    tmpWall.setCell2(tmpVertex1->cell(c1I));
	  }
	  ++numFoundCell;
	}
    if (numFoundCell==1) {
      tmpWall.setCell2(background());
      ++numFoundCell;
    }
    if (numFoundCell!=2) {
      std::cerr << "Tissue::readInitMerryProj() Found " << numFoundCell
		<< " cells for wall " << i << std::endl;
      std::cerr << "Vertices: " << tmpVertex1->index() << " (" << v1Name << ") " 
		<< tmpVertex2->index() << " (" << v2Name << ")" << std::endl;
      for (size_t c1I=0; c1I<tmpVertex1->numCell(); ++c1I)
	for (size_t c2I=0; c2I<tmpVertex2->numCell(); ++c2I)
	  std::cerr << tmpVertex1->cell(c1I)->index() << " " << tmpVertex2->cell(c2I)->index()
		    << std::endl;
      exit(-1);
    } 
    // Add the wall to tissue (including pointers from cells and vertices)
    addWall(tmpWall);
    Wall *tmpWallP=wallP(numWall()-1);
    tmpVertex1->addWall(tmpWallP);
    tmpVertex2->addWall(tmpWallP);
    if (tmpWallP->cell1()!=background())
      tmpWallP->cell1()->addWall(tmpWallP);
    if (tmpWallP->cell2()!=background())
      tmpWallP->cell2()->addWall(tmpWallP);
  }
  IN.close();
  
  
  if( verbose>1 ) {
    std::cerr << "Vertices:" << std::endl;
    for (size_t i=0; i<numVertex(); ++i) {
      std::cerr << vertex(i).index() << "\t";
      for (size_t dim=0; dim<vertex(i).numPosition(); ++dim)
	std::cerr << vertex(i).position(dim) << " ";
      std::cerr << "\t";
      for (size_t k=0; k<vertex(i).numCell(); ++k)
	std::cerr << vertex(i).cell(k)->index() << " ";
      std::cerr << std::endl;
    }
    
    std::cerr << "Cells:" << std::endl;
    for (size_t i=0; i<numCell(); ++i) {
      std::cerr << cell(i).index() << " (" << cellName[i] << ")\t";
      for (size_t k=0; k<cell(i).numVertex(); ++k)
	std::cerr << cell(i).vertex(k)->index() << " ";
      std::cerr << std::endl;
    }
    std::cerr << "Walls:" << std::endl;
    for (size_t i=0; i<numWall(); ++i) {
      std::cerr << wall(i).index() << "\t";
      std::cerr << wall(i).cell1()->index() << " " << wall(i).cell2()->index()<< " "
		<< wall(i).vertex1()->index() << " " << wall(i).vertex2()->index()
		<< std::endl;
    }
    
  }
  
  if (verbose)
    std::cerr << numCell() << " cells and " << numVertex() 
	      << " vertices and " << numWall() << " walls extracted by "
	      << "readInitMerryProj()" << std::endl;
  sortCellWallAndCellVertex();
  checkConnectivity(verbose);
  
  return;
  
  // Extract internal walls
  for (size_t i=0; i<numCell(); ++i) {
    for (size_t ii=i+1; ii<numCell(); ++ii) {
      std::vector<Vertex*> commonVertex;
      for (size_t j=0; j<cell(i).numVertex(); ++j)
	for (size_t jj=0; jj<cell(ii).numVertex(); ++jj)
	  if (cell(i).vertex(j)==cell(ii).vertex(jj))
	    commonVertex.push_back(cell(i).vertex(j));
      if (commonVertex.size()==2) {
	//cell i and ii joined by wall connected to vertices j and jj
	size_t numWallBefore=numWall();
	Wall wallTmp;
	wallTmp.setIndex(numWallBefore);
	wallTmp.setCell(&cell(i),&cell(ii));
	wallTmp.setVertex(commonVertex[0],commonVertex[1]);
	double length=0.0;
	for (size_t dim=0; dim<commonVertex[0]->numPosition(); ++dim)
	  length += ( commonVertex[0]->position(dim)-commonVertex[1]->position(dim) ) *
	    ( commonVertex[0]->position(dim)-commonVertex[1]->position(dim) );
	length = std::sqrt(length);
	wallTmp.setLength(length);
	if( verbose>1 )
	  std::cerr << "wall " << numWallBefore << " added with length "
		    << wallTmp.length() << " and connected to cells"
		    << cell(i).index() << "," << cell(ii).index() << " and vertices "
		    << commonVertex[0]->index() << "," << commonVertex[1]->index()
		    << std::endl; 
	addWall(wallTmp);
	// Add wall to cells and vertices
	cell(i).addWall( &(wall(numWallBefore)) );
	cell(ii).addWall( &(wall(numWallBefore)) );
	commonVertex[0]->addWall( &(wall(numWallBefore)) );
	commonVertex[1]->addWall( &(wall(numWallBefore)) );
      }
    }
  }
  // Extract walls towards boundary
  //assert( c1==static_cast<size_t>(-1) || c1<numCell() );
  for (size_t i=0; i<numCell(); ++i) {
    if (cell(i).numVertex() != cell(i).numWall() && cell(i).numVertex()>2 ) {
      //Sort vertices
      size_t Nv = cell(i).numVertex();
      std::vector<double> cellCenter = cell(i).positionFromVertex();
      std::vector<double> theta(Nv);
      std::vector<size_t> vI(Nv);
      for (size_t k=0; k<Nv; ++k) {				
	vI[k]=k;
	theta[k] = std::atan2( cell(i).vertex(k)->position(1)-cellCenter[1],
			       cell(i).vertex(k)->position(0)-cellCenter[0] );
      }
      for (size_t k=0; k<Nv; ++k)	
	for (size_t kk=k+1; kk<Nv; ++kk)				
	  if (theta[kk]<theta[k]) {
	    size_t tmpvI=vI[kk];
	    double tmpTheta=theta[kk];
	    vI[kk]=vI[k];
	    theta[kk]=theta[k];
	    vI[k]=tmpvI;
	    theta[k]=tmpTheta;
	  }
      for (size_t k=0; k<Nv; ++k) {
	size_t vI1 = vI[k];
	size_t vI2 = vI[(k+1)%Nv];
	//check if neigh vertices has wall inbetween
	size_t hasWall=0;
	for (size_t v=0; v<cell(i).vertex(vI1)->numWall(); ++v)
	  for (size_t vv=0; vv<cell(i).vertex(vI2)->numWall(); ++vv)
	    if (cell(i).vertex(vI1)->wall(v) == cell(i).vertex(vI2)->wall(vv) )
	      hasWall++;
	if (hasWall>1) {
	  std::cerr << "Tissue::readInitMerryProj() More than one vertex pair are connected "
		    << "for cell " << i << " (" << cellName[i] << ")" << std::endl;
	  std::cerr << "Vertices: ";
	  for (size_t kk=0; kk<cell(i).numVertex(); ++kk)
	    std::cerr << cell(i).vertex(kk)->index() << " ";
	  std::cerr << std::endl;
	  std::cerr << "Walls: ";
	  for (size_t kk=0; kk<cell(i).numWall(); ++kk)
	    std::cerr << cell(i).wall(kk)->index() << " ";
	  std::cerr << std::endl;					
	  for (size_t v=0; v<cell(i).vertex(vI1)->numWall(); ++v)
	    for (size_t vv=0; vv<cell(i).vertex(vI2)->numWall(); ++vv)
	      if (cell(i).vertex(vI1)->wall(v) == cell(i).vertex(vI2)->wall(vv) )
		std::cerr << cell(i).vertex(vI1)->wall(v)->index()<< " " 
			  <<  cell(i).vertex(vI2)->wall(vv)->index() << std::endl;
	  exit(-1);
	}
	if (!hasWall) {
	  //Add wall between vertices and add cell and bg
	  size_t numWallBefore=numWall();
	  Wall tmpWall;
	  tmpWall.setIndex(numWallBefore);
	  tmpWall.setCell(&cell(i),background());
	  tmpWall.setVertex(cell(i).vertex(vI1),cell(i).vertex(vI2));
	  double length=0.0;
	  for (size_t dim=0; dim<cell(i).vertex(vI1)->numPosition(); ++dim)
	    length += ( cell(i).vertex(vI1)->position(dim)-
			cell(i).vertex(vI2)->position(dim) ) *
	      ( cell(i).vertex(vI1)->position(dim)-
		cell(i).vertex(vI2)->position(dim) );
	  length = std::sqrt(length);
	  tmpWall.setLength(length);
	  if (verbose>1)
	    std::cerr << "wall " << numWallBefore << " added with length "
		      << tmpWall.length() << " and connected to cells "
		      << cell(i).index() << "," << background()->index() 
		      << " and vertices "
		      << cell(i).vertex(vI1)->index() << "," 
		      << cell(i).vertex(vI2)->index()
		      << std::endl; 
	  addWall(tmpWall);
	  // Add wall to cells and vertices
	  cell(i).addWall( &(wall(numWallBefore)) );
	  cell(i).vertex(vI1)->addWall( &(wall(numWallBefore)) );
	  cell(i).vertex(vI2)->addWall( &(wall(numWallBefore)) );					
	}
      }
    }
    else {
      removeCell(i);
    }
  }
  if (verbose)
    std::cerr << numCell() << " cells and " << numVertex() 
	      << " vertices and " << numWall() << " walls extracted by "
	      << "Tissue::readInitMerryProj()" << std::endl;
  checkConnectivity(verbose);
  sortCellWallAndCellVertex();
  checkConnectivity(verbose);
}

void Tissue::readInitMGXTriCell( const char *initFile, int verbose ) 
{
  std::ifstream IN(initFile);
  if( !IN ) {
    std::cerr << "Tissue::readInitMGXTriCell(char*) - "
	      << "Cannot open file " << initFile 
	      << std::endl; 
    exit(EXIT_FAILURE);
  }
  unsigned int numVertexVal,dimension=3;//assuming always 3
  IN >> numVertexVal;
  setNumVertex( numVertexVal );
  for( size_t i=0 ; i<numVertexVal ; ++i )
    vertex(i).setIndex(i);
  
  std::vector<size_t> vLabel(numVertexVal);
  std::vector<double> pos(dimension);
  std::vector<size_t> cellName,vertexName;
  // Read information about vertices
  for( size_t i=0 ; i<numVertexVal ; ++i ) {
    size_t tmp;
    IN >> tmp;
    vertexName.push_back(tmp);
    if (tmp!=vertex(i).index()) {
      std::cerr << "Tissue::readInitMGXTriCell() Expecting consecutive indices"
		<< " in file." << std::endl;
      exit(EXIT_FAILURE);
    }
    for( size_t dim=0 ; dim<dimension ; ++dim )
      IN >> pos[dim];
    vertex(i).setPosition(pos);
    IN >> vLabel[i];
    std::string sTmp;
    IN >> sTmp;
    //Check sTmp==j?
  }
  
  // Read vertex connectivity and create walls
  size_t wallIndex=0;
  for( size_t i=0 ; i<numVertexVal ; ++i ) {
    size_t indexVal;
    IN >> indexVal;
    if (i!=indexVal) {
      std::cerr << "Tissue::readInitMGXTriCell() Expecting consecutive indices"
		<< " in file when reading connectivity." << std::endl;
      exit(EXIT_FAILURE);
    }
    size_t numVertexNeigh;
    IN >> numVertexNeigh;
    for( size_t k=0 ; k<numVertexNeigh ; ++k ) {
      size_t j;
      IN >> j;
      if (i<j) {
	Wall tmpWall;
	tmpWall.setIndex(wallIndex);
	tmpWall.setVertex(vertexP(i),vertexP(j));
	tmpWall.setLength(tmpWall.lengthFromVertexPosition());
	tmpWall.setCell(background(),background());//temporary, generate cells below
	addWall(tmpWall);
	vertex(i).addWall(wallP(wallIndex));
	vertex(j).addWall(wallP(wallIndex));
      }
    }
  }
  IN.close();
  
  // Generate cells (assuming triangeles)
  
  
  // Mark wall boundaries and indices
  //for (size_t i=0; i<numWall(); ++i) {
  //if (wall(i).cell1().variable(0)==wall.cell2().variable(0)) {
  //	wall(i).addVariable(1);
  //}
  //else {
  //	wall(i).addVariable(0);
  //}
  //}
  if (verbose) {
    std::cerr << numCell() << " cells and " << numVertex() 
	      << " vertices and " << numWall() << " walls extracted by "
	      << "Tissue::readInitMGXTriCell()" << std::endl;
  }
  sortCellWallAndCellVertex();
  checkConnectivity(verbose);
  
  return;
}

void Tissue::readInitMGXTriVtu( const char *initFile, int verbose ) 
{
  std::ifstream IN(initFile);
  if( !IN ) {
    std::cerr << "Tissue::readInitMGXTriVtu(char*) - "
	      << "Cannot open file " << initFile 
	      << std::endl; 
    exit(EXIT_FAILURE);
  }
  unsigned int numVertexVal,dimension=3,dim;//assuming always 3
  IN >> numVertexVal;
  IN >> dim;
  assert (dim==dimension);
  
  setNumVertex( numVertexVal );
  
  std::vector<double> pos(dimension);
  //std::vector<size_t> cellName,vertexName;
  // Read information about vertices, positions first and label next
  for( size_t i=0 ; i<numVertexVal ; ++i ) {
    vertex(i).setIndex(i);
    //vertex(i).setId(i);
    for( size_t dim=0 ; dim<dimension ; ++dim ) {
      IN >> pos[dim];
    }
    vertex(i).setPosition(pos);
  }
  size_t numVtmp,numLtmp;
  IN >> numVtmp;
  if (numVtmp != numVertexVal) {
    std::cerr << "Tissue::readInitMGXTriVtu(char*) - "
	      << "Number of vertex positions not same as vertex labels." 
	      << std::endl; 
    exit(EXIT_FAILURE);
  }
  IN >> numLtmp;
  assert (numLtmp==1);
  std::vector<int> vLabel(numVertexVal);
  for( size_t i=0 ; i<numVertexVal ; ++i ) {
    IN >> vLabel[i];
  }
  // Read cell connectivity
  size_t numCellVal,numCellVertex;
  IN >> numCellVal;
  
  setNumCell( numCellVal );
  
  IN >> numCellVertex;
  if (numCellVertex != 3) {
    std::cerr << "Tissue::readInitMGXTriVtu(char*) - "
	      << "Number of vertex per cell expected to be 3."
	      << std::endl; 
    exit(EXIT_FAILURE);
  }	
  for( size_t i=0 ; i<numCellVal ; ++i ) {
    std::vector<size_t> cellV(numCellVertex);
    cell(i).setIndex(i);
    //cell(i).setId(i);
    for( size_t k=0 ; k<numCellVertex ; ++k ) {
      IN >> cellV[k];
      cell(i).addVertex(vertexP(cellV[k]));
      vertex(cellV[k]).addCell(cellP(i));
    }
  }  
  // Read cell label
  size_t numCtmp;
  IN >> numCtmp;
  if (numCtmp != numCellVal) {
    std::cerr << "Tissue::readInitMGXTriVtu(char*) - "
	      << "Number of cell connections (vertex) not same as cell labels." 
	      << std::endl; 
    exit(EXIT_FAILURE);
  }
  IN >> numLtmp;
  assert (numLtmp==1);
  for( size_t i=0 ; i<numCellVal ; ++i ) {
    size_t cellLabel;
    IN >> cellLabel;
    // Add direction in front of variables
    cell(i).addVariable(1);
    cell(i).addVariable(0);
    cell(i).addVariable(0);
    cell(i).addVariable(1);
    cell(i).addVariable(cellLabel);
  }
  IN.close();
  
  size_t wallIndex=0;
  // Create walls from cells and their vertex connections
  for( size_t i=0 ; i<numCellVal ; ++i ) {    
    // Check potential wall for all permutations of vertex pairs (assuming triangles)
    assert (cell(i).numVertex()==3);
    for (size_t k=0; k<cell(i).numVertex(); ++k) {
      for (size_t kk=k+1; kk<cell(i).numVertex(); ++kk) {
	int vertexWallK = cell(i).vertex(k)->isNeighborViaWall(cell(i).vertex(kk));
	if ( vertexWallK != -1 ) {
	  // Wall exist, just check it has an appropriate cell connected and add the cell
	  size_t vertexWallI = static_cast<size_t>(vertexWallK);
	  cell(i).vertex(k)->wall(vertexWallI)->setCell2( cellP(i) );
	  cell(i).addWall( cell(i).vertex(k)->wall(vertexWallI) );
	}
	else {
	  // Create a new wall
	  Wall tmpWall;
	  tmpWall.setIndex(wallIndex);
	  tmpWall.setVertex(vertexP(cell(i).vertex(k)->index()),
			    vertexP(cell(i).vertex(kk)->index()));
	  tmpWall.setLength(tmpWall.lengthFromVertexPosition());
	  tmpWall.setCell(cellP(i),background());
	  // Add variable indicating if boundary between cell labels (proper wall=1)
	  //if (vLabel[tmpWall.vertex1()->index()] == -1 && vLabel[tmpWall.vertex1()->index()] == -1) {
	  //tmpWall.addVariable(1.0);
	  //}
	  //else {
	  //tmpWall.addVariable(0.0);
	  //}
	  addWall(tmpWall);
	  vertex(cell(i).vertex(k)->index()).addWall(wallP(wallIndex));
	  vertex(cell(i).vertex(kk)->index()).addWall(wallP(wallIndex));
	  cell(i).addWall(wallP(wallIndex));
	  ++wallIndex;
	}
      }
    }
  }
  // Add variable indicating if boundary between cell labels (proper wall=1)
  for (size_t i=0; i<numWall(); ++i) {
    if (wall(i).cell1() == background() || wall(i).cell2() == background() 
	|| wall(i).cell1()->variable(4) != wall(i).cell2()->variable(4)) {
      wall(i).addVariable(1.0);
    }
    else {
      wall(i).addVariable(0.0);
    }
  }
  

  if (verbose) {
    std::cerr << numCell() << " cells and " << numVertex() 
	      << " vertices and " << numWall() << " walls extracted by "
	      << "readMGXTriInit()" << std::endl;
  }
  sortCellWallAndCellVertex();
  checkConnectivity(verbose);
  
  return;
}

void Tissue::readInitMGXTriMesh( const char *initFile, int verbose ) 
{
  std::ifstream IN(initFile);
  if( !IN ) {
    std::cerr << "Tissue::readInitMGXTriMesh(char*) - "
	      << "Cannot open file " << initFile 
	      << std::endl; 
    exit(EXIT_FAILURE);
  }
  // Read and check header
  //
  std::string tmpString;
  int tmpInt;
  IN >> tmpString;
  IN >> tmpInt;
  if (tmpString != "MeshVersionFormatted" || tmpInt!=1) { 
    std::cerr << "Tissue::readInitMGXTriMesh Expecting 'MeshVersionFormatted 1' on first line of file, "
	      << "not '" << tmpString << " " << tmpInt << "'." << std::endl;
    exit(EXIT_FAILURE);
  }
  IN >> tmpString;
  size_t dimension;
  IN >> dimension;
  if (tmpString != "Dimension") { 
    std::cerr << "Tissue::readInitMGXTriMesh Expecting 'Dimension d' on second line of file, "
	      << "not '" << tmpString << " " << dimension << "'." << std::endl;
    exit(EXIT_FAILURE);
  }
  IN >> tmpString;
  if (tmpString != "Vertices") { 
    std::cerr << "Tissue::readInitMGXTriMesh Expecting 'Vertices' on third line of file, "
	      << "not '" << tmpString << "'." << std::endl;
    exit(EXIT_FAILURE);
  }

  // Read and create vertices and store labels
  //
  unsigned int numVertexVal;
  IN >> numVertexVal;
  setNumVertex( numVertexVal );
  
  std::vector<double> pos(dimension);
  std::vector<int> vLabel(numVertexVal);
  for( size_t i=0 ; i<numVertexVal ; ++i ) {
    vertex(i).setIndex(i);
    //vertex(i).setId(i);
    for( size_t dim=0 ; dim<dimension ; ++dim ) {
      IN >> pos[dim];
    }
    vertex(i).setPosition(pos);
    IN >> vLabel[i];
  }

  IN >> tmpString;
  if (tmpString != "Triangles") { 
    std::cerr << "Tissue::readInitMGXTriMesh Expecting 'Triangles' on line after reading vertices, "
	      << "not '" << tmpString << "'." << std::endl;
    exit(EXIT_FAILURE);
  }

  // Read cell (triangle) connectivity and cell labels
  // NOTE: In this format the vertex indices start at 1 (not zero).
  size_t numCellVal,numCellVertex=3;// Assuming triangles in this file
  IN >> numCellVal;
  setNumCell( numCellVal );
  
  for (size_t i=0; i<numCellVal; ++i) {
    std::vector<size_t> cellV(numCellVertex);
    cell(i).setIndex(i);
    //cell(i).setId(i);
    for( size_t k=0 ; k<numCellVertex ; ++k ) {
      IN >> cellV[k];
      cellV[k]--;//Since vertex indices in this list starts from one.
      cell(i).addVertex(vertexP(cellV[k]));
      vertex(cellV[k]).addCell(cellP(i));
    }

    int cellLabel;
    IN >> cellLabel;
    // Add direction in front of variables
    cell(i).addVariable(1);
    cell(i).addVariable(0);
    cell(i).addVariable(0);
    cell(i).addVariable(1);
    cell(i).addVariable(cellLabel);
  }

  IN >> tmpString;
  if (tmpString != "End") { 
    std::cerr << "Tissue::readInitMGXTriMesh Expecting 'End' on final line after reading cells, "
	      << "not '" << tmpString << "'." << std::endl;
    //exit(EXIT_FAILURE);
  }

  IN.close();
  
  // Create walls from cells and their vertex connections
  //
  size_t wallIndex=0;
  for( size_t i=0 ; i<numCellVal ; ++i ) {    
    // Check potential wall for all permutations of vertex pairs (assuming triangles)
    assert (cell(i).numVertex()==3);
    for (size_t k=0; k<cell(i).numVertex(); ++k) {
      for (size_t kk=k+1; kk<cell(i).numVertex(); ++kk) {
	int vertexWallK = cell(i).vertex(k)->isNeighborViaWall(cell(i).vertex(kk));
	if ( vertexWallK != -1 ) {
	  // Wall exist, just check it has an appropriate cell connected and add the cell
	  size_t vertexWallI = static_cast<size_t>(vertexWallK);
	  cell(i).vertex(k)->wall(vertexWallI)->setCell2( cellP(i) );
	  cell(i).addWall( cell(i).vertex(k)->wall(vertexWallI) );
	}
	else {
	  // Create a new wall
	  Wall tmpWall;
	  tmpWall.setIndex(wallIndex);
	  tmpWall.setVertex(vertexP(cell(i).vertex(k)->index()),
			    vertexP(cell(i).vertex(kk)->index()));
	  tmpWall.setLength(tmpWall.lengthFromVertexPosition());
	  tmpWall.setCell(cellP(i),background());
	  // Add variable indicating if boundary between cell labels (proper wall=1)
	  //if (vLabel[tmpWall.vertex1()->index()] == -1 && vLabel[tmpWall.vertex1()->index()] == -1) {
	  //tmpWall.addVariable(2.0);
	  //}
	  //else if (vLabel[tmpWall.vertex1()->index()] == -1 && vLabel[tmpWall.vertex1()->index()] == -1) {
	  //tmpWall.addVariable(1.0);
	  //}
	  //else {
	  //tmpWall.addVariable(0.0);
	  //}
	  addWall(tmpWall);
	  vertex(cell(i).vertex(k)->index()).addWall(wallP(wallIndex));
	  vertex(cell(i).vertex(kk)->index()).addWall(wallP(wallIndex));
	  cell(i).addWall(wallP(wallIndex));
	  ++wallIndex;
	}
      }
    }
  }
  // Add variable indicating if boundary between cell labels (proper wall=1)
  for (size_t i=0; i<numWall(); ++i) {
    if (wall(i).cell1() == background() || wall(i).cell2() == background() 
	|| wall(i).cell1()->variable(4) != wall(i).cell2()->variable(4)) {
      wall(i).addVariable(1.0);
    }
    else {
      wall(i).addVariable(0.0);
    }
  }
  
  if (verbose) {
    std::cerr << numCell() << " cells and " << numVertex() 
	      << " vertices and " << numWall() << " walls extracted by "
	      << "Tissue::readInitMGXTriMesh()" << std::endl;
  }
  sortCellWallAndCellVertex();
  checkConnectivity(verbose);
  
  return;
}

void Tissue::readInitPLYMesh(const char *initFile, int verbose)
{
  //PLY_file ply_file(initFile);
  //PLY_reader ply_reader;
  //ply_reader.read(ply_file, *this);
}

void Tissue::readModel(std::ifstream &IN,int verbose) {
  
  unsigned int numReactionVal,numCompartmentChangeVal,numDirection;
  
  if( verbose )
    std::cerr << "Reading model file:\n";
  IN >> numReactionVal;
  IN >> numCompartmentChangeVal;
  IN >> numDirection;
  assert(numDirection==0 || numDirection==1);
  
  //Read Reactions
  //
  //Remove any present reactions before adding
  if( numReaction() ) 
    reaction_.resize(0);
  
  if( verbose )
    std::cerr << "reactions...\n"; 
  for( size_t i=0 ; i<numReactionVal ; i++ ) {
    if( addReaction(IN) )
      std::cerr << "Tissue::ReadModel(ifstream) "
		<< "Warning Adding reaction failed for "
		<< "tissue " << id() << " (index " << i << ")\n";
    else if( verbose )
      std::cerr << reaction(numReaction()-1)->id() << std::endl;
  }
  //Read compartmentChanges
  //
  //Remove any present compartmentChanges before adding
  if( numCompartmentChange() ) 
    compartmentChange_.resize(0);
  
  if( verbose )
    std::cerr << "compartment changes...\n";
  for( size_t i=0 ; i<numCompartmentChangeVal ; i++ ) {
    if( addCompartmentChange(IN) ) 
      std::cerr << "Tissue::ReadModel(ifstream) "
		<< "Warning Adding compartmentChange failed for "
		<< "tissue " << id() << " (index " << i << ")\n";  
    else if( verbose )
      std::cerr << compartmentChange(numCompartmentChange()-1)->id() 
		<< std::endl;
  }
  // Read direction if applicable
  if( verbose )
    std::cerr << "direction...\n";
  if( numDirection ) {
    if( direction()->readDirection(IN) ) {
      std::cerr << "Tissue::ReadModel(ifstream) "
		<< "Adding direction failed." << std::endl;
      exit(-1);
    }
    else if( verbose )
      std::cerr << direction()->directionUpdate()->id() << std::endl
		<< direction()->directionDivision()->id() << std::endl;
  }			
  if( verbose )
    std::cerr << "Done\n\n";
}

void Tissue::readModel(const char *fileName, int verbose) 
{
  std::string tmp(fileName);
  readModel(tmp,verbose);
}

void Tissue::readModel(std::string fileName, int verbose) 
{  
  const char* fName = fileName.c_str();
  std::istream *IN = myFiles::openFile(fName);
  if( !IN ) {
    std::cerr << "Tissue::readModel(std::string) - "
	      << "Cannot open file " << fileName << "\n\n\7";
    exit(-1);
  }
  readModel((std::ifstream &) *IN,verbose);
}

size_t wallFromCellPair(std::vector< std::pair<size_t,size_t> > &wallCell,
			size_t c1,size_t c2) {
  
  for( size_t i=0 ; i<wallCell.size() ; ++i )
    if( (wallCell[i].first==c1 && wallCell[i].second==c2) ||
	(wallCell[i].first==c2 && wallCell[i].second==c1) )
      return i;
  return static_cast<size_t>(-1);
}

void Tissue::readSphereInit( const char *initFile, int verbose ) 
{
  //Read the sphere data 
  std::ifstream IN(initFile);
  if( !IN ) {
    std::cerr << "Tissue::readSphereInit() Cannot open file " << initFile << std::endl; 
    exit(EXIT_FAILURE);
  }
  size_t numCell,numCol;
  IN >> numCell;
  IN >> numCol;
  //If not given asssume only positions and radii in file
  size_t dimension = numCol-1;
  if( dimension<2 || dimension>3 ) {
    std::cerr << "Tissue::readSphereInit() Dimension "
	      << "(read by number of columns-1 in the file) needs to be 2 or 3." 
	      << std::endl;
    exit(EXIT_FAILURE);
  }
  std::vector< std::vector<double> > data(numCell);
  double tmp;
  for( size_t i=0 ; i<numCell ; ++i ) {
    data[i].resize( dimension+1 );
    for( size_t j=0 ; j<data[i].size() ; ++j )
      IN >> data[i][j];
    for( size_t j=data[i].size() ; j<numCol ; ++j )
      IN >> tmp;		
  }  
  //Create the tissue from the sphere data
  double rFac=1.0;
  createTissueFromSpheres(data,rFac,verbose);  
}

void Tissue::
createTissueFromSpheres(DataMatrix &y,
			double rFac, int verbose) 
{  
  size_t N = y.size();
  if( !N ) return;
  std::vector< std::vector<size_t> > cellCellNeighbor(N),
    cellWallNeighbor(N);
  double d,r;
  unsigned int numWall=0,numVertex=0,num2Vertex=0,num3Vertex=0,
    num4Vertex=0;
  
  size_t rIndex = y[0].size()-1;
  if( rIndex<2 || rIndex>3 ) {
    std::cerr << "Tissue::createTissueFromSpheres()"
	      << "Only allowed for 2D and 3D.\n";
    exit(0);
  }
  //Get cell-cell neighbors which also defines the walls
  //
  std::vector< std::pair<size_t,size_t> > wallCell;
  std::vector< std::vector<size_t> > cellWall(N);
  for(size_t i=0 ; i<N ; ++i ) {
    for(size_t j=i+1 ; j<N ; ++j ) {
      r = y[i][rIndex] + y[j][rIndex];
      d=0.;
      for( size_t dim=0 ; dim<rIndex ; dim++ )
	d += (y[i][dim]-y[j][dim])*(y[i][dim]-y[j][dim]);
      d = std::sqrt( d );
      
      if( d<=r*rFac ) {
        //Add cell-cell neighbors
        cellCellNeighbor[i].push_back(j);
        cellCellNeighbor[j].push_back(i);
	cellWall[i].push_back( numWall );
	cellWall[j].push_back( numWall );
	wallCell.push_back( std::pair<size_t,size_t>(i,j) );
        numWall++;
	if( verbose>1 )
	  std::cerr << "Wall between cells " << i << "," << j << std::endl;
      }
    }
  }
  
  //Get verteces by checking cell neighbor relationships
  //
  std::vector< std::pair<size_t,size_t> > 
    wallVertex(numWall,std::pair<size_t,size_t> (static_cast<size_t>(-1),
						 static_cast<size_t>(-1)));
  std::vector< std::vector<size_t> > cellVertex(N),vertexCell,vertexWall;
  for(size_t i1=0 ; i1<N ; ++i1 ) {
    for(size_t k1=0 ; k1<cellCellNeighbor[i1].size() ; ++k1 ) {
      size_t i2=cellCellNeighbor[i1][k1]; 
      if( true ) { //if( i2>i1 ) { //For not double checking
	//get vector of common neighbors
	std::vector<size_t> common;
	for(size_t k2=k1+1 ; k2<cellCellNeighbor[i1].size() ; ++k2 ) {
	  size_t i3=cellCellNeighbor[i1][k2]; 
	  std::vector<size_t>::iterator it32 = 
	    find(cellCellNeighbor[i2].begin(),
		 cellCellNeighbor[i2].end(),i3);
	  //if( i3>i2 && it32 != cellCellNeighbor[i2].end() ) //common found
	  if( true && it32 != cellCellNeighbor[i2].end() ) //common found
	    common.push_back(i3);
	}
	//std::cerr << common.size() << " common indeces found for cells "
	//  << i1 << " " << i2 << " ( ";
	//for(size_t tmp=0 ; tmp<common.size() ; ++tmp )
	//std::cerr << common[tmp] << " ";
	//std::cerr << ")" << std::endl;
	//Check pairwise for neighbors (to identify 4-verteces)
	std::vector< std::pair<size_t,size_t> > pairs;
	std::vector<size_t> commonMarkedForPair(common.size());	
	for(size_t c1=0 ; c1<common.size() ; ++c1 )
	  for(size_t c2=c1+1 ; c2<common.size() ; ++c2 ) {
	    std::vector<size_t>::iterator cit = 
	      find(cellCellNeighbor[common[c1]].begin(),
		   cellCellNeighbor[common[c1]].end(),common[c2]);
	    if( cit != cellCellNeighbor[common[c1]].end() ) {//pair found
	      pairs.push_back(std::pair<size_t,size_t>(common[c1],common[c2]));
	      commonMarkedForPair[c1]++;
	      commonMarkedForPair[c2]++;
	    }
	  }	
	//Add all 3-vertex
	for(size_t c=0 ; c<common.size() ; ++c )
	  if( commonMarkedForPair[c]==0 ) {
	    std::vector<size_t> tmpVec(3);
	    tmpVec[0]=i1;tmpVec[1]=i2;tmpVec[2]=common[c];
	    size_t addFlag=1;
	    for( size_t v=0 ; v<vertexCell.size() ; ++v ) {
	      std::vector<size_t>::iterator cit0 = 
		find(vertexCell[v].begin(),
		     vertexCell[v].end(),tmpVec[0]);
	      std::vector<size_t>::iterator cit1 = 
		find(vertexCell[v].begin(),
		     vertexCell[v].end(),tmpVec[1]);
	      std::vector<size_t>::iterator cit2 = 
		find(vertexCell[v].begin(),
		     vertexCell[v].end(),tmpVec[2]);
	      if( cit0 != vertexCell[v].end() &&
		  cit1 != vertexCell[v].end() &&
		  cit2 != vertexCell[v].end() ) {//vertex already found
		addFlag=0;
		break;
	      }
	    }
	    if( addFlag ) {
	      vertexCell.push_back(tmpVec);
	      vertexWall.push_back(std::vector<size_t>(0));
	      for( size_t iTmp=0 ; iTmp<tmpVec.size() ; ++iTmp )
		cellVertex[tmpVec[iTmp]].push_back(numVertex);
	      for( size_t iTmp=0 ; iTmp<tmpVec.size() ; ++iTmp )
		for( size_t jTmp=iTmp+1 ; jTmp<tmpVec.size() ; ++jTmp ) {
		  size_t wallTmp = wallFromCellPair(wallCell,tmpVec[iTmp],tmpVec[jTmp]);
		  if( wallTmp<numWall ) {
		    if( wallVertex[wallTmp].first>numWall )
		      wallVertex[wallTmp].first = numVertex;
		    else if( wallVertex[wallTmp].second>numWall )
		      wallVertex[wallTmp].second = numVertex;
		    else
		      std::cerr << "Warning, trying to add a third "
				<< "vertex to wall " << wallTmp 
				<< std::endl; 
		    vertexWall[vertexWall.size()-1].push_back(wallTmp);
		  }
		}
	      if( verbose>1 )
		std::cerr << "3-vertex for " << i1 << " " << i2 << " "
			  << common[c] << " added" << std::endl;
	      numVertex++;
	      num3Vertex++;
	    }
	  }	
	  else if( commonMarkedForPair[c]==1 ) {
	    size_t i3Tmp=common[c],i4Tmp=common[c],tmpCount=0;
	    for(size_t p=0 ; p<pairs.size() ; ++p )
	      if( pairs[p].first==i3Tmp ) {
		i4Tmp=pairs[p].second;
		tmpCount++;
	      }
	      else if( pairs[p].second==i3Tmp ) {
		i4Tmp=pairs[p].first;
		tmpCount++;
	      }
	    //if( tmpCount==1 && i3Tmp<i4Tmp ) {
	    if( tmpCount==1 && true ) {
	      std::vector<size_t> tmpVec(4);
	      tmpVec[0]=i1;tmpVec[1]=i2;tmpVec[2]=i3Tmp;tmpVec[3]=i4Tmp;
	      size_t addFlag=1;
	      for( size_t v=0 ; v<vertexCell.size() ; ++v ) {
		std::vector<size_t>::iterator cit0 = 
		  find(vertexCell[v].begin(),
		       vertexCell[v].end(),tmpVec[0]);
		std::vector<size_t>::iterator cit1 = 
		  find(vertexCell[v].begin(),
		       vertexCell[v].end(),tmpVec[1]);
		std::vector<size_t>::iterator cit2 = 
		  find(vertexCell[v].begin(),
		       vertexCell[v].end(),tmpVec[2]);
		std::vector<size_t>::iterator cit3 = 
		  find(vertexCell[v].begin(),
		       vertexCell[v].end(),tmpVec[3]);
		if( cit0 != vertexCell[v].end() &&
		    cit1 != vertexCell[v].end() &&
		    cit2 != vertexCell[v].end() &&
		    cit3 != vertexCell[v].end() ) {//vertex already found
		  addFlag=0;
		  break;
		}
	      }
	      if( addFlag ) {
		vertexCell.push_back(tmpVec);
		vertexWall.push_back(std::vector<size_t>(0));
		for( size_t iTmp=0 ; iTmp<tmpVec.size() ; ++iTmp )
		  cellVertex[tmpVec[iTmp]].push_back(numVertex);
		for( size_t iTmp=0 ; iTmp<tmpVec.size() ; ++iTmp )
		  for( size_t jTmp=iTmp+1 ; jTmp<tmpVec.size() ; ++jTmp ) {
		    size_t wallTmp = wallFromCellPair(wallCell,tmpVec[iTmp],tmpVec[jTmp]);
		    if( wallTmp<numWall ) {
		      if( wallVertex[wallTmp].first>numWall )
			wallVertex[wallTmp].first = numVertex;
		      else if( wallVertex[wallTmp].second>numWall )
			wallVertex[wallTmp].second = numVertex;
		      else
			std::cerr << "Warning, trying to add a third "
				  << "vertex to wall " << wallTmp 
				  << std::endl; 
		      vertexWall[vertexWall.size()-1].push_back(wallTmp);
		    }
		    else
		      std::cerr << "Warning found cell pair without wall"
				<< std::endl;
		  }
		
		if( verbose>1 )
		  std::cerr << "4-vertex for " << i1 << " " << i2 << " "
			    << i3Tmp << " " << i4Tmp << " added" 
			    << std::endl;
		numVertex++;
		num4Vertex++;
		//std::cerr << "vertexCell:" << std::endl;
		//for( size_t iTmp=0 ; iTmp<vertexCell.size() ; ++iTmp ) {
		//for( size_t jTmp=0 ; jTmp<vertexCell[iTmp].size() ; ++jTmp )
		//  std::cerr << vertexCell[iTmp][jTmp] << " ";
		//std::cerr << std::endl;
		//}
	      }
	    }
	    else if( tmpCount != 1 ) {
	      std::cerr << "Plausibel 4-vertex does not match pairs"
			<< std::endl;
	    }
	  }
	  else {
	    std::cerr << "Warning: possible 5-vertex found..." 
		      << std::endl;
	  }
      } 
    }
  }
  //std::cerr << "WallCell:" << std::endl;
  //for(size_t i=0 ; i<wallCell.size() ; ++i )
  //std::cerr << i << " - " << wallCell[i].first << " " 
  //      << wallCell[i].second << std::endl;
  //std::cerr << "WallVertex:" << std::endl;
  //for(size_t i=0 ; i<wallVertex.size() ; ++i )
  //std::cerr << i << " - " << wallVertex[i].first << " " 
  //      << wallVertex[i].second << std::endl;
  
  
  //Add vertecis to empty walls
  for( size_t wallI=0 ; wallI<numWall ; ++wallI ) {
    if( wallVertex[wallI].first>numVertex ) {
      wallVertex[wallI].first = numVertex;
      cellVertex[wallCell[wallI].first].push_back(numVertex);
      cellVertex[wallCell[wallI].second].push_back(numVertex);
      std::vector<size_t> tmpVec(1);tmpVec[0]=wallI;
      vertexWall.push_back(tmpVec);
      tmpVec.resize(2);
      tmpVec[0]=wallCell[wallI].first;
      tmpVec[1]=wallCell[wallI].second;
      vertexCell.push_back(tmpVec);
      if( verbose>1 )
	std::cerr << "Adding vertex " << numVertex << " to wall (first) "
		  << wallI << std::endl;
      numVertex++;
      num2Vertex++;
    }
    if( wallVertex[wallI].second>numVertex ) {
      wallVertex[wallI].second = numVertex;
      cellVertex[wallCell[wallI].first].push_back(numVertex);
      cellVertex[wallCell[wallI].second].push_back(numVertex);
      std::vector<size_t> tmpVec(1);tmpVec[0]=wallI;
      vertexWall.push_back(tmpVec);
      tmpVec.resize(2);
      tmpVec[0]=wallCell[wallI].first;
      tmpVec[1]=wallCell[wallI].second;
      vertexCell.push_back(tmpVec);
      if( verbose>1 )
	std::cerr << "Adding vertex " << numVertex << " to wall (second) "
		  << wallI << std::endl;
      numVertex++;
      num2Vertex++;
    }
  }
  //Add walls between two twoVertices in a cell (only if there are two)
  //
  for( size_t cellI=0 ; cellI<cellVertex.size() ; ++cellI ) {
    std::vector<size_t> twoVertexList;
    size_t numTwoVertex=0;
    for( size_t vertexI=0 ; vertexI<cellVertex[cellI].size() ; ++vertexI ) {
      if( vertexCell[ cellVertex[cellI][vertexI] ].size() == 2 ) {
	numTwoVertex++;
	twoVertexList.push_back( cellVertex[cellI][vertexI] );
      }
    }
    if( numTwoVertex == 2 && twoVertexList[0] != twoVertexList[1] ) {
      //Add new wall between these vertices and assume wall boundary towards
      //outside (ie one of the wallCell index is -1...
      size_t wallI = wallCell.size();
      numWall++;
      cellWall[cellI].push_back( wallI );
      vertexWall[twoVertexList[0]].push_back( wallI );
      vertexWall[twoVertexList[1]].push_back( wallI );
      wallCell.resize(wallI+1);
      wallVertex.resize(wallI+1);
      wallCell[wallI].first = cellI;
      wallCell[wallI].second = static_cast<size_t>(-1);
      wallVertex[wallI].first = twoVertexList[0];
      wallVertex[wallI].second = twoVertexList[1];
      if( verbose )
	std::cerr << "Wall " << wallI << " added in cell " << cellI 
		  << " (" << wallCell[wallI].second << ") between vertices " 
		  << twoVertexList[0] << " and "
		  << twoVertexList[1] << std::endl;
    }
  }
  
  std::cerr << N << " cells " << numWall << " walls and " << numVertex
	    << " vertices defined (" << num2Vertex << " 2v, " 
	    << num3Vertex << " 3v, " << num4Vertex << " 4v)"
	    << std::endl << std::endl;
  
  //Fill the tissue with cells, walls and verteces
  //
  if( verbose )
    std::cerr << "Tissue::createTissueFromSpheres() "
	      << "Creating the tissue." 
	      << std::endl;    
  setNumCell( N );
  setNumWall( numWall );
  setNumVertex( numVertex );
  for( size_t i=0 ; i<numCell() ; i++ ) {
    cell(i).setIndex(i);
    cell(i).addVariable( y[i][rIndex] );
    for( size_t j=0 ; j<cellWall[i].size() ; j++ )
      cell(i).addWall( &wall( cellWall[i][j] ) ); 
    for( size_t j=0 ; j<cellVertex[i].size() ; j++ )
      cell(i).addVertex( &vertex( cellVertex[i][j] ) ); 
  }
  for( size_t i=0 ; i<(this->numWall()) ; i++ ) {
    wall(i).setIndex(i);
    Cell *cell1,*cell2;
    //Check if wallCell is background
    if( wallCell[i].first<numCell() )
      cell1 = &cell(wallCell[i].first);
    else
      cell1 = &background_;
    if( wallCell[i].second<numCell() )
      cell2 = &cell(wallCell[i].second);
    else
      cell2 = &background_;
    
    wall(i).setCell(cell1,cell2);
    wall(i).setVertex(&vertex(wallVertex[i].first),
		      &vertex(wallVertex[i].second));
  }
  for( size_t i=0 ; i<(this->numVertex()) ; i++ ) {
    vertex(i).setIndex(i);
    for( size_t j=0 ; j<vertexCell[i].size() ; j++ )
      vertex(i).addCell( &cell( vertexCell[i][j] ) ); 
    for( size_t j=0 ; j<vertexWall[i].size() ; j++ )
      vertex(i).addWall( &wall( vertexWall[i][j] ) ); 
    size_t dimension = y[0].size()-1;
    std::vector<double> pos(dimension);
    for( size_t j=0 ; j<vertexCell[i].size() ; j++ )
      for( size_t d=0 ; d<dimension ; d++ )
	pos[d] += y[vertexCell[i][j]][d];
    for( size_t d=0 ; d<dimension ; d++ )
      pos[d] /= static_cast<double>( vertexCell[i].size() );
    vertex(i).setPosition(pos);
  }
  //Set the wall lengths from the vertices positions
  //
  if( verbose )
    std::cerr << "Tissue::createTissueFromSpheres() "
	      << "Setting wall lengths from vertex positions." 
	      << std::endl;
  setWallLengthFromVertexPosition();
  
  //Check that the tissue is ok
  checkConnectivity(verbose);
  
  //std::cerr << "CellWall:" << std::endl;
  //for(size_t i=0 ; i<cellWall.size() ; ++i ) {
  //std::cerr << i << " - ";
  //for(size_t j=0 ; j<cellWall[i].size() ; ++j )
  //  std::cerr << cellWall[i][j] << " ";
  //std::cerr << std::endl;
  //}
  //std::cerr << "CellVertex:" << std::endl;
  //for(size_t i=0 ; i<cellVertex.size() ; ++i ) {
  //std::cerr << i << " - ";
  //for(size_t j=0 ; j<cellVertex[i].size() ; ++j )
  //  std::cerr << cellVertex[i][j] << " ";
  //std::cerr << std::endl;
  //}
  //std::cerr << "WallCell:" << std::endl;
  //for(size_t i=0 ; i<wallCell.size() ; ++i )
  //std::cerr << i << " - " << wallCell[i].first << " " 
  //      << wallCell[i].second << std::endl;
  //std::cerr << "WallVertex:" << std::endl;
  //for(size_t i=0 ; i<wallVertex.size() ; ++i )
  //std::cerr << i << " - " << wallVertex[i].first << " " 
  //      << wallVertex[i].second << std::endl;
  //std::cerr << "VertexCell:" << std::endl;
  //for(size_t i=0 ; i<vertexCell.size() ; ++i ) {
  //std::cerr << i << " - ";
  //for(size_t j=0 ; j<vertexCell[i].size() ; ++j )
  //  std::cerr << vertexCell[i][j] << " ";
  //std::cerr << std::endl;
  //}
  //std::cerr << "VertexWall:" << std::endl;
  //for(size_t i=0 ; i<vertexWall.size() ; ++i ) {
  //std::cerr << i << " - ";
  //for(size_t j=0 ; j<vertexWall[i].size() ; ++j )
  //  std::cerr << vertexWall[i][j] << " ";
  //std::cerr << std::endl;
  //}  
}

void Tissue::readVoronoiInit( const char *initFile, int verbose ) 
{
  //Read the voronoi output data
  std::ifstream IN(initFile);
  if( !IN ) {
    std::cerr << "Tissue::readVoronoiInit() Cannot open file " << initFile << std::endl; 
    exit(EXIT_FAILURE);
  }
  size_t dimension, numVertex,numCell;
  int tmpI=0;
  double tmpD=0.0;
  //Read header information
  IN >> dimension;
  assert(dimension==2 || dimension==3);
  IN >> numVertex;
  --numVertex;
  IN >> numCell;
  IN >> tmpI;
  for( size_t j=0 ; j<dimension ; ++j )
    IN >> tmpD;
  
  //Read vertexPositions
  std::vector< std::vector<double> > vertexPos(numVertex);
  for( size_t i=0 ; i<vertexPos.size() ; ++i ) {
    vertexPos[i].resize(dimension);
    for( size_t j=0 ; j<dimension ; ++j )
      IN >> vertexPos[i][j];
  }
  
  //Read cellVertex connections
  std::vector< std::vector<size_t> > cellVertex(numCell);
  size_t numCellVertex=0;
  for( size_t i=0 ; i<cellVertex.size() ; ++i ) {
    IN >> numCellVertex;
    cellVertex[i].resize(numCellVertex);
    for( size_t j=0 ; j<numCellVertex ; ++j )
      IN >> cellVertex[i][j];
  }  
  // Create the tissue from the data and print the init
  //
  Tissue T;
  createTissueFromVoronoi(vertexPos,cellVertex,verbose);  
}

void Tissue::createTissueFromVoronoi(DataMatrix &vertexPos,
				     std::vector< std::vector<size_t> > &cellVertexTmp,
				     int verbose) 
{
  std::vector< std::vector<size_t> > cellVertex(cellVertexTmp.size()),
    cellWall(cellVertexTmp.size()),vertexCell(vertexPos.size());
  size_t boundaryIndex = static_cast<size_t>(-1);
  std::set<size_t> boundaryNeighVertex;
  //Convert cellVertexTmp to cellVertex by adding additional vertices
  //at boundary and lower each index by one
  //
  for( size_t i=0 ; i<cellVertexTmp.size() ; ++i ) {
    for( size_t k=0 ; k<cellVertexTmp[i].size() ; ++k ) {
      size_t vI1=static_cast<size_t>(cellVertexTmp[i][k]-1);
      size_t vI2=static_cast<size_t>(cellVertexTmp[i][(k+1)%cellVertexTmp[i].size()]-1);
      cellVertex[i].push_back( vI1 );
      if( vI2==boundaryIndex && vI1 != boundaryIndex ) {
	boundaryNeighVertex.insert(vI1);
      }				
      if( vI1==boundaryIndex ) {
	cellVertex[i].push_back( vI1 );
	if( vI2 != boundaryIndex )
	  boundaryNeighVertex.insert(vI2);
      }
      else
	vertexCell[vI1].push_back(i);
    }
  }
  std::cerr << "cellVertexTmp converted\n";
  //Add new vertices at the boundary
  //
  size_t numOldVertex=vertexPos.size();
  size_t numNewVertex=boundaryNeighVertex.size();
  vertexPos.resize(vertexPos.size()+numNewVertex,vertexPos[0]);
  vertexCell.resize(vertexPos.size());
  size_t kCount=0;
  for( std::set<size_t>::iterator k=boundaryNeighVertex.begin() ; 
       k!=boundaryNeighVertex.end() ; ++k ) {
    for( size_t i=0 ; i<cellVertex.size() ; ++i ) {
      for( size_t j=0 ; j<cellVertex[i].size() ; ++j ) {
	if( *k==cellVertex[i][j] ) {
	  size_t jPlus=(j+1)%cellVertex[i].size();
	  size_t jMinus = j!=0 ? j-1 : cellVertex[i].size()-1; 
	  if( cellVertex[i][jPlus]==boundaryIndex ) {
	    cellVertex[i][jPlus]=numOldVertex+kCount;
	    vertexCell[numOldVertex+kCount].push_back(i);
	  }
	  else if( cellVertex[i][jMinus]==boundaryIndex ) {
	    cellVertex[i][jMinus]=numOldVertex+kCount;
	    vertexCell[numOldVertex+kCount].push_back(i);
	  }
	}
      }
    }
    ++kCount;
  }
  std::cerr << "New vertices added\n";
  //Create walls
  //
  std::vector< std::vector<size_t> > vertexWall( vertexCell.size() );
  std::vector< std::pair<size_t,size_t> > wallCell,wallVertex;
  for( size_t i=0 ; i<cellVertex.size() ; ++i ) {
    for( size_t j=0 ; j<cellVertex[i].size() ; ++j ) {
      size_t jPlus=(j+1)%cellVertex[i].size();
      size_t vI1=cellVertex[i][j];
      size_t vI2=cellVertex[i][jPlus];
      //First check if it is a boundary wall (new wall)
      if( vI1>=numOldVertex && vI2>=numOldVertex ) {
	cellWall[i].push_back(wallCell.size());
	std::pair<size_t,size_t> tmpPair(i,boundaryIndex);
	wallCell.push_back(tmpPair);
	vertexWall[vI1].push_back(wallVertex.size());
	vertexWall[vI2].push_back(wallVertex.size());
	tmpPair.first=vI1;
	tmpPair.second=vI2;
	wallVertex.push_back(tmpPair);
      }
      else {
	//Check if already in wall
	size_t inWallFlag=0;
	for( size_t k=0 ; k<wallVertex.size() ; ++k ) {
	  if( ( wallVertex[k].first==vI1 && wallVertex[k].second==vI2 ) ||
	      ( wallVertex[k].first==vI2 && wallVertex[k].second==vI1 ) ) {
	    ++inWallFlag;
	    if( wallCell[k].first != wallCell[k].second ) {
	      std::cerr << "Tissue::createInitFromVoronoi() Wall not"
			<< " marked for additional cell." << std::endl;
	      exit(-1);
	    }
	    wallCell[k].second=i;
	    cellWall[i].push_back(k);
	  }
	}
	if( !inWallFlag ) {
	  //New wall, set one cell and two vertices
	  cellWall[i].push_back(wallCell.size());
	  std::pair<size_t,size_t> tmpPair(i,i);
	  wallCell.push_back(tmpPair);
	  vertexWall[vI1].push_back(wallVertex.size());
	  vertexWall[vI2].push_back(wallVertex.size());
	  tmpPair.first=vI1;
	  tmpPair.second=vI2;
	  wallVertex.push_back(tmpPair);
	}
	else if( inWallFlag>1 ) {
	  std::cerr << "Tissue::createInitFromVoronoi() Vertices in"
		    << " multiple (>2) walls." << std::endl;
	  exit(-1);
	}
      }
    }
  }//for i (adding walls)
  std::cerr << "New walls added\n";
  assert( cellWall.size() == cellVertex.size() );
  assert( wallCell.size() == wallVertex.size() );
  assert( vertexCell.size() == vertexWall.size() );
  assert( vertexCell.size() == vertexPos.size() );
  
  // 	std::cerr << "cellVertex and cellWall:" << std::endl;
  // 	for( size_t i=0 ; i<cellVertex.size() ; ++i ) {
  // 		std::cerr << i << " (" << cellVertex[i].size() << ") ";
  // 		for( size_t k=0 ; k<cellVertex[i].size() ; ++k )
  // 			std::cerr << cellVertex[i][k] << " ";
  // 		std::cerr << " (" << cellWall[i].size() << ") ";		
  // 		for( size_t k=0 ; k<cellWall[i].size() ; ++k )
  // 			std::cerr << cellWall[i][k] << " ";
  // 		std::cerr << std::endl;
  // 	}
  // 	std::cerr << "vertexCell and vertexWall:" << std::endl;
  // 	for( size_t i=0 ; i<vertexCell.size() ; ++i ) {
  // 		std::cerr << i << " (" << vertexCell[i].size() << ") "; 
  // 		for( size_t k=0 ; k<vertexCell[i].size() ; ++k )
  // 			std::cerr << vertexCell[i][k] << " ";
  // 		std::cerr << " (" << vertexWall[i].size() << ") ";
  // 		for( size_t k=0 ; k<vertexWall[i].size() ; ++k )
  // 			std::cerr << vertexWall[i][k] << " ";
  // 		std::cerr << std::endl;
  // 	}
  // 	std::cerr << "wallCell and wallVertex:" << std::endl;
  // 	for( size_t i=0 ; i<wallCell.size() ; ++i ) {
  // 		std::cerr << i << "  "; 
  // 		std::cerr << wallCell[i].first << " " << wallCell[i].second 
  // 							<< "  ";
  // 		std::cerr << wallVertex[i].first << " " << wallVertex[i].second;
  // 		std::cerr << std::endl;
  // 	}
  
  //Extract possible positions for the new vertices
  for( size_t i=numOldVertex ; i<vertexPos.size() ; ++i ) {
    std::cerr << i << std::endl;
    size_t foundWallFlag=0;
    size_t falseCell1=0,falseCell2=0;
    size_t wallI=0,cellI=0;
    //Identify wall
    for( size_t k=0 ; k<vertexWall[i].size() ; ++k )
      if( wallCell[ vertexWall[i][k] ].first != boundaryIndex &&
	  wallCell[ vertexWall[i][k] ].second != boundaryIndex ) {
	++foundWallFlag;
	wallI = vertexWall[i][k];
	falseCell1 = wallCell[ vertexWall[i][k] ].first;
	falseCell2 = wallCell[ vertexWall[i][k] ].second;
      }
    if( foundWallFlag != 1 ) {
      std::cerr << "Tissue::createInitFromVoronoi() Multiple ("
		<< foundWallFlag << ") walls found for vertex" 
		<< std::endl;
      exit(-1);
    }
    //Extract vertex at other end
    assert( wallVertex[wallI].first==i || wallVertex[wallI].second==i );
    size_t vertexI = wallVertex[wallI].first==i ? 
      wallVertex[wallI].second : wallVertex[wallI].first; 
    
    //Find cell not connected to previous wall
    size_t foundCellFlag=0;
    for( size_t k=0 ; k<vertexCell[vertexI].size() ; ++k )
      if( vertexCell[vertexI][k] != falseCell1 && 
	  vertexCell[vertexI][k] != falseCell2 ) {
	++foundCellFlag;				
	cellI = vertexCell[vertexI][k];
      }
    if( foundCellFlag != 1 ) {
      std::cerr << "Tissue::createInitFromVoronoi() Multiple ("
		<< foundCellFlag << ") cells found for vertex " 
		<< vertexI << std::endl;
      std::cerr << falseCell1 << " " << falseCell2 << std::endl;
      exit(-1);
    }
    std::vector<double> cellPos( vertexPos[0].size() );
    if( !cellVertex[cellI].size() ) {
      std::cerr << "No vertices defined for chosen cell" << std::endl;
      exit(-1);
    }
    for( size_t k=0 ; k<cellVertex[cellI].size() ; ++k )
      for( size_t d=0 ; d<cellPos.size() ; ++d )
	cellPos[d] += vertexPos[cellVertex[cellI][k]][d];
    for( size_t d=0 ; d<cellPos.size() ; ++d )
      cellPos[d] /= cellVertex[cellI].size();
    
    //Extract direction and normalize
    std::vector<double> direction(cellPos.size());
    double norm=0.0;
    for( size_t d=0 ; d<cellPos.size() ; ++d ) {
      direction[d] = vertexPos[vertexI][d]-cellPos[d];
      norm += direction[d]*direction[d];
    }
    if( norm<=0.0 ) {
      std::cerr << "Tissue::createInitFromVoronoi() Direction "
		<< "without length (" << norm << ")" << std::endl; 
      exit(-1);
    }
    norm = std::sqrt(norm);
    if( norm<=0.0 ) {
      std::cerr << "Tissue::createInitFromVoronoi() Direction "
		<< "without length (" << norm << ")" << std::endl; 
      exit(-1);
    }
    for( size_t d=0 ; d<cellPos.size() ; ++d )
      direction[d] /= norm;
    //Set new vertex position
    double length=1.0;
    for( size_t d=0 ; d<cellPos.size() ; ++d )
      vertexPos[i][d] = vertexPos[vertexI][d] + length*direction[d];
  }
  
  //Create the tissue
  assert( cellWall.size() == cellVertex.size() );
  assert( wallCell.size() == wallVertex.size() );
  assert( vertexCell.size() == vertexWall.size() );
  assert( vertexCell.size() == vertexPos.size() );
  setNumCell( cellWall.size() );
  setNumWall( wallCell.size() );
  setNumVertex( vertexCell.size() );
  for( size_t i=0 ; i<numCell() ; ++i ) {
    cell(i).setIndex(i);
    for( size_t k=0 ; k<cellWall[i].size() ; ++k )
      cell(i).addWall( &wall(cellWall[i][k]) );
    for( size_t k=0 ; k<cellVertex[i].size() ; ++k )
      cell(i).addVertex( &vertex(cellVertex[i][k]) );
  }
  for( size_t i=0 ; i<numVertex() ; ++i ) {
    vertex(i).setIndex(i);
    for( size_t k=0 ; k<vertexCell[i].size() ; ++k )
      vertex(i).addCell( &cell(vertexCell[i][k]) );
    for( size_t k=0 ; k<vertexWall[i].size() ; ++k )
      vertex(i).addWall( &wall(vertexWall[i][k]) );
  }
  for( size_t i=0 ; i<numWall() ; ++i ) {
    wall(i).setIndex(i);
    if( wallCell[i].first == static_cast<size_t>(-1) )
      wall(i).setCell(background(),&cell(wallCell[i].second));
    else if( wallCell[i].second == static_cast<size_t>(-1) )
      wall(i).setCell(&cell(wallCell[i].first),background());
    else
      wall(i).setCell(&cell(wallCell[i].first),&cell(wallCell[i].second));
    wall(i).setVertex( &vertex(wallVertex[i].first),&vertex(wallVertex[i].second) );
  }
  
  // 	std::cerr << "cellVertex and cellWall:" << std::endl;
  // 	for( size_t i=0 ; i<cellVertex.size() ; ++i ) {
  // 		std::cerr << i << "," << cell(i).index() << " (" << cellVertex[i].size() << ","
  // 							<< cell(i).numVertex() << ") ";
  // 		for( size_t k=0 ; k<cellVertex[i].size() ; ++k )
  // 			std::cerr << cellVertex[i][k] << "," << cell(i).vertex(k)->index()
  // 								<< " ";
  // 		std::cerr << " (" << cellWall[i].size() << ","
  // 							<< cell(i).numWall() << ") ";
  // 		for( size_t k=0 ; k<cellWall[i].size() ; ++k )
  // 			std::cerr << cellWall[i][k] << "," << cell(i).wall(k)->index()
  // 								<< " ";
  // 		std::cerr << std::endl;
  // 	}
  // 	std::cerr << "vertexCell and vertexWall:" << std::endl;
  // 	for( size_t i=0 ; i<vertexCell.size() ; ++i ) {
  // 		std::cerr << i << "," << vertex(i).index() << " (" 
  // 							<< vertexCell[i].size() << "," 
  // 							<< vertex(i).numCell() << ") ";
  // 		for( size_t k=0 ; k<vertexCell[i].size() ; ++k )
  // 			std::cerr << vertexCell[i][k] << "," << vertex(i).cell(k)->index()
  // 								<< " ";
  // 		std::cerr << " (" << vertexWall[i].size() << "," 
  // 							<< vertex(i).numWall() << ") ";
  // 		for( size_t k=0 ; k<vertexWall[i].size() ; ++k )
  // 			std::cerr << vertexWall[i][k] << "," << vertex(i).wall(k)->index()
  // 								<< " ";
  // 		std::cerr << std::endl;
  // 	}
  // 	std::cerr << "wallCell and wallVertex:" << std::endl;
  // 	for( size_t i=0 ; i<wallCell.size() ; ++i ) {
  // 		std::cerr << i << "," << wall(i).index() << "  "; 
  // 		std::cerr << wallCell[i].first << "," << wall(i).cell1()->index()
  // 							<< " " << wallCell[i].second << "," 
  // 							<< wall(i).cell2()->index() << "  ";
  // 		std::cerr << wallVertex[i].first << "," << wall(i).vertex1()->index()
  // 							<< " " << wallVertex[i].second << "," 
  // 							<< wall(i).vertex2()->index() << "  ";
  
  // 		std::cerr << std::endl;
  // 	}
  
  //Set positions for the vertices
  assert( numVertex()==vertexPos.size() );
  for( size_t i=0 ; i<numVertex() ; ++i )
    vertex(i).setPosition(vertexPos[i]);
  
  //Get the wall lengths from the vertex positions
  setWallLengthFromVertexPosition();
  
  std::cerr << "Checking tissue" << std::endl;
  checkConnectivity(verbose);
  std::cerr << "Tissue created" << std::endl;
  DataMatrix cellData( numCell() ),
    cellDeriv( numCell() ),wallData( numWall() ),wallDeriv( numWall() ),
    vertexData( numVertex() ),vertexDeriv( numVertex() );
  removeEpidermalCells(cellData,wallData,vertexData,cellDeriv,wallDeriv,
		       vertexDeriv);
  std::cerr << "Checking tissue after first removal" << std::endl;
  checkConnectivity(verbose);
  std::cerr << "Tisue created" << std::endl;
  removeEpidermalCells(cellData,wallData,vertexData,cellDeriv,wallDeriv,
		       vertexDeriv);
  std::cerr << "Checking tissue after second removal" << std::endl;
  checkConnectivity(verbose);
  std::cerr << "Tisue created" << std::endl;
  removeEpidermalCells(cellData,wallData,vertexData,cellDeriv,wallDeriv,
		       vertexDeriv);
  std::cerr << "Checking tissue after third removal" << std::endl;
  checkConnectivity(verbose);
  std::cerr << "Tisue created" << std::endl;
  
}

void Tissue::copyState(DataMatrix &cellData, DataMatrix &wallData, DataMatrix &vertexData)
{
  if (cellData.size()!=numCell()) {
    std::cerr << "Tissue::copyState() Not the same number of cells in Tissue as in data." << std::endl;
    exit(EXIT_FAILURE);
  }
  if (wallData.size()!=numWall()) {
    std::cerr << "Tissue::copyState() Not the same number of walls in Tissue as in data." << std::endl;
    exit(EXIT_FAILURE);
  }
  if (vertexData.size()!=numVertex()) {
    std::cerr << "Tissue::copyState() Not the same number of vertices in Tissue as in data." << std::endl;
    exit(EXIT_FAILURE);
  }
  for (size_t i=0; i<numCell(); ++i) {
    assert( cell(i).numVariable()==cellData[i].size() );
    cell(i).setVariable(cellData[i]);
  }
  for (size_t i=0; i<numWall(); ++i) {
    wall(i).setLength(wallData[i][0]);
    assert( wall(i).numVariable()==wallData[i].size()-1 ); //wallData also stores length
    for (size_t j=0; j<wall(i).numVariable(); ++j)
      wall(i).setVariable(j,wallData[i][j+1]);
  }  
  for (size_t i=0; i<numVertex(); ++i) {
    vertex(i).setPosition(vertexData[i]);
  }
}

void Tissue::derivs( DataMatrix &cellData,
		     DataMatrix &wallData,
		     DataMatrix &vertexData,
		     DataMatrix &cellDeriv,
		     DataMatrix &wallDeriv,
		     DataMatrix &vertexDeriv ) 
{  
  //Set all derivatives to zero
  for( size_t i=0 ; i<cellDeriv.size() ; ++i )
    std::fill(cellDeriv[i].begin(),cellDeriv[i].end(),0.0);
  for( size_t i=0 ; i<wallDeriv.size() ; ++i )
    std::fill(wallDeriv[i].begin(),wallDeriv[i].end(),0.0);
  for( size_t i=0 ; i<vertexDeriv.size() ; ++i )
    std::fill(vertexDeriv[i].begin(),vertexDeriv[i].end(),0.0);
  
  //Calculate derivative contributions from all reactions
  for( size_t r=0 ; r<numReaction() ; ++r )
    reaction(r)->derivs(*this,cellData,wallData,vertexData,
			cellDeriv,wallDeriv,vertexDeriv);
}

void Tissue::derivsWithAbs( DataMatrix &cellData,
			    DataMatrix &wallData,
			    DataMatrix &vertexData,
			    DataMatrix &cellDeriv,
			    DataMatrix &wallDeriv,
			    DataMatrix &vertexDeriv,
			    DataMatrix &sdydtCell,
			    DataMatrix &sdydtWall,
			    DataMatrix &sdydtVertex ) 
{  
  //Set all derivatives to zero
  for( size_t i=0 ; i<cellDeriv.size() ; ++i ) {
    std::fill(cellDeriv[i].begin(),cellDeriv[i].end(),0.0);
    std::fill(sdydtCell[i].begin(),sdydtCell[i].end(),0.0);
  }
  for( size_t i=0 ; i<wallDeriv.size() ; ++i ) {
    std::fill(wallDeriv[i].begin(),wallDeriv[i].end(),0.0);
    std::fill(sdydtWall[i].begin(),sdydtWall[i].end(),0.0);
  }
  for( size_t i=0 ; i<vertexDeriv.size() ; ++i ) {
    std::fill(vertexDeriv[i].begin(),vertexDeriv[i].end(),0.0);
    std::fill(sdydtVertex[i].begin(),sdydtVertex[i].end(),0.0);
  }  

  //Calculate derivative contributions from all reactions
  for( size_t r=0 ; r<numReaction() ; ++r ){
    reaction(r)->derivsWithAbs(*this,cellData,wallData,vertexData,
			       cellDeriv,wallDeriv,vertexDeriv,
			       sdydtCell,sdydtWall,sdydtVertex);}
}

void::Tissue::initiateReactions(DataMatrix &cellData,
				DataMatrix &wallData,
				DataMatrix &vertexData,
				DataMatrix &cellDeriv,
				DataMatrix &wallDeriv,
				DataMatrix &vertexDeriv ) 
{
  for (size_t i=0; i<numReaction(); ++i) {
    reaction(i)->initiate(*this,cellData,wallData,vertexData,cellDeriv,wallDeriv,vertexDeriv);
  }
}

void::Tissue::updateReactions(DataMatrix &cellData,
			      DataMatrix &wallData,
			      DataMatrix &vertexData,
			      double step) 
{
  for (size_t i=0; i<numReaction(); ++i)
    reaction(i)->update(*this,cellData,wallData,vertexData,step);	
}

void::Tissue::
initiateDirection(DataMatrix &cellData,
		  DataMatrix &wallData,
		  DataMatrix &vertexData,
		  DataMatrix &cellDerivs,
		  DataMatrix &wallDerivs,
		  DataMatrix &vertexDerivs ) 
{
  direction()->initiate(*this,cellData,wallData,vertexData,cellDerivs,
			wallDerivs,vertexDerivs);
}

void::Tissue::
updateDirection(double step,
		DataMatrix &cellData,
		DataMatrix &wallData,
		DataMatrix &vertexData,
		DataMatrix &cellDerivs,
		DataMatrix &wallDerivs,
		DataMatrix &vertexDerivs) 
{
  direction()->update(*this,step,cellData,wallData,vertexData,cellDerivs,
		      wallDerivs,vertexDerivs);	
}

void::Tissue::
updateDirectionDivision(size_t cellI,
			DataMatrix &cellData,
			DataMatrix &wallData,
			DataMatrix &vertexData,
			DataMatrix &cellDerivs,
			DataMatrix &wallDerivs,
			DataMatrix &vertexDerivs) 
{
  direction()->divide(*this,cellI,cellData,wallData,vertexData,
		      cellDerivs,wallDerivs,vertexDerivs);	
}

void Tissue::
checkCompartmentChange( DataMatrix &cellData,
			DataMatrix &wallData,
			DataMatrix &vertexData,
			DataMatrix &cellDeriv,
			DataMatrix &wallDeriv,
			DataMatrix &vertexDeriv ) {
  
  unsigned int uglyHackCounter = 0;
  
  for( size_t l=0 ; l<numCompartmentChange() ; ++l ) {
    for( size_t i=0 ; i<numCell() ; ++i ) {
      ++uglyHackCounter;
      
      if (uglyHackCounter > 1000000) {
	// Time to bail out.
	std::cerr << "Ugly hack counter lager than a million!\n";
	std::exit(EXIT_FAILURE);
      }
      
      if( compartmentChange(l)->flag(this,i,cellData,wallData,vertexData,cellDeriv,wallDeriv,vertexDeriv) ) {
	compartmentChange(l)->update(this,i,cellData,wallData,vertexData,cellDeriv,wallDeriv,vertexDeriv);
	//If cell division, sort walls and vertices for cell plus 
	//divided cell plus their neighbors
	//Get list of potential cells to be sorted
	//Also add division rule for directions
	if( compartmentChange(l)->numChange()==1 ) {
	  std::set<size_t> sortCell;
	  sortCell.insert(i);
	  size_t ii=numCell()-1;
	  sortCell.insert(ii);
	  for( size_t k=0 ; k<cell(i).numWall() ; ++k ) {
	    if( cell(i).wall(k)->cell1()->index() == i )
	      sortCell.insert(cell(i).wall(k)->cell2()->index());
	    else
	      sortCell.insert(cell(i).wall(k)->cell1()->index());
	  }
	  for( size_t k=0 ; k<cell(ii).numWall() ; ++k ) {
	    if( cell(ii).wall(k)->cell1()->index() == ii )
	      sortCell.insert(cell(ii).wall(k)->cell2()->index());
	    else
	      sortCell.insert(cell(ii).wall(k)->cell1()->index());
	  }									
	  //Remove if background within the list
	  sortCell.erase( static_cast<size_t>(-1) );
	  //Sort the cells
	  //for( std::set<size_t>::iterator k=sortCell.begin() ; 
	  //	 k!=sortCell.end() ; ++k )
	  //std::cerr << *k << " ";
	  //std::cerr << "to be sorted" << std::endl;
	  
	  // If one of the daughter
	  // cells is on the edge and
	  // only has the other daughter
	  // cell as its neighbor the
	  // other daughter cell needs
	  // to be sorted
	  // first. Therefore we save
	  // all cells with only one
	  // neighbor in a second set of
	  // indexes and sort them in a
	  // second round.
	  std::set<size_t> oneNeighborCells;
	  
	  for (std::set<size_t>::iterator k = sortCell.begin(); k != sortCell.end(); ++k)
	    {
	      Cell &cellToSort = cell(*k);
	      
	      int counter = 0;
	      
	      for (size_t wallIndex = 0; wallIndex < cellToSort.numWall(); ++wallIndex)
		{
		  if (cellToSort.cellNeighbor(wallIndex) != background())
		    {
		      ++counter;
		    }
		}
	      
	      if (counter == 1)
		{
		  oneNeighborCells.insert(cellToSort.index());
		}
	      else
		{
		  cellToSort.sortWallAndVertex(*this);
		}
	    }
	  
	  for (std::set<size_t>::iterator k = oneNeighborCells.begin(); k != oneNeighborCells.end(); ++k)
	    {
	      Cell &cellToSort = cell(*k);
	      
	      cellToSort.sortWallAndVertex(*this);
	    }
	}	
	else if( compartmentChange(l)->numChange()==-1 )
	  --i;
	else if( compartmentChange(l)->numChange()<-1 )
	  i=numCell()+1;
      }
    }
  }
}

void Tissue::removeCell(size_t cellIndex,
                        DataMatrix &cellData,
                        DataMatrix &wallData,
                        DataMatrix &vertexData,
                        DataMatrix &cellDeriv,
                        DataMatrix &wallDeriv,
                        DataMatrix &vertexDeriv ) 
{
  assert(cellIndex<numCell());
  std::vector<size_t> wallRemove;
  //Mark walls for removal via index or change wallCell to background
  //To be removed if connected to removed cell(by default) and background
  for( size_t k=0 ; k<cell(cellIndex).numWall() ; ++k )
    if( ( cell(cellIndex).wall(k)->cell1()->index() == cellIndex &&
          cell(cellIndex).wall(k)->cell2() == background() ) ||
        ( cell(cellIndex).wall(k)->cell2()->index() == cellIndex &&
          cell(cellIndex).wall(k)->cell1() == background() ) ) {
      wallRemove.push_back(cell(cellIndex).wall(k)->index());
      //cell(cellIndex).wall(k)->setCell1( background() );
      //cell(cellIndex).wall(k)->setCell2( background() );
    }
    else if( cell(cellIndex).wall(k)->cell1()->index() == cellIndex ) {
      cell(cellIndex).wall(k)->setCell1( background() );
      //std::cerr << "Cell " << cellIndex << " switched to bg for wall "
      //				<< cell(cellIndex).wall(k)->index() << std::endl;
    }
    else if( cell(cellIndex).wall(k)->cell2()->index() == cellIndex ) {
      cell(cellIndex).wall(k)->setCell2( background() );
      //std::cerr << "Cell " << cellIndex << " switched to bg for wall "
      //				<< cell(cellIndex).wall(k)->index() << std::endl;
    }
    else {
      std::cerr << "Tissue::removeCell() wall not connected to cell"
                << std::endl;
      exit(-1);
    }
  //Remove cell and potential wall connections from vertices
  //Caveat: Also remove background...
  for( size_t k=0 ; k<cell(cellIndex).numVertex() ; ++k ) {
    //cell(cellIndex).vertex(k)->removeCell( &cell(cellIndex) );
    //if( cell(cellIndex).vertex(k)->removeCell( &cell(cellIndex) ) )
    //std::cerr << "Cell " << cellIndex << " removed from vertex "
    //				<< cell(cellIndex).vertex(k)->index() << std::endl;
    //cell(cellIndex).vertex(k)->removeCell( background() );
    //if( cell(cellIndex).vertex(k)->removeCell( background() ) )
    //std::cerr << "Background(cell) removed from vertex "
    //				<< cell(cellIndex).vertex(k)->index() << std::endl;
    
    //for( size_t w=0 ; w<wallRemove.size() ; ++w )
    //cell(cellIndex).vertex(k)->removeWall( &wall(wallRemove[w]) );
    
    //if( cell(cellIndex).vertex(k)->removeWall( &wall(wallRemove[w]) ) )
    //	std::cerr << "Wall " << wallRemove[w] << " removed from vertex "
    //						<< cell(cellIndex).vertex(k)->index() << std::endl;		
  }
  
  static size_t numCR=0,numWR=0,numVR=0; 
  //Remove vertices without connection to cells or walls
  for( size_t k=0 ; k<cell(cellIndex).numVertex() ; ++k ) {
    //Remove cell from vertex
    cell(cellIndex).vertex(k)->removeCell( &cell(cellIndex) );
    //Remove walls from vertex
    for( size_t w=0 ; w<wallRemove.size() ; ++w )
      cell(cellIndex).vertex(k)->removeWall( &wall(wallRemove[w]) );
    
    if( cell(cellIndex).vertex(k)->numCell() == 0 &&
        cell(cellIndex).vertex(k)->numWall() == 0 ) {
      //remove vertex
      size_t vI=cell(cellIndex).vertex(k)->index();
      if( vI>=vertexData.size() ) {
        std::cerr << "Tissue::removeCell() wrong in index " << std::endl
                  << numCell() << " " << numWall() << " " << numVertex()
                  << std::endl
                  << cellIndex << " " << cell(cellIndex).index() << " "
                  << cell(cellIndex).numVertex() << " "
                  << cell(cellIndex).numWall() << std::endl;
        for( size_t kk=0 ; kk<cell(cellIndex).numVertex() ; ++kk )
          std::cerr << cell(cellIndex).vertex(kk)->index() << " ";
        std::cerr << std::endl;				
      } 
      assert( vI<vertexData.size() );
      vertexData[vI] = vertexData[vertexData.size()-1];
      vertexDeriv[vI] = vertexDeriv[vertexDeriv.size()-1];
      vertexData.pop_back();
      vertexDeriv.pop_back();
      removeVertex(vI);
      std::cerr << "Vertex " << vI << " removed" << std::endl;
      numVR++;
    }
    else if( cell(cellIndex).vertex(k)->numCell() == 0 ||
             cell(cellIndex).vertex(k)->numWall() == 0 ) {
      std::cerr << "Tissue::removeCell() strange vertex." << std::endl;
      std::cerr << "It has " << cell(cellIndex).vertex(k)->numCell() 
                << " cells and " << cell(cellIndex).vertex(k)->numWall()
                << " walls." << std::endl;
      std::cerr << "Cells: ";
      for( size_t kk=0 ; kk<cell(cellIndex).vertex(k)->numCell() ; ++kk )
        std::cerr << cell(cellIndex).vertex(k)->cell(kk)->index() << " ";
      std::cerr << "\nWalls: ";
      for( size_t kk=0 ; kk<cell(cellIndex).vertex(k)->numWall() ; ++kk )
        std::cerr << cell(cellIndex).vertex(k)->wall(kk)->index() << " ";			
      exit(-1);
    }
  }
  //Remove walls connected to cellIndex and background
  for( size_t k=0 ; k<cell(cellIndex).numWall() ; ++k ) {
    if( ( cell(cellIndex).wall(k)->cell1() == background() &&
          cell(cellIndex).wall(k)->cell2() == &cell(cellIndex) ) ||
        ( cell(cellIndex).wall(k)->cell2() == background() &&
          cell(cellIndex).wall(k)->cell1() == &cell(cellIndex) ) ) {
      size_t wI=cell(cellIndex).wall(k)->index();
      //std::cerr << wI << " " << numWall() << " " << wallData.size() << std::endl;
      assert( wI<wallData.size() );
      wallData[wI] = wallData[wallData.size()-1];
      wallDeriv[wI] = wallDeriv[wallDeriv.size()-1];
      wallData.pop_back();
      wallDeriv.pop_back();
      removeWall(wI);
      //std::cerr << "Wall " << wI << " removed." << std::endl;
      numWR++;
      //wall(wI).setIndex(wI);
      //std::cerr << wI << " " << numWall() << " " << wallData.size() << std::endl;
    }
  }
  // 	//Old malfunctional version
  // 	for( size_t k=0 ; k<wallRemove.size() ; ++k ) {
  // 		size_t wI=wallRemove[k];
  // 		std::cerr << wI << " " << numWall() << " " << wallData.size() << std::endl;
  // 		assert( wI<wallData.size() );
  // 		wallData[wI] = wallData[wallData.size()-1];
  // 		wallDeriv[wI] = wallDeriv[wallDeriv.size()-1];
  // 		wallData.pop_back();
  // 		wallDeriv.pop_back();
  // 		removeWall(wI);
  // 		std::cerr << "Wall " << wI << " removed." << std::endl;
  // 		numWR++;
  // 		wall(wI).setIndex(wI);
  // 	}
  //Remove cell
  //std::cerr << cellIndex << " " << numCell() << " " << cellData.size()
  //				<< std::endl;
  assert( cellIndex<cellData.size() );
  cellData[cellIndex] = cellData[cellData.size()-1];
  cellDeriv[cellIndex] = cellDeriv[cellDeriv.size()-1];
  cellData.pop_back();
  cellDeriv.pop_back();
  removeCell(cellIndex);
  //std::cerr << "Cell " << cellIndex << " removed." << std::endl;
  //cell(cellIndex).setIndex(cellIndex);
  numCR++;
  //std::cerr << cellIndex << " " << numCell() << " " << cellData.size()
  //				<< std::endl;
  
  assert( cellData.size() == numCell() );
  assert( wallData.size() == numWall() );
  assert( vertexData.size() == numVertex() );	
  //checkConnectivity(1);
  std::cerr << numCR << " cells, " << numWR << " walls, and "
            << numVR << " vertices removed in total" << std::endl;
}

void Tissue::
removeCells(std::vector<size_t> &cellIndex,
            DataMatrix &cellData,
            DataMatrix &wallData,
            DataMatrix &vertexData,
            DataMatrix &cellDeriv,
            DataMatrix &wallDeriv,
            DataMatrix &vertexDeriv ) 
{
  // Sort the indices to make sure highest indices are removed first 
  // (since removed index is occupied with the last one)
  sort(cellIndex.begin(),cellIndex.end());
  
  size_t numRemove = cellIndex.size();
  for (size_t ii=0; ii<numRemove; ++ii) {
    size_t i = numRemove-(ii+1);
    removeCell(cellIndex[i],cellData,wallData,vertexData,cellDeriv,wallDeriv,
               vertexDeriv);
  }
}

void Tissue::
removeEpidermalCells(DataMatrix &cellData,
                     DataMatrix &wallData,
                     DataMatrix &vertexData,
                     DataMatrix &cellDeriv,
                     DataMatrix &wallDeriv,
                     DataMatrix &vertexDeriv,
                     double radialThreshold,
                     const bool checkBackground) 
{
  size_t dimension=vertexData[0].size();
  std::vector<size_t> cellR;
  //Mark cells for removal (sorted with highest index first
  for( size_t i=0 ; i<numCell() ; ++i ) {
    size_t cellI=numCell()-1-i;
    if (!checkBackground || cell(cellI).isNeighbor(background())) {
      if( radialThreshold>0.0 ) {//check that cell is outside
        std::vector<double> cellPos;
        cellPos = cell(cellI).positionFromVertex(vertexData);
        double r=0.0;
        for( size_t d=0 ; d<dimension ; ++d )
          r += cellPos[d]*cellPos[d];
        if( r>0.0 )
          r = std::sqrt(r);
        if( r>0.0 && r>radialThreshold )					
          cellR.push_back( cellI );
      }
      else
        cellR.push_back( cellI );
    }
  }
  
  if (cellR.size() > 0) {
    std::cerr << "Removing " << cellR.size() << " epidermal cells:\n" ;
    
    for (size_t i = 0; i < cellR.size(); ++i) {
      std::cerr << cellR[i] << " ";
    }
    
    std::cerr << "\n";
  }
  
  // Remove cells
  for (size_t i = 0; i < cellR.size(); ++i) {
    removeCell(cellR[i], cellData, wallData, vertexData, cellDeriv, wallDeriv, vertexDeriv);
    //checkConnectivity(1);
  }
}

void Tissue::removeEpidermalCellsMk2(DataMatrix &cellData,
                                     DataMatrix &wallData,
                                     DataMatrix &vertexData,
                                     DataMatrix &cellDeriv,
                                     DataMatrix &wallDeriv,
                                     DataMatrix &vertexDeriv,
                                     double radialThreshold) 
{
  size_t dimensions = vertexData[0].size();
  std::vector<size_t> cellR;
  
  //Mark cells for removal (sorted with highest index first)
  for (size_t i = 0; i < numCell(); ++i) {
    size_t cellI = numCell() - 1 - i;
    
    if (cell(cellI).isNeighbor(background())) {
      if (radialThreshold > 0.0) { //check that cell is outside
        Cell &c = cell(cellI);
        
        bool marked = true;
        
        for (size_t j = 0; j < c.numVertex(); ++j) {
          Vertex *vertex = c.vertex(j);
          
          double r = 0.0;
          
          for (size_t d = 0; d < dimensions; ++d) {
            r += vertexData[vertex->index()][d] * vertexData[vertex->index()][d];
          }
          
          if (r < radialThreshold * radialThreshold) {
            marked = false;
            break;
          }
        }
        
        if (marked == true) {
          cellR.push_back(cellI);
        }
      }
      else {
        cellR.push_back(cellI);
      }
    }
  }
  
  if (cellR.size() > 0) {
    std::cerr << "Removing " << cellR.size() << " epidermal cells:\n" ;
    
    for (size_t i = 0; i < cellR.size(); ++i) {
      std::cerr << cellR[i] << " ";
    }
    
    std::cerr << "\n";
  }
  
  // Remove cells
  for (size_t i = 0; i < cellR.size(); ++i) {
    removeCell(cellR[i], cellData, wallData, vertexData, cellDeriv, wallDeriv, vertexDeriv);
    //checkConnectivity(1);
  }
}

void Tissue::
removeEpidermalCellsAtDistance(DataMatrix &cellData,
                               DataMatrix &wallData,
                               DataMatrix &vertexData,
                               DataMatrix &cellDeriv,
                               DataMatrix &wallDeriv,
                               DataMatrix &vertexDeriv,
                               double distanceThreshold,double max,
                               size_t direction ) 
{
  size_t dimension;
  dimension = vertexData[0].size();
  assert( direction<dimension );
  std::vector<size_t> cellR;
  //Mark cells for removal (sorted with highest index first
  for( size_t i=0 ; i<numCell() ; ++i ) {
    size_t cellI=numCell()-1-i;
    if( cell(cellI).isNeighbor( background() ) ) {
      std::vector<double> cellPos;
      cellPos = cell(cellI).positionFromVertex(vertexData);
      double dist = std::fabs( cellPos[direction]-max );
      if( dist>distanceThreshold )					
        cellR.push_back( cellI );
    }
  }
  
  if (cellR.size() > 0) {
    std::cerr << "Removing " << cellR.size() << " epidermal cells:\n" ;
    
    for (size_t i = 0; i < cellR.size(); ++i) {
      std::cerr << cellR[i] << " ";
    }
    
    std::cerr << "\n";
  }
  
  // Remove cells
  for (size_t i = 0; i < cellR.size(); ++i) {
    removeCell(cellR[i], cellData, wallData, vertexData, cellDeriv, wallDeriv, vertexDeriv);
    //checkConnectivity(1);
  }
}

void Tissue::divideCell( Cell *divCell, size_t wI, size_t w3I, 
			 std::vector<double> &v1Pos,
			 std::vector<double> &v2Pos,
			 DataMatrix &cellData,
			 DataMatrix &wallData,
			 DataMatrix &vertexData,
			 DataMatrix &cellDeriv,
			 DataMatrix &wallDeriv,
			 DataMatrix &vertexDeriv,
			 std::vector<size_t> &volumeChangeList,
			 double threshold) 
  
{	
  size_t Nc=numCell(),Nw=numWall(),Nv=numVertex();
  size_t i = divCell->index();
  size_t dimension = vertexData[0].size();
  
  //Move new vertices if closer than threshold to old vertex
  if( threshold>=0.0 ) {
    Wall *w1 = divCell->wall(wI), *w2 = divCell->wall(w3I);
    double w1L=0.0,w2L=0.0,t1=0.0,t2=0.0;
    for( size_t dim=0; dim<dimension; ++dim ) {
      w1L += (vertexData[w1->vertex1()->index()][dim]-
	      vertexData[w1->vertex2()->index()][dim]) *
	(vertexData[w1->vertex1()->index()][dim] -
	 vertexData[w1->vertex2()->index()][dim]);
      w2L += (vertexData[w2->vertex1()->index()][dim]-
	      vertexData[w2->vertex2()->index()][dim]) *
	(vertexData[w2->vertex1()->index()][dim] -
	 vertexData[w2->vertex2()->index()][dim]);
      t1 += (v1Pos[dim]-vertexData[w1->vertex2()->index()][dim])*
	(v1Pos[dim]-vertexData[w1->vertex2()->index()][dim]);
      t2 += (v2Pos[dim]-vertexData[w2->vertex2()->index()][dim])*
	(v2Pos[dim]-vertexData[w2->vertex2()->index()][dim]);
    }
    w1L = std::sqrt(w1L);
    w2L = std::sqrt(w2L);
    t1 = std::sqrt(t1)/w1L;
    t2 = std::sqrt(t2)/w2L;
    assert( t1>=0.0 && t1<=1.0 );
    assert( t2>=0.0 && t2<=1.0 );
    if( t1<threshold ) {
      std::cerr << "Tissue::divideCell() Moving vertex 1 from "
		<< t1 << " to " << threshold << std::endl;
      t1=threshold;
      for( size_t dim=0; dim<dimension; ++dim ) {
	v1Pos[dim] = vertexData[w1->vertex2()->index()][dim] +
	  t1*(vertexData[w1->vertex1()->index()][dim]-
	      vertexData[w1->vertex2()->index()][dim]);
      }
    }
    else if( t1>(1.0-threshold) ) {
      std::cerr << "Tissue::divideCell() Moving vertex 1 from "
		<< t1 << " to " << 1.0-threshold << std::endl;
      t1 = threshold;
      for( size_t dim=0; dim<dimension; ++dim ) {
	v1Pos[dim] = vertexData[w1->vertex1()->index()][dim] +
	  t1*(vertexData[w1->vertex2()->index()][dim]-
	      vertexData[w1->vertex1()->index()][dim]);
      }
    }
    if( t2<threshold ) {
      std::cerr << "Tissue::divideCell() Moving vertex 2 from "
		<< t2 << " to " << threshold << std::endl;
      t2=threshold;
      for( size_t dim=0; dim<dimension; ++dim ) {
	v2Pos[dim] = vertexData[w2->vertex2()->index()][dim] +
	  t2*(vertexData[w2->vertex1()->index()][dim]-
	      vertexData[w2->vertex2()->index()][dim]);
      }
    }
    else if( t2>(1.0-threshold) ) {
      std::cerr << "Tissue::divideCell() Moving vertex 2 from "
		<< t2 << " to " << 1.0-threshold << std::endl;
      t2 = threshold;
      for( size_t dim=0; dim<dimension; ++dim ) {
	v2Pos[dim] = vertexData[w2->vertex1()->index()][dim] +
	  t2*(vertexData[w2->vertex2()->index()][dim]-
	      vertexData[w2->vertex1()->index()][dim]);
      }
    }		
  }
  
  //Create the new data structure and set indices in the tissue vectors
  //
  //Add the new cell
  addCell( cell(i) );
  cell(Nc).setIndex(Nc);
  cellData.resize(Nc+1,cellData[i]);
  cellDeriv.resize(Nc+1,cellDeriv[0]);
  
  //Add the two new vertices
  Vertex tmpVertex;
  tmpVertex.setPosition(v1Pos);
  tmpVertex.setIndex(Nv);
  addVertex(tmpVertex);
  vertexData.resize(Nv+1,v1Pos);  
  tmpVertex.setPosition(v2Pos);
  tmpVertex.setIndex(Nv+1);
  addVertex(tmpVertex);
  vertexData.resize(Nv+2,v2Pos);  
  vertexDeriv.resize(Nv+2,vertexDeriv[0]);  
  
  //Add the three new walls
  Wall tmpWall;
  //New wall dividing old cell into two
  tmpWall.setIndex(Nw);
  addWall(tmpWall);
  wallData.resize(Nw+1,wallData[0]);
  double tmpLength = 0.0;
  for( size_t d=0 ; d<dimension ; ++d )
    tmpLength += (v1Pos[d]-v2Pos[d])*(v1Pos[d]-v2Pos[d]);
  wallData[Nw][0] = std::sqrt( tmpLength );
  
  //Wall continuing the first selected wall
  addWall( *(cell(i).wall(wI)) );
  wall(Nw+1).setIndex(Nw+1);
  //Set new lengths as fractions of the old determined from the new 
  //vertex position
  double oldL = wallData[cell(i).wall(wI)->index()][0];
  size_t v1w = cell(i).wall(wI)->vertex1()->index();
  size_t v2w = cell(i).wall(wI)->vertex2()->index();
  tmpLength=0.0; 
  double tmpLengthFrac=0.0; 
  for( size_t d=0 ; d<dimension ; ++d ) {
    tmpLengthFrac += (v1Pos[d]-vertexData[v1w][d])*
      (v1Pos[d]-vertexData[v1w][d]);
    tmpLength += (vertexData[v2w][d]-vertexData[v1w][d])*
      (vertexData[v2w][d]-vertexData[v1w][d]);
  }
  tmpLength = std::sqrt( tmpLength );
  tmpLengthFrac = std::sqrt( tmpLengthFrac );
  double lengthFrac = tmpLengthFrac/tmpLength;
  wallData.resize(Nw+2,wallData[cell(i).wall(wI)->index()]);
  wallData[cell(i).wall(wI)->index()][0] = lengthFrac*oldL;
  wallData[Nw+1][0] = oldL-wallData[cell(i).wall(wI)->index()][0];
  
  //Wall continuing the second selected wall
  addWall(*(cell(i).wall(w3I)));
  wall(Nw+2).setIndex(Nw+2);
  //Set new lengths as fractions of the old determined from the new
  //vertex position
  oldL = wallData[cell(i).wall(w3I)->index()][0];
  v1w = cell(i).wall(w3I)->vertex1()->index();
  v2w = cell(i).wall(w3I)->vertex2()->index();
  tmpLength=0.0; 
  tmpLengthFrac=0.0; 
  for( size_t d=0 ; d<dimension ; ++d ) {
    tmpLengthFrac += (v2Pos[d]-vertexData[v1w][d])*
      (v2Pos[d]-vertexData[v1w][d]);
    tmpLength += (vertexData[v2w][d]-vertexData[v1w][d])*
      (vertexData[v2w][d]-vertexData[v1w][d]);
  }
  tmpLength = std::sqrt( tmpLength );
  tmpLengthFrac = std::sqrt( tmpLengthFrac );
  lengthFrac = tmpLengthFrac/tmpLength;
  wallData.resize(Nw+3,wallData[cell(i).wall(w3I)->index()]);
  wallData[cell(i).wall(w3I)->index()][0] = lengthFrac*oldL;
  wallData[Nw+2][0] = oldL-wallData[cell(i).wall(w3I)->index()][0];
  
  //Resize derivative matrix as well
  wallDeriv.resize(Nw+3,wallDeriv[0]);
  
  // Create connection matrix by first selecting vertices and walls for
  // the new and old cells
  //
  std::vector<size_t> oldVIndex,newVIndex,oldWIndex,newWIndex,
    usedWIndex( cell(i).numWall() );
  
  //Extract walls and vertices for 'old' cell
  size_t tmpWIndex = cell(i).wall(wI)->index(); 
  size_t tmpVIndex = cell(i).wall(wI)->vertex1()->index();
  size_t nextW = wI;
  //size_t count=0;
  //std::cerr << "wI=" << wI << "(" << cell(i).wall(wI)->index() 
  //				<< ") w3I=" << w3I << "(" << cell(i).wall(w3I)->index() 
  //				<< ")" << std::endl;
  do {
    //std::cerr << count << " " << nextW << "  " << tmpWIndex << " " 
    //				<< tmpVIndex << std::endl;
    //add to new vectors
    oldWIndex.push_back( tmpWIndex );
    oldVIndex.push_back( tmpVIndex );
    //mark as used
    usedWIndex[nextW] = 1;
    //usedVIndex[xxx] = 1;
    //find next wall
    nextW = cell(i).numWall();
    size_t flag=0;
    for( size_t w=0 ; w<cell(i).numWall() ; ++w ) {
      if( !usedWIndex[w] && (cell(i).wall(w)->vertex1()->index()==tmpVIndex ||
			     cell(i).wall(w)->vertex2()->index()==tmpVIndex ) ) {
	nextW = w;
	flag++;
      }
    }
    if( flag != 1 ) {
      std::cerr << "Tissue::divideCell() " << flag  
		<< " walls marked for next wall..." << std::endl;
      std::cerr << tmpVIndex << std::endl;
      for( size_t w=0 ; w<cell(i).numWall() ; ++w ) {
	std::cerr << w << " " << usedWIndex[w] << " " 
		  << cell(i).wall(w)->vertex1()->index() << " "
		  << cell(i).wall(w)->vertex2()->index() << std::endl;
      }
    }
    assert( flag==1 );
    tmpWIndex = cell(i).wall(nextW)->index();
    if( cell(i).wall(nextW)->vertex1()->index()==tmpVIndex )
      tmpVIndex = cell(i).wall(nextW)->vertex2()->index();
    else if( cell(i).wall(nextW)->vertex2()->index()==tmpVIndex )
      tmpVIndex = cell(i).wall(nextW)->vertex1()->index();
    else {
      std::cerr << "Tissue::DivideCell() " 
		<< "Wrong vertex indices for chosen wall " 
		<< cell(i).wall(nextW)->vertex1()->index() << " "
		<< cell(i).wall(nextW)->vertex2()->index() << std::endl;
      exit(-1);
    }
    //std::cerr << count++ << " " << nextW << "  " << tmpWIndex << " " 
    //				<< tmpVIndex << std::endl;
  } while( nextW != w3I ); 
  if( cell(i).wall(nextW)->vertex1()->index()==oldVIndex[oldVIndex.size()-1] )
    oldWIndex.push_back(cell(i).wall(w3I)->index());
  else if( cell(i).wall(nextW)->vertex2()->index()==oldVIndex[oldVIndex.size()-1] )
    oldWIndex.push_back(Nw+2);
  else {
    std::cerr << "Wrong last index for old cell (not in w3I)" << std::endl;
    exit(-1);
  }
  //Extract walls and vertices for 'new' cell
  tmpWIndex = Nw+1;//new copy of wI; 
  tmpVIndex = wall(Nw+1).vertex2()->index();
  //std::cerr << "Old cell done..." << std::endl;
  usedWIndex[wI]=0;
  usedWIndex[w3I]=0;
  nextW=wI;
  do {
    //std::cerr << count << " " << nextW << "  " << tmpWIndex << " " 
    //				<< tmpVIndex << std::endl;
    //add to new vectors
    newWIndex.push_back( tmpWIndex );
    newVIndex.push_back( tmpVIndex );
    //mark as used
    usedWIndex[nextW] = 1;
    //usedVIndex[xxx] = 1;
    //find next wall
    nextW = cell(i).numWall();
    size_t flag=0;
    for( size_t w=0 ; w<cell(i).numWall() ; ++w ) {
      if( !usedWIndex[w] && ( cell(i).wall(w)->vertex1()->index()==tmpVIndex ||
			      cell(i).wall(w)->vertex2()->index()==tmpVIndex ) ) {
	nextW = w;
	flag++;
      }
    }
    if( flag != 1 )
      std::cerr << flag << " walls marked for next wall..." << std::endl;
    assert( flag==1 );
    tmpWIndex = cell(i).wall(nextW)->index();
    if( cell(i).wall(nextW)->vertex1()->index()==tmpVIndex )
      tmpVIndex = cell(i).wall(nextW)->vertex2()->index();
    else if( cell(i).wall(nextW)->vertex2()->index()==tmpVIndex )
      tmpVIndex = cell(i).wall(nextW)->vertex1()->index();
    else {
      std::cerr << "Tissue::DivideCell() " 
		<< "Wrong vertex indices for chosen wall " 
		<< cell(i).wall(nextW)->vertex1() << " "
		<< cell(i).wall(nextW)->vertex2() << std::endl;
      exit(-1);
    }
    //std::cerr << count++ << " " << nextW << "  " << tmpWIndex << " " 
    //				<< tmpVIndex << std::endl;
  } while( nextW != w3I ); 
  if( cell(i).wall(nextW)->vertex1()->index()==newVIndex[newVIndex.size()-1] )
    newWIndex.push_back(cell(i).wall(w3I)->index());
  else if( cell(i).wall(nextW)->vertex2()->index()==newVIndex[newVIndex.size()-1] )
    newWIndex.push_back(Nw+2);
  else {
    std::cerr << "Wrong last index for new cell (not in w3I)" << std::endl;
    exit(-1);
  }
  
  //std::cerr << "New cell done..." << std::endl;
  
  //Add new vertices
  oldVIndex.push_back(Nv);
  oldVIndex.push_back(Nv+1);
  newVIndex.push_back(Nv);
  newVIndex.push_back(Nv+1);
  
  //Add new walls
  oldWIndex.push_back(Nw);
  newWIndex.push_back(Nw);
  //newWIndex.push_back(Nw+1);
  //newWIndex.push_back(Nw+2);
  
  //std::cerr << "Cell " << i << " has " << cell(i).numWall() << " walls and " 
  //				<< cell(i).numVertex() << " vertices" << std::endl;
  // 	for( size_t w=0 ; w<cell(i).numWall() ; ++w )
  // 		std::cerr << cell(i).wall(w)->index() << " ";
  // 	std::cerr << std::endl;
  // 	for( size_t v=0 ; v<cell(i).numVertex() ; ++v )
  // 		std::cerr << cell(i).vertex(v)->index() << " ";
  // 	std::cerr << std::endl;
  
  // 	for( size_t k=0 ; k<oldWIndex.size() ; ++k )
  // 		std::cerr << oldWIndex[k] << " ";
  // 	std::cerr << "  ";
  // 	for( size_t k=0 ; k<newWIndex.size() ; ++k )
  // 		std::cerr << newWIndex[k] << " ";
  // 	std::cerr << std::endl;
  
  // 	for( size_t k=0 ; k<oldVIndex.size() ; ++k )
  // 		std::cerr << oldVIndex[k] << " ";
  // 	std::cerr << "  ";
  // 	for( size_t k=0 ; k<newVIndex.size() ; ++k )
  // 		std::cerr << newVIndex[k] << " ";
  // 	std::cerr << std::endl;
  //exit(0);
  
  //Set vertices and cells for the walls
  wall(Nw).setVertex( &(vertex(Nv)),&(vertex(Nv+1)) );
  wall(Nw).setCell( &(cell(i)),&(cell(Nc)) );
  assert( wall(Nw).cell1()->index()==i );
  assert( wall(Nw).cell2()->index()==Nc );
  
  wall(cell(i).wall(wI)->index()).setVertex2( &(vertex(Nv)) );
  int vInCellFlag=0;
  for( size_t k=0 ; k<oldVIndex.size() ; ++k )
    if( oldVIndex[k]==cell(i).wall(wI)->vertex1()->index() ) 
      vInCellFlag++;
  assert( vInCellFlag==1 );
  
  wall(Nw+1).setVertex1( &(vertex(Nv)) );
  vInCellFlag=0;
  for( size_t k=0 ; k<newVIndex.size() ; ++k )
    if( newVIndex[k]==wall(Nw+1).vertex2()->index() ) 
      vInCellFlag++;
  assert( vInCellFlag==1 );
  
  if( wall(Nw+1).cell1()->index() == i ) {
    wall(Nw+1).setCell1( &(cell(Nc)) );
    if( wall(Nw+1).cell2()->index() < Nc ) {
      wall(Nw+1).cell2()->addWall( &(wall(Nw+1)) );
      wall(Nw+1).cell2()->addVertex( &(vertex(Nv)) );
    }
  }
  else if( wall(Nw+1).cell2()->index() == i ) {
    wall(Nw+1).setCell2( &(cell(Nc)) );
    if( wall(Nw+1).cell1()->index() < Nc ) {
      wall(Nw+1).cell1()->addWall( &(wall(Nw+1)) );
      wall(Nw+1).cell1()->addVertex( &(vertex(Nv)) );
    }
  }
  else {
    std::cerr << "Tissue::divideCell() "
	      << "First wall not connected to dividing cell" << std::endl;
    std::cerr << i << " " << cell(i).index() << "\t" 
	      << wall( cell(i).wall(wI)->index() ).cell1()->index() << " "
	      << wall(Nw+1).cell1()->index() << "\t"
	      << wall(cell(i).wall(wI)->index()).cell2()->index() << " "
	      << wall(Nw+1).cell2()->index() << "\n";		
    exit(-1);
  }
  
  vInCellFlag=0;
  size_t wInOldCellFlag=0;
  for( size_t k=0 ; k<oldWIndex.size() ; ++k )
    if( cell(i).wall(w3I)->index()==oldWIndex[k] )
      wInOldCellFlag=1;
  if( wInOldCellFlag ) {
    for( size_t k=0 ; k<oldVIndex.size() ; ++k )
      if( oldVIndex[k]==wall(cell(i).wall(w3I)->index()).vertex1()->index() ) 
	vInCellFlag++;
    if( vInCellFlag == 1 ) {
      wall(cell(i).wall(w3I)->index()).setVertex2( &(vertex(Nv+1)) );
      wall(Nw+2).setVertex1( &(vertex(Nv+1)) );
    }
    else {
      wall(cell(i).wall(w3I)->index()).setVertex1( &(vertex(Nv+1)) );
      wall(Nw+2).setVertex2( &(vertex(Nv+1)) );
    }
  }
  else {		
    for( size_t k=0 ; k<oldVIndex.size() ; ++k )
      if( oldVIndex[k]==wall(Nw+2).vertex1()->index() ) 
	vInCellFlag++;
    if( vInCellFlag == 1 ) {
      wall(Nw+2).setVertex2( &(vertex(Nv+1)) );
      wall(cell(i).wall(w3I)->index()).setVertex1( &(vertex(Nv+1)) );
    }
    else {
      wall(Nw+2).setVertex1( &(vertex(Nv+1)) );
      wall(cell(i).wall(w3I)->index()).setVertex2( &(vertex(Nv+1)) );
    }
  }
  
  //Change cell connection for new wall (w3I or Nw+2)
  size_t newWallIndex = cell(i).wall(w3I)->index();
  size_t w3IInNewCellFlag = 0;
  for( size_t w=0 ; w<newWIndex.size() ; ++w )
    if( newWIndex[w]==newWallIndex )
      w3IInNewCellFlag++;
  if( !w3IInNewCellFlag++ )
    newWallIndex=Nw+2;
  
  if( wall(newWallIndex).cell1()->index() == i ) {
    wall(newWallIndex).setCell1( &(cell(Nc)) );
    if( wall(newWallIndex).cell2()->index() < Nc ) {
      wall(newWallIndex).cell2()->addWall( &(wall(Nw+2)) );
      wall(newWallIndex).cell2()->addVertex( &(vertex(Nv+1)) );
    }
  }    
  else if( wall(newWallIndex).cell2()->index() == i ) {
    wall(newWallIndex).setCell2( &(cell(Nc)) );
    if( wall(newWallIndex).cell1()->index() < Nc ) {
      wall(newWallIndex).cell1()->addWall( &(wall(Nw+2)) );
      wall(newWallIndex).cell1()->addVertex( &(vertex(Nv+1)) );
    }    
  }
  else {
    std::cerr << "Tissue::divideCell() "
	      << "Second wall not connected to dividing cell" << std::endl;
    for( size_t k=0 ; k<cell(i).numWall() ; ++k ) {
      size_t v1w3Itmp = cell(i).wall(k)->vertex1()->index();
      size_t v2w3Itmp = cell(i).wall(k)->vertex2()->index();
      std::cerr << vertexData[v1w3Itmp][0] << " " 
		<< vertexData[v1w3Itmp][1] << " 4\n"
		<< vertexData[v2w3Itmp][0] << " "
		<< vertexData[v2w3Itmp][1] << " 4\n\n\n";
    }
    exit(-1);
  }
  
  //Set cells and walls for the vertices
  //
  for( size_t v=0 ; v<newVIndex.size() ; ++v ) {
    std::vector<Cell*> tmpCell;
    std::vector<Wall*> tmpWall;
    if( newVIndex[v]==Nv ) {
      //First new vertex
      tmpCell.push_back( &(cell(i)) );
      tmpCell.push_back( &(cell(Nc)) );
      //plus one more
      if( cell(i).wall(wI)->cell1()->index()==i &&
	  cell(i).wall(wI)->cell2() != background() )
	tmpCell.push_back( cell(i).wall(wI)->cell2() );
      else if( cell(i).wall(wI)->cell2()->index()==i &&
	       cell(i).wall(wI)->cell1() != background() )
	tmpCell.push_back( cell(i).wall(wI)->cell1() );
      else if( cell(i).wall(wI)->cell1()->index() != i &&
	       cell(i).wall(wI)->cell2()->index() != i ){
	std::cerr << "Tissue::divideCell() "
		  << "Wall wI not connected to dividing cell" 
		  << std::endl;
	exit(-1);
      }
      tmpWall.push_back( cell(i).wall(wI) );
      tmpWall.push_back( &(wall(Nw)) );
      tmpWall.push_back( &(wall(Nw+1)) );
    }
    else if( newVIndex[v]==Nv+1 ) {
      //Second new vertex
      tmpCell.push_back( &(cell(i)) );
      tmpCell.push_back( &(cell(Nc)) );
      //plus one more
      if( cell(i).wall(w3I)->cell1()->index()==i ||
	  cell(i).wall(w3I)->cell1()->index()==Nc ) {
	if( cell(i).wall(w3I)->cell2() != background() )
	  tmpCell.push_back( cell(i).wall(w3I)->cell2() );
      }
      else if( cell(i).wall(w3I)->cell2()->index()==i ||
	       cell(i).wall(w3I)->cell2()->index()==Nc ) {
	if( cell(i).wall(w3I)->cell1() != background() )
	  tmpCell.push_back( cell(i).wall(w3I)->cell1() );
      }
      else {
	std::cerr << "Tissue::divideCell() "
		  << "Wall w3I not connected to dividing cell" 
		  << std::endl;
	exit(-1);
      }
      tmpWall.push_back( cell(i).wall(w3I) );
      tmpWall.push_back( &(wall(Nw)) );
      tmpWall.push_back( &(wall(Nw+2)) );
    }
    else {
      //For other vertices
      for( size_t c=0 ; c<vertex(newVIndex[v]).numCell() ; ++c )
	if( vertex(newVIndex[v]).cell(c)->index() == i )
	  vertex(newVIndex[v]).setCell(c,&(cell(Nc)));
      for( size_t w=0 ; w<vertex(newVIndex[v]).numWall() ; ++w ) {
	if( vertex(newVIndex[v]).wall(w)->index() == 
	    cell(i).wall(wI)->index() )
	  vertex(newVIndex[v]).setWall(w,&(wall(Nw+1)));
      	else if( vertex(newVIndex[v]).wall(w)->index() == 
		 cell(i).wall(w3I)->index() ) {
	  size_t w3IInNewFlag=0;
	  for( size_t ww=0 ; ww<newWIndex.size() ; ++ww )
	    if( newWIndex[ww] == cell(i).wall(w3I)->index() )
	      w3IInNewFlag++;
	  if( !w3IInNewFlag )
	    vertex(newVIndex[v]).setWall(w,&(wall(Nw+2)));
	}
      }
    }
    if( newVIndex[v]>=Nv ) {
      vertex(newVIndex[v]).setCell(tmpCell);
      vertex(newVIndex[v]).setWall(tmpWall);
    }
  }
  
  //Finally check if vertices in old cell connected to w3I should be changed
  //to Nw+2
  for( size_t v=0 ; v<oldVIndex.size() ; ++v ) {
    if( oldVIndex[v] < Nv ) {
      for( size_t w=0 ; w<vertex(oldVIndex[v]).numWall() ; ++w ) {
      	if( vertex(oldVIndex[v]).wall(w)->index() == 
	    cell(i).wall(w3I)->index() ) {
	  size_t w3IInOldFlag=0;
	  for( size_t ww=0 ; ww<oldWIndex.size() ; ++ww )
	    if( oldWIndex[ww] == cell(i).wall(w3I)->index() )
	      w3IInOldFlag++;
	  if( !w3IInOldFlag )
	    vertex(oldVIndex[v]).setWall(w,&(wall(Nw+2)));
	}
      }
    }
  }
  
  //Set walls and vertices for the cells
  //This should be done after the vertex and wall insertion not to
  //mess up wI and w3I
  std::vector<Wall*> tmpW(oldWIndex.size());
  std::vector<Vertex*> tmpV(oldVIndex.size());
  for( size_t w=0 ; w<oldWIndex.size() ; ++w )
    tmpW[w] = &(wall(oldWIndex[w]));
  for( size_t v=0 ; v<oldVIndex.size() ; ++v )
    tmpV[v] = &(vertex(oldVIndex[v]));
  cell(i).setWall(tmpW);
  cell(i).setVertex(tmpV);
  
  tmpW.resize( newWIndex.size() );
  tmpV.resize( newVIndex.size() );
  for( size_t w=0 ; w<newWIndex.size() ; ++w ) {
    tmpW[w] = &(wall(newWIndex[w]));
    if( tmpW[w]->index() != Nw && tmpW[w]->cell1()->index()==i )
      tmpW[w]->setCell1( &(cell(Nc)) );
    else if( tmpW[w]->index() != Nw && tmpW[w]->cell2()->index()==i )
      tmpW[w]->setCell2( &(cell(Nc)) );
  }
  for( size_t v=0 ; v<newVIndex.size() ; ++v )
    tmpV[v] = &(vertex(newVIndex[v]));
  cell(Nc).setWall(tmpW);
  cell(Nc).setVertex(tmpV);
  assert( wall(Nw).cell1()->index()==i );
  assert( wall(Nw).cell2()->index()==Nc );
  
  //Update directional wall vector if applicable
  updateDirectionDivision(i,cellData,wallData,vertexData,
			  cellDeriv,wallDeriv,vertexDeriv);
  
  // Update the volume dependent variables for each cell variable index 
  // given in volumeChangeList
  if (volumeChangeList.size()) {
    cell(i).sortWallAndVertex(*this);
    cell(Nc).sortWallAndVertex(*this);
    double Vi = cell(i).calculateVolume(vertexData);
    double Vn = cell(Nc).calculateVolume(vertexData);
    double fi = Vi/(Vi+Vn);
    double fn = Vn/(Vi+Vn);
    
    for (size_t k=0; k<volumeChangeList.size(); ++k) {
      cellData[i][volumeChangeList[k]] *= fi;
      cellData[Nc][volumeChangeList[k]] *= fn;
    }
  }		
  //checkConnectivity(1);


  
  //    cell(Nc).sortWallAndVertex(*this);


  
  // for( size_t cl=0 ; cl<Nc+1 ; cl++ ){
  //   cell(cl).sortWallAndVertex(*this);
  //   std::cerr<<"cell "<< cl<<std::endl;    
  //   for( size_t wl=0 ; wl<cell(cl).numWall() ; wl++ ){
  //     std::cerr<<"wall "<< wl<<"("<<cell(cl).wall(wl)-> index()<<")"
  //              <<" between nodes " 
  //              <<cell(cl).wall(wl)->vertex1()-> index()<<","
  //              <<cell(cl).wall(wl)->vertex2()-> index()
  //              <<std::endl;
  //   }
    
  // }
  // std::cerr<<std::endl;
  // for( size_t cl=0 ; cl<Nc+1 ; cl++ ){
  //   std::cerr<<"cell "<< cl<<" vertices    :";    
  //   for( size_t wl=0 ; wl<cell(cl).numWall() ; wl++ ){
  //     std::cerr << wl<<"("<<cell(cl).vertex(wl)-> index()<<")  ";
      
  //   }
  //   std::cerr<<std::endl;
  // }
  

}


void Tissue::branchCell( Cell *brCell, size_t wI, 
			 std::vector<double> &v1Pos,
                         std::vector<double> &v2Pos,
                         std::vector<double> &v3Pos,
                         std::vector<double> &v4Pos,
			 DataMatrix &cellData,
			 DataMatrix &wallData,
			 DataMatrix &vertexData,
			 DataMatrix &cellDeriv,
			 DataMatrix &wallDeriv,
			 DataMatrix &vertexDeriv,
                         double threshold)
  
{	
  size_t Nc=numCell(),Nw=numWall(),Nv=numVertex();
  size_t i = brCell->index();
  size_t dimension = vertexData[0].size();
  
  //Create the new data structure and set indices in the tissue vectors
  //
  //Add the new cell  
 
  addCell( cell(i) );
  cell(Nc).setIndex(Nc);
  cellData.resize(Nc+1,cellData[1]);
  cellDeriv.resize(Nc+1,cellDeriv[0]);

 
  
  //Add the four new vertices
  Vertex tmpVertex;
  tmpVertex.setPosition(v1Pos);
  tmpVertex.setIndex(Nv);
  addVertex(tmpVertex);
  vertexData.resize(Nv+1,v1Pos);
  
  tmpVertex.setPosition(v2Pos);
  tmpVertex.setIndex(Nv+1);
  addVertex(tmpVertex);
  vertexData.resize(Nv+2,v2Pos);  

  tmpVertex.setPosition(v3Pos);
  tmpVertex.setIndex(Nv+2);
  addVertex(tmpVertex);
  vertexData.resize(Nv+3,v3Pos);  
  
  tmpVertex.setPosition(v4Pos);
  tmpVertex.setIndex(Nv+3);
  addVertex(tmpVertex);
  vertexData.resize(Nv+4,v4Pos);  
  
  //Resize derivative matrix for vertices    
  vertexDeriv.resize(Nv+4,vertexDeriv[0]);  
  
  //Add the five new walls
  
  //New wall connecting v1 to v2
  addWall(*(cell(i).wall(wI)));
  wall(Nw).setIndex(Nw);
  wallData.resize(Nw+1,wallData[0]);
  double tmpLength = 0.0;
  for( size_t d=0 ; d<dimension ; ++d )
    tmpLength += (v1Pos[d]-v2Pos[d])*(v1Pos[d]-v2Pos[d]);
  wallData[Nw][0] = std::sqrt( tmpLength );
  
  //New wall connecting v2 to v3
  addWall(*(cell(i).wall(wI)));
  wall(Nw+1).setIndex(Nw+1);
  wallData.resize(Nw+2,wallData[0]);
  tmpLength = 0.0;
  for( size_t d=0 ; d<dimension ; ++d )
    tmpLength += (v2Pos[d]-v3Pos[d])*(v2Pos[d]-v3Pos[d]);
  wallData[Nw+1][0] = std::sqrt( tmpLength );
  
  //New wall connecting v3 to v4
  addWall(*(cell(i).wall(wI)));
  wall(Nw+2).setIndex(Nw+2);
  wallData.resize(Nw+3,wallData[0]);
  tmpLength = 0.0;
  for( size_t d=0 ; d<dimension ; ++d )
    tmpLength += (v3Pos[d]-v4Pos[d])*(v3Pos[d]-v4Pos[d]);
  wallData[Nw+2][0] = std::sqrt( tmpLength );
  
  //New wall connecting v4 to v1
  addWall(*(cell(i).wall(wI)));
  wall(Nw+3).setIndex(Nw+3);
  wallData.resize(Nw+4,wallData[0]);
  tmpLength = 0.0;
  for( size_t d=0 ; d<dimension ; ++d )
    tmpLength += (v4Pos[d]-v1Pos[d])*(v4Pos[d]-v1Pos[d]);
  wallData[Nw+3][0] = std::sqrt( tmpLength );
  
  
 //Wall continuing the  selected wall
  addWall( *(cell(i).wall(wI)) );
  wall(Nw+4).setIndex(Nw+4);
  //Set new lengths as fractions of the old determined from the new 
  //vertex position
  double oldL = wallData[cell(i).wall(wI)->index()][0];
  size_t v1w = cell(i).wall(wI)->vertex1()->index();
  size_t v2w = cell(i).wall(wI)->vertex2()->index();
  
  tmpLength=0;
  for( size_t d=0 ; d<dimension ; ++d ) {
    tmpLength += (vertexData[v2w][d]-vertexData[v1w][d])*
      (vertexData[v2w][d]-vertexData[v1w][d]);
  }
  tmpLength = std::sqrt( tmpLength ); // current length for the choosen wall for branching
  
  double tmpLengthFrac1=0.0;
  double tmpLengthFrac2=0.0; 
  bool v1posWithv1w;
  // which is closer to which: v1Pos and v2Pos and v1w and v2w
  for( size_t d=0 ; d<dimension ; ++d ) {
    tmpLengthFrac1 += (v1Pos[d]-vertexData[v1w][d])*
      (v1Pos[d]-vertexData[v1w][d]);
    tmpLengthFrac2 += (v2Pos[d]-vertexData[v1w][d])*
      (v2Pos[d]-vertexData[v1w][d]);
  } 
  if (tmpLengthFrac1<tmpLengthFrac2){  
    v1posWithv1w=1;
    tmpLengthFrac1=std::sqrt(tmpLengthFrac1);
    tmpLengthFrac2=0;
    for( size_t d=0 ; d<dimension ; ++d ) {
      tmpLengthFrac2 += (v2Pos[d]-vertexData[v2w][d])*
        (v2Pos[d]-vertexData[v2w][d]);
    }  
    tmpLengthFrac2=std::sqrt(tmpLengthFrac1);
  }
  else {
    v1posWithv1w=0;  
    tmpLengthFrac2=std::sqrt(tmpLengthFrac2);
    tmpLengthFrac1=0;
    for( size_t d=0 ; d<dimension ; ++d ) {
      tmpLengthFrac1 += (v1Pos[d]-vertexData[v2w][d])*
        (v1Pos[d]-vertexData[v2w][d]);
    }  
    tmpLengthFrac1=std::sqrt(tmpLengthFrac1);
  }
         
  double lengthFrac1 = tmpLengthFrac1/tmpLength;
  double lengthFrac2 = tmpLengthFrac2/tmpLength;

  size_t wIG=cell(i).wall(wI)-> index();
  wallData.resize(Nw+5,wallData[wIG]); // adding five new walls
  wallData[wIG][0] = lengthFrac1*oldL; //wall on the same side as v1Pos keeps the old index
  wallData[Nw+4][0] = lengthFrac2*oldL;
  wallData[Nw][0] = oldL*(1-lengthFrac1-lengthFrac2);


  //Resize derivative matrix for walls
  wallDeriv.resize(Nw+5,wallDeriv[0]);



  // Create connection matrix
  // setting up the connectivity

  //extracting the wall and vertex indices for the old cell (cell i)
  std::vector<size_t> oldVIndex,newVIndex,oldWIndex,newWIndex;
  
  for(size_t wIndex=0; wIndex<cell(i).numWall(); wIndex++){
    oldWIndex.push_back( cell(i).wall(wIndex)->index() );
  }
  
  for(size_t vIndex=0; vIndex<cell(i).numVertex(); vIndex++){
    oldVIndex.push_back( cell(i).vertex(vIndex)->index() );
  }

  
  cell(i).addWall( &(wall(0)));
  cell(i).addWall( &(wall(0)));
  cell(i).addVertex( &(vertex(0)));
  cell(i).addVertex( &(vertex(0)));


  vertex(Nv).addCell( &(cell(Nc)));
  vertex(Nv).addCell( &(cell(i)));
  vertex(Nv).addWall( &(wall(Nw)));
  vertex(Nv).addWall( &(wall(Nw+3)));
  vertex(Nv).addWall( &(wall(wIG)));


  vertex(Nv+1).addCell( &(cell(Nc)));
  vertex(Nv+1).addCell( &(cell(i)));
  vertex(Nv+1).addWall( &(wall(Nw)));
  vertex(Nv+1).addWall( &(wall(Nw+1)));
  vertex(Nv+1).addWall( &(wall(Nw+4)));

  vertex(Nv+2).addCell( &(cell(Nc)));
  vertex(Nv+2).addWall( &(wall(Nw+1)));
  vertex(Nv+2).addWall( &(wall(Nw+2)));
  
  vertex(Nv+3).addCell( &(cell(Nc)));
  vertex(Nv+3).addWall( &(wall(Nw+2)));
  vertex(Nv+3).addWall( &(wall(Nw+3)));
    

  
  size_t wIPlus1, wIPlus2;
  wIPlus1= wI+1;
  wIPlus2= wI+2;
    
  if (wIPlus1>cell(i).numWall()-1){
    wIPlus1-= cell(i).numWall();
    wIPlus2-= cell(i).numWall();
  }
  if (wIPlus2 > cell(i).numWall()-1){
    wIPlus2-= cell(i).numWall();
  }
 
  if ( ((cell(i).wall(wI) -> vertex1() -> index())==v1w && v1posWithv1w==1)||
       ((cell(i).wall(wI) -> vertex2() -> index())==v2w && v1posWithv1w==0) ){

    cell(i).setVertex(wIPlus1, &(vertex(Nv+1)));
    cell(i).setVertex(wIPlus2, &(vertex(Nv)));

    cell(i).setWall(wIPlus2, &(wall(wIG)));
    cell(i).setWall(wIPlus1, &(wall(Nw)));
    cell(i).setWall(wI, &(wall(Nw+4)));



    cell(Nc).setWall(0, &(wall(Nw)));
    cell(Nc).setWall(1, &(wall(Nw+3)));
    cell(Nc).setWall(2, &(wall(Nw+2)));
    cell(Nc).setWall(3, &(wall(Nw+1)));

    cell(Nc).setVertex(0, &(vertex(Nv)));
    cell(Nc).setVertex(1, &(vertex(Nv+1)));
    cell(Nc).setVertex(2, &(vertex(Nv+2)));
    cell(Nc).setVertex(3, &(vertex(Nv+3)));
  
    
    
  }
  else{


    cell(i).setVertex(wIPlus1, &(vertex(Nv)));
    cell(i).setVertex(wIPlus2, &(vertex(Nv+1)));
 

    cell(i).setWall(wI, &(wall(wIG)));
    cell(i).setWall(wIPlus1, &(wall(Nw)));
    cell(i).setWall(wIPlus2, &(wall(Nw+4)));


   
    cell(Nc).setWall(0, &(wall(Nw)));
    cell(Nc).setWall(1, &(wall(Nw+1)));
    cell(Nc).setWall(2, &(wall(Nw+2)));
    cell(Nc).setWall(3, &(wall(Nw+3)));


    cell(Nc).setVertex(0, &(vertex(Nv+1)));
    cell(Nc).setVertex(1, &(vertex(Nv)));
    cell(Nc).setVertex(2, &(vertex(Nv+3)));
    cell(Nc).setVertex(3, &(vertex(Nv+2)));

  }
 
  
  for(size_t wIndex=wIPlus2+1; wIndex<cell(i).numWall(); wIndex++){ 
    cell(i).setVertex(wIndex, &(vertex(oldVIndex[wIndex-2])));
    cell(i).setWall(wIndex, &(wall(oldWIndex[wIndex-2])));
  }



  // setting vertices for walls 




  if (v1posWithv1w){
    
    wall(wIG).setVertex1( &(vertex(v1w)) );
    wall(wIG).setVertex2( &(vertex(Nv)) );

    wall(Nw).setVertex1( &(vertex(Nv)) );
    wall(Nw).setVertex2( &(vertex(Nv+1)) );

    wall(Nw+4).setVertex1( &(vertex(Nv+1)) );
    wall(Nw+4).setVertex2( &(vertex(v2w)) );

    wall(Nw+1).setVertex1( &(vertex(Nv+1)) );
    wall(Nw+1).setVertex2( &(vertex(Nv+2)) );

    wall(Nw+2).setVertex1( &(vertex(Nv+2)) );
    wall(Nw+2).setVertex2( &(vertex(Nv+3)) );

    wall(Nw+3).setVertex1( &(vertex(Nv+3)) );
    wall(Nw+3).setVertex2( &(vertex(Nv)) );
  }
  else{

    wall(wIG).setVertex1( &(vertex(Nv)) );
    wall(wIG).setVertex2( &(vertex(v2w)) );

    wall(Nw).setVertex1( &(vertex(Nv+1)) );
    wall(Nw).setVertex2( &(vertex(Nv)) );

    wall(Nw+4).setVertex1( &(vertex(v1w)) );   
    wall(Nw+4).setVertex2( &(vertex(Nv+1)) );

    wall(Nw+1).setVertex1( &(vertex(Nv+2)) );
    wall(Nw+1).setVertex2( &(vertex(Nv+1)) );
    
    wall(Nw+2).setVertex1( &(vertex(Nv+3)) );
    wall(Nw+2).setVertex2( &(vertex(Nv+2)) );
    
    wall(Nw+3).setVertex1( &(vertex(Nv)) );
    wall(Nw+3).setVertex2( &(vertex(Nv+3)) );
   
  }
  
    wall(Nw+1).setCell( &(cell(Nc)), background() );
    wall(Nw+2).setCell( &(cell(Nc)), background() );
    wall(Nw+3).setCell( &(cell(Nc)), background() );

  if (  (wall(wIG).cell1() -> index())==i ){

    wall(Nw+4).setCell( &(cell(i)), background() );
    wall(Nw).setCell( &(cell(i)), &(cell(Nc)) );  
    wall(wIG).setCell( &(cell(i)), background() );
    
    wall(Nw).setCellSort1( wall(wIG).cellSort1() );
    wall(Nw+4).setCellSort1( wall(wIG).cellSort1() );
    wall(Nw).setCellSort2( -(wall(wIG).cellSort1()) );


    wall(Nw+1).setCellSort1( -(wall(wIG).cellSort1()) );
    wall(Nw+2).setCellSort1( -(wall(wIG).cellSort1()) );
    wall(Nw+3).setCellSort1( -(wall(wIG).cellSort1()) );


  }
  else{
   
    wall(Nw+4).setCell(  background(), &(cell(i)) );
    wall(Nw).setCell( &(cell(Nc)), &(cell(i)) );  
    wall(wIG).setCell( background(), &(cell(i)) );


    wall(Nw).setCellSort2( wall(wIG).cellSort2() );
    wall(Nw+4).setCellSort2( wall(wIG).cellSort2() );
    wall(Nw).setCellSort1( -(wall(wIG).cellSort2()) );


    wall(Nw+1).setCellSort1( -(wall(wIG).cellSort2()) );
    wall(Nw+2).setCellSort1( -(wall(wIG).cellSort2()) );
    wall(Nw+3).setCellSort1( -(wall(wIG).cellSort2()) );


  }




//   std::cerr<<"wall 0      "<<wall(0).cell1() ->index()<<"    "<<wall(0).cell2() ->index()<<std::endl;
//   std::cerr<<"wall 1      "<<wall(1).cell1() ->index()<<"    "<<wall(1).cell2() ->index()<<std::endl;
//   std::cerr<<"wall 2      "<<wall(2).cell1() ->index()<<"    "<<wall(2).cell2() ->index()<<std::endl;
//   std::cerr<<"wall 3      "<<wall(3).cell1() ->index()<<"    "<<wall(3).cell2() ->index()<<std::endl;
//   std::cerr<<"wall 4      "<<wall(4).cell1() ->index()<<"    "<<wall(4).cell2() ->index()<<std::endl;
//   std::cerr<<"wall 5      "<<wall(5).cell1() ->index()<<"    "<<wall(5).cell2() ->index()<<std::endl;
//   std::cerr<<"wall 6      "<<wall(6).cell1() ->index()<<"    "<<wall(6).cell2() ->index()<<std::endl;
    
//   std::cerr<<"wall 11     "<<wall(11).cell1() ->index()<<"   "<<wall(11).cell2() ->index()<<std::endl;
//   std::cerr<<"wall 7      "<<wall(7).cell1() ->index()<<"    "<<wall(7).cell2() ->index()<<std::endl;
//   std::cerr<<"wall 8      "<<wall(8).cell1() ->index()<<"    "<<wall(8).cell2() ->index()<<std::endl;
//   std::cerr<<"wall 9      "<<wall(9).cell1() ->index()<<"    "<<wall(9).cell2() ->index()<<std::endl;
//   std::cerr<<"wall 10     "<<wall(10).cell1() ->index()<<"   "<<wall(10).cell2() ->index()<<std::endl;
  


//   std::cerr<<"wall 0  vertices    "<<wall(0).vertex1() ->index()<<"    "<<wall(0).vertex2() ->index()<<std::endl;
//   std::cerr<<"wall 1  vertices    "<<wall(1).vertex1() ->index()<<"    "<<wall(1).vertex2() ->index()<<std::endl;
//   std::cerr<<"wall 2  vertices    "<<wall(2).vertex1() ->index()<<"    "<<wall(2).vertex2() ->index()<<std::endl;
//   std::cerr<<"wall 3  vertices    "<<wall(3).vertex1() ->index()<<"    "<<wall(3).vertex2() ->index()<<std::endl;
//   std::cerr<<"wall 4  vertices    "<<wall(4).vertex1() ->index()<<"    "<<wall(4).vertex2() ->index()<<std::endl;
//   std::cerr<<"wall 5  vertices    "<<wall(5).vertex1() ->index()<<"    "<<wall(5).vertex2() ->index()<<std::endl;
//   std::cerr<<"wall 6  vertices    "<<wall(6).vertex1() ->index()<<"    "<<wall(6).vertex2() ->index()<<std::endl;
//   std::cerr<<"wall 7  vertices    "<<wall(7).vertex1() ->index()<<"    "<<wall(7).vertex2() ->index()<<std::endl;
//   std::cerr<<"wall 8  vertices    "<<wall(8).vertex1() ->index()<<"    "<<wall(8).vertex2() ->index()<<std::endl;
//   std::cerr<<"wall 9  vertices    "<<wall(9).vertex1() ->index()<<"    "<<wall(9).vertex2() ->index()<<std::endl;
//   std::cerr<<"wall 10 vertices    "<<wall(10).vertex1() ->index()<<"   "<<wall(10).vertex2() ->index()<<std::endl;
//   std::cerr<<"wall 11 vertices    "<<wall(11).vertex1() ->index()<<"   "<<wall(11).vertex2() ->index()<<std::endl;

//   std::cerr<<"wall 0  vertices sort   "<<wall(0).cellSort1()<<std::endl;
//   std::cerr<<"wall 1  vertices sort   "<<wall(1).cellSort1()<<std::endl;
//   std::cerr<<"wall 2  vertices sort   "<<wall(2).cellSort1()<<std::endl;
//   std::cerr<<"wall 3  vertices sort   "<<wall(3).cellSort1()<<std::endl;
//   std::cerr<<"wall 4  vertices sort   "<<wall(4).cellSort1()<<std::endl;
//   std::cerr<<"wall 5  vertices sort   "<<wall(5).cellSort1()<<std::endl;
//   std::cerr<<"wall 6  vertices sort   "<<wall(6).cellSort1()<<std::endl;
//   std::cerr<<"wall 7  vertices sort   "<<wall(7).cellSort1()<<std::endl;
//   std::cerr<<"wall 8  vertices sort   "<<wall(8).cellSort1()<<std::endl;
//   std::cerr<<"wall 9  vertices sort   "<<wall(9).cellSort1()<<std::endl;
//   std::cerr<<"wall 10 vertices sort   "<<wall(10).cellSort1()<<std::endl;
//   std::cerr<<"wall 11 vertices sort   "<<wall(11).cellSort1()<<std::endl;
 
// std::cerr<<"size of cell 2  variables...................."<<cellData[2].size();


	
  //checkConnectivity(1);
}







void Tissue::
divideCellCenterTriangulation( Cell *divCell, size_t wI, size_t w3I, 
                               size_t centerIndex,size_t wallLengthIndex, 
                               std::vector<double> &v1Pos,
                               std::vector<double> &v2Pos,
                               DataMatrix &cellData,
                               DataMatrix &wallData,
                               DataMatrix &vertexData,
                               DataMatrix &cellDeriv,
                               DataMatrix &wallDeriv,
                               DataMatrix &vertexDeriv,
                               std::vector<size_t> &volumeChangeList,
                               double threshold,
                               size_t doubleFlag) 

{	
  size_t Nc=numCell(),Nw=numWall(),Nv=numVertex();
  size_t cellIndex = divCell->index();
  size_t dimension = vertexData[0].size();
  size_t numCellWalls=cell(cellIndex).numWall();
  size_t lengthInternalIndex=centerIndex+dimension;
  size_t wIindex=cell(cellIndex).wall(wI) -> index();
  size_t w3Iindex=cell(cellIndex).wall(w3I) -> index();



  // storing internal walls of the neghbor cells if they exist
  size_t cellWI, cellW3I;
  bool WIhasNeighbor,W3IhasNeighbor;
  WIhasNeighbor=false;
  W3IhasNeighbor=false;
 
  if( cell(cellIndex).wall(wI)->cell1()->index()==cellIndex &&
      cell(cellIndex).wall(wI)->cell2() != background() ){
    cellWI=cell(cellIndex).wall(wI)->cell2()->index();
    WIhasNeighbor=true;
  }
  else if( cell(cellIndex).wall(wI)->cell2()->index()==cellIndex &&
           cell(cellIndex).wall(wI)->cell1() != background() ){
    cellWI=cell(cellIndex).wall(wI)->cell1()->index();
    WIhasNeighbor=true;
  }
  else if( cell(cellIndex).wall(wI)->cell1()->index() != cellIndex &&
           cell(cellIndex).wall(wI)->cell2()->index() != cellIndex ){
    std::cerr << "Tissue::divideCell() "
              << "Wall wI not connected to dividing cell" 
              << std::endl;
    exit(-1);
  }
   
  if(cell(cellIndex).wall(w3I)->cell1()->index()==cellIndex &&
     cell(cellIndex).wall(w3I)->cell2() != background() ){
    cellW3I=cell(cellIndex).wall(w3I)->cell2()->index();
    W3IhasNeighbor=true;   
  }
  else if( cell(cellIndex).wall(w3I)->cell2()->index()==cellIndex &&
           cell(cellIndex).wall(w3I)->cell1() != background() ){
    cellW3I=cell(cellIndex).wall(w3I)->cell1()->index();           
    W3IhasNeighbor=true;   
  }
  else if( cell(cellIndex).wall(w3I)->cell1()->index() != cellIndex &&
           cell(cellIndex).wall(w3I)->cell2()->index() != cellIndex ){
    std::cerr << "Tissue::divideCell() "
              << "Wall w3I not connected to dividing cell" 
              << std::endl;
    exit(-1);
  }

 
  // storing growing component of the main walls of the dividing cell and its neighbors 
  // and internal length of the neighbors.
  std::vector<std::vector<double> > cellL, cellWIL, cellWIintL, cellW3IL, cellW3IintL;
  if (doubleFlag==1){ // if double resting length
    cellL.resize(cell(cellIndex).numWall());
    for(size_t kk=0; kk<cell(cellIndex).numWall(); kk++){
      cellL[kk].resize(2);
      cellL[kk][0]=cell(cellIndex).wall(kk) -> index();
      cellL[kk][1]=cellData[cellIndex][centerIndex+dimension+2*cell(cellIndex).numWall()+kk];
    }  
    
    if(WIhasNeighbor){
      size_t numWal=cell(cellWI).numVertex();
      cellWIL.resize(numWal);      
      cellWIintL.resize(2*numWal);
      for(size_t kk=0; kk<numWal; ++kk){
        
        cellWIintL[2*kk].resize(3); // vertex and wall indices + internal wall length
        cellWIintL[2*kk+1].resize(3);
        cellWIL[kk].resize(2); // wall index and wall length

        cellWIintL[2*kk][1]=cell(cellWI).wall((kk+numWal-1)%numWal) -> index();
        cellWIintL[2*kk+1][1]=cell(cellWI).wall(kk) -> index();
       
        cellWIintL[2*kk][0]=cell(cellWI).vertex(kk) -> index();
        cellWIintL[2*kk+1][0]=cell(cellWI).vertex(kk) -> index();

        cellWIintL[2*kk][2]=cellData[cellWI][centerIndex+dimension+2*kk];
        cellWIintL[2*kk+1][2]=cellData[cellWI][centerIndex+dimension+2*kk+1];

        cellWIL[kk][0]=cell(cellWI).wall(kk) -> index();
        cellWIL[kk][1]=cellData[cellWI][centerIndex+dimension+2*numWal+kk];

      }
    }

    if(W3IhasNeighbor){
      size_t numWal=cell(cellW3I).numVertex();
      cellW3IL.resize(numWal);      
      cellW3IintL.resize(2*numWal);
      
      for(size_t kk=0; kk<numWal; ++kk){
        
        cellW3IintL[2*kk].resize(3); // vertex and wall indices + internal wall length
        cellW3IintL[2*kk+1].resize(3);
        cellW3IL[kk].resize(2); // wall index and wall length
        
        cellW3IintL[2*kk][1]=cell(cellW3I).wall((kk+numWal-1)%numWal) -> index();
        cellW3IintL[2*kk+1][1]=cell(cellW3I).wall(kk) -> index();
        
        cellW3IintL[2*kk][0]=cell(cellW3I).vertex(kk) -> index();
        cellW3IintL[2*kk+1][0]=cell(cellW3I).vertex(kk) -> index();
        cellW3IintL[2*kk][2]=cellData[cellW3I][centerIndex+dimension+2*kk];
        cellW3IintL[2*kk+1][2]=cellData[cellW3I][centerIndex+dimension+2*kk+1];

        cellW3IL[kk][0]=cell(cellW3I).wall(kk) -> index();
        cellW3IL[kk][1]=cellData[cellW3I][centerIndex+dimension+2*numWal+kk];

      }
    }
    // for(size_t z=0; z<cellW3IL.size(); z++) // go through the cellW3IL vector
    //   std::cerr<<cellData[cellW3I][centerIndex+dimension+2*cell(cellW3I).numVertex()+z] <<"  ";      

  }
  else if(doubleFlag==0){// if single resting length with averaging
    if(WIhasNeighbor){
      cellWIintL.resize(cell(cellWI).numVertex()); 
      for(size_t kk=0; kk<cell(cellWI).numVertex(); kk++){
        cellWIintL[kk].resize(2);
        cellWIintL[kk][0]=cell(cellWI).vertex(kk) -> index();
        cellWIintL[kk][1]=cellData[cellWI][centerIndex+dimension+kk];
      }
    }
    if(W3IhasNeighbor){
      cellW3IintL.resize(cell(cellW3I).numVertex()); 
      for(size_t kk=0; kk<cell(cellW3I).numVertex(); kk++){
        cellW3IintL[kk].resize(2);
        cellW3IintL[kk][0]=cell(cellW3I).vertex(kk) -> index();
        cellW3IintL[kk][1]=cellData[cellW3I][centerIndex+dimension+kk];
      }
    }
  }

  
  // calculating the averaged strain tensor of the cell
  ////////////begin strain
  double StrainCellGlobal[3][3]={{0,0,0},{0,0,0},{0,0,0}};  
  double TotalCellArea=0;
  for (size_t wallindex=0; wallindex<numCellWalls; ++wallindex) { 
    size_t kPlusOneMod = (wallindex+1)%numCellWalls;
    //size_t v1 = com;
    size_t v2 = cell(cellIndex).vertex(wallindex)->index();
    size_t v3 = cell(cellIndex).vertex(kPlusOneMod)->index();
    //size_t w1 = internal wallindex
    size_t w2 = cell(cellIndex).wall(wallindex)->index();
    //size_t w3 = internal wallindex+1
    
    // Position matrix holds in rows positions for com, vertex(wallindex), vertex(wallindex+1)
    DataMatrix position(3,vertexData[v2]);
    for (size_t d=0; d<dimension; ++d)
      position[0][d] = cellData[cellIndex][centerIndex+d]; // com position
    //position[1] = vertexData[v2]; // given by initiation
    position[2] = vertexData[v3];
    //position[0][2] z for vertex 1 of the current element
    
    
    
    // Resting lengths are from com-vertex(wallindex), 
    // vertex(wallindex)-vertex(wallindex+1) (wall(wallindex)), com-vertex(wallindex+1)
    std::vector<double> restingLength(3);
    if(doubleFlag==1){ // double resting length
      restingLength[0] = cellData[cellIndex][lengthInternalIndex + 2*wallindex+1];
      //restingLength[1] = wallData[w2][wallLengthIndex]+wallData[w2][wallLengthIndex+restAddIndex];
      restingLength[1] = wallData[w2][wallLengthIndex]+
        cellData[cellIndex][lengthInternalIndex+2*numCellWalls+wallindex];
      restingLength[2] = cellData[cellIndex][lengthInternalIndex + 2*kPlusOneMod];      
    }
    else if(doubleFlag==0){    // single resting length
      restingLength[0] = cellData[cellIndex][lengthInternalIndex + wallindex];
      restingLength[1] = wallData[w2][wallLengthIndex];
      restingLength[2] = cellData[cellIndex][lengthInternalIndex + kPlusOneMod];
    }

    // Lengths are from com-vertex(wallindex), 
    // vertex(wallindex)-vertex(wallindex+1) (wall(wallindex)), com-vertex(wallindex+1)
    std::vector<double> length(3);
    length[0] = std::sqrt( (position[0][0]-position[1][0])*(position[0][0]-position[1][0]) +
                           (position[0][1]-position[1][1])*(position[0][1]-position[1][1]) +
                           (position[0][2]-position[1][2])*(position[0][2]-position[1][2]) );
    
    length[1] = std::sqrt( (position[2][0]-position[1][0])*(position[2][0]-position[1][0]) +
                           (position[2][1]-position[1][1])*(position[2][1]-position[1][1]) +
                           (position[2][2]-position[1][2])*(position[2][2]-position[1][2]) );
    
    length[2] = std::sqrt( (position[0][0]-position[2][0])*(position[0][0]-position[2][0]) +
                           (position[0][1]-position[2][1])*(position[0][1]-position[2][1]) +
                           (position[0][2]-position[2][2])*(position[0][2]-position[2][2]) );
    
    
    // restingArea of the element (using Heron's formula)                                      
    double restingArea=std::sqrt( ( restingLength[0]+restingLength[1]+restingLength[2])*
                                  (-restingLength[0]+restingLength[1]+restingLength[2])*
                                  ( restingLength[0]-restingLength[1]+restingLength[2])*
                                  ( restingLength[0]+restingLength[1]-restingLength[2])  )*0.25;
    TotalCellArea +=restingArea;
    
    //Current shape local coordinate of the element  (counterclockwise ordering of nodes/edges)
    double CurrentAngle1=std::acos(  (length[0]*length[0]+length[1]*length[1]-length[2]*length[2])/
                                     (length[0]*length[1]*2)    );
    
    double Qa=std::cos(CurrentAngle1)*length[0];
    double Qc=std::sin(CurrentAngle1)*length[0];
    double Qb=length[1];
    
    double positionLocal[3][2]={ {Qa , Qc}, 
                                 {0  , 0 },  
                                 {Qb , 0 }  };
    
    // Local coordinates of the resting shape ( counterclockwise )
    double RestingAngle1=std::acos(  (restingLength[0]*restingLength[0]
                                      +restingLength[1]*restingLength[1]
                                      -restingLength[2]*restingLength[2])/
                                     (restingLength[0]*restingLength[1]*2)    );
    
    double Pa=std::cos(RestingAngle1)*restingLength[0];
    double Pc=std::sin(RestingAngle1)*restingLength[0];
    double Pb=restingLength[1];
    
    // shape vector matrix in resting shape in local coordinate system  
    //= inverse of coordinate matrix ( only first two elements i.e. ShapeVectorResting[3][2] )      
    double ShapeVectorResting[3][3]={ {  0   ,       1/Pc      , 0 }, 
                                      {-1/Pb , (Pa-Pb)/(Pb*Pc) , 1 },       
                                      { 1/Pb ,     -Pa/(Pb*Pc) , 0 }  };
    
    // Rotation Matrix for changing coordinate systems for both Local to Global( Strain Tensor) 
    // and Global to Local( Aniso Vector in the current shape)
    double rotation[3][3];  
    
    double tempA=std::sqrt((position[2][0]-position[1][0])*(position[2][0]-position[1][0])+
                           (position[2][1]-position[1][1])*(position[2][1]-position[1][1])+
                           (position[2][2]-position[1][2])*(position[2][2]-position[1][2])  );
    
    double tempB=std::sqrt((position[0][0]-position[1][0])*(position[0][0]-position[1][0])+
                         (position[0][1]-position[1][1])*(position[0][1]-position[1][1])+
                           (position[0][2]-position[1][2])*(position[0][2]-position[1][2])  );
    
    
    double Xcurrent[3];      
    Xcurrent[0]= (position[2][0]-position[1][0])/tempA;
    Xcurrent[1]= (position[2][1]-position[1][1])/tempA;
    Xcurrent[2]= (position[2][2]-position[1][2])/tempA;
    
    double Bcurrent[3];      
    Bcurrent[0]= (position[0][0]-position[1][0])/tempB;
    Bcurrent[1]= (position[0][1]-position[1][1])/tempB;
    Bcurrent[2]= (position[0][2]-position[1][2])/tempB;
    
    double Zcurrent[3];      
    Zcurrent[0]= Xcurrent[1]*Bcurrent[2]-Xcurrent[2]*Bcurrent[1];
    Zcurrent[1]= Xcurrent[2]*Bcurrent[0]-Xcurrent[0]*Bcurrent[2];
    Zcurrent[2]= Xcurrent[0]*Bcurrent[1]-Xcurrent[1]*Bcurrent[0];
    
    tempA=std:: sqrt(Zcurrent[0]*Zcurrent[0]+Zcurrent[1]*Zcurrent[1]+Zcurrent[2]*Zcurrent[2]);
    Zcurrent[0]=Zcurrent[0]/tempA;
    Zcurrent[1]=Zcurrent[1]/tempA;
    Zcurrent[2]=Zcurrent[2]/tempA;
    
    double Ycurrent[3];      
    Ycurrent[0]= Zcurrent[1]*Xcurrent[2]-Zcurrent[2]*Xcurrent[1];
    Ycurrent[1]= Zcurrent[2]*Xcurrent[0]-Zcurrent[0]*Xcurrent[2];
    Ycurrent[2]= Zcurrent[0]*Xcurrent[1]-Zcurrent[1]*Xcurrent[0];
    
    
    rotation[0][0]=Xcurrent[0];
    rotation[1][0]=Xcurrent[1];
    rotation[2][0]=Xcurrent[2];
    
    rotation[0][1]=Ycurrent[0];
    rotation[1][1]=Ycurrent[1];
    rotation[2][1]=Ycurrent[2];
    
    rotation[0][2]=Zcurrent[0];
    rotation[1][2]=Zcurrent[1];
    rotation[2][2]=Zcurrent[2];  
    
    
    // STRAIN:
    
    double DeformGrad[2][2]={{0,0},{0,0}}; // F=Sigma i Qi x Di ----------------------
    for ( int ii=0 ; ii<3 ; ++ii ) {
      DeformGrad[0][0]=DeformGrad[0][0]+positionLocal[ii][0]*ShapeVectorResting[ii][0];
      DeformGrad[1][0]=DeformGrad[1][0]+positionLocal[ii][1]*ShapeVectorResting[ii][0];
      DeformGrad[0][1]=DeformGrad[0][1]+positionLocal[ii][0]*ShapeVectorResting[ii][1];
      DeformGrad[1][1]=DeformGrad[1][1]+positionLocal[ii][1]*ShapeVectorResting[ii][1];
    } 
    
    
    double LeftCauchy[2][2]; // B=FFt
    LeftCauchy[0][0]=DeformGrad[0][0]*DeformGrad[0][0]+DeformGrad[0][1]*DeformGrad[0][1];
    LeftCauchy[1][0]=DeformGrad[1][0]*DeformGrad[0][0]+DeformGrad[1][1]*DeformGrad[0][1];
    LeftCauchy[0][1]=DeformGrad[0][0]*DeformGrad[1][0]+DeformGrad[0][1]*DeformGrad[1][1];
    LeftCauchy[1][1]=DeformGrad[1][0]*DeformGrad[1][0]+DeformGrad[1][1]*DeformGrad[1][1];
    
    double StrainAlmansi[2][2]; // e=0.5(1-B^-1)  True strain tensor
    double temp=LeftCauchy[0][0]*LeftCauchy[1][1]-LeftCauchy[1][0]*LeftCauchy[0][1]; // det(B)
    StrainAlmansi[0][0]=0.5*(1-(LeftCauchy[1][1]/temp));
    StrainAlmansi[1][0]=0.5*LeftCauchy[1][0]/temp;
    StrainAlmansi[0][1]=0.5*LeftCauchy[0][1]/temp;  
    StrainAlmansi[1][1]=0.5*(1-(LeftCauchy[0][0]/temp));
  
    
    double StrainTensor[3][3]; // there are other alternatives than StrainAlmansi for strain tensor
    StrainTensor[0][0]=StrainAlmansi[0][0];
    StrainTensor[1][0]=StrainAlmansi[1][0];
    StrainTensor[0][1]=StrainAlmansi[0][1];
    StrainTensor[1][1]=StrainAlmansi[1][1];
    
    
    StrainTensor[0][2]=0;  // adding 3rd dimension which is zero, the tensor is still in element plane
    StrainTensor[1][2]=0;
    StrainTensor[2][2]=0;
    StrainTensor[2][0]=0;
    StrainTensor[2][1]=0;
    
    
    // rotating strain tensor to the global coordinate system
    double tempR[3][3]={{0,0,0},{0,0,0},{0,0,0}};
    for (int r=0 ; r<3 ; r++) 
      for (int s=0 ; s<3 ; s++) 
        for(int w=0 ; w<3 ; w++) 
          tempR[r][s] += rotation[r][w]*StrainTensor[w][s];
    
    for (int r=0 ; r<3 ; r++) 
    for (int s=0 ; s<3 ; s++) {
      StrainTensor[r][s]=0;
      for(int w=0 ; w<3 ; w++) 
        StrainTensor[r][s] += tempR[r][w]*rotation[s][w]; 
    }
    
    for (int r=0 ; r<3 ; r++) 
      for (int s=0 ; s<3 ; s++)    
        StrainCellGlobal[r][s] += restingArea*StrainTensor[r][s];
  }
  
  for (int r=0 ; r<3 ; r++) 
    for (int s=0 ; s<3 ; s++)    
      StrainCellGlobal[r][s]= StrainCellGlobal[r][s]/TotalCellArea; 
  


  
  
  double strainEpcilon =0.000001;
  
  double eigenVectorStrain[3][3]={{1,0,0},{0,1,0},{0,0,1}};
  double pivot=1;
  double pi=3.1415;
  int I,J;
  double RotAngle,Si,Co;
  
  
  
  pivot=std::abs(StrainCellGlobal[1][0]);
  I=1;
  J=0;
  if (std::abs(StrainCellGlobal[2][0])>pivot) {
    pivot=std::abs(StrainCellGlobal[2][0]);
    I=2;
    J=0;
  }
  if (std::abs(StrainCellGlobal[2][1])>pivot) {
    pivot=std::abs(StrainCellGlobal[2][1]);
    I=2;
    J=1;
  }
  
  
  while (pivot> strainEpcilon) {
    
    if (std::abs(StrainCellGlobal[I][I]-StrainCellGlobal[J][J])< strainEpcilon ) {
      RotAngle=pi/4;
    }            
    else {
      RotAngle=0.5*std::atan((2*StrainCellGlobal[I][J])/(StrainCellGlobal[J][J]-StrainCellGlobal[I][I]));
    }
    Si=std::sin(RotAngle);
    Co=std::cos(RotAngle);
    double tempRot[3][3]={{1,0,0},{0,1,0},{0,0,1}};
    tempRot[I][I]=Co;
    tempRot[J][J]=Co;
    tempRot[I][J]=Si;
    tempRot[J][I]=-Si;
    double tempStrain[3][3]={{0,0,0},{0,0,0},{0,0,0}};
    for (int r=0 ; r<3 ; r++) 
      for (int s=0 ; s<3 ; s++) 
        for(int w=0 ; w<3 ; w++) 
          tempStrain[r][s]=tempStrain[r][s]+StrainCellGlobal[r][w]*tempRot[w][s];
    
    for (int r=0 ; r<3 ; r++) 
      for (int s=0 ; s<3 ; s++) 
        StrainCellGlobal[r][s]=0;
    
    for (int r=0 ; r<3 ; r++) 
      for (int s=0 ; s<3 ; s++) 
        for(int w=0 ; w<3 ; w++) 
          StrainCellGlobal[r][s]=StrainCellGlobal[r][s]+tempRot[w][r]*tempStrain[w][s];
    
    for (int r=0 ; r<3 ; r++) 
      for (int s=0 ; s<3 ; s++) 
        tempStrain[r][s]=eigenVectorStrain[r][s];
    
    for (int r=0 ; r<3 ; r++) 
      for (int s=0 ; s<3 ; s++) 
        eigenVectorStrain[r][s]=0;
    
    for (int r=0 ; r<3 ; r++) 
      for (int s=0 ; s<3 ; s++) 
        for(int w=0 ; w<3 ; w++) 
          eigenVectorStrain[r][s]=eigenVectorStrain[r][s]+tempStrain[r][w]*tempRot[w][s];
    
    
    pivot=std::abs(StrainCellGlobal[1][0]);
    I=1;
    J=0;
    if (std::abs(StrainCellGlobal[2][0])>pivot) {
      pivot=std::abs(StrainCellGlobal[2][0]);
      I=2;
      J=0;
    }
    if (std::abs(StrainCellGlobal[2][1])>pivot) {
      pivot=std::abs(StrainCellGlobal[2][1]);
      I=2;
      J=1;
    }
  }
  
  // normalizing eigenvectors (remove if not necessary)  
  double temp=std::sqrt(eigenVectorStrain[0][0]*eigenVectorStrain[0][0] +
                 eigenVectorStrain[1][0]*eigenVectorStrain[1][0] +
                 eigenVectorStrain[2][0]*eigenVectorStrain[2][0] );
  if(temp>0){
    eigenVectorStrain[0][0]/=temp;
    eigenVectorStrain[1][0]/=temp;
    eigenVectorStrain[2][0]/=temp;
  }
  temp=std::sqrt(eigenVectorStrain[0][1]*eigenVectorStrain[0][1] +
                 eigenVectorStrain[1][1]*eigenVectorStrain[1][1] +
                 eigenVectorStrain[2][1]*eigenVectorStrain[2][1] );
  if(temp>0){
    eigenVectorStrain[0][1]/=temp;
    eigenVectorStrain[1][1]/=temp;
    eigenVectorStrain[2][1]/=temp;
  }
  temp=std::sqrt(eigenVectorStrain[0][2]*eigenVectorStrain[0][2] +
                 eigenVectorStrain[1][2]*eigenVectorStrain[1][2] +
                 eigenVectorStrain[2][2]*eigenVectorStrain[2][2] );
  if(temp>0){
    eigenVectorStrain[0][2]/=temp;
    eigenVectorStrain[1][2]/=temp;
    eigenVectorStrain[2][2]/=temp;
  }
  
  // maximal strain direction
  double maximalStrainValue=StrainCellGlobal[0][0];
  int Istrain=0;
  if (StrainCellGlobal[1][1]>maximalStrainValue) 
    {
      maximalStrainValue=StrainCellGlobal[1][1];
      Istrain=1;
    }
  if (StrainCellGlobal[2][2]>maximalStrainValue) 
    {
      maximalStrainValue=StrainCellGlobal[2][2];
      Istrain=2;
    }
  
  
  // 2nd maximalstrain direction/value
  double maximalStrainValue2,maximalStrainValue3;
  int Istrain2,Istrain3;
  if (Istrain==0) {
    Istrain2=1;
    Istrain3=2;
  }
  if (Istrain==1) {
    Istrain2=0;
    Istrain3=2;
  }
  if (Istrain==2) {
    Istrain2=0;
    Istrain3=1;
  }
  if(StrainCellGlobal[Istrain3][Istrain3] > StrainCellGlobal[Istrain2][Istrain2] ) {
    temp=Istrain2;
    Istrain2=Istrain3;
    Istrain3=temp;
  }
  maximalStrainValue2=StrainCellGlobal[Istrain2][Istrain2];
  maximalStrainValue3=StrainCellGlobal[Istrain3][Istrain3];

  // if(std::abs(normalGlob[0]*eigenVectorStrain[0][Istrain]+   //<<<<<<<<<<<<<<<<<<
  //             normalGlob[1]*eigenVectorStrain[1][Istrain]+
  //             normalGlob[2]*eigenVectorStrain[2][Istrain]) > .1) {
    
  //   std::cerr << "max strain vector is out of cell plane in the cell " <<cellIndex <<std::endl;
  //   Istrain=Istrain2; 
  //   Istrain2=Istrain3;
  //   maximalStrainValue=maximalStrainValue2;
  //   maximalStrainValue2=maximalStrainValue3; 
  //}
  //std::cerr<<std::endl;
  std::cerr<<" in the cell "<<cellIndex<<std::endl
           <<"  strain before "
           <<maximalStrainValue <<"  "
           <<maximalStrainValue2<<"  "
           <<maximalStrainValue3<<"  "<<std::endl;

  //////////// end strain
  
  std::vector<double> centroid(3);
  for( size_t ind=0; ind<cell(cellIndex).numWall(); ++ind )
    for( size_t d=0; d<dimension; ++d ){
      centroid[d]+=vertexData[cell(cellIndex).vertex(ind) -> index()][d];
    }
  
  for( size_t d=0; d<dimension; ++d ){
    centroid[d]/=cell(cellIndex).numWall();
  }
  

  std::vector<double> centroidNorm(3);
  for( size_t d=0; d<dimension; ++d ){
    centroidNorm[d]=cellData[cellIndex][centerIndex+d]-centroid[d];
  }

  //Move new vertices if closer than threshold to old vertex
  if( threshold>=0.0 ) {
    Wall *w1 = divCell->wall(wI), *w2 = divCell->wall(w3I);
    double w1L=0.0,w2L=0.0,t1=0.0,t2=0.0;
    for( size_t dim=0; dim<dimension; ++dim ) {
      w1L += (vertexData[w1->vertex1()->index()][dim]-
	      vertexData[w1->vertex2()->index()][dim]) *
	(vertexData[w1->vertex1()->index()][dim] -
	 vertexData[w1->vertex2()->index()][dim]);
      w2L += (vertexData[w2->vertex1()->index()][dim]-
	      vertexData[w2->vertex2()->index()][dim]) *
	(vertexData[w2->vertex1()->index()][dim] -
	 vertexData[w2->vertex2()->index()][dim]);
      t1 += (v1Pos[dim]-vertexData[w1->vertex2()->index()][dim])*
	(v1Pos[dim]-vertexData[w1->vertex2()->index()][dim]);
      t2 += (v2Pos[dim]-vertexData[w2->vertex2()->index()][dim])*
	(v2Pos[dim]-vertexData[w2->vertex2()->index()][dim]);
    }
    w1L = std::sqrt(w1L);
    w2L = std::sqrt(w2L);
    t1 = std::sqrt(t1)/w1L;
    t2 = std::sqrt(t2)/w2L;
    assert( t1>=0.0 && t1<=1.0 );
    assert( t2>=0.0 && t2<=1.0 );
    if( t1<threshold ) {
      std::cerr << "Tissue::divideCell() Moving vertex 1 from "
		<< t1 << " to " << threshold << std::endl;
      t1=threshold;
      for( size_t dim=0; dim<dimension; ++dim ) {
	v1Pos[dim] = vertexData[w1->vertex2()->index()][dim] +
	  t1*(vertexData[w1->vertex1()->index()][dim]-
	      vertexData[w1->vertex2()->index()][dim]);
      }
    }
    else if( t1>(1.0-threshold) ) {
      std::cerr << "Tissue::divideCell() Moving vertex 1 from "
		<< t1 << " to " << 1.0-threshold << std::endl;
      t1 = threshold;
      for( size_t dim=0; dim<dimension; ++dim ) {
	v1Pos[dim] = vertexData[w1->vertex1()->index()][dim] +
	  t1*(vertexData[w1->vertex2()->index()][dim]-
	      vertexData[w1->vertex1()->index()][dim]);
      }
    }
    if( t2<threshold ) {
      std::cerr << "Tissue::divideCell() Moving vertex 2 from "
		<< t2 << " to " << threshold << std::endl;
      t2=threshold;
      for( size_t dim=0; dim<dimension; ++dim ) {
	v2Pos[dim] = vertexData[w2->vertex2()->index()][dim] +
	  t2*(vertexData[w2->vertex1()->index()][dim]-
	      vertexData[w2->vertex2()->index()][dim]);
      }
    }
    else if( t2>(1.0-threshold) ) {
      std::cerr << "Tissue::divideCell() Moving vertex 2 from "
		<< t2 << " to " << 1.0-threshold << std::endl;
      t2 = threshold;
      for( size_t dim=0; dim<dimension; ++dim ) {
	v2Pos[dim] = vertexData[w2->vertex1()->index()][dim] +
	  t2*(vertexData[w2->vertex2()->index()][dim]-
	      vertexData[w2->vertex1()->index()][dim]);
      }
    }		
  }
  
  //Create the new data structure and set indices in the tissue vectors
  //
  //Add the new cell
  addCell( cell(cellIndex) );
  cell(Nc).setIndex(Nc);
  cellData.resize(Nc+1,cellData[cellIndex]);
  cellDeriv.resize(Nc+1,cellDeriv[0]);
  
  //Add the two new vertices
  Vertex tmpVertex;
  tmpVertex.setPosition(v1Pos);
  tmpVertex.setIndex(Nv);
  addVertex(tmpVertex);
  vertexData.resize(Nv+1,v1Pos);  
  tmpVertex.setPosition(v2Pos);
  tmpVertex.setIndex(Nv+1);
  addVertex(tmpVertex);
  vertexData.resize(Nv+2,v2Pos);  
  vertexDeriv.resize(Nv+2,vertexDeriv[0]);  
  
  //Add the three new walls
  Wall tmpWall;
  //New wall dividing old cell into two
  tmpWall.setIndex(Nw);
  addWall(tmpWall);
  wallData.resize(Nw+1,wallData[0]);
  // the new wall calculated based on strain
  std::vector<double> comWallVector(3);
  for( size_t d=0 ; d<dimension ; ++d )
    comWallVector[d]= (v1Pos[d]-v2Pos[d]);
  wallData[Nw][wallLengthIndex] =std::sqrt( 
                                           ( comWallVector[0]*eigenVectorStrain[0][Istrain]
                                             +comWallVector[1]*eigenVectorStrain[1][Istrain]
                                             +comWallVector[2]*eigenVectorStrain[2][Istrain])
                                           *( comWallVector[0]*eigenVectorStrain[0][Istrain]
                                              +comWallVector[1]*eigenVectorStrain[1][Istrain]
                                              +comWallVector[2]*eigenVectorStrain[2][Istrain])
                                           *(1-2*maximalStrainValue)
			    
                                           +( comWallVector[0]*eigenVectorStrain[0][Istrain2]
                                              +comWallVector[1]*eigenVectorStrain[1][Istrain2]
                                              +comWallVector[2]*eigenVectorStrain[2][Istrain2])
                                           *( comWallVector[0]*eigenVectorStrain[0][Istrain2]
                                              +comWallVector[1]*eigenVectorStrain[1][Istrain2]
                                              +comWallVector[2]*eigenVectorStrain[2][Istrain2])
                                           *(1-2*maximalStrainValue2)
                                           
                                           + ( comWallVector[0]*eigenVectorStrain[0][Istrain3]
                                               +comWallVector[1]*eigenVectorStrain[1][Istrain3]
                                               +comWallVector[2]*eigenVectorStrain[2][Istrain3])
                                           *( comWallVector[0]*eigenVectorStrain[0][Istrain3]
                                              +comWallVector[1]*eigenVectorStrain[1][Istrain3]
                                              +comWallVector[2]*eigenVectorStrain[2][Istrain3])
                                           *(1-2*maximalStrainValue3)
                                            );
  if (doubleFlag==1)// for the wall springs
    wallData[Nw][wallLengthIndex+1]=std::sqrt(comWallVector[0]*comWallVector[0]+
                                              comWallVector[1]*comWallVector[1]+
                                              comWallVector[2]*comWallVector[2]);

  
  //Wall continuing the first selected wall
  addWall( *(cell(cellIndex).wall(wI)) );
  wall(Nw+1).setIndex(Nw+1);
  //Set new lengths as fractions of the old determined from the new 
  //vertex position
  double oldL = wallData[cell(cellIndex).wall(wI)->index()][0];
  size_t v1w = cell(cellIndex).wall(wI)->vertex1()->index();
  size_t v2w = cell(cellIndex).wall(wI)->vertex2()->index();
  double tmpLength=0.0; 
  double tmpLengthFrac=0.0; 
  for( size_t d=0 ; d<dimension ; ++d ) {
    tmpLengthFrac += (v1Pos[d]-vertexData[v1w][d])*
      (v1Pos[d]-vertexData[v1w][d]);
    tmpLength += (vertexData[v2w][d]-vertexData[v1w][d])*
      (vertexData[v2w][d]-vertexData[v1w][d]);
  }
  tmpLength = std::sqrt( tmpLength );
  tmpLengthFrac = std::sqrt( tmpLengthFrac );
  double lengthFrac = tmpLengthFrac/tmpLength;
  wallData.resize(Nw+2,wallData[cell(cellIndex).wall(wI)->index()]);
  wallData[cell(cellIndex).wall(wI)->index()][wallLengthIndex] = lengthFrac*oldL;
  wallData[Nw+1][wallLengthIndex] = oldL-wallData[cell(cellIndex).wall(wI)->index()][wallLengthIndex];
  if (doubleFlag==1){// for the wall springs
    wallData[cell(cellIndex).wall(wI)->index()][wallLengthIndex+1] = lengthFrac*oldL;
    wallData[Nw+1][wallLengthIndex+1] = oldL-wallData[cell(cellIndex).wall(wI)->index()][wallLengthIndex];
  }

  //Wall continuing the second selected wall
  addWall(*(cell(cellIndex).wall(w3I)));
  wall(Nw+2).setIndex(Nw+2);
  //Set new lengths as fractions of the old determined from the new
  //vertex position
  oldL = wallData[cell(cellIndex).wall(w3I)->index()][0];
  v1w = cell(cellIndex).wall(w3I)->vertex1()->index();
  v2w = cell(cellIndex).wall(w3I)->vertex2()->index();
  tmpLength=0.0; 
  tmpLengthFrac=0.0; 
  for( size_t d=0 ; d<dimension ; ++d ) {
    tmpLengthFrac += (v2Pos[d]-vertexData[v1w][d])*
      (v2Pos[d]-vertexData[v1w][d]);
    tmpLength += (vertexData[v2w][d]-vertexData[v1w][d])*
      (vertexData[v2w][d]-vertexData[v1w][d]);
  }
  tmpLength = std::sqrt( tmpLength );
  tmpLengthFrac = std::sqrt( tmpLengthFrac );
  lengthFrac = tmpLengthFrac/tmpLength;
  wallData.resize(Nw+3,wallData[cell(cellIndex).wall(w3I)->index()]);
  wallData[cell(cellIndex).wall(w3I)->index()][wallLengthIndex] = lengthFrac*oldL;
  wallData[Nw+2][wallLengthIndex] = oldL-wallData[cell(cellIndex).wall(w3I)->index()][wallLengthIndex];
  if (doubleFlag==1){// for the wall springs
    wallData[cell(cellIndex).wall(w3I)->index()][wallLengthIndex+1] = lengthFrac*oldL;
    wallData[Nw+2][wallLengthIndex+1] = oldL-wallData[cell(cellIndex).wall(w3I)->index()][wallLengthIndex];
  }

  
  //Resize derivative matrix as well
  wallDeriv.resize(Nw+3,wallDeriv[0]);
  
  // Create connection matrix by first selecting vertices and walls for
  // the new and old cells
  //
  std::vector<size_t> oldVIndex,newVIndex,oldWIndex,newWIndex,
    usedWIndex( cell(cellIndex).numWall() );
  
  //Extract walls and vertices for 'old' cell
  size_t tmpWIndex = cell(cellIndex).wall(wI)->index(); 
  size_t tmpVIndex = cell(cellIndex).wall(wI)->vertex1()->index();
  size_t nextW = wI;
  do {
    //add to new vectors
    oldWIndex.push_back( tmpWIndex );
    oldVIndex.push_back( tmpVIndex );
    //mark as used
    usedWIndex[nextW] = 1;
    //usedVIndex[xxx] = 1;
    //find next wall
    nextW = cell(cellIndex).numWall();
    size_t flag=0;
    for( size_t w=0 ; w<cell(cellIndex).numWall() ; ++w ) {
      if( !usedWIndex[w] && (cell(cellIndex).wall(w)->vertex1()->index()==tmpVIndex ||
			     cell(cellIndex).wall(w)->vertex2()->index()==tmpVIndex ) ) {
	nextW = w;
	flag++;
      }
    }
    if( flag != 1 ) {
      std::cerr << "Tissue::divideCellCentertriangulation() " << flag  
		<< " walls marked for next wall..." << std::endl;
      std::cerr << tmpVIndex << std::endl;
      for( size_t w=0 ; w<cell(cellIndex).numWall() ; ++w ) {
	std::cerr << w << " " << usedWIndex[w] << " " 
		  << cell(cellIndex).wall(w)->vertex1()->index() << " "
		  << cell(cellIndex).wall(w)->vertex2()->index() << std::endl;
      }
    }
    assert( flag==1 );
    tmpWIndex = cell(cellIndex).wall(nextW)->index();
    if( cell(cellIndex).wall(nextW)->vertex1()->index()==tmpVIndex )
      tmpVIndex = cell(cellIndex).wall(nextW)->vertex2()->index();
    else if( cell(cellIndex).wall(nextW)->vertex2()->index()==tmpVIndex )
      tmpVIndex = cell(cellIndex).wall(nextW)->vertex1()->index();
    else {
      std::cerr << "Tissue::DivideCellCentertriangulation() " 
		<< "Wrong vertex indices for chosen wall " 
		<< cell(cellIndex).wall(nextW)->vertex1()->index() << " "
		<< cell(cellIndex).wall(nextW)->vertex2()->index() << std::endl;
      exit(-1);
    }
    //std::cerr << count++ << " " << nextW << "  " << tmpWIndex << " " 
    //				<< tmpVIndex << std::endl;
  } while( nextW != w3I ); 
  if( cell(cellIndex).wall(nextW)->vertex1()->index()==oldVIndex[oldVIndex.size()-1] )
    oldWIndex.push_back(cell(cellIndex).wall(w3I)->index());
  else if( cell(cellIndex).wall(nextW)->vertex2()->index()==oldVIndex[oldVIndex.size()-1] )
    oldWIndex.push_back(Nw+2);
  else {
    std::cerr << "Wrong last index for old cell (not in w3I)" << std::endl;
    exit(-1);
  }
  //Extract walls and vertices for 'new' cell
  tmpWIndex = Nw+1;//new copy of wI; 
  tmpVIndex = wall(Nw+1).vertex2()->index();
  //std::cerr << "Old cell done..." << std::endl;
  usedWIndex[wI]=0;
  usedWIndex[w3I]=0;
  nextW=wI;
  do {
    //add to new vectors
    newWIndex.push_back( tmpWIndex );
    newVIndex.push_back( tmpVIndex );
    //mark as used
    usedWIndex[nextW] = 1;
    //usedVIndex[xxx] = 1;
    //find next wall
    nextW = cell(cellIndex).numWall();
    size_t flag=0;
    for( size_t w=0 ; w<cell(cellIndex).numWall() ; ++w ) {
      if( !usedWIndex[w] && ( cell(cellIndex).wall(w)->vertex1()->index()==tmpVIndex ||
			      cell(cellIndex).wall(w)->vertex2()->index()==tmpVIndex ) ) {
	nextW = w;
	flag++;
      }
    }
    if( flag != 1 )
      std::cerr << flag << " walls marked for next wall..." << std::endl;
    assert( flag==1 );
    tmpWIndex = cell(cellIndex).wall(nextW)->index();
    if( cell(cellIndex).wall(nextW)->vertex1()->index()==tmpVIndex )
      tmpVIndex = cell(cellIndex).wall(nextW)->vertex2()->index();
    else if( cell(cellIndex).wall(nextW)->vertex2()->index()==tmpVIndex )
      tmpVIndex = cell(cellIndex).wall(nextW)->vertex1()->index();
    else {
      std::cerr << "Tissue::DivideCell() " 
		<< "Wrong vertex indices for chosen wall " 
		<< cell(cellIndex).wall(nextW)->vertex1() << " "
		<< cell(cellIndex).wall(nextW)->vertex2() << std::endl;
      exit(-1);
    }
  } while( nextW != w3I ); 
  if( cell(cellIndex).wall(nextW)->vertex1()->index()==newVIndex[newVIndex.size()-1] )
    newWIndex.push_back(cell(cellIndex).wall(w3I)->index());
  else if( cell(cellIndex).wall(nextW)->vertex2()->index()==newVIndex[newVIndex.size()-1] )
    newWIndex.push_back(Nw+2);
  else {
    std::cerr << "Wrong last index for new cell (not in w3I)" << std::endl;
    exit(-1);
  }
  
  //std::cerr << "New cell done..." << std::endl;
  
  //Add new vertices
  oldVIndex.push_back(Nv);
  oldVIndex.push_back(Nv+1);
  newVIndex.push_back(Nv);
  newVIndex.push_back(Nv+1);
  
  //Add new walls
  oldWIndex.push_back(Nw);
  newWIndex.push_back(Nw);
  //newWIndex.push_back(Nw+1);
  //newWIndex.push_back(Nw+2);
  
  
  
  //Set vertices and cells for the walls
  wall(Nw).setVertex( &(vertex(Nv)),&(vertex(Nv+1)) );
  wall(Nw).setCell( &(cell(cellIndex)),&(cell(Nc)) );
  assert( wall(Nw).cell1()->index()==cellIndex );
  assert( wall(Nw).cell2()->index()==Nc );
  
  wall(cell(cellIndex).wall(wI)->index()).setVertex2( &(vertex(Nv)) );
  int vInCellFlag=0;
  for( size_t k=0 ; k<oldVIndex.size() ; ++k )
    if( oldVIndex[k]==cell(cellIndex).wall(wI)->vertex1()->index() ) 
      vInCellFlag++;
  assert( vInCellFlag==1 );
  
  wall(Nw+1).setVertex1( &(vertex(Nv)) );
  vInCellFlag=0;
  for( size_t k=0 ; k<newVIndex.size() ; ++k )
    if( newVIndex[k]==wall(Nw+1).vertex2()->index() ) 
      vInCellFlag++;
  assert( vInCellFlag==1 );
  
  if( wall(Nw+1).cell1()->index() == cellIndex ) {
    wall(Nw+1).setCell1( &(cell(Nc)) );
    if( wall(Nw+1).cell2()->index() < Nc ) {
      wall(Nw+1).cell2()->addWall( &(wall(Nw+1)) );
      wall(Nw+1).cell2()->addVertex( &(vertex(Nv)) );
    }
  }
  else if( wall(Nw+1).cell2()->index() == cellIndex ) {
    wall(Nw+1).setCell2( &(cell(Nc)) );
    if( wall(Nw+1).cell1()->index() < Nc ) {
      wall(Nw+1).cell1()->addWall( &(wall(Nw+1)) );
      wall(Nw+1).cell1()->addVertex( &(vertex(Nv)) );
    }
  }
  else {
    std::cerr << "Tissue::divideCellCentertriangulation() "
	      << "First wall not connected to dividing cell" << std::endl;
    std::cerr << cellIndex << " " << cell(cellIndex).index() << "\t" 
	      << wall( cell(cellIndex).wall(wI)->index() ).cell1()->index() << " "
	      << wall(Nw+1).cell1()->index() << "\t"
	      << wall(cell(cellIndex).wall(wI)->index()).cell2()->index() << " "
	      << wall(Nw+1).cell2()->index() << "\n";		
    exit(-1);
  }
  
  vInCellFlag=0;
  size_t wInOldCellFlag=0;
  for( size_t k=0 ; k<oldWIndex.size() ; ++k )
    if( cell(cellIndex).wall(w3I)->index()==oldWIndex[k] )
      wInOldCellFlag=1;
  if( wInOldCellFlag ) {
    for( size_t k=0 ; k<oldVIndex.size() ; ++k )
      if( oldVIndex[k]==wall(cell(cellIndex).wall(w3I)->index()).vertex1()->index() ) 
	vInCellFlag++;
    if( vInCellFlag == 1 ) {
      wall(cell(cellIndex).wall(w3I)->index()).setVertex2( &(vertex(Nv+1)) );
      wall(Nw+2).setVertex1( &(vertex(Nv+1)) );
    }
    else {
      wall(cell(cellIndex).wall(w3I)->index()).setVertex1( &(vertex(Nv+1)) );
      wall(Nw+2).setVertex2( &(vertex(Nv+1)) );
    }
  }
  else {		
    for( size_t k=0 ; k<oldVIndex.size() ; ++k )
      if( oldVIndex[k]==wall(Nw+2).vertex1()->index() ) 
	vInCellFlag++;
    if( vInCellFlag == 1 ) {
      wall(Nw+2).setVertex2( &(vertex(Nv+1)) );
      wall(cell(cellIndex).wall(w3I)->index()).setVertex1( &(vertex(Nv+1)) );
    }
    else {
      wall(Nw+2).setVertex1( &(vertex(Nv+1)) );
      wall(cell(cellIndex).wall(w3I)->index()).setVertex2( &(vertex(Nv+1)) );
    }
  }
  
  //Change cell connection for new wall (w3I or Nw+2)
  size_t newWallIndex = cell(cellIndex).wall(w3I)->index();
  size_t w3IInNewCellFlag = 0;
  for( size_t w=0 ; w<newWIndex.size() ; ++w )
    if( newWIndex[w]==newWallIndex )
      w3IInNewCellFlag++;
  if( !w3IInNewCellFlag++ )
    newWallIndex=Nw+2;
  
  if( wall(newWallIndex).cell1()->index() == cellIndex ) {
    wall(newWallIndex).setCell1( &(cell(Nc)) );
    if( wall(newWallIndex).cell2()->index() < Nc ) {
      wall(newWallIndex).cell2()->addWall( &(wall(Nw+2)) );
      wall(newWallIndex).cell2()->addVertex( &(vertex(Nv+1)) );
    }
  }    
  else if( wall(newWallIndex).cell2()->index() == cellIndex ) {
    wall(newWallIndex).setCell2( &(cell(Nc)) );
    if( wall(newWallIndex).cell1()->index() < Nc ) {
      wall(newWallIndex).cell1()->addWall( &(wall(Nw+2)) );
      wall(newWallIndex).cell1()->addVertex( &(vertex(Nv+1)) );
    }    
  }
  else {
    std::cerr << "Tissue::divideCell() "
	      << "Second wall not connected to dividing cell" << std::endl;
    for( size_t k=0 ; k<cell(cellIndex).numWall() ; ++k ) {
      size_t v1w3Itmp = cell(cellIndex).wall(k)->vertex1()->index();
      size_t v2w3Itmp = cell(cellIndex).wall(k)->vertex2()->index();
      std::cerr << vertexData[v1w3Itmp][0] << " " 
		<< vertexData[v1w3Itmp][1] << " 4\n"
		<< vertexData[v2w3Itmp][0] << " "
		<< vertexData[v2w3Itmp][1] << " 4\n\n\n";
    }
    exit(-1);
  }
  
  //Set cells and walls for the vertices

  for( size_t v=0 ; v<newVIndex.size() ; ++v ) {
    std::vector<Cell*> tmpCell;
    std::vector<Wall*> tmpWall;
    if( newVIndex[v]==Nv ) {
      //First new vertex
      tmpCell.push_back( &(cell(cellIndex)) );
      tmpCell.push_back( &(cell(Nc)) );
      //plus one more
      if( cell(cellIndex).wall(wI)->cell1()->index()==cellIndex &&
	  cell(cellIndex).wall(wI)->cell2() != background() ){
	tmpCell.push_back( cell(cellIndex).wall(wI)->cell2() );
      }
      else if( cell(cellIndex).wall(wI)->cell2()->index()==cellIndex &&
               cell(cellIndex).wall(wI)->cell1() != background() ){
	tmpCell.push_back( cell(cellIndex).wall(wI)->cell1() );
      }
      else if( cell(cellIndex).wall(wI)->cell1()->index() != cellIndex &&
	       cell(cellIndex).wall(wI)->cell2()->index() != cellIndex ){
	std::cerr << "Tissue::divideCell() "
		  << "Wall wI not connected to dividing cell" 
		  << std::endl;
	exit(-1);
      }
      tmpWall.push_back( cell(cellIndex).wall(wI) );
      tmpWall.push_back( &(wall(Nw)) );
      tmpWall.push_back( &(wall(Nw+1)) );
    }
    else if( newVIndex[v]==Nv+1 ) {
      //Second new vertex
      tmpCell.push_back( &(cell(cellIndex)) );
      tmpCell.push_back( &(cell(Nc)) );
      //plus one more
      if( cell(cellIndex).wall(w3I)->cell1()->index()==cellIndex ||
	  cell(cellIndex).wall(w3I)->cell1()->index()==Nc ) {
	if( cell(cellIndex).wall(w3I)->cell2() != background() ){
	  tmpCell.push_back( cell(cellIndex).wall(w3I)->cell2() );
        }
      }
      else if( cell(cellIndex).wall(w3I)->cell2()->index()==cellIndex ||
	       cell(cellIndex).wall(w3I)->cell2()->index()==Nc ) {
	if( cell(cellIndex).wall(w3I)->cell1() != background() ){
	  tmpCell.push_back( cell(cellIndex).wall(w3I)->cell1() );
        }
      }
      else {
	std::cerr << "Tissue::divideCell() "
		  << "Wall w3I not connected to dividing cell" 
		  << std::endl;
	exit(-1);
      }
      tmpWall.push_back( cell(cellIndex).wall(w3I) );
      tmpWall.push_back( &(wall(Nw)) );
      tmpWall.push_back( &(wall(Nw+2)) );
    }
    else {
      //For other vertices
      for( size_t c=0 ; c<vertex(newVIndex[v]).numCell() ; ++c )
	if( vertex(newVIndex[v]).cell(c)->index() == cellIndex )
	  vertex(newVIndex[v]).setCell(c,&(cell(Nc)));
      for( size_t w=0 ; w<vertex(newVIndex[v]).numWall() ; ++w ) {
	if( vertex(newVIndex[v]).wall(w)->index() == 
	    cell(cellIndex).wall(wI)->index() )
	  vertex(newVIndex[v]).setWall(w,&(wall(Nw+1)));
      	else if( vertex(newVIndex[v]).wall(w)->index() == 
		 cell(cellIndex).wall(w3I)->index() ) {
	  size_t w3IInNewFlag=0;
	  for( size_t ww=0 ; ww<newWIndex.size() ; ++ww )
	    if( newWIndex[ww] == cell(cellIndex).wall(w3I)->index() )
	      w3IInNewFlag++;
	  if( !w3IInNewFlag )
	    vertex(newVIndex[v]).setWall(w,&(wall(Nw+2)));
	}
      }
    }
    if( newVIndex[v]>=Nv ) {
      vertex(newVIndex[v]).setCell(tmpCell);
      vertex(newVIndex[v]).setWall(tmpWall);
    }
  }
  
  //Finally check if vertices in old cell connected to w3I should be changed
  //to Nw+2
  for( size_t v=0 ; v<oldVIndex.size() ; ++v ) {
    if( oldVIndex[v] < Nv ) {
      for( size_t w=0 ; w<vertex(oldVIndex[v]).numWall() ; ++w ) {
      	if( vertex(oldVIndex[v]).wall(w)->index() == 
	    cell(cellIndex).wall(w3I)->index() ) {
	  size_t w3IInOldFlag=0;
	  for( size_t ww=0 ; ww<oldWIndex.size() ; ++ww )
	    if( oldWIndex[ww] == cell(cellIndex).wall(w3I)->index() )
	      w3IInOldFlag++;
	  if( !w3IInOldFlag )
	    vertex(oldVIndex[v]).setWall(w,&(wall(Nw+2)));
	}
      }
    }
  }

  //Set walls and vertices for the cells
  //This should be done after the vertex and wall insertion not to
  //mess up wI and w3I
  std::vector<Wall*> tmpW(oldWIndex.size());
  std::vector<Vertex*> tmpV(oldVIndex.size());
  for( size_t w=0 ; w<oldWIndex.size() ; ++w )
    tmpW[w] = &(wall(oldWIndex[w]));
  for( size_t v=0 ; v<oldVIndex.size() ; ++v )
    tmpV[v] = &(vertex(oldVIndex[v]));
  cell(cellIndex).setWall(tmpW);
  cell(cellIndex).setVertex(tmpV);
  
  tmpW.resize( newWIndex.size() );
  tmpV.resize( newVIndex.size() );
  for( size_t w=0 ; w<newWIndex.size() ; ++w ) {
    tmpW[w] = &(wall(newWIndex[w]));
    if( tmpW[w]->index() != Nw && tmpW[w]->cell1()->index()==cellIndex )
      tmpW[w]->setCell1( &(cell(Nc)) );
    else if( tmpW[w]->index() != Nw && tmpW[w]->cell2()->index()==cellIndex )
      tmpW[w]->setCell2( &(cell(Nc)) );
  }
  for( size_t v=0 ; v<newVIndex.size() ; ++v )
    tmpV[v] = &(vertex(newVIndex[v]));
  cell(Nc).setWall(tmpW);
  cell(Nc).setVertex(tmpV);
  assert( wall(Nw).cell1()->index()==cellIndex );
  assert( wall(Nw).cell2()->index()==Nc );
  
  //Update directional wall vector if applicable
  updateDirectionDivision(cellIndex,cellData,wallData,vertexData,
			  cellDeriv,wallDeriv,vertexDeriv);
  
  // Update the volume dependent variables for each cell variable index 
  // given in volumeChangeList
  double fi;
  double fn;
  if (volumeChangeList.size()) {

    cell(cellIndex).sortWallAndVertex(*this);
    cell(Nc).sortWallAndVertex(*this);
    if(WIhasNeighbor)
      cell(cellWI).sortWallAndVertex(*this);
    if(W3IhasNeighbor)
      cell(cellW3I).sortWallAndVertex(*this);

    
    double Vi = cell(cellIndex).calculateVolume(vertexData);
    double Vn = cell(Nc).calculateVolume(vertexData);
    fi = Vi/(Vi+Vn);
    fn = Vn/(Vi+Vn);
    
    for (size_t k=0; k<volumeChangeList.size(); ++k) {
      cellData[cellIndex][volumeChangeList[k]] *= fi;
      cellData[Nc][volumeChangeList[k]] *= fn;
    }
  }		
  //checkConnectivity(1);
  
  

  // sorting all the cells involved in division
  cell(cellIndex).sortWallAndVertex(*this);
  cell(Nc).sortWallAndVertex(*this);
  if(WIhasNeighbor)
    cell(cellWI).sortWallAndVertex(*this);
  if(W3IhasNeighbor)
    cell(cellW3I).sortWallAndVertex(*this);
  
  
  
 

  //setting cell centers and internal walls in cellData for cells cellIndex and Nc
  if (doubleFlag==1){
    cellData[cellIndex].resize(centerIndex+dimension+3*cell(cellIndex).numWall());
    cellData[Nc].resize(centerIndex+dimension+3*cell(Nc).numWall());
    
    cellDeriv[cellIndex].resize(centerIndex+dimension+3*cell(cellIndex).numWall());
    cellDeriv[Nc].resize(centerIndex+dimension+3*cell(Nc).numWall());
  }
  if (doubleFlag==0){
    cellData[cellIndex].resize(centerIndex+dimension+cell(cellIndex).numWall());
    cellData[Nc].resize(centerIndex+dimension+cell(Nc).numWall());
    
    cellDeriv[cellIndex].resize(centerIndex+dimension+cell(cellIndex).numWall());
    cellDeriv[Nc].resize(centerIndex+dimension+cell(Nc).numWall());
  }
  // cell centroids 
  for (size_t d=0; d< dimension; ++d){
    cellData[cellIndex][centerIndex+d]=0;
    cellData[Nc][centerIndex+d]=0;
  }
  
  for (size_t k=0; k<cell(cellIndex).numWall(); ++k)
    for (size_t d=0; d< dimension; ++d)
      cellData[cellIndex][centerIndex+d]+=vertexData[cell(cellIndex).vertex(k)->index()][d];
  
  for (size_t k=0; k<cell(Nc).numWall(); ++k)
    for (size_t d=0; d< dimension; ++d)
      cellData[Nc][centerIndex+d]+=vertexData[cell(Nc).vertex(k)->index()][d];
  
  for (size_t d=0; d< dimension; ++d){
    cellData[cellIndex][centerIndex+d]/=cell(cellIndex).numWall();
    cellData[Nc][centerIndex+d]/=cell(Nc).numWall();
  }
  
  //centroidNorm
  for (size_t d=0; d< dimension; ++d){
    cellData[cellIndex][centerIndex+d]+=fi*centroidNorm[d];
    cellData[Nc][centerIndex+d]+=fn*centroidNorm[d];
  }
  
  std::vector<std::vector<double> > 
    internalVectorI(cell(cellIndex).numWall(), std::vector<double>(3)), 
    internalVectorNc(cell(Nc).numWall(), std::vector<double>(3)) ;
  
  for (size_t k=0; k<cell(cellIndex).numWall(); ++k) 
    for (size_t d=0; d< dimension; ++d)
      internalVectorI[k][d]=
        vertexData[cell(cellIndex).vertex(k)->index()][d]
        -cellData[cellIndex][centerIndex+d];
  for (size_t k=0; k<cell(Nc).numWall(); ++k) 
    for (size_t d=0; d< dimension; ++d)
      internalVectorNc[k][d]=
        vertexData[cell(Nc).vertex(k)->index()][d]
        -cellData[Nc][centerIndex+d];
  

  
  //strain
  if (doubleFlag==1){
    for (size_t k=0; k<cell(cellIndex).numWall(); ++k){ 
      cellData[cellIndex][centerIndex+dimension+2*k]=
        std::sqrt( 
                  ( internalVectorI[k][0]*eigenVectorStrain[0][Istrain]
                    +internalVectorI[k][1]*eigenVectorStrain[1][Istrain]
                    +internalVectorI[k][2]*eigenVectorStrain[2][Istrain])
                  *( internalVectorI[k][0]*eigenVectorStrain[0][Istrain]
                     +internalVectorI[k][1]*eigenVectorStrain[1][Istrain]
                     +internalVectorI[k][2]*eigenVectorStrain[2][Istrain])
                  *(1-2*maximalStrainValue)
                  
                  +( internalVectorI[k][0]*eigenVectorStrain[0][Istrain2]
                     +internalVectorI[k][1]*eigenVectorStrain[1][Istrain2]
                     +internalVectorI[k][2]*eigenVectorStrain[2][Istrain2])
                  *( internalVectorI[k][0]*eigenVectorStrain[0][Istrain2]
                     +internalVectorI[k][1]*eigenVectorStrain[1][Istrain2]
                     +internalVectorI[k][2]*eigenVectorStrain[2][Istrain2])
                  *(1-2*maximalStrainValue2)
                  
                  +( internalVectorI[k][0]*eigenVectorStrain[0][Istrain3]
                     +internalVectorI[k][1]*eigenVectorStrain[1][Istrain3]
                     +internalVectorI[k][2]*eigenVectorStrain[2][Istrain3])
                  *( internalVectorI[k][0]*eigenVectorStrain[0][Istrain3]
                     +internalVectorI[k][1]*eigenVectorStrain[1][Istrain3]
                     +internalVectorI[k][2]*eigenVectorStrain[2][Istrain3])
                  *(1-2*maximalStrainValue3)
                   );
      cellData[cellIndex][centerIndex+dimension+2*k+1]=cellData[cellIndex][centerIndex+dimension+2*k];
    }
    for (size_t k=0; k<cell(Nc).numWall(); ++k){ 
      cellData[Nc][centerIndex+dimension+2*k]=
        std::sqrt( 
                  ( internalVectorNc[k][0]*eigenVectorStrain[0][Istrain]
                    +internalVectorNc[k][1]*eigenVectorStrain[1][Istrain]
                    +internalVectorNc[k][2]*eigenVectorStrain[2][Istrain])
                  *( internalVectorNc[k][0]*eigenVectorStrain[0][Istrain]
                     +internalVectorNc[k][1]*eigenVectorStrain[1][Istrain]
                     +internalVectorNc[k][2]*eigenVectorStrain[2][Istrain])
                  *(1-2*maximalStrainValue)
                  
                  +( internalVectorNc[k][0]*eigenVectorStrain[0][Istrain2]
                     +internalVectorNc[k][1]*eigenVectorStrain[1][Istrain2]
                     +internalVectorNc[k][2]*eigenVectorStrain[2][Istrain2])
                  *( internalVectorNc[k][0]*eigenVectorStrain[0][Istrain2]
                     +internalVectorNc[k][1]*eigenVectorStrain[1][Istrain2]
                     +internalVectorNc[k][2]*eigenVectorStrain[2][Istrain2])
                  *(1-2*maximalStrainValue2)
                  
                  +( internalVectorNc[k][0]*eigenVectorStrain[0][Istrain3]
                     +internalVectorNc[k][1]*eigenVectorStrain[1][Istrain3]
                     +internalVectorNc[k][2]*eigenVectorStrain[2][Istrain3])
                  *( internalVectorNc[k][0]*eigenVectorStrain[0][Istrain3]
                     +internalVectorNc[k][1]*eigenVectorStrain[1][Istrain3]
                     +internalVectorNc[k][2]*eigenVectorStrain[2][Istrain3])
                  *(1-2*maximalStrainValue3)
                   );
      cellData[Nc][centerIndex+dimension+2*k+1]=cellData[Nc][centerIndex+dimension+2*k];
    }
  }
  

  if (doubleFlag==0){
    for (size_t k=0; k<cell(cellIndex).numWall(); ++k) 
      cellData[cellIndex][centerIndex+dimension+k]=
        std::sqrt( 
                  ( internalVectorI[k][0]*eigenVectorStrain[0][Istrain]
                    +internalVectorI[k][1]*eigenVectorStrain[1][Istrain]
                    +internalVectorI[k][2]*eigenVectorStrain[2][Istrain])
                  *( internalVectorI[k][0]*eigenVectorStrain[0][Istrain]
                     +internalVectorI[k][1]*eigenVectorStrain[1][Istrain]
                     +internalVectorI[k][2]*eigenVectorStrain[2][Istrain])
                  *(1-2*maximalStrainValue)
                  
                  +( internalVectorI[k][0]*eigenVectorStrain[0][Istrain2]
                     +internalVectorI[k][1]*eigenVectorStrain[1][Istrain2]
                     +internalVectorI[k][2]*eigenVectorStrain[2][Istrain2])
                  *( internalVectorI[k][0]*eigenVectorStrain[0][Istrain2]
                     +internalVectorI[k][1]*eigenVectorStrain[1][Istrain2]
                     +internalVectorI[k][2]*eigenVectorStrain[2][Istrain2])
                  *(1-2*maximalStrainValue2)
                  
                  +( internalVectorI[k][0]*eigenVectorStrain[0][Istrain3]
                     +internalVectorI[k][1]*eigenVectorStrain[1][Istrain3]
                     +internalVectorI[k][2]*eigenVectorStrain[2][Istrain3])
                  *( internalVectorI[k][0]*eigenVectorStrain[0][Istrain3]
                     +internalVectorI[k][1]*eigenVectorStrain[1][Istrain3]
                     +internalVectorI[k][2]*eigenVectorStrain[2][Istrain3])
                  *(1-2*maximalStrainValue3)
                   );
    
    for (size_t k=0; k<cell(Nc).numWall(); ++k) 
      cellData[Nc][centerIndex+dimension+k]=
        std::sqrt( 
                  ( internalVectorNc[k][0]*eigenVectorStrain[0][Istrain]
                    +internalVectorNc[k][1]*eigenVectorStrain[1][Istrain]
                    +internalVectorNc[k][2]*eigenVectorStrain[2][Istrain])
                  *( internalVectorNc[k][0]*eigenVectorStrain[0][Istrain]
                     +internalVectorNc[k][1]*eigenVectorStrain[1][Istrain]
                     +internalVectorNc[k][2]*eigenVectorStrain[2][Istrain])
                  *(1-2*maximalStrainValue)
                  
                  +( internalVectorNc[k][0]*eigenVectorStrain[0][Istrain2]
                     +internalVectorNc[k][1]*eigenVectorStrain[1][Istrain2]
                     +internalVectorNc[k][2]*eigenVectorStrain[2][Istrain2])
                  *( internalVectorNc[k][0]*eigenVectorStrain[0][Istrain2]
                     +internalVectorNc[k][1]*eigenVectorStrain[1][Istrain2]
                     +internalVectorNc[k][2]*eigenVectorStrain[2][Istrain2])
                  *(1-2*maximalStrainValue2)
                  
                  +( internalVectorNc[k][0]*eigenVectorStrain[0][Istrain3]
                     +internalVectorNc[k][1]*eigenVectorStrain[1][Istrain3]
                     +internalVectorNc[k][2]*eigenVectorStrain[2][Istrain3])
                  *( internalVectorNc[k][0]*eigenVectorStrain[0][Istrain3]
                     +internalVectorNc[k][1]*eigenVectorStrain[1][Istrain3]
                     +internalVectorNc[k][2]*eigenVectorStrain[2][Istrain3])
                  *(1-2*maximalStrainValue3)
                   );
  }

  // setting growing components of main walls in cellIndex and Nc
  if (doubleFlag==1){// double resting length
    size_t numWallCellIndex=cell(cellIndex).numWall();
    size_t numWallNc=cell(Nc).numWall();
    double wIcomp, w3Icomp;
    
    for(size_t z=0; z<numWallCellIndex+numWallNc-4; z++){ // go through the cellL vector
      if(cellL[z][0]==wIindex)
        wIcomp=cellL[z][1];
      if(cellL[z][0]==w3Iindex)
        w3Icomp=cellL[z][1];
    }
    


    for(size_t yy=0; yy<numWallCellIndex; yy++){// for each wall in the new cellIndex
      size_t wallyy=cell(cellIndex).wall(yy) -> index(); // get the index of the wall

      for(size_t z=0; z<numWallCellIndex+numWallNc-4; z++) // go through the cellL vector
        if(cellL[z][0]==wallyy)
          cellData[cellIndex][centerIndex+dimension+2*numWallCellIndex+yy]=cellL[z][1];
      
      if(wallyy==wIindex)
        cellData[cellIndex][centerIndex+dimension+2*numWallCellIndex+yy]=
          (wallData[wIindex][0]/(wallData[wIindex][0]+wallData[Nw+1][0]))*wIcomp;
      if(wallyy==w3Iindex)
        cellData[cellIndex][centerIndex+dimension+2*numWallCellIndex+yy]=
          (wallData[w3Iindex][0]/(wallData[w3Iindex][0]+wallData[Nw+2][0]))*w3Icomp;
      if(wallyy==Nw)//New wall dividing old cell is new with zero growing component
        cellData[cellIndex][centerIndex+dimension+2*numWallCellIndex+yy]=0;
      if(wallyy==Nw+1)
        cellData[cellIndex][centerIndex+dimension+2*numWallCellIndex+yy]=
          (wallData[Nw+1][0]/(wallData[wIindex][0]+wallData[Nw+1][0]))*wIcomp;
      if(wallyy==Nw+2)
        cellData[cellIndex][centerIndex+dimension+2*numWallCellIndex+yy]=
          (wallData[Nw+2][0]/(wallData[w3Iindex][0]+wallData[Nw+2][0]))*w3Icomp;
    }



    for(size_t yy=0; yy<numWallNc; yy++){// for each wall in the new Nc
      size_t wallyy=cell(Nc).wall(yy) -> index(); // get the index of the wall
      
      for(size_t z=0; z<numWallCellIndex+numWallNc-4; z++) // go through the cellL vector
        if(cellL[z][0]==wallyy)
          cellData[Nc][centerIndex+dimension+2*numWallNc+yy]=cellL[z][1];
      
      if(wallyy==wIindex)
        cellData[Nc][centerIndex+dimension+2*numWallNc+yy]=
          (wallData[wIindex][0]/(wallData[wIindex][0]+wallData[Nw+1][0]))*wIcomp;
      if(wallyy==w3Iindex)
        cellData[Nc][centerIndex+dimension+2*numWallNc+yy]=
          (wallData[w3Iindex][0]/(wallData[w3Iindex][0]+wallData[Nw+2][0]))*w3Icomp;
      if(wallyy==Nw)//New wall dividing old cell is new with zero growing component
        cellData[Nc][centerIndex+dimension+2*numWallNc+yy]=0;
      if(wallyy==Nw+1)
        cellData[Nc][centerIndex+dimension+2*numWallNc+yy]=
          (wallData[Nw+1][0]/(wallData[wIindex][0]+wallData[Nw+1][0]))*wIcomp;
      if(wallyy==Nw+2)
        cellData[Nc][centerIndex+dimension+2*numWallNc+yy]=
          (wallData[Nw+2][0]/(wallData[w3Iindex][0]+wallData[Nw+2][0]))*w3Icomp;
      
    }
   
    if(WIhasNeighbor){

      cellData[cellWI].resize(cellData[cellWI].size()+3);
      size_t numWallCellWI=cell(cellWI).numWall();

      for(size_t z=0; z<numWallCellWI-1; z++){ // go through the cellWIL vector
        if(cellWIL[z][0]==wIindex)
          wIcomp=cellWIL[z][1];
      }

      for(size_t yy=0; yy<numWallCellWI; yy++){// for each wall in the cellWI
        size_t wallyy=cell(cellWI).wall(yy) -> index(); // get the index of the wall
        
        for(size_t z=0; z<numWallCellWI-1; z++) // go through the cellWIL vector
          if(cellWIL[z][0]==wallyy)
            cellData[cellWI][centerIndex+dimension+2*numWallCellWI+yy]=cellWIL[z][1];
        
        if(wallyy==wIindex)
          cellData[cellWI][centerIndex+dimension+2*numWallCellWI+yy]=
            (wallData[wIindex][0]/(wallData[wIindex][0]+wallData[Nw+1][0]))*wIcomp;

        if(wallyy==Nw+1)
          cellData[cellWI][centerIndex+dimension+2*numWallCellWI+yy]=
            (wallData[Nw+1][0]/(wallData[wIindex][0]+wallData[Nw+1][0]))*wIcomp;
        
      }
    }
    
    if(W3IhasNeighbor){
      
      cellData[cellW3I].resize(cellData[cellW3I].size()+3);
      size_t numWallCellW3I=cell(cellW3I).numWall();
      
      for(size_t z=0; z<numWallCellW3I-1; z++){ // go through the cellW3IL vector
        if(cellW3IL[z][0]==w3Iindex)
          w3Icomp=cellW3IL[z][1];
      }

      for(size_t yy=0; yy<numWallCellW3I; yy++){// for each wall in the cellW3IL
        size_t wallyy=cell(cellW3I).wall(yy) -> index(); // get the index of the wall
        
        for(size_t z=0; z<numWallCellW3I-1; z++) // go through the cellW3IL vector
          if(cellW3IL[z][0]==wallyy)
            cellData[cellW3I][centerIndex+dimension+2*numWallCellW3I+yy]=cellW3IL[z][1];
        
        if(wallyy==w3Iindex)
          cellData[cellW3I][centerIndex+dimension+2*numWallCellW3I+yy]=
            (wallData[w3Iindex][0]/(wallData[w3Iindex][0]+wallData[Nw+2][0]))*w3Icomp;
                
        if(wallyy==Nw+2)
          cellData[cellW3I][centerIndex+dimension+2*numWallCellW3I+yy]=
            (wallData[Nw+2][0]/(wallData[w3Iindex][0]+wallData[Nw+2][0]))*w3Icomp;

               
      }
    }
  }
  
  //Update directional wall vector if applicable 
  updateDirectionDivision(cellIndex,cellData,wallData,vertexData,
			  cellDeriv,wallDeriv,vertexDeriv);
  
  
  //setting internal length of neighbors if applicable
  if (WIhasNeighbor){
    size_t numVcellWI=cell(cellWI).numVertex();    
    if(doubleFlag==1){
      //cellData[cellWI].resize(cellData[cellWI].size()+3);  // it is done before
      for(size_t kk=0; kk<numVcellWI; kk++){                                    // go through vertices
        size_t kkPlusOne = (kk+1)%numVcellWI;
        size_t VcellWI=cell(cellWI).vertex(kk)-> index();                       // pick the vertex index
        size_t WcellWI=cell(cellWI).wall(kk)-> index();                         // pick the wall index
        size_t VcellWIPlusOne=cell(cellWI).vertex(kkPlusOne)-> index();         // pick the next vertex index
        for(size_t yy=0; yy<2*numVcellWI-2; yy++){                              // for each vertex go through  

          if(VcellWI==cellWIintL[yy][0] && (WcellWI==cellWIintL[yy][1] || WcellWI==Nw+1) )
            cellData[cellWI][centerIndex+dimension+2*kk+1]=cellWIintL[yy][2];
          
          if(VcellWIPlusOne==cellWIintL[yy][0] && (WcellWI==cellWIintL[yy][1] || WcellWI==Nw+1))
            cellData[cellWI][centerIndex+dimension+2*kkPlusOne]=cellWIintL[yy][2];
        }
      }
    }
    if(doubleFlag==0){
      cellData[cellWI].resize(cellData[cellWI].size()+1);
      for(size_t kk=0; kk<numVcellWI; kk++){
        size_t VcellWI=cell(cellWI).vertex(kk)-> index();
        for(size_t yy=0; yy<numVcellWI-1; yy++)
          if(VcellWI==cellWIintL[yy][0])
            cellData[cellWI][centerIndex+dimension+kk]=cellWIintL[yy][1];
      }
    }
    
    size_t WI1, WI2, InsrtLind;
    double newLength;
    bool WI1Chosen=false, WI2Chosen=false;
    size_t vl1,vl2;

    for( size_t wl=0 ; wl<numVcellWI ; wl++ ) {
      if((cell(cellWI).wall(wl)->vertex1()-> index()==Nv ||cell(cellWI).wall(wl)->vertex2()-> index()==Nv) 
         && !WI1Chosen){
        WI1=cell(cellWI).wall(wl)-> index();
        WI1Chosen=true;
        if(cell(cellWI).wall(wl)->vertex1()-> index()==Nv)
          vl1=cell(cellWI).wall(wl)->vertex2()-> index();
        else
          vl1=cell(cellWI).wall(wl)->vertex1()-> index();
      }
      if((cell(cellWI).wall(wl)->vertex1()-> index()==Nv ||cell(cellWI).wall(wl)->vertex2()-> index()==Nv) 
         && WI1Chosen){
        WI2=cell(cellWI).wall(wl)-> index();
        WI2Chosen=true;
        if(cell(cellWI).wall(wl)->vertex1()-> index()==Nv)
          vl2=cell(cellWI).wall(wl)->vertex2()-> index();
        else
          vl2=cell(cellWI).wall(wl)->vertex1()-> index();
      }
    }
    if(WI1Chosen && WI2Chosen){
      size_t vl1temp,vl2temp;      
      for(size_t iii=0; iii<numVcellWI; iii++){
        if(cell(cellWI).vertex(iii)-> index()==vl1)
          vl1temp=iii;
        if(cell(cellWI).vertex(iii)-> index()==vl2)
          vl2temp=iii;
        if(cell(cellWI).vertex(iii)-> index()==Nv){
          InsrtLind=iii;
        }
      }
      vl1=vl1temp;
      vl2=vl2temp;
      double l1,l2,l3,l4;
      if(doubleFlag==1){
        if(InsrtLind==0){
          l1=cellData[cellWI][centerIndex+dimension+2*vl1];
          l2=cellData[cellWI][centerIndex+dimension+2*vl2+1];
          l3=wallData[WI1][wallLengthIndex]+
            cellData[cellWI][centerIndex+dimension+2*numVcellWI+InsrtLind];
          l4=wallData[WI2][wallLengthIndex]+
            cellData[cellWI][centerIndex+dimension+2*numVcellWI+vl2];           }
        else{
          l1=cellData[cellWI][centerIndex+dimension+2*vl1+1];
          l2=cellData[cellWI][centerIndex+dimension+2*vl2];
          l3=wallData[WI1][wallLengthIndex]+
            cellData[cellWI][centerIndex+dimension+2*numVcellWI+vl1];
          l4=wallData[WI2][wallLengthIndex]+
            cellData[cellWI][centerIndex+dimension+2*numVcellWI+InsrtLind];
        }
      }
      if(doubleFlag==0){
        l1=cellData[cellWI][centerIndex+dimension+vl1];
        l2=cellData[cellWI][centerIndex+dimension+vl2];   
        l3=wallData[WI1][wallLengthIndex];
        l4=wallData[WI2][wallLengthIndex];
      }      
      newLength=std::sqrt(((l2*l2*l3+l1*l1*l4)/(l3+l4))-l3*l4);
      if(doubleFlag==1){
        cellData[cellWI][centerIndex+dimension+2*InsrtLind]=newLength;
        cellData[cellWI][centerIndex+dimension+2*InsrtLind+1]=newLength;
      }
      if(doubleFlag==0){
        cellData[cellWI][centerIndex+dimension+InsrtLind]=newLength;
      }
    }
    else{
      std::cerr<<"in celldivision something goes wrong"<<std::endl;
      exit(-1);      
    }
    
    cellDeriv[cellWI].resize(cellData[cellWI].size());
  }

 
  if (W3IhasNeighbor){
    size_t numVcellW3I=cell(cellW3I).numVertex();
    if(doubleFlag==1){
      //cellData[cellW3I].resize(cellData[cellW3I].size()+3);  // it is done before
      for(size_t kk=0; kk<numVcellW3I; kk++){
        size_t kkPlusOne = (kk+1)%numVcellW3I;
        size_t VcellW3I=cell(cellW3I).vertex(kk)-> index();
        size_t WcellW3I=cell(cellW3I).wall(kk)-> index();
        size_t VcellW3IPlusOne=cell(cellW3I).vertex(kkPlusOne)-> index();
        for(size_t yy=0; yy<2*numVcellW3I-2; yy++){

          if(VcellW3I==cellW3IintL[yy][0] && (WcellW3I==cellW3IintL[yy][1] ||WcellW3I==Nw+2) )
            cellData[cellW3I][centerIndex+dimension+2*kk+1]=cellW3IintL[yy][2];

          if(VcellW3IPlusOne==cellW3IintL[yy][0] && (WcellW3I==cellW3IintL[yy][1] ||WcellW3I==Nw+2) )
            cellData[cellW3I][centerIndex+dimension+2*kkPlusOne]=cellW3IintL[yy][2];
        }
      }
    }
    
    if(doubleFlag==0){
      cellData[cellW3I].resize(cellData[cellW3I].size()+1);
      for(size_t kk=0; kk<numVcellW3I; kk++){
        size_t VcellW3I=cell(cellW3I).vertex(kk)-> index();
        for(size_t yy=0; yy<numVcellW3I-1; yy++)
          if(VcellW3I==cellW3IintL[yy][0])
            cellData[cellW3I][centerIndex+dimension+kk]=cellW3IintL[yy][1];
      }
    }
    
    size_t W3I1, W3I2, InsrtLind;
    double newLength;
    bool W3I1Chosen=false, W3I2Chosen=false;
    size_t vl1,vl2;

    for( size_t wl=0 ; wl<numVcellW3I ; wl++ ) {
      if((cell(cellW3I).wall(wl)->vertex1()-> index()==Nv+1 || cell(cellW3I).wall(wl)->vertex2()-> index()==Nv+1) 
         && !W3I1Chosen){
        W3I1=cell(cellW3I).wall(wl)-> index();
        W3I1Chosen=true;
        if(cell(cellW3I).wall(wl)->vertex1()-> index()==Nv+1)
          vl1=cell(cellW3I).wall(wl)->vertex2()-> index();
        else
          vl1=cell(cellW3I).wall(wl)->vertex1()-> index();
      }
      if((cell(cellW3I).wall(wl)->vertex1()-> index()==Nv+1 || cell(cellW3I).wall(wl)->vertex2()-> index()==Nv+1) 
         && W3I1Chosen){
        W3I2=cell(cellW3I).wall(wl)-> index();
        W3I2Chosen=true;
        if(cell(cellW3I).wall(wl)->vertex1()-> index()==Nv+1)
          vl2=cell(cellW3I).wall(wl)->vertex2()-> index();
        else
          vl2=cell(cellW3I).wall(wl)->vertex1()-> index();
      }
    }
    if(W3I1Chosen && W3I2Chosen){
      size_t vl1temp,vl2temp;      
      for(size_t iii=0; iii<numVcellW3I; iii++){
        if(cell(cellW3I).vertex(iii)-> index()==vl1)
          vl1temp=iii;
        if(cell(cellW3I).vertex(iii)-> index()==vl2)
          vl2temp=iii;
        if(cell(cellW3I).vertex(iii)-> index()==Nv+1)
          InsrtLind=iii;
      }
      vl1=vl1temp;
      vl2=vl2temp;
      double l1,l2,l3,l4;
      if(doubleFlag==1){
        if(InsrtLind==0){
          l1=cellData[cellW3I][centerIndex+dimension+2*vl1];
          l2=cellData[cellW3I][centerIndex+dimension+2*vl2+1];
          l3=wallData[W3I1][wallLengthIndex]+
            cellData[cellW3I][centerIndex+dimension+2*numVcellW3I+InsrtLind];
          l4=wallData[W3I2][wallLengthIndex]+
            cellData[cellW3I][centerIndex+dimension+2*numVcellW3I+vl2];           }
        else{
          l1=cellData[cellW3I][centerIndex+dimension+2*vl1+1];
          l2=cellData[cellW3I][centerIndex+dimension+2*vl2];
          l3=wallData[W3I1][wallLengthIndex]+
            cellData[cellW3I][centerIndex+dimension+2*numVcellW3I+vl1];
          l4=wallData[W3I2][wallLengthIndex]+
            cellData[cellW3I][centerIndex+dimension+2*numVcellW3I+InsrtLind];
        }
      }
      if(doubleFlag==0){
        l1=cellData[cellW3I][centerIndex+dimension+vl1];
        l2=cellData[cellW3I][centerIndex+dimension+vl2];   
        l3=wallData[W3I1][wallLengthIndex];
        l4=wallData[W3I2][wallLengthIndex];
      }      
      newLength=std::sqrt(((l2*l2*l3+l1*l1*l4)/(l3+l4))-l3*l4);
      if(doubleFlag==1){
        cellData[cellW3I][centerIndex+dimension+2*InsrtLind]=newLength;
        cellData[cellW3I][centerIndex+dimension+2*InsrtLind+1]=newLength;
      }
      if(doubleFlag==0){
        cellData[cellW3I][centerIndex+dimension+InsrtLind]=newLength;
      }
    }
    else{
      std::cerr<<"in celldivision something goes wrong"<<std::endl;
      exit(-1);      
    }
    
    cellDeriv[cellW3I].resize(cellData[cellW3I].size());
  }

  
  
 ///////////////////////////////////////// adhoc for comparing strain before and after division
 numCellWalls=cell(cellIndex).numWall();
  // calculating the averaged strain tensor of the cell
  ////////////begin strain
 for(size_t gg=0; gg<3;gg++)  
   for(size_t ggk=0; ggk<3;ggk++)  
     StrainCellGlobal[gg][ggk]=0;
 
  TotalCellArea=0;
  for (size_t wallindex=0; wallindex<numCellWalls; ++wallindex) { 
    size_t kPlusOneMod = (wallindex+1)%numCellWalls;
    //size_t v1 = com;
    size_t v2 = cell(cellIndex).vertex(wallindex)->index();
    size_t v3 = cell(cellIndex).vertex(kPlusOneMod)->index();
    //size_t w1 = internal wallindex
    size_t w2 = cell(cellIndex).wall(wallindex)->index();
    //size_t w3 = internal wallindex+1
    
    // Position matrix holds in rows positions for com, vertex(wallindex), vertex(wallindex+1)
    DataMatrix position(3,vertexData[v2]);
    for (size_t d=0; d<dimension; ++d)
      position[0][d] = cellData[cellIndex][centerIndex+d]; // com position
    //position[1] = vertexData[v2]; // given by initiation
    position[2] = vertexData[v3];
    //position[0][2] z for vertex 1 of the current element
    
    
      
    // Resting lengths are from com-vertex(wallindex), 
    // vertex(wallindex)-vertex(wallindex+1) (wall(wallindex)), com-vertex(wallindex+1)
    std::vector<double> restingLength(3);
    if(doubleFlag==1){ // double resting length
      restingLength[0] = cellData[cellIndex][lengthInternalIndex + 2*wallindex+1];
      //restingLength[1] = wallData[w2][wallLengthIndex]+wallData[w2][wallLengthIndex+restAddIndex];
      restingLength[1] = wallData[w2][wallLengthIndex]+
        cellData[cellIndex][lengthInternalIndex+2*numCellWalls+wallindex];
      restingLength[2] = cellData[cellIndex][lengthInternalIndex + 2*kPlusOneMod];      
    }
    else if(doubleFlag==0){    // single resting length
      restingLength[0] = cellData[cellIndex][lengthInternalIndex + wallindex];
      restingLength[1] = wallData[w2][wallLengthIndex];
      restingLength[2] = cellData[cellIndex][lengthInternalIndex + kPlusOneMod];
    }
    // Lengths are from com-vertex(wallindex), 
    // vertex(wallindex)-vertex(wallindex+1) (wall(wallindex)), com-vertex(wallindex+1)
    std::vector<double> length(3);
    length[0] = std::sqrt( (position[0][0]-position[1][0])*(position[0][0]-position[1][0]) +
                           (position[0][1]-position[1][1])*(position[0][1]-position[1][1]) +
                           (position[0][2]-position[1][2])*(position[0][2]-position[1][2]) );
    
    length[1] = std::sqrt( (position[2][0]-position[1][0])*(position[2][0]-position[1][0]) +
                           (position[2][1]-position[1][1])*(position[2][1]-position[1][1]) +
                           (position[2][2]-position[1][2])*(position[2][2]-position[1][2]) );
    
    length[2] = std::sqrt( (position[0][0]-position[2][0])*(position[0][0]-position[2][0]) +
                           (position[0][1]-position[2][1])*(position[0][1]-position[2][1]) +
                           (position[0][2]-position[2][2])*(position[0][2]-position[2][2]) );
    
    
    // restingArea of the element (using Heron's formula)                                      
    double restingArea=std::sqrt( ( restingLength[0]+restingLength[1]+restingLength[2])*
                                  (-restingLength[0]+restingLength[1]+restingLength[2])*
                                  ( restingLength[0]-restingLength[1]+restingLength[2])*
                                  ( restingLength[0]+restingLength[1]-restingLength[2])  )*0.25;
    TotalCellArea +=restingArea;
    
    //Current shape local coordinate of the element  (counterclockwise ordering of nodes/edges)
    double CurrentAngle1=std::acos(  (length[0]*length[0]+length[1]*length[1]-length[2]*length[2])/
                                     (length[0]*length[1]*2)    );
    
    double Qa=std::cos(CurrentAngle1)*length[0];
    double Qc=std::sin(CurrentAngle1)*length[0];
    double Qb=length[1];
    
    double positionLocal[3][2]={ {Qa , Qc}, 
                                 {0  , 0 },  
                                 {Qb , 0 }  };
    
    // Local coordinates of the resting shape ( counterclockwise )
    double RestingAngle1=std::acos(  (restingLength[0]*restingLength[0]
                                      +restingLength[1]*restingLength[1]
                                      -restingLength[2]*restingLength[2])/
                                     (restingLength[0]*restingLength[1]*2)    );
    
    double Pa=std::cos(RestingAngle1)*restingLength[0];
    double Pc=std::sin(RestingAngle1)*restingLength[0];
    double Pb=restingLength[1];
    
    // shape vector matrix in resting shape in local coordinate system  
    //= inverse of coordinate matrix ( only first two elements i.e. ShapeVectorResting[3][2] )      
    double ShapeVectorResting[3][3]={ {  0   ,       1/Pc      , 0 }, 
                                      {-1/Pb , (Pa-Pb)/(Pb*Pc) , 1 },       
                                      { 1/Pb ,     -Pa/(Pb*Pc) , 0 }  };
    
    // Rotation Matrix for changing coordinate systems for both Local to Global( Strain Tensor) 
    // and Global to Local( Aniso Vector in the current shape)
    double rotation[3][3];  
    
    double tempA=std::sqrt((position[2][0]-position[1][0])*(position[2][0]-position[1][0])+
                           (position[2][1]-position[1][1])*(position[2][1]-position[1][1])+
                           (position[2][2]-position[1][2])*(position[2][2]-position[1][2])  );
    
    double tempB=std::sqrt((position[0][0]-position[1][0])*(position[0][0]-position[1][0])+
                         (position[0][1]-position[1][1])*(position[0][1]-position[1][1])+
                           (position[0][2]-position[1][2])*(position[0][2]-position[1][2])  );
    
    
    double Xcurrent[3];      
    Xcurrent[0]= (position[2][0]-position[1][0])/tempA;
    Xcurrent[1]= (position[2][1]-position[1][1])/tempA;
    Xcurrent[2]= (position[2][2]-position[1][2])/tempA;
    
    double Bcurrent[3];      
    Bcurrent[0]= (position[0][0]-position[1][0])/tempB;
    Bcurrent[1]= (position[0][1]-position[1][1])/tempB;
    Bcurrent[2]= (position[0][2]-position[1][2])/tempB;
    
    double Zcurrent[3];      
    Zcurrent[0]= Xcurrent[1]*Bcurrent[2]-Xcurrent[2]*Bcurrent[1];
    Zcurrent[1]= Xcurrent[2]*Bcurrent[0]-Xcurrent[0]*Bcurrent[2];
    Zcurrent[2]= Xcurrent[0]*Bcurrent[1]-Xcurrent[1]*Bcurrent[0];
    
    tempA=std:: sqrt(Zcurrent[0]*Zcurrent[0]+Zcurrent[1]*Zcurrent[1]+Zcurrent[2]*Zcurrent[2]);
    Zcurrent[0]=Zcurrent[0]/tempA;
    Zcurrent[1]=Zcurrent[1]/tempA;
    Zcurrent[2]=Zcurrent[2]/tempA;
    
    double Ycurrent[3];      
    Ycurrent[0]= Zcurrent[1]*Xcurrent[2]-Zcurrent[2]*Xcurrent[1];
    Ycurrent[1]= Zcurrent[2]*Xcurrent[0]-Zcurrent[0]*Xcurrent[2];
    Ycurrent[2]= Zcurrent[0]*Xcurrent[1]-Zcurrent[1]*Xcurrent[0];
    
    
    rotation[0][0]=Xcurrent[0];
    rotation[1][0]=Xcurrent[1];
    rotation[2][0]=Xcurrent[2];
    
    rotation[0][1]=Ycurrent[0];
    rotation[1][1]=Ycurrent[1];
    rotation[2][1]=Ycurrent[2];
    
    rotation[0][2]=Zcurrent[0];
    rotation[1][2]=Zcurrent[1];
    rotation[2][2]=Zcurrent[2];  
    
    
    // STRAIN:
    
    double DeformGrad[2][2]={{0,0},{0,0}}; // F=Sigma i Qi x Di ----------------------
    for ( int ii=0 ; ii<3 ; ++ii ) {
      DeformGrad[0][0]=DeformGrad[0][0]+positionLocal[ii][0]*ShapeVectorResting[ii][0];
      DeformGrad[1][0]=DeformGrad[1][0]+positionLocal[ii][1]*ShapeVectorResting[ii][0];
      DeformGrad[0][1]=DeformGrad[0][1]+positionLocal[ii][0]*ShapeVectorResting[ii][1];
      DeformGrad[1][1]=DeformGrad[1][1]+positionLocal[ii][1]*ShapeVectorResting[ii][1];
    } 
    
    
    double LeftCauchy[2][2]; // B=FFt
    LeftCauchy[0][0]=DeformGrad[0][0]*DeformGrad[0][0]+DeformGrad[0][1]*DeformGrad[0][1];
    LeftCauchy[1][0]=DeformGrad[1][0]*DeformGrad[0][0]+DeformGrad[1][1]*DeformGrad[0][1];
    LeftCauchy[0][1]=DeformGrad[0][0]*DeformGrad[1][0]+DeformGrad[0][1]*DeformGrad[1][1];
    LeftCauchy[1][1]=DeformGrad[1][0]*DeformGrad[1][0]+DeformGrad[1][1]*DeformGrad[1][1];
    
    double StrainAlmansi[2][2]; // e=0.5(1-B^-1)  True strain tensor
    double temp=LeftCauchy[0][0]*LeftCauchy[1][1]-LeftCauchy[1][0]*LeftCauchy[0][1]; // det(B)
    StrainAlmansi[0][0]=0.5*(1-(LeftCauchy[1][1]/temp));
    StrainAlmansi[1][0]=0.5*LeftCauchy[1][0]/temp;
    StrainAlmansi[0][1]=0.5*LeftCauchy[0][1]/temp;  
    StrainAlmansi[1][1]=0.5*(1-(LeftCauchy[0][0]/temp));
  
    
    double StrainTensor[3][3]; // there are other alternatives than StrainAlmansi for strain tensor
    StrainTensor[0][0]=StrainAlmansi[0][0];
    StrainTensor[1][0]=StrainAlmansi[1][0];
    StrainTensor[0][1]=StrainAlmansi[0][1];
    StrainTensor[1][1]=StrainAlmansi[1][1];
    
    
    StrainTensor[0][2]=0;  // adding 3rd dimension which is zero, the tensor is still in element plane
    StrainTensor[1][2]=0;
    StrainTensor[2][2]=0;
    StrainTensor[2][0]=0;
    StrainTensor[2][1]=0;
    
    
    // rotating strain tensor to the global coordinate system
    double tempR[3][3]={{0,0,0},{0,0,0},{0,0,0}};
    for (int r=0 ; r<3 ; r++) 
      for (int s=0 ; s<3 ; s++) 
        for(int w=0 ; w<3 ; w++) 
          tempR[r][s] += rotation[r][w]*StrainTensor[w][s];
    
    for (int r=0 ; r<3 ; r++) 
    for (int s=0 ; s<3 ; s++) {
      StrainTensor[r][s]=0;
      for(int w=0 ; w<3 ; w++) 
        StrainTensor[r][s] += tempR[r][w]*rotation[s][w]; 
    }
    
    for (int r=0 ; r<3 ; r++) 
      for (int s=0 ; s<3 ; s++)    
        StrainCellGlobal[r][s] += restingArea*StrainTensor[r][s];
  }
  
  for (int r=0 ; r<3 ; r++) 
    for (int s=0 ; s<3 ; s++)    
      StrainCellGlobal[r][s]= StrainCellGlobal[r][s]/TotalCellArea; 
  


  
  
  strainEpcilon =0.000001;
  for(size_t gg=0; gg<3;gg++)  
    for(size_t ggk=0; ggk<3;ggk++)
      if(gg==ggk)  
        eigenVectorStrain[gg][ggk]=1;
      else
        eigenVectorStrain[gg][ggk]=0;
  
  pivot=1;
  pi=3.1415;
  
   //double RotAngle,Si,Co;
  
  
  
  pivot=std::abs(StrainCellGlobal[1][0]);
  I=1;
  J=0;
  if (std::abs(StrainCellGlobal[2][0])>pivot) {
    pivot=std::abs(StrainCellGlobal[2][0]);
    I=2;
    J=0;
  }
  if (std::abs(StrainCellGlobal[2][1])>pivot) {
    pivot=std::abs(StrainCellGlobal[2][1]);
    I=2;
    J=1;
  }
  
  
  while (pivot> strainEpcilon) {
    
    if (std::abs(StrainCellGlobal[I][I]-StrainCellGlobal[J][J])< strainEpcilon ) {
      RotAngle=pi/4;
    }            
    else {
      RotAngle=0.5*std::atan((2*StrainCellGlobal[I][J])/(StrainCellGlobal[J][J]-StrainCellGlobal[I][I]));
    }
    Si=std::sin(RotAngle);
    Co=std::cos(RotAngle);
    double tempRot[3][3]={{1,0,0},{0,1,0},{0,0,1}};
    tempRot[I][I]=Co;
    tempRot[J][J]=Co;
    tempRot[I][J]=Si;
    tempRot[J][I]=-Si;
    double tempStrain[3][3]={{0,0,0},{0,0,0},{0,0,0}};
    for (int r=0 ; r<3 ; r++) 
      for (int s=0 ; s<3 ; s++) 
        for(int w=0 ; w<3 ; w++) 
          tempStrain[r][s]=tempStrain[r][s]+StrainCellGlobal[r][w]*tempRot[w][s];
    
    for (int r=0 ; r<3 ; r++) 
      for (int s=0 ; s<3 ; s++) 
        StrainCellGlobal[r][s]=0;
    
    for (int r=0 ; r<3 ; r++) 
      for (int s=0 ; s<3 ; s++) 
        for(int w=0 ; w<3 ; w++) 
          StrainCellGlobal[r][s]=StrainCellGlobal[r][s]+tempRot[w][r]*tempStrain[w][s];
    
    for (int r=0 ; r<3 ; r++) 
      for (int s=0 ; s<3 ; s++) 
        tempStrain[r][s]=eigenVectorStrain[r][s];
    
    for (int r=0 ; r<3 ; r++) 
      for (int s=0 ; s<3 ; s++) 
        eigenVectorStrain[r][s]=0;
    
    for (int r=0 ; r<3 ; r++) 
      for (int s=0 ; s<3 ; s++) 
        for(int w=0 ; w<3 ; w++) 
          eigenVectorStrain[r][s]=eigenVectorStrain[r][s]+tempStrain[r][w]*tempRot[w][s];
    
    
    pivot=std::abs(StrainCellGlobal[1][0]);
    I=1;
    J=0;
    if (std::abs(StrainCellGlobal[2][0])>pivot) {
      pivot=std::abs(StrainCellGlobal[2][0]);
      I=2;
      J=0;
    }
    if (std::abs(StrainCellGlobal[2][1])>pivot) {
      pivot=std::abs(StrainCellGlobal[2][1]);
      I=2;
      J=1;
    }
  }
  
  // normalizing eigenvectors (remove if not necessary)  
  temp=std::sqrt(eigenVectorStrain[0][0]*eigenVectorStrain[0][0] +
                 eigenVectorStrain[1][0]*eigenVectorStrain[1][0] +
                 eigenVectorStrain[2][0]*eigenVectorStrain[2][0] );
  if(temp>0){
    eigenVectorStrain[0][0]/=temp;
    eigenVectorStrain[1][0]/=temp;
    eigenVectorStrain[2][0]/=temp;
  }
  temp=std::sqrt(eigenVectorStrain[0][1]*eigenVectorStrain[0][1] +
                 eigenVectorStrain[1][1]*eigenVectorStrain[1][1] +
                 eigenVectorStrain[2][1]*eigenVectorStrain[2][1] );
  if(temp>0){
    eigenVectorStrain[0][1]/=temp;
    eigenVectorStrain[1][1]/=temp;
    eigenVectorStrain[2][1]/=temp;
  }
  temp=std::sqrt(eigenVectorStrain[0][2]*eigenVectorStrain[0][2] +
                 eigenVectorStrain[1][2]*eigenVectorStrain[1][2] +
                 eigenVectorStrain[2][2]*eigenVectorStrain[2][2] );
  if(temp>0){
    eigenVectorStrain[0][2]/=temp;
    eigenVectorStrain[1][2]/=temp;
    eigenVectorStrain[2][2]/=temp;
  }
  
  // maximal strain direction
  maximalStrainValue=StrainCellGlobal[0][0];
  Istrain=0;
  if (StrainCellGlobal[1][1]>maximalStrainValue) 
    {
      maximalStrainValue=StrainCellGlobal[1][1];
      Istrain=1;
    }
  if (StrainCellGlobal[2][2]>maximalStrainValue) 
    {
      maximalStrainValue=StrainCellGlobal[2][2];
      Istrain=2;
    }
  
  
  // 2nd maximalstrain direction/value
  //double maximalStrainValue2,maximalStrainValue3;
  //int Istrain2,Istrain3;
  if (Istrain==0) {
    Istrain2=1;
    Istrain3=2;
  }
  if (Istrain==1) {
    Istrain2=0;
    Istrain3=2;
  }
  if (Istrain==2) {
    Istrain2=0;
    Istrain3=1;
  }
  if(StrainCellGlobal[Istrain3][Istrain3] > StrainCellGlobal[Istrain2][Istrain2] ) {
    temp=Istrain2;
    Istrain2=Istrain3;
    Istrain3=temp;
  }
  maximalStrainValue2=StrainCellGlobal[Istrain2][Istrain2];
  maximalStrainValue3=StrainCellGlobal[Istrain3][Istrain3];

  std::cerr<<"  strain after 1 "
           <<maximalStrainValue <<"  "
           <<maximalStrainValue2<<"  "
           <<maximalStrainValue3<<"  "<<std::endl;


  cellIndex=Nc;
 // adhoc for comparing strain before and after division
 numCellWalls=cell(cellIndex).numWall();
  // calculating the averaged strain tensor of the cell
  ////////////begin strain
 for(size_t gg=0; gg<3;gg++)  
   for(size_t ggk=0; ggk<3;ggk++)  
     StrainCellGlobal[gg][ggk]=0;
 
  TotalCellArea=0;
  for (size_t wallindex=0; wallindex<numCellWalls; ++wallindex) { 
    size_t kPlusOneMod = (wallindex+1)%numCellWalls;
    //size_t v1 = com;
    size_t v2 = cell(cellIndex).vertex(wallindex)->index();
    size_t v3 = cell(cellIndex).vertex(kPlusOneMod)->index();
    //size_t w1 = internal wallindex
    size_t w2 = cell(cellIndex).wall(wallindex)->index();
    //size_t w3 = internal wallindex+1
    
    // Position matrix holds in rows positions for com, vertex(wallindex), vertex(wallindex+1)
    DataMatrix position(3,vertexData[v2]);
    for (size_t d=0; d<dimension; ++d)
      position[0][d] = cellData[cellIndex][centerIndex+d]; // com position
    //position[1] = vertexData[v2]; // given by initiation
    position[2] = vertexData[v3];
    //position[0][2] z for vertex 1 of the current element
    
    
    // Resting lengths are from com-vertex(wallindex), 
    // vertex(wallindex)-vertex(wallindex+1) (wall(wallindex)), com-vertex(wallindex+1)
    std::vector<double> restingLength(3);
    if(doubleFlag==1){ // double resting length
      restingLength[0] = cellData[cellIndex][lengthInternalIndex + 2*wallindex+1];
      //restingLength[1] = wallData[w2][wallLengthIndex]+wallData[w2][wallLengthIndex+restAddIndex];
      restingLength[1] = wallData[w2][wallLengthIndex]+
        cellData[cellIndex][lengthInternalIndex+2*numCellWalls+wallindex];
      restingLength[2] = cellData[cellIndex][lengthInternalIndex + 2*kPlusOneMod];      
    }
    else if(doubleFlag==0){    // single resting length
      restingLength[0] = cellData[cellIndex][lengthInternalIndex + wallindex];
      restingLength[1] = wallData[w2][wallLengthIndex];
      restingLength[2] = cellData[cellIndex][lengthInternalIndex + kPlusOneMod];
    }
    // Lengths are from com-vertex(wallindex), 
    // vertex(wallindex)-vertex(wallindex+1) (wall(wallindex)), com-vertex(wallindex+1)
    std::vector<double> length(3);
    length[0] = std::sqrt( (position[0][0]-position[1][0])*(position[0][0]-position[1][0]) +
                           (position[0][1]-position[1][1])*(position[0][1]-position[1][1]) +
                           (position[0][2]-position[1][2])*(position[0][2]-position[1][2]) );
    
    length[1] = std::sqrt( (position[2][0]-position[1][0])*(position[2][0]-position[1][0]) +
                           (position[2][1]-position[1][1])*(position[2][1]-position[1][1]) +
                           (position[2][2]-position[1][2])*(position[2][2]-position[1][2]) );
    
    length[2] = std::sqrt( (position[0][0]-position[2][0])*(position[0][0]-position[2][0]) +
                           (position[0][1]-position[2][1])*(position[0][1]-position[2][1]) +
                           (position[0][2]-position[2][2])*(position[0][2]-position[2][2]) );
    
    
    // restingArea of the element (using Heron's formula)                                      
    double restingArea=std::sqrt( ( restingLength[0]+restingLength[1]+restingLength[2])*
                                  (-restingLength[0]+restingLength[1]+restingLength[2])*
                                  ( restingLength[0]-restingLength[1]+restingLength[2])*
                                  ( restingLength[0]+restingLength[1]-restingLength[2])  )*0.25;
    TotalCellArea +=restingArea;
    
    //Current shape local coordinate of the element  (counterclockwise ordering of nodes/edges)
    double CurrentAngle1=std::acos(  (length[0]*length[0]+length[1]*length[1]-length[2]*length[2])/
                                     (length[0]*length[1]*2)    );
    
    double Qa=std::cos(CurrentAngle1)*length[0];
    double Qc=std::sin(CurrentAngle1)*length[0];
    double Qb=length[1];
    
    double positionLocal[3][2]={ {Qa , Qc}, 
                                 {0  , 0 },  
                                 {Qb , 0 }  };
    
    // Local coordinates of the resting shape ( counterclockwise )
    double RestingAngle1=std::acos(  (restingLength[0]*restingLength[0]
                                      +restingLength[1]*restingLength[1]
                                      -restingLength[2]*restingLength[2])/
                                     (restingLength[0]*restingLength[1]*2)    );
    
    double Pa=std::cos(RestingAngle1)*restingLength[0];
    double Pc=std::sin(RestingAngle1)*restingLength[0];
    double Pb=restingLength[1];
    
    // shape vector matrix in resting shape in local coordinate system  
    //= inverse of coordinate matrix ( only first two elements i.e. ShapeVectorResting[3][2] )      
    double ShapeVectorResting[3][3]={ {  0   ,       1/Pc      , 0 }, 
                                      {-1/Pb , (Pa-Pb)/(Pb*Pc) , 1 },       
                                      { 1/Pb ,     -Pa/(Pb*Pc) , 0 }  };
    
    // Rotation Matrix for changing coordinate systems for both Local to Global( Strain Tensor) 
    // and Global to Local( Aniso Vector in the current shape)
    double rotation[3][3];  
    
    double tempA=std::sqrt((position[2][0]-position[1][0])*(position[2][0]-position[1][0])+
                           (position[2][1]-position[1][1])*(position[2][1]-position[1][1])+
                           (position[2][2]-position[1][2])*(position[2][2]-position[1][2])  );
    
    double tempB=std::sqrt((position[0][0]-position[1][0])*(position[0][0]-position[1][0])+
                         (position[0][1]-position[1][1])*(position[0][1]-position[1][1])+
                           (position[0][2]-position[1][2])*(position[0][2]-position[1][2])  );
    
    
    double Xcurrent[3];      
    Xcurrent[0]= (position[2][0]-position[1][0])/tempA;
    Xcurrent[1]= (position[2][1]-position[1][1])/tempA;
    Xcurrent[2]= (position[2][2]-position[1][2])/tempA;
    
    double Bcurrent[3];      
    Bcurrent[0]= (position[0][0]-position[1][0])/tempB;
    Bcurrent[1]= (position[0][1]-position[1][1])/tempB;
    Bcurrent[2]= (position[0][2]-position[1][2])/tempB;
    
    double Zcurrent[3];      
    Zcurrent[0]= Xcurrent[1]*Bcurrent[2]-Xcurrent[2]*Bcurrent[1];
    Zcurrent[1]= Xcurrent[2]*Bcurrent[0]-Xcurrent[0]*Bcurrent[2];
    Zcurrent[2]= Xcurrent[0]*Bcurrent[1]-Xcurrent[1]*Bcurrent[0];
    
    tempA=std:: sqrt(Zcurrent[0]*Zcurrent[0]+Zcurrent[1]*Zcurrent[1]+Zcurrent[2]*Zcurrent[2]);
    Zcurrent[0]=Zcurrent[0]/tempA;
    Zcurrent[1]=Zcurrent[1]/tempA;
    Zcurrent[2]=Zcurrent[2]/tempA;
    
    double Ycurrent[3];      
    Ycurrent[0]= Zcurrent[1]*Xcurrent[2]-Zcurrent[2]*Xcurrent[1];
    Ycurrent[1]= Zcurrent[2]*Xcurrent[0]-Zcurrent[0]*Xcurrent[2];
    Ycurrent[2]= Zcurrent[0]*Xcurrent[1]-Zcurrent[1]*Xcurrent[0];
    
    
    rotation[0][0]=Xcurrent[0];
    rotation[1][0]=Xcurrent[1];
    rotation[2][0]=Xcurrent[2];
    
    rotation[0][1]=Ycurrent[0];
    rotation[1][1]=Ycurrent[1];
    rotation[2][1]=Ycurrent[2];
    
    rotation[0][2]=Zcurrent[0];
    rotation[1][2]=Zcurrent[1];
    rotation[2][2]=Zcurrent[2];  
    
    
    // STRAIN:
    
    double DeformGrad[2][2]={{0,0},{0,0}}; // F=Sigma i Qi x Di ----------------------
    for ( int ii=0 ; ii<3 ; ++ii ) {
      DeformGrad[0][0]=DeformGrad[0][0]+positionLocal[ii][0]*ShapeVectorResting[ii][0];
      DeformGrad[1][0]=DeformGrad[1][0]+positionLocal[ii][1]*ShapeVectorResting[ii][0];
      DeformGrad[0][1]=DeformGrad[0][1]+positionLocal[ii][0]*ShapeVectorResting[ii][1];
      DeformGrad[1][1]=DeformGrad[1][1]+positionLocal[ii][1]*ShapeVectorResting[ii][1];
    } 
    
    
    double LeftCauchy[2][2]; // B=FFt
    LeftCauchy[0][0]=DeformGrad[0][0]*DeformGrad[0][0]+DeformGrad[0][1]*DeformGrad[0][1];
    LeftCauchy[1][0]=DeformGrad[1][0]*DeformGrad[0][0]+DeformGrad[1][1]*DeformGrad[0][1];
    LeftCauchy[0][1]=DeformGrad[0][0]*DeformGrad[1][0]+DeformGrad[0][1]*DeformGrad[1][1];
    LeftCauchy[1][1]=DeformGrad[1][0]*DeformGrad[1][0]+DeformGrad[1][1]*DeformGrad[1][1];
    
    double StrainAlmansi[2][2]; // e=0.5(1-B^-1)  True strain tensor
    double temp=LeftCauchy[0][0]*LeftCauchy[1][1]-LeftCauchy[1][0]*LeftCauchy[0][1]; // det(B)
    StrainAlmansi[0][0]=0.5*(1-(LeftCauchy[1][1]/temp));
    StrainAlmansi[1][0]=0.5*LeftCauchy[1][0]/temp;
    StrainAlmansi[0][1]=0.5*LeftCauchy[0][1]/temp;  
    StrainAlmansi[1][1]=0.5*(1-(LeftCauchy[0][0]/temp));
  
    
    double StrainTensor[3][3]; // there are other alternatives than StrainAlmansi for strain tensor
    StrainTensor[0][0]=StrainAlmansi[0][0];
    StrainTensor[1][0]=StrainAlmansi[1][0];
    StrainTensor[0][1]=StrainAlmansi[0][1];
    StrainTensor[1][1]=StrainAlmansi[1][1];
    
    
    StrainTensor[0][2]=0;  // adding 3rd dimension which is zero, the tensor is still in element plane
    StrainTensor[1][2]=0;
    StrainTensor[2][2]=0;
    StrainTensor[2][0]=0;
    StrainTensor[2][1]=0;
    
    
    // rotating strain tensor to the global coordinate system
    double tempR[3][3]={{0,0,0},{0,0,0},{0,0,0}};
    for (int r=0 ; r<3 ; r++) 
      for (int s=0 ; s<3 ; s++) 
        for(int w=0 ; w<3 ; w++) 
          tempR[r][s] += rotation[r][w]*StrainTensor[w][s];
    
    for (int r=0 ; r<3 ; r++) 
    for (int s=0 ; s<3 ; s++) {
      StrainTensor[r][s]=0;
      for(int w=0 ; w<3 ; w++) 
        StrainTensor[r][s] += tempR[r][w]*rotation[s][w]; 
    }
    
    for (int r=0 ; r<3 ; r++) 
      for (int s=0 ; s<3 ; s++)    
        StrainCellGlobal[r][s] += restingArea*StrainTensor[r][s];
  }
  
  for (int r=0 ; r<3 ; r++) 
    for (int s=0 ; s<3 ; s++)    
      StrainCellGlobal[r][s]= StrainCellGlobal[r][s]/TotalCellArea; 
  


  
  
  strainEpcilon =0.000001;
  for(size_t gg=0; gg<3;gg++)  
    for(size_t ggk=0; ggk<3;ggk++)
      if(gg==ggk)  
        eigenVectorStrain[gg][ggk]=1;
      else
        eigenVectorStrain[gg][ggk]=0;
  
  pivot=1;
  pi=3.1415;
  
   //double RotAngle,Si,Co;
  
  
  
  pivot=std::abs(StrainCellGlobal[1][0]);
  I=1;
  J=0;
  if (std::abs(StrainCellGlobal[2][0])>pivot) {
    pivot=std::abs(StrainCellGlobal[2][0]);
    I=2;
    J=0;
  }
  if (std::abs(StrainCellGlobal[2][1])>pivot) {
    pivot=std::abs(StrainCellGlobal[2][1]);
    I=2;
    J=1;
  }
  
  
  while (pivot> strainEpcilon) {
    
    if (std::abs(StrainCellGlobal[I][I]-StrainCellGlobal[J][J])< strainEpcilon ) {
      RotAngle=pi/4;
    }            
    else {
      RotAngle=0.5*std::atan((2*StrainCellGlobal[I][J])/(StrainCellGlobal[J][J]-StrainCellGlobal[I][I]));
    }
    Si=std::sin(RotAngle);
    Co=std::cos(RotAngle);
    double tempRot[3][3]={{1,0,0},{0,1,0},{0,0,1}};
    tempRot[I][I]=Co;
    tempRot[J][J]=Co;
    tempRot[I][J]=Si;
    tempRot[J][I]=-Si;
    double tempStrain[3][3]={{0,0,0},{0,0,0},{0,0,0}};
    for (int r=0 ; r<3 ; r++) 
      for (int s=0 ; s<3 ; s++) 
        for(int w=0 ; w<3 ; w++) 
          tempStrain[r][s]=tempStrain[r][s]+StrainCellGlobal[r][w]*tempRot[w][s];
    
    for (int r=0 ; r<3 ; r++) 
      for (int s=0 ; s<3 ; s++) 
        StrainCellGlobal[r][s]=0;
    
    for (int r=0 ; r<3 ; r++) 
      for (int s=0 ; s<3 ; s++) 
        for(int w=0 ; w<3 ; w++) 
          StrainCellGlobal[r][s]=StrainCellGlobal[r][s]+tempRot[w][r]*tempStrain[w][s];
    
    for (int r=0 ; r<3 ; r++) 
      for (int s=0 ; s<3 ; s++) 
        tempStrain[r][s]=eigenVectorStrain[r][s];
    
    for (int r=0 ; r<3 ; r++) 
      for (int s=0 ; s<3 ; s++) 
        eigenVectorStrain[r][s]=0;
    
    for (int r=0 ; r<3 ; r++) 
      for (int s=0 ; s<3 ; s++) 
        for(int w=0 ; w<3 ; w++) 
          eigenVectorStrain[r][s]=eigenVectorStrain[r][s]+tempStrain[r][w]*tempRot[w][s];
    
    
    pivot=std::abs(StrainCellGlobal[1][0]);
    I=1;
    J=0;
    if (std::abs(StrainCellGlobal[2][0])>pivot) {
      pivot=std::abs(StrainCellGlobal[2][0]);
      I=2;
      J=0;
    }
    if (std::abs(StrainCellGlobal[2][1])>pivot) {
      pivot=std::abs(StrainCellGlobal[2][1]);
      I=2;
      J=1;
    }
  }
  
  // normalizing eigenvectors (remove if not necessary)  
  temp=std::sqrt(eigenVectorStrain[0][0]*eigenVectorStrain[0][0] +
                 eigenVectorStrain[1][0]*eigenVectorStrain[1][0] +
                 eigenVectorStrain[2][0]*eigenVectorStrain[2][0] );
  if(temp>0){
    eigenVectorStrain[0][0]/=temp;
    eigenVectorStrain[1][0]/=temp;
    eigenVectorStrain[2][0]/=temp;
  }
  temp=std::sqrt(eigenVectorStrain[0][1]*eigenVectorStrain[0][1] +
                 eigenVectorStrain[1][1]*eigenVectorStrain[1][1] +
                 eigenVectorStrain[2][1]*eigenVectorStrain[2][1] );
  if(temp>0){
    eigenVectorStrain[0][1]/=temp;
    eigenVectorStrain[1][1]/=temp;
    eigenVectorStrain[2][1]/=temp;
  }
  temp=std::sqrt(eigenVectorStrain[0][2]*eigenVectorStrain[0][2] +
                 eigenVectorStrain[1][2]*eigenVectorStrain[1][2] +
                 eigenVectorStrain[2][2]*eigenVectorStrain[2][2] );
  if(temp>0){
    eigenVectorStrain[0][2]/=temp;
    eigenVectorStrain[1][2]/=temp;
    eigenVectorStrain[2][2]/=temp;
  }
  
  // maximal strain direction
  maximalStrainValue=StrainCellGlobal[0][0];
  Istrain=0;
  if (StrainCellGlobal[1][1]>maximalStrainValue) 
    {
      maximalStrainValue=StrainCellGlobal[1][1];
      Istrain=1;
    }
  if (StrainCellGlobal[2][2]>maximalStrainValue) 
    {
      maximalStrainValue=StrainCellGlobal[2][2];
      Istrain=2;
    }
  
  
  // 2nd maximalstrain direction/value
  //double maximalStrainValue2,maximalStrainValue3;
  //int Istrain2,Istrain3;
  if (Istrain==0) {
    Istrain2=1;
    Istrain3=2;
  }
  if (Istrain==1) {
    Istrain2=0;
    Istrain3=2;
  }
  if (Istrain==2) {
    Istrain2=0;
    Istrain3=1;
  }
  if(StrainCellGlobal[Istrain3][Istrain3] > StrainCellGlobal[Istrain2][Istrain2] ) {
    temp=Istrain2;
    Istrain2=Istrain3;
    Istrain3=temp;
  }
  maximalStrainValue2=StrainCellGlobal[Istrain2][Istrain2];
  maximalStrainValue3=StrainCellGlobal[Istrain3][Istrain3];
  std::cerr<<"  strain after 2 "<<maximalStrainValue<<"  "<<maximalStrainValue2<<"  "<<maximalStrainValue3<<"  "<<std::endl;

  ///////////////////////////////////////////////////////// end ad-hoc

}



void Tissue::removeTwoVertex( size_t index ) 
{
	if (vertex(index).numWall() != 2) {
		std::cerr << "Tissue::removeTwoVertex() Vertex not a two-vertex, not removed!" 
							<< std::endl;
		return;
	}
	Vertex* v=vertexP(index);
	Wall* w1=vertex(index).wall(0);
	Wall* w2=vertex(index).wall(1);
	assert(v->numCell()==1 || v->numCell()==2);
	Cell* c1=vertex(index).cell(0);
	Cell* c2;
	if (v->numCell()==2)
		c2=vertex(index).cell(1);
	else
		c2=background();
	// Get 'neighboring' vertices from walls
	Vertex* v1=w1->vertex1();
	if (v1==v)
		v1=w1->vertex2();
	Vertex* v2=w2->vertex1();
	if (v2==v)
		v2=w2->vertex2();
	
	// Remove v,w2 from c1 and c2 (if not background)
	if (c1!=background()) {
		size_t numW=c1->numWall();
		size_t numV=c1->numVertex();
		assert(numW>3 && numW==numV);
		std::vector<Vertex*> newV(numV-1);
		std::vector<Wall*> newW(numV-1);
		size_t vI=0;
		size_t wI=0;
		for (size_t i=0; i<numV; ++i) {
			if (c1->vertex(i) != v)
				newV[vI++]=c1->vertex(i);
			if (c1->wall(i) != w2)
				newW[wI++]=c1->wall(i);
		}
		if (vI!=numV-1 || wI!=numW-1) {
			std::cerr << "Tissue::removeTwoVertex Vertex or wall to be removed not found in cell."
								<< std::endl;
			std::cerr << "Cell: " << c1->index() << " NumVertex(wall): " << c1->numVertex() << " ("
								<< c1->numWall() << ") Vertex: " << v->index() << " Wall: "
								<< w2->index() << " vI: " << vI << " wI: " << wI << std::endl;
			exit(-1);
		}
		c1->setVertex(newV);
		c1->setWall(newW);
	}		
	if (c2!=background()) {
		size_t numW=c2->numWall();
		size_t numV=c2->numVertex();
		assert(numW>3 && numW==numV);
		std::vector<Vertex*> newV(numV-1);
		std::vector<Wall*> newW(numV-1);
		size_t vI=0;
		size_t wI=0;
		for (size_t i=0; i<numV; ++i) {
			if (c2->vertex(i) != v)
				newV[vI++]=c2->vertex(i);
			if (c2->wall(i) != w2)
				newW[wI++]=c2->wall(i);
		}
		if (vI!=numV-1 || wI!=numW-1) {
			std::cerr << "Tissue::removeTwoVertex Vertex or wall to be removed not found in cell."
								<< std::endl;
			std::cerr << "Cell: " << c2->index() << " NumVertex(wall): " << c2->numVertex() << " ("
								<< c2->numWall() << ") Vertex: " << v->index() << " Wall: "
								<< w2->index() << " vI: " << vI << " wI: " << wI << std::endl;
			exit(-1);
		}
		c2->setVertex(newV);
		c2->setWall(newW);
	}		
	// Update w1
	if (w1->vertex1()==v && w1->vertex2()==v1)
		w1->setVertex1(v2);
	else if (w1->vertex2()==v && w1->vertex1()==v1)
		w1->setVertex2(v2);
	else {
		std::cerr << "Tissue::removeTwoVertex() Wrong in updating wall." << std::endl; 
		exit(-1);
	}
	// update w1 length?
	// Update v2 to connect to w1 instead of w2
	size_t numW=v2->numWall();
	size_t updateFlag=0;
	for (size_t i=0; i<numW; ++i)
		if (v2->wall(i)==w2) {
			v2->setWall(i,w1);
			++updateFlag;
		}
	if (updateFlag!=1) {
		std::cerr << "Tissue::removeTwoVertex() Update of v2 wrong." << std::endl;
		exit(-1);
	}
	// Remove v and w2
	removeVertex(v->index());
	removeWall(w2->index());
}

void Tissue::sortCellWallAndCellVertex(Cell *cell) 
{	
	std::vector<size_t> sortedFlag(numCell());
	size_t numSorted=0;
	if( !cell ) { //sort all		
	  size_t startSortIndex=0;
	  do {
	    if (!sortedFlag[startSortIndex]) {
	      cell = this->cellP(startSortIndex);	  
	      sortCellRecursive(cell, sortedFlag, numSorted);
	    }
	    startSortIndex++;
	    //std::cerr << "Tissue::sortCellWallAndCellVertex() numSorted = " << numSorted 
	    //	      << " startSortIndex = " << startSortIndex << std::endl;
	  } while (numSorted<numCell() && startSortIndex<numCell());
	}
	else { //sort recursively around given cell (at division?)
	  sortCellRecursive(cell, sortedFlag, numSorted);
	}
	std::cerr << "Tissue::sortCellWallAndCellVertex() " << numSorted << " of " << numCell()
		  << " faces (cells) sorted recursively." << std::endl;
}

void Tissue::sortCellRecursive( Cell* cell, std::vector<size_t> &sortedFlag, size_t &numSorted)
{
  if (sortedFlag[cell->index()])
    return;
  cell->sortWallAndVertex(*this);
  sortedFlag[cell->index()]++;
  numSorted++;
  for (size_t k=0; k<cell->numWall(); ++k) {
    Cell *cellNext = cell->cellNeighbor(k);
    if (cellNext!=background())
      sortCellRecursive(cellNext,sortedFlag,numSorted);
  }
  return;
}

void Tissue::checkConnectivity(size_t verbose) 
{	
  int exitFlag=0;
  size_t numC = numCell();
  // Check if all indices are used
  //
  for (size_t i=0; i<numC; ++i) {
    if( verbose ) {
      if( cell(i).index() != i ) {
	std::cerr << "Tissue::checkConnectivity() "
		  << "Cell " << i << " has index " << cell(i).index()
		  << std::endl;
	exitFlag++;
      }
    }
    else
      assert( cell(i).index() == i );
  }
  for (size_t i=0; i<numWall(); ++i) {
    if( verbose ) {
      if( wall(i).index() != i ) {
	std::cerr << "Tissue::checkConnectivity() "
		  << "Wall " << i << " has index " << wall(i).index()
		  << std::endl;
	exitFlag++;
      }
    }
    else
      assert( wall(i).index() == i );
  }
  for (size_t i=0; i<numVertex(); ++i) {
    if( verbose ) {
      if( vertex(i).index() != i ) {
	std::cerr << "Tissue::checkConnectivity() "
		  << "Vertex " << i << " has index " << vertex(i).index()
		  << std::endl;
	exitFlag++;
      }
    }
    else
      assert( vertex(i).index() == i );
  }
  // Make sure all cellVertex(Wall) are real vertices(walls) via index
  //
  for( size_t k=0 ; k<numC ; ++k ) {
    for( size_t l=0 ; l<cell(k).numWall() ; ++l ) { 
      if( verbose ) {
	if( cell(k).wall(l)->index()>=numWall() ) {
	  std::cerr << "Tissue::checkConnectivity() " << "Cell " << k 
		    << " is connected to wall "
		    << cell(k).wall(l)->index() << "("
		    << numWall() << " walls in total)" << std::endl;
	  exitFlag++;
	}
      }
      else {
	assert( cell(k).wall(l)->index()<numWall() );
      }
    }
    for( size_t l=0 ; l<cell(k).numVertex() ; ++l ) { 
      if( verbose ) {
	if( cell(k).vertex(l)->index()>=numVertex() ) {
	  std::cerr << "Tissue::checkConnectivity() " << "Cell " << k 
		    << " is connected to vertex "
		    << cell(k).vertex(l)->index() << "("
		    << numVertex() << " vertices in total)" << std::endl;
	  exitFlag++;
	}
	for( size_t ll=l+1 ; ll<cell(k).numVertex() ; ++ll ) { 
	  if( cell(k).vertex(l)==cell(k).vertex(ll) ) {
	    std::cerr << "Tissue::checkConnectivity() " << "Cell " << k 
		      << " is connected to vertex "
		      << cell(k).vertex(l)->index() << "(twice)"
		      << std::endl;
	    exitFlag++;
	  }
	}
      }
      else {
	assert( cell(k).vertex(l)->index()<numVertex() );
      }
    }
  }
  //Make sure all wallVertex(Cell) are real vertices(cells) via index
  //
  for( size_t k=0 ; k<numWall() ; ++k ) {
    if( verbose ) {
      if( ( wall(k).cell1()->index()>=numCell() &&
	    wall(k).cell1() != background() ) ||
	  ( wall(k).cell2()->index()>=numCell() &&
	    wall(k).cell2() != background() ) ) {
	std::cerr << "Tissue::checkConnectivity() " << "Wall " << k 
		  << " is connected to cell "
		  << wall(k).cell1()->index() << " and "
		  << wall(k).cell2()->index() << " ("
		  << numCell() << " cells in total)" << std::endl;
	exitFlag++;
      }
      if( wall(k).cell1() == wall(k).cell2() ) {
	std::cerr << "Tissue::checkConnectivity() " << "Wall " << k 
		  << " is connected to cell "
		  << wall(k).cell1()->index() << " and "
		  << wall(k).cell2()->index() << " (same cell)" 
		  << std::endl;
	exitFlag++;
      }
    }
    else {
      assert( ( wall(k).cell1()->index()<numCell() ||
		wall(k).cell1() == background() ) &&
	      ( wall(k).cell2()->index()<numCell() ||
		wall(k).cell2() == background() ) );
      assert( wall(k).cell1() != wall(k).cell2() );
    }
    if( verbose ) {
      if( wall(k).vertex1()->index()>=numVertex() ||
	  wall(k).vertex2()->index()>=numVertex() ) {
	std::cerr << "Tissue::checkConnectivity() " << "Wall " << k 
		  << " is connected to vertex "
		  << wall(k).vertex1()->index() << " and "
		  << wall(k).vertex2()->index() << " ("
		  << numVertex() << " vertices in total)" << std::endl;
	exitFlag++;
      }
      if( wall(k).vertex1() == wall(k).vertex2() ) {
	std::cerr << "Tissue::checkConnectivity() " << "Wall " << k 
		  << " is connected to vertex "
		  << wall(k).vertex1()->index() << " and "
		  << wall(k).vertex2()->index() << " (same vertex)" 
		  << std::endl;
	exitFlag++;
      }
    }
    else {
      assert( wall(k).vertex1()->index()<numVertex() &&
	      wall(k).vertex2()->index()<numVertex() );			
      assert( wall(k).vertex1() != wall(k).vertex2() );
    }
  }
  
  //Make sure all vertexCell(Wall) are real cells(walls) via index
  //
  for( size_t k=0 ; k<numVertex() ; ++k ) {
    for( size_t l=0 ; l<vertex(k).numCell() ; ++l ) { 
      if( verbose ) {
	if( vertex(k).cell(l)->index()>=numCell() ) {
	  std::cerr << "Tissue::checkConnectivity() " << "Vertex " << k 
		    << " is connected to cell "
		    << vertex(k).cell(l)->index() << "("
		    << numCell() << " cells in total)" << std::endl;
	  exitFlag++;
	}
	if( vertex(k).cell(l) == background() ) {
	  std::cerr << "Tissue::checkConnectivity() " 
		    << "Vertex " << k << " is connected to background"
		    << std::endl;
	  exitFlag++;
	}
	for( size_t ll=l+1 ; ll<vertex(k).numCell() ; ++ll ) { 
	  if( vertex(k).cell(l)==vertex(k).cell(ll) ) {
	    std::cerr << "Tissue::checkConnectivity() " << "Vertex " << k 
		      << " is connected to cell "
		      << vertex(k).cell(l)->index() << " twice."
		      << std::endl;
	    exitFlag++;
	  }
	}		
      }
      else {
	assert( vertex(k).cell(l)->index()<numCell() );
	assert( vertex(k).cell(l) != background() );
      }
    }
    for( size_t l=0 ; l<vertex(k).numWall() ; ++l ) { 
      if( verbose ) {
	if( vertex(k).wall(l)->index()>=numWall() ) {
	  std::cerr << "Tissue::checkConnectivity() " << "Vertex " << k 
		    << " is connected to wall "
		    << vertex(k).wall(l)->index() << "("
		    << numWall() << " walls in total)" << std::endl;
	  exitFlag++;
	}
	for( size_t ll=l+1 ; ll<vertex(k).numWall() ; ++ll ) { 
	  if( vertex(k).wall(l)==vertex(k).wall(ll) ) {
	    std::cerr << "Tissue::checkConnectivity() " << "Vertex " << k 
		      << " is connected to wall "
		      << vertex(k).wall(l)->index() << " twice."
		      << std::endl;
	    exitFlag++;
	  }
	}		
      }
      else {
	assert( vertex(k).wall(l)->index()<numWall() );
      }
    }
  }
  
  // Make sure all compartments include same set of variables
  //
  size_t numVarTmp=cell(0).numVariable();
  for (size_t i=1; i<numC; ++i) {
    if( verbose ) {
      if( numVarTmp != cell(i).numVariable() ) {
	std::cerr << "Tissue::checkConnectivity() " << "Cell " << i 
		  << " has " << cell(i).numVariable() << " variables" 
		  << " while cell 0 has " << numVarTmp << std::endl;
	exitFlag++;
      }
    }
    else
      assert( numVarTmp==cell(i).numVariable() );
  }
  // Do checks on connectivity from cells
  //
  for (size_t i=0; i<numC; ++i) {
    
    if( verbose ) {
      if ( cell(i).numWall() != cell(i).numVertex() ) {
	std::cerr << "Tissue::checkConnectivity() "
		  << "Cell " << i << " has " << cell(i).numWall()
		  << " walls and " << cell(i).numVertex() 
		  << " vertices!" << std::endl; 
	exitFlag++;
      }
    }
    else
      assert( cell(i).numWall() == cell(i).numVertex() );
    
    //Make sure that vertecis in all cell-walls are cell-vertices
    for (size_t w=0 ; w<cell(i).numWall(); ++w) {
      if ( verbose ) {
	if ( !cell(i).hasVertex( cell(i).wall(w)->vertex1() ) ) {
	  std::cerr << "Tissue::checkConnectivity() "
		    << "Cell " << i << " has wall " << cell(i).wall(w)->index()
		    << " with vertex " << cell(i).wall(w)->vertex1()->index()
		    << " but is not connected to the vertex!"
		    << std::endl;
	  exitFlag++;
	}
	if ( !cell(i).hasVertex( cell(i).wall(w)->vertex2() ) ) {
	  std::cerr << "Tissue::checkConnectivity() "
		    << "Cell " << i << " has wall " << cell(i).wall(w)->index() 
		    << " with vertex " << cell(i).wall(w)->vertex2()->index()
		    << " but is not connected to the vertex!"
		    << std::endl;
	  exitFlag++;
	}				
      }
      else {
	assert( cell(i).hasVertex( cell(i).wall(w)->vertex1() ) );
	assert( cell(i).hasVertex( cell(i).wall(w)->vertex2() ) );
      }
    }
    //Make sure that two walls in all cell-vertices are cell-walls
    for (size_t v=0 ; v<cell(i).numVertex(); ++v) {
      int numWall=0;
      for (size_t w=0 ; w<cell(i).vertex(v)->numWall(); ++w) {
	numWall += cell(i).hasWall( cell(i).vertex(v)->wall(w) );
      }
      if ( verbose ) {
	if( numWall != 2 ) {
	  std::cerr << "Tissue::checkConnectivity() "
		    << "Cell " << i << " has vertex " << cell(i).vertex(v)->index() 
		    << " with " << cell(i).vertex(v)->numWall() 
		    << " walls but " << numWall 
		    << " walls are connected to the cell!"
		    << std::endl;
	  exitFlag++;
	}								
      }
      else {
	assert( numWall==2 );
      }
    }
  }//for i
  
  // Check that walls and vertices are properly sorted in cells
  //
  for (size_t i=0; i<numC; ++i) {
    Cell* cP = cellP(i);
    size_t numW=cP->numWall();
    for (size_t k=0; k<numW; ++k) {
      size_t kPlus = (k+1)%numW;
      if (cP->wall(k)->cell1()==cP) {
	if (cP->wall(k)->cellSort1()==-1) {
	  if (cP->vertex(kPlus) != cP->wall(k)->vertex1() ||
	      cP->vertex(k) != cP->wall(k)->vertex2() ) {
	    std::cerr << "Tissue::checkConnectivity() "
		      << "1: vertices and walls not sorted correctly in cell " 
		      << i << " wall " << cP->wall(k)->index() << std::endl;
	    ++exitFlag;
	  }
	}
	else { //cellSort=1 or cellSort=0
	  if (cP->vertex(k) != cP->wall(k)->vertex1() ||
	      cP->vertex(kPlus) != cP->wall(k)->vertex2() ) {
	    std::cerr << "Tissue::checkConnectivity() "
		      << "2: vertices and walls not sorted correctly in cell "
		      << i << " wall " << cP->wall(k)->index() << std::endl;
	    ++exitFlag;
	  }
	}
      }
      else if (cP->wall(k)->cell2()==cP) {
	if (cP->wall(k)->cellSort2()==-1) {
	  if (cP->vertex(kPlus) != cP->wall(k)->vertex1() ||
	      cP->vertex(k) != cP->wall(k)->vertex2() ) {
	    std::cerr << "Tissue::checkConnectivity() "
		      << "3: vertices and walls not sorted correctly in cell "
		      << i << " wall " << cP->wall(k)->index() << std::endl;
	    ++exitFlag;
	  }
	}
	else { //cellSort=1 or cellSort=0
	  if (cP->vertex(k) != cP->wall(k)->vertex1() ||
	      cP->vertex(kPlus) != cP->wall(k)->vertex2() ) {
	    std::cerr << "Tissue::checkConnectivity() "
		      << "4: vertices and walls not sorted correctly in cell "
		      << i << " wall " << cP->wall(k)->index() << std::endl;
	    ++exitFlag;
	  }
	}
      }
      else {
	std::cerr << "Tissue::checkConnectivity() "
		  << "cellWall not connected to cell." << std::endl; 
	++exitFlag;
      }
    }
  }
  
  if ( exitFlag ) {
    std::cerr << "Tissue::checkConnectivity() "
	      << exitFlag << " errors found in tissue." << std::endl;
    exit(-1);
  }
}

unsigned int Tissue::
findPeaksGradientAscent( DataMatrix &cellData, 
			 size_t col, std::vector<size_t> &cellMax,
			 std::vector<size_t> &flag )
{
  assert(cellData.size() == numCell() );
  assert( cellData[0].size()>col );
  
  if( cellMax.size() )
    cellMax.resize(0);
  if( flag.size() != numCell() ) flag.resize(numCell());
  for( size_t i=0 ; i<numCell() ; ++i )
    flag[i]=0;
  
  std::vector<size_t> cellTmp;//Values before threshold check
  std::vector<unsigned int> numTmp;//times cellTmp been visited
  std::vector<size_t> walkTmp;//positions for a walk (start point)
  
  size_t count=1;
  //Find the maxima from each cell
  for( size_t iStart=0 ; iStart<numCell() ; ++iStart ) {
    size_t i=iStart;
    double value,newValue;
    walkTmp.resize(1);
    walkTmp[0]=i;
    //find the max by walking uphill (greedy)
    if( !flag[i] ) {
      do {
	newValue=value=cellData[i][col];
	size_t newI=i;
	//Check all neighboring cells
	for(size_t k=0 ; k<cell(i).numWall() ; k++ ) {
	  size_t j = cell(i).wall(k)->cell1()->index();
	  if( j==i )
	    j = cell(i).wall(k)->cell2()->index();
	  if( j != background()->index() ) {
	    if( cellData[j][col]>newValue ) {
	      newValue=cellData[j][col];
	      newI=j;
	    }
	  }
	}
	i=newI;
	walkTmp.push_back( i );
      } while( newValue>value && !flag[i] );
    }
    //Collect the path data and add one visit for the maximum
    if( !flag[i] ) { //new maximum
      cellTmp.push_back( i );
      numTmp.push_back(1);
      unsigned int n=count++;//cellTmp.size();
      for( size_t a=0 ; a<walkTmp.size() ; a++ )
	flag[ walkTmp[a] ] = n;
    }
    else { //old maximum or background
      size_t n = flag[i];
      for( size_t a=0 ; a<walkTmp.size() ; a++ )
	flag[ walkTmp[a] ] = n;
      if( flag[i]>0 )//old maxima
	numTmp[n-1]++;
    }
  }
  //No threshold checking...
  unsigned int threshold = 1;
  double valThreshold = 0.0;
  //Get the maxima visited more than threshold times and with an intensity
  //value higher than threshold
  std::vector<int> clusterNum;
  for( size_t n=0 ; n<cellTmp.size() ; n++ )
    if( numTmp[n]>=threshold && cellData[ cellTmp[n] ][col]>valThreshold ) {
      cellMax.push_back( cellTmp[n] ); 
      clusterNum.push_back( n+1 );
    }
  
  //Save the basins of attraction
  //   boa.resize( cellMax.size() );
  //   for( int i=0 ; i<numCompartment() ; i++ )
  //     for( int n=0 ; n<cellMax.size() ; n++ )
  //       if( flag[i] == clusterNum[n] )
  // 	boa[n].push_back(i);
  
  return static_cast<unsigned int>(cellMax.size());
}

void Tissue::printInit(std::ostream &os) const {
  
  // Increase resolution to max for doubles
  unsigned int oldPrecision = os.precision(); 
  os.precision(15);
  //std::cerr << "Tissue::prinitInit(): old precision: " << oldPrecision << " new " 
  //				<< os.precision() << std::endl;
  
  os << numCell() << " " << numWall() << " " << numVertex() << std::endl;
  
  //Print the connectivity from walls
  for( size_t i=0 ; i<numWall() ; ++i ) {
    os << i << " ";
    if( wall(i).cell1()->index()<numCell() )
      os << wall(i).cell1()->index() << " " ;
    else
      os << "-1 ";
    if( wall(i).cell2()->index()<numCell() )
      os << wall(i).cell2()->index() << " ";
    else
      os << "-1 ";
    os << wall(i).vertex1()->index() 
       << " " << wall(i).vertex2()->index() << std::endl;
  }
  os << std::endl;
  
  //Print the vertex positions
  os << numVertex() << " " << vertex(0).numPosition() << std::endl;
  for( size_t i=0 ; i<numVertex() ; ++i ) {
    for( size_t j=0 ; j<vertex(i).numPosition() ; ++j )
      os << vertex(i).position(j) << " ";
    os << std::endl;
  }
  os << std::endl;
  
  //Print wall data
  os << numWall() << " 1 " << wall(0).numVariable() << std::endl;
  for( size_t i=0 ; i<numWall() ; ++i ) {
    os << wall(i).length() << " ";
    for( size_t j=0 ; j<wall(i).numVariable() ; ++j )
      os << wall(i).variable(j) << " ";
    os << std::endl;
  }
  os << std::endl;
  
  //Print cell data
  os << numCell() << " " << cell(0).numVariable() << std::endl;
  if( cell(0).numVariable() ) {
    for( size_t i=0 ; i<numCell() ; ++i ) {
      for( size_t j=0 ; j<cell(i).numVariable() ; ++j )
	os << cell(i).variable(j) << " ";
      os << std::endl;
    }
    os << std::endl;
  }  
  os.precision(oldPrecision);		
}

void Tissue::printInit(DataMatrix &cellData,
		       DataMatrix &wallData,
		       DataMatrix &vertexData,
		       std::ostream &os) {
  
  assert( numCell()==cellData.size() && 
	  numWall()==wallData.size() &&
	  numVertex()==vertexData.size() );
  
  // Increase resolution to max for doubles
  unsigned int oldPrecision = os.precision(); 
  os.precision(15);
  //std::cerr << "Tissue::prinitInit(): old precision: " << oldPrecision << " new " 
  //				<< os.precision() << std::endl;
  
  os << numCell() << " " << numWall() << " " << numVertex() << std::endl;
  
  // Print the connectivity from walls
  for( size_t i=0 ; i<numWall() ; ++i ) {
    os << i << " ";
    if( wall(i).cell1()->index()<numCell() )
      os << wall(i).cell1()->index() << " " ;
    else
      os << "-1 ";
    if( wall(i).cell2()->index()<numCell() )
      os << wall(i).cell2()->index() << " ";
    else
      os << "-1 ";
    os << wall(i).vertex1()->index() 
       << " " << wall(i).vertex2()->index() << std::endl;
  }
  os << std::endl;
  
  //Print the vertex positions
  os << numVertex() << " " << vertex(0).numPosition() << std::endl;
  for( size_t i=0 ; i<numVertex() ; ++i ) {
    assert( vertex(i).numPosition()==vertexData[i].size() );
    for( size_t j=0 ; j<vertex(i).numPosition() ; ++j )
      os << vertexData[i][j] << " ";
    os << std::endl;
  }
  os << std::endl;
  
  //Print wall data
  os << numWall() << " 1 " << wall(0).numVariable() << std::endl;
  for( size_t i=0 ; i<numWall() ; ++i ) {
    assert( wallData[i].size() );
    for( size_t j=0 ; j<wallData[i].size() ; ++j )
      os << wallData[i][j] << " ";
    os << std::endl;
  }
  os << std::endl;
  
  //Print cell data
  os << numCell() << " " << cell(0).numVariable() << std::endl;
  if( cell(0).numVariable() ) {
    for( size_t i=0 ; i<numCell() ; ++i ) {
      assert( cellData[i].size() );
      for( size_t j=0 ; j<cellData[i].size() ; ++j )
	os << cellData[i][j] << " ";
      os << std::endl;
    }
    os << std::endl;
  }  
  os.precision(oldPrecision);	
}

void Tissue::printInitFem(std::ostream &os) const
{
  // Increase resolution 
  unsigned int oldPrecision = os.precision(); 
  os.precision(20);
  //std::cerr << "Tissue::printInit(): old precision: " << oldPrecision << " new " 
  //				<< os.precision() << std::endl;	
  
  //Print the vertex positions
  size_t numV=numVertex();
  os << numV << " nodes" << std::endl;
  for (size_t i=0; i<numV; ++i) {
    os << i << " : ";
    for (size_t j=0; j<vertex(i).numPosition(); ++j )
      os << vertex(i).position(j) << " ";
    os << std::endl;
  }
  //os << std::endl;
  
  //Print cell connection data
  size_t numC=numCell();
  os << numC << " faces" << std::endl;
  for (size_t i=0; i<numC; ++i) {
    os << i << " : ";
    size_t numV=cell(i).numVertex();
    os << numV << ", ";
    for (size_t k=0; k<numV; ++k)
      os << cell(i).vertex(k)->index() << " ";
    os << std::endl;
  }
  os.precision(oldPrecision);
}

void Tissue::printInitTri(std::ostream &os)
{
  // Increase resolution 
  unsigned int oldPrecision = os.precision(); 
  os.precision(20);
  //std::cerr << "Tissue::printInitTri(): old precision: " << oldPrecision << " new " 
  //	    << os.precision() << std::endl;	

  // Create data structure for triangulated tissue.
  //
  size_t numC=0; //added below
  size_t numW=numWall(); //appended below
  size_t numV=numVertex()+numCell(); //all
  std::vector<size_t> cellIndexStart(numCell());
  std::vector<size_t> wallIndexStart(numCell());
  for( size_t i=0 ; i<numCell() ; ++i ) {
    numC += cell(i).numWall();
    numW += cell(i).numVertex();
    if (i==0) {
      cellIndexStart[i] = numCell();
      wallIndexStart[i] = numWall();
    }
    else {
      cellIndexStart[i] = cellIndexStart[i-1]+cell(i-1).numWall()-1;
      wallIndexStart[i] = wallIndexStart[i-1]+cell(i-1).numWall();
    }
  }
  DataMatrix c(numC);
  DataMatrix w(numW);
  DataMatrix v(numV);
  for (size_t i=0; i<numVertex(); ++i) {
    v[i] = vertex(i).position();
  }
  std::vector< std::pair<size_t,size_t> > cellNeigh(numW); // Wall connections to cells
  std::vector< std::pair<size_t,size_t> > vertexNeigh(numW); // Wall connections to vertices

  std::vector<double> wallTmpData(1+wall(0).numVariable(),0.0);
  //size_t D=vertex(0).numPosition(); //dimension

  for (size_t i=0; i<numCell(); ++i) {
    size_t numCellVar = cell(i).numVariable();
    // Add central vertex with position
    size_t vI = numVertex()+i;
    v[vI].resize(v[0].size());
    for (size_t d=0; d<v[vI].size(); ++d) {
      v[vI] = cell(i).positionFromVertex(); // Getting the central position of the cell
    }
    for (size_t k=0; k<cell(i).numWall(); ++k) {
      // add cell data in correct cell      
      if (k==0) {
	// old index used
	c[i].resize(numCellVar); 
	for( size_t j=0; j<numCellVar; ++j) {
	  c[i][j] = cell(i).variable(j);
	}
      }
      else {
	// new index for the rest (and store same cell data)
	size_t ii = cellIndexStart[i]+k-1; 
	c[ii].resize(numCellVar); 
	for( size_t j=0; j<numCellVar; ++j) {
	  c[ii][j] = cell(i).variable(j);
	}
      }
      // update current walls including neighborhood
      size_t wI = cell(i).wall(k)->index();
      // temporary wall data should store length and variables 
      w[wI].resize( 1+wall(wI).numVariable() );
      w[wI][0] = wall(wI).length();
      for (size_t j=0; j<wall(wI).numVariable(); ++j) {
	w[wI][1+j] = wall(wI).variable(j);
      }
      // vertex neighbors stay the same
      vertexNeigh[wI].first = cell(i).vertex(k)->index();
      if (k<cell(i).numVertex()-1) {
	vertexNeigh[wI].second = cell(i).vertex(k+1)->index();
      }
      else {
	vertexNeigh[wI].second = cell(i).vertex(0)->index();
      }
      // first find current cell as neighbor
      if (cell(i).wall(k)->cell1()->index() == i) {
	if (k==0) {
	  cellNeigh[wI].first = i;
	}
	else {
	  cellNeigh[wI].first = cellIndexStart[i]+k-1;
	}
      }
      else if (cell(i).wall(k)->cell2()->index() == i) {
	if (k==0) {
	  cellNeigh[wI].second = i;
	}
	else {
	  cellNeigh[wI].second = cellIndexStart[i]+k-1;
	}
      }
      else {
	std::cerr << "Tissue::printInitTri() Error: Cell wall not"
		  << "connected to cell." << std::endl;
	exit(EXIT_FAILURE);
      }
      // also recognize if the wall is connected to background
      if (cell(i).wall(k)->cell1() == background() ) {
	cellNeigh[wI].first = size_t(-1);
      }
      else if (cell(i).wall(k)->cell2() == background() ) {
	cellNeigh[wI].second = size_t(-1);
      }

      // add new wall between vertex and central vertex
      wI = wallIndexStart[i]+k;
      vertexNeigh[wI].first = cell(i).vertex(k)->index();
      vertexNeigh[wI].second = vI;
      // Calculate length of new wall from vertex positions
      double length=0.0;
      size_t vII = vertexNeigh[wI].first;
      for (size_t d=0; d<v[vI].size(); ++d) {
	length += (v[vI][d]-v[vII][d])*(v[vI][d]-v[vII][d]);
      }
      length = std::sqrt(length);
      wallTmpData[0] = length;
      w[wI] = wallTmpData;// length + 0s for the rest of the wall variables.      
      if (k==0) {
	cellNeigh[wI].first = cellIndexStart[i]+
	  cell(i).numWall()-2; // last from added list of new cells
	cellNeigh[wI].second = cell(i).index(); //i
      }
      else {
	cellNeigh[wI].first = cellNeigh[wI-1].second; // second from prev wall (going around)
	cellNeigh[wI].second = cellIndexStart[i]+k-1; // from added list of new cells
      }
    }
  }
  
  // Print the init file with the new data
  //
  os << numC << " " << numW << " " << numV << std::endl;
  
  // Print the connectivity from walls
  for( size_t i=0 ; i<numW ; ++i ) {
    os << i << " "; // index
    if (cellNeigh[i].first<numC) { // cell neighbors
      os << cellNeigh[i].first << " ";
    }
    else {
      os << "-1 ";
    }
    if (cellNeigh[i].second<numC) {
      os << cellNeigh[i].second << " "; 
    }
    else {
      os << "-1 ";
    }
    os << vertexNeigh[i].first << " " << vertexNeigh[i].second << std::endl; // vertex neighbors
  }
  os << std::endl;
  
  // Print the vertex positions
  os << numV << " " << v[0].size() << std::endl;
  for( size_t i=0 ; i<numV ; ++i ) {
    for( size_t j=0 ; j<v[i].size() ; ++j )
      os << v[i][j] << " ";
    os << std::endl;
  }
  os << std::endl;
  
  // Print wall data
  os << numW << " 1 " << w[0].size()-1 << std::endl;
  for( size_t i=0 ; i<numW ; ++i ) {
//     assert( wa[i].size() );
    for( size_t j=0 ; j<w[i].size() ; ++j )
      os << w[i][j] << " ";
    os << std::endl;
  }
  os << std::endl;
  
  // Print cell data
  os << numC << " " << c[0].size() << std::endl;
  if( c[0].size() ) {
    for( size_t i=0 ; i<numC ; ++i ) {
      assert( c[i].size() );
      for( size_t j=0 ; j<c[i].size() ; ++j )
	os << c[i][j] << " ";
      os << std::endl;
    }
    os << std::endl;
  }  
  os.precision(oldPrecision);
}

void Tissue::printInitOrganism(std::ostream &os)
{
  // Increase resolution 
  unsigned int oldPrecision = os.precision(); 
  os.precision(20);
  //std::cerr << "Tissue::printInit(): old precision: " << oldPrecision << " new " 
  //				<< os.precision() << std::endl;	
  
  //Print the cell
  size_t numC=numCell();
  size_t numVar = cell(0).numVariable();
  size_t dimension = vertex(0).numPosition();
  os << numC << " " << numVar+dimension+1 << std::endl;
  for (size_t i=0; i<numC; ++i) {
    std::vector<double> posTmp = cell(i).positionFromVertex();
    if (posTmp.size()!=dimension) {
      std::cerr << "Tissue::printInitOrganism() Wrong number of positional variables: "
		<< posTmp.size() << ", expecting " << dimension << std::endl;
      exit(EXIT_FAILURE);
    }
    for (size_t d=0; d<dimension; ++d) {
      os << posTmp[d] << " ";
    }
    os << cell(i).volume() << " ";
    for (size_t j=0; j<cell(i).numVariable(); ++j )
      os << cell(i).variable(j) << " ";
    os << std::endl;
  }
  os << std::endl;
  
  //Print cell connection data
  os << numC << " 1" << std::endl;//assuming connections and areas stored
  for (size_t i=0; i<numC; ++i) {
    if (cell(i).index() != i) {
      std::cerr << "Tissue::printInitOrganism() Assumes cell indices are same as position in cell vector."
		<< std::endl << "For cell at position " << i << " the index is " << cell(i).index() 
		<< "." << std::endl; 
      exit(EXIT_FAILURE);
    }
    os << i << " ";
    size_t numW=cell(i).numWall();
    // Collect real neighbors
    std::vector<size_t> kList,jList;
    for (size_t k=0; k<numW; ++k) {
      if (cell(i).wall(k)->cell1()->index()==cell(i).index()) {
	if (cell(i).wall(k)->cell2()!=background()) {
	  kList.push_back(k);
	  jList.push_back(cell(i).wall(k)->cell2()->index());
	}
      }
      else if (cell(i).wall(k)->cell2()->index()==cell(i).index()) {
	if (cell(i).wall(k)->cell1()!=background()) {
	  kList.push_back(k);
	  jList.push_back(cell(i).wall(k)->cell1()->index());
	}
      }
    }
    // Print
    os << kList.size() << " ";
    for (size_t k=0; k<jList.size(); ++k) {
      os << jList[k] << " ";
    }
    for (size_t k=0; k<kList.size(); ++k) {
      os << cell(i).wall(kList[k])->length() << " ";
    }    
    os << std::endl;
  }
  os.precision(oldPrecision);
}

void Tissue::printInitPly(std::ostream &os) {

  // Print header
  os << "ply" << std::endl
     << "format ascii 1.0" << std::endl
     << "comment Geometry (simplified) generated by tissue/organism" << std::endl
     << "comment http://www.thep.lu.se/~henrik/, henrik@thep.lu.se" << std::endl;
  // vertex element definition
  os << "element vertex " << numVertex() << std::endl
     << "property float x" << std::endl
     << "property float y" << std::endl;
  if (vertex(0).numPosition()==3) {
    os << "property float z" << std::endl;
  }
  os << "property uint index" << std::endl;
  // face element definition
  os << "element face " << numCell() << std::endl
     << "property list uchar uint vertex_index" << std::endl
     << "property uint index" << std::endl
     << "property list uchar float var" << std::endl;
  // edge element definition
  os << "element edge " << numWall() << std::endl
     << "property uint source" << std::endl
     << "property uint target" << std::endl
     << "property uint index" << std::endl
     << "property list uchar float var" << std::endl;
  os << "end_header" << std::endl;

  // Print vertex positions and index
  //os << "comment start of vertex list" << std::endl;
  for (size_t i=0; i<numVertex(); ++i) {
    for (size_t d=0; d<vertex(i).numPosition(); ++d) {
      os << vertex(i).position(d) << " ";
    }
    os << vertex(i).index() << std::endl;
  }
  // Print face information, vertex_list, index  and variables
  //os << "comment start of face list" << std::endl;
  for (size_t i=0; i<numCell(); ++i) {
    os << cell(i).numVertex() << " ";
    for (size_t k=0; k<cell(i).numVertex(); ++k) {
      os << cell(i).vertex(k)->index() << " ";
    }
    os << cell(i).index() << " ";
    os << cell(i).numVariable() << " ";
    for (size_t k=0; k<cell(i).numVariable(); ++k) {
      os << cell(i).variable(k) << " ";
    }
     os << std::endl;
  }
  // Print edge information, source, target, index, and variables
  //os << "comment start of edge list" << std::endl;
  for (size_t i=0; i<numWall(); ++i) {
    os << wall(i).vertex1()->index() << " " << wall(i).vertex2()->index() << " ";
    os << wall(i).index() << " ";
    os << wall(i).numVariable() << " ";
    for (size_t k=0; k<wall(i).numVariable(); ++k) {
      os << wall(i).variable(k) << " ";
    }
    os << std::endl;
  }
}

void Tissue::printVertex(std::ostream &os) {
  
  for( size_t i=0 ; i<numVertex() ; ++i )
    os << vertex(i).position(0) << " "
       << vertex(i).position(1) << std::endl;
}

void Tissue::printWall(std::ostream &os) {
  
  for( size_t i=0 ; i<numWall() ; ++i ) {
    for( size_t d=0 ; d<wall(i).vertex1()->numPosition() ; ++d )
      os << wall(i).vertex1()->position(d) << " "; 
    os << std::endl;
    for( size_t d=0 ; d<wall(i).vertex1()->numPosition() ; ++d )
      os << wall(i).vertex2()->position(d) << " "; 
    os << std::endl << std::endl << std::endl;
  }
}

void Tissue::printVertexAndCell(std::ostream &os) {
  
  size_t Nv = numVertex(); 
  if( !Nv ) {
    os << "0 0" << std::endl << "0 0" << std::endl;
    return;
  }
  //Print the vertex positions
  size_t dimension = vertex(0).numPosition();
  os << numVertex() << " " << dimension << std::endl;
  for( size_t i=0 ; i<Nv ; ++i ) {
    for( size_t d=0 ; d<dimension ; ++d )
      os << vertex(i).position(d) << " ";
    os << std::endl;
  }
  os << std::endl;
  //Print the cells, first connected vertecis and then variables
  size_t Nc = numCell();
  static std::vector<size_t> randomIndex(Nc);
  static size_t flag=1;
  if( flag ) {
    for( size_t i=0 ; i<Nc ; ++i ) randomIndex[i]=i;
    random_shuffle(randomIndex.begin(),randomIndex.end());
    flag=0;
  }
  os << Nc << " 3" << std::endl;
  for( size_t i=0 ; i<Nc ; ++i ) {
    size_t Ncv = cell(i).numVertex(); 
    os << Ncv << " ";
    for( size_t k=0 ; k<Ncv ; ++k )
      os << cell(i).vertex(k)->index() << " ";
    os << i << " " << randomIndex[i] << " " << cell(i).volume() 
       << std::endl;
  }
}

void Tissue::
printVertexAndCell(DataMatrix &cellData,
		   DataMatrix &vertexData,
		   std::ostream &os) {
  
  size_t Nv = vertexData.size(); 
  if( !Nv ) {
    os << "0 0" << std::endl << "0 0" << std::endl;
    return;
  }
  //Print the vertex positions
  size_t dimension = vertexData[0].size();
  os << numVertex() << " " << dimension << std::endl;
  for( size_t i=0 ; i<Nv ; ++i ) {
    for( size_t d=0 ; d<dimension ; ++d )
      os << vertexData[i][d] << " ";
    os << std::endl;
  }
  os << std::endl;
  //Print the cells, first connected vertecis and then variables
  size_t Nc = cellData.size();
  //For normal printing
  //
  int numPrintVar=cell(0).numVariable()+3;
  os << Nc << " " << numPrintVar << std::endl;
  for( size_t i=0 ; i<Nc ; ++i ) {
    size_t Ncv = cell(i).numVertex(); 
    os << Ncv << " ";
    for( size_t k=0 ; k<Ncv ; ++k )
      os << cell(i).vertex(k)->index() << " ";
    
    for( size_t k=0 ; k<cellData[i].size() ; ++k )
      os << cellData[i][k] << " ";
    os << i << " " << cell(i).calculateVolume(vertexData) << " " 
       << cell(i).numWall() << std::endl;
  }
  //
  //End, normal printing
  
  //For membrane-PIN printing (version1)
  //
  // 	os << Nc << std::endl;
  // 	size_t aI=0,pI=1;
  //   for( size_t i=0 ; i<Nc ; ++i ) {
  //     size_t Ncv = cell(i).numVertex(); 
  //     os << Ncv << " ";
  //     for( size_t k=0 ; k<Ncv ; ++k )
  //       os << cell(i).vertex(k)->index() << " ";
  
  // 		double p3=0.01;
  // 		double sum=p3;
  // 		for( size_t n=0 ; n<Ncv ; ++n ) {
  // 			if( !cell(i).isNeighbor(background()) ) { 
  // 				if( cell(i).wall(n)->cell1()->index()==i )
  // 					sum += cellData[ cell(i).wall(n)->cell2()->index() ][ aI ];
  // 				else
  // 					sum += cellData[ cell(i).wall(n)->cell1()->index() ][ aI ];
  // 			}
  // 		}
  
  // 		os << p3*cellData[i][pI]/sum << " ";		
		
  // 		for( size_t k=0 ; k<Ncv ; ++k ) {
  // 			if( !cell(i).isNeighbor(background()) ) { 
  // 				size_t neighIndex; 
  // 				if( cell(i).wall(k)->cell1()->index()==i )
  // 					neighIndex = cell(i).wall(k)->cell2()->index();				
  // 				else
  // 					neighIndex = cell(i).wall(k)->cell1()->index();				
  // 				os << cellData[i][pI] * cellData[neighIndex][aI] / sum << " ";
  // 			}
  // 			else
  // 				os << 0.0 << " ";
  // 		}
  // 		os << std::endl;
  //   }	
  
  //   //For membrane PIN1 printing (version3)
  //   //
  //  	os << Nc << std::endl;
  // 	size_t auxinI=1,pinI=2,auxI=3,pidI=4,xI=5;
  //  	std::vector<double> parameter(8);
  // 	parameter[0]=0.1;
  // 	parameter[1]=1.0;
  // 	parameter[2]=0.9;
  // 	parameter[3]=0.1;
  // 	parameter[4]=1.0;
  // 	parameter[5]=1.0;
  // 	parameter[6]=0.0001;
  // 	parameter[7]=4;
  
  // 	for( size_t i=0 ; i<Nc ; ++i ) {
  // 		size_t Ncv = cell(i).numVertex(); 
  // 		os << Ncv << " ";
  // 		for( size_t k=0 ; k<Ncv ; ++k )
  // 			os << cell(i).vertex(k)->index() << " ";
  
  // 		//Transport
  // 		//
  // 		size_t numWalls=cell(i).numWall();
  // 		//PID factor
  // 		double tmpPow = std::pow(cellData[i][pidI],parameter[7]);
  // 		double Ci = tmpPow/(tmpPow+std::pow(parameter[6],parameter[7]));
  
  // 		//Polarization coefficient normalization constant
  // 		double sum=0.0;
  // 		std::vector<double> Pij(numWalls);
  // 		for( size_t n=0 ; n<numWalls ; ++n ) {
  // 			if( cell(i).wall(n)->cell1() != background() &&
  // 					cell(i).wall(n)->cell2() != background() ) { 
  // 				size_t neighI;
  // 				if( cell(i).wall(n)->cell1()->index()==i )
  // 					neighI = cell(i).wall(n)->cell2()->index();
  // 				else
  // 					neighI = cell(i).wall(n)->cell1()->index();
  // 				//double powX = std::pow(cellData[ neighI ][ xI ],parameter[5));
  // 				//double Cij = powX/(std::pow(parameter[4),parameter[5))+powX);
  // 				double Cij = cellData[ neighI ][ xI ];
  // 				sum += Pij[n] = (1.0-parameter[2]) + 
  // 					parameter[2]*(Ci*Cij+(1.0-Ci)*(1.0-Cij));
  // 				//sum += Pij[n] = (1.0-parameter[2]) + 
  // 				//parameter[2]*cellData[ neighI ][xI];
  // 			}
  // 			else 
  // 				sum += Pij[n] = (1.0-parameter[2]);
  // 		}
  // 		//sum /= numWalls;//For adjusting for different num neigh
  // 		sum += parameter[3];
  
  // 		if( sum >= 0.0 )
  // 			std::cout << parameter[3]*cellData[i][pinI] / sum << " ";
  // 		else
  // 			std::cout << "0.0 ";
  
  // 		for( size_t n=0 ; n<numWalls ; ++n ) {
  // 			double pol=0.0;
  // 			if( sum != 0.0 )
  // 				pol = cellData[i][pinI] * Pij[n] / sum;
  // 			std::cout << pol << " ";
  // 		}
  // 		std::cout << std::endl;
  // 	}
  //
  //End Pij printing version3
  
  // 	//Transport
  // 	//
  // 	size_t numWalls=T.cell(i).numWall();
  // 	//PID factor
  // 	double tmpPow = std::pow(cellData[i][pidI],parameter(7));
  // 	double Ci = tmpPow/(tmpPow+std::pow(parameter(29),parameter(30)));
  //     //Polarization coefficient normalization constant
  //     double sum=0.0;
  //     std::vector<double> Pij(numWalls);
  //     for( size_t n=0 ; n<numWalls ; ++n ) {
  //       if( T.cell(i).wall(n)->cell1() != T.background() &&
  // 					T.cell(i).wall(n)->cell2() != T.background() ) { 
  // 				size_t neighI;
  // 				if( T.cell(i).wall(n)->cell1()->index()==i )
  // 					neighI = T.cell(i).wall(n)->cell2()->index();
  // 				else
  // 					neighI = T.cell(i).wall(n)->cell1()->index();
  // 				double powX = std::pow(cellData[ neighI ][ xI ],parameter(28));
  // 				double Cij = powX/(std::pow(parameter(27),parameter(28))+powX);
  // 				sum += Pij[n] = (1.0-parameter(25)) + 
  // 					parameter(25)*(Ci*Cij+(1.0-Ci)*(1.0-Ci));
  //       }
  //       else 
  // 				sum += Pij[n] = (1.0-parameter(25));
  //     }
  //     //sum /= numWalls;//For adjusting for different num neigh
  //     sum += parameter(26);
  
  //     for( size_t n=0 ; n<numWalls ; ++n ) {
  //       //if( !T.cell(i).isNeighbor(T.background()) ) { 
  //       if( T.cell(i).wall(n)->cell1() != T.background() &&
  // 					T.cell(i).wall(n)->cell2() != T.background() ) { 
  // 				size_t neighI; 
  // 				if( T.cell(i).wall(n)->cell1()->index()==i )
  // 					neighI = T.cell(i).wall(n)->cell2()->index();				
  // 				else
  // 					neighI = T.cell(i).wall(n)->cell1()->index();				
  // 				double pol=0.0;
  // 				if( sum != 0.0 )
  // 					pol = cellData[i][pinI] * Pij[n] / sum;
  // 				double transportRate = parameter(23)*cellData[neighI][auxI]*
  // 					pol*cellData[i][auxinI] /
  // 					( (parameter(24)+cellData[i][auxinI])*
  // 						(cellData[i][auxI]+cellData[neighI][auxI]) );
  // 				cellDerivs[i][auxinI] -= transportRate + 
  // 					parameter(31)*cellData[i][auxinI];
  // 				cellDerivs[neighI][auxinI] += transportRate +
  // 					parameter(31)*cellData[i][auxinI];
  //       }
  //     }
  
  //
  //End Pij printing version2
  
  
  //   for( size_t i=0 ; i<Nc ; ++i ) {
  //     size_t Ncv = cell(i).numVertex(); 
  //     os << Ncv << " ";
  //     for( size_t k=0 ; k<Ncv ; ++k )
  //       os << cell(i).vertex(k)->index() << " ";
  
  // 		double p3=0.01;
  // 		double sum=p3;
  // 		for( size_t n=0 ; n<Ncv ; ++n ) {
  // 			if( !cell(i).isNeighbor(background()) ) { 
  // 				if( cell(i).wall(n)->cell1()->index()==i )
  // 					sum += cellData[ cell(i).wall(n)->cell2()->index() ][ aI ];
  // 				else
  // 					sum += cellData[ cell(i).wall(n)->cell1()->index() ][ aI ];
  // 			}
  // 		}
  
  // 		os << p3*cellData[i][pI]/sum << " ";		
  
  // 		for( size_t k=0 ; k<Ncv ; ++k ) {
  // 			if( !cell(i).isNeighbor(background()) ) { 
  // 				size_t neighIndex; 
  // 				if( cell(i).wall(k)->cell1()->index()==i )
  // 					neighIndex = cell(i).wall(k)->cell2()->index();				
  // 				else
  // 					neighIndex = cell(i).wall(k)->cell1()->index();				
  // 				os << cellData[i][pI] * cellData[neighIndex][aI] / sum << " ";
  // 			}
  // 			else
  // 				os << 0.0 << " ";
  // 		}
  // 		os << std::endl;
  //   }	
  
}

void Tissue::
printVertexAndWall(DataMatrix &wallData,
		   DataMatrix &vertexData,
		   std::ostream &os) {
  
  size_t Nv = vertexData.size(); 
  if( !Nv ) {
    os << "0 0" << std::endl << "0 0" << std::endl;
    return;
  }
  //Print the vertex positions
  size_t dimension = vertexData[0].size();
  os << numVertex() << " " << dimension << std::endl;
  for( size_t i=0 ; i<Nv ; ++i ) {
    for( size_t d=0 ; d<dimension ; ++d )
      os << vertexData[i][d] << " ";
    os << std::endl;
  }
  os << std::endl;
  // Print the walls, first connected vertecis and then variables
  size_t Nw = wallData.size();
  //
  // Print wall variables
  //
  int numPrintVar=wallData[0].size()+2;
  os << Nw << " " << numPrintVar << std::endl;
  for( size_t i=0 ; i<Nw ; ++i ) {
    os << "2 ";
    os << wall(i).vertex1()->index() << " " 
       << wall(i).vertex2()->index() << " ";
    for( size_t k=0 ; k<wallData[i].size() ; ++k )
      os << wallData[i][k] << " ";
    os << i << " " << wall(i).lengthFromVertexPosition(vertexData)
       << std::endl;
  }
}	
