/**
 * Filename     : adhocReaction.cc
 * Description  : Classes describing some ad hoc updates
 * Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
 * Created      : September 2007
 * Revision     : $Id:$
 */

#include"adhocReaction.h"
#include"tissue.h"
#include"baseReaction.h"
#include<cmath>

VertexNoUpdateFromPosition::
VertexNoUpdateFromPosition(std::vector<double> &paraValue, 
			   std::vector< std::vector<size_t> > 
			   &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //////////////////////////////////////////////////////////////////////
  if( paraValue.size()!=2 ) {
    std::cerr << "VertexNoUpdateFromPosition::"
	      << "VertexNoUpdateFromPosition() "
	      << "Uses two parameters threshold and direction "
	      << "(-1 -> less than)" 
	      << std::endl;
    exit(0);
  }
  if( paraValue[1] != 1.0 && paraValue[1] != -1.0 ) {
    std::cerr << "VertexNoUpdateFromPosition::"
	      << "VertexNoUpdateFromPosition() "
	      << "direction (second parameter) need to be 1 (greater than) "
	      << "or -1 (less than)." << std::endl;
    exit(0);
  }
  if( indValue.size() != 1 || indValue[0].size() != 1 ) {
    std::cerr << "VertexNoUpdateFromPosition::"
	      << "VertexNoUpdateFromPosition() "
	      << "Pos index in first level." << std::endl;
    exit(0);
  }
  //Set the variable values
  //////////////////////////////////////////////////////////////////////
  setId("VertexNoUpdateFromPosition");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //////////////////////////////////////////////////////////////////////
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "threshold";
  tmp[1] = "direction";
  setParameterId( tmp );
}

void VertexNoUpdateFromPosition::
derivs(Tissue &T,
       std::vector< std::vector<double> > &cellData,
       std::vector< std::vector<double> > &wallData,
       std::vector< std::vector<double> > &vertexData,
       std::vector< std::vector<double> > &cellDerivs,
       std::vector< std::vector<double> > &wallDerivs,
       std::vector< std::vector<double> > &vertexDerivs ) {
  
  //Check the cancelation for every vertex
  size_t numVertices = T.numVertex();
  size_t posIndex = variableIndex(0,0);
  size_t dimension = vertexData[0].size();
  assert( posIndex < vertexData[0].size() );
  
  for( size_t i=0 ; i<numVertices ; ++i )
    if( (parameter(1)>0 && vertexData[i][posIndex]>parameter(0)) || 
	(parameter(1)<0 && vertexData[i][posIndex]<parameter(0)) )
      for( size_t d=0 ; d<dimension ; ++d )
	vertexDerivs[i][d] = 0.0;
}

VertexNoUpdateBoundary::
VertexNoUpdateBoundary(std::vector<double> &paraValue, 
		       std::vector< std::vector<size_t> > 
		       &indValue ) {
  
  //Do some checks on the parameters and variable indices
  //
  if (paraValue.size()) {
    std::cerr << "VertexNoUpdateBoundary::"
	      << "VertexNoUpdateBoundary() "
	      << "Uses no parameters."
	      << std::endl;
    exit(0);
  }
  if (indValue.size()!=0 && indValue.size()!=1) {
    std::cerr << "VertexNoUpdateBoundary::"
	      << "VertexNoUpdateBoundary() "
	      << "Either no variable indices is used (all directions fixed), " 
	      << "or fixed directions are given in first level." << std::endl;
    exit(0);
  }
  //Set the variable values
  //
  setId("VertexNoUpdateBoundary");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  setParameterId( tmp );
}

void VertexNoUpdateBoundary::
derivs(Tissue &T,
       std::vector< std::vector<double> > &cellData,
       std::vector< std::vector<double> > &wallData,
       std::vector< std::vector<double> > &vertexData,
       std::vector< std::vector<double> > &cellDerivs,
       std::vector< std::vector<double> > &wallDerivs,
       std::vector< std::vector<double> > &vertexDerivs ) 
{  
  //Check the cancelation for every vertex
  size_t numVertices = T.numVertex();
  size_t dimension = vertexData[0].size();
  if (numVariableIndexLevel())
    for (size_t dI=0; dI<numVariableIndex(0); ++dI)
      assert(variableIndex(0,dI)<dimension);

  for (size_t i=0; i<numVertices; ++i)
    if (T.vertex(i).isBoundary(T.background())) {
      if (numVariableIndexLevel()) {
	for (size_t dI=0; dI<numVariableIndex(0); ++dI)
	  vertexDerivs[i][variableIndex(0,dI)] = 0.0;
      }
      else {
	for (size_t d=0; d<dimension; ++d)
	  vertexDerivs[i][d] = 0.0;
      }
    }
}

VertexTranslateToMax::
VertexTranslateToMax(std::vector<double> &paraValue, 
		     std::vector< std::vector<size_t> > 
		     &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //////////////////////////////////////////////////////////////////////
  if( paraValue.size()!=1 ) {
    std::cerr << "VertexTranslateToMax::VertexTranslateToMax() "
	      << "Uses one parameter, maxPos " << std::endl;
    exit(0);
  }
  if( indValue.size() != 1 || indValue[0].size() != 1 ) {
    std::cerr << "VertexTranslateToMax::"
	      << "VertexTranslateToMax() "
	      << "Pos index in first level." << std::endl;
    exit(0);
  }
  //Set the variable values
  //////////////////////////////////////////////////////////////////////
  setId("VertexTranslateToMax");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //////////////////////////////////////////////////////////////////////
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "maxPos";
  setParameterId( tmp );
}

void VertexTranslateToMax::
initiate(Tissue &T,
	 std::vector< std::vector<double> > &cellData,
	 std::vector< std::vector<double> > &wallData,
	 std::vector< std::vector<double> > &vertexData)
{
  size_t numVertices = T.numVertex();
  size_t posIndex = variableIndex(0,0);
  assert( posIndex < vertexData[0].size() );
  
  double max = vertexData[0][posIndex];
  for( size_t i=1 ; i<numVertices ; ++i )
    if (vertexData[i][posIndex]>max)
      max = vertexData[i][posIndex];
  double delta = max-parameter(0);
  for( size_t i=0 ; i<numVertices ; ++i )
    vertexData[i][posIndex] -= delta;
}

void VertexTranslateToMax::
derivs(Tissue &T,
       std::vector< std::vector<double> > &cellData,
       std::vector< std::vector<double> > &wallData,
       std::vector< std::vector<double> > &vertexData,
       std::vector< std::vector<double> > &cellDerivs,
       std::vector< std::vector<double> > &wallDerivs,
       std::vector< std::vector<double> > &vertexDerivs ) 
{
}

void VertexTranslateToMax::
update(Tissue &T,
       std::vector< std::vector<double> > &cellData,
       std::vector< std::vector<double> > &wallData,
       std::vector< std::vector<double> > &vertexData,
			 double h)
{
  size_t numVertices = T.numVertex();
  size_t posIndex = variableIndex(0,0);
  assert( posIndex < vertexData[0].size() );
  
	double max = vertexData[0][posIndex];
  for( size_t i=1 ; i<numVertices ; ++i )
    if (vertexData[i][posIndex]>max)
			max = vertexData[i][posIndex];
	double delta = max-parameter(0);
  for( size_t i=0 ; i<numVertices ; ++i )
		vertexData[i][posIndex] -= delta;
}

CenterCOM::CenterCOM(std::vector<double> &paraValue, std::vector< std::vector<size_t> > &indValue) {
	//Do some checks on the parameters and variable indeces
	//////////////////////////////////////////////////////////////////////
	if (paraValue.size() != 0) {
		std::cerr << "CenterCOM::CenterCOM() Uses no parameters.\n";
		std::exit(EXIT_FAILURE);
	}

	if (indValue.size() != 0) {
		std::cerr << "CenterCOM::CenterCOM() No variable indices used.\n";
		std::exit(EXIT_FAILURE);
	}

	//Set the variable values
	//////////////////////////////////////////////////////////////////////
	setId("CenterCOM");
	setParameter(paraValue);  
	setVariableIndex(indValue);
	
	//Set the parameter identities
	//////////////////////////////////////////////////////////////////////
	std::vector<std::string> tmp(numParameter());
	setParameterId(tmp);
}

void CenterCOM::initiate(Tissue &T,
	std::vector< std::vector<double> > &cellData,
	std::vector< std::vector<double> > &wallData,
	std::vector< std::vector<double> > &vertexData)
{
	update(T, cellData, wallData, vertexData, 0.0);
}

void CenterCOM::derivs(Tissue &T,
	std::vector< std::vector<double> > &cellData,
	std::vector< std::vector<double> > &wallData,
	std::vector< std::vector<double> > &vertexData,
	std::vector< std::vector<double> > &cellDerivs,
	std::vector< std::vector<double> > &wallDerivs,
	std::vector< std::vector<double> > &vertexDerivs) 
{

}

void CenterCOM::update(Tissue &T,
	std::vector< std::vector<double> > &cellData,
	std::vector< std::vector<double> > &wallData,
	std::vector< std::vector<double> > &vertexData,
	double h)
{
	size_t dimension = vertexData[0].size();
	size_t numVertices = T.numVertex();
	
	std::vector<double> com(dimension, 0.0);

	for (size_t i = 0; i < numVertices; ++i) {
		for (size_t d = 0; d < dimension; ++d) {
			com[d] += vertexData[i][d];
		}
	}
	
	for (size_t d = 0; d < dimension; ++d) {
		com[d] /= numVertices;
	}
	
	for (size_t i = 0; i < numVertices; ++i) {
		for (size_t d = 0; d < dimension; ++d) {
			vertexData[i][d] -= com[d];
		}
	}
}

CalculatePCAPlane::
CalculatePCAPlane(std::vector<double> &paraValue, 
									std::vector< std::vector<size_t> > 
						 &indValue ) {
  
	//
  // Do some checks on the parameters and variable indeces
  //
  if (paraValue.size()!=1) {
    std::cerr << "CalculatePCAPlane::"
							<< "CalculatePCAPlane() "
							<< "Uses one parameter: onlyInUpdateFlag "
							<< std::endl;
    exit(0);
  }
  if (indValue.size()!=0) {
    std::cerr << "CalculatePCAPlane::"
							<< "CalculatePCAPlane() "
							<< "No indices used." << std::endl;
    exit(0);
  }
  //Set the variable values
  //////////////////////////////////////////////////////////////////////
  setId("CalculatePCAPlane");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //////////////////////////////////////////////////////////////////////
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "onlyInUpdateFlag";
  setParameterId( tmp );
}

void CalculatePCAPlane::
initiate(Tissue &T,
				 std::vector< std::vector<double> > &cellData,
				 std::vector< std::vector<double> > &wallData,
				 std::vector< std::vector<double> > &vertexData)
{
	if (parameter(0)==1.0) {
		size_t numCell = T.numCell();
		for (size_t i=0; i<numCell; ++i)
			T.cell(i).calculatePCAPlane(vertexData);
	}
}

void CalculatePCAPlane::
derivs(Tissue &T,
       std::vector< std::vector<double> > &cellData,
       std::vector< std::vector<double> > &wallData,
       std::vector< std::vector<double> > &vertexData,
       std::vector< std::vector<double> > &cellDerivs,
       std::vector< std::vector<double> > &wallDerivs,
       std::vector< std::vector<double> > &vertexDerivs ) 
{
	if (parameter(0)!=1.0) {
		size_t numCell = T.numCell();
		for (size_t i=0; i<numCell; ++i)
			T.cell(i).calculatePCAPlane(vertexData);
	}
}

void CalculatePCAPlane::
update(Tissue &T,
       std::vector< std::vector<double> > &cellData,
       std::vector< std::vector<double> > &wallData,
       std::vector< std::vector<double> > &vertexData,
			 double h)
{
	if (parameter(0)==1.0) {
		size_t numCell = T.numCell();
		for (size_t i=0; i<numCell; ++i)
			T.cell(i).calculatePCAPlane(vertexData);
	}
}

InitiateWallLength::
InitiateWallLength(std::vector<double> &paraValue, 
									 std::vector< std::vector<size_t> > 
									 &indValue ) 
{
  //
  //Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=1 ) {
    std::cerr << "InitiateWallLength::InitiateWallLength() "
							<< "Uses one parameter, factor. " << std::endl;
    exit(0);
  }
  if( indValue.size() != 0 ) {
    std::cerr << "InitiateWallLength::"
							<< "InitiateWallLength() "
							<< "No variable indices used." << std::endl;
    exit(0);
  }
	//
  //Set the variable values
  //
  setId("InitiateWallLength");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  //
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
	tmp[0] = "factor";
  setParameterId( tmp );
}

void InitiateWallLength::
initiate(Tissue &T,
				 std::vector< std::vector<double> > &cellData,
				 std::vector< std::vector<double> > &wallData,
				 std::vector< std::vector<double> > &vertexData)
{
	size_t dimension = vertexData[0].size();
  size_t numWall = T.numWall();
  
  for (size_t i=0; i<numWall; ++i) {
		double distance=0.0;
		size_t v1I=T.wall(i).vertex1()->index();
		size_t v2I=T.wall(i).vertex2()->index();
		for (size_t d=0; d<dimension; ++d )
			distance += (vertexData[v2I][d]-vertexData[v1I][d])*(vertexData[v2I][d]-vertexData[v1I][d]);
		distance = std::sqrt(distance);
		wallData[i][0] = parameter(0)*distance;
	}
}

void InitiateWallLength::
derivs(Tissue &T,
       std::vector< std::vector<double> > &cellData,
       std::vector< std::vector<double> > &wallData,
       std::vector< std::vector<double> > &vertexData,
       std::vector< std::vector<double> > &cellDerivs,
       std::vector< std::vector<double> > &wallDerivs,
       std::vector< std::vector<double> > &vertexDerivs ) 
{
}

InitiateWallMesh::
InitiateWallMesh(std::vector<double> &paraValue, 
		 std::vector< std::vector<size_t> > 
		 &indValue ) 
{
  //
  // Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=1 ) {
    std::cerr << "InitiateWallMesh::InitiateWallMesh() "
	      << "Uses one parameter, #of added vertices per wall." << std::endl;
    exit(0);
  }
  if( indValue.size() != 0 ) {
    std::cerr << "InitiateWallMesh::"
	      << "InitiateWallMesh() "
	      << "No variable indices used." << std::endl;
    exit(0);
  }
  //
  //Set the variable values
  //
  setId("InitiateWallMesh");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  //
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "numAddedVerticesPerWall";
  setParameterId( tmp );
}

void InitiateWallMesh::
initiate(Tissue &T,
	 std::vector< std::vector<double> > &cellData,
	 std::vector< std::vector<double> > &wallData,
	 std::vector< std::vector<double> > &vertexData)
{
  size_t dimension = vertexData[0].size();
  size_t numWallOld = T.numWall();
  size_t numWall = numWallOld;
  size_t numVertex = T.numVertex();
  size_t numMesh = (size_t) parameter(0);
  if (!numMesh) 
    return;
  double scaleFactor = 1.0/(1.0+numMesh); 
  
  for (size_t i=0; i<numWallOld; ++i) {		
    // Add the new walls and vertices, and set indices, 
    // and wall lengths and vertex positions (into the data matrices)
    wallData[i][0] *= scaleFactor;
    size_t v1I = T.wall(i).vertex1()->index();
    size_t v2I = T.wall(i).vertex2()->index();
    std::vector<double> n(dimension);
    for (size_t d=0; d<dimension; ++d)
      n[d] = scaleFactor*(vertexData[v2I][d]-vertexData[v1I][d]);
    for (size_t k=0; k<numMesh; ++k) {
      T.addWall(T.wall(i));
      T.wall(numWall+k).setIndex(numWall+k);
      wallData.push_back(wallData[i]);//with the new correct length
      T.addVertex(T.vertex(v1I));
      T.vertex(numVertex+k).setIndex(numVertex+k);
      vertexData.push_back(vertexData[v1I]);
      for (size_t d=0; d<dimension; ++d)
	vertexData[numVertex+k][d] += (k+1)*n[d];
    }
    // Adjust cell connections
    size_t c1I = T.wall(i).cell1()->index();
    size_t c2I = T.wall(i).cell2()->index();
    for (size_t k=0; k<numMesh; ++k) {
      T.cell(c1I).addWall(T.wallP(numWall+k));
      T.cell(c2I).addWall(T.wallP(numWall+k));
      T.cell(c1I).addVertex(T.vertexP(numVertex+k));
      T.cell(c2I).addVertex(T.vertexP(numVertex+k));
    }
    // Adjust wall connections
    T.wall(i).setVertex2(T.vertexP(numVertex));
    for (size_t k=0; k<numMesh-1; ++k) {
      T.wall(numWall+k).setVertex1(T.vertexP(numVertex+k));
      T.wall(numWall+k).setVertex2(T.vertexP(numVertex+k+1));
    }
    T.wall(numWall+numMesh-1).setVertex1(T.vertexP(numVertex+numMesh-1));
    // Adjust vertex connections
    std::vector<Cell*> tmpCell;
    if (T.wall(i).cell1()!=T.background())
      tmpCell.push_back(T.wall(i).cell1());
    if (T.wall(i).cell2()!=T.background())
      tmpCell.push_back(T.wall(i).cell2());
    std::vector<Wall*> tmpWall(2);
    tmpWall[1] = T.wallP(i);
    for (size_t k=0; k<numMesh; ++k) {
      T.vertex(numVertex+k).setCell(tmpCell);
      tmpWall[0] = tmpWall[1];
      tmpWall[1] = T.wallP(numWall+k);
      T.vertex(numVertex+k).setWall(tmpWall);
    }
    size_t numVertexWall = T.vertex(v2I).numWall();
    size_t updatedFlag=0;
    for (size_t k=0; k<numVertexWall; ++k) {
      if (T.vertex(v2I).wall(k)==T.wallP(i)) {
	T.vertex(v2I).setWall(k,T.wallP(numWall+numMesh-1));
	++updatedFlag;
      }
    }
    assert(updatedFlag==1);
    numWall += numMesh;
    numVertex += numMesh;
  }
  // Sort the walls and vertices in the cells again since things have been added
  T.sortCellWallAndCellVertex();
  T.checkConnectivity(1);
}

void InitiateWallMesh::
derivs(Tissue &T,
       std::vector< std::vector<double> > &cellData,
       std::vector< std::vector<double> > &wallData,
       std::vector< std::vector<double> > &vertexData,
       std::vector< std::vector<double> > &cellDerivs,
       std::vector< std::vector<double> > &wallDerivs,
       std::vector< std::vector<double> > &vertexDerivs ) 
{
}

StrainTest::
StrainTest(std::vector<double> &paraValue, 
								 std::vector< std::vector<size_t> > 
								 &indValue ) 
{
  //
  //Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()<2 ) {
    std::cerr << "StrainTest::StrainTest() "
							<< "Uses at least two parameters, version, angle and optional parameters." << std::endl;
    exit(0);
  }
  if( indValue.size() != 0 ) {
    std::cerr << "StrainTest::"
							<< "StrainTest() "
							<< "No variable indices used." << std::endl;
    exit(0);
  }
	//
  //Set the variable values
  //
  setId("StrainTest");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  //
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
	tmp[0] = "version";
	tmp[1] = "angle";
  setParameterId( tmp );
}

void StrainTest::
initiate(Tissue &T,
				 std::vector< std::vector<double> > &cellData,
				 std::vector< std::vector<double> > &wallData,
				 std::vector< std::vector<double> > &vertexData)
{
}

void StrainTest::
derivs(Tissue &T,
       std::vector< std::vector<double> > &cellData,
       std::vector< std::vector<double> > &wallData,
       std::vector< std::vector<double> > &vertexData,
       std::vector< std::vector<double> > &cellDerivs,
       std::vector< std::vector<double> > &wallDerivs,
       std::vector< std::vector<double> > &vertexDerivs ) 
{
	size_t dimension=vertexData[0].size();
	double angleRad = 2*3.14159*parameter(1)/360;
	assert(dimension==2);
	std::vector<double> n(dimension);
	n[0] = std::cos(angleRad);
	n[1] = std::sin(angleRad);
	if (parameter(0)==0.0) {
		if (numParameter()!=2) {
			std::cerr << "StrainTest::derivs() version " << parameter(0) << " requires " << 2
								<< " parameters." << std::endl;
			exit(-1);
		}
		for (size_t d=0; d<dimension; ++d) {
			double sign = vertexData[3][d]>0.0 ? 1.0 : -1.0;
			vertexDerivs[3][d] = n[d]*sign;
		}
	}
	else if (parameter(0)==1.0) {
		if (numParameter()!=2) {
			std::cerr << "StrainTest::derivs() version " << parameter(0) << " requires " << 2
								<< " parameters." << std::endl;
			exit(-1);
		}
		for (size_t vI=0; vI<vertexData.size(); ++vI) {
			for (size_t d=0; d<dimension; ++d) {
				double sign = vertexData[vI][d]>0.0 ? 1.0 : -1.0;
				vertexDerivs[vI][d] = n[d]*sign;
			}
		}
	}	
}

CalculateVertexStressDirection::CalculateVertexStressDirection(std::vector<double> &paraValue, 
	std::vector< std::vector<size_t> > &indValue)
{
	if (paraValue.size() != 1) {
		std::cerr << "CalculateVertexStressDirection::CalculateVertexStressDirection() "
		<< "Uses one parameter: onlyInUpdateFlag\n";
		exit(0);
	}
	
	if (indValue.size() != 1) {
		std::cerr << "CalculateVertexStressDirection::CalculateVertexStressDirection() "
		<< "Wall force indexes.\n";
		exit(0);
	}
	
	setId("CalculateVertexStressDirection");
	setParameter(paraValue);  
	setVariableIndex(indValue);
	
	std::vector<std::string> tmp(numParameter());
	tmp[0] = "onlyInUpdateFlag";
	setParameterId(tmp);

	for (size_t i = 0; i < indValue[0].size(); ++i) {
		wallForceIndexes_.push_back(indValue[0][i]);
	}
}

void CalculateVertexStressDirection::initiate(Tissue &T,
	std::vector< std::vector<double> > &cellData,
	std::vector< std::vector<double> > &wallData,
	std::vector< std::vector<double> > &vertexData)
{
	if (parameter(0) == 1.0) {
		size_t numVertex = T.numVertex();
		
		for (size_t i = 0; i < numVertex; ++i) {
			T.vertex(i).calculateStressDirection(vertexData, wallData, wallForceIndexes_);
		}
	}
}

void CalculateVertexStressDirection::derivs(Tissue &T,
	std::vector< std::vector<double> > &cellData,
	std::vector< std::vector<double> > &wallData,
	std::vector< std::vector<double> > &vertexData,
	std::vector< std::vector<double> > &cellDerivs,
	std::vector< std::vector<double> > &wallDerivs,
	std::vector< std::vector<double> > &vertexDerivs) 
{
	if (parameter(0) != 1.0) {
		size_t numVertex = T.numVertex();
		for (size_t i=0; i<numVertex; ++i) {
			T.vertex(i).calculateStressDirection(vertexData, wallData, wallForceIndexes_);
		}
	}
}

void CalculateVertexStressDirection::update(Tissue &T,
	std::vector< std::vector<double> > &cellData,
	std::vector< std::vector<double> > &wallData,
	std::vector< std::vector<double> > &vertexData,
	double h)
{
	if (parameter(0) == 1.0) {
		size_t numVertex = T.numVertex();
		for (size_t i=0; i<numVertex; ++i) {
			T.vertex(i).calculateStressDirection(vertexData, wallData, wallForceIndexes_);
		}
	}
}
