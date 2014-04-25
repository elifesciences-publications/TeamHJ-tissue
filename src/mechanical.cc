//
// Filename     : mechanical.cc
// Description  : Classes describing updates due to mechanical interaction
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : April 2006
// Revision     : $Id:$
//
#include <utility>
#include <vector>
#include "baseReaction.h"
#include "mechanical.h"
#include "tissue.h"

VertexFromCellPressure::
VertexFromCellPressure(std::vector<double> &paraValue, 
		       std::vector< std::vector<size_t> > 
		       &indValue ) 
{  
  //Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=2 || (paraValue[1] != 0.0 && paraValue[1] != 1.0) ) {
    std::cerr << "VertexFromCellPressure::"
	      << "VertexFromCellPressure() "
	      << "Uses two parameters K_force and normalizeVolumeFlag (= 0 or 1).\n";
    exit(0);
  }
  if( indValue.size() != 0 ) {
    std::cerr << "VertexFromCellPressure::"
	      << "VertexFromCellPressure() "
	      << "No index given.\n";
    exit(0);
  }
  //Set the variable values
  //
  setId("VertexFromCellPressure");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "K_force";
  tmp[1] = "f_V_norm";
  setParameterId( tmp );
}

void VertexFromCellPressure::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) {
  
  //Do the update for each vertex via each wall in each cell
  size_t numCells = T.numCell();
  size_t dimension;
  dimension = T.vertex(0).numPosition(); 
  
  //Assumming vertices and walls are sorted and 2 dimensions
  //
  assert( dimension==2 );
  //For each cell
  for (size_t cellI = 0; cellI < numCells; ++cellI) {
    Cell &tmpCell = T.cell(cellI);
    
    double factor = 0.5 * parameter(0);
    
    if (parameter(1) == 1)
      {
	double cellVolume = tmpCell.calculateVolume(vertexData);
	factor /= std::fabs(cellVolume);
      }
    
    for (size_t k = 0; k < tmpCell.numVertex(); ++k) {
      size_t v1I = tmpCell.vertex(k)->index();
      size_t v1PlusI = tmpCell.vertex((k + 1) % (tmpCell.numVertex()))->index();
      size_t v1MinusK = k > 0 ? k - 1 : tmpCell.numVertex() - 1;
      size_t v1MinusI = tmpCell.vertex(v1MinusK)->index();
      
      vertexDerivs[v1I][0] += factor * (vertexData[v1PlusI][1] - vertexData[v1MinusI][1]);
      vertexDerivs[v1I][1] += factor * (vertexData[v1MinusI][0]- vertexData[v1PlusI][0]);
    }
  }
}

namespace CenterTriangulation {
  VertexFromCellPressure::
  VertexFromCellPressure(std::vector<double> &paraValue, 
			 std::vector< std::vector<size_t> > 
			 &indValue ) 
  {  
    //Do some checks on the parameters and variable indeces
  //
    if( paraValue.size()!=2 || (paraValue[1] != 0.0 && paraValue[1] != 1.0) ) {
      std::cerr << "CenterTriangulation::VertexFRomCellPressure::"
		<< "VertexFromCellPressure() " << std::endl
		<< "Uses two parameters K_force and normalizeVolumeFlag (= 0 or 1).\n";
      exit(EXIT_FAILURE);
    }
    if( indValue.size() != 1 || indValue[0].size() != 2 ) {
      std::cerr << "CenterTriangulation::VertexFRomCellPressure::"
		<< "VertexFromCellPressure() " << std::endl
		<< "Start of additional Cell variable indices (center(x,y,z) "
		<< "L_1,...,L_n, n=num vertex) is given in first level." 
		<< "concentration index is given in second level. "
		<< std::endl;    
      exit(EXIT_FAILURE);
    }
    //Set the variable values
    //
    setId("CenterTriangulation::VertexFromCellPressure");
    setParameter(paraValue);  
    setVariableIndex(indValue);
    
    //Set the parameter identities
    //
    std::vector<std::string> tmp( numParameter() );
    tmp[0] = "K_force";
    tmp[1] = "f_V_norm";
    setParameterId( tmp );
  }
  
  void VertexFromCellPressure::
  derivs(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs ) {
  
    //Do the update for each vertex via each wall in each cell
    size_t numCells = T.numCell();
    size_t dimension = T.vertex(0).numPosition(); 
    std::vector<double> cellCenter(dimension);
    
  //Assumming vertices and walls are sorted
  //
  //For each cell
    for (size_t cellI = 0; cellI < numCells; ++cellI) {
      Cell &tmpCell = T.cell(cellI);
      
      for (size_t d=0; d<dimension; ++d) {
	cellCenter[d] = cellData[cellI][variableIndex(0,0)+d];
      }
      
      for (size_t k = 0; k < tmpCell.numVertex(); ++k) {
	size_t v1I = tmpCell.vertex(k)->index();
	size_t v2I = tmpCell.vertex((k + 1) % (tmpCell.numVertex()))->index();

	// std::vector<double> nCell(dimension),nWall(dimension);
	// //normal to the triangle plane
	// nCell[0]=(vertexData[v1I][1]-cellCenter[1])*(vertexData[v2I][2]-cellCenter[2])-(vertexData[v1I][2]-cellCenter[2])*(vertexData[v2I][1]-cellCenter[1]);
	// nCell[1]=(vertexData[v1I][2]-cellCenter[2])*(vertexData[v2I][0]-cellCenter[0])-(vertexData[v1I][0]-cellCenter[0])*(vertexData[v2I][2]-cellCenter[2]);
	// nCell[2]=(vertexData[v1I][0]-cellCenter[0])*(vertexData[v2I][1]-cellCenter[1])-(vertexData[v1I][1]-cellCenter[1])*(vertexData[v2I][0]-cellCenter[0]);
	// //normal to the wall outward
	// nWall[0]=(vertexData[v2I][1]-vertexData[v1I][1])*nCells[2]-(vertexData[v2I][2]-vertexData[v1I][2])*nCells[1];
	// nWall[1]=(vertexData[v2I][2]-vertexData[v1I][2])*nCells[0]-(vertexData[v2I][0]-vertexData[v1I][0])*nCells[2];
	// nWall[2]=(vertexData[v2I][0]-vertexData[v1I][0])*nCells[1]-(vertexData[v2I][1]-vertexData[v1I][1])*nCells[0];

	// double dxNorm = 0.0;
        // dxNorm =nWall[0]*nWall[0]+nWall[1]*nWall[1]+nWall[2]*nWall[2]; 
	// if (dxNorm>0.0) {
	//   dxNorm = std::sqrt(dxNorm);
	//   for( size_t d=0 ; d<dimension ; ++d ) {
	//     nWall[d] /=dxNorm;
	//   }
	// }
	// else {
	//   std::cerr << "CellTriangulationVertexFromCellPressure::derivs() "
	// 	    << "strange wall length or direction." << std::endl;
	// }




	std::vector<double> x0(dimension),dx(dimension);
	double dxNorm = 0.0;
	for( size_t d=0 ; d<dimension ; ++d ) {
	  x0[d] = 0.5*(vertexData[v1I][d] + vertexData[v2I][d]);
	  // Caveat: This is NOT  normal to the wall in the triangle plane
	  dx[d] = x0[d]-cellCenter[d];
	}




	double wallLength = 0.0;
	double wallFactor = 0.0; // dx.wallVector/wallVector^2 for making dx perpendicular to the wall 
	std::vector<double> wallVector(dimension);
	
	for( size_t d=0 ; d<dimension ; ++d ) {
	  wallVector[d]=vertexData[v1I][d] - vertexData[v2I][d];
	  wallLength +=wallVector[d]*wallVector[d];
	}
	if (wallLength>0.0) {
	  wallLength = std::sqrt(wallLength);
	}
	else {
	  std::cerr << "CenterTriangulation::VertexFromCellPressure::derivs() "
		    << "Strange wall length." << std::endl;
	}
	wallFactor= ( dx[0]*wallVector[0]
		      +dx[1]*wallVector[1]
		      +dx[2]*wallVector[2])/(wallLength*wallLength);
	for( size_t d=0 ; d<dimension ; ++d ) {
	  
	  dx[d]=dx[d]-wallFactor*wallVector[d];
	  dxNorm += dx[d]*dx[d];
	}
	if (dxNorm>0.0) {
	  dxNorm = std::sqrt(dxNorm);
	  double dxInv = 1.0/dxNorm;
	  for( size_t d=0 ; d<dimension ; ++d ) {
	    dx[d] = dx[d]*dxInv;
	  }
	}
	else {
	  std::cerr << "CellTriangulationVertexFromCellPressure::derivs() "
		    << "Force direction undetermined." << std::endl;
	}
	double factor = 0.5 * parameter(0);
	if (parameter(1) == 1)
	  {
	    //NOTE maybe this one should be calculated using the central mesh vertex?
	    double cellVolume = tmpCell.calculateVolume(vertexData);                                  
	    factor /= std::fabs(cellVolume);         
	  }
	if(variableIndex(0,1)!=0) {
	  factor*=cellData[cellI][variableIndex(0,1)];
	}

	factor *= wallLength;    
	for( size_t d=0 ; d<dimension ; ++d ) {
	  vertexDerivs[v1I][d] += factor * dx[d];
	  vertexDerivs[v2I][d] += factor * dx[d];	
	}
      }
    }
  }
  
  VertexFromCellPressureLinear::
  VertexFromCellPressureLinear(std::vector<double> &paraValue, 
			       std::vector< std::vector<size_t> > 
			       &indValue ) 
  {  
    //Do some checks on the parameters and variable indeces
    //
    if( paraValue.size()!=3 || (paraValue[1] != 0.0 && paraValue[1] != 1.0) ) {
      std::cerr << "CenterTriangulation::VertexFromCellPressureLinear::"
		<< "VertexFromCellPressureLinear() " << std::endl
		<< "Uses three parameters K_force and normalizeVolumeFlag (= 0 or 1)" 
		<< "and deltaT that sets the time the linear increase is applied.\n";
      exit(0);
    }
    if( indValue.size() != 1 || indValue[0].size() != 1 ) {
      std::cerr << "CenterTriangulation::VertexFromCellPressureLinear::"
		<< "VertexFromCellPressureLinear() " << std::endl
		<< "Start of additional Cell variable indices (center(x,y,z) "
		<< "L_1,...,L_n, n=num vertex) is given in first level." 
		<< std::endl;    
      exit(0);
    }
    //Set the variable values
    //
    setId("CenterTriangulation::VertexFromCellPressureLinear");
    setParameter(paraValue);  
    setVariableIndex(indValue);
    
    //Set the parameter identities
    //
    std::vector<std::string> tmp( numParameter() );
    tmp[0] = "K_force";
    tmp[1] = "f_V_norm";
    tmp[1] = "deltaT";
    timeFactor_=0.0;
    setParameterId( tmp );
  }
  
  void VertexFromCellPressureLinear::
  derivs(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs ) {
    
    //Do the update for each vertex via each wall in each cell
    size_t numCells = T.numCell();
    size_t dimension = T.vertex(0).numPosition(); 
    std::vector<double> cellCenter(dimension);
    
    //Assumming vertices and walls are sorted
    //
    //For each cell
    for (size_t cellI = 0; cellI < numCells; ++cellI) {
      Cell &tmpCell = T.cell(cellI);
      
      for (size_t d=0; d<dimension; ++d) {
	cellCenter[d] = cellData[cellI][variableIndex(0,0)+d];
      }
      
      for (size_t k = 0; k < tmpCell.numVertex(); ++k) {
	size_t v1I = tmpCell.vertex(k)->index();
	size_t v2I = tmpCell.vertex((k + 1) % (tmpCell.numVertex()))->index();
	std::vector<double> x0(dimension),dx(dimension);
	double dxNorm = 0.0;
	for( size_t d=0 ; d<dimension ; ++d ) {
	  x0[d] = 0.5*(vertexData[v1I][d] + vertexData[v2I][d]);
	  // Caveat: This is  normal to the wall in the triangle plane
	  dx[d] = x0[d]-cellCenter[d];
	}
	double wallLength = 0.0;
	double wallFactor = 0.0; // dx.wallVector/wallVector^2 for making dx perpendicular to the wall 
	std::vector<double> wallVector(dimension);
	
	for( size_t d=0 ; d<dimension ; ++d ) {
	  wallVector[d]=vertexData[v1I][d] - vertexData[v2I][d];
	  wallLength +=wallVector[d]*wallVector[d];
	}
	if (wallLength>0.0) {
	  wallLength = std::sqrt(wallLength);
	}
	else {
	  std::cerr << "CenterTriangulation::VertexFromCellPressureLinear::derivs() "
		    << "Strange wall length." << std::endl;
	}
	wallFactor= ( dx[0]*wallVector[0]
		      +dx[1]*wallVector[1]
		      +dx[2]*wallVector[2])/(wallLength*wallLength);
	for( size_t d=0 ; d<dimension ; ++d ) {
	  
	  dx[d]=dx[d]-wallFactor*wallVector[d];
	  dxNorm += dx[d]*dx[d];
	}
	if (dxNorm>0.0) {
	  dxNorm = std::sqrt(dxNorm);
	  double dxInv = 1.0/dxNorm;
	  for( size_t d=0 ; d<dimension ; ++d ) {
	    dx[d] = dx[d]*dxInv;
	  }
	}
	else {
	  std::cerr << "CenterTriangulation::VertexFromCellPressureLinear::derivs() "
		    << "Force direction undetermined." << std::endl;
	}
	double factor = 0.5 * timeFactor_ * parameter(0);
	if (parameter(1) == 1)
	  {
	    //NOTE maybe this one should be calculated using the central mesh vertex?
	    double cellVolume = tmpCell.calculateVolume(vertexData);
	    factor /= std::fabs(cellVolume);         
	  }
	factor *= wallLength;    
	for( size_t d=0 ; d<dimension ; ++d ) {
	  vertexDerivs[v1I][d] += factor * dx[d];
	  vertexDerivs[v2I][d] += factor * dx[d];	
	}
      }
    }
  }
  
  void VertexFromCellPressureLinear::
  update(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 double h)
  {
    if (timeFactor_ < 1.0 ) {
      timeFactor_ += h/parameter(numParameter()-1);
    }
    if (timeFactor_ >1.0)
      timeFactor_=1.0;
    //cellData[0][12]=timeFactor_*parameter(0);
  }
}// end namespace CenterTriangulation
  

VertexFromCellPressureVolumeNormalized::
VertexFromCellPressureVolumeNormalized(std::vector<double> &paraValue, 
			   std::vector< std::vector<size_t> > 
			   &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //////////////////////////////////////////////////////////////////////
  if( paraValue.size()!=1 ) {
    std::cerr << "VertexFromCellPressureVolumeNormalized::"
	      << "VertexFromCellPressureVolumeNormalized() "
	      << "Uses one parameter K_force.\n";
    exit(0);
  }
  if( indValue.size() != 0 ) {
    std::cerr << "VertexFromCellPressureVolumeNormalized::"
	      << "VertexFromCellPressureVolumeNormalized() "
	      << "No index given.\n";
    exit(0);
  }
  //Set the variable values
  //////////////////////////////////////////////////////////////////////
  setId("VertexFromCellPressureVolumeNormalized");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //////////////////////////////////////////////////////////////////////
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "K_force";
  setParameterId( tmp );
}

void VertexFromCellPressureVolumeNormalized::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) {
  
  //Do the update for each vertex via each wall in each cell
  size_t numCells = T.numCell();
  size_t dimension = T.vertex(0).numPosition(); 
  
  //For each cell
  for( size_t cellI=0 ; cellI<numCells ; ++cellI ) {
    
    Cell &tmpCell = T.cell(cellI);
    //Calculate cell position from vertices
    std::vector<double> xCenter = tmpCell.positionFromVertex(vertexData);
    assert( xCenter.size()==dimension );
    double cellVolume = tmpCell.calculateVolume(vertexData);
    
    //Calculate derivative contributions to vertices from each wall
    for( size_t k=0 ; k<tmpCell.numWall() ; ++k ) {
      Wall &tmpWall = tmpCell.wallRef(k);
      size_t v1I = tmpWall.vertex1()->index();
      size_t v2I = tmpWall.vertex2()->index();
      std::vector<double> n(dimension),dx(dimension), x0(dimension);
      double b=0;
      for( size_t d=0 ; d<dimension ; ++d ) {
	n[d] = vertexData[v2I][d]-vertexData[v1I][d];
	b += n[d]*n[d];
	x0[d] = 0.5*(vertexData[v1I][d] + vertexData[v2I][d]);
	dx[d] = xCenter[d]-x0[d];
      }
      assert( b>0.0 );
      b = std::sqrt(b);
      for( size_t d=0 ; d<dimension ; ++d )
	n[d] /= b;
      double bInv = 1.0/b;
      double h = dx[0]*dx[0] + dx[1]*dx[1]
	-(n[0]*dx[0]+n[1]*dx[1])*(n[0]*dx[0]+n[1]*dx[1]);
      assert( h>0.0 );
      h = std::sqrt(h);
      double hInv = 1.0/h;
      double fac = parameter(0)*0.5/cellVolume;
      
      vertexDerivs[v1I][0] += fac*( 0.5*b*hInv*
				    ( -dx[0] - 2.0*(n[0]*dx[0]+n[1]*dx[1])
				      *(-n[1]*n[1]*bInv*dx[0]
					-0.5*n[0]
					+n[0]*n[1]*bInv*dx[1]) ) 
				    - h*n[0] ); 
      vertexDerivs[v2I][0] += fac*( 0.5*b*hInv*
				    ( -dx[0] - 2.0*(n[0]*dx[0]+n[1]*dx[1])
				      *(n[1]*n[1]*bInv*dx[0]
					-0.5*n[0]
					-n[0]*n[1]*bInv*dx[1]) )
				    + h*n[0] );
      vertexDerivs[v1I][1] += fac*( 0.5*b*hInv*
				    ( -dx[1] - 2.0*(n[0]*dx[0]+n[1]*dx[1])
				      *(-n[0]*n[0]*bInv*dx[1]
					-0.5*n[1]
					+n[0]*n[1]*bInv*dx[0]) ) 
				    - h*n[1] ); 
      vertexDerivs[v2I][1] += fac*( 0.5*b*hInv*
				    ( -dx[1] - 2.0*(n[0]*dx[0]+n[1]*dx[1])
				      *(n[0]*n[0]*bInv*dx[1]
					-0.5*n[1]
					-n[0]*n[1]*bInv*dx[0]) )
				    + h*n[1] );      
    }
  }
}

VertexFromCellPressureThresholdFromMaxPos::
VertexFromCellPressureThresholdFromMaxPos(std::vector<double> &paraValue, 
			   std::vector< std::vector<size_t> > 
			   &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //////////////////////////////////////////////////////////////////////
  if( paraValue.size()!=2 ) {
    std::cerr << "VertexFromCellPressureThresholdFromMaxPos::"
	      << "VertexFromCellPressureThresholdFromMaxPos() "
	      << "Uses two parameters K_force and X_th.\n";
    exit(0);
  }
  if( indValue.size() != 1 || indValue[0].size() != 1 ) {
    std::cerr << "VertexFromCellPressureThresholdFromMaxPos::"
	      << "VertexFromCellPressureThresholdFromMaxPos() "
	      << "One index given (direction).\n";
    exit(0);
  }
  //Set the variable values
  //////////////////////////////////////////////////////////////////////
  setId("VertexFromCellPressureThresholdFromMaxPos");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //////////////////////////////////////////////////////////////////////
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "K_force";
  tmp[1] = "X_th";
  setParameterId( tmp );
}

void VertexFromCellPressureThresholdFromMaxPos::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) {
  
  //Do the update for each vertex via each wall in each cell
  size_t numCells = T.numCell();
  size_t dimension = T.vertex(0).numPosition(); 
  
  //For each cell
  for( size_t cellI=0 ; cellI<numCells ; ++cellI ) {
    
    Cell &tmpCell = T.cell(cellI);
    //Calculate cell position from vertices
    std::vector<double> xCenter = tmpCell.positionFromVertex(vertexData);
    assert( xCenter.size()==dimension );
    assert( variableIndex(0,0)<dimension );
    //Find max pos in given direction;
    double max=vertexData[0][variableIndex(0,0)];
    for (size_t i=1; i<vertexData.size(); ++i ) {
      if ( vertexData[i][variableIndex(0,0)]>max )
	max = vertexData[i][variableIndex(0,0)];
    }
    //Only if close to apex
    if ( max-xCenter[variableIndex(0,0)]<parameter(1) ) {
      
      //Calculate derivative contributions to vertices from each wall
      for( size_t k=0 ; k<tmpCell.numWall() ; ++k ) {
	Wall &tmpWall = tmpCell.wallRef(k);
	size_t v1I = tmpWall.vertex1()->index();
	size_t v2I = tmpWall.vertex2()->index();
	std::vector<double> n(dimension),dx(dimension), x0(dimension);
	double b=0;
	for( size_t d=0 ; d<dimension ; ++d ) {
	  n[d] = vertexData[v2I][d]-vertexData[v1I][d];
	  b += n[d]*n[d];
	  x0[d] = 0.5*(vertexData[v1I][d] + vertexData[v2I][d]);
	  dx[d] = xCenter[d]-x0[d];
	}
	assert( b>0.0 );
	b = std::sqrt(b);
	for( size_t d=0 ; d<dimension ; ++d )
	  n[d] /= b;
	double bInv = 1.0/b;
	double h = dx[0]*dx[0] + dx[1]*dx[1]
	  -(n[0]*dx[0]+n[1]*dx[1])*(n[0]*dx[0]+n[1]*dx[1]);
	assert( h>0.0 );
	h = std::sqrt(h);
	double hInv = 1.0/h;
	double fac = parameter(0)*0.5;
	
	vertexDerivs[v1I][0] += fac*( 0.5*b*hInv*
				      ( -dx[0] - 2.0*(n[0]*dx[0]+n[1]*dx[1])
					*(-n[1]*n[1]*bInv*dx[0]
					  -0.5*n[0]
					  +n[0]*n[1]*bInv*dx[1]) ) 
				      - h*n[0] ); 
	vertexDerivs[v2I][0] += fac*( 0.5*b*hInv*
				      ( -dx[0] - 2.0*(n[0]*dx[0]+n[1]*dx[1])
					*(n[1]*n[1]*bInv*dx[0]
					  -0.5*n[0]
					  -n[0]*n[1]*bInv*dx[1]) )
				      + h*n[0] );
	vertexDerivs[v1I][1] += fac*( 0.5*b*hInv*
				      ( -dx[1] - 2.0*(n[0]*dx[0]+n[1]*dx[1])
					*(-n[0]*n[0]*bInv*dx[1]
					  -0.5*n[1]
					  +n[0]*n[1]*bInv*dx[0]) ) 
				      - h*n[1] ); 
	vertexDerivs[v2I][1] += fac*( 0.5*b*hInv*
				      ( -dx[1] - 2.0*(n[0]*dx[0]+n[1]*dx[1])
					*(n[0]*n[0]*bInv*dx[1]
					  -0.5*n[1]
					  -n[0]*n[1]*bInv*dx[0]) )
				      + h*n[1] );
	
	//vertexDerivs[v1I][0] += fac*( b*hInv*(bInv*(n[0]*dx[0]+n[1]*dx[1])*
	//				    (n[1]*n[1]*dx[0]+b*n[0]-n[0]*n[1]*dx[1]) + dx[0] ) - h*n[0] );
	//vertexDerivs[v1I][1] += fac*( b*hInv*(bInv*(n[1]*dx[1]+n[0]*dx[0])*
	//				    (-n[0]*n[0]*dx[1]+b*n[1]+n[1]*n[0]*dx[0]) + dx[1] ) - h*n[1] );
	//vertexDerivs[v2I][0] += fac*( n[0]*h + n[1]*hInv*(n[0]*dx[0]+n[1]*dx[1])*(n[0]*dx[1]-n[1]*dx[0]) );
	//vertexDerivs[v2I][1] += fac*( n[1]*h + n[0]*hInv*(n[1]*dx[1]+n[0]*dx[0])*(n[1]*dx[0]-n[0]*dx[1]) );      
      }
    }
  }
}

VertexFromCellPowerdiagram::
VertexFromCellPowerdiagram(std::vector<double> &paraValue, 
			   std::vector< std::vector<size_t> > 
			   &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //////////////////////////////////////////////////////////////////////
  if( paraValue.size()!=1 ) {
    std::cerr << "VertexFromCellPowerdiagram::"
	      << "VertexFromCellPowerdiagram() "
	      << "Uses one parameter K_force.\n";
    exit(0);
  }
  if( indValue.size() != 1 || indValue[0].size() != 1 ) {
    std::cerr << "VertexFromCellPowerdiagram::"
	      << "VertexFromCellPowerdiagram() "
	      << "Cell radius index given.\n";
    exit(0);
  }
  //Set the variable values
  //////////////////////////////////////////////////////////////////////
  setId("VertexFromCellPowerdiagram");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //////////////////////////////////////////////////////////////////////
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "K_force";
  setParameterId( tmp );
}

void VertexFromCellPowerdiagram::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) {
  
  //Do the update for each vertex
  //size_t numCells = T.numCell();
  size_t numVertices = T.numVertex();
  size_t cellRIndex = variableIndex(0,0);
  size_t dimension = T.vertex(0).numPosition(); 
  
  for( size_t i=0 ; i<numVertices ; ++i ) {
    //Only calculate if vertex connected to three cells
    size_t numCellForVertex = T.vertex(i).numCell();
    if(  numCellForVertex == 3 ) {
      //Calculate position for all three cells
      DataMatrix cellPos(numCellForVertex);
      std::vector<double> cellR(numCellForVertex);
      for( size_t cellI=0 ; cellI<numCellForVertex ; ++cellI ) {
	cellR[cellI] = cellData[ T.vertex(i).cell(cellI)->index() ][ cellRIndex ];
	//cellR[cellI] = T.cell( T.vertex(i).cell(cellI)->index() ).variable(0);
	cellPos[cellI] = T.cell( T.vertex(i).cell(cellI)->index() ).
	  positionFromVertex(vertexData);
      }
      //std::cerr << "* " << cellPos[0][0] << " " << cellPos[0][1] << "  "
      //	<< cellPos[1][0] << " " << cellPos[1][1] << "  "
      //	<< cellPos[2][0] << " " << cellPos[2][1] << std::endl;
      
      //Calculate optimal position from cell positions and raddii
      //according to power diagram formulation
      if( dimension != 2 ) {
	std::cerr << "VertexFromCellPowerdiagram::derivs() Only defined "
		  << "for two dimensions so far..." << std::endl;
	exit(-1);
      }
      double Kji = cellPos[1][0]*cellPos[1][0]+cellPos[1][1]*cellPos[1][1]-
	cellPos[0][0]*cellPos[0][0]-cellPos[0][1]*cellPos[0][1]-
	(cellR[1]*cellR[1]-cellR[0]*cellR[0]);
      double Kki = cellPos[2][0]*cellPos[2][0]+cellPos[2][1]*cellPos[2][1]-
	cellPos[0][0]*cellPos[0][0]-cellPos[0][1]*cellPos[0][1]-
	(cellR[2]*cellR[2]-cellR[0]*cellR[0]);
      std::vector<double> powPos(dimension);
      powPos[1] = 0.5*( Kji*(cellPos[2][0]-cellPos[0][0]) -
			Kki*(cellPos[1][0]-cellPos[0][0]) ) /
	( (cellPos[1][1]-cellPos[0][1])*(cellPos[2][0]-cellPos[0][0]) -
	  (cellPos[2][1]-cellPos[0][1])*(cellPos[1][0]-cellPos[0][0]) );
      powPos[0] = 0.5*Kji/(cellPos[1][0]-cellPos[0][0]) -
	powPos[1]*(cellPos[1][1]-cellPos[0][1])/(cellPos[1][0]-cellPos[0][0]);
      //std::cerr << i << " " << vertexData[i][0] << " " << vertexData[i][1] << "\t"
      //	<< powPos[0] << " " << powPos[1] << "\t"
      //	<< cellPos[0][0] << " " << cellPos[0][1] << "  "
      //	<< cellPos[1][0] << " " << cellPos[1][1] << "  "
      //	<< cellPos[2][0] << " " << cellPos[2][1] << "\t"
      //	<< cellR[0] << " " << cellR[1] << " " << cellR[2] << std::endl;
      
      //Update vertex for each dimension
      for(size_t d=0 ; d<dimension ; d++ )
	vertexDerivs[i][d] -= parameter(0)*(vertexData[i][d]-powPos[d]);
    }
  }
}

VertexFromCellInternalPressure::
VertexFromCellInternalPressure(std::vector<double> &paraValue, 
			       std::vector< std::vector<size_t> > 
			       &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //////////////////////////////////////////////////////////////////////
  if( paraValue.size()!=1 ) {
    std::cerr << "VertexFromCellInternalPressure::"
	      << "VertexFromCellInternalPressure() "
	      << "Uses one parameter K_force.\n";
    exit(0);
  }
  if( indValue.size() != 0 ) {
    std::cerr << "VertexFromCellInternalPressure::"
	      << "VertexFromCellInternalPressure() "
	      << "No index given.\n";
    exit(0);
  }
  //Set the variable values
  //////////////////////////////////////////////////////////////////////
  setId("VertexFromCellInternalPressure");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //////////////////////////////////////////////////////////////////////
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "K_force";
  setParameterId( tmp );
}

void VertexFromCellInternalPressure::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) {
  
  //Do the update for each vertex via each wall in each cell
  size_t numCells = T.numCell();
  size_t dimension = T.vertex(0).numPosition(); 
  
  //For each cell
  for( size_t cellI=0 ; cellI<numCells ; ++cellI ) {
    Cell &tmpCell = T.cell(cellI);
    if( !(tmpCell.isNeighbor(T.background())) ) {
      //Calculate cell position from vertices
      std::vector<double> xCenter = tmpCell.positionFromVertex(vertexData);
      assert( xCenter.size()==dimension );
      
      //Calculate derivative contributions to vertices from each wall
      for( size_t k=0 ; k<tmpCell.numWall() ; ++k ) {
	Wall &tmpWall = tmpCell.wallRef(k);
	size_t v1I = tmpWall.vertex1()->index();
	size_t v2I = tmpWall.vertex2()->index();
	std::vector<double> n(dimension),dx(dimension), x0(dimension);
	double b=0;
	for( size_t d=0 ; d<dimension ; ++d ) {
	  n[d] = vertexData[v2I][d]-vertexData[v1I][d];
	  b += n[d]*n[d];
	  x0[d] = 0.5*(vertexData[v1I][d] + vertexData[v2I][d]);
	  dx[d] = xCenter[d]-x0[d];
	}
	assert( b>0.0 );
	b = std::sqrt(b);
	for( size_t d=0 ; d<dimension ; ++d )
	  n[d] /= b;
	double bInv = 1.0/b;
	double h = dx[0]*dx[0] + dx[1]*dx[1]
	  -(n[0]*dx[0]+n[1]*dx[1])*(n[0]*dx[0]+n[1]*dx[1]);
	assert( h>0.0 );
	h = std::sqrt(h);
	double hInv = 1.0/h;
	double fac = parameter(0)*0.5;
	
	vertexDerivs[v1I][0] += fac*( 0.5*b*hInv*
				      ( -dx[0] - 2.0*(n[0]*dx[0]+n[1]*dx[1])
					*(-n[1]*n[1]*bInv*dx[0]
					  -0.5*n[0]
					  +n[0]*n[1]*bInv*dx[1]) ) 
				      - h*n[0] ); 
	vertexDerivs[v2I][0] += fac*( 0.5*b*hInv*
				      ( -dx[0] - 2.0*(n[0]*dx[0]+n[1]*dx[1])
					*(n[1]*n[1]*bInv*dx[0]
					  -0.5*n[0]
					  -n[0]*n[1]*bInv*dx[1]) )
				      + h*n[0] );
	vertexDerivs[v1I][1] += fac*( 0.5*b*hInv*
				      ( -dx[1] - 2.0*(n[0]*dx[0]+n[1]*dx[1])
					*(-n[0]*n[0]*bInv*dx[1]
					  -0.5*n[1]
					  +n[0]*n[1]*bInv*dx[0]) ) 
				      - h*n[1] ); 
	vertexDerivs[v2I][1] += fac*( 0.5*b*hInv*
				      ( -dx[1] - 2.0*(n[0]*dx[0]+n[1]*dx[1])
					*(n[0]*n[0]*bInv*dx[1]
					  -0.5*n[1]
					  -n[0]*n[1]*bInv*dx[0]) )
				      + h*n[1] );
	
	//vertexDerivs[v1I][0] += fac*( b*hInv*(bInv*(n[0]*dx[0]+n[1]*dx[1])*
	//				    (n[1]*n[1]*dx[0]+b*n[0]-n[0]*n[1]*dx[1]) + dx[0] ) - h*n[0] );
	//vertexDerivs[v1I][1] += fac*( b*hInv*(bInv*(n[1]*dx[1]+n[0]*dx[0])*
	//				    (-n[0]*n[0]*dx[1]+b*n[1]+n[1]*n[0]*dx[0]) + dx[1] ) - h*n[1] );
	//vertexDerivs[v2I][0] += fac*( n[0]*h + n[1]*hInv*(n[0]*dx[0]+n[1]*dx[1])*(n[0]*dx[1]-n[1]*dx[0]) );
	//vertexDerivs[v2I][1] += fac*( n[1]*h + n[0]*hInv*(n[1]*dx[1]+n[0]*dx[0])*(n[1]*dx[0]-n[0]*dx[1]) );      
      }
    }
  }
}

VertexForceOrigoFromIndex::
VertexForceOrigoFromIndex(std::vector<double> &paraValue, 
			  std::vector< std::vector<size_t> > 
			  &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //////////////////////////////////////////////////////////////////////
  if( paraValue.size()!=2 ) {
    std::cerr << "VertexForceOrigoFromIndex::"
	      << "VertexForceOrigoFromIndex() "
	      << "Uses two parameters K_force direction(-1 -> inwards)" 
	      << std::endl;
    exit(0);
  }
  if( paraValue[1] != 1.0 && paraValue[1] != -1.0 ) {
    std::cerr << "VertexForceOrigoFromIndex::"
	      << "VertexForceOrigoFromIndex() "
	      << "direction (second parameter) need to be 1 (outward) "
	      << "or -1 (inwards)." << std::endl;
    exit(0);
  }
  if( indValue.size() != 1 || indValue[0].size()<1 ) {
    std::cerr << "VertexForceOrigoFromIndex::"
	      << "VertexForceOrigoFromIndex() "
	      << "Vertex indices in first level." << std::endl;
    exit(0);
  }
  //Set the variable values
  //////////////////////////////////////////////////////////////////////
  setId("VertexForceOrigoFromIndex");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //////////////////////////////////////////////////////////////////////
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "K_force";
  tmp[1] = "direction";
  setParameterId( tmp );
}

void VertexForceOrigoFromIndex::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) {
  
  size_t dimension = vertexData[variableIndex(0,0)].size();
  double coeff = parameter(0)*parameter(1);
  //For each vertex in list
  for( size_t k=0 ; k<numVariableIndex(0) ; ++k ) {
    size_t i=variableIndex(0,k);
    for( size_t d=0 ; d<dimension ; d++ )
      vertexDerivs[i][d] -= coeff*vertexData[i][d];
  }
}

CellForceOrigoFromIndex::
CellForceOrigoFromIndex(std::vector<double> &paraValue, 
			  std::vector< std::vector<size_t> > 
			  &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //////////////////////////////////////////////////////////////////////
  if( paraValue.size()!=2 ) {
    std::cerr << "CellForceOrigoFromIndex::"
	      << "CellForceOrigoFromIndex() "
	      << "Uses two parameters K_force direction(-1 -> inwards)" 
	      << std::endl;
    exit(0);
  }
  if( paraValue[1] != 1.0 && paraValue[1] != -1.0 ) {
    std::cerr << "CellForceOrigoFromIndex::"
	      << "CellForceOrigoFromIndex() "
	      << "direction (second parameter) needs to be 1 (outward) "
	      << "or -1 (inwards)." << std::endl;
    exit(0);
  }
  if( indValue.size() != 1 || indValue[0].size()<1 ) {
    std::cerr << "CellForceOrigoFromIndex::"
	      << "CellForceOrigoFromIndex() "
	      << "Cell indices in first level." << std::endl;
    exit(0);
  }
  //Set the variable values
  //////////////////////////////////////////////////////////////////////
  setId("CellForceOrigoFromIndex");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //////////////////////////////////////////////////////////////////////
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "K_force";
  tmp[1] = "direction";
  setParameterId( tmp );
}

void CellForceOrigoFromIndex::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) {
  
  double coeff = parameter(0)*parameter(1);
  size_t dimension = vertexData[0].size();
  //For each vertex in cell list
  for( size_t k=0 ; k<numVariableIndex(0) ; ++k ) {
    size_t cellI=variableIndex(0,k);
    for( size_t l=0 ; l<T.cell(cellI).numVertex() ; ++l ) {
      size_t i = T.cell(cellI).vertex(l)->index();      
      for( size_t d=0 ; d<dimension ; d++ )
	vertexDerivs[i][d] += coeff*vertexData[i][d];
    }
  }
}

CylinderForce::
CylinderForce(std::vector<double> &paraValue, 
	      std::vector< std::vector<size_t> > 
	      &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //////////////////////////////////////////////////////////////////////
  if( paraValue.size()!=2 ) {
    std::cerr << "CylinderForce::"
	      << "CylinderForce() "
	      << "Uses two parameters K_force direction(-1 -> inwards)" 
	      << std::endl;
    exit(0);
  }
  if( paraValue[1] != 1.0 && paraValue[1] != -1.0 ) {
    std::cerr << "CylinderForce::"
	      << "CylinderForce() "
	      << "direction (second parameter) needs to be 1 (outward) "
	      << "or -1 (inwards)." << std::endl;
    exit(0);
  }
  if( indValue.size() != 0 ) {
    std::cerr << "CylinderForce::"
	      << "CylinderForce() "
	      << "No indices used." << std::endl;
    exit(0);
  }
  //Set the variable values
  //////////////////////////////////////////////////////////////////////
  setId("CylinderForce");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //////////////////////////////////////////////////////////////////////
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "K_force";
  tmp[1] = "direction";
  setParameterId( tmp );
}

void CylinderForce::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) {
  
  double coeff = parameter(0)*parameter(1);
  size_t dimension = vertexData[0].size();
  size_t lastPosIndex = dimension-1;
  //For each vertex 
  for( size_t i=0 ; i<T.numVertex() ; ++i ) {
		//On cylinder
		double norm = 0.0;
		for( size_t d=0 ; d<lastPosIndex ; d++ )
			norm += vertexData[i][d]*vertexData[i][d];
		if( norm>0.0 )
			norm = 1.0/std::sqrt(norm);
		else norm=0.0;
		for( size_t d=0 ; d<lastPosIndex ; d++ )
			vertexDerivs[i][d] += coeff*norm*vertexData[i][d];
	}
}

SphereCylinderForce::
SphereCylinderForce(std::vector<double> &paraValue, 
			  std::vector< std::vector<size_t> > 
			  &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //////////////////////////////////////////////////////////////////////
  if( paraValue.size()!=2 ) {
    std::cerr << "SphereCylinderForce::"
	      << "SphereCylinderForce() "
	      << "Uses two parameters K_force direction (-1 -> inwards)" 
	      << std::endl;
    exit(0);
  }
  if( paraValue[1] != 1.0 && paraValue[1] != -1.0 ) {
    std::cerr << "SphereCylinderForce::"
	      << "SphereCylinderForce() "
	      << "direction (second parameter) needs to be 1 (outward) "
	      << "or -1 (inwards)." << std::endl;
    exit(0);
  }
  if( indValue.size() != 0 ) {
    std::cerr << "SphereCylinderForce::"
	      << "SphereCylinderForce() "
	      << "Cell indices in first level." << std::endl;
    exit(0);
  }
  //Set the variable values
  //////////////////////////////////////////////////////////////////////
  setId("SphereCylinderForce");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //////////////////////////////////////////////////////////////////////
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "K_force";
  tmp[1] = "direction";
  setParameterId( tmp );
}

void SphereCylinderForce::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) {
  
  double coeff = parameter(0)*parameter(1);
  size_t dimension = vertexData[0].size();
  size_t lastPosIndex = dimension-1;
  //For each vertex 
  for( size_t i=0 ; i<T.numVertex() ; ++i ) {
    if( vertexData[i][lastPosIndex]>0.0 ) {
      //On sphere
      double norm = 0.0;
      for( size_t d=0 ; d<dimension ; d++ )
				norm += vertexData[i][d]*vertexData[i][d];
      if( norm>0.0 )
				norm = 1.0/std::sqrt(norm);
      else norm=0.0;
      for( size_t d=0 ; d<dimension ; d++ )
				vertexDerivs[i][d] += coeff*norm*vertexData[i][d];
    }
    else {
      //On cylinder
      double norm = 0.0;
      for( size_t d=0 ; d<dimension-1 ; d++ )
				norm += vertexData[i][d]*vertexData[i][d];
      if( norm>0.0 )
				norm = 1.0/std::sqrt(norm);
      else norm=0.0;
      for( size_t d=0 ; d<dimension-1 ; d++ )
				vertexDerivs[i][d] += coeff*norm*vertexData[i][d];
    }
  }
}

SphereCylinderForceFromRadius::
SphereCylinderForceFromRadius(std::vector<double> &paraValue, 
			      std::vector< std::vector<size_t> > 
			      &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //////////////////////////////////////////////////////////////////////
  if( paraValue.size()!=3 ) {
    std::cerr << "SphereCylinderForceFromRadius::"
	      << "SphereCylinderForceFromRadius() "
	      << "Uses three parameters F_out, F_in and radius." 
	      << std::endl;
    exit(0);
  }
  if( paraValue[0] < 0.0 || paraValue[1] < 0.0 || paraValue[2] < 0.0 ) {
    std::cerr << "SphereCylinderForceFromRadius::"
	      << "SphereCylinderForceFromRadius() "
	      << "Parameters need to be positive." << std::endl;
    exit(0);
  }
  if( indValue.size() != 0 ) {
    std::cerr << "SphereCylinderForceFromRadius::"
	      << "SphereCylinderForceFromRadius() "
	      << "No indices used." << std::endl;
    exit(0);
  }
  //Set the variable values
  //////////////////////////////////////////////////////////////////////
  setId("SphereCylinderForceFromRadius");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //////////////////////////////////////////////////////////////////////
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "F_out";
  tmp[1] = "F_in";
  tmp[2] = "R";
  setParameterId( tmp );
}

void SphereCylinderForceFromRadius::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) {
  
  size_t dimension = vertexData[0].size();
  size_t lastPosIndex = dimension-1;
  double F_out = parameter(0);
  double F_in = parameter(1);
  double R = parameter(2);
  //For each vertex 
  for( size_t i=0 ; i<T.numVertex() ; ++i ) {
    if( vertexData[i][lastPosIndex]>0.0 ) {
      //On sphere
      double norm = 0.0;
      for( size_t d=0 ; d<dimension ; d++ )
				norm += vertexData[i][d]*vertexData[i][d];
      norm = std::sqrt(norm);
      if( norm<R ) {
				double coeff = F_out*(R-norm)/norm;
				for( size_t d=0 ; d<dimension ; d++ )
					vertexDerivs[i][d] += coeff*vertexData[i][d];
      }
      else {
				double coeff = F_in*(R-norm)/norm;
				for( size_t d=0 ; d<dimension ; d++ )
					vertexDerivs[i][d] += coeff*vertexData[i][d];
      }
    }
    else {
      //On cylinder
      double norm = 0.0;
      for( size_t d=0 ; d<dimension-1 ; d++ )
				norm += vertexData[i][d]*vertexData[i][d];
      norm = std::sqrt(norm);
      if( norm<R ) {
				double coeff = F_out*(R-norm)/norm;
				for( size_t d=0 ; d<dimension-1 ; d++ )
					vertexDerivs[i][d] += coeff*vertexData[i][d];
      }
      else {
				double coeff = F_in*(R-norm)/norm;
				for( size_t d=0 ; d<dimension-1 ; d++ )
					vertexDerivs[i][d] += coeff*vertexData[i][d];
      }      
    }
  }
}

InfiniteWallForce::
InfiniteWallForce(std::vector<double> &paraValue, 
		  std::vector< std::vector<size_t> > 
		  &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //////////////////////////////////////////////////////////////////////
  if( paraValue.size()!=3 ) {
    std::cerr << "InfiniteWallForce::"
	      << "InfiniteWallForce() "
	      << "Uses three parameters k_spring threshold and direction "
	      << "(-1 -> less than)" 
	      << std::endl;
    exit(0);
  }
  if( paraValue[2] != 1.0 && paraValue[2] != -1.0 ) {
    std::cerr << "InfiniteWallForce::"
	      << "InfiniteWallForce() "
	      << "direction (second parameter) need to be 1 (greater than) "
	      << "or -1 (less than)." << std::endl;
    exit(0);
  }
  if( indValue.size() != 1 || indValue[0].size() != 1 ) {
    std::cerr << "InfiniteWallForce::"
	      << "InfiniteWallForce() "
	      << "Pos index in first level." << std::endl;
    exit(0);
  }
  //Set the variable values
  //////////////////////////////////////////////////////////////////////
  setId("InfiniteWallForce");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //////////////////////////////////////////////////////////////////////
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "k_spring";
  tmp[1] = "threshold";
  tmp[2] = "direction";
  setParameterId( tmp );
}

void InfiniteWallForce::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) {
  
  //Check the cancelation for every vertex
  size_t numVertices = T.numVertex();
  size_t posIndex = variableIndex(0,0);
  //size_t dimension = vertexData[0].size();
  assert( posIndex < vertexData[0].size() );
  
  for( size_t i=0 ; i<numVertices ; ++i ) {
    if( parameter(2)>0 && vertexData[i][posIndex]>parameter(1) ) {
      vertexDerivs[i][posIndex] -= parameter(0)*
	(vertexData[i][posIndex]-parameter(1));
    }
    else if( parameter(2)<0 && vertexData[i][posIndex]<parameter(1) ) {
	    vertexDerivs[i][posIndex] -= parameter(0)*
		    (vertexData[i][posIndex]-parameter(1));
    }
  }
}

EpidermalVertexForce::
EpidermalVertexForce(std::vector<double> &paraValue, 
		  std::vector< std::vector<size_t> > 
		  &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //////////////////////////////////////////////////////////////////////
  if( paraValue.size()!=2 ) {
    std::cerr << "EpidermalVertexForce::"
	      << "EpidermalVertexForce() "
	      << "Uses two parameters k_strength direction "
	      << "(-1 -> down)" 
	      << std::endl;
    exit(0);
  }
  if( paraValue[1] != 1.0 && paraValue[1] != -1.0 ) {
    std::cerr << "EpidermalVertexForce::"
	      << "EpidermalVertexForce() "
	      << "direction (second parameter) need to be 1 (pos) "
	      << "or -1 (neg direction)." << std::endl;
    exit(0);
  }
  if( indValue.size() != 1 || indValue[0].size() != 1 ) {
    std::cerr << "EpidermalVertexForce::"
	      << "EpidermalVertexForce() "
	      << "Pos index in first level." << std::endl;
    exit(0);
  }
  //Set the variable values
  //////////////////////////////////////////////////////////////////////
  setId("EpidermalVertexForce");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //////////////////////////////////////////////////////////////////////
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "k_spring";
  tmp[1] = "direction";
  setParameterId( tmp );
}

void EpidermalVertexForce::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) {
  
  //Check the cancelation for every vertex
  size_t numVertices = T.numVertex();
  size_t posIndex = variableIndex(0,0);
  //size_t dimension = vertexData[0].size();
  assert( posIndex < vertexData[0].size() );
  
  for( size_t i=0 ; i<numVertices ; ++i ) {
		int epidermisFlag=0;
		for( size_t k=0 ; k<T.vertex(i).numWall() ; ++k ) {
			if( T.vertex(i).wall(k)->cell1()==T.background() ||
					T.vertex(i).wall(k)->cell1()==T.background() ) {
				epidermisFlag++;
			}
		}		
    if( epidermisFlag ) {
      vertexDerivs[i][posIndex] += parameter(0)*parameter(1);
    }
  }
}

VertexFromPressureExperimental::VertexFromPressureExperimental(std::vector<double> &paraValue, 
												   std::vector< std::vector<size_t> > &indValue)
{
	if (paraValue.size() != 2) {
	  std::cerr << "VertexFromPressureExperimental::VertexFromPressureExperimental() "
		    << "Uses two parameter: k no_contraction_flag" << std::endl;
	  exit(EXIT_FAILURE);
	}
	
	if (indValue.size() != 1 || indValue[0].size() != 1) {
	  std::cerr << "VertexFromPressureExperimental::VertexFromPressureExperimental() "
		    << "Water volume index given.\n";
	  exit(EXIT_FAILURE);
	}
	
	setId("VertexFromPressureExperimental");
	setParameter(paraValue);  
	setVariableIndex(indValue);
	
	std::vector<std::string> tmp(numParameter());
	tmp[0] = "k";
	tmp[1] = "no_contraction_flag";
	setParameterId(tmp);
}

void VertexFromPressureExperimental::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs)
{
  // ATTENTION! This function assume two dimensions (x, y) and that
  // the vertices are sorted in a clock-wise order.
  
  assert(T.vertex(0).numPosition());
  
  double area;
  for (size_t n = 0; n < T.numCell(); ++n) {
    Cell cell = T.cell(n);
    std::vector< std::pair<double, double> > vertices;
    for (size_t i = 0; i < cell.numVertex(); ++i) {
      std::pair<double, double> vertex;
      vertex.first = vertexData[cell.vertex(i)->index()][0];
      vertex.second = vertexData[cell.vertex(i)->index()][1];
      vertices.push_back(vertex);
    }
    area = polygonArea(vertices);
    double sa = area < 0 ? -1 : +1;
    area = std::fabs(area);
    
    if (parameter(1) != 0.0 && cellData[n][variableIndex(0, 0)] - area < 0) {
      //Std::cerr << "Lower water volume than volume: " << n << " " 
      //				<< area << " " << cellData[n][variableIndex(0,0)]
      //				<< std::endl;
      continue;
    }
    for (size_t i = 1; i < (cell.numVertex() + 1); ++i) {
      Vertex *vertex = cell.vertex(i % cell.numVertex()); // Current vertex in polygon.
      Vertex *pvertex = cell.vertex((i - 1) % cell.numVertex()); // Previous vertex in polygon.
      Vertex *nvertex = cell.vertex((i + 1) % cell.numVertex()); // Next vertex in polygon.
      
      double px = vertexData[pvertex->index()][0];
      double py = vertexData[pvertex->index()][1];

      double nx = vertexData[nvertex->index()][0];
      double ny = vertexData[nvertex->index()][1];
      
      double dAdx = sa * 0.5 * (-py + ny);
      double dAdy = sa * 0.5 * (px - nx);

      
      vertexDerivs[vertex->index()][0] += parameter(0) * (1 - area / cellData[n][variableIndex(0, 0)]) * dAdx;
      vertexDerivs[vertex->index()][1] += parameter(0) * (1 - area / cellData[n][variableIndex(0, 0)]) * dAdy;

    }
  }
}

double VertexFromPressureExperimental::polygonArea(std::vector< std::pair<double, double> > vertices)
{
  double area = 0.0;
  size_t N = vertices.size();
  for (size_t n = 0; n < N; ++n) {
    area += vertices[n].first * vertices[(n + 1) % N].second;
    area -= vertices[(n + 1) % N].first * vertices[n].second;
  }
  area *= 0.5;
  return area;
}

CellVolumeExperimental::
CellVolumeExperimental(std::vector<double> &paraValue,
		       std::vector< std::vector<size_t> > &indValue)
{
  if (paraValue.size() != 4) {
    std::cerr << "CellVolumeExperimental::CellVolumeExperimental() "
	      << "Uses four parameters: k_p, P_max, k_pp and "
	      << "allowShrink_flag." << std::endl;
    exit(EXIT_FAILURE);
  }
  
  if (indValue.size() < 2 || indValue.size() > 3
      || indValue[0].size() != 2 
      || (indValue.size()==3 && indValue[2].size() != 1 ) ) {
    std::cerr << "CellVolumeExperimental::CellVolumeExperimental() "
	      << "Wall length index and cell volume index must be given in "
	      << "first level.\n"
	      << "Force indices must be given in second level. "
	      << "Optionally index for saving the pressure can be"
	      << " given at third level." << std::endl; 		
    exit(EXIT_FAILURE);
  }
  
  setId("CellVolumeExperimental");
  setParameter(paraValue);
  setVariableIndex(indValue);
  
  std::vector<std::string> tmp(numParameter());
  tmp[0] = "k_p";
  tmp[1] = "P_max";
  tmp[2] = "k_pp";
  tmp[3] = "allowShrink_flag";
  
  setParameterId(tmp);
}

void CellVolumeExperimental::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs)
{
  for (size_t n = 0; n < T.numCell(); ++n) {
    Cell cell = T.cell(n);
    
    double P = 0.0;
    double sum = 0.0;
    for (size_t i = 0; i < cell.numWall(); ++i) {
      size_t vertex1Index = cell.wall(i)->vertex1()->index();
      size_t vertex2Index = cell.wall(i)->vertex2()->index();
      size_t dimensions = vertexData[vertex1Index].size();
      
      double distance = 0.0;
      for (size_t d = 0; d < dimensions; ++d) {
	distance += (vertexData[vertex1Index][d] - vertexData[vertex2Index][d])
	  * (vertexData[vertex1Index][d] - vertexData[vertex2Index][d]);
      }
      distance = std::sqrt(distance);
      
      for (size_t j = 0; j < numVariableIndex(1); ++j)
	P += wallData[cell.wall(i)->index()][variableIndex(1, j)]/distance;
      sum += distance;
    }
    P *= parameter(2);
    
    if (numVariableIndexLevel()==3)
      cellData[n][variableIndex(2,0)]=P;
    
    if( parameter(3) || parameter(1)-P>0.0 )
      cellDerivs[cell.index()][variableIndex(0, 1)] += 
	parameter(0) * (parameter(1) - P) * sum;
  }
}

EpidermalRadialForce::EpidermalRadialForce(std::vector<double> &paraValue,
					   std::vector< std::vector<size_t> > &indValue)
{
  if (paraValue.size() != 1) {
    std::cerr << "EpidermalRadialForce::EpidermalRadialForce() " 
	      << "Uses one parameter: force" << std::endl;
    exit(EXIT_FAILURE);
  }
  
  if (indValue.size() != 0) {
    std::cerr << "EpidermalRadialForce::EpidermalRadialForce() "
	      << "No indices are given." << std::endl;
    exit(EXIT_FAILURE);
  }
  
  setId("EpidermalRadialForce");
  setParameter(paraValue);
  setVariableIndex(indValue);
  
  std::vector<std::string> tmp(numParameter());
  tmp[0] = "force";
  
  setParameterId(tmp);
}

void EpidermalRadialForce::derivs(Tissue &T,
				  DataMatrix &cellData,
				  DataMatrix &wallData,
				  DataMatrix &vertexData,
				  DataMatrix &cellDerivs,
				  DataMatrix &wallDerivs,
				  DataMatrix &vertexDerivs)
{
  for (size_t n = 0; n < T.numVertex(); ++n) {
    Vertex vertex = T.vertex(n);
    
    // If vertex is not in the epidermal layer then skip.
    bool isEpidermalVertex = false;
    for (size_t i = 0; i < vertex.numWall(); ++i) {
      Wall *wall = vertex.wall(i);
      if ((wall->cell1()->index() == (size_t) -1) || (wall->cell2()->index() == (size_t) -1)) {
	isEpidermalVertex = true;
	break;
      }
    }
    if (isEpidermalVertex == false) {
      continue;
      
    }
    
    double x = vertex.position(0);
    double y = vertex.position(1);
    double A = std::sqrt(x * x + y * y);
    if (A == 0)
      continue;
    x /= A;
    y /= A;
    
    // 		std::cerr << "Vertex " << vertex.index() << std::endl;
    // 		std::cerr << " x = " << x << std::endl;
    // 		std::cerr << " y = " << y << std::endl;
    
    vertexDerivs[vertex.index()][0] += -parameter(0) * x;
    vertexDerivs[vertex.index()][1] += -parameter(0) * y;
  }
}

PerpendicularWallPressure::PerpendicularWallPressure(std::vector<double> &paraValue,
						     std::vector< std::vector<size_t> > &indValue)
{
  if (paraValue.size() != 1) {
    std::cerr << "PerpendicularWallPressure::PerpendicularWallPressure() " 
	      << "Uses one parameter: k_force" << std::endl;
		exit(EXIT_FAILURE);
  }
  
  if (indValue.size() != 1 || indValue[0].size() != 1) {
    std::cerr << "PerpendicularWallPressure::PerpendicularWallPressure() "
	      << "First level gives pressure index in cell." << std::endl;
    exit(EXIT_FAILURE);
  }
  
  setId("PerpendicularWallPressure");
  setParameter(paraValue);
  setVariableIndex(indValue);
  
  std::vector<std::string> tmp(numParameter());
  tmp[0] = "k_force";
  
  setParameterId(tmp);
}

void PerpendicularWallPressure::derivs(Tissue &T,
				       DataMatrix &cellData,
				       DataMatrix &wallData,
				       DataMatrix &vertexData,
				       DataMatrix &cellDerivs,
				       DataMatrix &wallDerivs,
				       DataMatrix &vertexDerivs)
{
  size_t dimension = vertexData[0].size();
  if (dimension==2) {
    for (size_t n = 0; n < T.numWall(); ++n) {
      Wall wall = T.wall(n);
      std::vector<Cell *> cells(2);
      cells[0] = wall.cell1();
      cells[1] = wall.cell2();
      std::vector<Vertex *> vertices(2);
      vertices[0] = wall.vertex1();
      vertices[1] = wall.vertex2();
      
      double nx = vertexData[vertices[1]->index()][0] - vertexData[vertices[0]->index()][0];
      double ny = vertexData[vertices[1]->index()][1] - vertexData[vertices[0]->index()][1];
      double A = std::sqrt(nx * nx + ny * ny);
      nx = nx / A;
      ny = ny / A;
      
      for (size_t i = 0; i < cells.size(); ++i) {
	Cell *cell = cells[i];
	
	if (cell == T.background()) {
	  continue;
	}
	
	double force = parameter(0) * 
	  cellData[cell->index()][variableIndex(0, 0)] *  
	  wall.lengthFromVertexPosition(vertexData);
	
	double x = 0.0;
	double y = 0.0;
	for (size_t j = 0; j < cell->numVertex(); ++j) {
	  x += vertexData[cell->vertex(j)->index()][0];
	  y += vertexData[cell->vertex(j)->index()][1];
	}
	x = (x / cell->numVertex()) - vertexData[vertices[0]->index()][0];
	y = (y / cell->numVertex()) - vertexData[vertices[0]->index()][1];
	
	int sign = (nx * y - ny * x) >= 0 ? 1 : -1;
	
	for (size_t k = 0; k < vertices.size(); ++k) {
	  Vertex *vertex = vertices[k];
	  vertexDerivs[vertex->index()][0] += ny * sign * force;
	  vertexDerivs[vertex->index()][1] += - nx * sign * force;
	}
      }
    }
  }
  else if (dimension==3) {
    
    size_t numCell = T.numCell();
    for (size_t i=0; i<numCell; ++i) {
      Cell* cell1 = &(T.cell(i));
      size_t c1I = cell1->index();
      std::vector<double> n1 = cell1->getNormalToPCAPlane();
      int n1Sign = cell1->vectorSignFromSort(n1,vertexData);
      std::vector<double> n = n1;
      std::vector<double> x1 = cell1->positionFromVertex(vertexData);
      
      size_t numNeigh = cell1->numWall();
      for (size_t k=0; k<numNeigh; ++k) {
	Cell* cell2 = cell1->cellNeighbor(k);
	size_t c2I = cell2->index();
	// Get normal to second cell (if not background) and calculate combined normal (n)
	if (cell2!=T.background() && c2I<c1I) {
	  continue;
	}
	else if (cell2!=T.background()) {
	  std::vector<double> n2(dimension);
	  n2 = cell2->getNormalToPCAPlane();
	  int n2Sign = cell2->vectorSignFromSort(n2,vertexData);
	  for (size_t d=0; d<dimension; ++d) 
	    n[d] = n1Sign*n1[d]+n2Sign*n2[d];
	}
	// Get wall direction
	size_t v1I = cell1->wall(k)->vertex1()->index();
	size_t v2I = cell1->wall(k)->vertex2()->index();
	std::vector<double> nw(dimension),xw(dimension);
	for (size_t d=0; d<dimension; ++d) {
	  nw[d] = vertexData[v2I][d]-vertexData[v1I][d];
	  xw[d] = 0.5*(vertexData[v2I][d]+vertexData[v1I][d]);					
	}
	// Normalize n and nw
	double norm=0.0,normW=0.0;
	for (size_t d=0; d<dimension; ++d) { 
	  norm += n[d]*n[d];
	  normW += nw[d]*nw[d];
	}
	assert(norm>0.0 && normW>0.0);
	norm = std::sqrt(norm);
	normW = std::sqrt(normW);
	double normFac=1.0/norm,normWFac=1.0/normW;
	for (size_t d=0; d<dimension; ++d) {
	  n[d] *= normFac;
	  nw[d] *= normWFac;
	}
	
	// Calculate force direction perpendicular to n and nw, plus sign
	std::vector<double> nF(dimension);
	nF[0] = n[1]*nw[2]-n[2]*nw[1];
	nF[1] = n[2]*nw[0]-n[0]*nw[2];
	nF[2] = n[0]*nw[1]-n[1]*nw[0];
	norm=0.0;
	for (size_t d=0; d<dimension; ++d) { 
	  norm += nF[d]*nF[d];
	}
	assert(norm>0.0);
	norm = std::sqrt(norm);
	normFac=1.0/norm;
	for (size_t d=0; d<dimension; ++d) {
	  nF[d] *= normFac;
	}
	
	// sign from scalar product?
	double scalar=0.0;
	for (size_t d=0; d<dimension; ++d)
	  scalar += (xw[d]-x1[d])*nF[d];
	int sign=1;
	if (scalar<0.0)
	  sign=-1;
	double forceFactor = parameter(0)*cellData[c1I][variableIndex(0,0)]*normW;
	if (cell2!=T.background())
	  forceFactor -= parameter(0)*cellData[c2I][variableIndex(0,0)]*normW;
	
	// Print result (wall and nF) for debugging
	//double fac=0.5*sign;
	//std::cerr << c1I << " " << vertexData[v1I][0] << " " 
	//				<< vertexData[v1I][1] << " " 
	//				<< vertexData[v1I][2] << " " 
	//				<< xw[0] << " " << xw[1] << " " << xw[2] << std::endl;
	//std::cerr << c1I << " " << vertexData[v2I][0] << " " 
	//				<< vertexData[v2I][1] << " " 
	//				<< vertexData[v2I][2] << " " 
	//				<< xw[0]+fac*nF[0] << " " << xw[1]+fac*nF[1] << " " << xw[2]+fac*nF[2] 
	//				<< std::endl << std::endl << std::endl;
	
	for (size_t d=0; d<dimension; ++d) {
	  vertexDerivs[v1I][d] += forceFactor*sign*nF[d];
	  vertexDerivs[v2I][d] += forceFactor*sign*nF[d];
	}				
      }
    }
  }
  else {
    std::cerr << "PerpendicularWallPressure::derivs() Only applicable for two or three"
	      << " dimensions." << std::endl;
    exit(-1);
  }
}

VertexFromCellPlane::
VertexFromCellPlane(std::vector<double> &paraValue,
		    std::vector< std::vector<size_t> > &indValue)
{
  if (paraValue.size() != 2) {
    std::cerr << "VertexFromCellPlane::VertexFromCellPlane() " 
	      << "Uses two parameters: k_force and areaFlag" << std::endl;
    exit(EXIT_FAILURE);
  }
  if (paraValue[1]!=0.0 && paraValue[1]!=1.0) {
    std::cerr << "VertexFromCellPlane::VertexFromCellPlane() " 
	      << "areaFlag must be zero (no area included) or one (area included)." << std::endl;
    exit(EXIT_FAILURE);
  }
  if (indValue.size() != 0) {
    std::cerr << "VertexFromCellPlane::VertexFromCellPlane() " 
	      << std::endl
	      << "No variable index used." << std::endl;
    exit(EXIT_FAILURE);
  }
  
  setId("VertexFromCellPlane");
  setParameter(paraValue);
  setVariableIndex(indValue);
  
  std::vector<std::string> tmp(numParameter());
  tmp[0] = "k_force";
  tmp[1] = "areaFlag";
  
  setParameterId(tmp);
}

void VertexFromCellPlane::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs)
{
  size_t dimension = vertexData[0].size();
  if (dimension!=3) {
    std::cerr << "VertexFromCellPlane::VertexFromCellPlane() " 
	      << "Only implemented for three dimensions." << std::endl;
    exit(EXIT_FAILURE);
  }
  unsigned int numFlipNormal=0;
  for (size_t n = 0; n < T.numCell(); ++n) {
    Cell cell = T.cell(n);
    unsigned int flipFlag=0;
    std::vector<double> normal = cell.getNormalToPCAPlane();
    double norm=0.0;
    for (size_t d=0; d<dimension; ++d)
      norm += normal[d]*normal[d];
    if (norm != 1.0) {
      norm = std::sqrt(norm);
      assert(norm>0.0);
      double normFac = 1.0/norm; 
      for (size_t d=0; d<dimension; ++d)
	normal[d] *= normFac;
    }
    
    std::vector<int> scalarProdSign(cell.numVertex());
    std::vector<double> scalarProdVal(cell.numVertex());
    double scalarProdSum=0.0;
    for (size_t k=0; k<cell.numVertex(); ++k) {
      size_t k2=(k+1)%cell.numVertex();
      size_t k3=(k+2)%cell.numVertex();
      //Make sure ((v2-v1)x(v3-v2))n has same sign for all cells
      assert(cell.numVertex()>2);
      std::vector<double> nw1(dimension),nw2(dimension);
      for (size_t d=0; d<dimension; ++d) {
	nw1[d] = vertexData[cell.vertex(k2)->index()][d]-vertexData[cell.vertex(k)->index()][d];
	nw2[d] = vertexData[cell.vertex(k3)->index()][d]-vertexData[cell.vertex(k2)->index()][d];
      }
      //cross product
      double scalarProd=0.0;
      for (size_t d1=0; d1<dimension; ++d1) {
	size_t d2=(d1+1)%dimension;
	size_t d3=(d1+2)%dimension;
	scalarProd += (nw1[d1]*nw2[d2]-nw1[d2]*nw2[d1])*normal[d3];
      }
      scalarProdVal[k] = scalarProd;
      scalarProdSum += scalarProd;
      if (scalarProd>0.0)
	scalarProdSign[k]=1;
      else
	scalarProdSign[k]=-1;
    }
    // for (size_t k=1; k<cell.numVertex(); ++k)
    //   if (scalarProdSign[k]!=scalarProdSign[0]) {
    //     std::cerr << "Cell " << n << " has diverging signs on scalar product." << std::endl;
    //     break;
    //   }
    int scalarProdSignSum=0;
    for (size_t k=0; k<scalarProdSign.size(); ++k)
      scalarProdSignSum += scalarProdSign[k];
    
    if (scalarProdSignSum<0) {
      numFlipNormal++;
      flipFlag=1;
      for (size_t d=0; d<dimension; ++d)
	normal[d] = -normal[d];
    }
    else if (scalarProdSignSum==0) {
      std::cerr << "Cell " << n << " has no majority sign in right hand rule expression." 
		<< std::endl;
      if (std::fabs(scalarProdSum)>0.01) {
	if (scalarProdSum<0.0) {
	  numFlipNormal++;
	  flipFlag=1;
	  for (size_t d=0; d<dimension; ++d)
	    normal[d] = -normal[d];
	}
      }
      else {
	std::vector<double> center = cell.positionFromVertex(vertexData);
	// Print all walls
	for (size_t k=0; k<cell.numWall(); ++k) {
	  std::cerr << "0 "; 
	  for (size_t d=0; d<dimension; ++d)
	    std::cerr << vertexData[cell.wall(k)->vertex1()->index()][d] << " ";
	  std::cerr << std::endl << "0 "; 
	  for (size_t d=0; d<dimension; ++d)
	    std::cerr << vertexData[cell.wall(k)->vertex2()->index()][d] << " ";
	  std::cerr << std::endl << std::endl;
	}		 
	std::cerr << "1 "; 
	for (size_t d=0; d<dimension; ++d)
	  std::cerr << vertexData[cell.wall(0)->vertex1()->index()][d] << " ";
	std::cerr << std::endl << std::endl;
	std::cerr << "2 "; 
	for (size_t d=0; d<dimension; ++d)
	  std::cerr << center[d] << " ";
	std::cerr << std::endl << "2 ";
	for (size_t d=0; d<dimension; ++d)
	  std::cerr << center[d]+normal[d] << " ";
	std::cerr << std::endl << std::endl;
	std::cerr << "3 "; 
	for (size_t d=0; d<dimension; ++d)
	  std::cerr << normal[d] << " ";
	std::cerr << std::endl << std::endl;
	std::cerr << "4 "; 
	for (size_t i=0; i<scalarProdVal.size(); ++i)
	  std::cerr << scalarProdVal[i] << " ";
	std::cerr << std::endl << std::endl;
	std::cerr << "5 "; 
	for (size_t i=0; i<scalarProdSign.size(); ++i)
	  std::cerr << scalarProdSign[i] << " ";
	std::cerr << std::endl << std::endl;
	
	exit(-1);
      }
    }
    // Get the cell size
    double A=1.0;
    if (parameter(1)==1.0)
      A = cell.calculateVolume(vertexData)/cell.numVertex();
    
    double coeff = parameter(0) * A;
    //update the vertex derivatives
    for (size_t k=0; k<cell.numVertex(); ++k) {
      double vCoeff=coeff;
      if (cell.vertex(k)->isBoundary(T.background()))
	vCoeff *= 1.5;
      for (size_t d=0; d<dimension; ++d) {
	vertexDerivs[cell.vertex(k)->index()][d] += vCoeff * normal[d];
      }
    }	
    // For saving normals in direction used for test plotting
    //for (size_t d=0; d<dimension; ++d)
    //cellData[cell.index()][d] = normal[d];
  }
  //std::cerr << numFlipNormal << " cells out of " << T.numCell() << " has flipped normal."
  //	      << std::endl;
}
///////////////////////////////////////////////////////////////////////////////////////////////77

VertexFromCellPlaneLinear::
VertexFromCellPlaneLinear(std::vector<double> &paraValue,
		    std::vector< std::vector<size_t> > &indValue)
{
  if (paraValue.size() != 3) {
    std::cerr << "VertexFromCellPlaneLinear::VertexFromCellPlaneLinear() " 
	      << "Uses two parameters: k_force and areaFlag and time span" << std::endl;
    exit(EXIT_FAILURE);
  }
  if (paraValue[1]!=0.0 && paraValue[1]!=1.0) {
    std::cerr << "VertexFromCellPlaneLinear::VertexFromCellPlaneLinear() " 
	      << "areaFlag must be zero (no area included) or one (area included)." << std::endl;
    exit(EXIT_FAILURE);
  }
  if (paraValue[1]<=0.0 ) {
    std::cerr << "VertexFromCellPlaneLinear::VertexFromCellPlaneLinear() " 
	      << "time span must be positive." << std::endl;
    exit(EXIT_FAILURE);
  }
  if (indValue.size() != 0) {
    std::cerr << "VertexFromCellPlaneLinear::VertexFromCellPlane() " 
	      << std::endl
	      << "No variable index used." << std::endl;
    exit(EXIT_FAILURE);
  }
  
  setId("VertexFromCellPlaneLinear");
  setParameter(paraValue);
  setVariableIndex(indValue);
  
  std::vector<std::string> tmp(numParameter());
  tmp[0] = "k_force";
  tmp[1] = "areaFlag";
  tmp[2] = "deltaT";

  timeFactor_=0.0;
  setParameterId(tmp);
}

void VertexFromCellPlaneLinear::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs)
{
  size_t dimension = vertexData[0].size();
  if (dimension!=3) {
    std::cerr << "VertexFromCellPlaneLinear::VertexFromCellPlaneLinear() " 
	      << "Only implemented for three dimensions." << std::endl;
    exit(EXIT_FAILURE);
  }
  unsigned int numFlipNormal=0;
  for (size_t n = 0; n < T.numCell(); ++n) {
    Cell cell = T.cell(n);
    unsigned int flipFlag=0;
    std::vector<double> normal = cell.getNormalToPCAPlane();
    double norm=0.0;
    for (size_t d=0; d<dimension; ++d)
      norm += normal[d]*normal[d];
    if (norm != 1.0) {
      norm = std::sqrt(norm);
      assert(norm>0.0);
      double normFac = 1.0/norm; 
      for (size_t d=0; d<dimension; ++d)
	normal[d] *= normFac;
    }
    
    std::vector<int> scalarProdSign(cell.numVertex());
    std::vector<double> scalarProdVal(cell.numVertex());
    double scalarProdSum=0.0;
    for (size_t k=0; k<cell.numVertex(); ++k) {
      size_t k2=(k+1)%cell.numVertex();
      size_t k3=(k+2)%cell.numVertex();
      //Make sure ((v2-v1)x(v3-v2))n has same sign for all cells
      assert(cell.numVertex()>2);
      std::vector<double> nw1(dimension),nw2(dimension);
      for (size_t d=0; d<dimension; ++d) {
	nw1[d] = vertexData[cell.vertex(k2)->index()][d]-vertexData[cell.vertex(k)->index()][d];
	nw2[d] = vertexData[cell.vertex(k3)->index()][d]-vertexData[cell.vertex(k2)->index()][d];
      }
      //cross product
      double scalarProd=0.0;
      for (size_t d1=0; d1<dimension; ++d1) {
	size_t d2=(d1+1)%dimension;
	size_t d3=(d1+2)%dimension;
	scalarProd += (nw1[d1]*nw2[d2]-nw1[d2]*nw2[d1])*normal[d3];
      }
      scalarProdVal[k] = scalarProd;
      scalarProdSum += scalarProd;
      if (scalarProd>0.0)
	scalarProdSign[k]=1;
      else
	scalarProdSign[k]=-1;
    }
    // for (size_t k=1; k<cell.numVertex(); ++k)
    //   if (scalarProdSign[k]!=scalarProdSign[0]) {
    //     std::cerr << "Cell " << n << " has diverging signs on scalar product." << std::endl;
    //     break;
    //   }
    int scalarProdSignSum=0;
    for (size_t k=0; k<scalarProdSign.size(); ++k)
      scalarProdSignSum += scalarProdSign[k];
    
    if (scalarProdSignSum<0) {
      numFlipNormal++;
      flipFlag=1;
      for (size_t d=0; d<dimension; ++d)
	normal[d] = -normal[d];
    }
    else if (scalarProdSignSum==0) {
      std::cerr << "Cell " << n << " has no majority sign in right hand rule expression." 
		<< std::endl;
      if (std::fabs(scalarProdSum)>0.01) {
	if (scalarProdSum<0.0) {
	  numFlipNormal++;
	  flipFlag=1;
	  for (size_t d=0; d<dimension; ++d)
	    normal[d] = -normal[d];
	}
      }
      else {
	std::vector<double> center = cell.positionFromVertex(vertexData);
	// Print all walls
	for (size_t k=0; k<cell.numWall(); ++k) {
	  std::cerr << "0 "; 
	  for (size_t d=0; d<dimension; ++d)
	    std::cerr << vertexData[cell.wall(k)->vertex1()->index()][d] << " ";
	  std::cerr << std::endl << "0 "; 
	  for (size_t d=0; d<dimension; ++d)
	    std::cerr << vertexData[cell.wall(k)->vertex2()->index()][d] << " ";
	  std::cerr << std::endl << std::endl;
	}		 
	std::cerr << "1 "; 
	for (size_t d=0; d<dimension; ++d)
	  std::cerr << vertexData[cell.wall(0)->vertex1()->index()][d] << " ";
	std::cerr << std::endl << std::endl;
	std::cerr << "2 "; 
	for (size_t d=0; d<dimension; ++d)
	  std::cerr << center[d] << " ";
	std::cerr << std::endl << "2 ";
	for (size_t d=0; d<dimension; ++d)
	  std::cerr << center[d]+normal[d] << " ";
	std::cerr << std::endl << std::endl;
	std::cerr << "3 "; 
	for (size_t d=0; d<dimension; ++d)
	  std::cerr << normal[d] << " ";
	std::cerr << std::endl << std::endl;
	std::cerr << "4 "; 
	for (size_t i=0; i<scalarProdVal.size(); ++i)
	  std::cerr << scalarProdVal[i] << " ";
	std::cerr << std::endl << std::endl;
	std::cerr << "5 "; 
	for (size_t i=0; i<scalarProdSign.size(); ++i)
	  std::cerr << scalarProdSign[i] << " ";
	std::cerr << std::endl << std::endl;
	
	exit(-1);
      }
    }
    // Get the cell size
    double A=1.0;
    if (parameter(1)==1.0)
      A = cell.calculateVolume(vertexData)/cell.numVertex();
    
    double coeff =timeFactor_*parameter(0) * A;
    //update the vertex derivatives
    for (size_t k=0; k<cell.numVertex(); ++k) {
      double vCoeff=coeff;
      if (cell.vertex(k)->isBoundary(T.background()))
	vCoeff *= 1.5;
      for (size_t d=0; d<dimension; ++d) {
	vertexDerivs[cell.vertex(k)->index()][d] += vCoeff * normal[d];
      }
      // std::cerr << timeFactor_ * parameter(0)  << std::endl;
    }	
    // For saving normals in direction used for test plotting
    //for (size_t d=0; d<dimension; ++d)
    //cellData[cell.index()][d] = normal[d];
  }
  //std::cerr << numFlipNormal << " cells out of " << T.numCell() << " has flipped normal."
  //	      << std::endl;
}


void VertexFromCellPlaneLinear::update(Tissue &T,
                                       DataMatrix &cellData,
                                       DataMatrix &wallData,
                                       DataMatrix &vertexData,
                                       double h)
{
  if (timeFactor_ < 1.0 ) {
    timeFactor_ += h/parameter(numParameter()-1);
  }
  if (timeFactor_ >1.0)
    timeFactor_=1.0;
  //cellData[0][12]=timeFactor_*parameter(0);
}



VertexFromCellPlaneLinearCenterTriangulation::
VertexFromCellPlaneLinearCenterTriangulation(std::vector<double> &paraValue,
		    std::vector< std::vector<size_t> > &indValue)
{
  if (paraValue.size() != 3) {
    std::cerr << "VertexFromCellPlaneLinearCenterTriangulation::VertexFromCellPlaneLinearCenterTriangulation() " 
	      << "Uses two parameters: k_force and areaFlag and time span" << std::endl;
    exit(EXIT_FAILURE);
  }
  if (paraValue[1]!=0.0 && paraValue[1]!=1.0) {
    std::cerr << "VertexFromCellPlaneLinearCenterTriangulation::VertexFromCellPlaneLinearCenterTriangulation() " 
	      << "areaFlag must be zero (no area included) or one (area included)." << std::endl;
    exit(EXIT_FAILURE);
  }
  if (paraValue[1]<=0.0 ) {
    std::cerr << "VertexFromCellPlaneLinearCenterTriangulation::VertexFromCellPlaneLinearCenterTriangulation() " 
	      << "time span must be positive." << std::endl;
    exit(EXIT_FAILURE);
  }
  if (indValue.size() != 1) {
    std::cerr << "VertexFromCellPlaneLinearCenterTriangulation::VertexFromCellPlane() " 
	      << std::endl
	       << "Start of additional Cell variable indices (center(x,y,z) " << std::endl;
    exit(EXIT_FAILURE);
  }
  
  setId("VertexFromCellPlaneLinearCenterTriangulation");
  setParameter(paraValue);
  setVariableIndex(indValue);
  
  std::vector<std::string> tmp(numParameter());
  tmp[0] = "k_force";
  tmp[1] = "areaFlag";
  tmp[2] = "deltaT";

  timeFactor_=0.0;
  setParameterId(tmp);
}

void VertexFromCellPlaneLinearCenterTriangulation::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs)
{
  size_t dimension = vertexData[0].size();
  if (dimension!=3) {
    std::cerr << " VertexFromCellPlaneLinearCenterTriangulation::VertexFromCellPlaneLinearCenterTriangulation() " 
	      << " Only implemented for three dimensions ." << std::endl;
    exit(EXIT_FAILURE);
  }
  
 size_t numCells = T.numCell();
 size_t comIndex = variableIndex(0,0);
 size_t lengthInternalIndex = comIndex+dimension;

  for (size_t cellIndex= 0; cellIndex< numCells; ++cellIndex) {  

    size_t numWalls = T.cell(cellIndex).numWall();
    //Cell cell = T.cell(n);        
    // double normal[3]={0,0,0};
    // normal[0]=cellData[n][normalVectorIndex  ]; 
    // normal[1]=cellData[n][normalVectorIndex+1]; 
    // normal[2]=cellData[n][normalVectorIndex+2]; 

     if(  T.cell(cellIndex).numVertex()!= numWalls ) {
       std::cerr << "VertexFromTRBCellPlaneLinearcenterTriangulationMT::derivs() same number of vertices and walls."
		<< " Not for cells with " << T.cell(cellIndex).numWall() << " walls and "
		<< T.cell(cellIndex).numVertex() << " vertices!"	
		<< std::endl;
      exit(-1);
    }
     // One triangle per 'vertex' in cyclic order
     for (size_t wallindex=0; wallindex<numWalls; ++wallindex) { 
       size_t kPlusOneMod = (wallindex+1)%numWalls;
       //size_t v1 = com;
       size_t v2 = T.cell(cellIndex).vertex(wallindex)->index();
       size_t v3 = T.cell(cellIndex).vertex(kPlusOneMod)->index();

       //size_t w1 = internal wallindex
      size_t w2 = T.cell(cellIndex).wall(wallindex)->index();
      //size_t w3 = internal wallindex+1
      // Position matrix holds in rows positions for com, vertex(wallindex), vertex(wallindex+1)
      DataMatrix position(3,vertexData[v2]);
      for (size_t d=0; d<dimension; ++d)
	position[0][d] = cellData[cellIndex][comIndex+d]; // com position
      //position[1] = vertexData[v2]; // given by initiation
      position[2] = vertexData[v3];
      //position[0][2] z for vertex 1 of the current element

      // Lengths are from com-vertex(wallindex), vertex(wallindex)-vertex(wallindex+1) (wall(wallindex)), com-vertex(wallindex+1)
      std::vector<double> length(numWalls);
      length[0] = std::sqrt( (position[0][0]-position[1][0])*(position[0][0]-position[1][0]) +
			     (position[0][1]-position[1][1])*(position[0][1]-position[1][1]) +
			     (position[0][2]-position[1][2])*(position[0][2]-position[1][2]) );
      
      length[1] = T.wall(w2).lengthFromVertexPosition(vertexData);
        
      length[2] = std::sqrt( (position[0][0]-position[2][0])*(position[0][0]-position[2][0]) +
			     (position[0][1]-position[2][1])*(position[0][1]-position[2][1]) +
			     (position[0][2]-position[2][2])*(position[0][2]-position[2][2]) );

      // current Area of the element (using Heron's formula)                                      
      double Area=std::sqrt( ( length[0]+length[1]+length[2])*
			     (-length[0]+length[1]+length[2])*
			     ( length[0]-length[1]+length[2])*
			     ( length[0]+length[1]-length[2])  )*0.25;

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

      double normal[3];      
      normal[0]= Xcurrent[1]*Bcurrent[2]-Xcurrent[2]*Bcurrent[1];
      normal[1]= Xcurrent[2]*Bcurrent[0]-Xcurrent[0]*Bcurrent[2];
      normal[2]= Xcurrent[0]*Bcurrent[1]-Xcurrent[1]*Bcurrent[0];

      tempA=std:: sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
      normal[0]=normal[0]/tempA;
      normal[1]=normal[1]/tempA;
      normal[2]=normal[2]/tempA;

      // normal[0]=cellData[cellIndex][32];  
      // normal[1]=cellData[cellIndex][33];  
      // normal[2]=cellData[cellIndex][34];  

       // Get the cell size
       double A=1.0/3;
       if (parameter(1)==1.0)
	 A = Area/3;
       double coeff =timeFactor_*parameter(0) * A;
       //update the vertex derivatives
                        
       cellDerivs[cellIndex][comIndex  ] +=   coeff * normal[0];  
       cellDerivs[cellIndex][comIndex+1] +=  coeff * normal[1];  
       cellDerivs[cellIndex][comIndex+2] +=  coeff * normal[2];  
       
       vertexDerivs[v2][0] +=  coeff * normal[0];  
       vertexDerivs[v2][1] +=  coeff * normal[1];  
       vertexDerivs[v2][2] +=  coeff * normal[2];  
       
       vertexDerivs[v3][0] +=  coeff * normal[0];  
       vertexDerivs[v3][1] +=  coeff * normal[1];  
       vertexDerivs[v3][2] +=  coeff * normal[2];         
			
     }// for over walls
  }// for over cells  
}




void VertexFromCellPlaneLinearCenterTriangulation::update(Tissue &T,
                                       DataMatrix &cellData,
                                       DataMatrix &wallData,
                                       DataMatrix &vertexData,
                                       double h)
{
  if (timeFactor_ < 1.0 ) {
    timeFactor_ += h/parameter(numParameter()-1);
  }
  if (timeFactor_ >1.0)
    timeFactor_=1.0;
  //cellData[0][12]=timeFactor_*parameter(0);
}



VertexFromCellPlaneSpatial::
VertexFromCellPlaneSpatial(std::vector<double> &paraValue,
			   std::vector< std::vector<size_t> > &indValue)
{
  if (paraValue.size() != 5) {
    std::cerr << "VertexFromCellPlaneSpatial::VertexFromCellPlaneSpatial() " 
	      << "Uses five parameters: k_min, k_max, K_spatial, n_spatial and areaFlag." 
	      << std::endl;
    exit(EXIT_FAILURE);
  }
  if (paraValue[4]!=0.0 && paraValue[4]!=1.0) {
    std::cerr << "VertexFromCellPlaneSpatial::VertexFromCellPlaneSpatial() " 
	      << "areaFlag must be zero (no area included) or one (area included)." << std::endl;
    exit(EXIT_FAILURE);
  }
  if (indValue.size() != 1 || indValue[0].size() != 1) {
    std::cerr << "VertexFromCellPlaneSpatial::VertexFromCellPlaneSpatial() " 
	      << std::endl
	      << "Variable index for spatial max index is used." << std::endl;
    exit(EXIT_FAILURE);
  }
  
  setId("VertexFromCellPlaneSpatial");
  setParameter(paraValue);
  setVariableIndex(indValue);
  
  std::vector<std::string> tmp(numParameter());
  tmp[0] = "k_min";
  tmp[1] = "k_max";
  tmp[2] = "K_spatial";
  tmp[3] = "n_spatial";
  tmp[4] = "areaFlag";
  
  setParameterId(tmp);
  Kpow_= std::pow(paraValue[2],paraValue[3]);
  
}

void VertexFromCellPlaneSpatial::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs)
{
  size_t dimension = vertexData[0].size();
  if (dimension!=3) {
    std::cerr << "VertexFromCellPlaneSpatial::VertexFromCellPlaneSpatial() " 
	      << "Only implemented for three dimensions." << std::endl;
    exit(EXIT_FAILURE);
  }
  size_t sI=variableIndex(0,0);
  assert (sI<vertexData[0].size());
  double max = vertexData[0][sI];
  size_t maxI=0;
  size_t numVertices = vertexData.size();
  for (size_t i=1; i<numVertices; ++i)
    if (vertexData[i][sI]>max) {
      max=vertexData[i][sI];
      maxI=i;
    }
  std::vector<double> maxPos(dimension);
  for (size_t d=0; d<dimension; ++d)
    maxPos[d] = vertexData[maxI][d];
  
  unsigned int numFlipNormal=0;
  for (size_t n = 0; n < T.numCell(); ++n) {
    Cell cell = T.cell(n);
    unsigned int flipFlag=0;
    std::vector<double> normal = cell.getNormalToPCAPlane();
    double norm=0.0;
    for (size_t d=0; d<dimension; ++d)
      norm += normal[d]*normal[d];
    if (norm != 1.0) {
      norm = std::sqrt(norm);
      assert(norm>0.0);
      double normFac = 1.0/norm; 
      for (size_t d=0; d<dimension; ++d)
	normal[d] *= normFac;
    }
    
    std::vector<int> scalarProdSign(cell.numVertex());
    std::vector<double> scalarProdVal(cell.numVertex());
    double scalarProdSum=0.0;
    for (size_t k=0; k<cell.numVertex(); ++k) {
      size_t k2=(k+1)%cell.numVertex();
      size_t k3=(k+2)%cell.numVertex();
      //Make sure ((v2-v1)x(v3-v2))n has same sign for all cells
      assert(cell.numVertex()>2);
      std::vector<double> nw1(dimension),nw2(dimension);
      for (size_t d=0; d<dimension; ++d) {
	nw1[d] = vertexData[cell.vertex(k2)->index()][d]-vertexData[cell.vertex(k)->index()][d];
	nw2[d] = vertexData[cell.vertex(k3)->index()][d]-vertexData[cell.vertex(k2)->index()][d];
      }
      //cross product
      double scalarProd=0.0;
      for (size_t d1=0; d1<dimension; ++d1) {
	size_t d2=(d1+1)%dimension;
	size_t d3=(d1+2)%dimension;
	scalarProd += (nw1[d1]*nw2[d2]-nw1[d2]*nw2[d1])*normal[d3];
      }
      scalarProdVal[k] = scalarProd;
      scalarProdSum += scalarProd;
      if (scalarProd>0.0)
	scalarProdSign[k]=1;
      else
	scalarProdSign[k]=-1;
    }
    //  		for (size_t k=1; k<cell.numVertex(); ++k)
    //  			if (scalarProdSign[k]!=scalarProdSign[0]) {
    //  				std::cerr << "Cell " << n << " has diverging signs on scalar product." << std::endl;
    //  				break;
    //  			}
    int scalarProdSignSum=0;
    for (size_t k=0; k<scalarProdSign.size(); ++k)
      scalarProdSignSum += scalarProdSign[k];
    
    if (scalarProdSignSum<0) {
      numFlipNormal++;
      flipFlag=1;
      for (size_t d=0; d<dimension; ++d)
	normal[d] = -normal[d];
    }
    else if (scalarProdSignSum==0) {
      std::cerr << "Cell " << n << " has no majority sign in right hand rule expression." 
		<< std::endl;
      if (std::fabs(scalarProdSum)>0.01) {
	if (scalarProdSum<0.0) {
	  numFlipNormal++;
	  flipFlag=1;
	  for (size_t d=0; d<dimension; ++d)
	    normal[d] = -normal[d];
	}
      }
      else {
	std::vector<double> center = cell.positionFromVertex(vertexData);
	// Print all walls
	for (size_t k=0; k<cell.numWall(); ++k) {
	  std::cerr << "0 "; 
	  for (size_t d=0; d<dimension; ++d)
	    std::cerr << vertexData[cell.wall(k)->vertex1()->index()][d] << " ";
	  std::cerr << std::endl << "0 "; 
	  for (size_t d=0; d<dimension; ++d)
	    std::cerr << vertexData[cell.wall(k)->vertex2()->index()][d] << " ";
	  std::cerr << std::endl << std::endl;
	}		 
	std::cerr << "1 "; 
	for (size_t d=0; d<dimension; ++d)
	  std::cerr << vertexData[cell.wall(0)->vertex1()->index()][d] << " ";
	std::cerr << std::endl << std::endl;
	std::cerr << "2 "; 
	for (size_t d=0; d<dimension; ++d)
	  std::cerr << center[d] << " ";
	std::cerr << std::endl << "2 ";
	for (size_t d=0; d<dimension; ++d)
	  std::cerr << center[d]+normal[d] << " ";
	std::cerr << std::endl << std::endl;
	std::cerr << "3 "; 
	for (size_t d=0; d<dimension; ++d)
	  std::cerr << normal[d] << " ";
	std::cerr << std::endl << std::endl;
	std::cerr << "4 "; 
	for (size_t i=0; i<scalarProdVal.size(); ++i)
	  std::cerr << scalarProdVal[i] << " ";
	std::cerr << std::endl << std::endl;
	std::cerr << "5 "; 
	for (size_t i=0; i<scalarProdSign.size(); ++i)
	  std::cerr << scalarProdSign[i] << " ";
	std::cerr << std::endl << std::endl;
	
	exit(EXIT_FAILURE);
      }
    }
    // Get the cell size
    double A=1.0;
    if (parameter(4)==1.0)
      A = cell.calculateVolume(vertexData)/cell.numVertex();
    // Calculate the spatial factor
    
    double maxDistance=0.0;
    std::vector<double> cellPos = cell.positionFromVertex(vertexData);
    for (size_t d=0; d<dimension; ++d)
      maxDistance += (maxPos[d]-cellPos[d])*(maxPos[d]-cellPos[d]);
    maxDistance = std::sqrt(maxDistance);
    double sFactor = std::pow(maxDistance,parameter(3));
    sFactor = parameter(1)*Kpow_/(Kpow_+sFactor);
    
    double coeff = (parameter(0)+sFactor) * A;
    
    
    //update the vertex derivatives
    for (size_t k=0; k<cell.numVertex(); ++k) {
      double vCoeff=coeff;
      if (cell.vertex(k)->isBoundary(T.background()))
	vCoeff *= 1.5;
      for (size_t d=0; d<dimension; ++d) {
	vertexDerivs[cell.vertex(k)->index()][d] += vCoeff * normal[d];
      }
    }	
  }
  //std::cerr << numFlipNormal << " cells out of " << T.numCell() << " has flipped normal."
  //				<< std::endl;
}

VertexFromCellPlaneConcentrationHill::
VertexFromCellPlaneConcentrationHill(std::vector<double> &paraValue,
				     std::vector< std::vector<size_t> > &indValue)
{
  if (paraValue.size() != 5) {
    std::cerr << "VertexFromCellPlaneConcentrationHill::"
	      << "VertexFromCellPlaneConcentrationHill() " 
	      << "Uses five parameters: k_min, k_max, K_H, n_H and areaFlag." 
	      << std::endl;
    exit(EXIT_FAILURE);
  }
  if (paraValue[4]!=0.0 && paraValue[4]!=1.0) {
    std::cerr << "VertexFromCellPlaneConcentrationHill::"
	      << "VertexFromCellPlaneConcentrationHill() " 
	      << "areaFlag must be zero (no area included) or one (area included)." << std::endl;
    exit(EXIT_FAILURE);
  }
  if (indValue.size() != 1 || indValue[0].size() != 1) {
    std::cerr << "VertexFromCellPlaneSpatial::VertexFromCellPlaneSpatial() " 
	      << std::endl
	      << "Variable index for concentration needed." << std::endl;
    exit(EXIT_FAILURE);
  }
  
  setId("VertexFromCellPlaneConcentrationHill");
  setParameter(paraValue);
  setVariableIndex(indValue);
  
  std::vector<std::string> tmp(numParameter());
  tmp[0] = "k_min";
  tmp[1] = "k_max";
  tmp[2] = "K_H";
  tmp[3] = "n_H";
  tmp[4] = "areaFlag";
  
  setParameterId(tmp);
  Kpow_= std::pow(paraValue[2],paraValue[3]);
}

void VertexFromCellPlaneConcentrationHill::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs)
{
  size_t dimension = vertexData[0].size();
  if (dimension!=3) {
    std::cerr << "VertexFromCellPlaneConcentrationHill::"
	      << "VertexFromCellPlaneConcentrationHill() " 
	      << "Only implemented for three dimensions." << std::endl;
    exit(EXIT_FAILURE);
  }
  size_t cI=variableIndex(0,0);
  assert (cI<cellData[0].size());

  unsigned int numFlipNormal=0;
  for (size_t n = 0; n < T.numCell(); ++n) {
    Cell cell = T.cell(n);
    unsigned int flipFlag=0;
    std::vector<double> normal = cell.getNormalToPCAPlane();
    double norm=0.0;
    for (size_t d=0; d<dimension; ++d)
      norm += normal[d]*normal[d];
    if (norm != 1.0) {
      norm = std::sqrt(norm);
      assert(norm>0.0);
      double normFac = 1.0/norm; 
      for (size_t d=0; d<dimension; ++d)
	normal[d] *= normFac;
    }
    
    std::vector<int> scalarProdSign(cell.numVertex());
    std::vector<double> scalarProdVal(cell.numVertex());
    double scalarProdSum=0.0;
    for (size_t k=0; k<cell.numVertex(); ++k) {
      size_t k2=(k+1)%cell.numVertex();
      size_t k3=(k+2)%cell.numVertex();
      //Make sure ((v2-v1)x(v3-v2))n has same sign for all cells
      assert(cell.numVertex()>2);
      std::vector<double> nw1(dimension),nw2(dimension);
      for (size_t d=0; d<dimension; ++d) {
	nw1[d] = vertexData[cell.vertex(k2)->index()][d]-vertexData[cell.vertex(k)->index()][d];
	nw2[d] = vertexData[cell.vertex(k3)->index()][d]-vertexData[cell.vertex(k2)->index()][d];
      }
      //cross product
      double scalarProd=0.0;
      for (size_t d1=0; d1<dimension; ++d1) {
	size_t d2=(d1+1)%dimension;
	size_t d3=(d1+2)%dimension;
	scalarProd += (nw1[d1]*nw2[d2]-nw1[d2]*nw2[d1])*normal[d3];
      }
      scalarProdVal[k] = scalarProd;
      scalarProdSum += scalarProd;
      if (scalarProd>0.0)
	scalarProdSign[k]=1;
      else
	scalarProdSign[k]=-1;
    }
    // for (size_t k=1; k<cell.numVertex(); ++k)
    //   if (scalarProdSign[k]!=scalarProdSign[0]) {
    //     std::cerr << "Cell " << n << " has diverging signs on scalar product." << std::endl;
    //     break;
    //   }
    int scalarProdSignSum=0;
    for (size_t k=0; k<scalarProdSign.size(); ++k)
      scalarProdSignSum += scalarProdSign[k];
    
    if (scalarProdSignSum<0) {
      numFlipNormal++;
      flipFlag=1;
      for (size_t d=0; d<dimension; ++d)
	normal[d] = -normal[d];
    }
    else if (scalarProdSignSum==0) {
      std::cerr << "Cell " << n << " has no majority sign in right hand rule expression." 
		<< std::endl;
      if (std::fabs(scalarProdSum)>0.01) {
	if (scalarProdSum<0.0) {
	  numFlipNormal++;
	  flipFlag=1;
	  for (size_t d=0; d<dimension; ++d)
	    normal[d] = -normal[d];
	}
      }
      else {
	std::vector<double> center = cell.positionFromVertex(vertexData);
	// Print all walls
	for (size_t k=0; k<cell.numWall(); ++k) {
	  std::cerr << "0 "; 
	  for (size_t d=0; d<dimension; ++d)
	    std::cerr << vertexData[cell.wall(k)->vertex1()->index()][d] << " ";
	  std::cerr << std::endl << "0 "; 
	  for (size_t d=0; d<dimension; ++d)
	    std::cerr << vertexData[cell.wall(k)->vertex2()->index()][d] << " ";
	  std::cerr << std::endl << std::endl;
	}		 
	std::cerr << "1 "; 
	for (size_t d=0; d<dimension; ++d)
	  std::cerr << vertexData[cell.wall(0)->vertex1()->index()][d] << " ";
	std::cerr << std::endl << std::endl;
	std::cerr << "2 "; 
	for (size_t d=0; d<dimension; ++d)
	  std::cerr << center[d] << " ";
	std::cerr << std::endl << "2 ";
	for (size_t d=0; d<dimension; ++d)
	  std::cerr << center[d]+normal[d] << " ";
	std::cerr << std::endl << std::endl;
	std::cerr << "3 "; 
	for (size_t d=0; d<dimension; ++d)
	  std::cerr << normal[d] << " ";
	std::cerr << std::endl << std::endl;
	std::cerr << "4 "; 
	for (size_t i=0; i<scalarProdVal.size(); ++i)
	  std::cerr << scalarProdVal[i] << " ";
	std::cerr << std::endl << std::endl;
	std::cerr << "5 "; 
	for (size_t i=0; i<scalarProdSign.size(); ++i)
	  std::cerr << scalarProdSign[i] << " ";
	std::cerr << std::endl << std::endl;
	exit(EXIT_FAILURE);
      }
    }
    // Get the cell size
    double A=1.0;
    if (parameter(4)==1.0)
      A = cell.calculateVolume(vertexData)/cell.numVertex();
    // Calculate the concentration dependent factor    
    double cFactor = std::pow(cellData[n][cI],parameter(3));
    cFactor = parameter(1)*cFactor/(Kpow_+cFactor);    
    double coeff = (parameter(0)+cFactor) * A;
    
    //update the vertex derivatives
    for (size_t k=0; k<cell.numVertex(); ++k) {
      double vCoeff=coeff;
      if (cell.vertex(k)->isBoundary(T.background()))
	vCoeff *= 1.5;
      for (size_t d=0; d<dimension; ++d) {
	vertexDerivs[cell.vertex(k)->index()][d] += vCoeff * normal[d];
      }
    }	
  }
  //std::cerr << numFlipNormal << " cells out of " << T.numCell() << " has flipped normal."
  //				<< std::endl;
}

VertexFromCellPlaneNormalized::
VertexFromCellPlaneNormalized(std::vector<double> &paraValue,
			      std::vector< std::vector<size_t> > &indValue)
{
  if (paraValue.size() != 2) {
    std::cerr << "VertexFromCellPlaneNormalized::VertexFromCellPlaneNormalized() " 
	      << "Uses two parameters: k_force and areaFlag" << std::endl;
    exit(EXIT_FAILURE);
  }
  if (paraValue[1]!=0.0 && paraValue[1]!=1.0) {
    std::cerr << "VertexFromCellPlaneNormalized::VertexFromCellPlaneNormalized() " 
	      << "areaFlag (p1) needs to be zero or one (1=area factor included)." << std::endl;
    exit(EXIT_FAILURE);
  }
  
  if (indValue.size() != 0) {
    std::cerr << "VertexFromCellPlaneNormalized::VertexFromCellPlaneNormalized() " 
	      << std::endl
	      << "No variable index used." << std::endl;
    exit(EXIT_FAILURE);
  }
  
  setId("VertexFromCellPlaneNormalized");
  setParameter(paraValue);
  setVariableIndex(indValue);
  
  std::vector<std::string> tmp(numParameter());
  tmp[0] = "k_force";
  tmp[1] = "areaFlag";
  
  setParameterId(tmp);
}

void VertexFromCellPlaneNormalized::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs)
{
  size_t dimension = vertexData[0].size();
  if (dimension!=3) {
    std::cerr << "VertexFromCellPlaneNormalized::VertexFromCellPlaneNormalized() " 
	      << "Only implemented for three dimensions." << std::endl;
    exit(EXIT_FAILURE);
  }
  std::vector<double> tmpD(dimension);
  DataMatrix vertexDerivsTmp(vertexDerivs.size(),tmpD);
  unsigned int numFlipNormal=0;
  for (size_t n = 0; n < T.numCell(); ++n) {
    Cell cell = T.cell(n);
    // This calculation should now be done in reaction CalculatePCAPlane
    //cell.calculatePCAPlane(vertexData);
    unsigned int flipFlag=0;
    
    std::vector<double> normal = cell.getNormalToPCAPlane();
    double norm=0.0;
    for (size_t d=0; d<dimension; ++d)
      norm += normal[d]*normal[d];
    if (norm != 1.0) {
      norm = std::sqrt(norm);
      assert(norm>0.0);
      double normFac = 1.0/norm; 
      for (size_t d=0; d<dimension; ++d)
	normal[d] *= normFac;
    }
    
    std::vector<int> scalarProdSign(cell.numVertex());
    std::vector<double> scalarProdVal(cell.numVertex());
    double scalarProdSum=0.0;
    for (size_t k=0; k<cell.numVertex(); ++k) {
      size_t k2=(k+1)%cell.numVertex();
      size_t k3=(k+2)%cell.numVertex();
      //Make sure ((v2-v1)x(v3-v2))n has same sign for all cells
      assert(cell.numVertex()>2);
      std::vector<double> nw1(dimension),nw2(dimension);
      for (size_t d=0; d<dimension; ++d) {
	nw1[d] = vertexData[cell.vertex(k2)->index()][d]-vertexData[cell.vertex(k)->index()][d];
	nw2[d] = vertexData[cell.vertex(k3)->index()][d]-vertexData[cell.vertex(k2)->index()][d];
      }
      //cross product
      double scalarProd=0.0;
      for (size_t d1=0; d1<dimension; ++d1) {
	size_t d2=(d1+1)%dimension;
	size_t d3=(d1+2)%dimension;
	scalarProd += (nw1[d1]*nw2[d2]-nw1[d2]*nw2[d1])*normal[d3];
      }
      scalarProdVal[k] = scalarProd;
      scalarProdSum += scalarProd;
      if (scalarProd>0.0)
	scalarProdSign[k]=1;
      else
	scalarProdSign[k]=-1;
    }
    // for (size_t k=1; k<cell.numVertex(); ++k)
    //   if (scalarProdSign[k]!=scalarProdSign[0]) {
    //     std::cerr << "Cell " << n << " has diverging signs on scalar product." << std::endl;
    //     break;
    //   }
    int scalarProdSignSum=0;
    for (size_t k=0; k<scalarProdSign.size(); ++k)
      scalarProdSignSum += scalarProdSign[k];
    
    if (scalarProdSignSum<0) {
      numFlipNormal++;
      flipFlag=1;
      for (size_t d=0; d<dimension; ++d)
	normal[d] = -normal[d];
    }
    else if (scalarProdSignSum==0) {
      //std::cerr << "Cell " << n << " has no majority sign in right hand rule expression." 
      //				<< std::endl;
      if (std::fabs(scalarProdSum)>0.01) {
	if (scalarProdSum<0.0) {
	  numFlipNormal++;
	  flipFlag=1;
	  for (size_t d=0; d<dimension; ++d)
	    normal[d] = -normal[d];
	}
      }
      else {
	std::vector<double> center = cell.positionFromVertex(vertexData);
	// Print all walls
	for (size_t k=0; k<cell.numWall(); ++k) {
	  std::cerr << "0 "; 
	  for (size_t d=0; d<dimension; ++d)
	    std::cerr << vertexData[cell.wall(k)->vertex1()->index()][d] << " ";
	  std::cerr << std::endl << "0 "; 
	  for (size_t d=0; d<dimension; ++d)
	    std::cerr << vertexData[cell.wall(k)->vertex2()->index()][d] << " ";
	  std::cerr << std::endl << std::endl;
	}		 
	std::cerr << "1 "; 
	for (size_t d=0; d<dimension; ++d)
	  std::cerr << vertexData[cell.wall(0)->vertex1()->index()][d] << " ";
	std::cerr << std::endl << std::endl;
	std::cerr << "2 "; 
	for (size_t d=0; d<dimension; ++d)
	  std::cerr << center[d] << " ";
	std::cerr << std::endl << "2 ";
	for (size_t d=0; d<dimension; ++d)
	  std::cerr << center[d]+normal[d] << " ";
	std::cerr << std::endl << std::endl;
	std::cerr << "3 "; 
	for (size_t d=0; d<dimension; ++d)
	  std::cerr << normal[d] << " ";
	std::cerr << std::endl << std::endl;
	std::cerr << "4 "; 
	for (size_t i=0; i<scalarProdVal.size(); ++i)
	  std::cerr << scalarProdVal[i] << " ";
	std::cerr << std::endl << std::endl;
	std::cerr << "5 "; 
	for (size_t i=0; i<scalarProdSign.size(); ++i)
	  std::cerr << scalarProdSign[i] << " ";
	std::cerr << std::endl << std::endl;
	
	exit(-1);
      }
    }
    // Get the cell size
    double A=1.0;
    if (parameter(1)!=0.0)
      A = cell.calculateVolume(vertexData)/cell.numVertex();
    
    //update the vertex derivatives
    for (size_t d=0; d<dimension; ++d) {
      //escellData[n][d]=normal[d];
      for (size_t k=0; k<cell.numVertex(); ++k)
	vertexDerivsTmp[cell.vertex(k)->index()][d] += A * normal[d];
    }
  }	
  // Normalize and update all vertices
  size_t nVertex = vertexDerivs.size();
  for (size_t i=0; i<nVertex; ++i) {
    double norm=0.0;
    for (size_t d=0; d<dimension; ++d)
      norm += vertexDerivsTmp[i][d]*vertexDerivsTmp[i][d];
    double normFac = 1.0/std::sqrt(norm);
    for (size_t d=0; d<dimension; ++d) {
      vertexDerivsTmp[i][d] *= normFac;
      vertexDerivs[i][d] += parameter(0) * vertexDerivsTmp[i][d];
    }
  }
  //std::cerr << numFlipNormal << " cells out of " << T.numCell() << " has flipped normal."
  //				<< std::endl;
}

VertexFromCellPlaneNormalizedSpatial::
VertexFromCellPlaneNormalizedSpatial(std::vector<double> &paraValue,
				     std::vector< std::vector<size_t> > &indValue)
{
  if (paraValue.size() != 5) {
    std::cerr << "VertexFromCellPlaneNormalizedSpatial::VertexFromCellPlaneNormalizedSpatial() " 
	      << "Uses four parameters: F_min F_max K_spatial n_spatial areaFlag" << std::endl;
    exit(EXIT_FAILURE);
  }
  if (paraValue[4]!=0.0 && paraValue[4]!=1.0) {
    std::cerr << "VertexFromCellPlaneNormalizedSpatial::VertexFromCellPlaneNormalizedSpatial() " 
	      << "areaFlag needs to be zero or one (1=area factor included)." << std::endl;
    exit(EXIT_FAILURE);
  }
  
  if (indValue.size() != 1 || indValue[0].size() != 1) {
    std::cerr << "VertexFromCellPlaneNormalizedSpatial::VertexFromCellPlaneNormalizedSpatial() " 
	      << std::endl
	      << "Variable index for spatial coordinate used." << std::endl;
    exit(EXIT_FAILURE);
  }
  
  setId("VertexFromCellPlaneNormalizedSpatial");
  setParameter(paraValue);
  setVariableIndex(indValue);
  
  std::vector<std::string> tmp(numParameter());
  tmp[0] = "F_min";
  tmp[1] = "F_max";
  tmp[2] = "K_spatial";
  tmp[3] = "n_spatial";
  tmp[4] = "areaFlag";
  
  setParameterId(tmp);
  Kpow_= std::pow(paraValue[2],paraValue[3]);
}

void VertexFromCellPlaneNormalizedSpatial::
derivs(Tissue &T,
			 DataMatrix &cellData,
			 DataMatrix &wallData,
			 DataMatrix &vertexData,
			 DataMatrix &cellDerivs,
			 DataMatrix &wallDerivs,
			 DataMatrix &vertexDerivs)
{
	size_t dimension = vertexData[0].size();
	if (dimension!=3) {
		std::cerr << "VertexFromCellPlaneNormalizedSpatial::VertexFromCellPlaneNormalizedSpatial() " 
							<< "Only implemented for three dimensions." << std::endl;
		exit(EXIT_FAILURE);
	}
	// Initiate spatial factor
	size_t sI=variableIndex(0,0);
	size_t numVertices=vertexData.size();
	assert (sI<dimension);
	double sMax=vertexData[0][sI];
	double maxI=0;
	for (size_t i=1; i<numVertices; ++i)
		if (vertexData[i][sI]>sMax) {
			sMax=vertexData[i][sI];
			maxI=i;
		}

	std::vector<double> maxPos(dimension);
	for (size_t d=0; d<dimension; ++d)
		maxPos[d] = vertexData[maxI][d];


	std::vector<double> tmpD(dimension);
	DataMatrix vertexDerivsTmp(vertexDerivs.size(),tmpD);
	unsigned int numFlipNormal=0;
	size_t numCells=T.numCell();
	for (size_t n=0; n<numCells; ++n) {
		Cell cell = T.cell(n);
		// This calculation should now be done in reaction CalculatePCAPlane
		//cell.calculatePCAPlane(vertexData);
		unsigned int flipFlag=0;
		
		std::vector<double> normal = cell.getNormalToPCAPlane();
		double norm=0.0;
		for (size_t d=0; d<dimension; ++d)
			norm += normal[d]*normal[d];
		if (norm != 1.0) {
			norm = std::sqrt(norm);
			assert(norm>0.0);
			double normFac = 1.0/norm; 
			for (size_t d=0; d<dimension; ++d)
				normal[d] *= normFac;
		}
		
		std::vector<int> scalarProdSign(cell.numVertex());
		std::vector<double> scalarProdVal(cell.numVertex());
		double scalarProdSum=0.0;
		for (size_t k=0; k<cell.numVertex(); ++k) {
			size_t k2=(k+1)%cell.numVertex();
			size_t k3=(k+2)%cell.numVertex();
			//Make sure ((v2-v1)x(v3-v2))n has same sign for all cells
			assert(cell.numVertex()>2);
			std::vector<double> nw1(dimension),nw2(dimension);
			for (size_t d=0; d<dimension; ++d) {
				nw1[d] = vertexData[cell.vertex(k2)->index()][d]-vertexData[cell.vertex(k)->index()][d];
				nw2[d] = vertexData[cell.vertex(k3)->index()][d]-vertexData[cell.vertex(k2)->index()][d];
			}
			//cross product
			double scalarProd=0.0;
			for (size_t d1=0; d1<dimension; ++d1) {
				size_t d2=(d1+1)%dimension;
				size_t d3=(d1+2)%dimension;
				scalarProd += (nw1[d1]*nw2[d2]-nw1[d2]*nw2[d1])*normal[d3];
			}
			scalarProdVal[k] = scalarProd;
			scalarProdSum += scalarProd;
			if (scalarProd>0.0)
				scalarProdSign[k]=1;
			else
				scalarProdSign[k]=-1;
		}
//  		for (size_t k=1; k<cell.numVertex(); ++k)
//  			if (scalarProdSign[k]!=scalarProdSign[0]) {
//  				std::cerr << "Cell " << n << " has diverging signs on scalar product." << std::endl;
//  				break;
//  			}
		int scalarProdSignSum=0;
		for (size_t k=0; k<scalarProdSign.size(); ++k)
			scalarProdSignSum += scalarProdSign[k];
		
		if (scalarProdSignSum<0) {
			numFlipNormal++;
			flipFlag=1;
			for (size_t d=0; d<dimension; ++d)
				normal[d] = -normal[d];
		}
		else if (scalarProdSignSum==0) {
			//std::cerr << "Cell " << n << " has no majority sign in right hand rule expression." 
			//				<< std::endl;
			if (std::fabs(scalarProdSum)>0.01) {
				if (scalarProdSum<0.0) {
					numFlipNormal++;
					flipFlag=1;
					for (size_t d=0; d<dimension; ++d)
						normal[d] = -normal[d];
				}
			}
			else {
				std::vector<double> center = cell.positionFromVertex(vertexData);
				// Print all walls
				for (size_t k=0; k<cell.numWall(); ++k) {
					std::cerr << "0 "; 
					for (size_t d=0; d<dimension; ++d)
						std::cerr << vertexData[cell.wall(k)->vertex1()->index()][d] << " ";
					std::cerr << std::endl << "0 "; 
					for (size_t d=0; d<dimension; ++d)
						std::cerr << vertexData[cell.wall(k)->vertex2()->index()][d] << " ";
					std::cerr << std::endl << std::endl;
				}		 
				std::cerr << "1 "; 
				for (size_t d=0; d<dimension; ++d)
					std::cerr << vertexData[cell.wall(0)->vertex1()->index()][d] << " ";
				std::cerr << std::endl << std::endl;
				std::cerr << "2 "; 
				for (size_t d=0; d<dimension; ++d)
					std::cerr << center[d] << " ";
				std::cerr << std::endl << "2 ";
				for (size_t d=0; d<dimension; ++d)
					std::cerr << center[d]+normal[d] << " ";
				std::cerr << std::endl << std::endl;
				std::cerr << "3 "; 
				for (size_t d=0; d<dimension; ++d)
					std::cerr << normal[d] << " ";
				std::cerr << std::endl << std::endl;
				std::cerr << "4 "; 
				for (size_t i=0; i<scalarProdVal.size(); ++i)
					std::cerr << scalarProdVal[i] << " ";
				std::cerr << std::endl << std::endl;
				std::cerr << "5 "; 
				for (size_t i=0; i<scalarProdSign.size(); ++i)
					std::cerr << scalarProdSign[i] << " ";
				std::cerr << std::endl << std::endl;
				
				exit(-1);
			}
		}
		// Get the cell size
		double A=1.0;
		if (parameter(4)!=0.0)
			A = cell.calculateVolume(vertexData)/cell.numVertex();
			
		//update the vertex derivatives
		for (size_t d=0; d<dimension; ++d) {
			//escellData[n][d]=normal[d];
			for (size_t k=0; k<cell.numVertex(); ++k)
				vertexDerivsTmp[cell.vertex(k)->index()][d] += A * normal[d];
		}
	}	
	// Normalize and update all vertices
	size_t nVertex = vertexDerivs.size();
	for (size_t i=0; i<nVertex; ++i) {
		double norm=0.0;
		for (size_t d=0; d<dimension; ++d)
			norm += vertexDerivsTmp[i][d]*vertexDerivsTmp[i][d];
		double normFac = 1.0/std::sqrt(norm);
		// Calculate spatial factor
		double maxDistance=0.0;
		for (size_t d=0; d<dimension; ++d)
			maxDistance += (maxPos[d]-vertexData[i][d])*(maxPos[d]-vertexData[i][d]);
		maxDistance = std::sqrt(maxDistance);
		double sFactor = std::pow(maxDistance,parameter(3));
		sFactor = parameter(1)*Kpow_/(Kpow_+sFactor);
		for (size_t d=0; d<dimension; ++d) {
			vertexDerivsTmp[i][d] *= normFac;
			vertexDerivs[i][d] += (parameter(0)+sFactor) * vertexDerivsTmp[i][d];
		}
	}
	//std::cerr << numFlipNormal << " cells out of " << T.numCell() << " has flipped normal."
	//				<< std::endl;
}

VertexFromCellPlaneSphereCylinder::
VertexFromCellPlaneSphereCylinder(std::vector<double> &paraValue,
																	std::vector< std::vector<size_t> > &indValue)
{
	if (paraValue.size() != 2) {
		std::cerr << "VertexFromCellPlaneSphereCylinder::VertexFromCellPlaneSphereCylinder() " 
							<< "Uses two parameters: k_force and areaFlag" << std::endl;
		exit(EXIT_FAILURE);
	}
	if (paraValue[1]!=0.0 && paraValue[1]!=1.0) {
		std::cerr << "VertexFromCellPlaneSphereCylinder::VertexFromCellPlaneSphereCylinder() " 
							<< "areaFlag (p1) needs to be zero or one (1=area factor included)." << std::endl;
		exit(EXIT_FAILURE);
	}
	
	if (indValue.size() != 0) {
		std::cerr << "VertexFromCellPlaneSphereCylinder::VertexFromCellPlaneSphereCylinder() " 
							<< std::endl
							<< "No variable index used." << std::endl;
		exit(EXIT_FAILURE);
	}
	
	setId("VertexFromCellPlaneSphereCylinder");
	setParameter(paraValue);
	setVariableIndex(indValue);
	
	std::vector<std::string> tmp(numParameter());
	tmp[0] = "k_force";
	tmp[1] = "areaFlag";
	
	setParameterId(tmp);
}

void VertexFromCellPlaneSphereCylinder::
derivs(Tissue &T,
			 DataMatrix &cellData,
			 DataMatrix &wallData,
			 DataMatrix &vertexData,
			 DataMatrix &cellDerivs,
			 DataMatrix &wallDerivs,
			 DataMatrix &vertexDerivs)
{
	size_t dimension = vertexData[0].size();
	if (dimension!=3) {
		std::cerr << "VertexFromCellPlaneSphereCylinder::VertexFromCellPlaneSphereCylinder() " 
							<< "Only implemented for three dimensions." << std::endl;
		exit(EXIT_FAILURE);
	}
	for (size_t n = 0; n < T.numCell(); ++n) {
		Cell cell = T.cell(n);
		// This calculation should now be done in CalculatePCAPlane reaction.
		//cell.calculatePCAPlane(vertexData);
		std::vector<double> normal = cell.getNormalToPCAPlane();
		double norm=0.0;
		for (size_t d=0; d<dimension; ++d)
			norm += normal[d]*normal[d];
		if (norm != 1.0) {
			norm = std::sqrt(norm);
			assert(norm>0.0);
			for (size_t d=0; d<dimension; ++d)
				normal[d] /= norm;
		}
		std::vector<double> cellPos = cell.positionFromVertex(vertexData);

		std::vector<double> scNorm(dimension);
		for (size_t d=0; d<dimension; ++d)
			scNorm[d] = cellPos[d];
		if( cellPos[2]<0.0 )
			scNorm[2]=0.0;
		double scalarProd=0.0;
		for (size_t d=0; d<dimension; ++d)
			scalarProd += scNorm[d]*normal[d];
		if (scalarProd<0.0) {
			for (size_t d=0; d<dimension; ++d)
				normal[d] = -normal[d];
		}
		// Introduce cell size
		
		double A = 1.0;
		if (parameter(1)!=0.0)
			A = cell.calculateVolume(vertexData)/cell.numVertex(); 

		double coeff = parameter(0) * A; 
		for (size_t k=0; k<cell.numVertex(); ++k) {
			double vCoeff=coeff;
			if (cell.vertex(k)->isBoundary(T.background()))
				vCoeff *= 1.5;
			for (size_t d=0; d<dimension; ++d)
				vertexDerivs[cell.vertex(k)->index()][d] += vCoeff * normal[d];
		}
	}	
}

VertexFromCellPlaneSphereCylinderConcentrationHill::
VertexFromCellPlaneSphereCylinderConcentrationHill(std::vector<double> &paraValue,
						   std::vector< std::vector<size_t> > &indValue)
{
  if (paraValue.size() != 5) {
    std::cerr << "VertexFromCellPlaneSphereCylinderConcentrationHill::"
	      << "VertexFromCellPlaneSphereCylinderConcentrationHill() " 
	      << "Uses five parameters: k_forceConst, k_forceHill, K_Hill, n_Hill areaFlag" 
	      << std::endl;
    exit(EXIT_FAILURE);
  }
  if (paraValue[4]!=0.0 && paraValue[4]!=1.0) {
    std::cerr << "VertexFromCellPlaneSphereCylinderConcentrationHill::"
	      << "VertexFromCellPlaneSphereCylinderConcentrationHill() " 
	      << "areaFlag (p4) needs to be zero or one (1=area factor included)." 
	      << std::endl;
    exit(EXIT_FAILURE);
  }
  
  if (indValue.size() != 1 || indValue[0].size() != 1) {
    std::cerr << "VertexFromCellPlaneSphereCylinderConcentrationHill::"
	      << "VertexFromCellPlaneSphereCylinderConcentrationHill() " 
	      << std::endl
	      << "One variable index for concentration index is used." << std::endl;
    exit(EXIT_FAILURE);
  }
  
  setId("VertexFromCellPlaneSphereCylinderConcentrationHill");
  setParameter(paraValue);
  setVariableIndex(indValue);
  
  std::vector<std::string> tmp(numParameter());
  tmp[0] = "k_forceConst";
  tmp[1] = "k_forceHill";
  tmp[2] = "K_Hill";
  tmp[3] = "n_Hill";
  tmp[4] = "areaFlag";
  
  setParameterId(tmp);
}

void VertexFromCellPlaneSphereCylinderConcentrationHill::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs)
{
  size_t dimension = vertexData[0].size();
  if (dimension!=3) {
    std::cerr << "VertexFromCellPlaneSphereCylinderConcentrationHill::"
	      << "VertexFromCellPlaneSphereCylinderConcentrationHill() " 
	      << "Only implemented for three dimensions." << std::endl;
    exit(EXIT_FAILURE);
  }
  double Kpow = std::pow(parameter(2),parameter(3));
  for (size_t n = 0; n < T.numCell(); ++n) {
    Cell cell = T.cell(n);
    // This calculation should now be done in reaction CalculatePCAPlane
    //cell.calculatePCAPlane(vertexData);
    double concpow = std::pow(cellData[cell.index()][variableIndex(0,0)],parameter(3));
    std::vector<double> normal = cell.getNormalToPCAPlane();
    double norm=0.0;
    for (size_t d=0; d<dimension; ++d)
      norm += normal[d]*normal[d];
    if (norm != 1.0) {
      norm = std::sqrt(norm);
      assert(norm>0.0);
      for (size_t d=0; d<dimension; ++d)
	normal[d] /= norm;
    }
    std::vector<double> cellPos = cell.positionFromVertex(vertexData);
    
    std::vector<double> scNorm(dimension);
    for (size_t d=0; d<dimension; ++d)
      scNorm[d] = cellPos[d];
    if( cellPos[2]<0.0 )
      scNorm[2]=0.0;
    double scalarProd=0.0;
    for (size_t d=0; d<dimension; ++d)
      scalarProd += scNorm[d]*normal[d];
    if (scalarProd<0.0) {
      for (size_t d=0; d<dimension; ++d)
	normal[d] = -normal[d];
    }
    double coeff = parameter(0) + parameter(1)*concpow/(concpow+Kpow);
    double A=1.0;
    if (parameter(4)!=0.0) {
      A = cell.calculateVolume(vertexData)/cell.numVertex();
      coeff *= A;
    }
    for (size_t k=0; k<cell.numVertex(); ++k) {
      double vCoeff=coeff;
      if (cell.vertex(k)->isBoundary(T.background()))
	vCoeff *= 1.5;
      for (size_t d=0; d<dimension; ++d)
	vertexDerivs[cell.vertex(k)->index()][d] += vCoeff * normal[d];
    }
  }
}


VertexFromCellPlaneTriangular::
VertexFromCellPlaneTriangular(std::vector<double> &paraValue,
			      std::vector< std::vector<size_t> > &indValue)
{
  if (paraValue.size() != 3 && paraValue.size() !=6) {
    std::cerr << "VertexFromCellPlaneTriangular::VertexFromCellPlaneTriangular()" 
	      << " Uses three parameters: k_force, areaFlag(0/1 without/with area-constant"
              << " pressure) and Zplane or six parameters: k_force, areaFlag (2/3 without/with area-pressure"
              << " increase with area decrease), Zplane, equilibrium volume,volume" 
              << " decrease factor and pressure increase factor." << std::endl;
    exit(EXIT_FAILURE);
  }
  if (paraValue[1]!=0.0 && paraValue[1]!=1.0 && paraValue[1]!=2.0 && paraValue[1]!=3.0) {
    std::cerr << "VertexFromCellPlaneTriangular::VertexFromCellPlaneTriangular() " 
	      << "areaFlag must be zero/two (no area included) or one/three (area included)." << std::endl;
    exit(EXIT_FAILURE);
  }
  if (indValue.size() != 0) {
    std::cerr << "VertexFromCellPlaneTriangular::VertexFromCellPlaneTriangular() " 
	      << std::endl
	      << "No variable index used." << std::endl;
    exit(EXIT_FAILURE);
  }
  
  setId("VertexFromCellPlaneTriangular");
  setParameter(paraValue);
  setVariableIndex(indValue);
  
  std::vector<std::string> tmp(numParameter());
  tmp[0] = "k_force";
  tmp[1] = "areaFlag";
  tmp[2] = "Zplane";

  if (paraValue.size() == 6){
  tmp[3] = "V0";
  tmp[4] = "Vfactor";
  tmp[5] = "Pfactor";
  }
  setParameterId(tmp);
}

void VertexFromCellPlaneTriangular::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs)
{
  size_t dimension = vertexData[0].size();
  if (dimension!=3) {
    std::cerr << "VertexFromCellPlaneTriangular::VertexFromCellPlaneTriangular() " 
	      << "Only implemented for three dimensions." << std::endl;
    exit(EXIT_FAILURE);
  }

  if (parameter(1) == 2 || parameter(1) == 3){
    
    
    double totalVolume=0;  
    size_t numCell = cellData.size();
    assert (numCell==T.numCell());  
    
    
    
    for (size_t cellIndex=0; cellIndex<numCell; ++cellIndex) { 
      Cell cell = T.cell(cellIndex);
      if (cell.numVertex()!=3) {
        std::cerr << "VertexFromCellPlaneTriangular::VertexFromCellPlaneTriangular() " 
                  << "Only implemented for triangular cells." << std::endl;
        exit(EXIT_FAILURE);
      }
      size_t v1 = T.cell(cellIndex).vertex(0)->index();
      size_t v2 = T.cell(cellIndex).vertex(1)->index();
      size_t v3 = T.cell(cellIndex).vertex(2)->index();
      DataMatrix position(3,vertexData[v1]);
      position[1] = vertexData[v2];
      position[2] = vertexData[v3];
      
      std::vector<double> normal = cell.getNormalTriangular(vertexData);
      double norm=0.0;
      for (size_t d=0; d<dimension; ++d)
        norm += normal[d]*normal[d];
      if (norm != 1.0) {
        norm = std::sqrt(norm);
        assert(norm>0.0);
        double normFac = 1.0/norm; 
        for (size_t d=0; d<dimension; ++d)
          normal[d] *= normFac;
      }
      double A= cell.calculateVolume(vertexData)*normal[2];
      totalVolume += A*(((position[0][2]+position[1][2]+position[2][2])/3)-parameter(2));
    }
    

    for (size_t cellIndex=0; cellIndex<numCell; ++cellIndex) { 
      Cell cell = T.cell(cellIndex);
      
      std::vector<double> normal = cell.getNormalTriangular(vertexData);
      double norm=0.0;
      for (size_t d=0; d<dimension; ++d)
        norm += normal[d]*normal[d];
      if (norm != 1.0) {
        norm = std::sqrt(norm);
        assert(norm>0.0);
        double normFac = 1.0/norm; 
        for (size_t d=0; d<dimension; ++d)
          normal[d] *= normFac;
      }
      
      // Get the cell size
      double A=1.0;
      if (parameter(1)==3.0)
        A = cell.calculateVolume(vertexData)/cell.numVertex();
      
      double coeff =((parameter(5)-1)*parameter(0)*totalVolume/((parameter(4)-1)*parameter(3))
                     +(parameter(4)-parameter(5))*parameter(0)/(parameter(4)-1))* A;
      //update the vertex derivatives
      for (size_t k=0; k<cell.numVertex(); ++k) {
        double vCoeff=coeff;
        //if (cell.vertex(k)->isBoundary(T.background()))
        //vCoeff *= 1.5;
        for (size_t d=0; d<dimension; ++d) {
          vertexDerivs[cell.vertex(k)->index()][d] += vCoeff * normal[d];
        }
      }	
    }
  }
  
  if (parameter(1) == 0 || parameter(1) ==1 ){
    double totalVolume=0;  
    size_t numCell = cellData.size();
    assert (numCell==T.numCell());  
    for (size_t cellIndex=0; cellIndex<numCell; ++cellIndex) { 
      Cell cell = T.cell(cellIndex);
      if (cell.numVertex()!=3) {
        std::cerr << "VertexFromCellPlaneTriangular::VertexFromCellPlaneTriangular() " 
                  << "Only implemented for triangular cells." << std::endl;
        exit(EXIT_FAILURE);
      }
      size_t v1 = T.cell(cellIndex).vertex(0)->index();
      size_t v2 = T.cell(cellIndex).vertex(1)->index();
      size_t v3 = T.cell(cellIndex).vertex(2)->index();
      DataMatrix position(3,vertexData[v1]);
      position[1] = vertexData[v2];
      position[2] = vertexData[v3];
      
      std::vector<double> normal = cell.getNormalTriangular(vertexData);
      double norm=0.0;
      for (size_t d=0; d<dimension; ++d)
        norm += normal[d]*normal[d];
      if (norm != 1.0) {
        norm = std::sqrt(norm);
        assert(norm>0.0);
        double normFac = 1.0/norm; 
        for (size_t d=0; d<dimension; ++d)
          normal[d] *= normFac;
      }
      double A= cell.calculateVolume(vertexData)*normal[2];
      totalVolume += A*(((position[0][2]+position[1][2]+position[2][2])/3)-parameter(2));
    }
    
    //  std:: cerr<<"total Volume:  "<<totalVolume<<std::endl;

    for (size_t n = 0; n < T.numCell(); ++n) { 
      Cell cell = T.cell(n);
      if (cell.numVertex()!=3) {
        std::cerr << "VertexFromCellPlaneTriangular::VertexFromCellPlaneTriangular() " 
                  << "Only implemented for triangular cells." << std::endl;
        exit(EXIT_FAILURE);
      }
      std::vector<double> normal = cell.getNormalTriangular(vertexData);
      double norm=0.0;
      for (size_t d=0; d<dimension; ++d)
        norm += normal[d]*normal[d];
      if (norm != 1.0) {
        norm = std::sqrt(norm);
        assert(norm>0.0);
        double normFac = 1.0/norm; 
        for (size_t d=0; d<dimension; ++d)
          normal[d] *= normFac;
      }
      
      // Get the cell size
      double A=1.0;
      if (parameter(1)==1.0)
        A = cell.calculateVolume(vertexData)/cell.numVertex();
      
      double coeff = parameter(0) * A;
      //update the vertex derivatives
      for (size_t k=0; k<cell.numVertex(); ++k) {
          double vCoeff=coeff;
          //if (cell.vertex(k)->isBoundary(T.background()))
          //vCoeff *= 1.5;
          for (size_t d=0; d<dimension; ++d) {
            vertexDerivs[cell.vertex(k)->index()][d] += vCoeff * normal[d];
          }
      }	
      // For saving normals in direction used for test plotting
      //for (size_t d=0; d<dimension; ++d)
      //cellData[cell.index()][d] = normal[d];
    }    
  }
}



VertexFromForce::
VertexFromForce(std::vector<double> &paraValue, 
		       std::vector< std::vector<size_t> > 
		       &indValue ) 
{  
  //Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()>3 ) {
    std::cerr << "VertexFromForce::"
	      << "VertexFromForce() "
	      << "Uses a force vector that should be in one (x), two (x,y) or three (x,y,z) "
	      << "dimensions." << std::endl;
    exit(0);
  }
  if( indValue.size() != 1 || indValue[0].size() < 1 ) {
    std::cerr << "VertexFromForce::"
	      << "VertexFromForce() "
	      << "List of vertex indices given in first level." << std::endl;
    exit(0);
  }
  //Set the variable values
  //
  setId("VertexFromForce");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "F_x";
  if (numParameter()>1)
    tmp[1] = "F_y";
  if (numParameter()==3)
    tmp[2] = "F_z";
  setParameterId( tmp );
}

void VertexFromForce::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) {
  
  //Do the update for each vertex in the list given.
  for (size_t k=0 ; k<numVariableIndex(0); ++k) {
    size_t i = variableIndex(0,k);
    assert(i<vertexData.size());
    for (size_t d=0; d<vertexData[i].size(); ++d)
      if( numParameter()>d )
	vertexDerivs[i][d] += parameter(d);
  }
}

VertexFromForceLinear::
VertexFromForceLinear(std::vector<double> &paraValue, 
		       std::vector< std::vector<size_t> > 
		       &indValue ) 
{  
  //Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()<2 || paraValue.size()>4 ) {
    std::cerr << "VertexFromForceLinear::"
	      << "VertexFromForceLinear() "
	      << "Uses a force vector that should be in one (x), two (x,y) or three (x,y,z) "
	      << "dimensions plus a deltaT that sets the time the linear increase is applied." 
	      << std::endl;
    exit(0);
  }
  if( indValue.size() != 1 || indValue[0].size() < 1 ) {
    std::cerr << "VertexFromForceLinear::"
	      << "VertexFromForceLinear() "
	      << "List of vertex indices given in first level." << std::endl;
    exit(0);
  }
  //Set the variable values
  //
  setId("VertexFromForceLinear");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "F_x";
  if (numParameter()==2)
    tmp[1] = "deltaT";
  else if (numParameter()==3) {
    tmp[1] = "F_y";
    tmp[2] = "deltaT";
  }
  else {
    tmp[1] = "F_y";
    tmp[2] = "F_z";
    tmp[3] = "deltaT";
  }
  timeFactor_=0.0;
  setParameterId( tmp );
}

void VertexFromForceLinear::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) {
  
  //Do the update for each vertex in the list given.
  for (size_t k=0 ; k<numVariableIndex(0); ++k) {
    size_t i = variableIndex(0,k);
    assert(i<vertexData.size());
 
    for (size_t d=0; d<vertexData[i].size(); ++d)
    if( numParameter()>d )
    	vertexDerivs[i][d] += timeFactor_*parameter(d);
   
   
  }
 
}

void VertexFromForceLinear::update(Tissue &T,
				   DataMatrix &cellData,
				   DataMatrix &wallData,
				   DataMatrix &vertexData,
				   double h)
{
  if (timeFactor_ < 1.0 ) {
    timeFactor_ += h/parameter(numParameter()-1);
  }
  if (timeFactor_ >1.0)
    timeFactor_=1.0;
 
}




VertexFromBall::
VertexFromBall(std::vector<double> &paraValue, 
		       std::vector< std::vector<size_t> > 
		       &indValue ) 
{  
  //Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=5 && paraValue.size()!=8 ) {
    std::cerr << "VertexFromBall::"
	      << "VertexFromBall() "
	      << "Puts a ball(sphere) of a given radius (radius) in a given "
	      << "position (x,y,z) around meriestem or moves it toward meristem by a given velocity vector "
	      << "5 parameters for static : radius, x, y, z. Kforce" << std::endl
	      << "8 parameters for dynamic : radius, x, y, z, Kforce, dx, dy, dz."
	      << std::endl;
    exit(0);
  }
  if( indValue.size() != 0 ) {
    std::cerr << "VertexFromBall::"
	      << "VertexFromBall() "
	      << "No variable indices used." << std::endl;
    exit(0);
  }
  //Set the variable values
  //
  setId("VertexFromBall");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "Radius";
  tmp[1] = "xc";
  tmp[2] = "yc";
  tmp[3] = "zc";
  tmp[4] = "Kforce";

  setParameterId( tmp );
}

void VertexFromBall::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) {
  
  //Do the update for each vertex .
  size_t numVertex = T.numVertex();
  for (size_t vertexIndex=0 ; vertexIndex<numVertex; ++vertexIndex) {
    double Radius=parameter(0);
    double Xc=parameter(1);
    double Yc=parameter(2);
    double Zc=parameter(3);
    double Kforce=parameter(4);
    DataMatrix position(1,vertexData[vertexIndex]);
    double d2=(position[0][0]-Xc)*(position[0][0]-Xc)+
      (position[0][1]-Yc)*(position[0][1]-Yc)+
      (position[0][2]-Zc)*(position[0][2]-Zc);
    if( d2 < Radius*Radius ){
      vertexDerivs[vertexIndex][0]+=Kforce*(Radius-std::sqrt(d2))*std::sqrt(Radius-std::sqrt(d2))*(position[0][0]-Xc)/std::sqrt(d2);
      vertexDerivs[vertexIndex][1]+=Kforce*(Radius-std::sqrt(d2))*std::sqrt(Radius-std::sqrt(d2))*(position[0][1]-Yc)/std::sqrt(d2);
      vertexDerivs[vertexIndex][2]+=Kforce*(Radius-std::sqrt(d2))*std::sqrt(Radius-std::sqrt(d2))*(position[0][2]-Zc)/std::sqrt(d2);
    }
  }
}


void VertexFromBall::update(Tissue &T,
			    DataMatrix &cellData,
			    DataMatrix &wallData,
			    DataMatrix &vertexData,
			    double h)
{
  if( numParameter()>5 ) {
    setParameter(1,parameter(1)+h*parameter(5));
    setParameter(2,parameter(2)+h*parameter(6));
    setParameter(3,parameter(3)+h*parameter(7));
  }
}










VertexFromParabolid::
VertexFromParabolid(std::vector<double> &paraValue, 
		       std::vector< std::vector<size_t> > 
		       &indValue ) 
{  
  //Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=5 && paraValue.size()!=6 ) {
    std::cerr << "VertexFromParabolid::"
	      << "VertexFromParabolid() "
	      << "Puts a vertical(z) parabolid z=a((x-xc)2+(y-yc)2)+b with given properties in a given "
	      << "position (x,y) and b below or above template or  moves it upward or downwars by a given velocity"
	      << "5 parameters for static : a, xc, yc, b. Kforce" << std::endl
	      << "6 parameters for dynamic : a, xc, yc, b, Kforce, v."
	      << std::endl;
    exit(0);
  }
  if( indValue.size() != 0 ) {
    std::cerr << "VertexFromParabolid::"
	      << "VertexFromParabolid() "
	      << "No variable indices used." << std::endl;
    exit(0);
  }
  //Set the variable values
  //
  setId("VertexFromParabolid");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "a";
  tmp[1] = "xc";
  tmp[2] = "yc";
  tmp[3] = "b";
  tmp[4] = "Kforce";

  setParameterId( tmp );
}

void VertexFromParabolid::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) {
  
  //Do the update for each vertex .
  size_t numVertex = T.numVertex();
  for (size_t vertexIndex=0 ; vertexIndex<numVertex; ++vertexIndex) {
    double a=parameter(0);
    double Xc=parameter(1);
    double Yc=parameter(2);
    double b=parameter(3);
    double Kforce=parameter(4);
    DataMatrix position(1,vertexData[vertexIndex]);
    double d=position[0][2]-a*(position[0][0]-Xc)*(position[0][0]-Xc)-a*(position[0][1]-Yc)*(position[0][1]-Yc)-b;
    if( d < 0 ){

      double m=2*a*(position[0][0]-Xc);
      double n=2*a*(position[0][1]-Yc);

      // vertexDerivs[vertexIndex][0]+=-Kforce*(m/std::sqrt(1+m*m+n*n));
      // vertexDerivs[vertexIndex][1]+=-Kforce*(n/std::sqrt(1+m*m+n*n));
      // vertexDerivs[vertexIndex][2]+=Kforce*(1/std::sqrt(1+m*m+n*n));
      vertexDerivs[vertexIndex][2]+=-Kforce*(d);
    }
  }
}


void VertexFromParabolid::update(Tissue &T,
			    DataMatrix &cellData,
			    DataMatrix &wallData,
			    DataMatrix &vertexData,
			    double h)
{
  if( numParameter()>5 ) 
    setParameter(3,parameter(3)+h*parameter(5));
   
}
























VertexFromExternalWall::
VertexFromExternalWall(std::vector<double> &paraValue, 
		       std::vector< std::vector<size_t> > 
		       &indValue ) 
{  
  //Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=12 ) {
    std::cerr << "VertexFromExternalWall::"
	      << "VertexFromExternalWall()"
	      << "Puts a wall in a given point(x,y,z) with a given normal vector(nx,ny,nz)"
              << "with a given height beween zmin and zmax around meriestem and  moves it"
              << "toward meristem with a given velocity vector(dx,dy,dz)"
	      << "12 parameters used:x0, y0, z0, nx, ny, nz, zmin, zmax, dx, dy, dz, Kforce" << std::endl
              << std::endl;
    exit(0);
  }
  if( indValue.size() != 0 ) {
    std::cerr << "VertexFromExternalWall::"
	      << "VertexFromExternalWall() "
	      << "No variable indices used." << std::endl;
    exit(0);
  }
  //Set the variable values
  //
  setId("VertexFromExternalWall");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "x0";
  tmp[1] = "y0";
  tmp[2] = "z0";
  tmp[3] = "nx";
  tmp[4] = "ny";
  tmp[5] = "nz";
  tmp[6] = "zmin";
  tmp[7] = "zmax";
  tmp[8] = "dx";
  tmp[9] = "dy";
  tmp[10] = "dz";
  tmp[11] = "Kforce";

  setParameterId( tmp );
}

void VertexFromExternalWall::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) {
  
  //Do the update for each vertex .
  size_t numVertex = T.numVertex();
  for (size_t vertexIndex=0 ; vertexIndex<numVertex; ++vertexIndex) {
    double X0=parameter(0);
    double Y0=parameter(1);
    double Z0=parameter(2);
    double nx=parameter(3);
    double ny=parameter(4);
    double nz=parameter(5);
    double Zmin=parameter(6);
    double Zmax=parameter(7);
    double Kforce=parameter(11);
    DataMatrix position(1,vertexData[vertexIndex]);
    double d=(nx*position[0][0]+ny*position[0][1]+nz*position[0][2]-nx*X0-ny*Y0-nz*Z0)/std::sqrt(nx*nx+ny*ny+nz*nz);
    if( d < 0 && position[0][2] > Zmin && position[0][2] < Zmax ){
      vertexDerivs[vertexIndex][0]+=Kforce*(-d)*std::sqrt(-d)*nx/std::sqrt(nx*nx+ny*ny+nz*nz);
      vertexDerivs[vertexIndex][1]+=Kforce*(-d)*std::sqrt(-d)*ny/std::sqrt(nx*nx+ny*ny+nz*nz);
      vertexDerivs[vertexIndex][2]+=Kforce*(-d)*std::sqrt(-d)*nz/std::sqrt(nx*nx+ny*ny+nz*nz);
    }
  }
}


void VertexFromExternalWall::update(Tissue &T,
			    DataMatrix &cellData,
			    DataMatrix &wallData,
			    DataMatrix &vertexData,
			    double h)
{
  if( numParameter()>5 ) {
    setParameter(0,parameter(0)+h*parameter(8));
    setParameter(1,parameter(1)+h*parameter(9));
    setParameter(2,parameter(2)+h*parameter(10));
  }
}


TemplateVolumeChange::TemplateVolumeChange(std::vector<double> &paraValue,
                                           std::vector< std::vector<size_t> > &indValue )
{ // Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=0 ) {
    std::cerr << "TemplateVolumeChange:: "
	      << "TemplateVolumeChange() "
	      << "Calculates change in template volume and its time derivative " 
              << "and total Derivative and stores them in the given indices in cellData vector "
              << "uses no parameter "
              << std::endl;
    exit(0);
  }
 if( indValue.size() != 1 || indValue[0].size() != 6 ) {
    std::cerr << "TemplateVolumeChange:: "
	      << "TemplateVolumeChange() "
	      << "1 level with 6 variable indices used: "
              << "cell-index-VolumeChange       component-index-VolumeChange "
              << "cell-index-deltaVolumeChange  component-index-deltaVolumeChange "
              << "cell-index-totalDerivative    component-index-totalDerivative "
              << std::endl;
    exit(0);
  }
  //Set the variable values
  //
  setId("TemplateVolumeChange");
  setParameter(paraValue);  
  setVariableIndex(indValue);

}

void TemplateVolumeChange:: initiate(Tissue &T,
                                     DataMatrix &cellData,
                                     DataMatrix &wallData,
                                     DataMatrix &vertexData,
                                     DataMatrix &cellDerivs,
                                     DataMatrix &wallDerivs,
                                     DataMatrix &vertexDerivs)
{
  size_t dimension=3; //Only implemented for 3D models
  assert (dimension==vertexData[0].size());
  vertexDataRest.resize(vertexData.size());
  for(size_t i=0 ; i<vertexDataRest.size() ;++i)
    vertexDataRest[i].resize(vertexData[i].size());

  size_t numVertices = T.numVertex();

  for (size_t vertexIndex=0; vertexIndex<numVertices; ++vertexIndex) 
    for (size_t d=0; d<dimension; ++d) {
      vertexDataRest[vertexIndex][d]=vertexData[vertexIndex][d];
    }
}

void TemplateVolumeChange::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs )
{
  size_t numVertices =T.numVertex();
  VolumeChange=0;
  for(size_t vertexIndex=0; vertexIndex<numVertices; ++vertexIndex) 
    { VolumeChange += std::sqrt( (vertexData[vertexIndex][0]-vertexDataRest[vertexIndex][0])*
                                 (vertexData[vertexIndex][0]-vertexDataRest[vertexIndex][0])+
                                 (vertexData[vertexIndex][1]-vertexDataRest[vertexIndex][1])*
                                 (vertexData[vertexIndex][1]-vertexDataRest[vertexIndex][1])+
                                 (vertexData[vertexIndex][2]-vertexDataRest[vertexIndex][2])*
                                 (vertexData[vertexIndex][2]-vertexDataRest[vertexIndex][2])  );
    }
  deltaVolumeChange = VolumeChange - cellData[variableIndex(0,0)][variableIndex(0,1)];
  totalDerivative=0;
  for(size_t vertexIndex=0; vertexIndex<numVertices; ++vertexIndex) 
    { 
      if (vertexDerivs[0].size()>2)
        totalDerivative += std::sqrt( vertexDerivs[vertexIndex][0]*vertexDerivs[vertexIndex][0]+
                                      vertexDerivs[vertexIndex][1]*vertexDerivs[vertexIndex][1]+
                                      vertexDerivs[vertexIndex][2]*vertexDerivs[vertexIndex][2]  );
    }
  cellData[variableIndex(0,0)][variableIndex(0,1)]=VolumeChange;
  cellData[variableIndex(0,2)][variableIndex(0,3)]=deltaVolumeChange;
  cellData[variableIndex(0,4)][variableIndex(0,5)]=totalDerivative;
  //std:: cerr << VolumeChange<<" "<<deltaVolumeChange<<" "<<VolumeChangeDerT<<std::endl;
}

void TemplateVolumeChange:: update(Tissue &T,
                                   DataMatrix &cellData,
                                   DataMatrix &wallData,
                                   DataMatrix &vertexData,
                                   DataMatrix &vertexDerivs,
                                   double h)
{ 
}




CalculateAngleVectors::CalculateAngleVectors(std::vector<double> &paraValue,
                                           std::vector< std::vector<size_t> > &indValue )
{ // Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=0 ) {
    std::cerr << "CalculateAngleVectors:: "
	      << "CalculateAngleVectors() "
	      << "Calculates abs(cos(...)) of angle between two 3d vectors(starting from given indices) in cellData vector "
              << "and stores it in the given index in cellData vector, uses no parameter "
              << std::endl;
    exit(0);
  }
 if( indValue.size() != 2 || indValue[0].size() != 2 || indValue[1].size() != 1 ) {
    std::cerr << "CalculateAngleVectors:: "
	      << "CalculateAngleVectors() "
	      << "1st level with 2 variable indices used: "
              << "start index for 1st vector   start index for 2nd vector "
              << "2nd level with 1 variable index used: "
              << "store index for angle(deg) "
              << std::endl;
    exit(0);
  }
  //Set the variable values
  //
  setId("CalculateAngleVectors");
  setParameter(paraValue);  
  setVariableIndex(indValue);

}

// void CalculateAngleVectors:: initiate(Tissue &T,
//                                      DataMatrix &cellData,
//                                      DataMatrix &wallData,
//                                      DataMatrix &vertexData,
//                                      DataMatrix &cellDerivs,
//                                      DataMatrix &wallDerivs,
//                                      DataMatrix &vertexDerivs)
// {
//   size_t dimension=3; //Only implemented for 3D models
//   assert (dimension==vertexData[0].size());
//   vertexDataRest.resize(vertexData.size());
//   for(size_t i=0 ; i<vertexDataRest.size() ;++i)
//     vertexDataRest[i].resize(vertexData[i].size());

//   size_t numVertices = vertexData.size();

//   for (size_t vertexIndex=0; vertexIndex<numVertices; ++vertexIndex) 
//     for (size_t d=0; d<dimension; ++d) {
//       vertexDataRest[vertexIndex][d]=vertexData[vertexIndex][d];
//   }
// }

void CalculateAngleVectors::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) {
  size_t numCells =T.numCell();
  for(size_t cellIndex=0; cellIndex<numCells; ++cellIndex) 
    { double teta=0;
      double tmp=std::sqrt( cellData[cellIndex][variableIndex(0,0)  ]*cellData[cellIndex][variableIndex(0,0)  ] +
                            cellData[cellIndex][variableIndex(0,0)+1]*cellData[cellIndex][variableIndex(0,0)+1] +
                            cellData[cellIndex][variableIndex(0,0)+2]*cellData[cellIndex][variableIndex(0,0)+2]   );
      cellData[cellIndex][variableIndex(0,0)  ]=cellData[cellIndex][variableIndex(0,0)  ]/tmp;
      cellData[cellIndex][variableIndex(0,0)+1]=cellData[cellIndex][variableIndex(0,0)+1]/tmp;
      cellData[cellIndex][variableIndex(0,0)+2]=cellData[cellIndex][variableIndex(0,0)+2]/tmp;

      tmp=std::sqrt( cellData[cellIndex][variableIndex(0,1)  ]*cellData[cellIndex][variableIndex(0,1)  ] +
                     cellData[cellIndex][variableIndex(0,1)+1]*cellData[cellIndex][variableIndex(0,1)+1] +
                     cellData[cellIndex][variableIndex(0,1)+2]*cellData[cellIndex][variableIndex(0,1)+2]   );
      cellData[cellIndex][variableIndex(0,1)  ]=cellData[cellIndex][variableIndex(0,1)  ]/tmp;
      cellData[cellIndex][variableIndex(0,1)+1]=cellData[cellIndex][variableIndex(0,1)+1]/tmp;
      cellData[cellIndex][variableIndex(0,1)+2]=cellData[cellIndex][variableIndex(0,1)+2]/tmp;

      teta=std::abs( cellData[cellIndex][variableIndex(0,0)  ]*cellData[cellIndex][variableIndex(0,1)  ] +
                     cellData[cellIndex][variableIndex(0,0)+1]*cellData[cellIndex][variableIndex(0,1)+1] +
                     cellData[cellIndex][variableIndex(0,0)+2]*cellData[cellIndex][variableIndex(0,1)+2]   );
  
      cellData[cellIndex][variableIndex(1,0)]=teta;
    }
}



CalculateAngleVectorXYplane::CalculateAngleVectorXYplane(std::vector<double> &paraValue,
                                           std::vector< std::vector<size_t> > &indValue )
{ // Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=0 ) {
    std::cerr << "CalculateAngleVectorXYplane:: "
	      << "CalculateAngleVectorXYplane() "
	      << "Calculates abs(cos(...)) of angle between a vector and XY plane "
              << "and stores it in the given index in cellData vector, uses no parameter "
              << std::endl;
    exit(0);
  }
 if( indValue.size() != 2 || indValue[0].size() != 1 || indValue[1].size() != 1 ) {
    std::cerr << "CalculateAngleVectorXYplane:: "
	      << "CalculateAngleVectorXYplane() "
	      << "1st level with 1 variable index used: "
              << "start index for the vector  "
              << "2nd level with 1 variable index used: "
              << "store index for cos(angle) "
              << std::endl;
    exit(0);
  }
  //Set the variable values
  //
  setId("CalculateAngleVectorXYplane");
  setParameter(paraValue);  
  setVariableIndex(indValue);

}


void CalculateAngleVectorXYplane::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) {
  size_t numCells =T.numCell();
  for(size_t cellIndex=0; cellIndex<numCells; ++cellIndex) 
    { double teta=0;
      double temp=std::sqrt( cellData[cellIndex][variableIndex(0,0)  ]*cellData[cellIndex][variableIndex(0,0)  ] +
                            cellData[cellIndex][variableIndex(0,0)+1]*cellData[cellIndex][variableIndex(0,0)+1]  );
      double pi=3.14159265;
      if(temp<0.00000001)
	teta=pi/2;
      else
	teta=std::atan(cellData[cellIndex][variableIndex(0,0)+2]/temp);
       
      cellData[cellIndex][variableIndex(1,0)]=teta;
    }
}



AngleVector::AngleVector(std::vector<double> &paraValue,
                                           std::vector< std::vector<size_t> > &indValue )
{ // Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=1 ) {
    std::cerr << "AngleVector:: "
	      << "AngleVector() "
	      << "Calculates  angle between a 3d vector(starting from given index) in cellData vector and "
              << "a given axes(x,y,z) and stores it in the given index in cellData vector, uses one parameter "
	      << "for specifying the axes "
              << std::endl;
    exit(0);
  }
 if( indValue.size() != 2 || indValue[0].size() != 1 || indValue[1].size() != 1 ) {
    std::cerr << "AngleVector:: "
	      << "AngleVector() "
	      << "1st level with 1 variable index used: "
              << "start index for the vector  "
              << "2nd level with 1 variable index used: "
              << "store index for angle(deg) "
              << std::endl;
    exit(0);
  }
  //Set the variable values
  //
  setId("AngleVector");
  setParameter(paraValue);  
  setVariableIndex(indValue);

 // Set the parameter identities
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "axes";   // axes specifier

  setParameterId( tmp );

  if( parameter(0)!=0 &&  parameter(0)!=1 &&  parameter(0)!=2 ) {
    std::cerr << " AngleVector:: AngleVector() "
              << " 0nly 0(X), 1(Y) or 2(Z) " << std::endl;
    exit(0);
  }
  if( parameter(0)==1 ||  parameter(0)==2 ) {
    std::cerr << " AngleVector:: AngleVector() "
              << " The code should be modified for 3d " << std::endl;
    exit(0);
  }

}



void AngleVector::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) {
  size_t numCells =T.numCell();
  for(size_t cellIndex=0; cellIndex<numCells; ++cellIndex) 
    { double teta=0;
      double tmp=std::sqrt( cellData[cellIndex][variableIndex(0,0)  ]*cellData[cellIndex][variableIndex(0,0)  ] +
                            cellData[cellIndex][variableIndex(0,0)+1]*cellData[cellIndex][variableIndex(0,0)+1] +
                            cellData[cellIndex][variableIndex(0,0)+2]*cellData[cellIndex][variableIndex(0,0)+2]   );
      cellData[cellIndex][variableIndex(0,0)  ]=cellData[cellIndex][variableIndex(0,0)  ]/tmp;
      cellData[cellIndex][variableIndex(0,0)+1]=cellData[cellIndex][variableIndex(0,0)+1]/tmp;
      cellData[cellIndex][variableIndex(0,0)+2]=cellData[cellIndex][variableIndex(0,0)+2]/tmp;

      double pi= 3.141592;

      if(parameter(0)==0){    
	teta=180*std::acos(cellData[cellIndex][variableIndex(0,0)  ])/pi;
	if(cellData[cellIndex][variableIndex(0,0)+1]<0)	{
	  teta=180-teta;
	}
	cellData[cellIndex][variableIndex(1,0)]=teta;

      }

      // should be modified for 3d
      if(parameter(0)==1){    
	teta=180*std::acos(cellData[cellIndex][variableIndex(0,0)+1])/pi;
	
	cellData[cellIndex][variableIndex(1,0)]=teta;
      }

      if(parameter(0)==2){    
	teta=180*std::acos(cellData[cellIndex][variableIndex(0,0)+2])/pi;
	
	cellData[cellIndex][variableIndex(1,0)]=teta;
      }

    }
}

VertexFromHypocotylGrowth::
VertexFromHypocotylGrowth(std::vector<double> &paraValue, 
			  std::vector< std::vector<size_t> > 
			  &indValue ) 
{  
  // Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=4 ) {
    std::cerr << "VertexFromHypocotylGrowth::"
	      << "VertexFromHypocotylGrowth()"
	      << "Applies forces F axially to the vertices in the region between Y0-a and Y0+a"
              << "uses 4 parameters: Y0, a, d and F" << std::endl
              << std::endl;
    exit(EXIT_FAILURE);
  }
  if( indValue.size() != 0 ) {
    std::cerr << "VertexFromHypocotylGrowth::"
	      << "VertexFromHypocotylGrowth() "
	      << "No variable indices used." << std::endl;
    exit(EXIT_FAILURE);
  }
  //Set the variable values
  //
  setId("VertexFromHypocotylGrowth");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "y0";      // the center of the growing region
  tmp[1] = "a";       // growing region interwall 
  tmp[2] = "d";       // force application interwall
  tmp[3] = "F";       // Force for each node
 
  setParameterId( tmp );
}

void VertexFromHypocotylGrowth::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) {
  
  //Do the update for each vertex .
  size_t numVertex = T.numVertex();
  for (size_t vertexIndex=0 ; vertexIndex<numVertex; ++vertexIndex) {
    double Y0=parameter(0);
    double a=parameter(1);
    double d=parameter(2);
    double F=parameter(3);
    
    DataMatrix position(1,vertexData[vertexIndex]);
    

    if( position[0][2] > Y0+a && position[0][2] < Y0+a+d ){

      //vertexDerivs[vertexIndex][0]+=Kforce*(-d)*std::sqrt(-d)*nx/std::sqrt(nx*nx+ny*ny+nz*nz);
      //vertexDerivs[vertexIndex][1]+=Kforce*(-d)*std::sqrt(-d)*ny/std::sqrt(nx*nx+ny*ny+nz*nz);
      vertexDerivs[vertexIndex][2]+=F;
    }
    if( position[0][2] < Y0-a && position[0][2] > Y0-a-d ){
      
      //vertexDerivs[vertexIndex][0]+=Kforce*(-d)*std::sqrt(-d)*nx/std::sqrt(nx*nx+ny*ny+nz*nz);
      //vertexDerivs[vertexIndex][1]+=Kforce*(-d)*std::sqrt(-d)*ny/std::sqrt(nx*nx+ny*ny+nz*nz);
      vertexDerivs[vertexIndex][2]-=F;
    }
  }
}




void VertexFromHypocotylGrowth::update(Tissue &T,
			    DataMatrix &cellData,
			    DataMatrix &wallData,
			    DataMatrix &vertexData,
			    double h)
{
 
}


// VertexFromHypocotylGrowth::
// VertexFromHypocotylGrowth(std::vector<double> &paraValue, 
// 		       std::vector< std::vector<size_t> > 
// 		       &indValue ) 
// {  
//   //Do some checks on the parameters and variable indeces
//   //
//   if( paraValue.size()!=8 ) {
//     std::cerr << "VertexFromHypocotylGrowth::"
// 	      << "VertexFromHypocotylGrowth()"
// 	      << "Applies forces by simple springs from a growing center in the Hypocotyl"
//               << "uses 7 parameters: Y0, a, d, delta, lambda, b, K, epsilon" << std::endl
//               << std::endl;
//     exit(0);
//   }
//   if( indValue.size() != 0 ) {
//     std::cerr << "VertexFromHypocotylGrowth::"
// 	      << "VertexFromHypocotylGrowth() "
// 	      << "No variable indices used." << std::endl;
//     exit(0);
//   }
//   //Set the variable values
//   //
//   setId("VertexFromHypocotylGrowth");
//   setParameter(paraValue);  
//   setVariableIndex(indValue);
  
//   //Set the parameter identities
//   //
//   std::vector<std::string> tmp( numParameter() );
//   tmp[0] = "y0";      // the center of the growing region
//   tmp[1] = "a";       // growing region interwall 
//   tmp[2] = "d";       // distance between the boundary of the growing region and epidermis
//   tmp[3] = "delta";   // maximum displacement ( at the center of the growing region)
//   tmp[4] = "lambda";  // Gaussian distribution of displacement of the growing center
//   tmp[5] = "b";       // interwall on epidermis recieving the force
//   tmp[6] = "Ks";       // spring constant
//   tmp[7] = "epsilon"; // deltaX for evaluating the integral force

//   setParameterId( tmp );
// }

// void VertexFromHypocotylGrowth::
// derivs(Tissue &T,
//        DataMatrix &cellData,
//        DataMatrix &wallData,
//        DataMatrix &vertexData,
//        DataMatrix &cellDerivs,
//        DataMatrix &wallDerivs,
//        DataMatrix &vertexDerivs ) {
  
//   //Do the update for each vertex .
//   size_t numVertex = T.numVertex();
//   for (size_t vertexIndex=0 ; vertexIndex<numVertex; ++vertexIndex) {
//     double Y0=parameter(0);
//     double a=parameter(1);
//     double d=parameter(2);
//     double delta=parameter(3);
//     double lambda=parameter(4);
//     double b=parameter(5);
//     double Ks=parameter(6);
//     double epsilon=parameter(7);
   
//     DataMatrix position(1,vertexData[vertexIndex]);
    

//     if( position[0][2] > Y0-b && position[0][2] < Y0+b ){

//       double Faxial=0;
//       double Fradial=0;
//       double X=Y0-a;
//       double Y=position[0][2];
//       while (X < Y0+a){
// 	int sgnx=0;
//         if ( X-Y0 > 0 ) sgnx=1;
// 	else sgnx=-1;

// 	Faxial -=(epsilon/(2*a))*Ks 
// 	  *(    (  std::sqrt( (Y-X)*(Y-X)+d*d )
// 		   /std::sqrt((Y-X-sgnx*delta*exp(-lambda*X*X))*(Y-X-sgnx*delta*exp(-lambda*X*X))+d*d )   
// 		   )       
// 		-1 )
//           *(Y-X-sgnx*delta*exp(-lambda*X*X));    
	 	  
// 	Fradial =(epsilon/(2*a))*Ks 
// 	  *(    (  std::sqrt( (Y-X)*(Y-X)+d*d )
// 		   /std::sqrt((Y-X-sgnx*delta*exp(-lambda*X*X))*(Y-X-sgnx*delta*exp(-lambda*X*X))+d*d )   
// 		   )       
// 		-1 )
//           *(d);    
// 	  X+=epsilon;	 
//       }
//       //vertexDerivs[vertexIndex][0]+=Kforce*(-d)*std::sqrt(-d)*nx/std::sqrt(nx*nx+ny*ny+nz*nz);
//       //vertexDerivs[vertexIndex][1]+=Kforce*(-d)*std::sqrt(-d)*ny/std::sqrt(nx*nx+ny*ny+nz*nz);
//       vertexDerivs[vertexIndex][2]+=Faxial;
//     }
//   }
// }




// void VertexFromHypocotylGrowth::update(Tissue &T,
// 			    DataMatrix &cellData,
// 			    DataMatrix &wallData,
// 			    DataMatrix &vertexData,
// 			    double h)
// {
 
// }






DebugReaction::DebugReaction(std::vector<double> &paraValue,
			     std::vector< std::vector<size_t> > &indValue)
{
}

void DebugReaction::derivs(Tissue &T,
			   DataMatrix &cellData,
			   DataMatrix &wallData,
			   DataMatrix &vertexData,
			   DataMatrix &cellDerivs,
			   DataMatrix &wallDerivs,
			   DataMatrix &vertexDerivs)
{
	for (size_t i = 0; i < T.numCell(); ++i) {
		Cell cell = T.cell(i);
		
// 		std::cerr << "Cell " << cell.index() << std::endl;

// 		for (size_t i = 0; i < cell.numVertex(); ++i) {
// 			Vertex *vertex = cell.vertex(i);
// 			std::cerr << "   Vertex: " << i << std::endl;
// 			std::cerr << "      x = " << vertexData[vertex->index()][0] << std::endl;
// 			std::cerr << "      y = " << vertexData[vertex->index()][1] << std::endl;
// 			std::cerr << "      z = " << vertexData[vertex->index()][2] << std::endl;
// 		}

		// This calculation should now be done in reaction CalculatePCAPlane
		//cell.calculatePCAPlane(vertexData);
		
		std::vector<double> N = cell.getNormalToPCAPlane();

		std::cerr << "x = " << N[0] << std::endl;
		std::cerr << "y = " << N[1] << std::endl;
		std::cerr << "z = " << N[2] << std::endl;
	}
}
