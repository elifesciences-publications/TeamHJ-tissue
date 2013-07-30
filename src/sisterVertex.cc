//
// Filename     : sisterVertex.cc
// Description  : Classes describing reactions related to sisterVertices
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : July 2013
// Revision     : $Id:$
//
#include<cmath>
#include"baseReaction.h"
#include"sisterVertex.h"
#include"tissue.h"

namespace SisterVertex {
  
  InitiateFromFile::
  InitiateFromFile(std::vector<double> &paraValue, 
		   std::vector< std::vector<size_t> > 
		   &indValue )
  {
    //Do some checks on the parameters and variable indices
    //
    if( paraValue.size()!=0 ) {
      std::cerr << "SisterVertex::InitiateFromFile::"
		<< "InitiateFromFile() "
		<< "No parameters should be provided." << std::endl;
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
    setId("SisterVertex::InitiateFromFile");
    setParameter(paraValue);  
    setVariableIndex(indValue);
  }
    
  void InitiateFromFile::
  initiate(Tissue &T,
	   DataMatrix &cellData,
	   DataMatrix &wallData,
	   DataMatrix &vertexData,
	   DataMatrix &cellDerivs,
	   DataMatrix &wallDerivs,
	   DataMatrix &vertexDerivs)
  {
    //open the file for reading
    std::ifstream IN("sister");
    if (!IN) {
      std::cerr << "SisterVertex::InitiateFromFile::initiate() "
		<< "Cannot open file sister." << std::endl;
      exit(EXIT_FAILURE);
    }
    size_t N=0;
    IN >> N;
    T.setNumSisterVertex(N);
    for (size_t i=0; i<N; ++i) {
      size_t v1,v2;
      IN >> v1;
      IN >> v2;
      T.setSisterVertex(i,0,v1);
      T.setSisterVertex(i,1,v2);
    }
  }

  InitiateFromDistance::
  InitiateFromDistance(std::vector<double> &paraValue, 
		       std::vector< std::vector<size_t> > 
		       &indValue )
  {
    //Do some checks on the parameters and variable indices
    //
    if( paraValue.size()!=1 ) {
      std::cerr << "SisterVertex::InitiateFromDistance::"
		<< "InitiateFromDistance() "
		<< "Uses one parameter, d_max "
		<< "(maximal distance to considered sister)." << std::endl;
      exit(0);
    }
    if( indValue.size() != 0 ) {
      std::cerr << "SisterVertex::InitiateFromDistance::"
		<< "InitiateFromDistance() "
		<< "No index given.\n";
      exit(0);
    }
    //Set the variable values
    //
    setId("SisterVertex::InitiateFromDistance");
    setParameter(paraValue);  
    setVariableIndex(indValue);
    
    //Set the parameter identities
    //
    std::vector<std::string> tmp( numParameter() );
    tmp[0] = "d_max";
    setParameterId( tmp );
  }
  
  void InitiateFromDistance::
  initiate(Tissue &T,
	   DataMatrix &cellData,
	   DataMatrix &wallData,
	   DataMatrix &vertexData,
	   DataMatrix &cellDerivs,
	   DataMatrix &wallDerivs,
	   DataMatrix &vertexDerivs)
  {
    // Go through all pairs of vertices and put close ones into the sisterVertex vector
    size_t N = T.numVertex();
    double dimension = T.vertex(0).numPosition();

    for (size_t i=0; i<N; ++i) {
      for (size_t j=i+1; j<N; ++j) {
	double distance = 0.0;
	for (size_t d=0; d<dimension; ++d) {
	  distance += (vertexData[i][d]-vertexData[j][d])*(vertexData[i][d]-vertexData[j][d]);
	}
	distance = std::sqrt(distance);
	if (distance<=parameter(0)) {//add pair to sisterVertex vector
	  T.addSisterVertex(i,j);
	}
      }
    }
  }    
  
  Spring::
  Spring(std::vector<double> &paraValue, 
	 std::vector< std::vector<size_t> > 
	 &indValue )
  {
    //Do some checks on the parameters and variable indices
    //
    if( paraValue.size()!=1 ) {
      std::cerr << "SisterVertex::Spring::"
		<< "Spring() "
		<< "Uses one parameter, k_spring." << std::endl;
      exit(0);
    }
    if( indValue.size() != 0 ) {
      std::cerr << "SisterVertex::Spring::"
		<< "Spring() "
		<< "No index given.\n";
      exit(0);
    }
    //Set the variable values
    //
    setId("SisterVertex::Spring");
    setParameter(paraValue);  
    setVariableIndex(indValue);
    
    //Set the parameter identities
    //
    std::vector<std::string> tmp( numParameter() );
    tmp[0] = "K_spring";
    setParameterId( tmp );
  }

  void Spring::
  derivs(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs )
  {
    size_t N = T.numSisterVertex();
    size_t dimension = T.vertex(0).numPosition();
    for (size_t i=0; i<N; ++i) {
      for (size_t d=0; d<dimension; ++d) {
	double derivs = -parameter(0)*(vertexData[T.sisterVertex(i,0)][d] - vertexData[T.sisterVertex(i,1)][d]);
	vertexDerivs[T.sisterVertex(i,0)][d] += derivs; 
	  vertexDerivs[T.sisterVertex(i,1)][d] -= derivs;	
      }
    }  
  }
  void Spring::
  update(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &vertexData,
	 DataMatrix &wallData,
	 DataMatrix &vertexDerivs,
	 double h)
  {
    // Maybe add a possibility to disconnect sisters that are moving too far from each other?
  }  

  CombineDerivatives::
  CombineDerivatives(std::vector<double> &paraValue, 
		     std::vector< std::vector<size_t> > 
		     &indValue )
  {
    //Do some checks on the parameters and variable indices
    //
    if( paraValue.size()!=0 ) {
      std::cerr << "SisterVertex::CombineDerivatives::"
		<< "CombineDerivatives() "
		<< "Does not use any parameters." << std::endl;
      exit(0);
    }
    if( indValue.size() != 0 ) {
      std::cerr << "SisterVertex::CombineDerivatives::"
		<< "CombineDerivatives() "
		<< "No index given.\n";
      exit(0);
    }
    //Set the variable values
    //
    setId("SisterVertex::CombineDerivatives");
    setParameter(paraValue);  
    setVariableIndex(indValue);
    
    //Set the parameter identities
    //
    std::vector<std::string> tmp( numParameter() );
    tmp[0] = "K_spring";
    setParameterId( tmp );
  }
  
  void CombineDerivatives::
  derivs(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs )
  {
    size_t N = T.numSisterVertex();
    size_t dimension = T.vertex(0).numPosition();
    for (size_t i=0; i<N; ++i) {
      for (size_t d=0; d<dimension; ++d) {
	double sum = vertexDerivs[T.sisterVertex(i,0)][d] + vertexDerivs[T.sisterVertex(i,1)][d];
	vertexDerivs[T.sisterVertex(i,0)][d] = vertexDerivs[T.sisterVertex(i,1)][d] = sum;	
      }
    }
  }
}
