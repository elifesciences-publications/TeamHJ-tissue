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

  void InitiateFromFile::
  derivs(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs )
  {
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
	  std::cerr << "SisterVertex::InitiateFromDistance::initiate() added sisters "
		    << i << " " << j << std::endl;
	}
      }
    }
  }    

  void InitiateFromDistance::
  derivs(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs )
  {
  }
  
  Spring::
  Spring(std::vector<double> &paraValue, 
	 std::vector< std::vector<size_t> > 
	 &indValue )
  {
    //Do some checks on the parameters and variable indices
    //
    if (paraValue.size()!=1 && paraValue.size()!=2) {
      std::cerr << "SisterVertex::Spring::"
		<< "Spring() "
		<< "Uses one or two parameters, k_spring and [lengthFractionBreak]." 
		<< std::endl;
      exit(EXIT_FAILURE);
    }
    if( indValue.size() != 0 ) {
      std::cerr << "SisterVertex::Spring::"
		<< "Spring() "
		<< "No index given.\n";
      exit(EXIT_FAILURE);
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
    if (numParameter()==2)
      tmp[1] = "BreakLength";
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
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 double h)
  {
    if (numParameter()==2) {
      std::vector<size_t> remove;
      size_t N=T.numSisterVertex();
      size_t dimension = vertexData[0].size();
      // Mark sister vertex pairs for removal
      for (size_t i=0; i<N; ++i) {
	double distance = 0.0;
	for (size_t d=0; d<dimension; ++d) {
	  distance += (vertexData[T.sisterVertex(i,0)][d] - vertexData[T.sisterVertex(i,1)][d])*
	    (vertexData[T.sisterVertex(i,0)][d] - vertexData[T.sisterVertex(i,1)][d]);
	}
	distance = std::sqrt(distance);
	if (distance>parameter(1)) {
	  remove.push_back(i);
	  //std::cerr << "SisterVertex::Spring::update() Marked for removal sisterVertex between "
	  //	    << T.sisterVertex(i,0) << " " << T.sisterVertex(i,1) << " " << i << std::endl;

	}
      }
      // Remove sistervertex pairs marked
      if (remove.size()) {
	for (int k=remove.size()-1; (k>=0 && k<remove.size()); --k) {
	  size_t i = remove[k];
	  size_t NN = T.numSisterVertex()-1;
	  // If last element, remove it
	  if (i==NN) {
	    //std::cerr << "SisterVertex::Spring::update() Removed (last) sisterVertex between "
	    //	      << T.sisterVertex(i,0) << " " << T.sisterVertex(i,1) << std::endl;
	    T.sisterVertexPopBack();	    
	  }
	  else {
	    // If not last element, copy last element to i and then remove last element
	    //std::cerr << "SisterVertex::Spring::update() Remove sisterVertex between "
	    //	      << T.sisterVertex(i,0) << " " << T.sisterVertex(i,1) << std::endl;
	    T.setSisterVertex(i,0,T.sisterVertex(NN,0));
	    T.setSisterVertex(i,1,T.sisterVertex(NN,1));
	    T.sisterVertexPopBack();	    
	  }
	}
      }
    }
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
    //std::vector<std::string> tmp( numParameter() );
    //setParameterId( tmp );
  }
  
  void CombineDerivatives::
  initiate(Tissue &T,
                DataMatrix &cellData,
                DataMatrix &wallData,
                DataMatrix &vertexData,
                DataMatrix &cellDerivs,
                DataMatrix &wallDerivs,
                DataMatrix &vertexDerivs)
  {
    // Go through all pairs of sistervertices and put the pairs with 
    // common nodes together in the private vector sisters
    size_t N=T.numSisterVertex();

    std::vector<std::vector<double>> tmpsisters;
    tmpsisters.resize(N);
    for(size_t i=0; i<N; ++i){
      tmpsisters[i].resize(3);
      tmpsisters[i][0]=T.sisterVertex(i,0);
      tmpsisters[i][1]=T.sisterVertex(i,1);
      tmpsisters[i][2]=0;
    }
    
    size_t counter=1;
    for(size_t i=0; i<N; ++i)
      if(tmpsisters[i][2]==0){
        tmpsisters[i][2]=counter;
        sisters.resize(counter);
        sisters[counter-1].push_back(tmpsisters[i][0]);
        sisters[counter-1].push_back(tmpsisters[i][1]);
        for(size_t j=i+1; j<N; j++)
          if(tmpsisters[j][2]==0){
            size_t M=sisters[counter-1].size();
            size_t ww=0;
            for(size_t k=0; k<M; ++k){
              if(tmpsisters[j][0]==sisters[counter-1][k])
                ww+=1;
              if(tmpsisters[j][1]==sisters[counter-1][k])
                ww+=2;
            }
            if(ww==1)
              sisters[counter-1].push_back(tmpsisters[j][1]);
            if(ww==2)
              sisters[counter-1].push_back(tmpsisters[j][0]);
            if(ww!=0)
              tmpsisters[j][2]=counter;
          }
        counter++;
      }
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
    size_t N = sisters.size();
    size_t dimension = T.vertex(0).numPosition();
    
    for (size_t i=0; i<N; ++i) {
      size_t M=sisters[i].size();
      for (size_t d=0; d<dimension; ++d) {
        double sum=0;
        for(size_t j=0; j<M; ++j) 
          sum += vertexDerivs[sisters[i][j]][d];
        
        for(size_t j=0; j<M; ++j)
          vertexDerivs[sisters[i][j]][d] = sum;	  
      }
    }
    
  }
  
}
