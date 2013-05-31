//
// Filename     : growth.cc
// Description  : Classes describing growth updates
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : April 2006
// Revision     : $Id:$
//
#include"growth.h"
#include"baseReaction.h"

namespace WallGrowth {
  Constant::
  Constant(std::vector<double> &paraValue, 
	   std::vector< std::vector<size_t> > 
	   &indValue ) {
    
    //Do some checks on the parameters and variable indeces
    //
    if( paraValue.size()!=2 && paraValue.size() !=3) {
      std::cerr << "WallGrowth::Constant::"
		<< "Constant() "
		<< "Two or three parameters used  k_growth, linearFlag, [L_trunc]"
		<< std::endl;
      exit(EXIT_FAILURE);
    }
    if( indValue.size() != 1 || indValue[0].size() != 1 ) {
      std::cerr << "WallGrowth::Const::"
		<< "Const() "
		<< "One variable index is used for specifying resting length."
		<< std::endl;
      exit(EXIT_FAILURE);
    }
    //Set the variable values
    //
    setId("WallGrowth::Constant");
    setParameter(paraValue);  
    setVariableIndex(indValue);
    
    //Set the parameter identities
    //
    std::vector<std::string> tmp( numParameter() );
    tmp.resize( numParameter() );
    tmp[0] = "k_growth";
    tmp[1] = "linearFlag";
    if (numParameter()>2) {
      tmp[1] = "L_trunc";
    }
    setParameterId( tmp );
  }
  
  void Constant::
  derivs(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs ) {
    
    size_t numWalls = T.numWall();
    size_t lengthIndex = variableIndex(0,0);
    
    for( size_t i=0 ; i<numWalls ; ++i ) {
      double arg = parameter(0);
      if (parameter(1)==1) {//linearFlag
	arg *= wallData[i][lengthIndex];
      }
      if (numParameter()>2) {//truncated at maximal length
	arg *= (1 - wallData[i][lengthIndex]/parameter(2));
      }
      wallDerivs[i][lengthIndex] += arg;
    }
  }

  Stress::
  Stress(std::vector<double> &paraValue, 
	 std::vector< std::vector<size_t> > 
	 &indValue ) {
    
    // Do some checks on the parameters and variable indeces
    //
    if( paraValue.size()!=4 && paraValue.size()!=5 ) {
      std::cerr << "WallGrowth::Stress::"
		<< "Stress() "
		<< "Uses four or five parameters k_growth, stress_threshold "
		<< "stretch_flag (0 for stress (read from wall variable),1 for strain),"
		<< " linear_flag (0 const, 1 prop to wall length), and [L_threshold]." 
		<< std::endl;
      exit(EXIT_FAILURE);
    }
    if( paraValue[2] != 0.0 && paraValue[2] != 1.0 ) {
      std::cerr << "WallGrowth::Stress::"
		<< "Stress() "
		<< "stretch_flag parameter must be 0 (stress used) or " 
		<< "1 (stretch used)." << std::endl;
      exit(EXIT_FAILURE);
    }
    if( paraValue[3] != 0.0 && paraValue[3] != 1.0 ) {
      std::cerr << "WallGrowth::Stress::"
		<< "Stress() "
		<< "linear_flag parameter must be 0 (constant growth) or " 
		<< "1 (length dependent growth)." << std::endl;
      exit(EXIT_FAILURE);
    }
    
    if( (indValue.size()!=1 && indValue.size()!=2) || indValue[0].size() != 1 
	|| (paraValue[2]==0 && (indValue.size()!=2 || !indValue[1].size())) ) {
      std::cerr << "WallGrowth::Stress::"
		<< "Stress() "
		<< "One variable index is used (wall length index) at first "
		<< "level, and stress variable index/indices at second (if strain_flag not set)."
		<< std::endl;
      exit(EXIT_FAILURE);
    }
    //Set the variable values
    //
    setId("WallGrowth::Stress");
    setParameter(paraValue);  
    setVariableIndex(indValue);
    
    //Set the parameter identities
    //
    std::vector<std::string> tmp( numParameter() );
    tmp.resize( numParameter() );
    tmp[0] = "k_growth";
    tmp[1] = "s_threshold";
    tmp[2] = "strain_flag";
    tmp[3] = "linear_flag";
    if (numParameter()>4) {
      tmp[4] = "L_threshold";
    }
    setParameterId( tmp );
  }
  
  void Stress::
  derivs(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs ) {
    
    size_t numWalls = T.numWall();
    size_t lengthIndex = variableIndex(0,0);
    
    for( size_t i=0 ; i<numWalls ; ++i ) {
      size_t v1 = T.wall(i).vertex1()->index();
      size_t v2 = T.wall(i).vertex2()->index();
      double stress=0.0;
      if (!parameter(2)) {//Stress used, read from saved data in the wall
	for (size_t k=0; k<numVariableIndex(1); ++k)
	  stress += wallData[i][variableIndex(1,k)];
      }
      else { //Strain/stretch used
	double distance=0.0;
	for( size_t d=0 ; d<vertexData[v1].size() ; d++ )
	  distance += (vertexData[v1][d]-vertexData[v2][d])*
	    (vertexData[v1][d]-vertexData[v2][d]);
	distance = std::sqrt(distance);
	stress = (distance-wallData[i][lengthIndex]) /
	  wallData[i][lengthIndex];
      }
      if (stress > parameter(1)) {
	double growthRate = parameter(0)*(stress - parameter(1));
	if (parameter(3))
	  growthRate *= wallData[i][lengthIndex];
	if (numParameter()>4) {
	  growthRate *= (1.0 - wallData[i][lengthIndex]/parameter(4));
	}
	wallDerivs[i][lengthIndex] += growthRate;
      }
    }
  }
  
  namespace CenterTriangulation {
    Constant::
    Constant(std::vector<double> &paraValue, 
	     std::vector< std::vector<size_t> > 
	     &indValue ) {
      
      //Do some checks on the parameters and variable indeces
      //
      if( paraValue.size()!=2 && paraValue.size() !=3) {
	std::cerr << "WallGrowth::CenterTriangulation::Constant::"
		  << "Constant() "
		  << "Two or three parameters used  k_growth, linearFlag, [L_trunc]"
		  << std::endl;
	exit(EXIT_FAILURE);
      }
      if( indValue.size() != 1 || indValue[0].size() != 1 ) {
	std::cerr << "WallGrowth::CenterTriangulation::Const::"
		  << "Const() "
		  << "Start of additional Cell variable indices (center(x,y,z) "
		  << "L_1,...,L_n, n=num vertex) is given in first level." 
		  << std::endl;
	exit(EXIT_FAILURE);
      }
      //Set the variable values
      //
      setId("WallGrowth::CenterTriangulation::Constant");
      setParameter(paraValue);  
      setVariableIndex(indValue);
      
      //Set the parameter identities
      //
      std::vector<std::string> tmp( numParameter() );
      tmp.resize( numParameter() );
      tmp[0] = "k_growth";
      tmp[1] = "linearFlag";
      if (numParameter()>2) {
	tmp[1] = "L_trunc";
      }
      setParameterId( tmp );
    }
    
    void Constant::
    derivs(Tissue &T,
	   DataMatrix &cellData,
	   DataMatrix &wallData,
	   DataMatrix &vertexData,
	   DataMatrix &cellDerivs,
	   DataMatrix &wallDerivs,
	   DataMatrix &vertexDerivs ) {
      
      size_t numCells = T.numCell();
      size_t lengthIndex = variableIndex(0,0);
      size_t lengthStartIndex = lengthIndex+3;//assuming 3D
      
      for (size_t i=0; i<numCells; ++i) {
				for (size_t k=0; k<T.cell(i).numVertex(); ++k) {
					double arg = parameter(0);
	  if (parameter(1)==1) {//linearFlag (prop to length)
	    arg *= cellData[i][k+lengthStartIndex];
	  }
	  if (numParameter()>2) {//truncated at maximal length
	    arg *= (1 - cellData[i][k+lengthStartIndex]/parameter(2));
	  }
	  cellDerivs[i][k+lengthStartIndex] += arg;
	}
      }
    }
    
    Stress::
    Stress(std::vector<double> &paraValue, 
	   std::vector< std::vector<size_t> > 
	   &indValue ) {
      
      //Do some checks on the parameters and variable indeces
      //
      if( paraValue.size()!=4 ) {
	std::cerr << "WallGrowthStresscenterTriangulation::"
		  << "WallGrowthStresscenterTriangulation() "
		  << "Uses four parameters k_growth, stress_threshold "
		  << "stretch_flag and linear_flag (0 const, 1 prop to wall length)" 
		  << std::endl;
	exit(0);
      }
      if( paraValue[2] != 0.0 && paraValue[2] != 1.0 ) {
	std::cerr << "WallGrowthStresscenterTriangulation::"
		  << "WallGrowthStresscenterTriangulation() "
		  << "stretch_flag parameter must be 0 (stress used) or " 
		  << "1 (stretch used)." << std::endl;
	exit(0);
      }
      if( paraValue[2] == 0.0 ) {
	std::cerr << "WallGrowthStresscenterTriangulation::"
		  << "WallGrowthStresscenterTriangulation() "
		  << "stretch_flag parameter must be 1 (stretch used) (not implemented for" 
		  << " stress yet..." 
		  << std::endl;
	exit(0);
      }
      if( paraValue[3] != 0.0 && paraValue[3] != 1.0 ) {
	std::cerr << "WallGrowthStresscenterTriangulation::"
		  << "WallGrowthStresscenterTriangulation() "
		  << "linear_flag parameter must be 0 (constant growth) or " 
		  << "1 (length dependent growth)." << std::endl;
	exit(0);
      }
      
      if( (indValue.size()!=1 && indValue.size()!=2) || indValue[0].size() != 1 
	  || (paraValue[2]==0 && (indValue.size()!=2 || !indValue[1].size())) ) {
	std::cerr << "WallGrowthStresscenterTriangulation::"
		  << "WallGrowthStresscenterTriangulation() "
		  << "Start of additional Cell variable indices (center(x,y,z) "
		  << "L_1,...,L_n, n=num vertex) is given in first level, " 
		  << "and stress variable indices at second (if strain_flag not set)."
		  << std::endl;
	exit(0);
      }
      //Set the variable values
      //
      setId("WallGrowthStresscenterTriangulation");
      setParameter(paraValue);  
      setVariableIndex(indValue);
      
      //Set the parameter identities
      //
      std::vector<std::string> tmp( numParameter() );
      tmp.resize( numParameter() );
      tmp[0] = "k_growth";
      tmp[1] = "s_threshold";
      tmp[2] = "strain_flag";
      tmp[3] = "linear_flag";
      setParameterId( tmp );
    }
    
    void Stress::
    derivs(Tissue &T,
	   DataMatrix &cellData,
	   DataMatrix &wallData,
	   DataMatrix &vertexData,
	   DataMatrix &cellDerivs,
	   DataMatrix &wallDerivs,
	   DataMatrix &vertexDerivs ) {
      
      size_t numCells = T.numCell();
      size_t posStartIndex = variableIndex(0,0);
      size_t lengthStartIndex = posStartIndex+3;
      
      for (size_t i=0; i<numCells; ++i) {
	for (size_t k=0; k<T.cell(i).numVertex(); ++k) {
	  size_t v = T.cell(i).vertex(k)->index();
	  double stress=0.0;
	  if (!parameter(2)) {//Stress used, read from saved data in the wall
	    std::cerr << "WallGrowthStresscenterTriangulation::derivs() " << std::endl
		      << "Strain (and not stress) is the only implemented version sofar."
		      << std::endl;
	    //for (size_t k=0; k<numVariableIndex(1); ++k)
	    //stress += wallData[i][variableIndex(1,k)];
	  }
	  else { //Strain/stretch used
	    double distance=0.0;
	    for( size_t d=0 ; d<vertexData[v].size() ; d++ )
	      distance += (vertexData[v][d]-cellData[i][d+posStartIndex])*
		(vertexData[v][d]-cellData[i][d+posStartIndex]);
	    distance = std::sqrt(distance);
	    stress = (distance-cellData[i][k+lengthStartIndex]) /
	      cellData[i][k+lengthStartIndex];
	  }
	  if (stress > parameter(1)) {
	    double growthRate = parameter(0)*(stress - parameter(1));
	    if (parameter(3))
	      growthRate *= cellData[i][k+lengthStartIndex];
	    cellDerivs[i][k+lengthStartIndex] += growthRate;
	  }
	}
      }
    }
  }

  StressSpatial::
  StressSpatial(std::vector<double> &paraValue, 
		std::vector< std::vector<size_t> > 
		&indValue ) {
    
    // Do some checks on the parameters and variable indeces
    //
    if( paraValue.size()!=6 ) {
      std::cerr << "WallGrowth::StressSpatial::"
		<< "StressSpatial() "
		<< "Uses six parameters k_growth, stress(stretch)_threshold "
		<< "K_hill n_Hill "
		<< "stretch_flag and linear_flag" << std::endl;
      exit(0);
    }
    if( paraValue[4] != 0.0 && paraValue[4] != 1.0 ) {
      std::cerr << "WallGrowth::StressSpatial::"
		<< "StressSpatial() "
		<< "stretch_flag parameter must be 0 (stress used) or " 
		<< "1 (stretch used)." << std::endl;
      exit(0);
    }
    if( paraValue[5] != 0.0 && paraValue[5] != 1.0 ) {
      std::cerr << "WallGrowth::StressSpatial::"
		<< "StressSpatial() "
		<< "linear_flag parameter must be 0 (constant growth) or " 
		<< "1 (length dependent growth)." << std::endl;
      exit(0);
    }
    
    if( indValue.size() != 2 || indValue[0].size() != 2 ) {
      std::cerr << "WallGrowth::StressSpatial::"
		<< "StressSpatial() "
		<< "Two variable index is used (wall length,spatial coordinate) at first "
		<< "level, and force variable index at second."
		<< std::endl;
      exit(0);
    }
    // Set the variable values
    //
    setId("WallGrowth::StressSpatial");
    setParameter(paraValue);  
    setVariableIndex(indValue);
    Kpow_=std::pow(paraValue[2],paraValue[3]);
    
    // Set the parameter identities
    //
    std::vector<std::string> tmp( numParameter() );
    tmp.resize( numParameter() );
    tmp[0] = "k_growth";
    tmp[1] = "stress_threshold";
    tmp[2] = "K_Hill";
    tmp[3] = "n_Hill";
    tmp[4] = "stretch_flag";
    tmp[5] = "linear_flag";
    setParameterId( tmp );
  }
  
  void StressSpatial::
  derivs(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs ) {
    
    size_t numWalls = T.numWall();
    size_t lengthIndex = variableIndex(0,0);
    size_t dimension = vertexData[0].size();
    
    // Prepare spatial factor
    size_t sI=variableIndex(0,1);
    assert (sI<vertexData[0].size());
    size_t numVertices = vertexData.size();
    double sMax= vertexData[0][sI];
    size_t maxI=0;
    for (size_t i=1; i<numVertices; ++i)
      if (vertexData[i][sI]>sMax) {
	sMax=vertexData[i][sI];
	maxI=i;
      }
    std::vector<double> maxPos(dimension);
    for (size_t d=0; d<dimension; ++d)
      maxPos[d] = vertexData[maxI][d];
    
    for( size_t i=0 ; i<numWalls ; ++i ) {
      size_t v1 = T.wall(i).vertex1()->index();
      size_t v2 = T.wall(i).vertex2()->index();
      double stress=0.0;
      if (!parameter(4)) {
	for (size_t k=0; k<numVariableIndex(1); ++k)
	  stress += wallData[i][variableIndex(1,k)];
      }
      else {
	double distance=0.0;
	for( size_t d=0 ; d<dimension ; ++d )
	  distance += (vertexData[v1][d]-vertexData[v2][d])*
	    (vertexData[v1][d]-vertexData[v2][d]);
	distance = std::sqrt(distance);
	stress = (distance-wallData[i][lengthIndex]) /
	  wallData[i][lengthIndex];
      }
      if (stress > parameter(1)) {
	// Calculate spatial factor
	double maxDistance = 0.0;
	for (size_t d=0; d<dimension; ++d) {
	  double pos = 0.5*(vertexData[v1][d]+vertexData[v2][d]);
	  maxDistance += (maxPos[d]-pos)*(maxPos[d]-pos);
	}
	maxDistance = std::sqrt(maxDistance);
	double spatialFactor = Kpow_/(Kpow_+std::pow(maxDistance,parameter(3)));
	
	double growthRate = parameter(0)*(stress - parameter(1))*spatialFactor;
		
	if (parameter(5))
	  growthRate *= wallData[i][lengthIndex];
	wallDerivs[i][lengthIndex] += growthRate;
      }
    }
  }
  
  StressSpatialSingle::
  StressSpatialSingle(std::vector<double> &paraValue, 
		      std::vector< std::vector<size_t> > 
		      &indValue ) {
    
    // Do some checks on the parameters and variable indeces
    //
    if( paraValue.size()!=6 ) {
      std::cerr << "WallGrowth::StressSpatialSingle::"
		<< "StressSpatialSingle() "
		<< "Uses six parameters k_growth, stress(stretch)_threshold "
		<< "K_hill n_Hill "
		<< "stretch_flag and linear_flag" << std::endl;
      exit(0);
    }
    if( paraValue[4] != 0.0 && paraValue[4] != 1.0 ) {
      std::cerr << "WallGrowth::StressSpatialSingle::"
		<< "StressSpatialSingle() "
		<< "stretch_flag parameter must be 0 (stress used) or " 
		<< "1 (stretch used)." << std::endl;
      exit(0);
    }
    if( paraValue[5] != 0.0 && paraValue[5] != 1.0 ) {
      std::cerr << "WallGrowth::StressSpatialSingle::"
		<< "StressSpatialSingle() "
		<< "linear_flag parameter must be 0 (constant growth) or " 
		<< "1 (length dependent growth)." << std::endl;
      exit(0);
    }
    
    if( indValue.size() != 2 || indValue[0].size() != 2 ) {
      std::cerr << "WallGrowth::StressSpatialSingle::"
		<< "StressSpatialSingle() "
		<< "Two variable index is used (wall length,spatial coordinate) at first "
		<< "level, and force variable index at second."
		<< std::endl;
      exit(0);
    }
    // Set the variable values
    //
    setId("WallGrowth::StressSpatialSingle");
    setParameter(paraValue);  
    setVariableIndex(indValue);
    Kpow_=std::pow(paraValue[2],paraValue[3]);
    
    // Set the parameter identities
    //
    std::vector<std::string> tmp( numParameter() );
    tmp.resize( numParameter() );
    tmp[0] = "k_growth";
    tmp[1] = "stress_threshold";
    tmp[2] = "K_Hill";
    tmp[3] = "n_Hill";
    tmp[4] = "stretch_flag";
    tmp[5] = "linear_flag";
    setParameterId( tmp );
  }
  
  void StressSpatialSingle::
  derivs(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs ) {
    
    size_t numWalls = T.numWall();
    size_t lengthIndex = variableIndex(0,0);
    size_t dimension = vertexData[0].size();
    
    // Prepare spatial factor
    size_t sI=variableIndex(0,1);
    assert (sI<vertexData[0].size());
    size_t numVertices = vertexData.size();
    double sMax= vertexData[0][sI];
    size_t maxI=0;
    for (size_t i=1; i<numVertices; ++i)
      if (vertexData[i][sI]>sMax) {
	sMax=vertexData[i][sI];
	maxI=i;
      }
    
    for( size_t i=0 ; i<numWalls ; ++i ) {
      size_t v1 = T.wall(i).vertex1()->index();
      size_t v2 = T.wall(i).vertex2()->index();
      double stress=0.0;
      if (!parameter(4)) {
	for (size_t k=0; k<numVariableIndex(1); ++k)
	  stress += wallData[i][variableIndex(1,k)];
      }
      else {
	double distance=0.0;
	for( size_t d=0 ; d<dimension ; ++d )
	  distance += (vertexData[v1][d]-vertexData[v2][d])*
	    (vertexData[v1][d]-vertexData[v2][d]);
	distance = std::sqrt(distance);
	stress = (distance-wallData[i][lengthIndex]) /
	  wallData[i][lengthIndex];
      }
      if (stress > parameter(1)) {
	// Calculate spatial factor
	double maxDistance = sMax - 0.5*(vertexData[v1][sI]+vertexData[v2][sI]);;
	double spatialFactor = Kpow_/(Kpow_+std::pow(maxDistance,parameter(3)));
	
	double growthRate = parameter(0)*(stress - parameter(1))*spatialFactor;
	
	if (parameter(5))
	  growthRate *= wallData[i][lengthIndex];
	wallDerivs[i][lengthIndex] += growthRate;
      }
    }
  }
  
  StressConcentrationHill::
  StressConcentrationHill(std::vector<double> &paraValue, 
			  std::vector< std::vector<size_t> > 
			  &indValue ) {
    
    // Do some checks on the parameters and variable indeces
    //
    if( paraValue.size()!=7 ) {
      std::cerr << "WallGrowth::StressConcentrationHill::"
		<< "StressConcentrationHill() "
		<< "Uses seven parameters k_growthConst, k_growthHill, K_Hill, n_Hill,"
		<< " stretch_threshold stretch_flag and linear_flag" << std::endl;
      exit(0);
    }
    if( paraValue[5] != 0.0 && paraValue[5] != 1.0 ) {
      std::cerr << "WallGrowth::StressConcentrationHill::"
		<< "StressConcentrationHill() "
		<< "stretch_flag parameter must be 0 (stress used) or " 
		<< "1 (stretch used)." << std::endl;
      exit(0);
    }
    if( paraValue[6] != 0.0 && paraValue[6] != 1.0 ) {
      std::cerr << "WallGrowth::StressConcentrationHill::"
		<< "StressConcentrationHill() "
		<< "linear_flag parameter must be 0 (constant growth) or " 
		<< "1 (length dependent growth)." << std::endl;
      exit(0);
    }
    
    if( indValue.size() != 2 || indValue[0].size() != 2 ) {
      std::cerr << "WallGrowth::StressConcentrationHill::"
		<< "StressConcentrationHill() "
		<< "wall length index and concentration index at first "
		<< "level, and spring constant variable indices at second"
		<< std::endl;
      exit(0);
    }
    // Set the variable values
    //
    setId("WallGrowth::StressConcentrationHill");
    setParameter(paraValue);  
    setVariableIndex(indValue);
    
    // Set the parameter identities
    //
    std::vector<std::string> tmp( numParameter() );
    tmp.resize( numParameter() );
    tmp[0] = "k_growth";
    tmp[1] = "stress_threshold";
    tmp[2] = "stretch_flag";
    tmp[3] = "linear_flag";
    setParameterId( tmp );
  }

  void StressConcentrationHill::
  derivs(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs ) {
    
    size_t numWalls = T.numWall();
    size_t lengthIndex = variableIndex(0,0);
    size_t concIndex = variableIndex(0,1);
    double Kpow = std::pow(parameter(2),parameter(3));
    
    for( size_t i=0 ; i<numWalls ; ++i ) {
      size_t v1 = T.wall(i).vertex1()->index();
      size_t v2 = T.wall(i).vertex2()->index();
      double stress=0.0;
      if (!parameter(5)) {
	for (size_t k=0; k<numVariableIndex(1); ++k)
	  stress += wallData[i][variableIndex(1,k)];
      }
      else {
	double distance=0.0;
	for( size_t d=0 ; d<vertexData[v1].size() ; d++ )
	  distance += (vertexData[v1][d]-vertexData[v2][d])*
	    (vertexData[v1][d]-vertexData[v2][d]);
	distance = std::sqrt(distance);
	stress = (distance-wallData[i][lengthIndex]) /
	  wallData[i][lengthIndex];
      }
      if (stress > parameter(4)) {
	// Get the Hill factor from the two cells
	double hillFactor=0.0;
	if (T.wall(i).cell1() != T.background()) {
	  double concpow = std::pow(cellData[T.wall(i).cell1()->index()][concIndex],parameter(3));
	  hillFactor += concpow/(Kpow+concpow); 
	}
	if (T.wall(i).cell2() != T.background()) {
	  double concpow = std::pow(cellData[T.wall(i).cell2()->index()][concIndex],parameter(3));
	  hillFactor += concpow/(Kpow+concpow); 
	}
	double growthRate = (parameter(0)+hillFactor*parameter(1))*(stress - parameter(4));
	if (parameter(6))
	  growthRate *= wallData[i][lengthIndex];
	wallDerivs[i][lengthIndex] += growthRate;
      }
    }
  }
  
  ConstantStressEpidermalAsymmetric::
  ConstantStressEpidermalAsymmetric(std::vector<double> &paraValue, 
				    std::vector< std::vector<size_t> > 
				    &indValue ) {
    
    // Do some checks on the parameters and variable indeces
    //
    if( paraValue.size()!=2 ) {
      std::cerr << "WallGrowth::ConstantStressEpidermalAsymmetric::"
		<< "ConstantStressEpidermalAsymmetric() "
		<< "Uses two parameters k_growth and frac_epi\n";
      exit(0);
    }  
    if( indValue.size() != 1 || indValue[0].size() != 1 ) {
      std::cerr << "WallGrowth::ConstantStressEpidermalAsymmetric::"
		<< "ConstantStressEpidermalAsymmetric() "
		<< "One variable index is used.\n";
      exit(0);
    }
    //Set the variable values
    //
    setId("WallGrowth::ConstantStressEpidermalAsymmetric");
    setParameter(paraValue);  
    setVariableIndex(indValue);
    
    //Set the parameter identities
    //
    std::vector<std::string> tmp( numParameter() );
    tmp.resize( numParameter() );
    tmp[0] = "k_growth";
    tmp[0] = "frac_epi";
    setParameterId( tmp );
  }
  
  void ConstantStressEpidermalAsymmetric::
  derivs(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs ) {
    
    size_t numWalls = T.numWall();
    size_t lengthIndex = variableIndex(0,0);
    
    for( size_t i=0 ; i<numWalls ; ++i ) {
      size_t v1 = T.wall(i).vertex1()->index();
      size_t v2 = T.wall(i).vertex2()->index();
      double kGrowth = parameter(0);
      if( T.wall(i).cell1() == T.background() || 
	  T.wall(i).cell2() == T.background() )
	kGrowth *= parameter(1);
      double distance=0.0;
      for( size_t d=0 ; d<vertexData[v1].size() ; d++ )
	distance += (vertexData[v1][d]-vertexData[v2][d])*
	  (vertexData[v1][d]-vertexData[v2][d]);
      distance = std::sqrt(distance);
      if( distance>wallData[i][lengthIndex] )
	wallDerivs[i][lengthIndex] += kGrowth*
	  (distance-wallData[i][lengthIndex]);
    }
  }
  
  Force::Force(std::vector<double> &paraValue,
	       std::vector< std::vector<size_t> > &indValue)
  {
    if (paraValue.size() != 2) {
      std::cerr << "WallGrowth::Force::Force() "
		<< "Uses two parameters: k_growth and Force_threshold" << std::endl;
      exit(EXIT_FAILURE);
    }
    
    if (indValue.size() != 2 || indValue[0].size() != 1) {
      std::cerr << "WallGrowth::Force::Force() "
		<< "Wall length index must be given in first level.\n"
		<< "Wall force index/indices must be given in second level.\n";
      exit(EXIT_FAILURE);
    }
    
    setId("VertexFromWallSpringExperimental");
    setParameter(paraValue);
    setVariableIndex(indValue);
    
    std::vector<std::string> tmp(numParameter());
    tmp[0] = "k_L";
    tmp[1] = "phi";
    
    setParameterId(tmp);
  }

  void Force::
  derivs(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs)
  {
    for (size_t i = 0; i < T.numWall(); ++i) {
      size_t vertex1Index = T.wall(i).vertex1()->index();
      size_t vertex2Index = T.wall(i).vertex2()->index();
      size_t dimensions = vertexData[vertex1Index].size();
      
      double distance = 0.0;
      for (size_t d = 0; d < dimensions; ++d) {
	distance += (vertexData[vertex1Index][d] - vertexData[vertex2Index][d])
	  * (vertexData[vertex1Index][d] - vertexData[vertex2Index][d]);
      }
      distance = std::sqrt(distance);
      
      double F = 0.0;
      for (size_t j = 0; j < numVariableIndex(1); ++j)
	F += wallData[T.wall(i).index()][variableIndex(1, j)];
      
      double arg = F - parameter(1);
      if (arg > 0)
	wallDerivs[T.wall(i).index()][variableIndex(0, 0)] 
	  += parameter(0) * arg; 
      //*wallData[T.wall(i).index()][variableIndex(0, 0)];
    }
  }
}

MoveVertexRadially::
MoveVertexRadially(std::vector<double> &paraValue, 
		   std::vector< std::vector<size_t> > 
		   &indValue ) {
  
  // Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=2 || ( paraValue[1]!=0 && paraValue[1]!=1) ) {
    std::cerr << "MoveVertexRadially::"
	      << "MoveVertexRadially() "
	      << "Uses two parameters k_growth and r_pow (0,1)\n";
    exit(0);
  }  
  if( indValue.size() != 0 ) {
    std::cerr << "MoveVertexRadially::"
	      << "MoveVertexRadially() "
	      << "No variable index is used.\n";
    exit(0);
  }
  // Set the variable values
  //
  setId("MoveVertexRadially");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  // Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "k_growth";
  tmp[0] = "r_pow";
  setParameterId( tmp );
}

void MoveVertexRadially::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) {
  
  size_t numVertices = T.numVertex();
  size_t dimension=vertexData[0].size();
  
  for( size_t i=0 ; i<numVertices ; ++i ) {
    double fac=parameter(0);
    if( parameter(1)==0.0 ) {
      double r=0.0;
      for( size_t d=0 ; d<dimension ; ++d )
	r += vertexData[i][d]*vertexData[i][d];
      if( r>0.0 )
	r = std::sqrt(r);
      if( r>0.0 )
	fac /= r;
      else
	fac=0.0;
    }
    for( size_t d=0 ; d<dimension ; ++d )
      vertexDerivs[i][d] += fac*vertexData[i][d];
  }
}

MoveVertexRadiallycenterTriangulation::
MoveVertexRadiallycenterTriangulation(std::vector<double> &paraValue, 
				      std::vector< std::vector<size_t> > 
				      &indValue ) {
  
  // Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=2 || ( paraValue[1]!=0 && paraValue[1]!=1) ) {
    std::cerr << "MoveVertexRadiallycenterTriangulation::"
	      << "MoveVertexRadiallycenterTriangulation() "
	      << "Uses two parameters k_growth and r_pow (0,1)\n";
    exit(0);
  }  
  if( indValue.size() != 1 || indValue[0].size() != 1 ) {
    std::cerr << "MoveVertexRadiallycenterTriangulation::"
	      << "MoveVertexRadiallycenterTriangulation() " << std::endl
	      << "Start of additional Cell variable indices (center(x,y,z) "
	      << "L_1,...,L_n, n=num vertex) is given in first level." 
	      << std::endl;
    exit(0);
  }
  // Set the variable values
  //
  setId("MoveVertexRadiallycenterTriangulation");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  // Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "k_growth";
  tmp[0] = "r_pow";
  setParameterId( tmp );
}

void MoveVertexRadiallycenterTriangulation::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) {
  
  size_t numVertices = T.numVertex();
  size_t numCells = T.numCell();
  size_t dimension=vertexData[0].size();
  
  // Move vertices
  for( size_t i=0 ; i<numVertices ; ++i ) {
    double fac=parameter(0);
    if( parameter(1)==0.0 ) {
      double r=0.0;
      for( size_t d=0 ; d<dimension ; ++d )
	r += vertexData[i][d]*vertexData[i][d];
      if( r>0.0 )
	r = std::sqrt(r);
      if( r>0.0 )
	fac /= r;
      else
	fac=0.0;
    }
    for( size_t d=0 ; d<dimension ; ++d )
      vertexDerivs[i][d] += fac*vertexData[i][d];
  }
  // Move vertices defined in cell centers
  for( size_t i=0 ; i<numCells ; ++i ) {
    double fac=parameter(0);
    if( parameter(1)==0.0 ) {
      double r=0.0;
      for( size_t d=variableIndex(0,0) ; d<variableIndex(0,0)+dimension ; ++d )
	r += cellData[i][d]*cellData[i][d];
      if( r>0.0 )
	r = std::sqrt(r);
      if( r>0.0 )
	fac /= r;
      else
	fac=0.0;
    }
    for( size_t d=variableIndex(0,0) ; d<variableIndex(0,0)+dimension ; ++d )
      cellDerivs[i][d] += fac*cellData[i][d];
  }
}

MoveVertexSphereCylinder::
MoveVertexSphereCylinder(std::vector<double> &paraValue, 
			 std::vector< std::vector<size_t> > 
			 &indValue ) {
  
  // Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=2 || ( paraValue[1]!=0 && paraValue[1]!=1) ) {
    std::cerr << "MoveVertexSphereCylinder::"
	      << "MoveVertexSphereCylinder() "
	      << "Uses two parameters k_growth and r_pow (0,1)\n";
    exit(0);
  }  
  if( indValue.size() != 0 ) {
    std::cerr << "MoveVertexSphereCylinder::"
	      << "MoveVertexSphereCylinder() "
	      << "No variable index is used.\n";
    exit(0);
  }
  // Set the variable values
  //
  setId("MoveVertexSphereCylinder");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  // Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "k_growth";
  tmp[0] = "r_pow";
  setParameterId( tmp );
}

void MoveVertexSphereCylinder::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) {
  
  size_t numVertices = T.numVertex();
  if (vertexData[0].size()!=3) {
    std::cerr << "MoveVertexSphereCylinder:: Only works for 3 dimensions." << std::endl;
    exit(-1);
  }
  size_t xI=0;
  size_t yI=1;
  size_t zI=2;
 
  for( size_t i=0 ; i<numVertices ; ++i ) {
    if (vertexData[i][zI]<0.0) { // on cylinder
      if( parameter(1)==0.0 ) {
	vertexDerivs[i][zI] -= parameter(0);
      }
      else {
	double r = std::sqrt(vertexData[i][xI]*vertexData[i][xI]+
			     vertexData[i][yI]*vertexData[i][yI]);
	vertexDerivs[i][zI] -= parameter(0)*(3.14159265*0.5*r-vertexData[i][zI]);
      }
    }
    else { // on half sphere
      double r = std::sqrt(vertexData[i][xI]*vertexData[i][xI]+
			   vertexData[i][yI]*vertexData[i][yI]+
			   vertexData[i][zI]*vertexData[i][zI]);
      double rPrime = std::sqrt(vertexData[i][xI]*vertexData[i][xI]+
				vertexData[i][yI]*vertexData[i][yI]);
      double theta = std::asin(rPrime/r);

      double fac=parameter(0)*theta;
      if (parameter(0)==1) {
	fac *= r;
      }
      vertexDerivs[i][xI] += fac*vertexData[i][xI]*vertexData[i][zI]/rPrime;
      vertexDerivs[i][yI] += fac*vertexData[i][yI]*vertexData[i][zI]/rPrime;
      vertexDerivs[i][zI] -= fac*rPrime;
    }
  }
}

WaterVolumeFromTurgor::
WaterVolumeFromTurgor(std::vector<double> &paraValue,
		      std::vector< std::vector<size_t> > &indValue)
{
  if (paraValue.size() != 5) {
    std::cerr << "WaterVolumeFromTurgor::WaterVolumeFromTurgor() "
	      << "Uses five parameters: k_p, P_max, k_pp and "
	      << "denyShrink_flag allowNegTurgor_flag." << std::endl;
    exit(EXIT_FAILURE);
  }
  
  if (indValue.size() < 1 || indValue.size() > 2
      || indValue[0].size() != 1 
      || (indValue.size()==2 && indValue[1].size() != 1 ) ) {
    std::cerr << "WaterVolumeFromTurgor::WaterVolumeFromTurgor() "
	      << "Water volume index must be given in "
	      << "first level.\n"
	      << "Optionally index for saving the turgor pressure can be"
	      << " given at second level." << std::endl; 		
    exit(EXIT_FAILURE);
  }
  
  setId("WaterVolumeFromTurgor");
  setParameter(paraValue);
  setVariableIndex(indValue);
  
  std::vector<std::string> tmp(numParameter());
  tmp[0] = "k_p";
  tmp[1] = "P_max";
  tmp[2] = "k_pp";
  tmp[3] = "denyShrink_flag";
  tmp[4] = "allowNegTurgor_flag";
  setParameterId(tmp);
}

void WaterVolumeFromTurgor::
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
    double totalLength = 0.0;
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
      totalLength += distance;
      
      // Old turgor measure from wall extensions
      //for (size_t j = 0; j < numVariableIndex(1); ++j)
      //P += wallData[cell.wall(i)->index()][variableIndex(1, j)]/distance;
    }
    
    // Calculate turgor measure from volume and 'water volume'
    // P ~ p_2(V_w-V)/V
    double cellVolume =cell.calculateVolume(vertexData);
    P = (cellData[cell.index()][variableIndex(0,0)]-cellVolume) / cellVolume;
    if (P<0.0 && !parameter(4))
      P=0.0;
    
    P *= parameter(2);
    
    if (numVariableIndexLevel()==2)
      cellData[n][variableIndex(1,0)]=P;
    
    if( !parameter(3) || parameter(1)-P>0.0 )
      cellDerivs[cell.index()][variableIndex(0,0)] += 
	parameter(0) * (parameter(1) - P) * totalLength;
  }
}

DilutionFromVertexDerivs::
DilutionFromVertexDerivs(std::vector<double> &paraValue,
			 std::vector< std::vector<size_t> > &indValue)
{
  if (paraValue.size()) {
    std::cerr << "DilutionFromVertexDerivs::DilutionFromVertexDerivs() "
	      << "Uses no parameters." << std::endl;
    exit(EXIT_FAILURE);
  }
  
  if (indValue.size() != 1 || indValue[0].size() < 1) {
    std::cerr << "DilutionFromVertexDerivs::DilutionFromVertexDerivs() "
	      << "List of concentration variable index must be given in "
	      << "first level." << std::endl;
    exit(EXIT_FAILURE);
  }
  
  setId("DilutionFromVertexDerivs");
  setParameter(paraValue);
  setVariableIndex(indValue);
  
  std::vector<std::string> tmp(numParameter());
  setParameterId(tmp);
}

void DilutionFromVertexDerivs::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs)
{
  size_t dimension;
  dimension = vertexData[0].size();
  assert(dimension==2);
  
  for (size_t n = 0; n < T.numCell(); ++n) {
    Cell cell = T.cell(n);
    double area = cell.calculateVolume(vertexData,1);
    
    double areaDerivs=0.0;
    for( size_t k=0 ; k<cell.numVertex() ; ++k ) {
      size_t vI = cell.vertex(k)->index();
      size_t vIPlus = cell.vertex((k+1)%(cell.numVertex()))->index();
      areaDerivs += vertexData[vIPlus][1]*vertexDerivs[vI][0] - 
	vertexData[vI][1]*vertexDerivs[vIPlus][0] -
	vertexData[vIPlus][0]*vertexDerivs[vI][1] +
	vertexData[vI][0]*vertexDerivs[vIPlus][1];
    }
    
    double fac = areaDerivs/area;
    for (size_t k=0; k<numVariableIndex(0); ++k)
      cellDerivs[n][variableIndex(0,k)] -= cellData[n][variableIndex(0,k)]*
	fac;
  }
}
