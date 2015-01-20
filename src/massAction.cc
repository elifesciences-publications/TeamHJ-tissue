//
// Filename     : massAction.cc
// Description  : Classes describing mass action reactions
// Author(s)    : Laura Brown (laura.brown@slcu.cam.ac.uk)
// Created      : August 2014
// Revision     : $Id:$
//
#include "massAction.h"
#include "baseReaction.h"


namespace MassAction {
  
  General::
  General(std::vector<double> &paraValue, 
	  std::vector< std::vector<size_t> > &indValue ) 
  {  
    //
    // Do some checks on the parameters and variable indeces
    //
    if( paraValue.size()!=1 ) {
      std::cerr << "MassAction::General::General() "
		<< "Uses only one parameter k_f" << std::endl;
      exit(0);
    }
    if( indValue.size() !=2 ) {
      std::cerr << "MassAction::General::General() "
		<< "Two levels of variable indices are used." << std::endl
		<< "One for reactants and one for products." << std::endl;
      exit(0);
    }
    if( indValue[0].size()<1 ) {
      std::cerr << "MassAction::General::General() "
		<< "If no reactants given, there will be no reaction..." << std::endl;
      exit(0);
    }
    //
    // Set the variable values
    //
    setId("MassAction::General");
    setParameter(paraValue);  
    setVariableIndex(indValue);
    //
    // Set the parameter identities
    //
    std::vector<std::string> tmp( numParameter() );
    tmp.resize( numParameter() );
    tmp[0] = "k_f";
    setParameterId( tmp );
  }
  
  void General::
  derivs(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs ) 
  {  
    if( numVariableIndex(0) ) {
      for( size_t cellIndex=0 ; cellIndex<cellData.size() ; cellIndex++ ) {
	double rate = parameter(0);
	for( size_t i=0 ; i< numVariableIndex(0) ; i++ )
	  rate *= cellData[cellIndex][variableIndex(0,i)];
    
	if (rate>0.0) {
	  for( size_t i=0 ; i< numVariableIndex(0) ; ++i )
	    cellDerivs[cellIndex][variableIndex(0,i)] -= rate;
	  
	  for( size_t i=0 ; i< numVariableIndex(1) ; ++i )
	    cellDerivs[cellIndex][variableIndex(1,i)] += rate;
	}
      }
    }
    else { //No reaction defined...
      return;
    }
  }

  void General::
  derivsWithAbs(Tissue &T,
		DataMatrix &cellData,
		DataMatrix &wallData,
		DataMatrix &vertexData,
		DataMatrix &cellDerivs,
		DataMatrix &wallDerivs,
		DataMatrix &vertexDerivs,
		DataMatrix &sdydtCell,
		DataMatrix &sdydtWall,
		DataMatrix &sdydtVertex ) 
  {  
    if( numVariableIndex(0) ) {
      for( size_t cellIndex=0 ; cellIndex<cellData.size() ; cellIndex++ ) {
	double rate = parameter(0);
	for( size_t i=0 ; i< numVariableIndex(0) ; i++ )
	  rate *= cellData[cellIndex][variableIndex(0,i)];
	
	if (rate>0.0) {
	  for( size_t i=0 ; i< numVariableIndex(0) ; ++i ) {
	    cellDerivs[cellIndex][variableIndex(0,i)] -= rate;
	    sdydtCell[cellIndex][variableIndex(0,i)] += rate;
	  }
	  
	  for( size_t i=0 ; i< numVariableIndex(1) ; ++i ) {
	    cellDerivs[cellIndex][variableIndex(1,i)] += rate;
	    sdydtCell[cellIndex][variableIndex(1,i)] += rate;
	  }
	}
      }
    }
    else { //No reaction defined...
      return;
    }
  }

  OneToTwo::
  OneToTwo(std::vector<double> &paraValue, 
	   std::vector< std::vector<size_t> > 
	   &indValue ) {
    
    //Do some checks on the parameters and variable indeces
    //
    if( paraValue.size()!=1 ) {
      std::cerr << "OneToTwo::"
		<< "OneToTwo() "
		<< "1 parameters used (see MembraneCycling.h)\n";
      exit(0);
    }
    if( indValue.size() != 2 || indValue[0].size() != 3 || indValue[1].size() != 0 ) {
      std::cerr << "OneToTwo::"
		<< "OneToTwo() "
		<< "Three cell variable indices (first row) and Zero wall variable"
		<< " indices are used." << std::endl;
      exit(0);
    }
    //Set the variable values
    //
    setId("OneToTwo");
    setParameter(paraValue);  
    setVariableIndex(indValue);
    
    //Set the parameter identities
    //
    std::vector<std::string> tmp( numParameter() );
    tmp.resize( numParameter() );
    tmp[0] = "k1";
    setParameterId( tmp );
  }
  
  void OneToTwo::
  derivs(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs ) 
  {  
    size_t numCells = T.numCell();
    size_t rI = variableIndex(0,0);//reactant
    size_t p1I = variableIndex(0,1);//product1
    size_t p2I = variableIndex(0,2);//product2
    
    assert( rI<cellData[0].size() &&
	    p1I<cellData[0].size() &&
	    p2I<cellData[0].size() );
    
    for (size_t i=0; i<numCells; ++i) {
      
      double fac = parameter(0)*cellData[i][rI];
      
      cellDerivs[i][p1I] += fac;
      cellDerivs[i][p2I] += fac;
      cellDerivs[i][rI] -= fac;
      
    }
  }
  
  TwoToOne::
  TwoToOne(std::vector<double> &paraValue, 
	   std::vector< std::vector<size_t> > 
	   &indValue ) {
    
    //Do some checks on the parameters and variable indeces
    //
    if( paraValue.size()!=1 ) {
      std::cerr << "TwoToOne::"
		<< "TwoToOne() "
		<< "1 parameters used (see MembraneCycling.h)\n";
      exit(0);
    }
    if( indValue.size() != 2 || indValue[0].size() != 3 || indValue[1].size() != 0 ) {
      std::cerr << "TwoToOne::"
		<< "TwoToOne() "
		<< "Three cell variable indices (first row) and Zero wall variable"
		<< " indices are used." << std::endl;
    exit(0);
    }
    //Set the variable values
    //
    setId("TwoToOne");
    setParameter(paraValue);  
    setVariableIndex(indValue);
    
    //Set the parameter identities
    //
    std::vector<std::string> tmp( numParameter() );
    tmp.resize( numParameter() );
    tmp[0] = "k1";
    setParameterId( tmp );
  }
  
  void TwoToOne::
  derivs(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs ) 
  {  
    size_t numCells = T.numCell();
    size_t r1I = variableIndex(0,0);//reactant1
    size_t r2I = variableIndex(0,1);//reactant2
    size_t pI = variableIndex(0,2);//product
    
    
    
    assert( r1I<cellData[0].size() &&
	    r2I<cellData[0].size() &&
	    pI<cellData[0].size() );
    
    for (size_t i=0; i<numCells; ++i) {
      
      double fac = parameter(0)*cellData[i][r1I]*cellData[i][r2I];
      
      cellDerivs[i][r1I] -= fac;
      cellDerivs[i][r2I] -= fac;
      cellDerivs[i][pI] += fac;
      
    }
  }

  GeneralWall::
  GeneralWall(std::vector<double> &paraValue, 
	  std::vector< std::vector<size_t> > &indValue ) 
  {  
    //
    // Do some checks on the parameters and variable indeces
    //
    if( paraValue.size()!=1 ) {
      std::cerr << "MassAction::GeneralWall::GeneralWall() "
		<< "Uses only one parameter k_f" << std::endl;
      exit(0);
    }
    if( indValue.size() !=2 ) {
      std::cerr << "MassAction::GeneralWall::GeneralWall() "
		<< "Two levels of variable indices are used." << std::endl
		<< "One for reactants and one for products." << std::endl;
      exit(0);
    }
    if( indValue[0].size()<1 ) {
      std::cerr << "MassAction::GeneralWall::GeneralWall() "
		<< "If no reactants given, there will be no reaction..." << std::endl;
      exit(0);
    }
    //
    // Set the variable values
    //
    setId("MassAction::GeneralWall");
    setParameter(paraValue);  
    setVariableIndex(indValue);
    //
    // Set the parameter identities
    //
    std::vector<std::string> tmp( numParameter() );
    tmp.resize( numParameter() );
    tmp[0] = "k_f";
    setParameterId( tmp );
  }
  
  void GeneralWall::
  derivs(Tissue &T,
	 DataMatrix &cellData,
	 DataMatrix &wallData,
	 DataMatrix &vertexData,
	 DataMatrix &cellDerivs,
	 DataMatrix &wallDerivs,
	 DataMatrix &vertexDerivs ) 
  {  
    if( numVariableIndex(0) ) {
      for( size_t wallIndex=0 ; wallIndex<wallData.size() ; wallIndex++ ) {
	double rate = parameter(0);
	for( size_t i=0 ; i< numVariableIndex(0) ; i++ )
	  rate *= wallData[wallIndex][variableIndex(0,i)];
    
	if (rate>0.0) {
	  for( size_t i=0 ; i< numVariableIndex(0) ; ++i )
	    wallDerivs[wallIndex][variableIndex(0,i)] -= rate;
	  
	  for( size_t i=0 ; i< numVariableIndex(1) ; ++i )
	    wallDerivs[wallIndex][variableIndex(1,i)] += rate;
	}
      }
    }
    else { //No reaction defined...
      return;
    }
  }

  void GeneralWall::
  derivsWithAbs(Tissue &T,
		DataMatrix &cellData,
		DataMatrix &wallData,
		DataMatrix &vertexData,
		DataMatrix &cellDerivs,
		DataMatrix &wallDerivs,
		DataMatrix &vertexDerivs,
		DataMatrix &sdydtCell,
		DataMatrix &sdydtWall,
		DataMatrix &sdydtVertex ) 
  {  
    if( numVariableIndex(0) ) {
      for( size_t wallIndex=0 ; wallIndex<wallData.size() ; wallIndex++ ) {
	double rate = parameter(0);
	for( size_t i=0 ; i< numVariableIndex(0) ; i++ )
	  rate *= wallData[wallIndex][variableIndex(0,i)];
	
	if (rate>0.0) {
	  for( size_t i=0 ; i< numVariableIndex(0) ; ++i ) {
	    wallDerivs[wallIndex][variableIndex(0,i)] -= rate;
	    sdydtWall[wallIndex][variableIndex(0,i)] += rate;
	  }
	  
	  for( size_t i=0 ; i< numVariableIndex(1) ; ++i ) {
	    wallDerivs[wallIndex][variableIndex(1,i)] += rate;
	    sdydtWall[wallIndex][variableIndex(1,i)] += rate;
	  }
	}
      }
    }
    else { //No reaction defined...
      return;
    }
  }
  
  OneToTwoWall::
  OneToTwoWall(std::vector<double> &paraValue,
	       std::vector< std::vector<size_t> >
	       &indValue ) {
    
    //Do some checks on the parameters and variable indeces
    //
    if( paraValue.size()!=1 ) {
      std::cerr << "OneToTwo::"
		<< "OneToTwo() "
		<< "1 parameters used (see MembraneCycling.h)\n";
      exit(0);
    }
    if( indValue.size() != 2 || indValue[0].size() != 0 || indValue[1].size() != 3 ) {
      std::cerr << "OneToTwo::"
		<< "OneToTwo() "
		<< "Zero cell variable indices (first row) and three wall variable"
		<< " indices are used." << std::endl;
      exit(0);
    }
    //Set the variable values
    //
    setId("OneToTwoWall");
    setParameter(paraValue);
    setVariableIndex(indValue);
    
    //Set the parameter identities
    //
    std::vector<std::string> tmp( numParameter() );
    tmp.resize( numParameter() );
    tmp[0] = "k1";
    setParameterId( tmp );
  }
  
    void OneToTwoWall::
    derivs(Tissue &T,
           DataMatrix &cellData,
           DataMatrix &wallData,
           DataMatrix &vertexData,
           DataMatrix &cellDerivs,
           DataMatrix &wallDerivs,
           DataMatrix &vertexDerivs )
    {
        size_t numWalls = T.numWall();
        size_t rwI = variableIndex(1,0);//reactant
        size_t pw1I = variableIndex(1,1);//product1
        size_t pw2I = variableIndex(1,2);//product2
             
        assert( rwI<wallData[0].size() &&
               pw1I<wallData[0].size() &&
               pw2I<wallData[0].size() );
        
        for (size_t j=0; j<numWalls; ++j) {
            
            double fac = parameter(0)*wallData[j][rwI];
            double fac2 = parameter(0)*wallData[j][rwI+1];
            
            wallDerivs[j][pw1I] += fac;
            wallDerivs[j][pw2I] += fac;
            wallDerivs[j][rwI] -= fac;
            wallDerivs[j][pw1I+1] += fac2;
            wallDerivs[j][pw2I+1] += fac2;
            wallDerivs[j][rwI+1] -= fac2;
            
        }
    }
   
    TwoToOneWall::
    TwoToOneWall(std::vector<double> &paraValue,
             std::vector< std::vector<size_t> >
             &indValue ) {
        
        //Do some checks on the parameters and variable indeces
        //
        if( paraValue.size()!=1 ) {
            std::cerr << "TwoToOne::"
            << "TwoToOne() "
            << "1 parameters used (see MembraneCycling.h)\n";
            exit(0);
        }
        if( indValue.size() != 2 || indValue[0].size() != 0 || indValue[1].size() != 3 ) {
            std::cerr << "TwoToOne::"
            << "TwoToOne() "
            << "Zero cell variable indices (first row) and Three wall variable"
            << " indices are used." << std::endl;
            exit(0);
        }
        //Set the variable values
        //
        setId("TwoToOneWall");
        setParameter(paraValue);
        setVariableIndex(indValue);
        
        //Set the parameter identities
        //
        std::vector<std::string> tmp( numParameter() );
        tmp.resize( numParameter() );
        tmp[0] = "k1";
        setParameterId( tmp );
    }
    
    void TwoToOneWall::
    derivs(Tissue &T,
           DataMatrix &cellData,
           DataMatrix &wallData,
           DataMatrix &vertexData,
           DataMatrix &cellDerivs,
           DataMatrix &wallDerivs,
           DataMatrix &vertexDerivs ) 
    {  
        size_t numWalls = T.numWall();
        size_t rw1I = variableIndex(1,0);//reactant1
        size_t rw2I = variableIndex(1,1);//reactant2
        size_t pwI = variableIndex(1,2);//product
        
        assert( rw1I<wallData[0].size() &&
               rw2I<wallData[0].size() &&
               pwI<wallData[0].size() );
        
        for (size_t j=0; j<numWalls; ++j) {
            
            double fac = parameter(0)*wallData[j][rw1I]*wallData[j][rw2I];
            double fac2 = parameter(0)*wallData[j][rw1I+1]*wallData[j][rw2I+1];
            
            wallDerivs[j][rw1I] -= fac;
            wallDerivs[j][rw2I] -= fac;
            wallDerivs[j][pwI] += fac;
            wallDerivs[j][rw1I+1] -= fac2;
            wallDerivs[j][rw2I+1] -= fac2;
            wallDerivs[j][pwI+1] += fac2;       
        }
    }
}

