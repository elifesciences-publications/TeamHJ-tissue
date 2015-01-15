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

