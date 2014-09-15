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
 







}

