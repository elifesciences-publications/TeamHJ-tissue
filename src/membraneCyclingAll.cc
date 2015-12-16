//
// Filename     : MembraneCyclingAll.cc
// Description  : Classes describing cycling to/from mebrane
// Author(s)    : Laura Brown (laura.brown@slcu.cam.ac.uk)
// Created      : August 2014
// Revision     : $Id:$
//

// This file is a modification of MembraneCycling.cc file by Pau

#include "membraneCyclingAll.h"
#include "baseReaction.h"


namespace MembraneCyclingAll {
  


Constant::
Constant(std::vector<double> &paraValue, 
	      std::vector< std::vector<size_t> > 
	      &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=2 ) {
    std::cerr << "Constant::"
	      << "Constant() "
	      << "2 parameters used (see MembraneCycling.h)\n";
    exit(0);
  }
  if( indValue.size() != 2 || indValue[0].size() != 1 || indValue[1].size() != 1 ) {
    std::cerr << "Constant::"
	      << "Constant() "
	      << "One cell variable indices (first row) and One wall variable"
	      << " indices are used (PIN)." << std::endl;
    exit(0);
  }
  //Set the variable values
  //
  setId("Constant");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "k_on";
  tmp[1] = "k_off";
    setParameterId( tmp );
}

void Constant::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{  
  size_t numCells = T.numCell();
  size_t pI = variableIndex(0,0);//pin
  size_t pwI = variableIndex(1,0);//pin (membrane/wall)


  assert( pI<cellData[0].size() &&
	  pwI<wallData[0].size() );

 for (size_t i=0; i<numCells; ++i) {
	  
       
    //protein cycling
    size_t numWalls = T.cell(i).numWall();
    for (size_t k=0; k<numWalls; ++k) {
      size_t j = T.cell(i).wall(k)->index();
      if( T.cell(i).wall(k)->cell1()->index() == i) {
	//PIN cycling
	double fac = -parameter(0)*cellData[i][pI]+parameter(1)*wallData[j][pwI];

	wallDerivs[j][pwI] -= fac;
	cellDerivs[i][pI] += fac;

      }
      else if( T.cell(i).wall(k)->cell2()->index() == i ) {
	//PIN cycling
	
	double fac = -parameter(0)*cellData[i][pI]+parameter(1)*wallData[j][pwI+1];

	wallDerivs[j][pwI+1] -= fac;
	cellDerivs[i][pI] += fac;
      }
    }
  }
 }





LocalWallFeedbackNonLinear::
LocalWallFeedbackNonLinear(std::vector<double> &paraValue, 
	      std::vector< std::vector<size_t> > 
	      &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=4 ) {
    std::cerr << "LocalWallFeedbackNonLinear::"
	      << "LocalWallFeedbackNonLinear() "
	      << "4 parameters used (see MembraneCycling.h)\n";
    exit(0);
  }
  if( indValue.size() != 2 || indValue[0].size() != 1 || indValue[1].size() != 2 ) {
    std::cerr << "LocalWallFeedbackNonLinear::"
	      << "LocalWallFeedbackNonLinear() "
	      << "One cell variable indices (first row) and One wall variable"
	      << " indices are used (PIN)." << std::endl;
    exit(0);
  }
  //Set the variable values
  //
  setId("LocalWallFeedbackNonLinear");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "k_on";
  tmp[1] = "k_off";
  tmp[2] = "k";
  tmp[3] = "n";
    setParameterId( tmp );
}

void LocalWallFeedbackNonLinear::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{  
  size_t numCells = T.numCell();
  size_t pI = variableIndex(0,0);//pin
  size_t xwI = variableIndex(1,0);//x (membrane/wall)
  size_t pwI = variableIndex(1,1);//pin (membrane/wall)


  assert(  pI<cellData[0].size() &&
	   xwI<wallData[0].size() &&
	   pwI<wallData[0].size() );

 for (size_t i=0; i<numCells; ++i) {
	         
    //protein cycling
    size_t numWalls = T.cell(i).numWall();
    for (size_t k=0; k<numWalls; ++k) {
      size_t j = T.cell(i).wall(k)->index();
     
      if( T.cell(i).wall(k)->cell1()->index() == i ) {
	//PIN cycling
	double fac = parameter(0)*cellData[i][pI]*(std::pow(wallData[j][xwI],parameter(3)))/(std::pow(wallData[j][xwI],parameter(3))+std::pow(parameter(2),parameter(3)))
                    -parameter(1)*wallData[j][pwI]*(std::pow(wallData[j][xwI],parameter(3)))/(std::pow(wallData[j][xwI],parameter(3))+std::pow(parameter(2),parameter(3))) ;
	wallDerivs[j][pwI] += fac;
	cellDerivs[i][pI] -= fac;
      }

      else if( T.cell(i).wall(k)->cell2()->index() == i ) {
  //PIN cycling
	double fac = parameter(0)*cellData[i][pI]*(std::pow(wallData[j][xwI+1],parameter(3)))/(std::pow(wallData[j][xwI+1],parameter(3))+std::pow(parameter(2),parameter(3)))
                    -parameter(1)*wallData[j][pwI+1]*(std::pow(wallData[j][xwI+1],parameter(3)))/(std::pow(wallData[j][xwI+1],parameter(3))+std::pow(parameter(2),parameter(3)));
	wallDerivs[j][pwI+1] += fac;
	cellDerivs[i][pI] -= fac;
      }
    }
  }
 }



LocalWallFeedbackNonLinearInhibition::
LocalWallFeedbackNonLinearInhibition(std::vector<double> &paraValue, 
        std::vector< std::vector<size_t> > 
        &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=4 ) {
    std::cerr << "LocalWallFeedbackNonLinearInhibition::"
        << "LocalWallFeedbackNonLinearInhibition() "
        << "4 parameters used (see MembraneCycling.h)\n";
    exit(0);
  }
  if( indValue.size() != 2 || indValue[0].size() != 1 || indValue[1].size() != 2 ) {
    std::cerr << "LocalWallFeedbackNonLinearInhibition::"
        << "LocalWallFeedbackNonLinearInhibition() "
        << "One cell variable indices (first row) and One wall variable"
        << " indices are used (PIN)." << std::endl;
    exit(0);
  }
  //Set the variable values
  //
  setId("LocalWallFeedbackNonLinearInhibition");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "k_on";
  tmp[1] = "k_off";
  tmp[2] = "k";
  tmp[3] = "n";
    setParameterId( tmp );
}

void LocalWallFeedbackNonLinearInhibition::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{  
  size_t numCells = T.numCell();
  size_t pI = variableIndex(0,0);//pin
  size_t xwI = variableIndex(1,0);//x (membrane/wall)
  size_t pwI = variableIndex(1,1);//pin (membrane/wall)


  assert(  pI<cellData[0].size() &&
     xwI<wallData[0].size() &&
     pwI<wallData[0].size() );

 for (size_t i=0; i<numCells; ++i) {
           
    //protein cycling
    size_t numWalls = T.cell(i).numWall();
    for (size_t k=0; k<numWalls; ++k) {
      size_t j = T.cell(i).wall(k)->index();
     
      if( T.cell(i).wall(k)->cell1()->index() == i ) {
  //PIN cycling
  double fac = parameter(0)*cellData[i][pI]/(std::pow(wallData[j][xwI],parameter(3))+std::pow(parameter(2),parameter(3)))
                    -parameter(1)*wallData[j][pwI]/(1.0+std::pow(wallData[j][xwI],parameter(3))/std::pow(parameter(2),parameter(3))) ;
  wallDerivs[j][pwI] += fac;
  cellDerivs[i][pI] -= fac;
      }

      else if( T.cell(i).wall(k)->cell2()->index() == i ) {
  //PIN cycling
  double fac = parameter(0)*cellData[i][pI]/(std::pow(wallData[j][xwI+1],parameter(3))+std::pow(parameter(2),parameter(3)))
                    -parameter(1)*wallData[j][pwI+1]/(1.0+std::pow(wallData[j][xwI],parameter(3))/std::pow(parameter(2),parameter(3)));
  wallDerivs[j][pwI+1] += fac;
  cellDerivs[i][pI] -= fac;
      }
    }
  }
 }




}



