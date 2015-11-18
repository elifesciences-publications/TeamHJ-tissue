//
// Filename     : MembraneCycling.cc
// Description  : Classes describing cycling to/from mebrane
// Author(s)    : Laura Brown (laura.brown@slcu.cam.ac.uk)
// Created      : August 2014
// Revision     : $Id:$
//
#include "membraneCycling.h"
#include "baseReaction.h"


namespace MembraneCycling {
  


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
      if( T.cell(i).wall(k)->cell1()->index() == i && T.cell(i).wall(k)->cell2() != T.background() ) {
	//PIN cycling
	double fac = -parameter(0)*cellData[i][pI]+parameter(1)*wallData[j][pwI];

	wallDerivs[j][pwI] -= fac;
	cellDerivs[i][pI] += fac;

      }
      else if( T.cell(i).wall(k)->cell2()->index() == i && T.cell(i).wall(k)->cell1() != T.background() ) {
	//PIN cycling
	
	double fac = -parameter(0)*cellData[i][pI]+parameter(1)*wallData[j][pwI+1];

	wallDerivs[j][pwI+1] -= fac;
	cellDerivs[i][pI] += fac;
      }
    }
  }
 }






CrossMembraneNonLinear::
CrossMembraneNonLinear(std::vector<double> &paraValue, 
	      std::vector< std::vector<size_t> > 
	      &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=4 ) {
    std::cerr << "CrossMembraneNonLinear::"
	      << "CrossMembraneNonLinear() "
	      << "4 parameters used (see MembraneCycling.h)\n";
    exit(0);
  }
  if( indValue.size() != 2 || indValue[0].size() != 1 || indValue[1].size() != 1 ) {
    std::cerr << "CrossMembraneNonLinear::"
	      << "CrossMembraneNonLinear() "
	      << "One cell variable indices (first row) and One wall variable"
	      << " indices are used (PIN)." << std::endl;
    exit(0);
  }
  //Set the variable values
  //
  setId("CrossMembraneNonLinear");
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

void CrossMembraneNonLinear::
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
      if( T.cell(i).wall(k)->cell1()->index() == i && T.cell(i).wall(k)->cell2() != T.background() ) {
	//PIN cycling
	double fac =-parameter(0)*cellData[i][pI]*std::pow(wallData[j][pwI+1],parameter(3))/(std::pow(parameter(2),parameter(3))+std::pow(wallData[j][pwI+1],parameter(3)))
                    +parameter(1)*wallData[j][pwI]*std::pow(wallData[j][pwI+1],parameter(3))/(std::pow(parameter(2),parameter(3))+std::pow(wallData[j][pwI+1],parameter(3))) ;
	wallDerivs[j][pwI] -= fac;
	cellDerivs[i][pI] += fac;

      }
      else if( T.cell(i).wall(k)->cell2()->index() == i && T.cell(i).wall(k)->cell1() != T.background() ) {
	//PIN cycling
	double fac =  -parameter(0)*std::pow(wallData[j][pwI],parameter(3))/(std::pow(parameter(2),parameter(3))+std::pow(wallData[j][pwI],parameter(3)))*cellData[i][pI]
                      +parameter(1)*std::pow(wallData[j][pwI],parameter(3))/(std::pow(parameter(2),parameter(3))+std::pow(wallData[j][pwI],parameter(3)))*wallData[j][pwI+1];
	wallDerivs[j][pwI+1] -= fac;
	cellDerivs[i][pI] += fac;
      }
    }
  }
 }







CrossMembraneLinear::
CrossMembraneLinear(std::vector<double> &paraValue, 
	      std::vector< std::vector<size_t> > 
	      &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=2 ) {
    std::cerr << "CrossMembraneLinear::"
	      << "CrossMembraneLinear() "
	      << "2 parameters used (see MembraneCycling.h)\n";
    exit(0);
  }
  if( indValue.size() != 2 || indValue[0].size() != 1 || indValue[1].size() != 1 ) {
    std::cerr << "CrossMembraneLinear::"
	      << "CrossMembraneLinear() "
	      << "One cell variable indices (first row) and One wall variable"
	      << " indices are used (PIN)." << std::endl;
    exit(0);
  }
  //Set the variable values
  //
  setId("CrossMembraneLinear");
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

void CrossMembraneLinear::
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
      if( T.cell(i).wall(k)->cell1()->index() == i && T.cell(i).wall(k)->cell2() != T.background() ) {
	//PIN cycling
	double fac =-parameter(0)*cellData[i][pI]*wallData[j][pwI+1]
                    +parameter(1)*wallData[j][pwI]*wallData[j][pwI+1] ;
	wallDerivs[j][pwI] -= fac;
	cellDerivs[i][pI] += fac;

      }
      else if( T.cell(i).wall(k)->cell2()->index() == i && T.cell(i).wall(k)->cell1() != T.background() ) {
	//PIN cycling
	double fac =  -parameter(0)*wallData[j][pwI]*cellData[i][pI]
                      +parameter(1)*wallData[j][pwI]*wallData[j][pwI+1];
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
     
      if( T.cell(i).wall(k)->cell1()->index() == i && T.cell(i).wall(k)->cell2() != T.background() ) {
	//PIN cycling
	double fac = parameter(0)*cellData[i][pI]*(std::pow(wallData[j][xwI],parameter(3)))/(std::pow(wallData[j][xwI],parameter(3))+std::pow(parameter(2),parameter(3)))
                    -parameter(1)*wallData[j][pwI]*(std::pow(wallData[j][xwI],parameter(3)))/(std::pow(wallData[j][xwI],parameter(3))+std::pow(parameter(2),parameter(3))) ;
	wallDerivs[j][pwI] += fac;
	cellDerivs[i][pI] -= fac;
      }

      else if( T.cell(i).wall(k)->cell2()->index() == i && T.cell(i).wall(k)->cell1() != T.background() ) {
       	//PIN cycling
	double fac = parameter(0)*cellData[i][pI]*(std::pow(wallData[j][xwI+1],parameter(3)))/(std::pow(wallData[j][xwI+1],parameter(3))+std::pow(parameter(2),parameter(3)))
                    -parameter(1)*wallData[j][pwI+1]*(std::pow(wallData[j][xwI+1],parameter(3)))/(std::pow(wallData[j][xwI+1],parameter(3))+std::pow(parameter(2),parameter(3)));
	wallDerivs[j][pwI+1] += fac;
	cellDerivs[i][pI] -= fac;
      }
    }
  }
 }





LocalWallFeedbackLinear::
LocalWallFeedbackLinear(std::vector<double> &paraValue, 
	      std::vector< std::vector<size_t> > 
	      &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=2 ) {
    std::cerr << "LocalWallFeedbackLinear::"
	      << "LocalWallFeedbackLinear() "
	      << "2 parameters used (see MembraneCycling.h)\n";
    exit(0);
  }
  if( indValue.size() != 2 || indValue[0].size() != 1 || indValue[1].size() != 2 ) {
    std::cerr << "LocalWallFeedbackLinear::"
	      << "LocalWallFeedbackLinear() "
	      << "One cell variable indices (first row) and One wall variable"
	      << " indices are used (PIN)." << std::endl;
    exit(0);
  }
  //Set the variable values
  //
  setId("LocalWallFeedbackLinear");
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

void LocalWallFeedbackLinear::
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
  size_t xwI = variableIndex(1,0);//pin (membrane/wall)
  size_t pwI = variableIndex(1,1);//x (membrane/wall)

  assert(  pI<cellData[0].size() &&
	   pwI<wallData[0].size() &&
	   xwI<wallData[0].size() );

 for (size_t i=0; i<numCells; ++i) {
	         
    //protein cycling
    size_t numWalls = T.cell(i).numWall();
    for (size_t k=0; k<numWalls; ++k) {
      size_t j = T.cell(i).wall(k)->index();
     
      if( T.cell(i).wall(k)->cell1()->index() == i && T.cell(i).wall(k)->cell2() != T.background() ) {
	//PIN cycling
	double fac = parameter(0)*cellData[i][pI]*wallData[j][xwI]
                    - parameter(1)*wallData[i][pwI]*wallData[j][xwI] ;
	wallDerivs[j][pwI] += fac;
	cellDerivs[i][pI] -= fac;
      }

      else if( T.cell(i).wall(k)->cell2()->index() == i && T.cell(i).wall(k)->cell1() != T.background() ) {
       	//PIN cycling
	double fac = parameter(0)*wallData[j][xwI+1]*cellData[i][pI]
	  -parameter(1)*wallData[j][xwI+1]*wallData[j][pwI+1];
	wallDerivs[j][pwI+1] += fac;
	cellDerivs[i][pI] -= fac;
      }
    }
  }
 }



CellUpTheGradientNonLinear::
CellUpTheGradientNonLinear(std::vector<double> &paraValue, 
	      std::vector< std::vector<size_t> > 
	      &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=4 ) {
    std::cerr << "CellUpTheGradientNonLinear::"
	      << "CellUpTheGradientNonLinear() "
	      << "4 parameters used (see MembraneCycling.h)\n";
    exit(0);
  }
  if( indValue.size() != 2 || indValue[0].size() != 2 || indValue[1].size() != 1 ) {
    std::cerr << "CellUpTheGradientNonlinear::"
	      << "CellUpTheGradientNonLinear() "
	      << "One cell variable indices (first row) and One wall variable"
	      << " indices are used (PIN)." << std::endl;
    exit(0);
  }
  //Set the variable values
  //
  setId("CellUpTheGradientNonLinear");
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

void CellUpTheGradientNonLinear::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{  
  size_t numCells = T.numCell();
 size_t aI = variableIndex(0,0);//auxin
  size_t pI = variableIndex(0,1);//pin
  size_t pwI = variableIndex(1,0);//pin (membrane/wall)


  assert( aI<cellData[0].size() &&
	  pI<cellData[0].size() &&
	  pwI<wallData[0].size() );

 for (size_t i=0; i<numCells; ++i) {
	  
       
    //protein cycling
    size_t numWalls = T.cell(i).numWall();
    

     for (size_t k=0; k<numWalls; ++k) {
      size_t j = T.cell(i).wall(k)->index();

      if( T.cell(i).wall(k)->cell1()->index() == i && T.cell(i).wall(k)->cell2() != T.background() ) {
	size_t cellNeigh = T.cell(i).wall(k)->cell2()->index();
	//PIN cycling
	double fac = parameter(0)*cellData[i][pI]*(std::pow(cellData[cellNeigh][aI],parameter(3)))/(std::pow(cellData[cellNeigh][aI],parameter(3))+std::pow(parameter(2),parameter(3)))
                    -parameter(1)*wallData[j][pwI]*(std::pow(cellData[cellNeigh][aI],parameter(3)))/(std::pow(cellData[cellNeigh][aI],parameter(3))+std::pow(parameter(2),parameter(3))) ;
	wallDerivs[j][pwI] += fac;
	cellDerivs[i][pI] -= fac;
      }

      else if( T.cell(i).wall(k)->cell2()->index() == i && T.cell(i).wall(k)->cell1() != T.background() ) {
	size_t cellNeigh = T.cell(i).wall(k)->cell1()->index();

	//PIN cycling
	double fac = parameter(0)*(std::pow(cellData[cellNeigh][aI],parameter(3)))/(std::pow(cellData[cellNeigh][aI],parameter(3))+std::pow(parameter(2),parameter(3)))*cellData[i][pI]
                    -parameter(1)*(std::pow(cellData[cellNeigh][aI],parameter(3)))/(std::pow(cellData[cellNeigh][aI],parameter(3))+std::pow(parameter(2),parameter(3)))*wallData[j][pwI+1];
	wallDerivs[j][pwI+1] += fac;
	cellDerivs[i][pI] -= fac;
      }
    }
  }
 }





CellUpTheGradientLinear::
CellUpTheGradientLinear(std::vector<double> &paraValue, 
	      std::vector< std::vector<size_t> > 
	      &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=2 ) {
    std::cerr << "CellUpTheGradientLinear::"
	      << "CellUpTheGradientLinear() "
	      << "4 parameters used (see MembraneCycling.h)\n";
    exit(0);
  }
  if( indValue.size() != 2 || indValue[0].size() != 2 || indValue[1].size() != 1 ) {
    std::cerr << "CellUpTheGradientLinear::"
	      << "CellUpTheGradientLinear() "
	      << "One cell variable indices (first row) and One wall variable"
	      << " indices are used (PIN)." << std::endl;
    exit(0);
  }
  //Set the variable values
  //
  setId("CellUpTheGradientLinear");
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

void CellUpTheGradientLinear::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{  
  size_t numCells = T.numCell();
 size_t aI = variableIndex(0,0);//auxin
  size_t pI = variableIndex(0,1);//pin
  size_t pwI = variableIndex(1,0);//pin (membrane/wall)


  assert( aI<cellData[0].size() &&
	  pI<cellData[0].size() &&
	  pwI<wallData[0].size() );

 for (size_t i=0; i<numCells; ++i) {
	  
       
    //protein cycling
    size_t numWalls = T.cell(i).numWall();
    

     for (size_t k=0; k<numWalls; ++k) {
      size_t j = T.cell(i).wall(k)->index();

      if( T.cell(i).wall(k)->cell1()->index() == i && T.cell(i).wall(k)->cell2() != T.background() ) {
	size_t cellNeigh = T.cell(i).wall(k)->cell2()->index();
	//PIN cycling
	double fac = parameter(0)*cellData[i][pI]*cellData[cellNeigh][aI]
                    -parameter(1)*wallData[j][pwI]*cellData[cellNeigh][aI] ;
	wallDerivs[j][pwI] += fac;
	cellDerivs[i][pI] -= fac;
      }

      else if( T.cell(i).wall(k)->cell2()->index() == i && T.cell(i).wall(k)->cell1() != T.background() ) {
	size_t cellNeigh = T.cell(i).wall(k)->cell1()->index();

	//PIN cycling
	double fac = parameter(0)*cellData[cellNeigh][aI]*cellData[i][pI]
                    -parameter(1)*cellData[cellNeigh][aI]*wallData[j][pwI+1];
	wallDerivs[j][pwI+1] += fac;
	cellDerivs[i][pI] -= fac;
      }
    }
  }
 }




CellFluxExocytosis::
CellFluxExocytosis(std::vector<double> &paraValue, 
	      std::vector< std::vector<size_t> > 
	      &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=3 ) {
    std::cerr << "CellFluxExocytosis::"
	      << "CellFluxExocytosis() "
	      << "3 parameters used (see network.h)\n";
    exit(0);
  }
  if( indValue.size() != 2 || indValue[0].size() != 2 || indValue[1].size() != 1 ) {
    std::cerr << "CellFluxExocytosis::"
	      << "CellFluxExocytosis() "
	      << "Two cell variable indices (first row) and two wall variable"
	      << " indices are used (auxin,PIN)." << std::endl;
    exit(0);
  }
  //Set the variable values
  //
  setId("CellFluxExocytosis");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "D";
  tmp[1] = "T";
  tmp[2] = "exo_PIN";
  setParameterId( tmp );
}

void CellFluxExocytosis::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{  
  size_t numCells = T.numCell();
  size_t aI = variableIndex(0,0);//auxin
  size_t pI = variableIndex(0,1);//pin
  size_t pwI = variableIndex(1,0);//pin (membrane/wall)


  assert( aI<cellData[0].size() &&
	  pI<cellData[0].size() &&
	  pwI<wallData[0].size() );

 for (size_t i=0; i<numCells; ++i) {
	  
   

    
    //Auxin transport and protein cycling
    size_t numWalls = T.cell(i).numWall();
    for (size_t k=0; k<numWalls; ++k) {
      size_t j = T.cell(i).wall(k)->index();
      if( T.cell(i).wall(k)->cell1()->index() == i && T.cell(i).wall(k)->cell2() != T.background() ) {
	// cell-cell transport
	size_t iNeighbor = T.cell(i).wall(k)->cell2()->index();
	double fac = parameter(0)*cellData[i][aI] + parameter(1)*cellData[i][aI]*wallData[j][pwI]; //p_2 A_i + p_4 A_i P_ij
	
	cellDerivs[i][aI] -= fac;
	cellDerivs[iNeighbor][aI] += fac;
	if (parameter(0)*cellData[i][aI]+parameter(1)*cellData[i][aI]*wallData[j][pwI]-parameter(0)*cellData[iNeighbor][aI]-parameter(1)*cellData[iNeighbor][aI]*wallData[j][pwI+1]>=0) {
	//PIN cycling
	double fac =  - parameter(2)*cellData[i][pI]* std::pow(parameter(0)*cellData[i][aI]+parameter(1)*cellData[i][aI]*wallData[j][pwI]-parameter(0)*cellData[iNeighbor][aI]-parameter(1)*cellData[iNeighbor][aI]*wallData[j][pwI+1],2);
	wallDerivs[j][pwI] -= fac;
	cellDerivs[i][pI] += fac;}

	else if (parameter(0)*cellData[i][aI]+parameter(1)*cellData[i][aI]*wallData[j][pwI]-parameter(0)*cellData[iNeighbor][aI]-parameter(1)*cellData[iNeighbor][aI]*wallData[j][pwI+1]<0){
	//PIN cycling
        double	fac =0;
	wallDerivs[j][pwI] -= 0;
	cellDerivs[i][pI] += 0;}




      }
      else if( T.cell(i).wall(k)->cell2()->index() == i && T.cell(i).wall(k)->cell1() != T.background() ) {
	// cell-cell transport
	size_t iNeighbor = T.cell(i).wall(k)->cell1()->index();
	// cell-cell transport
	double fac = parameter(0)*cellData[i][aI] + parameter(1)*cellData[i][aI]*wallData[j][pwI+1]; //p_2 A_i + p_4 A_i P_ij
	
	cellDerivs[i][aI] -= fac;
	cellDerivs[iNeighbor][aI] += fac;

      
	//PIN cycling
	if (parameter(0)*cellData[i][aI]+parameter(1)*cellData[i][aI]*wallData[j][pwI+1]-parameter(0)*cellData[iNeighbor][aI]-parameter(1)*cellData[iNeighbor][aI]*wallData[j][pwI]>=0) {

	  double fac = -parameter(2)*cellData[i][pI]*std::pow(parameter(0)*cellData[i][aI]+parameter(1)*cellData[i][aI]*wallData[j][pwI+1]-parameter(0)*cellData[iNeighbor][aI]-parameter(1)*cellData[iNeighbor][aI]*wallData[j][pwI],2);
	wallDerivs[j][pwI+1] -= fac;
	cellDerivs[i][pI] += fac;}
	else if (parameter(0)*cellData[i][aI]+parameter(1)*cellData[i][aI]*wallData[j][pwI+1]-parameter(0)*cellData[iNeighbor][aI]-parameter(1)*cellData[iNeighbor][aI]*wallData[j][pwI]<0) {

	double fac = 0;
	wallDerivs[j][pwI+1] -= 0;
	cellDerivs[i][pI] += 0;}
      
      }
    }
  }
}
























InternalCellNonLinear::
InternalCellNonLinear(std::vector<double> &paraValue, 
	      std::vector< std::vector<size_t> > 
	      &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=4 ) {
    std::cerr << "InternalCellNonLinear::"
	      << "InternalCellNonLinear() "
	      << "4 parameters used (see MembraneCycling.h)\n";
    exit(0);
  }
  if( indValue.size() != 2 || indValue[0].size() != 2 || indValue[1].size() != 1 ) {
    std::cerr << "InternalCellNonlinear::"
	      << "InternalCellNonLinear() "
	      << "One cell variable indices (first row) and One wall variable"
	      << " indices are used (PIN)." << std::endl;
    exit(0);
  }
  //Set the variable values
  //
  setId("InternalCellNonLinear");
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

void InternalCellNonLinear::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{  
  size_t numCells = T.numCell();
 size_t aI = variableIndex(0,0);//auxin
  size_t pI = variableIndex(0,1);//pin
  size_t pwI = variableIndex(1,0);//pin (membrane/wall)


  assert( aI<cellData[0].size() &&
	  pI<cellData[0].size() &&
	  pwI<wallData[0].size() );

 for (size_t i=0; i<numCells; ++i) {
	  
       
    //protein cycling
    size_t numWalls = T.cell(i).numWall();
    

     for (size_t k=0; k<numWalls; ++k) {
      size_t j = T.cell(i).wall(k)->index();

      if( T.cell(i).wall(k)->cell1()->index() == i && T.cell(i).wall(k)->cell2() != T.background() ) {

	//PIN cycling
	double fac = parameter(0)*cellData[i][pI]*(std::pow(cellData[i][aI],parameter(3)))/(std::pow(cellData[i][aI],parameter(3))+std::pow(parameter(2),parameter(3)))
                    -parameter(1)*wallData[j][pwI]*(std::pow(cellData[i][aI],parameter(3)))/(std::pow(cellData[i][aI],parameter(3))+std::pow(parameter(2),parameter(3))) ;
	wallDerivs[j][pwI] += fac;
	cellDerivs[i][pI] -= fac;
      }

      else if( T.cell(i).wall(k)->cell2()->index() == i && T.cell(i).wall(k)->cell1() != T.background() ) {


	//PIN cycling
	double fac = parameter(0)*(std::pow(cellData[i][aI],parameter(3)))/(std::pow(cellData[i][aI],parameter(3))+std::pow(parameter(2),parameter(3)))*cellData[i][pI]
                    -parameter(1)*(std::pow(cellData[i][aI],parameter(3)))/(std::pow(cellData[i][aI],parameter(3))+std::pow(parameter(2),parameter(3)))*wallData[j][pwI+1];
	wallDerivs[j][pwI+1] += fac;
	cellDerivs[i][pI] -= fac;
      }
    }
  }
 }





InternalCellLinear::
InternalCellLinear(std::vector<double> &paraValue, 
	      std::vector< std::vector<size_t> > 
	      &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=2 ) {
    std::cerr << "InternalCellLinear::"
	      << "InternalCellLinear() "
	      << "4 parameters used (see MembraneCycling.h)\n";
    exit(0);
  }
  if( indValue.size() != 2 || indValue[0].size() != 2 || indValue[1].size() != 1 ) {
    std::cerr << "InternalCellLinear::"
	      << "InternalCellLinear() "
	      << "One cell variable indices (first row) and One wall variable"
	      << " indices are used (PIN)." << std::endl;
    exit(0);
  }
  //Set the variable values
  //
  setId("InternalCellLinear");
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

void InternalCellLinear::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{  
  size_t numCells = T.numCell();
 size_t aI = variableIndex(0,0);//auxin
  size_t pI = variableIndex(0,1);//pin
  size_t pwI = variableIndex(1,0);//pin (membrane/wall)


  assert( aI<cellData[0].size() &&
	  pI<cellData[0].size() &&
	  pwI<wallData[0].size() );

 for (size_t i=0; i<numCells; ++i) {
	  
       
    //protein cycling
    size_t numWalls = T.cell(i).numWall();
    

     for (size_t k=0; k<numWalls; ++k) {
      size_t j = T.cell(i).wall(k)->index();

      if( T.cell(i).wall(k)->cell1()->index() == i && T.cell(i).wall(k)->cell2() != T.background() ) {

	//PIN cycling
	double fac = parameter(0)*cellData[i][pI]*cellData[i][aI]
                    -parameter(1)*wallData[j][pwI]*cellData[i][aI] ;
	wallDerivs[j][pwI] += fac;
	cellDerivs[i][pI] -= fac;
      }

      else if( T.cell(i).wall(k)->cell2()->index() == i && T.cell(i).wall(k)->cell1() != T.background() ) {
		//PIN cycling
	double fac = parameter(0)*cellData[i][aI]*cellData[i][pI]
                    -parameter(1)*cellData[i][aI]*wallData[j][pwI+1];
	wallDerivs[j][pwI+1] += fac;
	cellDerivs[i][pI] -= fac;
      }
    }
  }
 }





PINFeedbackNonLinear::
PINFeedbackNonLinear(std::vector<double> &paraValue, 
	      std::vector< std::vector<size_t> > 
	      &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=4 ) {
    std::cerr << "PINFeedbackNonLinear::"
	      << "PINFeedbackNonLinear() "
	      << "4 parameters used (see MembraneCycling.h)\n";
    exit(0);
  }
  if( indValue.size() != 2 || indValue[0].size() != 1 || indValue[1].size() != 1 ) {
    std::cerr << "PINFeedbackNonLinear::"
	      << "PINFeedbackNonLinear() "
	      << "One cell variable indices (first row) and One wall variable"
	      << " indices are used (PIN)." << std::endl;
    exit(0);
  }
  //Set the variable values
  //
  setId("PINFeedbackNonLinear");
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

void PINFeedbackNonLinear::
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

  assert(  pI<cellData[0].size() &&
	   pwI<wallData[0].size());

 for (size_t i=0; i<numCells; ++i) {
	         
    //protein cycling
    size_t numWalls = T.cell(i).numWall();
    for (size_t k=0; k<numWalls; ++k) {
      size_t j = T.cell(i).wall(k)->index();
     
      if( T.cell(i).wall(k)->cell1()->index() == i && T.cell(i).wall(k)->cell2() != T.background() ) {
	//PIN cycling
	double fac = parameter(0)*cellData[i][pI]*(std::pow(wallData[j][pwI],parameter(2)))/(std::pow(wallData[j][pwI],parameter(2))+std::pow(parameter(2),parameter(3)))
                    - parameter(1)*wallData[j][pwI]*(std::pow(wallData[j][pwI],parameter(2)))/(std::pow(wallData[j][pwI],parameter(2))+std::pow(parameter(2),parameter(3))) ;
	wallDerivs[j][pwI] += fac;
	cellDerivs[i][pI] -= fac;
      }

      else if( T.cell(i).wall(k)->cell2()->index() == i && T.cell(i).wall(k)->cell1() != T.background() ) {
       	//PIN cycling
	double fac = parameter(0)*(std::pow(wallData[j][pwI+1],parameter(3)))/(std::pow(wallData[j][pwI+1],parameter(3))+std::pow(parameter(2),parameter(3)))*cellData[i][pI]
	  -parameter(1)*(std::pow(wallData[j][pwI+1],parameter(3)))/(std::pow(wallData[j][pwI+1],parameter(3))+std::pow(parameter(2),parameter(3)))*wallData[j][pwI+1];
	wallDerivs[j][pwI+1] += fac;
	cellDerivs[i][pI] -= fac;
      }
    }
  }
 }





PINFeedbackLinear::
PINFeedbackLinear(std::vector<double> &paraValue, 
	      std::vector< std::vector<size_t> > 
	      &indValue ) {
  
  //Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=2 ) {
    std::cerr << "PINFeedbackLinear::"
	      << "PINFeedbackLinear() "
	      << "2 parameters used (see MembraneCycling.h)\n";
    exit(0);
  }
  if( indValue.size() != 2 || indValue[0].size() != 1 || indValue[1].size() != 1 ) {
    std::cerr << "PINFeedbackLinear::"
	      << "PINFeedbackLinear() "
	      << "One cell variable indices (first row) and One wall variable"
	      << " indices are used (PIN)." << std::endl;
    exit(0);
  }
  //Set the variable values
  //
  setId("PINFeedbackLinear");
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

void PINFeedbackLinear::
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
 

  assert(  pI<cellData[0].size() &&
	   pwI<wallData[0].size() );

 for (size_t i=0; i<numCells; ++i) {
	         
    //protein cycling
    size_t numWalls = T.cell(i).numWall();
    for (size_t k=0; k<numWalls; ++k) {
      size_t j = T.cell(i).wall(k)->index();
     
      if( T.cell(i).wall(k)->cell1()->index() == i && T.cell(i).wall(k)->cell2() != T.background() ) {
	//PIN cycling
	double fac = parameter(0)*cellData[i][pI]*wallData[j][pwI]
                    - parameter(1)*wallData[i][pwI]*wallData[j][pwI] ;
	wallDerivs[j][pwI] += fac;
	cellDerivs[i][pI] -= fac;
      }

      else if( T.cell(i).wall(k)->cell2()->index() == i && T.cell(i).wall(k)->cell1() != T.background() ) {
       	//PIN cycling
	double fac = parameter(0)*wallData[j][pwI+1]*cellData[i][pI]
	  -parameter(1)*wallData[j][pwI+1]*wallData[j][pwI+1];
	wallDerivs[j][pwI+1] += fac;
	cellDerivs[i][pI] -= fac;
      }
    }
  }
 }








}



