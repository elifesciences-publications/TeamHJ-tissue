//
// Filename     : grn.cc
// Description  : Classes describing gene regulatory updates. Converted from organism.
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : February 2013
// Revision     : $Id: grn.cc 556 2012-11-15 15:30:35Z henrik $
//

#include"tissue.h"
#include "baseReaction.h"
#include "grn.h"
#include<cstdlib>

Hill::
Hill(std::vector<double> &paraValue, 
     std::vector< std::vector<size_t> > 
     &indValue ) 
{
  // Do some checks on the parameters and variable indeces
  //
  if( indValue.size() != 3 || indValue[0].size() != 1 ) {
    std::cerr << "Hill::Hill() "
	      << "Three levels of variable indices are used, "
	      << "the first for the updated molecule, 2nd for activators and 3rd for repressors"
	      << std::endl;
    exit(EXIT_FAILURE);
  }
  if( 1 + 2*(indValue[1].size()+indValue[2].size()) != paraValue.size() ) {
    std::cerr << "Hill::Hill() "
	      << "Number of parameters does not agree with number of "
	      << "activators/repressors.\n"
	      << indValue[1].size() << " activators, " 
	      << indValue[2].size() << " repressors, and "
	      << paraValue.size() << " parameters\n"
 	      << "One plus pairs of parameters (V_max,K_half1,n_Hill1,"
 	      << "K2,n2,...) must be given.\n";
    exit(EXIT_FAILURE);
  }
  
  // Set the variable values
  //
  setId("Hill");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  // Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "V_max";
  size_t count=1;
  for (size_t i=0; i<indValue[1].size(); ++i) {
    tmp[count++] = "K_A";
    tmp[count++] = "n_A";
  }
  for (size_t i=0; i<indValue[2].size(); ++i) {
    tmp[count++] = "K_R";
    tmp[count++] = "n_R";
  }
  setParameterId( tmp );
}

void Hill::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs )
{
  // Do the update for each cell
  size_t numCells = T.numCell();
  size_t cIndex = variableIndex(0,0);
  //For each cell
  for (size_t cellI = 0; cellI < numCells; ++cellI) {      
  //
    double contribution=parameter(0);
    size_t parameterIndex=1;
    // Activator contributions
    for( size_t i=0 ; i<numVariableIndex(1) ; i++ ) {
      double c = std::pow(cellData[cellI][variableIndex(1,i)], 
      	   parameter(parameterIndex+1));
      contribution *= c
      / ( std::pow(parameter(parameterIndex),parameter(parameterIndex+1)) + c );
      parameterIndex+=2;
    }
    // Repressor contributions
   for( size_t i=0 ; i<numVariableIndex(2) ; i++ ) {
     double c = std::pow(parameter(parameterIndex),parameter(parameterIndex+1));
     contribution *= c /
       ( c + std::pow(cellData[cellI][variableIndex(2,i)],
		      parameter(parameterIndex+1)) );
     parameterIndex+=2;
   }   
    cellDerivs[cellI][cIndex] += contribution; 
  }
}


void Hill::
derivsWithAbs(Tissue &T,
        DataMatrix &cellData,
        DataMatrix &wallData,
        DataMatrix &vertexData,
        DataMatrix &cellDerivs,
        DataMatrix &wallDerivs,
        DataMatrix &vertexDerivs,
        DataMatrix &sdydtCell,
        DataMatrix &sdydtWall,
        DataMatrix &sdydtVertex ){
  // Do the update for each cell
  size_t numCells = T.numCell();
  size_t cIndex = variableIndex(0,0);
  //For each cell
  for (size_t cellI = 0; cellI < numCells; ++cellI) {      
  //
    double contribution=parameter(0);
    size_t parameterIndex=1;
    // Activator contributions
    for( size_t i=0 ; i<numVariableIndex(1) ; i++ ) {
      double c = std::pow(cellData[cellI][variableIndex(1,i)], 
           parameter(parameterIndex+1));
      contribution *= c
      / ( std::pow(parameter(parameterIndex),parameter(parameterIndex+1)) + c );
      parameterIndex+=2;
    }
    // Repressor contributions
   for( size_t i=0 ; i<numVariableIndex(2) ; i++ ) {
     double c = std::pow(parameter(parameterIndex),parameter(parameterIndex+1));
     contribution *= c /
       ( c + std::pow(cellData[cellI][variableIndex(2,i)],
          parameter(parameterIndex+1)) );
     parameterIndex+=2;
   }   
    cellDerivs[cellI][cIndex] += contribution; 
    sdydtCell[cellI][cIndex] += contribution;    

  }
}


HillGeneralOne::HillGeneralOne(std::vector<double> &paraValue, 
			       std::vector< std::vector<size_t> > 
			       &indValue ) 
{
  // Do some checks on the parameters and variable indeces
  //
  if (indValue.size() != 2 || indValue[0].size() != 1 || indValue[1].size() != 1) {
    std::cerr << "HillGeneralOne::HillGeneralOne() "
	      << "The variable to be updated is used from first level and one variable"
	      << " index (TF) is used from second." << std::endl;
    exit(EXIT_FAILURE);
  }
  if (paraValue.size() != 4) {
    std::cerr << "HillGeneralOne::HillGeneralOne() "
	      << "Uses four parameters (Vunbound, Vbound, K_H, n_H)."
	      << std::endl;
    exit(EXIT_FAILURE);
  }
  
  // Set the variable values
  //
  setId("HillGeneralOne");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  // Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "V_unbound";
  tmp[1] = "V_bound";
  tmp[2] = "K_half";
  tmp[3] = "n_Hill";
  setParameterId( tmp );
}

void HillGeneralOne::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs )
{
  double KPow = std::pow(parameter(2),parameter(3));   
  size_t cIndex = variableIndex(0,0);
  
  // Do the update for each cell
  size_t numCells = T.numCell();
  //For each cell
  for (size_t cellI = 0; cellI < numCells; ++cellI) {      
    double tfPow = std::pow(cellData[cellI][variableIndex(1,0)],parameter(3)); 

    cellDerivs[cellI][cIndex] += (parameter(0)*KPow + parameter(1)*tfPow) /
      (KPow+tfPow);
  }
}


void HillGeneralOne::
derivsWithAbs(Tissue &T,
        DataMatrix &cellData,
        DataMatrix &wallData,
        DataMatrix &vertexData,
        DataMatrix &cellDerivs,
        DataMatrix &wallDerivs,
        DataMatrix &vertexDerivs,
        DataMatrix &sdydtCell,
        DataMatrix &sdydtWall,
        DataMatrix &sdydtVertex) 
{
  double KPow = std::pow(parameter(2),parameter(3));   
  size_t cIndex = variableIndex(0,0);
  
  // Do the update for each cell
  size_t numCells = T.numCell();
  //For each cell
  for (size_t cellI = 0; cellI < numCells; ++cellI) {      
    double tfPow = std::pow(cellData[cellI][variableIndex(1,0)],parameter(3)); 

    cellDerivs[cellI][cIndex] += (parameter(0)*KPow + parameter(1)*tfPow) /
      (KPow+tfPow);
    sdydtCell[cellI][cIndex]+=(parameter(0)*KPow + parameter(1)*tfPow) /
      (KPow+tfPow);
  }
}



HillGeneralOne_TwoInputs::HillGeneralOne_TwoInputs(std::vector<double> &paraValue, 
             std::vector< std::vector<size_t> > 
             &indValue ) 
{
  // Do some checks on the parameters and variable indeces
  //
  if (indValue.size() != 2 || indValue[0].size() != 1 || indValue[1].size() != 2) {
    std::cerr << "HillGeneralOne_TwoInputs::HillGeneralOne_TwoInputs() "
        << "The variable to be updated is used from first level and two variables"
        << " index (two inputs acting as a single TF) are used from second." << std::endl;
    exit(EXIT_FAILURE);
  }
  if (paraValue.size() != 4) {
    std::cerr << "HillGeneralOne_TwoInputs::HillGeneralOne_TwoInputs() "
        << "Uses four parameters (Vunbound, Vbound, K_H, n_H)."
        << std::endl;
    exit(EXIT_FAILURE);
  }
  
  // Set the variable values
  //
  setId("HillGeneralOne_TwoInputs");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  // Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "V_unbound";
  tmp[1] = "V_bound";
  tmp[2] = "K_half";
  tmp[3] = "n_Hill";
  setParameterId( tmp );
}

void HillGeneralOne_TwoInputs::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs )
{
  double KPow = std::pow(parameter(2),parameter(3));   
  size_t cIndex = variableIndex(0,0);
  
  // Do the update for each cell
  size_t numCells = T.numCell();
  //For each cell
  for (size_t cellI = 0; cellI < numCells; ++cellI) {      
    double tfPow = std::pow(cellData[cellI][variableIndex(1,0)]+cellData[cellI][variableIndex(1,1)],parameter(3)); 

    cellDerivs[cellI][cIndex] += (parameter(0)*KPow + parameter(1)*tfPow) /
      (KPow+tfPow);
  }
}


void HillGeneralOne_TwoInputs::
derivsWithAbs(Tissue &T,
        DataMatrix &cellData,
        DataMatrix &wallData,
        DataMatrix &vertexData,
        DataMatrix &cellDerivs,
        DataMatrix &wallDerivs,
        DataMatrix &vertexDerivs,
        DataMatrix &sdydtCell,
        DataMatrix &sdydtWall,
        DataMatrix &sdydtVertex) 
{
  double KPow = std::pow(parameter(2),parameter(3));   
  size_t cIndex = variableIndex(0,0);
  
  // Do the update for each cell
  size_t numCells = T.numCell();
  //For each cell
  for (size_t cellI = 0; cellI < numCells; ++cellI) {      
    double tfPow = std::pow((cellData[cellI][variableIndex(1,0)]+cellData[cellI][variableIndex(1,1)]),parameter(3)); 

    cellDerivs[cellI][cIndex] += (parameter(0)*KPow + parameter(1)*tfPow) /
      (KPow+tfPow);
    sdydtCell[cellI][cIndex]+=(parameter(0)*KPow + parameter(1)*tfPow) /
      (KPow+tfPow);
  }
}



HillGeneralTwo::
HillGeneralTwo(std::vector<double> &paraValue, 
	       std::vector< std::vector<size_t> > 
	       &indValue ) 
{
  // Do some checks on the parameters and variable indeces
  //
  if (indValue.size() != 2 || indValue[0].size() != 1 || indValue[1].size() != 2) {
    std::cerr << "HillGeneralTwo::HillGeneralTwo() "
	      << "Two levels of indices are used, first with variable to be updated and "
	      << "and second with two variable indices (TFs)." << std::endl;
    exit(EXIT_FAILURE);
  }
  if (paraValue.size() != 8) {
    std::cerr << "HillGeneralTwo::HillGeneralTwo() "
	      << "Uses eight parameters (Vunbound, Vfirstbound, Vsecondbound, Vbothbound, "
	      << "K1_H, n1_H, K2_H, n2_H)." << std::endl;
    exit(EXIT_FAILURE);
  }
  
  // Set the variable values
  //
  setId("HillGeneralTwo");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  // Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "V_unbound";
  tmp[1] = "V_firstbound";
  tmp[2] = "V_secondbound";
  tmp[3] = "V_bothbound";
  tmp[4] = "K1_half";
  tmp[5] = "n1_Hill";
  tmp[6] = "K2_half";
  tmp[7] = "n2_Hill";
  setParameterId( tmp );
}

void HillGeneralTwo::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs ) 
{
  double K1Pow = std::pow(parameter(4),parameter(5));   
  double K2Pow = std::pow(parameter(6),parameter(7));   
  
  size_t cIndex = variableIndex(0,0);
  
  // Do the update for each cell
  size_t numCells = T.numCell();
  //For each cell
  for (size_t cellI = 0; cellI < numCells; ++cellI) {      
    double tf1Pow = std::pow(cellData[cellI][variableIndex(1,0)],parameter(5)); 
    double tf2Pow = std::pow(cellData[cellI][variableIndex(1,1)],parameter(7)); 
    
    
    cellDerivs[cellI][cIndex] += (parameter(0)*K1Pow*K2Pow + parameter(1)*tf1Pow*K2Pow + 
				  parameter(2)*K1Pow*tf2Pow + parameter(3)*tf1Pow*tf2Pow)
      / ( (K1Pow+tf1Pow)*(K2Pow+tf2Pow) );
  }
}

HillGeneralThree::
HillGeneralThree(std::vector<double> &paraValue, 
		 std::vector< std::vector<size_t> > 
		 &indValue ) 
{
  // Do some checks on the parameters and variable indeces
  //
  if (indValue.size() != 2 || indValue[0].size() != 1 || indValue[1].size() != 3) {
    std::cerr << "HillGeneralThree::HillGeneralThree() "
	      << "Two levels of indices are used, first with variable to be updated and "
	      << "and second with three variable indices (TFs)." << std::endl;
    exit(EXIT_FAILURE);
  }
  if (paraValue.size() != 14) {
    std::cerr << "HillGeneralThree::HillGeneralThree() "
	      << "Uses fourteen parameters (V000, V100, V010, V001, V110, V101, V011, V111, "
	      << "K1_H, n1_H, K2_H, n2_H, K3_H, n3_H)." << std::endl;
    exit(EXIT_FAILURE);
  }
  
  // Set the variable values
  //
  setId("HillGeneralThree");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  // Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "V_000";
  tmp[1] = "V_100";
  tmp[2] = "V_010";
  tmp[3] = "V_001";
  tmp[4] = "V_110";
  tmp[5] = "V_101";
  tmp[6] = "V_011";
  tmp[7] = "V_111";
  tmp[8] = "K1_half";
  tmp[9] = "n1_Hill";
  tmp[10] = "K2_half";
  tmp[11] = "n2_Hill";
  tmp[12] = "K3_half";
  tmp[13] = "n3_Hill";
  setParameterId( tmp );
}

void HillGeneralThree::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs )
{
  double K1Pow = std::pow(parameter(8),parameter(9));   
  double K2Pow = std::pow(parameter(10),parameter(11));   
  double K3Pow = std::pow(parameter(12),parameter(13));   

  size_t cIndex = variableIndex(0,0);
  
  // Do the update for each cell
  size_t numCells = T.numCell();
  //For each cell
  for (size_t cellI = 0; cellI < numCells; ++cellI) {      
    double tf1Pow = std::pow(cellData[cellI][variableIndex(1,0)],parameter(9)); 
    double tf2Pow = std::pow(cellData[cellI][variableIndex(1,1)],parameter(11)); 
    double tf3Pow = std::pow(cellData[cellI][variableIndex(1,2)],parameter(13)); 

    cellDerivs[cellI][cIndex] += (parameter(0)*K1Pow*K2Pow*K3Pow + parameter(1)*tf1Pow*K2Pow*K3Pow
				  + parameter(2)*K1Pow*tf2Pow*K3Pow + parameter(3)*K1Pow*K2Pow*tf3Pow
				  + parameter(4)*tf1Pow*tf2Pow*K3Pow + parameter(5)*tf1Pow*K2Pow*tf3Pow
				  + parameter(6)*K1Pow*tf2Pow*tf3Pow + parameter(7)*tf1Pow*tf2Pow*tf3Pow )
      / ( (K1Pow+tf1Pow)*(K2Pow+tf2Pow)*(K3Pow+tf3Pow) );
  }
}

Grn::Grn(std::vector<double> &paraValue, 
	 std::vector< std::vector<size_t> > 
	 &indValue ) 
{
  // Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()<2) {
    std::cerr << "Grn::Grn() At least two parameters (h,tau) must be given.\n";
    exit(0);
  }
  if( indValue.size() != 2 || indValue[0].size() !=1 ) {
    std::cerr << "Grn::Grn() "
	      << "Two levels of variable index can be used." 
	      << " The first is the cell variable to be updated, and the second"
	      << " is the list of input variables." << std::endl;
    exit(0);
  }
  if( ( paraValue.size() != indValue[1].size()+2 ) ) {
    std::cerr << "Grn::Grn() "
	      << "Wrong number of parameters/variable indeces given (see documentation)."
	      << std::endl;
    std::cerr << paraValue.size() << " " << indValue.size();
    std::cerr << " " << indValue[1].size() << std::endl;
    exit(0);
  }
  
  // Set the variable values
  //
  setId("Grn");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  // Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "h";
  tmp[1] = "tau";
  for(size_t i=2; i<numParameter(); i++ )
    tmp[i] = "T_i";
  setParameterId( tmp );
}

void Grn::
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs )
{
  for (size_t cellI=0; cellI<T.numCell(); ++cellI) { 
    // Threshold
    double u = parameter(0);//h
    
    // Internal contribution
    size_t add=2;
    for(size_t b=0; b<numVariableIndex(1); b++ )
      u += parameter(b + add)*cellData[cellI][variableIndex(1,b)];
    
  // Apply sigmoid and tau parameter
    if(parameter(1)>0.)
      cellDerivs[cellI][variableIndex(0,0)] += sigmoid(u)/parameter(1);
    else {
      std::cerr << "Grn::derivs Division by tau=0." << std::endl;
      exit(-1);
    }
  }
}


Gsrn2::Gsrn2(std::vector<double> &paraValue, 
	     std::vector< std::vector<size_t> > 
	     &indValue ) 
{  
  //Do some checks on the parameters and variable indeces
  //
  if( paraValue.size() !=4 ) {
    std::cerr << "Grns2::Grns2() "
              << "four parameters (p0,... p3) must be given.  reaction 4 3 1 1 1\n";
    exit(0);
  }
  if( indValue.size()!=3 ) {
    std::cerr << "Grns2::Grns2() "
              << " three levels of variable indeces must be used.\n";
    exit(0);
  }
  
  //Set the variable values
  //
  setId("gsrn2");
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  tmp[0] = "h";
  tmp[1] = "tau";
  tmp[2] = "p2";// interacellular
  tmp[3] = "p3";// intercellular

  
  setParameterId( tmp );
}

void Gsrn2::						       
derivs(Tissue &T,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs )		
  
{
  for (size_t cellI=0; cellI<T.numCell(); ++cellI) { 
    
    Cell* cell1 = &(T.cell(cellI));
    size_t numWalls = cell1->numWall();
    
    //Threshold
    double u = parameter(0);//h
    
    //size_t addIndex=2;//avoiding h and tau
    
    // Internal contribution
    u += parameter(2)*cellData[cellI][variableIndex(1,0)];
    
    // Neighbor contributions (direct)
    for (size_t k=0; k<numWalls; ++k){
      Cell* cell2 = cell1->cellNeighbor(k);
      if(cell2!=T.background()){	
        u -= parameter(3)
          *cellData[cell2->index()][variableIndex(2,0)];
      }
    }     
    // Apply sigmoid and tau parameter
    cellDerivs[cellI][variableIndex(0,0)] += sigmoid(u)/parameter(1);
  }
}


// Gsrn2::Gsrn2(std::vector<double> &paraValue, 
// 	     std::vector< std::vector<size_t> > 
// 	     &indValue ) 
// {  
//   //Do some checks on the parameters and variable indeces
//   //
//   if( paraValue.size()<2 ) {
//     std::cerr << "Grns2::Grns2() "
//               << "At least two parameters (h,tau) must be given.\n";
//     exit(0);
//   }
//   if( indValue.size()>3 ) {
//     std::cerr << "Grns2::Grns2() "
//               << "At most three levels of variable indeces can be used.\n";
//     exit(0);
//   }
//   if( 
//      ( indValue.size()==3 && 
//        paraValue.size() != indValue[0].size()+indValue[1].size()+
//        0.5*indValue[2].size()+2 ) ||
//      ( indValue.size()==2 && 
//        paraValue.size() != indValue[0].size()+indValue[1].size()+2 ) ||
//      ( indValue.size()==1 && paraValue.size() != indValue[0].size()+2 ) ||
//      ( indValue.size()==0 && paraValue.size() != 2 ) 
//       ) {
//     std::cerr << "Grns2::Grns2() "
//               << "Wrong number of parameters/variable indeces given.\n";
//     std::cerr << paraValue.size() << " " << indValue.size() << " ";
//     for( size_t i=0 ; i<indValue.size() ; i++ )
//       std::cerr << indValue[i].size() << " ";
//     std::cerr << "\n";
//     exit(0);
//   }
//   if( indValue.size()==3 && indValue[2].size()%2 ) {
//     std::cerr << "Grns2::Grns2() "
//               << "Ligand/Receptor indeces must come in pairs.\n";
//     exit(0);
//   }
  
//   //Set the variable values
//   //
//   setId("gsrn2");
//   setParameter(paraValue);  
//   setVariableIndex(indValue);
  
//   //Set the parameter identities
//   //
//   std::vector<std::string> tmp( numParameter() );
//   tmp[0] = "h";
//   tmp[1] = "tau";
//   size_t add=2;
//   for( size_t i=0 ; i<numVariableIndex(0) ; i++ )
//     tmp[i+add] = "T";
//   if( indValue.size()>1 ) {
//     add += numVariableIndex(0);
//     for( size_t i=0 ; i<numVariableIndex(1) ; i++ )
//       tmp[i+add] = "T_hat";
//   }
//   if( indValue.size()>1 ) {
//     add += numVariableIndex(1);
//     for( size_t i=add ; i<tmp.size() ; i++ )
//       tmp[i] = "T_T";
//   }
//   setParameterId( tmp );
// }

// void Gsrn2::						       
// derivs(Tissue &T,
//        DataMatrix &cellData,
//        DataMatrix &wallData,
//        DataMatrix &vertexData,
//        DataMatrix &cellDerivs,
//        DataMatrix &wallDerivs,
//        DataMatrix &vertexDerivs )		
  
// {
//   //Threshold
//   double u = parameter(0);//h
  
//   size_t addIndex=2;//avoiding h and tau
//   // Internal contribution
//   if( numVariableIndexLevel() )
//     for( size_t b=0 ; b<numVariableIndex(0) ; b++ )
//       u += parameter(addIndex+b)*y[compartment.index()][variableIndex(0,b)];
  
//   // Neighbor contributions (direct)
//   if( numVariableIndexLevel()>1 ) {
//     addIndex += numVariableIndex(0);//after local T parameters
//     for( size_t n=0 ; n<compartment.numNeighbor() ; n++ )
//       for( size_t b=0 ; b<numVariableIndex(1) ; b++ )
//         u += parameter(addIndex+b)
//           *y[compartment.neighbor(n)][variableIndex(1,b)];
//   }
//   // Neighbor contributions ligand receptor version
//   if( numVariableIndexLevel()>2 ) {
//     addIndex += numVariableIndex(1);//after first neighbor T parameters
//     for( size_t n=0 ; n<compartment.numNeighbor() ; n++ ) {
//       size_t parIndexCount=0;
//       for( size_t b=0 ; b<numVariableIndex(2) ; b+=2 )
//         u += parameter(addIndex+parIndexCount++)
//           *y[compartment.index()][variableIndex(2,b)]
//           *y[compartment.neighbor(n)][variableIndex(2,b+1)];
//     }
//   }
//   // Apply sigmoid and tau parameter
//   if( parameter(1)>0. )
//     dydt[compartment.index()][species] += sigmoid(u)/parameter(1);
//   else {
//     std::cerr << "Gsrn2::derivs Division by tau=0.\n";
//     exit(-1);
//   }
// }
