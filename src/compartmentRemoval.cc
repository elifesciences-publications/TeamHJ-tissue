/**
 * Filename     : compartmentRemoval.cc
 * Description  : Classes describing compartment removal updates
 * Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
 * Created      : November 2006
 * Revision     : $Id:$
 */
#include"compartmentRemoval.h"
#include"baseCompartmentChange.h"

//!Constructor
RemovalOutsideRadius::
RemovalOutsideRadius(std::vector<double> &paraValue, 
										 std::vector< std::vector<size_t> > 
										 &indValue ) 
{
	//Do some checks on the parameters and variable indeces
  //////////////////////////////////////////////////////////////////////
  if( paraValue.size()!=1 ) {
    std::cerr << "RemovalOutsideRadius::"
							<< "RemovalOutsideRadius() "
							<< "One parameter used R_threshold\n";
    exit(0);
  }
  if( indValue.size() != 0 ) {
    std::cerr << "RemovalOutsideRadius::"
							<< "RemovalOutsideRadius() "
							<< "No variable index is used.\n";
    exit(0);
  }
  //Set the variable values
  //////////////////////////////////////////////////////////////////////
  setId("RemovalOutsideRadius");
	setNumChange(-1);
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //////////////////////////////////////////////////////////////////////
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "R_threshold";
  setParameterId( tmp );
}

//! Flags a cell for division if its position is outside threshold radius
/*! 
 */
int RemovalOutsideRadius::
flag(Tissue *T,size_t i,
     std::vector< std::vector<double> > &cellData,
     std::vector< std::vector<double> > &wallData,
     std::vector< std::vector<double> > &vertexData,
     std::vector< std::vector<double> > &cellDerivs,
     std::vector< std::vector<double> > &wallDerivs,
     std::vector< std::vector<double> > &vertexDerivs ) {
	
	//Calculate cell center from vertices positions
	std::vector<double> cellCenter;
	cellCenter = T->cell(i).positionFromVertex(vertexData);
	assert( cellCenter.size() == vertexData[0].size() );
	double R=0.0;
	for( size_t d=0 ; d<cellCenter.size() ; ++d )
		R += cellCenter[d]*cellCenter[d];
	R = std::sqrt(R);
  if( R > parameter(0) ) {
    std::cerr << "Cell " << i 
							<< " marked for removal at radial distance " 
							<< R << std::endl;
    return 1;
  } 
  return 0;
}

//! Updates the dividing cell by adding a prependicular wall from the longest
/*! 
 */
void RemovalOutsideRadius::
update(Tissue *T,size_t i,
       std::vector< std::vector<double> > &cellData,
       std::vector< std::vector<double> > &wallData,
       std::vector< std::vector<double> > &vertexData,
       std::vector< std::vector<double> > &cellDeriv,
       std::vector< std::vector<double> > &wallDeriv,
       std::vector< std::vector<double> > &vertexDeriv ) {
  
	//Remove cell and adjust its neighboring walls and vertices
	T->removeCell(i,cellData,wallData,vertexData,cellDeriv,wallDeriv,
								vertexDeriv);
	
	//Check that the removal did not mess up the data structure
	//T->checkConnectivity(1);	
}

//!Constructor
RemovalOutsideRadiusEpidermis::
RemovalOutsideRadiusEpidermis(std::vector<double> &paraValue, 
										 std::vector< std::vector<size_t> > 
										 &indValue ) 
{
	//Do some checks on the parameters and variable indeces
  //////////////////////////////////////////////////////////////////////
  if( paraValue.size()!=2 ) {
    std::cerr << "RemovalOutsideRadiusEpidermis::"
							<< "RemovalOutsideRadiusEpidermis() "
							<< "Two parameters used, R_threshold1 and R_threshold2\n";
    exit(0);
  }
  if( indValue.size() != 0 ) {
    std::cerr << "RemovalOutsideRadiusEpidermis::"
							<< "RemovalOutsideRadiusEpidermis() "
							<< "No variable index is used.\n";
    exit(0);
  }
  //Set the variable values
  //////////////////////////////////////////////////////////////////////
  setId("RemovalOutsideRadiusEpidermis");
	setNumChange(-2);
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //////////////////////////////////////////////////////////////////////
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "R_threshold1";
  tmp[0] = "R_threshold2";
  setParameterId( tmp );
}

//! Flags a cell for removal if its position is outside threshold radius
/*! 
 */
int RemovalOutsideRadiusEpidermis::
flag(Tissue *T,size_t i,
     std::vector< std::vector<double> > &cellData,
     std::vector< std::vector<double> > &wallData,
     std::vector< std::vector<double> > &vertexData,
     std::vector< std::vector<double> > &cellDerivs,
     std::vector< std::vector<double> > &wallDerivs,
     std::vector< std::vector<double> > &vertexDerivs ) {
	
	//Calculate cell center from vertices positions
	std::vector<double> cellCenter;
	cellCenter = T->cell(i).positionFromVertex(vertexData);
	assert( cellCenter.size() == vertexData[0].size() );
	double R=0.0;
	for( size_t d=0 ; d<cellCenter.size() ; ++d )
		R += cellCenter[d]*cellCenter[d];
	R = std::sqrt(R);
  if( R > parameter(0) ) {
    std::cerr << "Epidermal cells outside " << parameter(1) 
							<< " marked for removal due to cell " << i 
							<< " at radial distance " << R << std::endl;
    return 1;
  } 
  return 0;
}

//! Updates the dividing cell by adding a prependicular wall from the longest
/*! 
 */
void RemovalOutsideRadiusEpidermis::
update(Tissue *T,size_t i,
       std::vector< std::vector<double> > &cellData,
       std::vector< std::vector<double> > &wallData,
       std::vector< std::vector<double> > &vertexData,
       std::vector< std::vector<double> > &cellDeriv,
       std::vector< std::vector<double> > &wallDeriv,
       std::vector< std::vector<double> > &vertexDeriv ) {
  
	//Remove cell and adjust its neighboring walls and vertices
	T->removeEpidermalCells(cellData,wallData,vertexData,cellDeriv,
													wallDeriv,vertexDeriv,parameter(1));
	
	//Check that the removal did not mess up the data structure
	//T->checkConnectivity(1);	
}

//!Constructor
RemovalOutsideMaxDistanceEpidermis::
RemovalOutsideMaxDistanceEpidermis(std::vector<double> &paraValue, 
																	 std::vector< std::vector<size_t> > 
																	 &indValue ) 
{
	//Do some checks on the parameters and variable indeces
  //////////////////////////////////////////////////////////////////////
  if( paraValue.size()!=2 ) {
    std::cerr << "RemovalOutsideMaxDistanceEpidermis::"
							<< "RemovalOutsideMaxDistanceEpidermis() "
							<< "Two parameters used, R_threshold1 and R_threshold2\n";
    exit(0);
  }
  if( indValue.size() != 1 || indValue[0].size() != 1 ) {
    std::cerr << "RemovalOutsideMaxDistanceEpidermis::"
							<< "RemovalOutsideMaxDistanceEpidermis() "
							<< "One variable index is used.\n";
    exit(0);
  }
  //Set the variable values
  //////////////////////////////////////////////////////////////////////
  setId("RemovalOutsideMaxDistanceEpidermis");
	setNumChange(-2);
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //////////////////////////////////////////////////////////////////////
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "R_threshold1";
  tmp[0] = "R_threshold2";
  setParameterId( tmp );
}

//! Flags a cell for division if its position is outside threshold radius
/*! 
 */
int RemovalOutsideMaxDistanceEpidermis::
flag(Tissue *T,size_t i,
     std::vector< std::vector<double> > &cellData,
     std::vector< std::vector<double> > &wallData,
     std::vector< std::vector<double> > &vertexData,
     std::vector< std::vector<double> > &cellDerivs,
     std::vector< std::vector<double> > &wallDerivs,
     std::vector< std::vector<double> > &vertexDerivs ) 
{	
	size_t d = variableIndex(0,0);
	assert( d<vertexData[0].size() );

	//Find max position in the specified direction
	max_=vertexData[0][d];
	for( size_t j=1 ; j<vertexData.size() ; ++j )
		if( vertexData[j][d]>max_ )
			max_=vertexData[j][d];
	
	//Calculate cell center from vertex positions
	std::vector<double> cellCenter;
	cellCenter = T->cell(i).positionFromVertex(vertexData);
	assert( cellCenter.size() == vertexData[0].size() );
	double dist=std::fabs( cellCenter[d]-max_ );
  if( dist > parameter(0) ) {
    std::cerr << "Epidermal cells outside distance " << parameter(1) 
							<< " marked for removal due to cell " << i 
							<< " at distance from max " << dist << std::endl;
    return 1;
  } 
  return 0;
}

//! Removes cells 'outside' a distance in a specific direction
/*! 
 */
void RemovalOutsideMaxDistanceEpidermis::
update(Tissue *T,size_t i,
       std::vector< std::vector<double> > &cellData,
       std::vector< std::vector<double> > &wallData,
       std::vector< std::vector<double> > &vertexData,
       std::vector< std::vector<double> > &cellDeriv,
       std::vector< std::vector<double> > &wallDeriv,
       std::vector< std::vector<double> > &vertexDeriv ) {
  
	//Remove cell and adjust its neighboring walls and vertices
	T->removeEpidermalCellsAtDistance(cellData,wallData,vertexData,
																		cellDeriv,wallDeriv,vertexDeriv,
																		parameter(1),max_,variableIndex(0,0));
	
	//Check that the removal did not mess up the data structure
	//T->checkConnectivity(1);	
}

RemovalOutsidePosition::
RemovalOutsidePosition(std::vector<double> &paraValue, 
											 std::vector< std::vector<size_t> > 
											 &indValue ) 
{
	// Do some checks on the parameters and variable indeces
  //
  if (paraValue.size()!=2) {
    std::cerr << "RemovalOutsidePosition::"
							<< "RemovalOutsidePosition() "
							<< "Two parameters used, x_threshold and sign "
							<< "+1(-1) for removal of cells above (below) threshold." << std::endl; 
    exit(0);
  }
  if (paraValue[1]!=1 && paraValue[1]!=-1) {
    std::cerr << "RemovalOutsidePosition::"
							<< "RemovalOutsidePosition() "
							<< "sign (second) parameter has to be +-1." << std::endl;
    exit(0);
  }
	
  if( indValue.size() != 1 || indValue[0].size() != 1 ) {
    std::cerr << "RemovalOutsidePosition::"
							<< "RemovalOutsidePosition() "
							<< "One variable index is used for determining the threshold variable."
							<< std::endl;
    exit(0);
  }
  //Set the variable values
  //////////////////////////////////////////////////////////////////////
  setId("RemovalOutsidePosition");
	setNumChange(-1);
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //////////////////////////////////////////////////////////////////////
  std::vector<std::string> tmp( numParameter() );
  tmp.resize( numParameter() );
  tmp[0] = "threshold";
  tmp[0] = "sign";
  setParameterId( tmp );
}

int RemovalOutsidePosition::
flag(Tissue *T,size_t i,
     std::vector< std::vector<double> > &cellData,
     std::vector< std::vector<double> > &wallData,
     std::vector< std::vector<double> > &vertexData,
     std::vector< std::vector<double> > &cellDerivs,
     std::vector< std::vector<double> > &wallDerivs,
     std::vector< std::vector<double> > &vertexDerivs ) 
{	
	size_t d = variableIndex(0,0);
	assert( d<vertexData[0].size() );
	
	//Calculate cell center from vertex positions
	std::vector<double> cellCenter;
	cellCenter = T->cell(i).positionFromVertex(vertexData);
	assert( cellCenter.size() == vertexData[0].size() );
	double dist = parameter(1)*(cellCenter[d]-parameter(0));
  if( dist > 0.0 ) {
    return 1;
  } 
  return 0;
}

void RemovalOutsidePosition::
update(Tissue *T,size_t i,
       std::vector< std::vector<double> > &cellData,
       std::vector< std::vector<double> > &wallData,
       std::vector< std::vector<double> > &vertexData,
       std::vector< std::vector<double> > &cellDeriv,
       std::vector< std::vector<double> > &wallDeriv,
       std::vector< std::vector<double> > &vertexDeriv ) {
  
	//Remove cell and adjust its neighboring walls and vertices
	//Remove cell and adjust its neighboring walls and vertices
	T->removeCell(i,cellData,wallData,vertexData,cellDeriv,wallDeriv,
								vertexDeriv);
}

