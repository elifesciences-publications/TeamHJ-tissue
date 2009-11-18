/**
 * Filename     : compartmentRemoval.cc
 * Description  : Classes describing compartment removal updates
 * Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
 * Created      : November 2006
 * Revision     : $Id:$
 */
#include"compartmentRemoval.h"
#include"baseCompartmentChange.h"

RemovalIndex::
RemovalIndex(std::vector<double> &paraValue, 
						 std::vector< std::vector<size_t> > 
						 &indValue ) 
{
	//Do some checks on the parameters and variable indeces
  //
  if( paraValue.size()!=0 ) {
    std::cerr << "RemovalIndex::"
							<< "RemovalIndex() "
							<< "No parameters used." << std::endl;
    exit(0);
  }
  if( indValue.size() != 1 ) {
    std::cerr << "RemovalIndex::"
							<< "RemovalIndex() "
							<< "List of cell indices to be removed is used.\n";
    exit(0);
  }
  //Set the variable values
  //
  setId("RemovalIndex");
	setNumChange(-1);
  setParameter(paraValue);  
  setVariableIndex(indValue);
  
  //Set the parameter identities
  //
  std::vector<std::string> tmp( numParameter() );
  setParameterId( tmp );
}

int RemovalIndex::
flag(Tissue *T,size_t i,
     std::vector< std::vector<double> > &cellData,
     std::vector< std::vector<double> > &wallData,
     std::vector< std::vector<double> > &vertexData,
     std::vector< std::vector<double> > &cellDerivs,
     std::vector< std::vector<double> > &wallDerivs,
     std::vector< std::vector<double> > &vertexDerivs ) 
{	
	// Should only be done once!
	static int updateFlag=1;
	if (updateFlag) {
		updateFlag=0;
		return 1;
	}
  return 0;
}

void RemovalIndex::
update(Tissue *T,size_t i,
       std::vector< std::vector<double> > &cellData,
       std::vector< std::vector<double> > &wallData,
       std::vector< std::vector<double> > &vertexData,
       std::vector< std::vector<double> > &cellDeriv,
       std::vector< std::vector<double> > &wallDeriv,
       std::vector< std::vector<double> > &vertexDeriv ) {
  
	//Remove all cells listed and adjust its neighboring walls and vertices
	T->removeCells(variableIndex(0),cellData,wallData,vertexData,cellDeriv,wallDeriv,
								 vertexDeriv);	
}

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



RemovalOutsideRadiusEpidermisMk2::RemovalOutsideRadiusEpidermisMk2(std::vector<double> &paraValue, 
	std::vector< std::vector<size_t> > &indValue) 
{
	//Do some checks on the parameters and variable indeces
	//////////////////////////////////////////////////////////////////////
	if (paraValue.size() != 2) {
		std::cerr << "RemovalOutsideRadiusEpidermisMk2::RemovalOutsideRadiusEpidermis() "
		<< "Two parameters used, R_threshold1 and R_threshold2\n";

		exit(EXIT_FAILURE);
	}
	
	if (indValue.size() != 0) {
		std::cerr << "RemovalOutsideRadiusEpidermis::"
		<< "RemovalOutsideRadiusEpidermis() "
		<< "No variable index is used.\n";
		
		exit(EXIT_FAILURE);
	}

	//Set the variable values
	//////////////////////////////////////////////////////////////////////
	setId("RemovalOutsideRadiusEpidermisMk2");
	setNumChange(-2);
	setParameter(paraValue);  
	setVariableIndex(indValue);
	
	//Set the parameter identities
	//////////////////////////////////////////////////////////////////////
	std::vector<std::string> tmp(numParameter());
	tmp.resize(numParameter());
	tmp[0] = "R_threshold1";
	tmp[0] = "R_threshold2";
	setParameterId(tmp);
}

int RemovalOutsideRadiusEpidermisMk2::flag(Tissue *T, size_t i,
	std::vector< std::vector<double> > &cellData,
	std::vector< std::vector<double> > &wallData,
	std::vector< std::vector<double> > &vertexData,
	std::vector< std::vector<double> > &cellDerivs,
	std::vector< std::vector<double> > &wallDerivs,
	std::vector< std::vector<double> > &vertexDerivs)
{
	//Calculate cell center from vertices positions
	std::vector<double> cellCenter;
	cellCenter = T->cell(i).positionFromVertex(vertexData);
	assert( cellCenter.size() == vertexData[0].size() );
	double R = 0.0;
	for (size_t d = 0; d < cellCenter.size(); ++d) {
		R += cellCenter[d]*cellCenter[d];
	}
	R = std::sqrt(R);
	if (R > parameter(0)) {
		std::cerr << "Epidermal cells outside " << parameter(1) 
		<< " marked for removal due to cell " << i 
		<< " at radial distance " << R << std::endl;

		return 1;
	} 
	return 0;
}

void RemovalOutsideRadiusEpidermisMk2::update(Tissue *T, size_t i,
	std::vector< std::vector<double> > &cellData,
	std::vector< std::vector<double> > &wallData,
	std::vector< std::vector<double> > &vertexData,
	std::vector< std::vector<double> > &cellDeriv,
	std::vector< std::vector<double> > &wallDeriv,
	std::vector< std::vector<double> > &vertexDeriv)
{
	//Remove cell and adjust its neighboring walls and vertices
	T->removeEpidermalCellsMk2(cellData, wallData, vertexData, cellDeriv, wallDeriv, vertexDeriv, parameter(1));
	
	//Check that the removal did not mess up the data structure
	//T->checkConnectivity(1);	
}






RemovalWholeCellOutsideRadiusEpidermis::RemovalWholeCellOutsideRadiusEpidermis(std::vector<double> &paraValue, std::vector< std::vector<size_t> > &indValue) 
{
	if (paraValue.size() != 2)
	{
		std::cerr << "RemovalWholeCellOutsideRadiusEpidermis::RemovalWholeCellOutsideRadiusEpidermis() "
		<< "Two parameter used R_threshold and R_threshold2\n";
		std::exit(EXIT_FAILURE);
	}
	if (indValue.size() != 0)
	{
		std::cerr << "RemovalWholeCellOutsideRadiusEpidermis::RemovalWholeCellOutsideRadiusEpidermis() "
		<< "No variable index is used.\n";
		std::exit(EXIT_FAILURE);
	}

	setId("RemovalWholeCellOutsideRadiusEpidermis");
	setNumChange(-1);
	setParameter(paraValue);  
	setVariableIndex(indValue);
  
	std::vector<std::string> tmp(numParameter());
	tmp.resize(numParameter());
	tmp[0] = "R_threshold";
	tmp[1] = "R_threshold2";
	setParameterId( tmp );
}

int RemovalWholeCellOutsideRadiusEpidermis::flag(Tissue *T, size_t i, std::vector< std::vector<double> > &cellData,
	std::vector< std::vector<double> > &wallData,
	std::vector< std::vector<double> > &vertexData,
	std::vector< std::vector<double> > &cellDerivs,
	std::vector< std::vector<double> > &wallDerivs,
	std::vector< std::vector<double> > &vertexDerivs)
{
	Cell &cell = T->cell(i);

	if (checkIfCellIsOutside(cell, vertexData, parameter(0)))
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

void RemovalWholeCellOutsideRadiusEpidermis::update(Tissue *T, size_t i, std::vector< std::vector<double> > &cellData,
	std::vector< std::vector<double> > &wallData,
	std::vector< std::vector<double> > &vertexData,
	std::vector< std::vector<double> > &cellDeriv,
	std::vector< std::vector<double> > &wallDeriv,
	std::vector< std::vector<double> > &vertexDeriv)
{
	T->removeEpidermalCells(cellData, wallData, vertexData, cellDeriv, wallDeriv, vertexDeriv, parameter(1), false);
}


bool RemovalWholeCellOutsideRadiusEpidermis::checkIfCellIsOutside(Cell &cell, std::vector< std::vector<double> > &vertexData, const double radius) const
{	
	const size_t dimensions = vertexData[0].size();
	const double radius2 = std::pow(radius, 2.0);

	for (size_t i = 0; i < cell.numVertex(); ++i)
	{
		const Vertex &vertex = *cell.vertex(i);
		const size_t vertexIndex = vertex.index();

		double sum2 = 0.0;
		
		for (size_t dimension = 0; dimension < dimensions; ++dimension)
		{
			sum2 += std::pow(vertexData[vertexIndex][dimension], 2.0);
		}

		if (sum2 < radius2)
		{
			return false;
		}
	}

	return true;
}



RemovalConcaveCellsAtEpidermis::RemovalConcaveCellsAtEpidermis(std::vector<double> &paraValue, std::vector< std::vector<size_t> > &indValue) 
{
	if (paraValue.size() != 0)
	{
		std::cerr << "RemovalConcaveCellsAtEpidermis::RemovalConcaveCellsAtEpidermis() uses no parameters.\n";
		std::exit(EXIT_FAILURE);
	}
	if (indValue.size() != 0)
	{
		std::cerr << "RemovalConcaveCellsAtEpidermis::RemovalConcaveCellsAtEpidermis() "
		<< "No variable index is used.\n";
		std::exit(EXIT_FAILURE);
	}

	setId("RemovalConcaveCellsAtEpidermis");
	setNumChange(-1);
	setParameter(paraValue);  
	setVariableIndex(indValue);
  
	std::vector<std::string> tmp(numParameter());
	tmp.resize(numParameter());
	setParameterId(tmp);
}

int RemovalConcaveCellsAtEpidermis::flag(Tissue *T, size_t i, std::vector< std::vector<double> > &cellData,
	std::vector< std::vector<double> > &wallData,
	std::vector< std::vector<double> > &vertexData,
	std::vector< std::vector<double> > &cellDerivs,
	std::vector< std::vector<double> > &wallDerivs,
	std::vector< std::vector<double> > &vertexDerivs)
{
	Cell &cell = T->cell(i);

	for (size_t k = 0; k < cell.numWall(); ++k)
	{
		if (cell.cellNeighbor(k) == T->background() && cell.isConcave(vertexData))
		{
			return 1;
		}
	}

	return 0;
}

void RemovalConcaveCellsAtEpidermis::update(Tissue *T, size_t i, std::vector< std::vector<double> > &cellData,
	std::vector< std::vector<double> > &wallData,
	std::vector< std::vector<double> > &vertexData,
	std::vector< std::vector<double> > &cellDeriv,
	std::vector< std::vector<double> > &wallDeriv,
	std::vector< std::vector<double> > &vertexDeriv)
{
	std::cerr << "Removing concave cell " << i << " at boundary.\n";
	T->removeCell(i, cellData, wallData, vertexData, cellDeriv, wallDeriv, vertexDeriv);	
}


RemoveIsolatedCells::RemoveIsolatedCells(std::vector<double> &paraValue, std::vector< std::vector<size_t> > &indValue) 
{
	if (paraValue.size() != 0)
	{
		std::cerr << "RemoveIsolatedCells::RemoveIsolatedCells() uses no parameters.\n";
		std::exit(EXIT_FAILURE);
	}
	if (indValue.size() != 0)
	{
		std::cerr << "RemoveIsolatedCells::RemoveIsolatedCells() "
		<< "No variable index is used.\n";
		std::exit(EXIT_FAILURE);
	}

	setId("RemoveIsolatedCells");
	setNumChange(-1);
	setParameter(paraValue);  
	setVariableIndex(indValue);
  
	std::vector<std::string> tmp(numParameter());
	tmp.resize(numParameter());
	setParameterId(tmp);
}

int RemoveIsolatedCells::flag(Tissue *T, size_t i, std::vector< std::vector<double> > &cellData,
	std::vector< std::vector<double> > &wallData,
	std::vector< std::vector<double> > &vertexData,
	std::vector< std::vector<double> > &cellDerivs,
	std::vector< std::vector<double> > &wallDerivs,
	std::vector< std::vector<double> > &vertexDerivs)
{
	if (T->numCell() == 1)
	{
		return 0;
	}

	Cell &cell = T->cell(i);

	for (size_t k = 0; k < cell.numWall(); ++k)
	{
		if (cell.cellNeighbor(k) != T->background())
		{
			return 0;
		}
	}

	return 1;
}

void RemoveIsolatedCells::update(Tissue *T, size_t i, std::vector< std::vector<double> > &cellData,
	std::vector< std::vector<double> > &wallData,
	std::vector< std::vector<double> > &vertexData,
	std::vector< std::vector<double> > &cellDeriv,
	std::vector< std::vector<double> > &wallDeriv,
	std::vector< std::vector<double> > &vertexDeriv)
{
	std::cerr << "Removing isolated cell " << i << " at boundary.\n";
	T->removeCell(i, cellData, wallData, vertexData, cellDeriv, wallDeriv, vertexDeriv);	
}



RemoveFoldedCells::RemoveFoldedCells(std::vector<double> &paraValue, std::vector< std::vector<size_t> > &indValue) 
{
	if (paraValue.size() != 1)
	{
		std::cerr << "RemoveFoldedCells::RemoveFoldedCells() uses one parameter.\n";
		std::cerr << "1 = Only cells on boundary, 0 = All folded cells.\n";
		std::exit(EXIT_FAILURE);
	}
	if (indValue.size() != 0)
	{
		std::cerr << "RemoveFoldedCells::RemoveFoldedCells() No variable index is used.\n";
		std::exit(EXIT_FAILURE);
	}
	
	setId("RemoveFoldedCells");
	setNumChange(-1);
	setParameter(paraValue);  
	setVariableIndex(indValue);
	
	std::vector<std::string> tmp(numParameter());
	tmp.resize(numParameter());
	tmp[0] = "boundary_flag";
	setParameterId(tmp);
}

int RemoveFoldedCells::flag(Tissue *T, size_t i, std::vector< std::vector<double> > &cellData,
	std::vector< std::vector<double> > &wallData,
	std::vector< std::vector<double> > &vertexData,
	std::vector< std::vector<double> > &cellDerivs,
	std::vector< std::vector<double> > &wallDerivs,
	std::vector< std::vector<double> > &vertexDerivs)
{
	Cell &cell = T->cell(i);

	if (parameter(0) && !cell.isNeighbor(T->background()))
	{
		return 0;
	}

	if (cell.isFolded(vertexData))
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

void RemoveFoldedCells::update(Tissue *T, size_t i, std::vector< std::vector<double> > &cellData,
	std::vector< std::vector<double> > &wallData,
	std::vector< std::vector<double> > &vertexData,
	std::vector< std::vector<double> > &cellDeriv,
	std::vector< std::vector<double> > &wallDeriv,
	std::vector< std::vector<double> > &vertexDeriv)
{
	Cell &cell = T->cell(i);

	if (cell.isNeighbor(T->background()))
	{
		std::cerr << "Removing folded cell on boundary.\n";
	}
	else
	{
		std::cerr << "Removing folded cell not on boundary.\n";
	}

	T->removeCell(i, cellData, wallData, vertexData, cellDeriv, wallDeriv, vertexDeriv);	
}


