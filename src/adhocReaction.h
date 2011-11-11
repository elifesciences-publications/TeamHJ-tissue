/**
 * Filename     : adhocReaction.h
 * Description  : Classes describing some ad hoc updates
 * Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
 * Created      : September 2007
 * Revision     : $Id:$
 */
#ifndef ADHOCREACTION_H
#define ADHOCREACTION_H

#include"tissue.h"
#include"baseReaction.h"
#include<cmath>

//!Sets positional derivatives to zero for vertices in specified region  
class VertexNoUpdateFromPosition : public BaseReaction {
  
 public:
  
  VertexNoUpdateFromPosition(std::vector<double> &paraValue, 
			     std::vector< std::vector<size_t> > 
			     &indValue );
  
  void derivs(Tissue &T,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      DataMatrix &cellDerivs,
	      DataMatrix &wallDerivs,
	      DataMatrix &vertexDerivs );
};

///
/// @brief Sets positional derivatives to zero for vertices at boundary  
///
class VertexNoUpdateBoundary : public BaseReaction {
  
 public:
  
  VertexNoUpdateBoundary(std::vector<double> &paraValue, 
			 std::vector< std::vector<size_t> > 
			 &indValue );
  
  void derivs(Tissue &T,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      DataMatrix &cellDerivs,
	      DataMatrix &wallDerivs,
	      DataMatrix &vertexDerivs );
};

///
/// @brief Moves the complete tissue such that the maximal value in specified direction is constant
///
/// The translation is done at each update (i.e. after each ODE integration step)
///
/// p_0 maxPos
/// vI_00 positional index
///  
class VertexTranslateToMax : public BaseReaction {
  
 public:
  
  VertexTranslateToMax(std::vector<double> &paraValue, 
		       std::vector< std::vector<size_t> > 
		       &indValue );
  
  void initiate(Tissue &T,
		DataMatrix &cellData,
		DataMatrix &wallData,
		DataMatrix &vertexData);
  
  void derivs(Tissue &T,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      DataMatrix &cellDerivs,
	      DataMatrix &wallDerivs,
	      DataMatrix &vertexDerivs );
  
  void update(Tissue &T,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      double h);
};

///
/// @brief Centers the tissue such that the center of mass is in origo
///
/// The translation is done at each update (i.e. after each ODE integration step)
///  
class CenterCOM : public BaseReaction
{
public:
  
  CenterCOM(std::vector<double> &paraValue, std::vector< std::vector<size_t> > &indValue);
	
  void initiate(Tissue &T,
		DataMatrix &cellData,
		DataMatrix &wallData,
		DataMatrix &vertexData);
  
  void derivs(Tissue &T,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      DataMatrix &cellDerivs,
	      DataMatrix &wallDerivs,
	      DataMatrix &vertexDerivs);
  
  void update(Tissue &T,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      double h);
};

///
/// @brief Helper reaction that only calculates the PCA plane for every cell.
///
/// Since calculating the PCA plane is computationally expensive, this
/// reaction is calculating it once for all, such that other reaction in need
/// of the PCA plane does not have to do the calculation. 
/// 
/// Note that this reaction has to be before the reactions using the PCA
/// plane.
///
class CalculatePCAPlane : public BaseReaction {
  
public:
  
  CalculatePCAPlane(std::vector<double> &paraValue, 
		    std::vector< std::vector<size_t> > 
		    &indValue );
  
  void initiate(Tissue &T,
		DataMatrix &cellData,
		DataMatrix &wallData,
		DataMatrix &vertexData);
  
  void derivs(Tissue &T,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      DataMatrix &cellDerivs,
	      DataMatrix &wallDerivs,
	      DataMatrix &vertexDerivs );
  
  void update(Tissue &T,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      double h);
};

///
/// @brief Initiate the wall length variables to a factor times the distance between the two vertices
///
/// Sets the wall length variables to be the distance between the two vertices times a factor.
/// 
class InitiateWallLength : public BaseReaction {
  
public:
  
  InitiateWallLength(std::vector<double> &paraValue, 
		     std::vector< std::vector<size_t> > 
		     &indValue );
  
  void initiate(Tissue &T,
		DataMatrix &cellData,
		DataMatrix &wallData,
		DataMatrix &vertexData);
  
  void derivs(Tissue &T,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      DataMatrix &cellDerivs,
	      DataMatrix &wallDerivs,
	      DataMatrix &vertexDerivs );
};

///
/// @brief Adds additional vertices to all walls
///
/// Wall 'meshing' is applied by inserting additional vertices for all walls.
/// This creates new vertices connected to one of the previous walls that are 
/// now divided into several walls. Only a single parameter is given which
/// sets the number of new vertices per wall. The 'remeshing' is done in the 
/// initiation step and there is no contribution to the derivative for this 
/// reaction.
///
/// In the model file the reaction is defined as:
///
/// @verbatim
/// InitiateWallMesh 1 0
/// numVertex
/// @endverbatim 
/// 
class InitiateWallMesh : public BaseReaction {
  
 public:
  
  InitiateWallMesh(std::vector<double> &paraValue, 
		   std::vector< std::vector<size_t> > 
		   &indValue );
  
  void initiate(Tissue &T,
		DataMatrix &cellData,
		DataMatrix &wallData,
		DataMatrix &vertexData);
  
  void derivs(Tissue &T,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      DataMatrix &cellDerivs,
	      DataMatrix &wallDerivs,
	      DataMatrix &vertexDerivs );
};

///
/// @brief Different well-defined derivatives used to check the strain direction update
///
class StrainTest : public BaseReaction {
  
public:
  
  StrainTest(std::vector<double> &paraValue, 
	     std::vector< std::vector<size_t> > 
	     &indValue );
  
  void initiate(Tissue &T,
		DataMatrix &cellData,
		DataMatrix &wallData,
		DataMatrix &vertexData);
	
  void derivs(Tissue &T,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      DataMatrix &cellDerivs,
	      DataMatrix &wallDerivs,
	      DataMatrix &vertexDerivs );
};

class CalculateVertexStressDirection : public BaseReaction
{
public:
	
	CalculateVertexStressDirection(std::vector<double> &paraValue, 
		std::vector< std::vector<size_t> > 
		&indValue );
	
	void initiate(Tissue &T,
		DataMatrix &cellData,
		DataMatrix &wallData,
		DataMatrix &vertexData);
	
	void derivs(Tissue &T,
		DataMatrix &cellData,
		DataMatrix &wallData,
		DataMatrix &vertexData,
		DataMatrix &cellDerivs,
		DataMatrix &wallDerivs,
		DataMatrix &vertexDerivs );
	
	void update(Tissue &T,
		DataMatrix &cellData,
		DataMatrix &wallData,
		DataMatrix &vertexData,
		double h);

private:
	std::vector<size_t> wallForceIndexes_;
};

#endif //ADHOCREACTION_H
