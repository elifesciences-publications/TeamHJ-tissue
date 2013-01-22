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

///
/// @brief Sets positional derivatives to zero for vertices in specified region  
///
/// A threshold is specified in a specific dimension and vertices on either side of this
/// is not updated.
///
/// In the model file, the reaction is specified as:
///
/// @verbatim
/// VertexNoUpdateFromPosition 2 1 1 
/// threshold directionFlag
/// vertexPositionIndex
/// @endverbatim
///
/// where threshold sets the threshold in the dimension (x,y,z) specified with vertexPositionIndex
/// and directionFlag is set to 1 if vertices above the threshold are to be kept fixed and -1
/// for vertices below the threshold.
///
/// @note This function sets the derivatives to zero, which means it has to be provided after
/// reactions that update the vertex derivatives.
/// 
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
/// @brief Sets positional derivatives to zero for vertices with listed indices
///
/// A list of vertex indices are specified for which vertex positions are not updated.
///
/// In the model file, the reaction is specified as:
///
/// @verbatim
/// VertexNoUpdateFromIndex 0 1 N 
/// vertexIndex1 [vertexIndex2...vertexIndexN]
/// @endverbatim
///
/// where the list if indices are the vertices not to be updated.
///
/// @note This function sets the derivatives to zero, which means it has to be provided after
/// reactions that update the vertex derivatives.
/// 
class VertexNoUpdateFromIndex : public BaseReaction {
  
 public:
  
  VertexNoUpdateFromIndex(std::vector<double> &paraValue, 
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
/// @brief Sets positional derivatives to zero for vertices at boundary so that boundary vertices would be restricted from moving in x and/or y ... direction(s) 
///
/// In the model file, the reaction is specified as:
///
/// @verbatim
/// for holding the boundary vertices in all directions:
/// VertexNoUpdateBoundary 0 0
/// 
/// or for holding them only in x and z directions
/// VertexNoUpdateBoundary 0 1 2
/// 0 2 
///
/// @endverbatim
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
		DataMatrix &vertexData,
		DataMatrix &cellDerivs,
		DataMatrix &wallDerivs,
		DataMatrix &vertexDerivs );
  
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
/// The translation is done at each update (i.e. after each ODE integration step), mainly for
/// plotting by keeping the tissue centered.
///
/// In a model file the reaction is given by:
///
/// @verbatim
/// CenterCOM 0 0
/// @endverbatim
///  
class CenterCOM : public BaseReaction
{
public:
  
  CenterCOM(std::vector<double> &paraValue, std::vector< std::vector<size_t> > &indValue);
	
  void initiate(Tissue &T,
		DataMatrix &cellData,
		DataMatrix &wallData,
		DataMatrix &vertexData,
		DataMatrix &cellDerivs,
		DataMatrix &wallDerivs,
		DataMatrix &vertexDerivs);
  
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
/// @brief Centers the tissue such that the center of mass is in origo including for central mesh points
///
/// The translation is done at each update (i.e. after each ODE integration step).
/// This function should be used in connection with mechanical TRBScenterTriangulation
/// since those reactions adds an extra vertex in the cell data vector, which also 
/// have to be moved. In a model file the reaction is defined as
///
/// @verbatim
/// CenterCOMcenterTriangulation 0 1 1 
/// cellDataPositionIndex
/// @endverbatim
///
/// where the index have to be matched with what is given in the mechanical TRBS reaction.
///  
class CenterCOMcenterTriangulation : public BaseReaction
{
public:
  
  CenterCOMcenterTriangulation(std::vector<double> &paraValue, std::vector< std::vector<size_t> > &indValue);
	
  void initiate(Tissue &T,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      DataMatrix &cellDerivs,
	      DataMatrix &wallDerivs,
	      DataMatrix &vertexDerivs);
  
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
		DataMatrix &vertexData,
		DataMatrix &cellDerivs,
		DataMatrix &wallDerivs,
		DataMatrix &vertexDerivs );
  
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
		DataMatrix &vertexData,
		DataMatrix &cellDerivs,
		DataMatrix &wallDerivs,
		DataMatrix &vertexDerivs );
	
  
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
		DataMatrix &vertexData,
		DataMatrix &cellDerivs,
		DataMatrix &wallDerivs,
		DataMatrix &vertexDerivs );
  
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
		DataMatrix &vertexData,
		DataMatrix &cellDerivs,
		DataMatrix &wallDerivs,
		DataMatrix &vertexDerivs );
	
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
		DataMatrix &vertexData,
		DataMatrix &cellDerivs,
		DataMatrix &wallDerivs,
		DataMatrix &vertexDerivs );
  
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




class MoveVerticesRandomlyCapCylinder : public BaseReaction
{
public:
	
	MoveVerticesRandomlyCapCylinder(std::vector<double> &paraValue, 
		std::vector< std::vector<size_t> > 
		&indValue );
	
	void initiate(Tissue &T,
		      DataMatrix &cellData,
		      DataMatrix &wallData,
		      DataMatrix &vertexData,
		      DataMatrix &cellDerivs,
		      DataMatrix &wallDerivs,
		      DataMatrix &vertexDerivs );
	
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










#endif //ADHOCREACTION_H
