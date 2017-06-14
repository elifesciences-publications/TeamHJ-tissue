//
// Filename     : adhocReaction.h
// Description  : Classes describing some ad hoc updates
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : September 2007
// Revision     : $Id:$
//
#ifndef ADHOCREACTION_H
#define ADHOCREACTION_H

#include"tissue.h"
#include"baseReaction.h"
#include "myRandom.h"

#include<cmath>

///
/// @brief Sets positional derivatives to zero for vertices in specified region  
///
/// @details A threshold is specified in a specific dimension and vertices on either side of this
/// is not updated.
///
/// In the model file, the reaction is specified as:
/// @verbatim
/// VertexNoUpdateFromPosition 2 1 1 
/// threshold directionFlag
/// vertexPositionIndex
/// @endverbatim
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
/// @details A list of vertex indices are specified for which vertex positions are not 
/// In the model file, the reaction is specified as:
/// @verbatim
/// VertexNoUpdateFromIndex 0 1 N 
/// vertexIndex1 [vertexIndex2...vertexIndexN]
/// @endverbatim
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
/// @brief Sets positional derivatives to zero for vertices with listed indices
///
/// @details Freezes all of the vertices but the leading ones(at the tip). 
/// In the model file, the reaction is specified as:
/// @verbatim
/// VertexNoUpdateFromList 0 0
/// @endverbatim
/// where the list if indices are the vertices not to be updated.
///
/// @note This function sets the derivatives to zero, which means it has to be provided after
/// reactions that update the vertex derivatives.
/// 
class VertexNoUpdateFromList : public BaseReaction {
  
 public:
  
  VertexNoUpdateFromList(std::vector<double> &paraValue, 
			     std::vector< std::vector<size_t> > 
			     &indValue );
  
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

  std::vector<size_t> updateVertices;
};


///
/// @brief randomizes the growth direction of tip cells(for two dim only)
///
/// @details Freezes all of the vertices but the leading ones(at the tip). 
/// In the model file, the reaction is specified as:
/// @verbatim
/// VertexRandTip 1 1 2
/// max_randomization_angle
/// tip signal index
/// mt vector starting index(2d)
/// @endverbatim
/// where the list if indices are the vertices not to be updated.
///
/// @note This function sets the derivatives to zero, which means it has to be provided after
/// reactions that update the vertex derivatives.
/// 
class VertexRandTip : public BaseReaction {
  
 public:
  
  VertexRandTip(std::vector<double> &paraValue, 
			     std::vector< std::vector<size_t> > 
			     &indValue );
  
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
/// @brief Sets positional derivatives to zero for vertices at boundary so that boundary 
/// vertices would be restricted from moving in x and/or y ... direction(s) 
///
/// @details In the model file, the reaction is specified as:
/// @verbatim
/// for holding the boundary vertices in all directions:
/// VertexNoUpdateBoundary 0 0
/// 
/// or for holding them only in x and z directions
/// VertexNoUpdateBoundary 0 1 2
/// 0 2 
///
/// @endverbatim
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
/// @brief Sets positional derivatives to zero for vertices at boundary so that boundary 
/// vertices would be restricted from moving parallel to the template edges. The boundary 
/// condition gets updated in each step.(not completely debugged)
///
/// @details In the model file, the reaction is specified as:
/// @verbatim
/// for holding the boundary vertices in all directions:
/// VertexNoUpdateBoundaryPtemplate 0 1 1
/// comIndex-cellVector
/// 
/// @endverbatim
///
class VertexNoUpdateBoundaryPtemplate : public BaseReaction { // BB
  
 public:
  
  VertexNoUpdateBoundaryPtemplate(std::vector<double> &paraValue, 
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
/// @brief Sets positional derivatives to zero for vertices at boundary so that boundary 
/// vertices would be restricted from moving parallel to the initial template edges 
/// the boundary is static. 
///
/// @details In the model file, the reaction is specified as:
/// @verbatim
/// for holding the boundary vertices in all directions:
/// VertexNoUpdateBoundaryPtemplateStatic 0 1 1
/// comIndex-cellVector
/// 
/// @endverbatim
///
class VertexNoUpdateBoundaryPtemplateStatic : public BaseReaction { // BB
  
public:
  
  VertexNoUpdateBoundaryPtemplateStatic(std::vector<double> &paraValue, 
                                        std::vector< std::vector<size_t> > 
                                        &indValue );
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
	      DataMatrix &vertexDerivs );
private:
  std::vector<size_t> boundaryVertices;
  std::vector< std::vector<double> > boundaryNormal;
  size_t numBoundaryVertices;
};





///  
/// @brief Sets positional derivatives to zero for vertices at boundary so that boundary 
/// vertices would be restricted from moving parallel to the initial template edges 
/// the boundary is static. 
///
/// @details In the model file, the reaction is specified as:
/// @verbatim
/// for holding the boundary vertices in all directions:
/// VertexNoUpdateBoundaryPtemplateStatic3D 0 1 1
/// comIndex-cellVector
/// 
/// @endverbatim
///
class VertexNoUpdateBoundaryPtemplateStatic3D : public BaseReaction { // BB
  
public:
  
  VertexNoUpdateBoundaryPtemplateStatic3D(std::vector<double> &paraValue, 
                                        std::vector< std::vector<size_t> > 
                                        &indValue );
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
	      DataMatrix &vertexDerivs );
private:
  //  std::vector<size_t> bottomVertices;
  //  std::vector<size_t> sideVertices;
 
  std::vector< std::vector<double> > bottomNormals;
  std::vector< std::vector<double> > sideNormals;
  size_t numBottomCells;
  size_t numSideCells;
};




///  
/// @brief Sets positional derivatives to zero for vertices at boundary so that boundary 
/// vertices would be restricted from moving parallel to the initial template edges 
/// the boundary is static. 
///
/// @details In the model file, the reaction is specified as:
/// @verbatim
/// for holding the boundary vertices in all directions:
/// VertexNoUpdateBoundary3D 0 1 2
/// comIndex-cellVector
/// neighborIndex
/// 
/// @endverbatim
///
class VertexNoUpdateBoundary3D : public BaseReaction { // BB
  
public:
  
  VertexNoUpdateBoundary3D(std::vector<double> &paraValue, 
                                        std::vector< std::vector<size_t> > 
                                        &indValue );
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
	      DataMatrix &vertexDerivs );
private:
   
  std::vector< double > bottomVertices, sideVertices;
  
};




///  
/// @brief Sets positional derivatives to zero for vertices at boundary so that boundary 
/// vertices would be restricted from moving parallel to the initial template edges 
/// the boundary is static. 
///
/// @details In the model file, the reaction is specified as:
/// @verbatim
/// for holding the boundary vertices in all directions:
/// VertexFromConstStressBoundary 8 0 
/// stress x
/// stress y
/// dx
/// dy
/// boundary right  x 
/// boundary left   x 
/// boundary top    y
/// boundary bottom y
/// 
/// @endverbatim
///
class VertexFromConstStressBoundary : public BaseReaction { // BB
  
public:
  
  VertexFromConstStressBoundary(std::vector<double> &paraValue, 
                                        std::vector< std::vector<size_t> > 
                                        &indValue );
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
	      DataMatrix &vertexDerivs );
  void update(Tissue &T,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      double h);
private:
  
  std::vector< std::vector<double> > rightVertices, leftVertices, topVertices, bottomVertices;
  size_t numOldVertices;
  double totaltime;
};





///  
/// @brief for ad-hoc manipulation of the templates

/// @details In the model file, the reaction is specified as:
/// @verbatim
/// for holding the boundary vertices in all directions:
///
/// manipulate 0 0 
/// 
/// @endverbatim
///
class manipulate : public BaseReaction { // BB
  
public:
  
  manipulate(std::vector<double> &paraValue, 
             std::vector< std::vector<size_t> > 
             &indValue );
  
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
	      DataMatrix &vertexDerivs );
  void update(Tissue &T,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      double h);

};



///  
/// @brief Calculates cell polarity vector based on a vector and a measure for each 
/// cell wall and places it in the cell center. The vector field could be printed 
/// via a seperate vtk file using its own printState function.
///
/// @details In the model file, the reaction is specified as:
/// @verbatim
/// for holding the boundary vertices in all directions:
/// cellPolarity3D 0 2 2 2
///
/// CELL_IDENTITY_INDEX
/// comIndex
///
/// polarity_vector_index
/// polarity_measure_INDEX
/// 
/// @endverbatim
///
class cellPolarity3D : public BaseReaction { // BB
  
public:
  
  cellPolarity3D(std::vector<double> &paraValue, 
                 std::vector< std::vector<size_t> > 
                 &indValue );
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
	      DataMatrix &vertexDerivs );
  void update(Tissue &T,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      double h);
  void printState(Tissue *T,
                  DataMatrix &cellData,
                  DataMatrix &wallData,
                  DataMatrix &vertexData, 
                  std::ostream &os); 
    private:
   
  std::vector< std::vector<double> > cellFaces;
  std::vector< std::vector<double> > cellCentPol;
  
};


///  
/// @brief Diffusion in 3 dimensional templates assuming the same volume of the cells (for now ad-hoc) 
///
/// @details In the model file, the reaction is specified as:
/// @verbatim
/// 
/// diffusion3D 1 1 3
/// diff-constant
/// 
/// conc_index
/// neighborWall_index
/// 3Dcell_index
/// 
/// @endverbatim
///
class diffusion3D : public BaseReaction { // BB
  
public:
  
  diffusion3D(std::vector<double> &paraValue, 
                                        std::vector< std::vector<size_t> > 
                                        &indValue );
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
	      DataMatrix &vertexDerivs );
private:
   
  std::vector< std::vector<double> > Cells3d; //holds the wall_indices and neighbohrhood info
  std::vector< std::vector<double> > faceArea; //holds the area of the faces
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
/// @details The translation is done at each update (i.e. after each ODE integration step), 
/// mainly for plotting by keeping the tissue centered.
///
/// In a model file the reaction is given by:
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
/// @brief Centers the tissue such that the center of mass is in origo including for central mesh points.
///
/// @details The translation is done at each update (i.e. after each ODE integration step).
/// This function should be used in connection with mechanical TRBScenterTriangulation
/// since those reactions adds an extra vertex in the cell data vector, which also 
/// have to be moved. In a model file the reaction is defined as
/// @verbatim
/// CenterCOMcenterTriangulation 0 1 1 
/// cellDataPositionIndex
/// @endverbatim
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
/// @details Sets the wall length variables to be the distance between the two vertices times a factor.
/// In the model file the reaction is defined as:
/// @verbatim
/// InitiateWallLength 1 0
/// factor
/// @endverbatim 
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
/// @details Wall 'meshing' is applied by inserting additional vertices for all walls.
/// This creates new vertices connected to one of the previous walls that are 
/// now divided into several walls. Only a single parameter is given which
/// sets the number of new vertices per wall. The 'remeshing' is done in the 
/// initiation step and there is no contribution to the derivative for this 
/// reaction.
///
/// In the model file the reaction is defined as:
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

///
/// @brief scales the template by a factor via Initiate 
/// position of vertices and wall length variables will be scaled.
/// 
/// @details In the model file the reaction is defined as:
/// @verbatim
/// scaleTemplate 1 0
/// factor
/// @endverbatim 
/// 
class scaleTemplate : public BaseReaction
{
public:
	
	scaleTemplate(std::vector<double> &paraValue, 
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
/// @brief scales the template by a factor via Initiate 
/// copies vectors from one index to another in the cell vector
/// ( 4 component after the indices will be copied)
///
/// @details In the model file the reaction is defined as:
/// @verbatim
/// copyCellVector 0 1 2
/// copy_from_index
/// copy_to_index
/// @endverbatim 
/// 
class copyCellVector : public BaseReaction
{
public:
	
	copyCellVector(std::vector<double> &paraValue, 
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
/// @brief scales the template by a factor via Initiate 
/// copies vectors from one index to another in the cell vector
/// ( 4 component after the indices will be copied)
///
/// @details In the model file the reaction is defined as:
/// @verbatim
/// limitZdis 0 1 2
/// copy_from_index
/// copy_to_index
/// @endverbatim 
/// 
class limitZdis : public BaseReaction
{
public:
  
  limitZdis(std::vector<double> &paraValue, 
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
   
  std::vector< size_t > topVertices; //holds the vertex indices at the top
  std::vector< size_t > bottomVertices; //holds the vertex indices at the bottom
};



///
/// @brief scales the template by a factor via Initiate 
/// randomizes the MT direction of the cells within the cell plane
///
/// @details In the model file the reaction is defined as:
/// @verbatim
/// randomizeMT 0 1 1
/// MT_index
/// @endverbatim 
/// 
class randomizeMT : public BaseReaction
{
public:
  
  randomizeMT(std::vector<double> &paraValue, 
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

class restrictVertexRadially : public BaseReaction
{
public:
	
	restrictVertexRadially(std::vector<double> &paraValue, 
		std::vector< std::vector<size_t> > 
		&indValue );
	
	
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
/// @brief A molecule is produced/created with constant rate in the cells 
/// located around CZ to creat phyllotactic pattern.
///
/// 
///
/// @f[ \frac{dc}{dt} = k_c/Area @f] if cell index is in the given list
/// Area will be  one if constant concentration(instead of amount) is needed to be produced
///
/// In a model file the reaction is defined as
///
/// @verbatim
/// creationPrimordiaTime 5 1 2 
/// k_c # creation rate
/// tp # time delay
/// teta # delta teta
/// z # vertical distance from max
/// constant-amount-flag
///
/// c_index
/// com_index
/// @endverbatim
///
class CreationPrimordiaTime: public BaseReaction {
  
  
private: 
  std::vector<size_t> proCells;
  
public:
  
  ///
  /// @brief Main constructor
  ///
  /// This is the main constructor which sets the parameters and variable
  /// indices that defines the reaction.
  ///
  /// @param paraValue vector with parameters
  ///
  /// @param indValue vector of vectors with variable indices
  ///
  /// @see BaseReaction::createReaction(std::vector<double> &paraValue,...)
  ///
  CreationPrimordiaTime(std::vector<double> &paraValue, 
                   std::vector< std::vector<size_t> > &indValue );
  
  ///
  /// @brief Derivative function for this reaction class
  ///
  /// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
  ///
  void derivs(Tissue &T,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      DataMatrix &cellDerivs,
	      DataMatrix &wallDerivs,
	      DataMatrix &vertexDerivs );
  ///
  /// @brief Update function for this reaction class
  ///
  /// @see BaseReaction::update(Tissue &T,...)
  ///
  void update(Tissue &T,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      double h);

};




///
/// @brief Updates list of vertices with a given force applied where the force is 
/// linearly increased from zero across a given time span (deltaT).
///
/// @details In a model file the reaction is defined as
/// @verbatim
/// VertexFromRotationalForceLinear 1/2/3(dimension+1) 1 (no of vertices)
/// Force component(s) deltaT
/// 1st vertex index
/// 2nd vertex index
/// ...
/// @endverbatim
///
class VertexFromRotationalForceLinear : public BaseReaction {
  
private:
  
  double timeFactor_;
  std::vector<std::vector<double> > boundVerticesUp, boundVerticesDn;

public:
  
  ///
  /// @brief Main constructor
  ///
  /// This is the main constructor which sets the parameters and variable
  /// indices that defines the reaction.
  ///
  /// @param paraValue vector with parameters
  ///
  /// @param indValue vector of vectors with variable indices
  ///
  /// @see BaseReaction::createReaction(std::vector<double> &paraValue,...)
  ///
  VertexFromRotationalForceLinear(std::vector<double> &paraValue, 
                                  std::vector< std::vector<size_t> > &indValue );
  

  void initiate(Tissue &T,
                DataMatrix &cellData,
                DataMatrix &wallData,
                DataMatrix &vertexData,
                DataMatrix &cellDerivs,
                DataMatrix &wallDerivs,
                DataMatrix &vertexDerivs );
  
  ///
  /// @brief Derivative function for this reaction class
  ///
  /// @see BaseReaction::derivs(Compartment &compartment,size_t species,...)
  ///
  void derivs(Tissue &T,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      DataMatrix &cellDerivs,
	      DataMatrix &wallDerivs,
	      DataMatrix &vertexDerivs );

  ///
  /// @brief Update function for this reaction class
  ///
  /// @see BaseReaction::update(Tissue &T,...)
  ///
  void update(Tissue &T,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      double h);
};









#endif //ADHOCREACTION_H
