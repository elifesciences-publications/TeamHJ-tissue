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


  void derivsWithAbs(Tissue &T,
     DataMatrix &cellData,
     DataMatrix &wallData,
     DataMatrix &vertexData,
     DataMatrix &cellDerivs,
     DataMatrix &wallDerivs,
     DataMatrix &vertexDerivs,
     DataMatrix &sdydtCell,
     DataMatrix &sdydtWall,
     DataMatrix &sdydtVertex );

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







class ThresholdSwitch : public BaseReaction {
  
 public:
  ///
  /// @brief Switches an output variable from 0 to 1 if an input variable is above a threshold. 
  ///
  /// This is the main constructor which checks and sets the parameters and
  /// variable indices that defines the reaction.
  ///
  /// @param paraValue vector with parameters
  ///
  /// @param indValue vector of vectors with variable indices
  ///

/// In the model file, the reaction is specified as:
/// @verbatim
/// ThresholdSwitch 2 2 1 1   # number of parameters is set to two (threshold and switch_type)
/// threshold		 		  # threshold above which a variable is reset to zero.  
/// switch_type		 		  # the switch_type parameter takes the values 0 and 1 for defining the reversible and irreversible switch, respectively.  
/// index_var   	 		  # index of the index variable upstream the switch.
/// index_var_out  			  # list of updated indices - for the moment it can be just one index - where the output of the switch is written. 
/// @endverbatim

/// @note This function makes a downstream species reversibly or irreversibly  
/// switch from 0 to 1, upon being above a certain threshold of an upstream variable.

  /// @see BaseReaction::createReaction(std::vector<double> &paraValue,...)
  ThresholdSwitch(std::vector<double> &paraValue, 
			     std::vector< std::vector<size_t> > &indValue );
		
  ///
  /// @brief This class does not use derivatives for updates.
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

  void derivsWithAbs(Tissue &T,
     DataMatrix &cellData,
     DataMatrix &wallData,
     DataMatrix &vertexData,
     DataMatrix &cellDerivs,
     DataMatrix &wallDerivs,
     DataMatrix &vertexDerivs,
     DataMatrix &sdydtCell,
     DataMatrix &sdydtWall,
     DataMatrix &sdydtVertex );
  ///
  /// @brief Update function for this reaction class
  ///
  /// @see BaseReaction::update(double h, double t, ...)
  ///
	void update(Tissue &T,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      double h);
};


class AndGate : public BaseReaction {
  
 public:
  ///
  /// @brief Main constructor
  ///
  /// This is the main constructor which checks and sets the parameters and
  /// variable indices that defines the reaction.
  ///
  /// @param paraValue vector with parameters
  ///
  /// @param indValue vector of vectors with variable indices
  ///

/// In the model file, the reaction is specified as:
/// @verbatim
/// AndGate 1 2 2 1   	 	 # number of parameters is set to one (gate_type)
/// gate_type		 		  # the gate_type parameter takes the values 0 and 1 for defining the reversible and irreversible gate, respectively.  
/// index_var1   	 		  # index of the fist variable upstream the gate.
/// index_var2   	 		  # index of the second variable upstream the gate.
/// index_var_out  			  # updated index where the output of the gate is written. 
/// @endverbatim

/// @note This logical gate function makes a downstream species reversibly or irreversibly  
/// switch from 0 to 1 if the two input variables are 1.

  /// @see BaseReaction::createReaction(std::vector<double> &paraValue,...)
  AndGate(std::vector<double> &paraValue, 
			     std::vector< std::vector<size_t> > &indValue );
		
  ///
  /// @brief This class does not use derivatives for updates.
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

  void derivsWithAbs(Tissue &T,
     DataMatrix &cellData,
     DataMatrix &wallData,
     DataMatrix &vertexData,
     DataMatrix &cellDerivs,
     DataMatrix &wallDerivs,
     DataMatrix &vertexDerivs,
     DataMatrix &sdydtCell,
     DataMatrix &sdydtWall,
     DataMatrix &sdydtVertex );
  ///
  /// @brief Update function for this reaction class
  ///
  /// @see BaseReaction::update(double h, double t, ...)
  ///
	void update(Tissue &T,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      double h);
};



class AndNotGate : public BaseReaction {
  
 public:
  ///
  /// @brief Main constructor
  ///
  /// This is the main constructor which checks and sets the parameters and
  /// variable indices that defines the reaction.
  ///
  /// @param paraValue vector with parameters
  ///
  /// @param indValue vector of vectors with variable indices
  ///

/// In the model file, the reaction is specified as:
/// @verbatim
/// AndNotGate 1 2 2 1   	  # number of parameters is set to one (gate_type)
/// gate_type		 		  # the gate_type parameter takes the values 0 and 1 for defining the reversible and irreversible gate, respectively.  
/// index_var1   	 		  # index of the fist variable upstream the gate.
/// index_var2   	 		  # index of the second variable upstream the gate.
/// index_var_out  			  # updated index where the output of the gate is written. 
/// @endverbatim

/// @note This logical gate function makes a downstream species reversibly or irreversibly  
/// switch from 0 to 1 if a first input variable is 1 and a second input variable is 0.

  /// @see BaseReaction::createReaction(std::vector<double> &paraValue,...)
  AndNotGate(std::vector<double> &paraValue, 
			     std::vector< std::vector<size_t> > &indValue );
		
  ///
  /// @brief This class does not use derivatives for updates.
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

  void derivsWithAbs(Tissue &T,
     DataMatrix &cellData,
     DataMatrix &wallData,
     DataMatrix &vertexData,
     DataMatrix &cellDerivs,
     DataMatrix &wallDerivs,
     DataMatrix &vertexDerivs,
     DataMatrix &sdydtCell,
     DataMatrix &sdydtWall,
     DataMatrix &sdydtVertex );
  ///
  /// @brief Update function for this reaction class
  ///
  /// @see BaseReaction::update(double h, double t, ...)
  ///
	void update(Tissue &T,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      double h);
};



class AndSpecialGate : public BaseReaction {
  
 public:
  ///
  /// @brief Main constructor
  ///
  /// This is the main constructor which checks and sets the parameters and
  /// variable indices that defines the reaction.
  ///
  /// @param paraValue vector with parameters
  ///
  /// @param indValue vector of vectors with variable indices
  ///

/// In the model file, the reaction is specified as:
/// @verbatim
/// AndSpecialGate 0 2 3 1    # number of parameters is set to zero
/// index_var1   	 		  # index of the fist variable upstream the gate.
/// index_var2   	 		  # index of the second variable upstream the gate.
/// index_var3   	 		  # index of the third variable upstream the gate.
/// index_var_out  			  # updated index where the output of the gate is written. 
/// @endverbatim

/// @note This logical gate function makes a downstream species   
/// switch from 0 to 1 if the folowing conditions are met: 
 	// the first input variable is 1 
 	// the second input variable is 0 
 	// the third input variable is 0.

  /// @see BaseReaction::createReaction(std::vector<double> &paraValue,...)
  AndSpecialGate(std::vector<double> &paraValue, 
			     std::vector< std::vector<size_t> > &indValue );
		
  ///
  /// @brief This class does not use derivatives for updates.
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

  void derivsWithAbs(Tissue &T,
     DataMatrix &cellData,
     DataMatrix &wallData,
     DataMatrix &vertexData,
     DataMatrix &cellDerivs,
     DataMatrix &wallDerivs,
     DataMatrix &vertexDerivs,
     DataMatrix &sdydtCell,
     DataMatrix &sdydtWall,
     DataMatrix &sdydtVertex );
  ///
  /// @brief Update function for this reaction class
  ///
  /// @see BaseReaction::update(double h, double t, ...)
  ///
	void update(Tissue &T,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      double h);
};


class AndGateCount : public BaseReaction {
  
 public:
  ///
  /// @brief Main constructor
  ///
  /// This is the main constructor which checks and sets the parameters and
  /// variable indices that defines the reaction.
  ///
  /// @param paraValue vector with parameters
  ///
  /// @param indValue vector of vectors with variable indices
  ///

/// In the model file, the reaction is specified as:
/// @verbatim
/// AndGateCount 0 2 2 1   	  # number of parameters is set to zero
/// index_var1   	 		  # index of the fist variable upstream the gate.
/// index_var2   	 		  # index of the second variable upstream the gate.
/// index_var_out  			  # updated index where the output of the gate is written. 
/// @endverbatim

/// @note This logical gate function makes a downstream species add +1  
/// if the two input variables are 1.

  /// @see BaseReaction::createReaction(std::vector<double> &paraValue,...)
  AndGateCount(std::vector<double> &paraValue, 
			     std::vector< std::vector<size_t> > &indValue );
		
  ///
  /// @brief This class does not use derivatives for updates.
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

  void derivsWithAbs(Tissue &T,
     DataMatrix &cellData,
     DataMatrix &wallData,
     DataMatrix &vertexData,
     DataMatrix &cellDerivs,
     DataMatrix &wallDerivs,
     DataMatrix &vertexDerivs,
     DataMatrix &sdydtCell,
     DataMatrix &sdydtWall,
     DataMatrix &sdydtVertex );
  ///
  /// @brief Update function for this reaction class
  ///
  /// @see BaseReaction::update(double h, double t, ...)
  ///
	void update(Tissue &T,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      double h);
};



class OrGateCount : public BaseReaction {
  
 public:
  ///
  /// @brief Main constructor
  ///
  /// This is the main constructor which checks and sets the parameters and
  /// variable indices that defines the reaction.
  ///
  /// @param paraValue vector with parameters
  ///
  /// @param indValue vector of vectors with variable indices
  ///

/// In the model file, the reaction is specified as:
/// @verbatim
/// AndGateCount 0 2 2 1   	  # number of parameters is set to one (switch_type)
/// index_var1   	 		  # index of the fist variable upstream the gate.
/// index_var2   	 		  # index of the second variable upstream the gate.
/// index_var_out  			  # updated index where the output of the gate is written. 
/// @endverbatim

/// @note This logical gate function makes a downstream species add +1  
/// if one of the two input variables is 1.

  /// @see BaseReaction::createReaction(std::vector<double> &paraValue,...)
  OrGateCount(std::vector<double> &paraValue, 
			     std::vector< std::vector<size_t> > &indValue );
		
  ///
  /// @brief This class does not use derivatives for updates.
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

  void derivsWithAbs(Tissue &T,
     DataMatrix &cellData,
     DataMatrix &wallData,
     DataMatrix &vertexData,
     DataMatrix &cellDerivs,
     DataMatrix &wallDerivs,
     DataMatrix &vertexDerivs,
     DataMatrix &sdydtCell,
     DataMatrix &sdydtWall,
     DataMatrix &sdydtVertex );
  ///
  /// @brief Update function for this reaction class
  ///
  /// @see BaseReaction::update(double h, double t, ...)
  ///
	void update(Tissue &T,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      double h);
};


class OrSpecialGateCount : public BaseReaction {
  
 public:
  ///
  /// @brief Main constructor
  ///
  /// This is the main constructor which checks and sets the parameters and
  /// variable indices that defines the reaction.
  ///
  /// @param paraValue vector with parameters
  ///
  /// @param indValue vector of vectors with variable indices
  ///

/// In the model file, the reaction is specified as:
/// @verbatim
/// AndGateCount 0 2 2 1   	  # number of parameters is set to one (switch_type)
/// index_var1   	 		  # index of the fist variable upstream the gate.
/// index_var2   	 		  # index of the second variable upstream the gate.
/// index_var_out  			  # updated index where the output of the gate is written. 
/// @endverbatim

/// @note This logical gate function makes a downstream species add +1  
/// if the first input variables is 1 or if the second input variable is larger than 0. 

  /// @see BaseReaction::createReaction(std::vector<double> &paraValue,...)
  OrSpecialGateCount(std::vector<double> &paraValue, 
			     std::vector< std::vector<size_t> > &indValue );
		
  ///
  /// @brief This class does not use derivatives for updates.
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

  void derivsWithAbs(Tissue &T,
     DataMatrix &cellData,
     DataMatrix &wallData,
     DataMatrix &vertexData,
     DataMatrix &cellDerivs,
     DataMatrix &wallDerivs,
     DataMatrix &vertexDerivs,
     DataMatrix &sdydtCell,
     DataMatrix &sdydtWall,
     DataMatrix &sdydtVertex );
  ///
  /// @brief Update function for this reaction class
  ///
  /// @see BaseReaction::update(double h, double t, ...)
  ///
	void update(Tissue &T,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      double h);
};


class Count : public BaseReaction {
  
 public:
  ///
  /// @brief Main constructor
  ///
  /// This is the main constructor which checks and sets the parameters and
  /// variable indices that defines the reaction.
  ///
  /// @param paraValue vector with parameters
  ///
  /// @param indValue vector of vectors with variable indices
  ///

/// In the model file, the reaction is specified as:
/// @verbatim
/// Count 0 1 1 	  	      # number of parameters is set to zero
/// index_var_out  			  # updated index where the output of the gate is written. 
/// @endverbatim

/// @note This function makes a downstream species add +1. 
/// 

  /// @see BaseReaction::createReaction(std::vector<double> &paraValue,...)
  Count(std::vector<double> &paraValue, 
			     std::vector< std::vector<size_t> > &indValue );
		
  ///
  /// @brief This class does not use derivatives for updates.
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

  void derivsWithAbs(Tissue &T,
     DataMatrix &cellData,
     DataMatrix &wallData,
     DataMatrix &vertexData,
     DataMatrix &cellDerivs,
     DataMatrix &wallDerivs,
     DataMatrix &vertexDerivs,
     DataMatrix &sdydtCell,
     DataMatrix &sdydtWall,
     DataMatrix &sdydtVertex );
  ///
  /// @brief Update function for this reaction class
  ///
  /// @see BaseReaction::update(double h, double t, ...)
  ///
	void update(Tissue &T,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      double h);
};


class FlagCount : public BaseReaction {
  
 public:
  ///
  /// @brief Main constructor
  ///
  /// This is the main constructor which checks and sets the parameters and
  /// variable indices that defines the reaction.
  ///
  /// @param paraValue vector with parameters
  ///
  /// @param indValue vector of vectors with variable indices
  ///

/// In the model file, the reaction is specified as:
/// @verbatim
/// FlagCount 0 2 1 1	  	  # number of parameters is set to zero
/// index_var   	 		  # index of the variable upstream the gate.
/// index_var_out  			  # updated index where the output of the gate is written. 
/// @endverbatim

/// @note This logical gate function makes a downstream species add +1  
/// if the input variable is 1.

  /// @see BaseReaction::createReaction(std::vector<double> &paraValue,...)
  FlagCount(std::vector<double> &paraValue, 
			     std::vector< std::vector<size_t> > &indValue );
		
  ///
  /// @brief This class does not use derivatives for updates.
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

  void derivsWithAbs(Tissue &T,
     DataMatrix &cellData,
     DataMatrix &wallData,
     DataMatrix &vertexData,
     DataMatrix &cellDerivs,
     DataMatrix &wallDerivs,
     DataMatrix &vertexDerivs,
     DataMatrix &sdydtCell,
     DataMatrix &sdydtWall,
     DataMatrix &sdydtVertex );
  ///
  /// @brief Update function for this reaction class
  ///
  /// @see BaseReaction::update(double h, double t, ...)
  ///
	void update(Tissue &T,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      double h);
};

class ThresholdReset : public BaseReaction {
  
 public:
  ///
  /// @brief Main constructor
  ///
  /// This is the main constructor which checks and sets the parameters and
  /// variable indices that defines the reaction.
  ///
  /// @param paraValue vector with parameters
  ///
  /// @param indValue vector of vectors with variable indices
  ///

/// In the model file, the reaction is specified as:
/// @verbatim
/// ThresholdReset 2 2 1 1   # number of parameters is set to two (threshold and switch_type)
/// threshold		 		  # threshold above which a variable is reset to zero.  
/// switch_type		 		  # the switchtype parameter takes the values 0 and 1 for defining the reversible and irreversible switch, respectively.  
/// index_var   	 		  # index of the index variable upstream the switch.
/// index_var_out  			  # list of updated indices - for the moment it can be just one index - where the output of the switch is written. 
/// @endverbatim

/// @note This function makes a downstream species reversibly or irreversibly  
/// switch from 1 to 0, upon being above a certain threshold of an upstream variable.

  /// @see BaseReaction::createReaction(std::vector<double> &paraValue,...)
  ThresholdReset(std::vector<double> &paraValue, 
			     std::vector< std::vector<size_t> > &indValue );
		
  ///
  /// @brief This class does not use derivatives for updates.
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
  
  void derivsWithAbs(Tissue &T,
     DataMatrix &cellData,
     DataMatrix &wallData,
     DataMatrix &vertexData,
     DataMatrix &cellDerivs,
     DataMatrix &wallDerivs,
     DataMatrix &vertexDerivs,
     DataMatrix &sdydtCell,
     DataMatrix &sdydtWall,
     DataMatrix &sdydtVertex );
  ///
  /// @brief Update function for this reaction class
  ///
  /// @see BaseReaction::update(double h, double t, ...)
  ///
	void update(Tissue &T,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      double h);
};


#endif //ADHOCREACTION_H
