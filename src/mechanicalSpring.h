/**
 * Filename     : mechanicalSpring.h
 * Description  : Classes describing mechanical updates for wall springs
 * Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
 * Created      : September 2007
 * Revision     : $Id:$
 */
#ifndef MECHANICALSPRING_H
#define MECHANICALSPRING_H

#include"tissue.h"
#include"baseReaction.h"
#include<cmath>

//!Updates vertices from an asymmetric wall spring potential
class VertexFromWallSpring : public BaseReaction {
  
 public:
  
  VertexFromWallSpring(std::vector<double> &paraValue, 
		       std::vector< std::vector<size_t> > 
		       &indValue );
  
  void derivs(Tissue &T,
	      std::vector< std::vector<double> > &cellData,
	      std::vector< std::vector<double> > &wallData,
	      std::vector< std::vector<double> > &vertexData,
	      std::vector< std::vector<double> > &cellDerivs,
	      std::vector< std::vector<double> > &wallDerivs,
	      std::vector< std::vector<double> > &vertexDerivs );
};

///
/// @brief Updates vertices from an asymmetric wall spring potential
/// where the spring constants are taken from two wall variables
///
/// This class updates according to a spring potential where two
/// spring constants are set from elsewhere, i.e. values are taken
/// from two variables representing the values for the TWO parts of
/// the wall connected to the two neighboring cells. 
///
class VertexFromDoubleWallSpring : public BaseReaction {
  
 public:
  
  VertexFromDoubleWallSpring(std::vector<double> &paraValue, 
			     std::vector< std::vector<size_t> > 
			     &indValue );
  
  void derivs(Tissue &T,
	      std::vector< std::vector<double> > &cellData,
	      std::vector< std::vector<double> > &wallData,
	      std::vector< std::vector<double> > &vertexData,
	      std::vector< std::vector<double> > &cellDerivs,
	      std::vector< std::vector<double> > &wallDerivs,
	      std::vector< std::vector<double> > &vertexDerivs );
};

///
/// @brief Updates vertices from an asymmetric wall spring potential with a spatial factor
///
class VertexFromWallSpringSpatial : public BaseReaction {
  
 private:
  
  double Kpow_;
  
 public:
  
  VertexFromWallSpringSpatial(std::vector<double> &paraValue, 
			      std::vector< std::vector<size_t> > 
			      &indValue );
  
  void derivs(Tissue &T,
	      std::vector< std::vector<double> > &cellData,
	      std::vector< std::vector<double> > &wallData,
	      std::vector< std::vector<double> > &vertexData,
	      std::vector< std::vector<double> > &cellDerivs,
	      std::vector< std::vector<double> > &wallDerivs,
	      std::vector< std::vector<double> > &vertexDerivs );
};

///
/// @brief Updates vertices from an asymmetric wall spring potential with a spatial MT factor
///
/// This function has a base level for the spring constant and then adds a
/// factor dependent on the MT direction in the cells as well as a spatial
/// factor. The spring constant of the wall is determined by:
///
/// @f[ k = p_0 + (p_1 + p_2*spatialFactor)*mtFactor @f]
///
/// where the spatial factor is determined by a Hill function and the MT
/// factor comes from the absolute value of the scalar products of the wall
/// vector and cell direction vectors.
///
class VertexFromWallSpringMTSpatial : public BaseReaction {
  
 private:
  
  double Kpow_;
  
 public:
  
  VertexFromWallSpringMTSpatial(std::vector<double> &paraValue, 
				std::vector< std::vector<size_t> > 
				&indValue );
  
  void derivs(Tissue &T,
	      std::vector< std::vector<double> > &cellData,
	      std::vector< std::vector<double> > &wallData,
	      std::vector< std::vector<double> > &vertexData,
	      std::vector< std::vector<double> > &cellDerivs,
	      std::vector< std::vector<double> > &wallDerivs,
	      std::vector< std::vector<double> > &vertexDerivs );
};

///
/// @brief Updates vertices from an spatially asymmetric wall spring potential
/// given by microtubule directions
///
class VertexFromWallSpringMT : public BaseReaction {
  
 public:
  
  VertexFromWallSpringMT(std::vector<double> &paraValue, 
			 std::vector< std::vector<size_t> > 
			 &indValue );
  
  void derivs(Tissue &T,
	      std::vector< std::vector<double> > &cellData,
	      std::vector< std::vector<double> > &wallData,
	      std::vector< std::vector<double> > &vertexData,
	      std::vector< std::vector<double> > &cellDerivs,
	      std::vector< std::vector<double> > &wallDerivs,
	      std::vector< std::vector<double> > &vertexDerivs );
  void initiate(Tissue &T,
		std::vector< std::vector<double> > &cellData,
		std::vector< std::vector<double> > &wallData,
		std::vector< std::vector<double> > &vertexData );
};

///
/// @brief Updates vertices from an spatially asymmetric wall spring potential
/// given by microtubule directions and also updates the spring constants
/// slowly.
///
class VertexFromWallSpringMTHistory : public BaseReaction {
  
 public:
  
  VertexFromWallSpringMTHistory(std::vector<double> &paraValue, 
				std::vector< std::vector<size_t> > 
				&indValue );
  
  void derivs(Tissue &T,
	      std::vector< std::vector<double> > &cellData,
	      std::vector< std::vector<double> > &wallData,
	      std::vector< std::vector<double> > &vertexData,
	      std::vector< std::vector<double> > &cellDerivs,
	      std::vector< std::vector<double> > &wallDerivs,
	      std::vector< std::vector<double> > &vertexDerivs );
  void initiate(Tissue &T,
		std::vector< std::vector<double> > &cellData,
		std::vector< std::vector<double> > &wallData,
		std::vector< std::vector<double> > &vertexData);
};

//!Updates vertices from an asymmetric epidermal wall spring potential
class VertexFromEpidermalWallSpring : public BaseReaction {
  
 public:
  
  VertexFromEpidermalWallSpring(std::vector<double> &paraValue, 
				std::vector< std::vector<size_t> > 
				&indValue );
  
  void derivs(Tissue &T,
	      std::vector< std::vector<double> > &cellData,
	      std::vector< std::vector<double> > &wallData,
	      std::vector< std::vector<double> > &vertexData,
	      std::vector< std::vector<double> > &cellDerivs,
	      std::vector< std::vector<double> > &wallDerivs,
	      std::vector< std::vector<double> > &vertexDerivs );
};

//!Updates vertices from an asymmetric epidermal wall spring potential
class VertexFromEpidermalCellWallSpring : public BaseReaction {
  
 public:
  
  VertexFromEpidermalCellWallSpring(std::vector<double> &paraValue, 
				    std::vector< std::vector<size_t> > 
				    &indValue );
  
  void derivs(Tissue &T,
	      std::vector< std::vector<double> > &cellData,
	      std::vector< std::vector<double> > &wallData,
	      std::vector< std::vector<double> > &vertexData,
	      std::vector< std::vector<double> > &cellDerivs,
	      std::vector< std::vector<double> > &wallDerivs,
	      std::vector< std::vector<double> > &vertexDerivs );
};

class VertexFromWallSpringExperimental : public BaseReaction
{
 public:
  VertexFromWallSpringExperimental(std::vector<double> &paraValue,
				   std::vector< std::vector<size_t> > &indValue);
  
  void derivs(Tissue &T,
	      std::vector< std::vector<double> > &cellData,
	      std::vector< std::vector<double> > &wallData,
	      std::vector< std::vector<double> > &vertexData,
	      std::vector< std::vector<double> > &cellDerivs,
	      std::vector< std::vector<double> > &wallDerivs,
	      std::vector< std::vector<double> > &vertexDerivs);
};

///
/// @brief Updates vertices from an asymmetric wall spring potential depending
/// on cell concentrations
///
class VertexFromWallSpringConcentrationHill : public BaseReaction {
  
 public:
  
  VertexFromWallSpringConcentrationHill(std::vector<double> &paraValue, 
					std::vector< std::vector<size_t> > 
					&indValue );
  
  void derivs(Tissue &T,
	      std::vector< std::vector<double> > &cellData,
	      std::vector< std::vector<double> > &wallData,
	      std::vector< std::vector<double> > &vertexData,
	      std::vector< std::vector<double> > &cellDerivs,
	      std::vector< std::vector<double> > &wallDerivs,
	      std::vector< std::vector<double> > &vertexDerivs );
};

class VertexFromWallSpringMTConcentrationHill : public BaseReaction {
  
 public:
  
  VertexFromWallSpringMTConcentrationHill(std::vector<double> &paraValue, 
					  std::vector< std::vector<size_t> > 
					  &indValue );
  
  void derivs(Tissue &T,
	      std::vector< std::vector<double> > &cellData,
	      std::vector< std::vector<double> > &wallData,
	      std::vector< std::vector<double> > &vertexData,
	      std::vector< std::vector<double> > &cellDerivs,
	      std::vector< std::vector<double> > &wallDerivs,
	      std::vector< std::vector<double> > &vertexDerivs );
};

#endif
