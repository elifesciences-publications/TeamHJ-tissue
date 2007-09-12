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
class VertexFromWallSpringAsymmetric : public BaseReaction {
  
 public:
  
  VertexFromWallSpringAsymmetric(std::vector<double> &paraValue, 
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
class VertexFromEpidermalWallSpringAsymmetric : public BaseReaction {
  
 public:
  
  VertexFromEpidermalWallSpringAsymmetric(std::vector<double> &paraValue, 
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
class VertexFromEpidermalCellWallSpringAsymmetric : public BaseReaction {
  
 public:
  
  VertexFromEpidermalCellWallSpringAsymmetric(std::vector<double> &paraValue, 
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

#endif
