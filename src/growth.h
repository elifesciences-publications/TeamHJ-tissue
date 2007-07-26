/**
 * Filename     : growth.h
 * Description  : Classes describing growth updates
 * Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
 * Created      : April 2006
 * Revision     : $Id:$
 */
#ifndef GROWTH_H
#define GROWTH_H

#include<cmath>

#include"tissue.h"
#include"baseReaction.h"

//!Exponential wall growth which truncates at threshold length
/*!Exponential growth truncated at threshold length. The wall lengths,
  L, are updated according to

  dL/dt = p_0*L*(1-L/p_1)

  p_0 is the growth rate, p_1 is the threshold for growth truncation.
  In addition, the column index for the wall length should be given.
 */
class WallGrowthExponentialTruncated : public BaseReaction {
  
 public:
  
  WallGrowthExponentialTruncated(std::vector<double> &paraValue, 
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

//!Exponential strech-driven wall growth which truncates at threshold length
/*!Exponential growth driven by a streched wall and truncated at
  threshold length. The wall lengths, L, are updated only if the
  length is shorter than the distance between the vertices of the
  wall and then according to

  dL/dt = p_0*(d_v-L)*L*(1-L/p_1).

  p_0 is the growth rate, p_1 is the threshold for growth truncation.
  d_v is the distance between the two wall vertices.

  In addition, the column index for the wall length should be given.
 */
class WallGrowthExponentialStressTruncated : public BaseReaction {
  
 public:
  
  WallGrowthExponentialStressTruncated(std::vector<double> &paraValue, 
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

//!Constant strech-driven wall growth 
/*!Constant growth driven by a streched wall. The wall lengths, L, are
  updated only if the length is shorter than the distance between the
  vertices of the wall and then according to

  dL/dt = p_0*(d_v-L)

  p_0 is the growth rate. 
  d_v is the distance between the two wall vertices.

  In addition, the column index for the wall length should be given.
 */
class WallGrowthConstantStress : public BaseReaction {
  
 public:
  
  WallGrowthConstantStress(std::vector<double> &paraValue, 
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

//!Constant strech-driven wall growth 
/*!Constant growth driven by a streched wall. The wall lengths, L, are
  updated only if the length is shorter than the distance between the
  vertices of the wall and then according to

  dL/dt = p_0*(d_v-L)*f_e

  p_0 is the growth rate. 
  f_e fraction for epidermal walls
  d_v is the distance between the two wall vertices.

  In addition, the column index for the wall length should be given.
 */
class WallGrowthConstantStressEpidermalAsymmetric : public BaseReaction {
  
 public:
  
  WallGrowthConstantStressEpidermalAsymmetric(std::vector<double> &paraValue, 
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

//!Growth via vertex movement radially outwards
/*!The tissue grows from vertex movement radially outwards. The update can be
	
  dr/dt = p_0 ( p_1=0) or 
	dr/dt = c*r (p_1=1) 

	p_0 is the rate, 
	p_1 is a flag determining function. 
*/
class MoveVertexRadially : public BaseReaction {
  
 public:
  
  MoveVertexRadially(std::vector<double> &paraValue, 
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

class WallLengthGrowExperimental : public BaseReaction
{
public:
	WallLengthGrowExperimental(std::vector<double> &paraValue,
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
/// @brief Updates the water volume variable given osmotic and turgor
/// potentials
///
/// This function uses a constant osmotic potential and calculates the turgor
/// potential to calculate water intake into the cell according to
///
/// @f \frac{V_w}{dt} = p_0 A (p_1-p_2T) @/f
///
/// where V_w is the water volume, T is the turgor,p_0 is the rate, p_1 is the
/// osmotic potential and p_2 is an scaling factor. Also p_3=denyShrink_flag
/// and p_4=allowNegTurgor_flag can be set to restrict the behavior.
///
/// The turgor, T, is calculated as T=V_w-V.
///
class WaterVolumeFromTurgor : public BaseReaction
{
public:
	WaterVolumeFromTurgor(std::vector<double> &paraValue,
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
/// @brief Updates 'concentration' variables according to volume changes from
/// derivatives of the vertex positions.
///
/// The dilution of 'concentration' variables, C_i are calculated according to
///
/// @f$ \frac{C_i}{dt} = - \frac{C_i}{V}\frac{dV}{dt} @f$
///
/// where V is the volume, and dV/dt is calculated from vertex position
/// derivatives.
///
/// @note Since this function uses the derivatives of the vertex positions it
/// needs to be applied after all other derivatives applied to the vertices.
///
class DilutionFromVertexDerivs : public BaseReaction
{
public:
	DilutionFromVertexDerivs(std::vector<double> &paraValue,
													 std::vector< std::vector<size_t> > &indValue);
	
	void derivs(Tissue &T,
							std::vector< std::vector<double> > &cellData,
							std::vector< std::vector<double> > &wallData,
							std::vector< std::vector<double> > &vertexData,
							std::vector< std::vector<double> > &cellDerivs,
							std::vector< std::vector<double> > &wallDerivs,
							std::vector< std::vector<double> > &vertexDerivs);
};

#endif
