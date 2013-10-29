//
// Filename     : directionReaction.h
// Description  : Classes describing some reaction updates related to directions
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : May 2008
// Revision     : $Id:$
//
#ifndef DIRECTIONREACTION_H
#define DIRECTIONREACTION_H

#include"tissue.h"
#include"baseReaction.h"
#include<cmath>

///
/// @brief Updates a direction continuosly towards a direction updated for the cell
///
/// @details This function applies the update within the derivative and moves the
/// direction as a function of the differnce in angle.
///
/// In a model file the reaction is defined as
/// @verbatim
/// ContinousMTDirection 1 2 1 1
/// k_rate
/// target index
/// MT index
/// @endverbatim
///
/// @note It is only implemented for two dimensions.
///
class ContinousMTDirection : public BaseReaction
{
public:
  ContinousMTDirection(std::vector<double> &paraValue,
		       std::vector< std::vector<size_t> > &indValue);
  
  void derivs(Tissue &T,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      DataMatrix &cellDerivs,
	      DataMatrix &wallDerivs,
	      DataMatrix &vertexDerivs);
};


///
/// @brief Updates a direction continuosly towards a direction updated for the cell
///
/// @details This function applies the update within the derivative and moves the
/// direction based on velocity vector proportional to the subtraction of 
/// target vector and the direction which is going to be updated.
/// Note: it is only implemented for three dimensions.
///
/// In a model file the reaction is defined as
/// @verbatim
/// ContinousMTDirection3d 1 2 1 1
///
/// k_rate
///
/// target index
///
/// MT index
///
/// @endverbatim
///
class ContinousMTDirection3d : public BaseReaction
{
public:
  ContinousMTDirection3d(std::vector<double> &paraValue,
		       std::vector< std::vector<size_t> > &indValue);
  
  void derivs(Tissue &T,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      DataMatrix &cellDerivs,
	      DataMatrix &wallDerivs,
	      DataMatrix &vertexDerivs);
};

///
/// @brief Updates a direction continuosly towards a direction updated for the cell
///
/// @details This function applies the update within the update and moves the
/// direction toward a given direction using an Euler step directly on the direction vector.
/// It normalizes the vector after the update.
///
/// In a model file the reaction is defined as
/// @verbatim
/// UpdateMTDirection 1 2 1 1
/// k_rate
/// target index
/// MT index
/// @endverbatim
///
class UpdateMTDirection : public BaseReaction
{
 public:
  UpdateMTDirection(std::vector<double> &paraValue,
		    std::vector< std::vector<size_t> > &indValue);
  
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
/// @brief Updates a direction continuosly towards a direction updated for the cell
///
/// @details This function applies the update within the "update" and moves the
/// direction toward a given direction using an Euler step directly on the direction vector.
/// The update is done only in the cells that are close "enough" to equilibrium. Equilibrium 
/// is considered as a state in which the addition of absolute values of velocity vectors of 
/// vertices of a cell is less than a user defined threshold. 
/// It normalizes the vector after the update.
///
/// In a model file the reaction is defined as
/// @verbatim
/// UpdateMTDirectionEquilibrium 3 4 1 1 1 2 
/// k_rate
/// equilibrium velocity threshold
/// stress difference threshold
/// target index
/// MT index
/// max-stress index
/// MT-stress index
/// velocity store index
/// @endverbatim
///
class UpdateMTDirectionEquilibrium : public BaseReaction
{
 public:
  UpdateMTDirectionEquilibrium(std::vector<double> &paraValue,
		    std::vector< std::vector<size_t> > &indValue);
  
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
/// @brief Updates a direction continuosly towards a direction updated for the cell
///
/// @details This function applies the update within the update and moves the
/// direction toward a given direction using an Euler step directly on the direction vector.
/// the update rate depends on both a user defined rate and a Hill function of a concentration 
/// that can be stress or strain anisotropy. It normalizes the vector after the update.
///
/// In a model file the reaction is defined as
/// @verbatim
/// UpdateMTDirectionConcenHill 3 3 1 1 1
/// k_rate
/// k_Hill
/// n_hill
/// target index
/// MT index
/// concentration(anisotropy) index
/// @endverbatim
///
class UpdateMTDirectionConcenHill : public BaseReaction
{
 public:
  UpdateMTDirectionConcenHill(std::vector<double> &paraValue,
		    std::vector< std::vector<size_t> > &indValue);
  
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




class RotatingDirection : public BaseReaction
{
 public:
  RotatingDirection(std::vector<double> &paraValue,
		    std::vector< std::vector<size_t> > &indValue);
  
  void derivs(Tissue &T,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      DataMatrix &cellDerivs,
	      DataMatrix &wallDerivs,
	      DataMatrix &vertexDerivs);
};

#endif //DIRECTIONREACTION_H
