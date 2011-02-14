//
// Filename     : baseDirectionUpdate.h
// Description  : The common base for classes describing direction updates
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : June 2007
// Revision     : $Id:$
//
#ifndef BASEDIRECTIONUPDATE_H
#define BASEDIRECTIONUPDATE_H

#include<vector>
#include<string>
#include<iostream>
#include<fstream>
//#include"tissue.h"

class Tissue;

///
/// @brief Base class for rules on how defined Directions are updated during simulations.
///
/// @see Direction
///
class BaseDirectionUpdate {
  
 private:
  
  std::string id_;
  std::vector<double> parameter_;           
  std::vector<std::string> parameterId_;           
  std::vector< std::vector<size_t> > variableIndex_;
  
 public:
  
  static BaseDirectionUpdate* createDirectionUpdate(std::vector<double> &paraValue, 
						    std::vector< std::vector<size_t> > 
						    &indValue,
						    std::string idValue );
  static BaseDirectionUpdate* createDirectionUpdate( std::istream &IN ); 
  
  //Constructor/destructor not defined!
  //BaseDirectionUpdate();
  virtual ~BaseDirectionUpdate();
  
  BaseDirectionUpdate & operator=( const BaseDirectionUpdate & baseDirectionUpdateCopy );
  
  // Get values
  inline std::string id() const;
  inline size_t numParameter() const;  
  inline size_t numVariableIndexLevel() const;
  inline size_t numVariableIndex(size_t level) const;
  
  inline double parameter(size_t i) const;
  inline double& parameterAddress(size_t i);
  inline std::string parameterId(size_t i) const;
  inline size_t variableIndex(size_t i,size_t j) const;
  
  // Set values
  inline void setId(std::string value);
  inline void setParameter(size_t i,double value);
  inline void setParameter(std::vector<double> &value);
  inline void setParameterId(size_t i,std::string value);
  inline void setParameterId(std::vector<std::string> &value);
  inline void setVariableIndex(size_t i, size_t j,size_t value);
  inline void setVariableIndex(size_t i, std::vector<size_t> &value);
  inline void setVariableIndex(std::vector< std::vector<size_t> > &value);
  
  virtual void initiate(Tissue &T,
			std::vector< std::vector<double> > &cellData,
			std::vector< std::vector<double> > &wallData,
			std::vector< std::vector<double> > &vertexData,
			std::vector< std::vector<double> > &cellDerivs,
			std::vector< std::vector<double> > &wallDerivs,
			std::vector< std::vector<double> > &vertexDerivs );
  virtual void update(Tissue &T, double h,
		      std::vector< std::vector<double> > &cellData,
		      std::vector< std::vector<double> > &wallData,
		      std::vector< std::vector<double> > &vertexData,
		      std::vector< std::vector<double> > &cellDerivs,
		      std::vector< std::vector<double> > &wallDerivs,
		      std::vector< std::vector<double> > &vertexDerivs );
  virtual void print( std::ofstream &os );
};

//!Returns the id string
inline std::string BaseDirectionUpdate::id() const 
{
  return id_;
}

//!Returns the number of parameters
inline size_t BaseDirectionUpdate::numParameter() const 
{
  return parameter_.size();
}

//!Returns the number of variable-index levels
inline size_t BaseDirectionUpdate::numVariableIndexLevel() const 
{
  return variableIndex_.size();
}

//!Returns the number of variable-index within one level 
inline size_t BaseDirectionUpdate::numVariableIndex(size_t level) const 
{
  return variableIndex_[level].size();
}

//!Returns a parameter
inline double BaseDirectionUpdate::parameter(size_t i) const 
{
  return parameter_[i];
}

//!Returns a reference of a parameter
inline double& BaseDirectionUpdate::parameterAddress(size_t i) 
{
  return parameter_[i];
}

//!Returns the id for a parameter
inline std::string BaseDirectionUpdate::parameterId(size_t i) const 
{
  return parameterId_[i];
}

//!Returns the index of variable used in the update equation
inline size_t BaseDirectionUpdate::variableIndex(size_t i,size_t j) const 
{
  return variableIndex_[i][j];
}

//!Sets the id string
inline void BaseDirectionUpdate::setId(std::string value) 
{
  id_=value;
}

//!Sets the parameter value with index i
inline void BaseDirectionUpdate::setParameter(size_t i,double value) 
{
  parameter_[i]=value;
}

//!Sets all parameters within a directionUpdate
inline void BaseDirectionUpdate::setParameter(std::vector<double> &value) 
{
  parameter_ = value;
}

//!Sets the parameter id value with index i
inline void BaseDirectionUpdate::setParameterId(size_t i,std::string value) 
{
  parameterId_[i]=value;
}

//!Sets all parameter id's within a directionUpdate
inline void BaseDirectionUpdate::setParameterId(std::vector<std::string> &value) 
{
  parameterId_=value;
}

//!Sets the variable index at specific level and place
inline void BaseDirectionUpdate::
setVariableIndex(size_t i, size_t j,size_t value) 
{
  variableIndex_[i][j] = value;
}

//!Sets the variable indeces for a given level
inline void BaseDirectionUpdate::
setVariableIndex(size_t i, std::vector<size_t> &value) 
{
  variableIndex_[i] = value;
}

//!Sets the variable indeces used in a directionUpdate
inline void BaseDirectionUpdate::
setVariableIndex(std::vector< std::vector<size_t> > &value) 
{
  variableIndex_ = value;
}

#endif


