//
// Filename     : baseDirectionDivision.h
// Description  : The common base for classes describing direction division updates
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : June 2007
// Revision     : $Id:$
//
#ifndef BASEDIRECTIONDIVISION_H
#define BASEDIRECTIONDIVISION_H

#include<vector>
#include<string>
#include<iostream>
#include<fstream>
#include"myTypedefs.h"

class Tissue;

///
/// @brief Base class for rules on how Directions are updated at cell divisions.
///
/// @see Direction
///
class BaseDirectionDivision {
  
 private:
  
  std::string id_;
  std::vector<double> parameter_;           
  std::vector<std::string> parameterId_;           
  std::vector< std::vector<size_t> > variableIndex_;
  
 public:
  
  static BaseDirectionDivision* createDirectionDivision(std::vector<double> &paraValue, 
																										std::vector< std::vector<size_t> > 
																										&indValue,
																										std::string idValue );
  static BaseDirectionDivision* createDirectionDivision( std::istream &IN ); 
  
  //Constructor/destructor not defined!
  //BaseDirectionDivision();
  virtual ~BaseDirectionDivision();
  
  BaseDirectionDivision & operator=( const BaseDirectionDivision & baseDirectionDivisionCopy );
  
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
  
  virtual void update(Tissue &T,size_t cellI,
		      DataMatrix &cellData,
		      DataMatrix &wallData,
		      DataMatrix &vertexData,
		      DataMatrix &cellDerivs,
		      DataMatrix &wallDerivs,
		      DataMatrix &vertexDerivs );
  virtual void print( std::ofstream &os );
};

//!Returns the id string
inline std::string BaseDirectionDivision::id() const 
{
  return id_;
}

//!Returns the number of parameters
inline size_t BaseDirectionDivision::numParameter() const 
{
  return parameter_.size();
}

//!Returns the number of variable-index levels
inline size_t BaseDirectionDivision::numVariableIndexLevel() const 
{
  return variableIndex_.size();
}

//!Returns the number of variable-index within one level 
inline size_t BaseDirectionDivision::numVariableIndex(size_t level) const 
{
  return variableIndex_[level].size();
}

//!Returns a parameter
inline double BaseDirectionDivision::parameter(size_t i) const 
{
  return parameter_[i];
}

//!Returns a reference of a parameter
inline double& BaseDirectionDivision::parameterAddress(size_t i) 
{
  return parameter_[i];
}

//!Returns the id for a parameter
inline std::string BaseDirectionDivision::parameterId(size_t i) const 
{
  return parameterId_[i];
}

//!Returns the index of variable used in the division equation
inline size_t BaseDirectionDivision::variableIndex(size_t i,size_t j) const 
{
  return variableIndex_[i][j];
}

//!Sets the id string
inline void BaseDirectionDivision::setId(std::string value) 
{
  id_=value;
}

//!Sets the parameter value with index i
inline void BaseDirectionDivision::setParameter(size_t i,double value) 
{
  parameter_[i]=value;
}

//!Sets all parameters within a directionDivision
inline void BaseDirectionDivision::setParameter(std::vector<double> &value) 
{
  parameter_ = value;
}

//!Sets the parameter id value with index i
inline void BaseDirectionDivision::setParameterId(size_t i,std::string value) 
{
  parameterId_[i]=value;
}

//!Sets all parameter id's within a directionDivision
inline void BaseDirectionDivision::setParameterId(std::vector<std::string> &value) 
{
  parameterId_=value;
}

//!Sets the variable index at specific level and place
inline void BaseDirectionDivision::
setVariableIndex(size_t i, size_t j,size_t value) 
{
  variableIndex_[i][j] = value;
}

//!Sets the variable indeces for a given level
inline void BaseDirectionDivision::
setVariableIndex(size_t i, std::vector<size_t> &value) 
{
  variableIndex_[i] = value;
}

//!Sets the variable indeces used in a directionDivision
inline void BaseDirectionDivision::
setVariableIndex(std::vector< std::vector<size_t> > &value) 
{
  variableIndex_ = value;
}

#endif


