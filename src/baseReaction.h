//
// Filename     : baseReaction.h
// Description  : The common base for classes describing reaction updates
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : April 2006
// Revision     : $Id:$
//
#ifndef BASEREACTION_H
#define BASEREACTION_H

#include<vector>
#include<string>
#include<iostream>
#include<fstream>
//#include"tissue.h"

class Tissue;

//!A base class for describing diff equation updates for dynamical variables
/*! The BaseReaction class is a base class used when defining
  different types of "reaction" classes. Each reaction class uses a
  vector of parameters and variable indeces to calculate a derivative
  of model variables. The variable indeces are divided into multiple
  layers to allow for different types of contributions for different
  variables. The baseReaction class can be seen as a 'factory'
  creating different reactions of different types.
*/ 
class BaseReaction {
  
 private:
  
  std::string id_;
  std::vector<double> parameter_;           
  std::vector<std::string> parameterId_;           
  std::vector< std::vector<size_t> > variableIndex_;
  
 public:
  
  static BaseReaction* createReaction(std::vector<double> &paraValue, 
				      std::vector< std::vector<size_t> > 
				      &indValue,
				      std::string idValue );
  static BaseReaction* createReaction( std::istream &IN ); 
  
  //Constructor/destructor not defined!
  //BaseReaction();
  virtual ~BaseReaction();
  
  BaseReaction & operator=( const BaseReaction & baseReactionCopy );
  
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
												std::vector< std::vector<double> > &vertexData);
  virtual void derivs(Tissue &T,
											std::vector< std::vector<double> > &cellData,
											std::vector< std::vector<double> > &wallData,
											std::vector< std::vector<double> > &vertexData,
											std::vector< std::vector<double> > &cellDerivs,
											std::vector< std::vector<double> > &wallDerivs,
											std::vector< std::vector<double> > &vertexDerivs );
  virtual void update(Tissue &T,
											std::vector< std::vector<double> > &cellData,
											std::vector< std::vector<double> > &wallData,
											std::vector< std::vector<double> > &vertexData,
											double h);
  virtual void print( std::ofstream &os );
};

//!Returns the id string
inline std::string BaseReaction::id() const {
  return id_;
}

//!Returns the number of parameters
inline size_t BaseReaction::numParameter() const {
  return parameter_.size();
}

//!Returns the number of variable-index levels
inline size_t BaseReaction::numVariableIndexLevel() const {
  return variableIndex_.size();
}

//!Returns the number of variable-index within one level 
inline size_t BaseReaction::numVariableIndex(size_t level) const {
  return variableIndex_[level].size();
}

//!Returns a parameter
inline double BaseReaction::parameter(size_t i) const {
  return parameter_[i];
}

//!Returns a reference of a parameter
inline double& BaseReaction::parameterAddress(size_t i) {
  return parameter_[i];
}

//!Returns the id for a parameter
inline std::string BaseReaction::parameterId(size_t i) const {
  return parameterId_[i];
}

//!Returns the index of variable used in the update equation
inline size_t BaseReaction::variableIndex(size_t i,size_t j) const {
  return variableIndex_[i][j];
}

//!Sets the id string
inline void BaseReaction::setId(std::string value) {
  id_=value;
}

//!Sets the parameter value with index i
inline void BaseReaction::setParameter(size_t i,double value) {
  parameter_[i]=value;
}

//!Sets all parameters within a reaction
inline void BaseReaction::setParameter(std::vector<double> &value) {
  parameter_ = value;
}

//!Sets the parameter id value with index i
inline void BaseReaction::setParameterId(size_t i,std::string value) {
  parameterId_[i]=value;
}

//!Sets all parameter id's within a reaction
inline void BaseReaction::setParameterId(std::vector<std::string> &value) {
  parameterId_=value;
}

//!Sets the variable index at specific level and place
inline void BaseReaction::setVariableIndex(size_t i, size_t j,size_t value) {
  variableIndex_[i][j] = value;
}

//!Sets the variable indeces for a given level
inline void BaseReaction::setVariableIndex(size_t i, std::vector<size_t> &value) {
  variableIndex_[i] = value;
}

//!Sets the variable indeces used in a reaction
inline void BaseReaction::setVariableIndex(std::vector< std::vector<size_t> > 
					     &value) {
  variableIndex_ = value;
}

#endif


