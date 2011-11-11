/**
 * Filename     : baseCompartmentChange.h
 * Description  : The common base for classes describing compartmentChange updates
 * Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
 * Created      : April 2006
 * Revision     : $Id:$
 */
#ifndef BASECOMPARTMENTCHANGE_H
#define BASECOMPARTMENTCHANGE_H

#include<vector>
#include<string>
#include<iostream>
#include<fstream>
#include "myTypedefs.h"

class Tissue;
class Cell;

///
/// @brief A base class for classes defining updates relating to changes in tissue size, e.g. cell division
///
/// The BaseCompartmentChange class is a base class used when defining
/// different types of "compartmentChange" classes. Each compartmentChange class uses a
/// vector of parameters and variable indeces to calculate a derivative
/// of model variables. The variable indeces are divided into multiple
/// layers to allow for different types of contributions for different
/// variables. The baseCompartmentChange class can be seen as a 'factory'
/// creating different compartmentChanges of different types.
/// In addition to the constructor, new compartmentChanges classes should
/// define a flag function for determining the rule for when the update 
/// (e.g. cell division) should
/// be applied, and an update function defining what to do at the update.
///
/// @see Cell
///
class BaseCompartmentChange {
  
 private:
  
  std::string id_;
  int numChange_;
  std::vector<double> parameter_;           
  std::vector<std::string> parameterId_;           
  std::vector< std::vector<size_t> > variableIndex_;
  
 public:
  
  static BaseCompartmentChange* createCompartmentChange(std::vector<double> &paraValue, 
							std::vector< std::vector<size_t> > 
							&indValue,
							std::string idValue );
  static BaseCompartmentChange* createCompartmentChange( std::istream &IN ); 
  
  //Constructor/destructor not defined!
  //BaseCompartmentChange();
  virtual ~BaseCompartmentChange();
  
  BaseCompartmentChange & operator=( const BaseCompartmentChange & baseCompartmentChangeCopy );
  
  // Get values
  inline std::string id() const;
  inline int numChange() const;
  inline size_t numParameter() const;  
  inline size_t numVariableIndexLevel() const;
  inline size_t numVariableIndex(size_t level) const;
  
  inline double parameter(size_t i) const;
  inline double& parameterAddress(size_t i);
  inline std::string parameterId(size_t i) const;
  inline size_t variableIndex(size_t i,size_t j) const;
  inline std::vector<size_t>& variableIndex(size_t i);
  
  // Set values
  inline void setId(std::string value);
  inline void setNumChange(int val);
  inline void setParameter(size_t i,double value);
  inline void setParameter(std::vector<double> &value);
  inline void setParameterId(size_t i,std::string value);
  inline void setParameterId(std::vector<std::string> &value);
  inline void setVariableIndex(size_t i, size_t j,size_t value);
  inline void setVariableIndex(size_t i, std::vector<size_t> &value);
  inline void setVariableIndex(std::vector< std::vector<size_t> > &value);
  
  ///
  /// @brief Defines a rule for when an update will take place.
  ///
  virtual int flag(Tissue* T,size_t i,
		   DataMatrix &cellData,
		   DataMatrix &wallData,
		   DataMatrix &vertexData,
		   DataMatrix &cellDerivs,
		   DataMatrix &wallDerivs,
		   DataMatrix &vertexDerivs );
  ///
  /// @brief Defines a rule for how an update will update the Tissue
  ///
  virtual void update(Tissue* T,size_t i,
		      DataMatrix &cellData,
		      DataMatrix &wallData,
		      DataMatrix &vertexData,
		      DataMatrix &cellDerivs,
		      DataMatrix &wallDerivs,
		      DataMatrix &vertexDerivs );
  
  ///
  /// @brief Prints the wall positions and different other stuff for
  /// error checking in gnuplot.
  ///
  void printCellWallError(DataMatrix &vertexData,
			  Cell *divCell, 
			  std::vector<size_t> &w3Tmp, 
			  size_t &wI, 
			  size_t &w3I,
			  std::vector<double> &point,
			  std::vector<double> &normal,
			  std::ostream &os=std::cerr);
  
  ///
  /// @brief Delivers the wall indices and vertex positions of walls for division
  ///
  /// It uses a position and a direction to find two walls
  /// vertex positions for a division.
  ///
  int findTwoDivisionWalls(DataMatrix &vertexData, 
			   Cell *divCell, std::vector<size_t> &wI,
			   std::vector<double> &point, 
			   std::vector<double> &direction, 
			   std::vector<double> &v1Pos, 
			   std::vector<double> &v2Pos);
  
  ///
  /// @brief Delivers the index and position of a second wall for division
  ///
  /// It uses a position on a wall and a direction to find a second
  /// vertex position for a division.
  ///
  int findSecondDivisionWall(DataMatrix &vertexData, 
			     Cell *divCell, size_t &wI, size_t &w3I, 
			     std::vector<double> &v1Pos, 
			     std::vector<double> &nW2, 
			     std::vector<double> &v2Pos);
};

//!Returns the id string
inline std::string BaseCompartmentChange::id() const {
  return id_;
}

//!Returns the cell number change
inline int BaseCompartmentChange::numChange() const {
	return numChange_;
}

//!Returns the number of parameters
inline size_t BaseCompartmentChange::numParameter() const {
  return parameter_.size();
}

//!Returns the number of variable-index levels
inline size_t BaseCompartmentChange::numVariableIndexLevel() const {
  return variableIndex_.size();
}

//!Returns the number of variable-index within one level 
inline size_t BaseCompartmentChange::numVariableIndex(size_t level) const {
  return variableIndex_[level].size();
}

//!Returns a parameter
inline double BaseCompartmentChange::parameter(size_t i) const {
  return parameter_[i];
}

//!Returns a reference of a parameter
inline double& BaseCompartmentChange::parameterAddress(size_t i) {
  return parameter_[i];
}

//!Returns the id for a parameter
inline std::string BaseCompartmentChange::parameterId(size_t i) const {
  return parameterId_[i];
}

//!Returns the index of variable used in the update equation
inline size_t BaseCompartmentChange::variableIndex(size_t i,size_t j) const {
  return variableIndex_[i][j];
}

//!Returns the index vector of variables used in the update equation
inline std::vector<size_t>& BaseCompartmentChange::variableIndex(size_t i) {
  return variableIndex_[i];
}

//!Sets the id string
inline void BaseCompartmentChange::setId(std::string value) {
  id_=value;
}

//!Sets the numChange variable
inline void BaseCompartmentChange::setNumChange(int val) {
	numChange_=val;
}

//!Sets the parameter value with index i
inline void BaseCompartmentChange::setParameter(size_t i,double value) {
  parameter_[i]=value;
}

//!Sets all parameters within a compartmentChange
inline void BaseCompartmentChange::setParameter(std::vector<double> &value) {
  parameter_ = value;
}

//!Sets the parameter id value with index i
inline void BaseCompartmentChange::setParameterId(size_t i,std::string value) {
  parameterId_[i]=value;
}

//!Sets all parameter id's within a compartmentChange
inline void BaseCompartmentChange::setParameterId(std::vector<std::string> &value) {
  parameterId_=value;
}

//!Sets the variable index at specific level and place
inline void BaseCompartmentChange::setVariableIndex(size_t i, size_t j,size_t value) {
  variableIndex_[i][j] = value;
}

//!Sets the variable indeces for a given level
inline void BaseCompartmentChange::setVariableIndex(size_t i, std::vector<size_t> &value) {
  variableIndex_[i] = value;
}

//!Sets the variable indeces used in a compartmentChange
inline void BaseCompartmentChange::setVariableIndex(std::vector< std::vector<size_t> > 
					     &value) {
  variableIndex_ = value;
}

#endif


