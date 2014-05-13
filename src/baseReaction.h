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
#include"myTypedefs.h"

class Tissue;

///
/// @brief A factory class for classes describing differential equation
/// updates for dynamical variables
///
/// The BaseReaction class is a base (factory) class used when defining
/// different types of "reaction" classes. Each reaction class uses a vector
/// of parameters and variable indices to calculate a derivative of model
/// variables. The variable indices are divided into multiple layers to allow
/// for different types of contributions for different variables. The
/// BaseReaction class can be seen as a 'factory' creating reactions of
/// different types.
///
class BaseReaction {
  
 private:
  
  std::string id_;
  std::vector<double> parameter_;           
  std::vector<std::string> parameterId_;           
  std::vector< std::vector<size_t> > variableIndex_;
  
 public:
  
  ///
  /// @brief Main factory creator, all creation should be mapped onto this one
  ///
  /// Given the idValue a reaction of the defined type is returned (using new
  /// Class). It chooses from the user defined list of possible reactions, and
  /// returns the correct one if defined. If the idValue given is not a defined
  /// reaction it exits.
  ///
  /// @param paraValue vector with parameters
  ///
  /// @param indValue vector of vectors with variable indices
  ///
  /// @param idValue identification of which reaction that should be created
  ///
  /// @return Returns a pointer to an instance of a reaction class defined by
  /// the idValue string.
  ///
  static BaseReaction* createReaction(std::vector<double> &paraValue, 
				      std::vector< std::vector<size_t> > 
				      &indValue,
				      std::string idValue );
  ///
  /// @brief Reads reaction data from an open file into data structures and 
  /// calls createReaction(paraValue,...)
  ///
  /// @see createReaction(std::vector<double>&,...)
  /// 
  static BaseReaction* createReaction( std::istream &IN ); 
  
  //Constructor/destructor not defined!
  //BaseReaction();
  virtual ~BaseReaction();
  
  BaseReaction & operator=( const BaseReaction & baseReactionCopy );
  
  // Get values
  ///
  /// @brief Returns the id string
  ///
  /// The id (name) string is used to define the reaction in the model file
  /// and it has to be a perfect match.
  ///
  /// @see createReaction(std::istream&)
  inline std::string id() const;
  ///
  /// @brief Returns the number of parameters
  ///
  /// @see parameter(size_t)
  ///
  inline size_t numParameter() const;  
  ///
  /// @brief Returns the number of variable-index levels
  ///
  /// @see variableIndex(size_t,size_t)
  ///
  inline size_t numVariableIndexLevel() const;
  ///
  /// @brief Returns the number of variable-index within the given level 
  ///
  /// @see variableIndex(size_t,size_t)
  ///
  inline size_t numVariableIndex(size_t level) const;  
  ///
  /// @brief Returns a parameter
  ///
  /// The developer has to define how many parameters that are needed for a specific reaction.
  /// Then it is checked such that the user provides this many parameters in the model file.
  ///
  inline double parameter(size_t i) const;
  ///
  /// @brief Returns a reference to a specific parameter
  ///
  /// Used such that the parameter can be updated from outside the class.
  ///
  /// @see parameter(size_t)
  ///
  inline double& parameterRef(size_t i);
  ///
  /// @brief Returns the id for a parameter
  ///
  /// Each parameter has a name defined, only used for information.
  ///
  /// @see parameter(size_t)
  ///
  inline std::string parameterId(size_t i) const;
  ///
  /// @brief Returns the index of variable used in the update equation
  ///
  /// A reaction may depend on the state of different (user defined) variables.
  /// The variable indices are stored in multiple vectors (if needed). The 
  /// indices relates to the indices in the data structure of Cells, Walls and
  /// Vertices, and the user provides these in the model file when creating a
  /// reaction. The developer of a new reaction defines which variable indices 
  /// that are needed for the reaction and in which order they should be given.
  ///
  /// @see Cell
  /// @see Wall
  /// @see Vertex
  ///
  inline size_t variableIndex(size_t i,size_t j) const;
  
  // Set values
  ///
  /// @brief Sets the id string
  ///
  /// @see id()
  ///
  inline void setId(std::string value);
  ///
  /// @brief Sets the parameter value with index i
  ///
  /// @see parameter(size_t)
  ///
  inline void setParameter(size_t i,double value);
  ///
  /// @brief Sets all parameters within a reaction
  ///
  /// @see parameter(size_t)
  ///
  inline void setParameter(std::vector<double> &value);
  ///
  /// @brief Sets the parameter id value with index i
  ///
  /// @see parameterId(size_t)
  ///
  inline void setParameterId(size_t i,std::string value);
  ///
  /// @brief Sets all parameter id's within a reaction.
  ///
  /// @see parameterId(size_t)
  ///
  inline void setParameterId(std::vector<std::string> &value);
  ///
  /// @brief Sets the variable index at specific level and place.
  ///
  /// @see variableIndex(size_t,size_t)
  ///
  inline void setVariableIndex(size_t i, size_t j,size_t value);
  ///
  /// @brief Sets the variable indices for a given level.
  ///
  /// @see variableIndex(size_t,size_t)
  ///
  inline void setVariableIndex(size_t i, std::vector<size_t> &value);
  ///
  /// @brief Sets all the variable indices used in a reaction.
  ///
  /// @see variableIndex(size_t,size_t)
  ///
  inline void setVariableIndex(std::vector< std::vector<size_t> > &value);
  ///
  /// @brief For reactions that need to initiate its parameters (state) before
  /// the simulation
  ///
  /// Occasionally, a reaction have 'states' included as parameters. If these
  /// need to be initiated before the simulation this is done by this
  /// function. Initiation of the state variables (overriding the read initial 
  /// state) can also be done. 
  /// Note: if the initiate function is not defined for a specific reaction
  /// this virtual function will be called and it does nothing.
  ///
  /// @param T The Tissue is provided such that connection information and other
  /// functions defined in the Tissue, Cells, Walls, and Vertices can be used.
  ///
  /// @param cellData The current state for the cells used for the update.
  /// @param wallData The current state for the walls used for the update.
  /// @param vertexData The current state for the vertices used for the update.
  /// @param cellDerivs The cell derivs to which the output of the reaction is added.
  /// @param wallDerivs The wall derivs to which the output of the reaction is added.
  /// @param vertexDerivs The vertex derivs to which the output of the reaction is added.
  ///
  /// @see Tissue::initiateReactions()
  /// @see Tissue::readInit(std::istream &,int)
  /// @note This is typically only included for quite specific reactions.
  ///  
  virtual void initiate(Tissue &T,
			DataMatrix &cellData,
			DataMatrix &wallData,
			DataMatrix &vertexData,
			DataMatrix &cellDerivs,
			DataMatrix &wallDerivs,
			DataMatrix &vertexDerivs );
  ///
  /// @brief Calculates the derivative given the state in
  /// cell[wall,vertex]Data and adds it to cell[wall,vertex]Derivs.
  ///
  /// This function, together with the constructor are the main functions
  /// needed to create a new reaction class. It calculates the time derivative
  /// contribution from a reaction given the state in cellData, wallData, and vertexData
  /// via the indices given
  /// in the member variableIndex and parameters given in the member parameter
  /// and adds the result to cellDerivs, wallDerivs, and vertexDerivs. 
  /// The individual reaction derivs() functions are called from Tissue::derivs().
  /// Note: The reactions need to loop over the tissue elements themselves, since this is
  /// not done in Tissue::derivs() (unlike the similar setup in Organism::derivs()).
  /// Note: If the derivs function is not defined for a specific reaction, this virtual
  /// function will be called, which will exit the program.
  ///
  /// @see Tissue::derivs()
  ///
  /// @param T The Tissue is provided such that connection information and other
  /// functions defined in the Tissue, Cells, Walls, and Vertices can be used.
  ///
  /// @param cellData The current state for the cells used for the update.
  /// @param wallData The current state for the walls used for the update.
  /// @param vertexData The current state for the vertices used for the update.
  /// @param cellDerivs The cell derivs to which the output of the reaction is added.
  /// @param wallDerivs The wall derivs to which the output of the reaction is added.
  /// @param vertexDerivs The vertex derivs to which the output of the reaction is added.
  ///
  virtual void derivs(Tissue &T,
		      DataMatrix &cellData,
		      DataMatrix &wallData,
		      DataMatrix &vertexData,
		      DataMatrix &cellDerivs,
		      DataMatrix &wallDerivs,
		      DataMatrix &vertexDerivs );
  ///
  /// @brief For reactions that need to update its parameters (state) during
  /// the simulation
  ///
  /// Occasionally, a reaction have 'states' included as parameters. If these
  /// need to be updated during the simulation this is done by this
  /// function. Update of the state variables can also be done. Also for
  /// example some approximate stochastic updates can be done in this
  /// function. This is done in between normal derivative steps in the
  /// solver. Note: if the update function is not defined for a specific reaction
  /// this virtual function will be called and it does nothing.
  ///
  /// @param h The current time step taken by the solver.
  /// @param T The Tissue is provided such that connection information and other
  /// functions defined in the Tissue, Cells, Walls, and Vertices can be used.
  ///
  /// @param cellData The current state for the cells used for the update.
  /// @param wallData The current state for the walls used for the update.
  /// @param vertexData The current state for the vertices used for the update.
  ///
  /// @see Tissue::updateReactions()
  /// @note This is typically only included for very specific reactions.
  ///
  virtual void update(Tissue &T,
		      DataMatrix &cellData,
		      DataMatrix &wallData,
		      DataMatrix &vertexData,
		      double h);
  ///
  /// @brief Prints the data structure of a reaction.
  ///
  /// Prints the data structure in a format readable for (re)creating a reaction.
  /// It is defined as virtual for the BaseReaction class (and will exit if not
  /// defined for the specific reaction) and relies upon that it is properly
  /// defined for the specific reaction to give a good output.
  ///
  /// @note In the current implementation, the print function is typically not used.
  ///
  virtual void print( std::ofstream &os );

  ///
  /// @brief Prints internal variables stored by a reaction
  ///
  /// This allows for printing internal variable/parameter values as specified
  /// by individual reactions. The BaseReaction version does not do anything.
  /// 
  virtual void printState(Tissue *T,
			  DataMatrix &cellData,
			  DataMatrix &wallData,
			  DataMatrix &vertexData, 
			  std::ostream &os=std::cout);
};

inline std::string BaseReaction::id() const {
  return id_;
}

inline size_t BaseReaction::numParameter() const {
  return parameter_.size();
}

inline size_t BaseReaction::numVariableIndexLevel() const {
  return variableIndex_.size();
}

inline size_t BaseReaction::numVariableIndex(size_t level) const {
  return variableIndex_[level].size();
}

inline double BaseReaction::parameter(size_t i) const {
  return parameter_[i];
}

inline double& BaseReaction::parameterRef(size_t i) {
  return parameter_[i];
}

inline std::string BaseReaction::parameterId(size_t i) const {  return parameterId_[i];
}

inline size_t BaseReaction::variableIndex(size_t i,size_t j) const {
  return variableIndex_[i][j];
}

inline void BaseReaction::setId(std::string value) {
  id_=value;
}

inline void BaseReaction::setParameter(size_t i,double value) {
  parameter_[i]=value;
}

inline void BaseReaction::setParameter(std::vector<double> &value) {
  parameter_ = value;
}

inline void BaseReaction::setParameterId(size_t i,std::string value) {
  parameterId_[i]=value;
}

inline void BaseReaction::setParameterId(std::vector<std::string> &value) {
  parameterId_=value;
}

inline void BaseReaction::setVariableIndex(size_t i, size_t j,size_t value) {
  variableIndex_[i][j] = value;
}

inline void BaseReaction::setVariableIndex(size_t i, std::vector<size_t> &value) {
  variableIndex_[i] = value;
}

inline void BaseReaction::setVariableIndex(std::vector< std::vector<size_t> > 
					   &value) {
  variableIndex_ = value;
}

#endif


