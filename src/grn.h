//
// Filename     : grn.h
// Description  : Classes describing gene regulatory network updates
//                Converted from Organism
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : February 2013
// Revision     : $Id: grn.h 561 2013-01-30 13:26:17Z henrik $
//
#ifndef GRN_H
#define GRN_H

#include <cmath>

#include "baseReaction.h"

///
/// @brief This class describes a Hill-type production with activators and
/// repressors.
///
/// @details This reaction is a Michaelis-Menten formalism for production with a
/// restricting variable given by
///
/// @f[ \frac{dy_{ij}}{dt} = p_0 \frac{y_{ik}^{p_2}}{p_1^{p_2}+y_{ik}^{p_2}}...
/// \frac{p_{1'}^{p_{2'}}}{p_{1'}^{p_{2'}}+y_{ik'}^{p_{2'}}}@f]
///
/// where p_0 is the maximal rate (@f$V_{max}@f$), and @f$p_1@f$ is the Hill
/// constant (@f$K_{half}@f$), and @f$p_2@f$ is the Hill coefficient (n). The
/// k index is given in the first level of varIndex and corresponds to the
/// activators, and the k' in the second level which corresponds to the
/// repressors. In both layers each index corresponds to a pair of parameters
/// (K,n) that are preceded by the @f$V_{max}@f$ parameter.
///
/// In the model file the reaction will be defined (somewhat different for different 
/// number of activators and repressors) e.g. as
///
/// One activator:
/// @verbatim
/// hill 3 3 1 1 0
/// V K n
/// produced_index
/// activator_index
/// @endverbatim
/// One repressor:
/// @verbatim
/// hill 3 3 1 0 1
/// V K n
/// produced_index
/// repressor_index
/// @endverbatim
/// Two activators and one repressor:
/// @verbatim
/// hill 5 3 1 2 1
/// V K_A1 n_A1 K_A2 n_A2 K_R n_R
/// produced_index
/// A1_index A2_index 
/// R_index
/// @endverbatim
///
class Hill : public BaseReaction {
  
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
  Hill(std::vector<double> &paraValue, 
       std::vector< std::vector<size_t> > &indValue );
  
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
        
};

///
/// @brief This class describes a simple Hill-type production where each bound/unbound state contributes
/// and for a single transcription factor
///
/// @details This reaction uses a single input and describes a Hill production (activation or repression) by
///
/// @f[ \frac{dy_{ij}}{dt} =
/// \frac{p_0 p_2^{p_3} + p_1 y_{ik}^{p_3}}{p_2^{p_3}+y_{ik}^{p_3}}@f]
///
/// where @f$ p_0 @f$ is the maximal unbound rate @f$ p_1 @f$ is the maximal bound rate
/// (@f$V_{max}@f$ if the other V is set to zero), and @f$ p_2 @f$ is the Hill
/// constant (@f$K_{half}@f$), and @f$ p_3 @f$ is the Hill coefficient. The
/// single k index is given in the first level of varIndex.
///
/// In a model file the reaction looks like
/// @verbatim
/// hillGeneralOne 4 2 1 1
/// Vunbound Vbound K_H n_H
/// j (produced_index)
/// k (activator/repressor_index)
/// @endverbatim
///
class HillGeneralOne : public BaseReaction {
  
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
  HillGeneralOne(std::vector<double> &paraValue, 
		 std::vector< std::vector<size_t> > &indValue );
  
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
};

///
/// @brief This class describes a Hill-type production where each bound/unbound state contributes
/// and for two transcription factors
///
/// @details This reaction uses two inputs and describes a Hill production (activation and repression) by
///
/// @f[ \frac{dy_{ij}}{dt} =
/// \frac{p_0 p_4^{p_5} p_6^{p_7} + p_1 y_{ik}^{p_5} p_6^{p_7} + p_2 p_4^{p_5} y_{il}^{p_7} + 
/// p_3 y_{ik}^{p_5} y_{il}^{p_7} }{ (p_4^{p_5}+y_{ik}^{p_5})(p_6^{p_7}+y_{ik}^{p_7})}@f]
///
/// where @f$ p_0 @f$ is the maximal unbound rate, @f$ p_1 @f$ is the maximal rate when the first TF is bound
/// , @f$ p_2 @f$ is the maximal rate when the second TF is bound, and @f$ p_3 @f$ is the maximal rate when
/// both TFs are bound (@f$V_{max}@f$ if the other Vs are set to zero). @f$ p_4,p_6 @f$ are the Hill
/// constants for the individual TFs (@f$K_{half}@f$), and @f$ p_5,p_7 @f$ are the Hill coefficients. The
/// two indices (@f$ k,l @f$) are given in the first level of varIndex.
///
/// In a model file the reaction looks like
/// @verbatim
/// hillGeneralTwo 8 2 1 2
/// V_unbound V_firstbound V_secondbound V_bothbound K_H1 n_H1 K_H2 n_H2 
/// j
/// k l
/// @endverbatim
///
class HillGeneralTwo : public BaseReaction {
  
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
  HillGeneralTwo(std::vector<double> &paraValue, 
		 std::vector< std::vector<size_t> > &indValue );
  
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
};

///
/// @brief This class describes a Hill-type production where each bound/unbound state contributes
/// and for three transcription factors
///
/// @details This reaction uses three inputs and describes a Hill production (activation and repression) by
///
/// @f[ \frac{dy_{ij}}{dt} =
/// \frac{p_0 p_8^{p_9} p_{10}^{p_{11}} p_{12}^{p_{13}} + p_1 y_{ik}^{p_9} p_{10}^{p_{11}} p_{12}^{p_{13}}
/// + p_2 p_8^{p_9} y_{jl}^{p_{11}} p_{12}^{p_{13}} + p_3 p_8^{p_9} p_{10}^{p_{11}} y_{jm}^{p_{13}}
/// + p_4 y_{jk}^{p_9} y_{il}^{p_{11}} p_{12}^{p_{13}} + p_5 y_{ik}^{p_9} p_{10}^{p_{11}} y_{im}^{p_{13}}
/// + p_6 p_8^{p_9} y_{il}^{p_{11}} y_{im}^{p_{13}} + p_7 y_{ik}^{p_9} y_{il}^{p_{11}} y_{im}^{p_{13}}  }
/// { (p_8^{p_9}+y_{ik}^{p_9})(p_{10}^{p_{11}}+y_{il}^{p_{11}})(p_{12}^{p_13}+y_{im}^{p_{13}})}@f]
///
/// where @f$ p_0 @f$ is the maximal unbound rate, @f$ p_1 @f$ is the maximal rate when the first TF is bound
/// , @f$ p_2 @f$ is the maximal rate when the second TF is bound, and @f$ p_3 @f$ is the maximal rate when
/// the third TF is bound. For two bound TFs the maximal rates are given by @f$ p_4 @f$ (first,second), 
/// @f$ p_5 @f$ (first,third), and @f$ p_6 @f$ (second,third), and maximal rate or production for all 
/// three TFs bound is given by @f$ p_7 @f$. 
/// @f$ p_8,p_{10},p_{12} @f$ are the Hill
/// constants for the individual TFs (@f$K_{half}@f$), and @f$ p_9,p_{11},p_{}13 @f$ are the Hill coefficients. The
/// three indices (@f$ k,l,m @f$) are given in the first level of varIndex.
///
/// In a model file the reaction looks like
/// @verbatim
/// hillGeneralThree 14 2 1 3
/// V_000 V_100 V_010 V_001 V_110 V_101 V_011 V_111 K_H1 n_H1 K_H2 n_H2 K_H3 n_H3 
/// j
/// k l m
/// @endverbatim
///
class HillGeneralThree : public BaseReaction {
  
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
  HillGeneralThree(std::vector<double> &paraValue, 
		   std::vector< std::vector<size_t> > &indValue );
  
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
};

///
/// @brief The class Grn use a neural network inspired mathematics for gene
/// regulation
///
/// This class uses a neural network inspired update with a sigmoidal function
/// for gene regulation as defined in Mjolsness et al (1991). The update is
/// done according to:
///
/// @f[ \frac{dy_{ij}}{dt} = \frac{1}{p_1}g\left( p_0 + \sum_k p_k y_{ik} \right) @f]
///
/// where g() is a sigmoidal function, the @f$p_k@f$ are the parameters given
/// from position 2 and forward, and the k indices for @f$y_{ik}@f$ are given
/// at the first level of variableIndex. Hence the number of parameters must
/// be two more than the number of indices given at the first level.
///
/// @see Grn::sigmoid(double x) for implementation of sigmoid function.
///
class Grn : public BaseReaction {
  
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
  /// @see BaseReaction::createReaction(std::vector<double> &paraValue,...)
 
  Grn(std::vector<double> &paraValue, 
      std::vector< std::vector<size_t> > &indValue );
  
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
  /// @brief Sigmoidal function used in the derivs function
  ///
  /// Sigmoidal function defined by
  ///
  /// @f[ g(x) = \frac{1}{2}\left( 1 + \frac{x}{\sqrt{1+x^2}}\right)@f]
  ///
  inline double sigmoid( double value );  
};



///
/// @brief The class Grn use a neural network inspired mathematics for gene
/// regulation where neighbor input is accounted for.
///
class Gsrn2 : public BaseReaction {
  
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
  /// @see BaseReaction::createReaction(std::vector<double> &paraValue,...)
  
  Gsrn2(std::vector<double> &paraValue, 
        std::vector< std::vector<size_t> > &indValue );
  
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
  /// @brief Sigmoidal function used in the derivs function
  ///
  /// Sigmoidal function defined by
  ///
  /// @f[ g(x) = \frac{1}{2}\left( 1 + \frac{x}{\sqrt{1+x^2}}\right)@f]
  ///
  inline double sigmoid( double value );  
};


inline double Grn::sigmoid(double x) 
{ 
  return 0.5*(1 + x/sqrt(1+x*x) );
}
inline double Gsrn2::sigmoid(double x) 
{ 
  return 0.5*(1 + x/sqrt(1+x*x) );
}



#endif
