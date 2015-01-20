// Filename     : MassAction.h
// Description  : Classes describing mass action reactions
// Author(s)    : Laura Brown (laura.brown@slcu.cam.ac.uk))
// Created      : August 2014
// Revision     : $Id:$
//
#ifndef MASSACTION_H
#define  MASSACTION_H

#include<cmath>

#include"tissue.h"
#include"baseReaction.h"



///
/// @brief A collect of mass action reactions including those acting on wall variables
///
/// These reactions follows mass action kinetics, i.e. the reaction rate is proportional to the reactants
/// concentrations. The general reaction handles any number of reactants and products, while specialized
/// versions a restricted to specific numbers of reactants/products. 
///
namespace MassAction {

  ///
  /// @brief A one way mass action reaction applied in all cells
  ///
  /// The reactant indices are in level zero and products in level one in
  /// variableIndex. One parameter (rate constant) is needed.
  ///
  /// Reactants [R] and products [P] stored in variableIndexLevel 0/1 are
  /// updated. The parameter is k_f in 
  ///
  /// @f[ \frac{d[P]}{dt} = k_f * \prod [R] @f]
  /// @f[ \frac{d[R]}{dt} = - k_f \prod [R] @f]
  ///
  /// In a model file the reaction is defined as:
  ///
  /// @verbatim
  /// MassAction::General 1 2 N_R N_P
  /// k_f
  /// R_1 ... R_{N_R}
  /// P_1 ... P_{N_P}
  /// @endverbatim
  ///
  class General : public BaseReaction {
    
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
    General(std::vector<double> &paraValue, 
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
    /// @brief Derivative function for this reaction class calculating the absolute value for noise solvers
    ///
    /// @see BaseReaction::derivsWithAbs(Compartment &compartment,size_t species,...)
    ///
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
/// @brief A one way mass action reaction assume one reactants and two products.
///
///
///  
///
/// @f[ \frac{dP}{dt} = k_1 R_1 R_2 @f] 
///  
/// @f[ \frac{dR_i}{dt} = - k_1 R_1 R_2 @f]
///
///
///  
/// In the model file the reaction is given by:
/// @verbatim
/// MassAction::OneToTwo 1 2 3 0
/// p_0
/// r_cell p1_cell P2_cell 
/// @endverbatim
///
class OneToTwo : public BaseReaction {
  
 public:
  
 OneToTwo(std::vector<double> &paraValue, 
		std::vector< std::vector<size_t> > 
		&indValue );
  
  void derivs(Tissue &T,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      DataMatrix &cellDerivs,
	      DataMatrix &wallDerivs,
	      DataMatrix &vertexDerivs );
};



///
/// @brief A one way mass action reaction assume two reactants and a single product
///
///
///  
///
/// @f[ \frac{dP}{dt} = k_1 R_1 R_2 @f] 
///  
/// @f[ \frac{dR_i}{dt} = - k_1 R_1 R_2 @f]
///
///
///  
/// In the model file the reaction is given by:
/// @verbatim
/// MassAction::TwoToOne 1 2 3 0
/// p_0
/// r1_cell r2_cell P_cell 
/// @endverbatim
///
class TwoToOne : public BaseReaction {
  
 public:
  
 TwoToOne(std::vector<double> &paraValue, 
		std::vector< std::vector<size_t> > 
		&indValue );
  
  void derivs(Tissue &T,
	      DataMatrix &cellData,
	      DataMatrix &wallData,
	      DataMatrix &vertexData,
	      DataMatrix &cellDerivs,
	      DataMatrix &wallDerivs,
	      DataMatrix &vertexDerivs );
};

  ///
  /// @brief A one way mass action reaction applied in all walls
  ///
  /// The reactant indices are in level zero and products in level one in
  /// variableIndex. One parameter (rate constant) is needed.
  ///
  /// Reactants [R] and products [P] stored in variableIndexLevel 0/1 are
  /// updated. The parameter is k_f in 
  ///
  /// @f[ \frac{d[P]}{dt} = k_f * \prod [R] @f]
  /// @f[ \frac{d[R]}{dt} = - k_f \prod [R] @f]
  ///
  /// In a model file the reaction is defined as:
  ///
  /// @verbatim
  /// MassAction::GeneralWall 1 2 N_R N_P
  /// k_f
  /// R_1 ... R_{N_R}
  /// P_1 ... P_{N_P}
  /// @endverbatim
  ///
  class GeneralWall : public BaseReaction {
    
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
    GeneralWall(std::vector<double> &paraValue, 
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
    /// @brief Derivative function for this reaction class calculating the absolute value for noise solvers
    ///
    /// @see BaseReaction::derivsWithAbs(Compartment &compartment,size_t species,...)
    ///
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
    /// @brief A one way mass action reaction assume one reactants and two products in walls.
    ///
    ///
    ///
    ///
    /// @f[ \frac{dP_ij}{dt} = k_1 R_1 R_2 @f]
    ///
    /// @f[ \frac{dR_ij}{dt} = - k_1 R_1 R_2 @f]
    ///
    ///
    ///
    /// In the model file the reaction is given by:
    /// @verbatim
    /// MassAction::OneToTwo 1 2 3 0
    /// p_0
    /// r_wall p1_wall P2_wall
    /// @endverbatim
    ///
    class OneToTwoWall : public BaseReaction {
        
    public:
        
        OneToTwoWall(std::vector<double> &paraValue,
                 std::vector< std::vector<size_t> >
                 &indValue );
        
        void derivs(Tissue &T,
                    DataMatrix &cellData,
                    DataMatrix &wallData,
                    DataMatrix &vertexData,
                    DataMatrix &cellDerivs,
                    DataMatrix &wallDerivs,
                    DataMatrix &vertexDerivs );
    };
    
    
    
    ///
    /// @brief A one way mass action reaction assume two reactants and a single product in walls
    ///
    ///
    ///
    ///
    /// @f[ \frac{dP}{dt} = k_1 R_1 R_2 @f]
    ///
    /// @f[ \frac{dR_i}{dt} = - k_1 R_1 R_2 @f]
    ///
    ///
    ///
    /// In the model file the reaction is given by:
    /// @verbatim
    /// MassAction::TwoToOneWall 1 2 3 0
    /// p_0
    /// r1_wall r2_wall P_wall
    /// @endverbatim
    ///
    class TwoToOneWall : public BaseReaction {
        
    public:
        
        TwoToOneWall(std::vector<double> &paraValue,
                 std::vector< std::vector<size_t> > 
                 &indValue );
        
        void derivs(Tissue &T,
                    DataMatrix &cellData,
                    DataMatrix &wallData,
                    DataMatrix &vertexData,
                    DataMatrix &cellDerivs,
                    DataMatrix &wallDerivs,
                    DataMatrix &vertexDerivs );
    };
    
    



}


#endif

