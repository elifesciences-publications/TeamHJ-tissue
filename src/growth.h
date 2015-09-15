//
// Filename     : growth.h
// Description  : Classes describing growth updates
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : April 2006
// Revision     : $Id:$
//
#ifndef WALLGROWTH_H
#define WALLGROWTH_H

#include<cmath>

#include"tissue.h"
#include"baseReaction.h"

/// 
/// @brief Reactions describing wall growth
///
/// Different reactions for wall growth implemented as increase of the resting length of the walls.
/// Either a constant or strain/stress based approach is used, and in some cases the possibility to
/// truncate the growth at specific lengths can be applied. Also versions reading cell variables
/// (typically a concentration) that affects the growth rate are implemented and named ConcentrationHill.
/// In all cases the index of the resting length is given and this is typically index 0.
///
/// In some cases called stress, strain or stress (read from wall variables) can be used and selected by a flag.
/// In some cases one can choose between const (no) or linear dependence on the current length of the wall,
/// again given by a flag.
///
namespace WallGrowth {

  ///
  /// @brief Constant/exponential wall growth which can be truncated at threshold length
  ///
  /// Constant or exponential (proportional to length) growth that can be truncated at 
  /// a threshold length. 
  /// The wall lengths, L, are updated according to
  ///
  /// @f[ \frac{dL}{dt} = p_0 L (1-\frac{L}{p_2}) @f]
  ///
  /// if linearFlag (@f$p_1@f$) is set, or 
  ///
  /// @f[ \frac{dL}{dt} = p_0 (1-\frac{L}{p_2}) @f]
  ///
  /// otherwise. @f$ p_0 @f$ is the growth rate (@f$ k_{growth} @f$), @f$ p_1 @f$ is the flag setting
  /// the update version (0 for constant and 1 for linear dependence on the current length), and 
  /// the optional parameter @f$ p_2 @f$ is the threshold for growth truncation (@f$ L_{trunc} @f$).
  /// If this parameter is not given, no truncation is applied.
  /// The column index for the wall length should be given at the first level.
  ///
  /// In a model file the reaction is defined as
  ///
  /// @verbatim
  /// WallGrowth::Constant 2/3 1 1
  /// k_growth linear_flag [L_trunc]
  /// L
  /// @endverbatim
  ///
  class Constant : public BaseReaction {
  
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
  Constant(std::vector<double> &paraValue, 
	   std::vector< std::vector<size_t> > 
	   &indValue );
  ///
  /// @brief Derivative function for this reaction class
  ///
  /// @see BaseReaction::derivs(Tissue &T,...)
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
  /// @brief Stress/strain-driven wall growth dependent on a threshold.
  ///
  /// Constant (constant mode) or exponential (linear mode) growth driven by a 
  /// streched/stressed wall. The wall 
  /// lengths are updated only if the length plus a threshold value is shorter 
  /// than the distance between the vertices of the wall (in strain mode), and if 
  /// the total stress is above the threshold (in stress mode). Hence the update follows:
  ///
  ///  @f[ \frac{dL}{dt} = p_0 L (S-p_1) @f] if @f$ S > p_1 @f$ and 0 otherwise.
  ///
  /// where S is the stress/strain and L is the wall length.
  /// @f$ p_0 @f$ is the growth rate (@f$ k_{growth} @f$). 
  /// @f$ p_1 @f$ is a threshold (@f$ s_{threshold} @f$) in stress or strain depending 
  /// on @f$ p_2 @f$. If set to zero, shinkage is allowed. 
  /// @f$ p_2 @f$ is a flag (@f$ strain_{flag} @f$) for using stretch/strain (@f$p_2=1 @f$) 
  /// or stress (@f$p_2=0 @f$). Strain is calculated by @f$ S=(d-L)/L @f$, where d is the distance
  /// between the vertices. Stress is read from the second layer of variable indices.  
  /// @f$ p_3 @f$ is a flag (@f$ linear_{flag} @f$) for using growth proportional to 
  /// wall length (@f$ p_3=1 @f$, as in equation above. If @f$ p_3=0 @f$ the above equation
  /// will not include the length L on the right-hand side.
  /// If a 5th parameter is given, the growth is truncated at a maximal length by multiplying the 
  /// rate with a factor @f$ (1 - L/p_4) @f$.
  ///  The column index for the wall length should be given in the first layer, and if stress
  ///  is to be used, a second layer with calculated wall stress variables has to be given.
  ///
  /// In a model file the reaction is defined as
  ///
  /// @verbatim
  /// WallGrowthStress 4/5 1/2 1 [N]
  /// k_growth s_threshold stretch_flag linear_flag [L_trunc]  
  /// L
  /// [stress1 ... stressN]
  /// @endverbatim
  ///
  /// If stress is used (stretch_flag=0) a second level of wall stresses has to be read
  /// (calculated and updated from other (mechanical) reactions).
  ///
  /// @note If s_threshold is set to zero, also shrinkage is allowed. To avoid shrinkage set small value.
  ///
  class Stress : public BaseReaction {
    
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
    Stress(std::vector<double> &paraValue, 
	   std::vector< std::vector<size_t> > 
	   &indValue );
    ///
    /// @brief Derivative function for this reaction class
    ///
    /// @see BaseReaction::derivs(Tissue &T,...)
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
  /// @brief almansi-strain-driven wall growth dependent on a threshold via "update".
  ///
  /// In a model file the reaction is defined as
  ///
  /// @verbatim
  /// WallGrowthStrain 4 1 1 
  /// k_growth s_threshold strain_flag linear_flag   
  /// L0-index
  /// @endverbatim
  ///
  /// If Almansi strain is used (strain_flag=0) .
  ///
  /// @note If s_threshold is set to zero, also shrinkage is allowed. To avoid shrinkage set small value.
  ///
  class Strain : public BaseReaction {
    
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
    Strain(std::vector<double> &paraValue, 
	   std::vector< std::vector<size_t> > 
	   &indValue );
    ///
    /// @brief Derivative function for this reaction class
    ///
    /// @see BaseReaction::derivs(Tissue &T,...)
    ///  
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
                double h );
  };
  
    
    
    
    
    ///
    /// @brief Constant stress/strain-driven wall growth dependent on a
    /// distance to an maximal coordinate (e.g. tip)
    ///
    /// Constant growth driven by a streched wall. The wall lengths, L, are
    /// updated only if the length is shorter than the distance between the
    ///  vertices of the wall and then according to
    ///
    ///  dL/dt = p_0*(d_v-L-p_1)*p_2^p_3/(p_2^p_3+d^p_3) iff (d_v-L) > p_1
    ///
    /// p_0 is the growth rate.
    /// p_1 is a stress/strain threshold
    /// p_2 is the K_Hill of the spatial factor
    /// p_3 is the n_Hill of the spatial factor
    /// p_4 is a flag for using stretch/strain instead of stress
    /// p_5 is a flag for using growth proportional to wall length (not constant)
    /// d_v is the distance between the two wall vertices.
    /// d is the distance between the max value and wall
    ///
    ///  In addition, the column index for the wall length, the distance
    ///  coordinate should be given at first level and stress index in second.
    ///
    class StressSpatial : public BaseReaction {
        
    private:
        
        double Kpow_;
        
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
        StressSpatial(std::vector<double> &paraValue,
                      std::vector< std::vector<size_t> >
                      &indValue );
        ///
        /// @brief Derivative function for this reaction class
        ///
        /// @see BaseReaction::derivs(Tissue &T,...)
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
    /// @brief Constant stress/strain-driven wall growth dependent on a
    /// distance to an maximal coordinate (e.g. tip)
    ///
    /// Constant growth driven by a streched wall. The wall lengths, L, are
    /// updated only if the length is shorter than the distance between the
    ///  vertices of the wall and then according to
    ///
    ///  dL/dt = p_0*(d_v-L-p_1)*p_2^p_3/(p_2^p_3+d^p_3) iff (d_v-L) > p_1
    ///
    /// p_0 is the growth rate.
    /// p_1 is a stress/strain threshold
    /// p_2 is the K_Hill of the spatial factor
    /// p_3 is the n_Hill of the spatial factor
    /// p_4 is a flag for using stretch/strain instead of stress
    /// p_5 is a flag for using growth proportional to wall length (not constant)
    /// d_v is the distance between the two wall vertices.
    /// d is the distance between the max value and wall
    ///
    ///  In addition, the column index for the wall length, the distance
    ///  coordinate should be given at first level and stress index in second.
    ///
    class StressSpatialSingle : public BaseReaction {
        
    private:
        
        double Kpow_;
        
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
        StressSpatialSingle(std::vector<double> &paraValue,
                            std::vector< std::vector<size_t> >
                            &indValue );
        ///
        /// @brief Derivative function for this reaction class
        ///
        /// @see BaseReaction::derivs(Tissue &T,...)
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
    /// @brief Constant stress/strain-driven wall growth dependent on a
    /// concentration level in the cell
    ///
    /// Constant growth driven by a streched wall. The wall lengths, L, are
    /// updated only if the length is shorter than the distance between the
    ///  vertices of the wall and then according to
    ///
    ///  dL/dt = (p_0+p_1*f(c,p_2,p_3)) * (d_v-L) if (d_v-L) > p_4
    ///
    /// p_0 is the constant growth rate.
    /// p_1 is the maximal added growth rate (V_max in the Hill function)
    /// p_2 is the K_Hill
    /// p_3 is the n_Hill
    /// p_4 is a threshold
    /// p_5 is a flag for using stretch/strain instead of stress
    /// p_6 is a flag for using growth proportional to wall length (not constant)
    /// f is the hill function (increasing)
    /// c is the concentration
    /// d_v is the distance between the two wall vertices.
    ///
    ///
    ///  In addition, the column index for the wall length and the cell
    ///  concentration should be given.
    ///
    class StressConcentrationHill : public BaseReaction {
        
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
        StressConcentrationHill(std::vector<double> &paraValue,
                                std::vector< std::vector<size_t> >
                                &indValue );
        ///
        /// @brief Derivative function for this reaction class
        ///
        /// @see BaseReaction::derivs(Tissue &T,...)
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
    /// @brief Constant strech-driven wall growth with epidermal walls treated specially
    ///
    ///Constant growth driven by a streched wall. The wall lengths, L, are
    ///updated only if the length is shorter than the distance between the
    ///vertices of the wall and then according to
    ///
    ///dL/dt = p_0*(d_v-L)*f_e
    ///
    ///p_0 is the growth rate.
    ///f_e fraction for epidermal walls
    ///d_v is the distance between the two wall vertices.
    ///
    ///In addition, the column index for the wall length should be given.
    ///
    class ConstantStressEpidermalAsymmetric : public BaseReaction {
        
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
        ConstantStressEpidermalAsymmetric(std::vector<double> &paraValue,
                                          std::vector< std::vector<size_t> >
                                          &indValue );
        ///
        /// @brief Derivative function for this reaction class
        ///
        /// @see BaseReaction::derivs(Tissue &T,...)
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
    /// @brief Reads Force(s) from wall variables and increase resting length if larger than
    /// threshold
    ///
    /// Implements growth as
    ///
    ///  @f[ \frac{dL}{dt} = p_0 (F-p_1) @f]  if @f$ F > p_1 @f$ and 0 otherwise.
    ///
    /// where F is added up from wall variables given as indices. In am model file the reaction
    /// is given as:
    ///
    /// @verbatim
    /// WallGrowth::Force 2 2 1 N
    /// k_growth F_threshold
    /// L
    /// F_1 ... F_N
    /// @endverbatim
    ///
    /// This is an old reaction that is replaced by the more general WallGrowth::Stress(), but
    /// hangs around since it was used in publications.
    ///
    /// @note used to be called WallLengthGrowExperimental.
    /// @see WallGrowth::Stress()
    ///
    class Force : public BaseReaction {
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
        Force(std::vector<double> &paraValue,
              std::vector< std::vector<size_t> >
              &indValue );
        /// @brief Derivative function for this reaction class
        ///
        /// @see BaseReaction::derivs(Tissue &T,...)
        ///
        void derivs(Tissue &T,
                    DataMatrix &cellData,
                    DataMatrix &wallData,
                    DataMatrix &vertexData,
                    DataMatrix &cellDerivs,
                    DataMatrix &wallDerivs,
                    DataMatrix &vertexDerivs);
    };
    
    ///
    /// @brief nested name space
    ///

  namespace CenterTriangulation {
    ///
    /// @brief Constant internal edge growth which can be truncated at threshold length
    ///
    /// Constant (constant mode) or exponential (proportional to length, linear mode) 
    /// growth of the internal edges 
    /// in a central cell vertex meshed description. The internal edges 
    /// lengths are updated only if the length plus a threshold value is shorter 
    /// than the distance between the vertices (vertex and central cell vertex). 
    ///
    /// The internal edge lengths, L, are updated according to
    ///
    /// @f[ \frac{dL}{dt} = p_0 L (1-\frac{L}{p_2}) @f]
    ///
    /// if linearFlag (@f$p_1@f$) is set, or 
    ///
    /// @f[ \frac{dL}{dt} = p_0 (1-\frac{L}{p_2}) @f]
    ///
    /// otherwise. @f$ p_0 @f$ is the growth rate (@f$ k_{growth} @f$), @f$ p_1 @f$ is the flag setting
    /// the update version (0 for constant and 1 for linear dependence on the current length), and 
    /// the optional parameter @f$ p_2 @f$ is the threshold for growth truncation (@f$ L_{trunc} @f$).
    /// If this parameter is not given, no truncation is applied.
    ///
    /// The column index for the cell additional variables of the central mesh 
    /// (x,y,z,L_1,...,L_n) should be given in the first level of indices.
    ///
    /// In a model file the reaction is defined as
    ///
    /// @verbatim
    /// WallGrowth::Constant 2/3 1 1
    /// k_growth linear_flag [L_trunc]
    /// index
    /// @endverbatim
    ///
    /// @see WallGrowth::Constant (for same update of 1D walls)
    ///
    class Constant : public BaseReaction {
      
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
      Constant(std::vector<double> &paraValue, 
	       std::vector< std::vector<size_t> > 
	       &indValue );
      ///
      /// @brief Derivative function for this reaction class
      ///
      /// @see BaseReaction::derivs(Tissue &T,...)
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
    /// @brief Constant stress/strain-driven internal edge growth dependent on a threshold
    ///
    /// Constant (constant mode) or exponential (linear mode) growth of the internal edges 
    /// in a central cell vertex meshed description and driven by a 
    /// streched/stressed wall. The internal edges 
    /// lengths are updated only if the length plus a threshold value is shorter 
    /// than the distance between the vertices (vertex and central cell vertex) 
    /// (in strain mode), and if 
    /// the total stress is above the threshold (in stress mode). Hence the update follows:
    ///
    ///  @f[ \frac{dL}{dt} = p_0 L (S-p_1) @f] if @f$ S > p_1 @f$ and 0 otherwise.
    ///
    /// where S is the stress/strain and L is the wall length.
    /// @f$ p_0 @f$ is the growth rate (@f$ k_{growth} @f$). 
    /// @f$ p_1 @f$ is a threshold (@f$ s_{threshold} @f$) in stress or strain depending 
    /// on @f$ p_2 @f$. 
    /// @f$ p_2 @f$ is a flag (@f$ strain_{flag} @f$) for using stretch/strain (@f$p_2=1 @f$) 
    /// or stress (@f$p_2=0 @f$). Strain is calculated by @f$ S=(d-L)/L @f$, where d is the distance
    /// between the vertices. Stress is read from the second layer of variable indices.  
    /// @f$ p_3 @f$ is a flag (@f$ linear_{flag} @f$) for using growth proportional to 
    /// wall length (@f$ p_3=1 @f$, as in equation above. If @f$ p_3=0 @f$ the above equation
    /// will not include the length L on the right-hand side.
    /// If a 5th parameter is given, the growth is truncated at a maximal length by multiplying the 
    /// rate with a factor @f$ (1 - L/p_4) @f$.
    ///
    /// The column index for the cell additional variables of the central mesh (x,y,z,L_1,...,L_n) 
    /// should be given in the first layer, and if stress
    /// is to be used, a second layer with calculated wall stress variables has to be given.
    ///
    /// In a model file the reaction is defined as
    ///
    /// @verbatim
    /// CenterTriangulation::WallGrowth::Stress 4/5 1/2 1 [N]
    /// k_growth s_threshold stretch_flag linear_flag [L_trunc] 
    /// L
    /// [stress1 ... stressN]
    /// @endverbatim
    ///
    /// If stress is used (stretch_flag=0) a second level of wall stresses has to be read
    /// (calculated and updated from other (mechanical) reactions).
    ///
    /// @see WallGrowthStress (for same updatre of 1D walls)
    ///
    class Stress : public BaseReaction {
      
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
      Stress(std::vector<double> &paraValue, 
	     std::vector< std::vector<size_t> > 
	     &indValue );
      ///
      /// @brief Derivative function for this reaction class
      ///
      /// @see BaseReaction::derivs(Tissue &T,...)
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
    /// @brief 
    ///
    /// In a model file the reaction is defined as
    ///
    /// @verbatim
    /// CenterTriangulation::WallGrowth::StrainTRBS 2 3 1 1 3
    /// k_growth s_threshold  
    ///    
    /// L_ij-index
    /// InternalVarStartIndex 
    ///
    /// strain1_index
    /// strain2_index
    /// strain_vector_index
    /// @endverbatim
    ///
    /// (strain value and direction calculated and updated from other (mechanical) reactions).
    ///
    ///

    class StrainTRBS : public BaseReaction {
      
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
      StrainTRBS(std::vector<double> &paraValue, 
	     std::vector< std::vector<size_t> > 
	     &indValue );
      ///
      /// @brief Derivative function for this reaction class
      ///
      /// @see BaseReaction::derivs(Tissue &T,...)
      ///  
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
                  double h );
    };

  

  }
 
}

///
/// @brief Growth via vertex movement radially outwards
///
/// The tissue grows from vertex movement radially outwards. The update can be
///
/// @f[ \frac{dr}{dt} = p_{0} @f] (if @f$ p_1=0 @f$) or
/// @f[ \frac{dr}{dt} = p_{0} r @f] (if @f$ p_{1}=1 @f$)
///
/// @f$ p_{0} @f$ is the rate (@f$ k_{growth} @f$),
/// @f$ p_{1} @f$ {0,1} is a flag determining which function to be used (@f$ r_{pow} @f$).
///
/// In a model file the reaction is defined as
///
/// @verbatim
/// MoveVertexRadially 2 0
/// p_0 p_1
/// @endverbatim
///
class MoveVertexRadially : public BaseReaction {
    
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
    MoveVertexRadially(std::vector<double> &paraValue,
                       std::vector< std::vector<size_t> >
                       &indValue );
    ///
    /// @brief Derivative function for this reaction class
    ///
    /// @see BaseReaction::derivs(Tissue &T,...)
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
/// @brief Growth via vertex movement radially outwards
///
/// The tissue grows the epidermal cells (cells bordering to the background)
/// from vertex movement radially outwards. The update can be
///
/// @f[ \frac{dr}{dt} = p_{0} @f] (if @f$ p_1=0 @f$) or
/// @f[ \frac{dr}{dt} = p_{0} r @f] (if @f$ p_{1}=1 @f$)
///
/// @f$ p_{0} @f$ is the rate (@f$ k_{growth} @f$),
/// @f$ p_{1} @f$ {0,1} is a flag determining which function to be used (@f$ r_{pow} @f$).
///
/// In a model file the reaction is defined as
///
/// @verbatim
/// MoveVertexRadially 2 0
/// p_0 p_1
/// @endverbatim
///
class MoveEpidermalVertexRadially : public BaseReaction {
    
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
    MoveEpidermalVertexRadially(std::vector<double> &paraValue,
                       std::vector< std::vector<size_t> >
                       &indValue );
    ///
    /// @brief Derivative function for this reaction class
    ///
    /// @see BaseReaction::derivs(Tissue &T,...)
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
/// @brief Growth via vertex movement in the x-direction
///
/// The tissue grows from vertex movement outwards in the x-direction. The update can be
///
/// @f[ \frac{dx}{dt} = p_{0} @f] (if @f$ p_{1}=0 @f$) or
/// @f[ \frac{dx}{dt} = p_{0} x @f] (if @f$ p_{1}=1 @f$)
///
/// @f$ p_{0} @f$ is the rate (@f$ k_{growth} @f$),
/// @f$ p_{1} @f$ {0,1} is a flag determining which function to be used.
///
/// In a model file the reaction is defined as
///
/// @verbatim
/// MoveVerteX 2 0
/// p_0 p_1
/// @endverbatim
///
class MoveVerteX : public BaseReaction {
    
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
    MoveVerteX(std::vector<double> &paraValue,
               std::vector< std::vector<size_t> >
               &indValue );
    ///
    /// @brief Derivative function for this reaction class
    ///
    /// @see BaseReaction::derivs(Tissue &T,...)
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
/// @brief Growth via vertex movement in the y-direction
///
/// The tissue grows from vertex movement outwards in the x-direction. The update can be
///
/// @f[ \frac{dy}{dt} = p_{0} @f] (if @f$ p_{1}=0 @f$) or
/// @f[ \frac{dy}{dt} = p_{0} x @f] (if @f$ p_{1}=1 @f$)
///
/// @f$ p_{0} @f$ is the rate (@f$ k_{growth} @f$),
/// @f$ p_{1} @f$ {0,1} is a flag determining which function to be used.
///
/// In a model file the reaction is defined as
///
/// @verbatim
/// MoveVertexY 2 0
/// p_0 p_1
/// @endverbatim
///
class MoveVertexY : public BaseReaction {
    
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
    MoveVertexY(std::vector<double> &paraValue,
               std::vector< std::vector<size_t> >
               &indValue );
    ///
    /// @brief Derivative function for this reaction class
    ///
    /// @see BaseReaction::derivs(Tissue &T,...)
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
/// @brief Growth via vertex movement radially outwards
///
///  The tissue grows from vertex movement radially outwards,  and also
/// includes moving the vertex defining the 'center' of the cells in the
/// center triangulated mesh. The update is given by
///
/// @f[ \frac{dr}{dt} = p_{0} @f] (if @f$ p_1=0 @f$) or
/// @f[ \frac{dr}{dt} = p_{0} r @f] (if @f$ p_{1}=1 @f$)
///
/// @f$ p_{0} @f$ is the rate (@f$ k_{growth} @f$),
/// @f$ p_{1} @f$ {0,1} is a flag determining which function to be used (@f$ r_{pow} @f$).
///
/// In a model file the reaction is defined as
///
/// @verbatim
/// MoveVertexRadially 2 0
/// p_0 p_1
/// InternalVarStartIndex
/// @endverbatim
///
/// @see MoveVertexRadially
///
class MoveVertexRadiallycenterTriangulation : public BaseReaction {
    
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
    MoveVertexRadiallycenterTriangulation(std::vector<double> &paraValue,
                                          std::vector< std::vector<size_t> >
                                          &indValue );
    ///
    /// @brief Derivative function for this reaction class
    ///
    /// @see BaseReaction::derivs(Tissue &T,...)
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
/// @brief Growth via vertex movement along sphereCylinder
///
/// The tissue grows from vertex movement outwards (from apex) along sphereCylinder.
/// The update can be described by the angular (v) movement
///
/// dv/dt = p_0 (p_1=0) or
/// dv/dt = p_0*r (p_1=1)
///
/// where
///
/// v is the angle from the apex
/// r is the sphere radius
/// p_0 is the rate,
/// p_1 is a flag determining function.
///
/// On the cylinder the vertex is moved downwards (in -z direction).
///
class MoveVertexSphereCylinder : public BaseReaction {
    
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
    MoveVertexSphereCylinder(std::vector<double> &paraValue,
                             std::vector< std::vector<size_t> >
                             &indValue );
    ///
    /// @brief Derivative function for this reaction class
    ///
    /// @see BaseReaction::derivs(Tissue &T,...)
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
/// @brief Updates the water volume variable given osmotic and turgor
/// potentials
///
/// This function uses a constant osmotic potential and calculates the turgor
/// potential to calculate water intake into the cell according to
///
/// @f[ \frac{V_w}{dt} = p_0 A (p_1-p_2T) @f]
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
    WaterVolumeFromTurgor(std::vector<double> &paraValue,
                          std::vector< std::vector<size_t> > &indValue);
    ///
    /// @brief Derivative function for this reaction class
    ///
    /// @see BaseReaction::derivs(Tissue &T,...)
    ///
    void derivs(Tissue &T,
                DataMatrix &cellData,
                DataMatrix &wallData,
                DataMatrix &vertexData,
                DataMatrix &cellDerivs,
                DataMatrix &wallDerivs,
                DataMatrix &vertexDerivs);
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
    DilutionFromVertexDerivs(std::vector<double> &paraValue,
                             std::vector< std::vector<size_t> > &indValue);
    ///
    /// @brief Derivative function for this reaction class
    ///
    /// @see BaseReaction::derivs(Tissue &T,...)
    ///  
    void derivs(Tissue &T,
                DataMatrix &cellData,
                DataMatrix &wallData,
                DataMatrix &vertexData,
                DataMatrix &cellDerivs,
                DataMatrix &wallDerivs,
                DataMatrix &vertexDerivs);
};

#endif
