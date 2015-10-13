//
// Filename     : compartmentDivision.h
// Description  : Classes describing compartmentDivision updates
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : July 2006
// Revision     : $Id:$
//
#ifndef COMPARTMENTDIVISION_H
#define COMPARTMENTDIVISION_H

#include <cmath>
#include <vector>

#include "tissue.h"
#include "baseCompartmentChange.h"

namespace Division {
  ///
  /// @brief Divides a cell when volume above a threshold, with new wall perpendicular to the longest wall segment.
  /// Divides a cell when volume above a threshold. New wall is created
  /// prependicular to the longest cell wall. In a model file it is defined as
  /// 
  /// @verbatim
  /// Division::VolumeViaLongestWall 3 1 [1]
  /// V_th //L^{wall}_{frac} L^{wall}_{threshold}
  /// I1
  /// @endverbatim
  ///
  /// where @f$V_{th}@f$ is the cell volume threshold, @f$L^{wall}_{frac}@f$ is the resting length 
  /// of the new wall (1.0 sets it to the distance between the vertices, and @f$L^{wall}_{threshold}@f$ 
  /// is the smallest (relative) length of the new subwalls (i.e. if closer than this to an existing vertex
  /// it will be moved to this distance from the old vertex).
  ///
  /// The list of indices given are for those variables that need to be updated due to the division,
  /// e.g. concentrations do not, the volume itself (if stored) needs to as well as molecular numbers.
  /// 
  class VolumeViaLongestWall : public BaseCompartmentChange {
    
  public:
    
    VolumeViaLongestWall(std::vector<double> &paraValue, 
				 std::vector< std::vector<size_t> > 
				 &indValue );
    int flag(Tissue *T,size_t i,
	     DataMatrix &cellData,
	     DataMatrix &wallData,
	     DataMatrix &vertexData,
	     DataMatrix &cellDerivs,
	     DataMatrix &wallDerivs,
	     DataMatrix &vertexDerivs );
    void update(Tissue* T,size_t i,
		DataMatrix &cellData,
		DataMatrix &wallData,
		DataMatrix &vertexData,
		DataMatrix &cellDerivs,
		DataMatrix &wallDerivs,
		DataMatrix &vertexDerivs );  
  };




  ///
  /// @brief Divides a cell when volume above a threshold, with new wall perpendicular to the longest wall segment.
  ///
  /// Divides a cell when volume above a threshold. New wall is created
  /// prependicular to the longest cell wall. In a model file it is defined as
  /// @verbatim
  /// Division::VolumeViaLongestWallCenterTriangulation 3 2 [1] 2
  /// V_th L^{wall}_{frac} L^{wall}_{threshold}
  /// I1
  /// com index, restinglengthIndex
  /// @endverbatim
  ///
  /// where @f$V_{th}@f$ is the cell volume threshold, @f$L^{wall}_{frac}@f$ is the resting length 
  /// of the new wall (1.0 sets it to the distance between the vertices, and @f$L^{wall}_{threshold}@f$ 
  /// is the smallest (relative) length of the new subwalls (i.e. if closer than this to an existing vertex
  /// it will be moved to this distance from the old vertex).
  ///
  /// The list of indices given are for those variables that need to be updated due to the division,
  /// e.g. concentrations do not, the volume itself (if stored) needs to as well as molecular numbers.
  /// 
  class VolumeViaLongestWallCenterTriangulation : public BaseCompartmentChange {
    
  public:
    
    VolumeViaLongestWallCenterTriangulation(std::vector<double> &paraValue, 
						    std::vector< std::vector<size_t> > 
						    &indValue );
    int flag(Tissue *T,size_t i,
	     DataMatrix &cellData,
	     DataMatrix &wallData,
	     DataMatrix &vertexData,
	     DataMatrix &cellDerivs,
	     DataMatrix &wallDerivs,
	     DataMatrix &vertexDerivs );
    void update(Tissue* T,size_t i,
		DataMatrix &cellData,
		DataMatrix &wallData,
		DataMatrix &vertexData,
		DataMatrix &cellDerivs,
		DataMatrix &wallDerivs,
		DataMatrix &vertexDerivs );  
  };
  
  
  
  
  ///
  /// @brief Divides a cell when volume above a threshold, with new wall perpendicular to the longest wall segment.
  ///
  /// Divides a cell when volume above a threshold in 3D with centerTriangulation. New wall is created
  /// prependicular to the longest cell wall. In a model file it is defined as
  /// @verbatim
  /// Division::VolumeViaLongestWall3DCenterTriangulation 3 2 [1] 2
  /// V_th L^{wall}_{frac} L^{wall}_{threshold}
  /// I1
  /// com index, restinglengthIndex
  /// @endverbatim
  ///
  /// where @f$V_{th}@f$ is the cell volume threshold, @f$L^{wall}_{frac}@f$ is the resting length 
  /// of the new wall (1.0 sets it to the distance between the vertices, and @f$L^{wall}_{threshold}@f$ 
  /// is the smallest (relative) length of the new subwalls (i.e. if closer than this to an existing vertex
  /// it will be moved to this distance from the old vertex).
  ///
  /// The list of indices given are for those variables that need to be updated due to the division,
  /// e.g. concentrations do not, the volume itself (if stored) needs to as well as molecular numbers.
  /// 
  class VolumeViaLongestWall3DCenterTriangulation : public BaseCompartmentChange {
    
  public:
    
    VolumeViaLongestWall3DCenterTriangulation(std::vector<double> &paraValue, 
						      std::vector< std::vector<size_t> > 
						      &indValue );
    int flag(Tissue *T,size_t i,
	     DataMatrix &cellData,
	     DataMatrix &wallData,
	     DataMatrix &vertexData,
	     DataMatrix &cellDerivs,
	     DataMatrix &wallDerivs,
	     DataMatrix &vertexDerivs );
    void update(Tissue* T,size_t i,
		DataMatrix &cellData,
		DataMatrix &wallData,
		DataMatrix &vertexData,
		DataMatrix &cellDerivs,
		DataMatrix &wallDerivs,
		DataMatrix &vertexDerivs );  
  };
  
  ///
  /// @brief Divides a cell when volume above a threshold, with new wall perpendicular to the longest wall segment.
  ///
  /// Creates a branch as
  /// @verbatim
  /// Branching 3 1 [1]
  /// V_th L^{wall}_{frac} L^{wall}_{threshold}
  /// I1
  /// @endverbatim
  ///
  /// where @f$V_{th}@f$ is the cell volume threshold, @f$L^{wall}_{frac}@f$ is the resting length 
  /// of the new wall (1.0 sets it to the distance between the vertices, and @f$L^{wall}_{threshold}@f$ 
  /// is the smallest (relative) length of the new subwalls (i.e. if closer than this to an existing vertex
  /// it will be moved to this distance from the old vertex).
  ///
  /// The list of indices given are for those variables that need to be updated due to the division,
  /// e.g. concentrations do not, the volume itself (if stored) needs to as well as molecular numbers.
  /// 
  class Branching : public BaseCompartmentChange {
    
  public:
    
    Branching(std::vector<double> &paraValue, 
	      std::vector< std::vector<size_t> > 
	      &indValue );
    int flag(Tissue *T,size_t i,
	     DataMatrix &cellData,
	     DataMatrix &wallData,
	     DataMatrix &vertexData,
	     DataMatrix &cellDerivs,
	     DataMatrix &wallDerivs,
	     DataMatrix &vertexDerivs );
    void update(Tissue* T,size_t i,
		DataMatrix &cellData,
		DataMatrix &wallData,
		DataMatrix &vertexData,
		DataMatrix &cellDerivs,
		DataMatrix &wallDerivs,
		DataMatrix &vertexDerivs );  
  };
  
  
  ///
  /// @brief Divides a cell when volume above a threshold and cell close enough to 'apex'
  ///
  /// Divides a cell when volume above a threshold and cell within a distance
  /// from max. New wall is created prependicular to the longest cell wall.
  ///
  class VolumeViaLongestWallSpatial : public BaseCompartmentChange {
    
  private:
    
    double sMax_;
    
  public:
    
    VolumeViaLongestWallSpatial(std::vector<double> &paraValue, 
					std::vector< std::vector<size_t> > 
					&indValue );
    
    int flag(Tissue *T,size_t i,
	     DataMatrix &cellData,
	     DataMatrix &wallData,
	     DataMatrix &vertexData,
	     DataMatrix &cellDerivs,
	     DataMatrix &wallDerivs,
	     DataMatrix &vertexDerivs );
    void update(Tissue* T,size_t i,
		DataMatrix &cellData,
		DataMatrix &wallData,
		DataMatrix &vertexData,
		DataMatrix &cellDerivs,
		DataMatrix &wallDerivs,
		DataMatrix &vertexDerivs );  
  };
  
  /// @brief Divides a cell when volume above a threshold in 3D
  /// Divides a cell when volume above a threshold. Same as
  /// Division::VolumeViaLongestWall but used for surfaces in 3D.
  
  class VolumeViaLongestWall3D : public BaseCompartmentChange {
    
  public:
    
    VolumeViaLongestWall3D(std::vector<double> &paraValue, 
				   std::vector< std::vector<size_t> > 
				   &indValue );
    
    int flag(Tissue *T,size_t i,
	     DataMatrix &cellData,
	     DataMatrix &wallData,
	     DataMatrix &vertexData,
	     DataMatrix &cellDerivs,
	     DataMatrix &wallDerivs,
	     DataMatrix &vertexDerivs );
    void update(Tissue* T,size_t i,
		DataMatrix &cellData,
		DataMatrix &wallData,
		DataMatrix &vertexData,
		DataMatrix &cellDerivs,
		DataMatrix &wallDerivs,
		DataMatrix &vertexDerivs );  
  };
  
  ///
  /// @brief Divides a cell when volume above a threshold in 3D and distance to apex smaller than th
  ///
  /// Divides a cell when volume above a threshold. Same as
  /// Division::VolumeViaLongestWall but used for surfaces in 3D.
  ///
  class VolumeViaLongestWall3DSpatial : public BaseCompartmentChange {
    
  private:
    
    double sMax_;
    
  public:
    
    VolumeViaLongestWall3DSpatial(std::vector<double> &paraValue, 
					  std::vector< std::vector<size_t> > 
					  &indValue );
    
    int flag(Tissue *T,size_t i,
	     DataMatrix &cellData,
	     DataMatrix &wallData,
	     DataMatrix &vertexData,
	     DataMatrix &cellDerivs,
	     DataMatrix &wallDerivs,
	     DataMatrix &vertexDerivs );
    void update(Tissue* T,size_t i,
		DataMatrix &cellData,
		DataMatrix &wallData,
		DataMatrix &vertexData,
		DataMatrix &cellDerivs,
		DataMatrix &wallDerivs,
		DataMatrix &vertexDerivs );  
  };
  
  /// @brief Divides a cell when volume above a threshold
  /// Divides a cell when volume above a threshold. New wall is created
  ///  prependicular to maximal strain rate.
  
  class VolumeViaStrain : public BaseCompartmentChange {
    
  public:
    
    VolumeViaStrain(std::vector<double> &paraValue, 
			    std::vector< std::vector<size_t> > 
			    &indValue );
    
    int flag(Tissue *T,size_t i,
	     DataMatrix &cellData,
	     DataMatrix &wallData,
	     DataMatrix &vertexData,
	     DataMatrix &cellDerivs,
	     DataMatrix &wallDerivs,
	     DataMatrix &vertexDerivs );
    void update(Tissue* T,size_t i,
		DataMatrix &cellData,
		DataMatrix &wallData,
		DataMatrix &vertexData,
		DataMatrix &cellDerivs,
		DataMatrix &wallDerivs,
		DataMatrix &vertexDerivs );  
  };
  
  /// @brief Divides a cell when volume above a threshold
  /// Divides a cell when volume above a threshold. New wall is created
  /// prependicular to direction given as cell variable.
  
  class VolumeViaDirection : public BaseCompartmentChange {
    
  public:
    
    VolumeViaDirection(std::vector<double> &paraValue, 
			       std::vector< std::vector<size_t> > 
			       &indValue );
    
    int flag(Tissue *T,size_t i,
	     DataMatrix &cellData,
	     DataMatrix &wallData,
	     DataMatrix &vertexData,
	     DataMatrix &cellDerivs,
	     DataMatrix &wallDerivs,
	     DataMatrix &vertexDerivs );
    void update(Tissue* T,size_t i,
		DataMatrix &cellData,
		DataMatrix &wallData,
		DataMatrix &vertexData,
		DataMatrix &cellDerivs,
		DataMatrix &wallDerivs,
		DataMatrix &vertexDerivs );  
  };
  
  /// @brief Divides a cell when volume above a threshold
  /// Divides a cell when volume above a threshold. New wall is created
  ///  in a random direction through center of mass.
  
  class VolumeRandomDirection : public BaseCompartmentChange {
    
  public:
    
    VolumeRandomDirection(std::vector<double> &paraValue, 
				  std::vector< std::vector<size_t> > 
				  &indValue );
    
    int flag(Tissue *T,size_t i,
	     DataMatrix &cellData,
	     DataMatrix &wallData,
	     DataMatrix &vertexData,
	     DataMatrix &cellDerivs,
	     DataMatrix &wallDerivs,
	     DataMatrix &vertexDerivs );
    void update(Tissue* T,size_t i,
		DataMatrix &cellData,
		DataMatrix &wallData,
		DataMatrix &vertexData,
		DataMatrix &cellDerivs,
		DataMatrix &wallDerivs,
		DataMatrix &vertexDerivs );  
  };
  
  ///
  /// @brief  Divides a cell when volume above a threshold in a random direction
  ///
  /// Divides a cell when volume above a threshold. New wall is created
  /// from a random (wall central) vertex and the second vertex is
  /// chosen such that the two (daughter) areas are as equal as
  /// possible.
  ///
  /// @note This class relies on that a central vertex is defined.
  ///
  class VolumeRandomDirectionCenterTriangulation : public BaseCompartmentChange {
    
  public:
    
    VolumeRandomDirectionCenterTriangulation(std::vector<double> &paraValue, 
						     std::vector< std::vector<size_t> > 
						     &indValue );
    
    int flag(Tissue *T,size_t i,
	     DataMatrix &cellData,
	     DataMatrix &wallData,
	     DataMatrix &vertexData,
	     DataMatrix &cellDerivs,
	     DataMatrix &wallDerivs,
	     DataMatrix &vertexDerivs );
    void update(Tissue* T,size_t i,
		DataMatrix &cellData,
		DataMatrix &wallData,
		DataMatrix &vertexData,
		DataMatrix &cellDerivs,
		DataMatrix &wallDerivs,
		DataMatrix &vertexDerivs );  
  };
  
  class ForceDirection : public BaseCompartmentChange
  {
  public:
    ForceDirection(std::vector<double> &paraValue, 
			   std::vector< std::vector<size_t> > &indValue);
    
    int flag(Tissue *T, size_t i,
	     DataMatrix &cellData,
	     DataMatrix &wallData,
	     DataMatrix &vertexData,
	     DataMatrix &cellDerivs,
	     DataMatrix &wallDerivs,
	     DataMatrix &vertexDerivs);
    void update(Tissue* T, size_t i,
		DataMatrix &cellData,
		DataMatrix &wallData,
		DataMatrix &vertexData,
		DataMatrix &cellDerivs,
		DataMatrix &wallDerivs,
		DataMatrix &vertexDerivs);  
  };
  
  /// @brief Divides a cell when volume above a threshold
  /// Divides a cell when volume above a threshold. New wall is created at shortest
  /// path that divides the volume in equal parts. 
  
  class VolumeViaShortestPath : public BaseCompartmentChange {
    
  public:
    
    VolumeViaShortestPath(std::vector<double> &paraValue, 
				  std::vector< std::vector<size_t> > 
				  &indValue );
    
    int flag(Tissue *T,size_t i,
	     DataMatrix &cellData,
	     DataMatrix &wallData,
	     DataMatrix &vertexData,
	     DataMatrix &cellDerivs,
	     DataMatrix &wallDerivs,
	     DataMatrix &vertexDerivs );
    void update(Tissue* T,size_t i,
		DataMatrix &cellData,
		DataMatrix &wallData,
		DataMatrix &vertexData,
		DataMatrix &cellDerivs,
		DataMatrix &wallDerivs,
		DataMatrix &vertexDerivs );  
  };
  
  ///
  /// @brief Divides a cell when volume above a threshold, with New wall created at shortest
  ///  path that divides the volume (not!) in equal parts. Using centerTriangulation and doubleLength 
  /// formats are optional and can be done by setting the coresponding flags. 
  /// 
  ///
  /// @verbatim
  ///
  /// Division::ShortestPath 4 2 0/1 1 
  /// V_{threshold} 
  /// L^{wall}_{frac} (relative of new wall)
  /// L^{wall}_{threshold} (disallowed closeness)
  /// centerCom flag(0:random, 1:COM)
  ///
  /// I1 (optional volume related index to be updated)
  ///
  /// cell time index(optional)
  ///
  /// @endverbatim
  ///
  /// or
  ///
  /// @verbatim
  ///
  /// Division::ShortestPath 6 3 0/1 1 2 
  /// V_{threshold} 
  /// L^{wall}_{frac} (relative of new wall)
  /// L^{wall}_{threshold} (disallowed closeness)
  /// centerCom flag(0:random, 1:COM)
  /// centerTriangulation flag (0/1)
  /// double length flag (0/1)
  ///
  /// I1 (optional volume related index to be updated)
  ///
  /// cell time index(optional)
  ///
  /// com index 
  /// restinglengthIndex
  ///
  /// @endverbatim
  
  class ShortestPath : public BaseCompartmentChange
  {
  public:
    struct Candidate {
      double distance;
      size_t wall1;
      size_t wall2;
      double px, py;
      double qx, qy;
    };
    
    ShortestPath(std::vector<double> &paraValue, 
			 std::vector< std::vector<size_t> > &indValue);
    
    int flag(Tissue *T, size_t i,
	     DataMatrix &cellData,
	     DataMatrix &wallData,
	     DataMatrix &vertexData,
	     DataMatrix &cellDerivs,
	     DataMatrix &wallDerivs,
	     DataMatrix &vertexDerivs);
    void update(Tissue* T, size_t i,
		DataMatrix &cellData,
		DataMatrix &wallData,
		DataMatrix &vertexData,
		DataMatrix &cellDerivs,
		DataMatrix &wallDerivs,
		DataMatrix &vertexDerivs);  
    
    std::vector<ShortestPath::Candidate> 
      getCandidates(Tissue* T, size_t i,
		    DataMatrix &cellData,
		    DataMatrix &wallData,
		    DataMatrix &vertexData,
		    DataMatrix &cellDerivs,
		    DataMatrix &wallDerivs,
		    DataMatrix &vertexDerivs);
    
    double astar(double sigma, double A, double B);
    double f(double a, double sigma, double A, double B);
    int sign(double a);
  };

  ///
  /// @brief Divides a cell when volume above a sizer/timer/adder (STA) rule, with new wall created at shortest
  ///  path through center-of-mass (or with random fluctuation outside of com).
  ///
  /// This division rule use a sizer/timer/adder rule for deciding when to divide, by dividing when
  ///
  /// @f[ V_{Division} > (2-p_{0}) + p_{0} V_{Birth} @f]
  ///
  /// i.e. relating the division volume to the size of the birth (size after previous division) of the cell.
  /// The parameter f sets the rule to:
  /// p0 = 2 -> timer
  /// p0 = 0 -> sizer (scaled to division at volume 2)
  /// p0 = 1 -> adder
  ///
  /// At division the division plane is chosen as the shortest path through the center of mass of the cell.
  /// The com can be replaced with a random point close to the com (setting flag in p3).
  /// There is a restriction to not be too close to a vertex (p2). 
  /// A parameter will set the resting length of the new wall (p1).
  ///
  /// In addition, formats are optional and can be done by setting the coresponding flags. 
  ///
  /// @verbatim
  /// Division::ShortestPath 4 2 0/1 1 
  /// cellDivisionRule (p0)
  /// L^{wall}_{frac} (relative of new wall)
  /// L^{wall}_{threshold} (disallowed closeness)
  /// centerCom flag(0:random, 1:COM)
  ///
  /// I1 (optional volume related index to be updated)
  ///
  /// cell time index(optional)
  ///
  /// @endverbatim
  ///
  /// or
  ///
  /// @verbatim
  ///
  /// Division::ShortestPath 6 3 0/1 1 2 
  /// V_{threshold} 
  /// L^{wall}_{frac} (relative of new wall)
  /// L^{wall}_{threshold} (disallowed closeness)
  /// centerCom flag(0:random, 1:COM)
  /// centerTriangulation flag (0/1)
  /// double length flag (0/1)
  ///
  /// I1 (optional volume related index to be updated)
  ///
  /// cell time index(optional)
  ///
  /// com index 
  /// restinglengthIndex
  ///
  /// @endverbatim  
  class STAViaShortestPath : public BaseCompartmentChange
  {
  public:
    struct Candidate {
      double distance;
      size_t wall1;
      size_t wall2;
      double px, py;
      double qx, qy;
    };
    
    STAViaShortestPath(std::vector<double> &paraValue,
		       std::vector< std::vector<size_t> > &indValue);
    
    int flag(Tissue *T, size_t i,
	     DataMatrix &cellData,
	     DataMatrix &wallData,
	     DataMatrix &vertexData,
	     DataMatrix &cellDerivs,
	     DataMatrix &wallDerivs,
	     DataMatrix &vertexDerivs);
    void update(Tissue* T, size_t i,
		DataMatrix &cellData,
		DataMatrix &wallData,
		DataMatrix &vertexData,
		DataMatrix &cellDerivs,
		DataMatrix &wallDerivs,
		DataMatrix &vertexDerivs);  
    
    std::vector<STAViaShortestPath::Candidate> 
      getCandidates(Tissue* T, size_t i,
		    DataMatrix &cellData,
		    DataMatrix &wallData,
		    DataMatrix &vertexData,
		    DataMatrix &cellDerivs,
		    DataMatrix &wallDerivs,
		    DataMatrix &vertexDerivs);
    
    double astar(double sigma, double A, double B);
    double f(double a, double sigma, double A, double B);
    int sign(double a);
  };


  /// @brief UNDER CONSTRUCTION, DO NOT USE YET!!!  Divides a cell when a certain condition (flag=1) is fulfilled, with New wall created at shortest
  ///  path that divides the volume (not!) in equal parts. Using centerTriangulation and doubleLength 
  ///  formats are optional and can be done by setting the coresponding flags. 
  /// 
  ///
  /// @verbatim
  ///
  /// Division::FlagResetShortestPath 4 2 0/1 1 
  /// V_{threshold} 
  /// L^{wall}_{frac} (relative of new wall)
  /// L^{wall}_{threshold} (disallowed closeness)
  /// centerCom flag(0:random, 1:COM)
  ///
  /// I1 (optional volume related index to be updated)
  ///
  /// cell time index(optional)
  ///
  /// @endverbatim
  ///
  /// or
  ///
  /// @verbatim
  ///
  /// Division::FlagResetShortestPath 6 4 0/1 0/1 2 1 
  /// V_{threshold} 
  /// L^{wall}_{frac} (relative of new wall)
  /// L^{wall}_{threshold} (disallowed closeness)
  /// centerCom flag(0:random, 1:COM)
  /// centerTriangulation flag (0/1)
  /// double length flag (0/1)
  /// I1 (optional volume related index to be updated)
  ///
  /// cell time index(optional)
  ///
  /// com index 
  /// restinglengthIndex
  /// flag index, indicating the variable-flag that is responsible for division (when flag=1 a cell divides, when flag=0 a cell does not divide)

  ///
  /// @endverbatim
  class FlagResetShortestPath : public BaseCompartmentChange
  {
  public:
    struct Candidate {
      double distance;
      size_t wall1;
      size_t wall2;
      double px, py;
      double qx, qy;
    };
    
    FlagResetShortestPath(std::vector<double> &paraValue, 
       std::vector< std::vector<size_t> > &indValue);
    
    int flag(Tissue *T, size_t i,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs);
    void update(Tissue* T, size_t i,
    DataMatrix &cellData,
    DataMatrix &wallData,
    DataMatrix &vertexData,
    DataMatrix &cellDerivs,
    DataMatrix &wallDerivs,
    DataMatrix &vertexDerivs);  
    
    std::vector<FlagResetShortestPath::Candidate> 
      getCandidates(Tissue* T, size_t i,
        DataMatrix &cellData,
        DataMatrix &wallData,
        DataMatrix &vertexData,
        DataMatrix &cellDerivs,
        DataMatrix &wallDerivs,
        DataMatrix &vertexDerivs);
    
    double astar(double sigma, double A, double B);
    double f(double a, double sigma, double A, double B);
    int sign(double a);
  };

  class ShortestPathGiantCells : public BaseCompartmentChange
  {
  public:
    struct Candidate {
      double distance;
      size_t wall1;
      size_t wall2;
      double px, py;
      double qx, qy;
    };
    
    ShortestPathGiantCells(std::vector<double> &paraValue, 
				   std::vector< std::vector<size_t> > &indValue);
    
    int flag(Tissue *T, size_t i,
	     DataMatrix &cellData,
	     DataMatrix &wallData,
	     DataMatrix &vertexData,
	     DataMatrix &cellDerivs,
	     DataMatrix &wallDerivs,
	     DataMatrix &vertexDerivs);
    
    void update(Tissue* T, size_t i,
		DataMatrix &cellData,
		DataMatrix &wallData,
		DataMatrix &vertexData,
		DataMatrix &cellDerivs,
		DataMatrix &wallDerivs,
		DataMatrix &vertexDerivs);  
    
    std::vector<ShortestPathGiantCells::Candidate> 
      getCandidates(Tissue* T, size_t i,
		    DataMatrix &cellData,
		    DataMatrix &wallData,
		    DataMatrix &vertexData,
		    DataMatrix &cellDerivs,
		    DataMatrix &wallDerivs,
		    DataMatrix &vertexDerivs);
    
    double astar(double sigma, double A, double B);
    double f(double a, double sigma, double A, double B);
    int sign(double a);
  };

  
  class Random : public BaseCompartmentChange
  {
  public:
    Random(std::vector<double> &paraValue, std::vector< std::vector<size_t> > &indValue);
    
    int flag(Tissue *T, size_t i,
	     DataMatrix &cellData,
	     DataMatrix &wallData,
	     DataMatrix &vertexData,
	     DataMatrix &cellDerivs,
	     DataMatrix &wallDerivs,
	     DataMatrix &vertexDerivs);
    
    void update(Tissue* T, size_t i,
		DataMatrix &cellData,
		DataMatrix &wallData,
		DataMatrix &vertexData,
		DataMatrix &cellDerivs,
		DataMatrix &wallDerivs,
		DataMatrix &vertexDerivs);  
    
    // Returns an integer between 0 and n - 1. 
    int random(int n);
  };
  
  class VolumeRandomDirectionGiantCells : public BaseCompartmentChange
  {
  public:
    
    VolumeRandomDirectionGiantCells(std::vector<double> &paraValue, 
					    std::vector< std::vector<size_t> > &indValue);
    
    int flag(Tissue *T, size_t i,
	     DataMatrix &cellData,
	     DataMatrix &wallData,
	     DataMatrix &vertexData,
	     DataMatrix &cellDerivs,
	     DataMatrix &wallDerivs,
	     DataMatrix &vertexDerivs);
    
    void update(Tissue* T, size_t i,
		DataMatrix &cellData,
		DataMatrix &wallData,
		DataMatrix &vertexData,
		DataMatrix &cellDerivs,
		DataMatrix &wallDerivs,
		DataMatrix &vertexDerivs);  
  };
  
  class MainAxis : public BaseCompartmentChange
  {
  public:
    
    MainAxis(std::vector<double> &paraValue, 
		     std::vector< std::vector<size_t> > &indValue);
    
    int flag(Tissue *T, size_t i,
	     DataMatrix &cellData,
	     DataMatrix &wallData,
	     DataMatrix &vertexData,
	     DataMatrix &cellDerivs,
	     DataMatrix &wallDerivs,
	     DataMatrix &vertexDerivs);
    
    void update(Tissue* T, size_t i,
		DataMatrix &cellData,
		DataMatrix &wallData,
		DataMatrix &vertexData,
		DataMatrix &cellDerivs,
		DataMatrix &wallDerivs,
		DataMatrix &vertexDerivs);  
    
    std::vector<double> getMainAxis(Cell &cell, DataMatrix &vertexData);
    
  private:
    struct Candidate
    {
      double s;
      size_t index;
      std::vector<double> p;
    };
    
    class CompareCandidate
    {
    public:
      bool operator()(const Candidate &a, const Candidate &b)
      {
	return std::abs(a.s) > std::abs(b.s);
      }
    };
  };


  /// @brief UNDER CONSTRUCTION, DO NOT USE YET!!!

    class FlagResetViaLongestWall : public BaseCompartmentChange {
    
  public:
    
    FlagResetViaLongestWall(std::vector<double> &paraValue, 
         std::vector< std::vector<size_t> > 
         &indValue );
    int flag(Tissue *T,size_t i,
       DataMatrix &cellData,
       DataMatrix &wallData,
       DataMatrix &vertexData,
       DataMatrix &cellDerivs,
       DataMatrix &wallDerivs,
       DataMatrix &vertexDerivs );
    void update(Tissue* T,size_t i,
    DataMatrix &cellData,
    DataMatrix &wallData,
    DataMatrix &vertexData,
    DataMatrix &cellDerivs,
    DataMatrix &wallDerivs,
    DataMatrix &vertexDerivs );  
  };



} //end namespace Division
#endif
  
  
