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

#include "tissue.h"
#include "baseCompartmentChange.h"




///
/// @brief Divides a cell when volume above a threshold, with new wall perpendicular to the longest wall segment.
/// Divides a cell when volume above a threshold. New wall is created
/// prependicular to the longest cell wall. In a model file it is defined as
/// 
/// @verbatim
/// DivisionVolumeViaLongestWall 3 1 [1]
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
class DivisionVolumeViaLongestWall : public BaseCompartmentChange {
  
 public:
  
  DivisionVolumeViaLongestWall(std::vector<double> &paraValue, 
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
/// DivisionVolumeViaLongestWallCenterTriangulation 3 2 [1] 2
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
class DivisionVolumeViaLongestWallCenterTriangulation : public BaseCompartmentChange {
  
 public:
  
  DivisionVolumeViaLongestWallCenterTriangulation(std::vector<double> &paraValue, 
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
/// DivisionVolumeViaLongestWall3DCenterTriangulation 3 2 [1] 2
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
class DivisionVolumeViaLongestWall3DCenterTriangulation : public BaseCompartmentChange {
  
 public:
  
  DivisionVolumeViaLongestWall3DCenterTriangulation(std::vector<double> &paraValue, 
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
class DivisionVolumeViaLongestWallSpatial : public BaseCompartmentChange {
  
 private:
  
  double sMax_;
  
 public:
  
  DivisionVolumeViaLongestWallSpatial(std::vector<double> &paraValue, 
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
/// DivisionVolumeViaLongestWall but used for surfaces in 3D.

class DivisionVolumeViaLongestWall3D : public BaseCompartmentChange {
  
 public:
  
  DivisionVolumeViaLongestWall3D(std::vector<double> &paraValue, 
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
/// DivisionVolumeViaLongestWall but used for surfaces in 3D.
///
class DivisionVolumeViaLongestWall3DSpatial : public BaseCompartmentChange {
  
 private:
  
  double sMax_;
  
 public:
  
  DivisionVolumeViaLongestWall3DSpatial(std::vector<double> &paraValue, 
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

class DivisionVolumeViaStrain : public BaseCompartmentChange {
  
 public:
  
  DivisionVolumeViaStrain(std::vector<double> &paraValue, 
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

class DivisionVolumeViaDirection : public BaseCompartmentChange {
  
 public:
  
  DivisionVolumeViaDirection(std::vector<double> &paraValue, 
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

class DivisionVolumeRandomDirection : public BaseCompartmentChange {
  
 public:
  
  DivisionVolumeRandomDirection(std::vector<double> &paraValue, 
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
class DivisionVolumeRandomDirectionCenterTriangulation : public BaseCompartmentChange {
  
 public:
  
  DivisionVolumeRandomDirectionCenterTriangulation(std::vector<double> &paraValue, 
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

class DivisionForceDirection : public BaseCompartmentChange
{
 public:
  DivisionForceDirection(std::vector<double> &paraValue, 
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

class DivisionVolumeViaShortestPath : public BaseCompartmentChange {
  
 public:
  
  DivisionVolumeViaShortestPath(std::vector<double> &paraValue, 
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
/// @verbatim
///
/// DivisionShortestPath 4 2 0/1 1 
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
/// DivisionShortestPath 6 3 0/1 1 2 
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

class DivisionShortestPath : public BaseCompartmentChange
{
 public:
  struct Candidate {
    double distance;
    size_t wall1;
    size_t wall2;
    double px, py;
    double qx, qy;
  };
  
  DivisionShortestPath(std::vector<double> &paraValue, 
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
  
  std::vector<DivisionShortestPath::Candidate> 
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

class DivisionShortestPathGiantCells : public BaseCompartmentChange
{
 public:
  struct Candidate {
    double distance;
    size_t wall1;
    size_t wall2;
    double px, py;
    double qx, qy;
  };
  
  DivisionShortestPathGiantCells(std::vector<double> &paraValue, 
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
  
  std::vector<DivisionShortestPathGiantCells::Candidate> 
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

class DivisionRandom : public BaseCompartmentChange
{
 public:
  DivisionRandom(std::vector<double> &paraValue, std::vector< std::vector<size_t> > &indValue);
  
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

class DivisionVolumeRandomDirectionGiantCells : public BaseCompartmentChange
{
 public:
  
  DivisionVolumeRandomDirectionGiantCells(std::vector<double> &paraValue, 
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

class DivisionMainAxis : public BaseCompartmentChange
{
 public:
  
  DivisionMainAxis(std::vector<double> &paraValue, 
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

#endif

