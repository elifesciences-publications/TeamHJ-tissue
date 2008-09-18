/**
 * Filename     : compartmentDivision.h
 * Description  : Classes describing compartmentDivision updates
 * Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
 * Created      : July 2006
 * Revision     : $Id:$
 */
#ifndef COMPARTMENTDIVISION_H
#define COMPARTMENTDIVISION_H

#include <cmath>

#include "tissue.h"
#include "baseCompartmentChange.h"

//!Divides a cell when volume above a threshold
/*!Divides a cell when volume above a threshold. New wall is created
  prependicular to the longest cell wall.
 */
class DivisionVolumeViaLongestWall : public BaseCompartmentChange {
  
 public:
  
  DivisionVolumeViaLongestWall(std::vector<double> &paraValue, 
															 std::vector< std::vector<size_t> > 
															 &indValue );
  
  int flag(Tissue *T,size_t i,
					 std::vector< std::vector<double> > &cellData,
					 std::vector< std::vector<double> > &wallData,
					 std::vector< std::vector<double> > &vertexData,
					 std::vector< std::vector<double> > &cellDerivs,
					 std::vector< std::vector<double> > &wallDerivs,
					 std::vector< std::vector<double> > &vertexDerivs );
  void update(Tissue* T,size_t i,
							std::vector< std::vector<double> > &cellData,
							std::vector< std::vector<double> > &wallData,
							std::vector< std::vector<double> > &vertexData,
							std::vector< std::vector<double> > &cellDerivs,
							std::vector< std::vector<double> > &wallDerivs,
							std::vector< std::vector<double> > &vertexDerivs );  
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
					 std::vector< std::vector<double> > &cellData,
					 std::vector< std::vector<double> > &wallData,
					 std::vector< std::vector<double> > &vertexData,
					 std::vector< std::vector<double> > &cellDerivs,
					 std::vector< std::vector<double> > &wallDerivs,
					 std::vector< std::vector<double> > &vertexDerivs );
  void update(Tissue* T,size_t i,
							std::vector< std::vector<double> > &cellData,
							std::vector< std::vector<double> > &wallData,
							std::vector< std::vector<double> > &vertexData,
							std::vector< std::vector<double> > &cellDerivs,
							std::vector< std::vector<double> > &wallDerivs,
							std::vector< std::vector<double> > &vertexDerivs );  
};

//!Divides a cell when volume above a threshold in 3D
/*!Divides a cell when volume above a threshold. Same as
  DivisionVolumeViaLongestWall but used for surfaces in 3D.
*/
class DivisionVolumeViaLongestWall3D : public BaseCompartmentChange {
  
 public:
  
  DivisionVolumeViaLongestWall3D(std::vector<double> &paraValue, 
																 std::vector< std::vector<size_t> > 
																 &indValue );
  
  int flag(Tissue *T,size_t i,
					 std::vector< std::vector<double> > &cellData,
					 std::vector< std::vector<double> > &wallData,
					 std::vector< std::vector<double> > &vertexData,
					 std::vector< std::vector<double> > &cellDerivs,
					 std::vector< std::vector<double> > &wallDerivs,
					 std::vector< std::vector<double> > &vertexDerivs );
  void update(Tissue* T,size_t i,
							std::vector< std::vector<double> > &cellData,
							std::vector< std::vector<double> > &wallData,
							std::vector< std::vector<double> > &vertexData,
							std::vector< std::vector<double> > &cellDerivs,
							std::vector< std::vector<double> > &wallDerivs,
							std::vector< std::vector<double> > &vertexDerivs );  
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
					 std::vector< std::vector<double> > &cellData,
					 std::vector< std::vector<double> > &wallData,
					 std::vector< std::vector<double> > &vertexData,
					 std::vector< std::vector<double> > &cellDerivs,
					 std::vector< std::vector<double> > &wallDerivs,
					 std::vector< std::vector<double> > &vertexDerivs );
  void update(Tissue* T,size_t i,
							std::vector< std::vector<double> > &cellData,
							std::vector< std::vector<double> > &wallData,
							std::vector< std::vector<double> > &vertexData,
							std::vector< std::vector<double> > &cellDerivs,
							std::vector< std::vector<double> > &wallDerivs,
							std::vector< std::vector<double> > &vertexDerivs );  
};

//!Divides a cell when volume above a threshold
/*!Divides a cell when volume above a threshold. New wall is created
  prependicular to maximal strain rate.
 */
class DivisionVolumeViaStrain : public BaseCompartmentChange {
  
 public:
  
  DivisionVolumeViaStrain(std::vector<double> &paraValue, 
													std::vector< std::vector<size_t> > 
													&indValue );
  
  int flag(Tissue *T,size_t i,
					 std::vector< std::vector<double> > &cellData,
					 std::vector< std::vector<double> > &wallData,
					 std::vector< std::vector<double> > &vertexData,
					 std::vector< std::vector<double> > &cellDerivs,
					 std::vector< std::vector<double> > &wallDerivs,
					 std::vector< std::vector<double> > &vertexDerivs );
  void update(Tissue* T,size_t i,
							std::vector< std::vector<double> > &cellData,
							std::vector< std::vector<double> > &wallData,
							std::vector< std::vector<double> > &vertexData,
							std::vector< std::vector<double> > &cellDerivs,
							std::vector< std::vector<double> > &wallDerivs,
							std::vector< std::vector<double> > &vertexDerivs );  
};

//!Divides a cell when volume above a threshold
/*!Divides a cell when volume above a threshold. New wall is created
  prependicular to direction given as cell variable.
 */
class DivisionVolumeViaDirection : public BaseCompartmentChange {
  
 public:
  
  DivisionVolumeViaDirection(std::vector<double> &paraValue, 
														 std::vector< std::vector<size_t> > 
														 &indValue );
  
  int flag(Tissue *T,size_t i,
					 std::vector< std::vector<double> > &cellData,
					 std::vector< std::vector<double> > &wallData,
					 std::vector< std::vector<double> > &vertexData,
					 std::vector< std::vector<double> > &cellDerivs,
					 std::vector< std::vector<double> > &wallDerivs,
					 std::vector< std::vector<double> > &vertexDerivs );
  void update(Tissue* T,size_t i,
							std::vector< std::vector<double> > &cellData,
							std::vector< std::vector<double> > &wallData,
							std::vector< std::vector<double> > &vertexData,
							std::vector< std::vector<double> > &cellDerivs,
							std::vector< std::vector<double> > &wallDerivs,
							std::vector< std::vector<double> > &vertexDerivs );  
};

//!Divides a cell when volume above a threshold
/*!Divides a cell when volume above a threshold. New wall is created
  in a random direction through center of mass.
 */
class DivisionVolumeRandomDirection : public BaseCompartmentChange {
  
 public:
  
  DivisionVolumeRandomDirection(std::vector<double> &paraValue, 
														 std::vector< std::vector<size_t> > 
														 &indValue );
  
  int flag(Tissue *T,size_t i,
					 std::vector< std::vector<double> > &cellData,
					 std::vector< std::vector<double> > &wallData,
					 std::vector< std::vector<double> > &vertexData,
					 std::vector< std::vector<double> > &cellDerivs,
					 std::vector< std::vector<double> > &wallDerivs,
					 std::vector< std::vector<double> > &vertexDerivs );
  void update(Tissue* T,size_t i,
							std::vector< std::vector<double> > &cellData,
							std::vector< std::vector<double> > &wallData,
							std::vector< std::vector<double> > &vertexData,
							std::vector< std::vector<double> > &cellDerivs,
							std::vector< std::vector<double> > &wallDerivs,
							std::vector< std::vector<double> > &vertexDerivs );  
};


class DivisionForceDirection : public BaseCompartmentChange
{
 public:
	DivisionForceDirection(std::vector<double> &paraValue, 
					   std::vector< std::vector<size_t> > &indValue);
	
	int flag(Tissue *T, size_t i,
		    std::vector< std::vector<double> > &cellData,
		    std::vector< std::vector<double> > &wallData,
		    std::vector< std::vector<double> > &vertexData,
		    std::vector< std::vector<double> > &cellDerivs,
		    std::vector< std::vector<double> > &wallDerivs,
		    std::vector< std::vector<double> > &vertexDerivs);
	void update(Tissue* T, size_t i,
			  std::vector< std::vector<double> > &cellData,
			  std::vector< std::vector<double> > &wallData,
			  std::vector< std::vector<double> > &vertexData,
			  std::vector< std::vector<double> > &cellDerivs,
			  std::vector< std::vector<double> > &wallDerivs,
			  std::vector< std::vector<double> > &vertexDerivs);  
};

//!Divides a cell when volume above a threshold
/*!Divides a cell when volume above a threshold. New wall is created at shortest
	path that divides the volume in equal parts. 
 */
class DivisionVolumeViaShortestPath : public BaseCompartmentChange {
  
 public:
  
  DivisionVolumeViaShortestPath(std::vector<double> &paraValue, 
																std::vector< std::vector<size_t> > 
																&indValue );
  
  int flag(Tissue *T,size_t i,
					 std::vector< std::vector<double> > &cellData,
					 std::vector< std::vector<double> > &wallData,
					 std::vector< std::vector<double> > &vertexData,
					 std::vector< std::vector<double> > &cellDerivs,
					 std::vector< std::vector<double> > &wallDerivs,
					 std::vector< std::vector<double> > &vertexDerivs );
  void update(Tissue* T,size_t i,
							std::vector< std::vector<double> > &cellData,
							std::vector< std::vector<double> > &wallData,
							std::vector< std::vector<double> > &vertexData,
							std::vector< std::vector<double> > &cellDerivs,
							std::vector< std::vector<double> > &wallDerivs,
							std::vector< std::vector<double> > &vertexDerivs );  
};

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
		 std::vector< std::vector<double> > &cellData,
		 std::vector< std::vector<double> > &wallData,
		 std::vector< std::vector<double> > &vertexData,
		 std::vector< std::vector<double> > &cellDerivs,
		    std::vector< std::vector<double> > &wallDerivs,
		 std::vector< std::vector<double> > &vertexDerivs);
	void update(Tissue* T, size_t i,
		    std::vector< std::vector<double> > &cellData,
		    std::vector< std::vector<double> > &wallData,
		    std::vector< std::vector<double> > &vertexData,
		    std::vector< std::vector<double> > &cellDerivs,
		    std::vector< std::vector<double> > &wallDerivs,
		    std::vector< std::vector<double> > &vertexDerivs);  
	
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
		std::vector< std::vector<double> > &cellData,
		std::vector< std::vector<double> > &wallData,
		std::vector< std::vector<double> > &vertexData,
		std::vector< std::vector<double> > &cellDerivs,
		std::vector< std::vector<double> > &wallDerivs,
		std::vector< std::vector<double> > &vertexDerivs);
	
	void update(Tissue* T, size_t i,
		std::vector< std::vector<double> > &cellData,
		std::vector< std::vector<double> > &wallData,
		std::vector< std::vector<double> > &vertexData,
		std::vector< std::vector<double> > &cellDerivs,
		std::vector< std::vector<double> > &wallDerivs,
		std::vector< std::vector<double> > &vertexDerivs);  
	
	double astar(double sigma, double A, double B);
	double f(double a, double sigma, double A, double B);
	int sign(double a);
};

class DivisionRandom : public BaseCompartmentChange
{
public:
	DivisionRandom(std::vector<double> &paraValue, std::vector< std::vector<size_t> > &indValue);
  
	int flag(Tissue *T, size_t i,
		std::vector< std::vector<double> > &cellData,
		std::vector< std::vector<double> > &wallData,
		std::vector< std::vector<double> > &vertexData,
		std::vector< std::vector<double> > &cellDerivs,
		std::vector< std::vector<double> > &wallDerivs,
		std::vector< std::vector<double> > &vertexDerivs);
	
	void update(Tissue* T, size_t i,
		std::vector< std::vector<double> > &cellData,
		std::vector< std::vector<double> > &wallData,
		std::vector< std::vector<double> > &vertexData,
		std::vector< std::vector<double> > &cellDerivs,
		std::vector< std::vector<double> > &wallDerivs,
		std::vector< std::vector<double> > &vertexDerivs);  

	/** Returns an integer between 0 and n - 1. */
	int random(int n);
};

class DivisionVolumeRandomDirectionGiantCells : public BaseCompartmentChange
{
public:
	
	DivisionVolumeRandomDirectionGiantCells(std::vector<double> &paraValue, 
		std::vector< std::vector<size_t> > &indValue);
  
	int flag(Tissue *T, size_t i,
		std::vector< std::vector<double> > &cellData,
		std::vector< std::vector<double> > &wallData,
		std::vector< std::vector<double> > &vertexData,
		std::vector< std::vector<double> > &cellDerivs,
		std::vector< std::vector<double> > &wallDerivs,
		std::vector< std::vector<double> > &vertexDerivs);

	void update(Tissue* T, size_t i,
		std::vector< std::vector<double> > &cellData,
		std::vector< std::vector<double> > &wallData,
		std::vector< std::vector<double> > &vertexData,
		std::vector< std::vector<double> > &cellDerivs,
		std::vector< std::vector<double> > &wallDerivs,
		std::vector< std::vector<double> > &vertexDerivs);  
};

class DivisionMainAxis : public BaseCompartmentChange
{
public:
	
	DivisionMainAxis(std::vector<double> &paraValue, 
		std::vector< std::vector<size_t> > &indValue);
  
	int flag(Tissue *T, size_t i,
		std::vector< std::vector<double> > &cellData,
		std::vector< std::vector<double> > &wallData,
		std::vector< std::vector<double> > &vertexData,
		std::vector< std::vector<double> > &cellDerivs,
		std::vector< std::vector<double> > &wallDerivs,
		std::vector< std::vector<double> > &vertexDerivs);
	
	void update(Tissue* T, size_t i,
		std::vector< std::vector<double> > &cellData,
		std::vector< std::vector<double> > &wallData,
		std::vector< std::vector<double> > &vertexData,
		std::vector< std::vector<double> > &cellDerivs,
		std::vector< std::vector<double> > &wallDerivs,
		std::vector< std::vector<double> > &vertexDerivs);  

	std::vector<double> getMainAxis(Cell &cell, std::vector< std::vector<double> > &vertexData);

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

