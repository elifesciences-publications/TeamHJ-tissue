/**
 * Filename     : tissue.h
 * Description  : A class describing a two-dimensional tissue of cells
 * Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
 * Created      : April 2006
 * Revision     : $Id$
 */
#ifndef TISSUE_H
#define TISSUE_H

#include <assert.h>
#include <fstream>
#include <iostream>
#include <list>
#include <string>
#include <vector>
#include "baseReaction.h"
#include "baseCompartmentChange.h"
#include "cell.h"
#include "direction.h"
#include "myTypedefs.h"
#include "vertex.h"
#include "wall.h"

///
/// @brief Defines the properties of a two-dimensional cell tissue model
///
/// The Tissue handles the update and information of the cells, walls,
/// vertices and species (molecules). It is
/// the 'top' class for the defined model. It includes reactions
/// updating more than one species simultaneously.
/// Note: currently the tissue is restricted to hold up to 10000 cells, 
/// walls and vertices, such that performance increases as memory is reserved.
///
class Tissue {
  
 private:
  
  std::string id_;           
  
  std::vector<Cell> cell_;
  std::vector<Wall> wall_;
  std::vector<Vertex> vertex_;
  Cell background_;
  Direction direction_;
  std::vector<BaseReaction*> reaction_;
  std::vector<BaseCompartmentChange*> compartmentChange_;
  //DataMatrix tmpCellData_;
  
  std::vector<size_t> directionalWall_;
	
 public:
  
  ///
  /// @brief Empty constructor
  ///
  /// The empty constructor reserves memory for the cells, walls
  /// and vertices, and sets the background 'cell'.
  ///
  Tissue();
  ///
  /// @brief Copy constructor
  ///
  /// Is not yet properly defined and should not be used. It does the
  /// same thing as the empty constructor.
  ///
  Tissue( const Tissue & tissueCopy );
  ///
  /// @brief Constructor from vectors of cells, walls and vertices
  ///
  /// Reserves memory, creates a background 'cell' and copy the vectors of 
  /// cells, walls, and vertices into the tissue.
  ///
  Tissue( const std::vector<Cell> &cellVal,
	  const std::vector<Wall> &wallVal,
	  const std::vector<Vertex> &vertexVal );
  ///
  /// @brief Constructor from reading an init file
  ///
  /// Reserves memory, creates a background 'cell', and then calls the 
  /// readInit() which reads the cell,
  /// wall and vertex information and store it in the tissue.
  ///
  /// @see Tissue::readInit(const char*,int)
  ///
  Tissue( const char *initFile, int verbose=0 );
  ///
  /// @brief Constructor from reading an init file
  ///
  /// Reserves memory, creates a background 'cell', and then calls the 
  /// readInit() which reads the cell,
  /// wall and vertex information and store it in the tissue.
  ///
  /// @see Tissue::readInit(std::string,int)
  ///
  Tissue( std::string initFile, int verbose=0 );
  ///
  /// @brief Constructor from data extracted from an output data file
  ///
  /// Reserves memory, creates a background 'cell', and then use the data
  /// stored in the vectors to create cells, walls and vertices for the 
  /// tissue, including connection (neighborhood) information.
  /// Also checks that the connectivity is proper.
  ///
  /// @see Tissue::checkConnectivity()
  ///
  Tissue( DataMatrix &cellData,
	  DataMatrix &wallData,
	  DataMatrix &vertexData,
	  std::vector< std::vector<size_t> > &cellVertex,
	  std::vector< std::vector<size_t> > &wallVertex,
	  int verbose=0);
  ///
  /// @brief Destructor (currently empty)
  ///
  ~Tissue();
  ///
  /// @brief Opens the file initFile and then calls readInit(std::istream&,int)
  ///
  /// @see Tissue::readInit(std::istream&,int)
  ///
  void readInit(const char *initFile,int verbose=0);
  ///
  /// @brief Opens the file initFile and then calls readInit(std::istream&,int)
  ///
  /// @see Tissue::readInit(std::istream&,int)
  ///
  void readInit(std::string initFile,int verbose=0);
  ///
  /// @brief Reads an initial tissue configuration from an open file
  ///
  /// This function implements the reading of an init file. It first reads
  /// some size information, then the connectivity for each wall. Then follows
  /// the vertex positions, wall length and variables, and cell variables.
  /// An init file has the format:
  ///
  /// @verbatim
  /// N_cell N_wall N_vertex
  /// 
  /// w_i c_i1 c_i2 v_i1 v_i2
  /// ...
  ///
  /// N_vertex dimension
  /// x_i y_i (z_i)
  /// ...
  ///
  /// N_wall N_length N_wallvar
  /// l_i v_i1 [v_i2] ...
  /// ...
  ///
  /// N_cell N_cellvar
  /// v_i1 [v_i2] ...
  /// ...
  ///
  /// @endverbatim
  ///
  /// N_cell, N_wall, N_vertex - number of cells, walls, vertices.
  ///
  /// w_i - wall index (0,1,...), c_i1,c_i2 - the two cells connected to wall i 
  /// (the real indices starts with 0,1,... and if connected to the background
  /// the index is -1),
  /// v_i1, v_i2 - the two vertices connected to wall i (indices 0,1,...).
  ///
  /// N_vertex dimension - number of vertices (again) and dimension. The
  /// dimension has to be two or three.
  /// 
  /// x_i, y_i, (z_i) - the positions for all vertices (z is given if 
  /// dimension=3).
  ///
  /// N_wall, N_length, N_wallvar - Number of walls (again), Number of
  /// (resting) lengths for each wall (in the current implementation
  /// it has to be 1), number of variables for each wall.
  ///
  /// l_i, v_i1, [v_in] - for each wall a length and a variable vector is read.
  ///
  /// N_cell, N_cellvar - Number of cells (again), number of variables for
  /// each cell.
  ///
  /// v_i1, [v_in] - for each cell a variable vector is to be provided.
  ///
  /// Comments can be included in the file by starting a row with # (not yet). 
  /// The connectivity/neighborhood is checked after the tissue is read.
  /// Caveat:
  /// No check on the validity of the initial variable values provided is
  /// is done.
  /// 
  /// @see Tissue::checkConnectivity(int)
  ///
  void readInit(std::istream &IN,int verbose=0);
  ///
  /// @brief Reads an initial tissue configuration in merryproj format
  ///
  /// Reads an initial tissue configuration from the file initFile. It 
  /// assumes that the format is the format provided by the merryproj software.
  /// See the function implementation for the format.
  ///
  void readMerryInit(const char *initFile,int verbose=0);
  ///
  /// @brief Reads an initial tissue configuration from MGX triangular mesh format
  ///
  /// Reads an initial tissue configuration from the file initFile. It 
  /// assumes that the format is the format provided by the MGX software,
  /// by saving the cell format after segmentation on triangulat mesh
  /// (before makeCell).
  /// See the function implementation for the format.
  ///
  void readMGXTriCellInit(const char *initFile,int verbose=0);
  ///
  /// @brief Reads an initial tissue configuration from MGX triangular mesh vtu format
  ///
  /// Reads an initial tissue configuration from the file initFile. It 
  /// assumes that the format is the format provided by the MGX software,
  /// by saving the vtu ascii format after segmentation on triangulated mesh
  /// (before makeCell).
  /// See the function implementation for the format.
  ///
  void readMGXTriVtuInit(const char *initFile,int verbose=0);
  ///
  /// @brief Reads an initial tissue configuration from MGX triangular mesh (mesh) format
  ///
  /// Reads an initial tissue configuration from the file initFile. It 
  /// assumes that the format is the format provided by the MGX software,
  /// by saving the mesh format after segmentation on triangulated mesh
  /// (before makeCell).
  /// See the function implementation for the format.
  ///
  void readMGXTriMeshInit( const char *initFile, int verbose=0);
  ///
  /// @brief Reads an initial tissue configuration from MGX cell mesh format
  ///
  /// Reads an initial tissue configuration from the file initFile. It 
  /// assumes that the format is the format provided by the MGX software,
  /// by exporting to the mesh format after segmentation on triangulated mesh 
  /// and after makeCell.
  /// See the function implementation for the format.
  ///
  void readMGXCellMeshInit(const char *initFile,int verbose=0, int skipCenter=1);
  ///
  /// @brief Opens the file modelFile and then calls readModel(std::ifstream&,int)
  ///
  /// @see Tissue::readModel(std::ifstream&,int)
  ///
  void readModel(const char *modelFile,int verbose=0);
  ///
  /// @brief Opens the file modelFile and then calls readModel(std::ifstream&,int)
  ///
  /// @see Tissue::readModel(std::ifstream&,int)
  ///
  void readModel(std::string modelFile,int verbose=0);
  ///
  /// @brief Reads a tissue model from an open file
  ///
  /// This is the readModel function implementing the actual
  /// reading. First it reads a number of size parameters. Then it 
  /// reads all reactions (rules for updating variables), followed by
  /// all compartmentChanges (rules for adding/dividing and removing cells).
  /// Finally it reads directions (rules for time update and division update
  /// for directional vectors defined for cells). It uses addReaction(),
  /// addCompartmentChange(), and Direction::readDirection() respectively
  /// to read the information for each entity.
  ///
  /// Here follows a detailed discussion on the format of the model file,
  /// divided into reasonable sub-parts. First, the number of entities are read:
  ///
  /// @verbatim
  /// N_reaction N_compChange N_direction
  /// @endverbatim
  ///
  /// N_reaction - number of reactions to be read. N_compChange - number
  /// of compartment changes to be read. N_direction - number of directions 
  /// in the model file (can be 0 or 1 in the current implementation). 
  ///
  /// Then N_reaction reactions are read in the format:
  ///
  /// @verbatim
  /// reaction1 N_p^{r1} N_{il}^{r1} N_{i1}^{r1} ...
  /// p_1^{r1} p_2^{r1} ...
  /// i_{11}^{r1} i_{12}^{r1} ...
  /// i_{21}^{r1} i_{22}^{r1} ...
  /// ...
  /// @endverbatim
  ///
  /// Possible reactions can be found among the classes inheriting
  /// the class BaseReaction.
  ///
  /// Then N_compChange compartment changes are read in the format:
  ///
  /// @verbatim
  /// compartmentChange1 N_p^{r1} N_{il}^{r1} N_{i1}^{r1} ...
  /// p_1^{r1} p_2^{r1} ...
  /// i_{11}^{r1} i_{12}^{r1} ...
  /// i_{21}^{r1} i_{22}^{r1} ...
  /// ...
  /// @endverbatim
  ///
  /// Possible compartment changes can be found among the classes inheriting
  /// the class BaseCompartmentChange.
  ///
  /// Finally, if N_direction=1 rules for defining a cell direction is read:
  ///
  /// @verbatim
  ///  ...
  /// ...
  /// @endverbatim
  ///
  /// Available (time) update rules and division rules for directions can be 
  /// found in classes inheriting BaseDirectionUpdate and 
  /// BaseDirectionDivision.
  ///
  /// Comments can be included in the file by starting a row with # (not 
  /// yet). Caveat:
  /// No complete check on the validity of the provided values are given.
  /// 
  /// @see addReaction(std::ifstream&)
  /// @see addCompartmentChange(std::ifstream&)
  /// @see Direction::readReaction(std::ifstream&)
  /// @see BaseReaction
  /// @see BaseCompartmentChange
  /// @see Direction
  /// @see BaseDirectionUpdate
  /// @see BaseDirectionDivision
  ///
  void readModel(std::ifstream &IN,int verbose=0);
  ///
  /// @brief Reads data from organism sphere format and creates a tissue
  ///
  /// Reads data in an organism sphere format and generates a tissue via createTissueFromSpheres()
  /// The assumed format is:
  /// @verbatim
  /// N_cell N_var
  /// x y [z] r
  /// ...
  /// @endverbatim
  ///
  /// where the data must only contain the positional information at the moment (N_var=3,4).
  /// Also the function assumes the rFac to be set to 1.
  ///
  /// @see createTissueFromSpheres()
  ///
  void readSphereInit( const char *initFile, int verbose=0); 
  ///
  /// @brief Creates a tissue from a vector of sphere data read by readSphereInit()
  ///
  /// @note Caveat: not fully tested. 
  ///
  /// @see readSphereInit(const char*,int)
  ///
  void createTissueFromSpheres(DataMatrix &y,
			       double rFac=1.0, int verbose=0);
  ///
  /// @brief Reads data in voronoi format and creates a tissue
  ///
  /// Reads voronoi output from qhull? and generates a tissue via createTissueFromVoronoi()
  /// The assumed format is:
  /// @verbatim
  /// to be completed...
  /// @endverbatim
  ///
  /// @see createTissueFromVoronoi()
  ///
  void readVoronoiInit( const char *initFile, int verbose=0); 
  ///
  /// @brief Creates a tissue from Voronoi data read by readVoronoiInit()
  ///
  /// @note Caveat: not fully tested.
  ///
  /// @see readVoronoiInit(const char*,int verbose)
  ///
  void createTissueFromVoronoi(DataMatrix &vertexPos,
			       std::vector< std::vector<size_t> > &cellVertex,
			       int verbose=0);
  ///
  /// @brief Sets all cell, wall, and vertex values to their updated values
  ///
  /// This function is used to update the Tissue variables (cell, wall, and vertex variables)
  /// to the state given. This is used to refresh the state since while updating, only
  /// the state matrix variables are updated. At the same time the number of elements in
  /// T is checked.
  ///
  void copyState(DataMatrix &cellData, DataMatrix &wallData, DataMatrix &vertexData);
  ///
  /// @brief The tissue name
  ///
  inline std::string id() const;
  //inline double volume() const;
  //inline size_t mitosisFlag() const;
  ///
  /// @brief Returns the number of cells in the tissue
  ///
  inline size_t numCell() const;
  ///
  /// @brief Returns the number of walls in the tissue
  ///
  inline size_t numWall() const;
  ///
  /// @brief Returns the number of vertices in the tissue
  ///
  inline size_t numVertex() const;
  ///
  /// @brief Returns the number of reactions in the tissue model
  ///
  inline size_t numReaction() const;
  ///
  /// @brief Returns the number of compartmentChanges in the tissue model
  ///
  inline size_t numCompartmentChange() const;
  ///
  /// @brief Returns the number of directionalWalls in the tissue model
  ///
  inline size_t numDirectionalWall() const;
  ///
  /// @brief Returns the directionalWall with index i
  ///
  inline size_t directionalWall(size_t i) const;
  ///
  /// @brief Returns a (const) reference to the tissue cell vector
  ///
  inline const std::vector<Cell> & cell() const;
  ///
  /// @brief Returns a (const) reference to cell i of the tissue
  ///
  inline const Cell & cell(size_t i) const;
  ///
  /// @brief Returns a reference to cell i of the tissue
  ///
  inline Cell & cell(size_t i);
  ///
  /// @brief Returns a pointer to cell i of the tissue
  ///
  inline Cell * cellP(size_t i);
  ///
  /// @brief Adds a cell at the end of the tissue cell vector
  ///
  inline void addCell( Cell val );
  ///
  /// @brief Removes cell i from the tissue and updates the connections
  ///
  inline void removeCell( size_t index );
  ///
  /// @brief Returns a pointer to the background 'cell'
  ///
  /// Every tissue has a background, defined as a cell with index=-1. This is used for knowing
  /// weather an element is at the boundary of the tissue, especially the (outmost) walls that 
  /// need to be connected to two cells use this.
  ///
  inline Cell* background();
  ///
  /// @brief Returns a pointer to the tissue direction
  ///
  inline Direction* direction();
  ///
  /// @brief Returns a reference to the tissue wall vector
  ///
  inline const std::vector<Wall> & wall() const;
  ///
  /// @brief Returns a (const) reference to wall i of the tissue
  ///
  inline const Wall & wall(size_t i) const;
  ///
  /// @brief Returns a reference to wall i of the tissue
  ///
  inline Wall & wall(size_t i);
  ///
  /// @brief Returns a pointer to wall i of the tissue
  ///
  inline Wall * wallP(size_t i);
  ///
  /// @brief Adds a wall at the end of the tissue wall vector
  ///
  inline void addWall( Wall val );
  ///
  /// @brief Removes wall i from the tissue and updates connectivity
  ///
  inline void removeWall( size_t index );
  ///
  /// @brief Returns a reference to the tissue vertex vector
  ///
  inline const std::vector<Vertex> & vertex() const;
  ///
  /// @brief Returns a (const) reference to vertex i of the tissue
  ///
  inline const Vertex & vertex(size_t i) const;
  ///
  /// @brief Returns a reference to vertex i of the tissue
  ///
  inline Vertex & vertex(size_t i);
  ///
  /// @brief Returns a pointer to vertex i of the tissue
  ///
  inline Vertex * vertexP(size_t i);
  ///
  /// @brief returns the number of spatial dimensions of the tissue.
  ///
  /// This function gets the number of spatial dimensions for the
  /// first vertex in the tissue (which is the same for all vertices). 
  /// The number of spatial dimensions for the tissue has to be 2 or 3.
  ///
  inline size_t numDimension();
  ///
  /// @brief Adds a vertex to the tissue at the end of the vector.
  /// 
  inline void addVertex( Vertex val );
  ///
  /// @brief Removes vertex index from the tissue and updates connectivity
  ///
  inline void removeVertex( size_t index );
  ///
  /// @brief Removes a two-vertex and updates connecting walls and cells.
  ///
  /// A two-vertex has only two connected walls (cells). This function removes
  /// such a vertex, and merge the two connecting walls into one, and removes
  /// the connections to the vertex from connecting cells. This is a way to
  /// remove ill-defined vertices from an init file, or after division.
  ///
  void removeTwoVertex( size_t index );
  ///
  /// @brief Returns a pointer to reaction i in the tissue model
  ///
  inline BaseReaction* reaction(size_t i) const;
  ///
  /// @brief Returns a pointer to compartmentChange i in the tissue model
  ///
  inline BaseCompartmentChange* compartmentChange(size_t i) const;
  ///
  /// @brief Sets the name of the tissue.
  ///
  inline void setId(std::string value);
  ///
  /// @brief Resizes the cell vector in the tissue.
  ///
  inline void setNumCell(size_t val);
  ///
  /// @brief Resizes the wall vector in the tissue.
  ///
  inline void setNumWall(size_t val);
  ///
  /// @brief Resizes the vertex vector in the tissue.
  ///
  inline void setNumVertex(size_t val);
  ///
  /// @brief Resizes the directionalWall vector in the tissue.
  ///
  inline void setNumDirectionalWall(size_t val);
  
  //inline const DataMatrix & tmpCellData() const;
  //inline void setTmpCellData(DataMatrix &val);
  ///
  /// @brief Sets all wall length variables from the two vertices positions
  ///
  /// Sets all wall (resting) lengths to the distance between the two
  /// vertices connected to the respective wall.
  ///
  inline void setWallLengthFromVertexPosition();
  ///
  /// @brief Sets directionalWall i of th etissue to the value val
  ///
  inline void setDirectionalWall(size_t i,size_t val);
  ///
  /// @brief Adds a reaction to the list of reactions from an open file 
  ///
  /// Adds a reaction by creating a new reaction via BaseReaction::createReaction
  /// and adds it to the tissue list of reactions.Used when a model is read 
  /// from a file.
  ///
  /// @see BaseReaction::createReaction(std::ifstream&)
  /// @see Tissue::readModel(std::ifstream&,int)
  /// @see BaseReaction
  ///
  int addReaction( std::istream &IN );
  ///
  /// @brief Adds a compartmentChange to the list of from an open file 
  ///
  /// Adds a compartmentChange to the tissue by creating a new 
  /// compartmentChange via BaseCompartmentChange::createCompartmentChange.
  /// A compartmentChange defines a rule for dividing cells or remove 
  /// cells (any rule that change the number of cells). Used when a 
  /// model is read from a file.
  ///
  /// @see BaseCompartmentChange::createReaction(std::ifstream&)
  /// @see Tissue::readModel(std::ifstream&,int)
  /// @see BaseCompartmentChange
  ///
  int addCompartmentChange( std::istream &IN );
  ///
  /// @brief Calculates the derivatives given the state provided 
  ///
  /// This is the main derivatives function used when numerically 
  /// integrating the system. It calls the derivs functions for all 
  /// reactions defined in the tissue model. The current (variable) 
  /// state is given in the
  /// cellData, wallData and vertexData matrices and the derivatives 
  /// are stored in the *Derivs matrices.
  ///
  /// @see BaseReaction::derivs()
  ///
  void derivs( DataMatrix &cellData,
	       DataMatrix &wallData,
	       DataMatrix &vertexData,
	       DataMatrix &cellDeriv,
	       DataMatrix &wallDeriv,
	       DataMatrix &vertexDeriv );
  ///
  /// @brief Initiates the variables via reactions
  ///
  /// This function is called before numerical integration of the
  /// tissue. It loops over reactions and initiate the state variables
  /// where such rules have been defined in the reactions.
  ///
  /// @see BaseReaction::initiate()
  ///
  void initiateReactions(DataMatrix &cellData,
			 DataMatrix &wallData,
			 DataMatrix &vertexData,
			 DataMatrix &cellDeriv,
			 DataMatrix &wallDeriv,
			 DataMatrix &vertexDeriv );
  ///
  /// @brief Takes care of 'discrete' updates during simulations
  ///
  /// Some reactions has discrete updates to be applied in between
  /// steps of the ODE solvers. This functions loop through all reactions
  /// and apply such updates.
  /// 
  /// @see BaseReaction::update()
  ///
  void updateReactions(DataMatrix &cellData,
		       DataMatrix &wallData,
		       DataMatrix &vertexData,
		       double step);
  ///
  /// @brief Initiates direction variables before simulation
  ///
  /// @see Direction::initiate()
  ///
  void initiateDirection(DataMatrix &cellData,
			 DataMatrix &wallData,
			 DataMatrix &vertexData,
			 DataMatrix &cellDerivs,
			 DataMatrix &wallDerivs,
			 DataMatrix &vertexDerivs );
  ///
  /// @brief Updates the directions during simulations according to defined rules.
  ///
  /// @see Direction::update()
  ///
  void updateDirection(double step,							
		       DataMatrix &cellData,
		       DataMatrix &wallData,
		       DataMatrix &vertexData,
		       DataMatrix &cellDerivs,
		       DataMatrix &wallDerivs,
		       DataMatrix &vertexDerivs );
  ///
  /// @brief Updates cell directions at cell divisions
  ///
  /// @see Direction::divide()
  ///
  void updateDirectionDivision(size_t cellI,
			       DataMatrix &cellData,
			       DataMatrix &wallData,
			       DataMatrix &vertexData,
			       DataMatrix &cellDerivs,
			       DataMatrix &wallDerivs,
			       DataMatrix &vertexDerivs);
  ///
  /// @brief Checks for and updates the tissue according to CompartmentChange rules 
  ///
  /// @see BaseCompartmentChange
  ///
  void checkCompartmentChange(DataMatrix &cellData,
			      DataMatrix &wallData,
			      DataMatrix &vertexData,
			      DataMatrix &cellDeriv,
			      DataMatrix &wallDeriv,
			      DataMatrix &vertexDeriv );
  ///
  /// @brief Updates topology and variables for a cell removal
  ///
  void removeCell(size_t cellIndex,
		  DataMatrix &cellData,
		  DataMatrix &wallData,
		  DataMatrix &vertexData,
		  DataMatrix &cellDeriv,
		  DataMatrix &wallDeriv,
		  DataMatrix &vertexDeriv );			
  ///
  /// @brief Calls removeCell(index,...) for all indices given in the vector
  ///
  void removeCells(std::vector<size_t> &cellIndex,
		   DataMatrix &cellData,
		   DataMatrix &wallData,
		   DataMatrix &vertexData,
		   DataMatrix &cellDeriv,
		   DataMatrix &wallDeriv,
		   DataMatrix &vertexDeriv );			
  ///
  /// @brief Calls removeCell(index,...) for all cells that are at 
  /// boundary and outside a radial threshold
  ///
  void removeEpidermalCells(DataMatrix &cellData,
			    DataMatrix &wallData,
			    DataMatrix &vertexData,
			    DataMatrix &cellDeriv,
			    DataMatrix &wallDeriv,
			    DataMatrix &vertexDeriv,
			    double radialThreshold = 0.0,
			    const bool checkBackground = true);
  
  ///
  /// @brief Calls removeCell(index,...) for all cells that are at 
  /// boundary and outside (all vertices) a radial threshold
  ///
  void removeEpidermalCellsMk2(DataMatrix &cellData,
			       DataMatrix &wallData,
			       DataMatrix &vertexData,
			       DataMatrix &cellDeriv,
			       DataMatrix &wallDeriv,
			       DataMatrix &vertexDeriv,
			       double radialThreshold = 0.0);
  
  ///
  /// @brief Calls removeCell(index,...) for all cells that are at the boundary and away from the max
  ///
  void removeEpidermalCellsAtDistance(DataMatrix &cellData,
				      DataMatrix &wallData,
				      DataMatrix &vertexData,
				      DataMatrix &cellDeriv,
				      DataMatrix &wallDeriv,
				      DataMatrix &vertexDeriv,
				      double radialThreshold,double max,
				      size_t direction);
  ///
  /// @brief Updates topology and variables at a cell division
  ///
  void divideCell( Cell *divCell, size_t w1, size_t w2, 
		   std::vector<double> &v1Pos,
		   std::vector<double> &v2Pos,
		   DataMatrix &cellData,
		   DataMatrix &wallData,
		   DataMatrix &vertexData,
		   DataMatrix &cellDeriv,
		   DataMatrix &wallDeriv,
		   DataMatrix &vertexDeriv,
		   std::vector<size_t> &volumeChangeList,
		   double threshold=0.0);

  ///
  /// @brief Updates topology and variables at a cell division assuming trangulation with center
  ///
  /// ...
  /// ...
  ///
  /// @note Requires that the cell is triangulated with the center position stored in the celldata
  ///
  void divideCellCenterTriangulation( Cell *divCell, size_t v1, size_t v2, size_t centerIndex,
				      DataMatrix &cellData,
				      DataMatrix &wallData,
				      DataMatrix &vertexData,
				      DataMatrix &cellDeriv,
				      DataMatrix &wallDeriv,
				      DataMatrix &vertexDeriv,
				      std::vector<size_t> &volumeChangeList );
  
  ///
  /// @brief Sorts cell.wall and cell.vertex vectors to be cyclic 
  ///
  void sortCellWallAndCellVertex(Cell* cell=NULL);
  ///
  /// @brief Recursive algorithm to sort all cells in a tissue
  ///
  void sortCellRecursive( Cell* cell, std::vector<size_t> &sortedFlag, size_t &numSorted);
  ///
  /// @brief Checks all connectivities as well as cell sort for inconsistencies
  ///
  void checkConnectivity(size_t verbose=0);  
  ///
  /// @brief Finds maxima in a variable column for the cells 
  ///
  /// Finds maxima in cells for a specific variable via a local search, and stores
  /// infomration on which cells 'belong' to the different maxima in flag.
  ///
  unsigned int findPeaksGradientAscent( DataMatrix &cellData, 
					size_t col, std::vector<size_t> &cellMax,
					std::vector<size_t> &flag );
  
  // Print functions ----------------------------------------
  ///
  /// @brief print standard tissue init format
  ///  
  void printInit(std::ostream &os=std::cout) const;
  ///
  /// @brief print standard tissue init format using variable values from provided data matrices
  ///  
  void printInit(DataMatrix &cellData,
		 DataMatrix &wallData,
		 DataMatrix &vertexData,
		 std::ostream &os);
  /// 
  /// @brief Prints init in Pawels FEM format
  ///
  /// Prints the current state in Pawels FEM init format using vertex data and
  /// tissue connections.
  ///
  void printInitFem(std::ostream &os) const;
  /// 
  /// @brief Prints init in a triangulated tissue format
  ///
  /// Prints the current state in a triangulated tissue format using cell, wall, and vertex data
  /// and tissue connections. It will do this in a triangulation with a vertex at the center of each
  /// cell.
  /// It triangulates using the cell center calculated by the vertex positions.
  ///
  void printInitTri(std::ostream &os);
  ///
  /// @brief Prints init state in organism init format together with wall areas 
  ///
  /// Creates a file in organism init format together
  /// with cell neighbor information (usually .neigh file for organism) in a single file.
  ///
  /// @see http://dev.thep.lu.se/organism
  /// 
  void printInitOrganism(std::ostream &os);
  void printVertex(std::ostream &os=std::cout);
  void printWall(std::ostream &os=std::cout);
  void printVertexAndCell(std::ostream &os=std::cout);
  void printVertexAndCell(DataMatrix &cellData,
			  DataMatrix &vertexData,
			  std::ostream &os=std::cout);  
  void printVertexAndWall(DataMatrix &wallData,
			  DataMatrix &vertexData,
			  std::ostream &os=std::cout);  
};

inline std::string Tissue::id() const { return id_; }

inline size_t Tissue::numCell() const { return cell_.size(); }

inline size_t Tissue::numWall() const { return wall_.size(); }

inline size_t Tissue::numVertex() const { return vertex_.size(); }

inline size_t Tissue::numReaction() const { return reaction_.size(); }

inline size_t Tissue::numCompartmentChange() const 
{ return compartmentChange_.size(); }

inline size_t Tissue::numDirectionalWall() const 
{ 
  return directionalWall_.size(); 
}

inline size_t Tissue::directionalWall(size_t i) const 
{ 
  return directionalWall_[i];
}

inline const std::vector<Cell> & Tissue::cell() const { return cell_; }

inline const Cell & Tissue::cell(size_t i) const { return cell_[i]; }

inline Cell & Tissue::cell(size_t i) { return cell_[i]; }

inline Cell* Tissue::cellP(size_t i) {return &cell_[i];}

inline void Tissue::addCell( Cell val ) { cell_.push_back(val);}

inline void Tissue::removeCell( size_t index ) {
  assert(index<numCell());
  Cell *cp = &cell_[numCell()-1];
  Cell *cpNew = &cell_[index];
  if( cp != cpNew ) {
    cell_[index] = cell_[numCell()-1];
    cell_[index].setIndex(index);
    //Reconnect wall neighbors
    for( size_t k=0 ; k<cell_[index].numWall() ; ++k )
      if( cell_[index].wall(k)->cell1() == cp )
	cell_[index].wall(k)->setCell1(cpNew);
      else if( cell_[index].wall(k)->cell2() == cp )
	cell_[index].wall(k)->setCell2(cpNew);
    //Reconnect vertex neighbors
    for( size_t k=0 ; k<cell_[index].numVertex() ; ++k )
      for( size_t l=0 ; l<cell_[index].vertex(k)->numCell() ; ++l )
	if( cell_[index].vertex(k)->cell(l) == cp )
	  cell_[index].vertex(k)->setCell(l,cpNew);
  }
  cell_.pop_back();
}

inline Cell* Tissue::background() { return &background_;}

inline Direction* Tissue::direction() { return &direction_;}

inline const std::vector<Wall> & Tissue::wall() const { return wall_; }

inline const Wall & Tissue::wall(size_t i) const { return wall_[i]; }

inline Wall & Tissue::wall(size_t i) { return wall_[i]; }

inline Wall * Tissue::wallP(size_t i) { return &wall_[i]; }

inline void Tissue::addWall( Wall val ) { wall_.push_back(val);}

inline void Tissue::removeWall( size_t index ) {
  assert(index<numWall());
  Wall *wp = &wall_[numWall()-1];
  Wall *wpNew = &wall_[index];
  if( wp != wpNew ) {
    wall_[index] = wall_[numWall()-1];
    wall_[index].setIndex(index);
    //Reconnect cell neighbors
    for( size_t l=0 ; l<wall_[index].cell1()->numWall() ; ++l )
      if( wall_[index].cell1()->wall(l) == wp )
	wall_[index].cell1()->setWall(l,wpNew);
    for( size_t l=0 ; l<wall_[index].cell2()->numWall() ; ++l )
      if( wall_[index].cell2()->wall(l) == wp )
	wall_[index].cell2()->setWall(l,wpNew);
    //Reconnect vertex neighbors
    for( size_t l=0 ; l<wall_[index].vertex1()->numWall() ; ++l )
      if( wall_[index].vertex1()->wall(l) == wp )
	wall_[index].vertex1()->setWall(l,wpNew);
    for( size_t l=0 ; l<wall_[index].vertex2()->numWall() ; ++l )
      if( wall_[index].vertex2()->wall(l) == wp )
	wall_[index].vertex2()->setWall(l,wpNew);	
  }
  wall_.pop_back();
}

inline const std::vector<Vertex> & Tissue::vertex() const { return vertex_; }

inline const Vertex & Tissue::vertex(size_t i) const { return vertex_[i]; }

inline Vertex & Tissue::vertex(size_t i) { return vertex_[i]; }

inline Vertex * Tissue::vertexP(size_t i) { return &vertex_[i]; }

inline size_t Tissue::numDimension() 
{
  return vertex(0).numPosition();
}

inline void Tissue::addVertex( Vertex val ) { vertex_.push_back(val);}

inline void Tissue::removeVertex( size_t index ) {
  assert(index<numVertex());
  Vertex *vp = &vertex_[numVertex()-1];
  Vertex *vpNew = &vertex_[index];
  std::cerr << "Vertex " << vpNew->index() << " removed and "
	    << vp->index() << " moved" << std::endl;
  if( vp != vpNew ) {
    vertex_[index] = vertex_[numVertex()-1];
    vertex_[index].setIndex(index);
    std::cerr << "In cell (" << vertex_[index].numCell() << ") ";
    //Reconnect cell neighbors
    for( size_t k=0 ; k<vertex_[index].numCell() ; ++k ) {
      size_t moved=0;
      for( size_t l=0 ; l<vertex_[index].cell(k)->numVertex() ; ++l )
	if( vertex_[index].cell(k)->vertex(l) == vp ) {
	  vertex_[index].cell(k)->setVertex(l,vpNew);
	  std::cerr << vertex_[index].cell(k)->index() << " ";
	  ++moved;
	}
      if( !moved )
	std::cerr << vertex_[index].cell(k)->index() << "(not) ";
    }
    std::cerr << std::endl;
    //Reconnect wall neighbors
    for( size_t k=0 ; k<vertex_[index].numWall() ; ++k )
      if( vertex_[index].wall(k)->vertex1() == vp )
	vertex_[index].wall(k)->setVertex1(vpNew);
      else if( vertex_[index].wall(k)->vertex2() == vp )
	vertex_[index].wall(k)->setVertex2(vpNew);
  }
  vertex_.pop_back();
}

inline BaseReaction* Tissue::reaction(size_t i) const { 
  assert( i<reaction_.size() );
  return reaction_[i];
}

inline BaseCompartmentChange* Tissue::compartmentChange(size_t i) const { 
  assert( i<compartmentChange_.size() );
  return compartmentChange_[i];
}

inline void Tissue::setId(std::string value) {
  id_ = value;
}

inline void Tissue::setNumCell(size_t val) {
  cell_.resize(val);
}

inline void Tissue::setNumWall(size_t val) {
  wall_.resize(val);
}

inline void Tissue::setNumVertex(size_t val) {
  vertex_.resize(val);
}

inline void Tissue::setNumDirectionalWall(size_t val) 
{
  directionalWall_.resize(val);
}

inline void Tissue::setDirectionalWall(size_t i,size_t val) 
{
  directionalWall_[i]=val;
}

//inline const DataMatrix & Tissue::
//tmpCellData() const 
//{
//	return tmpCellData_;
//}

//inline void Tissue::
//setTmpCellData(DataMatrix &val)
//{
//	tmpCellData_ = val;
//}
#endif
