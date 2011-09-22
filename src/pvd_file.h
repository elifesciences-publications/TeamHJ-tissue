/**
 * File:    VTUostream.h
 * Author:  pkrupinski
 * Created: September 15, 2011, 3:20 PM
 */
#ifndef _PVD_FILE_H_
#define _PVD_FILE_H_
#include <fstream>
#include <vector>

class Tissue;
class VTUostream;
//-----------------------------------------------------------------------------
// example of usage:
// Tissue T;
// //... initialize T ... 
// PVD_file pvdfile("output.pvd", "VTK_cells.vtu", "VTK_walls.vtu");
// pvdfile << T;
/// @brief Handles a pvd file for writing a Tissue. Uses multiblock dataset for the cell and wall data
class PVD_file
{
public:
  /// @brief Opens a pvd file for writing just a cell geometry and data
  PVD_file(const std::string filename);
  /// @brief Opens a pvd file for writing both cell and wall geometry and data as a multiblock
  PVD_file(const std::string filename, const std::string vtu_filename1, const std::string vtu_filename2);
  /// @brief Destructor
  ~PVD_file();
  /// @brief Outputs current Tissue  
  void operator<<(Tissue const& t);
  /// @brief Outputs current Tissue with a timestamp
  void write(Tissue const& t, double time);

  /// @brief Get the filename of an i'th vtu file associated with this pvd
  std::string const& get_vtu_filename(int i = 0) const { return vtu_filenames[i]; }
  /// @brief Get current counter of time steps already written
  size_t counter() const { return m_counter; }

private:
  typedef void (VTUostream::*CellOutFunPtr)(Tissue const& t);
  /// @brief Append pvd file with new time step vtu members
  void pvdFileWrite(size_t num, double time = -1.0);
  /// @brief Update names of vtu files for a given step
  void vtuNameUpdate(const int counter);

  /// @brief Marks a possition in the pvd file for appending the next step
  std::ios::pos_type mark;
  /// @brief NAme of the pvd file
  std::string filename;
  /// @brief Names of the vtu member files which will be extended with numerical step anotation
  std::vector<std::string> vtu_basenames;
  /// @brief Full current names of the vtu member files
  std::vector<std::string> vtu_filenames;
  /// @brief Counter of the timesteps written so far
  size_t m_counter;
  /// @brief Function pointer to the VTUostream cell output function for the choice of the output mode
  CellOutFunPtr cell_out_fun_ptr;
};

//-----------------------------------------------------------------------------
#endif
