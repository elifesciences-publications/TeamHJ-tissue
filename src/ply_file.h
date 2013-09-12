//
// File:    ply_file.h
// Author:  pkrupinski
// Created: May 06 2013, 3:20 PM
//
#ifndef _PLY_FILE_H_
#define _PLY_FILE_H_
#include <fstream>
#include <vector>

class Tissue;
// example of usage:
// Tissue T;
// //... initialize T ...
// PLY_file plyfile("output.ply");
// plyfile << T;
//-----------------------------------------------------------------------------
class PLY_output_options
{
public:
  
    PLY_output_options() : m_center_triangulation_output(false), m_bare_geometry_output(false)
    {
    }
    /// @brief Sets output option for center triangulation
    bool& center_triangulation_output() 
    {
        return m_center_triangulation_output;
    }
    bool center_triangulation_output() const 
    {
        return m_center_triangulation_output;
    }
    /// @brief Sets output option for bare geometry
    bool& bare_geometry_output() 
    {
        return m_bare_geometry_output;
    }
    bool bare_geometry_output() const 
    {
        return m_bare_geometry_output;
    }
private:
    /// @brief Flag for center triangulation output
    bool m_center_triangulation_output;
    /// @brief Flag for bare geometry output
    bool m_bare_geometry_output;
};
//-----------------------------------------------------------------------------
/// @brief Handles a ply file for writing a Tissue.
class PLY_file : public PLY_output_options
{
public:
    /// @brief Empty constructor
    PLY_file();
    /// @brief Opens a ply file for writing cell geometry and data
    PLY_file(const std::string filename);
    /// @brief Destructor
    ~PLY_file();
    /// @brief open
    PLY_file& open(const std::string filename);
    /// @brief Close the file
    void close();
    /// @brief Outputs current Tissue
    void operator<<(Tissue const& t);

private:
    /// @brief Name of the ply file
    std::string filename;
    friend class PLY_reader;
};
//-----------------------------------------------------------------------------
class PLY_ostream : public PLY_output_options
{
public:
    PLY_ostream(std::ostream &os);
    /// @brief PLY_ostream output of current Tissue
    PLY_ostream& operator<<(Tissue const& t);
private:
    /// @brief associated ostream
    std::ostream &m_os;
};
//-----------------------------------------------------------------------------
#endif
