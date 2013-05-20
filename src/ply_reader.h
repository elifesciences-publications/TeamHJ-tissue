//
// File:    ply_file.h
// Author:  pkrupinski
// Created: May 06 2013, 3:20 PM
//
#ifndef _PLY_READER_H_
#define _PLY_READER_H_
#include <fstream>
#include <vector>
#include "ply_file.h"
#include "ply.h"
#include <tr1/functional>

class Tissue;
class Cell;
class Wall;
class Vertex;
//-----------------------------------------------------------------------------
/// @brief Reades a ply file to a Tissue.
class PLY_reader
{
public:
    /// @brief Constructor
    PLY_reader();
    /// @brief Destructor
    ~PLY_reader();
    /// @brief Reads the ply file to the Tissue
    void read(PLY_file const&f, Tissue &t);

private:
    /// @brief PLY input file stream associated with parser

    const PLY_file *m_ply_file;
    Tissue *m_tissue;

    //array to hold vertex x, y, z values
    std::vector<double> m_position;
    //space to hold index of different parsed quanities
    size_t m_index;
    //space for vertices
    std::vector<Vertex*> m_vertices;
    //space for variables
    std::vector<double> m_variables;
    std::vector<double> m_ceter_vertex;
    std::vector<double> m_edges_lengths;
    
    /// @brief Set connectivity information between tissue Walls and Cells given that their Vertex information is alredy specified
    void set_cell_wall_connectivity(Tissue &t);

    void vertex_begin();
    void vertex_x_callback(ply::float32 x);
    void vertex_y_callback(ply::float32 y);
    void vertex_z_callback(ply::float32 z);
    void vertex_index_callback(ply::uint32 i);
    void vertex_end();

    void face_begin();
    void face_vertex_indices_begin(ply::uint8 size);
    void face_vertex_indices_element(ply::uint32 vertex_index);
    void face_vertex_indices_end();
    void face_index_callback(ply::uint32 i);
    void face_variables_begin(ply::uint8 size);
    void face_variables_element(ply::float32 variable);
    void face_variables_end();
    void face_center_vertex_begin(ply::uint8 size);
    void face_center_vertex_element(ply::float32 variable);
    void face_center_vertex_end();
    void face_edges_lengths_begin(ply::uint8 size);
    void face_edges_lengths_element(ply::float32 variable);
    void face_edges_lengths_end();
    void face_end();

    void edge_begin();
    void edge_source_callback(ply::uint32 i);
    void edge_target_callback(ply::uint32 i);
    void edge_index_callback(ply::uint32 i);
    void edge_variables_begin(ply::uint8 size);
    void edge_variables_element(ply::float32 variable);
    void edge_variables_end();
    void edge_end();



    void info_callback(const std::string& filename, std::size_t line_number, const std::string& message);
    void warning_callback(const std::string& filename, std::size_t line_number, const std::string& message);
    void error_callback(const std::string& filename, std::size_t line_number, const std::string& message);

    std::tr1::tuple<std::tr1::function<void()>, std::tr1::function<void()> > element_definition_callback(const std::string& element_name, std::size_t count);
    template <typename ScalarType> std::tr1::function<void (ScalarType)> scalar_property_definition_callback(const std::string& element_name, const std::string& property_name);
    template <typename SizeType, typename ScalarType> std::tr1::tuple<std::tr1::function<void (SizeType)>, std::tr1::function<void (ScalarType)>, std::tr1::function<void ()> > list_property_definition_callback(const std::string& element_name, const std::string& property_name);
};

//-----------------------------------------------------------------------------
template  <>
inline std::tr1::function <void (ply::float32)> PLY_reader::scalar_property_definition_callback(const std::string& element_name, const std::string& property_name)
{
    if (element_name == "vertex") {
        if (property_name == "x") {
            return std::tr1::bind(&PLY_reader::vertex_x_callback, this, std::tr1::placeholders::_1);
        }
        else if (property_name == "y") {
            return std::tr1::bind(&PLY_reader::vertex_y_callback, this, std::tr1::placeholders::_1);
        }
        else if (property_name == "z") {
            return std::tr1::bind(&PLY_reader::vertex_z_callback, this, std::tr1::placeholders::_1);
        }
        else {
            return 0;
        }
    }
    else {
        return 0;
    }
}
//-----------------------------------------------------------------------------
template  <>
inline std::tr1::function <void (ply::uint32)> PLY_reader::scalar_property_definition_callback(const std::string& element_name, const std::string& property_name)
{
    if (element_name == "vertex") {
        if (property_name == "index") {
            return std::tr1::bind(&PLY_reader::vertex_index_callback, this, std::tr1::placeholders::_1);
        }
        else {
            return 0;
        }
    }
    if (element_name == "face") {
        if (property_name == "index") {
            return std::tr1::bind(&PLY_reader::face_index_callback, this, std::tr1::placeholders::_1);
        }
        else {
            return 0;
        }
    }
    if (element_name == "edge") {
        if (property_name == "index") {
            return std::tr1::bind(&PLY_reader::edge_index_callback, this, std::tr1::placeholders::_1);
        }
        else if (property_name == "source") {
            return std::tr1::bind(&PLY_reader::edge_source_callback, this, std::tr1::placeholders::_1);
        }
        else if (property_name == "target") {
            return std::tr1::bind(&PLY_reader::edge_target_callback, this, std::tr1::placeholders::_1);
        }
        else {
            return 0;
        }
    }
    else {
        return 0;
    }
}
//-----------------------------------------------------------------------------
template <>
inline std::tr1::tuple<std::tr1::function<void (ply::uint8)>, std::tr1::function<void (ply::uint32)>, std::tr1::function<void ()> > PLY_reader::list_property_definition_callback(const std::string& element_name, const std::string& property_name)
{
    if ((element_name == "face") && (property_name == "vertex_index")) {
        return std::tr1::tuple<std::tr1::function<void (ply::uint8)>, std::tr1::function<void (ply::uint32)>, std::tr1::function<void ()> >(
                   std::tr1::bind(&PLY_reader::face_vertex_indices_begin, this, std::tr1::placeholders::_1),
                   std::tr1::bind(&PLY_reader::face_vertex_indices_element, this, std::tr1::placeholders::_1),
                   std::tr1::bind(&PLY_reader::face_vertex_indices_end, this)
               );
    }
    else {
        return std::tr1::tuple<std::tr1::function<void (ply::uint8)>, std::tr1::function<void (ply::uint32)>, std::tr1::function<void ()> >(0, 0, 0);
    }
}
//-----------------------------------------------------------------------------
template <>
inline std::tr1::tuple<std::tr1::function<void (ply::uint8)>, std::tr1::function<void (ply::float32)>, std::tr1::function<void ()> > PLY_reader::list_property_definition_callback(const std::string& element_name, const std::string& property_name)
{
    if (element_name == "face") {
        if (property_name == "var") {
            return std::tr1::tuple<std::tr1::function<void (ply::uint8)>, std::tr1::function<void (ply::float32)>, std::tr1::function<void ()> >(
                       std::tr1::bind(&PLY_reader::face_variables_begin, this, std::tr1::placeholders::_1),
                       std::tr1::bind(&PLY_reader::face_variables_element, this, std::tr1::placeholders::_1),
                       std::tr1::bind(&PLY_reader::face_variables_end, this));
        }
        else if (property_name == "center_vertex") {
            return std::tr1::tuple<std::tr1::function<void (ply::uint8)>, std::tr1::function<void (ply::float32)>, std::tr1::function<void ()> >(
                       std::tr1::bind(&PLY_reader::face_center_vertex_begin, this, std::tr1::placeholders::_1),
                       std::tr1::bind(&PLY_reader::face_center_vertex_element, this, std::tr1::placeholders::_1),
                       std::tr1::bind(&PLY_reader::face_center_vertex_end, this));
        }
        else if (property_name == "internal_edges_lengths") {
            return std::tr1::tuple<std::tr1::function<void (ply::uint8)>, std::tr1::function<void (ply::float32)>, std::tr1::function<void ()> >(
                       std::tr1::bind(&PLY_reader::face_edges_lengths_begin, this, std::tr1::placeholders::_1),
                       std::tr1::bind(&PLY_reader::face_edges_lengths_element, this, std::tr1::placeholders::_1),
                       std::tr1::bind(&PLY_reader::face_edges_lengths_end, this));
        }
    }
    else if ((element_name == "edge") && (property_name == "var")) {
        return std::tr1::tuple<std::tr1::function<void (ply::uint8)>, std::tr1::function<void (ply::float32)>, std::tr1::function<void ()> >(
                   std::tr1::bind(&PLY_reader::edge_variables_begin, this, std::tr1::placeholders::_1),
                   std::tr1::bind(&PLY_reader::edge_variables_element, this, std::tr1::placeholders::_1),
                   std::tr1::bind(&PLY_reader::edge_variables_end, this));
    }
    return std::tr1::tuple<std::tr1::function<void (ply::uint8)>, std::tr1::function<void (ply::float32)>, std::tr1::function<void ()> >(0, 0, 0);
}
//-----------------------------------------------------------------------------
#endif
