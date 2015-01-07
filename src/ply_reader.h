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
#include "tissue.h"

#if C11NSPACE == std
#include <functional>
#elif C11NSPACE == boost
#include <boost/functional.hpp>
#endif

class Tissue;
class Cell;
class Wall;
class Vertex;
size_t const INVALID_SIZE = -1;
//----------------------------------------------------------------------------
struct Edge_element
{
    typedef std::pair<size_t, size_t> Edge;

    Edge_element ( const Edge& e ) : edge ( e ), index ( INVALID_SIZE )
    {

    }

    bool operator< ( const Edge_element& ee ) const
    {
        return edge < ee.edge;
    }

    Edge edge;
    mutable size_t index;

};
//----------------------------------------------------------------------------
class PLY_reader_options
{
public:
    PLY_reader_options() : m_index_base ( 0 ) {}
    int& index_base()
    {
        return m_index_base;
    }
    int index_base() const
    {
        return m_index_base;
    }
private:
    int m_index_base;
};
//-----------------------------------------------------------------------------
/// @brief Reades a ply file to a Tissue.
class PLY_reader: public PLY_reader_options
{
public:
    /// @brief Constructor
    PLY_reader();
    /// @brief Destructor
    ~PLY_reader();
    /// @brief Reads the ply file to the Tissue
    void read ( PLY_file const&f, Tissue &t );

private:
    /// @brief PLY input file stream associated with parser

    const PLY_file *m_ply_file;
    Tissue *m_tissue;

    //array to hold vertex x, y, z values
    std::vector<double> m_position;
    //space to hold index of different parsed quanities
    size_t m_index;
    //vertex counter
    size_t m_v_counter;
    //face counter
    size_t m_f_counter;
    //edge counter
    size_t m_e_counter;
    //space for vertices
    std::vector<Vertex*> m_vertices;
    //space for variables
    std::vector<double> m_variables;
    std::vector<double> m_center_vertex;
    std::vector<double> m_edges_lengths;

    /// @brief Set connectivity information between tissue Walls and Cells given that their Vertex information is alredy specified
    void set_cell_wall_connectivity ( Tissue &t );
    /// @brief adds an edge information to the tissue given the cyclic ordering of vertices in each cell
    void infer_walls_from_cells ( Tissue &t );

    void vertex_begin();
    void vertex_x_callback ( ply::float32 x );
    void vertex_y_callback ( ply::float32 y );
    void vertex_z_callback ( ply::float32 z );
    void vertex_index_callback ( ply::uint32 i );
    void vertex_end();

    void face_begin();
    void face_vertex_indices_end();
    void face_index_callback ( ply::uint32 i );
    void face_variables_begin ( ply::uint8 size );
    void face_variables_element ( ply::float32 variable );
    void face_variables_end();
    void face_center_vertex_begin ( ply::uint8 size );
    void face_center_vertex_element ( ply::float32 variable );
    void face_center_vertex_end();
    void face_edges_lengths_begin ( ply::uint8 size );
    void face_edges_lengths_element ( ply::float32 variable );
    void face_edges_lengths_end();
    void face_end();

    void edge_begin();
    void edge_source_callback ( ply::uint32 i );
    void edge_target_callback ( ply::uint32 i );
    void edge_index_callback ( ply::uint32 i );
    void edge_variables_begin ( ply::uint8 size );
    void edge_variables_element ( ply::float32 variable );
    void edge_variables_end();
    void edge_end();

    void info_callback ( const std::string& filename, std::size_t line_number, const std::string& message );
    void warning_callback ( const std::string& filename, std::size_t line_number, const std::string& message );
    void error_callback ( const std::string& filename, std::size_t line_number, const std::string& message );

    C11NSPACE::tuple<C11NSPACE::function<void() >, C11NSPACE::function<void() > > element_definition_callback ( const std::string& element_name, std::size_t count );
    template <typename ScalarType> C11NSPACE::function<void ( ScalarType ) > scalar_property_definition_callback ( const std::string& element_name, const std::string& property_name );
    template <typename SizeType, typename ScalarType> C11NSPACE::tuple<C11NSPACE::function<void ( SizeType ) >, C11NSPACE::function<void ( ScalarType ) >, C11NSPACE::function<void () > > list_property_definition_callback ( const std::string& element_name, const std::string& property_name );

    template <typename SizeType> void face_vertex_indices_begin ( SizeType size );
    template <typename ScalarType>void face_vertex_indices_element ( ScalarType vertex_index );
};
//----------------------------------------------------------------------------
template <typename SizeType>
void PLY_reader::face_vertex_indices_begin ( SizeType size )
{
    m_vertices.reserve ( size );
}
//----------------------------------------------------------------------------
template <typename ScalarType>
void PLY_reader::face_vertex_indices_element ( ScalarType vertex_index )
{
    assert ( vertex_index >= 0 && vertex_index < ScalarType(m_tissue->numVertex()) );
    m_vertices.push_back ( m_tissue->vertexP ( vertex_index ) );
}
//-----------------------------------------------------------------------------
template  <>
inline C11NSPACE::function <void ( ply::float32 ) > PLY_reader::scalar_property_definition_callback ( const std::string& element_name, const std::string& property_name )
{
    if ( element_name == "vertex" )
    {
        if ( property_name == "x" )
        {
            return C11NSPACE::bind ( &PLY_reader::vertex_x_callback, this, C11NSPACE::placeholders::_1 );
        }
        else if ( property_name == "y" )
        {
            return C11NSPACE::bind ( &PLY_reader::vertex_y_callback, this, C11NSPACE::placeholders::_1 );
        }
        else if ( property_name == "z" )
        {
            return C11NSPACE::bind ( &PLY_reader::vertex_z_callback, this, C11NSPACE::placeholders::_1 );
        }
        else
        {
            return 0;
        }
    }
    else
    {
        return 0;
    }
}
//-----------------------------------------------------------------------------
template  <>
inline C11NSPACE::function <void ( ply::uint32 ) > PLY_reader::scalar_property_definition_callback ( const std::string& element_name, const std::string& property_name )
{
    if ( element_name == "vertex" )
    {
        if ( property_name == "index" )
        {
            return C11NSPACE::bind ( &PLY_reader::vertex_index_callback, this, C11NSPACE::placeholders::_1 );
        }
        else
        {
            return 0;
        }
    }
    if ( element_name == "face" )
    {
        if ( property_name == "index" )
        {
            return C11NSPACE::bind ( &PLY_reader::face_index_callback, this, C11NSPACE::placeholders::_1 );
        }
        else
        {
            return 0;
        }
    }
    if ( element_name == "edge" )
    {
        if ( property_name == "index" )
        {
            return C11NSPACE::bind ( &PLY_reader::edge_index_callback, this, C11NSPACE::placeholders::_1 );
        }
        else if ( property_name == "source" )
        {
            return C11NSPACE::bind ( &PLY_reader::edge_source_callback, this, C11NSPACE::placeholders::_1 );
        }
        else if ( property_name == "target" )
        {
            return C11NSPACE::bind ( &PLY_reader::edge_target_callback, this, C11NSPACE::placeholders::_1 );
        }
        else
        {
            return 0;
        }
    }
    else
    {
        return 0;
    }
}
//-----------------------------------------------------------------------------
template <typename SizeType, typename ScalarType>
inline C11NSPACE::tuple<C11NSPACE::function<void ( SizeType ) >, C11NSPACE::function<void ( ScalarType ) >, C11NSPACE::function<void () > > PLY_reader::list_property_definition_callback ( const std::string& element_name, const std::string& property_name )
{
    if ( ( element_name == "face" ) && ( property_name == "vertex_index" || property_name == "vertex_indices" ) )
    {
        return C11NSPACE::tuple<C11NSPACE::function<void ( SizeType ) >, C11NSPACE::function<void ( ScalarType ) >, C11NSPACE::function<void () > > (
                   C11NSPACE::bind ( &PLY_reader::face_vertex_indices_begin<SizeType>, this, C11NSPACE::placeholders::_1 ),
                   C11NSPACE::bind ( &PLY_reader::face_vertex_indices_element<ScalarType>, this, C11NSPACE::placeholders::_1 ),
                   C11NSPACE::bind ( &PLY_reader::face_vertex_indices_end, this )
               );
    }
    else
    {
        return C11NSPACE::tuple<C11NSPACE::function<void ( SizeType ) >, C11NSPACE::function<void ( ply::uint32 ) >, C11NSPACE::function<void () > > ( 0, 0, 0 );
    }
}
//-----------------------------------------------------------------------------
template <>
inline C11NSPACE::tuple<C11NSPACE::function<void ( ply::uint8 ) >, C11NSPACE::function<void ( ply::float32 ) >, C11NSPACE::function<void () > > PLY_reader::list_property_definition_callback ( const std::string& element_name, const std::string& property_name )
{
    if ( element_name == "face" )
    {
        if ( property_name == "var" )
        {
            return C11NSPACE::tuple<C11NSPACE::function<void ( ply::uint8 ) >, C11NSPACE::function<void ( ply::float32 ) >, C11NSPACE::function<void () > > (
                       C11NSPACE::bind ( &PLY_reader::face_variables_begin, this, C11NSPACE::placeholders::_1 ),
                       C11NSPACE::bind ( &PLY_reader::face_variables_element, this, C11NSPACE::placeholders::_1 ),
                       C11NSPACE::bind ( &PLY_reader::face_variables_end, this ) );
        }
        else if ( property_name == "center_vertex" )
        {
            return C11NSPACE::tuple<C11NSPACE::function<void ( ply::uint8 ) >, C11NSPACE::function<void ( ply::float32 ) >, C11NSPACE::function<void () > > (
                       C11NSPACE::bind ( &PLY_reader::face_center_vertex_begin, this, C11NSPACE::placeholders::_1 ),
                       C11NSPACE::bind ( &PLY_reader::face_center_vertex_element, this, C11NSPACE::placeholders::_1 ),
                       C11NSPACE::bind ( &PLY_reader::face_center_vertex_end, this ) );
        }
        else if ( property_name == "internal_edges_lengths" )
        {
            return C11NSPACE::tuple<C11NSPACE::function<void ( ply::uint8 ) >, C11NSPACE::function<void ( ply::float32 ) >, C11NSPACE::function<void () > > (
                       C11NSPACE::bind ( &PLY_reader::face_edges_lengths_begin, this, C11NSPACE::placeholders::_1 ),
                       C11NSPACE::bind ( &PLY_reader::face_edges_lengths_element, this, C11NSPACE::placeholders::_1 ),
                       C11NSPACE::bind ( &PLY_reader::face_edges_lengths_end, this ) );
        }
    }
    else if ( ( element_name == "edge" ) && ( property_name == "var" ) )
    {
        return C11NSPACE::tuple<C11NSPACE::function<void ( ply::uint8 ) >, C11NSPACE::function<void ( ply::float32 ) >, C11NSPACE::function<void () > > (
                   C11NSPACE::bind ( &PLY_reader::edge_variables_begin, this, C11NSPACE::placeholders::_1 ),
                   C11NSPACE::bind ( &PLY_reader::edge_variables_element, this, C11NSPACE::placeholders::_1 ),
                   C11NSPACE::bind ( &PLY_reader::edge_variables_end, this ) );
    }
    return C11NSPACE::tuple<C11NSPACE::function<void ( ply::uint8 ) >, C11NSPACE::function<void ( ply::float32 ) >, C11NSPACE::function<void () > > ( 0, 0, 0 );
}
//-----------------------------------------------------------------------------
#endif
