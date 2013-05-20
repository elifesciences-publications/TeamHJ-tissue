//
// File:    ply_file.cc
// Author:  pkrupinski
// Created: May 06 2013, 3:20 PM
//
#include "ply_file.h"
#include "tissue.h"
// #include <sstream>
#include <stdio.h>
//----------------------------------------------------------------------------
PLY_file::PLY_file() : m_center_triangulation(false)
{
}
//----------------------------------------------------------------------------
PLY_file::PLY_file(const std::string filename) : filename(filename), m_center_triangulation(false)
{
}
//----------------------------------------------------------------------------
PLY_file::~PLY_file()
{
}
//----------------------------------------------------------------------------
PLY_file& PLY_file::open(const std::string filename)
{
    this->filename = filename;
    return *this;
}
//----------------------------------------------------------------------------
void PLY_file::close()
{
    filename = "";
}
//----------------------------------------------------------------------------
void PLY_file::operator<<(Tissue const& t)
{
    write(t);
}
//----------------------------------------------------------------------------
void PLY_file::center_triangulation_output(bool b)
{
  m_center_triangulation = b;
}
//----------------------------------------------------------------------------
void PLY_file::write(Tissue const& t)
{
    std::ofstream m_os(filename.c_str());
    // Print header
    m_os << "ply" << std::endl
    << "format ascii 1.0" << std::endl;
//      << "comment Geometry (simplified) generated by tissue/organism" << std::endl
//      << "comment http://www.thep.lu.se/~henrik/, henrik@thep.lu.se" << std::endl;
    // vertex element definition
    m_os << "element vertex " << t.numVertex() << std::endl
    << "property float x" << std::endl
    << "property float y" << std::endl
    << "property float z" << std::endl;
    m_os << "property uint index" << std::endl;

    // face element definition
    m_os << "element face " << t.numCell() << std::endl
    << "property list uchar uint vertex_index" << std::endl
    << "property uint index" << std::endl
    << "property list uchar float var" << std::endl;
    //include center triangulation definitions if required
    if (m_center_triangulation)
    {
        m_os << "property list uchar float center_vertex" << std::endl;
        m_os << "property list uchar float internal_edges_lengths" << std::endl;
    }
    // edge element definition
    m_os << "element edge " << t.numWall() << std::endl
    << "property uint source" << std::endl
    << "property uint target" << std::endl
    << "property uint index" << std::endl
    << "property list uchar float var" << std::endl;
    m_os << "end_header" << std::endl;

    // Print vertex positions and index
    for (size_t i=0; i<t.numVertex(); ++i) {
        for (size_t d=0; d<t.vertex(i).numPosition(); ++d) {
            m_os << t.vertex(i).position(d) << " ";
        }
        if (t.vertex(0).numPosition()<3) {
            m_os << "0 ";
        }
        m_os << t.vertex(i).index() << std::endl;
    }
    // Print face information, vertex_list, index  and variables
    for (size_t i=0; i<t.numCell(); ++i) {
        m_os << t.cell(i).numVertex() << " ";
        for (size_t k=0; k<t.cell(i).numVertex(); ++k) {
            m_os << t.cell(i).vertex(k)->index() << " ";
        }
        m_os << t.cell(i).index() << " ";
        m_os << t.cell(i).numVariable() << " ";
        for (size_t k=0; k<t.cell(i).numVariable(); ++k) {
            m_os << t.cell(i).variable(k) << " ";
        }
        m_os << std::endl;
        //Print center triangulation variables if required
        if (m_center_triangulation)
        {
            m_os << t.cell(i).numCenterPosition() << " ";
            for (size_t k=0; k<t.cell(i).numCenterPosition(); ++k) {
                m_os << t.cell(i).centerPosition(k) << " ";
            }
            m_os << std::endl;
            
            m_os << t.cell(i).numEdgeLength() << " ";
            for (size_t k=0; k<t.cell(i).numEdgeLength(); ++k) {
                m_os << t.cell(i).edgeLength(k) << " ";
            }
            m_os << std::endl;
        }
    }
    // Print edge information, source, target, index, and variables
    for (size_t i=0; i<t.numWall(); ++i) {
        m_os << t.wall(i).vertex1()->index() << " " << t.wall(i).vertex2()->index() << " ";
        m_os << t.wall(i).index() << " ";
        m_os << t.wall(i).numVariable() << " ";
        for (size_t k=0; k<t.wall(i).numVariable(); ++k) {
            m_os << t.wall(i).variable(k) << " ";
        }
        m_os << std::endl;
    }
    m_os.close();
}
// //----------------------------------------------------------------------------
// PLY_parser::PLY_parser(PLY_file const& f): m_input(f.filename.c_str())m_sep(" "), m_tok(std::string(), m_sep)
// {
//
// }
// //----------------------------------------------------------------------------
// void PLY_parser::parse_header()
// {
//   std::string line;
//   while (!m_input.eof())
//   {
//     getline(m_input, line);
//     if()
//     end_header
//   }
// }
//----------------------------------------------------------------------------
