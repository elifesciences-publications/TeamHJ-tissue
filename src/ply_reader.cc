//
// File:    ply_file.cc
// Author:  pkrupinski
// Created: May 06 2013, 3:20 PM
//
#include "ply_reader.h"
#include "tissue.h"
#include "cell.h"
#include "wall.h"

using namespace std::tr1::placeholders;
//----------------------------------------------------------------------------
PLY_reader::PLY_reader() : m_ply_file(NULL), m_tissue(NULL), m_position(3)
{
}
//----------------------------------------------------------------------------
PLY_reader::~PLY_reader()
{
}
//----------------------------------------------------------------------------
void PLY_reader::read(PLY_file const&f, Tissue &t)
{
    ply::ply_parser parser;

    m_ply_file = &f;
    m_tissue = &t;
    std::string filename  = m_ply_file->filename;

    parser.info_callback(std::tr1::bind(&PLY_reader::info_callback, this, std::tr1::ref(filename), _1, _2));
    parser.warning_callback(std::tr1::bind(&PLY_reader::warning_callback, this, std::tr1::ref(filename), _1, _2));
    parser.error_callback(std::tr1::bind(&PLY_reader::error_callback, this, std::tr1::ref(filename), _1, _2));

    parser.element_definition_callback(std::tr1::bind(&PLY_reader::element_definition_callback, this, _1, _2));

    ply::ply_parser::scalar_property_definition_callbacks_type scalar_property_definition_callbacks;
    ply::at<ply::float32>(scalar_property_definition_callbacks) = std::tr1::bind(&PLY_reader::scalar_property_definition_callback<ply::float32>, this, _1, _2);
    ply::at<ply::uint32>(scalar_property_definition_callbacks) = std::tr1::bind(&PLY_reader::scalar_property_definition_callback<ply::uint32>, this, _1, _2);
    parser.scalar_property_definition_callbacks(scalar_property_definition_callbacks);

    ply::ply_parser::list_property_definition_callbacks_type list_property_definition_callbacks;
    ply::at<ply::uint8, ply::uint32>(list_property_definition_callbacks) = std::tr1::bind(&PLY_reader::list_property_definition_callback<ply::uint8, ply::uint32>, this, _1, _2);
    ply::at<ply::uint8, ply::float32>(list_property_definition_callbacks) = std::tr1::bind(&PLY_reader::list_property_definition_callback<ply::uint8, ply::float32>, this, _1, _2);
    parser.list_property_definition_callbacks(list_property_definition_callbacks);

    parser.parse(filename);
    
    set_cell_wall_connectivity(t);
}
//----------------------------------------------------------------------------
std::tr1::tuple<std::tr1::function<void()>, std::tr1::function<void()> > PLY_reader::element_definition_callback(const std::string& element_name, std::size_t count)
{
    if (element_name == "vertex") {
        m_tissue->setNumVertex(count);
        return std::tr1::tuple<std::tr1::function<void()>, std::tr1::function<void()> >(
                   std::tr1::bind(&PLY_reader::vertex_begin, this),
                   std::tr1::bind(&PLY_reader::vertex_end, this)
               );
    }
    else if (element_name == "face") {
        m_tissue->setNumCell(count);
        return std::tr1::tuple<std::tr1::function<void()>, std::tr1::function<void()> >(
                   std::tr1::bind(&PLY_reader::face_begin, this),
                   std::tr1::bind(&PLY_reader::face_end, this)
               );
    }
    else if (element_name == "edge") {
        m_tissue->setNumWall(count);
        return std::tr1::tuple<std::tr1::function<void()>, std::tr1::function<void()> >(
                   std::tr1::bind(&PLY_reader::edge_begin, this),
                   std::tr1::bind(&PLY_reader::edge_end, this)
               );
    }
    else {
        return std::tr1::tuple<std::tr1::function<void()>, std::tr1::function<void()> >(0, 0);
    }
}
//----------------------------------------------------------------------------
void PLY_reader::vertex_begin()
{

}
//----------------------------------------------------------------------------
void PLY_reader::vertex_x_callback(ply::float32 x)
{
    m_position[0] = x;
}
//----------------------------------------------------------------------------
void PLY_reader::vertex_y_callback(ply::float32 y)
{
    m_position[1] = y;
}
//----------------------------------------------------------------------------
void PLY_reader::vertex_z_callback(ply::float32 z)
{
    m_position[2] = z;
}
//----------------------------------------------------------------------------
void PLY_reader::vertex_index_callback(ply::uint32 i)
{
    m_index = i;
}
//----------------------------------------------------------------------------
void PLY_reader::vertex_end()
{
    Vertex * v = m_tissue->vertexP(m_index);
    v->setIndex(m_index);
    v->setPosition(m_position);
}
//----------------------------------------------------------------------------
void PLY_reader::face_begin()
{
    m_vertices.clear();
    m_variables.clear();
}
//----------------------------------------------------------------------------
void PLY_reader::face_vertex_indices_begin(ply::uint8 size)
{
    m_vertices.reserve(size);
}
//----------------------------------------------------------------------------
void PLY_reader::face_vertex_indices_element(ply::uint32 vertex_index)
{
    m_vertices.push_back( m_tissue->vertexP(vertex_index) );
}
//----------------------------------------------------------------------------
void PLY_reader::face_vertex_indices_end()
{
  
}
//----------------------------------------------------------------------------
void PLY_reader::face_index_callback(ply::uint32 i)
{
    m_index = i;
}
//----------------------------------------------------------------------------
void PLY_reader::face_variables_begin(ply::uint8 size)
{
    m_variables.reserve(size);
}
//----------------------------------------------------------------------------
void PLY_reader::face_variables_element(ply::float32 variable)
{
    m_variables.push_back(variable);
}
//----------------------------------------------------------------------------
void PLY_reader::face_variables_end()
{

}
//----------------------------------------------------------------------------
void PLY_reader::face_center_vertex_begin(ply::uint8 size)
{
    m_ceter_vertex.reserve(size);
}
//----------------------------------------------------------------------------
void PLY_reader::face_center_vertex_element(ply::float32 variable)
{
    m_ceter_vertex.push_back(variable);
}
//----------------------------------------------------------------------------
void PLY_reader::face_center_vertex_end()
{

}
//----------------------------------------------------------------------------
void PLY_reader::face_edges_lengths_begin(ply::uint8 size)
{
    m_edges_lengths.reserve(size);
}
//----------------------------------------------------------------------------
void PLY_reader::face_edges_lengths_element(ply::float32 variable)
{
    m_edges_lengths.push_back(variable);
}
//----------------------------------------------------------------------------
void PLY_reader::face_edges_lengths_end()
{

}
//----------------------------------------------------------------------------
void PLY_reader::face_end()
{
    Cell * c = m_tissue->cellP(m_index);
    c->setIndex(m_index);
    c->setVertex(m_vertices);
    c->setNumVariable(m_variables.size());
    c->setVariable(m_variables);
    c->setCenterPosition(m_ceter_vertex);
    c->setEdgeLength(m_edges_lengths);
    std::vector<Vertex*>::iterator it, v_end;
    for (it = m_vertices.begin(), v_end = m_vertices.end(); it != v_end; ++it)
    {
        (*it)->addCell(c);
    }
}
//----------------------------------------------------------------------------
void PLY_reader::edge_begin()
{
    m_vertices.clear();
    m_vertices.resize(2);
    m_variables.clear();
}
//----------------------------------------------------------------------------
void PLY_reader::edge_source_callback(ply::uint32 size)
{
    m_vertices[0] = m_tissue->vertexP(size);
}
//----------------------------------------------------------------------------
void PLY_reader::edge_target_callback(ply::uint32 size)
{
    m_vertices[1] = m_tissue->vertexP(size);
}
//----------------------------------------------------------------------------
void PLY_reader::edge_index_callback(ply::uint32 i)
{
    m_index = i;
}
//----------------------------------------------------------------------------
void PLY_reader::edge_variables_begin(ply::uint8 size)
{
    m_variables.reserve(size);
}
//----------------------------------------------------------------------------
void PLY_reader::edge_variables_element(ply::float32 variable)
{
    m_variables.push_back(variable);
}
//----------------------------------------------------------------------------
void PLY_reader::edge_variables_end()
{

}
//----------------------------------------------------------------------------
void PLY_reader::edge_end()
{
    Wall * w = m_tissue->wallP(m_index);
    w->setIndex(m_index);
    w->setVertex(m_vertices[0], m_vertices[1]);
    w->setNumVariable(m_variables.size());
    w->setVariable(m_variables);
    std::vector<Vertex*>::iterator it, v_end;
    for (it = m_vertices.begin(), v_end = m_vertices.end(); it != v_end; ++it)
    {
        (*it)->addWall(w);
    }
}
//----------------------------------------------------------------------------
void PLY_reader::info_callback(const std::string& filename, std::size_t line_number, const std::string& message)
{
    std::cerr << filename << ":" << line_number << ": " << "info: " << message << std::endl;
}
//----------------------------------------------------------------------------
void PLY_reader::warning_callback(const std::string& filename, std::size_t line_number, const std::string& message)
{
    std::cerr << filename << ":" << line_number << ": " << "warning: " << message << std::endl;
}
//----------------------------------------------------------------------------
void PLY_reader::error_callback(const std::string& filename, std::size_t line_number, const std::string& message)
{
    std::cerr << filename << ":" << line_number << ": " << "error: " << message << std::endl;
}
//----------------------------------------------------------------------------
void PLY_reader::set_cell_wall_connectivity(Tissue &t)
{
    std::vector<size_t> v1_cells, v2_cells;
    v1_cells.reserve(10);
    v2_cells.reserve(10);

    //iterate all walls
    for (size_t i = 0; i < t.numWall(); ++i)
    {
        Wall * w = t.wallP(i);
        v1_cells.clear();
        v2_cells.clear();
        
        //for each vertex of the wall find intersection of its cells
        Vertex *v1 = w->vertex1();
        Vertex *v2 = w->vertex2();
        for (size_t j = 0; j < v1->numCell(); ++j)
            v1_cells.push_back(v1->cell(j)->index());
        for (size_t j = 0; j < v2->numCell(); ++j)
            v2_cells.push_back(v2->cell(j)->index());
        std::sort(v1_cells.begin(), v1_cells.end());
        std::sort(v2_cells.begin(), v2_cells.end());
        std::vector<size_t>::iterator it=std::set_intersection (v1_cells.begin(), v1_cells.end(), v2_cells.begin(), v2_cells.end(), v1_cells.begin());
        v1_cells.resize(it-v1_cells.begin());

        //set wall cells and cell walls acordingly
        assert(v1_cells.size() <= 2);
        std::vector<Cell*> cells(2, NULL);
        std::vector<Cell*>::iterator cell_it = cells.begin();
        for (it = v1_cells.begin(); it != v1_cells.end(); ++it)
        {
            Cell *c = m_tissue->cellP(*it);
            c->addWall(w);
            *cell_it++ = c;
        }
        w->setCell(cells[0], cells[1]);
    }
}
//----------------------------------------------------------------------------
