//
// File:    VTUostream.cc
// Author:  pkrupinski
// Created: September 15, 2011, 3:20 PM
//
#include "VTUostream.h"
#include "tissue.h"
#include "vertex.h"
// #include <unordered_set>
const double WALL_RELATIVE_THICKNESS = 0.1;
using namespace IO;
//-----------------------------------------------------------------------------

Point::Point ( Vertex const& vv ) : x ( 0.0 ), y ( 0.0 ), z ( 0.0 )
{
    Vertex &v = const_cast<Vertex &> ( vv );
    int n = v.numPosition();
    if ( n > 2 )
    {
        x = v.position ( 0 );
        y = v.position ( 1 );
        z = v.position ( 2 );
    }
    else if ( n > 1 )
    {
        x = v.position ( 0 );
        y = v.position ( 1 );
    }
    else if ( n > 0 )
        x = v.position ( 0 );
}
//-----------------------------------------------------------------------------

Point const& Point::displace ( Vector const&v, double d )
{
    x += d * v.x;
    y += d * v.y;
    z += d * v.z;
    return *this;
}
//-----------------------------------------------------------------------------

VTUostream::VTUostream() : m_os ( 0 ), D ( WALL_RELATIVE_THICKNESS )
{
}
//-----------------------------------------------------------------------------

VTUostream::VTUostream ( std::ostream& o ) : m_os ( &o ), D ( WALL_RELATIVE_THICKNESS )
{
    header();
}
//-----------------------------------------------------------------------------
//BEGIN cell and wall geometry for walls as line segments
void VTUostream::write_cells ( Tissue const& t )
{
    typedef std::vector<Cell>::const_iterator CellIter;
    std::vector<Cell> const& cells = t.cell();
    CellIter cit, cend;
    bool triangles = true;
    for ( cit = cells.begin(), cend = cells.end(); cit != cend; ++cit )
    {
        if ( cit->numVertex() != 3 )
            triangles = false;
    }
    write_piece_header ( t.numVertex(), t.numCell() );
    write_cell_point_geometry ( t );
    if ( triangles )
        write_cell_geometry ( t, VTUostream::TRIANGLE );
    else
        write_cell_geometry ( t, VTUostream::POLYGON );
    write_cell_data_header ( "Scalars=\"cell variable 4\" Vectors=\"cell vector\"" );
    //write_cell_data3V ( t );  //writes three cell vectors
    write_cell_data ( t );  //writes one cell vectors
    write_cell_data_footer();
    write_piece_footer();



}
//-----------------------------------------------------------------------------
void VTUostream::write_walls ( Tissue const& t )
{
    typedef std::vector<Cell>::const_iterator CellIter;
    std::vector<Cell> const& cells = t.cell();
    int npts = 0;
    CellIter cit, cend;
    for ( cit = cells.begin(), cend = cells.end(); cit != cend; ++cit )
        npts += cit->numWall();

    write_piece_header ( npts, npts );
    write_wall_point_geometry ( t );
    write_wall_geometry ( t );
    write_wall_data_header ( "Scalars=\"wall variable 0\"" );
    write_wall_data ( t );
    write_wall_data_footer();
    write_piece_footer();
}
//-----------------------------------------------------------------------------
void VTUostream::write_cell_point_geometry ( Tissue const& t )
{
    typedef std::vector<Vertex>::const_iterator VertexIter;
    std::vector<Vertex> const& vertices = t.vertex();

    VertexIter viter, vend;
    *m_os << "<Points>\n"
          << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for ( viter = vertices.begin(), vend = vertices.end(); viter != vend; ++viter )
    {
        std::vector<double> const& p = viter->position();
        size_t npos = viter->numPosition();
        for ( size_t j = 0; j < npos; ++j )
            *m_os << p[j] << " ";
        if ( npos < 3 )
            *m_os << "0 ";
        *m_os << "\n";
    }
    *m_os << "</DataArray>\n"
          << "</Points>\n";
}
//-----------------------------------------------------------------------------
void VTUostream::write_cell_geometry ( Tissue const& t, Cell_type ct )
{
    typedef std::vector<Cell>::const_iterator CellIter;
    std::vector<Cell> const& cells = t.cell();

    *m_os << "<Cells>\n"
          << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    CellIter cit, cend;
    for ( cit = cells.begin(), cend = cells.end(); cit != cend; ++cit )
    {
        // this part relays on the experimental 'unordered_set' but does not assume that vertices are cyclically ordered vertices in the cell
        //     typedef std::vector<Wall*>::const_iterator WallIter;
        //     std::vector<Wall*> const& walls = cit->wall();
        //     std::unordered_set<size_t> indices;
        //     typedef std::unordered_set<size_t>::iterator IndexIter;
        //     WallIter wit, wend;
        //     wit = walls.begin();
        //     Vertex const *v1 = (*wit)->vertex1(), *v2 = (*wit)->vertex2();
        //     if(const_cast<Cell&>(*cit).hasVertex(const_cast<Vertex*>(v1)))
        //       indices.insert( v1->index() );
        //     else if(const_cast<Cell&>(*cit).hasVertex(const_cast<Vertex*>(v2)))
        //       indices.insert( v2->index() );
        //     else{
        //       std::cerr << "Cell wall vertives do not agree with cell vertices. Writing to file aborted!\n";
        //       return;
        //     }
        //     ++wit;
        //     for (wend = walls.end(); wit != wend; ++wit)
        //     {
        //       indices.insert( (*wit)->vertex1()->index() );
        //       indices.insert( (*wit)->vertex2()->index() );
        //     }
        //     IndexIter iiter, iend;
        //     for (iiter = indices.begin(), iend = indices.end(); iiter != iend; ++iiter)
        //     {
        //       *m_os << *iiter << " ";
        //     }
        //     *m_os << "\n";
        // this part needs cyclically ordered vertices in the cell
        typedef std::vector<Vertex*>::const_iterator CVIter;
        std::vector<Vertex*> const& verts = cit->vertex();
        CVIter cviter, cvend;
        for ( cviter = verts.begin(), cvend = verts.end(); cviter != cvend; ++cviter )
        {
            *m_os << ( *cviter )->index() << " ";
        }
        *m_os << "\n";
    }

    *m_os << "</DataArray>\n"
          << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    int total_offset = 0;
    for ( cit = cells.begin(), cend = cells.end(); cit != cend; ++cit )
    {
        total_offset += cit->numWall();
        *m_os << total_offset << " ";
    }

    *m_os << "\n"
          << "</DataArray>\n"
          << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for ( size_t i = 0; i < t.numCell(); ++i )
    {
        *m_os << ct << " ";
    }
    *m_os << "\n";
    //*/
    *m_os << "</DataArray>\n"
          << "</Cells>\n";
}
//-----------------------------------------------------------------------------
void VTUostream::write_wall_point_geometry ( Tissue const& t )
{
    const double D = 0.05;
    typedef std::vector<Cell>::const_iterator CellIter;
    std::vector<Cell> const& cells = t.cell();
    std::vector<Point> points;
    points.reserve ( t.numWall() *2 );
    CellIter cit, cend;
    for ( cit = cells.begin(), cend = cells.end(); cit != cend; ++cit )
    {
        Cell &c = const_cast<Cell&> ( *cit );
        typedef std::vector<Wall*>::const_iterator WallIter;
        std::vector<Wall*> const& walls = c.wall();

        std::vector<double> cent ( 3 );
        cent = c.positionFromVertex();
        Point center ( cent[0], cent[1], cent[2] );

        WallIter wit, wend;
        wit = walls.begin();
        Vertex const *w1v1 = ( *wit )->vertex1(), *w1v2 = ( *wit )->vertex2();
        ++wit;
        Vertex const *w2v1 = ( *wit )->vertex1(), *w2v2 = ( *wit )->vertex2();
        if ( w1v1 == w2v1 ) //w1v2, (w1v1 = w2v1), w2v2
        {
            points.push_back ( Point ( *w1v2 ).displace_towards ( center, D ) );
            points.push_back ( Point ( *w1v1 ).displace_towards ( center, D ) );
            w1v2 = w2v2;
        }
        else if ( w1v1 == w2v2 ) //w1v2, (w1v1 = w2v2), w2v1
        {
            points.push_back ( Point ( *w1v2 ).displace_towards ( center, D ) );
            points.push_back ( Point ( *w1v1 ).displace_towards ( center, D ) );
            w1v2 = w2v1;
        }
        else if ( w1v2 == w2v1 ) //w1v1, (w1v2 = w2v1), w2v2
        {
            points.push_back ( Point ( *w1v1 ).displace_towards ( center, D ) );
            points.push_back ( Point ( *w1v2 ).displace_towards ( center, D ) );
            w1v2 = w2v2;
        }
        else if ( w1v2 == w2v2 ) //w1v1, (w1v2 = w2v2), w2v1
        {
            points.push_back ( Point ( *w1v1 ).displace_towards ( center, D ) );
            points.push_back ( Point ( *w1v2 ).displace_towards ( center, D ) );
            w1v2 = w2v1;
        }
        else
        {
            std::cerr << "VTUostream::write_wall_point_geometry(Tissue const& t); Consecutive cell walls do not share a vertex.\n";
            exit ( EXIT_FAILURE );
        }
        ++wit;
        for ( wend = walls.end(); wit != wend; ++wit )
        {
            w2v1 = ( *wit )->vertex1();
            w2v2 = ( *wit )->vertex2();
            if ( w1v2 == w2v1 )
            {
                points.push_back ( Point ( *w2v1 ).displace_towards ( center, D ) );
                w1v2 = w2v2;
            }
            else if ( w1v2 == w2v2 )
            {
                points.push_back ( Point ( *w2v2 ).displace_towards ( center, D ) );
                w1v2 = w2v1;
            }
            else
            {
                std::cerr << "VTUostream::write_wall_point_geometry(Tissue const& t); Consecutive cell walls do not share a vertex.\n";
                exit ( EXIT_FAILURE );
            }
        }
    }
    typedef std::vector<Point>::const_iterator PointIter;

    PointIter piter, pend;
    *m_os << "<Points>\n"
          << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for ( piter = points.begin(), pend = points.end(); piter != pend; ++piter )
    {
        *m_os << piter->x << " " << piter->y << " " << piter->z << "\n";
    }
    *m_os << "</DataArray>\n"
          << "</Points>\n";
}
//-----------------------------------------------------------------------------
void VTUostream::write_wall_geometry ( Tissue const& t )
{
    typedef std::vector<Cell>::const_iterator CellIter;
    std::vector<Cell> const& cells = t.cell();

    *m_os << "<Cells>\n"
          << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    int count = 0;
    CellIter cit, cend;
    for ( cit = cells.begin(), cend = cells.end(); cit != cend; ++cit )
    {
        int nwall = cit->numWall() - 1;
        for ( int i = 0; i < nwall; ++i )
        {
            *m_os << count + i << " " << count + i + 1 << " ";
        }
        *m_os << count + nwall << " " << count << " ";
        count += cit->numWall();

    }

    *m_os << "</DataArray>\n"
          << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    int total_offset = 0;
    for ( int i = 0; i < count; ++i )
    {
        total_offset += 2;
        *m_os << total_offset << " ";
    }
    *m_os << "\n"

          << "</DataArray>\n"
          << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for ( int i = 0; i < count; ++i )
    {
        *m_os << 3 << " ";
    }
    *m_os << "\n";
    //*/
    *m_os << "</DataArray>\n"
          << "</Cells>\n";
}
//END cell and wall geometry for walls as line segments
//-----------------------------------------------------------------------------
//BEGIN cell and wall geometry for 2D walls
void VTUostream::write_cells2 ( Tissue const& t )
{
    typedef std::vector<Cell>::const_iterator CellIter;
    std::vector<Cell> const& cells = t.cell();
    int npts = 0;
    CellIter cit, cend;
    bool triangles = true;
    for ( cit = cells.begin(), cend = cells.end(); cit != cend; ++cit )
    {
        int nvrt = cit->numVertex();
        if ( nvrt != 3 )
            triangles = false;
        npts += nvrt;
    }

    write_piece_header ( npts, t.numCell() );
    write_cell_point_geometry2 ( t );
    if ( triangles )
        write_cell_geometry2 ( t, VTUostream::TRIANGLE );
    else
        write_cell_geometry2 ( t, VTUostream::POLYGON );
    write_cell_data_header ( "Scalars=\"cell variable 4\" Vectors=\"cell vector\"" );
  
    //  write_cell_data ( t );
    write_cell_data3V ( t );
    
    write_cell_data_footer();
    write_piece_footer();
}
//-----------------------------------------------------------------------------
// single walls
void VTUostream::write_walls2 ( Tissue const& t )
{
    typedef std::vector<Cell>::const_iterator CellIter;
    std::vector<Cell> const& cells = t.cell();
    int ncell = 0;
    CellIter cit, cend;
    for ( cit = cells.begin(), cend = cells.end(); cit != cend; ++cit )
        ncell += cit->numWall();
    ;
    std::vector<Vertex*> verts;
    verts.reserve ( t.numWall() *2 );

    write_piece_header ( t.numVertex() + ncell, ncell );
    write_wall_point_geometry2 ( t, verts );
    write_wall_geometry2 ( t, verts );
    write_wall_data_header ( "Scalars=\"wall variable 0\"" );
    write_wall_data ( t );
    write_wall_data_footer();
    write_piece_footer();
}
//-----------------------------------------------------------------------------
// double walls
void VTUostream::write_walls3 ( Tissue const& t )
{
    typedef std::vector<Cell>::const_iterator CellIter;
    std::vector<Cell> const& cells = t.cell();
    int ncell = 0;
    CellIter cit, cend;
    for ( cit = cells.begin(), cend = cells.end(); cit != cend; ++cit )
        ncell += cit->numWall();

    std::vector<Vertex*> verts;
    verts.reserve ( t.numWall() *2 );
    write_piece_header ( t.numVertex() + ncell, ncell );
    write_wall_point_geometry2 ( t, verts );
    write_wall_geometry2 ( t, verts );
    write_wall_data_header ( "Scalars=\"wall variable 0\"" );
    //  std::cout << "write_wall data2\n";
    write_wall_data2 ( t );
    write_wall_data_footer();
    write_piece_footer();
}
//-----------------------------------------------------------------------------
void VTUostream::write_cell_point_geometry2 ( Tissue const& t )
{
    typedef std::vector<Cell>::const_iterator CellIter;
    std::vector<Cell> const& cells = t.cell();

    *m_os << "<Points>\n"
          << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    CellIter cit, cend;
    for ( cit = cells.begin(), cend = cells.end(); cit != cend; ++cit )
    {
        Cell &c = const_cast<Cell&> ( *cit );
	//HJ: removed due to unused variable warning
        //typedef std::vector<Wall*>::const_iterator WallIter;
        //std::vector<Wall*> const& walls = c.wall();

        std::vector<double> cent ( 3 );
        cent = c.positionFromVertex();
        Point center ( cent[0], cent[1], cent[2] );

        typedef std::vector<Vertex*>::const_iterator CVIter;
        std::vector<Vertex*> const& verts = c.vertex();
        CVIter cviter, cvend;
        for ( cviter = verts.begin(), cvend = verts.end(); cviter != cvend; ++cviter )
        {
            Point p = Point ( **cviter ).displace_towards ( center, D );
            *m_os << p.x << " " << p.y << " " << p.z << "\n";
        }
    }
    *m_os << "</DataArray>\n"
          << "</Points>\n";
}
//-----------------------------------------------------------------------------
void VTUostream::write_cell_geometry2 ( Tissue const& t, Cell_type ct )
{
    typedef std::vector<Cell>::const_iterator CellIter;
    std::vector<Cell> const& cells = t.cell();

    *m_os << "<Cells>\n"
          << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    CellIter cit, cend;
    int count = 0;
    for ( cit = cells.begin(), cend = cells.end(); cit != cend; ++cit )
    {
        int nverts = cit->numVertex();
        for ( int i = 0; i < nverts; ++i, ++count )
        {
            *m_os << count << " ";
        }
        *m_os << "\n";
    }

    *m_os << "</DataArray>\n"
          << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    int total_offset = 0;
    for ( cit = cells.begin(), cend = cells.end(); cit != cend; ++cit )
    {
        total_offset += cit->numVertex();
        *m_os << total_offset << " ";
    }
    *m_os << "\n"

          << "</DataArray>\n"
          << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for ( size_t i = 0; i < t.numCell(); ++i )
    {
        *m_os << ct << " ";
    }
    *m_os << "\n";
    *m_os << "</DataArray>\n"
          << "</Cells>\n";
}
//-----------------------------------------------------------------------------
void VTUostream::write_wall_point_geometry2 ( Tissue const& t, std::vector<Vertex*> &verts )
{
    //write all the vertices first in the order of their indices
    typedef std::vector<Vertex>::const_iterator VertexIter;
    std::vector<Vertex> const& vertices = t.vertex();

    VertexIter viter, vend;
    *m_os << "<Points>\n"
          << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for ( viter = vertices.begin(), vend = vertices.end(); viter != vend; ++viter )
    {
        std::vector<double> const& p = viter->position();
        size_t npos = viter->numPosition();
        for ( size_t j = 0; j < npos; ++j )
            *m_os << p[j] << " ";
        if ( npos < 3 )
            *m_os << "0 ";
        *m_os << "\n";
    }

    //write the displaced vertices of each cell
    typedef std::vector<Cell>::const_iterator CellIter;
    std::vector<Cell> const& cells = t.cell();

    CellIter cit, cend;
    for ( cit = cells.begin(), cend = cells.end(); cit != cend; ++cit )
    {
        Cell &c = const_cast<Cell&> ( *cit );
        typedef std::vector<Wall*>::const_iterator WallIter;
        std::vector<Wall*> const& walls = c.wall();

        WallIter wit, wend;
        wit = walls.begin();
        Vertex *w1v1 = ( *wit )->vertex1(), *w1v2 = ( *wit )->vertex2();
        ++wit;
        Vertex *w2v1 = ( *wit )->vertex1(), *w2v2 = ( *wit )->vertex2();
        if ( w1v1 == w2v1 ) //w1v2, (w1v1 = w2v1), w2v2
        {
            verts.push_back ( w1v2 );
            verts.push_back ( w1v1 );
            w1v2 = w2v2;
        }
        else if ( w1v1 == w2v2 ) //w1v2, (w1v1 = w2v2), w2v1
        {
            verts.push_back ( w1v2 );
            verts.push_back ( w1v1 );
            w1v2 = w2v1;
        }
        else if ( w1v2 == w2v1 ) //w1v1, (w1v2 = w2v1), w2v2
        {
            verts.push_back ( w1v1 );
            verts.push_back ( w1v2 );
            w1v2 = w2v2;
        }
        else if ( w1v2 == w2v2 ) //w1v1, (w1v2 = w2v2), w2v1
        {
            verts.push_back ( w1v1 );
            verts.push_back ( w1v2 );
            w1v2 = w2v1;
        }
        else
        {
            std::cerr << "VTUostream::write_wall_point_geometry(Tissue const& t); Consecutive cell walls do not share a vertex.\n";
            exit ( EXIT_FAILURE );
        }
        ++wit;
        for ( wend = walls.end(); wit != wend; ++wit )
        {
            w2v1 = ( *wit )->vertex1();
            w2v2 = ( *wit )->vertex2();
            if ( w1v2 == w2v1 )
            {
                verts.push_back ( w2v1 );
                w1v2 = w2v2;
            }
            else if ( w1v2 == w2v2 )
            {
                verts.push_back ( w2v2 );
                w1v2 = w2v1;
            }
            else
            {
                std::cerr << "VTUostream::write_wall_point_geometry(Tissue const& t); Consecutive cell walls do not share a vertex.\n";
                exit ( EXIT_FAILURE );
            }
        }
    }
    typedef std::vector<Vertex*>::const_iterator VertexPIter;
    VertexPIter vpit = verts.begin(); //, vpend = verts.end();
    for ( cit = cells.begin(), cend = cells.end(); cit != cend; ++cit )
    {
        Cell &c = const_cast<Cell&> ( *cit );
        std::vector<double> cent ( 3 );
        cent = c.positionFromVertex();
        Point center ( cent[0], cent[1], cent[2] );
        
        int nwall = cit->numWall();
        for ( int i = 0; i < nwall; ++i )

        {
            Point p = Point ( **vpit++ ).displace_towards ( center, D );
            *m_os << p.x << " " << p.y << " " << p.z << "\n";
        }
    }
    *m_os << "</DataArray>\n"
          << "</Points>\n";
}
//-----------------------------------------------------------------------------
void VTUostream::write_wall_geometry2 ( Tissue const& t, std::vector<Vertex*> const& verts )
{
    typedef std::vector<Cell>::const_iterator CellIter;
    std::vector<Cell> const& cells = t.cell();
    *m_os << "<Cells>\n"
          << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    int count = 0, offset = t.numVertex();
    //construct quad walls by ordering vertices circularly: first 2 original verices of the wall then 2 displaced vertices
    std::vector<Vertex*> ::const_iterator vit = verts.begin(), vstart; //, vend = verts.end();
    CellIter cit, cend;
    for ( cit = cells.begin(), cend = cells.end(); cit != cend; ++cit )
    {
        vstart = vit;
        int nwall = cit->numWall() - 1, temp = offset + count;
        for ( int i = 0; i < nwall; ++i )
        {
            *m_os << ( *vit++ )->index() << " ";
            *m_os << ( *vit )->index() << " " << temp + i + 1 << " " << temp + i << " ";
        }
        *m_os << ( *vit++ )->index() << " " << ( *vstart )->index() << " " << temp << " " << temp + nwall << " ";
        count += cit->numWall();
    }
    *m_os << "\n";
    *m_os << "</DataArray>\n"
          << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    int total_offset = 0;
    for ( int i = 0; i < count; ++i )
    {
        total_offset += 4;
        *m_os << total_offset << " ";
    }
    *m_os << "\n"

          << "</DataArray>\n"
          << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for ( int i = 0; i < count; ++i )
    {
        *m_os << 7 << " ";
    }
    *m_os << "\n";
    //*/
    *m_os << "</DataArray>\n"
          << "</Cells>\n";
}
//END cell and wall geometry for 2D walls
//-----------------------------------------------------------------------------
//BEGIN cell and wall geometry for 2D walls printing inner and outer cell walls separately based on the flag in the last wall variable
void VTUostream::write_cells3 ( Tissue const& t )
{
    typedef std::vector<Cell>::const_iterator CellIter;
    std::vector<Cell> const& cells = t.cell();
    std::vector<IO::Point> disp_points;
    std::vector<char> vertex_flag;
    std::vector< std::map<size_t,size_t> > cvp_map;

    int npts = 0;
    CellIter cit, cend;
    bool triangles = true;
    for ( cit = cells.begin(), cend = cells.end(); cit != cend; ++cit )
    {
        int nvrt = cit->numVertex();
        if ( nvrt != 3 )
            triangles = false;
        npts += nvrt;
    }

    size_t flag_pos = t.wall(0).numVariable()-1;
    write_piece_header ( npts, t.numCell() );
    write_cell_point_geometry3 ( t, disp_points, vertex_flag, cvp_map, flag_pos );

    if ( triangles )
        write_cell_geometry3 ( t, VTUostream::TRIANGLE );
    else
        write_cell_geometry3 ( t, VTUostream::POLYGON );
    write_cell_data_header ( "Scalars=\"cell variable 4\" Vectors=\"cell vector\"" );
    write_cell_data ( t );
    write_cell_data_footer();
    write_piece_footer();
}
//-----------------------------------------------------------------------------
void VTUostream::write_inner_walls ( Tissue const& t )
{
    std::vector<char> vertex_flag;
    std::vector<uint> index_map;
    size_t flag_pos = t.wall(0).numVariable()-1;
        double FLAG_VAL = 0.0;
    size_t num_cell = prepare_marked_vertices ( t, vertex_flag, flag_pos, FLAG_VAL )*2;
    size_t num_verts = std::count ( vertex_flag.begin(), vertex_flag.end(), 1 );

    write_piece_header ( num_verts, num_cell );
    write_inner_wall_point_geometry3 ( t, vertex_flag, index_map, num_cell );
    write_inner_wall_geometry3 ( t, index_map, flag_pos , FLAG_VAL );
    write_wall_data_header ( "Scalars=\"wall variable 0\"" );
    write_wall_data3 ( t, flag_pos, FLAG_VAL, true );
    write_wall_data_footer();
    write_piece_footer();
}
//-----------------------------------------------------------------------------
void VTUostream::write_outer_walls ( Tissue const& t )
{
    std::vector<Point> disp_points;
    std::vector<char> vertex_flag;
    std::vector< std::map<size_t,size_t> > cvp_map;
    std::vector<uint> index_map;
    size_t flag_pos = t.wall(0).numVariable()-1;
    double FLAG_VAL = 1.0;
    
    std::pair<size_t, size_t> num = prepare_wall_point_geometry3 ( t, disp_points,  vertex_flag,  cvp_map, flag_pos, FLAG_VAL );
    write_piece_header ( disp_points.size() + num.first, num.second );
    write_outer_wall_point_geometry3 ( t, disp_points,  vertex_flag,  index_map );
    write_outer_wall_geometry3 ( t, index_map, cvp_map, num.first, flag_pos , FLAG_VAL );
    write_wall_data_header ( "Scalars=\"wall variable 0\"" );
    write_wall_data3 ( t, flag_pos, FLAG_VAL, false );
    write_wall_data_footer();
    write_piece_footer();
}
//-----------------------------------------------------------------------------
void VTUostream::write_outer_line_walls ( Tissue const& t )
{
    std::vector<Point> disp_points;
    std::vector<char> vertex_flag;
    std::vector< std::map<size_t,size_t> > cvp_map;
    std::vector<uint> index_map;
    size_t flag_pos = t.wall(0).numVariable()-1;
    double FLAG_VAL = 1.0;
    
    std::pair<size_t, size_t> num = prepare_wall_point_geometry3 ( t, disp_points,  vertex_flag,  cvp_map, flag_pos, FLAG_VAL );
    write_piece_header ( disp_points.size() + num.first, num.second );
    write_outer_wall_point_geometry3 ( t, disp_points,  vertex_flag,  index_map );
    write_outer_wall_geometry3 ( t, index_map, cvp_map, num.first, flag_pos , FLAG_VAL );
    write_wall_data_header ( "Scalars=\"wall variable 0\"" );
    write_wall_data3 ( t, flag_pos, FLAG_VAL, false );
    write_wall_data_footer();
    write_piece_footer();
}
//-----------------------------------------------------------------------------
void VTUostream::write_cell_point_geometry3 ( Tissue const& t, std::vector<Point>& disp_points, std::vector<char>& vertex_flag, std::vector< std::map<size_t,size_t> >& cvp_map, size_t flag_pos )
{
    typedef std::vector<Cell>::const_iterator CellIter;
    std::vector<Cell> const& cells = t.cell();
    cvp_map.resize ( t.numCell() ); //maps for each cell displaced vertices point indices

    size_t counter = prepare_marked_vertices ( t, vertex_flag, flag_pos, 1.0 );
    disp_points.reserve ( 2*counter ); //contains displaced points for the vertices that needed displacements indexed by vertex index

    *m_os << "<Points>\n"
          << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    CellIter cit, cend;
    for ( cit = cells.begin(), cend = cells.end(); cit != cend; ++cit )
    {
        Cell &c = const_cast<Cell&> ( *cit );
        std::vector<double> cent ( 3 );
        typedef std::vector<Vertex*>::const_iterator CVIter;
        std::vector<Vertex*> const& verts = c.vertex();
        CVIter cviter, cvend, cvprev = --verts.end(), cvnext = verts.begin();
        for ( cviter = verts.begin(), cvend = verts.end(); cviter != cvend; ++cviter )
        {
            ++cvnext;
            if ( cvnext == cvend )
                cvnext = verts.begin();

            Point p = Point ( **cviter );
            size_t index = ( *cviter )->index();
            if ( vertex_flag[ index ] )
            {
                //figure out vector to displace point along
                Vector prev_vec ( p, Point ( **cvprev ) );
                Vector next_vec ( p, Point ( **cvnext ) );
                if ( vertex_flag[ ( *cvprev )->index() ] )
                {
                    if ( vertex_flag[ ( *cvnext )->index() ] )
                        p.displace ( prev_vec+next_vec, D );
                    else
                        p.displace ( next_vec, D );
                }
                else
                {
                    if ( vertex_flag[ ( *cvnext )->index() ] )
                        p.displace ( prev_vec, D );
                }
                cvp_map[cit->index()].insert ( std::pair<size_t, size_t> ( index, disp_points.size() ) ); //maps vertex index to index in disp_points
                disp_points.push_back ( p );
            }
            *m_os << p.x << " " << p.y << " " << p.z << "\n";
            cvprev = cviter;
        }
    }

    *m_os << "</DataArray>\n"
          << "</Points>\n";
}
//-----------------------------------------------------------------------------
void VTUostream::write_cell_geometry3 ( Tissue const& t, Cell_type ct )
{
    typedef std::vector<Cell>::const_iterator CellIter;
    std::vector<Cell> const& cells = t.cell();

    *m_os << "<Cells>\n"
          << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    CellIter cit, cend;
    int count = 0;
    for ( cit = cells.begin(), cend = cells.end(); cit != cend; ++cit )
    {
        int nverts = cit->numVertex();
        for ( int i = 0; i < nverts; ++i, ++count )
        {
            *m_os << count << " ";
        }
        *m_os << "\n";
    }

    *m_os << "</DataArray>\n"
          << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    int total_offset = 0;
    for ( cit = cells.begin(), cend = cells.end(); cit != cend; ++cit )
    {
        total_offset += cit->numVertex();
        *m_os << total_offset << " ";
    }
    *m_os << "\n"

          << "</DataArray>\n"
          << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for ( size_t i = 0; i < t.numCell(); ++i )
    {
        *m_os << ct << " ";
    }
    *m_os << "\n";
    *m_os << "</DataArray>\n"
          << "</Cells>\n";
}
//-----------------------------------------------------------------------------
std::pair<size_t, size_t> VTUostream::prepare_wall_point_geometry3 ( Tissue const& t, std::vector<Point>& disp_points, std::vector<char>& vertex_flag, std::vector< std::map<size_t,size_t> >& cvp_map, size_t flag_pos, double flag_val )
{
    typedef std::vector<Cell>::const_iterator CellIter;
    std::vector<Cell> const& cells = t.cell();
    cvp_map.resize ( t.numCell() ); //maps for each cell displaced vertices point indices
    size_t counter = prepare_marked_vertices ( t, vertex_flag, flag_pos, flag_val );
    disp_points.reserve ( 2*counter ); //contains displaced points. Indices of those points are mapped later to cell and vertex indices

    counter = 0;
    CellIter cit, cend;
    for ( cit = cells.begin(), cend = cells.end(); cit != cend; ++cit )
    {
        Cell &c = const_cast<Cell&> ( *cit );
        //count the walls that need displacement
        typedef std::vector<Wall*>::const_iterator CWIter;
        CWIter cwiter, cwend;
        for ( cwiter = c.wall().begin(), cwend = c.wall().end(); cwiter != cwend; ++cwiter )
        {
            Wall const& w = **cwiter;
            double flag = w.variable ( flag_pos ); //flag for inner and outer walls
            if ( flag == flag_val )
                ++counter;
        }
        //displace vertices of these walls
        typedef std::vector<Vertex*>::const_iterator CVIter;
        std::vector<Vertex*> const& verts = c.vertex();
        CVIter cviter, cvend, cvprev = --verts.end(), cvnext = verts.begin();
        for ( cviter = verts.begin(), cvend = verts.end(); cviter != cvend; ++cviter )
        {
            ++cvnext;
            if ( cvnext == cvend )
                cvnext = verts.begin();

            Point p = Point ( **cviter );
            size_t index = ( *cviter )->index();
            if ( vertex_flag[ index ] )
            {
                //figure out vector to displace point along
                Vector prev_vec ( p, Point ( **cvprev ) );
                Vector next_vec ( p, Point ( **cvnext ) );
                if ( vertex_flag[ ( *cvprev )->index() ] )
                {
                    if ( vertex_flag[ ( *cvnext )->index() ] )
                        p.displace ( prev_vec+next_vec, D );
                    else
                        p.displace ( next_vec, D );
                }
                else
                {
                    if ( vertex_flag[ ( *cvnext )->index() ] )
                        p.displace ( prev_vec, D );
                }

                cvp_map[cit->index()].insert ( std::pair<size_t, size_t> ( index, disp_points.size() ) ); //maps vertex index to index in disp_points
                disp_points.push_back ( p );
            }
            cvprev = cviter;
        }
    }
    size_t nverts = std::count ( vertex_flag.begin(), vertex_flag.end(), 1 );
    return std::make_pair ( nverts, counter );
}
//-----------------------------------------------------------------------------
void VTUostream::write_outer_wall_point_geometry3 ( Tissue const& t, std::vector<Point>& disp_points, std::vector<char>& vertex_flag, std::vector<uint>& index_map )
{
    //write all the outer vertices first in the order of their indices and create the map from tissue index to file index
    typedef std::vector<Vertex>::const_iterator VertexIter;
    std::vector<Vertex> const& vertices = t.vertex();
    index_map.resize ( t.numVertex(), 0 ); //index map of verices to ordering indices in the file

    VertexIter viter, vend;
    *m_os << "<Points>\n"
          << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    uint counter = 0;
    for ( viter = vertices.begin(), vend = vertices.end(); viter != vend; ++viter )
    {
        size_t index = viter->index();
        if ( vertex_flag[index] == 1 ) //outer vertex
        {
            std::vector<double> const& p = viter->position();
            size_t npos = viter->numPosition();
            for ( size_t j = 0; j < npos; ++j )
                *m_os << p[j] << " ";
            if ( npos < 3 )
                *m_os << "0 ";
            *m_os << "\n";

            index_map[index] = counter;
            ++counter;
        }
    }

    //write the displaced vertices of each cell using precomputed points
    typedef std::vector<Point>::const_iterator PIter;
    PIter piter, pend;
    for ( piter = disp_points.begin(), pend = disp_points.end(); piter != pend; ++piter )
    {
        *m_os << piter->x << " " << piter->y << " " << piter->z << "\n";
    }

    *m_os << "</DataArray>\n"
          << "</Points>\n";
}
//-----------------------------------------------------------------------------
void VTUostream::write_outer_wall_geometry3 ( Tissue const& t, std::vector<uint>& index_map, std::vector< std::map<size_t,size_t> >& cvp_map, size_t offset, size_t flag_pos, double flag_val )
{
    typedef std::vector<Cell>::const_iterator CellIter;
    std::vector<Cell> const& cells = t.cell();
    *m_os << "<Cells>\n"
          << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";

    CellIter cit, cend;
    int count = 0;
    for ( cit = cells.begin(), cend = cells.end(); cit != cend; ++cit )
    {
        Cell &c = const_cast<Cell&> ( *cit );

        typedef std::vector<Wall*>::const_iterator CWIter;
        CWIter cwiter, cwend;
        for ( cwiter = c.wall().begin(), cwend = c.wall().end(); cwiter != cwend; ++cwiter )
        {
            Wall const& w = **cwiter;
            double flag = w.variable ( flag_pos ); //flag for inner and outer walls
            if ( flag == flag_val )
            {
                std::map<size_t,size_t> &map = cvp_map[c.index()];
                size_t i1 = w.vertex1()->index();
                size_t i2 = w.vertex2()->index();
                std::map<size_t,size_t>::iterator it1 = map.find ( i1 );
                std::map<size_t,size_t>::iterator it2 = map.find ( i2 );
                if ( it1 == map.end() || it2 == map.end() )
                {
                    std::cerr << "VTUostream::write_outer_wall_geometry3( Tissue const& t, std::vector<uint>& index_map, std::vector< std::map<size_t,size_t> >& cvp_map, size_t flag_pos ); Outer wall vertex index not found in displaced vertices map.\n";
                    exit ( EXIT_FAILURE );
                }
                //write indicess of points in the file constituing the wall quadrilateral
                *m_os << index_map[i1] << " " << index_map[i2] << " " << offset + it2->second << " " << offset + it1->second << " ";
                ++count;
            }
        }
    }
    *m_os << "\n";
    *m_os << "</DataArray>\n"
          << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    int total_offset = 0;
    for ( int i = 0; i < count; ++i )
    {
        total_offset += 4;
        *m_os << total_offset << " ";
    }
    *m_os << "\n"

          << "</DataArray>\n"
          << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for ( int i = 0; i < count; ++i )
    {
        *m_os << 7 << " ";
    }
    *m_os << "\n";
    //*/
    *m_os << "</DataArray>\n"
          << "</Cells>\n";
}
//-----------------------------------------------------------------------------
void VTUostream::write_inner_wall_point_geometry3 ( Tissue const& t, std::vector<char>& vertex_flag, std::vector<uint>& index_map, size_t counter )
{
  //HJ: removed due to unused variable warning
  //typedef std::vector<Cell>::const_iterator CellIter;
    typedef std::vector<Vertex>::const_iterator VertexIter;
    std::vector<Vertex> const& vertices = t.vertex();
    index_map.resize ( t.numVertex(), 0 ); //index map of verices to ordering indices in the file

    VertexIter viter, vend;
    *m_os << "<Points>\n"
          << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    counter = 0;
    for ( viter = vertices.begin(), vend = vertices.end(); viter != vend; ++viter )
    {
        size_t index = viter->index();
        if ( vertex_flag[index] == 1 ) //inner vertex
        {
            std::vector<double> const& p = viter->position();
            size_t npos = viter->numPosition();
            for ( size_t j = 0; j < npos; ++j )
                *m_os << p[j] << " ";
            if ( npos < 3 )
                *m_os << "0 ";
            *m_os << "\n";

            index_map[index] = counter;
            ++counter;
        }
    }
    *m_os << "</DataArray>\n"
          << "</Points>\n";
}
//-----------------------------------------------------------------------------
void VTUostream::write_inner_wall_geometry3 ( Tissue const& t, std::vector<uint>& index_map, size_t flag_pos, double flag_val )
{
    typedef std::vector<Cell>::const_iterator CellIter;
    std::vector<Cell> const& cells = t.cell();
    *m_os << "<Cells>\n"
          << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";

    CellIter cit, cend;
    int count = 0;
    for ( cit = cells.begin(), cend = cells.end(); cit != cend; ++cit )
    {
        Cell &c = const_cast<Cell&> ( *cit );

        typedef std::vector<Wall*>::const_iterator CWIter;
        CWIter cwiter, cwend;
        for ( cwiter = c.wall().begin(), cwend = c.wall().end(); cwiter != cwend; ++cwiter )
        {
            Wall const& w = **cwiter;
            double flag = w.variable ( flag_pos ); //flag for inner and outer walls
            if ( flag == flag_val )
            {

                size_t i1 = w.vertex1()->index();
                size_t i2 = w.vertex2()->index();

                //write indicess of points in the file constituing the wall line
                *m_os << index_map[i1] << " " << index_map[i2] << " ";
                ++count;
            }
        }
    }
    *m_os << "\n";

    *m_os << "</DataArray>\n"
          << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    int total_offset = 0;
    for ( int i = 0; i < count; ++i )
    {
        total_offset += 2;
        *m_os << total_offset << " ";
    }
    *m_os << "\n"

          << "</DataArray>\n"
          << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for ( int i = 0; i < count; ++i )
    {
        *m_os << 3 << " ";
    }
    *m_os << "\n";
    *m_os << "</DataArray>\n"
          << "</Cells>\n";
}
//-----------------------------------------------------------------------------
size_t VTUostream::prepare_marked_vertices ( Tissue const& t, std::vector<char>& vertex_flag, size_t flag_pos, double flag_val )
{
    vertex_flag.resize ( t.numVertex(), 0 ); //flag for which vertices to displace
    typedef std::vector<Wall>::const_iterator WallIter;
    std::vector<Wall> const& walls = t.wall();
    WallIter witer, wend;
    size_t counter = 0;
    for ( witer = walls.begin(), wend = walls.end(); witer != wend; ++witer )
    {
        Wall const& w = *witer;
        double flag = w.variable ( flag_pos );
        if ( flag == flag_val )
        {
            vertex_flag[w.vertex1()->index()] = 1;
            vertex_flag[w.vertex2()->index()] = 1;
            ++counter;
        }
    }
    return counter;
}
//END cell and wall geometry for 2D walls printing inner and outer cell walls separately based on the flag in the last wall variable
//-----------------------------------------------------------------------------
void VTUostream::write_piece_header ( int n_pts, int n_cell )
{
    //    mark = m_os.tellp();
    *m_os << "<Piece  NumberOfPoints=\"";
    m_os->width ( 8 );
    *m_os << n_pts;
    *m_os << "\" NumberOfCells=\"";
    m_os->width ( 8 );
    *m_os << n_cell << "\">";
    *m_os << "\n";
}
//-----------------------------------------------------------------------------
void VTUostream::write_piece_footer()
{
    *m_os << "</Piece>\n";
}
//-----------------------------------------------------------------------------
void VTUostream::write_cell_data ( Tissue const& t )
{
    typedef std::vector<Cell>::const_iterator CellIter;
    std::vector<Cell> const& cells = t.cell();

    CellIter cit = cells.begin(), cend;
    int nvars = cit->numVariable();
    //write a cell vector data assuming first 3 cell variables are vector components and 4th is a length
    *m_os << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"cell vector\" format=\"ascii\">\n";
    for ( cit = cells.begin(), cend = cells.end(); cit != cend; ++cit )
    {
        Cell &c = const_cast<Cell&> ( *cit );
        *m_os << c.variable ( 0 ) << " " << c.variable ( 1 ) << " " << c.variable ( 2 ) << "\n";
    }
    *m_os << "</DataArray>\n";
    *m_os << "<DataArray type=\"Float64\" Name=\"cell vector length\" format=\"ascii\">\n";
    for ( cit = cells.begin(), cend = cells.end(); cit != cend; ++cit )
    {
        Cell &c = const_cast<Cell&> ( *cit );
        *m_os << c.variable ( 3 ) << " ";
    }
    *m_os << "\n"<< "</DataArray>\n";
    //write rest of the cell data
    for ( int i = 4; i < nvars; ++i )
    {
        *m_os << "<DataArray type=\"Float64\" Name=\"cell variable " << i << "\" format=\"ascii\">\n";
        for ( cit = cells.begin(), cend = cells.end(); cit != cend; ++cit )
        {
            Cell &c = const_cast<Cell&> ( *cit );
            *m_os << c.variable ( i ) << " ";
        }
        *m_os << "\n"
              << "</DataArray>\n";
    }
}
//----------------------------------------------------------for having 3 cell vectors
void VTUostream::write_cell_data3V ( Tissue const& t )
{
    typedef std::vector<Cell>::const_iterator CellIter;
    std::vector<Cell> const& cells = t.cell();

    CellIter cit = cells.begin(), cend;
    int nvars = cit->numVariable();
    //write 3 cell vector data assuming 
    // 0,1,2 cell variables are 1st vector components and 3 is a length
    // 4,5,6 cell variables are 1st vector components and 7 is a length
    // 8,9,10 cell variables are 1st vector components and 11 is a length
    *m_os << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"cell vector1\" format=\"ascii\">\n";
    for ( cit = cells.begin(), cend = cells.end(); cit != cend; ++cit )
    {
        Cell &c = const_cast<Cell&> ( *cit );
        *m_os << c.variable ( 0 ) << " " << c.variable ( 1 ) << " " << c.variable ( 2 ) << "\n";
    }
    *m_os << "</DataArray>\n";
    *m_os << "<DataArray type=\"Float64\" Name=\"cell vector1 length\" format=\"ascii\">\n";
    for ( cit = cells.begin(), cend = cells.end(); cit != cend; ++cit )
    {
        Cell &c = const_cast<Cell&> ( *cit );
        *m_os << c.variable ( 3 ) << " ";
    }
    *m_os << "\n"<< "</DataArray>\n";

    *m_os << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"cell vector2\" format=\"ascii\">\n";
    for ( cit = cells.begin(), cend = cells.end(); cit != cend; ++cit )
    {
        Cell &c = const_cast<Cell&> ( *cit );
        *m_os << c.variable ( 4 ) << " " << c.variable ( 5 ) << " " << c.variable ( 6 ) << "\n";
    }
    *m_os << "</DataArray>\n";
    *m_os << "<DataArray type=\"Float64\" Name=\"cell vector2 length\" format=\"ascii\">\n";
    for ( cit = cells.begin(), cend = cells.end(); cit != cend; ++cit )
    {
        Cell &c = const_cast<Cell&> ( *cit );
        *m_os << c.variable ( 7 ) << " ";
    }
    *m_os << "\n"<< "</DataArray>\n";

    *m_os << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"cell vector3\" format=\"ascii\">\n";
    for ( cit = cells.begin(), cend = cells.end(); cit != cend; ++cit )
    {
        Cell &c = const_cast<Cell&> ( *cit );
        *m_os << c.variable ( 8 ) << " " << c.variable ( 9 ) << " " << c.variable ( 10 ) << "\n";
    }
    *m_os << "</DataArray>\n";
    *m_os << "<DataArray type=\"Float64\" Name=\"cell vector3 length\" format=\"ascii\">\n";
    for ( cit = cells.begin(), cend = cells.end(); cit != cend; ++cit )
    {
        Cell &c = const_cast<Cell&> ( *cit );
        *m_os << c.variable ( 11 ) << " ";
    }
    *m_os << "\n"<< "</DataArray>\n";

    //write rest of the cell data
    for ( int i = 12; i < nvars; ++i )
    {
        *m_os << "<DataArray type=\"Float64\" Name=\"cell variable " << i << "\" format=\"ascii\">\n";
        for ( cit = cells.begin(), cend = cells.end(); cit != cend; ++cit )
        {
            Cell &c = const_cast<Cell&> ( *cit );
            *m_os << c.variable ( i ) << " ";
        }
        *m_os << "\n"
              << "</DataArray>\n";
    }
}
//-----------------------------------------------------------------------------
void VTUostream::write_wall_data ( Tissue const& t )
{
    typedef std::vector<Cell>::const_iterator CellIter;
    std::vector<Cell> const& cells = t.cell();
    CellIter cit = cells.begin(), cend;
    Cell &c = const_cast<Cell&> ( *cit );

    //Print wall lengths
    *m_os << "<DataArray type=\"Float64\" Name=\"wall length\" format=\"ascii\">\n";
    for ( cit = cells.begin(), cend = cells.end(); cit != cend; ++cit )
    {
        Cell &c = const_cast<Cell&> ( *cit );
        typedef std::vector<Wall*>::const_iterator WallIter;
        std::vector<Wall*> const& walls = c.wall();
        WallIter wit, wend;
        for ( wit = walls.begin(), wend = walls.end(); wit != wend; ++wit )
        {
            Wall &w = const_cast<Wall&> ( **wit );
            *m_os << w.length() << " ";
        }
    }
    *m_os << "\n"
          << "</DataArray>\n";

    // Print wall variables
    int nvars = c.wall ( 0 )->numVariable();
    for ( int i = 0; i < nvars; ++i )
    {
        *m_os << "<DataArray type=\"Float64\" Name=\"wall variable " << i << "\" format=\"ascii\">\n";
        for ( cit = cells.begin(), cend = cells.end(); cit != cend; ++cit )
        {
            Cell &c = const_cast<Cell&> ( *cit );
            typedef std::vector<Wall*>::const_iterator WallIter;
            std::vector<Wall*> const& walls = c.wall();
            WallIter wit, wend;
            for ( wit = walls.begin(), wend = walls.end(); wit != wend; ++wit )
            {
                Wall &w = const_cast<Wall&> ( **wit );
                *m_os << w.variable ( i ) << " ";
            }
        }
        *m_os << "\n"
              << "</DataArray>\n";
    }
}
//-----------------------------------------------------------------------------
void VTUostream::write_wall_data2 ( Tissue const& t )
{
    typedef std::vector<Cell>::const_iterator CellIter;
    std::vector<Cell> const& cells = t.cell();
    CellIter cit = cells.begin(), cend;
    Cell &c = const_cast<Cell&> ( *cit );
    //Print wall lengths
    *m_os << "<DataArray type=\"Float64\" Name=\"wall length\" format=\"ascii\">\n";
    for ( cit = cells.begin(), cend = cells.end(); cit != cend; ++cit )
    {
        Cell &c = const_cast<Cell&> ( *cit );
        typedef std::vector<Wall*>::const_iterator WallIter;
        std::vector<Wall*> const& walls = c.wall();
        WallIter wit, wend;
        for ( wit = walls.begin(), wend = walls.end(); wit != wend; ++wit )
        {
            Wall &w = const_cast<Wall&> ( **wit );
            *m_os << w.length() << " ";
        }
    }
    *m_os << "\n"
          << "</DataArray>\n";
    // Print wall variables assuming paired structure
    //total number of variables in the wall
    int nvars = c.wall ( 0 )->numVariable();
    //check if nvars is odd so it can be of the form (length, v1c1, v1c2, v2c1, v2c2, ...)
    //    if (!(nvars & 0x1))
    //    {
    //        std::cerr << "VTUostream::write_wall_data2(Tissue const& t); number of wall variables does not fit requested format\n";
    //        exit(EXIT_FAILURE);
    //    }
    for ( int i = 0; i < nvars; i+=2 )
    {
        *m_os << "<DataArray type=\"Float64\" Name=\"wall variable " << i/2 << "\" format=\"ascii\">\n";
        for ( cit = cells.begin(), cend = cells.end(); cit != cend; ++cit )
        {
            Cell &c = const_cast<Cell&> ( *cit );
            size_t c_id = c.index();
            typedef std::vector<Wall*>::const_iterator WallIter;
            std::vector<Wall*> const& walls = c.wall();
            WallIter wit, wend;
            for ( wit = walls.begin(), wend = walls.end(); wit != wend; ++wit )
            {
                Wall &w = const_cast<Wall&> ( **wit );
                size_t c1id = w.cell1()->index();
                size_t c2id = w.cell2()->index();
                int j = -1;
                if ( c_id == c1id )
                    j = i;
                else if ( c_id == c2id )
                    j = i+1;
                else
                {
                    std::cerr << "VTUostream::write_wall_data2(Tissue const& t); "
                              << "Wall does not report connection to the cell "
                              << "which was accessed through\n";
                    exit ( EXIT_FAILURE );
                }
                *m_os << w.variable ( j ) << " ";
            }
        }
        *m_os << "\n"
              << "</DataArray>\n";
    }
}
//-----------------------------------------------------------------------------
void VTUostream::write_wall_data3 ( Tissue const& t, size_t flag_pos, double flag_val, bool project_cell_variables )
{
    typedef std::vector<Cell>::const_iterator CellIter;
    std::vector<Cell> const& cells = t.cell();
    CellIter cit = cells.begin(), cend;
    Cell &c = const_cast<Cell&> ( *cit );
    //Print wall lengths
    *m_os << "<DataArray type=\"Float64\" Name=\"wall length\" format=\"ascii\">\n";
    for ( cit = cells.begin(), cend = cells.end(); cit != cend; ++cit )
    {
        Cell &c = const_cast<Cell&> ( *cit );
        typedef std::vector<Wall*>::const_iterator WallIter;
        std::vector<Wall*> const& walls = c.wall();
        WallIter wit, wend;
        for ( wit = walls.begin(), wend = walls.end(); wit != wend; ++wit )
        {
            Wall &w = const_cast<Wall&> ( **wit );
            double flag = w.variable ( flag_pos ); //flag for inner and outer walls
            if ( flag == flag_val )
                *m_os << w.length() << " ";
        }
    }
    *m_os << "\n"
          << "</DataArray>\n";
    // Print wall variables assuming paired structure
    //total number of variables in the wall
    int nvars = c.wall ( 0 )->numVariable();
// //     check if nvars is odd so it can be of the form (length, v1c1, v1c2, v2c1, v2c2, ...)
//     if ( ! ( nvars & 0x1 ) )
//     {
//         std::cerr << "VTUostream::write_wall_data3 ( Tissue const& t, size_t flag_pos, double flag_val ); number of wall variables does not fit requested format\n";
//         exit ( EXIT_FAILURE );
//     }
    for ( int i = 0; i < nvars; i+=2 )
    {
        *m_os << "<DataArray type=\"Float64\" Name=\"wall variable " << i/2 << "\" format=\"ascii\">\n";
        for ( cit = cells.begin(), cend = cells.end(); cit != cend; ++cit )
        {
            Cell &c = const_cast<Cell&> ( *cit );
            size_t c_id = c.index();
            typedef std::vector<Wall*>::const_iterator WallIter;
            std::vector<Wall*> const& walls = c.wall();
            WallIter wit, wend;
            for ( wit = walls.begin(), wend = walls.end(); wit != wend; ++wit )
            {
                Wall &w = const_cast<Wall&> ( **wit );
                double flag = w.variable ( flag_pos ); //flag for inner and outer walls
                if ( flag == flag_val )
                {
                    size_t c1id = w.cell1()->index();
                    size_t c2id = w.cell2()->index();
                    int j = -1;
                    if ( c_id == c1id )
                        j = i;
                    else if ( c_id == c2id )
                        j = i+1;
                    else
                    {
                        std::cerr << "VTUostream::write_wall_data3(Tissue const& t, size_t flag_pos, double flag_val ); "
                                  << "Wall does not report connection to the cell which it was accessed through\n";
                        exit ( EXIT_FAILURE );
                    }
                    *m_os << w.variable ( j ) << " ";
                }
            }
        }
        *m_os << "\n"
              << "</DataArray>\n";
    }

    //project cell variable on the wall
    if ( project_cell_variables )
    {
        //write a cell vector data assuming first 3 cell variables are vector components and 4th is a length
        *m_os << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"cell vector\" format=\"ascii\">\n";
        for ( cit = cells.begin(), cend = cells.end(); cit != cend; ++cit )
        {
            Cell &c = const_cast<Cell&> ( *cit );
            typedef std::vector<Wall*>::const_iterator WallIter;
            std::vector<Wall*> const& walls = c.wall();
            WallIter wit, wend;
            for ( wit = walls.begin(), wend = walls.end(); wit != wend; ++wit )
            {
                Wall &w = const_cast<Wall&> ( **wit );
                double flag = w.variable ( flag_pos ); //flag for inner and outer walls
                if ( flag == flag_val )
                    *m_os << c.variable ( 0 ) << " " << c.variable ( 1 ) << " " << c.variable ( 2 ) << "\n";
            }
        }
        *m_os << "</DataArray>\n";
        *m_os << "<DataArray type=\"Float64\" Name=\"cell vector length\" format=\"ascii\">\n";
        for ( cit = cells.begin(), cend = cells.end(); cit != cend; ++cit )
        {
            Cell &c = const_cast<Cell&> ( *cit );
            typedef std::vector<Wall*>::const_iterator WallIter;
            std::vector<Wall*> const& walls = c.wall();
            WallIter wit, wend;
            for ( wit = walls.begin(), wend = walls.end(); wit != wend; ++wit )
            {
                Wall &w = const_cast<Wall&> ( **wit );
                double flag = w.variable ( flag_pos ); //flag for inner and outer walls
                if ( flag == flag_val )
                    *m_os << c.variable ( 3 ) << " ";
            }
        }
        *m_os << "\n"<< "</DataArray>\n";
        //write rest of the cell data
        nvars = c.numVariable();
        for ( int i = 4; i < nvars; ++i )
        {
            *m_os << "<DataArray type=\"Float64\" Name=\"cell variable " << i << "\" format=\"ascii\">\n";
            for ( cit = cells.begin(), cend = cells.end(); cit != cend; ++cit )
            {
                Cell &c = const_cast<Cell&> ( *cit );
                typedef std::vector<Wall*>::const_iterator WallIter;
                std::vector<Wall*> const& walls = c.wall();
                WallIter wit, wend;
                for ( wit = walls.begin(), wend = walls.end(); wit != wend; ++wit )
                {
                    Wall &w = const_cast<Wall&> ( **wit );
                    double flag = w.variable ( flag_pos ); //flag for inner and outer walls
                    if ( flag == flag_val )
                        *m_os << c.variable ( i ) << " ";
                }
            }
            *m_os << "\n"
                  << "</DataArray>\n";
        }
    }
}
//-----------------------------------------------------------------------------
