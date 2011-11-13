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

Point::Point(Vertex const& vv) : x(0.0), y(0.0), z(0.0)
{
    Vertex &v = const_cast<Vertex &> (vv);
    int n = v.numPosition();
    if (n > 2)
    {
        x = v.position(0);
        y = v.position(1);
        z = v.position(2);
    }
    else if (n > 1)
    {
        x = v.position(0);
        y = v.position(1);
    }
    else if (n > 0)
        x = v.position(0);
}
//-----------------------------------------------------------------------------

Point const& Point::displace(Vector const&v, double d)
{
    x += d * v.x;
    y += d * v.y;
    z += d * v.z;
    return *this;
}
//-----------------------------------------------------------------------------

VTUostream::VTUostream() : m_os(0), D(WALL_RELATIVE_THICKNESS)
{
}
//-----------------------------------------------------------------------------

VTUostream::VTUostream(std::ostream& o) : m_os(&o), D(WALL_RELATIVE_THICKNESS)
{
    header();
}
//-----------------------------------------------------------------------------

void VTUostream::write_cells(Tissue const& t)
{
    write_piece_header(t.numVertex(), t.numCell());
    write_cell_point_geometry(t);
    write_cell_geometry(t);
    write_cell_data_header("Scalars=\"cell variable 0\"");
    write_cell_data(t);
    write_cell_data_footer();
    write_piece_footer();
}
//-----------------------------------------------------------------------------

void VTUostream::write_cells2(Tissue const& t)
{
    typedef std::vector<Cell>::const_iterator CellIter;
    std::vector<Cell> const& cells = t.cell();
    int npts = 0;
    CellIter cit, cend;
    for (cit = cells.begin(), cend = cells.end(); cit != cend; ++cit)
        npts += cit->numVertex();

    write_piece_header(npts, t.numCell());
    write_cell_point_geometry2(t);
    write_cell_geometry2(t);
    write_cell_data_header("Scalars=\"cell variable 0\"");
    write_cell_data(t);
    write_cell_data_footer();
    write_piece_footer();
}
//-----------------------------------------------------------------------------

void VTUostream::write_walls(Tissue const& t)
{
    typedef std::vector<Cell>::const_iterator CellIter;
    std::vector<Cell> const& cells = t.cell();
    int npts = 0;
    CellIter cit, cend;
    for (cit = cells.begin(), cend = cells.end(); cit != cend; ++cit)
        npts += cit->numWall();

    write_piece_header(npts, npts);
    write_wall_point_geometry(t);
    write_wall_geometry(t);
    write_wall_data_header("Scalars=\"wall variable 0\"");
    write_wall_data(t);
    write_wall_data_footer();
    write_piece_footer();
}
//-----------------------------------------------------------------------------
void VTUostream::write_walls2(Tissue const& t)
{
	typedef std::vector<Cell>::const_iterator CellIter;
	std::vector<Cell> const& cells = t.cell();
	int ncell = 0;
	CellIter cit, cend;
	for (cit = cells.begin(), cend = cells.end(); cit != cend; ++cit)
		ncell += cit->numWall();
	;
	std::vector<Vertex*> verts;
	verts.reserve(t.numWall()*2);
	
	write_piece_header(t.numVertex() + ncell, ncell);
	write_wall_point_geometry2(t, verts);
	write_wall_geometry2(t, verts);
	write_wall_data_header("Scalars=\"wall variable 0\"");
	write_wall_data(t);
	write_wall_data_footer();
	write_piece_footer();
}
//-----------------------------------------------------------------------------

void VTUostream::write_walls3(Tissue const& t)
{
    typedef std::vector<Cell>::const_iterator CellIter;
    std::vector<Cell> const& cells = t.cell();
    int ncell = 0;
    CellIter cit, cend;
    for (cit = cells.begin(), cend = cells.end(); cit != cend; ++cit)
        ncell += cit->numWall();
    ;
    std::vector<Vertex*> verts;
    verts.reserve(t.numWall()*2);

    write_piece_header(t.numVertex() + ncell, ncell);
    write_wall_point_geometry2(t, verts);
    write_wall_geometry2(t, verts);
    write_wall_data_header("Scalars=\"wall variable 0\"");
    //  std::cout << "write_wall data2\n";
    write_wall_data2(t);
    write_wall_data_footer();
    write_piece_footer();
}
//-----------------------------------------------------------------------------

void VTUostream::write_piece_header(int n_pts, int n_cell)
{
    //    mark = m_os.tellp();
    *m_os << "<Piece  NumberOfPoints=\"";
    m_os->width(8);
    *m_os << n_pts;
    *m_os << "\" NumberOfCells=\"";
    m_os->width(8);
    *m_os << n_cell << "\">";
    *m_os << "\n";
}
//-----------------------------------------------------------------------------

void VTUostream::write_piece_footer()
{
    *m_os << "</Piece>\n";
}
//-----------------------------------------------------------------------------

void VTUostream::write_cell_point_geometry(Tissue const& t)
{
    typedef std::vector<Vertex>::const_iterator VertexIter;
    std::vector<Vertex> const& vertices = t.vertex();

    VertexIter viter, vend;
    *m_os << "<Points>\n"
            << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (viter = vertices.begin(), vend = vertices.end(); viter != vend; ++viter)
    {
        std::vector<double> const& p = viter->position();
        size_t npos = viter->numPosition();
        for (size_t j = 0; j < npos; ++j)
            *m_os << p[j] << " ";
        if (npos < 3)
            *m_os << "0 ";
        *m_os << "\n";
    }
    *m_os << "</DataArray>\n"
            << "</Points>\n";
}
//-----------------------------------------------------------------------------

void VTUostream::write_cell_geometry(Tissue const& t)
{
    typedef std::vector<Cell>::const_iterator CellIter;
    std::vector<Cell> const& cells = t.cell();

    *m_os << "<Cells>\n"
            << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    CellIter cit, cend;
    for (cit = cells.begin(), cend = cells.end(); cit != cend; ++cit)
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
        for (cviter = verts.begin(), cvend = verts.end(); cviter != cvend; ++cviter)
        {
            *m_os << (*cviter)->index() << " ";
        }
        *m_os << "\n";
    }

    *m_os << "</DataArray>\n"
            << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    int total_offset = 0;
    for (cit = cells.begin(), cend = cells.end(); cit != cend; ++cit)
    {
        total_offset += cit->numWall();
        *m_os << total_offset << " ";
    }

    *m_os << "\n"
            << "</DataArray>\n"
            << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (int i = 0; i < t.numCell(); ++i)
    {
        *m_os << 7 << " ";
    }
    *m_os << "\n";
    //*/
    *m_os << "</DataArray>\n"
            << "</Cells>\n";
}
//-----------------------------------------------------------------------------

void VTUostream::write_wall_point_geometry(Tissue const& t)
{
    const double D = 0.05;
    typedef std::vector<Cell>::const_iterator CellIter;
    std::vector<Cell> const& cells = t.cell();
    std::vector<Point> points;
    points.reserve(t.numWall()*2);
    CellIter cit, cend;
    for (cit = cells.begin(), cend = cells.end(); cit != cend; ++cit)
    {
        Cell &c = const_cast<Cell&> (*cit);
        typedef std::vector<Wall*>::const_iterator WallIter;
        std::vector<Wall*> const& walls = c.wall();

        std::vector<double> cent(3);
        cent = c.positionFromVertex();
        Point center(cent[0], cent[1], cent[2]);

        WallIter wit, wend;
        wit = walls.begin();
        Vertex const *w1v1 = (*wit)->vertex1(), *w1v2 = (*wit)->vertex2();
        ++wit;
        Vertex const *w2v1 = (*wit)->vertex1(), *w2v2 = (*wit)->vertex2();
        if (w1v1 == w2v1) //w1v2, (w1v1 = w2v1), w2v2
        {
            points.push_back(Point(*w1v2).displace_towards(center, D));
            points.push_back(Point(*w1v1).displace_towards(center, D));
            w1v2 = w2v2;
        }
        else if (w1v1 == w2v2) //w1v2, (w1v1 = w2v2), w2v1
        {
            points.push_back(Point(*w1v2).displace_towards(center, D));
            points.push_back(Point(*w1v1).displace_towards(center, D));
            w1v2 = w2v1;
        }
        else if (w1v2 == w2v1) //w1v1, (w1v2 = w2v1), w2v2
        {
            points.push_back(Point(*w1v1).displace_towards(center, D));
            points.push_back(Point(*w1v2).displace_towards(center, D));
            w1v2 = w2v2;
        }
        else if (w1v2 == w2v2) //w1v1, (w1v2 = w2v2), w2v1
        {
            points.push_back(Point(*w1v1).displace_towards(center, D));
            points.push_back(Point(*w1v2).displace_towards(center, D));
            w1v2 = w2v1;
        }
        else
        {
            std::cerr << "VTUostream::write_wall_point_geometry(Tissue const& t); Consecutive cell walls do not share a vertex.\n";
            exit(EXIT_FAILURE);
        }
        ++wit;
        for (wend = walls.end(); wit != wend; ++wit)
        {
            w2v1 = (*wit)->vertex1();
            w2v2 = (*wit)->vertex2();
            if (w1v2 == w2v1)
            {
                points.push_back(Point(*w2v1).displace_towards(center, D));
                w1v2 = w2v2;
            }
            else if (w1v2 == w2v2)
            {
                points.push_back(Point(*w2v2).displace_towards(center, D));
                w1v2 = w2v1;
            }
            else
            {
                std::cerr << "VTUostream::write_wall_point_geometry(Tissue const& t); Consecutive cell walls do not share a vertex.\n";
                exit(EXIT_FAILURE);
            }
        }
    }
    typedef std::vector<Point>::const_iterator PointIter;

    PointIter piter, pend;
    *m_os << "<Points>\n"
            << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (piter = points.begin(), pend = points.end(); piter != pend; ++piter)
    {
        *m_os << piter->x << " " << piter->y << " " << piter->z << "\n";
    }
    *m_os << "</DataArray>\n"
            << "</Points>\n";
}
//-----------------------------------------------------------------------------

void VTUostream::write_wall_geometry(Tissue const& t)
{
    typedef std::vector<Cell>::const_iterator CellIter;
    std::vector<Cell> const& cells = t.cell();

    *m_os << "<Cells>\n"
            << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    int count = 0;
    CellIter cit, cend;
    for (cit = cells.begin(), cend = cells.end(); cit != cend; ++cit)
    {
        int nwall = cit->numWall() - 1;
        for (int i = 0; i < nwall; ++i)
        {
            *m_os << count + i << " " << count + i + 1 << " ";
        }
        *m_os << count + nwall << " " << count << " ";
        count += cit->numWall();

    }

    *m_os << "</DataArray>\n"
            << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    int total_offset = 0;
    for (int i = 0; i < count; ++i)
    {
        total_offset += 2;
        *m_os << total_offset << " ";
    }
    *m_os << "\n"

            << "</DataArray>\n"
            << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (int i = 0; i < count; ++i)
    {
        *m_os << 3 << " ";
    }
    *m_os << "\n";
    //*/
    *m_os << "</DataArray>\n"
            << "</Cells>\n";
}
//-----------------------------------------------------------------------------

void VTUostream::write_cell_data(Tissue const& t)
{
    typedef std::vector<Cell>::const_iterator CellIter;
    std::vector<Cell> const& cells = t.cell();

    CellIter cit = cells.begin(), cend;
    int nvars = cit->numVariable();
    for (int i = 0; i < nvars; ++i)
    {
        *m_os << "<DataArray type=\"Float64\" Name=\"cell variable " << i << "\" format=\"ascii\">\n";
        for (cit = cells.begin(), cend = cells.end(); cit != cend; ++cit)
        {
            Cell &c = const_cast<Cell&> (*cit);
            *m_os << c.variable(i) << " ";
        }
        *m_os << "\n"
                << "</DataArray>\n";
    }
}
//-----------------------------------------------------------------------------

void VTUostream::write_wall_data(Tissue const& t)
{
    typedef std::vector<Cell>::const_iterator CellIter;
    std::vector<Cell> const& cells = t.cell();
    CellIter cit = cells.begin(), cend;
    Cell &c = const_cast<Cell&> (*cit);

		//Print wall lengths
    *m_os << "<DataArray type=\"Float64\" Name=\"wall length\" format=\"ascii\">\n";
    for (cit = cells.begin(), cend = cells.end(); cit != cend; ++cit)
			{
				Cell &c = const_cast<Cell&>(*cit);
				typedef std::vector<Wall*>::const_iterator WallIter;
				std::vector<Wall*> const& walls = c.wall();
				WallIter wit, wend;
				for (wit = walls.begin(), wend = walls.end(); wit != wend; ++wit)
					{
						Wall &w = const_cast<Wall&>(**wit);
						*m_os << w.length() << " ";
					}
			}
    *m_os << "\n"
					<< "</DataArray>\n";

		// Print wall variables
    int nvars = c.wall(0)->numVariable();
    for (int i = 0; i < nvars; ++i)
    {
        *m_os << "<DataArray type=\"Float64\" Name=\"wall variable " << i << "\" format=\"ascii\">\n";
        for (cit = cells.begin(), cend = cells.end(); cit != cend; ++cit)
        {
            Cell &c = const_cast<Cell&> (*cit);
            typedef std::vector<Wall*>::const_iterator WallIter;
            std::vector<Wall*> const& walls = c.wall();
            WallIter wit, wend;
            for (wit = walls.begin(), wend = walls.end(); wit != wend; ++wit)
            {
                Wall &w = const_cast<Wall&> (**wit);
                *m_os << w.variable(i) << " ";
            }
        }
        *m_os << "\n"
                << "</DataArray>\n";
    }
}
//-----------------------------------------------------------------------------

void VTUostream::write_wall_data2(Tissue const& t)
{
    typedef std::vector<Cell>::const_iterator CellIter;
    std::vector<Cell> const& cells = t.cell();
    CellIter cit = cells.begin(), cend;
    Cell &c = const_cast<Cell&> (*cit);

		//Print wall lengths
    *m_os << "<DataArray type=\"Float64\" Name=\"wall length\" format=\"ascii\">\n";
    for (cit = cells.begin(), cend = cells.end(); cit != cend; ++cit)
			{
				Cell &c = const_cast<Cell&>(*cit);
				typedef std::vector<Wall*>::const_iterator WallIter;
				std::vector<Wall*> const& walls = c.wall();
				WallIter wit, wend;
				for (wit = walls.begin(), wend = walls.end(); wit != wend; ++wit)
					{
						Wall &w = const_cast<Wall&>(**wit);
						*m_os << w.length() << " ";
					}
			}
    *m_os << "\n"
					<< "</DataArray>\n";
		
		// Print wall variables assuming paired structure
    //total number of variables in the wall
    int nvars = c.wall(0)->numVariable();
    //check if nvars is odd so it can be of the form (length, v1c1, v1c2, v2c1, v2c2, ...)
//    if (!(nvars & 0x1))
//    {
//        std::cerr << "VTUostream::write_wall_data2(Tissue const& t); number of wall variables does not fit requested format\n";
//        exit(EXIT_FAILURE);
//    }
    for (int i = 0; i < nvars; i+=2)
			{
				*m_os << "<DataArray type=\"Float64\" Name=\"wall variable " << i << "\" format=\"ascii\">\n";
				for (cit = cells.begin(), cend = cells.end(); cit != cend; ++cit)
					{
            Cell &c = const_cast<Cell&> (*cit);
            size_t c_id = c.index();
            typedef std::vector<Wall*>::const_iterator WallIter;
            std::vector<Wall*> const& walls = c.wall();
            WallIter wit, wend;
            for (wit = walls.begin(), wend = walls.end(); wit != wend; ++wit)
							{
                Wall &w = const_cast<Wall&> (**wit);
                size_t c1id = w.cell1()->index();
                size_t c2id = w.cell2()->index();
                int j = -1;
								if (c_id == c1id)
									j = i;
								else if (c_id == c2id)
									j = i+1;
								else
									{
										std::cerr << "VTUostream::write_wall_data2(Tissue const& t); Wall does not report connection to the cell which was accessed through\n";
										exit(EXIT_FAILURE);
									}
							}
						*m_os << w.variable(j) << " ";
					}
			}
		*m_os << "\n"
					<< "</DataArray>\n";
}

//-----------------------------------------------------------------------------
void VTUostream::write_cell_point_geometry2(Tissue const& t)
{
    typedef std::vector<Cell>::const_iterator CellIter;
    std::vector<Cell> const& cells = t.cell();

    *m_os << "<Points>\n"
            << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    CellIter cit, cend;
    for (cit = cells.begin(), cend = cells.end(); cit != cend; ++cit)
    {
        Cell &c = const_cast<Cell&> (*cit);
        typedef std::vector<Wall*>::const_iterator WallIter;
        std::vector<Wall*> const& walls = c.wall();

        std::vector<double> cent(3);
        cent = c.positionFromVertex();
        Point center(cent[0], cent[1], cent[2]);

        typedef std::vector<Vertex*>::const_iterator CVIter;
        std::vector<Vertex*> const& verts = c.vertex();
        CVIter cviter, cvend;
        for (cviter = verts.begin(), cvend = verts.end(); cviter != cvend; ++cviter)
        {
            Point p = Point(**cviter).displace_towards(center, D);
            *m_os << p.x << " " << p.y << " " << p.z << "\n";
        }
    }
    *m_os << "</DataArray>\n"
            << "</Points>\n";
}
//-----------------------------------------------------------------------------

void VTUostream::write_cell_geometry2(Tissue const& t)
{
    typedef std::vector<Cell>::const_iterator CellIter;
    std::vector<Cell> const& cells = t.cell();

    *m_os << "<Cells>\n"
            << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    CellIter cit, cend;
    int count = 0;
    for (cit = cells.begin(), cend = cells.end(); cit != cend; ++cit)
    {
        int nverts = cit->numVertex();
        for (int i = 0; i < nverts; ++i, ++count)
        {
            *m_os << count << " ";
        }
        *m_os << "\n";
    }

    *m_os << "</DataArray>\n"
            << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    int total_offset = 0;
    for (cit = cells.begin(), cend = cells.end(); cit != cend; ++cit)
    {
        total_offset += cit->numVertex();
        *m_os << total_offset << " ";
    }
    *m_os << "\n"

            << "</DataArray>\n"
            << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (int i = 0; i < t.numCell(); ++i)
    {
        *m_os << 7 << " ";
    }
    *m_os << "\n";
    *m_os << "</DataArray>\n"
            << "</Cells>\n";
}
//-----------------------------------------------------------------------------

void VTUostream::write_wall_point_geometry2(Tissue const& t, std::vector<Vertex*> &verts)
{
    //write all the vertices first in the order of their indices
    typedef std::vector<Vertex>::const_iterator VertexIter;
    std::vector<Vertex> const& vertices = t.vertex();

    VertexIter viter, vend;
    *m_os << "<Points>\n"
            << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (viter = vertices.begin(), vend = vertices.end(); viter != vend; ++viter)
    {
        std::vector<double> const& p = viter->position();
        size_t npos = viter->numPosition();
        for (size_t j = 0; j < npos; ++j)
            *m_os << p[j] << " ";
        if (npos < 3)
            *m_os << "0 ";
        *m_os << "\n";
    }

    //write the displaced vertices of each cell
    typedef std::vector<Cell>::const_iterator CellIter;
    std::vector<Cell> const& cells = t.cell();

    CellIter cit, cend;
    for (cit = cells.begin(), cend = cells.end(); cit != cend; ++cit)
    {
        Cell &c = const_cast<Cell&> (*cit);
        typedef std::vector<Wall*>::const_iterator WallIter;
        std::vector<Wall*> const& walls = c.wall();

        WallIter wit, wend;
        wit = walls.begin();
        Vertex *w1v1 = (*wit)->vertex1(), *w1v2 = (*wit)->vertex2();
        ++wit;
        Vertex *w2v1 = (*wit)->vertex1(), *w2v2 = (*wit)->vertex2();
        if (w1v1 == w2v1) //w1v2, (w1v1 = w2v1), w2v2
        {
            verts.push_back(w1v2);
            verts.push_back(w1v1);
            w1v2 = w2v2;
        }
        else if (w1v1 == w2v2) //w1v2, (w1v1 = w2v2), w2v1
        {
            verts.push_back(w1v2);
            verts.push_back(w1v1);
            w1v2 = w2v1;
        }
        else if (w1v2 == w2v1) //w1v1, (w1v2 = w2v1), w2v2
        {
            verts.push_back(w1v1);
            verts.push_back(w1v2);
            w1v2 = w2v2;
        }
        else if (w1v2 == w2v2) //w1v1, (w1v2 = w2v2), w2v1
        {
            verts.push_back(w1v1);
            verts.push_back(w1v2);
            w1v2 = w2v1;
        }
        else
        {
            std::cerr << "VTUostream::write_wall_point_geometry(Tissue const& t); Consecutive cell walls do not share a vertex.\n";
            exit(EXIT_FAILURE);
        }
        ++wit;
        for (wend = walls.end(); wit != wend; ++wit)
        {
            w2v1 = (*wit)->vertex1();
            w2v2 = (*wit)->vertex2();
            if (w1v2 == w2v1)
            {
                verts.push_back(w2v1);
                w1v2 = w2v2;
            }
            else if (w1v2 == w2v2)
            {
                verts.push_back(w2v2);
                w1v2 = w2v1;
            }
            else
            {
                std::cerr << "VTUostream::write_wall_point_geometry(Tissue const& t); Consecutive cell walls do not share a vertex.\n";
                exit(EXIT_FAILURE);
            }
        }
    }
    typedef std::vector<Vertex*>::const_iterator VertexPIter;
    VertexPIter vpit = verts.begin(), vpend = verts.end();
    for (cit = cells.begin(), cend = cells.end(); cit != cend; ++cit)
    {
        Cell &c = const_cast<Cell&> (*cit);
        std::vector<double> cent(3);
        cent = c.positionFromVertex();
        Point center(cent[0], cent[1], cent[2]);

        int nwall = cit->numWall();
        for (int i = 0; i < nwall; ++i)
        {
            Point p = Point(**vpit++).displace_towards(center, D);
            *m_os << p.x << " " << p.y << " " << p.z << "\n";
        }

    }
    *m_os << "</DataArray>\n"
            << "</Points>\n";
}
//-----------------------------------------------------------------------------

void VTUostream::write_wall_geometry2(Tissue const& t, std::vector<Vertex*> const& verts)
{
    typedef std::vector<Cell>::const_iterator CellIter;
    std::vector<Cell> const& cells = t.cell();
    *m_os << "<Cells>\n"
            << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    int count = 0, offset = t.numVertex();
    //construct quad walls by ordering vertices circularly: first 2 original verices of the wall then 2 displaced vertices
    std::vector<Vertex*> ::const_iterator vit = verts.begin(), vend = verts.end(), vstart;
    CellIter cit, cend;
    for (cit = cells.begin(), cend = cells.end(); cit != cend; ++cit)
    {
        vstart = vit;
        int nwall = cit->numWall() - 1, temp = offset + count;
        for (int i = 0; i < nwall; ++i)
        {
            *m_os << (*vit++)->index() << " ";
            *m_os << (*vit)->index() << " " << temp + i + 1 << " " << temp + i << " ";
        }
        *m_os << (*vit++)->index() << " " << (*vstart)->index() << " " << temp << " " << temp + nwall << " ";
        count += cit->numWall();
    }
    *m_os << "\n";
    *m_os << "</DataArray>\n"
            << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    int total_offset = 0;
    for (int i = 0; i < count; ++i)
    {
        total_offset += 4;
        *m_os << total_offset << " ";
    }
    *m_os << "\n"

            << "</DataArray>\n"
            << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (int i = 0; i < count; ++i)
    {
        *m_os << 7 << " ";
    }
    *m_os << "\n";
    //*/
    *m_os << "</DataArray>\n"
            << "</Cells>\n";
}
//-----------------------------------------------------------------------------
