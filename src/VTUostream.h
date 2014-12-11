//
// File:    VTUostream.h
// Author:  pkrupinski
// Created: September 15, 2011, 3:20 PM
//
#ifndef _VTUOSTREAM_H_
#define	_VTUOSTREAM_H_

#include <iostream>
#include <vector>
#include <map>
class Tissue;
class Vertex;
//-----------------------------------------------------------------------------
namespace IO
{
  class Vector;
  /// @brief 3D point
  class Point
  {
  public:
    Point():x(0.0), y(0.0), z(0.0) {}
    Point(double x, double y):x(x), y(y), z(0.0) {}
    Point(double x, double y, double z):x(x), y(y), z(z) {}
    Point(Vertex const& v);

    double sqr_length_to(Point const& p) const
    {
      double d = p.x-x;
      double res = d*d;
      d = p.y-y;
      res += d*d;
      d = p.z-z;
      res += d*d;
      return res;
    }

    Point const& displace_towards(Point const&p, double d)
    {
      x += d*(p.x - x);
      y += d*(p.y - y);
      z += d*(p.z - z);
      return *this;
    }

    Point const& displace(Vector const&v, double d);

    double x,y,z;
  };
  //-----------------------------------------------------------------------------
  /// @brief 3D vector
  class Vector
  {
  public:
    Vector():x(0.0), y(0.0), z(0.0) {}
    Vector(double x, double y):x(x), y(y), z(0.0) {}
    Vector(double x, double y, double z):x(x), y(y), z(z) {}
    Vector(Point const& p1, Point const& p2):x(p2.x-p1.x), y(p2.y-p1.y), z(p2.z-p1.z) {}

    double sqr_length() const
    {
      return x*x+y*y+z*z;
    }
    Vector operator+(const Vector& v)
    {
      return Vector(x+v.x,y+v.y,z+v.z);
    }
    double x,y,z;
  };
}
//-----------------------------------------------------------------------------
/// @brief Output stream for VTU format. Supports cell and wall output separately
class VTUostream
{
public:

  enum Cell_type{ POLYGON = 7, TRIANGLE = 5};
  
  VTUostream();
  VTUostream(std::ostream& o);
  ~VTUostream() {close();}

  void open(std::ostream& o)
  {
    m_os = &o;
    header();
  }

  void close()
  {
    if (m_os)
      footer();
    m_os = 0;
  }

  typedef const void* Const_void_ptr;

  operator Const_void_ptr() const
  {
    if (m_os)
      return m_os;
    return NULL;
  }

  std::ostream& os()
  {
    //        assert(m_os);
    return *m_os;
  }

  VTUostream& seekp(std::streampos pos)
  {
    m_os->seekp(pos);
    return *this;
  }

  std::streampos tellp()
  {
    return m_os->tellp();
  }

  std::streamsize width(std::streamsize wide)
  {
    return m_os->width(wide);
  }
  /// @brief Write cells using geometry directly from tissue without making room for the walls display
  void write_cells(Tissue const& t);
  /// @brief Write cells with shrinked geometry leaving space for walls display
  void write_cells2(Tissue const& t);
  /// @brief Write cells with shrinked geometry leaving space only for outer walls display
  void write_cells3(Tissue const& t);
  /// @brief Write walls using lines
  void write_walls(Tissue const& t);
  /// @brief Write walls using 2D elements and assuming single wall between cells
  void write_walls2(Tissue const& t);
  /// @brief Write walls using 2D elements and assuming double wall between cells
  void write_walls3(Tissue const& t);
  /// @brief Write inner walls only using lines and assuming double wall between cells
  void write_inner_walls(Tissue const& t);
  /// @brief Write outer walls using 2D elements and assuming double wall between cells
  void write_outer_walls(Tissue const& t);
protected:
  // cell and wall geometry for walls as line segments
  void write_cell_point_geometry(Tissue const& t);
  void write_cell_geometry(Tissue const& t, Cell_type ct = POLYGON);
  void write_wall_point_geometry(Tissue const& t);
  void write_wall_geometry(Tissue const& t);
  // cell and wall geometry for 2D walls
  void write_cell_point_geometry2(Tissue const& t);
  void write_cell_geometry2(Tissue const& t, Cell_type ct = POLYGON);
  void write_wall_point_geometry2(Tissue const& t, std::vector<Vertex*> & verts);
  void write_wall_geometry2(Tissue const& t, std::vector<Vertex*> const& verts);
  // cell and wall geometry for 2D walls printing inner and outer cell walls separately based on the flag in the last wall variable
  void write_cell_point_geometry3(Tissue const& t, std::vector<IO::Point>& disp_points, std::vector<char>& vertex_flag, std::vector< std::map<size_t,size_t> >& cvp_map, size_t flag_pos);
  void write_cell_geometry3(Tissue const& t, Cell_type ct = POLYGON);
  size_t prepare_marked_vertices( Tissue const& t, std::vector<char>& vertex_flag, size_t flag_pos, double flag_val );
  std::pair<size_t, size_t> prepare_wall_point_geometry3(Tissue const& t, std::vector<IO::Point>& disp_points, std::vector<char>& vertex_flag, std::vector< std::map<size_t,size_t> >& cvp_map, size_t flag_pos, double flag_val);
  void write_outer_wall_point_geometry3 ( Tissue const& t, std::vector<IO::Point>& disp_points, std::vector<char>& vertex_flag, std::vector<uint>& index_map ); 
  void write_outer_wall_geometry3 ( Tissue const& t, std::vector<uint>& index_map, std::vector< std::map<size_t,size_t> >& cvp_map, size_t offset, size_t flag_pos, double flag_val );
  void write_inner_wall_point_geometry3 ( Tissue const& t, std::vector<char>& vertex_flag, std::vector<uint>& index_map, size_t counter );
  void write_inner_wall_geometry3( Tissue const& t, std::vector<uint>& index_map, size_t flag_pos, double flag_val );
  void write_piece_header(int n_pts, int n_cell);
  void write_piece_footer();

  void write_point_data_header(std::string s){ *m_os << "<PointData " << s << ">\n";  }
  void write_point_data_footer(){ *m_os << "</PointData>\n"; }
  void write_point_data(Tissue const& t);

  void write_cell_data_header(std::string s){ *m_os << "<CellData " << s << ">\n"; }
  void write_cell_data_footer(){ *m_os << "</CellData>\n"; }
  void write_cell_data(Tissue const& t);

    void write_wall_data_header ( std::string s )
    {
        *m_os << "<CellData " << s << ">\n";
    }
    void write_wall_data_footer()
    {
        *m_os << "</CellData>\n";
    }
    // wall data for single wall between cells
    void write_wall_data ( Tissue const& t );
    /// @brief Writes the data from the wall variables assuming paired structure of composite double wall between cells
    void write_wall_data2 ( Tissue const& t );
    /// @brief Writes the data only from the walls matched by the value of the flag. Variables are assumed to have paired structure of composite double wall between cells
    void write_wall_data3 ( Tissue const& t, size_t flag_pos, double flag_val, bool project_cell_variables = false );

  void header()
  {
    os() << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n"
    << "<UnstructuredGrid>\n";
  }

  void footer()
  {

    os() << "</UnstructuredGrid>\n"
    << "</VTKFile>\n";
  }
  //    std::ios::pos_type mark;
  std::ostream* m_os;
  const double D;
};
//-----------------------------------------------------------------------------
inline VTUostream& operator<<(VTUostream& os, const char* s)
{
  os.os() << s;
  return os;
}
//-----------------------------------------------------------------------------
inline VTUostream& operator<<(VTUostream& os, double d)
{
  os.os() << d;
  return os;
}
//-----------------------------------------------------------------------------
inline VTUostream& operator<<(VTUostream& os, int i)
{
  os.os() << i;
  return os;
}
//-----------------------------------------------------------------------------
#endif	/* VTUOSTREAM_H */
