//
// File:    VTUostream.h
// Author:  pkrupinski
// Created: September 15, 2011, 3:20 PM
//
#ifndef _VTUOSTREAM_H_
#define	_VTUOSTREAM_H_

#include <iostream>
#include <vector>
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

    double x,y,z;
  };
}
//-----------------------------------------------------------------------------
/// @brief Output stream for VTU format. Supports cell and wall output separately
class VTUostream
{
public:

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
      return *m_os;
    return 0;
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

  void write_cells(Tissue const& t);
  void write_cells2(Tissue const& t);
  void write_walls(Tissue const& t);
  void write_walls2(Tissue const& t);
  
protected:
  void write_cell_point_geometry(Tissue const& t);
  void write_cell_geometry(Tissue const& t);
  void write_wall_point_geometry(Tissue const& t);
  void write_wall_geometry(Tissue const& t);

  void write_cell_point_geometry2(Tissue const& t);
  void write_cell_geometry2(Tissue const& t);
  void write_wall_point_geometry2(Tissue const& t, std::vector<Vertex*> & verts);
  void write_wall_geometry2(Tissue const& t, std::vector<Vertex*> const& verts);

  void write_piece_header(int n_pts, int n_cell);
  void write_piece_footer();

  void write_point_data_header(std::string s){ *m_os << "<PointData " << s << ">\n";  }
  void write_point_data_footer(){ *m_os << "</PointData>\n"; }
  void write_point_data(Tissue const& t);

  void write_cell_data_header(std::string s){ *m_os << "<CellData " << s << ">\n"; }
  void write_cell_data_footer(){ *m_os << "</CellData>\n"; }
  void write_cell_data(Tissue const& t);

  void write_wall_data_header(std::string s){ *m_os << "<CellData " << s << ">\n"; }
  void write_wall_data_footer(){ *m_os << "</CellData>\n"; }
  void write_wall_data(Tissue const& t);

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
