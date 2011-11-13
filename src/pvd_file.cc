//
// File:    pvd_file.cc
// Author:  pkrupinski
// Created: September 15, 2011, 3:20 PM
//
#include "pvd_file.h"
#include "VTUostream.h"
#include <fstream>
#include <sstream>
#include <stdio.h>
//----------------------------------------------------------------------------

PVD_file::PVD_file() : cell_out_fun_ptr(NULL)
{
}
//----------------------------------------------------------------------------

PVD_file::PVD_file(const std::string filename) : filename(filename), vtu_basenames(1), vtu_filenames(1), cell_out_fun_ptr(&VTUostream::write_cells)
{
    vtu_basenames[0] = filename;
}
//----------------------------------------------------------------------------

PVD_file::PVD_file(const std::string filename, const std::string vtu_filename1, const std::string vtu_filename2) : filename(filename), m_counter(0), vtu_basenames(2), vtu_filenames(2), cell_out_fun_ptr(&VTUostream::write_cells2)
{
    vtu_basenames[0] = vtu_filename1;
    vtu_basenames[1] = vtu_filename2;
}
//----------------------------------------------------------------------------

//PVD_file::PVD_file(const std::string filename, const std::string vtu_filename1, const std::string vtu_filename2, size_t n) : filename(filename), m_counter(0), vtu_basenames(2), cell_out_fun_ptr(NULL)
//{
//    vtu_basenames[0] = vtu_filename1;
//    vtu_basenames[1] = vtu_filename2;
//    pvdFileWriteFull(n);
//}
//----------------------------------------------------------------------------

PVD_file::~PVD_file()
{
}
//----------------------------------------------------------------------------
PVD_file& PVD_file::open(const std::string filename, const std::string vtu_filename1, const std::string vtu_filename2)
{
    this->filename = filename;
    m_counter = 0;
    vtu_basenames.resize(2);
    vtu_filenames.resize(2);
    vtu_basenames[0] = vtu_filename1;
    vtu_basenames[1] = vtu_filename2;   
}
//----------------------------------------------------------------------------
void PVD_file::operator<<(Tissue const& t)
{
    write(t, -1.0);
}
//----------------------------------------------------------------------------

void PVD_file::write(Tissue const& t, double time)
{
    // Update vtu file name and clear file
    vtuNameUpdate(m_counter);
    // Write pvd file
    pvdFileWrite(m_counter, time);
    std::ofstream co(vtu_filenames[0].c_str());
    VTUostream out(co);
    //call approporiate cell output function depending on the mode defined by the choice of constructor
    (out.*cell_out_fun_ptr)(t);
    out.close();
    if (vtu_filenames.size() > 1)
    {
        std::ofstream wo(vtu_filenames[1].c_str());
        out.open(wo);
        out.write_walls2(t);
        out.close();
    }
    // Increase the number of times we have saved the mesh
    m_counter++;
}
//----------------------------------------------------------------------------

void PVD_file::write(Tissue const& t, const std::string vtu_filename1, const std::string vtu_filename2, size_t count)
{
    std::vector<std::string> basenames(2), filenames(2);
    basenames[0] = vtu_filename1;
    basenames[1] = vtu_filename2;
    // Update vtu file name
    vtuNameUpdate(count, basenames, filenames);
    std::ofstream co(filenames[0].c_str());
    VTUostream out(co);
    out.write_cells2(t);
    out.close();
    std::ofstream wo(filenames[1].c_str());
    out.open(wo);
    out.write_walls2(t);
    out.close();
}
//----------------------------------------------------------------------------
void PVD_file::writeTwoWall(Tissue const& t, const std::string vtu_filename1, const std::string vtu_filename2, size_t count)
{
    std::vector<std::string> basenames(2), filenames(2);
    basenames[0] = vtu_filename1;
    basenames[1] = vtu_filename2;
    // Update vtu file name
    vtuNameUpdate(count, basenames, filenames);
    std::ofstream co(filenames[0].c_str());
    VTUostream out(co);
    out.write_cells2(t);
    out.close();
    std::ofstream wo(filenames[1].c_str());
    out.open(wo);
//    std::cout << "write_walls3\n";
    out.write_walls3(t);
    out.close();
}
//----------------------------------------------------------------------------

void PVD_file::pvdFileWrite(size_t num, double time)
{
    std::fstream pvdFile;

    if (num == 0)
    {
        // Open pvd file
        pvdFile.open(filename.c_str(), std::ios::out | std::ios::trunc);
        // Write header
        pvdFile << "<?xml version=\"1.0\"?> " << std::endl;
        pvdFile << "<VTKFile type=\"Collection\" version=\"0.1\" > " << std::endl;
        pvdFile << "<Collection> " << std::endl;
    }
    else
    {
        // Open pvd file
        pvdFile.open(filename.c_str(), std::ios::out | std::ios::in);
        pvdFile.seekp(mark);
    }
    if (time < 0.0)
        time = num;
    for (int i = 0; i < vtu_filenames.size(); ++i)
    {
        // Remove directory path from name for pvd file
        std::string fname;
        fname.assign(vtu_filenames[i], vtu_filenames[i].find_last_of("/") + 1, vtu_filenames[i].size());
        // Data file name
        pvdFile << "<DataSet timestep=\"" << time << "\" part=\"" << i << "\" file=\"" << fname << "\"/>" << std::endl;
    }
    mark = pvdFile.tellp();

    // Close headers
    pvdFile << "</Collection> " << std::endl;
    pvdFile << "</VTKFile> " << std::endl;

    // Close file
    pvdFile.close();
}
//----------------------------------------------------------------------------
void PVD_file::write_full_pvd(const std::string filename, const std::string vtu_filename1, const std::string vtu_filename2, size_t n)
{
    std::vector<std::string> basenames(2);
    basenames[0] = vtu_filename1;
    basenames[1] = vtu_filename2;
    std::fstream pvdFile;
    size_t n_files = basenames.size();
    std::vector<std::string> filestart(n_files), extension(n_files), fname(n_files);

    for (int i = 0; i < n_files; ++i)
    {
        std::string fname = basenames[i];
        filestart[i].assign(fname, 0, fname.find("."));
        extension[i].assign(fname, fname.find("."), fname.size());
    }
    // Open pvd file
    //    pvdFile.open(filename.c_str(), std::ios::out | std::ios::trunc);
    pvdFile.open(filename.c_str(), std::ios::out);
    // Write header
    pvdFile << "<?xml version=\"1.0\"?> " << std::endl;
    pvdFile << "<VTKFile type=\"Collection\" version=\"0.1\" > " << std::endl;
    pvdFile << "<Collection> " << std::endl;

//    std::cout << "num = " << num << "\n";
    for (int i = 0; i < n; ++i)
    {
//        std::cout << "i = " << i << "\n";
        std::ostringstream fileid;
        fileid.fill('0');
        fileid.width(6);
        fileid << i;
        for (int j = 0; j < n_files; ++j)
        {
            std::string fname = filestart[j] + fileid.str() + extension[j];
//            std::cout << "j = " << j << "; fname = \'" << fname << "\'\n";
            // Data file name
            pvdFile << "<DataSet timestep=\"" << i << "\" part=\"" << j << "\" file=\"" << fname << "\"/>" << std::endl;
        }
    }
    // Close headers
    pvdFile << "</Collection> " << std::endl;
    pvdFile << "</VTKFile> " << std::endl;

    // Close file
    pvdFile.close();
}

void PVD_file::pvdFileWriteFull(size_t num)
{
    std::fstream pvdFile;
    size_t n_files = vtu_basenames.size();
    std::vector<std::string> filestart(n_files), extension(n_files), fname(n_files);

    for (int i = 0; i < n_files; ++i)
    {
        std::string fname = vtu_basenames[i];
        filestart[i].assign(fname, 0, fname.find("."));
        extension[i].assign(fname, fname.find("."), fname.size());
    }
    // Open pvd file
    //    pvdFile.open(filename.c_str(), std::ios::out | std::ios::trunc);
    pvdFile.open(filename.c_str(), std::ios::out);
    // Write header
    pvdFile << "<?xml version=\"1.0\"?> " << std::endl;
    pvdFile << "<VTKFile type=\"Collection\" version=\"0.1\" > " << std::endl;
    pvdFile << "<Collection> " << std::endl;

//    std::cout << "num = " << num << "\n";
    for (int i = 0; i < num; ++i)
    {
//        std::cout << "i = " << i << "\n";
        std::ostringstream fileid;
        fileid.fill('0');
        fileid.width(6);
        fileid << i;
        for (int j = 0; j < n_files; ++j)
        {
            std::string fname = filestart[j] + fileid.str() + extension[j];
//            std::cout << "j = " << j << "; fname = \'" << fname << "\'\n";
            // Data file name
            pvdFile << "<DataSet timestep=\"" << i << "\" part=\"" << j << "\" file=\"" << fname << "\"/>" << std::endl;
        }
    }
    // Close headers
    pvdFile << "</Collection> " << std::endl;
    pvdFile << "</VTKFile> " << std::endl;

    // Close file
    pvdFile.close();
}
//----------------------------------------------------------------------------

void PVD_file::vtuNameUpdate(const int m_counter, std::vector<std::string>const& basenames, std::vector<std::string>& filenames)
{
    std::string filestart, extension;
    std::ostringstream fileid;
    fileid.fill('0');
    fileid.width(6);
    fileid << m_counter;

    for (int i = 0; i < basenames.size(); ++i)
    {
        std::ostringstream newfilename;
        std::string fname = basenames[i];
        filestart.assign(fname, 0, fname.find("."));
        extension.assign(fname, fname.find("."), fname.size());
        newfilename << filestart << fileid.str() << extension;
        filenames[i] = newfilename.str();
        //       // Make sure file is empty
        //       FILE* fp = fopen(vtu_filenames[i].c_str(), "w");
        //       fclose(fp);
    }
}
//-----------------------------------------------------------------------------
