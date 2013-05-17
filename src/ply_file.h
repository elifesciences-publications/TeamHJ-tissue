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
// class PLY_reader;
//-----------------------------------------------------------------------------
// example of usage:
// Tissue T;
// //... initialize T ...
// PLY_file plyfile("output.ply");
// plyfile << T;
/// @brief Handles a ply file for writing a Tissue.
class PLY_file
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
    /// @brief Outputs current Tissue
    void write(Tissue const& t);
    friend class PLY_reader;
};
//-----------------------------------------------------------------------------
// class PLY_parser
// {
// public:
//     typedef boost::char_separator<char> Separator;
//     typedef boost::tokenizer<Separator> Tokenizer;
//     typedef Tokenizer::iterator Token_iterator;
//     /// @brief Constructor
//     PLY_parser(PLY_file const& f);
//     /// @brief Extracts information from header
//     void parse_header();
// private:
//     /// @brief PLY input file stream associated with parser
//     std::ifstream m_input;
//     /// @brief
//     int m_num_vert;
//     /// @brief
//     int m_num_wall;
//     /// @brief
//     int m_num_cell;
//     
//     Separator m_sep;
//     Tokenizer m_tok;
// };
//-----------------------------------------------------------------------------
#endif
