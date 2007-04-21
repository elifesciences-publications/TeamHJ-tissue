#include "myFiles.h"

#include <fstream>
#include <sstream>

std::istream *myFiles::openFile(const std::string &file)
{
	std::ifstream IN(file.c_str());
	if (!IN)
		return NULL;
	
	std::ostringstream sout;
	std::string tmp;
	std::string::size_type pos;
	
	while (getline(IN, tmp)) {
		pos = tmp.find("#", 0);
		if (pos == std::string::npos)
			sout << tmp << std::endl;
		else
			sout << tmp.erase(pos) << std::endl;
	}
	IN.close();
	return new std::istringstream(sout.str());
}
