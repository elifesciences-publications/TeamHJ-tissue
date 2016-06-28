#ifndef MYFILES_H
#define MYFILES_H

#include <iostream>
#include <string>

namespace myFiles {
  std::istream *openFile(const std::string &fileName);
  void createDir(const char path);
}

#endif /* MYFILES_H */
