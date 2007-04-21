#ifndef MYFILES_H
#define MYFILES_H

#include <iostream>
#include <string>

namespace myFiles {
  std::istream *openFile(const std::string &fileName);
}

#endif /* MYFILES_H */
