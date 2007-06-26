//
// Filename     : myTimes.h
// Description  : Defining namespace for timing functions
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
// Created      : 
// Revision     : $Id: myTimes.h 189 2007-02-23 01:18:49Z henrik $
//
#ifndef MYTIMES_H
#define MYTIMES_H

#include <stdio.h>
#include <sys/times.h>
#include <unistd.h>

///
/// @brief Namespace with functions for measuring time aspects
///
namespace myTimes {
  
  long clocktick();
  
  ///
  /// @brief Returns user time since start.
  ///
  double getTime(void);
  
  ///
  /// @brief Returns system time since start
  ///
  double getSysTime(void);
  
  ///
  /// @brief Prints result from calling getTime() and getSysTime()
  ///
  void showTime(void);
  
  ///
  /// @brief Returns user time since last call to getDiffTime()
  ///
  double getDiffTime(void);
  
}
#endif
