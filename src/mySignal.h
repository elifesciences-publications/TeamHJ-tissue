#ifndef MYSIGNAL_H
#define MYSIGNAL_H

#include <vector>
#include <signal.h>
#include "baseSolver.h"

namespace mySignal {
  void myExit();
	void addSolver(BaseSolver *S);
  void solverSignalHandler(int signal);
}

#endif /* MYSIGNAL_H */
