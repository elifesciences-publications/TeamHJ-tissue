#include <signal.h>

#include "mySignal.h"
#include "myConfig.h"

namespace mySignal {
  std::vector<BaseSolver *> solvers;
}

void mySignal::myExit()
{
	solverSignalHandler(-1);
}
void mySignal::addSolver(BaseSolver *S)
{
  signal(SIGINT,  solverSignalHandler); // ctrl-c
  signal(SIGHUP,  solverSignalHandler); // Close terminal window
  //  signal(SIGQUIT, solverSignalHandler);
  signal(SIGTERM, solverSignalHandler); // kill  
  solvers.push_back(S); 
}

void mySignal::solverSignalHandler(int signal)
{
  size_t i;
	
  std::cerr << "mySignal::solverSignalHandler - "
						<< "Recieved signal ";
  switch (signal) {
		case SIGINT:
			std::cerr << "SIGINT";
			break;
		case SIGHUP:
			std::cerr << "SIGHUP";
			break;
		case SIGQUIT:
			std::cerr << "SIGQUIT";
			break;
		case SIGTERM:
			std::cerr << "SIGTERM";
			break;
		default:
			std::cerr << "default";	
  }
  std::cerr << std::endl;
	// Bad design. You can add more than one solver, but they will all
	// overwrite the same file.
  for (i=0; i<solvers.size(); i++) {
		std::string fileName;
		fileName = myConfig::getValue("init_output", 0);
		if(!fileName.empty()) {
			std::ofstream OUT(fileName.c_str());
			if (!OUT) {
				std::cerr << "Warning: mySignal::solverSignalHandler() -"
									<< "Cant open file for init output.\n";
			} else {
				solvers[i]->printInit(OUT);					
				OUT.close();					
			}
		}
		if (solvers[i]->debugFlag()) {
			fileName = myConfig::getValue("debug_output", 0);
			if(!fileName.empty()) {
				std::ofstream OUT(fileName.c_str());
				if (!OUT) {
					std::cerr << "Warning: mySignal::solverSignalHandler() -"
										<< "Cant open file for debug output.\n";
				} else {
					solvers[i]->printDebug(OUT);
				}
			}
		}
	}
	exit(EXIT_FAILURE);
}
