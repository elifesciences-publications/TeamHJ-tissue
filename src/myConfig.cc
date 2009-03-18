//
// Filename     : myConfig.cc
// Description  : Singleton handling user configuration
// Author(s)    : Patrik Sahlin (sahlin@thep.lu.se)
// Created      : March 2007
// Revision     : $Id: myConfig.cc 217 2007-03-30 11:53:16Z henrik $
//
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>

#include "myConfig.h"
#include "myFiles.h"

std::map< std::string, std::pair<std::vector<std::string>, bool> > myConfig::configure_;
std::vector<std::string> myConfig::argv_;
std::vector<ConfigOption> myConfig::options_;
myConfig *myConfig::instance_ = NULL;

myConfig::myConfig(int argc, char *argv[], const std::string &file)
{
	parseFile(file);
	parseCommandLine(argc, argv);
}

std::vector<std::string> myConfig::initConfig(int argc, char *argv[], 
																							const std::string &file)
{
	if (instance_ == NULL)
		instance_ = new myConfig(argc, argv, file);
	return argv_;
}

void myConfig::registerOption(const std::string &key, size_t arguments)
{
	// Add option to options_ vector.
	ConfigOption option = { key, arguments };
	options_.push_back(option);

	// Initialise configure_ entry.
	configure_[key].second = false;
	std::vector<std::string> tmpVector(arguments);
	configure_[key].first =  tmpVector;
}

std::string myConfig::getValue(const std::string &key, size_t index)
{
	std::map< std::string, std::pair<std::vector<std::string>, bool> >::iterator i;
	i = configure_.find(key);
	if (i == configure_.end()) {
		std::cerr << "Internal Error: myConfig::getValue() - key: " << key
							<< " not defined." << std::endl;
		exit(EXIT_FAILURE);
	}
	if (index > i->second.first.size() - 1) {
		std::cerr << "Internal Error: myConfig::getValue() - Index out of bounds." 
							<< std::endl;
		exit(EXIT_FAILURE);
	}
	return configure_[key].first[index];
}

bool myConfig::getBooleanValue(const std::string &key)
{
	std::map< std::string, std::pair<std::vector<std::string>, bool> >::iterator i;
	i = configure_.find(key);
	if (i == configure_.end()) {
		std::cerr << "Error: myConfig::getValue() - key: " << key
							<< " not defined." << std::endl;
		exit(EXIT_FAILURE);
	}
	return configure_[key].second;
}

int myConfig::argc(void)
{
	return argv_.size();
}

std::string myConfig::argv(int index)
{
	return argv_[index];
}

void myConfig::parseFile(const std::string &file)
{
	std::istream *IN = myFiles::openFile(file);
	if (IN != NULL) {
		while (!IN->eof()) {
			// Read line from file.
			std::string tmpLine;
			getline(*IN, tmpLine);
			
			// Return if EOF is reached.
			if (IN->eof())
				break;

			// Make temporary streamstream of line read from file.
			std::stringstream LINE;
			LINE.str(tmpLine);

			// Read key name.
			std::string tmpKey;
			LINE >> tmpKey;
			if (tmpKey.empty())
				continue;

			// Search for option with name equal to key read from file.
			std::vector<ConfigOption>::iterator iter;
			for (iter = options_.begin(); iter != options_.end(); ++iter)
				if (iter->key.compare(tmpKey) == 0)
					break;

			// Warn the user if key is not registered.
			if (iter == options_.end()) {
				std::cerr << "Warning: Configuration option " << tmpKey 
						<< " is unsupported." << std::endl; 
				continue;
			}

			// Parse arguments if option requires a set of arguments.
			// Mark option's boolean value to be set.
			if (iter->indices > 0) {
				for (size_t i = 0; i < iter->indices; ++i) {
					std::string tmpValue;
					LINE >> tmpValue;
					// Make safety check.
					if (tmpValue.empty()) {
						std::cerr << "Error: Option " << iter->key << " requires "
											<< iter->indices << " arguments." << std::endl;
						exit(EXIT_FAILURE);
					}

					configure_[tmpKey].first[i] = tmpValue;
					configure_[tmpKey].second = true;
				}
			} else if (iter->indices == 0)
				configure_[tmpKey].second = true;
		}
	}
	delete IN;
}

void myConfig::parseCommandLine(int argc, char *argv[])
{
	// Clear argv_ vector.
	argv_.clear();

	// Loop over command line arguments.
	for (size_t i = 0; i < (size_t) argc; ++i) {
		// If an argument starts with character '-' it is considered to be
		// an option and the rest of the argument is used as the key.  If
		// it is not an option the argument is pushed to the back of the
		// argv_ vector.
		if (argv[i][0] != '-')
			argv_.push_back(argv[i]);
		else {
			// Search for option with name equal to key parsed from command line.
			std::vector<ConfigOption>::iterator iter;
			for (iter = options_.begin(); iter != options_.end(); ++iter)
				if (iter->key.compare(argv[i] + 1) == 0)
					break;

			// Warn the user if key is not registered.
			if (iter == options_.end()) {
				std::cerr << "Warning: Configuration option " << argv[i] + 1
									<< " is unsupported." << std::endl; 
				continue;
			}

			// Parse arguments if option requires a set of arguments.
			// Mark option's boolean value to be set.
			if (iter->indices > 0) {
				for (size_t j = 0; j < iter->indices; ++j) {
					// Make safety check.
					if (i + 1 + j > (size_t) (argc - 1)) {
						std::cerr << "Error: Option -" << iter->key << " requires "
											<< iter->indices << " arguments." << std::endl;
						exit(EXIT_FAILURE);
					}
					configure_[iter->key].first[j] = argv[i + 1 + j];
					configure_[iter->key].second = true;
				}
			} else {
				configure_[iter->key].second = true;
			}
			i += iter->indices;
		}
	}

}

