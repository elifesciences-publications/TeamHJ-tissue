//
// Filename     : myConfig.h
// Description  : Singleton handling user configuration
// Author(s)    : Patrik Sahlin (sahlin@thep.lu.se)
// Created      : March 2007
// Revision     : $Id: myConfig.h 217 2007-03-30 11:53:16Z henrik $
//

#ifndef MYCONFIG_H
#define MYCONFIG_H

#include <string>
#include <vector>
#include <map>
#include <utility>

typedef struct {
	std::string key;
	size_t indices;
} ConfigOption;

///
/// @brief A singleton class to handle user configuration.
///
/// The myConfig class handles user configuration. User
/// configuration is read from the command line and from the file
/// $HOME/.organism. 
/// The following keys are supported: 
///
/// <ul> 
///
/// <li> @b -init_output file - Set file name for
/// output of final state in init file format.
///
/// <li> @b -debug_output file - Saves the last ten configurations
/// before exiting (note you have to exit with mySignal::myExit() in
/// order for this to work)
///
/// <li> @b -parameter_input file - Reads parameter values from file replacing
/// those values read in the model file.
///
/// <li> @b -verbose flag - Set flag for verbose (flag=1) or silent (0) output
/// mode to stderr.
///
/// </ul>
///
class myConfig
{
 public:
	///
	/// @brief Initializes user configuration.
	///
	/// The function takes argc and argv from main() as arguments and
	/// initializes user configuration. The function first reads the
	/// file given by the third argument (in organism it is set to
	/// $HOME/.organism and will then parse command line
	/// arguments. In this way the user can override configuration
	/// from .organism by command line arguments.
	///
	/// For the command line the function is looking for arguments
	/// following this pattern:
	/// @verbatim
	/// ./binary -key1 arg11 arg12 -key2 arg21 arg22 ...
	/// @endverbatim
	///
	/// The .organism configuration file located at $HOME/.organism
	/// should be constructed according to this pattern:
	/// @verbatim
	/// key1 arg11 arg12 ...
	/// key2 arg21 arg22 ...
	/// ...  ...
	/// @endverbatim
	static std::vector<std::string> initConfig(int argc, char *argv[], const std::string &file);
	static std::vector<std::string> initConfig(int argc, char *argv[]);

	///
	/// @brief Registers an option to the application.
	///
	/// The first argument is the name of the option. The second
	/// arguments sepcifies the option's number of arguments.
	///
	static void registerOption(const std::string &key, size_t arguments);

	///
	/// @brief Get value from user configuration.
	///
	/// Returns the value for the argument given by index belonging to
	/// the option given by key.
	///
	static std::string getValue(const std::string &key, size_t index);

	///
	/// @brief Get boolean value from user configuration.
	///
	/// This function is mainly used for options with no arguments, but
	/// can also be used for options with a non-zero amount of
	/// arguments.
	///
	static bool getBooleanValue(const std::string &key);

	///
	/// @brief Returns number of unparsed arguments from command line.
	///
	static int argc(void);

	///
	/// @brief Returns unparsed argument from command line.
	///
	static std::string argv(int index);

 private:
	// Private constructor and destructor.
	myConfig();
	myConfig(int argc, char *argv[], const std::string &file);
	myConfig(int argc, char *argv[]);

	// Private functions hidden for readability.
	void parseFile(const std::string &file);
	void parseCommandLine(int argc, char *argv[]);

 	// Contains user configuration.
	static std::map< std::string, std::pair<std::vector<std::string>, bool> > configure_;
	// Contains registered options.
	static std::vector<ConfigOption> options_;
	// Contains unparsed arguments from the command line.
	static std::vector<std::string> argv_;

	// The singleton instance.
	static myConfig *instance_;
};

#endif /* MYCONFIG_H */

/*
namespace myConfig {

	///
	/// @brief Initializes user configuration.
	///
	/// The function takes argc and argv from main() as arguments and
	/// initializes user configuration. The function first reads the
	/// file .organism from the user's home directory and will then
	/// parse command line arguments. In this way the user can override
	/// configuration from .organism by command line arguments.
	///
	/// For the command line the function is looking for arguments
	/// following this pattern:
	/// @verbatim
	/// ./binary -key1 value1 -key2 value2 ...
	/// @endverbatim
	///
	/// The .organism configuration file located at $HOME/.organism
	/// should be constructed according to this pattern:
	/// @verbatim
	/// key1 value1
	/// key2 value2
	/// ...  ...
	/// @endverbatim
	///
	/// @return A std::string with unparsed arguments from the command
	/// line.
	///

	///
	/// @brief Gives value from user configuration.
	///
	/// getValue() takes a key as a string as it's argument.
	///
	/// @return The value of given key.
	///
	/// @warning If the key has no value the returned string is
	/// empty. Do a check before using it.  
	///

	/// 
	/// @brief Changes user configuration.
	///
	/// The function is used to change user configuration. The first
	/// argument is the key and the second argument is the value.
	///

}

*/
