#ifndef MYRANDOM_H
#define MYRANDOM_H
//
// Filename     : myRandom.h
// Description  : Defining namespace for random functions
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
//              : Bo Soderberg (bosse@thep.lu.se)
// Created      : 
// Revision     : $Id: myRandom.h 215 2007-03-30 10:29:03Z henrik $
//

#include <climits>
#include <cmath>
#include <cstdlib>
#include <sys/time.h>

///
/// @brief Namespace including functions for generating random numbers
///
/// This namespace includes several functions for generating random
/// numbers of various kinds. It also includes Randomize functions for
/// seeding the random number generators. Some of these functions are
/// inherited from Bo Soderberg.
///
namespace myRandom {
	
	///
	/// @brief generates a random number in [0:1) by using ran3()
	///
	/// This function is provided for backward compability.
	/// 
	/// @see ran3()
	///
	double Rnd( void );
	
	///
	/// @brief Returns a Gaussian random number using Rnd()
	///
	/// @see Rnd()
	///
	double Grand( void );
	
	///
	/// @brief Calls the ran3Randomize() function
	///
	/// This function is provided for backward compability.
	/// 
	/// @see ran3Randomize()
	///
	long int Randomize( void );
	
	///
	/// @brief (0,1) random number generator
	///
	/// This random number generator returns a number in [0:1). It is
	/// adopted from Numerical Recipes and is credited to Knuth. This
	/// function is supposed to create better random properties compared
	/// to the standard random().
	///
	double ran3( void );
	//double ran3(long *idum);

	///
	/// @brief Generatess an idum based on time for the ran3() function
	///
	long int ran3Randomize( void ); 

	///
	/// @brief Sets the idum variable for the ran3 random sequence
	///
	void sran3(long idumVal);

	///
	/// @brief Generates a random number between 0 and 1 with the
	/// resolution of step
	///
	/// The step sets the number of possible levels.
	///
	double Rndstep(int step);
}

#endif

