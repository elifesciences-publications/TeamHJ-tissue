//
// Filename     : myRandom.cc
// Description  : Defining the functions for random number generators
// Author(s)    : Henrik Jonsson (henrik@thep.lu.se)
//              : Bo Soderberg (bosse@thep.lu.se)
// Created      : 
// Revision     : $Id: myRandom.cc 215 2007-03-30 10:29:03Z henrik $
//

#include <iostream>

#include "myRandom.h"

namespace myRandom {

	double Rnd( void )
	{
		return ran3();
		//return random() / double(INT_MAX);
	}
	
	double Grand( void )
	{
		return sqrt( -2 * log(Rnd())) * cos( M_PI * Rnd() );
	}
		
	double Rndstep(int step)
	{
		return (double)((int)((step+1)*Rnd()))/step;
	}
	
	long int Randomize( void ) 
	{
		return ran3Randomize();
	}
		
#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)
//double ran3(long *idum)
	static long idum = 1;
	double ran3( void )
	{
		static int inext,inextp;
		static long ma[56];
		static int iff=0;
		long mj,mk;
		int i,ii,k;
		double ret_val;
		
		//if (*idum < 0 || iff == 0) {
		if (myRandom::idum < 0 || iff == 0) {
			iff=1;
			//mj=MSEED-(*idum < 0 ? -*idum : *idum);
			mj=MSEED-(myRandom::idum < 0 ? -myRandom::idum : myRandom::idum);
			mj %= MBIG;
			ma[55]=mj;
			mk=1;
			for (i=1;i<=54;i++) {
				ii=(21*i) % 55;
				ma[ii]=mk;
				mk=mj-mk;
				if (mk < MZ) mk += MBIG;
				mj=ma[ii];
			}
			for (k=1;k<=4;k++)
				for (i=1;i<=55;i++) {
					ma[i] -= ma[1+(i+30) % 55];
					if (ma[i] < MZ) ma[i] += MBIG;
				}
			inext=0;
			inextp=31;
			//*idum=1;
			myRandom::idum=1;
		}
		if (++inext == 56) inext=1;
		if (++inextp == 56) inextp=1;
		mj=ma[inext]-ma[inextp];
		if (mj < MZ) mj += MBIG;
		ma[inext]=mj;
		ret_val = mj*FAC;
		if (mj == 0) ret_val = FAC;
		return ret_val;
	}
#undef MBIG
#undef MSEED
#undef MZ
#undef FAC
	
	void sran3( long idumVal )
	{
		std::cerr << "sran3()\n";
		// Make sure seed in [1:10000000] 
		// The ran3 function has proven to be instable otherwise
		if( idumVal<1 || idumVal>10000000 ) {
			std::cerr << "myRandom::sran3() Seed (idum) provided in dangerous "
								<< "region. Usevalues in [1:10000000]." << std::endl;
			exit(-1);
		}
		idum=-idumVal;
	}
	
	long int ran3Randomize( void ) 
	{
		std::cerr << "ran3Randomize()\n";
		struct timeval tp; // { long tv_sec; long tv_usec; }
		long int seed;
		
		gettimeofday(&tp, 0);
		//seed = 1024 * tp.tv_sec + (int) (.001024 * tp.tv_usec)
		seed = 1+tp.tv_usec%9999999;//Put seed in [1:10000000]
		idum = -seed;
		return -seed;
	}
} //end namespace myRandom

