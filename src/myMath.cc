/**
 * Filename     : myMath.cc
 * Description  : A set of functions dealing with various mathematical operations.
 * Author(s)    : Patrik Sahlin (sahlin@thep.lu.se)
 * Created      : August 2008
 * Revision     : $Id$
 */

#include <cmath>
#include <cstdlib>
#include <iostream>

#include "myMath.h"

namespace myMath
{
  size_t jacobiTransformation(std::vector< std::vector<double> > A, 
			      std::vector< std::vector<double > > &V, std::vector<double> &d)
  {
    size_t n = A.size();
    
    // Check input and initialize vectors.
	
    for (size_t i = 0; i < n; ++i) {
      if (A[i].size() != n) {
	std::cerr << "jacobiTransformation(): Size of row " << i 
		  << " in matrix A is not equal to number of rows (" << n << ")." << std::endl;
	exit(EXIT_FAILURE);
      }
    }
    
    V.resize(n);
    for (size_t i = 0; i < n; ++i) {
      V[i].resize(n);
      for (size_t j = 0; j < n; ++j) {
	if (i == j) {
	  V[i][i] = 1.0;
	} else {
	  V[i][j] = 0.0;
	}
      }
    }
    
    // Start of iterations.
    
    size_t iterations = 0;
    
    while (true) {
      // Measure sum of off-diagonal elements.
      
      double sum = 0.0;
      for (size_t i = 0; i < n; ++i) {
	for (size_t j = 0; j < n; ++j) {
	  if (i != j) {
	    sum += std::abs(A[i][j]);
	  }
	}
      }
      
      if (sum < 0.00001) {
	break;
      }
      
      // Loop over all off-diagonal elements
      
      for (size_t p = 0; p < n; ++p) {
	for (size_t q = 0; q < n; ++q) {
	  if (p == q) {
	    continue;
	  }
	  
	  if (A[p][q] == 0.0) {
	    continue;
	  }
	  
	  // Get theta and t.
	  
	  double t;
	  
	  double g = 100.0 * std::abs(A[p][q]);
	  double h = A[q][q] - A[p][p];
	  
	  if (std::abs(h) + g == std::abs(h)) {
	    t = A[p][q] / h;
	  } else {
	    double theta = h / (2.0 * A[p][q]);
	    
	    t = (theta > 0) ? +1 : -1;
	    t /= std::abs(theta) + std::sqrt(theta * theta + 1.0);
	  }
	  
	  double c = 1.0 / std::sqrt(t * t + 1.0);
	  double s = t * c;
	  
	  // Calculate the rotation matrix.
	  
	  // Caveat: is n=3 always true?
	  std::vector< std::vector<double> > P(3);
	  for (size_t i = 0; i < 3; ++i) {
	    P[i].resize(3);
	  }
	  for (size_t i = 0; i < n; ++i) {
	    for (size_t j = 0; j < n; ++j) {
	      if (i == j) {
		if (i == p || i == q) {
		  P[i][j] = c;
		} else {
		  P[i][j] = 1.0;
		}
	      } else {
		if (i == p && j == q) {
		  P[i][j] = s;
		} else if (i == q && j == p) {
		  P[i][j] = -s;
		} else {
		  P[i][j] = 0.0;
		}
	      }
	    }
	  }
	  
	  // Calculate Aprime and set A = Aprime.
	  
	  std::vector< std::vector<double> > Aprime(n);
	  
	  for (size_t i = 0; i < n; ++i) {
	    Aprime[i].resize(n);
	    for (size_t j = 0; j < n; ++j) {
	      Aprime[i][j] = 0.0;
	      
	      for (size_t k = 0; k < n; ++k) {
		for (size_t l = 0; l < n; ++l) {
		  Aprime[i][j] += P[k][i] * A[k][l] * P[l][j];
		}
	      }
	    }
	  }
	  
	  A = Aprime;
	  
	  // Calculate Vprime and set V = Vprime
	  
	  std::vector< std::vector<double> > Vprime(n);
	  
	  for (size_t i = 0; i < n; ++i) {
	    Vprime[i].resize(n);
	    for (size_t j = 0; j < n; ++j) {
	      Vprime[i][j] = 0.0;
	      
	      for (size_t k = 0; k < n; ++k) {
		Vprime[i][j] += V[i][k] * P[k][j];
	      }
	    }
	  }
	  
	  V = Vprime;
	}
      }
      
      ++iterations;
      if (iterations>100 && !iterations%100)
	std::cerr << "jacobiTransform iterations: " << iterations << std::endl;
    }
    
    // Fix diagonal output vector.
    
    d.resize(n);
    for (size_t i = 0; i < n; ++i) {
      d[i] = A[i][i];
    }
    
    return iterations;
  }
  
  signed int sign(const double argument)
  {
    return (argument >= 0) ? +1 : -1;
  }

  double areaTriangle(const std::vector<double> &p1,
		      const std::vector<double> &p2,
		      const std::vector<double> &p3)
  {
    double a= std::sqrt( (p1[0]-p2[0])*(p1[0]-p2[0]) +   
			 (p1[1]-p2[1])*(p1[1]-p2[1]) +
			 (p1[2]-p2[2])*(p1[2]-p2[2]) );
    
    double  b= std::sqrt( (p1[0]-p3[0])*(p1[0]-p3[0]) +
			  (p1[1]-p3[1])*(p1[1]-p3[1]) +
			  (p1[2]-p3[2])*(p1[2]-p3[2]) );
    
    double  c= std::sqrt( (p2[0]-p3[0])*(p2[0]-p3[0]) +
			  (p2[1]-p3[1])*(p2[1]-p3[1]) +
			  (p2[2]-p3[2])*(p2[2]-p3[2]) );
    
    return std::sqrt( ( a+b+c)*(-a+b+c)*( a-b+c)*( a+b-c) )*0.25;
  }

  double pi()
  {
    return 3.14159265;
  }

}
