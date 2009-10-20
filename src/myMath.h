/**
 * Filename     : myMath.h
 * Description  : A set of functions dealing with various mathematical operations.
 * Author(s)    : Patrik Sahlin (sahlin@thep.lu.se)
 * Created      : August 2008
 * Revision     : $Id$
 */
#ifndef MYMATH_H
#define MYMATH_H

#include <vector>

// This function takes an n x n real symmetric matrix and finds the
// eigenvalues and eigenvectors. The algorithm is based on Numerical
// Recipes, p. 463-466. The implementation is loosely based on
// p. 467-469. It might need a better check for convergence.

size_t jacobiTransformation(std::vector< std::vector<double> > A, std::vector< std::vector<double > > &V, std::vector<double> &d);

#endif /* MYMATH_H */
