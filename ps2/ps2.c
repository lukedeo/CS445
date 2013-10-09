//-----------------------------------------------------------------------------
//  Yale CPSC 445a, Problem Set 2
//  Implementation of Steepest Descent
//  By: Luke de Oliveira
//  Yale College `14
//  10/9/13
//-----------------------------------------------------------------------------

#include <math.h>
#include <stdio.h>
#include <stdlib.h>


void dumb_solve(double* a, double* y, int n, double eps, int numit, double* x, int* niter, double* discreps);
