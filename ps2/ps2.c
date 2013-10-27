//-----------------------------------------------------------------------------
//  Yale CPSC 445a, Problem Set 2
//  Implementation of Steepest Descent
//  By: Luke de Oliveira
//  Yale College `14
//  10/9/13
//-----------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


void dumb_solve(double *a, double *y, int n, double eps, int numit, double *x, int *niter, double *discreps) 
{
	int i;
	double *grad = NULL; // init ptr to carry the gradient of the objective
	double step_size = 0.0, error = 0.0;
	double *Ax = (double *) malloc (sizeof (double) * n);
	double *Ax_b = (double *) malloc (sizeof (double) * n);
	double *x_n = (double *) malloc (sizeof (double) * n);


	//We let the initial guess be zero
	for (i = 0; i < n; i++) 
	{
		x_n[i] = 0;
	}

	for (*niter = 0; (*niter < numit); (*niter)++) 
	{
		general_multiply (a, x_n, Ax, n, 1); // Calculates Ax
		subtract (Ax, y, Ax_b, n); // Calculates Ax - b, sticks into Ax_b

		error = residual_norm(Ax_b, n);

		if (error < eps) // terminate if our error is small enough
		{
			break;
		}
		discreps[*niter] = error; // dump the 

		//Calculare the new solution iterate.

		grad = gradient(a, Ax_b, n);
		step_size = stepSize (a, Ax_b, n);

		for (i = 0; i < n; i++)
		{
			grad[i] = step_size * grad[i];
		}
		subtract (x_n, grad, x_n, n);
		free (grad);
	}

	/* Copy over answer to solution vector x */
	for (i = 0; i < n; i++)
	{
		x[i] = x_n[i];
	}

	free (x_n); 
	free (Ax); 
	free (Ax_b);
	return;
}
//----------------------------------------------------------------------------
double stepSize (double *a, double *Ax_b, int n) {
	int i;
	double *at = (double *) malloc (sizeof (double) * (n * n));
	double *aat = (double *) malloc (sizeof (double) * (n * n));
	double *nom = (double *) malloc (sizeof (double) * n);
	double *denom = (double *) malloc (sizeof (double) * n);
	double nominator = 0.0;
	double denominator = 0.0;

	/* Calculate A transpose */
	for (i = 0; i < (n * n); i++) 
	{
		at[i] = a[i];
	}
	transpose (at, n);

	/* Calculate At * (Ax - b) */
	general_multiply (at, Ax_b, nom, n, 1);

	/* Calculate norm squared of At * (Ax - b), i.e. nominator */
	for (i = 0; i < n; i++) 
	{
		nominator += nom[i] * nom[i];
	}

	/* Calculate A*At */
	general_multiply (a, at, aat, n, n);

	/* Calculate A*At * (Ax - b) */
	general_multiply (aat, Ax_b, denom, n, 1);

	/* Calculate the norm squared of A*At * (Ax - b), i.e. denominator */
	for (i = 0; i < n; i++) 
	{
		denominator += denom[i] * denom[i];
	}

	free (at); free (aat); free (nom); free (denom);
	return (nominator / (2.0 * denominator));
}
//----------------------------------------------------------------------------
double residual_norm(double *Ax_b, int n) 
{
	int i;
	double out = 0.0;
	/* Dot (Ax - b) with itself which is equal to || Ax - b || ^ 2 */
	for (i = 0; i < n; i++) 
	{
		out += Ax_b[i] * Ax_b[i];
	}
	return out;
}
//----------------------------------------------------------------------------
void general_multiply (double *left, double *right, double *result, int m, int n) 
{
	int i, k, j;
	double temp = 0;

	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n; j++)
		{
			for (k = 0; k < m; k++)
			{
				temp += left[i * m + k] * right[k * n + j];
			}	
			result[i * n + j] = temp;
			temp = 0;
		}
	}
	return;
}

void column_gram_matrix(double *A, double *result, int n)
{
	int i, j, k;
	for (i = 0; i < n; ++i)
	{
		for (j = i; j < n; ++j)
		{
			double temp = 0.0;
			for (k = 0; k < n; ++k)
			{
				temp += A[i * n + k] * A[k * n + k];
			}
			result[i * n + j] = temp;
			if (i != j)
			{
				result[j * n + i] = temp;
			}
		}
	}
	return;
}



//----------------------------------------------------------------------------
double *gradient (double *a, double *Ax_b, int n) 
{
	int i;
	double *at = (double *) malloc (sizeof (double) * (n * n));
	double *gradient = (double *) malloc (sizeof (double) * n);

	/* Calculate 2 * At */
	for (i = 0; i < (n * n); i++) 
	{
		at[i] = 2 * a[i];
	}
	transpose (at, n);

	/* Calculate 2At(Ax - b) */
	general_multiply (at, Ax_b, gradient, n, 1);

	free (at);
	return gradient;
}
//----------------------------------------------------------------------------
void print (double *matrix, int m, int n) 
{
	int i, j;
	printf ("\n");
	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++)
			printf ("%.9f\t", matrix[i * n + j]);
		printf ("\n");
	}
	printf ("\n");
	return;
}
//----------------------------------------------------------------------------
void subtract (double *left, double *right, double *result, int n) 
{
	int i;
	for (i = 0; i < n; i++)
	{
		result[i] = left[i] - right[i];
	}
	return;
}
//----------------------------------------------------------------------------
double norm (double *vector, int n) 
{
	double result = 0.0;
	int i = 0;
	for (i = 0; i < n; i++)
	{
		result += (vector[i] * vector[i]);
	}		
	return sqrt (result);
}
//----------------------------------------------------------------------------
void transpose (double *matrix, int n) 
{
	double temp;
	int i, j;
	for (i = 0; i < n; i++) 
	{
		for (j = i + 1; j < n; j++) 
		{
			if (j == i) 
			{
				continue;
			}
			temp = matrix[i * n + j];
			matrix[i * n + j] = matrix[j * n + i];
			matrix[j * n + i] = temp;
		}
	}
	return;
}
