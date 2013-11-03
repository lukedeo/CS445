//-----------------------------------------------------------------------------
//  Yale CPSC 445a, Problem Set 2
//  Implementation of Steepest Descent on Ax = y
//  By: Luke de Oliveira
//  Yale College `14
//  10/9/13
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
//	Standard include files
//-----------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


//-----------------------------------------------------------------------------
//	Forward-Declarations of functions, including dumb_solve
//-----------------------------------------------------------------------------

void dumb_solve(double *a, double *y, int n, double eps, int numit, double *x, int *niter, double *discreps);
double line_search(double *a, double *Ax_b, int n);
double residual_norm(double *Ax_b, int n);
void general_multiply(double *left, double *right, double *result, int m, int n);
void matrix_vector(double *matrix, double *vector, double *result, int n);
double *gradient(double *a, double *Ax_b, int n);
void print(double *matrix, int m, int n);
void vector_subtract(double *left, double *right, double *result, int n);
double norm(double *vector, int n);
void matrix_transpose(double *matrix, int n);
void pair_gen(double *A, double *y, int n);

//-----------------------------------------------------------------------------
//	main routine which runs through the tests.
//-----------------------------------------------------------------------------
int main(int argc, char const *argv[])
{
	int n = atoi(argv[1]);
	int numit = 10000, niter, i;
	double eps = 10e-6;
	double *matrix, *vector, *discreps, *x, *temp, *x_real;
	matrix = (double *) malloc (sizeof (double) * (n * n));
	vector = (double *) malloc (sizeof (double) * n);
	discreps = (double *) malloc (sizeof (double) * numit);
	x = (double *) malloc (sizeof (double) * n);
	x_real = (double *) malloc (sizeof (double) * n);
	temp = (double *) malloc (sizeof (double) * n);

	pair_gen(matrix, vector, n);

	for (i = 0; i < n; ++i)
	{
		x_real[i] = vector[i] / matrix[i * n + i];
	}

	//solve the system!
	dumb_solve(matrix, vector, n, eps, numit, x, &niter, discreps);

	//Print the results
	printf("Original A:");
	print(matrix, n, n);

	printf("Solution to linear system");
	print(x, n, 1);

	printf("Analytical Solution to linear system");
	print(x_real, n, 1);

	printf("Multiply solution by A to see if we get the original vector");
	general_multiply(matrix, x, temp, n, 1);        
	print(temp, n, 1);
	printf("Number of iterations: %d\n", niter);
	printf("Discrepancy values:\n");
	for (i = 0; i < niter; i++) 
	{
		// if (i == (numit - 1)) 
		// {
		// 	break;
		// }

		printf("%i, %.8f\n", i, discreps[i]);
	}
	printf("\n");
	free(matrix); 
	free(vector); 
	free(discreps); 
	free(x_real); 
	free(x); 
	free(temp);
	return 0;
}

//-----------------------------------------------------------------------------
//	pair_gen: for any n, generates the A and y pair as specified in ps2
//-----------------------------------------------------------------------------
void pair_gen(double *A, double *y, int n)
{
	int i;
	for (i = 0; i < n; ++i)
	{
		double ip1 = (double)(i + 1);
		A[i * n + i] = 1 / ((ip1) * (ip1));
		y[i] = 1.0;
	}
}
//-----------------------------------------------------------------------------
//	dumb_solve: implementation of steepest descent on the Ax = y system
//-----------------------------------------------------------------------------
void dumb_solve(double *a, double *y, int n, double eps, int numit, double *x, int *niter, double *discreps) 
{
	int i;
	double *grad = NULL; // init ptr to carry the gradient of the objective function.
	double step_size = 0.0, error = 0.0;
	double *Ax = (double *) malloc (sizeof (double) * n);
	double *Ax_b = (double *) malloc (sizeof (double) * n);
	double *x_n = (double *) malloc (sizeof (double) * n);

	//We let the initial guess be the vector of ones
	for (i = 0; i < n; i++) 
	{
		x_n[i] = 1;
	}

	for (*niter = 0; (*niter < numit); (*niter)++) 
	{
		matrix_vector(a, x_n, Ax, n); // Calculates Ax
		vector_subtract(Ax, y, Ax_b, n); // Calculates Ax - b, sticks into Ax_b
		error = residual_norm(Ax_b, n);
		if (error < eps) // terminate if our error is small enough
		{
			break;
		}
		discreps[*niter] = error; // dump the error here

		//Calculare the new solution iterate.
		grad = gradient(a, Ax_b, n);
		step_size = line_search(a, grad, n);
		for (i = 0; i < n; i++)
		{
			grad[i] = step_size * grad[i];
		}
		vector_subtract(x_n, grad, x_n, n);
		free(grad);
	}
	/* Copy over answer to solution vector x */
	for (i = 0; i < n; i++)
	{
		x[i] = x_n[i];
	}
	free(x_n);
	free(Ax); 
	free(Ax_b);
	return;
}
//-----------------------------------------------------------------------------
//	line_search: given A and grad(f(x)), calculates optimal step size.
//-----------------------------------------------------------------------------
double line_search (double *a, double *grad, int n) {
	int i;
	double *denom = (double *) malloc (sizeof (double) * n);
	double numerator = 0.0;
	double denominator = 0.0;

	/* Calculate norm squared of gradient */
	for (i = 0; i < n; i++) 
	{
		numerator += grad[i] * grad[i];
	}

	matrix_vector(a, grad, denom, n);

	/* Calculate the norm squared of A*gradient, i.e. denominator */
	for (i = 0; i < n; i++) 
	{
		denominator += denom[i] * denom[i];
	}

	free(denom);
	return (numerator / (2.0 * denominator));
}
//----------------------------------------------------------------------------
double residual_norm(double *Ax_b, int n) 
{
	int i;
	double out = 0.0;
	// Dot (Ax - b) with itself which equals ||Ax - b|| ^ 2
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
//----------------------------------------------------------------------------
void matrix_vector(double *matrix, double *vector, double *result, int n)
{
	int i, j;
	for (i = 0; i < n; ++i)
	{
		double temp = 0.0;
		for (j = 0; j < n; ++j)
		{
			temp += matrix[i * n + j] * vector[j];
		}
		result[i] = temp;
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
	matrix_transpose(at, n);

	// Calculate 2At(Ax - b) 
	matrix_vector(at, Ax_b, gradient, n);
	free (at);
	return gradient;
}
//----------------------------------------------------------------------------
void print(double *matrix, int m, int n) 
{
	int i, j;
	printf("\n");
	for (i = 0; i < m; i++) 
	{
		for (j = 0; j < n; j++)
		{
			printf("%.9f\t", matrix[i * n + j]);
		}
		printf("\n");
	}
	printf("\n");
	return;
}
//----------------------------------------------------------------------------
void vector_subtract(double *left, double *right, double *result, int n) 
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
	return sqrt(result);
}
//----------------------------------------------------------------------------
void matrix_transpose (double *matrix, int n) 
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
