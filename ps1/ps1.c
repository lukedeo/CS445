//-----------------------------------------------------------------------------
//  Yale CPSC 445a, Problem Set 1
//  Implementation of SVD (One-Sided Jacobi Rotations)
//  By: Luke de Oliveira
//  Yale College `14
//  9/30/13
//-----------------------------------------------------------------------------

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
//-----------------------------------------------------------------------------
//	Stuff for sorting the columns in decreasing order of singular values
//-----------------------------------------------------------------------------
typedef struct singular_value
{
	int index;
	double value;
} singular_value;
int compare (const void * a, const void * b);
//-----------------------------------------------------------------------------
//	Pre declarations of some functions, and the main Jacobi one.
//-----------------------------------------------------------------------------
void jacobi(double* a, int n, double *s, double *u, double *v);
void eyeify(double* a, int n);
void mat_print(double* a, int n);
void vec_print(double* a, int n);
void rotate(double* a, int n, int p, int q, double* u);
double sign(double val);
double gram_entry(double* a, int n, int row, int col);
double max_entry(double* a, int n, int *row, int *col);
singular_value* singular_value_accumulate(double* a, int n);
double* matrix_part_b(int n);
double* matrix_part_c(int n); 
double* eye(int n);
double off_norm(double* a, int n);
//-----------------------------------------------------------------------------
//	Main function, to be deleted.
//-----------------------------------------------------------------------------
int main(int argc, char const *argv[])
{
	if (argc == 1)
	{
		printf("CPSC 445 PSET 1. Invoke as ./ps1 [b|c] [n] \n");
		return 0;
	}
	char part = (char)(*(argv[1]));
	int n = atoi(argv[2]);
	double *a, *u, *v, *s;
	v = eye(n);
	u = eye(n);
	s = malloc(n * sizeof(double));
	if ((part == 'b') || (part == 'B'))
	{
		a = matrix_part_b(n);
	}
	else if ((part == 'C') || (part == 'c'))
	{
		a = matrix_part_c(n);
	}
	else
	{
		printf("Error: invalid part.\n");
		return -1;
	}

	jacobi(a, n, s, u, v);

	a = matrix_part_b(n);

	printf("\na = \n");
	mat_print(a, n);
	printf("\nu = \n");
	mat_print(u, n);
	printf("\nv = \n");
	mat_print(v, n);
	printf("\ns = \n");
	vec_print(s, n);
	printf("\n");

	free(s);
	free(a);
	free(u);
	free(v);
	return 0;
}
//-----------------------------------------------------------------------------
//	Implementations of functions.
//-----------------------------------------------------------------------------
double gram_entry(double* a, int n, int row, int col)
{
	int i;
	double sum = 0;
	for(i = 0; i < n; ++i)
	{
		sum += a[i * n + row] * a[i * n + col];
	}
	return sum;
}
//----------------------------------------------------------------------------
void rotate(double *a, int n, int p, int q, double* u)
{
	int i;
	double a_pp, a_pq, 
           a_qq, t, 
           c, s, a_temp;
	a_pp = gram_entry(a, n, p, p);
	a_pq = gram_entry(a, n, p, q);
	a_qq = gram_entry(a, n, q, q);
	t = (a_pp - a_qq) / (2 *(a_pq));
	t = sign(t) / (fabs(t) + sqrt(1 + t * t));
	c = 1 / sqrt(1 + t * t);
	s = (c * t);
	for (i = 0; i < n; ++i)
	{
		a_temp = a[i * n + p];
		a[i * n + p] = s * a[i * n + q] + c * a_temp;
		a[i * n + q] = c * a[i * n + q] - s * a_temp;

		a_temp = u[i * n + p];
		u[i * n + p] = s * u[i * n + q] + c * a_temp;
		u[i * n + q] = c * u[i * n + q] - s * a_temp;
	}
}
//----------------------------------------------------------------------------
void jacobi(double* a, int n, double *s, double *u, double *v)
{
	int i, j, nullspace = 0;
	double epsilon = 1e-15, comparison;
	int iter_max = 10000, iters = 0, count = 1;

	double* T;
	T = eye(n);

	while ((count != 0) && (iters < iter_max))
	{
		count = 0;
		for (i = 0; i < (n - 1); ++i)
		{
			for (j = i + 1; j < n; ++j)
			{
				comparison = epsilon * (sqrt(gram_entry(a, n, i, i) * gram_entry(a, n, j, j)));
				if (fabs(gram_entry(a, n, i, j)) > comparison)
				{
					rotate(a, n, i, j, T);
					++count;
				}
			}
		}
		++iters;
	}
	if (iters == iter_max)
	{
		printf("\nJacobi's Method failed to converge. Terminating\n");
		return;
	}
	else
	{
		printf("\nJacobi's method converged in %i iterations.\n", iters);
	}
	singular_value *sigma;
	sigma = singular_value_accumulate(a, n);
	printf("\nBelow are the singular values > 1.0E-13:\n\n");
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < n; ++j)
		{
			u[i * n + j] = a[i * n + sigma[j].index];
			v[i * n + j] = T[i * n + sigma[j].index];
		}
		s[i] = sigma[i].value;
		if (sigma[i].value > 1.0e-13)
		{
			printf("%e  ", sigma[i].value);
		}
		else
		{
			++nullspace;
		}
	}
	if (nullspace > 0)
	{
		printf("\n\nThe matrix is probably not of full rank. n = %i, and rank(A) is approximately %i.\n", n, n - nullspace);
	}
	else
	{
		printf("\n\nThe matrix is probably of full rank.\n");
	}
	printf("\nBelow is the column vector associated with the largest singular value (eigenvalue of the implicit Gram Matrix):\n");
	for (i = 0; i < n; ++i)
	{
		printf("%e\n", u[i * n]);
	}

	free(T);
	free(sigma);
}
//----------------------------------------------------------------------------
double sign(double val)
{
	if(val < 0)
	{
		return -1;
	}
	return 1;
}
//----------------------------------------------------------------------------
double* eye(int n)
{
	double* I = malloc(n * n * sizeof(double));
	int i, j;
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < n; ++j)
		{
			if (i != j)
			{
				I[i * n + j] = 0.0;
			}
			else
			{
				I[i * n + j] = 1.0;
			}
		}
	}
	return I;
}
//----------------------------------------------------------------------------
void eyeify(double* a, int n)
{
	int i, j;
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < n; ++j)
		{
			if (i != j)
			{
				a[i * n + j] = 0.0;
			}
			else
			{
				a[i * n + j] = 1.0;
			}
		}
	}
}
//----------------------------------------------------------------------------
void mat_print(double* a, int n)
{
	int i, j;
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < n; ++j)
		{
			printf("%f  ", a[i * n + j]);
		}
		printf("\n");
	}
}
//----------------------------------------------------------------------------
void vec_print(double* a, int n)
{
	int j;
	for (j = 0; j < n; ++j)
	{
		printf("%e  ", a[j]);
	}
}
//----------------------------------------------------------------------------
double* matrix_part_b(int n)
{
	double* a = malloc(n * n * sizeof(double));
	int i, j;
	double i_p_1, j_p_1;
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < n; ++j)
		{
			i_p_1 = ((double) i) + 1;
			j_p_1 = ((double) j) + 1;
			a[i * n + j] = sqrt(i_p_1 * i_p_1 + j_p_1 * j_p_1);
		}
	}
	return a;
}
//----------------------------------------------------------------------------
double* matrix_part_c(int n)
{
	double* a = malloc(n * n * sizeof(double));
	int i, j;
	double i_p_1, j_p_1;
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < n; ++j)
		{
			i_p_1 = ((double) i) + 1;
			j_p_1 = ((double) j) + 1;
			a[i * n + j] = i_p_1 * i_p_1 + j_p_1 * j_p_1;
		}
	}
	return a;
}
//----------------------------------------------------------------------------
double max_entry(double* a, int n, int *row, int *col)
{
	int j, i;
	double val, max_val = 0;
	for (i = 0; i < (n - 1); ++i)
	{
		for (j = i + 1; j < n; ++j)
		{
			val = fabs(gram_entry(a, n, i, j));
			if (val > max_val)
			{
				max_val = val;
				(*row) = i;
				(*col) = j;
			}
		}
	}
	return max_val;
}
//----------------------------------------------------------------------------
double off_norm(double* a, int n)
{
	int j, i;
	double val, sum = 0;

	for (i = 0; i < (n - 1); ++i)
	{
		for (j = i + 1; j < n; ++j)
		{
			val = fabs(gram_entry(a, n, i, j));
			sum += (val * val);
		}
	}
	return sqrt(2 * sum);
}
//----------------------------------------------------------------------------
int compare (const void * a, const void * b)
{
	if ((*(singular_value*)a).value > (*(singular_value*)b).value)
	{
		return -1;
	}
	if ((*(singular_value*)a).value == (*(singular_value*)b).value)
	{
		return 0;
	}
	if ((*(singular_value*)a).value <  (*(singular_value*)b).value)
	{
		return 1;
	}
}
//----------------------------------------------------------------------------
singular_value* singular_value_accumulate(double* a, int n)
{
	singular_value* data;
	data = malloc(n * sizeof(singular_value));
	int i, j;
	double sum, sing_val;
	for (j = 0; j < n; ++j)
	{
		sum = 0;
		for (i = 0; i < n; ++i)
		{
			sum += (a[i * n + j] * a[i * n + j]);
		}
		sing_val = sqrt(sum);
		data[j].value = sing_val;
		data[j].index = j;
		for (i = 0; i < n; ++i)
		{
			a[i * n + j] /= sing_val;
		}
	}
	qsort(data, n, sizeof(singular_value), compare);
	return data;
}



