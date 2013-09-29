#include <math.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct singular_value
{
	int index;
	double value;
} singular_value;

int compare (const void * a, const void * b);
void eyeify(double** A, int n);
void mat_print(double** A, int n);
void vec_print(double* A, int n);
void mat_free(double** A, int n);
void jacobi(double** A, int n, double *s, double **U, double **V);
void rotate(double** A, int n, int p, int q, double** U);
double sign(double val);
double A_transpose_A(double* A, int n, int row, int col);
double max_entry(double** A, int n, int *row, int *col);
singular_value* singular_value_accumulate(double** A, int n);
double** matrix_part_b(int n);
double** matrix_part_c(int n);
double** eye(int n);

int main(int argc, char const *argv[])
{
	if (argc == 1)
	{
		printf("CPSC 445 PSET 1. Invoke as ./ps1 [b|c] [n] \n");
		return 0;
	}
	char part = (char)(*(argv[1]));
	int n = atoi(argv[2]);
	double **A, **U, **V, *s;
	V = eye(n);
	U = eye(n);
	s = malloc(n * sizeof(double));
	if ((part == 'b') || (part == 'B'))
	{
		A = matrix_part_b(n);
	}
	else if ((part == 'C') || (part == 'c'))
	{
		A = matrix_part_c(n);
	}
	else
	{
		printf("Error: invalid part.\n");
		return -1;

	}
	jacobi(A, n, s, U, V);
	printf("\nU = \n");
	mat_print(U, n);
	printf("\nV = \n");
	mat_print(V, n);
	printf("\ns = \n");
	vec_print(s, n);
	printf("\n");

	free(s);
	mat_free(A, n);
	mat_free(U, n);
	mat_free(V, n);
	return 0;
}
//----------------------------------------------------------------------------
double A_transpose_A(double* A, int n, int row, int col)
{
	int i;
	double sum = 0;
	for(i = 0; i < n; ++i)
	{
		sum += A[i][row] * A[i][col];
	}
	return sum;
}
//----------------------------------------------------------------------------
void rotate(double **A, int n, int p, int q, double** U)
{
	int i;
	double a_pp, a_pq, 
	       a_qq, t, 
	       c, s, A_temp;

	a_pp = A_transpose_A((double(*)[n])A, n, p, p);
	a_pq = A_transpose_A((double(*)[n])A, n, p, q);
	a_qq = A_transpose_A((double(*)[n])A, n, q, q);

	t = (a_pp - a_qq) / (2 *(a_pq));
	t = sign(t) / (fabs(t) + sqrt(1 + t * t));
	c = 1 / sqrt(1 + t * t);
	s = (c * t);
	for (i = 0; i < n; ++i)
	{
		A_temp = A[i][p];
		A[i][p] = s * A[i][q] + c * A_temp;
		A[i][q] = c * A[i][q] - s * A_temp;

		A_temp = U[i][p];
		U[i][p] = s * U[i][q] + c * A_temp;
		U[i][q] = c * U[i][q] - s * A_temp;
	}
}
//----------------------------------------------------------------------------
void jacobi(double** A, int n, double *s, double **U, double **V)
{
	int i, j;
	double epsilon = 1e-10, temp;
	int iter_max = 1000000, iters = 0;

	double** T;
	T = eye(n);

	while ((max_entry(A, n, &i, &j) > epsilon * (sqrt(A_transpose_A(A, n, i, i) * A_transpose_A(A, n, j, j)))) && (iters < iter_max))
	{
		rotate(A, n, i, j, T);
		++iters;
	}

	singular_value *sigma;
	sigma = singular_value_accumulate(A, n);
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < n; ++j)
		{
			U[i][j] = A[i][sigma[j].index];
			V[i][j] = T[i][sigma[j].index];
		}
		s[i] = sigma[i].value;
		printf("%f  ", sigma[i].value);
	}
	mat_free(T, n);
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
double** eye(int n)
{
	double** I = malloc(n * sizeof(double));
	int i, j;
	for (i = 0; i < n; ++i)
	{
		I[i] = malloc(n * sizeof(double));
	}
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < n; ++j)
		{
			if (i != j)
			{
				I[i][j] = 0.0;
			}
			else
			{
				I[i][j] = 1.0;
			}
		}
	}
	return I;
}
//----------------------------------------------------------------------------
void eyeify(double** A, int n)
{
	int i, j;
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < n; ++j)
		{
			if (i != j)
			{
				A[i][j] = 0.0;
			}
			else
			{
				A[i][j] = 1.0;
			}
		}
	}
}
//----------------------------------------------------------------------------
void mat_print(double** A, int n)
{
	int i, j;
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < n; ++j)
		{
			printf("%f  ", A[i][j]);
		}
		printf("\n");
	}
}
//----------------------------------------------------------------------------
void vec_print(double* A, int n)
{
	int j;
	for (j = 0; j < n; ++j)
	{
		printf("%e  ", A[j]);
	}
}

//----------------------------------------------------------------------------
void mat_free(double** A, int n)
{
	int i;
	for (i = 0; i < n; ++i)
	{
		free(A[i]);
	}
	free(A);
}
//----------------------------------------------------------------------------
double** matrix_part_b(int n)
{
	double** A = malloc(n * sizeof(double));
	int i, j;
	for (i = 0; i < n; ++i)
	{
		A[i] = malloc(n * sizeof(double));
	}
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < n; ++j)
		{
			A[i][j] = sqrt(i * i + j * j);
		}
	}
	return A;
}
//----------------------------------------------------------------------------
double** matrix_part_c(int n)
{
	double** A = malloc(n * sizeof(double));
	int i, j;
	for (i = 0; i < n; ++i)
	{
		A[i] = malloc(n * sizeof(double));
	}
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < n; ++j)
		{
			A[i][j] = i * i + j * j;
		}
	}
	return A;
}
//----------------------------------------------------------------------------
double max_entry(double** A, int n, int *row, int *col)
{
	int j, i;
	double val, max_val = 0;
	for (i = 0; i < (n - 1); ++i)
	{
		for (j = i + 1; j < n; ++j)
		{
			val = fabs(A_transpose_A(A, n, i, j));
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
int compare (const void * a, const void * b)
{
	if ( (*(singular_value*)a).value >  (*(singular_value*)b).value )
	{
		return -1;
	}
	if ( (*(singular_value*)a).value == (*(singular_value*)b).value )
	{
		return 0;
	}
	if ( (*(singular_value*)a).value <  (*(singular_value*)b).value )
	{
		return 1;
	}
}
//----------------------------------------------------------------------------
singular_value* singular_value_accumulate(double** A, int n)
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
			sum += (A[i][j] * A[i][j]);
		}
		sing_val = sqrt(sum);
		data[j].value = sing_val;
		data[j].index = j;
		for (i = 0; i < n; ++i)
		{
			A[i][j] /= sing_val;
		}
	}
	qsort(data, n, sizeof(singular_value), compare);
	return data;
}



