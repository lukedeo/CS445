#include <math.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct singular_value
{
	int index;
	double value;
} singular_value;

int compare (const void * a, const void * b);
void eyeify(double* a, int n);
void mat_print(double* a, int n);
void vec_print(double* a, int n);
// void mat_free(double* a, int n);
void jacobi(double* a, int n, double *s, double *u, double *v);
void rotate(double* a, int n, int p, int q, double* u);
double sign(double val);
double a_transpose_a(double* a, int n, int row, int col);
double max_entry(double* a, int n, int *row, int *col);
singular_value* singular_value_accumulate(double* a, int n);
double* matrix_part_b(int n);
double* matrix_part_c(int n);
double* eye(int n);

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
//----------------------------------------------------------------------------
double a_transpose_a(double* a, int n, int row, int col)
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

	a_pp = a_transpose_a(a, n, p, p);
	a_pq = a_transpose_a(a, n, p, q);
	a_qq = a_transpose_a(a, n, q, q);

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
	int i, j;
	double epsilon = 1e-10, temp;
	int iter_max = 1000000, iters = 0;

	double* T;
	T = eye(n);

	while ((max_entry(a, n, &i, &j) > epsilon * (sqrt(a_transpose_a(a, n, i, i) * a_transpose_a(a, n, j, j)))) && (iters < iter_max))
	{
		rotate(a, n, i, j, T);
		++iters;
	}

	singular_value *sigma;
	sigma = singular_value_accumulate(a, n);
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < n; ++j)
		{
			u[i * n + j] = a[i * n + sigma[j].index];
			v[i * n + j] = T[i * n + sigma[j].index];
		}
		s[i] = sigma[i].value;
		printf("%f  ", sigma[i].value);
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
// void mat_free(double* a, int n)
// {
// 	int i;
// 	for (i = 0; i < n; ++i)
// 	{
// 		free(a[i]);
// 	}
// 	free(a);
// }
//----------------------------------------------------------------------------
double* matrix_part_b(int n)
{
	double* a = malloc(n * n * sizeof(double));
	int i, j;
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < n; ++j)
		{
			a[i * n + j] = sqrt(i * i + j * j);
		}
	}
	return a;
}
//----------------------------------------------------------------------------
double* matrix_part_c(int n)
{
	double* a = malloc(n * n * sizeof(double));
	int i, j;
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < n; ++j)
		{
			a[i * n + j] = i * i + j * j;
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
			val = fabs(a_transpose_a(a, n, i, j));
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



