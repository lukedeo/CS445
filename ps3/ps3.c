#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <time.h>

//----------------------------------------------------------------------------
typedef struct pair
{
	int idx;
	double value;
} pair;

//----------------------------------------------------------------------------
int compare (const void * a, const void * b)
{
	if ((*(pair*)a).value > (*(pair*)b).value)
	{
		return 1;
	}
	if ((*(pair*)a).value == (*(pair*)b).value)
	{
		return 0;
	}
	if ((*(pair*)a).value < (*(pair*)b).value)
	{
		return -1;
	}
}

//----------------------------------------------------------------------------
double dist(double x1, double y1, double x2, double y2)
{
	double d1 = x1 - x2;
	double d2 = y1 - y2;
	return sqrt(d1 * d1 + d2 * d2);
}

//----------------------------------------------------------------------------
void dist_matrix(double *D, double *a, int n)
{
	int i, j;
	for (i = 0; i < (n - 1); ++i)
	{	
		D[i * n + i] = 0;
		for (j = (i + 1); j < n; ++j)
		{
			D[i * n + j] = dist(a[2 * i], a[2 * i + 1], a[2 * j], a[2 * j + 1]);
			D[j * n + i] = D[i * n + j];
		}
	}
}

//----------------------------------------------------------------------------
void naive_knn(double *D, int i, int n, int k, int *iz)
{
	pair* pairs;
	pairs = malloc((n - 1) * sizeof(pair));
	int ix = 0, j;
	for (j = 0; j < n; ++j)
	{
		if (j == i)
		{
			continue;
		}
		pairs[ix].idx = j;
		pairs[ix].value = D[i * n + j];
		++ix;
	}
	qsort(pairs, n - 1, sizeof(pair), compare);
	for (j = 0; j < k; ++j)
	{
		iz[i * k + j] = pairs[j].idx;
	}
	free(pairs);
}

//----------------------------------------------------------------------------
void seek_naive(double *a, int n, int k, int *iz)
{
	int i;
	double *D;
	D = malloc(n * n * sizeof(double));
	dist_matrix(D, a, n);
	for (i = 0; i < n; ++i)
	{
		naive_knn(D, i, n, k, iz);
	}
	free(D);
}

//----------------------------------------------------------------------------
void print_matrix(double *a, int row, int col) 
{
	int i, j;
	for (i = 0; i < row; ++i) 
	{
		for (j = 0; j < col; j++)
		{
			printf("%10.4f", a[i*col+j]);
		}
		printf("\n");
	}
	return;
}
void print_int_matrix(int *a, int row, int col) 
{
	int i, j;
	for (i = 0; i < row; ++i) 
	{
		for (j = 0; j < col; ++j)
		{
			printf("%d\t", a[i*col+j]);		
		}
		printf("\n");
	}
	return;
}

int main(int argc, char const *argv[])
{
	int n = 5, i, j, k = 3, *iz;
	double *a, *D;
	a = malloc(2 * n * sizeof(double));
	D = malloc(n * n * sizeof(double));
	iz = malloc(n * k * sizeof(int));


	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < 2; ++j)
		{
			a[i * 2 + j] = i / (j + 1) *( j + 1) + 1;
		}
	}
	dist_matrix(D, a, n);
	seek_naive(a, n, k, iz);
	printf("a matrix\n");
	printDoubleMatrix(a, n, 2);
	printf("D matrix\n");
	printDoubleMatrix(D, n, n);
	printf("iz matrix\n");
	printMatrix(iz, n, k);







	return 0;
}