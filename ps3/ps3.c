#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <time.h>

const int BUF = 100000;


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
// SEEK NAIVE
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
void int_swap(int *a, int *b) 
{
	int tmp = *a;
	*a = *b;
	*b = tmp;
	return;
}
//----------------------------------------------------------------------------

typedef struct point
{
	double x;
	double y;
} point;

typedef struct box
{
	point b_left;
	point t_right;
} box;

typedef struct circle
{
	point center;
	double radius;
} circle;

typedef struct control_object
{
	int is_processed;
	int start;
	int end;
	int parent;
	int children[4];
	box Box;
} control_object;

int is_in_box(point p, box B)
{
	if ((p.x > B.Box.t_right.x) || (p.x < B.Box.b_left.x))
	{
		return 0;
	}
	if ((p.y > B.Box.t_right.y) || (p.y < B.Box.b_left.y))
	{
		return 0;
	}
	return 1;
}

void arr_to_points(double *a, point *points, int n)
{
	for (int i = 0; i < n; ++i)
	{
		points[i].x = a[2 * i];
		points[i].y = a[2 * i + 1];
	}
}


void print_point(point p)
{
	printf("(%f, %f)\n", p.x, p.y);
}


void gen_sub_boxes(box X, box *a, box *b, box *c, box *d)
{
	a->b_left.x = X.b_left.x;
	a->b_left.y = (X.t_right.y - X.b_left.y) / 2;

	a->t_right.y = X.t_right.y;
	a->t_right.x = (X.t_right.x - X.b_left.x) / 2;

	b->b_left.x = (X.t_right.x - X.b_left.x) / 2;
	b->b_left.y = (X.t_right.y - X.b_left.y) / 2;

	b->t_right.y = X.t_right.y;
	b->t_right.x = X.t_right.x;


	c->b_left.x = X.b_left.x;
	c->b_left.y = X.b_left.y;

	c->t_right.y = b->b_left.y;
	c->t_right.x = b->b_left.x;

	d->b_left.x = (X.t_right.x - X.b_left.x) / 2;
	d->b_left.y = X.b_left.y;

	d->t_right.y = (X.t_right.y - X.b_left.y) / 2;
	d->t_right.x = X.t_right.x;
}


//----------------------------------------------------------------------------


void add_control_entry(control_object *control, int *permutation, int current, int next, point *Data, int n)
{
	int head = control[current].start;
	int tail = control[current].end;

	int i;
	for (i = 0; i < 4; ++i)
	{
		control[current].children[i] = next + i;
		control[next + i].parent = current;
	}
	gen_sub_boxes(control[current].Box, &control[next].Box, 
		                                &control[next + 1].Box, 
		                                &control[next + 2].Box, 
		                                &control[next + 3].Box);
	printf("%s\n", "parent box:");
	print_point(control[current].Box.b_left);
	print_point(control[current].Box.t_right);

	int count = 0;

	while(count < 4)
	{
		int s = head;
		int e = tail;
		while(s <= e)
		{
			if (is_in_box(Data[permutation[s]], control[next+count].Box))
			{
				s++;
			}
			else
			{
				int_swap(&permutation[s], &permutation[e]);
			}
		}
		control[next+count].start = head;
		control[next+count].end = s - 1;
		count++;
		head = s;
	}	
}

//----------------------------------------------------------------------------
void seek(double *a, int n, int k, int *iz) 
{
	point *Data;
	Data = malloc(n * sizeof(point));
	arr_to_points(a, Data, n);


	int i, j;
	int *permutation = (int*) malloc(n * sizeof(int));
	for (i = 0; i < n; i++)
	{
		permutation[i] = i;
	}

	control_object *Control;
	Control = malloc(BUF * sizeof (control_object));
	int current = 0;
	Control[0].is_processed = 3, Control[0].start = 0, Control[0].end = n - 1;
	Control[0].parent = -8, 
	Control[0].Box.b_left.x = 0, Control[0].Box.b_left.y = 0;
	Control[0].Box.t_right.x = 1, Control[0].Box.t_right.y = 1;
	print_point(Control[current].Box.b_left);
	print_point(Control[current].Box.t_right);

	
	int next = 1;

	while(current < next)
	{
		int e = Control[current].end;
		int s = Control[current].start;
		if (e - s > k)
		{
			add_control_entry(Control, permutation, current, next, Data, n);
			next += 4;
		}
		Control[current].is_processed;
		current++;
	}


	for (i = 0; i < 4; ++i)
	{
		printf("\nChild box #%d:\n", i+1);
		print_point(Control[next+i].Box.b_left);
		print_point(Control[next+i].Box.t_right);
	}

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
	seek(a, n, k, iz);
	printf("a matrix\n");
	print_matrix(a, n, 2);
	printf("D matrix\n");
	print_matrix(D, n, n);
	printf("iz matrix\n");
	print_int_matrix(iz, n, k);







	return 0;
}