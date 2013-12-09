//-----------------------------------------------------------------------------
//  Yale CPSC 445a, Problem Set 3
//  Quad-tree k-NN
//  By: Luke de Oliveira
//  Yale College `14
//  12/30/13
//-----------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define BUFFER 1000000

//-----------------------------------------------------------------------------
//	Utility Structures
//-----------------------------------------------------------------------------
typedef struct pair
{
	int idx;
	double value;
} pair;
//----------------------------------------------------------------------------
typedef struct point
{
	double x;
	double y;
} point;
//----------------------------------------------------------------------------
typedef struct box
{
	point b_left, b_right;
	point t_left, t_right;
} box;
//----------------------------------------------------------------------------
typedef struct circle
{
	point center;
	double radius;
} circle;
//----------------------------------------------------------------------------
typedef struct control_object
{
	int is_processed;
	int start;
	int end;
	int parent;
	int children[4];
	box Box;
} control_object;
//----------------------------------------------------------------------------
typedef struct perf_t
{
	int n;
	int k;
	int passed;
	double fast_time;
	double naive_time;
} perf_t;
//-----------------------------------------------------------------------------
//	Predeclarations
//-----------------------------------------------------------------------------

// Utility functions
void print_int_matrix(int *a, int row, int col);
void arr_to_points(double *a, point *points, int n);
double distance(point p, point q);
double maximum_real(double a, double b);
int compare(const void * a, const void * b);
void int_swap(int *a, int *b);
int is_in_box(point p, box B);
void print_point(point p);
void print_matrix(double *a, int row, int col);
void print_perf_array(perf_t *perf, int n);
void print_perf(perf_t perf);
perf_t test(int n, int k);
double get_radius(point p, box B);
int intersect(box B, circle C);

//Functions for seek_naive
void naive_knn(int idx, point *Data, int n, int k, int *iz, int verbose);
void seek_naive(double *a, int n, int k, int *iz);

//functions for seek
void gen_sub_boxes(box X, box *a, box *b, box *c, box *d, int verbose);
void add_control_entry(control_object *control, int *permutation, int current, int next, point *Data, int n, int verbose);
int find_parent_idx(point p, control_object *control);
void knn(int idx, point *Data, int points[], int length, int n, int k, int *iz, int verbose);
void seek(double *a, int n, int k, int *iz, int verbose);

//-----------------------------------------------------------------------------
//	INSERT MAIN FUNCTION HERE
//-----------------------------------------------------------------------------

int main(int argc, char const *argv[])
{
	// int n = atoi(argv[1]), k = atoi(argv[2]);
	// srand(2342342341);
	// perf_t t = test(n, k);
	// print_perf(t);
	// 
	// 
	perf_t t;
	int n, k;
	for (n = 20; n <= 1000; n += 10)
	{
		for (k = 1; k <= 15; k += 4)
		{
			t = test(n, k);
			print_perf(t);
		}
	}


	return 0;
}
//-----------------------------------------------------------------------------
//	Functions that do useful stuff
//-----------------------------------------------------------------------------
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
//----------------------------------------------------------------------------
void arr_to_points(double *a, point *points, int n)
{
	int i;
	for (i = 0; i < n; ++i)
	{
		points[i].x = a[2 * i];
		points[i].y = a[2 * i + 1];
	}
}
//----------------------------------------------------------------------------
double distance(point p, point q)
{
	double dx = p.x - q.x;
	double dy = p.y - q.y;
	return sqrt((dx * dx) + (dy * dy));
}
//----------------------------------------------------------------------------
double maximum_real(double a, double b)
{
	if (a >= b)
	{
		return a;
	}
	return b;
}
//----------------------------------------------------------------------------
int compare(const void * a, const void * b)
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
void int_swap(int *a, int *b) 
{
	int _tmp = *a;
	*a = *b;
	*b = _tmp;
	return;
}
//----------------------------------------------------------------------------
int is_in_box(point p, box B)
{
	if ((p.x > B.t_right.x) || (p.x < B.b_left.x))
	{
		return 0;
	}
	if ((p.y > B.t_right.y) || (p.y < B.b_left.y))
	{
		return 0;
	}
	return 1;
}
//----------------------------------------------------------------------------
void print_point(point p)
{
	printf("(%f, %f)\n", p.x, p.y);
}
//----------------------------------------------------------------------------
double get_radius(point p, box B) 
{
	double result = distance(p, B.b_left);
	result = maximum_real(result, distance(p, B.b_right));
	result = maximum_real(result, distance(p, B.t_right));
	return maximum_real(result, distance(p, B.t_left));
}
//----------------------------------------------------------------------------
void print_matrix(double *a, int row, int col) 
{
	int i, j;
	for (i = 0; i < row; ++i) 
	{
		for (j = 0; j < col; j++)
		{
			printf("%16.10f", a[i*col+j]);
		}
		printf("\n");
	}
	return;
}
//----------------------------------------------------------------------------
void print_perf_array(perf_t *perf, int n)
{
	int i;
	printf("N, k, Passed, Fast.time, Naive.time\n");
	for (i = 0; i < n; ++i)
	{
		printf("%d, %d, %d, %f, %f\n", perf[i].n, 
			                           perf[i].k, 
			                           perf[i].passed, 
			                           perf[i].fast_time, 
			                           perf[i].naive_time);
	}
}
//----------------------------------------------------------------------------
void print_perf(perf_t perf)
{
	printf("%d, %d, %d, %f, %f\n", perf.n, 
		                           perf.k, 
		                           perf.passed, 
		                           perf.fast_time, 
		                           perf.naive_time);
}
//----------------------------------------------------------------------------
perf_t test(int n, int k)
{
	int i, j, *iz, *iz2, t1, t2, t3, ok = 1;
	
	double *a;
	a = (double*) malloc(2 * n * sizeof(double));
	iz = (int*) malloc(n * k * sizeof(int));
	iz2 = (int*) malloc(n * k * sizeof(int));
	// srand((unsigned)time(NULL));
	
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < 2; ++j)
		{
			a[i * 2 + j] = ((double)rand()/(double)RAND_MAX);
		}
	}
	t1 = clock();
	seek_naive(a, n, k, iz);
	// print_matrix(a, n, 2);
	t2 = clock();
	seek(a, n, k, iz2, 0);
	// printf("\n\n\n\n\n");
	// print_matrix(a, n, 2);
	t3 = clock();
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < k; ++j)
		{
			if (iz[i * k + j] != iz2[i * k + j])
			{
				ok = 0;
				// printf("Diff: (%d, %d) = %d, %d\n", i, j, iz[i * k + j], iz2[i * k + j]);
				// break;
			}
		}
	}
	free(a);
	a = NULL;
	free(iz);
	iz = NULL;
	free(iz2);
	iz2 = NULL;

	perf_t _tmp = {.n = n, .k = k, .passed = ok, 
		.fast_time = (double)(t3-t2)/(double)CLOCKS_PER_SEC, 
		.naive_time = (double)(t2-t1)/(double)CLOCKS_PER_SEC};
	return _tmp;
}
//-----------------------------------------------------------------------------
//	Functions for and including seek_naive
//-----------------------------------------------------------------------------
void naive_knn(int idx, point *Data, int n, int k, int *iz, int verbose)
{
	pair* pairs;
	pairs = (pair*)malloc((n-1) * sizeof(pair));
	int ix = 0, j;
	for (j = 0; j < (n); ++j)
	{
		if (j == idx)
		{
			continue;
		}
		pairs[ix].idx = j;
		pairs[ix].value = distance(Data[idx], Data[j]);
		++ix;
	}
	qsort(pairs, n - 1, sizeof(pair), compare);
	for (j = 0; j < k; ++j)
	{
		iz[idx * k + j] = pairs[j].idx;
	}
	free(pairs);
	pairs = NULL;
}
//----------------------------------------------------------------------------
void seek_naive(double *a, int n, int k, int *iz)
{
	int i;
	point *Data;
	Data = (point*) malloc(n * sizeof(point));
	arr_to_points(a, Data, n);
	for (i = 0; i < n; ++i)
	{
		naive_knn(i, Data, n, k, iz, 0);
	}
	free(Data);
	Data = NULL;
}
//-----------------------------------------------------------------------------
//	Functions for and including the real seek
//-----------------------------------------------------------------------------
void gen_sub_boxes(box X, box *a, box *b, box *c, box *d, int verbose)
{
	if (verbose)
	{
		printf("parent box:\n");

		print_point(X.b_left);
		print_point(X.b_right);
		print_point(X.t_right);
		print_point(X.t_left);
	}

	a->b_left.x = X.b_left.x;
	a->b_left.y = (X.t_right.y + X.b_left.y) / 2;

	a->t_right.y = X.t_right.y;
	a->t_right.x = (X.t_right.x + X.b_left.x) / 2;

	b->b_left.x = (X.t_right.x + X.b_left.x) / 2;
	b->b_left.y = (X.t_right.y + X.b_left.y) / 2;

	b->t_right.y = X.t_right.y;
	b->t_right.x = X.t_right.x;

	c->b_left.x = X.b_left.x;
	c->b_left.y = X.b_left.y;

	c->t_right.y = b->b_left.y;
	c->t_right.x = b->b_left.x;

	d->b_left.x = (X.t_right.x + X.b_left.x) / 2;
	d->b_left.y = X.b_left.y;

	d->t_right.y = (X.t_right.y + X.b_left.y) / 2;
	d->t_right.x = X.t_right.x;




	a->b_right = c->t_right;
	a->t_left = X.t_left;

	b->t_left = a->t_right;
	b->b_right = d->t_right;

	c->b_right = d->b_left;
	c->t_left = a->b_left; 

	d->b_right = X.b_right;
	d->t_left = b->b_left;
	
	if (verbose)
	{
		printf("box a:\n----------\n");
		print_point(a->b_left);
		print_point(a->b_right);
		print_point(a->t_right);
		print_point(a->t_left);
		printf("box b:\n----------\n");
		print_point(b->b_left);
		print_point(b->b_right);
		print_point(b->t_right);
		print_point(b->t_left);
		printf("box c:\n----------\n");
		print_point(c->b_left);
		print_point(c->b_right);
		print_point(c->t_right);
		print_point(c->t_left);
		printf("box d:\n----------\n");
		print_point(d->b_left);
		print_point(d->b_right);
		print_point(d->t_right);
		print_point(d->t_left);
	}
}
//----------------------------------------------------------------------------
void add_control_entry(control_object *control, int *permutation, int current, int next, point *Data, int n, int verbose)
{

	if (verbose)
	{
		printf("looking at current = %d\n", current);
	}
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
		                                &control[next + 3].Box, verbose);

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
				e--;
			}
		}
		if (verbose)
		{
			printf("start = %d, end = %d\n", s, e);
		}
		control[next+count].start = head;
		control[next+count].end = s - 1;
		count++;
		head = s;
	}	
}
//----------------------------------------------------------------------------
int intersect(box B, circle C)
{
	if (is_in_box(C.center, B))
	{
		return 1;
	}
	if ((distance(B.t_right, C.center) < C.radius) || 
		(distance(B.t_left, C.center) < C.radius) ||
		(distance(B.b_right, C.center) < C.radius) ||
		(distance(B.b_left, C.center) < C.radius))
	{
		return 1;
	}
	point p = C.center;
	p.y += C.radius;
	if (is_in_box(p, B))
	{
		return 1;
	}
	p.y -= 2 * C.radius;
	if (is_in_box(p, B))
	{
		return 1;
	}
	p.y += C.radius;
	p.x += C.radius;
	if (is_in_box(p, B))
	{
		return 1;
	}
	p.x -= 2 * C.radius;
	if (is_in_box(p, B))
	{
		return 1;
	}
	return 0;	
}
//----------------------------------------------------------------------------
int find_parent_idx(point p, control_object *control)
{
	int current = 0, i, j;
	int leaf_idx = control[0].children[0];
	while(leaf_idx != 0)
	{
		if (is_in_box(p, control[leaf_idx].Box))
		{
			current = leaf_idx;
		}
		else if (is_in_box(p, control[leaf_idx + 1].Box))
		{
			current = leaf_idx + 1;
		}
		else if (is_in_box(p, control[leaf_idx + 2].Box))
		{
			current = leaf_idx + 2;
		}
		else
		{
			current = leaf_idx + 3;
		}
		leaf_idx = control[current].children[0];
	}
	return control[current].parent;
}
//----------------------------------------------------------------------------
void knn(int idx, point *Data, int points[], int length, int n, int k, int *iz, int verbose)
{
	pair* pairs;
	pairs = (pair*)malloc((length - 1) * sizeof(pair));
	// pair pairs[length];
	int ix = 0, j;
	for (j = 0; j < (length); ++j)
	{
		if (points[j] == idx)
		{
			continue;
		}
		pairs[ix].idx = points[j];
		pairs[ix].value = distance(Data[idx], Data[points[j]]);
		++ix;
	}
	if (verbose)
	{
		printf("%s\n", "accumulated distances.");
	}
	qsort(pairs, ix, sizeof(pair), compare);
	if (verbose)
	{
		printf("%s\n", "values sorted!");
	}
	for (j = 0; j < k; ++j)
	{
		iz[idx * k + j] = pairs[j].idx;
	}
	if (verbose)
	{
		printf("%s\n", "values imputed.");
	}
	free(pairs);
	pairs = NULL;
}
//----------------------------------------------------------------------------
void seek(double *a, int n, int k, int *iz, int verbose) 
{
	point *Data;
	Data = (point*) malloc(n * sizeof(point));
	arr_to_points(a, Data, n);
	// getchar();
	// print_point(Data[0]);
	// print_point(Data[1]);
	// printf("test distance = %f", distance(Data[0], Data[1]));
	// getchar();


	int i, j, l;
	int *permutation = (int*) malloc(n * sizeof(int));
	for (i = 0; i < n; i++)
	{
		permutation[i] = i;
	}

	control_object *Control;
	Control = (control_object*) malloc(BUFFER * sizeof (control_object));
	// control_object Control[BUFFER];
	
	Control[0].is_processed = 0, Control[0].start = 0, Control[0].end = (n - 1);
	Control[0].parent = -8, 
	Control[0].Box.b_left.x = 0, Control[0].Box.b_left.y = 0;
	Control[0].Box.t_left.x = 0, Control[0].Box.t_left.y = 1;
	Control[0].Box.t_right.x = 1, Control[0].Box.t_right.y = 1;
	Control[0].Box.b_right.x = 1, Control[0].Box.b_right.y = 0;
	if(verbose)
	{
		print_point(Control[0].Box.b_left);
	}
	if(verbose)
	{
		print_point(Control[0].Box.t_right);
	}

	int current = 0;
	int next = 1;

	while(current < next)
	{
		int e = Control[current].end;
		int s = Control[current].start;
		if(verbose)
		{
			printf("-------\ne = %d, s = %d, num in node = %d\n", e, s, e - s + 1);
		}
		if (((e - s + 1) > k)/* && (!Control[current].is_processed)*/) // if there are more than k points at this node.
		{	
			if(verbose)
			{
				printf("splitting a node!\n");
			}
			add_control_entry(Control, permutation, current, next, Data, n, verbose);
			next += 4;
		};
		Control[current].is_processed = 1;
		current++;
	}
	if (verbose)
	{
		getchar();
		printf("where do things start and end?\n");
		// for (i = 0; i < next; ++i)
		// {
		// 	printf("node %d, starts at %d, ends at %d.\n", i+1, Control[i].start, Control[i].end);
		// }
		getchar();
	}
	for (i = 0; i < n; ++i)
	{
		if(verbose)
		{
			printf("-------------------------\n%s%d\n", "Searching for neighbors of point ", i+1);
		}
		int parent_idx = find_parent_idx(Data[i], Control);
		circle C;
		C.center = Data[i];
		C.radius = get_radius(Data[i], Control[parent_idx].Box);

		int points[n];
		int to_consider[n];
		for (j = 0; j < n; ++j)
		{	
			to_consider[j] = 0;
		}
		int length = 0;
		if(verbose)
		{
			printf("parent idx = %d, index starts: %d, index ends: %d\n", parent_idx,  Control[parent_idx].start, Control[parent_idx].end);
			// printf("Child indices of: ");
			printf("we've considered points:\n");
		}
		int h;
		// for (h = 0; h < 4; ++h)
		// {
		// 	if(verbose)
		// 	{
		// 		printf("%d  ", Control[parent_idx].children[h]);
		// 	}
		// }
		
		for (j = Control[parent_idx].start; j <= Control[parent_idx].end; ++j)
		{
			points[length++] = permutation[j];

			if (verbose)
			{
				printf("id: %d\n", permutation[j]);
			}
			to_consider[permutation[j]] = 1;
		}
		// getchar();
		// print_int_matrix(permutation, 1, n);
		// getchar();


		for (j = 0; j < next; ++j)
		{
			if ((Control[j].children[0] != 0) ||
				(Control[j].parent == parent_idx) ||
				(j == parent_idx))
			{
				continue;
			}
			// if(verbose)
			// {
			// 	printf("search index j = %d\n", j);
			// }
			if (intersect(Control[j].Box, C))
			{
				for (l = Control[j].start; l <= Control[j].end; l++)
				{
					if (to_consider[permutation[l]] == 0)
					{
						if (verbose)
						{
							printf("id: %d\n", permutation[l]);
						}
						
						points[length++] = permutation[l];
						to_consider[permutation[l]] = 1;
					}
				}
			}
		}
		if(verbose)
		{
			printf("length is %d\n", length);
		}
		knn(i, Data, points, length, n, k, iz, verbose);

	}
	free(Data);
	Data = NULL;
	free(permutation);
	permutation = NULL;
	// free(Control);
	// Control = NULL;
}
//----------------------------------------------------------------------------