#include <stdio.h>
#include <stdlib.h>

struct xy_point {
	int x;
	int y;
};

int k1 = 256, k2 = 257, a1 = 57, a2 = 73, c1 = 1, c2 = 1;
int xn1 = 0, xn2 = 0;

int point_cmp(const struct xy_point *x1, const struct xy_point *x2)
{
	return (x1->x >= x2->x ? (x1->x == x2->x ? 0 : 1) : -1);
}

int bad_rng1();
int bad_rng2();
int improved_rng();

#define N_CHOICES 5
#define N_ELEMENTS 1000

int main(int argc, char **argv)
{
	int tests[] = { 2, 3, 5, 7, 13 };
	struct xy_point *list = NULL;
	list = malloc(N_ELEMENTS * sizeof(struct xy_point));
	if (!list) {
		fprintf(stderr, "can't allocate memory, exit\n");
		return -1;
	}
	int i = 0, j = 0, n = 0, k, tmp[2], tmp2[2];
	char buffer[512];
	FILE *fout;
	xn1 = 10;
	xn2 = 42;
	for (n = 0; n < N_CHOICES; n++) {
		for (i = 0; i < N_ELEMENTS; i++) {
			list[i].x = improved_rng(tests[n]);
			list[i].y = improved_rng(tests[n]);

			// For exercise (a) uncomment the following two lines, comment
			// out the above two lines and set N_CHOICES to 1

			/*(list[i]).x=bad_rng1();
			   (list[i]).y=bad_rng1(); */
		}

		qsort(list, N_ELEMENTS, sizeof(struct xy_point), &point_cmp);

		sprintf(buffer, "rng1_%d.dat", n);
		fout = fopen(buffer, "we");

		// Gnuplot
		// print 'rng1_i.dat' using 1:2
		for (i = 0; i < N_ELEMENTS; i++) {
			fprintf(fout, "%d %d \n", list[i].x, list[i].y);
		}
		fprintf(fout, "\n");
		fclose(fout);
	}

	return 0;
}

int bad_rng1()
{
	xn1 = ((a1 * xn1 + c1) % k1);
	return xn1;
}

int bad_rng2()
{
	xn2 = ((a2 * xn2 + c2) % k2);
	return xn2;
}

int improved_rng(int n)
{
	int x1 = 0, x2 = 1;
	if (n < 2)
		n = 2;
	while (x2 % n) {
		x2 = bad_rng2();
		x1 = bad_rng1();
	}
	return x1;
}
