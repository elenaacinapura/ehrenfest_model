#include "functions.h"

#include <linear_algebra/blas_wrappers.h>

#include <assert.h>
#include <cblas.h>
#include <gsl/gsl_rng.h>
#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>

void initialize_simulation(int N, int x0, bool box[], int *x, double distribution[]) {
	for (int i = 0; i <= N; i++) {
		distribution[i] = 0.0;
	}
	for (int i = 0; i < N; i++) {
		if (i < x0) {
			box[i] = 0;
		} else {
			box[i] = 1;
		}
	}
	*x = x0;
	distribution[x0] += 1.0;
}

void simulation_step(int N, bool box[], int *x, gsl_rng *r) {
	int random_part = gsl_rng_uniform_int(r, N); /* Pick a particle at random */
	int random_box = gsl_rng_uniform_int(r, 2);  /* Pick a box at random */
	// printf("part = %d, box = %d, box(p) = %d, x = %d\n", random_part, random_box, box[random_part], *x);
	if (box[random_part] != random_box) {
		if (random_box == 0) {
			(*x)++;
		} else {
			(*x)--;
		}
		box[random_part] = !box[random_part]; /* Change box */
	}
}

double distance(int N, double d[], double ld[]) {
	double sum = 0.0;
	for (int i = 0; i <= N; i++) {
		sum += fabs(d[i] - ld[i]);
	}
	return sum;
}

void normalize(int N, double d[]) {
	double sum = 0.0;
	for (int i = 0; i <= N; i++) {
		sum += d[i];
	}
	for (int i = 0; i <= N; i++) {
		d[i] /= sum;
	}
}

void initialize_P(int N, double P[N + 1][N + 1], bool modified) {
	int DIM = N + 1;
	for (int i = 0; i < DIM; i++) {
		for (int j = 0; j < DIM; j++) {
			if (modified) {
				if (i == j - 1) {
					P[i][j] = 0.5 * (double)(j) / N; /* Probability of transition j->j-1 */
				} else if (i == j + 1) {
					P[i][j] = 0.5 * (double)(N - j) / N; /* Probability of transition j->j+1 */
				} else if (i == j) {
					P[i][j] = 0.5 * (double)(j) / N + 0.5 * (double)(N - j) / N;
				} else {
					P[i][j] = 0.0;
				}
			} else {
				if (i == j - 1) {
					P[i][j] = (double)(j) / N; /* Probability of transition j->j-1 */
				} else if (i == j + 1) {
					P[i][j] = (double)(N - j) / N; /* Probability of transition j->j+1 */
				} else {
					P[i][j] = 0.0;
				}
			}
		}
	}
}

unsigned long binomial(unsigned long n, unsigned long k) {
	unsigned long c = 1, i;
	if (k > n - k) { // take advantage of symmetry
		k = n - k;
	}
	for (i = 1; i <= k; i++, n--) {
		if (c / i > UINT_MAX / n) { // return 0 on overflow
			return 0;
		}
		c = c / i * n + c % i * n / i; // split c * n / i into (c / i * i + c % i) * n / i
	}
	return c;
}

void fill_limiting_distribution(int N, double ld[]) {
	for (int i = 0; i <= N; i++) {
		ld[i] = pow(0.5, N) * binomial(N, i);
	}
}

void theoretical_prediction(int N, double f[N + 1], double P[N + 1][N + 1], int num_steps) {
	int DIM = N + 1;
	double res[DIM];
	for (int k = 0; k < num_steps; k++) {
		mat_vec_mul(DIM, DIM, (double *)P, f, res);
		vec_copy(DIM, res, f);
	}
}