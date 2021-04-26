#include "functions.h"

#include <linear_algebra/blas_wrappers.h>

#include <assert.h>
#include <cblas.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>

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

void theoretical_prediction(int N, double f[N + 1], double P[N + 1][N + 1], int num_steps) {
	int DIM = N + 1;
	double res[DIM];
	for (int k = 0; k < num_steps; k++) {
		mat_vec_mul(DIM, DIM, (double *)P, f, res);
		vec_copy(DIM, res, f);
	}
}