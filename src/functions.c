#include "functions.h"

#include <linear_algebra/blas_wrappers.h>

#include <assert.h>
#include <cblas.h>
#include <math.h>
#include <stdio.h>



void initialize_P(int N, double P[N + 1][N + 1]) {
    int DIM = N + 1;
	for (int i = 0; i < DIM; i++) {
		if (i > 0) {
            P[i-1][i] = (double)(i)/N;  /* Probability of transition i->i-1 */
        }
        if (i < DIM - 1) {
            P[i + 1][i] = (double) (N - i)/N;   /* Probability of transition i->i+1 */
        }
	}
}

void theoretical_prediction (int N, double f[N+1], double P[N+1][N+1], int num_steps) {
    int DIM = N + 1;
    for (int k = 0; k < num_steps; k++) {
        double res[N+1];
        mat_vec_mul(DIM, DIM, (double *)P, f, res);
        vec_copy(DIM, res, f);
    }
}