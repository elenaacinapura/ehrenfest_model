#include "functions.h"

#include <linear_algebra/array_routines.h>

#include <assert.h>
#include <cblas.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <time.h>

int main () {
    /*************************************************************************
	 * PARAMETERS
	 * ***********************************************************************
     * N = number of particles
     * T = number of timesteps of the simulation
     * N_REP = number of repetitions of the same simulation
     * x0 = number of particles in box 0 at t = 0
     * ***********************************************************************
	 */
    int N = 10;	
	int T = 1000;
	int N_REP = 10000;
    int x0 = 1;

    assert(N > 0 && T >= 0 && N_REP > 0);
	assert(x0 >= 0 && x0 <= N);

    /*************************************************************************
     * SIMULATION
     * ***********************************************************************
     * box[i] = box where particle i is
     * occupation_cnt = current number of particle in box 0
     * count[i] = how many times box 0 had i particles after T timesteps
     * ***********************************************************************
     */
    bool box[N];
	int occupation_cnt;
    int count[N + 1];
    /**--------------------
     * Initialize count[]
     * --------------------
     */
	for (int i = 0; i <= N; i++) {
		count[i] = 0;
	}
    /**-----------------------------------
     * Initialize random number generator
     * ----------------------------------- 
     */
	gsl_rng *r = gsl_rng_alloc(gsl_rng_taus);
	assert(r != NULL);
	gsl_rng_set(r, time(NULL)); /* Set seed */
    /**---------------------------
     * Run simulation N_REP times
     * ---------------------------
     */
    for (int rep = 0; rep < N_REP; rep++) {
        /**-----------------------
         * Set initial conditions 
         * -----------------------
         */
        occupation_cnt = x0;
        for (int i = 0; i < x0; i++) {
			box[i] = 0;
		}
		for (int i = x0; i < N; i++) {
			box[i] = 1;
		}
        /**--------------------
         * Run for T timesteps 
         * --------------------
         */
        int chosen;
		for (int t = 0; t < T; t++) {
			chosen = gsl_rng_uniform_int(r, N); /* Pick a particle at random */
			if (box[chosen] == 0) {
				occupation_cnt--;
			} else {
				occupation_cnt++;
			}
			assert(occupation_cnt >= 0 && occupation_cnt <= N);
			box[chosen] = !box[chosen]; /* Change box */
		}
		count[occupation_cnt]++;
    }

    /* Always remember to free allocated memory!! */
	gsl_rng_free(r);

    /*************************************************************************
     * THEORETICAL PREDICTION
     * ***********************************************************************
     * f[i] = probability that box 0 has i particles at the given time
     * P[i][j] = probability of the transition from state j to state i
     * ***********************************************************************
     */
    double f[N + 1];
	double P[N + 1][N + 1];

	for (int i = 0; i <= N; i++) {
		f[i] = 0.0;
	}
	f[x0] = 1.0;

	initialize_P(N, P, 0);

	theoretical_prediction(N, f, P, T);

    /*************************************************************************
     * Print to file
     * ***********************************************************************
     */
    /* Parameters */
    FILE *f_parameters;
	f_parameters = fopen("basic_parameters.csv", "w");
	assert(f_parameters != NULL);
	fprintf(f_parameters, "N\tx0\tT\tn_rep\n%d\t%d\t%d\t%d\n", N, x0, T, N_REP);
	fclose(f_parameters);

    /* Count */
    FILE * f_count;
    f_count = fopen("basic_count.csv", "w");
    assert(f_count != NULL);
    fprintf(f_count, "count\n");
    for ( int i = 0; i <= N; i++) {
		fprintf(f_count, "%d\n", count[i]);
	}
    fclose(f_count);

    /* Prediction */
	FILE *f_theo;
	f_theo = fopen("basic_prediction.csv", "w");
	assert(f_theo != NULL);
	fprintf(f_theo, "prob\n");
	for (int i = 0; i <= N; i++) {
		fprintf(f_theo, "%e\n", f[i]);
	}
	fclose(f_theo);
}