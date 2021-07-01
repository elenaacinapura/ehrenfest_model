#include "functions.h"

#include <linear_algebra/array_routines.h>
#include <linear_algebra/blas_wrappers.h>

#include <assert.h>
#include <cblas.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <time.h>

#define EPS 5e-4

int main() {
	/*************************************************************************
	 * PARAMETERS
	 * ***********************************************************************
	 * N = number of particles
	 * SAMPLING_STEPS = number of time steps of the sampling phase
	 * x0 = number of particles in box 0 at t = 0
	 * ***********************************************************************
	 */
	int N = 100;
	int x0 = N;

	/*************************************************************************
	 * VARIABLES
	 * ***********************************************************************
	 * box[i] 				box where particle i is
	 * x 					state of the system, number of particle in box 0
	 * distribution[i] 		number of instants spent in state i
	 * normalized_dist[i]	same as distribution[i], but normalized
	 * last_hit_time[i] 	last hitting time of state i
	 * avg_ret_time[i] 		average return time to state i
	 * times_returned[i]	number of times that system returned to state i
	 * ***********************************************************************
	 */
	bool box[N];
	int x;
	double distribution[N + 1];
	double normalized_dist[N + 1];
	double limiting_dist[N + 1];
	int last_hit_time[N + 1];
	double avg_return_time[N + 1];
	int times_returned[N + 1];

	FILE *f;

	/*************************************************************************/
	/* Initialize random number generator */
	/*************************************************************************/
	gsl_rng *r = gsl_rng_alloc(gsl_rng_taus);
	assert(r != NULL);
	gsl_rng_set(r, time(NULL)); /* Set seed */

	/*************************************************************************/
	/* Set initial conditions */
	/*************************************************************************/
	initialize_simulation(N, x0, box, &x, distribution, last_hit_time, avg_return_time, times_returned);
	fill_limiting_distribution(N, limiting_dist);

	/*************************************************************************/
	/* Welcome print */
	/*************************************************************************/
	printf("\n=========================================================\n");
	printf("EHRENFEST URN PROBLEM SIMULATION\n");
	printf("=========================================================\n");
	printf("Starting the simulation.\n");

	/*************************************************************************
	 * SIMULATION
	 * **********************************************************************/
	int t = -1;
	double distance_from_limit;

	/*********************************
	 * Converging phase
	 *********************************
	 */
	printf("**************************************\n");
	printf("Converging phase\n");

	f = fopen("distance.csv", "w");
	fprintf(f, "t\tdist\n");

	do {
		t++;

		/* Do the Ehrenfest step */
		simulation_step(N, box, &x, r);

		/* Update distributions*/
		distribution[x] += 1.0;
		vec_copy(N + 1, distribution, normalized_dist);
		normalize(N, normalized_dist);

		/* Update hitting time */
		if (last_hit_time[x] != -1) {
			avg_return_time[x] += t - last_hit_time[x];
			times_returned[x]++;
		}
		last_hit_time[x] = t;

		/* Calculate the distance from the limit distribution */
		distance_from_limit = distance(N, normalized_dist, limiting_dist);

		/* Print info about the current iteration */
		if (t % 100 == 0) {
			printf("\r\tt = %d delta = %.3lf", t, distance_from_limit);
		}
		fprintf(f, "%d\t%lf\n", t, distance_from_limit);
		

	} while (distance_from_limit > EPS); 
	fclose(f);

	for (int i = 0; i <= N; i++) {
		if (times_returned[i] == 0) {
			avg_return_time[i] = __DBL_MAX__;
		}
		avg_return_time[i] = (double)(avg_return_time[i]) / times_returned[i];
	}

	
	f = fopen("limiting_dist.csv", "w");
	fprintf(f, "i\tval\n");
	for (int i = 0; i <= N; i++) {
		fprintf(f, "%d\t%lf\n", i, limiting_dist[i]);
	}
	fclose(f);

	f = fopen("distribution.csv", "w");
	fprintf(f, "i\tval\n");
	for (int i = 0; i <= N; i++) {
		fprintf(f, "%d\t%lf\n", i, normalized_dist[i]);
	}
	fclose(f);

	f = fopen("avg_return_time_sampled.csv", "w");
	fprintf(f, "i\tval\n");
	for (int i = 0; i <= N; i++) {
		fprintf(f, "%d\t%lf\n", i, avg_return_time[i]);
	}
	fclose(f);

	f = fopen("avg_return_time_th.csv", "w");
	fprintf(f, "i\tval\n");
	for (int i = 0; i <= N; i++) {
		fprintf(f, "%d\t%lf\n", i, 1.0 / limiting_dist[i]);
	}
	fclose(f);

	printf("\n**************************************\n");
	printf("End of the simulation.\nShowing the plots.\n");
	printf("=========================================================\n");
}