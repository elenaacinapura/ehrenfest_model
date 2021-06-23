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

int main() {
	/*************************************************************************
	 * PARAMETERS
	 * ***********************************************************************
	 * N = number of particles
	 * CONVERGENCE_REP = number of timesteps of the convergence phase
	 * SAMPLING_REP = number of repetitions of the sampling phase
	 * x0 = number of particles in box 0 at t = 0
	 * ***********************************************************************
	 */
	int N = 10;
	int CONVERGENCE_REP = 100;
	int SAMPLING_REP = 10000;
	int x0 = 0;

	/*************************************************************************
	 * VARIABLES
	 * ***********************************************************************
	 * box[i] 				box where particle i is
	 * x 					state of the system, number of particle in box 0
	 * distribution[i] 		time spent with i particles in box 0
	 * last_hit_time[i] 	last hitting time of state i
	 * avg_ret_time[i] 		average return time to state i
	 * hitting_cnt[i]		number of times the system hit the state i
	 * ***********************************************************************
	 */
	bool box[N];
	int x;
	double distribution[N + 1];
	double normalized_dist[N + 1];
	double limiting_dist[N + 1];
	
	/* Initialize random number generator */
	gsl_rng *r = gsl_rng_alloc(gsl_rng_taus);
	assert(r != NULL);
	gsl_rng_set(r, time(NULL)); /* Set seed */

	/* Set initial conditions */
	initialize_simulation(N, x0, box, &x, distribution);
	fill_limiting_distribution(N, limiting_dist);

	/* Welcome */
	printf("\n***********************************************************\n");
	printf("EHRENFEST URN PROBLEM SIMULATION\n");
	printf("***********************************************************\n");
	printf("Simulating...\n\n");

	/*************************************************************************
	 * START THE SIMULATION
	 * **********************************************************************/
	int t = 0;
	double distance_from_limit;
	do {
		printf("\rt = %d", t);
		simulation_step(N, box, &x, r);
		distribution[x] += 1.0;
		vec_copy(N+1, distribution, normalized_dist);
		normalize(N, normalized_dist);
		distance_from_limit = distance(N, normalized_dist, limiting_dist);
		t++;
		printf(", dist = %lf", distance_from_limit);
	} while (t < 100);

	printf("\nConvergence time = %d", t);

	FILE *f;
	f = fopen("limiting_dist.csv", "w");
	fprintf(f, "i\tval\n");
	for (int i = 0; i <= N; i++) {
		fprintf(f, "%d\t%lf\n", i, limiting_dist[i]);
	}

	f = fopen("distribution.csv", "w");
	fprintf(f, "i\tval\n");
	for (int i = 0; i <= N; i++) {
		fprintf(f, "%d\t%lf\n", i, normalized_dist[i]);
	}
}