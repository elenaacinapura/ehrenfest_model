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
	 * T = number of time steps
	 * x0 = number of particles in box 0 at t = 0
	 * ***********************************************************************
	 */
	int N = 30;
	int x0 = N;
	int T = 1e7;

	/*************************************************************************
	 * VARIABLES
	 * ***********************************************************************
	 * box[i] 				box where particle i is
	 * x 					state of the system, number of particle in box 0
	 * visits[i] 			number of instants spent in state i
	 * visit_freq[i]		same as visits[i], but normalized
	 * last_visit[i] 		last hitting time of state i
	 * recurrence_time[i] 	average return time to state i
	 * times_returned[i]	number of times that system returned to state i 
	 * 						(basically visit[i]-1)
	 * ***********************************************************************
	 */
	bool box[N];
	int x;
	double visits[N + 1];
	double visit_freq[N + 1];
	double limiting_dist[N + 1];
	int last_visit[N + 1];
	double recurrence_time[N + 1];
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
	initialize_simulation(N, x0, box, &x, visits, last_visit, recurrence_time, times_returned);
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

	f = fopen("distance.csv", "w");
	fprintf(f, "t\tdist\n");

	do {
		t++;

		/* Do the Ehrenfest step */
		simulation_step(N, box, &x, r);

		/* Update distributions*/
		visits[x] += 1.0;
		vec_copy(N + 1, visits, visit_freq);
		normalize(N, visit_freq);

		/* Update hitting time */
		if (last_visit[x] != -1) {
			recurrence_time[x] += t - last_visit[x];
			times_returned[x]++;
		}
		last_visit[x] = t;

		/* Calculate the distance from the limit distribution */
		distance_from_limit = distance(N, visit_freq, limiting_dist);

		/* Print info about the current iteration */
		if (t % 100 == 0) {
			printf("\r\tt = %d D = %.3lf", t, distance_from_limit);
		}
		fprintf(f, "%d\t%lf\n", t, distance_from_limit);

	} while (t < T); // distance_from_limit > EPS
	fclose(f);

	/* Final calculation of the mean recurrence time */
	for (int i = 0; i <= N; i++) {
		if (times_returned[i] < 1) {
			recurrence_time[i] = -1;
		} else {
			recurrence_time[i] = (double)(recurrence_time[i]) / times_returned[i];
		}
	}

	f = fopen("limiting_dist.csv", "w");
	fprintf(f, "i\tval\n");
	for (int i = 0; i <= N; i++) {
		fprintf(f, "%d\t%lf\n", i, limiting_dist[i]);
	}
	fclose(f);

	f = fopen("visit_freq.csv", "w");
	fprintf(f, "i\tval\n");
	for (int i = 0; i <= N; i++) {
		fprintf(f, "%d\t%lf\n", i, visit_freq[i]);
	}
	fclose(f);

	f = fopen("recurrence_time_sampled.csv", "w");
	fprintf(f, "i\tval\n");
	for (int i = 0; i <= N; i++) {
		fprintf(f, "%d\t%lf\n", i, recurrence_time[i]);
	}
	fclose(f);

	f = fopen("recurrence_time_th.csv", "w");
	fprintf(f, "i\tval\n");
	for (int i = 0; i <= N; i++) {
		fprintf(f, "%d\t%lf\n", i, 1.0 / limiting_dist[i]);
	}
	fclose(f);

	printf("\n**************************************\n");
	printf("End of the simulation.\nShowing the plots.\n");
	printf("=========================================================\n");
}