#include "../functions.h"

#include <linear_algebra/array_routines.h>

#include <assert.h>
#include <cblas.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <time.h>

int main() {
	/* Welcome */
	printf("\n***********************************************************\n");
	printf("EHRENFEST URN PROBLEM SIMULATION\n");
	printf("***********************************************************\n");
	printf("Simulating...\n\n");

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
	int T_single = 100000000;
	int T_mult = 100;
	int N_REP = 10000;
	int x0 = 1;

	assert(N > 0 && T_mult >= 0 && T_single >= 0 && N_REP > 0);
	assert(x0 >= 0 && x0 <= N);

	/*************************************************************************
	 * SIMULATION
	 * ***********************************************************************
	 * box[i] = box where particle i is
	 * occupation_cnt = current number of particle in box 0
	 * count[i] = how many times box 0 had i particles after T timesteps
	 * time_spent[i] = time spent with i particles in box 0
	 * last_hit_time[i] = last hitting time of state i
	 * avg_ret_time[i] = average return time to state i
	 * hitting_cnt = number of times the system hit the state i
	 * ***********************************************************************
	 */
	bool box[N];
	int occupation_cnt;
	int count[N + 1];
	int time_spent[N + 1];
	int last_hit_time[N + 1];
	double avg_ret_time[N + 1];
	int hitting_cnt[N + 1];
	/**-----------------------------------
	 * Initialize random number generator
	 * -----------------------------------
	 */
	gsl_rng *r = gsl_rng_alloc(gsl_rng_taus);
	assert(r != NULL);
	gsl_rng_set(r, time(NULL)); /* Set seed */

	/**************************************************************
	 * ONE-REP SIMULATION
	 **************************************************************
	 */
	/**-----------------------
	 * Set initial conditions
	 * -----------------------
	 */
	for (int i = 0; i <= N; i++) {
		count[i] = 0;
		time_spent[i] = 0;
		avg_ret_time[i] = 0;
		hitting_cnt[i] = 0;
		last_hit_time[i] = 0;
	}
	for (int i = 0; i < x0; i++) {
		box[i] = 0;
	}
	for (int i = x0; i < N; i++) {
		box[i] = 1;
	}
	occupation_cnt = x0;
	time_spent[x0]++;
	last_hit_time[x0] = 0;
	hitting_cnt[x0] = 1;
	/**--------------------
	 * Run for T timesteps
	 * --------------------
	 */
	int random_part, random_box;
	for (int t = 0; t < T_single; t++) {
		random_part = gsl_rng_uniform_int(r, N); /* Pick a particle at random */
		random_box = gsl_rng_uniform_int(r, 2);	 /* Pick a box at random */
		if (box[random_part] != random_box) {
			if (random_box == 0) {
				occupation_cnt++;
			} else {
				occupation_cnt--;
			}
			box[random_part] = !box[random_part]; /* Change box */
		}
		assert(occupation_cnt >= 0 && occupation_cnt <= N);
		time_spent[occupation_cnt]++;
		hitting_cnt[occupation_cnt]++;

		if (hitting_cnt[occupation_cnt] > 1) {
			avg_ret_time[occupation_cnt] += (double)(t - last_hit_time[occupation_cnt]);
		}
		last_hit_time[occupation_cnt] = t;
	}

	for (int i = 0; i < N + 1; i++) {
		if (hitting_cnt[i] > 1) {
			avg_ret_time[i] /= hitting_cnt[i] - 1;
		}
	}

	/**************************************************************
	 * MULTIPLE-REP SIMULATION
	 **************************************************************
	 */
	/**--------------------
	 * Initialize vectors
	 * --------------------
	 */
	for (int i = 0; i <= N; i++) {
		count[i] = 0;
	}
	for (int rep = 0; rep < N_REP; rep++) {
		/**-----------------------
		 * Set initial conditions
		 * -----------------------
		 */
		for (int i = 0; i < x0; i++) {
			box[i] = 0;
		}
		for (int i = x0; i < N; i++) {
			box[i] = 1;
		}
		for (int i = 0; i < N + 1; i++) {
			hitting_cnt[i] = 0;
			last_hit_time[i] = 0;
		}
		occupation_cnt = x0;

		/**--------------------
		 * Run for T timesteps
		 * --------------------
		 */
		int random_part, random_box;
		for (int t = 0; t < T_mult; t++) {
			random_part = gsl_rng_uniform_int(r, N); /* Pick a particle at random */
			random_box = gsl_rng_uniform_int(r, 2);	 /* Pick a box at random */
			if (box[random_part] != random_box) {
				if (random_box == 0) {
					occupation_cnt++;
				} else {
					occupation_cnt--;
				}
				box[random_part] = !box[random_part]; /* Change box */
			}
			assert(occupation_cnt >= 0 && occupation_cnt <= N);
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

	initialize_P(N, P, 1);

	theoretical_prediction(N, f, P, T_mult);

	/*************************************************************************
	 * Print to file
	 * ***********************************************************************
	 */

	/* Parameters */
	FILE *f_parameters;
	f_parameters = fopen("basic_parameters.csv", "w");
	assert(f_parameters != NULL);
	fprintf(f_parameters, "N\tx0\tT_single\tT_mult\tn_rep\n%d\t%d\t%d\t%d\t%d\n", N, x0, T_single, T_mult, N_REP);
	fclose(f_parameters);

	/* Count */
	FILE *f_count;
	f_count = fopen("basic_count.csv", "w");
	assert(f_count != NULL);
	fprintf(f_count, "count\n");
	for (int i = 0; i <= N; i++) {
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

	/* Time spent */
	FILE *f_time;
	f_time = fopen("time_spent.csv", "w");
	assert(f_time != NULL);
	fprintf(f_time, "times\n");
	for (int i = 0; i <= N; i++) {
		fprintf(f_time, "%d\n", time_spent[i]);
	}
	fclose(f_time);

	/* Return time */
	FILE *f_ret;
	f_ret = fopen("return.csv", "w");
	assert(f_ret != NULL);
	fprintf(f_ret, "time\n");
	for (int i = 0; i <= N; i++) {
		fprintf(f_ret, "%e\n", avg_ret_time[i]);
	}
	fclose(f_ret);

	/* Bye-bye */
	printf("Simulation ended successfully.\n");
	printf("Now showing the plots.\n");
}