#include "functions.h"

#include <linear_algebra/array_routines.h>

#include <assert.h>
#include <cblas.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <time.h>

int main() {
	/**-----------
	 * Parameters
	 * -----------
	 */
	int N = 10;	 /* Number of particles */
	int T = 100;	 /* Number of time steps of the simulation */
	int n_rep = 10000; /* Number of repetitions of the simulation */

	int x0 = 1;	 /* Number of particles in box 0 at t = 0 */
	assert(x0 >= 0 && x0 <= N);

	FILE *f_parameters;
	f_parameters = fopen("parameters.csv", "w");
	assert(f_parameters != NULL);
	fprintf(f_parameters, "N\tx0\tT\tn_rep\n%d\t%d\t%d\t%d\n", N, x0, T, n_rep);
	fclose(f_parameters);

	/****************************************************************************************
	 * SIMULATION
	 ***************************************************************************************/
	bool box[N];		/* box[i] = box to wich particle i belongs to */
	int occ_freq[N + 1];	/* occ_freq[i] = num of times box 0 has had i particles */
	int cnt_occupation = x0; /* Current number of particles in box 0 */

	int final_state_counting[N + 1]; /* [i] = how many times box 0 had i particles after T timesteps */

	/**-----------------------------------
	 * Initialize final_state_counting
	 * -----------------------------------
	 */
	for (int i = 0; i <= N; i++) {
		final_state_counting[i] = 0;
	}

	/**----------------------------------
	 * Initialize random number generator
	 * ----------------------------------
	 */
	gsl_rng *r = gsl_rng_alloc(gsl_rng_taus);
	assert(r != NULL);
	gsl_rng_set(r, time(NULL)); /* Set seed */

	/**---------------
	 * Run rimulation
	 * ---------------
	 */
	for (int rep = 0; rep < n_rep; rep++) {
		/**-----------------------------------
		 * Initialize box and occ_freq vector
		 * -----------------------------------
		 */
		cnt_occupation = x0;
		for (int i = 0; i < N; i++) {
			if (i < x0) {
				box[i] = 0;
			} else {
				box[i] = 1;
			}
		}
		for (int i = 0; i <= N; i++) {
			occ_freq[i] = 0;
		}
		occ_freq[x0] = 1;

		/**---------------------
		 * Run round simulation
		 * ---------------------
		 */
		int chosen;
		for (int t = 0; t < T; t++) {
			chosen = gsl_rng_uniform_int(r, N); /* Pick a particle at random */
			if (box[chosen] == 0) {
				cnt_occupation--;
			} else {
				cnt_occupation++;
			}
			assert(cnt_occupation >= 0 && cnt_occupation <= N);
			box[chosen] = !box[chosen]; /* Change box */
			occ_freq[cnt_occupation]++;
		}

		final_state_counting[cnt_occupation]++;
	}

	/* Always remember to free allocated memory!! */
	gsl_rng_free(r);

	/**-------------------------------------------
	 * Print frequencies and final state counting
	 * -------------------------------------------
	 */
	FILE *f_freq, *f_final_state;
	f_freq = fopen("occ_freq.csv", "w");
	f_final_state = fopen("final_state.csv", "w");
	assert(f_freq != NULL && f_final_state != NULL);

	fprintf(f_freq, "freq\n");
	for (int i = 0; i <= N; i++) {
		fprintf(f_freq, "%lf\n", (double)(occ_freq[i])/n_rep);
	}

	fprintf(f_final_state, "count\n");
	for ( int i = 0; i <= N; i++) {
		fprintf(f_final_state, "%d\n", final_state_counting[i]);
	}

	fclose(f_freq);



	/****************************************************************************************
	 * THEORETICAL PREDICTION
	 ***************************************************************************************/

	/* Fill the initial distribution */
	double f0[N + 1];
	double P[N + 1][N + 1];

	for (int i = 0; i <= N; i++) {
		f0[i] = 0.0;
	}
	f0[x0] = 1.0;

	initialize_P(N, P, 0);

	theoretical_prediction(N, f0, P, T);

	/* Print */
	FILE *f_theo;
	f_theo = fopen("theoretical.csv", "w");
	assert(f_theo != NULL);
	fprintf(f_theo, "freq\n");
	for (int i = 0; i <= N; i++) {
		fprintf(f_theo, "%e\n", f0[i]);
	}
	fclose(f_theo);
}
