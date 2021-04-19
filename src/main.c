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
	/**------------------------------
	 * Parameters and data structures
	 * ------------------------------
	 */
	int N = 100;   /* Number of particles */
	int x0 = 1;  /* Number of particles in box 0 at t = 0 */
	int T = 10000; /* Number of time steps of the simulation */

	bool box[N]; /* box[i] = box to wich particle i belongs to */
	int frequency[N + 1]; 	/* frequency[i] = num of times box 0 has had i particles */
	int cnt_occupancy = x0;

	/**-----------------------------------
	 * Initialize box and frequency vector
	 * -----------------------------------
	 */
	for (int i = 0; i < N; i++) {
		if (i < x0) {
			box[i] = 0;
		} else {
			box[i] = 1;
		}
	}
	for (int i = 0; i <= N; i++) {
		frequency[i] = 0;
	}
	assert(x0 >= 0 && x0 <= N);
	frequency[x0] = 1;

	/**----------------------------------
	 * Initialize random number generator
	 * ----------------------------------
	 */
	gsl_rng *r = gsl_rng_alloc(gsl_rng_taus);
	assert(r != NULL);
	gsl_rng_set(r, time(NULL));		/* Set seed */
	
	/**-----------------------------------
	 * Run simulation
	 * -----------------------------------
	 */
	int chosen;
	for (int t = 0; t < T; t++) {
		chosen = gsl_rng_uniform_int(r, N);	/* Pick a particle at random */
		if (box[chosen] == 0) {
			cnt_occupancy--;
		} else {
			cnt_occupancy++;
		}
		assert(cnt_occupancy >= 0 && cnt_occupancy <= N);
		box[chosen] = !box[chosen]; 	/* Change box */
		frequency[cnt_occupancy]++;
	}
	/* Always remember to free allocated memory!! */
	gsl_rng_free(r);

	FILE * f_freq;
	f_freq = fopen("frequency.csv", "w");
	assert(f_freq != NULL);
	for (int i = 0; i <= N; i++) {
		fprintf(f_freq, "%d\n", frequency[i]);
	}
	fclose(f_freq);
}
