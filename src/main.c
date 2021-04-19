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
	int N = 3;   /* Number of particles */
	int x0 = 1;  /* Number of particles in box 0 at t = 0 */
	int T = 100; /* Number of time steps of the simulation */

	bool box[N]; /* box[i] = box to wich particle i belongs to */
	int cnt_occupancy = 0;

	/**----------------------
	 * Initialize box vector
	 * ----------------------
	 */
	for (int i = 0; i < N; i++) {
		if (i < x0) {
			box[i] = 0;
		} else {
			box[i] = 1;
		}
	}

	/**-----------------------------------
	 * Initialize random number generator
	 * -----------------------------------
	 */
	gsl_rng *r = gsl_rng_alloc(gsl_rng_taus);
	gsl_rng_set(r, time(NULL));		/* Set seed */
	
	/**-----------------------------------
	 * Run simulation
	 * -----------------------------------
	 */
	for (int t = 0; t < T; t++) {
		int chosen = gsl_rng_uniform_int(r, N);	/* Pick a particle at random */
		if (box[chosen] == 0) {
			cnt_occupancy--;
		} else {
			cnt_occupancy++;
		}
		box[chosen] = !box[chosen]; 	/* Change box */
	}

	/* Always remember to free allocated memory!! */
	gsl_rng_free(r);
}
