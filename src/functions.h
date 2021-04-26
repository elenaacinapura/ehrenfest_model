#ifndef __FUNCTIONS_H__
#define __FUNCTIONS_H__

#include <stdbool.h>
/** 
 * @brief Fill stochastic matrix P for an Ehrenfest problem.
 * 
 * @param N Number of particles in the model
 * @param P Stochastic matrix (N + 1)*(N + 1) to be initialized
 * @param modified Flag to indicate whether to consider the modified Ehrenfest problem
 * 
 * @return
 */
void initialize_P(int N, double P[N + 1][N + 1], bool modified);

/** 
 * @brief Evaluate the probability distribution after given timesteps.
 * 
 * @param N Number of particles in the model
 * @param f Initial distribution that will be let evolve
 * @param P Stochastic matrix (N + 1)*(N + 1)
 * @param num_steps Number of timesteps for the evolution
 * 
 * @return
 */
void theoretical_prediction(int N, double f[N+1], double P[N+1][N+1], int num_steps);

#endif
