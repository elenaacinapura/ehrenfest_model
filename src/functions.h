#ifndef __FUNCTIONS_H__
#define __FUNCTIONS_H__

/** 
 * @brief Fill stochastic matrix P
 *
 */
void initialize_P(int N, double P[N + 1][N + 2]);

void theoretical_prediction(int N, double f[N+1], double P[N+1][N+1], int num_steps);

#endif