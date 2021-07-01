#ifndef __FUNCTIONS_H__
#define __FUNCTIONS_H__

#include <stdbool.h>
#include <gsl/gsl_rng.h>

void initialize_simulation(int N, int x0, bool box[], int *x, double distribution[], int last_hit_time[], double avg_return_time [], int times_returning[]);

void simulation_step(int N, bool box[], int *x, gsl_rng *r);

double distance(int N, double d[], double ld[]);

void normalize (int N, double d[]);

void fill_limiting_distribution(int N, double ld[]);

void initialize_P(int N, double P[N + 1][N + 1], bool modified);

void theoretical_prediction(int N, double f[N+1], double P[N+1][N+1], int num_steps);

#endif
