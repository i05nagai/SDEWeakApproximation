#ifndef _MLMC_H_
#define _MLMC_H_

#include <sde_wa/sde_wa.h>

typedef struct mlmc_l{
	int l;

	int N;
	int newN;
	double sumY;
	double sumYY;
	double V;

	struct mlmc_l *prev;
	struct mlmc_l *next;
} MLMC_L;

typedef struct sde_mlmc{
	SDE_WA_SLTN *sl;
	MLMC_L *list;

	double *init_y;
	double T;
	double (*payoff_function)(double *x, void *params);

	int (*level_l_estimation)(struct sde_mlmc *mlmc, MLMC_L *estimator, int l, int M, int N);
	double P;
} SDE_MLMC;


void trans_rand_unif_to_normal(double *unif_rand, double *normal_rand, int n);
SDE_MLMC *alloc_SDE_MLMC(SDE_WA_SLTN *sl);
void free_SDE_MLMC(SDE_MLMC *mlmc);
MLMC_L *alloc_MLMC_L(void);
void free_MLMC_L(MLMC_L *mlmc_l);
void init_level_l_estimator(MLMC_L *mlmc_l);
void add_level_l_estimator(MLMC_L *mlmc_l);

int mlmc_algorithm(SDE_MLMC *mlmc, int M, double eps, int extraporation_flag);
#endif
