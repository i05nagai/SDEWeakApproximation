/**
 * @file	mlmc_algo.c
 * @brief	an algorithm for MLMC estimation.

 *			MLMC Algorithm proposed by Giles. 
 * @date	09/Mar/2015 (Mon) 17:32 
 * @author	i05nagai
 */

#include <stdlib.h>
#include <math.h>
#include <mlmc.h>

/**
 * MLMC Algorithm proposed by Mike Giles.
 * MLMC Algorithm proposed by Mike Giles.
 * @param mlmc a pointer of a structure for MLMC.
 * @param M a base of the number of timesteps.
 * @param eps Mean Square Error.
 * @param extraporation_flag 0 or 1. Set 1 to use an extrapolation.
 * @return none.
 */
int mlmc_algorithm(SDE_MLMC *mlmc, int M, double eps, int extraporation_flag) {
	int converged_flag=0;
	int l;
	int L=-1;
	int n;
	int dNl=0;
	double con=0.0;
	double sum=0.0;
	double hl=0.0;
	MLMC_L *mlmc_ests;
	MLMC_L *mlmc_ests1;

	mlmc_ests = mlmc->list;
	while (!converged_flag) {
		L++;
		n = (int)pow(M,L);

		mlmc->level_l_estimation(mlmc, mlmc_ests, L, M, mlmc_ests->N);

		//estimate Vl
		sum = 0.0;
		mlmc_ests1 = mlmc->list;
		for (l=0; l<=L; l++) {
			mlmc_ests1->V = mlmc_ests1->sumYY/mlmc_ests1->N - 
				(mlmc_ests1->sumY/mlmc_ests1->N)*(mlmc_ests1->sumY/mlmc_ests1->N);
			hl = mlmc->T/pow(M,l);
			sum += sqrt(mlmc_ests1->V/hl);

			mlmc_ests1 = mlmc_ests1->next;
		}
		
		//define optimal Nl
		mlmc_ests1 = mlmc->list;
		for (l=0; l<=L; l++) {
			hl = mlmc->T/pow(M,l);
			mlmc_ests1->newN = ceil(2.0 * sqrt(mlmc_ests1->V * hl) * sum / (eps * eps));

			mlmc_ests1 = mlmc_ests1->next;
		}

		//update sample sum
		mlmc_ests1 = mlmc->list;
		for (l=0; l<=L; l++) {
			dNl = mlmc_ests1->newN - mlmc_ests1->N;
			if (dNl > 0) {
				mlmc_ests1->N += dNl;
				mlmc->level_l_estimation(mlmc, mlmc_ests1, l, M, dNl);
			}

			mlmc_ests1 = mlmc_ests1->next;
		}

		//test for convergence
		if (extraporation_flag) {
			if (L > 1) {
				con = mlmc_ests->sumY/mlmc_ests->N - (mlmc->T/(double)n) * mlmc_ests->prev->sumY/mlmc_ests->prev->N;
				if (fabs(con) < 1.0*(M*M - 1.0)*eps/sqrt(2.0)) {
					converged_flag = 1;
				}
			}
		} else {
			if (L > 1) {
				con = fmax(fabs(mlmc_ests->sumY/mlmc_ests->N), fabs(mlmc_ests->prev->sumY/mlmc_ests->prev->N/M));
				if (con < (M-1)*eps/sqrt(2.0)) {
					converged_flag = 1;
				}
			}
		}

		if (mlmc_ests->next == NULL) {
			add_level_l_estimator(mlmc_ests);
		}
		mlmc_ests = mlmc_ests->next;
	}

	sum = 0.0;
	mlmc_ests1 = mlmc->list;
	while (mlmc_ests1->next != NULL) {
		sum += mlmc_ests1->sumY / mlmc_ests1->N;
		mlmc_ests1 = mlmc_ests1->next;
	}
	if (extraporation_flag) {
		sum += mlmc_ests->prev->sumY/mlmc_ests->prev->N/(M-1);
	}
	mlmc->P = sum;

	return 0;
}

