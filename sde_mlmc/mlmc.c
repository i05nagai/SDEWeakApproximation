/**
 * @file	mlmc.c
 * @brief	some methods required MLMC algorithms & level l estimators.
 *
 * 			There are some utility functions.
 * @date	09/Mar/2015 (Mon) 17:32 
 * @author	i05nagai
 */
#include <stdlib.h>
#include <math.h>
#include <sde_mlmc/mlmc.h>

/**
 * transform uniform random number to normal random number.
 * transform uniform random number to normal random number.
 * @param unif_rand n dimension uniform random number.
 * @param normal_rand n dimension normal random number.
 * @param n dimension of randm number. n must be even number.
 * @return none.
 */
void trans_rand_unif_to_normal(double *unif_rand, double *normal_rand, int n) {
	int i;
	double *u_seq1, *u_seq2;
	double *n_seq1, *n_seq2;
	int mid = n/2;

	if (n%2 == 0) {
		u_seq1 = unif_rand;
		u_seq2 = unif_rand + mid;
		n_seq1 = normal_rand;
		n_seq2 = normal_rand + mid;
		for (i=0; i<mid; i++) {
			n_seq1[i] = sqrt(-2.0*log(u_seq1[i]))*cos(2.0*M_PI*u_seq2[i]);
			n_seq2[i] = sqrt(-2.0*log(u_seq1[i]))*sin(2.0*M_PI*u_seq2[i]);
		}
	} else {
		u_seq1 = unif_rand;
		u_seq2 = unif_rand + mid + n%2;
		n_seq1 = normal_rand;
		n_seq2 = normal_rand + mid;
		for (i=0; i<mid; i++) {
			n_seq1[i] = sqrt(-2.0*log(u_seq1[i]))*cos(2.0*M_PI*u_seq2[i]);
			n_seq2[i] = sqrt(-2.0*log(u_seq1[i]))*sin(2.0*M_PI*u_seq2[i]);
		}
		n_seq2[i] = sqrt(-2.0*log(u_seq1[i]))*sin(2.0*M_PI*u_seq2[i]);
	}
}


/**
 * level l estimator allocation.
 * level l estimator allocation.
 * @return allocated level l estimator.
 */
SDE_MLMC *alloc_SDE_MLMC(SDE_WA_SLTN *sl) {
	SDE_MLMC *mlmc;

	mlmc = (SDE_MLMC *)malloc(sizeof(SDE_MLMC));
	mlmc->sl = sl;
	mlmc->list = alloc_MLMC_L();
	mlmc->list->l = 0;
	mlmc->init_y = (double *)malloc(sizeof(double)*sl->sde->dim_y);
	mlmc->P = 0.0;

	return mlmc;
}

/**
 * freeing the SDE_MLMC.
 * The function frees the pointer of SDE_MLMC.
 * You must free SDE_WA_SLTN in the pointer by calling free_SDE_WA_SLTN.
 * @param mlmc_l This parameter will be freed.
 */
void free_SDE_MLMC(SDE_MLMC *mlmc){
	free_MLMC_L(mlmc->list);
	free(mlmc->init_y);
	free(mlmc);
}

/**
 * level l estimator allocation.
 * level l estimator allocation.
 * @return allocated level l estimator.
 */
MLMC_L *alloc_MLMC_L(void) {
	MLMC_L *mlmc_l;

	mlmc_l = (MLMC_L *)malloc(sizeof(MLMC_L));
	mlmc_l->sumY = 0.0;
	mlmc_l->sumYY = 0.0;
	mlmc_l->N = 10000;
	mlmc_l->newN = 0;

	mlmc_l->prev = NULL;
	mlmc_l->next = NULL;

	return mlmc_l;
}

/**
 * freeing the list of level l estimators.
 * The function frees the list of level l estimators correctly.
 * @param mlmc_l This parameter will be freed.
 */
void free_MLMC_L(MLMC_L *mlmc_l){
	MLMC_L *temp = mlmc_l->next;
	while (mlmc_l->next != NULL) {
		free(mlmc_l);
		mlmc_l = temp;
		temp = temp->next;
	}
	free(mlmc_l);
}


/**
 * initializing level l estimator.
 * initializing level l estimator.
 * @param mlmc_l This parameter will be initilized.
 */
void init_level_l_estimator(MLMC_L *mlmc_l) {
    mlmc_l->sumY = 0.0;
    mlmc_l->sumYY = 0.0;
    mlmc_l->V = 0.0;
    mlmc_l->N = 10000;
    mlmc_l->newN = 0;
}              

/**
 * adding a new estimator at the end of the list.
 * adding a new estimator at the end of the list.
 * @param mlmc_l the list added a new new estimator.
 */
void add_level_l_estimator(MLMC_L *mlmc_l) {
	if (mlmc_l == NULL) {
		return;
	}

	MLMC_L *temp = mlmc_l;
	while (temp->next != NULL) {
		temp = temp->next;
	}
	temp->next = alloc_MLMC_L();
	temp->next->prev = temp;
	temp->next->l = temp->l + 1;
}

