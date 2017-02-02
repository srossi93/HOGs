 /***************************************************************************
 *   Copyright (C) 2009 by Gabriele Moser                                  *
 *                                                                         *
 *   This code has been written for the OPERA "Civil protection from       *
 *   flooding events" project, funded by the Italian Space Agency (ASI)    *
 *   by modifying the LIBSVM 2.89 package by C.-C. Chang and C.-J. Lin.    *
 *   Copyright information for LIBSVM 2.89 is below.                       *
 *                                                                         *
 *                                 *******                                 *
 *                                                                         *
 *					          MODIFIED LIBSVM                              *
 *                                                                         *
 ***************************************************************************
 *                                                                         *
 * Copyright (c) 2000-2009 Chih-Chung Chang and Chih-Jen Lin               *
 * All rights reserved.                                                    *
 *                                                                         *
 * Redistribution and use in source and binary forms, with or without      *
 * modification, are permitted provided that the following conditions      *
 * are met:                                                                *
 *                                                                         *
 * 1. Redistributions of source code must retain the above copyright       *
 * notice, this list of conditions and the following disclaimer.           *
 *                                                                         *
 * 2. Redistributions in binary form must reproduce the above copyright    *
 * notice, this list of conditions and the following disclaimer in the     *
 * documentation and/or other materials provided with the distribution.    *
 *                                                                         *
 * 3. Neither name of copyright holders nor the names of its contributors  *
 * may be used to endorse or promote products derived from this software   *
 * without specific prior written permission.                              *
 *                                                                         *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS     *
 * ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT     *
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR   *
 * A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE REGENTS OR  *
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,   *
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,     *
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR      *
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF  *
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING    *
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS      *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.            *
 *                                                                         *
 ***************************************************************************/


#ifndef _LIBSVM_H
#define _LIBSVM_H
#define LIBSVM_VERSION 289


enum { C_SVC, EPSILON_SVR};
enum { LINEAR, POLY, RBF, SIGMOID, MRF, GUSTAVO};
const unsigned short powelliter=10;
const unsigned short brentiter=100;
const unsigned short nUSVmin=2;	// era 0
const double param_min[3]={0.000001,0.00001,0.000001},param_max[3]={100000,10,0.1},param_init[3]={1,2,0.01},SigmaThr=DBL_MAX;
const short theta_min[3]={-10,-9,-10},theta_max[3]={10,5,-2};
const bool powell_svr_tanh=false;
enum svmoption { GIVEN,PSB};


struct svm_node				// vettore delle feature
{
	unsigned short dim;		// n. di feature
	double *values;			// v. delle f.
};


struct svm_problem			// training set
{
	unsigned long l;		// n. campioni di training
	double *y;				// puntatore per allocazione dinamica delle etichette di classe
	struct svm_node *x;		// e quello per i v. delle f.
};


struct svm_parameter		// configurazione dei parametri
{
	int svm_type;			// class? altro? per noi è C_SVC
	int kernel_type;		// gaussiano (rbf = radial basis function)? lineare? altro?
	double gamma;			// = 0.5/var nel kernel gaussiano
	double lambda;			// per svm mischiata con mrf
	double C;

	// parametri di configurazione del programmatore quadratico
	int shrinking;
	double cache_size;	// [MB]
	double eps;			// stop threshold

	// flag per ottimizzazione automatica dei parametri
	bool PSB;			// svm "pura"
	bool HK;			// svm mischiata con mrf

	// inutili per noi
	bool GGGH;
	int degree;			// for poly
	double p;			// epsilon for EPSILON_SVR
	double coef0;		// for poly, sigmoid
};


struct svm_model			// svm addestrata
{
	svm_parameter param;	// parameter
	int nr_class;			// number of classes, = 2 in regression/one class svm
	unsigned long l;		// total #SV
	svm_node *SV;			// SVs (SV[l])
	double **sv_coef;		// coefficients for SVs in decision functions (sv_coef[k-1][l]) ---> gli alpha
	double *rho;			// constants in decision functions (rho[k*(k-1)/2])
	int *label;				// label of each class (label[k])
	int *nSV;				// number of SVs for each class (nSV[k]), nSV[0] + nSV[1] + ... + nSV[k-1] = l
	int free_sv;			// 1 if svm_model is created by svm_load_model; 0 if svm_model is created by svm_train
	double SB;				// (average) span bound
	double noise_var;
	Matrix L;
	vector<longword> usv;
};


struct svm_model *svm_train(const struct svm_problem *prob, const struct svm_parameter *param);
double svm_cross_validation(const struct svm_problem *prob, const struct svm_parameter *param, int nr_fold);
void svm_predict_values(const struct svm_model *model, const struct svm_node *x, double* dec_values);
double svm_predict(const struct svm_model *model, const struct svm_node *x);
double svm_predict(const struct svm_model *model, const struct svm_node *x, double standard_deviation[3]);
void init_default_param(svm_parameter *param, int type);
void svm_destroy_prob(struct svm_problem *prob);
void svm_destroy_model(struct svm_model *model);
void Powell(struct svm_problem *prob, struct svm_parameter *param);
vector<bool> PDDP(Matrix X, longword target);
void svm_regression_read_ground(struct svm_problem *prob, Vector gain, Vector offset, char trainname[80]);
void svm_regression_write_test(struct svm_model *model, char testname[80], char estimatename[80]);
Vector svm_regression_holdout(struct svm_model *model, struct svm_problem *prob, double gain, double offset, char estimatename[80]);
Matrix PCA(struct svm_problem *prob);
void svm_regression_bias_training(struct svm_problem *prob, struct svm_model *model);

Vector QPparamMRF(Matrix energy_difference_transposed);



#endif /* _LIBSVM_H */
