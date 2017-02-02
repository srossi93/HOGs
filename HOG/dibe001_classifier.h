/***************************************************************************
 *   Copyright (C) 2008 by Gabriele Moser                                  *
 *                                                                         *
 *   This code has been written for the OPERA "Civil protection from       *
 *   flooding events" project, funded by the Italian Space Agency (ASI).   *
 *                                                                         *
 *                                 *******                                 *
 *                                                                         *
 *                        CLASSIFIER CLASS (DIBE001)                       *
 *                                                                         *
 ***************************************************************************/

#include "kdtree.h"
#include "svm.h"
#include "graphcut.h"

enum bayesoption { ML,MAP};
enum knnoption { INPUT,CV};
enum multiresoption { SOLBERG,GEEXP,GESWAP,GELBP};
const word svm_sizeratio_threshold=10;
const word svm_small_subsample_threshold=3000;
const word svm_large_subsample_threshold=10000;
const word svm_min_class_size_threshold=10;
const word graph_cut_iterations=100;
const bool graph_cut_first_order=false;
const bool lowmem=true;


svm_problem *InitProblemPDDP(Image *img, Image *train, word C, longword svm_subsample_threshold, bool overwrite_train);

class Classifier: public Image
{
public:

	void GaussBayes(Image *image, Image *training_map, Image *posterior, bayesoption option);						// MAP o ML gaussiano
	word KNN(Image *image, Image *training_map, Image *posterior, word K_neighbor_or_folds, knnoption option);		// k punti vicini
	void MarkovPottsICM(Image *posterior, Image *training_map, Image *init_map, double &beta, mrfoption option);	// classificatore su MRF
	void SVC(Image *image, Image *training_map, svm_parameter *param);												// SVM per classificazione
	void MarkovSVC(Image *image, Image *training_map, svm_parameter *param);										// combinazione SVM-MRF
	void GraphCutOptimizePotts(Image *pixelwise_post, Vector training_set_size, float beta, GraphCutOption option, NeighborhoodOrder order);	// sempre MRF
	// ma con energia minimizzata mediante graph cut


protected:

};


