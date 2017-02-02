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
#include "omp.h"
#include "image.h"
#include "library.h"
//#include "kdtree.h"
//#include "svm.h"
//#include "graphcut.h"
#include "dibe001_classifier.h"


using namespace std;

void Classifier::GaussBayes(Image *img, Image *train, Image *postprob, bayesoption option)
{
	word w = img->get_width(), h = img->get_height(), b = img->get_bands(), C = word(train->maxrange(0)), r, s, m, n, c;
	if ((w != train->get_width()) || (h != train->get_height()))
		errorlog("Classifier->GaussBayes. Inconsistent image sizes.");
	if (C < 2)
		errorlog("Classifier->GaussBayes. Inconsistent numbers of classes.");
	this->create(w, h, 1);
	postprob->create(w, h, C);
	Vector size = Vector(C, 0);
	Vector detcov = Vector(C, 1);
	Matrix mean = Matrix(C, b, 0);
	Tensor cov = Tensor(C, b, b, 0);
	Tensor invcov = Tensor(C, b, b, 0);
	Vector post = Vector(C, 0);
	for (m = 0; m < w; m++)
	for (n = 0; n < h; n++)
	if (c = word(train->get(m, n, 0)))
	{
		size[c - 1]++;
		for (r = 0; r < b; r++)
		{
			mean[c - 1][r] += img->get(m, n, r);
			for (s = r; s < b; s++) cov[c - 1][r][s] += img->get(m, n, r)*img->get(m, n, s);
		}
	}
	for (c = 0; c<C; c++)
	if (size[c]>0)
	{
		for (r = 0; r < b; r++)
		{
			mean[c][r] /= size[c];
			cov[c][r][r] = cov[c][r][r] / size[c] - mean[c][r] * mean[c][r];
		}
		for (r = 0; r < b; r++)
		for (s = r + 1; s < b; s++)
		{
			cov[c][r][s] = cov[c][r][s] / size[c] - mean[c][r] * mean[c][s];
			cov[c][s][r] = cov[c][r][s];
		}
		if (option == ML) size[c] = 1;
	}
	else
		errorlog("Classifier->GaussBayes. Empty class.");
	for (c = 0; c < C; c++)
	{
		detcov[c] = InvDetCholesky(cov[c], invcov[c]);
		if (detcov[c] <= 0)
			errorlog("Classifier->GaussBayes. Nonpositive-definite sample covariance.");
	}
	for (m = 0; m < w; m++)
	for (n = 0; n < h; n++)
	{
		word label = 0;
		double postmax = -1, posttot = 0;
		for (c = 0; c < C; c++)
		{
			double bhatta = 0;
			for (r = 0; r < b; r++)
			{
				bhatta += invcov[c][r][r] * (img->get(m, n, r) - mean[c][r])*(img->get(m, n, r) - mean[c][r]);
				for (s = r + 1; s<b; s++) bhatta += 2 * invcov[c][r][s] * (img->get(m, n, r) - mean[c][r])*(img->get(m, n, s) - mean[c][s]);
			}
			post[c] = size[c] * exp(-0.5*bhatta) / sqrt(detcov[c]);
			if (post[c]>postmax)
			{
				postmax = post[c];
				label = c + 1;
			}
			posttot += post[c];
		}
		this->put(m, n, 0, label);
		for (c = 0; c<C; c++)
		{
			if (posttot>0)
				post[c] /= posttot;
			else
				post[c] = 1 / double(C);
			postprob->put(m, n, c, post[c]);
		}
	}
}


/***************************************************************************************************************************************/
double MarkovPottsHK(Image *post, Image *train, Image *map, Vector size)
{
	word w = post->get_width(), h = post->get_height(), C = post->get_bands(), k, m, n, c;
	longword trainsize = 0, constraints = 0;
	for (c = 0; c < C; c++)
		trainsize += longword(size[c]);
	Matrix XTR = Matrix(2, trainsize*(C - 1), 0);
	for (m = 1; m < w - 1; m++)
	for (n = 1; n<h - 1; n++)
	if ((c = word(train->get(m, n, 0))) && (post->get(m, n, c - 1)>0))
	{
		Vector UX = Vector(C, 0);
		Vector US = Vector(C, 0);
		for (k = 0; k < C; k++)
		{
			UX[k] = -log(post->get(m, n, k)) + log(size[k]);
			US[k] = 0;
		}
		short u, v;
		for (u = -1; u<2; u++)
		for (v = -1; v<2; v++)
		if ((u*u + v*v>0) && (k = word(map->get(m + u, n + v, 0)))) US[k - 1]--;
		for (k = 0; k<C; k++)
			//if((k!=c-1)&&(post->get(m,n,k)>0)&&((UX[k]-UX[c-1])*(US[k]-US[c-1])<0))
		if ((k != c - 1) && (post->get(m, n, k)>0) && ((UX[k] - UX[c - 1]>0) || (US[k] - US[c - 1] > 0)))
		{
			XTR[0][constraints] = UX[k] - UX[c - 1];
			XTR[1][constraints] = US[k] - US[c - 1];
			constraints++;
		}
	}
	if (constraints == 0)
		errorlog("Classifier->MarkovPottsHK. No linear equations to be solved.");
	XTR[0].resize(constraints);
	XTR[1].resize(constraints);
	//Vector alpha=HoKashyap(XTR,hkiter);
	Vector alpha = QPparamMRF(XTR);
	double beta = 1;
	if ((alpha[0] > 0) && (alpha[1] > 0)) beta = alpha[1] / alpha[0];
	//if((beta!=beta)||(beta<=std::numeric_limits<double>::min())||(beta>=std::numeric_limits<double>::max()))
	//	beta=1;		// controllo per verificare che beta non sia NaN ne +-oo
	cout << endl << beta << endl;
	return beta;
}


void Classifier::MarkovPottsICM(Image *post, Image *train, Image *map, double &beta, mrfoption option)
{
	word w = post->get_width(), h = post->get_height(), C = post->get_bands(), m, n, c, iter;
	if ((map->get_bands() > 1) || (train->get_bands() > 1))
		errorlog("Classifier->MarkovPottsICM. Erroneously multichannel input map.");
	if ((w != train->get_width()) || (h != train->get_height()) || (w != map->get_width()) || (h != map->get_height()))
		errorlog("Classifier->MarkovPottsICM. Inconsistent image sizes.");
	if ((C < 2) || (C != word(train->maxrange(0))) || (C < map->maxrange(0)))
		errorlog("Classifier->MarkovPottsICM. Inconsistent numbers of classes.");
	Vector size = Vector(C, 0);
	this->create(w, h, 1);
	for (m = 0; m < w; m++)
	for (n = 0; n < h; n++)
	{
		this->put(m, n, 0, map->get(m, n, 0));
		if (c = word(train->get(m, n, 0))) size[c - 1]++;
	}
	for (c = 0; c < C; c++)
	if (size[c] == 0)
		errorlog("Classifier->MarkovPottsICM. Empty class.");
	if (option == HK)
		beta = MarkovPottsHK(post, train, map, size);
	Image *buffer = new Image;
	buffer->create(w, h, 1);
	for (iter = 0; iter < icmiter; iter++)
	{
		cout << iter + 1 << " ";
		flush(cout);
		for (m = 0; m < w; m++)
		for (n = 0; n < h; n++)
			buffer->put(m, n, 0, this->get(m, n, 0));
		for (m = 1; m < w - 1; m++)
		for (n = 1; n < h - 1; n++)
		{
			Vector US = Vector(C, 0);
			short u, v;
			for (u = -1; u < 2; u++)
			for (v = -1; v<2; v++)
			if ((u*u + v*v>0) && (c = word(buffer->get(m + u, n + v, 0)))) US[c - 1]--;
			word label = 0;
			double maxlocalpost = -DBL_MAX;
			for (c = 0; c<C; c++)
			{
				double localpost = post->get(m, n, c)*exp(-beta*US[c]) / size[c];
				if (localpost>maxlocalpost)
				{
					maxlocalpost = localpost;
					label = c + 1;
				}
			}
			this->put(m, n, 0, label);
		}
	}
	buffer->kill();
	delete buffer;
}


/***************************************************************************************************************************************/
word CrossValidation(record *training, longword numrec, word folds, word b, word C)
{
	const word Klow = 5, Kup = 50, Kstep = 5;
	const int bucket_size = 1;
	vector<word> partition = vector<word>(numrec, 0);
	word *counter = new word[C], f, k, K = 1, r, c;
	longword i, j;
	srand(unsigned int(time(NULL)));
	for (i = 0; i < numrec; i++)
		partition[i] = word(floor(folds*float(rand()) / RAND_MAX));
	unsigned int *neighbor;
	double minCV = DBL_MAX;
	for (k = Klow; k < Kup + 1; k += Kstep)
	{
		double CV = 0;
		for (f = 0; f < folds; f++)
		{
			longword subnumrec = 0, I = 0, suberrors = 0, subtest = 0;
			for (i = 0; i < numrec; i++)
			if (partition[i] != f)
				subnumrec++;
			if (subnumrec == 0)
				errorlog("CrossValidation. Empty training subset.");
			record *subtraining = new record[subnumrec];
			for (i = 0; i < numrec; i++)
			if (partition[i] != f)
			{
				subtraining[I].label = training[i].label;
				subtraining[I].data = new float[b];
				subtraining[I].index = I;
				for (r = 0; r < b; r++)
					subtraining[I].data[r] = training[i].data[r];
				I++;
			}
			kdtree *tree = new kdtree(subtraining, subnumrec, bucket_size, b);
			for (i = 0; i < numrec; i++)
			if (partition[i] == f)
			{
				subtest++;
				neighbor = tree->search(&training[i], k);
				for (c = 0; c < C; c++)
					counter[c] = 0;
				for (j = 1; j < k; j++)
				if (neighbor[j] != inf)
					counter[subtraining[neighbor[j]].label - 1]++;
				word maxcounter = 0, Label = 0;
				for (c = 0; c<C; c++)
				if (counter[c]>maxcounter)
				{
					maxcounter = counter[c];
					Label = c + 1;
				}
				if (Label != training[i].label) suberrors++;
			}
			if (subtest == 0)
				errorlog("CrossValidation. Empty test subset.");
			CV += double(suberrors) / (subtest*folds);
			tree->delete_pqr();
			delete tree;
		}
		cout << k << "	" << CV << endl;
		flush(cout);
		if (CV < minCV)
		{
			minCV = CV;
			K = k;
		}
	}
	delete[] counter;
	return K;
}


word Classifier::KNN(Image *img, Image *train, Image *postprob, word fold_neighbor, knnoption option)
{
	word w = img->get_width(), h = img->get_height(), b = img->get_bands(), C = word(train->maxrange(0)), m, n, r, c, K;
	const double eps = 0.01;
	if ((w != train->get_width()) || (h != train->get_height()))
		errorlog("Classifier->KNN. Inconsistent image sizes.");
	if (C < 2)
		errorlog("Classifier->KNN. Inconsistent numbers of classes.");
	this->create(w, h, 1);
	postprob->create(w, h, C);
	longword numrec = 0, i = 0;
	for (m = 0; m < w; m++)
	for (n = 0; n<h; n++)
	if (train->get(m, n, 0)>0) numrec++;
	if (numrec == 0)
		errorlog("Classifier->KNN. Empty training set.");
	record *training = new record[numrec], *query = new record;
	query->label = 0;
	query->index = 0;
	query->data = new float[b];
	for (m = 0; m < w; m++)
	for (n = 0; n < h; n++)
	if (c = word(train->get(m, n, 0)))
	{
		training[i].label = c;
		training[i].data = new float[b];
		training[i].index = i;
		for (r = 0; r < b; r++)
			training[i].data[r] = float(img->get(m, n, r));
		i++;
	}

	if (option == INPUT)
		K = fold_neighbor;
	else
		K = CrossValidation(training, numrec, fold_neighbor, b, C);

	const int bucket_size = 1;
	kdtree *tree = new kdtree(training, numrec, bucket_size, b);
	unsigned int *neighbor;
	vector<word> counter = vector<word>(C, 0);
	for (m = 0; m < w; m++)
	{
		cout << m + 1 << " "; flush(cout);
		for (n = 0; n < h; n++)
		{
			for (r = 0; r < b; r++)
				query->data[r] = float(img->get(m, n, r));
			neighbor = tree->search(query, K);
			for (c = 0; c < C; c++)
				counter[c] = 0;
			for (i = 0; i < K; i++)
			if (neighbor[i] != inf)
				counter[training[neighbor[i]].label - 1]++;
			word maxcounter = 0, Label = 0;
			for (c = 0; c<C; c++)
			{
				postprob->put(m, n, c, double(counter[c]) / K);
				if (counter[c]>maxcounter)
				{
					maxcounter = counter[c];
					Label = c + 1;
				}
			}
			this->put(m, n, 0, Label);
		}
	}
	tree->delete_pqr();
	//for(i=0;i<numrec;i++) delete[] training[i].data;
	delete[] training, query->data;
	delete query;
	return (option == INPUT ? 0 : K);
}


/***************************************************************************************************************************************/
vector<longword> TargetSizes(Image *train, vector<bool> &subsample, longword &sizetot, longword svm_subsample_threshold)
{
	word w = train->get_width(), h = train->get_height(), C = subsample.size(), m, n, c;
	longword balancedsizetot = 0, sizemin = w*h + 1;
	vector<longword> size = vector<longword>(C, 0);
	for (m = 0; m < w; m++)
	for (n = 0; n < h; n++)
	if (c = word(train->get(m, n, 0)))
		size[c - 1]++;
	for (c = 0; c < C; c++)
	if (size[c] == 0)
		errorlog("BalanceClasses. Empty class.");
	else if (size[c] < sizemin)
		sizemin = size[c];
	for (c = 0; c<C; c++)
	{
		if (size[c]>svm_sizeratio_threshold*sizemin)
		{
			size[c] = svm_sizeratio_threshold*sizemin;
			subsample[c] = true;
		}
		else
			subsample[c] = false;
		balancedsizetot += size[c];
	}
	if (balancedsizetot > svm_subsample_threshold)
	for (c = 0; c < C; c++)
	{
		longword sizecandidate = longword(size[c] * float(svm_subsample_threshold) / balancedsizetot);
		size[c] = (sizecandidate < svm_min_class_size_threshold ? svm_min_class_size_threshold : sizecandidate);
		subsample[c] = true;
	}
	sizetot = 0;
	for (c = 0; c < C; c++)
		sizetot += size[c];
	return size;
}


svm_problem *InitProblemPDDP(Image *img, Image *train, word C, longword svm_subsample_threshold, bool overwrite_train)
{
	word w = img->get_width(), h = img->get_height(), b = img->get_bands(), m, n, r, c;
	longword i = 0, sizetot = 0;
	vector<bool> subsample = vector<bool>(C, false);
	vector<longword> target = TargetSizes(train, subsample, sizetot, svm_subsample_threshold);
	struct svm_problem *prob = new svm_problem;
	prob->l = sizetot;
	prob->y = new double[prob->l];
	prob->x = new svm_node[prob->l];
	for (c = 0; c < C; c++)
	{
		if (subsample[c])
		{
			cout << "PDDP over class " << c + 1 << endl;
			flush(cout);
			Matrix trainsamples;
			for (m = 0; m < w; m++)
			for (n = 0; n < h; n++)
			if (word(train->get(m, n, 0)) == c + 1)
			{
				Vector x = Vector(b);
				for (r = 0; r < b; r++)
					x[r] = img->get(m, n, r);
				trainsamples.push_back(x);
			}
			vector<bool> selected = PDDP(trainsamples, target[c]);
			longword j = 0;
			for (m = 0; m < w; m++)
			for (n = 0; n < h; n++)
			if (word(train->get(m, n, 0)) == c + 1)
			{
				if (selected[j])
				{
					prob->x[i].dim = b;
					prob->x[i].values = new double[b];
					for (r = 0; r < b; r++)
						prob->x[i].values[r] = img->get(m, n, r);
					prob->y[i] = c + 1;
					i++;
				}
				else if (overwrite_train)
					train->put(m, n, 0, 0);
				j++;
			}
		}
		else
		{
			for (m = 0; m < w; m++)
			for (n = 0; n < h; n++)
			if (word(train->get(m, n, 0)) == c + 1)
			{
				prob->x[i].dim = b;
				prob->x[i].values = new double[b];
				for (r = 0; r < b; r++)
					prob->x[i].values[r] = img->get(m, n, r);
				prob->y[i] = c + 1;
				i++;
			}
		}
	}
	return prob;
}


svm_problem *InitSVC(Image *img, Image *train, svm_parameter *param, word C)
{
	init_default_param(param, C_SVC);
	struct svm_problem *prob = InitProblemPDDP(img, train, C, (param->PSB ? svm_small_subsample_threshold : svm_large_subsample_threshold), !param->PSB);
	if (param->PSB)
	{
		Powell(prob, param);
		svm_destroy_prob(prob);
		prob = InitProblemPDDP(img, train, C, svm_large_subsample_threshold, true);
		train->store("TrainMapSubsample.raw", BYTE);
		param->PSB = false;
	}
	return prob;
}


void Classifier::SVC(Image *img, Image *train, svm_parameter *param)
{
	word w = img->get_width(), h = img->get_height(), b = img->get_bands(), C = word(train->maxrange(0)), m, n, r;
	if ((w != train->get_width()) || (h != train->get_height()) || (train->get_bands() != 1))
		errorlog("Classifier->SVC. Inconsistent image sizes.");
	if (C < 2)
		errorlog("Classifier->SVC. Inconsistent number of classes.");
	this->create(w, h, 1, 0);
	struct svm_problem *prob = InitSVC(img, train, param, C);
	struct svm_model *model = svm_train(prob, param);
	if (model->nSV == 0)
		errorlog("Classifier->SVC. No support vectors.");
	svm_node *x = new svm_node;
	x->dim = b;
	x->values = new double[b];
	cout << "column ";

	//Image *test = new Image;
	//test->load("TestMap.raw",w,h,1,BYTE);
	

	omp_set_num_threads(3);
	//double t = (double)cv::getTickCount();
	int IterPerCycle = 500;
	int outer = 0, totalNo = ceil((float)w / IterPerCycle);
	for (; outer<totalNo; outer++) {
		long loc_max = (outer + 1)*IterPerCycle;
		if (loc_max>w) loc_max = w;
		cout << "Processing columns " << outer*IterPerCycle + 1 << " thru " << loc_max << ";\n";
#pragma omp parallel for
		for (int m = outer*IterPerCycle; m<loc_max; m++) {
			svm_node *x = new svm_node;
			x->dim = b;
			x->values = new double[b];
			for (int n = 0; n<h; n++){
				//if (mask->get(m, n, 0) == 0) continue;
				for (int r = 0; r<b; r++)
					x->values[r] = img->get(m, n, r);
				//if(x->values[0]>0) this->put(m,n,0,svm_predict(model,x)); else this->put(m,n,0,0);
				this->put(m, n, 0, svm_predict(model, x));
			}
			delete[] x->values;
			delete x;
		}
	}
	//t = ((double)cv::getTickCount() - t) / cv::getTickFrequency();
	//std::cout << "\nTime elapsed is " << t << " sec;\n";

	
	////////////NON PARALLEL SVM
	//for (m = 0; m < w; m++)
	//{
	//	cout << m + 1 << " ";
	//	flush(cout);
	//	//for(n=0;n<h;n++) if((test->get(m,n,0)>0)||param->HK)
	//	for (n = 0; n < h; n++)
	//	{
	//		for (r = 0; r < b; r++)
	//			x->values[r] = img->get(m, n, r);
	//		//if(x->values[0]>0) this->put(m,n,0,svm_predict(model,x)); else this->put(m,n,0,0);
	//		this->put(m, n, 0, svm_predict(model, x));
	//	}
	//}
	svm_destroy_prob(prob);
	svm_destroy_model(model);
	delete[] x->values;
	delete x, prob, model;
}


/***************************************************************************************************************************************/
double EpsilonCompute(Image *map, word m, word n)
{
	short u, v, epsilon = 0, counter = 0, w = map->get_width(), h = map->get_height();
	for (u = -1; u < 2; u++)
	for (v = -1; v<2; v++)
	if ((u*u + v*v>0) && (m + u >= 0) && (m + u < w) && (n + v >= 0) && (n + v < h))
	if (word(map->get(m + u, n + v, 0)) == 1)
	{
		epsilon++;
		counter++;
	}
	else if (word(map->get(m + u, n + v, 0)) == 2)
	{
		epsilon--;
		counter++;
	}
	else if (word(map->get(m + u, n + v, 0)) > 2)
		errorlog("EpsilonCompute. Unacceptable binary class label.");
	return (counter > 0 ? float(epsilon) / counter : 0);
}


double EpsilonCompute(Image *img, svm_model *model, word m, word n)
{
	short u, v, epsilon = 0, counter = 0, w = img->get_width(), h = img->get_height(), b = img->get_bands(), r;
	svm_node *x = new svm_node;
	x->dim = b + 1;
	x->values = new double[b + 1];
	x->values[b] = m + n*w;
	for (u = -1; u < 2; u++)
	for (v = -1; v<2; v++)
	if ((u*u + v*v>0) && (m + u >= 0) && (m + u < w) && (n + v >= 0) && (n + v < h))
	{
		for (r = 0; r < b; r++)
			x->values[r] = img->get(m + u, n + v, r);
		word label = word(svm_predict(model, x));
		if (label == 1)
		{
			epsilon++;
			counter++;
		}
		else if (label == 2)
		{
			epsilon--;
			counter++;
		}
	}
	delete[] x->values;
	delete x;
	return (counter > 0 ? float(epsilon) / counter : 0);
}


void MarkovSVCHK(Image *img, Image *train, svm_parameter *param, word C)
{
	word w = img->get_width(), h = img->get_height(), b = img->get_bands(), m, n, r, c, d;
	longword i;
	Matrix X;
	for (c = 0; c < C - 1; c++)
	for (d = c + 1; d < C; d++)
	{
		cout << endl << "pair (" << c + 1 << ", " << d + 1 << ")";
		svm_problem *prob = new svm_problem;
		prob->l = i = 0;
		for (m = 0; m < w; m++)
		for (n = 0; n < h; n++)
		if ((word(train->get(m, n, 0)) == c + 1) || (word(train->get(m, n, 0)) == d + 1))
			prob->l++;
		prob->y = new double[prob->l];
		prob->x = new svm_node[prob->l];
		for (m = 0; m < w; m++)
		for (n = 0; n < h; n++)
		if ((word(train->get(m, n, 0)) == c + 1) || (word(train->get(m, n, 0)) == d + 1))
		{
			prob->x[i].dim = b + 1;
			prob->x[i].values = new double[b + 1];
			for (r = 0; r < b; r++)
				prob->x[i].values[r] = img->get(m, n, r);
			prob->x[i].values[b] = m + n*w;
			prob->y[i] = ((word(train->get(m, n, 0)) == c + 1) ? 1 : 2);
			i++;
		}
		param->lambda = 0;	// giunto qui, ho gia' param->kernel_type = MRF
		struct svm_model *model = svm_train(prob, param);
		double sumepsilon = 0;
		for (i = 0; i<model->l; i++)
			sumepsilon += model->sv_coef[0][i] * EpsilonCompute(img, model, word(model->SV[i].values[b]) % w, word(model->SV[i].values[b]) / w);
		for (i = 0; i<prob->l; i++)
		{
			Vector column = Vector(2, model->rho[0]);
			double *dec_values = new double;
			svm_predict_values(model, &prob->x[i], dec_values);
			column[0] = dec_values[0];
			delete dec_values;
			column[1] = EpsilonCompute(img, model, word(prob->x[i].values[b]) % w, word(prob->x[i].values[b]) / w)*sumepsilon;
			if ((column[0]>0) || (column[1]>0))
				X.push_back(column);
		}
		svm_destroy_prob(prob);
		svm_destroy_model(model);
		delete prob, model;
	}
	Matrix XTR = Matrix(X.cols(), X.rows());
	for (i = 0; i<X.rows(); i++)
	for (r = 0; r<2; r++)
		XTR[r][i] = X[i][r];
	Vector Lambda = HoKashyap(XTR, hkiter);
	param->lambda = (Lambda[0]>0 && Lambda[1]>0 ? Lambda[1] / Lambda[0] : 0.1);
	param->HK = false; // giunto qui, ho ottimizzato automaticamente tutto quanto poteva essere ottimizzato automaticamente
}


void Classifier::MarkovSVC(Image *img, Image *train, svm_parameter *param)
{
	word w = img->get_width(), h = img->get_height(), b = img->get_bands(), C = word(train->maxrange(0)), m, n, r, c, d, iter, maxiter = 5;
	if ((w != train->get_width()) || (h != train->get_height()) || (train->get_bands() != 1))
		errorlog("Classifier->MarkovSVC. Inconsistent image sizes.");
	if (C < 2)
		errorlog("Classifier->MarkovSVC. Inconsistent number of classes.");
	InitSVC(img, train, param, C);	// qui dentro: sottocampionamento train + eventuale Powell
	Image *extimg = new Image;
	extimg->create(w, h, b + 1);
	for (m = 0; m < w; m++)
	for (n = 0; n < h; n++)
	for (r = 0; r < b; r++)
		extimg->put(m, n, r, img->get(m, n, r));
	param->kernel_type = MRF;
	if (param->HK)
		MarkovSVCHK(img, train, param, C);
	Image *votemap = new Image;
	votemap->create(w, h, C, 0);
	for (c = 0; c < C - 1; c++)
	for (d = c + 1; d < C; d++)
	{
		Image *bintrain = new Image;
		bintrain->create(w, h, 1, 0);
		for (m = 0; m < w; m++)
		for (n = 0; n < h; n++)
		{
			extimg->put(m, n, b, 0);
			if (word(train->get(m, n, 0)) == c + 1)
				bintrain->put(m, n, 0, 1);
			else if (word(train->get(m, n, 0)) == d + 1)
				bintrain->put(m, n, 0, 2);
		}
		for (iter = 0; iter < maxiter; iter++)
		{
			cout << endl << "ICM iteration " << iter + 1 << "; pair (" << c + 1 << ", " << d + 1 << ")";
			flush(cout);
			Classifier *binmap = new Classifier;
			binmap->SVC(extimg, bintrain, param);
			for (m = 0; m < w; m++)
			for (n = 0; n < h; n++)
			{
				extimg->put(m, n, b, EpsilonCompute(binmap, m, n));
				if ((iter == maxiter - 1) && (word(binmap->get(m, n, 0)) == 1))
					votemap->put(m, n, c, votemap->get(m, n, c) + 1);
				else if ((iter == maxiter - 1) && (word(binmap->get(m, n, 0)) == 2))
					votemap->put(m, n, d, votemap->get(m, n, d) + 1);
			}
			binmap->kill();
			delete binmap;
		}
		bintrain->kill();
		delete bintrain;
	}
	this->create(w, h, 1, 0);
	for (m = 0; m < w; m++)
	for (n = 0; n < h; n++)
	{
		word maxvote = 0, label = 0;
		for (c = 0; c<C; c++)
		{
			word vote = word(votemap->get(m, n, c));
			if (vote>maxvote)
			{
				maxvote = vote;
				label = c + 1;
			}
		}
		this->put(m, n, 0, label);
	}
	votemap->kill();
	extimg->kill();
	delete extimg, votemap;
}


/***************************************************************************************************************************************/
void Classifier::GraphCutOptimizePotts(Image *post, Vector size, float beta, GraphCutOption option, NeighborhoodOrder order)
{
	word w = post->get_width(), h = post->get_height(), C = post->get_bands(), m, n, c;
	float *UX = new float[w*h*C], tol = float(1e-6);
	for (m = 0; m < w; m++)
	for (n = 0; n < h; n++)
	for (c = 0; c < C; c++)
		UX[(m + n*w)*C + c] = float(log(size[c]) - log(post->get(m, n, c) + tol));
	post->kill();
	DataCost *data = new DataCost(UX);
	SmoothnessCost *smooth = new SmoothnessCost(1, 1, beta);
	EnergyFunction *eng = new EnergyFunction(data, smooth, order);
	MRFmodel *mrf;
	switch (option)
	{
	case GCEXP:
		mrf = new Expansion(w, h, C, eng);
		break;
	case GCSWAP:
		mrf = new Swap(w, h, C, eng);
		break;
	case MPBP:
		mrf = new MaxProdBP(w, h, C, eng);
		break;
	}
	mrf->initialize();
	mrf->clearAnswer();
	mrf->optimize(graph_cut_iterations);
	delete[] UX;
	this->create(w, h, 1);
	for (m = 0; m < w; m++)
	for (n = 0; n < h; n++)
		this->put(m, n, 0, word(mrf->getLabel(m + n*w) + 1));
	delete mrf, data, smooth, eng;
}


/*void ICMstepMultiRes(Image *pan, Image *ms, Classifier *map, Image *energymap, word ratio, Vector panmean, Vector var, Matrix msmean, Tensor cov)
{
word w=pan->get_width(),h=pan->get_height(),b=ms->get_bands(),C=panmean.size(),i,j,m,n,r,s,c,l;
for(i=0;i<w/ratio;i++)
{
cout << i+1 << " "; flush(cout);
for(j=0;j<h/ratio;j++)
{
double totlocalsize=0;
labelconfig baseconfig;
baseconfig.localsize = Vector(C,0);
baseconfig.localmean = Vector(b,0);
baseconfig.localcov = Matrix(b,b,0);
baseconfig.invlocalcov = Matrix(b,b);
for(m=i*ratio;m<(i+1)*ratio;m++)
for(n=j*ratio;n<(j+1)*ratio;n++)
if(c=word(map->get(m,n,0)))
{
baseconfig.localsize[c-1]++;
totlocalsize++;
}
if(totlocalsize>0)
for(c=0;c<C;c++)
for(r=0;r<b;r++)
{
baseconfig.localmean[r]+=baseconfig.localsize[c]*msmean[c][r]/totlocalsize;
baseconfig.localcov[r][r]+=baseconfig.localsize[c]*cov[c][r][r]/(totlocalsize*totlocalsize);
for(s=r+1;s<b;s++)
{
baseconfig.localcov[r][s]+=baseconfig.localsize[c]*cov[c][r][s]/(totlocalsize*totlocalsize);
baseconfig.localcov[s][r]=baseconfig.localcov[r][s];
}
}
else
errorlog("ICMstepMultiRes. Empty square.");
baseconfig.localdet=InvDetCholesky(baseconfig.localcov,baseconfig.invlocalcov);
vector<labelconfig> listconfig = vector<labelconfig>(1,baseconfig);
for(m=(i*ratio>1 ? i*ratio : 1);m<((i+1)*ratio<w-1 ? (i+1)*ratio : w-1);m++)
for(n=(j*ratio>1 ? j*ratio : 1);n<((j+1)*ratio<h-1 ? (j+1)*ratio : h-1);n++)
{
word current=word(map->get(m,n,0)),label=0;
for(c=0;c<C;c++)
if(var[c]>0)
{
labelconfig currentconfig;
currentconfig.localsize = Vector(baseconfig.localsize);
currentconfig.localsize[c]++;
currentconfig.localsize[current-1]--;
bool go_on=true;
for(l=0;(l<listconfig.size())&&go_on;l++)
if(listconfig[l].localsize==currentconfig.localsize)
{
currentconfig=listconfig[l];
go_on=false;
}
if(go_on)
{
currentconfig.localmean = Vector(baseconfig.localmean);
currentconfig.localcov = Matrix(baseconfig.localcov);
currentconfig.invlocalcov = Matrix(b,b);
for(r=0;r<b;r++)
{
currentconfig.localmean[r]+=(msmean[c][r]-msmean[current-1][r])/totlocalsize;
currentconfig.localcov[r][r]+=(cov[c][r][r]-cov[current-1][r][r])/(totlocalsize*totlocalsize);
for(s=r+1;s<b;s++)
{
currentconfig.localcov[r][s]+=(cov[c][r][s]-cov[current-1][r][s])/(totlocalsize*totlocalsize);
currentconfig.localcov[s][r]=currentconfig.localcov[r][s];
}
}
currentconfig.localdet=InvDetCholesky(currentconfig.localcov,currentconfig.invlocalcov);
listconfig.push_back(currentconfig);
}
double UPan=0.5*(pan->get(m,n,0)-panmean[c])*(pan->get(m,n,0)-panmean[c])/var[c]+0.5*log(var[c]),UMs=0.5*log(currentconfig.localdet);
for(r=0;r<b;r++)
{
UMs+=0.5*currentconfig.invlocalcov[r][r]*(ms->get(i,j,r)-currentconfig.localmean[r])*(ms->get(i,j,r)-currentconfig.localmean[r]);
for(s=r+1;s<b;s++)
UMs+=currentconfig.invlocalcov[r][s]*(ms->get(i,j,r)-currentconfig.localmean[r])*(ms->get(i,j,s)-currentconfig.localmean[s]);
}
energymap->put(m,n,c,UPan);
energymap->put(m,n,c+C,UMs);
}
else
errorlog("ICMstepMultiRes. Unacceptable zero variance.");
}
}
}
}


Vector HKloopMultiRes(Image *energymap, Image *train, Image *map, Vector size)
{
word w=energymap->get_width(), h=energymap->get_height(), C=size.size(),k,m,n,c;
longword trainsize=0,constraints=0;
for(c=0;c<C;c++)
trainsize+=longword(size[c]);
Matrix XTR = Matrix(3,trainsize*(C-1),0);
for(m=1;m<w-1;m++)
for(n=1;n<h-1;n++)
if(c=word(train->get(m,n,0)))
{
Vector UPan = Vector(C), UMs = Vector(C), US = Vector(C,0);
for(k=0;k<C;k++)
{
UPan[k]=energymap->get(m,n,k);
UMs[k]=energymap->get(m,n,k+C);
//US[k]=-(word(map->get(m-1,n,0))==k+1 ? 1 : 0)-(word(map->get(m+1,n,0))==k+1 ? 1 : 0)-(word(map->get(m,n-1,0))==k+1 ? 1 : 0)-(word(map->get(m,n+1,0))==k+1 ? 1 : 0);
}
short u,v;
for(u=-1;u<2;u++)
for(v=-1;v<2;v++)
if((u*u+v*v>0)&&(k=word(map->get(m+u,n+v,0)))) US[k-1]--;
for(k=0;k<C;k++)
if((k!=c-1)&&((UPan[k]-UPan[c-1]>0)||(UMs[k]-UMs[c-1]>0)||(US[k]-US[c-1]>0)))
{
XTR[0][constraints]=UPan[k]-UPan[c-1];
XTR[1][constraints]=UMs[k]-UMs[c-1];
XTR[2][constraints]=US[k]-US[c-1];
constraints++;
}
}
if(constraints==0)
errorlog("Classifier->MarkovPottsHK. No linear equations to be solved.");
XTR[0].resize(constraints);
XTR[1].resize(constraints);
XTR[2].resize(constraints);
Vector beta=HoKashyap(XTR,hkiter), alpha = Vector(2,1);
for(c=0;c<2;c++)
if((beta[c]>0)&&(beta[2]>0)) alpha[c]=beta[c]/beta[2];
return alpha;
}*/

