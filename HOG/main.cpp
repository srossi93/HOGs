/***************************************************************************
 *   Copyright (C) 2008 by Gabriele Moser                                  *
 *                                                                         *
 *   This code has been written for the OPERA "Civil protection from       *
 *   flooding events" project, funded by the Italian Space Agency (ASI).   *
 *                                                                         *
 ***************************************************************************/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif




#include "image.h"
#include "library.h"
using namespace std;
#include "dibe001_classifier.h"

#include "hog.h"


#define bins 9
#define cellSize 8
#define blockSize 2

#if !defined(ARRAY_SIZE)
#define ARRAY_SIZE(x) (int)(sizeof((x)) / sizeof((x)[0]))
#endif


int main(void)
{
	//Caricamento in memoria l'immagine da classificare (I) e l'immagine di training (TR)
	Image *I = new Image, *TR = new Image;

	I->load("Amiens_2006_SPOT_5m.raw", 2000, 2200, 3, BYTE);
	TR->load("GT_Amiens2006_5m_10classes_TR.raw", 2000, 2200, 1, BYTE);
	TR->stdmaker();

	/***************************************************************/

	Image *HOG0 = new Image, *HOG1 = new Image, *HOG2 = new Image;
	HOG0->create(I->get_width(),I->get_height(),bins);
	HOG1->create(I->get_width(), I->get_height(), bins);
	HOG2->create(I->get_width(), I->get_height(), bins);
	

	HOG(I, cellSize, blockSize, bins, 0, HOG0);
	
	HOG(I, cellSize, blockSize, bins, 1, HOG1);
	HOG(I, cellSize, blockSize, bins, 2, HOG2);
	HOG0->store("HOG0.raw", BYTE);
	HOG1->store("HOG1.raw", BYTE);
	HOG2->store("HOG2.raw", BYTE);
	
	cout << "HOGs created" << endl << flush;

	Image *finalHog = new Image;
	finalHog = MergeImg(I, HOG0, HOG1, HOG2);
	finalHog->store("FinalHog.raw", BYTE);

  finalHog->rescale(-1, 1);

	svm_parameter *param = new svm_parameter;
	param->kernel_type = RBF;
	
	Classifier *map = new Classifier;
	map->SVC(finalHog, TR, param);
	map->store("ClassMap.raw", BYTE);

	finalHog->kill();
	TR->kill();
	map->kill();
	delete TR, finalHog, map, param;

	cout << "Finito" << endl;

	getchar();
	exit(0);



	return 0;
}
