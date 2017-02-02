 /***************************************************************************
 *																           *
 *                                                                         *
 *                                                                         *
 *                LIBRARY OF SUPPORT FOR REALIZING HOG                     *
 *                                                                         *
 *                                                                         *
 *                                                                         *
 *                                                                         *
 *                                                                         *
 ***************************************************************************/

////#include <cstdio>
////#include <iostream>
////#include <fstream>
////#include <cmath>
////#include <cfloat>
////#include <cstdlib>
////#include <ctime>
////#include <vector>
////#include <limits>

#ifndef HOG_library
#define HOG_library

#include "image.h"
using namespace std;

extern Image* MergeImg(Image *I, Image *hog0, Image *hog1, Image *hog2);



struct Point
{
	int x;
	int y;
};


extern Image *Gradient(Image *img, int band);
extern Image *Transpose(Image *img);
extern void MagnitudePh(Image *imgX, Image *imgY, Image *magnI, Image *phI);
extern Image * CreateCell(Image * img, int size, int row, int column);
extern double * CreateIstogram(Image *cellPh, Image *cellMagn, int bins);
extern double Norm2Calculator(double arr[], int dim);
extern void NormalizeBlock(double ***cellMatrix, int FirstCellRow, int FirstCellColumn, int blockSize, int bins);
extern int CalculateIndexTrueWindow(int windowCenter, int sizeImage, bool upper);

extern void InsertHistogramFeaturesVector(Image* HOG, double ***cellMatrix, int cellSize, int bins);
extern void HOG(Image *img, int cellSize, int blockSize, int bins, int band, Image *HOG);


//static void CreateHOGFeaturesVector(Image *HOG, double ***cellMatrix, int cellSize, int bins, Point windowCenter, Point WUL, Point WDR);
//static void InsertHistogramInFeatureVector(Image *HOG, double ***cellMatrix, Point cell, int bins, Point windowCenter, int index);

#endif