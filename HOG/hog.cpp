/***************************************************************************
 *                                                                         *
 *                                                                         *
 *                                                                         *
 *                LIBRARY OF SUPPORT FOR REALIZING HOG                     *
 *                                                                         *
 *                      Simone & Eugenio & Margherita                      *
 *                                                                         *
 *                                                                         *
 *                                                                         *
 ***************************************************************************/

//#include <cstdio>
//#include <iostream>
//#include <fstream>
//#include <cmath>
//#include <cfloat>
//#include <cstdlib>
//#include <ctime>
//#include <vector>
//#include <limits>

#include "image.h"

#include "library.h"

#include "hog.h"

using namespace std;

#define quad(a) (double)a*a
#define epsilon 0.00000001


// Funzione per unire tre immagini a singola banda in un unica immagine a tre bande
extern Image* MergeImg(Image *I, Image *hog0, Image *hog1, Image *hog2)
{
	Image *mImg = new Image;
	//creo la nuova immagine a tre bande
	int band = I->get_bands() + hog0->get_bands() + hog1->get_bands() + hog2->get_bands();

	mImg->create(I->get_width(), I->get_height(), band);
	                  
	//copio l'Immagine (prime tre bande)
	for (int b = 0; b < I->get_bands(); b++)
	{
		for (int w = 0; w < I->get_width(); w++)
		{
			for (int h = 0; h < I->get_height(); h++)
			{
				mImg->put(w, h, b, I->get(w, h, b));
			}
		}
	}
	Image *arrBands[3] = { hog0, hog1, hog2 };
	for (int k = 0; k < 3; k++)
	{
		for (int b = 0; b < hog0->get_bands(); b++)                                              //scorro le bande della nuova immagine
		{
			for (int i = 0; i < I->get_width(); i++)                                       //scorro width
			{
				for (int j = 0; j < I->get_height(); j++)                                   //scorro height
				{
					mImg->put(i, j, 3+b*k, arrBands[k]->get(i, j, b));                //setto il corrispondente pixel della mergeImg copiando banda per banda
				}
			}
		}
	}
	return mImg;
}



//Calcolo del gradiente di un'immagine con kernel [-1 0 1]
//Parametri
//-Image *img: riferimento all'immagine sulla quale calcolare il gradiente 
//-int band: indice della banda su cui calcolare il gradiente
//Return
//-Image *: puntantore all'immagine creata internamente con le informazioni
//          estratte dall'immagine di partenza
extern Image *Gradient(Image *img, int band)
{
	Image *nImg = new Image;
	nImg->create(img->get_width(), img->get_height(), 1, 0);  ///la creo con le dimensioni dell'immagine passata per farla funzionare con la trasposta
	float val = 0;  //il valore deve essere float
	int ker[] = { -1, 0, 1 };

	for (int i = 0; i < img->get_height(); i++)
	{
		for (int j = 0; j < img->get_width(); j++)
		{
			for (int k = 0; k < 3; k++)
			{
				if (j == 0 && k == 0)
				{
					//j++;
					continue;
				}
				else if (j == (img->get_width() - 1) && k == 2)
				{
					//j++;
					continue;
				}
				else
				{
					//j++;
					val += (float)img->get(j + ker[k], i, band) * ker[k];
				}
			}
			nImg->put(j, i, 0, val);

			val = 0;
		}

	}

	//nImg->store("Gradient.raw", BYTE);
	return nImg;
}

//Calcolo della trasposta di un'immagine
//Parametri
//-Image *img: riferimento all'immagine di partenza
//Return
//-Image *: puntatore ad un'immagine creata internamente  
extern Image *Transpose(Image *img)
{
	Image *tImg = new Image;
	tImg->create(img->get_height(), img->get_width(), img->get_bands()); ///creo immagine con x e y invertiti (rif. immagine passata)

	for (int b = 0; b < img->get_bands(); b++)//scorro le bande
	{
		for (int i = 0; i < img->get_width(); i++)//scorro width
		{
			for (int j = 0; j < img->get_height(); j++) //scorro height
			{
				tImg->put(j, i, b, img->get(i, j, b)); //copio il pixel nella nuova img con banda b e indici invertiti
			}
		}
	}
	return tImg;
}

//Calcolo della magnitudine e fase del gradiente a partire dai valori dei 
//gradienti su x e y
//Parametri
//-Image *imgX e Image *imgY: riferimenti ai gradienti su righe e colonne
//-Image *magnI e Imgae *phI: immagini create internamente per lo storage 
//di magnitudine e fase. Al termine del lavoro verrano modificate direttamente
//queste immagini
//Reminder: la fase è calcolata, in prima battuta, in gradi da -90° a +90°
//successivamente, solo per quelle negative, si opera una traslazione di 180°
extern void MagnitudePh(Image *imgX, Image *imgY, Image *magnI, Image *phI)
{
	magnI->create(imgX->get_width(), imgX->get_height(), 1);
	phI->create(imgX->get_width(), imgX->get_height(), 1);

	for (int i = 0; i < imgX->get_height(); i++)
	{
		for (int j = 0; j < imgX->get_width(); j++)
		{
			double value = sqrt(quad(imgX->get(j, i, 0)) + quad(imgY->get(j, i, 0)));
			magnI->put(j, i, 0, value);
			//value = atan2((double)imgX->get(j, i, 0),(double) imgY->get(j, i, 0));
			value = (atan((double)imgY->get(j, i, 0) / (double)imgX->get(j, i, 0)))*(180 / PI);
			//value = value + 90;
			if (value < 0)
			{
				value = value + 2 * 90;
			}
			phI->put(j, i, 0, value/* *(255/180) or (255/3.14)*/); // per visualizzarla in modo corretto se no nera, valori di radianti troppo piccoli
		}
	}
	return;

}

//Creazione di una cella
//Parametri
//-Image *img: immagine di partenza sulla quale si estrae la cella
//-int size: dimensione in pixel della cella (ipotesi: cella quadrata)
//-int cellRow, cellColumn: riferimento a quale cella si vuole estrarre
//Return
//-Image *: riferimento alla cella creata internamente  
extern Image *CreateCell(Image * img, int size, int cellRow, int cellColumn)
{
	Image *cell = new Image;
	cell->create(size, size, 1);

	int indexImageX = size*cellColumn;
	int indexImageY = size*cellRow;

	for (int y = 0; y < size; y++)
	{
		for (int x = 0; x < size; x++)
		{
			cell->put(x, y, 0, img->get(indexImageX + x, indexImageY + y, 0));
		}
	}
	return cell;
}

extern double * CreateIstogram(Image *cellPh, Image *cellMagn, int bins)
{
	//const int i = bins;

	double *istogram;// [9];// [i];
	istogram = (double*)malloc(sizeof(double)*bins); //struttura dinamica dell'array. non si puo fare con variabile. serve la free?
	//****** Simo: fare controllo su malloc

	for (int i = 0; i < bins; i++)
	{
		istogram[i] = 0;
	}


	for (int y = 0; y < cellPh->get_height(); y++)
	{
		for (int x = 0; x < cellPh->get_width(); x++)
		{
			for (int k = 0; k < bins; k++)
			{
				double angle = cellPh->get(x, y, 0);//la funzione atan restituisce radianti

				if (angle < (180 / bins) * (k + 1) && angle >= (180 / bins)*k)  //change order 
				{
					istogram[k] += cellMagn->get(x, y, 0);
					break;
				}
			}
		}
	}

	return istogram;
}

extern void NormalizeBlock(double ***cellMatrix, int FirstCellRow, int FirstCellColumn, int blockSize, int bins)
{
	double normVal = 0;
	for (int x = 0; x < blockSize; x++)
	{
		for (int y = 0; y < blockSize; y++)
		{
			normVal += sqrt(Norm2Calculator(cellMatrix[FirstCellRow + x][FirstCellColumn + y], bins) + quad(epsilon));			//Ho calcolato il denominatore di normalizzazione
		}
	}

	//divido ogni elemento degli istogrammi delle celle del blocco per il valore di normalizzazione calcolato su (L2-NORM)

	for (int x = 0; x < blockSize; x++)
	{
		for (int y = 0; y < blockSize; y++)
		{
			for (int z = 0; z < bins; z++)			{
				//cout << cellMatrix[x][y][z] << endl;
				cellMatrix[FirstCellRow + x][FirstCellColumn + y][z] /= normVal;
				//cout << cellMatrix[x][y][z] << endl;
			}
			//cout << endl;
		}
	}

}

extern double Norm2Calculator(double arr[], int dim)
{
	double norm2 = 0;
	/*int zise = arr.size()*/;
	for (int i = 0; i < dim; i++)	{
    norm2 += quad(arr[i]);
	}
	return norm2;
}



extern void InsertHistogramFeaturesVector(Image* HOG, double ***cellMatrix, int cellSize, int bins)
{
	int countrow = 0;
	int countcolumn = 0;
	for (int row = 0; row < HOG->get_height(); row++)
	{
		for (int column = 0; column < HOG->get_width(); column++)
		{
			for (int b = 0; b < bins; b++)
			{
				HOG->put(column, row, b, cellMatrix[countrow][countcolumn][b]*1000);
			}
			if (column != 0 && (column+1)%cellSize == 0)
				countcolumn++;
		}
		countcolumn = 0;

		if (row != 0 && (row+1)%cellSize == 0)
			countrow++;
	}
}

extern void HOG(Image *I, int cellSize, int blockSize, int B, int band, Image *HOG)
{
	//IGX0 contiene il gradiente su x dell'immagine I
	Image *IGX0 = new Image();

	//IGY0 contiene il gradiente su y dell'immagine I
	Image *IGY0 = new Image();

	//I_TRASP è un'immagine temporanea che contiene la trasposta dell'immagine originale (serve per 
	//calcolare il gradiente su y usando la stessa funzione di prima sapendo che (I^T*K)^T=I*K^T
	//dove I è l'immagine di partenza e K è il Kernel base [-1 0 1] per il calcolo del gradiente
	Image *I_TRASP = new Image();

	///Calcolo del gradiente su x
	IGX0 = Gradient(I, band);

	//Calcolo del gradiente su y
	I_TRASP = Transpose(I);
	IGY0 = Gradient(I_TRASP, band);
	IGY0 = Transpose(IGY0);
	I_TRASP->kill();
	delete I_TRASP;


	//gradientMagn e gradientPhase conterranno, rispettivamente, la magnitudine e la phase del 
	//gradiente appena calcolato
	Image *gradientMagn = new Image(), *gradientPhase = new Image();

	//Calcolo della magnitudine e fase del gradiente. Al termine della computazione le due immagini
	//IGX0 e IGY0 sono cancellate, in quanto non più utili
	MagnitudePh(IGX0, IGY0, gradientMagn, gradientPhase);
	IGX0->kill();
	IGY0->kill();
	delete IGX0, IGY0;


	//Dichiarazione e inizializzazione di alcuni parametri per il calcolo dell'HOG 

	int N = I->get_height() / cellSize;	// Numero di righe della matrice, Numero di cella su una colonna
	int M = I->get_width() / cellSize;	// Numero di colonne della matrice, Numero di celle su una riga 



	//Simone alloca magicamente un cubo <3
	//Allocazione della matrice di istogrammi
	//Uso: cellMatrix[x][y][]
	//Per ogni coordinata sulla matrice viene restituito un puntatore a double (double*) 
	//ovvero l'istogramma per quella determinata cella
	double ***cellMatrix;

	//Allocazione dinamica della matrice di istogrammi 
	//con l'aggiunta di controlli su effettiva riuscita dell'allocazione
	cellMatrix = (double***)malloc(N*sizeof(double**));
	if (cellMatrix == NULL)
	{
		errorlog("Allocation Error -> CellMatric can't be created");
	}
	for (int k = 0; k < N; k++)
	{
		cellMatrix[k] = (double**)malloc(M*sizeof(double*));
		if (cellMatrix[k] == NULL)
		{
			errorlog("Allocation Error -> *CellMatric can't be created");
		}
		for (int i = 0; i < M; i++)
		{
			cellMatrix[k][i] = (double*)malloc(B*sizeof(double));
			if (cellMatrix[k][i] == NULL)
			{
				errorlog("Allocation Error -> **CellMatric can't be created");
			}
		}
	}

	//Riempimento della matrice di array con gli istogrammi.
	//Per ogni cella, si demanda il calcolo dell'istogramma 
	//a una specifica funzione 
	//Reminder: k scorre tra su x (in orizzontale), i scorre su y (verticale)
	for (int i = 0; i < N; i++)
	{
		for (int k = 0; k < M; k++)
		{
			//Dichiarazione e inizializzazione di due immagini ausiliarie per il calcolo
			//della cella (una per la Magnitudine e una per la Fase)
			Image *cellM = new Image(), *cellP = new Image();
			cellM = CreateCell(gradientMagn, cellSize, i, k);
			cellP = CreateCell(gradientPhase, cellSize, i, k);

			//Calcolo dell'istogramma per la cella k-i
			cellMatrix[i][k] = CreateIstogram(cellP, cellM, B);

			//Eliminazione delle due immagini temporanee
			cellM->kill();
			cellP->kill();
			delete cellM, cellP;
		}
	}

	//Normalizzazione di tutti gli istogrammi blocco per blocco
	//Per evitare di accedere a zone di memoria fuori da quella allocata per la matrice,
	//a causa di dimensioni della stessa non definibili a priori, ferminamo il calcolo 
	//della normalizzazione una riga prima e una colonna prima. 
	//Reminder: l'ultima riga e l'ultima colonna della matrice non sono normalizzate
	for (int x = 0; x < N - 1; x += blockSize)
	{
		for (int y = 0; y < M - 1; y += blockSize)
		{
			NormalizeBlock(cellMatrix, x, y, blockSize, B);
		}
	}
 
	//Image *HOG = new Image(I->get_width(), I->get_height(), 9);
	InsertHistogramFeaturesVector(HOG, cellMatrix, 8, B);



	gradientMagn->kill();
	gradientPhase->kill();


	delete gradientMagn, gradientPhase;


}

//static int CalculateIndexTrueWindow(int windowCenter, int sizeImage, bool upper)
//{
//	int windowIndex = windowCenter;
//	for (int i = 0; i < (WSIZE - 1) / 2; i++)
//	{
//		if (windowIndex - (int)upper*sizeImage == 0)
//		{
//			break;
//		}
//		windowIndex = windowIndex - 1 + 2 * (int)upper;
//	}
//	return windowIndex;
//}


// riempimento descrittori basato su finestra mobile
//static void CreateHOGFeaturesVector(Image *HOG, double ***cellMatrix, int cellSize, int bins, Point windowCenter, Point WUL, Point WDR)
//{
//	Point startCell = { WUL.x, WUL.y };
//	while (startCell.x%cellSize != 0)
//	{
//		startCell.x++;
//	}
//	while (startCell.y%cellSize != 0)
//	{
//		startCell.y++;
//	}
//	InsertHistogramInFeatureVector(HOG, cellMatrix, startCell, bins, windowCenter, 0);
//
//	Point cellUp = startCell;
//	Point cellDown = { startCell.x + 8, startCell.y + 8 };
//	
//	int counter = 1;
//	while (cellUp.y > WUL.y && cellDown.y < WDR.y)
//	{
//		while (cellUp.x > WUL.x && cellDown.x < WDR.x)
//		{
//			InsertHistogramInFeatureVector(HOG, cellMatrix, cellUp, bins, windowCenter, counter*bins);
//			counter++;
//			cellUp.x += cellSize;
//			cellDown.x += cellSize;
//		}
//		cellUp.y += cellSize;
//		cellDown.y += cellSize;
//	}
//}

//static void InsertHistogramInFeatureVector(Image *HOG, double ***cellMatrix, Point cell, int bins, Point windowCenter, int index)
//{
//	for (int i = 0; i < bins; i++)
//	{
//		double value = cellMatrix[cell.x / 8][cell.y / 8][i];
//		HOG->put(windowCenter.x, windowCenter.y, index + i, value);
//	}
//
//}