/***************************************************************************
 *   Copyright (C) 2008 by Gabriele Moser                                  *
 *                                                                         *
 *   This code has been written for the OPERA "Civil protection from       *
 *   flooding events" project, funded by the Italian Space Agency (ASI).   *
 *                                                                         *
 *                                 *******                                 *
 *                                                                         *
 *                               IMAGE CLASS                               *
 *                                                                         *
 ***************************************************************************/
#ifndef IMG_LIB_MOSER
#define IMG_LIB_MOSER

typedef unsigned char byte; //8 bit, unsigned
typedef unsigned short word; //16 bit, unsigned
typedef unsigned long longword;
enum basetype { ASCII, BYTE, WORD, LONG, FLOAT, DOUBLE };
enum netpbmformat { PBM, PGM, PPM };

class Image
{

protected:
	float *pixel;	///< Pointer to the image in the memory
	word width;		///< Width of the image, in pixels 
	word height;	///< Height of the image, in pixels
	word bands;		///< Number of channels  

public:
	Image();									///< Defalut costructor 
	Image(word width, word height, word bands); ///< Basic Costructor 

	//~Image(); ///< Basic Destructor

	void create(word Width, word Height, word Bands);									// alloca spazio in memoria
	void create(word Width, word Height, word Bands, double init_value);				// idem ma inizializza anche tutto ad un certo valore
	void kill();																		// disalloca
	void load(char filename[256], word Width, word Height, word Bands, basetype type);	// carica da disco in formato RAW-BSQ
	void store(char filename[80], basetype type);										// salva su disco, stesso formato
	void store(char filename[80], word channel, basetype type);
	void put(word x, word y, word band, double z);
	double get(word x, word y, word band);
	word get_width() const { return width; }
	word get_height() const { return height; }
	word get_bands() const { return bands; }

	double minrange(word band);
	double maxrange(word band);
	void rescale(double outmin, double outmax, word band);
	void rescale(double outmin, double outmax);

	word stdmaker();
	double ConfusionMatrix(Image *test_map, char filename[80]);
	double ConfusionMatrix(Image *test_map, char filename[80], char headername[256]);
	void BinaryError(Image *test_map, char filename[80]);
	void McNemar(Image *ref_map, Image *test_map, char filename[80]);

	void mapstore(char filename[80], netpbmformat format);
	void threshold(Image *image, double given_threshold);


};

#endif