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

#include "image.h"
#include "library.h"

using namespace std;

///Default costructor with dimensions of image set to 0 and the poiter to image to NULL. Edited by Simone
Image::Image()
{
	width = 0;
	height = 0;
	bands = 0;
	pixel = NULL;
}

///Costructor with the dimensions of the image and the number of bands. Edited by Simone
Image::Image(word width, word height, word bands)
{
	if (width > 0 && height > 0 && bands > 0)
	{
		this->width = width;
		this->height = height;
		this->bands = bands;

		pixel = new float[width * height * bands];
		if (pixel == NULL)
		{
			errorlog("The image can't be create");
		}
		return;
	}
	errorlog("Image dimensions should be grater then 0");
}


///Defaul Destructor: set to 0 every pixel and to NULL the pointer to the image. Edited by Simone
//Image::~Image()
//{
//	if (pixel != NULL)
//	{
//		for (int i = 0; i < height*width*bands; i++)
//		{
//			pixel[i] = 0;
//		}
//		delete pixel;
//		pixel = 0;
//	}
//}


void Image::create(word w, word h, word b)
{
	width = w;
	height = h;
	bands = b;
	pixel = new float[w * h * b];
}


void Image::create(word w, word h, word b, double z)
{
	width=w;
	height=h;
	bands=b;
	pixel = new float[w*h*b];
	longword k;
	for (k = 0; k < longword(w*h*b); k++) pixel[k] = float(z);
}


void Image::kill()
{
	delete[] pixel;
}


void Image::load(char name[256], word w, word h, word b, basetype t)
{
	word m,n,r;
	this->create(w,h,b);
	if(t==ASCII)
	{
		ifstream in(name);
		if(in)
		{
			double Z;
			for(r=0;r<b;r++)
				for(n=0;n<h;n++)
					for(m=0;m<w;m++)
					{
						in >> Z;
						pixel[m+n*width+r*width*height]=float(Z);
					}
			in.close();
		}
		else
		{
			char errormessage[256]="Image->load(). Missing input file ";
			strcat(errormessage,name);
			errorlog(errormessage);
		}
	}
	else
	{
		ifstream in(name,ios::binary);
		if(in)
		{
			byte B=0;
			word W=0;
			float F=0;
			double D=0;
			for(r=0;r<b;r++)
				for(n=0;n<h;n++)
					for(m=0;m<w;m++)
					{
						switch(t)
						{
						case BYTE:
							in.read((char*)&B,sizeof(byte));
							pixel[m+n*width+r*width*height]=B;
							break;
						case WORD:
							in.read((char*)&W,sizeof(word));
							pixel[m+n*width+r*width*height]=W;
							break;
						case FLOAT:
							in.read((char*)&F,sizeof(float));
							pixel[m+n*width+r*width*height]=F;
							break;
						case DOUBLE:
							in.read((char*)&D,sizeof(double));
							pixel[m+n*width+r*width*height]=float(D);
							break;
						}
					}
			in.close();
		}
		else
		{
			char errormessage[256]="Image->load(). Missing input file ";
			strcat(errormessage,name);
			errorlog(errormessage);
		}
	}
}


void Image::store(char name[80], basetype t)
{
	word m,n,r;
	if(t==ASCII)
	{
		ofstream out(name);
		if(out)
		{
			for(r=0;r<bands;r++)
				for(n=0;n<height;n++)
					for(m=0;m<width;m++)
						out << pixel[m+n*width+r*width*height] << " ";
			out.close();
		}
		else
		{
			char errormessage[256]="Image->load(). Error opening output file ";
			strcat(errormessage,name);
			errorlog(errormessage);
		}
	}
	else
	{
		ofstream out(name,ios::binary);
		if(out)
		{
			for(r=0;r<bands;r++)
				for(n=0;n<height;n++)
					for(m=0;m<width;m++)
					{
						double D=pixel[m+n*width+r*width*height];
						byte B=byte(D);
						word W=word(D);
						float F=float(D);
						switch(t)
						{
						case BYTE:
							out.write((char*)&B,sizeof(byte));
							break;
						case WORD:
							out.write((char*)&W,sizeof(word));
							break;
						case FLOAT:
							out.write((char*)&F,sizeof(float));
							break;
						case DOUBLE:
							out.write((char*)&D, sizeof(double));
							break;
						}
				}
			out.close();
		}
		else
		{
			char errormessage[256]="Image->load(). Error opening output file ";
			strcat(errormessage,name);
			errorlog(errormessage);
		}
	}
}


void Image::mapstore(char name[80], netpbmformat format)
{
	if(bands>1)
		errorlog("Image->mapstore. Unsupported multichannel map.");
	ofstream out(name,ios::binary);
	word m,n,M=(width%8>0 ? width/8+1 : width/8);
	switch(format)
	{
	case PBM:
		out << "P4 " << width << " " << height << " ";
		for(n=0;n<height;n++)
			for(m=0;m<M;m++)
			{
				word B=0,i,power=128;
				for(i=0;(i<8)&&(8*m+i<width);i++)
				{
					if(pixel[8*m+i+n*width]>1) B+=power;
					power/=2;
				}
				out.write((char*)&B,sizeof(byte));
			}
		break;
	case PGM:
		out << "P5 " << width << " " << height << " " << this->maxrange(0) << " ";
		for(n=0;n<height;n++)
			for(m=0;m<width;m++)
			{
				byte B=byte(pixel[m+n*width]);
				out.write((char*)&B,sizeof(byte));
			}
		break;
	}
	out.close();
}


void Image::store(char name[80], word channel, basetype t)
{
	word m,n;
	Image *buffer = new Image;
	buffer->create(width,height,1);
	for(m=0;m<width;m++)
		for(n=0;n<height;n++)
			buffer->put(m,n,0,pixel[m+n*width+channel*width*height]);
	buffer->store(name,t);
	buffer->kill();
	delete buffer;
}


void Image::put(word x, word y, word band, double z)
{
	if((x>=width)||(y>=height)||(band>=bands))
		errorlog("Image->put. Out of range.");
	pixel[x+y*width+band*width*height]=float(z);
}


double Image::get(word x, word y, word band)
{
	if ((x >= width) || (y >= height) || (band >= bands))
		errorlog("Image->get. Out of range.");
	return pixel[x + y*width + band*width*height];
}



double Image::minrange(word r)
{
	word m,n;
	if(r>=bands)
		errorlog("Image->minrange. Out of range.");
	double M=DBL_MAX;
	for(m=0;m<width;m++)
		for(n=0;n<height;n++)
		if (pixel[m + n*width + r*width*height] < M) M = pixel[m + n*width + r*width*height];
	if(M==DBL_MAX)
		errorlog("Image->minrange. Invalid result.");
	return M;
}


double Image::maxrange(word r)
{
	word m,n;
	if(r>=bands)
		errorlog("Image->maxrange. Out of range.");
	double M=-DBL_MAX;
	for(m=0;m<width;m++)
		for(n=0;n<height;n++)
			if(pixel[m+n*width+r*width*height]>M) M=pixel[m+n*width+r*width*height];
	if(M==-DBL_MAX)
		errorlog("Image->minrange. Invalid result.");
	return M;
}


void Image::rescale(double A, double B, word r)
{
	word m,n;
	if((A>=B)||(r>=bands))
		errorlog("Image->rescale. Out of range.");
	double a=this->minrange(r),b=this->maxrange(r);
	for(m=0;m<width;m++)
		for(n=0;n<height;n++)
			if(b>a)
				pixel[m+n*width+r*width*height]=float((B-A)*(pixel[m+n*width+r*width*height]-a)/(b-a)+A);
			else
				pixel[m+n*width+r*width*height]=float(0.5*(A+B));
}


void Image::rescale(double A, double B)
{
	if(A>=B)
		errorlog("Image->rescale. Out of range.");
	word r;
	for(r=0;r<bands;r++)
		this->rescale(A,B,r);
}



void Image::threshold(Image *img, double T)
{
	if(img->get_bands()>1)
		errorlog("Image->threshold. Thresholding undefined for multichannel images.");
	word w=img->get_width(),h=img->get_height(),m,n;
	this->create(w,h,1);
	for(m=0;m<w;m++)
		for(n=0;n<h;n++)
			if(img->get(m,n,0)>T)
				this->put(m,n,0,2);
			else
				this->put(m,n,0,1);
}


word Image::stdmaker()
{
	// occhio! assume che non si siano più di 255 classi (ragionevole!)
	if(bands>1)
		errorlog("Image->stdmaker. Map standardization undefined on multichannel images.");
	word m,n,i;
	vector<bool> nonzero = vector<bool>(255,false);
	for(m=0;m<width;m++)
		for(n=0;n<height;n++)
			if(i=word(pixel[m+n*width]))
			{
				if(i<256)
					nonzero[i-1]=true;
				else
					errorlog("Image->stdmaker. Pixel values out of range for a map.");
			}
	vector<byte> newldg = vector<byte>(255,0);
	byte b=0;
	for(i=0;i<255;i++)
		if(nonzero[i])
		{
			b++;
			newldg[i]=b;
		}
	for(m=0;m<width;m++)
		for(n=0;n<height;n++)
			if(i=word(pixel[m+n*width]))
				pixel[m+n*width]=newldg[i-1];
	return word(b);
}
