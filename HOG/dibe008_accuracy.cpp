/***************************************************************************
 *   Copyright (C) 2008 by Gabriele Moser                                  *
 *                                                                         *
 *   This code has been written for the OPERA "Civil protection from       *
 *   flooding events" project, funded by the Italian Space Agency (ASI).   *
 *                                                                         *
 *                                 *******                                 *
 *                                                                         *
 *                            ACCURACY (DIBE008)                           *
 *                                                                         *
 ***************************************************************************/

#include "image.h"
#include "library.h"

using namespace std;


double roundperc(double p)
{
	if(10000*p-floor(10000*p)>0.5)
		return 0.0001*ceil(10000*p);
	else
		return 0.0001*floor(10000*p);
}


double Image::ConfusionMatrix(Image *test, char name[80])
{
	word w=this->get_width(),h=this->get_height(),K=word(this->maxrange(0)),C=word(test->maxrange(0)),c,k,m,n;
	bool valid=true;
	if((this->get_bands()>1)||(test->get_bands()>1))
		errorlog("Image->ConfusionMatrix. Multichannel format unsupported for a map.");
	if((w!=test->get_width())||(h!=test->get_height()))
		errorlog("Image->ConfusionMatrix. Inconsistent image sizes.");
	if(K>C)
		valid=false;
	else if(K<C)
		K=C;
	Matrix H = Matrix(C,C+1,0);
	for(m=0;m<w;m++)
		for(n=0;n<h;n++)
			if(c=word(test->get(m,n,0)))
			{
				if(k=word(this->get(m,n,0)))
					H[c-1][k-1]++;
				else
					H[c-1][C]++;
			}
	//ofstream out(name);
	ofstream out(name,ios::app);
	out << endl << name;
	if(!valid) out << "casino!";
	out << endl;
	Vector PA = Vector(C,0);
	Vector UA = Vector(C,0);
	double OA=0,AA=0,kappa=0,tot=0;
	for(c=0;c<C;c++)
	{
		PA[c]=0;
		UA[c]=0;
		for(k=0;k<C+1;k++)
		{
			tot+=H[c][k];
			PA[c]+=H[c][k];
			if(k<C) UA[c]+=H[k][c];
			if(c==k) OA+=H[c][k];
			out << H[c][k] << ";";
		}
		kappa+=PA[c]*UA[c];
		if(PA[c]>0)
			PA[c]=H[c][c]/PA[c];
		if(UA[c]>0)
			UA[c]=H[c][c]/UA[c];
		AA+=PA[c];
		out << endl;
	}
	for(c=0;c<C;c++)
	{
		PA[c]=roundperc(PA[c]);
		UA[c]=roundperc(UA[c]);
	}
	if(tot*tot>kappa)
		kappa=roundperc((tot*OA-kappa)/(tot*tot-kappa));
	AA=roundperc(AA/C);
	if(tot>0)
		OA=roundperc(OA/tot);
	for(c=0;c<C;c++)
		out << PA[c] << ";";
	out << endl;
	for(c=0;c<C;c++)
		out << UA[c] << ";";
	out << endl << OA << ";" << AA << ";" << kappa << endl;
	out.close();

	ofstream outOA("OA.csv",ios::app),outAA("AA.csv",ios::app),outKA("Kappa.csv",ios::app);
	outOA << name << ";" << OA << endl;
	outAA << name << ";" << AA << endl;
	outKA << name << ";" << kappa << endl;
	outOA.close();
	outAA.close();
	outKA.close();
	return OA;
}


double Image::ConfusionMatrix(Image *test, char name[80], char header[256])
{
	word w=this->get_width(),h=this->get_height(),K=word(this->maxrange(0)),C=word(test->maxrange(0)),c,k,m,n;
	bool valid=true;
	if((this->get_bands()>1)||(test->get_bands()>1))
		valid=false;
	if((w!=test->get_width())||(h!=test->get_height()))
		valid=false;
	if(K>C)
	{
		valid=false;
		C=K;
	}
	else if(K<C)
		K=C;
	Matrix H = Matrix(C,C+1,0);
	for(m=0;m<w;m++)
		for(n=0;n<h;n++)
			if(c=word(test->get(m,n,0)))
			{
				if(k=word(this->get(m,n,0)))
					H[c-1][k-1]++;
				else
					H[c-1][C]++;
			}
	//ofstream out(name);
	ofstream out(name,ios::app);
	out << endl << header;
	if(!valid) out << ";casino!";
	out << endl;
	Vector PA = Vector(C,0);
	Vector UA = Vector(C,0);
	double OA=0,AA=0,kappa=0,tot=0;
	for(c=0;c<C;c++)
	{
		PA[c]=0;
		UA[c]=0;
		for(k=0;k<C+1;k++)
		{
			tot+=H[c][k];
			PA[c]+=H[c][k];
			if(k<C) UA[c]+=H[k][c];
			if(c==k) OA+=H[c][k];
			out << H[c][k] << ";";
		}
		kappa+=PA[c]*UA[c];
		if(PA[c]>0)
			PA[c]=H[c][c]/PA[c];
		if(UA[c]>0)
			UA[c]=H[c][c]/UA[c];
		AA+=PA[c];
		out << endl;
	}
	for(c=0;c<C;c++)
	{
		PA[c]=roundperc(PA[c]);
		UA[c]=roundperc(UA[c]);
	}
	if(tot*tot>kappa)
		kappa=roundperc((tot*OA-kappa)/(tot*tot-kappa));
	AA=roundperc(AA/C);
	if(tot>0)
		OA=roundperc(OA/tot);
	for(c=0;c<C;c++)
		out << PA[c] << ";";
	out << endl;
	for(c=0;c<C;c++)
		out << UA[c] << ";";
	out << endl << OA << ";" << AA << ";" << kappa << endl;
	out.close();

	ofstream outACC("Accuracy.csv",ios::app);
	outACC << header << ";" << OA << ";" << AA << ";" << kappa;
	if(!valid) outACC << ";casino!";
	outACC << endl;
	outACC.close();
	return OA;
}


void Image::BinaryError(Image *test, char name[80])
{
	word w=this->get_width(), h=this->get_height(),m,n,labelmap,labeltest;
	if((this->get_bands()>1)||(test->get_bands()>1))
		errorlog("Image->BinaryError. Multichannel format unsupported for a map.");
	if((w!=test->get_width())||(h!=test->get_height()))
		errorlog("Image->BinaryError. Inconsistent image sizes.");
	if((this->maxrange(0)>2)||(test->maxrange(0)>2))
		errorlog("Image->BinaryError. Inconsistent numbers of classes.");
	double PF=0,PD=0,Perr=0,testH1=0,testH0=0;
	for(m=0;m<w;m++)
		for(n=0;n<h;n++)
			if((labelmap=word(this->get(m,n,0)))&&(labeltest=word(test->get(m,n,0))))
			{
				if(labeltest==1)		// no-change in test map
				{
					testH0++;
					if(labelmap==2)
					{
						PF++;
						Perr++;
					}
				}
				else					// change in test map
				{
					testH1++;
					if(labelmap==2)
						PD++;
					else
						Perr++;
				}
			}
	if(testH0>0)
		PF=roundperc(PF/testH0);
	if(testH1>0)
		PD=roundperc(PD/testH1);
	if(testH0+testH1>0)
		Perr=roundperc(Perr/(testH0+testH1));
	ofstream out(name);
	out << PF << "	" << PD << "	" << Perr;
	out.close();
}


void Image::McNemar(Image *map, Image *test, char name[80])
{
	word w=this->get_width(),h=this->get_height(),C=word(test->maxrange(0)),c,m,n;
	if((this->get_bands()>1)||(test->get_bands()>1)||(map->get_bands()>1))
		errorlog("Image->McNemar. Multichannel format unsupported for a map.");
	if((w!=test->get_width())||(h!=test->get_height())||(w!=map->get_width())||(h!=map->get_height()))
		errorlog("Image->McNemar. Inconsistent image sizes.");
	if((C!=word(this->maxrange(0)))||(C!=word(map->maxrange(0))))
		errorlog("Image->McNemar. Unequal numbers of classes.");
	double thiserror=0,maperror=0;
	for(m=0;m<w;m++)
		for(n=0;n<h;n++)
			if(c=word(test->get(m,n,0)))
			{
				word thislabel=word(this->get(m,n,0)),maplabel=word(map->get(m,n,0));
				if((thislabel!=c)&&(maplabel==c)) thiserror++;
				else if((thislabel==c)&&(maplabel!=c)) maperror++;
			}
	double Z1=(thiserror-maperror)/sqrt(thiserror+maperror);
	double Z2=(fabs(thiserror-maperror)-1)*(fabs(thiserror-maperror)-1)/(thiserror+maperror);
	ofstream out(name);
	out << 100*roundperc(Z1*0.01) << "	" << 100*roundperc(Z2*0.01) << endl
		<< Queue(fabs(Z1)) << "	" << 1-IncGammaP(0.5,0.5*Z2);
	out.close();
}
