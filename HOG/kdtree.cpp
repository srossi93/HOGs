/***************************************************************************
 *   Copyright (C) 2008 by Gabriele Moser                                  *
 *                                                                         *
 *   This code has been written for the OPERA "Civil protection from       *
 *   flooding events" project, funded by the Italian Space Agency (ASI).   *
 *                                                                         *
 *                                 *******                                 *
 *                                                                         *
 *                   RECORD, KDNODE, AND KDTREE CLASSES                    *
 *                                                                         *
 ***************************************************************************/

#include "image.h"
#include "library.h"
#include "kdtree.h"

using namespace std;


//#define SUBOPT

//piu' veloce che chiamare la pow2 di math.h visto che dobbiamo fare solo il quadrato
inline float pow2 (float a) {return a*a;}

/*********************************************
Funzioni per la ricerca in un kdtree
**********************************************/
//ritorna true se ricerca e' finita (il true si propaga a catena)
bool kdtree::search2(kdnode* node)
{
	float temp;
	float sum;
	uint i;
	int j,l;

	if (node->left == NULL && node->right == NULL)		//se nodo e' terminale
	{
		//fnode_examinated++;			//counter for performance
		//examine records in bucket, updating pqd, pqr
		for (i=node->l; i<=node->u; i++)
		{
			sum=0.0;
			for (j=0; j<dim; j++)
			{
				sum += pow2(rset[perm[i]].data[j] - Xq->data[j]);
				//rec_examinated++;			//counter for performance
				if (sum > pqd[0])
					break;				//se sono gia' oltre l'm-esimo record non calcolo le altre features
			}
			if (j==dim)					//se sono piu' vicino dell'm-esimo
			{							//aggiorno pqd e pqr
				for (l=0;l<m;l++)			//cerco dove inserire la nuova distanza
				{
					if (l==m-1 || sum > pqd[l+1])
					{	//inserisco qui
						pqd[l]=sum;
						pqr[l]=rset[perm[i]].index;
						break;
					}
					else
					{	//shifto
						pqd[l]=pqd[l+1];
						pqr[l]=pqr[l+1];
					}
				}
				if (pqd_m>0)					//aggiorno indice
					pqd_m--;
				if (pqd[pqd_m]==inf)
					pqd_m++;
			}
		}

#ifdef SUBOPT
        return true;
#else
		return ball_within_bounds();
#endif
	}

	//recursive call on closer son
	if ( Xq->data[node->d] <= node->p )
	{
		temp = Bplus[node->d];
		Bplus[node->d]=node->p;
		if ( search2(node->left) )
			return true;
		Bplus[node->d] = temp;
	}
	else
	{
		temp = Bminus[node->d];
		Bminus[node->d]=node->p;
		if ( search2(node->right) )
			return true;
		Bminus[node->d] = temp;
	}

	//recursive call on farther son, if necessary
	if ( Xq->data[node->d] <= node->p )
	{
		temp = Bminus[node->d];
		Bminus[node->d] = node->p;
		if ( bounds_overlap_ball() )
		{
			if ( search2(node->right) )
				return true;
		}
		Bminus[node->d] = temp;
	}
	else
	{
		temp = Bplus[node->d];
		Bplus[node->d] = node->p;
		if ( bounds_overlap_ball() )
		{
			if ( search2(node->left) )
				return true;
		}
		Bplus[node->d] = temp;
	}

	//see if we should return or terminate
	return ball_within_bounds();
}

//ritorna il pqr, un array con gli indici dei record piu' vicini alla query
uint* kdtree::search(record* query, int k)
{
	int i;

	m=k;

	if (pqd!=NULL)
	{
		delete[] pqd;
		pqd=NULL;
	}
	pqd = new float[m];

	if (pqr!=NULL)
	{
		delete[] pqr;
		pqr=NULL;
	}
	pqr = new uint[m];

	//inizializzazione
	for (i=0;i<m;i++)
	{
		pqd[i]=inf;
		pqr[i]=(uint)inf;
	}
	for (i=0;i<dim;i++)
	{
		Bplus[i]=inf;
		Bminus[i]=-inf;
	}
	fnode_examinated=0;
	rec_examinated=0;
	pqd_m=m-1;

	Xq=query;

	search2(root);


	delete[] pqd;
	pqd=NULL;
	return pqr;
}



bool kdtree::ball_within_bounds()
{
	for (int d=0;d<dim;d++)
	{
		if ( Xq->data[d] < Bminus[d] || pow2(Xq->data[d] - Bminus[d]) < pqd[pqd_m] || 
				Xq->data[d] > Bplus[d] || pow2(Xq->data[d] - Bplus[d]) < pqd[pqd_m])
					return false;
	}
	return true;
}

bool kdtree::bounds_overlap_ball()
{
	float sum =0.0;
	for (int d=0;d<dim;d++)
	{
		if (Xq->data[d] < Bminus[d])
		{
			//lower than low boundary
			sum += pow2(Xq->data[d] - Bminus[d]);
			if (sum >= pqd[pqd_m])
				return false;
		}
		else if (Xq->data[d] > Bplus[d])
		{
			//higher than high boundary
			sum += pow2(Xq->data[d] - Bplus[d]);
			if (sum >= pqd[pqd_m])
				return false;
		}
	}
	return true;
}


/*********************************************
costruttore di kdnode
*********************************************/
kdnode::kdnode (uint l, uint u)
{this->l=l;
 this->u=u;
 this->d=0;
 this->p=0;
 this->left=NULL;
 this->right=NULL; 
}

/*********************************************
distruttore di kdnode 
*********************************************/
kdnode::~kdnode ()
{}

/*********************************************
kdtree::median 
*********************************************/
float kdtree::median (int feature, int l, int u)
{return rset[perm[(l+u)/2]].data[feature];}


/*********************************************
kdtree::spreadest 
*********************************************/
int kdtree::spreadest (int l, int u)
{int feature,j,maxdim;
float max       =-999999999.0f,
      min       = 999999999.0f,
      maxspread =-999999999.0f;

//trova la feature con massima dispersione
for (feature=0; feature < dim; feature++)
    {max =-999999999.0f;
     min = 999999999.0f;
     for (j=l; j <= u; j++)
         {if (max < rset[perm[j]].data[feature]) max =   rset[perm[j]].data[feature];
          if (min > rset[perm[j]].data[feature]) min =   rset[perm[j]].data[feature];
    
          if (maxspread < fabs(max-min))
             {maxspread = float(fabs(max-min));
              maxdim = feature;}
         }
    }
return(maxdim);
}

/*********************************************
Costruttore di Kdtree 
*********************************************/
kdtree::kdtree(record* recset, int rec_num, int bucket_size, int feat_num)
{if (recset) rset=recset;
 else {cout << "ERROR: non è possibile costruire il KDtree a partire da un recordset nullo";
       exit(1);}
 rnum   = rec_num;
 dim    = feat_num;
 b      = bucket_size;
 Bplus  = new float[dim];
 Bminus = new float[dim];
 pqd=NULL;
 pqr=NULL;
 //inizializza perm
 perm = new uint[rnum];
 int i;
 for (i=0; i<rnum; i++)
     {perm[i]=i;}
 //costruisce l'albero
 root = build_tree(0, rnum-1);
}

/**********************************************
Distruttore di Kdtree
**********************************************/
kdtree::~kdtree()
{delete_node (root);
 root = NULL;
}


void kdtree::delete_node (kdnode* node)
{if (node == NULL) return;
 if (node->left != NULL ) 
     delete_node (node->left);
 if (node->right != NULL ) 
     delete_node (node->right);
 delete node;
 node = NULL;
}


/**********************************************
Funzione di utilita'
**********************************************/
void swap (unsigned int &a, unsigned int& b)
{unsigned int temp;
 temp = a;
 a = b;
 b = temp;
}


/*********************************************
kdtree::build_tree 
*********************************************/
kdnode* kdtree::build_tree(uint l, uint u)
{kdnode* node = new kdnode (l ,u);

 //check per nodi terminali
 if (int(u-l+1) <= b)  	// GAB-VC2008
    {//d, p unused on terminal nodes
     return node;}
 
 //trova la discriminator key (feature con massima dispersione)
 //e la partition             (mediana della feature con massima dispersione)
 int    d=0;    //feature con massima dispersione
 float p=0;    //mediana della feature con massima dispersione	// GAB-VC2008

 //calcola feature con max dispersione
 node->d = d = spreadest (l, u);

 //ordinamento dei dati
 float v;
 int i,j;
 int right=u;
 int left =l;
 while(right>left) {
    v=rset[perm[right]].data[d]; i=left-1; j=right;
    for (;;) {
       while (rset[perm[++i]].data[d] < v);
       while (rset[perm[--j]].data[d] > v && j>left); 
       if (i >= j) break;
       swap (perm[i], perm[j]);
  }
    swap (perm[i], perm [right]); 
    if (i>=int((l+u)/2)) right=i-1;	// GAB-VC2008
    if (i<=int((l+u)/2)) left=i+1;	// GAB-VC2008
 }
 
 //calcola mediana
 node->p = p = median(d,l,u);
 
 //costruzione dei nodi successori
 node -> left = build_tree  (l,(l+u)/2);
 node -> right = build_tree ((l+u)/2+1,u);
 return node;
}


void kdtree::delete_pqr()
{
	delete[] pqr;
	pqr=NULL;
}

