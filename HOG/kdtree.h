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

#ifndef KDTREE_H
#define KDTREE_H

typedef unsigned int uint;
const float inf=1e+08;


class record
{
public:
	float*	data;
	int		label;
	int		index;
	record() { data=0;}
	~record() { if(data) delete data;}
};


class kdnode
{
public:
	uint	d;
	float	p;
	uint	l;
	uint	u;
	kdnode*	left;
	kdnode*	right;
	kdnode(uint l, uint u);
	~kdnode();
};


class kdtree
{
	float*  pqd;
	uint*   pqr;
	record* Xq;
	float*  Bplus;
	float*  Bminus;
	uint*   perm;
	record* rset;
	int	    b;
	int     m;
	int	    dim;
    int     rnum;
	int 	fnode_examinated;
	int     rec_examinated;
	int     pqd_m;
	kdnode* root;

public:
	kdtree(record* recset, int rec_num, int bucket_size, int feat_num);
	~kdtree();
	uint*   search(record* query, int k);			//restituisce il pqr
	void	delete_pqr();
private:
	bool    search2(kdnode* node);
	bool    ball_within_bounds();
	bool    bounds_overlap_ball();
	void    delete_node (kdnode* node);
	kdnode* build_tree(uint l, uint u);
	float   median    (int feature, int l, int u);
	int     spreadest (int l, int u);
};

#endif
