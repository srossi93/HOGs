 /**************************************************************************
 *   Copyright (C) 2011 by Gabriele Moser                                  *
 *                                                                         *
 *   This code has been written for the AO-CSK "Development and validation *
 *   of multitemporal image analysis methodologies for multirisk           *
 *   monitoring of critical structures and infrastructures" project,       *
 *   funded by the Italian Space Agency (ASI) by modifying the "MRF energy *
 *   minimization software," Version 1.6, May 5, 2006, by O. Veksler,      *
 *   R. Zabih, Y. Boykov, and V. Kolmogorov. Copyright information for     *
 *   this package is below.                                                *
 *                                                                         *
 *                                 *******                                 *
 *                                                                         *
 *			   MODIFIED MRF ENERGY MINIMIZATION BY GRAPH-CUTS              *
 *                                                                         *
 ***************************************************************************
 * Copyright Olga Veksler, Ramin Zabih, and Vladimir Kolmogorov            *
 * Send any questions to olga@csd.uwo.ca, rdz@cs.cornell.edu,              *
 * vnk@microsoft.com                                                       *
 ***************************************************************************
 * (C) 2002 Marshall Tappen, MIT AI Lab                                    *
 ***************************************************************************
 * Copyright 2001 Vladimir Kolmogorov (vnk@cs.cornell.edu), Yuri Boykov    *
 * (yuri@csd.uwo.ca).                                                      *
 * This program is free software; you can redistribute it and/or modify    *
 * it under the terms of the GNU General Public License as published by    *
 * the Free Software Foundation; either version 2 of the License, or       *
 * (at your option) any later version.                                     *
 * This program is distributed in the hope that it will be useful,         *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 * GNU General Public License for more details.                            *
 * You should have received a copy of the GNU General Public License       *
 * along with this program; if not, write to the Free Software             *
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA *
 *                                                                         *
 ***************************************************************************/


#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <float.h>
#include "graphcut.h"


/***************************************************************************************************************************************/
void MRFmodel::initialize()
{
    m_initialized = 1;
    if ( m_dataType == ARRAY )
        setData(m_e->m_dataCost->m_costArray);
    else  setData(m_e->m_dataCost->m_costFn);
    if ( m_smoothType == FUNCTION )
        setSmoothness(m_e->m_smoothCost->m_costFn);
    else 
    {
        if ( m_smoothType == ARRAY )
        {
            checkArray(m_e->m_smoothCost->m_V);
            setSmoothness(m_e->m_smoothCost->m_V);
        }
        else
        {
            int smoothExp = m_e->m_smoothCost->m_smoothExp;
            if ( smoothExp != 1 && smoothExp != 2) 
            { 
                fprintf(stderr, "Wrong exponent in setSmoothness()!\n"); 
                exit(1); 
            }
            setSmoothness(smoothExp,m_e->m_smoothCost->m_smoothMax,m_e->m_smoothCost->m_lambda);
        }
        if (m_e->m_smoothCost->m_varWeights )
        {
            if (!m_grid_graph) 
            { 
                fprintf(stderr, "Edge multiplier cannot be used with non-grid graphs\n"); 
                exit(1); 
            }
            setCues(m_e->m_smoothCost->m_hWeights,m_e->m_smoothCost->m_vWeights);
        }
    }
}


void MRFmodel::commonInitialization(EnergyFunction *e)
{
    m_dataType    = e->m_dataCost->m_type;
    m_smoothType  = e->m_smoothCost->m_type;
    m_varWeights  = e->m_smoothCost->m_varWeights;
    m_initialized = 0;
	m_second_order = (e->m_order==SECOND_ORDER);
    m_e = e;
}


MRFmodel::MRFmodel(int width, int height, int nLabels, EnergyFunction *e)
{
    m_width        = width;
    m_height       = height;
    m_nLabels      = nLabels;
    m_nPixels      = width*height;
    m_grid_graph   = 1;
    commonInitialization(e);
}


MRFmodel::MRFmodel(int nPixels, int nLabels, EnergyFunction *e)
{
    m_nLabels     = nLabels;
    m_nPixels     = nPixels;
    m_grid_graph  = 0;
    commonInitialization(e);
}


char MRFmodel::checkEnergy()
{ 
    if ( !m_initialized) {printf("Call initialize() first,exiting!");exit(1);}
    return(2);
}


MRFmodel::EnergyVal MRFmodel::totalEnergy()
{
    if (!isValid()) { fprintf(stderr, "totalEnergy() cannot be called for invalid energy!\n"); exit(1); }
    if (!m_initialized) { fprintf(stderr, "Call initialize() first!\n"); exit(1); }
    return dataEnergy()+smoothnessEnergy();
}


void MRFmodel::optimize(int nIterations, float& time)
{
    if (!isValid()) { fprintf(stderr, "optimize() cannot be called for invalid energy!\n"); exit(1); }
    if (!m_initialized ) { fprintf(stderr, "run initialize() first!\n"); exit(1); }
    // start timer
    clock_t start = clock();
    optimizeAlg(nIterations);
    // stop timer
    clock_t finish = clock();
    time = (float) (((double)(finish - start)) / CLOCKS_PER_SEC);
}


void MRFmodel::optimize(int nIterations)
{
	float dummy_time;
	this->optimize(nIterations,dummy_time);
}


void MRFmodel::checkArray(CostVal *V)
{
    int i, j;
    for (i=0; i< m_nLabels; i++)
        for (j=i; j<m_nLabels; j++)
        if (V[i*m_nLabels + j] != V[j*m_nLabels + i])
        { 
            fprintf(stderr, "Error in setSmoothness(): V is not symmetric!\n"); 
            exit(1); 
        }
}


/***************************************************************************************************************************************/
Graph::Graph(void (*err_function)(char *))
{
    error_function = err_function;
    node_block_first = NULL;
    arc_for_block_first = NULL;
    arc_rev_block_first = NULL;
    flow = 0;
}


Graph::~Graph()
{
    while (node_block_first)
    {
        node_block *next = node_block_first -> next;
        delete node_block_first;
        node_block_first = next;
    }

    while (arc_for_block_first)
    {
        arc_for_block *next = arc_for_block_first -> next;
        delete arc_for_block_first -> start;
        arc_for_block_first = next;
    }

    while (arc_rev_block_first)
    {
        arc_rev_block *next = arc_rev_block_first -> next;
        delete arc_rev_block_first -> start;
        arc_rev_block_first = next;
    }
}


Graph::node_id Graph::add_node()
{
    node *i;
    if (!node_block_first || node_block_first->current+1 > &node_block_first->nodes[NODE_BLOCK_SIZE-1])
    {
        node_block *next = node_block_first;
        node_block_first = (node_block *) new node_block;
        if (!node_block_first) { if (error_function) (*error_function)("Not enough memory!"); exit(1); }
        node_block_first -> current = & ( node_block_first -> nodes[0] );
        node_block_first -> next = next;
    }
    i = node_block_first -> current ++;
    i -> first_out = (arc_forward *) 0;
    i -> first_in = (arc_reverse *) 0;
    i -> tr_cap = 0;
    return (node_id) i;
}


void Graph::add_edge(node_id from, node_id to, captype cap, captype rev_cap)
{
    arc_forward *a_for;
    arc_reverse *a_rev;
    if (!arc_for_block_first || arc_for_block_first->current+1 > &arc_for_block_first->arcs_for[ARC_BLOCK_SIZE])
    {
        arc_for_block *next = arc_for_block_first;
        char *ptr = new char[sizeof(arc_for_block)+1];
        if (!ptr) { if (error_function) (*error_function)("Not enough memory!"); exit(1); }
        if ((int)ptr & 1) arc_for_block_first = (arc_for_block *) (ptr + 1);
        else              arc_for_block_first = (arc_for_block *) ptr;
        arc_for_block_first -> start = ptr;
        arc_for_block_first -> current = & ( arc_for_block_first -> arcs_for[0] );
        arc_for_block_first -> next = next;
    }
    if (!arc_rev_block_first || arc_rev_block_first->current+1 > &arc_rev_block_first->arcs_rev[ARC_BLOCK_SIZE])
    {
        arc_rev_block *next = arc_rev_block_first;
        char *ptr = new char[sizeof(arc_rev_block)+1];
        if (!ptr) { if (error_function) (*error_function)("Not enough memory!"); exit(1); }
        if ((int)ptr & 1) arc_rev_block_first = (arc_rev_block *) (ptr + 1);
        else              arc_rev_block_first = (arc_rev_block *) ptr;
        arc_rev_block_first -> start = ptr;
        arc_rev_block_first -> current = & ( arc_rev_block_first -> arcs_rev[0] );
        arc_rev_block_first -> next = next;
    }
    a_for = arc_for_block_first -> current ++;
    a_rev = arc_rev_block_first -> current ++;
    a_rev -> sister = (arc_forward *) from;
    a_for -> shift  = (int) to;
    a_for -> r_cap = cap;
    a_for -> r_rev_cap = rev_cap;
    ((node *)from) -> first_out =
        (arc_forward *) ((int)(((node *)from) -> first_out) + 1);
    ((node *)to) -> first_in =
        (arc_reverse *) ((int)(((node *)to) -> first_in) + 1);
}


void Graph::set_tweights(node_id i, captype cap_source, captype cap_sink)
{
    flow += (cap_source < cap_sink) ? cap_source : cap_sink;
    ((node*)i) -> tr_cap = cap_source - cap_sink;
}


void Graph::add_tweights(node_id i, captype cap_source, captype cap_sink)
{
    register captype delta = ((node*)i) -> tr_cap;
    if (delta > 0) cap_source += delta;
    else           cap_sink   -= delta;
    flow += (cap_source < cap_sink) ? cap_source : cap_sink;
    ((node*)i) -> tr_cap = cap_source - cap_sink;
}


/*
    Converts arcs added by 'add_edge()' calls
    to a forward star graph representation.

    Linear time algorithm.
    No or little additional memory is allocated
    during this process
    (it may be necessary to allocate additional
    arc blocks, since arcs corresponding to the
    same node must be contiguous, i.e. be in one
    arc block.)
*/
void Graph::prepare_graph()
{
    node *i;
    arc_for_block *ab_for, *ab_for_first;
    arc_rev_block *ab_rev, *ab_rev_first, *ab_rev_scan;
    arc_forward *a_for;
    arc_reverse *a_rev, *a_rev_scan, a_rev_tmp;
    node_block *nb;
    bool for_flag = false, rev_flag = false;
    int k;
    if (!arc_rev_block_first)
    {
        node_id from = add_node(), to = add_node();
        add_edge(from, to, 1, 0);
    }
    /* FIRST STAGE */
    a_rev_tmp.sister = NULL;
    for (a_rev=arc_rev_block_first->current; a_rev<&arc_rev_block_first->arcs_rev[ARC_BLOCK_SIZE]; a_rev++)
    {
        a_rev -> sister = NULL;
    }
    ab_for = ab_for_first = arc_for_block_first;
    ab_rev = ab_rev_first = ab_rev_scan = arc_rev_block_first;
    a_for = &ab_for->arcs_for[0];
    a_rev = a_rev_scan = &ab_rev->arcs_rev[0];
    for (nb=node_block_first; nb; nb=nb->next)
    {
        for (i=&nb->nodes[0]; i<nb->current; i++)
        {
            /* outgoing arcs */
            k = (int) i -> first_out;
            if (a_for + k > &ab_for->arcs_for[ARC_BLOCK_SIZE])
            {
                if (k > ARC_BLOCK_SIZE) { if (error_function) (*error_function)("# of arcs per node exceeds block size!"); exit(1); }
                if (for_flag) ab_for = NULL;
                else          { ab_for = ab_for -> next; ab_rev_scan = ab_rev_scan -> next; }
                if (ab_for == NULL)
                {
                    arc_for_block *next = arc_for_block_first;
                    char *ptr = new char[sizeof(arc_for_block)+1];
                    if (!ptr) { if (error_function) (*error_function)("Not enough memory!"); exit(1); }
                    if ((int)ptr & 1) arc_for_block_first = (arc_for_block *) (ptr + 1);
                    else              arc_for_block_first = (arc_for_block *) ptr;
                    arc_for_block_first -> start = ptr;
                    arc_for_block_first -> current = & ( arc_for_block_first -> arcs_for[0] );
                    arc_for_block_first -> next = next;
                    ab_for = arc_for_block_first;
                    for_flag = true;
                }
                else a_rev_scan = &ab_rev_scan->arcs_rev[0];
                a_for = &ab_for->arcs_for[0];
            }
            if (ab_rev_scan)
            {
                a_rev_scan += k;
                i -> parent = (arc_forward *) a_rev_scan;
            }
            else i -> parent = (arc_forward *) &a_rev_tmp;
            a_for += k;
            i -> first_out = a_for;
            ab_for -> last_node = i;

            /* incoming arcs */
            k = (int) i -> first_in;
            if (a_rev + k > &ab_rev->arcs_rev[ARC_BLOCK_SIZE])
            {
                if (k > ARC_BLOCK_SIZE) { if (error_function) (*error_function)("# of arcs per node exceeds block size!"); exit(1); }
                if (rev_flag) ab_rev = NULL;
                else          ab_rev = ab_rev -> next;
                if (ab_rev == NULL)
                {
                    arc_rev_block *next = arc_rev_block_first;
                    char *ptr = new char[sizeof(arc_rev_block)+1];
                    if (!ptr) { if (error_function) (*error_function)("Not enough memory!"); exit(1); }
                    if ((int)ptr & 1) arc_rev_block_first = (arc_rev_block *) (ptr + 1);
                    else              arc_rev_block_first = (arc_rev_block *) ptr;
                    arc_rev_block_first -> start = ptr;
                    arc_rev_block_first -> current = & ( arc_rev_block_first -> arcs_rev[0] );
                    arc_rev_block_first -> next = next;
                    ab_rev = arc_rev_block_first;
                    rev_flag = true;
                }
                a_rev = &ab_rev->arcs_rev[0];
            }
            a_rev += k;
            i -> first_in = a_rev;
            ab_rev -> last_node = i;
        }
        /* i is the last node in block */
        i -> first_out = a_for;
        i -> first_in  = a_rev;
    }
    /* SECOND STAGE */
    for (ab_for=arc_for_block_first; ab_for; ab_for=ab_for->next)
    {
        ab_for -> current = ab_for -> last_node -> first_out;
    }
    for ( ab_for=ab_for_first, ab_rev=ab_rev_first;
          ab_for;
          ab_for=ab_for->next, ab_rev=ab_rev->next )
    for ( a_for=&ab_for->arcs_for[0], a_rev=&ab_rev->arcs_rev[0];
          a_for<&ab_for->arcs_for[ARC_BLOCK_SIZE];
          a_for++, a_rev++ )
    {
        arc_forward *af;
        arc_reverse *ar;
        node *from;
        int shift = 0, shift_new;
        captype r_cap = 0, r_rev_cap = 0, r_cap_new, r_rev_cap_new;
        if (!(from=(node *)(a_rev->sister))) continue;
        af = a_for;
        ar = a_rev;
        do
        {
            ar -> sister = NULL;
            shift_new = ((char *)(af->shift)) - (char *)from;
            r_cap_new = af -> r_cap;
            r_rev_cap_new = af -> r_rev_cap;
            if (shift)
            {
                af -> shift = shift;
                af -> r_cap = r_cap;
                af -> r_rev_cap = r_rev_cap;
            }
            shift = shift_new;
            r_cap = r_cap_new;
            r_rev_cap = r_rev_cap_new;
            af = -- from -> first_out;
            if ((arc_reverse *)(from->parent) != &a_rev_tmp)
            {
                from -> parent = (arc_forward *)(((arc_reverse *)(from -> parent)) - 1);
                ar = (arc_reverse *)(from -> parent);
            }
        } while ((from=(node *)(ar->sister)));
        af -> shift = shift;
        af -> r_cap = r_cap;
        af -> r_rev_cap = r_rev_cap;
    }
    for (ab_for=arc_for_block_first; ab_for; ab_for=ab_for->next)
    {
        i = ab_for -> last_node;
        a_for = i -> first_out;
        ab_for -> current -> shift     = a_for -> shift;
        ab_for -> current -> r_cap     = a_for -> r_cap;
        ab_for -> current -> r_rev_cap = a_for -> r_rev_cap;
        a_for -> shift = (int) (ab_for -> current + 1);
        i -> first_out = (arc_forward *) (((char *)a_for) - 1);
    }
    /* THIRD STAGE */
    for (ab_rev=arc_rev_block_first; ab_rev; ab_rev=ab_rev->next)
    {
        ab_rev -> current = ab_rev -> last_node -> first_in;
    }
    for (nb=node_block_first; nb; nb=nb->next)
    for (i=&nb->nodes[0]; i<nb->current; i++)
    {
        arc_forward *a_for_first, *a_for_last;
        a_for_first = i -> first_out;
        if (IS_ODD(a_for_first))
        {
            a_for_first = (arc_forward *) (((char *)a_for_first) + 1);
            a_for_last = (arc_forward *) ((a_for_first ++) -> shift);
        }
        else a_for_last = (i + 1) -> first_out;
        for (a_for=a_for_first; a_for<a_for_last; a_for++)
        {
            node *to = NEIGHBOR_NODE(i, a_for -> shift);
            a_rev = -- to -> first_in;
            a_rev -> sister = a_for;
        }
    }
    for (ab_rev=arc_rev_block_first; ab_rev; ab_rev=ab_rev->next)
    {
        i = ab_rev -> last_node;
        a_rev = i -> first_in;
        ab_rev -> current -> sister = a_rev -> sister;
        a_rev -> sister = (arc_forward *) (ab_rev -> current + 1);
        i -> first_in = (arc_reverse *) (((char *)a_rev) - 1);
    }
}


/***************************************************************************************************************************************/
void LinkedBlockList::addFront(ListType item)
{
    if ( m_head_block_size == GCLL_BLOCK_SIZE )
    {
        LLBlock *tmp      = (LLBlock *) new LLBlock;
        tmp -> m_next     = m_head;
        m_head            = tmp;
        m_head_block_size = 0;
    }
    m_head ->m_item[m_head_block_size] = item;
    m_head_block_size++;
}


ListType LinkedBlockList::next()
{
    ListType toReturn = m_cursor -> m_item[m_cursor_ind];
    m_cursor_ind++;
    if ( m_cursor == m_head && m_cursor_ind >= m_head_block_size )
    {
        m_cursor     = m_cursor ->m_next;
        m_cursor_ind = 0;
    }
    else if ( m_cursor_ind == GCLL_BLOCK_SIZE )
    {
        m_cursor = m_cursor ->m_next;
        m_cursor_ind = 0;
    }
    return(toReturn);
}


bool LinkedBlockList::hasNext()
{
    if ( m_cursor != 0 ) return (true);
    else return(false);
}


LinkedBlockList::~LinkedBlockList()
{
    LLBlock *tmp;
    while ( m_head != 0 ) 
    {
        tmp = m_head;
        m_head = m_head->m_next;
        delete tmp;
    }
};


/*void LinkedBlockList::clearList()
{
    LLBlock *tmp;
    while ( m_head != 0 ) 
    {
        tmp = m_head;
        m_head = m_head->m_next;
        delete tmp;
    }
};*/


/***************************************************************************************************************************************/
#define MAX_INTT 1000000000

void GCoptimization::initialize_memory()
{
    m_lookupPixVar = (PixelType *) new PixelType[m_nPixels];
    m_labelTable   = (LabelType *) new LabelType[m_nLabels];
    terminateOnError( !m_lookupPixVar || !m_labelTable ,"Not enough memory");
    for ( int i = 0; i < m_nLabels; i++ )
        m_labelTable[i] = i;
}


void GCoptimization::commonGridInitialization(PixelType width, PixelType height, int nLabels)
{
    terminateOnError( (width < 0) || (height <0) || (nLabels <0 ),"Illegal negative parameters");
    m_width       = width;
    m_height      = height;
    m_nPixels     = width*height;
    m_nLabels     = nLabels;
    m_grid_graph  = 1;
}


void GCoptimization::setParameters(int numParam, void *param)
{
    if (numParam != 1 ) 
        printf("\nInvalid number of parameters, can only set one parameter\\that is boolean label order\n");
    else
    {
        m_random_label_order = *((bool *) param);
    }
}


void GCoptimization::commonNonGridInitialization(PixelType nupixels, int num_labels)
{
    terminateOnError( (nupixels <0) || (num_labels <0 ),"Illegal negative parameters");
    m_nLabels        = num_labels;
    m_nPixels        = nupixels;
    m_grid_graph         = 0;
    m_neighbors = (LinkedBlockList *) new LinkedBlockList[nupixels];
    terminateOnError(!m_neighbors,"Not enough memory");
}


void GCoptimization::commonInitialization()
{
    m_needToFreeV        = 0;
    m_random_label_order = 1;
    initialize_memory();
}

/* Use this constructor only for grid graphs                                          */
GCoptimization::GCoptimization(PixelType width,PixelType height,LabelType nLabels,EnergyFunction *eng):MRFmodel(width,height,nLabels,eng)
{
    commonGridInitialization(width,height,nLabels);
    m_labeling           = (LabelType *) new LabelType[m_nPixels];
    terminateOnError(!m_labeling,"out of memory");
    for ( int i = 0; i < m_nPixels; i++ ) m_labeling[i] = (LabelType) 0;
    commonInitialization(); 
}

/* Use this constructor for general graphs                                            */
GCoptimization::GCoptimization(PixelType nPixels,int nLabels,EnergyFunction *eng):MRFmodel(nPixels,nLabels,eng) 
{
    commonNonGridInitialization(nPixels, nLabels);
    m_labeling           = (LabelType *) new LabelType[m_nPixels];
    terminateOnError(!m_labeling,"out of memory");
    for ( int i = 0; i < nPixels; i++ ) m_labeling[i] = (LabelType) 0;
    commonInitialization(); 
}


void GCoptimization::setData(EnergyTermType* dataArray)
{
    m_datacost  = dataArray;
}


void GCoptimization::setData(DataCostFn dataFn)
{
    m_dataFnPix = dataFn;
}


void GCoptimization::setSmoothness(EnergyTermType* V)
{
    m_smoothcost = V;
}


void GCoptimization::setSmoothness(int smoothExp,CostVal smoothMax, CostVal lambda)
{
    int i, j;
    CostVal cost;
    m_needToFreeV = 1;
    m_smoothcost = (CostVal *) new CostVal[m_nLabels*m_nLabels*sizeof(CostVal)];
    if (!m_smoothcost) { fprintf(stderr, "Not enough memory!\n"); exit(1); }
    for (i=0; i<m_nLabels; i++)
        for (j=i; j<m_nLabels; j++)
        {
            cost = (CostVal) ((smoothExp == 1) ? j - i : (j - i)*(j - i));
            if (cost > smoothMax) cost = smoothMax;
            m_smoothcost[i*m_nLabels + j] = m_smoothcost[j*m_nLabels + i] = cost*lambda;
        }
}


void GCoptimization::setCues(EnergyTermType* hCue, EnergyTermType* vCue)
{
    m_horizWeights    = hCue;
    m_vertWeights     = vCue;
}


void GCoptimization::setSmoothness(SmoothCostGeneralFn cost)
{
    m_smoothFnPix   = cost;
}


GCoptimization::EnergyType GCoptimization::dataEnergy()
{
    if ( m_dataType == ARRAY) 
        return(giveDataEnergyArray());
    else return(giveDataEnergyFnPix());
    return(0);
}

    
GCoptimization::EnergyType GCoptimization::giveDataEnergyArray()
{
    EnergyType eng = (EnergyType) 0;
    for ( int i = 0; i < m_nPixels; i++ )
        eng = eng + m_datacost(i,m_labeling[i])/m_nPixels;
    return(eng);
}

    
GCoptimization::EnergyType GCoptimization::giveDataEnergyFnPix()
{
    EnergyType eng = (EnergyType) 0;
    for ( int i = 0; i < m_nPixels; i++ )
        eng = eng + m_dataFnPix(i,m_labeling[i])/m_nPixels;
    return(eng);
}


GCoptimization::EnergyType GCoptimization::smoothnessEnergy()
{
    if ( m_grid_graph )
    {
        if ( m_smoothType != FUNCTION )
        {
			if(m_second_order) return(giveSmoothEnergy_G_ARRAY_SECOND_ORDER());
            else if(m_varWeights) return(giveSmoothEnergy_G_ARRAY_VW());
            else return(giveSmoothEnergy_G_ARRAY());
        }
        else return(giveSmoothEnergy_G_FnPix());
    }
    else
    {
        if ( m_smoothType != FUNCTION ) return(giveSmoothEnergy_NG_ARRAY());
        else  return(giveSmoothEnergy_NG_FnPix());
    }
    return(0);
}


GCoptimization::EnergyType GCoptimization::giveSmoothEnergy_NG_FnPix()
{
    EnergyType eng = (EnergyType) 0;
    int i;
    Neighbor *temp; 
    for ( i = 0; i < m_nPixels; i++ )
        if ( !m_neighbors[i].isEmpty() )
        {
            m_neighbors[i].setCursorFront();
            while ( m_neighbors[i].hasNext() )
            {
                temp = (Neighbor *) m_neighbors[i].next();
                if ( i < temp->to_node )
                    eng = eng + m_smoothFnPix(i,temp->to_node, m_labeling[i],m_labeling[temp->to_node])/m_nPixels;
            }
        }
    return(eng);
}


GCoptimization::EnergyType GCoptimization::giveSmoothEnergy_NG_ARRAY()
{
    EnergyType eng = (EnergyType) 0;
    int i;
    Neighbor *temp; 
    for ( i = 0; i < m_nPixels; i++ )
        if ( !m_neighbors[i].isEmpty() )
        {
            m_neighbors[i].setCursorFront();
            while ( m_neighbors[i].hasNext() )
            {
                temp = (Neighbor *) m_neighbors[i].next();

                if ( i < temp->to_node )
                    eng = eng + m_smoothcost(m_labeling[i],m_labeling[temp->to_node])*(temp->weight)/m_nPixels;
            }
        }
    return(eng);
}


GCoptimization::EnergyType GCoptimization::giveSmoothEnergy_G_ARRAY_VW()
{
    EnergyType eng = (EnergyType) 0;
    int x,y,pix;
    for ( y = 0; y < m_height; y++ )
        for ( x = 1; x < m_width; x++ )
        {
            pix = x+y*m_width;
            eng = eng + m_smoothcost(m_labeling[pix],m_labeling[pix-1])*m_horizWeights[pix-1]/m_nPixels;
        }
    for ( y = 1; y < m_height; y++ )
        for ( x = 0; x < m_width; x++ )
        {
            pix = x+y*m_width;
            eng = eng + m_smoothcost(m_labeling[pix],m_labeling[pix-m_width])*m_vertWeights[pix-m_width]/m_nPixels;
        }
    return(eng);
}


GCoptimization::EnergyType GCoptimization::giveSmoothEnergy_G_ARRAY()
{
    EnergyType eng = (EnergyType) 0;
    int x,y,pix;
    for ( y = 0; y < m_height; y++ )
        for ( x = 1; x < m_width; x++ )
        {
            pix = x+y*m_width;
            eng = eng + m_smoothcost(m_labeling[pix],m_labeling[pix-1])/m_nPixels;
        }
    for ( y = 1; y < m_height; y++ )
        for ( x = 0; x < m_width; x++ )
        {
            pix = x+y*m_width;
            eng = eng + m_smoothcost(m_labeling[pix],m_labeling[pix-m_width])/m_nPixels;
        }
    return(eng);
}


GCoptimization::EnergyType GCoptimization::giveSmoothEnergy_G_ARRAY_SECOND_ORDER()
{
    EnergyType eng = (EnergyType) 0;
    int x,y,pix;
    for ( y = 0; y < m_height; y++ )
        for ( x = 1; x < m_width; x++ )
        {
            pix = x+y*m_width;
            eng = eng + m_smoothcost(m_labeling[pix],m_labeling[pix-1])/m_nPixels;
        }
    for ( y = 1; y < m_height; y++ )
        for ( x = 0; x < m_width; x++ )
        {
            pix = x+y*m_width;
            eng = eng + m_smoothcost(m_labeling[pix],m_labeling[pix-m_width])/m_nPixels;
        }
    for ( y = 1; y < m_height; y++ )
        for ( x = 1; x < m_width; x++ )
        {
            pix = x+y*m_width;
            eng = eng + m_smoothcost(m_labeling[pix],m_labeling[pix-m_width-1])/m_nPixels;
        }
    for ( y = 1; y < m_height; y++ )
        for ( x = 0; x < m_width-1; x++ )
        {
            pix = x+y*m_width;
            eng = eng + m_smoothcost(m_labeling[pix],m_labeling[pix-m_width+1])/m_nPixels;
        }
    return(2*eng);
}


GCoptimization::EnergyType GCoptimization::giveSmoothEnergy_G_FnPix()
{
    EnergyType eng = (EnergyType) 0;
    int x,y,pix;
    for ( y = 0; y < m_height; y++ )
        for ( x = 1; x < m_width; x++ )
        {
            pix = x+y*m_width;
            eng = eng + m_smoothFnPix(pix,pix-1,m_labeling[pix],m_labeling[pix-1])/m_nPixels;
        }
    for ( y = 1; y < m_height; y++ )
        for ( x = 0; x < m_width; x++ )
        {
            pix = x+y*m_width;
            eng = eng + m_smoothFnPix(pix,pix-m_width,m_labeling[pix],m_labeling[pix-m_width])/m_nPixels;
        }
    return(eng);
}


void GCoptimization::setLabelOrder(bool RANDOM_LABEL_ORDER)
{
    m_random_label_order = RANDOM_LABEL_ORDER;
}


void GCoptimization::terminateOnError(bool error_condition,const char *message)
{ 
   if  (error_condition) 
   {
      printf("\n %s \n", message);
      exit(1);
    }
}


void GCoptimization::clearAnswer()
{
    if (!m_labeling ) {printf("First initialize algorithm" );exit(0);}
    memset(m_labeling, 0, m_nPixels*sizeof(Label));
}


void GCoptimization::scramble_label_table()
{
   LabelType r1,r2,temp;
   int num_times,cnt;
   num_times = m_nLabels*2;
   srand(clock());
   for ( cnt = 0; cnt < num_times; cnt++ )
   {
      r1 = rand()%m_nLabels;  
      r2 = rand()%m_nLabels;  
      temp             = m_labelTable[r1];
      m_labelTable[r1] = m_labelTable[r2];
      m_labelTable[r2] = temp;
   }
}


void GCoptimization::add_t_links_ARRAY(Energy *e,Energy::Var *variables,int size,LabelType alpha_label)
{
    for ( int i = 0; i < size; i++ )
        e -> add_term1(variables[i], m_datacost(m_lookupPixVar[i],alpha_label),
                                     m_datacost(m_lookupPixVar[i],m_labeling[m_lookupPixVar[i]]));

}


void GCoptimization::add_t_links_FnPix(Energy *e,Energy::Var *variables,int size,LabelType alpha_label)
{
    for ( int i = 0; i < size; i++ )
        e -> add_term1(variables[i], m_dataFnPix(m_lookupPixVar[i],alpha_label),
                                     m_dataFnPix(m_lookupPixVar[i],m_labeling[m_lookupPixVar[i]]));

}


void GCoptimization::setNeighbors(PixelType pixel1, int pixel2, EnergyTermType weight)
{
    assert(pixel1 < m_nPixels && pixel1 >= 0 && pixel2 < m_nPixels && pixel2 >= 0);
    assert(m_grid_graph == 0);
    Neighbor *temp1 = (Neighbor *) new Neighbor;
    Neighbor *temp2 = (Neighbor *) new Neighbor;
    temp1->weight  = weight;
    temp1->to_node = pixel2;
    temp2->weight  = weight;
    temp2->to_node = pixel1;
    m_neighbors[pixel1].addFront(temp1);
    m_neighbors[pixel2].addFront(temp2);
}


/*void GCoptimization::clearNeighbors()
{
	m_neighbors->clearList();
}*/


GCoptimization::~GCoptimization()
{
    delete [] m_labeling;
    if ( ! m_grid_graph ) delete [] m_neighbors;            
    delete [] m_labelTable;
    delete [] m_lookupPixVar;
    if (m_needToFreeV) delete [] m_smoothcost;
}


Swap::Swap(PixelType width,PixelType height,int num_labels, EnergyFunction *eng):GCoptimization(width,height,num_labels,eng)
{
    m_pixels = new PixelType[m_nPixels];
    terminateOnError( !m_pixels ,"Not enough memory");
}


Swap::Swap(PixelType nPixels, int num_labels, EnergyFunction *eng):GCoptimization(nPixels,num_labels,eng)
{
    m_pixels = new PixelType[m_nPixels];
    terminateOnError( !m_pixels ,"Not enough memory");
}


Swap::~Swap()
{
    delete [] m_pixels;
}


GCoptimization::EnergyType Swap::swap(int max_num_iterations)
{
    return(start_swap(max_num_iterations)); 
}


GCoptimization::EnergyType Swap::swap()
{
    return(start_swap(MAX_INTT));
}


GCoptimization::EnergyType Swap::start_swap(int max_num_iterations )
{
    int curr_cycle = 1;
    EnergyType new_energy,old_energy;
    new_energy = dataEnergy()+smoothnessEnergy();
    old_energy = new_energy+1;
    while ( old_energy > new_energy  && curr_cycle <= max_num_iterations)
    {
		printf("%f\n",new_energy);
        old_energy = new_energy;
        new_energy = oneSwapIteration();
        curr_cycle++;   
    }
    return(new_energy);
}


GCoptimization::EnergyType Swap::oneSwapIteration()
{
    int next,next1;
    if (m_random_label_order) scramble_label_table();
    for (next = 0;  next < m_nLabels;  next++ )
        for (next1 = m_nLabels - 1;  next1 >= 0;  next1-- )
            if ( m_labelTable[next] < m_labelTable[next1] )
            {
                perform_alpha_beta_swap(m_labelTable[next],m_labelTable[next1]);
            }
    return(dataEnergy()+smoothnessEnergy());
}


GCoptimization::EnergyType Swap::alpha_beta_swap(LabelType alpha_label, LabelType beta_label)
{
    terminateOnError( alpha_label < 0 || alpha_label >= m_nLabels || beta_label < 0 || beta_label >= m_nLabels,
        "Illegal Label to Expand On");
    perform_alpha_beta_swap(alpha_label,beta_label);
    return(dataEnergy()+smoothnessEnergy());
}


void Swap::add_t_links_ARRAY_swap(Energy *e,Energy::Var *variables,int size, LabelType alpha_label, LabelType beta_label, PixelType *pixels)
{
    for ( int i = 0; i < size; i++ )
        e -> add_term1(variables[i], m_datacost(pixels[i],alpha_label),
                                     m_datacost(pixels[i],beta_label));
}
    

void Swap::add_t_links_FnPix_swap(Energy *e,Energy::Var *variables,int size, LabelType alpha_label, LabelType beta_label, PixelType *pixels)
{
    for ( int i = 0; i < size; i++ )
        e -> add_term1(variables[i], m_dataFnPix(pixels[i],alpha_label),
                                     m_dataFnPix(pixels[i],beta_label));
}


void Swap::perform_alpha_beta_swap(LabelType alpha_label, LabelType beta_label)
{
    PixelType i,size = 0;
    Energy *e = new Energy();
    for ( i = 0; i < m_nPixels; i++ )
    {
        if ( m_labeling[i] == alpha_label || m_labeling[i] == beta_label)
        {
            m_pixels[size]    = i;
            m_lookupPixVar[i] = size;
            size++;
        }
    }
    if ( size == 0 ) return;
    Energy::Var *variables = (Energy::Var *) new Energy::Var[size];
    if (!variables) { fprintf(stderr, "Not enough memory!\n"); exit(1); }
    for ( i = 0; i < size; i++ )
        variables[i] = e ->add_variable();
    if ( m_dataType == ARRAY ) add_t_links_ARRAY_swap(e,variables,size,alpha_label,beta_label,m_pixels);
    else  add_t_links_FnPix_swap(e,variables,size,alpha_label,beta_label,m_pixels);
    if ( m_grid_graph )
    {
        if ( m_smoothType != FUNCTION ) //GAB
        {
			if(m_second_order) set_up_swap_energy_G_ARRAY_SECOND_ORDER(size,alpha_label,beta_label,m_pixels,e,variables);
            else if (m_varWeights) set_up_swap_energy_G_ARRAY_VW(size,alpha_label,beta_label,m_pixels,e,variables);
            else set_up_swap_energy_G_ARRAY(size,alpha_label,beta_label,m_pixels,e,variables);
        }
        else  set_up_swap_energy_G_FnPix(size,alpha_label,beta_label,m_pixels,e,variables);
    }
    else
    {
        if ( m_smoothType  != FUNCTION  ) set_up_swap_energy_NG_ARRAY(size,alpha_label,beta_label,m_pixels,e,variables);
        else set_up_swap_energy_NG_FnPix(size,alpha_label,beta_label,m_pixels,e,variables);
    }
    e -> minimize();
    for ( i = 0; i < size; i++ )
        if ( e->get_var(variables[i]) == 0 )
            m_labeling[m_pixels[i]] = alpha_label;
        else m_labeling[m_pixels[i]] = beta_label;
    delete [] variables;
    delete e;
}


void Swap::set_up_swap_energy_NG_ARRAY(int size,LabelType alpha_label,LabelType beta_label, PixelType *pixels,Energy* e, Energy::Var *variables)
{
    PixelType nPix,pix,i;
    EnergyTermType weight;
    Neighbor *tmp;
    for ( i = 0; i < size; i++ )
    {
        pix = pixels[i];
        if ( !m_neighbors[pix].isEmpty() )
        {
            m_neighbors[pix].setCursorFront();
            while ( m_neighbors[pix].hasNext() )
            {
                tmp = (Neighbor *) (m_neighbors[pix].next());
                nPix   = tmp->to_node;
                weight = tmp->weight;
                
                if ( m_labeling[nPix] == alpha_label || m_labeling[nPix] == beta_label)
                {
                    if ( pix < nPix )
                        e ->add_term2(variables[i],variables[m_lookupPixVar[nPix]],
                                      m_smoothcost(alpha_label,alpha_label)*weight,
                                      m_smoothcost(alpha_label,beta_label)*weight,
                                      m_smoothcost(beta_label,alpha_label)*weight,
                                      m_smoothcost(beta_label,beta_label)*weight);
                }
                else
                    e ->add_term1(variables[i],m_smoothcost(alpha_label,m_labeling[nPix])*weight,
                                               m_smoothcost(beta_label,m_labeling[nPix])*weight);
            }
        }
    }
}


void Swap::set_up_swap_energy_NG_FnPix(int size,LabelType alpha_label,LabelType beta_label, PixelType *pixels,Energy* e, Energy::Var *variables)
{
    PixelType nPix,pix,i;
    Neighbor *tmp;
    for ( i = 0; i < size; i++ )
    {
        pix = pixels[i];
        if ( !m_neighbors[pix].isEmpty() )
        {
            m_neighbors[pix].setCursorFront();
            
            while ( m_neighbors[pix].hasNext() )
            {
                tmp = (Neighbor *) (m_neighbors[pix].next());
                nPix   = tmp->to_node;
                
                if ( m_labeling[nPix] == alpha_label || m_labeling[nPix] == beta_label)
                {
                    if ( pix < nPix )
                        e ->add_term2(variables[i],variables[m_lookupPixVar[nPix]],
                                      m_smoothFnPix(pix,nPix,alpha_label,alpha_label),
                                      m_smoothFnPix(pix,nPix,alpha_label,beta_label),
                                      m_smoothFnPix(pix,nPix,beta_label,alpha_label),
                                      m_smoothFnPix(pix,nPix,beta_label,beta_label) );
                }
                else
                    e ->add_term1(variables[i],m_smoothFnPix(pix,nPix,alpha_label,m_labeling[nPix]),
                                               m_smoothFnPix(pix,nPix,beta_label,m_labeling[nPix]));
            }
        }
    }
}


void Swap::set_up_swap_energy_G_FnPix(int size,LabelType alpha_label,LabelType beta_label, PixelType *pixels,Energy* e, Energy::Var *variables)
{
    PixelType nPix,pix,i,x,y;
    for ( i = 0; i < size; i++ )
    {
        pix = pixels[i];
        y = pix/m_width;
        x = pix - y*m_width;
        if ( x > 0 )
        {
            nPix = pix - 1;
    
            if ( m_labeling[nPix] == alpha_label || m_labeling[nPix] == beta_label)
                e ->add_term2(variables[i],variables[m_lookupPixVar[nPix]],
                              m_smoothFnPix(pix,nPix,alpha_label,alpha_label),
                              m_smoothFnPix(pix,nPix,alpha_label,beta_label),
                              m_smoothFnPix(pix,nPix,beta_label,alpha_label),
                              m_smoothFnPix(pix,nPix,beta_label,beta_label) );
                else
                    e ->add_term1(variables[i],m_smoothFnPix(pix,nPix,alpha_label,m_labeling[nPix]),
                                           m_smoothFnPix(pix,nPix,beta_label,m_labeling[nPix]));
        }   
        if ( y > 0 )
        {
            nPix = pix - m_width;
            if ( m_labeling[nPix] == alpha_label || m_labeling[nPix] == beta_label)
                e ->add_term2(variables[i],variables[m_lookupPixVar[nPix]],
                              m_smoothFnPix(pix,nPix,alpha_label,alpha_label),
                              m_smoothFnPix(pix,nPix,alpha_label,beta_label),
                              m_smoothFnPix(pix,nPix,beta_label,alpha_label),
                              m_smoothFnPix(pix,nPix,beta_label,beta_label) );
    
                else
                    e ->add_term1(variables[i],m_smoothFnPix(pix,nPix,alpha_label,m_labeling[nPix]),
                                           m_smoothFnPix(pix,nPix,beta_label,m_labeling[nPix]));
        }   
        if ( x < m_width - 1 )
        {
            nPix = pix + 1;
    
            if ( !(m_labeling[nPix] == alpha_label || m_labeling[nPix] == beta_label) )
                    e ->add_term1(variables[i],m_smoothFnPix(pix,nPix,alpha_label,m_labeling[nPix]),
                                               m_smoothFnPix(pix,nPix,beta_label,m_labeling[nPix]));
        }   
        if ( y < m_height - 1 )
        {
            nPix = pix + m_width;
    
            if ( !(m_labeling[nPix] == alpha_label || m_labeling[nPix] == beta_label) )
                e ->add_term1(variables[i],m_smoothFnPix(pix,nPix,alpha_label,m_labeling[nPix]),
                                           m_smoothFnPix(pix,nPix,beta_label,m_labeling[nPix]));
        }
    }
}



void Swap::set_up_swap_energy_G_ARRAY_VW(int size,LabelType alpha_label,LabelType beta_label, PixelType *pixels,Energy* e, Energy::Var *variables)
{
    PixelType nPix,pix,i,x,y;
    EnergyTermType weight;  
    for ( i = 0; i < size; i++ )
    {
        pix = pixels[i];
        y = pix/m_width;
        x = pix - y*m_width;
        if ( x > 0 )
        {
            nPix = pix - 1;
            weight = m_horizWeights[nPix];
    
            if ( m_labeling[nPix] == alpha_label || m_labeling[nPix] == beta_label)
                e ->add_term2(variables[i],variables[m_lookupPixVar[nPix]],
                              m_smoothcost(alpha_label,alpha_label)*weight,
                              m_smoothcost(alpha_label,beta_label)*weight,
                              m_smoothcost(beta_label,alpha_label)*weight,
                              m_smoothcost(beta_label,beta_label)*weight );
                else
                    e ->add_term1(variables[i],m_smoothcost(alpha_label,m_labeling[nPix])*weight,
                                           m_smoothcost(beta_label,m_labeling[nPix])*weight);
        }   
        if ( y > 0 )
        {
            nPix = pix - m_width;
            weight = m_vertWeights[nPix];

            if ( m_labeling[nPix] == alpha_label || m_labeling[nPix] == beta_label)
                e ->add_term2(variables[i],variables[m_lookupPixVar[nPix]],
                              m_smoothcost(alpha_label,alpha_label)*weight,
                              m_smoothcost(alpha_label,beta_label)*weight,
                              m_smoothcost(beta_label,alpha_label)*weight,
                              m_smoothcost(beta_label,beta_label)*weight );
                else
                    e ->add_term1(variables[i],m_smoothcost(alpha_label,m_labeling[nPix])*weight,
                                           m_smoothcost(beta_label,m_labeling[nPix])*weight);
        }   
        if ( x < m_width - 1 )
        {
            nPix = pix + 1;
    
            if ( !(m_labeling[nPix] == alpha_label || m_labeling[nPix] == beta_label) )
                    e ->add_term1(variables[i],m_smoothcost(alpha_label,m_labeling[nPix])*m_horizWeights[pix],
                                               m_smoothcost(beta_label,m_labeling[nPix])*m_horizWeights[pix]);
        }   
        if ( y < m_height - 1 )
        {
            nPix = pix + m_width;
    
            if ( !(m_labeling[nPix] == alpha_label || m_labeling[nPix] == beta_label) )
                e ->add_term1(variables[i],m_smoothcost(alpha_label,m_labeling[nPix])*m_vertWeights[pix],
                                           m_smoothcost(beta_label,m_labeling[nPix])*m_vertWeights[pix]);
        }
    }
}


void Swap::set_up_swap_energy_G_ARRAY(int size,LabelType alpha_label,LabelType beta_label, PixelType *pixels,Energy* e, Energy::Var *variables)
{
    PixelType nPix,pix,i,x,y;
    for ( i = 0; i < size; i++ )
    {
        pix = pixels[i];
        y = pix/m_width;
        x = pix - y*m_width;
        if ( x > 0 )
        {
            nPix = pix - 1;
            if ( m_labeling[nPix] == alpha_label || m_labeling[nPix] == beta_label)
                e ->add_term2(variables[i],variables[m_lookupPixVar[nPix]],
                              m_smoothcost(alpha_label,alpha_label),
                              m_smoothcost(alpha_label,beta_label),
                              m_smoothcost(beta_label,alpha_label),
                              m_smoothcost(beta_label,beta_label) );
                else
                    e ->add_term1(variables[i],m_smoothcost(alpha_label,m_labeling[nPix]),
                                           m_smoothcost(beta_label,m_labeling[nPix]));
        }   
        if ( y > 0 )
        {
            nPix = pix - m_width;
            if ( m_labeling[nPix] == alpha_label || m_labeling[nPix] == beta_label)
                e ->add_term2(variables[i],variables[m_lookupPixVar[nPix]],
                              m_smoothcost(alpha_label,alpha_label),
                              m_smoothcost(alpha_label,beta_label),
                              m_smoothcost(beta_label,alpha_label),
                              m_smoothcost(beta_label,beta_label) );
                else
                    e ->add_term1(variables[i],m_smoothcost(alpha_label,m_labeling[nPix]),
                                           m_smoothcost(beta_label,m_labeling[nPix]));
        }   
        if ( x < m_width - 1 )
        {
            nPix = pix + 1;
            if ( !(m_labeling[nPix] == alpha_label || m_labeling[nPix] == beta_label) )
                    e ->add_term1(variables[i],m_smoothcost(alpha_label,m_labeling[nPix]),
                                               m_smoothcost(beta_label,m_labeling[nPix]));
        }   
        if ( y < m_height - 1 )
        {
            nPix = pix + m_width;
    
            if ( !(m_labeling[nPix] == alpha_label || m_labeling[nPix] == beta_label))
                e ->add_term1(variables[i],m_smoothcost(alpha_label,m_labeling[nPix]),
                                           m_smoothcost(beta_label,m_labeling[nPix]));
        }
    }
}


void Swap::set_up_swap_energy_G_ARRAY_SECOND_ORDER(int size,LabelType alpha_label,LabelType beta_label, PixelType *pixels,Energy* e, Energy::Var *variables)
{
    PixelType nPix,pix,i,x,y;
    for ( i = 0; i < size; i++ )
    {
        pix = pixels[i];
        y = pix/m_width;
        x = pix - y*m_width;
        if ( x > 0 )
        {
            nPix = pix - 1;
            if ( m_labeling[nPix] == alpha_label || m_labeling[nPix] == beta_label)
                e ->add_term2(variables[i],variables[m_lookupPixVar[nPix]],
                              m_smoothcost(alpha_label,alpha_label),
                              m_smoothcost(alpha_label,beta_label),
                              m_smoothcost(beta_label,alpha_label),
                              m_smoothcost(beta_label,beta_label) );
                else
                    e ->add_term1(variables[i],m_smoothcost(alpha_label,m_labeling[nPix]),m_smoothcost(beta_label,m_labeling[nPix]));
        }   
        if ( y > 0 )
        {
            nPix = pix - m_width;
            if ( m_labeling[nPix] == alpha_label || m_labeling[nPix] == beta_label)
                e ->add_term2(variables[i],variables[m_lookupPixVar[nPix]],
                              m_smoothcost(alpha_label,alpha_label),
                              m_smoothcost(alpha_label,beta_label),
                              m_smoothcost(beta_label,alpha_label),
                              m_smoothcost(beta_label,beta_label) );
                else
                    e ->add_term1(variables[i],m_smoothcost(alpha_label,m_labeling[nPix]),m_smoothcost(beta_label,m_labeling[nPix]));
        }   
        if ( x > 0 && y > 0 )
        {
            nPix = pix - m_width - 1;
            if ( m_labeling[nPix] == alpha_label || m_labeling[nPix] == beta_label)
                e ->add_term2(variables[i],variables[m_lookupPixVar[nPix]],
                              m_smoothcost(alpha_label,alpha_label),
                              m_smoothcost(alpha_label,beta_label),
                              m_smoothcost(beta_label,alpha_label),
                              m_smoothcost(beta_label,beta_label) );
                else
                    e ->add_term1(variables[i],m_smoothcost(alpha_label,m_labeling[nPix]),m_smoothcost(beta_label,m_labeling[nPix]));
        }
        if ( x < m_width - 1 && y > 0 )
        {
            nPix = pix - m_width + 1;
            if ( m_labeling[nPix] == alpha_label || m_labeling[nPix] == beta_label)
                e ->add_term2(variables[i],variables[m_lookupPixVar[nPix]],
                              m_smoothcost(alpha_label,alpha_label),
                              m_smoothcost(alpha_label,beta_label),
                              m_smoothcost(beta_label,alpha_label),
                              m_smoothcost(beta_label,beta_label) );
                else
                    e ->add_term1(variables[i],m_smoothcost(alpha_label,m_labeling[nPix]),m_smoothcost(beta_label,m_labeling[nPix]));
        }
        if ( x < m_width - 1 )
        {
            nPix = pix + 1;
            if ( !(m_labeling[nPix] == alpha_label || m_labeling[nPix] == beta_label) )
                    e ->add_term1(variables[i],m_smoothcost(alpha_label,m_labeling[nPix]),m_smoothcost(beta_label,m_labeling[nPix]));
        }   
        if ( y < m_height - 1 )
        {
            nPix = pix + m_width;
            if ( !(m_labeling[nPix] == alpha_label || m_labeling[nPix] == beta_label))
                e ->add_term1(variables[i],m_smoothcost(alpha_label,m_labeling[nPix]),m_smoothcost(beta_label,m_labeling[nPix]));
        }
        if ( x < m_width - 1 && y < m_height - 1 )
        {
            nPix = pix + m_width + 1;
            if ( !(m_labeling[nPix] == alpha_label || m_labeling[nPix] == beta_label))
                e ->add_term1(variables[i],m_smoothcost(alpha_label,m_labeling[nPix]),m_smoothcost(beta_label,m_labeling[nPix]));
        }
        if ( x > 0 && y < m_height - 1 )
        {
            nPix = pix + m_width - 1;
            if ( !(m_labeling[nPix] == alpha_label || m_labeling[nPix] == beta_label))
                e ->add_term1(variables[i],m_smoothcost(alpha_label,m_labeling[nPix]),m_smoothcost(beta_label,m_labeling[nPix]));
        }
    }
}


void Swap::optimizeAlg(int nIterations)
{
    swap(nIterations);
}


GCoptimization::EnergyType Expansion::expansion(int max_num_iterations)
{
    return(start_expansion(max_num_iterations)); 
}


GCoptimization::EnergyType Expansion::expansion()
{
    return(start_expansion(MAX_INTT));
}


GCoptimization::EnergyType Expansion::start_expansion(int max_num_iterations )
{
    int curr_cycle = 1;
    EnergyType new_energy,old_energy;
    new_energy = dataEnergy()+smoothnessEnergy();
    old_energy = new_energy+1;
    while ( old_energy > new_energy  && curr_cycle <= max_num_iterations)
    {
		printf("%f\n",new_energy);
        old_energy = new_energy;
        new_energy = oneExpansionIteration();
        curr_cycle++;   
    }
    return(new_energy);
}


GCoptimization::EnergyType Expansion::alpha_expansion(LabelType label)
{
    terminateOnError( label < 0 || label >= m_nLabels,"Illegal Label to Expand On");
    perform_alpha_expansion(label);
    return(dataEnergy()+smoothnessEnergy());
}


void Expansion::perform_alpha_expansion(LabelType alpha_label)
{
    PixelType i,size = 0; 
    Energy *e = new Energy();
    for ( i = 0; i < m_nPixels; i++ )
    {
        if ( m_labeling[i] != alpha_label )
        {
            m_lookupPixVar[size] = i;
            size++;
        }
    }
    if ( size > 0 ) 
    {
        Energy::Var *variables = (Energy::Var *) new Energy::Var[size];
        if ( !variables) {printf("\nOut of memory, exiting");exit(1);}
        for ( i = 0; i < size; i++ )
            variables[i] = e ->add_variable();
        if ( m_dataType == ARRAY ) add_t_links_ARRAY(e,variables,size,alpha_label);
        else  add_t_links_FnPix(e,variables,size,alpha_label);
        if ( m_grid_graph )
        {
            if ( m_smoothType != FUNCTION )
            {
				if(m_second_order) set_up_expansion_energy_G_ARRAY_SECOND_ORDER(size,alpha_label,e,variables);
                else if (m_varWeights) set_up_expansion_energy_G_ARRAY_VW(size,alpha_label,e,variables);
                else set_up_expansion_energy_G_ARRAY(size,alpha_label,e,variables);
            }
            else set_up_expansion_energy_G_FnPix(size,alpha_label,e,variables);
        }
        else
        {
            if ( m_smoothType != FUNCTION ) set_up_expansion_energy_NG_ARRAY(size,alpha_label,e,variables);
            else if ( m_smoothType == FUNCTION) set_up_expansion_energy_NG_FnPix(size,alpha_label,e,variables);
        }
        e -> minimize();
        for ( i = 0,size = 0; i < m_nPixels; i++ )
        {
            if ( m_labeling[i] != alpha_label )
            {
                if ( e->get_var(variables[size]) == 0 )
                    m_labeling[i] = alpha_label;

                size++;
            }
        }
        delete [] variables;
    }
    delete e;
}


/* Performs alpha-expansion for non regular grid graph for case when energy terms are NOT     */
/* specified by a function */
void Expansion::set_up_expansion_energy_NG_ARRAY(int size, LabelType alpha_label,Energy *e,Energy::Var *variables )
{
    EnergyTermType weight;
    Neighbor *tmp;
    int i,nPix,pix;;
    for ( i = size - 1; i >= 0; i-- )
    {
        pix = m_lookupPixVar[i];
        m_lookupPixVar[pix] = i;
        if ( !m_neighbors[pix].isEmpty() )
        {
            m_neighbors[pix].setCursorFront();
            while ( m_neighbors[pix].hasNext() )
            {
                tmp = (Neighbor *) (m_neighbors[pix].next());
                nPix   = tmp->to_node;
                weight = tmp->weight;
                if ( m_labeling[nPix] != alpha_label )
                {
                    if ( pix < nPix )
                        e ->add_term2(variables[i],variables[m_lookupPixVar[nPix]],
                                      m_smoothcost(alpha_label,alpha_label)*weight,
                                      m_smoothcost(alpha_label,m_labeling[nPix])*weight,
                                      m_smoothcost(m_labeling[pix],alpha_label)*weight,
                                      m_smoothcost(m_labeling[pix],m_labeling[nPix])*weight);
                }
                else
                    e ->add_term1(variables[i],m_smoothcost(alpha_label,alpha_label)*weight,
                                           m_smoothcost(m_labeling[pix],alpha_label)*weight);
            }
        }
    }
}


/* Performs alpha-expansion for non regular grid graph for case when energy terms are        */
/* specified by a function */
void Expansion::set_up_expansion_energy_NG_FnPix(int size, LabelType alpha_label,Energy *e,Energy::Var *variables )
{
    Neighbor *tmp;
    int i,nPix,pix;
    for ( i = size - 1; i >= 0; i-- )
    {
        pix = m_lookupPixVar[i];
        m_lookupPixVar[pix] = i;
        if ( !m_neighbors[pix].isEmpty() )
        {
            m_neighbors[pix].setCursorFront();
            while ( m_neighbors[pix].hasNext() )
            {
                tmp = (Neighbor *) (m_neighbors[pix].next());
                nPix   = tmp->to_node;
                if ( m_labeling[nPix] != alpha_label )
                {
                    if ( pix < nPix )
                        e ->add_term2(variables[i],variables[m_lookupPixVar[nPix]],
                                      m_smoothFnPix(pix,nPix,alpha_label,alpha_label),
                                      m_smoothFnPix(pix,nPix,alpha_label,m_labeling[nPix]),
                                      m_smoothFnPix(pix,nPix,m_labeling[pix],alpha_label),
                                      m_smoothFnPix(pix,nPix,m_labeling[pix],m_labeling[nPix]));
                }
                else
                    e ->add_term1(variables[i],m_smoothFnPix(pix,nPix,alpha_label,m_labeling[nPix]),
                                               m_smoothFnPix(pix,nPix,m_labeling[pix],alpha_label));
            }
        }
    }
}


/* Performs alpha-expansion for  regular grid graph for case when energy terms are NOT        */
/* specified by a function */
void Expansion::set_up_expansion_energy_G_ARRAY_VW(int size, LabelType alpha_label,Energy *e, Energy::Var *variables )
{
    int i,nPix,pix,x,y;
    EnergyTermType weight;
    for ( i = size - 1; i >= 0; i-- )
    {
        pix = m_lookupPixVar[i];
        y = pix/m_width;
        x = pix - y*m_width;
        m_lookupPixVar[pix] = i;
        if ( x < m_width - 1 )
        {
            nPix = pix + 1;
            weight = m_horizWeights[pix];
            if ( m_labeling[nPix] != alpha_label )
                e ->add_term2(variables[i],variables[m_lookupPixVar[nPix]],
                              m_smoothcost(alpha_label,alpha_label)*weight,
                              m_smoothcost(alpha_label,m_labeling[nPix])*weight,
                              m_smoothcost(m_labeling[pix],alpha_label)*weight,
                              m_smoothcost(m_labeling[pix],m_labeling[nPix])*weight);
            else   e ->add_term1(variables[i],m_smoothcost(alpha_label,m_labeling[nPix])*weight,
                                 m_smoothcost(m_labeling[pix],alpha_label)*weight);
        }   
        if ( y < m_height - 1 )
        {
            nPix = pix + m_width;
            weight = m_vertWeights[pix];
            if ( m_labeling[nPix] != alpha_label )
                e ->add_term2(variables[i],variables[m_lookupPixVar[nPix]],
                              m_smoothcost(alpha_label,alpha_label)*weight ,
                              m_smoothcost(alpha_label,m_labeling[nPix])*weight,
                              m_smoothcost(m_labeling[pix],alpha_label)*weight ,
                              m_smoothcost(m_labeling[pix],m_labeling[nPix])*weight );
            else   e ->add_term1(variables[i],m_smoothcost(alpha_label,m_labeling[nPix])*weight,
                                 m_smoothcost(m_labeling[pix],alpha_label)*weight);
        }   
        if ( x > 0 )
        {
            nPix = pix - 1;
    
            if ( m_labeling[nPix] == alpha_label )
               e ->add_term1(variables[i],m_smoothcost(alpha_label,m_labeling[nPix])*m_horizWeights[nPix],
                                 m_smoothcost(m_labeling[pix],alpha_label)*m_horizWeights[nPix]);
        }   
        if ( y > 0 )
        {
            nPix = pix - m_width;
            if ( m_labeling[nPix] == alpha_label )
               e ->add_term1(variables[i],m_smoothcost(alpha_label,alpha_label)*m_vertWeights[nPix],
                                 m_smoothcost(m_labeling[pix],alpha_label)*m_vertWeights[nPix]);
        }   
    }
}


/* Performs alpha-expansion for  regular grid graph for case when energy terms are NOT        */
/* specified by a function */
void Expansion::set_up_expansion_energy_G_ARRAY(int size, LabelType alpha_label,Energy *e, Energy::Var *variables )
{
    int i,nPix,pix,x,y;
    for ( i = size - 1; i >= 0; i-- )
    {
        pix = m_lookupPixVar[i];
        y = pix/m_width;
        x = pix - y*m_width;
        m_lookupPixVar[pix] = i;
        if ( x < m_width - 1 )
        {
            nPix = pix + 1;
            if ( m_labeling[nPix] != alpha_label )
                e ->add_term2(variables[i],variables[m_lookupPixVar[nPix]],
                              m_smoothcost(alpha_label,alpha_label),
                              m_smoothcost(alpha_label,m_labeling[nPix]),
                              m_smoothcost(m_labeling[pix],alpha_label),
                              m_smoothcost(m_labeling[pix],m_labeling[nPix]));
            else   e ->add_term1(variables[i],m_smoothcost(alpha_label,m_labeling[nPix]),
                                 m_smoothcost(m_labeling[pix],alpha_label));
        }   
        if ( y < m_height - 1 )
        {
            nPix = pix + m_width;
            if ( m_labeling[nPix] != alpha_label )
                e ->add_term2(variables[i],variables[m_lookupPixVar[nPix]],
                              m_smoothcost(alpha_label,alpha_label) ,
                              m_smoothcost(alpha_label,m_labeling[nPix]),
                              m_smoothcost(m_labeling[pix],alpha_label) ,
                              m_smoothcost(m_labeling[pix],m_labeling[nPix]) );
            else   e ->add_term1(variables[i],m_smoothcost(alpha_label,m_labeling[nPix]),
                                 m_smoothcost(m_labeling[pix],alpha_label));
        }   
        if ( x > 0 )
        {
            nPix = pix - 1;
            if ( m_labeling[nPix] == alpha_label )
               e ->add_term1(variables[i],m_smoothcost(alpha_label,m_labeling[nPix]),
                                 m_smoothcost(m_labeling[pix],alpha_label) );
        }   
        if ( y > 0 )
        {
            nPix = pix - m_width;
            if ( m_labeling[nPix] == alpha_label )
               e ->add_term1(variables[i],m_smoothcost(alpha_label,alpha_label),
                                 m_smoothcost(m_labeling[pix],alpha_label));
        }   
    }
}


void Expansion::set_up_expansion_energy_G_ARRAY_SECOND_ORDER(int size, LabelType alpha_label,Energy *e, Energy::Var *variables )
{
    int i,nPix,pix,x,y;
    for ( i = size - 1; i >= 0; i-- )
    {
        pix = m_lookupPixVar[i];
        y = pix/m_width;
        x = pix - y*m_width;
        m_lookupPixVar[pix] = i;
        if ( x < m_width - 1 )
        {
            nPix = pix + 1;
            if ( m_labeling[nPix] != alpha_label )
                e ->add_term2(variables[i],variables[m_lookupPixVar[nPix]],
                              m_smoothcost(alpha_label,alpha_label),
                              m_smoothcost(alpha_label,m_labeling[nPix]),
                              m_smoothcost(m_labeling[pix],alpha_label),
                              m_smoothcost(m_labeling[pix],m_labeling[nPix]));
            else   e ->add_term1(variables[i],m_smoothcost(alpha_label,m_labeling[nPix]),m_smoothcost(m_labeling[pix],alpha_label));
        }   
        if ( y < m_height - 1 )
        {
            nPix = pix + m_width;
            if ( m_labeling[nPix] != alpha_label )
                e ->add_term2(variables[i],variables[m_lookupPixVar[nPix]],
                              m_smoothcost(alpha_label,alpha_label) ,
                              m_smoothcost(alpha_label,m_labeling[nPix]),
                              m_smoothcost(m_labeling[pix],alpha_label) ,
                              m_smoothcost(m_labeling[pix],m_labeling[nPix]) );
            else   e ->add_term1(variables[i],m_smoothcost(alpha_label,m_labeling[nPix]),m_smoothcost(m_labeling[pix],alpha_label));
        }
        if ( x < m_width - 1 && y < m_height - 1 )
        {
            nPix = pix + m_width + 1;
            if ( m_labeling[nPix] != alpha_label )
                e ->add_term2(variables[i],variables[m_lookupPixVar[nPix]],
                              m_smoothcost(alpha_label,alpha_label) ,
                              m_smoothcost(alpha_label,m_labeling[nPix]),
                              m_smoothcost(m_labeling[pix],alpha_label) ,
                              m_smoothcost(m_labeling[pix],m_labeling[nPix]) );
            else   e ->add_term1(variables[i],m_smoothcost(alpha_label,m_labeling[nPix]),m_smoothcost(m_labeling[pix],alpha_label));
        }
        if ( x > 0 && y < m_height - 1 )
        {
            nPix = pix + m_width - 1;
            if ( m_labeling[nPix] != alpha_label )
                e ->add_term2(variables[i],variables[m_lookupPixVar[nPix]],
                              m_smoothcost(alpha_label,alpha_label) ,
                              m_smoothcost(alpha_label,m_labeling[nPix]),
                              m_smoothcost(m_labeling[pix],alpha_label) ,
                              m_smoothcost(m_labeling[pix],m_labeling[nPix]) );
            else   e ->add_term1(variables[i],m_smoothcost(alpha_label,m_labeling[nPix]),m_smoothcost(m_labeling[pix],alpha_label));
        }
        if ( x > 0 )
        {
            nPix = pix - 1;
            if ( m_labeling[nPix] == alpha_label )
               e ->add_term1(variables[i],m_smoothcost(alpha_label,m_labeling[nPix]),m_smoothcost(m_labeling[pix],alpha_label) );
        }   
        if ( y > 0 )
        {
            nPix = pix - m_width;
            if ( m_labeling[nPix] == alpha_label )
               e ->add_term1(variables[i],m_smoothcost(alpha_label,alpha_label),m_smoothcost(m_labeling[pix],alpha_label));
        }
        if ( x > 0 && y > 0 )
        {
            nPix = pix - m_width - 1;
            if ( m_labeling[nPix] == alpha_label )
               e ->add_term1(variables[i],m_smoothcost(alpha_label,alpha_label),m_smoothcost(m_labeling[pix],alpha_label));
        }
        if ( x < m_width - 1 && y > 0 )
        {
            nPix = pix - m_width + 1;
            if ( m_labeling[nPix] == alpha_label )
               e ->add_term1(variables[i],m_smoothcost(alpha_label,alpha_label),m_smoothcost(m_labeling[pix],alpha_label));
        }
    }
}


/* Performs alpha-expansion for  regular grid graph for case when energy terms are NOT        */
/* specified by a function */
void Expansion::set_up_expansion_energy_G_FnPix(int size, LabelType alpha_label,Energy *e, Energy::Var *variables )
{
    int i,nPix,pix,x,y;
    for ( i = size - 1; i >= 0; i-- )
    {
        pix = m_lookupPixVar[i];
        y = pix/m_width;
        x = pix - y*m_width;
        m_lookupPixVar[pix] = i;
        if ( x < m_width - 1 )
        {
            nPix = pix + 1;
            if ( m_labeling[nPix] != alpha_label )
                e ->add_term2(variables[i],variables[m_lookupPixVar[nPix]],
                              m_smoothFnPix(pix,nPix,alpha_label,alpha_label),
                              m_smoothFnPix(pix,nPix,alpha_label,m_labeling[nPix]),
                              m_smoothFnPix(pix,nPix,m_labeling[pix],alpha_label),
                              m_smoothFnPix(pix,nPix,m_labeling[pix],m_labeling[nPix]));
            else   e ->add_term1(variables[i],m_smoothFnPix(pix,nPix,alpha_label,m_labeling[nPix]),
                                 m_smoothFnPix(pix,nPix,m_labeling[pix],alpha_label));
        }   
        if ( y < m_height - 1 )
        {
            nPix = pix + m_width;
            if ( m_labeling[nPix] != alpha_label )
                e ->add_term2(variables[i],variables[m_lookupPixVar[nPix]],
                              m_smoothFnPix(pix,nPix,alpha_label,alpha_label) ,
                              m_smoothFnPix(pix,nPix,alpha_label,m_labeling[nPix]),
                              m_smoothFnPix(pix,nPix,m_labeling[pix],alpha_label) ,
                              m_smoothFnPix(pix,nPix,m_labeling[pix],m_labeling[nPix]) );
            else   e ->add_term1(variables[i],m_smoothFnPix(pix,nPix,alpha_label,m_labeling[nPix]),
                                 m_smoothFnPix(pix,nPix,m_labeling[pix],alpha_label));
        }   
        if ( x > 0 )
        {
            nPix = pix - 1;
            if ( m_labeling[nPix] == alpha_label )
               e ->add_term1(variables[i],m_smoothFnPix(pix,nPix,alpha_label,m_labeling[nPix]),
                                 m_smoothFnPix(pix,nPix,m_labeling[pix],alpha_label) );
        }   
        if ( y > 0 )
        {
            nPix = pix - m_width;
            if ( m_labeling[nPix] == alpha_label )
               e ->add_term1(variables[i],m_smoothFnPix(pix,nPix,alpha_label,alpha_label),
                                 m_smoothFnPix(pix,nPix,m_labeling[pix],alpha_label));
        }   
    }
}


GCoptimization::EnergyType Expansion::oneExpansionIteration()
{
    int next;
    terminateOnError( m_dataType == NONE,"You have to set up the data cost before running optimization");
    terminateOnError( m_smoothType == NONE,"You have to set up the smoothness cost before running optimization");
    if (m_random_label_order) scramble_label_table();
    for (next = 0;  next < m_nLabels;  next++ )
        perform_alpha_expansion(m_labelTable[next]);
    return(dataEnergy()+smoothnessEnergy());
}


void Expansion::optimizeAlg(int nIterations)
{
    expansion(nIterations);
}


/***************************************************************************************************************************************/
#define TERMINAL ( (arc_forward *) 1 )      /* to terminal */
#define ORPHAN   ( (arc_forward *) 2 )      /* orphan */
#define INFINITE_D 1000000000       /* infinite distance to the terminal */

/*
    Functions for processing active list.
    i->next points to the next node in the list
    (or to i, if i is the last node in the list).
    If i->next is NULL iff i is not in the list.
    There are two queues. Active nodes are added
    to the end of the second queue and read from
    the front of the first queue. If the first queue
    is empty, it is replaced by the second queue
    (and the second queue becomes empty).
*/
inline void Graph::set_active(node *i)
{
    if (!i->next)
    {
        /* it's not in the list yet */
        if (queue_last[1]) queue_last[1] -> next = i;
        else               queue_first[1]        = i;
        queue_last[1] = i;
        i -> next = i;
    }
}

/*
    Returns the next active node.
    If it is connected to the sink, it stays in the list,
    otherwise it is removed from the list
*/
inline Graph::node * Graph::next_active()
{
    node *i;
    while ( 1 )
    {
        if (!(i=queue_first[0]))
        {
            queue_first[0] = i = queue_first[1];
            queue_last[0]  = queue_last[1];
            queue_first[1] = NULL;
            queue_last[1]  = NULL;
            if (!i) return NULL;
        }
        /* remove it from the active list */
        if (i->next == i) queue_first[0] = queue_last[0] = NULL;
        else              queue_first[0] = i -> next;
        i -> next = NULL;
        /* a node in the list is active iff it has a parent */
        if (i->parent) return i;
    }
}


void Graph::maxflow_init()
{
    node *i;
    node_block *nb;
    queue_first[0] = queue_last[0] = NULL;
    queue_first[1] = queue_last[1] = NULL;
    orphan_first = NULL;
    for (nb=node_block_first; nb; nb=nb->next)
    for (i=&nb->nodes[0]; i<nb->current; i++)
    {
        i -> next = NULL;
        i -> TS = 0;
        if (i->tr_cap > 0)
        {
            /* i is connected to the source */
            i -> is_sink = 0;
            i -> parent = TERMINAL;
            set_active(i);
            i -> TS = 0;
            i -> DIST = 1;
        }
        else if (i->tr_cap < 0)
        {
            /* i is connected to the sink */
            i -> is_sink = 1;
            i -> parent = TERMINAL;
            set_active(i);
            i -> TS = 0;
            i -> DIST = 1;
        }
        else
        {
            i -> parent = NULL;
        }
    }
    TIME = 0;
}


void Graph::augment(node *s_start, node *t_start, captype *cap_middle, captype *rev_cap_middle)
{
    node *i;
    arc_forward *a;
    captype bottleneck;
    nodeptr *np;
    /* 1. Finding bottleneck capacity */
    /* 1a - the source tree */
    bottleneck = *cap_middle;
    for (i=s_start; ; )
    {
        a = i -> parent;
        if (a == TERMINAL) break;
        if (IS_ODD(a))
        {
            a = MAKE_EVEN(a);
            if (bottleneck > a->r_cap) bottleneck = a -> r_cap;
            i = NEIGHBOR_NODE_REV(i, a -> shift);
        }
        else
        {
            if (bottleneck > a->r_rev_cap) bottleneck = a -> r_rev_cap;
            i = NEIGHBOR_NODE(i, a -> shift);
        }
    }
    if (bottleneck > i->tr_cap) bottleneck = i -> tr_cap;
    /* 1b - the sink tree */
    for (i=t_start; ; )
    {
        a = i -> parent;
        if (a == TERMINAL) break;
        if (IS_ODD(a))
        {
            a = MAKE_EVEN(a);
            if (bottleneck > a->r_rev_cap) bottleneck = a -> r_rev_cap;
            i = NEIGHBOR_NODE_REV(i, a -> shift);
        }
        else
        {
            if (bottleneck > a->r_cap) bottleneck = a -> r_cap;
            i = NEIGHBOR_NODE(i, a -> shift);
        }
    }
    if (bottleneck > - i->tr_cap) bottleneck = - i -> tr_cap;
    /* 2. Augmenting */
    /* 2a - the source tree */
    *rev_cap_middle += bottleneck;
    *cap_middle -= bottleneck;
    for (i=s_start; ; )
    {
        a = i -> parent;
        if (a == TERMINAL) break;
        if (IS_ODD(a))
        {
            a = MAKE_EVEN(a);
            a -> r_rev_cap += bottleneck;
            a -> r_cap -= bottleneck;
            if (!a->r_cap)
            {
                /* add i to the adoption list */
                i -> parent = ORPHAN;
                np = nodeptr_block -> New();
                np -> ptr = i;
                np -> next = orphan_first;
                orphan_first = np;
            }
            i = NEIGHBOR_NODE_REV(i, a -> shift);
        }
        else
        {
            a -> r_cap += bottleneck;
            a -> r_rev_cap -= bottleneck;
            if (!a->r_rev_cap)
            {
                /* add i to the adoption list */
                i -> parent = ORPHAN;
                np = nodeptr_block -> New();
                np -> ptr = i;
                np -> next = orphan_first;
                orphan_first = np;
            }
            i = NEIGHBOR_NODE(i, a -> shift);
        }
    }
    i -> tr_cap -= bottleneck;
    if (!i->tr_cap)
    {
        /* add i to the adoption list */
        i -> parent = ORPHAN;
        np = nodeptr_block -> New();
        np -> ptr = i;
        np -> next = orphan_first;
        orphan_first = np;
    }
    /* 2b - the sink tree */
    for (i=t_start; ; )
    {
        a = i -> parent;
        if (a == TERMINAL) break;
        if (IS_ODD(a))
        {
            a = MAKE_EVEN(a);
            a -> r_cap += bottleneck;
            a -> r_rev_cap -= bottleneck;
            if (!a->r_rev_cap)
            {
                /* add i to the adoption list */
                i -> parent = ORPHAN;
                np = nodeptr_block -> New();
                np -> ptr = i;
                np -> next = orphan_first;
                orphan_first = np;
            }
            i = NEIGHBOR_NODE_REV(i, a -> shift);
        }
        else
        {
            a -> r_rev_cap += bottleneck;
            a -> r_cap -= bottleneck;
            if (!a->r_cap)
            {
                /* add i to the adoption list */
                i -> parent = ORPHAN;
                np = nodeptr_block -> New();
                np -> ptr = i;
                np -> next = orphan_first;
                orphan_first = np;
            }
            i = NEIGHBOR_NODE(i, a -> shift);
        }
    }
    i -> tr_cap += bottleneck;
    if (!i->tr_cap)
    {
        /* add i to the adoption list */
        i -> parent = ORPHAN;
        np = nodeptr_block -> New();
        np -> ptr = i;
        np -> next = orphan_first;
        orphan_first = np;
    }
    flow += bottleneck;
}


void Graph::process_source_orphan(node *i)
{
    node *j;
    arc_forward *a0_for, *a0_for_first, *a0_for_last;
    arc_reverse *a0_rev, *a0_rev_first, *a0_rev_last;
    arc_forward *a0_min = NULL, *a;
    nodeptr *np;
    int d, d_min = INFINITE_D;
    /* trying to find a new parent */
    a0_for_first = i -> first_out;
    if (IS_ODD(a0_for_first))
    {
        a0_for_first = (arc_forward *) (((char *)a0_for_first) + 1);
        a0_for_last = (arc_forward *) ((a0_for_first ++) -> shift);
    }
    else a0_for_last = (i + 1) -> first_out;
    a0_rev_first = i -> first_in;
    if (IS_ODD(a0_rev_first))
    {
        a0_rev_first = (arc_reverse *) (((char *)a0_rev_first) + 1);
        a0_rev_last  = (arc_reverse *) ((a0_rev_first ++) -> sister);
    }
    else a0_rev_last = (i + 1) -> first_in;
    for (a0_for=a0_for_first; a0_for<a0_for_last; a0_for++)
    if (a0_for->r_rev_cap)
    {
        j = NEIGHBOR_NODE(i, a0_for -> shift);
        if (!j->is_sink && (a=j->parent))
        {
            /* checking the origin of j */
            d = 0;
            while ( 1 )
            {
                if (j->TS == TIME)
                {
                    d += j -> DIST;
                    break;
                }
                a = j -> parent;
                d ++;
                if (a==TERMINAL)
                {
                    j -> TS = TIME;
                    j -> DIST = 1;
                    break;
                }
                if (a==ORPHAN) { d = INFINITE_D; break; }
                if (IS_ODD(a))
                    j = NEIGHBOR_NODE_REV(j, MAKE_EVEN(a) -> shift);
                else
                    j = NEIGHBOR_NODE(j, a -> shift);
            }
            if (d<INFINITE_D) /* j originates from the source - done */
            {
                if (d<d_min)
                {
                    a0_min = a0_for;
                    d_min = d;
                }
                /* set marks along the path */
                for (j=NEIGHBOR_NODE(i, a0_for->shift); j->TS!=TIME; )
                {
                    j -> TS = TIME;
                    j -> DIST = d --;
                    a = j->parent;
                    if (IS_ODD(a))
                        j = NEIGHBOR_NODE_REV(j, MAKE_EVEN(a) -> shift);
                    else
                        j = NEIGHBOR_NODE(j, a -> shift);
                }
            }
        }
    }
    for (a0_rev=a0_rev_first; a0_rev<a0_rev_last; a0_rev++)
    {
        a0_for = a0_rev -> sister;
        if (a0_for->r_cap)
        {
            j = NEIGHBOR_NODE_REV(i, a0_for -> shift);
            if (!j->is_sink && (a=j->parent))
            {
                /* checking the origin of j */
                d = 0;
                while ( 1 )
                {
                    if (j->TS == TIME)
                    {
                        d += j -> DIST;
                        break;
                    }
                    a = j -> parent;
                    d ++;
                    if (a==TERMINAL)
                    {
                        j -> TS = TIME;
                        j -> DIST = 1;
                        break;
                    }
                    if (a==ORPHAN) { d = INFINITE_D; break; }
                    if (IS_ODD(a))
                        j = NEIGHBOR_NODE_REV(j, MAKE_EVEN(a) -> shift);
                    else
                        j = NEIGHBOR_NODE(j, a -> shift);
                }
                if (d<INFINITE_D) /* j originates from the source - done */
                {
                    if (d<d_min)
                    {
                        a0_min = MAKE_ODD(a0_for);
                        d_min = d;
                    }
                    /* set marks along the path */
                    for (j=NEIGHBOR_NODE_REV(i,a0_for->shift); j->TS!=TIME; )
                    {
                        j -> TS = TIME;
                        j -> DIST = d --;
                        a = j->parent;
                        if (IS_ODD(a))
                            j = NEIGHBOR_NODE_REV(j, MAKE_EVEN(a) -> shift);
                        else
                            j = NEIGHBOR_NODE(j, a -> shift);
                    }
                }
            }
        }
    }
    if ((i->parent = a0_min))
    {
        i -> TS = TIME;
        i -> DIST = d_min + 1;
    }
    else
    {
        /* no parent is found */
        i -> TS = 0;
        /* process neighbors */
        for (a0_for=a0_for_first; a0_for<a0_for_last; a0_for++)
        {
            j = NEIGHBOR_NODE(i, a0_for -> shift);
            if (!j->is_sink && (a=j->parent))
            {
                if (a0_for->r_rev_cap) set_active(j);
                if (a!=TERMINAL && a!=ORPHAN && IS_ODD(a) && NEIGHBOR_NODE_REV(j, MAKE_EVEN(a)->shift)==i)
                {
                    /* add j to the adoption list */
                    j -> parent = ORPHAN;
                    np = nodeptr_block -> New();
                    np -> ptr = j;
                    if (orphan_last) orphan_last -> next = np;
                    else             orphan_first        = np;
                    orphan_last = np;
                    np -> next = NULL;
                }
            }
        }
        for (a0_rev=a0_rev_first; a0_rev<a0_rev_last; a0_rev++)
        {
            a0_for = a0_rev -> sister;
            j = NEIGHBOR_NODE_REV(i, a0_for -> shift);
            if (!j->is_sink && (a=j->parent))
            {
                if (a0_for->r_cap) set_active(j);
                if (a!=TERMINAL && a!=ORPHAN && !IS_ODD(a) && NEIGHBOR_NODE(j, a->shift)==i)
                {
                    /* add j to the adoption list */
                    j -> parent = ORPHAN;
                    np = nodeptr_block -> New();
                    np -> ptr = j;
                    if (orphan_last) orphan_last -> next = np;
                    else             orphan_first        = np;
                    orphan_last = np;
                    np -> next = NULL;
                }
            }
        }
    }
}


void Graph::process_sink_orphan(node *i)
{
    node *j;
    arc_forward *a0_for, *a0_for_first, *a0_for_last;
    arc_reverse *a0_rev, *a0_rev_first, *a0_rev_last;
    arc_forward *a0_min = NULL, *a;
    nodeptr *np;
    int d, d_min = INFINITE_D;
    /* trying to find a new parent */
    a0_for_first = i -> first_out;
    if (IS_ODD(a0_for_first))
    {
        a0_for_first = (arc_forward *) (((char *)a0_for_first) + 1);
        a0_for_last = (arc_forward *) ((a0_for_first ++) -> shift);
    }
    else a0_for_last = (i + 1) -> first_out;
    a0_rev_first = i -> first_in;
    if (IS_ODD(a0_rev_first))
    {
        a0_rev_first = (arc_reverse *) (((char *)a0_rev_first) + 1);
        a0_rev_last  = (arc_reverse *) ((a0_rev_first ++) -> sister);
    }
    else a0_rev_last = (i + 1) -> first_in;
    for (a0_for=a0_for_first; a0_for<a0_for_last; a0_for++)
    if (a0_for->r_cap)
    {
        j = NEIGHBOR_NODE(i, a0_for -> shift);
        if (j->is_sink && (a=j->parent))
        {
            /* checking the origin of j */
            d = 0;
            while ( 1 )
            {
                if (j->TS == TIME)
                {
                    d += j -> DIST;
                    break;
                }
                a = j -> parent;
                d ++;
                if (a==TERMINAL)
                {
                    j -> TS = TIME;
                    j -> DIST = 1;
                    break;
                }
                if (a==ORPHAN) { d = INFINITE_D; break; }
                if (IS_ODD(a))
                    j = NEIGHBOR_NODE_REV(j, MAKE_EVEN(a) -> shift);
                else
                    j = NEIGHBOR_NODE(j, a -> shift);
            }
            if (d<INFINITE_D) /* j originates from the sink - done */
            {
                if (d<d_min)
                {
                    a0_min = a0_for;
                    d_min = d;
                }
                /* set marks along the path */
                for (j=NEIGHBOR_NODE(i, a0_for->shift); j->TS!=TIME; )
                {
                    j -> TS = TIME;
                    j -> DIST = d --;
                    a = j->parent;
                    if (IS_ODD(a))
                        j = NEIGHBOR_NODE_REV(j, MAKE_EVEN(a) -> shift);
                    else
                        j = NEIGHBOR_NODE(j, a -> shift);
                }
            }
        }
    }
    for (a0_rev=a0_rev_first; a0_rev<a0_rev_last; a0_rev++)
    {
        a0_for = a0_rev -> sister;
        if (a0_for->r_rev_cap)
        {
            j = NEIGHBOR_NODE_REV(i, a0_for -> shift);
            if (j->is_sink && (a=j->parent))
            {
                /* checking the origin of j */
                d = 0;
                while ( 1 )
                {
                    if (j->TS == TIME)
                    {
                        d += j -> DIST;
                        break;
                    }
                    a = j -> parent;
                    d ++;
                    if (a==TERMINAL)
                    {
                        j -> TS = TIME;
                        j -> DIST = 1;
                        break;
                    }
                    if (a==ORPHAN) { d = INFINITE_D; break; }
                    if (IS_ODD(a))
                        j = NEIGHBOR_NODE_REV(j, MAKE_EVEN(a) -> shift);
                    else
                        j = NEIGHBOR_NODE(j, a -> shift);
                }
                if (d<INFINITE_D) /* j originates from the sink - done */
                {
                    if (d<d_min)
                    {
                        a0_min = MAKE_ODD(a0_for);
                        d_min = d;
                    }
                    /* set marks along the path */
                    for (j=NEIGHBOR_NODE_REV(i,a0_for->shift); j->TS!=TIME; )
                    {
                        j -> TS = TIME;
                        j -> DIST = d --;
                        a = j->parent;
                        if (IS_ODD(a))
                            j = NEIGHBOR_NODE_REV(j, MAKE_EVEN(a) -> shift);
                        else
                            j = NEIGHBOR_NODE(j, a -> shift);
                    }
                }
            }
        }
    }
    if ((i->parent = a0_min))
    {
        i -> TS = TIME;
        i -> DIST = d_min + 1;
    }
    else
    {
        /* no parent is found */
        i -> TS = 0;
        /* process neighbors */
        for (a0_for=a0_for_first; a0_for<a0_for_last; a0_for++)
        {
            j = NEIGHBOR_NODE(i, a0_for -> shift);
            if (j->is_sink && (a=j->parent))
            {
                if (a0_for->r_cap) set_active(j);
                if (a!=TERMINAL && a!=ORPHAN && IS_ODD(a) && NEIGHBOR_NODE_REV(j, MAKE_EVEN(a)->shift)==i)
                {
                    /* add j to the adoption list */
                    j -> parent = ORPHAN;
                    np = nodeptr_block -> New();
                    np -> ptr = j;
                    if (orphan_last) orphan_last -> next = np;
                    else             orphan_first        = np;
                    orphan_last = np;
                    np -> next = NULL;
                }
            }
        }
        for (a0_rev=a0_rev_first; a0_rev<a0_rev_last; a0_rev++)
        {
            a0_for = a0_rev -> sister;
            j = NEIGHBOR_NODE_REV(i, a0_for -> shift);
            if (j->is_sink && (a=j->parent))
            {
                if (a0_for->r_rev_cap) set_active(j);
                if (a!=TERMINAL && a!=ORPHAN && !IS_ODD(a) && NEIGHBOR_NODE(j, a->shift)==i)
                {
                    /* add j to the adoption list */
                    j -> parent = ORPHAN;
                    np = nodeptr_block -> New();
                    np -> ptr = j;
                    if (orphan_last) orphan_last -> next = np;
                    else             orphan_first        = np;
                    orphan_last = np;
                    np -> next = NULL;
                }
            }
        }
    }
}


Graph::flowtype Graph::maxflow()
{
    node *i, *j, *current_node = NULL, *s_start, *t_start = NULL;
    captype *cap_middle = NULL, *rev_cap_middle = NULL;
    arc_forward *a_for, *a_for_first, *a_for_last;
    arc_reverse *a_rev, *a_rev_first, *a_rev_last;
    nodeptr *np, *np_next;
    prepare_graph();
    maxflow_init();
    nodeptr_block = new DBlock<nodeptr>(NODEPTR_BLOCK_SIZE, error_function);
    while ( 1 )
    {
        if ((i=current_node))
        {
            i -> next = NULL; /* remove active flag */
            if (!i->parent) i = NULL;
        }
        if (!i)
        {
            if (!(i = next_active())) break;
        }
        /* growth */
        s_start = NULL;
        a_for_first = i -> first_out;
        if (IS_ODD(a_for_first))
        {
            a_for_first = (arc_forward *) (((char *)a_for_first) + 1);
            a_for_last = (arc_forward *) ((a_for_first ++) -> shift);
        }
        else a_for_last = (i + 1) -> first_out;
        a_rev_first = i -> first_in;
        if (IS_ODD(a_rev_first))
        {
            a_rev_first = (arc_reverse *) (((char *)a_rev_first) + 1);
            a_rev_last = (arc_reverse *) ((a_rev_first ++) -> sister);
        }
        else a_rev_last = (i + 1) -> first_in;
        if (!i->is_sink)
        {
            /* grow source tree */
            for (a_for=a_for_first; a_for<a_for_last; a_for++)
            if (a_for->r_cap)
            {
                j = NEIGHBOR_NODE(i, a_for -> shift);
                if (!j->parent)
                {
                    j -> is_sink = 0;
                    j -> parent = MAKE_ODD(a_for);
                    j -> TS = i -> TS;
                    j -> DIST = i -> DIST + 1;
                    set_active(j);
                }
                else if (j->is_sink)
                {
                    s_start = i;
                    t_start = j;
                    cap_middle     = & ( a_for -> r_cap );
                    rev_cap_middle = & ( a_for -> r_rev_cap );
                    break;
                }
                else if (j->TS <= i->TS &&
                         j->DIST > i->DIST)
                {
                    /* heuristic - trying to make the distance from j to the source shorter */
                    j -> parent = MAKE_ODD(a_for);
                    j -> TS = i -> TS;
                    j -> DIST = i -> DIST + 1;
                }
            }
            if (!s_start)
            for (a_rev=a_rev_first; a_rev<a_rev_last; a_rev++)
            {
                a_for = a_rev -> sister;
                if (a_for->r_rev_cap)
                {
                    j = NEIGHBOR_NODE_REV(i, a_for -> shift);
                    if (!j->parent)
                    {
                        j -> is_sink = 0;
                        j -> parent = a_for;
                        j -> TS = i -> TS;
                        j -> DIST = i -> DIST + 1;
                        set_active(j);
                    }
                    else if (j->is_sink)
                    {
                        s_start = i;
                        t_start = j;
                        cap_middle     = & ( a_for -> r_rev_cap );
                        rev_cap_middle = & ( a_for -> r_cap );
                        break;
                    }
                    else if (j->TS <= i->TS &&
                             j->DIST > i->DIST)
                    {
                        /* heuristic - trying to make the distance from j to the source shorter */
                        j -> parent = a_for;
                        j -> TS = i -> TS;
                        j -> DIST = i -> DIST + 1;
                    }
                }
            }
        }
        else
        {
            /* grow sink tree */
            for (a_for=a_for_first; a_for<a_for_last; a_for++)
            if (a_for->r_rev_cap)
            {
                j = NEIGHBOR_NODE(i, a_for -> shift);
                if (!j->parent)
                {
                    j -> is_sink = 1;
                    j -> parent = MAKE_ODD(a_for);
                    j -> TS = i -> TS;
                    j -> DIST = i -> DIST + 1;
                    set_active(j);
                }
                else if (!j->is_sink)
                {
                    s_start = j;
                    t_start = i;
                    cap_middle     = & ( a_for -> r_rev_cap );
                    rev_cap_middle = & ( a_for -> r_cap );
                    break;
                }
                else if (j->TS <= i->TS &&
                         j->DIST > i->DIST)
                {
                    /* heuristic - trying to make the distance from j to the sink shorter */
                    j -> parent = MAKE_ODD(a_for);
                    j -> TS = i -> TS;
                    j -> DIST = i -> DIST + 1;
                }
            }
            for (a_rev=a_rev_first; a_rev<a_rev_last; a_rev++)
            {
                a_for = a_rev -> sister;
                if (a_for->r_cap)
                {
                    j = NEIGHBOR_NODE_REV(i, a_for -> shift);
                    if (!j->parent)
                    {
                        j -> is_sink = 1;
                        j -> parent = a_for;
                        j -> TS = i -> TS;
                        j -> DIST = i -> DIST + 1;
                        set_active(j);
                    }
                    else if (!j->is_sink)
                    {
                        s_start = j;
                        t_start = i;
                        cap_middle     = & ( a_for -> r_cap );
                        rev_cap_middle = & ( a_for -> r_rev_cap );
                        break;
                    }
                    else if (j->TS <= i->TS &&
                             j->DIST > i->DIST)
                    {
                        /* heuristic - trying to make the distance from j to the sink shorter */
                        j -> parent = a_for;
                        j -> TS = i -> TS;
                        j -> DIST = i -> DIST + 1;
                    }
                }
            }
        }
        TIME ++;
        if (s_start)
        {
            i -> next = i; /* set active flag */
            current_node = i;
            /* augmentation */
            augment(s_start, t_start, cap_middle, rev_cap_middle);
            /* augmentation end */
            /* adoption */
            while ((np=orphan_first))
            {
                np_next = np -> next;
                np -> next = NULL;
                while ((np=orphan_first))
                {
                    orphan_first = np -> next;
                    i = np -> ptr;
                    nodeptr_block -> Delete(np);
                    if (!orphan_first) orphan_last = NULL;
                    if (i->is_sink) process_sink_orphan(i);
                    else            process_source_orphan(i);
                }
                orphan_first = np_next;
            }
            /* adoption end */
        }
        else current_node = NULL;
    }
    delete nodeptr_block;
    return flow;
}


Graph::termtype Graph::what_segment(node_id i)
{
    if (((node*)i)->parent && !((node*)i)->is_sink) return SOURCE;
    return SINK;
}


/***************************************************************************************************************************************/
int numIterRun;
// Some of the GBP code has been disabled here

#ifdef __GNUC__
#define _finite(n) finite(n)
#endif

#define mexPrintf printf
#define mexErrMsgTxt printf
OneNodeCluster::OneNodeCluster()
{
}


int OneNodeCluster::numStates;
float int_exp(float val, int exp)
{
  float tmp=1.0;
  for(int i =0; i < exp; i++)
    tmp *=val;
  return tmp;
}

void getPsiMat(OneNodeCluster & /*cluster*/, float *&destMatrix, int r, int c, MaxProdBP *mrf, int direction)
{
  int mrfHeight = mrf->getHeight();
  int mrfWidth = mrf->getWidth();
  int numLabels = mrf->getNLabels();

  float *currMatrix = mrf->getScratchMatrix();

  if(((direction==/*OneNodeCluster::*/UP) &&(r==0)) ||
     ((direction==/*OneNodeCluster::*/DOWN) &&(r==(mrfHeight-1)))||
     ((direction==/*OneNodeCluster::*/LEFT) &&(c==0))||
     ((direction==/*OneNodeCluster::*/RIGHT) &&(c==(mrfWidth-1))))
  {
    for(int i=0; i < numLabels * numLabels; i++)
    {
      currMatrix[i] = 1.0;
    }
  }
  else
  {
    int weight_mod = 1;
    if(mrf->varWeights())
/*    {

      switch(direction)
      {
      case OneNodeCluster::LEFT:
        weight_mod = mrf->getHorizWeight(r,c-1);
        break;
      case OneNodeCluster::RIGHT:
        weight_mod = mrf->getHorizWeight(r,c);
        break;
      case OneNodeCluster::UP:
        weight_mod = mrf->getVertWeight(r-1,c);
        break;
      case OneNodeCluster::DOWN:
        weight_mod = mrf->getVertWeight(r,c);
        break;
      }
    }
        */
        {
                if ( direction == /*OneNodeCluster::*/LEFT ) weight_mod = int(mrf->getHorizWeight(r,c-1));
                else if ( direction == /*OneNodeCluster::*/RIGHT ) weight_mod = int(mrf->getHorizWeight(r,c));
                else if ( direction == /*OneNodeCluster::*/UP ) weight_mod = int(mrf->getVertWeight(r-1,c));
                else    weight_mod = int(mrf->getVertWeight(r,c));
        
        }
    if(weight_mod!=1)
    {
      for(int i = 0; i < numLabels*numLabels; i++)
      {
        currMatrix[i] = int_exp(mrf->getExpV(i),weight_mod);
      }
      destMatrix = currMatrix;
    }
    else
    {
      destMatrix = mrf->getExpV();
    }
  }
}
  
void initOneNodeMsgMem(OneNodeCluster *nodeArray, float *memChunk, 
                       const int numNodes, const int msgChunkSize)
{
  float *currPtr = memChunk;
  OneNodeCluster *currNode = nodeArray;
  for(int i = 0; i < numNodes; i++)
  {

    currNode->receivedMsgs[0] = currPtr; currPtr+=msgChunkSize;
    currNode->receivedMsgs[1] = currPtr; currPtr+=msgChunkSize;
    currNode->receivedMsgs[2] = currPtr; currPtr+=msgChunkSize;
    currNode->receivedMsgs[3] = currPtr; currPtr+=msgChunkSize;
    currNode->nextRoundReceivedMsgs[0] = currPtr; currPtr +=msgChunkSize;
    currNode->nextRoundReceivedMsgs[1] = currPtr; currPtr +=msgChunkSize;
    currNode->nextRoundReceivedMsgs[2] = currPtr; currPtr +=msgChunkSize;
    currNode->nextRoundReceivedMsgs[3] = currPtr; currPtr +=msgChunkSize;
  
    currNode++;
  }
}


void OneNodeCluster::ComputeMsgRight(float *msgDest, int r, int c, MaxProdBP *mrf)
{
  float *nodeLeftMsg = receivedMsgs[LEFT],
    *nodeDownMsg = receivedMsgs[DOWN],
    *nodeUpMsg =    receivedMsgs[UP];

  float *psiMat;

  getPsiMat(*this,psiMat,r,c,mrf,/*OneNodeCluster::*/RIGHT);
  
  float *cmessage = msgDest,
    total = 0;
  
  for(int rightNodeInd = 0; rightNodeInd < numStates; rightNodeInd++)
  {

    *cmessage = 0;
    for(int leftNodeInd = 0; leftNodeInd < numStates; leftNodeInd++)
    {
      float tmp =   nodeLeftMsg[leftNodeInd] * 
        nodeUpMsg[leftNodeInd] * 
        nodeDownMsg[leftNodeInd] * 
        localEv[leftNodeInd] *
        psiMat[leftNodeInd * numStates + rightNodeInd];

      if(tmp > *cmessage)
        *cmessage = tmp;

      if ((*cmessage != *cmessage)||!(_finite(*cmessage))||(*cmessage < 1e-10))
    {
      mexPrintf("%f %f %f %f %e\n ",nodeLeftMsg[leftNodeInd],nodeUpMsg[leftNodeInd],nodeDownMsg[leftNodeInd] ,
                psiMat[leftNodeInd * numStates + rightNodeInd], localEv[leftNodeInd]);
      
      mexErrMsgTxt("Break Here\n");

      assert(0);
    } 
    }
//      if (*cmessage < 0.000001)
//        *cmessage=0.000001;
    total += *cmessage;
    cmessage++;
  }

  int errFlag = 0;
  for(int i = 0; i < numStates; i++)
  {
    
    msgDest[i] /= total;
    if (msgDest[i] != msgDest[i])
      errFlag=1;
//      if (msgDest[i] < 0.000001)
//        msgDest[i] = 0.000001;
  }
  if(errFlag)
  {
    mexPrintf("%f |",total);
    for(int i = 0; i < numStates; i++)
    {
      mexPrintf("%f ",msgDest[i]);
    }
    mexErrMsgTxt(" ");
    assert(0);
  }
}

// This means, "Compute the message to send left."

void OneNodeCluster::ComputeMsgLeft(float *msgDest, int r, int c, MaxProdBP *mrf)
{
  float *nodeRightMsg = receivedMsgs[RIGHT],
    *nodeDownMsg = receivedMsgs[DOWN],
    *nodeUpMsg =    receivedMsgs[UP];

  float *psiMat;

  getPsiMat(*this,psiMat,r,c,mrf,/*OneNodeCluster::*/LEFT);
  
  float *cmessage = msgDest,
    total = 0;
  
  for(int leftNodeInd = 0; leftNodeInd < numStates; leftNodeInd++)
  {

    *cmessage = 0;
    for(int rightNodeInd = 0; rightNodeInd < numStates; rightNodeInd++)
    {
       float tmp =  nodeRightMsg[rightNodeInd] * 
        nodeUpMsg[rightNodeInd] * 
        nodeDownMsg[rightNodeInd] *
        localEv[rightNodeInd] *
        psiMat[leftNodeInd * numStates + rightNodeInd] ; 
       
       if(tmp > *cmessage)
         *cmessage = tmp;

      if ((*cmessage != *cmessage)||!(_finite(*cmessage)))
      {
        mexPrintf("%f %f %f %f \n",nodeRightMsg[rightNodeInd] ,nodeUpMsg[rightNodeInd],nodeDownMsg[rightNodeInd],
                  psiMat[leftNodeInd * numStates + rightNodeInd]);
        
        mexErrMsgTxt("Break Here\n");
      assert(0);
      }
    }
//        if (*cmessage < 0.000001)
//      *cmessage=0.000001;
    total += *cmessage;
    cmessage++;
  }

  for(int i = 0; i < numStates; i++)
  {
    msgDest[i] /= total;
//      if (msgDest[i] < 0.000001)
//        msgDest[i] = 0.000001;
  }
}

void OneNodeCluster::ComputeMsgUp(float *msgDest, int r, int c, MaxProdBP *mrf)
{
  float *nodeRightMsg = receivedMsgs[RIGHT],
    *nodeDownMsg = receivedMsgs[DOWN],
    *nodeLeftMsg =    receivedMsgs[LEFT];
 
  float *psiMat;

  getPsiMat(*this,psiMat,r,c,mrf,/*OneNodeCluster::*/UP);
  
  float *cmessage = msgDest,
    total = 0;
  
  for(int upNodeInd = 0; upNodeInd < numStates; upNodeInd++)
  {

    *cmessage = 0;
    for(int downNodeInd = 0; downNodeInd < numStates; downNodeInd++)
    {
      float tmp = nodeRightMsg[downNodeInd] * 
        nodeLeftMsg[downNodeInd] * 
        nodeDownMsg[downNodeInd] * 
        localEv[downNodeInd] *
        psiMat[upNodeInd * numStates + downNodeInd] ; 

       if(tmp > *cmessage)
         *cmessage = tmp;


      if ((*cmessage != *cmessage)||!_finite(*cmessage))
    {
      mexPrintf("%f %f %f %f\n ",nodeRightMsg[downNodeInd],nodeLeftMsg[downNodeInd],nodeDownMsg[downNodeInd],psiMat[upNodeInd * numStates + downNodeInd]);
      
      mexErrMsgTxt("Break Here\n");
      assert(0);
    } 
    }
//      if (*cmessage < 0.000001)
//        *cmessage=0.000001;
    total += *cmessage;
    cmessage++;
  }

  for(int i = 0; i < numStates; i++)
  {
    msgDest[i] /= total;
//      if (msgDest[i] < 0.000001)
//        msgDest[i] = 0.000001;
  }
}

void OneNodeCluster::ComputeMsgDown(float *msgDest, int r, int c, MaxProdBP *mrf)
{
  float *nodeRightMsg = receivedMsgs[RIGHT],
    *nodeUpMsg = receivedMsgs[UP],
    *nodeLeftMsg =    receivedMsgs[LEFT];
 
  float *psiMat;

  getPsiMat(*this,psiMat,r,c,mrf,/*OneNodeCluster::*/DOWN);
  
  float *cmessage = msgDest,
    total = 0;
  
  for(int downNodeInd = 0; downNodeInd < numStates; downNodeInd++)
  {

    *cmessage = 0;
    for(int upNodeInd = 0; upNodeInd < numStates; upNodeInd++)
    {
      float tmp =   nodeRightMsg[upNodeInd] * 
        nodeLeftMsg[upNodeInd] * 
        nodeUpMsg[upNodeInd] * 
        localEv[upNodeInd] *
        psiMat[upNodeInd * numStates + downNodeInd] ; 

       if(tmp > *cmessage)
         *cmessage = tmp;

    }
    if (*cmessage != *cmessage)
    {      printf("Break Here\n");
//      if (*cmessage < 0.000001)
//        *cmessage=0.000001;
      assert(0);
    }
   total += *cmessage;
    cmessage++;
  }

  for(int i = 0; i < numStates; i++)
  {
    msgDest[i] /= total;
//      if (msgDest[i] < 0.000001)
//        msgDest[i] = 0.000001;
  }

}




//  void TwoNodeCluster::deliverMsgs()
//  {
//    float *tmp;
  
//    tmp = nextRoundReceivedMsgs[0];
//    nextRoundReceivedMsgs[0] = receivedMsgs[0];
//    receivedMsgs[0] = tmp;

//    tmp = nextRoundReceivedMsgs[1];
//    nextRoundReceivedMsgs[1] = receivedMsgs[1];
//    receivedMsgs[1] = tmp;


//  }

//  void OneNodeCluster::deliverMsgs()
//  {
//      float *tmp;

//      tmp =nextRoundReceivedMsgs[UP];
//      nextRoundReceivedMsgs[UP] =receivedMsgs[UP];
//      receivedMsgs[UP] = tmp;

//      tmp = nextRoundReceivedMsgs[DOWN];
//      nextRoundReceivedMsgs[DOWN] = receivedMsgs[DOWN];
//      receivedMsgs[DOWN] = tmp;

//      tmp = nextRoundReceivedMsgs[LEFT];
//      nextRoundReceivedMsgs[LEFT] = receivedMsgs[LEFT];
//      receivedMsgs[LEFT] = tmp;

//      tmp = nextRoundReceivedMsgs[RIGHT];
//      nextRoundReceivedMsgs[RIGHT] = receivedMsgs[RIGHT];
//      receivedMsgs[RIGHT] = tmp;


//  }


void OneNodeCluster::deliverMsgs()
{
  const double alpha = 0.8,omalpha = 1-alpha;
  int i;
        for( i = 0; i < numStates; i++)
        {
          receivedMsgs[UP][i] = float(alpha * receivedMsgs[UP][i] + omalpha * nextRoundReceivedMsgs[UP][i]);
        }
        for( i = 0; i < numStates; i++)
        {
          receivedMsgs[DOWN][i] = float(alpha * receivedMsgs[DOWN][i] + omalpha * nextRoundReceivedMsgs[DOWN][i]);
        }
        for( i = 0; i < numStates; i++)
        {
          receivedMsgs[LEFT][i] = float(alpha * receivedMsgs[LEFT][i] + omalpha * nextRoundReceivedMsgs[LEFT][i]);
        }
        for( i = 0; i < numStates; i++)
        {
          receivedMsgs[RIGHT][i] = float(alpha * receivedMsgs[RIGHT][i] + omalpha * nextRoundReceivedMsgs[RIGHT][i]);
        }

}

void OneNodeCluster::getBelief(float *beliefVec)
{
        float sum=0;
	int i;
        for(i = 0; i < numStates; i++)
        {
                beliefVec[i] = receivedMsgs[UP][i] * receivedMsgs[DOWN][i] *
                  receivedMsgs[LEFT][i] * receivedMsgs[RIGHT][i] * localEv[i];
                if(!_finite(beliefVec[i]))
                {
                  mexPrintf("SPROBLEM!\n %f %f %f  %f\n",receivedMsgs[UP][i] , receivedMsgs[DOWN][i] ,
                         receivedMsgs[LEFT][i] , receivedMsgs[RIGHT][i] );
                  mexErrMsgTxt(" ");
                  assert(0);

                }
                sum += beliefVec[i];
        }

        for( i = 0; i < numStates; i++)
        {
                beliefVec[i] /= sum;
        }
}


void computeMessagesLeftRight(OneNodeCluster *nodeArray, const int numCols, const int /*numRows*/, const int currRow, const float alpha, MaxProdBP *mrf)
{
  const int numStates = OneNodeCluster::numStates;
  const float omalpha = float(1.0 - alpha);
  int col, i;
  for(col = 0; col < numCols-1; col++)
  {
    nodeArray[currRow * numCols + col].ComputeMsgRight(nodeArray[currRow * numCols + col+1].nextRoundReceivedMsgs[/*OneNodeCluster::*/LEFT],currRow, col, mrf);
    for(i = 0; i < numStates; i++)
    {
      nodeArray[currRow * numCols + col+1].receivedMsgs[/*OneNodeCluster::*/LEFT][i] = 
        omalpha * nodeArray[currRow * numCols + col+1].receivedMsgs[/*OneNodeCluster::*/LEFT][i] + 
        alpha * nodeArray[currRow * numCols + col+1].nextRoundReceivedMsgs[/*OneNodeCluster::*/LEFT][i];
    }
  } 
  for( col = numCols-1; col > 0; col--)
  {
    nodeArray[currRow * numCols + col].ComputeMsgLeft(nodeArray[currRow * numCols + col-1].nextRoundReceivedMsgs[/*OneNodeCluster::*/RIGHT],currRow, col, mrf);
    for(i = 0; i < numStates; i++)
    {
      nodeArray[currRow * numCols + col-1].receivedMsgs[/*OneNodeCluster::*/RIGHT][i] = 
        omalpha * nodeArray[currRow * numCols + col-1].receivedMsgs[/*OneNodeCluster::*/RIGHT][i] + 
        alpha * nodeArray[currRow * numCols + col-1].nextRoundReceivedMsgs[/*OneNodeCluster::*/RIGHT][i];
    }
  } 

}

void computeMessagesUpDown(OneNodeCluster *nodeArray, const int numCols, const int numRows, const int currCol, const float alpha, MaxProdBP *mrf)
{
  const int numStates = OneNodeCluster::numStates;
  const float omalpha = float(1.0 - alpha);
  int row;
  for(row = 0; row < numRows-1; row++)
  {
    nodeArray[row * numCols + currCol].ComputeMsgDown(nodeArray[(row+1) * numCols + currCol].nextRoundReceivedMsgs[/*OneNodeCluster::*/UP],row, currCol, mrf);
    for(int i = 0; i < numStates; i++)
    {
      nodeArray[(row+1) * numCols + currCol].receivedMsgs[/*OneNodeCluster::*/UP][i] = 
        omalpha * nodeArray[(row+1) * numCols + currCol].receivedMsgs[/*OneNodeCluster::*/UP][i] + 
        alpha * nodeArray[(row+1) * numCols + currCol].nextRoundReceivedMsgs[/*OneNodeCluster::*/UP][i];
    }
  } 
  for(row = numRows-1; row > 0; row--)
  {
    nodeArray[row * numCols + currCol].ComputeMsgUp(nodeArray[(row-1) * numCols + currCol].nextRoundReceivedMsgs[/*OneNodeCluster::*/DOWN], row, currCol, mrf);
    for(int i = 0; i < numStates; i++)
    {
      nodeArray[(row-1) * numCols + currCol].receivedMsgs[/*OneNodeCluster::*/DOWN][i] = 
        omalpha * nodeArray[(row-1) * numCols + currCol].receivedMsgs[/*OneNodeCluster::*/DOWN][i] + 
        alpha * nodeArray[(row-1) * numCols + currCol].nextRoundReceivedMsgs[/*OneNodeCluster::*/DOWN][i];
    }
  } 

}


/***************************************************************************************************************************************/
#define m_D(pix,l)  m_D[(pix)*m_nLabels+(l)]
#define m_V(l1,l2)  m_V[(l1)*m_nLabels+(l2)]


MaxProdBP::MaxProdBP(int width, int height, int nLabels,EnergyFunction *eng):MRFmodel(width,height,nLabels,eng)
{
    m_needToFreeV = 0;
    initializeAlg();
}


MaxProdBP::MaxProdBP(int nPixels, int nLabels,EnergyFunction *eng):MRFmodel(nPixels,nLabels,eng)
{
    m_needToFreeV = 0;
    initializeAlg();
}


MaxProdBP::~MaxProdBP()
{ 
    delete[] m_answer;
    if (m_message_chunk) delete[] m_message_chunk;
    if (!m_grid_graph) delete[] m_neighbors;
    if ( m_needToFreeV ) delete[] m_V;
}


void MaxProdBP::initializeAlg()
{
  m_exp_scale = 0;
    m_answer = (Label *) new Label[m_nPixels];
    if ( !m_answer ){printf("\nNot enough memory, exiting");exit(0);}
    m_scratchMatrix = new float[m_nLabels * m_nLabels];
    nodeArray =     new OneNodeCluster[m_nPixels];
    OneNodeCluster::numStates = m_nLabels;
    if (!m_grid_graph)
    {
      assert(0);
      // Only Grid Graphs are supported
        m_neighbors = (LinkedBlockList *) new LinkedBlockList[m_nPixels];
        if (!m_neighbors) {printf("Not enough memory,exiting");exit(0);};
    }
    else
    {
      m_message_chunk = (float *) new float[2*4*m_nPixels * m_nLabels];
      if ( !m_message_chunk ){printf("\nNot enough memory for messages, exiting");exit(0);}
      for(int i = 0; i < 2*4*m_nPixels*m_nLabels; i++)
        m_message_chunk[i]=1.0;
      initOneNodeMsgMem(nodeArray, m_message_chunk, m_nPixels, m_nLabels);
    }
}


int MaxProdBP::getNLabels()
{
  return m_nLabels;
}


int MaxProdBP::getWidth()
{
  return m_width;
}


int MaxProdBP::getHeight()
{
  return m_height;
}


float MaxProdBP::getExpV(int i)
{
  return m_ExpData[i];
}


float *MaxProdBP::getExpV()
{
  return m_ExpData;
}


MRFmodel::CostVal MaxProdBP::getHorizWeight(int r, int c)
{
  int x = c;
  int y = r;
  int pix    = x+y*m_width;
  return m_varWeights ? m_horizWeights[pix] :  1;
}


MRFmodel::CostVal MaxProdBP::getVertWeight(int r, int c)
{
  int x = c;
  int y = r;
  int pix    = x+y*m_width;
  return  m_varWeights ? m_vertWeights[pix] :  1;
}


bool MaxProdBP::varWeights()
{
  return m_varWeights;
}


float *MaxProdBP::getScratchMatrix()
{
  return m_scratchMatrix;
}


void MaxProdBP::clearAnswer()
{
    memset(m_answer, 0, m_nPixels*sizeof(Label));
}


void MaxProdBP::setNeighbors(int pixel1, int pixel2, CostVal weight)
{
	assert(0);
	//Only Grid Graphs are supported
    assert(!m_grid_graph);
    assert(pixel1 < m_nPixels && pixel1 >= 0 && pixel2 < m_nPixels && pixel2 >= 0);
    Neighbor *temp1 = (Neighbor *) new Neighbor;
    Neighbor *temp2 = (Neighbor *) new Neighbor;
    if ( !temp1 || ! temp2 ) {printf("\nNot enough memory, exiting");exit(0);}
    temp1->weight  = weight;
    temp1->to_node = pixel2;
    temp2->weight  = weight;
    temp2->to_node = pixel1;
    m_neighbors[pixel1].addFront(temp1);
    m_neighbors[pixel2].addFront(temp2);
}


/*void MaxProdBP::clearNeighbors()
{
	m_neighbors->clearList();
}*/


MRFmodel::EnergyVal MaxProdBP::smoothnessEnergy()
{
    EnergyVal eng = (EnergyVal) 0;
    EnergyVal weight;
    int x,y,pix;
    if ( m_grid_graph )
    {
        if ( m_smoothType != FUNCTION  )
        {
            for ( y = 0; y < m_height; y++ )
              for ( x = 1; x < m_width; x++ )
              {
                pix    = x+y*m_width;
                weight = m_varWeights ? m_horizWeights[pix-1] :  1;
                eng = eng + m_V(m_answer[pix],m_answer[pix-1])*weight;
              }
            for ( y = 1; y < m_height; y++ )
              for ( x = 0; x < m_width; x++ )
              {
                pix = x+y*m_width;
                weight = m_varWeights ? m_vertWeights[pix-m_width] :  1;
                eng = eng + m_V(m_answer[pix],m_answer[pix-m_width])*weight;
              }
        }
        else
        {
          for ( y = 0; y < m_height; y++ )
            for ( x = 1; x < m_width; x++ )
            {
              pix = x+y*m_width;
              eng = eng + m_smoothFn(pix,pix-1,m_answer[pix],m_answer[pix-1]);
            }
          for ( y = 1; y < m_height; y++ )
            for ( x = 0; x < m_width; x++ )
            {
              pix = x+y*m_width;
              eng = eng + m_smoothFn(pix,pix-m_width,m_answer[pix],m_answer[pix-m_width]);
            }
        }
    }
    else
    {
      assert(0);
    }
    return(eng);
}


MRFmodel::EnergyVal MaxProdBP::dataEnergy()
{
    EnergyVal eng = (EnergyVal) 0;
    if ( m_dataType == ARRAY) 
    {
        for ( int i = 0; i < m_nPixels; i++ )
            eng = eng + m_D(i,m_answer[i]);
    }
    else
    {
        for ( int i = 0; i < m_nPixels; i++ )
            eng = eng + m_dataFn(i,m_answer[i]);
    }
    return(eng);
}


void MaxProdBP::setData(DataCostFn dcost)
{
    m_dataFn = dcost;
}


void MaxProdBP::setData(CostVal* data)
{
    m_D = data;
    m_ExpData = new float[m_nPixels * m_nLabels];
    if(!m_ExpData)
    {
      exit(0);
    }
    CostVal cmax = -10;
    int i;
    for (i= 0; i < m_nPixels; i++)
    {
      for(int j = 0; j < m_nLabels; j++)
      {
        if(m_D(i,j) > cmax)
          cmax = m_D(i,j);
      }
    }
    m_exp_scale = float(float(cmax)*4.0);
    float *cData = m_ExpData;
    for (i= 0; i < m_nPixels; i++)
    {
      nodeArray[i].localEv = cData;
      for(int j = 0; j < m_nLabels; j++)
      {
        *cData = float(exp(-1.0*m_D(i,j)/m_exp_scale));
//      printf("%e %e %e| \n",float(m_D(i,j)), *cData, m_exp_scale);
        cData++;
      }
    }
}


void MaxProdBP::setSmoothness(SmoothCostGeneralFn cost)
{
    m_smoothFn = cost;
}


void MaxProdBP::setSmoothness(CostVal* V)
{
    m_V = V;
    m_ExpData = new float[m_nLabels * m_nLabels];
    float *cptr = m_ExpData;
    for (int i = 0; i < m_nLabels; i++)
    {
      for(int j = 0; j < m_nLabels; j++)
      {
        *cptr = float(exp(-1.0*m_V(i,j)));
//      printf("%f %f | ", *cptr, m_V(i,j));
        cptr++;
      }
    }
    printf("\n");
}


void MaxProdBP::setSmoothness(int smoothExp,CostVal smoothMax, CostVal lambda)
{
    int i, j;
    CostVal cost;
    m_needToFreeV = 1;
    m_V = (CostVal *) new CostVal[m_nLabels*m_nLabels*sizeof(CostVal)];
    if (!m_V) { fprintf(stderr, "Not enough memory!\n"); exit(1); }
    for (i=0; i<m_nLabels; i++)
        for (j=i; j<m_nLabels; j++)
        {
			cost = MRFmodel::CostVal((smoothExp == 1) ? j - i : (j - i)*(j - i));
            if (cost > smoothMax) cost = smoothMax;
            m_V[i*m_nLabels + j] = m_V[j*m_nLabels + i] = cost*lambda;
        }
    m_ExpData = new float[m_nLabels * m_nLabels];
    float *cptr = m_ExpData;
//  m_exp_scale = lambda;
    for ( i = 0; i < m_nLabels; i++)
    {
      for( j = 0; j < m_nLabels; j++)
      {
//      printf("%e  ",  float(m_V(i,j)));
        *cptr = float(expf(float(-1.0*float(m_V(i,j))/m_exp_scale)));
//      printf("%e %i | ", *cptr, m_V(i,j));
        cptr++;
      }
    }
//  printf("SmoothMax: %e lambda: %e\n",float(smoothMax),float(lambda));
}


void MaxProdBP::setCues(CostVal* hCue, CostVal* vCue)
{
    m_horizWeights = hCue;
    m_vertWeights  = vCue;
}


void MaxProdBP::optimizeAlg(int nIterations)
{
    //int x, y, i, j, n;
    //Label* l;
    //CostVal* dataPtr;
    CostVal *D = (CostVal *) new CostVal[m_nLabels];
    if ( !D ) {printf("\nNot enough memory, exiting");exit(0);}
    if ( !m_grid_graph) {printf("\nMaxProdBP is not implemented for nongrids yet!");exit(1);}
    int numRows = getHeight();
    int numCols = getWidth();
    const float alpha = float(0.95);
    for (int niter=0; niter < nIterations; niter++)
    {
      for(int r = 0; r < numRows; r++)
      {
        computeMessagesLeftRight(nodeArray, numCols, numRows, r, alpha, this);
      }
      for(int c = 0; c < numCols; c++)
      {
        computeMessagesUpDown(nodeArray, numCols, numRows, c, alpha, this);
      }
      //printf("Iter: %d\n",niter);
    }
      //float tmpBeliefVec[m_nLabels];
     float *tmpBeliefVec = new  float[m_nLabels];
      Label *currAssign = m_answer;
      for(int m = 0; m < numRows; m++)
      {
    for(int n = 0; n < numCols; n++)
    {
      nodeArray[m*numCols+n].getBelief(tmpBeliefVec);
      float currMax = 0;
      int maxInd = -100;
      for (int i = 0; i < m_nLabels; i++)
      {
        if(tmpBeliefVec[i] > currMax)
      {
        currMax = tmpBeliefVec[i];
        maxInd = i;
      }
      }
      currAssign[m * numCols +n] = maxInd;
    }
      }
    delete [] tmpBeliefVec;
}


/***************************************************************************************************************************************/
#define m_D(pix,l)  m_D[(pix)*m_nLabels+(l)]
#define m_V(l1,l2)  m_V[(l1)*m_nLabels+(l2)]


ICM::ICM(int width, int height, int nLabels,EnergyFunction *eng):MRFmodel(width,height,nLabels,eng)
{
    m_needToFreeV = 0;
    initializeAlg();
}


ICM::ICM(int nPixels, int nLabels,EnergyFunction *eng):MRFmodel(nPixels,nLabels,eng)
{
    m_needToFreeV = 0;
    initializeAlg();
}


ICM::~ICM()
{ 
    delete[] m_answer;
    if (!m_grid_graph) delete[] m_neighbors;
    if ( m_needToFreeV ) delete[] m_V;
}


void ICM::initializeAlg()
{
    m_answer = (Label *) new Label[m_nPixels];
    if ( !m_answer ){printf("\nNot enough memory, exiting");exit(0);}
    if (!m_grid_graph)
    {
        m_neighbors = (LinkedBlockList *) new LinkedBlockList[m_nPixels];
        if (!m_neighbors) {printf("Not enough memory,exiting");exit(0);};
    }
}


void ICM::clearAnswer()
{
    memset(m_answer, 0, m_nPixels*sizeof(Label));
}


void ICM::setNeighbors(int pixel1, int pixel2, CostVal weight)
{
    assert(!m_grid_graph);
    assert(pixel1 < m_nPixels && pixel1 >= 0 && pixel2 < m_nPixels && pixel2 >= 0);
    Neighbor *temp1 = (Neighbor *) new Neighbor;
    Neighbor *temp2 = (Neighbor *) new Neighbor;
    if ( !temp1 || ! temp2 ) {printf("\nNot enough memory, exiting");exit(0);}
    temp1->weight  = weight;
    temp1->to_node = pixel2;
    temp2->weight  = weight;
    temp2->to_node = pixel1;
    m_neighbors[pixel1].addFront(temp1);
    m_neighbors[pixel2].addFront(temp2);
}


/*void ICM::clearNeighbors()
{
	m_neighbors->clearList();
}*/


MRFmodel::EnergyVal ICM::smoothnessEnergy()
{
    EnergyVal eng = (EnergyVal) 0;
    EnergyVal weight;
    int x,y,pix,i;
    if ( m_grid_graph )
    {
        if ( m_smoothType != FUNCTION  )
        {
            for ( y = 0; y < m_height; y++ )
                for ( x = 1; x < m_width; x++ )
                {
                    pix    = x+y*m_width;
                    weight = m_varWeights ? m_horizWeights[pix-1] :  1;
                    eng = eng + m_V(m_answer[pix],m_answer[pix-1])*weight;
                }
            for ( y = 1; y < m_height; y++ )
                for ( x = 0; x < m_width; x++ )
                {
                    pix = x+y*m_width;
                    weight = m_varWeights ? m_vertWeights[pix-m_width] :  1;
                    eng = eng + m_V(m_answer[pix],m_answer[pix-m_width])*weight;
                }
        }
        else
        {
            for ( y = 0; y < m_height; y++ )
                for ( x = 1; x < m_width; x++ )
                {
                    pix = x+y*m_width;
                    eng = eng + m_smoothFn(pix,pix-1,m_answer[pix],m_answer[pix-1]);
                }
            for ( y = 1; y < m_height; y++ )
                for ( x = 0; x < m_width; x++ )
                {
                    pix = x+y*m_width;
                    eng = eng + m_smoothFn(pix,pix-m_width,m_answer[pix],m_answer[pix-m_width]);
                }
        }
    }
    else
    {
        Neighbor *temp; 
        if ( m_smoothType != FUNCTION  )
        {
            for ( i = 0; i < m_nPixels; i++ )
                if ( !m_neighbors[i].isEmpty() )
                {
                    m_neighbors[i].setCursorFront();
                    while ( m_neighbors[i].hasNext() )
                    {
                        temp = (Neighbor *) m_neighbors[i].next();
                        if ( i < temp->to_node )
                            eng = eng + m_V(m_answer[i],m_answer[temp->to_node])*(temp->weight);
                    }
                }
        }
        else
        {
            for ( i = 0; i < m_nPixels; i++ )
                if ( !m_neighbors[i].isEmpty() )
                {
                    m_neighbors[i].setCursorFront();
                    while ( m_neighbors[i].hasNext() )
                    {
                        temp = (Neighbor *) m_neighbors[i].next();
                        if ( i < temp->to_node )
                            eng = eng + m_smoothFn(i,temp->to_node, m_answer[i],m_answer[temp->to_node]);
                }
            }
        }
        
    }
    return(eng);
}


MRFmodel::EnergyVal ICM::dataEnergy()
{
    EnergyVal eng = (EnergyVal) 0;
    if ( m_dataType == ARRAY) 
    {
        for ( int i = 0; i < m_nPixels; i++ )
            eng = eng + m_D(i,m_answer[i]);
    }
    else
    {
        for ( int i = 0; i < m_nPixels; i++ )
            eng = eng + m_dataFn(i,m_answer[i]);
    }
    return(eng);
}


void ICM::setData(DataCostFn dcost)
{
    m_dataFn = dcost;
}


void ICM::setData(CostVal* data)
{
    m_D = data;
}


void ICM::setSmoothness(SmoothCostGeneralFn cost)
{
    m_smoothFn = cost;
}


void ICM::setSmoothness(CostVal* V)
{
    m_V = V;
}


void ICM::setSmoothness(int smoothExp,CostVal smoothMax, CostVal lambda)
{
    int i, j;
    CostVal cost;
    m_needToFreeV = 1;
    m_V = (CostVal *) new CostVal[m_nLabels*m_nLabels*sizeof(CostVal)];
    if (!m_V) { fprintf(stderr, "Not enough memory!\n"); exit(1); }
    for (i=0; i<m_nLabels; i++)
        for (j=i; j<m_nLabels; j++)
        {
            cost = (CostVal) ((smoothExp == 1) ? j - i : (j - i)*(j - i));
            if (cost > smoothMax) cost = smoothMax;
            m_V[i*m_nLabels + j] = m_V[j*m_nLabels + i] = cost*lambda;
        }
}


void ICM::setCues(CostVal* hCue, CostVal* vCue)
{
    m_horizWeights = hCue;
    m_vertWeights  = vCue;
}


void ICM::optimizeAlg(int nIterations)
{
    int x, y, i, j, n;
    Label* l;
    CostVal* dataPtr;
    CostVal *D = (CostVal *) new CostVal[m_nLabels];
    if ( !D ) {printf("\nNot enough memory, exiting");exit(0);}
    if ( !m_grid_graph) {printf("\nICM is not implemented for nongrids yet!");exit(1);}
    for ( ; nIterations > 0; nIterations --)
    {
        n = 0;
        l = m_answer;
        dataPtr = m_D;
        for (y=0; y<m_height; y++)
        for (x=0; x<m_width; x++, l++, dataPtr+=m_nLabels, n++)
        {
            // set array D
            if (m_dataType == FUNCTION)
            {
                for (i=0; i<m_nLabels; i++)
                {
                    D[i] = m_dataFn(x+y*m_width, i);
                }
            }
            else memcpy(D, dataPtr, m_nLabels*sizeof(CostVal));
            // add smoothness costs
            if (m_smoothType == FUNCTION)
            {
                if (x > 0)
                {
                    j = *(l-1);
                    for (i=0; i<m_nLabels; i++) D[i] += m_smoothFn(x+y*m_width-1, x+y*m_width, j, i);
                }
                if (y > 0)
                {
                    j = *(l-m_width);
                    for (i=0; i<m_nLabels; i++) D[i] += m_smoothFn(x+y*m_width-m_width,x+y*m_width , j, i);
                }
                if (x < m_width-1)
                {
                    j = *(l+1);
                    for (i=0; i<m_nLabels; i++) D[i] += m_smoothFn(x+y*m_width+1, x+y*m_width, i, j);
                }
                if (y < m_height-1)
                {
                    j = *(l+m_width);
                    for (i=0; i<m_nLabels; i++) D[i] += m_smoothFn(x+y*m_width+m_width, x+y*m_width, i, j);
                }
            }
            else
            {
                if (x > 0)
                {
                    j = *(l-1);
                    CostVal lambda = (m_varWeights) ? m_horizWeights[n-1] : 1;
                    for (i=0; i<m_nLabels; i++) D[i] += lambda * m_V[j*m_nLabels + i];
                }
                if (y > 0)
                {
                    j = *(l-m_width);
                    CostVal lambda = (m_varWeights) ? m_vertWeights[n-m_width] : 1;
                    for (i=0; i<m_nLabels; i++) D[i] += lambda * m_V[j*m_nLabels + i];
                }
                if (x < m_width-1)
                {
                    j = *(l+1);
                    CostVal lambda = (m_varWeights) ? m_horizWeights[n] : 1;
                    for (i=0; i<m_nLabels; i++) D[i] += lambda * m_V[j*m_nLabels + i];
                }
                if (y < m_height-1)
                {
                    j = *(l+m_width);
                    CostVal lambda = (m_varWeights) ? m_vertWeights[n] : 1;
                    for (i=0; i<m_nLabels; i++) D[i] += lambda * m_V[j*m_nLabels + i];
                }
            }
            // compute minimum of D, set new label for (x,y)
            CostVal D_min = D[0];
            *l = 0;
            for (i=1; i<m_nLabels; i++)
            {
                if (D_min > D[i])
                {
                    D_min = D[i];
                    *l = i;
                }
            }
        }
    }
    delete[] D;
}
