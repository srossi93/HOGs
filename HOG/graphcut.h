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


#ifndef __MRF_H__
#define __MRF_H__
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>


enum NeighborhoodOrder { FIRST_ORDER, SECOND_ORDER};
const NeighborhoodOrder default_graph_cut_order = SECOND_ORDER;


class EnergyFunction;

class MRFmodel
{
public:
    // *********** CONSTRUCTORS/DESTRUCTOR
    // Constructor. After you call this, you must call setData and setSmoothness            
    // Use this constructor for 2D grid graphs of size width by height Standard 4-connected 
    // neighborhood system is assumed. Labels are in the range 0,1,...nLabels - 1           
    // Width is in  the range 0,1,...width-1 and height is in the range 0,1,...height-1     
    // Input parameter eng specifies the data and smoothness parts of the energy
    // For 2D grids, since 4 connected neighborhood structure is assumed, this 
    // fully specifies the energy
    MRFmodel(int width, int height, int nLabels, EnergyFunction *eng);
    // Use this constructor for a general neighborhood system. Pixels are in the range      
    // 0,1,..nPixels-1, and labels are in the range 0,1,...,nLabels-1 
    // Input parameter eng specifies the data and smoothness parts of the energy
    // after this constructor you need to call setNeighbors() to specify the heighborhood system
    MRFmodel(int nPixels, int nLabels, EnergyFunction *eng);
    virtual ~MRFmodel() { }
    // Returns true if energy function has been specified, returns false otherwise 
    // By default, it always returns true. Can be modified by the suppler of       
    // optimization algorithm                                                      
    virtual int isValid(){return true;};  
    // *********** EVALUATING THE ENERGY
    typedef int Label;
    typedef float EnergyVal;        /* The total energy of a labeling */
    typedef float CostVal;          /* costs of individual terms of the energy */
    EnergyVal totalEnergy();      /* returns energy of current labeling */
    virtual EnergyVal dataEnergy() = 0;        /* returns the data part of the energy */
    virtual EnergyVal smoothnessEnergy() = 0;  /* returns the smoothness part of the energy */
    //Functional representation for data costs
    typedef CostVal (*DataCostFn)(int pix, Label l); 
    // Functional representation for the general cost function type 
    typedef CostVal (*SmoothCostGeneralFn)(int pix1, int pix2,  Label l1, Label l2); 
    // Use this function only for non-grid graphs. Sets pix1 and pix2 to be neighbors 
    // with the specified weight.  Can be called ONLY once for each pair of pixels    
    // That is if pixel1 and pixel2 are neihbors, call  either setNeighbors(pixel1,pixel2,weight) 
    // or setNeighbors(pixel2,pixel1,weight), but NOT BOTH                                     
    virtual void setNeighbors(int pix1, int pix2, CostVal weight)= 0;
    //virtual void clearNeighbors()= 0;
    void initialize();
    // Runs optimization for nIterations. Input parameter time returns the time it took to
    // perform nIterations of optimization
    void optimize(int nIterations, float& time);
    void optimize(int nIterations);
    virtual void optimizeAlg(int nIterations)=0;
    // *********** ACCESS TO SOLUTION
    // Returns pointer to array of size nPixels. Client may then read/write solution (but not deallocate array).
    virtual Label* getAnswerPtr()= 0;
    // returns the label of the input pixel
    virtual Label getLabel(int pixel)= 0;
    // sets label of a pixel
    virtual void setLabel(int pixel,Label label)= 0;
    // sets all the labels to zero
    virtual void clearAnswer()  = 0;
    // use this function to pass any parameters to optimization algorithm. 
    // The first argument is the number of passed, parameters  and
    // the second argument is the pointer to the array of parameters
    virtual void setParameters(int numParam, void *param) = 0;
    // This function returns lower bound computed by the algorithm (if any)
    // By default, it returns 0.
    virtual EnergyVal lowerBound(){return((EnergyVal) 0);};
    // Returns 0 if the energy is not suitable for current optimization algorithm 
    // Returns 1 if the energy is suitable for current optimization algorithm
    // Returns 2 if current optimizaiton algorithm does not check the energy 
    virtual char checkEnergy();
    typedef enum
        {
            FUNCTION,
            ARRAY,
            THREE_PARAM,
            NONE
        } InputType;
protected:
    int  m_width, m_height;  // width and height of a grid,if graph is a grid
    int  m_nPixels;          // number of pixels, for both grid and non-grid graphs
    int  m_nLabels;          // number of labels, for both grid and non-grid graphs
    bool m_grid_graph;   // true if the graph is a 2D grid
	bool m_second_order;
    bool m_varWeights;   // true if weights are spatially varying. To be used only with 2D grids
    bool m_initialized;  // true if array m_V is allocated memory.  
    EnergyFunction *m_e;
    InputType m_dataType;     
    InputType m_smoothType;
    // *********** SET THE DATA COSTS 
    // Following 2 functions set the data costs
    virtual void setData(DataCostFn dcost)=0; 
    virtual void setData(CostVal* data)=0;   
    // *********** SET THE SMOOTHNESS COSTS 
    // following 3 functions set the smoothness costs 
    // there are 2 ways to represent the smoothness costs, one with array, one with function
    // In addition, for 2D grid graphs spacially varying weights can be specified by 2 arrays
    // Smoothness cost depends on labels  V(l1,l2) for all edges (except for a multiplier - see setCues() ).
    // V must be symmetric: V(l1,l2) = V(l2,l1)
    // V must be an array of size nLabels*nLabels. It is NOT copied into internal memory
    virtual void setSmoothness(CostVal* V)=0;
    // General smoothness cost can be specified by passing pointer to a function 
    virtual  void setSmoothness(SmoothCostGeneralFn cost)=0;
    // Use if the smoothness is V(l1,l2) = lambda * min ( |l1-l2|^m_smoothExp, m_smoothMax )
    // Can also add optional spatially varying weights for 2D grid graphs using setCues()
    virtual void setSmoothness(int smoothExp,CostVal smoothMax, CostVal lambda)=0;
    // You are not required to call setCues, in which case there is no multiplier.
    // Function below cannot be called for general cost function.
    // This function can be only used for a 2D grid graph
    // hCue and vCue must be arrays of size width*height in row major order. 
    // They are NOT copied into internal memory.
    // hCue(x,y) holds the variable weight for edge between pixels (x+1,y) and (x,y)
    // vCue(x,y) holds the variable weight for edge between pixels (x,y+1) and (x,y)
    virtual void setCues(CostVal* hCue, CostVal* vCue)=0; 
    void commonInitialization(EnergyFunction *e);
    void checkArray(CostVal *V);
};


// *********** This class is for data costs
// Data costs can be specified eithe by an array or by a pointer to a function 
// If specified by an array, use constructor DataCost(cost) where 
// cost is the array of type CostVal. The cost of pixel p and label l is
// stored at cost[p*nLabels+l] where nLabels is the number of labels
// If data costs are to be specified by a function, pass
// a pointer to a function 
// CostVal costFn(int pix, Label lab)
// which returns the
// data cost of pixel pix to be assigned label lab
class DataCost
{
friend class MRFmodel;
public:
    typedef MRFmodel::CostVal CostVal;
    typedef MRFmodel::DataCostFn DataCostFn;
    DataCost(CostVal *cost){m_costArray = cost;m_type = MRFmodel::ARRAY; };
    DataCost(DataCostFn costFn){m_costFn = costFn;m_type = MRFmodel::FUNCTION;};
private:
    MRFmodel::CostVal *m_costArray;
    MRFmodel::DataCostFn m_costFn;
    MRFmodel::InputType m_type;     
};


// ***************** This class represents smoothness costs 
// If the smoothness is V(l1,l2) = lambda * min ( |l1-l2|^m_smoothExp, m_smoothMax )
// use constructor SmoothnessCost(smoothExp,smoothMax,lambda)
// If, in addition,  there are spacially varying weights use constructor
// SmoothnessCost(smoothExp,smoothMax,lambda,hWeights,vWeights)
// hWeights and vWeights can be only used for a 2D grid graph
// hWeights and vWeights must be arrays of size width*height in row major order. 
// They are NOT copied into internal memory.
// hWeights(x,y) holds the variable weight for edge between pixels (x+1,y) and (x,y)
// vWeights(x,y) holds the variable weight for edge between pixels (x,y+1) and (x,y)
//  If the smoothness costs are specified by input array V of type CostVal and
// size nLabels*nLabels, use consructor SmoothnessCost(V). 
// If in addition, there are 
// are spacially varying weights use constructor SmoothnessCost(V,hWeights,vWeights)
// Note that array V must be of size nLabels*nLabels, and be symmetric. 
// That is V[i*nLabels+j] = V[j*nLabels+i]
// Finally, if the smoothness term is specified by a general function, use
// constructor SmoothnessCost(costFn)
class SmoothnessCost
{
    friend class MRFmodel;
public:
    typedef MRFmodel::CostVal CostVal;
    // Can be used for  2D grids and for general graphs
    // In case if used for 2D grids, the smoothness term WILL NOT be spacially varying
    SmoothnessCost(int smoothExp,CostVal smoothMax,CostVal lambda)
        {m_type=MRFmodel::THREE_PARAM;m_smoothMax = smoothMax;m_smoothExp = smoothExp;m_lambda=lambda;m_varWeights=false;};
    // Can be used only for 2D grids 
    // the smoothness term WILL BE be spacially varying
    SmoothnessCost(int smoothExp,CostVal smoothMax,CostVal lambda,CostVal *hWeights, CostVal *vWeights)
        {m_type=MRFmodel::THREE_PARAM;m_smoothMax = smoothMax;m_smoothExp = smoothExp;m_lambda=lambda;
        m_varWeights = true;m_hWeights = hWeights; m_vWeights = vWeights;};
    // Can be used 2D grids and for general graphs
    // In case if used for 2D grids, the smoothness term WILL NOT be spacially varying
    SmoothnessCost(CostVal *V){m_V = V;m_type = MRFmodel::ARRAY;m_varWeights=false;};
    // Can be used only for 2D grids 
    // the smoothness term WILL BE be spacially varying
    SmoothnessCost(CostVal *V,CostVal *hWeights, CostVal *vWeights )
        {m_V = V;m_hWeights = hWeights; m_vWeights = vWeights; m_varWeights = true; m_type=MRFmodel::ARRAY;};
    // Can be used 2D grids and for general graphs
    SmoothnessCost(MRFmodel::SmoothCostGeneralFn costFn){m_costFn = costFn;m_type = MRFmodel::FUNCTION;};
private:
    CostVal *m_V,*m_hWeights, *m_vWeights;
    MRFmodel::SmoothCostGeneralFn m_costFn;
    MRFmodel::InputType m_type;
    int m_smoothExp;
    CostVal m_smoothMax,m_lambda;
    bool m_varWeights;
    EnergyFunction *m_eng;
};


class EnergyFunction
{
public:
    EnergyFunction(DataCost *dataCost,SmoothnessCost *smoothCost)
        {m_dataCost = dataCost;m_smoothCost = smoothCost;};
    EnergyFunction(DataCost *dataCost,SmoothnessCost *smoothCost, NeighborhoodOrder order)
        {m_dataCost = dataCost;m_smoothCost = smoothCost;m_order=order;};
    DataCost *m_dataCost;
    SmoothnessCost *m_smoothCost;
	NeighborhoodOrder m_order;
};


#endif /*  __MRF_H__ */

/*
    virtual EnergyVal dataEnergy() = 0;        
    virtual EnergyVal smoothnessEnergy() = 0;  
    virtual void setNeighbors(int pix1, int pix2, CostVal weight)= 0;
    virtual void optimizeAlg(int nIterations)=0;
    virtual Label* getAnswerPtr()= 0;
    virtual Label getLabel(int pixel)= 0;
    virtual void setLabel(int pixel,Label label)= 0;
    virtual void clearAnswer()  = 0;
    virtual void setParameters(int numParam, void *param) = 0;
    virtual void setData(DataCostFn dcost)=0; 
    virtual void setData(CostVal* data)=0;   
    virtual void setSmoothness(CostVal* V)=0;
    virtual void setSmoothness(SmoothCostGeneralFn cost)=0;
    virtual void setCues(CostVal* hCue, CostVal* vCue)=0; 
    virtual void setSmoothness(int smoothExp,CostVal smoothMax, CostVal lambda);
    virtual EnergyVal lowerBound(){return((EnergyVal) 0);};
*/
/***************************************************************************************************************************************/
#ifndef __BLOCK_H__
#define __BLOCK_H__

template <class Type> class Block
{
public:
    /* Constructor. Arguments are the block size and
       (optionally) the pointer to the function which
       will be called if allocation failed; the message
       passed to this function is "Not enough memory!" */
    Block(int size, void (*err_function)(char *) = NULL) { first = last = NULL; block_size = size; error_function = err_function; }
    /* Destructor. Deallocates all items added so far */
    ~Block() { while (first) { block *next = first -> next; delete first; first = next; } }
    /* Allocates 'num' consecutive items; returns pointer
       to the first item. 'num' cannot be greater than the
       block size since items must fit in one block */
    Type *New(int num = 1)
    {
        Type *t;

        if (!last || last->current + num > last->last)
        {
            if (last && last->next) last = last -> next;
            else
            {
                block *next = (block *) new char [sizeof(block) + (block_size-1)*sizeof(Type)];
                if (!next) { if (error_function) (*error_function)("Not enough memory!"); exit(1); }
                if (last) last -> next = next;
                else first = next;
                last = next;
                last -> current = & ( last -> data[0] );
                last -> last = last -> current + block_size;
                last -> next = NULL;
            }
        }
        t = last -> current;
        last -> current += num;
        return t;
    }
    /* Returns the first item (or NULL, if no items were added) */
    Type *ScanFirst()
    {
        scan_current_block = first;
        if (!scan_current_block) return NULL;
        scan_current_data = & ( scan_current_block -> data[0] );
        return scan_current_data ++;
    }
    /* Returns the next item (or NULL, if all items have been read)
       Can be called only if previous ScanFirst() or ScanNext()
       call returned not NULL. */
    Type *ScanNext()
    {
        if (scan_current_data >= scan_current_block -> current)
        {
            scan_current_block = scan_current_block -> next;
            if (!scan_current_block) return NULL;
            scan_current_data = & ( scan_current_block -> data[0] );
        }
        return scan_current_data ++;
    }
    /* Marks all elements as empty */
    void Reset()
    {
        block *b;
        if (!first) return;
        for (b=first; ; b=b->next)
        {
            b -> current = & ( b -> data[0] );
            if (b == last) break;
        }
        last = first;
    }

private:
    typedef struct block_st
    {
        Type                    *current, *last;
        struct block_st         *next;
        Type                    data[1];
    } block;
    int     block_size;
    block   *first;
    block   *last;
    block   *scan_current_block;
    Type    *scan_current_data;
    void    (*error_function)(char *);
};


template <class Type> class DBlock
{
public:
    /* Constructor. Arguments are the block size and
       (optionally) the pointer to the function which
       will be called if allocation failed; the message
       passed to this function is "Not enough memory!" */
    DBlock(int size, void (*err_function)(char *) = NULL) { first = NULL; first_free = NULL; block_size = size; error_function = err_function; }
    /* Destructor. Deallocates all items added so far */
    ~DBlock() { while (first) { block *next = first -> next; delete first; first = next; } }
    /* Allocates one item */
    Type *New()
    {
        block_item *item;
        if (!first_free)
        {
            block *next = first;
            first = (block *) new char [sizeof(block) + (block_size-1)*sizeof(block_item)];
            if (!first) { if (error_function) (*error_function)("Not enough memory!"); exit(1); }
            first_free = & (first -> data[0] );
            for (item=first_free; item<first_free+block_size-1; item++)
                item -> next_free = item + 1;
            item -> next_free = NULL;
            first -> next = next;
        }
        item = first_free;
        first_free = item -> next_free;
        return (Type *) item;
    }
    /* Deletes an item allocated previously */
    void Delete(Type *t)
    {
        ((block_item *) t) -> next_free = first_free;
        first_free = (block_item *) t;
    }

private:
    typedef union block_item_st
    {
        Type            t;
        block_item_st   *next_free;
    } block_item;
    typedef struct block_st
    {
        struct block_st         *next;
        block_item              data[1];
    } block;
    int         block_size;
    block       *first;
    block_item  *first_free;
    void    (*error_function)(char *);
};

#endif


/***************************************************************************************************************************************/
#ifndef __GRAPH_H__
#define __GRAPH_H__
/*
    Nodes, arcs and pointers to nodes are
    added in blocks for memory and time efficiency.
    Below are numbers of items in blocks
*/
#define NODE_BLOCK_SIZE 512
#define ARC_BLOCK_SIZE 1024
#define NODEPTR_BLOCK_SIZE 128

class Graph
{
public:
    typedef enum
    {
        SOURCE  = 0,
        SINK    = 1
    } termtype; /* terminals */
    typedef MRFmodel::CostVal captype;
    /* Type of total flow */
    typedef MRFmodel::EnergyVal flowtype;
    typedef void * node_id;
    /* interface functions */
    /* Constructor. Optional argument is the pointer to the
       function which will be called if an error occurs;
       an error message is passed to this function. If this
       argument is omitted, exit(1) will be called. */
    Graph(void (*err_function)(char *) = NULL);
    /* Destructor */
    ~Graph();
    /* Adds a node to the graph */
    node_id add_node();
    /* Adds a bidirectional edge between 'from' and 'to'
       with the weights 'cap' and 'rev_cap' */
    void add_edge(node_id from, node_id to, captype cap, captype rev_cap);
    /* Sets the weights of the edges 'SOURCE->i' and 'i->SINK'
       Can be called at most once for each node before any call to 'add_tweights'.
       Weights can be negative */
    void set_tweights(node_id i, captype cap_source, captype cap_sink);
    /* Adds new edges 'SOURCE->i' and 'i->SINK' with corresponding weights
       Can be called multiple times for each node.
       Weights can be negative */
    void add_tweights(node_id i, captype cap_source, captype cap_sink);
    /* After the maxflow is computed, this function returns to which
       segment the node 'i' belongs (Graph::SOURCE or Graph::SINK) */
    termtype what_segment(node_id i);
    /* Computes the maxflow. Can be called only once. */
    flowtype maxflow();

private:
    /* internal variables and functions */
    struct arc_forward_st;
    struct arc_reverse_st;

#define IS_ODD(a) ((int)(a) & 1)
#define MAKE_ODD(a)  ((arc_forward *) ((int)(a) | 1))
#define MAKE_EVEN(a) ((arc_forward *) ((int)(a) & (~1)))
#define MAKE_ODD_REV(a)  ((arc_reverse *) ((int)(a) | 1))
#define MAKE_EVEN_REV(a) ((arc_reverse *) ((int)(a) & (~1)))

    /* node structure */
    typedef struct node_st
    {
        /*
            Usually i->first_out is the first outgoing
            arc, and (i+1)->first_out-1 is the last outgoing arc.
            However, it is not always possible, since
            arcs are allocated in blocks, so arcs corresponding
            to two consecutive nodes may be in different blocks.
            If outgoing arcs for i are last in the arc block,
            then a different mechanism is used. i->first_out
            is odd in this case; the first outgoing arc
            is (a+1), and the last outgoing arc is
            ((arc_forward *)(a->shift))-1, where
            a = (arc_forward *) (((char *)(i->first_out)) + 1);
            Similar mechanism is used for incoming arcs.
        */
        arc_forward_st  *first_out; /* first outcoming arc */
        arc_reverse_st  *first_in;  /* first incoming arc */
        arc_forward_st  *parent;    /* describes node's parent
                                       if IS_ODD(parent) then MAKE_EVEN(parent) points to 'arc_reverse',
                                       otherwise parent points to 'arc_forward' */
        node_st         *next;      /* pointer to the next active node
                                       (or to itself if it is the last node in the list) */
        int             TS;         /* timestamp showing when DIST was computed */
        int             DIST;       /* distance to the terminal */
        short           is_sink;    /* flag showing whether the node is in the source or in the sink tree */
        captype         tr_cap;     /* if tr_cap > 0 then tr_cap is residual capacity of the arc SOURCE->node
                                       otherwise         -tr_cap is residual capacity of the arc node->SINK */
    } node;
    /* arc structures */
#define NEIGHBOR_NODE(i, shift) ((node *) ((char *)(i) + (shift)))
#define NEIGHBOR_NODE_REV(i, shift) ((node *) ((char *)(i) - (shift)))
    typedef struct arc_forward_st
    {
        int             shift;      /* node_to = NEIGHBOR_NODE(node_from, shift) */
        captype         r_cap;      /* residual capacity */
        captype         r_rev_cap;  /* residual capacity of the reverse arc*/
    } arc_forward;
    typedef struct arc_reverse_st
    {
        arc_forward     *sister;    /* reverse arc */
    } arc_reverse;
    /* 'pointer to node' structure */
    typedef struct nodeptr_st
    {
        node_st         *ptr;
        nodeptr_st      *next;
    } nodeptr;
    typedef struct node_block_st
    {
        node                    *current;
        struct node_block_st    *next;
        node                    nodes[NODE_BLOCK_SIZE];
    } node_block;

#define last_node LAST_NODE.LAST_NODE

    typedef struct arc_for_block_st
    {
        char                    *start;     /* the actual start address of this block.
                                               May be different from 'this' since 'this'
                                               must be at an even address. */
        arc_forward             *current;
        struct arc_for_block_st *next;
        arc_forward             arcs_for[ARC_BLOCK_SIZE]; /* all arcs must be at even addresses */
        union
        {
            arc_forward         dummy;
            node                *LAST_NODE; /* used in graph consruction */
        }                       LAST_NODE;
    } arc_for_block;

    typedef struct arc_rev_block_st
    {
        char                    *start;     /* the actual start address of this block.
                                               May be different from 'this' since 'this'
                                               must be at an even address. */
        arc_reverse             *current;
        struct arc_rev_block_st *next;
        arc_reverse             arcs_rev[ARC_BLOCK_SIZE]; /* all arcs must be at even addresses */
        union
        {
            arc_reverse         dummy;
            node                *LAST_NODE; /* used in graph consruction */
        }                       LAST_NODE;
    } arc_rev_block;

    node_block          *node_block_first;
    arc_for_block       *arc_for_block_first;
    arc_rev_block       *arc_rev_block_first;
    DBlock<nodeptr>     *nodeptr_block;
    void    (*error_function)(char *);  /* this function is called if a error occurs,
                                           with a corresponding error message
                                           (or exit(1) is called if it's NULL) */
    flowtype            flow;       /* total flow */
    node                *queue_first[2], *queue_last[2];    /* list of active nodes */
    nodeptr             *orphan_first, *orphan_last;        /* list of pointers to orphans */
    int                 TIME;                               /* monotonically increasing global counter */
    /* functions for processing active list */
    void set_active(node *i);
    node *next_active();
    void prepare_graph();
    void maxflow_init();
    void augment(node *s_start, node *t_start, captype *cap_middle, captype *rev_cap_middle);
    void process_source_orphan(node *i);
    void process_sink_orphan(node *i);
};

#endif

/***************************************************************************************************************************************/
class Energy : Graph
{
public:
    typedef node_id Var;

    /* Types of energy values.
       Value is a type of a value in a single term
       TotalValue is a type of a value of the total energy.
       By default Value = short, TotalValue = int.
       To change it, change the corresponding types in graph.h */
    typedef captype Value;
    typedef flowtype TotalValue;
    /* interface functions */
    /* Constructor. Optional argument is the pointer to the
       function which will be called if an error occurs;
       an error message is passed to this function. If this
       argument is omitted, exit(1) will be called. */
    Energy(void (*err_function)(char *) = NULL);
    /* Destructor */
    ~Energy();
    /* Adds a new binary variable */
    Var add_variable();
    /* Adds a constant E to the energy function */
    void add_constant(Value E);
    /* Adds a new term E(x) of one binary variable
       to the energy function, where
           E(0) = E0, E(1) = E1
       E0 and E1 can be arbitrary */
    void add_term1(Var x,
                   Value E0, Value E1);
    /* Adds a new term E(x,y) of two binary variables
       to the energy function, where
           E(0,0) = E00, E(0,1) = E01
           E(1,0) = E10, E(1,1) = E11
       The term must be submodular, i.e. E00 + E11 <= E01 + E10 */
    void add_term2(Var x, Var y,
                   Value E00, Value E01,
                   Value E10, Value E11);
    /* Adds a new term E(x,y,z) of three binary variables
       to the energy function, where
           E(0,0,0) = E000, E(0,0,1) = E001
           E(0,1,0) = E010, E(0,1,1) = E011
           E(1,0,0) = E100, E(1,0,1) = E101
           E(1,1,0) = E110, E(1,1,1) = E111
       The term must be submodular. It means that if one
       of the variables is fixed (for example, y=1), then
       the resulting function of two variables must be submodular.
       Since there are 6 ways to fix one variable
       (3 variables times 2 binary values - 0 and 1),
       this is equivalent to 6 inequalities */
    void add_term3(Var x, Var y, Var z,
                   Value E000, Value E001,
                   Value E010, Value E011,
                   Value E100, Value E101,
                   Value E110, Value E111);
    /* After the energy function has been constructed,
       call this function to minimize it.
       Returns the minimum of the function */
    TotalValue minimize();
    /* After 'minimize' has been called, this function
       can be used to determine the value of variable 'x'
       in the optimal solution.
       Returns either 0 or 1 */
    int get_var(Var x);

private:
    /* internal variables and functions */
    TotalValue  Econst;
    void        (*error_function)(char *);  /* this function is called if a error occurs,
                                            with a corresponding error message
                                            (or exit(1) is called if it's NULL) */
};

/************************  Implementation ******************************/
#ifndef __ENERGY_H__
#define __ENERGY_H__


inline Energy::Energy(void (*err_function)(char *)) : Graph(err_function)
{
    Econst = 0;
    error_function = err_function;
}

inline Energy::~Energy() {}

inline Energy::Var Energy::add_variable() { return add_node(); }

inline void Energy::add_constant(Value A) { Econst += A; }

inline void Energy::add_term1(Var x, Value A, Value B)
{
    add_tweights(x, B, A);
}

inline void Energy::add_term2(Var x, Var y, Value A, Value B, Value C, Value D)
{
    if ( A+D > C+B) 
    {
        Value delta = A+D-C-B;
        Value subtrA = delta/3;

        A = A-subtrA;
        C = C+subtrA;
        B = B+(delta-subtrA*2);
    }
    /* 
       E = A A  +  0   B-A
           D D     C-D 0
       Add edges for the first term
    */
    add_tweights(x, D, A);
    B -= A; C -= D;
    /* now need to represent
       0 B
       C 0
    */
    assert(B + C >= 0); /* check regularity */
    if (B < 0)
    {
        /* Write it as
           B B  +  -B 0  +  0   0
           0 0     -B 0     B+C 0
        */
        add_tweights(x, 0, B); /* first term */
        add_tweights(y, 0, -B); /* second term */
        add_edge(x, y, 0, B+C); /* third term */
    }
    else if (C < 0)
    {
        /* Write it as
           -C -C  +  C 0  +  0 B+C
            0  0     C 0     0 0
        */
        add_tweights(x, 0, -C); /* first term */
        add_tweights(y, 0, C); /* second term */
        add_edge(x, y, B+C, 0); /* third term */
    }
    else /* B >= 0, C >= 0 */
    {
        add_edge(x, y, B, C);
    }
}

inline void Energy::add_term3(Var x, Var y, Var z, Value E000, Value E001, Value E010, Value E011, Value E100, Value E101, Value E110, Value E111)
{
    register Value pi = (E000 + E011 + E101 + E110) - (E100 + E010 + E001 + E111);
    register Value delta;
    register Var u;
    if (pi >= 0)
    {
        Econst += E111 - (E011 + E101 + E110);
        add_tweights(x, E101, E001);
        add_tweights(y, E110, E100);
        add_tweights(z, E011, E010);
        delta = (E010 + E001) - (E000 + E011); /* -pi(E[x=0]) */
        assert(delta >= 0); /* check regularity */
        add_edge(y, z, delta, 0);
        delta = (E100 + E001) - (E000 + E101); /* -pi(E[y=0]) */
        assert(delta >= 0); /* check regularity */
        add_edge(z, x, delta, 0);
        delta = (E100 + E010) - (E000 + E110); /* -pi(E[z=0]) */
        assert(delta >= 0); /* check regularity */
        add_edge(x, y, delta, 0);
        if (pi > 0)
        {
            u = add_variable();
            add_edge(x, u, pi, 0);
            add_edge(y, u, pi, 0);
            add_edge(z, u, pi, 0);
            add_tweights(u, 0, pi);
        }
    }
    else
    {
        Econst += E000 - (E100 + E010 + E001);
        add_tweights(x, E110, E010);
        add_tweights(y, E011, E001);
        add_tweights(z, E101, E100);
        delta = (E110 + E101) - (E100 + E111); /* -pi(E[x=1]) */
        assert(delta >= 0); /* check regularity */
        add_edge(z, y, delta, 0);
        delta = (E110 + E011) - (E010 + E111); /* -pi(E[y=1]) */
        assert(delta >= 0); /* check regularity */
        add_edge(x, z, delta, 0);
        delta = (E101 + E011) - (E001 + E111); /* -pi(E[z=1]) */
        assert(delta >= 0); /* check regularity */
        add_edge(y, x, delta, 0);
        u = add_variable();
        add_edge(u, x, -pi, 0);
        add_edge(u, y, -pi, 0);
        add_edge(u, z, -pi, 0);
        add_tweights(u, -pi, 0);
    }
}

inline Energy::TotalValue Energy::minimize() { return Econst + maxflow(); }

inline int Energy::get_var(Var x) { return (int)what_segment(x); }

#endif


/***************************************************************************************************************************************/
/* Singly Linked List of Blocks */
// This data structure should be used only for the GCoptimization class implementation
// because it lucks some important general functions for general list, like remove_item()
// The head block may be not full 
// For regular 2D grids, it's better to set GCLL_BLOCK_SIZE to 2
// For other graphs, it should be set to the average expected number of neighbors
// Data in linked list for the neighborhood system is allocated in blocks of size GCLL_BLOCK_SIZE 

#ifndef __LINKEDBLOCKLIST_H__
#define __LINKEDBLOCKLIST_H__
#define GCLL_BLOCK_SIZE 4  
// GCLL_BLOCKSIZE should "fit" into the type BlockType. That is 
// if GCLL_BLOCKSIZE is larger than 255 but smaller than largest short integer
// then  BlockType should be set to short
typedef unsigned char BlockType;
//The type of data stored in the linked list
typedef void * ListType;


class LinkedBlockList{
public: 
    void addFront(ListType item);
    inline bool isEmpty(){if (m_head == 0) return(true); else return(false);};
    inline LinkedBlockList(){m_head = 0; m_head_block_size = GCLL_BLOCK_SIZE;}; 
    ~LinkedBlockList();
    // Next three functins are for the linked list traversal
    inline void setCursorFront(){m_cursor = m_head; m_cursor_ind = 0;};
    ListType next();
    bool hasNext();
	//void clearList();

private:
    typedef struct LLBlockStruct{
        ListType m_item[GCLL_BLOCK_SIZE];
        struct LLBlockStruct *m_next;
    } LLBlock;
    LLBlock *m_head;
    // Remembers the number of elements in the head block, since it may not be full
    BlockType m_head_block_size;
    // For block traversal, points to current element in the current block
    BlockType m_cursor_ind;
    // For block traversal, points to current block in the linked list
    LLBlock *m_cursor;
};

#endif


/***************************************************************************************************************************************/
#ifndef __GCOPTIMIZATION_H__
#define __GCOPTIMIZATION_H__

#define m_datacost(pix,lab)     (m_datacost[(pix)*m_nLabels+(lab)] )
#define m_smoothcost(lab1,lab2) (m_smoothcost[(lab1)+(lab2)*m_nLabels] )
#define USE_MEMBER_FUNCTION 0
#define PASS_AS_PARAMETER   1

class GCoptimization:public MRFmodel
{
public:
    ///////////////////////////////////////////////////////////////////////////////////////
    //    First define needed data types                                                 //
    ///////////////////////////////////////////////////////////////////////////////////////
    /* Type for the total energy calculation. Change it in Graph.h   */
    typedef Graph::flowtype EnergyType;
    /* Type for the individual terms in the energy function.Change it in Graph.h */
    typedef Graph::captype EnergyTermType;
    /* Type of label. Can be set to char, short, int, long */
    typedef int LabelType;
    /* Type for pixel. Can be set to  short, int, long */ 
    typedef int PixelType;
    
    ///////////////////////////////////////////////////////////////////////////////////////
    //    Next define constructors                                                       //
    ///////////////////////////////////////////////////////////////////////////////////////
    /* Use this constructor only for grid of size width by height.  If you use this constructor,  */
    /* 4 connected grid neigbhorhood structure is assumed.  num_labels is the number of labels,   */
    /* Labels are assumed to be in  {0, 1, 2,....num_labels-1}                                    */
    /* dataSetup and smoothSetup can take only 2 values, namely USE_MEMBER_FUNCTION (0) and       */
    /* PASS_AS_PARAMETER(1)                                                                       */ 
    /* In case dataSetup = USE_MEMBER_FUNCTION, to set the data cost you must use the             */
    /* member function setDataCost(pixel,label,cost).  If dataSetup = PASS_AS_PARAMETER,          */
    /* to set the data costs you must use any of the setData() functions.                         */
    /* In case smoothSetup = USE_MEMBER_FUNCTION, to set the smooth cost you must use the         */
    /* member function setSmoothCost(label1,label2,cost).  If smoothSetup = PASS_AS_PARAMETER,    */
    /* to set the smoothness costs you must use any of the setSmoothness() functions.             */
    GCoptimization(PixelType width,PixelType height,int num_labels,EnergyFunction *eng);
    /* This is the constructor for non-grid graphs. Everything is the same as in the constructor  */
    /* above, but neighborhood structure must  be specified by any of the two setNeighbors()      */
    /* functions */
    GCoptimization(PixelType num_pixels,int nLabels, EnergyFunction *eng);
    ~GCoptimization();
    /* Returns current label assigned to input pixel */
    inline LabelType getLabel(PixelType pixel){assert(pixel >= 0 && pixel < m_nPixels);return(m_labeling[pixel]);};
    /* Sets data cost of pixel to label */
    /* returns pointer to the answer */
    Label* getAnswerPtr(){return(m_labeling);}
    // Sets answer to zeros
    void clearAnswer();
    /* Makes pixel1 and pixel2 neighbors of each other. Can be called only 1 time for each         */
    /* unordered pair of pixels. Parameter weight can be used to set spacially varying terms       */
    /* If the desired penalty for neighboring pixels pixel1 and pixel2 is                          */
    /*  V(label1,label2) = weight*SmoothnessPenalty(label1,label2), then                           */
    /* member function setLabel should be called as: setLabel(pixel1,pixel2,weight)                */
    void setNeighbors(PixelType pixel1, PixelType pixel2, EnergyTermType weight);
    //void clearNeighbors();
    /* This function can be used to change the label of any pixel at any time      */
    inline void setLabel(PixelType pixel, LabelType label){
        assert(label >= 0 && label < m_nLabels && pixel >= 0 && pixel < m_nPixels);m_labeling[pixel] = label;};
    /* By default, the labels are visited in random order for both the swap and alpha-expansion moves */
    /* Use this function with boolean argument 0 to fix the order to be not random                    */
    /* Use this function with argumnet 1 to fix the order back to random                              */
    void setLabelOrder(bool RANDOM_LABEL_ORDER);
	void setParameters(int numParam, void *param);
    /* Returns Data Energy of current labeling */
    EnergyType dataEnergy();
    /* Returns Smooth Energy of current labeling */
    EnergyType smoothnessEnergy();

protected:
    /* This function is used to set the data term, and it can be used only if dataSetup = PASS_AS_PARAMETER */
    /* DataCost is an array s.t. the data cost for pixel p and  label l is stored at                        */
    /* DataCost[pixel*num_labels+l].  If the current neighborhood system is a grid, then                    */
    /* the data term for label l and pixel with coordinates (x,y) is stored at                              */ 
    /* DataCost[(x+y*width)*num_labels + l]. Thus the size of array DataCost is num_pixels*num_labels       */
     void setData(EnergyTermType *DataCost);
    /* This function is used to set the data term, and it can be used only if dataSetup = PASS_AS_PARAMETER */
    /* dataFn is a pointer to a function  f(Pixel p, Label l), s.t. the data cost for pixel p to have       */
    /* label l  is given by f(p,l) */
     void setData(DataCostFn dataFn);
	 /* This function is used to set the smoothness term, and it can be used only if                         */
     /* smoothSetup = PASS_AS_PARAMETER                                                                      */
     /*  V is an array of costs, such that V(label1,label2)  is stored at V[label1+num_labels*label2]        */
     /* If graph is a grid, then using this  function only if the smooth costs are not spacially varying     */
     /* that is the smoothness penalty V depends only on labels, but not on pixels                           */
     void setSmoothness(EnergyTermType* V);
	 void setSmoothness(int smoothExp,CostVal smoothMax, CostVal lambda);
     /* This function is used to set the smoothness term, and it can be used only if                         */
     /* smoothSetup = PASS_AS_PARAMETER AND the graph is a grid                                              */
     /* hCue(x,y,label1,label2) is the smoothness penalty for pixels (x,y) and (x+1,y) to have labels        */
     /* label1 and label2, correspondingly.                                                                  */
     /* vCue(x,y,label1,label2) is the smoothness penalty for pixels (x,y) and (x,y+1) to have labels        */
     /* label1 and label2, correspondingly.                                                                  */
     void setCues(CostVal* hCue, CostVal* vCue); // hCue and vCue must be arrays of size width*height in row major order. 
     /* This function is used to set the smoothness term, and it can be used only if                         */
     /* smoothSetup = PASS_AS_PARAMETER AND the graph is a grid                                              */
     /* horz_cost is a function f(x,y,label1,label2) such that smoothness penalty for neigboring pixels      */
     /* (x,y) and (x+1,y) to have labels, respectively, label1 and label2 is f(x,y,label1,label2)            */
     /* vert_cost is a function f(x,y,label1,label2) such that smoothness penalty for neigboring pixels      */
     /* (x,y) and (x,y+1) to have labels, respectively, label1 and label2 is f(x,y,label1,label2)            */
     void setSmoothness(SmoothCostGeneralFn cost);
	 virtual void optimizeAlg(int nIterations) = 0;
	 typedef struct NeighborStruct{
        PixelType  to_node;
        EnergyTermType weight;
	 } Neighbor;
    LabelType *m_labeling;
    bool m_random_label_order;
    bool m_needToFreeV;
    EnergyTermType *m_datacost;
    EnergyTermType *m_smoothcost;
    EnergyTermType *m_vertWeights;
    EnergyTermType *m_horizWeights;
    LinkedBlockList *m_neighbors;
    LabelType *m_labelTable;
    PixelType *m_lookupPixVar;
    EnergyTermType m_weight;
    /* Pointers to function for energy terms */
    DataCostFn m_dataFnPix;
    SmoothCostGeneralFn m_smoothFnPix;
    void commonGridInitialization( PixelType width, PixelType height, int nLabels);
    void commonNonGridInitialization(PixelType num_pixels, int num_labels);
    void commonInitialization();    
    void scramble_label_table();
    EnergyType giveDataEnergyArray();
    EnergyType giveDataEnergyFnPix();
    EnergyType giveSmoothEnergy_G_ARRAY_VW();
    EnergyType giveSmoothEnergy_G_ARRAY();
    EnergyType giveSmoothEnergy_G_ARRAY_SECOND_ORDER();
    EnergyType giveSmoothEnergy_G_FnPix();
    EnergyType giveSmoothEnergy_NG_ARRAY();
    EnergyType giveSmoothEnergy_NG_FnPix();
    void add_t_links_ARRAY(Energy *e,Energy::Var *variables,int size,LabelType alpha_label);
    void add_t_links_FnPix(Energy *e,Energy::Var *variables,int size,LabelType alpha_label);
    void initialize_memory();
    void terminateOnError(bool error_condition,const char *message);
};


class Swap: public GCoptimization
{
public:
    Swap(PixelType width,PixelType height,int num_labels, EnergyFunction *eng);
    Swap(PixelType nPixels, int num_labels, EnergyFunction *eng);
    ~Swap();
    /* Peforms swap algorithm. Runs it the specified number of iterations */
    EnergyType swap(int max_num_iterations);
    /* Peforms swap algorithm. Runs it until convergence */
    EnergyType swap();
    /* Peforms one iteration (one pass over all labels)  of swap algorithm.*/
    EnergyType oneSwapIteration();
    /* Peforms  swap on a pair of labels, specified by the input parameters alpha_label, beta_label */
    EnergyType alpha_beta_swap(LabelType alpha_label, LabelType beta_label);

protected:
    void optimizeAlg(int nIterations);

private:
    PixelType *m_pixels;
    void set_up_swap_energy_G_ARRAY_VW(int size, LabelType alpha_label,LabelType beta_label,PixelType *pixels,Energy* e, Energy::Var *variables);
    void set_up_swap_energy_G_ARRAY(int size, LabelType alpha_label,LabelType beta_label,PixelType *pixels,Energy* e, Energy::Var *variables);
	void set_up_swap_energy_G_ARRAY_SECOND_ORDER(int size, LabelType alpha_label,LabelType beta_label,PixelType *pixels,Energy* e, Energy::Var *variables);
    void set_up_swap_energy_G_FnPix(int size, LabelType alpha_label,LabelType beta_label,PixelType *pixels,Energy* e, Energy::Var *variables);
    void set_up_swap_energy_NG_ARRAY(int size, LabelType alpha_label,LabelType beta_label,PixelType *pixels,Energy* e, Energy::Var *variables);     
    void set_up_swap_energy_NG_FnPix(int size, LabelType alpha_label,LabelType beta_label,PixelType *pixels,Energy* e, Energy::Var *variables);     
    EnergyType start_swap(int max_iterations);
    void add_t_links_ARRAY_swap(Energy *e,Energy::Var *variables,int size,LabelType alpha_label,LabelType beta_label,PixelType *pixels);
    void add_t_links_FnPix_swap(Energy *e,Energy::Var *variables,int size,LabelType alpha_label,LabelType beta_label,PixelType *pixels);
    void perform_alpha_beta_swap(LabelType alpha_label, LabelType beta_label);
};


class Expansion: public GCoptimization
{
public:
    Expansion(PixelType width,PixelType height,int num_labels,EnergyFunction *eng):GCoptimization(width,height,num_labels,eng){};
    Expansion(PixelType nPixels, int num_labels,EnergyFunction *eng):GCoptimization(nPixels,num_labels,eng){};
    /* Peforms expansion algorithm. Runs the number of iterations specified by max_num_iterations */
    /* Returns the total energy of labeling   */
    EnergyType expansion(int max_num_iterations);
    /* Peforms expansion algorithm. Runs it until convergence */
    /* Returns the total energy of labeling   */
    EnergyType expansion();
    /* Peforms one iteration (one pass over all labels)  of expansion algorithm.*/
    EnergyType oneExpansionIteration();
    /* Peforms  expansion on one label, specified by the input parameter alpha_label */
    EnergyType alpha_expansion(LabelType alpha_label);
    
protected:
    void optimizeAlg(int nIterations);

private:
    void set_up_expansion_energy_G_ARRAY_VW(int size, LabelType alpha_label,Energy* e, Energy::Var *variables);
    void set_up_expansion_energy_G_ARRAY(int size, LabelType alpha_label,Energy* e, Energy::Var *variables);
	void set_up_expansion_energy_G_ARRAY_SECOND_ORDER(int size, LabelType alpha_label,Energy* e, Energy::Var *variables);
    void set_up_expansion_energy_G_FnPix(int size, LabelType alpha_label,Energy* e, Energy::Var *variables);
    void set_up_expansion_energy_NG_ARRAY(int size, LabelType alpha_label,Energy* e, Energy::Var *variables);       
    void set_up_expansion_energy_NG_FnPix(int size, LabelType alpha_label,Energy* e, Energy::Var *variables);       
    void perform_alpha_expansion(LabelType label);  
    EnergyType start_expansion(int max_iterations);
};

#endif


/***************************************************************************************************************************************/
#ifndef _reg_h
#define _reg_h

#define UP 0
#define DOWN 1
#define LEFT 2
#define RIGHT 3

class MaxProdBP;
class OneNodeCluster
{
public:
  OneNodeCluster();
  //static const int UP = 0, DOWN = 1, LEFT = 2, RIGHT = 3;
  static int numStates;
  float   *receivedMsgs[4],
              *nextRoundReceivedMsgs[4],
              *localEv;
//   float *psi[4];
  void ComputeMsgRight(float *msgDest, int r, int c, MaxProdBP *mrf);
  void ComputeMsgUp(float *msgDest, int r, int c, MaxProdBP *mrf);
  void ComputeMsgLeft(float *msgDest, int r, int c, MaxProdBP *mrf);
  void ComputeMsgDown(float *msgDest, int r, int c, MaxProdBP *mrf);
  void getBelief(float *beliefVec);
  int msgChange(float thresh);
  void deliverMsgs();
};

void initOneNodeMsgMem(OneNodeCluster *nodeArray, float *memChunk, const int numNodes, const int msgChunkSize);
void computeMessagesUpDown(OneNodeCluster *nodeArray, const int numCols, const int numRows, const int currCol, const float alpha, MaxProdBP *mrf);
void computeMessagesLeftRight(OneNodeCluster *nodeArray, const int numCols, const int numRows, const int currRow, const float alpha, MaxProdBP *mrf);
void computeOneNodeMessagesPeriodic(OneNodeCluster *nodeTopArray, OneNodeCluster *nodeBotArray, const int numCols, const float alpha);

#endif


/***************************************************************************************************************************************/
#ifndef __MAXPRODBP_H__
#define __MAXPRODBP_H__

class MaxProdBP : public MRFmodel{
public:
    MaxProdBP(int width, int height, int nLabels, EnergyFunction *eng);
    MaxProdBP(int nPixels, int nLabels,EnergyFunction *eng);
    ~MaxProdBP();
    void setNeighbors(int pix1, int pix2, CostVal weight);
    //void clearNeighbors();
    Label getLabel(int pixel){return(m_answer[pixel]);};
    void setLabel(int pixel,Label label){m_answer[pixel] = label;};
    Label* getAnswerPtr(){return(m_answer);};
    void clearAnswer();
    void setParameters(int /*numParam*/, void * /*param*/){printf("No optional parameters to set");}
    EnergyVal smoothnessEnergy();
    EnergyVal dataEnergy();
    int getWidth();
    int getHeight();
    float *getScratchMatrix();
    int getNLabels();
    bool varWeights();
    float getExpV(int i);
    float *getExpV();
    CostVal getHorizWeight(int r, int c);
    CostVal getVertWeight(int r, int c);
protected:
    void setData(DataCostFn dcost); 
    void setData(CostVal* data);    
    void setSmoothness(SmoothCostGeneralFn cost);
    void setSmoothness(CostVal* V);
    void setSmoothness(int smoothExp,CostVal smoothMax, CostVal lambda);
    void setCues(CostVal* hCue, CostVal* vCue); 
    void initializeAlg();
    void optimizeAlg(int nIterations);
private:
    Label *m_answer;
    CostVal *m_V;
    CostVal *m_D;
    CostVal *m_horizWeights;
    CostVal *m_vertWeights;
    float m_exp_scale;
    DataCostFn m_dataFn;
    SmoothCostGeneralFn m_smoothFn;
    bool m_needToFreeV;
    float *m_scratchMatrix;
    float *m_ExpData;
    float *m_message_chunk;
    OneNodeCluster *nodeArray;
    typedef struct NeighborStruct
	{
		int     to_node;
		CostVal weight;
    } Neighbor;
    LinkedBlockList *m_neighbors;
};

#endif


/***************************************************************************************************************************************/
#ifndef __ICM_H__
#define __ICM_H__

class ICM : public MRFmodel{
public:
    ICM(int width, int height, int nLabels, EnergyFunction *eng);
    ICM(int nPixels, int nLabels,EnergyFunction *eng);
    ~ICM();
    void setNeighbors(int pix1, int pix2, CostVal weight);
    //void clearNeighbors();
    Label getLabel(int pixel){return(m_answer[pixel]);};
    void setLabel(int pixel,Label label){m_answer[pixel] = label;};
    Label* getAnswerPtr(){return(m_answer);};
    void clearAnswer();
    void setParameters(int /*numParam*/, void * /*param*/){printf("No optional parameters to set");}
    EnergyVal smoothnessEnergy();
    EnergyVal dataEnergy();
protected:
    void setData(DataCostFn dcost); 
    void setData(CostVal* data);    
    void setSmoothness(SmoothCostGeneralFn cost);
    void setSmoothness(CostVal* V);
    void setSmoothness(int smoothExp,CostVal smoothMax, CostVal lambda);
    void setCues(CostVal* hCue, CostVal* vCue); 
    void initializeAlg();
    void optimizeAlg(int nIterations);
private:
    Label *m_answer;
    CostVal *m_V;
    CostVal *m_D;
    CostVal *m_horizWeights;
    CostVal *m_vertWeights;
    DataCostFn m_dataFn;
    SmoothCostGeneralFn m_smoothFn;
    bool m_needToFreeV;
    typedef struct NeighborStruct{
        int     to_node;
        CostVal weight;
    } Neighbor;
    LinkedBlockList *m_neighbors;
};


#endif


/***************************************************************************************************************************************/
enum GraphCutOption { GCEXP,GCSWAP,MPBP};
const GraphCutOption GraphCutChoice[3]={GCSWAP,GCEXP,MPBP};


/*static MRFmodel* optimizePottsFirst(unsigned short w, unsigned short h, unsigned short C, float *UX, float beta, int iter, GraphCutOption option)
{
	DataCost *data = new DataCost(UX);
	SmoothnessCost *smooth = new SmoothnessCost(1,1,beta);
	EnergyFunction *eng = new EnergyFunction(data,smooth);
	MRFmodel *mrf;
	switch(option)
	{
	case GCEXP:
		mrf = new Expansion(w,h,C,eng);
		break;
	case GCSWAP:
		mrf = new Swap(w,h,C,eng);
		break;
	case MPBP:
		mrf = new MaxProdBP(w,h,C,eng);
		break;
	}
	mrf->initialize();
	mrf->clearAnswer();
	mrf->optimize(iter);
	delete[] UX,eng,smooth,data;
	return mrf;
}


static MRFmodel* optimizePottsSecond(unsigned short w, unsigned short h, unsigned short C, float *UX, float beta, int iter, GraphCutOption option)
{
	DataCost *data = new DataCost(UX);
	SmoothnessCost *smooth = new SmoothnessCost(1,1,beta);
	EnergyFunction *eng = new EnergyFunction(data,smooth);
	MRFmodel *mrf;
	switch(option)
	{
	case GCEXP:
		mrf = new Expansion(w*h,C,eng);
		break;
	case GCSWAP:
		mrf = new Swap(w*h,C,eng);
		break;
	case MPBP:
		mrf = new MaxProdBP(w*h,C,eng);
		break;
	}
	unsigned short m,n;
	short u,v;
	for(m=1;m<w-1;m++)
		for(n=1;n<h-1;n++)
			for(u=-1;u<2;u++)
				for(v=-1;v<2;v++)
					if(u*u+v*v>0)
						mrf->setNeighbors(m+n*w,m+u+(n+v)*w,1);
	mrf->initialize();
	mrf->clearAnswer();
	mrf->optimize(iter);
	//delete eng,smooth,data;
	return mrf;
}*/
