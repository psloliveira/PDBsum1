// ==================================================================
// KdTree.h
// Copyright (c) Jonathan Barker, 2002
// ==================================================================

// ==================================================================
// Declaration of abstract type Region.
// ==================================================================

#ifndef REGION_H
#define REGION_H

// ==================================================================
// Forward declarations
// ==================================================================
// Region				An abstraction of a region
// ==================================================================

typedef struct _Region Region;

// ==================================================================
// type Region
// ==================================================================
// intersectionQ		The intersection oracle
// inclusionQ			The inclusion oracle
// free(R)				Destroy region (free memory)
// ==================================================================

struct _Region
{
	int (*intersectionQ)(Region*,double*,double*,int);
	int (*inclusionQ)(Region*,double*,int);
	void (*free)(Region*);
};

// ==================================================================

#endif

// ==================================================================
// Declaration of type KdTree, related types and methods.
// ==================================================================

#ifndef KDTREE_H
#define KDTREE_H

// ==================================================================
// Forward declarations
// ==================================================================
// KdTree					A kd-tree
// KdTreeQuery				Result of query on a KdTree
// ==================================================================

typedef struct _KdTree KdTree;
typedef struct _KdTreeQuery KdTreeQuery;

// ==================================================================
// Methods of type KdTree
// ==================================================================
// create(u,n,k)			Create kd-tree on u[0],...,u[n-1]
// free(K)					Free the kd-tree K
// query(K,R)				Initialise a query object (see code)
// ==================================================================

extern KdTree *KdTree_create(double**,int,int);
extern void KdTree_free(KdTree*);
extern KdTreeQuery *KdTree_query(KdTree*,Region*);

// ==================================================================
// Methods of type KdTreeQuery
// ==================================================================
// free(Q)					Free memory for query object and Region
// next(Q)					Return next result in Q (-1 if no more)
// ==================================================================

extern void KdTreeQuery_free(KdTreeQuery*);
extern int KdTreeQuery_next(KdTreeQuery*);

// ==================================================================

#endif

// ==================================================================
// Declaration of Annulus construction and destruction.
// ==================================================================

#ifndef ANNULUS_H
#define ANNULUS_H

// ==================================================================
// Methods for Annulus manipulation
// ==================================================================
// create(u,a,b,d)		Make region {x in R^d : a <= |x-u| <= b }.
// free(A)				Free the region given (or use R->free)
// ==================================================================

extern Region *Annulus_create(double*,double,double,int);
extern void Annulus_free(Region*);

// ==================================================================

#endif

