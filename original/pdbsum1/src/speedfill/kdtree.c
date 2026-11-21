// ==================================================================
// KdTree.c
// Copyright (c) Jonathan Barker, 2002
// ==================================================================
// Implementation of type KdTree, related types and methods.
// ==================================================================

#include "kdtree.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

// ==================================================================
// Forward declarations of local types
// ==================================================================
// KdTreeNode                        One node of a KdTree
// ==================================================================

typedef struct _KdTreeNode KdTreeNode;

// ==================================================================
// Local type KdTreeNode
// ==================================================================
// type                                        One of 0,1,2 for a node or -1 for a leaf
// index                                Index of splitting coordinate or point
// left,right                        Branches of tree (unless this is a leaf)
// min,max                                The box region encompassed by this node
// depth                                The depth of the tree below this node
// ==================================================================

struct _KdTreeNode
{
  int type;
  int index;
  KdTreeNode *left;
  KdTreeNode *right;
  double *min;
  double *max;
  int depth;
};

// ==================================================================
// Declaration of methods of local type KdTreeNode
// ==================================================================
// create(...)                        Creates a new node and its descendants
// free(N)                                Frees a node and all its descendants
// ==================================================================

static KdTreeNode *KdTreeNode_create(int*,int,int,double**,int);
static void KdTreeNode_free(KdTreeNode*);

// ==================================================================
// type KdTree
// ==================================================================
// root                                        The root of the tree
// dim                                        Dimension of the points
// ==================================================================

struct _KdTree
{
  KdTreeNode *root;
  int dim;
};

// ==================================================================
// type KdTreeQuery
// ==================================================================
// tree                                        The tree which the query relates to
// region                                The region being queried (see Region.h)
// count                                Number of nodes on the stack
// stack                                The stack of nodes in the query
// ==================================================================

struct _KdTreeQuery
{
  KdTree *tree;
  Region *region;
  int count;
  KdTreeNode *stack[0];
};


// ==================================================================
// Declaration of private methods/members of type KdTree
// ==================================================================
// compare(pa,pb)                Used during node creation (qsort)
// data,index                        Static globals (see KdTreeNode_create)
// ==================================================================

static int KdTree_compare(const void*,const void*);
static double **KdTree_data;
static int KdTree_index;

// ==================================================================
// Local functions
// ==================================================================

static double dmin(double x, double y)
{
  return x<y ? x:y;
}

static double dmax(double x, double y)
{
  return x>y ? x:y;
}

static int imin(int x, int y)
{
  return x<y ? x:y;
}

static int imax(int x, int y)
{
  return x>y ? x:y;
}

// ==================================================================
// Public methods of type KdTree
// ==================================================================

KdTree *KdTree_create(double **u, int n, int d)
{
  KdTree *K;
  int *tmp;
    int i;

    if(n<1 || d<1 || !u) return NULL;

    // 1. Create memory for the object.

    K = (KdTree*)calloc(1,sizeof(KdTree));
    K->dim=d;

    // 3a. Create a temporary array to hold indices

    tmp = (int*)calloc(n,sizeof(int));
    for(i=0; i<n; i++) tmp[i]=i;

    // 3b. Create the tree recursively. This takes time
    // of order at most n.log(n)^2, assuming that qsort
    // always manages n.log(n) and d is constant.

    K->root = KdTreeNode_create(tmp,n,0,u,d);
    free(tmp);

    // 4. Return the result!

    return K;
}

void KdTree_free(KdTree *K)
{
  if(K)
    {
      KdTreeNode_free(K->root);
      free(K);
    }
}

KdTreeQuery *KdTree_query(KdTree *K, Region *R)
{
  KdTreeQuery *Q;
  int rq;

  rq = sizeof(KdTreeQuery)+K->root->depth*sizeof(KdTreeNode*);

  Q = (KdTreeQuery*)calloc(1,rq);
  Q->tree=K;
  Q->region=R;
  Q->count=1;
  Q->stack[0]=K->root;

  return Q;
}

// ==================================================================
// Methods of type KdTreeQuery
// ==================================================================

int KdTreeQuery_next(KdTreeQuery *Q)
{
  KdTreeNode *N;
  Region *R = Q->region;
  KdTreeNode **stack=&(Q->stack[0]);
  int dim = Q->tree->dim;
  int *count = &(Q->count);

  // Until the stack is empty (or we return inside
  // the while loop...

  while(*count>0)
    {
      // Pull the top node off the stack.

      N = stack[--(*count)];

      // If the node is a leaf we simply test it and
      // remove it from the stack. If the point is in
      // the query region, return it; otherwise continue
      // with the rest of the stack.

      if(N->type<0)
        {
          if(R->inclusionQ(R,N->min,dim))
            {
              return N->index;
            }
          else
            {
              continue;
            }
        }

      // So the node is internal (ie not a leaf).
      // If the query region does not intersect
      // the node's region then we can remove it
      // and continue with the rest of the stack.

      if(!R->intersectionQ(R,N->min,N->max,dim))
        {
          continue;
        }

      // The query region *does* intersect the node's
      // region. So now we must place the child nodes
      // onto the stack.

      stack[(*count)++]=N->left;
      stack[(*count)++]=N->right;
    }

  return -1;
}

void KdTreeQuery_free(KdTreeQuery *Q)
{
  if(Q) free(Q);
}

// ==================================================================
// Private methods of type KdTree
// ==================================================================

static int KdTree_compare(const void *pa, const void *pb)
{
  const int a = *((const int*)pa);
  const int b = *((const int*)pb);
  double c = KdTree_data[a][KdTree_index];
  double d = KdTree_data[b][KdTree_index];

  return c<d ? -1 : c>d ? 1 : 0;
}

// ==================================================================
// Methods of local type KdTreeNode
// ==================================================================

static KdTreeNode *KdTreeNode_create(int *idx,int n,int type,double **u,int dim)
{
  KdTreeNode *N;
  int split;
  int i,rq;

  // 1. The really easy case. If n is 0 do nothing!

  if(n<=0) return NULL;

  // 1.5. We'll need to create a node in all other cases.

  rq = sizeof(KdTreeNode)+dim*2*sizeof(double);
  N = (KdTreeNode*)calloc(1,rq);
  N->min=(double*)&N[1];
  N->max=&N->min[dim];

  // 2. The easy case. If n is 1, create a leaf.

  if(n==1)
    {
      N->type=-1;
      N->index=idx[0];
      N->depth=1;
      memcpy(N->min,u[idx[0]],sizeof(double)*dim);
      memcpy(N->max,u[idx[0]],sizeof(double)*dim);

      return N;
    }

  // 2.5. Now we need to order the indices by coordinate
  // numbered type. THIS IS NOT THREAD-SAFE. This kludge
  // is used because it's not possible to pass extra
  // parameters to qsort.

  KdTree_data=u;
  KdTree_index=type;
  qsort(idx,n,sizeof(int),KdTree_compare);

  // 3. The recursive case. Find [n/2] and split the array into
  // two pieces. Create a node whose splitting value is the median.
  // But make sure that if there are several entries with the same
  // coordinate that we take the right-most.

  split = n/2;
  N->index=idx[split-1];
  while(split<n-1 && u[split+1][type]==u[split][type]) split++;
  N->type=type;

  // Now create the left and right branches of the node.

  type = (type+1)%dim;
  N->left = KdTreeNode_create(idx,split,type,u,dim);
  N->right = KdTreeNode_create(&idx[split],n-split,type,u,dim);

  // Compute max,min and depth...

  N->depth = imax(N->left->depth,N->right->depth)+1;

  for(i=0; i<dim; i++)
    {
      N->min[i]=dmin(N->left->min[i],N->right->min[i]);
      N->max[i]=dmax(N->left->max[i],N->right->max[i]);
    }

  // We're done...

  return N;
}

static void KdTreeNode_free(KdTreeNode *N)
{
  if(N)
    {
      KdTreeNode_free(N->left);
      KdTreeNode_free(N->right);
      free(N);
    }
}

// ==================================================================

// ==================================================================
// Annulus.c
// Copyright (c) Jonathan Barker, 2002
// ==================================================================
// Implementation of the Annulus region and its oracles
// ==================================================================

// ==================================================================
// Local "functions"
// ==================================================================

#define min(x,y) (x<y ? x:y)
#define max(x,y) (x>y ? x:y) 

// ==================================================================
// type Annulus
// ==================================================================
// centre                                The center of the annulus
// min,max                                Limits of the radii
// dim                                        The dimension of the space
// ==================================================================

typedef struct _Annulus
{
  double *centre;
  double min;
  double max;
  int dim;
}
Annulus;

// ==================================================================
// Oracles for type Annulus
// ==================================================================

static int Annulus_po(Region *vA, double *x, int d)
{
  Annulus *A=(Annulus*)&vA[1];
  double tmp,sum;
  int i;

  // Does x lie within annulus A?

  if(A->dim!=d) return 0;

  for(sum=0.0,i=0; i<d; i++)
    {
      tmp = A->centre[i]-x[i];
      sum += tmp*tmp;
    }

  return sum<A->min || sum>A->max ? 0:1;
}

static int Annulus_ro(Region *vA, double *minBox, double *maxBox, int d)
{
  Annulus *A=(Annulus*)&vA[1];
  double minSum;
  double maxSum;
  double t1,t2;
  int i;

  if(d!=A->dim) return 0;

  // Does the box region [minBox,maxBox] intersect the annulus A?

  minSum=0.0;
  maxSum=0.0;
  for(i=0; i<d; i++)
    {
      t1 = A->centre[i]-minBox[i];
      t2 = A->centre[i]-maxBox[i];
      t1 *= t1;
      t2 *= t2;

      if(minBox[i]>A->centre[i] || maxBox[i]<A->centre[i])
        {
          minSum += min(t1,t2);
        }

      maxSum += max(t1,t2);
    }

  return minSum>A->max || maxSum<A->min ? 0:1;
}


// ==================================================================
// Methods of for regions of type Annulus
// ==================================================================

Region *Annulus_create(double *u, double a, double b, int d)
{
  Annulus *A;
  Region *R;
  int rq;

  rq = sizeof(Region)+sizeof(Annulus)+d*sizeof(double);
  R = (Region*)calloc(1,rq);
  A = (Annulus*)&R[1];
  R->intersectionQ=Annulus_ro;
  R->inclusionQ=Annulus_po;
  R->free=Annulus_free;

  A->centre=(double*)&A[1];
  memcpy(A->centre,u,sizeof(double)*d);
  A->min=a*a;
  A->max=b*b;
  A->dim=d;

  return R;
}

void Annulus_free(Region *R)
{
  if(R) free(R);
}

// ==================================================================

