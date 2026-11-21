/***********************************************************************

  ras3d.h - Include file for common Raster3D routines in ras3d.c

***********************************************************************/

/*------------------ Standard includes ----------------------*/

#ifndef STANDARD_INCLUDES
#define STANDARD_INCLUDES

/* Standard system include files */
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <ctype.h>

#endif

/*-----------------------------------------------------------*/

/*------------------ Common ---------------------------------*/

#ifndef COMMON
#define COMMON

/* Flag-values */
#define FALSE                0
#define TRUE                 1
#define OK                   0

//* Array parameters */
#define COLNAM_LEN          12
#define COLNO_LEN           15
#define MAX_OBJ_COL          8
#define NAMINO              20
#define NATOM_TYPES          8
#define NCONSERV            10
#define NRESTYPES            7

/* Residue conservation ranges */
static double RANGE[NCONSERV + 1]
    = { 0.00, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70,  0.80, 0.90,
	1.000001 };

#endif

/*-----------------------------------------------------------*/

/*----------------- Structures ------------------------------*/

#ifndef STRUCTURES
#define STRUCTURES

/* Colouring schemes for ligands and buried vertices */
#define NONE                 -1
#define STANDARD              0
#define CPK                   1
#define RESIDUE_TYPE          2
#define CONSERVATION          3
#define BURIED_EXPOSED        4
#define BY_CHAIN              5
#define ELECTROSTATIC         6
#define BY_BVALUE             7

#define NCOLOUR_SCHEMES       8

/* Electrostatic potential colouring parameters */
#define NEG_UPPER          -60.0
#define NEG_LOWER           -5.0
#define POS_UPPER           60.0
#define POS_LOWER            5.0

/* Structure for coordinate limits
   ------------------------------- */
struct limits
{
  int                     first;
  double                  valmin[3];
  double                  valmax[3];
  double                  c_of_g[3];
  double                  max_radius;
  double                  min_radius;
  int                     npoints;
};

/* Structure for grid limits
   ------------------------- */
struct grid_limits
{
  int                     minpos[3];
  int                     maxpos[3];
};

/* Structure for each atom record to be written out
   ------------------------------------------------ */
struct atom
{
  int                     atom_number;
  int                     type;
  char                    atom_name[5];
  char                    element[3];
  double                  radius;
  double                  coords[3];
  double                  bvalue;
  int                     checked;
  int                     icluster;
  double                  depth;
  int                     ndepth;
  struct residue          *residue_ptr;
  struct atom             *next_atom_ptr;
};

/* Structure for each residue
   -------------------------- */
struct residue
{
  char                    res_name[4];
  char                    res_num[6];
  char                    chain;
  int                     type;
  int                     natoms;
  int                     in_range;
  double                  conservation_score;
  int                     checked;
  struct atom             *first_atom_ptr;
  struct rasdef           *rasdef_ptr;
  struct residue          *next_residue_ptr;
};

/* Structure for each vertex record
   -------------------------------- */
struct vertex
{
  int                     vertex_no;
  int                     edge_no;
  double                  coords[3];
  double                  ratio;
  double                  bvalue;
  double                  normal[3];
  int                     cluster;
  int                     deleted;
  int                     exposed;
  int                     npolygons;
  int                     sort_pos[3];
  struct atom             *atom_ptr;
  struct polist           *first_polist_ptr;
  struct polist           *last_polist_ptr;
  struct vertex           *next_vertex_list_ptr;
  struct vertex           *next_sorted_vertex_ptr[3];
  struct vertex           *next_vertex_ptr;
};

#endif

/*-----------------------------------------------------------*/

/*----------------- Raster3D data items and functions -------*/

#ifndef RAS3D_H
#define RAS3D_H

/* Data types */

/* Structure for each polist record giving each vertex's list of polygons
   ---------------------------------------------------------------------- */
struct polist
{
  struct polygon          *polygon_ptr;
  struct polist           *next_polist_ptr;
};

/* Structure for each polygon record
   --------------------------------- */
struct polygon
{
  int                     polygon_no;
  int                     polygon_colour;
  int                     outside_point;
  struct vertex           *vertex_ptr[3];
  struct polygon          *next_polygon_ptr;
};

/* Routines */
extern double calculate_rascale(struct limits vertex_limits,double border);
extern struct vertex *binary_edge_search(int edge_no,
					 struct vertex **vertex_index_ptr,
					 int nvertices);
extern struct polist *create_polist_record(struct polist **fst_polist_ptr,
					   struct polist **lst_polist_ptr,
					   struct polygon *polygon_ptr);
extern struct polygon *create_polygon_record(struct polygon **fst_polygon_ptr,
					     struct polygon **lst_polygon_ptr,
					     int polygon_no,
					     struct vertex *vertex1_ptr,
					     struct vertex *vertex2_ptr,
					     struct vertex *vertex3_ptr,
					     int low_point);
extern int make_line(int ncross_face,int nlines,int cell_line[12][2],
		     struct vertex *line_vertex_ptr[12][2],
		     struct vertex *side_vertex_ptr[4],int on_edge[4],
		     int last_value,int contour_level);
extern int generate_triangles(struct polygon **fst_polygon_ptr,
			      struct polygon **lst_polygon_ptr,int npol,
			      struct vertex *polygon_vertex_ptr[12],int nvert,
			      int store_low_point[12],double corner_coords[3],
			      int corner_value[8],double grid_sep,
			      int contour_level,int npoints[3],
			      double valmin[3],double valmax[3]);
extern int make_polygons(int nlines,struct vertex *line_vertex_ptr[12][2],
			 int cell_line[12][2],int contour_level,
			 int low_point[12],double corner_coords[3],
			 int corner_value[8],double grid_sep,
			 int npoints[3],double valmin[3],double valmax[3],
			 struct polygon **fst_polygon_ptr,
			 struct polygon **lst_polygon_ptr,int *npol);
extern int generate_polygons(int *grid,int npoints[3],double valmin[3],
			     double valmax[3],double grid_sep,
			     int contour_level,
			     struct grid_limits glimits,
			     struct vertex **vertex_index_ptr,int nvertices,
			     struct polygon **fst_polygon_ptr,
			     struct atom *dummy_atom_ptr,
			     int buried_vertices_only);
extern void vecprd(double px,double py,double pz,double qx,double qy,double qz,
		   double *vx, double *vy, double *vz);
extern double dot3(double x1,double y1,double z1,double x2,double y2,
		   double z2);
extern int orien3(double x1,double y1,double z1,double x2,double y2,double z2,
		  double x3,double y3,double z3,double ex,double ey,double ez);
extern void calculate_normal(double x1,double y1,double z1,
			     double x2,double y2,double z2,
			     double x3,double y3,double z3,double normal[3]);
extern void order_polygon_vertices(struct polygon *first_polygon_ptr,
				   int npoints[3],double valmin[3],
				   double valmax[3],double grid_sep);
extern void calculate_vertex_normals(struct vertex *first_vertex_ptr,
				     struct limits *vertex_limits);
extern void get_colour_rgb(char *c_name,double rgb[3]);
extern void write_raster3d_headers(FILE *outfile_ptr,char *prog_name,
				   char *background_colour,int sizex,
				   int sizey);
extern int get_colour_cpk(char *element,double rgb[3]);
extern int get_colour_restype(char *res_name,double rgb[3]);
extern int get_colour_conserv(double conservation,double rgb[3],
			      char *conservation_colours);
extern int get_chain_colour(char chain,char *colour_name,
			    char *colour_numbers);
extern void write_raster3d_polygons(FILE *outfile_ptr,
                             struct polygon *first_polygon_ptr,
                             double c_of_g[3],double rascale,char chain,
                             int vertex_colouring,
			     char *conservation_colours,
			     char *col_name);
extern int get_grid_value(int grid_value,int *masked,int contour_level);
extern int extract_grid_value(int grid_value,int *masked);

#endif

/*-----------------------------------------------------------*/

