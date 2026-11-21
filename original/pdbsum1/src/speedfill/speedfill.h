/***********************************************************************

  speedfill.h - Include file for speedfill.c

***********************************************************************/

/* Include files for Jonathan Barker's KD-tree routines */
#include "kdtree.h"

/* Include files for common Raster3D routines */
#include "ras3d.h"

/*------------------ Common ---------------------------------*/

#ifndef COMMON
#define COMMON

/* Array parameters */
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

/* Maths parameters */
#define SQRT_THREE   (1.732050808)

/* Parameters for exponentials lookup table */
#define RADSTEP           (0.1)
#define DISTEP            (0.1)

/* Array parameters */
#define FILENAME_LEN       512
#define LINELEN            132
#define LISTLEN            512
#define MAX_ARRAY_SIZE 10000000
#define MAX_CHAINS         100
#define MAX_LOCAL_PDB       27
#define MAX_LOOPS           10
#define MAX_NEIGHBOURS      40
#define MAX_PATTN_LEN      300
#define MAX_SPHERE_POINTS  800
#define MAX_STORE          500
#define MAX_TOKEN           20
#define MAX_TOP_GAPS        50
#define MIN_DELETE           5
#define NAME_LEN           512
#define NAME_LINES           6
#define NHASH                2
#define NTOP_BS             10
#define NAMLEN             600
#define PCDIFF            5.0
#define STRING_LEN        1000
#define TOKEN_LEN           50

/* Definition of operating system */
#define UNIX               0
#define WINDOWS            1
#define OP_SYS          UNIX

/* Commands set according to operating system */
char CP_COMMAND[6];
char MV_COMMAND[6];
char RM_COMMAND[6];
char SLEEP_COMMAND[30];

/* Grid-point values */
#define EMPTY                 0
#define MASKED                1
#define SPHERE_LOC            2
#define MAX_SPHERE            3
#define SMALLER_SPHERE        4

/* Hash types */
#define ATOM_ATOM             0
#define ATOM_TILE             1

/* Atom types */
#define UNKNOWN               0
#define NITROGEN              1
#define CARBON                2
#define OXYGEN                3
#define SULPHUR               4
#define PHOSPHORUS            5
#define IRON                  6

/* Record types defined in the rasmol.dfn file */
#define OTHER                -1
#define CA_ONLY               0
#define LIGAND                1
#define METAL                 2
#define NUCLEIC_ACID          3
#define PROTEIN               4
#define LIGAND_RANGE          5

/* Vertex types */
#define OK                    0
#define DELETED               1

/* Sphere clustering and vertex-filtering methods */
#define OLD                   0
#define NEW                   1

/* Cluster sort methods */
#define SORT_BY_SPHERES       0
#define SORT_BY_VOLUME        1

/* Plot parameters */
#define BOND_RADIUS          (0.25)
#define BORDER               (0.05)

/* Numerical parameters */
#define CONTOUR_LEVEL        (100)
#define COV_DIST2            (3.6)
#define GRID_SCALE_FACTOR    (10000)
#define NEIGHBOUR_DIST2      (16.0)
#define OVERLAP_DIST         (2.00)
#define OVERLAP_DIST2       (OVERLAP_DIST * OVERLAP_DIST)
#define TOO_CLOSE_DIST       (0.95)
#define TOO_CLOSE_SQRD      (TOO_CLOSE_DIST * TOO_CLOSE_DIST)
#define TOO_CLOSE_BV         (0.70)
#define TOO_CLOSE_BV_SQRD   (TOO_CLOSE_BV * TOO_CLOSE_BV)
#define VERTEX_DIST          (2.00)
#define VERTEX_DIST_SQRD    (VERTEX_DIST * VERTEX_DIST)

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
  double                  normal[3];
  double                  ratio;
  double                  bvalue;
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

/* Structure for input parameters
   ------------------------------ */
struct parameters
{
  double                  min_gap_radius;
  double                  max_gap_radius;
  double                  max_gap_span;
  int                     ntop_clusters;
  int                     indiv_atom_radii;
  int                     include_H_atoms;
  int                     include_waters;
  int                     select;
  int                     show_gap[MAX_TOP_GAPS];
  int                     no_gaps;
  int                     no_gif;
  int                     no_protein;
  int                     no_script;
  int                     named_file;
  int                     from_www;
  int                     generate_html;
  int                     generate_r3d;
  int                     output_vrml;
  int                     no_rasmol;
  int                     buried_vertices_only;
  int                     vertex_colouring;
  int                     write_surfnet_input;
  int                     read_vertices;
  int                     save_vertices;
  char                    user_id[LINELEN];
  int                     profunc;
  int                     centroid;
  char                    delphi_file[FILENAME_LEN];
  int                     store_atom_pointers;
  int                     col_by_conservation;
  int                     use_consurf_data;
  int                     have_rescons;
  int                     check_all_spheres;
  double                  background_rgb[3];
  char                    background_colour[COLNAM_LEN];
  char                    conservation_colours[NCONSERV + 1];
  double                  rotation_angle[3];
  int                     have_transformations;
  int                     r3d_sizex;
  int                     r3d_sizey;
  int                     windows;
  char                    ligand_colour[COLNAM_LEN];
  double                  metatom_radius;
  double                  ligatom_radius;
  double                  ligbond_radius;
  int                     depth_analysis;
  int                     save_grid;
  int                     filter_vertices;
  int                     raster3d_size;
  int                     delete_temp_files;
  int                     jmol_resnames;
  int                     flip;
  int                     pdb_type;
  int                     sph_resno;
  int                     pymol;
  int                     relinks;
  int                     pymol_www;
  int                     run_64;
};

/* Structure for each RasMol definition entry
   ------------------------------------------ */
struct rasdef
{
  int                type;
  int                number;
  char               res_name[4];
  char               res_num[6];
  char               chain;
  int                count;
  char               colour[COLNAM_LEN];
  char               *rasmol_ligand_def;
  char               ligand_from[7];
  char               ligand_to[7];
  struct atom        *first_atom_ptr;
  struct atom        *last_atom_ptr;
  struct rasdef      *next_rasdef_ptr;
};

/* Structure for each vjoin record
   ------------------------------- */
struct vjoin
{
  struct vertex           *vertex1_ptr;
  struct vertex           *vertex2_ptr;
  struct vjoin            *next_vjoin_ptr;
};

/* Structure for each gap sphere
   ----------------------------- */
struct sphere
{
  double                  radius;
  double                  coords[3];
  int                     cluster;
  struct sphere           *next_sphere_ptr;
};

/* Structure for each binding site stored for comparison
   ----------------------------------------------------- */
struct site
{
  char                    pdb_code[5];
  char                    chain;
  int                     with_ligand;
  double                  volume;
  double                  percen;
  double                  percen2;
  double                  depth;
  double                  eigen_value2;
  double                  eigen_value3;
  int                     count_type[NRESTYPES];
  double                  score;
};

/* Structure for each GROW data item
   ------------------------------------- */
struct grow
{
  char               data_type;
  char               res_name1[4];
  char               res_num1[6];
  char               chain1;
  char               res_name2[4];
  char               res_num2[6];
  char               chain2;
  struct grow        *next_grow_ptr;
};

/* Structure for each chain data item
   ---------------------------------- */
struct chain
{
  char               chain_id;
  char               copy_of;
  struct chain       *next_chain_ptr;
};


