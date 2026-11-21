/***********************************************************************

  ligplot.h - Include file for ligplot.c

***********************************************************************/

#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <ctype.h>

/* Flag-values */
#define FALSE           0
#define TRUE            1
#define FREE         TRUE
#define NOT_FREE    FALSE
#define OFF         FALSE
#define ON           TRUE
#define RESTORE     FALSE
#define SAVE         TRUE
#define PLUS            1
#define BOTH            0
#define MINUS          -1
#define ALL_ATOMS      -1

/* v.5.4.4 --> */
static char *T_F[3] = { "FALSE", "TRUE", "*UNSET*" };
/* <-- v.5.4.4 */

/* v.5.0--> */
/* Run-types: standard or PostScript generation for LigPlot+ */
#define STANDARD_RUN    0
#define LIGPLUS_PLOT    1
/* <--v.5.0 */

/* Object types */
#define LIGAND          1
#define HGROUP          2
#define HYDROPHOBIC     3
#define WATER           4
#define SIMPLE_HGROUP   5
/* v.5.0--> */
#define DUMMY           6
/* <--v.5.0 */

/* Residue types */
#define STANDARD        1
#define SIMPLE_LIGAND   2
/* v.5.0--> */
#define DUMMY_RESIDUE   3
/* <--v.5.0 */

/* Bond types */
#define COVALENT        0
#define HBOND           1
#define CONTACT         2
#define INTERNAL        3
/* v.5.4--> */
#define SALT_BRIDGE     4
#define DISULPHIDE      5
/* <--v.5.4 */
#define DELETED       -99

/* Bond order types */
#define SINGLE          1
#define DOUBLE          2
#define TRIPLE          3
#define AROMATIC        4

/* Bond sources */
#define CONECT          0
#define CALCULATED      1
#define HBPLUS          2
#define DRW_FILE        3

/* Hhb_info sources */
/* v.3.1--> */
/*#define EXTRA          -1 */
/* <--v.3.1 */
/* v.4.2--> */
#define RESET          -1
/* <--v.4.2 */
#define CONECT          0
#define HHB_FILE        1
#define NNB_FILE        2

/* Bond ends */
#define EITHER          0
#define FIRST           1
#define SECOND          2

/* Move types */
#define MOVE_OBJECT     0
#define SWIVEL_OBJECT   1
#define SWING_OBJECT    2
#define FLIP_OBJECT     3
#define FLIP_ABOUT_ATOM 4

/* v.5.3.4--> */
//#define NMOVE_TYPES     4
#define NMOVE_TYPES     5

//static char *MOVE_TYPE[NMOVE_TYPES] = {
//  "MOVE_OBJECT", "SWIVEL_OBJECT", "SWING_OBJECT", "FLIP_OBJECT",
//  "FLIP_ABOUT_ATOM"
//};
/* <--v.5.3.4 */

/* v.5.1--> */
/* Maximum number of rotatable bonds for all-conformation test */
/* v.5.2.1 --> */
// #define MAX_ROTATABLE  20
#define MAX_ROTATABLE   12
#define ROTBOND_OVERLAP 10
#define MAX_INTERNAL_ENERGY  100.0
/* <-- v.5.2.1 */
/* <--v.5.1 */

/* v.3.2--> */
/* Non-bonded contact types */
/* v.5.3 --> */
#define NOT_DEFINED     -1
/* <--v.5.3 */
#define HYDROPHOBIC_ONLY 0
#define HYDROPHOBIC_ANY  1
#define ANY_ANY          2
/* <--v.3.2 */

/* v.4.1.1--> */
/* Interaction types */
#define MAIN_MAIN       0
#define SIDE_SIDE       1
#define MAIN_SIDE       2
/* <--v.4.1.1 */

/* v.4.2--> */
/* Record types to be picked up from the XMAS file */
#define REC_ATOM        0
#define REC_CONECT      1
#define NTYPES          2
/* <--v.4.2 */

/* Loop maxima */
#define MAXCYCLE       50
#define START_SWIVELS  (MAXCYCLE / 3)

/* Array parameters */
/* v.4.1--> */
/* #define COL_NAME_LEN   15 */
#define PRM_COL_NAME_LEN 15
#define COL_NAME_LEN     18
/* <--v.4.1 */
#define CONNECT        20
#define CONNECTIONS    30
/* v.5.0--> */
/* #define FILENAME_LEN   80 */
#define FILENAME_LEN   1024
/* <--v.5.0 */
#define HHB            24
/* v.5.1.2--> */
/* #define LINELEN       100 */
#define LINELEN       1024
/* <--v.5.1.2 */
/* v.4.0--> */
/* #define MAXATOMS     1200 */
#define MAXATOMS     5000
/* <--v.4.0 */
/* v.5.4.4 --> */
#define LARGE_STRUCTURE  1000
/* <-- v.5.4.4 */
#define MAXBONDS      600
/* v.5.3.8--> */
#define MAX_CHAINS    100
/* <--v.5.3.8 */
#define MAX_COLOURS    20
#define MAXHBONDS     800
#define MAXLIGRES    1000
#define MAXNATM      1000
#define MAXNRES       150
#define MAX_PARAMETERS  4
#define MAXRESIDUES 10000
#define MAXLABELS  (2 * MAXRESIDUES)
/* v.4.0--> */
/* #define MAXSIDESTACK (MAXATOMS * MAXBONDS) */
/* <--v.4.0 */
/* v.3.2--> */
/* v.4.0.2--> */
/* #define MAXTOKENS       8 */
#define MAXTOKENS       20
/* <--v.4.0.2 */
/* <--v.3.2 */
#define MAX_FLAT_LOOP  10
#define MAXELBONDS     50
/* v.4.5.2--> */
#define MAXPATHS      100
#define MAXPATH_LEN     6
/* <--v.4.5.2 */
#define MAXSPECIAL_RES 10
#define MAXSTACK      600
#define POINTER         3
/* v.5.1.2--> */
/* #define TITLE_LEN      80 */
#define TITLE_LEN      LINELEN
/* <--v.5.1.2 */
/* v.4.3--> */
#define TOKENLEN       20
/* <--v.4.3 */

/* Limits for using anchor coords directly */
#define MIN_ANCHOR_LEN 1.0
#define MAX_ANCHOR_LEN 2.5
#define MIN_ANCHOR_LEN_SQRD (MIN_ANCHOR_LEN * MIN_ANCHOR_LEN)
#define MAX_ANCHOR_LEN_SQRD (MAX_ANCHOR_LEN * MAX_ANCHOR_LEN)

/* v.5.1--> */
#define RING_ANCHOR_WEIGHT  30.0
/* <--v.5.1 */

/* Minimization parameters */
/* v.4.0--> */
/* #define NTERMS          5 */
#define NTERMS          7
/* <--v.4.0 */
#define MIN_ANGLE     ((float) 10.0)
#define MIN_MOVE      ((float)  1.5)
/* v.4.4--> */
#define TOO_CLOSE     ((float)  1.0)
/* <--v.4.4 */

/* Miscellaneous parameters */
/* v.4.0--> */
/* #define INT_MAX       32767 */
/* <--v.4.0 */
#define SWIVEL_ANGLE   ((float) 10.0)

/* v.4.4.4--> */
/* Fixed molecule size scaling parameters */
#define CHEM_WIDTH     ((float)  10.0)
#define CHEM_HEIGHT    ((float)   7.2)
/* <--v.4.4.4 */

/* Plot parameters */
#define BBOXX1         ((float)  30.0)
#define BBOXX2         ((float) 550.0)
#define BBOXY1         ((float)  50.0)
#define BBOXY2         ((float) 780.0)
#define BORDER_MARGIN  ((float)  0.93)
#define CHAR_ASPECT    ((float)  0.484)
#define SPOKE_EXTENT   ((float)  1.4)
#define SPREAD_FACTOR  ((float)  1.0)
#define TITLE_FRACTION ((float)  0.04)
#define TITLE_SIZE     ((float)  20.0)

/* Parameter arrays */
/* v.3.2--> */
/* #define PICTURE_OPTIONS 2
   #define INCLUDE_PARAMETERS      19 */
#define PICTURE_OPTIONS          3
/* v.4.0--> */
/* #define INCLUDE_PARAMETERS      21 */
/* v.4.3--> */
/* #define INCLUDE_PARAMETERS      23 */
#define INCLUDE_PARAMETERS      24
/* <--v.4.3 */
/* <--v.4.0 */
/* <--v.3.2 */
#define OBJECT_COLOURS          18
#define OBJECT_SIZES             9
#define TEXT_COLOURS             9
#define TEXT_SIZES               8
/* v.4.0--> */
/* #define MINIMIZATION_PARAMETERS 11 */
#define MINIMIZATION_PARAMETERS 15
/* <--v.4.0 */

/* Parameters defining layout of Key to symbols */
#define KEY_SCALE          ((float)  0.15)

/* Hydrophobic symbol parameters */
#define HYDROPHOBIC_RADIUS ((float)  1.1)
#define SPOKE_ANGLE        ((float) 12.0)

#define NON_FLAT           ((float)  0.2)
#define KEY_HEIGHT         ((float)  0.14)
#define MARGIN_ACC         ((float)  0.04)
#define MARGIN_NO_HPH      ((float)  0.04)
#define IDEAL_ENERGY       ((float) 10.0)
/* v.5.3.6 --> */
#define LABEL_GAP          ((float)  0.2)
/* <-- v.5.3.6 */
#define SWING_STEPS        45
#define MOVE_STEPS         30
/* v.4.4--> */
#define ABORT_ENERGY 1000000.0
/* <--v.4.4 */

#define ATOM_INTERACT_DIST      ((float) 1.0)
#define ATOM_DIST2      (ATOM_INTERACT_DIST * ATOM_INTERACT_DIST)
#define BOND_INTERACT_DIST      ((float) 0.6)
#define INTERACT_DIST2  (BOND_INTERACT_DIST * BOND_INTERACT_DIST)
#define INTERNAL_INTERACT_DIST  ((float) 0.3)
#define INTERNAL_DIST2  (INTERNAL_INTERACT_DIST * INTERNAL_INTERACT_DIST)
#define PS_STRING_LENGTH               132
#define MAXH                             4
#define HSCALE_FACTOR           ((float) 0.0)

/* Angle constants */
#define PI         ((float)  3.141592654)
#define RADDEG     ((float)  (180.0 / PI))
#define PI_BY_2    ((float)  (PI / 2.0))

/* debug */
int debug;
/* debug */

/* v.4.3--> */
/* Line-types in .drw file */
#define GLOBAL_PARAMS   0
#define COLOURMAP_DATA  1
#define COLOUR_DATA     2
#define SIZE_DATA       3
#define OBJECT_DATA     4
#define RES_DATA        5
#define ATOM_DATA       6
#define BOND_DATA       7
#define PLOT_TITLE      8
static char hash_codes[] = "PMCSORABT";
/* <--v.4.3 */

/* v.5.0--> */
/* Constants */
#define BOLTZMANN_CONST ((float) 0.1)
#define STAGE1        500
#define STAGE2       1000
/* v.5.1--> */
#define STAGE3       1500
#define STAGE4       2000
/* <--v.5.1 */

/* v.5.2 --> */
/* Residue similarities */
#define RES_IDENTICAL         0
#define RES_SIMILAR           1
#define RES_DIFFERENT         2
#define RES_MISSING           3
#define RES_ANY               4

#define NRES_SIMTYPES         5

static double SIMILARITY_COLOUR[NRES_SIMTYPES][3] = {
  {1.000,  0.412,  0.706},  /* pink */
  {1.000,  0.640,  0.000},  /* orange */
  {0.500,  0.500,  0.500},  /* grey */
  {1.000,  1.000,  1.000},  /* white */
  {1.000,  1.000,  1.000}   /* white */
};

static char *SIMILARITY_DESC[NRES_SIMTYPES] = {
  "Identical_colour",
  "Similar_colour  ",
  "Different_colour",
  "Missing_colour  ",
  "Any_colour      ",
};

/* Interaction types */
/* v.5.4--> */
#define INVALID              -1
/* <--v.5.4 */
#define S_HBONDS              0
#define S_CONTACTS            1
#define S_SALT_BRIDGES        2
#define S_DISULPHIDES         3
/* v.5.4--> */
#define S_COVALENT            4

//#define S_NTYPES              4
#define S_NTYPES              5

/* Residue charges */
#define POSITIVE              1                  
#define NEGATIVE              2                  

/* Distance parameters */
#define SALT_BRIDGE_DIST    4.0
#define SBRIDGE_DIST_SQRD   (SALT_BRIDGE_DIST * SALT_BRIDGE_DIST)
#define DISULPHIDE_DIST      2.5
#define DISULPH_DIST_SQRD   (DISULPHIDE_DIST * DISULPHIDE_DIST)
/* <--v.5.4 */

//static char S_INT_TYPE[S_NTYPES] = { 'H', 'N', 'S', 'D' };
/* <-- v.5.2 */

/* v.5.3.4--> */
/* Ideal gap between the interfaces in a DIMPLOT diagram */
#define IDEAL_GAP         ( 3.0 )
#define IDEAL_SEPARATION  ( 0.4 )
#define MIN_OVERLAP       ( 1.0 )
#define LEFT_FACE             0
#define RIGHT_FACE            1
#define INTERACTING_WATER     2
/* <--v.5.3.4 */

/* v.5.3.8--> */
/* Ellipse paramters for highlighted residue */
#define ANGLE_STEP           10
#define ELLIPSE_MARGIN    ( 1.0 )
#define ELLIPSE_LINEWIDTH ( 0.1 )
#define NEAR_ZERO         ( 0.000001 )

/* The DisaStr variant types */
#define NO_VARIANT           -2
#define NAT_VARIANT          -1
#define DISEASE               0
#define DISEASE_NOTE          1
#define MUTAGEN_EXPT          2
#define DISEASE_DDD           3
#define CLINVAR_DATA          4
#define ADDED_VAR             5

#define NVAR_TYPES            6

static char *ELLIPSE_COL_DEF[NVAR_TYPES + 2] = { 
  "none", "natvar", "disease", "note", "mutagen", "DDD", "clinvar",
  "added"
};

static char *ELLIPSE_COLOUR[NVAR_TYPES + 2] = {
  "0.0 0.0 0.0", "0.0 0.0 0.0", "1.0 0.0 0.0", "0.7 0.7 0.7",
  "0.0 1.0 0.0", "0.0 0.0 1.0", "1.0 0.0 0.0", "0.5 0.0 1.0"
};

/* Interaction type */
#define NOT_WANTED           -9
#define TO_PROTEIN            0
#define TO_DNA                1
#define TO_LIGAND             2
#define TO_METAL              3
/* <--v.5.3.8 */

/* Working directory */
char wkdir[FILENAME_LEN];
char ligplot_prm[FILENAME_LEN];
char ligplot_drw_name[FILENAME_LEN];
int nensemble;
/* v.5.3.8--> */
char highlight_res[16];
char pdb_code[5];
char uniprot_acc[20];
char pdbsum_dir[FILENAME_LEN];
/* <--v.5.3.8 */
int ligfit;
/* v.5.3.4--> */
int dimfit;
/* <--v.5.3.4 */
/* v.5.2 --> */
char ligsearch_id[FILENAME_LEN];
int ligsearch_annot;
/* v.5.3 --> */
int contact_type;
/* <--v.5.3 */
/* <-- v.5.2 */
/* <--v.5.0 */

/* v.5.0.1--> */
int conect;
/* <--v.5.0.1 */

/* File pointers */
FILE *fil_bonds;
FILE *fil_hbplus;
FILE *fil_pdb;
FILE *ligplot_bonds_out;
FILE *ligplot_frm;
FILE *ligplot_hhb_out;
FILE *ligplot_nnb_out;
FILE *ligplot_par;
FILE *ligplot_pdb;
/* v.3.2--> */
FILE *ligplot_res;
/* <--v.3.2 */
/* v.4.0--> */
FILE *ligplot_rcm;
/* <--v.4.0 */
/* v.4.1--> */
FILE *ligplot_drw;
/* <--v.4.1 */
/* v.4.1.1--> */
FILE *ligplot_sum;
/* <--v.4.1.1 */

FILE *ps_file;

/* Miscellaneous counters, flags and arrays */
int Nobjects;
int Water_as_Ligand;
/* v.3.1.2--> */
int Metal_as_Ligand;
/* <--v.3.1.2 */
/* v.3.2--> */
int Write_Res_File;
/* <--v.3.2 */
/* v.4.3--> */
int Write_Lig_File;
/* <--v.4.3 */
int xsite_file;
char xsite_probe[7];
/* v.4.0--> */
int Interface_Plot;
/* <--v.4.0 */
/* v.5.3.4--> */
int From_Dimer;
/* <--v.5.3.4 */
/* v.4.3--> */
int Plot_from_drw;
/* <--v.4.3 */
/* v.4.4.4--> */
int Scale_Plot;
/* <--v.4.4.4 */
/* v.4.5.4--> */
int Write_Bond_Types;
/* <--v.4.5.4 */
/* v.5.4.5 --> */
int PDBsum1;
/* <-- v.5.4.5 */

float Maximum_accessibility, label_coord[MAXLABELS][2];
int n_labels;

char root_name[TITLE_LEN];

char lig_res_name_store[MAXLIGRES][4];
char lig_res_num_store[MAXLIGRES][6];
char lig_chain_store[MAXLIGRES];
int ligand_residues;
/* v.4.1--> */
int ncols;
/* <--v.4.1 */

/* v.5.0--> */
/* Structure for parameter keys (from LigPlus) and mappings (to LIGPLOT)
   --------------------------------------------------------------------- */
  struct pkey
  {
    char                *key;
    int                 group;
    int                 iparam;
  };

#define NPARAM_KEYS     2

/* Parameter key mappings (LigPlus to LIGPLOT */
  static const struct pkey PARAM_KEY[NPARAM_KEYS] = {
    { "IN_COLOUR", 0, 0 },
    { "PORTRAIT", 0, 1 }
  };

/* Structure for input parameters from LigPlot+
   -------------------------------------------- */
struct parameter
{
  char                *key;
  char                *value;
  struct parameter    *next_parameter_ptr;
};

/* Structure for ensemble (when plotting from LigPlot+)
   ---------------------------------------------------- */
struct ensemble
{
  struct ensemble    *next_ensemble_ptr;
};
struct ensemble *first_ensemble_ptr;
/* <--v.5.0 */

/* Structure for information read in from HBPLUS files
   --------------------------------------------------- */
struct hhb_info
{
  int                source;
/* v.5.4--> */
  int                bond_type;
/* <--v.5.4 */
  int                atom_number1;
  char               atom_name1[5];
  char               res_name1[4];
  char               res_num1[6];
  char               chain1;
  float              x1, y1, z1;
  int                atom_number2;
  char               atom_name2[5];
  char               res_name2[4];
  char               res_num2[6];
  char               chain2;
  float              x2, y2, z2;
  float              bond_length;
  struct hhb_info    *next_hhb_info_ptr;
};
struct hhb_info *first_hhb_info_ptr, *last_hhb_info_ptr;

/* Structure for objects to be drawn on the final LIGPLOT picture
   -------------------------------------------------------------- */
struct object
{
  int                object_type;
  int                nresidues;
  int                nbonds;
  int                nrot_bonds;
  int                nmain_chain;
  int                n_ca;
  int                weight;
  float              minx, miny, maxx, maxy;
/* v.4.0--> */
  int                interface;
  int                have_anchors;
/* <--v.4.0 */
/* v.5.0--> */
  int                fixed;
/* <--v.5.0 */
/* v.5.3.4--> */
  int                coords_fixed;
/* <--v.5.3.4 */
  float              max_atom_size;
  float              place_x, place_y;
  float              internal_energy;
  float              total_energy;
  struct residue     *first_residue_ptr;
  struct object_bond *first_object_bond_ptr;
  struct object      *next_object_ptr;
};
struct object *first_object_ptr;

/* Structure for residue data read in from PDB file
   ------------------------------------------------ */
struct residue
{
  char               res_name[4];
  char               res_num[6];
  char               chain;
  int                residue_type;
  int                natoms;
  int                inligand;
  float              minx, miny, maxx, maxy;
  float              max_atom_size;
/* v.3.2--> */
  int                deleted;
/* <--v.3.2 */
/* v.4.0--> */
  float              original_mean_x, original_mean_y, original_mean_z;
  float              flattened_mean_x, flattened_mean_y;
  int                have_anchor;
  float              anchor_pstn_x, anchor_pstn_y;
  int                nanchored_atoms;
/* <--v.4.0 */
/* v.4.1--> */
  float              label_x;
  float              label_y;
/* <--v.4.1 */
  struct coordinate  *first_atom_ptr;
  struct object      *object_ptr;
  struct residue     *next_residue_ptr;
};
struct residue *first_residue_ptr;

/* Structure for atom data
   ----------------------- */
struct coordinate
{
  char               atom_number[6];
  char               atom_name[5];
  char               print_name[5];
/* v.4.3--> */
  char               element[3];
/* <--v.4.3 */
  struct residue     *residue_ptr;
  float              original_x, original_y, original_z;
  float              x, y, z;
  float              fit_x, fit_y, fit_z;
  float              save_x, save_y;
  char               occupancy[7];
  char               bvalue[7];
  float              accessibility;
  float              atom_size;
  int                side_chain;
  int                checked;
  int                deleted;
  int                natom_links;
  int                plot_atom;
/* v.4.0--> */
  int                have_anchor;
  float              anchor_pstn_x, anchor_pstn_y;
/* <--v.4.0 */
/* v.5.0--> */
  float              anchor_weight;
/* <--v.5.0 */
/* v.5.1--> */
  int                ring_atom;
/* <--v.5.1 */
/* v.4.1--> */
  int                atom_count;
  float              label_x;
  float              label_y;
/* <--v.4.1 */
/* v.5.2 --> */
  int                score[S_NTYPES];
  int                have_score;
/* <-- v.5.2 */
  struct coordinate  *next_stack_ptr;
  struct atom_link   *first_atom_link_ptr;
  struct coordinate  *next;
};
struct coordinate *first_atom_ptr;

/* Structure containing info about each bond
   ---------------------------------------- */
struct bond
{
  struct coordinate  *first_atom_ptr;
  struct coordinate  *second_atom_ptr;
  int                bond_type;
  int                elastic;
  int                rotatable_bond;
  int                ring_bond;
  int                end_bond;
  int                bond_order;
  int                bond_source;
  float              bond_length;
  int                checked;
  int                flattened;
  int                nfirst_atom_links;
  int                nsecond_atom_links;
  int                nbonds_from_first_atom;
  int                nbonds_from_second_atom;
/* v.4.1--> */
  char               label[5];
/* <--v.4.1 */
/* v.4.4--> */
  int                plotted;
/* <--v.4.4 */
/* v.4.5.2--> */
  int                ring_no;
/* <--v.4.5.2 */
  struct bond_link   *first_bond_link_ptr;
  struct bond        *next_stack_ptr;
  struct bond        *next_link_ptr;
  struct bond        *next_bond_ptr;
};
struct bond *first_bond_ptr;


/* Structure giving atoms covalently bonded to current atom
   -------------------------------------------------------- */
struct atom_link
{
  struct coordinate  *atom_ptr;
  struct atom_link   *next_atom_link_ptr;
};

/* Structure giving bonds linked to either end of each bond
   -------------------------------------------------------- */
struct bond_link
{
  struct bond        *bond_ptr;
  int                bond_end;
  struct bond_link   *next_bond_link_ptr;
};

/* Structure giving bonds belonging to each object
   ----------------------------------------------- */
struct object_bond
{
  struct bond        *bond_ptr;
  struct object_bond *next_object_bond_ptr;
};

/* v.4.1.1--> */
/* Structure giving list of interacting residue pairs and the
   interactions between them
   ---------------------------------------------------------- */
struct pair
{
  struct residue     *residue1_ptr;
  struct residue     *residue2_ptr;
  int                nhbonds[3];
  int                ncontacts;
  int                duplicate;
  struct pair        *next_pair_ptr;
};

/* <--v.4.1.1 */

/* v.4.5.2--> */
/* Structure for tracing a path from a given bond
   ---------------------------------------------- */
struct path
{
  struct bond        *bond_ptr;
  int                end;
  struct path        *next_path_ptr;
};
/* <--v.4.5.2 */

/* Plot parameters defined by user
   ------------------------------- */
struct Picture_Options
{
int In_Colour, Portrait;
/* v.3.2--> */
float Rotation_Angle;
/* <--v.3.2 */
/* v.5.3.4--> */
int Allow_Rotation;
/* <--v.5.3.4 */
};
struct Picture_Options *Picture;

struct Include_Parameters
{
int Hydrophobics, Waters, Mainchain_Atoms, Linked_Residues, Hbonds,
    Internal_Hbonds, External_Bonds, Hydrophobic_Bonds,
    Simple_Ligand_Residues, Simple_Nonligand_Residues, Accessibilities,
    Ligand_Atoms, Nonligand_Atoms, Double_Bonds, Key,
    Residue_Names, Atom_Names, Hbond_Lengths, Filename_for_Title;
/* v.3.2--> */
int External_Bonds_Solid, Contact_Type;
/* v.4.0--> */
int Ligand_Accessibilities_Only, Water_Atoms;
/* <--v.4.0 */
/* <--v.3.2 */
/* v.4.3--> */
int Chemical_Notation;
/* <--v.4.3 */
};
struct Include_Parameters *Include;

int Include_Hydrogens;
/* v.4.4.1--> */
int Greek_Names;
/* <--v.4.4.1 */

/* v.5.0--> */
int Include_Waters;
/* <--v.5.0 */

/* v.5.1--> */
int Filter_Waters;
/* <--v.5.1 */

/* v.5.1.4--> */
int Clip_Diagram;
/* <--v.5.1.4 */

struct Sizes
{
char Ligand_Atoms[COL_NAME_LEN + 1],
     Nonligand_Atoms[COL_NAME_LEN + 1],
     Waters[COL_NAME_LEN + 1],
     Hydrophobics[COL_NAME_LEN + 1],
     Simple_Residues[COL_NAME_LEN + 1],
     Ligand_Bonds[COL_NAME_LEN + 1],
     Nonligand_Bonds[COL_NAME_LEN + 1],
     Hydrogen_Bonds[COL_NAME_LEN + 1],
     External_Bonds[COL_NAME_LEN + 1];
};
struct Sizes *Size;

struct Size_Vals
{
float Ligand_Atoms, Nonligand_Atoms, Waters, Hydrophobics, Simple_Residues,
      Ligand_Bonds, Nonligand_Bonds, Hydrogen_Bonds, External_Bonds;
};
struct Size_Vals *Size_Val;

/* v.4.2--> */
/* float object_size[OBJECT_SIZES]; */
float object_size[OBJECT_SIZES], param_object_size[OBJECT_SIZES];
/* <--v.4.2 */
/* v.4.1--> */
char object_size_desc[OBJECT_SIZES][COL_NAME_LEN + 1];
/* <--v.4.1 */

struct Text_Sizes
{
char Ligand_Residue_Names[COL_NAME_LEN + 1],
     Nonligand_Residue_Names[COL_NAME_LEN + 1],
     Water_Names[COL_NAME_LEN + 1],
     Hydrophobic_Names[COL_NAME_LEN + 1],
     Simple_Residue_Names[COL_NAME_LEN + 1],
     Ligand_Atom_Names[COL_NAME_LEN + 1],
     Nonligand_Atom_Names[COL_NAME_LEN + 1],
     Hbond_Lengths[COL_NAME_LEN + 1];
};
struct Text_Sizes *Text_Size;

struct Text_Size_Vals
{
float Ligand_Residue_Names, Nonligand_Residue_Names, Water_Names,
      Hydrophobic_Names, Simple_Residue_Names, Ligand_Atom_Names, 
      Nonligand_Atom_Names, Hbond_Lengths;
};
struct Text_Size_Vals *Text_Size_Val;

/* v.4.2--> */
/* float text_size[TEXT_SIZES]; */
float text_size[TEXT_SIZES], param_text_size[TEXT_SIZES];
/* <--v.4.2 */
/* v.4.1--> */
char text_size_desc[TEXT_SIZES][COL_NAME_LEN + 1];
/* <--v.4.1 */

struct Colours
{
char Background[COL_NAME_LEN + 1], Ligand_Bonds[COL_NAME_LEN + 1],
    Nonligand_Bonds[COL_NAME_LEN + 1], Hydrogen_Bonds[COL_NAME_LEN + 1],
    External_Bonds[COL_NAME_LEN + 1], Hydrophobics[COL_NAME_LEN + 1],
    Accessibility_Min[COL_NAME_LEN + 1], Accessibility_Max[COL_NAME_LEN + 1],
    Nitrogen[COL_NAME_LEN + 1], Oxygen[COL_NAME_LEN + 1],
    Carbon[COL_NAME_LEN + 1], Sulphur[COL_NAME_LEN + 1],
    Water[COL_NAME_LEN + 1], Phosphorus[COL_NAME_LEN + 1],
    Iron[COL_NAME_LEN + 1], Other[COL_NAME_LEN + 1],
    Atom_Edges[COL_NAME_LEN + 1], Simple_Residues[COL_NAME_LEN + 1];
};
struct Colours *Colour;

char object_colour[OBJECT_COLOURS][COL_NAME_LEN + 1];
/* v.4.1--> */
int object_grey[OBJECT_COLOURS];
int object_inverse[OBJECT_COLOURS];
/* <--v.4.1 */

struct Text_Colours
{
char Title[COL_NAME_LEN + 1],
    Key_Text[COL_NAME_LEN + 1],
    Ligand_Residue_Names[COL_NAME_LEN + 1],
    Nonligand_Residue_Names[COL_NAME_LEN + 1],
    Water_Names[COL_NAME_LEN + 1],
    Hydrophobic_Names[COL_NAME_LEN + 1],
    Ligand_Atom_Names[COL_NAME_LEN + 1],
    Nonligand_Atom_Names[COL_NAME_LEN + 1],
    Hbond_Lengths[COL_NAME_LEN + 1];
};
struct Text_Colours *Text_Colour;

/* v.5.3.8--> */
/* Structure for list of residues to be highlighted on the plot
   ------------------------------------------------------------ */
struct highres
{
  int                     type;
  int                     ref;
  char                    res_name[4];
  char                    res_num[6];
  char                    chain;
  struct highres          *next_highres_ptr;
};
/* <--v.5.3.8 */

/* Various global parameters */
/* v.5.0--> */
int run_type;
/* <--v.5.0 */
char text_colour[TEXT_COLOURS][COL_NAME_LEN + 1];
/* v.4.1--> */
int text_grey[TEXT_COLOURS];
int text_inverse[TEXT_COLOURS];
/* <--v.4.1 */

char Colour_Table_Name[MAX_COLOURS][COL_NAME_LEN + 1];
float Colour_Table[MAX_COLOURS][3];

char special_res[MAXSPECIAL_RES][15];

char Print_Title[TITLE_LEN + 1];
/* v.5.4.5 --> */
char Drw_Title[TITLE_LEN + 1];
/* <-- v.5.4.5 */

/* Minimization parameters */
int Max_Loops, Random_Start;
float Atom_Atom_Clash, Bond_Atom_Clash, Bond_Stretch, Energy_Drop_Maximum;
float HB_Weight, Internal_Energy_Weight, Max_HBdist, Nonbond_Weight;
float Overlap_Score;
/* v.4.0--> */
float Anchor_Weight, Boundary_Energy_Weight, Min_Boundary_Dist;
float Rel_Pstn_Energy_Weight;
/* <--v.4.0 */
/* v.5.3.4--> */
float Closeup_Energy_Weight;
/* <--v.5.3.4 */

/* Global Plot variables
   --------------------- */
int Nwarnings;
int Print_as_is, Split_Colour_Ligand_Bonds, Split_Colour_Nonligand_Bonds;
float Margin_y;
float Page_Max_x, Page_Min_x, Page_Max_y, Page_Min_y;
float Plot_Centre_x, Plot_Centre_y;
float Accessibility_Max[3], Accessibility_Min[3], Max_atom_radius,
      Max_object_size, Max_object_height, Max_object_width,
      Max_Textsize, Max_Textwidth;
/* v.5.0--> */
int verbose;
/* <--v.5.0 */
/* v.4.0--> */
int Have_Anchors;
/* <--v.4.0 */
/* v.4.2--> */
int Ligand_Molecule_id;
int Output_drw_File_Only, Have_Xmas, Plot_All_Xmas_Ligands, Read_Xmas;
/* <--v.4.2 */

/*************************************************************************

Include variables for Andrew Martin's least-squares fitting routines,
matfit and qikfit

*************************************************************************/

#define SMALL  1.0e-20     /* Convergence cutoffs                       */
#define SMALSN 1.0e-10

typedef short  BOOL;
typedef double REAL;
typedef struct
{  REAL x, y, z;
}  VEC3F;

typedef VEC3F COOR;

#define ABS(x)   (((x)<0)   ? (-(x)) : (x))

/* Prototypes for functions defined in fit.c */
/* BOOL matfit(COOR *x1, COOR *x2, REAL rm[3][3], int n, REAL *wt1, 
            BOOL column);

static void qikfit(REAL umat[3][3], REAL rm[3][3], BOOL column);  */

/***********************************************************************

Include variables for the Genetic Algorithm routines

***********************************************************************/

/* GA parameters */
#define EPSILON 1e-40
#define MAXFUN 1000000  // Upper bound for number of function evaluations
                        // in case best fitness doesn't reach the desired
                        // value
#define MAXP 100        // Population size (best if chosen between 100-150)
#define KIDS 2          // Pool size of kids to be formed (use 2,3 or 4)
#define M 1             // M+2 is the number of parents participating in
                        // xover (use 1)
#define FAMILY 2        // Number of parents to be replaced by good
                        // individuals(use 1 or 2)
#define SIGMA_ZETA 0.1
#define SIGMA_ETA 0.1   // Variances used in PCX (best if fixed at these
                        // values)

/* Other run parameters */
#define MINIMIZE 1      // Set 1 to minimize and -1 to maximize number
#define RANDPARENT M+2  // of parents participating in PCX

/* Can change these */
#define MAXRUN 10       // Number of runs each with different random
                        // initial population
#define LIMIT 100.0     // Accuracy of best solution fitness desired


/* Structure for each population
   ----------------------------- */
struct pop
{
  double *vari;
  double obj;
};

/**************************************************************/


