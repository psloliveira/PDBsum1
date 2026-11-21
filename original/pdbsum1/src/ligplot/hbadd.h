/**************************************************************************
 *
 * hbadd.h - Include file for program hbadd.c
 *
 *************************************************************************/

#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <ctype.h>

/* Flag-values */
#define FALSE           0
#define TRUE            1
#define UNSET           2

static char *T_F[3] = { "FALSE", "TRUE", "*UNSET*" };

/* Array parameters */
#define COLNAM_LEN     11
/* v.5.0--> */
/* #define FILENAME_LEN  130         Max filename-length */
/* #define LINE_LEN      130         Max line-length */
#define FILENAME_LEN  512         /* Max filename-length */
#define LINE_LEN     1026         /* Max line-length */
/* <--v.5.0 */
#define MAXCONNECT      9         /* Max no. of connections per atom */
#define MAXABORT        5         /* Max no. of aborts on graph-matching */
#define MAXANGLES      10         /* Max no. of bond angles per atom */
#define MAX_INT_CHAINS 10
#define MAXPATHS     1000         /* Max no. of paths to generate before 
                                     giving up */
/* v.5.1--> */
#define MAX_TOKEN     100
#define MAX_TYPES     100

/* Scale factor for coords in input .sdf file */
#define MOL_SCALE_FACTOR    1.0
/* <--v.5.1 */

/* v.5.0.1--> */
/* #define MIN_MAPPED      5 */
#define MIN_MAPPED      3
#define MAX_MISSING     4
#define MIN_PERCENT ( 40.0 )
/* <--v.5.0.1 */
#define REF_LEN        20
#define RESNAM_LEN      8         /* Maximum length for dic residue name */
#define SHOWPATHS     100
#define TOKEN_LEN     100

/* Directory separator */
static char DIR_SEP[] = "/";

/* NULL character */
static char null = '\0';

/* Coordinate subscripts */
#define X                     0
#define Y                     1
#define Z                     2

/* Record types defined in the rasmol.dfn file */
#define OTHER                -1
#define CA_ONLY               0
#define LIGAND                1
#define METAL                 2
#define NUCLEIC_ACID          3
#define PROTEIN               4
#define LIGAND_RANGE          5

/* Miscellaneous parameters */
#define COVAL_DIST    2.0
#define COVAL_DIST2 (COVAL_DIST * COVAL_DIST)
#define HATOM_DIST    1.42
#define HATOM_DIST2 (HATOM_DIST * HATOM_DIST)
#ifdef INT_MAX
#else
#define INT_MAX     32767
#endif

/* Angle constants */
#define  PI          3.141592654
#define  RADDEG      (180.0 / PI)
#define  PI_BY_2     (PI / 2.0)
#define  ANGLE_ERROR    4.0       /* Allowed bond-angle error */

/* Dictionary types */
#define TEXT            0
#define CIF             1

/* Bond order types */
#define SINGLE          1
#define DOUBLE          2
#define TRIPLE          3
#define AROMATIC        4

/* Match-scores for comparing subgraphs */
/* v.5.0--> */
/* #define MATCH_SCORE     2
   #define MAP_SCORE       1 */
#define MATCH_SCORE     5
#define MAP_SCORE       4

#define VERBOSE         FALSE
/* <--v.5.0 */

static char *bond_type[] = { "UNKN", "SING", "DOUB", "TRIP", "AROM" };

/* Error flag */
int error_flag;

/* v.5.0--> */
/* Match and map scores */
int Match_Score;
/* <--v.5.0 */

/* File pointers */
FILE *fil_dic;
FILE *fil_map;
FILE *fil_out;
FILE *fil_pdb;
FILE *fil_sum;

/* Structure for residue data read in from PDB file
   ------------------------------------------------ */
struct residue
{
  char              res_name[4];
  char              res_num[6];
  char              chain;
  int               natoms;
/* v.5.0--> */
  int               bonds;
/* <--v.5.0 */
  int               nconnect;
  int               non_hydrogen;
  int               nmatch;
  int               fake;
  int               nmatched;
  int               nmapped;
  int               nmissed;
  int               nmissing;
  int               nbonds_matched;
  int               liguniq_no;
  int               ligand_no;
  struct atom       *first_atom_ptr;
  struct dictionary *dictionary_ptr;
  struct rasdef     *rasdef_ptr;
  struct residue    *next_residue_ptr;
};

/* Structure for atom coordinates data read in from PDB file
   --------------------------------------------------------- */
struct atom
{
  int               atom_number;
  char              atom_name[5];
  char              element[3];
  int               hydrogen;
/* v.5.1--> */
  int               type;
/* <--v.5.1 */
  float             x, y, z;
  struct residue    *residue_ptr;
  struct dic_data   *dic_data_ptr;
  int               checked;
  int               nequiv;
  struct atom       *connect_atom_ptr[MAXCONNECT];
  struct atom       *next_atom_ptr;
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
  char               interacts_with[MAX_INT_CHAINS];
  struct rasdef      *next_rasdef_ptr;
};

/* Structure for atom connectivities, as read in from PDB file
   ----------------------------------------------------------- */
struct connect
{
  int               atom_number1;
  int               atom_number2;
  struct connect    *next_connect_ptr;
};

/* Structure for HET group entries from HET Group Dictionary
   --------------------------------------------------------- */
struct dictionary
{
  char              res_name[RESNAM_LEN];
  char              *mol_name;
  int               natoms;
/* v.5.0--> */
  int               nbonds;
/* <--v.5.0 */
  int               nmatch;
  int               non_hydrogen;
  int               have_hydrogens;
  int               fake_entry;
  int               nconnect;
  struct residue    *residue_ptr;
  struct dic_data   *first_dic_data_ptr;
  struct dictionary *next_dictionary_ptr;
};


/* Structure for atom data read in from HET Group Dictionary
   --------------------------------------------------------- */
struct dic_data
{
  char              atom_name[5];
  char              output_atom_name[5];
  char              element[3];
  int               hydrogen;
  int               number;
  int               nbound;
  char              boundto[MAXCONNECT][5];   
  int               boundto_num[MAXCONNECT];
  int               bond_order[MAXCONNECT];
  int               checked;
  struct atom       *atom_ptr;
  struct dic_data   *connect_atom_ptr[MAXCONNECT];
  struct dic_data   *next_dic_data_ptr;
};

/* Structure for possible atom equivalences
   ---------------------------------------- */
struct equiv_atom
{
  int                equiv_atom_number;
  int                checked;
  int                in_best_clique;
  int                nlinks;
  int                tree_level;
  int                match_score;
  int                links_combination;
  struct equiv_bond  *first_equiv_bond_ptr;
  struct atom        *atom_ptr;
  struct dic_data    *dic_data_ptr;
  struct equiv_atom  *next_equiv_atom_ptr;
};

/* Structure for links between graph nodes (ie possibly equivalenced bonds
   given possible atom equivalences)
   ----------------------------------------------------------------------- */
struct equiv_bond
{
  int                link_number;
  int                possible_start;
  int                checked;
  int                on_stack;
  int                match_score;
  struct equiv_atom  *equiv_atom_ptr;
  struct equiv_atom  *other_equiv_atom_ptr;
  struct equiv_bond  *reverse_equiv_bond_ptr;
  struct equiv_bond  *ea_next_equiv_bond_ptr;
  struct equiv_bond  *next_equiv_bond_ptr;
};

/* Structure for paths tree
   ------------------------ */
struct tree
{
  int                branches_checked;
  int                score;
  struct equiv_bond  *equiv_bond_ptr;
  struct tree        *prev_level_tree_ptr;
  struct tree        *next_level_tree_ptr;
  struct tree        *link_tree_ptr;
  struct tree        *first_end_tree_ptr;
  struct tree        *next_end_tree_ptr;
  struct tree        *prev_tree_ptr;
  struct tree        *next_tree_ptr;
  struct tree        *all_next_tree_ptr;
};

/* v.5.1--> */
/* Structure for each bond record
   ------------------------------ */
struct bond
{
  int                     number;
  int                     type;
  struct atom             *atom_ptr[2];
  struct bond             *next_bond_ptr;
};

/* <--v.5.1 */
