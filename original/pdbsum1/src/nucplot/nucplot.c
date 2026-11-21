/**************************************************************************
 *                               NUCPLOT
 *
 * Program to read a PDB file for a protein-DNA complex and draw a
 * schematic diagram of the interactions.
 * 
 * Written by Nicholas M Luscombe. Modified by Roman A Laskowski.
 * Biomolecular Structures & Modelling Unit
 * Department of Biochemistry & Molecular Biology
 * University College London
 * Gower Street, London WC1E 6BT, UK
 *
 * Written 18 November 1996
 *
 * Usage  : nucplot
 * Input  : pdb format file of a protein-DNA complex
 *          .hb2 format HBPLUS output file between protein and DNA:
 * 
 *
 * 3. Modify the links for the list of DNA strands
 *    Link the paired strands together and break all other links
 *    beg_chain[] points to first atom of first strand
 *    second_strand[] points to the first atom of the second strand
 *    end_chain[] points to the last atom of the second strand
 *
 * 4. Store pairing information in an array of type PAIRS:
 *    beg_pair[] points to the first pair
 *    end_pair[] points to the last pair
 *
 * 5. Drawing:
 *    Calculate scaling and positioning 
 *    Draw nucleic acids (select by choosing parameters:
 *       text version
 *       box version
 *
 * Modifications
 * -------------
 *
 * v.1.1   Fixes required for certain compilers
 *                                                      (RAL 14 Aug 2001)
 * v.1.1.1 Bug fixes for bridge-interactions to last phosphate on each page.
 *         Errors on plot when waters labelled as belonging to more than
 *         one chain (as in 1qna)
 *                                                      (RAL  3 Mar 2005)
 *
 * v.1.1.2 Additional info on H-bond warning. Fix to bug introduced in
 *         v.1.1.1 concerning bridging waters of chain blank.
 *         Fix to sep_bridged_waters which was failing to put bridging
 *         waters at the end of the interactions list.
 *                                                      (RAL 24 Mar 2005)
 *
 * v.1.1.3 Change to list of acceptable nucleic acid names to match new
 *         format introduced by PDB remediation project.
 *                                                      (RAL  6 Aug 2007)
 *
 * v.1.1.4 Fix to deal with first residue number being -1.
 *                                                      (RAL  7 Nov 2010)
 * v.1.1.5 Various minor bug-fixes to unitialised variables.
 *                                                      (RAL  3 Jul 2012)
 * v.1.1.6 Minor bug-fixes to reset dashed-line type after key.
 *                                                      (RAL 20 Aug 2013)
 * v.1.1.7 Addition of plot of red ellipse around selected residue.
 *                                                      (RAL 19 Feb 2019)
 * v.1.1.8 Additions for running in PDBsum1.
 *                                                      (RAL 23 Mar 2022)
 * v.1.1.9 Addition of "DA.", "DC.", etc base names.
 *                                                      (RAL 31 May 2023)
 *
 * v.1.1.10 Bug fix atom number read and terminal interactions lost
 *          from plot.
 *                                                      (RAL 15 Jan 2024)
 *
 * Compilation
 * -----------
 *
 * cc -o nucplot nucplot.c -lm
 *
 * Optimization flags:-  cc -O2 -funroll-loops -o nucplot nucplot.c -lm
 *
 *
 * Function calling tree
 *
 * main
 *    -> initialize_parameters
 *          -> verify_colours
 *                -> verify_colour_name
 *    -> read_parameters
 *          -> update_current_parameter
 *    -> prepare_colour_table
 *    -> verify_colours
 *          -> verify_colour_name
 *    -> com_line (command_line)
 *          -> print_usage
 *    -> read_pdb
 *          -> open_file
 *          -> is_hydrogen_atom
 *          -> get_atom_info
 *    -> split
 *          -> is_nucleic_acid
 *    -> baselink
 *    -> print 
 *    -> print_rev 
 *    -> readhbplus
 *          -> open_file
 *          -> check_hbplus
 *          -> dna_base_pair
 *          -> prot_dna_hbond
 *                -> get_atom_types
 *                      -> is_hydrogen
 *                -> find_hbond_res
 *                -> find_hbond_res
 *                -> hbond_warning
 *                -> hbond_list
 *          -> water_dna_hbond
 *                -> wd_hbond_type
 *                -> find_hbond_water
 *                -> find_hbond_res
 *                -> hbond_warning
 *                -> hbond_list
 *                -> mark_water_hbond
 *          -> prot_water_hbond
 *                -> get_atom_types
 *                      -> is_hydrogen
 *                -> find_hbond_water
 *                -> find_hbond_res
 *                -> hbond_warning
 *                -> check_water_bridge
 *                -> hbond_list
 *    -> read_connect_records
 *          -> get_connect_records
 *          -> connect_list
 *          -> check_connect_records
 *          -> mark_first_connect
 *    -> write_bonds
 *    -> paircatalog
 *          -> create_nodes
 *          -> link_nodes
 *          -> assign_sides
 *                -> check_levels
 *                -> number_levels
 *          -> sort_nodes
 *          -> create_pairs_list
 *    -> num_base
 *    -> store_variants
 *          -> extract_residue_name
 *                -> to_upper
 *                -> format_res_number
 *          -> get_other_chains
 *                -> string_truncate
 *          -> create_highres_entry
 *          -> get_variants
 *                -> string_truncate
 *                -> format_res_number
 *                -> find_highres_entry
 *                -> create_highres_entry
 *    -> draw_picture
 *          -> ss_or_ds
 *          -> find_bonded_base
 *                -> store_interaction
 *          -> find_first_pw
 *          -> sort_interactions
 *                -> interact_resnum_order
 *                      -> swap_interactions
 *                -> sep_bridged_waters
 *                -> find_bridge_bonds
 *                -> bridgeint_resnum_order
 *                -> rm_dup_int
 *                      -> move_up_one
 *          -> open_psfile
 *                -> open_file
 *          -> calc_no_of_pages
 *          -> define_object_sizes
 *                -> calc_no_of_pages
 *          -> psopen_
 *          -> base_pos
 *          -> pspage_
 *          -> text_nucleic_acids
 *                -> pslwid_
 *                -> pscomm_
 *                -> psends_
 *                      -> plot_chain_id
 *                -> pstrans_
 *                -> psctxt_
 *                -> psbasenum_
 *          -> text_interactions
 *                -> find_bonded_base
 *                -> interact_pos
 *                      -> set_pair_flag
 *                      -> count_bond_num
 *                      -> init_int_pos
 *                      -> nonphos_int_pos
 *                      -> draw_interaction_group
 *                            -> psphostxt_    
 *                            -> get_interact_det
 *                                  -> match_res_code
 *                            -> pstrans_
 *                            -> psctxt_
 *                            -> psline_
 *                            -> psctxt_
 *                      -> phos_int_pos
 *                      -> last_phos_int_pos
 *          -> text_backbone
 *                -> is_sugar_atom
 *                -> psugar_
 *                      -> get_pentagon
 *                      -> plot_pentagon
 *                -> psbasenum_
 *                -> is_phosphate_atom
 *                -> psphostxt_
 *          -> text_bridge_interactions
 *          -> psendp_
 *          -> psclos_
 *    -> draw_nucleic_acids
 *
 *************************************************************************/
/* header files =========================================================*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/* Define flag parameters ===============================================*/
#define TRUE               0
#define FALSE              1
#define LEFT               0
#define RIGHT              1
#define DONOR              0
#define ACCEPTOR           1
#define HBOND              0
#define COVALENT           1
#define NBOND              2
#define BOXES              1
#define TOP                0
#define BOTTOM             1
#define X                  0
#define Y                  1
#define ONE                0
#define TWO                1
#define NONPHOS            0 
#define PHOS               1 
#define LAST_PHOS          2
#define PD                 0 
#define WD                 1
#define PW                 2
#define BACKBONE           0
#define BASE               1
#define PHOSPHATE          2
#define FORWARDS           0
#define BACKWARDS          1
#define CANT_TELL          0
#define RISING             1
#define FALLING            2

/* v.1.1.7--> */
/* Special characters */
static char null = '\0';
/* <--v.1.1.7 */

/* No of parameters =====================================================*/
#define PRINT_OPTION       5
#define PLOT_PARAM         7
#define NUM_PARAM          3
#define COL_PARAM          3

/* Define array parameters ==============================================*/
#define COLNAME_LENGTH    15           /* String-length for colour names */
/* v.1.1.7--> */
//#define FILENAME_LEN     130
//#define LINE_LEN         130                           /* Max file width */
//#define PS_STRING_LENGTH 132
#define FILENAME_LEN     512
#define LINE_LEN        1024                           /* Max file width */
#define PS_STRING_LENGTH 512
/* <--v.1.1.7 */
#define MAXALIGN         500 /* Maximum number of base-pairs in alignment */
#define MAXARRAY          30 /* Max no. of protein chains or DNA strands */
#define MAXCHAINS        100                         /* Half of MAXARRAY */
#define MEDIARRAY         15                         /* Half of MAXARRAY */
#define MAX_COLOURS       20
#define MAX_NON_STANDARD  50
#define NUMRES            20
/* v.1.1.3--> */
//#define NUMBASES          14      /* Number of standard bases recognized */
/* v.1.1.9--> */
//#define NUMBASES          20      /* Number of standard bases recognized */
#define NUMBASES          26      /* Number of standard bases recognized */
/* <--v.1.1.9 */
/* <--v.1.1.3 */
#define OBJECT_COLOURS    26
#define MAXCON             5             /* No of allowed covalent bonds */
#define TITLE_LENGTH      80

/* Define plot parameters ===============================================*/
#define BBOXX1           50.0
#define BBOXX2          550.0
#define BBOXY1           50.0
#define BBOXY2          800.0
#define BOTTOM_MARGIN_Y  10.0  
#define LEFT_MARGIN_X    10.0
#define RIGHT_MARGIN_X   10.0
#define TOP_MARGIN_Y     10.0

#define LPAGE_POSX     (BBOXY2 - 40.0)
#define LPAGE_POSY     (BBOXX2 - 20.0)
#define PAGE_POSX      (BBOXX2 - 40.0)
#define PAGE_POSY      (BBOXY2 - 20.0)
#define PAGE_TXT         12.0

/* Define atom types */
#define DNA                 0
#define PROTEIN             1
#define WATER               2
#define HYDROGEN            3

/* Define widths and heights of picture objects =========================*/
#define BASE_TO_BASE_WIDTH    100
#define CHAIN_POS_HEIGHT       20
#define END_DASH_LENGTH        40
#define KEY_HEIGHT_LANDSCAPE  100
#define KEY_WIDTH_LANDSCAPE  (10.0 * KEY_HEIGHT_LANDSCAPE)
#define KEY_HEIGHT_PORTRAIT   120
#define KEY_WIDTH_PORTRAIT   (6.8 * KEY_HEIGHT_PORTRAIT)
#define KEY_SCALE             0.5
#define SUGAR_TO_BASE_WIDTH    60
#define P_RADIUS               12
#define P_TO_P_HEIGHT         100
#define P_TO_SUGAR_WIDTH       40
#define SUGAR_SIZE             20
#define WATER_RADIUS           12
#define WATER_TO_P_WIDTH      100
#define WATER_TO_SUGAR_WIDTH   77

/* Text sizes */
#define BASENAME_TXT      55
#define BASENUM_TXT       22
#define CHAINAME_TXT      30
#define END_TXT           40
#define KEY_TEXT_TXT      20
#define KEY_TXT           40
#define PHOS_TXT          16 
#define RESNAME_TXT       22
#define TITLE_TXT         22
#define WATER_TXT         15

/* Line widths */
#define BSLNWID          3.5  /* Bonds between bases */
#define PSLNWID          2.3  /* Phosphate lines */
#define BDLNWID          0.8  /* External bond lines */
#define EDLNWID          2.0  /* Dashed lines at either end */

/* Maths constants */
#define PI               3.141592654

/* v.1.1.10--> */
static char *T_F[3] = { "FALSE", "TRUE", "*UNSET*" };
/* <--v.1.1.10 */

/* Declare global variables =============================================*/
char  pdbfilename[FILENAME_LEN];
char  hb2filename[FILENAME_LEN];
char  nb2filename[FILENAME_LEN];
char  bondfilename[FILENAME_LEN];
char  psfilename[FILENAME_LEN];
FILE  (*ps_file);
FILE  (*bondfile);
FILE  (*parfile);


char resname[NUMRES][4] = {
  "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
  "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL" };
char oneletter[NUMRES][2] = {
  "A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
  "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V" };

/* Colour tables */
char Colour_Table_Name[MAX_COLOURS][COLNAME_LENGTH + 1];
float Colour_Table[MAX_COLOURS][3];

/* v.1.1.7--> */
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
/* <--v.1.1.7 */

/* v.1.1.8--> */
int PDBsum1;
/* <--v.1.1.8 */

typedef struct param
{
   /* print options */
   int               landscape;
   int               colour;
   int               plot_title;
   int               plot_key;
   int               write_to_nucplot_ps;

   /* plot parameters */
   int               rd_bond_file;
   int               numbering;
   int               highlight;
   int               prot_res;
   int               non_bridge_water;
   int               bridge_water;
   int               non_bond_contact;
   int               one_letter_codes;
   int               nonbond_by_residue;

   /* numerical parameters */
   int               numpairs;
   double            hb_distance;
   double            nb_distance;

   /* Colours */
   char              background_colour[COLNAME_LENGTH + 1];
   char              res_name_colour[COLNAME_LENGTH + 1];
   char              text_colour[COLNAME_LENGTH + 1];
   char              hbond_colour[COLNAME_LENGTH + 1];
   char              nbond_colour[COLNAME_LENGTH + 1];
   char              adenine_colour[COLNAME_LENGTH + 1];
   char              cytosine_colour[COLNAME_LENGTH + 1];
   char              thymine_colour[COLNAME_LENGTH + 1];
   char              guanine_colour[COLNAME_LENGTH + 1];
   char              uracil_colour[COLNAME_LENGTH + 1];
   char              sugar_colour[COLNAME_LENGTH + 1];
   char              inside_sugar_colour[COLNAME_LENGTH + 1];
   char              phosphate_colour[COLNAME_LENGTH + 1];
   char              inside_phosphate_colour[COLNAME_LENGTH + 1];
   char              backbone_colour[4][COLNAME_LENGTH + 1];
   char              bbline_colour[COLNAME_LENGTH + 1];
   char              nitrogen_colour[COLNAME_LENGTH + 1];
   char              oxygen_colour[COLNAME_LENGTH + 1];
   char              carbon_colour[COLNAME_LENGTH + 1];
   char              sulphur_colour[COLNAME_LENGTH + 1];
   char              water_colour[COLNAME_LENGTH + 1];
   char              inside_water_colour[COLNAME_LENGTH + 1];
   char              unknown_atom_colour[COLNAME_LENGTH + 1];

   /* Non-standard bases */
   int               n_nonstan;
   char              nonstandard_name[MAX_NON_STANDARD][4];

   /* Plot title */
   char              title[TITLE_LENGTH];

/* v.1.1.7--> */
   /* Additional info for highlighting residues by DisaStr variants */
   char              pdb_code[5];
   char              highlight_res[16];
   char              *disastr_dir;
   char              *pdbsum_dir;
   int               disastr;
/* <--v.1.1.7 */

} PARAM;

typedef struct atoms
{
   /* atom details */
   int               number;
   char              name[5];
   char              resname[4];
   char              chain;
   int               resnum;
   double            coord[3];
   int               type;        /* Whether DNA, PROTEIN or WATER */

   /* base pairing details */
   char              pair_name[5];
   int               pair_resnum;
   char              pair_resname[4];
   char              pair_chain;
   struct atoms      *pair_dna;
   int               dna_chain;
   int               npair_bond;

   int               numbering;   /* Number assigned to the base by 
                                     numbering system */
   /* linking pointers */
   struct atoms      *next_atom;
   struct atoms      *prev_atom;
   struct atoms      *next_res;
   struct atoms      *prev_res;

   /* hbonding information */
   int               nbonds;
   int               bridge_flag;
   struct atoms      *bondto[MAXARRAY];

   /* dna base-pair to which atom belongs */
   struct pairs      *pair;

   int               donor_acceptor[MAXARRAY];
   struct atoms      *bonded_pair[MAXARRAY]; /* the atoms to which the res
                                                is hbonded */
   struct atoms      *bonded_atom[MAXARRAY]; /* the atoms from which the
                                                res is hbonded */
   struct node       *node_ptr;
} ATOMS;

typedef struct bonds
{
   int               type;
   int               bridge_flag;
   int               entity;
   double            distance;
   struct atoms      *bond_atom[2];
   struct bonds      *prev_bond;
   struct bonds      *next_bond;
}BONDS;

typedef struct pairs
{
   /* X, Y pos of strand 1, 2 */
   double            x_pos[2], y_pos[2];         
   double            bondx_pos[2][MAXARRAY], pbondx_pos[2][MAXARRAY];
   double            bondy_pos[2][MAXARRAY], pbondy_pos[2][MAXARRAY];
   double            bdgex_pos[2][MAXARRAY][MAXARRAY];
   double            pbdgex_pos[2][MAXARRAY][MAXARRAY];
   double            bdgey_pos[2][MAXARRAY][MAXARRAY];
   double            pbdgey_pos[2][MAXARRAY][MAXARRAY];

   /* number position */
   double            x_numpos[2];                     

   int               colour[2];
   int               first_base[2];
   int               last_base[2];
   int               nbond[2];
   int               npbond[2];
   int               ndups[2][MAXARRAY], npdups[2][MAXARRAY];
   int               bridge_flag[2][MAXARRAY], pbridge_flag[2][MAXARRAY];
   int               nbridge[2][MAXARRAY], npbridge[2][MAXARRAY];
   int               type[2][MAXARRAY];  /* bond type */
   int               ptype[2][MAXARRAY]; /* bond type */
   int               paired;
   struct atoms      *bondto[2][MAXARRAY];   /* interacting atom bonded to */
   struct atoms      *pbondto[2][MAXARRAY];
   struct atoms      *bridgeto[2][MAXARRAY][MAXARRAY];
   struct atoms      *pbridgeto[2][MAXARRAY][MAXARRAY];

   struct atoms      *combine[2][MAXARRAY];   /* interacting atom bonded to */
   struct atoms      *pcombine[2][MAXARRAY];
   struct atoms      *bridge_combine[2][MAXARRAY][MAXARRAY];
   struct atoms      *pbridge_combine[2][MAXARRAY][MAXARRAY];

   struct atoms      *bondfrom[2][MAXARRAY]; /* atom in DNA bonded to */
   struct atoms      *pbondfrom[2][MAXARRAY];

   struct atoms      *strand[2];
   struct pairs      *next_pair;
   struct pairs      *prev_pair;
}PAIRS;

typedef struct sizes
{
   int               npages;
   int               total_pages;

   /* text version */
   double            txtscale;                           /* scale factor */
   double            x_centre, y_centre;                  /* page centre */
   double            x_txtstart, y_txtstart;           /* start position */
   double            y_keystart;               /* start position for key */
   double            y_title_pos;             /* position for plot title */

}SIZES;


/* Structure for linked nodes giving the connectivities and pairing of
   the DNA chains found in the structure
   ------------------------------------------------------------------- */
struct node
{
  ATOMS             *atom_ptr;
  float             level;
  struct node       *link_node_ptr[2];
  struct node       *base_pair_node_ptr;
  int               linked;
  int               paired;
  int               nhbond;
  int               side;
  int               colour;
  int               done;
  int               first_base;
  int               last_base;
  struct node       *next_node_ptr;
 };

/* v.1.1.7--> */
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
/* <--v.1.1.7 */

/* Function prototypes ==================================================*/
int
initialize_parameters(PARAM *par);

int
read_parameters(PARAM *par);

int
update_current_parameter(PARAM *par, int parameter_group, int parameter, 
                         int true_false, double number,
                         char colour[COLNAME_LENGTH + 1], float red,
                         float green, float blue);
void
prepare_colour_table(void);

void 
verify_colours(PARAM *par);

void 
verify_colour_name(char *colour, int in_colour,char *default_colour);
int
com_line(PARAM *par, int argc, char *argv[]);

void
print_usage(void);

int
read_pdb(ATOMS **first_prot, ATOMS **last_prot, int *nprot, 
         ATOMS **first_water, ATOMS **last_water, int *nwater);
int
is_hydrogen_atom(char atname[5]);

int
get_atom_info(char input_line[LINE_LEN], ATOMS *atom);

int
split(ATOMS *first_atom, ATOMS *first_prot[MAXARRAY],
      ATOMS *last_prot[MAXARRAY], ATOMS *first_DNA[MAXARRAY],
      ATOMS *last_dna[MAXARRAY], int *nprot, int *ndna, PARAM *par);
int
is_nucleic_acid(char resname[4], char *code, PARAM *par);

int
baselink(ATOMS *first_atom, ATOMS *last_atom);

int
get_bond_info(ATOMS *first_dna[MAXARRAY], ATOMS *last_dna[MAXARRAY],
              ATOMS *first_prot[MAXARRAY], ATOMS *last_prot[MAXARRAY],
              ATOMS *first_water, BONDS **first_bond,
              BONDS **first_nbond, 
              BONDS **first_connect, BONDS **last_bond,
              PARAM *par, int ndna, int nprot, int argc, char *argv[]);
int
readhbplus(ATOMS *first_dna[MAXARRAY], ATOMS *first_prot[MAXARRAY], 
           ATOMS *first_water, BONDS **first_bond, BONDS **last_bond,
           PARAM *par, int ndna, int nprot, int bond_flag);
int
getnam(char pdbfil[FILENAME_LEN],int *istart,int *iend);

void
get_filename(char *pdbfilename, char *file_name, char *root_name, 
             char *extension);

int
check_hbplus(FILE *hbfile);

int
dna_base_pair(char input_line[LINE_LEN], ATOMS *first_dna[MAXARRAY],
              int ndna, PARAM *par);
int
prot_dna_hbond(char input_line[LINE_LEN], ATOMS *first_prot[MAXARRAY], 
               ATOMS *first_dna[MAXARRAY], BONDS **first_bond,
               BONDS **last_bond, int nprot, int ndna, PARAM *par);
int
get_atom_types(char input_line[LINE_LEN], int *dna_flag, int *prot_flag,
               int *water_flag, PARAM *par);
int
water_dna_hbond(char input_line[LINE_LEN], ATOMS *first_water,
                ATOMS *first_dna[MAXARRAY],  
                BONDS **first_bond, BONDS **last_bond, int ndna, PARAM *par);
int
mark_water_hbond(ATOMS *water_atom, ATOMS *dna_atom);

int
prot_water_hbond(char input_line[LINE_LEN], ATOMS *first_prot[MAXARRAY],
                 ATOMS *first_water, BONDS **first_bond,
                 BONDS **last_bond, int nprot, PARAM *par);
int
find_hbond_res(ATOMS *first_atom[MAXARRAY], int nchains, ATOMS **hbres, 
               ATOMS **hbatom, char input_line[LINE_LEN], int flag);
int
find_hbond_water(ATOMS *first_water, ATOMS **water_res,
                 ATOMS **water_atom, 
                 char input_line[LINE_LEN], int flag);
int
hbond_warning(char input_line[LINE_LEN], ATOMS *res1, ATOMS *atom1, 
/* v.1.1.2--> */
//              ATOMS *res2, ATOMS *atom2, int *flag);
              ATOMS *res2, ATOMS *atom2, int *flag, int call_no);
/* <--v.1.1.2 */
int
check_water_bridge(ATOMS *water_atom, int *bridge_flag, PARAM *par);

int
hbond_list(BONDS **first_bond, BONDS **last_bond, ATOMS *atom1,
           ATOMS *atom2, 
           int flag1, int flag2, int bondtype, char input_line[LINE_LEN],
           int entity_flag);
int
exclude_bonds(BONDS **first_bond, BONDS **last_bond, PARAM *par);

int
read_bond_file(BONDS **first_bond, BONDS **first_nbond, 
               BONDS **first_connect, BONDS **last_bond,
               ATOMS *first_dna[MAXARRAY], ATOMS *first_prot[MAXARRAY], 
               ATOMS *first_water, PARAM *par, int ndna, int nprot);
int
get_hbond_info(BONDS **first_bond, BONDS **last_bond,
               ATOMS *first_dna[MAXARRAY], ATOMS *first_prot[MAXARRAY], 
               ATOMS *first_water, PARAM *par, int ndna, int nprot);
int
get_nbond_info(BONDS **first_bond, BONDS **first_nbond, BONDS **last_bond,
               ATOMS *first_dna[MAXARRAY], ATOMS *first_prot[MAXARRAY], 
               ATOMS *first_water, PARAM *par, int ndna, int nprot);
int
get_cbond_info(BONDS **first_bond, BONDS **first_connect,
               BONDS **last_bond,
               ATOMS *first_dna[MAXARRAY], ATOMS *first_prot[MAXARRAY], 
               ATOMS *first_water, PARAM *par, int ndna, int nprot);
int
check_bond_file(FILE *bondfile);

int
convert_bond_to_hb2(char input_line[LINE_LEN], char hb2_line[LINE_LEN]);

int
readnb2(ATOMS *first_dna[MAXARRAY], ATOMS *first_prot[MAXARRAY], 
        ATOMS *first_water, BONDS **first_bond, BONDS **last_bond,
        BONDS **first_nbond, PARAM *par, int ndna, int nprot, int argc,
        char *argv[]);
int
prot_dna_nbond(char input_line[LINE_LEN], ATOMS *first_prot[MAXARRAY], 
               ATOMS *first_dna[MAXARRAY], BONDS **first_bond,
               BONDS **last_bond,
               BONDS **first_nbond, int nprot, int ndna, PARAM *par);
int
accept_reject_nbond(BONDS *first_bond, ATOMS *prot_atom);

int
nbond_list(BONDS **first_bond, BONDS **last_bond, BONDS **first_nbond, 
           ATOMS *prot_atom, ATOMS *dna_atom, int bondtype, 
           char input_line[LINE_LEN], int entity_flag);
int
prot_dna_cbond(char input_line[LINE_LEN], ATOMS *first_prot[MAXARRAY], 
               ATOMS *first_dna[MAXARRAY], BONDS **first_bond,
               BONDS **last_bond,
               BONDS **first_connect, int nprot, int ndna, PARAM *par);
int
cbond_list(BONDS **first_bond, BONDS **last_bond, BONDS **first_connect, 
           ATOMS *prot_atom, ATOMS *dna_atom, int bondtype,
           char input_line[LINE_LEN], int entity_flag);
int
read_connect_records(ATOMS *first_prot[MAXARRAY],
                     ATOMS *last_prot[MAXARRAY], 
                     ATOMS *first_dna[MAXARRAY],
                     ATOMS *last_dna[MAXARRAY],
                     BONDS **first_bond, BONDS **last_bond,
                     BONDS **first_connect,
                     int nprot, int ndna);
int
get_connect_numbers(char input_line[LINE_LEN],
                    ATOMS *first_prot[MAXARRAY],
                    ATOMS *last_prot[MAXARRAY], int nprot,
                    ATOMS *first_dna[MAXARRAY],
                    ATOMS *last_dna[MAXARRAY], int ndna,
                    int numbers[MAXCON][2], int store_prot[MAXCON],
                    int store_dna[MAXCON],
                    int *inprot, int *indna, int *ncon);
int
connect_list(ATOMS *first_prot[MAXARRAY], ATOMS *first_dna[MAXARRAY], 
             BONDS **first_bond, BONDS **last_bond,
             BONDS **first_connect, int numbers[MAXCON][2],
             int store_prot[MAXCON], int store_dna[MAXCON], 
             int inprot, int indna, int ncon);
int
check_connect_record(BONDS **first_bond, BONDS **last_bond,
                     double *distance);
int
mark_first_connect(BONDS **first_bond, BONDS **first_connect);

int
write_bonds(BONDS *first_bond, BONDS *first_nbond, BONDS *first_connect);

void
paircatalog(ATOMS *first_dna[MAXARRAY], ATOMS *last_dna[MAXARRAY],
            ATOMS *second_strand[MEDIARRAY], 
            PAIRS *beg_pair[MEDIARRAY], PAIRS *end_pair[MEDIARRAY],
            int ndna, int *nchain, ATOMS *first_atom);
void
create_nodes(ATOMS *first_dna[MAXARRAY], int ndna,
             struct node **first_node_ptr, struct node **last_node_ptr);
void
link_nodes(struct node *first_node_ptr, ATOMS *first_atom);

void
assign_sides(ATOMS *first_dna[MAXARRAY], int ndna,
             struct node *first_node_ptr);
void
check_levels(struct node *first_node_ptr, int dna_chain, float max_level,
             float *start_level, int *direction);
void
number_levels(struct node *first_node_ptr, int dna_chain, float start_level,
              int direction, float *max_level);
void
sort_nodes(struct node **first_node_ptr);

void
create_pairs_list(struct node *first_node_ptr, PAIRS *beg_pair[MEDIARRAY],
                  PAIRS *end_pair[MEDIARRAY], int *nchain);
int
num_base(int num_system, ATOMS *first_dna[MAXARRAY], int ndna);

/* Picture drawing functions */ 
int
draw_picture(PAIRS *beg_pair[MEDIARRAY], PAIRS *end_pair[MEDIARRAY],
             int nchain, BONDS *first_bond, BONDS *first_nbond,
             BONDS *first_connect, PARAM *par);
int
ss_or_ds(PAIRS *first_pair, int ds_flag[2]);

int
find_bonded_base(PAIRS *first_pair, PAIRS *last_pair, 
                 BONDS *first_bond, PARAM *par, int *max_nbond,
                 int *max_npbond);
void
store_interaction(PAIRS *pair,PAIRS *next_pair,PAIRS *prev_pair,
                  BONDS *bond,ATOMS *bond_atom[2],int i,int j,int k);
int
find_first_pw(BONDS *first_bond, BONDS **first_pw_bond);

int
sort_interactions(PAIRS *first_pair, PAIRS *last_pair,
                  BONDS *first_bond, 
                  PARAM *par, int *max_nbridge, int *max_npbridge);
int
interact_resnum_order(PAIRS *pair, int i, int *nprot_int,
                      int *nwater_int, int phos_flag);
void
swap_interactions(PAIRS *pair, int i, int j, int k, int phos_flag);

int
sep_bridged_waters(PAIRS *pair, int i, int nprot_int, int nwater_int,
                   int phos_flag);
int
find_bridge_bonds(PAIRS *pair, BONDS *first_bond, int i, int nprot_int,
                  int nwater_int, int *max_nbridge, int *max_npbridge,
                  int phos_flag);
int
bridgeint_resnum_order(PAIRS *pair, int i, int nprot_int,
                       int nwater_int, int phos_flag);
int
rm_dup_int(PAIRS *pair, int i, int phos_flag, int nonbond_by_residue);

void
move_up_one(PAIRS *pair, int i, int j, int nbonds, int phos_flag);

int
define_object_sizes(PAIRS *beg_pair[MEDIARRAY], int ichain, SIZES *size,
                    PARAM *par,
                    int max_nbond, int max_npbond, 
                    int max_nbridge, int max_npbridge, int ds_flag[2]);
int
calc_no_of_pages(PAIRS  *pair, int numpairs, int *npair);

int
open_psfile(PAIRS *beg_pair[MEDIARRAY], int i, PARAM *par);

void
psopen_(SIZES *size, PARAM *par);

int
base_pos(PAIRS *beg_pair, SIZES *size, PARAM *par, int ds_flag[2]);

PAIRS*
text_nucleic_acids(PAIRS *beg_pair[MEDIARRAY],
                    PAIRS *end_pair[MEDIARRAY], 
                    PAIRS *first_pair, int pairno, SIZES *size,
                    PARAM *par, int ds_flag[2]);
int
text_backbone(PAIRS *first_pair, PAIRS *last_pair, SIZES *size,
              PARAM *par);
int
is_sugar_atom(char atname[5]);

int
is_phosphate_atom(char atname[5]);

int
text_interactions(PAIRS *first_pair, PAIRS *last_pair, 
                  PAIRS *first_page_pair, PAIRS *last_page_pair, 
                  SIZES *size, PARAM *par,
                  int max_nbond, int max_npbond);
int
set_pair_flag(PAIRS *pair, PAIRS *first_pair, PAIRS *last_pair, 
              int *first_flag, int *last_flag);
int
count_bond_num(PAIRS *pair, int i, int nbond[3], int nslots[3]);

int
init_int_pos(PAIRS *pair, SIZES *size, PARAM *par, int i,
             int nslots[3], double y_spacing[3], double text_size[3],
             double bondx_pos[3], double bondy_pos[3]);
int
nonphos_int_pos(PAIRS *pair, SIZES *size, PARAM *par, int i,
                double bondx_pos, double bondy_pos, double y_spacing,
                double text_size);
int
phos_int_pos(PAIRS *pair, SIZES *size, PARAM *par, int first_flag, int i,
             double bondx_pos, double bondy_pos, double y_spacing,
             double text_size);
int
last_phos_int_pos(PAIRS *pair, SIZES *size, PARAM *par, int last_flag, int i, 
                  double bondx_pos, double bondy_pos, double y_spacing,
                  double text_size);
double
calc_bridge_dist(double text_size, PARAM *par);

int
draw_interaction_group(PAIRS *first_pair, PAIRS *last_pair, PAIRS *pair,
                       PARAM *par, SIZES *size, int i, int last_flag,
                       int bond_point, double y_spacing, double text_size);
int
get_interact_det(PAIRS *first_pair, PAIRS *last_pair, PAIRS *pair,
                 PARAM *par, ATOMS *interaction, ATOMS *int_combine,
                 char string[LINE_LEN], int i, int bond_point,
                 char atom_name[5], char atom_colour[COLNAME_LENGTH + 1],
                 char atom_name_combine[11],
                 char atom_colour_combine[COLNAME_LENGTH + 1],
                 int *water_flag, int bond_type);
void
get_atom_name(PARAM *par, ATOMS *atom_ptr, char atom_name[5],
              char atom_colour[COLNAME_LENGTH + 1]);
int
check_for_double(PAIRS *first_pair, PAIRS *last_pair, PAIRS *current_pair,
                 ATOMS *current_interaction, int current_strand,
                 int current_bond_group);
int
draw_bond_line(PAIRS *pair, PARAM *par, SIZES *size, int i, int ibond,
               int bond_point, double y_spacing);
int
match_res_code(PARAM *par, char *three_letter, char *string);

int
text_draw_bdge_int(PAIRS *first_pair, PAIRS *last_pair, PAIRS *pair,
                   PAIRS *next_pair, SIZES *size,
                   PARAM *par, int i, int flag, int last_flag,
                   double text_size[3]);
int
bridge_bond_lines(PAIRS *pair, PARAM *par, SIZES *size, 
                  int i, int phos_flag, double text_size);
void
set_dash_type(int line_type);

/* v.1.1.7--> */
struct highres
*find_highres_entry(struct highres *first_highres_ptr,
                    char *res_name,char *res_num,
                    char chain);

struct highres
*create_highres_entry(struct highres **fst_highres_ptr,
                      struct highres **lst_highres_ptr,
                      int type,int ref,char *res_name,
                      char *res_num,char chain);

void
format_res_number(char *tmpstr,char *res_num);

void
extract_residue_name(char *highlight_res,char *res_name,
                     char *res_num,char *chn);

int
get_variants(struct highres **fst_highres_ptr,
             struct highres **lst_highres_ptr,int nhigh,
             char *chains,int nchains,PARAM *par);

int
store_variants(struct highres **fst_highres_ptr,PARAM *par);

void
store_string(char **string_ptr,char *string);
/* <--v.1.1.7 */

void
psbasecol_(ATOMS *atom, PARAM *par, char *code);

void
psbasenum_(PAIRS *pair, SIZES *size, PARAM *par, int i, double x,
           double y);
void 
psccol_(char *colour);

void 
pscirc_(double x, double y, double radius);

void
pscolb_(char *text_string);

void
pscomm_(char *text_string);

void
psltxt_(double x,double y,int size, char *text_string);
void
psrtxt_(double x,double y,int size, char *text_string);

void
psrltx_(double x,double y,int size, char *text_string);

void
psrotx_(double x,double y,int size, char *text_string);

void
psctxt_(double x,double y,int size,
        char *text_string);
void
psdash_(int ONOFF1, int ONOFF2);

void
psendp_(void);

void
psends_(PAIRS *edge_pair[MEDIARRAY], PAIRS *pair, SIZES *size, 
        PARAM *par, int pairno, int which_end, int ds_flag[2],
        char chain_top[2], char chain_bottom[2]);
void
plot_key(PARAM *par, SIZES *size, double x_start, double y_start,
         double key_height, double key_width);
void
plot_chain_id(char chain, double x, double y, double text_size,
              PARAM *par, SIZES *size);
void
psline_(double x1,double y1,double x2,double y2);

void
pslwid_(double width);

void
pspage_(int page_no, SIZES *size, PARAM *par);

void
psphostxt_(PAIRS *pair, SIZES *size, PARAM *par, double xstart,
           double x_line[2], double y_line[2], 
           double prev_xline[2], double prev_yline[2], 
           int flag, int end_flag, int first_base, int last_base,
           int top_of_page);
void
psrest_(void);

void
psrot_(float new_origin_x, float new_origin_y);

void
pssave_(void);

void
pselip_(float x,float y,float width,float ratio,float angle);

void
psetgray_(double grayscale);

void
pshade_(char *colour);

void
pstrans_(double x, double y, double xo, double yo, 
         double *x_new, double *y_new, PARAM *par);
void
psubox_(double x1, double y1, double x2, double y2,
        double x3, double y3, double x4, double y4);
void
psugar_(PAIRS *pair, SIZES *size, PARAM *par, int i, double xstart,
        double ystart, double x_line[2], double y_line[2]);
void
get_pentagon(double x[5], double y[5], double xstart, double ystart,
             double sugar_size, int i);
void
plot_pentagon(double x[5], double y[5], PARAM *par, double line_width);

void
pswatertxt_(SIZES *size, PARAM *par, double x1, double y1,
            char int_det[LINE_LEN], int i, double text_size,
            int ibond);
void
psclos_(void);

FILE*
open_file(char filename[FILENAME_LEN], char *open_type);




/**************************************************************************
 *
 *
 * Main program
 *
 *
 *************************************************************************/
int main(int argc, char *argv[])
{
/* Declare variables ====================================================*/
  int   i;
  int   num_system;                      /* Numbering system identifier */
  int   natoms, natoms_p[MAXARRAY], natoms_d[MAXARRAY]; /* No. of atoms */
  int   nprot,  ndna;          /* No. of protein chains and dna strands */
  int   nchain;                              /* No. of pairs of strands */
  int   nwaters;
  PARAM parameters;
  PARAM *par;
   
  ATOMS *first_atom, *last_atom;
  ATOMS *first_water, *last_water;
  ATOMS *first_prot[MAXARRAY], *last_prot[MAXARRAY];
  ATOMS *first_dna[MAXARRAY], *last_dna[MAXARRAY];

  ATOMS *beg_chain[MEDIARRAY], *end_chain[MEDIARRAY];
  ATOMS *second_strand[MEDIARRAY];
  BONDS *first_bond, *last_bond, *first_nbond, *first_connect;
  PAIRS *pair[MEDIARRAY], *beg_pair[MEDIARRAY], *end_pair[MEDIARRAY];
  SIZES size;
  SIZES *ptsize;
   
/* Initialize variables =================================================*/
  num_system = 0;
  natoms     = nprot     = ndna          = nwaters     = nchain = 0;
  first_atom = last_atom = first_water   = last_water  = NULL;
  first_bond = last_bond = first_connect = first_nbond = NULL;
  pdbfilename[0] = hb2filename[0] = nb2filename[0] = bondfilename[0] = '\0';

  for (i = 0; i < MAXARRAY; i++)
    {
      natoms_p[i]   = 0;
      first_prot[i] = NULL;
      natoms_d[i]   = 0;
      first_dna[i]  = NULL;
    }

  for (i = 0; i < MEDIARRAY; i++)
    {
      beg_chain[i] = NULL;
      end_chain[i] = NULL;
      second_strand[i] = NULL;
      beg_pair[i] = NULL;
      end_pair[i] = NULL;
      pair[i]     = NULL;
    }
  ptsize = &size;
  par    = &parameters;
/* v.1.1.8--> */
  PDBsum1 = FALSE;
/* <--v.1.1.8 */
   
/* Initialize parameters ================================================*/
  initialize_parameters(par);

  /* Check to see if have parameter file */
  if ((parfile = fopen("nucplot.par", "r")) == NULL)
    {
      printf("\n");
      printf("Note: No parameter file supplied - default settings");
      printf(" will be used\n");
      printf("\n");
    }

  /* Read in the parameters */
  else
    {
      printf("\n");
      printf("Reading parameter file ...\n");
      read_parameters(par);
    }

  /* Prepare the colour tables */
  prepare_colour_table();

  /* Check that the entered colours are valid */
  verify_colours(par);

  /* Get the command-line parameters, including input file name */
  num_system = com_line(par, argc, argv);

/* Check usage ==========================================================*/
  if (num_system == 0)
    printf("default numbering...\n");

/* Read PDB file ========================================================*/
   read_pdb(&first_atom, &last_atom, &natoms, &first_water, &last_water,
            &nwaters);

/* Segment the stored information into protein and DNA parts ============*/
   split(first_atom, first_prot, last_prot, first_dna, last_dna, 
         &nprot, &ndna, par); 
   printf("Number of DNA strands identified:     %4d\n",ndna);
   printf("Number of protein chains identified:  %4d\n",nprot);

/* Check that have at least one DNA chain */
  if (ndna == 0)
    {
      printf("\n");
      printf("*** No nucleic acid chains found in the PDB file.\n");
      printf("*** Nothing to plot.\n");
      printf("\n");
      exit(1);
    }

/* Link the first atom for each base ====================================*/
   baselink(first_atom, last_atom);

/* Bonding information ==================================================*/
   get_bond_info(first_dna, last_dna, first_prot, last_prot, first_water, 
                 &first_bond, &first_nbond, &first_connect, &last_bond, 
                 par, ndna, nprot, argc, argv);

/* Store strand pairing info in an array ================================*/
   paircatalog(first_dna, last_dna, second_strand, beg_pair, end_pair,
               ndna, &nchain, first_atom);

/* Assign the numbering system to each base =============================*/
   num_base(num_system, first_dna, ndna);


/* @@@ */
   
/*=========================================================================
 * Draw pictures 
 *=======================================================================*/
   draw_picture(beg_pair, end_pair, nchain, first_bond, first_nbond,
                first_connect, par);

/* End of program =======================================================*/
   return(0);
}


   
/**************************************************************************
 *
 * initialize_parameters - Function to initialize the drawing parameters
 *         
 *************************************************************************/
int
initialize_parameters(PARAM *par)
{
  int ibase, icolour;

/* Print options ========================================================*/
  par->landscape   = FALSE;
  par->colour      = TRUE;
  par->plot_title  = TRUE;
  par->plot_key    = TRUE;
  par->write_to_nucplot_ps = TRUE;

/* Plot parameters ======================================================*/
  par->rd_bond_file     = FALSE;
  par->numbering        = TRUE;
  par->highlight        = FALSE;
  par->prot_res         = TRUE;
  par->non_bridge_water = TRUE;
  par->bridge_water     = TRUE;
  par->non_bond_contact = TRUE;
  par->one_letter_codes = TRUE;
  par->nonbond_by_residue = TRUE;

/* Numerical parameters =================================================*/
  par->numpairs    = 10;
  par->hb_distance = 3.20;
  par->nb_distance = 3.50;

/* Colours ==============================================================*/
  par->colour      = FALSE;
  verify_colours(par);
  par->colour      = TRUE;

  /* Colour definitions */
  for (icolour = 0; icolour < MAX_COLOURS; icolour++)
    {
      strcpy(Colour_Table_Name[icolour],"Black          ");
      Colour_Table[icolour][0] = 0.0000;
      Colour_Table[icolour][1] = 0.0000;
      Colour_Table[icolour][2] = 0.0000;
    }

  /* Non-standard bases */
  par->n_nonstan = 0;
  for (ibase = 0; ibase < MAX_NON_STANDARD; ibase++)
    par->nonstandard_name[ibase][0] = '\0';
  
/* Return to main =======================================================*/
  return(0);
}



/**************************************************************************
 *
 * read_parameters - Function to read nucplot.par
 *         
 *************************************************************************/
int
read_parameters(PARAM *par)
{
/* Declare variables ====================================================*/
  int     nline;
   
  int     true_false;
  double  number;
  char    colour_name[COLNAME_LENGTH + 1];

  int     parameter;
  int     parameter_group, last_group;
  int     nparams, no_of_parameters[] = {
    PRINT_OPTION, PLOT_PARAM, NUM_PARAM, OBJECT_COLOURS, MAX_COLOURS,
    MAX_NON_STANDARD + 1 };
  char    parameter_type, parameter_types[] = { 'P', 'Y', 'N', 'C', 'D', 'R' };
  char    group_names[6][20] = {
    "PRINT OPTIONS", "PLOT PARAMETERS", "NUMERICAL PARAMETERS",
    "COLOUR PARAMETERS", "COLOUR DEFINITIONS", "NON-STANDARD BASES" };
  char    input_line[LINE_LEN];
  
  float blue, green, red;

/* Initialize variables =================================================*/
  nline = 0;
  parameter_group = last_group = -1;
/* v.1.1.5--> */
  parameter = -1;
  nparams = 0;
/* <--v.1.1.5 */

/*=========================================================================
 * Read file
 *=======================================================================*/
  while (fgets(input_line, LINE_LEN, parfile) != NULL)
    {
      
/* Check that it is a parameter file ====================================*/
      if ((strncmp(input_line, "NUCPLOT v.1.0  -  Parameter file", 32))
          && (nline == 0))
        {
          printf("*** ERROR - Incorrect parameter file.\n");
          printf("            Must contain \"NUCPLOT v.1.0  -  Parameter file\" in the first line.\n");
          exit(-1);
        }
      
/* Get the parameter group ==============================================*/
      if (!strncmp(input_line, group_names[0], 13))
        parameter_group = 0;
      else if (!strncmp(input_line, group_names[1], 15))
        parameter_group = 1;
      else if (!strncmp(input_line, group_names[2], 20))
        parameter_group = 2;
      else if (!strncmp(input_line, group_names[3], 16))
        parameter_group = 3;
      else if (!strncmp(input_line, group_names[4], 18))
        parameter_group = 4;
      else if (!strncmp(input_line, group_names[5], 18))
        parameter_group = 5;

/* Initialize line_no if new parameter group ============================*/
      if (parameter_group != last_group)
        {
          parameter      = -2;
          last_group     = parameter_group;
          nparams        = no_of_parameters[parameter_group];
          parameter_type = parameter_types[parameter_group];
        }

/* Extract the value of the parameter ===================================*/
      if (parameter > -1 && parameter < nparams)
        {
          /* Yes/No response */
          if (parameter_type == 'Y')
            {
              if (input_line[0] == 'Y' || input_line[0] == 'y')
                true_false = TRUE;
              else if (input_line[0] == 'N' || input_line[0] == 'n')
                true_false = FALSE;
              else
                {
                  true_false = FALSE;
                  printf("*** Warning - Error in parameter file \"%s\" ",
                         group_names[parameter_group]);
                  printf("at line number %2d.\n", parameter+1);
                  printf("              Set to default\n");
                }
            }

          /* If this is one of the print options, interpret it */
          else if (parameter_type == 'P')
            {
              if (parameter != 1 && (input_line[0] == 'Y' ||
                                     input_line[0] == 'y'))
                true_false = TRUE;
              else if (parameter == 1 && (input_line[0] == 'L' ||
                                          input_line[0] == 'l'))
                true_false = TRUE;
              else
                true_false = FALSE;
            }

         /* If a number parameter */
         else if (parameter_type  == 'N')
            sscanf(input_line,"%lf ",&number);
        
         /* If a colour parameter */
         else if (parameter_type  == 'C' && par->colour == TRUE)
         {
            strncpy(colour_name,input_line,COLNAME_LENGTH);
            colour_name[COLNAME_LENGTH] = '\0';
         }

         /* If name of non-standard base */
         else if (parameter_type  == 'R')
         {
           /* Check that maximum number of nonstandard bases not
              exceeded */
           if (par->n_nonstan < MAX_NON_STANDARD)
             {
               strncpy(par->nonstandard_name[par->n_nonstan],input_line,3);
               par->nonstandard_name[par->n_nonstan][3] = '\0';
               par->n_nonstan++;
             }

           /* Otherwise, show a warning message */
           else
             {
               printf("*** Warning. Maximum number of non-standard bases");
               printf(" exceeded (%d)\n",MAX_NON_STANDARD);
             }
         }
        
         /* If this is a colour definition line, read in the RGB
            values and the colour name */
         else if (parameter_type == 'D')
           {
             /* Read in the RGB values */
             sscanf(input_line,"%f %f %f ",&red,&green,&blue);
              
             /* Read in the colour name */
             strncpy(colour_name,input_line+22,COLNAME_LENGTH);
             colour_name[COLNAME_LENGTH] = '\0';
           }

/* Update the current parameter =========================================*/
          if (parameter_type  != 'R')
            update_current_parameter(par, parameter_group, parameter, 
                                     true_false, number, colour_name,
                                     red, green, blue);
        }

/* Increment the parameter count ========================================*/
      nline++;
      parameter++;
    }

/* Return to main =======================================================*/
  return(0);
}



/**************************************************************************
 *
 * update_current_parameter - Function to record the details to the 
 *                            current parameter 
 *         
 *************************************************************************/
int
update_current_parameter(PARAM *par, int parameter_group, int parameter, 
                         int true_false, double number,
                         char colour[COLNAME_LENGTH + 1], float red,
                         float green, float blue)
{
/* Declare variables ====================================================*/

/* Initialize variables =================================================*/

/* Process according to parameter group =================================*/
  switch (parameter_group)
    {
      
/* Print option =========================================================*/
    case 0:
      switch (parameter)
        {
        case 0:  par->colour              = true_false; break;
        case 1:  par->landscape           = true_false; break;
        case 2:  par->plot_title          = true_false; break;
        case 3:  par->plot_key            = true_false; break;
        case 4:  par->write_to_nucplot_ps = true_false; break;
        default: printf("*** ERROR - Invalid print option in nucplot.par\n");
        }
      break;
      
/* Plot parameters ======================================================*/
    case 1:
      switch (parameter)
        {
        case 0: par->rd_bond_file       = true_false; break;
        case 1: par->prot_res           = true_false; break;
        case 2: par->non_bridge_water   = true_false; break;
        case 3: par->bridge_water       = true_false; break;
        case 4: par->non_bond_contact   = true_false; break;
        case 5: par->one_letter_codes   = true_false; break;
        case 6: par->nonbond_by_residue = true_false; break;
        default: printf("*** ERROR - Invalid plot parameter in nucplot.par\n");
        } 
      break;
      
/* Numerical parameters =================================================*/
    case 2:
      switch (parameter)
        {
        case 0: par->numpairs = (int) number; break;
        case 1: par->hb_distance = number; break;
        case 2: par->nb_distance = number; break;
        default: printf("*** ERROR - Invalid numerical parameter");
          printf(" in nucplot.par\n");
        }
      break;

/* Colour parameters ====================================================*/
    case 3:
      switch (parameter)
        {
        case  0: strcpy(par->background_colour,       colour); break;
        case  1: strcpy(par->res_name_colour,         colour); break;
        case  2: strcpy(par->text_colour,             colour); break;
        case  3: strcpy(par->hbond_colour,            colour); break;
        case  4: strcpy(par->nbond_colour,            colour); break;
        case  5: strcpy(par->adenine_colour,          colour); break;
        case  6: strcpy(par->cytosine_colour,         colour); break;
        case  7: strcpy(par->thymine_colour,          colour); break;
        case  8: strcpy(par->guanine_colour,          colour); break;
        case  9: strcpy(par->uracil_colour,           colour); break;
        case 10: strcpy(par->sugar_colour,            colour); break;
        case 11: strcpy(par->inside_sugar_colour,     colour); break;
        case 12: strcpy(par->phosphate_colour,        colour); break;
        case 13: strcpy(par->inside_phosphate_colour, colour); break;
        case 14: strcpy(par->backbone_colour[0],      colour); break;
        case 15: strcpy(par->backbone_colour[1],      colour); break;
        case 16: strcpy(par->backbone_colour[2],      colour); break;
        case 17: strcpy(par->backbone_colour[3],      colour); break;
          case 18: strcpy(par->bbline_colour,           colour); break;
          case 19: strcpy(par->nitrogen_colour,         colour); break;
        case 20: strcpy(par->oxygen_colour,           colour); break;
        case 21: strcpy(par->carbon_colour,           colour); break;
        case 22: strcpy(par->sulphur_colour,          colour); break;
        case 23: strcpy(par->water_colour,            colour); break;
        case 24: strcpy(par->inside_water_colour,     colour); break;
        case 25: strcpy(par->unknown_atom_colour,     colour); break;
        default: printf("*** ERROR - Invalid colour parameter ");
          printf(" in nucplot.par\n");
        }
      break;

      /* Colour definitions */
    case 4 :
      /* Update the appropriate parameter */
      Colour_Table[parameter][0] = red;
      Colour_Table[parameter][1] = green;
      Colour_Table[parameter][2] = blue;
      strcpy(Colour_Table_Name[parameter],colour);
      break;
    }
   
/* Return to read_parameters ============================================*/
  return(0);
}

/**************************************************************************
 *
 * prepare_colour_table - Modify colour-names in the table such that any
 *                        inter-word spaces are replaced by underscores
 *         
 *************************************************************************/
void
prepare_colour_table(void)
{
  int icolour, iletter, last_letter;

  /* Loop through all the colour entries */
  for (icolour = 0; icolour < MAX_COLOURS; icolour++)
    {
      /* Initialise position of last non-blank letter in the name */
      last_letter = 0;

      /* Process the name to replace any inter-word spaces by underscores */
      for (iletter = 0; iletter < COLNAME_LENGTH; iletter++)
        {
          /* Replace blanks by underscores */
          if (Colour_Table_Name[icolour][iletter] == ' ')
            Colour_Table_Name[icolour][iletter] = '_';
          else
            last_letter = iletter;
        }

      /* Blank out any trailing spaces replace by underscores */
      for (iletter = last_letter + 1; iletter < COLNAME_LENGTH; iletter++)
        Colour_Table_Name[icolour][iletter] = ' ';
    }
}


/**************************************************************************
 *
 * verify_colours - Verify that the colour names in the parameter file
 *                  are valid
 *         
 *************************************************************************/
void verify_colours(PARAM *par)
{
  int i;

  /* Background colour */
  verify_colour_name(par->background_colour, par->colour, "White");

  /* Text colours */
  verify_colour_name(par->res_name_colour, par->colour, "Black");
  verify_colour_name(par->text_colour, par->colour, "Black");
  
  /* H-bond colour */
  verify_colour_name(par->hbond_colour, par->colour, "Black");
  
  /* Non-bonded contact colour */
  verify_colour_name(par->nbond_colour, par->colour, "Black");

  /* Nucleic acid base colours */
  verify_colour_name(par->adenine_colour, par->colour, "Black");
  verify_colour_name(par->cytosine_colour, par->colour, "Black");
  verify_colour_name(par->thymine_colour, par->colour, "Black");
  verify_colour_name(par->guanine_colour, par->colour, "Black");
  verify_colour_name(par->uracil_colour, par->colour, "Black");

  /* Sugar colours */
  verify_colour_name(par->sugar_colour, par->colour, "Black");
  verify_colour_name(par->inside_sugar_colour, par->colour, "White");

  /* Phosphate colours */
  verify_colour_name(par->phosphate_colour, par->colour, "Black");
  verify_colour_name(par->inside_phosphate_colour, par->colour, "White");

  /* Backbone colour */
  for (i = 0; i < 4; i++)
    verify_colour_name(par->backbone_colour[i], par->colour, "Black");

  /* Base-to-base and continuation lines */
  verify_colour_name(par->bbline_colour, par->colour, "Black");

  /* Atom colours */
  verify_colour_name(par->nitrogen_colour, par->colour, "Black");
  verify_colour_name(par->oxygen_colour, par->colour, "Black");
  verify_colour_name(par->carbon_colour, par->colour, "Black");
  verify_colour_name(par->sulphur_colour, par->colour, "Black");

  /* Water colours */
  verify_colour_name(par->water_colour, par->colour, "Black");
  verify_colour_name(par->inside_water_colour, par->colour, "White");

  /* Other atoms */
  verify_colour_name(par->unknown_atom_colour, par->colour, "Black");
}



/**************************************************************************
 *
 * verify_colour_name - Verify that given colour name is in the Colour Table
 *         
 *************************************************************************/
void verify_colour_name(char *colour, int in_colour,
                        char *default_colour)
{
  char new_name[COLNAME_LENGTH + 1];
  int icolour, iletter, last_letter;
  int found;

  /* Initialise variables */
  found = FALSE;

  /* If plotting in black-and-white, transfer the default colour across */
  if (in_colour == FALSE)
    {
      strcpy(colour,default_colour);
      return;
    }

  /* Save the colour name */
  strcpy(new_name,colour);

  /* Initialise position of last non-blank letter in the name */
  last_letter = 0;

  /* Process the name to replace any inter-word spaces by underscores */
  for (iletter = 0; iletter < COLNAME_LENGTH; iletter++)
    {
      /* Replace blanks by underscores */
      if (new_name[iletter] == ' ')
        new_name[iletter] = '_';
      else
        last_letter = iletter;
    }

  /* Blank out any trailing spaces replace by underscores */
  for (iletter = last_letter + 1; iletter < COLNAME_LENGTH; iletter++)
    new_name[iletter] = ' ';

  /* Loop through the colour definitions until match is found */
  for (icolour = 0; icolour < MAX_COLOURS && found == FALSE; icolour++)
    {
      /* If colour matches then set flag to TRUE */
      if (!strncmp(Colour_Table_Name[icolour],new_name,COLNAME_LENGTH))
        found = TRUE;
    }

  /* If colour not matched, then set to default colour */
  if (found == FALSE)
    {
      if (strcmp(new_name,"White") && strcmp(new_name,"Black"))
        {
          printf("*** Invalid colour name [%s]. Set to default colour [%s]\n",
                 colour,default_colour);
          strcpy(new_name,default_colour);
        }
    }

  /* Save the new colour name */
  strcpy(colour,new_name);
}


/**************************************************************************
 *
 * com_line - Function to check the command line arguments and extract 
 *            information from them
 *         
 *************************************************************************/
int com_line(PARAM *par, int argc, char *argv[])
{
/* Declare variables ====================================================*/
  char root_name[FILENAME_LEN];
/* v.1.1.8--> */
  char tmp[FILENAME_LEN];
/* <--v.1.1.8 */
  int  i;
  int  recog;                        /* See if flag has been recognized */
  int  flag;      /* See if numbering system has already been specified */
  int  gotpdb;                      /* See if have got the PDB filename */
  int  num_system;                       /* numbering system identifier */


/* Initialize variables =================================================*/
  i          = 0;
  num_system = 0;
  recog      = FALSE; 
  flag       = FALSE;
  gotpdb     = FALSE;
  par->title[0] = '\0';
/* v.1.1.7--> */
  par->pdb_code[0] = '\0';
  par->highlight_res[0] = '\0';
  par->disastr_dir = &null;
  par->pdbsum_dir = &null;
  par->disastr = FALSE;
/* <--v.1.1.7 */

/* Are there any command line arguments? - if no print usage ============*/
/* v.1.1.7--> */
//  if (argc < 2 || argc > 2)
  if (argc < 2)
/* <--v.1.1.7 */
    print_usage();

/* v.1.1.7--> */
  /* Get the name of the input PDB file */
  strcpy(pdbfilename, argv[1]);
/* <--v.1.1.7 */
      
/*=========================================================================   
 * Loop through the remaining command line arguments
 *=======================================================================*/
/* v.1.1.7--> */
  for (i = 2; i < argc; i++)
    {
      /* Store the input PDB file name */
//      strcpy(pdbfilename, argv[i]);

      /* Check for the PDB code */
      if (!strcmp(argv[i],"-pdb"))
        {
          /* Get the PDB code from the next parameter */
          if (i < argc - 1)
            {
              strcpy(par->pdb_code,argv[i + 1]);
              i++;
            }
        }
      
      /* Get the DisaStr residue to be highlighted */
      else if (!strcmp(argv[i],"-res"))
        {
          /* Get the residue from the next parameter */
          if (i < argc - 1)
            {
              strcpy(par->highlight_res,argv[i + 1]);
              i++;
            }
        }
      
      /* Check for the PDBsum directory path */
      else if (!strcmp(argv[i],"-pdbsum"))
        {
          /* Get the PDBsum directory from the next parameter */
          if (i < argc - 1)
            {
              store_string(&(par->pdbsum_dir),argv[i + 1]);
              i++;
            }
        }

/* v.1.1.8--> */
      /* Check for the PDBsum1 flag */
      else if (!strcmp(argv[i],"-pdbsum1"))
        PDBsum1 = TRUE;
/* <--v.1.1.8 */

      /* Check for the DisaStr directory path */
      else if (!strcmp(argv[i],"-disastr"))
        {
          /* Get the DisaStr directory from the next parameter */
          if (i < argc - 1)
            {
              store_string(&(par->disastr_dir),argv[i + 1]);
              i++;
            }
        }
      
    }

  /* See if we have all the required parameters for a DisaStr run */
  if (par->pdb_code[0] != '\0' && par->highlight_res[0] != '\0' &&
      par->disastr_dir != &null && par->pdbsum_dir != &null)
    par->disastr = TRUE;
/* <--v.1.1.7 */

/* Get the other filenames and print them out ===========================*/
  printf("Input file names:-\n");
  printf("   PDB file              %s\n",pdbfilename);
  get_filename(pdbfilename,hb2filename,root_name,".hb2");
/* v.1.1.8--> */
  if (PDBsum1 == TRUE)
    {
      sprintf(tmp,"../hbplus/%s",hb2filename);
      strcpy(hb2filename,tmp);
    }
/* <--v.1.1.8 */
  printf("   H-bonds file          %s\n",hb2filename);
  get_filename(pdbfilename,nb2filename,root_name,".nb2");
/* v.1.1.8--> */
  if (PDBsum1 == TRUE)
    {
      sprintf(tmp,"../hbplus/%s",nb2filename);
      strcpy(nb2filename,tmp);
    }
/* <--v.1.1.8 */
  printf("   Non-bonded contacts   %s\n",nb2filename);
  get_filename(pdbfilename,bondfilename,root_name,".bond");
  printf("   Bonds file            %s\n",bondfilename);
  printf("\n");

/* Initialise the Postscript file name as the root name */
  strcpy(psfilename,root_name);

/* Save the root-name as the title for the plot */
  if (par->title[0] == '\0')
    strcpy(par->title,root_name);

/* Return to main =======================================================*/
  return(num_system);
}



/**************************************************************************
 *
 * print_usage - Function to print the command line usage
 *
 *************************************************************************/
void print_usage(void)
{
   printf("\nUSAGE:\n");
   printf("      nucplot PDB_file\n");
   printf("\n");
   printf("   where:-   PDB_file = name of input PDB file\n");
   printf("\n");
/*   printf("numbering SYSTEMS:\n");
   printf("      See manual for a full explanation of numbering systems.\n");
   printf("      Default numbering is the one found in the PDB file.\n");
   printf("      -a                   antiparallel\n");
   printf("      -b                   primed antiparallel\n");
   printf("      -h                   hairpin\n");
   printf("      -p                   parallel\n");
   printf("      -q                   primed parallel\n"); */
   printf("\n");
   exit(1); 
}



/**************************************************************************
 * 
 * read_pdb - Function to read the PDB file 
 *
 *************************************************************************/
int read_pdb(ATOMS **first_prot, ATOMS **last_prot, int *natoms, 
             ATOMS **first_water, ATOMS **last_water, int *nwaters)
{
/* Declare variables ====================================================*/
  int done, hydrogen;
  char    atom_name[5], input_line[LINE_LEN];
  ATOMS   *atom, *water;
  FILE    (*pdbfile);

/* Initialize variables =================================================*/
  atom = water = NULL;
  done = FALSE;

/* Open files ===========================================================*/
  pdbfile = open_file(pdbfilename, "r");

/* Read PDB file ========================================================*/
  while (fgets(input_line, LINE_LEN, pdbfile) != NULL && done == FALSE)
    {
      /* Check whether this is a hydrogen atom */
      if (!strncmp(input_line, "ATOM", 4) ||
          !strncmp(input_line, "HETATM", 6))
        {
          strncpy(atom_name, input_line+12, 4);
          atom_name[4] = '\0';
          hydrogen = is_hydrogen_atom(atom_name);
        }
      else hydrogen = FALSE;

      /* Get water information */
      if ((!strncmp(input_line, "ATOM", 4) ||
           !strncmp(input_line, "HETATM", 6)) && 
          !strncmp(input_line+17, "HOH", 3) && hydrogen == FALSE)
        {
          /* memory allocation */
          water = (ATOMS *) malloc(sizeof(ATOMS));
          if (water == NULL)
            {
              printf("*** ERROR - Unable to store water atom information\n");
              exit(1);
            }
         
          /* Initialize or update pointers */
          if (*first_water == NULL)
            {
              *first_water = water;
              water->prev_atom = NULL;
            }
          else
            {
              (*last_water)->next_atom = water;
              water->prev_atom = *last_water;
            }

          get_atom_info(input_line, water);
          (*nwaters)++;
         
          /* Update pointers */
          water->next_atom = NULL;
          *last_water = water;
        }

      /* Get protein and DNA information */
      else if ((!strncmp(input_line, "ATOM", 4) ||
                !strncmp(input_line, "HETATM", 6)) && hydrogen == FALSE)
        {
          /* Memory allocation */
          atom = (ATOMS *) malloc(sizeof(ATOMS));
          if (atom == NULL)
            {
              printf("*** ERROR - Unable to store protein atom information\n");
              exit(1);
            }

          /* Initialize or update pointers */
          if (*first_prot == NULL)
            {
              *first_prot = atom;
              atom->prev_atom = NULL;
            }
          else
            {
              (*last_prot)->next_atom = atom;
              atom->prev_atom = *last_prot;
            }

          /* Get the info */
          get_atom_info(input_line, atom);
          (*natoms)++;

          /* Update pointers */
          atom->next_atom = NULL;
          *last_prot = atom;
        }

      /* Check for end of first model in an NMR structure */
      else if (!strncmp(input_line, "ENDMDL", 6))
        done = TRUE;
    }

/* Check that ATOM records were found ===================================*/
  if (atom == NULL)
    {
      printf("*** ERROR - No ATOM records found in PDB file\n");
      exit(1);
    }
   
/* Close files ==========================================================*/
  fclose(pdbfile);

/* Return to main =======================================================*/
  return(0);
}



/***********************************************************************

is_hydrogen_atom  -  Check whether this is a hydrogen or a pseudo atom
                     from any of the many possible representations of
                     such in the PDB

***********************************************************************/

int
is_hydrogen_atom(char atname[5])
{
  char atom_name[5];

  int is_hydrogen;

  /* Initialise */
  is_hydrogen = FALSE;
  strncpy(atom_name,atname,4);
  atom_name[4] = '\0';

  /* Check the atom name for any indicators to suggest this is either
     a hydrogen or a pseudo atom */
  if (atom_name[0] == 'H' || atom_name[1] == 'H' || atom_name[1] == 'Q' ||
      (atom_name[1] == 'D' && (atom_name[0] == ' ' || atom_name[0] == '1' ||
                               atom_name[0] == '2' || atom_name[0] == '3')))
    is_hydrogen = TRUE;

  /* Check for one of the possible metals */
  if (!strncmp(atom_name,"HG  ",4) || !strncmp(atom_name,"HO  ",4))
    is_hydrogen = FALSE;

  return is_hydrogen;
}



/**************************************************************************
 *
 * get_atom_info - Function to get all the atom details 
 *
 *************************************************************************/
int get_atom_info(char input_line[LINE_LEN], ATOMS *atom)
{
/* Declare variables ====================================================*/
   int i;
   int atom_number;
   double  x, y, z;

/* Store atom details ===================================================*/
   /* Atom number */
/* v.1.1.10--> */
//   sscanf(input_line+7, "%d", &atom_number);
   sscanf(input_line+6, "%d", &atom_number);
/* <--v.1.1.10 */
   atom->number = atom_number;
   
   /* Atom name */
   strncpy(atom->name, input_line+12, 4);
   atom->name[4] = '\0';
   
   /* Residue name */
   strncpy(atom->resname, input_line+17, 3);
   atom->resname[3] = '\0';

   /* Chain ID */
   atom->chain = input_line[21];
   
   /* Sequence ID */
   sscanf(input_line+22, "%4d", &(atom->resnum));
   
   /* Coordinates */
   x = y = z = 0;
   sscanf(input_line+30, "%lf %lf %lf", &x, &y, &z);
   atom->coord[0] = x;
   atom->coord[1] = y;
   atom->coord[2] = z;
   
/* Initialize unused variables ==========================================*/
   atom->pair_name[0] = '\0';
   atom->pair_resnum = -1;
   atom->pair_resname[0] = '\0';
   atom->pair_chain  = '\0';
   atom->numbering   = -1;
   atom->nbonds      = 0;
   atom->dna_chain   = FALSE;
   atom->dna_chain   = -1;
   atom->npair_bond  =  0;
   atom->bridge_flag = FALSE;
   atom->pair        = NULL;
   atom->pair_dna    = NULL;
   atom->next_res    = NULL;
   atom->node_ptr    = NULL;
   for (i = 0; i < MAXARRAY; i++)
     {
       atom->bondto[i] = NULL;
       atom->bonded_atom[i] = NULL;
       atom->bonded_pair[i] = NULL;
     }
   
/* Return to read_pdb ===================================================*/
   return(0);
}



/**************************************************************************
 *
 * split - Function to distinguish between the atom types into proteins and
 *         each DNA strand
 *
 *************************************************************************/
int split(ATOMS *first_atom, ATOMS *first_prot[MAXARRAY],
          ATOMS *last_prot[MAXARRAY], ATOMS *first_dna[MAXARRAY],
          ATOMS *last_dna[MAXARRAY], int *nprot, int *ndna, PARAM *par)
{
/* Declare variables ====================================================*/
  int   last_resnum;
  int   dna_chain, dna_base, new_chain;
  int   nphos, non_phos;
  char  chain, code, last_resname[4];
  ATOMS *atom;
   
/* Initialize variables =================================================*/
  *nprot = *ndna = 0;
  atom   = first_atom;
  chain = '\0';
  dna_chain = FALSE;
  new_chain = TRUE;
  last_resnum = 0;
  last_resname[0] = '\0';
  nphos = 0;
  non_phos = 0;

/* Determine whether first residue is a nucleic acid or not */
  dna_chain = is_nucleic_acid(atom->resname, &code, par);
  dna_base = dna_chain;

/* Loop through atoms ===================================================*/
  while (atom != NULL)
    {
      /* If chain has not changed, then check that residue-type is
         consistent with the current chain type */
      if (atom->chain == chain)
        {
          /* Initialise chain-break flag */
          new_chain = FALSE;

          /* Check whether this is a DNA base */
          dna_base = is_nucleic_acid(atom->resname, &code, par);

          /* Check whether residue-type is consistent */
          if (dna_chain != dna_base)
            {
              /* Show warning message */
              if (dna_chain == TRUE)
                {
                  printf("*** Warning. Unusual DNA base encountered %s\n",
                         atom->resname);
                }
              else
                {
                  printf("*** Warning. DNA base encountered in protein %s\n",
                         atom->resname);
                }

              /* Create a chain-break here */
              new_chain = TRUE;
            }

          /* If the base has changed, check that starts with a phosphor */
          if (dna_chain == TRUE && (atom->resnum != last_resnum ||
                   strncmp(atom->resname,last_resname,3)))
            {
              if (atom->name[1] != 'P')
                new_chain = TRUE;
            }
        }

      /* Assign new starting pointer when chain ID changes */
      if (atom->chain != chain || new_chain == TRUE)
        {
          /* Determine whether this is a DNA base */
          dna_chain = is_nucleic_acid(atom->resname, &code, par);

          /* Create a new DNA chain-start */
          if (dna_chain == TRUE)
            {
              /* Check maximum not exceeded */
              if (*ndna >= MAXARRAY - 1)
                {
                  printf("*** ERROR - Exceeded maximum number\n");
                  printf("            of DNA strands: %d\n",MAXARRAY);
                  exit(1);
                }

              /* Otherwise, store chain-start */
              else
                {
                  first_dna[*ndna] = atom;
                  (*ndna)++;
                }
            }

          /* Create a new protein chain-start */
          else
            {
              /* Check maximum not exceeded */
              if (*nprot >= MAXARRAY - 1)
                {
                  printf("*** ERROR - Exceeded maximum number\n");
                  printf("            of protein chains: %d\n",MAXARRAY);
                  exit(1);
                }

              /* Otherwise, store chain-start */
              else
                {
                  first_prot[*nprot] = atom;
                  (*nprot)++;
                }
            }

          /* Store the current chain */
          chain = atom->chain;
        }

      /* Store the current pointer as the last of the current chain */
      if (dna_chain == TRUE)
        {
          last_dna[(*ndna) - 1] = atom;
          atom->type = DNA;
          atom->dna_chain = (*ndna) - 1;

          /* If DNA, then increment counts of atom-types */
          if (atom->name[1] == 'P')
            nphos++;
          else
            non_phos++;
        }
      else
        {
          last_prot[(*nprot) - 1] = atom;
          atom->type = PROTEIN;

          /* Check whether this is a water */
          if (!strncmp(atom->resname, "HOH", 3))
            atom->type = WATER;
        }

      /* Store the residue name and number */
      last_resnum = atom->resnum;
      strncpy(last_resname,atom->resname,3);
      last_resname[3] = '\0';

      /* Get the next atom */
      atom = atom->next_atom;
    }
   
  /* Check that DNA strands are phosphate-only */
  if ((*ndna) > 0 && non_phos < 6 * nphos)
    {
      printf("\n");
      printf("*** Incomplete nucleic acid chains found in the PDB file.\n");
      printf("*** No plot will be produced.\n");
      printf("\n");
      exit(1);
    }

/* Return to main =======================================================*/
  return(0);
}



/***********************************************************************

is_nucleic_acid  -  Check whether the given 3-letter residue name is a
                    nucleic acid

***********************************************************************/

int
is_nucleic_acid(char resname[4], char *code, PARAM *par)
{
  int ibase;
  int nucleic_acid;

  char nucname[NUMBASES][4] = {
    "  A", " +A", "ADE", "  C", " +C", "CYT", "  G", " +G", "GUA",
/* v.1.1.3--> */
//    "  T", " +T", "THY", "  U", "  I" }; 
    "  T", " +T", "THY", "  U", "  I",
/* v.1.1.9--> */
//  " DA", " DC", " DG", " DI", " DT", " DU" }; 
    " DA", " DC", " DG", " DI", " DT", " DU",
    "DA ", "DC ", "DG ", "DI ", "DT ", "DU " }; 
/* <--v.1.1.9 */
/* <--v.1.1.3 */
  char one_letter[NUMBASES] = {
      'A',   'A',   'A',   'C',   'C',   'C',   'G',   'G',   'G',
/* v.1.1.3--> */
//      'T',   'T',   'T',   'U',   'I' }; 
      'T',   'T',   'T',   'U',   'I',
/* v.1.1.9--> */
//      'A',   'C',   'G',   'I',   'T',   'U', }; 
      'A',   'C',   'G',   'I',   'T',   'U',
      'A',   'C',   'G',   'I',   'T',   'U' }; 
/* <--v.1.1.9 */
/* <--v.1.1.3 */

  /* Initialise */
  nucleic_acid = FALSE;
  *code = '\0';

  /* Loop through all the nucleic acid bases */
  for (ibase = 0; ibase < NUMBASES && nucleic_acid == FALSE; ibase++)
    if (!strncmp(resname, nucname[ibase], 3))
      {
        /* Match found. Get equivalent one-letter code, if there is one */
        nucleic_acid = TRUE;
        if (one_letter[ibase] != ' ')
          *code = one_letter[ibase];
      }

  /* If not found in the standard names, check the non-standard ones */
  if (nucleic_acid == FALSE)
    for (ibase = 0; ibase < MAX_NON_STANDARD && nucleic_acid == FALSE; ibase++)
      if (!strncmp(resname, par->nonstandard_name[ibase], 3))
        nucleic_acid = TRUE;

  /* Return answer */
  return nucleic_acid;
}



/**************************************************************************
 *
 * baselink - Function link to the previous and next base.
 *            next_res - link to the first atom of the next base
 *            prev_res - link to the last atom of the previous base
 *
 *************************************************************************/
int
baselink(ATOMS *first_atom, ATOMS *last_atom)
{
/* Declare variables ====================================================*/
   int   lastresnum;
   ATOMS *atom, *last_res, *this_res;

/* Initialize variables =================================================*/
/* v.1.1.4--> */
/*   lastresnum = -1; */
   lastresnum = -999999;
/* <--v.1.1.4 */
   last_res   = this_res = NULL;
   atom       = first_atom;
   

/* Loop through the strand ==============================================*/
   while (atom != NULL)
   {
      if (atom->resnum != lastresnum) 
      {
         /* Point first atom of last res to first atom of new res */
         if (last_res != NULL)
            last_res->next_res = atom;
         
         /* Update pointer for first atom of last res */
         last_res = atom;
         lastresnum = atom->resnum;
      }
      atom = atom->next_atom;
   }

/* v.1.1.4--> */
/*   lastresnum = -1; */
   lastresnum = -999999;
/* <--v.1.1.4 */
   atom       = first_atom;
   while (atom != NULL)
   {
      /* If meet new base */
      if (atom->resnum != lastresnum)
      {
         this_res = atom;
         lastresnum = atom->resnum;
      }
      
      atom->next_res = this_res->next_res;
      
      atom = atom->next_atom;
   }

/* Loop through the strand backwards ====================================*/
/* v.1.1.4--> */
/*   lastresnum = -1; */
   lastresnum = -999999;
/* <--v.1.1.4 */
   atom = last_res = last_atom;
   
   while (atom != NULL)
   {
      /* If meet new base */
      if (atom->resnum != lastresnum) 
      {
         /* Point first atom of last res to first atom of new res */
         last_res->prev_res = atom;
            
         /* Update pointer for first atom of last res */
         last_res = atom;
         lastresnum = atom->resnum;
      }
      atom = atom->prev_atom;
   }

/* v.1.1.4--> */
/*   lastresnum = -1; */
   lastresnum = -999999;
/* <--v.1.1.4 */
   atom       = last_atom;
   while (atom != NULL)
   {
      if (atom->resnum != lastresnum)
      {
         this_res = atom;
         lastresnum = atom->resnum;
      }
      
      atom->prev_res = this_res->prev_res;
      atom = atom->prev_atom;
   }
   
/* Return to main =======================================================*/
   return(0);
}



/**************************************************************************
 *
 * get_bond_info - Function to get bonding information 
 *
 *************************************************************************/
int get_bond_info(ATOMS *first_dna[MAXARRAY], ATOMS *last_dna[MAXARRAY],
                  ATOMS *first_prot[MAXARRAY], ATOMS *last_prot[MAXARRAY], 
                  ATOMS *first_water, BONDS **first_bond, BONDS **first_nbond, 
                  BONDS **first_connect, BONDS **last_bond,
                  PARAM *par, int ndna, int nprot, int argc, char *argv[])
{
/* Declare variables ====================================================*/
   int bond_flag;

/* Initialize variables =================================================*/
   bond_flag = FALSE;

/* Check if have bond file or not =======================================*/
   if ((bondfile = fopen(bondfilename, "r")) != NULL &&
       par->rd_bond_file == TRUE)
      bond_flag = TRUE;
   else
      bond_flag = FALSE;
/* v.1.1--> */
   if (bondfile != NULL)
/* <--v.1.1 */
     fclose(bondfile);

/* Read .hb2 file =======================================================*/
   readhbplus(first_dna, first_prot, first_water, &(*first_bond),
              &(*last_bond), par, ndna, nprot, bond_flag);

/* If have .bond file ===================================================*/
   if (bond_flag == TRUE)
   {
      printf("Bond file supplied - read bonding information from %s\n",
             bondfilename);
      read_bond_file(&(*first_bond), &(*first_nbond), &(*first_connect), 
                     &(*last_bond), first_dna, first_prot, first_water,
                     par, ndna, nprot);
   }
   
/* If not have .bond file ===============================================*/
   else if (bond_flag == FALSE)
   {
      printf("No bond file supplied:- \n");
      printf("    Reading external interaction information from PDB CONNECT records\n");
      printf("    and HBPLUS output files\n");      

      /* Get the non bonded contact information */
      if (par->non_bond_contact == TRUE)
         readnb2(first_dna, first_prot, first_water, &(*first_bond),
                 &(*last_bond), &(*first_nbond), par, ndna, nprot,
                 argc, argv);

      /* Get covalent bonding information */
      read_connect_records(first_prot, last_prot, first_dna, last_dna, 
                           &(*first_bond), &(*last_bond), &(*first_connect),
                           nprot, ndna);

      /* Write bond interactions file */
      write_bonds(*first_bond, *first_nbond, *first_connect);
   }

/* Exclude unwanted bonds ===============================================*/
   if (par->prot_res == FALSE
       || par->non_bridge_water == FALSE 
       || par->bridge_water == FALSE)
      exclude_bonds(&(*first_bond), &(*last_bond), par);

/* Return to main =======================================================*/
   return(0);
}



/**************************************************************************
 *
 * readhbplus - Open and read the HBPLUS output file 
 *
 *************************************************************************/
int readhbplus(ATOMS *first_dna[MAXARRAY], ATOMS *first_prot[MAXARRAY], 
               ATOMS *first_water, BONDS **first_bond, BONDS **last_bond,
               PARAM *par, int ndna, int nprot, int bond_flag)
{
/* Declare variables ====================================================*/
   int    i, j, k;
   int    flag;
   int    nchar;
/* v.1.1.5--> */
   int headers = TRUE;
/* <--v.1.1.5 */
   char   input_line[LINE_LEN];
   char   lastres;
   double distance;
   ATOMS  *atom;
   FILE   (*hbfile);

/* Initialize variables =================================================*/
   i = j   = k = 0;
   flag    = FALSE;
   nchar   = 0;
   atom    = NULL;
   lastres = '\0';


/* Open the HBPLUS file =================================================*/
   /* Open file */
   hbfile = open_file(hb2filename, "r");

/* Read the HBPLUS file =================================================*/
   /* Check the HBPLUS file */
   check_hbplus(hbfile);

   /* Get DNA base pair and protein-DNA hbond info first */
   while (fgets(input_line, LINE_LEN, hbfile) != NULL)
   {
/* v.1.1.5--> */
     /* If reading in header records, check whether this is the last one */
     if (headers == TRUE)
       {
         /* Check for last of the headers */
         if (!strncmp(input_line,"n    s   type  num",18))
           headers = FALSE;
       }

     /* Otherwise, process this data line */
     else
       {
/* <--v.1.1.5 */

         sscanf(input_line+28, "%lf", &distance);

         /* Get the DNA base pairing information ****See lab notes 2**** */
         dna_base_pair(input_line, first_dna, ndna, par);

         /* Get Protein-DNA Hbonding information - only if no .bond file */
         if (bond_flag == FALSE)
           if (distance <= par->hb_distance)
             prot_dna_hbond(input_line, first_prot, first_dna, &(*first_bond), 
                            &(*last_bond), nprot, ndna, par);
/* v.1.1.5--> */
       }
/* <--v.1.1.5 */
   }

   /* Get other Hbond interactions if no .bond file */
   if (bond_flag == FALSE)
   {
      /* Stack DNA-water interactions below protein-DNA interactions */
      rewind(hbfile);
      while (fgets(input_line, LINE_LEN, hbfile) != NULL)
      {
         sscanf(input_line+28, "%lf", &distance);
         if (distance <= par->hb_distance)
            water_dna_hbond(input_line, first_water, first_dna,
                            &(*first_bond), &(*last_bond), ndna, par);
      }

      /* Stack Protein-Water H-bonding information bridging waters below
         DNA_Water interactions */
      rewind(hbfile);
      while (fgets(input_line, LINE_LEN, hbfile) != NULL)
      {
         sscanf(input_line+28, "%lf", &distance);
         if (distance <= par->hb_distance)
            prot_water_hbond(input_line, first_prot, first_water,
                             &(*first_bond), &(*last_bond), nprot, par);
      }   
   
      /* Loop through bonds to delete water interactions as specified
         in parameter file */
      /*if (par->prot_res == FALSE
          || par->non_bridge_water == FALSE 
          || par->bridge_water == FALSE)
         exclude_bonds(&(*first_bond), &(*last_bond), par);*/
   }
   
/* Close HBPLUS output file =============================================*/
   fclose(hbfile);

/* Return to main =======================================================*/
   return(0);
}



/***********************************************************************

getnam  -  Peel off the directory path and extension from the full name
           of the .pdb file

***********************************************************************/

int getnam(char pdbfil[FILENAME_LEN],int *istart,int *iend)
{
    char pchar;
    int  finish, gotdot, ipos, istate;


    /* Initialise variables */
    finish = FALSE;
    *iend = -1;
    *istart = 0;
    istate = 1;
    ipos = strlen(pdbfil) - 1;
    gotdot = FALSE;

    /* Check through the filename from right to left */
    while (finish == FALSE && ipos > -1)
    {

        /* Pick off next character */
        pchar = pdbfil[ipos];


        /* State 1: Searching for first non-blank character */
        if (istate == 1)
        {
            if (pchar == '/' || pchar == '\\' || pchar == ']')
            {
                printf("*** ERROR in supplied name of file: [%s\n",
                       pdbfil);
                return(-1);
            }

            if (pchar != ' ' && pchar != '.')
            {
                *iend = ipos;
                istate = 2;
            }
        }

        /* State 2: Searching for end of extension, or end of
           directory path */
        else if (istate == 2)
        {
            /* If character is a dot, and is the first dot, then
               note position */
            if (pchar == '.' && gotdot == FALSE)
            {
                *iend = ipos - 1;
                gotdot = TRUE;
            }

            /* If character signifies the end of a directory path,
               note pstn */
            else if (pchar == '/' || pchar == '\\' || pchar == ']')
            {
                *istart = ipos + 1;
                finish = TRUE;
            }
        }

        /* Step back a character */
        ipos = ipos - 1;
    }

    /* Check whether file name is sensible */
    if (*istart > *iend)
    {
        printf("*** ERROR in supplied name of file: %s\n", pdbfil);
        printf("    No name found\n");
        return(-1);
    }
    return(0);
}
/***********************************************************************

get_filename  -  Strip off directory path of PDB filename
                 and append required extention to root

***********************************************************************/

void get_filename(char *pdbfilename, char *file_name, char *root_name,
                  char *extension)
{
    int      iend, ierror, inew, ipos, istart;

    /* Initialise returned filename */
    file_name[0] = '\0';

    /* Check that PDB filename is valid */
    if (pdbfilename[0] != '\0' && pdbfilename[0] != '\n')
    {
        /* Peel off the directory path and extension, returning just
           the root file name */
        ierror = getnam(pdbfilename,&istart,&iend);
        if (ierror != -1)
        {
            /* Transfer just the root name into new file name */
            for (ipos = istart, inew = 0; ipos < iend + 1;
                 ipos++, inew++)
                file_name[inew] = pdbfilename[ipos];
            file_name[inew] = '\0';

            /* Store the root name */
            strcpy(root_name,file_name);

            /* Append the extension */
            strcat(file_name,extension);
        }
    }
}


/**************************************************************************
 *
 * check_hbplus - Function to check the hbplus output file
 *
 *************************************************************************/
int check_hbplus(FILE *hbfile)
{
/* Declare variables ====================================================*/
   char  input_line[LINE_LEN];
   
/* Check 1st line to see it is a hb2 file ===============================*/
   if (fgets(input_line, LINE_LEN, hbfile) != NULL)
   {
      if (strncmp(input_line, "HBPLUS Hydrogen Bond Calculator", 31))
      {
         printf("\n*** ERROR - %s not recognized as a HBPLUS output file\n",
                hb2filename); 
         printf("            File must contain the following phrase\n");
         printf("            in line 1, position 1:\n");
         printf("                \"HBPLUS Hydrogen Bond Calculator\"\n");
         exit(1);
      }
   }

/* Check it isn't an empty file =========================================*/
   else if (fgets(input_line, LINE_LEN, hbfile) == NULL)
   {
      printf("\n*** ERROR - No data found in HBPLUS output file\n");
      exit(1);   
   }

/* Return to readhbplus =================================================*/
   return(0);
}



/**************************************************************************
 *
 * dna_base_pair - Function to obtain the DNA base pairing information
 *
 *************************************************************************/
int dna_base_pair(char input_line[LINE_LEN], ATOMS *first_dna[MAXARRAY],
                  int ndna, PARAM *par)
{
/* Declare variables ====================================================*/
   int i, j, k;
   int dna_baseA, dna_baseB, skip;
   int   resA_no, resB_no;
   char code;
   char atomnameA[5], atomnameB[5], chainA, chainB, resnameA[4], resnameB[4];

   ATOMS *atom;
   
/* Initialize variables =================================================*/
   i = j = k = 0;
   atom  = NULL;

/* Get the relevant details from the input line */
   strncpy(atomnameA,input_line+9,4);
   atomnameA[4] = '\0';
   strncpy(atomnameB,input_line+23,4);
   atomnameB[4] = '\0';
   strncpy(resnameA,input_line+6,3);
   resnameA[3] = '\0';
   strncpy(resnameB,input_line+20,3);
   resnameB[3] = '\0';
   sscanf(input_line+1, "%4d", &resA_no);
   sscanf(input_line+15, "%4d", &resB_no);
   chainA = input_line[0];
   chainB = input_line[14];

   /* Check whether this is a DNA base */
   dna_baseA = is_nucleic_acid(resnameA, &code, par);

   /* Check whether other one is a DNA base */
   dna_baseB = is_nucleic_acid(resnameB, &code, par);

   /* Check that H-bond not between backbone atoms */
   skip = FALSE;
   if (atomnameA[3] == '*' || atomnameB[3] == '*')
     skip = TRUE;

/* Check that H-bond is between DNA bases ================================
   & it is between different strands ===================================*/
   if (dna_baseA == TRUE && dna_baseB == TRUE && skip == FALSE)
     {
       /* Loop through the DNA strands to mark the pairing */
       for (k = 0; k < ndna; k++)
         {
           /* Get the pointer to the first atom of this strand */
           atom = first_dna[k];

           /* Loop through all the atoms in this strand */
           while (atom != NULL && atom->dna_chain == k)
             {
               /* If atom's base matches first base of pair, then
                  store other base's details */
               if (!strncmp(atom->name, atomnameA, 4) &&
                   !strncmp(atom->resname, resnameA, 3) &&
                   atom->resnum == resA_no && atom->chain == chainA)
                 {
                   strncpy(atom->pair_name,atomnameB,4);
                   atom->pair_name[4] = '\0';
                   strncpy(atom->pair_resname,resnameB,3);
                   atom->pair_resname[3] = '\0';
                   atom->pair_resnum   = resB_no;
                   atom->pair_chain    = chainB;
                 }

               /* If atom's base matches second base of pair, then
                  store first one's details */
               if (!strncmp(atom->name, atomnameB, 4) &&
                   !strncmp(atom->resname, resnameB, 3) &&
                   atom->resnum == resB_no && atom->chain == chainB)
                 {
                   strncpy(atom->pair_resname,resnameA,3);
                   atom->pair_resname[3] = '\0';
                   strncpy(atom->pair_name,atomnameA,4);
                   atom->pair_name[4] = '\0';
                   atom->pair_resnum   = resA_no;
                   atom->pair_chain    = chainA;
                 }

               /* Get pointer to the next atom record in this strand */
               atom = atom->next_atom;
             }                              
         }
     }

/* Return to readhbplus =================================================*/
   return(0);
}



/**************************************************************************
 *
 * exclude_bonds - Function to delete water interactions from the
 *                 bonds linked list as specified in the parameter 
 *                 file
 *
 *************************************************************************/
int
exclude_bonds(BONDS **first_bond, BONDS **last_bond, PARAM *par)
{
/* Declare variables ====================================================*/
  ATOMS *atom[2];
  BONDS *bond, *next_bond, *prev_bond;
   
/* Initialize variables =================================================*/
  bond = *first_bond;
   
/* Loop through bonds ===================================================*/
  while (bond != NULL)
    {
      atom[0]   = bond->bond_atom[0];
      atom[1]   = bond->bond_atom[1];
      next_bond = bond->next_bond;
      prev_bond = bond->prev_bond;
             
/* Only for Hydrogen bonds ==============================================*/
      if (bond->type == HBOND)
        {
/*=========================================================================
 * For protein-DNA interactions
 *=======================================================================*/
          if (par->prot_res == FALSE)
            {
              /* Check if one atom is a DNA atom and the other a protein
                 atom */
              if ((atom[0]->type == DNA && atom[1]->type == PROTEIN) ||
                  (atom[1]->type == DNA && atom[0]->type == PROTEIN))
                {
                  if (prev_bond != NULL)
                    prev_bond->next_bond = next_bond;

                  if (next_bond != NULL)
                    {
                      next_bond->prev_bond = prev_bond;
                      if (bond == *first_bond)
                        *first_bond = next_bond;
                    }

                  if (next_bond == NULL)
                    if (bond == *first_bond)
                      *first_bond = NULL;

                  if (*first_bond == NULL)
                    *last_bond = NULL;
                }
            }
/*=========================================================================      
 * For water interactions
 *=======================================================================*/
          if (!strncmp(atom[0]->resname, "HOH", 3)
              || !strncmp(atom[1]->resname, "HOH", 3))
            {
              
/* Non-bridging interactions ============================================*/ 
              if (((!strncmp(atom[0]->resname, "HOH", 3) &&
                    atom[0]->bridge_flag == FALSE)
                   || (!strncmp(atom[1]->resname, "HOH", 3) &&
                       atom[1]->bridge_flag == FALSE))
                      && par->non_bridge_water == FALSE)
                {
                  if (prev_bond != NULL)
                    prev_bond->next_bond = next_bond;
                  if (next_bond != NULL)
                    {
                      next_bond->prev_bond = prev_bond;
                          if (bond == *first_bond)
                            *first_bond = next_bond;
                    }
                  if (next_bond == NULL)
                    if (bond == *first_bond)
                      *first_bond = NULL;
                  if (*first_bond == NULL)
                    *last_bond = NULL;
                }
            
/* Bridging water interactions ==========================================*/
              if ((atom[0]->bridge_flag == TRUE ||
                   atom[1]->bridge_flag == TRUE)
                  && par->bridge_water == FALSE)
                {
                  if (prev_bond != NULL)
                    prev_bond->next_bond = next_bond;
                  if (next_bond != NULL)
                    {
                      next_bond->prev_bond = prev_bond;
                      if (bond == *first_bond)
                        *first_bond = next_bond;
                    }
                  if (next_bond == NULL)
                    if (bond == *first_bond)
                      *first_bond = NULL;
                  if (*first_bond == NULL)
                    *last_bond = NULL;
                }
            }
        }
      bond = bond->next_bond;
    }

/* Return to readhbplus =================================================*/
  return(0);
}



/**************************************************************************
 *
 * prot_dna_hbond - Function to get the Hbonds between the protein and DNA
 *
 *************************************************************************/
int prot_dna_hbond(char input_line[LINE_LEN], ATOMS *first_prot[MAXARRAY], 
                   ATOMS *first_dna[MAXARRAY], BONDS **first_bond,
                   BONDS **last_bond, int nprot, int ndna, PARAM *par)
{
/* Declare variables ====================================================*/ 
   int   dna_flag, prot_flag, water_flag;
   int   found_flag;
   int   bondtype_flag;
   ATOMS *prot_res, *prot_atom; /* First atom of right res and right atom */
   ATOMS *dna_res, *dna_atom;  /* First atom of right base and right atom */
   

/* Initialize variables =================================================*/
   dna_flag      = prot_flag = -1;
   found_flag    = TRUE;
   bondtype_flag = HBOND;
   prot_res      = prot_atom  = dna_res = dna_atom = NULL; 

/* Check which types of atoms the bond is between */
   get_atom_types(input_line, &dna_flag, &prot_flag, &water_flag, par);

/* Record the H-bond if it is ===========================================*/
   if ((prot_flag == LEFT && dna_flag == RIGHT)
       || (prot_flag == RIGHT && dna_flag == LEFT))
   {
      /* Find the right protein residue */
      find_hbond_res(first_prot, nprot, &prot_res, &prot_atom, input_line,
                     prot_flag);
      
      /* Find the right DNA base */
      find_hbond_res(first_dna, ndna, &dna_res, &dna_atom, input_line,
                     dna_flag);

      /* Warning if HBPLUS atoms not found in PDB */
      hbond_warning(input_line, prot_res, prot_atom, dna_res, dna_atom,
/* v.1.1.2--> */
//                    &found_flag);
                    &found_flag,1);
/* <--v.1.1.2 */

      /* Proceed if found_flag = TRUE */
      if (found_flag == TRUE)
      {
         /* Pair up the hbonded atoms */
         hbond_list(&(*first_bond), &(*last_bond),prot_atom, dna_atom, 
                    prot_flag, dna_flag, bondtype_flag, input_line, PD);
      }
   }
   
/* Return to readhbplus =================================================*/
   return(0);
}



/**************************************************************************
 *
 * get_atom_types - Function to get the two types of atoms involved in the
 *                  given bond
 *
 *************************************************************************/
int
get_atom_types(char input_line[LINE_LEN], int *dna_flag, int *prot_flag,
               int *water_flag, PARAM *par)
{
  int left_type, right_type;
  int is_hydrogen;
   
  char atom_type[5], code, resname[4];

  /* Check whether left-hand atom is a hydrogen */
  strncpy(atom_type, input_line+9, 4);
  atom_type[4] = '\0';
  is_hydrogen = is_hydrogen_atom(atom_type);

  /* If atom is a hydrogen, then don't want it */
  if (is_hydrogen == TRUE)
    left_type = HYDROGEN;

  /* Otherwise, check residue-type the atom belongs to */
  else
    {
      /* Check left-hand atom for DNA base */
      strncpy(resname, input_line+6, 3);
      resname[3] = '\0';
      *dna_flag = is_nucleic_acid(resname, &code, par);

      /* If atom is a DNA atom, the store this */
      if (*dna_flag == TRUE)
        left_type = DNA;

      /* Otherwise is either a protein or a water atom */
      else
        if (!strncmp(input_line+6, "HOH", 3))
          left_type = WATER;
        else
          left_type = PROTEIN;
    }

  /* Repeat for right-hand atom */
  strncpy(atom_type, input_line+23, 4);
  atom_type[4] = '\0';
  is_hydrogen = is_hydrogen_atom(atom_type);

  /* If atom is a hydrogen, then don't want it */
  if (is_hydrogen == TRUE)
    right_type = HYDROGEN;

  /* Otherwise, check residue-type the atom belongs to */
  else
    {
      /* Get residue name */
      strncpy(resname, input_line+20, 3);
      resname[3] = '\0';
      *dna_flag = is_nucleic_acid(resname, &code, par);

      /* If atom is a DNA atom, the store this */
      if (*dna_flag == TRUE)
        right_type = DNA;

      /* Otherwise is either a protein or a water atom */
      else
        if (!strncmp(input_line+20, "HOH", 3))
          right_type = WATER;
        else
          right_type = PROTEIN;
    }

  /* Initialise the flags to be returned */
  *dna_flag = *prot_flag = *water_flag = -1;

/* Determine which side the different atoms are on */
  if (left_type == DNA)
    *dna_flag = LEFT;
  else if (left_type == PROTEIN)
    *prot_flag = LEFT;
  else if (left_type == WATER)
    *water_flag = LEFT;
  if (right_type == DNA)
    *dna_flag = RIGHT;
  else if (right_type == PROTEIN)
    *prot_flag = RIGHT;
  else if (right_type == WATER)
    *water_flag = RIGHT;

/* Return to prot_dna_hbond =============================================*/
   return(0);
}



/**************************************************************************
 *
 * water_dna_hbond - Function to get Hbonds between the DNA and water 
 *                   molecules 
 *
 *************************************************************************/
int water_dna_hbond(char input_line[LINE_LEN], ATOMS *first_water, 
                    ATOMS *first_dna[MAXARRAY], 
                    BONDS **first_bond, BONDS **last_bond, int ndna,
                    PARAM *par)
{
/* Declare variables ====================================================*/ 
   int   dna_flag, prot_flag, water_flag;
   int   found_flag;
   int   bondtype_flag;
   ATOMS *water_res, *water_atom; /* Dummy variable and right atom */
   ATOMS *dna_res, *dna_atom; /* First atom of right base and right atom */
   

/* Initialize variables =================================================*/
   dna_flag      = water_flag = -1;
   found_flag    = TRUE;
   bondtype_flag = HBOND;
   water_res     = water_atom  = dna_res = dna_atom = NULL; 

/* Check that H-bond is between water and DNA ===========================*/ 
   get_atom_types(input_line, &dna_flag, &prot_flag, &water_flag, par);
   
/* Record the H-bond if it is ===========================================*/
   if ((water_flag == LEFT && dna_flag == RIGHT)
       || (water_flag == RIGHT && dna_flag == LEFT))
   {
      /* Find right water molecule */
      find_hbond_water(first_water, &water_res, &water_atom, input_line,
                       water_flag);

      /* Find the right DNA base */
      find_hbond_res(first_dna, ndna, &dna_res, &dna_atom, input_line,
                     dna_flag);

      /* Warning if HBPLUS atoms not found in PDB */
      hbond_warning(input_line, water_res, water_atom, dna_res, dna_atom,
/* v.1.1.2--> */
//                    &found_flag);
                    &found_flag,2);
/* <--v.1.1.2 */

      /* Proceed if found_flag = TRUE */
      if (found_flag == TRUE)
      {
         /* Pair up the hbonded atoms */
         hbond_list(&(*first_bond), &(*last_bond), water_atom, dna_atom, 
                    water_flag, dna_flag, bondtype_flag, input_line, WD);

         /* Mark the details of the water-DNA Hbond for later check on
            bridging water */
         mark_water_hbond(water_atom, dna_atom);
      }
   }
   
  
/* Return to readhbplus =================================================*/
   return(0);
}



/**************************************************************************
 *
 * mark_water_bond - Function to mark the water hbond for later checking 
 *                   for bridging
 *
 *************************************************************************/
int mark_water_hbond(ATOMS *water_atom, ATOMS *dna_atom)
{
/* Record ===============================================================*/
   water_atom->bondto[water_atom->nbonds] = dna_atom;
   water_atom->nbonds++;

/* Return to water_dna_hbond ============================================*/
   return(0);
}



/**************************************************************************
 *
 * prot_water_hbond - Function to get Hbonds between the protein and 
 *                    bridging water molecules 
 *
 *************************************************************************/
int
prot_water_hbond(char input_line[LINE_LEN], ATOMS *first_prot[MAXARRAY],
                 ATOMS *first_water, BONDS **first_bond,
                 BONDS **last_bond, int nprot, PARAM *par)
{
/* Declare variables ====================================================*/ 
  int   dna_flag, prot_flag, water_flag;
  int   found_flag;
  int   bridge_flag;
   int   bondtype_flag;
  ATOMS *prot_res, *prot_atom; /* First atom of right res and right atom */
  ATOMS *water_res, *water_atom; /* Dummy variable and right atom */

/* Initialize variables =================================================*/
  prot_flag     = water_flag = -1;
  found_flag    = TRUE;
  bridge_flag   = FALSE;
  bondtype_flag = HBOND;
  prot_res      = prot_atom = water_res = water_atom = NULL; 

/* Check that H-bond is between water and protein =======================*/ 
  get_atom_types(input_line, &dna_flag, &prot_flag, &water_flag, par);

  if ((water_flag == LEFT && prot_flag == RIGHT)
      || (water_flag == RIGHT && prot_flag == LEFT))
    {
      /* Find right water molecule */
      find_hbond_water(first_water, &water_res, &water_atom, input_line,
                       water_flag); 
  
      /* Find the right protein residue */
      find_hbond_res(first_prot, nprot, &prot_res, &prot_atom, input_line,
                     prot_flag);
      
      /* Warning if HBPLUS atoms not found in PDB */
      hbond_warning(input_line, prot_res, prot_atom, water_res, water_atom,
/* v.1.1.2--> */
//                    &found_flag);
                    &found_flag,3);
/* <--v.1.1.2 */

      /* Proceed if found_flag = TRUE */
      if (found_flag == TRUE)
        {
          /* Check if the water is bridging */
         /*printf("%s\n", input_line);*/
          check_water_bridge(water_atom, &bridge_flag, par);
          
          /* Record the H-bond if it is bridging */
          if (bridge_flag == TRUE)
            {
              /* Pair up the hbonded atoms */
              hbond_list(&(*first_bond), &(*last_bond), prot_atom, water_atom, 
                       prot_flag, water_flag, bondtype_flag, input_line, PW);
            }
        }
    }
/* Return to readhbplus =================================================*/
  return(0);
}



/**************************************************************************
 *
 * find_hbond_res - Function to find the protein residues and DNA bases 
 *                  involved in the H-bond
 *
 *************************************************************************/
int
find_hbond_res(ATOMS *first_atom[MAXARRAY], int nchains, ATOMS **hbres, 
               ATOMS **hbatom, char input_line[LINE_LEN], int flag)
{
/* Declare variables ====================================================*/
  int   i;
  int   resnum[2];
  ATOMS *atom, *other_atom;

/* Get residue numbers from HBPLUS output file ==========================*/
  sscanf(input_line+1, "%4d", &resnum[LEFT]);
  sscanf(input_line+15, "%4d", &resnum[RIGHT]);

/* Check for blank chain-id */
  if (input_line[0] == '-')
    input_line[0] = ' ';
  if (input_line[14] == '-')
    input_line[14] = ' ';

/* Loop through chains ==================================================*/
  for (i = 0; i < nchains; i++)
    {
      if ((flag == LEFT && first_atom[i]->chain == input_line[0])
          || (flag == RIGHT && first_atom[i]->chain == input_line[14]))
        {
          /* If right chain, loop through residues */
          atom = first_atom[i];
          while ((atom != NULL) && (atom->chain == first_atom[i]->chain))
            {
              if ((flag == LEFT && atom->resnum == resnum[LEFT])
                  || (flag == RIGHT && atom->resnum == resnum[RIGHT]))
                {        
                  /* If right residue record to protres */
                  *hbres = atom;

                  /* Initialise pointer to first atom of this residue */
                  other_atom = atom;

                  /* Loop through atoms of the residue */
                  while (other_atom != NULL &&
                         other_atom->resnum == (*hbres)->resnum)
                    {
                      /* If right atom record to protatom */
                      if ((flag == LEFT &&
                           !strncmp(other_atom->name, input_line+9,4))
                          || (flag == RIGHT &&
                              !strncmp(other_atom->name, input_line+23, 4)))
                        *hbatom = other_atom;

                      /* Get pointer to the next atom */
                      other_atom = other_atom->next_atom;
                    }
                }
              if (atom != NULL)
                atom = atom->next_res;
            }
        }
    }

/* Return to prot_dna_hbond =============================================*/
  return(0);
}



/**************************************************************************
 *
 * find_hbond_water - Function to find the water molecule involved in 
 *                    the H-bond
 *
 *************************************************************************/
int find_hbond_water(ATOMS *first_water, ATOMS **water_res, ATOMS **water_atom,
                      char input_line[LINE_LEN], int flag)
{
/* Declare variables ====================================================*/
/* v.1.1.1--> */
  char chain[2];
/* <--v.1.1.1 */
   int   resnum[2];
   ATOMS *atom;

/* Get residue numbers from HBPLUS output file ==========================*/
   sscanf(input_line+1, "%4d", &resnum[LEFT]);
   sscanf(input_line+15, "%4d", &resnum[RIGHT]);
/* v.1.1.1--> */
   chain[LEFT] = input_line[0];
   chain[RIGHT] = input_line[14];
/* v.1.1.2--> */
   if (chain[LEFT] == '-')
     chain[LEFT] = ' ';
   if (chain[RIGHT] == '-')
     chain[RIGHT] = ' ';
/* <--v.1.1.2 */
/* <--v.1.1.1 */

/* Loop through chains ==================================================*/
   /* If right chain, loop through residues */
   atom = first_water;
   while (atom != NULL)
   {
/* v.1.1.1--> */
     // if ((flag == LEFT && atom->resnum == resnum[LEFT])
     //    || (flag == RIGHT && atom->resnum == resnum[RIGHT]))
     if ((flag == LEFT && atom->chain == chain[LEFT] &&
          atom->resnum == resnum[LEFT])
         || (flag == RIGHT && atom->chain == chain[RIGHT] &&
             atom->resnum == resnum[RIGHT]))
/* <--v.1.1.1 */
      {        
         /* If right molecule record to water_res & water_atom */
         *water_res = atom;
         *water_atom = atom;
      }
      atom = atom->next_atom;
   }

/* Return to prot_dna_hbond =============================================*/
   return(0);
}



/**************************************************************************
 *
 * hbond_warning - Function to produce a warning if the atom contained in 
 *                 the HBPLUS file is not found in the PDB
 *
 *************************************************************************/
int hbond_warning(char input_line[LINE_LEN], ATOMS *res1, ATOMS *atom1, 
/* v.1.1.2--> */
//                  ATOMS *res2, ATOMS *atom2, int *flag)
                  ATOMS *res2, ATOMS *atom2, int *flag, int call_no)
/* <--v.1.1.2 */
{
   if ((res1 == NULL) || (atom1 == NULL) ||
       (res2 == NULL) || (atom2 == NULL))
   {
      /* Warning */
      printf("*** Warning: Protein, DNA or water atom in the following");
      printf(" HPLUS file output\n");
/* v.1.1.2--> */
      printf("             line not found in PDB file.");
      printf(" Ignoring HBPLUS output for this line.\n"); 
      printf("             ");
      if (res1 == NULL)
        printf("1st res missing ");
      if (res2 == NULL)
        printf("2nd res missing ");
      if (atom1 == NULL)
        printf("1st atom missing ");
      if (atom2 == NULL)
        printf("2nd atom missing ");
      printf("[Call-pt: %d]\n",call_no);
/* <--v.1.1.2 */
      input_line[strlen(input_line)] = '\0';
      printf("%s\n", input_line);

      /* Update flag */
      *flag = FALSE;
   }
   return(0);
}



/**************************************************************************
 *
 * check_water_bridge - Function to check if the water bridges the protein 
 *                      and DNA 
 *
 *************************************************************************/
int check_water_bridge(ATOMS *water_atom, int *bridge_flag, PARAM *par)
{
/* Declare variables ====================================================*/ 
   int   i;
   int dna_flag;
   char code;
   ATOMS *bondto;

/* Initialize variables =================================================*/

/* Check bridging =======================================================*/
   for (i = 0; i < water_atom->nbonds; i++)
     {
       /* Get pointer to residue water is interacting with */
       bondto = water_atom->bondto[i];

       /* Check if the residue is a nucleic acid */
       dna_flag = is_nucleic_acid(bondto->resname, &code, par);

       /* If bond to nucleic acid, then store as a bridging water */
       if (dna_flag == TRUE)
         {
           (*bridge_flag) = TRUE;
           water_atom->bridge_flag = TRUE;
         }
     }
   
/* Return to prot_water_hbond ===========================================*/
   return(0);
}



/**************************************************************************
 *
 * hbond_list - Function to record the Hbond in a linked list
 *
 *************************************************************************/
int hbond_list(BONDS **first_bond, BONDS **last_bond, ATOMS *atom1,
               ATOMS *atom2, int flag1, int flag2, int bondtype,
               char input_line[LINE_LEN], int entity_flag)
{
/* Declare variables ====================================================*/ 
   double distance;
   BONDS  *bond;

   /* Check that we don't already have this bond */
   bond = *first_bond;
   while (bond != NULL)
     {
       /* If have this bond already, then return */
       if ((bond->bond_atom[DONOR] == atom1 &&
            bond->bond_atom[ACCEPTOR] == atom2) ||
           (bond->bond_atom[DONOR] == atom2 &&
            bond->bond_atom[ACCEPTOR] == atom1))
         return(0);

       bond = bond->next_bond;
     }

/* Memory allocation ====================================================*/
   bond = (BONDS *) malloc(sizeof(BONDS));
   if (bond == NULL)
   {
      printf("*** ERROR - Unable to store bonding information\n");
      exit(1);
   }
   
/* Initialize or update pointers ========================================*/
   if (*first_bond == NULL)
   {
      *first_bond = bond;
      bond->prev_bond = NULL;
   }
   else
   {
      (*last_bond)->next_bond = bond;
      bond->prev_bond = *last_bond;
   }
   
/* Record bonding information ===========================================*/
   bond->type = bondtype;
   if ((flag1 == LEFT) && (flag2 == RIGHT))
   {
      bond->bond_atom[DONOR]    = atom1;
      bond->bond_atom[ACCEPTOR] = atom2;
  }
   else if ((flag1 == RIGHT) && (flag2 == LEFT))
   {
      bond->bond_atom[DONOR]    = atom2;
      bond->bond_atom[ACCEPTOR] = atom1;
   }

/* Record bond distance =================================================*/
   sscanf(input_line+28, "%lf", &distance);
   bond->distance = distance;

/* Record entity of interacting groups ==================================*/
   bond->entity = entity_flag;

/* Update pointers ======================================================*/
   bond->next_bond = NULL;
   *last_bond = bond;
   
   return(0);
}



/**************************************************************************
 *
 * read_bond_file - Function to read the .bond file
 *
 *************************************************************************/
int read_bond_file(BONDS **first_bond, BONDS **first_nbond, 
                   BONDS **first_connect, BONDS **last_bond,
                   ATOMS *first_dna[MAXARRAY], ATOMS *first_prot[MAXARRAY], 
                   ATOMS *first_water, PARAM *par, int ndna, int nprot)
{
/* Declare variables ====================================================*/
   int  bond_type;
   
/* Initialize variables =================================================*/
   bond_type = -1;

/* Open .bond file ======================================================*/
   bondfile = open_file(bondfilename, "r");

/*========================================================================
 * Read through file 
 *======================================================================*/
/* Check bond file =====================================================*/
   check_bond_file(bondfile);

/* Get Hbonding info ===================================================*/
   get_hbond_info(&(*first_bond), &(*last_bond), first_dna, first_prot, 
                  first_water, par, ndna, nprot);

/* Get non bonded contact info ==========================================*/
   if (par->non_bond_contact == TRUE)
      get_nbond_info(&(*first_bond), &(*first_nbond), &(*last_bond),
                     first_dna,  first_prot, first_water, par, ndna,
                     nprot);
   
/* Get covalent bond info ===============================================*/
   get_cbond_info(&(*first_bond), &(*first_connect), &(*last_bond), first_dna, 
                  first_prot, first_water, par, ndna, nprot);

/* Close .bond file =====================================================*/
   fclose(bondfile);
   
/* Return ===============================================================*/
   return(0);
}



/**************************************************************************
 *
 * get_hbond_info - Function to get Hbond interaction information from the 
 *                  bond file 
 *
 *************************************************************************/
int get_hbond_info(BONDS **first_bond, BONDS **last_bond,
                   ATOMS *first_dna[MAXARRAY], ATOMS *first_prot[MAXARRAY], 
                   ATOMS *first_water, PARAM *par, int ndna, int nprot)
{
/* Declare variables ====================================================*/
   int    bond_type;
   float  distance;
   char   input_line[LINE_LEN];
   char   hb2_line[LINE_LEN];

/* Initialize variables =================================================*/   
   bond_type = HBOND;
   
/* Read through file ====================================================*/
   /* Get DNA base pair and protein-DNA hbond info first */
   while ((fgets(input_line, LINE_LEN, bondfile) != NULL)
          && bond_type == HBOND)
   {
      /* Update bond_type when get to other bond type data */
      if (!strncmp(input_line, "**** Non Bonded Contacts", 24)
          || !strncmp(input_line, "**** Covalent Bonds", 19))
         bond_type = -1;

      /* Convert input line into HBPLUS format */
      convert_bond_to_hb2(input_line, hb2_line);

      /* Record bond details if within distance constraints */
      sscanf(hb2_line+28, "%f", &distance);
      if (distance <= par->hb_distance)
         prot_dna_hbond(hb2_line, first_prot, first_dna, &(*first_bond), 
                        &(*last_bond), nprot, ndna, par);
   }

   /* Stack DNA-water interactions below protein-DNA interactions */
   bond_type = HBOND;
   rewind(bondfile);
   while (fgets(input_line, LINE_LEN, bondfile) != NULL)
   {
      /* Update bond_type when get to other bond type data */
      if (!strncmp(input_line, "**** Non Bonded Contacts", 24)
          || !strncmp(input_line, "**** Covalent Bonds", 19))
         bond_type = -1;

      /* Convert input line into HBPLUS format */
      convert_bond_to_hb2(input_line, hb2_line);

      sscanf(hb2_line+28, "%f", &distance);
      if (distance <= par->hb_distance)
         water_dna_hbond(hb2_line, first_water, first_dna, &(*first_bond),
                         &(*last_bond), ndna, par);
   }

   /* Stack Protein-Bridging Water interactions below DNA-Water interactions */
   rewind(bondfile);
   while (fgets(input_line, LINE_LEN, bondfile) != NULL)
   {
      /* Update bond_type when get to other bond type data */
      if (!strncmp(input_line, "**** Non Bonded Contacts", 24)
          || !strncmp(input_line, "**** Covalent Bonds", 19))
         bond_type = -1;

      /* Convert input line into HBPLUS format */
      convert_bond_to_hb2(input_line, hb2_line);

      sscanf(hb2_line+28, "%f", &distance);
      if (distance <= par->hb_distance)
         prot_water_hbond(hb2_line, first_prot, first_water, &(*first_bond),
                          &(*last_bond), nprot, par);
   }

   /* Loop through bonds to delete water interactions as specified in
      parameter file */
   /*if (par->prot_res == FALSE
       || par->non_bridge_water == FALSE 
       || par->bridge_water == FALSE)
      exclude_bonds(&(*first_bond), &(*last_bond), par);*/

/* Return ===============================================================*/
   return(0);
}



/**************************************************************************
 *
 * get_nbond_info - Function to get non bonded contact information from the 
 *                  bond file 
 *
 *************************************************************************/
int get_nbond_info(BONDS **first_bond, BONDS **first_nbond, BONDS **last_bond,
                   ATOMS *first_dna[MAXARRAY], ATOMS *first_prot[MAXARRAY], 
                   ATOMS *first_water, PARAM *par, int ndna, int nprot)
{
/* Declare variables ====================================================*/
   int    bond_type;
   double distance;
   char   input_line[LINE_LEN];
   char   hb2_line[LINE_LEN];

/* Initialize variables =================================================*/
   bond_type = -1;

/* Read through file ====================================================*/
   rewind(bondfile);
   while (fgets(input_line, LINE_LEN, bondfile) != NULL)
   {
      /* Update bond_type when get to other bond type data */
      if (!strncmp(input_line, "**** Hydrogen Bonds", 19)
          || !strncmp(input_line, "**** Covalent Bonds", 19))
         bond_type = -1;
      else if (!strncmp(input_line, "**** Non Bonded Contacts", 24))
         bond_type = NBOND;
      
      if (bond_type == NBOND)
      {
         /* Convert input line into HBPLUS format */
         convert_bond_to_hb2(input_line, hb2_line);

         /* Record bond details if within distance constraints */
         sscanf(hb2_line+28, "%lf", &distance);
         if (distance <= par->nb_distance)
            prot_dna_nbond(hb2_line, first_prot, first_dna, &(*first_bond), 
                           &(*last_bond), &(*first_nbond), nprot, ndna, par);
      }
   }

/* Return ===============================================================*/
   return(0); 
}



/**************************************************************************
 *
 * get_cbond_info - Function to get covalent bond information from the 
 *                  bond file 
 *
 *************************************************************************/
int get_cbond_info(BONDS **first_bond, BONDS **first_connect,
                   BONDS **last_bond, ATOMS *first_dna[MAXARRAY],
                   ATOMS *first_prot[MAXARRAY], ATOMS *first_water,
                   PARAM *par, int ndna, int nprot)
{
/* Declare variables ====================================================*/
   int    bond_type;
   char   input_line[LINE_LEN];

/* Initialize variables =================================================*/
   bond_type = -1;

/* Read through file ====================================================*/
   while (fgets(input_line, LINE_LEN, bondfile) != NULL)
   {
      /* Update bond_type when get to other bond type data */
      if (!strncmp(input_line, "**** Hydrogen Bonds", 19)
          || !strncmp(input_line, "**** Non Bonded Contacts", 24))
         bond_type = -1;
      else if (!strncmp(input_line, "**** Covalent Bonds", 19))
         bond_type = COVALENT;
      
      if (bond_type == COVALENT)
         prot_dna_cbond(input_line, first_prot, first_dna, &(*first_bond), 
                        &(*last_bond), &(*first_connect), nprot, ndna, par);
   }
/* Return ===============================================================*/
   return(0);
}



/**************************************************************************
 *
 * check_bond_file - Function to check the bond file 
 *
 *************************************************************************/
int check_bond_file(FILE *bondfile)
{
/* Declare variables ====================================================*/
   char input_line[LINE_LEN];

/* Read first line of bond file =========================================*/
   if (fgets(input_line, LINE_LEN, bondfile) != NULL)
   {
      /* Check version */ 
      if (strncmp(input_line, "NUCPLOT v.1.0  -  Bond file", 27))
      {
         printf("*** ERROR - Incorrect bond file.\n");
         printf("            Must contain \"NUCPLOT v.1.0  -  Bond file\" in the first line.\n");
         exit(-1);
      }
      
      /* Check file name */
      if (strncmp(input_line+29, bondfilename, strlen(bondfilename)))
      {
         printf("*** ERROR - Incorrect bond file.\n");
         printf("            The bond file is not for the pdbfile %s.\n", pdbfilename);
         exit(-1);
      }
   }
   
/* Return ===============================================================*/
   return(0);
}



/**************************************************************************
 *
 * convert_bond_to_hb2 - Function to convert the bond file input line to 
 *                       hbplus output file format
 *
 *************************************************************************/
int convert_bond_to_hb2(char input_line[LINE_LEN], char hb2_line[LINE_LEN])
{
/* Declare variables ====================================================*/
   int i, len;
   char   resname[2][4];
   char   chain[2];
   char   resnum[2][6];
   char   atomname[2][5];
   char   distance[5];

   /* Initialise variables */
   for (i = 0; i < LINE_LEN; i++)
     hb2_line[i] = ' ';
   hb2_line[33] = '\0';
   len = strlen(input_line);

   /* Proceed if this looks like a sensible line */
   if (len > 43)
     {
/* Get details ==========================================================*/
       /* Resname */
       strncpy(resname[DONOR], input_line, 3);
       resname[DONOR][3] = '\0';
       strncpy(resname[ACCEPTOR], input_line+21, 3);
       resname[ACCEPTOR][3] = '\0';
   
       /* Chain ID */
       chain[0] = input_line[4];
       chain[1] = input_line[25];

       /* Resnum */
       strncpy(resnum[DONOR], input_line+6, 5);
       resnum[DONOR][5] = '\0';
       strncpy(resnum[ACCEPTOR], input_line+27, 5);
       resnum[ACCEPTOR][5] = '\0';
   
       /* Atom name */
       strncpy(atomname[DONOR], input_line+13, 4);
       atomname[DONOR][4] = '\0';
       strncpy(atomname[ACCEPTOR], input_line+34, 4);
       atomname[ACCEPTOR][4] = '\0';

       /* Distance */
       strncpy(distance, input_line+41, 4);
       distance[4] = '\0';

/* Store in hb2_line ====================================================*/
       hb2_line[0] = chain[0];
       strncpy(hb2_line+1, resnum[DONOR], 5);
       strncpy(hb2_line+6, resname[DONOR], 3);
       strncpy(hb2_line+9, atomname[DONOR], 4);
       hb2_line[14] = chain[1];
       strncpy(hb2_line+15, resnum[ACCEPTOR], 5);
       strncpy(hb2_line+20, resname[ACCEPTOR], 3);
       strncpy(hb2_line+23, atomname[ACCEPTOR], 4);
       strncpy(hb2_line+28, distance, 4);
       hb2_line[33] = '\0';
     }

   /* Test print */

/* Return ===============================================================*/
   return(0);
}



/**************************************************************************
 *
 * readnb2 - Function to read the nb2 file
 *
 *************************************************************************/
int readnb2(ATOMS *first_dna[MAXARRAY], ATOMS *first_prot[MAXARRAY], 
            ATOMS *first_water, BONDS **first_bond, BONDS **last_bond,
            BONDS **first_nbond, PARAM *par, int ndna, int nprot, int argc,
            char *argv[])
{
/* Declare variables ====================================================*/
   char   input_line[LINE_LEN];
/* v.1.1.5--> */
   int headers = TRUE;
/* <--v.1.1.5 */
   double distance;
   FILE   (*nbfile);

/* Initialize variables =================================================*/

/* Open file ============================================================*/ 
   nbfile = open_file(nb2filename, "r");

/* Read the HBPLUS file =================================================*/
   /* Check the HBPLUS file */
   check_hbplus(nbfile);

   while (fgets(input_line, LINE_LEN, nbfile) != NULL)
   {
/* v.1.1.5--> */
     /* If reading in header records, check whether this is the last one */
     if (headers == TRUE)
       {
         /* Check for last of the headers */
         if (!strncmp(input_line,"n    s   type  num",18))
           headers = FALSE;
       }

     /* Otherwise, process this data line */
     else
       {
/* <--v.1.1.5 */
         sscanf(input_line+28, "%lf", &distance);
      
         /* Get Protein-DNA non bonding contact information */
         if (distance <= par->nb_distance)
           prot_dna_nbond(input_line, first_prot, first_dna, &(*first_bond), 
                          &(*last_bond), &(*first_nbond), nprot, ndna, par);
/* v.1.1.5--> */
       }
/* <--v.1.1.5 */
   }

/* Close HBPLUS output file =============================================*/
   fclose(nbfile);

/* Return to main =======================================================*/
   return(0);
}



/**************************************************************************
 *
 * prot_dna_nbond - Function to get the non bonded contacts between the 
 *                  protein and DNA
 *
 *************************************************************************/
int prot_dna_nbond(char input_line[LINE_LEN], ATOMS *first_prot[MAXARRAY], 
                   ATOMS *first_dna[MAXARRAY], BONDS **first_bond,
                   BONDS **last_bond, BONDS **first_nbond, int nprot,
                   int ndna, PARAM *par)
{
/* Declare variables ====================================================*/ 
  int   dna_flag, prot_flag, water_flag;
  int   found_flag;
  int   bondtype_flag, wanted;
  ATOMS *prot_res, *prot_atom; /* First atom of right res and right atom */
  ATOMS *dna_res, *dna_atom;  /* First atom of right base and right atom */
   

/* Initialize variables =================================================*/
  dna_flag         = prot_flag = -1;
  found_flag       = TRUE;
  bondtype_flag    = NBOND;
  prot_res         = prot_atom  = dna_res = dna_atom = NULL; 

/* Check that N-bond is between protein and DNA =========================*/ 
  get_atom_types(input_line, &dna_flag, &prot_flag, &water_flag, par);

/* Record the N-bond if it is ===========================================*/
  if ((prot_flag == LEFT && dna_flag == RIGHT)
      || (prot_flag == RIGHT && dna_flag == LEFT))
    {
      /* Find the right protein residue */
      find_hbond_res(first_prot, nprot, &prot_res, &prot_atom, input_line,
                     prot_flag);
      
      /* Find the right DNA base */
      find_hbond_res(first_dna, ndna, &dna_res, &dna_atom, input_line,
                     dna_flag);

      /* Warning if HBPLUS atoms not found in PDB */
      hbond_warning(input_line, prot_res, prot_atom, dna_res, dna_atom,
/* v.1.1.2--> */
//                    &found_flag);
                    &found_flag,4);
/* <--v.1.1.2 */
      
      /* Proceed if found_flag = TRUE */
      if (found_flag == TRUE)
        {
          /* Check whether residue already involved in a H-bond */
          wanted = accept_reject_nbond(*first_bond, prot_atom);

          /* Pair up the non-bonded atoms */
          if (wanted == TRUE)
            nbond_list(&(*first_bond), &(*last_bond), &(*first_nbond),
                       prot_atom, dna_atom,  bondtype_flag, input_line, PD);
        }
    }
  
/* Return to readnb =====================================================*/
   return(0);
}



/**************************************************************************
 *
 * accept_reject_nbond - Function to determine whether to accept or reject
 *                       this non-bonded contact according to whether
 *                       protein residue already recorded as having a H-bond
 *                       to the DNA
 *
 *************************************************************************/
int
accept_reject_nbond(BONDS *first_bond, ATOMS *prot_atom)
{
/* Declare variables ====================================================*/ 
  int wanted;
  BONDS  *bond;
  ATOMS *atom[2];
   
  /* Initialise variables */
  wanted = TRUE;

  /* Get pointer to the first bond */
  bond = first_bond;

  /* Loop through bonds */
  while (bond != NULL && wanted == TRUE)
    {
      /* Get the two atoms involved in this bond */
      atom[0] = bond->bond_atom[0];
      atom[1] = bond->bond_atom[1];

      /* Check whether this residue is involved in an H-bond */
      if (bond->type == HBOND)
        {
          /* Check residue and chain-id */
          if ((prot_atom->resnum == atom[0]->resnum &&
              prot_atom->chain == atom[0]->chain) ||
              (prot_atom->resnum == atom[1]->resnum &&
              prot_atom->chain == atom[1]->chain))
            {
              /* Already have this residue, so don't want to store
                 non-bonded contacts */
              wanted = FALSE;
            }
        }
      
      /* Get pointer to next bond */
      bond = bond->next_bond;
   }
   
  /* Return wanted flag */
  return(wanted);
}



/**************************************************************************
 *
 * nbond_list - Function to record the non bonded contacts in a linked list
 *
 *************************************************************************/
int nbond_list(BONDS **first_bond, BONDS **last_bond, BONDS **first_nbond, 
               ATOMS *prot_atom, ATOMS *dna_atom, int bondtype,
               char input_line[LINE_LEN], int entity_flag)
{
/* Declare variables ====================================================*/ 
   double distance;
   BONDS  *bond;

/* Memory allocation ====================================================*/
   bond = (BONDS *) malloc(sizeof(BONDS));
   if (bond == NULL)
   {
      printf("*** ERROR - Unable to store bonding information\n");
      exit(1);
   }
   
/* Initialize or update pointers ========================================*/
   if (*first_bond == NULL)
   {
      *first_bond = bond;
      bond->prev_bond = NULL;
   }
   else
   {
      (*last_bond)->next_bond = bond;
      bond->prev_bond = *last_bond;
   }
   
   if (*first_nbond == NULL)
      *first_nbond = bond;
   
/* Record bonding information ===========================================*/
   bond->type = bondtype;
   bond->bond_atom[DONOR]    = prot_atom;
   bond->bond_atom[ACCEPTOR] = dna_atom;

/* Record bond distance =================================================*/
   sscanf(input_line+28, "%lf", &distance);
   bond->distance = distance;

/* Record entity of interacting groups ==================================*/
   bond->entity = entity_flag;

/* Update pointers ======================================================*/
   bond->next_bond = NULL;
   *last_bond = bond;

   return(0);
}



/**************************************************************************
 *
 * prot_dna_cbond - Function to get the covalent bonds between the 
 *                  protein and DNA
 *
 *************************************************************************/
int prot_dna_cbond(char input_line[LINE_LEN], ATOMS *first_prot[MAXARRAY], 
                   ATOMS *first_dna[MAXARRAY], BONDS **first_bond,
                   BONDS **last_bond,
                   BONDS **first_connect, int nprot, int ndna, PARAM *par)
{
/* Declare variables ====================================================*/ 
   int    dna_flag, prot_flag, water_flag;
   int    found_flag;
   int    bondtype_flag;
   double distance;
   ATOMS  *prot_res, *prot_atom; /* First atom of right res and right atom */
   ATOMS  *dna_res, *dna_atom;  /* First atom of right base and right atom */

/* Initialize variables =================================================*/
   dna_flag         = prot_flag = -1;
   found_flag       = TRUE;
   bondtype_flag    = COVALENT;
   prot_res         = prot_atom  = dna_res = dna_atom = NULL; 

/* Check which types of atoms the bond is between */
   get_atom_types(input_line, &dna_flag, &prot_flag, &water_flag, par); 

/* Record the N-bond if it is ===========================================*/
   if ((prot_flag == LEFT && dna_flag == RIGHT)
       || (prot_flag == RIGHT && dna_flag == LEFT))
   {
      /* Find the right protein residue */
      find_hbond_res(first_prot, nprot, &prot_res, &prot_atom, input_line, prot_flag);
      
      /* Find the right DNA base */
      find_hbond_res(first_dna, ndna, &dna_res, &dna_atom, input_line, dna_flag);
      
      /* Warning if HBPLUS atoms not found in PDB */
      hbond_warning(input_line, prot_res, prot_atom, dna_res, dna_atom,
/* v.1.1.2--> */
//                    &found_flag);
                    &found_flag,5);
/* <--v.1.1.2 */

      /* Proceed if found_flag = TRUE */
      if (found_flag == TRUE)
      {
         /* Pair up the hbonded atoms */
         cbond_list(&(*first_bond), &(*last_bond), &(*first_connect), prot_atom, dna_atom, 
                    bondtype_flag, input_line, PD);
      }
   }

/* Check the connect record =============================================*/
   check_connect_record(&(*first_connect), &(*last_bond), &distance);
   
   return(0);

}



/**************************************************************************
 *
 * cbond_list - Function to record the non bonded contacts in a linked list
 *
 *************************************************************************/
int cbond_list(BONDS **first_bond, BONDS **last_bond, BONDS **first_connect, 
               ATOMS *prot_atom, ATOMS *dna_atom, int bondtype,
               char input_line[LINE_LEN], int entity_flag)
{
/* Declare variables ====================================================*/ 
   double distance;
   BONDS  *bond;

/* Memory allocation ====================================================*/
   bond = (BONDS *) malloc(sizeof(BONDS));
   if (bond == NULL)
   {
      printf("*** ERROR - Unable to store bonding information\n");
      exit(1);
   }
   
/* Initialize or update pointers ========================================*/
   if (*first_bond == NULL)
   {
      *first_bond = bond;
      bond->prev_bond = NULL;
   }
   else
   {
      (*last_bond)->next_bond = bond;
      bond->prev_bond = *last_bond;
   }
   
   if (*first_connect == NULL)
      *first_connect = bond;
   
/* Record bonding information ===========================================*/
   bond->type = bondtype;
   bond->bond_atom[DONOR]    = prot_atom;
   bond->bond_atom[ACCEPTOR] = dna_atom;

/* Record bond distance =================================================*/
   sscanf(input_line+28, "%lf", &distance);
   bond->distance = distance;

/* Record entity of interacting groups ==================================*/
   bond->entity = entity_flag;

/* Update pointers ======================================================*/
   bond->next_bond = NULL;
   *last_bond = bond;

/* Return ===============================================================*/
   return(0);
}



/**************************************************************************
 *
 * read_connect_records - Function to read the PDB file and extract the 
 *                        covalent bonding information
 *
 *************************************************************************/
int
read_connect_records(ATOMS *first_prot[MAXARRAY], ATOMS *last_prot[MAXARRAY], 
                     ATOMS *first_dna[MAXARRAY], ATOMS *last_dna[MAXARRAY], 
                     BONDS **first_bond, BONDS **last_bond,
                     BONDS **first_connect, int nprot, int ndna)
{
/* Declare variables ====================================================*/
   int    ncon;
   int    numbers[MAXCON][2], store_dna[MAXCON], store_prot[MAXCON];
   int    indna, inprot;
   double distance;
   char   input_line[LINE_LEN];
   FILE   (*pdbfile);

/* Open PDB file ========================================================*/
   pdbfile = open_file(pdbfilename, "r");

/*=========================================================================
 * Read CONECT records
 *=======================================================================*/
   while (fgets(input_line, LINE_LEN, pdbfile) != NULL)
   {
      if (!strncmp(input_line,"CONECT",6))
      {
/* Loop through atom numbers on line ====================================*/
         get_connect_numbers(input_line, first_prot, last_prot, nprot,
                             first_dna, last_dna, ndna, numbers, store_prot,
                             store_dna, &inprot, &indna, &ncon);
         
/* Proceed if there is a bond between the protein and DNA ===============*/
         connect_list(first_prot, first_dna, &(*first_bond), &(*last_bond), 
                      &(*first_connect), numbers, store_prot, store_dna, 
                      inprot, indna, ncon);
 
/* Check the connect record =============================================*/
         check_connect_record(&(*first_bond), &(*last_bond), &distance);
      }
   }
   
/* Mark the first covalent bond =========================================*/
   mark_first_connect(&(*first_bond), &(*first_connect));
   
   /* test */
   /*bond = *first_bond;
   printf("bond type      donor      accepter\n");
   while (bond != NULL)
   {
      bond1 = bond->bond_atom[DONOR];
      bond2 = bond->bond_atom[ACCEPTOR];
      printf("%d     %s %c  %3d  %s      %s %c  %3d  %s\n",
             bond->type,
             bond1->resname, bond1->chain, bond1->resnum, bond1->name,
             bond2->resname, bond2->chain, bond2->resnum, bond2->name);
      bond = bond->next_bond;
   }*/
   
/* Return to main =======================================================*/
   return(0);
}



/**************************************************************************
 *
 * get_connect_numbers - Function to get the connectivity information from 
 *                       the PDB file
 *
 *************************************************************************/
int get_connect_numbers(char input_line[LINE_LEN],
                        ATOMS *first_prot[MAXARRAY], ATOMS *last_prot[MAXARRAY], int nprot,
                        ATOMS *first_dna[MAXARRAY], ATOMS *last_dna[MAXARRAY], int ndna,
                        int numbers[MAXCON][2], int store_prot[MAXCON], int store_dna[MAXCON],
                        int *inprot, int *indna, int *ncon)
{
/* Declare variables ===================================================*/
   int  i;
   int  line_pos;
   int  line_done;
   char atmnum[6];

/* Initialize variables ================================================*/
   *ncon     = 0;
   line_pos  = 6;
   line_done = FALSE;
   *inprot   = *indna = FALSE;
   for (i = 0; i < MAXCON; i++)
   {
      numbers[i][0]   = numbers[i][1]   = -1;
      store_dna[i] = store_dna[i] = store_prot[i] = store_prot[i] = FALSE;
   }

/* Loop through atom numbers on line ====================================*/
   while (line_done == FALSE && *ncon < MAXCON)
   {
      /* Extract this atom number and store */
      strncpy(atmnum, input_line+line_pos, 5);
      atmnum[5] = '\0';
      
      /* If have blank, then end of connections reached */
      if (!strncmp(atmnum,"     ",5))
         line_done = TRUE;
      
      /* Store the number as protein or DNA atom */
      else
      {
         numbers[*ncon][0] = atoi(atmnum);

         /* Check whether the atom belongs to the DNA */
         /* Loop through DNA strands */
         for (i = 0; i < ndna; i++)
         {
            if ((numbers[*ncon][0] >= first_dna[i]->number) &&
                (numbers[*ncon][0] <= last_dna[i]->number))
            {
               *indna = TRUE;
               store_dna[*ncon] = TRUE;

               /* Store chain number */
               numbers[*ncon][1] = i;
            }
         }
         
         for (i = 0; i < nprot; i++)
         {
            if ((numbers[*ncon][0] >= first_prot[i]->number) &&
                (numbers[*ncon][0] <= last_prot[i]->number))
            {
               *inprot = TRUE;
               store_prot[*ncon] = TRUE;

               /* Store chain number */
               numbers[*ncon][1] = i;
            }
         }
         
         /* Increment variables */
         (*ncon)++;
         line_pos = line_pos + 5;
      }
   }

/* Return to read_connect_records =======================================*/
   return(0);
}



/**************************************************************************
 *
 * connect_list - Function to record the covalent bond int the same list 
 *                as the H-bonds
 *
 *************************************************************************/
int
connect_list(ATOMS *first_prot[MAXARRAY], ATOMS *first_dna[MAXARRAY], 
             BONDS **first_bond, BONDS **last_bond, BONDS **first_connect, 
             int numbers[MAXCON][2], int store_prot[MAXCON],
             int store_dna[MAXCON], int inprot, int indna, int ncon)
{
/* Declare variables ====================================================*/
   int   i, j;
   int   have_bond;
   int   dna_chain;
   char  prot_chain;
   ATOMS *prot_atom, *dna_atom;
   ATOMS *bond1, *bond2;
   BONDS *bond;
   
/* Initialize variables =================================================*/
   j = 0;
   
/* Check if have connection between protein and DNA =====================*/   
   if (indna == TRUE && inprot == TRUE)
     {

/* Store bond details btw prot and DNA ==================================*/
      for (i = 1; i < ncon; i++)
      {
         if ((store_dna[0] == TRUE && store_prot[i] == TRUE)
             || (store_prot[0] == TRUE && store_dna[i] == TRUE))
         {

/* Check if this bond is already recorded ===============================*/
            /* Intialize variables */
            have_bond = FALSE;
            bond = *first_bond;

            while (bond != NULL)
            {
               bond1 = bond->bond_atom[DONOR];
               bond2 = bond->bond_atom[ACCEPTOR];
               
               if ((bond1->number == numbers[0][0] &&
                    bond2->number == numbers[i][0])
                   || (bond1->number == numbers[i][0] &&
                       bond2->number == numbers[0][0]))
                 have_bond = TRUE;

               bond = bond->next_bond;
            }

/* If not recorded, record details in bond ==============================*/
            if (have_bond == FALSE)
            {
               /* Memory allocation */
               bond = (BONDS *) malloc(sizeof(BONDS));
               if (bond == NULL)
               {
                  printf("*** ERROR - Unable to store bonding information\n");
                  exit(1);
               }
               
               /* Record bond type */
               bond->type = COVALENT;

               /* Find DNA atom */
               if (store_dna[0] == TRUE)
                 {
                   dna_atom = first_dna[numbers[0][1]];
                   j = 0;
                 }
               else if (store_dna[i] == TRUE)
                 {
                   dna_atom = first_dna[numbers[i][1]];
                   j = i;
                 }
               dna_chain = dna_atom->dna_chain;
               while (dna_atom != NULL && 
                      dna_atom->dna_chain == dna_chain)
               {
                  if (numbers[j][0] == dna_atom->number)
                     bond->bond_atom[ACCEPTOR] = dna_atom;
                  dna_atom = dna_atom->next_atom;
               }
               
               /* Find protein atom */
               if (store_prot[0] == TRUE)
                 {
                   prot_atom = first_prot[numbers[0][1]];
                   j = 0;
                 }
               else if (store_prot[i] == TRUE)
                 {
                   prot_atom = first_prot[numbers[i][1]];
                   j = i;
                 }
               prot_chain = prot_atom->chain;
               while ((prot_atom != NULL) && 
                      (prot_atom->chain == prot_chain))
               {
                  if (numbers[j][0] == prot_atom->number)
                     bond->bond_atom[DONOR] = prot_atom;
                  prot_atom = prot_atom->next_atom;
               }

               /* Update links */
               bond->prev_bond = *last_bond;
               bond->next_bond = NULL;
               (*last_bond)->next_bond = bond;
               *last_bond = bond; 
            }
            
         }            
      }            
   }
   
/* Return to read_connect_records =======================================*/
   return(0);
}



/**************************************************************************
 *
 * check_connect_record - Function to check the validity of the connect 
 *                        record 
 *
 *************************************************************************/
int
check_connect_record(BONDS **first_bond, BONDS **last_bond, double *distance)
{
/* Define variables =====================================================*/
   int    wanted;
   ATOMS  *bond1, *bond2;
   BONDS  *bond, *prev_bond, *next_bond;
   
/* Initialize variables =================================================*/
   bond = *first_bond;
   
/* Loop through the covalent bonds ======================================*/
   while (bond != NULL)
   {
      if (bond->type == COVALENT)
      {
         bond1 = bond->bond_atom[ACCEPTOR];
         bond2 = bond->bond_atom[DONOR];

         /* Check that neither is null */
         if (bond1 != NULL && bond2 != NULL)
           {

/* Calculate the distance between them ==================================*/
             *distance = (bond1->coord[0] - bond2->coord[0])
               * (bond1->coord[0] - bond2->coord[0])
                 + (bond1->coord[1] - bond2->coord[1])
                   * (bond1->coord[1] - bond2->coord[1])
                     + (bond1->coord[2] - bond2->coord[2])
                       * (bond1->coord[2] - bond2->coord[2]);
             *distance = sqrt(*distance);

/* Record distance of bond ==============================================*/
             bond->distance = *distance;

/* Check the bond length ================================================*/
             wanted = TRUE;
             if (*distance > 3.0)
               {
                 printf("*** Warning. Long bond from CONECT ");
                 printf("[%5d ->%5d] Length %7.2f",
                        bond1->number, bond2->number, *distance);
                 if (*distance > 4.0)
                   {
                     printf(" - discarded\n");
                     wanted = FALSE;
                   }
                 else
                   printf("\n");
               }
           }
         else
           wanted = FALSE;

         /* If link to be deleted, do so by making pointers by-pass it */
         if (wanted == FALSE)
         {
            /* Set up pointers */
            if (bond->prev_bond != NULL)
               prev_bond = bond->prev_bond;

            if (bond->next_bond != NULL)
               next_bond = bond->next_bond;
            
            /* If this is the first link, then set the next link as the
               first */
            if ( bond->prev_bond == NULL)
               *first_bond = bond->next_bond;
            
            /* Otherwise, make pointers bypass the current link */
            else
            {
               if (bond->prev_bond != NULL)
               {
                  prev_bond = bond->prev_bond;
                  prev_bond->next_bond = bond->next_bond;
               }
               if (bond->next_bond != NULL)
               {
                  next_bond = bond->next_bond;
                  next_bond->prev_bond = prev_bond;
               }
            }
            
            /* If this was the last in the list, then set the previous
               pointer as the new terminus */
            if (bond->next_bond == NULL)
               *last_bond = prev_bond;
         }
      }
      bond = bond->next_bond;
   }
   
/* Return to read_connect_records =======================================*/
   return(0);
}



/**************************************************************************
 *
 * mark_first_connect - Function to find the first covalent bond 
 *
 *************************************************************************/
int mark_first_connect(BONDS **first_bond, BONDS **first_connect)
{
/* Declare variables ====================================================*/
   int   found_bond;
   BONDS *bond;
   
/* Initialize variables =================================================*/
   found_bond = FALSE;
   bond = *first_bond;
   
/* Loop through bonds ===================================================*/
   while (bond != NULL)
   {
      if (bond->type == COVALENT && found_bond == FALSE)
      {
         *first_connect = bond;
         found_bond = TRUE;
      }
      bond = bond->next_bond;
   }
   
/* Return to read_connect_records =======================================*/
   return(0);
}



/**************************************************************************
 *
 * write_bonds - Function to write the bond to a .bond file
 *
 *************************************************************************/
int write_bonds(BONDS *first_bond, BONDS *first_nbond, BONDS *first_connect)
{
/* Declare variables ====================================================*/
   BONDS *hbond, *nbond, *cbond;
   ATOMS *acceptor, *donor;
   ATOMS *protein, *dna, *swap;

/* Initialize variables =================================================*/
   hbond = first_bond;
   nbond = first_nbond;
   cbond = first_connect;

/* Open bond file =======================================================*/
   bondfile = open_file(bondfilename, "w");
   
/* Initialize file ======================================================*/
   fprintf(bondfile, "NUCPLOT v.1.0  -  Bond file (%s)\n", bondfilename);
   fprintf(bondfile, "-----------------------------------------------\n\n");
   
/* Write H-bonds ========================================================*/
   fprintf(bondfile, "**** Hydrogen Bonds ***************************\n");
   fprintf(bondfile, "      Donor               Acceptor     Distance\n");

   while (hbond != NULL && hbond->type == HBOND)
   {
      donor    = hbond->bond_atom[DONOR];
      acceptor = hbond->bond_atom[ACCEPTOR];
      
      fprintf(bondfile, "%s %c %4d   %s    %s %c %4d   %s   %3.2f\n", 
              donor->resname, donor->chain, donor->resnum, donor->name,
              acceptor->resname, acceptor->chain, acceptor->resnum,
              acceptor->name,
              hbond->distance);
      
      hbond = hbond->next_bond;
   }

/* Write non bonded contacts ============================================*/
   fprintf(bondfile, "\n\n**** Non Bonded Contacts **********************\n");
   fprintf(bondfile, "     protein                DNA        Distance\n");
   
   while (nbond != NULL && nbond->type == NBOND)
   {
      donor    = nbond->bond_atom[DONOR];
      acceptor = nbond->bond_atom[ACCEPTOR];

      fprintf(bondfile, "%s %c %4d   %s    %s %c %4d   %s   %3.2f\n", 
              donor->resname, donor->chain, donor->resnum, donor->name,
              acceptor->resname, acceptor->chain, acceptor->resnum,
              acceptor->name,
              nbond->distance);

      nbond = nbond->next_bond;
   }

/* Write covalent bonds =================================================*/
   fprintf(bondfile, "\n\n**** Covalent Bonds ***************************\n");
   fprintf(bondfile, "     protein                DNA        Distance\n");

   while (cbond != NULL && cbond->type == COVALENT)
   {
      protein = cbond->bond_atom[DONOR];
      dna     = cbond->bond_atom[ACCEPTOR];

      /* If atom marked as protein is not a protein, then swap */
      if (protein->type != PROTEIN)
        {
          swap    = dna;
          dna     = protein;
          protein = swap;
        }
      
      fprintf(bondfile, "%s %c %4d   %s    %s %c %4d   %s   %3.2f\n", 
              protein->resname, protein->chain, protein->resnum, protein->name,
              dna->resname, dna->chain, dna->resnum, dna->name,
              hbond->distance);
      
      cbond = cbond->next_bond;
   }
   

/* Close bond file ======================================================*/
   fclose(bondfile);
   
/* Return to main =======================================================*/
   return(0);
}



/**************************************************************************
 *
 * paircatalog - Function to pair up the chains from the base-pairings
 *               previously calculated
 *               See lab notes 4
 *
 *************************************************************************/
void
paircatalog(ATOMS *first_dna[MAXARRAY], ATOMS *last_dna[MAXARRAY],
            ATOMS *second_strand[MEDIARRAY], 
            PAIRS *beg_pair[MEDIARRAY], PAIRS *end_pair[MEDIARRAY],
            int ndna, int *nchain, ATOMS *first_atom)
{
  int nalign, next_align, nchains;
  int flag;

  PAIRS *pair;
  struct node *first_node_ptr, *last_node_ptr;
   
  /* Initialize variables */
  flag = FALSE;
  pair = NULL;
  *nchain = 0;
  nalign = 0;
  next_align = 0;

  /* Create a linked list of nodes, one node per base */
  create_nodes(first_dna, ndna, &first_node_ptr, &last_node_ptr);

  /* Use the base-pairing information to link all the base-paired nodes
     together */
  link_nodes(first_node_ptr, first_atom);

  /* Assign the chains to the left- and right-hand sides of the page */
  assign_sides(first_dna, ndna, first_node_ptr);

  /* Sort the nodes according to their level numbers */
  sort_nodes(&first_node_ptr);

  /* Use the list of sorted nodes to create the list of pairs (matched
     as well as unmatched) to be plotted */
  create_pairs_list(first_node_ptr, beg_pair, end_pair, &nchains);
  *nchain = nchains;
}



/**************************************************************************
 *
 * create_nodes - Function to create a node for each base
 *
 *************************************************************************/
void
create_nodes(ATOMS *first_dna[MAXARRAY], int ndna,
             struct node **first_node_ptr, struct node **last_node_ptr)
{
  int colour, istrand;
  int first_base;
  int dna_chain;
  ATOMS *atom_ptr;
  struct node *node_ptr;

  /* Initialise node pointers */
  *last_node_ptr = NULL;
  *first_node_ptr = NULL;
  node_ptr = NULL;
  colour = 0;

  /* Loop through the separate DNA strands, creating a node for each
     base  */
  for (istrand = 0; istrand < ndna; istrand++)
    {
      /* Get first atom of this DNA strand */
      atom_ptr = first_dna[istrand];

      /* Initialise flag indicating this is the first base of the
         current strand */
      first_base = TRUE;
      dna_chain = atom_ptr->dna_chain;

      /* Loop through the residues for this chain */
      while (atom_ptr != NULL && atom_ptr->dna_chain == dna_chain)
        {
          /* Create a node for this base */
          node_ptr = (struct node*) malloc(sizeof(struct node));
          if (node_ptr == NULL)
            {
              printf("*** ERROR. Unable to allocate memory for struct");
              printf(" node\n");
              printf("***        Program nucplot terminated with error.\n");
              exit (1);
            }

          /* If this is the first node stored, then update
             pointer to head of linked list */
          if (*first_node_ptr == NULL)
            *first_node_ptr = node_ptr;
              
          /* Otherwise, make previous node point to the current one */
          else
            (*last_node_ptr)->next_node_ptr = node_ptr;

          /* Store the details for this node */
          node_ptr->atom_ptr = atom_ptr;
          node_ptr->level = (float) -99999.99;
          node_ptr->link_node_ptr[0] = NULL;
          node_ptr->link_node_ptr[1] = NULL;
          node_ptr->base_pair_node_ptr = NULL;
          node_ptr->paired = FALSE;
          node_ptr->linked = FALSE;
          node_ptr->nhbond = 0;
          node_ptr->colour = colour;
          node_ptr->first_base = first_base;
          node_ptr->last_base = FALSE;
          node_ptr->done = FALSE;
          node_ptr->next_node_ptr = NULL;

          /* Store link from base to node */
          atom_ptr->node_ptr = node_ptr;

          /* If this is not the first base of this chain, then connect
             previous base to this one */
          if (first_base == FALSE)
            {
              node_ptr->link_node_ptr[0] = *last_node_ptr;
              (*last_node_ptr)->link_node_ptr[1] = node_ptr;
            }
          first_base = FALSE;

          /* Save the current record's pointer */
          *last_node_ptr = node_ptr;

          /* Get pointer to atom representing the next base */
          atom_ptr = atom_ptr->next_res;
        }

      /* Mark last node as being the last base */
      if (node_ptr != NULL)
        node_ptr->last_base = TRUE;

      /* Shift chain colour by 1 */
      colour ++;
      if (colour > 3)
        colour = 0;
    }
}



/**************************************************************************
 *
 * link_nodes - Function to link the nodes representing the bases according
 *              to the base-pairing information
 *
 *************************************************************************/
void
link_nodes(struct node *first_node_ptr, ATOMS *first_atom)
{
  int resnum, other_resnum;
  int nhbond, max_hbond;
  int dna_chain, other_dna_chain;
  char chain, other_chain, resname[4], other_resname[4];
  ATOMS *atom_ptr, *other_atom_ptr, *scan_atom_ptr;
  struct node *node_ptr, *best_node_ptr, *other_node_ptr;

  /* Get pointer to the first node */
  node_ptr = first_node_ptr;

  /* Loop through all the nodes, pairing them up according to the base-
     pairing information */
  while (node_ptr != NULL)
    {
      /* If node not already linked, then check for base-pairing
         possibilities with all other nodes */
      if (node_ptr->linked == FALSE)
        {
          /* Get the base corresponding to this node */
          atom_ptr = node_ptr->atom_ptr;
          resnum = atom_ptr->resnum;
          strncpy(resname,atom_ptr->resname,3);
          resname[3] = '\0';
          chain = atom_ptr->chain;
          dna_chain = atom_ptr->dna_chain;

          /* Initialise counts of H-bond interactions */
          max_hbond = 0;
          best_node_ptr = NULL;

          /* Get pointer to the first node */
          other_node_ptr = first_node_ptr;

          /* Loop through all the other nodes counting numbers of H-bonds
             between each one and the current one */
          while (other_node_ptr != NULL)
            {
              /* Get the base corresponding to this node */
              other_atom_ptr = other_node_ptr->atom_ptr;
              other_resnum = other_atom_ptr->resnum;
              strncpy(other_resname,other_atom_ptr->resname,3);
              other_resname[3] = '\0';
              other_chain = other_atom_ptr->chain;
              other_dna_chain = other_atom_ptr->dna_chain;

              /* If node not already linked, then check for base-pairing
                 possibilities with all other nodes */
              if (other_node_ptr->linked == FALSE &&
                  other_node_ptr != node_ptr && other_dna_chain != dna_chain)
                {
                  /* Initialise count of H-bonds */
                  nhbond = 0;

                  /* Get pointer to the first atom */
                  scan_atom_ptr = first_atom;

                  /* Loop through all the atoms to count H-bonds between
                     these bases */
                  while (scan_atom_ptr != NULL)
                    {
                      /* If atom belongs to first base and is H-bonded
                         to second, then have a hit */
                      if (!strncmp(resname,scan_atom_ptr->resname,3) &&
                          resnum == scan_atom_ptr->resnum &&
                          chain ==  scan_atom_ptr->chain &&
                          !strncmp(other_resname,
                                   scan_atom_ptr->pair_resname,3) &&
                          other_resnum == scan_atom_ptr->pair_resnum &&
                          other_chain ==  scan_atom_ptr->pair_chain)
                        {
                          /* Increment count */
                          nhbond++;
                        }

                      /* Get pointer to next atom */
                      scan_atom_ptr = scan_atom_ptr->next_atom;
                    }

                  /* Check if number of H-bonds is greater than for
                     previous possible pairing */
                  if (nhbond > max_hbond)
                    {
                      max_hbond = nhbond;
                      best_node_ptr = other_node_ptr;
                    }

                  /* If same number of H-bonds, check for complementarity */
                  else if (nhbond == max_hbond && nhbond > 1)
                    {
                      /* Print warning */
                      printf("*** Warning. Ambiguous base-pairing:-\n");
                      printf("***          %s %d %c has %d H-bonds to ",
                             atom_ptr->resname,
                             atom_ptr->resnum,
                             atom_ptr->chain,nhbond);
                      printf("%s %d %c and %s %d %c\n",
                             best_node_ptr->atom_ptr->resname,
                             best_node_ptr->atom_ptr->resnum,
                             best_node_ptr->atom_ptr->chain,
                             other_node_ptr->atom_ptr->resname,
                             other_node_ptr->atom_ptr->resnum,
                             other_node_ptr->atom_ptr->chain);
                    }                  

                }
              /* Get pointer to the next node */
              other_node_ptr = other_node_ptr->next_node_ptr;
            }

          /* If have found a pairing, then link the nodes */
          if (max_hbond > 1)
            {
              /* Link the two nodes */
              node_ptr->base_pair_node_ptr = best_node_ptr;
              node_ptr->paired = TRUE;
              node_ptr->linked = TRUE;
              node_ptr->nhbond = max_hbond;
              best_node_ptr->base_pair_node_ptr = node_ptr;
              best_node_ptr->paired = TRUE;
              best_node_ptr->linked = TRUE;
              best_node_ptr->nhbond = max_hbond;
            }
        }

      /* Get pointer to the next node */
      node_ptr = node_ptr->next_node_ptr;
    }
}



/**************************************************************************
 *
 * assign_sides - Function to assign each chain to either the left or
 *                right-hand sides
 *
 *************************************************************************/
void
assign_sides(ATOMS *first_dna[MAXARRAY], int ndna,
             struct node *first_node_ptr)
{
  int nchain, next_chain, side, side_stack[MAXCHAINS];
  int all_chains[MAXCHAINS], chain_stack[MAXCHAINS];
  int dna_chain, ichain, other_dna_chain;
  int direction, ipos;
  float level, max_level;
  int added, done;

  float start_level;

  struct node *node_ptr, *other_node_ptr;


  /* Store all the different chains that we have */
  for (ichain = 0; ichain < ndna; ichain++)
    all_chains[ichain] = first_dna[ichain]->dna_chain;

  /* Initialise the side to the LEFT */
  side = LEFT;
  nchain = 0;
  next_chain = 0;

  /* Get starting chain and place on stack */
  chain_stack[nchain] = first_node_ptr->atom_ptr->dna_chain;
  side_stack[nchain] = side;
  nchain++;

  /* Initialise level-numbering at zero */
  level = -1.0;
  max_level = level;

  /* Loop while there are still chains on the stack */
  while (next_chain < nchain)
    {
      /* Get the next chain off the stack and its corresponding side */
      dna_chain = chain_stack[next_chain];
      side = side_stack[next_chain];

      /* Check whether this chain has any levels already assigned and
         determine */
      check_levels(first_node_ptr, dna_chain, max_level, &start_level,
                   &direction);

      /* Fill in all the level numbers for this chain */
      number_levels(first_node_ptr, dna_chain, start_level, direction,
                    &max_level);

      /* Get pointer to the first node */
      node_ptr = first_node_ptr;

      /* Loop through all the nodes, assigning them to the side on which
         this chain is to be plotted */
      while (node_ptr != NULL)
        {
          /* If the base is of the current chain, then update its side */
          if (node_ptr->atom_ptr->dna_chain == dna_chain)
            {
              node_ptr->side = side;

              /* If this is a linked node, then assign opposite side to its
                 partner chain */
              if (node_ptr->linked == TRUE)
                {
                  /* Get the linked node */
                  other_node_ptr = node_ptr->base_pair_node_ptr;

                  /* Give it the same level-number */
                  other_node_ptr->level = node_ptr->level;

                  /* Get its chain ID */
                  other_dna_chain = other_node_ptr->atom_ptr->dna_chain;

                  /* Loop through all the chains on the stack to see if
                     this one has been done */
                  done = FALSE;
                  for (ipos = 0; ipos < nchain; ipos++)
                    if (chain_stack[ipos] == other_dna_chain)
                      done = TRUE;

                  /* If chain hasn't been processed, then add to stack */
                  if (done == FALSE)
                    {
                      chain_stack[nchain] = other_dna_chain;
                      if (side == LEFT)
                        side_stack[nchain] = RIGHT;
                      else
                        side_stack[nchain] = LEFT;
                      nchain++;
                    }
                }
            }

          /* Get pointer to the next node */
          node_ptr = node_ptr->next_node_ptr;
        }

      /* Go to the next position in the stack */
      next_chain++;

      /* If have reached the end of the stack, check whether there are
         any chains that haven't been processed yet and add the first of
         these to the stack */
      if (next_chain == nchain)
        {
          /* Loop through each of the known chains */
          added = FALSE;
          for (ichain = 0; ichain < ndna && added == FALSE; ichain++)
            {
              /* Get this chain */
              dna_chain = all_chains[ichain];

              /* Check whether it has been processed */
              done = FALSE;
              for (ipos = 0; ipos < nchain; ipos++)
                if (chain_stack[ipos] == dna_chain)
                  done = TRUE;

              /* If chain hasn't been processed, then add to stack */
              if (done == FALSE)
                {
                  chain_stack[nchain] = dna_chain;
                  side_stack[nchain] = LEFT;
                  nchain++;

                  /* Set marker that one chain added */
                  added = TRUE;
                }
            }
        }
    }
}



/**************************************************************************
 *
 * check_levels - Function to check whether the given strand already has
 *                level numbers assigned and, if so, where they start and
 *                in which direction they run
 *
 *************************************************************************/
void
check_levels(struct node *first_node_ptr, int dna_chain, float max_level,
             float *start_level, int *direction)
{
  int unnumbered_levels;
  int done, first;
  float first_level, level;
  ATOMS *atom_ptr;
  struct node *node_ptr;

  /* Initialise variables */
  *direction = CANT_TELL;
  level = 0.0;
  unnumbered_levels = 0;
  first = TRUE;
  done = FALSE;

  /* Get pointer to the first node */
  node_ptr = first_node_ptr;

  /* Loop through all the nodes, numbering them according to their
     level in the final diagram */
  while (node_ptr != NULL && done == FALSE)
    {
      /* Get get the corresponding base */
      atom_ptr = node_ptr->atom_ptr;

      /* If this base belongs to the right chain, then process */
      if (atom_ptr->dna_chain == dna_chain)
        {
          /* If have a level number, then process */
          if (node_ptr->level > -99999.0)
            {
              /* Get the level number */
              level = node_ptr->level;

              /* If this is the first level-number encountered, then
                 store */
              if (first == TRUE)
                {
                  /* Store the level and change the flag */
                  first_level = level;
                  first = FALSE;
                }

              /* Otherwise, check whether numbers are rising or falling */
              else
                {
                  /* If level has increased, then levels are rising */
                  if (level > first_level)
                    *direction = RISING;
                  else
                    *direction = FALLING;

                  /* Have all the information we need */
                  done = TRUE;
                }

            }

          /* If this node has no level-number, increment count of levels
             without numbers */
          else if (first == TRUE)
            unnumbered_levels++;
        }

      /* Get pointer to the next node */
      node_ptr = node_ptr->next_node_ptr;
    }

  /* Determine what the starting level should be */

  /* If no level numbers found, then set starting level to the maximum */
  if (*direction == CANT_TELL)
    *start_level = max_level + 1;

  /* Otherwise, calculate starting level from information obtained */
  else
    {
      /* If direction is rising, extrapolate back to get starting level */
      if (*direction == RISING)
        *start_level = first_level - unnumbered_levels;

      /* Otherwise, extrapolate the other way */
      else
        *start_level = first_level + unnumbered_levels;
    }
}



/**************************************************************************
 *
 * number_levels - Function to number the node levels of the given chain
 *
 *************************************************************************/
void
number_levels(struct node *first_node_ptr, int dna_chain, float start_level,
              int direction, float *max_level)
{
  int unnumbered_nodes;
  int first, found;
  float level, level_step, next_level;
  ATOMS *atom_ptr, *other_atom_ptr;
  struct node *node_ptr, *other_node_ptr;

  /* Initialise variables */
  level = start_level;
  level_step = 1.0;
  first = TRUE;

  /* Get pointer to the first node */
  node_ptr = first_node_ptr;

  /* Loop through all the nodes */
  while (node_ptr != NULL)
    {
      /* Get get the corresponding base */
      atom_ptr = node_ptr->atom_ptr;

      /* If this base belongs to the right chain, then process */
      if (atom_ptr->dna_chain == dna_chain)
        {
          /* If have a level number, then store value */
          if (node_ptr->level > -99999.0)
            {
              /* Store current level_number */
              level = node_ptr->level;
              first = FALSE;

              /* Reset the level increment to 1.0 */
              level_step = 1.0;
            }

          /* Otherwise, have to assign a level-number */
          else
            {
              /* If this is an unnumbered node that is an insertion
                 (ie is not one of the ones leading up to the first
                 numbered one in the chain) then locate the next
                 numbered level in the chain */
              if (first == FALSE)
                {
                  /* Initialise variables for search for numbered level */
                  next_level = (float) -99999.9;
                  found = FALSE;
                  unnumbered_nodes = 1;

                  /* Get pointer to the next node */
                  other_node_ptr = node_ptr->next_node_ptr;

                  /* Search until get to another level-number in this
                     chain or come off the end */
                  while (other_node_ptr != NULL && found == FALSE)
                    {
                      /* Process only if node belongs to correct chain */
                      other_atom_ptr = other_node_ptr->atom_ptr;
                      if (other_atom_ptr->dna_chain == dna_chain)
                        {
                          /* If node has a level number then end here */
                          if (other_node_ptr->level > -99999.0)
                            {
                              next_level = other_node_ptr->level;
                              found = TRUE;
                            }

                          /* Otherwise, increment count of unnumbered
                             nodes */
                          else
                            unnumbered_nodes++;
                        }

                      /* Get pointer to the next node */
                      other_node_ptr = other_node_ptr->next_node_ptr;
                    }

                  /* If have found a number, then calculate the
                     level increment to be applied to the unnumbered
                     nodes in between */
                  if (found == TRUE)
                    {
                      level_step
                        = (float) (fabs(next_level - level) / (unnumbered_nodes + 1.0));
                    }

                  /* Otherwise, have come off the end of the chain, so
                     continue numbering in current direction with
                     standard increment of 1.0 */
                  else
                    {
                      /* Set increment to 1.0 */
                      level_step = 1.0;
                    }

                  /* Calculate next level from known direction */
                  if (direction == FALLING)
                    level = level - level_step;
                  else
                    level = level + level_step;
                }

              /* Assign a level-number to the current node */
              node_ptr->level = level;
              
              /* Increment or decrement accordingly */
              if (direction == FALLING)
                level = level - level_step;
              else
                level = level + level_step;

              /* Set flag as though first time, so level-numbering
                 will continue until next numbered node encountered */
              first = TRUE;
            }
        }

      /* If this is the highest level encountered so far, then store */
      if (node_ptr->level > *max_level)
        *max_level = node_ptr->level;

      /* Get pointer to the next node */
      node_ptr = node_ptr->next_node_ptr;
    }
}



/**************************************************************************
 *
 * sort_nodes - Function to sort the nodes according to their level numbers
 *
 *************************************************************************/
void
sort_nodes(struct node **first_node_ptr)
{
  int swapped;
  struct node *node_ptr, *last_node_ptr, *next_node_ptr;

  /* Initialise swap-pointer */
  swapped = TRUE;

  /* Loop until no more node-positions have been swapped */
  while (swapped == TRUE)
    {
      /* Reset swap-pointer */
      swapped = FALSE;

      /* Get pointer to the first node and the one after it */
      node_ptr = *first_node_ptr;
      if (node_ptr != NULL)
        next_node_ptr = node_ptr->next_node_ptr;
      else
        next_node_ptr = NULL;

      /* Initialise the last node pointer */
      last_node_ptr = NULL;

      /* Loop through all the nodes */
      while (node_ptr != NULL && next_node_ptr != NULL)
        {
          /* If level number drops, then swap the pointers */
          if (next_node_ptr->level < node_ptr->level)
            {
              /* Swap the pointers round */
              if (last_node_ptr == NULL)
                *first_node_ptr = next_node_ptr;
              else
                last_node_ptr->next_node_ptr = next_node_ptr;
              node_ptr->next_node_ptr = next_node_ptr->next_node_ptr;
              next_node_ptr->next_node_ptr = node_ptr;

              /* Set the current pointer as the one that's just come
                 into this position */
              node_ptr = next_node_ptr;

              /* Set flag showing pointers still out of order */
              swapped = TRUE;
            }

          /* Store the current node as the last one encountered */
          last_node_ptr = node_ptr;

          /* Get pointer to the next node and the one after */
          node_ptr = node_ptr->next_node_ptr;
          next_node_ptr = node_ptr->next_node_ptr;
        }
    }
}



/**************************************************************************
 *
 * create_pairs_list - Function to generate a list of pairs to be plotted
 *                     according to the sorted level numbers
 *
 *************************************************************************/
void
create_pairs_list(struct node *first_node_ptr, PAIRS *beg_pair[MEDIARRAY],
                  PAIRS *end_pair[MEDIARRAY], int *nchain)
{
  int i, ichain, j, k, l, side, other_side;
  int continues;
  float diff;
  ATOMS *atom_ptr, *other_atom_ptr;
  PAIRS *pair, *last_pair;
  struct node *node_ptr, *other_node_ptr, *next_node_ptr;

  /* Initialise pointers for pair lists */
  for (i = 0; i < MEDIARRAY; i++)
    {
      beg_pair[i] = NULL;
      end_pair[i] = NULL;
    }
  last_pair = NULL;
  *nchain = 0;
  ichain = 0;

  /* Get pointer to the first node and the one after it */
  node_ptr = first_node_ptr;

  /* Loop through all the nodes */
  while (node_ptr != NULL)
    {
      /* If this is an unpaired node, see if it can be paired up with
         another unpaired node */
      if (node_ptr->linked == FALSE)
        {
          /* Check if the next node has the same level number */
          next_node_ptr = node_ptr->next_node_ptr;

          /* Check the difference in level numbers */
          if (next_node_ptr != NULL)
            {
              diff = (float) fabs(node_ptr->level - next_node_ptr->level);

              /* If the levels are the same, then can combine the
                 nodes into a single (unconnected) base pair */
              if (diff < 0.00001)
                {
                  /* Link the two nodes */
                  node_ptr->base_pair_node_ptr = next_node_ptr;
                  next_node_ptr->base_pair_node_ptr = node_ptr;

                  /* Set link-flags */
                  node_ptr->linked = TRUE;
                  next_node_ptr->linked = TRUE;
                }
            }
        }

      /* Create a pair record for the given pair of nodes, unless already
         done */
      if (node_ptr->done == FALSE)
        {
          /* Get the linked base-pair (if any) */
          other_node_ptr = node_ptr->base_pair_node_ptr;

          /* Create a pair record */
          pair = (PAIRS *) malloc (sizeof(PAIRS));
          if (pair == NULL)
            {
              printf("*** ERROR - Unable to store pairing information\n");
              printf("    Cannot allocate enough memory for pair");
              printf("%c %2d--%2d %c\n",
                     node_ptr->atom_ptr->chain,
                     node_ptr->atom_ptr->resnum,
                     node_ptr->atom_ptr->pair_resnum,
                     node_ptr->atom_ptr->pair_chain);
              exit(1);
            }

          /* Initialise various fields */
          for (i = 0; i < 2; i++)
            {
              pair->strand[i] = NULL;
              pair->first_base[i] = FALSE;
              pair->last_base[i] = FALSE;
            }
          pair->paired = FALSE;

          /* Get this node's side */
          side = node_ptr->side;
          pair->first_base[side] = node_ptr->first_base;
          pair->last_base[side] = node_ptr->last_base;
          pair->colour[side] = node_ptr->colour;

          /* Store pairing information */
          atom_ptr = node_ptr->atom_ptr;
          pair->strand[side] = atom_ptr;

          /* Point from atom-record to the pair record */
          atom_ptr->pair = pair;

          /* Update node as done */
          node_ptr->done = TRUE;

          /* If have a paired base, then get its details */
          if (other_node_ptr != NULL)
            {
              /* Store other atom's details */
              other_atom_ptr = other_node_ptr->atom_ptr;
              other_side = other_node_ptr->side;
              pair->strand[other_side] = other_atom_ptr;
              pair->first_base[other_side] = other_node_ptr->first_base;
              pair->last_base[other_side] = other_node_ptr->last_base;
              pair->colour[other_side] = other_node_ptr->colour;

              /* Point from atom-record to the pair record */
              other_atom_ptr->pair = pair;

              /* If the bases are paired, then set paired flag */
              if (node_ptr->paired == TRUE)
                pair->paired = TRUE;

              /* Update other node as done */
              other_node_ptr->done = TRUE;
            }

          /* Initialize other information */
          pair->x_pos[ONE] = pair->x_pos[TWO] = 0;
          pair->y_pos[ONE] = pair->y_pos[TWO] = 0;
          for (j = 0; j < 2; j++)
            {
              pair->nbond[j] = pair->npbond[j] = 0;
              for (k = 0; k < MAXARRAY; k++)
                {
                  pair->type[j][k] = pair->ptype[j][k] = 0;
                  pair->nbridge[j][k] = pair->npbridge[j][k] = 0;
                  pair->bridge_flag[j][k] = pair->pbridge_flag[j][k]
                    = FALSE;
                  pair->bondto[j][k] = NULL;
                  pair->pbondto[j][k] = NULL;
                  pair->combine[j][k] = NULL;
                  pair->pcombine[j][k] = NULL;
                  pair->bondfrom[j][k] = NULL;
                  pair->pbondfrom[j][k] = NULL;
                  pair->ndups[j][k] = pair->npdups[j][k] = 1;
                  for (l = 0; l < MAXARRAY; l++)
                    {
                      pair->bridgeto[j][k][l] = NULL;
                      pair->pbridgeto[j][k][l] = NULL;
                      pair->bridge_combine[j][k][l] = NULL;
                      pair->pbridge_combine[j][k][l] = NULL;
                    }
                }
            }

          /* If this pair is not connected to the previous one, then
             have a break and need to start a new chain */
          continues = FALSE;
          if (last_pair != NULL)
            {
              /* Check continuity of chains in either strand */
              for (i = 0; i < 2; i++)
                {
                  if (last_pair->strand[i] != NULL &&
                      pair->strand[i] != NULL)
                    {
                      /* Get the two atom pointers */
                      atom_ptr = pair->strand[i];
                     other_atom_ptr = last_pair->strand[i];

                      /* If the chain doesn't change, have continuity */
                      if (atom_ptr->chain == other_atom_ptr->chain)
                        continues = TRUE;
                    }
                }
            }

          /* If chain continuity broken, start a new chain */
          if (continues == FALSE)
            {
              /* Start a new chain */
              ichain = *nchain;
              beg_pair[ichain] = pair;
              (*nchain)++;
              last_pair = NULL;
            }

          /* Save this position as the end of the current chain */
          end_pair[ichain] = pair;

          /* Update pointers */
          pair->next_pair = NULL;
          pair->prev_pair = last_pair;

          /* Get last pair to point to this pair */
          if (last_pair != NULL)
            last_pair->next_pair = pair;

          /* Store the current pair */
          last_pair = pair;
        }

      /* Get pointer to the next node */
      node_ptr = node_ptr->next_node_ptr;
    }
}



/**************************************************************************
 *
 * num_base - Function to number each base according to the numbering 
 *            system
 *
 *************************************************************************/
int num_base(int num_system, ATOMS *first_dna[MAXARRAY], int ndna)
{
/* Declare variables ====================================================*/
   int   i;
   int dna_chain;
   int   num;                          /* Counter for numbering residues */
   ATOMS *atom, *last_atom, *pair;
   
/* Initialize variables =================================================*/
   i    = 0;
   atom = NULL;

/*=========================================================================
 * Assign numbering to the P atoms
 *=======================================================================*/

/* Default numbering ====================================================*/
   if (num_system == 0)
     {
       /* Loop through all the DNA atoms */
       atom = first_dna[0];
       while (atom != NULL)
         {
           /* Number only the P atoms */
           if (atom->name[1] == 'P')
             
             /* numbering  number is the same as the resnum */
             atom->numbering = atom->resnum;
           atom = atom->next_atom;
         }
     }

/* Antiparallel numbering system ========================================*/
/*   if ((num_system == 1) || (num_system == 2))
   {
      for (i = 0; i < ndna; i++)
      {
         atom = first_dna[i];
         while ((atom != NULL) && (atom->chain == first_dna[i]->chain))
         {
            if (atom->resnum != last_atom->resnum)
               num++;
            
            atom->numbering = num;
            
            last_atom = atom;
            atom = atom->next_atom;
         }
      }
   }
*/ 

/* Parallel numbering system ============================================*/
   if ((num_system == 4) || (num_system == 5))
   {
      for (i = 0; i < ndna; i++)
      {
         num = 1;
         atom = first_dna[i];
         last_atom = atom;
         dna_chain = atom->dna_chain;

         /* Loop through the atoms of the strand */
         while (atom != NULL && atom->dna_chain == dna_chain)
         {
            /* Increment number count if new base */
            if (atom->resnum != last_atom->resnum)
               num++;
            
            /* Number the P atom only if not previously numbered */
            if (atom->numbering == -1 && atom->name[1] == 'P')
            {
               /* Number the first and second strands identically */
               atom->numbering = num;
               pair = atom->pair_dna;
               pair->numbering = num;
            }
            last_atom = atom; 
            atom = atom->next_atom;
         }
      }
   }   



/*=========================================================================
 * Print base pairing
 *=======================================================================*/
/*   for (i = 0; i < ndna; i++)
   {
      atom = first_dna[i];
      printf("\n\n");

      while ((atom != NULL) && 
                (atom->chain == first_dna[i]->chain))
      {
         if (atom->pair_dna != NULL)
         {
            pair = atom->pair_dna;
            printf("%4d   (%c)%s  --%s  (%c) %4d\n", 
                   atom->numbering, atom->chain, atom->resname,
                   pair->resname, pair->chain, pair->numbering);
         }
         
         atom = atom->next_res;
      }
   }
*/
/* Return to main =======================================================*/
   return(0);
}



/***********************************************************************
 *
 * draw_picture - Function to draw the postscript picture
 *
 **********************************************************************/
int draw_picture(PAIRS *beg_pair[MEDIARRAY], PAIRS *end_pair[MEDIARRAY],
                 int nchain, BONDS *first_bond, BONDS *first_nbond,
                 BONDS *first_connect, PARAM *par)
{
/* Declare variables =================================================*/
   int    ichain;
   int    npages, npairs, page_no, total_pages;
   int    max_nbond, max_npbond, max_nbridge, max_npbridge;
   int    ds_flag[2];  /* elements for strand 1 and 2 - TRUE if exist */ 
   char   page_number[10];
   BONDS  *first_pw_bond;
   PAIRS  *pair, *last_pair, *first_page_pair, *last_page_pair;
   SIZES  sz;
   SIZES  *size;
   
/* Initialize variables ==============================================*/
   max_nbond = max_npbond = max_nbridge = max_npbridge = 0;
   size      = &sz;
   page_no    = 1;
   total_pages = 0;

   /* Calculate total number of pages for the whole plot */
   for (ichain = 0; ichain < nchain; ichain++)
     {
       npages = calc_no_of_pages(beg_pair[ichain], par->numpairs,
                                 &npairs);
       total_pages = total_pages + npages;
     }
   size->total_pages = total_pages;

/* Loop through DNA chains ===========================================*/
   for (ichain = 0; ichain < nchain; ichain++)
   {
      /* Initialize */
      if (par->write_to_nucplot_ps == FALSE)
        page_no    = 1;
      pair       = beg_pair[ichain];
      last_pair  = end_pair[ichain];
      ds_flag[0] = ds_flag[1] = FALSE;

      /* Find out whether the DNA is single stranded or double stranded */
      ss_or_ds(pair, ds_flag);

      /* Find the right DNA residues for bonds */
      find_bonded_base(pair, last_pair, first_bond, par, &max_nbond,
                       &max_npbond);

      /* Find first protein-water bond in list */
      find_first_pw(first_bond, &first_pw_bond);

      /* Sort interactions in order for drawing */
      sort_interactions(pair, last_pair, first_pw_bond, par, 
                        &max_nbridge, &max_npbridge);

      /* Open new psfile for each DNA chain */
      if ((par->write_to_nucplot_ps == TRUE && ichain == 0) ||
          par->write_to_nucplot_ps == FALSE)
        open_psfile(beg_pair, ichain, par);
      
      /* Scale the drawing accordingly */
      define_object_sizes(beg_pair, ichain, size, par, max_nbond,
                          max_npbond, max_nbridge, max_npbridge, ds_flag);
      
      /* Opening lines of postscript file */
/* v.1.1.1--> */
      if ((par->write_to_nucplot_ps == TRUE && ichain == 0) ||
          par->write_to_nucplot_ps == FALSE)
/* <--v.1.1.1 */
        psopen_(size, par);

/*======================================================================
 * Text version picture
 *====================================================================*/
      while (pair != NULL)
        {
          /* Record x-y coordinates of the positions of bases */
          base_pos(pair, size, par, ds_flag);
          
          /* Write starting lines for page */
          pspage_(page_no, size, par);
            
          /* Draw the base-pairs and bond between them */
          first_page_pair = pair;
          pair = text_nucleic_acids(beg_pair, end_pair, pair, ichain, size,
                                    par, ds_flag);
          last_page_pair = pair;

          /* Draw all the residues and waters interacting directly
             with the DNA */
          text_interactions(beg_pair[ichain], end_pair[ichain],
                            first_page_pair, last_page_pair,
                            size, par, max_nbond, max_npbond);
          
          /* Draw the nucleic acid backbone, icluding phosphates
             and sugars */
          text_backbone(first_page_pair, last_page_pair, size, par);

          /* Write out the page-number */
          sprintf(page_number,"Page %d",page_no);
          if (par->landscape == TRUE)
            psrtxt_(LPAGE_POSX,LPAGE_POSY,(int) PAGE_TXT,page_number);
          else
            psrtxt_(PAGE_POSX,PAGE_POSY,(int) PAGE_TXT,page_number);
              
          /* Close page */
          psendp_();

          /* Increment page count */
          page_no++;
        }

      /* Close file */
      if (par->write_to_nucplot_ps == FALSE)
        psclos_();
   }

   /* Close the nucplot.ps file, if there is one */
   if (par->write_to_nucplot_ps == TRUE)
     psclos_();

/* Return to main ====================================================*/
   return(0);
}



/***********************************************************************
 *
 * ss_or_ds - Function to find out whether the DNA is double or single 
 *            stranded 
 *
 **********************************************************************/
int ss_or_ds(PAIRS *first_pair, int ds_flag[2])
{
/* Declare variables =================================================*/
   int   strand_flag[2];
   PAIRS *pair;
   
/* Initialize variables ==============================================*/
   strand_flag[ONE] = strand_flag[TWO] = FALSE;
   pair             = first_pair;

/* Loop throught DNA =================================================*/
   while (pair != NULL)
   {
      if (pair->strand[ONE] != NULL)
         ds_flag[ONE] = TRUE;
      if (pair->strand[TWO] != NULL)
         ds_flag[TWO] = TRUE;
      
      pair = pair->next_pair;
   }

/* Return to draw_picture ============================================*/
   return(0);
}



/***********************************************************************
 *
 * find_first_pw - Function to find the first protein-water bond in the 
 *                 bond list
 *
 **********************************************************************/
int find_first_pw(BONDS *first_bond, BONDS **first_pw_bond)
{
/* Declare variables =================================================*/
   int   found_flag;
   ATOMS *atom[2];
   BONDS *bond;
   
/* Initialize variables ==============================================*/
   found_flag = FALSE;
   bond       = first_bond;
   
/* Loop through bonds ================================================*/
   while (bond != NULL)
   {
      atom[0] = bond->bond_atom[0];
      atom[1] = bond->bond_atom[1];
      
/* Check to see if protein-water bond ================================*/
      if (found_flag == FALSE)
         if (bond->entity == PW)
         {
            *first_pw_bond = bond;
            found_flag = TRUE;
         }
      bond = bond->next_bond;
   }
   
/* Return to draw_pictures ===========================================*/
   return(0);
}



/***********************************************************************
 *
 * sort_interactions - Function to sort interactions in order for pairing 
 *
 **********************************************************************/
int sort_interactions(PAIRS *first_pair, PAIRS *last_pair, BONDS *first_bond, 
                      PARAM *par, int *max_nbridge, int *max_npbridge)
{
/* Declare variables =================================================*/
   int   i;
   int   nprot_int, nwater_int;
   PAIRS *pair;

/* Initialize variables ==============================================*/
   pair = first_pair;

/* Loop through pairs ================================================*/
   while (pair != NULL)
   {

/* Loop through strands ==============================================*/
      for (i = 0; i < 2; i++)
      {

/* Non-phosphate bonds ===============================================*/
         if (pair->nbond[i] != 0)
         { 
            nprot_int = nwater_int = 0;

            /* Place interactions in resnum order */
            interact_resnum_order(pair, i, &nprot_int, &nwater_int, NONPHOS);
            
            /* Separate out bridged and non-bridged waters */
            sep_bridged_waters(pair, i, nprot_int, nwater_int, NONPHOS);

            /* Record the bridging bonds PAIRS */
            find_bridge_bonds(pair, first_bond, i, nprot_int, 
                              nwater_int, &(*max_nbridge), &(*max_npbridge),
                              NONPHOS);

            /* Place bridging int in resnum order */
            bridgeint_resnum_order(pair, i, nprot_int, nwater_int, NONPHOS);

            /* Remove duplicate interactions */
            rm_dup_int(pair, i, NONPHOS, par->nonbond_by_residue);
         }
         
/* Phosphate bonds ===================================================*/
         if (pair->npbond[i] != 0)
         {
            nprot_int = nwater_int = 0;

            /* Place interactions in resnum order */
            interact_resnum_order(pair, i, &nprot_int, &nwater_int, PHOS);

            /* Separate out bridged and non-bridged waters */
            sep_bridged_waters(pair, i, nprot_int, nwater_int, PHOS);

            /* Record the bridging bonds PAIRS */
            find_bridge_bonds(pair, first_bond, i, nprot_int, 
                              nwater_int, &(*max_nbridge), &(*max_npbridge),
                              PHOS);

            /* Place bridging int in resnum order */
            bridgeint_resnum_order(pair, i, nprot_int, nwater_int, PHOS);

            /* Remove duplicate interactions */
            rm_dup_int(pair, i, PHOS, par->nonbond_by_residue);
         }
      }
      pair = pair->next_pair;
   }

/* Return to draw_pictures ===========================================*/
   return(0);
}



/***********************************************************************
 *
 * interact_resnum_order - Function to place interacting groups in 
 *                         resnum order
 *
 **********************************************************************/
int 
interact_resnum_order(PAIRS *pair, int i, int *nprot_int,
                      int *nwater_int, int phos_flag)
{
/* Declare variables =================================================*/
   int   j;
   int   nbonds;
   int   swap_flag;
   ATOMS *interaction[MAXARRAY];
   ATOMS *swap_int;

/* Initialize variables ==============================================*/
   if (phos_flag == NONPHOS)
   {
      nbonds = pair->nbond[i];
      if (nbonds > MAXARRAY)
        nbonds = MAXARRAY;
      for (j = 0; j < nbonds; j++)
         interaction[j] = pair->bondto[i][j];
   }
   else if (phos_flag == PHOS)
   {
      nbonds = pair->npbond[i];
      if (nbonds > MAXARRAY)
        nbonds = MAXARRAY;
      for (j = 0; j < nbonds; j++)
         interaction[j] = pair->pbondto[i][j];
   }

/* Count up no. of protein interactions and water interactions =======*/
   for (j = 0; j < nbonds; j++)
   {
      if (strncmp(interaction[j]->resname, "HOH", 3))
         (*nprot_int)++;
      else if (!strncmp(interaction[j]->resname, "HOH", 3))
         (*nwater_int)++;
   }

/* Sort in resnum order ==============================================*/
   /* For protein interactions */
   swap_flag = TRUE;
   while (swap_flag == TRUE)
   {
      swap_flag = FALSE;
      for (j = 1; j < *nprot_int; j++)
      {
         if (interaction[j]->resnum < interaction[j-1]->resnum)
         {
            swap_int         = interaction[j];
            interaction[j]   = interaction[j-1];
            interaction[j-1] = swap_int;
            
            /* Swap these two interactions */
            swap_interactions(pair, i, j, j - 1, phos_flag);

            /* Set swap flag */
            swap_flag = TRUE;
         }
      }
   }

   /* For water interactions */
   swap_flag = TRUE;
   while (swap_flag == TRUE)
   {
      swap_flag = FALSE;
      for (j = (*nprot_int)+1; j < nbonds; j++)
      {
         if (interaction[j]->resnum < interaction[j-1]->resnum)
         {
            swap_int         = interaction[j];
            interaction[j]   = interaction[j-1];
            interaction[j-1] = swap_int;
               
            /* Swap these two interactions */
            swap_interactions(pair, i, j, j - 1, phos_flag);

            /* Set swap flag */
            swap_flag = TRUE;
         }
      }
   }

/* Return to sort_interactions =======================================*/
   return(0);
}



/***********************************************************************
 *
 * swap_interactions - Function to swap the given pair of interactions
 *
 **********************************************************************/
void
swap_interactions(PAIRS *pair, int i, int j, int k, int phos_flag)
{
  int m;
  int swap_bridge_flag, swap_nbridge, swap_type;
  ATOMS *swap_bondfrom, *swap_bondto, *swap_bridgeto;

  /* Perform the swap according to which set of interactions require
     swapping */
  if (phos_flag == NONPHOS)
    {
      /* Swap the pointers to the interacting atom on protein */
      swap_bondto        = pair->bondto[i][j];
      pair->bondto[i][j] = pair->bondto[i][k];
      pair->bondto[i][k] = swap_bondto;

      /* Swap the pointers to the combine-atoms */
      swap_bondto        = pair->combine[i][j];
      pair->combine[i][j] = pair->combine[i][k];
      pair->combine[i][k] = swap_bondto;

      /* Swap the pointers to the interacting atom on base */
      swap_bondfrom        = pair->bondfrom[i][j];
      pair->bondfrom[i][j] = pair->bondfrom[i][k];
      pair->bondfrom[i][k] = swap_bondfrom;

      /* Swap the associated bridge flags */
      swap_bridge_flag        = pair->bridge_flag[i][j];
      pair->bridge_flag[i][j] = pair->bridge_flag[i][k];
      pair->bridge_flag[i][k] = swap_bridge_flag;

      /* Swap number of bridging interactions */
      swap_nbridge        = pair->nbridge[i][j];
      pair->nbridge[i][j] = pair->nbridge[i][k];
      pair->nbridge[i][k] = swap_nbridge;

      /* Swap bond types */
      swap_type        = pair->type[i][j];
      pair->type[i][j] = pair->type[i][k];
      pair->type[i][k] = swap_type;

      /* Swap all the bridged interactions */
      for (m = 0; m < MAXARRAY; m++)
        {
          swap_bridgeto           = pair->bridgeto[i][j][m];
          pair->bridgeto[i][j][m] = pair->bridgeto[i][k][m];
          pair->bridgeto[i][k][m] = swap_bridgeto;

          /* Swap combined bridge interactions */
          swap_bridgeto           = pair->bridge_combine[i][j][m];
          pair->bridge_combine[i][j][m] = pair->bridge_combine[i][k][m];
          pair->bridge_combine[i][k][m] = swap_bridgeto;
        }
    }
  else if (phos_flag == PHOS)
    {
      /* Swap the pointers to the interacting atom */
      swap_bondto         = pair->pbondto[i][j];
      pair->pbondto[i][j] = pair->pbondto[i][k];
      pair->pbondto[i][k] = swap_bondto;

      /* Swap the pointers to the combine-atoms */
      swap_bondto        = pair->pcombine[i][j];
      pair->pcombine[i][j] = pair->pcombine[i][k];
      pair->pcombine[i][k] = swap_bondto;

      /* Swap the pointers to the interacting atom on base */
      swap_bondfrom         = pair->pbondfrom[i][j];
      pair->pbondfrom[i][j] = pair->pbondfrom[i][k];
      pair->pbondfrom[i][k] = swap_bondfrom;

      /* Swap the associated bridge flags */
      swap_bridge_flag         = pair->pbridge_flag[i][j];
      pair->pbridge_flag[i][j] = pair->pbridge_flag[i][k];
      pair->pbridge_flag[i][k] = swap_bridge_flag;

      /* Swap number of bridging interactions */
      swap_nbridge         = pair->npbridge[i][j];
      pair->npbridge[i][j] = pair->npbridge[i][k];
      pair->npbridge[i][k] = swap_nbridge;

      /* Swap bond types */
      swap_type         = pair->ptype[i][j];
      pair->ptype[i][j] = pair->ptype[i][k];
      pair->ptype[i][k] = swap_type;

      /* Swap all the bridged interactions */
      for (m = 0; m < MAXARRAY; m++)
        {
          swap_bridgeto            = pair->pbridgeto[i][j][m];
          pair->pbridgeto[i][j][m] = pair->pbridgeto[i][k][m];
          pair->pbridgeto[i][k][m] = swap_bridgeto;

          /* Swap combined bridge interactions */
          swap_bridgeto           = pair->pbridge_combine[i][j][m];
          pair->pbridge_combine[i][j][m] = pair->pbridge_combine[i][k][m];
          pair->pbridge_combine[i][k][m] = swap_bridgeto;
        }
    }
}



/***********************************************************************
 *
 * sep_bridged_waters - Function to place the bridged waters at the 
 *                      bottom of the pile
 *
 **********************************************************************/
int sep_bridged_waters(PAIRS *pair, int i, int nprot_int, int nwater_int,
                       int phos_flag)
{
/* Declare variables =================================================*/
   int   j, k, l;
   int   nbonds;
   ATOMS *interaction[MAXARRAY];
   ATOMS *swap_int;

/* Initialize variables ==============================================*/
   if (phos_flag == NONPHOS)
   {
      nbonds = pair->nbond[i];
      for (j = 0; j < nbonds; j++)
         interaction[j] = pair->bondto[i][j];
   }
   else if (phos_flag == PHOS)
   {
      nbonds = pair->npbond[i];
      for (j = 0; j < nbonds; j++)
         interaction[j] = pair->pbondto[i][j];
   }

/* Check bridge flag is working */
   /*for (j = nprot_int; j < nbonds; j++)
   {
      if (interaction[j]->bridge_flag == TRUE)
         printf("Found bridging water: HOH %d\n", interaction[j]->resnum);
   }*/
   
/* Stack bridging waters at the bottom of the pile ===================*/
/* v.1.1.2--> */
//   j = l = nprot_int;
//   while (j < nbonds)
//   {
//      if (interaction[l]->bridge_flag == TRUE)
//      {
//         /* Place at bottom of stack */
//         swap_int = interaction[l];
//         interaction[nbonds] = swap_int;
//
//         /* Swap with position off the end of the stack */
//         swap_interactions(pair, i, l, nbonds, phos_flag);
//         
//         /* Move everything up one */
//         move_up_one(pair, i, l, nbonds+1, phos_flag);
//
//         /* Move corresponding record of interaction */
//         for (k = l; k < nbonds+1; k++)
//           interaction[k] = interaction[k+1];
//      }
//      else if (interaction[nprot_int]->bridge_flag == FALSE)
//         l++;
//      j++;
//   }
   for (j = 0; j < nbonds; j++)
     {
       if (interaction[j]->bridge_flag == TRUE)
         {
           /* Place at bottom of stack */
           swap_int = interaction[j];
           interaction[nbonds] = swap_int;

           /* Swap with position off the end of the stack */
           swap_interactions(pair, i, j, nbonds, phos_flag);
         
           /* Move everything up one */
           move_up_one(pair, i, j, nbonds+1, phos_flag);

           /* Move corresponding record of interaction */
           for (k = j; k < nbonds+1; k++)
             interaction[k] = interaction[k+1];
         }
     }
/* <--v.1.1.2 */
   
/* Return to sort_interactions =======================================*/
   return(0);
}



/***********************************************************************
 *
 * find_bridge_bonds - Function to find the protein residues to which 
 *                     the briding waters are bonded
 *
 **********************************************************************/
int find_bridge_bonds(PAIRS *pair, BONDS *first_bond, int i, int nprot_int, 
                      int nwater_int, int *max_nbridge, int *max_npbridge, 
                      int phos_flag)
{
/* Declare variables =================================================*/
   int   j, k, l;
   int   nbonds;
   int   total_bonds;
   int   met_flag;
   ATOMS *interaction[MAXARRAY];
   ATOMS *atom[2];
   ATOMS *water_atom;
   ATOMS *prot_atom; 
   ATOMS *check_atom;
   BONDS *bond;
   
/* Initialize variables ==============================================*/
   total_bonds = 0;
   if (phos_flag == NONPHOS)
   {
      nbonds = pair->nbond[i];
      for (j = 0; j < nbonds; j++)
         interaction[j] = pair->bondto[i][j];
   }
   else if (phos_flag == PHOS)
   {
      nbonds = pair->npbond[i];
      for (j = 0; j < nbonds; j++)
         interaction[j] = pair->pbondto[i][j];
   }

/* Loop through the waters bonded to the base ========================*/
   for (j = nprot_int; j < nbonds; j++)
   {
      if (interaction[j]->bridge_flag == TRUE)
      {
         bond = first_bond;
         if (phos_flag == NONPHOS)
            pair->bridge_flag[i][j] = TRUE;
         else if (phos_flag == PHOS)
            pair->pbridge_flag[i][j] = TRUE;

/* Loop through protein water bonds ==================================*/
         while (bond != NULL)
         {
            atom[0] = bond->bond_atom[0];
            atom[1] = bond->bond_atom[1];
            l       = 1;

            for (k = 0; k < 2; k++)
               if (!strncmp(atom[k]->resname, "HOH", 3))
               {
                  l = l - k;
                  water_atom = atom[k];
                  prot_atom  = atom[l];
               }
            
/* Compare bond list atom with pair list atom - are they the same water? =*/
            if (interaction[j]->resnum == water_atom->resnum)
            {
               if (phos_flag == NONPHOS)
               {
                  /* Check that the prot atom has not been recorded before */
                  met_flag = FALSE;
                  for (k = 0; k < pair->nbridge[i][j]; k++)
                  {
                     check_atom = pair->bridgeto[i][j][k];
                     if (prot_atom->resnum == check_atom->resnum
                         && !strncmp(prot_atom->name, check_atom->name, 4))
                        met_flag = TRUE;
                  }
                  
                  if (met_flag == FALSE)
                  {
                     pair->bridgeto[i][j][pair->nbridge[i][j]] = prot_atom;
                     pair->nbridge[i][j]++;
                     total_bonds++;
                  }
               }
               else if (phos_flag == PHOS)
               {        
                  met_flag = FALSE;
                  for (k = 0; k < pair->npbridge[i][j]; k++)
                  {
                     check_atom = pair->pbridgeto[i][j][k];
                     if (prot_atom->resnum == check_atom->resnum
                         && !strncmp(prot_atom->name, check_atom->name, 4))
                        met_flag = TRUE;
                  }

                  if (met_flag == FALSE)
                  {
                     pair->pbridgeto[i][j][pair->npbridge[i][j]] = prot_atom;
                     pair->npbridge[i][j]++;
                     total_bonds++;
                  }
               }
            }
            bond = bond->next_bond;
         }
       }
   }      
   
/* Update max_nbridge and max_npbridge ===============================*/
   if (phos_flag == NONPHOS)
   {
      if (total_bonds > *max_nbridge)
         *max_nbridge = total_bonds;
   }
   else if (phos_flag == PHOS)
      if (total_bonds > *max_npbridge)
         *max_npbridge = total_bonds;
   

/* Return to sort_interactions =======================================*/
   return(0);
}



/***********************************************************************
 *
 * bridgeint_resnum_order - Function to place the bridging interactions 
 *                          in resnum order
 *
 **********************************************************************/
int bridgeint_resnum_order(PAIRS *pair, int i, int nprot_int, int nwater_int, int phos_flag)
{
/* Declare variables =================================================*/
   int   j, k;
   int   nbonds;
   int   swap_flag;
   ATOMS *bridge_int[2], *swap_bondto;

/* Initialize variables ==============================================*/   
   if (phos_flag == NONPHOS)
      nbonds = pair->nbond[i];
   else if (phos_flag == PHOS)
      nbonds = pair->npbond[i];

/* Loop through bridging interactions ================================*/
   for (j = nprot_int; j < nbonds; j++)
   {
      swap_flag = TRUE;
      while (swap_flag == TRUE)
      {
         swap_flag = FALSE;
         if (phos_flag == NONPHOS)
         {
            for (k = 1; k < pair->nbridge[i][j]; k++)
            {
               bridge_int[0] = pair->bridgeto[i][j][k-1];
               bridge_int[1] = pair->bridgeto[i][j][k];         
      
               if (bridge_int[1]->resnum < bridge_int[0]->resnum)
               {
                  swap_bondto               = pair->bridgeto[i][j][k];
                  pair->bridgeto[i][j][k]   = pair->bridgeto[i][j][k-1];
                  pair->bridgeto[i][j][k-1] = swap_bondto;
                  swap_flag = TRUE;
               }
            }
         }
         else if (phos_flag == PHOS)
         {
            for (k = 1; k < pair->npbridge[i][j]; k++)
            {
               bridge_int[0] = pair->pbridgeto[i][j][k-1];
               bridge_int[1] = pair->pbridgeto[i][j][k];         
               
               if (bridge_int[1]->resnum < bridge_int[0]->resnum)
               {
                  swap_bondto                = pair->pbridgeto[i][j][k];
                  pair->pbridgeto[i][j][k]   = pair->pbridgeto[i][j][k-1];
                  pair->pbridgeto[i][j][k-1] = swap_bondto;
                  swap_flag = TRUE;
               } 
            }
         }
      }
   }
   

/* Return to sort_interactions =======================================*/
   return(0);
}



/***********************************************************************
 *
 * rm_dup_int - Function to remove duplicate interactions
 *
 **********************************************************************/
int
rm_dup_int(PAIRS *pair, int i, int phos_flag, int nonbond_by_residue)
{
/* Declare variables =================================================*/
  int   j, k;
  int bond_type;
  int duplicate, merge_flag;
  ATOMS *interaction[2];

/* Initialize variables ==============================================*/   

/* Loop through first level interactions =============================*/
  if (phos_flag == NONPHOS)
    {
      /* Initialise merge flag */
      merge_flag = TRUE;

      /* Loop while duplicates being encountered */
      while (merge_flag == TRUE)
        {
          /* Initialise merge flag */
           merge_flag = FALSE;

           /* Loop over all the interactions */
           for (j = 0; j < pair->nbond[i]; j++)
             {
               /* Save pointer to current interaction */
               interaction[0] = pair->bondto[i][j];
         
               /* Loop over all the other interactions */
               for (k = 0; k < j; k++)
                 {
                   /* Save pointer */
                   interaction[1] = pair->bondto[i][k];
                   bond_type = pair->type[i][k];

                   /* If interacting atom is the same as the one being
                      checked, then have a duplicate */
                   duplicate = FALSE;

                   /* Check whether the interactions belong to the
                      same residue */
                   if (interaction[0]->resnum == interaction[1]->resnum &&
                       interaction[0]->chain == interaction[1]->chain)
                     {
                       /* If atom is the same, then have a simple duplicate
                          and can combine the interactions */
                       if (!strncmp(interaction[0]->name,
                                interaction[1]->name, 4))
                         duplicate = TRUE;

                       /* Otherwise, have different atoms belonging to
                          the same residue */
                       else
                         {
                           /* If non-bonded interactions and these are being
                              combined, then have duplicate */
                           if (bond_type == NBOND &&
                               nonbond_by_residue == TRUE)
                             duplicate = TRUE;

                           /* If hydrogen bonds, then see if can be
                              combined */
                           else if (bond_type == HBOND)
                             {
                               /* Check that interactions haven't already
                                  been combined */
                               if (pair->combine[i][k] == NULL)
                                 {
                                   pair->combine[i][k] = interaction[0];
                                   duplicate = TRUE;
                                 }
                             }
                         }
                     }

                   /* If have a duplicate, then merge */
                   if (duplicate == TRUE)
                     {
                       /* Merge the bridging interactions, if any */
                       if (pair->nbridge[i][j] > 0 &&
                           pair->nbridge[i][k] > 0)
                         {
/* debug 
                           printf("*** HAVE TO MERGE BRIDGING INTERACTIONS!!\n");
 debug */
                         }

                       /* Move everything up one */
                       move_up_one(pair, i, j, pair->nbond[i], phos_flag);
               
                       /* Reduce total number of bonds and increment count
                          of dupliucates */
                       pair->nbond[i]--;
                       pair->ndups[i][k]++;
                       merge_flag = TRUE;
                     }
                 }
             }
         }
    }

  /* Process interactions to phosphate group */
  else if (phos_flag == PHOS)
    {
      /* Initialise merge flag */
      merge_flag = TRUE;

      /* Loop while duplicates being encountered */
      while (merge_flag == TRUE)
        {
          /* Initialise merge flag */
           merge_flag = FALSE;

           /* Loop over all the interactions to the phosphate backbone */
           for (j = 0; j < pair->npbond[i]; j++)
             {
               /* Save pointer to current interaction */
               interaction[0] = pair->pbondto[i][j];
         
               /* Loop over all the other interactions */
               for (k = 0; k < j; k++)
                 {
                   /* Save pointer */
                   interaction[1] = pair->pbondto[i][k];
                   bond_type = pair->ptype[i][k];

                   /* If interacting atom is the same as the one being checked,
                      then have a duplicate */
                   duplicate = FALSE;

                   /* Check whether the interactions belong to the
                      same residue */
                   if (interaction[0]->resnum == interaction[1]->resnum &&
                       interaction[0]->chain == interaction[1]->chain)
                     {
                       /* If atom is the same, then have a simple duplicate
                          and can combine the interactions */
                       if (!strncmp(interaction[0]->name,
                                interaction[1]->name, 4))
                         duplicate = TRUE;

                       /* Otherwise, have different atoms belonging to
                          the same residue */
                       else
                         {
                           /* If non-bonded interactions and these are being
                              combined, then have duplicate */
                           if (bond_type == NBOND &&
                               nonbond_by_residue == TRUE)
                             duplicate = TRUE;

                           /* If hydrogen bonds, then see if can be
                              combined */
                           else if (bond_type == HBOND)
                             {
                               /* Check that interactions haven't already
                                  been combined */
                               if (pair->pcombine[i][k] == NULL)
                                 {
                                   pair->pcombine[i][k] = interaction[0];
                                   duplicate = TRUE;
                                 }
                             }
                         }
                     }

                   /* If have a duplicate, then merge */
                   if (duplicate == TRUE)
                     {
                       /* Merge the bridging interactions, if any */
                       if (pair->npbridge[i][j] > 0 &&
                           pair->npbridge[i][k] > 0)
                         {
/* debug 
                           printf("*** HAVE TO MERGE BRIDGING INTERACTIONS!!\n");
 debug */
                         }

                       /* Move everything up one */
                       move_up_one(pair, i, j, pair->npbond[i], phos_flag);
                       
                       /* Reduce total number of bonds and increment count
                          of dupliucates */
                       pair->npbond[i]--;
                       pair->npdups[i][k]++;
                       merge_flag = TRUE;
                     }
                }
            }
        }
    }

/* Return to sort_interactions =======================================*/
  return(0);
}


/***********************************************************************
 *
 * move_up_one - Function to move a give set of interactions up a slot
 *
 **********************************************************************/
void
move_up_one(PAIRS *pair, int i, int j, int nbonds, int phos_flag)
{
  int l, m;

  /* Perform the move according to which set of interactions require
     moving */
  for (l = j; l < nbonds; l++)
    {
      if (phos_flag == NONPHOS)
        {
          pair->ndups[i][l] = pair->ndups[i][l+1];
          pair->bridge_flag[i][l] = pair->bridge_flag[i][l+1];
          pair->nbridge[i][l] = pair->nbridge[i][l+1];
          pair->type[i][l] = pair->type[i][l+1];
          pair->bondto[i][l] = pair->bondto[i][l+1];
          pair->combine[i][l] = pair->combine[i][l+1];
          for (m = 0; m < MAXARRAY; m++)
            {
              pair->bridgeto[i][l][m] = pair->bridgeto[i][l+1][m];
              pair->bridge_combine[i][l][m]
                = pair->bridge_combine[i][l+1][m];
            }
          pair->bondfrom[i][l] = pair->bondfrom[i][l+1];
        }
      else
        {
          pair->npdups[i][l] = pair->npdups[i][l+1];
          pair->pbridge_flag[i][l] = pair->pbridge_flag[i][l+1];
          pair->npbridge[i][l] = pair->npbridge[i][l+1];
          pair->ptype[i][l] = pair->ptype[i][l+1];
          pair->pbondto[i][l] = pair->pbondto[i][l+1];
          pair->pcombine[i][l] = pair->pcombine[i][l+1];
          for (m = 0; m < MAXARRAY; m++)
            {
              pair->pbridgeto[i][l][m] = pair->pbridgeto[i][l+1][m];
              pair->pbridge_combine[i][l][m]
                = pair->pbridge_combine[i][l+1][m];
            }
          pair->pbondfrom[i][l] = pair->pbondfrom[i][l+1];
        }
    }
}


/***********************************************************************
 *
 * define_object_sizes - Function to scale the picture to the page
 *
 **********************************************************************/
int define_object_sizes(PAIRS *beg_pair[MEDIARRAY], int ichain, SIZES *size,
                        PARAM *par, int max_nbond, int max_npbond, 
                        int max_nbridge, int max_npbridge, int ds_flag[2])
{

  /* Definitions of size-variables:-
                                        +-+
   +                  RESID(C) ATOM --- |P|
   |                                    +-+
   |                                       \
   |                                        \___
 A*|                              +-+        |  \___ BASE ___|___ BASE
   |                              |w| ------ |__/     G            C
   |                              +-+       /
   |                                       /
   |                            +-+     +-+
   +      RESID(C) ATOM --- NUM |w| --- |P|
                                +-+     +-+
          +------------+    +--+               +-------+
          resname_width  water_no_width    SUGAR_TO_BASE_WIDTH
         10 * text_width 3 * text_width
                       +----+                          +-----------+
                 resname_interact_width             BASE_TO_BASE_WIDTH
                    = 3 * text_width
                               +-+
                      water_no_to_water_width
                       = 1.5 * WATER_RADIUS
                                  +------+
                              WATER_TO_P_WIDTH
                                         +-----+
                                     P_TO_SUGAR_WIDTH

                                   +-----------+
                                WATER_TO_SUGAR_WIDTH

     where text_width = 0.71 * RESNAME_TXT

     A* is P_TO_P_HEIGHT

     Sizes:-

     WATER_RADIUS - Radius for water
     P_RADIUS     - Radius for phosphate group
     SUGAR_SIZE   - Parameter defining size of sugar
     RESNAME_TXT  - Size of protein residue name
     BASENAME_TXT - Size of single letter representing nucleic acid base

     Page definitions:-
                                      (BBOXX2,BBOXY2)
       +------------------------------------+  +
       |                                    |  |  TOP_MARGIN_Y
       |  +------------------------------+  |  +
       |  |                              |  |  |
       |  |                              |  |  |
       |  |                              |  |  |
       |  |                              |  |  |
       |  |                              |  |  |
       |  |                              |  |  |
       |  |                              |  |  |
       |  |                              |  |  |
       |  |                              |  |  |
       |  |                              |  |  |  y_height
       |  |              X               |  |  |
       |  |                              |  |  |
       |  |                              |  |  |
       |  |                              |  |  |
       |  |                              |  |  |
       |  |                              |  |  |
       |  |                              |  |  |
       |  |                              |  |  |
       |  |                              |  |  |
       |  |                              |  |  |
       |  |                              |  |  |
       |  * - - - - - - - - - - - - - - -|  |  +
       |  |                              |  |  |  key_height
       |  |                              |  |  |  (if key plotted)
       |  | - - - - - - - - - - - - - - -|  |  +
       |  |                              |  |  |  title_size (if used)
       |  +------------------------------+  |  +
       |                                    |  |  BOTTOM_MARGIN_Y
       +------------------------------------+  +
     (BBOXX1,BBOXY1)

       +--+                              +--+
     LEFT_MARGIN_X                     RIGHT_MARGIN_X

          +------------------------------+
                      x_width

     * is given by (ps_origin_x,ps_origin_y) and all coords are scaled to
       PostScript coords by ps_scale_factor
     X is given by (ps_centre_x,ps_centre_y)

     where title_size = TITLE_TXT * 1.5 (if plot title required)

     ----------------------------------------------------------------------*/



/* Declare variables =================================================*/
  int    npairs, pairs_per_page;
  int    single_strand;
  float  base_width, resname_width, text_width, title_size;
  float  total_height, total_width;
  float  resname_interact_width, water_no_width, water_no_to_water_width;
  float  x, y;
  float  key_width, key_height;
  float  ps_centre_x, ps_centre_y, ps_origin_x, ps_origin_y;
  float  picture_width, picture_height;
  float  ps_scale_factor;

  PAIRS  *pair;

/* Initialize variables ==============================================*/
  npairs = 0;
  pair   = beg_pair[ichain];
  single_strand = FALSE;
  if (par->landscape == TRUE)
    {
      key_width = KEY_WIDTH_LANDSCAPE * KEY_SCALE;
      key_height = KEY_HEIGHT_LANDSCAPE * KEY_SCALE;
    }
  else
    {
      key_width = KEY_WIDTH_PORTRAIT * KEY_SCALE;
      key_height = KEY_HEIGHT_PORTRAIT * KEY_SCALE;
    }

/* Initialize parameters =============================================*/

  /* Count the no. of pages the picture should be on ===================*/
  size->npages  = calc_no_of_pages(pair, par->numpairs, &npairs);

  /* Determine how many base pairs will be plotted on each page */
  if (size->npages == 1)
    pairs_per_page = npairs;
  else
    pairs_per_page = par->numpairs;

  /* If writing a separate plot for each groups of chains, then adjust
     page-total */
  if (par->write_to_nucplot_ps == FALSE)
    size->total_pages = size->npages;

/*======================================================================
 * Determine scale factor and centering position 
 *====================================================================*/

/*======================================================================
 * Text version 
 *====================================================================*/  
/* Precalculate widths that depend on text-sizes */
  text_width = (float) 0.71 * RESNAME_TXT;
  resname_width = 10 * text_width;
  water_no_width = 3 * text_width;
  if (par->landscape == TRUE)
    water_no_width = 1.5 * RESNAME_TXT;
  resname_interact_width = 3 * text_width;
  water_no_to_water_width = 1.5 * WATER_RADIUS;
  title_size = 2.0 * TITLE_TXT;

  /* Calculate width of base with no interactions */
  base_width = P_TO_SUGAR_WIDTH + SUGAR_TO_BASE_WIDTH
    + BASE_TO_BASE_WIDTH / 2.0;

  /* Add in width of residue name and its interaction distance */
  total_width = base_width + resname_width + resname_interact_width;

  /* If have bridging waters, then add additional distance required to
     allow for those */
  if (max_npbridge != 0 || max_nbridge != 0)
    total_width = total_width + water_no_width
      + water_no_to_water_width + WATER_RADIUS / (float) 2.0 + WATER_TO_P_WIDTH;

  /* If have two strands of DNA, then double the width */
  if (ds_flag[ONE] == TRUE && ds_flag[TWO] == TRUE)
    total_width = (float) 2.0 * total_width;
  else
    single_strand = TRUE;

/* Calculate total height required for each page =====================*/
  total_height = (float) (3.0 * CHAINAME_TXT + (pairs_per_page + 1) * P_TO_P_HEIGHT);

/* Determine extent of PostScript page that is available for the
   picture area ======================================================*/

  /* Get the picture width and height */
  picture_width = BBOXX2 - BBOXX1 - LEFT_MARGIN_X - RIGHT_MARGIN_X;
  picture_height = BBOXY2 - BBOXY1 - BOTTOM_MARGIN_Y - TOP_MARGIN_Y;

  /* Calculate reference origin in PostScript coordinates */
  ps_origin_x = BBOXX1 + LEFT_MARGIN_X;
  ps_origin_y = BBOXY1 + BOTTOM_MARGIN_Y;

  /* If plotting a title, adjust height to allow room for this */
  if (par->plot_title == TRUE)
    {
      /* Adjust height or width according to orientation */
      if (par->landscape == TRUE)
        {
          picture_width = picture_width - title_size;
          ps_origin_x = ps_origin_x + title_size;
        }
      else
        {
          picture_height = picture_height - title_size;
          ps_origin_y = ps_origin_y + title_size;
        }
    }

  /* If plotting a key, adjust height to allow room for this */
  if (par->plot_key == TRUE)
    {
      /* Adjust height or width according to orientation */
      if (par->landscape == TRUE)
        {
          picture_width = (float) (picture_width - key_height
            - KEY_TXT * KEY_SCALE);
          ps_origin_x = ps_origin_x + key_height;
        }
      else
        {
          picture_height = (float) (picture_height - key_height
            - 2.0 * KEY_TXT * KEY_SCALE);
          ps_origin_y = ps_origin_y + key_height;
        }
    }

  /* Get the PostScript coords of the centre-point on the page */
  ps_centre_x = (float) (ps_origin_x + picture_width / 2.0);
  ps_centre_y = (float) (ps_origin_y + picture_height / 2.0);

  /* Determine whether picture bounded in the x- or y-direction */
  if ((total_height / total_width) < (picture_height / picture_width))
    {
      ps_scale_factor = picture_width / total_width;
    }
  else
    {
      ps_scale_factor = picture_height / total_height;
    }

  /* Save the scaling factor */
  size->txtscale = ps_scale_factor;

/* Centering =========================================================*/
  size->x_centre = ps_centre_x;
  size->y_centre = ps_centre_y;

  /* Calculate the starting x-coordinate of the first base, relative
     to the centre of the picture */
  x = - BASE_TO_BASE_WIDTH / 2.0;

  /* If have only a single strand, then centre according to the width
     of the strand plus interactions */
  if (single_strand == TRUE)
    {
      if (ds_flag[ONE] == TRUE)
        x = (float) (total_width / 2.0 - BASE_TO_BASE_WIDTH / 2.0);
      else
        x = - (float) (total_width / 2.0 + BASE_TO_BASE_WIDTH / 2.0);
    }

  /* Calculate the starting y-coordinate of the first base, relative
     to the centre of the picture */
  y = (float) (((float) pairs_per_page / 2.0) * P_TO_P_HEIGHT - 2.0 * CHAINAME_TXT);

  /* Convert starting position of first base into PostScript coordinates
     and store */
  size->x_txtstart = size->x_centre + x * size->txtscale;
  size->y_txtstart = size->y_centre + y * size->txtscale;

  /* Calculate position for key and title, if required */
  size->y_title_pos = BBOXY1 + BOTTOM_MARGIN_Y + title_size / 2.0;
  size->y_keystart = BBOXY1 + BOTTOM_MARGIN_Y + key_height;
  if (par->plot_title == TRUE)
    size->y_keystart = size->y_keystart + 1.2 * TITLE_TXT;

  /* Store the filename as the plot title */
  if (par->plot_title == FALSE)
    par->title[0] = '\0';


/* Return to main ====================================================*/
  return(0);
}


/***********************************************************************
 *
 * calc_no_of_pages - Function to calculate the number of pages for the
 *                    given chain
 *
 **********************************************************************/
int
calc_no_of_pages(PAIRS  *pair, int numpairs, int *npair)
{
  int npairs, npages;

  /* Initialize variables */
  npairs = 0;
  npages = 1;

  /* Count the no. of pages the picture should be on */
  while (pair != NULL)
    {
      npairs++;
      if (npairs > numpairs)
        {
          npages++;
          npairs = 1;
        }
      pair = pair->next_pair;
    }

  /* Return the number of pairs and pages */
  *npair = npairs;
  return(npages);
}


/***********************************************************************
 *
 * open_psfile - Function to open PostScript files
 *
 **********************************************************************/
int
open_psfile(PAIRS *beg_pair[MEDIARRAY], int i, PARAM *par)
{
/* Declare variables =================================================*/
  char file_name[FILENAME_LEN];
  char chain[2];
  int   j;
  int   dot_pos;
  int   have_dot;
  ATOMS *atom;

/* Initialize variables ==============================================*/
  j             = 0;
  dot_pos       = 0;
  have_dot      = FALSE;
  atom          = beg_pair[i]->strand[ONE];

/* Name the Postscript file ==========================================*/

  /* If writing to nucplot.ps, set this as the filename */
  if (par->write_to_nucplot_ps == TRUE)
    strcpy(file_name,"nucplot.ps");
   
  /* Otherwise, form filename from name of input PDB file and the chain-ID
     of the first DNA chain */
  else
    {
      /* If have a chain-id, add that to the filename */
      strcpy(file_name,psfilename);
      if (atom != NULL && atom->chain != ' ')
        {
          strcat(file_name,"_");
          chain[0] = atom->chain;
          chain[1] = '\0';
          strcat(file_name,chain);
        }

      /* Add the .ps extension */
      strcat(file_name,".ps");
    }

/* Open the Postscript file ==========================================*/
  ps_file = open_file(file_name, "w");
  printf("\n");
  printf("Opened PostScript file %s ...\n",file_name);
  printf("\n");

/* Return to main ====================================================*/
  return(0);
}



/***********************************************************************
 *
 * psopen_ - Function to write out header records
 *
 **********************************************************************/
void psopen_(SIZES *size, PARAM *par)
{
/* Declare variables =================================================*/
  char   colno[3], exclamation, percent;
  double grey_shade, xx1, xx2, xy1, xy2;
  int    icolour, igrey;
/* v.1.1.7--> */
  int i;
/* <--v.1.1.7 */

/* Initialize variables ==============================================*/
  percent = '%';
  exclamation = '!';

/* Define border around picture ======================================*/
  if (par->landscape == FALSE)
    {
      xx1 = -1.0;
      xx2 = 650.0;
      xy1 = -1.0;
      xy2 = 951.0;
    }
  else if (par->landscape == TRUE)
    {
      xy1 = -1.0;
      xy2 = 650.0;
      xx1 = -1.0;
      xx2 = 951.0;
    }
  
/*======================================================================
 * Write out headings to PostScript file
 *====================================================================*/
  fprintf(ps_file,"%c%cPS-Adobe-3.0\n",percent,exclamation);
  fprintf(ps_file,"%%%%Creator: Nucplot\n");
  fprintf(ps_file,"%%%%DocumentNeededResources: font Times-Bold Symbol\n");
  fprintf(ps_file,"%%%%BoundingBox: (atend)\n");
  fprintf(ps_file,"%%%%Pages: %d\n", size->total_pages);
  if (par->title[0] != '\0')
    fprintf(ps_file,"%%%%Title: %s\n",par->title);
  fprintf(ps_file,"%%%%EndComments\n");
  fprintf(ps_file,"%%%%BeginProlog\n");
  fprintf(ps_file,"/L { moveto lineto stroke } bind def\n");
  fprintf(ps_file,"/G { gsave } bind def\n");
  fprintf(ps_file,"/W { setlinewidth } bind def\n");
  fprintf(ps_file,"/D { setdash } bind def\n");
  fprintf(ps_file,"/Col { setrgbcolor } bind def\n");
  fprintf(ps_file,"/Zero_linewidth { 0.0 } def\n");
  fprintf(ps_file,"/Sphcol { 1 setgray } def\n");
  fprintf(ps_file,"/Sphere {");
  fprintf(ps_file,"  newpath 3 copy 0 360 arc gsave Sphcol fill 0\n");
  fprintf(ps_file,"  setgray 0.5 setlinewidth\n");
  fprintf(ps_file,"  3 copy 0.94 mul 260 350 arc stroke 3 copy 0.87");
  fprintf(ps_file," mul 275 335 arc stroke\n");
  fprintf(ps_file,"  3 copy 0.79 mul 295 315 arc stroke 3 copy 0.8");
  fprintf(ps_file," mul 115 135 arc\n");
  fprintf(ps_file,"  3 copy 0.6 mul 135 115 arcn closepath gsave 1");
  fprintf(ps_file," setgray fill grestore stroke\n");
  fprintf(ps_file,"  3 copy 0.7 mul 115 135 arc stroke 3 copy 0.6");
  fprintf(ps_file," mul 124.9 125 arc\n");
  fprintf(ps_file,"  0.8 mul 125 125.1 arc stroke grestore stroke");
  fprintf(ps_file," } bind def\n");

/* Loop to write out the different grey-shades =======================*/
  for (igrey = 0; igrey < 11; igrey++)
  {
      grey_shade = (double)igrey / 10.0;
      sprintf(colno,"%2d",igrey);
      if (colno[0] == ' ')
          colno[0] = '0';
      fprintf(ps_file,"/Grey%s { %6.2f setgray } def\n",colno,grey_shade);
  }

/* Write out default black and white colours */
  fprintf(ps_file,
          "/Black           {   0.0000   0.0000   0.0000 setrgbcolor } def\n");
  fprintf(ps_file,
          "/White           {   1.0000   1.0000   1.0000 setrgbcolor } def\n");

/* Write out colour definitions ======================================*/
  if (par->colour == TRUE)
    for (icolour = 0; icolour < MAX_COLOURS; icolour++)
      {
        if (Colour_Table_Name[icolour][0] != '\0')
          {
            fprintf(ps_file,"/%s { %8.4f %8.4f %8.4f",
                    Colour_Table_Name[icolour],Colour_Table[icolour][0],
                    Colour_Table[icolour][1],Colour_Table[icolour][2]);
            fprintf(ps_file," setrgbcolor } def\n");
          }
      }

/* Write out background colour =======================================*/
  fprintf(ps_file, "/Background_col  { %s } def\n",par->background_colour);

/* Print the remainder of the definitions ============================*/
  fprintf(ps_file,"/Poly3 { moveto lineto lineto fill grestore } bind ");
  fprintf(ps_file,"def\n");
  fprintf(ps_file,"/Pl3 { 6 copy Poly3 moveto moveto moveto closepath ");
  fprintf(ps_file,"stroke } bind def\n");
  fprintf(ps_file,"/Pline3 { 6 copy Poly3 moveto lineto lineto closepa");
  fprintf(ps_file,"th stroke } bind def\n");
  fprintf(ps_file,"/LPhos3 { moveto -15.0 -20.00 rlineto 15.0 -20.00 rlineto");
  fprintf(ps_file," stroke } bind def\n");
  fprintf(ps_file,"/RPhos3 { moveto 15.0 -20.00 rlineto -15.0 -20.00 rlineto");
  fprintf(ps_file," stroke } bind def\n");
  fprintf(ps_file,"/Poly4 { moveto lineto lineto lineto fill grestore } ");
  fprintf(ps_file,"bind def\n");
  fprintf(ps_file,"/Pl4 { 8 copy Poly4 moveto moveto moveto moveto ");
  fprintf(ps_file,"closepath stroke } bind def\n");
  fprintf(ps_file,"/Pline4 { 8 copy Poly4 moveto lineto lineto lineto");
  fprintf(ps_file," closepath stroke } bind def\n");
  fprintf(ps_file,"/Pupoly4 { moveto 0.0 -23.0 rlineto 60.0 0.0 rlineto"); 
  fprintf(ps_file," 0.0 23.0 rlineto fill grestore } bind def\n");
  fprintf(ps_file,"/Pubox4 { 2 copy Pupoly4 moveto 0.0 -23.0 rlineto");
  fprintf(ps_file," 60.0 0.0 rlineto 0.0 23.0 rlineto closepath stroke } bind def\n");
  fprintf(ps_file,"/Pypoly4 { moveto 0.0 -23.0 rlineto 45.0 0.0 rlineto");
  fprintf(ps_file," 0.0 23.0 rlineto fill grestore } bind def\n");
  fprintf(ps_file,"/Pybox4 { 2 copy Pypoly4 moveto 0.0 -23.0 rlineto");
  fprintf(ps_file," 45.0 0.0 rlineto 0.0 23.0 rlineto closepath stroke } bind def\n");
  fprintf(ps_file,"/Circol { 1 setgray } def\n");
  fprintf(ps_file,"/Circle { gsave newpath 0 360 arc gsave Circol fill ");
  fprintf(ps_file,"grestore stroke grestore } bind def\n");
  fprintf(ps_file,"/Ucircle { gsave newpath 0 360 arc gsave Circol fill ");
  fprintf(ps_file,"grestore grestore } bind def\n");
  fprintf(ps_file,"/Arc { newpath arc stroke newpath } bind def\n");
/* v.1.1.7--> */
  fprintf(ps_file,"/Ocircle { gsave newpath 0 360 arc stroke grestore } ");
  fprintf(ps_file,"bind def\n");

  /* Loop over the ellipse colours used for highlighting variants */
  for (i = 0; i < NVAR_TYPES + 2; i++)
    fprintf(ps_file,"/Ellipse_col_%s { %s setrgbcolor } def\n",
            ELLIPSE_COL_DEF[i],ELLIPSE_COLOUR[i]);
  fprintf(ps_file,"/Ellipse { gsave translate scale Circle grestore } ");
  fprintf(ps_file,"bind def\n");
  fprintf(ps_file,"/Oellipse { gsave translate scale Ocircle grestore } ");
  fprintf(ps_file,"bind def\n");
/* <--v.1.1.7 */
  fprintf(ps_file,"/Print { /Times-Bold findfont exch scalefont setfont");
  fprintf(ps_file," show } bind def\n");
  fprintf(ps_file,"/Gprint { /Symbol findfont exch scalefont setfont show");
  fprintf(ps_file," } bind def\n");
  fprintf(ps_file,"/Lalign {\n");
  fprintf(ps_file,"  dup /Times-Bold findfont exch scalefont setfont\n");
  fprintf(ps_file,"  exch stringwidth pop pop -3 div 0.0 exch rmoveto\n");
  fprintf(ps_file," } bind def\n");
  fprintf(ps_file,"/Ralign {\n");
  fprintf(ps_file,"  dup /Times-Bold findfont exch scalefont setfont\n");
  fprintf(ps_file,"  exch stringwidth pop -1 mul exch -3 div rmoveto\n");
  fprintf(ps_file," } bind def\n");
  fprintf(ps_file,"/Center {\n");
  fprintf(ps_file,"  dup /Times-Bold findfont exch scalefont setfont\n");
  fprintf(ps_file,"  exch stringwidth pop -2 div exch -3 div rmoveto\n");
  fprintf(ps_file," } bind def\n");
  fprintf(ps_file,"/LalignRot90 {\n");
  fprintf(ps_file,"  dup /Times-Bold findfont exch scalefont setfont\n");
  fprintf(ps_file,"  exch stringwidth pop pop 3 div 0.0 exch exch rmoveto\n");
  fprintf(ps_file," } bind def\n");
  fprintf(ps_file,"/RalignRot90 {\n");
  fprintf(ps_file,"  dup /Times-Bold findfont exch scalefont setfont\n");
  fprintf(ps_file,"  exch stringwidth pop -1 mul exch 3 div exch rmoveto\n");
  fprintf(ps_file," } bind def\n");
  fprintf(ps_file,"/Center {\n");
  fprintf(ps_file,"  dup /Times-Bold findfont exch scalefont setfont\n");
  fprintf(ps_file,"  exch stringwidth pop -2 div exch -3 div rmoveto\n");
  fprintf(ps_file," } bind def\n");
  fprintf(ps_file,"/CenterRot90 {\n");
  fprintf(ps_file,"  dup /Times-Bold findfont exch scalefont setfont\n");
  fprintf(ps_file,"  exch stringwidth pop -2 div exch 3 div exch rmoveto\n");
  fprintf(ps_file," } bind def\n");
  fprintf(ps_file,"/UncenterRot90 {\n");
  fprintf(ps_file,"  dup /Times-Bold findfont exch scalefont setfont\n");
  fprintf(ps_file,"  exch stringwidth } bind def\n");
  fprintf(ps_file,"/Rot90 { gsave currentpoint translate 90 rotate }");
  fprintf(ps_file," bind def\n");
/* v.1.1.7--> */
  fprintf(ps_file,"/RotAngle { gsave currentpoint translate rotate } bind def\n");
/* <--v.1.1.7 */
  fprintf(ps_file,"%%%%EndProlog\n");
  fprintf(ps_file,"%%%%BeginSetup\n");
  fprintf(ps_file,"1 setlinecap 1 setlinejoin 1 setlinewidth 0 setgray");
  fprintf(ps_file," [ ] 0 setdash newpath\n");
  fprintf(ps_file,"%%%%EndSetup\n");

  /*fprintf(ps_file,"%%%%Page:  1 %2d\n", size->npages);
  fprintf(ps_file,"/NucplotSave save def\n");
  fprintf(ps_file,"%6.2f%6.2f moveto %7.2f%7.2f lineto%7.2f%7.2f lineto\n",
          xx1,xy1,xx2,xy1,xx2,xy2);
  fprintf(ps_file,"%7.2f%7.2f lineto closepath\n",xx1,xy2);
  fprintf(ps_file,"gsave 1.0000 setgray fill grestore\n");
  fprintf(ps_file,"stroke gsave\n");*/

/*Draw in the background colour ======================================*/
  /*pscomm_("Background");
  pshade_("Background_col");
  psubox_(BBOXX1, BBOXY1, BBOXX2, BBOXY1, BBOXX2, BBOXY2, BBOXX1, BBOXY2);*/
  
  /*if (Picture->In_Colour == TRUE)
    {
      pscomm_("Background");

      pshade_(Colour->Background);
      psubox_(BBOXX1, BBOXY1, BBOXX2, BBOXY1, BBOXX2, BBOXY2, BBOXX1, BBOXY2);
    }*/
}



/***********************************************************************
 *
 * base_pos - Function to assign the positions of the bases in the 
 *            postscript picture
 *
 **********************************************************************/
int base_pos(PAIRS *beg_pair, SIZES *size, PARAM *par, int ds_flag[2])
{
/* Declare variables =================================================*/
   int    npairs_done;
   double y_pos;
   PAIRS  *pair;
   
/* Initialize variables ==============================================*/
   npairs_done = 0;
   y_pos       = size->y_txtstart;
   pair        = beg_pair;

/* Loop through pairs ===============================================*/
   while (pair != NULL && npairs_done < par->numpairs)
   {

/* Calculate positions ==============================================*/
      /* y positions */
      pair->y_pos[ONE] = y_pos;
      pair->y_pos[TWO] = y_pos;

      /* Increment the y-position */
      y_pos = y_pos - P_TO_P_HEIGHT * size->txtscale;

      /* x positions */
      pair->x_pos[ONE] = size->x_txtstart;
      pair->x_pos[TWO] = size->x_txtstart
        + BASE_TO_BASE_WIDTH * size->txtscale;

      /* If only the second strand is present, then set its x-position
         accordingly */
      if (ds_flag[ONE] == FALSE)
        pair->x_pos[TWO] = size->x_txtstart;

      /* numbering position */
      if (par->numbering == FALSE)
      {
         pair->x_numpos[ONE] = pair->x_pos[ONE];
         pair->x_numpos[TWO] = pair->x_pos[TWO];
      }
      else if (par->numbering == TRUE)
      {
         pair->x_numpos[ONE] = pair->x_pos[ONE];
         pair->x_numpos[TWO] = pair->x_pos[TWO];
      }
         
      pair = pair->next_pair;
      npairs_done++;
   }
   
/* Return to draw_pictures ===========================================*/
   return(0);
}



/***********************************************************************
 *
 * text_nucleic_acids - Function to draw the nucleic acid pictures out to
 *                      Postscript file (text version)
 *
 **********************************************************************/
PAIRS *text_nucleic_acids(PAIRS *beg_pair[MEDIARRAY], PAIRS *end_pair[MEDIARRAY], 
                          PAIRS *first_pair, int pairno, SIZES *size, PARAM *par, int ds_flag[2])
{
/* Declare variables =================================================*/
   char   chain_top[2], chain_bottom[2], code[4];
   int    npairs;
   double text_size, x, y, x1, y1, x2, y2;
   double cx, cy;
   int    save_x1, save_x2, x_centre;
   PAIRS  *pair, *last_pair;
   ATOMS  *base1, *base2;

/* Initialize variables ==============================================*/
   cx = size->x_centre;
   cy = size->y_centre;
   npairs    = 0;
   last_pair = NULL;
   pair      = first_pair;
   chain_top[ONE] = '\0';
   chain_top[TWO] = '\0';
   chain_bottom[ONE] = '\0';
   chain_bottom[TWO] = '\0';

/* Comment ===========================================================*/
   pscomm_("Nucleic Acids");

/*======================================================================
 * Write out bases
 *====================================================================*/
/* Set up font sizes, positions etc ==================================*/
   
/* Loop through base pairs ===========================================*/      
   while ((pair != NULL) && (npairs < par->numpairs))
   {
      /* Draw strand one */
      if (pair->strand[ONE] != NULL)
      {
         pssave_();

         if (pair->nbond[ONE] != 0 && par->highlight == TRUE)
         {
            fprintf(ps_file, "0.0 setgray\n");
         }
            
         /* Resname */
         base1 = pair->strand[ONE];
         pstrans_(pair->x_pos[ONE], pair->y_pos[ONE], cx, cy, &x1, &y1, par);
         psbasecol_(base1, par, code);
         text_size = BASENAME_TXT*size->txtscale;
         if ((int) strlen(code) > 1)
           text_size = text_size / 3.0;
         psctxt_(x1, y1, (int) text_size, code);
         save_x1 = (int) pair->x_pos[ONE];
         if (chain_top[ONE] == '\0')
           chain_top[ONE] = base1->chain;
         chain_bottom[ONE] = base1->chain;

         psrest_();
      }
      
      /* Draw strand two */
      if (pair->strand[TWO] != NULL)
      {
         pssave_();
         if (pair->nbond[TWO] != 0 && par->highlight == TRUE)
         {
            fprintf(ps_file, "0.0 setgray\n");
         }

         /* Resname */
         base2 = pair->strand[TWO];
         pstrans_(pair->x_pos[TWO], pair->y_pos[TWO], cx, cy, &x1, &y1, par);
         psbasecol_(base2, par, code);
         text_size = BASENAME_TXT*size->txtscale;
         if ((int) strlen(code) > 1)
           text_size = text_size / 3.0;
         psctxt_(x1, y1, (int) text_size, code);
         save_x2 = (int) pair->x_pos[TWO];
         if (chain_top[TWO] == '\0')
           chain_top[TWO] = base2->chain;
         chain_bottom[TWO] = base2->chain;
 
         /* Res number */
         /*if (par->numbering == TRUE)
           psbasenum_(pair, size, par, RIGHT);*/
         
         psrest_();
      }
      
      /* Draw the bond line between the bases */
      if ((pair->strand[ONE] != NULL)
          && (pair->strand[TWO] != NULL))
      {
         pssave_();
         pslwid_(BSLNWID * size->txtscale);

         /* Calculate start- and end-coords of line between the two
            bases */
         x_centre = (int) ((save_x1 + save_x2) / 2.0);
         x = size->txtscale * (BASE_TO_BASE_WIDTH - BASENAME_TXT) / 2.0;
         y = pair->y_pos[ONE];

         /* Calculate the transformed coordinates (for landscape version) */
         pstrans_(x_centre - x, y, cx, cy, &x1, &y1, par);
         pstrans_(x_centre + x, y, cx, cy, &x2, &y2, par);
         if (par->highlight == TRUE)
           psetgray_(0.0);

         /* Draw the bond line, unless bases are mismatched */
         if (pair->paired == TRUE)
           {
             pscolb_(par->bbline_colour);
             psline_(x1, y1, x2, y2);
           }
         psrest_();
      }
      
      npairs++;
      last_pair = pair;
      pair = pair->next_pair;
   }
   
/* Write out 5' end or dashed lines for top end ======================*/
   psends_(beg_pair, first_pair, size, par, pairno, TOP, ds_flag,
           chain_top, chain_bottom);
  
/* Write out 5' end or dashed lines for bottom ends ==================*/
   psends_(end_pair, last_pair, size, par, pairno, BOTTOM, ds_flag,
           chain_top, chain_bottom);

/* Return to main ====================================================*/
   return(pair);
}

/* v.1.1.7--> */
/***********************************************************************

string_truncate  -  Truncate the given string at the last non-blank
                    character, or at its maximum length

***********************************************************************/

int string_truncate(char *string,int max_length)
{
  char ch;

  int iend, ipos, len;

  /* Initialise variables */
  iend = -1;
  ch = string[0];
  for (ipos = 0; ipos < max_length && ch != '\0'; ipos++)
    {
      /* Get the current character and check for end of string */
      ch = string[ipos];
      if (ch == '\n' || ch == 13)
        ch = '\0';
      if (ch != '\0')
        {
          /* If not a space, then mark end of string */
          if (ch != ' ')
            iend = ipos;
        }
    }

  /* End the string after the last non-blank character */
  if (iend + 1 < max_length)
    {
      string[iend + 1] = '\0';
      len = iend + 1;
    }
  else
    {
      string[max_length - 1] = '\0';
      len = max_length - 1;
    }

  /* Return the string length */
  return(len);
}
/***********************************************************************

to_upper  -  Convert given character string to upper case

***********************************************************************/

void to_upper(char *search_string,int length)
{
  int ichar;

  int ipos;

  /* Loop through all the character positions, converting any letters
     to upper case */
  for (ipos = 0; ipos < length; ipos++)
    {
      /* Check the current character */
      ichar = search_string[ipos] - 'a';
      if (ichar >= 0 && ichar < 26)
        {
          ichar = ichar + 'A';
          search_string[ipos] = ichar;
        }
    }
}
/***********************************************************************

get_other_chains  -  Read in the chain equivalences from the chains.dat
                   file and return the other chain in the interface if
                   it is equivalent to the one with the highlighted
                   residue

***********************************************************************/

int get_other_chains(char chain,char *chains,PARAM *par)
{
  char file_name[FILENAME_LEN + 1], input_line[LINE_LEN + 1];
  char chain_id, copy_of, other_chain, store_chain;

  int ichain, len, nchains, nlines;

  FILE *file_ptr;

  /* Initialise variables */
  nchains = 0;
  nlines = 0;
  other_chain = '\0';

  /* Store our chain */
  chains[nchains] = chain;
  nchains++;

  /* Form the name of the chains.dat file in the PDBsum directory */
  strcpy(file_name,par->pdbsum_dir);
  strcat(file_name,"chains.dat");

  /* Open the chains.dat input file */
  if ((file_ptr = fopen(file_name,"r")) !=  NULL)
    {
      /* Loop while reading in records from the input file */
      while (fgets(input_line,LINE_LEN,file_ptr)!= NULL)
        {
          /* Increment count of lines read in */
          nlines++;

          /* Truncate the string to strip off the newline */
          len = string_truncate(input_line,LINE_LEN);

          /* If this is a copy of our chain, store it */
          if (len > 2)
            {
              /* Get the two chain ids */
              chain_id = input_line[0];
              copy_of = input_line[2];
              store_chain = '\0';

              /* If either chain corresponds to our chain, store
                 the other one as one of the equivalents */
              if (chain_id == chain)
                store_chain = copy_of;
              else if (copy_of == chain)
                store_chain = chain_id;

              /* If have a chain to store, do so */
              if (store_chain != '\0' && nchains < MAXCHAINS)
                {
                  /* Store this chain */
                  chains[nchains] = store_chain;

                  /* Increment count of stored chains */
                  nchains++;
                }
            }
        }

      /* Close the input file */
      fclose(file_ptr);
    }

  /* Show warning tnat chain.dat file not found */
  else
    printf("*** Warning. File not found: %s\n",file_name);

  /* Return number of chains */
  return(nchains);
}
/***********************************************************************

find_highres_entry  -  Find the given record

***********************************************************************/

struct highres *find_highres_entry(struct highres *first_highres_ptr,
                                   char *res_name,char *res_num,
                                   char chain)
{
  struct highres *highres_ptr, *found_highres_ptr;

  /* Initialise variables */
  found_highres_ptr = NULL;

  /* Get pointer to the first record */
  highres_ptr = first_highres_ptr;

  /* Loop over all the records to find the one required */
  while (highres_ptr != NULL && found_highres_ptr == NULL)
    {
      /* If have a match, store the pointer */
      if (chain == highres_ptr->chain &&
          !strcmp(res_num,highres_ptr->res_num) &&
          !strcmp(res_name,highres_ptr->res_name))
        found_highres_ptr = highres_ptr;

      /* Get pointer to the next domain */
      highres_ptr = highres_ptr->next_highres_ptr;
    }

  /* Return the matching record pointer */
  return(found_highres_ptr);
}
/***********************************************************************

create_highres_entry  -  Create a record holding residue to be
                         highlighted

***********************************************************************/

struct highres *create_highres_entry(struct highres **fst_highres_ptr,
                                     struct highres **lst_highres_ptr,
                                     int type,int ref,char *res_name,
                                     char *res_num,char chain)
{
  struct highres *highres_ptr;
  struct highres *first_highres_ptr, *last_highres_ptr;

  /* Initialise pfam pointers */
  highres_ptr = NULL;
  first_highres_ptr = *fst_highres_ptr;
  last_highres_ptr = *lst_highres_ptr;

  /* Create pfam entry */
  highres_ptr = (struct highres *) malloc(sizeof(struct highres));
  if (highres_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for struct highres\n");
      printf("***        Program terminated with error.\n");
      exit (1);
    }

  /* Store the known data */
  highres_ptr->type = type;
  highres_ptr->ref = ref;
  highres_ptr->chain = chain;
  strcpy(highres_ptr->res_name,res_name);
  strcpy(highres_ptr->res_num,res_num);

  /* Initialise pointer to next entry */
  highres_ptr->next_highres_ptr = NULL;

  /* If this is the first entry stored, then update
     pointer to head of linked list */
  if (first_highres_ptr == NULL)
    first_highres_ptr = highres_ptr;

  /* Otherwise, make previous entry point to the current
     one */
  else
    last_highres_ptr->next_highres_ptr = highres_ptr;
                      
  /* Save the current residue's pointer */
  last_highres_ptr = highres_ptr;

  /* Return current pointers */
  *fst_highres_ptr = first_highres_ptr;
  *lst_highres_ptr = last_highres_ptr;
  return(highres_ptr);
}
/***********************************************************************

format_res_number  -  Interpret the residue-number from the input
                      parameter for the residue to be highlighted

***********************************************************************/

void format_res_number(char *tmpstr,char *res_num)
{
  char ch, insertion_code, new_number[6], res_number[6];
  int cpos, idigit, ipos, ndigits;
  int nchars;
  int at_end, got_number;

  /* Initialise variables */
  at_end = FALSE;
  got_number = FALSE;
  nchars = 0;
  ndigits = 0;
  strcpy(new_number,"     ");
  strcpy(res_number,"     ");

  /* Initialise position of insertion code */
  cpos = 0;
  insertion_code = ' ';

  /* Loop through the characters searching for end of number or
     insertion-code */
  for (ipos = 0; ipos < 5 && at_end == FALSE; ipos++)
    {
      /* Get this character */
      ch = tmpstr[ipos];

      /* If this is the end of the string, store end-position */
      if (ch == '\0')
        {
          cpos = ipos;
          at_end = TRUE;
        }

      /* If this is a digit, then have part of a residue number */
      else if ((ch >= '0' && ch <= '9') ||
               ch == '-')
        {
          res_number[ndigits] = ch;
          ndigits++;
          got_number = TRUE;
        }

      /* If this is a character it might be the insertion code */
      else if ((ch >= 'A' && ch <= 'Z') ||
               ch == '\'')
        {
          /* Store position of insertion code */
          cpos = ipos;
          insertion_code = ch;

          /* Increment count of characters */
          nchars++;
        }

      /* If this is a space, have reached the end of the number */
      else if (ch == ' ' || ch == '\0')
        {
          if (got_number == TRUE)
            cpos = ipos;
        }

      /* Otherwise, have an invalid character in the residue
         number */
      else
        cpos = -1;
    }

  /* If have too many digits, then last one is an insertion code */
  if (ndigits > 4)
    {
      insertion_code = res_number[4];
      ndigits--;
    }

  /* Reformat the number */
  if (cpos == 0)
    cpos = 4;

  /* Now re-form the number such the first four characters
     contain the digits and the 5th contains the insertion
     code (if any) */
  for (idigit = 0; idigit < ndigits; idigit++)
    {
      ipos = idigit + (4 - ndigits);
      new_number[ipos] = res_number[idigit];
    }

  /* Add the insertion code */
  new_number[4] = insertion_code;
  new_number[5] = '\0';

  /* Transfer the new number across */
  strncpy(res_num,new_number,5);
  res_num[5] = '\0';
}
/***********************************************************************

extract_residue_name  -  Extract resisdue name, number and chain id

***********************************************************************/

void extract_residue_name(char *highlight_res,char *res_name,
                          char *res_num,char *chn)
{
  char chain, number_string[20], tmp_string[20];

  char *string_ptr;

  int name_end;

  /* Initialise variables */
  res_name[0] = res_num[0] = '\0';
  chain = ' ';
  
  /* Check for square brackets defining residue name */
  if (highlight_res[0] == '[')
    {
      /* Get the name bounded by square brackets */
      strcpy(tmp_string,highlight_res + 1);

      /* Check for close-bracket */
      string_ptr = strstr(tmp_string,"]");
      if (string_ptr != NULL)
        {
          string_ptr[0] = '\0';
          name_end = (string_ptr - tmp_string) + 1;
        }
      else
        {
          tmp_string[3] = '\0';
          name_end = 3;
        }

      /* Store the residue name */
      strcpy(res_name,tmp_string);
    }

  /* Otherwise, take residue name to be first 3 characters */
  else
    {
      /* Transfer name */
      strncpy(res_name,highlight_res,3);
      res_name[3] = '\0';
      name_end = 2;
    }

  /* Convert residue name to upper case */
  to_upper(res_name,3);

  /* Get the residue number */
  strcpy(number_string,highlight_res + name_end + 1);
  string_ptr = strstr(number_string,":");
  if (string_ptr != NULL)
    {
      /* Extract the chain identifier after the colon */
      chain = string_ptr[1];
      string_ptr[0] = '\0';
    }

  /* Convert the residue number into the standard fixed format */
  format_res_number(number_string,res_num);

  /* Return the chain id */
  *chn = chain;
}
/***********************************************************************

get_variants  -  Read in the DisaStr veriants on the current LIGPLOT
                 diagram

***********************************************************************/

int get_variants(struct highres **fst_highres_ptr,
                 struct highres **lst_highres_ptr,int nhigh,
                 char *chains,int nchains,PARAM *par)
{
  char file_name[FILENAME_LEN], input_line[LINE_LEN + 1];
  char chain, res_name[4], res_num[6];
  char number_string[20], tmp[LINE_LEN + 1];

  char *string_ptr;

  int ichain, len, line, type;
  int ref, wanted;

  struct highres *highres_ptr, *first_highres_ptr, *last_highres_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  line = 0;
  ref = FALSE;
  type = NO_VARIANT;
  wanted = FALSE;

  /* Initialise pointers */
  first_highres_ptr = *fst_highres_ptr;
  last_highres_ptr = *lst_highres_ptr;
  highres_ptr = NULL;

  /* Form name of input annotations.dat file */
  strcpy(file_name,par->disastr_dir);
  strcat(file_name,"/wiring/annotations.dat");

  /* If file not found, abort */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    {
      printf("*** ERROR. File not found: %s\n",file_name);
      return(nhigh);
    }

  /* Loop while reading in records from input file */
  while (fgets(input_line,LINE_LEN,file_ptr) != NULL)
    {
      /* Truncate the input line */
      len = string_truncate(input_line,LINE_LEN);

      /* Get the position of the space */
      string_ptr = strchr(input_line,' ');

      /* If this is the start of a new residue, then re-initialise */
      if (!strncmp(input_line,"KEY: ",5))
        {
          highres_ptr = NULL;
          type = NO_VARIANT;
          wanted = FALSE;
        }

      /* If interaction type, save the type */
      else if (!strncmp(input_line,"INT_VAR_TYPE ",13))
        {
          /* Get the variant type */
          if (string_ptr != NULL)
            {
              /* Get the type */
              strcpy(number_string,string_ptr + 1);
              type = atoi(number_string);
              if (type == -9)
                type = NO_VARIANT;
            }
        }

      /* If this is the start of a new residue, then re-initialise */
      else if (!strncmp(input_line,"INT_PDB_CODE[",13))
        {
          /* Reinitialise flag */
          wanted = FALSE;

          /* Check if this is the right PDB code */
          if (string_ptr != NULL && !strcmp(string_ptr + 1,par->pdb_code))
            wanted = TRUE;
        }

      /* Only consider line if we have the right type of interaction */
      if (wanted == TRUE && string_ptr != NULL)
        {
          /* If this is the residue name, then save */
          if (!strncmp(input_line,"INT_RES_NAME[",13))
            strcpy(res_name,string_ptr + 1);

          /* If this is the residue number, then save */
          else if (!strncmp(input_line,"INT_RES_NUM[",12))
            {
              /* Get the raw number */
              strcpy(tmp,string_ptr + 1);
              
              /* Format the number */
              format_res_number(tmp,res_num);
            }
              
          /* If this is the chain, then save */
          else if (!strncmp(input_line,"INT_CHAIN[",10))
            chain = string_ptr[1];

          /* If this is the LIGPLOT identifer, save the residue */
          else if (!strncmp(input_line,"INT_LIGPLOT[",12))
            {
              /* See if the current chain is one of the wanted ones */
              wanted = FALSE;
              for (ichain = 0; ichain < nchains; ichain++)
                {
                  /* If the right chain, set the flag */
                  if (chains[ichain] == chain)
                    wanted = TRUE;
                }
              
              /* Loop over all the equivalent chains, to create a
                 highres record for each */
              for (ichain = 0; ichain < nchains && wanted == TRUE;
                   ichain++)
                {
                  /* See if we already have this residue */
                  highres_ptr
                    = find_highres_entry(first_highres_ptr,res_name,
                                         res_num,chains[ichain]);

                  /* If not the reference residue, then create a record
                     for it */
                  if (highres_ptr == NULL)
                    {
                      /* Create record */
                      highres_ptr
                        = create_highres_entry(&first_highres_ptr,
                                               &last_highres_ptr,type,ref,
                                               res_name,res_num,
                                               chains[ichain]);

                      /* Increment count of variants stored */
                      nhigh++;
                    }

                  /* If this is the reference residue, then save the type */
                  else if (highres_ptr->ref == TRUE &&
                           highres_ptr->type == NO_VARIANT)
                    highres_ptr->type = type;
                }
            }
        }
    }

  /* Close the input file */
  fclose(file_ptr);

  /* Return the pointers */
  *fst_highres_ptr = first_highres_ptr;
  *lst_highres_ptr = last_highres_ptr;

  /* Return the number of highlighted residues */
  return(nhigh);
}
/***********************************************************************

store_variants  -  If called from DisaStr, pick up all the disease-
                   associated variants for the given UniProt entry

***********************************************************************/

int store_variants(struct highres **fst_highres_ptr,PARAM *par)
{
  char chain, chains[MAXCHAINS], res_name[4], res_num[6];

  int ichain, nchains, nhigh, var_type;
  int ref;

  struct highres *highres_ptr, *first_highres_ptr, *last_highres_ptr;
    
  /* Initialise variables */
  nhigh = 0;
  ref = TRUE;
  var_type = NO_VARIANT;

  /* Initialise pointers */
  highres_ptr = first_highres_ptr = last_highres_ptr = NULL;

  /* If we have a residue to highlight, then process variants */
  if (par->highlight_res[0] != '\0')
    {
      /* Unpack the name of the highlighted residue */
      extract_residue_name(par->highlight_res,res_name,res_num,
                           &chain);

      /* Pick up the chains from the chains.dat file, if present */
      nchains = get_other_chains(chain,chains,par);

      /* Loop over all the equivalent chains, to create a highres record
         for each */
      for (ichain = 0; ichain < nchains; ichain++)
        {
          /* Create a highlight record for this residue */
          highres_ptr
            = create_highres_entry(&first_highres_ptr,&last_highres_ptr,
                                   var_type,ref,res_name,res_num,chain);
          nhigh++;

          /* Save the variant type */
          highres_ptr->type = var_type;
        }

      /* Get all the other residues to be highlighted */
      nhigh = get_variants(&first_highres_ptr,&last_highres_ptr,nhigh,
                           chains,nchains,par);
    }

  /* Return current pointers */
  *fst_highres_ptr = first_highres_ptr;

  /* Return number of variant residues stored for highlighting */
  return(nhigh);
}  
/* <--v.1.1.7 */


/***********************************************************************
 *
 * text_backbone - Function to draw the DNA backbone
 *
 **********************************************************************/
int text_backbone(PAIRS *first_pair, PAIRS *last_pair, SIZES *size, PARAM *par)
{
/* Declare variables =================================================*/
   int    i, j;
   int    proceed_sugar, proceed_phos, top_of_page;
   int    drawn_flag, drawn_sugar, first_base, last_base;
   int    is_sugar, is_phosphate;
   double xstart, xstart_left, ystart;
   double x_line[2][2], y_line[2][2];
   double prev_xline[2][2], prev_yline[2][2];
   ATOMS  *atom, *first_atom, *next_base_atom, *first_n_b_atom;
   PAIRS  *pair, *next_pair;

/* Initialize variables ==============================================*/
/* v.1.1.5--> */
   proceed_sugar = FALSE;
/* <--v.1.1.5 */
   pair = first_pair;
   pscomm_("Backbone");
   for (i = 0; i < 2; i++)
     for (j = 0; j < 2; j++)
       {
         x_line[i][j] = 0.0;
         y_line[i][j] = 0.0;
       }
   top_of_page = TRUE;

/* Loop through pairs ================================================*/   
   while (pair != last_pair)
   {

/* Do left and right strands =========================================*/
      for (i = 0; i < 2; i++)
      {
         atom = first_atom = pair->strand[i];
         next_pair = pair->next_pair;
         drawn_flag = drawn_sugar = FALSE;

/* Loop through DNA residue ==========================================*/
         while (atom != NULL && atom->resnum == first_atom->resnum)
         {

/* Sugar =============================================================*/
            is_sugar = is_sugar_atom(atom->name);
            if (drawn_sugar == FALSE && is_sugar == TRUE)
              proceed_sugar = TRUE;

            if (proceed_sugar == TRUE)
            {
               /* Left */
               if (i == ONE)
               {
                  /* Save previous sugar's coordinates */
                  for (j = 0; j < 2; j++)
                  {
                     prev_xline[LEFT][j] = x_line[LEFT][j];
                     prev_yline[LEFT][j] = y_line[LEFT][j];
                  }
                  xstart = pair->x_numpos[ONE]
                    - SUGAR_TO_BASE_WIDTH * size->txtscale;
                  ystart = pair->y_pos[ONE];
                  psugar_(pair, size, par, i, xstart, ystart, x_line[LEFT],
                          y_line[LEFT]);
                  psbasenum_(pair, size, par, i, xstart, ystart);

               }

               /* Right */
               else if (i == TWO)
               {
                  /* Save previous sugar's coordinates */
                  for (j = 0; j < 2; j++)
                  {
                     prev_xline[RIGHT][j] = x_line[RIGHT][j];
                     prev_yline[RIGHT][j] = y_line[RIGHT][j];
                  }
                  xstart = pair->x_numpos[TWO]
                    + SUGAR_TO_BASE_WIDTH * size->txtscale;
                  ystart  = pair->y_pos[TWO];
                  psugar_(pair, size, par, i, xstart, ystart, x_line[RIGHT],
                          y_line[RIGHT]);
                  psbasenum_(pair, size, par, i, xstart, ystart);
               }
               drawn_sugar = TRUE;
               proceed_sugar = FALSE;
            }
 
/* Phosphate =========================================================*/
            is_phosphate = is_phosphate_atom(atom->name);
            if (drawn_flag == FALSE && is_phosphate == TRUE)
              proceed_phos = TRUE;
            
            if (proceed_phos == TRUE && drawn_sugar == TRUE) 
            {
              /* Determine whether this is the first or last base of a
                 chain */
              first_base = pair->first_base[i];
              last_base = pair->last_base[i];

               /* Left */
               if (i == ONE)
               {
                  xstart = pair->x_numpos[ONE]
                    - (SUGAR_TO_BASE_WIDTH + P_TO_SUGAR_WIDTH)
                      * size->txtscale;
                  xstart_left = xstart;
                  if (pair != first_pair)
                     psphostxt_(pair, size, par, xstart, x_line[LEFT],
                                y_line[LEFT], prev_xline[LEFT],
                                prev_yline[LEFT], LEFT, FALSE,
                                first_base, last_base, top_of_page);
                  else if (pair == first_pair)
                     psphostxt_(pair, size, par, xstart, x_line[LEFT],
                                y_line[LEFT], prev_xline[RIGHT],
                                prev_yline[RIGHT], LEFT, FALSE,
                                first_base, last_base, top_of_page);
               }

               /* Right */
               else if (i == TWO)
               {
                  xstart = pair->x_numpos[TWO]
                    + (SUGAR_TO_BASE_WIDTH + P_TO_SUGAR_WIDTH)
                      * size->txtscale;
                  if (pair != last_pair)
                     psphostxt_(pair, size, par, xstart, x_line[RIGHT],
                                y_line[RIGHT], prev_xline[RIGHT],
                                prev_yline[RIGHT], RIGHT, FALSE,
                                first_base, last_base, top_of_page);
               }
               drawn_flag = TRUE;
               proceed_phos = FALSE;
            }
            atom = atom->next_atom;
         }

         /* Draw next phosphate if last base on the page and 
            the next base has a phosphate */
         if (last_pair != NULL && i == ONE)
            if (pair == last_pair->prev_pair)
            {
              /* Determine whether this is the first base of a chain */
              first_base = FALSE;
              last_base = pair->last_base[i];

               next_base_atom = first_n_b_atom = next_pair->strand[i];
               drawn_flag = FALSE;
               while(drawn_flag == FALSE && next_base_atom != NULL 
                     && next_base_atom->resnum == first_n_b_atom->resnum)
               {
                 is_phosphate = is_phosphate_atom(next_base_atom->name);
                 if (is_phosphate == TRUE)
                   {
                      xstart = xstart_left;
                      psphostxt_(pair, size, par, xstart, x_line[LEFT],
                                 y_line[LEFT], prev_xline[LEFT],
                                 prev_yline[LEFT], LEFT, TRUE,
                                 first_base, last_base, top_of_page);
                      drawn_flag = TRUE;
                    }
                  next_base_atom = next_base_atom->next_atom;
               }
            }
       }         
      top_of_page = FALSE;
      pair = pair->next_pair;
   }
   
/* Return to draw_picture ============================================*/
   return(0);
}



/***********************************************************************

is_sugar_atom  -  Check whether this is an atom belonging to the ring
                  of the sugar backbone

***********************************************************************/

int
is_sugar_atom(char atname[5])
{
  int is_sugar;

  /* Initialise */
  is_sugar = FALSE;

  /* Check the atom name for any indicators suggesting this might belong
     to the sugar ring */
  if (atname[3] == '*' || atname[3] == '\'')
    {
      if (!strncmp(atname, " O3", 3) || !strncmp(atname, " O5", 3))
        is_sugar = FALSE;
      else
        is_sugar = TRUE;
    }

  /* Return the answer */
  return is_sugar;
}



/***********************************************************************

is_phosphate_atom  -  Check whether this is an atom belonging to the
                      phosphate backbone

***********************************************************************/

int
is_phosphate_atom(char atname[5])
{
  int is_phosphate;

  /* Initialise */
  is_phosphate = FALSE;

  /* Check the atom name for any indicators suggesting this might belong
     to the phosphate backbone */
  if ((atname[3] == '*' || atname[3] == '\'') &&
      (!strncmp(atname, " O3", 3) || !strncmp(atname, " O5", 3)))
    is_phosphate = TRUE;
/* v.1.1.3--> */
//  if (atname[1] == 'P' || atname[3] == 'P')
  if (atname[1] == 'P' || atname[2] == 'P' || atname[3] == 'P')
/* <--v.1.1.3 */
    is_phosphate = TRUE;

  /* Return the answer */
  return is_phosphate;
}



/***********************************************************************
 *
 * find_bonded_base - Function to find the right DNA base for each bond
 *
 **********************************************************************/
int
find_bonded_base(PAIRS *first_pair, PAIRS *last_pair, BONDS *first_bond, 
                 PARAM *par, int *max_nbond, int *max_npbond)
{
/* Declare variables =================================================*/
  int   i, j, k, l;
  ATOMS *strand[2];
  ATOMS *bond_atom[2];
  PAIRS *pair, *next_pair, *prev_pair;
  BONDS *bond;

/* Initialize variables ==============================================*/
  bond = first_bond;
  l = 0;
   
/* Loop through bonds ================================================*/
  while (bond != NULL)
    {
      /* Check if a water or not */
      bond_atom[0] = bond->bond_atom[0];
      bond_atom[1] = bond->bond_atom[1];

      /* Loop over the base-pairs */
      pair = first_pair;
/* v.1.1.10--> */
//      while (pair != last_pair && pair != NULL)
      while (pair != NULL)
/* <--v.1.1.10 */
        {      
          /* Get pointers to next and previous base pairs */
          next_pair = pair->next_pair;
          prev_pair = pair->prev_pair;
         
          /* Loop through strands */
          for (i = 0; i < 2; i++)
            {
              /* Compare with each strand */
              if (pair->strand[i] != NULL)
                {
                  strand[i] = pair->strand[i];
                  k = 1;
                  
                  /* Compare with each bond_atom */
                  for (j = 0; j < 2; j++)
                    {
                      /* Check that residue name, number and chain-id match */
                      if (strand[i]->chain == bond_atom[j]->chain
                          && strand[i]->resnum == bond_atom[j]->resnum
                          && !strncmp(strand[i]->resname,
                                      bond_atom[j]->resname,3))
                        {
                          /* Store the details in the appropriate structure */
                          store_interaction(pair,next_pair,prev_pair,
                                            bond,bond_atom,i,j,k);
                        }
                      k = 0;
                    }

                  /* Record if largest number of non-phosphate bonds found */
                  if (pair->nbond[i] > *max_nbond)
                    *max_nbond = pair->nbond[i];
               
                  /* Record if largest number of phosphate bonds found */
                  if (pair->npbond[i] > *max_npbond)
                    *max_npbond = pair->npbond[i];
                }
            }
          pair = pair->next_pair;
        }
      bond = bond->next_bond;
    }

/* Return to text_interactions =======================================*/
  return(0);
}



/***********************************************************************
 *
 * store_interaction - Store the current protein-DNA interaction in the
 *                     appropriate base-pair record
 *
 **********************************************************************/
void
store_interaction(PAIRS *pair,PAIRS *next_pair,PAIRS *prev_pair,
                  BONDS *bond,ATOMS *bond_atom[2],int i,int j,int k)
{
/* Declare variables =================================================*/
  int is_phosphate;

  /* Store the details in the appropriate structure */

  /* Determine whether the interacting atom on the DNA belongs to a
     phosphate group */
  is_phosphate = is_phosphate_atom(bond_atom[j]->name);

  /* Phosphate bonds */
  if (is_phosphate == TRUE)
    {
      /* If this is any atom other than for the next base, then
         process */
      if (strncmp(bond_atom[j]->name+1, "O3", 2))
        {
          /* Check that maximum number of interactions with this
             phosphate group not exceeded */
          if (pair->npbond[i] > MAXARRAY - 1)
            {
              printf("*** WARNING. Maximum number of interactions (%d) ",
                     MAXARRAY);
              printf("exceeded for phosphate group of [%s %4d %c]\n",
                     pair->strand[i]->resname,
                     pair->strand[i]->resnum,
                     pair->strand[i]->chain);
              printf("***          Interaction between ");
              printf("atom %s and %s %s %4d %c lost\n",
                     bond_atom[j]->name,
                     bond_atom[k]->name,
                     bond_atom[k]->resname,
                     bond_atom[k]->resnum,
                     bond_atom[k]->chain);
            }
      
          /* Otherwise, store the interaction */
          else
            {
              pair->pbondto[i][pair->npbond[i]] = bond_atom[k];
              pair->pbondfrom[i][pair->npbond[i]] = bond_atom[j];
              pair->ptype[i][pair->npbond[i]] = bond->type;
              pair->npdups[i][pair->npbond[i]] = 1;
              (pair->npbond[i])++;
            }
        }
      else
        {
          if (i == ONE && next_pair != NULL)
            {
              /* Check that maximum number of interactions with this
                 phosphate group not exceeded */
              if (next_pair->npbond[i] > MAXARRAY - 1)
                {
                  printf("*** WARNING. Maximum number of interactions (%d) ",
                         MAXARRAY);
                  printf("exceeded for phosphate group of [%s %4d %c]\n",
                         next_pair->strand[i]->resname,
                         next_pair->strand[i]->resnum,
                         next_pair->strand[i]->chain);
                  printf("***          Interaction between ");
                  printf("atom %s and %s %s %4d %c lost\n",
                         bond_atom[j]->name,
                         bond_atom[k]->name,
                         bond_atom[k]->resname,
                         bond_atom[k]->resnum,
                         bond_atom[k]->chain);
                }
      
              /* Otherwise, store the interaction */
              else
                {
                  next_pair->pbondto[i][next_pair->npbond[i]] = bond_atom[k];
                  next_pair->pbondfrom[i][next_pair->npbond[i]] = bond_atom[j];
                  next_pair->ptype[i][next_pair->npbond[i]] = bond->type;
                  next_pair->npdups[i][next_pair->npbond[i]] = 1;
                  (next_pair->npbond[i])++;
                }
            }
          else if (i == TWO && prev_pair != NULL)
            {
              /* Check that maximum number of interactions with this
                 phosphate group not exceeded */
              if (prev_pair->npbond[i] > MAXARRAY - 1)
                {
                  printf("*** WARNING. Maximum number of interactions (%d) ",
                         MAXARRAY);
                  printf("exceeded for phosphate group of [%s %4d %c]\n",
                         prev_pair->strand[i]->resname,
                         prev_pair->strand[i]->resnum,
                         prev_pair->strand[i]->chain);
                  printf("***          Interaction between ");
                  printf("atom %s and %s %s %4d %c lost\n",
                         bond_atom[j]->name,
                         bond_atom[k]->name,
                         bond_atom[k]->resname,
                         bond_atom[k]->resnum,
                         bond_atom[k]->chain);
                }
      
              /* Otherwise, store the interaction */
              else
                {
                  prev_pair->pbondto[i][prev_pair->npbond[i]]
                    = bond_atom[k];
                  prev_pair->pbondfrom[i][prev_pair->npbond[i]]
                    = bond_atom[j];
                  prev_pair->ptype[i][prev_pair->npbond[i]]
                    = bond->type;
                  prev_pair->npdups[i][prev_pair->npbond[i]]
                    = 1;
                  (prev_pair->npbond[i])++;
                }
            } 
        }
    }
  
  /* Non phosphate bonds */
  else
    {
      /* Check that maximum number of interactions with this base
         not exceeded */
      if (pair->nbond[i] > MAXARRAY - 1)
        {
          printf("*** WARNING. Maximum number of interactions (%d) ",
                 MAXARRAY);
          printf("exceeded for base [%s %4d %c]\n",
                 pair->strand[i]->resname,
                 pair->strand[i]->resnum,
                 pair->strand[i]->chain);
          printf("***          Interaction between ");
          printf("atom %s and %s %s %4d %c lost\n",
                 bond_atom[j]->name,
                 bond_atom[k]->name,
                 bond_atom[k]->resname,
                 bond_atom[k]->resnum,
                 bond_atom[k]->chain);
        }
      
      /* Otherwise, store the interaction */
      else
        {
          pair->bondto[i][pair->nbond[i]] = bond_atom[k];
          pair->bondfrom[i][pair->nbond[i]] = bond_atom[j];
          pair->type[i][pair->nbond[i]] = bond->type;
          pair->ndups[i][pair->nbond[i]] = 1;
          (pair->nbond[i])++;
        }
    }
}



/***********************************************************************
 *
 * text_interactions - Function to draw the interactions to the 
 *                     Postscript file
 *
 **********************************************************************/
int
text_interactions(PAIRS *first_pair, PAIRS *last_pair, 
                  PAIRS *first_page_pair, PAIRS *last_page_pair, 
                  SIZES *size, PARAM *par,
                  int max_nbond, int max_npbond)
{
/* Declare variables =================================================*/
   int    bond_point, i;
   int    nbond[3], nslots[3];
   int    first_flag, last_flag;
   double bondx_pos[3], bondy_pos[3];
   double text_size[3], y_spacing[3];
   PAIRS  *pair, *next_pair;

/* Initialize variables ==============================================*/
   pair = first_page_pair;

/* Loop throught the pairs ===========================================*/
   while (pair != last_page_pair)
   {
      next_pair = pair->next_pair;
      
/* Calculate the positions ===========================================*/      
      for (i = 0; i < 2; i++)
      {

/* Set flags =========================================================*/
         set_pair_flag(pair, first_page_pair, last_page_pair,
                       &first_flag, &last_flag);
         
/* Count number of interactions to sugar, base, phosphate and next
   phosphate */
         count_bond_num(pair, i, nbond, nslots);

/* Initial values ====================================================*/
         init_int_pos(pair, size, par, i, nslots, y_spacing,
                      text_size, bondx_pos, bondy_pos);

/* Non phosphate bonds (ie to base and sugar) ========================*/
         bond_point = NONPHOS;
         nonphos_int_pos(pair, size, par, i, bondx_pos[bond_point],
                         bondy_pos[bond_point], y_spacing[bond_point],
                         text_size[bond_point]);

/* Draw non-phosphate bonds ==========================================*/
          draw_interaction_group(first_pair, last_pair, pair, par, size,
                                i, last_flag, bond_point,
                                 y_spacing[bond_point],
                                text_size[bond_point]);

/* Phosphate bonds ===================================================*/
         bond_point = PHOS;
         phos_int_pos(pair, size, par, first_flag, i, bondx_pos[bond_point],
                      bondy_pos[bond_point], y_spacing[bond_point],
                         text_size[bond_point]);
         
/* Draw phosphate bonds ==============================================*/
         if (((i == ONE && first_flag == FALSE) || i == TWO) &&
             pair->strand[i] != NULL)
           draw_interaction_group(first_pair, last_pair, pair, par, size,
                                  i, last_flag, bond_point,
                                  y_spacing[bond_point],
                                  text_size[bond_point]);

/* Bonds to last phosphate ===========================================*/
/* v.1.1.1--> */
         // if (pair->next_pair == last_pair && i == ONE)
         if (pair->next_pair == last_page_pair && last_page_pair != NULL
             && i == ONE)
/* <--v.1.1.1 */
           { 
            bond_point = LAST_PHOS;
             last_phos_int_pos(pair, size, par, last_flag, i, 
                               bondx_pos[bond_point], bondy_pos[bond_point],
                               y_spacing[bond_point],
                               text_size[bond_point]);

/* Draw bonds ========================================================*/
             draw_interaction_group(first_pair, last_pair, pair, par, size,
                                    i, last_flag, bond_point,
                                    y_spacing[bond_point],
                                    text_size[bond_point]);
           }

/* Draw bridging interactions ========================================== */
         text_draw_bdge_int(first_pair, last_pair, pair, next_pair, size,
                            par, i, first_flag, last_flag, text_size); 
       }
      pair = pair->next_pair;      
   }

/* Return to text_interactions =======================================*/
   return(0);
}



/***********************************************************************
 *
 * set_pair_flag - Function to return apropriate flags if the pair is 
 *                 the first or the last on the page 
 *
 **********************************************************************/
int
set_pair_flag(PAIRS *pair, PAIRS *first_pair, PAIRS *last_pair, 
              int *first_flag, int *last_flag)
{
/* Initialize flags ==================================================*/
   *first_flag = *last_flag = FALSE;

/* Set flags accordingly =============================================*/
   if (pair == first_pair && pair->prev_pair != NULL)
      *first_flag = TRUE;

   if (last_pair != NULL)
      if (pair == last_pair->prev_pair && pair->next_pair != NULL)
         *last_flag = TRUE;

   return(0);
}



/***********************************************************************
 *
 * count_bond_num - Function to count the number of each type of bond
 *                  attached the the 3 bond_points for this base
 *                  on this strand
 *
 **********************************************************************/
int
count_bond_num(PAIRS *pair, int i, int nbond[3], int nslots[3])
{
/* Declare variables =================================================*/
  int   bond_point, j, k, nbridge;
  PAIRS *next_pair;
      
/* Initialize variables ==============================================*/
  for (j = 0; j < 3; j++)
    {
      nbond[j] = 0;
      nslots[j] = 0;
    }
  next_pair = pair->next_pair;
   
/* Counts for interactions to base and sugar =========================*/
  bond_point = NONPHOS;
  for (j = 0; j < pair->nbond[i]; j++)
    {
      /* Increment count of direct bonds */
      nbond[bond_point]++;

      /* Check whether there are any bridging interactions here */
      if (pair->bridge_flag[i][j] == TRUE)
        {
          /* Count the number of interactions */
          nbridge = pair->nbridge[i][j];

          /* Store the slot positions for each interaction */
          for (k = 0; k < pair->nbridge[i][j]; k++)
            pair->bdgey_pos[i][j][k] = nslots[bond_point] + k;

          /* Calculate position for the direct interaction */
          pair->bondy_pos[i][j] = nslots[bond_point]
            + (float) (nbridge - 1) / 2.0;

          /* Increment count of available slots */
          nslots[bond_point] = nslots[bond_point] + nbridge;
        }

      /* If no bridging interactions, then place interaction in next
         slot */
      else
        {
          /* Get the relative position for this interaction */
          pair->bondy_pos[i][j] = nslots[bond_point];
          nslots[bond_point]++;
        }
    }

  /* Repeat for interactions to phosphate  ===========================*/
  bond_point = PHOS;
  for (j = 0; j < pair->npbond[i]; j++)
    {
      /* Increment count of direct bonds */
      nbond[bond_point]++;

      /* Check whether there are any bridging interactions here */
      if (pair->pbridge_flag[i][j] == TRUE)
        {
          /* Count the number of interactions */
          nbridge = pair->npbridge[i][j];

          /* Store the slot positions for each interaction */
          for (k = 0; k < pair->npbridge[i][j]; k++)
            pair->pbdgey_pos[i][j][k] = nslots[bond_point] + k;

          /* Calculate position for the direct interaction */
          pair->pbondy_pos[i][j] = nslots[bond_point]
            + (float) (nbridge - 1) / 2.0;

          /* Increment count of available slots */
          nslots[bond_point] = nslots[bond_point] + nbridge;
        }

      /* If no bridging interactions, then place interaction in next
         slot */
      else
        {
          /* Get the relative position for this interaction */
          pair->pbondy_pos[i][j] = nslots[bond_point];
          nslots[bond_point]++;
        }
    }

  /* Repeat for interactions to last phosphate on page ===========*/
  bond_point = LAST_PHOS;
  if (next_pair != NULL)
    for (j = 0; j < next_pair->npbond[i]; j++)
      {
        /* Increment count of direct bonds */
        nbond[bond_point]++;

        /* Check whether there are any bridging interactions here */
/* v.1.1.1--> */
        // if (next_pair->bridge_flag[i][j] == TRUE)
        if (next_pair->pbridge_flag[i][j] == TRUE)
/* <--v.1.1.1 */
          {
            /* Count the number of interactions */
/* v.1.1.1--> */
            // nbridge = next_pair->nbridge[i][j];
            nbridge = next_pair->npbridge[i][j];
/* <--v.1.1.1 */

            /* Store the slot positions for each interaction */
/* v.1.1.1--> */
            // for (k = 0; k < next_pair->nbridge[i][j]; k++)
            //   next_pair->bdgey_pos[i][j][k] = nslots[bond_point] + k;
            for (k = 0; k < next_pair->npbridge[i][j]; k++)
              next_pair->pbdgey_pos[i][j][k] = nslots[bond_point] + k;
/* <--v.1.1.1 */

            /* Calculate position for the direct interaction */
/* v.1.1.1--> */
            // next_pair->bondy_pos[i][j] = nslots[bond_point]
            //  + (float) (nbridge - 1) / 2.0;
            next_pair->pbondy_pos[i][j] = nslots[bond_point]
              + (float) (nbridge - 1) / 2.0;
/* <--v.1.1.1 */

            /* Increment count of available slots */
            nslots[bond_point] = nslots[bond_point] + nbridge;
          }

        /* If no bridging interactions, then place interaction in next
           slot */
        else
          {
            /* Get the relative position for this interaction */
/* v.1.1.1--> */
            // next_pair->bondy_pos[i][j] = nslots[bond_point];
            next_pair->pbondy_pos[i][j] = nslots[bond_point];
/* <--v.1.1.1 */
            nslots[bond_point]++;
          }
      }

/* Return ============================================================*/
  return(0);
}


/***********************************************************************
 *
 * init_int_pos - Function to calculate interaction position starting 
 *                points for each of the 3 different bond-point types
 *
 **********************************************************************/
int
init_int_pos(PAIRS *pair, SIZES *size, PARAM *par, int i,
             int nslots[3], double y_spacing[3], double text_size[3],
             double bondx_pos[3], double bondy_pos[3])
{
/* Declare variables =================================================*/
  int bond_point, side;
  double x_diff;
  float  rel_pos[3], text_height;
  PAIRS  *next_pair;
   
/* Initialize variables ==============================================*/
  next_pair = pair->next_pair;
  for (bond_point = 0; bond_point < 3; bond_point++)
    {
      y_spacing[bond_point] = 0.0;
      bondx_pos[bond_point] = 0.0;
      bondy_pos[bond_point] = 0.0;
    }

  /* Get the appropriate side (left = -ve, right = +ve) */
  side = i * 2 - 1;

  /* Get distance from sugar/phosphate to interacting group */
  x_diff = SUGAR_TO_BASE_WIDTH + WATER_TO_P_WIDTH;

/* Calculate the relative mid-points for each bond-group =============*/
  /* Interactions to sugar/base are the defining level */
  rel_pos[0] = 0.0;

  /* Interactions to phosphate are higher by half the phosphate-to-
     phosphate height */
  rel_pos[1] = - (float) (side * (P_TO_P_HEIGHT * size->txtscale) / 2.0);
      
  /* Interactions to next phosphate are lower by half the phosphate-
     to-phosphate height */
  rel_pos[2] = (float) (side * (P_TO_P_HEIGHT * size->txtscale) / 2.0);

/* Check whether the text-sizes for the residue names need to be
   reduced at all and get the starting y-positions ===================*/
  /* Loop over the bond-groups (ie interactions to: the sugar/base,
     the phospate and the next phosphate) */
  for (bond_point = 0; bond_point < 3; bond_point++)
    {
      /* Calculate the overall text height */
      text_size[bond_point] = RESNAME_TXT;
      text_height = (float) (text_size[bond_point] * (nslots[bond_point] - 1));

      /* If interactions will take up too much space, reduce the text
         size */
      if (text_height > 0.5 * P_TO_P_HEIGHT)
        text_size[bond_point]
          = (0.5 * P_TO_P_HEIGHT) / (nslots[bond_point] - 1);

      /* Calculate the y-spacing according to the text-size */
      y_spacing[bond_point] = text_size[bond_point] * size->txtscale;

      /* Calculate the starting y-position for this group, given the
         y-spacing and number of interactions */
      bondy_pos[bond_point] = pair->y_pos[i] + rel_pos[bond_point]
        + (y_spacing[bond_point] * (nslots[bond_point] - 1)) / 2.0;
    }

/* Starting x-positions ==============================================*/
  /* Bonds to sugar or base */
  bondx_pos[0] = pair->x_pos[i]
    + side * (SUGAR_TO_BASE_WIDTH + WATER_TO_SUGAR_WIDTH) * size->txtscale;

  /* Bonds to phosphate */
  bondx_pos[1] = pair->x_pos[i]
    + side * (SUGAR_TO_BASE_WIDTH + P_TO_SUGAR_WIDTH + WATER_TO_P_WIDTH)
      * size->txtscale;

  /* Bonds to next phosphate */
  bondx_pos[2] = bondx_pos[1];

/* Return to interact_pos ============================================*/   
   return(0);
}



/***********************************************************************
 *
 * nonphos_int_pos - Function to determine the non-phosphate bond 
 *                   positions 
 *
 **********************************************************************/
int
nonphos_int_pos(PAIRS *pair, SIZES *size, PARAM *par, int i,
                double bondx_pos,double bondy_pos, double y_spacing,
                double text_size)
{
/* Declare variables =================================================*/
  int j, k, side;
  double bridgex_pos, bridging_dist; 

  /* Get the appropriate side (left = -ve, right = +ve) */
  side = i * 2 - 1;

/* Calculate the bridging distance from the water position to the
   interacting residue position */
  bridging_dist = calc_bridge_dist(text_size, par);

  /* Loop over all the interactions with the current sugar/base unit */
  for (j = 0; j < pair->nbond[i]; j++)
    {
      /* Store the y-position for this interaction */
      pair->bondy_pos[i][j] = bondy_pos - pair->bondy_pos[i][j] * y_spacing;

      /* Store the x-position for this interaction */
      pair->bondx_pos[i][j] = bondx_pos;

      /* Check for any bridging interactions */
      if (pair->bridge_flag[i][j] == TRUE)
        {
          /* Get the x-coordinate for this interaction */
          bridgex_pos = bondx_pos + side * bridging_dist * size->txtscale;

          /* Loop over all the bridging positions and calculate x- and
             y-coords */
          for (k = 0; k < pair->nbridge[i][j]; k++)
            {
              /* X-pos */
              pair->bdgex_pos[i][j][k] = bridgex_pos;

              /* Y-pos */
              pair->bdgey_pos[i][j][k] = bondy_pos
                - pair->bdgey_pos[i][j][k] * y_spacing;
            }  
        }
    }

/* Return to interact_pos ============================================*/
  return(0);
}



/***********************************************************************
 *
 * phos_int_pos - Function to determine the phosphate bond positions 
 *
 **********************************************************************/
int
phos_int_pos(PAIRS *pair, SIZES *size, PARAM *par, int first_flag, int i,
             double bondx_pos, double bondy_pos, double y_spacing,
             double text_size)
{
/* Declare variables =================================================*/
  int j, k, side;
  double bridgex_pos, bridging_dist;

  /* Get the appropriate side (left = -ve, right = +ve) */
  side = i * 2 - 1;

/* Calculate the bridging distance from the water position to the
   interacting residue position */
  bridging_dist = calc_bridge_dist(text_size, par);

  if ((i == ONE && first_flag == FALSE) || (i == TWO))
    for (j = 0; j < pair->npbond[i]; j++)
      {
        /* Y-pos */
        pair->pbondy_pos[i][j]
          = bondy_pos - pair->pbondy_pos[i][j] * y_spacing;
            
        /* X-pos */
        pair->bondx_pos[i][j] = bondx_pos;

        /* Check for any bridging interactions */
        if (pair->pbridge_flag[i][j] == TRUE)
          {
            /* Get the x-coordinate for this interaction */
            bridgex_pos = bondx_pos + side * bridging_dist * size->txtscale;

            /* Loop over all the bridging positions and calculate x- and
               y-coords */
            for (k = 0; k < pair->npbridge[i][j]; k++)
              {
                /* X-pos */
                pair->pbdgex_pos[i][j][k] = bridgex_pos;

                /* Y-pos */
                pair->pbdgey_pos[i][j][k] = bondy_pos
                  - pair->pbdgey_pos[i][j][k] * y_spacing;
              }  
          }
      }

/* Return to interact_pos ============================================*/
   return(0);
}



/***********************************************************************
 *
 * last_phos_int_pos - Function to determine the phosphate bond positions
 *                     for next pair if last pair on the page 
 *
 **********************************************************************/
int
last_phos_int_pos(PAIRS *pair, SIZES *size, PARAM *par, int last_flag, int i, 
                  double bondx_pos, double bondy_pos, double y_spacing,
                  double text_size)
{
/* Declare variables =================================================*/
  int j, k, side;
  double bridgex_pos, bridging_dist;
  PAIRS  *next_pair;

  /* Get the appropriate side (left = -ve, right = +ve) */
  side = i * 2 - 1;

/* Calculate the bridging distance from the water position to the
   interacting residue position */
  bridging_dist = calc_bridge_dist(text_size, par);

  if (pair != NULL)
    {
      /* Get the next pair */
      next_pair = pair->next_pair;

      if (next_pair != NULL)
        {
          for (j = 0; j < next_pair->npbond[i]; j++)
            {
/* v.1.1.1--> */
              // pair->bondx_pos[i][j] = bondx_pos;
              next_pair->bondx_pos[i][j] = bondx_pos;
/* <--v.1.1.1 */
              next_pair->pbondy_pos[i][j]
                = bondy_pos - next_pair->pbondy_pos[i][j] * y_spacing;;

              /* Check for any bridging interactions */
/* v.1.1.1--> */
              // if (next_pair->bridge_flag[i][j] == TRUE)
              if (next_pair->pbridge_flag[i][j] == TRUE)
/* <--v.1.1.1 */
                {
                  /* Get the x-coordinate for this interaction */
                  bridgex_pos = bondx_pos
                    + side * bridging_dist * size->txtscale;

                  /* Loop over all the bridging positions and calculate x-
                     and y-coords */
/* v.1.1.1--> */
                  // for (k = 0; k < next_pair->nbridge[i][j]; k++)
                  for (k = 0; k < next_pair->npbridge[i][j]; k++)
/* <--v.1.1.1 */
                    {
                      /* X-pos */
/* v.1.1.1--> */
                      // next_pair->bdgex_pos[i][j][k] = bridgex_pos;
                      next_pair->pbdgex_pos[i][j][k] = bridgex_pos;
/* <--v.1.1.1 */

                      /* Y-pos */
/* v.1.1.1--> */
                      // next_pair->bdgey_pos[i][j][k] = bondy_pos
                      //  - next_pair->bdgey_pos[i][j][k] * y_spacing;
                      next_pair->pbdgey_pos[i][j][k] = bondy_pos
                        - next_pair->pbdgey_pos[i][j][k] * y_spacing;
/* <--v.1.1.1 */
                    }  
                }
            }
        }
    }

/* Return to interact_pos ============================================*/
   return(0);
}



/***********************************************************************
 *
 * calc_bridge_dist - Function to calculate the bridging distance given
 *                    the text size
 *
 **********************************************************************/
double
calc_bridge_dist(double text_size, PARAM *par)
{
/* Declare variables =================================================*/
  double bridging_dist, resname_interact_width, text_width, water_no_width; 

/* Calculate the bridging distance from the water position to the
   interacting residue position */
  text_width = 0.71 * text_size;
  resname_interact_width = 3 * text_width;
  water_no_width = 3 * text_width;
  if (par->landscape == TRUE)
    water_no_width = 1.5 * text_size;
  bridging_dist = 2.5 * WATER_RADIUS + water_no_width
    + resname_interact_width;

/* Return bridging distance ==========================================*/   
   return(bridging_dist);
}



/***********************************************************************
 *
 * draw_interaction_group - Function to draw the interactions to the given
 *                          object: sugar/base, phosphate, or next
 *                          phosphate
 *
 **********************************************************************/
int
draw_interaction_group(PAIRS *first_pair, PAIRS *last_pair, PAIRS *pair,
                       PARAM *par, SIZES *size, int i, int last_flag,
                       int bond_point, double y_spacing, double text_size)
{
/* Declare variables =================================================*/
  int    ibond, nbonds, side, wbond;
  int    on_left, water_flag;
  int    bond_type;
  double x1, y1;
  double x_init, y_init;
  double cx, cy;
  double text_width;
  char   int_det[LINE_LEN], atom_name[5], atom_colour[COLNAME_LENGTH + 1];
  char   atom_name_combine[11], atom_colour_combine[COLNAME_LENGTH + 1];
  ATOMS  *interaction, *int_combine;
  PAIRS  *next_pair;

/* Initialize variables ==============================================*/
  cx     = size->x_centre;
  cy     = size->y_centre;
  next_pair = pair->next_pair;
  text_width = 0.71 * text_size;
  wbond = 0;
  y_init = 0.0;

  /* Get the appropriate side (left = -ve, right = +ve) */
  side = i * 2 - 1;

  /* Determine the number of bonds to be plotted */
  if (bond_point == NONPHOS)
    nbonds = pair->nbond[i];
  else if (bond_point == PHOS)
    nbonds = pair->npbond[i];
  else
    {
      /* If there is no next phosphate, then return */
      if (next_pair == NULL)
        return(0);

      /* Otherwise, get the number of interactions */
      else
        nbonds = next_pair->npbond[i];
    }

/* Draw ==============================================================*/
  for (ibond = 0; ibond < nbonds; ibond++)
    {
      water_flag = FALSE;

      /* Get the pointer to the atom involved in this interaction */
      if (bond_point == NONPHOS)
        {
          interaction = pair->bondto[i][ibond];
          int_combine = pair->combine[i][ibond];
          bond_type = pair->type[i][ibond];
        }
      else if (bond_point == PHOS)
        {
          interaction = pair->pbondto[i][ibond];
          int_combine = pair->pcombine[i][ibond];
          bond_type = pair->ptype[i][ibond];
        }
      else
        {
          interaction = next_pair->pbondto[i][ibond];
          int_combine = next_pair->pcombine[i][ibond];
          bond_type = next_pair->ptype[i][ibond];
        }

      /* get interaction details */
      get_interact_det(first_pair, last_pair, pair, par, interaction,
                       int_combine, int_det, i, bond_point, atom_name,
                       atom_colour, atom_name_combine, atom_colour_combine,
                       &water_flag, bond_type);

      /* Write to postscript file */
      x_init = pair->bondx_pos[i][ibond];
/* v.1.1.2--> */
      if (bond_point == LAST_PHOS && next_pair != NULL)
        x_init = next_pair->bondx_pos[i][ibond];
/* <--v.1.1.2 */

      /* If not a water, shift name slightly */
      if (water_flag == FALSE)
        x_init = x_init + side * (text_width * size->txtscale) / 2.0;

      /* Get appropriate y-position */
      if (bond_point == NONPHOS)
        y_init = pair->bondy_pos[i][ibond];
      else if (bond_point == PHOS)
        y_init = pair->pbondy_pos[i][ibond];
      else if (bond_point == LAST_PHOS && pair->next_pair != NULL)
        y_init = next_pair->pbondy_pos[i][ibond];
      
      /* Transform for landscape mode */
      pstrans_(x_init, y_init, cx, cy, &x1, &y1, par);
      
      /* Draw water interaction */
      if (water_flag == TRUE)
        {
          pswatertxt_(size, par, x_init, y_init, int_det, i, text_size,
                      wbond);
          wbond++;
        }

      /* Draw protein interaction */
      else if (water_flag == FALSE)
        {
          /* Reset water-counter */
          wbond = 0;

          /* Determine whether doing left- or right-hand strand */
          if ((par->landscape == FALSE && i == ONE) ||
              (par->landscape == TRUE && i == TWO))
            on_left = TRUE;
          else
            on_left = FALSE;

          /* On left-hand strand */
          if (on_left == TRUE)
            {
              /* Print full atom and residue details */
              if (par->landscape == TRUE)
                psrotx_(x1, y1, (int) (text_size * size->txtscale), int_det);
              else
                psrtxt_(x1, y1, (int) (text_size * size->txtscale), int_det);
            
              /* If have two atom names combined, print both in the
                 colour of the left-hand atom */
              if (atom_name_combine[0] != '\0')
                {
                  pssave_();
                  pscolb_(atom_colour_combine);
                  if (par->landscape == TRUE)
                    psrotx_(x1, y1, (int) (text_size * size->txtscale),
                            atom_name_combine);
                  else
                    psrtxt_(x1, y1, (int) (text_size * size->txtscale),
                            atom_name_combine);
                  psrest_();
                }
            
              /* Print just the atom name in the given colour */
              if (atom_name[0] != '\0')
                {
                  pssave_();
                  pscolb_(atom_colour);
                  if (par->landscape == TRUE)
                    psrotx_(x1, y1, (int) (text_size * size->txtscale), atom_name);
                  else
                    psrtxt_(x1, y1, (int) (text_size * size->txtscale), atom_name);
                  psrest_();
                }
            }

          /* Right-hand strand */
          else
            {
              /* Print full atom and residue details */
              if (par->landscape == TRUE)
                psrltx_(x1, y1, (int) (text_size * size->txtscale), int_det);
              else
                psltxt_(x1, y1, (int) (text_size * size->txtscale), int_det);
              
              /* If have two atom names combined, print both in the
                 colour of the right-hand atom */
              if (atom_name_combine[0] != '\0')
                {
                  pssave_();
                  pscolb_(atom_colour_combine);
                  if (par->landscape == TRUE)
                    psrltx_(x1, y1, (int) (text_size * size->txtscale),
                            atom_name_combine);
                  else
                    psltxt_(x1, y1, (int) (text_size * size->txtscale),
                            atom_name_combine);
                  psrest_();
                }
              
              /* Print just the atom name in the given colour */
              if (atom_name[0] != '\0')
                {
                  pssave_();
                  pscolb_(atom_colour);
                  if (par->landscape == TRUE)
                    psrltx_(x1, y1, (int) (text_size * size->txtscale), atom_name);
                  else
                    psltxt_(x1, y1, (int) (text_size * size->txtscale), atom_name);
                  psrest_();
                }
            }
        }
      
      /* Draw bond line(s) for this interaction */
      draw_bond_line(pair, par, size, i, ibond, bond_point, y_spacing);
    }

/* Return ===========================================================*/
   return(0);
}

/***********************************************************************
 *
 * get_interact_det - Function to get the interaction information 
 *                    drawing
 *
 **********************************************************************/
int
get_interact_det(PAIRS *first_pair, PAIRS *last_pair, PAIRS *pair,
                 PARAM *par, ATOMS *interaction, ATOMS *int_combine,
                 char string[LINE_LEN], int i, int bond_point,
                 char atom_name[5], char atom_colour[COLNAME_LENGTH + 1],
                 char atom_name_combine[11],
                 char atom_colour_combine[COLNAME_LENGTH + 1],
                 int *water_flag, int bond_type)
{
/* Declare variables =================================================*/
  char asterisk[2], res_name[4], temp_atom_name[5], temp_name[30];
  char swap_colour[COLNAME_LENGTH + 1];
  int flag, on_left, is_double;

/* Initialize variables ==============================================*/
  flag = FALSE;
  atom_name_combine[0] = '\0';
  atom_colour_combine[0] = '\0';

/* Check whether this residue appears anywhere else on the plot */
  is_double = check_for_double(first_pair, last_pair, pair, interaction,
                               i, bond_point);
  if (is_double == TRUE)
    strcpy(asterisk,"*");
  else
    asterisk[0] = '\0';

/* Get information ===================================================*/
  /* If this is a water, then set water flag */
  if (!strncmp(interaction->resname, "HOH", 3))
    {
      *water_flag = TRUE;
      sprintf(string, "%d", interaction->resnum);
      if (interaction->chain != ' ')
        {
          sprintf(temp_name,"(%c)", interaction->chain);
          strcat(string,temp_name);
        }
      strcat(string,asterisk);
    }

  /* Otherwise, for a residue name, string the residue and atom
     details together */
  else
    {
      /* Determine whether atom-name to be on left or right of output
         string */
      if ((par->landscape == FALSE && i == ONE) ||
          (par->landscape == TRUE && i == TWO))
        on_left = FALSE;
      else
        on_left = TRUE;

      /* Get 1-letter code for this residue */
      match_res_code(par, interaction->resname, res_name);

      /* Get the atom name and colour */
      get_atom_name(par, interaction, atom_name, atom_colour);

      /* If there is another atom combined with the current one, get
         its name and colour */
      if (int_combine != NULL)
        {
          /* Get combined atom name */
          get_atom_name(par, int_combine, temp_atom_name,
                        atom_colour_combine);

          /* Combine the two names */
          if (on_left == TRUE)
            {
              /* For left-hand strand put first atom name second */
              strcpy(atom_name_combine, temp_atom_name);
              strcat(atom_name_combine, ", ");
              strcat(atom_name_combine, atom_name);
              strncpy(atom_name,temp_atom_name,4);
              atom_name[4] = '\0';
            }
          else
            {
              /* For right-hand strand put first atom name first */
              strcpy(atom_name_combine, atom_name);
              strcat(atom_name_combine, ", ");
              strcat(atom_name_combine, temp_atom_name);
              strncpy(atom_name,temp_atom_name,4);
              atom_name[4] = '\0';
            }
        }

      /* If atom name is on the left arrange output string accordingly */
      if (on_left == TRUE)
        {
          /* Start with atom name (unless only showing residue name) */
          if (par->nonbond_by_residue == TRUE && bond_type == NBOND)
            {
              string[0] = '\0';
              atom_name[0] = '\0';
            }
          else
            {
              /* Use combined name where relevant */
              if (atom_name_combine[0] != '\0')
                {
                  strcpy(string,atom_name_combine);
                  strcat(atom_name,",");

                  /* Swap the atom colours */
                  strcpy(swap_colour,atom_colour_combine);
                  strcpy(atom_colour_combine,atom_colour);
                  strcpy(atom_colour,swap_colour);
                }
              else
                strcpy(string,atom_name);
              strcat(string," ");
            }

          /* Form residue name and number */
          sprintf(temp_name, "%s%d",res_name,interaction->resnum); 
          strcat(string,temp_name);

          /* Add chain-id if not blank */
          if (interaction->chain != ' ')
            {
              sprintf(temp_name,"(%c)", interaction->chain);
              strcat(temp_name,asterisk);
              strcat(temp_name," ");
              strcat(string,temp_name);
            }

          /* Otherwise, just add an asterisk if relevant */
          else
            { 
              if (asterisk[0] != '\0')
                strcat(string,"*");
              strcat(string," ");
            }
        }

      /* Otherwise, put atom name on the right */
      else
        {
          /* Form the output string */
          sprintf(string, "%s%d",res_name,interaction->resnum); 

          /* Add chain-id if not blank */
          if (interaction->chain != ' ')
            {
              sprintf(temp_name,"(%c)", interaction->chain);
              strcat(temp_name,asterisk);
              strcat(string,temp_name);
            }

          /* Otherwise, just add an asterisk if relevant */
          else
            { 
              if (asterisk[0] != '\0')
                strcat(string,"*");
            }

          /* Add atom name (unless only showing residue name) */
          if (par->nonbond_by_residue == TRUE && bond_type == NBOND)
            atom_name[0] = '\0';
          else
            {
              strcat(string," ");

              /* Use combined name where relevant */
              if (atom_name_combine[0] != '\0')
                {
                  strcat(string,atom_name_combine);

                  /* Swap the atom colours */
                  strcpy(swap_colour,atom_colour_combine);
                  strcpy(atom_colour_combine,atom_colour);
                  strcpy(atom_colour,swap_colour);
                }
              else
                strcat(string,atom_name);
            }
        }
   }
            
/* Return to text_draw_interaction ===================================*/
   return(0);
}



/***********************************************************************
 *
 * get_atom_name - Function to get the name and colour for the given
 *                 atom
 *
 **********************************************************************/
void
get_atom_name(PARAM *par, ATOMS *atom_ptr, char atom_name[5],
              char atom_colour[COLNAME_LENGTH + 1])
{
/* Declare variables =================================================*/
  char temp_name[20];
  int ifirst, ilast, j, len;

  /* Determine atom-type to get appropriate colour */
  strncpy(atom_name,atom_ptr->name,4);
  atom_name[4] = '\0';

  /* Nitrogen */
  if (atom_name[1] == 'N')
    strcpy(atom_colour, par->nitrogen_colour);
  else if (atom_name[1] == 'O')
    strcpy(atom_colour, par->oxygen_colour);
  else if (atom_name[1] == 'S')
    strcpy(atom_colour, par->sulphur_colour);
  else if (atom_name[1] == 'C')
    strcpy(atom_colour, par->carbon_colour);
  else
    strcpy(atom_colour, par->unknown_atom_colour);

  /* Remove trailing blanks from atom name */
  ifirst = -1;
  ilast = 0;
  for (j = 0; j < 4; j++)
    if (atom_name[j] != ' ')
      {
        if (ifirst < 0)
          ifirst = j;
        ilast = j;
      }
  if (ifirst != -1)
    {
      len = ilast - ifirst + 1;
      strncpy(temp_name,atom_name+ifirst,len);
      temp_name[len] = '\0';
      strcpy(atom_name,temp_name);
    }
}



/***********************************************************************
 *
 * check_for_double - Function to check whether the given residue appears
 *                    anywhere else on the plot
 *
 **********************************************************************/
int
check_for_double(PAIRS *first_pair, PAIRS *last_pair, PAIRS *current_pair,
                 ATOMS *current_interaction, int current_strand,
                 int current_bond_point)
{
  int in_current, is_double;
  int i, ibond, jbond, nbonds;
  ATOMS *interaction;
  PAIRS  *pair;

  /* Initilaise variables */
  is_double = FALSE;
  pair = first_pair;

  /* Loop over all the base pairs */
/* v.1.1.10--> */
//  while (pair != last_pair && is_double == FALSE)
  while (pair != NULL && is_double == FALSE)
/* <--v.1.1.10 */
    {
      /* Loop over the two strands */
      for (i = 0; i < 2 && is_double == FALSE; i++)
        {
          /* Check whether we are in the current interaction group */
          in_current = FALSE;
          if (pair == current_pair && i == current_strand &&
              current_bond_point == NONPHOS)
            in_current = TRUE;

          /* Process if not in the current interaction group */
          if (in_current == FALSE)
            {
              /* Get number of bond to the base/sugar */
              nbonds = pair->nbond[i];

              /* Loop over all the appropriate bonds */
              for (ibond = 0; ibond < nbonds && is_double == FALSE; ibond++)
                {
                  /* Get the pointer to the atom involved in
                     this interaction */
                  interaction = pair->bondto[i][ibond];

                  /* Check for match on residue number and chain-id */
                  if (interaction->resnum == current_interaction->resnum &&
                      !strncmp(interaction->resname,
                               current_interaction->resname,3) &&
                      interaction->chain == current_interaction->chain)
                    is_double = TRUE;

                  /* Check for any bridging interactions */
                  for (jbond = 0; jbond < pair->nbridge[i][ibond]; jbond++)
                    {
                      /* Get the pointer to the atom involved in
                         this interaction */
                      interaction = pair->bridgeto[i][ibond][jbond];

                      /* Check for match on residue number and chain-id */
                      if (interaction->resnum == current_interaction->resnum &&
                          !strncmp(interaction->resname,
                                   current_interaction->resname,3) &&
                          interaction->chain == current_interaction->chain)
                        is_double = TRUE;
                    }
                }

            }

          /* Check whether we are in the current interaction group */
          in_current = FALSE;
/* v.1.1.2--> */
//          if (pair == current_pair && i == current_strand &&
//              (current_bond_point == PHOS ||
//                current_bond_point == LAST_PHOS))
          if ((pair == current_pair && i == current_strand &&
               current_bond_point == PHOS) ||
              (pair == current_pair->next_pair && i == current_strand &&
               current_bond_point == LAST_PHOS))
/* <--v.1.1.2 */
            in_current = TRUE;

          /* Process if not in the current interaction group */
          if (in_current == FALSE)
            {
              /* Get number of bond to the base/sugar */
              nbonds = pair->npbond[i];

              /* Loop over all the appropriate bonds */
              for (ibond = 0; ibond < nbonds && is_double == FALSE; ibond++)
                {
                  /* Get the pointer to the atom involved in
                     this interaction */
                  interaction = pair->pbondto[i][ibond];

                  /* Check for match on residue number and chain-id */
                  if (interaction->resnum == current_interaction->resnum &&
                      !strncmp(interaction->resname,
                               current_interaction->resname,3) &&
                      interaction->chain == current_interaction->chain)
                    is_double = TRUE;

                  /* Check for any bridging interactions */
                  for (jbond = 0; jbond < pair->npbridge[i][ibond]; jbond++)
                    {
                      /* Get the pointer to the atom involved in
                         this interaction */
                      interaction = pair->pbridgeto[i][ibond][jbond];

                      /* Check for match on residue number and chain-id */
                      if (interaction->resnum == current_interaction->resnum &&
                          !strncmp(interaction->resname,
                                   current_interaction->resname,3) &&
                          interaction->chain == current_interaction->chain)
                        is_double = TRUE;
                    }
                }
            }
        }

      /* Get the next base-pair */
      pair = pair->next_pair;
    }

  /* Return flag indicating whether residue found elsewhere */
  return(is_double);
}


/***********************************************************************
 *
 * draw_bond_line - Function to draw bond line (and any duplicates) for
 *                  the given interactions
 *
 **********************************************************************/
int
draw_bond_line(PAIRS *pair, PARAM *par, SIZES *size, int i, int ibond,
               int bond_point, double y_spacing)
{
/* Declare variables =================================================*/
   int    k, nbond, side;
   int    nduplicates;
   int    line_type[MAXARRAY];
   int    bond_type[MAXARRAY];
   int    is_sugar;
   double x_start, y_start, x_end, y_end;
   double y1[MAXARRAY], y2[MAXARRAY], x1[MAXARRAY], x2[MAXARRAY];
   double y1_init;
   double cx, cy;
   double bond_gap, drop_dist, rel_pos, stop_dist;
   ATOMS  *atom, *other_atom;
   PAIRS  *next_pair;

/* Initialize variables ==============================================*/
   nbond = 0;
   bond_gap = 0.0;
   cx = size->x_centre;
   cy = size->y_centre;
   next_pair = pair->next_pair;
   y1_init = 0.0;

   /* Get the appropriate side (left = -ve, right = +ve) */
   side = i * 2 - 1;

   /* Get the relative y-positions according to the where the bond-group
      bonds */
   if (bond_point == NONPHOS)
     rel_pos = 0.0;
   else if (bond_point == PHOS)
     rel_pos = - side * (P_TO_P_HEIGHT * size->txtscale) / 2.0;
   else
     rel_pos = side * (P_TO_P_HEIGHT * size->txtscale) / 2.0;

   /* Adjustments for bond-lines to the base (the drop-distance from
      the centre-point and the stopping distance short of the base) */
   drop_dist = 0.3 * y_spacing;
   stop_dist = 0.4 * BASENAME_TXT;

   /* Get the number of duplicates for this interaction */
   if (bond_point == NONPHOS)
     nduplicates = pair->ndups[i][ibond];
   else if (bond_point == PHOS)
     nduplicates = pair->npdups[i][ibond];
   else
/* v.1.1.2--> */
//     nduplicates = pair->npdups[i][ibond];
     nduplicates = next_pair->npdups[i][ibond];
/* <--v.1.1.2 */

   /* Calculate the starting y-position */
   if (bond_point == NONPHOS)
     y1_init = pair->bondy_pos[i][ibond];
   else if (bond_point == PHOS)
     y1_init = pair->pbondy_pos[i][ibond];
   else if (bond_point == LAST_PHOS && next_pair != NULL)
     y1_init = next_pair->pbondy_pos[i][ibond];

   /* Calculate the bond-gap according to the number of duplicates */
   if (nduplicates > 1)
     {
       /* Calculate spacing between duplicate bonds */
       bond_gap = y_spacing / (nduplicates + 1);

       /* Adjust starting y-position */
       y1_init = y1_init + ((nduplicates - 1) * bond_gap) / 2.0;
     }

   /* Loop over all the duplicates */
   for (k = 0; k < nduplicates; k++)
     {
       /* Get the line type depending on whether it is to a base or
          backbone */
       if (bond_point == NONPHOS)
         {
           /* Get the pointer to this interaction */
           atom = pair->bondfrom[i][ibond];
           other_atom = pair->bondto[i][ibond];

           /* Check whether interaction is to sugar or backbone */
           is_sugar = is_sugar_atom(atom->name);
           if (is_sugar == TRUE)
             line_type[nbond] = BACKBONE;
           else 
             line_type[nbond] = BASE;
         }
       else
         line_type[nbond] = PHOSPHATE;

       /* Get the y-positions for this bond */
       y1[nbond] = y1_init;
       y2[nbond] = pair->y_pos[i] + rel_pos;

       /* Increment the y-position */
       y1_init = y1_init - bond_gap;
            
       /* Get the x-positions for this bond */
       x1[nbond] = pair->bondx_pos[i][ibond];
/* v.1.1.2--> */
       if (bond_point == LAST_PHOS && next_pair != NULL)
         x1[nbond] = next_pair->bondx_pos[i][ibond];
/* <--v.1.1.2 */
       if (line_type[nbond] == BASE)
         x2[nbond] = pair->x_numpos[i] + side * stop_dist * size->txtscale;
       else if (line_type[nbond] == BACKBONE)
         x2[nbond] = pair->x_numpos[i]
           + side * SUGAR_TO_BASE_WIDTH * size->txtscale;
       else if (line_type[nbond] == PHOSPHATE)
         x2[nbond] = pair->x_numpos[i]
           + side * (SUGAR_TO_BASE_WIDTH + P_TO_SUGAR_WIDTH) * size->txtscale;

       /* For interactions with phosphates, check whether the
          interaction is with an end-oxygen */
       if (line_type[nbond] == PHOSPHATE)
         {
           /* If this is either the first or last base of the chain, then
              need to check which phosphate atom is involved */
           if (pair->first_base[i] == TRUE || pair->last_base[i] == TRUE)
             {
               /* Get the pointers to the two atoms involved in this
                  interaction */
               atom = NULL;
               other_atom = NULL;
               if (bond_point == PHOS)
                 {
                   atom = pair->pbondfrom[i][ibond];
                   other_atom = pair->pbondto[i][ibond];
                 }
               else if (bond_point == LAST_PHOS && pair->next_pair != NULL)
                 {
                   atom = next_pair->pbondfrom[i][ibond];
                   other_atom = next_pair->pbondto[i][ibond];
                 }

               /* Proceed only if both atoms are valid */
               if (atom != NULL && other_atom != NULL)
                 {
                   /* If this is the first base, and atom is an O5*,
                      then have an end */
                   if (pair->first_base[i] == TRUE &&
                       !strncmp(atom->name+1,"O5",2))
                     {
                       /* If this is first strand, then move y downwards by
                          phosphate radius, otherwise move it up */
                       y2[nbond] = y2[nbond]
                         + side * P_RADIUS * size->txtscale;
                     }

                   /* If this is the last base, and atom is an O3*,
                      then have an end */
                   if (pair->last_base[i] == TRUE &&
                       !strncmp(atom->name+1,"O3",2))
                     {
                       /* If this is second strand, then move y upwards by
                          phosphate radius, otherwise move it down */
                       y2[nbond] = y2[nbond]
                         - side * P_RADIUS * size->txtscale;
                     }
                 }
             }
         }

       /* For lines to the base, draw straight across, and slightly below
          centre */
       if (line_type[nbond] == BASE)
         {
           y1[nbond] = y1[nbond] - drop_dist * size->txtscale;
           y2[nbond] = y1[nbond];
         }

       /* Bond type */
       bond_type[nbond] = HBOND;
       if (bond_point == NONPHOS)
         bond_type[nbond] = pair->type[i][ibond];
       else if (bond_point == PHOS)
         bond_type[nbond] = pair->ptype[i][ibond];
       else if (bond_point == LAST_PHOS && pair->next_pair != NULL)
/* v.1.1.2--> */
//         bond_type[nbond] = pair->next_pair->type[i][ibond];
         bond_type[nbond] = pair->next_pair->ptype[i][ibond];
/* <--v.1.1.2 */

       /* Increment bond number */
       nbond++;
     }


/*======================================================================
 * Draw
 *====================================================================*/
   pssave_();
   pslwid_(BDLNWID * size->txtscale);
   for (k = 0; k < nbond; k++)
   {
     if (bond_type[k] == HBOND)
       pscolb_(par->hbond_colour);
     else if (bond_type[k] == NBOND)
       pscolb_(par->nbond_colour);
     set_dash_type(bond_type[k]);
     pstrans_(x1[k], y1[k], cx, cy, &x_start, &y_start, par);
     pstrans_(x2[k], y2[k], cx, cy, &x_end, &y_end, par); 
     psline_(x_start, y_start, x_end, y_end);
   }  
   psdash_(0, 0);
   psrest_();

/* Return ============================================================*/
   return(0);
}





/***********************************************************************
 *
 * match_res_code - Function to find the matching one letter AA code 
 *
 **********************************************************************/
int match_res_code(PARAM *par, char *three_letter, char *res_name)
{
/* Declare variables =================================================*/
  int  i, istart, j;
  int found;
  char upper[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
  char lower[] = "abcdefghijklmnopqrstuvwxyz";

/* Copy full 3-letter code into output */
  strncpy(res_name,three_letter,3);
  res_name[3] = '\0';
   
/* Find match ========================================================*/
  /* If showing as single-letter amino acid codes, then get the
     corresponding code */
  if (par->one_letter_codes == TRUE)
    {
      for (i = 0; i < NUMRES; i++)
        if (!strncmp(three_letter, resname[i], 3))
          {
            res_name[0] = oneletter[i][0];
            res_name[1] = '\0';
          }
    }

  /* Otherwise, convert 2nd and 3rd letters to lower-case */
  else
    {
      /* If first letter is a space, then convert 3rd letter only */
      istart = 1;
      if (res_name[0] == ' ')
        istart = 2;

      for (i = istart; i < 3; i++)
        {
          found = FALSE;
          for (j = 0; j < 26 && found == FALSE; j++)
            {
              if (res_name[i] == upper[j])
                {
                  found = TRUE;
                  res_name[i] = lower[j];
                }
            }
        }
    }

/* Return to get_interact_det ========================================*/
   return(0);
}



/***********************************************************************
 *
 * text_draw_bdge_int - Function to draw bridging interactions to the 
 *                      Postscript file
 *
 **********************************************************************/
int 
text_draw_bdge_int(PAIRS *first_pair, PAIRS *last_pair, PAIRS *pair,
                   PAIRS *next_pair, SIZES *size, PARAM *par, int i,
                   int flag, int last_flag, double text_size[3])
{
/* Declare variables =================================================*/
   int    bond_point, j, k;
   int    on_left, water_flag;
   int    bond_type;
   double tsize, x1, y1;
   double x_init, y_init;
   double cx, cy;
   char   int_det[LINE_LEN], atom_name[5], atom_colour[COLNAME_LENGTH + 1];
   char   atom_name_combine[11], atom_colour_combine[COLNAME_LENGTH + 1];
   
/* Initialize variables ==============================================*/
   water_flag = FALSE;
   cx         = size->x_centre;
   cy         = size->y_centre;
   bond_type = HBOND;
   
/*======================================================================
  * Non-phosphate bonds 
  *===================================================================*/
/* Loop through the first level interactions =========================*/
   tsize = text_size[0];
   bond_point = NONPHOS;
   for (j = 0; j < pair->nbond[i]; j++)
   {
      
/* Loop through bridging interactions ================================*/
      for (k = 0; k < pair->nbridge[i][j]; k++)
      {
         get_interact_det(first_pair, last_pair, pair, par,
                          pair->bridgeto[i][j][k],
                          pair->bridge_combine[i][j][k], int_det, i,
                          bond_point, atom_name, atom_colour,
                          atom_name_combine, atom_colour_combine,
                          &water_flag, bond_type);
         
         /* Write to Postscript file */
         x_init = pair->bdgex_pos[i][j][k];
         y_init = pair->bdgey_pos[i][j][k];
         pstrans_(x_init, y_init, cx, cy, &x1, &y1, par);

          /* Determine whether doing left- or right-hand strand */
          if ((par->landscape == FALSE && i == ONE) ||
              (par->landscape == TRUE && i == TWO))
            on_left = TRUE;
          else
            on_left = FALSE;

          /* On left-hand strand */
          if (on_left == TRUE)
           {
             /* Print full atom and residue details */
             if (par->landscape == TRUE)
               psrotx_(x1, y1, (int) (tsize * size->txtscale), int_det);
             else
               psrtxt_(x1, y1, (int) (tsize * size->txtscale), int_det);

             /* Print just the atom name in the given colour */
             pssave_();
             pscolb_(atom_colour);
             if (par->landscape == TRUE)
               psrotx_(x1, y1, (int) (tsize * size->txtscale), atom_name);
             else
               psrtxt_(x1, y1, (int) (tsize * size->txtscale), atom_name);
             psrest_();
           }
         else
           {
             /* Print full atom and residue details */
             if (par->landscape == TRUE)
               psrltx_(x1, y1, (int) (tsize * size->txtscale), int_det);
             else
               psltxt_(x1, y1, (int) (tsize * size->txtscale), int_det);

             /* Print just the atom name in the given colour */
             pssave_();
             pscolb_(atom_colour);
             if (par->landscape == TRUE)
               psrltx_(x1, y1, (int) (tsize * size->txtscale), atom_name);
             else
               psltxt_(x1, y1, (int) (tsize * size->txtscale), atom_name);
             psrest_();
         }
      }
   }
   bridge_bond_lines(pair, par, size, i, NONPHOS, tsize);

/*======================================================================
 * Present pair phosphate bonds 
 *====================================================================*/
/* Loop through first_level interactions =============================*/
   tsize = text_size[1];
   bond_point = PHOS;
   if (!(flag == TRUE && i == ONE))
   {
      for (j = 0; j < pair->npbond[i]; j++)
      {
         
/* Loop through bridging interactions ================================*/
         for (k = 0; k < pair->npbridge[i][j]; k++)
         {
           /* Get the atom and residue names */
            get_interact_det(first_pair, last_pair, pair, par,
                             pair->pbridgeto[i][j][k],
                             pair->pbridge_combine[i][j][k], int_det, i,
                             bond_point, atom_name, atom_colour,
                             atom_name_combine, atom_colour_combine,
                             &water_flag,
                             bond_type);
            
            /* Write to postscript file */
            x_init = pair->pbdgex_pos[i][j][k];
            y_init = pair->pbdgey_pos[i][j][k];
            pstrans_(x_init, y_init, cx, cy, &x1, &y1, par);
            
            /* Determine whether doing left- or right-hand strand */
            if ((par->landscape == FALSE && i == ONE) ||
                (par->landscape == TRUE && i == TWO))
              on_left = TRUE;
            else
              on_left = FALSE;

            /* On left-hand strand */
            if (on_left == TRUE)
              {
                /* Print full atom and residue details */
                if (par->landscape == TRUE)
                  psrotx_(x1, y1, (int) (tsize * size->txtscale), int_det);
                else
                  psrtxt_(x1, y1, (int) (tsize * size->txtscale), int_det);

                /* Print just the atom name in the given colour */
                pssave_();
                pscolb_(atom_colour);
                if (par->landscape == TRUE)
                  psrotx_(x1, y1, (int) (tsize * size->txtscale), atom_name);
                else
                  psrtxt_(x1, y1, (int) (tsize * size->txtscale), atom_name);
                psrest_();
              }
            else
              {
                /* Print full atom and residue details */
                if (par->landscape == TRUE)
                  psrltx_(x1, y1, (int) (tsize * size->txtscale), int_det);
                else
                  psltxt_(x1, y1, (int) (tsize * size->txtscale), int_det);
                
                /* Print just the atom name in the given colour */
                pssave_();
                pscolb_(atom_colour);
                if (par->landscape == TRUE)
                  psrltx_(x1, y1, (int) (tsize * size->txtscale), atom_name);
                else
                  psltxt_(x1, y1, (int) (tsize * size->txtscale), atom_name);
                psrest_();
            }
         }
      }
      bridge_bond_lines(pair, par, size, i, PHOS, tsize);
   }
   
/*======================================================================
 * Next pair phosphate bonds when last pair on the page
 *====================================================================*/
   bond_point = LAST_PHOS;
/* v.1.1.1--> */
   tsize = text_size[LAST_PHOS];
/* <--v.1.1.1 */
   if (next_pair != NULL && last_flag == TRUE && i == ONE)
   {
      for (j = 0; j < next_pair->npbond[i]; j++)
      {
         for (k = 0; k < next_pair->npbridge[i][j]; k++)
         {
            get_interact_det(first_pair, last_pair, next_pair, par,
                             next_pair->pbridgeto[i][j][k],
                             next_pair->pbridge_combine[i][j][k], int_det,
                             i, bond_point, atom_name, atom_colour,
                             atom_name_combine, atom_colour_combine,
                             &water_flag, bond_type);

            /* Write to postscript file */
            x_init = next_pair->pbdgex_pos[i][j][k];
            y_init = next_pair->pbdgey_pos[i][j][k];
            pstrans_(x_init, y_init, cx, cy, &x1, &y1, par);
            
/* v.1.1.1--> */
//            /* Print full atom and residue details */
//            if (par->landscape == TRUE)
//              psrotx_(x1, y1, (int) (tsize * size->txtscale), atom_name);
//            else
//              psrtxt_(x1, y1, (int) (tsize * size->txtscale), atom_name);
//
//            /* Print just the atom name in the given colour */
//            pssave_();
//            pscolb_(atom_colour);
//            if (par->landscape == TRUE)
//              psrotx_(x1, y1, (int) (tsize * size->txtscale), atom_name);
//            else
//              psrtxt_(x1, y1, (int) (tsize * size->txtscale), atom_name);
//            psrest_();

            /* Determine whether doing left- or right-hand strand */
            if ((par->landscape == FALSE && i == ONE) ||
                (par->landscape == TRUE && i == TWO))
              on_left = TRUE;
            else
              on_left = FALSE;

            /* On left-hand strand */
            if (on_left == TRUE)
              {
                /* Print full atom and residue details */
                if (par->landscape == TRUE)
                  psrotx_(x1, y1, (int) (tsize * size->txtscale), int_det);
                else
                  psrtxt_(x1, y1, (int) (tsize * size->txtscale), int_det);

                /* Print just the atom name in the given colour */
                pssave_();
                pscolb_(atom_colour);
                if (par->landscape == TRUE)
                  psrotx_(x1, y1, (int) (tsize * size->txtscale), atom_name);
                else
                  psrtxt_(x1, y1, (int) (tsize * size->txtscale), atom_name);
                psrest_();
              }
            else
              {
                /* Print full atom and residue details */
                if (par->landscape == TRUE)
                  psrltx_(x1, y1, (int) (tsize * size->txtscale), int_det);
                else
                  psltxt_(x1, y1, (int) (tsize * size->txtscale), int_det);
                
                /* Print just the atom name in the given colour */
                pssave_();
                pscolb_(atom_colour);
                if (par->landscape == TRUE)
                  psrltx_(x1, y1, (int) (tsize * size->txtscale), atom_name);
                else
                  psltxt_(x1, y1, (int) (tsize * size->txtscale), atom_name);
                psrest_();
            }
/* <--v.1.1.1 */
         }
      }
      bridge_bond_lines(next_pair, par, size, i, PHOS, tsize);
   }
   

/* Return to text_bridge_interactions ================================*/
   return(0);
}



/***********************************************************************
 *
 * bridge_bond_lines - Function to draw bond lines to bridging 
 *                     interactions
 *
 **********************************************************************/
int bridge_bond_lines(PAIRS *pair, PARAM *par, SIZES *size, 
                      int i, int phos_flag, double text_size)
{
/* Declare variables =================================================*/
   int    j, k;
   int    is_sugar;
   int    line_type[MAXARRAY];
   double x1, x2, y1, y2;
   double xstart, xend, ystart, yend;
   double cx, cy;
   double extra;
   double bridge_length, resname_interact_width, text_width;
   ATOMS  *atom, *bridge_atom;

/* Initialize variables ==============================================*/
   cx       = size->x_centre;
   cy       = size->y_centre;
/*   bond_gap = BONDGAP * size->txtscale; */

/* Calculate minimum bridge length */
   text_width = 0.71 * text_size;
   resname_interact_width = 3 * text_width;
   bridge_length = resname_interact_width;
   
/*======================================================================
 * Portrait
 *====================================================================*/
/* Calculate positions ===============================================*/
/* Non-phosphate bonds ===============================================*/
   if (phos_flag == NONPHOS)
     {
/* Loop through first-level and bridging interactions ================*/
       for (j = 0; j < pair->nbond[i]; j++)
         {
           atom = pair->bondfrom[i][j];
           bridge_atom = pair->bondto[i][j];
            
           for (k = 0; k < pair->nbridge[i][j]; k++)
             {
               y1 = pair->bondy_pos[i][j];
               y2 = pair->bdgey_pos[i][j][k];
               is_sugar = is_sugar_atom(atom->name);
               if (is_sugar == TRUE)
                 line_type[0] = BACKBONE;
               else 
                 line_type[0] = BASE;

               /* Calculate extra length of bridging line to take into
                  account length of water's residue number */
               extra = 0.0;
               if (bridge_atom->resnum < 10)
                 extra = 2 * text_width;
               else if (bridge_atom->resnum < 100)
                  extra = 1 * text_width;
               else 
                 extra = 0.0;

               /* Have a fixed start-position for landscape mode */
               if (par->landscape == TRUE)
                 extra = 0.7 * text_width;
               
               if (i == ONE)
                 {
                   x2 = pair->bdgex_pos[i][j][k]
                     + WATER_RADIUS * size->txtscale;
                   x1 = pair->bdgex_pos[i][j][k]
                     + (bridge_length + extra) * size->txtscale;
                 }
               else if (i == TWO)
                 {
                   x2 = pair->bdgex_pos[i][j][k]
                     - WATER_RADIUS * size->txtscale;
                   x1 = pair->bdgex_pos[i][j][k]
                     - (bridge_length + extra) * size->txtscale;
                 }
               pssave_();
               pscolb_(par->hbond_colour);
               pslwid_(BDLNWID * size->txtscale);
               set_dash_type(HBOND);
               pstrans_(x1, y1, cx, cy, &xstart, &ystart, par);
               pstrans_(x2, y2, cx, cy, &xend, &yend, par); 
               psline_(xstart, ystart, xend, yend);
               psdash_(0, 0);
               psrest_();
             }
         }
     }         
/* Phosphate bonds ===================================================*/
   else if (phos_flag == PHOS)
     {
       for (j = 0; j < pair->npbond[i]; j++)
         {
           bridge_atom = pair->pbondto[i][j];

           for (k = 0; k < pair->npbridge[i][j]; k++)
             {
               y1 = pair->pbondy_pos[i][j];
               y2 = pair->pbdgey_pos[i][j][k];

               /* Calculate extra length of bridging line to take into
                  account length of water's residue number */
               extra = 0.0;
               if (bridge_atom->resnum < 10)
                 extra = 2 * text_width;
               else if (bridge_atom->resnum < 100)
                 extra = 1 * text_width;
               else 
                 extra = 0.0;

               /* Have a fixed start-position for landscape mode */
               if (par->landscape == TRUE)
                 extra = 0.7 * text_width;
               
               if (i == ONE)
                 {
                   x2 = pair->pbdgex_pos[i][j][k]
                     + WATER_RADIUS * size->txtscale;
                  x1 = pair->pbdgex_pos[i][j][k]
                    + (bridge_length + extra) * size->txtscale;
                 }
               else if (i == TWO)
                 {
                   x2 = pair->pbdgex_pos[i][j][k]
                     - WATER_RADIUS * size->txtscale;
                  x1 = pair->pbdgex_pos[i][j][k]
                    - (bridge_length + extra) * size->txtscale;
                 }
               
               pssave_();
               pscolb_(par->hbond_colour);
               pslwid_(BDLNWID * size->txtscale);
               set_dash_type(HBOND);
               pstrans_(x1, y1, cx, cy, &xstart, &ystart, par);
               pstrans_(x2, y2, cx, cy, &xend, &yend, par); 
               psline_(xstart, ystart, xend, yend);
               psdash_(0, 0);
               psrest_();
             }
         }
     }
   
/* Return to text_draw_bdge_int ======================================*/
   return(0);
}

/***********************************************************************
 *
 * set_dash_type - Set the dash-patter according to bond type
 *
 **********************************************************************/
void set_dash_type(int line_type)
{
/* Set the dash type =================================================*/
  if (line_type == HBOND)
    psdash_(1, 1);
  else if (line_type == NBOND)
    psdash_(2, 3);
}



/***********************************************************************
 *
 * psbasecol_ - Function to assign the NDB colouring to bases
 *
 **********************************************************************/
void psbasecol_(ATOMS *atom, PARAM *par, char *code)
{
/* Declare variables =================================================*/
  char one_letter_code;

/* Get the one-letter code for this base-name */
  is_nucleic_acid(atom->resname, &one_letter_code, par);

/* Assign colours ====================================================*/
  if (one_letter_code == 'A')
    pscolb_(par->adenine_colour);
  else if (one_letter_code == 'G')
    pscolb_(par->guanine_colour);
  else if (one_letter_code == 'T')
    pscolb_(par->thymine_colour);
  else if (one_letter_code == 'C')
    pscolb_(par->cytosine_colour);
  else if (one_letter_code == 'U')
    pscolb_(par->uracil_colour);
  else
    pscolb_(par->unknown_atom_colour);

  /* Return the one-letter code */
  if (one_letter_code == '\0')
    strcpy(code,atom->resname);
  else
    {
      code[0] = one_letter_code;
      code[1] = '\0';
    }
}



/***********************************************************************
 *
 * psbasenum_ - Function to write the base numbering 
 *
 **********************************************************************/
void psbasenum_(PAIRS *pair, SIZES *size, PARAM *par, int i, double x,
                double y)
{
/* Declare variables =================================================*/
   double text_size, x1, y1;
   double cx, cy;
   int base_number, ipos, first_char;
   char   resnum_char[5], temp_string[5];
   ATOMS *base1, *base2;
   
/* Initialize variables ==============================================*/
   x1    = x;
   y1    = y;
   cx    = size->x_centre;
   cy    = size->y_centre;
   base1 = pair->strand[ONE];
   base2 = pair->strand[TWO];

/* Get the character version of the number */
   if (i == ONE)
     base_number = base1->resnum;
   else
     base_number = base2->resnum;
   sprintf(temp_string,"%4d",base_number);
   temp_string[4] = '\0';

   /* Get rid of leading spaces */
   first_char = -1;
   for (ipos = 0; ipos < 4 && first_char == -1; ipos++)
     if (temp_string[ipos] != ' ')
       first_char = ipos;
   strcpy(resnum_char,temp_string+first_char);

/* Write number ======================================================*/
   /* Transform for landscape mode */
   pstrans_(x, y, cx, cy, &x1, &y1, par);

   /* Print the base number */
   text_size = BASENUM_TXT * size->txtscale;
   if (first_char == 1)
     text_size = text_size * 0.75;
   else if (first_char == 0)
     text_size = text_size * 0.5;
   psctxt_(x1, y1, (int) text_size, resnum_char);
}



/***********************************************************************
 *
 * psccol_  -  Function to define interior circle colour
 *
 **********************************************************************/

void 
psccol_(char *colour)
{
  fprintf(ps_file,"/Circol { %s } def\n",colour);
}


/***********************************************************************
 *
 * pscirc_  -  Function to draw a circle
 *
 **********************************************************************/

void 
pscirc_(double x, double y, double radius)
{
  fprintf(ps_file,"%7.2f %7.2f %7.2f Circle\n",x,y,radius);
}



/***********************************************************************
 *
 * pscolb_ - Function to write background colour level (for text, etc)
 *
 **********************************************************************/
void
pscolb_(char *text_string)
{
    fprintf(ps_file,"%s\n",text_string);
}



/***********************************************************************
 *
 * pscomm_ - Function to write out a comment to the PostScript file
 *
 **********************************************************************/
void pscomm_(char *text_string)
{
  fprintf(ps_file,"\n");
  fprintf(ps_file,"%% %s\n",text_string);
  fprintf(ps_file,"\n");
}



/***********************************************************************
 *
 * psltxt_ - Function to write out left aligned text to PostScript file
 *
 **********************************************************************/
void psltxt_(double x,double y,int size,
             char *text_string)
{
   fprintf(ps_file,"%7.2f %7.2f moveto\n", x, y);
   fprintf(ps_file,"(%s) %d Lalign\n",text_string,size);
   fprintf(ps_file,"(%s) %d Print\n",text_string,size);
}



/***********************************************************************
 *
 * psrtxt_ - Function to write out right aligned text to PostScript file
 *
 **********************************************************************/
void psrtxt_(double x,double y,int size,char *text_string)
{
   fprintf(ps_file,"%7.2f %7.2f moveto\n", x, y);
   fprintf(ps_file,"(%s) %d Ralign\n",text_string,size);
   fprintf(ps_file,"(%s) %d Print\n",text_string,size);
}



/***********************************************************************
 *
 * psrltx_ - Function to write out left-aligned text, rotated through
 *           90 degrees
 *
 **********************************************************************/
void
psrltx_(double x,double y,int size, char *text_string)
{
   fprintf(ps_file,"%7.2f %7.2f moveto\n", x, y);
   fprintf(ps_file,"(%s) %d LalignRot90 Rot90\n",text_string,size);
   fprintf(ps_file,"(%s) %d Print\n",text_string,size);
   fprintf(ps_file,"grestore\n");
}



/***********************************************************************
 *
 * psrotx_ - Function to write out right-aligned text, rotated through
 *           90 degrees
 *
 **********************************************************************/
void
psrotx_(double x,double y,int size, char *text_string)
{
   fprintf(ps_file,"%7.2f %7.2f moveto\n", x, y);
   fprintf(ps_file,"(%s) %d RalignRot90 Rot90\n",text_string,size);
   fprintf(ps_file,"(%s) %d Print\n",text_string,size);
   fprintf(ps_file,"grestore\n");
}



/***********************************************************************
 *
 * psctxt_ - Function to write out centred text to PostScript file
 *
 **********************************************************************/
void psctxt_(double x,double y,int size,char *text_string)
{
   fprintf(ps_file,"%7.2f %7.2f moveto\n", x, y);
   fprintf(ps_file,"(%s) %d Center\n",text_string,size);
   fprintf(ps_file,"(%s) %d Print\n",text_string,size);
}



/***********************************************************************
 *
 * psdash_  -  Switch dashed lines on/off
 *
 **********************************************************************/
void psdash_(int ONOFF1, int ONOFF2)
{
   static int last_ONOFF = -9;
   
   if (ONOFF1 != last_ONOFF && ONOFF2 != last_ONOFF)
   {
      if (ONOFF1 == 0)
         fprintf(ps_file,"[]  0 setdash\n");
      else
         fprintf(ps_file,"[ %2d %2d ] 0 setdash\n",ONOFF1,ONOFF2);
   }
   last_ONOFF = ONOFF1;
}



/***********************************************************************
 *
 * psendp_ - Function to write final lines for current page to 
 *           Postscript file
 *
 **********************************************************************/
void psendp_(void)
{
   fprintf(ps_file,"\n\n");
   fprintf(ps_file,"NucplotSave restore\n");
   fprintf(ps_file,"showpage\n");
}



/***********************************************************************
 *
 * psends_ - Function to draw the top and bottom sections of the DNA 
 *           chain
 *
 **********************************************************************/
void psends_(PAIRS *edge_pair[MEDIARRAY], PAIRS *pair, SIZES *size, 
             PARAM *par, int pairno, int which_end, int ds_flag[2],
             char chain_top[2], char chain_bottom[2])
{
/* Declare variables =================================================*/
  double chain_pos_x, chain_pos_y;
  double x, y, x1, y1, x2, y2;
  double cx, cy;
  double key_height, key_width;

/*Initialize variables ===============================================*/
   cx = size->x_centre;
   cy = size->y_centre;

/* Calculate position for placement of chain name */
   chain_pos_x = SUGAR_TO_BASE_WIDTH + P_TO_SUGAR_WIDTH;
   chain_pos_y = P_TO_P_HEIGHT + CHAIN_POS_HEIGHT;

/* Set text colour */
   pscolb_(par->text_colour);

/* Write top end of strands ==========================================*/
   if (which_end == TOP)
   {
      /* 5' end strand 1 */
      if (pair == edge_pair[pairno])
      {
         if (ds_flag[ONE] == TRUE)
         {
           pstrans_(pair->x_pos[ONE], 
                    pair->y_pos[ONE]
                    + P_TO_P_HEIGHT * size->txtscale * 2.0 / 3.0,
                    cx, cy, &x1, &y1, par);
           psctxt_(x1, y1, (int) (END_TXT * size->txtscale), "5'");
         }
      }

      /* Dashed lines */
      else if (pair != edge_pair[pairno])
      {
         pssave_();
         pslwid_(EDLNWID * size->txtscale);
         psdash_(5, 5);
         pscolb_(par->bbline_colour);
         
         /* Strand 1 */
         pstrans_(pair->x_pos[ONE],
                  pair->y_pos[ONE] + P_TO_P_HEIGHT * size->txtscale,
                  cx, cy, &x1, &y1, par);
         pstrans_(pair->x_pos[ONE],
                  pair->y_pos[ONE] + END_DASH_LENGTH * size->txtscale,
                  cx, cy, &x2, &y2, par);
         if (ds_flag[ONE] == TRUE)
            psline_(x1, y1, x2, y2);

         /* Strand 2 */
         pstrans_(pair->x_pos[TWO],
                  pair->y_pos[ONE] + P_TO_P_HEIGHT * size->txtscale,
                  cx, cy, &x1, &y1, par);
         pstrans_(pair->x_pos[TWO],
                  pair->y_pos[ONE] + END_DASH_LENGTH * size->txtscale,
                  cx, cy, &x2, &y2, par);
         if (ds_flag[TWO] == TRUE)
            psline_(x1, y1, x2, y2);

         psdash_(0, 0);
         psrest_();
      }

      /* Show chain id of strand 1 */
      plot_chain_id(chain_top[ONE],
                    pair->x_pos[ONE] - chain_pos_x * size->txtscale, 
                    pair->y_pos[ONE] + chain_pos_y * size->txtscale,
                    CHAINAME_TXT * size->txtscale, par, size);

      /* Show chain id of strand 2 */
      plot_chain_id(chain_top[TWO],
                    pair->x_pos[TWO] + chain_pos_x * size->txtscale, 
                    pair->y_pos[TWO] + chain_pos_y * size->txtscale,
                    CHAINAME_TXT * size->txtscale, par, size);

    }

/* Write bottom end of strands =======================================*/
   else if (which_end == BOTTOM)
     {
       /* 5' end strand 2 */
       if (pair == edge_pair[pairno])
         {
           if (ds_flag[TWO] == TRUE)
             {
               pstrans_(pair->x_pos[TWO],
                        pair->y_pos[TWO]
                        - P_TO_P_HEIGHT * size->txtscale * 2.0 / 3.0,
                        cx, cy, &x1, &y1, par);
               psctxt_(x1, y1, (int) (END_TXT * size->txtscale), "5'");
             }
         }

       /* Dashed lines */ 
       else if (pair != edge_pair[pairno])
         {
           pssave_();
           pslwid_(EDLNWID*size->txtscale);
           psdash_(5, 5);
           pscolb_(par->bbline_colour);

           /* Strand 1 */
           pstrans_(pair->x_pos[ONE],
                    pair->y_pos[ONE] - P_TO_P_HEIGHT * size->txtscale,
                    cx, cy, &x1, &y1, par);
           pstrans_(pair->x_pos[ONE],
                    pair->y_pos[ONE] - END_DASH_LENGTH * size->txtscale,
                    cx, cy, &x2, &y2, par);
           if (ds_flag[ONE] == TRUE)
             psline_(x1, y1, x2, y2);

           /* Strand 2 */
           pstrans_(pair->x_pos[TWO],
                    pair->y_pos[ONE] - P_TO_P_HEIGHT * size->txtscale,
                    cx, cy, &x1, &y1, par);
           pstrans_(pair->x_pos[TWO],
                    pair->y_pos[ONE] - END_DASH_LENGTH * size->txtscale,
                    cx, cy, &x2, &y2, par);
           if (ds_flag[TWO] == TRUE)
             psline_(x1, y1, x2, y2);
           
           psdash_(0, 0);
           psrest_();
         }

      /* Show chain id of strand 1 only if it is different from that
         printed at the top */
      if (chain_bottom[ONE] != chain_top[ONE])
        plot_chain_id(chain_bottom[ONE],
                      pair->x_pos[ONE] - chain_pos_x * size->txtscale, 
                      pair->y_pos[ONE] - chain_pos_y * size->txtscale, 
                      CHAINAME_TXT * size->txtscale, par, size);

      /* Show chain id of strand 2 only if it is different from that
         printed at the top */
       if (chain_bottom[TWO] != chain_top[TWO])
         plot_chain_id(chain_bottom[TWO],
                       pair->x_pos[TWO] + chain_pos_x * size->txtscale,
                       pair->y_pos[TWO] - chain_pos_y * size->txtscale,
                       CHAINAME_TXT * size->txtscale, par, size);

      /* Write out the plot title */
      if (par->plot_title == TRUE && par->title[0] != '\0')
        {
          /* Get the position of the title */
          if (par->landscape == TRUE)
            {
              x = size->y_centre;
              y = size->y_title_pos;
            }
          else
            {
              x = size->x_centre;
              y = size->y_title_pos;
            }
          pscomm_("Title");
          psctxt_(x, y, TITLE_TXT, par->title);
        }

       /* Write out the key, if required */
       if (par->plot_key == TRUE)
         {
          /* Get the position of the key-start */
           if (par->landscape == TRUE)
             {
               x = size->y_centre;
               y = size->y_keystart;
               key_height = KEY_HEIGHT_LANDSCAPE * KEY_SCALE;
               key_width = KEY_WIDTH_LANDSCAPE * KEY_SCALE;
             }
           else
             {
               x = size->x_centre;
               y = size->y_keystart;
               key_height = KEY_HEIGHT_PORTRAIT * KEY_SCALE;
               key_width = KEY_WIDTH_PORTRAIT * KEY_SCALE;
             }

           /* Plot the key */
           plot_key(par, size, x, y, key_height, key_width);
         }
     }
}



/***********************************************************************
 *
 * plot_key - Function to plot the key to the diagram
 *
 **********************************************************************/
void
plot_key(PARAM *par, SIZES *size, double x_start, double y_start,
         double key_height, double key_width)
{
  char nb_string[40], number_string[6];
  int ilimit, one_place, two_places;
  double radius, x[5], xl1, xl2, x_left, x1, x2, y[5], y1;

  /* Initialise variables */

  /* Plot start of key */
  pscomm_("Key");
  x_left = x_start - key_width / 2.0;
  if (par->landscape == TRUE)
    {
      x1 = x_left - 3.0 * KEY_TXT * KEY_SCALE;
      y1 = y_start + 0.5 * KEY_TXT * KEY_SCALE;
    }
  else
    {
      x1 = x_left;
      y1 = y_start + 1.5 * KEY_TXT * KEY_SCALE;
    }
  pscolb_(par->text_colour);
  psctxt_(x1, y1, (int) (KEY_TXT * KEY_SCALE), "Key");

  /* Plot sugar */
  x1 = x_left;
  y1 = y_start;
  get_pentagon(x, y, x1, y1, SUGAR_SIZE * KEY_SCALE, 0);
  plot_pentagon(x, y, par, PSLNWID * KEY_SCALE);

  /* Plot a number inside the sugar */
  psctxt_(x1, y1, (int) (BASENUM_TXT * KEY_SCALE), "3");

  /* Print description */
  x1 = x1 + 1.8 * SUGAR_SIZE * KEY_SCALE;
  psltxt_(x1, y1, (int) (KEY_TEXT_TXT * KEY_SCALE), "Backbone sugar and base-number");

  /* Get position and radius of phosphate */
  x1 = x_left;
  y1 = y_start - 2.0 * SUGAR_SIZE * KEY_SCALE;
  radius = P_RADIUS * KEY_SCALE;

  /* Draw the phosphate */
  pssave_();
  pscolb_(par->phosphate_colour);
  psccol_(par->inside_phosphate_colour);
  pscirc_(x1, y1, radius);
  psctxt_(x1, y1, (int) (PHOS_TXT * KEY_SCALE), "P");
  psrest_();

  /* Print explanation */
  x1 = x1 + 1.8 * SUGAR_SIZE * KEY_SCALE;
  psltxt_(x1, y1, (int) (KEY_TEXT_TXT * KEY_SCALE), "Phosphate group");

  /* Print asterisk and explanation */
  if (par->landscape == TRUE)
    {
      x1 = x_start + 300.0 * KEY_SCALE;
      y1 = y_start - 1.3 * KEY_TEXT_TXT * KEY_SCALE;
    }
  else
    {
      x1 = x_left + 1.5 * SUGAR_SIZE * KEY_SCALE;
      y1 = y1 - 1.5 * KEY_TEXT_TXT * KEY_SCALE;
    }
  pssave_();
  pscolb_(par->res_name_colour);
  psrtxt_(x1, y1, (int) (RESNAME_TXT * KEY_SCALE), "*");
  psrest_();
  x1 = x1 + 0.3 * SUGAR_SIZE * KEY_SCALE;
  psltxt_(x1, y1, (int) (KEY_TEXT_TXT * KEY_SCALE),
          "Residue/water on plot more than once");

  /* Plot items on the right-hand side of the key */
  if (par->landscape == TRUE)
    x2 = x_start - 50.0 * KEY_SCALE;
  else
    x2 = x_start + 100.0 * KEY_SCALE;
  y1 = y_start;

  /* Plot H-bond */
  radius = WATER_RADIUS * KEY_SCALE;
  x1 = x2;
  xl2 = x1;
  xl1 = xl2 - 60.0 * KEY_SCALE;
  pssave_();
  pscolb_(par->hbond_colour);
  pslwid_(BDLNWID * KEY_SCALE);
  set_dash_type(HBOND);
  psline_(xl1, y1, xl2, y1);
  psrest_();

  /* Print explanation */
  x1 = x2 + radius;
  psltxt_(x1, y1, (int) (KEY_TEXT_TXT * KEY_SCALE), "Hydrogen bond to DNA");

  /* Plot line representing non-bonded contact */
  y1 = y1 - 1.3 * KEY_TEXT_TXT * KEY_SCALE;
  pssave_();
  pscolb_(par->nbond_colour);
  pslwid_(BDLNWID * KEY_SCALE);
  set_dash_type(NBOND);
  psline_(xl1, y1, xl2, y1);
/* v.1.1.6--> */
  psdash_(0, 0);
/* <--v.1.1.6 */
  psrest_();

  /* Print explanation */
  ilimit = (int) (par->nb_distance * 10.0);
  one_place = (int) (ilimit * 10 + 0.5);
  two_places = (int) (par->nb_distance * 100.0 + 0.5);
  if (one_place == two_places)
    sprintf(number_string,"%3.1f",par->nb_distance);
  else
    sprintf(number_string,"%4.2f",par->nb_distance);
  sprintf(nb_string,"Nonbonded contact to DNA (< %sA)",
          number_string);
  psltxt_(x1, y1, (int) (KEY_TEXT_TXT * KEY_SCALE), nb_string);

  /* Shift water by its radius */
  radius = WATER_RADIUS * KEY_SCALE;
  if (par->landscape == TRUE)
    {
      x2 = x_start + 300.0 * KEY_SCALE;
      x1 = x2 - radius;
      y1 = y_start;
    }
  else
    {
      x1 = x2 - radius;
      y1 = y1 - 1.3 * KEY_TEXT_TXT * KEY_SCALE;
    }

  /* Draw water */
  pssave_();
  pscolb_(par->water_colour);
  psccol_(par->inside_water_colour);
  pscirc_(x1, y1, radius);
  psctxt_(x1, y1, (int) (WATER_TXT * KEY_SCALE), "W");

  /* Print a water number by its side */
  x1 = x1 - 1.6 * radius;
  pscolb_(par->res_name_colour);
  psrtxt_(x1, y1, (int) (RESNAME_TXT * KEY_SCALE), "88");
  psrest_();

  /* Print explanation */
  x1 = x2 + radius;
  psltxt_(x1, y1, (int) (KEY_TEXT_TXT * KEY_SCALE),
          "Water molecule and number");
}



/***********************************************************************
 *
 * plot_chain_id - Print the chain id of the given strand at the top or
 *                 bottom of the plot
 *
 **********************************************************************/
void 
plot_chain_id(char chain, double x, double y, double text_size,
              PARAM *par, SIZES *size)

{
/* Declare variables =================================================*/
  char   chain_string[8];
  double x1, y1;
  double cx, cy;

/*Initialize variables ===============================================*/
  cx = size->x_centre;
  cy = size->y_centre;

  /* Plot is chain not null */
  if (chain != '\0')
    {
      /* Show chain id of given strand */
      strcpy(chain_string,"Chain ");
      if (chain == ' ')
        chain_string[6] = '-';
      else
        chain_string[6] = chain;
      chain_string[7] = '\0';

      /* Transform for landscape mode */
      pstrans_(x, y, cx, cy, &x1, &y1, par);

      /* Plot the chain-id */
      psctxt_(x1, y1, (int) text_size, chain_string);
    }
}



/***********************************************************************
 *
 * psline_ - Function to write line out to PostScript file
 *
 **********************************************************************/
void psline_(double x1,double y1,double x2,double y2)
{
  fprintf(ps_file,"%7.2f %7.2f %7.2f %7.2f L\n",x1, y1, x2, y2);
}



/***********************************************************************
 *
 * pslwid_ - Function to write line-width out to PostScript file
 *
 **********************************************************************/
void pslwid_(double width)
{
    fprintf(ps_file,"%7.2f W\n",width);
}



/***********************************************************************
 *
 * pspage_ - Function to write starting lines for new page to 
 *           Postscript file
 *
 **********************************************************************/
void pspage_(int page_no, SIZES *size, PARAM *par)
{
/* Declare variables =================================================*/
   double xx1, xx2, xy1, xy2;
   
/* initialize variables ==============================================*/


/* Define border around picture ======================================*/
  xx1 = -1.0;
  xx2 = 650.0;
  xy1 = -1.0;
  xy2 = 951.0;
   
/*======================================================================
 * Write out headings to PostScript file
 *====================================================================*/
  fprintf(ps_file,"\n\n%%%%Page: %d %d\n", page_no, page_no);
  fprintf(ps_file,"/NucplotSave save def\n");
  fprintf(ps_file,"%6.2f%6.2f moveto %7.2f%7.2f lineto%7.2f%7.2f lineto\n",
          xx1,xy1,xx2,xy1,xx2,xy2);
  fprintf(ps_file,"%7.2f%7.2f lineto closepath\n",xx1,xy2);
  fprintf(ps_file,"gsave 1.0000 setgray fill grestore\n");
  fprintf(ps_file,"stroke gsave\n");

/* Draw in the background colour =====================================*/
   pscomm_("Background");
   pshade_("Background_col");
   psubox_(BBOXX1, BBOXY1, BBOXX2, BBOXY1, BBOXX2, BBOXY2, BBOXX1, BBOXY2);

/* Rotate page =======================================================*/
   pscomm_("Rotate");
   if (par->landscape == TRUE)
      psrot_(BBOXX2 + BBOXX1, 0);

}



/***********************************************************************
 *
 * psphostxt_ - Function to write out text type phosphate to 
 *              Postscript file
 *
 **********************************************************************/
void
psphostxt_(PAIRS *pair, SIZES *size, PARAM *par, double xstart,
           double x_line[2], double y_line[2], 
           double prev_xline[2], double prev_yline[2], 
           int flag, int end_flag, int first_base, int last_base,
           int top_of_page)
{
/* Declare variables =================================================*/
  int    bold_flag;
  int colour;
  int side;
  double x1, x2, x3, y1, y2, y3;
  double y_diff;
  double cx, cy;
  double radius;

/* Initialize variables ==============================================*/
  y_diff    = y_line[0] - y_line[1];
  cx        = size->x_centre;
  cy        = size->y_centre;
  radius    = P_RADIUS * size->txtscale;
   
  /* Get the appropriate side (left = -ve, right = +ve) */
  side = flag * 2 - 1;

  if (par->highlight == TRUE)
    bold_flag = FALSE;
  else if (par->highlight == FALSE)
    bold_flag = TRUE;

  /* Get the backbone colour */
  colour = pair->colour[flag];
  
/* Calculate coordinates =============================================*/
  if (end_flag == FALSE)
    {
      /* Left coordinates */
      if (flag == LEFT)
        {         
          x1 = x_line[0];
          x2 = xstart;
          x3 = x_line[1];
          y1 = y_line[0];
          y2 = pair->y_pos[ONE] + (P_TO_P_HEIGHT / 2.0) * size->txtscale;
          y3 = y1 + P_TO_P_HEIGHT * size->txtscale - y_diff;
          if (pair->npbond[ONE] != 0)
            bold_flag = TRUE;
        }
      
      /* Right coordinates */
      else if (flag == RIGHT)
        {
          x1 = x_line[0];
          x2 = xstart;   
          x3 = x_line[0];
          y1 = y_line[1];
          y2 = pair->y_pos[TWO] - (P_TO_P_HEIGHT / 2.0) * size->txtscale;
          y3 = y_line[1] - P_TO_P_HEIGHT * size->txtscale + y_diff;
          if (pair->npbond[TWO] != 0)
            bold_flag = TRUE;
        }
    }
  else if (end_flag == TRUE)
    {
      /* Left coordinates */
      if (flag == LEFT)
        {         
          x1 = x_line[0];
          x2 = xstart;
          x3 = x_line[1];
          y1 = y_line[1];
          y2 = pair->y_pos[ONE] - (P_TO_P_HEIGHT / 2.0) * size->txtscale;
          y3 = y_line[1] - P_TO_P_HEIGHT * size->txtscale + y_diff;
          if (pair->npbond[ONE] != 0)
            bold_flag = TRUE;
        }
    }
   
  /* Translate and rotate for Landscape */
  pstrans_(x1, y1, cx, cy, &x1, &y1, par);
  pstrans_(x2, y2, cx, cy, &x2, &y2, par);
  pstrans_(x3, y3, cx, cy, &x3, &y3, par);

/* Draw phosphate ====================================================*/
  pscomm_("phosphate");
  pssave_();
  pslwid_(PSLNWID * size->txtscale);
  if (bold_flag == FALSE)
    psetgray_(0.2);
  else if (bold_flag == TRUE)
    psetgray_(0.0);

  /* If this is the first base, then plot oxygen only */
  if (first_base == TRUE)
    {
      x3 = x2;
      y3 = y2 + side * P_RADIUS * size->txtscale;
      pscolb_(par->backbone_colour[colour]);
      psline_(x1, y1, x3, y3);
      pscolb_(par->oxygen_colour);
      psccol_(par->oxygen_colour);
      pscirc_(x3, y3, radius / 4.0);
    }
  else if (flag == LEFT && top_of_page == TRUE)
    {
      x3 = x2;
      y3 = y2 + side * P_RADIUS * size->txtscale;
      pscolb_(par->backbone_colour[colour]);
      psline_(x1, y1, x3, y3);
    }
  else
    {
      pscolb_(par->backbone_colour[colour]);
      psline_(x1, y1, x2, y2);
      psline_(x2, y2, x3, y3);
      pscolb_(par->phosphate_colour);
      psccol_(par->inside_phosphate_colour);
      pscirc_(x2, y2, radius);
      psctxt_(x2, y2, (int) (PHOS_TXT * size->txtscale), "P");
    }

  /* If this is a last-base of a chain, then plot oxygen on the
     other side */
  if (last_base == TRUE)
    {
      /* Left coordinates */
      if (flag == LEFT)
        {         
          x1 = x_line[0];
          x2 = xstart;
          x3 = x_line[1];
          y1 = y_line[1];
          y2 = pair->y_pos[ONE] - (P_TO_P_HEIGHT / 2.0) * size->txtscale;
          y3 = y1 - P_TO_P_HEIGHT * size->txtscale - y_diff;
        }
      
      /* Right coordinates */
      else if (flag == RIGHT)
        {
          x1 = x_line[0];
          x2 = xstart;   
          x3 = x_line[0];
          y1 = y_line[0];
          y2 = pair->y_pos[TWO] + (P_TO_P_HEIGHT / 2.0) * size->txtscale;
          y3 = y1 + P_TO_P_HEIGHT * size->txtscale + y_diff;
        }

      /* Plot the oxygen */
      x3 = x2;
      y3 = y2 - side * P_RADIUS * size->txtscale;
      pscolb_(par->backbone_colour[colour]);
      psline_(x1, y1, x3, y3);
      pscolb_(par->oxygen_colour);
      psccol_(par->oxygen_colour);
      pscirc_(x3, y3, radius / 4.0);
    }

  /* If this is the top of the page, and on the right, draw the
     extension */
  else if (flag == RIGHT && top_of_page == TRUE)
    {
      x1 = x_line[0];
      x2 = xstart;   
      x3 = x_line[0];
      y1 = y_line[0];
      y2 = pair->y_pos[TWO] + (P_TO_P_HEIGHT / 2.0) * size->txtscale;
      y3 = y1 + P_TO_P_HEIGHT * size->txtscale + y_diff;

      /* Plot the oxygen */
      x3 = x2;
      y3 = y2 - side * P_RADIUS * size->txtscale;
      pscolb_(par->backbone_colour[colour]);
      psline_(x1, y1, x3, y3);
    }

  psrest_();
}



/***********************************************************************
 *
 * psrest_ - Function to restore graphics state
 *
 **********************************************************************/
void psrest_(void)
{
    fprintf(ps_file,"grestore\n");
}



/***********************************************************************
 *
 * psrot_  -  Function to rotate graphics page through 90 degrees
 *
 **********************************************************************/
void psrot_(float new_origin_x, float new_origin_y)
{
   fprintf(ps_file," %7.2f %7.2f moveto Rot90\n",new_origin_x,
           new_origin_y);
}



/***********************************************************************
 *
 * pssave_ - Function to save graphics state
 *
 **********************************************************************/
void pssave_(void)
{
    fprintf(ps_file,"gsave\n");
}


/* v.1.1.7--> */
/***********************************************************************

pselip_  -  Write out a coloured ellipse to PostScript file

***********************************************************************/

void pselip_(float x,float y,float width,float ratio,float angle)
{
  fprintf(ps_file,"%7.2f %7.2f moveto\n",x,y);
  fprintf(ps_file,"%7.2f RotAngle\n",-angle);
  fprintf(ps_file,"0.0 0.0 %7.2f 1.0 %7.2f 0 0 Oellipse\n",width,ratio);
  fprintf(ps_file,"grestore\n");
}
/* v.1.1.7--> */


/***********************************************************************
 *
 * psetgray_ - Function to set the gray scale 
 *
 **********************************************************************/
void psetgray_(double grayscale)
{
   fprintf(ps_file,"%4.2f setgray\n", grayscale);
}



/***********************************************************************
 *
 * pshade_ - Function to write colour/shade out to PostScript file
 *
 **********************************************************************/
void pshade_(char *colour)
{
    fprintf(ps_file,"G %s\n",colour);
}



/***********************************************************************
 *
 * psugar_ - Function to draw the sugar
 *
 **********************************************************************/
void psugar_(PAIRS *pair, SIZES *size, PARAM *par, int i, 
             double xstart, double ystart, double x_line[2],
             double y_line[2])
{
/* Declare variables =================================================*/
   int    colour, j;
   double x[5], y[5], x1, y1;
   float  x_diff, x_end, y_end;
   double cx, cy;
   double sugar_size;

/* Initialize variables ==============================================*/
   sugar_size = SUGAR_SIZE * size->txtscale;
   cx        = size->x_centre;
   cy        = size->y_centre;
   colour = pair->colour[i];

/* Calculate coordinates of the pentagon =============================*/
   get_pentagon(x, y, xstart, ystart, sugar_size, i);

/* Record positions to which the phosphate lines must go =============*/
   x_line[0] = x[0];
   x_line[1] = x[1];
   y_line[0] = y[0];
   y_line[1] = y[1];

   /* Calculate length of bond from sugar to base */
   x_diff = (float) (fabs(pair->x_numpos[i] - x[3])
     - size->txtscale * BASENAME_TXT / 2.0);
   if (pair->x_numpos[i] > x[3])
     x_end = (float) (x[3] + x_diff);
   else
     x_end = (float) (x[3] - x_diff);

   /* Transform for landscape mode */
   pstrans_(x_end, y[3], cx, cy, &x1, &y1, par);
   x_end = (float) x1;
   y_end = (float) y1;
   for (j = 0; j < 5; j++)
     {
       pstrans_(x[j], y[j], cx, cy, &x1, &y1, par);
       x[j] = x1;
       y[j] = y1;
     }

  /* Draw bond from base to sugar */
  pslwid_(PSLNWID * size->txtscale);
  pscolb_(par->backbone_colour[colour]);
  if (par->landscape == TRUE)
    psline_(x[3], y[3], x[3], y_end);
  else
    psline_(x[3], y[3], x_end, y[3]);

/* Draw pentagon =====================================================*/   
   plot_pentagon(x, y, par, PSLNWID * size->txtscale);
}




/***********************************************************************
 *
 * get_pentagon - Get the coordinates of the five corners of the petagon
 *                representing the sugar, given the coords of its centre
 *                and its relative size
 *
 **********************************************************************/
void
get_pentagon(double x[5], double y[5], double xstart, double ystart,
             double sugar_size, int i)
{
/* Declare variables =================================================*/
   int    j;
   double xrel[5], yrel[5];
   double theta;

/* Initialize variables ==============================================*/
   theta  = 72.5 * (PI / 180);

/* Calculate coordinates of the pentagon =============================*/
   xrel[0] = sugar_size * cos(theta/2); 
   yrel[0] = sugar_size * sin(theta/2);

   xrel[1] = xrel[0]*cos(theta) + yrel[0]*sin(theta);
   yrel[1] = -xrel[0]*sin(theta) + yrel[0]*cos(theta);

   xrel[2] = xrel[1]*cos(theta) + yrel[1]*sin(theta);
   yrel[2] = -xrel[1]*sin(theta) + yrel[1]*cos(theta);

   xrel[3] = xrel[2]*cos(theta) + yrel[2]*sin(theta);
   yrel[3] = -xrel[2]*sin(theta) + yrel[2]*cos(theta);

   xrel[4] = xrel[3]*cos(theta) + yrel[3]*sin(theta);
   yrel[4] = -xrel[3]*sin(theta) + yrel[3]*cos(theta);

   /* If sugar is on strand 1, then flip to face the other way */
   if (i == ONE)
      for (j = 0; j < 5; j++)
      {
         xrel[j] = -xrel[j];
         yrel[j] = yrel[j];
      }
   
   /* Shift pentagon to required coordinates */
   for (j = 0; j < 5; j++)
   {
      x[j] = xrel[j] + xstart;
      y[j] = yrel[j] + ystart;
   }
 }


   
/***********************************************************************
 *
 * plot_pentagon - Plot the given pentagon
 *
 **********************************************************************/
void
plot_pentagon(double x[5], double y[5], PARAM *par, double line_width)
{
/* Declare variables =================================================*/


  pscomm_("sugar");
  pssave_();
  pslwid_(line_width);

  /* Draw a pentagon representing the sugar */
  pscolb_(par->sugar_colour);
  fprintf(ps_file, "/Sugarcol { %s } def\n", par->inside_sugar_colour);
  fprintf(ps_file, "%7.2f %7.2f moveto %7.2f %7.2f lineto\n", 
          x[0], y[0], x[1], y[1]);
  fprintf(ps_file,
          "%7.2f %7.2f lineto %7.2f %7.2f lineto %7.2f %7.2f lineto\n", 
          x[2], y[2], x[3], y[3], x[4], y[4]);
  fprintf(ps_file, "gsave Sugarcol fill grestore closepath stroke\n");
  psrest_();
 }


   
/***********************************************************************
 *
 * pstrans_ - Function to apply a transformation of a 90 degree 
 *            anticlockwise rotation followed by a flip about the 
 *            x-axis IF NEEDED
 *
 **********************************************************************/
void pstrans_(double x, double y, double xo, double yo, double *x_new, double *y_new, PARAM *par)
{
/* Declare variables =================================================*/
   double xr, yr;
   double xr_new, yr_new;
   double xo_new, yo_new;
   
/* Check whether transformation required or not ======================*/
   if (par->landscape == FALSE)
   {
      *x_new = x;
      *y_new = y;
   }
   
   else if (par->landscape == TRUE)
   {
/* Transformation ====================================================*/
      /* Set new centre */
      xo_new = yo;
      yo_new = xo;
      
      /* Set coordinates relative to centre */
      xr = x - xo;
      yr = y - yo;
      
      /* Transform */
      xr_new = -yr;
      yr_new = -xr;
      
      /* Return to real coordinates */
      *x_new = xr_new + xo_new;
      *y_new = yr_new + yo_new;
   }
   
/* Return ============================================================*/
}



/***********************************************************************
 *
 * psubox_ - Function to write out unbounded box to PostScript file
 *
 **********************************************************************/
void
psubox_(double x1,double y1,double x2,double y2,double x3,double y3,
        double x4,double y4)
{
  fprintf(ps_file,"%7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f Pl4\n",
          x1, y1, x2, y2, x3, y3, x4, y4);
}



/***********************************************************************
 *
 * pswatertxt_ - Function to draw the water to the Postscript file 
 *
 **********************************************************************/
void 
pswatertxt_(SIZES *size, PARAM *par, double x, double y,
            char int_det[LINE_LEN], int i, double text_size,
            int ibond)
{
/* Declare variables =================================================*/
  int side;
  double radius, fraction, w_size;
  double cx, cy;
  double xnum_pos, ynum_pos;
  double x1, y1;

  /* Get the appropriate side (left = -ve, right = +ve) */
  side = i * 2 - 1;

/* Initialize variables ==============================================*/
  cx        = size->x_centre;
  cy        = size->y_centre;
  fraction = text_size / RESNAME_TXT;
  radius = fraction * WATER_RADIUS * size->txtscale;
  w_size = fraction * WATER_TXT * size->txtscale;

  /* Shift water by its radius */
  if (i == ONE)
    x = x - radius;
  else
    x = x + radius;

  /* Transform for landscape mode */
  pstrans_(x, y, cx, cy, &x1, &y1, par);

  /* Draw water and circle with colour option */
  pssave_();
  pslwid_(1.0 * size->txtscale);
  pscolb_(par->water_colour);
  psccol_(par->inside_water_colour);
  pscirc_(x1, y1, radius);
  psctxt_(x1, y1, (int) w_size, "W");

  /* If landscape mode, reduce text size */
  if (par->landscape == TRUE)
    text_size = 0.7 * text_size;

  /* Calculate position of water-number */
  xnum_pos = x + side * 1.6 * radius;
  if (par->landscape == TRUE && ibond % 2 != 0)
    xnum_pos = xnum_pos + side * 0.6 * text_size * size->txtscale;
  ynum_pos = y;
  
  /* Transform for landscape mode */
  pstrans_(xnum_pos, ynum_pos, cx, cy, &x1, &y1, par);
  pscolb_(par->res_name_colour);

   /* Write number */
   if (i == ONE)
     {
       if (par->landscape == TRUE)
         psctxt_(x1, y1, (int) (text_size * size->txtscale), int_det);
       else
         psrtxt_(x1, y1, (int) (text_size * size->txtscale), int_det);
     }
   else if (i == TWO)
     {
       if (par->landscape == TRUE)
         psctxt_(x1, y1, (int) (text_size * size->txtscale), int_det);
       else
         psltxt_(x1, y1, (int) (text_size * size->txtscale), int_det);
     } 

   psrest_();
}



/***********************************************************************
 * 
 * psclos_ - Function to write final lines to PostScript file and close
 *
 **********************************************************************/
void
psclos_(void)
{
/* Write out closing lines to PostScript file ========================*/
    fprintf(ps_file,"%%%%Trailer\n");
    fprintf(ps_file,"%%%%BoundingBox: %d %d %d %d\n",(int) BBOXX1 - 1,
        (int) BBOXY1 - 1, (int) BBOXX2 + 1, (int) BBOXY2 + 1);
    fprintf(ps_file,"%%%%EOF\n");

/* Close postscript file =============================================*/
    fclose(ps_file);
}



/**************************************************************************
 *
 * open_file - Function to open files with appropriate error messages
 *
 *************************************************************************/
FILE*
open_file(char filename[FILENAME_LEN], char *open_type)
{
   FILE *file;
   
   if ((file = fopen(filename, open_type)) == NULL)
   {
      printf("\n*** ERROR - Unable to open file %s for reading\n",
             filename);
      exit(1);
   }
   return(file);
}

/***********************************************************************

store_string  -  Allocate memory for given string and store pointer to
                 it

***********************************************************************/

void store_string(char **string_ptr,char *string)
{
  int len;

  /* Get the length of the string */
  len = strlen(string);

  /* Allocate memory for the string */
  *string_ptr = (char *) malloc(sizeof(char)*(len + 1));

  /* Store the string */
  strcpy(*string_ptr,string);
}

