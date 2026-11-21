/***********************************************************************

  dimhtml.h - Include file for dimhtml.c

***********************************************************************/

/* Array parameters */
#define COLNAMLEN           12
#define FILENAME_LEN       512
#define LINELEN           1024
#define MAX_CHAIN_COL        8
#define MAX_TOKEN            8
#define NCONSERV            10
#define STRING_LEN        1000
#define TOKEN_LEN           50

/* Random number generator parameters */
/* debug */
#define RANDOM_START       FALSE
/* debug */
#define LEFT_SEQ              0
#define RIGHT_SEQ             1
#define BOTH_SEQS             2
#define LEFT_SEQ_CUTOFF    (0.4)
#define RIGHT_SEQ_CUTOFF   (0.8)

/* Interaction types */
#define HBONDS                0
#define CONTACTS              1
#define SALT_BRIDGES          2
#define DISULPHIDES           3

#define NTYPES                4

/* Non-bonded contact types */
#define HYDROPHOBIC_ONLY      0
#define HYDROPHOBIC_ANY       1
#define ANY_ANY               2

/* Maximum number of residues that can be plotted */
/* v.5.4.4 --> */
//#define MAX_RESIDUES        300
#define MAX_RESIDUES       1000
/* <-- v.5.4.4 */

/* Colouring schemes */
#define RESIDUE_TYPE          0
#define CONSERVATION          1

/* Aligment parameters */
#define GAP_PENALTY           5
#define GAP_SCORE            20
#define HBOND_WEIGHT         12
#define DISULPHIDE_WEIGHT    10
#define SALT_BRIDGE_WEIGHT    5
#define WATER_BRIDGE_WEIGHT   5
#define CROSS_PENALTY         5
/* v.5.3.8--> */
//#define MAX_ALIGN_LEN        40
#define MAX_ALIGN_LEN        20
/* <--v.5.3.8 */

/* Alignment positions between added gaps, and margin either side of
   alignment to improve minimization */
#define ALIGN_FRAG_LEN        5
#define ALIGN_GAP             3
#define ALIGN_MARGIN          5

/* Minimization parameters */
#define MAX_STEPS           100
/* v.5.4.4 --> */
#define MIN_STEPS            50
/* <-- v.5.4.4 */
#define STARTING_TEMP       (100.0)
#define BOLTZMANN_CONSTANT  ( 1.38062E-23 )

/* Data types read in from dimhtmlXY.dat file */
#define NOTHING              -1
#define RESIDUE_DATA          0
#define ALIGNMENT_DATA        1
#define INTERACTION_DATA      2

/* Residue conservation ranges */
static double RANGE[NCONSERV + 1]
    = {  0.00,  0.10,  0.20,  0.30, 
           0.40,  0.50,  0.60,  0.70, 
           0.80,  0.90,  1.000001 };

/* Upper domain colour */
#define UPPER_DOMAIN_COL      7

/* v.5.4.4 --> */
//FILE *dump_file_ptr;
/* <-- v.5.4.4 */

/* Special files */
#define STANDARD_PDB          0
#define UPLOADED_PDB          1
#define MCSG_PDB              2

/* Order residues in first, second or neither sequence */
#define NEITHER               0
#define FIRST_SEQ             1
#define SECOND_SEQ            2

/* C of M separations for dimplot.rcm file (no longer used) */
#define RCM_SEP_X           8.0
#define RCM_WATER_SEP_X     2.0
#define RCM_SEP_Y           4.0

/* v.5.3.4--> */
/* CofM parameters for dimplot.rcm file */
#define GAP_Y             ( 1.0 )
#define IDEAL_GAP         ( 3.0 )
#define SIZE_HYDROPHOBIC  ( 3.0 )
#define SIZE_HBOND        ( 5.0 )
/* <--v.5.3.4 */

/* v.5.3.8--> */
/* Ellipse paramters for highlighted residue */
#define ELLIPSE_FACTOR    ( 1.2 )
#define ELLIPSE_LINEWIDTH ( 4.0 )
#define NEAR_ZERO         ( 0.000001 )

/* The DisaStr variant types */
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

/* Structure for input parameters
   ------------------------------ */
struct parameters
{
  char                    pdb_code[5];
  char                    chain[2];
  int                     domain[2];
  int                     chains;
  int                     pdb_type;
  int                     colour_scheme;
  int                     dimplot;
  int                     use_data;
  int                     save_data;
  char                    work_dir[FILENAME_LEN];
  int                     labels_outside;
  int                     show_chains;
  int                     order;
  int                     flip_seqs;
  int                     contact_type;
  int                     include_waters;
  int                     have_waters;
/* v.5.3.6--> */
  float                   size_hydrophobic;
  int                     portrait;
/* <--v.5.3.6 */
/* v.5.3.8--> */
  char                    highlight_res[16];
  char                    *disastr_dir;
/* <--v.5.3.8 */
  int                     run_64;
};

/* Structure for each residue
   -------------------------- */
struct residue
{
  char               res_name[4];
  char               res_num[6];
  char               chain;
  int                chain_no;
  double             ycoord;
  int                number;
  double             conservation;
  char               colour[COLNAMLEN + 1];
  int                ref_pos;
  int                nhbonds;
/* v.5.3.8--> */
  int                deleted;
/* <--v.5.3.8 */
  struct link        *first_link_ptr;
  struct link        *last_link_ptr;
  struct residue     *align_residue_ptr;
  struct residue     *next_residue_ptr;
};

/* Structure for each residue link
   ------------------------------- */
struct link
{
  int                count[3];
  int                salt_bridge;
  int                disulphide;
  int                water_bridge;
  struct residue     *residue_ptr[2];
  struct link        *nextres_link_ptr[2];
  struct link        *next_link_ptr;
};

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
