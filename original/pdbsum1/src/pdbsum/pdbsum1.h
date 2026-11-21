/***********************************************************************

  pdbsum1.h - Include file for pdbsum1.c

***********************************************************************/

/* v.1.0.1 --> */
/* Program version */
static char *VERSION = "v.1.0.1";
/* <-- v.1.0.1 */

/* Array parameters */
#define FILENAME_LEN       512
#define LINELEN           1024
#define LONG_LINELEN     10000
#define MAX_CHAINS         100
#define MAX_CLEFT_CHAINS     4
#define MAX_RUN_CODES       10
#define MAX_TOKEN           10

/* Promotif plots */
#define BETA_BARRELS           0
#define SHEETS                 1
#define BAB_MOTIFS             2
#define BETA_HAIRPINS          3
#define PSI_LOOPS              4
#define BETA_BULGES            5
#define STRANDS                6
#define HELICES                7
#define HELIX_INTERACTIONS     8
#define BETA_TURNS             9
#define GAMMA_TURNS           10
#define DISULPHIDES           11

#define NPLOTS                12

static char *PLOT_DESC[NPLOTS] = {
  "beta_barrel", "sheet", "beta_alpha_beta_unit", "beta_hairpin",
  "psi_loop", "beta_bulge", "strands", "helices", "helint",
  "bturns", "gturns", "disulphide"
};

static char *PLOT_NAME[NPLOTS] = {
  "BETA_BARRELS", "SHEETS", "BAB_MOTIFS", "BETA_HAIRPINS", "PSI_LOOPS",
  "BETA_BULGES", "STRANDS", "HELICES", "HELIX_INTERACTIONS",
  "BETA_TURNS", "GAMMA_TURNS", "DISULPHIDES"
};

static char *HAVE_PS[NPLOTS] = {
  "No", "No", "No", "Yes", "No",
  "Yes", "No", "Yes", "No",
  "Yes", "Yes", "No"
};

static int MOTIF_NO[NPLOTS] = { -1, 7, 5, 3, 6, 2, -1, 4, 8, 0, 1, -1 };

/* Interaction types */
#define PP_HBONDS             0
#define PP_CONTACTS           1
#define PP_SALT_BRIDGES       2
#define PP_DISULPHIDES        3

#define NPP_TYPES             4

static char *LETTERS = "abcdefghijklmnopqrstuvwxyz";
static char *NUMBERS = "0123456789";

/* Commands set according to operating system */
char CP_COMMAND[6];
char MV_COMMAND[6];
char RM_COMMAND[6];
char SLEEP_COMMAND[30];

/* Structure for input parameters
   ------------------------------ */
struct parameters
{
/* v.1.0.1 --> */
  char               *version;
/* <-- v.1.0.1 */
  char               *file;
  char               pdb_code[5];
  int                single_file;
  int                max_cleft_chains;
  int                run_romlas;
  int                relinks;
};

/* Structure for PDB records
   ------------------------- */
struct pdb
{
  char              pdb_code[5];
  char              *file_name;
  char              *title;
  struct pdb        *next_pdb_ptr;
};

/* Structure for Promotif plots
   ---------------------------- */
struct plot
{
  char              chain;
  int               chain_no;
  int               count[NPLOTS];
  struct plot       *next_plot_ptr;
};

/* Structure for each chain record
   ------------------------------- */
struct chain
{
  char                    chain;
  int                     number;
  int                     chain_no;
  int                     nifaces;
  struct iface            *first_iface_ptr;
  struct iface            *last_iface_ptr;
  struct chain            *next_chain_ptr;
};

/* Structure for each iface record
   ----------------------------------- */
struct iface
{
  int                     iface_no;
  int                     nresid[2];
  double                  chain_area[2];
  double                  iface_area[2];
  int                     count[NPP_TYPES];
  int                     ntypes;
  struct chain            *chain_ptr;
  struct iface            *next_iface_ptr;
};

/* Structure for each ligmet record
   ----------------------------------- */
struct ligmet
{
  char                    *id;
  char                    *ligand_seq;
  char                    *het_name;
  char                    *resname1;
  char                    *res1;
  char                    *resname2;
  char                    *res2;
  char                    *chain;
  char                    *colour;
  int                     from_dic;
  int                     natoms;
  struct resid            *first_resid_ptr;
  struct ligmet           *next_ligmet_ptr;
};

/* Structure for each resid record
   ------------------------------- */
struct resid
{
  char                    *id;
  struct resid            *next_resid_ptr;
};

/* Structure for run record
   ------------------------ */
struct run
{
  char              *input_line;
  struct run        *next_run_ptr;
};

