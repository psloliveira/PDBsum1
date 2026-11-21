/***********************************************************************

  seqdata.h - Include file for seqdata.c

***********************************************************************/

/* Array parameters */
#define FILENAME_LEN       512
#define GO_LENGTH            7
#define LINELEN           1024
#define MAX_ANNOT           10
#define MAX_CHAINS         100
#define MAX_CODES          100
#define MAX_COPIES          10
#define MAX_LINES           10
#define MAX_TOKEN           30
#define NAMINO              20
#define NAMLEN             512
#define TOKEN_LEN          512
#define TUPLEN              20
#define UNIPROT_LEN         20

FILE *dump_file_ptr;

/* Code types */
#define UNKNOWN             -1
#define SWISSPROT_ID         0
#define SWISSPROT_CODE       1
#define EC_NUMBER            2
#define GO_COMPONENT         3
#define GO_FUNCTION          4
#define GO_PROCESS           5

#define NTYPES               6

/* Code sources */
#define MSD                  0
#define PDB_FILE             1
#define GO_FILE              2

/* Structure for input parameters
   ------------------------------ */
struct parameters
{
  char                    pdb_code[5];
  int                     uploaded_pdb;
  char                    file_name[FILENAME_LEN + 1];
  int                     pdb_type;
  int                     run_64;
};

/* Structure for each PDB entry
   ---------------------------- */
struct pdb
{
  char               pdb_code[5];
  int                nchains;
  int                nunique;
  int                ncodes[NTYPES];
  int                have_ec;
  struct chain       *first_chain_ptr;
  struct code        *first_code_ptr;
  struct code        *last_code_ptr;
  struct pdb         *next_pdb_ptr;
};

/* Structure for each chain data item
   ---------------------------------- */
struct chain
{
  char               chain_id;
  char               copy_of;
  int                ncodes;
  struct code        *code_ptr[MAX_CODES];
  struct chain       *next_chain_ptr;
};

/* Structure for each code data item
   --------------------------------- */
struct code
{
  int                type;
  char               *value;
  char               chain_id;
  int                source;
  int                count;
  char               *annotation[MAX_ANNOT];
  int                nannot;
  struct go          *go_ptr;
  struct code        *next_code_ptr;
};

/* Structure for each go data item
   --------------------------------- */
struct go
{
  char               go_id[GO_LENGTH + 1];
  int                type;
  char               *name;
  float              ave_level;
  int                level;
  struct go          *next_go_ptr;
};


