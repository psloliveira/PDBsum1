/***********************************************************************

  promfilt.h - Include file for promfilt.c

***********************************************************************/

/* Array parameters */
#define FILENAME_LEN       512
#define LINELEN           1024
#define MAX_CHAINS         100
#define MAX_COPIES          10
#define NAMLEN             512

#define NMOTIF_TYPES        12


/* Structure for input parameters
   ------------------------------ */
struct parameters
{
  char                    pdb_code[5];
  int                     pdb_type;
  int                     pdbsum1;
  int                     run_64;
};

/* Structure for each chain data item
   ---------------------------------- */
struct chain
{
  char               chain_id;
  int                chain_no;
  char               copy_of;
  struct chain       *next_chain_ptr;
};

/* Structure for each line data item
   ---------------------------------- */
struct line
{
  char               *line;
  struct line        *next_line_ptr;
};
