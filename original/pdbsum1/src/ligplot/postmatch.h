/***********************************************************************

  postmatch.h - Include file for postmatch.c

***********************************************************************/

#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <ctype.h>

/* Flag-values */     
#define FALSE                0
#define TRUE                 1

/* Array parameters */
#define FILENAME_LEN       512
#define LINELEN           1024
#define MAX_TOKEN          100
#define MAX_TYPES          100
#define STRING_LEN        1000
#define TOKEN_LEN          100

/* Minimum similarity for match to be considered a match */
#define MIN_SIMILARITY       40.0

/* Minimum number of atoms for molecule to be a significant one */
#define MIN_ATOMS             6

/* Molecule types */
#define REACTANT              0
#define PRODUCT               1
#define COFACTOR              2

#define NTYPES                3

FILE *dump_file_ptr;

/* Structure for file and directory names
   -------------------------------------- */
struct file_data
{
  char              *pdbsum_dir;
};

/* Structure for each match record
   ------------------------------- */
struct match
{
  char              pdb_code[5];
  int               ec_no;
  int               rstep_no;
  int               component_no;
  char              molec_id[7];
  int               mol_type;
  char              *mol_name;
  char              res_name[4];
  char              res_num[6];
  char              chain;
  double            similarity;
  int               nmatched;
  int               nmissed;
  int               nmissing;
  int               liguniq_no;
  int               ligand_no;
  int               wanted;
  int               written;
  struct match      *next_match_ptr;
};
