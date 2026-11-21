/***********************************************************************

  seqdata.c - Program to get SWISS-PROT code and E.C. number from data
              extracted from MSD database, or to use values in the PDB
              file if not present

************************************************************************

Date:-         12 Feb 2004
Dir update:-   19 Nov 2012

Written by:-   Roman Laskowski

----------------------------------------------------------------------*/

/* Datafiles
   ---------

Actual name                  Description
-----------                  -----------
AFanalyses.tsv               Input file of analyses of AlphaFold models
CATHPARAM                    Input parameter file
chains.dat                   Input chains file listing unique and
                             duplicate protein chains
pdb_chain_sp_ec              Input data extracted from the MSD
                             database giving PDB code, chain-id
                             and corresponding SWISS-PROT codes
                             and E.C. numbers (where present)
gene_association.goa_pdb.gz  Input PDB->GO annotations file
pdb_ec                       Output file listing E.C. number for
                             each chain for use by program enzyme
seqinfo.dat                  Output database file containing
                             extracted SWISS-PROT codes and E.C.
                             numbers
seqdata.dat                  Output file of summary information

Compilation
-----------

cc -c seqdata.c
cc -c common.c
cc -c reapar.c
cc -c readpdb.c
cc -o seqdata seqdata.o common.o reapar.o readpdb.o -lm -lz

Optimization flags:-  

cc -O2 -funroll-loops -c seqdata.c
cc -O2 -funroll-loops -c common.c
cc -O2 -funroll-loops -c reapar.c
cc -O2 -funroll-loops -c readpdb.c
cc -o seqdata seqdata.o common.o reapar.o readpdb.o -lm -lz

For testing on laptop, without read of gzipped files, use:
gcc -DTEST -c seqdata.c


*/

/* Function calling-tree
   ---------------------

main
   -> get_command_arguments
   -> c_get_paths (in common.c)
         -> get_param (in reapar.c)
         -> store_string (in common.c)
   -> get_pdb_codes
         -> string_truncate (in common.c)
         -> create_pdb_entry
   -> create_pdb_index
         -> pdb_sort
   -> get_msd_data
         -> read_chain_equivs
               -> string_chomp (in common.c)
               -> create_chain_record
         -> read_msd_data
               -> string_chomp (in common.c)
               -> get_tokens (in common.c)
               -> binary_pdb_search
               -> string_strip
               -> store_string (in common.c)
               -> create_code_record
                     -> store_string (in common.c)
               -> free_tokens (in common.c)
         -> get_pdb_data
               -> r_find_pdb_file (in readpdb.c)
               -> r_open_pdb_file (in readpdb.c)
               -> string_chomp (in common.c)
               -> store_string (in common.c)
               -> create_code_record
                     -> store_string (in common.c)
         -> get_af_ec
               -> string_chomp (in common.c)
               -> get_tokens (in common.c)
               -> store_string (in common.c)
               -> create_code_record
                     -> store_string (in common.c)
               -> free_tokens (in common.c)
   -> get_go_data
         -> to_lower (in common.c)
         -> to_upper (in common.c)
         -> get_tokens (in common.c)
         -> create_code_record
               -> store_string (in common.c)
         -> free_tokens (in common.c)
   -> get_go_names
         -> string_truncate (in common.c)
         -> create_go_record
   -> create_go_index
         -> go_sort
   -> match_go
         -> binary_go_search
   -> get_go_levels
          -> read_ontology_file
                -> binary_go_search
          -> sort_go_terms
          -> renumber_levels
   -> write_seqinfo
         -> open_output_file (in common.c)

*/


#include "common.h"
#include "seqdata.h"


/* Prototypes
   ---------- */

/* In common.c */
void c_get_paths(struct c_file_data *file_name_ptr,int run_64);
int find_pdbsum_dir(char *pdb_code,char *pdbsum_template[NPDBSUM_DIR],
                    int pdb_type,char *dir_name);
void free_tokens(char **token,int ntokens);
int get_tokens(char *line,char **token,int max_token,char token_sep);
FILE *open_output_file(char *out_name,int terminate);
void store_string(char **string_ptr,char *string);
int string_chomp(char *string);
int string_truncate(char *string,int max_length);
int tokenize(char *string,char **token,int max_tokens,char token_sep);
void to_lower(char *search_string,int length);
void to_upper(char *search_string,int length);

/* Function prototypes from reapar.c */
int get_param(char *token_name,char *file_name,int filename_len,
              int exit_on_error);

/* Function prototypes from readpdb.c */
int r_find_pdb_file(char *pdb_code,char *pdb_template[NPDB_DIR],
                    char *file_name,int pdb_type);
FILE *r_open_pdb_file(char *file_name,int *gz);

/***********************************************************************

get_command_arguments  -  Interpret the command line arguments

***********************************************************************/

void get_command_arguments(char *strptrs[],int num,
                           struct parameters *params)
{
  int itoken;

  /* Initialise parameters */
  params->pdb_code[0] = '\0';
  params->pdb_type = STANDARD_PDB;
  params->run_64 = FALSE;

  /* If -halp option entered, show options available */
  if (num == 1 && (!strcmp(strptrs[1],"-help") ||
                   !strcmp(strptrs[1],"-h")))
    {
      printf("\n");
      printf("Usage is ...\n"
             "\n"
             "  seqdata  [UPLOAD | MCSG |  AF]  [-pdbsum1 pdb_code]  [-64]\n"
             "\n"
             "where\n"
             "   * UPLOAD            = Uploaded PDB file for PDBsum "
             "generation\n"
             "   * MCSG              = MCSG structure\n"
             "   * AF                = Alpha Fold structure\n"
             "   * -pdbsum1 pdb_code = PDBsum1 structure\n"
             "   * -64               = run on 64-bit processor\n"
             "\n"
             "\n");
      exit(1);
    }

  /* Loop through the command arguments */
  for (itoken = 1; itoken < num + 1; itoken++)
    {
      /* Check for the UPLOAD option */
      if (!strcmp(strptrs[itoken],"UPLOAD"))
        params->pdb_type = UPLOAD_PDB;

      /* Check for the MCSG option */
      else if (!strcmp(strptrs[itoken],"MCSG"))
        params->pdb_type = MCSG_PDB;

      /* Check for the Alpha Fold option */
      else if (!strcmp(strptrs[itoken],"AF"))
        params->pdb_type = ALPHA_FOLD_PDB;

      /* Check for the PDBsum1 option */
      else if (!strcmp(strptrs[itoken],"-pdbsum1"))
        {
          /* Get the PDB code from the next token */
          if (itoken < num)
            {
              /* Get the PDB code */
              strcpy(params->pdb_code,strptrs[itoken + 1]);

              /* Increment token count */
              itoken++;
            }

          /* Set flag */
          params->pdb_type = PDBSUM1;
        }

      /* Check for the 64-bit option */
      else if (!strcmp(strptrs[itoken],"-64"))
        params->run_64 = TRUE;
    }
}
/***********************************************************************

create_pdb_entry  -  Create a pdb entry

***********************************************************************/

struct pdb *create_pdb_entry(struct pdb **fst_pdb_ptr,
                             struct pdb **lst_pdb_ptr,char *pdb_code)
{
  int i, isort;

  struct pdb *pdb_ptr;
  struct pdb *first_pdb_ptr, *last_pdb_ptr;

  /* Initialise pdb pointers */
  pdb_ptr = NULL;
  first_pdb_ptr = *fst_pdb_ptr;
  last_pdb_ptr = *lst_pdb_ptr;

  /* Create pdb entry */
  pdb_ptr = (struct pdb *) malloc(sizeof(struct pdb));
  if (pdb_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for struct pdb\n");
      printf("***        Program terminated with error.\n");
      exit (1);
    }

  /* Store the PDB code */
  strcpy(pdb_ptr->pdb_code,pdb_code);

  /* Initialise remaining fields */
  pdb_ptr->nchains = 0;
  pdb_ptr->nunique = 0;
  for (i = 0; i < NTYPES; i++)
    pdb_ptr->ncodes[i] = 0;
  pdb_ptr->have_ec = FALSE;
  pdb_ptr->first_chain_ptr = NULL;
  pdb_ptr->first_code_ptr = NULL;
  pdb_ptr->last_code_ptr = NULL;

  /* Initialise pointer to next entry */
  pdb_ptr->next_pdb_ptr = NULL;

  /* If this is the first entry stored, then update
     pointer to head of linked list */
  if (first_pdb_ptr == NULL)
    first_pdb_ptr = pdb_ptr;

  /* Otherwise, make previous entry point to the current
     one */
  else
    last_pdb_ptr->next_pdb_ptr = pdb_ptr;
                      
  /* Save the current residue's pointer */
  last_pdb_ptr = pdb_ptr;

  /* Return current pointers */
  *fst_pdb_ptr = first_pdb_ptr;
  *lst_pdb_ptr = last_pdb_ptr;
  return(pdb_ptr);
}
/***********************************************************************

get_pdb_codes  -  Read through the dataset.lst file to read in the PDB
                  codes to be processed

***********************************************************************/

int get_pdb_codes(char *dataset_name,struct pdb **fst_pdb_ptr,
                  char *single_pdb_code)
{
  char pdb_code[5];
  char file_name[FILENAME_LEN], input_line[LINELEN + 1];

  int len, ncodes;

  struct pdb *pdb_ptr;
  struct pdb *first_pdb_ptr, *last_pdb_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  pdb_ptr = NULL;
  first_pdb_ptr = last_pdb_ptr = NULL;
  ncodes = 0;

  /* If have a single PDB code, then save that and return */
  if (single_pdb_code[0] != '\0')
    {
      /* Create a new pdb entry */
      pdb_ptr = create_pdb_entry(&first_pdb_ptr,&last_pdb_ptr,
                                 single_pdb_code);

      /* Increment count */
      ncodes++;

      /* Return pointer to first pdb entry */
      *fst_pdb_ptr = first_pdb_ptr;

      /* Return the number of sequences read in */
      return(ncodes);
    }

  /* Form name of dataset.lst file */
  strcpy(file_name,dataset_name);

  /* Open the dataset.lst file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    {
      printf("\n*** ERROR. Unable to open file [%s]\n",file_name);
      exit(1);
    }

  /* Loop through the file, line by line */
  while (fgets(input_line,LINELEN,file_ptr) != NULL)
    {
      /* Truncate the input line */
      len = string_truncate(input_line,LINELEN);

      /* Check line-length */
      if (len >= 4)
        {
          /* Get the PDB code */
          strncpy(pdb_code,input_line,4);
          pdb_code[4] = '\0';

          /* Create a new pdb entry */
          pdb_ptr = create_pdb_entry(&first_pdb_ptr,&last_pdb_ptr,
                                     pdb_code);

          /* Increment count */
          ncodes++;
        }
    }

  /* Close file */
  fclose(file_ptr);

  /* Return pointer to first pdb entry */
  *fst_pdb_ptr = first_pdb_ptr;

  /* Return the number of sequences read in */
  printf("Number of PDB codes read in: %d\n",ncodes);
  return(ncodes);
}
/***********************************************************************

pdb_sort  -  Shell-sort which sorts the UniProt pdbs into ascending
              order

***********************************************************************/

void pdb_sort(struct pdb **pdb_index_ptr,int nindex)
{
  char pdb_name1[5], pdb_name2[5];

  int endloop, i, indi, indl, j, k, l, lognb2, m, n, nn;
  float aln2i = (float) 1.4426950, tiny = (float) 1.e-5;
  float calc;
  double dble, dble1;

  struct pdb *pdb1_ptr, *pdb2_ptr, *swap_pdb_ptr;

  /*--Initialise values */
  n = nindex;
  dble = n;
  dble1 = log(dble);
  calc = (float) dble1;
  calc = calc * aln2i + tiny;
  lognb2 = (int) calc;

  /* Perform sort */
  m = n;
  for (nn = 1; nn <= lognb2; nn++)
    {
      m = m / 2;
      k = n - m;
      for (j = 1; j <= k; j++)
        {
          i = j;
          endloop = FALSE;
          while (endloop == FALSE)
            {
              l = i + m;
              indi = i - 1;
              indl = l - 1;

              /* Get the two pdb records */
              pdb1_ptr = pdb_index_ptr[indi];
              pdb2_ptr = pdb_index_ptr[indl];

              /* Form full pdb names on which sort is based */
              strcpy(pdb_name1,pdb1_ptr->pdb_code);
              strcpy(pdb_name2,pdb2_ptr->pdb_code);

              /* If in wrong order, the swap sort-pointers */
              if (strcmp(pdb_name1,pdb_name2) > 0)
                {
                  swap_pdb_ptr = pdb_index_ptr[indi];
                  pdb_index_ptr[indi] = pdb_index_ptr[indl];
                  pdb_index_ptr[indl] = swap_pdb_ptr;

                  i = i - m;
                  if (i < 1)
                    endloop = TRUE;
                }
              else
                endloop = TRUE;
            }
        }
    }
}
/***********************************************************************

create_pdb_index  -  Create index of PDB records for use in binary
                     search

***********************************************************************/

void create_pdb_index(struct pdb *first_pdb_ptr,
                      struct pdb ***pdb_index_ptr,int npdb)
{
  int nentry;

  struct pdb *pdb_ptr;

  struct pdb **pdb_ind_ptr;

  /* Initialise variables */

  /* Create the index array as one large array */
  pdb_ind_ptr
    = (struct pdb **) malloc(npdb * sizeof(struct pdb *));

  /* Check that sufficient memory was allocated */
  if (pdb_ind_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for index array\n");
      exit(1);
    }

  /* Get pointer to the first pdb entry */
  pdb_ptr = first_pdb_ptr;
  nentry = 0;

  /* Loop through all entries to place each one in the sorted index
     array */
  while (pdb_ptr != NULL)
    {
      /* Store current pointer at end of the index array */
      if (nentry < npdb)
        pdb_ind_ptr[nentry] = pdb_ptr;

      /* Increment count of stored entries */
      nentry++;

      /* Get pointer to the next NIST entry in the linked list */
      pdb_ptr = pdb_ptr->next_pdb_ptr;
    }

  /* Sort the code pointers into order of UniProt code */
  if (npdb > 1)
    pdb_sort(pdb_ind_ptr,npdb);

  /* Return the array pointer */
  *pdb_index_ptr = pdb_ind_ptr;
}
/***********************************************************************

binary_search  -  Search for the given PDB record

***********************************************************************/

struct pdb *binary_search(struct pdb **pdb_index_ptr,int npdb,
                          char *pdb_code)
{
  int bottom_point, ipoint, new_point, top_point;
  int found;

  struct pdb *pdb_ptr, *found_pdb_ptr, *prev_pdb_ptr;

  /* Initialise variables */
  found = FALSE;
  found_pdb_ptr = NULL;

  /* Make starting point half-way into the index */
  ipoint = npdb / 2;
  bottom_point = 0;
  top_point = npdb - 1;

  /* Loop until required index position found */
  while (found == FALSE)
    {
      /* Get the associated GO pointer */
      pdb_ptr = pdb_index_ptr[ipoint];

      /* Check if index value is equal to the value we want */
      if (!strcmp(pdb_code,pdb_ptr->pdb_code))
        {
          /* Set flag to indicate that match found and store pointer */
          found = TRUE;
          found_pdb_ptr = pdb_ptr;
        }

      /* If index value is higher than the one we want, then need
         to search lower */
      else if (strcmp(pdb_code,pdb_ptr->pdb_code) < 0)
        {
          /* If this is the first go in the sorted list, then can't
             go any lower */
          if (ipoint == 0)
            {
              /* Set flag that we are done */
              found = TRUE;
            }

          /* Otherwise, check whether the previous entry's name
             is below the one we're searching for */
          else
            {
              /* Get the pointer to the previous entry */
              prev_pdb_ptr = pdb_index_ptr[ipoint - 1];

              /* If value is below the search value, we've found the
                 position in the index we're after, but don't have an exact
                 match*/
              if (strcmp(prev_pdb_ptr->pdb_code,pdb_code) < 0)
                found = TRUE;

              /* Otherwise, need to search earlier in the list */
              else
                {
                  /* Save the current position as the uppermost position
                     to search on */
                  top_point = ipoint;

                  /* Set pointer to halfway between the last search
                     low-search and here */
                  ipoint = (ipoint + bottom_point) / 2;
                }
            }
        }

      /* Otherwise, have to look further up the list */
      else
        {
          /* Save the current position as the lowest position to search
             on */
          bottom_point = ipoint;

          /* If have reached the top of the list, then our value is
             higher than any value in the list */
          if (ipoint >= npdb - 1)
            {
              /* Set flag as done */
              found = TRUE;
            }

          /* Otherwise, can search further up the list */
          else
            {
              /* Set pointer to halfway between the last search high-search
                 and here */
              new_point = (ipoint + top_point) / 2;

              /* Adjust if going back to the same point */
              if (new_point == ipoint)
                new_point++;
              ipoint = new_point;
            }
        }
    }

  /* Return the GO pointer */
  return(found_pdb_ptr);
}
/***********************************************************************

string_strip  -  Remove any leading and trailing spaces from the given
                 string

***********************************************************************/

int string_strip(char *string,int max_length)
{
  char ch;

  int iend, ipos, nblanks;
  int blank;

  /* Initialise variables */
  blank = TRUE;
  nblanks = 0;
  iend = -1;
  ch = string[0];
  for (ipos = 0; ipos < max_length && ch != '\0'; ipos++)
    {
      /* Get the current character and check for end of string */
      ch = string[ipos];
      if (ch != '\0')
        {
          /* If not a space, then mark end of string */
          if (ch != ' ')
            {
              iend = ipos;
              blank = FALSE;
            }

          /* If a leading space, then shift whole string */
          else if (blank == TRUE)
            nblanks++;
        }
    }

  /* End the string after the last non-blank character */
  if (iend + 1 < max_length)
    string[iend + 1] = '\0';
  else
    string[max_length - 1] = '\0';

  /* If have any leading blanks, remove these */
  if (nblanks > 0)
    strcpy(string,string+nblanks);

  /* Return the string length */
  return(strlen(string));
}
/***********************************************************************

create_chain_record  -  Create and initialise a new chain record

***********************************************************************/

struct chain *create_chain_record(struct chain **fst_chain_ptr,
                                  struct chain **lst_chain_ptr,
                                  char chain_id,char copy_of)
{
  int i;

  struct chain *chain_ptr, *first_chain_ptr, *last_chain_ptr;

  /* Initialise chain pointers */
  chain_ptr = NULL;
  first_chain_ptr = *fst_chain_ptr;
  last_chain_ptr = *lst_chain_ptr;

  /* Allocate memory for structure to hold chain info */
  chain_ptr = (struct chain *) malloc(sizeof(struct chain));
  if (chain_ptr == NULL)
    {
      printf("*** Can't allocate memory for struct chain\n");
      exit (1);
    }

  /* If this is the very first chain, then store pointer */
  if (first_chain_ptr == NULL)
    first_chain_ptr = chain_ptr;

  /* Add link from previous chain to the current one */
  if (last_chain_ptr != NULL)
    last_chain_ptr->next_chain_ptr = chain_ptr;
  last_chain_ptr = chain_ptr;

  /* Store the chain details */
  chain_ptr->chain_id = chain_id;
  chain_ptr->copy_of = copy_of;
  chain_ptr->next_chain_ptr = NULL;

  /* Initialise code pointers */
  chain_ptr->ncodes = 0;
  for (i = 0; i < MAX_CODES; i++)
    chain_ptr->code_ptr[i] = NULL;

  /* Return current pointers */
  *fst_chain_ptr = first_chain_ptr;
  *lst_chain_ptr = last_chain_ptr;

  /* Return new chain pointer */
  return(chain_ptr);
}
/***********************************************************************

read_chain_equivs  -  Read in the chain equivalences from the chains.dat
                      file

***********************************************************************/

int read_chain_equivs(char *pdb_code,char *pdbsum_dir,
                      struct chain **fst_chain_ptr)
{
  char file_name[FILENAME_LEN + 1], input_line[LINELEN + 1];
  char chain_id, copy_of;

  int len, nchains, nlines;

  struct chain *chain_ptr, *first_chain_ptr, *last_chain_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  nchains = 0;
  nlines = 0;

  /* Initialise chain pointers */
  chain_ptr = first_chain_ptr = last_chain_ptr = NULL;

  /* Form the name of the chains.dat file */
  strcpy(file_name,pdbsum_dir);
  strcat(file_name,"chains.dat");
  // printf("Reading chains.dat file (%s) ...\n",file_name);

  /* Open the chains.dat input file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    {
      //      printf("\n*** Warning. Unable to open input file [%s]\n",file_name);
    }

  /* If OK, read through the file v*/
  else
    {
      /* Loop while reading in records from the input file */
      while (fgets(input_line,LINELEN,file_ptr)!= NULL)
        {
          /* Increment count of lines read in */
          nlines++;

          /* Truncate the string to strip off the newline */
          // string_truncate(input_line,LINELEN);
          len = string_chomp(input_line);

          /* If line is blank, then take to be chain blank */
          if (len == 0)
            chain_id = ' ';
          else
            chain_id = input_line[0];

          /* If this is a copy of another chain, store that */
          if (len > 2)
            copy_of = input_line[2];
          else
            copy_of = '\0';

          /* Create a chain record */
          chain_ptr
            = create_chain_record(&first_chain_ptr,&last_chain_ptr,
                                  chain_id,copy_of);

          /* Increment chain count */
          nchains++;
        }

      /* Close the input file */
      fclose(file_ptr);
    }

  /* Sho stats */
  // printf("  Number of chains defined           %8d\n",nchains);
  // printf("\n");

  /* Return pointer to first chain record */
  *fst_chain_ptr = first_chain_ptr;

  /* Return the number of chains stored */
  return(nchains);
}
/***********************************************************************

create_code_record  -  Create a new code record or update existing one
                       if we already have it

***********************************************************************/

struct code *create_code_record(struct code **fst_code_ptr,
                                struct code **lst_code_ptr,
                                int type,char *value,int source,
                                char *annotation,int *new,char chain)
{
  int iannot;
  int found;

  struct code *code_ptr, *match_code_ptr;
  struct code *first_code_ptr, *last_code_ptr;

  /* Initialise variables */
  found = FALSE;
  *new = FALSE;

  /* Initialise code pointers */
  last_code_ptr = *lst_code_ptr;
  first_code_ptr = *fst_code_ptr;
  match_code_ptr = NULL;

  /* Get pointer to first entry */
  code_ptr = first_code_ptr;

  /* Loop through all these to determine whether we already have this
     entry */
  while (code_ptr != NULL && found == FALSE)
    {
      /* If have match on name, then already have this entry */
      if (code_ptr->type == type && !strcmp(code_ptr->value,value) &&
          (code_ptr->chain_id == chain ||
           (type != SWISSPROT_ID && type != SWISSPROT_CODE)))
        {
          /* Store pointer to located entry */
          match_code_ptr = code_ptr;

          /* Increment count of chains belonging to this code */
          code_ptr->count++;

          /* Set flag that found */
          found = TRUE;
        }

      /* Get next code entry */
      code_ptr = code_ptr->next_code_ptr;
    }

  /* If match found, then retrieve pointer */
  if (found == TRUE)
    code_ptr = match_code_ptr;

  /* Otherwise, create a new code entry */
  else
    {
      code_ptr = (struct code*)malloc(sizeof(struct code));
      if (code_ptr == NULL)
        {
          printf("*** ERROR. Unable to allocate memory for");
          printf(" struct code\n");
          printf("***        Program seqdata terminated with error.\n");
          exit (1);
        }

      /* Set flag to indicate that this is a new entry */
      *new = TRUE;

      /* Store the details for this entry */
      code_ptr->type = type;
      code_ptr->chain_id = chain;
      store_string(&(code_ptr->value),value);
      code_ptr->source = source;
      for (iannot = 0; iannot < MAX_ANNOT; iannot++)
        code_ptr->annotation[iannot] = &null;
      code_ptr->nannot = 0;
      code_ptr->go_ptr = NULL;

      /* Initialise count */
      code_ptr->count = 1;

      /* Initialise pointer to next code record */
      code_ptr->next_code_ptr = NULL;

      /* If this is the first entry stored, then update pointer to
         head of linked list */
      if (first_code_ptr == NULL)
        first_code_ptr = code_ptr;

      /* Otherwise, make previous entry point to the current
         one */
      else
        last_code_ptr->next_code_ptr = code_ptr;
                      
      /* Save the current code's pointer */
      last_code_ptr = code_ptr;
    }

  /* If have an annotation, add to list */
  if (annotation[0] != '\0' && code_ptr->nannot < MAX_ANNOT)
    {
      /* Initialise flag */
      found = FALSE;

      /* Loop to see whather we already have this annotation */
      for (iannot = 0; iannot < code_ptr->nannot; iannot++)
        if (!strcmp(code_ptr->annotation[iannot],annotation))
          found = TRUE;

      /* If don't already have it, store in next available slot */
      if (found == FALSE)
        {
          iannot = code_ptr->nannot;
          store_string(&(code_ptr->annotation[iannot]),annotation);
          code_ptr->nannot++;
        }
    }

  /* Return current pointers */
  *fst_code_ptr = first_code_ptr;
  *lst_code_ptr = last_code_ptr;
  return(code_ptr);
}
/***********************************************************************

binary_pdb_search  -  Search for the given PDB file

***********************************************************************/

struct pdb *binary_pdb_search(struct pdb **pdb_index_ptr,int ncodes,
                              char *pdb_code)
{
  int bottom_point, ipoint, new_point, top_point;
  int found;

  struct pdb *pdb_ptr, *found_pdb_ptr, *prev_pdb_ptr;

  /* Initialise variables */
  found = FALSE;
  found_pdb_ptr = NULL;

  /* Make starting point half-way into the index */
  ipoint = ncodes / 2;
  bottom_point = 0;
  top_point = ncodes - 1;

  /* Loop until required index position found */
  while (found == FALSE)
    {
      /* Get the associated GO pointer */
      pdb_ptr = pdb_index_ptr[ipoint];

      /* Check if index value is equal to the value we want */
      if (!strcmp(pdb_code,pdb_ptr->pdb_code))
        {
          /* Set flag to indicate that match found and store pointer */
          found = TRUE;
          found_pdb_ptr = pdb_ptr;
        }

      /* If index value is higher than the one we want, then need
         to search lower */
      else if (strcmp(pdb_code,pdb_ptr->pdb_code) < 0)
        {
          /* If this is the first go in the sorted list, then can't
             go any lower */
          if (ipoint == 0)
            {
              /* Set flag that we are done */
              found = TRUE;
            }

          /* Otherwise, check whether the previous entry's name
             is below the one we're searching for */
          else
            {
              /* Get the pointer to the previous entry */
              prev_pdb_ptr = pdb_index_ptr[ipoint - 1];

              /* If value is below the search value, we've found the
                 position in the index we're after, but don't have an exact
                 match*/
              if (strcmp(prev_pdb_ptr->pdb_code,pdb_code) < 0)
                found = TRUE;

              /* Otherwise, need to search earlier in the list */
              else
                {
                  /* Save the current position as the uppermost position
                     to search on */
                  top_point = ipoint;

                  /* Set pointer to halfway between the last search
                     low-search and here */
                  ipoint = (ipoint + bottom_point) / 2;
                }
            }

        }

      /* Otherwise, have to look further up the list */
      else
        {
          /* Save the current position as the lowest position to search
             on */
          bottom_point = ipoint;

          /* If have reached the top of the list, then our value is
             higher than any value in the list */
          if (ipoint >= ncodes - 1)
            {
              /* Set flag as done */
              found = TRUE;
            }

          /* Otherwise, can search further up the list */
          else
            {
              /* Set pointer to halfway between the last search high-search
                 and here */
              new_point = (ipoint + top_point) / 2;

              /* Adjust if going back to the same point */
              if (new_point == ipoint)
                new_point++;
              ipoint = new_point;
            }
        }
    }

  /* Return the matched pointer */
  return(found_pdb_ptr);
}
/***********************************************************************

read_msd_data  -  Read in the given PDB file

***********************************************************************/

int read_msd_data(struct pdb **pdb_index_ptr,int npdb,
                  char *pdbsum_data_dir)
{
  char file_name[LINELEN + 1], input_line[LINELEN + 1], pdb_code[5];
  char chain, code[3][TOKEN_LEN], copy_of;
  char annotation[TOKEN_LEN];

  char *token[MAX_TOKEN];

  int i, icode, len, ncodes, nitems, nlines, ntokens, type;
  int found, have_code, new;

  struct chain *chain_ptr, *match_chain_ptr;
  struct code *code_ptr;
  struct pdb *pdb_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  annotation[0] = '\0';
  found = FALSE;
  nitems = 0;
  nlines = 0;

  /* Form file name */
  strcpy(file_name,pdbsum_data_dir);
  strcat(file_name,DIR_SEP);
  strcat(file_name,"pdb_chain_sp_ec");

  /* Open the MSD data file */
  if ((file_ptr = fopen(file_name,"r")) !=  NULL)
    {
      /* Loop while reading in records from the input PDB file */
      while (fgets(input_line,LINELEN,file_ptr)!= NULL)
        {
          /* Truncate the string to strip off the newline */
          // len = string_truncate(input_line,LINELEN);
          len = string_chomp(input_line);

          /* Increment count of lines read in */
          nlines++;

          /* Retrieve all the tokens from the current line */
          // ntokens = get_tokens(input_line,token,MAX_TOKEN,'"');
          ntokens = tokenize(input_line,token,MAX_TOKEN,'"');

          /* Get the PDB code */
          strcpy(pdb_code,token[1]);

          /* Find this code */
          pdb_ptr = binary_pdb_search(pdb_index_ptr,npdb,pdb_code);

          /* If PDB code found, extract info on this line */
          if (ntokens > 8 && pdb_ptr != NULL)
            {
              /* Get the chain identifier */
              chain = token[3][0];
              if (chain == '\0' || chain == '@')
                chain = ' ';

              /* Get pointer to the first chain record */
              chain_ptr = pdb_ptr->first_chain_ptr;
              match_chain_ptr = NULL;

              /* Loop through all the chain records until find the correct
                 one */
              while (chain_ptr != NULL && match_chain_ptr == NULL)
                {
                  /* If this is the correct record, store the pointer */
                  if (chain_ptr->chain_id == chain)
                    match_chain_ptr = chain_ptr;

                  /* Get pointer to the next chain record */
                  chain_ptr = chain_ptr->next_chain_ptr;
                }

              /* If have a match, store the SWISS-PROT and E.C. codes */
              if (match_chain_ptr != NULL)
                {
                  /* Strip off the leading and trailing blanks, if any */
                  string_strip(token[5],TOKEN_LEN);
                  string_strip(token[7],TOKEN_LEN);
                  if (ntokens > 9)
                    string_strip(token[9],TOKEN_LEN);

                  /* Pull out the three codes of interest from the current
                     line's tokens */
                  strcpy(code[SWISSPROT_ID],token[5]);
                  strcpy(code[SWISSPROT_CODE],token[7]);
                  if (ntokens > 9)
                    strcpy(code[EC_NUMBER],token[9]);
                  else
                    code[EC_NUMBER][0] = '\0';

                  /* If this is a copy chain, determine which chain it is a
                     copy of */
                  if (match_chain_ptr->copy_of != '\0')
                    copy_of = match_chain_ptr->copy_of;
                  else
                    copy_of = match_chain_ptr->chain_id;

                  /* Loop over the different code types to store */
                  for (type = 0; type < 3; type++)
                    {
                      /* If have this code, then store */
                      if (code[type][0] != '\0')
                        {
                          /* Store the code itself */
                          code_ptr
                            = create_code_record(&(pdb_ptr->first_code_ptr),
                                                 &(pdb_ptr->last_code_ptr),
                                                 type,code[type],MSD,
                                                 annotation,&new,copy_of);

                          /* If this is a new code, increment count of codes */
                          if (new == TRUE)
                            pdb_ptr->ncodes[type]++;

                          /* If this is an E.C. number, set flag */
                          if (type == EC_NUMBER)
                            pdb_ptr->have_ec = TRUE;

                          /* Store for current chain */
                          if (match_chain_ptr->ncodes < MAX_CODES)
                            {
                              /* Determine whether this chain already has this
                                 code */
                              have_code = FALSE;

                              /* Loop over the stored codes */
                              for (icode = 0; icode < match_chain_ptr->ncodes;
                                   icode++)
                                {
                                  /* If have code, then set flag */
                                  if (match_chain_ptr->code_ptr[icode]
                                      == code_ptr)
                                    have_code = TRUE;
                                }

                              /* If this chain doesn't already have this code,
                                 then store */
                              if (have_code == FALSE)
                                {
                                  ncodes = match_chain_ptr->ncodes;
                                  match_chain_ptr->code_ptr[ncodes]
                                    = code_ptr;
                                  match_chain_ptr->ncodes++;
                                }
                            }
                        }
                    }

                  /* Increment item count */
                  nitems++;
                }

              /* Set flag that correct code found */
              found = TRUE;
            }

          /* Free the tokens */
          // free_tokens(token,ntokens);
        }

      /* Close the file */
      fclose(file_ptr);
    }
  else
    printf("\n*** Warning. Unable to open input file [%s]\n",
           file_name);

  /* Return the number of items found */
  return(nitems);
}
/***********************************************************************

get_pdb_data  -  Obtain SWISS-PROT code(s) and E.C. number(s) from the
                 PDB file

***********************************************************************/

int get_pdb_data(char *pdb_code,struct chain *first_chain_ptr,
                 struct code **fst_code_ptr,struct code **lst_code_ptr,
                 int *ncodes,int nitems,struct c_file_data file_name_ptr,
                 struct parameters params)
{
  char file_name[FILENAME_LEN], input_line[LINELEN + 1];
  char annotation[TOKEN_LEN], tmp_string[LINELEN + 1];
  char chain, code_value[2][UNIPROT_LEN], copy_of;

  char *string_ptr, *token[MAX_TOKEN];

  int dir, i, icode, len, nlines, type[2], ntokens, ntypes;
  int done, eof, gzipped, have_code, have_ec_data, have_sws_data, new;
  int next_sws_entry, unp_continued;

  struct code *code_ptr, *first_code_ptr, *last_code_ptr;
  struct chain *chain_ptr, *match_chain_ptr;

#ifdef TEST
  FILE *file_ptr;
#else
  FILE *file_ptr;
  gzFile gzfile_ptr;
#endif

  /* Initialise counts */
  annotation[0] = '\0';
  done = FALSE;
  eof = FALSE;
  next_sws_entry = FALSE;
  unp_continued = FALSE;
  nlines = 0;
  chain = '\0';
  if (ncodes[SWISSPROT_ID] + ncodes[SWISSPROT_CODE] > 0)
    have_sws_data = TRUE;
  else
    have_sws_data = FALSE;
  if (ncodes[EC_NUMBER] > 0)
    have_ec_data = TRUE;
  else
    have_ec_data = FALSE;

  /* Initialise pointers */
  code_ptr = NULL;
  first_code_ptr = *fst_code_ptr;
  last_code_ptr = *lst_code_ptr;

  /* Form the name of the PDB file */
  if (params.pdb_type == PDBSUM1)
    sprintf(file_name,"../%s.pdb",pdb_code);
  else
    dir = r_find_pdb_file(pdb_code,file_name_ptr.pdb_template,
                          file_name,params.pdb_type);

  /* If file not found, then return */
  if (file_name[0] == '\0' || dir == NOT_FOUND)
      return(nitems);

  /* Determine whether the PDB file is a gzipped one */
  gzipped = FALSE;
  len = strlen(file_name);
  if (!strncmp(file_name + len - 3,".gz",3))
    gzipped = TRUE;
#ifdef TEST
  gzipped = FALSE;
#endif

  /* Open PDB file */
#ifdef TEST
  file_ptr = fopen(file_name,"r");
  if (file_ptr ==  NULL)
    eof = TRUE;
#else
  if (gzipped == TRUE)
    {
      gzfile_ptr = gzopen(file_name,"r");
      if (gzfile_ptr ==  NULL)
        eof = TRUE;
    }
  else
    {
      file_ptr = fopen(file_name,"r");
      if (file_ptr ==  NULL)
        eof = TRUE;
    }
#endif

  /* If file not opened, then abort */
  if (eof == TRUE)
    {
      printf("\n*** Warning. Unable to open PDB file [%s]\n",file_name);
      return(nitems);
    }

  /* Loop through the PDB file, picking up all the atoms that are
     required */
  while (eof == FALSE && done == FALSE)
    {
      /* Read in the next record */
      if (gzipped == TRUE)
        {
#ifdef TEST
#else
          if (gzgets(gzfile_ptr,input_line,LINELEN) == NULL)
            eof = TRUE;
#endif
        }
      else
        {
          if (fgets(input_line,LINELEN,file_ptr) == NULL)
            eof = TRUE;
        }

      /* If record read in OK, then process it */
      if (eof == FALSE)
        {
          /* Initialise */
          ntypes = 0;

          /* Truncate the string to strip off the newline */
          // len = string_truncate(input_line,LINELEN);
          len = string_chomp(input_line);

          /* Increment count of lines read in */
          nlines++;

          /* If COMPND record, check for E.C. number */
          if (!strncmp("COMPND",input_line,6) && have_ec_data == FALSE)
            {
              /* Old-style entry: (E.C.i.j.k.l)*/
              if ((string_ptr = strstr(input_line,"(E.C.")) != NULL)
                {
                  /* Splice out the number */
                  strcpy(tmp_string,string_ptr+5);
                  string_ptr = strstr(tmp_string,")");
                  if (string_ptr != NULL)
                    string_ptr[0] = '\0';

                  /* Check for there being at least one dot in the number
                     and it not being too long */
                  len = strlen(tmp_string);
                  if (len < 20 && strstr(tmp_string,".") != NULL)
                    {
                      /* Store the E.C. number */
                      strcpy(code_value[0],tmp_string);

                      /* Update count */
                      ntypes = 1;
                      type[0] = EC_NUMBER;
                    }
                }

              /* Otherwise, check for new-style entry:  EC: i.j.k.l; */
              else if ((string_ptr = strstr(input_line,"EC: ")) != NULL)
                {
                  /* Splice out the number */
                  strcpy(tmp_string,string_ptr+4);
                  string_ptr = strstr(tmp_string,";");
                  if (string_ptr != NULL)
                    string_ptr[0] = '\0';

                  /* Check for there being at least one dot in the number
                     and it not being too long */
                  len = strlen(tmp_string);
                  if (len < 20 && strstr(tmp_string,".") != NULL)
                    {
                      /* Store the E.C. number */
                      strcpy(code_value[0],tmp_string);

                      /* Update count */
                      ntypes = 1;
                      type[0] = EC_NUMBER;
                    }
                }

              /* Otherwise, check for chain identifier */
              else if ((string_ptr = strstr(input_line," CHAIN: ")) != NULL)
                {
                  /* Store the first chain identifier */
                  chain = string_ptr[8];
                }
            }

          /* If DBREF record, check for SWISS-PROT code */
          if (!strncmp("DBREF",input_line,5) && have_sws_data == FALSE)
            {
              /* Check if this line contains a link to a SWISS-PROT code */
              if (((!strncmp(input_line+26,"SWS",3) ||
                   !strncmp(input_line+26,"UNP",3)) &&
                   strstr(input_line," SPAC ") == NULL) ||
                  (input_line[5] == '2' && unp_continued == TRUE))
                {
                  /* Get the chain-id */
                  chain = input_line[12];

                  /* Extract the SWISS-PROT id and accession number */
                  ntokens
                    = get_tokens(input_line + 29,token,MAX_TOKEN,' ');

                  /* If only have one token, then have long UniProt id,
                     with accession on next line */
                  if (ntokens == 1 && input_line[5] == '1')
                    {
                      /* Save the UniProt id */
                      len = strlen(token[0]);
                      if (len > 2 && len < UNIPROT_LEN)
                        {
                          /* Store code */
                          strcpy(code_value[ntypes],token[0]);
                          type[ntypes] = SWISSPROT_CODE;
                          ntypes++;
                        }

                      /* Set flag that accession is on the next line */
                      unp_continued = TRUE;
                    }

                  /* If this is the accesson on the next line, then save */
                  else if (input_line[5] == '2' && unp_continued == TRUE)
                    {
                      /* Free the tokens */
                      free_tokens(token,ntokens);

                      /* Extract the SWISS-PROT accession number */
                      ntokens
                        = get_tokens(input_line + 14,token,MAX_TOKEN,' ');

                      /* Save the UniProt accession */
                      len = strlen(token[0]);
                      if (len > 2 && len < UNIPROT_LEN )
                        {
                          /* Store code */
                          strcpy(code_value[ntypes],token[0]);
                          type[ntypes] = SWISSPROT_ID;
                          ntypes++;
                        }

                      /* Unset flag */
                      unp_continued = FALSE;
                    }

                  /* If codes not too long (or too short), save them */
                  else if (ntokens > 1 && input_line[5] == ' ')
                    {
                      len = strlen(token[0]);
                      if (len > 2 && len < UNIPROT_LEN)
                        {
                          /* Store code */
                          strcpy(code_value[ntypes],token[0]);
                          type[ntypes] = SWISSPROT_ID;
                          ntypes++;
                        }

                      /* Repeat for the SWISS-PROT code */
                      len = strlen(token[1]);
                      if (len > 2 && len < UNIPROT_LEN )
                        {
                          /* Store code */
                          strcpy(code_value[ntypes],token[1]);
                          type[ntypes] = SWISSPROT_CODE;
                          ntypes++;
                        }

                      /* Unset flag */
                      unp_continued = FALSE;
                    }

                  /* Free the tokens */
                  free_tokens(token,ntokens);
                }
            }

          /* Check for old-style SWISS-PROT code */
          else if (have_sws_data == FALSE &&
                   !strncmp("SWISS-PROT ENTRY NAME",input_line+11,21))
            {
              /* Next line contains the SWISS-PROT code and chain-id */
              next_sws_entry = TRUE;
            }

          /* Check for line containing old-style SWISS-PROT code */
          else if (next_sws_entry == TRUE)
            {
              /* Extract the SWISS-PROT code */
              /* Get the chain-id and the SWISS-PROT id */
              chain = input_line[47];
              strncpy(tmp_string,input_line+14,12);
              tmp_string[12] = '\0';
              string_ptr = strstr(tmp_string," ");
              if (string_ptr != NULL)
                string_ptr[0] = '\0';

              /* If id not too long (or too short), save it */
              len = strlen(tmp_string);
              if (len > 2 && len < 15)
                {
                  /* Store code */
                  strcpy(code_value[ntypes],tmp_string);
                  type[ntypes] = SWISSPROT_CODE;
                  ntypes++;
                }

              /* Unset the flag */
              next_sws_entry = FALSE;
            }

          /* If hit first ATOM or HETATM record can stop reading */
          else if (!strncmp("ATOM  ",input_line,6) ||
                   !strncmp("HETATM",input_line,6))
            done = TRUE;

          /* If have one or more codes to store, then do so */
          if (ntypes > 0)
            {
              /* Get pointer to the first chain record */
              chain_ptr = first_chain_ptr;
              match_chain_ptr = NULL;
              copy_of = '\0';

              /* Loop through all the chain records until find the correct
                 one */
              while (chain_ptr != NULL && match_chain_ptr == NULL)
                {
                  /* If this is the correct record, store the pointer */
                  if (chain_ptr->chain_id == chain ||
                      (chain == ' ' && chain_ptr->chain_id == 'A'))
                    {
                      /* Store pointer */
                      match_chain_ptr = chain_ptr;

                      /* If this is a copy chain, determine which chain it
                         is a copy of */
                      if (chain_ptr->copy_of != '\0')
                        copy_of = chain_ptr->copy_of;
                      else
                        copy_of = chain_ptr->chain_id;
                    }

                  /* Get pointer to the next chain record */
                  chain_ptr = chain_ptr->next_chain_ptr;
                }

              /* Loop over the types */
              for (i = 0; i < ntypes; i++)
                {
                  /* Store the code itself */
                  code_ptr
                    = create_code_record(&first_code_ptr,&last_code_ptr,
                                         type[i],code_value[i],PDB_FILE,
                                         annotation,&new,copy_of);

                  /* If this is a new code, increment count of codes */
                  if (new == TRUE)
                    ncodes[type[i]]++;
                  nitems++;

                  /* Adjust chains count */
                  code_ptr->count--;

                  /* Get pointer to the first chain */
                  chain_ptr = first_chain_ptr;

                  /* Loop through all the chain records to find all those
                     for which this code applies */
                  while (chain_ptr != NULL)
                    {
                      /* If this is the right chain, or a copy of the right
                         one, store pointer to code */
                      if (chain_ptr->chain_id == chain ||
                          chain_ptr->copy_of == chain || chain == '\0' ||
                          (chain == ' ' &&
                           (chain_ptr->chain_id == 'A' ||
                            chain_ptr->copy_of == 'A')))
                        {
                          /* Check whether we already have this code
                             stored for this chain */
                          have_code = FALSE;
                          for (icode = 0; icode < chain_ptr->ncodes; icode++)
                            {
                              /* If have code, then don't want to store
                                 again */
                              if (chain_ptr->code_ptr[icode] == code_ptr)
                                have_code = TRUE;
                            }

                          /* If code not yet stored, then store it now */
                          if (have_code == FALSE)
                            {
                              /* Store for current chain */
                              if (chain_ptr->ncodes < MAX_CODES)
                                {
                                  chain_ptr->code_ptr[chain_ptr->ncodes]
                                    = code_ptr;
                                  chain_ptr->ncodes++;
                                }

                              /* Increment chain count for current code */
                              code_ptr->count++;
                            }
                        }

                      /* Get the next chain */
                      chain_ptr = chain_ptr->next_chain_ptr;
                    }
                }
            }
        }
    }

  /* Close PDB file */
  if (gzipped == TRUE)
    {
#ifdef TEST
#else
      gzclose(gzfile_ptr);
#endif
    }
  else
    fclose(file_ptr);

  /* Return whether E.C. number read in */

  /* Return the first and last code-pointers */
  *fst_code_ptr = first_code_ptr;
  *lst_code_ptr = last_code_ptr;

  /* Return the number of codes retained */
  return(nitems);
}
/***********************************************************************

get_af_ec  -  Get any EC numbers for this protein from the AlphaFold
              analysis file

***********************************************************************/

int get_af_ec(char *pdb_code,struct pdb *pdb_ptr,int nitems,
              struct c_file_data file_name_ptr)
{
  char file_name[FILENAME_LEN], input_line[LINELEN + 1];
  char annotation[TOKEN_LEN], ec_number[20];

  char *string_ptr, *subtoken[MAX_TOKEN], *token[MAX_TOKEN];

  int i, len, ncodes, nsubtokens, ntokens, type;
  int copy_of, done, new;

  struct chain *chain_ptr;
  struct code *code_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  annotation[0] = '\0';
  copy_of = FALSE;
  done = FALSE;
  type = EC_NUMBER;

  /* Initialise pointers */
  code_ptr = NULL;

  /* Form the name of the analysis file */
  strcpy(file_name,file_name_ptr.other_dir[ALPHAFOLD_PDBSUM_DATA_DIR]);
  strcat(file_name,"/AFanalyses.tsv");

  /* Open the file */
  file_ptr = fopen(file_name,"r");
  if (file_ptr ==  NULL)
    {
      printf("*** Warning. Unable to open %s\n",file_name);
      return(nitems);
    }

  /* Read through the file until we find our entry */
  while (fgets(input_line,LINELEN,file_ptr) != NULL && done == FALSE)
    {
      /* Truncate the input line */
      len = string_chomp(input_line);

      /* Extract all the tokens */
      ntokens = get_tokens(input_line,token,MAX_TOKEN,'\t');

      /* See if this is our PDB entry */
      if (ntokens > 12 && !strcmp(token[2],pdb_code))
        {
          /* See if there are any EC numbers */
          if (strcmp(token[11],"-"))
            {
              /* Strip out the individual EC numbers */
              nsubtokens
                = get_tokens(token[11],subtoken,MAX_TOKEN,' ');

              /* Loop over the subtokens to store each EC number */
              for (i = 0; i < nsubtokens; i++)
                {
                  /* If not too long, then store */
                  len = strlen(subtoken[i]);
                  if (len < 20)
                    {
                      /* Get the EC number and remove final comma,
                         if there is one */
                      strcpy(ec_number,subtoken[i]);
                      string_ptr = strchr(ec_number,',');
                      if (string_ptr != NULL)
                        string_ptr[0] = '\0';

                      /* Create a code record */
                      code_ptr
                        = create_code_record(&(pdb_ptr->first_code_ptr),
                                             &(pdb_ptr->last_code_ptr),
                                             type,ec_number,PDB_FILE,
                                             annotation,&new,copy_of);

                      /* If this is a new code, increment count of codes */
                      if (new == TRUE)
                        pdb_ptr->ncodes[type]++;

                      /* Get the chain record */
                      chain_ptr = pdb_ptr->first_chain_ptr;

                      /* Save pointer to code record */
                      ncodes = chain_ptr->ncodes;
                      chain_ptr->code_ptr[ncodes] = code_ptr;
                      chain_ptr->ncodes++;
                    }
                }

              /* Free the subtokens */
              free_tokens(subtoken,nsubtokens);

              /* Increment item count */
              nitems++;
            }

          /* Set flag that we're done */
          done = TRUE;
        }

      /* Free the tokens */
      free_tokens(token,ntokens);
    }

  /* Close the file */
  fclose(file_ptr);

  /* Return the number of codes retained */
  return(nitems);
}
/***********************************************************************

get_msd_data  -  Read through the pdb_chain_sp_ec file to pick up the
                 UniProt codes and EC numbers associated with each
                 PDB file

***********************************************************************/

int get_msd_data(struct pdb **pdb_index_ptr,int npdb,
                 struct c_file_data file_name_ptr,
                 struct parameters params)
{
  char pdb_code[5], file_name[FILENAME_LEN], pdbsum_dir[FILENAME_LEN];

  int dir_type, ipdb, nitems;

  struct code *first_code_ptr, *last_code_ptr;
  struct pdb *pdb_ptr;

  FILE *auditfile_ptr;

  /* Initialise variables */
  nitems = 0;

  /* Initialise pointers */
  first_code_ptr = NULL;
  last_code_ptr = NULL;

  /* Loop through the PDB entries */
  for (ipdb = 0; ipdb < npdb; ipdb++)
    {
      /* Get this PDB record */
      pdb_ptr = pdb_index_ptr[ipdb];

      /* Get the PDB code */
      strcpy(pdb_code,pdb_ptr->pdb_code);
      if (ipdb % 1000 == 0)
	{
	  printf("%d. Getting PDBsum data for %s\n",ipdb,pdb_code);
	  fflush(stdout);
	}

      /* Form the name of the PDBsum directory */
      if (params.pdb_type == PDBSUM1)
        {
          dir_type = PDBSUM1;
          strcpy(pdbsum_dir,"../newsec/");
        }
      else
        {
          dir_type
            = find_pdbsum_dir(pdb_code,file_name_ptr.pdbsum_template,
                              params.pdb_type,pdbsum_dir);
        }
      if (dir_type != NOT_FOUND)
        {
          /* Read in all the chain-equivalences from the chains.dat file */
          pdb_ptr->nchains
            = read_chain_equivs(pdb_code,pdbsum_dir,
                                &(pdb_ptr->first_chain_ptr));
        }
      else
        printf("\n*** Warning. Unable to find PDBsum directory "
               "for PDB code [%s]\n",pdb_code);
    }

  /* Read through the MSD data and match up to our PDB codes */
  if (params.pdb_type == STANDARD_PDB)
    nitems
      = read_msd_data(pdb_index_ptr,npdb,
                      file_name_ptr.other_dir[PDBSUM_REAL_DATA_DIR]);
  
  /* Loop through all the PDB entries again to pick up data from PDB
     file from any not in the MSD */
  for (ipdb = 0; ipdb < npdb; ipdb++)
    {
      /* Get this PDB record */
      pdb_ptr = pdb_index_ptr[ipdb];

      /* Get the PDB code */
      strcpy(pdb_code,pdb_ptr->pdb_code);
      if (ipdb % 1000 == 0)
	{
	  printf("%d. Getting PDB data for %s\n",ipdb,pdb_code);
	  fflush(stdout);
	}

      /* If no data from the MSD database, then read the PDB file to see
         if the data can be extracted from it */
      if (pdb_ptr->have_ec == FALSE)
        {
          /* Get the data from the PDB file */
          nitems = get_pdb_data(pdb_code,pdb_ptr->first_chain_ptr,
                                &(pdb_ptr->first_code_ptr),
                                &(pdb_ptr->last_code_ptr),pdb_ptr->ncodes,
                                nitems,file_name_ptr,params);

          /* If this is an AlphaFold model, see if we have EC numbers
             for it in the analysis file */
          if (params.pdb_type == ALPHA_FOLD_PDB)
            nitems = get_af_ec(pdb_code,pdb_ptr,nitems,file_name_ptr);
        }
    }

  /* Return number of items picked up */
  return(nitems);
}
/***********************************************************************

get_go_data  -  Get the PDB->GO annotations from the gzipped GO file

***********************************************************************/

int get_go_data(struct pdb **pdb_index_ptr,int npdb,
                struct c_file_data file_name_ptr)
{
  char file_name[FILENAME_LEN + 1], input_line[LINELEN + 1];
  char zfile_name[FILENAME_LEN + 1], go_names_file[FILENAME_LEN + 1];
  char command[FILENAME_LEN + 1];
  char brcode[5], chain, copy_of, pdb_code[5], prefix;
  char code[TOKEN_LEN];
  char annotation[TOKEN_LEN];

  char *token[MAX_TOKEN];

  int i, icode, ipdb, len, ngo, ngo_types[3], nlines, ntokens, source;
  int go_type, type;
  int err_no, have_code, new;

  struct code *code_ptr, *first_code_ptr, *last_code_ptr;
  struct chain *chain_ptr, *match_chain_ptr;
  struct pdb *pdb_ptr;

#ifdef TEST
  FILE *file_ptr;
#else
  FILE *file_ptr;
  gzFile gzfile_ptr;
#endif

  /* Initialise counts */
  ngo = 0;
  for (i = 0; i < 3; i++)
    ngo_types[i] = 0;
  nlines = 0;
  chain = '\0';
  match_chain_ptr = NULL;
  type = UNKNOWN;
  source = GO_FILE;

  /* Initialise pointers */
  code_ptr = NULL;

  /* Initialise file names */
  strcpy(go_names_file,file_name_ptr.other_dir[ZPDB_TO_GO_FILE]);
#ifdef TEST
  strcpy(file_name,"gene_association.goa_pdb");

  /* Open the file, if it exists */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    {
      printf("\n*** Warning. Unable to open file [%s]\n",file_name);
      return(ngo);
    }

  /* Loop while reading in records from input file */
  while (fgets(input_line,LINELEN,file_ptr)!= NULL)
#else
    /* Form name of input file */
    strcpy(zfile_name,file_name_ptr.other_dir[ZPDB_TO_GO_FILE]);

  /* Open the file, if it exists */
  gzfile_ptr = gzopen(zfile_name,"rb");
  if (gzfile_ptr ==  NULL)
    {
      /* Check for insufficiency of memory */
      if(err_no)
        {
          printf("*** ERROR. Insufficient memory opening file [%s]\n",
                 zfile_name);
          return(ngo);
        }

      /* Otherwise, error must be due to missing file */
      printf("\n*** Warning. Unable to open file [%s]\n",zfile_name);
      return(ngo);
    }

  /* Loop while reading in records from the pdb2go file */
  while (gzgets(gzfile_ptr,input_line,LINELEN)!= NULL)
#endif
    {
      /* Truncate the string to strip off the newline */
      len = string_truncate(input_line,LINELEN);

      /* Increment count of lines read in */
      nlines++;

      /* Extract all the tokens on this line */
      ntokens = get_tokens(input_line,token,MAX_TOKEN,'\t');

      /* Get the GO prefix (P,F or C) and store */
      type = UNKNOWN;
      if (ntokens > 8 && token[0][0] != '!')
        {
          prefix = token[8][0];
          if (prefix == 'C')
            type = GO_COMPONENT;
          else if (prefix == 'F')
            type = GO_FUNCTION;
          else if (prefix == 'P')
            type = GO_PROCESS;
        }

      /* If have sufficient tokens get the corresponding PDB code */
      if (ntokens > 11 && type != UNKNOWN)
        {
          /* Get the PDB code and convert to lower case */
          strncpy(brcode,token[1],4);
          brcode[4] = '\0';
          to_lower(brcode,4);

          /* See if this corresponds to one of our PDB codes */
          pdb_ptr = binary_pdb_search(pdb_index_ptr,npdb,brcode);

          /* If found, then extract and store GO annotation data */
          if (pdb_ptr != NULL)
            {
              /* Get the chain identifier */
              chain = token[1][5];
              if (chain == '@')
                chain = ' ';

              /* Get pointer to the first chain record */
              chain_ptr = pdb_ptr->first_chain_ptr;
              match_chain_ptr = NULL;

              /* Loop through all the chain records until find the
                 correct one */
              while (chain_ptr != NULL && match_chain_ptr == NULL)
                {
                  /* If this is the correct record, store the pointer */
                  if (chain_ptr->chain_id == chain)
                    match_chain_ptr = chain_ptr;

                  /* Get pointer to the next chain record */
                  chain_ptr = chain_ptr->next_chain_ptr;
                }

              /* If have a match, store the GO code */
              if (match_chain_ptr != NULL)
                {
                  /* Get the GO code and store */
                  strcpy(code,token[4]+3);

                  /* Get the source of the annotation */
                  strcpy(annotation,token[6]);

                  /* If this is a copy chain, determine which chain
                     it is a copy of */
                  if (match_chain_ptr->copy_of != '\0')
                    copy_of = match_chain_ptr->copy_of;
                  else
                    copy_of = match_chain_ptr->chain_id;

                  /* Store the code itself */
                  code_ptr
                    = create_code_record(&(pdb_ptr->first_code_ptr),
                                         &(pdb_ptr->last_code_ptr),
                                         type,code,source,annotation,
                                         &new,copy_of);
                  ngo++;
                  go_type = type - GO_COMPONENT;
                  ngo_types[go_type]++;

                  /* If this is a new code, increment count of codes */
                  if (new == TRUE)
                    pdb_ptr->ncodes[type]++;

                  /* Store for current chain */
                  if (match_chain_ptr->ncodes < MAX_CODES)
                    {
                      /* Determine whether this chain already has this
                         code */
                      have_code = FALSE;

                      /* Loop over the stored codes */
                      for (icode = 0; icode < match_chain_ptr->ncodes;
                           icode++)
                        {
                          /* If have code, then set flag */
                          if (match_chain_ptr->code_ptr[icode]
                              == code_ptr)
                            have_code = TRUE;
                        }

                      /* If this chain doesn't already have this code,
                         then store */
                      if (have_code == FALSE)
                        {
                          match_chain_ptr->code_ptr[match_chain_ptr->ncodes]
                            = code_ptr;
                          match_chain_ptr->ncodes++;
                        }
                    }
                }
            }
        }

      /* Free the tokens */
      free_tokens(token,ntokens);
    }

  /* Close the input file */
#ifdef TEST
  fclose(file_ptr);
#else
  gzclose(gzfile_ptr);
#endif

  /* Return the number of GO terms encountered */
  return(ngo);
}
/***********************************************************************

create_go_record  -  Create a new go record

***********************************************************************/

struct go *create_go_record(struct go **fst_go_ptr,
                            struct go **lst_go_ptr,char *code,
                            int type,char *name)
{
  struct go *go_ptr;
  struct go *first_go_ptr, *last_go_ptr;

  /* Initialise go pointers */
  last_go_ptr = *lst_go_ptr;
  first_go_ptr = *fst_go_ptr;

  /* Create a new go entry */
  go_ptr = (struct go*)malloc(sizeof(struct go));
  if (go_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for");
      printf(" struct go\n");
      printf("***        Program seqdata terminated with error.\n");
      exit (1);
    }

  /* Store the details for this entry */
  strcpy(go_ptr->go_id,code);
  go_ptr->type = type;
  store_string(&(go_ptr->name),name);

  /* Initialise remaining fields */
  go_ptr->ave_level = 0.0;
  go_ptr->level = 0;

  /* Initialise pointer to next go record */
  go_ptr->next_go_ptr = NULL;

  /* If this is the first entry stored, then update pointer to
     head of linked list */
  if (first_go_ptr == NULL)
    first_go_ptr = go_ptr;

  /* Otherwise, make previous entry point to the current
     one */
  else
    last_go_ptr->next_go_ptr = go_ptr;
                      
  /* Save the current go's pointer */
  last_go_ptr = go_ptr;

  /* Return current pointers */
  *fst_go_ptr = first_go_ptr;
  *lst_go_ptr = last_go_ptr;
  return(go_ptr);
}
/***********************************************************************

get_go_names  -  Read through the GO terms and ids file to pick up the
                 description of each GO number we have

***********************************************************************/

int get_go_names(struct go **fst_go_ptr,struct c_file_data file_name_ptr)
{
  char file_name[FILENAME_LEN], input_line[LINELEN + 1];
  char code[TOKEN_LEN];

  char *token[MAX_TOKEN];

  int len, nnames, ntokens, type;
  int found;

  struct code *code_ptr;
  struct go *go_ptr, *first_go_ptr, *last_go_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  nnames = 0;
  strcpy(file_name,file_name_ptr.other_dir[GO_NAMES_FILE]);

  /* Initialise pointers */
  first_go_ptr = *fst_go_ptr = NULL;
  last_go_ptr = NULL;

  /* Open file */
  printf("Reading GO names file (%s) ...\n",file_name);
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    {
      printf("\n*** Warning. Unable to open GO names file [%s]\n",
             file_name);
      return(nnames);
    }

  /* Loop while reading in records from the input file */
  while (fgets(input_line,LINELEN,file_ptr)!= NULL)
    {
      /* Truncate the string to strip off the newline */
      len = string_truncate(input_line,LINELEN);

      /* Extract all the tokens on this line */
      ntokens = get_tokens(input_line,token,MAX_TOKEN,'\t');

      /* Check that this is a valid line */
      if (ntokens > 2 && !strncmp(token[0],"GO:",3))
        {
          /* Extract the GO number */
          strcpy(code,token[0]+3);

          /* Get the ontology */
          if (token[2][0] == 'C')
            type = GO_COMPONENT;
          else if (token[2][0] == 'F')
            type = GO_FUNCTION;
          else if (token[2][0] == 'P')
            type = GO_PROCESS;
          else
            type = UNKNOWN;

          /* Create a GO record to store this description */
          if (type != UNKNOWN)
            {
              go_ptr = create_go_record(&first_go_ptr,&last_go_ptr,
                                        code,type,token[1]);

              /* Increment count of names read in */
              nnames++;
            }
        }

      /* Free the tokens */
      free_tokens(token,ntokens);
    }

  /* Close the file */
  fclose(file_ptr);

  /* Return first of the GO records */
  *fst_go_ptr = first_go_ptr;

  /* Return number of names read in */
  return(nnames);
}
/***********************************************************************

go_sort  -  Shell-sort which sorts the UniProt gos into ascending
              order

***********************************************************************/

void go_sort(struct go **go_index_ptr,int nindex)
{
  char go_id1[GO_LENGTH + 1], go_id2[GO_LENGTH + 1];

  int endloop, i, indi, indl, j, k, l, lognb2, m, n, nn;
  float aln2i = (float) 1.4426950, tiny = (float) 1.e-5;
  float calc;
  double dble, dble1;

  struct go *go1_ptr, *go2_ptr, *swap_go_ptr;

  /*--Initialise values */
  n = nindex;
  dble = n;
  dble1 = log(dble);
  calc = (float) dble1;
  calc = calc * aln2i + tiny;
  lognb2 = (int) calc;

  /* Perform sort */
  m = n;
  for (nn = 1; nn <= lognb2; nn++)
    {
      m = m / 2;
      k = n - m;
      for (j = 1; j <= k; j++)
        {
          i = j;
          endloop = FALSE;
          while (endloop == FALSE)
            {
              l = i + m;
              indi = i - 1;
              indl = l - 1;

              /* Get the two go records */
              go1_ptr = go_index_ptr[indi];
              go2_ptr = go_index_ptr[indl];

              /* Form full go names on which sort is based */
              strcpy(go_id1,go1_ptr->go_id);
              strcpy(go_id2,go2_ptr->go_id);

              /* If in wrong order, the swap sort-pointers */
              if (strcmp(go_id1,go_id2) > 0)
                {
                  swap_go_ptr = go_index_ptr[indi];
                  go_index_ptr[indi] = go_index_ptr[indl];
                  go_index_ptr[indl] = swap_go_ptr;

                  i = i - m;
                  if (i < 1)
                    endloop = TRUE;
                }
              else
                endloop = TRUE;
            }
        }
    }
}
/***********************************************************************

create_go_index  -  Create index of GO records for use in binary
                     search

***********************************************************************/

void create_go_index(struct go *first_go_ptr,
                      struct go ***go_index_ptr,int ngo)
{
  int nentry;

  struct go *go_ptr;

  struct go **go_ind_ptr;

  /* Initialise variables */

  /* Create the index array as one large array */
  go_ind_ptr
    = (struct go **) malloc(ngo * sizeof(struct go *));

  /* Check that sufficient memory was allocated */
  if (go_ind_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for index array\n");
      exit(1);
    }

  /* Get pointer to the first go entry */
  go_ptr = first_go_ptr;
  nentry = 0;

  /* Loop through all entries to place each one in the sorted index
     array */
  while (go_ptr != NULL)
    {
      /* Store current pointer at end of the index array */
      if (nentry < ngo)
        go_ind_ptr[nentry] = go_ptr;

      /* Increment count of stored entries */
      nentry++;

      /* Get pointer to the next NIST entry in the linked list */
      go_ptr = go_ptr->next_go_ptr;
    }

  /* Sort the code pointers into order of UniProt code */
  if (ngo > 1)
    go_sort(go_ind_ptr,ngo);

  /* Return the array pointer */
  *go_index_ptr = go_ind_ptr;
}
/***********************************************************************

binary_go_search  -  Search for the given GO record

***********************************************************************/

struct go *binary_go_search(struct go **go_index_ptr,int ngo,
                            char *go_id)
{
  int bottom_point, ipoint, new_point, top_point;
  int found;

  struct go *go_ptr, *found_go_ptr, *prev_go_ptr;

  /* Initialise variables */
  found = FALSE;
  found_go_ptr = NULL;

  /* Make starting point half-way into the index */
  ipoint = ngo / 2;
  bottom_point = 0;
  top_point = ngo - 1;

  /* Loop until required index position found */
  while (found == FALSE)
    {
      /* Get the associated GO pointer */
      go_ptr = go_index_ptr[ipoint];

      /* Check if index value is equal to the value we want */
      if (!strcmp(go_id,go_ptr->go_id))
        {
          /* Set flag to indicate that match found and store pointer */
          found = TRUE;
          found_go_ptr = go_ptr;
        }

      /* If index value is higher than the one we want, then need
         to search lower */
      else if (strcmp(go_id,go_ptr->go_id) < 0)
        {
          /* If this is the first go in the sorted list, then can't
             go any lower */
          if (ipoint == 0)
            {
              /* Set flag that we are done */
              found = TRUE;
            }

          /* Otherwise, check whether the previous entry's name
             is below the one we're searching for */
          else
            {
              /* Get the pointer to the previous entry */
              prev_go_ptr = go_index_ptr[ipoint - 1];

              /* If value is below the search value, we've found the
                 position in the index we're after, but don't have an exact
                 match*/
              if (strcmp(prev_go_ptr->go_id,go_id) < 0)
                found = TRUE;

              /* Otherwise, need to search earlier in the list */
              else
                {
                  /* Save the current position as the uppermost position
                     to search on */
                  top_point = ipoint;

                  /* Set pointer to halfway between the last search
                     low-search and here */
                  ipoint = (ipoint + bottom_point) / 2;
                }
            }
        }

      /* Otherwise, have to look further up the list */
      else
        {
          /* Save the current position as the lowest position to search
             on */
          bottom_point = ipoint;

          /* If have reached the top of the list, then our value is
             higher than any value in the list */
          if (ipoint >= ngo - 1)
            {
              /* Set flag as done */
              found = TRUE;
            }

          /* Otherwise, can search further up the list */
          else
            {
              /* Set pointer to halfway between the last search high-search
                 and here */
              new_point = (ipoint + top_point) / 2;

              /* Adjust if going back to the same point */
              if (new_point == ipoint)
                new_point++;
              ipoint = new_point;
            }
        }
    }

  /* Return the GO pointer */
  return(found_go_ptr);
}
/***********************************************************************

match_go  -  Match the GO records to those associated with each PDB
             entry

***********************************************************************/

void match_go(struct pdb **pdb_index_ptr,int npdb,
              struct go **go_index_ptr,int ngo)
{
  int icode, ipdb, type;

  struct chain *chain_ptr;
  struct code *code_ptr;
  struct pdb *pdb_ptr;

  /* Loop through the PDB entries */
  for (ipdb = 0; ipdb < npdb; ipdb++)
    {
      /* Get this PDB record */
      pdb_ptr = pdb_index_ptr[ipdb];

      /* Get pointer to the first chain record */
      chain_ptr = pdb_ptr->first_chain_ptr;

      /* Loop through all the chain records to write out SWISS-PROT codes
         and E.C. numbers */
      while (chain_ptr != NULL)
        {
          /* If this is a unique chain, write its details out */
          if (chain_ptr->copy_of == '\0')
            {
              /* Loop over the stored codes */
              for (icode = 0; icode < chain_ptr->ncodes; icode++)
                {
                  /* Get this code */
                  code_ptr = chain_ptr->code_ptr[icode];

                  /* Get this code's type */
                  type = code_ptr->type;

                  /* If this is a GO code, find its corresponding
                     name record */
                  if (type == GO_COMPONENT || type == GO_FUNCTION ||
                      type == GO_PROCESS)
                    {
                      /* Find this GO record */
                      code_ptr->go_ptr
                        = binary_go_search(go_index_ptr,ngo,
                                           code_ptr->value);
                    }
                }
            }

          /* Get pointer to the next chain record */
          chain_ptr = chain_ptr->next_chain_ptr;
        }
    }
}
/***********************************************************************

read_ontology_file  -  Read through the given .ontology file to pick up
                       all the names and levels for each annotation of
                       this GO type

***********************************************************************/

void read_ontology_file(char *file_name,struct go **go_index_ptr,
                        int ngo,int type)
{
  char input_line[LINELEN + 1];
  char description[LINELEN + 1], go_number[GO_LENGTH + 1];

  char *string_ptr, *end_ptr;

  int go_type, len, level, nlevel;

  struct go *go_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  printf("\n");
  printf("Reading .ontology file (%s) ...\n",file_name);
  printf("\n");
  go_type = type + GO_COMPONENT;
  level = 0;
  nlevel = -1;

  /* Open file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    {
      printf("\n*** ERROR. Unable to open .ontology file [%s]\n",file_name);
      return;
    }

  /* Loop while reading in records from the input file */
  while (fgets(input_line,LINELEN,file_ptr)!= NULL)
    {
      /* Truncate the string to strip off the newline */
      len = string_truncate(input_line,LINELEN);

      /* Check whether this is a definition line */
      string_ptr = strstr(input_line,"  %");
      if (string_ptr == NULL)
        string_ptr = strstr(input_line,"  <");

      /* If line contains a definition, determine what level it lies at */
      if (string_ptr != NULL)
        {
          /* Get definition's level */
          level = string_ptr - input_line;

          /* Get the end of the name, plus the GO id */
          end_ptr = strstr(string_ptr + 3," ;");
          if (end_ptr != NULL && !strncmp(end_ptr + 3,"GO:",3))
            {
              /* Extract the GO number */
              strncpy(go_number,end_ptr + 6,7);
              go_number[7] = '\0';
            }

          /* If have a terminator, truncate string */
          if (end_ptr != NULL)
            end_ptr[0] = '\0';

          /* Get the length of the name */
          len = strlen(string_ptr + 3);

          /* Store the name */
          strcpy(description,string_ptr + 3);

          /* Find this GO record */
          go_ptr = binary_go_search(go_index_ptr,ngo,go_number);

          /* If found, store the details */
          if (go_ptr != NULL)
            {
              /* Store level */
              go_ptr->ave_level = go_ptr->ave_level + level;
              go_ptr->level++;
            }
        }
    }

  /* Close the file */
  fclose(file_ptr);
}
/***********************************************************************

sort_go_terms  -  Sort current set of GO terms according to their
                  average number and renumber the levels

***********************************************************************/

void sort_go_terms(struct code **fst_code_ptr,int total_codes)
{
  int nswap, maxswap;
  int swapped;

  struct code *first_code_ptr, *code_ptr, *next_code_ptr, *prev_code_ptr;

  /* Initialise variables */
  nswap = 0;
  maxswap = total_codes * total_codes;

  /* Store starting code */
  first_code_ptr = *fst_code_ptr;

  /* Initialise variables */
  prev_code_ptr = NULL;
  swapped = TRUE;

  /* Loop until chains in required order */
  while (swapped == TRUE && nswap < maxswap)
    {
      /* Initialise swapped flag */
      swapped = FALSE;
      nswap++;

      /* Get pointer to the first code record */
      code_ptr = first_code_ptr;
      prev_code_ptr = NULL;

      /* Loop through the codes to swap any that are in the wrong order */
      while (code_ptr != NULL)
        {
          /* Get the next code record */
          next_code_ptr = code_ptr->next_code_ptr;

          /* If have a next record of the right type, check the
             order */
          if (next_code_ptr != NULL)
            {
              /* If codes in wrong order, need to swap */
              if (code_ptr->type > next_code_ptr->type ||
                  (code_ptr->type == next_code_ptr->type &&
                   code_ptr->go_ptr != NULL &&
                   next_code_ptr->go_ptr != NULL &&
                   code_ptr->go_ptr->ave_level >
                   next_code_ptr->go_ptr->ave_level))
                {
                  /* Perform the swap */

                  /* If have a previous record, get it to point to
                     the next entry */
                  if (prev_code_ptr != NULL)
                    {
                      /* Change the pointer */
                      prev_code_ptr->next_code_ptr = next_code_ptr;
                    }

                  /* Otherwise, if have no previous record, it must
                     be the first one */
                  else
                    first_code_ptr = next_code_ptr;

                  /* Transfer next record's next pointer to current
                     record */
                  code_ptr->next_code_ptr
                    = next_code_ptr->next_code_ptr;

                  /* Get the next record to point back to this one */
                  next_code_ptr->next_code_ptr = code_ptr;
                  next_code_ptr = code_ptr->next_code_ptr;
                  
                  /* Set flag that swap made */
                  swapped = TRUE;
                }
            }

          /* Store current code pointer */
          prev_code_ptr = code_ptr;

          /* Get pointer to the next code record */
          code_ptr = next_code_ptr;
        }
    }

  /* Return starting code pointer */
  *fst_code_ptr = first_code_ptr;
}
/***********************************************************************

renumber_levels  -  Renumber the GO levels using the average levels as a
                    guide

***********************************************************************/

void renumber_levels(struct code *first_code_ptr)
{
  int ilevel, level, next_level, type;

  struct code *code_ptr;

  /* Initialise variables */
  level = 0;
  next_level = 0;
  type = 0;

  /* Get pointer to first entry */
  code_ptr = first_code_ptr;

  /* Loop through all the codes to assign a level number */
  while (code_ptr != NULL)
    {
      /* Consider only if type is a GO annotation */
      if (code_ptr->type > 2 && code_ptr->go_ptr != NULL)
        {
          /* If type has changed, reset levels */
          if (code_ptr->type != type)
            {
              /* Re-initialise levels */
              level = 0;
              next_level = 0;

              /* Store the type */
              type = code_ptr->type;
            }

          /* Convert average level to an integer */
          ilevel = (int) code_ptr->go_ptr->ave_level;

          /* If don't know which level to assign this to, increment
             level and store mapping */
          if (ilevel != level)
            {
              /* Increment level and store mapping */
              next_level++;
              level = ilevel;
            }

          /* Store the current level */
          code_ptr->go_ptr->level = next_level;
        }

      /* Get next code entry */
      code_ptr = code_ptr->next_code_ptr;
    }
}
/***********************************************************************

get_go_levels  -  Read through the three ontology files to pick up the
                  go names and levels

***********************************************************************/

void get_go_levels(struct go **go_index_ptr,int ngo,
                   struct pdb **pdb_index_ptr,int npdb,
                   struct c_file_data file_name_ptr)
{
  char file_name[FILENAME_LEN + 1];

  static char *ontology[3] = { "component", "function", "process" };

  int i, igo, ipdb, total_codes, type;

  struct go *go_ptr;
  struct pdb *pdb_ptr;

  /* Loop over the GO types */
  for (type = 0; type < 3; type++)
    {
      /* Form name of this .ontology file */
      strcpy(file_name,file_name_ptr.other_dir[GO_ONTOLOGY_DIR]);
      strcat(file_name,DIR_SEP);
      strcat(file_name,ontology[type]);
      strcat(file_name,".ontology");

      /* Read through the file to pick up all the names and levels
         for each annotation of this GO type */
      read_ontology_file(file_name,go_index_ptr,ngo,type);
    }

  /* Loop through all the GO terms to calculate the average GO level
     for each term  */
  for (igo = 0; igo < ngo; igo++)
    {
      /* Get this GO record */
      go_ptr = go_index_ptr[igo];

      /* Calculate average level */
      if (go_ptr->level > 0)
        go_ptr->ave_level = go_ptr->ave_level / go_ptr->level;
    }

  /* Loop over all PDB entries to sort their GO terms and number their
     levels */
  for (ipdb = 0; ipdb < npdb; ipdb++)
    {
      /* Get this PDB record */
      pdb_ptr = pdb_index_ptr[ipdb];

      /* Count the total number of codes store for this PDB entry */
      total_codes = 0;
      for (i = 0; i < NTYPES; i++)
        total_codes = total_codes + pdb_ptr->ncodes[i];

      /* Sort the GO terms according to their average level and
         renumber the levels */
      if (total_codes > 1)
        sort_go_terms(&(pdb_ptr->first_code_ptr),total_codes);

      /* Renumber the levels */
      renumber_levels(pdb_ptr->first_code_ptr);
    }
}
/***********************************************************************

write_seqinfo  -  Write the output seqinfo.dat files

***********************************************************************/

void write_seqinfo(struct pdb **pdb_index_ptr,int npdb,
                   struct c_file_data file_name_ptr,
                   struct parameters params)
{
  char file_name[LINELEN + 1], pdb_code[5];
  char pdbsum_dir[FILENAME_LEN];
  char desc[LINELEN + 1], value[40];

  static char *code_desc[NTYPES]
    = { "SPROT_ID", "SPROT_CODE", "EC_NO", "GO_ID_COMP", "GO_ID_FUNC",
        "GO_ID_PROC" };

  static char *code_type[NTYPES]
    = { "SWISSPROT_ID", "SWISSPROT_CODE", "EC_NUMBER", "GO_NUMBER_COMP",
        "GO_NUMBER_FUNC", "GO_NUMBER_PROC" };

  static char *description[NTYPES]
    = { " ", " ", " ", "COMP", "FUNC", "PROC" };

  static char *code[NTYPES]
    = { "SWS_ID:  ", "SWS_CODE:", "EC_NO:   ", "GO_ID_C: ", "GO_ID_F: ",
        "GO_ID_P: " };

  static char *source[3] = { "MSD", "PDB", "GO " };

  int chain_no, count, iannot, ichain, icode, ipdb, iunique, ncopy;
  int dir_type, level, ntypes[NTYPES], type;
  int first, have_go, any_go, ok;

  struct chain *chain_ptr, *other_chain_ptr;
  struct code *code_ptr;
  struct go *go_ptr;
  struct pdb *pdb_ptr;

  FILE *file_ptr, *outfile_ptr, *seqdata_ptr;

  /* Open summary data file */
  strcpy(file_name,"seqdata.dat");

  /* Open the output file */
  seqdata_ptr = open_output_file(file_name,TRUE);

  /* Form name of output pdb_EC file */
  strcpy(file_name,"pdb_EC");

  /* Open the output file */
  outfile_ptr = open_output_file(file_name,TRUE);

  /* Loop through the PDB entries */
  for (ipdb = 0; ipdb < npdb; ipdb++)
    {
      /* Initialise variables */
      ichain = 0;
      any_go = FALSE;

      /* Get this PDB record */
      pdb_ptr = pdb_index_ptr[ipdb];

      /* Get the PDB code */
      strcpy(pdb_code,pdb_ptr->pdb_code);

      /* Form the name of the PDBsum directory */
      if (params.pdb_type == PDBSUM1)
        {
          strcpy(pdbsum_dir,"./");
          dir_type = PDBSUM1;
        }
      else
        dir_type
          = find_pdbsum_dir(pdb_code,file_name_ptr.pdbsum_template,
                            params.pdb_type,pdbsum_dir);
      if (dir_type != NOT_FOUND)
        {
          /* Form output file name */
          strcpy(file_name,pdbsum_dir);
          strcat(file_name,"seqinfo.dat");

          /* Open the output file */
          file_ptr = open_output_file(file_name,FALSE);

          /* If file opened, set flag */
          if (file_ptr != NULL)
            ok = TRUE;
          else
            {
              /* Set flag */
              ok = FALSE;

              /* Show warning */
              printf("*** Warning. Unable to open output file: %s\n",
                     file_name);
            }

          /* Get pointer to the first chain record */
          chain_ptr = pdb_ptr->first_chain_ptr;

          /* Loop through all the chain records to write out
             SWISS-PROT codes and E.C. numbers */
          while (chain_ptr != NULL)
            {
              /* If this is a unique chain, write its details out */
              if (chain_ptr->copy_of == '\0')
                {
                  /* Initialise counts */
                  for (type = 0; type < NTYPES; type++)
                    ntypes[type] = 0;

                  /* Increment unique chain counter */
                  ichain++;

                  /* Loop over the stored codes */
                  for (icode = 0; icode < chain_ptr->ncodes; icode++)
                    {
                      /* Get this code's type */
                      type = chain_ptr->code_ptr[icode]->type;

                      /* Increment types count */
                      ntypes[type]++;

                      /* Write out the code */
                      if (ok == TRUE) 
                        fprintf(file_ptr,"%s[%d][%d] %s\n",
                                code_type[type],ichain,ntypes[type],
                                chain_ptr->code_ptr[icode]->value);

                      /* Write out to the summary file */
                      code_ptr = chain_ptr->code_ptr[icode];
                      strcpy(value,code_ptr->value);
                      strcat(value,"            ");
                      value[12] = '\0';
                      if (ok == TRUE)
                        fprintf(seqdata_ptr,"%s%c %s %s %s",pdb_code,
                                chain_ptr->chain_id,code[type],value,
                                source[code_ptr->source]);

                      /* Get pointer to the first chain record */
                      other_chain_ptr = pdb_ptr->first_chain_ptr;
                      ncopy = 0;

                      /* Loop through all the chains to find copies of
                         the current one */
                      while (other_chain_ptr != NULL)
                        {
                          /* If this is a copy of the current chain,
                             write out the chain identifier */
                          if (other_chain_ptr->copy_of
                              == chain_ptr->chain_id)
                            {
                              /* If this is the first copy, write out
                                 space */
                              if (ncopy == 0 && ok == TRUE)
                                fprintf(seqdata_ptr,"   ");

                              /* Otherwise, write out a comma */
                              else if (ok == TRUE)
                                fprintf(seqdata_ptr,",");

                              /* Write out the copy chain */
                              if (ok == TRUE) 
                                fprintf(seqdata_ptr,"%c",
                                        other_chain_ptr->chain_id);

                              /* Increment count of copies */
                              ncopy++;
                            }

                          /* Get pointer to the next chain record */
                          other_chain_ptr = other_chain_ptr->next_chain_ptr;
                        }

                      /* Close off the current seqdata.dat line */
                      if (ok == TRUE) 
                        fprintf(seqdata_ptr,"\n");
                    }

                  /* Initialise GO flag */
                  have_go = FALSE;

                  /* Write out the counts of each code type for this chain */
                  for (type = 0; type < NTYPES; type++)
                    {
                      /* Write out if non-zero */
                      if (ntypes[type] > 0)
                        {
                          /* Write out the count */
                          if (ok == TRUE) 
                            fprintf(file_ptr,"N%s[%d] %d\n",code_type[type],
                                    ichain,ntypes[type]);

                          /* If this is a GO annotation, set flag */
                          if (type > 2)
                            have_go = TRUE;
                        }
                    }

                  /* If have any GO annotation, write out flag */
                  if (have_go == TRUE)
                    {
                      /* Write out flag indicating that this chain has GO
                         annotation */
                      if (ok == TRUE) 
                        fprintf(file_ptr,"GO_ANNOT[%d] TRUE\n",ichain);

                      /* Set overall flag */
                      any_go = TRUE;
                    }
                }

              /* Write out any E.C. numbers to the pdb_EC file */
              first = TRUE;

              /* Loop over the stored codes */
              for (icode = 0; icode < chain_ptr->ncodes; icode++)
                {
                  /* Get this code's type */
                  type = chain_ptr->code_ptr[icode]->type;

                  /* If it's an enzyme number, write out */
                  if (type == EC_NUMBER)
                    {
                      /* If this is the first one for this chain, write
                         out PDB code and chain identifier */
                      if (first == TRUE)
                        {
                          /* Write out PDB code and chain identifier */
                          if (ok == TRUE) 
                            fprintf(outfile_ptr,"%s%c ",pdb_code,
                                    chain_ptr->chain_id);

                          /* Reset flag */
                          first = FALSE;
                        }

                      /* Write out the E.C number */
                      if (ok == TRUE) 
                        fprintf(outfile_ptr," %s",
                                chain_ptr->code_ptr[icode]->value);
                    }
                }

              /* If any E.C. numbers written to the pdb_EC file, close off
                 this line */
              if (ok == TRUE && first == FALSE)
                fprintf(outfile_ptr,"\n");

              /* Get pointer to the next chain record */
              chain_ptr = chain_ptr->next_chain_ptr;
            }

          /* Loop over the different code types */
          for (type = 0; type < NTYPES; type++)
            {
              /* Initialise counter */
              count = 0;

              /* If have any of this code type, write out */
              if (pdb_ptr->ncodes[type] > 0)
                {
                  /* Write out the number of codes of this type */
                  if (ok == TRUE) 
                    fprintf(file_ptr,"N%s %d\n",code_desc[type],
                            pdb_ptr->ncodes[type]);

                  /* Get pointer to the first code record */
                  code_ptr = pdb_ptr->first_code_ptr;

                  /* Loop through all the code records to write out the
                     corresponding SWISS-PROT codes or E.C. numbers and
                     which chains they represent */
                  while (code_ptr != NULL)
                    {
                      /* If this is the right code type, write out */
                      if (code_ptr->type == type)
                        {
                          /* Increment count of this code type */
                          count++;

                          /* Write out this code */
                          if (ok == TRUE) 
                            fprintf(file_ptr,"%s[%d] %s\n",code_desc[type],
                                    count,code_ptr->value);

                          /* If this is a GO number, write out its
                             description and GO level */
                          if (type > 2)
                            {
                              /* Get description, replacing by question
                                 mark if blank */
                              go_ptr = code_ptr->go_ptr;
                              if (go_ptr != NULL)
                                {
                                  strcpy(desc,go_ptr->name);
                                  level = go_ptr->level;
                                }
                              else
                                {
                                  strcpy(desc,"?");
                                  level = 0;
                                }

                              /* Write out the description */
                              if (ok == TRUE) 
                                {
                                  fprintf(file_ptr,"GO_DESC_%s[%d] %s\n",
                                          description[type],count,desc);

                                  /* Write out the level */
                                  fprintf(file_ptr,"GO_LEVEL_%s[%d] %d\n",
                                          description[type],count,level);
                                }

                              /* If have any annotations giving origin
                                 of this GO assignment, write these out */
                              if (ok == TRUE && code_ptr->nannot > 0)
                                {
                                  /* Show how many we have */
                                  fprintf(file_ptr,"GO_NANNOT_%s[%d] %d\n",
                                          description[type],count,
                                          code_ptr->nannot);

                                  /* Loop through them all to write out */
                                  for (iannot = 0;
                                       iannot < code_ptr->nannot;
                                       iannot++)
                                    fprintf(file_ptr,
                                            "GO_ANNOT_%s[%d][%d] %s\n",
                                            description[type],count,
                                            iannot+1,
                                            code_ptr->annotation[iannot]);
                                }
                            }

                          /* Write out count of chains to which this code
                             applies */
                          if (ok == TRUE) 
                            {
                              fprintf(file_ptr,"%s_CHAINS[%d] ",
                                      code_desc[type],count);
                            }

                          /* Get pointer to the first chain record */
                          chain_ptr = pdb_ptr->first_chain_ptr;
                          ichain = 0;
                          iunique = 0;
                          chain_no = 0;

                          /* Loop through all the chain records to write
                             the chains to which this code applies */
                          while (chain_ptr != NULL)
                            {
                              /* If this is a unique chain, increment
                                 counter */
                              if (chain_ptr->copy_of == '\0')
                                iunique++;

                              /* Loop over the stored codes */
                              for (icode = 0; icode < chain_ptr->ncodes;
                                   icode++)
                                {
                                  /* If this chain has this code, then
                                     write out chain id */
                                  if (chain_ptr->code_ptr[icode]
                                      == code_ptr)
                                    {
                                      /* Increment chain counter */
                                      ichain++;

                                      /* Write out this chain */
                                      if (ok == TRUE)
                                        {
                                          if (ichain > 1)
                                            fprintf(file_ptr,",");
                                          fprintf(file_ptr," %c",
                                                  chain_ptr->chain_id);
                                        }

                                      /* If this is a unique chain, then
                                         store chain number */
                                      if (chain_ptr->copy_of == '\0')
                                        chain_no = iunique;
                                    }
                                }

                              /* Get pointer to the next chain record */
                              chain_ptr = chain_ptr->next_chain_ptr;
                            }

                          /* Close off current line */
                          if (ok == TRUE)
                            fprintf(file_ptr,"\n");

                          /* Write out which unique chain this
                             corresponds to */
                          if (ok == TRUE && chain_no > 0)
                            fprintf(file_ptr,"%s_UNIQUE_NO[%d] %d\n",
                                    code_desc[type],count,chain_no);

                          /* Write out the number of chains corresponding
                             to the current code */
                          if (ok == TRUE && ichain > 0)
                            fprintf(file_ptr,"N%s_CHAINS[%d] %d\n",
                                    code_desc[type],count,ichain);
                        }

                      /* Get pointer to the next code record */
                      code_ptr = code_ptr->next_code_ptr;
                    }
                }
            }

          /* If have any GO annotation, write out flag */
          if (ok == TRUE && any_go == TRUE)
            fprintf(file_ptr,"HAVE_GO TRUE\n");

          /* Close the seqinfo.dat file */
          if (ok == TRUE) fclose(file_ptr);
        }
    }

  /* Close the output files */
  fclose(seqdata_ptr);
  fclose(outfile_ptr);
}

/***********************************************************************

                             M   A   I   N

***********************************************************************/

int main(int argc,char *argv[])
{
  int i, nitems, nchains, ngo, nnames, npdb;

  struct c_file_data file_name_ptr;
  struct parameters params;

  struct go *first_go_ptr;
  struct pdb *first_pdb_ptr;

  struct go **go_index_ptr;
  struct pdb **pdb_index_ptr;

  /* Initialise variables */
  nitems = 0;
  ngo = 0;

  /* Read in the command-line parameters and store */
  get_command_arguments(argv,argc - 1,&params);

  /* Get directories holding PDBsum and PDB data from CATHPARAM file */
  c_get_paths(&file_name_ptr,params.run_64);

  /* Read in the PDB codes in the dataset.lst file */
  npdb = get_pdb_codes("dataset.lst",&first_pdb_ptr,params.pdb_code);

  /* If have any PDB entries, process */
  if (npdb > 0)
    {
      /* Create an index to allow binary search through PDB records */
      create_pdb_index(first_pdb_ptr,&pdb_index_ptr,npdb);

      /* Get the UniProt codes and EC numbers associated with each
         PDB code from the MSD data file pdb_chain_sp_ec */
      nitems = get_msd_data(pdb_index_ptr,npdb,file_name_ptr,params);

      /* Check for PDB to GO annotations */
/* debug */
//      if (params.pdb_type == STANDARD_PDB)
//        ngo = get_go_data(pdb_index_ptr,npdb,file_name_ptr);
      ngo = 0;
/* debug */

      /* If have read in any GO annotations, read the GO terms and
         ids file to pick up the descriptions */
      if (ngo > 0)
        {
          /* Read in the GO names */
          nnames = get_go_names(&first_go_ptr,file_name_ptr);

          /* Create an index for the GO terms just read in */
          create_go_index(first_go_ptr,&go_index_ptr,nnames);

          /* Match the GO records to those associated with each
             PDB entry */
          match_go(pdb_index_ptr,npdb,go_index_ptr,nnames);

          /* Read through the three .ontology files to calculate an
             approximate level for each GO term */
          get_go_levels(go_index_ptr,nnames,pdb_index_ptr,npdb,
                        file_name_ptr);
        }
    }

  /* If have any data from either source, then write out to the 
     seqinfo.dat and pdb_EC files */
  if (nitems > 0 || ngo > 0)
    {
      /* Write the output seqinfo.dat file */
      write_seqinfo(pdb_index_ptr,npdb,file_name_ptr,params);
    }

  return(0);
}

