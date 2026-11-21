/***********************************************************************

  promfilt.c - Program to filter out all PROMOTIF data relating to
               non-unique chains for a given PDB entry and write to
               promotif.dat for storage in the PDBsum directory

************************************************************************

Date:-         31 Mar 2005
Dir update:-   14 Nov 2012

Written by:-   Roman Laskowski

----------------------------------------------------------------------*/

/* Datafiles
   ---------

Actual name      Description
-----------      -----------
CATHPARAM        Input parameter file
chains.dat       Input chains file listing unique and duplicate
                 protein chains
prepromotif.dat  Input promotif.dat file to be filtered for unique
                 chains only
promotif.dat     Output file containing PROMOTIF data for unique chains
                 only
motifs.dat       Output file containing summary PROMOTIF data for use
                 by showmain

Compilation
-----------

cc -c promfilt.c
cc -c common.c
cc -c reapar.c
cc -o promfilt promfilt.o common.o reapar.o -lm -lz

Optimization flags:-  

cc -O2 -funroll-loops -c promfilt.c
cc -O2 -funroll-loops -c common.c
cc -O2 -funroll-loops -c reapar.c
cc -o promfilt promfilt.o common.o reapar.o -lm -lz


*/

/* Function calling-tree
   ---------------------

main
   -> get_command_arguments
   -> c_get_paths (in common.c)
         -> get_param (in reapar.c)
   -> find_pdbsum_dir (in common.c)
   -> read_chain_equivs
         -> string_truncate (in common.c)
         -> create_chain_record
   -> sort_promotif_data
         -> open_output_file (in common.c)
         -> string_truncate (in common.c)
         -> create_line_record
               -> store_string (in common.c)
         -> sort_lines
         -> free_lines
   -> filter_data
         -> string_truncate (in common.c)

*/


#include "common.h"
#include "promfilt.h"


/* Prototypes
   ---------- */

/* In common.c */
void c_get_paths(struct c_file_data *file_name_ptr,int run_64);
int find_pdbsum_dir(char *pdb_code,char *pdbsum_template[NPDBSUM_DIR],
                    int pdb_type,char *dir_name);
FILE *open_output_file(char *out_name,int terminate);
void store_string(char **string_ptr,char *string);
int string_truncate(char *string,int max_length);

/* Function prototypes from reapar.c */
int get_param(char *token_name,char *file_name,int filename_len,
              int exit_on_error);

/***********************************************************************

get_command_arguments  -  Interpret the command line arguments

***********************************************************************/

void get_command_arguments(char *strptrs[],int ntoken,
                           struct parameters *params)
{
  int itoken;

  /* Initialise parameters */
  params->pdb_code[0] = '\0';
  params->pdb_type = STANDARD_PDB;
  params->pdbsum1 = FALSE;
  params->run_64 = FALSE;

  /* If -help entered on the command line, show options available */
  if (ntoken < 1 || (ntoken == 1 && (!strcmp(strptrs[1],"-help") ||
                                     !strcmp(strptrs[1],"-h"))))
    {
      printf("\n");
      printf("Correct usage is ...\n");
      printf("\n");
      printf("    promfilt  pdb_code  [UPLOAD | MCSG | AF]  [-pdbsum1]  "
	     "[-64]\n");
      printf("\n");
      printf("where\n");
      printf("     * pdb_code = PDB code for file to be processed or uplo"
	     " if files in current dir\n");
      printf("     * UPLOAD   = Uploaded PDB file for PDBsum generation\n");
      printf("     * MCSG     = MCSG structure\n");
      printf("     * AF       = Alpha Fold structure\n");
      printf("     * -pdbsum1 = PDBsum1 structure, so sort required\n");
      printf("     * -64      = run on 64-bit processor\n");
      printf("\n");
      exit(1);
    }

  /* Get the PDB code from the first token */
  strcpy(params->pdb_code,strptrs[1]);
  params->pdb_code[4] = '\0';

  /* Loop through any additional command arguments */
  for (itoken = 2; itoken < ntoken + 1; itoken++)
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
        params->pdbsum1 = TRUE;

      /* Check for the 64-bit option */
      else if (!strcmp(strptrs[itoken],"-64"))
        params->run_64 = TRUE;
    }
}
/***********************************************************************

create_chain_record  -  Create and initialise a new chain record

***********************************************************************/

struct chain *create_chain_record(struct chain **fst_chain_ptr,
                                  struct chain **lst_chain_ptr,
                                  char chain_id,char copy_of)
{
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
  chain_ptr->chain_no = 0;
  chain_ptr->copy_of = copy_of;
  chain_ptr->next_chain_ptr = NULL;

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

int read_chain_equivs(char *pdbsum_dir,struct chain **fst_chain_ptr)
{
  char file_name[FILENAME_LEN + 1], input_line[LINELEN + 1];
  char chain_id, copy_of;

  int len, nchains, nlines, nunique;

  struct chain *chain_ptr, *first_chain_ptr, *last_chain_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  nchains = 0;
  nlines = 0;
  nunique = 0;

  /* Initialise chain pointers */
  chain_ptr = NULL;
  first_chain_ptr = NULL;
  last_chain_ptr = NULL;

  /* Form the name of the chains.dat file */
#ifdef TEST
  /* Form name of input file */
  strcpy(file_name,"chains.dat");
#else
  strcpy(file_name,pdbsum_dir);
  strcat(file_name,"chains.dat");
#endif
  printf("Reading chains.dat file (%s) ...\n",file_name);

  /* Open the chains.dat input file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    {
      printf("\n*** Warning. Unable to open input file [%s]\n",file_name);
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
          string_truncate(input_line,LINELEN);

          /* Get the line length */
          len = strlen(input_line);

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

          /* If this is a unique chain, assign a number to it */
          if (chain_ptr->copy_of == '\0')
            {
              /* Increment count of unique chains and store */
              nunique++;
              chain_ptr->chain_no = nunique;
            }

          /* Increment chain count */
          nchains++;
        }

      /* Close the input file */
      fclose(file_ptr);
    }

  /* Sho stats */
  printf("  Number of chains defined           %8d\n",nchains);
  printf("\n");

  /* Return pointer to first chain record */
  *fst_chain_ptr = first_chain_ptr;

  /* Return the number of chains stored */
  return(nunique);
}
/***********************************************************************

create_line_record  -  Create and initialise a new line record

***********************************************************************/

struct line *create_line_record(struct line **fst_line_ptr,
                                  struct line **lst_line_ptr,
                                  char *input_line)
{
  struct line *line_ptr, *first_line_ptr, *last_line_ptr;

  /* Initialise line pointers */
  line_ptr = NULL;
  first_line_ptr = *fst_line_ptr;
  last_line_ptr = *lst_line_ptr;

  /* Allocate memory for structure to hold line info */
  line_ptr = (struct line *) malloc(sizeof(struct line));
  if (line_ptr == NULL)
    {
      printf("*** Can't allocate memory for struct line\n");
      exit (1);
    }

  /* If this is the very first line, then store pointer */
  if (first_line_ptr == NULL)
    first_line_ptr = line_ptr;

  /* Add link from previous line to the current one */
  if (last_line_ptr != NULL)
    last_line_ptr->next_line_ptr = line_ptr;
  last_line_ptr = line_ptr;

  /* Store the line details */
  store_string(&(line_ptr->line),input_line);

  /* Initialise pointer to the next record */
  line_ptr->next_line_ptr = NULL;

  /* Return current pointers */
  *fst_line_ptr = first_line_ptr;
  *lst_line_ptr = last_line_ptr;

  /* Return new line pointer */
  return(line_ptr);
}
/***********************************************************************

line_sort  -  Shell-sort which sorts the line types into ascending
               order ot type, copy and duplicate

***********************************************************************/

void line_sort(int *index,struct line **line_index_ptr,int nindex)
{
  int endloop, i, indi, indl, iswap, j, k, l, lognb2, m, n, nn;
  float aln2i = (float) 1.4426950, tiny = (float) 1.e-5;
  float calc;
  double dble, dble1;

  struct line *line1_ptr, *line2_ptr;

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

              /* Get the two line records */
              line1_ptr = line_index_ptr[index[indi]];
              line2_ptr = line_index_ptr[index[indl]];

              /* If in wrong order, the swap sort-pointers */
              if (strcmp(line1_ptr->line,line2_ptr->line) > 0)
                {
                  iswap = index[indi];
                  index[indi] = index[indl];
                  index[indl] = iswap;

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

sort_lines  -  Sort the lines read in from the Promotif file

***********************************************************************/

int sort_lines(struct line **fst_line_ptr,int nlines)
{
  int iline, iorder, ipos;
  int *index;

  struct line *first_line_ptr, *last_line_ptr;
  struct line *line_ptr;
  struct line **line_index_ptr;

  /* Initialise pointers */
  first_line_ptr = *fst_line_ptr;
  last_line_ptr = NULL;

  /* Grab memory for temporary sort arrays */
  index = (int *) malloc(nlines * sizeof(int));
  line_index_ptr
    = (struct line **) malloc(nlines * sizeof(struct line *));

  /* Get pointer to first line record */
  line_ptr = first_line_ptr;
  ipos = 0;

  /* Loop through line records to initialise the sort index */
  while (line_ptr != NULL)
    {
      /* Initialise sort index */
      index[ipos] = ipos;
      line_index_ptr[ipos] = line_ptr;

      /* Increment line counter */
      ipos++;

      /* Get next line record */
      line_ptr = line_ptr->next_line_ptr;
    }

  /* Initialise line counter */
  iline = 0;
  nlines = ipos;

  /* Sort the lines */
  line_sort(index,line_index_ptr,nlines);

  /* Having sorted the lines, rearrange order of the line-record
     pointers accordingly */
  for (iorder = 0; iorder < nlines; iorder++)
    {
      /* Get the array position of this line record */
      ipos = index[iorder];

      /* Get the current line record */
      line_ptr = line_index_ptr[ipos];

      /* If this is the first, then make the first record */
      if (iorder == 0)
        first_line_ptr = line_ptr;

      /* Otherwise, get the previous line record to point to this one */
      else
        last_line_ptr->next_line_ptr = line_ptr;

      /* Blank out current record's next line pointer */
      line_ptr->next_line_ptr = NULL;

      /* Store current record as last encountered */
      last_line_ptr = line_ptr;
    }

  /* Free up memory used for temporary sort arrays */
  free(index);
  free(line_index_ptr);

  /* Return current pointer */
  *fst_line_ptr = first_line_ptr;

  /* Return number of sorted lines */
  return(nlines);
}
/***********************************************************************

free_line_memory  -  Free up memory used by residue line

***********************************************************************/

void free_line_memory(struct line *first_line_ptr)
{
  struct line *line_ptr, *next_line_ptr;

  /* Get first line */
  line_ptr = first_line_ptr;

  /* Loop through the line records until done */
  while (line_ptr != NULL)
    {
      /* Get pointer to the next line record */
      next_line_ptr = line_ptr->next_line_ptr;

      /* If line not blank, then free the string */
      if (line_ptr->line != &null)
	free(line_ptr->line);

      /* Free up memory taken by current line */
      free(line_ptr);

      /* Go to the next line */
      line_ptr = next_line_ptr;
    }
}
/***********************************************************************

sort_promotif_data  -  Sort all the individual Promotif files and output
                       a prepromotif.dat file

***********************************************************************/

void sort_promotif_data(void)
{
  char file_name[FILENAME_LEN + 1], input_line[LINELEN + 1];

  static char *fname[NMOTIF_TYPES + 1] = { 
    "barrels", "sheets", "bab", "hairpins", "psi", "bulges", "strands",
    "helices", "helint", "bturns", "gturns", "disulph", "summary" };

  int ifile, iline, nlines, nchains;

  struct line *line_ptr, *first_line_ptr, *last_line_ptr;

  FILE *file_ptr, *outfile_ptr;

  /* Initialise variables */
  nchains = 0;

  /* Open the output prepromotif.dat file */
  strcpy(file_name,"prepromotif.dat");
  outfile_ptr = open_output_file(file_name,TRUE);

  /* Loop through all the files */
  for (ifile = 0; ifile < NMOTIF_TYPES + 1; ifile++)
    {
      /* Form the name of this file */
      sprintf(file_name,"%s.promotif",fname[ifile]);

      /* Initialise line count */
      nlines = 0;

      /* Initialise pointers */
      line_ptr = first_line_ptr = last_line_ptr = NULL;

      /* Open the input file */
      if ((file_ptr = fopen(file_name,"r")) !=  NULL)
	{
	  /* Loop while reading the file */
	  while (fgets(input_line,LINELEN,file_ptr)!= NULL)
	    {
	      /* Truncate the string to strip off the newline */
	      string_truncate(input_line,LINELEN);

	      /* Increment the line count */
	      nlines++;

	      /* If one of the header lines, write straight out */
	      if (nlines < 3)
		fprintf(outfile_ptr,"%s\n",input_line);

	      /* Otherwise, store the line for sorting */
	      else
		{
		  /* Create a line record */
		  line_ptr
		    = create_line_record(&first_line_ptr,&last_line_ptr,
					 input_line);

		  /* Increment count of lines */
		  nlines++;
		}
	    }

	  /* Close the input file */
	  fclose(file_ptr);

	  /* Sort all the stored lines */
	  sort_lines(&first_line_ptr,nlines);

	  /* Get the first line record */
	  line_ptr = first_line_ptr;

	  /* Loop to write out */
	  while (line_ptr != NULL)
	    {
	      /* Write out this line */
	      fprintf(outfile_ptr,"%s\n",line_ptr->line);

	      /* Get pointer to the next line record */
	      line_ptr = line_ptr->next_line_ptr;
	    }

	  /* Free memory */
	  free_line_memory(first_line_ptr);
	}
    }

  /* Close the output file */
  fclose(outfile_ptr);
}
/***********************************************************************

filter_data  -  Write the output seqinfo.dat file

***********************************************************************/

void filter_data(struct chain *first_chain_ptr,int nunique)
{
  char ch, chain, last_chain;
  char file_name[LINELEN + 1], input_line[LINELEN + 1];
  char motif_name[30];

  int chain_no, i, last_no, nmotifs, ntypes[MAX_CHAINS], type;
  int found, file_open, new_chain, unique, wanted;

  struct chain *chain_ptr;

  static char *motif_names[NMOTIF_TYPES] = { 
    "BETA_BARRELS", "SHEETS", "BAB_MOTIFS", "BETA_HAIRPINS", "PSI_LOOPS",
    "BETA_BULGES", "STRANDS", "HELICES", "HELIX_INTERACTIONS",
    "BETA_TURNS", "GAMMA_TURNS", "DISULPHIDES" };

  FILE *file_ptr, *outfile_ptr, *motifile_ptr;

  /* Initialise variables */
  file_open = FALSE;
  last_chain = '\0';
  last_no = 0;
  chain_no = 0;
  new_chain = TRUE;
  nmotifs = 0;
  for (i = 0; i < MAX_CHAINS; i++)
    ntypes[i] = 0;
  motif_name[0] = '\0';
  type = -1;
  unique = FALSE;

  /* Form input file name */
  strcpy(file_name,"prepromotif.dat");

  /* Open the input file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    {
      printf("\n*** Warning. Unable to open input file [%s]\n",file_name);
      exit(0);
    }

  /* Open the output promotif.dat file */
  strcpy(file_name,"promotif.dat");
  outfile_ptr = open_output_file(file_name,TRUE);

  /* Loop while reading in records from the input file */
  while (fgets(input_line,LINELEN,file_ptr)!= NULL)
    {
      /* Truncate the string to strip off the newline */
      string_truncate(input_line,LINELEN);

      /* Assume line not wanted */
      wanted = FALSE;

      /* Get the first character on the line */
      ch = input_line[0];

      /* If this is a START record or a new chain record, write out
         number of motifs for previous chain */
      if ((!strncmp(input_line,"START:",6) ||
          (!strncmp(input_line+1,"::",2) && ch != last_chain) ||
           !strcmp(input_line,"EOF")) &&
          nmotifs > 0)
        {
          /* If a type to be listed on the PDBsum page, write out */
          if (chain_no > 0 && type + 1 > 0)
            {
              /* If not yet open, open summary data file */
              if (file_open == FALSE)
                {
                  strcpy(file_name,"motifs.dat");
                  motifile_ptr = open_output_file(file_name,TRUE);

                  /* Write start record to motifs.dat file */
                  fprintf(motifile_ptr,"NEW_PROMOTIF TRUE\n");

                  /* Write start record to motifs.dat file */
                  fprintf(motifile_ptr,"PROMOTIF TRUE\n");

                  /* Set flag that file open */
                  file_open = TRUE;
                }

              /* Increment count of this chain's motif types */
              if (chain_no < MAX_CHAINS + 1)
                ntypes[chain_no - 1]++;

              /* Get the motif's place in this chain's list of motifs */
              i = ntypes[chain_no - 1];

              /* Determine the motif type and write out */
              fprintf(motifile_ptr,"P_TYPE[%d][%d] %d\n",chain_no,i,
                      type + 1);

              /* Write out number of motifs for previous chain */
              fprintf(motifile_ptr,"P_COUNT[%d][%d] %d\n",chain_no,i,
                      nmotifs);
            }

          /* Re-initialise count */
          nmotifs = 0;
        }

      /* If this is a START record, pick up the new motif name */
      if (!strncmp(input_line,"START:",6))
        {
          /* Get the motif name */
          strcpy(motif_name,input_line+6);

          /* Determine motif type */
          type = -1;
          for (i = 0; i < NMOTIF_TYPES; i++)
            if (!strcmp(motif_name,motif_names[i]))
              type = i;

          /* Re-initialise motif count */
          nmotifs = 0;

          /* Flag as wanted */
          wanted = TRUE;
        }

      /* If this is a data record, check whether chain is a wanted one */
      else if (!strncmp(input_line+1,"::",2))
        {
          /* Get the chain */
          chain = ch;

          /* If this is a new chain, need to determine whether it is
             a new chain or a copy chain */
          if (chain != last_chain)
            {
              /* Assume is unique */
              found = FALSE;
              unique = FALSE;

              /* Get pointer to the first chain record */
              chain_ptr = first_chain_ptr;

              /* Loop through all the chain records to determine
                 whether this is  a unique or a copy chain */
              while (chain_ptr != NULL && found == FALSE)
                {
                  /* If this is our chain, check whether it is a copy
                     chain */
                  if (chain_ptr->chain_id == chain)
                    {
                      /* If this is a copy chain flag it as non-unique */
                      if (chain_ptr->copy_of != '\0')
                        unique = FALSE;
		      else
                        unique = TRUE;

                      /* Mark as found */
                      found = TRUE;

                      /* Store the chain number */
                      chain_no = chain_ptr->chain_no;
                    }

                  /* Get pointer to the next chain record */
                  chain_ptr = chain_ptr->next_chain_ptr;
                }
            }

          /* If record corresponds to a unique chains, then want to
             keep it */
          if (unique == TRUE)
            {
              /* Set flag */
              wanted = TRUE;

              /* Increment count of motifs written out */
              nmotifs++;
            }

          /* For disulphides and helix interactions need to keep all
             records for the protein-protein interface analyses */
          if (!strcmp(motif_name,"DISULPHIDES") ||
              !strcmp(motif_name,"HELIX_INTERACTIONS"))
            wanted = TRUE;

          /* Store chain */
          last_chain = chain;
        }

      /* Otherwise, take this to be the field names line */
      else
        wanted = TRUE;

      /* If line is wanted, write out to the promotif.dat file */
      if (wanted == TRUE)
        {
          /* Write out */
          fprintf(outfile_ptr,"%s\n",input_line);
        }
    }

  /* Write the motif counts for each chain */
  if (file_open == TRUE)
    {
      for (i = 0; i < nunique; i++)
        fprintf(motifile_ptr,"NMOTIFS[%d] %d\n",i + 1,ntypes[i]);

      /* Close file */
      fclose(motifile_ptr);
    }

  /* Close the other files */
  fclose(file_ptr);
  fclose(outfile_ptr);
}

/***********************************************************************

                             M   A   I   N

***********************************************************************/

int main(int argc,char *argv[])
{
  char pdbsum_dir[FILENAME_LEN + 1];

  int nunique;
  int dir_type;

  struct c_file_data file_name_ptr;
  struct parameters params;

  struct chain *first_chain_ptr;

  /* Interpret the command-line parameters */
  get_command_arguments(argv,argc - 1,&params);

  /* Get directories holding PDBsum and PDB data from CATHPARAM file */
  c_get_paths(&file_name_ptr,params.run_64);

  /* Form the name of the PDBsum directory */
  if (params.pdbsum1 == TRUE)
    strcpy(pdbsum_dir,"../newsec/");
  else if (strcmp(params.pdb_code,"uplo"))
    {
      dir_type
        = find_pdbsum_dir(params.pdb_code,file_name_ptr.pdbsum_template,
                          params.pdb_type,pdbsum_dir);
      if (dir_type == NOT_FOUND)
        {
          printf("*** ERROR. Unable to locate PDBsum directory for %s\n",
                 params.pdb_code);
          exit(-1);
        }
    }
  else
    pdbsum_dir[0] = '\0';

  /* Read in all the chain-equivalences from the chains.dat file */
  nunique = read_chain_equivs(pdbsum_dir,&first_chain_ptr);

  /* If we have a PDBsum1 structure, need to sort all the individual
     Promotif files and output a prepromotif.dat file */
  if (params.pdbsum1 == TRUE)
    sort_promotif_data();

  /* Process the prepromotif.dat file to extract all the required data */
  filter_data(first_chain_ptr,nunique);

  return(0);
}

