/***********************************************************************

  resdata.c - Program to pick up various residue annotations (such as
              PROSITE patterns, PROMOTIF motifs, CATH domains, etc)
              required by SAS and stored in the entry's PDBsum directory
              as resdata.dat

If C-alpha only, but have secondary structure, show in grey
                 if no sst assignments, don't plot secondary structure
                 part

Get protein1.html to show sequence insetad!

************************************************************************

Date:-         7 Apr 2004
Dir update:-   16 Oct 2012

Written by:-   Roman Laskowski

----------------------------------------------------------------------*/

/* Datafiles
   ---------

Actual name      Description
-----------      -----------
CATHPARAM        Input parameter file
cath.domains     Input CATH domains file giving domain definitions for
                 all PDB chains classified in CATH
chains.dat       Input chains file listing unique and duplicate
                 protein chains
grow.out         Input file listing H-bonds and non-bonded interactions
                 that all protein residues make with any other molecules
                 in the structure
ligands.dat      Input file listing all ligands and metal ions
pdbXXXX.ent      Input PDB file
promotif.dat     Input .dat file holding structure's PROMOTIF data
rasmol.dfn       Input file holding the RasMol script definitions from
                 which the PDB file's contents are determined
resdata.dat      Output file listing all the protein residues and
                 their annotations
3dp.matches.list Input file listing all PROSITE matches to PDB chains
3dp.tof.list     Input file marking PROSITE matches as true or false


Compilation
-----------

cc -c resdata.c
cc -c common.c
cc -c reapar.c
cc -c readpdb.c
cc -o resdata resdata.o common.o reapar.o readpdb.o -lm -lz

Optimization flags:-  

cc -O2 -funroll-loops -c resdata.c
cc -O2 -funroll-loops -c common.c
cc -O2 -funroll-loops -c reapar.c
cc -O2 -funroll-loops -c readpdb.c
cc -o resdata resdata.o common.o reapar.o readpdb.o -lm -lz

*/

/* Function calling-tree
   ---------------------

main
   -> get_command_arguments
   -> c_get_paths (in common.c)
         -> get_param (in reapar.c)
         -> store_string (in common.c)
   -> find_pdbsum_dir (in common.c)
   -> get_dfn_data
         -> string_truncate (in common.c)
         -> create_rasdef_entry
               -> store_string (in common.c)
         -> get_ligand_residues
               -> get_token
               -> format_res_number (in common.c)
               -> create_rasdef_entry
                     -> store_string (in common.c)
   -> get_ligands
         -> string_truncate (in common.c)
   -> define_rep_ligands
   -> get_hets
         -> string_truncate (in common.c)
         -> create_het_record
   -> read_chain_equivs
         -> string_truncate (in common.c)
         -> create_chain_record
   -> read_pdb_file
         -> form_filename (in common.c)
         -> string_truncate (in common.c)
         -> get_atom_details
         -> assign_to_object
         -> create_residue_record
               -> get_aa_code (in common.c)
         -> create_atom_record
         -> scan_sites
         -> create_site_record
   -> get_promotif_data 
         -> read_promotif_dat
               -> string_truncate (in common.c)
               -> create_promotif_record
         -> read_in_promotif_data
               -> string_truncate (in common.c)
               -> create_motif_record
         -> motif_annotate
         -> free_motif_memory
   -> get_prosite_patterns
         -> string_truncate (in common.c)
         -> get_prosite_strings
         -> check_for_duplicate_pattern
         -> create_pattern_record
               -> store_string (in common.c)
         -> pattern_annotate
   -> verify_prosite_patterns
         -> string_truncate (in common.c)
         -> pattern_annotate
   -> read_in_domain_data
         -> string_truncate (in common.c)
         -> extract_domain_info (in common.c)
               -> get_tokens (in common.c)
               -> format_residue_number (in common.c)
               -> free_tokens (in common.c)
         -> extract_domain_data -- No longer used
               -> get_domain_token
               -> format_residue_number
         -> domain_annotate
   -> get_grow_data
         -> read_in_grow_data
               -> string_truncate (in common.c)
               -> create_grow_record
         -> add_rasdef
   -> write_resdata
         -> open_output_file (in common.c)
   -> write_domains
         -> open_output_file (in common.c)

*/

#include "common.h"
#include "resdata.h"

/* Prototypes
   ---------- */

/* In common.c */
void c_get_paths(struct c_file_data *file_name_ptr,int run_64);
int extract_domain_info(char *input_line,char *wolf_code,char *cath_code,
                        char domain_details[C_MAXDOMCHAIN][16]);
int find_pdbsum_dir(char *pdb_code,char *pdbsum_template[NPDBSUM_DIR],
                    int pdb_type,char *dir_name);
void format_residue_number(char residue_number[6],char *token);
char get_aa_code(char *res_name);
FILE *open_output_file(char *out_name,int terminate);
void store_string(char **string_ptr,char *string);
int string_truncate(char *string,int max_length);

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

void get_command_arguments(char *strptrs[],int ntoken,
                           struct parameters *params)
{
  int itoken;

  /* Initialise parameters */
  params->pdb_code[0] = '\0';
  params->pdb_type = STANDARD_PDB;
  params->run_64 = FALSE;

  /* If -help entered on the command line, show options available */
  if (ntoken < 1 || (ntoken == 1 && (!strcmp(strptrs[1],"-help") ||
                                     !strcmp(strptrs[1],"-h"))))
    {
      printf("\n");
      printf("Correct usage is ...\n");
      printf("\n");
      printf("    resdata  {pdb_code | -uplo}  [UPLOAD | MCSG | AF]  [-64]"
             "  [-pdbsum1]\n");
      printf("where\n");
      printf("     * pdb_code = PDB code for file to be processed\n");
      printf("     * -uplo    = File is an uploaded PDB file\n");
      printf("     * UPLOAD   = Uploaded PDB file for PDBsum generation\n");
      printf("     * MCSG     = MCSG structure\n");
      printf("     * AF       = Alpha Fold structure\n");
      printf("     * -pdbsum1 = PDBsum1 structure\n");
      printf("     * -64      = run on 64-bit processor\n");
      printf("\n");
      exit(1);
    }

  /* Check for uploaded PDB file flag */
  if (!strcmp(strptrs[1],"-uplo"))
    strcpy(params->pdb_code,"uplo");

  /* Otherwise, take first parameter to be the PDB code */
  else
    strcpy(params->pdb_code,strptrs[1]);

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
        params->pdb_type = PDBSUM1;

      /* Check for the 64-bit option */
      else if (!strcmp(strptrs[itoken],"-64"))
        params->run_64 = TRUE;
    }
}
/***********************************************************************

create_rasdef_entry  -  Create a rasdef entry

***********************************************************************/

struct rasdef *create_rasdef_entry(struct rasdef **fst_rasdef_ptr,
                                   struct rasdef **lst_rasdef_ptr,
                                   int type,int number,char *res_name,
                                   char *res_num,char chain,
                                   char *colour,char *ligand_def,
                                   char *ligand_from,char *ligand_to,
                                   int count)
{
  int i;

  struct rasdef *rasdef_ptr;
  struct rasdef *first_rasdef_ptr, *last_rasdef_ptr;

  /* Initialise rasdef pointers */
  rasdef_ptr = NULL;
  first_rasdef_ptr = *fst_rasdef_ptr;
  last_rasdef_ptr = *lst_rasdef_ptr;

  /* Create rasdef entry */
  rasdef_ptr =
    (struct rasdef*) malloc(sizeof(struct rasdef));
  if (rasdef_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for");
      printf(" struct rasdef\n");
      printf("***        Program postssm terminated with");
      printf(" error.\n");
      exit (1);
    }

  /* Store the details */
  rasdef_ptr->type = type;
  rasdef_ptr->number = number;
  strcpy(rasdef_ptr->res_name,res_name);
  strcpy(rasdef_ptr->res_num,res_num);
  rasdef_ptr->chain = chain;
  strcpy(rasdef_ptr->colour,colour);
  rasdef_ptr->count = count;

  /* If this is a ligand entry, store the RasMol definition of the
     ligand */
  if (type == LIGAND)
    {
      /* Store the string */
      store_string(&(rasdef_ptr->rasmol_ligand_def),ligand_def);
    }
  else
    rasdef_ptr->rasmol_ligand_def = &null;

  /* Initialise other fields */
  strcpy(rasdef_ptr->ligand_from,ligand_from);
  strcpy(rasdef_ptr->ligand_to,ligand_to);
  for (i = 0; i < MAX_INT_CHAINS; i++)
    rasdef_ptr->interacts_with[i] = ' ';
  rasdef_ptr->nhets = 0;
  rasdef_ptr->first_het_ptr = NULL;
  rasdef_ptr->last_het_ptr = NULL;
  rasdef_ptr->range_rasdef_ptr = NULL;

  /* Initialise pointer to next entry */
  rasdef_ptr->next_rasdef_ptr = NULL;

  /* If this is the first entry stored, then update
     pointer to head of linked list */
  if (first_rasdef_ptr == NULL)
    first_rasdef_ptr = rasdef_ptr;

  /* Otherwise, make previous entry point to the current
     one */
  else
    last_rasdef_ptr->next_rasdef_ptr = rasdef_ptr;
                      
  /* Save the current residue's pointer */
  last_rasdef_ptr = rasdef_ptr;

  /* Return current pointers */
  *fst_rasdef_ptr = first_rasdef_ptr;
  *lst_rasdef_ptr = last_rasdef_ptr;
  return(rasdef_ptr);
}
/***********************************************************************

get_token - Get the next space- or tab-delimited token from the given
            input line, starting at the given position

**********************************************************************/  

int get_token(char *line,int line_len,int *ipos,char *token,
              int max_token_len,int *done)
{
  char ch;

  int i, token_len;
  int have_token;

  /* Initialise variables */
  *done = FALSE;
  have_token = FALSE;
  i = 0;
  token[0] = '\0';

  /* Loop until token-string found or end of line reached */
  while (*done == FALSE && have_token == FALSE && *ipos < line_len)
    {
      /* Get the next character in the string */
      ch = line[*ipos];

      /* Check for end of line */
      if (ch == '\0' || ch == '\n')
        *done = TRUE;

      /* If have a space, then take to be the end of the token if have
         already got something in it */
      else if (ch == ' ')
        {
          /* If have a token, then we are done */
          if (token[0] != '\0')
            have_token = TRUE;
        }

      /* Otherwise, add character to token-string provided not too long */
      else
        {
          if (i < max_token_len)
            {
              token[i] = ch;
              i++;
            }
        }

      /* Go to the next character in the line */
      (*ipos)++;
    }

  /* Terminate the token-string */
  token[i] = '\0';
  token_len = i;
  if (i > 0)
    *done = FALSE;

  /* Return the token */
  return(token_len);
} 
/***********************************************************************

format_res_number  -  Interpret the residue-number entered by the
                       user

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
  for (idigit = 0; idigit < 5; idigit++)
    {
      new_number[idigit] = ' ';
      res_number[idigit] = ' ';
    }

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

get_ligand_residues  -  Extract the ligand details on the current line

***********************************************************************/

struct rasdef *get_ligand_residues(char *input_line,
                                   struct rasdef **fst_rasdef_ptr,
                                   struct rasdef **lst_rasdef_ptr,int number,
                                   int *nd)
{
  char chain, chain_str[2], res_name[4], res_num[6];
  char colour[COLNAM_LEN], number_string[10];
  char token_string[TOKEN_LEN + 1], tmp_string[TOKEN_LEN + 1];
  char *string_ptr;
  char ligand_from[7], ligand_to[7];
  char first_residue[7], last_residue[7];

  int count, len, lpos, ndef, token_len, type;
  int name_end;
  int done;

  struct rasdef *rasdef_ptr, *range_rasdef_ptr, *start_rasdef_ptr;
  struct rasdef *first_rasdef_ptr, *last_rasdef_ptr;

  /* Initialise variables */
  chain = ' ';
  count = 0;
  ndef = *nd;
  done = FALSE;
  type = LIGAND;
  colour[0] = '\0';
  first_residue[0] = '\0';
  last_residue[0] = '\0';

  /* Initialise rasdef pointers */
  rasdef_ptr = range_rasdef_ptr = start_rasdef_ptr = NULL;
  first_rasdef_ptr = *fst_rasdef_ptr;
  last_rasdef_ptr = *lst_rasdef_ptr;

  /* Get the line length */
  len = strlen(input_line);
  lpos = 0;

  /* Loop until line has been processed */
  while (done == FALSE)
    {
      /* Get the next token on the line */
      token_len = get_token(input_line,LINELEN,&lpos,token_string,
                            TOKEN_LEN,&done);

      /* If have token, then interpret it */
      if (token_len > 1)
        {
          /* Check for square brackets defining residue name */
          if (token_string[0] == '[')
            {
              /* Get the name bounded by square brackets */
              strcpy(tmp_string,token_string+1);

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
              strncpy(res_name,token_string,3);
              res_name[3] = '\0';
              name_end = 2;
            }

          /* Get the residue number */
          strcpy(number_string,token_string + name_end + 1);
          string_ptr = strstr(number_string,",");
          if (string_ptr != NULL)
            string_ptr[0] = '\0';
          string_ptr = strstr(number_string,":");
          if (string_ptr != NULL)
            {
              /* Extract the chain identifier after the colon */
              chain = string_ptr[1];
              string_ptr[0] = '\0';
            }

          /* Convert the residue number into the standard fixed format */
          format_res_number(number_string,res_num);

          /* Copy chain into chain-string */
          chain_str[0] = chain;
          chain_str[1] = '\0';

          /* If this is the first residue of this ligand, store as first */
          if (first_residue[0] == '\0')
            {
              /* Copy across and add chain-id */
              strcpy(first_residue,number_string);
              strcat(first_residue,chain_str);
            }

          /* Copy residue number and chain-id across as last residue of
             current ligand */
          strcpy(last_residue,number_string);
          strcat(last_residue,chain_str);

          /* Increment count of residues belonging to this ligand */
          count++;

          /* Blank out unwanted fields */
          ligand_from[0] = '\0';
          ligand_to[0] = '\0';

          /* Create a rasdef entry for this ligand residue */
          rasdef_ptr
            = create_rasdef_entry(&first_rasdef_ptr,
                                  &last_rasdef_ptr,type,number,
                                  res_name,res_num,chain,colour,
                                  input_line,ligand_from,ligand_to,
                                  count);

          /* If this is the first record, store its pointer */
          if (start_rasdef_ptr == NULL)
            start_rasdef_ptr = rasdef_ptr;
        }
    }

  /* If have a valid residue range for this ligand, create an extra
     rasdef record to store the range */
  if (first_residue[0] != '\0')
    {
      /* Define type as a range of residues representing a ligand */
      type = LIGAND_RANGE;

      /* Blank out any unwanted data */
      res_name[0] = '\0';
      res_num[0] = '\0';
      colour[0] = '\0';

      /* Create a rasdef entry for this ligand residue */
      range_rasdef_ptr
        = create_rasdef_entry(&first_rasdef_ptr,
                              &last_rasdef_ptr,type,number,
                              res_name,res_num,chain,colour,
                              input_line,ligand_from,ligand_to,
                              count);
      ndef++;

      /* Get pointer to first of the residue records defined above */
      rasdef_ptr = start_rasdef_ptr;

      /* Loop through the latest set of residue records to get them to
         pount to the range record */
      while (rasdef_ptr != NULL)
        {
          /* Store pointer */
          rasdef_ptr->range_rasdef_ptr = range_rasdef_ptr;

          /* Get next definition record */
          rasdef_ptr = rasdef_ptr->next_rasdef_ptr;
        }
    }

  /* Return current pointers */
  *fst_rasdef_ptr = first_rasdef_ptr;
  *lst_rasdef_ptr = last_rasdef_ptr;

  /* Return updated number of rasdef entries */
  *nd = ndef;

  /* Return pointer to residue range record */
  return(range_rasdef_ptr);
}
/***********************************************************************

get_dfn_data  -  Read in and store the RasMol definitions from the
                 rasmol.dfn file

***********************************************************************/

int get_dfn_data(char *pdbsum_dir,struct rasdef **fst_rasdef_ptr,
                 struct parameters params)
{
  char input_line[LINELEN + 1], file_name[FILENAME_LEN];
  char chain, res_name[4], res_num[6];
  char colour[COLNAM_LEN], tmp_line[LINELEN + 1];
  char ligand_from[7], ligand_to[7];

  char *string_ptr;

  int count, len, line, nchain, ndef, nligand, nmetal, type;
  int first;

  struct rasdef *rasdef_ptr, *range_rasdef_ptr;
  struct rasdef *first_rasdef_ptr, *last_rasdef_ptr;

  FILE *file_ptr;

  /* Initialise */
  first = TRUE;
  count = 0;
  ligand_from[0] = '\0';
  ligand_to[0] = '\0';
  line = 0;
  nchain = 0;
  ndef = 0;
  nligand = 0;
  nmetal = 0;

  /* Form the file name of the rasmol.dfn file */
  strcpy(file_name,pdbsum_dir);
  if (params.pdb_type == PDBSUM1)
    strcat(file_name,"newsec/");
  strcat(file_name,"rasmol.dfn");

  /* If have an uploaded PDB file, the file will be in the current directory */
  if (!strcmp(params.pdb_code,"uplo"))
    strcpy(file_name,"rasmol.dfn");

  /* Initialise rasdef pointers */
  rasdef_ptr = NULL;
  first_rasdef_ptr = NULL;
  last_rasdef_ptr = NULL;
  range_rasdef_ptr = NULL;

  /* Open the rasmol.dfn input file */
  if ((file_ptr = fopen(file_name,"r")) !=  NULL)
    {
      /* Loop through the file processing each record */
      while (fgets(input_line,LINELEN,file_ptr)!= NULL)
        {
          /* Increment line-count */
          line++;

          /* Truncate the string to strip off the newline */
          string_truncate(input_line,LINELEN);

          /* Check line-type according to 1st character: C, L, M,
             N or P */

          /* If protein CA-only chain, then flag as such */
          if (input_line[0] == 'C')
            {
              /* Define type */
              type = CA_ONLY;
              nchain++;

              /* Define residue details */
              res_name[0] = '\0';
              res_num[0] = '\0';
              chain = input_line[1];
              if (chain == '\0' || chain == '\n')
                chain = ' ';
              colour[0] = '\0';

              /* Create a rasdef entry */
              rasdef_ptr
                = create_rasdef_entry(&first_rasdef_ptr,
                                      &last_rasdef_ptr,type,nchain,
                                      res_name,res_num,chain,colour,
                                      input_line,ligand_from,ligand_to,
                                      count);
              ndef++;
            }

          /* Check for ligand identifier */
          else if (input_line[0] == 'L')
            {
              /* Increment ligand-count */
              nligand++;

              /* Extract the ligand residues */
              range_rasdef_ptr
                = get_ligand_residues(input_line+4,&first_rasdef_ptr,
                                      &last_rasdef_ptr,nligand,&ndef);
            }


          /* Check for ligand range */
          else if (input_line[0] == 'l' && range_rasdef_ptr != NULL)
            {
              /* Store the ligand residue range */
              store_string(&(range_rasdef_ptr->rasmol_ligand_def),
                           input_line+5);

              /* Re-initialise pointer */
              range_rasdef_ptr = NULL;
            }

          /* Check for metal */
          else if (input_line[0] == 'M')
            {
              /* Define type */
              type = METAL;
              nmetal++;
              res_name[0] = '\0';

              /* Store the metal name, right-justified */
              strcpy(tmp_line,input_line+1);
              string_ptr = strchr(tmp_line,' ');
              if (string_ptr != NULL)
                string_ptr[0] = '\0';
              len = strlen(tmp_line);
              if (len == 1)
                {
                  strcpy(res_name,"  ");
                  strcat(res_name,tmp_line);
                }
              else if (len == 2)
                {
                  strcpy(res_name," ");
                  strcat(res_name,tmp_line);
                }
              else
                {
                  strncpy(res_name,tmp_line,3);
                  res_name[3] = '\0';
                }

              /* Extract the metal name into the residue number field */
              strncpy(res_num,input_line+1,3);
              res_num[3] = ' ';
              res_num[4] = '\0';

              /* Get the metal colour */
              if (input_line[4] == ' ')
                strcpy(colour,input_line+5);
              else
                strcpy(colour,input_line+4);

              /* Create a rasdef entry */
              rasdef_ptr
                = create_rasdef_entry(&first_rasdef_ptr,
                                      &last_rasdef_ptr,type,nmetal,
                                      res_name,res_num,chain,colour,
                                      input_line,ligand_from,ligand_to,
                                      count);
              ndef++;
            }
          
          /* Check for nucleic acid chain */
          else if (input_line[0] == 'N')
            {
              /* Define type */
              type = NUCLEIC_ACID;
              nchain++;

              /* Define residue details */
              res_name[0] = '\0';
              res_num[0] = '\0';
              chain = input_line[1];
              if (chain == '\0' || chain == '\n')
                chain = ' ';
              colour[0] = '\0';

              /* Create a rasdef entry */
              rasdef_ptr
                = create_rasdef_entry(&first_rasdef_ptr,
                                      &last_rasdef_ptr,type,nchain,
                                      res_name,res_num,chain,colour,
                                      input_line,ligand_from,ligand_to,
                                      count);
              ndef++;
            }

          /* If protein chain, then write the script records to define
             and display the chain */
          else if (input_line[0] == 'P')
            {
              /* Define type */
              type = PROTEIN;
              nchain++;

              /* Define residue details */
              res_name[0] = '\0';
              res_num[0] = '\0';
              chain = input_line[1];
              if (chain == '\0' || chain == '\n')
                chain = ' ';
              colour[0] = '\0';

              /* Create a rasdef entry */
              rasdef_ptr
                = create_rasdef_entry(&first_rasdef_ptr,
                                      &last_rasdef_ptr,type,nchain,
                                      res_name,res_num,chain,colour,
                                      input_line,ligand_from,ligand_to,
                                      count);
              ndef++;
            }
        }

      /* Close the rasdef file */
      fclose(file_ptr);
    }

  /* If no rasmol.dfn file, then can't continue */
  else
    {
      /* Show error message and stop */
      printf("*** ERROR. Unable to open file: %s\n",file_name);
      exit(-1);
    }

  /* Return pointer to first RasMol definition */
  *fst_rasdef_ptr = first_rasdef_ptr;

  /* Show number of records read in */
  printf("Number of definitions read in from rasmol.dfn file:  %8d\n",
         ndef);

  /* Return number of definition entries read in */
  return(ndef);
}
/***********************************************************************

get_ligands  -  Read the ligands.dat file to pick up the ligand names
                and replacing current ligand ranges by the ligand
                sequences

***********************************************************************/

int get_ligands(char *pdbsum_dir,struct rasdef *first_rasdef_ptr,
                struct parameters params)
{
  char input_line[LINELEN + 1], file_name[FILENAME_LEN];
  char desc[LINELEN + 1], seq_name[LINELEN + 1];

  char *string_ptr;

  int len, line, nligands;

  struct rasdef *rasdef_ptr;

  FILE *file_ptr;

  /* Initialise */
  seq_name[0] = '\0';
  line = 0;
  nligands = 0;

  /* Form the file name of the rasmol.dfn file */
  strcpy(file_name,pdbsum_dir);
  if (params.pdb_type == PDBSUM1)
    strcat(file_name,"newsec/");
  strcat(file_name,"ligands.dat");

  /* If have an uploaded PDB file, the file will be in the current
     directory */
  if (!strcmp(params.pdb_code,"uplo"))
    strcpy(file_name,"ligands.dat");

  /* Open the ligands.dat input file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    return(nligands);

  /* Loop through the file processing each record */
  while (fgets(input_line,LINELEN,file_ptr)!= NULL)
    {
      /* Increment line-count */
      line++;

      /* Truncate the string to strip off the newline */
      len = string_truncate(input_line,LINELEN);

      /* Separate the field name from the field value */
      string_ptr = strstr(input_line," ");

      /* If this is the ligand sequence, store it */
      if (!strncmp(input_line,"LIGAND_SEQ[",11))
        {
          /* Save the name */
          strcpy(seq_name,string_ptr+1);
        }

      /* If this is the ligand range check it against the stored
         ligand ranges */
      else if ((!strncmp(input_line,"LIG_RANGE[",10) ||
                !strncmp(input_line,"LIG_DUPLIC[",11)) &&
               seq_name[0] != '\0')
        {
          /* Get the description */
          strcpy(desc,string_ptr+1);

          /* Get pointer to the first RasMol definition record */
          rasdef_ptr = first_rasdef_ptr;

          /* Loop through all the definition records */
          while (rasdef_ptr != NULL)
            {
              /* Check if this is a ligand range */
              if (rasdef_ptr->type == LIGAND_RANGE)
                {
                  /* Check the range */
                  if (!strcmp(rasdef_ptr->rasmol_ligand_def,desc))
                    {
                      /* Replace current range/name by sequence
                         name */
                      store_string(&(rasdef_ptr->rasmol_ligand_def),
                                   seq_name);

                      /* Increment ligand count */
                      nligands++;
                    }
                }

              /* Get next definition record */
              rasdef_ptr = rasdef_ptr->next_rasdef_ptr;
            }
        }
    }

  /* Close the rasdef file */
  fclose(file_ptr);

  /* Show number of records read in */
  printf("Number of ligands read in from ligands.dat file:     %8d\n",
         nligands);

  /* Return number of definition entries read in */
  return(nligands);
}
/***********************************************************************

define_rep_ligands  -  Get the ligand residue records to point to a
                       representative sequence record

***********************************************************************/

void define_rep_ligands(struct rasdef *first_rasdef_ptr)
{
  int number;

  struct rasdef *rasdef_ptr, *other_rasdef_ptr, *update_rasdef_ptr;

  /* Initialise variables */
  number = 0;

  /* Get pointer to the first RasMol definition record */
  rasdef_ptr = first_rasdef_ptr;

  /* Loop through all the definition records */
  while (rasdef_ptr != NULL)
    {
      /* Check if this is a ligand range */
      if (rasdef_ptr->type == LIGAND_RANGE)
        {
          /* Increment ligand counter and save number */
          number++;
          rasdef_ptr->number = number;

          /* Get pointer to the first RasMol definition record again */
          other_rasdef_ptr = first_rasdef_ptr;

          /* Loop through all the definition records to find other
             ligand range records having the same ligand sequence */
          while (other_rasdef_ptr != NULL)
            {
              /* Check if this is a ligand range and matches the
                 ligand sequence of our current range */
              if (other_rasdef_ptr->type == LIGAND_RANGE &&
                  other_rasdef_ptr != rasdef_ptr &&
                  !strcmp(other_rasdef_ptr->rasmol_ligand_def,
                          rasdef_ptr->rasmol_ligand_def))
                {
                  /* Get first rasdef record again */
                  update_rasdef_ptr = first_rasdef_ptr;

                  /* Loop through all the ligand residue records to
                     get them to point to the representative ligand
                     sequence record */
                  while (update_rasdef_ptr != NULL)
                    {
                      /* If this is a residue record that points to
                         the duplicate range record, set it to point
                         to the representative ligand sequence record */
                      if (update_rasdef_ptr->type == LIGAND &&
                          update_rasdef_ptr->range_rasdef_ptr
                          == other_rasdef_ptr)
                        {
                          /* Change the pointer */
                          update_rasdef_ptr->range_rasdef_ptr
                            = rasdef_ptr;
                        }

                      /* Get next definition record */
                      update_rasdef_ptr = update_rasdef_ptr->next_rasdef_ptr;
                    }

                  /* Delete the duplicate ligand range record */
                  other_rasdef_ptr->type = DELETED;
                }

              /* Get next definition record */
              other_rasdef_ptr = other_rasdef_ptr->next_rasdef_ptr;
            }
        }

      /* Get next definition record */
      rasdef_ptr = rasdef_ptr->next_rasdef_ptr;
    }
}
/***********************************************************************

create_het_record  -  Create a new het record

***********************************************************************/

struct het *create_het_record(struct het **fst_het_ptr,
                              struct het **lst_het_ptr,
                              char *name,char *desc)
{
  struct het *het_ptr;
  struct het *first_het_ptr, *last_het_ptr;

  /* Initialise het pointers */
  last_het_ptr = *lst_het_ptr;
  first_het_ptr = *fst_het_ptr;

  /* Create a new het entry */
  het_ptr = (struct het*)malloc(sizeof(struct het));
  if (het_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for");
      printf(" struct het\n");
      printf("***        Program resdata terminated with error.\n");
      exit (1);
    }

  /* Store the details for this entry */
  strcpy(het_ptr->name,name);
  store_string(&(het_ptr->desc),desc);

  /* Initialise pointer to next het record */
  het_ptr->next_het_ptr = NULL;

  /* If this is the first entry stored, then update pointer to
     head of linked list */
  if (first_het_ptr == NULL)
    first_het_ptr = het_ptr;

  /* Otherwise, make previous entry point to the current
     one */
  else
    last_het_ptr->next_het_ptr = het_ptr;
                      
  /* Save the current het's pointer */
  last_het_ptr = het_ptr;

  /* Return current pointers */
  *fst_het_ptr = first_het_ptr;
  *lst_het_ptr = last_het_ptr;

  return(het_ptr);
}
/***********************************************************************

get_hets  -  Read the ligands.dat file to pick up the Het Group
             names

***********************************************************************/

int get_hets(char *pdbsum_dir,struct rasdef *first_rasdef_ptr,
             struct parameters params)
{
  char input_line[LINELEN + 1], file_name[FILENAME_LEN];
  char het_name[4], desc[LINELEN + 1], metal[4], seq_name[LINELEN + 1];

  char *string_ptr;

  int len, line, nhets;

  struct het *het_ptr;
  struct rasdef *rasdef_ptr;

  FILE *file_ptr;

  /* Initialise */
  het_name[0] = '\0';
  metal[0] = '\0';
  seq_name[0] = '\0';
  line = 0;
  nhets = 0;

  /* Form the file name of the rasmol.dfn file */
  strcpy(file_name,pdbsum_dir);
  if (params.pdb_type == PDBSUM1)
    strcat(file_name,"newsec/");
  strcat(file_name,"ligands.dat");

  /* If have an uploaded PDB file, the file will be in the current
     directory */
  if (!strcmp(params.pdb_code,"uplo"))
    strcpy(file_name,"ligands.dat");

  /* Open the ligands.dat input file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    return(nhets);

  /* Loop through the file processing each record */
  while (fgets(input_line,LINELEN,file_ptr)!= NULL)
    {
      /* Increment line-count */
      line++;

      /* Truncate the string to strip off the newline */
      len = string_truncate(input_line,LINELEN);

      /* Separate the field name from the field value */
      string_ptr = strstr(input_line," ");
      if (string_ptr == NULL)
        {
          strcat(input_line," ");
          string_ptr = strstr(input_line," ");
        }

      /* If this is the ligand sequence, store it */
      if (!strncmp(input_line,"LIGAND_SEQ[",11))
        {
          /* Save the name */
          strcpy(seq_name,string_ptr+1);
        }

      /* If this is a Het Group name, store it */
      else if (!strncmp(input_line,"HET_NAME[",9))
        {
          /* Save the name */
          strcpy(het_name,string_ptr+1);
        }

      /* If this is a Het Group description, apply it to all ligands
         containing this Het Group */
      else if (!strncmp(input_line,"HET_DESC[",9) &&
               het_name[0] != '\0')
        {
          /* Get the description */
          strcpy(desc,string_ptr+1);

          /* If first character is a space, shift along by one */
          if (desc[0] == ' ')
            strcpy(desc,string_ptr+2);

          /* Get pointer to the first RasMol definition record */
          rasdef_ptr = first_rasdef_ptr;

          /* Loop through all the definition records */
          while (rasdef_ptr != NULL)
            {
              /* Check if this is a ligand */
              if (rasdef_ptr->type == LIGAND_RANGE &&
                  !strcmp(rasdef_ptr->rasmol_ligand_def,seq_name))
                {
                  /* Create and store a Het Group record */
                  het_ptr
                    = create_het_record(&(rasdef_ptr->first_het_ptr),
                                        &(rasdef_ptr->last_het_ptr),
                                        het_name,desc);

                  /* Increment count of Het Groups */
                  rasdef_ptr->nhets++;
                  nhets++;
                }

              /* Get next definition record */
              rasdef_ptr = rasdef_ptr->next_rasdef_ptr;
            }

          /* Re-initialise Het Group name */
          het_name[0] = '\0';
        }

      /* If this is a metal name, store it */
      else if (!strncmp(input_line,"MET_NAME[",9))
        {
          /* Save the name */
          strcpy(metal,string_ptr+1);
        }

      /* If this is a metal description, store it */
      else if (!strncmp(input_line,"MET_DESC[",9) && metal[0] != '\0')
        {
          /* Get the description */
          strcpy(desc,string_ptr+1);

          /* If first character is a space, shift along by one */
          if (desc[0] == ' ')
            strcpy(desc,string_ptr+2);

          /* Get pointer to the first RasMol definition record */
          rasdef_ptr = first_rasdef_ptr;

          /* Loop through all the definition records */
          while (rasdef_ptr != NULL)
            {
              /* Check if this is a metal */
              if (rasdef_ptr->type == METAL &&
                  !strcmp(rasdef_ptr->res_name,metal))
                {
                  /* Store metal name */
                  store_string(&(rasdef_ptr->rasmol_ligand_def),desc);
                }

              /* Get next definition record */
              rasdef_ptr = rasdef_ptr->next_rasdef_ptr;
            }

          /* Re-initialise metal name */
          metal[0] = '\0';
        }
    }

  /* Close the rasdef file */
  fclose(file_ptr);

  /* Show number of records read in */
  printf("Number of Het Groups read in from ligands.dat file:  %8d\n",
         nhets);

  /* Return number of definition entries read in */
  return(nhets);
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
  chain_ptr->copy_of = copy_of;
  chain_ptr->nresid = 0;
  chain_ptr->first_residue_ptr = NULL;
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

int read_chain_equivs(char *pdbsum_dir,struct chain **fst_chain_ptr,
                      struct parameters params)
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
  strcpy(file_name,pdbsum_dir);
  if (params.pdb_type == PDBSUM1)
    strcat(file_name,"newsec/");
  strcat(file_name,"chains.dat");

  /* If have an uploaded PDB file, the file will be in the current directory */
  if (!strcmp(params.pdb_code,"uplo"))
    strcpy(file_name,"chains.dat");

  /* Open the chains.dat input file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    {
      printf("\n*** Warning. Unable to chains.dat file [%s]\n",file_name);
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

          /* Otherwise, chain is unique */
          else
            {
              copy_of = '\0';
              nunique++;
            }

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

  /* Return pointer to first chain record */
  *fst_chain_ptr = first_chain_ptr;

  /* Show number of records read in */
  printf("Number of chains read in from chains.dat file:       %8d\n",
         nchains);
  printf("Number of unique chains:                             %8d\n",
         nunique);

  /* Return the number of chains stored */
  return(nchains);
}
/***********************************************************************

get_atom_details  -  Extract the atom details from the given line

***********************************************************************/

void get_atom_details(char input_line[LINELEN + 1],int *atom_number,
                      char atom_name[5],
                      char *chain,char res_name[4],char res_num[6],
                      double *xval,double *yval,double *zval,
                      double *bval,double *occup)
{
  char number_string[20];

  int len;

  double bvalue, occupancy, x, y, z;

  /* Get the length of the input line */
  len = strlen(input_line);

  /* Get the residue details */
  strncpy(res_name,input_line+17,3);
  res_name[3] = '\0';
  strncpy(res_num,input_line+22,5);
  res_num[5] = '\0';
  *chain = input_line[21];

  /* Get the atom-number */
  strncpy(number_string,input_line+6,5);
  number_string[5] = '\0';
  *atom_number = atoi(number_string);

  /* Get the atom type */
  strncpy(atom_name,input_line+12,4);
  atom_name[4] = '\0';

  /* Retrieve the atomic coordinates, occupancy and B-value */
  strncpy(number_string,input_line+30,8);
  number_string[8] = '\0';
  x = atof (number_string);
  strncpy(number_string,input_line+38,8);
  number_string[8] = '\0';
  y = atof (number_string);
  strncpy(number_string,input_line+46,8);
  number_string[8] = '\0';
  z = atof (number_string);

  /* Retrieve the occupancy and B-value, if they are present */
  if (len > 59)
    {
      strncpy(number_string,input_line+54,6);
      number_string[6] = '\0';
      occupancy = atof(number_string);
    }
  else
    occupancy = 0.0;
  if (len > 65)
    {
      strncpy(number_string,input_line+60,6);
      number_string[6] = '\0';
      bvalue = atof(number_string);
    }
  else
    bvalue = 0.0;

  /* Store the retrieved values */
  *xval = x;
  *yval = y;
  *zval = z;
  *bval = bvalue;
  *occup = occupancy;
}
/***********************************************************************

check_rasdef  -  Check rasdef pointers

***********************************************************************/

void check_rasdef(struct rasdef *first_rasdef_ptr)
{
  int have_exact;

  struct rasdef *rasdef_ptr;

  /* Get pointer to the first RasMol definition record */
  rasdef_ptr = first_rasdef_ptr;

  /* Loop through all the definition records */
  while (rasdef_ptr != NULL)
    {
      if (rasdef_ptr->type == -9)
        {
          printf("GONE WRONG\n");
          exit(1);
        }

      /* Get next definition record */
      rasdef_ptr = rasdef_ptr->next_rasdef_ptr;
    }
}
/***********************************************************************

assign_to_object  -  Check which of the elements in the rasmol.dfn file
                     the given residue should be assigned to

***********************************************************************/

struct rasdef *assign_to_object(struct rasdef *first_rasdef_ptr,
                                char *atom_name,char *res_name,
                                char *res_num,char chain)
{
  int have_exact;

  struct rasdef *rasdef_ptr, *assigned_rasdef_ptr;

  /* Initialise variables */
  assigned_rasdef_ptr = NULL;
  have_exact = FALSE;

  /* Get pointer to the first RasMol definition record */
  rasdef_ptr = first_rasdef_ptr;

  /* Loop through all the definition records */
  while (rasdef_ptr != NULL && have_exact == FALSE)
    {
      /* Process if not a deleted entry */
      if (rasdef_ptr->type != DELETED)
        {
          /* Check for a match to a metal */
          if (rasdef_ptr->type == METAL)
            {
              /* Check that atom name matches that stored in the
                 residue number */
              if (!strncmp(atom_name,rasdef_ptr->res_num,4))
                {
                  /* Store the pointer and set flag */
                  have_exact = TRUE;
                  assigned_rasdef_ptr = rasdef_ptr;
                }
            }

          /* Check for exact match */
          else if (!strcmp(rasdef_ptr->res_name,res_name) &&
                   !strcmp(rasdef_ptr->res_num,res_num) &&
                   rasdef_ptr->chain == chain)
            {
              /* Store the pointer and set flag */
              have_exact = TRUE;
              assigned_rasdef_ptr = rasdef_ptr;
            }

          /* Otherwise, check for a match on the chain-id only */
          else if (rasdef_ptr->chain == chain &&
                   rasdef_ptr->res_name[0] == '\0' &&
                   rasdef_ptr->type != LIGAND_RANGE)
            assigned_rasdef_ptr = rasdef_ptr;
        }

      /* Get next definition record */
      rasdef_ptr = rasdef_ptr->next_rasdef_ptr;
    }

  /* Return the object to which this residue is to be assigned */
  return(assigned_rasdef_ptr);
}
/***********************************************************************

create_residue_record  -  Create and initialise a new residue record

***********************************************************************/

struct residue *create_residue_record(struct residue **fst_residue_ptr,
                                      struct residue **lst_residue_ptr,
                                      char chain,char *res_name,
                                      char *res_num,
                                      struct chain *chain_ptr)
{
  int i, type;

  struct residue *residue_ptr, *first_residue_ptr, *last_residue_ptr;

  /* Initialise residue pointers */
  residue_ptr = NULL;
  first_residue_ptr = *fst_residue_ptr;
  last_residue_ptr = *lst_residue_ptr;

  /* Allocate memory for structure to hold residue info */
  residue_ptr = (struct residue *) malloc(sizeof(struct residue));
  if (residue_ptr == NULL)
    {
      printf("*** Can't allocate memory for struct residue\n");
      exit (1);
    }

  /* If this is the very first, save its pointer */
  if (first_residue_ptr == NULL)
    first_residue_ptr = residue_ptr;

  /* Add link from previous residue to the current one */
  if (last_residue_ptr != NULL)
    last_residue_ptr->next_residue_ptr = residue_ptr;
  last_residue_ptr = residue_ptr;

  /* Get the single-letter amino acid code */
  residue_ptr->aa_code = get_aa_code(res_name);

  /* Initialise the elements of the structure */
  strcpy(residue_ptr->res_name,res_name);
  strcpy(residue_ptr->res_num,res_num);
  residue_ptr->chain = chain;
  residue_ptr->sec_struc = '.';
  residue_ptr->beta_turn = FALSE;
  residue_ptr->gamma_turn = FALSE;
  residue_ptr->disulphide_id = 0;
  for (type = 0; type < CONTACT_TYPES; type++)
    {
      residue_ptr->nhbonds[type] = 0;
      residue_ptr->ncontacts[type] = 0;
    }
  residue_ptr->domain = -1;
  residue_ptr->nsites = 0;
  residue_ptr->npatterns = 0;
  for (i = 0; i < MAX_SITES; i++)
    residue_ptr->site_ptr[i] = NULL;

  /* Store chain, and update chain record */
  residue_ptr->chain_ptr = chain_ptr;
  if (chain_ptr->first_residue_ptr == NULL)
    chain_ptr->first_residue_ptr = residue_ptr;
  chain_ptr->nresid++;
  residue_ptr->first_atom_ptr = NULL;
  residue_ptr->first_rlink_ptr = NULL;

  /* Initialise pointer to the next residue */
  residue_ptr->next_residue_ptr = NULL;

  /* Return current pointers */
  *fst_residue_ptr = first_residue_ptr;
  *lst_residue_ptr = last_residue_ptr;

  /* Return new residue pointer */
  return(residue_ptr);
}
/***********************************************************************

create_atom_record  -  Create a new atom record

***********************************************************************/

struct atom *create_atom_record(struct atom **fst_atom_ptr,
                                struct atom **lst_atom_ptr,
                                int atom_number,char *atom_name,
                                double x,double y,double z,
                                double bvalue,double occupancy,
                                struct residue *residue_ptr)
{
  struct atom *atom_ptr;
  struct atom *first_atom_ptr, *last_atom_ptr;

  /* Initialise atom pointers */
  atom_ptr = NULL;
  last_atom_ptr = *lst_atom_ptr;
  first_atom_ptr = *fst_atom_ptr;

  /* Create atom entry */
  atom_ptr = (struct atom*)malloc(sizeof(struct atom));
  if (atom_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for struct atom\n");
      printf("***        Program resdata terminated with error.\n");
      exit (1);
    }

  /* Store the details for this entry */
  atom_ptr->atom_number = atom_number;
  strcpy(atom_ptr->atom_name,atom_name);
  atom_ptr->coords[0] = x;
  atom_ptr->coords[1] = y;
  atom_ptr->coords[2] = z;
  atom_ptr->bvalue = bvalue;
  atom_ptr->occupancy = occupancy;

  /* If residue's first atom, store pointer to it */
  if (residue_ptr->first_atom_ptr == NULL)
    residue_ptr->first_atom_ptr = atom_ptr;

  /* Initialise pointer to next atom record */
  atom_ptr->next_atom_ptr = NULL;

  /* If this is the first entry stored, then update
     pointer to head of linked list */
  if (first_atom_ptr == NULL)
    first_atom_ptr = atom_ptr;

  /* Otherwise, make previous entry point to the current
     one */
  else
    last_atom_ptr->next_atom_ptr = atom_ptr;
                      
  /* Save the current atom's pointer */
  last_atom_ptr = atom_ptr;

  /* Return current pointers */
  *fst_atom_ptr = first_atom_ptr;
  *lst_atom_ptr = last_atom_ptr;
  return(atom_ptr);
}
/***********************************************************************

scan_sites  -  Determine whether given residue is part of any active
               site

***********************************************************************/

int scan_sites(struct residue *residue_ptr,struct site *first_site_ptr)
{
  int i, nsite_res;
  int done;

  struct site *site_ptr;

  /* Initialise variables */
  done = FALSE;
  nsite_res = 0;

  /* Get pointer to the first site */
  site_ptr = first_site_ptr;

  /* Loop over all the sites */
  while (site_ptr != NULL && done == FALSE)
    {
      /* Check residue name, number and chain id */
      if (!strcmp(site_ptr->res_num,residue_ptr->res_num) &&
          !strcmp(site_ptr->res_name,residue_ptr->res_name) &&
          site_ptr->chain == residue_ptr->chain)
        {
          /* Store the site pointer, if can */
          if (residue_ptr->nsites < MAX_SITES)
            {
              i = residue_ptr->nsites;
              residue_ptr->site_ptr[i] = site_ptr;
              residue_ptr->nsites++;
            }

          /* Increment count of site residues */
          nsite_res++;

          /* Set flag that we're done */
          done = TRUE;
        }

      /* Get the next site */
      site_ptr = site_ptr->next_site_ptr;
    }

  /* Return number of identified site residues */
  return(nsite_res);
}
/***********************************************************************

create_site_record  -  Create and initialise a new site record

***********************************************************************/

struct site *create_site_record(struct site **fst_site_ptr,
                                struct site **lst_site_ptr,char *name,
                                int number,char *res_name,char *res_num,
                                char chain,int type)
{
  struct site *site_ptr, *first_site_ptr, *last_site_ptr;

  /* Initialise site pointers */
  site_ptr = NULL;
  first_site_ptr = *fst_site_ptr;
  last_site_ptr = *lst_site_ptr;

  /* Allocate memory for structure to hold site info */
  site_ptr = (struct site *) malloc(sizeof(struct site));
  if (site_ptr == NULL)
    {
      printf("*** Can't allocate memory for struct site\n");
      exit (1);
    }

  /* If this is the very first site, then store pointer */
  if (first_site_ptr == NULL)
    first_site_ptr = site_ptr;

  /* Add link from previous site to the current one */
  if (last_site_ptr != NULL)
    last_site_ptr->next_site_ptr = site_ptr;
  last_site_ptr = site_ptr;

  /* Store the site details */
  store_string(&(site_ptr->name),name);
  site_ptr->number = number;
  strcpy(site_ptr->res_name,res_name);
  strcpy(site_ptr->res_num,res_num);
  site_ptr->chain = chain;
  site_ptr->type = type;

  /* Initialise remaining fields */
  site_ptr->desc = &null;
  site_ptr->next_site_ptr = NULL;

  /* Return current pointers */
  *fst_site_ptr = first_site_ptr;
  *lst_site_ptr = last_site_ptr;

  /* Return new site pointer */
  return(site_ptr);
}
/***********************************************************************

create_motif_record  -  Create and initialise a new motif record

***********************************************************************/

struct motif *create_motif_record(struct motif **fst_motif_ptr,
                                        struct motif **lst_motif_ptr,
                                        int data_type,char sec_struc,
                                        char *res_num1,char chain1,
                                        char *res_num2,char chain2)
{
  struct motif *motif_ptr, *first_motif_ptr, *last_motif_ptr;

  /* Initialise motif pointers */
  motif_ptr = NULL;
  first_motif_ptr = *fst_motif_ptr;
  last_motif_ptr = *lst_motif_ptr;

  /* Allocate memory for structure to hold motif info */
  motif_ptr = (struct motif *) malloc(sizeof(struct motif));
  if (motif_ptr == NULL)
    {
      printf("*** Can't allocate memory for struct motif\n");
      exit (1);
    }

  /* If this is the very first, save its pointer */
  if (first_motif_ptr == NULL)
    first_motif_ptr = motif_ptr;

  /* Add link from previous motif to the current one */
  if (last_motif_ptr != NULL)
    last_motif_ptr->next_motif_ptr = motif_ptr;
  last_motif_ptr = motif_ptr;

  /* Initialise the elements of the structure */
  motif_ptr->data_type = data_type;
  motif_ptr->sec_struc = sec_struc;
  strcpy(motif_ptr->res_num1,res_num1);
  strcpy(motif_ptr->res_num2,res_num2);
  motif_ptr->chain1 = chain1;
  motif_ptr->chain2 = chain2;

  /* Initialise pointer to the next motif */
  motif_ptr->next_motif_ptr = NULL;

  /* Return current pointers */
  *fst_motif_ptr = first_motif_ptr;
  *lst_motif_ptr = last_motif_ptr;

  /* Return new motif pointer */
  return(motif_ptr);
}
/***********************************************************************

read_pdb_file  -  Read in the given structure from the PDB file

***********************************************************************/

int read_pdb_file(char *pdbsum_dir,struct c_file_data file_name_ptr,
                  struct residue **fst_residue_ptr,
                  struct rasdef *first_rasdef_ptr,
                  struct motif **fst_author_motif_ptr,
                  struct site **fst_site_ptr,
                  struct chain *first_chain_ptr,
                  struct parameters params)
{
  char input_line[LINELEN + 1];
  char file_name[FILENAME_LEN];
  char atom_name[5], chain, res_name[4], res_num[6];
  char res_name1[4], res_name2[4], res_num1[6], res_num2[6];
  char last_chain, last_res_name[6], last_res_num[6], site_name[4];
  char sec_struc;

  int atom_number, ifile, ipos, ires, len, natoms, nresid, nsite_res;
  int dir, done, eof, gzipped, keep_reading, wanted;

  double bvalue, occupancy, x, y, z;

  struct atom *atom_ptr, *last_atom_ptr;
  struct chain *chain_ptr, *match_chain_ptr;
  struct motif *first_motif_ptr, *last_motif_ptr, *motif_ptr;
  struct rasdef *rasdef_ptr;
  struct residue *residue_ptr, *first_residue_ptr, *last_residue_ptr;
  struct site *site_ptr, *first_site_ptr, *last_site_ptr;

#ifdef TEST
  FILE *file_ptr;
#else
  FILE *file_ptr;
  gzFile gzfile_ptr;
#endif

  /* Initialise variables */
  eof = FALSE;
  nresid = 0;
  keep_reading = TRUE;
  last_res_num[0] = '\0';
  last_res_name[0] = '\0';
  last_chain = '\0';
  natoms = 0;
  nsite_res = 0;
  site_name[0] = '\0';

  /* Initialise pointers */
  atom_ptr = NULL;
  last_atom_ptr = NULL;
  residue_ptr = NULL;
  first_residue_ptr = NULL;
  last_residue_ptr = NULL;
  site_ptr = NULL;
  first_site_ptr = NULL;
  last_site_ptr = NULL;
  match_chain_ptr = NULL;
  motif_ptr = NULL;
  first_motif_ptr = NULL;
  last_motif_ptr = NULL;

  /* Form the name of the PDB file */
  if (params.pdb_type == PDBSUM1)
    sprintf(file_name,"%s/%s.new",pdbsum_dir,params.pdb_code);
  else
    dir = r_find_pdb_file(params.pdb_code,file_name_ptr.pdb_template,
                          file_name,params.pdb_type);

  /* If file not found, use the file.pdb in the PDBsum directory */
  if (dir == NOT_FOUND)
    {
      strcpy(file_name,pdbsum_dir);
      strcat(file_name,"/file.pdb");
    }

  /* If have an uploaded PDB file, the file will be in the current
     directory */
  if (!strcmp(params.pdb_code,"uplo"))
    strcpy(file_name,"pdb1xxx.new");

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
      printf("echo \"*** ERROR. Unable to open PDB file ");
      printf("[%s] \"\n",file_name);
      printf("echo \" \"\n");
      exit(-1);
    }

  /* Initialise flag */
  eof = FALSE;

  /* Loop through the PDB file, picking up all the atoms that are
     required */
  while (eof == FALSE && keep_reading == TRUE)
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
          /* Check for the end of an NMR structure model or start of CONECT
             records */
          if (!strncmp("ENDMDL",input_line,6) ||
              !strncmp("CONECT",input_line,6))
            keep_reading = FALSE;

          /* If this is an ATOM or HETATM record, then process */
          else if (!strncmp(input_line,"ATOM  ",6) ||
                   !strncmp(input_line,"HETATM",6))
            {
              /* Get all the atom details */
              get_atom_details(input_line,&atom_number,atom_name,
                               &chain,res_name,res_num,&x,&y,&z,
                               &bvalue,&occupancy);

              /* Check whether this is a new chain */
              if (chain != last_chain)
                {
                  /* New chain, so check whether this chain is wanted */
                  match_chain_ptr = NULL;
                  chain_ptr = first_chain_ptr;

                  /* Loop over the stored protein chain identifiers to see
                     if this is one of the wanted chains */
                  while (chain_ptr != NULL)
                    {
                      /* If chain matches and is a unique chain, then
                         want this chain */
                      if (chain == chain_ptr->chain_id &&
                          chain_ptr->copy_of == '\0')
                        match_chain_ptr = chain_ptr;

                      /* Get the next chain */
                      chain_ptr = chain_ptr->next_chain_ptr;
                    }
                }

              /* Set flag to indicate if this residue is wanted */
              if (match_chain_ptr != NULL)
                {
                  wanted = TRUE;
                  chain_ptr = match_chain_ptr;
                }
              else
                wanted = FALSE;

              /* Ignore certain residue types */
              if (!strncmp(res_name,"ACE",3) || !strncmp(res_name,"HOH",3))
                wanted = FALSE;

              /* If wanted and is a new residue, then store */
              if (wanted == TRUE &&
                  (chain != last_chain || strcmp(res_num,last_res_num) ||
                   strcmp(res_name,last_res_name)))
                {
                  /* Determine whether this residuse is a protein residue */
                  rasdef_ptr
                    = assign_to_object(first_rasdef_ptr,atom_name,res_name,
                                       res_num,chain);

                  /* If this residue is wanted, then save */
                  if (rasdef_ptr != NULL &&
                      (rasdef_ptr->type == PROTEIN ||
                       rasdef_ptr->type == CA_ONLY))
                    {
                      residue_ptr
                        = create_residue_record(&first_residue_ptr,
                                                &last_residue_ptr,chain,
                                                res_name,res_num,chain_ptr);
                      last_atom_ptr = NULL;

                      /* Increment count of residues stored */
                      nresid++;

                      /* Determine whether this residue belongs to an
                         active site (as defined by the previously
                         read in SITE records) */
                      scan_sites(residue_ptr,first_site_ptr);
                    }
                }

              /* If atom wanted, store it */
              if (wanted == TRUE && residue_ptr != NULL)
                {
                  /* Create atom record */
                  atom_ptr
                    = create_atom_record(&(residue_ptr->first_atom_ptr),
                                         &last_atom_ptr,atom_number,
                                         atom_name,x,y,z,bvalue,occupancy,
                                         residue_ptr);
                  natoms++;
                }

              /* Store the current residue details */
              last_chain = chain;
              strcpy(last_res_num,res_num);
              strcpy(last_res_name,res_name);
            }

          /* If this is a SITE identifier, save the code */
          else if (!strncmp("REMARK 800 SITE_IDENTIFIER:",input_line,27))
            {
              /* Save the site name */
              strncpy(site_name,input_line+28,3);
              site_name[3] = '\0';
            }

          /* If this is the SITE description, create a new SITE record */
          else if (!strncmp("REMARK 800 SITE_DESCRIPTION:",input_line,28) &&
                   site_name[0] != '\0')
            {
              /* Initialise unwanted fields */
              nsite_res = 0;
              res_name[0] = res_num[0] = chain = '\0';

              /* Create the record */
              site_ptr
                = create_site_record(&first_site_ptr,&last_site_ptr,site_name,
                                     nsite_res,res_name,res_num,chain,
                                     SITE_DESCRIPTION);

              /* Save the description */
              store_string(&(site_ptr->desc),input_line+29);

              /* Re-initialise site name */
              site_name[0] = '\0';
            }

          /* If this is a SITE record, store the site details */
          else if (!strncmp("SITE  ",input_line,6))
            {
              /* Add blank to end of line to ensure that truncated lines
                 are correctly interpreted */
              strcat(input_line,"                                      ");

              /* Get the site name */
              strncpy(site_name,input_line+11,3);
              site_name[3] = '\0';

              /* Initialise the record position for reading in the residues
                 in this record */
              ipos = 18;
              done = FALSE;

              /* Process each of the potential site residues in turn */
              for (ires = 0; ires < 4 && done == FALSE; ires++)
                {
                  /* Get the chain-ID, residue-name and number */
                  strncpy(res_name,input_line+ipos,3);
                  res_name[3] = '\0';
                  strncpy(res_num,input_line+ipos+5,5);
                  res_num[5] = '\0';
                  chain = input_line[ipos+4];

                  /* If all is blank, then we are done */
                  if (!strcmp(res_name,"   ") || !strcmp(res_num,"     "))
                    done = TRUE;

                  /* Otherwise, create a site record to store these details */
                  else
                    {
                      /* Increment count of site residues */
                      nsite_res++;

                      /* Create the record */
                      site_ptr
                        = create_site_record(&first_site_ptr,
                                             &last_site_ptr,site_name,
                                             nsite_res,res_name,res_num,
                                             chain,RESIDUE_DETAILS);
                    }

                  /* Increment position for next residue */
                  ipos = ipos + 11;
                }
            }

          /* Check for secondary structure definition */
          else if (!strncmp("HELIX ",input_line,6) ||
                   !strncmp("SHEET ",input_line,6))
            {
              /* Read in the helix details */
              if (!strncmp("HELIX ",input_line,6))
                {
                  ifile = 0;
                  sec_struc = 'H';
                  chain = input_line[19];

                  /* Start residue */
                  strncpy(res_name1,input_line+15,3);
                  res_name1[3] = '\0';
                  strncpy(res_num1,input_line+21,5);
                  res_num1[5] = '\0';

                  /* End residue */
                  strncpy(res_name2,input_line+27,3);
                  res_name2[3] = '\0';
                  strncpy(res_num2,input_line+33,5);
                  res_num2[5] = '\0';

                }

              /* Read in the strand details */
              else
                {
                  ifile = 1;
                  sec_struc = 'E';
                  chain = input_line[21];

                  /* Start residue */
                  strncpy(res_name1,input_line+17,3);
                  res_name1[3] = '\0';
                  strncpy(res_num1,input_line+22,5);
                  res_num1[5] = '\0';

                  /* End residue */
                  strncpy(res_name2,input_line+28,3);
                  res_name2[3] = '\0';
                  strncpy(res_num2,input_line+33,5);
                  res_num2[5] = '\0';
                }

              /* Create the motif record */
              motif_ptr
                = create_motif_record(&first_motif_ptr,&last_motif_ptr,
                                      ifile,sec_struc,res_num1,chain,
                                      res_num2,chain);
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

  /* Show number of records read in */
  printf("Number of residues stored:                           %8d\n",
         nresid);
  printf("Number of atoms stored:                              %8d\n",
         natoms);

  /* Return pointer to first residue */
  *fst_residue_ptr = first_residue_ptr;

  /* Return starting motif pointer */
  *fst_author_motif_ptr = first_motif_ptr;

  /* Return first SITE record */
  *fst_site_ptr = first_site_ptr;

  /* Return number of residues read in */
  return(nresid);
}
/***********************************************************************

read_promotif_dat  -  Read in the PROMOTIF analyses from the promotif.dat
                      file, if present

***********************************************************************/

int read_promotif_dat(char *pdb_code,struct motif **fst_promotif_ptr,
                      char *pdbsum_dir,char *chains,int nchains,
                      struct parameters params)
{
#define NMOTIFS 5

  char file_name[FILENAME_LEN], input_line[LINELEN];
  char chain, chain1, chain2, res_num1[6], res_num2[6];
  char res_name1[4], res_name2[4];

  int class, data_end, data_start, i, ichain, ifile, line, motif, nmotif;
  int have_file;

  struct motif *last_promotif_ptr, *first_promotif_ptr;
  struct motif *promotif_ptr;

  FILE *file_ptr;

  static char *motif_name[NMOTIFS] = {
    "HELICES", "STRANDS", "BETA_TURNS", "GAMMA_TURNS", "DISULPHIDES"
  };
  static char sec_struc[]={"HEbgD"};

  /* Initialise pointers */
  first_promotif_ptr = NULL;
  last_promotif_ptr = NULL;

  /* Initialise variables */
  class = 1;
  have_file = TRUE;
  motif = -1;
  nmotif = 0;
  strcpy(res_name1,"XXX");
  strcpy(res_name2,"XXX");

  /* Form name of promotif.dat file */
  strcpy(file_name,pdbsum_dir);
  if (params.pdb_type == PDBSUM1)
    strcat(file_name,"promotif/");
  strcat(file_name,"promotif.dat");

  /* If have an uploaded PDB file, the file will be in the current directory */
  if (!strcmp(pdb_code,"uplo"))
    strcpy(file_name,"promotif.dat");

  /* Open the PROMOTIF file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    have_file = FALSE;

  /* If have the file, read in all its data */
  if (have_file == TRUE)
    {
      /* Loop while reading in records from the file */
      while (fgets(input_line,LINELEN,file_ptr) != NULL)
        {
          /* Truncate the string to strip off the newline */
          string_truncate(input_line,LINELEN);

          /* Check whether this is a START record */
          if (!strncmp(input_line,"START:",6))
            {
              /* Determine which motif this is */
              motif = -1;
              for (i = 0; i < NMOTIFS; i++)
                if (!strcmp(input_line+6,motif_name[i]))
                  motif = i;
            }

          /* Otherwise, if this is a motif record and this is one of
             the motifs we need, store the details */
          else if (motif > -1 && !strncmp(input_line+1,"::",2))
            {
              /* Get the chain that this motif belongs to */
              chain = input_line[0];

              /* Initialise data */
              chain1 = chain;
              chain2 = chain;
              res_num1[0] = '\0';
              res_num2[0] = '\0';

              /* Process according to the different record formats */
              switch (motif)
                {
                  /* a. Helix */
                case 0 :
                  strncpy(res_num1,input_line+8,5);
                  res_num1[5] = '\0';
                  strncpy(res_num2,input_line+20,5);
                  res_num2[5] = '\0';
                  break;

                  /* b. Beta strand */
                case 1 :
                  strncpy(res_num1,input_line+8,5);
                  res_num1[5] = '\0';
                  strncpy(res_num2,input_line+20,5);
                  res_num2[5] = '\0';
                  break;

                  /* c. Beta turns */
                case 2 :
                  strncpy(res_num1,input_line+3,5);
                  res_num1[5] = '\0';
                  strncpy(res_num2,input_line+39,5);
                  res_num2[5] = '\0';
                  break;

                  /* D. Gamma turns */
                case 3 :
                  strncpy(res_num1,input_line+3,5);
                  res_num1[5] = '\0';
                  strncpy(res_num2,input_line+27,5);
                  res_num2[5] = '\0';
                  break;

                  /* E. Disulphide bridge */
                case 4 :
                  chain1 = input_line[3];
                  strncpy(res_num1,input_line+6,5);
                  res_num1[5] = '\0';
                  chain2 = input_line[13];
                  strncpy(res_num2,input_line+16,5);
                  res_num2[5] = '\0';

                  /* Check which of these belong to the same chain */
                  if (chain1 != chain2)
                    {
                      res_num1[0] = '\0';
                      res_num2[0] = '\0';
                    }
                  break;

                  /* F. Beta hairpin */
                case 5 :
                  strncpy(res_num1,input_line+20,5);
                  res_num1[5] = '\0';
                  strncpy(res_num2,input_line+32,5);
                  res_num2[5] = '\0';
                  break;

                default :
                  break;
                }

              /* If have a valid entry, add it to the linked
                 list of this structures PROMOTIF data items */
              if (res_num1[0] != '\0' && res_num2[0] != '\0')
                {
                  /* Create a PROMOTIF record */
                  promotif_ptr =
                    create_motif_record(&first_promotif_ptr,
                                        &last_promotif_ptr,
                                        ifile,sec_struc[motif],
                                        res_num1,chain1,
                                        res_num2,chain2);

                  /* Increment count of motifs stored */
                  nmotif++;
                }
            }
        }

      /* Close the file */
      fclose(file_ptr);
    }

  /* Return pointer to first PROMOTIF entry */
  *fst_promotif_ptr = first_promotif_ptr;

  /* Return number of motifs stored */
  return(nmotif);
}
/***********************************************************************

read_in_promotif_data  -  Read in the PROMOTIF analyses for the given
                          structure, if present

***********************************************************************/

int read_in_promotif_data(char *pdb_code,char *pdbsum_dir,char chain,
                          struct motif **fst_promotif_ptr,
                          struct chain *chain_ptr,
                          struct parameters params)
{
  char chain_str[2], file_name[FILENAME_LEN], input_line[LINELEN];
  char chain1, chain2, res_num1[6], res_num2[6];

  int data_end, data_start, have_file, ifile, line;
  int nmotif, npromotif_files;

  struct motif *first_promotif_ptr, *last_promotif_ptr, *promotif_ptr;

  FILE *file_ptr;

  static char *file_extension[] = {
    "hlx", "str", "btn", "gtn", "dsf" };
  static char sec_struc[]={"HEbgD"};

  /* Initialise variables */
  chain_str[0] = chain;
  chain_str[1] = '\0';
  npromotif_files = strlen(sec_struc);
  nmotif = 0;

  /* Initialise pointers */
  promotif_ptr = NULL;
  first_promotif_ptr = NULL;
  last_promotif_ptr = NULL;

  /* Loop over the different type of PROMOTIF data that may be present */
  for (ifile = 0; ifile < npromotif_files; ifile++)
    {
      /* Initialise flag indicating whether PROMOTIF file present */
      have_file = TRUE;
      line = 0;

      /* Form the current filename */
      strcpy(file_name,pdbsum_dir);
      if (params.pdb_type == PDBSUM1)
        strcat(file_name,"promotif/");
      strcat(file_name,pdb_code);

      /* If have an uploaded PDB file, the file will be in the current
         directory */
      if (!strcmp(pdb_code,"uplo"))
        strcpy(file_name,"1xxx");

      /* Add chain identifier and extension */
      if (chain != ' ')
        strcat(file_name,chain_str);
      strcat(file_name,"_");
      strcat(file_name,file_extension[ifile]);
      strcat(file_name,".html");

      /* Open the PROMOTIF file */
      if ((file_ptr = fopen(file_name,"r")) ==  NULL)
        have_file = FALSE;

      /* If have the file, read in all its data */
      if (have_file == TRUE)
        {
          /* Initialise variables */
          data_start = FALSE;
          data_end = FALSE;

          /* Loop while reading in records from the file */
          while (fgets(input_line,LINELEN,file_ptr) != NULL &&
                 data_end == FALSE)
            {
              /* If waiting for the start of the data, check for dashes */
              if (data_start == FALSE)
                {
                  if (!strncmp(input_line,"----------",10))
                    data_start = TRUE;
                }

              /* Otherwise, get this line of data */
              else
                {
                  /* Check for end of data */
                  if (!strncmp(input_line,"</PRE>",6) ||
                      !strncmp(input_line," <p>",4))
                    data_end = TRUE;
                      
                  /* Read in the data according to the file type */
                  else if ((int) strlen(input_line) > 10)
                    {
                      /* Initialise data */
                      chain1 = chain;
                      chain2 = chain;
                      res_num1[0] = '\0';
                      res_num2[0] = '\0';

                      /* Process according to the different file formats */
                      switch (ifile)
                        {
                          /* a. Helix */
                        case 0 :
                          strncpy(res_num1,input_line+10,5);
                          res_num1[5] = '\0';
                          strncpy(res_num2,input_line+16,5);
                          res_num2[5] = '\0';
                          break;

                          /* b. Beta strand */
                        case 1 :
                          strncpy(res_num1,input_line+11,5);
                          res_num1[5] = '\0';
                          strncpy(res_num2,input_line+17,5);
                          res_num2[5] = '\0';
                          break;

                          /* c. Beta or gamma turns */
                        case 2 :
                        case 3 :
                          strncpy(res_num1,input_line,5);
                          res_num1[5] = '\0';
                          strncpy(res_num2,input_line+6,5);
                          res_num2[5] = '\0';
                          break;

                          /* d. Disulphide bridge */
                        case 4 :
                          chain1 = input_line[5];
                          strncpy(res_num1,input_line+6,5);
                          res_num1[5] = '\0';
                          chain2 = input_line[18];
                          strncpy(res_num2,input_line+19,5);
                          res_num2[5] = '\0';

                          /* Check which of these belong to the
                             same chain */
                          if (chain1 != ' ' && chain1 != chain)
                            res_num1[0] = '\0';
                          if (chain2 != ' ' && chain2 != chain)
                            res_num2[0] = '\0';
                          break;

                          /* e. Beta hairpin */
                        case 5 :
                          strncpy(res_num1,input_line+6,5);
                          res_num1[5] = '\0';
                          strncpy(res_num2,input_line+21,5);
                          res_num2[5] = '\0';
                          break;

                        default :
                          break;
                        }

                      /* If have a valid entry, create a PROMOTIF
                         record */
                      if (res_num1[0] != '\0' && res_num2[0] != '\0')
                        {
                          /* Create the record */
                          promotif_ptr
                            = create_motif_record(&first_promotif_ptr,
                                                  &last_promotif_ptr,
                                                  ifile,sec_struc[ifile],
                                                  res_num1,chain1,
                                                  res_num2,chain2);

                          /* Increment count of records created */
                          nmotif++;
                        }
                    }
                }
            }

          /* Close the file */
          fclose(file_ptr);
        }
    }

  /* Return starting motif pointer */
  *fst_promotif_ptr = first_promotif_ptr;

  /* Return number of motifs stored */
  return(nmotif);
}
/***********************************************************************

motif_annotate  -  Apply the motif annotations to all the residues
                   concerned

***********************************************************************/

void motif_annotate(struct chain *chain_ptr,struct motif *first_motif_ptr)
{
  int iresid, ndisulphide, nresid;
  int in_range;

  struct motif *motif_ptr;
  struct residue *residue_ptr, *next_residue_ptr;

  /* Initialise variables */
  ndisulphide = 0;

  /* Get pointer to the first motif */
  motif_ptr = first_motif_ptr;

  /* Loop over all the motifs */
  while (motif_ptr != NULL)
    {
      /* If this PROMOTIF item corresponds to a disulphide bond, then
         increment count of these */
      if (motif_ptr->sec_struc == 'D')
        ndisulphide++;

      /* Check that this motif is for the correct chain */
      if (motif_ptr->chain1 == chain_ptr->chain_id)
        {
          /* Initialise variables */
          in_range = FALSE;

          /* Get this chain's first residue */
          residue_ptr = chain_ptr->first_residue_ptr;
          nresid = chain_ptr->nresid;
          iresid = 0;

          /* Loop through this chain's residues, to store the MOTIF
             annotations */
          while (residue_ptr != NULL && iresid < nresid)
            {
              /* Check whether this is the start of the residue
                 range */
              if (!strcmp(motif_ptr->res_num1,
                          residue_ptr->res_num) &&
                  motif_ptr->chain1 == residue_ptr->chain)
                in_range = TRUE;

              /* For helix or strand, mark residue if in range */
              if (motif_ptr->sec_struc == 'H' ||
                  motif_ptr->sec_struc == 'E')
                {
                  /* If within the residue range, then update the
                     residue's secondary structure */
                  if (in_range == TRUE)
                    residue_ptr->sec_struc = motif_ptr->sec_struc;

                  /* Check whether this is the end of the
                     residue-range */
                  if (!strcmp(motif_ptr->res_num2,
                              residue_ptr->res_num) &&
                      motif_ptr->chain2 == residue_ptr->chain)
                    in_range = FALSE;
                }

              /* For turns, shift start-position by 1 residue */
              else if (motif_ptr->sec_struc == 'b' ||
                       motif_ptr->sec_struc == 'g')
                {
                  /* If within the residue range, then update the
                     next residue's secondary structure */
                  if (in_range == TRUE)
                    {
                      /* Get the next residue */
                      next_residue_ptr = residue_ptr->next_residue_ptr;
                      if (next_residue_ptr != NULL)
                        {
                          next_residue_ptr->sec_struc = 'T';
                          if (motif_ptr->sec_struc == 'b')
                            next_residue_ptr->beta_turn = TRUE;
                          else
                            next_residue_ptr->gamma_turn = TRUE;
                          in_range = FALSE;
                        }
                    }
                }

              /* For disulphides, check for a residue match */
              else if (motif_ptr->sec_struc == 'D')
                {
                  /* If either end of disulphide bond is matched
                     then store the disulphide number */
                  if ((!strncmp(motif_ptr->res_num1,
                                residue_ptr->res_num,5) &&
                       motif_ptr->chain1 == residue_ptr->chain) ||
                      (!strncmp(motif_ptr->res_num2,
                                residue_ptr->res_num,5) &&
                       motif_ptr->chain1 == residue_ptr->chain))
                    {
                      residue_ptr->disulphide_id = ndisulphide;
                    }
                }

              /* Get the next residue */
              residue_ptr = residue_ptr->next_residue_ptr;
              iresid++;
            }
        }

      /* Get the next motif */
      motif_ptr = motif_ptr->next_motif_ptr;
    }
}
/***********************************************************************

free_motif_memory  -  Free up memory used by promotif records

***********************************************************************/

void free_motif_memory(struct motif *first_motif_ptr)
{
  struct motif *motif_ptr, *next_motif_ptr;

  /* Get first motif */
  motif_ptr = first_motif_ptr;

  /* Loop through the motif records until done */
  while (motif_ptr != NULL)
    {
      /* Get pointer to the next motif record */
      next_motif_ptr = motif_ptr->next_motif_ptr;

      /* Free up memory taken by current motif */
      free(motif_ptr);

      /* Go to the next motif */
      motif_ptr = next_motif_ptr;
    }
}
/***********************************************************************

get_promotif_data  -  Pick up Gail Hutchinson's PROMOTIF assignments from
                      the appropriate WWW directories

***********************************************************************/

void get_promotif_data(char *pdbsum_dir,struct chain *first_chain_ptr,
                       struct motif *first_author_motif_ptr,
                       struct parameters params)
{
  char chains[2];

  int nmotif;

  struct chain *chain_ptr;
  struct motif *first_promotif_ptr;

  /* Get pointer to the first chain */
  chain_ptr = first_chain_ptr;

  /* Loop over all the chains */
  while (chain_ptr != NULL)
    {
      /* If this is a unique chain, get all the PROMOTIF annotations
         for its residues */
      if (chain_ptr->copy_of == '\0')
        {
          /* Get the current chain */
          chains[0] = chain_ptr->chain_id;
          chains[1] = '\0';

          /* Read in all the PROMOTIF annotations from new promotif.dat
             file, if present */
          nmotif
            = read_promotif_dat(params.pdb_code,&first_promotif_ptr,
                                pdbsum_dir,chains,1,params);

          /* If no motifs found, read in all the PROMOTIF annotations for
             the current chain from old-style HTML PROMOTIF files */
          if (nmotif == 0)
            nmotif
              = read_in_promotif_data(params.pdb_code,pdbsum_dir,
                                      chain_ptr->chain_id,
                                      &first_promotif_ptr,chain_ptr,
                                      params);

          /* If have any motifs defined for this chain, identify which
             residues are involved */
          if (nmotif > 0)
            motif_annotate(chain_ptr,first_promotif_ptr);

          /* Otherwise, use any author-defined secondary structure
             assignments read in from the PDB file's header records */
          else
            motif_annotate(chain_ptr,first_author_motif_ptr);

          /* Free memory taken by PROMOTIF pointers */
          free_motif_memory(first_promotif_ptr);
        }

      /* Get the next chain */
      chain_ptr = chain_ptr->next_chain_ptr;
    }
}
/***********************************************************************

get_prosite_strings  -  Get the residue- and colour-strings from the
                        PROSITE pattern line

***********************************************************************/

void get_prosite_strings(char input_line[LINELEN + 1],
                         char residue_string[STRING_LEN],
                         char colours_string[STRING_LEN])
{
  char ch;

  int ipos, ncols, nres, spos;
  int colours, done;

  /* Initialise variables */
  colours = FALSE;
  done = FALSE;
  ipos = 63;
  ncols = 0;
  nres = 0;
  spos = 0;
  while (done == FALSE)
    {
      /* Get the next character from the input line */
      ch = input_line[ipos];

      /* If this is the end of the line, then end here */
      if (ipos > LINELEN || ch == '\0' || ch == '\n')
        done = TRUE;

      /* Otherwise, if non-blank, then store */
      else if (ch != ' ' && ch != '\t' && ch != '?' && ch != '>'
               && ch != '<' && ch != '!' && ch != '|')
        {
          /* Check that string not already too long */
          if (spos < STRING_LEN - 1)
            {
              /* Store one-letter residue code, or colour */
              if (colours == FALSE)
                {
                  residue_string[spos] = ch;
                  nres++;
                }
              else
                {
                  colours_string[spos] = ch;
                  ncols++;
                }

              /* Increment position in output string */
              spos++;
            }
        }

      /* If it is an exclamation mark, need to delete previous
         character as well */
      else if (ch == '!')
        {
          /* If have a previous character, then delete */
          if (spos > 0)
            spos--;
        }

      /* If it is a blank or a TAB, determine what to do about it */
      else if (ch == ' ' || ch == '\t')
        {
          /* If have been reading in the amino acid sequence,
             then flip flag to expect the colours sequence next */
          if (colours == FALSE && nres > 0)
            {
              colours = TRUE;
              spos = 0;
            }
        }

      /* Go to the next input character position */
      ipos++;
    }

  /* Terminate both strings */
  residue_string[nres] = '\0';
  colours_string[ncols] = '\0';
}
/***********************************************************************

check_for_duplicate_pattern  -  Check that we don't already have this
                                PROSITE pattern

***********************************************************************/

int check_for_duplicate_pattern(struct pattern 
                                *first_pattern_ptr,
                                char prosite_id[PROSITE_PATLEN],
                                char chain,char first_res_num[6],
                                char last_res_num[6],
                                char residue_string[STRING_LEN],
                                char colours_string[STRING_LEN])
{
  int wanted;

  struct pattern *pattern_ptr;

  /* Initialise variables */
  wanted = TRUE;

  /* Start at first pattern */
  pattern_ptr = first_pattern_ptr;

  /* Loop through all the patterns */
  while (pattern_ptr != NULL && wanted == TRUE)
    {
      /* Check if have a duplicate */
      if (!strcmp(pattern_ptr->prosite_id,prosite_id) &&
          pattern_ptr->chain == chain &&
          !strcmp(pattern_ptr->first_res_num,first_res_num) &&
          !strcmp(pattern_ptr->last_res_num,last_res_num) &&
          !strcmp(pattern_ptr->residue_string,residue_string) &&
          !strcmp(pattern_ptr->colours_string,colours_string))
        wanted = FALSE;

      /* Go to next pattern */
      pattern_ptr = pattern_ptr->next_pattern_ptr;
    }

  /* Return whether pattern is wanted */
  return(wanted);
}
/***********************************************************************

create_pattern_record  -  Create and initialise a new pattern record

***********************************************************************/

struct pattern *create_pattern_record(struct pattern **fst_pattern_ptr,
                                      struct pattern **lst_pattern_ptr,
                                      char *prosite_id,int occurrence,
                                      int position,char *first_res_num,
                                      char *last_res_num,char chain,
                                      char *residue_string,
                                      char *colours_string)
{
  struct pattern *pattern_ptr, *first_pattern_ptr, *last_pattern_ptr;

  /* Initialise pattern pointers */
  pattern_ptr = NULL;
  first_pattern_ptr = *fst_pattern_ptr;
  last_pattern_ptr = *lst_pattern_ptr;

  /* Allocate memory for structure to hold pattern info */
  pattern_ptr = (struct pattern *) malloc(sizeof(struct pattern));
  if (pattern_ptr == NULL)
    {
      printf("*** Can't allocate memory for struct pattern\n");
      exit (1);
    }

  /* If this is the very first, save its pointer */
  if (first_pattern_ptr == NULL)
    first_pattern_ptr = pattern_ptr;

  /* Add link from previous pattern to the current one */
  if (last_pattern_ptr != NULL)
    last_pattern_ptr->next_pattern_ptr = pattern_ptr;
  last_pattern_ptr = pattern_ptr;

  /* Initialise the elements of the structure */
  strcpy(pattern_ptr->prosite_id,prosite_id);
  pattern_ptr->occurrence = occurrence;
  pattern_ptr->chain = chain;
  pattern_ptr->pattern_no = 0;
  pattern_ptr->true_false = FALSE;

  /* Get start- and end-residue numbers */
  strcpy(pattern_ptr->first_res_num,first_res_num);
  strcpy(pattern_ptr->last_res_num,last_res_num);

  /* Get sequential residue position */
  pattern_ptr->position = position;

  /* Get the strings giving the one-letter residue codes and their
     colours */
  store_string(&(pattern_ptr->residue_string),residue_string);
  store_string(&(pattern_ptr->colours_string),colours_string);

  /* Initialise pointer to the next pattern */
  pattern_ptr->next_pattern_ptr = NULL;

  /* Return current pointers */
  *fst_pattern_ptr = first_pattern_ptr;
  *lst_pattern_ptr = last_pattern_ptr;

  /* Return new pattern pointer */
  return(pattern_ptr);
}
/***********************************************************************

pattern_annotate  -  Apply the pattern annotations to all the residues
                   concerned

***********************************************************************/

void pattern_annotate(struct residue *first_residue_ptr,
                      struct pattern *pattern_ptr)
{
  int i, iresid;
  int done, in_range;

  struct residue *residue_ptr;

  /* Initialise variables */
  done = FALSE;
  in_range = FALSE;
  iresid = 0;

  /* Get the first residue */
  residue_ptr = first_residue_ptr;

  /* Loop through the residues to marke those within the current PROSITE
     pattern */
  while (residue_ptr != NULL && done == FALSE)
    {
      /* Check chain identifier */
      if (pattern_ptr->chain == residue_ptr->chain)
        {
          /* Check whether this is the start of the residue range */
          if (!strcmp(pattern_ptr->first_res_num,residue_ptr->res_num))
            in_range = TRUE;

          /* If in-range, annotate this residue */
          if (in_range == TRUE)
            {
              /* Store the pattern pointer, if can */
              if (residue_ptr->npatterns < MAX_PATTERNS)
                {
                  i = residue_ptr->npatterns;
                  residue_ptr->pattern_ptr[i] = pattern_ptr;
                  residue_ptr->npatterns++;

                  /* Store residue's position in the current pattern */
                  residue_ptr->pattern_pos[i] = iresid;
                }
              
              /* Increment residue position */
              iresid++;
            }

          /* Check whether this is the start of the residue range */
          if (!strcmp(pattern_ptr->last_res_num,residue_ptr->res_num))
            done = TRUE;
        }

      /* Switch in-range flag off */
      else if (in_range == TRUE)
        done = TRUE;

      /* Get the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
    }
}
/***********************************************************************

get_prosite_patterns  -  Pick up all this sequence's PROSITE required
                         patterns

***********************************************************************/

int get_prosite_patterns(char *pdb_code,struct c_file_data file_name_ptr,
                         struct pattern **fst_pattern_ptr,
                         struct residue *first_residue_ptr,
                         struct parameters params)
{
  char file_name[FILENAME_LEN], input_line[LINELEN + 1];
  char chain, last_chain;
  char first_res_num[6], last_res_num[6];
  char last_prosite_id[PROSITE_PATLEN], prosite_id[PROSITE_PATLEN];
  char colours_string[LINELEN], residue_string[LINELEN];

  int npattn, occurrence, position;
  int done, wanted;

  struct pattern *pattern_ptr, *first_pattern_ptr, *last_pattern_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  done = FALSE;
  last_chain = '\0';
  occurrence = 0;
  last_prosite_id[0] = '\0';
  npattn = 0;

  /* Initialise pointers */
  pattern_ptr = NULL;
  first_pattern_ptr = NULL;
  last_pattern_ptr = NULL;

  /* If have an uploaded PDB file, the file will be in the current directory */
  if (!strcmp(pdb_code,"uplo"))
    strcpy(file_name,"3dp.matches.list");
  else if (params.pdb_type == STANDARD_PDB)
    strcpy(file_name,file_name_ptr.other_dir[PROSITE_FILE]);
  else
    strcpy(file_name,file_name_ptr.other_dir[UPLOAD_PROSITE_FILE]);

  /* Open the PROSITE matches file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    {
      /* If not found, show warning */
      printf("*** Warning. PROSITE matches file not found: %s\n",file_name);
      return(0);
    }

  /* Loop while reading in records from the file */
  while (fgets(input_line,LINELEN,file_ptr) != NULL && done == FALSE)
    {
      /* Truncate the string to strip off the newline */
      string_truncate(input_line,LINELEN);

      /* Check whether this is the correct PDB code */
      if (!strncmp(input_line,pdb_code,4) ||
          (!strcmp(pdb_code,"uplo") && !strncmp(input_line,"1xxx",4)))
        {
          /* Get the chain-id */
          chain = input_line[6];

          /* Get this pattern-id */
          strncpy(prosite_id,input_line+10,7);
          prosite_id[7] = '\0';

          /* If this is not the same chain, or the same pattern, as the
             previous pattern for this structure, then re-initialise
             the occurrence count */
          if (chain != last_chain || strncmp(prosite_id,last_prosite_id,7))
            occurrence = 0;

          /* Otherwise, increment the occurrence count */
          occurrence++;

          /* Save the pattern-id and chain */
          strncpy(last_prosite_id,prosite_id,7);
          last_prosite_id[7] = '\0';
          last_chain = chain;

          /* Get start- and end-residue numbers */
          strncpy(first_res_num,input_line+51,5);
          first_res_num[5] = '\0';
          strncpy(last_res_num,input_line+57,5);
          last_res_num[5] = '\0';

          /* Get sequential residue position */
          sscanf(input_line+46,"%d",&position);

          /* Get the strings giving the one-letter residue codes and
             their colours */
          get_prosite_strings(input_line,residue_string,colours_string);

          /* Check that we don't already have this pattern */
          wanted = 
            check_for_duplicate_pattern(first_pattern_ptr,
                                        prosite_id,chain,first_res_num,
                                        last_res_num,residue_string,
                                        colours_string);

          /* If wanted, then store details */
          if (wanted == TRUE)
            {
              /* Create a new pattern record */
              pattern_ptr
                = create_pattern_record(&first_pattern_ptr,
                                        &last_pattern_ptr,prosite_id,
                                        occurrence,position,
                                        first_res_num,last_res_num,chain,
                                        residue_string,colours_string);

              /* Increment count of records created */
              npattn++;

              /* If this is an uploaded PDB file want to show all matching
                 patterns (as can't verify them anyway), so annotate
                 residues now */
              if (!strcmp(pdb_code,"uplo"))
                {
                  /* Mark pattern as wanted */
                  pattern_ptr->true_false = TRUE;

                  /* Annotate all residues in this pattern */
                  pattern_annotate(first_residue_ptr,pattern_ptr);

                  /* Store the pattern number */
                  pattern_ptr->pattern_no = npattn;
                 }
            }
        }

      /* If have gone past the PDB code we're after, then can end here */
      else if (strncmp(input_line,pdb_code,4) > 0)
        done = TRUE;
    }

  /* Close the file */
  fclose(file_ptr);

  /* Show number of unique patterns stored */
  printf("Number of PROSITE patterns read in:                  %8d\n",
         npattn);

  /* Return pointer to first PROSITE pattern */
  *fst_pattern_ptr = first_pattern_ptr;

  /* Return number of patterns stored */
  return(npattn);
}
/***********************************************************************

verify_prosite_patterns  -  Verify the PROSITE patterns just read in
                            using the true/false file

***********************************************************************/

int verify_prosite_patterns(char *pdb_code,struct c_file_data file_name_ptr,
                            struct pattern *first_pattern_ptr,
                            struct residue *first_residue_ptr,
                            struct parameters params)
{
  char file_name[FILENAME_LEN], input_line[LINELEN + 1];
  char chain, true_false;
  char prosite_id[PROSITE_PATLEN];

  int npattn;
  int done, wanted;

  struct pattern *pattern_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  done = FALSE;
  npattn = 0;

  /* Get the name of the appropriate PROSITE true/false file */
  if (params.pdb_type == STANDARD_PDB)
    strcpy(file_name,file_name_ptr.other_dir[TOF_MATCHES_FILE]);
  else
    strcpy(file_name,file_name_ptr.other_dir[UPLOAD_TOF_MATCHES_FILE]);

  /* Open the PROSITE true/false file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    {
      /* If not found, show warning */
      printf("*** Warning. PROSITE truefalse file not found: %s\n",
             file_name);
      return(0);
    }

  /* Loop while reading in records from the file */
  while (fgets(input_line,LINELEN,file_ptr) != NULL && done == FALSE)
    {
      /* Truncate the string to strip off the newline */
      string_truncate(input_line,LINELEN);

      /* Check whether this is the correct PDB code */
      if (!strncmp(input_line,pdb_code,4))
        {
          /* Get the chain-id */
          chain = input_line[6];

          /* Get this pattern-id */
          strncpy(prosite_id,input_line+11,7);
          prosite_id[7] = '\0';

          /* Get true/false identifier */
          true_false = input_line[43];
          if (true_false == 't' || params.pdb_type != STANDARD_PDB)
            true_false = 'T';

          /* Set flag */
          wanted = FALSE;

          /* Get pointer to the first pattern */
          pattern_ptr = first_pattern_ptr;

          /* Loop over all the patterns to find this one */
          while (pattern_ptr != NULL && wanted == FALSE && true_false == 'T')
            {
              /* Check whether this is the right pattern */
              if (!strcmp(prosite_id,pattern_ptr->prosite_id) &&
                  chain == pattern_ptr->chain)
                {
                  /* Mark pattern as wanted */
                  pattern_ptr->true_false = TRUE;
                  wanted = TRUE;

                  /* Annotate all residues in this pattern */
                  pattern_annotate(first_residue_ptr,pattern_ptr);

                  /* Increment count of patterns retained */
                  npattn++;

                  /* Store the pattern number */
                  pattern_ptr->pattern_no = npattn;
                }

              /* Get the next pattern */
              pattern_ptr = pattern_ptr->next_pattern_ptr;
            }
        }

      /* If have gone past the PDB code we're after, then can end here */
      else if (strncmp(input_line,pdb_code,4) > 0)
        done = TRUE;
    }

  /* Close the file */
  fclose(file_ptr);

  /* Show number of unique patterns stored */
  printf("Number of valid PROSITE patterns retained:           %8d\n",
         npattn);

  /* Return number of patterns stored */
  return(npattn);
}
/***********************************************************************

get_domain_token - Get the next space- or tab-delimited token from the
                   given input line, starting at the given position

***********************************************************************/

int get_domain_token(char *line,int *ipos,char *token,int *done)
{
  char ch;

  int i, token_len;
  int have_token;

  /* Initialise variables */
  *done = FALSE;
  have_token = FALSE;
  i = 0;
  token[0] = '\0';

  /* Loop until token-string found or end of line reached */
  while (*done == FALSE && have_token == FALSE && *ipos < LINELEN)
    {
      /* Get the next character in the string */
      ch = line[*ipos];

      /* Check for end of line */
      if (ch == '\0' || ch == '\n')
        *done = TRUE;

      /* If have a space, then take to be the end of the token if have
         already got something in it */
      else if (ch == ' ')
        {
          /* If have a token, then we are done */
          if (token[0] != '\0')
            have_token = TRUE;
        }

      /* Otherwise, add character to token-string provided not too long */
      else
        {
          if (i < TOKEN_LEN)
            {
              token[i] = ch;
              i++;
            }
        }

      /* Go to the next character in the line */
      (*ipos)++;
    }

  /* Terminate the token-string */
  token[i] = '\0';
  token_len = i;
  if (i > 0)
    *done = FALSE;

  /* Return the token */
  return(token_len);
} 
/***********************************************************************

extract_domain_data - Extract the domain data from the free-format line
                      of the CATH domains file and store in a fixed-
                      format array

**********************************************************************/

int extract_domain_data(char *input_line,char *wolf_code,
                        char domain_details[C_MAXDOMCHAIN][16],
                        int *ndomains,int *nfrags,char *error_message,
                        int verbose)
{
  char chain, number_string[20], res_num[6], token_string[TOKEN_LEN + 1];

  int idomain, ifrag, isegment, lpos, ndetails, nsegments;
  int len, start_pos, token_len;
  int done, have_hyphen, next_is_frag, state;

  /* Initialise variables */
  chain = input_line[4];
  error_message[0] = '\0';
  idomain = 0;
  ifrag = 0;
  ndetails = 0;
  *ndomains = 0;
  *nfrags = 0;
  lpos = 0;
  state = WOLF_CODE;
  token_len = 0;
  token_string[0] = '\0';
  done = FALSE;

  /* Loop through line, character by character, until end reached */
  while (done == FALSE)
    {
      /* Get next token string */
      start_pos = lpos;
      token_len = get_domain_token(input_line,&lpos,token_string,&done);

      /* Process the token-string according to what it is */

      /* If token is empty, then skip */
      if (token_len == 0 || (token_len > 1 && token_string[0] == '('))
        continue;

      else if (state == WOLF_CODE)
        {
          /* Store the wolf code */
          strcpy(wolf_code,token_string);

          /* Switch state */
          state = NDOMAINS;
        }

      /* Number of domains */
      else if (state == NDOMAINS)
        {
          /* Get the number of domains */
          if (token_string[0] == 'D')
            {
              /* Extract the number */
              strcpy(number_string,token_string+1);
              *ndomains = atoi(number_string);
            }

          /* If have wrong token, then generate error message */
          else
            {
              sprintf(error_message,"Missing Dxx token. Instead have [%s]",
                      token_string);
              done = TRUE;
            }

          /* Switch state */
          state = NFRAGS;
        }

      /* Number of fragments */
      else if (state == NFRAGS)
        {
          /* Get the number of inter-domain fragments */
          if (token_string[0] == 'F')
            {
              /* Extract the number */
              strcpy(number_string,token_string+1);
              *nfrags = atoi(number_string);
            }

          /* If have wrong token, then generate error message */
          else
            {
              sprintf(error_message,"Missing Fxx token. Instead have [%s]",
                      token_string);
              done = TRUE;
            }

          /* Switch state */
          state = NSEGMENTS;
        }

      /* Number of segments */
      else if (state == NSEGMENTS)
        {
          /* Get the number of segments */
          strcpy(number_string,token_string);
          nsegments = atoi(number_string);

          /* Initialise the details line */
          sprintf(domain_details[ndetails],"D%2d",idomain + 1);
          isegment = 0;

          /* Switch state */
          state = CHAIN1;
        }

      /* Chain identifier */
      else if (state == CHAIN1 || state == CHAIN2)
        {
          /* If too long to be a chain-id, then show error message */
          if (token_len > 1)
            {
              sprintf(error_message,"Expected chain-id is too long [%s]",
                      token_string);
              done = TRUE;
            }

          /* Otherwise, store the chain-id */
          else 
            strcat(domain_details[ndetails],token_string);

          /* Check whether the chain agrees with that at the start of
             the record */
          if (token_string[0] != chain)
            {
              sprintf(error_message,"Chain mismatch [%c] != [%c]",
                      token_string[0],chain);
              done = TRUE;
            }

          /* Switch state */
          if (state == CHAIN1)
            state = RES1;
          else
            state = RES2;
        }

      /* Residue number */
      else if (state == RES1 || state == RES2)
        {
          /* Check whether have a hyphen on the end of the residue number */
          if (token_string[token_len - 1] == '-')
            {
              /* Have hyphen, so strip it off */
              have_hyphen = TRUE;
              token_len = token_len - 1;
              token_string[token_len] = '\0';
            }
          else
            have_hyphen = FALSE;

          /* Format the residue number as a 5-character string */
          format_residue_number(res_num,token_string);

          /* Store the residue number */
          strcat(domain_details[ndetails],res_num);

          /* Switch state */
          if (state == RES1)
            {
              if (have_hyphen == TRUE)
                state = CHAIN2;
              else
                state = HYPHEN1;
            }
          else
            {
              if (have_hyphen == TRUE)
                {
                  strcpy(token_string,"-");
                  token_len = 1;
                  state = DONE_SEG;
                }
              else
                state = HYPHEN2;
            }
        }

      /* Hyphen/insertion-code between residues */
      else if (state == HYPHEN1)
        {
          /* Check that is really a hyphen, then switch state */
          if (token_string[0] == '-')
            state = CHAIN2;

          /* Otherwise, take this to be the insertion code */
          else if (res_num[4] == ' ')
            {
              domain_details[ndetails][8] = token_string[0];
              state = CHAIN2;
            }

          /* Otherwise, if already have an insertion code, then there's
             something wrong */
          else
            {
              sprintf(error_message,
                      "Duplicate insertion code on 1st residue [%s] - [%s]",
                      res_num,token_string);
              done = TRUE;
            }
        }

      /* Hyphen/insertion-code after second residue */
      else if (state == HYPHEN2)
        {
          /* Check that is really a hyphen, then switch state */
          if (token_string[0] == '-')
            state = DONE_SEG;

          /* Otherwise, take this to be the insertion code */
          else if (res_num[4] == ' ')
            {
              domain_details[ndetails][14] = token_string[0];
              state = DONE_SEG;
            }

          /* Otherwise, if already have an insertion code, then there's
             something wrong */
          else
            {
              sprintf(error_message,
                      "Duplicate insertion code on 2nd residue [%s] - [%s]",
                      res_num,token_string);
              done = TRUE;
            }
        }

      /* If details complete, prepare for the next */
      if (state == DONE_SEG)
        {
          /* Check that length exactly 15 characters */
          len = strlen(domain_details[ndetails]);
          if (len == 15)
            {
              /* Determine what's coming next */

              /* Increment segment count */
              isegment++;

              /* If have finished segments for this domain/fragment
                 then expect number of segments next */
              if (isegment >= nsegments)
                {
                  /* Initialise flag */
                  next_is_frag = FALSE;

                  /* If domain just completed, the increment count */
                  if (domain_details[ndetails][0] == 'D')
                    {
                      /* Increment domains count */
                      idomain++;

                      /* If all domains just done, next info will be
                         about the fragments */
                      if (idomain >= *ndomains)
                        next_is_frag = TRUE;

                      /* Otherwise, expecting number of segments making
                         up the domain next */
                      else
                        state = NSEGMENTS;
                    }

                  /* If domain just completed was a fragment, increment
                     count */
                  else if (domain_details[ndetails][0] == 'F')
                    {
                      /* Increment fragment count and set flag */
                      ifrag++;
                      next_is_frag = TRUE;

                      /* If all fragments done, then can stop here */
                      if (ifrag >= *nfrags)
                        {
                          done = TRUE;
                          next_is_frag = FALSE;
                        }
                    }

                  /* For a fragment, initialise the details */
                  if (next_is_frag == TRUE)
                    {
                      /* Start the next details line */
                      if (ndetails < C_MAXDOMCHAIN - 1)
                        sprintf(domain_details[ndetails + 1],"F%2d",ifrag);

                      /* Set state to look for a chain-id */
                      state = CHAIN1;
                      isegment = 0;
                    }
                }

              /* Otherwise, expect chain-id of first residue next */
              else
                {
                  /* Copy the domain details for next residue range */
                  if (ndetails < C_MAXDOMCHAIN - 1)
                    {
                      strncpy(domain_details[ndetails + 1],
                              domain_details[ndetails],3);
                      domain_details[ndetails + 1][3] = '\0';
                    }

                  /* Set state to look for chain-id */
                  state = CHAIN1;
                }

              /* Increment count of details stored */
              if (ndetails < C_MAXDOMCHAIN - 1)
                ndetails++;
              else
                {
                  sprintf(error_message,
                          "Maximum number of domain details exceeded %d",
                          ndetails);
                  done = TRUE;
                }
            }

          /* Otherwise, give an error message */
          else
            {
              sprintf(error_message,
                      "Domain details wrong length [%s] from pstn %d",
                      token_string,start_pos);
              done = TRUE;
            }
        }
    }

  /* Return numbers of domain segments stored */
  return(ndetails);
}  
/***********************************************************************

domain_annotate  -  Apply the domain annotations to all the residues
                    concerned

***********************************************************************/

void domain_annotate(struct residue *first_residue_ptr,char chain,
                     char *first_res_num,char *last_res_num,int domain)
{
  int iresid;
  int done, in_range;

  struct residue *residue_ptr;

  /* Initialise variables */
  done = FALSE;
  in_range = FALSE;
  iresid = 0;

  /* Get the first residue */
  residue_ptr = first_residue_ptr;

  /* Loop through the residues to marke those within the current PROSITE
     domain */
  while (residue_ptr != NULL && done == FALSE)
    {
      /* Check chain identifier */
      if (chain == residue_ptr->chain)
        {
          /* Check whether this is the start of the residue range */
          if (!strcmp(first_res_num,residue_ptr->res_num))
            in_range = TRUE;

          /* If in-range, annotate this residue */
          if (in_range == TRUE)
            {
              /* Store the domain number */
              residue_ptr->domain = domain;
            }

          /* Check whether this is the start of the residue range */
          if (!strcmp(last_res_num,residue_ptr->res_num))
            done = TRUE;
        }

      /* Switch in-range flag off */
      else if (in_range == TRUE)
        done = TRUE;

      /* Get the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
    }
}
/***********************************************************************

read_in_domain_data  -  Read in the domain data for the given PDB file

***********************************************************************/

int read_in_domain_data(char *pdb_code,struct c_file_data file_name_ptr,
                        struct residue *first_residue_ptr)
{
  char file_name[FILENAME_LEN], input_line[LINELEN + 1];
  char chain, first_res_num[6], last_res_num[6];
  char cath_code[C_CATH_LEN], number_string[20], wolf_code[10];
  char domain_details[C_MAXDOMCHAIN][16], error_message[LINELEN + 1];

  int domain, idetail, line, len;
  int ndetails, ndomains, nfrags;
  int verbose;

  FILE *file_ptr;

  /* Initialise variables */
  domain = 0;
  error_message[0] = '\0';
  line = 0;
  ndetails = 0;
  verbose = TRUE;

  /* Form the name of the domains file */
  // strcpy(file_name,file_name_ptr.other_dir[DOMALL_FILE]);
  strcpy(file_name,file_name_ptr.other_dir[CATH_DOMAINS]);

  /* Open the domains file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    {
      /* Show warning message about domains file being missing */
      printf("*** Warning. CATH domains file not found: %s\n",file_name);
      return(0);
    }

  /* Loop while reading in records from the file */
  while (fgets(input_line,LINELEN,file_ptr) != NULL)
    {
      /* Truncate the string to strip off the newline */
      len = string_truncate(input_line,LINELEN);

      /* Increment line-count */
      line++;

      /* Check whether this record applies to our PDB code */
      if (len > 14 && !strncmp(input_line,pdb_code,4))
        {
          /* Get chain identifier */
          chain = input_line[4];
          if (chain == '0')
            chain = 'A';

          /* Extract the domain data from the free-format line and
             store in a fixed-format array for ease of use */
          /* Old version for old format CATH files
          ndetails
            = extract_domain_data(input_line,wolf_code,domain_details,
                                  &ndomains,&nfrags,
                                  error_message,verbose);
          */
          ndetails
            = extract_domain_info(input_line,wolf_code,cath_code,
                                  domain_details);

          /* Increment domain count */
          domain++;

          /* If have an error message, display it */
          if (error_message[0] != '\0')
            {
              /* Print the error message */
              printf("*** ERROR in domains file at line %d for ",line);
              printf("entry [%s]\n",wolf_code);
              printf("***       %s\n",error_message);
            }

          /* Otherwise, annotate the residues belonging to this domain */
          else
            {
              /* Loop through the domain detail entries */
              for (idetail = 0; idetail < ndetails; idetail++)
                {
                  /* Process if not a fragment */
                  if (domain_details[idetail][0] == 'D')
                    {
                      /* Get the start and end residue numbers for this
                         domain fragment */
                      strncpy(first_res_num,domain_details[idetail]+4,5);
                      first_res_num[5] = '\0';
                      strncpy(last_res_num,domain_details[idetail]+10,5);
                      last_res_num[5] = '\0';

                      /* Get the domain number */
                      strncpy(number_string,domain_details[idetail]+1,2);
                      number_string[2] = '\0';
                      //domain = atoi(number_string);
 
                      /* Annotate residues belonging to this domain */
                      domain_annotate(first_residue_ptr,chain,
                                      first_res_num,last_res_num,domain);
                    }
                }
            }
        }
    }

  /* Close the file */
  fclose(file_ptr);

  /* Return number of detail items found */
  return(ndetails);
}
/***********************************************************************

create_grow_record  -  Create and initialise a new grow record

***********************************************************************/

struct grow *create_grow_record(struct grow **fst_grow_ptr,
                                struct grow **lst_grow_ptr)
{
  struct grow *grow_ptr, *first_grow_ptr, *last_grow_ptr;

  /* Initialise grow pointers */
  grow_ptr = NULL;
  first_grow_ptr = *fst_grow_ptr;
  last_grow_ptr = *lst_grow_ptr;

  /* Allocate memory for structure to hold grow info */
  grow_ptr = (struct grow *) malloc(sizeof(struct grow));
  if (grow_ptr == NULL)
    {
      printf("*** Can't allocate memory for struct grow\n");
      exit (1);
    }

  /* If this is the very first, save its pointer */
  if (first_grow_ptr == NULL)
    first_grow_ptr = grow_ptr;

  /* Add link from previous grow to the current one */
  if (last_grow_ptr != NULL)
    last_grow_ptr->next_grow_ptr = grow_ptr;
  last_grow_ptr = grow_ptr;

  /* Initialise the elements of the structure */
  grow_ptr->res_name[0] = '\0';
  grow_ptr->res_num[0] = '\0';
  grow_ptr->chain = ' ';
  grow_ptr->data_type = ' ';
  grow_ptr->other_molecule_type = ' ';
  grow_ptr->rasdef_ptr = NULL;

  /* Initialise pointer to the next grow */
  grow_ptr->next_grow_ptr = NULL;

  /* Return current pointers */
  *fst_grow_ptr = first_grow_ptr;
  *lst_grow_ptr = last_grow_ptr;

  /* Return new grow pointer */
  return(grow_ptr);
}
/***********************************************************************

read_in_grow_data  -  Read in the GROW data on contacts with ligands for
                      the given structure, if present

***********************************************************************/

void read_in_grow_data(char *pdbsum_dir,struct rasdef *first_rasdef_ptr,
                       struct grow **fst_grow_ptr,
                       struct parameters params)
{
  char file_name[FILENAME_LEN], input_line[LINELEN + 1];
  char atom_name[5], chain, res_name[4], res_num[6];

  int i;
  int found, have_file, line;

  struct grow *first_grow_ptr, *last_grow_ptr, *grow_ptr;
  struct rasdef *rasdef_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  have_file = TRUE;
  line = 0;

  /* Initialise pointers */
  grow_ptr = NULL;
  first_grow_ptr = NULL;
  last_grow_ptr = NULL;

  /* Form the name of the grow.out */
  strcpy(file_name,pdbsum_dir);
  if (params.pdb_type == PDBSUM1)
    strcat(file_name,"newsec/");
  strcat(file_name,"grow.out");

  /* If have an uploaded PDB file, the file will be in the current directory */
  if (!strcmp(params.pdb_code,"uplo"))
    strcpy(file_name,"grow.out");

  /* Open the GROW file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    have_file = FALSE;

  /* If have the file, read in all its data */
  if (have_file == TRUE)
    {
      /* Loop while reading in records from the file */
      while (fgets(input_line,LINELEN,file_ptr) != NULL)
        {
          /* Truncate the string to strip off the newline */
          string_truncate(input_line,LINELEN);

          /* If this is not a comment line, process its data */
          if (input_line[0] != '#')
            {
              /* Create a grow record */
              grow_ptr
                = create_grow_record(&first_grow_ptr,&last_grow_ptr);

              /* Get the residue name, number and chain */
              strncpy(grow_ptr->res_name,input_line+14,3);
              grow_ptr->res_name[3] = '\0';
              strncpy(grow_ptr->res_num,input_line+9,5);
              grow_ptr->res_num[5] = '\0';
              grow_ptr->chain = input_line[6];
              grow_ptr->data_type = input_line[0];
              grow_ptr->other_molecule_type = input_line[4];

              /* If this is a contact/H-bond to a ligand or metal,
                 determine which ligand/metal it is */
              if (grow_ptr->other_molecule_type == 'L' ||
                  grow_ptr->other_molecule_type == 'M')
                {
                  /* Get the atom name */
                  strncpy(atom_name,input_line+33,4);
                  atom_name[4] = '\0';

                  /* Get the residue name */
                  strncpy(res_name,input_line+44,3);
                  res_name[3] = '\0';
                  strncpy(res_num,input_line+49,5);
                  res_num[5] = '\0';
                  chain = input_line[54];

                  /* Determine which ligand this is */
                  rasdef_ptr
                    = assign_to_object(first_rasdef_ptr,atom_name,
                                       res_name,res_num,chain);

                  /* Store the ligand reference */
                  grow_ptr->rasdef_ptr = rasdef_ptr;

                  /* If this residue is wanted, then save */
                  if (rasdef_ptr != NULL && (rasdef_ptr->type == LIGAND ||
                                             rasdef_ptr->type == METAL))
                    {
                      /* Loop through the chains this ligand interacts
                         with to see if we already have this one */
                      i = 0;
                      found = FALSE;

                      /* Update the chains this ligand interacts
                         with */
                      while (i < MAX_INT_CHAINS && found == FALSE &&
                             rasdef_ptr->interacts_with[i] != ' ')
                        {
                          /* If have this chain, then set flag */
                          if (rasdef_ptr->interacts_with[i] ==
                              grow_ptr->chain)
                            found = TRUE;

                          /* Increment chain position */
                          i++;
                        }

                      /* If chain not found, add it to the current list */
                      if (found == FALSE && i < MAX_INT_CHAINS)
                        rasdef_ptr->interacts_with[i] = grow_ptr->chain;
                    }
                }
            }
        }

      /* Close the file */
      fclose(file_ptr);
    }

  /* Return first grow pointers */
  *fst_grow_ptr = first_grow_ptr;
}
/***********************************************************************

add_rasdef  -  See if we already have this range record for the current
               residue, and add if we don't

***********************************************************************/

void add_rasdef(struct residue *residue_ptr,struct rasdef *rasdef_ptr,
                int hbond,int contact)
{
  struct rlink *rlink_ptr, *match_rlink_ptr, *last_rlink_ptr;

  /* Initialise rasdef pointers */
  last_rlink_ptr = NULL;
  match_rlink_ptr = NULL;

  /* Get pointer to first entry */
  rlink_ptr = residue_ptr->first_rlink_ptr;

  /* Loop through all these to determine whether we already have this
     entry */
  while (rlink_ptr != NULL)
    {
      /* If have match, store pointer */
      if (rlink_ptr->rasdef_ptr == rasdef_ptr)
        match_rlink_ptr = rlink_ptr;

      /* Store pointer as last record */
      last_rlink_ptr = rlink_ptr;

      /* Get next rasdef entry */
      rlink_ptr = rlink_ptr->next_rlink_ptr;
    }

  /* If match found, then retrieve pointer */
  if (match_rlink_ptr != NULL)
    rlink_ptr = match_rlink_ptr;

  /* Otherwise, create a new rasdef entry */
  else
    {
      rlink_ptr = (struct rlink*) malloc(sizeof(struct rlink));
      if (rlink_ptr == NULL)
        {
          printf("*** ERROR. Unable to allocate memory for");
          printf(" struct rasdef\n");
          printf("***        Program terminated with error.\n");
          exit (1);
        }

      /* Store the details for this entry */
      rlink_ptr->rasdef_ptr = rasdef_ptr;

      /* Initialise remaining fields */
      rlink_ptr->nhbonds = 0;
      rlink_ptr->ncontacts = 0;

      /* Initialise pointer to next rasdef record */
      rlink_ptr->next_rlink_ptr = NULL;

      /* If this is the first entry stored, then update pointer to
         head of linked list */
      if (residue_ptr->first_rlink_ptr == NULL)
        residue_ptr->first_rlink_ptr = rlink_ptr;

      /* Otherwise, make previous entry point to the current
         one */
      else
        last_rlink_ptr->next_rlink_ptr = rlink_ptr;
    }

  /* Update this link's interactions */
  if (hbond == TRUE)
    rlink_ptr->nhbonds++;
  if (contact == TRUE)
    rlink_ptr->ncontacts++;
}
/***********************************************************************

get_grow_data  -  Pick up all the H-bonds and non-bonded contacts in this
                  PDB entry's grow.out file

***********************************************************************/

void get_grow_data(char *pdbsum_dir,struct rasdef *first_rasdef_ptr,
                   struct residue *first_residue_ptr,
                   struct parameters params)
{
  int ninter, type;
  int contact, hbond;
  int interact, wanted;

  struct grow *grow_ptr, *first_grow_ptr;
  struct rasdef *rasdef_ptr, *range_rasdef_ptr;
  struct residue *residue_ptr;

  /* Initialise variables */
  ninter = 0;

  /* Read in the GROW analyses for this file */
  read_in_grow_data(pdbsum_dir,first_rasdef_ptr,&first_grow_ptr,params);

  /* Set pointer to the first residue */
  residue_ptr = first_residue_ptr;

  /* Loop over all the residues to store any interaction data just read
     in from the grow.out file */
  while (residue_ptr != NULL)
    {
      /* Initialise flag */
      interact = FALSE;

      /* Get pointer to the first of the GROW data items */
      grow_ptr = first_grow_ptr;

      /* Loop through the grow items */
      while (grow_ptr != NULL)
        {
          /* If the GROW data item corresponds to the current residue,
             then store the interaction data */
          if (!strcmp(grow_ptr->res_num,residue_ptr->res_num) &&
              grow_ptr->chain == residue_ptr->chain)
            {
              /* Get contact type */
              hbond = contact = FALSE;
              type = NOT_WANTED;
              if (grow_ptr->other_molecule_type == 'M')
                type = TO_METAL;
              else if (grow_ptr->other_molecule_type == 'N')
                type = TO_DNA;
              else if (grow_ptr->other_molecule_type == 'L')
                type = TO_LIGAND;
              else
                type = NOT_WANTED;

              /* If one of the wanted types, increment counts */
              if (type != NOT_WANTED)
                {
                  /* Initialise flag */
                  wanted = FALSE;

                  /* Hydrogen bonds */
                  if (grow_ptr->data_type == 'H')
                    {
                      residue_ptr->nhbonds[type]++;
                      interact = TRUE;
                      wanted = TRUE;
                      hbond = TRUE;
                   }

                  /* Covalent and non-bonded contacts */
                  else if (grow_ptr->data_type == 'N' ||
                           grow_ptr->data_type == 'C' ||
                           grow_ptr->data_type == 'M')
                    {
                      residue_ptr->ncontacts[type]++;
                      interact = TRUE;
                      wanted = TRUE;
                      contact = TRUE;
                    }

                  /* For metal or ligand interaction, determine
                     which metal/ligand this is */
                  if ((type == TO_METAL || type == TO_LIGAND) &&
                      wanted == TRUE)
                    {
                      /* Get the ligand/metal residue */
                      rasdef_ptr = grow_ptr->rasdef_ptr;

                      /* Get the corresponding ligand range record */
                      if (type == TO_METAL)
                        range_rasdef_ptr = rasdef_ptr;
                      else if (rasdef_ptr != NULL)
                        range_rasdef_ptr = rasdef_ptr->range_rasdef_ptr;
                      else
                        range_rasdef_ptr = NULL;

                      /* See if we already have this range record for
                         the current residue, and add if we don't */
                      if (range_rasdef_ptr != NULL)
                        add_rasdef(residue_ptr,range_rasdef_ptr,
                                   hbond,contact);
                    }
                }
            }

          /* Get pointer to the next GROW data item */
          grow_ptr = grow_ptr->next_grow_ptr;
        }

      /* If residue involved in interaction, increment count of
         interacting residues */
      if (interact == TRUE)
        ninter++;

      /* Get the next residue-pointer */
      residue_ptr = residue_ptr->next_residue_ptr;
    }

  /* Show number of residues making interactions */
  printf("Number of residues making interactions (from grow)   %8d\n",
         ninter);
}
/***********************************************************************

write_resdata  -  Write out the residue annotations to the resdata.dat
                  file

***********************************************************************/

void write_resdata(struct chain *first_chain_ptr,
                   struct pattern *first_pattern_ptr,int npattn,
                   struct site *first_site_ptr,
                   struct residue *first_residue_ptr,
                   struct rasdef *first_rasdef_ptr,
                   struct parameters params)
{
  char file_name[FILENAME_LEN];

  static char *subscript[] = { "_D", "_L", "_M" };

  int ctype, isite, ipattn, nresid, position, type;

  struct chain *chain_ptr;
  struct het *het_ptr;
  struct pattern *pattern_ptr;
  struct rasdef *rasdef_ptr;
  struct residue *residue_ptr;
  struct rlink *rlink_ptr;
  struct site *site_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  nresid = 0;

  /* Form the file name for the output resdata.dat file */
  strcpy(file_name,"resdata.dat");

  /* Open the output resdata.dat file */
  file_ptr = open_output_file(file_name,TRUE);

  /* Write out the PDB code */
  fprintf(file_ptr,"PDB %s\n",params.pdb_code);

  /* Get pointer to the first chain */
  chain_ptr = first_chain_ptr;

  /* Loop over all the chains to write out the chain equivalences */
  while (chain_ptr != NULL)
    {
      /* Write out this chain */
      fprintf(file_ptr,"CHAIN %c",chain_ptr->chain_id);

      /* If this is a copy chain, show which chain it equivalences to */
      if (chain_ptr->copy_of != '\0')
        fprintf(file_ptr," = %c",chain_ptr->copy_of);

      /* Close off the line */
      fprintf(file_ptr,"\n");

      /* Get the next chain */
      chain_ptr = chain_ptr->next_chain_ptr;
    }

  /* Get pointer to the first site */
  site_ptr = first_site_ptr;

  /* Loop over all the sites to write out their descriptions */
  while (site_ptr != NULL)
    {
      /* If this is a site description record, write out */
      if (site_ptr->type == SITE_DESCRIPTION)
        fprintf(file_ptr,"SITE %s %s",site_ptr->name,site_ptr->desc);

      /* Get the next site */
      site_ptr = site_ptr->next_site_ptr;
    }

  /* Get pointer to the first RasMol definition record */
  rasdef_ptr = first_rasdef_ptr;

  /* Loop through all the definition records to identify ligands */
  while (rasdef_ptr != NULL)
    {
      /* Get the record type */
      type = rasdef_ptr->type;

      /* Check for a ligand-range or metal record */
      if (type == LIGAND_RANGE || type == METAL)
        {
          /* Write out the ligand sequence */
          if (type == METAL)
            fprintf(file_ptr,"METAL %d %s %s\n",rasdef_ptr->number,
                    rasdef_ptr->res_name,rasdef_ptr->rasmol_ligand_def);
          else
            fprintf(file_ptr,"LIGAND %d %s\n",rasdef_ptr->number,
                    rasdef_ptr->rasmol_ligand_def);

          /* Get the first of the Het Groups */
          het_ptr = rasdef_ptr->first_het_ptr;

          /* Loop over all the Het Groups to write out their
             descriptions */
          while (het_ptr != NULL)
            {
              /* Write the Het Group and its description */
              fprintf(file_ptr," HET_GRP: %s %s\n",het_ptr->name,
                      het_ptr->desc);

              /* Get the next Het Group */
              het_ptr = het_ptr->next_het_ptr;
            }
        }

      /* Get next definition record */
      rasdef_ptr = rasdef_ptr->next_rasdef_ptr;
    }

  /* If have any PROSITE patterns, write these out */
  if (npattn > 0)
    {
      /* Write out the number of patterns */
      fprintf(file_ptr,"NPROSITE %2d\n",npattn);

      /* Initialise count */
      ipattn = 0;

      /* Get pointer to the first pattern */
      pattern_ptr = first_pattern_ptr;

      /* Loop over all the patterns to write out their details */
      while (pattern_ptr != NULL)
        {
          /* If this is a valid pattern, write out */
          if (pattern_ptr->true_false == TRUE)
            {
              /* Increment count */
              ipattn++;

              /* Write out this pattern */
              fprintf(file_ptr," PROSITE %2d. %s %s to %s '%c'\n",
                      ipattn,pattern_ptr->prosite_id,
                      pattern_ptr->first_res_num,
                      pattern_ptr->last_res_num,pattern_ptr->chain);
            }

          /* Get the next pattern */
          pattern_ptr = pattern_ptr->next_pattern_ptr;
        }
    }

  /* Line prior to residue list */
  fprintf(file_ptr,"---------------------------------------\n");

  /* Get pointer to the first residue record */
  residue_ptr = first_residue_ptr;

  /* Loop over the residues to write out the annotations */
  while (residue_ptr != NULL)
    {
      /* Write out this residue */
      fprintf(file_ptr,"%c %s %s %c ",residue_ptr->aa_code,
              residue_ptr->res_name,residue_ptr->res_num,
              residue_ptr->chain);

      /* Write out secondary structure annotation and domain number */
      fprintf(file_ptr,"%c %2d",residue_ptr->sec_struc,
              residue_ptr->domain);

      /* Close record off */
      fprintf(file_ptr,"\n");

      /* If residue belongs to a disulphide bond, write out its id */
      if (residue_ptr->disulphide_id != 0)
        fprintf(file_ptr," DISUL %2d\n",residue_ptr->disulphide_id);

      /* If residue is in a beta turn, write out */
      if (residue_ptr->beta_turn == TRUE)
        fprintf(file_ptr," BTURN\n");

      /* Ditto for gamma turn */
      if (residue_ptr->gamma_turn == TRUE)
        fprintf(file_ptr," GTURN\n");

      /* Write out any active sites to which this residue belongs */
      for (isite = 0; isite < residue_ptr->nsites; isite++)
        {
          /* Get this site record */
          site_ptr = residue_ptr->site_ptr[isite];

          /* Write out the site name */
          fprintf(file_ptr," SITE %s\n",site_ptr->name);
        }

      /* Write out any PROSITE patterns to which this residue belongs */
      for (ipattn = 0; ipattn < residue_ptr->npatterns; ipattn++)
        {
          /* Get this pattern record */
          pattern_ptr = residue_ptr->pattern_ptr[ipattn];
          position = residue_ptr->pattern_pos[ipattn];

          /* Write out the pattern number and residue's colour in
             the pattern */
          fprintf(file_ptr," PROSITE %2d %c\n",pattern_ptr->pattern_no,
                  pattern_ptr->colours_string[position]);
        }

      /* Write out numbers of H-bonds, if any, that this residue makes
         with other molecules (DNA, ligands and metals) */
      for (ctype = 0; ctype < CONTACT_TYPES; ctype++)
        {
          /* If have any H-bonds, write out */
          if (residue_ptr->nhbonds[ctype] > 0)
            fprintf(file_ptr," HB%s %3d\n",subscript[ctype],
                    residue_ptr->nhbonds[ctype]);

          /* Repeat for non-bonded contacts */
          if (residue_ptr->ncontacts[ctype] > 0)
            fprintf(file_ptr," NB%s %3d\n",subscript[ctype],
                    residue_ptr->ncontacts[ctype]);
        }

      /* Get pointer to this residue's first link record */
      rlink_ptr = residue_ptr->first_rlink_ptr;

      /* List the ligands/metals involved in the contacts and H-bonds */
      while (rlink_ptr != NULL)
        {
          /* Get the corresponding rasdef record */
          rasdef_ptr = rlink_ptr->rasdef_ptr;

          /* If metal, write out number */
          if (rasdef_ptr->type == METAL)
            fprintf(file_ptr," MET %d\n",rasdef_ptr->number);

          /* Otherwise, write out ligand number */
          else
            {
              /* Write out the ligand number */
              fprintf(file_ptr," LIG %d\n",rasdef_ptr->number);

              /* Write out numbers of H-bonds and non-bonded contacts
                 between the current residue and this ligand */
              fprintf(file_ptr,"   Counts: %d %d\n",rlink_ptr->nhbonds,
                     rlink_ptr->ncontacts);
            }

          /* Get next definition record */
          rlink_ptr = rlink_ptr->next_rlink_ptr;
        }

      /* Get the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
    }

  /* Close the file */
  fclose(file_ptr);
}
/***********************************************************************

write_domains  -  Write out the coords of each domain or of the whole
                  chain if no domains found

***********************************************************************/

void write_domains(struct chain *first_chain_ptr,
                   struct residue *first_residue_ptr,
                   struct c_file_data file_name_ptr,
                   struct parameters params)
{
  char file_name[FILENAME_LEN], input_line[LINELEN + 1];
  char chain;

  int domain, domain_no, ndomains, start_domain, end_domain;
  int have_domains;

  struct atom *atom_ptr;
  struct chain *chain_ptr;
  struct residue *residue_ptr;

  FILE *file_ptr;

  /* Get pointer to the first chain */
  chain_ptr = first_chain_ptr;

  /* Loop over all the chains to write out the domain coords for each
     unique chain */
  while (chain_ptr != NULL)
    {
      /* If this is a unique chain, break it up into its domains
         (if any) */
      if (chain_ptr->copy_of == '\0')
        {
          /* Set flag */
          end_domain = start_domain = -1;
          ndomains = 1;
          have_domains = FALSE;

          /* Get the first of this chain's residues */
          residue_ptr = chain_ptr->first_residue_ptr;

          /* Loop through this chain's residues to determine how many
             domains we have */
          while (residue_ptr != NULL && residue_ptr->chain_ptr == chain_ptr)
            {
              /* If have a domain annotation, store largest domain
                 number */
              if (start_domain == -1)
                start_domain = residue_ptr->domain;
              if (residue_ptr->domain > end_domain)
                end_domain = residue_ptr->domain;

              if (residue_ptr->domain > 0)
                have_domains = TRUE;

              residue_ptr = residue_ptr->next_residue_ptr;
            }

          /* If no domains defined, take whole chain to be a single
             domain */
          if (have_domains == FALSE)
            start_domain = end_domain = ndomains = 1;

          /* Initialise output domain number */
          domain_no = 0;

          /* Loop over the domains to be written out */
          for (domain = start_domain; domain <= end_domain; domain++)
            {
              /* Define chain identifier */
              chain = chain_ptr->chain_id;
              if (chain == ' ')
                chain = '0';

              /* Define file name */
              domain_no++;
              sprintf(file_name,"dom%c%2.2d.pdb",chain,domain_no);

              /* Open the domain output file */
              file_ptr = open_output_file(file_name,TRUE);

              /* Write out the chain identifier */
              fprintf(file_ptr,"set chain = '%c'\n",chain);

              /* Get the first of this chain's residues again */
              residue_ptr = chain_ptr->first_residue_ptr;

              /* Loop through this chain's residues to determine how many
                 domains we have */
              while (residue_ptr != NULL &&
                     residue_ptr->chain_ptr == chain_ptr)
                {
                  /* If this is the correct domain, or we have no domains
                     defined, write out this residue's atom records */
                  if (residue_ptr->domain == domain ||
                      have_domains == FALSE)
                    {
                      /* Get this residue's first atom */
                      atom_ptr = residue_ptr->first_atom_ptr;

                      /* Loop while writing out this residue's atoms */
                      while (atom_ptr != NULL)
                        {
                          /* Write atom to file */
                          fprintf(file_ptr,
                                  "ATOM  %5d %s %s %c%s   %8.3f%8.3f%8.3f",
                                  atom_ptr->atom_number,atom_ptr->atom_name,
                                  residue_ptr->res_name,residue_ptr->chain,
                                  residue_ptr->res_num,atom_ptr->coords[0],
                                  atom_ptr->coords[1],atom_ptr->coords[2]);
                          fprintf(file_ptr,"%6.2f%6.2f\n",atom_ptr->occupancy,
                                  atom_ptr->bvalue);

                          /* Get the next atom */
                          atom_ptr = atom_ptr->next_atom_ptr;
                        }
                    }

                  /* Get the next residue */
                  residue_ptr = residue_ptr->next_residue_ptr;
                }

              /* Close the output domains file */
              fclose(file_ptr);
            }
        }

      /* Get the next chain */
      chain_ptr = chain_ptr->next_chain_ptr;
    }
}

/***********************************************************************

                             M   A   I   N

***********************************************************************/

int main(int argc,char *argv[])
{
  char pdbsum_dir[FILENAME_LEN + 1];

  int nchains, ndef, ndomains, npattn, nresid;
  int dir_type;

  struct c_file_data file_name_ptr;

  struct parameters params;

  struct chain *first_chain_ptr;
  struct motif *first_author_motif_ptr;
  struct pattern *first_pattern_ptr;
  struct rasdef *first_rasdef_ptr;
  struct residue *first_residue_ptr;
  struct site *first_site_ptr;

  /* Initialise variables */
  ndomains = npattn = 0;

  /* Pick up the PDB code from the command-line arguments */
  get_command_arguments(argv,argc - 1,&params);

  /* Get required paths and file names */
  c_get_paths(&file_name_ptr,params.run_64);

  /* Form the name of the PDBsum directory */
  dir_type
    = find_pdbsum_dir(params.pdb_code,file_name_ptr.pdbsum_template,
                      params.pdb_type,pdbsum_dir);
  if (dir_type == NOT_FOUND)
    {
      printf("*** ERROR. Unable to locate PDBsum directory for %s\n",
             params.pdb_code);
      exit(-1);
    }

  /* Read in the definitions from the rasmol.dfn file */
  ndef = get_dfn_data(pdbsum_dir,&first_rasdef_ptr,params);

  /* Read in the ligands.dat file to pick up ligand names */
  get_ligands(pdbsum_dir,first_rasdef_ptr,params);

  /* Define a single residue range to represent each unique ligand
     sequence and get all ligand residue records to point to it */
  define_rep_ligands(first_rasdef_ptr);

  /* Read in the ligands.dat file to pick up Het Group names */
  get_hets(pdbsum_dir,first_rasdef_ptr,params);

  /* Read in all the chain-equivalences from the chains.dat file */
  nchains
    = read_chain_equivs(pdbsum_dir,&first_chain_ptr,params);

  /* Read in the protein sequence(s) from the corresponding PDB file */
  nresid
    = read_pdb_file(pdbsum_dir,file_name_ptr,&first_residue_ptr,
                    first_rasdef_ptr,&first_author_motif_ptr,
                    &first_site_ptr,first_chain_ptr,params);

  /* Read in the PROMOTIF data and annotate the residues accordingly */
  get_promotif_data(pdbsum_dir,first_chain_ptr,first_author_motif_ptr,
                    params);

  /* Read in any PROSITE data for this entry */
  if (params.pdb_type != ALPHA_FOLD_PDB && params.pdb_type != PDBSUM1)
    npattn
      = get_prosite_patterns(params.pdb_code,file_name_ptr,
                             &first_pattern_ptr,first_residue_ptr,
                             params);

  /* Verify the PROSITE patterns just read in (unless have an uploaded
     pdb file) */
  if (strcmp(params.pdb_code,"uplo") &&
      params.pdb_type != ALPHA_FOLD_PDB &&
      params.pdb_type != PDBSUM1)
    npattn
      = verify_prosite_patterns(params.pdb_code,file_name_ptr,
                                first_pattern_ptr,first_residue_ptr,
                                params);

  /* Read in the CATH domain data and annotate residues */
  if (strcmp(params.pdb_code,"uplo") &&
      params.pdb_type != ALPHA_FOLD_PDB &&
      params.pdb_type != PDBSUM1)
    ndomains
      = read_in_domain_data(params.pdb_code,file_name_ptr,
                            first_residue_ptr);

  /* Read in the residue interactions from the grow.out file */
  get_grow_data(pdbsum_dir,first_rasdef_ptr,first_residue_ptr,
                params);

  /* Write out the residue annotations to the resdata.dat file */
  write_resdata(first_chain_ptr,first_pattern_ptr,npattn,
                first_site_ptr,first_residue_ptr,first_rasdef_ptr,
                params);

  /* Write out coords of each of the different domains present in the
     structure */
  write_domains(first_chain_ptr,first_residue_ptr,file_name_ptr,
                params);

  return(0);
}
