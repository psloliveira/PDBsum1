/***********************************************************************

  fixpdb.c - Program to fix common errors in PDB file (missing chain
             ids,incorrect atom and residue naming, etc)

************************************************************************

Date:-         19 Apr 2022

Written by:-   Roman Laskowski

The fixes:

     1. Replace blank chain id by next available
     2. Change zero occupancy to 1.0
     3. Fix non-standard water names (SOL, WAT, H2O, etc) and atom names
     4. Fix duplicate ligand atom names to make them unique
     5. Change DNA base names to standard PDB format
     6. Make all residue names upper case
     7. Remove all ATOM and HETATM records from models > 1.
     8. Set B-value to maximum of 999.99 if over

----------------------------------------------------------------------*/

/* Datafiles
   ---------

Actual name      Description
-----------      -----------
file-name        Input PDB file
file-name.fix    Output PDB file with all fixes

Compilation
-----------

cc -c fixpdb.c
cc -c common.c
cc -c reapar.c
cc -c readpdb.c
cc -o fixpdb fixpdb.o common.o reapar.o readpdb.o -lm -lz

Optimization flags:-  

cc -O2 -funroll-loops -c fixpdb.c
cc -O2 -funroll-loops -c common.c
cc -O2 -funroll-loops -c reapar.c
cc -O2 -funroll-loops -c readpdb.c
cc -o fixpdb fixpdb.o common.o reapar.o readpdb.o -lm -lz

*/

/* Function calling-tree
   ---------------------

main
   -> get_command_arguments
   -> c_get_paths (in common.c)
         -> get_param (in reapar.c)
         -> store_string (in common.c)
   -> check_chains
         -> store_string (in common.c)
         -> append_string (in common.c)
   -> rewrite_pdb_file
         -> string_chomp (in common.c)
         -> r_get_atom_details (in readpdb.c)
         -> to_upper (in common.c)

*/

#include "common.h"
#include "fixpdb.h"

/* Prototypes
   ---------- */

/* In common.c */
void append_string(char **string_ptr,char *string);
void c_get_paths(struct c_file_data *file_name_ptr,int run_64);
void free_tokens(char **token,int ntokens);
int get_tokens(char *line,char **token,int max_token,char token_sep);
FILE *open_output_file(char *out_name,int terminate);
void store_string(char **string_ptr,char *string);
int string_chomp(char *string);
void to_upper(char *search_string,int length);

/* Function prototypes from reapar.c */
int get_param(char *token_name,char *file_name,int filename_len,
              int exit_on_error);

/* Function prototypes from readpdb.c */
int r_get_atom_details(char *input_line,int *atom_number,char *atom_name,
                       char *chain,char *res_name,char *res_num,
                       double *xval,double *yval,double *zval,
                       double *bval,double *occup,char *element,
                       int *rec_type,char wanted_chain);

/***********************************************************************

get_command_arguments  -  Interpret the command line arguments

***********************************************************************/

void get_command_arguments(char *strptrs[],int num,
                           struct parameters *params)
{
  int itoken;

  /* Initialise variables */
  params->file_name = &null;
  params->run_64 = TRUE;

  /* If nothing entered on the command line, show options available */
  if (num == 0 || (num == 1 && !strcmp(strptrs[1],"-help")))
    {
      printf("\n");
      printf("Correct usage is:\n");
      printf("\n");
      printf("    fixpdb  filename\n");
      printf("\n");
      printf("where\n");
      printf("   * filename  = Name of PDB file to be fixed\n");
      printf("   * -help     = Show this help message\n");
      printf("\n");
      printf("\n");
      exit(1);
    }

  /* Get the file name from the first argument */
  store_string(&(params->file_name),strptrs[1]);
}
/***********************************************************************

check_chains  -  Read in the given structure from the PDB file

***********************************************************************/

char check_chains(char *file_name,char **chn_type)
{
  char input_line[LINELEN + 1];
  char chain, last_chain, new_chain;
  char res_name[4], last_res_name[4];
  char ctype[2], curr_chain_type;

  char *chain_list, *chain_type, *string_ptr;

  int blank_chain, chain_list_len, i, base, len, nchains;
  int done;

  double bvalue, occupancy, x, y, z;

  FILE *file_ptr;

  /* Initialise variables */
  blank_chain = -1;
  done = FALSE;
  last_chain = '\0';
  last_res_name[0] = '\0';
  nchains = 0;
  new_chain = ' ';

  /* Prepare the chain list */
  store_string(&chain_list,CHAIN_LIST);
  chain_list_len = strlen(chain_list);

  /* Create array of chain types */
  chain_type = &null;
  curr_chain_type = UNKNOWN;
  
  /* Open PDB file */
  file_ptr = fopen(file_name,"r");
  if (file_ptr ==  NULL)
    {
      printf("*** ERROR. Unable to open PDB file [%s] \"\n",file_name);
      exit(-1);
    }

  /* Loop through the PDB file, picking up all the atoms */
  while (fgets(input_line,LINELEN,file_ptr) != NULL && done == FALSE)
    {
      /* Check for the end of an NMR structure model or start of CONECT
         records */
      if (!strncmp("ENDMDL",input_line,6) ||
          !strncmp("CONECT",input_line,6))
        done = TRUE;

      /* If this is an ATOM record, then get chain */
      else if (!strncmp(input_line,"ATOM  ",6) ||
               !strncmp(input_line,"HETATM",6))
        {
          /* Get the chain identifier */
          chain = input_line[21];

          /* Check whether this is a new chain */
          if (chain != last_chain)
            {
              /* If this chain is blank, then set flag */
              if (chain == ' ' && blank_chain == -1)
                blank_chain = nchains;

              /* Otherwise, mark this chain id as in use */
              else
                {
                  /* Find this chain in the chain list */
                  string_ptr = strchr(chain_list,chain);

                  /* If found, mark this chain as in use */
                  if (string_ptr != NULL)
                    string_ptr[0] = ' ';
                }

              /* Store the current residue chain */
              last_chain = chain;

              /* Initialise residue name */
              last_res_name[0] = '\0';

              /* Increment chain count */
              nchains++;

              /* Append type to list of chain types */
              ctype[0] = UNKNOWN;
              ctype[1] = '\0';
              curr_chain_type = UNKNOWN;
              if (chain_type == &null)
                store_string(&chain_type,ctype);
              else
                append_string(&chain_type,ctype);
            }

          /*Get the residue name */
          strncpy(res_name,input_line+17,3);
          res_name[3] = '\0';

          /* If a new name, what residue type is */
          if (strcmp(res_name,last_res_name))
            {
              /* If not yet assigned see if it is a protein */
              if (curr_chain_type == UNKNOWN)
                {
                  /* Check if this is an amino acid */
                  for (i = 0;
                       i < NAMINO && curr_chain_type == UNKNOWN;i++)
                    {
                      /* If match, save type as protein */
                      if (!strcmp(res_name,AMINO_CODE[i]))
                        curr_chain_type = PROTEIN;
                    }
                }
              
              /* If still unassigned, check for nucleic acid */
              if (curr_chain_type == UNKNOWN)
                {
                  /* Check for possible base names */
                  string_ptr = strstr(BASE_NAMES,res_name);

                  /* If found, mark as nucleic acid */
                  if (string_ptr != NULL)
                    curr_chain_type = NUCLEIC_ACID;
                }

              /* If marked as nucleic acid, check for uracil or
                 thymine */
              if (curr_chain_type == NUCLEIC_ACID)
                {
                  /* Check base names */
                  string_ptr = strstr(BASE_NAMES,res_name);

                  /* If found, see if either uracil or thymine */
                  if (string_ptr != NULL)
                    {
                      /* Determine match location */
                      base = ((string_ptr - BASE_NAMES) / 4) % 5;

                      /* If thymine, then we have DNA */
                      if (base == 3)
                        curr_chain_type = DNA;
                      else if (base == 4)
                        curr_chain_type = RNA;
                    }
                }

              /* Save the current chain type */
              if (chain_type[nchains - 1] == UNKNOWN ||
                  chain_type[nchains - 1] == NUCLEIC_ACID)
                chain_type[nchains - 1] = curr_chain_type;

              /* Save the current residue name */
              strcpy(last_res_name,res_name);
            }
        }
    }

  /* Close PDB file */
  fclose(file_ptr);

  /* If we have a blank chain, decide what to rename it to */
  if (blank_chain > -1 && chain_type[blank_chain] != UNKNOWN)
    {
      /* Initialise flag */
      done = FALSE;

      /* Loop through all the chains, to find next unused one */
      for (i = 0; i < chain_list_len && done == FALSE; i++)
        {
          /* If chain id not been blanked out, can use it */
          if (chain_list[i] != ' ')
            {
              /* Save this chain */
              new_chain = chain_list[i];

              /* Set flag */
              done = TRUE;
            }
        }
    }

  /* Return the chain types */
  *chn_type = chain_type;

  /* Return the new name for any blank chains */
  return(new_chain);
}
/***********************************************************************

rewrite_pdb_file  -  Rewrite the PDB file, making all the necessary
                     fixes

***********************************************************************/

int rewrite_pdb_file(char *file_name,char new_chain,char *chain_type)
{
  char out_name[FILENAME_LEN], input_line[LINELEN + 1];
  char atom_name[5], atom_type[7], chain, res_name[4], res_num[6];
  char res_name1[4], res_name2[4], res_num1[6], res_num2[6];
  char last_chain, last_res_name[6], last_res_num[6], use_res_name[6];
  char element[3], number_string[20], wanted_chain;
  char ch, curr_chain_type;

  char *atom_names, *string_ptr;

  int atom_number, base, i, len, nchains, number, rec_type;
  int natnames, natoms, nchange_chain, nmodels, nnucleic, nwaters;
  int done, last_ter, water;

  double bvalue, occupancy, x, y, z;

  FILE *file_ptr, *outfile_ptr;

  /* Initialise variables */
  atom_names = &null;
  last_chain = '\0';
  last_res_num[0] = '\0';
  last_res_name[0] = '\0';
  use_res_name[0] = '\0';
  last_ter = FALSE;
  natoms = nchains = nmodels = 0;
  natnames = nchange_chain = nnucleic = nwaters = 0;
  wanted_chain = '\0';

  /* Open PDB file */
  file_ptr = fopen(file_name,"r");
  if (file_ptr ==  NULL)
    {
      printf("*** ERROR. Unable to open PDB file [%s] \"\n",file_name);
      exit(-1);
    }

  /* Form name of the output file */
  strcpy(out_name,file_name);
  string_ptr = strchr(out_name,'.');
  if (string_ptr != NULL)
    string_ptr[0] = '\0';
  strcat(out_name,".fix");

  /* Open the output file */
  outfile_ptr = open_output_file(out_name,TRUE);

  /* Loop through the PDB file, picking up all the atoms */
  while (fgets(input_line,LINELEN,file_ptr) != NULL)
    {
      /* Truncate the string to strip off the newline */
      string_chomp(input_line);

      /* Check for the end of an ensemble or NMR structure model */
      if (!strncmp("ENDMDL",input_line,6))
        {
          /* Re-initialise for next model */
          atom_names = &null;
          last_chain = '\0';
          last_res_num[0] = '\0';
          last_res_name[0] = '\0';
          use_res_name[0] = '\0';
          nchains = 0;

          /* If model had any atoms, then increment count */
          if (natoms > 0)
            nmodels++;
        }

      /* If this is an ATOM or HETATM record, then process it */
      if (!strncmp(input_line,"ATOM  ",6) ||
          !strncmp(input_line,"HETATM",6))
        {
          if (!strncmp(input_line,"ATOM  ",6))
            strcpy(atom_type,"ATOM  ");
          else
            strcpy(atom_type,"HETATM");

          /* Get the chain identifier */
          r_get_atom_details(input_line,&atom_number,atom_name,
                             &chain,res_name,res_num,&x,&y,&z,
                             &bvalue,&occupancy,element,&rec_type,
                             wanted_chain);

          /* Initialise flags */
          water = FALSE;

          /* Check whether this is a new chain */
          if (chain != last_chain)
            {
              /* Store the current chain */
              last_chain = chain;

              /* Get the current chain type */
              curr_chain_type = chain_type[nchains];

              /* Increment chain count */
              nchains++;

              /* Re-initialise residue name */
              last_res_name[0] = '\0';
              last_res_num[0] = '\0';
              use_res_name[0] = '\0';
            }

          /* If chain is blank, replace by new chain */
          if (chain == ' ')
            {
              chain = new_chain;
              nchange_chain++;
            }

          /* If a new residue name, process it */
          if (strcmp(res_name,last_res_name))
            {
              /* Take the input name as the output name */
              strcpy(use_res_name,res_name);

              /* Convert to upper case */
              to_upper(use_res_name,3);

              /* If a Het Group, then check for water */
              if (!strcmp(atom_type,"HETATM"))
                {
                  /* Check for one of the non-standard water names */
                  for (i = 0; i < NWATER_NAMES && water == FALSE; i++)
                    {
                      /* If this is a water name, then replace */
                      if (!strcmp(res_name,WATER_NAME[i]))
                        {
                          /* Replace residue name */
                          strcpy(use_res_name,"HOH");

                          /* Replace atom name */
                          strcpy(atom_name," O  ");

                          /* Set flag that we have water */
                          water = TRUE;

                          /* Increment count of water names fixed */
                          nwaters++;
                        }
                    }
                }

              /* If not Het Group, check for DNA/RNA base */
              else if (curr_chain_type == NUCLEIC_ACID ||
                       curr_chain_type == DNA || curr_chain_type == RNA)
                {
                  /* Determine which base this is */
                  string_ptr = strstr(BASE_NAMES,res_name);

                  /* If found, determine which base this is */
                  if (string_ptr != NULL)
                    {
                      /* Get the base */
                      base = ((string_ptr - BASE_NAMES) / 4) % 5;

                      /* Replace the residue name */
                      if (curr_chain_type == RNA)
                        strcpy(use_res_name,RNA_NAME[base]);
                      else
                        strcpy(use_res_name,DNA_NAME[base]);

                      /* If change made, then increment count */
                      if (strcmp(use_res_name,res_name))
                        nnucleic++;
                    }
                }

              /* Save the current residue name */
              strcpy(last_res_name,res_name);
            }

          /* If a new residue number, initialise atom names */
          if (strcmp(res_num,last_res_num))
            {
              /* If already have a atom names string, then free */
              if (atom_names != &null)
                free(atom_names);

              /* Re-initialise */
              atom_names = &null;

              /* Save the current residue number */
              strcpy(last_res_num,res_num);
            }

          /* If a Het Group, then check for duplicate atom names */
          if (!strcmp(atom_type,"HETATM") &&
              strcmp(use_res_name,"HOH"))
            {
              /* Check if we have already seen this atom
                 name in this residue */
              string_ptr = strstr(atom_names,atom_name);

              /* If we have a duplicate, then need to rename
                 the current atom */
              if (string_ptr != NULL)
                {
                  /* Initialise trial number */
                  done = FALSE;
                  number = 0;

                  /* Loop until a suitable number has been generated */
                  while (number < 99 && done == FALSE)
                    {
                      /* Increment the number */
                      number++;

                      /* Format it */
                      sprintf(number_string,"%d ",number);
                      number_string[2] = '\0';

                      /* Add to atom number */
                      strcpy(atom_name + 2,number_string);

                      /* See if we have this name */
                      string_ptr = strstr(atom_names,atom_name);

                      /* If not, then we're done */
                      if (string_ptr == NULL)
                        done = TRUE;
                    }

                  /* If atom name changed, increment count */
                  if (done == TRUE)
                    natnames++;
                }

              /* Append atom name to string */
              if (atom_names == &null)
                store_string(&atom_names,"|");
              else
                append_string(&atom_names,"|");
              append_string(&atom_names,atom_name);
            }

          /* If occupancy is zero, set to 1.0 */
          if (occupancy == 0)
            occupancy = 1.0;

          /* Write out atom, unless belongs to model > 1 */
          if (nmodels == 0)
            {
              /* Write this line out */
              fprintf(outfile_ptr,"%s%5d %s %s %c%s   %8.3f%8.3f%8.3f",
                      atom_type,atom_number,atom_name,use_res_name,chain,
                      res_num,x,y,z);
              fprintf(outfile_ptr,"%6.2f%6.2f          %s\n",occupancy,
                      bvalue,element);

              /* Increment count of atoms written out */
              natoms++;
            }
        }


      /* Otherwise, write straight out */
      else
        fprintf(outfile_ptr,"%s\n",input_line);
    }

  /* Close both files */
  fclose(file_ptr);
  fclose(outfile_ptr);

  /* Show any changes made */
  if (new_chain != ' ')
    printf("Blank chain changed to [%c] for %d atoms\n",new_chain,
           nchange_chain);
  if (natnames > 0)
    printf("%8d atom names corrected\n",natnames);
  if (nwaters > 0)
    printf("%8d water names corrected\n",nwaters);
  if (nnucleic > 0)
    printf("%8d DNA/RNA names corrected\n",nnucleic);

  /* Return number of atoms */
  return(natoms);
}

/***********************************************************************

                             M   A   I   N

***********************************************************************/

int main(int argc,char *argv[])
{
  char new_chain;

  char *chain_type;

  struct c_file_data file_name_ptr;
  struct parameters params;

  /* Initialise variables */

  /* Interpret the command-line parameters */
  get_command_arguments(argv,argc - 1,&params);

  /* Read in all the directory paths from the parameter file, CATHPARAM */
  c_get_paths(&file_name_ptr,params.run_64);

  /* Read the PDB file to check for missing chains */
  new_chain = check_chains(params.file_name,&chain_type);

  /* Rewrite the PDB file, fixing errors */
  rewrite_pdb_file(params.file_name,new_chain,chain_type);

  return(0);
}
