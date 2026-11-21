/***********************************************************************

  postmatch.c - Program to analyse HBADD matches to between enzyme
                reaction/cofactor molecules and the ligands in the
                PDB file. Writes out a molec.dat file to the PDB entry's
                PDBsum directory

************************************************************************

Date:-         16 Jul 2004

Written by:-   Roman Laskowski

----------------------------------------------------------------------*/

/* Datafiles
   ---------

Actual name      Description
-----------      -----------
match_all.dat    Input file containing the molecule matches
molecules.dat    Output file giving matches between enzyme reaction/
                 cofactor molecules and ligands in the structure

Compilation
-----------

cc -c postmatch.c
cc -c common.c
cc -c reapar.c
cc -o postmatch postmatch.o common.o reapar.o -lm -lz

Optimization flags:-  

cc -O2 -funroll-loops -c postmatch.c
cc -O2 -funroll-loops -c common.c
cc -O2 -funroll-loops -c reapar.c
cc -o postmatch postmatch.o common.o reapar.o -lm -lz

*/

/* Function calling-tree
   ---------------------

main
   -> c_get_paths (in common.c)
         -> get_param (in reapar.c)
         -> store_string (in common.c)
   -> read_matches_dat
         -> string_chomp (in common.c)
         -> get_tokens (in common.c)
         -> count_atoms
         -> create_match_entry
         -> free_tokens (in common.c)
   -> write_molecules_dat
         -> find_pdbsum_dir (in common.c)
         -> open_output_file (in common.c)

*/


#include "common.h"
#include "postmatch.h"


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

/* In reapar.c */
int get_param(char *token_name,char *file_name,int filename_len,
              int exit_on_error);

/***********************************************************************

create_match_entry  -  Create a match entry

***********************************************************************/

struct match *create_match_entry(struct match **fst_match_ptr,
                                 struct match **lst_match_ptr)
{
  struct match *match_ptr;
  struct match *first_match_ptr, *last_match_ptr;

  /* Initialise match pointers */
  match_ptr = NULL;
  first_match_ptr = *fst_match_ptr;
  last_match_ptr = *lst_match_ptr;

  /* Create match entry */
  match_ptr =
    (struct match*) malloc(sizeof(struct match));
  if (match_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for");
      printf(" struct match\n");
      printf("***        Program postmatch terminated with");
      printf(" error.\n");
      exit (1);
    }

  /* Initialise all the fields */
  match_ptr->pdb_code[0] = '\0';
  match_ptr->ec_no = 0;
  match_ptr->rstep_no = 0;
  match_ptr->component_no = 0;
  match_ptr->molec_id[0] = '\0';
  match_ptr->mol_type = -1;
  match_ptr->mol_name = &null;
  match_ptr->res_name[0] = '\0';
  match_ptr->res_num[0] = '\0';
  match_ptr->chain = '\0';
  match_ptr->similarity = 0.0;
  match_ptr->nmatched = 0;
  match_ptr->nmissed = 0;
  match_ptr->nmissing = 0;
  match_ptr->liguniq_no = 0;
  match_ptr->ligand_no = 0;
  match_ptr->wanted = TRUE;
  match_ptr->written = FALSE;

  /* Initialise pointer to next entry */
  match_ptr->next_match_ptr = NULL;

  /* If this is the first entry stored, then update
     pointer to head of linked list */
  if (first_match_ptr == NULL)
    first_match_ptr = match_ptr;

  /* Otherwise, make previous entry point to the current
     one */
  else
    last_match_ptr->next_match_ptr = match_ptr;
                      
  /* Save the current residue's pointer */
  last_match_ptr = match_ptr;

  /* Return current pointers */
  *fst_match_ptr = first_match_ptr;
  *lst_match_ptr = last_match_ptr;
  return(match_ptr);
}
/***********************************************************************

format_res_num  -  Interpret the residue-number as read in from file

***********************************************************************/

void format_res_num(char *res_num)
{
  char insertion_code, res_number[6], new_number[6];

  int cpos, got_number, idigit, ipos, ndigits;
  int at_end;
  int nchars;

  /* Initialise position of insertion code */
  cpos = 0;
  insertion_code = ' ';
  got_number = FALSE;
  ndigits = 0;
  for (idigit = 0; idigit < 5; idigit++)
    new_number[idigit] = ' ';
  at_end = FALSE;
  nchars = 0;

  /* Loop through the characters in the entered residue-number
     searching for end of number or insertion-code */
  for (ipos = 0; ipos < 5 && at_end == FALSE; ipos++)
    {
      /* If this is the end of the string, store end-position */
      if (res_num[ipos] == '\0')
        {
          cpos = ipos;
          at_end = TRUE;
        }

      /* If this is a digit, then have part of a residue number */
      else if ((res_num[ipos] >= '0' && res_num[ipos] <= '9') ||
               res_num[ipos] == '-')
        {
          res_number[ndigits] = res_num[ipos];
          ndigits++;
          got_number = TRUE;
        }

      /* If this is a character it might be a prefix to the
         residue number, or the insertion code */
      else if ((res_num[ipos] >= 'A' && res_num[ipos] <= 'Z') ||
               res_num[ipos] == '\'')
        {
          /* If we already have a number, this must be an
             insertion code */
          if (got_number == TRUE)
            {
              cpos = ipos;
              insertion_code = res_num[ipos];
            }

          /* Increment count of characters */
          nchars++;
        }

      /* If this is a space, have reached the end of the number */
      else if (res_num[ipos] == ' ')
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

count_atoms  -  If this is a chemical formula, count its number of atoms

***********************************************************************/

int count_atoms(char *token)
{
  char name[LINELEN + 1];

  int i, len;
  int small;

  char *string_ptr;

#define NMOLECS     27

  static char *small_molecule[NMOLECS] = {
    "Ag(+)", "Cl(-)", "CO", "CO2", "CO(2)", "Cu(+)", "Fe(2+)", "Fe(3+)",
    "H(+)", "H(2)", "HBr", "HCl", "Hg", "Hg(2+)", "H(2)O", "H(2)O(2)",
    "H(2)S", "K(+)", "Mg(2+)", "Mn(2+)", "N(2)", "Na(+)", "NH(3)",
    "Ni(2+)", "NO", "O(2)", "Ni(2+)"
  };

  /* Initialise variables */
  small = FALSE;
  strcpy(name,token);

  /* Remove (Out) and (In) at end, if present */
  if ((string_ptr = strstr(name,"(In)")) != NULL)
    string_ptr[0] = '\0';
  else if ((string_ptr = strstr(name,"(Out)")) != NULL)
    string_ptr[0] = '\0';

  /* Check for common cases */
  for (i = 0; i < NMOLECS && small == FALSE; i++)
    {
      if (!strcmp(name,small_molecule[i]))
        small = TRUE;
    }

  /* Return whether this is a small molecule */
  return(small);
}
/***********************************************************************

read_matches_dat  -  Read in and store the molecule matches for the given
                     PDB entry

***********************************************************************/

int read_matches_dat(struct match **fst_match_ptr)
{
  char input_line[LINELEN + 1], file_name[FILENAME_LEN];
  char mol_name[LINELEN + 1], number_string[20];

  char *string_ptr, *token[MAX_TOKEN];

  int ipos, itoken, line, nmatch, ntokens;
  int small;

  struct match *match_ptr;
  struct match *first_match_ptr, *last_match_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  file_name[0] = '\0';
  line = 0;
  nmatch = 0;

  /* Form the file name of the match_all.dat file */
  strcat(file_name,"match_all.dat");

  /* Initialise match pointers */
  match_ptr = NULL;
  first_match_ptr = NULL;
  last_match_ptr = NULL;

  /* Open the match_all.dat input file */
  if ((file_ptr = fopen(file_name,"r")) !=  NULL)
    {
      /* Loop through the file processing each record */
      while (fgets(input_line,LINELEN,file_ptr) != NULL)
        {
          /* Increment line-count */
          line++;

          /* Truncate the string to strip off the newline */
          string_chomp(input_line);

          /* Create a match entry */
          match_ptr
            = create_match_entry(&first_match_ptr,&last_match_ptr);

          /* Extract the space-separated tokens */
          ntokens = get_tokens(input_line,token,MAX_TOKEN,' ');

          /* PDB code and molecule id */
          strcpy(match_ptr->pdb_code,token[0]);
          match_ptr->pdb_code[4] = '\0';
          strcpy(match_ptr->molec_id,token[1]);
          match_ptr->molec_id[6] = '\0';

          /* Molecule typre (reactant, product or cofactor) */
          if (token[2][0] == 'R')
            match_ptr->mol_type = REACTANT;
          else if (token[2][0] == 'P')
            match_ptr->mol_type = PRODUCT;
          else if (token[2][0] == 'C')
            match_ptr->mol_type = COFACTOR;

          /* E.C. ref no, reaction step no and component number */
          strcpy(number_string,token[3]);
          match_ptr->ec_no = atoi(number_string);
          strcpy(number_string,token[4]);
          match_ptr->rstep_no = atoi(number_string);
          strcpy(number_string,token[5]);
          match_ptr->component_no = atoi(number_string);

          /* If molecule not matched, then check if it is a significant one */
          if (!strcmp(token[6],"UNMATCHED"))
            {
              /* If have the molecule name, can check for it */
              if (ntokens > 7)
                {
                  /* If molecule name consists of more than one token,
                     concatenate them */
                  strcpy(mol_name,token[7]);
                  if (ntokens > 8)
                    {
                      mol_name[0] = '\0';
                      for (itoken = 7; itoken < ntokens; itoken++)
                        {
                          if (itoken > 7)
                            strcat(mol_name," ");
                          strcat(mol_name,token[itoken]);
                        }
                    }
              
                  /* See if molecule is something insignificant like
                     water */
                  small = count_atoms(mol_name);

                  /* If this is a significant molecule, then count as
                     missed */
                  if (small == TRUE)
                    strcpy(match_ptr->res_name,"...");
                  else
                    {
                      /* Create a dummy match record for this one */
                      strcpy(match_ptr->res_name,"xxx");

                      /* Increment match count */
                      nmatch++;
                    }

                  /* Molecule name */
                  store_string(&(match_ptr->mol_name),token[7]);
                  match_ptr->chain = ' ';
                  match_ptr->similarity = 0;
                }
            }

          /* Otherwise, retrieve the rest of the data on the line */
          else
            {
              /* Get the residue name */
              strcpy(match_ptr->res_name,token[6]);

              /* Get the residue number and format it */
              format_res_num(token[8]);

              /* Get the chain, checking whether it might be blank */
              ipos = 8;
              if (!strcmp(token[13],"Lig:"))
                {
                  match_ptr->chain = token[8][0];
                  ipos++;
                }
              else
                match_ptr->chain = ' ';

              /* Percentage similarity */
              strcpy(number_string,token[ipos]);
              string_ptr = strchr(number_string,'%');
              if (string_ptr != NULL)
                string_ptr[0] = '\0';
              match_ptr->similarity = atof(number_string);
              ipos++;

              /* Numbers of matches and misses */
              strcpy(number_string,token[ipos]);
              match_ptr->nmatched = atoi(number_string);
              ipos++;
              strcpy(number_string,token[ipos]);
              match_ptr->nmissed = atoi(number_string);
              ipos++;
              strcpy(number_string,token[ipos]);
              match_ptr->nmissing = atoi(number_string);
              ipos++;
              ipos++;

              /* PDBsum ligand reference */
              strcpy(number_string,token[ipos]);
              match_ptr->liguniq_no = atoi(number_string);
              ipos++;
              strcpy(number_string,token[ipos]);
              match_ptr->ligand_no = atoi(number_string);
              ipos++;

              /* Form the molecule name */
              input_line[0] = '\0';
              for (itoken = ipos; itoken < ntokens; itoken++)
                {
                  if (itoken > ipos)
                    strcat(input_line," ");
                  strcat(input_line,token[itoken]);
                }
              store_string(&(match_ptr->mol_name),input_line);

              /* Increment count of matches read in */
              nmatch++;
            }

          /* Free the tokens */
          free_tokens(token,ntokens);
        }

      /* Close the match file */
      fclose(file_ptr);
    }

  /* If no match_all.dat file, then can't continue */
  else
    {
      /* Show error message and stop */
      printf("*** ERROR. Unable to open file: %s\n",file_name);
      exit(-1);
    }

  /* Return pointer to first match */
  *fst_match_ptr = first_match_ptr;

  /* Show number of records read in */
  printf("Number of matches read in from match_all.dat file:  %8d\n",
         nmatch);

  /* Return number of matches read in */
  return(nmatch);
}
/***********************************************************************

write_molecules_dat  -  Write out the molecules.dat file

***********************************************************************/

void write_molecules_dat(struct c_file_data file_name_ptr,
                         struct match *first_match_ptr)
{
  char file_name[FILENAME_LEN + 1], pdb_code[5];
  char pdbsum_dir[FILENAME_LEN + 1];
  char mol_type[20], number_string[20], type[10], subscript[10];
  char cofactors[LINELEN], reactants[LINELEN], products[LINELEN];
  char coverage[NTYPES + 1];

  int ec_no, i, rstep_no, nec_no, nrstep_no, nwritten;
  int count[NTYPES], dir_type, npresent[NTYPES], itype, pdb_type;
  int wanted;

  struct match *match_ptr, *other_match_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  nec_no = nwritten = 0;
  pdb_type = STANDARD_PDB;
  strcpy(pdb_code,first_match_ptr->pdb_code);
  for (i = 0; i < NTYPES; i++)
    npresent[i] = 0;
  
  /* Form the name of the PDBsum directory */
  dir_type = find_pdbsum_dir(pdb_code,file_name_ptr.pdbsum_template,
                             pdb_type,pdbsum_dir);
  if (dir_type == NOT_FOUND)
    {
      printf("*** ERROR. Unable to locate PDBsum directory for %s\n",
             pdb_code);
      exit(-1);
    }

  /* Form output file names */
  strcpy(file_name,pdbsum_dir);
  strcat(file_name,DIR_SEP);
  strcat(file_name,"molecules.dat");

  /* Open the molecules.dat output file */
  file_ptr = open_output_file(file_name,TRUE);

  /* Get pointer to the first match record */
  match_ptr = first_match_ptr;

  /* Loop over the matches */
  while (match_ptr != NULL)
    {
      /* Save the PDB code */
      strcpy(pdb_code,match_ptr->pdb_code);
      
      /* Save the E.C. counter */
      if (match_ptr->ec_no > nec_no)
        nec_no = match_ptr->ec_no;

      /* Get pointer to the next match record */
      other_match_ptr = match_ptr->next_match_ptr;

      /* Initialise flag */
      wanted = TRUE;

      /* If current Het Group has already been assigned, then don't
         write out */
      wanted = match_ptr->wanted;

      /* Loop over the other matches for this reaction step and see if
         there is a better reaction match than this one */
      while (other_match_ptr != NULL)
        {
          /* Check whether this is the same Het Group and it belongs to
             the same reaction step or is cofactor in the same reaction */
          if (!strcmp(other_match_ptr->res_name,match_ptr->res_name) &&
              other_match_ptr->ec_no == match_ptr->ec_no &&
              (other_match_ptr->rstep_no == match_ptr->rstep_no ||
               other_match_ptr->rstep_no == 0))
            {
              /* If current other molecule has a better similarity score
                 then don't want this one */
              if (other_match_ptr->similarity > match_ptr->similarity)
                wanted = FALSE;

              /* Otherwise mark that one as done */
              else
                other_match_ptr->wanted = FALSE;
            }

          /* Get the next match */
          other_match_ptr = other_match_ptr->next_match_ptr;
        }

      /* Check whether similarity score is above the minimum */
      if (match_ptr->similarity < MIN_SIMILARITY)
        wanted = FALSE;

      /* If a dummy match record, then don't want it written out */
      if (!strcmp(match_ptr->res_name,"xxx"))
        {
          wanted = FALSE;
          match_ptr->written = TRUE;
        }

      /* If match is wanted, write out to molecules.dat file */
      if (wanted == TRUE)
        {
          /* Mark this match as written */
          match_ptr->written = TRUE;

          /* Increment count of matches written out for this E.C.
             number */
          nwritten++;

          /* Write comment */
          fprintf(file_ptr,"# Ordered by entry's E.C. code\n");

          /* If match is to a cofactor, subscript identifies E.C. number
             and component */
          if (match_ptr->mol_type == COFACTOR)
            {
              strcpy(type,"COF");
              sprintf(subscript,"[%d][%d]",match_ptr->ec_no,
                      match_ptr->component_no);
            }

          /* Otherwise, for a reaction product or reactant, need to
             include the reaction step as well */
          else
            {
              strcpy(type,"REAC");
              sprintf(subscript,"[%d][%d][%d]",match_ptr->ec_no,
                      match_ptr->rstep_no,match_ptr->component_no);
            }

          /* Write out the Het Group name */
          fprintf(file_ptr,"%s_HETNAME%s %s\n",type,subscript,
                  match_ptr->res_name);

          /* Write out the matching molecule identifier */
          fprintf(file_ptr,"%s_MOLECULE%s %s\n",type,subscript,
                  match_ptr->molec_id);

          /* Write out the molecule name */
          fprintf(file_ptr,"%s_MOLNAME%s %s\n",type,subscript,
                  match_ptr->mol_name);

          /* Write out the ligand reference: unique ligand number */
          fprintf(file_ptr,"%s_LIG%s %d\n",type,subscript,
                  match_ptr->liguniq_no);

          /* Write out the ligand number */
          fprintf(file_ptr,"%s_LIGNUM%s %d\n",type,subscript,
                  match_ptr->ligand_no);

          /* Write out the percentage similarity */
          fprintf(file_ptr,"%s_MATCH%s %6.2f\n",type,subscript,
                  match_ptr->similarity);

          /* Write comment */
          fprintf(file_ptr,"# Ordered by ligand reference\n");

          /* Form subscript */
          strcpy(type,"LIG");
          sprintf(subscript,"[%d][%d]",match_ptr->liguniq_no,
                  match_ptr->ligand_no);

          /* Write corresponding E.C. ref */
          fprintf(file_ptr,"%s_EC_NO%s %d\n",type,subscript,
                  match_ptr->ec_no);

          /* Write out the molecule this ligand matches */
          fprintf(file_ptr,"%s_MOL%s %s\n",type,subscript,
                  match_ptr->molec_id);

          /* Write out the molecule type */
          if (match_ptr->mol_type == REACTANT)
            strcpy(mol_type,"reactant");
          else if (match_ptr->mol_type == PRODUCT)
            strcpy(mol_type,"product");
          else if (match_ptr->mol_type == COFACTOR)
            strcpy(mol_type,"cofactor");
          fprintf(file_ptr,"%s_MOLTYPE%s %s\n",type,subscript,
                  mol_type);

          /* Write out the molecule name */
          fprintf(file_ptr,"%s_MOLNAME%s %s\n",type,subscript,
                  match_ptr->mol_name);

          /* Write out the Het Group that matches (in case part of a
             larger ligand molecule) */
          fprintf(file_ptr,"%s_HETGROUP%s %s\n",type,subscript,
                  match_ptr->res_name);

          /* Write out the percentage similarity */
          fprintf(file_ptr,"%s_MATCH%s %6.2f\n",type,subscript,
                  match_ptr->similarity);
        }

      /* Get the next match */
      match_ptr = match_ptr->next_match_ptr;
    }

  /* Loop over the different enzyme classes */
  for (ec_no = 0; ec_no < nec_no; ec_no++)
    {
      /* Initialise the counts */
      for (i = 0; i < NTYPES; i++)
        count[i] = npresent[i];
      nrstep_no = 0;
      cofactors[0] = '\0';

      /* Get pointer to the first match record */
      match_ptr = first_match_ptr;

      /* Loop over the matches */
      while (match_ptr != NULL)
        {
          /* If this is the right enzyme class, then proceed */
          if (match_ptr->ec_no == ec_no + 1 &&
              strcmp(match_ptr->res_name,"..."))
            {
              /* If this is a reaction step, then increment count */
              if (match_ptr->rstep_no > nrstep_no)
                nrstep_no = match_ptr->rstep_no;

              /* If this is a cofactor, then increment counts */
              if (match_ptr->rstep_no == 0)
                {
                  /* Increment cofactor count for this EC class */
                  count[COFACTOR]++;

                  /* If matched, increment count of cofactor matches */
                  if (match_ptr->written == TRUE &&
                      strcmp(match_ptr->res_name,"xxx"))
                    {
                      npresent[COFACTOR]++;
                      sprintf(number_string," %8.3f",
                              match_ptr->similarity);
                      strcat(cofactors,number_string);
                    }
                  else
                    strcat(cofactors," ........");
                }
            }

          /* Get the next match */
          match_ptr = match_ptr->next_match_ptr;
        }

      /* Loop over the reaction steps for the current E.C. class and 
         count reactant and product molecules present in the PDB file */
      for (rstep_no = 0; rstep_no < nrstep_no; rstep_no++)
        {
          /* Initialise counts of product and reactant molecules */
          count[PRODUCT] = count[REACTANT] = 0;
          npresent[PRODUCT] = npresent[REACTANT] = 0;
          reactants[0] = products[0] = '\0';

          /* Get pointer to the first match record */
          match_ptr = first_match_ptr;

          /* Loop over the matches */
          while (match_ptr != NULL)
            {
              /* If this is the right enzyme class and right reaction
                 step, then proceed */
              if (match_ptr->ec_no == ec_no + 1 &&
                  match_ptr->rstep_no == rstep_no + 1 &&
                  strcmp(match_ptr->res_name,"..."))
                {
                  /* Get the molecule type */
                  itype = match_ptr->mol_type;

                  /* Increment appropriate count */
                  count[itype]++;

                  /* If matched, increment count of matches */
                  if (match_ptr->written == TRUE &&
                      strcmp(match_ptr->res_name,"xxx"))
                    {
                      npresent[itype]++;
                      sprintf(number_string," %8.3f",
                              match_ptr->similarity);
                    }
                  else
                    strcpy(number_string," ........");

                  /* Append similarity to appropriate string */
                  if (itype == REACTANT)
                    strcat(reactants,number_string);
                  else
                    strcat(products,number_string);
                }

              /* Get the next match */
              match_ptr = match_ptr->next_match_ptr;
            }

          /* Initialise coverage string */
          coverage[0] = '\0';

          /* Calculate the coverage flags */
          for (itype = 0; itype < NTYPES; itype++)
            {
              /* Determine flag according to coverage */
              if (count[itype] == 0)
                strcat(coverage,"-");
              else if (npresent[itype] == 0)
                strcat(coverage,".");
              else if (npresent[itype] == count[itype])
                strcat(coverage,"*");
              else
                strcat(coverage,"x");
            }

          /* Write out the counts and flags for this reaction step */
          fprintf(file_ptr,"COVERAGE[%d][%d] %s  %s  R: %s P: %s C: %s\n",
                  ec_no + 1,rstep_no + 1,pdb_code,coverage,
                  reactants,products,cofactors);
        }
    }

  /* Close the molecules.dat file */
  fclose(file_ptr);
}

/***********************************************************************

                             M   A   I   N

***********************************************************************/

int main(int argc,char *argv[])
{
  int nmatch;

  struct c_file_data file_name_ptr;

  struct match *first_match_ptr;

  /* Get required paths and file names */
  c_get_paths(&file_name_ptr,TRUE);

  /* Read in all the molecule matches */
  nmatch = read_matches_dat(&first_match_ptr);

  /* Write the molecules.dat file to this PDB entry's PDBsum directory */
  write_molecules_dat(file_name_ptr,first_match_ptr);

  return(0);
}
