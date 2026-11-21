/***********************************************************************

  readpdb.c - Program to read in data from PDB file

************************************************************************

Date:-         12 Oct 2009

Written by:-   Roman Laskowski

----------------------------------------------------------------------*/

/* Datafiles
   ---------

Actual name  Description
-----------  -----------
<file_name>  Input PDB file

Compilation
-----------

cc -c readpdb.c

*/

/* Function calling-tree
   ---------------------

main
   -> r_find_pdb_file
         -> form_filename (in common.c)
               -> store_string (in common.c)
   -> r_find_residue
   -> r_get_dfn_data
         -> form_filename (in common.c)
               -> store_string (in common.c)
         -> string_chop (in common.c)
         -> create_r_rasdef_entry
         -> r_get_ligand_residues
               -> r_get_token
               -> r_format_res_number
               -> r_remove_leading_blanks
               -> create_r_rasdef_entry
               -> store_string (in common.c)
         -> r_to_lower
         -> store_string (in common.c)
         -> r_number_ligands
         -> map_to_ligplots
               -> string_chomp (in common.c)
               -> get_tokens (in common.c)
               -> split_metal_range
               -> create_r_rasdef_entry
               -> r_remove_leading_blanks (in common.c)
               -> store_string (in common.c)
   -> r_get_pdb_data
         -> r_get_dfn_data
               -> create_r_ligand_record
         -> r_find_pdb_file
               -> form_filename (in common.c)
                     -> store_string (in common.c)
         -> r_open_pdb_file
         -> r_read_pdb_file
               -> create_r_pdb_entry
               -> string_chop (in common.c)
               -> r_store_data
                     -> to_lower (in common.c)
                     -> create_r_info_record
                     -> get_tokens (in common.c)
                     -> store_string (in common.c)
                     -> r_get_reference
                           -> string_chop (in common.c)
                           -> remove_marks (in common.c)
                           -> type_set (in common.c)
                     -> remove_leading_blanks
                     -> create_r_ref_record
                     -> type_set (in common.c)
                     -> fix_upper (in common.c)
               -> process_line
               -> get_atom_data
                     -> r_get_atom_details
                           -> remove_leading_blanks
                     -> r_check_for_duplicate
                     -> create_r_chain_record
                     -> r_assign_to_object
                     -> create_r_residue_record
                           -> r_get_aa_code
                     -> create_r_atom_record
                           -> get_element
                                 -> is_hydrogen_atom
                                       -> r_check_for_metal
                                 -> r_check_for_metal
               -> get_conect_data
                     -> r_find_atom
                     -> r_check_for_metal
                     -> create_r_conect_record
                           -> create_r_rbond_record
         -> identify_residues
               -> get_residue_type
         -> r_assign_to_molecules
               -> create_r_molec_record
         -> r_calc_limits
               -> r_get_max_min_coords
   -> r_get_pdb_file
         -> r_open_pdb_file (as above)
         -> r_read_pdb_file (as above)
         -> r_res_assignments
               -> classify_residues
                     -> get_residue_type
                     -> get_het_desc
                     -> excise_metal
                          -> create_r_chain_record
                          -> update_residues
               -> shift_metals
               -> adjust_residues
               -> classify_chains
               -> join_chains
               -> reclassify_chains
               -> split_ligands
                     -> find_joins
                           -> check_for_bond
                     -> identify_molecule
                           -> create_r_list_record
                           -> free_r_list_memory
                     -> create_r_chain_record
                     -> update_residues
               -> single_residue_ligands
                     -> r_check_for_metal
               -> form_ligands
                     -> create_r_ligand_record
                     -> add_res_name
                     -> r_format_residue_number
                     -> store_string
               -> identify_unique_ligands
               -> r_sort_ligands
                     -> r_shell_sort
               -> assign_chains
                     -> create_r_rasdef_record
                     -> assign_rasdefs
         -> r_assign_to_molecules
               -> create_r_molec_record
         -> r_calc_limits
               -> r_get_max_min_coords
   -> r_get_pdb_chain
         -> r_get_pdb_data (as above)
   -> r_get_res_name
   -> r_free_all_data
         -> r_free_pdb_data
   -> r_write_pdb_coords
   -> r_list_rasdef_entries
*/


#include "common.h"
#include "readpdb.h"

/* Prototypes
   ---------- */

/* In common.c */
void fix_upper(char *name);
char *form_filename(char *pdb_code,char *filename_template);
int get_tokens(char *line,char **token,int max_token,char token_sep);
void to_lower(char *search_string,int length);
void remove_leading_blanks(char *string);
void remove_marks(char *string);
void right_justify(char *string,int length);
void store_string(char **string_ptr,char *string);
int string_chomp(char *string);
int string_chop(char *string);
void type_set(char *name,int cap_start);


/***********************************************************************

r_find_pdb_file  -  Locate the PDB file from its code, trying the
                    standard PDB directories, or MCSG and Uploads
                    directories, as appropriate

***********************************************************************/

int r_find_pdb_file(char *pdb_code,char *pdb_template[NPDB_DIR],
                    char *file_name,int pdb_type)
{
  char *fname;

  int start, trial, trying, dir;

  FILE *file_ptr;

  /* Initialise */
  file_name[0] = '\0';
  dir = NOT_FOUND;

  /* Determine the starting PDB directory to search */
  start = 0;
  if (pdb_type == UPLOAD_PDB)
    start = PDB_UPLOAD;
  else if (pdb_type == MCSG_PDB)
    start = PDB_MCSG;
  else if (pdb_type == ALPHA_FOLD_PDB)
    start = PDB_ALPHA_FOLD;

  /* Loop over the possible directories the file might be in, until we
     find it */
  for (trial = start; trial < NPDB_DIR && dir == NOT_FOUND; trial++)
    {
      /* Form the trial name */
      fname = form_filename(pdb_code,pdb_template[trial]);

      /* Open the file, if it exists */
      if ((file_ptr = fopen(fname,"r")) !=  NULL)
        {
          /* Set the appropriate flag */
          dir = trial;

          /* Save the file name */
          strcpy(file_name,fname);

          /* Close the file */
          fclose(file_ptr);
        }

      /* Free memory of file name */
      free(fname);
    }

  /* Rule out chance hits in the wrong directory */
  if (pdb_type != FIND_FIRST)
    {
      if (dir == PDB_UPLOAD)
        {
          if (pdb_type != UPLOAD_PDB)
            dir = NOT_FOUND;
        }
      else if (dir == PDB_MCSG)
        {
          if (pdb_type != MCSG_PDB)
            dir = NOT_FOUND;
        }
      else
        {
          if (pdb_type == UPLOAD_PDB || pdb_type == MCSG_PDB)
            dir = NOT_FOUND;
        }
    }
  if (dir == NOT_FOUND)
    file_name[0] = '\0';

  /* Return the file type */
  return(dir);
}
/***********************************************************************

r_find_residue  -  Find the given residue entry

***********************************************************************/

struct r_residue *r_find_residue(struct r_residue *first_r_residue_ptr,
                                 char *res_name,char *res_num,char chain)
{
  struct r_residue *r_residue_ptr, *found_r_residue_ptr;

  /* Initialise variables */
  found_r_residue_ptr = NULL;

  /* Get pointer to the first residue */
  r_residue_ptr = first_r_residue_ptr;

  /* Loop over all the residues to find the one of interest */
  while (r_residue_ptr != NULL && found_r_residue_ptr == NULL)
    {
      /* If this is the right residue, store pointer */
      if (chain == r_residue_ptr->chain &&
          !strcmp(res_num,r_residue_ptr->res_num) &&
          (!strcmp(res_name,r_residue_ptr->res_name) || res_name[0] == '\0'))
        found_r_residue_ptr = r_residue_ptr;

      /* Get pointer to the next residue */
      r_residue_ptr = r_residue_ptr->next_r_residue_ptr;
    }

  /* Return the residue pointer */
  return(found_r_residue_ptr);
}
/***********************************************************************

create_r_rasdef_entry  -  Create a r_rasdef entry

***********************************************************************/

struct r_rasdef *create_r_rasdef_entry(struct r_rasdef **fst_r_rasdef_ptr,
                                       struct r_rasdef **lst_r_rasdef_ptr,
                                       int type,int number,char *res_name,
                                       char *res_num,char chain,
                                       char *colour,char *ligand_def,
                                       char *ligand_from,char *ligand_to,
                                       int nresid)
{
  struct r_rasdef *r_rasdef_ptr;
  struct r_rasdef *first_r_rasdef_ptr, *last_r_rasdef_ptr;

  /* Initialise r_rasdef pointers */
  r_rasdef_ptr = NULL;
  first_r_rasdef_ptr = *fst_r_rasdef_ptr;
  last_r_rasdef_ptr = *lst_r_rasdef_ptr;

  /* Create r_rasdef entry */
  r_rasdef_ptr =
    (struct r_rasdef*) malloc(sizeof(struct r_rasdef));
  if (r_rasdef_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for");
      printf(" struct r_rasdef\n");
      printf("***        Program postssm terminated with");
      printf(" error.\n");
      exit (1);
    }

  /* Store the details */
  r_rasdef_ptr->type = type;
  r_rasdef_ptr->number = number;
  strcpy(r_rasdef_ptr->res_name,res_name);
  strcpy(r_rasdef_ptr->res_num,res_num);
  r_rasdef_ptr->chain = chain;
  strcpy(r_rasdef_ptr->colour,colour);
  r_rasdef_ptr->nresid = nresid;

  /* Initialise other fields */
  strcpy(r_rasdef_ptr->ligand_from,ligand_from);
  strcpy(r_rasdef_ptr->ligand_to,ligand_to);
  r_rasdef_ptr->ligplot_id = &null;
  r_rasdef_ptr->ligand_name = &null;
  r_rasdef_ptr->ligand_range = &null;
  r_rasdef_ptr->copy = 0;
  r_rasdef_ptr->ncopy = 1;
  r_rasdef_ptr->wanted = TRUE;
  r_rasdef_ptr->first_r_residue_ptr = NULL;

  /* Initialise pointer to next entry */
  r_rasdef_ptr->next_r_rasdef_ptr = NULL;

  /* If this is the first entry stored, then update
     pointer to head of linked list */
  if (first_r_rasdef_ptr == NULL)
    first_r_rasdef_ptr = r_rasdef_ptr;

  /* Otherwise, make previous entry point to the current
     one */
  else
    last_r_rasdef_ptr->next_r_rasdef_ptr = r_rasdef_ptr;
                      
  /* Save the current residue's pointer */
  last_r_rasdef_ptr = r_rasdef_ptr;

  /* Return current pointers */
  *fst_r_rasdef_ptr = first_r_rasdef_ptr;
  *lst_r_rasdef_ptr = last_r_rasdef_ptr;
  return(r_rasdef_ptr);
}
/***********************************************************************

r_format_res_number  -  Interpret the residue-number from the rasmol.def
                        or catres.dat line

***********************************************************************/

void r_format_res_number(char *tmpstr,char *res_num)
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
  strcpy(new_number,"     ");
  strcpy(res_number,"     ");

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

r_get_token - Get the next space- or tab-delimited token from the given
              input line, starting at the given position

**********************************************************************/  

int r_get_token(char *line,int line_len,int *ipos,char *token,
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

r_remove_leading_blanks  -  Remove leading blank spaces from the given
                            string

***********************************************************************/

void r_remove_leading_blanks(char *string)
{
  char *tmp;

  int ipos, fpos, len;

  /* If no leading blanks, then return */
  if (string[0] != ' ' && string[0] != '_')
    return;

  /* Get the length of the string */
  len = strlen(string);

  /* Allocate memory for temporary string */
  tmp = (char *) malloc(sizeof(char)*(len + 1));

  /* Initialise position of first non-blank character */
  fpos = 0;

  /* Loop through to first non-blank character */
  for (ipos = 0; ipos < len && fpos == 0; ipos++)
    if (string[ipos] != ' ' && string[ipos] != '_')
      fpos = ipos;

  /* Shift whole string across */
  if (fpos > 0)
    {
      strcpy(tmp,string+fpos);
      strcpy(string,tmp);
    }
  else
    string[0] = '\0';

  /* Free up memory taken by the temporary string */
  free(tmp);
}
/***********************************************************************

r_get_ligand_residues  -  Extract the ligand details on the current line

***********************************************************************/

int r_get_ligand_residues(char *input_line,struct r_rasdef **fst_r_rasdef_ptr,
                          struct r_rasdef **lst_r_rasdef_ptr,int number,
                          int ndef,struct r_ligand *r_ligand_ptr,
                          char *ligand_name)
{
  char chain, chain_str[2], res_name[4], residue_name[4], res_num[6];
  char colour[R_COLNAM_LEN], number_string[10];
  char token_string[R_TOKEN_LEN + 1], tmp_string[R_TOKEN_LEN + 1];
  char ligand_from[7], ligand_to[7];
  char first_residue[7], last_residue[7];

  char *string_ptr;

  int len, lpos, nresid, token_len, type;
  int name_end;
  int done;

  struct r_rasdef *r_rasdef_ptr;
  struct r_rasdef *first_r_rasdef_ptr, *last_r_rasdef_ptr, *start_r_rasdef_ptr;

  /* Initialise variables */
  chain = ' ';
  nresid = 0;
  done = FALSE;
  ligand_name[0] = '\0';
  type = R_LIGAND;
  colour[0] = '\0';
  first_residue[0] = '\0';
  last_residue[0] = '\0';

  /* Initialise r_rasdef pointers */
  r_rasdef_ptr = NULL;
  first_r_rasdef_ptr = *fst_r_rasdef_ptr;
  last_r_rasdef_ptr = *lst_r_rasdef_ptr;
  start_r_rasdef_ptr = NULL;

  /* Get the line length */
  len = strlen(input_line);
  lpos = 0;

  /* Loop until line has been processed */
  while (done == FALSE)
    {
      /* Get the next token on the line */
      token_len = r_get_token(input_line,R_LINELEN,&lpos,token_string,
                              R_TOKEN_LEN,&done);

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

          /* Format residue name and add to current ligand name */
          if (ligand_name[0] != '\0')
            strcat(ligand_name,"-");

          /* Remove leading blanks */
          strcpy(residue_name,res_name);
          r_remove_leading_blanks(residue_name);
          strcat(ligand_name,residue_name);

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
          r_format_res_number(number_string,res_num);

          /* Copy chain into chain-string */
          chain_str[0] = chain;
          chain_str[1] = '\0';

          /* If this is the first residue of this ligand, store as first */
          if (first_residue[0] == '\0')
            {
              /* Copy across and add chain-id */
              strcpy(first_residue,number_string);
              strcat(first_residue,chain_str);

              /* Store in ligand record */
              strcpy(r_ligand_ptr->first_residue,res_name);
              strcat(r_ligand_ptr->first_residue,res_num);
              strcat(r_ligand_ptr->first_residue,chain_str);
            }

          /* Store as this ligand's last residue */
          strcpy(r_ligand_ptr->last_residue,res_name);
          strcat(r_ligand_ptr->last_residue,res_num);
          strcat(r_ligand_ptr->last_residue,chain_str);

          /* Copy residue number and chain-id across as last residue of
             current ligand */
          strcpy(last_residue,number_string);
          strcat(last_residue,chain_str);

          /* Increment count of residues belonging to this ligand */
          nresid++;
          r_ligand_ptr->nresid++;

          /* Blank out unwanted fields */
          ligand_from[0] = '\0';
          ligand_to[0] = '\0';

          /* Create a r_rasdef entry for this ligand residue */
          r_rasdef_ptr
            = create_r_rasdef_entry(&first_r_rasdef_ptr,
                                  &last_r_rasdef_ptr,type,number,
                                  res_name,res_num,chain,colour,
                                  input_line,ligand_from,ligand_to,
                                  nresid);

          /* If token contains a comma, remove */
          string_ptr = strchr(token_string,',');
          if (string_ptr != NULL)
            string_ptr[0] = '\0';

          /* Store the token in the range field */
          store_string(&(r_rasdef_ptr->ligand_range),token_string);

          /* If this is the first residue record for this ligand,
             store starting pointer */
          if (start_r_rasdef_ptr == NULL)
            start_r_rasdef_ptr = r_rasdef_ptr;
        }
    }

  /* If have a ligand name, store in all r_rasdef records relating to this
     ligand */
  if (ligand_name[0] != '\0')
    {
      /* Get pointer to the starting RasMol definition record */
      r_rasdef_ptr = start_r_rasdef_ptr;

      /* Loop through all the definitions to store current ligand
         name */
      while (r_rasdef_ptr != NULL)
        {
          /* Store ligand name */
          store_string(&(r_rasdef_ptr->ligand_name),ligand_name);

          /* Get next definition record */
          r_rasdef_ptr = r_rasdef_ptr->next_r_rasdef_ptr;
        }

      /* Store name as residue sequence in the ligand record */
      store_string(&(r_ligand_ptr->residue_sequence),ligand_name);
    }

  /* Return current pointers */
  *fst_r_rasdef_ptr = first_r_rasdef_ptr;
  *lst_r_rasdef_ptr = last_r_rasdef_ptr;

  /* Return updated number of r_rasdef entries */
  return(ndef);
}
/***********************************************************************

create_r_ligand_record  -  Create and initialise a new ligand record

***********************************************************************/

struct r_ligand *create_r_ligand_record(struct r_ligand **fst_r_ligand_ptr,
                                        struct r_ligand **lst_r_ligand_ptr,
                                        struct r_chain *r_chain_ptr,
                                        int number,int metal)
{
  int i;

  struct r_ligand *r_ligand_ptr, *first_r_ligand_ptr, *last_r_ligand_ptr;

  /* Initialise ligand pointers */
  r_ligand_ptr = NULL;
  first_r_ligand_ptr = *fst_r_ligand_ptr;
  last_r_ligand_ptr = *lst_r_ligand_ptr;

  /* Allocate memory for structure to hold ligand info */
  r_ligand_ptr = (struct r_ligand *) malloc(sizeof(struct r_ligand));
  if (r_ligand_ptr == NULL)
    {
      printf("*** Can't allocate memory for struct r_ligand\n");
      exit (1);
    }

  /* If this is the very first ligand, then store pointer */
  if (first_r_ligand_ptr == NULL)
    first_r_ligand_ptr = r_ligand_ptr;

  /* Add link from previous ligand to the current one */
  if (last_r_ligand_ptr != NULL)
    last_r_ligand_ptr->next_r_ligand_ptr = r_ligand_ptr;
  last_r_ligand_ptr = r_ligand_ptr;

  /* Store the ligand details */
  r_ligand_ptr->first_residue[0] = '\0';
  r_ligand_ptr->last_residue[0] = '\0';
  r_ligand_ptr->number = number;
  r_ligand_ptr->metal = metal;
  r_ligand_ptr->metal_colour = 0;
  r_ligand_ptr->metal_name[0] = '\0';
  r_ligand_ptr->type = 0;
  r_ligand_ptr->count = 1;
  r_ligand_ptr->unique = TRUE;
  r_ligand_ptr->unique_count = 0;
  r_ligand_ptr->ncopy = 0;
  r_ligand_ptr->duplicate_count = 0;
  r_ligand_ptr->nduplicates = 0;
  r_ligand_ptr->residue_list = &null;
  r_ligand_ptr->residue_sequence = &null;
  r_ligand_ptr->residue_range = &null;
  r_ligand_ptr->natoms = 0;
  r_ligand_ptr->sort_score = 0;
  r_ligand_ptr->got_bonds = FALSE;
  r_ligand_ptr->r_rasdef_ptr = NULL;
  if (r_chain_ptr != NULL)
    {
      r_ligand_ptr->nresid = r_chain_ptr->nresid;
      r_ligand_ptr->first_r_residue_ptr
        = r_chain_ptr->first_r_residue_ptr;
    }
  else
    {
      r_ligand_ptr->nresid = 0;
      r_ligand_ptr->first_r_residue_ptr = NULL;
    }
  r_ligand_ptr->first_r_list_ptr = NULL;
  r_ligand_ptr->duplicate_r_ligand_ptr = NULL;
  r_ligand_ptr->next_r_ligand_ptr = NULL;

  /* Return current pointers */
  *fst_r_ligand_ptr = first_r_ligand_ptr;
  *lst_r_ligand_ptr = last_r_ligand_ptr;

  /* Return new ligand pointer */
  return(r_ligand_ptr);
}
/***********************************************************************

r_to_lower  -  Convert given character string to lower case

***********************************************************************/

void r_to_lower(char *search_string,int length)
{
  int ichar;

  int ipos;

  /* Loop through all the character positions, converting any letters
     to upper case */
  for (ipos = 0; ipos < length; ipos++)
    {
      /* Check the current character */
      ichar = search_string[ipos] - 'A';
      if (ichar >= 0 && ichar < 26)
        {
          ichar = ichar + 'a';
          search_string[ipos] = ichar;
        }
    }
}
/***********************************************************************

r_number_ligands  -  Identify and number unique ligands and copies

***********************************************************************/

int r_number_ligands(struct r_rasdef *first_r_rasdef_ptr)
{
  int copy, nligands, number;

  struct r_rasdef *r_rasdef_ptr, *r_other_rasdef_ptr;

  /* Initialise variables */
  copy = 0;
  nligands = 0;
  number = -1;

  /* Get pointer to the first RasMol definition record */
  r_rasdef_ptr = first_r_rasdef_ptr;

  /* Loop through all the definition records */
  while (r_rasdef_ptr != NULL)
    {
      /* If this is a ligand record, then check against all others */
      if (r_rasdef_ptr->type == R_LIGAND && r_rasdef_ptr->copy == 0)
        {
          /* Initialise the count for this ligand */
          copy = r_rasdef_ptr->copy = 1;
          number = r_rasdef_ptr->number;

          /* Increment ligand count */
          nligands++;

          /* Save ligand number */
          r_rasdef_ptr->number = nligands;

          /* Get pointer to the next RasMol definition record */
          r_other_rasdef_ptr = r_rasdef_ptr->next_r_rasdef_ptr;

          /* Loop through all other definition records to find ligands
             matching the current one */
          while (r_other_rasdef_ptr != NULL)
            {
              /* If this is a matching ligand record, then store count */
              if (r_other_rasdef_ptr->type == R_LIGAND)
                {
                  /* If this is same ligand, then save copy number */
                  if (r_other_rasdef_ptr->number == number)
                    {
                      /* Save copy and ligand numbers */
                      r_other_rasdef_ptr->copy = copy;
                      r_other_rasdef_ptr->number = nligands;
                    }

                  /* Otherwise, if this is a copy of our ligand, increment
                     the copy count and save */
                  else if (!strcmp(r_other_rasdef_ptr->ligand_name,
                                   r_rasdef_ptr->ligand_name))
                    {
                      /* Increment copy count */
                      copy++;

                      /* Save the copy count */
                      r_other_rasdef_ptr->copy = copy;

                      /* Save the current number */
                      number = r_other_rasdef_ptr->number;

                      /* Save ligand number */
                      r_other_rasdef_ptr->number = nligands;
                    }
                }

              /* Get next definition record */
              r_other_rasdef_ptr = r_other_rasdef_ptr->next_r_rasdef_ptr;
            }
        }

      /* Get next definition record */
      r_rasdef_ptr = r_rasdef_ptr->next_r_rasdef_ptr;
    }

   /* Return number of ligands */
   return(nligands);
}
/***********************************************************************

 split_metal_range -  Split the given metal range into residue name,
                      number and chain id

***********************************************************************/

void split_metal_range(char *range,char *res_name,char *res_num,
                       char *chn)
{
  char ch, chain, number_string[20];

  char *string_ptr, *token[R_MAX_TOKEN];

  int ipos, len, npos, ntokens, res_no;

  /* Initialise */
  chain = ' ';
  res_name[0] = res_num[0] = '\0';

  /* Find the chain id */
  string_ptr = strchr(range,'(');
  if (string_ptr != NULL)
    {
      /* Get the chain id */
      chain = string_ptr[1];

      /* Strip it off */
      string_ptr[0] = '\0';
    }

  /* Get the tokens making up this range */
  ntokens = get_tokens(range,token,R_MAX_TOKEN,' ');

  /* If have 2 tokens, then residue name and number are separate */
  if (ntokens == 2)
    {
      /* Save the residue name and number */
      strcpy(res_name,token[0]);
      strcpy(res_num,token[1]);
    }

  /* Otherwise, we have to separate the residue name from the number */
  else
    {
      /* Save the first token */
      strcpy(number_string,token[0]);

      /* Get its length */
      len = strlen(number_string);

      /* Find the location of the first numeral */
      npos = -1;
      for (ipos = 0; ipos < len && npos == -1; ipos++)
        {
          /* Get this character */
          ch = number_string[ipos];

          /* If numeric, then we're done */
          if (ch >= '0' && ch <= '9')
            npos = ipos;
        }

      /* Split the residue name and number */
      if (npos > -1)
        {
          /* Get the residue name */
          strncpy(res_name,number_string,npos);
          res_name[npos] = '\0';

          /* Get the residue number */
          strcpy(res_num,number_string + npos);
        }
    }

  /* Remove trailing spaces from residue name */
  string_ptr = strchr(res_name,' ');
  if (string_ptr != NULL)
    string_ptr[0] = '\0';

  /* Right-justify it */
  right_justify(res_name,4);

  /* If last character of residue number is numeric,
     append a space */
  len = strlen(res_num);
  ch = res_num[len - 1];
  if (ch >= '0' && ch <= '9')
    {
      /* Format the number */
      res_no = atoi(res_num);
      sprintf(res_num,"%4d ",res_no);
    }
  else
    {
      /* Format the number */
      ch = res_num[len - 1];
      strcpy(number_string,res_num);
      number_string[len - 1] = '\0';
      res_no = atoi(number_string);
      sprintf(res_num,"%4d%c",res_no,ch);
    }

  /* Return the chain id */
  *chn = chain;
}
/***********************************************************************

map_to_ligplots -  Process all the ligand entries to map them onto the
                   PDBsum ligplot in which they appear

***********************************************************************/

int map_to_ligplots(char *pdbsum_dir,struct r_rasdef **fst_r_rasdef_ptr,
                    struct r_rasdef **lst_r_rasdef_ptr,
                    struct r_ligand **fst_r_ligand_ptr,
                    struct r_ligand **lst_r_ligand_ptr,int *nmet)
{
  char file_name[R_FILENAME_LEN], input_line[R_LONG_LINELEN + 1];
  char colour[R_COLNAM_LEN], number_string[20], range[R_LONG_LINELEN];
  char chain, ligplot_id[8], metal_name[20], res_name[4], res_num[6];
  char tmp[40];

  char *string_ptr, *token[R_MAX_TOKEN];

  int copy, len, loop, ndef, nmetal, number, old_number, res_no, type;
  int ntokens, num1, num2, num3;
  int done, metal, metal_duplic, search;

  struct r_ligand *r_ligand_ptr, *first_r_ligand_ptr, *last_r_ligand_ptr;
  struct r_rasdef *r_rasdef_ptr, *first_r_rasdef_ptr, *last_r_rasdef_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  colour[0] = '\0';
  copy = 0;
  metal = FALSE;
  ndef = 0;
  nmetal = *nmet;
  number = old_number = 0;
  range[0] = tmp[0] = '\0';

  /* Initialise pointers */
  r_ligand_ptr = NULL;
  first_r_ligand_ptr = *fst_r_ligand_ptr;
  last_r_ligand_ptr = *lst_r_ligand_ptr;
  r_rasdef_ptr = NULL;
  first_r_rasdef_ptr = *fst_r_rasdef_ptr;
  last_r_rasdef_ptr = *lst_r_rasdef_ptr;

  /* Loop through all the definition records */
  while (r_rasdef_ptr != NULL)
    {
      /* Save this as the last rasdef record */
      last_r_rasdef_ptr = r_rasdef_ptr;
      
      /* Get next definition record */
      r_rasdef_ptr = r_rasdef_ptr->next_r_rasdef_ptr;
    }
  
  /* Form the name of the ligands.dat file */
  strcpy(file_name,pdbsum_dir);
  strcat(file_name,"/ligands.dat");

  /* Open the file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    return(ndef);

  /* Loop through the file, line by line */
  while (fgets(input_line,R_LONG_LINELEN,file_ptr) != NULL)
    {
      /* Truncate the input line */
       len = string_chomp(input_line);

       /* Initialise flag */
       metal = metal_duplic = search = FALSE;

       /* If ligand range, then store its LIGPLOT number */
       if (!strncmp(input_line,"LIG_RANGE[",10) ||
           !strncmp(input_line,"MET_RANGE[",10) ||
           !strncmp(input_line,"LIG_DUPLIC[",11) ||
           !strncmp(input_line,"MET_DUPLIC[",11))
         {
           /* Initialise */
           num1 = num2 = num3 = 0;
           ligplot_id[0] = '\0';
           tmp[0] = '\0';

           /* Get the ligand range */
           string_ptr = strchr(input_line,' ');
           if (string_ptr != NULL)
             {
               /* Save the ligand range and delete */
               strcpy(range,string_ptr + 1);
               string_ptr[0] = '\0';
             }

           /* Get the first open bracket */
           string_ptr = strchr(input_line,'[');
           if (string_ptr != NULL)
             strcpy(tmp,string_ptr + 1);

           /* Extract the ligplot numbers */
           ntokens = get_tokens(tmp,token,R_MAX_TOKEN,']');

           /* If have at least two numbers, then extract */
           if (ntokens > 2)
             {
               /* Get the first and second numbers */
               number = num1 = atoi(token[0]);
               copy = num2 = atoi(token[1] + 1);
               if (ntokens > 3)
                 copy = num3 = atoi(token[2] + 1) + 1;

               /* Form the LIGPLOT id */
               sprintf(ligplot_id,"%2.2d_%2.2d",num1,num2);

               /* Set flag to search for this ligand range */
               search = TRUE;

               /* If metal, set flag */
               if (!strncmp(input_line,"MET_RANGE[",10))
                 metal = TRUE;
               else if (!strncmp(input_line,"MET_DUPLIC[",11))
                 metal_duplic = TRUE;
             }
         }

       /* If this is the metal colour, then store */
       else if (!strncmp(input_line,"METAL_COLOUR[",11))
         {
           /* Get the start of the colour name */
           string_ptr = strchr(input_line,' ');
           if (string_ptr != NULL)
             strcpy(colour,string_ptr + 1);
         }

       /* For metal, extract the residue name and number */
       if (metal == TRUE || metal_duplic == TRUE)
         {
           /* Define type */
           type = R_METAL;

           /* Split the metal range into residue name, number and
              chain id */
           split_metal_range(range,res_name,res_num,&chain);
         }

       /* If this is a metal, create special rasdef records for it */
       if (metal == TRUE)
         {
           /* If we have the residue number, create a rasdef record */
           if (chain != '\0' && res_num[0] != '\0')
             {
               /* Create a r_rasdef entry */
               r_rasdef_ptr
                 = create_r_rasdef_entry(&first_r_rasdef_ptr,
                                         &last_r_rasdef_ptr,type,number,
                                         res_name,res_num,chain,colour,
                                         input_line,range,range,copy);

               /* Update counts */
               r_rasdef_ptr->number = number;
               r_rasdef_ptr->copy = copy;
               ndef++;

               /* Format the metal name */
               strcpy(metal_name,res_name);
               r_remove_leading_blanks(metal_name);

               /* Store metal name as ligand name */
               store_string(&(r_rasdef_ptr->ligand_name),metal_name);

               /* Add the residue number and chain id */
               strcpy(res_name,metal_name);
               res_no = atoi(res_num);
               sprintf(metal_name,"%s %d(%c)",res_name,res_no,chain);

               /* Store as the ligand range */
               store_string(&(r_rasdef_ptr->ligand_range),metal_name);

               /* Store the LIGPLOT id */
               if (ligplot_id[0] != '\0')
                 store_string(&(r_rasdef_ptr->ligplot_id),ligplot_id);
             }
         }
 
       /* If search required, loop through the rasdef records to
          annotate the ones belonging to the current range */
       else if (search == TRUE)
         {
           /* Loop twice - first to locate the LIG_RANGE records, and
              then update all the component ligand residues */
           for (loop = 0; loop < 2; loop++)
             {
               /* Initialise flag */
               done = FALSE;
               
               /* Get pointer to the first RasMol definition record */
               r_rasdef_ptr = first_r_rasdef_ptr;
               
               /* Loop through all the definition records */
               while (r_rasdef_ptr != NULL && done == FALSE)
                 {
                   /* If this is the correct ligand-range record,
                      then get the ligand number and copy number */
                   if ((r_rasdef_ptr->type == R_LIGAND_RANGE &&
			!strcmp(r_rasdef_ptr->ligand_range,range)) ||
                       (r_rasdef_ptr->type == R_METAL &&
			!strcmp(r_rasdef_ptr->res_name,res_name)))
                     {
                       /* Update this record's ligand details */
                       r_rasdef_ptr->number = number;
                       r_rasdef_ptr->copy = copy;

                       /* Save the ligand number */
                       old_number = r_rasdef_ptr->number;

                       /* Store the LIGPLOT id */
                       if (ligplot_id[0] != '\0')
                         store_string(&(r_rasdef_ptr->ligplot_id),
                                      ligplot_id);

                       /* Set flag that we are done */
                       done = TRUE;
                     }

                   /* If this is a ligand belonging to our ligand
                      range, then update its numbers */
                   if (loop == 1 && r_rasdef_ptr->type == R_LIGAND &&
                       r_rasdef_ptr->number == old_number &&
                       r_rasdef_ptr->ligplot_id == &null)
                     {
                       /* Update this record's ligand details */
                       r_rasdef_ptr->number = number;
                       r_rasdef_ptr->copy = copy;

                       /* Store the LIGPLOT id */
                       if (ligplot_id[0] != '\0')
                         store_string(&(r_rasdef_ptr->ligplot_id),
                                      ligplot_id);
                     }

                   /* Get next definition record */
                   r_rasdef_ptr = r_rasdef_ptr->next_r_rasdef_ptr;
                 }
             }
        }
    }

  /* Return current pointers */
  *fst_r_ligand_ptr = first_r_ligand_ptr;
  *lst_r_ligand_ptr = last_r_ligand_ptr;
  *fst_r_rasdef_ptr = first_r_rasdef_ptr;
  *lst_r_rasdef_ptr = last_r_rasdef_ptr;

  /* Return number of metals */
  *nmet = nmetal;

  /* Close the file */
  fclose(file_ptr);

  /* Return number of new definition records created */
  return(ndef);
}
/***********************************************************************

r_get_dfn_data  -  Read in and store the RasMol definitions from the
                   rasmol.dfn file

***********************************************************************/

int r_get_dfn_data(char *pdb_code,int pdb_type,char *pdbsum_template,
                   struct r_rasdef **fst_r_rasdef_ptr,
                   struct r_ligand **fst_r_ligand_ptr,int *nlig,
                   int *nmet)
{
  char input_line[R_LINELEN + 1], file_name[R_FILENAME_LEN];
  char chain, metal_name[4], res_name[4], res_num[6];
  char col_name[R_COLNAM_LEN], ligand_name[R_LINELEN];
  char ligand_from[7], ligand_to[7], number_string[20];

  char *pdbsum_dir, *string_ptr;

  int count, i, line, nchain, ndef, nligand, nmetal, number;
  int colour, type;
  int first, metal, ok;

  struct r_chain *r_chain_ptr;
  struct r_ligand *r_ligand_ptr, *first_r_ligand_ptr, *last_r_ligand_ptr;
  struct r_rasdef *r_rasdef_ptr;
  struct r_rasdef *first_r_rasdef_ptr, *last_r_rasdef_ptr;

  FILE *file_ptr;

  /* Initialise */
  first = TRUE;
  count = 0;
  ligand_from[0] = '\0';
  ligand_to[0] = '\0';
  ligand_name[0] = '\0';
  line = 0;
  metal = FALSE;
  nchain = 0;
  ndef = 0;
  nligand = 0;
  nmetal = 0;
  ok = FALSE;

  /* Initialise r_rasdef pointers */
  r_chain_ptr = NULL;
  r_rasdef_ptr = NULL;
  *fst_r_rasdef_ptr = first_r_rasdef_ptr = last_r_rasdef_ptr = NULL;

  /* Initialise ligand pointers */
  r_ligand_ptr = first_r_ligand_ptr = last_r_ligand_ptr = NULL;

  /* Form the file name of the rasmol.dfn file */
  pdbsum_dir = form_filename(pdb_code,pdbsum_template);
  strcpy(file_name,pdbsum_dir);
  if (pdb_type == PDBSUM1)
    strcat(file_name,"/newsec");
  strcat(file_name,"/rasmol.dfn");

  /* Open the rasmol.dfn input file */
  if ((file_ptr = fopen(file_name,"r")) !=  NULL)
    {
      /* Set flag */
      ok = TRUE;

      /* Loop through the file processing each record */
      while (fgets(input_line,R_LINELEN,file_ptr)!= NULL)
        {
          /* Increment line-count */
          line++;

          /* Truncate the string to strip off the newline */
          string_chop(input_line);

          /* Check line-type according to 1st character: C, L, M,
             N or P */

          /* If protein CA-only chain, then flag as such */
          if (input_line[0] == 'C')
            {
              /* Define type */
              type = R_CALPHA_ONLY;
              nchain++;

              /* Define residue details */
              res_name[0] = '\0';
              res_num[0] = '\0';
              chain = input_line[1];
              if (chain == '\0' || chain == '\n')
                chain = ' ';
              col_name[0] = '\0';

              /* Create a r_rasdef entry */
              r_rasdef_ptr
                = create_r_rasdef_entry(&first_r_rasdef_ptr,
                                        &last_r_rasdef_ptr,type,nchain,
                                        res_name,res_num,chain,col_name,
                                        input_line,ligand_from,ligand_to,
                                        count);
              ndef++;
            }

          /* Check for ligand identifier */
          else if (input_line[0] == 'L')
            {
              /* Increment ligand-count */
              nligand++;

              /* Set flag that not metal */
              metal = FALSE;
              
              /* Create ligand record */
              r_ligand_ptr
                = create_r_ligand_record(&first_r_ligand_ptr,
                                         &last_r_ligand_ptr,r_chain_ptr,
                                         nligand,metal);

              /* Save the residue list */
              store_string(&(r_ligand_ptr->residue_list),input_line + 5);

              /* Extract the ligand residues */
              ndef
                = r_get_ligand_residues(input_line+4,&first_r_rasdef_ptr,
                                        &last_r_rasdef_ptr,nligand,ndef,
                                        r_ligand_ptr,ligand_name);
            }

          /* Check for ligand range */
          else if (input_line[0] == 'l')
            {
              /* Extract the number */
              strncpy(number_string,input_line + 1,3);
              number_string[3] = '\0';
              number = atoi(number_string);

              /* Define type as a range of residues representing a ligand */
              type = R_LIGAND_RANGE;

              /* Blank out any unwanted data */
              chain = ' ';
              res_name[0] = '\0';
              res_num[0] = '\0';
              col_name[0] = '\0';
              ligand_from[0] = '\0';
              ligand_to[0] = '\0';

              /* Get the chain */
              string_ptr = strchr(input_line,'(');
              if (string_ptr != NULL)
                chain = string_ptr[1];
              
              /* Create a r_rasdef entry for this ligand range */
              r_rasdef_ptr
                = create_r_rasdef_entry(&first_r_rasdef_ptr,
                                        &last_r_rasdef_ptr,type,number,
                                        res_name,res_num,chain,col_name,
                                        input_line,ligand_from,
                                        ligand_to,count);

              /* Increment count of definitions */
              ndef++;

              /* Save the ligand range */
              store_string(&(r_rasdef_ptr->ligand_range),input_line + 5);

              /* Save the ligand range */
              store_string(&(r_ligand_ptr->residue_range),input_line + 5);

              /* Store the corresponding link to the rasdef range record */
              r_ligand_ptr->r_rasdef_ptr = r_rasdef_ptr;

              /* Save the ligand name */
              if (ligand_name[0] != '\0')
                store_string(&(r_rasdef_ptr->ligand_name),ligand_name);
            }

          /* Check for metal */
          else if (input_line[0] == 'M')
            {
              /* Define type */
              type = R_METAL;
              nmetal++;

              /* Extract the metal name into the residue number field */
              strncpy(res_num,input_line+1,3);
              res_num[3] = ' ';
              res_num[4] = '\0';

              /* Repeat for the residue name field */
              strncpy(res_name,input_line,3);
              res_name[0] = ' ';
              res_name[3] = '\0';

              /* Get the metal colour */
              if (input_line[4] == ' ')
                strcpy(col_name,input_line+5);
              else
                strcpy(col_name,input_line+4);

              /* Create a r_rasdef entry */
              r_rasdef_ptr
                = create_r_rasdef_entry(&first_r_rasdef_ptr,
                                      &last_r_rasdef_ptr,type,nmetal,
                                      res_name,res_num,chain,col_name,
                                      input_line,ligand_from,ligand_to,
                                      count);
              ndef++;

              /* Format the metal name */
              strcpy(metal_name,res_name);
              r_remove_leading_blanks(metal_name);
              if (strlen(metal_name) > 1 && metal_name[0] != ' ')
                r_to_lower(metal_name+1,strlen(metal_name+1));

              /* Store metal name as ligand name */
              store_string(&(r_rasdef_ptr->ligand_name),metal_name);

              /* Set flag that metal */
              metal = TRUE;
              
              /* Create ligand record */
              r_ligand_ptr
                = create_r_ligand_record(&first_r_ligand_ptr,
                                         &last_r_ligand_ptr,r_chain_ptr,
                                         nmetal,metal);

              /* Save the metal name */
              strcpy(r_ligand_ptr->metal_name,metal_name);

              /* Get the metal colour */
              colour = 0;
              for (i = 0; i < R_NMETAL_COLOURS; i++)
                if (!strcmp(col_name,R_METAL_COLOUR[i]))
                  colour = i;

              /* Save the metal colour */
              r_ligand_ptr->metal_colour = colour;

              /* Store the corresponding link to the rasdef range record */
              r_ligand_ptr->r_rasdef_ptr = r_rasdef_ptr;
            }
          
          /* Check for nucleic acid chain */
          else if (input_line[0] == 'N')
            {
              /* Define type */
              type = R_NUCLEIC_ACID;
              nchain++;

              /* Define residue details */
              res_name[0] = '\0';
              res_num[0] = '\0';
              chain = input_line[1];
              if (chain == '\0' || chain == '\n')
                chain = ' ';
              col_name[0] = '\0';

              /* Create a r_rasdef entry */
              r_rasdef_ptr
                = create_r_rasdef_entry(&first_r_rasdef_ptr,
                                      &last_r_rasdef_ptr,type,nchain,
                                      res_name,res_num,chain,col_name,
                                      input_line,ligand_from,ligand_to,
                                      count);
              ndef++;

              /* Store "DNA/RNA" as the ligand name */
              store_string(&(r_rasdef_ptr->ligand_name),"DNA/RNA");
            }

          /* If protein chain, then write the script records to define
             and display the chain */
          else if (input_line[0] == 'P')
            {
              /* Define type */
              type = R_PROTEIN;
              nchain++;

              /* Define residue details */
              res_name[0] = '\0';
              res_num[0] = '\0';
              chain = input_line[1];
              if (chain == '\0' || chain == '\n')
                chain = ' ';
              col_name[0] = '\0';

              /* Create a r_rasdef entry */
              r_rasdef_ptr
                = create_r_rasdef_entry(&first_r_rasdef_ptr,
                                      &last_r_rasdef_ptr,type,nchain,
                                      res_name,res_num,chain,col_name,
                                      input_line,ligand_from,ligand_to,
                                      count);
              ndef++;
            }
        }

      /* Close the r_rasdef file */
      fclose(file_ptr);
    }

  /* If no rasmol.dfn file, then can't continue */
  else
    {
      /* Show error message and stop */
      printf("*** Warning. Unable to open file: %s\n",file_name);
      ok = FALSE;
      return(ok);
    }

  /* If we have any ligands or metals, process them */
  if (nligand > 0 || nmetal > 0)
    {
      /* Identify and number unique ligands and copies */
      r_number_ligands(first_r_rasdef_ptr);

      /* Map to corresponding LIGPLOT diagrams in PDBsum */
      map_to_ligplots(pdbsum_dir,&first_r_rasdef_ptr,&last_r_rasdef_ptr,
                      &first_r_ligand_ptr,&last_r_ligand_ptr,&nmetal);
    }
  
  /* Return pointer to first RasMol definition */
  *fst_r_rasdef_ptr = first_r_rasdef_ptr;

  /* Return pointer to the first ligand */
  *fst_r_ligand_ptr = first_r_ligand_ptr;

  /* Return total numbers of ligand and metals */
  *nlig = nligand;
  *nmet = nmetal;

  /* Return flag */
  return(ok);
}
/***********************************************************************

create_r_pdb_entry  -  Create a PDB entry

***********************************************************************/

struct r_pdb *create_r_pdb_entry(char *pdb_code)
{
  struct r_pdb *r_pdb_ptr;

  /* Initialise r_pdb pointers */
  r_pdb_ptr = NULL;

  /* Create r_pdb entry */
  r_pdb_ptr = (struct r_pdb *) malloc(sizeof(struct r_pdb));
  if (r_pdb_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for struct r_pdb\n");
      printf("***        Program terminated with error.\n");
      exit (1);
    }

  /* Store the R_PDB code */
  store_string(&(r_pdb_ptr->pdb_code),pdb_code);

  /* Initialise remaining fields */
  r_pdb_ptr->protein_name = &null;
  r_pdb_ptr->ninfo = 0;
  r_pdb_ptr->nmodels = 0;
  r_pdb_ptr->nchains = 0;
  r_pdb_ptr->nmolecs = 0;
  r_pdb_ptr->nligands = 0;
  r_pdb_ptr->nmetals = 0;
  r_pdb_ptr->nresid = 0;
  r_pdb_ptr->natoms = 0;
  r_pdb_ptr->first_r_info_ptr = NULL;
  r_pdb_ptr->first_r_ref_ptr = NULL;
  r_pdb_ptr->first_r_hetdesc_ptr = NULL;
  r_pdb_ptr->first_r_chain_ptr = NULL;
  r_pdb_ptr->first_r_molec_ptr = NULL;
  r_pdb_ptr->first_r_ligand_ptr = NULL;
  r_pdb_ptr->first_r_residue_ptr = NULL;
  r_pdb_ptr->first_r_atom_ptr = NULL;
  r_pdb_ptr->last_r_info_ptr = NULL;
  r_pdb_ptr->last_r_ref_ptr = NULL;
  r_pdb_ptr->last_r_chain_ptr = NULL;
  r_pdb_ptr->last_r_residue_ptr = NULL;
  r_pdb_ptr->last_r_atom_ptr = NULL;

  /* Return current pointers */
  return(r_pdb_ptr);
}
/***********************************************************************

create_r_info_record  -  Create a new data record for an item of PDB
                         header information

***********************************************************************/

struct r_info *create_r_info_record(struct r_info **fst_r_info_ptr,
                                    struct r_info **lst_r_info_ptr,
                                    char *name,char *value)
{
  int level;

  struct r_info *r_info_ptr;
  struct r_info *first_r_info_ptr, *last_r_info_ptr;

  /* Initialise data pointers */
  r_info_ptr = NULL;
  first_r_info_ptr = *fst_r_info_ptr;
  last_r_info_ptr = *lst_r_info_ptr;

  /* Create data entry */
  r_info_ptr = (struct r_info*)malloc(sizeof(struct r_info));
  if (r_info_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for");
      printf(" struct r_info\n");
      printf("***        Program terminated with error.\n");
      exit (1);
    }

  /* Store the known data items */
  store_string(&(r_info_ptr->name),name);
  store_string(&(r_info_ptr->value),value);

  /* Initialise remaining fields */
  r_info_ptr->nlevel = 0;
  for (level = 0; level < R_MAX_LEVEL; level++)
    r_info_ptr->instance[level] = 0;

  /* Initialise pointer to next data record */
  r_info_ptr->next_r_info_ptr = NULL;

  /* If this is the first entry stored, then update pointer to head
     of linked list */
  if (first_r_info_ptr == NULL)
    first_r_info_ptr = r_info_ptr;

  /* Otherwise, make previous entry point to the current one */
  else
    last_r_info_ptr->next_r_info_ptr = r_info_ptr;
                      
  /* Save the current data's pointer */
  last_r_info_ptr = r_info_ptr;

  /* Return current pointers */
  *fst_r_info_ptr = first_r_info_ptr;
  *lst_r_info_ptr = last_r_info_ptr;
  return(r_info_ptr);
}
/***********************************************************************

r_assign_to_object  -  Check which of the elements in the rasmol.dfn file
                       the given residue should be assigned to

***********************************************************************/

struct r_rasdef *r_assign_to_object(struct r_rasdef *first_r_rasdef_ptr,
                                    char *atom_name,char *res_name,
                                    char *res_num,char chain)
{
  int type;
  int have_exact;

  struct r_rasdef *r_rasdef_ptr, *assigned_r_rasdef_ptr;

  /* Initialise variables */
  assigned_r_rasdef_ptr = NULL;
  have_exact = FALSE;
  type = R_UNKNOWN;

  /* If this is a water, then return */
  if (!strcmp(res_name,"HOH"))
    return(assigned_r_rasdef_ptr);

  /* Get pointer to the first RasMol definition record */
  r_rasdef_ptr = first_r_rasdef_ptr;

  /* Loop through all the definition records */
  while (r_rasdef_ptr != NULL && have_exact == FALSE)
    {
      /* Check for a match to a metal */
      if (r_rasdef_ptr->type == R_METAL)
        {
          /* Check that atom name matches that stored in the
             residue number */
          if (!strncmp(atom_name,r_rasdef_ptr->res_num,4))
            {
              /* Store the pointer and set flag */
              have_exact = TRUE;
              assigned_r_rasdef_ptr = r_rasdef_ptr;
            }
        }

      /* Check for exact match */
      else if (!strcmp(r_rasdef_ptr->res_name,res_name) &&
               !strcmp(r_rasdef_ptr->res_num,res_num) &&
               r_rasdef_ptr->chain == chain)
        {
          /* Store the pointer and set flag */
          have_exact = TRUE;
          assigned_r_rasdef_ptr = r_rasdef_ptr;
        }

      /* Otherwise, check for a match on the chain-id only */
      else if (r_rasdef_ptr->chain == chain &&
               r_rasdef_ptr->res_name[0] == '\0' &&
               r_rasdef_ptr->type != R_LIGAND_RANGE)
        assigned_r_rasdef_ptr = r_rasdef_ptr;

      /* Get next definition record */
      r_rasdef_ptr = r_rasdef_ptr->next_r_rasdef_ptr;
    }

  /* Return the matched rasdef pointer */
  return(assigned_r_rasdef_ptr);
}
/***********************************************************************

r_process_line  -  Process the current PDB header line

***********************************************************************/

int r_process_line(char *input_line,char *string1,char *string2,
                   char *string3,char *string4,char *string5,
                   char *string6,int *synth,int keyword,int molecule,
                   int *protein_no,int *np,char *reference,int *flag,
                   struct r_rasdef *first_r_rasdef_ptr)
{
  char atom_name[5], chain, res_name[4], res_num[6], tmp[R_LINELEN + 1];
  char number_string[20];

  char *string_ptr;

  int imolec, len, nprot, number;
  int is_protein, synthetic;

  struct r_rasdef *r_rasdef_ptr;

  /* Initialise variables */
  chain = ' ';
  atom_name[0] = '\0';
  nprot = *np;
  number_string[0] = '\0';
  res_name[0] = '\0';
  res_num[0] = '\0';
  synthetic = *synth;

  /* Check for continuation character */
  if (input_line[3] == ' ')
    strcpy(tmp,input_line+4);
  else
    strcpy(tmp,input_line+5);
  strcpy(input_line,tmp);

  /* Remove leading blanks */
  if (input_line[0] != '\0' && input_line[0] == ' ')
    remove_leading_blanks(input_line);

  /* Truncate the string after last non-blank character */
  len = string_chop(input_line);

  /* If this is the molecule id line, get the number */
  if (!strncmp(input_line,"MOL_ID: ",8))
    {
      /* Get the molecule number */
      strcpy(number_string,input_line+8);
      string_ptr = strstr(number_string,";");
      if (string_ptr != NULL)
        string_ptr[0] = '\0';

      /* Get number */
      molecule = atoi(number_string);
    }

  /* Process line according to current keyword */
  switch (keyword)
    {
      /* a. Title */
    case TITLE :

      /* If no title yet, copy current title across */
      if (string1[0] == '\0')
        strcat(string1,input_line);

      /* Otherwise, add onto end of existing string, if have room */
      else if (len + strlen(string1) < R_STRING_LEN - 2)
        {
          /* Add space and then current line */
          strcat(string1," ");
          strcat(string1,input_line);
        }

      break;

      /* b. Compound info */
    case COMPOUND :

      /* Check for molecule name */
      if (!strncmp(input_line,"MOLECULE: ",10))
        {
          /* Save the molecule name */
          strcpy(string1,input_line+10);

          /* Set flag that we're in the molecule name */
          flag[keyword] = MOLECULE_NAME;
        }

      /* Check for list of synonyms */
      else if (!strncmp(input_line,"SYNONYM: ",9))
        {
          /* Save the synonym list */
          strcpy(string2,input_line+9);

          /* Set flag that we're in the synonyms */
          flag[keyword] = SYNONYMS;
        }

      /* Check for chain identifier(s) */
      else if (!strncmp(input_line,"CHAIN: ",7))
        {
          /* Extract the chain identifier(s) */
          strcpy(tmp,input_line+7);
          strcpy(input_line,tmp);

          /* Check for blank chain-id */
          if (!strcmp(input_line,"NULL;"))
            strcpy(input_line," ;");

          /* Check first chain to determine whether this molecule
             is a protein or not */
          chain = input_line[0];
          r_rasdef_ptr = r_assign_to_object(first_r_rasdef_ptr,atom_name,
                                            res_name,res_num,chain);

          /* If we don't have the rasmol.dfn data, then use the molecule's
             name instead */
          is_protein = FALSE;
          if (r_rasdef_ptr == NULL)
            {
              if (!strncmp(string1,"DNA;",4) ||
                  !strncmp(string1,"RNA;",4) ||
                  !strncmp(string1,"DNA (",5) ||
                  !strncmp(string1,"RNA (",5) || 
                  !strncmp(string1,"DNA/RNA",7) ||
                  strstr(string1," DNA;") != NULL ||
                  strstr(string1," RNA;") != NULL)
                is_protein = FALSE;
              else
                is_protein = TRUE;
            }

          /* If this is a protein chain, set flag */
          if ((r_rasdef_ptr != NULL &&
              (r_rasdef_ptr->type == R_PROTEIN ||
               r_rasdef_ptr->type == R_CALPHA_ONLY)) || is_protein == TRUE)
            {
              /* If have a valid molecule number, set flag */
              if (molecule > 0 && molecule < R_MAX_MOLECS)
                {
                  /* Increment count of proteins and store number */
                  nprot++;
                  protein_no[molecule - 1] = nprot;
                }
            }

          /* Unset flag */
          flag[keyword] = NONE;
        }

      /* Check whether this is some other keyword */
      else
        {
          /* Look for colon, which signifies a key word */
          string_ptr = strstr(input_line,":");

          /* If no keyword, then this might be a continuation of a
             previous term */
          if (string_ptr == NULL)
            {
              /* If previous keyword is the molecule name or list of
                 synonyms, then add to current string */
              if (flag[keyword] == MOLECULE_NAME)
                {
                  /* Add if not too long */
                  if (len + strlen(string1) < R_STRING_LEN - 2)
                    {
                      /* Add space and then current line */
                      strcat(string1," ");
                      strcat(string1,input_line);
                    }
                }

              /* Check for coninuation of list of synonyms */
              else if (flag[keyword] == SYNONYMS)
                {
                  /* Add if not too long */
                  if (len + strlen(string2) < R_STRING_LEN - 2)
                    {
                      /* Add space and then current line */
                      strcat(string2," ");
                      strcat(string2,input_line);
                    }
                }
            }

          /* Otherwise, reset flag */
          else
            flag[keyword] = NONE;
        }

      break;

      /* c. Species name */
    case SOURCE :

      /* Check for scientific name */
      if (!strncmp(input_line,"ORGANISM_SCIENTIFIC: ",21))
        {
          /* Save the name */
          strcpy(string1,input_line+21);

          /* Set flag that we're in the scientific name */
          flag[keyword] = SCIENTIFIC_NAME;
        }

      /* Check for common species name */
      else if (!strncmp(input_line,"ORGANISM_COMMON: ",17))
        {
          /* Save the common species name  */
          strcpy(string2,input_line+17);

          /* Set flag that we're in the common name */
          flag[keyword] = COMMON_NAME;
        }

      /* Check for common species name */
      else if (!strncmp(input_line,"ORGANISM_TAXID: ",16))
        {
          /* Save the taxonomy id  */
          strcpy(string3,input_line+16);

          /* Set flag that we're in the common name */
          flag[keyword] = TAXONOMY_ID;
        }

      /* Special case of synthetic constructs */
      else if (!strncmp(input_line,"SYNTHETIC: YES",14))
        synthetic = TRUE;

      /* Check whether this is some other keyword */
      else
        {
          /* Look for colon, which signifies a key word */
          string_ptr = strstr(input_line,":");

          /* If no keyword, then this might be a continuation of a
             previous term */
          if (string_ptr == NULL)
            {
              /* If previous keyword is the scientific name, then add
                 to current string */
              if (flag[keyword] == SCIENTIFIC_NAME)
                {
                  /* Add if not too long */
                  if (len + strlen(string1) < R_STRING_LEN - 2)
                    {
                      /* Add space and then current line */
                      strcat(string1," ");
                      strcat(string1,input_line);
                    }
                }

              /* Check for continuation of common name */
              else if (flag[keyword] == COMMON_NAME)
                {
                  /* Add if not too long */
                  if (len + strlen(string2) < R_STRING_LEN - 2)
                    {
                      /* Add space and then current line */
                      strcat(string2," ");
                      strcat(string2,input_line);
                    }
                }
            }

          /* Otherwise, reset flag */
          else
            flag[keyword] = NONE;
        }

      break;

      /* d. Key words */
    case KEY_WORDS :

      /* If no key words yet, copy current line across */
      if (string1[0] == '\0')
        strcat(string1,input_line);

      /* Otherwise, add onto end of existing string, if have room */
      else if (len + strlen(string1) < R_STRING_LEN - 2)
        {
          /* Add space and then current line */
          strcat(string1," ");
          strcat(string1,input_line);
        }

      break;

      /* e. Authors */
    case AUTHOR :

      /* If no authors yet, copy current line across */
      if (string1[0] == '\0')
        strcat(string1,input_line);

      /* Otherwise, add onto end of existing string, if have room */
      else if (len + strlen(string1) < R_STRING_LEN - 2)
        {
          /* Add current line */
          strcat(string1,input_line);
        }

      break;

      /* g. Additional reference */
    case REFERENCE :

      /* If this is the reference identifier, determine reference number */
      if (!strncmp(input_line,"REFERENCE ",10))
        {
          /* Extract the number */
          strcpy(number_string,input_line+10);
          number = atoi(number_string);
          sprintf(reference,"%d",number);
        }

      /* f. Key reference */
    case JOURNAL :

      /* Author line */
      if (!strncmp(input_line,"AUTH ",5))
        {
          /* If no author list yet, copy current line across */
          if (string1[0] == '\0')
            strcat(string1,input_line+7);

          /* Otherwise, add onto end of existing string, if have room */
          else if (len + strlen(string1) < R_STRING_LEN - 2)
            {
              /* Add space and then current line */
              strcat(string1,input_line+7);
            }
        }

      /* Title line */
      else if (!strncmp(input_line,"TITL ",5))
        {
          /* If no title yet, copy current line across */
          if (string2[0] == '\0')
            strcat(string2,input_line+7);

          /* Otherwise, add onto end of existing string, if have room */
          else if (len + strlen(string2) < R_STRING_LEN - 2)
            {
              /* Add space and then current line */
              strcat(string2," ");
              strcat(string2,input_line+7);
            }
        }

      /* Reference */
      else if (!strncmp(input_line,"REF  ",5))
        {
          /* If no reference yet, copy current line across */
          if (string3[0] == '\0')
            strcat(string3,input_line+7);

          /* Otherwise, add onto end of existing string, if have room */
          else if (len + strlen(string3) < R_STRING_LEN - 2)
            {
              /* Add space and then current line */
              strcat(string3," ");
              strcat(string3,input_line+7);
            }
        }

      /* Reference number */
      else if (!strncmp(input_line,"REFN ",5))
        {
          /* If no reference yet, copy current line across */
          if (string4[0] == '\0')
            strcat(string4,input_line+7);

          /* Otherwise, add onto end of existing string, if have room */
          else if (len + strlen(string4) < R_STRING_LEN - 2)
            {
              /* Add space and then current line */
              strcat(string4," ");
              strcat(string4,input_line+7);
            }
        }

      /* PubMed id */
      else if (!strncmp(input_line,"PMID ",5))
        {
          /* If no id yet, copy current line across */
          if (string5[0] == '\0')
            strcat(string5,input_line+7);
        }

      /* DOI number */
      else if (!strncmp(input_line,"DOI  ",5))
        {
          /* If no reference yet, copy current line across */
          if (string6[0] == '\0')
            strcat(string6,input_line+7);
        }

      break;

    default :
      break;
    }

  /* Return current number of proteins */
  *np = nprot;

  /* Return the synthetic flag */
  *synth = synthetic;

  /* Return current molecule number */
  return(molecule);
}
/***********************************************************************

create_r_ref_record  -  Create a reference entry

***********************************************************************/

struct r_ref *create_r_ref_record(struct r_ref **fst_r_ref_ptr,
                                  struct r_ref **lst_r_ref_ptr,
                                  int type,char *number)
{
  int i;

  struct r_ref *r_ref_ptr;
  struct r_ref *first_r_ref_ptr, *last_r_ref_ptr;

  /* Initialise pointers */
  r_ref_ptr = NULL;
  first_r_ref_ptr = *fst_r_ref_ptr;
  last_r_ref_ptr = *lst_r_ref_ptr;

  /* Create ref entry */
  r_ref_ptr = (struct r_ref *) malloc(sizeof(struct r_ref));
  if (r_ref_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for struct r_ref\n");
      printf("***        Program terminated with error.\n");
      exit (1);
    }

  /* Store the known data */
  r_ref_ptr->type = type;
  strcpy(r_ref_ptr->number,number);

  /* Initialise other fields */
  r_ref_ptr->nauthors = 0;
  for (i = 0; i < R_MAX_AUTHORS; i++)
    r_ref_ptr->author[i] = &null;
  r_ref_ptr->title = &null;
  r_ref_ptr->journal = &null;
  r_ref_ptr->doi_number = &null;
  r_ref_ptr->year = 0;
  r_ref_ptr->volume = 0;
  r_ref_ptr->start_page = &null;
  r_ref_ptr->end_page = &null;
  r_ref_ptr->to_be_published = FALSE;
  r_ref_ptr->new = FALSE;

  /* Initialise pointer to next entry */
  r_ref_ptr->next_r_ref_ptr = NULL;

  /* If this is the first entry stored, then update pointer to head
     of linked list */
  if (first_r_ref_ptr == NULL)
    first_r_ref_ptr = r_ref_ptr;

  /* Otherwise, make previous entry point to the current one */
  else
    last_r_ref_ptr->next_r_ref_ptr = r_ref_ptr;
                      
  /* Save the current residue's pointer */
  last_r_ref_ptr = r_ref_ptr;

  /* Return current pointers */
  *fst_r_ref_ptr = first_r_ref_ptr;
  *lst_r_ref_ptr = last_r_ref_ptr;
  return(r_ref_ptr);
}
/***********************************************************************

r_get_reference  -  Break up the journal reference into its components

***********************************************************************/

void r_get_reference(char *string,char *journal,int *year,int *volume,
                     char *page)
{
  char ch, number_string[20];

  char *string_ptr, *tmp;

  int i, ipos, len;

  /* Initialise numbers */
  page[0] = '\0';
  *year = 0;
  *volume = 0;

  /* Get the length of the string */
  len = strlen(string);
  if (len == 0)
    return;

  /* Allocate memory for temporary string */
  tmp = (char *) malloc(sizeof(char)*(len + 1));

  /* If article not marked as "to be published" then extract info */
  if (strcmp(string,"to be published"))
    {
      /* Get just the journal name */
      strcpy(journal,string);
      string_ptr = strstr(journal," v.");
      if (string_ptr != NULL)
        string_ptr[0] = '\0';

      /* If details not given, get just the journal name */
      else
        journal[30] = '\0';

      /* Remove trailing blank spaces */
      len = string_chop(journal);

      /* Remove any funny characters */
      remove_marks(journal);

      /* Format the journal name */
      type_set(journal,TRUE);

      /* Retrieve the volume, page and year */
      string_ptr = strstr(string," v.");
      if (string_ptr != NULL)
        {
          /* Shift string along */
          strcpy(tmp,string_ptr+3);
          strcpy(string,tmp);

          /* Get the volume, check for cases where first character is
             a letter rather than a number (eg in the Acta Cryst jnls) */
          strncpy(number_string,string,4);
          number_string[4] = '\0';
          ch = ' ';
          ipos = 0;
          for (i = 0; i < 3 && ch == ' '; i++)
            {
              /* Get this character */
              ch = number_string[i];

              /* If first non-blank character, check if letter or number */
              if (ch != ' ')
                {
                  /* If not number set number position one higher */
                  if (ch < '0' || ch > '9')
                    ipos = i + 1;
                }
            }

          /* Get the voulme number */
          *volume = atoi(number_string+ipos);

          /* Get the page */
          strncpy(number_string,string+4,6);
          number_string[6] = '\0';
          strcpy(page,number_string);

          /* Get the year */
          strncpy(number_string,string+10,5);
          number_string[5] = '\0';
          *year = atoi(number_string);
        }
    }

  /* Copy "TO BE PUBLISHED" to journal name */
  else
    strcpy(journal,"TO BE PUBLISHED");

  /* Free up memory taken by the temporary string */
  free(tmp);
}
/***********************************************************************

r_store_data  -  Store the current PDB header item

***********************************************************************/

void r_store_data(int keyword,char *string1,char *string2,char *string3,
                  char *string4,char *string5,char *string6,
                  int *protein_no,int molecule,char *reference,
                  struct r_info **fst_r_info_ptr,
                  struct r_info **lst_r_info_ptr,
                  struct r_ref **fst_r_ref_ptr,
                  struct r_ref **lst_r_ref_ptr)
{
  char journal[R_LINELEN+1];
  char name[R_TOKEN_LEN], number_string[20], page[20], tmp[R_TOKEN_LEN];

  char *string_ptr, *token[R_MAX_TOKEN];
  int itoken, ntokens, type, volume, year;

  struct r_info *first_r_info_ptr, *last_r_info_ptr, *r_info_ptr;
  struct r_ref *first_r_ref_ptr, *last_r_ref_ptr, *r_ref_ptr;

  /* Initialise pointers */
  first_r_info_ptr = *fst_r_info_ptr;
  last_r_info_ptr = *lst_r_info_ptr;
  first_r_ref_ptr = *fst_r_ref_ptr;
  last_r_ref_ptr = *lst_r_ref_ptr;

  /* Convert strings to lower case */
  to_lower(string1,strlen(string1));
  to_lower(string2,strlen(string2));
  to_lower(string3,strlen(string3));
  to_lower(string4,strlen(string4));

  /* Remove any semi-colons */
  string_ptr = strstr(string1,";");
  if (string_ptr != NULL)
    string_ptr[0] = '\0';
  string_ptr = strstr(string2,";");
  if (string_ptr != NULL)
    string_ptr[0] = '\0';
  string_ptr = strstr(string3,";");
  if (string_ptr != NULL)
    string_ptr[0] = '\0';
  string_ptr = strstr(string4,";");
  if (string_ptr != NULL)
    string_ptr[0] = '\0';

  /* Store according to current keyword */
  switch (keyword)
    {
      /* a. Title */
    case TITLE :

      /* Format the title */
      type_set(string1,TRUE);
      fix_upper(string1);

      /* Store title */
      r_info_ptr
        = create_r_info_record(&first_r_info_ptr,&last_r_info_ptr,
                               "TITLE",string1);

      break;

      /* b. Compound */
    case COMPOUND :

      /* If this is a protein molecule, store */
      if (molecule > 0 && protein_no[molecule - 1] > 0)
        {
          /* Store the molecule name */
          r_info_ptr
            = create_r_info_record(&first_r_info_ptr,&last_r_info_ptr,
                                   "PROTEIN",string1);
          r_info_ptr->nlevel = 1;
          r_info_ptr->instance[0] = protein_no[molecule - 1];

          /* Break up the synonyms line into the individual synonyms */
          ntokens = get_tokens(string2,token,R_MAX_TOKEN,',');

          /* Store the number of synonyms */
          if (ntokens > 0)
            {
              sprintf(number_string,"%d",ntokens);
              r_info_ptr
                = create_r_info_record(&first_r_info_ptr,&last_r_info_ptr,
                                       "NSYNONYMS",number_string);
              r_info_ptr->nlevel = 1;
              r_info_ptr->instance[0] = protein_no[molecule - 1];
            }

          /* Loop through the tokens to store each one separately */
          for (itoken = 0; itoken < ntokens; itoken++)
            {
              /* Remove leading blank space, if there is one */
              if (token[itoken][0] == ' ')
                {
                  strcpy(tmp,token[itoken]+1);
                  strcpy(token[itoken],tmp);
                }

              /* Store this synonym */
              r_info_ptr
                = create_r_info_record(&first_r_info_ptr,&last_r_info_ptr,
                                       "SYNONYM",token[itoken]);
              r_info_ptr->nlevel = 2;
              r_info_ptr->instance[0] = protein_no[molecule - 1];
              r_info_ptr->instance[1] = itoken + 1;
            }
        }

      break;

      /* c. Species */
    case SOURCE :

      /* If this is a protein molecule, store species */
      if (molecule > 0 && protein_no[molecule - 1] > 0)
        {
          /* Store the organism name */
          r_info_ptr
            = create_r_info_record(&first_r_info_ptr,&last_r_info_ptr,
                                   "ORGANISM_SCIENTIFIC",string1);
          r_info_ptr->nlevel = 1;
          r_info_ptr->instance[0] = protein_no[molecule - 1];

          /* Store the common name */
          r_info_ptr
            = create_r_info_record(&first_r_info_ptr,&last_r_info_ptr,
                                   "ORGANISM_COMMON",string2);
          r_info_ptr->nlevel = 1;
          r_info_ptr->instance[0] = protein_no[molecule - 1];

          /* Store the taxonomy id */
          r_info_ptr
            = create_r_info_record(&first_r_info_ptr,&last_r_info_ptr,
                                   "ORGANISM_TAXID",string3);
          r_info_ptr->nlevel = 1;
          r_info_ptr->instance[0] = protein_no[molecule - 1];
        }

      break;

      /* d. Key words */
    case KEY_WORDS :

      /* Break up the key words line into the individual key words or
         phrases */
      ntokens = get_tokens(string1,token,R_MAX_TOKEN,',');

      /* Store the number of key words */
      if (ntokens > 0)
        {
          sprintf(number_string,"%d",ntokens);
          r_info_ptr
            = create_r_info_record(&first_r_info_ptr,&last_r_info_ptr,
                                   "NKEY_WORDS",number_string);
        }

      /* Loop through the tokens to store */
      for (itoken = 0; itoken < ntokens; itoken++)
        {
          /* Remove leading blank spaces */
          if (token[itoken][0] == ' ')
            remove_leading_blanks(token[itoken]);

          /* Store the key word */
          r_info_ptr
            = create_r_info_record(&first_r_info_ptr,&last_r_info_ptr,
                                   "KEY",token[itoken]);
          r_info_ptr->nlevel = 1;
          r_info_ptr->instance[0] = itoken + 1;
        }

      break;

      /* e. Authors */
    case AUTHOR :

      /* Break up the author line into the individual author names */
      ntokens = get_tokens(string1,token,R_MAX_TOKEN,',');

      /* Store the number of authors */
      if (ntokens > 0)
        {
          sprintf(number_string,"%d",ntokens);
          r_info_ptr
            = create_r_info_record(&first_r_info_ptr,&last_r_info_ptr,
                                   "NAUTHORS",number_string);
        }

      /* Loop through the tokens to store */
      for (itoken = 0; itoken < ntokens; itoken++)
        {
          /* Format the author name */
          strcpy(name,token[itoken]);
          type_set(name,TRUE);

          /* Write out this author */
          r_info_ptr
            = create_r_info_record(&first_r_info_ptr,&last_r_info_ptr,
                                   "AUTHOR",name);
          r_info_ptr->nlevel = 1;
          r_info_ptr->instance[0] = itoken + 1;
        }

      break;

      /* f. Key reference */
    case JOURNAL :

      /* Define reference type */
      type = KEY_REF;

      /* Break up the author line into the individual author names */
      ntokens = get_tokens(string1,token,R_MAX_TOKEN,',');

      /* Create a reference record and store the number of authors */
      r_ref_ptr
        = create_r_ref_record(&first_r_ref_ptr,&last_r_ref_ptr,type,"0");
      if (ntokens > R_MAX_AUTHORS)
        r_ref_ptr->nauthors = R_MAX_AUTHORS;
      else
        r_ref_ptr->nauthors = ntokens;

      /* Flag reference as a new one */
      r_ref_ptr->new = TRUE;

      /* Loop through the tokens to format and store author names */
      for (itoken = 0; itoken < ntokens && itoken < R_MAX_AUTHORS; itoken++)
        {
          /* Format the author name */
          strcpy(name,token[itoken]);
          type_set(name,TRUE);

          /* Store this author */
          store_string(&(r_ref_ptr->author[itoken]),name);
        }

      /* format the title */
      type_set(string2,TRUE);
      fix_upper(string2);

      /* Store the article title */
      store_string(&(r_ref_ptr->title),string2);

      /* Break up the journal reference into its components */
      r_get_reference(string3,journal,&year,&volume,page);

      /* Store the name of the journal */
      store_string(&(r_ref_ptr->journal),journal);

      /* If reference annotated as "TO BE PUBLISHED", set flag */
      if (!strcmp(journal,"TO BE PUBLISHED"))
        r_ref_ptr->to_be_published = TRUE;
      else
        r_ref_ptr->to_be_published = FALSE;

      /* If have the PubMed id, store that */
      if (string5[0] != '\0' && strlen(string5) < R_PUBMED_LEN)
        strcpy(r_ref_ptr->pubmed_id,string5);

      /* If have the DOI number, store that */
      if (string6[0] != '\0')
        store_string(&(r_ref_ptr->doi_number),string6);

      /* Store the volume, year and page */
      r_ref_ptr->year = year;
      r_ref_ptr->volume = volume;

      /* Remove leading blanks from page number and store */
      remove_leading_blanks(page);
      store_string(&(r_ref_ptr->start_page),page);

      break;

      /* g. Additional references */
    case REFERENCE :

      /* Define reference type */
      type = SECONDARY_REF;

      /* Break up the author line into the individual author names */
      ntokens = get_tokens(string1,token,R_MAX_TOKEN,',');

      /* Create a reference record and store the number of authors */
      r_ref_ptr = create_r_ref_record(&first_r_ref_ptr,&last_r_ref_ptr,type,
                                      reference);
      if (ntokens > R_MAX_AUTHORS)
        r_ref_ptr->nauthors = R_MAX_AUTHORS;
      else
        r_ref_ptr->nauthors = ntokens;

      /* Flag reference as a new one */
      r_ref_ptr->new = TRUE;

      /* Loop through the tokens to store author names */
      for (itoken = 0; itoken < ntokens && itoken < R_MAX_AUTHORS; itoken++)
        {
          /* Format the author name */
          strcpy(name,token[itoken]);
          type_set(name,TRUE);

          /* Store this author */
          store_string(&(r_ref_ptr->author[itoken]),name);
        }

      /* format the title */
      type_set(string2,TRUE);
      fix_upper(string2);

      /* Store the article title */
      store_string(&(r_ref_ptr->title),string2);

      /* Break up the journal reference into its components */
      r_get_reference(string3,journal,&year,&volume,page);

      /* Store the name of the journal */
      store_string(&(r_ref_ptr->journal),journal);

      /* If reference annotated as "TO BE PUBLISHED", set flag */
      if (!strcmp(journal,"TO BE PUBLISHED"))
        r_ref_ptr->to_be_published = TRUE;

      /* Store the volume, year and page */
      r_ref_ptr->year = year;
      r_ref_ptr->volume = volume;
      store_string(&(r_ref_ptr->start_page),page);

      break;

    default :
      break;
    }

  /* Reinitialise all strings */
  string1[0] = '\0';
  string2[0] = '\0';
  string3[0] = '\0';
  string4[0] = '\0';
  string5[0] = '\0';
  string6[0] = '\0';

  /* Return pointers */
  *fst_r_info_ptr = first_r_info_ptr;
  *lst_r_info_ptr = last_r_info_ptr;
  *fst_r_ref_ptr = first_r_ref_ptr;
  *lst_r_ref_ptr = last_r_ref_ptr;
}
/***********************************************************************

r_check_for_metal  -  Check whether the given atom is a metal ion

***********************************************************************/

int r_check_for_metal(char *atom_name,char *element)
{
  int imetal;
  int metal;

  /* Initialise variables */
  metal = FALSE;

  /* Loop through all the metals */
  for (imetal = 0; imetal < R_NMETALS && metal == FALSE; imetal++)
    {
      /* If name matches in first two characters, and third character
         is either a space or a number, then have a metal */
      if (!strncmp(atom_name,R_METAL_NAME[imetal],2) &&
          (atom_name[2] == ' ' ||
           (atom_name[2] >= '0' && atom_name[2] <= '9' &&
            (atom_name[3] == ' ' || atom_name[3] == '\'' ||
            (atom_name[3] >= '0' && atom_name[3] <= '9')))))
        {
          /* Set flag */
          metal = TRUE;

          /* Return the metal name matched */
          strncpy(element,R_METAL_NAME[imetal],2);
          element[2] = '\0';
        }
    }

  /* Return the answer */
  return(metal);
}
/***********************************************************************

is_hydrogen_atom  -  Check whether this is a hydrogen or a pseudo atom
                     from any of the many possible representations of
                     such in the PDB

***********************************************************************/

int is_hydrogen_atom(char *atname)
{
  char atom_name[5], element[3];

  int is_hydrogen, is_metal;

  /* Initialise */
  is_hydrogen = FALSE;
  strncpy(atom_name,atname,4);
  atom_name[4] = '\0';

  /* Check the atom name for any indicators to suggest this is either
     a hydrogen or a pseudo atom */
  if (atom_name[0] == 'H' || atom_name[1] == 'H' || atom_name[1] == 'Q' ||
      (atom_name[1] == 'D' && (atom_name[0] == ' ' || atom_name[0] == '1' ||
                               atom_name[0] == '2' || atom_name[0] == '3')))
    is_hydrogen = TRUE;

  /* Check for possible strange hydrogens */
  if (atom_name[0] >= '0' && atom_name[0] <= '9' && atom_name[2] == 'H')
    is_hydrogen = TRUE;

  /* If this looks likely to be a hydrogen, check that it isn't in
     fact a metal */
  if (is_hydrogen == TRUE)
    {
      /* Check whether this is a metal */
      is_metal = r_check_for_metal(atom_name,element);

      /* If a metal, then isn't a hydrogen */
      if (is_metal == TRUE
          && atom_name[3] != '\'' && atom_name[3] != '*')
        is_hydrogen = FALSE;
    }

  return is_hydrogen;
}
/***********************************************************************

get_element  -  Determine element from atom's name

***********************************************************************/

int get_element(char *atom_name,char *element)
{
  char metal_name[3];

  int ielem;
  int found, is_hydrogen, is_metal;

  static char *valid_elem[] = { " C", " N", " O", " P", " S" };

  static int nelem = 5;

  /* Initialise variables */
  is_metal = FALSE;
  strcpy(element,"  ");

  /* Check whether this is a hydrogen atom */
  is_hydrogen = is_hydrogen_atom(atom_name);

  /* If atom is a hydrogen, then have its element type */
  if (is_hydrogen == TRUE)
    strcpy(element," H");

  /* Otherwise, identify non-hydrogen element */
  else
    {
      /* First, whether this might be a metal */
      is_metal = r_check_for_metal(atom_name,metal_name);

      /* If it is, then store the element */
      if (is_metal == TRUE)
        strcpy(element,metal_name);

      /* Otherwise, take character in 2nd position to be representative
         of the element */
      else
        {
          /* Form the proposed element name */
          element[0] = ' ';
          element[1] = atom_name[1];
          element[2] = '\0';

          /* Check that this is a valid element name */
          found = FALSE;
          for (ielem = 0; ielem < nelem && found == FALSE; ielem++)
            {
              /* If have a match, then can stop looking */
              if (!strcmp(element,valid_elem[ielem]))
                found = TRUE;
            }

          /* If not found, then show warning message */
          if (found == FALSE)
            {
              /* Show warning */
              printf("*** Warning. Can't identify element for atom");
              printf(" [%s]\n",atom_name);
            }
        }
    }

  /* Return whether atom is a metal */
  return(is_metal);
}
/***********************************************************************

r_get_atom_details  -  Extract the atom details from the given line

***********************************************************************/

int r_get_atom_details(char *input_line,int *atom_number,char *atom_name,
                       char *chain,char *res_name,char *res_num,
                       double *xval,double *yval,double *zval,
                       double *bval,double *occup,char *element,
                       int *rec_type,char wanted_chain)
{
  char number_string[20];

  int len;
  int wanted;

  double bvalue, occupancy, x, y, z;

  /* Initialise flag */
  wanted = FALSE;

  /* Get the length of the input line */
  len = strlen(input_line);

  /* Get record type */
  if (!strncmp(input_line,"ATOM  ",6))
    *rec_type = R_ATOM_REC;
  else
    *rec_type = R_HETATM_REC;

  /* Get the residue details */
  strncpy(res_name,input_line+17,3);
  res_name[3] = '\0';
  strncpy(res_num,input_line+22,5);
  res_num[5] = '\0';
  *chain = input_line[21];

  /* If we only want a specific chain, check whether this is it */
  if (wanted_chain != '\0' && wanted_chain != *chain)
    return(wanted);

  /* Set flag as wanted */
  wanted = TRUE;

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

  /* If B-value greater than maximum, set to maximum */
  if (bvalue > R_MAX_BVALUE)
    bvalue = R_MAX_BVALUE;

  /* Retrieve the element, if present */
  if (len > 76)
    {
      strncpy(element,input_line+76,2);
      element[2] = '\0';

      /* Check whether this is a valid element */
      if ((element[0] >= '0' && element[0] <= '9') ||
          (element[1] >= '0' && element[1] <= '9'))
        element[0] = '\0';
    }
  else
    element[0] = '\0';

  /* Get the element if we don't already have it */
  if (element[0] == '\0')
    get_element(atom_name,element);

  /* Store the retrieved values */
  *xval = x;
  *yval = y;
  *zval = z;
  *bval = bvalue;
  *occup = occupancy;

  /* Return whether this atom is wanted */
  return(wanted);
}
/***********************************************************************

create_r_chain_record  -  Create and initialise a new r_chain record

***********************************************************************/

struct r_chain *create_r_chain_record(struct r_chain **fst_r_chain_ptr,
                                      struct r_chain **lst_r_chain_ptr,
                                      struct r_pdb *r_pdb_ptr,
                                      char chain_id)
{
  struct r_chain *r_chain_ptr, *first_r_chain_ptr, *last_r_chain_ptr;

  /* Initialise r_chain pointers */
  r_chain_ptr = NULL;
  first_r_chain_ptr = *fst_r_chain_ptr;
  last_r_chain_ptr = *lst_r_chain_ptr;

  /* Allocate memory for structure to hold r_chain info */
  r_chain_ptr = (struct r_chain *) malloc(sizeof(struct r_chain));
  if (r_chain_ptr == NULL)
    {
      printf("*** Can't allocate memory for struct r_chain\n");
      exit (1);
    }

  /* If this is the very first r_chain, then store pointer */
  if (first_r_chain_ptr == NULL)
    first_r_chain_ptr = r_chain_ptr;

  /* Add link from previous r_chain to the current one */
  if (last_r_chain_ptr != NULL)
    last_r_chain_ptr->next_r_chain_ptr = r_chain_ptr;
  last_r_chain_ptr = r_chain_ptr;

  /* Store the r_chain details */
  r_chain_ptr->chain_id = chain_id;
  r_chain_ptr->r_pdb_ptr = r_pdb_ptr;
  r_chain_ptr->next_r_chain_ptr = NULL;

  /* Initialise remaining fields */
  r_chain_ptr->type = R_CHAIN_UNKNOWN;
  r_chain_ptr->nresid = 0;
  r_chain_ptr->nhelices = 0;
  r_chain_ptr->nstrands = 0;
  r_chain_ptr->nstand = 0;
  r_chain_ptr->nbase = 0;
  r_chain_ptr->natom_recs = 0;
  r_chain_ptr->dominant_sst = R_HELIX;
  r_chain_ptr->symm_chain = FALSE;
  r_chain_ptr->nhetatm_recs = 0;
  r_chain_ptr->link_r_chain_ptr = NULL;
  r_chain_ptr->first_r_residue_ptr = NULL;

  /* Return current pointers */
  *fst_r_chain_ptr = first_r_chain_ptr;
  *lst_r_chain_ptr = last_r_chain_ptr;

  /* Return new r_chain pointer */
  return(r_chain_ptr);
}
/***********************************************************************

r_get_aa_code  -  Get single-letter code from amino acid name

***********************************************************************/

char r_get_aa_code(char *res_name)
{
  char aa_code;
  int iamino, icode;
  int namino;

  static char amino[] = "ACDEFGHIKLMNPQRSTVWYXBZMXCMYSX";
  static char *amino_name[] = {
    "ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU",
    "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR",
    "PCA ","ASX ","GLX ", "UNK", "XXX", "CSS", "MIL", "PTR", "SEP", "MSE"
  };

  /* Initialise variables */
  namino = strlen(amino) - 1;

  /* Get the corresponding three-letter code */
  iamino = namino;
  for (icode = 0; icode < namino; icode++)
    if (!strcmp(res_name,amino_name[icode]))
      iamino = icode;
  aa_code = amino[iamino];

  /* Return the single-letter code */
  return(aa_code);
}
/***********************************************************************

get_residue_type  -  Get the residue type from its name

***********************************************************************/

void get_residue_type(struct r_residue *r_residue_ptr)
{
  int iamino, ires, itype, namino, nbase, nonamino;
  int metal;

  static char amino[] = "ACDEFGHIKLMNPQRSTVWY";
  static char *amino_name[] = {
    "ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU",
    "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR"
  };

  static char non_amino[] = "CMMYSX";
  static char *non_amino_name[] = {
    "CSS", "MIL", "MSE", "PTR", "SEP", "UNK"
  };

  static char base[] = "ACGITUACGITUACGITU";
  static char *base_name[] = {
    "  A", "  C", "  G", "  I", "  T", "  U",
    " +A", " +C", " +G", " +I", " +T", " +U",
    " DA", " DC", " DG", " DI", " DT", " DU"
  };

  /* Initialise variables */
  iamino = -1;
  metal = FALSE;
  namino = strlen(amino);
  nonamino = strlen(non_amino);
  nbase = strlen(base);

  /* Loop through all the residues until find a match */
  for (ires = 0; ires < namino && iamino == -1; ires++)
    if (!strcmp(r_residue_ptr->res_name,amino_name[ires]))
      {
        /* Found residue name, so retrieve single letter code and
           set as standard amino acid */
        iamino = ires;
        r_residue_ptr->aa_code = amino[iamino];
        r_residue_ptr->type = R_TYPE_STANDARD;
      }

  /* If still not identified, repeat for the non-standard amino acids */
  for (ires = 0; ires < nonamino && iamino == -1; ires++)
    if (!strcmp(r_residue_ptr->res_name,non_amino_name[ires]))
      {
        /* Found residue name, so retrieve single letter code */
        iamino = ires;
        r_residue_ptr->aa_code = non_amino[iamino];
        r_residue_ptr->type = R_TYPE_NON_STANDARD;
      }

  /* Loop through all the residues until find a match */
  for (ires = 0; ires < nbase && iamino == -1; ires++)
    if (!strcmp(r_residue_ptr->res_name,base_name[ires]))
      {
        /* Found residue name, so retrieve single letter code */
        iamino = ires;
        r_residue_ptr->aa_code = base[iamino];
        r_residue_ptr->type = R_TYPE_NUCLEIC_ACID;
      }

  /* If still not found, check whether this might be a non-standard
     amino acid (ie one containing all three mainchain atoms N, CA
     and C) */
  if (iamino == -1 && (r_residue_ptr->nmain_chain == 3 ||
                       r_residue_ptr->have_calpha == TRUE))
    {
      r_residue_ptr->type = R_TYPE_NON_STANDARD;
      iamino = 0;
    }

  /* Check whether this might be a metal ion */
  if (iamino == -1 && r_residue_ptr->natoms == 1)
    {
      /* Check whether atom is a metal */
      metal = TRUE;
      for (itype = 1; itype < R_NATOM_TYPES && metal == TRUE; itype++)
        {
         /* If one of the standard atoms, then not a metal */
          if (!strcmp(r_residue_ptr->first_r_atom_ptr->element,
                      R_ELEMENT[itype]))
            metal = FALSE;
        }

      if (metal == TRUE)
        {
          /* Set flags to mark residue and chain as a single
             metal ion*/
          r_residue_ptr->type = R_TYPE_METAL_ION;
        }
    }
}
/***********************************************************************

create_r_residue_record  -  Create and initialise a new r_residue record

***********************************************************************/

struct r_residue *create_r_residue_record(struct r_residue **fst_r_residue_ptr,
                                          struct r_residue **lst_r_residue_ptr,
                                          struct r_pdb *r_pdb_ptr,char athet,
                                          struct r_chain *r_chain_ptr,
                                          char *res_name,char *res_num,
                                          struct r_rasdef *r_rasdef_ptr)
{
  int i;

  struct r_residue *r_residue_ptr, *first_r_residue_ptr, *last_r_residue_ptr;

  /* Initialise r_residue pointers */
  r_residue_ptr = NULL;
  first_r_residue_ptr = *fst_r_residue_ptr;
  last_r_residue_ptr = *lst_r_residue_ptr;

  /* Allocate memory for structure to hold r_residue info */
  r_residue_ptr = (struct r_residue *) malloc(sizeof(struct r_residue));
  if (r_residue_ptr == NULL)
    {
      printf("*** Can't allocate memory for struct r_residue\n");
      exit (1);
    }

  /* If this is the very first, save its pointer */
  if (first_r_residue_ptr == NULL)
    first_r_residue_ptr = r_residue_ptr;

  /* Add link from previous r_residue to the current one */
  if (last_r_residue_ptr != NULL)
    last_r_residue_ptr->next_r_residue_ptr = r_residue_ptr;
  last_r_residue_ptr = r_residue_ptr;

  /* Get the single-letter amino acid code */
  r_residue_ptr->aa_code = r_get_aa_code(res_name);

  /* Save known elements of the structure */
  r_residue_ptr->r_pdb_ptr = r_pdb_ptr;
  r_residue_ptr->athet = athet;
  strcpy(r_residue_ptr->res_name,res_name);
  strcpy(r_residue_ptr->res_num,res_num);
  r_residue_ptr->r_rasdef_ptr = r_rasdef_ptr;
  if (r_chain_ptr != NULL)
    r_residue_ptr->chain = r_chain_ptr->chain_id;
  else
    r_residue_ptr->chain = ' ';

  /* Initialise remaining elements */
  r_residue_ptr->seq_no = -1;
  r_residue_ptr->ref_resid = -1;
  r_residue_ptr->type = R_TYPE_UNKNOWN;
  r_residue_ptr->sec_struc = R_OTHER;
  r_residue_ptr->domain = 0;
  r_residue_ptr->natoms = 0;
  r_residue_ptr->nmain_chain = 0;
  r_residue_ptr->have_calpha = FALSE;
  r_residue_ptr->connected = FALSE;
  r_residue_ptr->accessibility = 0;
  r_residue_ptr->wanted = TRUE;
  r_residue_ptr->site = &null;
  r_residue_ptr->r_chain_ptr = r_chain_ptr;
  r_residue_ptr->r_ligand_ptr = NULL;
  r_residue_ptr->r_hetdesc_ptr = NULL;
  r_residue_ptr->first_r_atom_ptr = NULL;
  r_residue_ptr->first_r_rbond_ptr = NULL;

  /* Update chain's residue links */
  if (r_chain_ptr != NULL)
    {
      if (r_chain_ptr->first_r_residue_ptr == NULL)
        r_chain_ptr->first_r_residue_ptr = r_residue_ptr;
      if (strcmp(res_name,"HOH") && strcmp(res_name,"SOL") &&
          strcmp(res_name,"WAT"))
        r_chain_ptr->nresid++;
    }

  /* Initialise the limits */
  r_residue_ptr->r_limits_ptr.first = TRUE;
  for (i = 0; i < 3; i++)
    {
      r_residue_ptr->r_limits_ptr.valmin[i] = 0;
      r_residue_ptr->r_limits_ptr.valmax[i] = 0;
      r_residue_ptr->r_limits_ptr.c_of_g[i] = 0;
    }
  r_residue_ptr->r_limits_ptr.npoints = 0;

  /* Initialise pointer to the next r_residue */
  r_residue_ptr->next_r_residue_ptr = NULL;

  /* Return current pointers */
  *fst_r_residue_ptr = first_r_residue_ptr;
  *lst_r_residue_ptr = last_r_residue_ptr;

  /* Return new r_residue pointer */
  return(r_residue_ptr);
}
/***********************************************************************

r_check_for_duplicate  -  Check for an earlier copy of the current atom

***********************************************************************/

int r_check_for_duplicate(struct r_residue *r_residue_ptr,char *atom_name,
                        char *res_name,char *res_num,char chain)
{
  int iatom, natoms;
  int have_atom;

  struct r_atom *r_atom_ptr;

  /* Initialise variables */
  have_atom = FALSE;

  /* If have a current r_residue, and its name matches, check its atoms
     for an existing copy of the current atom */
  if (r_residue_ptr != NULL &&
      !strcmp(res_name,r_residue_ptr->res_name) &&
      !strcmp(res_num,r_residue_ptr->res_num) &&
      chain == r_residue_ptr->chain)
    {
      /* Get number of atoms */
      natoms = r_residue_ptr->natoms;

      /* Initialise pointer to first atom of this r_residue */
      r_atom_ptr = r_residue_ptr->first_r_atom_ptr;
      iatom = 0;
  
      /* Loop through all the residue's atoms to look for match */
      while (iatom < natoms && r_atom_ptr != NULL && have_atom == FALSE)
        {
          /* Check for duplication */
          if (!strcmp(atom_name,r_atom_ptr->atom_name))
            have_atom = TRUE;

          /* Get pointer to next atom in linked-list */
          r_atom_ptr = r_atom_ptr->next_r_atom_ptr;
          iatom++;
        }
    }

  /* Return flag */
  return(have_atom);
}
/***********************************************************************

create_r_atom_record  -  Create a new atom record

***********************************************************************/

struct r_atom *create_r_atom_record(struct r_atom **fst_r_atom_ptr,
                                    struct r_atom **lst_r_atom_ptr,
                                    struct r_residue *r_residue_ptr,
                                    int atom_number,char *atom_name,
                                    double x,double y,double z,
                                    double bvalue,double occupancy,
                                    char *element,int rec_type)
{
  int i;

  struct r_atom *r_atom_ptr;
  struct r_atom *first_r_atom_ptr, *last_r_atom_ptr;

  /* Initialise atom pointers */
  r_atom_ptr = NULL;
  last_r_atom_ptr = *lst_r_atom_ptr;
  first_r_atom_ptr = *fst_r_atom_ptr;

  /* Create atom entry */
  r_atom_ptr = (struct r_atom*)malloc(sizeof(struct r_atom));
  if (r_atom_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for");
      printf(" struct r_atom\n");
      printf("***        Program terminated with error.\n");
      exit (1);
    }

  /* Store the details for this entry */
  r_atom_ptr->atom_number = atom_number;
  strcpy(r_atom_ptr->atom_name,atom_name);
  r_atom_ptr->orig_coords[X] = r_atom_ptr->coords[X] = x;
  r_atom_ptr->orig_coords[Y] = r_atom_ptr->coords[Y] = y;
  r_atom_ptr->orig_coords[Z] = r_atom_ptr->coords[Z] = z;
  r_atom_ptr->bvalue = bvalue;
  r_atom_ptr->occupancy = occupancy;
  r_atom_ptr->rec_type = rec_type;
  r_atom_ptr->r_residue_ptr = r_residue_ptr;
  r_atom_ptr->type = 0;

  /* Get the element if we don't already have it */
  if (element[0] == '\0')
    get_element(atom_name,element);
  strcpy(r_atom_ptr->element,element);

  /* Increment count of this residue's atoms */
  r_residue_ptr->natoms++;

  /* If residue's first atom, store pointer to it */
  if (r_residue_ptr->first_r_atom_ptr == NULL)
    {
      /* Store atom pointer */
      r_residue_ptr->first_r_atom_ptr = r_atom_ptr;

      /* Store atom's coords as the max- and min-coords for the residue */
      for (i = 0; i < 3; i++)
        {
          r_residue_ptr->coords_min[i] = r_atom_ptr->coords[i];
          r_residue_ptr->coords_max[i] = r_atom_ptr->coords[i];
        }
    }

  /* Otherwise, update residue's min and max coordinates */
  else
    {
      for (i = 0; i < 3; i++)
        {
          if (r_atom_ptr->coords[i] < r_residue_ptr->coords_min[i])
            r_residue_ptr->coords_min[i] = r_atom_ptr->coords[i];
          if (r_atom_ptr->coords[i] > r_residue_ptr->coords_max[i])
            r_residue_ptr->coords_max[i] = r_atom_ptr->coords[i];
        }
    }

  /* Recompute residue's centroid */
  for (i = 0; i < 3; i++)
    r_residue_ptr->centroid[i]
      = (r_residue_ptr->coords_max[i] + r_residue_ptr->coords_min[i]) / 2;

  /* If atom is one of the 3 mainchain atoms (CA, N or C), increment
     residue's count of these */
  if (!strcmp(atom_name," CA ") || !strcmp(atom_name," N  ") ||
      !strcmp(atom_name," C  "))
    r_residue_ptr->nmain_chain++;

  /* Check for C-alpha */
  if (!strcmp(atom_name," CA "))
    r_residue_ptr->have_calpha = TRUE;

  /* Initialise pointers */
  r_atom_ptr->first_r_conect_ptr = NULL;
  r_atom_ptr->last_r_conect_ptr = NULL;

  /* Initialise pointer to next atom record */
  r_atom_ptr->next_r_atom_ptr = NULL;

  /* If this is the first entry stored, then update
     pointer to head of linked list */
  if (first_r_atom_ptr == NULL)
    first_r_atom_ptr = r_atom_ptr;

  /* Otherwise, make previous entry point to the current
     one */
  else
    last_r_atom_ptr->next_r_atom_ptr = r_atom_ptr;
                      
  /* Save the current atom's pointer */
  last_r_atom_ptr = r_atom_ptr;

  /* Return current pointers */
  *fst_r_atom_ptr = first_r_atom_ptr;
  *lst_r_atom_ptr = last_r_atom_ptr;
  return(r_atom_ptr);
}
/***********************************************************************

get_atom_data  -  Read in the data from the current ATOM/HETATM line

***********************************************************************/

int get_atom_data(char *input_line,struct r_pdb *r_pdb_ptr,char *lst_chn,
                  char *last_res_name,char *last_res_num,
                  struct r_rasdef *first_r_rasdef_ptr,char wanted_chain)
{
  char athet, atom_name[5], element[3], last_chain;
  char chain, res_name[4], res_num[6], occup_bvalue[13];

  int atom_number;
  int rec_type;
  int have_atom, wanted;

  double bvalue, occupancy, x, y, z;

  struct r_atom *r_atom_ptr = NULL;
  struct r_chain *r_chain_ptr = NULL;
  struct r_rasdef *r_rasdef_ptr = NULL;
  static struct r_rasdef *start_r_rasdef_ptr = NULL;
  struct r_residue *r_residue_ptr = NULL;

  /* Initialise variables */
  last_chain = *lst_chn;

  /* Initialise pointers */
  r_chain_ptr = r_pdb_ptr->last_r_chain_ptr;
  r_residue_ptr = r_pdb_ptr->last_r_residue_ptr;

  /* Retrieve the atom details from the input line */
  wanted
    = r_get_atom_details(input_line,&atom_number,atom_name,&chain,
                         res_name,res_num,&x,&y,&z,&bvalue,&occupancy,
                         element,&rec_type,wanted_chain);

  /* Check for duplicate records */
  if (wanted == TRUE && occupancy > 0.0 && occupancy < 1.0)
    {
      /* If have partial occupancy atom, accept only the first instance
         of it */
      have_atom = r_check_for_duplicate(r_residue_ptr,atom_name,res_name,
                                        res_num,chain);
      if (have_atom == TRUE)
        wanted = FALSE;
    }

  /* If atom wanted, then save it */
  if (wanted == TRUE)
    {
      /* Check whether atom belongs to a new chain */
      if (chain != last_chain)
        {
          /* Create a new chain record */
          r_chain_ptr
            = create_r_chain_record(&(r_pdb_ptr->first_r_chain_ptr),
                                    &(r_pdb_ptr->last_r_chain_ptr),
                                    r_pdb_ptr,chain);

          /* Increment count of chains */
          r_pdb_ptr->nchains++;
        }

      /* Check whether atom belongs to a new residue */
      if (chain != last_chain || strcmp(res_num,last_res_num) ||
          strcmp(res_name,last_res_name))
        {
          /* Save whether ATOM or HETATM record */
          athet = input_line[0];

          /* Determine residue type using the definitions read in
             from the rasmol.dfn file, if we have them */
          r_rasdef_ptr
            = r_assign_to_object(first_r_rasdef_ptr,atom_name,res_name,
                                 res_num,chain);

          /* If chain type not yet defined, store its type */
          if (r_rasdef_ptr != NULL && r_chain_ptr->type == R_UNKNOWN)
            r_chain_ptr->type = r_rasdef_ptr->type;

          /* Create a new residue record */
          r_residue_ptr
            = create_r_residue_record(&(r_pdb_ptr->first_r_residue_ptr),
                                      &(r_pdb_ptr->last_r_residue_ptr),
                                      r_pdb_ptr,athet,r_chain_ptr,res_name,
                                      res_num,r_rasdef_ptr);

          /* If this residue is a continuation of a ligand, then get
             ligand's first rasdef record */
          if (r_rasdef_ptr != NULL && r_rasdef_ptr->type == R_LIGAND &&
              start_r_rasdef_ptr != NULL &&
              r_rasdef_ptr->number == start_r_rasdef_ptr->number)
            r_rasdef_ptr = start_r_rasdef_ptr;
          else
            start_r_rasdef_ptr = NULL;

          /* If we have a rasdef definition, update its start residue
             or number of residues */
          if (r_rasdef_ptr != NULL)
            {
              /* If we don't yet have a first residue, save this as the
                 first */
              if (r_rasdef_ptr->first_r_residue_ptr == NULL)
                r_rasdef_ptr->first_r_residue_ptr = r_residue_ptr;

              /* Increment count of residues belonging to this object */
              r_rasdef_ptr->nresid++;

              /* Save this as a possible start */
              if (r_rasdef_ptr->type == R_LIGAND)
                start_r_rasdef_ptr = r_rasdef_ptr;
            }

          /* Save these residue details */
          last_chain = chain;
          strcpy(last_res_name,res_name);
          strcpy(last_res_num,res_num);

          /* Increment residue count */
          r_pdb_ptr->nresid++;
        }

      /* Create a new atom record */
      r_atom_ptr
        = create_r_atom_record(&(r_pdb_ptr->first_r_atom_ptr),
                               &(r_pdb_ptr->last_r_atom_ptr),
                               r_residue_ptr,atom_number,atom_name,x,
                               y,z,bvalue,occupancy,element,rec_type);

      /* Increment count of atoms */
      r_pdb_ptr->natoms++;
    }

  /* Save the last chain id */
  *lst_chn = last_chain;

   /* Return last atom number */
   return(atom_number);
}
/***********************************************************************

r_find_atom  -  Find the given atom record

***********************************************************************/

struct r_atom *r_find_atom(struct r_atom *first_r_atom_ptr,
                           int atom_number)
{
  struct r_atom *r_atom_ptr, *match_r_atom_ptr;

  /* Initialise variables */
  match_r_atom_ptr = NULL;

  /* Get pointer to the first atom */
  r_atom_ptr = first_r_atom_ptr;

  /* Loop over the atoms to write out */
  while (r_atom_ptr != NULL && match_r_atom_ptr == NULL)
    {
      /* Check if this is the required record */
      if (r_atom_ptr->atom_number == atom_number)
        {
          /* Save the atom pointer */
          match_r_atom_ptr = r_atom_ptr;
        }

      /* Get pointer to next atom record */
      r_atom_ptr = r_atom_ptr->next_r_atom_ptr;
    }
      
  /* Return match */
  return(match_r_atom_ptr);
}
/***********************************************************************

create_r_rbond_record  -  Create a connection from the first residue to
                          the second

***********************************************************************/

void create_r_rbond_record(struct r_residue *r_residue1_ptr,
                           struct r_residue *r_residue2_ptr)
{
  struct r_rbond *r_rbond_ptr;

  /* Create r_rbond entry */
  r_rbond_ptr = (struct r_rbond*)malloc(sizeof(struct r_rbond));
  if (r_rbond_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for");
      printf(" struct r_rbond\n");
      printf("***        Program terminated with error.\n");
      exit (1);
    }

  /* Store the details for this entry */
  r_rbond_ptr->to_r_residue_ptr = r_residue2_ptr;
  r_rbond_ptr->next_r_rbond_ptr = NULL;

  /* If this is the first entry stored, then update pointer to
     head of linked list */
  if (r_residue1_ptr->first_r_rbond_ptr == NULL)
    r_residue1_ptr->first_r_rbond_ptr = r_rbond_ptr;

  /* Otherwise, make current entry point to the first one and make it
     the new first */
  else
    {
      r_rbond_ptr->next_r_rbond_ptr = r_residue1_ptr->first_r_rbond_ptr;
      r_residue1_ptr->first_r_rbond_ptr = r_rbond_ptr;
    }                 
}
/***********************************************************************

create_r_conect_record  -  Create a connection from the first atom to
                           the second

***********************************************************************/

void create_r_conect_record(struct r_atom *r_atom1_ptr,
                            struct r_atom *r_atom2_ptr)
{
  int have_conect;

  struct r_conect *r_conect_ptr;

  /* Initialise flag */
  have_conect = FALSE;

  /* Get pointer to the first of this atom's CONECT record */
  r_conect_ptr = r_atom1_ptr->first_r_conect_ptr;

  /* Loop over the CONECT records */
  while (r_conect_ptr != NULL && have_conect == FALSE)
    {
      /* If the CONECT links to the second atom, then we already have
         this bond */
      if (r_conect_ptr->to_r_atom_ptr == r_atom2_ptr)
        have_conect = TRUE;

      /* Get pointer to next CONECT record */
      r_conect_ptr = r_conect_ptr->next_r_conect_ptr;
    }

  /* If we already have a CONECT record between these two atoms,
     then return */
  if (have_conect == TRUE)
    return;

  /* Create conect entry */
  r_conect_ptr = (struct r_conect*) malloc(sizeof(struct r_conect));
  if (r_conect_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for");
      printf(" struct conect\n");
      printf("***        Program terminated with error.\n");
      exit (1);
    }

  /* Store the details for this entry */
  r_conect_ptr->to_r_atom_ptr = r_atom2_ptr;
  r_conect_ptr->next_r_conect_ptr = NULL;

  /* If this is the first entry stored, then update pointer to
     head of linked list */
  if (r_atom1_ptr->first_r_conect_ptr == NULL)
    r_atom1_ptr->first_r_conect_ptr = r_conect_ptr;

  /* Otherwise, make previous entry point to the current one */
  else
    r_atom1_ptr->last_r_conect_ptr->next_r_conect_ptr = r_conect_ptr;

  /* If atoms belong to different residues, create a connect record
     between the two residues */
  if (r_atom1_ptr->r_residue_ptr != r_atom2_ptr->r_residue_ptr)
    create_r_rbond_record(r_atom1_ptr->r_residue_ptr,
                          r_atom2_ptr->r_residue_ptr);
                      
  /* Save the current conect's pointer */
  r_atom1_ptr->last_r_conect_ptr = r_conect_ptr;
}
/***********************************************************************

get_conect_data  -  Extract the atom connectivity in the given CONECT
                    record

***********************************************************************/

void get_conect_data(char *input_line,struct r_atom *first_r_atom_ptr,
                     int natoms)
{
  char atmnum[10], element[3];

  int atom_number, iatom, iconect, ipos, len;
  int done, ok;

  double dist_sqrd, x1, x2, y1, y2, z1, z2;

  struct r_atom *r_atom_ptr, *to_r_atom_ptr;

  /* Initialise variables */
  element[0] = '\0';
  iatom = 0;
  done = FALSE;
  ipos = 6;

  /* Initialise pointers */
  r_atom_ptr = NULL;
  to_r_atom_ptr = NULL;

  /* Get line length */
  len = strlen(input_line);

  /* Loop through all the possible atom-number positions */
  for (iconect = 0; iconect < R_MAX_CONECT && done == FALSE; iconect++)
    {
      /* If this position is blank, then have reached end of atom-numbers
         on this line */
      if (ipos + 5 > len || !strncmp(input_line+ipos,"     ",5))
        done = TRUE;

      /* Otherwise, check for this atom number */
      else
        {
          /* Get the atom-number at this position */
          strncpy(atmnum,input_line+ipos,5);
          atmnum[5] = '\0';
          atom_number = atoi(atmnum);

          /* If the first atom number, store as first pointer */
          if (iatom == 0)
            {
              /* Get the corresponding atom pointer */
              r_atom_ptr = r_find_atom(first_r_atom_ptr,atom_number);

              /* If this is null, then can't make use of the
                 current CONECT line */
              if (r_atom_ptr == NULL)
                done = TRUE;

              /* Otherwise, retrieve atom's coordinates */
              else
                {
                  x1 = r_atom_ptr->coords[0];
                  y1 = r_atom_ptr->coords[1];
                  z1 = r_atom_ptr->coords[2];
                }
            }

          /* Otherwise, get which atom first one is bonded to and
             check out this bond */
          else
            {
              /* If first position is null, then no point looking
                 at this, or any further, positions */
              if (r_atom_ptr == NULL)
                done = TRUE;

              /* Otherwise, get the bonded atom */
              else
                {
                  /* Get the corresponding atom pointer */
                  to_r_atom_ptr
                    = r_find_atom(first_r_atom_ptr,atom_number);

                  /* If this is not null, then check out the
                     connection */
                  if (to_r_atom_ptr != NULL)
                    {
                      /* Retrieve atom's coordinates */
                      x2 = to_r_atom_ptr->coords[0];
                      y2 = to_r_atom_ptr->coords[1];
                      z2 = to_r_atom_ptr->coords[2];

                      /* Calculate squared distance between
                         these two atoms */
                      dist_sqrd
                        = (x2 - x1) * (x2 - x1)
                        + (y2 - y1) * (y2 - y1)
                        + (z2 - z1) * (z2 - z1);

                      /* Initialise flag */
                      ok = TRUE;

                      /* If distance greater than standard cut-off,
                         check if either atom is a metal */
                      if (dist_sqrd > R_MAX_CONECT_LEN_SQRD)
                        {
                          /* If either atom is a metal, allow
                             extra distance in the bond */
                          if (r_check_for_metal(r_atom_ptr->atom_name,
                                                element) == TRUE ||
                              r_check_for_metal(to_r_atom_ptr->atom_name,
                                                element) == TRUE)
                            {
                              /* If bond too long even for metals,
                                 then discard */
                              if (dist_sqrd > R_MAX_METAL_LEN_SQRD)
                                ok = FALSE;
                            }

                          /* If neither is metal, then bond too long */
                          else
                            ok = FALSE;
                        }

                      /* If distance is a reasonable one, store
                         a connection between these two atoms */
                      if (ok == TRUE)
                        {
                          /* Create a connection from the first
                             atom to the second */
                          create_r_conect_record(r_atom_ptr,
                                                 to_r_atom_ptr);

                          /* Create a connection the other way
                             round */
                          create_r_conect_record(to_r_atom_ptr,
                                                 r_atom_ptr);
                        }

                      /* If distance too great, show warning message */
                      else
                        {
                          printf("* Warning. Long bond from CONECT: "
                                 "%5d. %s %s %s %c -> %5d. %s %s %s %c"
                                 "%8.3f\n",
                                 r_atom_ptr->atom_number,
                                 r_atom_ptr->atom_name,
                                 r_atom_ptr->r_residue_ptr->res_name,
                                 r_atom_ptr->r_residue_ptr->res_num,
                                 r_atom_ptr->r_residue_ptr->chain,
                                 to_r_atom_ptr->atom_number,
                                 to_r_atom_ptr->atom_name,
                                 to_r_atom_ptr->r_residue_ptr->res_name,
                                 to_r_atom_ptr->r_residue_ptr->res_num,
                                 to_r_atom_ptr->r_residue_ptr->chain,
                                 sqrt(dist_sqrd));
                        }
                    }
                }
            }

          /* Increment count of atom position */
          iatom++;
        }

      /* Shift position along input line */
      ipos = ipos + 5;
    }
}
/***********************************************************************

r_read_pdb_file  -  Read in the data from the PDB file

***********************************************************************/

#ifdef TEST
struct r_pdb *r_read_pdb_file(FILE *file_ptr,int gzipped,char *pdb_code,
#else
struct r_pdb *r_read_pdb_file(gzFile file_ptr,int gzipped,char *pdb_code,
#endif
                              char wanted_chain,int data_rqrd,
                              struct r_rasdef *first_r_rasdef_ptr)
{
  char input_line[R_LINELEN + 1];
  char last_chain, last_res_num[6], last_res_name[4];
  char last_reference[5], number_string[20], reference[5];
  char string1[R_STRING_LEN], string2[R_STRING_LEN];
  char string3[R_STRING_LEN], string4[R_STRING_LEN];
  char string5[R_STRING_LEN], string6[R_STRING_LEN];

  char *string_ptr;

  int ikey, keylen[NKEY_WORDS], keyword, last_keyword, len, nlines;
  int molecule, nprotein, nref[2];
  int done, eof, keep_reading, flag[NKEY_WORDS];
  int protein_no[R_MAX_MOLECS];
  int first_compnd, have_error, have_header, have_title, old_style;
  int synthetic;

  struct r_info *first_r_info_ptr, *last_r_info_ptr, *r_info_ptr;
  struct r_pdb *r_pdb_ptr;
  struct r_ref *first_r_ref_ptr, *last_r_ref_ptr, *r_ref_ptr;

  /* Initialise variables */
  keep_reading = FALSE;
  first_compnd = TRUE;
  have_error = FALSE;
  have_title = FALSE;
  have_header = FALSE;
  string1[0] = '\0';
  string2[0] = '\0';
  string3[0] = '\0';
  string4[0] = '\0';
  string5[0] = '\0';
  string6[0] = '\0';
  synthetic = FALSE;
  for (molecule = 0; molecule < R_MAX_MOLECS; molecule++)
    protein_no[molecule] = 0;
  molecule = 0;
  nlines = 0;
  nprotein = 0;
  old_style = FALSE;
  last_reference[0] = '\0';
  reference[0] = '\0';
  last_keyword = NONE;
  keyword = NONE;

  /* Initialise variables */
  last_chain = '\0';
  last_res_name[0] = '\0';
  last_res_num[0] = '\0';

  /* Initialise pointers */
  r_info_ptr = first_r_info_ptr = last_r_info_ptr = NULL;
  r_ref_ptr = first_r_ref_ptr = last_r_ref_ptr = NULL;
  nref[0] = 0;
  nref[1] = 0;

  /* Create a new PDB entry */
  r_pdb_ptr = create_r_pdb_entry(pdb_code);

  /* Precalculate the lengths of the key words and initialise flags */
  for (ikey = 0; ikey < NKEY_WORDS; ikey++)
    {
      keylen[ikey] = strlen(key[ikey]);
      flag[ikey] = NONE;
    }

  /* Initialise flag */
  eof = done = FALSE;
  keep_reading = TRUE;

  /* Loop through the PDB file, picking up all the atoms that are
     required */
  while (eof == FALSE && done == FALSE)
    {
      /* Read in the next record */
#ifdef TEST
      if (fgets(input_line,R_LINELEN,file_ptr) == NULL)
        eof = TRUE;
#else
      if (gzgets(file_ptr,input_line,R_LINELEN) == NULL)
        eof = TRUE;
#endif

      /* If record read in OK, then process it */
      if (eof == FALSE)
        {
          /* Truncate the string to strip off the newline */
          len = string_chop(input_line);
          nlines++;

          /* Check determine which keyword we have */
          keyword = NONE;
          for (ikey = 0; ikey < NKEY_WORDS && keyword == NONE; ikey++)
            {
              /* If have match, set flag */
              if (!strncmp(input_line,key[ikey],keylen[ikey]))
                keyword = ikey;
            }

          /* If old-style COMPND, set keyword to be TITLE */
          if (keyword == COMPOUND && old_style == TRUE)
            keyword = TITLE;

          /* If key word has changed, Store the previous data if
             any */
          if ((strstr(input_line,"MOL_ID: ") != NULL ||
               keyword != last_keyword || strcmp(reference,last_reference))
              && last_keyword != NONE && have_header == TRUE)
            {
              /* Store the data */
              if (string1[0] != '\0' || string2[0] != '\0' ||
                  string3[0] != '\0' || string4[0] != '\0' ||
                  string5[0] != '\0' || string6[0] != '\0' ||
                  synthetic == TRUE)
                {
                  /* If this is a synthetic construct, and we don't have
                     any source organism, set as synthetic */
                  if (synthetic == TRUE && string1[0] == '\0')
                    {
                      strcpy(string1,"SYNTHETIC CONSTRUCT");
                      strcpy(string2,"SYNTHETIC");
                      strcpy(string3,"32630");
                    }

                  /* Store the data */
                  r_store_data(last_keyword,string1,string2,
                               string3,string4,string5,string6,protein_no,
                               molecule,last_reference,&first_r_info_ptr,
                               &last_r_info_ptr,&first_r_ref_ptr,
                               &last_r_ref_ptr);

                  /* Reset the synthetic flag */
                  synthetic = FALSE;

                  /* If this is the TITLE line, set flag */
                  if (last_keyword == TITLE)
                    have_title = TRUE;
                }

              /* Unset header flag */
              have_header = FALSE;
            }

          /* If this is the first COMPND record, determine whether it
             is old- or new-style */
          if (keyword == COMPOUND && first_compnd == TRUE)
            {
              /* If line doesn't contain MOL_ID then take it to be an
                 old-style entry and take this to be the title */
              if (strstr(input_line,"MOL_ID: ") == NULL)
                {
                  /* Set all the flags */
                  keyword = TITLE;
                  old_style = TRUE;
                }

              /* Unset first-compound flag */
              first_compnd = FALSE;
            }

          /* Save the current key word and reference */
          last_keyword = keyword;
          strcpy(last_reference,reference);

          /* Process line according to current keyword */
          switch (keyword)
            {
              /* a. Header records */
            case TITLE :
            case COMPOUND :
            case SOURCE :
            case KEY_WORDS :
            case AUTHOR :
            case JOURNAL :
            case REFERENCE :

              /* If we're not interested in the header info, break here */
              if (data_rqrd == COORDS_ONLY)
                break;

              /* Process current line */
              molecule
                = r_process_line(input_line+6,string1,string2,string3,
                                 string4,string5,string6,&synthetic,keyword,
                                 molecule,protein_no,&nprotein,
                                 reference,flag,first_r_rasdef_ptr);

              /* If this is a COMPND record, unset first flag */
              if (keyword == COMPOUND)
                first_compnd = FALSE;

              /* Set flag that have a header record */
              have_header = TRUE;

              break;

              /* b. Start of new model */
            case MODEL :

              /* Increment count of models */
              r_pdb_ptr->nmodels++;

              /* c. End of model */
            case ENDMDL :

              /* If already have a model, don't need to read in any more atom
                 records */
              if (r_pdb_ptr->nmodels > 1)
                keep_reading = FALSE;

              break;

              /* d. Atom data */
            case ATOM :
            case HETATM :

              /* If only interested in the header records, can stop
                 reading here */
              if (data_rqrd == HEADERS_ONLY)
                {
                  keep_reading = FALSE;
                  done = TRUE;
                }

              /* Pick up the atom data */
              if (keep_reading == TRUE)
                get_atom_data(input_line,r_pdb_ptr,&last_chain,last_res_name,
                              last_res_num,first_r_rasdef_ptr,wanted_chain);

              break;

              /* e. TER records */
            case TER :

              /* Re-initialise last chain id */
              last_chain = '\0';

              break;

              /* f. Connect records */
            case CONECT :

              /* Extract the atom connectivities from this line */
              get_conect_data(input_line,r_pdb_ptr->first_r_atom_ptr,
                              r_pdb_ptr->natoms);

              break;

            default :
              break;
            }

          /* Save the current key word and reference */
          last_keyword = keyword;
        }
    }

  /* Save the pointers */
  r_pdb_ptr->first_r_info_ptr = first_r_info_ptr;
  r_pdb_ptr->last_r_info_ptr = last_r_info_ptr;
  r_pdb_ptr->first_r_ref_ptr = first_r_ref_ptr;
  r_pdb_ptr->last_r_ref_ptr = last_r_ref_ptr;

  /* Return the PDB record */
  return(r_pdb_ptr);
}
/***********************************************************************

r_open_pdb_file  -  Open the given PDB file

***********************************************************************/

#ifdef TEST
FILE *r_open_pdb_file(char *file_name,int *gz)
#else
gzFile r_open_pdb_file(char *file_name,int *gz)
#endif
{
  int len;
  int eof, gzipped;

#ifdef TEST
  FILE *file_ptr;
#else
  gzFile file_ptr;
#endif

  /* Initialise variables */
  file_ptr = NULL;
  if (file_name == &null)
    return(file_ptr);

  /* Determine whether the PDB file is a gzipped one */
  eof = gzipped = FALSE;
  len = strlen(file_name);
  if (!strncmp(file_name + len - 3,".gz",3))
    gzipped = TRUE;
#ifdef TEST
  gzipped = FALSE;
#endif

  /* Open PDB file */
#ifdef TEST
  file_ptr = fopen(file_name,"r");
#else
  file_ptr = gzopen(file_name,"r");
#endif

  /* Return whether file is gzipped */
  *gz = gzipped;

  /* Return the file pointer */
  return(file_ptr);
}
/***********************************************************************

identify_residues  -  Identify residue types

***********************************************************************/

void identify_residues(struct r_pdb *r_pdb_ptr)
{
  int i;

  struct r_residue *r_residue_ptr;

  /* Get pointer to the first residue */
  r_residue_ptr = r_pdb_ptr->first_r_residue_ptr;

  /* Loop through all the residues to classify */
  while (r_residue_ptr != NULL)
    {
      /* Get the residue type and amino acid code from its name */
      get_residue_type(r_residue_ptr);

      /* Get the next residue */
      r_residue_ptr = r_residue_ptr->next_r_residue_ptr;
    }
}
/***********************************************************************

create_r_molec_record  -  Create and initialise a new r_molec record

***********************************************************************/

struct r_molec *create_r_molec_record(struct r_molec **fst_r_molec_ptr,
                                      struct r_molec **lst_r_molec_ptr,
                                      struct r_pdb *r_pdb_ptr,
                                      struct r_residue *r_residue_ptr,
                                      int type)
{
  int i;

  struct r_molec *r_molec_ptr, *first_r_molec_ptr, *last_r_molec_ptr;

  /* Initialise r_molec pointers */
  r_molec_ptr = NULL;
  first_r_molec_ptr = *fst_r_molec_ptr;
  last_r_molec_ptr = *lst_r_molec_ptr;

  /* Allocate memory for structure to hold r_molec info */
  r_molec_ptr = (struct r_molec *) malloc(sizeof(struct r_molec));
  if (r_molec_ptr == NULL)
    {
      printf("*** Can't allocate memory for struct r_molec\n");
      exit (1);
    }

  /* If this is the very first r_molec, then store pointer */
  if (first_r_molec_ptr == NULL)
    first_r_molec_ptr = r_molec_ptr;

  /* Add link from previous r_molec to the current one */
  if (last_r_molec_ptr != NULL)
    last_r_molec_ptr->next_r_molec_ptr = r_molec_ptr;
  last_r_molec_ptr = r_molec_ptr;

  /* Store the r_molec details */
  r_molec_ptr->type = type;
  r_molec_ptr->r_pdb_ptr = r_pdb_ptr;
  r_molec_ptr->next_r_molec_ptr = NULL;
  r_molec_ptr->first_r_residue_ptr = r_residue_ptr;

  /* Initialise remaining fields */
  r_molec_ptr->nresid = 0;

  /* Initialise the limits */
  r_molec_ptr->r_limits_ptr.first = TRUE;
  for (i = 0; i < 3; i++)
    {
      r_molec_ptr->r_limits_ptr.valmin[i] = 0;
      r_molec_ptr->r_limits_ptr.valmax[i] = 0;
      r_molec_ptr->r_limits_ptr.c_of_g[i] = 0;
    }
  r_molec_ptr->r_limits_ptr.npoints = 0;

  /* Return current pointers */
  *fst_r_molec_ptr = first_r_molec_ptr;
  *lst_r_molec_ptr = last_r_molec_ptr;

  /* Return new r_molec pointer */
  return(r_molec_ptr);
}
/***********************************************************************

r_assign_to_molecules  -  Assign residues to their individual molecules

***********************************************************************/

int r_assign_to_molecules(struct r_pdb *r_pdb_ptr)
{
  int i, last_number, r_molec_type, ncalpha, nmolecs, type;
  int new_molec;

  struct r_molec *r_molec_ptr, *first_r_molec_ptr, *last_r_molec_ptr;
  struct r_rasdef *r_rasdef_ptr, *last_r_rasdef_ptr;
  struct r_residue *r_residue_ptr;

  /* Initialise variables */
  nmolecs = 0;
  last_r_rasdef_ptr = NULL;

  /* Initialise pointers */
  r_molec_ptr = first_r_molec_ptr = last_r_molec_ptr = NULL;

  /* Get pointer to the first residue read in from the PDB file */
  r_residue_ptr = r_pdb_ptr->first_r_residue_ptr;

  /* Loop through all the residues to pick out the data we need */
  while (r_residue_ptr != NULL)
    {
      /* Get this residue's definition */
      r_rasdef_ptr = r_residue_ptr->r_rasdef_ptr;

      /* Get the molecule type */
      type = R_UNKNOWN;
      new_molec = FALSE;

      /* See if this is a new rasdef definition */
      if (r_rasdef_ptr != NULL)
        {
          /* If this is a new definition, then want to create a new
             molecule */
          if (r_rasdef_ptr->type == R_LIGAND)
            {
              /* Check if this residue is a continuation of the same
                 ligand */
              if (last_r_rasdef_ptr == NULL ||
                  last_r_rasdef_ptr->type != R_LIGAND ||
                  last_r_rasdef_ptr->number != r_rasdef_ptr->number)
                new_molec = TRUE;
            }
          else if (r_rasdef_ptr != last_r_rasdef_ptr ||
                   r_rasdef_ptr->type == R_METAL)
            new_molec = TRUE;

          /* Save the type and the definition */
          type = r_rasdef_ptr->type;
          last_r_rasdef_ptr = r_rasdef_ptr;
        }

      else if (!strcmp(r_residue_ptr->res_name,"HOH"))
        {
          type = R_WATER;
          new_molec = TRUE;
        }

      else if (r_residue_ptr->natoms == 1)
        {
          type = R_METAL;
          new_molec = TRUE;
        }

      /* If this is a new molecule, create new record */
      if (new_molec == TRUE)
        {
          /* Create a molecule record */
          r_molec_ptr
            = create_r_molec_record(&first_r_molec_ptr,
                                    &last_r_molec_ptr,r_pdb_ptr,
                                    r_residue_ptr,type);

          /* Increment count of residues belonging to this molecule */
          r_molec_ptr->nresid++;

          /* Increment residue count */
          nmolecs++;
        }

      /* If this is the same molecule as before, increment residue count */
      else if (r_rasdef_ptr != NULL)
        r_molec_ptr->nresid++;

      else
        printf("*** Residue %s %s %c unassigned!\n",
               r_residue_ptr->res_name,
               r_residue_ptr->res_num,
               r_residue_ptr->chain);

      /* Get the next r_residue */
      r_residue_ptr = r_residue_ptr->next_r_residue_ptr;
    }

  /* Show number of molecule records created */
  printf("Number of molecules created:      %8d\n",nmolecs);

  /* Save the first molecule pointer */
  r_pdb_ptr->first_r_molec_ptr = first_r_molec_ptr;

  /* Return number of molecules created */
  return(nmolecs);
}
/***********************************************************************

r_get_max_min_coords  -  Get the maximum and minimum x-, y-, and z-coords
                       of all atoms

***********************************************************************/

void r_get_max_min_coords(struct r_limits *limit,double x,double y,double z)
{
  /* If this is the first atom encountered, store its
     coords as the present max and min values */
  if (limit->first == TRUE)
    {
      limit->c_of_g[X] = x;
      limit->valmin[X] = x;
      limit->valmax[X] = x;
      limit->c_of_g[Y] = y;
      limit->valmin[Y] = y;
      limit->valmax[Y] = y;
      limit->c_of_g[Z] = z;
      limit->valmin[Z] = z;
      limit->valmax[Z] = z;
      limit->first = FALSE;
      limit->npoints = 1;
    }

  /* Otherwise, update the limits if coordinates are outside them */
  else
    {
      limit->c_of_g[X] = limit->c_of_g[X] + x;
      limit->c_of_g[Y] = limit->c_of_g[Y] + y;
      limit->c_of_g[Z] = limit->c_of_g[Z] + z;
      if (x < limit->valmin[0])
        limit->valmin[0] = x;
      if (y < limit->valmin[1])
        limit->valmin[1] = y;
      if (z < limit->valmin[2])
        limit->valmin[2] = z;
      if (x > limit->valmax[0])
        limit->valmax[0] = x;
      if (y > limit->valmax[1])
        limit->valmax[1] = y;
      if (z > limit->valmax[2])
        limit->valmax[2] = z;
      limit->npoints++;
    }
}
/***********************************************************************

r_calc_limits  -  Calculate the x-, y- and z- atom limits for each
                  molecule and each residue

***********************************************************************/

void r_calc_limits(struct r_pdb *r_pdb_ptr)
{
  int iatom, icoord, iresid, natoms, nmolec, nresid, type;
  int wanted;

  double coords[3], margin;

  struct r_atom *r_atom_ptr;
  struct r_molec *r_molec_ptr;
  struct r_residue *r_residue_ptr;

  /* Initialise variables */
  nmolec = 0;

  /* Get pointer to the first molecule */
  r_molec_ptr = r_pdb_ptr->first_r_molec_ptr;

  /* Loop through all the molecules to calc coord limits */
  while (r_molec_ptr != NULL)
    {
      /* Increment molecule count */
      nmolec++;

      /* Get pointer to this molecule's first residue */
      r_residue_ptr = r_molec_ptr->first_r_residue_ptr;
      nresid = r_molec_ptr->nresid;
      iresid = 0;

      /* Loop through all the r_residues to store the transformed coords */
      while (r_residue_ptr != NULL && iresid < nresid)
        {
          /* Initialise pointer to first r_atom of this r_residue */
          r_atom_ptr = r_residue_ptr->first_r_atom_ptr;
          iatom = 0;
          natoms = r_residue_ptr->natoms;

          /* Loop through all the r_residue's r_atoms to count */
          while (iatom < natoms && r_atom_ptr != NULL)
            {
              /* Update the residue maximum and minimum coordinates */
              r_get_max_min_coords(&(r_residue_ptr->r_limits_ptr),
                                   r_atom_ptr->coords[X],
                                   r_atom_ptr->coords[Y],
                                   r_atom_ptr->coords[Z]);

              /* Repeat for maximum and minimum molecule coordinates */
              r_get_max_min_coords(&(r_molec_ptr->r_limits_ptr),
                                   r_atom_ptr->coords[X],
                                   r_atom_ptr->coords[Y],
                                   r_atom_ptr->coords[Z]);

              /* Get pointer to next r_atom in linked-list */
              r_atom_ptr = r_atom_ptr->next_r_atom_ptr;
              iatom++;
            }

          /* Calculate residue's C of G */
          if (r_residue_ptr->r_limits_ptr.npoints > 0)
            {
              for (icoord = 0; icoord < 3; icoord++)
                r_residue_ptr->r_limits_ptr.c_of_g[icoord]
                  = r_residue_ptr->r_limits_ptr.c_of_g[icoord]
                  / r_residue_ptr->r_limits_ptr.npoints;
            }

          /* Get the next r_residue */
          r_residue_ptr = r_residue_ptr->next_r_residue_ptr;
          iresid++;
        }


      /* Calculate molecule's C of G */
      if (r_molec_ptr->r_limits_ptr.npoints > 0)
        {
          for (icoord = 0; icoord < 3; icoord++)
            r_molec_ptr->r_limits_ptr.c_of_g[icoord]
              = r_molec_ptr->r_limits_ptr.c_of_g[icoord]
              / r_molec_ptr->r_limits_ptr.npoints;
        }
  
      /* Get the next molecule */
      r_molec_ptr = r_molec_ptr->next_r_molec_ptr;
    }

  /* Use maximum radius to compute extra border to be added round
     r_atom r_limits */
  margin = R_MAX_ATOM_RADIUS;
}
/***********************************************************************

r_assign_to_ligands  -  Assign residues to their individual ligands

***********************************************************************/

int r_assign_to_ligands(struct r_pdb *r_pdb_ptr)
{
  char chain, res_name[4], res_num[6];

  int nligands, nmetals;
  int found;

  struct r_ligand *r_ligand_ptr;
  struct r_residue *r_residue_ptr;

  /* Initialise variables */
  nligands = 0;
  nmetals = 0;

  /* Get pointer to the first ligand in the PDB entry */
  r_ligand_ptr = r_pdb_ptr->first_r_ligand_ptr;

  /* Loop through all the ligands to calc coord limits */
  while (r_ligand_ptr != NULL)
    {
      /* If ligand, then identify the start and end residues */
      if (r_ligand_ptr->metal == FALSE)
        {
          /* Increment ligand count */
          nligands++;

          /* Extract the residue details of the first residue */
          chain = r_ligand_ptr->first_residue[8];
          strncpy(res_name,r_ligand_ptr->first_residue,3);
          res_name[3] ='\0';
          strncpy(res_num,r_ligand_ptr->first_residue + 3,5);
          res_num[5] ='\0';

          /* Get pointer to the first residue read in from the PDB file */
          r_residue_ptr = r_pdb_ptr->first_r_residue_ptr;
          found = FALSE;

          /* Loop through all the residues to pick up the first residue
             of this ligand */
          while (r_residue_ptr != NULL && found == FALSE)
            {
              /* If we have the first residue, store its pointer */
              if (r_residue_ptr->chain == chain &&
                  !strcmp(r_residue_ptr->res_num,res_num) &&
                  !strcmp(r_residue_ptr->res_name,res_name))
                {
                  /* Save the pointer */
                  r_ligand_ptr->first_r_residue_ptr = r_residue_ptr;

                  /* Set flag */
                  found = TRUE;

                  /* Increment count of ligands assigned */
                  r_pdb_ptr->nligands++;
                }

              /* Get the next r_residue */
              r_residue_ptr = r_residue_ptr->next_r_residue_ptr;
            }
        }
      
      /* Get the next ligand */
      r_ligand_ptr = r_ligand_ptr->next_r_ligand_ptr;
    }

   /* Return number of ligands */
   return(nligands);
}
/***********************************************************************

r_get_pdb_data  -  Open the PDB file and read in the data from it,
                   using the rasmol.dfn template to identify each
                   residue

***********************************************************************/

struct r_pdb *r_get_pdb_data(char *pdb_code,
                             struct c_file_data file_name_ptr,
                             int pdb_type,int data_rqrd)
{
  char file_name[R_FILENAME_LEN + 1], wanted_chain;

  int nligand, nmetal, pdir;
  int gzipped;

  struct r_ligand *first_r_ligand_ptr;
  struct r_pdb *r_pdb_ptr;
  struct r_rasdef *first_r_rasdef_ptr;

#ifdef TEST
  FILE *file_ptr;
#else
  gzFile file_ptr;
#endif

  /* Initialise variables */
  file_ptr = NULL;
  nligand = nmetal = 0;
  r_pdb_ptr = NULL;
  first_r_ligand_ptr = NULL;
  first_r_rasdef_ptr = NULL;
  wanted_chain = '\0';

  /* Determine which PDBsum directory we need */
  if (pdb_type == UPLOAD_PDB || pdb_type == MCSG_PDB)
    pdir = PDBSUM_UPLOAD;
  else
    pdir = PDBSUM_CURRENT;

  /* Read in the definitions from the rasmol.dfn file */
  if (data_rqrd != HEADERS_ONLY)
    r_get_dfn_data(pdb_code,pdb_type,file_name_ptr.pdbsum_template[pdir],
                   &first_r_rasdef_ptr,&first_r_ligand_ptr,&nligand,
                   &nmetal);

  /* Find the PDB file */
  if (pdb_type == PDBSUM1)
    sprintf(file_name,"../%s.new",pdb_code);
  else
    r_find_pdb_file(pdb_code,file_name_ptr.pdb_template,file_name,
                    pdb_type);

  /* Open the PDB file, if found */
  if (file_name[0] != '\0')
    file_ptr = r_open_pdb_file(file_name,&gzipped);

  /* If file found, then read in the data */
  if (file_ptr != NULL)
    {
      /* Read the PDB file */
      r_pdb_ptr = r_read_pdb_file(file_ptr,gzipped,pdb_code,wanted_chain,
                                  data_rqrd,first_r_rasdef_ptr);

      /* If have coordinate data, then process atom records */
      if (data_rqrd != HEADERS_ONLY)
        {
          /* Identify the residues */
          identify_residues(r_pdb_ptr);

          /* Save the first rasdef definition */
          r_pdb_ptr->first_r_rasdef_ptr = first_r_rasdef_ptr;

          /* Save the first ligand definition */
          r_pdb_ptr->first_r_ligand_ptr = first_r_ligand_ptr;

          /* Identify the separate molecules contained in the file and
             calculate their coordinate limits */
          r_pdb_ptr->nmolecs = r_assign_to_molecules(r_pdb_ptr);

          /* Calculate all the residue and molecule coordinate limits */
          r_calc_limits(r_pdb_ptr);

          /* Identify the residues making up each ligand */
          r_assign_to_ligands(r_pdb_ptr);
        }

      /* Close the PDB file */
#ifdef TEST
      fclose(file_ptr);
#else
      gzclose(file_ptr);
#endif
    }

  /* Return the PDB record */
  return(r_pdb_ptr);
}
/***********************************************************************

r_get_pdb_chain  -  Read in just the given chain from the PDB file

***********************************************************************/

struct r_pdb *r_get_pdb_chain(char *pdb_code,char wanted_chain,
                             struct c_file_data file_name_ptr,
                             int pdb_type,int data_rqrd)
{
  char file_name[R_FILENAME_LEN + 1];

  int nligand, nmetal, pdir;
  int gzipped;

  struct r_ligand *first_r_ligand_ptr;
  struct r_pdb *r_pdb_ptr;
  struct r_rasdef *first_r_rasdef_ptr;

#ifdef TEST
  FILE *file_ptr;
#else
  gzFile file_ptr;
#endif

  /* Initialise variables */
  file_ptr = NULL;
  r_pdb_ptr = NULL;
  first_r_rasdef_ptr = NULL;

  /* Determine which PDBsum directory we need */
  if (pdb_type == UPLOAD_PDB || pdb_type == MCSG_PDB)
    pdir = PDBSUM_UPLOAD;
  else
    pdir = PDBSUM_CURRENT;

  /* Read in the definitions from the rasmol.dfn file */
  if (data_rqrd != HEADERS_ONLY)
    r_get_dfn_data(pdb_code,pdb_type,file_name_ptr.pdbsum_template[pdir],
                   &first_r_rasdef_ptr,&first_r_ligand_ptr,&nligand,
                   &nmetal);

  /* Find the PDB file */
  r_find_pdb_file(pdb_code,file_name_ptr.pdb_template,file_name,
                  pdb_type);

  /* Open the PDB file, if found */
  if (file_name[0] != '\0')
    file_ptr = r_open_pdb_file(file_name,&gzipped);

  /* If file found, then read in the data */
  if (file_ptr != NULL)
    {
      /* Read the PDB file */
      r_pdb_ptr = r_read_pdb_file(file_ptr,gzipped,pdb_code,wanted_chain,
                                  data_rqrd,first_r_rasdef_ptr);

      /* If have coordinate data, then process atom records */
      if (data_rqrd != HEADERS_ONLY)
        {
          /* Identify the residues */
          identify_residues(r_pdb_ptr);

          /* Save the first rasdef definition */
          r_pdb_ptr->first_r_rasdef_ptr = first_r_rasdef_ptr;

          /* Save the first ligand definition */
          r_pdb_ptr->first_r_ligand_ptr = first_r_ligand_ptr;

          /* Identify the spearate molecules contained in the file and
             calculate their coordinate limits */
          r_pdb_ptr->nmolecs = r_assign_to_molecules(r_pdb_ptr);

          /* Calculate all the residue and molecule coordinate limits */
          r_calc_limits(r_pdb_ptr);

          /* Identify the residues making up each ligand */
          r_assign_to_ligands(r_pdb_ptr);
        }

      /* Close the PDB file */
#ifdef TEST
      fclose(file_ptr);
#else
      gzclose(file_ptr);
#endif
    }

  /* Return the PDB record */
  return(r_pdb_ptr);
}
/***********************************************************************

r_get_res_name  -  Get amino acid name from single-letter code

***********************************************************************/

void r_get_res_name(char aa_code,char *res_name)
{
  char *string_ptr;

  int iamino;
  int namino;

  static char amino[] = "ACDEFGHIKLMNPQRSTVWYXBZ";
  static char *amino_name[] = {
    "ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU",
    "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR",
    "UNK ","ASX ","GLX "
  };

  /* Initialise variables */
  namino = strlen(amino) - 1;

  /* Locate the single-letter amino acid code */
  string_ptr = strchr(amino,aa_code);

  /* If found, save corresponding residue name */
  if (string_ptr != NULL)
    {
      /* Get which code this is */
      iamino = string_ptr - amino;
      strcpy(res_name,amino_name[iamino]);
    }
  else
    strcpy(res_name,"UNK");
}
/***********************************************************************

get_het_desc  -  Get Het Group description for given residue

***********************************************************************/

void get_het_desc(struct r_residue *r_residue_ptr,
                  struct r_hetdesc *first_r_hetdesc_ptr)
{
  int found;

  struct r_hetdesc *r_hetdesc_ptr;

  /* Initialise variables */
  found = FALSE;

  /* Get pointer to the first hetdesc record */
  r_hetdesc_ptr = first_r_hetdesc_ptr;

  /* Loop over all the hetdesc records */
  while (r_hetdesc_ptr != NULL && found == FALSE)
    {
      /* Check residue name for match */
      if (!strcmp(r_hetdesc_ptr->res_name,r_residue_ptr->res_name))
        {
          /* Store the pounter to this hetdesc record */
          r_residue_ptr->r_hetdesc_ptr = r_hetdesc_ptr;

          /* Set flag that we're found */
          found = TRUE;
        }

      /* Get the next hetdesc record */
      r_hetdesc_ptr = r_hetdesc_ptr->next_r_hetdesc_ptr;
    }
}
/***********************************************************************

update_residues -  Update chain pointers for all residues in the current
                   chain

***********************************************************************/

void update_residues(struct r_chain *r_chain_ptr)
{
  int iresid, nresid;

  struct r_residue *r_residue_ptr;

  /* Get the first of this chain's residues */
  r_residue_ptr = r_chain_ptr->first_r_residue_ptr;
  nresid = r_chain_ptr->nresid;
  iresid = 0;

  /* Loop through the chain's residues, determining to which
     other residues in this chain they are covalently bonded to */
  while (r_residue_ptr != NULL && iresid < nresid)
    {
      /* Set chain pointer */
      r_residue_ptr->r_chain_ptr = r_chain_ptr;

      /* Get the next residue */
      r_residue_ptr = r_residue_ptr->next_r_residue_ptr;
      iresid++;
    }
}
/***********************************************************************

excise_metal  -  Place metal ion as a separate chain

***********************************************************************/

void excise_metal(struct r_pdb *r_pdb_ptr,struct r_residue *r_residue_ptr,
                  struct r_chain **lst_r_chain_ptr)
{
  int iresid, nresid, nbefore;

  struct r_chain *r_chain_ptr, *new_r_chain_ptr;
  struct r_chain *last_r_chain_ptr;
  struct r_chain *first_dummy_r_chain_ptr, *last_dummy_r_chain_ptr;
  struct r_residue *other_r_residue_ptr;

  /* Initialise end chain pointer */
  last_r_chain_ptr = *lst_r_chain_ptr;

  /* Get this residue's chain */
  r_chain_ptr = r_residue_ptr->r_chain_ptr;

  /* Define its chain as a metal chain */
  r_chain_ptr->type = R_CHAIN_METAL;

  /* If it is not the first residue of the chain start a new chain
     for it */
  if (r_residue_ptr != r_chain_ptr->first_r_residue_ptr)
    {
      /* Initialise count of preceding residues */
      other_r_residue_ptr = r_chain_ptr->first_r_residue_ptr;
      nresid = r_chain_ptr->nresid;
      iresid = 0;
      nbefore = 0;

      /* Loop through the residues preceding it in the current chain
         to determine how many there are */
      while (other_r_residue_ptr != NULL &&
             other_r_residue_ptr != r_residue_ptr && iresid < nresid)
        {
          /* Increment count */
          nbefore++;

          /* Get the next residue */
          other_r_residue_ptr = other_r_residue_ptr->next_r_residue_ptr;
          iresid++;
        }

      /* Initialise dummy pointers */
      first_dummy_r_chain_ptr = NULL;
      last_dummy_r_chain_ptr = NULL;

      /* Create a new chain record */
      new_r_chain_ptr = create_r_chain_record(&first_dummy_r_chain_ptr,
                                              &last_dummy_r_chain_ptr,
                                              r_pdb_ptr,
                                              r_chain_ptr->chain_id);

      /* Store current residue as the new chain's first */
      new_r_chain_ptr->first_r_residue_ptr = r_residue_ptr;

      /* Create links to the new chain */
      new_r_chain_ptr->next_r_chain_ptr = r_chain_ptr->next_r_chain_ptr;
      r_chain_ptr->next_r_chain_ptr = new_r_chain_ptr;

      /* If the current chain is the last in the linked list, make
         the new chain the new last */
      if (last_r_chain_ptr == r_chain_ptr)
        last_r_chain_ptr = new_r_chain_ptr;

      /* Update first residues and residue counts */
      new_r_chain_ptr->nresid = r_chain_ptr->nresid - nbefore;
      r_chain_ptr->nresid = nbefore;

      /* Copy across other parameters to the new chain */
      new_r_chain_ptr->type = r_chain_ptr->type;

      /* Update all this chain's residues so that they point to the
         new chain */
      update_residues(new_r_chain_ptr);

      /* Store new chain pointer */
      r_chain_ptr = new_r_chain_ptr;
    }

  /* Split off any residues after the current one, creating a new
     chain */
  if (r_chain_ptr->nresid > 1)
    {
      /* Initialise dummy pointers */
      first_dummy_r_chain_ptr = NULL;
      last_dummy_r_chain_ptr = NULL;

      /* Create a new chain record */
      new_r_chain_ptr = create_r_chain_record(&first_dummy_r_chain_ptr,
                                              &last_dummy_r_chain_ptr,
                                              r_pdb_ptr,
                                              r_chain_ptr->chain_id);

      /* Store residue after current residue as the new chain's first */
      new_r_chain_ptr->first_r_residue_ptr
        = r_residue_ptr->next_r_residue_ptr;

      /* Create links to the new chain */
      new_r_chain_ptr->next_r_chain_ptr = r_chain_ptr->next_r_chain_ptr;
      r_chain_ptr->next_r_chain_ptr = new_r_chain_ptr;

      /* If the current chain is the last in the linked list, make
         the new chain the new last */
      if (last_r_chain_ptr == r_chain_ptr)
        last_r_chain_ptr = new_r_chain_ptr;

      /* Update first residues and residue counts */
      new_r_chain_ptr->nresid = r_chain_ptr->nresid - 1;
      r_chain_ptr->nresid = 1;

      /* Copy across other parameters to the new chain */
      new_r_chain_ptr->type = r_chain_ptr->type;
      r_chain_ptr->type = R_CHAIN_METAL;

      /* Update all this chain's residues so that they point to the
         new chain */
      update_residues(new_r_chain_ptr);
    }

  /* Return last chain pointers */
  *lst_r_chain_ptr = last_r_chain_ptr;
}
/***********************************************************************

classify_residues  -  Classify all the residues according to type

***********************************************************************/

void classify_residues(struct r_pdb *r_pdb_ptr,
                       struct r_chain **lst_r_chain_ptr)
{
  int i;

  struct r_chain *last_r_chain_ptr;
  struct r_residue *r_residue_ptr;

  /* Initialise end chain pointer */
  last_r_chain_ptr = *lst_r_chain_ptr;

  /* Get pointer to the first residue */
  r_residue_ptr = r_pdb_ptr->first_r_residue_ptr;

  /* Loop through all the residues to classify */
  while (r_residue_ptr != NULL)
    {
      /* Get the residue type and amino acid code from its name */
      get_residue_type(r_residue_ptr);

      /* If residue non-standard, get its description */
      if (r_residue_ptr->type == R_TYPE_NON_STANDARD ||
          r_residue_ptr->type == R_TYPE_METAL_ION ||
          r_residue_ptr->type == R_TYPE_UNKNOWN)
        get_het_desc(r_residue_ptr,r_pdb_ptr->first_r_hetdesc_ptr);

      /* If this residue is a metal ion, then want to place it in a
         separate chain of its own */
      if (r_residue_ptr->type == R_TYPE_METAL_ION)
        excise_metal(r_pdb_ptr,r_residue_ptr,&last_r_chain_ptr);

      /* Get the next residue */
      r_residue_ptr = r_residue_ptr->next_r_residue_ptr;
    }

  /* Return last chain pointers */
  *lst_r_chain_ptr = last_r_chain_ptr;
}
/***********************************************************************

classify_chains  -  Initial analysis of each chain to determine what it
                    might be (protein/DNA/ligand)

***********************************************************************/

int classify_chains(struct r_chain *first_r_chain_ptr)
{
  int iatom, iresid, last_sec_struc, niatoms, natoms, nchains, nresid;
  int have_non_blank;

  struct r_atom *r_atom_ptr;
  struct r_chain *r_chain_ptr;
  struct r_residue *r_residue_ptr;

  /* Initialise variables */
  have_non_blank = FALSE;
  last_sec_struc = R_COIL;
  nchains = 0;

  /* Get pointer to the first chain */
  r_chain_ptr = first_r_chain_ptr;

  /* Loop over all the currently defined chains */
  while (r_chain_ptr != NULL)
    {
      /* If this is a non-blank chain identifier, then set flag */
      if (r_chain_ptr->chain_id != ' ')
        have_non_blank = TRUE;

      /* Initialise counts */
      r_chain_ptr->nhelices = 0;
      r_chain_ptr->nstrands = 0;
      r_chain_ptr->nstand = 0;
      r_chain_ptr->nbase = 0;
      natoms = 0;

      /* Get this chain's first residue */
      r_residue_ptr = r_chain_ptr->first_r_residue_ptr;
      nresid = r_chain_ptr->nresid;
      iresid = 0;

      /* Initialise chain counts */
      r_chain_ptr->natom_recs = 0;
      r_chain_ptr->nhetatm_recs = 0;

      /* Loop through this chain's residues, updating all the counts */
      while (r_residue_ptr != NULL && iresid < nresid)
        {
          /* Increment total atom count for this chain */
          natoms = natoms + r_residue_ptr->natoms;

          /* Get the first of the first residue's atoms */
          r_atom_ptr = r_residue_ptr->first_r_atom_ptr;
          iatom = 0;
          niatoms = r_residue_ptr->natoms;
  
          /* Loop over the residue's atoms to update the chain counts */
          while (r_atom_ptr != NULL && iatom < niatoms)
            {
              /* Update ATOM or HETATM count */
              if (r_atom_ptr->rec_type == R_ATOM_REC)
                r_chain_ptr->natom_recs++;
              else
                r_chain_ptr->nhetatm_recs++;

              /* Get pointer to next atom in linked-list */
              r_atom_ptr = r_atom_ptr->next_r_atom_ptr;
              iatom++;
            }

          /* Update chain counts depending on this residue's type */
          if (r_residue_ptr->type == R_TYPE_STANDARD ||
              r_residue_ptr->type == R_TYPE_NON_STANDARD)
            r_chain_ptr->nstand++;
          else if (r_residue_ptr->type == R_TYPE_NUCLEIC_ACID)
            r_chain_ptr->nbase++;

          /* Update secondary structure counts */
          if (r_residue_ptr->sec_struc != last_sec_struc)
            {
              /* Update chain's helix and strand counts */
              if (r_residue_ptr->sec_struc == R_HELIX)
                r_chain_ptr->nhelices++;
              else if (r_residue_ptr->sec_struc == R_STRAND)
                r_chain_ptr->nstrands++;

              /* Save the current secondary structure */
              last_sec_struc = r_residue_ptr->sec_struc;
            }

          /* Get the next residue */
          r_residue_ptr = r_residue_ptr->next_r_residue_ptr;
          iresid++;
        }

      /* Consider non-metal chains */
      if (r_chain_ptr->type != R_CHAIN_METAL)
        {
          /* Check for obvious ligands: no ATOM records, only HETATM ones */
          if (r_chain_ptr->natom_recs == 0)
            r_chain_ptr->type = R_CHAIN_LIGAND;

          /* If chain has a blank chain-id, and we already have see non-blank
             chain identifiers, then this is not a protein */
          else if (r_chain_ptr->chain_id == ' ' && have_non_blank == TRUE &&
                   r_chain_ptr->nhelices + r_chain_ptr->nstrands == 0)
            r_chain_ptr->type = R_CHAIN_LIGAND;

          /* If chain contains nucleic acids, then mark down as a nucleic
             acid chain */
          else if (r_chain_ptr->nbase > 0)
            r_chain_ptr->type = R_CHAIN_DNA;

          /* Otherwise, if contains standard, or unusual, amino acids, then
             mark down as a polypeptide (or all-CA) chain */
          else if (r_chain_ptr->nstand > 0)
            {
              /* Determine whether this looks like an all-CA chain */
              if (natoms < 2 * r_chain_ptr->nresid)
                r_chain_ptr->type = R_CHAIN_CALPHA_ONLY;

              /* Otherwise, mark this chain down as a polypeptide chain */
              else
                r_chain_ptr->type = R_CHAIN_PROTEIN;
            }

          /* If no standard amino acids, and no nucleic acids, then this
             must be a ligand chain */
          else
            r_chain_ptr->type = R_CHAIN_LIGAND;
        }

      /* Increment count of chains */
      nchains++;

      /* Get the next chain */
      r_chain_ptr = r_chain_ptr->next_r_chain_ptr;
    }

  /* Return number of chains */
  return(nchains);
}
/***********************************************************************

shift_metals  -  Shift the chains representing metals to the end of the
                 list

***********************************************************************/

int shift_metals(struct r_chain **fst_r_chain_ptr,
                 struct r_chain **lst_r_chain_ptr)
{
  int nshift;

  struct r_chain *r_chain_ptr, *first_r_chain_ptr, *last_r_chain_ptr;
  struct r_chain *done_r_chain_ptr, *next_r_chain_ptr, *prev_r_chain_ptr;

  /* Initialise variables */
  done_r_chain_ptr = NULL;
  next_r_chain_ptr = NULL;
  prev_r_chain_ptr = NULL;
  nshift = 0;

  /* Initialise start and end chain pointers */
  first_r_chain_ptr = *fst_r_chain_ptr;
  last_r_chain_ptr = *lst_r_chain_ptr;

  /* Get pointer to the first chain */
  r_chain_ptr = first_r_chain_ptr;

  /* Loop over all the chains */
  while (r_chain_ptr != done_r_chain_ptr)
    {
      /* Get the next chain */
      next_r_chain_ptr = r_chain_ptr->next_r_chain_ptr;

      /* If this is a metal chain, shift to end of linked list */
      if (r_chain_ptr->type == R_CHAIN_METAL)
        {
          /* If have a previous chain, make its next pointer point
             around the metal */
          if (prev_r_chain_ptr != NULL)
            prev_r_chain_ptr->next_r_chain_ptr = next_r_chain_ptr;

          /* Otherwise, this must be the first chain, so adjust
             pointer to first chain */
          else
            first_r_chain_ptr = next_r_chain_ptr;

          /* Place current metal chain at end of linked list */
          if (last_r_chain_ptr != NULL)
            {
              /* Change pointers */
              last_r_chain_ptr->next_r_chain_ptr = r_chain_ptr;
              r_chain_ptr->next_r_chain_ptr = NULL;

              /* Make current chain the new last chain */
              last_r_chain_ptr = r_chain_ptr;
            }

          /* If this is the first shifted metal, save its pointer */
          if (done_r_chain_ptr == NULL)
            done_r_chain_ptr = r_chain_ptr;

          /* Increment count of chains shifted */
          nshift++;
        }

      /* If not metal, save as previous chain */
      else
        prev_r_chain_ptr = r_chain_ptr;

      /* Get the next chain */
      r_chain_ptr = next_r_chain_ptr;
    }

  /* Return current chain pointers */
  *fst_r_chain_ptr = first_r_chain_ptr;
  *lst_r_chain_ptr = last_r_chain_ptr;

  /* Return number of chains shifted */
  return(nshift);
}
/***********************************************************************

adjust_residues  -  Adjust the residues so that their order corresponds
                    to the current chain order

***********************************************************************/

int adjust_residues(struct r_residue **fst_r_residue_ptr,
                    struct r_chain *first_r_chain_ptr)
{
  int iresid, nres, nresid;

  struct r_chain *r_chain_ptr;
  struct r_residue *r_residue_ptr, *first_r_residue_ptr, *last_r_residue_ptr;

  /* Initialise variables */
  first_r_residue_ptr = NULL;
  last_r_residue_ptr = NULL;
  nres = 0;

  /* Get pointer to the first chain */
  r_chain_ptr = first_r_chain_ptr;

  /* Loop over the chains */
  while (r_chain_ptr != NULL)
    {
      /* Get this chain's first residue */
      r_residue_ptr = r_chain_ptr->first_r_residue_ptr;
      nresid = r_chain_ptr->nresid;
      iresid = 0;
  
      /* If this is the very first residue, store it */
      if (first_r_residue_ptr == NULL)
        first_r_residue_ptr = r_residue_ptr;

      /* If have a last residue, get it to point to this one */
      if (last_r_residue_ptr != NULL)
        last_r_residue_ptr->next_r_residue_ptr = r_residue_ptr;

      /* Initialise pointer to this chain's last residue */
      last_r_residue_ptr = NULL;

      /* Loop through the chain's residues to find the last */
      while (r_residue_ptr != NULL && iresid < nresid)
        {
          /* Update residue pointer */
          last_r_residue_ptr = r_residue_ptr;

          /* Get the next residue */
          r_residue_ptr = r_residue_ptr->next_r_residue_ptr;
          iresid++;
          nres++;
        }

      /* Get the next chain */
      r_chain_ptr = r_chain_ptr->next_r_chain_ptr;
    }

  /* Make last residue terminate the linked list */
  if (last_r_residue_ptr != NULL)
    last_r_residue_ptr->next_r_residue_ptr = NULL;

  /* Return current residue pointer */
  *fst_r_residue_ptr = first_r_residue_ptr;

  /* Return number of residues */
  return(nres);
}
/***********************************************************************

join_chains  -  Join any adjacent chains which seem to belong to the
                same polypeptide

***********************************************************************/

void join_chains(struct r_chain *first_r_chain_ptr)
{
  int iresid, nresid;

  struct r_chain *r_chain_ptr, *next_r_chain_ptr;
  struct r_residue *r_residue_ptr, *last_r_residue_ptr;

  /* Get pointer to the first chain */
  r_chain_ptr = first_r_chain_ptr;
  if (r_chain_ptr != NULL)
    next_r_chain_ptr = r_chain_ptr->next_r_chain_ptr;

  /* Loop over all the currently defined chains */
  while (r_chain_ptr != NULL && next_r_chain_ptr != NULL)
    {
      /* If both chains have the same chain-id and are of the same
         type (ie protein, DNA or ligand) then join together. The
         second chain needs to be longer than standard ligand length */
      if (r_chain_ptr->chain_id == next_r_chain_ptr->chain_id &&
          ((r_chain_ptr->type == next_r_chain_ptr->type &&
            next_r_chain_ptr->nresid > R_LIGAND_LENGTH) ||
           (r_chain_ptr->type == R_CHAIN_PROTEIN &&
            next_r_chain_ptr->type == R_CHAIN_CALPHA_ONLY)))
        {
          /* Get the first of the second chain's residues */
          r_residue_ptr = next_r_chain_ptr->first_r_residue_ptr;
          nresid = next_r_chain_ptr->nresid;
          iresid = 0;

          /* Loop through the chain's residues, marking them as now
             belonging to the first chain */
          while (r_residue_ptr != NULL && iresid < nresid)
            {
              /* Update chain pointer */
              r_residue_ptr->r_chain_ptr = r_chain_ptr;

              /* Get the next residue */
              r_residue_ptr = r_residue_ptr->next_r_residue_ptr;
              iresid++;
            }

          /* Join the second chain onto the first */
          r_chain_ptr->nresid
            = r_chain_ptr->nresid + next_r_chain_ptr->nresid;
          r_chain_ptr->nhelices
            = r_chain_ptr->nhelices + next_r_chain_ptr->nhelices;
          r_chain_ptr->nstrands
            = r_chain_ptr->nstrands + next_r_chain_ptr->nstrands;
          r_chain_ptr->nstand
            = r_chain_ptr->nstand + next_r_chain_ptr->nstand;
          r_chain_ptr->nbase
            = r_chain_ptr->nbase + next_r_chain_ptr->nbase;

          /* Get elongated chain to point past the removed second
             fragment */
          r_chain_ptr->next_r_chain_ptr = next_r_chain_ptr->next_r_chain_ptr;

          /* Free up memory from the discarded chain */
          free(next_r_chain_ptr);

          /* Get the new next chain */
          next_r_chain_ptr = r_chain_ptr->next_r_chain_ptr;
        }

      /* Otherwise, go on to the next chain */
      else
        {
          r_chain_ptr = next_r_chain_ptr;
          if (r_chain_ptr != NULL)
            next_r_chain_ptr = r_chain_ptr->next_r_chain_ptr;
        }
    }
}
/***********************************************************************

reclassify_chains  -  Reclassify very short polypeptide chains or
                      fragments of DNA as ligands

***********************************************************************/

void reclassify_chains(struct r_chain *first_r_chain_ptr)
{
  struct r_chain *r_chain_ptr;

  /* Get pointer to the first chain */
  r_chain_ptr = first_r_chain_ptr;

  /* Loop over all the currently defined chains */
  while (r_chain_ptr != NULL)
    {
      /* If short polypeptide chain, class as a ligand - unless it is
         the very first chain */
      if (r_chain_ptr->type == R_CHAIN_PROTEIN)
        {
          /* If chain shorter than ligand-length cutoff, reclassify as
             ligand */
          if (r_chain_ptr->nresid <= R_LIGAND_LENGTH)
            r_chain_ptr->type = R_CHAIN_LIGAND;
        }

      /* Similarly, a fragment of DNA of fewer than 3 bases can be
         classed a ligand */
      else if (r_chain_ptr->type == R_CHAIN_DNA && r_chain_ptr->nresid < 3)
        r_chain_ptr->type = R_CHAIN_LIGAND;

      /* Get the next chain */
      r_chain_ptr = r_chain_ptr->next_r_chain_ptr;
    }
}
/***********************************************************************

check_for_bond  -  Determine whether the given two residues are
                   covalently bonded

***********************************************************************/

int check_for_bond(struct r_residue *r_residue_ptr,
                   struct r_residue *other_r_residue_ptr)
{
  int i, iatom, niatoms, jatom, njatoms;
  int bonded, toofar;

  double dist_sqrd;

  struct r_atom *r_atom_ptr, *other_r_atom_ptr;
  struct r_rbond *r_rbond_ptr;

  /* Initialise variables */
  bonded = FALSE;

  /* Get the first of the residue's CONECT records */
  r_rbond_ptr = r_residue_ptr->first_r_rbond_ptr;

  /* Check for a bond defined by the CONECT records */
  while (r_rbond_ptr != NULL && bonded == FALSE)
    {
      /* If this bond joins the current two residues, then we're done */
      if (r_rbond_ptr->to_r_residue_ptr == other_r_residue_ptr)
        bonded = TRUE;

      /* Get the next residue CONECT record */
      r_rbond_ptr = r_rbond_ptr->next_r_rbond_ptr;
    }

  /* If no existing bond found, then check residue limits to see if
     residues might be within covalent bonding distance */
  if (bonded == FALSE)
    {
      /* Check whether these residue are within reach */
      toofar = FALSE;
      for (i = 0; i < 3 && toofar == FALSE; i++)
        {
          if (r_residue_ptr->coords_min[i]
              - other_r_residue_ptr->coords_max[i] > R_TOO_FAR_DIST)
            toofar = TRUE;
          if (other_r_residue_ptr->coords_min[i]
              - r_residue_ptr->coords_max[i] > R_TOO_FAR_DIST)
            toofar = TRUE;
        }
      
      /* If residues within range calculate all atom-atom distances */
      if (toofar == FALSE)
        {
          /* Get the first of the first residue's atoms */
          r_atom_ptr = r_residue_ptr->first_r_atom_ptr;
          iatom = 0;
          niatoms = r_residue_ptr->natoms;
  
          /* Loop through all the residue's atoms */
          while (r_atom_ptr != NULL && iatom < niatoms && bonded == FALSE)
            {
              /* Check whether this atom within range of the other
                 residue */
              toofar = FALSE;
              for (i = 0; i < 3 && toofar == FALSE; i++)
                {
                  if (r_atom_ptr->coords[i]
                      - other_r_residue_ptr->coords_max[i] > R_TOO_FAR_DIST ||
                      other_r_residue_ptr->coords_min[i]
                      - r_atom_ptr->coords[i] > R_TOO_FAR_DIST)
                    toofar = TRUE;
                }

              /* If atom not too far from the other residue, check its
                 distance from all that residue's atoms */
              if (toofar == FALSE)
                {
                  /* Get the first of the other residue's atoms */
                  other_r_atom_ptr = other_r_residue_ptr->first_r_atom_ptr;
                  jatom = 0;
                  njatoms = other_r_residue_ptr->natoms;
  
                  /* Loop through all the residue's atoms */
                  while (jatom < njatoms && other_r_atom_ptr != NULL &&
                         bonded == FALSE)
                    {
                      /* Calculate the squared distance between these
                         two atoms */
                      dist_sqrd = 0.0;
                      for (i = 0; i < 3; i++)
                        dist_sqrd
                          = dist_sqrd
                          + (r_atom_ptr->coords[i]
                             - other_r_atom_ptr->coords[i])
                          * (r_atom_ptr->coords[i]
                             - other_r_atom_ptr->coords[i]);

                      /* If within covalent bonding distance, have a
                         bond between the two residues */
                      if (dist_sqrd < R_COVALENT_DIST_SQRD)
                        {
                          /* Create a bond between these two residues */
                          create_r_rbond_record(r_residue_ptr,
                                                other_r_residue_ptr);

                          /* Set flag that residues bonded */
                          bonded = TRUE;
                        }

                      /* Get pointer to next atom in linked-list */
                      other_r_atom_ptr = other_r_atom_ptr->next_r_atom_ptr;
                      jatom++;
                    }
                }

              /* Get pointer to next atom in linked-list */
              r_atom_ptr = r_atom_ptr->next_r_atom_ptr;
              iatom++;
            }
        }
    }

  /* Return whether residues are bonded */
  return(bonded);
}
/***********************************************************************

find_joins  -  Determine which residues in this chain are connected
               by covalent bonds

***********************************************************************/

void find_joins(struct r_chain *r_chain_ptr)
{
  int iresid, jresid, nresid;
  int bonded;

  struct r_residue *r_residue_ptr, *other_r_residue_ptr;

  /* Get the first of this chain's residues */
  r_residue_ptr = r_chain_ptr->first_r_residue_ptr;
  nresid = r_chain_ptr->nresid;
  iresid = 0;

  /* Loop through the chain's residues, determining to which
     other residues in this chain they are covalently bonded to */
  while (r_residue_ptr != NULL && iresid < nresid)
    {
      /* Get pointer to the next residue */
      other_r_residue_ptr = r_residue_ptr->next_r_residue_ptr;
      jresid = iresid + 1;

      /* Loop through each of the chain's other residues
         to calculate residue connectivities */
      while (other_r_residue_ptr != NULL && jresid < nresid)
        {
          /* Determine whether these two residues are covalently
             bonded */
          bonded = check_for_bond(r_residue_ptr,other_r_residue_ptr);

          /* Get the next residue */
          other_r_residue_ptr = other_r_residue_ptr->next_r_residue_ptr;
          jresid++;
        }

      /* Get the next residue */
      r_residue_ptr = r_residue_ptr->next_r_residue_ptr;
      iresid++;
    }
}
/***********************************************************************

create_r_list_record  -  Create and initialise a new list record

***********************************************************************/

struct r_list *create_r_list_record(struct r_list **fst_r_list_ptr,
                                    struct r_list **lst_r_list_ptr,
                                    struct r_residue *r_residue_ptr)
{
  struct r_list *r_list_ptr, *first_r_list_ptr, *last_r_list_ptr;

  /* Initialise list pointers */
  r_list_ptr = NULL;
  first_r_list_ptr = *fst_r_list_ptr;
  last_r_list_ptr = *lst_r_list_ptr;

  /* Allocate memory for structure to hold list info */
  r_list_ptr = (struct r_list *) malloc(sizeof(struct r_list));
  if (r_list_ptr == NULL)
    {
      printf("*** Can't allocate memory for struct r_list\n");
      exit (1);
    }

  /* If this is the very first list, then store pointer */
  if (first_r_list_ptr == NULL)
    first_r_list_ptr = r_list_ptr;

  /* Add link from previous list to the current one */
  if (last_r_list_ptr != NULL)
    last_r_list_ptr->next_r_list_ptr = r_list_ptr;
  last_r_list_ptr = r_list_ptr;

  /* Store the list details */
  r_list_ptr->r_residue_ptr = r_residue_ptr;
  r_list_ptr->next_r_list_ptr = NULL;

  /* Return current pointers */
  *fst_r_list_ptr = first_r_list_ptr;
  *lst_r_list_ptr = last_r_list_ptr;

  /* Return new list pointer */
  return(r_list_ptr);
}
/***********************************************************************

free_r_list_memory  -  Free up memory used by residue list

***********************************************************************/

void free_r_list_memory(struct r_list *first_r_list_ptr)
{
  struct r_list *r_list_ptr, *next_r_list_ptr;

  /* Get first list */
  r_list_ptr = first_r_list_ptr;

  /* Loop through the list records until done */
  while (r_list_ptr != NULL)
    {
      /* Get pointer to the next list record */
      next_r_list_ptr = r_list_ptr->next_r_list_ptr;

      /* Free up memory taken by current list */
      free(r_list_ptr);

      /* Go to the next list */
      r_list_ptr = next_r_list_ptr;
    }
}
/***********************************************************************

identify_molecule  -  Identify all residues in this chain that are in
                      some way connected to the first residue

***********************************************************************/

void identify_molecule(struct r_chain *r_chain_ptr)
{
  int iresid, nresid;

  struct r_list *r_list_ptr;
  struct r_list *first_r_list_ptr, *last_r_list_ptr;
  struct r_rbond *r_rbond_ptr;
  struct r_residue *r_residue_ptr, *to_r_residue_ptr;

  /* Initialise pointers */
  r_list_ptr = first_r_list_ptr = last_r_list_ptr = NULL;

  /* Get the first of this chain's residues */
  r_residue_ptr = r_chain_ptr->first_r_residue_ptr;
  nresid = r_chain_ptr->nresid;
  iresid = 0;

  /* Loop through the chain's residues to initialise connected flag */
  while (r_residue_ptr != NULL && iresid < nresid)
    {
      /* Initialise flag */
      r_residue_ptr->connected = FALSE;

      /* Get the next residue */
      r_residue_ptr = r_residue_ptr->next_r_residue_ptr;
      iresid++;
    }

  /* Get the first residue again */
  r_residue_ptr = r_chain_ptr->first_r_residue_ptr;
  iresid = 0;
  
  /* Create a list record for the first residue */
  r_list_ptr
    = create_r_list_record(&first_r_list_ptr,&last_r_list_ptr,
                           r_residue_ptr);

  /* Loop while there are list records on the stack */
  while (r_list_ptr != NULL)
    {
      /* Get the residue stored */
      r_residue_ptr = r_list_ptr->r_residue_ptr;

      /* Process only if not already marked */
      if (r_residue_ptr->connected == FALSE)
        {
          /* Mark this residue as connected */
          r_residue_ptr->connected = TRUE;

          /* Get its first residue bond */
          r_rbond_ptr = r_residue_ptr->first_r_rbond_ptr;

          /* Loop through the list to create a list record for each */
          while (r_rbond_ptr != NULL)
            {
              /* Get the connected residue */
              to_r_residue_ptr = r_rbond_ptr->to_r_residue_ptr;

              /* If this residue not already connected, add to list */
              if (to_r_residue_ptr->connected == FALSE &&
                  to_r_residue_ptr->r_chain_ptr == r_chain_ptr)
                {
                  /* Create a list record for the residue to which the current
                     one is bonded */
                  create_r_list_record(&first_r_list_ptr,&last_r_list_ptr,
                                     r_rbond_ptr->to_r_residue_ptr);
                }

              /* Get the next residue bond */
              r_rbond_ptr = r_rbond_ptr->next_r_rbond_ptr;
            }
        }

      /* Get the next list record */
      r_list_ptr = r_list_ptr->next_r_list_ptr;
    }

  /* Free up memory taken by the list pointers */
  free_r_list_memory(first_r_list_ptr);
}
/***********************************************************************

split_ligands  -  Calculate the connectivities between all residues in
                  each ligand chain to determine which are phyically
                  distinct molecules

***********************************************************************/

void split_ligands(struct r_chain *first_r_chain_ptr)
{
  int iresid, nresid, start_nresid;
  int done;

  struct r_chain *r_chain_ptr, *new_r_chain_ptr;
  struct r_chain *first_new_r_chain_ptr, *last_new_r_chain_ptr;
  struct r_residue *r_residue_ptr, *start_r_residue_ptr;

  /* Get pointer to the first chain */
  r_chain_ptr = first_r_chain_ptr;

  /* Loop over all the currently defined chains */
  while (r_chain_ptr != NULL)
    {
      /* Only process if this is a ligand chain */
      if (r_chain_ptr->type == R_CHAIN_LIGAND)
        {
          /* Determine which residues in the chain are connected by
             covalent bonds */
          find_joins(r_chain_ptr);

          /* Initialise flag */
          done = FALSE;

          /* Loop while splitting the chain up */
          while (done == FALSE)
            {
              /* Mark all residues that can be reached from the first
                 residue */
              identify_molecule(r_chain_ptr);

              /* Get the first of this chain's residues */
              r_residue_ptr = r_chain_ptr->first_r_residue_ptr;
              nresid = r_chain_ptr->nresid;
              iresid = 0;
              start_r_residue_ptr = NULL;

              /* Loop through the chain's residues to retain all those
                 flagged as connected */
              while (r_residue_ptr != NULL && iresid < nresid &&
                     start_r_residue_ptr == NULL)
                {
                  /* If not connected, need to terminate here */
                  if (r_residue_ptr->connected == FALSE)
                    {
                      /* Save this residue as being the first of
                         a new chain */
                      start_r_residue_ptr = r_residue_ptr;

                      /* Calculate how many residues the new chain
                         will have */
                      start_nresid = nresid - iresid;

                      /* Save the number of residues the truncated
                         chain has */
                      r_chain_ptr->nresid = iresid;
                    }

                  /* Get the next residue */
                  r_residue_ptr = r_residue_ptr->next_r_residue_ptr;
                  iresid++;
                }

              /* If any residues left over, create a new chain for them
                 and transfer residues to it */
              if (start_r_residue_ptr != NULL)
                {
                  /* Initialise pointers */
                  first_new_r_chain_ptr = NULL;
                  last_new_r_chain_ptr = NULL;

                  /* Create a new chain record */
                  new_r_chain_ptr
                    = create_r_chain_record(&first_new_r_chain_ptr,
                                            &last_new_r_chain_ptr,
                                            r_chain_ptr->r_pdb_ptr,
                                            r_chain_ptr->chain_id);

                  /* Store the chain's first residue and number of
                     residues */
                  new_r_chain_ptr->first_r_residue_ptr = start_r_residue_ptr;
                  new_r_chain_ptr->nresid = start_nresid;
                  new_r_chain_ptr->type = r_chain_ptr->type;

                  /* Update all this chain's residues so that they
                     point to the new chain */
                  update_residues(new_r_chain_ptr);

                  /* Adjust chain pointers to insert the new chain
                     into the linked list */
                  new_r_chain_ptr->next_r_chain_ptr
                    = r_chain_ptr->next_r_chain_ptr;
                  r_chain_ptr->next_r_chain_ptr = new_r_chain_ptr;
                  r_chain_ptr = new_r_chain_ptr;
                }

              /* Otherwise, if none left, then we are done */
              else
                done = TRUE;
            }
        }

      /* Get the next chain */
      r_chain_ptr = r_chain_ptr->next_r_chain_ptr;
    }
}
/***********************************************************************

single_residue_ligands  -  Identify which single-residue ligands are
                           metals and which are actually residue
                           attachments

***********************************************************************/

void single_residue_ligands(struct r_chain *first_r_chain_ptr)
{
  int iresid, itype, nresid;
  int bonded, found, metal;

  struct r_atom *r_atom_ptr;
  struct r_chain *r_chain_ptr, *other_r_chain_ptr;
  struct r_residue *r_residue_ptr, *other_r_residue_ptr;

  /* Get pointer to the first chain */
  r_chain_ptr = first_r_chain_ptr;

  /* Loop over all the currently defined chains */
  while (r_chain_ptr != NULL)
    {
      /* Only process if this is a single-residue ligand chain */
      if (r_chain_ptr->type == R_CHAIN_LIGAND && r_chain_ptr->nresid == 1)
        {
          /* Get the residue */
          r_residue_ptr = r_chain_ptr->first_r_residue_ptr;

          /* If only a single atom, check whether its is a metal */
          if (r_residue_ptr->natoms == 1)
            {
              /* Get the atom */
              r_atom_ptr = r_residue_ptr->first_r_atom_ptr;

              /* Check whether atom is a metal */
              metal = TRUE;
              for (itype = 1; itype < R_NATOM_TYPES && metal == TRUE;
                   itype++)
                {
                  /* If one of the standard atoms, then not a metal */
                  if (!strcmp(r_residue_ptr->first_r_atom_ptr->element,
                              R_ELEMENT[itype]))
                    metal = FALSE;
                }

              /* If atom is a metal, mark residue and chain as metals */
              if (metal == TRUE)
                {
                  /* Set flags to mark residue and chain as a single
                     metal ion*/
                  r_residue_ptr->type = R_TYPE_METAL_ION;
                  r_chain_ptr->type = R_CHAIN_METAL;
                }
            }

          /* Get first chain */
          other_r_chain_ptr = first_r_chain_ptr;

          /* Loop over all the protein and DNA chains to see if this
             might be a residue attachment */
          while (other_r_chain_ptr != NULL)
            {
              /* Only process if this is a protein, ligand or DNA chain
                 with the same chain-id as our current residue */
              if (other_r_chain_ptr != r_chain_ptr &&
                  (other_r_chain_ptr->type == R_CHAIN_PROTEIN ||
                   other_r_chain_ptr->type == R_CHAIN_DNA ||
                   other_r_chain_ptr->type == R_CHAIN_LIGAND) &&
                  other_r_chain_ptr->chain_id == r_residue_ptr->chain)
                {
                  /* Get the first residue */
                  other_r_residue_ptr = other_r_chain_ptr->first_r_residue_ptr;
                  nresid = other_r_chain_ptr->nresid;
                  iresid = 0;
                  found = FALSE;

                  /* Loop through the chain's residues to find a
                     match on the residue number */
                  while (other_r_residue_ptr != NULL && iresid < nresid &&
                         found == FALSE)
                    {
                      /* Check for match on residue number */
                      if (!strcmp(other_r_residue_ptr->res_num,
                                  r_residue_ptr->res_num))
                        {
                          /* Check for a covalent bond between these
                             two residues */
                          bonded
                            = check_for_bond(r_residue_ptr,other_r_residue_ptr);

                          /* If bonded, have a residue attachment, so set
                             the residue and chain types accordingly */
                          if (bonded == TRUE)
                            {
                              /* Set chain type */
                              r_chain_ptr->type = R_CHAIN_ATTACHMENT;

                              /* Store pointers between the two
                                 residue */
                              //r_residue_ptr->attach_r_residue_ptr
                              //  = other_r_residue_ptr;
                              //other_r_residue_ptr->attach_r_residue_ptr
                              //  = r_residue_ptr;
                            }

                          /* Set flag that found a match */
                          found = TRUE;
                        }

                      /* Get this chain's next residue */
                      other_r_residue_ptr
                        = other_r_residue_ptr->next_r_residue_ptr;
                      iresid++;
                    }
                }

              /* Get the next chain */
              other_r_chain_ptr = other_r_chain_ptr->next_r_chain_ptr;
            }
        }

      /* Get the next chain */
      r_chain_ptr = r_chain_ptr->next_r_chain_ptr;
    }
}
/***********************************************************************

add_res_name  -  Add given residue name to the residue list and sequence

***********************************************************************/

void add_res_name(char *residue_sequence,char *residue_list,
                  char *res_name,char *res_num,char chain)
{
  char ch;
  char residue_name[20];

  int first_pos, i, ipos, last_pos, len, numeric;

  /* Initialise variables */
  len = strlen(res_name);
  numeric = FALSE;
  first_pos = -1;
  last_pos = 0;

  /* Replace blanks by underscores in residue name used for residue
     sequence */
  strcpy(residue_name,res_name);
  for (i = 0; i < 3; i++)
    if (residue_name[i] == ' ')
      residue_name[i] = '_';

  /* If residue sequence too long, then abort this entry */
  if (len + strlen(residue_sequence) > R_LINELEN - 1)
    {
      /* See if we already have a ... in the string */
      if (strstr(residue_sequence," ...") == NULL)
        {
          /* Append ellipsis onto end of current residue list */
          strcat(residue_sequence," ...");

          /* Show warning */
          printf("*** Warning. Ligand residue sequence too long! Some "
                 "residues lost.\n");
          printf("***        Residue sequence %s\n",residue_sequence);
        }
    }

  /* Otherwise, add onto name of current residue sequence */
  else
    {
      if (residue_sequence[0] != '\0')
        strcat(residue_sequence,"-");
      strcat(residue_sequence,residue_name);
    }

  /* Locate first and last non-blank positions in the residue name */
  for (ipos = 0; ipos < len; ipos++)
    {
      /* Get first character */
      ch = res_name[ipos];

      /* Check whether non-blank */
      if (first_pos == -1 && ch != ' ')
        first_pos = ipos;
      if (ch != ' ')
        last_pos = ipos;

      /* Check whether numeric */
      if (ch < 'A' || ch > 'Z')
        numeric = TRUE;
    }

  /* Form modified residue name */
  if (first_pos > -1)
    {
      residue_name[0] = '\0';
      if (numeric == TRUE)
        strcat(residue_name,"[");
      strncat(residue_name,res_name+first_pos,last_pos - first_pos + 1);
      if (numeric == TRUE)
        strcat(residue_name,"]");
    }

  /* Re-initialise variables for pre=ocessing the residue number */
  len = strlen(res_num);
  first_pos = -1;
  last_pos = 0;

  /* Locate first and last non-blank positions in the residue number */
  for (ipos = 0; ipos < len; ipos++)
    {
      /* Get first character */
      ch = res_num[ipos];

      /* Check whether non-blank */
      if (first_pos == -1 && ch != ' ')
        first_pos = ipos;
      if (ch != ' ')
        last_pos = ipos;
    }

  /* Add residue number onto residue name */
  len = strlen(residue_name);
  if (first_pos > -1)
    {
      strncat(residue_name,res_num+first_pos,last_pos - first_pos + 1);
      residue_name[len + last_pos - first_pos + 1] = '\0';
    }

  /* If chain not blank, add that, too */
  if (chain != ' ')
    {
      /* Get current length of residue name */
      len = strlen(residue_name);

      /* Add colon and chain identifier */
      residue_name[len] = ':';
      residue_name[len + 1] = chain;
      residue_name[len + 2] = '\0';
    }

  /* If residue list too long, then show warning */
  if (strlen(residue_name) + strlen(residue_list) > R_LINELEN - 7)
    {
      /* See if we already have a ... in the string */
      if (strstr(residue_list," ...") == NULL)
        {
          /* Append ellipsis onto end of current residue list */
          strcat(residue_list," ...");

          /* Show warning */
          printf("*** Warning. Ligand residue list too long! Some  "
                 "residues lost.\n");
          printf("***        Residue list %s + residue name %s\n",
                 residue_list,residue_name);
        }
    }

  /* Otherwise, add onto end of current residue list */
  else
    {
      if (residue_list[0] != '\0')
        strcat(residue_list,", ");
      strcat(residue_list,residue_name);
    }
}
/***********************************************************************

r_format_residue_number  -  Format the ligand's residue range

***********************************************************************/

void r_format_residue_number(char *residue_number,
                           struct r_residue *r_residue_ptr)
{
  char chain_id[4], tmp_residue_number[15];
  char *string_ptr;

  int ipos, len;
  int done;

  /* Initialise variables */
  done = FALSE;
  residue_number[0] = '\0';

  /* Add residue name and number */
  strcat(residue_number,r_residue_ptr->res_name);
  strcat(residue_number," ");
  strcat(residue_number,r_residue_ptr->res_num);

  /* Loop while removing all runs of pairs of consecutive blank spaces */
  while (done == FALSE)
    {
      /* Check for a pair of adjacent spaces */
      string_ptr = strstr(residue_number,"  ");

      /* If found, then replace by a single space */
      if (string_ptr != NULL)
        {
          /* Shift second half of number left one space */
          ipos = string_ptr - residue_number;
          strcpy(tmp_residue_number,residue_number);
          tmp_residue_number[ipos] = '\0';
          strcat(tmp_residue_number,string_ptr+1);
          strcpy(residue_number,tmp_residue_number);
        }

      /* Otherwise we are done */
      else
        done = TRUE;
    }

  /* Remove initial space, if present */
  if (residue_number[0] == ' ')
    strcpy(residue_number,residue_number+1);

  /* Ditto for final space */
  len = strlen(residue_number);
  if (residue_number[len - 1] == ' ')
    residue_number[len - 1] = '\0';

  /* If chain-id not blank, then add */
  if (r_residue_ptr->chain != ' ')
    {
      /* Chain id enclosed in brackets */
      sprintf(chain_id,"(%c)",r_residue_ptr->chain);
      strcat(residue_number,chain_id);
    }
}
/***********************************************************************

form_ligands  -  Create all the ligand records for each of the ligand
                 chains

***********************************************************************/

int form_ligands(struct r_chain *first_r_chain_ptr,
                 struct r_ligand **fst_r_ligand_ptr,int *nlig,int *nmet,
                 int biomol)
{
  char residue_list[R_LINELEN + 1], residue_range[R_LINELEN + 1];
  char residue_sequence[R_LINELEN + 1];
  char start_residue[15], end_residue[15];
  char stored_metal[R_NMETALS][3];

  int colour, i, iresid, metal_type, nligands, nligmet, nmetals, nmetyp;
  int metal_count[R_NMETALS], nresid;
  int metal;

  struct r_atom *r_atom_ptr;
  struct r_chain *r_chain_ptr;
  struct r_residue *r_residue_ptr, *end_r_residue_ptr, *start_r_residue_ptr;
  struct r_ligand *r_ligand_ptr, *first_r_ligand_ptr, *last_r_ligand_ptr;

  /* Initialise variables */
  nligands = 0;
  nligmet = 0;
  nmetals = 0;

  /* Initialise metal colours */
  for (i = 0; i < R_NMETALS; i++)
    {
      if (i < R_NCOMMON_METALS)
        strcpy(stored_metal[i],R_COMMON_METAL[i]);
      else
        stored_metal[i][0] = '\0';
      metal_count[i] = 0;
    }
  nmetyp = R_NCOMMON_METALS;

  /* Initialise ligand pointers */
  r_ligand_ptr = first_r_ligand_ptr = last_r_ligand_ptr = NULL;

  /* Get pointer to the first chain */
  r_chain_ptr = first_r_chain_ptr;

  /* Loop over all the currently defined chains */
  while (r_chain_ptr != NULL)
    {
      /* Only process if this is a ligand or metal chain */
      if ((r_chain_ptr->type == R_CHAIN_LIGAND ||
          r_chain_ptr->type == R_CHAIN_METAL) && r_chain_ptr->nresid > 0)
        {
          /* Determine which type this is */
          if (r_chain_ptr->type == R_CHAIN_LIGAND)
            metal = FALSE;
          else
            metal = TRUE;

          /* Initialise variables */
          start_r_residue_ptr = NULL;
          end_r_residue_ptr = NULL;

          /* Create a ligand record */
          r_ligand_ptr
            = create_r_ligand_record(&first_r_ligand_ptr,
                                     &last_r_ligand_ptr,r_chain_ptr,
                                     nligands,metal);

          /* If this is a biomolecule ligand, assign number to type */
          if (biomol == TRUE && metal == FALSE)
            r_ligand_ptr->type = nligands;
          else if (biomol == TRUE && metal == TRUE)
            r_ligand_ptr->type = nmetals;

          /* Initialise residue list and sequence */
          residue_list[0] = '\0';
          residue_sequence[0] = '\0';

          /* Get the first of this chain's residues */
          r_residue_ptr = r_chain_ptr->first_r_residue_ptr;
          nresid = r_chain_ptr->nresid;
          iresid = 0;
  
          /* If a metal, get its name */
          if (r_ligand_ptr->metal == TRUE && r_residue_ptr != NULL )
            {
              /* Get the atom corresponding to this metal ion */
              r_atom_ptr = r_ligand_ptr->first_r_residue_ptr->first_r_atom_ptr;

              /* Store the metal name */
              strcpy(r_ligand_ptr->metal_name,r_atom_ptr->element);

              /* Get the corresponding colour */
              colour = -1;
              for (metal_type = 0; metal_type < nmetyp && colour == -1;
                   metal_type++)
                {
                  /* If element name matches, store colour */
                  if (!strcmp(stored_metal[metal_type],r_atom_ptr->element))
                    {
                      /* Retrieve the metal's colour */
                      colour = metal_type % R_NMETAL_COLOURS;

                      /* Update count of metal ions in structure */
                      metal_count[metal_type]
                        = metal_count[metal_type] + r_ligand_ptr->count;
                    }
                }

              /* Store the metal's colour */
              if (colour < 0)
                colour = 0;
              r_ligand_ptr->metal_colour = colour;
            }

          /* Loop through the chain's residues to add to current ligand */
          while (r_residue_ptr != NULL && iresid < nresid)
            {
              /* Add current residue name to the ligand residue list */
              add_res_name(residue_sequence,residue_list,
                           r_residue_ptr->res_name,r_residue_ptr->res_num,
                           r_residue_ptr->chain);

              /* If this is the first or last residue of the ligand,
                 store */
              if (iresid == 0)
                start_r_residue_ptr = r_residue_ptr;
              if (iresid == nresid - 1)
                end_r_residue_ptr = r_residue_ptr;

              /* Store ligand pointer */
              r_residue_ptr->r_ligand_ptr = r_ligand_ptr;

              /* Update atom count */
              r_ligand_ptr->natoms
                = r_ligand_ptr->natoms + r_residue_ptr->natoms;

              /* Get the next residue */
              r_residue_ptr = r_residue_ptr->next_r_residue_ptr;
              iresid++;
            }

          /* Increment count of ligand residues read in */
          if (metal == TRUE)
            nmetals++;
          else
            nligands++;
          nligmet++;

          /* Store the ligand's residue list and sequence */
          store_string(&(r_ligand_ptr->residue_list),residue_list);
          store_string(&(r_ligand_ptr->residue_sequence),residue_sequence);

          /* Format the ligand's residue range */
          r_format_residue_number(start_residue,start_r_residue_ptr);
          r_format_residue_number(end_residue,end_r_residue_ptr);

          /* Form the range */
          strcpy(residue_range,start_residue);
          if (nresid > 1)
            {
              strcat(residue_range," to ");
              strcat(residue_range,end_residue);
            }

          /* Store the ligand's residue range */
          store_string(&(r_ligand_ptr->residue_range),residue_range);
        }

      /* Get the next chain */
      r_chain_ptr = r_chain_ptr->next_r_chain_ptr;
    }

  /* Return pointer to the first ligand */
  *fst_r_ligand_ptr = first_r_ligand_ptr;

  /* Return total numbers of ligand and metals */
  *nlig = nligands;
  *nmet = nmetals;
  return(nligmet);
}
/***********************************************************************

identify_unique_ligands  -  Determine which ligands are unique

***********************************************************************/

int identify_unique_ligands(struct r_ligand *first_r_ligand_ptr,
                            int *nmetyp)
{
  char element[3];

  int nequiv, max_equiv, ninter, ninteract, nmetypes, ntypes, nunique;
  int new_type;

  struct r_atom *r_atom_ptr;
  struct r_ligand *r_ligand_ptr;
  struct r_ligand *copy_r_ligand_ptr, *equiv_r_ligand_ptr;
  struct r_ligand *other_r_ligand_ptr;

  /* Initialise variables */
  nmetypes = 0;
  ntypes = 0;
  nunique = 0;

  /* Get pointer to the first ligand record */
  r_ligand_ptr = first_r_ligand_ptr;

  if (r_ligand_ptr != NULL)
    {
      /* The first ligand is, by definition, unique */
      r_ligand_ptr->unique = TRUE;
      if (r_ligand_ptr->metal == TRUE)
        {
          nmetypes++;
          r_ligand_ptr->type = nmetypes;
        }
      else
        {
          ntypes++;
          r_ligand_ptr->type = ntypes;
        }
      nunique++;

      /* Get the next ligand */
      r_ligand_ptr = r_ligand_ptr->next_r_ligand_ptr;
    }

  /* Loop through all the ligands to determine which are unique
     and which are symmetry-related copies */
  while (r_ligand_ptr != NULL)
    {
      /* Initialise variables */
      max_equiv = 0;
      ninteract = 0;
      copy_r_ligand_ptr = NULL;
      equiv_r_ligand_ptr = NULL;
      new_type = TRUE;

      /* If a metal, get its name */
      if (r_ligand_ptr->metal == TRUE)
        {
          /* Get the atom corresponding to this metal ion */
          r_atom_ptr = r_ligand_ptr->first_r_residue_ptr->first_r_atom_ptr;

          /* Store the metal name */
          strcpy(r_ligand_ptr->metal_name,r_atom_ptr->element);
        }

      /* Get pointer to the first ligand */
      other_r_ligand_ptr = first_r_ligand_ptr;

      /* Loop through all the previous ligands to see if this is a
         copy of any one of them */
      while (other_r_ligand_ptr != NULL &&
             other_r_ligand_ptr != r_ligand_ptr)
        {
          /* Only consider this ligand if it is a unique one */
          if (other_r_ligand_ptr->unique == TRUE)
            {
              /* Check whether these two ligands have the same string
                 of residues */
              if (!strcmp(r_ligand_ptr->residue_sequence,
                          other_r_ligand_ptr->residue_sequence))
                {
                  /* Assign same ligand type */
                  r_ligand_ptr->type = other_r_ligand_ptr->type;
                  new_type = FALSE;

                  /* Save ligand that this is a copy of */
                  if (copy_r_ligand_ptr == NULL)
                    copy_r_ligand_ptr = other_r_ligand_ptr;
                }
            }

          /* Get next ligand record */
          other_r_ligand_ptr = other_r_ligand_ptr->next_r_ligand_ptr;
        }

      /* If this is a new ligand type, assign it the next type
         number */
      if (new_type == TRUE)
        {
          /* Increment type count and assign to current ligand */
          if (r_ligand_ptr->metal == TRUE)
            {
              nmetypes++;
              r_ligand_ptr->type = nmetypes;
            }
          else
            {
              ntypes++;
              r_ligand_ptr->type = ntypes;
            }
        }

      /* If have enough matching residues, then current ligand is
         a duplicate */
      if (max_equiv > ninteract / 2)
        {
          /* Mark as non-unique */
          r_ligand_ptr->unique = FALSE;

          /* Save ligand of which it is a duplicate */
          equiv_r_ligand_ptr->nduplicates++;

          /* Store this as the duplicate number of this ligand */
          r_ligand_ptr->unique_count = equiv_r_ligand_ptr->unique_count;
          r_ligand_ptr->duplicate_count = equiv_r_ligand_ptr->nduplicates;
          r_ligand_ptr->duplicate_r_ligand_ptr = equiv_r_ligand_ptr;
        }

      /* Otherwise, deem it to be unique */
      else
        {
          /* Mark as unique */
          r_ligand_ptr->unique = TRUE;

          /* Get its copy number */
          if (copy_r_ligand_ptr != NULL)
            {
              /* Increment count of copies of original ligand type */
              copy_r_ligand_ptr->ncopy++;

              /* Store this as the copy number of this ligand */
              r_ligand_ptr->unique_count = copy_r_ligand_ptr->ncopy;
            }
        }

      /* Increment the total count of this ligand type */
      if (copy_r_ligand_ptr != NULL)
        copy_r_ligand_ptr->count++;

      /* Get next ligand record */
      r_ligand_ptr = r_ligand_ptr->next_r_ligand_ptr;
    }

  /* Return the numbers of ligand and metal types */
  *nmetyp = nmetypes;
  return(ntypes);
}
/***********************************************************************

r_shell_sort  -  Shell-sort which sorts the ligand types into ascending
                 order ot type, copy and duplicate

***********************************************************************/

void r_shell_sort(int *index,struct r_ligand **r_ligand_index_ptr,int nindex)
{
  int endloop, i, indi, indl, iswap, j, k, l, lognb2, m, n, nn;
  float aln2i = (float) 1.4426950, tiny = (float) 1.e-5;
  float calc;
  double dble, dble1;

  struct r_ligand *r_ligand1_ptr, *r_ligand2_ptr;

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

              /* Get the two ligand records */
              r_ligand1_ptr = r_ligand_index_ptr[index[indi]];
              r_ligand2_ptr = r_ligand_index_ptr[index[indl]];

              /* If in wrong order, the swap sort-pointers */
              if (r_ligand1_ptr->sort_score > r_ligand2_ptr->sort_score)
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

r_sort_ligands  -  Sort the ligands by their score and delete any of
                   the lower-scoring ones that overlap the higher
                   scoring ones significantly

***********************************************************************/

int r_sort_ligands(struct r_ligand **fst_r_ligand_ptr,int nligands)
{
  int iligand, iorder, ipos;
  int *index;

  struct r_ligand *first_r_ligand_ptr, *last_r_ligand_ptr;
  struct r_ligand *r_ligand_ptr;
  struct r_ligand **r_ligand_index_ptr;

  /* Initialise pointers */
  first_r_ligand_ptr = *fst_r_ligand_ptr;
  last_r_ligand_ptr = NULL;

  /* Grab memory for temporary sort arrays */
  index = (int *) malloc(nligands * sizeof(int));
  r_ligand_index_ptr
    = (struct r_ligand **) malloc(nligands * sizeof(struct r_ligand *));

  /* Get pointer to first ligand record */
  r_ligand_ptr = first_r_ligand_ptr;
  ipos = 0;

  /* Loop through ligand records to initialise the sort index */
  while (r_ligand_ptr != NULL)
    {
      /* Initialise sort index */
      index[ipos] = ipos;
      r_ligand_index_ptr[ipos] = r_ligand_ptr;

      /* Calculate the score on which to sort, involving the type,
         instance and duplicate copy */
      r_ligand_ptr->sort_score
        = 1000 * 1000 * r_ligand_ptr->type
        + 1000 * r_ligand_ptr->unique_count
        + r_ligand_ptr->duplicate_count;

      /* Want ligands to come ahead of metals */
      if (r_ligand_ptr->metal == TRUE)
        r_ligand_ptr->sort_score = r_ligand_ptr->sort_score + 10000000;

      /* Increment ligand counter */
      ipos++;

      /* Get next ligand record */
      r_ligand_ptr = r_ligand_ptr->next_r_ligand_ptr;
    }

  /* Initialise ligand counter */
  iligand = 0;
  nligands = ipos;

  /* Sort the ligands */
  r_shell_sort(index,r_ligand_index_ptr,nligands);

  /* Having sorted the ligands, rearrange order of the ligand-record
     pointers in descending order of ligand */
  for (iorder = 0; iorder < nligands; iorder++)
    {
      /* Get the array position of this ligand record */
      ipos = index[iorder];

      /* Get the current ligand record */
      r_ligand_ptr = r_ligand_index_ptr[ipos];

      /* If this is the first, then make the first record */
      if (iorder == 0)
        first_r_ligand_ptr = r_ligand_ptr;

      /* Otherwise, get the previous ligand record to point to this one */
      else
        last_r_ligand_ptr->next_r_ligand_ptr = r_ligand_ptr;

      /* Blank out current record's next ligand pointer */
      r_ligand_ptr->next_r_ligand_ptr = NULL;

      /* Increment ligand counter and store */
      iligand++;
      r_ligand_ptr->number = iligand;

      /* Store current record as last encountered */
      last_r_ligand_ptr = r_ligand_ptr;
    }

  /* Free up memory used for temporary sort arrays */
  free(index);
  free(r_ligand_index_ptr);

  /* Return current pointer */
  *fst_r_ligand_ptr = first_r_ligand_ptr;

  /* Return number of sorted ligands */
  return(nligands);
}
/***********************************************************************

assign_rasdefs  -  Assign the given rasdef to all the residues in this
                   molecule

***********************************************************************/

void assign_rasdefs(struct r_residue *first_r_residue_ptr,int nresid,
                    struct r_rasdef *r_rasdef_ptr)
{
  int iresid;

  struct r_residue *r_residue_ptr;

  /* Initialise variables */
  iresid = 0;

  /* Get the first residue */
  r_residue_ptr = first_r_residue_ptr;

  /* Loop through all the residues to assign to the given rasdef record */
  while (r_residue_ptr != NULL && iresid < nresid)
    {
      /* Save the rasdef pointer */
      r_residue_ptr->r_rasdef_ptr = r_rasdef_ptr;

      /* Get the next r_residue */
      r_residue_ptr = r_residue_ptr->next_r_residue_ptr;
      iresid++;
    }
}
/***********************************************************************

assign_chains  -  Create the rasdef records

***********************************************************************/

struct r_rasdef *assign_chains(struct r_chain *first_r_chain_ptr,
                               struct r_ligand *first_r_ligand_ptr)
{
  char chain, ctype, input_line[R_LINELEN], last_sequence[10];
  char metal_name[4], res_name[4], res_num[6];
  char colour[R_COLNAM_LEN];
  char ligand_from[7], ligand_to[7];

  static char ctypes[3] = { 'C', 'P', 'N' };

  int count, loop, metal_type, nchain, ndef, nligands, type;
  int col_no;

  static int types[3] = { R_CALPHA_ONLY, R_PROTEIN, R_NUCLEIC_ACID };

  struct r_chain *r_chain_ptr;
  struct r_ligand *r_ligand_ptr;
  struct r_rasdef *r_rasdef_ptr, *first_r_rasdef_ptr, *last_r_rasdef_ptr;

  /* Initialise variables */
  count = nchain = ndef = nligands = 0;
  ligand_from[0] = '\0';
  ligand_to[0] = '\0';
  chain = ' ';
  res_name[0] = '\0';
  res_num[0] = '\0';
  colour[0] = '\0';

  /* Initialise pointers */
  r_rasdef_ptr = first_r_rasdef_ptr = last_r_rasdef_ptr = NULL;

  /* Loop over protein and DNA chains */
  for (loop = 0; loop < 3; loop++)
    {
      /* Get chain type marker */
      ctype = ctypes[loop];
      type = types[loop];
      nchain++;

      /* Get pointer to the first chain */
      r_chain_ptr = first_r_chain_ptr;

      /* Loop over all the chains */
      while (r_chain_ptr != NULL)
        {
          /* If appropriate chain type, then write out details */
          if ((loop == 0 && r_chain_ptr->type == R_CHAIN_CALPHA_ONLY) ||
              (loop == 1 && r_chain_ptr->type == R_CHAIN_PROTEIN) ||
              (loop == 2 && r_chain_ptr->type == R_CHAIN_DNA))
            {
              /* Save this chain type and chain id */
              sprintf(input_line,"%c%c",ctype,r_chain_ptr->chain_id);

              /* Create a r_rasdef entry */
              r_rasdef_ptr
                = create_r_rasdef_entry(&first_r_rasdef_ptr,
                                        &last_r_rasdef_ptr,type,nchain,
                                        res_name,res_num,
                                        r_chain_ptr->chain_id,colour,
                                        input_line,ligand_from,ligand_to,
                                        count);
              ndef++;

              /* Assign all this chain's residues to this rasdef record */
              assign_rasdefs(r_chain_ptr->first_r_residue_ptr,
                             r_chain_ptr->nresid,r_rasdef_ptr);
            }

          /* Get the next chain */
          r_chain_ptr = r_chain_ptr->next_r_chain_ptr;
        }
    }

  /* Get pointer to the first ligand record */
  r_ligand_ptr = first_r_ligand_ptr;

  /* Loop over the ligands */
  while (r_ligand_ptr != NULL)
    {
      /* Write out the details of ligand */
      if (r_ligand_ptr->metal == FALSE && r_ligand_ptr->natoms > 0)
        {
          /* Increment ligand counter */
          nligands++;

          /* Save the residue list */
          sprintf(input_line,"L%3.3d %s",nligands,
                  r_ligand_ptr->residue_list);

          /* Get the chain identifier */
          chain = r_ligand_ptr->first_r_residue_ptr->chain;
          type = R_LIGAND;

          /* Create a r_rasdef entry */
          r_rasdef_ptr
            = create_r_rasdef_entry(&first_r_rasdef_ptr,
                                    &last_r_rasdef_ptr,type,nchain,
                                    res_name,res_num,chain,colour,
                                    input_line,ligand_from,ligand_to,
                                    count);
          ndef++;

          /* Assign all this chain's residues to this rasdef record */
          assign_rasdefs(r_ligand_ptr->first_r_residue_ptr,
                         r_ligand_ptr->nresid,r_rasdef_ptr);

          /* Save the residue range */
          sprintf(input_line,"l%3.3d %s",nligands,
                  r_ligand_ptr->residue_range);
          type = R_LIGAND_RANGE;

          /* Create a r_rasdef entry */
          r_rasdef_ptr
            = create_r_rasdef_entry(&first_r_rasdef_ptr,
                                    &last_r_rasdef_ptr,type,nchain,
                                    res_name,res_num,chain,colour,
                                    input_line,ligand_from,ligand_to,
                                    count);
          ndef++;
        }

      /* Get the next ligand */
      r_ligand_ptr = r_ligand_ptr->next_r_ligand_ptr;
    }

  /* Get pointer to the first ligand record again */
  r_ligand_ptr = first_r_ligand_ptr;
  last_sequence[0] = '\0';
  metal_type = 0;
  type = R_METAL;

  /* Loop over the ligand records to write out the metals */
  while (r_ligand_ptr != NULL)
    {
      /* Write out the details of the metal */
      if (r_ligand_ptr->metal == TRUE &&
          strcmp(r_ligand_ptr->residue_sequence,last_sequence))
        {
          /* Increment metal counter */
          metal_type++;

          /* Get the metal colour */
          col_no = r_ligand_ptr->metal_colour;

          /* Assign all this chain's residues to this rasdef record */
          assign_rasdefs(r_ligand_ptr->first_r_residue_ptr,
                         r_ligand_ptr->nresid,r_rasdef_ptr);

          /* Save the residue list */
          sprintf(input_line,"M%s %s",r_ligand_ptr->metal_name,
                  R_METAL_COLOUR[col_no]);

          /* Get the chain identifier */
          chain = r_ligand_ptr->first_r_residue_ptr->chain;

          /* Create a r_rasdef entry */
          r_rasdef_ptr
            = create_r_rasdef_entry(&first_r_rasdef_ptr,
                                    &last_r_rasdef_ptr,type,nchain,
                                    res_name,res_num,chain,colour,
                                    input_line,ligand_from,ligand_to,
                                    count);
          ndef++;

          /* Store metal name */
          strcpy(last_sequence,r_ligand_ptr->residue_sequence);
        }

      /* Get the next ligand */
      r_ligand_ptr = r_ligand_ptr->next_r_ligand_ptr;
    }

  /* Return the first rasdef pointer */
  return(first_r_rasdef_ptr);
}
/***********************************************************************

r_res_assignments  -  Interpret the data read in from the PDB file to
                      identify protein chains, DNA, ligands, etc

***********************************************************************/

int r_res_assignments(struct r_pdb *r_pdb_ptr)
{
  int nchains, nligmet, nligands, nmetals, nresid, nshift, ntypes;
  int nmetypes;

  struct r_chain *r_chain_ptr, *last_r_chain_ptr;
  struct r_ligand *first_r_ligand_ptr;
  struct r_rasdef *r_rasdef_ptr, *first_r_rasdef_ptr, *last_r_rasdef_ptr;

  /* Initialise pointers */
  first_r_rasdef_ptr = NULL;
  last_r_chain_ptr = NULL;

  /* Get the first chain */
  r_chain_ptr = r_pdb_ptr->first_r_chain_ptr;

  /* Loop over all the currently defined chains to find pointer to
     last one */
  while (r_chain_ptr != NULL)
    {
      /* Store as the last chain */
      last_r_chain_ptr = r_chain_ptr;

      /* Get the next chain */
      r_chain_ptr = r_chain_ptr->next_r_chain_ptr;
    }

  /* Classify the residues according to type */
  classify_residues(r_pdb_ptr,&last_r_chain_ptr);

  /* Shift all chains containing single metal ions to the end of
     the chains linked list */
  nshift
    = shift_metals(&(r_pdb_ptr->first_r_chain_ptr),&last_r_chain_ptr);

  /* Adjust the residues so that their order corresponds to the current
     chain order */
  if (nshift > 0)
    nresid = adjust_residues(&(r_pdb_ptr->first_r_residue_ptr),
                             r_pdb_ptr->first_r_chain_ptr);

  /* Initial classification of the "chains" */
  nchains = classify_chains(r_pdb_ptr->first_r_chain_ptr);

  /* If have more than 1 chain, check for consecutive chains with the
     same chain-ID which can be combine together into a single chain
     (currently split into separate chains by spurious TER records in
     the PDB file identifying breaks in the chain) */
  if (nchains > 1)
    join_chains(r_pdb_ptr->first_r_chain_ptr);

  /* Loop through all the chains to reclassify all very short chains
     as ligands */
  reclassify_chains(r_pdb_ptr->first_r_chain_ptr);

  /* Split the ligand chains into physically distinct molecules */
  split_ligands(r_pdb_ptr->first_r_chain_ptr);

  /* Identify which single-residue ligands are metals and which are
     actually residue attachments */
  single_residue_ligands(r_pdb_ptr->first_r_chain_ptr);

  /* Create all the ligand records for each of the ligand chains */
  nligmet
    = form_ligands(r_pdb_ptr->first_r_chain_ptr,&first_r_ligand_ptr,
                   &nligands,&nmetals,FALSE);

  /* If have any ligands, process them to determine which are unique */
  if (nligmet > 0)
    ntypes = identify_unique_ligands(first_r_ligand_ptr,&nmetypes);

  /* Sort ligands by type, uniqueness and duplicate number */
  if (nligmet > 1)
    r_sort_ligands(&first_r_ligand_ptr,nligmet);

  /* Save all the molecule definitions */
  first_r_rasdef_ptr
    = assign_chains(r_pdb_ptr->first_r_chain_ptr,first_r_ligand_ptr);

  /* Save the first rasdef pointer */
  r_pdb_ptr->first_r_rasdef_ptr = first_r_rasdef_ptr;

  /* Return the number of ligands and metals */
  return(nligmet);
}
/***********************************************************************

r_get_pdb_file  -  Open the given PDB file and read in the data from it

***********************************************************************/

struct r_pdb *r_get_pdb_file(char *file_name)
{
  char pdb_code[5];
  char wanted_chain;

  int data_rqrd, gzipped;

  struct r_pdb *r_pdb_ptr;
  struct r_rasdef *first_r_rasdef_ptr;

#ifdef TEST
  FILE *file_ptr;
#else
  gzFile file_ptr;
#endif

  /* Initialise variables */
  data_rqrd = ALL_DATA;
  strcpy(pdb_code,"none");
  r_pdb_ptr = NULL;
  first_r_rasdef_ptr = NULL;
  wanted_chain = '\0';

  /* Open the PDB file */
  file_ptr = r_open_pdb_file(file_name,&gzipped);
  if (file_ptr != NULL)
    {
      /* Read the PDB file */
      r_pdb_ptr = r_read_pdb_file(file_ptr,gzipped,pdb_code,wanted_chain,
                                  data_rqrd,first_r_rasdef_ptr);

      /* Identify the different components of the PDB file */
      r_res_assignments(r_pdb_ptr);

      /* Assign to separate molecules and calculate coord limits */
      r_pdb_ptr->nmolecs = r_assign_to_molecules(r_pdb_ptr);

      /* Calculate all the residue and molecule coordinate limits */
      r_calc_limits(r_pdb_ptr);

      /* Close the PDB file */
#ifdef TEST
      fclose(file_ptr);
#else
      gzclose(file_ptr);
#endif
    }

  /* Return the PDB record */
  return(r_pdb_ptr);
}
/***********************************************************************

r_free_pdb_data  -  Free up memory taken by PDB data

***********************************************************************/

void r_free_pdb_data(struct r_pdb *r_pdb_ptr)
{
  int iatom, natoms;

  struct r_atom *r_atom_ptr, *next_r_atom_ptr;
  struct r_info *r_info_ptr, *next_r_info_ptr;
  struct r_residue *r_residue_ptr, *next_r_residue_ptr;

  /* Get pointer to the first residue read in from the PDB file */
  r_residue_ptr = r_pdb_ptr->first_r_residue_ptr;

  /* Loop through all the residues to free up memory taken by them and
     their atoms */
  while (r_residue_ptr != NULL)
    {
      /* Get the next residue */
      next_r_residue_ptr = r_residue_ptr->next_r_residue_ptr;

      /* Get the first of this residue's atoms */
      r_atom_ptr = r_residue_ptr->first_r_atom_ptr;
      natoms = r_residue_ptr->natoms;
      iatom = 0;

      /* Loop over this residue's atoms to free up memory */
      while (r_atom_ptr != NULL && iatom < natoms)
        {
          /* Get the next atom record */
          next_r_atom_ptr = r_atom_ptr->next_r_atom_ptr;

          /* Free up the memory taken up by this atom */
          free(r_atom_ptr);

          /* Get pointer to next atom record and increment
             atom counter */
          r_atom_ptr = next_r_atom_ptr;
          iatom++;
        }

      /* Free up this residue record */
      free(r_residue_ptr);

      /* Get the next r_residue */
      r_residue_ptr = next_r_residue_ptr;
    }

  /* Get pointer to the first info record */
  r_info_ptr = r_pdb_ptr->first_r_info_ptr;

  /* Loop through all the info records to free up memory taken by them */
  while (r_info_ptr != NULL)
    {
      /* Get the next info record */
      next_r_info_ptr = r_info_ptr->next_r_info_ptr;

      /* Free up the character fields */
      if (r_info_ptr->name != &null)
        free(r_info_ptr->name);
      if (r_info_ptr->value != &null)
        free(r_info_ptr->value);

     /* Free up this info record */
      free(r_info_ptr);

      /* Get the next r_info record */
      r_info_ptr = next_r_info_ptr;
    }

  /* Free up the character fields */
  if (r_pdb_ptr->pdb_code != &null)
    free(r_pdb_ptr->pdb_code);
  if (r_pdb_ptr->protein_name != &null)
    free(r_pdb_ptr->protein_name);

  /* Free the PDB record */
  free(r_pdb_ptr);
}
/***********************************************************************

r_free_all_data  -  Free up memory taken by all data, including rasdef
                    rescords and PDB records

***********************************************************************/

void r_free_all_data(struct r_rasdef *first_r_rasdef_ptr,
                     struct r_pdb *r_pdb_ptr)
{
  /* @@@ To be done */


  /* Free up the PDB data */
  r_free_pdb_data(r_pdb_ptr);
}
/***********************************************************************

r_write_pdb_coords  -  Write the atom coords to the pymol.pdb file

***********************************************************************/

int r_write_pdb_coords(FILE *file_ptr,
                       struct r_residue *first_r_residue_ptr,int nresid,
                       int atom_number,int omit_hydrogens)
{
  char atom_type[7], number_string[20], res_num[6];

  int i, iatom, iresid, natoms, type;
  int wanted;

  float bvalue, coords[3], occupancy;

  struct r_atom *r_atom_ptr;
  struct r_residue *r_residue_ptr;

  /* Initialise variables */
  strcpy(atom_type,"ATOM  ");
  bvalue = 10.0;
  iresid = 0;
  occupancy = 1.0;

  /* Get pointer to the first residue record */
  r_residue_ptr = first_r_residue_ptr;

  /* Loop through this chain's residues, to write out */
  while (r_residue_ptr != NULL)
    {
      /* Get the first of this residue's atoms */
      r_atom_ptr = r_residue_ptr->first_r_atom_ptr;
      natoms = r_residue_ptr->natoms;
      iatom = 0;

      /* Get residue number */
      strcpy(res_num,r_residue_ptr->res_num);

      /* Determine whether to write out as an ATOM or HETATM record */
      if (r_residue_ptr->r_rasdef_ptr != NULL &&
          (r_residue_ptr->r_rasdef_ptr->type == R_LIGAND ||
           r_residue_ptr->r_rasdef_ptr->type == R_METAL ||
           r_residue_ptr->r_rasdef_ptr->type == R_WATER))
        strcpy(atom_type,"HETATM");
      else
        strcpy(atom_type,"ATOM  ");

      /* Loop over this residue's atoms to write out */
      while (r_atom_ptr != NULL && iatom < natoms)
        {
          /* If atom wanted, then write out */
          if (omit_hydrogens == FALSE ||
              strcmp(r_atom_ptr->element," H"))
            {
              /* Increment atom count */
              atom_number++;
              if (atom_number > 99999)
                atom_number = 1;

              /* Get this atom's coords */
              for (i = 0; i < 3; i++)
                coords[i] = r_atom_ptr->coords[i];

              /* Write out the details */
              fprintf(file_ptr,"%s%5d %s %s %c%s   %8.3f%8.3f%8.3f",
                      atom_type,atom_number,r_atom_ptr->atom_name,
                      r_residue_ptr->res_name,
                      r_residue_ptr->chain,res_num,coords[0],
                      coords[1],coords[2]);
              fprintf(file_ptr,"%6.2f%6.2f\n",r_atom_ptr->occupancy,
                      r_atom_ptr->bvalue);
            }

          /* Get pointer to next atom record and increment atom
             counter */
          r_atom_ptr = r_atom_ptr->next_r_atom_ptr;
          iatom++;
        }

      /* Get the next residue */
      r_residue_ptr = r_residue_ptr->next_r_residue_ptr;
      iresid++;

      /* If we have written out the required number of residues, stop
         here */
      if (iresid == nresid)
        r_residue_ptr = NULL;
    }

  /* Write out a TER record */
  fprintf(file_ptr,"TER\n");

  /* Return the current atom number */
  return(atom_number);
}
/***********************************************************************

r_list_rasdef_entries  -  List all the rasdef entries (for debugging)

***********************************************************************/

int r_list_rasdef_entries(struct r_rasdef *first_r_rasdef_ptr)
{
  int nout;

  struct r_rasdef *r_rasdef_ptr;

  /* Initialise variables */
  nout = 0;
  printf("   Listing rasdef entries ...\n");

  /* Get pointer to the first RasMol definition record */
  r_rasdef_ptr = first_r_rasdef_ptr;

  /* Loop through all the definition records */
  while (r_rasdef_ptr != NULL)
    {
      /* Increment count */
      nout++;

      /* Print this entry */
      if (r_rasdef_ptr->type != R_UNKNOWN)
        printf("   %4d. %2d %s %s %s %c %d:%d [%s] [%s] [%s] %d [%s]\n",
               nout,
               r_rasdef_ptr->type,
               R_RESTYPE_NAME[r_rasdef_ptr->type],
               r_rasdef_ptr->res_name,
               r_rasdef_ptr->res_num,
               r_rasdef_ptr->chain,
               r_rasdef_ptr->number,
               r_rasdef_ptr->copy,
               r_rasdef_ptr->colour,
               r_rasdef_ptr->ligand_name,
               r_rasdef_ptr->ligand_range,
               r_rasdef_ptr->nresid,
               r_rasdef_ptr->ligplot_id);

      /* Get next definition record */
      r_rasdef_ptr = r_rasdef_ptr->next_r_rasdef_ptr;
    }

   return(nout);
}
/***********************************************************************

                             M   A   I   N

***********************************************************************/

int r_main_test(int argc,char *argv[])
//int main(int argc,char *argv[])
{
   exit(0);
}

