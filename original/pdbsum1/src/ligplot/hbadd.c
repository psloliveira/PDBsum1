/**************************************************************************

 hbadd.c - Program to prepare a set of hbplus.rc entries for a given PDB
           file. An hbplus.rc entry is made for every HET group in the
           PDB file that is correctly located in the HET group dictionary,
           het_dictionary.txt, obtained from the PDB (from the PDB
           documentation page at http://www.pdb.bnl.gov/doc_help.html).

           The hbplus.rc entry defines the connectivities of the atoms
           in the HET group and each atom's numbers of available donors
           and acceptors. It is used by program HBPLUS for calculating
           H-bonds.

*************************************************************************

 Version:-      v.1.0  -  30 January 1997
 Written by:-   Nicholas M Luscombe(1) and Roman A Laskowski(2)
 
            (1) Department of Biochemistry and Molecular Biology
                University College London
                Gower Street, London WC1E 6BT, UK.
 
            (2) Department of Crystallography
                Birkbeck College
                Malet Street,  London WC1E 7HX, UK.

 Requirements:- An input PDB file for which the hbplus.rc entries are
                required. The HET Group Dictionary, het_dictionary.txt,
                obtained from the PDB from:-

                http://www.pdb.bnl.gov/doc_help.html



 ------------------------------------------------------------------------

 Program description
 -------------------

 Program reads in the given PDB file and identifies all the HET groups
 involved. Computes connectivities using CONECT records and distances
 between atoms.

 Locates each HET group in the HET Group Dictionary. Connectivities are
 recorded. Atoms are matched by name and then by connectivity. Where the
 HET group in the PDB file successfully matches the dictionary definition
 on both these criteria, the relevant bond angles are computed where
 necessary.

 Rules for H-bond formation are applied to all O and N atoms according to 
 the following rules (John Mitchell, personal communication):-
 
 H-bond donors: 
    Any O, S or N potentially donates the number of Hs bound to it.

 H-bond acceptors:
    O/S: sp3 - can accept up to 2
         sp2 - can accept up to 2

    N:   planar   - no acceptance
         sp3      - can accept 1
         sp2      - can accept 1
         aromatic - can accept (3 - no. of bonds)
         amide    - no acceptance

 ------------------------------------------------------------------------

 Original version was part of v.3.0 of LIGPLOT. Amendments from v.4.3
 onwards have been labelled by v.m.n--> and <--v.m.n where m.n is the
 version number corresponding to the change

 v.3.1          Amendment to write out the hbplus.rc file for HET groups
                even if not in the Het Group Dictionary, using the
                connectivity information in the PDB file.
                                                        12 Apr 1997 (RAL)
                Bug-fix to write out 2-character residue-names in the
                format that HBPLUS requires.
                                                        15 Apr 1997 (RAL)
 v.3.1.2        Amendment to stop reading PDB file ATOM data when have
                hit an ENDMDL record.
                Amendment for metals to fool HBPLUS in recognizing bonds
                to metal ions.
                Bug-fix in routine validate_pdb_hets from v.3.1.
                                                        23 Apr 1997 (RAL)
 v.3.2          Addition of graph-matching routines to identify atom
                correspondences between the PDB het group and the entry
                in the Het Group Dictionary when the atom names don't
                match. Functions check_atoms_present and
                check_connectivities replaced by function graph_match and
                its subfunctions.
                Generation of fake entries in cases when no match with
                the Het Group Dictionary possible. Numbers of acceptors
                and donors set to maximum for given atom type and bond
                connectivity in these cases.
                                                         1 May 1997 (RAL)
                Amendment for atom names containing double-quotes (eg
                in 1gac). Double-quotes written out as @ so that HBPLUS
                doesn't fail on reading string in. Amendment made to
                HBPLUS to convert @'s back to double-quotes.
                                                        25 May 1997 (RAL)
                Amendment to open hbplus.rc file only the first time it
                is needed, so that an empty file is not created if no data
                needs to be written out (or this causes problems for
                HBPLUS)
                                                         8 Jun 1997 (RAL)
 v.3.3          Amendment to include ATOM-record residues as well as 
                HETATM ones (specifically for dealing with non-standard
                amino acids and nucleic acids).
                                                         6 Oct 1997 (RAL)

 v.3.4          Addition of hbadd.sum file summarising matches between
                residues in the PDB file and their definitions in the
                Het Group Dictionary. (Used by PDBsum programs to pick
                the most representative structure for any het group type).
                                                        15 Sep 1998 (RAL)

 v.4.3          Amendments to read CIF-format Het Group Dictionary.
                                                        17 Nov 2000 (RAL)

 v.4.4          Rewrite of clique-detection routine as was failing to
                search all possible subgraphs and hence missing the
                correct mapping for some molecules (eg CSM in 1chm or
                any of the FADs). Spotted by John Mitchell.
                                                         5 Dec 2000 (RAL)

                Addition of routines to write out coordinates of each
                Het Group encountered for use in generating pages in

                PDBsum.
                                                         6 Feb 2001 (RAL)

 v.4.4.1        Possible bug-fix to interpretation of old-style Het
                Group Dictionary.
                                                        25 May 2001 (RAL)

 v.4.4.2        Modification for new-format components.cif with new field
                in atom definitions.
                                                         9 Oct 2007 (RAL)

 v.4.4.3        Bug-fix for bromines (and other metals) that are covalently
                bonded.
                                                        21 Oct 2008 (RAL)

                Check for PDB line length when extracting element from
                end of line.
                                                         8 Jan 2009 (RAL)

v.5.0           Addition of -wkdir command-line argument to provide
                working directory
                                                         6 Feb 2009 (RAL)

v.5.0.1         Change to MIN_MAPPED parameter for mapping enzyme reaction
                molecules.
                                                        29 Nov 2010 (RAL)

v.5.1           Change to allow a Het Group Dictionary in SDF format.
                                                      17 Jan 2011 (RAL)

v.5.1.1         Bug-fix to prevent program exiting too early if it fails
                to find any of the Het groups in the dictionary.
                                                       2 May 2011 (RAL)
v.5.1.2         Modifications to strcpy(field,field + n) commands which
                fail under some linux operating systems.
                                                        RAL 24 Jan 2012
v.5.3.1         Bug-fix for faults when failing to find Het Group 
                dictionary.
                                                        RAL  8 Oct 2012
v.5.3.3         Bug-fix for cases where dictionary empty causing program
                to crash.
                                                        RAL 13 May 2013
v.5.3.7         Bug-fix to identify flourine as a metal.
                                                        RAL  4 Jan 2016
v.5.4.3         Change for new format of the Het Group Dictionary
                                                        RAL  8 Nov 2020
v.5.4.5         Modification for metals not being written out to the
                hbplus.rc file.
                                                        RAL 24 Mar 2022

 ------------------------------------------------------------------------

 Datafiles
 ---------

 Actual name          File pointer    Description
 -----------          ------------    -----------
 <filename>.pdb       fil_pdb         Input PDB file holding coords of
                                      protein and for which hbplus.rc 
                                      entries are required.
 het_dictionary.txt   fil_dic         The Het Group Dictionary.
 hbplus.rc            fil_out         File containing the hbplus.rc
                                      entries.
 hbadd.sum            fil_sum         Summary file giving matches between
                                      residue in PDB file and its
                                      definition in the Het Group
                                      Dictionary.
 hbadd.map            fil_map         Output file showing the mapping
                                      between the PDB atoms and the
                                      dictionary atoms
 hetxx.pdb               -            Output PDB file(s) for each Het
                                      Group encountered for use in PDBsum
 
 ------------------------------------------------------------------------
 
 Function calling-tree
 ---------------------
 
 main
    -> get_command_arguments
    -> read_pdbfile
          -> to_lower
          -> get_dfn_data
                -> string_truncate
                -> create_rasdef_entry
                      -> store_string
                -> get_ligand_residues
                      -> get_token
                      -> format_res_number
                      -> create_rasdef_entry
                            -> store_string
          -> assign_to_object
          -> check_residue_name
          -> create_residue_entry
          -> is_duplicate
          -> get_element
                -> is_hydrogen_atom
                      -> check_for_metal
                -> check_for_metal
          -> check_element
          -> store_connections
    -> match_up_ligands
          -> get_ligand
                -> extract_range
                      -> get_token
    -> read_sdf_dictionary
          -> string_truncate
          -> get_tokens
          -> to_upper
          -> get_attype
          -> find_atom
          -> create_atom_record
          -> create_bond_record
          -> check_hetname
          -> add_to_dictionary
                -> create_dictionary_entry
          -> create_dic_data_entry
          -> store_sdf_atoms
                -> create_dic_data_entry
          -> store_sdf_bonds
                -> create_dic_data_entry
          -> free_atom_memory
          -> free_bond_memory
    -> read_dictionary
          -> string_truncate
          -> add_to_dictionary
                -> create_dictionary_entry
          -> create_dic_data_entry
          -> interpret_conect_line
                -> get_element
                      -> is_hydrogen_atom
                            -> check_for_metal
                      -> check_for_metal
          -> interpret_cif_line
                -> tidy_line
                -> cif_bond_line
                      -> get_atom_name
                -> cif_atom_line
                      -> get_atom_name
          -> store_bond
                -> create_dic_data_entry
                -> is_hydrogen_atom
    -> generate_duplicate_entries
          -> create_dic_data_entry
          -> store_bond
                -> create_dic_data_entry
    -> dic_connectivities
    -> open_sumfile
    -> open_mapfile
    -> validate_pdb_hets
          -> graph_match
                -> calc_covalent_connectivities
                -> create_equiv_atoms
                      -> name_match_score
                -> calc_equiv_bonds
                      -> create_equiv_bonds
                -> check_for_metal
                -> initialise_atoms
                -> initialise_dic_data
                -> get_best_possible_score
                -> check_perfect_match
                      -> initialise_equiv_atoms
                -> search_paths
                      -> initialise_equiv_bonds
                      -> initialise_equiv_atoms
                      -> create_tree_entry
                      -> process_tree
                            -> check_compatibility
                            -> create_tree_entry
                            -> delete_branch
                            -> reset_path
                                  -> initialise_equiv_atoms
                            -> traverse_path
                                  -> check_branches
                            -> score_path
                                  -> get_score
                                  -> get_bond_score
                            -> free_paths
                -> process_best_path
                      -> get_bond_score
                      -> reset_path
                            -> initialise_equiv_atoms
          -> write_sumfile
          -> print_missing_hets
          -> generate_fake_entry
                -> generate_connectivities
                      -> check_stored_conecs
                -> create_dictionary_entry
                -> create_dic_data_entry
    -> generate_hbplusrc
          -> open_hbplusrc
          -> writename_fn
                -> convert_quotes
          -> writecon_fn
          -> writehbdon_fn
                -> check_for_metal
          -> writehbacc_fn
                -> check_atoms_present
    -> best_dic_match
    -> write_bonds
    -> write_hetpdb_files


 ----------------------------------------------------------------------

 Compiling
 ---------

cc -o hbadd hbadd.c -lm

Optimization flags:-  cc -O2 -funroll-loops -o hbadd hbadd.c -lm

***********************************************************************/


/* H E A D E R   F I L E */

#include "hbadd.h"



/***********************************************************************

get_command_arguments  -  Get the name of the input PDB file and Het Group
                          Dictionary from the command-line arguments

***********************************************************************/

void get_command_arguments(char *strptrs[],int num_args,
                           char *pdb_name,
                           char *dictionary_name,
                           int *fpdbsum,char *dic_res_name,
/* v.5.0--> */
/*                           char *molec_type,char *reaction_ref) */
                           char *molec_type,char *reaction_ref,
                           int *match_names,char *work_dir,
                           char *pdbsum_dir)
/* <--v.5.0 */
{
  int iarg, namlen;
  int for_pdbsum;

  /* Initialise variables */
  dic_res_name[0] = '\0';
  for_pdbsum = FALSE;
  *molec_type = ' ';
/* v.5.0--> */
  *match_names = TRUE;
  Match_Score = MATCH_SCORE;
  work_dir[0] = '\0';
/* <--v.5.0 */

  /* If no command-line arguments entered, then report error */
  if (num_args < 2)
    {
      printf("\n");
      printf("*** ERROR. Program hbadd requires at least 2 arguments:\n");
      printf("\n");
      printf("        hbadd  pdb_file  dic_file [-pdbsum]  "
             "[-c dic_name]  [-mtype xY.Y.Y]\n");
/* v.5.0--> */
      printf("               [-nomatch]  [-wkdir directory]\n");
/* <--v.5.0 */
      printf("\n");
      printf("    where:\n");
      printf("      pdb_file = name of PDB file to be processed\n");
      printf("      dic_file = name of the Het Group Dictionary\n");
      printf("      -pdbsum  = run for PDBsum, so hetxx.pdb files to "
             "be written\n");
      printf("      -c dic_name\n");
      printf("               = compare Het Groups in the PDB file "
             "vs given dictionary entry\n");
      printf("      -mtype x = molecule type (C=cofactor, R=reactant,"
             "P=product)\n");
      printf("      Y.Y.Y    = reaction reference (E.C. ref, reaction "
             "step, component no)\n");
/* v.5.0--> */
      printf("      -nomatch = flag indicating that atom name match not "
             "required\n");
      printf("      -wkdir directory = working directory for output "
             "files\n");
      printf("      -pdir directory = name of PDBsum directory for "
             "rasmol.dfn and ligands.dat files\n");
/* <--v.5.0 */
      printf("\n");
      printf("*** Program terminated with error.\n");
      exit(1);
    }

  /* Get the PDB file name */

  /* Check that length of filename not too long */
  namlen = strlen(strptrs[1]);
  if (namlen > FILENAME_LEN - 1)
    {
      printf("\n");
      printf("*** ERROR. Length of PDB filename longer than %d characters\n",
             FILENAME_LEN);
      printf("***        [%s]\n",strptrs[1]);
      printf("***        Program hbadd terminated with error.\n");
      printf("\n");
      exit(1);
    }

  /* Store the PDB filename */
  strcpy(pdb_name,strptrs[1]);

  /* Get the Het Group Dictionary file name */

  /* Check that length of filename not too long */
  namlen = strlen(strptrs[2]);
  if (namlen > FILENAME_LEN - 1)
    {
      printf("\n");
      printf("*** ERROR. Length of dictionary filename longer ");
      printf("than %d characters\n",FILENAME_LEN);
      printf("***        [%s]\n",strptrs[2]);
      printf("***        Program hbadd terminated with error.\n");
      printf("\n");
      exit(1);
    }

  /* Store the PDB filename */
  strcpy(dictionary_name,strptrs[2]);

  /* Loop through any additional parameters */
  for (iarg = 3; iarg < num_args + 1; iarg++)
    {
      /* Check for -pdbsum option */
      if (!strcmp(strptrs[iarg],"-pdbsum"))
        for_pdbsum = TRUE;

      /* Check for -c option signifying comparison of specifically
         named residues in the PDB file against entries in the Het
         Group Dictionary */
      else if (!strcmp(strptrs[iarg],"-c"))
        {
          /* Get the dictionary entry to be compared against */
          if (iarg + 1 < num_args + 1)
            {
              /* Store name */
              strcpy(dic_res_name,strptrs[iarg + 1]);
              iarg++;
            }
        }

      /* Check for -mtype option giving molecule type (C=cofactor,
         R=reactant, P=product) and enzyme reference */
      else if (!strcmp(strptrs[iarg],"-mtype"))
        {
          /* Get the molecule type */
          if (iarg + 1 < num_args + 1)
            {
              /* Store name and reaction reference */
              *molec_type = strptrs[iarg + 1][0];
              strcpy(reaction_ref,strptrs[iarg + 1]+1);
              iarg++;
            }
        }

/* v.5.0--> */
      /* If -nomatch flag entered, don't want to score similar atom
         names */
      else if (!strcmp(strptrs[iarg],"-nomatch"))
        {
          *match_names = FALSE;
          Match_Score = MAP_SCORE;
        }

      /* Get the working directory */
      else if (!strcmp(strptrs[iarg],"-wkdir"))
        {
          /* Get the directory from the next command-line parameter */
          if (iarg + 1 < num_args + 1)
            {
              /* Get the directory name */
              strcpy(work_dir,strptrs[iarg + 1]);
            }
          iarg++;
        }

      /* Get the PDBsum directory, if entered */
      else if (!strcmp(strptrs[iarg],"-pdir"))
        {
          /* Get the directory from the next command-line parameter */
          if (iarg + 1 < num_args + 1)
            {
              /* Get the directory name */
              strcpy(pdbsum_dir,strptrs[iarg + 1]);
            }
          iarg++;
        }
/* <--v.5.0 */
    }

  /* Return PDBsum flag */
  *fpdbsum = for_pdbsum;
}
/* v.5.0--> */
/***********************************************************************

prefix_file_names  -  Prefix all input and output files with the name of
                      the working directory 

***********************************************************************/

void prefix_file_names(char *work_dir,char *bond_name,char *match_file,
                       char *out_name,char *sum_name,char *map_name)
{
  char ch, tmp[FILENAME_LEN];
  int len;

  /* Get the length of the directory name */
  len = strlen(work_dir);

  /* If empty string, then return */
  if (len == 0) return;

  /* Check whether directory name is correctly terminated by a
     directory separator */
  ch = work_dir[len - 1];
  if (ch != '/' && ch != '\\')
    strcat(work_dir,DIR_SEP);

  /* Prefix all the file names */
  strcpy(tmp,work_dir);
  strcat(tmp,bond_name);
  strcpy(bond_name,tmp);

  strcpy(tmp,work_dir);
  strcat(tmp,match_file);
  strcpy(match_file,tmp);

  strcpy(tmp,work_dir);
  strcat(tmp,out_name);
  strcpy(out_name,tmp);

  strcpy(tmp,work_dir);
  strcat(tmp,sum_name);
  strcpy(sum_name,tmp);

  strcpy(tmp,work_dir);
  strcat(tmp,map_name);
  strcpy(map_name,tmp);
}
/* <--v.5.0 */
/***********************************************************************

to_lower  -  Convert given character string to lower case

***********************************************************************/

void to_lower(char *search_string,int length)
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

string_truncate  -  Truncate the given string at the last non-blank
                    character, or at its maximum length

***********************************************************************/

int string_truncate(char *string,int max_length)
{
  char ch;

  int iend, ipos, len;

  /* Initialise variables */
  iend = -1;
  ch = string[0];
  for (ipos = 0; ipos < max_length && ch != '\0'; ipos++)
    {
      /* Get the current character and check for end of string */
      ch = string[ipos];
      if (ch == '\n' || ch == 13)
        ch = '\0';
      if (ch != '\0')
        {
          /* If not a space, then mark end of string */
          if (ch != ' ')
            iend = ipos;
        }
    }

  /* End the string after the last non-blank character */
  if (iend + 1 < max_length)
    {
      string[iend + 1] = '\0';
      len = iend + 1;
    }
  else
    {
      string[max_length - 1] = '\0';
      len = max_length - 1;
    }

  /* Return the string length */
  return(len);
}
/***********************************************************************

store_string  -  Allocate memory for given string and store pointer to
                 it

***********************************************************************/

void store_string(char **string_ptr,char *string)
{
  int len;

  /* Get the length of the string */
  len = strlen(string);

  /* Allocate memory for the string */
  *string_ptr = (char *) malloc(sizeof(char)*(len + 1));

  /* Store the string */
  strcpy(*string_ptr,string);
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

int get_ligand_residues(char *input_line,struct rasdef **fst_rasdef_ptr,
                        struct rasdef **lst_rasdef_ptr,int number,
                        int ndef)
{
  char chain, chain_str[2], res_name[4], res_num[6];
  char colour[COLNAM_LEN], number_string[10];
  char token_string[TOKEN_LEN + 1], tmp_string[TOKEN_LEN + 1];
  char *string_ptr;
  char ligand_from[7], ligand_to[7];
  char first_residue[7], last_residue[7];

  int count, len, lpos, token_len, type;
  int name_end;
  int done;

  struct rasdef *rasdef_ptr;
  struct rasdef *first_rasdef_ptr, *last_rasdef_ptr;

  /* Initialise variables */
  chain = ' ';
  count = 0;
  done = FALSE;
  type = LIGAND;
  colour[0] = '\0';
  first_residue[0] = '\0';
  last_residue[0] = '\0';
  ligand_from[0] = '\0';
  ligand_to[0] = '\0';

  /* Initialise rasdef pointers */
  rasdef_ptr = NULL;
  first_rasdef_ptr = *fst_rasdef_ptr;
  last_rasdef_ptr = *lst_rasdef_ptr;

  /* Get the line length */
  len = strlen(input_line);
  lpos = 0;

  /* Loop until line has been processed */
  while (done == FALSE)
    {
      /* Get the next token on the line */
      token_len = get_token(input_line,LINE_LEN,&lpos,token_string,
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

          /* Right-justify residue name, if necessary */
          if (strlen(res_name) == 1)
            {
              strcpy(tmp_string,"  ");
              strcat(tmp_string,res_name);
              strcpy(res_name,tmp_string);
            }
          else if (strlen(res_name) == 2)
            {
              strcpy(tmp_string," ");
              strcat(tmp_string,res_name);
              strcpy(res_name,tmp_string);
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

          /* If this is the first residue of this ligand, store its
             number */
          if (ligand_from[0] == '\0')
            strcpy(ligand_from,res_num);

          /* Create a rasdef entry for this ligand residue */
          rasdef_ptr
            = create_rasdef_entry(&first_rasdef_ptr,
                                  &last_rasdef_ptr,type,number,
                                  res_name,res_num,chain,colour,
                                  input_line,ligand_from,ligand_to,
                                  count);
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
      rasdef_ptr
        = create_rasdef_entry(&first_rasdef_ptr,
                              &last_rasdef_ptr,type,number,
                              res_name,res_num,chain,colour,
                              input_line,ligand_from,ligand_to,
                              count);
      ndef++;
    }

  /* Return current pointers */
  *fst_rasdef_ptr = first_rasdef_ptr;
  *lst_rasdef_ptr = last_rasdef_ptr;

  /* Return updated number of rasdef entries */
  return(ndef);
}
/***********************************************************************

get_dfn_data  -  Read in and store the RasMol definitions from the
                 rasmol.dfn file

***********************************************************************/

int get_dfn_data(char *pdb_code,char *pdbsum_dir,
                 struct rasdef **fst_rasdef_ptr)
{
  char input_line[LINE_LEN + 1], file_name[FILENAME_LEN];
  char chain, res_name[4], res_num[6];
  char colour[COLNAM_LEN];
  char ligand_from[7], ligand_to[7];

  int count, line, nchain, ndef, nligand, nmetal, type;
  int first;

  struct rasdef *rasdef_ptr;
  struct rasdef *first_rasdef_ptr, *last_rasdef_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  first = TRUE;
  count = 0;
  file_name[0] = '\0';
  ligand_from[0] = '\0';
  ligand_to[0] = '\0';
  line = 0;
  nchain = 0;
  ndef = 0;
  nligand = 0;
  nmetal = 0;

  /* Form the file name of the rasmol.dfn file */
  if (pdbsum_dir[0] != '\0')
    {
      strcpy(file_name,pdbsum_dir);
      strcat(file_name,DIR_SEP);
      strcat(file_name,pdb_code);
      strcat(file_name,DIR_SEP);
    }
  strcat(file_name,"rasmol.dfn");

  /* If have an uploaded PDB file, the file will be in the current directory */
  if (!strcmp(pdb_code,"uplo"))
    strcpy(file_name,"rasmol.dfn");

  /* Initialise rasdef pointers */
  rasdef_ptr = NULL;
  first_rasdef_ptr = NULL;
  last_rasdef_ptr = NULL;

  /* Open the rasmol.dfn input file */
  if ((file_ptr = fopen(file_name,"r")) !=  NULL)
    {
      /* Loop through the file processing each record */
      while (fgets(input_line,LINE_LEN,file_ptr)!= NULL)
        {
          /* Increment line-count */
          line++;

          /* Truncate the string to strip off the newline */
          string_truncate(input_line,LINE_LEN);

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
              ndef
                = get_ligand_residues(input_line+4,&first_rasdef_ptr,
                                      &last_rasdef_ptr,nligand,ndef);
            }

          /* Check for metal */
          else if (input_line[0] == 'M')
            {
              /* Define type */
              type = METAL;
              nmetal++;

              /* Extract the metal name into the residue number field */
              strncpy(res_num,input_line+1,3);
              res_num[3] = ' ';
              res_num[4] = '\0';

              /* Get the metal colour */
              strncpy(colour,input_line+5,COLNAM_LEN - 1);
              colour[COLNAM_LEN - 1] = '\0';

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

store_connections  -  Process the current CONECT record and store its
                      connections if they relate to any of the HET groups

***********************************************************************/

void store_connections(char line[LINE_LEN + 1],int first_atom_number, 
                       int last_atom_number,
                       struct connect **first_connect_ptr,
                       struct connect **last_connect_ptr,int *nlinks)
{
  char atmnum[6];
  int done, inhetgroup;
  int iatom, icon, ipos, jatom, ncon;
  int inhetgroup_atom[8], store_atom[8];

  struct connect *connect_ptr;
  
  /* Initialise variables */
  inhetgroup = 0;
  ipos = 6;
  ncon = 0;
  done = FALSE;

  /* Loop through the character positions in the record to extract
     all the atom numbers in the CONECT record */
  while (done == FALSE && ncon < 5)
    {
      /* Extract this atom number and store */
      strncpy(atmnum,line+ipos,5);
      atmnum[5] = '\0';

      /* If have blank, then end of connections reached */
      if (!strncmp(atmnum,"     ",5))
        done = TRUE;
      
      /* Otherwise, store the number */
      else
        {
          store_atom[ncon] = atoi(atmnum);
          
          /* Check whether the atom might belong to one of the HET groups */
          if (store_atom[ncon] >= first_atom_number &&
              store_atom[ncon] <= last_atom_number)
            {
              inhetgroup++;
              inhetgroup_atom[ncon] = TRUE;
            }
          else
            {
              inhetgroup_atom[ncon] = FALSE;
            }

          /* Increment variables */
          ncon++;
          ipos = ipos + 5;
        }
    }

  /* If have possible covalent bonds between HET group atoms, then store */
  if (inhetgroup > 1)
    {
      /* Consider each connection in turn */
      iatom = store_atom[0];
      for (icon = 1; icon < ncon; icon++)
        {
          jatom = store_atom[icon];

          /* If both atoms are in a HET groups then store this connection */
          if ((inhetgroup_atom[0] == TRUE && inhetgroup_atom[icon] == TRUE))
            {
              /* Allocate memory for structure to hold current atom's
                 details */
              connect_ptr
                = (struct connect*)malloc(sizeof(struct connect));
              if (connect_ptr == NULL)
                {
                  printf("*** ERROR. Unable to allocate memory for struct");
                  printf(" connect\n");
                  printf("***        Program hbadd terminated with error.\n");
                  exit (1);
                }

              /* If this is the first connection stored, then update
                 pointer to head of linked list */
              if (*first_connect_ptr == NULL)
                *first_connect_ptr = connect_ptr;
              
              /* Otherwise, make previous connection point to the current
                 one */
              else
                (*last_connect_ptr)->next_connect_ptr = connect_ptr;

              /* Save the current record's pointer */
              *last_connect_ptr = connect_ptr;

              /* Store the atom numbers of the two atoms */
              connect_ptr->atom_number1 = iatom;
              connect_ptr->atom_number2 = jatom;
              connect_ptr->next_connect_ptr = NULL;

              /* Increment number of connections stored */
              (*nlinks)++;
            }
        }
    }
}
/***********************************************************************

check_for_metal  -  Check whether the given atom is a metal ion

***********************************************************************/

int check_for_metal(char atom_name[5],char element[3])
{
  int imetal;
  int metal;

  static int nmetals = 100;

  /* Note "AC  " changed to "ac  " below as carbon atoms in certain
     molecules labelled AC1, AC1, etc (eg in FAD). Similarly NP, NO
     and PO */
  static char *metal_name[] = {
    "HE  ", "LI  ", "BE  ", " B  ", " F  ",
    "NE  ", "NA  ", "MG  ", "AL  ", "SI  ",
    "CL  ", "AR  ", " K  ", "CA  ", "SC  ",
    "TI  ", " V  ", "CR  ", "MN  ", "FE  ",
    "CO  ", "NI  ", "CU  ", "ZN  ", "GA  ",
    "GE  ", "AS  ", "SE  ", "BR  ", "KR  ",
    "RB  ", "SR  ", " Y  ", "ZR  ", "NB  ",
    "MO  ", "TC  ", "RU  ", "RH  ", "PD  ",
    "AG  ", "CD  ", "IN  ", "SN  ", "SB  ",
    "TE  ", " I  ", "XE  ", "CS  ", "BA  ",
    "LA  ", "CE  ", "PR  ", "ND  ", "PM  ",
    "SM  ", "EU  ", "GD  ", "TB  ", "DY  ",
    "HO  ", "ER  ", "TM  ", "YB  ", "LU  ",
    "HF  ", "TA  ", " W  ", "RE  ", "OS  ",
    "IR  ", "PT  ", "AU  ", "HG  ", "TL  ",
    "PB  ", "BI  ", "po  ", "AT  ", "RN  ",
    "FR  ", "RA  ", "ac  ", "TH  ", "PA  ",
    " U  ", "np  ", "PU  ", "AM  ", "CM  ",
    "BK  ", "CF  ", "ES  ", "FM  ", "MD  ",
    "no  ", "LR  ", "....", "....", "...."
  };

  /* Initialise variables */
  metal = FALSE;

  /* Loop through all the metals */
  for (imetal = 0; imetal < nmetals && metal == FALSE; imetal++)
    {
      /* If name matches in first two characters, and third character
         is either a space or a number, then have a metal */
/* v.5.3.7 --> */
//      if (!strncmp(atom_name,metal_name[imetal],2) &&
      if (!strncmp(element,metal_name[imetal],2) ||
          (!strncmp(atom_name,metal_name[imetal],2) &&
/* <-- v.5.3.7 */
          (atom_name[2] == ' ' ||
           (atom_name[2] >= '0' && atom_name[2] <= '9' &&
/* v.5.3.7 --> */
//            atom_name[3] == ' ')))
            (atom_name[3] == ' ' ||
             (atom_name[3] >= '0' && atom_name[3] <= '9'))))))
/* <-- v.5.3.7 */
        {
          /* Set flag */
          metal = TRUE;

          /* Return the metal name matched */
/* v.5.3.7 --> */
          if (element[0] == '\0')
            {
/* <-- v.5.3.7 */
              strncpy(element,metal_name[imetal],2);
              element[2] = '\0';
/* v.5.3.7 --> */
            }
/* <-- v.5.3.7 */
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

int is_hydrogen_atom(char atname[5])
{
  char atom_name[5], element[3];

  int is_hydrogen, is_metal;

  /* Initialise */
  is_hydrogen = FALSE;
  strncpy(atom_name,atname,4);
  atom_name[4] = '\0';
/* v.5.3.7 --> */
  element[0] = '\0';
/* <-- v.5.3.7 */

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
      is_metal = check_for_metal(atom_name,element);

      /* If a metal, then isn't a hydrogen */
      if (is_metal == TRUE)
        is_hydrogen = FALSE;
    }

  return is_hydrogen;
}
/***********************************************************************

get_element  -  Determine element from atom's name

***********************************************************************/

void get_element(char atom_name[5],char element[3])
{
  char metal_name[3];

  int ielem;
  int found, is_hydrogen, is_metal;

  static char *valid_elem[] = { " C", " N", " O", " P", " S" };

  static int nelem = 5;

  /* Initialise variables */
  strcpy(element,"  ");
/* v.5.3.7 --> */
  metal_name[0] = '\0';
/* <-- v.5.3.7 */

  /* Check whether this is a hydrogen atom */
  is_hydrogen = is_hydrogen_atom(atom_name);

  /* If atom is a hydrogen, then have its element type */
  if (is_hydrogen == TRUE)
    strcpy(element," H");

  /* Otherwise, identify non-hydrogen element */
  else
    {
      /* First, whether this might be a metal */
      is_metal = check_for_metal(atom_name,metal_name);

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
}
/***********************************************************************

check_element  -  Match whether the atom name is consistent with the
                  element definition

***********************************************************************/

int check_element(char atom_name[5],char element[3],
                  char res_name[4],char res_num[6],char chain)
{
  int match;

  /* Initialise variables */
  match = TRUE;

  /* If have a single-character element name */
  if (element[0] == ' ')
    {
      /* If element matches atom's second character, then atom name
         OK */
      if (element[1] == atom_name[1])
        match = TRUE;

      /* Otherwise, must match at first character */
      else if (element[1] == atom_name[0])
        match = TRUE;

      /* For hydrogens, the H might be in the third character position */
      else if (!strncmp(element," H",2) && atom_name[2] == 'H')
        match = TRUE;

      /* Otherwise, have a mismatch */
      else
        match = FALSE;
    }

  /* If have a two-character element name (ie metal), then must
     match first two characters of the atom name */
  else if (!strncmp(element,atom_name,2))
    match = TRUE;

  /* Except, of course, in special circumstances when the match
     is in the third and fourth character positions (!) */
  else if (!strncmp(element,atom_name+1,2))
    match = TRUE;

  /* Otherwise, have another likely mismatch */
  else
    match = FALSE;

  /* If element not consistent with atom name, then show warning
     message */
  if (match == FALSE)
    {
      printf("*** Warning. Can't match element identifier [%s]",
             element);
      printf(" to atom name [%s]\n",atom_name);
      printf("***          for residue [%s %s %c]\n",
             res_name,res_num,chain);
    }

  /* Return whether we have a match */
  return(match);
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

      /* Get next definition record */
      rasdef_ptr = rasdef_ptr->next_rasdef_ptr;
    }

  /* Return the object to which this residue is to be assigned */
  return(assigned_rasdef_ptr);
}
/***********************************************************************

check_residue_name  -  Check whether residue name is one of the standard
                       amino acids or a nucleic acid base

***********************************************************************/

int check_residue_name(char res_name[4])
{
  int iamino;
  static int namino = 0;
  int done, standard;

  static char *amino[] = {
    "ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS",
    "ILE", "LYS", "LEU", "MET", "ASN", "PRO", "GLN",
    "ARG", "SER", "THR", "VAL", "TRP", "TYR",
    "  A", "  C", "  G", "  I", "  T", "  U",
    " +A", " +C", " +G", " +I", " +T", " +U",
    " DA", " DC", " DG", " DI", " DT", " DU", "XXX"
  };

  /* Initialise flag */
  standard = FALSE;

  /* If this is the first call to this routine, count the number of
     standard residue names */
  if (namino == 0)
    {
      /* Initialise flag */
      done = FALSE;

      /* Loop through the codes until reach the last one */
      while (done == FALSE)
        {
          /* Check for the last code */
          if (!strncmp(amino[namino],"XXX",3))
            done = TRUE;

          /* Otherwise, increment the count */
          else
            namino++;
        }
    }

  /* Check whether residue appears in the list */
  for (iamino = 0; iamino < namino && standard == FALSE; iamino++)
    if (!strncmp(res_name,amino[iamino],3))
      standard = TRUE;

  /* Return whether residue is one of the standard one */
  return(standard);
}
/***********************************************************************

is_duplicate  -  Check whether we already have an atom of this name for
                 the current residue

***********************************************************************/

int is_duplicate(char atom_name[5],struct atom *first_atom_ptr)
{
  int duplicate;

  struct atom *atom_ptr;

  /* Initialise variables */
  duplicate = FALSE;

  /* Initialise pointer to start of linked list of atom coords */
  atom_ptr = first_atom_ptr;

  /* Loop through all the atoms to check whether we
     already have this one */
  while (atom_ptr != NULL && duplicate == FALSE)
    {
      /* If atom names match, then have a duplicate */
      if (!strcmp(atom_ptr->atom_name,atom_name))
        duplicate = TRUE;

      /* Get pointer to next atom in linked-list */
      atom_ptr = atom_ptr->next_atom_ptr;
    }               

  /* Return whether we have a duplicate atom name */
  return(duplicate);
}
/***********************************************************************

create_residue_entry  -  Create a residue entry

***********************************************************************/

struct residue *create_residue_entry(struct residue **fst_residue_ptr,
                                     struct residue **lst_residue_ptr,
                                     char *res_name,char *res_num,
                                     char chain,struct atom *atom_ptr,
                                     struct rasdef *rasdef_ptr)
{
  struct residue *residue_ptr;
  struct residue *first_residue_ptr, *last_residue_ptr;

  /* Initialise residue pointers */
  residue_ptr = NULL;
  first_residue_ptr = *fst_residue_ptr;
  last_residue_ptr = *lst_residue_ptr;

  /* Create residue entry */
  residue_ptr =
    (struct residue*) malloc(sizeof(struct residue));
  if (residue_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for");
      printf(" struct residue\n");
      printf("***        Program postssm terminated with");
      printf(" error.\n");
      exit (1);
    }

  /* Store the details */
  strncpy(residue_ptr->res_name,res_name,3);
  residue_ptr->res_name[3] = '\0';
  strncpy(residue_ptr->res_num,res_num,5);
  residue_ptr->res_num[5] = '\0';
  residue_ptr->chain = chain;

  /* Initialise other fields */
  residue_ptr->natoms = 0;
  residue_ptr->nconnect = 0;
  residue_ptr->non_hydrogen = 0;
  residue_ptr->nmatch = 0;
  residue_ptr->fake = FALSE;
  residue_ptr->nmatched = 0;
  residue_ptr->nmapped = 0;
  residue_ptr->nmissed = 0;
  residue_ptr->nmissing = 0;
  residue_ptr->nbonds_matched = 0;
  residue_ptr->liguniq_no = 0;
  residue_ptr->ligand_no = 0;
  residue_ptr->first_atom_ptr = atom_ptr;
  residue_ptr->dictionary_ptr = NULL; 
  residue_ptr->rasdef_ptr = rasdef_ptr; 
  residue_ptr->next_residue_ptr = NULL;

  /* Initialise pointer to next entry */
  residue_ptr->next_residue_ptr = NULL;

  /* If this is the first entry stored, then update
     pointer to head of linked list */
  if (first_residue_ptr == NULL)
    first_residue_ptr = residue_ptr;

  /* Otherwise, make previous entry point to the current
     one */
  else
    last_residue_ptr->next_residue_ptr = residue_ptr;
                      
  /* Save the current residue's pointer */
  last_residue_ptr = residue_ptr;

  /* Return current pointers */
  *fst_residue_ptr = first_residue_ptr;
  *lst_residue_ptr = last_residue_ptr;
  return(residue_ptr);
}
/***********************************************************************

read_pdbfile  -  Read through PDB file to get pick up all the HET groups

***********************************************************************/

void read_pdbfile(char *pdb_name,char *pdb_code,char *pdbsum_dir,
                  char *dic_res_name,struct atom **first_atom_ptr,
                  struct residue **fst_residue_ptr,
                  struct connect **first_connect_ptr,
                  int *natoms,int *nhetgroups,int *have_hydrogens)
{
  char line[LINE_LEN + 1], pdb_upper[5];
  char chain, last_chain, last_resnum[6], resnum[6];
  char atmnum[6], atom_name[5], element[3], resnam[4];
  char last_resnam[4];

  int first_atom_number, last_atom_number, len, nconnect, nlinks;
  int atom_count, residue_count;
  int ndef, nhydrogens;
  int duplicate, keep_reading, new_resid, standard, want_atom, want_resid;
  int consistent;

  struct atom *atom_ptr, *last_atom_ptr;
  struct connect *last_connect_ptr;
  struct rasdef *first_rasdef_ptr, *rasdef_ptr;
  struct residue *residue_ptr, *first_residue_ptr, *last_residue_ptr;

  *have_hydrogens = FALSE;
  nhydrogens = 0;

  /* Initialise variables */
  first_rasdef_ptr = NULL;
  rasdef_ptr = NULL;
  pdb_upper[0] = '\0';
  pdb_code[0] = '\0';

  /* Open the pdb file */
  printf("  Opening PDB file [%s] ...\n",pdb_name);
  if ((fil_pdb = fopen(pdb_name,"r")) == NULL)
    {
      printf("\n*** Unable to open your PDB file [%s]\n",pdb_name);
      exit(1);
    }
    
  /* Prepare to start reading in the data from the PDB file */
  printf("\n");
  printf("    Searching for HET groups in PDB file ...\n");

  /* Initialise variables */
  atom_count = 0;
  residue_count = 0;
  *first_atom_ptr = NULL;
  first_residue_ptr = NULL;
  residue_ptr = NULL;
  *first_connect_ptr = NULL;
  first_atom_number = 0;
  keep_reading = TRUE;
  last_atom_number = 999999;
  last_atom_ptr = NULL;
  last_residue_ptr = NULL;
  last_connect_ptr = NULL;
  *natoms = 0;
  ndef = 0;
  *nhetgroups = 0;
  nconnect = 0;
  new_resid = TRUE;
  nlinks = 0;
  last_resnam[0] = '\0';
  last_resnum[0] = '\0';
  last_chain = '\0';
  want_atom = FALSE;
  want_resid = FALSE;
  
  /* Search through the PDB file to find all the HET groups */
  while (fgets(line,LINE_LEN,fil_pdb) != NULL)
    {
/* v.4.4.3--> */
      /* Truncate the string to strip off the newline */
      len = string_truncate(line,LINE_LEN);
/* <--v.4.4.3 */

      /* If this is the HEADER record, pick up the PDB code */
      if (!strncmp(line,"HEADER",6))
        {
          /* Extract the PDB code */
          strncpy(pdb_upper,line+62,4);
          pdb_upper[4] = '\0';

          /* Store the PDB code */
          strcpy(pdb_code,pdb_upper);
          to_lower(pdb_code,4);

          /* If comparing a single dictionary entry, get this PDB file's
             molecule definitions from the rasmol.dfn file */
          if (dic_res_name[0] != '\0')
            {
              /* Read in this PDB entry's rasmol.dfn file */
              ndef
                = get_dfn_data(pdb_code,pdbsum_dir,&first_rasdef_ptr);
            }
        }

      /* If this is an HETATM record, process it */
      else if ((!strncmp(line,"HETATM",6) || !strncmp(line,"ATOM  ",6)) &&
               keep_reading == TRUE)
        {
          /* Initialise flag */
          new_resid = FALSE;

          /* Get the residue name, number and chain-ID */
          strncpy(resnam,line+17,3);
          resnam[3] = '\0';
          strncpy(resnum,line+22,5);
          resnum[5] = '\0';
          chain = line[21];
          
          /* Get the atom name */
          strncpy(atom_name,line+12,4);
          atom_name[4] = '\0';

          /* If this is a new residue type, check whether it is a standard
             residue type */
          if (strcmp(resnam,last_resnam) || strcmp(resnum,last_resnum) ||
              chain != last_chain)
            {
              /* Initialise flags */
              want_resid = FALSE;
              new_resid = TRUE;

              /* If have molecule definitions, determine which molecule
                 this residue belongs to */
              if (ndef > 0)
                {
                  /* Get the definition record */
                  rasdef_ptr
                    = assign_to_object(first_rasdef_ptr,atom_name,
                                       resnam,resnum,chain);

                  /* If this is a ligand residue, then set flag */
                  if (rasdef_ptr != NULL && rasdef_ptr->type == LIGAND)
                    want_resid = TRUE;
                }

              /* Otherwise, check whether this is a standard amino acid
                 or nucleic acid base */
              else
                {
                  /* Check residue type */
                  standard = check_residue_name(resnam);

                  /* If non-standard residue, then want it */
                  if (standard == FALSE)
                    want_resid = TRUE;
                }
            }

          /* Check that this is not a water */
          if (!strncmp(resnam,"HOH",3) || !strncmp(resnam,"WAT",3))
            want_resid = FALSE;

          /* If residue is wanted, then probably want atom, too */
          if (want_resid == TRUE)
            want_atom = TRUE;
          else
            want_atom = FALSE;

          /* If looks like atom is wanted, then check that we don't already
             have an atom of that name */
          if (want_resid == TRUE && new_resid == FALSE &&
              residue_ptr != NULL)
            {
              duplicate =
                is_duplicate(atom_name,residue_ptr->first_atom_ptr);
              if (duplicate == TRUE)
                want_atom = FALSE;
              else
                want_atom = TRUE;
            }

          /* If atom is wanted, then store it */
          if (want_atom == TRUE)
            {
              /* Store the current atom's details */

              /* Allocate memory for structure to hold current atom's
                 details */
              atom_ptr
                = (struct atom*)malloc(sizeof(struct atom));
              if (atom_ptr == NULL)
                {
                  printf("*** ERROR. Unable to allocate memory for struct");
                  printf(" atom\n");
                  printf("***        Program hbadd terminated with error.\n");
                  exit (1);
                }

              /* If this is the first atom stored, then update
                 pointer to head of linked list */
              if (*first_atom_ptr == NULL)
                *first_atom_ptr = atom_ptr;
              
              /* Otherwise, make previous residue point to the current
                 one */
              else
                last_atom_ptr->next_atom_ptr = atom_ptr;

              /* Save the current atom's pointer */
              last_atom_ptr = atom_ptr;

              /* Extract the relevant atom data from the record */
              strncpy(atmnum,line+6,5);
              atmnum[5] = '\0';
              atom_ptr->atom_number = atoi(atmnum);
              strcpy(atom_ptr->atom_name,atom_name);
              sscanf(line+30," %f %f %f ",&atom_ptr->x,&atom_ptr->y,
                     &atom_ptr->z);

              /* Check for old-style PDB files containing the PDB code
                 in columns 73-76 */
              if (!strncmp(line+72,pdb_upper,4))
                strcpy(element,"  ");

              /* Otherwise, try to extract the element from the new-style
                 files */
/* v.4.4.3--> */
              //else
              else if (len > 76)
/* <--v.4.4.3 */
                {
                  /* Get the element identifier, if there is one */
                  strncpy(element,line+76,2);
                  element[2] = '\0';

                  /* Check whether this is a valid element */
                  if ((element[0] >= '0' && element[0] <= '9') ||
                      (element[1] >= '0' && element[1] <= '9'))
                    strcpy(element,"  ");
                }
/* v.4.4.3--> */
              else
                strcpy(element,"  ");
/* <--v.4.4.3 */

              /* If have an element identifier, check whether it is
                 consistent with the atom name */
              if (strcmp(element,"  "))
                {
                  /* Compare element with atom name */
                  consistent = check_element(atom_name,element,resnam,
                                             resnum,chain);

                  /* If not consistent, assume atom name gives the
                     correct value */
                  if (consistent == FALSE)
                    strcpy(element,"  ");
                }

              /* If don't have the element, have to extract it from
                 the atom name */
              if (!strcmp(element,"  "))
                get_element(atom_ptr->atom_name,element);

              /* Store element */
              strcpy(atom_ptr->element,element);

              /* If this is a hydrogen atom, increment count */
              if (!strcmp(element," H"))
                {
                  nhydrogens++;
                  atom_ptr->hydrogen = TRUE;
                }
              else
                atom_ptr->hydrogen = FALSE;

              /* Initialise other fields */
              atom_ptr->nequiv = 0;
              atom_ptr->next_atom_ptr = NULL;

              /* If this is the first stored atom, store its number */
              if (first_atom_number == 0)
                first_atom_number = atom_ptr->atom_number;

              /* Store number of last atom read in */
              last_atom_number = atom_ptr->atom_number;
              
              /* Check whether this is a new residue */
              if (new_resid == TRUE)
                {
                  /* Store this residue's details */
                  residue_ptr
                    = create_residue_entry(&first_residue_ptr,
                                           &last_residue_ptr,resnam,
                                           resnum,chain,atom_ptr,
                                           rasdef_ptr);

                  /* Save current residue */
                  strncpy(last_resnam,resnam,3);
                  last_resnam[3] = '\0';
                  strncpy(last_resnum,resnum,5);
                  last_resnum[5] = '\0';
                  last_chain = chain;
                  
                  /* Increment the residue-count */
                  residue_count++;
                }

              /* Save pointer to residue in current atom's record */
              atom_ptr->residue_ptr = residue_ptr;
              atom_ptr->dic_data_ptr = NULL;
              atom_count++;

              /* Increment count of atoms in this residue */
              residue_ptr->natoms++;

              /* If this is not a hydrogen, then increment count of
                 non-hydrogen atoms */
              if (atom_ptr->hydrogen == FALSE)
                residue_ptr->non_hydrogen++;
            }

          /* Save current residue */
          strncpy(last_resnam,resnam,3);
          last_resnam[3] = '\0';
        }

      /* If this is an ENDMDL record, then don't want to read in any more
         HET groups */
      else if (!strncmp(line,"ENDMDL",6))
        keep_reading = FALSE;

      /* If this is an CONECT record, process it */
      else if (!strncmp(line,"CONECT",6))
        {
          /* If this is the first CONECT record encountered, then
             show message */
          if (nconnect == 0)
            {
              printf("\n");
              printf("    Processing CONECT records ...\n");
            }
          nconnect++;

          /* Check whether first atom corresponds to any atom within
             the range of atom numbers stored */

          /* Process the current CONECT record */
          store_connections(line,first_atom_number,last_atom_number,
                            &(*first_connect_ptr),&last_connect_ptr,&nlinks);
        }
    }

  /* Print counts of data read in */
  printf("\n");
  printf("      Number of HET residues read in        = %7d\n",residue_count);
  *nhetgroups = residue_count;
  printf("      Number of atom records stored         = %7d\n",atom_count);
  *natoms = atom_count;
  printf("      Number of CONECT details stored       = %7d\n",nlinks);
  printf("\n");

  /* If have an adequate number of hydrogen atoms, then set flag
     accordingly */
  if (atom_count > 0)
    if ((float) nhydrogens / (float) atom_count > 0.2)
      *have_hydrogens = TRUE;

  /* Return the first residue pointer */
  *fst_residue_ptr = first_residue_ptr;

  /* Close the PDB file */
  fclose(fil_pdb);
}
/***********************************************************************

extract_range  -  Extract the ligand range from the ligands.dat line

***********************************************************************/

void extract_range(char *string,char *res_name1,char *res_num1,
                   char *chain1,char *res_name2,char *res_num2,
                   char *chain2)
{
  char number_string[20], token_string[TOKEN_LEN + 1];

  int itoken, len, lpos, token_len;
  int done;

  /* Initialise variables */
  res_name1[0] = '\0';
  res_num1[0] = '\0';
  *chain1 = ' ';
  res_name2[0] = '\0';
  res_num2[0] = '\0';
  *chain2 = ' ';
  done = FALSE;
  itoken = 0;

  /* Get the line length */
  len = strlen(string);
  lpos = 0;

  /* Loop until string has been processed */
  while (done == FALSE)
    {
      /* Get the next token on the line */
      token_len = get_token(string,LINE_LEN,&lpos,token_string,TOKEN_LEN,
                            &done);

      /* If have token, then interpret it */
      if (token_len > 0)
        {
          /* Process according to which token this is */
          switch (itoken)
            {
              /* First residue name */
            case 0 :

              /* Save as start and end of range */
              strcpy(res_name1,token_string);
              strcpy(res_name2,token_string);

              break;

              /* First residue number and chain if */
            case 1 :

              /* Check for appended chain identifier */
              if (token_len > 3 && token_string[token_len - 1] == ')')
                {
                  /* Get the chain identifier */
                  *chain1 = token_string[token_len - 2];
                  *chain2 = *chain1;

                  /* Remove chain identifier from string */
                  token_string[token_len - 3] = '\0';
                }

              /* Get the residue name */
              strcpy(number_string,token_string);

              /* Set in standard format */
              format_res_number(number_string,res_num1);

              /* Store as end residue, also */
              strcpy(res_num2,res_num1);
      
              break;

              /* The "to" separator */
            case 2 :

              /* Check that have correct text */
              if (strcmp(token_string,"to"))
                printf("*** Warning. Unexpected text in LIG_RANGE: [%s]\n",
                       token_string);
      
              break;

              /* Second residue name */
            case 3 :
      
              /* Save as end of range */
              strcpy(res_name2,token_string);
              break;

              /* Second residue number and chain if */
            case 4 :
      
              /* Check for appended chain identifier */
              if (token_len > 3 && token_string[token_len - 1] == ')')
                {
                  /* Get the chain identifier */
                  *chain2 = token_string[token_len - 2];

                  /* Remove chain identifier from string */
                  token_string[token_len - 3] = '\0';
                }

              /* Get the residue name */
              strcpy(number_string,token_string);

              /* Set in standard format */
              format_res_number(number_string,res_num2);

              break;

              /* Have unexpected token on the LIG_RANGE line */
            default :
              printf("*** Warning. Extra field in LIG_RANGE: [%s]\n",
                     token_string);
            }
        }

      /* Increment token count */
      itoken++;
    }
}
/***********************************************************************

get_ligands  -  Read in the ligands.dat file to identify which ligand
                each residue belongs to

***********************************************************************/

void get_ligands(char *pdb_code,struct residue *first_residue_ptr,
                 char *pdbsum_dir)
{
  char input_line[LINE_LEN + 1], file_name[FILENAME_LEN];
  char number_string[20];
  char chain1, res_name1[4], res_num1[6];
  char chain2, res_name2[4], res_num2[6];

  char *string_ptr;

  int liguniq_no, ligand_no;

  struct rasdef *rasdef_ptr;
  struct residue *residue_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  file_name[0] = '\0';

  /* Form the file name of the ligands.dat file */
  if (pdbsum_dir[0] != '\0')
    {
      strcpy(file_name,pdbsum_dir);
      strcat(file_name,DIR_SEP);
      strcat(file_name,pdb_code);
      strcat(file_name,DIR_SEP);
    }
  strcat(file_name,"ligands.dat");

  /* Open the ligands.dat input file */
  if ((file_ptr = fopen(file_name,"r")) !=  NULL)
    {
      /* Loop through the file processing each record */
      while (fgets(input_line,LINE_LEN,file_ptr)!= NULL)
        {
          /* Truncate the string to strip off the newline */
          string_truncate(input_line,LINE_LEN);

          /* Check for line containing ligand range */
          if (!strncmp(input_line,"LIG_RANGE",9))
            {
              /* Find the start of the ligand range */
              string_ptr = strstr(input_line," ");
              if (string_ptr != NULL)
                {
                  /* Extract the ligand range */
                  extract_range(string_ptr+1,res_name1,res_num1,
                                &chain1,res_name2,res_num2,&chain2);

                  /* Get the first stored residue */
                  residue_ptr = first_residue_ptr;

                  /* Loop through all the residues to find the one
                     belonging to the current ligand */
                  while (residue_ptr != NULL)
                    {
                      /* Get the corresponding molcule definition */
                      rasdef_ptr = residue_ptr->rasdef_ptr;

                      /* If valid, check residue range */
                      if (rasdef_ptr != NULL)
                        {
                          /* Check whether this is the correct ligand */
                          if (rasdef_ptr->chain == chain1 &&
                              !strcmp(rasdef_ptr->ligand_from,res_num1))
                            {
                              /* Initialise ligand reference */
                              liguniq_no = 0;
                              ligand_no = 0;

                              /* Have match, so extract the identifying
                                 ligand numbers */
                              strcpy(number_string,input_line+10);
                              string_ptr = strstr(number_string,"]");
                              if (string_ptr != NULL)
                                string_ptr[0] = '\0';

                              /* Get the unique ligand number */
                              liguniq_no = atoi(number_string);

                              /* Get the second ligand number */
                              string_ptr = strstr(input_line+10,"][");
                              if (string_ptr != NULL)
                                {
                                  strcpy(number_string,string_ptr+2);
                                  string_ptr = strstr(number_string,"]");
                                  if (string_ptr != NULL)
                                    string_ptr[0] = '\0';

                                  /* Get the second ligand number */
                                  ligand_no = atoi(number_string);
                                }

                              /* Store the full ligand ref */
                              residue_ptr->liguniq_no = liguniq_no;
                              residue_ptr->ligand_no = ligand_no;
                            }
                        }

                      /* Get the next residue */
                      residue_ptr = residue_ptr->next_residue_ptr;
                    }
                }
            }

        }

      /* Close the rasdef file */
      fclose(file_ptr);
    }

  /* If no ligands.dat file, then can't continue */
  else
    {
      /* Show error message and stop */
      printf("*** Warning. Unable to open file: %s\n",file_name);
    }
}
/***********************************************************************

match_up_ligands  -  Match the residues just read in to unique ligands
                     in the ligands.dat file deleting any that have no
                     match

***********************************************************************/

void match_up_ligands(char *pdb_code,struct residue **fst_residue_ptr,
                      char *pdbsum_dir)
{
  struct residue *residue_ptr, *next_residue_ptr;
  struct residue *first_residue_ptr, *prev_residue_ptr;

  /* Initialise residue pointers */
  residue_ptr = NULL;
  first_residue_ptr = *fst_residue_ptr;
  prev_residue_ptr = NULL;

  /* Read in the ligands from the ligands.dat file and match to the
     residues stored */
  get_ligands(pdb_code,first_residue_ptr,pdbsum_dir);

  /* Get pointer to the first residue */
  residue_ptr = first_residue_ptr;

  /* Loop through all the residues to delete any without a ligand
     match */
  while (residue_ptr != NULL)
    {
      /* Get the next residue */
      next_residue_ptr = residue_ptr->next_residue_ptr;

      /* If have no ligand match, delete this entry */
      if (residue_ptr->liguniq_no == 0)
        {
          /* If this is not the first residue, get previous residue
             to point round it */
          if (prev_residue_ptr != NULL)
            prev_residue_ptr->next_residue_ptr = next_residue_ptr;

          /* Otherwise, this must be the first residue, so make the
             next residue the new start-point */
          else
            first_residue_ptr = next_residue_ptr;

          /* Free memory taken by this residue */
          free(residue_ptr);
        }

      /* Otherwise, keeping this residue, so save its pointer */
      else
        prev_residue_ptr = residue_ptr;


      /* Get the next residue */
      residue_ptr = next_residue_ptr;
    }

  /* Return current pointers */
  *fst_residue_ptr = first_residue_ptr;
}
/* v.4.3--> */
/***********************************************************************

create_dictionary_entry  -  Create a dictionary entry

***********************************************************************/

struct dictionary 
*create_dictionary_entry(struct dictionary **fst_dictionary_ptr,
                         struct dictionary **lst_dictionary_ptr,
                         char res_name[4],int *ndic_entries)
{
  struct dictionary *dictionary_ptr;
  struct dictionary *first_dictionary_ptr, *last_dictionary_ptr;


  /* Initialise dictionary pointers */
  dictionary_ptr = NULL;
  first_dictionary_ptr = *fst_dictionary_ptr;
  last_dictionary_ptr = *lst_dictionary_ptr;

  /* Create dictionary entry */
  dictionary_ptr =
    (struct dictionary*)malloc(sizeof(struct dictionary));
  if (dictionary_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for");
      printf(" struct dictionary\n");
      printf("***        Program hbadd terminated with");
      printf(" error.\n");
                          exit (1);
    }

  /* Store the details for this entry */
  strcpy(dictionary_ptr->res_name,res_name);
  dictionary_ptr->mol_name = &null;
  dictionary_ptr->natoms = 0;
/* v.5.0--> */
  dictionary_ptr->nbonds = 0;
/* <--v.5.0 */
  dictionary_ptr->nmatch = 0;
  dictionary_ptr->non_hydrogen = 0;
  dictionary_ptr->have_hydrogens = FALSE;
  dictionary_ptr->fake_entry = FALSE;
  dictionary_ptr->nconnect = 0;
  dictionary_ptr->residue_ptr = NULL;
  dictionary_ptr->first_dic_data_ptr = NULL;
  dictionary_ptr->next_dictionary_ptr = NULL;

  /* If this is the first entry stored, then update
     pointer to head of linked list */
  if (first_dictionary_ptr == NULL)
    first_dictionary_ptr = dictionary_ptr;

  /* Otherwise, make previous entry point to the current
     one */
  else
    last_dictionary_ptr->next_dictionary_ptr = dictionary_ptr;
                      
  /* Save the current residue's pointer */
  last_dictionary_ptr = dictionary_ptr;
  (*ndic_entries)++;

  /* Return current pointers */
  *fst_dictionary_ptr = first_dictionary_ptr;
  *lst_dictionary_ptr = last_dictionary_ptr;
  return(dictionary_ptr);
}
/***********************************************************************

create_dic_data_entry  -  Create a dic_data entry

***********************************************************************/

struct dic_data *create_dic_data_entry(struct dic_data **fst_dic_data_ptr,
                                       struct dic_data **lst_dic_data_ptr,
                                       int natoms)
{
  int i;

  struct dic_data *dic_data_ptr;
  struct dic_data *first_dic_data_ptr, *last_dic_data_ptr;


  /* Initialise dic_data pointers */
  dic_data_ptr = NULL;
  first_dic_data_ptr = *fst_dic_data_ptr;
  last_dic_data_ptr = *lst_dic_data_ptr;

  /* Create dic_data entry */
  dic_data_ptr =
    (struct dic_data*)malloc(sizeof(struct dic_data));
  if (dic_data_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for");
      printf(" struct dic_data\n");
      printf("***        Program hbadd terminated with");
      printf(" error.\n");
      exit (1);
    }

  /* Store the data for this entry */
  dic_data_ptr->atom_name[0] = '\0';
  dic_data_ptr->hydrogen = FALSE;
  dic_data_ptr->output_atom_name[0] = '\0';
  dic_data_ptr->number = natoms;
  dic_data_ptr->nbound = 0;

  /* Initialize other variables */
  strcpy(dic_data_ptr->element,"  ");
  for (i = 0; i < MAXCONNECT; i++)
    {
      dic_data_ptr->connect_atom_ptr[i] = NULL;
      dic_data_ptr->boundto[i][0] = '\0';
      dic_data_ptr->boundto_num[i] = 0;
      dic_data_ptr->bond_order[i] = SINGLE;
    }
  
  /* If this is the first entry stored, then update
     pointer to head of linked list */
  if (first_dic_data_ptr == NULL)
    first_dic_data_ptr = dic_data_ptr;

  /* Otherwise, make previous entry point to the current
     one */
  else
    last_dic_data_ptr->next_dic_data_ptr = dic_data_ptr;
                      
  /* Save the current residue's pointer */
  last_dic_data_ptr = dic_data_ptr;

  /* Set pointer to next entry to null */
  dic_data_ptr->next_dic_data_ptr = NULL;

  /* Return current pointers */
  *fst_dic_data_ptr = first_dic_data_ptr;
  *lst_dic_data_ptr = last_dic_data_ptr;
  return(dic_data_ptr);
}
/***********************************************************************

interpret_conect_line  -  Extract the connected atoms from the given
                          CONECT line in the Het Group Dictionary

***********************************************************************/

void interpret_conect_line(char line[LINE_LEN + 1],char res_name[4],
                           struct dic_data *dic_data_ptr,
                           struct dictionary *dictionary_ptr)
{
  char atom_name[5], element[3], nbound_char[3];
  static char ctrl_M = 13;

  int i, j;

  /* Extract the atom name */ 
  strncpy(atom_name,line+11,4);
/* v.4.4.1--> */
  atom_name[4] = '\0';
/* <--v.4.4.1 */

  /* Remove any spurious characters present in the atom name  */
  for (i = 0; i < 4; i++)
    if (atom_name[i] == '\0' || atom_name[i] == '\n')
      atom_name[i] = ' ';

  /* Store the atom name */
  strcpy(dic_data_ptr->atom_name,atom_name);
  strcpy(dic_data_ptr->output_atom_name,atom_name);

  /* Extract element from the atom name */
  get_element(atom_name,element);

  /* Store element */
  strcpy(dic_data_ptr->element,element);

  /* Determine whether this is a hydrogen atom */
  if (!strcmp(dic_data_ptr->element," H"))
    dic_data_ptr->hydrogen = TRUE;

  /* Number of connections */ 
  if ((int) strlen(line) >= 19)
    {
      strncpy(nbound_char, line+18, 2);
      nbound_char[2] = '\0';
      dic_data_ptr->nbound = atoi(nbound_char);

      /* Check that maximum not exceeded */
      if (dic_data_ptr->nbound > MAXCONNECT - 1)
        {
          printf("*** ERROR. Maximum number of atom connections ");
          printf("exceeded: %d\n",MAXCONNECT);
          printf("           Dictionary atom [%s] [%s]\n",
                 res_name,dic_data_ptr->atom_name);
          exit(0);
        }
    }
  else
    dic_data_ptr->nbound = 0;

  /* Connectivity information */
  for (i = 0; i < dic_data_ptr->nbound; i++)
    {
/* v.4.4.1--> */
/*      strncpy(dic_data_ptr->boundto[i], line+(20+(i*5)), 5); */
      strncpy(dic_data_ptr->boundto[i], line+(20+(i*5)),4);
      dic_data_ptr->boundto[i][4] = '\0';
/* <--v.4.4.1 */

      /* Remove carriage returns and other stray characters */
      for (j = 0; j < 4; j++)
        if ((dic_data_ptr->boundto[i][j] == '\n') 
            || (dic_data_ptr->boundto[i][j] == ctrl_M)
            || (dic_data_ptr->boundto[i][j] == '\0'))
          dic_data_ptr->boundto[i][j] = ' ';
      dic_data_ptr->boundto[i][4] = '\0';

      /* Get element of the bonded atom */
      get_element(dic_data_ptr->boundto[i],element);

      /* If neither atom is a hydrogen, then increment count of covalent
         bonds for this dictionary entry */
      if (strcmp(dic_data_ptr->element," H") && strcmp(element," H"))
        {
          /* Increment count only if haven't already done so for this
             atom pair */
          if (strcmp(dic_data_ptr->atom_name,dic_data_ptr->boundto[i]) < 0)
            dictionary_ptr->nconnect++;
        }
    }
}
/***********************************************************************

tidy_line  -  Edit input line to remove troublesome quotes

***********************************************************************/

void tidy_line(char *input_line)
{
  char *string_ptr;

  int iend, ipos, istart, len;
  int done;

  /* Initialise variables */
  done = FALSE;

  /* Get the length of the line */
  len = strlen(input_line);

  /* Loop until all quotes have been dealt with */
  while (done == FALSE)
    {
      /* Check for start of quotes */
      string_ptr = strstr(input_line," '");

      /* If have a start-quote, then continue */
      if (string_ptr != NULL)
        {
          /* Get the first character position of the string in quotes */
          istart = (int) (string_ptr - input_line) + 2;

          /* Check for end-quotes */
          string_ptr = strstr(input_line,"' ");

          /* Knock out the opening-quotes */
          input_line[istart - 1] = '.';

          /* If have end-quote, then examine string within the quotes */
          if (string_ptr != NULL)
            {
              /* Get the first character position of the string in quotes */
              iend = (int) (string_ptr - input_line) - 1;

              /* Loop through the characters enclosed in quotes */
              for (ipos = istart; ipos < iend + 1; ipos++)
                {
                  /* If character is a space, replace by ampersand */
                  if (input_line[ipos] == ' ')
                    input_line[ipos] = '&';
                }
              
              /* Knock out the end-quotes */
              input_line[iend + 1] = '.';
            }

          /* Otherwise, show warning message and finish */
          else
            {
              /* Show warning */
              printf("*** Warning. Start-quotes but no end quotes in "
                     "CIF file ");
              printf("[%s]\n",input_line);
              done = TRUE;
            }
        }

      /* Otherwise, no more processing necessary */
      else
        done = TRUE;
    }

  /* Replace all dots by spaces */
  for (ipos = 0; ipos < len; ipos++)
    if (input_line[ipos] == '.')
      input_line[ipos] = ' ';
}
/***********************************************************************

get_atom_name  -  Extract the atom name from the given token

***********************************************************************/

int get_atom_name(char token[TOKEN_LEN],char atom_name[5])
{
/* v.4.4.2--> */
//  char quote, tmp[TOKEN_LEN];
  char dquote, quote, tmp[TOKEN_LEN];
/* <--v.4.4.2 */
  static char q_string[] = { "'" };
/* v.4.4.2--> */
  static char dq_string[] = { "\"" };
/* <--v.4.4.2 */

  int ipos, len;
  int valid;

  /* Initialise variables */
  quote = q_string[0];
/* v.4.4.2--> */
  dquote = dq_string[0];
/* <--v.4.4.2 */
  valid = TRUE;

  /* Get the length of the token */
  strcpy(tmp,token);
  len = strlen(tmp);

  /* Check for opening quotes in the first position */
  if (token[0] == quote)
    {
      /* Remove the start- and end-tokens */
      token[len - 1] = '\0';
      strcpy(tmp,token+1);

      /* Get new length of token */
      len = strlen(tmp);
    }

/* v.4.4.2--> */
  /* Repeat for double quotes */
  if (token[0] == dquote)
    {
      /* Remove the start- and end-tokens */
      token[len - 1] = '\0';
      strcpy(tmp,token+1);

      /* Get new length of token */
      len = strlen(tmp);
    }
/* <--v.4.4.2 */

  /* If token not too long, copy across to atom name */
  if (len < 5)
    {
      /* Copy name across */
      strcpy(atom_name,tmp);

      /* Pad out atom-name to 4 characters in length */
      strcpy(tmp,atom_name);
      len = strlen(tmp);
      for (ipos = len; ipos < 4; ipos++)
        tmp[ipos] = ' ';
      tmp[4] = '\0';
      strcpy(atom_name,tmp);

      /* Remove any ampersands */
      for (ipos = 0; ipos < 4; ipos++)
        if (atom_name[ipos] == '&')
          atom_name[ipos] = ' ';
    }

  /* Otherwise, atom-name too long! */
  else
    {
      /* Show warning message */
      printf("*** Warning. Atom name too long! [%s] ",tmp);
      valid = FALSE;
    }

  /* Return whether atom name is valid or not */
  return(valid);
}
/***********************************************************************

cif_bond_line  -  Interpret the current CIF bond line and determine
                  whether it holds valid bond information

***********************************************************************/

int cif_bond_line(char token[4][TOKEN_LEN],struct dic_data **dic_d1_ptr,
                  struct dic_data **dic_d2_ptr,int *bond_order,
                  struct dictionary *dictionary_ptr)
{
  char atom_name[5], tmp[TOKEN_LEN];

  int iatom, itoken;

  int found, valid;

  struct dic_data *dic_data_ptr, *save_dic_data_ptr[2];

  /* Initialise dic_data pointers */
  save_dic_data_ptr[0] = *dic_d1_ptr;
  save_dic_data_ptr[1] = *dic_d2_ptr;

  /* Set flag */
  valid = TRUE;

  /* Process the two atom names */
  for (iatom = 0; iatom < 2 && valid == TRUE; iatom++)
    {
      /* Get index of the appropriate token */
      itoken = iatom + 1;

      /* Get the atom name from the token */
      strcpy(tmp,token[itoken]);

      /* Extract the atom name */
      valid = get_atom_name(tmp,atom_name);

      /* If atom name is valid, locate it in the list of atom records
         stored */
      if (valid == TRUE)
        {
          /* Get pointer to the first of this dictionary entry's 
             dic_data records */
          dic_data_ptr = dictionary_ptr->first_dic_data_ptr;
          found = FALSE;

          /* Loop through the entries to find the one corresponding to
             our atom */
          while (dic_data_ptr != NULL && found == FALSE)
            {
              /* If atom names agree, then have a match */
              if (!strcmp(atom_name,dic_data_ptr->atom_name) ||
                  (atom_name[3] == ' ' && 
                   !strncmp(atom_name,dic_data_ptr->atom_name+1,3) &&
                   dic_data_ptr->atom_name[0] == ' '))
                {
                  /* Set flag */
                  found = TRUE;

                  /* Save the pointer to this record */
                  save_dic_data_ptr[iatom] = dic_data_ptr;
                }

              /* Get next entry */
              dic_data_ptr = dic_data_ptr->next_dic_data_ptr;
            }

          /* If entry not found, then show warning */
          if (found == FALSE)
            {
              /* Show warning message */
              printf("*** Warning. Atom type [%s] not defined in ",
                     atom_name);
              printf("CIF file for residue entry [%s]\n",
                     dictionary_ptr->res_name);

              /* Set flag to denote invalid entry */
              valid = FALSE;
            }
        }
    }

  /* Get the bond order */
  if (!strncmp(token[3],"SING",4))
    *bond_order = SINGLE;
  else if (!strncmp(token[3],"DOUB",4))
    *bond_order = DOUBLE;
  else if (!strncmp(token[3],"TRIP",4))
    *bond_order = TRIPLE;
  else if (!strncmp(token[3],"AROM",4))
    *bond_order = AROMATIC;
  else
    {
      printf("*** Warning. Unrecognised bond-order in dictionary: [%s]\n",
             token[3]);
      *bond_order = SINGLE;
    }

  /* Return the two dic_data pointers */
  *dic_d1_ptr = save_dic_data_ptr[0];
  *dic_d2_ptr = save_dic_data_ptr[1];

  /* Return whether this line is valid or not */
  return(valid);
}
/***********************************************************************

cif_atom_line  -  Interpret the current CIF line and determine
                       whether it holds valid bond information

***********************************************************************/

int cif_atom_line(char token[4][TOKEN_LEN],char atom_name[5],
                  struct dic_data **fst_dic_data_ptr,
                  struct dic_data **lst_dic_data_ptr,
/* v.4.4.2--> */
//                  struct dictionary *dictionary_ptr)
                  struct dictionary *dictionary_ptr,int new_format)
/* <--v.4.4.2 */
{
  char element[5];
  char tmp[TOKEN_LEN];

/* v.4.4.2--> */
//  int len;
  int element_field, len;
/* <--v.4.4.2 */
  int shift, valid;

  struct dic_data *dic_data_ptr;
  struct dic_data *first_dic_data_ptr, *last_dic_data_ptr;

  /* Initialise dic_data pointers */
  first_dic_data_ptr = *fst_dic_data_ptr;
  last_dic_data_ptr = *lst_dic_data_ptr;

  /* Initialise variables */
  valid = TRUE;
/* v.4.4.2--> */
  if (new_format == TRUE)
    element_field = 3;
  else
    element_field = 2;
/* <--v.4.4.2 */

  /* Get the length of the element item */
/* v.4.4.2--> */
//  len = strlen(token[2]);
  len = strlen(token[element_field]);
/* <--v.4.4.2 */

  /* Check that it is valid */
  if (len < 3)
    {
      /* If single character, then shift to second character position */
      if (len == 1)
        {
          element[0] = ' ';
/* v.4.4.2--> */
//          element[1] = token[2][0];
          element[1] = token[element_field][0];
/* <--v.4.4.2 */
          element[2] = '\0';
        }

      /* Otherwise, copy straight across */
      else
/* v.4.4.2--> */
//        strcpy(element,token[2]);
        strcpy(element,token[element_field]);
/* <--v.4.4.2 */
    }

  /* Otherwise, show warning */
  else
    {
      /* Show warning message */
      printf("*** Warning. Invalid element identifier in CIF file");
/* v.4.4.2--> */
//      printf(" [%s]\n",token[2]);
      printf(" [%s]\n",token[element_field]);
/* <--v.4.4.2 */
      return(FALSE);
    }

  /* Get the atom name */
  valid = get_atom_name(token[1],atom_name);

  /* If name valid, then format it correctly */
  if (valid == TRUE)
    {
      /* Get the current length of the atom name */
      len = strlen(atom_name);

      /* Use the element name to determine the correct format of the 
         atom name */
      shift = FALSE;

      /* If have a single-character element name */
      if (element[0] == ' ')
        {
          /* If element matches atom's second character, then no shift
             required */
          if (element[1] == atom_name[1])
            shift = FALSE;

          /* Otherwise, must match at first character, and need to
             shift name */
          else if (element[1] == atom_name[0])
            {
               if (atom_name[3] == ' ')
                 shift = TRUE;
               else
                 shift = FALSE;
             }

          /* For hydrogens, the H might be in the third character position */
          else if (!strncmp(element," H",2) && atom_name[2] == 'H')
            shift = FALSE;

          /* Otherwise, have a problem */
          else
            {
              /* Show warning message */
              printf("*** Warning. Can't match element identifier [%s]",
                     element);
              printf(" to atom name [%s] in CIF file\n",atom_name);
              return(FALSE);
            }
        }

      /* If have a two-character element name (ie metal), then must
         match first two characters of the atom name */
      else if (!strncmp(element,atom_name,2))
        shift = FALSE;

      /* Except, of course, in special circumstances when the match
         is in the third and fourth character positions (!) */
      else if (!strncmp(element,atom_name+1,2))
        shift = FALSE;

      /* Otherwise, have another problem */
      else
        {
          /* Show warning message */
          printf("*** Warning. Can't match element identifier [%s] to ",
                 element);
          printf("atom name [%s] in CIF file\n",atom_name);
          valid = FALSE;
        }

      /* If shift required, then make it */
      if (shift == TRUE)
        {
          strcpy(tmp,atom_name);
          len = strlen(tmp);
          atom_name[0] = ' ';
          strcpy(atom_name+1,tmp);
          atom_name[4] = '\0';
        }
    }

  /* If have a valid atom name, then create a new dic-data entry for
     it */
  if (valid == TRUE)
    {
      /* Create new entry */
      dic_data_ptr
        = create_dic_data_entry(&first_dic_data_ptr,
                                &last_dic_data_ptr,
                                dictionary_ptr->natoms);

      /* Store the atom name */
      strcpy(dic_data_ptr->atom_name,atom_name);
      strcpy(dic_data_ptr->output_atom_name,atom_name);

      /* Store the element */
      strcpy(dic_data_ptr->element,element);

      /* Determine whether this is a hydrogen atom */
      if (!strncmp(element," H",2))
        dic_data_ptr->hydrogen = TRUE;

      /* If this is the first dic_data record for this dictionary
         entry, then store pointer to it */
      if (dictionary_ptr->natoms == 0)
        dictionary_ptr->first_dic_data_ptr = dic_data_ptr;

      /* Increment atom-count for this dictionary entry */
      dictionary_ptr->natoms++;
      if (dic_data_ptr->hydrogen == FALSE)
        dictionary_ptr->non_hydrogen++;
    }

  /* Return current pointers */
  *fst_dic_data_ptr = first_dic_data_ptr;
  *lst_dic_data_ptr = last_dic_data_ptr;

  /* Return whether this line is valid or not */
  return(valid);
}
/***********************************************************************

interpret_cif_line  -  Interpret the current CIF line and determine
                       whether it holds valid atom or bond information

***********************************************************************/

int interpret_cif_line(char *line,char atom_names[2][5],
                       int *bnd_ord,struct dic_data **fst_dic_data_ptr,
                       struct dic_data **lst_dic_data_ptr,
/* v.4.4.2--> */
//                       struct dictionary *dictionary_ptr)
                       struct dictionary *dictionary_ptr,int new_format)
/* <--v.4.4.2 */
{
  char atom_name[5];
  char token[4][TOKEN_LEN];

  int bond_order;
  int bond_line, valid;

  struct dic_data *dic_data1_ptr, *dic_data2_ptr;
  struct dic_data *first_dic_data_ptr, *last_dic_data_ptr;

  /* Initialise dic_data pointers */
  first_dic_data_ptr = *fst_dic_data_ptr;
  last_dic_data_ptr = *lst_dic_data_ptr;

  /* Initialise variables */
  bond_line = FALSE;
  bond_order = *bnd_ord;
  valid = FALSE;

  /* Tidy up the line to get rid of any nastiness */
  tidy_line(line);

  /* Read in the four tokens on the current line */
  sscanf(line,"%s%s%s%s",token[0],token[1],token[2],token[3]);

  /* If fourth token indicates the bond order, then this is a valid
     bond record */
  if (!strncmp(token[3],"SING",4) || !strncmp(token[3],"DOUB",4) || 
      !strncmp(token[3],"TRIP",4) || !strncmp(token[3],"AROM",4) ||
      token[3][0] == '?')
    {
      /* Extract the bond information */
      valid = cif_bond_line(token,&dic_data1_ptr,&dic_data2_ptr,
                            &bond_order,dictionary_ptr);

      /* If a valid bond-line, then flag it as such */
      if (valid == TRUE)
        {
          /* Set flag */
          bond_line = TRUE;

          /* Store the two atom names */
          strcpy(atom_names[0],dic_data1_ptr->atom_name);
          strcpy(atom_names[1],dic_data2_ptr->atom_name);
        }
    }

  /* Otherwise, take it to be an atom-defining record */
  else
    {
      /* Extract the atom information */
      valid = cif_atom_line(token,atom_name,&first_dic_data_ptr,
/* v.4.4.2--> */
//                            &last_dic_data_ptr,dictionary_ptr);
                            &last_dic_data_ptr,dictionary_ptr,new_format);
/* <--v.4.4.2 */
    }

  /* Return bond order */
  *bnd_ord = bond_order;

  /* Return current pointers */
  *fst_dic_data_ptr = first_dic_data_ptr;
  *lst_dic_data_ptr = last_dic_data_ptr;

  /* Return whether this line is valid or not */
  return(bond_line);
}
/***********************************************************************

store_bond  -  Store the given bond as two dic_data entries

***********************************************************************/

void store_bond(struct dic_data **fst_dic_data_ptr,
                struct dic_data **lst_dic_data_ptr,
                struct dictionary *dictionary_ptr,char atom_names[2][5],
                int bond_order)
{
  char atom_name[5], other_atom_name[5];

  int iatom, nbound;
  int found;

  struct dic_data *dic_data_ptr, *save_dic_data_ptr;
  struct dic_data *first_dic_data_ptr, *last_dic_data_ptr;
  struct dic_data *store_dic_data_ptr[2];

  /* Initialise dic_data pointers */
  first_dic_data_ptr = *fst_dic_data_ptr;
  last_dic_data_ptr = *lst_dic_data_ptr;

  /* Loop through the two atoms */
  for (iatom = 0; iatom < 2; iatom++)
    {
      /* Get the atom names the way we want */
      strcpy(atom_name,atom_names[iatom]);
      strcpy(other_atom_name,atom_names[1 - iatom]);

      /* Get pointer to the first of this dictionary entry's dic_data
         records */
      dic_data_ptr = dictionary_ptr->first_dic_data_ptr;
      save_dic_data_ptr = NULL;
      found = FALSE;

      /* Loop through the entries to find the one corresponding to
         our atom */
      while (dic_data_ptr != NULL && found == FALSE)
        {
          /* If atom names agree, then have a match */
          if (!strcmp(atom_name,dic_data_ptr->atom_name) ||
              (atom_name[3] == ' ' && 
               !strncmp(atom_name,dic_data_ptr->atom_name+1,3)))
            {
              /* Set flag */
              found = TRUE;

              /* Save the pointer to this record */
              save_dic_data_ptr = dic_data_ptr;
            }

          /* Get next entry */
          dic_data_ptr = dic_data_ptr->next_dic_data_ptr;
        }

      /* If no match found, then need to create a dic_data entry for
         this atom */
      if (found == FALSE)
        {
          /* Create an entry for this atom's connectivities */
          dic_data_ptr
            = create_dic_data_entry(&first_dic_data_ptr,
                                    &last_dic_data_ptr,
                                    dictionary_ptr->natoms);

          /* Store the atom name */
          strcpy(dic_data_ptr->atom_name,atom_name);
          strcpy(dic_data_ptr->output_atom_name,atom_name);

          /* Determine whether this is a hydrogen atom */
          dic_data_ptr->hydrogen = is_hydrogen_atom(atom_name);

          /* If this is the first dic_data record for this dictionary
             entry, then store pointer to it */
          if (dictionary_ptr->natoms == 0)
            dictionary_ptr->first_dic_data_ptr = dic_data_ptr;

          /* Increment atom-count for this dictionary entry */
          dictionary_ptr->natoms++;
          if (dic_data_ptr->hydrogen == FALSE)
            dictionary_ptr->non_hydrogen++;
        }

      /* Otherwise, retrieve the pointer to the entry */
      else
        dic_data_ptr = save_dic_data_ptr;

      /* Store this pointer */
      store_dic_data_ptr[iatom] = dic_data_ptr;

      /* Store the connected atom */
      nbound = dic_data_ptr->nbound;
      if (nbound < MAXCONNECT)
        {
          /* Copy the connected atom across */
          strcpy(dic_data_ptr->boundto[nbound],other_atom_name);

          /* Store the bond order */
          dic_data_ptr->bond_order[nbound] = bond_order;

          /* Increment count of bonds */
          dic_data_ptr->nbound++;
        }

      /* Otherwise, show error message */
      else
        {
          printf("*** ERROR. Maximum number of atom connections ");
          printf("exceeded: %d\n",MAXCONNECT);
          printf("           Dictionary atom [%s] [%s]\n",
                 dictionary_ptr->res_name,dic_data_ptr->atom_name);
          exit(0);
        }
    }

  /* If neither atom is a hydrogen, then increment count of covalent
     bonds for this dictionary entry */
  if (store_dic_data_ptr[0]->hydrogen == FALSE &&
      store_dic_data_ptr[1]->hydrogen == FALSE)
    dictionary_ptr->nconnect++;

  /* Return current pointers */
  *fst_dic_data_ptr = first_dic_data_ptr;
  *lst_dic_data_ptr = last_dic_data_ptr;
}
/* v.5.1--> */
/***********************************************************************

add_to_dictionary  -  Add the current residue definition to the
                      dictionary and create link to all relevant PDB
                      residues

***********************************************************************/

struct dictionary *add_to_dictionary(char *res_name,char *dic_res_name,
                                     struct residue *first_residue_ptr,
                                     struct dictionary **fst_dictionary_ptr,
                                     struct dictionary **lst_dictionary_ptr,
                                     int *ndic,int *wntd)
{
  int ndic_entries, nmatch;
  int have_entry, wanted;

  struct dictionary *dictionary_ptr;
  struct dictionary *first_dictionary_ptr, *last_dictionary_ptr;
  struct residue *residue_ptr;

  /* Initialise variables */
  have_entry = FALSE;
  ndic_entries = *ndic;

  /* Initialise pointers */
  dictionary_ptr = NULL;
  first_dictionary_ptr = *fst_dictionary_ptr;
  last_dictionary_ptr = *lst_dictionary_ptr;

  /* If looking for a specific dictionary entry, see if this
     is the one */
  if (!strcmp(res_name,dic_res_name))
    have_entry = TRUE;

  /* Initialise variables for search through stored residues */
  residue_ptr = first_residue_ptr;
  nmatch = 0;
  wanted = FALSE;

  /* If looking for a specific dictionary entry and this isn't
     the one, don't need to process all the residues */
  if (dic_res_name[0] != '\0' && have_entry == FALSE)
    residue_ptr = NULL;

  /* Loop through all HET groups read in from the PDB file to
     see if this is amongst them */
  while (residue_ptr != NULL)
    {
      /* If the residue names match, then want to store this
         dictionary entry */
      if (!strncmp(res_name,residue_ptr->res_name,3) ||
          (have_entry == TRUE && residue_ptr->dictionary_ptr == NULL))
        {
          wanted = TRUE;

          /* If this is the first match for this het group,
             create an entry in the dictionary linked list */
          if (nmatch == 0 ||
              (dic_res_name[0] != '\0' && residue_ptr->natoms > 3))
            {
              /* Create a dictionary entry */
              dictionary_ptr
                = create_dictionary_entry(&first_dictionary_ptr,
                                          &last_dictionary_ptr,
                                          res_name,&ndic_entries);
            }

          /* Add link from residue to dictionary entry */
          residue_ptr->dictionary_ptr = dictionary_ptr;

          /* Increment count of matches and store */
          nmatch++;
          residue_ptr->nmatch = nmatch;
          if (dic_res_name[0] != '\0')
            residue_ptr->nmatch = 1;
        }

      /* Get pointer to the next residue in the list */
      residue_ptr = residue_ptr->next_residue_ptr;
    }

  /* Show number of matches found */
  if (wanted == TRUE)
    {
      printf("      Found match for: %s ",res_name);
      if (nmatch > 1)
        printf("x %d\n",nmatch);
      else
        printf("\n");

      /* If looking for a specific dictionary entry, work with the
         first dictionary entry */
      if (dic_res_name[0] != '\0')
        dictionary_ptr = first_dictionary_ptr;
    }

  /* Return current count of dictionary entries and whether entry
     wanted */
  *ndic = ndic_entries;
  *wntd = wanted;

  /* Return current pointers */
  *fst_dictionary_ptr = first_dictionary_ptr;
  *lst_dictionary_ptr = last_dictionary_ptr;

  return(dictionary_ptr);
}
/* <--v.5.1 */
/***********************************************************************

read_dictionary  -  Read through HET Group Dictionary to pick up all the
                    required HET group definitions

***********************************************************************/

int read_dictionary(char dictionary_name[],
                    struct residue *first_residue_ptr,char *dic_res_name,
                    struct dictionary **fst_dictionary_ptr,
                    struct dictionary **lst_dictionary_ptr,
                    struct dic_data **fst_dic_data_ptr,
                    struct dic_data **lst_dic_data_ptr)
{
  char atom_names[2][5];
  char res_name[RESNAM_LEN], match_res_name[4];
  char line[LINE_LEN + 1], tmp_line[LINE_LEN + 1];
/* v.5.4.3--> */
  char element[3];
/* <--v.5.4.3 */

  char *string_ptr;

  int bond_order, len, match_len, ndic_entries;
  int bond_line, dic_type, residue_start, wanted;
/* v.5.4.3--> */
  int i, start_pos;
  int atom_next, bond_next;
/* <--v.5.4.3 */
/* v.4.4.2--> */
//  int done, have_entry, name_next;
  int done, have_entry, name_next, new_format;
/* <--v.4.4.2 */

  struct dictionary *dictionary_ptr;
  struct dictionary *first_dictionary_ptr, *last_dictionary_ptr;
  struct dic_data *dic_data_ptr, *first_dic_data_ptr, *last_dic_data_ptr;

  /* Initialise variables */
/* v.5.4.3--> */
  atom_next = bond_next = FALSE;
/* <--v.5.4.3 */
  bond_order = SINGLE;
  dic_type = TEXT;
  dictionary_ptr = NULL;
  done = FALSE;
  have_entry = FALSE;
  first_dictionary_ptr = NULL;
  last_dictionary_ptr = NULL;
  dic_data_ptr = NULL;
  first_dic_data_ptr = NULL;
  last_dic_data_ptr = NULL;
  match_len = 0;
  match_res_name[0] = '\0';
  name_next = FALSE;
  ndic_entries = 0;
/* v.4.4.2--> */
  new_format = FALSE;
/* <--v.4.4.2 */
  residue_start = FALSE;
  res_name[0] = '\0';
  wanted = FALSE;
  
  /* Open the HET Group Dictionary */
  printf("  Opening HET Group Dictionary ...\n");
  if ((fil_dic = fopen(dictionary_name,"r")) == NULL)
    {
      printf("\n*** Unable to open HET Group Dictionary [%s]\n",
             dictionary_name);
      return(ndic_entries);
    }
    
  /* Prepare to start reading in the data from the file */
  printf("\n");
  printf("    Scanning HET Group Dictionary ...\n");

  /* Read through the dictionary */
  while (fgets(line,LINE_LEN,fil_dic) != NULL && done == FALSE)
    {
      /* Truncate the line */
      len = string_truncate(line,LINE_LEN);

      /* If this is the start of a new residue, then check whether
         it is one of the ones required */
      if (!strncmp(line,"RESIDUE",7))
        {
          /* Extract residue name */ 
          strncpy(res_name,line+10,3);
          res_name[3] = '\0';

          /* Set flag to indicate that we have the start of a residue's
             details */
          residue_start = TRUE;
          dic_type = TEXT;
        }

      /* Check for CIF-format residue-start */
      else if (!strncmp(line,"data_",5) && strncmp(line+5,"UNK",3))
        {
          /* Initialise residue name */
          res_name[0] = '\0';

          /* Extract the residue name */
          strcpy(tmp_line,line+5);

          /* Check for terminating space */
          string_ptr = strstr(tmp_line," ");

          /* If have a space, then residue name is the KEGG identifier,
             so extract */
          if (string_ptr != NULL)
            {
              /* Copy KEGG identifier */
              strcpy(res_name,string_ptr+1);

              /* Store the version of the residue name that will need to
                 be matched in the atom and connectivity records to follow */
              strncpy(match_res_name,tmp_line,3);
              match_res_name[3] = '\0';
              match_len = 3;
            }

          /* Otherwise, have standard residue name */
          else
            {
              /* Replace underscore(s) by space(s) */
              if (tmp_line[0] == '_')
                tmp_line[0] = ' ';
              if (tmp_line[1] == '_')
                tmp_line[1] = ' ';

              /* Check length of name, and extract accordingly */
              len = strlen(tmp_line);
              if (len == 1)
                {
                  strcat(res_name,"  ");
                  strncat(res_name,tmp_line,1);
                }
              else if (len == 2)
                {
                  strcat(res_name," ");
                  strncat(res_name,tmp_line,2);
                }
              else
                {
                  strncpy(res_name,tmp_line,3);
                }
              res_name[3] = '\0';


              /* Store the version of the residue name that will need to
                 be matched in the atom and connectivity records to follow */
              strcpy(match_res_name,tmp_line);
              match_len = len;
            }

          /* Set flag to indicate that we have the start of a residue's
             details */
          residue_start = TRUE;
          dic_type = CIF;
          name_next = FALSE;
/* v.5.4.3--> */
          atom_next = bond_next = FALSE;
/* <--v.5.4.3 */
        }

/* v.5.4.3--> */
      /* If hit loop or hash and reading atoms or bonds, unset flag */
      else if ((line[0] == '#' || !strncmp(line,"loop_",5)) &&
               (atom_next == TRUE || bond_next == TRUE))
        atom_next = bond_next = FALSE;

      /* Check if this is the start of the atoms list */
      else if (!strncmp(line,"_chem_comp_atom.comp_id",23) &&
               dictionary_ptr != NULL)
        {
          /* Sat flag */
          atom_next = TRUE;

          /* See if this is a single atom, with the name on this line */
          if (line[23] == ' ')
            {
              /* Initialise start position */
              start_pos = 0;

              /* Find start of name */
              for (i = 23; i < len && start_pos == 0; i++)
                {
                  /* If not a space, have the start of the name */
                  if (line[i] != ' ')
                    start_pos = i;
                }
              
              /* If name found, create a dummy version of this line */
              if (start_pos > 0)
                {
                  /* Assume this is a metal, so form the element name */
                  strcpy(tmp_line,line + start_pos);
                  len = strlen(tmp_line);
                  if (len > 2)
                    tmp_line[2] = '\0';
                  strcpy(element,tmp_line);

                  /* Form dummy line */
                  sprintf(line,"%s %s %s %s",match_res_name,
                          line + start_pos,line + start_pos,element);

                  /* Interpret the current CIF line */
                  interpret_cif_line(line,atom_names,&bond_order,
                                     &first_dic_data_ptr,&last_dic_data_ptr,
                                     dictionary_ptr,new_format);

                  /* Unset the flag */
                  atom_next = FALSE;
                }
            }
        }

      /* Check if this is the start of the bonds list */
      else if (!strcmp(line,"_chem_comp_bond.comp_id") &&
               dictionary_ptr != NULL)
        bond_next = TRUE;
/* <--v.5.4.3 */

      /* Check for CIF-format molecule name tag */
      else if (!strncmp(line,"_chem_comp.name",15))
        {
          /* Set flag */
          name_next = TRUE;

/* v.5.4.3--> */
          /* Check if the name is on this line */
          if (line[15] == ' ')
            {
              /* Initialise start position */
              start_pos = 0;

              /* Find start of name */
              for (i = 15; i < len && start_pos == 0; i++)
                {
                  /* If not a space, have the start of the name */
                  if (line[i] != ' ')
                    start_pos = i;
                }
              
              /* If name found, then store it and unset flag */
              if (start_pos > 0)
                {
                  /* Store the name */
                  if (dictionary_ptr != NULL)
                    store_string(&(dictionary_ptr->mol_name),
                                 line + start_pos);

                  /* Unset the flag */
                  name_next = FALSE;
                }
            }
/* <--v.5.4.3 */
        }

/* v.4.4.2--> */
      /* Check for new-format CIF file */
      else if (!strncmp(line,"_chem_comp_atom.alt_atom_id",27))
        {
          /* Set flag */
          new_format = TRUE;
        }
/* <--v.4.4.2 */

      /* Check for CIF-format molecule name */
      else if (line[0] == ';' && name_next == TRUE)
        {
          /* Store the molecule name */
          if (dictionary_ptr != NULL)
            store_string(&(dictionary_ptr->mol_name),line+1);

          /* Unset flag */
          name_next = FALSE;
        }

      /* Unset molecule name flag if still set */
      else if (name_next == TRUE)
        name_next = FALSE;

      /* If this is the start of a new entry and we already have the
         one we were after, can stop here */
      if (residue_start == TRUE && have_entry == TRUE)
        done = TRUE;

/* v.4.4.2--> */
      /* If have reached the _pdbx lines then don't want the rest of
         this entry */
/* v.5.4.3--> 
      else if (!strncmp(line,"_pdbx",5))
        wanted = FALSE;
<--v.5.4.3 */
/* <--v.4.4.2 */

      /* If have a new residue, then check if it's one we want */
      else if (residue_start == TRUE)
        {
          /* Unset the flag */
          residue_start = FALSE;

/* v.5.1--> */
          /* Add this residue to the dictionary, linking to its entry
             from all PDB residues for this residue type */
          dictionary_ptr
            = add_to_dictionary(res_name,dic_res_name,first_residue_ptr,
                                &first_dictionary_ptr,
                                &last_dictionary_ptr,&ndic_entries,
                                &wanted);
/* <--v.5.1 */
        }

      /* If this is a CONECT record for a wanted residue, then store
         it */
      else if (wanted == TRUE && dic_type == TEXT &&
               !strncmp(line,"CONECT",6))
        {
          /* Create an entry for this atom's connectivities */
          dic_data_ptr = create_dic_data_entry(&first_dic_data_ptr,
                                               &last_dic_data_ptr,
                                               dictionary_ptr->natoms);

          /* Interpret the connectivities on the given line */
          interpret_conect_line(line,res_name,dic_data_ptr,
                                dictionary_ptr);

          /* If this is the first CONECT record for this entry, then
             store pointer to the first atom */
          if (dictionary_ptr->natoms == 0)
            dictionary_ptr->first_dic_data_ptr = dic_data_ptr;

          /* Increment atom-count for this dictionary entry */
          dictionary_ptr->natoms++;
          if (dic_data_ptr->hydrogen == FALSE)
            dictionary_ptr->non_hydrogen++;
        }

      /* If this is a CIF-style connectivity record for a wanted residue,
         then store the bond details */
      else if (wanted == TRUE && dic_type == CIF &&
/* v.5.4.3--> */
//               !strncmp(line,match_res_name,match_len))
               !strncmp(line,match_res_name,match_len) &&
               (atom_next == TRUE || bond_next == TRUE))
/* <--v.5.4.3 */
        {
          /* Interpret the current CIF line */
          bond_line
            = interpret_cif_line(line,atom_names,&bond_order,
                                 &first_dic_data_ptr,&last_dic_data_ptr,
/* v.4.4.2--> */
//                                 dictionary_ptr);
                                 dictionary_ptr,new_format);
/* <--v.4.4.2 */

          /* Process only if this is a valid bond line */
          if (bond_line == TRUE)
            {
              /* Store this bond in the appropriate dic_data entries */
              store_bond(&first_dic_data_ptr,&last_dic_data_ptr,
                         dictionary_ptr,atom_names,bond_order);
            }
        }
    }

  /* Print counts of data read in */
  printf("\n");
  printf("      Number of dictionary entries read in  = %7d\n",ndic_entries);

  /* Close the dictionary file */
  fclose(fil_dic);

  /* Return current pointers */
  *fst_dictionary_ptr = first_dictionary_ptr;
  *lst_dictionary_ptr = last_dictionary_ptr;
  *fst_dic_data_ptr = first_dic_data_ptr;
  *lst_dic_data_ptr = last_dic_data_ptr;

  /* Return number of dictionary entries stored */
  return(ndic_entries);
}
/* <--v.4.3 */
/* v.5.1--> */
/***********************************************************************

get_tokens - Get the space- or tab-delimited token from the given
             input line

**********************************************************************/  

int get_tokens(char *line,char **token,int max_token,char token_sep)
{
  char ch, token_sep2;

/* v.5.1.2 --> */
  // char *tmp, *string;
  char *tmp, *tmp_cpy;
/* <-- v.5.1.2 */

  int intok, ipos, len, ntoken, tokenlen;
  int done, keep_token;

  /* Initialise variables */
  done = FALSE;
  keep_token = FALSE;
  ntoken = 0;
  intok = 0;
  ipos = 0;
  if (token_sep == '[')
    {
      token_sep2 = ']';
      keep_token = TRUE;
    }
  else
    token_sep2 = '\0';

  /* Get the length of the line */
  len = strlen(line);

  /* Allocate space for temporary string */
  tmp = (char *) malloc(sizeof(char)*(len + 1));
  tmp[0] = '\0';
/* v.5.1.2 --> */
  tmp_cpy = (char *) malloc(sizeof(char)*(len + 1));
/* <-- v.5.1.2 */

  /* Loop until end of line reached */
  while (done == FALSE && ipos < len + 1)
    {
      /* Get the next character in the string */
      ch = line[ipos];

      /* Check for end of line */
      if (ch == '\0' || ch == '\n')
        {
          done = TRUE;
          ch = token_sep;
        }

      /* If have a token separator, then take to be the end of the
         token if have already got something in it */
      if (ch == token_sep || ch == token_sep2)
        {
          /* If have a token, increment token count */
          if (tmp[0] != '\0' || token_sep == ',')
            {
              /* If keeping tokens and this is a token-end, add it
                 to the current string */
              if (keep_token == TRUE && ch == token_sep2)
                {
                  tmp[intok] = ch;
                  intok++;
                }

              /* Terminate current token */
              tmp[intok] = '\0';

              /* Check for bounding quotes */
              if (tmp[0] == '"' &&
                  tmp[intok - 1] == '"')
                {
                  /* Strip out the quotes */
/* v.5.1.2 --> */
                  // strcpy(tmp,tmp+1);
                  strcpy(tmp_cpy,tmp+1);
                  strcpy(tmp,tmp_cpy);
/* <-- v.5.1.2 */
                  tmp[intok - 2] = '\0';
                }

              /* Save the current token */
              tokenlen = strlen(tmp);
              store_string(&(token[ntoken]),tmp);

              /* Increment token count */
              ntoken++;

              /* Reinitialise position within current token */
              intok = 0;
              tmp[intok] = '\0';

              /* If keeping tokens and this is a token-start, add it
                 at the start of the current string */
              if (keep_token == TRUE && ch == token_sep)
                {
                  tmp[intok] = ch;
                  intok++;
                }
            }

          /* If keeping separator and this is a start-separator,
             start ne token with it */
          else if (keep_token == TRUE && ch == token_sep)
            {
              tmp[intok] = ch;
              intok++;
            }
        }

      /* Otherwise, add character to token-string provided not
         too long */
      else
        {
          if (ntoken < max_token && intok < len)
            {
              tmp[intok] = ch;
              intok++;
            }
        }

      /* Go to the next character in the line */
      ipos++;
    }

  /* Free up memory taken by temporary line */
  free(tmp);

  /* Return the number of tokens picken up */
  return(ntoken);
} 
/***********************************************************************

to_upper  -  Convert given character string to upper case

***********************************************************************/

void to_upper(char *search_string,int length)
{
  int ichar;

  int ipos;

  /* Loop through all the character positions, converting any letters
     to upper case */
  for (ipos = 0; ipos < length; ipos++)
    {
      /* Check the current character */
      ichar = search_string[ipos] - 'a';
      if (ichar >= 0 && ichar < 26)
        {
          ichar = ichar + 'A';
          search_string[ipos] = ichar;
        }
    }
}
/***********************************************************************

get_attype  -  Identify this atom type and increment count

***********************************************************************/

int get_attype(char *atom_name,char stored_name[MAX_TYPES][5],
               int *name_count,int *ntyp)
{
  int atom_type, itype, natom_types;
  int found;

  /* Initialise variables */
  atom_type = 0;
  found = FALSE;
  natom_types = *ntyp;

  /* Loop over the stored atom names */
  for (itype = 0; itype < natom_types && found == FALSE; itype++)
    {
      /* If have a match, update the count */
      if (!strcmp(atom_name,stored_name[itype]))
        {
          /* Increment count */
          name_count[itype]++;

          /* Set the atom type */
          atom_type = itype;

          /* Set flag */
          found = TRUE;
        }
    }

  /* If don't already have this one, create a new slot for it */
  if (found == FALSE)
    {
      /* If not off end of array, store */
      if (natom_types < MAX_TYPES)
        {
          /* Store name name and increment count */
          strcpy(stored_name[natom_types],atom_name);
                     
          /* Increment count */
          name_count[natom_types] = 1;

          /* Set the atom type */
          atom_type = natom_types;

          /* Increment count of name types */
          natom_types++;
        }
    }

  /* Return current number of atom types */
  *ntyp = natom_types;

  /* Return the atom type */
  return(atom_type);
}
/***********************************************************************

find_atom  -  Find the given atom record

***********************************************************************/

struct atom *find_atom(struct atom *first_atom_ptr,int atom_number)
{
  struct atom *atom_ptr, *match_atom_ptr;

  /* Initialise variables */
  match_atom_ptr = NULL;

  /* Get pointer to the first atom */
  atom_ptr = first_atom_ptr;

  /* Loop over the atoms to write out */
  while (atom_ptr != NULL && match_atom_ptr == NULL)
    {
      /* Check if this is the required record */
      if (atom_ptr->atom_number == atom_number)
        {
          /* Save the atom pointer */
          match_atom_ptr = atom_ptr;
        }

      /* Get pointer to next atom record */
      atom_ptr = atom_ptr->next_atom_ptr;
    }
      
  /* Return match */
  return(match_atom_ptr);
}
/***********************************************************************

create_atom_record  -  Create a atom record

***********************************************************************/

struct atom *create_atom_record(struct atom **fst_atom_ptr,
                                struct atom **lst_atom_ptr,
                                int atom_number,char *atom_name,
                                double *coords,char *element,int type)
{
  int i;

  struct atom *atom_ptr;
  struct atom *first_atom_ptr, *last_atom_ptr;

  /* Initialise atom pointers */
  atom_ptr = NULL;
  first_atom_ptr = *fst_atom_ptr;
  last_atom_ptr = *lst_atom_ptr;

  /* Create atom entry */
  atom_ptr =
    (struct atom*)malloc(sizeof(struct atom));
  if (atom_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for");
      printf(" struct atom\n");
      printf("***        Program hetindex terminated with");
      printf(" error.\n");
      exit (1);
    }

  /* Store the known data */
  atom_ptr->atom_number = atom_number;
  strcpy(atom_ptr->atom_name,atom_name);
  strcpy(atom_ptr->element,element);
  atom_ptr->x = coords[X];
  atom_ptr->y = coords[Y];
  atom_ptr->z = coords[Z];
  atom_ptr->type = type;

  /* Initialise remaining fields */
  atom_ptr->hydrogen = FALSE;
  atom_ptr->checked = FALSE;
  atom_ptr->nequiv = 0;
  atom_ptr->residue_ptr = NULL;
  atom_ptr->dic_data_ptr = NULL;
  for (i = 0; i < MAXCONNECT; i++)
    atom_ptr->connect_atom_ptr[i] = NULL;

  /* Initialise pointer to next entry */
  atom_ptr->next_atom_ptr = NULL;

  /* If this is the first entry stored, then update
     pointer to head of linked list */
  if (first_atom_ptr == NULL)
    first_atom_ptr = atom_ptr;

  /* Otherwise, make previous entry point to the current
     one */
  else
    last_atom_ptr->next_atom_ptr = atom_ptr;
                      
  /* Save the current residue's pointer */
  last_atom_ptr = atom_ptr;

  /* Return current pointers */
  *fst_atom_ptr = first_atom_ptr;
  *lst_atom_ptr = last_atom_ptr;
  return(atom_ptr);
}
/***********************************************************************

create_bond_record  -  Create a new bond record

***********************************************************************/

struct bond *create_bond_record(struct bond **fst_bond_ptr,
                                struct bond **lst_bond_ptr,
                                int number,struct atom *atom1_ptr,
                                struct atom *atom2_ptr,int bond_type)
{
  struct bond *bond_ptr;
  struct bond *first_bond_ptr, *last_bond_ptr;

  /* Initialise bond pointers */
  bond_ptr = NULL;
  last_bond_ptr = *lst_bond_ptr;
  first_bond_ptr = *fst_bond_ptr;

  /* Create bond entry */
  bond_ptr = (struct bond*)malloc(sizeof(struct bond));
  if (bond_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for");
      printf(" struct bond\n");
      printf("***        Program catalact terminated with error.\n");
      exit (1);
    }

  /* Store the details for this entry */
  bond_ptr->number = number;
  bond_ptr->type = bond_type;
  bond_ptr->atom_ptr[0] = atom1_ptr;
  bond_ptr->atom_ptr[1] = atom2_ptr;

  /* Initialise pointer to next bond record */
  bond_ptr->next_bond_ptr = NULL;

  /* If this is the first entry stored, then update pointer to head
     of linked list */
  if (first_bond_ptr == NULL)
    first_bond_ptr = bond_ptr;

  /* Otherwise, make previous entry point to the current
     one */
  else
    last_bond_ptr->next_bond_ptr = bond_ptr;
                      
  /* Save the current bond's pointer */
  last_bond_ptr = bond_ptr;

  /* Return current pointers */
  *fst_bond_ptr = first_bond_ptr;
  *lst_bond_ptr = last_bond_ptr;
  return(bond_ptr);
}
/***********************************************************************

free_atom_memory  -  Free up memory used by residue atom

***********************************************************************/

void free_atom_memory(struct atom *first_atom_ptr)
{
  struct atom *atom_ptr, *next_atom_ptr;

  /* Get first atom */
  atom_ptr = first_atom_ptr;

  /* Loop through the atom records until done */
  while (atom_ptr != NULL)
    {
      /* Get pointer to the next atom record */
      next_atom_ptr = atom_ptr->next_atom_ptr;

      /* Free up memory taken by current atom */
      free(atom_ptr);

      /* Go to the next atom */
      atom_ptr = next_atom_ptr;
    }
}
/***********************************************************************

free_bond_memory  -  Free up memory used by residue bond

***********************************************************************/

void free_bond_memory(struct bond *first_bond_ptr)
{
  struct bond *bond_ptr, *next_bond_ptr;

  /* Get first bond */
  bond_ptr = first_bond_ptr;

  /* Loop through the bond records until done */
  while (bond_ptr != NULL)
    {
      /* Get pointer to the next bond record */
      next_bond_ptr = bond_ptr->next_bond_ptr;

      /* Free up memory taken by current bond */
      free(bond_ptr);

      /* Go to the next bond */
      bond_ptr = next_bond_ptr;
    }
}
/***********************************************************************

check_hetname  -  Check residue names to see if this is among them

***********************************************************************/

int check_hetname(struct residue *first_residue_ptr,char *id)
{
  char res_name[4];

  int len;
  int match;

  struct residue *residue_ptr;

  /* Initialise variables */
  match = FALSE;

  /* Right-justify the residue name */
  strcpy(res_name,id);
  len = strlen(res_name);
  if (len == 2)
    {
      strcpy(res_name," ");
      strcat(res_name,id);
    }
  else if (len == 1)
    {
      strcpy(res_name,"  ");
      strcat(res_name,id);
    }

  /* Get pointer to the first residue */
  residue_ptr = first_residue_ptr;

  /* Loop through all the residues to delete any without a ligand
     match */
  while (residue_ptr != NULL && match == FALSE)
    {
      /* If residue matches, then set flag */
      if (!strcmp(residue_ptr->res_name,res_name))
        match = TRUE;

      /* Get the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
    }

  /* Return whether match found */
  return(match);
}
/***********************************************************************

store_sdf_atoms  -  Store all the atom details for current dictionary
                    entry

***********************************************************************/

int store_sdf_atoms(struct atom *first_atom_ptr,
                    struct dictionary *dictionary_ptr,
                    struct dic_data **fst_dic_data_ptr,
                    struct dic_data **lst_dic_data_ptr)
{
  int natoms;

  struct atom *atom_ptr;
  struct dic_data *dic_data_ptr;
  struct dic_data *first_dic_data_ptr, *last_dic_data_ptr;

  /*  Initialise variables */
  natoms = 0;

  /* Initialise dic_data pointers */
  first_dic_data_ptr = *fst_dic_data_ptr;
  last_dic_data_ptr = *lst_dic_data_ptr;

  /* Initialise pointer to first atom */
  atom_ptr = first_atom_ptr;

  /* Loop through all the atoms to add to dictionary */
  while (atom_ptr != NULL)
    {
      /* Increment atom count */
      natoms++;

      /* Create new entry for this atom*/
      dic_data_ptr
        = create_dic_data_entry(&first_dic_data_ptr,
                                &last_dic_data_ptr,
                                dictionary_ptr->natoms);

      /* Store the atom name */
      strcpy(dic_data_ptr->atom_name,atom_ptr->atom_name);
      strcpy(dic_data_ptr->output_atom_name,atom_ptr->atom_name);

      /* Store the element */
      strcpy(dic_data_ptr->element,atom_ptr->element);

      /* Determine whether this is a hydrogen atom */
      if (!strncmp(atom_ptr->element," H",2))
        dic_data_ptr->hydrogen = TRUE;

      /* If this is the first dic_data record for this dictionary
         entry, then store pointer to it */
      if (dictionary_ptr->natoms == 0)
        dictionary_ptr->first_dic_data_ptr = dic_data_ptr;

      /* Increment atom-count for this dictionary entry */
      dictionary_ptr->natoms++;
      if (dic_data_ptr->hydrogen == FALSE)
        dictionary_ptr->non_hydrogen++;

      /* Store link from atom to dictionary data entry */
      atom_ptr->dic_data_ptr = dic_data_ptr;

      /* Get pointer to next atom in linked-list */
      atom_ptr = atom_ptr->next_atom_ptr;
    }               

  /* Return current pointers */
  *fst_dic_data_ptr = first_dic_data_ptr;
  *lst_dic_data_ptr = last_dic_data_ptr;

  /* Return number of atoms stored */
  return(natoms);
}
/***********************************************************************

store_sdf_bonds  -  Store all the bond details for current dictionary
                    entry

***********************************************************************/

int store_sdf_bonds(struct bond *first_bond_ptr,
                    struct dictionary *dictionary_ptr,
                    struct dic_data **fst_dic_data_ptr,
                    struct dic_data **lst_dic_data_ptr)
{
  int iatom, nbonds, nbound;
  int have_hydrogen;

  struct atom *atom_ptr, *other_atom_ptr;
  struct bond *bond_ptr;
  struct dic_data *dic_data_ptr;
  struct dic_data *first_dic_data_ptr, *last_dic_data_ptr;

  /*  Initialise variables */
  nbonds = 0;

  /* Initialise dic_data pointers */
  first_dic_data_ptr = *fst_dic_data_ptr;
  last_dic_data_ptr = *lst_dic_data_ptr;

  /* Initialise pointer to first bond */
  bond_ptr = first_bond_ptr;

  /* Loop through all the bonds to add to dictionary */
  while (bond_ptr != NULL)
    {
      /* Increment bond count */
      have_hydrogen = FALSE;
      nbonds++;

      /* Loop over the two atoms forming this bond */
      for (iatom = 0; iatom < 2; iatom++)
        {
          /* Get the two atoms */
          atom_ptr = bond_ptr->atom_ptr[iatom];
          other_atom_ptr = bond_ptr->atom_ptr[1 - iatom];

          /* Get this atom's dictionary data item */
          dic_data_ptr = atom_ptr->dic_data_ptr;

          /* If this is a hydrogen, set flag */
          if (dic_data_ptr->hydrogen == TRUE)
            have_hydrogen = TRUE;

          /* Store the connected atom */
          nbound = dic_data_ptr->nbound;
          if (nbound < MAXCONNECT)
            {
              /* Copy the connected atom across */
              strcpy(dic_data_ptr->boundto[nbound],
                     other_atom_ptr->atom_name);

              /* Store the bond order */
              dic_data_ptr->bond_order[nbound] = bond_ptr->type;

              /* Increment count of bonds */
              dic_data_ptr->nbound++;
            }

          /* Otherwise, show error message */
          else
            {
              printf("*** ERROR. Maximum number of atom connections ");
              printf("exceeded: %d\n",MAXCONNECT);
              printf("           Dictionary atom [%s] [%s]\n",
                     dictionary_ptr->res_name,dic_data_ptr->atom_name);
              exit(0);
            }
        }

      /* If neither atom is a hydrogen, then increment count of covalent
         bonds for this dictionary entry */
      if (have_hydrogen == FALSE)
        dictionary_ptr->nconnect++;

      /* Get pointer to next bond in linked-list */
      bond_ptr = bond_ptr->next_bond_ptr;
    }               

  /* Return current pointers */
  *fst_dic_data_ptr = first_dic_data_ptr;
  *lst_dic_data_ptr = last_dic_data_ptr;

  /* Return number of bonds stored */
  return(nbonds);
}
/***********************************************************************

read_sdf_dictionary  -  Read through HET Group Dictionary in SDF format

***********************************************************************/

int read_sdf_dictionary(char dictionary_name[],
                        struct residue *first_residue_ptr,
                        char *dic_res_name,
                        struct dictionary **fst_dictionary_ptr,
                        struct dictionary **lst_dictionary_ptr,
                        struct dic_data **fst_dic_data_ptr,
                        struct dic_data **lst_dic_data_ptr)
{
  char id[20], input_line[LINE_LEN + 1], new_id[20], number_string[20];
  char atom_name[10], element[3], stored_name[MAX_TYPES][5], tmp[3];

  char *token[MAX_TOKEN];

  int i, iatom, ibond, itoken, len, len1, lno;
  int natoms, natom_types, nbonds, ndic_entries, ndist, ntokens;
  int atom_type, bond_type, count, name_count[MAX_TYPES];
  int atom_number1, atom_number2;
  int have_key, match, wanted;

  double ave_dist, coords[3], dist_sqrd;

  struct atom *atom_ptr, *first_atom_ptr, *last_atom_ptr;
  struct atom *atom1_ptr, *atom2_ptr;
  struct bond *bond_ptr, *first_bond_ptr, *last_bond_ptr;
  struct dictionary *dictionary_ptr;
  struct dictionary *first_dictionary_ptr, *last_dictionary_ptr;
  struct dic_data *dic_data_ptr, *first_dic_data_ptr, *last_dic_data_ptr;

  /* Initialise variables */
  ave_dist = 0;
  have_key = FALSE;
  lno = -1;
  iatom = ibond = natoms = nbonds = 0;
  ndic_entries = 0;

  /* Initialise pointers */
  first_atom_ptr = last_atom_ptr = NULL;
  first_bond_ptr = last_bond_ptr = NULL;
  first_dictionary_ptr = NULL;
  last_dictionary_ptr = NULL;
  dic_data_ptr = NULL;
  first_dic_data_ptr = NULL;
  last_dic_data_ptr = NULL;
  
  /* Open the HET Group Dictionary */
  printf("  Opening SDF-format HET Group Dictionary ...\n");
  if ((fil_dic = fopen(dictionary_name,"r")) == NULL)
    {
      printf("\n*** Unable to open HET Group Dictionary [%s]\n",
             dictionary_name);
      return(ndic_entries);
    }
    
  /* Prepare to start reading in the data from the file */
  printf("\n");
  printf("    Scanning HET Group Dictionary ...\n");

  /* Loop while reading through the file */
  while (fgets(input_line,LINE_LEN,fil_dic) != NULL)
    {
      /* Truncate the string */
      len = string_truncate(input_line,LINE_LEN);

      /* Extract the tab-separated tokens */
      ntokens = get_tokens(input_line,token,MAX_TOKEN,' ');

      /* If looking for first line of this molecule, check if this
         could be it */
      if (lno == -1)
        {
          /* If enough tokens, then take to be first line */
          if (ntokens > 4)
            {
              /* Get the number of atoms and bonds */
              natoms = atoi(token[0]);
              nbonds = atoi(token[1]);

              /* Adjust for case where over 99 bonds and hence the
                 number of atoms and bonds run together */
              if (natoms > 999)
                {
                  /* Get number of atoms */
                  strncpy(number_string,input_line,3);
                  number_string[3] = '\0';
                  natoms = atoi(number_string);

                  /* Get number of bonds */
                  strncpy(number_string,input_line+3,3);
                  number_string[3] = '\0';
                  nbonds = atoi(number_string);
                }

              /* Re-initialise variables */
              ave_dist = 0;
              iatom = 0;
              ibond = 0;
              lno = 0;
              id[0] = '\0';
              have_key = FALSE;
              natom_types = 0;
              ndist = 0;
              for (i = 0; i < MAX_TYPES; i++)
                name_count[i] = 0;

              /* Re-initialise atom and bond pointers */
              first_atom_ptr = last_atom_ptr = NULL;
              first_bond_ptr = last_bond_ptr = NULL;
            }
        }

      /* If reading in atom line, store atom type and coords */
      else if (lno < natoms)
        {
          /* Get atom type and atom coords */
          strcpy(element,token[3]);
          to_upper(element,2);
          coords[0] = MOL_SCALE_FACTOR * atof(token[0]);
          coords[1] = MOL_SCALE_FACTOR * atof(token[1]);

          /* Identify this atom type and increment count */
          atom_type
            = get_attype(element,stored_name,name_count,&natom_types);
          count = name_count[atom_type];

          /* Form this atom's name */
          if (strlen(element) == 2)
            sprintf(atom_name,"%s%d  ",element,count);
          else
            sprintf(atom_name," %s%d  ",element,count);
          if (strlen(element) == 1)
            {
              strcpy(tmp,element);
              strcpy(element," ");
              strcat(element,tmp);
            }
          atom_name[4] = '\0';

          /* Increment atom counter */
          iatom++;

          /* Create atom record */
          atom_ptr
            = create_atom_record(&first_atom_ptr,&last_atom_ptr,
                                 iatom,atom_name,coords,element,
                                 atom_type);

          /* Increment line count */
          lno++;
        }

      /* If bond line, identify the two atoms involved */
      else if (lno <= natoms + nbonds)
        {
          /* Get the first atom number for this bond */
          atom_number1 = atoi(token[0]);

          /* Repeat for the second atom */
          atom_number2 = atoi(token[1]);

          /* Get the bond type */
          bond_type = atoi(token[2]);

          /* Check for case where the two atom numbers run together */
          if (atom_number1 > 999)
            {
              /* Get the length of this token */
              len = strlen(token[0]);

              /* Determine length of first atom number */
              len1 = len - 3;

              /* Extract the two atom numbers */
              if (len1 > 0)
                {
                  /* Extract the first number */
                  strcpy(number_string,token[0]);
                  number_string[len1] = '\0';
                  atom_number1 = atoi(number_string);

                  /* Extract the second number */
                  strcpy(number_string,token[0] + len1);
                  atom_number2 = atoi(number_string);
                }

              /* Get the bond type from the second token */
              bond_type = atoi(token[1]);
            }

          /* Get this atom record */
          atom1_ptr = find_atom(first_atom_ptr,atom_number1);

          /* Get this atom record */
          atom2_ptr = find_atom(first_atom_ptr,atom_number2);

          /* If both atoms found, then store the bond */
          if (atom1_ptr != NULL && atom2_ptr != NULL)
            {
              /* Increment count of bonds stored */
              ibond++;

              /* Create a bond record */
              bond_ptr
                = create_bond_record(&first_bond_ptr,&last_bond_ptr,
                                     ibond,atom1_ptr,atom2_ptr,
                                     bond_type);

              /* If atoms not hydrogens, calculate bond length */
              if (strcmp(atom1_ptr->element,"H") &&
                  strcmp(atom2_ptr->element,"H"))
                {
                  /* Calculate bond length and add to accumulated average */
                  dist_sqrd =
                    (atom1_ptr->x - atom2_ptr->x)
                    * (atom1_ptr->x - atom2_ptr->x)
                    + (atom1_ptr->y - atom2_ptr->y)
                    * (atom1_ptr->y - atom2_ptr->y)
                    + (atom1_ptr->z - atom2_ptr->z)
                    * (atom1_ptr->z - atom2_ptr->z);

                  /* Add to average */
                  ave_dist = ave_dist + sqrt(dist_sqrd);
                  ndist++;
                }
            }

          /* Increment line count */
          lno++;
        }

      /* If have hit the end of the molecule, process it */
      else if (!strcmp(input_line,"$$$$"))
        {
          /* If have found a unique id corresponding to a Het name, then
             save this molecule in the dictionary */
          if (id[0] != '\0')
            {
              /* Add this residue to the dictionary. linking to its entry
                 from all PDB residues for this residue type */
              dictionary_ptr
                = add_to_dictionary(id,dic_res_name,first_residue_ptr,
                                    &first_dictionary_ptr,
                                    &last_dictionary_ptr,&ndic_entries,
                                    &wanted);

              /* If wanted, then store the atom and bond details */
              if (wanted == TRUE)
                {
                  /* Store all the atoms */
                  natoms
                    = store_sdf_atoms(first_atom_ptr,dictionary_ptr,
                                      &first_dic_data_ptr,
                                      &last_dic_data_ptr);

                  /* Store all the bonds */
                  nbonds
                    = store_sdf_bonds(first_bond_ptr,dictionary_ptr,
                                      &first_dic_data_ptr,
                                      &last_dic_data_ptr);
                }
            }


          /* Reset line number */
          lno = -1;

          /* Free up memory taken by atom and bond records */
          if (first_atom_ptr != NULL)
            free_atom_memory(first_atom_ptr);
          if (first_bond_ptr != NULL)
            free_bond_memory(first_bond_ptr);
        }

      /* Otherwise, look for id corresponding to Het id */
      else
        {
          /* Initialise identifier */
          new_id[0] = '\0';
          wanted = TRUE;

          /* Loop over the tokens in the current line */
          for (itoken = 0; itoken < ntokens && wanted == TRUE; itoken++)
            {
              /* Get the length of the token */
              len1 = strlen(token[itoken]);

              /* If previous token was a key token, see if current
                 token might be a residue name */
              if (have_key == TRUE)
                {
                  /* If token has the right length, then store */
                  if (len1 < 4)
                    {
                      /* Store the token */
                      strcpy(new_id,token[0]);

                      /* See if it matches any residue names */
                      match = check_hetname(first_residue_ptr,new_id);

                      /* If residue name matches, then save the id */
                      if (match == TRUE)
                        {
                          /* If already have an id, then cannot uniquely
                             identify this molecule */
                          if (id[0] != '\0')
                            {
                              /* Show warning message */
                              printf("*** Warning. Molecule id in SDF "
                                     "file is confusing: [%s or %s]. "
                                     "Ignoring\n",id,new_id);

                              /* Unset the id and wanted flag */
                              id[0] = '\0';
                              wanted = FALSE;
                            }

                          /* Otherwise, save the id */
                          strcpy(id,new_id);
                          wanted = TRUE;
                        }
                    }

                  /* Unset flag */
                  have_key = FALSE;
                }

              /* Otherwise, see if this is a token */
              else if (token[itoken][0] == '<' &&
                       token[itoken][len1 - 1] == '>')
                have_key = TRUE;
            }
        }
    }

  /* Print counts of data read in */
  printf("\n");
  printf("      Number of dictionary entries read in  = %7d\n",
         ndic_entries);

  /* Close the dictionary file */
  fclose(fil_dic);

  /* Return current pointers */
  *fst_dictionary_ptr = first_dictionary_ptr;
  *lst_dictionary_ptr = last_dictionary_ptr;
  *fst_dic_data_ptr = first_dic_data_ptr;
  *lst_dic_data_ptr = last_dic_data_ptr;

  /* Return number of dictionary entries stored */
  return(ndic_entries);
}
/* <--v.5.1 */
/***********************************************************************

generate_duplicate_entries  -  Generate dictionary duplicate entries
                               - one for each residue to be compared
                               against

***********************************************************************/

void generate_duplicate_entries(struct dictionary *first_dictionary_ptr)
{
  int i, iatom, natoms;

  struct residue *residue_ptr;
  struct dictionary *dictionary_ptr;
  struct dic_data *dic_data_ptr, *new_dic_data_ptr, *last_dic_data_ptr;

  /* Initialise pointer */
/* v.5.3.2--> */
  residue_ptr = NULL;
/* <--v.5.3.2 */

  /* Get first dictionary entry */
  dictionary_ptr = first_dictionary_ptr;

  /* Update residue pointer */
/* v.5.3.2--> */
  if (dictionary_ptr != NULL)
    {
/* <--v.5.3.2 */
      residue_ptr = dictionary_ptr->residue_ptr;
      if (residue_ptr != NULL)
        residue_ptr->dictionary_ptr = dictionary_ptr;
/* v.5.3.2--> */
    }
/* <--v.5.3.2 */

  /* Start at second dictionary entry */
  if (dictionary_ptr != NULL)
    dictionary_ptr = dictionary_ptr->next_dictionary_ptr;

  /* Loop through all subsequent dictionary entries to duplicate the
     first entry's atoms and bonds */
  while (dictionary_ptr != NULL)
    {
      /* Update residue pointer */
      residue_ptr = dictionary_ptr->residue_ptr;
      if (residue_ptr != NULL)
        residue_ptr->dictionary_ptr = dictionary_ptr;

      /* Transfer relevant atom and connectivity info */
      dictionary_ptr->natoms = first_dictionary_ptr->natoms;
      dictionary_ptr->non_hydrogen = first_dictionary_ptr->non_hydrogen;
      dictionary_ptr->nconnect = first_dictionary_ptr->nconnect;

      /* Transfer the molecule name */
      store_string(&(dictionary_ptr->mol_name),
                   first_dictionary_ptr->mol_name);

      /* Initialise pointers */
      last_dic_data_ptr = NULL;

      /* Get the first dictionary entry's first atom */
      dic_data_ptr = first_dictionary_ptr->first_dic_data_ptr;
      natoms = first_dictionary_ptr->natoms;

      /* Loop through each of the atoms in turn and process */
      for (iatom = 0; iatom < natoms && dic_data_ptr != NULL; iatom++)
        {
          /* Create an entry for this atom's connectivities */
          new_dic_data_ptr
            = create_dic_data_entry(&(dictionary_ptr->first_dic_data_ptr),
                                    &last_dic_data_ptr,iatom);

          /* Copy the relevant details across */
          strcpy(new_dic_data_ptr->atom_name,dic_data_ptr->atom_name);
          strcpy(new_dic_data_ptr->output_atom_name,dic_data_ptr->atom_name);
          strcpy(new_dic_data_ptr->element,dic_data_ptr->element);
          new_dic_data_ptr->hydrogen = dic_data_ptr->hydrogen;
          new_dic_data_ptr->nbound = dic_data_ptr->nbound;

          /* Transfer the connectivity information */
          for (i = 0; i < dic_data_ptr->nbound; i++)
            {
              strcpy(new_dic_data_ptr->boundto[i],dic_data_ptr->boundto[i]);
              new_dic_data_ptr->bond_order[i] = dic_data_ptr->bond_order[i];
            }

          /* Get pointer to the next HET group atom in the list */
          dic_data_ptr = dic_data_ptr->next_dic_data_ptr;
        }

      /* Get pointer to the next dictionary entry */
      dictionary_ptr = dictionary_ptr->next_dictionary_ptr;
    }
}
/***********************************************************************

dic_connectivities  -  Scan through the connectivity info for each stored
                       CONEC record from the HET Groups Dictionary and
                       match up the other atoms that each atom is bonded to

***********************************************************************/

void dic_connectivities(struct dictionary *first_dictionary_ptr)
{
  int found, iatom, iconnect, jatom, natoms, nbonds;

  struct dictionary *dictionary_ptr;
  struct dic_data *dic_data_ptr, *other_dic_data_ptr;

  /* Initialise HET group pointer to start of linked list */
  dictionary_ptr = first_dictionary_ptr;

  /* Loop through all the HET group residues read in from the dictionary */
  while (dictionary_ptr != NULL)
    {
      /* Initialise count of bonds in dictionary entry */
      nbonds = 0;

      /* Get the number of atoms making up this HET group and the pointer
         to the first one */
      dic_data_ptr = dictionary_ptr->first_dic_data_ptr;
      natoms = dictionary_ptr->natoms;

      /* Loop through each of the atoms in turn and process */
      for (iatom = 0; iatom < natoms && dic_data_ptr != NULL; iatom++)
        {
          /* Loop through the list of atoms bonded to this one */
          for (iconnect = 0; iconnect < dic_data_ptr->nbound; iconnect++)
            {
              /* Initialise search list to first of the atoms in this
                 HET group */
              other_dic_data_ptr = dictionary_ptr->first_dic_data_ptr;
              found = FALSE;

              /* Loop through all the atoms in the HET group until find
                 the atom we're after */
              for (jatom = 0; jatom < natoms && found == FALSE; jatom++)
                {
                  /* Check whether this is the atom we want */
                  if (!strncmp(dic_data_ptr->boundto[iconnect],
                               other_dic_data_ptr->atom_name,4))
                    {
                      /* Have a match, so store the pointer matched atom */
                      dic_data_ptr->connect_atom_ptr[iconnect]
                        = other_dic_data_ptr;
                      dic_data_ptr->boundto_num[iconnect]
                        = other_dic_data_ptr->number;
                      found = TRUE;

                          /* Increment count of bonds */
                          nbonds++;
                    }

                  /* Get pointer to the next atom in the list */
                  other_dic_data_ptr = other_dic_data_ptr->next_dic_data_ptr;
                }

              /* If the atom not found, then show warning */
              if (found == FALSE)
                {
                  printf("*** Warning. Error in HET Group Dictionary");
                  printf(" for residue-type [%s]\n",
                         dictionary_ptr->res_name);
                  printf("***          Connect atom %d [%s] undefined",
                         iconnect + 1,dic_data_ptr->boundto[iconnect]);
                  printf(" for atom type [%s]\n",dic_data_ptr->atom_name);
                }
            }

          /* Get pointer to the next HET group atom in the list */
          dic_data_ptr = dic_data_ptr->next_dic_data_ptr;
        }

    /* Halve the number of bonds */
     /* nbonds = nbonds / 2; */

/* v.5.0--> */
      /* Store the number of bonds */
      dictionary_ptr->nbonds = nbonds / 2;
/* <--v.5.0 */

      /* Get pointer to the next HET group entry in the list */
      dictionary_ptr = dictionary_ptr->next_dictionary_ptr;
    }
}
/***********************************************************************

open_sumfile  -  Open the hbadd.sum file

***********************************************************************/

void open_sumfile(char *sum_name,char *pdb_name)
{
  /* Open the output hbadd.sum file */
  if ((fil_sum = fopen(sum_name,"w")) == NULL)
    {
      printf("\n*** Unable to open summary file [%s] for output\n",sum_name);
      exit(1);
    }

  /* Write PDB filename to the first line of the summary file */
  fprintf(fil_sum,"%s\n",pdb_name);

  /* Print column headings */
  fprintf(fil_sum,"Residue    Atoms Dic Match Exact Mapped Missed");
  fprintf(fil_sum,"  Missing     Bonds Dic Match\n");
}    
/***********************************************************************

open_mapfile  -  Open the hbadd.map file

***********************************************************************/

void open_mapfile(char map_name[],char pdb_name[FILENAME_LEN])
{
  /* Open the output hbadd.map file */
  if ((fil_map = fopen(map_name,"w")) == NULL)
    {
      printf("\n*** Unable to open file [%s] for output\n",map_name);
      exit(1);
    }

  /* Write PDB filename to the first line of the mapmary file */
  fprintf(fil_map,"%s\n",pdb_name);

  /* Print column headings */
  fprintf(fil_map,"PDB file atom   ->  Dictionary atom\n");
}    
/***********************************************************************

print_missing_hets  -  Print this HET group as missing from the HET
                       Group Dictionary, unless a warning has already
                       been printed for this residue type

***********************************************************************/

void print_missing_hets(struct residue *residue_ptr)
{
  int nmatch;

  struct residue *residue_comp_ptr;


  /* Flag this as the first of its kind */
  nmatch = 1;
  residue_ptr->nmatch = nmatch;

  /* Initialise search of residues beyond the current one */
  residue_comp_ptr = residue_ptr->next_residue_ptr;

  /* Loop through all the following HET group residues to find any
     duplicates of the current one */
  while (residue_comp_ptr != NULL)
    {
      /* If residue name the same, then flag as an additional copy */
      if (!strncmp(residue_comp_ptr->res_name,residue_ptr->res_name,3))
        {
          /* Increment matches count and store */
          nmatch++;
          residue_comp_ptr->nmatch = nmatch;
        }

      /* Get pointer to the next residue in the list */
      residue_comp_ptr = residue_comp_ptr->next_residue_ptr;
    }
}
/***********************************************************************

check_stored_conecs  -  Loop through all the stored connections from the
                        CONECT records in the PDB file and see if there
                        is one for the two atoms in question

***********************************************************************/

int check_stored_conecs(struct connect *first_connect_ptr,
                        int atom_number1,int atom_number2)
{
  int connected;

  struct connect *connect_ptr;

  /* Initialise variables */
  connected = FALSE;

  /* Initialise connect-pointer to start of linked list */
  connect_ptr = first_connect_ptr;

  /* Loop through all the HET group connects */
  while (connect_ptr != NULL && connected == FALSE)
    {
      /* Check whether this bond corresponds to the two atoms in question */
      if ((atom_number1 == connect_ptr->atom_number1 &&
           atom_number2 == connect_ptr->atom_number2) ||
          (atom_number1 == connect_ptr->atom_number2 &&
           atom_number2 == connect_ptr->atom_number1))
        connected = TRUE;

      /* Get pointer to the next connect in the list */
      connect_ptr = connect_ptr->next_connect_ptr;
    }

  /* Return the answer */
  return(connected);
}
/***********************************************************************

generate_connectivities  -  Generate all the atom-atom connectivities
                            for the fake dictionary entry using the
                            atom-atom distances and CONECT records in
                            the PDB 

***********************************************************************/

void generate_connectivities(struct residue *residue_ptr,
                             struct dictionary *dictionary_ptr,
                             struct connect *first_connect_ptr)
{
  int iatom, iconnect, jatom, natoms;
  int connected, hydrogen1, hydrogen2;

  float connect_dist2, distance2, x1, x2, y1, y2, z1, z2;

  struct atom *atom_ptr, *other_atom_ptr;
  struct dic_data *dic_data_ptr, *other_dic_data_ptr;

  /* Initialise pointer to first atom in the dictionary list of
     atoms making up this HET group */
  dic_data_ptr = dictionary_ptr->first_dic_data_ptr;
  natoms = dictionary_ptr->natoms;

  /* Loop through each of the atoms in turn and check each of its bonds */
  for (iatom = 0; iatom < natoms && dic_data_ptr != NULL; iatom++)
    {
      /* Get the pointer to the atom coordinates */
      atom_ptr = dic_data_ptr->atom_ptr;
      if (atom_ptr != NULL)
        {
          /* Get the atom's coordinates */
          x1 = atom_ptr->x;
          y1 = atom_ptr->y;
          z1 = atom_ptr->z;

          /* Check if this atom is a hydrogen */
          hydrogen1 = FALSE;
          if (!strcmp(atom_ptr->element," H"))
            hydrogen1 = TRUE;

          /* Loop through all this residue's atoms to find all atoms
             bonded to this one */

          /* Get pointer to first of residue's atoms */
          other_atom_ptr = residue_ptr->first_atom_ptr;

          /* Loop through the residue's atoms */
          for (jatom = 0; jatom < natoms && other_atom_ptr != NULL;
               jatom++)
            {
              /* Process if not the same atom */
              if (other_atom_ptr != atom_ptr)
                {
                  /* Check whether the two atoms are covalently
                     bonded */

                  /* Get the second atom's coordinates */
                  x2 = other_atom_ptr->x;
                  y2 = other_atom_ptr->y;
                  z2 = other_atom_ptr->z;
                  
                  /* Check if this atom is a hydrogen */
                  hydrogen2 = FALSE;
                  if (!strcmp(other_atom_ptr->element," H"))
                    hydrogen2 = TRUE;

                  /* Calculate the distance between them */
                  distance2 = (x1 - x2) * (x1 - x2)
                    + (y1 - y2) * (y1 - y2)
                    + (z1 - z2) * (z1 - z2);

                  /* Get the other atom's dictionary entry */
                  other_dic_data_ptr = other_atom_ptr->dic_data_ptr;

                  /* If either atom is a hydrogen, then reduce
                     connect distance */
                  if (hydrogen1 == TRUE || hydrogen2 == TRUE)
                    connect_dist2 = (float) HATOM_DIST2;
                  else
                    connect_dist2 = COVAL_DIST2;

                  /* If bond distance is not within allowed tolerances,
                     then check the connectivities as read in from
                     the CONECT records in the PDB file */
                  connected = FALSE;
                  if (distance2 > connect_dist2)
                    {
                      /* Check for a CONECT record between these two
                         atoms */
                      connected
                        = check_stored_conecs(first_connect_ptr,
                                              atom_ptr->atom_number,
                                              other_atom_ptr->atom_number);
                    }
                  else
                    connected = TRUE;

                  /* If atoms are connected, then add this bond to the
                     dictionary entry */
                  if (connected == TRUE)
                    {
                      /* Get next available slot for connectivity data */
                      iconnect = dic_data_ptr->nbound;

                      /* Check that maximum number of possible connections
                         not exceeded */
                      if (iconnect > MAXCONNECT - 1)
                        {
                          printf("*** Warning. Maximum number of bonds");
                          printf(" exceeded for atom [%s %s %s %c] %d",
                                 atom_ptr->atom_name,
                                 atom_ptr->residue_ptr->res_name,
                                 atom_ptr->residue_ptr->res_num,
                                 atom_ptr->residue_ptr->chain,
                                 iconnect + 1);
                        }

                      /* Otherwise, store this bond */
                      else
                        {
                          /* Store the details of the other atom */
                          dic_data_ptr->connect_atom_ptr[iconnect]
                            = other_dic_data_ptr;
                          dic_data_ptr->boundto_num[iconnect]
                            = other_atom_ptr->atom_number;
                          strncpy(dic_data_ptr->boundto[iconnect],
                                  other_atom_ptr->atom_name,4);
                          dic_data_ptr->boundto[iconnect][4] = '\0';
                      
                          /* Increment count of bonded connections */
                          dic_data_ptr->nbound++;
                        }
                    }
                }

              /* Get pointer to next atom */
              other_atom_ptr = other_atom_ptr->next_atom_ptr;
            }
        }
      
      /* Get pointer to the next HET group atom in the list */
      dic_data_ptr = dic_data_ptr->next_dic_data_ptr;
    }
}
/***********************************************************************

generate_fake_entry  -  Check that the atoms in the HET group from the
                        PDB file match those in the dictionary definition

***********************************************************************/

void generate_fake_entry(struct residue *residue_ptr,
                         struct residue *first_residue_ptr,
                         struct dictionary **fst_dictionary_ptr,
                         struct dictionary **lst_dictionary_ptr,
                         int *ndic_entries,
                         struct dic_data **lst_dic_data_ptr,
                         struct connect *first_connect_ptr)
{
  int iatom, natoms, ndic;
  int already_got_it;

  struct atom *atom_ptr;
  struct dictionary *dictionary_ptr, *save_dictionary_ptr;
  struct dictionary *first_dictionary_ptr, *last_dictionary_ptr;
  struct dic_data *dic_data_ptr, *last_dic_data_ptr;
  struct residue *check_residue_ptr;

  /* Initialise variables */
  ndic = *ndic_entries;

  /* Initialise pointers */
  first_dictionary_ptr = *fst_dictionary_ptr;
  last_dictionary_ptr = *lst_dictionary_ptr;
  last_dic_data_ptr = *lst_dic_data_ptr;

  /* Show message */
  printf("    Generating fake HBPLUS entry from PDB file for residue ");
  printf("type [%s]\n",residue_ptr->res_name);

  /* Check that don't already have a dictionary entry for this residue
     type */
  already_got_it = FALSE;
  check_residue_ptr = first_residue_ptr;
  save_dictionary_ptr = NULL;

  /* Loop through all the residues */
  while (check_residue_ptr != NULL && already_got_it == FALSE)
    {
      /* If residue-type matches and have a dictionary pointer, then
         must already have an entry in the dictionary */
      if (!strncmp(check_residue_ptr->res_name,residue_ptr->res_name,3) &&
          check_residue_ptr->dictionary_ptr != NULL &&
          check_residue_ptr->nmatch == 1)
        {
          /* Store the pointer */
          save_dictionary_ptr = check_residue_ptr->dictionary_ptr;
          already_got_it = TRUE;
        }

      /* Get pointer to the next residue */
      check_residue_ptr = check_residue_ptr->next_residue_ptr;
    }

  /* If already have this residue type, then don't need to create another
     entry */
  if (already_got_it == TRUE)
    {
      /* Create the appropriate link */
      residue_ptr->dictionary_ptr = save_dictionary_ptr;
      return;
    }

  /* Create a new (fake) dictionary entry corresponding to the unmatched
     residue */
  dictionary_ptr
    = create_dictionary_entry(&first_dictionary_ptr,
                              &last_dictionary_ptr,residue_ptr->res_name,
                              &ndic);

 /* Store the details for this entry */
  dictionary_ptr->natoms = residue_ptr->natoms;
  dictionary_ptr->non_hydrogen = residue_ptr->non_hydrogen;
  dictionary_ptr->fake_entry = TRUE;
  dictionary_ptr->residue_ptr = residue_ptr;

  /* Add link from residue to dictionary entry */
  residue_ptr->dictionary_ptr = dictionary_ptr;

  /* Initialise search to first PDB atom in the current residue */
  atom_ptr = residue_ptr->first_atom_ptr;
  natoms = residue_ptr->natoms;

  /* Loop through all the atoms, creating a dictionary entry for 
     each one */
  for (iatom = 0; iatom < natoms && atom_ptr != NULL; iatom++)
    {
      /* Create entry in the atom-data linked list for this atom */
      dic_data_ptr
        = create_dic_data_entry(&dictionary_ptr->first_dic_data_ptr,
                                &last_dic_data_ptr,natoms);

      /* Store the atom name */
      strcpy(dic_data_ptr->atom_name,atom_ptr->atom_name);
      strcpy(dic_data_ptr->output_atom_name,atom_ptr->atom_name);

      /* Store the element */
      strcpy(dic_data_ptr->element,atom_ptr->element);

      /* Determine whether this is a hydrogen atom */
      if (!strcmp(atom_ptr->element," H"))
        {
          dic_data_ptr->hydrogen = TRUE;
          dictionary_ptr->have_hydrogens = TRUE;
        }

      /* Store pointer to equivalenced atom */
      dic_data_ptr->number = atom_ptr->atom_number;
      dic_data_ptr->atom_ptr = atom_ptr;

      /* Store pointer the other way */
      atom_ptr->dic_data_ptr = dic_data_ptr;

      /* Save pointer to the current entry */
      last_dic_data_ptr = dic_data_ptr;

      /* Get pointer to the next coordinate record in the list */
      atom_ptr = atom_ptr->next_atom_ptr;
    }

  /* Calculate all the connectivities using distances and CONECT data */
  generate_connectivities(residue_ptr,dictionary_ptr,first_connect_ptr);

  /* Mark the residue as having a fake dictionary entry */
  residue_ptr->fake = TRUE;

  /* Return pointers */
  *lst_dic_data_ptr = last_dic_data_ptr;
  *fst_dictionary_ptr = first_dictionary_ptr;
  *lst_dictionary_ptr = last_dictionary_ptr;

  /* Return number of dictionary entries */
  *ndic_entries = ndic;
}
/***********************************************************************

calc_covalent_connectivities  -  Calculate the covalent connectivities
                                 of all the atoms in the given residue

***********************************************************************/

int calc_covalent_connectivities(struct residue *residue_ptr,
                                 struct connect *first_connect_ptr)
{
  int connected, iatom, iconnect, jatom, natoms, nconnect;

  float distance2, x1, x2, y1, y2, z1, z2;

  struct atom *atom_ptr, *other_atom_ptr;

  /* Initialise variables */
  nconnect = 0;

  /* Initialise pointer to first atom in the residue */
  atom_ptr = residue_ptr->first_atom_ptr;
  natoms = residue_ptr->natoms;

  /* Loop through each of the atoms in turn and check each of its bonds */
  for (iatom = 0; iatom < natoms && atom_ptr != NULL; iatom++)
    {
      /* Get this atom's coordinates */
      x1 = atom_ptr->x;
      y1 = atom_ptr->y;
      z1 = atom_ptr->z;

      /* Initialise all its connectivity pointers */
      for (iconnect = 0; iconnect < MAXCONNECT; iconnect++)
        atom_ptr->connect_atom_ptr[iconnect] = NULL;

      /* Re-initialise connections counter */
      iconnect = 0;

      /* Loop over all the other atoms */
      other_atom_ptr = residue_ptr->first_atom_ptr;

      /* Loop through each of the atoms in turn and check each of
         its bonds */
      for (jatom = 0; jatom < natoms && other_atom_ptr != NULL; jatom++)
        {
          /* Proceed if not the same as the first atom */
          if (other_atom_ptr != atom_ptr)
            {
              /* Get this atom's coordinates */
              x2 = other_atom_ptr->x;
              y2 = other_atom_ptr->y;
              z2 = other_atom_ptr->z;

              /* Check whether the two atoms are covalently bonded in
                 the structure */

              /* Calculate the distance between them */
              distance2 = (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2)
                + (z1 - z2) * (z1 - z2);

              /* If bond distance is not within allowed tolerances, then
                 check the connectivities as read in from the CONECT
                 records in the PDB file */
              connected = FALSE;
              if (distance2 > COVAL_DIST2)
                {
                  /* Check for a CONECT record between these two
                     atoms */
                  connected
                    = check_stored_conecs(first_connect_ptr,
                                          atom_ptr->atom_number,
                                          other_atom_ptr->atom_number);
                }
              else
                connected = TRUE;

              /* If atoms are bonded, add to the first atom's list of
                 connections */
              if (connected == TRUE)
                {
                  if (iconnect > MAXCONNECT - 1)
                    {
                      printf("*** Warning. Maximum number of bonds");
                      printf(" exceeded for atom [%s %s %s %c] %d\n",
                             atom_ptr->atom_name,
                             atom_ptr->residue_ptr->res_name,
                             atom_ptr->residue_ptr->res_num,
                             atom_ptr->residue_ptr->chain,
                             iconnect + 1);
                    }

                  /* Otherwise, store this bond */
                  else
                    {
                      /* Store pointer to bonded atom */
                      atom_ptr->connect_atom_ptr[iconnect] = other_atom_ptr;
                      iconnect++;

                      /* If neither atom is a hydrogen, then increment
                         count of covalent bonds */
                      if (strcmp(atom_ptr->element," H") &&
                          strcmp(other_atom_ptr->element," H") &&
                          atom_ptr->atom_number < other_atom_ptr->atom_number)
                          nconnect++;
                    }
                }
            }

          /* Get pointer to the next atom in the list */
          other_atom_ptr = other_atom_ptr->next_atom_ptr;
        }
      
      /* Get pointer to the next atom in the list */
      atom_ptr = atom_ptr->next_atom_ptr;
    }

  /* Return number of covalent bonds found */
  return(nconnect);
}
/***********************************************************************

name_match_score  -  CCalculate a "match" score based on how similar
                     the PDB and dictionary atom names are in the
                     current equivalence 

***********************************************************************/

int name_match_score(char atom_name[5],char dic_atom_name[5])
{
  char ch, name1[5], name2[5];

  int score;
  int ipos, jpos;
  int found;

  /* Initialise variables */
  score = 0;

  /* Transfer names to temporary areas */
  strcpy(name1,atom_name);
  strcpy(name2,dic_atom_name);

  /* Check for an exact match */
  if (!strcmp(name1,name2))
    score = 8;

  /* Otherwise calculate similarity score from partial matches */
  else
    {
      /* Loop over characters in first name */
      for (ipos = 0; ipos < 4; ipos++)
        {
          /* Get the current character */
          ch = name1[ipos];

          /* Initialise flag */
          found = FALSE;

          /* Check whether the match occurs at the equivalent character
             position */
          if (ch == name2[ipos])
            {
              /* Score 2 */
              score = score + 2;

              /* Mark as found */
              found = TRUE;

              /* Blank out the character in the second name */
              name2[ipos] = '\0';
            }

          /* Loop over the characters in the second name */
          for (jpos = 0; jpos < 4 && found == FALSE; jpos++)
            {
              /* Check for match */
              if (ch == name2[jpos])
                {
                  /* Score 1 */
                  score = score + 1;

                  /* Mark as found */
                  found = TRUE;

                  /* Blank out the character in the second name */
                  name2[jpos] = '\0';
                }
            }
        }
    }

  /* Return the match score */
  return(score);
}
/***********************************************************************

create_equiv_atoms  -  Create equiv_atoms of the graph, each corresponding to a
                       pairing between similar atoms of the residue and from
                       the Het Group Dictionary entry

***********************************************************************/

void create_equiv_atoms(struct residue *residue_ptr,
                        struct equiv_atom **first_equiv_atom_ptr,
/* v.5.0--> */
/*                        int *nequiv_atoms,int *mx_eq) */
                        int *nequiv_atoms,int *mx_eq,int match_names)
/* <--v.5.0 */
{
  char atom_name[5], dic_atom_name[5];

  int iatom, jatom, max_equiv, ndic_atoms, natoms;

  struct atom *atom_ptr;
  struct dictionary *dictionary_ptr;
  struct dic_data *dic_data_ptr;
  struct equiv_atom *last_equiv_atom_ptr, *equiv_atom_ptr;

  /* Initialise variables */
  last_equiv_atom_ptr = NULL;
  *first_equiv_atom_ptr = NULL;
  *nequiv_atoms = 0;
  max_equiv = 0;

  /* Get dictionary entry corresponding to this residue */
  dictionary_ptr = residue_ptr->dictionary_ptr;

  /* Initialise pointer to first atom in the residue */
  atom_ptr = residue_ptr->first_atom_ptr;
  natoms = residue_ptr->natoms;

  /* Loop through each of the atoms in turn */
  for (iatom = 0; iatom < natoms && atom_ptr != NULL; iatom++)
    {
      /* Get the atom name */
      strncpy(atom_name,atom_ptr->atom_name,4);
      atom_name[4] = '\0';

      /* Skip atom if it is a hydrogen */
      if (atom_ptr->hydrogen == FALSE)
        {
          /* Get the first of the atoms corresponding to the dictionary
             entry */
          dic_data_ptr = dictionary_ptr->first_dic_data_ptr;
          ndic_atoms = dictionary_ptr->natoms;

          /* Loop through each of the atoms in the dictionary entry */
          for (jatom = 0; jatom < ndic_atoms && dic_data_ptr != NULL;
               jatom++)
            {
              /* Get this atom's name */
              strncpy(dic_atom_name,dic_data_ptr->atom_name,4);
              dic_atom_name[4] = '\0';

              /* If the two atoms are of the same type, create a
                 graph-equiv_atom to correspond to their potentially
                 being equivalent */
              if (!strcmp(atom_ptr->element,dic_data_ptr->element))
                {
                  /* Create a equiv_atom equivalencing these two atoms */
                  equiv_atom_ptr
                    = (struct equiv_atom*)malloc(sizeof(struct equiv_atom));
                  if (equiv_atom_ptr == NULL)
                    {
                      printf("*** ERROR. Unable to allocate memory");
                      printf(" for struct equiv_atom\n");
                      printf("***        Program hbadd terminated with");
                      printf(" error.\n");
                      exit (1);
                    }

                  /* If this is the first equiv_atom, then update pointer
                     to head of linked list */
                  if (*first_equiv_atom_ptr == NULL)
                    *first_equiv_atom_ptr = equiv_atom_ptr;

                  /* Otherwise, make previous equiv_atom point to the current
                     one */
                  else
                    last_equiv_atom_ptr->next_equiv_atom_ptr = equiv_atom_ptr;

                  /* Save the current record's pointer */
                  last_equiv_atom_ptr = equiv_atom_ptr;

                  /* Store the two potentially equivalent atoms */
                  equiv_atom_ptr->atom_ptr = atom_ptr;
                  equiv_atom_ptr->dic_data_ptr = dic_data_ptr;

                  /* Fill in the remaining data items */
                  equiv_atom_ptr->equiv_atom_number = *nequiv_atoms;
                  equiv_atom_ptr->checked = FALSE;
                  equiv_atom_ptr->in_best_clique = FALSE;
                  equiv_atom_ptr->tree_level = 0;
                  equiv_atom_ptr->links_combination = 0;
                  equiv_atom_ptr->nlinks = 0;
                  equiv_atom_ptr->first_equiv_bond_ptr = NULL;
                  equiv_atom_ptr->next_equiv_atom_ptr = NULL;

                  /* Calculate a "match" score based on how similar
                     the PDB and dictionary atom names are in this
                     equivalence */
                  equiv_atom_ptr->match_score
                    = name_match_score(atom_name,dic_atom_name);
/* v.5.0--> */
                  if (match_names == FALSE)
                    equiv_atom_ptr->match_score = 0;
/* <--v.5.0 */

                  /* Increment count of number of times this atom has
                     been equivalenced */
                  atom_ptr->nequiv++;

                  /* If this is the largest number of equivalences so
                     far, then store */
                  if (atom_ptr->nequiv > max_equiv)
                    max_equiv = atom_ptr->nequiv;

                  /* Increment count of graph equiv_atoms created */
                  (*nequiv_atoms)++;
                }

              /* Get pointer to the next HET group atom in the list */
              dic_data_ptr = dic_data_ptr->next_dic_data_ptr;
            }
        }
      
      /* Get pointer to the next atom in the list */
      atom_ptr = atom_ptr->next_atom_ptr;
    }

  /* Return maximum number of atom equivalences */
  *mx_eq = max_equiv;
}
/***********************************************************************

create_equiv_bond  -  Create a new equivalent bond record

***********************************************************************/

struct equiv_bond
*create_equiv_bond(struct equiv_atom *equiv_atom_ptr,
                   struct equiv_atom *other_equiv_atom_ptr,
                   struct equiv_bond **fst_equiv_bond_ptr,
                   struct equiv_bond **lst_equiv_bond_ptr,
                   int *nequiv_bonds,int max_equiv)
{
  int other_score, score;

  struct atom *atom_ptr;
  struct equiv_bond *equiv_bond_ptr;
  struct equiv_bond *first_equiv_bond_ptr, *last_equiv_bond_ptr;
  struct equiv_bond *other_equiv_bond_ptr;
  struct equiv_bond *before_equiv_bond_ptr, *after_equiv_bond_ptr;

  /* Initialise pointers */
  first_equiv_bond_ptr = *fst_equiv_bond_ptr;
  last_equiv_bond_ptr = *lst_equiv_bond_ptr;

  /* Allocate memory for structure to hold equiv_bond info */
  equiv_bond_ptr = (struct equiv_bond *) malloc(sizeof(struct equiv_bond));
  if (equiv_bond_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for struct equiv_bond\n");
      printf("***        Program hbadd terminated with error.\n");
      exit (1);
    }

  /* If this is the first equivalent bond then update pointer to head
     of linked list */
  if (first_equiv_bond_ptr == NULL)
    first_equiv_bond_ptr = equiv_bond_ptr;

  /* Otherwise, get the previous equivalent bond and get it to point to
     the current equivalent bond */
  else
    {
      /* Update its pointer from previous equivalent bond */
      last_equiv_bond_ptr->next_equiv_bond_ptr = equiv_bond_ptr;
    }

  /* Calculate a "match" score for this equivalence */
  score = equiv_atom_ptr->match_score + other_equiv_atom_ptr->match_score;

  /* Add in scores reflecting rarity of the atoms in question */
  atom_ptr = equiv_atom_ptr->atom_ptr;
/* v.5.0--> */
/*  score = score + (max_equiv - atom_ptr->nequiv); */
  if (!strcmp(atom_ptr->element," C"))
    score = score + 1;
  else if (!strcmp(atom_ptr->element," O"))
    score = score + 3;
  else if (!strcmp(atom_ptr->element," N"))
    score = score + 6;
  else
    score = score + 10;
/* <--v.5.0 */

  atom_ptr = other_equiv_atom_ptr->atom_ptr;
/* v.5.0--> */
/*  score = score + (max_equiv - atom_ptr->nequiv); */
  if (!strcmp(atom_ptr->element," C"))
    score = score + 1;
  else if (!strcmp(atom_ptr->element," O"))
    score = score + 3;
  else if (!strcmp(atom_ptr->element," N"))
    score = score + 6;
  else
    score = score + 10;
/* <--v.5.0 */

  /* If this is the first equivalent bond for the current equiv_atom,
     then update pointer to head of its linked list of links */
  if (equiv_atom_ptr->first_equiv_bond_ptr == NULL)
    {
      equiv_atom_ptr->first_equiv_bond_ptr = equiv_bond_ptr;
      equiv_bond_ptr->ea_next_equiv_bond_ptr = NULL;
    }

  /* Otherwise, place it in the current list in accordance with
     its match score */
  else
    {
      /* Initialise pointers to new location */
      before_equiv_bond_ptr = NULL;
      after_equiv_bond_ptr = equiv_atom_ptr->first_equiv_bond_ptr;

      /* Start at the first equiv_bond record */
      other_equiv_bond_ptr = equiv_atom_ptr->first_equiv_bond_ptr;

      /* Get its score */
      other_score = other_equiv_bond_ptr->match_score;

      /* Loop until hit end of list, or find the correct slot to
         put current equivalenced bond */
      while (other_equiv_bond_ptr != NULL && score < other_score)
        {
          /* Amend the location pointers */
          before_equiv_bond_ptr = other_equiv_bond_ptr;
          after_equiv_bond_ptr
            = other_equiv_bond_ptr->ea_next_equiv_bond_ptr;

          /* Go to next one in the list */
          other_equiv_bond_ptr = after_equiv_bond_ptr;

          /* Get its score */
          if (other_equiv_bond_ptr != NULL)
            other_score = other_equiv_bond_ptr->match_score;
          else
            other_score = 0;
        }

      /* Place the equivalent bond in its appropriate location */

      /* If it is to be first, then adjust pointers */
      if (before_equiv_bond_ptr == NULL)
        equiv_atom_ptr->first_equiv_bond_ptr = equiv_bond_ptr;

      /* Otherwise, get the previous entry to point to the current one */
      else
        before_equiv_bond_ptr->ea_next_equiv_bond_ptr = equiv_bond_ptr;

      /* Point to previous first in the list */
      equiv_bond_ptr->ea_next_equiv_bond_ptr = after_equiv_bond_ptr;
    }

  /* Store current equivalent bond as the current equiv_atom's last
     equivalent bond and increment count of links */
  last_equiv_bond_ptr = equiv_bond_ptr;
  equiv_atom_ptr->nlinks++;

  /* Initialise other fields */
  equiv_bond_ptr->possible_start = TRUE;
  equiv_bond_ptr->checked = FALSE;
  equiv_bond_ptr->on_stack = FALSE;
  equiv_bond_ptr->link_number = *nequiv_bonds;

  /* Store the match score */
  equiv_bond_ptr->match_score = score;

  /* Store the link details */
  equiv_bond_ptr->equiv_atom_ptr = equiv_atom_ptr;
  equiv_bond_ptr->other_equiv_atom_ptr = other_equiv_atom_ptr;
  equiv_bond_ptr->reverse_equiv_bond_ptr = NULL;
  equiv_bond_ptr->next_equiv_bond_ptr = NULL;
  (*nequiv_bonds)++;

  /* Return first and last pointers */
  *fst_equiv_bond_ptr = first_equiv_bond_ptr;
  *lst_equiv_bond_ptr = last_equiv_bond_ptr;

  /* Return the new pointer */
  return(equiv_bond_ptr);
}
/***********************************************************************

calc_equiv_bonds  -  Calculate all the possible links between possibly
                     equivalenced atoms

***********************************************************************/

int calc_equiv_bonds(struct equiv_atom *first_equiv_atom_ptr,
                     struct equiv_bond **fst_equiv_bond_ptr,
                     int max_equiv)
{
  int done, have_bond;
  int iconnect;
  int nequiv_bonds;

  struct atom *atom_ptr, *other_atom_ptr;
  struct dic_data *dic_data_ptr, *other_dic_data_ptr;
  struct equiv_atom *equiv_atom_ptr, *other_equiv_atom_ptr;
  struct equiv_bond *first_equiv_bond_ptr, *last_equiv_bond_ptr;
  struct equiv_bond *equiv_bond1_ptr, *equiv_bond2_ptr;


  /* Initialise variables */
  last_equiv_bond_ptr = NULL;
  nequiv_bonds = 0;

  /* Initialise pointers */
  first_equiv_bond_ptr = *fst_equiv_bond_ptr;

  /* Get pointer to first equiv_atom */
  equiv_atom_ptr = first_equiv_atom_ptr;

  /* Loop through all the equiv_atoms */
  while (equiv_atom_ptr != NULL)
    {
      /* Get the atom pointers to the corresponding atoms */
      atom_ptr = equiv_atom_ptr->atom_ptr;
      dic_data_ptr = equiv_atom_ptr->dic_data_ptr;

      /* Get pointer to first equiv_atom */
      other_equiv_atom_ptr = equiv_atom_ptr->next_equiv_atom_ptr;

      /* Loop through all other equiv_atoms */
      while (other_equiv_atom_ptr != NULL)
        {
          /* Get the atom pointers to the corresponding atoms */
          other_atom_ptr = other_equiv_atom_ptr->atom_ptr;
          other_dic_data_ptr = other_equiv_atom_ptr->dic_data_ptr;

          /* Skip if any of the atoms in the two equiv_atoms are the
             same */
          if (other_atom_ptr != atom_ptr &&
              other_dic_data_ptr != dic_data_ptr)
            {
              /* Check whether have a bond between the two atoms in
                 the residue */
              have_bond = FALSE;

              /* Loop through first atom's connections to see if
                 other atom is amongst them */
              done = FALSE;
              for (iconnect = 0; iconnect < MAXCONNECT && done == FALSE;
                   iconnect++)
                {
                  /* Check for end of the connections */
                  if (atom_ptr->connect_atom_ptr[iconnect] == NULL)
                    done = TRUE;

                  /* Otherwise, check for the second atom */
                  else if (atom_ptr->connect_atom_ptr[iconnect]
                           == other_atom_ptr)
                    {
                      have_bond = TRUE;
                      done = TRUE;
                    }
                }

              /* If have a bond, then repeat for the two atoms from
                 the dictionary */
              if (have_bond == TRUE)
                {
                  /* Initialise flag */
                  have_bond = FALSE;

                  /* Loop through first dictionary atom's connections
                     to see if other dictionary atom is amongst them */
                  done = FALSE;
                  for (iconnect = 0; iconnect < MAXCONNECT && done == FALSE;
                       iconnect++)
                    {
                      /* Check for end of the connections */
                      if (dic_data_ptr->connect_atom_ptr[iconnect] == NULL)
                        done = TRUE;

                      /* Otherwise, check for the second atom */
                      else if (dic_data_ptr->connect_atom_ptr[iconnect]
                               == other_dic_data_ptr)
                        {
                          have_bond = TRUE;
                          done = TRUE;
                        }
                    }

                  /* If have a bond here, too, then can connect the
                     two equiv_atoms */
                  if (have_bond == TRUE)
                    {
                      /* Create a connection between the two possibly
                       equivalenced atoms */
                      equiv_bond1_ptr
                        = create_equiv_bond(equiv_atom_ptr,
                                            other_equiv_atom_ptr,
                                            &first_equiv_bond_ptr,
                                            &last_equiv_bond_ptr,
                                            &nequiv_bonds,max_equiv);

                      /* If debugging, show equivalence */
                      if (VERBOSE == TRUE)
                        printf("%d. [%s] - [%s]  ==  [%s] - [%s]\n",
                               nequiv_bonds,
                               equiv_atom_ptr->atom_ptr->atom_name,
                               other_equiv_atom_ptr->atom_ptr->atom_name,
                               equiv_atom_ptr->dic_data_ptr->atom_name,
                               other_equiv_atom_ptr->dic_data_ptr->atom_name);

                      /* Repeat to create a connection going in the opposite
                         direction */
                      equiv_bond2_ptr
                        = create_equiv_bond(other_equiv_atom_ptr,
                                            equiv_atom_ptr,
                                            &first_equiv_bond_ptr,
                                            &last_equiv_bond_ptr,
                                            &nequiv_bonds,max_equiv);

                      /* Store pointers across between these two links */
                      equiv_bond1_ptr->reverse_equiv_bond_ptr
                        = equiv_bond2_ptr;
                      equiv_bond2_ptr->reverse_equiv_bond_ptr
                        = equiv_bond1_ptr;
                    }
                }
            }

          /* Get pointer to the next equiv_atom */
          other_equiv_atom_ptr = other_equiv_atom_ptr->next_equiv_atom_ptr;
        }

      /* Get pointer to the next equiv_atom */
      equiv_atom_ptr = equiv_atom_ptr->next_equiv_atom_ptr;
    }

  /* Return pointer to first of the equivalenced bonds */
  *fst_equiv_bond_ptr = first_equiv_bond_ptr;

  /* Return the number of equivalenced bonds */
  return(nequiv_bonds);
}
/***********************************************************************

initialise_equiv_atoms  -  Initialise all the equiv_atom flags

***********************************************************************/

void initialise_equiv_atoms(struct equiv_atom *first_equiv_atom_ptr)
{
  struct equiv_atom *equiv_atom_ptr;

  /* Get pointer to the start of linked list of equiv_atom coords */
  equiv_atom_ptr = first_equiv_atom_ptr;

  /* Loop through all the equiv_atoms to initialise the checked flag */
  while (equiv_atom_ptr != NULL)
    {
      /* Initialise the flag */
      equiv_atom_ptr->checked = FALSE;
      
      /* Get pointer to next equiv_atom in linked-list */
      equiv_atom_ptr = equiv_atom_ptr->next_equiv_atom_ptr;
    }               
}
/***********************************************************************

check_perfect_match  -  Check whether we already have a perfect match
                        between the ligand atom names in the PDB file
                        and those in the Het Group Dictionary definition

***********************************************************************/

int check_perfect_match(struct equiv_bond *first_equiv_bond_ptr,
                        struct equiv_atom *first_equiv_atom_ptr,int *nb)
{
  char atom_name[5], other_atom_name[5];
  char dic_atom_name[5], other_dic_atom_name[5];

  int nbonds, nfailed, nmatched;
  int perfect;

  struct atom *atom_ptr, *other_atom_ptr;
  struct dic_data *dic_data_ptr, *other_dic_data_ptr;
  struct equiv_atom *equiv_atom_ptr, *other_equiv_atom_ptr;
  struct equiv_bond *equiv_bond_ptr;

  /* Initialise variables */
  nbonds = 0;
  nmatched = 0;
  nfailed = 0;
  perfect = FALSE;

  /* Initialise all the equivalenced atom records */
  initialise_equiv_atoms(first_equiv_atom_ptr);

  /* Get pointer to the start of linked list of equiv_bond coords */
  equiv_bond_ptr = first_equiv_bond_ptr;

  /* Loop through all the equiv_bonds to initialise the checked flag */
  while (equiv_bond_ptr != NULL)
    {
      /* Get the two equivalences at either end of this potential
         bond match */
      equiv_atom_ptr = equiv_bond_ptr->equiv_atom_ptr;
      other_equiv_atom_ptr = equiv_bond_ptr->other_equiv_atom_ptr;

      /* Get the atom and dictionary entries corresponding to these
         equivalences */
      atom_ptr = equiv_atom_ptr->atom_ptr;
      dic_data_ptr = equiv_atom_ptr->dic_data_ptr;
      other_atom_ptr = other_equiv_atom_ptr->atom_ptr;
      other_dic_data_ptr = other_equiv_atom_ptr->dic_data_ptr;

      /* Get all the corresponding atom names */
      strcpy(atom_name,atom_ptr->atom_name);
      strcpy(dic_atom_name,dic_data_ptr->atom_name);
      strcpy(other_atom_name,other_atom_ptr->atom_name);
      strcpy(other_dic_atom_name,other_dic_data_ptr->atom_name);

      /* Check for an exact match */
      if (!strcmp(atom_name,dic_atom_name) &&
          !strcmp(other_atom_name,other_dic_atom_name))
        {
          /* Mark the equivalenced atom records as matched */
          equiv_atom_ptr->checked = TRUE;
          other_equiv_atom_ptr->checked = TRUE;

          /* Increment count of perfectly matching bonds */
          nbonds++;
        }

      /* Get pointer to next equiv_bond in linked-list */
      equiv_bond_ptr = equiv_bond_ptr->next_equiv_bond_ptr;
    }

  /* Get pointer to the first equiv_atom record */
  equiv_atom_ptr = first_equiv_atom_ptr;

  /* Loop through all the equiv_atoms to determine how many identically-
     named atoms (between PDB file and dictionary) have been checked
     using the connectivity information in the PDB file and dictionary
     entry */
  while (equiv_atom_ptr != NULL)
    {
      /* Get the atom names for this equivalence */
      atom_ptr = equiv_atom_ptr->atom_ptr;
      dic_data_ptr = equiv_atom_ptr->dic_data_ptr;
      strcpy(atom_name,atom_ptr->atom_name);
      strcpy(dic_atom_name,dic_data_ptr->atom_name);

      /* If atom names are identical, see if this equivalence is
         compatible with the connectivity info */
      if (!strcmp(atom_name,dic_atom_name))
        {
          /* Check the flag */
          if (equiv_atom_ptr->checked == TRUE)
            nmatched++;
          else
            nfailed++;
        }

      /* Get pointer to next equiv_atom in linked-list */
      equiv_atom_ptr = equiv_atom_ptr->next_equiv_atom_ptr;
    }               

  /* Get number of matches and failures */
  if (nmatched > 0 && nfailed == 0)
    {
      /* Set flag to indicate perfect match achieved */
      perfect = TRUE;
    }

  /* Return number of bonds matched */
  *nb = nbonds / 2;

  /* Return whether have a perfect match */
  return(perfect);
}
/***********************************************************************

get_best_possible_score  -  Calculate the best possible score that can
                            be achieved in matching the current residue
                            against its dictionary entry

***********************************************************************/

int get_best_possible_score(struct residue *residue_ptr,
                            int *perfect_possible)
{
  int best_score, iatom, jatom, natoms, ndic_atoms, nlink, nmatch;
  int got_match;

  struct atom *atom_ptr;
  struct dic_data *dic_data_ptr;
  struct dictionary *dictionary_ptr;

  /* Initialise best possible score */
  nlink = 0;
  nmatch = 0;
  *perfect_possible = FALSE;

  /* Get pointer to the start of linked list of atom coords */
  atom_ptr = residue_ptr->first_atom_ptr;
  natoms = residue_ptr->natoms;

  /* Get the residue's dictionary entry */
  dictionary_ptr = residue_ptr->dictionary_ptr;

  /* Loop through all the atoms */
  for (iatom = 0; iatom < natoms && atom_ptr != NULL; iatom++)
    {
      /* Get pointer to the start of linked list of dictionary atoms */
      dic_data_ptr = dictionary_ptr->first_dic_data_ptr;
      ndic_atoms = dictionary_ptr->natoms;

      /* Initialise flag */
      got_match = FALSE;

      /* If not a hydrogen, increment count */
      if (strcmp(atom_ptr->element," H"))
        nlink++;

      /* Loop through all the dictionary atoms */
      for (jatom = 0; jatom < ndic_atoms && dic_data_ptr != NULL &&
           got_match == FALSE; jatom++)
        {
          /* Check if the atom name matches this dictionary atom's name */
          if (!strncmp(atom_ptr->atom_name,dic_data_ptr->atom_name,4))
            {
              /* Verify that the elements also match */
              if (strcmp(atom_ptr->element,dic_data_ptr->element))
                {
                  /* Print warning */
                  printf("*** Warning. PDB and dictionary atom names");
                  printf(" match [%s] for [%s %s %c]\n",
                         atom_ptr->atom_name,residue_ptr->res_name,
                         residue_ptr->res_num,residue_ptr->chain);
                  printf("***          But element assignments do not: ");
                  printf("[%s] vs [%s]\n",atom_ptr->element,
                         dic_data_ptr->element);
                }

              /* Set flag to indicate we've found a match */
              got_match = TRUE;

              /* If atoms are not hydrogens, then increase count of
                 potentially matching atoms */
              if (strcmp(atom_ptr->element," H"))
                nmatch++;
            }

          /* Get pointer to next dic_data in linked-list */
          dic_data_ptr = dic_data_ptr->next_dic_data_ptr;
        }               
      
      /* Get pointer to next atom in linked-list */
      atom_ptr = atom_ptr->next_atom_ptr;
    }

  /* If have as many possible matches as there are atoms, then a
     perfect match may be possible */
  if (nlink == nmatch)
    *perfect_possible = TRUE;

  /* Add maximum score obtainable from the best possible match of the
     atoms */
  if (dictionary_ptr->non_hydrogen < nlink)
    nlink = dictionary_ptr->non_hydrogen;
/* v.5.0--> */
/*  best_score = nmatch * MATCH_SCORE + (nlink - nmatch) * MAP_SCORE; */
  best_score = nmatch * Match_Score + (nlink - nmatch) * MAP_SCORE;
/* <--v.5.0 */

  /* Return the computed best possible score */
  return(best_score);
}
/***********************************************************************

initialise_atoms  -  Initialise all the atom flags for the given residue

***********************************************************************/

void initialise_atoms(struct residue *residue_ptr)
{
  int iatom, natoms;

  struct atom *atom_ptr;

  /* Get pointer to the start of linked list of atom coords */
  atom_ptr = residue_ptr->first_atom_ptr;
  natoms = residue_ptr->natoms;

  /* Loop through all the atoms to initialise the checked flag */
  for (iatom = 0; iatom < natoms && atom_ptr != NULL; iatom++)
    {
      /* Initialise the flag */
      atom_ptr->checked = FALSE;

      /* Blank out the pointer to the dictionary entry */
      atom_ptr->dic_data_ptr = NULL;
      
      /* Get pointer to next atom in linked-list */
      atom_ptr = atom_ptr->next_atom_ptr;
    }               
}
/***********************************************************************

initialise_dic_data  -  Initialise all the dic_data flags for the given
                        dictionary entry

***********************************************************************/

void initialise_dic_data(struct dictionary *dictionary_ptr)
{
  int iatom, natoms;

  struct dic_data *dic_data_ptr;

  /* Get pointer to the start of linked list of dic_data coords */
  dic_data_ptr = dictionary_ptr->first_dic_data_ptr;
  natoms = dictionary_ptr->natoms;

  /* Loop through all the dic_datas to initialise the checked flag */
  for (iatom = 0; iatom < natoms && dic_data_ptr != NULL; iatom++)
    {
      /* Initialise the flag */
      dic_data_ptr->checked = FALSE;
      
      /* Blank out the pointer to the dictionary entry */
      dic_data_ptr->atom_ptr = NULL;

      /* Get pointer to next dic_data in linked-list */
      dic_data_ptr = dic_data_ptr->next_dic_data_ptr;
    }               
}
/***********************************************************************

count_atom_matches  -  Loop through all the equiv_atoms in the current clique
                       and count the number of times that an atom-atom
                       mapping, represented by a equiv_atom, gives two atoms
                       with identical names

***********************************************************************/

void count_atom_matches(struct equiv_atom *first_equiv_atom_ptr,int *nmatch_atoms)
{
  struct equiv_atom *equiv_atom_ptr;

  /* Get pointer to the start of linked list of equiv_atom coords */
  equiv_atom_ptr = first_equiv_atom_ptr;

  /* Loop through all the equiv_atoms to check the mappings */
  while (equiv_atom_ptr != NULL)
    {
      /* Check only if equiv_atom belongs to the current clique */
      if (equiv_atom_ptr->checked == TRUE)
        {
          /* If the atom-names are identical, increase the score */
          if (!strncmp(equiv_atom_ptr->atom_ptr->atom_name,
                       equiv_atom_ptr->dic_data_ptr->atom_name,4))
            (*nmatch_atoms)++;
        }

      /* Get pointer to next equiv_atom in linked-list */
      equiv_atom_ptr = equiv_atom_ptr->next_equiv_atom_ptr;
    }               
}
/***********************************************************************

initialise_equiv_bonds  -  Initialise all the equiv_bond flags

***********************************************************************/

void initialise_equiv_bonds(struct equiv_bond *first_equiv_bond_ptr)
{
  struct equiv_bond *equiv_bond_ptr;

  /* Get pointer to the start of linked list of equiv_bond coords */
  equiv_bond_ptr = first_equiv_bond_ptr;

  /* Loop through all the equiv_bonds to initialise the checked flag */
  while (equiv_bond_ptr != NULL)
    {
      /* Initialise the flags */
      equiv_bond_ptr->checked = FALSE;
      equiv_bond_ptr->on_stack = FALSE;

      /* Get pointer to next equiv_bond in linked-list */
      equiv_bond_ptr = equiv_bond_ptr->next_equiv_bond_ptr;
    }               
}
/***********************************************************************

create_tree_entry  -  Create a new entry in the tree of possible paths

***********************************************************************/

struct tree *create_tree_entry(struct equiv_bond *equiv_bond_ptr,
                               struct tree *prev_level_tree_ptr,
                               struct tree *prev_tree_ptr,
                               struct tree **fst_all_tree_ptr,
                               struct tree **lst_all_tree_ptr)
{
  struct tree *tree_ptr, *first_all_tree_ptr, *last_all_tree_ptr;

  /* Initialise pointers */
  first_all_tree_ptr = *fst_all_tree_ptr;
  last_all_tree_ptr = *lst_all_tree_ptr;

  /* Allocate memory for structure to hold tree info */
  tree_ptr = (struct tree *) malloc(sizeof(struct tree));
  if (tree_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for struct tree\n");
      printf("***        Program hbadd terminated with error.\n");
      exit (1);
    }

  /* If have a level above, then update its entry */
  if (prev_level_tree_ptr != NULL)
    {
      /* If this is the first entry from previous level, then update
         pointer down from previous level */
      if (prev_level_tree_ptr->next_level_tree_ptr == NULL)
        prev_level_tree_ptr->next_level_tree_ptr = tree_ptr;

      /* Otherwise, get the previous equivalent bond at this level and get
         it to point to the current equivalent bond */
      else
        {
          /* Update its pointer from previous equivalent bond */
          prev_tree_ptr->next_tree_ptr = tree_ptr;
        }
    }

  /* Update links across all tree entries */
  if (first_all_tree_ptr == NULL)
    first_all_tree_ptr = tree_ptr;
  else
    last_all_tree_ptr->all_next_tree_ptr = tree_ptr;
  last_all_tree_ptr = tree_ptr;

  /* Initialise all fields */
  tree_ptr->branches_checked = FALSE;
  tree_ptr->score = 0;

  /* Initialise pointers */
  tree_ptr->equiv_bond_ptr = equiv_bond_ptr;
  tree_ptr->prev_level_tree_ptr = prev_level_tree_ptr;
  tree_ptr->next_level_tree_ptr = NULL;
  tree_ptr->link_tree_ptr = NULL;
  tree_ptr->first_end_tree_ptr = NULL;
  tree_ptr->next_end_tree_ptr = NULL;
  tree_ptr->all_next_tree_ptr = NULL;

  /* Initialise forward and backward pointers */
  tree_ptr->prev_tree_ptr = prev_tree_ptr;
  tree_ptr->next_tree_ptr = NULL;

  /* Flag that bond equivalence is on tree-stack */
  equiv_bond_ptr->on_stack = TRUE;
  equiv_bond_ptr->reverse_equiv_bond_ptr->on_stack = TRUE;

  /* Return pointer to last tree entry created */
  *fst_all_tree_ptr = first_all_tree_ptr;
  *lst_all_tree_ptr = last_all_tree_ptr;

  /* Return the new pointer */
  return(tree_ptr);
}
/***********************************************************************

check_compatibility  -  Check whether the given equivalenced bond is
                        compatible with all on the path up to this point

***********************************************************************/

int check_compatibility(struct equiv_bond *current_equiv_bond_ptr,
                        struct equiv_bond *equiv_bond_ptr,
                        struct equiv_atom *first_equiv_atom_ptr)
{
  int wanted;

  struct atom *atom_ptr;
  struct dic_data *dic_data_ptr;
  struct equiv_atom *equiv_atom_ptr, *other_equiv_atom_ptr;

  /* Initialise variables */
  wanted = TRUE;

  /* Check that not the same bond and not already in the current path */
  if (current_equiv_bond_ptr != NULL)
    {
      if (equiv_bond_ptr == current_equiv_bond_ptr ||
          equiv_bond_ptr == current_equiv_bond_ptr->reverse_equiv_bond_ptr ||
          equiv_bond_ptr->on_stack == TRUE)
        wanted = FALSE;
    }

  /* Otherwise, check compatibility of atom equivalences */
  if (wanted == TRUE)
    {
      /* Get atom at other end of this bond and its corresponding
         dictionary entry */
      equiv_atom_ptr = equiv_bond_ptr->other_equiv_atom_ptr;
      atom_ptr = equiv_atom_ptr->atom_ptr;
      dic_data_ptr = equiv_atom_ptr->dic_data_ptr;

      /* Get pointer to first atom equivalence */
      other_equiv_atom_ptr = first_equiv_atom_ptr;

      /* Loop through all atom equivalences to check all those in the
         current path */
      while (other_equiv_atom_ptr != NULL && wanted == TRUE)
        {
          /* If thise equivalence is in the current path, then check
             it out */
          if (other_equiv_atom_ptr->checked == TRUE)
            {
              /* Check for compatibility */
              if ((atom_ptr == other_equiv_atom_ptr->atom_ptr &&
                   dic_data_ptr != other_equiv_atom_ptr->dic_data_ptr) ||
                  (atom_ptr != other_equiv_atom_ptr->atom_ptr &&
                   dic_data_ptr == other_equiv_atom_ptr->dic_data_ptr))
                wanted = FALSE;
            }

          /* Get pointer to next atom, equivalence */
          other_equiv_atom_ptr = other_equiv_atom_ptr->next_equiv_atom_ptr;
        }
    }

  /* Return whether bond is compatible or not */
  return(wanted);
}
/***********************************************************************

get_score  -  Get score for the given path

***********************************************************************/

int get_score(struct equiv_atom *first_equiv_atom_ptr,int *nmch,
                int *nmap)
{
  char atom_name[5], dic_atom_name[5];

  int abort;
  int nmatched, nmapped, score;

  struct atom *atom_ptr;
  struct dic_data *dic_data_ptr;
  struct equiv_atom *equiv_atom_ptr, *other_equiv_atom_ptr;

  /* Initialise score */
  abort = FALSE;
  nmatched = 0;
  nmapped = 0;
  score = 0;
  if (VERBOSE == TRUE)
    printf("Path:");

  /* Get pointer to the start of linked list of equiv_atom coords */
  equiv_atom_ptr = first_equiv_atom_ptr;

  /* Loop through all the equiv_atoms to build up total score */
  while (equiv_atom_ptr != NULL && abort == FALSE)
    {
      /* If this atom-equivalence is on the current path, then
         include in scoring function */
      if (equiv_atom_ptr->checked == TRUE)
        {
          /* Get atom at other end of this bond and its corresponding
             dictionary entry */
          atom_ptr = equiv_atom_ptr->atom_ptr;
          dic_data_ptr = equiv_atom_ptr->dic_data_ptr;

          /* Get the PDB and dictionary atom names */
          strcpy(atom_name,atom_ptr->atom_name);
          strcpy(dic_atom_name,dic_data_ptr->atom_name);

          /* If atom-names match exactly, then add match-score */
          if (!strcmp(atom_name,dic_atom_name))
            nmatched++;

          /* Otherwise, add map-score */
          else
            nmapped++;

          /* Show current mapping */
          if (VERBOSE == TRUE)
            printf(" [%s/%s]",atom_name,dic_atom_name);

          /* Check that not incompatible with any existing atom
             equivalence */
          other_equiv_atom_ptr = first_equiv_atom_ptr;

          /* Loop through all the equiv_atoms to check compatibility */
          while (other_equiv_atom_ptr != NULL && abort == FALSE)
            {
              /* If this atom-equivalence is on the current path, then
                 include in scoring function */
              if (other_equiv_atom_ptr->checked == TRUE)
                {
                  /* Check for compatibility */
                  if ((atom_ptr == other_equiv_atom_ptr->atom_ptr &&
                       dic_data_ptr != other_equiv_atom_ptr->dic_data_ptr) ||
                      (atom_ptr != other_equiv_atom_ptr->atom_ptr &&
                       dic_data_ptr == other_equiv_atom_ptr->dic_data_ptr))
                    {
                      /* Abort scoring of this path and set its score
                         to zero */
                      abort = TRUE;
                      nmatched = 0;
                      nmapped = 0;
                    }
                }

              /* Get pointer to next equiv_atom in linked-list */
              other_equiv_atom_ptr
                = other_equiv_atom_ptr->next_equiv_atom_ptr;
            }
        }

      /* Get pointer to next equiv_atom in linked-list */
      equiv_atom_ptr = equiv_atom_ptr->next_equiv_atom_ptr;
    }

  /* Calculate score */
/* v.5.0--> */
/*  score = nmatched * MATCH_SCORE + nmapped * MAP_SCORE; */
  score = nmatched * Match_Score + nmapped * MAP_SCORE;
/* <--v.5.0 */
  if (VERBOSE == TRUE)
    printf("  Score = %d",score);

  /* Return the match counts */
  *nmch = nmatched;
  *nmap = nmapped;

  /* Return the score */
  return(score);
}
/***********************************************************************

get_bond_score  -  Calculate score from the equivalenced bonds

***********************************************************************/

int get_bond_score(struct equiv_bond *first_equiv_bond_ptr)
{
  int score;

  struct equiv_bond *equiv_bond_ptr;

  /* Initialise score */
  score = 0;

  /* Get pointer to the start of linked list of equiv_bond coords */
  equiv_bond_ptr = first_equiv_bond_ptr;

  /* Loop through all the equiv_bonds to count the number of the
         current stack */
  while (equiv_bond_ptr != NULL)
    {
      /* If bond on stack, increment score */
      if (equiv_bond_ptr->on_stack == TRUE)
                  score++;

      /* Get pointer to next equiv_bond in linked-list */
      equiv_bond_ptr = equiv_bond_ptr->next_equiv_bond_ptr;
    }               

  /* Return bond score */
  return(score);
}
/***********************************************************************

reset_path  -  Re-initialise all the flags on the atom along the current
               path

***********************************************************************/

void reset_path(struct tree *first_tree_ptr,
                struct equiv_atom *first_equiv_atom_ptr,
                struct equiv_bond *first_equiv_bond_ptr,int trace)
{
  struct equiv_bond *equiv_bond_ptr;
  struct tree *tree_ptr;

  /* Re-set all flags on equivalenced atoms */
  initialise_equiv_atoms(first_equiv_atom_ptr);
  initialise_equiv_bonds(first_equiv_bond_ptr);

  /* Set pointer to last entry on path */
  tree_ptr = first_tree_ptr;
  if (trace == TRUE)
    {
      printf("\n");
      printf("Tracing path:\n");
    }

  /* Loop through current path */
  while (tree_ptr != NULL)
    {
      /* Get the equivalenced bond-pair represented by this path
         branch */
      equiv_bond_ptr = tree_ptr->equiv_bond_ptr;
      if (trace == TRUE)
        printf("[%s->%s/%s->%s]\n",
               equiv_bond_ptr->equiv_atom_ptr->atom_ptr->atom_name,
               equiv_bond_ptr->other_equiv_atom_ptr->atom_ptr->atom_name,
               equiv_bond_ptr->equiv_atom_ptr->dic_data_ptr->atom_name,
               equiv_bond_ptr->other_equiv_atom_ptr->dic_data_ptr->atom_name);

      /* Flag this equivalenced bond as on the stack */
      equiv_bond_ptr->on_stack = TRUE;
      equiv_bond_ptr->reverse_equiv_bond_ptr->on_stack = TRUE;

      /* Set the flags of both equivalenced atoms */
      equiv_bond_ptr->equiv_atom_ptr->checked = TRUE;
      equiv_bond_ptr->other_equiv_atom_ptr->checked = TRUE;

      /* If path has a link, follow that */
      if (tree_ptr->link_tree_ptr != NULL)
        tree_ptr = tree_ptr->link_tree_ptr;

      /* Otherwise, go up a level to get next entry in the path */
      else
        tree_ptr = tree_ptr->prev_level_tree_ptr;
    }
}
/***********************************************************************

check_branches  -  Check which, if any, of the previous branches are
                   consistent with the current one

***********************************************************************/

void check_branches(struct tree *path_tree_ptr,
                    struct equiv_atom *first_equiv_atom_ptr)
{
  int best_score, score;
  int compatible;

  struct equiv_atom *equiv_atom_ptr, *other_equiv_atom_ptr;
  struct equiv_bond *equiv_bond_ptr;
  struct tree *prev_level_tree_ptr;
  struct tree *end_tree_ptr, *path_end_tree_ptr, *trace_path_tree_ptr;

  /* Initialise variables */
  best_score = 0;

  /* Initialise pointers */
  path_end_tree_ptr = path_tree_ptr->first_end_tree_ptr;

  /* Get the pointer to the first of the path-ends to search */
  prev_level_tree_ptr = path_tree_ptr->prev_level_tree_ptr;
  end_tree_ptr = prev_level_tree_ptr->first_end_tree_ptr;

  /* Loop through all the paths leading to this branch */
  while (end_tree_ptr != NULL && end_tree_ptr != path_end_tree_ptr)
    {
      /* Process only if the score of this path is non-zero */
      if (end_tree_ptr->score > 0)
        {
          /* Get the equivalent bond */
          equiv_bond_ptr = end_tree_ptr->equiv_bond_ptr;

          /* Set flags */
          compatible = TRUE;
          score = 0;

          /* Get the path-end */
          trace_path_tree_ptr = end_tree_ptr;

          /* Loop through this path up to the branch-point */
          while (trace_path_tree_ptr != NULL && compatible == TRUE &&
                 trace_path_tree_ptr != prev_level_tree_ptr)
            {
              /* Get the equivalent bond */
              equiv_bond_ptr = trace_path_tree_ptr->equiv_bond_ptr;

              /* Get the two equivalenced atoms */
              equiv_atom_ptr = equiv_bond_ptr->equiv_atom_ptr;
              other_equiv_atom_ptr = equiv_bond_ptr->other_equiv_atom_ptr;

              /* Get the scores for the match */
              if (!strcmp(equiv_atom_ptr->atom_ptr->atom_name,
                          equiv_atom_ptr->dic_data_ptr->atom_name))
/* v.5.0--> */
/*                score = score + MATCH_SCORE; */
                score = score + Match_Score;
/* <--v.5.0 */
              else
                score = score + MAP_SCORE;
              if (!strcmp(other_equiv_atom_ptr->atom_ptr->atom_name,
                          other_equiv_atom_ptr->dic_data_ptr->atom_name))
/* v.5.0--> */
/*                score = score + MATCH_SCORE; */
                score = score + Match_Score;
/* <--v.5.0 */
              else
                score = score + MAP_SCORE;

              /* Check compatibility */
              compatible
                = check_compatibility(NULL,equiv_bond_ptr,
                                      first_equiv_atom_ptr);

              /* Get next tree entry along this path */
              if (trace_path_tree_ptr->link_tree_ptr != NULL)
                trace_path_tree_ptr = trace_path_tree_ptr->link_tree_ptr;
              else
                trace_path_tree_ptr = trace_path_tree_ptr->prev_level_tree_ptr;
            }

          /* If both paths compatible, and this is the largest score so
             far, save path-end as possible link-point */
          if (compatible == TRUE && score > best_score)
            {
              /* Link to this path-end */
              path_tree_ptr->link_tree_ptr = end_tree_ptr;

              /* Save the best score so far */
              best_score = score;
            }
        }

      /* Get next path-end pointer */
      end_tree_ptr = end_tree_ptr->next_end_tree_ptr;
    }

  /* Flag current entry as having had its branches checked */
  path_tree_ptr->branches_checked = TRUE;
}
/***********************************************************************

traverse_path  -  Traverse the given path to link in any compatible
                  branches

***********************************************************************/

void traverse_path(struct tree *end_tree_ptr,
                   struct equiv_atom *first_equiv_atom_ptr,
                   struct equiv_bond *first_equiv_bond_ptr)
{
/*  struct equiv_bond *equiv_bond_ptr, *prev_equiv_bond_ptr; */
  struct tree *tree_ptr, *path_tree_ptr;

  /* Start that the first of the end-path pointers */
  tree_ptr = end_tree_ptr;

  /* Re-initialise all equivalenced bonds and atoms */
  reset_path(tree_ptr,first_equiv_atom_ptr,first_equiv_bond_ptr,
             FALSE);

  /* Store path-end pointer */
  path_tree_ptr = tree_ptr;

  /* Loop while chasing up the path */
  while (path_tree_ptr != NULL)
    {
      /* If this is the first path-end that leads to this node,
         store its pointer */
      if (path_tree_ptr->first_end_tree_ptr == NULL)
        path_tree_ptr->first_end_tree_ptr = tree_ptr;

      /* Get the equivalenced bond represented by this path-node */
/*      equiv_bond_ptr = path_tree_ptr->equiv_bond_ptr; */

      /* Check whether have any branches off to the left */
      if (path_tree_ptr->branches_checked == FALSE &&
          path_tree_ptr->prev_tree_ptr != NULL)
        {
          /* Get the previous branch */
/*          prev_tree_ptr = path_tree_ptr->prev_tree_ptr; */

          /* Get the corresponding equivalenced bond */
/*          prev_equiv_bond_ptr = prev_tree_ptr->equiv_bond_ptr; */

          /* Check which, if any, of the previous branches are
             consistent with the current one */
          check_branches(path_tree_ptr,first_equiv_atom_ptr);
        }

      /* Go to the next step in this path */
      if (path_tree_ptr->link_tree_ptr != NULL)
        path_tree_ptr = path_tree_ptr->link_tree_ptr;
      else
        path_tree_ptr = path_tree_ptr->prev_level_tree_ptr;
    }
}
/***********************************************************************

score_path  -  Calculate scores for the given path

***********************************************************************/

int score_path(struct tree *tree_ptr,struct tree **bst_tree_ptr,
               struct equiv_atom *first_equiv_atom_ptr,
               struct equiv_bond *first_equiv_bond_ptr,
               int best_score)
{
  int bond_score, nmapped, nmatched, score;

  struct tree *best_tree_ptr;

  /* Initialise pointers */
  best_tree_ptr = *bst_tree_ptr;

  /* Re-initialise all equivalenced bonds and atoms */
  reset_path(tree_ptr,first_equiv_atom_ptr,first_equiv_bond_ptr,
             FALSE);

  /* Calculate the path-score */
  score = get_score(first_equiv_atom_ptr,&nmatched,&nmapped);

  /* Calculate score for the equivalenced bonds */
  bond_score = get_bond_score(first_equiv_bond_ptr);
  if (VERBOSE == TRUE)
    printf("  Bond score = %d\n",bond_score);

  /* Store this path's score */
/* v.5.0--> */
  score = score + bond_score;
/* <--v.5.0 */
  tree_ptr->score = score;

  /* If this is the best score so far, then store */
  if (score > best_score)
    {
      /* Store score and pointer to end of path */
      best_score = score;
      best_tree_ptr = tree_ptr;
  }

  /* Return start of best path */
  *bst_tree_ptr = best_tree_ptr;

  /* Return the best score encountered */
  return(best_score);
}
/***********************************************************************

delete_branch  -  Delete brother of upper level node if equals node
                  just created

***********************************************************************/

void delete_branch(struct equiv_bond *equiv_bond_ptr,
                   struct tree *tree_ptr)
{
  int done;

  struct equiv_bond *other_equiv_bond_ptr;
  struct tree *other_tree_ptr, *prev_tree_ptr;

  /* Initialise variables */
  done = FALSE;

  /* Get first brother of the given tree entry */
  other_tree_ptr = tree_ptr->next_tree_ptr;

  /* Loop through all the brothers until done */
  while (other_tree_ptr != NULL && done == FALSE)
    {
      /* Get this entry's equivalenced bond */
      other_equiv_bond_ptr = other_tree_ptr->equiv_bond_ptr;

      /* If have a match, then delete the entry */
      if (other_equiv_bond_ptr == equiv_bond_ptr)
        {
          /* Get the tree entry to the left of this one */
          prev_tree_ptr = other_tree_ptr->prev_tree_ptr;

          /* Get it to by-pass current entry */
          prev_tree_ptr->next_tree_ptr = other_tree_ptr->next_tree_ptr;

          /* Amend reverse link, too */
          if (other_tree_ptr->next_tree_ptr != NULL)
            other_tree_ptr->next_tree_ptr->prev_tree_ptr
              = prev_tree_ptr;

          /* Set flag as done */
          done = TRUE;
        }

      /* Get the next tree entry */
      if (done == FALSE)
        other_tree_ptr = other_tree_ptr->next_tree_ptr;
    }
}
/***********************************************************************

process_tree  -  Process the entries on the equivalent bonds tree-stack
                 until all possible paths have been traced

***********************************************************************/

int process_tree(struct tree *first_tree_ptr,
                 struct equiv_atom *first_equiv_atom_ptr,
                 struct equiv_bond *first_equiv_bond_ptr,
                 struct tree **fst_end_tree_ptr,
                 struct tree **lst_end_tree_ptr,int *npth,
                 struct tree **fst_all_tree_ptr,
                 struct tree **lst_all_tree_ptr,
                 int best_possible,int *bst_score,
                 struct tree **bst_tree_ptr)
{
  int best_score, npaths;
  int iatom, np, ntree;
  int abort, done, go_down, got_next, wanted;
  static int nabort = 0;

/*  struct atom *atom_ptr, *atom_pair_ptr[2]; */
/*  struct dic_data *dic_data_pair_ptr[2]; */
  struct equiv_atom *equiv_atom_ptr;
  struct equiv_atom *equiv_atom_pair_ptr[2];
  struct equiv_bond *equiv_bond_ptr, *other_equiv_bond_ptr;
  struct tree *tree_ptr, *best_tree_ptr, *new_tree_ptr, *prev_tree_ptr;
  struct tree *first_end_tree_ptr, *last_end_tree_ptr;
  struct tree *first_all_tree_ptr, *last_all_tree_ptr;

  /* Initialise variables */
  abort = FALSE;
  best_score = *bst_score;
  done = FALSE;
  np = 0;
  ntree = 0;
  npaths = *npth;

  /* Get pointer to starting entry in the tree */
  tree_ptr = first_tree_ptr;

  /* Initalise all other pointers */
  prev_tree_ptr = NULL;
  best_tree_ptr = *bst_tree_ptr;
  first_end_tree_ptr = *fst_end_tree_ptr;
  last_end_tree_ptr = *lst_end_tree_ptr;
  first_all_tree_ptr = *fst_all_tree_ptr;
  last_all_tree_ptr = *lst_all_tree_ptr;

  /* Keep path-searching until no more paths to be found from the
     current stack */
  while (done == FALSE)
    {
      /* Initialise flags */
      go_down = FALSE;

      /* Get the equivalenced bond pointer to be processed */
      equiv_bond_ptr = tree_ptr->equiv_bond_ptr;
      ntree++;

      /* Unset start-flag on this bond */
      equiv_bond_ptr->possible_start = FALSE;
      equiv_bond_ptr->reverse_equiv_bond_ptr->possible_start = FALSE;

      /* Set current bond as checked */
      equiv_bond_ptr->checked = TRUE;
      equiv_bond_ptr->reverse_equiv_bond_ptr->checked = TRUE;

      /* Get the two atom-pairs */
      equiv_atom_ptr = equiv_bond_ptr->equiv_atom_ptr;
      equiv_atom_pair_ptr[0] = equiv_atom_ptr;
/*      atom_pair_ptr[0] = equiv_atom_ptr->atom_ptr; */
/*      dic_data_pair_ptr[0] = equiv_atom_ptr->dic_data_ptr; */
      equiv_atom_ptr = equiv_bond_ptr->other_equiv_atom_ptr;
      equiv_atom_pair_ptr[1] = equiv_atom_ptr;
/*      atom_pair_ptr[1] = equiv_atom_ptr->atom_ptr; */
/*      dic_data_pair_ptr[1] = equiv_atom_pair_ptr[1]->dic_data_ptr; */

      /* Mark both atom equivalence pairs as checked */
      equiv_atom_pair_ptr[0]->checked = TRUE;
      equiv_atom_pair_ptr[1]->checked = TRUE;

      /* Loop over the two atom-pairs in the current bond-pair */
      for (iatom = 0; iatom < 2; iatom++)
        {
          /* Get this atom's pointer and equivalent dictionary entry */
/*          atom_ptr = atom_pair_ptr[iatom]; */
/*          dic_data_ptr = dic_data_pair_ptr[iatom]; */
          equiv_atom_ptr = equiv_atom_pair_ptr[iatom];

          /* Get pointer to the first of this atom's equivalent
             bonds */
          other_equiv_bond_ptr = equiv_atom_ptr->first_equiv_bond_ptr;

          /* Loop through all the adjacent bond equivalences to
             find those that are compatible with all those in the
             path to this point */
          while (other_equiv_bond_ptr != NULL)
            {
              /* Check the compatibility of the current equivalenced
                 bond with all those in the path up to this point */
              wanted = check_compatibility(equiv_bond_ptr,
                                           other_equiv_bond_ptr,
                                           first_equiv_atom_ptr);

              /* If compatibility checks OK, then add to tree */
              if (wanted == TRUE)
                {
                  /* Create an entry on the tree */
                  new_tree_ptr = create_tree_entry(other_equiv_bond_ptr,
                                                   tree_ptr,prev_tree_ptr,
                                                   &first_all_tree_ptr,
                                                   &last_all_tree_ptr);

                  /* Check to see if any brothers of the level above
                     are the same as the current bond-equivalence and
                     can thus be deleted */
                  delete_branch(other_equiv_bond_ptr,tree_ptr);

                  /* Save the pointer returned */
                  prev_tree_ptr = new_tree_ptr;

                  /* Set flag to indicate that can go down to the
                     next level */
                  go_down = TRUE;

                  /* Increment count of entries in the tree */
                  ntree++;
                }

              /* Go to next equivalenced bond */
              other_equiv_bond_ptr
                = other_equiv_bond_ptr->ea_next_equiv_bond_ptr;
            }
        }

      /* If have reached the last node on this branch-point, then
         calculate score for the current path */
      if (go_down == FALSE)
        {
          /* Add this end-point to the linked list of path-ends */
          if (first_end_tree_ptr == NULL)
            first_end_tree_ptr = tree_ptr;
          else
            last_end_tree_ptr->next_end_tree_ptr = tree_ptr;

          /* Save this as the last path-end pointer so far */
          last_end_tree_ptr = tree_ptr;

          /* Traverse this path to link in any branch-points */
          traverse_path(tree_ptr,first_equiv_atom_ptr,
                        first_equiv_bond_ptr);

          /* Score the current path */
          best_score
            = score_path(tree_ptr,&best_tree_ptr,first_equiv_atom_ptr,
                         first_equiv_bond_ptr,best_score);

          /* If this equals the best possible, then end here */
          if (best_score >= best_possible || abort == TRUE)
            done = TRUE;

          /* Increment count of path-ends */
          npaths++;
          np++;

          /* If been processing a long time, then show progress */
          if (npaths % SHOWPATHS == 0)
            {
              printf("%12d paths searched. Best score = %d,\n",
                     npaths,best_score);
              printf(" (best possible = %d)\n",best_possible);
            }

          /* If have exceeded MAXPATHS, then possibly in an infinite
             loop, so abort here */
          if (np > MAXPATHS)
            {
              /* If this is the first time round, then show warning
                 message */
              if (nabort == 0)
                {
                  printf("*** Warning. Possibly in an infinite loop. ");
                  printf("Aborting graph-match\n");
                }

              /* Set flag so that we abort the current loop */
              done = TRUE;

              /* Increment count of aborts */
              nabort++;

              /* If have exceeded maximum number of aborts, then
                 abandon graph-matching altogether */
              if (nabort >= MAXABORT)
                {
                  printf("*** Warning. Have now aborted on %d such",
                         nabort);
                  printf(" occasions\n");
                  printf("***          Am giving up altogether\n");

                  /* Set abort flag */
                  abort = TRUE;
                }
            }
        }

      /* Get next entry on the tree to be processed */
      got_next = FALSE;

      /* Loop until have determine which is the next tree entry to
         process */
      while (got_next == FALSE && done == FALSE)
        {
          /* If can go down a level, get the tree entry to go to */
          if (go_down == TRUE)
            {
              new_tree_ptr = tree_ptr->next_level_tree_ptr;
              got_next = TRUE;
            }

          /* Otherwise, move across at the same level */
          else
            {
              /* Check whether can move across */
              new_tree_ptr = tree_ptr->next_tree_ptr;

              /* If can go across at current level, then have next
                 entry to go to */
              if (new_tree_ptr != NULL)
                got_next = TRUE;

              /* Otherwise, have nothing to move to at this level, so
                 need to go up a level */
              else
                {
                  /* Get pointer to higher level */
                  new_tree_ptr = tree_ptr->prev_level_tree_ptr;
                }

              /* If already at top level, then we are done */
              if (new_tree_ptr == NULL)
                done = TRUE;
            }

          /* Store the next tree pointer to be processed */
          tree_ptr = new_tree_ptr;
        }

      /* Prepare next entry */
      if (done == FALSE)
        {
          /* Re-initialise pointer to previous tree-entry at
             this level */
          prev_tree_ptr = NULL;

          /* Re-initialise all the flags on the atom and bond
             equivalences */
          reset_path(new_tree_ptr,first_equiv_atom_ptr,
                     first_equiv_bond_ptr,FALSE);
        }
    }

  /* Return the pointers to the start and end of the list of
     path-ends */
  *fst_end_tree_ptr = first_end_tree_ptr;
  *lst_end_tree_ptr = last_end_tree_ptr;
  *fst_all_tree_ptr = first_all_tree_ptr;
  *lst_all_tree_ptr = last_all_tree_ptr;

  /* Return the number of path-ends */
  *npth = npaths;
  if (VERBOSE == TRUE)
    printf("BEST SCORE: %d\n",best_score);

  /* Return best score and pointer to corresponding path-end */
  *bst_score = best_score;
  *bst_tree_ptr = best_tree_ptr;

  /* Return whether to abort whole graph-matching process */
  return(abort);
}
/***********************************************************************

free_paths  -  Free memory used by all paths in the current tree

***********************************************************************/

void free_paths(struct tree *first_tree_ptr)
{
  struct tree *tree_ptr, *next_tree_ptr;

  /* Start at the first tree pointers */
  tree_ptr = first_tree_ptr;

  /* Loop through all the tree entries */
  while (tree_ptr != NULL)
    {
      /* Save the pointer to the next entry */
      next_tree_ptr = tree_ptr->all_next_tree_ptr;

      /* Free the memory for the current tree entry */
      free(tree_ptr);

      /* Go to the next entry */
      tree_ptr = next_tree_ptr;
    }
}
/***********************************************************************

search_paths  -  Find longest (and best) path through the network of
                 potentially equivalenced bonds

***********************************************************************/

int search_paths(struct equiv_atom *first_equiv_atom_ptr,
                 struct equiv_bond *first_equiv_bond_ptr,
                 struct tree **bst_tree_ptr,int best_possible,
                 struct tree **fst_all_tree_ptr,
                 struct tree **lst_all_tree_ptr)
{
  int best_match, best_score, npaths;
  int abort, done;

  struct equiv_bond *equiv_bond_ptr;
  struct equiv_bond *start_equiv_bond_ptr;
  struct tree *best_tree_ptr, *first_tree_ptr;
  struct tree *first_end_tree_ptr, *last_end_tree_ptr;
  struct tree *first_all_tree_ptr, *last_all_tree_ptr;

  /* Initialise variables */
  abort = FALSE;
  best_score = 0;
  done = FALSE;
  npaths = 0;

  /* Initialise pointers */
  best_tree_ptr = NULL;
  *bst_tree_ptr = NULL;
  first_all_tree_ptr = *fst_all_tree_ptr;
  last_all_tree_ptr = *lst_all_tree_ptr;

  /* Loop until all possible start-points have been exhausted, or
     best possible match achieved */
  while (done == FALSE)
    {
      /* Initialise best match score */
      best_match = 0;
      start_equiv_bond_ptr = NULL;

      /* Start at the first of the potentially equivalenced bonds */
      equiv_bond_ptr = first_equiv_bond_ptr;

      /* Loop through the equivalenced bonds to find best path-start */
      while (equiv_bond_ptr != NULL)
        {
          /* If this bond-pair can be a path-start, then check if the
             best so far */
          if (equiv_bond_ptr->possible_start == TRUE)
            {
              /*
              if (VERBOSE == TRUE)
                printf("Considering %4d. [%s]->[%s] == [%s]->[%s] %d\n",
                       equiv_bond_ptr->link_number,
                       equiv_bond_ptr->equiv_atom_ptr->atom_ptr->atom_name,
                       equiv_bond_ptr->other_equiv_atom_ptr->atom_ptr->atom_name,
                       equiv_bond_ptr->equiv_atom_ptr->dic_data_ptr->atom_name,
                       equiv_bond_ptr->other_equiv_atom_ptr->dic_data_ptr->atom_name,
                       equiv_bond_ptr->match_score);
              */

              /* Check if score is higher than best so far */
              if (equiv_bond_ptr->match_score > best_match)
                {
                  /* Save the pointer */
                  start_equiv_bond_ptr = equiv_bond_ptr;

                  /* Save the best score */
                  best_match = equiv_bond_ptr->match_score;
                }
            }

          /* Go to next equivalenced bond */
          equiv_bond_ptr = equiv_bond_ptr->next_equiv_bond_ptr;
        }

      /* If no more left, then we are done */
      if (start_equiv_bond_ptr == NULL)
        done = TRUE;

      /* Otherwise, generate tree of paths starting at the current 
         bond equivalence */
      else
        {
          /* Get the starting bond equivalence */
          equiv_bond_ptr = start_equiv_bond_ptr;

          /* Show starting bond equivalence */
          if (VERBOSE == TRUE)
            printf("       Start-bond %4d. [%s]->[%s] == [%s]->[%s] %d\n",
                   equiv_bond_ptr->link_number,
                   equiv_bond_ptr->equiv_atom_ptr->atom_ptr->atom_name,
                   equiv_bond_ptr->other_equiv_atom_ptr->atom_ptr->atom_name,
                   equiv_bond_ptr->equiv_atom_ptr->dic_data_ptr->atom_name,
                   equiv_bond_ptr->other_equiv_atom_ptr->dic_data_ptr->atom_name,
                   equiv_bond_ptr->match_score);

          /* Initialise pointers */
          first_end_tree_ptr = NULL;
          last_end_tree_ptr = NULL;

          /* Initialise all the equivalent atom and bond records */
          initialise_equiv_atoms(first_equiv_atom_ptr);
          initialise_equiv_bonds(first_equiv_bond_ptr);
          npaths = 0;

          /* Create an entry on the tree for this equivalent bond */
          first_tree_ptr = create_tree_entry(equiv_bond_ptr,NULL,NULL,
                                             &first_all_tree_ptr,
                                             &last_all_tree_ptr);

          /* Process the entries in the tree until all possible paths
             from this start-point have been explored */
          abort = process_tree(first_tree_ptr,first_equiv_atom_ptr,
                               first_equiv_bond_ptr,&first_end_tree_ptr,
                               &last_end_tree_ptr,&npaths,
                               &first_all_tree_ptr,
                               &last_all_tree_ptr,best_possible,
                               &best_score,&best_tree_ptr);

          /* If this equals the best possible, then end here */
          if (best_score >= best_possible || abort == TRUE)
            done = TRUE;
        }
    }

  /* Return pointer to last tree entry created */
  *fst_all_tree_ptr = first_all_tree_ptr;
  *lst_all_tree_ptr = last_all_tree_ptr;

  /* Return pointer to best path found */
  *bst_tree_ptr = best_tree_ptr;

  /* Return the number of path-ends */
  return(best_score);
}
/***********************************************************************

process_best_path  -  Use the best atom-atom match between the residue
                      and the corresponding Het Group Dictionary entry,
                      linking the dictionary atoms to the corresponding
                      residue atoms

***********************************************************************/

int process_best_path(struct residue *residue_ptr,
                      struct equiv_atom *first_equiv_atom_ptr,
                      struct equiv_bond *first_equiv_bond_ptr,
                      int *nmch,int *nmap,int *nmissd,int *nmissg,int *nbmatch)
{
  int iatom, natoms;
  int nbonds_matched, nmapped, nmatched, nmissed, nmissing;
  int ok;

  struct atom *atom_ptr;
  struct dic_data *dic_data_ptr;
  struct dictionary *dictionary_ptr;
  struct equiv_atom *equiv_atom_ptr;

  /* Get the dictionary pointer corresponding to this residue */
  dictionary_ptr = residue_ptr->dictionary_ptr;

  /* Initialise count of matching atoms */
  ok = TRUE;
  nbonds_matched = 0;
  nmatched = 0;
  nmapped = 0;
  nmissed = 0;
  nmissing = 0;

  /* Calculate score for the equivalenced bonds */
  nbonds_matched = get_bond_score(first_equiv_bond_ptr);

  /* Get pointer to the start of linked list of equiv_atom coords */
  equiv_atom_ptr = first_equiv_atom_ptr;

  /* Loop through all the equiv_atoms to pick up those that are valid */
  while (equiv_atom_ptr != NULL)
    {
      /* If this atom-equivalence is on the current path, then include */
      if (equiv_atom_ptr->checked == TRUE)
        {
          /* Increment count of matches */
          nmatched++;

          /* Get the mapping represented by this equivalence */
          atom_ptr = equiv_atom_ptr->atom_ptr;
          dic_data_ptr = equiv_atom_ptr->dic_data_ptr;

          /* Mark these as accounted for */
          atom_ptr->checked = TRUE;
          dic_data_ptr->checked = TRUE;

          /* Establish pointers between the dictionary atom and the
             residue atom */
          dic_data_ptr->atom_ptr = atom_ptr;
          atom_ptr->dic_data_ptr = dic_data_ptr;

          /* If atom names not identical, then show mapping applied */
          if (strncmp(atom_ptr->atom_name,dic_data_ptr->atom_name,4))
            {
              /* Increment count of mismatched mappings */
              nmapped++;
              printf("           * Res atom [%s] maps to dict atom [%s]\n",
                     equiv_atom_ptr->atom_ptr->atom_name,
                     equiv_atom_ptr->dic_data_ptr->atom_name);
            }

          /* Copy atom-name from residue into the output name to be
             written to the hbplus.rc file */
          strncpy(dic_data_ptr->output_atom_name,atom_ptr->atom_name,4);
          dic_data_ptr->output_atom_name[4] = '\0';

          /* Write mapping out to hbadd.map file */
          fprintf(fil_map,"[%s]   [%s]  ->  [%s]",
                  residue_ptr->res_name,atom_ptr->atom_name,
                  dic_data_ptr->atom_name);

          /* If atom names match, show asterisk */
          if (!strcmp(atom_ptr->atom_name,dic_data_ptr->atom_name))
            fprintf(fil_map,"  *\n");

          /* Otherwise, show that mapping occurred */
          else
            fprintf(fil_map,"  Mapped\n");
        }

      /* Get pointer to next atom equivalence */
      equiv_atom_ptr = equiv_atom_ptr->next_equiv_atom_ptr;
    }   

  /* Get pointer to the first of this residue's atoms */
  atom_ptr = residue_ptr->first_atom_ptr;
  natoms = residue_ptr->natoms;

  /* Loop through the input atoms to write out those that were not
     matched against dictionary atoms */
  for (iatom = 0; iatom < natoms && atom_ptr != NULL; iatom++)
    {
      /* If atom not matched and nopt hydrogen, then write out as
         missing from dictionary */
      if (atom_ptr->checked == FALSE && atom_ptr->hydrogen == FALSE)
        {
          /* Write out to the map file */
          fprintf(fil_map,"[%s]   [%s]  ->  xxxxxx  Missed\n",
                  residue_ptr->res_name,atom_ptr->atom_name);

          /* Increment count */
          nmissed++;
        }

      /* Go to the next atom */
      atom_ptr = atom_ptr->next_atom_ptr;
    }

  /* Get pointer to the first of this dictionary residue's entries */
  dic_data_ptr = dictionary_ptr->first_dic_data_ptr;
  natoms = dictionary_ptr->natoms;

  /* Loop through the input dictionary atoms to write out those missing
     from the residue in the PDB file */
  for (iatom = 0; iatom < natoms && dic_data_ptr != NULL; iatom++)
    {
      /* If entry not matched and not hydrogen, then write out as
         missing from residue */
      if (dic_data_ptr->checked == FALSE && dic_data_ptr->hydrogen == FALSE)
        {
          /* Write out to the map file */
          fprintf(fil_map,"[%s]   xxxxxx  ->  [%s]  Missing\n",
                  residue_ptr->res_name,dic_data_ptr->atom_name);

          /* Increment count */
          nmissing++;
        }

      /* Go to the next dic_data */
      dic_data_ptr = dic_data_ptr->next_dic_data_ptr;
    }

  /* Show number of atoms mapped onto atoms with a different name */
  if (nmapped > 0)
    {
      printf("         Number of atoms mapped to atoms with a different ");
      printf("name:  %5d\n",nmapped);
      printf("         Number of identical matches:                     ");
      printf("       %5d\n",nmatched - nmapped);
    }

  /* If too few atoms matched, then take this to be an unsuccessful
     match, and abort */
  if ((residue_ptr->non_hydrogen < dictionary_ptr->non_hydrogen &&
       nmatched < residue_ptr->non_hydrogen - 1) ||
      (residue_ptr->non_hydrogen >= dictionary_ptr->non_hydrogen &&
       nmatched < dictionary_ptr->non_hydrogen - 1))
    {
      printf("           *** Warning. Can only map %d of %d atoms to ",
             nmatched,residue_ptr->non_hydrogen);
      printf("dictionary entry\n");
      ok = FALSE;
    }

  /* If have atoms present in the dictionary but not in the residue,
     list them */
  if (ok == TRUE &&
      residue_ptr->non_hydrogen < dictionary_ptr->non_hydrogen)
    {
      /* Get pointer to dictionary's first atom entry */
      dic_data_ptr = dictionary_ptr->first_dic_data_ptr;
      natoms = dictionary_ptr->natoms;

      /* Write heading to hbadd.map file */
      fprintf(fil_map,"\nAtoms in dictionary but not in PDB entry:\n");

      /* Loop through all the dic_data atoms to pick up the missed
         ones */
      for (iatom = 0; iatom < natoms && dic_data_ptr != NULL; iatom++)
        {
          /* If entry has no pointer to a residue atom, then is one
             that is missing from the residue */
          if (dic_data_ptr->atom_ptr == NULL &&
              dic_data_ptr->hydrogen == FALSE)
            {
              printf("            * Warning. Atom [%s] missing ",
                     dic_data_ptr->atom_name);
              printf("from residue %s\n",
                     dictionary_ptr->res_name);

              /* Write out to hbadd.map file */
              fprintf(fil_map,"   [%s] [%s]\n",residue_ptr->res_name,
                      dic_data_ptr->atom_name);
            }

          /* Get pointer to next dic_data in linked-list */
          dic_data_ptr = dic_data_ptr->next_dic_data_ptr;
        }
    }

  /* Otherwise, is have atoms present in the residue but not in the
     dictionary, list these */
  else if (ok == TRUE &&
           residue_ptr->non_hydrogen > dictionary_ptr->non_hydrogen)
    {
      /* Write heading to hbadd.map file */
      fprintf(fil_map,"\nAtoms in PDB entry not mapped to dictionary:\n");

      /* Get pointer to residue's first atom */
      atom_ptr = residue_ptr->first_atom_ptr;
      natoms = residue_ptr->natoms;

      /* Loop through all the atoms to pick up the missed ones */
      for (iatom = 0; iatom < natoms && atom_ptr != NULL; iatom++)
        {
          /* If atom has no pointer to a dictionary entry, then is one
             that is extra to the residue */
          if (atom_ptr->dic_data_ptr == NULL &&
              atom_ptr->hydrogen == FALSE)
            {
              printf("         *** Warning. Extra atom, not present in");
              printf(" dictionary [%s]\n",atom_ptr->atom_name);

              /* Write to hbadd.map file */
              fprintf(fil_map,"     [%s] [%s]\n",residue_ptr->res_name,
                      dic_data_ptr->atom_name);
            }
      
          /* Get pointer to next atom in linked-list */
          atom_ptr = atom_ptr->next_atom_ptr;
        }               
    }

  /* Write out missed and missing atoms */
  if (nmissed + nmissing > 0)
    {
      printf("                Number of missed atoms: %3d",nmissed);
      printf("   Number of missing atoms: %3d\n",nmissing);

      /* Determine whether match to dictionary entry is close enough */
/* v.5.0.1--> */
/*      if (nmissed + nmissing > 4) */
      if (nmissed + nmissing > MAX_MISSING)
/* <--v.5.0.1 */
        ok = FALSE;
    }

  /* If satisfied that the residue matches the dictionary entry
     sufficiently well, then link the dictionary entry to this residue */
  if (ok == TRUE)
    dictionary_ptr->residue_ptr = residue_ptr;

  /* Return the match-counts */
  *nmch = nmatched;
  *nbmatch = nbonds_matched / 2;
  *nmap = nmapped;
  *nmissd = nmissed;
  *nmissg = nmissing;

  /* Return whether mapping was successful */
  return(ok);
}
/***********************************************************************

write_sumfile  -  Write the hit/miss statistics to the summary file,

***********************************************************************/

void write_sumfile(struct residue *residue_ptr,int nmatched,int nmapped,
                   int nmissed,int nmissing,int nbonds_matched,
                   char *sum_name,char *pdb_name)
{
  struct dictionary *dictionary_ptr;

  /* Get the dictionary pointer corresponding to this residue */
  dictionary_ptr = residue_ptr->dictionary_ptr;

  /* If file not open, then open it */
  if (fil_sum == NULL)
    open_sumfile(sum_name,pdb_name);

  /* Summarise match success/failure in the hbadd.sum file */
  fprintf(fil_sum,"%s %c%s ",residue_ptr->res_name,residue_ptr->chain,
          residue_ptr->res_num);

  /* Print numbers of atom in residue and number matched */
  fprintf(fil_sum," %3d  %3d  %3d  %3d    %3d    %3d      %3d ",
          residue_ptr->non_hydrogen,dictionary_ptr->non_hydrogen,
          nmatched,nmatched - nmapped,nmapped,nmissed,nmissing);

  /* Show asterisks if residue unmapped */
  if (nmissing > 0 && nmissed > 0)
    fprintf(fil_sum,"   XXX");
  else if (nmissing > 0)
    fprintf(fil_sum,"   ---");
  else if (nmissed > 0)
    fprintf(fil_sum,"   ***");
  else
        fprintf(fil_sum,"      ");

  /* Show numbers of matched bonds */
  fprintf(fil_sum," %3d  %3d  %3d",residue_ptr->nconnect,
          dictionary_ptr->nconnect,nbonds_matched);

  /* Close off line */
  fprintf(fil_sum,"\n");
}
/***********************************************************************

graph_match  -  Perform a simple graph-match between the atoms of the
                residue and those of the entry from the Het Group
                Dictionary to get the best atom-atom correspondences

***********************************************************************/

int graph_match(struct residue *residue_ptr,
                struct connect *first_connect_ptr,char *sum_name,
/* v.5.0--> */
/*                 char *pdb_name) */
                char *pdb_name,int match_names)
/* <--v.5.0 */
{
  char asterisk;
/* v.5.4.5--> */
  char element[3];
/* <--v.5.4.5 */

  int best_possible, max_equiv;
  int nconnect, nequiv_atoms, nequiv_bonds;
  int nbonds_matched, nmapped, nmatched, nmissed, nmissing;
/* v.5.0--> */
  int nbonds;
/* <--v.5.0 */
 int ok, perfect, perfect_possible, trace;
/* v.5.4.5--> */
 int is_metal;
/* <--v.5.4.5 */

  struct atom *atom_ptr;
  struct dic_data *dic_data_ptr;
  struct dictionary *dictionary_ptr;
  struct equiv_atom *first_equiv_atom_ptr;
  struct equiv_bond *first_equiv_bond_ptr;
  struct tree *best_tree_ptr;
  struct tree *first_all_tree_ptr, *last_all_tree_ptr;

  /* Initialise variables */
  nbonds_matched = 0;
  nmatched = 0;
  nmissed = 0;
  nmissing = residue_ptr->non_hydrogen;
  nmapped = 0;
  ok = FALSE;
  trace = FALSE;

  /* Initialise pointers */
  best_tree_ptr = NULL;
  first_equiv_atom_ptr = NULL;
  first_equiv_bond_ptr = NULL;
  first_all_tree_ptr = NULL;
  last_all_tree_ptr = NULL;

  /* Get the dictionary pointer corresponding to this residue */
  dictionary_ptr = residue_ptr->dictionary_ptr;

  /* Show numbers of atoms in PDB file and dictionary */
  asterisk = ' ';
  if (residue_ptr->non_hydrogen != dictionary_ptr->non_hydrogen)
    asterisk = '*';
  printf("    Residue %s. Atoms: %6d in PDB entry, %6d in dictionary  %c\n",
         residue_ptr->res_name,residue_ptr->non_hydrogen,
         dictionary_ptr->non_hydrogen,asterisk);

  /* Calculate all the covalent connectivities for the atoms of the
     residue */
  nconnect = calc_covalent_connectivities(residue_ptr,first_connect_ptr);

  /* Show numbers of bonds in PDB file and dictionary */
  asterisk = ' ';
  if (nconnect != dictionary_ptr->nconnect)
    asterisk = '*';
  printf("                 Bonds: %6d in PDB entry, %6d in dictionary  %c\n",
         nconnect,dictionary_ptr->nconnect,asterisk);

  /* Store count of this residue's bonds */
  residue_ptr->nconnect = nconnect;

  /* Create the equiv_atoms of the graph by pairing of like-minded atoms
     between the residue and the entry from the Het Group Dictionary */
  create_equiv_atoms(residue_ptr,&first_equiv_atom_ptr,&nequiv_atoms,
/* v.5.0--> */
/*                     &max_equiv); */
                     &max_equiv,match_names);
/* <--v.5.0 */
  printf("                    Number of equiv_atoms created = %d\n",
         nequiv_atoms);

  /* Calculate links between equiv_atoms */
  nequiv_bonds = calc_equiv_bonds(first_equiv_atom_ptr,
                                  &first_equiv_bond_ptr,max_equiv);
  printf("                    Number of equiv_bonds created = %d\n",
         nequiv_bonds);

  /* If there is only a single atom, then map directly */
  if (residue_ptr->non_hydrogen == dictionary_ptr->non_hydrogen &&
      residue_ptr->non_hydrogen == 1)
    {
      /* Increment count of matches */
      nmatched = 1;
      nmissed = 0;
      nmissing = 0;
      nmapped = 0;

      /* Show match */
      printf("                    Single-atom residue %s %s %c",
             residue_ptr->res_name,residue_ptr->res_num,residue_ptr->chain);
      printf("  matched vs %s\n",dictionary_ptr->res_name);


      /* Get the residue and dictionary atoms */
      atom_ptr = residue_ptr->first_atom_ptr;
      dic_data_ptr = dictionary_ptr->first_dic_data_ptr;

      /* If atom names not identical, then show mapping applied */
      if (strncmp(atom_ptr->atom_name,dic_data_ptr->atom_name,4))
        {
          /* Show the residue name */
          printf("                    *** Residue:-  %s %s %c ",
                 residue_ptr->res_name,
                 residue_ptr->res_num,
                 residue_ptr->chain);
          printf("Name mismatch:\n");

          /* Increment count of mismatched mappings */
          nmapped = 1;
          printf("           * Res atom ");
          printf("[%s] maps to dict atom [%s]\n",
                 atom_ptr->atom_name,dic_data_ptr->atom_name);
        }

      /* Copy atom-name from residue into the output name to be
         written to the hbplus.rc file */
      strncpy(dic_data_ptr->output_atom_name,atom_ptr->atom_name,4);
      dic_data_ptr->output_atom_name[4] = '\0';

/* v.5.4.5--> */
      /* See if this atom is a metal */
      is_metal = check_for_metal(atom_ptr->atom_name,element);

      /* Save the link to the residue */
      if (is_metal == TRUE)
        dictionary_ptr->residue_ptr = residue_ptr;
/* <--v.5.4.5 */

      /* Flag that match was alright */
      ok = TRUE;
    }

  /* Otherwise, perform a graph-matching between the connectivities
     of the atoms in the PDB file and those in the Het Group Dictionary */
  else
    {
      /* Initialise indicators of which atoms have been processed */
      initialise_atoms(residue_ptr);
      initialise_dic_data(dictionary_ptr);

      /* Pre-calculate the best possible score that can be achieved */
      best_possible
        = get_best_possible_score(residue_ptr,&perfect_possible);

/* v.5.0--> */
      /* Add in best possible bond score */
      nbonds = residue_ptr->nconnect;
      if (dictionary_ptr->nbonds < residue_ptr->nconnect)
        nbonds = dictionary_ptr->nbonds;
      best_possible = best_possible + 2 * nbonds;
/* <--v.5.0 */

      /* Check whether we already have a perfect match between the 
         ligand atom names in the PDB file and those in the Het Group
         Dictionary definition */
      if (perfect_possible == TRUE &&
          residue_ptr->non_hydrogen == dictionary_ptr->non_hydrogen)
        perfect = check_perfect_match(first_equiv_bond_ptr,
                                      first_equiv_atom_ptr,&nbonds_matched);
      else
        perfect = FALSE;

      /* If have a perfect match, then use that */
      if (perfect == TRUE)
        {
          /* Show message */
          printf("          Have a perfect match with the dictionary");
          printf(" for %s ...\n",residue_ptr->res_name);

          /* Retrieve atom-equivalences corresponding to the perfect match */
          ok = process_best_path(residue_ptr,first_equiv_atom_ptr,
                                 first_equiv_bond_ptr,&nmatched,
                                 &nmapped,&nmissed,&nmissing,&nbonds_matched);
        }

      /* Otherwise, need to perform a graph-match */
      else
        {
          /* Show message */
          printf("       Performing graph-match for %s ...\n",
                 residue_ptr->res_name);

          /* Calculate all possible paths linking the equivalent bonds */
          search_paths(first_equiv_atom_ptr,first_equiv_bond_ptr,
                       &best_tree_ptr,best_possible,
                       &first_all_tree_ptr,&last_all_tree_ptr);

          /* If have a best-path, then retrieve the atom-equivalences it
             represents */
          if (best_tree_ptr != NULL)
            {
              /* Initialise which atom-equivalences are included in
                 the best path */
              reset_path(best_tree_ptr,first_equiv_atom_ptr,
                         first_equiv_bond_ptr,trace);

              /* Process the best set of equivalences identified */
              ok = process_best_path(residue_ptr,first_equiv_atom_ptr,
                                     first_equiv_bond_ptr,&nmatched,
                                     &nmapped,&nmissed,&nmissing,
                                     &nbonds_matched);
            }

          /* If no good match obtained, then set counts */
          else
            {
              /* Calculate counts */
              nmissed = dictionary_ptr->non_hydrogen;
              nmissing = residue_ptr->non_hydrogen;

              /* Show warning */
              printf("*** Warning. Cannot match Het Group %s to ",
                     residue_ptr->res_name);
              printf("dictionary entry\n");
            }

          /* Free up the memory used by the tree entries */
          free_paths(first_all_tree_ptr);
        }
    }

  /* Write the hit/miss statistics to the summary file */
  write_sumfile(residue_ptr,nmatched,nmapped,nmissed,nmissing,
                nbonds_matched,sum_name,pdb_name);

  /* Store the match stats for this residue */
  residue_ptr->nmatched = nmatched;
  residue_ptr->nmapped = nmapped;
  residue_ptr->nmissed = nmissed;
  residue_ptr->nmissing = nmissing;
  residue_ptr->nbonds_matched = nbonds_matched;

  /* Return whether graph-match was successful */
  return(ok);
}
/***********************************************************************

validate_pdb_hets  -  Validate atom names and connectivities of the HET
                      groups in the PDB file according to the data in the
                      HET Group Dictionary, listing any groups not found
                      in the dictionary

***********************************************************************/

void validate_pdb_hets(struct residue *first_residue_ptr,
                       struct dictionary **first_dictionary_ptr,
                       struct dictionary **last_dictionary_ptr,
                       int *ndic_entries,
                       struct dic_data **last_dic_data_ptr,
                       struct connect *first_connect_ptr,
/* v.5.0--> */
/*                       char *sum_name,char *pdb_name,char *dic_res_name) */
                       char *sum_name,char *pdb_name,char *dic_res_name,
                       int match_names)
/* <--v.5.0 */
{
  int ok, printed_warning;

  struct residue *residue_ptr;
  struct dictionary *dictionary_ptr;

  /* Initialise variables */
  printed_warning = FALSE;
  ok = FALSE;

  /* Initialise residue-pointer to start of linked list */
  residue_ptr = first_residue_ptr;

  /* Loop through all the HET group residues */
  while (residue_ptr != NULL)
    {
      /* If have a matching entry in the PDB, then validate all the atom
         names and connectivities */
      if (residue_ptr->dictionary_ptr != NULL && residue_ptr->nmatch == 1)
        {
          /* Get the dictionary entry corresponding to this residue */
          dictionary_ptr = residue_ptr->dictionary_ptr;

          /* If this dictionary entry has not yet been matched against
             a residue, then perform the matching */
          if (dictionary_ptr->residue_ptr == NULL)
            {
              /* Set flag */
              ok = FALSE;

              /* Apply a simple graph-matching technique to match
                 atoms between the residue and the entry in the Het
                 Group Dictionary */
              if (dictionary_ptr->natoms > 0)
                ok = graph_match(residue_ptr,first_connect_ptr,sum_name,
/* v.5.0--> */
/*                                  pdb_name); */
                                 pdb_name,match_names);
/* <--v.5.0 */
              else
                printf("    Residue %s. No atoms defined in dictionary\n",
                       residue_ptr->res_name);

              /* If graph-match failed, then remove residue-to-dictionary
                 link */
              if (ok == FALSE && dic_res_name[0] == '\0')
                {
                  residue_ptr->dictionary_ptr = NULL;
                  residue_ptr->nmatch = 0;
                }
            }
        }

      /* Get pointer to the next residue in the list */
      residue_ptr = residue_ptr->next_residue_ptr;
    }

  /* Re-initialise residue-pointer to start of linked list */
  residue_ptr = first_residue_ptr;
  if (dic_res_name[0] != '\0')
    residue_ptr = NULL;

  /* Loop through all the HET group residues to pick up any for which
     no dictionary entry has been written out */
  while (residue_ptr != NULL)
    {
      /* If this het group has no pointer to a dictionary entry, then
         must be missing from the dictionary */
      if (residue_ptr->dictionary_ptr == NULL)
        {
          /* If we haven't come across this residue type before, print
             warning and flag all others of the same type as done */
          if (residue_ptr->nmatch == 0)
            {
              /* Print warning message about this HET group */
              print_missing_hets(residue_ptr);
              printed_warning = TRUE;

              /* Do the best that can be done with the information
                 available - ie create a dictionary entry from the
                 atom-names and their connectivities, as in the PDB file */
              generate_fake_entry(residue_ptr,first_residue_ptr,
                                  &(*first_dictionary_ptr),
                                  &(*last_dictionary_ptr),
                                  &(*ndic_entries),
                                  &(*last_dic_data_ptr),first_connect_ptr);
            }
        }

      /* Get pointer to the next residue in the list */
      residue_ptr = residue_ptr->next_residue_ptr;
    }

  /* If any warning messages printed, then throw a blank line */
  if (printed_warning == TRUE)
    printf("\n");
}
/***********************************************************************

convert_quotes - Convert any double-quotes in the atom name to @'s

***********************************************************************/

void convert_quotes(char atom_name[5])
{
  int ipos;

  /* Loop through all the characters, making the conversion if necessary */
  for (ipos = 0; ipos < 4; ipos++)
    if (atom_name[ipos] == '"')
      atom_name[ipos] = '@';
}
/***********************************************************************

writename_fn - Function to write atom and residue names to hbplus.rc file

***********************************************************************/

void writename_fn(char res_name[4],int natoms,
                  struct dic_data *first_dic_data_ptr)
{
   int i, iatom;
   struct dic_data *dic_data_ptr;

   /* Initialize variables */
   i = iatom = 0;

   /* Write out the residue name */
   fprintf(fil_out, "-u \'%s\'\n",res_name);

   /* Initialise pointer to the first of the atom names */
   dic_data_ptr = first_dic_data_ptr;

   /* Loop through all the atoms to write out to file */
   while (iatom < natoms && dic_data_ptr != NULL)
     {
       /* Convert any double-quotes in the atom name to @'s */
       convert_quotes(dic_data_ptr->output_atom_name);

       /* Write out only if a non-hydrogen atom */
/* v.5.1--> */
//       if (dic_data_ptr->atom_name[1] != 'H')
       if (dic_data_ptr->hydrogen == FALSE)
/* <--v.5.1 */
         {
           /* Start new line every 10 atoms */
           if (i >= 10)
             {
               fprintf(fil_out,"\"\n");
               i = 0;
             }
           if (i == 0)
             fprintf(fil_out,"-M \'%s\' \"",res_name);
           
           /* Write out the current atom name */
           fprintf(fil_out, "%s",dic_data_ptr->output_atom_name);
           i++;
         }

       /* Set pointer to the next atom and increment counter */
       dic_data_ptr = dic_data_ptr->next_dic_data_ptr;
       iatom++;
     }
   /* Write out the end part of the current line */
   fprintf(fil_out,"\"\n");
}
/***********************************************************************

writecon_fn - Function to write the connectivities to the hbplus.rc file

***********************************************************************/

void writecon_fn(char res_name[4],int natoms,
                 struct dic_data *first_dic_data_ptr)
{
  int i, iatom, jatom;
  int flag;
  int flagb;
  struct dic_data *dic_data_ptr, *other_dic_data_ptr;

  /* Initialize variables */
  i = iatom = 0;
  flag  = 0;
  flagb = 0;

  /* Initialise pointer to the first of the atom names */
  dic_data_ptr = first_dic_data_ptr;
  
  /* Loop through all the atoms to write out to file */
  while (iatom < natoms && dic_data_ptr != NULL)
    {
      /* Loop through the connected atoms */
      for (jatom = 0; jatom < dic_data_ptr->nbound; jatom++)
        {
          /* Write out only if both atoms are non-hydrogens */
/* v.5.1--> */
//          if (dic_data_ptr->atom_name[1] != 'H' &&
//              dic_data_ptr->boundto[jatom][1] != 'H')
          if (dic_data_ptr->hydrogen == FALSE &&
             is_hydrogen_atom(dic_data_ptr->boundto[jatom]) == FALSE)
/* <--v.5.1 */
            {
              flagb = 1;
              /* Start new line every 5 pairs */
          
              /* Flag is to prevent writing "-T" twice on same line */
              if (i >= 5)
                {
                  fprintf(fil_out,"\"\n");
                  i    = 0;
                  flag = 0;
                }
              if (i == 0 && flag == 0)
                {
                  fprintf(fil_out,"-T \'%s\' \"",res_name);
                  flag = 1;
                }
              
              /* Write connectivity information */
              if (dic_data_ptr->boundto_num[jatom] >= dic_data_ptr->number)
                {
                  other_dic_data_ptr = dic_data_ptr->connect_atom_ptr[jatom];
                  fprintf(fil_out,"%s%s:",dic_data_ptr->output_atom_name,
                          other_dic_data_ptr->output_atom_name);
                  i++;
                }
            }
        }

      /* Set pointer to next atom */
      dic_data_ptr = dic_data_ptr->next_dic_data_ptr;
      iatom++;
    }
  /* Close line */
  if (flagb == 1)
    fprintf(fil_out,"\"\n");
}
/***********************************************************************

writehbdon_fn - Function to write Hbond donor information to hbplus.rc 

***********************************************************************/

void writehbdon_fn(char res_name[4],int natoms,
                   struct dic_data *first_dic_data_ptr,int fake_entry,
                   int have_hydrogens)
{
  int iatom, jatom;
  int nhyd;                              /* Count no. of H atoms bonded */   
  int ncoval;
  int metal;

  struct dic_data *dic_data_ptr;

  /* Initialize variables */
  iatom = 0;
  metal = FALSE;
   
  /* Initialise pointer to the first of the atom names */
  dic_data_ptr = first_dic_data_ptr;

  /* Loop through all the atoms to write out to file */
  while (iatom < natoms && dic_data_ptr != NULL)
    {
      /* If residue consists of a single atom, check whether it is a metal */
/* v.5.3.7 --> */
//      metal = check_for_metal(dic_data_ptr->atom_name,element);
      metal = check_for_metal(dic_data_ptr->atom_name,
                              dic_data_ptr->element);
/* <-- v.5.3.7 */

      /* Check if this atom is O, S or N - these potentially donate
         the same number of hydrogens as H atoms bonded to them */
      if (!strcmp(dic_data_ptr->element," O") ||
          !strcmp(dic_data_ptr->element," S") ||
          !strcmp(dic_data_ptr->element," N") || metal == TRUE)
        {
          /* Initialise count of hydrogens attached to current atom */
          nhyd = 0;
          ncoval = 0;

          /* Count number of H atoms bonded to */
          for (jatom = 0; jatom < dic_data_ptr->nbound; jatom++)
            {
/* v.5.1--> */
//              if (dic_data_ptr->boundto[jatom][0] == 'H' ||
//                  dic_data_ptr->boundto[jatom][1] == 'H')
              if (is_hydrogen_atom(dic_data_ptr->boundto[jatom]) == TRUE)
/* <--v.5.1 */
                nhyd++;
              else
                ncoval++;
            }

          /* If this is a metal, then write out up to 6 potential donors
             so that HBPLUS picks up polar atoms coordinating this metal */
          if (metal == TRUE)
            {
/* v.4.4.3--> */
              if (ncoval == 0)
                {
/* <--v.4.4.3 */
                  nhyd = 8 - ncoval;
                  fprintf(fil_out,"-E \'%s\' \"%s\" %d\n",
                          res_name,dic_data_ptr->output_atom_name,nhyd);
/* v.4.4.3--> */
                }
/* <--v.4.4.3 */
            }

          /* Write out number of Hbond donors = number of hydrogens found */
          else if (nhyd > 0)
            fprintf(fil_out,"-E \'%s\' \"%s\" %d\n",
                    res_name,dic_data_ptr->output_atom_name,nhyd);

          /* Otherwise, if no hydrogens found, and this is a fake entry,
             then write out the maximum number of donors for this type
             of atom */
          else if (fake_entry == TRUE && have_hydrogens == FALSE)
            {
              /* Assign maximum number of hydrogens */
              if (!strcmp(dic_data_ptr->element," O"))
                nhyd = 1;
              else if (!strcmp(dic_data_ptr->element," S"))
                nhyd = 1;
              else if (!strcmp(dic_data_ptr->element," N"))
                nhyd = 3;

              /* Write out the maximum number of hydrogen donors */
              if (nhyd - ncoval > 0)
                fprintf(fil_out,"-E \'%s\' \"%s\" %d\n",
                        res_name,dic_data_ptr->output_atom_name,
                        nhyd - ncoval);
            }
        }
      /* Set pointer to next atom */
      dic_data_ptr = dic_data_ptr->next_dic_data_ptr;
      iatom++;
    }
}
/***********************************************************************

calc_angle  -  Calculate angles between the 3 sets of coordinates
               supplied

***********************************************************************/

float calc_angle(float coords1[3], float coords2[3], float coords3[3])
{
  int icoord;
  float angle, cosang, diff1, diff2, diff3, dist1, dist2, dist3;

  /* Initialise variables */
  angle = 0.0;
  dist1 = 0;
  dist2 = 0;
  dist3 = 0;

  /* Calculate distances between points 1->2, 2->3, 3->1 */
  for (icoord = 0; icoord < 3; icoord++)
    {
      diff1 = coords2[icoord] - coords1[icoord];
      diff2 = coords3[icoord] - coords2[icoord];
      diff3 = coords1[icoord] - coords3[icoord];
      dist1 = dist1 + diff1 * diff1;
      dist2 = dist2 + diff2 * diff2;
      dist3 = dist3 + diff3 * diff3;
    }

  /* Calculate angle */
  if (dist1 > 0.0 && dist2 > 0.0 && dist3 > 0.0)
    {
      dist1 = (float) sqrt((double)dist1);
      dist2 = (float) sqrt((double)dist2);
      dist3 = (float) sqrt((double)dist3);
      cosang = (dist3 * dist3 - dist1 * dist1 - dist2 * dist2)
        / ((float) -2.0 * dist1 * dist2);
      if (fabs(cosang) <= 1.0)
        angle = (float) (RADDEG * acos(cosang));
      else
        angle = (float) 0.0;
    }

  /* Return the angle */
  return(angle);
}
/***********************************************************************

writehbacc_fn - Function to write Hbond acceptor information to hbplus.rc 

***********************************************************************/

void writehbacc_fn(char res_name[4],int natoms,
                   struct dictionary *dictionary_ptr,
                   struct dic_data *first_dic_data_ptr)
{
  int iatom, istore, jatom, jstore, naccept, nangles;
  int nplanar, nstore, nsp3;
  int metal;

  float angle, coords_N[3], coords_other[MAXCONNECT][3];

  struct atom *atom_ptr;
  struct dic_data *dic_data_ptr, *other_dic_data_ptr;

  /* Initialize variables */
/*  printed_warning = FALSE; */

  /* Get the HET group from the PDB file that matched all the
     atom names and connectivities */
/*  residue_ptr = dictionary_ptr->residue_ptr; */

  /* If a valid pointer, then map all the atoms in the HET group to the
     corresponding coordinates in the PDB file */
/*  if (residue_ptr != NULL)
    check_atoms_present(residue_ptr,&printed_warning,&nmatch); */
  
  /* Re-initialise pointer to the first of the atom names */
  dic_data_ptr = first_dic_data_ptr;
  iatom = 0;  

  /* Loop through all the atoms to write out to file */
  while (iatom < natoms && dic_data_ptr != NULL)
    {
      /* If residue consists of a single atom, check whether it is a metal */
/* v.5.3.7 --> */
//      metal = check_for_metal(dic_data_ptr->atom_name,element);
      metal = check_for_metal(dic_data_ptr->atom_name,
                              dic_data_ptr->element);
/* <-- v.5.3.7 */

      /* Check if this atom is O or S - these potentially accept
         2 hydrogens each */
      if (dic_data_ptr->atom_name[1] == 'O' ||
          dic_data_ptr->atom_name[1] == 'S')
        {
          /* Write out 2 Hbond acceptors */
          fprintf(fil_out,"-e \'%s\' \"%s\" 2\n",
                  res_name,dic_data_ptr->output_atom_name);
        }

      /* If metal, then write out up to 6 potential acceptors */
      else if (metal == TRUE)
        {
          fprintf(fil_out,"-e \'%s\' \"%s\" 8\n",
                  res_name,dic_data_ptr->output_atom_name);
        }

      /* If this atom is N, then more complicated rules apply according
         to connectivities */
      else if (dic_data_ptr->atom_name[1] == 'N')
        {
          /* Case 1: sp2 or aromatic - always accepts 1 */
          if (dic_data_ptr->nbound == 2)
            fprintf(fil_out, "-e \'%s\' \"%s\" 1\n",
                    res_name,dic_data_ptr->output_atom_name);

          /* Case 2: sp3 or amide - need to check bond angles */
          else if (dic_data_ptr->nbound == 3)
            {
              /* Calculate all the bond-angles at this nitrogen position,
                 using the corresponding coordinates from the atoms in the
                 PDB file */

              /* Get the coordinates of the N atom */
              atom_ptr = dic_data_ptr->atom_ptr;

              /* If have the N coordinates, then retrieve them */
              if (atom_ptr != NULL)
                {
                  coords_N[0] = atom_ptr->x;
                  coords_N[1] = atom_ptr->y;
                  coords_N[2] = atom_ptr->z;

                  /* Initialise count of stored coordinates of atoms
                     bonded to this one */
                  nstore = 0;

                  /* Loop through the connected atoms and store all
                     their coordinates */
                  for (jatom = 0; jatom < dic_data_ptr->nbound; jatom++)
                    {
                      /* Get the pointer to the bonded atom */
                      other_dic_data_ptr
                        = dic_data_ptr->connect_atom_ptr[jatom];
                      if (other_dic_data_ptr != NULL)
                        atom_ptr = other_dic_data_ptr->atom_ptr;

                      /* If have a valid coordinate record, then store
                         the coordinates */
                      if (atom_ptr != NULL)
                        {
                          coords_other[nstore][0] = atom_ptr->x;
                          coords_other[nstore][1] = atom_ptr->y;
                          coords_other[nstore][2] = atom_ptr->z;
                          nstore++;
                        }
                    }

                  /* Calculate all possible angles and store */
                  nangles = 0;
                  nplanar = 0;
                  nsp3 = 0;

                  /* Loop over all pairs of stored coordinates */
                  for (istore = 0; istore < nstore - 1; istore++)
                    {
                      for (jstore = istore + 1; jstore < nstore; jstore++)
                        {
                          /* Calculate this angle */
                          angle = calc_angle(coords_other[istore],coords_N,
                                             coords_other[jstore]);
                          nangles++;

                          /* Angle close to 120 => planar, amide bond */
                          if (angle > (120 - ANGLE_ERROR))
                            nplanar++;
                          
                          /* Angle close to 110 => sp3 */
                          else if (angle < (110 + ANGLE_ERROR))
                            nsp3++;
                        }
                    }

                  /* Now use the computed angles to determine the
                     hybridisation state of the nitrogen */

                  /* Determine number of H-bond acceptors according to
                     the angles just calculated */
                  naccept = 0;

                  /* No bond angles within limits */
                  if (nangles == 0)
                    {
                      printf("*** Warning. Cannot determine hybridization");
                      printf(" of %s atom in residue %s\n",
                             dic_data_ptr->output_atom_name,
                             dictionary_ptr->res_name);
                      printf("***          H-bond acceptance set to 1\n");
                      naccept = 1;
                    }

                  /* More angles at 110 than at 120, so assume sp3 */
                  else if (nsp3 > nplanar)
                    naccept = 1;

                  /* Hybridization state ambiguous */
                  else if (nsp3 == nplanar)
                    {
                      printf("*** Warning. Hybridization ambiguous");
                      printf(" for %s atom in residue-type %s\n",
                             dic_data_ptr->output_atom_name,
                             dictionary_ptr->res_name);
                      printf("***          H-bond acceptance set to 1\n");
                      naccept = 1;
                    }

                  /* Write out the number of H-bond acceptors */
                  if (naccept > 0)
                    fprintf(fil_out,"-e \'%s\' \"%s\" %d\n",
                            res_name,dic_data_ptr->output_atom_name,
                            naccept);
                }

              /* Otherwise, print an error message */
              else
                printf("*** Warning. Coords for N missing\n");
            }
        }

      /* Set pointer to next atom */
      dic_data_ptr = dic_data_ptr->next_dic_data_ptr;
      iatom++;
    }
}
/***********************************************************************

open_hbplusrc  -  Open the hbplus.rc file

***********************************************************************/

void open_hbplusrc(char out_name[])
{
  /* Open the output hbplus.rc file */
  printf("\n");
  printf("  Opening HBPLUS input file, hbplus.rc ...\n");
  if ((fil_out = fopen(out_name,"w")) == NULL)
    {
      printf("\n*** Unable to open HBPLUS file [%s] for output\n",out_name);
      exit(1);
    }
}    
/***********************************************************************

generate_hbplusrc  -  Process each correctly-matched dictionary HET
                      group entry and write out a set of hbplus.rc
                      records for it

***********************************************************************/

void generate_hbplusrc(char out_name[],
                       struct dictionary *first_dictionary_ptr)
{
  static int file_open = FALSE;
  struct dictionary *dictionary_ptr;
    
  /* Initialise residue-pointer to start of linked list */
  dictionary_ptr = first_dictionary_ptr;

  /* Loop through all the HET group residues */
  while (dictionary_ptr != NULL)
    {
      /* If this het group has a pointer to a residue from the PDB, then
         PDB entry can write out to hbplus.rc file */
      if (dictionary_ptr->residue_ptr != NULL)
        {
          /* If this is the first entry to be written out, open the
             hbplus.rc file */
          if (file_open == FALSE)
            {
              open_hbplusrc(out_name);
              file_open = TRUE;
            }

          /* Write out the residue and atom names to the hbplus.rc file */
          printf("    Writing H-bonding definitions for residue-type");
          printf(" [%s] ...\n",dictionary_ptr->res_name);
          writename_fn(dictionary_ptr->res_name,dictionary_ptr->natoms,
                       dictionary_ptr->first_dic_data_ptr);

          /* Write out all the connectivities to the hbplus.rc file */
          writecon_fn(dictionary_ptr->res_name,dictionary_ptr->natoms,
                      dictionary_ptr->first_dic_data_ptr);

          /* Calculate numbers of hydrogen donors on each atom */
          writehbdon_fn(dictionary_ptr->res_name,dictionary_ptr->natoms,
                        dictionary_ptr->first_dic_data_ptr,
                        dictionary_ptr->fake_entry,
                        dictionary_ptr->have_hydrogens);

          /* Calculate numbers of hydrogen acceptors on each atom */
          writehbacc_fn(dictionary_ptr->res_name,dictionary_ptr->natoms,
                         dictionary_ptr,dictionary_ptr->first_dic_data_ptr);
        }

      /* Get pointer to the next residue in the list */
      dictionary_ptr = dictionary_ptr->next_dictionary_ptr;
    }

  /* Close the hbplus.rc file */
  if (file_open == TRUE)
    fclose(fil_out);
}
/***********************************************************************

write_bonds  -  Write out the bonds-list

***********************************************************************/

void write_bonds(char *out_name,struct residue *first_residue_ptr)
{
  int bond_order, iatom, ibound, natoms, nbound;
  int wanted;

  struct atom *atom_ptr, *other_atom_ptr;
  struct dic_data *dic_data_ptr, *other_dic_data_ptr;
  struct residue *residue_ptr;

  FILE *file_ptr;

  /* Open the output file */
  printf("  Opening output bonds file, %s ...\n",out_name);
  if ((file_ptr = fopen(out_name,"w")) == NULL)
    {
      printf("\n*** Unable to open output bonds file [%s]\n",out_name);
      exit(1);
    }

  /* Initialise residue-pointer to start of linked list */
  residue_ptr = first_residue_ptr;

  /* Loop through all the residues to write out */
  while (residue_ptr != NULL)
    {
      /* Only write out if this is the first of its kind */
      if (residue_ptr->nmatch == 1)
        {
          /* Get pointer to residue's first atom */
          atom_ptr = residue_ptr->first_atom_ptr;
          natoms = residue_ptr->natoms;

          /* Loop through all the atoms to pick up the missed ones */
          for (iatom = 0; iatom < natoms && atom_ptr != NULL; iatom++)
            {
              /* Get dictionary entry onto which this atom is mapped */
              dic_data_ptr = atom_ptr->dic_data_ptr;

              /* If atom mapped onto dictionary entry, get its bonds */
              if (dic_data_ptr != NULL)
                {
                  /* Get number of connections */
                  nbound = dic_data_ptr->nbound;

                  /* Loop over all the bonded atoms */
                  for (ibound = 0; ibound < nbound; ibound++)
                    {
                      /* Get pointer to connected dictionary entry */
                      other_dic_data_ptr
                        = dic_data_ptr->connect_atom_ptr[ibound];

                      /* Get the atom to which this maps */
                      other_atom_ptr = other_dic_data_ptr->atom_ptr;

                      /* Determine whether valid */
                      wanted = FALSE;
                      if (other_atom_ptr != NULL)
                        {
                          /* Check that atom number higher */
                          if (other_atom_ptr->atom_number
                              > atom_ptr->atom_number)
                            wanted = TRUE;
                        }

                      /* If valid, then write bond out */
                      if (wanted == TRUE)
                        {
                          /* Get the bond order */
                          bond_order = dic_data_ptr->bond_order[ibound];

                          /* Write first atom out */
                          fprintf(file_ptr,"%5d %s %s %s %c ->  ",
                                  atom_ptr->atom_number,
                                  atom_ptr->atom_name,
                                  atom_ptr->residue_ptr->res_name,
                                  atom_ptr->residue_ptr->res_num,
                                  atom_ptr->residue_ptr->chain);

                          /* Write second atom out */
                          fprintf(file_ptr,"%5d %s %s %s %c   ",
                                  other_atom_ptr->atom_number,
                                  other_atom_ptr->atom_name,
                                  other_atom_ptr->residue_ptr->res_name,
                                  other_atom_ptr->residue_ptr->res_num,
                                  other_atom_ptr->residue_ptr->chain);

                          /* Write out the bond order */
                          fprintf(file_ptr,"%s\n",bond_type[bond_order]);
                        }
                    }
                }

              /* Get pointer to next atom in linked-list */
              atom_ptr = atom_ptr->next_atom_ptr;
            }
        }

      /* Get pointer to the next residue in the list */
      residue_ptr = residue_ptr->next_residue_ptr;
    }

  /* Close the output file */
  fclose(file_ptr);
}    
/***********************************************************************

best_dic_match  -  Determine which residue most closely matches the
                   selected dictionary entry

***********************************************************************/

void best_dic_match(char *pdb_code,struct residue *first_residue_ptr,
                    char *out_name,char molec_type,char *reaction_ref,
                    char *pdbsum_dir)
{
  char number_string[20];

  char *string_ptr;

  int best_match, len, total;
  int ec_no, rstep_no, component_no;
/* v.5.0.1--> */
  int ok;
/* <--v.5.0.1 */

/* v.5.0.1--> */
/*  double match_percentage; */
  double match_percentage, percentage;
/* <--v.5.0.1 */

  struct dictionary *dictionary_ptr;
  struct residue *best_residue_ptr, *residue_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  ec_no = 0;
  rstep_no = 0;
  component_no = 0;

  /* Extract the reaction reference */
  strcpy(number_string,reaction_ref);
  string_ptr = strstr(number_string,".");
  if (string_ptr != NULL)
    {
      /* Extract E.C. reference */
      string_ptr[0] = '\0';
      ec_no = atoi(number_string);

      /* Get length of the first number */
      len = strlen(number_string) + 1;

      /* Get next number */
      strcpy(number_string,reaction_ref+len);
      string_ptr = strstr(number_string,".");
      if (string_ptr != NULL)
        {
          /* Extract reaction step number */
          string_ptr[0] = '\0';
          rstep_no = atoi(number_string);

          /* Get length of the first and second numbers */
          len = len + strlen(number_string) + 1;

          /* Get component number */
          strcpy(number_string,reaction_ref+len);
          component_no = atoi(number_string);
        }
    }

  /* Initialise best score */
  best_match = 0;
  best_residue_ptr = NULL;

  /* Initialise residue-pointer to start of linked list */
  residue_ptr = first_residue_ptr;

  /* Loop through all the residues to determine closest match */
  while (residue_ptr != NULL)
    {
      /* Get corresponding dictionary pointer */
      dictionary_ptr = residue_ptr->dictionary_ptr;

      /* If have a dictionary entry, calculate match score */
      if (dictionary_ptr != NULL)
        {
/* v.5.0.1--> */
          /* Determine whether this match is OK to consider */
          if (residue_ptr->nmatched >= MIN_MAPPED ||
              (residue_ptr->nmissed + residue_ptr->nmissing)
              <= MAX_MISSING)
            ok = TRUE;
          else
            ok = FALSE;

          /* Check have reasonable percentage */
          total = residue_ptr->nmatched + residue_ptr->nmissed
            + residue_ptr->nmissing;
          if (total > 0)
            percentage = 100.0 * (double) residue_ptr->nmatched
              / (double) total;
          else
            percentage = 0.0;
          if (percentage < MIN_PERCENT)
            ok = FALSE;
/* <--v.5.0.1 */

          /* If this is the best match score so far, store the pointer */
/* v.5.0.1--> */
/*          if (residue_ptr->nmapped > best_match &&
            residue_ptr->nmapped > MIN_MAPPED) */
          if (residue_ptr->nmatched > best_match && ok == TRUE)
/* <--v.5.0.1 */
            {
              /* Store pointer and best score */
              best_residue_ptr = residue_ptr;
/* v.5.0.1--> */
/*              best_match = residue_ptr->nmapped; */
              best_match = residue_ptr->nmatched;
/* <--v.5.0.1 */
            }
        }

      /* Get pointer to the next residue in the list */
      residue_ptr = residue_ptr->next_residue_ptr;
    }

  /* If have a good mapping between residue and dictionary then
     write out */
  if (best_residue_ptr != NULL)
    {
      /* Open the output file */
      printf("  Opening output map-scores file, %s ...\n",out_name);
      if ((file_ptr = fopen(out_name,"w")) == NULL)
        {
          printf("\n*** Unable to open output file [%s]\n",out_name);
          exit(1);
        }

      /* Calculate match percentage */
      total = best_residue_ptr->nmatched + best_residue_ptr->nmissed
        + best_residue_ptr->nmissing;
      match_percentage = 100.0 * (double) best_residue_ptr->nmatched
        / (double) total;

      /* Get corresponding dictionary pointer */
      dictionary_ptr = best_residue_ptr->dictionary_ptr;

      /* Write out the details */
      fprintf(file_ptr,
              "%s %s %c%3d%3d%3d %s %s %c %6.2f%% %3d %3d %3d "
              "Lig: %3d %3d %s\n",
              pdb_code,best_residue_ptr->dictionary_ptr->res_name,
              molec_type,ec_no,rstep_no,component_no,
              best_residue_ptr->res_name,best_residue_ptr->res_num,
              best_residue_ptr->chain,match_percentage,
              best_residue_ptr->nmatched,best_residue_ptr->nmissed,
              best_residue_ptr->nmissing,best_residue_ptr->liguniq_no,
              best_residue_ptr->ligand_no,dictionary_ptr->mol_name);

      /* Close the output file */
      fclose(file_ptr);
    }
}    
/***********************************************************************

write_hetpdb_files  -  Write out the hetxx.pdb files for PDBsum

***********************************************************************/

void write_hetpdb_files(struct residue *first_residue_ptr)
{
  char file_name[FILENAME_LEN], res_name[4];

  int iatom, iconec, ipos, natoms, nhets;
  int done;

  struct atom *atom_ptr, *other_atom_ptr;
  struct residue *residue_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  printf("  Writing out hetpdb files ...\n");
  nhets = 0;

  /* Initialise residue-pointer to start of linked list */
  residue_ptr = first_residue_ptr;

  /* Loop through all the residues to write out */
  while (residue_ptr != NULL)
    {
      /* Only write out if this is the first of its kind */
      if (residue_ptr->nmatch == 1)
        {
          /* Increment count of het groups */
          nhets++;

          /* Get residue name */
          strcpy(res_name,residue_ptr->res_name);

          /* Replace all blanks by underscores */
          for (ipos = 0; ipos < 3; ipos++)
            if (res_name[ipos] == ' ')
              res_name[ipos] = '_';

          /* Form name of output file */
          strcpy(file_name,"het");
          strcat(file_name,res_name);
          strcat(file_name,".pdb");

          /* Open the output file */
          printf("    Opening output PDB file, %s ...\n",file_name);
          if ((file_ptr = fopen(file_name,"w")) == NULL)
            {
              printf("\n*** Unable to open output PDB file [%s]\n",
                     file_name);
              exit(1);
            }

          /* Write out whether bond details came from the Het Group
             dictionary */
          if (residue_ptr->fake == TRUE)
            fprintf(file_ptr,"set from_dic = FALSE  ##run.scr\n");
          else
            fprintf(file_ptr,"set from_dic = TRUE   ##run.scr\n");

          /* Write out number of atoms */
          fprintf(file_ptr,"@ natoms = %d  ##run.scr\n",
                  residue_ptr->natoms);

          /* Write out residue name, for pick-up by the script */
          fprintf(file_ptr,"set het_name = \"%s\"  ##run.scr\n",res_name);

          /* Write out residue name, number and chain identifier */
          fprintf(file_ptr,"set resname = \"-n%s\"  ##run.scr\n",
                  residue_ptr->res_name);
          fprintf(file_ptr,"set resnum = \"%s\"  ##run.scr\n",
                  residue_ptr->res_num);
          if (residue_ptr->chain == '!')
            fprintf(file_ptr,"set chain = \"\\%c\"  ##run.scr\n",
                    residue_ptr->chain);
          else
            fprintf(file_ptr,"set chain = \"%c\"  ##run.scr\n",
                    residue_ptr->chain);

          /* Get pointer to residue's first atom */
          atom_ptr = residue_ptr->first_atom_ptr;
          natoms = residue_ptr->natoms;

          /* Loop through all the atoms to pick up the missed ones */
          for (iatom = 0; iatom < natoms && atom_ptr != NULL; iatom++)
            {
              /* Write atom out in PDB format */
              fprintf(file_ptr,"ATOM  %5d %4s %3s %c%5s   %8.3f%8.3f%8.3f",
                      atom_ptr->atom_number,
                      atom_ptr->atom_name,
                      atom_ptr->residue_ptr->res_name,
                      atom_ptr->residue_ptr->chain,
                      atom_ptr->residue_ptr->res_num,
                      atom_ptr->x,atom_ptr->y,atom_ptr->z);
              fprintf(file_ptr,"  1.00 51.54          %s\n",
                      atom_ptr->element);

              /* Get pointer to next atom in linked-list */
              atom_ptr = atom_ptr->next_atom_ptr;
            }

          /* Close off the molecule */
          fprintf(file_ptr,"TER   \n");

          /* Write out all the bonds as CONECT records */

          /* Get pointer to residue's first atom again*/
          atom_ptr = residue_ptr->first_atom_ptr;
          natoms = residue_ptr->natoms;

          /* Loop through all the atoms to write out theis connectivities */
          for (iatom = 0; iatom < natoms && atom_ptr != NULL; iatom++)
            {
              /* Initialise flag */
              done = FALSE;

              /* Loop through this atom's connectivities */
              for (iconec = 0; iconec < MAXCONNECT && done == FALSE; iconec++)
                {
                  /* Get pointer to bonded atom */
                  other_atom_ptr = atom_ptr->connect_atom_ptr[iconec];

                  /* If null, then end here */
                  if (other_atom_ptr == NULL)
                    done = TRUE;

                  /* Otherwise, write out CONECT record */
                  else
                    fprintf(file_ptr,"CONECT%5d%5d\n",atom_ptr->atom_number,
                            other_atom_ptr->atom_number);
                }

              /* Get pointer to next atom in linked-list */
              atom_ptr = atom_ptr->next_atom_ptr;
            }

          /* Close the output file */
          fclose(file_ptr);
        }

      /* Get pointer to the next residue in the list */
      residue_ptr = residue_ptr->next_residue_ptr;
    }
}    
/***********************************************************************

                             M   A   I   N

***********************************************************************/

int main(int argc,char *argv[])
{
  char dictionary_name[FILENAME_LEN];
  char bond_name[FILENAME_LEN] = "hbadd.bonds";
  char match_file[FILENAME_LEN] = "match.dat";
  char out_name[FILENAME_LEN] = "hbplus.rc";
  char pdb_name[FILENAME_LEN];
  char pdbsum_dir[FILENAME_LEN];
  char sum_name[FILENAME_LEN] = "hbadd.sum";
  char map_name[FILENAME_LEN] = "hbadd.map";
  char molec_type, dic_res_name[RESNAM_LEN];
  char pdb_code[5], reaction_ref[REF_LEN];
/* v.5.0--> */
  char work_dir[FILENAME_LEN];
/* <--v.5.0 */

  int natoms, nhetgroups, ndic_entries;
  int for_pdbsum, have_hydrogens, match_names;
/* v.5.1--> */
  int len;
  int sdf;
/* <--v.5.1 */

  struct atom *first_atom_ptr;
  struct residue *first_residue_ptr;
  struct connect *first_connect_ptr;
  struct dictionary *first_dictionary_ptr;
  struct dic_data *first_dic_data_ptr;
  struct dictionary *last_dictionary_ptr;
  struct dic_data *last_dic_data_ptr;


  /* Initialise variables */
  fil_sum = NULL;
  pdbsum_dir[0] = '\0';

  /* Initialise pointers */
  first_residue_ptr = NULL;
/* v.5.3.1--> */
  first_dictionary_ptr = last_dictionary_ptr = NULL;
  first_dic_data_ptr = last_dic_data_ptr = NULL;
/* <--v.5.3.1 */

  /* Get the PDB file name from the command-line arguments*/
  get_command_arguments(argv,argc - 1,pdb_name,dictionary_name,
                        &for_pdbsum,dic_res_name,&molec_type,
/* v.5.0--> */
/*                         reaction_ref); */
                        reaction_ref,&match_names,work_dir,pdbsum_dir);

  /* If have a working directory, modify the names of all the input
     and output files to prefix with this directory */
  if (work_dir[0] != '\0')
    prefix_file_names(work_dir,bond_name,match_file,out_name,
                      sum_name,map_name);
/* <--v.5.0 */

  /* Read in all the HET groups in the PDB file */
  read_pdbfile(pdb_name,pdb_code,pdbsum_dir,dic_res_name,&first_atom_ptr,
               &first_residue_ptr,&first_connect_ptr,&natoms,
               &nhetgroups,&have_hydrogens);

  /* If looking for specific dictionary entry match the residues just
     read in to unique ligands, deleting any that have no match */
  if (dic_res_name[0] != '\0')
    match_up_ligands(pdb_code,&first_residue_ptr,pdbsum_dir);

/* v.5.1--> */
  /* Check whether distionary is in SDF format */
  sdf = FALSE;
  len = strlen(dictionary_name);
  if (!strcmp(dictionary_name + len - 4,".sdf") ||
      !strcmp(dictionary_name + len - 4,".sdf"))
    sdf = TRUE;

  /* If dictionary is in SDF format, then read in */
  if (sdf == TRUE)
    ndic_entries
      = read_sdf_dictionary(dictionary_name,first_residue_ptr,dic_res_name,
                            &first_dictionary_ptr,&last_dictionary_ptr,
                            &first_dic_data_ptr,&last_dic_data_ptr);

  /* Otherwise, read in .txt or .cif format dictionary */
  else
/* <--v.5.1 */

    /* Read in any matching HET group definitions from the HET Group
       Dictionary */
    ndic_entries
      = read_dictionary(dictionary_name,first_residue_ptr,dic_res_name,
                        &first_dictionary_ptr,&last_dictionary_ptr,
                        &first_dic_data_ptr,&last_dic_data_ptr);

  /* If have no dictionary entries, then end here */
/* v.5.1.1-->
  if (ndic_entries == 0)
    exit(0);
<--v.5.1.1 */

  /* If looking for a specific dictionary entry to compare against
     all Het Groups, duplicate the first dictionary entry */
  if (dic_res_name[0] != '\0')
    generate_duplicate_entries(first_dictionary_ptr);

  /* Match up all the connected atoms from the data read in from the 
     HET group definitions */
  dic_connectivities(first_dictionary_ptr);

  /* Open the hbadd.map output file */
  open_mapfile(map_name,pdb_name);

  /* Print any HET groups in the PDB file that weren't found in the 
     HET Group Dictionary, or don't quite match */
  validate_pdb_hets(first_residue_ptr,&first_dictionary_ptr,
                    &last_dictionary_ptr,&ndic_entries,
                    &last_dic_data_ptr,first_connect_ptr,sum_name,
/* v.5.0--> */
/*                    pdb_name,dic_res_name); */
                    pdb_name,dic_res_name,match_names);
/* <--v.5.0 */

  /* Process each correctly-matched dictionary HET group entry and write
     out a set of hbplus.rc records for it */
  if (dic_res_name[0] == '\0')
    generate_hbplusrc(out_name,first_dictionary_ptr);

  /* If comparing all Het Groups against a single dictionary entry,
     determine which one gave the closest match */
  else
    best_dic_match(pdb_code,first_residue_ptr,match_file,molec_type,
                   reaction_ref,pdbsum_dir);

  /* Write out the file of bonds */
  write_bonds(bond_name,first_residue_ptr);

  /* If required, write out the hetxx.pdb files for PDBsum */
  if (for_pdbsum == TRUE)
    write_hetpdb_files(first_residue_ptr);

  /* Finish */
  return(0);
}
