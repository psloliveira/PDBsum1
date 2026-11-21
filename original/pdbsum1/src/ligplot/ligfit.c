/***********************************************************************

  ligfit.c - Program to compare set of ligand-protein interactions with
             those from an ensemble of LIGPLOTs (written out by LigPlus)
             and compute the constraint file (ligplot.rcm) to generate a
             LIGPLOT for this new ligand to fit to the existing LIGPLOTs
             in the pdbentry

************************************************************************

Date:-         31 Jul 2009

Written by:-   Roman Laskowski

               (Fitting routines qikfit and matfit written by Andrew Martin,
                based on FORTRAN versions written by Mike Sutcliffe.)

****** Testing

ligfit . -dimplot d1 d2 -antibody -exe ~/src/ligplot -prm /ebi/research/thornton/userapps/roman/LigPlus/lib/params/ligplot.prm > debug

 ligplot dimplot.pdb -norot -dimfit -ligfit -wkdir ./ -prm /c/roman/java/LigPlus/lib/params/ligplot.prm -ctype 1

----------------------------------------------------------------------*/

/* Version history (version numbers as for LigPlus)
   ===============

v.1.0.4  Changes to test whether alignment of binding site residues gives
         a good fit when the structures are superposed.
                                                        RAL 27 May 2010

v.1.0.5  Minor changes to structural superposition cut-off values.
                                                        RAL 14 Jul 2010

v.1.1    Bug fix to transfer HHB, NNB, -HB and -NB records from input PDB
         file to output ligplus.pdb file.
                                                        RAL 26 Jul 2010
v.1.2    Change to ligand fitting to use least-squares superposition to
         fit the flattened ligands on the basis of the equivalenced atoms.
                                                        RAL 17 Sep 2010
         Addition of CORA alignments.
                                                        RAL 29 Sep 2010
v.1.3    Fix so that waters interacting via another water to ligand are
         not excluded.
                                                        RAL  9 Feb 2011
         Change to read in ALL the original PDB entry's coords. For
         standard fitting, check sequence alignment and use whole sequence.
                                                        RAL 10 Feb 2011
v.1.3.1  Reinstate anchor weight for equivalenced atoms that are involved
         in H-bonds.
                                                        RAL 19 Aug 2011
v.1.3.6  Strip off double-quotes from residue number if required.
         Modifications to strcpy(field,field + n) commands which fail
         under some linux operating systems.
                                                        RAL 24 Jan 2012
v.1.3.7  Problem of aspartic protease dimers binding symmetrical (or
         nearly symmetrical) ligands.
         Bug-fix on read-in of FASTA alignments.
                                                        RAL  3 May 2012
v.1.4.1  Bug-fix to save and restore hbadd.bonds file which was being
         overwritten whenever hbadd was be run from within ligfit to
         perform graph-matching and all bond info about ligand was being
         lost.
                                                        RAL  8 Jul 2012
v.2.0    New multiple-dimplot routines to superpose two protein-protein
         interfaces, based on mapping one side onto the other.
                                                        RAL  2 Sep 2014
v.2.0.1  Addition of verbose option for debugging.
         Fix to problem of ligand being skipped if its first atom is a
         hydrogen.
                                                        RAL 16 Feb 2015

v.2.0.2  Bug-fixes for antibody plot crashing.
                                                        RAL 13 Apr 2015

v.2.0.3  Bug-fixes for ligand-fitting failing when ligand contains
         unusually long covalent bond (eg to a bromine).
                                                        RAL  5 Jun 2015

v.2.1    Bug-fix on residue placement not allowing enough space for the
         placed residue, resulting in overlaps.
                                                        RAL 13 Oct 2015

v.2.2.1  Dealing with cases where protein residue has same number and
         chain-id as the ligand.
                                                        RAL 20 Mar 2020

v.2.2.8  Change to system calls to allow for spaces in file names and
         paths. 
                                                        RAL  9 Nov 2022
         Bug-fix preventing antibody plots superposing on macs.
                                                        RAL  7 Jan 2023

*/

/* Datafiles
   ---------

Actual name      Description
-----------      -----------
align.aln        Input CORA file giving structural alignment to be used
                 when overlaying the ligand binding sites
dimhtml.rcm      Input .rcm file giving layout for current plot
dimplot.rcm      Output constraints file to DIMPLOT
ensemble.drw     Input .drw file from LigPlus holding the 2D and 3D
                 coordinates of the existing LIGPLOT(s) to fit new ligand
                 and protein sidechains to
ligand.pdb, ligand.hhb, ligand.nnb, ligand.rcm
                 Output files for input to LIGPLOT to position current
                 ligand relative to previous ligands
ligplot.drw      Input file, created by LIGPLOT for current ligand
ligplus.dom      Input file giving domain definition for DIMPLOT
ligplus.pdb      Input/output PDB file holding coords of the protein-ligand
                 complex to be processed. New version of this file, with
                 coordinates tranformed is written out by program
ligplus.hhb      Input HBPLUS H-bonds file
ligplus.lgf      Output file of residue equivalences
ligplus.nnb      Input HBPLUS non-bonded contacts file
ligplus.rcm      Output constraints file to LIGPLOT
shifts.pdb       Output file of atom coords of any sidechains from previous
                 plot that were shifted out of the way of the current
                 plot's ligand

*/

/* Function calling-tree
   ---------------------

main
   -> get_command_arguments
         -> store_string
         -> get_residue_number
   -> get_domains
         -> get_domain_number
         -> get_residue_range
               -> get_token
               -> format_number
         -> create_domain_record
   -> read_drw_file
         -> string_truncate
         -> create_pdbentry_record
         -> get_drw_residue
               -> create_residue_record
                     -> get_aa_code
         -> get_drw_atom
               -> get_tokens
               -> create_atom_record
                     -> get_element
                           -> is_hydrogen_atom
                                 -> check_for_metal
                           -> check_for_metal
         -> get_drw_bond
               -> get_tokens
               -> create_bond_record
   -> sort_residues
         -> perform_sort
   -> init_filenames
   -> read_pdb_file
         -> create_pdbentry_record
         -> string_truncate
         -> create_hhbnnb_record
         -> get_atom_details
         -> check_for_duplicate
         -> unmark_residues
         -> create_residue_record
               -> get_aa_code
         -> create_atom_record
               -> get_element
                     -> is_hydrogen_atom
                           -> check_for_metal
                     -> check_for_metal
         -> find_atom
         -> calc_distance
         -> create_bond_record
   -> replace_plot_residues
         -> find_residue
         -> find_atom_name
         -> replace_bond_atom
   -> assign_domains
   -> update_all_maxmin
         -> update_max_min
   -> calc_bonds
         -> check_bond
         -> create_bond_record
   -> read_hbplus_file
         -> get_hbplus_info
         -> match_residue
         -> match_atom
         -> create_bond_record
   -> check_interactions
   -> classify_residues
   -> get_cora_alignment
         -> string_truncate
         -> get_tokens
         -> to_lower
         -> get_aa_name
         -> create_resseq_record
         -> store_alignment
   -> get_caf_alignment
         -> string_truncate
         -> match_alignment
               -> get_sequence
               -> create_alignment_arrays
               -> perform_alignment
                     -> fill_array
                     -> nw_sweep
                     -> find_best_path
                     -> get_mapping
               -> match_up_residues
   -> get_fasta_alignment
         -> string_truncate
         -> to_lower
         -> match_alignment (as above)
   -> get_seq_alignment
         -> get_sequence
         -> get_dimplot_sequence
         -> create_alignment_arrays
         -> perform_alignment (as above)
   -> get_bs_alignment
         -> get_residue_names
               -> get_aa_name
         -> create_residue_record
               -> get_aa_code
         -> get_bs_sequence
         -> get_dimplot_sequence
         -> copy_sequence
         -> reverse_chains
         -> create_alignment_arrays
         -> perform_alignment (as above)
         -> get_residue_mappings
               -> create_resseq_record
               -> store_alignment
   -> fit_structures
         -> fit_coords
               -> get_fit_coords
                     -> check_backbone
               -> centre_coords
               -> matfit
                     -> qikfit
               -> calc_separation
                     -> transform_atom
   -> match_ligands
         -> perform_graph_match
               -> create_dummy_dic
                     -> open_file
                           -> get_file_name
                           -> open_output_file
                     -> get_attype
               -> write_dummy_pdb
                     -> open_file
                           -> get_file_name
                           -> open_output_file
                     -> get_attype
               -> copy_file
               -> run_hbadd
               -> read_hbadd_map
                     -> get_file_name
                     -> create_atomap_record
               -> fit_ligands
                     -> mark_ligands
               -> copy_file
   -> map_extras
         -> residue_overlap
               -> find_atomap
               -> create_atomap_record
   -> weigh_anchors
         -> check_equiv_hbonds
   -> write_coords
   -> read_rcm_file
         -> string_truncate (in common.c)
         -> find_residue
         -> create_layout_record
   -> add_dummy_layouts
         -> create_residue_record
         -> create_atom_record
         -> update_max_min
   -> arrange_dimplot
         -> find_layout
         -> place_residue
               -> get_no
         -> adjust_atoms
         -> shift_residues
         -> update_all_maxmin
         -> list_residues
   -> dimplot_restraints
   -> get_equivalences
         -> open_file
               -> open_output_file
         -> write_ligand_restraints
         -> write_equiv_residues
               -> check_swappable_atoms
                     -> score_interactions
         -> write_equiv_ligatoms
               -> find_bond
               -> create_bond_record
               -> find_atomap
               -> create_atomap_record
               -> free_bonds
         -> write_residue_cofm
   -> write_ligand_files
   -> run_ligplot
         -> repeat_splice
               -> perform_splice
   -> shift_sidechains
         -> update_all_maxmin
               -> update_max_min
         -> calc_cofm
         -> check_overlap
         -> make_shift
         -> write_ligfit
         -> write_blocking_ligands
               -> write_residue
   -> map_ligands
   -> write_shifts_pdb
         -> open_file
         -> write_residue

Compile as follows:-
------------------

cc -o ligfit ligfit.c -lm

Optimization flags:-  cc -O2 -funroll-loops -o ligfit ligfit.c -lm

For testing without system calls, use:
    gcc -DTEST -o ligfit ligfit.c -lm

*/

#include "ligfit.h"

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

get_residue_number  -  Interpret the residue-number entered by the
                       user

***********************************************************************/

int get_residue_number(char *string,char *res_num,char *res_name)
{
  char insertion_code, res_number[6], new_number[6];
  char temp_string[LINELEN];

  int cpos, got_number, idigit, ipos, ndigits;
  int at_end, residue_name, residue_number;
  int len, nchars;

  /* Initialise position of insertion code */
  cpos = 0;
  insertion_code = ' ';
  got_number = FALSE;
  ndigits = 0;
  for (idigit = 0; idigit < 5; idigit++)
    new_number[idigit] = ' ';
  at_end = FALSE;
  nchars = 0;
  residue_name = FALSE;
  residue_number = FALSE;

/* v.1.3.6 --> */
  /* Check for quotes around residue name and strip off if present */
  len = strlen(string);
  if (string[len - 1] == '"')
    string[len - 1] = '\0';
  if (string[0] == '"')
    {
      strcpy(temp_string,string+1);
      strcpy(string,temp_string);
    }
/* <-- v.1.3.6 */

  /* Check whether this has a -n prefix indicating a residue name */
  if (!strncmp(string,"-n",2))
    {
      residue_name = TRUE;

      /* Excise the -n */
      strcpy(temp_string,string+2);
      strcpy(string,temp_string);
    }

  /* Loop through the characters in the entered residue-number
     searching for end of number or insertion-code */
  for (ipos = 0; ipos < 5 && at_end == FALSE; ipos++)
    {
      /* If this is the end of the string, store end-position */
      if (string[ipos] == '\0')
        {
          cpos = ipos;
          at_end = TRUE;
        }

      /* If this is a digit, then have part of a residue number */
      else if ((string[ipos] >= '0' && string[ipos] <= '9') ||
               string[ipos] == '-')
        {
          res_number[ndigits] = string[ipos];
          ndigits++;
          got_number = TRUE;

          /* If already have a character, then this must be a residue name */
          if (nchars > 0)
            {
              got_number = FALSE;
              residue_name = TRUE;
            }
        }

      /* If this is a character it might be a prefix to the
         residue number, or the insertion code */
      else if ((string[ipos] >= 'A' && string[ipos] <= 'Z') ||
               string[ipos] == '\'')
        {
          /* If we already have a number, this must be an
             insertion code */
          if (got_number == TRUE)
            {
              cpos = ipos;
              insertion_code = string[ipos];
            }

          /* Otherwise it must be a residue name */
          else
            {
              residue_name = TRUE;
              got_number = FALSE;
            }
          /* Increment count of characters */
          nchars++;
        }

      /* If this is a space, have reached the end of the number */
      else if (string[ipos] == ' ')
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

  /* If entry doesn't resemble a residue number, then see whether it
     might be a residue name instead */
  if (residue_name == TRUE || cpos < 0 || got_number == FALSE ||
      nchars > 1)
    {
      /* Get length of residue name */
      len = strlen(string);

      /* If too long, then have an error */
      if (len > 3)
        {
          printf("*** Invalid residue-number or name:");
          printf(" [%s]\n",string);
          exit(1);
        }

      /* Store the name */
      strcpy(res_name,string);

      /* Check for blanks at end */
      if (len == 3 && !strncmp(string+1,"  ",2))
        {
          len = 1;
          string[1] = '\0';
        }
      if (len == 3 && string[2] == ' ')
        {
          len = 2;
          string[2] = '\0';
        }

      /* Check for short residue names */
      if (len == 1)
        {
          strcpy(res_name,"  ");
          strcat(res_name,string);
        }
      else if (len == 2)
        {
          strcpy(res_name," ");
          strcat(res_name,string);
        }

      residue_number = FALSE;
    }

  /* Otherwise, reformat the number */
  else
    {
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

      /* Set flag to say that residue number is OK */
      residue_number = TRUE;
    }

  /* Return flag indicating whether residue number identified */
  return(residue_number);
}
/***********************************************************************

get_command_arguments  -  Interpret the command line arguments

***********************************************************************/

void get_command_arguments(char *argv[],int num,struct parameters *params)
{
/* v.2.2.8 --> */
//  char token[LINELEN];
  char tmp[LINELEN], token[LINELEN];
/* <-- v.2.2.8 */

  int iarg, irange, res_num;
/* v.2.2.8 --> */
  int i;
  int done;
/* <-- v.2.2.8 */

  /* Initialise variables */
  irange = 0;

  /* Initialise parameters */
  params->contact_type = HYDROPHOBIC_ANY;
  params->include_waters = FALSE;
/* v.1.2 --> */
  params->fit_ligands = FALSE;
/* <-- v.1.2 */
  params->directory = &null;
  params->prmfile = &null;
  params->alnfile = &null;
  store_string(&(params->ligplot_exe_dir),".");
  params->ligplot_params[0] = '\0';
  params->res_num1[0] = params->res_num2[0] = '\0';
  params->res_name1[0] = params->res_name2[0] = '\0';
/* v.2.0 --> */
//  params->chain = ' ';
  params->chain[0] = params->chain[1] = ' ';
  params->have_domains = FALSE;
  params->antibody = FALSE;
  params->order = FALSE;
  params->dimplot = FALSE;
/* <-- v.2.0 */
/* v.2.0.1 --> */
  params->verbose = FALSE;
/* <-- v.2.0.1 */

  /* If nothing entered on the command line, show options available */
  if (num < 2)
    {
      printf("\n");
      printf("*** ERROR. Correct usage is:\n");
      printf("\n");
/* v.2.0 --> */
//      printf("    ligfit  directory  [-w]  ligand_defntn "
      printf("    ligfit  directory  [-w]  [-dimplot X Y]  ligand_defntn "
/* <-- v.2.0 */
             "[-exe exe_directory]\n");
/* v.2.0.1 --> */
//      printf("            [-prm prm_file]  [-aln aln_file]  [-fl]\n");
      printf("            [-prm prm_file]  [-aln aln_file]  [-fl]  [-v]\n");
/* <-- v.2.0.1 */
      printf("\n");
      printf("where\n");
      printf("     * directory     = Directory containing all the input "
             "and output files\n");
      printf("     * -w            = Flag indicating waters are to be "
             "included\n");
      printf("     * -fl           = Flag indicating plot to be generated "
             "by fitting the ligands\n");
      printf("     * ligand_defntn = Definition of ligand in first of "
             "the structures\n");
      printf("          - Format:\n");
      printf("                [resname1] resnum1 [resname2] resnum2 "
             "[chain]\n");
      printf("            If residue name is numeric, prefix with -n\n");
/* v.2.0 --> */
      printf("     * -dimplot X Y  = DIMPLOT of chains X and Y, or domains d1 "
             "and d2 (defined in ligplus.dom)\n");
      printf("     * -antibody     = Flag indicating this is for an "
             "antibody plot\n");
      printf("     * -order        = Flag indicating residues in order in "
             "interface plot\n");
/* <-- v.2.0 */
      printf("     * -exe exe_directory = Full pathname of LIGPLOT "
             "executables\n");
      printf("     * -prm prm_file = Full pathname to ligplot.prm file\n");
      printf("     * -aln aln_file = CORA-format structural alignment "
             "file\n");
/* v.2.0.1 --> */
      printf("     * -v            = Verbose option for debugging\n");
/* <-- v.2.0.1 */
      printf("\n");
      printf("\n");
      printf("Examples:-\n");
      printf("\n");
      printf("   ligfit /tmp/run ECA 1 GLU 3 I\n");
      printf("   ligfit /tmp/run -n888 1\n");
      printf("\n");
      exit(1);
    }

  /* Get the name of the working directory from the first paremeter */
  store_string(&(params->directory),argv[1]);

  /* Loop through the command arguments */
  for (iarg = 2; iarg < num + 1; iarg++)
    {
      /* If we have the -w option, set the include waters flag */
      if (!strcmp(argv[iarg],"-w"))
        params->include_waters = TRUE;

      /* If we have the -fl option, set the fit-ligands flag */
      else if (!strcmp(argv[iarg],"-fl"))
        params->fit_ligands = TRUE;

/* v.2.0 --> */
      /* If we have the -dimplot option, set flag */
      else if (!strcmp(argv[iarg],"-dimplot"))
        {
          params->dimplot = TRUE;

          /* Get the two chain-ids from the next two token, if they exist */
          if (iarg < num - 1)
            {
              strcpy(token,argv[iarg + 1]);
              if (strlen(token) == 1)
                params->chain[0] = token[0];
              else if (token[0] == 'd')
                {
                  params->chain[0] = token[1];
                  params->have_domains = TRUE;
                }
              iarg++;
              strcpy(token,argv[iarg + 1]);
              if (strlen(token) == 1)
                params->chain[1] = token[0];
              else if (token[0] == 'd')
                {
                  params->chain[1] = token[1];
                  params->have_domains = TRUE;
                }
              iarg++;
            }
        }

      /* If we have the -order option, set the flag */
      else if (!strcmp(argv[iarg],"-order"))
        params->order = TRUE;

      /* If we have the -antibody option, set the flags */
      else if (!strcmp(argv[iarg],"-antibody"))
        {
          params->antibody = TRUE;
          params->order = TRUE;
        }
/* <-- v.2.0 */

      /* If we have the -exe option, save the full path of the LIGPLOT
         executable directory */
      else if (!strcmp(argv[iarg],"-exe"))
        {
/* v.2.2.8 --> */
          /* Get the path from the next token, if there is one */
/*          if (iarg < num)
            {
              store_string(&(params->ligplot_exe_dir),argv[iarg + 1]);
              iarg++;
            }
*/
          /* Initialise flag */
          done = FALSE;
          tmp[0] = '\0';

          /* Loop until all parts of the path have been obtained */
          for (i = iarg + 1; i < num && done == FALSE; i++)
            {
              /* If not a parameter, add to string */
              if (!strcmp(argv[i],"-") || argv[i][0] != '-')
                {
                  /* Append or start */
                  if (tmp[0] == '\0')
                    strcpy(tmp,argv[i]);
                  else
                    {
                      strcat(tmp," ");
                      strcat(tmp,argv[i]);
                    }
                }

              /* Otherwise, we are done */
              else
                {
                  done = TRUE;
                  iarg = i - 1;
                }

              /* Save the path */
              store_string(&(params->ligplot_exe_dir),tmp);
            }
/* <-- v.2.2.8 */
        }

      /* If we have the -prm option, save the full path of the ligplot.prm
         parameter file */
      else if (!strcmp(argv[iarg],"-prm"))
        {
          /* Get the file name from the next token, if there is one */
          if (iarg < num)
            {
              store_string(&(params->prmfile),argv[iarg + 1]);
              iarg++;
            }
        }

      /* If we have the -aln option, save the full path of the CORA-format
         alignment file */
      else if (!strcmp(argv[iarg],"-aln"))
        {
          /* Get the file name from the next token, if there is one */
          if (iarg < num)
            {
              store_string(&(params->alnfile),argv[iarg + 1]);
              iarg++;
            }
        }

/* v.2.0.1 --> */
      /* If we have the -v option, set to verbose mode */
      else if (!strcmp(argv[iarg],"-v"))
        params->verbose = TRUE;
/* <-- v.2.0.1 */

      /* Otherwise, take this to be the ligand residue range */
      else
        {
          /* Store in LIGPLOT parameters string */
          if (params->ligplot_params[0] == '\0')
            strcpy(params->ligplot_params,argv[iarg]);
          else
            {
              strcat(params->ligplot_params," ");
              strcat(params->ligplot_params,argv[iarg]);
            }

          /* If first token of residue range get residue name/number */
          if (irange == 0)
            {
              /* See if this is a residue number */
              strcpy(token,argv[iarg]);
              res_num
                = get_residue_number(token,params->res_num1,
                                     params->res_name1);

              /* If not a residue number, then decrement range as
                 still need to pick up the residue number */
              if (res_num == FALSE)
                irange--;
            }

          /* If this is the second of the residue-range parameters,
             then take it to be the second residue in the range */
          else if (irange == 1)
            {
              /* Interpret the entered residue number */
              strcpy(token,argv[iarg]);
              res_num
                = get_residue_number(token,params->res_num2,
                                     params->res_name2);

              /* If not a residue number, then decrement range as
                 still need to pick up the residue number */
              if (res_num == FALSE)
                irange--;
            }

          /* The third residue-range parameter will be the chain-id */
          else if (irange == 2)
/* v.2.0 --> */
//            params->chain = argv[iarg][0];
            params->chain[0] = argv[iarg][0];
/* <-- v.2.0 */
          
          /* Increment count of residue-range parameters */
          irange++;
        }
    }

/* v.2.0 --> */
  /* For LIGPLOT, check the ligand residue range */
  if (params->dimplot == FALSE)
    {
/* <-- v.2.0 */
      /* If only one residue entered for the range, then make the
         end-residue number the same */
      if (params->res_num2[0] == '\0')
        {
          /* If have a chain-id in the residue name, store it */
          if (strlen(params->res_name2) == 3 &&
              !strncmp(params->res_name2,"  ",2))
/* v.2.0 --> */
//            params->chain = params->res_name2[2];
            params->chain[0] = params->res_name2[2];
/* <-- v.2.0 */

          /* Copy residue name and number */
          strcpy(params->res_num2,params->res_num1);
          strcpy(params->res_name2,params->res_name1);
        }

      /* If no ligand residue range, then abort */
      if (params->res_num1[0] == '\0')
        {
          printf("*** ERROR. No ligand range entered!\n");
          exit(-1);
        }

      /* Show entered parameters */
      printf("Ligand range:           [%s.%s.%c] - [%s.%s.%c]\n",
/* v.2.0 --> */
//             params->res_name1,params->res_num1,params->chain,
//             params->res_name2,params->res_num2,params->chain);
             params->res_name1,params->res_num1,params->chain[0],
             params->res_name2,params->res_num2,params->chain[1]);
/* <-- v.2.0 */
      printf("\n");
      printf("Ligand residue range:  ");
      if (params->res_name1[0] != '\0')
        printf("%s ",params->res_name1);
      printf("%s to ",params->res_num1);
      if (params->res_name2[0] != '\0')
        printf("%s ",params->res_name2);
      printf("%s   ",params->res_num2);
/* v.2.0 --> */
//      if (params->chain != ' ')
//        printf("Chain [%c]\n",params->chain);
      if (params->chain[0] != ' ')
        printf("Chain [%c]\n",params->chain[0]);
/* <-- v.2.0 */
      else
        printf("\n");
/* v.2.0 --> */
    }

  /* Check the two DIMPLOT chains/domains */
  else
    {
      /* If no chains defined, then abort */
      if (params->chain[0] == ' ' || params->chain[1] == ' ')
        {
          printf("*** ERROR. No chains/domains entered!\n");
          exit(-1);
        }

      /* Show the chains/domains selected */
      if (params->have_domains == TRUE)
        printf("DIMPLOT for domains %c and %c\n",params->chain[0],params->chain[1]);
      else
        printf("DIMPLOT for chains %c and %c\n",params->chain[0],params->chain[1]);
    }
/* <-- v.2.0 */

}
/***********************************************************************

string_truncate  -  Truncate the given string at the last non-blank
                    character, or at its maximum length

***********************************************************************/

int string_truncate(char *string,int max_length)
{
  char ch;

  int iend, ipos;

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
    string[iend + 1] = '\0';
  else
    string[max_length - 1] = '\0';

  /* Return the string length */
  return(strlen(string));
}
/* v.2.0 --> */
/***********************************************************************

get_domain_number - Get the domain number from the current line

**********************************************************************/         

int get_domain_number(char *line,int *ipos)
{
  char ch, number_string[NUMBER_LEN + 1];

  int domain, i, len;
  int done;

  /* Initialise domain number */
  domain = -1;
  done = FALSE;
  i = 0;
  len = strlen(line);

  /* Loop until domain number found */
  while (done == FALSE && *ipos < len)
    {
      /* Get the next character in the string */
      ch = line[*ipos];

      /* If have a number, then add it to the number-string */
      if (ch >= '0' && ch <= '9')
        {
          /* Add to number-string provided not too long */
          if (i < NUMBER_LEN)
            {
              number_string[i] = ch;
              i++;
            }
        }

      /* Otherwise, take this to be the end of the domain-number */
      else if (i > 0)
        {
          /* Extract the domain number */
          number_string[i] = '\0';
          domain = atoi(number_string);
          done = TRUE;
        }

      /* Go to the next character in the line */
      (*ipos)++;
    }

  /* Return the domain number */
  return(domain);
} 
/***********************************************************************

get_token - Get the domain number from the current line

**********************************************************************/         

int get_token(char *line,int *ipos,char *token,int *done)
{
  char ch;

  int i, len;
  int have_token;

  /* Initialise variables */
  *done = FALSE;
  have_token = FALSE;
  i = 0;
  len = strlen(line);
  token[0] = '\0';

  /* Loop until domain number found */
  while (*done == FALSE && have_token == FALSE && *ipos < len)
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
  if (i > 0)
    {
      have_token = TRUE;
      *done = FALSE;
    }

  /* Return flag */
  return(have_token);
} 
/***********************************************************************

format_number - Right-justify residue number, taking into account whether
                it includes an insertion code

**********************************************************************/         

void format_number(char *residue_number,char *token)
{
  char ch;

  int ipos, len, start_pos;

  /* Initialise variables */
  strncpy(residue_number,"     ",5);
  residue_number[5] = '\0';

  /* Get length of token-string */
  len = strlen(token);

  /* Proceed if have a token to process */
  if (len > 0)
    {
      /* If token is maximum length, then transfer as is */
      if (len >= 5)
        {
          /* Set start-position at start of token */
          start_pos = 0;
          len = 5;
        }

      /* Otherwise, determine how number is to be transferred */
      else
        {
          /* Check whether last character is a number or an insertion code */
          ch = token[len - 1];

          /* If it is a number, then add a space onto the end */
          if (ch >= '0' && ch <= '9')
            {
              token[len] = ' ';
              token[len + 1] = '\0';
              len++;
            }

          /* Calculate start-position for transfer */
          start_pos = 5 - len;
        }

      /* Perform the transfer of the residue number to justify it as
         required */
      for (ipos = 0; ipos < len; ipos++)
        residue_number[start_pos + ipos] = token[ipos];
      residue_number[5] = '\0';
    }
} 
/***********************************************************************

get_residue_range - Get the next residue range from the current line

**********************************************************************/         

int get_residue_range(char *line,int *ipos,char *start_residue,
                      char *end_residue,char *chain,int *end)
{
  char residue_number[6], token[TOKEN_LEN + 1];

  int ipstn, len;
  int done, have_range, have_token;

  /* Initialise chain-id */
  *chain = ' ';
  *end = FALSE;
  start_residue[0] = '\0';
  end_residue[0] = '\0';
  done = FALSE;
  have_range = FALSE;
  ipstn = *ipos;
  len = strlen(line);
  if (ipstn > len -1)
    {
      *end = TRUE;
      return(have_range);
    }

  /* Loop until domain number found */
  while (have_range == FALSE && done == FALSE && *ipos < len)
    {
      /* Get the next token from the line */
      have_token = get_token(line,&ipstn,token,&done);

      /* If have a valid token, then extract appropriate part of range */
      if (have_token == TRUE)
        {
          /* Check that haven't found a range separator */
          if (!strncmp(token,"&",1))
            have_range = TRUE;

          /* Otherwise, store appropriate part of range */
          else if (strncmp(token,"-",1))
            {
              /* If don't have a start of range, then take this to be it */
              if (start_residue[0] == '\0')
                {
                  /* Format the residue number and store */
                  format_number(residue_number,token);
                  strcpy(start_residue,residue_number);
                }

              /* If don't have an end of range, then take this to be it */
              else if (end_residue[0] == '\0')
                {
                  /* Format the residue number and store */
                  format_number(residue_number,token);
                  strcpy(end_residue,residue_number);
                }

              /* Otherwise, take this to be the chain-id */
              else
                {
                  *chain = token[0];
                  have_range = TRUE;
                }
            }
        }
    }

  /* Check whether have a range */
  if (start_residue[0] != '\0' && end_residue[0] != '\0')
    have_range = TRUE;

  /* Return current position in the line */
  *ipos = ipstn;
  *end = done;

  /* Return the done flag */
  return(have_range);
} 
/***********************************************************************

create_domain_record  -  Create a new domain record

***********************************************************************/

struct domain *create_domain_record(struct domain **fst_domain_ptr,
                                    struct domain **lst_domain_ptr,
                                    int domain_number,char chain,
                                    char *start_residue_number,
                                    char *end_residue_number)
{
  struct domain *domain_ptr;
  struct domain *first_domain_ptr, *last_domain_ptr;

  /* Initialise domain pointers */
  domain_ptr = NULL;
  last_domain_ptr = *lst_domain_ptr;
  first_domain_ptr = *fst_domain_ptr;

  /* Create domain entry */
  domain_ptr = (struct domain*)malloc(sizeof(struct domain));
  if (domain_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for");
      printf(" struct domain\n");
      printf("***        Program terminated with error.\n");
      exit (1);
    }

  /* Store the details for this entry */
  domain_ptr->domain_number = domain_number;
  domain_ptr->chain = chain;
  strcpy(domain_ptr->start_residue_number,start_residue_number);
  strcpy(domain_ptr->end_residue_number,end_residue_number);

  /* Initialise pointer to next domain record */
  domain_ptr->next_domain_ptr = NULL;

  /* If this is the first entry stored, then update
     pointer to head of linked list */
  if (first_domain_ptr == NULL)
    first_domain_ptr = domain_ptr;

  /* Otherwise, make previous entry point to the current
     one */
  else
    last_domain_ptr->next_domain_ptr = domain_ptr;
                      
  /* Save the current domain's pointer */
  last_domain_ptr = domain_ptr;

  /* Return current pointers */
  *fst_domain_ptr = first_domain_ptr;
  *lst_domain_ptr = last_domain_ptr;
  return(domain_ptr);
}
/***********************************************************************

get_domains - Read in the domain definitions from the .dom file

**********************************************************************/         

void get_domains(char *dom_name,int domain1,int domain2,
/* v.2.0 --> */
//                 struct domain **fst_domain_ptr)
                 struct domain **fst_domain_ptr,struct parameters params)
/* <-- v.2.0 */
{
/* v.2.0 --> */
//  char line[LINELEN + 1];
  char file_name[FILENAME_LEN], line[LINELEN + 1];
/* <-- v.2.0 */
  char chain, start_residue[6], end_residue[6];
/* v.2.0 --> */
  char start_chain, end_chain;
/* <-- v.2.0 */

  int domain, ipos;
  int done, have_domain1, have_domain2, have_range;
 
  struct domain *domain_ptr;
  struct domain *first_domain_ptr, *last_domain_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  *fst_domain_ptr = first_domain_ptr = last_domain_ptr = NULL;
  have_domain1 = FALSE;
  have_domain2 = FALSE;
/* v.2.0 --> */
  start_chain = end_chain = ' ';

  /* Form filename of the domain file */
  strcpy(file_name,params.directory);
  strcat(file_name,DIR_SEP);
  strcat(file_name,dom_name);
/* <-- v.2.0 */

  /* Open the .dom domains file in the current working directory */
/* v.2.0 --> */
//  if ((file_ptr = fopen(dom_name,"r")) == NULL)
//    {
//      printf("*** Warning. Unable to open domain file %s ...\n",dom_name);
//      return;
//    }
  if ((file_ptr = fopen(file_name,"r")) == NULL)
    {
      printf("*** Warning. Unable to open domain file %s ...\n",file_name);
      return;
    }
/* <-- v.2.0 */

  /* Prepare to start reading in the data from the DOM file */
  printf("\n");
  printf("  Reading through .dom file ...\n\n");
  printf("    FILE: %s\n",file_name);
  fflush(stdout);

  /* Loop while reading in the domain data */
  while (fgets(line,LINELEN,file_ptr) !=NULL)      
    {
      /* Check if valid record */ 
      if (!strncmp(line,"DOMAIN",6) || !strncmp(line,"Domain",6) ||
          !strncmp(line,"domain",6))
        {
          /* Truncate the string */
          string_truncate(line,LINELEN);

          /* Get the domain number */
          ipos = 6;
          domain = get_domain_number(line,&ipos);

          /* Initialise flag */
          done = FALSE;

          /* Check whether this is one of the domains required */
          if (domain == domain1)
            have_domain1 = TRUE;
          else if (domain == domain2)
            have_domain2 = TRUE;
          else
            done = TRUE;

          /* If have a valid domain number, then extract all the residue
             details */
          while (done == FALSE && domain > -1)
            {
              /* Get the next residue-range in the current line */
              have_range = get_residue_range(line,&ipos,start_residue,
                                             end_residue,&chain,&done);

              /* If have a residue-range, then store */
              if (have_range == TRUE)
                {
                  /* Create a new domain record */
                  domain_ptr
                    = create_domain_record(&first_domain_ptr,&last_domain_ptr,
                                           domain,chain,start_residue,
                                           end_residue);
                }
            }
        }
    }

  /* Close the .dom file */
  fclose(file_ptr);

  /* Check whether have definitions for both domains */
  if (have_domain1 == FALSE || have_domain2 == FALSE)
    {
      /* Print error message if first domain undefined */
      if (have_domain1 == FALSE)
        printf("*** Error. Domain number %d undefined in .dom file\n",
               domain1);

      /* Print error message if second domain undefined */
      if (have_domain2 == FALSE)
        printf("*** Error. Domain number %d undefined in .dom file\n",
               domain2);

      /* Terminate the program */
      printf("***        Program aborted\n");     
      exit(1);
    }

  /* Return pointer to the first domain */
  *fst_domain_ptr = first_domain_ptr;
}  
/* <-- v.2.0 */
/***********************************************************************

create_pdbentry_record  -  Create and initialise a new pdbentry record

***********************************************************************/

struct pdbentry *create_pdbentry_entry(struct pdbentry **fst_pdbentry_ptr,
                                       struct pdbentry **lst_pdbentry_ptr)
{
  struct pdbentry *pdbentry_ptr, *first_pdbentry_ptr, *last_pdbentry_ptr;

  /* Initialise pdbentry pointers */
  pdbentry_ptr = NULL;
  first_pdbentry_ptr = *fst_pdbentry_ptr;
  last_pdbentry_ptr = *lst_pdbentry_ptr;

  /* Allocate memory for structure to hold pdbentry info */
  pdbentry_ptr = (struct pdbentry *) malloc(sizeof(struct pdbentry));
  if (pdbentry_ptr == NULL)
    {
      printf("*** Can't allocate memory for struct pdbentry\n");
      exit (1);
    }

  /* If this is the very first, save its pointer */
  if (first_pdbentry_ptr == NULL)
    first_pdbentry_ptr = pdbentry_ptr;

  /* Add link from previous pdbentry to the current one */
  if (last_pdbentry_ptr != NULL)
    last_pdbentry_ptr->next_pdbentry_ptr = pdbentry_ptr;
  last_pdbentry_ptr = pdbentry_ptr;

  /* Initialise the elements of the structure */
  pdbentry_ptr->pdb_code[0] = '\0';
  pdbentry_ptr->chain = ' ';
/* v.1.3.7 --> */
  pdbentry_ptr->map_chain = ' ';
/* <-- v.1.3.7 */
  pdbentry_ptr->nresid = 0;
  pdbentry_ptr->nbonds = 0;
/* v.1.1 --> */
  pdbentry_ptr->first_hhbnnb_ptr = NULL;
/* <-- v.1.1 */
  pdbentry_ptr->first_residue_ptr = NULL;
  pdbentry_ptr->first_bond_ptr = NULL;
  pdbentry_ptr->pdb_id = -1;

  /* Initialise pointer to the next pdbentry */
  pdbentry_ptr->next_pdbentry_ptr = NULL;

  /* Return current pointers */
  *fst_pdbentry_ptr = first_pdbentry_ptr;
  *lst_pdbentry_ptr = last_pdbentry_ptr;

  /* Return new pdbentry pointer */
  return(pdbentry_ptr);
}
/***********************************************************************

get_line_type  -  Get the .drw line type according to the keyword given
                  by #x, where x is one of the hash codes: D, M, C, S,
                  R, A, B, E, O

***********************************************************************/

int get_line_type(char *input_line)
{
  char code;

  int ltype = NOT_WANTED;

  /* Get the hash-code character */
  code = input_line[1];

  /* Determine which hash code this is (only want molecules,
     residues, atoms and bonds and ensemble records) */
  if (code == 'P')
    {
      ltype = GLOBAL_PARAMS;
    }
  else if (code == 'E')
    {
      ltype = NEW_PDB_ENTRY;
    }
  else if (code == 'M')
    {
      ltype = COLOURMAP_DATA;
    }
  else if (code == 'D')
    {
      ltype = DIMPLOT_FLAG;
    }
  else if (code == 'C')
    {
      ltype = COLOUR_DATA;
    }
  else if (code == 'S')
    {
      ltype = SIZE_DATA;
    }
  else if (code == 'O')
    {
      ltype = MOLECULE_DATA;
    }
  else if (code == 'R')
    {
      ltype = RES_DATA;
    }
  else if (code == 'A')
    {
      ltype = ATOM_DATA;
    }
  else if (code == 'B')
    {
      ltype = BOND_DATA;
    }
  else if (code == 'T')
    {
      ltype = PLOT_TITLE;
    }

  /* Return the line type */
  return(ltype);
}
/***********************************************************************

get_aa_code  -  Get single-letter amino acid code from the residue-name 

***********************************************************************/

char get_aa_code(char res_name[4])
{
  char aa_code;

  int iamino, ires;
  int namino;

  static char amino[] = "ACDEFGHIKLMNPQRSTVWYXBZX?X";
  static char *amino_name[] = {
    "ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU",
    "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR",
    "PCA ","ASX ","GLX ", "UNK", "DUM", "XXX"
  };

  /* Initialise variables */
  namino = strlen(amino) - 1;

  /* Get the corresponding one-letter code */
  iamino = namino;
  for (ires = 0; ires < namino; ires++)
    if (!strncmp(res_name,amino_name[ires],3))
      iamino = ires;
  aa_code = amino[iamino];

  /* Return the one-letter code */
  return(aa_code);
}
/***********************************************************************

create_residue_record  -  Create and initialise a new residue record

***********************************************************************/

struct residue *create_residue_record(struct residue **fst_residue_ptr,
                                      struct residue **lst_residue_ptr,
                                      char *res_name,char *res_num,
                                      char chain,int type,char *label,
                                      float label_x,float label_y,
                                      char *residue_id,
                                      struct pdbentry *pdbentry_ptr)
{
  int i, j;

  struct residue *residue_ptr;
  struct residue *last_residue_ptr, *first_residue_ptr;

  /* Initialise residue pointers */
  residue_ptr = NULL;
  last_residue_ptr = *lst_residue_ptr;
  first_residue_ptr = *fst_residue_ptr;

  /* Allocate memory for structure to hold sequence info */
  residue_ptr = (struct residue *) malloc(sizeof(struct residue));
  if (residue_ptr == NULL)
    {
      printf("*** Can't allocate memory for struct");
      printf(" residue\n");
      exit (1);
    }

  /* Store the details for this entry */
  strcpy(residue_ptr->res_name,res_name);
  strcpy(residue_ptr->res_num,res_num);
  residue_ptr->chain = chain;
  residue_ptr->type = type;
  residue_ptr->pdbentry_ptr = pdbentry_ptr;

  /* Store label and coords, if we have one */
  if (label[0] != '\0')
    {
      /* Store residue label and coords */
      store_string(&(residue_ptr->label),label);
      residue_ptr->label_x = label_x;
      residue_ptr->label_y = label_y;
    }
  else
    {
      residue_ptr->label = &null;
      residue_ptr->label_x = 0.0;
      residue_ptr->label_y = 0.0;
    }

  /* Get the single-letter amino acid code */
  residue_ptr->aa_code = get_aa_code(res_name);

  /* Initialise the remaining residue details */
  for (i = 0; i < 2; i++)
    for (j = 0; j < 3; j++)
      residue_ptr->valmin[i][j] = residue_ptr->valmax[i][j] = 0;
  residue_ptr->interaction = NO_INTERACTION;
  residue_ptr->natoms = 0;
  residue_ptr->nmap = 0;
  residue_ptr->shifted = FALSE;
  store_string(&(residue_ptr->residue_id),residue_id);
/* v.1.3 --> */
  residue_ptr->in_plot = FALSE;
  residue_ptr->aln_pos = -1;
/* <-- v.1.3 */
/* v.2.0 --> */
  residue_ptr->chain_type = NON_INTERFACE;
  residue_ptr->placed = FALSE;
/* <-- v.2.0 */
  residue_ptr->first_atom_ptr = NULL;
  residue_ptr->map_residue_ptr = NULL;
  residue_ptr->next_residue_ptr = NULL;

  /* If this is the first entry stored, then update pointer to head
     of linked list */
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

get_drw_residue  -  Extract residue details from the given .drw line

***********************************************************************/

struct residue *get_drw_residue(char *input_line,char *residue_id,
                                struct residue **fst_residue_ptr,
                                struct residue **lst_residue_ptr,
/* v.2.0 --> */
//                                struct pdbentry *pdbentry_ptr)
                                struct pdbentry *pdbentry_ptr,
                                int molecule,struct parameters params)
/* <-- v.2.0 */
{
  char chain, res_name[4], res_num[6], type[2];
  char label[LINELEN];

  int len, res_type;

  float label_x, label_y;

  struct residue *residue_ptr, *first_residue_ptr, *last_residue_ptr;

  /* Initialise variables */
  chain = ' ';
  label[0] = '\0';
  label_x = 0.0;
  label_y = 0.0;

  /* Initialise pointers */
  residue_ptr = NULL;
  first_residue_ptr = *fst_residue_ptr;
  last_residue_ptr = *lst_residue_ptr;

  /* Get the line length */
  len = strlen(input_line);

  /* Proceed only if line is long enough to contain the required data */
  if (len > 9)
    {
      /* Get the residue type */
      type[0] = input_line[0];
      type[1] = '\0';
      res_type = atoi(type); 

/* v.2.0 --> */
      /* For DIMPLOT, reassign the 'ligand' chain */
      if (params.dimplot == TRUE && res_type == LIGAND)
        res_type = HGROUP;
/* <-- v.2.0 */

      /* Extract the residue details */
      strncpy(res_name,input_line+1,3);
      res_name[3] = '\0';
      strncpy(res_num,input_line+4,5);
      res_num[5] = '\0';

      /* Extract chain identifier */
      chain = input_line[9];
      if (chain == '-')
        chain = ' ';

      /* If have residue name details, then store */
      if (len > 15)
        {
          /* Extract x- and y-coordinates for residue label */
          sscanf(input_line+11,"%f %f %s",&label_x,&label_y,label);
        }

      /* Create a new residue record */
      residue_ptr
        = create_residue_record(&first_residue_ptr,&last_residue_ptr,
                                res_name,res_num,chain,res_type,label,
                                label_x,label_y,residue_id,pdbentry_ptr);

/* v.2.0 --> */
      /* For DIMPLOT, determine the chain type */
      if (params.dimplot == TRUE)
        {
          /* Check if this is a water, or either interface */
          if (!strcmp(res_name,"HOH"))
            residue_ptr->chain_type = INTERFACE_WATER;
          else if (molecule == 1)
            residue_ptr->chain_type = CHAIN1;
          else if (molecule == 2)
            residue_ptr->chain_type = CHAIN2;
        }
/* <-- v.2.0 */
    }

  /* Return current pointers */
  *fst_residue_ptr = first_residue_ptr;
  *lst_residue_ptr = last_residue_ptr;

  /* Return residue */
  return(residue_ptr);
}
/***********************************************************************

check_for_metal  -  Check whether the given atom is a metal ion

***********************************************************************/

int check_for_metal(char *atom_name,char *element)
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
      if (!strncmp(atom_name,metal_name[imetal],2) &&
          (atom_name[2] == ' ' ||
           (atom_name[2] >= '0' && atom_name[2] <= '9')))
        {
          /* Set flag */
          metal = TRUE;

          /* Return the metal name matched */
          strncpy(element,metal_name[imetal],2);
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

int is_hydrogen_atom(char atname[5])
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

int get_element(char atom_name[5],char element[3])
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

  /* Return whether atom is a metal */
  return(is_metal);
}
/***********************************************************************

create_atom_record  -  Create a new atom record

***********************************************************************/

struct atom *create_atom_record(struct atom **fst_atom_ptr,
                                struct atom **lst_atom_ptr,
                                struct residue *residue_ptr,
                                char *atom_name,int atom_number,
                                float bvalue,float occupancy)
{
  char element[3];

  int i;

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
      printf("*** ERROR. Unable to allocate memory for");
      printf(" struct atom\n");
      printf("***        Program terminated with error.\n");
      exit (1);
    }

  /* Store the details for this entry */
  atom_ptr->atom_number = atom_number;
  strcpy(atom_ptr->atom_name,atom_name);
/* v.2.0.2 --> */
//  strcpy(atom_ptr->element,element);
/* <-- v.2.0.2 */
  atom_ptr->bvalue = bvalue;
  atom_ptr->occupancy = occupancy;
  atom_ptr->residue_ptr = residue_ptr;

  /* Get the element corresponding to this atom name */
  get_element(atom_name,element);
  strcpy(atom_ptr->element,element);

  /* Initialise remaining fields */
  for (i = 0; i < 3; i++)
    atom_ptr->coords[i] = atom_ptr->original_coords[i] = 0.0;
  for (i = 0; i < 2; i++)
    atom_ptr->label_coords[i] = 0.0;
  atom_ptr->interacting = FALSE;
  atom_ptr->fake_name[0] = '\0';
/* v.1.3 --> */
  atom_ptr->anchor_weight = 1.0;
/* <-- v.1.3 */

  /* If the first atom of this residue, then store pointer */
  if (residue_ptr->first_atom_ptr == NULL)
    residue_ptr->first_atom_ptr = atom_ptr;

  /* Increment count of atoms in this residue */
  residue_ptr->natoms++;

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

get_tokens - Get the space- or tab-delimited token from the given
             input line

**********************************************************************/  

int get_tokens(char *line,char token[MAX_TOKEN][TOKEN_LEN],
               char token_sep)
{
/* v.1.3.6 --> */
  // char ch;
  char ch, tmp[TOKEN_LEN];
/* <-- v.1.3.6 */

  int intok, ipos, itoken, ntoken;
  int done;

  /* Initialise variables */
  done = FALSE;
  ntoken = 0;
  intok = 0;
  ipos = 0;
  for (itoken = 0; itoken < MAX_TOKEN; itoken++)
    token[itoken][0] = '\0';

  /* Loop until end of line reached */
  while (done == FALSE && ipos < LINELEN)
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
      if (ch == token_sep)
        {
          /* If have a token, increment token count */
          if (token[ntoken][0] != '\0')
            {
              /* Terminate current token */
              token[ntoken][intok] = '\0';

              /* Check for bounding quotes */
              if (token[ntoken][0] == '"' &&
                  token[ntoken][intok - 1] == '"')
                {
                  /* Strip out the quotes */
/* v.1.3.6 --> */
                  // strcpy(token[ntoken],token[ntoken]+1);
                  strcpy(tmp,token[ntoken]+1);
                  strcpy(token[ntoken],tmp);
/* <-- v.1.3.6 */
                  token[ntoken][intok - 2] = '\0';
                }

              /* Increment token count */
              ntoken++;

              /* Reinitialise position within current token */
              intok = 0;
            }
        }

      /* Otherwise, add character to token-string provided not too long */
      else
        {
          if (ntoken < MAX_TOKEN && intok < TOKEN_LEN)
            {
              token[ntoken][intok] = ch;
              intok++;
            }
        }

      /* Go to the next character in the line */
      ipos++;
    }

  /* Return the number of tokens picked up */
  return(ntoken);
} 
/***********************************************************************

get_drw_atom  -  Extract atom details from the given .drw line

***********************************************************************/

struct atom *get_drw_atom(char *input_line,struct atom **fst_atom_ptr,
                          struct atom **lst_atom_ptr,int natoms,
                          struct residue *residue_ptr)
{
  char atom_type[5];
  char label[LINELEN], token[MAX_TOKEN][TOKEN_LEN];

  int item, itoken, len, ntokens;

  float atom_x, atom_y, label_x, label_y, value;

  struct atom *atom_ptr, *first_atom_ptr, *last_atom_ptr;

  /* Initialise variables */
  atom_x = 0.0;
  atom_y = 0.0;
  label[0] = '\0';
  label_x = 0.0;
  label_y = 0.0;

  /* Initialise pointers */
  atom_ptr = NULL;
  first_atom_ptr = *fst_atom_ptr;
  last_atom_ptr = *lst_atom_ptr;

  /* Get the line length */
  len = strlen(input_line);

  /* Extract the atom name */
  strncpy(atom_type,input_line,4);
  atom_type[4] = '\0';

  /* Create a new atom record */
  atom_ptr
    = create_atom_record(&first_atom_ptr,&last_atom_ptr,residue_ptr,
                         atom_type,natoms,10.0,1.0);

  /* Read in all the tokens from the remainder of this line */
  ntokens = get_tokens(input_line + 4,token,' ');
  item = 0;

  /* Loop over the tokens to store the data */
  for (itoken = 0; itoken < ntokens; itoken++)
    {
      /* If this is the hash, reset token number */
      if (!strcmp(token[itoken],"#"))
        item = 5;

      /* If not the hash, convert string to floating point number */
      if (item < 4 || item > 5)
        value = atof(token[itoken]);

      /* Interpret the token according to its type */
      switch(item)
        {
          /* Store the atom's LIGPLOT x coordinate */
        case 0:
          atom_ptr->coords[X] = value;
          break;

          /* Store the atom's LIGPLOT y coordinate */
        case 1:
          atom_ptr->coords[Y] = value;
          break;

          /* Store the label x coordinate */
        case 2:
          atom_ptr->label_coords[X] = value;
          break;

          /* Store the label y coordinate */
        case 3:
          atom_ptr->label_coords[Y] = value;
          break;

          /* Store the atom's label */
        case 4:
          strcpy(label,token[itoken]);
          break;

          /* Store the atom's original x coordinate */
        case 6:
          atom_ptr->original_coords[X] = value;
          break;

          /* Store the atom's original y coordinate */
        case 7:
          atom_ptr->original_coords[Y] = value;
          break;

          /* Store the atom's original z coordinate */
        case 8:
          atom_ptr->original_coords[Z] = value;
          break;

          /* Default */
        default:
          break;
        }

      /* Increment the item number */
      item++;
    }

  /* Return current pointers */
  *fst_atom_ptr = first_atom_ptr;
  *lst_atom_ptr = last_atom_ptr;

  /* Return the create atom record */
  return(atom_ptr);
}
/***********************************************************************

create_bond_record  -  Create a new bond record

***********************************************************************/

struct bond *create_bond_record(struct bond **fst_bond_ptr,
                                struct bond **lst_bond_ptr,int type,
                                float length,int bond_order,
                                struct atom *at1_ptr,
                                struct atom *at2_ptr)
{
  struct atom *atom_ptr, *atom1_ptr, *atom2_ptr;
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
      printf("***        Program terminated with error.\n");
      exit (1);
    }

  /* Save atom pointers */
  atom1_ptr = at1_ptr;
  atom2_ptr = at2_ptr;

  /* Swap atoms if necessary such that first atom is the lower-numbered
     one */
  if (atom2_ptr->atom_number < atom1_ptr->atom_number)
    {
      atom_ptr = atom1_ptr;
      atom1_ptr = atom2_ptr;
      atom2_ptr = atom_ptr;
    }

  /* Mark both atoms as interacting */
  atom1_ptr->interacting = TRUE;
  atom2_ptr->interacting = TRUE;

  /* Store the details for this entry */
  bond_ptr->atom1_ptr = atom1_ptr;
  bond_ptr->atom2_ptr = atom2_ptr;
  bond_ptr->bond_order = bond_order;
  bond_ptr->length = length;
  bond_ptr->type = type;

  /* Initialise pointer to next bond record */
  bond_ptr->next_bond_ptr = NULL;

  /* If this is the first entry stored, then update
     pointer to head of linked list */
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

get_drw_bond  -  Extract bond details from the given .drw line

***********************************************************************/

struct bond *get_drw_bond(char *input_line,struct bond **fst_bond_ptr,
                          struct bond **lst_bond_ptr,
                          struct atom *first_atom_ptr)
{
  char type[2], token[MAX_TOKEN][TOKEN_LEN];

  int atom1, atom2, bond_order, bond_type, len, ntokens;

  float bond_length;

  struct bond *bond_ptr, *first_bond_ptr, *last_bond_ptr;
  struct atom *atom_ptr, *atom1_ptr, *atom2_ptr;

  /* Initialise variables */
  bond_length = 0.0;
  bond_order = SINGLE;

  /* Initialise pointers */
  atom1_ptr = NULL;
  atom2_ptr = NULL;
  bond_ptr = NULL;
  first_bond_ptr = *fst_bond_ptr;
  last_bond_ptr = *lst_bond_ptr;

  /* Get the line length */
  len = strlen(input_line);

  /* Proceed only if line is long enough to contain the required data */
  if (len > 4)
    {
      /* Extract the bond type */
      type[0] = input_line[0];
      type[1] = '\0';
      bond_type = atoi(type);
      bond_order = SINGLE;

      /* Extract all the tokens on the input line */
      ntokens = get_tokens(input_line,token,' ');

      /* Extract bond's two atoms */
      atom1 = atoi(token[1]);
      atom2 = atoi(token[2]);

      /* Get the pointers for these two atoms */
      atom_ptr = first_atom_ptr;

      /* Loop through all the atoms to find the matching ones */
      while (atom_ptr != NULL && (atom1_ptr == NULL || atom2_ptr == NULL))
        {
          /* If have a match for the first atom, store pointer */
          if (atom_ptr->atom_number == atom1)
            atom1_ptr = atom_ptr;

          /* Or for second atom */
          else if (atom_ptr->atom_number == atom2)
            atom2_ptr = atom_ptr;

          /* Get the next atom record */
          atom_ptr = atom_ptr->next_atom_ptr;
        }

      /* If line contains coords bond label, read it in */
      if (ntokens > 3)
        {
          /* Check whether label simply gives the bond order */
          if (bond_type != HBOND && bond_type != INTERNAL)
            {
              if (token[3][0] == 'd')
                bond_order = DOUBLE;
              else if (token[3][0] == 't')
                bond_order = TRIPLE;
              else if (token[3][0] == 'a')
                bond_order = AROMATIC;
            }

          /* Otherwise, label gives the H-bond length */
          else
            {
              /* Extract bond length */
              bond_length = (float) atof(token[3]);
            }
        }

      /* If have a bond between two atoms, store it */
      if (atom1_ptr != NULL && atom2_ptr != NULL)
        {
          /* Create a new bond record */
          bond_ptr = create_bond_record(&first_bond_ptr,&last_bond_ptr,
                                        bond_type,bond_length,bond_order,
                                        atom1_ptr,atom2_ptr);
        }
    }

  /* Return current pointers */
  *fst_bond_ptr = first_bond_ptr;
  *lst_bond_ptr = last_bond_ptr;

  /* Return the bond */
  return(bond_ptr);
}
/***********************************************************************

read_drw_file  -  Read in and the data from the given .drw file

***********************************************************************/

int read_drw_file(struct pdbentry **fst_pdbentry_ptr,char *fname,
                  int *nr,struct parameters params)
{
  char file_name[FILENAME_LEN], input_line[LINELEN + 1];
/* v.2.2.8 --> */
//  char residue_id[20];
  char residue_id[LINELEN];
/* <-- v.2.2.8 */
/* v.2.0 --> */
  char atom_coords[30], last_atom_coords[30], number_string[20];

  char *string_ptr;
/* <-- v.2.0 */

  int atom_number, len, line, ltype, natoms, nbonds, npdb, nresid;
/* v.2.0 --> */
  int molecule;
/* <-- v.2.0 */

  struct atom *atom_ptr, *first_atom_ptr, *last_atom_ptr;
  struct bond *bond_ptr, *first_bond_ptr, *last_bond_ptr;
  struct pdbentry *pdbentry_ptr;
  struct pdbentry *first_pdbentry_ptr, *last_pdbentry_ptr;
  struct residue *residue_ptr, *first_residue_ptr, *last_residue_ptr;

  FILE *file_ptr;

  /* Initialise */
  atom_number = 0;
  line = 0;
/* v.2.0 --> */
  last_atom_coords[0] = '\0';
  molecule = 0;
/* <-- v.2.0 */
  natoms = 0;
  nbonds = 0;
  npdb = 0;
  nresid = 0;
  residue_id[0] = '\0';

  /* Initialise pointers */
  pdbentry_ptr = first_pdbentry_ptr = last_pdbentry_ptr = NULL;
  *fst_pdbentry_ptr = NULL;
/* v.2.0 --> */
  residue_ptr = NULL;
/* <-- v.2.0 */
  first_residue_ptr = last_residue_ptr = NULL;
  first_atom_ptr = last_atom_ptr = NULL;
  first_bond_ptr = last_bond_ptr = NULL;

  /* Form filename of the ensemble.drw file */
  strcpy(file_name,params.directory);
  strcat(file_name,DIR_SEP);
  strcat(file_name,fname);

  /* Open the ensemble.drw input file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    {
      printf("*** ERROR. Unable to open input file: %s\n",file_name);
      return(0);
    }
  printf("  Reading .drw file ...\n\n");
  printf("    FILE: %s\n",file_name);

  /* Loop through the file processing each record */
  while (fgets(input_line,LINELEN,file_ptr)!= NULL)
    {
      /* Increment line-count */
      line++;

      /* Truncate the string to strip off the newline */
      len = string_truncate(input_line,LINELEN);

      /* If a new line-type hash code, determine line type */
      if (len > 1 && input_line[0] == '#')
        {
          /* Get the line type */
          ltype = get_line_type(input_line);

          /* If the start of a new PDB entry, create a new PDB entry
             record */
          if (ltype == NEW_PDB_ENTRY)
            {
              /* Create a pdbentry entry */
              pdbentry_ptr = create_pdbentry_entry(&first_pdbentry_ptr,
                                                   &last_pdbentry_ptr);

              /* Increment count of PDB entries */
              npdb++;

              /* If line contains a PDB identifier, store it */
              if (len > 3)
                {
                  /* Get the PDB id */
                  pdbentry_ptr->pdb_id = atoi(input_line + 3);
                }

              /* Initialise atom number */
              atom_number = 0;

              /* Re-initialise residue, atom and bond pointers */
              first_residue_ptr = last_residue_ptr = NULL;
              first_atom_ptr = last_atom_ptr = NULL;
              first_bond_ptr = last_bond_ptr = NULL;
            }

/* v.2.0 --> */
          /* If new molecule, get its molecule number */
          else if (ltype == MOLECULE_DATA)
            {
              /* Check for molecule identifier */
              if (strlen(input_line) > 2)
                {
                  /* Get the number */
                  strncpy(number_string,input_line+2,1);
                  number_string[1] = '\0';
                  molecule = atoi(number_string);
                }
              else
                molecule = 0;
            }
/* <-- v.2.0 */

          /* If new residue record, get its identifier */
          else if (ltype == RES_DATA)
            {
              /* If line contains a residue identifier, store it */
              if (len > 3)
                {
                  /* Save the residue identifier */
                  strcpy(residue_id,input_line + 3);

/* v.2.0 --> */
                  /* If we have the antibody identifier, strip it
                     off */
                  string_ptr = strchr(residue_id,']');
                  if (string_ptr != NULL)
                    string_ptr[1] = '\0';
/* <-- v.2.0 */
                }
            }
        }

      /* Otherwise, if have valid line, extract its data */
      else if (len > 0 && ltype != NOT_WANTED)
        {
          /* Interpret the current line according to its type */
          switch(ltype)
            {
              /* These line types not needed here */
            case PLOT_TITLE:
            case GLOBAL_PARAMS:
            case COLOURMAP_DATA:
            case COLOUR_DATA:
            case SIZE_DATA:
            case MOLECULE_DATA:
            case DIMPLOT_FLAG:
              break;

              /* New PDB entry - check for PDB code */
            case NEW_PDB_ENTRY:

              /* Get the PDB code */
              if (!strncmp(input_line,"PDB CODE: ",10))
                {
                  if (strlen(input_line+10) <= 4)
                    {
                      strcpy(pdbentry_ptr->pdb_code,input_line+10);
                    }
                  else
                    {
                      strncpy(pdbentry_ptr->pdb_code,input_line+10,4);
                      pdbentry_ptr->pdb_code[4] = '\0';
                    }
                }
              break;

              /* New residue */
            case RES_DATA:

              /* If haven't created a PDB entry yet, do so now */
              if (first_pdbentry_ptr == NULL)
                {
                  /* Create a pdbentry entry */
                  pdbentry_ptr = create_pdbentry_entry(&first_pdbentry_ptr,
                                                       &last_pdbentry_ptr);

                  /* Increment count of PDB entries */
                  npdb++;
                }

              /* Get the residue type, name and number */
              residue_ptr
                = get_drw_residue(input_line,residue_id,&first_residue_ptr,
/* v.2.0 --> */
//                                  &last_residue_ptr,pdbentry_ptr);
                                  &last_residue_ptr,pdbentry_ptr,
                                  molecule,params);
/* <-- v.2.0 */
              if (residue_ptr != NULL)
                {
                  /* If this is this PDB entry's first residue, store
                     its pointer */
                  if (pdbentry_ptr->nresid == 0)
                    pdbentry_ptr->first_residue_ptr = residue_ptr;

                  /* Increment residue counts */
                  pdbentry_ptr->nresid++;
                  nresid++;
                }

/* v.2.0 --> */
              /* Reinitialize formatted coords */
              last_atom_coords[0] = '\0';
/* <-- v.2.0 */

              break;

              /* New atom */
            case ATOM_DATA:

              /* Get the atom name and coords */
              atom_ptr
                = get_drw_atom(input_line,&first_atom_ptr,&last_atom_ptr,
                               atom_number,residue_ptr);
              if (atom_ptr != NULL)
                {
                  atom_number++;
                  natoms++;
/* v.2.0 --> */
                  /* Get this atoms' formatted coords */
                  if (params.dimplot == TRUE)
                    {
                      /* Format the coords */
                      sprintf(atom_coords,"%8.3f %8.3f",
                              atom_ptr->coords[X],atom_ptr->coords[Y]);

                      /* If the same coords as before, mark this
                         residue as a hydrophobic group */
                      if (!strcmp(atom_coords,last_atom_coords))
                        {
                          /* Mark the residue */
                          residue_ptr->type = HYDROPHOBIC;

                          /* Store the formatted coords */
                          strcpy(last_atom_coords,atom_coords);
                        }
                    }
/* <-- v.2.0 */
                }

              break;

              /* New bond */
            case BOND_DATA:

              /* Extract the bond details */
              bond_ptr
                = get_drw_bond(input_line,&first_bond_ptr,&last_bond_ptr,
                               first_atom_ptr);
              if (bond_ptr != NULL)
                {
                  /* If this is this PDB entry's first residue, store
                     it pointer */
                  if (pdbentry_ptr->nbonds == 0)
                    pdbentry_ptr->first_bond_ptr = bond_ptr;

                  /* Increment residue counts */
                  pdbentry_ptr->nbonds++;
                  nbonds++;
                }
              break;

              /* Default */
            default:
              break;
            }
        }
    }

  /* Close the file */
  fclose(file_ptr);

  /* Return pointer to first pdbentry record*/
  *fst_pdbentry_ptr = first_pdbentry_ptr;
  printf("Number of PDB entries in the ensemble: %8d\n",npdb);
  printf("Number of residues stored:             %8d\n",nresid);
  printf("Number of atoms stored:                %8d\n",natoms);
  printf("Number of bonds stored:                %8d\n",nbonds);

  /* Return the number of residues */
  *nr = nresid;

  /* Return number of definitions read in */
  return(npdb);
}
/***********************************************************************

perform_sort  -  Bubble sort the residue by residue number

***********************************************************************/

void perform_sort(struct residue **fst_residue_ptr)
{
  int swap, swapped;

  struct residue *residue_ptr, *first_residue_ptr;
  struct residue *next_residue_ptr, *prev_residue_ptr, *swap_residue_ptr;

  /* Get pointer to the first residueerence record */
  first_residue_ptr = *fst_residue_ptr;
  swapped = TRUE;

  /* Loop until residues in required order */
  while (swapped == TRUE)
    {
      /* Initialise swapped flag */
      swapped = FALSE;
          
      /* Get pointer to first residue record */
      residue_ptr = first_residue_ptr;
      prev_residue_ptr = NULL;

      /* Loop over all the residues, adjusting their order where
         necessary */
      while (residue_ptr != NULL)
        {
          /* Get pointer to the next record */
          next_residue_ptr = residue_ptr->next_residue_ptr;

          /* Determine whether current two records are in the right
             order according to the special marker */
          swap = FALSE;
          if (next_residue_ptr != NULL)
            {
              /* Check if chains are in right order (with blank chains
                 going to the back) */
              if ((residue_ptr->chain == ' ' &&
                   next_residue_ptr->chain != ' ') ||
                  (next_residue_ptr->chain != ' ' &&
                   residue_ptr->chain > next_residue_ptr->chain))
                swap = TRUE;

              /* If chains the same, check if residue numbers in
                 right order */
              else if (residue_ptr->chain == next_residue_ptr->chain &&
                       strcmp(residue_ptr->res_num,
                              next_residue_ptr->res_num) > 0)
                swap = TRUE;
            }

          /* If need to swap, perform the swap */
          if (swap == TRUE)
            {
              /* If this is the first residue record, adjust first pointer */
              if (residue_ptr == first_residue_ptr)
                first_residue_ptr = next_residue_ptr;

              /* Otherwise, get the previous residue pointer to point
                 to next one */
              else
                prev_residue_ptr->next_residue_ptr = next_residue_ptr;

              /* Swap pointers within the records */
              residue_ptr->next_residue_ptr
                = next_residue_ptr->next_residue_ptr;
              next_residue_ptr->next_residue_ptr = residue_ptr;

              /* Swap the pointers themselves */
              swap_residue_ptr = residue_ptr;
              residue_ptr = next_residue_ptr;
              next_residue_ptr = swap_residue_ptr;

              /* Set flag that swap made */
              swapped = TRUE;
            }

          /* Save current pointer */
          prev_residue_ptr = residue_ptr;

          /* Get the next residue record */
          residue_ptr = next_residue_ptr;
        }

      /* Save the first residue pointer */
      *fst_residue_ptr = first_residue_ptr;
    }
}
/***********************************************************************
  
sort_residues -  Sort the residues in each PDB entry in the ensemble

***********************************************************************/

void sort_residues(struct pdbentry *first_pdbentry_ptr)
{
  struct pdbentry *pdbentry_ptr;

  /* Get the first PDB entry */
  pdbentry_ptr = first_pdbentry_ptr;

  /* Loop through the PDB entries */
  while (pdbentry_ptr != NULL)
    {
      /* Sort the residues into residue number order */
      perform_sort(&(pdbentry_ptr->first_residue_ptr));

      /* Get the next entry */
      pdbentry_ptr = pdbentry_ptr->next_pdbentry_ptr;
    }
}
/***********************************************************************

init_filenames  -  Form the names of the input PDB and HBPLUS files

***********************************************************************/

/* v.1.3 --> */
//void init_filenames(char *hhb_name,char *nnb_name,char *pdb_name,
//                    struct parameters params)
void init_filenames(char *hhb_name,char *nnb_name,char *ligall_pdb,
                    char *ligplus_pdb,struct parameters params)
/* <-- v.1.3 */
{
  /* Form the names of the input PDB and HBPLUS files */
  if (params.directory != &null)
    {
      /* Form the file names in the working directory */
      strcpy(hhb_name,params.directory);
      strcat(hhb_name,DIR_SEP);
      strcat(hhb_name,"ligplus.hhb");
      strcpy(nnb_name,params.directory);
      strcat(nnb_name,DIR_SEP);
      strcat(nnb_name,"ligplus.nnb");
/* v.1.3 --> */
//      strcpy(pdb_name,params.directory);
//      strcat(pdb_name,DIR_SEP);
/* v.1.2 --> */
/*      strcat(pdb_name,"ligplus.pdb"); */
//      if (params.alnfile == &null)
//        strcat(pdb_name,"ligplus.pdb");
//      else
//        strcat(pdb_name,"ligall.pdb");
      strcpy(ligplus_pdb,params.directory);
      strcat(ligplus_pdb,DIR_SEP);
      strcat(ligplus_pdb,"ligplus.pdb");
      strcpy(ligall_pdb,params.directory);
      strcat(ligall_pdb,DIR_SEP);
      strcat(ligall_pdb,"ligall.pdb");
/* <-- v.1.3 */
/* <-- v.1.2 */
    }
  else
    {
      strcpy(hhb_name,"ligplus.hhb");
      strcpy(nnb_name,"ligplus.nnb");
/* v.1.2 --> */
/*      strcpy(pdb_name,"ligplus.pdb"); */
/* v.1.3 --> */
//      if (params.alnfile == &null)
//        strcpy(pdb_name,"ligplus.pdb");
//      else
//        strcpy(pdb_name,"ligall.pdb");
      strcpy(ligall_pdb,"ligall.pdb");
      strcpy(ligplus_pdb,"ligplus.pdb");
/* <-- v.1.3 */
/* <-- v.1.2 */
    }
}
/* v.1.1 --> */
/***********************************************************************

create_hhbnnb_record  -  Create and initialise a new hhbnnb record

***********************************************************************/

struct hhbnnb *create_hhbnnb_record(struct hhbnnb **fst_hhbnnb_ptr,
                                    struct hhbnnb **lst_hhbnnb_ptr,
                                    char *data)
{
  struct hhbnnb *hhbnnb_ptr;
  struct hhbnnb *last_hhbnnb_ptr, *first_hhbnnb_ptr;

  /* Initialise hhbnnb pointers */
  hhbnnb_ptr = NULL;
  last_hhbnnb_ptr = *lst_hhbnnb_ptr;
  first_hhbnnb_ptr = *fst_hhbnnb_ptr;

  /* Allocate memory for structure to hold sequence info */
  hhbnnb_ptr = (struct hhbnnb *) malloc(sizeof(struct hhbnnb));
  if (hhbnnb_ptr == NULL)
    {
      printf("*** Can't allocate memory for struct");
      printf(" hhbnnb\n");
      exit (1);
    }

  /* Store the details for this entry */
  store_string(&(hhbnnb_ptr->data),data);

  /* Set pointer to next record */
  hhbnnb_ptr->next_hhbnnb_ptr = NULL;

  /* If this is the first entry stored, then update pointer to head
     of linked list */
  if (first_hhbnnb_ptr == NULL)
    first_hhbnnb_ptr = hhbnnb_ptr;

  /* Otherwise, make previous entry point to the current
     one */
  else
    last_hhbnnb_ptr->next_hhbnnb_ptr = hhbnnb_ptr;
                      
  /* Save the current hhbnnb's pointer */
  last_hhbnnb_ptr = hhbnnb_ptr;

  /* Return current pointers */
  *fst_hhbnnb_ptr = first_hhbnnb_ptr;
  *lst_hhbnnb_ptr = last_hhbnnb_ptr;
  return(hhbnnb_ptr);
}
/* <-- v.1.1 */
/***********************************************************************

get_atom_details  -  Extract the atom details from the given line

***********************************************************************/

void get_atom_details(char *input_line,int *athet,int *atom_number,
                      char *atom_name,char *chain,char *res_name,
                      char *res_num,float *xval,float *yval,float *zval,
                      float *bval,float *occup,char *element)
{
  char number_string[20];

  int len;
  int is_metal;

  float bvalue, occupancy, x, y, z;

  /* Initialise variables */
  is_metal = FALSE;

  /* Get the residue details */
  strncpy(res_name,input_line+17,3);
  res_name[3] = '\0';
  strncpy(res_num,input_line+22,5);
  res_num[5] = '\0';
  *chain = input_line[21];

  /* Get the atom type */
  if (!strncmp(input_line,"ATOM",4))
    *athet = ATOM;
  else
    *athet = HETATM;

  /* Get the atom-number */
  strncpy(number_string,input_line+6,5);
  number_string[5] = '\0';
  *atom_number = atoi(number_string);

  /* Get the atom name */
  strncpy(atom_name,input_line+12,4);
  atom_name[4] = '\0';

  /* Retrieve the atomic coordinates, occupancy and B-value */
  sscanf(input_line+30," %f %f %f %f %f",&x,&y,&z,&occupancy,
         &bvalue);

  /* Store the retrieved values */
  *xval = x;
  *yval = y;
  *zval = z;
  *bval = bvalue;
  *occup = occupancy;

  /* Get the length of the input line */
  len = strlen(input_line);

  /* Get the element identifier, if there is one */
  if (len > 77)
    {
      strncpy(element,input_line+76,2);
      element[2] = '\0';
    }
  else
    strcpy(element,"  ");

  /* Check whether this is a valid element */
  if ((element[0] >= '0' && element[0] <= '9') ||
      (element[1] >= '0' && element[1] <= '9'))
    strcpy(element,"  ");

  /* If don't have the element, have to extract it from the atom name */
  if (!strcmp(element,"  "))
    get_element(atom_name,element);
}
/***********************************************************************

check_for_duplicate  -  Check for an earlier copy of the current atom

***********************************************************************/

int check_for_duplicate(struct atom *first_atom_ptr,char *atom_name,
                        char *res_name,char *res_num,char chain)
{
  int have_atom;

  struct atom *atom_ptr;

  /* Initialise variables */
  have_atom = FALSE;

  /* Get the first atom record */
  atom_ptr = first_atom_ptr;
      
  /* Loop through all the atoms */
  while (atom_ptr != NULL && have_atom == FALSE)
    {
      /* Check for duplication */
      if (!strcmp(atom_name,atom_ptr->atom_name) &&
          !strcmp(res_name,atom_ptr->residue_ptr->res_name) &&
          !strcmp(res_num,atom_ptr->residue_ptr->res_num) &&
          chain == atom_ptr->residue_ptr->chain)
        have_atom = TRUE;
                  
      /* Get next atom record */
      atom_ptr = atom_ptr->next_atom_ptr;
    }

  /* Return flag */
  return(have_atom);
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

calc_distance  -  Calculate the distance between the two points

***********************************************************************/

float calc_distance(double *coords1,double *coords2)
{
  float distance, dist2;

  /* Initialise variables */
  distance = 0.0;

  /* Calculate the squared distance */
  dist2 = (coords1[X] - coords2[X]) * (coords1[X] - coords2[X])
    + (coords1[Y] - coords2[Y]) * (coords1[Y] - coords2[Y])
    + (coords1[Z] - coords2[Z]) * (coords1[Z] - coords2[Z]);

  /* Take root */
  distance = sqrt(dist2);

  /* Return the distance */
  return(distance);
}
/***********************************************************************

unmark_residues  -  Go through the residues marked as ligands and unmark
                    them

***********************************************************************/

void unmark_residues(struct residue *first_residue_ptr)
{
  int done;

  struct residue *residue_ptr;


  /* Initialise variables */
  done = FALSE;
  
  /* Get pointer to the first residue */
  residue_ptr = first_residue_ptr;

  /* Loop over all the residues to change their type */
  while (residue_ptr != NULL && done == FALSE)
    {
      /* If residue marked as a ligand, change its type */
      if (residue_ptr->type == LIGAND)
        {
          /* If not a water, set the type */
          if (strcmp(residue_ptr->res_name,"HOH"))
            residue_ptr->type = HGROUP;
        }

      /* Otherwise, we're done */
      else
        done = TRUE;

      /* Get pointer to the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
    }
}
/***********************************************************************

read_pdb_file  -  Read in the ligand and protein data from the
                  ligplus.pdb file

***********************************************************************/

/* v.2.2.1 --> */
// struct pdbentry *read_pdb_file(char *file_name,
struct pdbentry *read_pdb_file(char *file_name,int file_type,
                               struct residue **fst_lig_res_ptr,
/* <-- v.2.2.1 */
                               struct parameters params)
{
  char input_line[LINELEN + 1];
  char atom_name[5], element[3], number_string[20];
  char res_name[4], res_num[6], last_res_name[4], last_res_num[6];
  char chain, label[2], last_chain, occup[4], residue_id[20];

  int athet, atom_number, iconect, ipos, len;
  int natoms, nbonds, nresid;
/* v.2.0 --> */
  int domain;
/* <-- v.2.0 */
  int bond_order, bond_type, res_type;
  int done, keep_reading, have_atom, hydrogen, wanted;
  int inligand;
/* v.2.2.1 --> */
  int ligand_start;
/* <-- v.2.2.1 */

  float bond_length, bvalue, label_x, label_y, occupancy, x, y, z;

  struct atom *atom_ptr, *first_atom_ptr, *last_atom_ptr;
  struct atom *atom1_ptr, *atom2_ptr;
  struct bond *bond_ptr, *first_bond_ptr, *last_bond_ptr;
  struct pdbentry *pdbentry_ptr;
  struct pdbentry *first_pdbentry_ptr, *last_pdbentry_ptr;
  struct residue *residue_ptr, *first_residue_ptr, *last_residue_ptr;
/* v.2.2.1 --> */
  struct residue *first_ligand_residue_ptr;
/* <-- v.2.2.1 */
/* v.1.1 --> */
  struct hhbnnb *hhbnnb_ptr, *first_hhbnnb_ptr, *last_hhbnnb_ptr;
/* <-- v.1.1 */
  
  FILE *file_ptr;

  /* Initialise variables */
  atom_number = 0;
/* v.2.0 --> */
  domain = 0;
/* <-- v.2.0 */
  inligand = FALSE;
  keep_reading = TRUE;
  label[0] = '\0';
  label_x = label_y = 0.0;
  last_res_name[0] = '\0';
  last_res_num[0] = '\0';
  last_chain = '\0';
/* v.2.2.1 --> */
  ligand_start = FALSE;
/* <-- v.2.2.1 */
  natoms = 0;
  nbonds = 0;
  nresid = 0;

  /* Initialise pointers */
  atom_ptr = first_atom_ptr = last_atom_ptr = NULL;
  bond_ptr = first_bond_ptr = last_bond_ptr = NULL;
  pdbentry_ptr = first_pdbentry_ptr = last_pdbentry_ptr = NULL;
  residue_ptr = first_residue_ptr = last_residue_ptr = NULL;
/* v.1.1 --> */
  hhbnnb_ptr = first_hhbnnb_ptr = last_hhbnnb_ptr = NULL;
/* <-- v.1.1 */
/* v.2.2.1 --> */
  first_ligand_residue_ptr = *fst_lig_res_ptr;
/* <-- v.2.2.1 */

  /* Create the PDB entry record */
  pdbentry_ptr
    = create_pdbentry_entry(&first_pdbentry_ptr,&last_pdbentry_ptr);

  /* Mark this entry as current */
  strcpy(pdbentry_ptr->pdb_code,"Curr");

  /* Open PDB file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    {
      printf("\n*** Warning. Unable to open PDB file [%s]\n",file_name);
      return(pdbentry_ptr);
    }

  /* Loop through the PDB file, to read in all the atoms */
  while (fgets(input_line,LINELEN,file_ptr)!= NULL)
    {
      /* Truncate the string to strip off the newline */
      len = string_truncate(input_line,LINELEN);

      /* Initialise flags */
      wanted = FALSE;
/* v.2.2.1 --> */
      ligand_start = FALSE;
/* <-- v.2.2.1 */

      /* Check for the end of an NMR structure model */
      if (!strncmp(input_line,"ENDMDL",6))
        keep_reading = FALSE;

      /* Check for the PDB code */
      else if (!strncmp(input_line,"REMARK PDB CODE: ",17))
        {
          input_line[21] = '\0';
          strcpy(pdbentry_ptr->pdb_code,input_line+17);
        }

/* v.2.0 --> */
      /* Check for the domain identifier */
      else if (!strncmp(input_line,"REMARK DOMAIN ",14))
        {
          strcpy(number_string,input_line+14);
          domain = atoi(number_string);
        }
/* <-- v.2.0 */

/* v.1.1 --> */
      /* If this is an HHB/NNB records, then save it for writing out
         later */
      else if (!strncmp(input_line,"HHB   ",6) ||
               !strncmp(input_line,"NNB   ",6) ||
               !strncmp(input_line,"-HB   ",6) ||
               !strncmp(input_line,"-NB   ",6))
        {
          // Create a hhbnnb record to store data */
          hhbnnb_ptr
            = create_hhbnnb_record(&first_hhbnnb_ptr,
                                    &last_hhbnnb_ptr,input_line);
        }
/* <-- v.1.1 */

      /* If this is an ATOM or HETATM record, process it */
      else if (keep_reading == TRUE &&
               (!strncmp(input_line,"ATOM",4) ||
                !strncmp(input_line,"HETATM",6)))
        {
          /* Retrieve the information from the line */
          get_atom_details(input_line,&athet,&atom_number,
                           atom_name,&chain,res_name,res_num,
                           &x,&y,&z,&bvalue,&occupancy,element);

          /* Determine whether this is a wanted ATOM or HETATM record */
          wanted = TRUE;

          /* Remove any hydrogens */
          if (!strcmp(element," H"))
            hydrogen = TRUE;
          else
            hydrogen = FALSE;
          if (hydrogen == TRUE)
            wanted = FALSE;

          /* If this is a water molecule, and we don't want waters,
             then flag as unwanted */
          if (params.include_waters == FALSE && !strcmp(res_name,"HOH"))
            wanted = FALSE;

          /* If this is an alternate-occupancy atom, check whether
             we already have it */
          if (len > 58)
            {
              strncpy(occup,input_line+56,3);
              occup[3] = '\0';
            }
          else
            occup[0] = '\0';
          if (wanted == TRUE &&
              (input_line[16] != ' ' || strncmp(occup,"1.0",3)))
            {
              /* Check for earlier copy of this atom */
              have_atom
                = check_for_duplicate(first_atom_ptr,atom_name,
                                      res_name,res_num,chain);

              /* If have this atom already, then don't want this
                 second copy */
              if (have_atom == TRUE)
                wanted = FALSE;
            }

          /* If atom still wanted, then prepare to store */
          if (wanted == TRUE)
            {
              /* Check whether this is a new residue */
              if (chain != last_chain ||
                  strncmp(res_num,last_res_num,5) ||
                  strncmp(res_name,last_res_name,3))
                {
                  /* Determine whether the residue belongs to the
                     ligand */
                  res_type = HGROUP;
                  if (inligand == FALSE)
                    {
                      /* Check for a match with the first ligand residue */
/* v.2.0 --> */
//                      if (chain == params.chain &&
                      if (params.dimplot == FALSE &&
                          chain == params.chain[0] &&
/* <-- v.2.0 */
                          !strcmp(res_num,params.res_num1) &&
                          (!strcmp(res_name,params.res_name1) ||
/* v.2.2.1 --> */
//                           params.res_name1[0] == '\0'))
                           params.res_name1[0] == '\0') &&
                          file_type == LIGPLUS_FILE)
/* <-- v.2.2.1 */
                        {
                          /* Have match, so set flag that this is
                             the ligand */
                          inligand = TRUE;
                          res_type = LIGAND;

/* v.2.2.1 --> */
                          /* If we already have "started" the ligand,
                             it is likely that a protein residue has
                             the same number and chain id as the first
                             ligand residue */
                          if (first_ligand_residue_ptr != NULL)
                            {
                              /* Go through the residues marked as ligands
                                 and unmark them */
                              unmark_residues(first_ligand_residue_ptr);
                            }

                          /* Set the flag that is is a ligand start */
                          first_ligand_residue_ptr = NULL;
                          ligand_start = TRUE;
/* <-- v.2.2.1 */
                        }
                    }
/* v.2.0 --> */
//                  else
                  else if (params.dimplot == FALSE)
/* <-- v.2.0 */
                    res_type = LIGAND;

                  /* If water, then mark as such */
                  if (!strcmp(res_name,"HOH"))
                    res_type = WATER;

                  /* Check whether this is the last residue of the ligand
                     range */
/* v.2.0 --> */
//                  if (chain == params.chain &&
                  if (chain == params.chain[0] &&
/* <-- v.2.0 */
                      !strcmp(res_num,params.res_num2) &&
                      (!strcmp(res_name,params.res_name2) ||
                       params.res_name2[0] == '\0'))
                    {
                      /* Have match, so unset ligand flag */
                      inligand = FALSE;
                    }

                  /* Create the residue identifier */
                  sprintf(residue_id,"[%s %s %c]",res_name,res_num,chain);

                  /* Create a new residue record */
                  residue_ptr
                    = create_residue_record(&first_residue_ptr,
                                            &last_residue_ptr,res_name,
                                            res_num,chain,res_type,
                                            label,label_x,label_y,
                                            residue_id,pdbentry_ptr);

                  /* Increment count of residues read in */
                  nresid++;
/* v.2.2.1 --> */
                  /* If this is the first residue of the ligand, save
                     its pointer */
                  if (ligand_start == TRUE)
                    first_ligand_residue_ptr = residue_ptr;
/* <-- v.2.2.1 */

/* v.2.0 --> */
                  /* Save the chain type for DIMPLOT fits */
                  if (!strcmp(res_name,"HOH"))
                    residue_ptr->chain_type = INTERFACE_WATER;
                  else if (domain == 1)
                    residue_ptr->chain_type = CHAIN1;
                  else if (domain == 2)
                    residue_ptr->chain_type = CHAIN2;
/* <-- v.2.0 */

/* v.2.0.1 --> */
                  /* Save the current residue name, number and chain */
                  strncpy(last_res_name,res_name,3);
                  last_res_name[3] = '\0';
                  strncpy(last_res_num,res_num,5);
                  last_res_num[5] = '\0';
                  last_chain = chain;
/* <-- v.2.0.1 */
                }

              /* Save current atom's details */
              atom_ptr
                = create_atom_record(&first_atom_ptr,&last_atom_ptr,
/* v.2.0.3 --> */
//                                     residue_ptr,atom_name,natoms,
                                     residue_ptr,atom_name,atom_number,
/* <-- v.2.0.3 */
                                     bvalue,occupancy);

              /* Save the x, y, z coordinates */
              atom_ptr->original_coords[X] = atom_ptr->coords[X] = x;
              atom_ptr->original_coords[Y] = atom_ptr->coords[Y] = y;
              atom_ptr->original_coords[Z] = atom_ptr->coords[Z] = z;

              /* Increment count of atoms read in */
              natoms++;
            }

/* v.2.0.1 --> */
//          /* Save the current residue name, number and chain */
//          strncpy(last_res_name,res_name,3);
//          last_res_name[3] = '\0';
//          strncpy(last_res_num,res_num,5);
//          last_res_num[5] = '\0';
//          last_chain = chain;
/* <-- v.2.0.1 */
        }

      /* If this is a CONECT record, then create a bond record if wanted */
      else if (!strncmp(input_line,"CONECT",6))
        {
          /* Initialise variables */
          atom1_ptr = atom2_ptr = NULL;
          done = FALSE;
          ipos = 6;

          /* Loop through all the possible atom-number positions */
          for (iconect = 0; iconect < MAXCONECT && done == FALSE;
               iconect++)
            {
              /* If this position is blank, then have reached end of
                 atom-numbers on this line */
              if (!strncmp(input_line+ipos,"     ",5))
                done = TRUE;

              /* Otherwise, store this atom number */
              else
                {
                  /* Get the atom-number at this position */
                  strncpy(number_string,input_line+ipos,5);
                  number_string[5] = '\0';

                  /* If first position, store first atom number */
                  if (iconect == 0)
                    {
                      /* Get the atom number */
                      atom_number = atoi(number_string);

                      /* Locate this atom */
                      atom1_ptr = find_atom(first_atom_ptr,atom_number);
                    }

                  /* Otherwise, create a new conect record */
                  else if (atom1_ptr != NULL)
                    {
                      /* Get the connected atom number */
                      atom_number = atoi(number_string);

                      /* Locate this atom */
                      atom2_ptr = find_atom(first_atom_ptr,atom_number);

                      /* If valid, then create a bond record */
                      if (atom2_ptr != NULL)
                        {
                          /* Define bond type, bond order and bond length */
                          bond_type = COVALENT;
                          bond_order = SINGLE;
                          bond_length
                            = calc_distance(atom1_ptr->original_coords,
                                            atom2_ptr->original_coords);

                          /* Create a bond record */
                          bond_ptr
                            = create_bond_record(&first_bond_ptr,
                                                 &last_bond_ptr,bond_type,
                                                 bond_length,bond_order,
                                                 atom1_ptr,atom2_ptr);

                          /* Increment count */
                          nbonds++;
                        }
                    }
                }

              /* Increment position in line */
              ipos = ipos + 5;
            }
        }
    }
  
  /* Close PDB file */
  fclose(file_ptr);

  /* Show stats */
  printf("Number of residues from PDB file:      %8d\n",nresid);

  /* Save the first residue and bond for this PDB entry */
/* v.1.1 --> */
  pdbentry_ptr->first_hhbnnb_ptr = first_hhbnnb_ptr;
/* <-- v.1.1 */
  pdbentry_ptr->first_residue_ptr = first_residue_ptr;
  pdbentry_ptr->nresid = nresid;
  pdbentry_ptr->first_bond_ptr = first_bond_ptr;
  pdbentry_ptr->nbonds = nbonds;

/* v.2.2.1 --> */
  /* Return the first ligand residue */
  *fst_lig_res_ptr = first_ligand_residue_ptr;
/* <-- v.2.2.1 */

  /* Return the PDB entry */
  return(pdbentry_ptr);
}
/* v.1.3 --> */
/***********************************************************************

find_residue  -  Find the given residue entry

***********************************************************************/

struct residue *find_residue(struct residue *first_residue_ptr,
                             char *res_name,char *res_num,char chain)
{
  struct residue *residue_ptr, *found_residue_ptr;

  /* Initialise variables */
  found_residue_ptr = NULL;

  /* Get pointer to the first residue */
  residue_ptr = first_residue_ptr;

  /* Loop over all the residues to find the one of interest */
  while (residue_ptr != NULL && found_residue_ptr == NULL)
    {
      /* If this is the right residue, store pointer */
      if (chain == residue_ptr->chain &&
          !strcmp(res_num,residue_ptr->res_num) &&
          !strcmp(res_name,residue_ptr->res_name))
        found_residue_ptr = residue_ptr;

      /* Get pointer to the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
    }

  /* Return the residueing residue pointer */
  return(found_residue_ptr);
}
/***********************************************************************

find_atom_name  -  Find given atom in given residue

***********************************************************************/

struct atom *find_atom_name(struct residue *residue_ptr,char *atom_name)
{
  int iatom, natoms;
  int found;

  struct atom *atom_ptr, *match_atom_ptr;

  /* Initialise variables */
  found = FALSE;
  match_atom_ptr = NULL;

  /* Get the first of the given residue's atoms */
  atom_ptr = residue_ptr->first_atom_ptr;
  iatom = 0;
  natoms = residue_ptr->natoms;
  
  /* Loop through all the residue's atoms to find the given atom */
  while (atom_ptr != NULL && iatom < natoms && found == FALSE)
    {
      /* If this is the right atom, store its pointer */
      if (!strcmp(atom_ptr->atom_name,atom_name))
        {
          /* Store pointer and set flag */
          match_atom_ptr = atom_ptr;
          found = TRUE;
        }

      /* Get pointer to next atom in linked-list */
      atom_ptr = atom_ptr->next_atom_ptr;
      iatom++;
    }

  /* Return matched atom */
  return(match_atom_ptr);
}
/***********************************************************************
  
replace_bond_atom  -  Replace the current atom in bonds featuring old
                      atom from .drw file

***********************************************************************/

void replace_bond_atom(struct pdbentry *pdbentry_ptr,
                       struct atom *match_atom_ptr,struct atom *atom_ptr)
{
  struct bond *bond_ptr;

  /* Get pointer to the first bond record */
  bond_ptr = pdbentry_ptr->first_bond_ptr;

  /* Loop over the bonds to find atom to be replaced */
  while (bond_ptr != NULL)
    {
      /* If either atom matches our current atom, then replace */
      if (bond_ptr->atom1_ptr == match_atom_ptr)
        bond_ptr->atom1_ptr = atom_ptr;
      if (bond_ptr->atom2_ptr == match_atom_ptr)
        bond_ptr->atom2_ptr = atom_ptr;
              
      /* Get pointer to next bond record */
      bond_ptr = bond_ptr->next_bond_ptr;
    }
}
/***********************************************************************
  
replace_plot_residues  -  Replace the residues in the full PDB file with
                          the corresponding details read in from the
                          ensemble.drw file

***********************************************************************/

void replace_plot_residues(struct pdbentry *full_pdbentry_ptr,
                           struct pdbentry *first_pdbentry_ptr)
{
  int iatom, natoms, nmatch;

  struct atom *atom_ptr, *match_atom_ptr;
  struct residue *match_residue_ptr, *residue_ptr;

  /* Initialise variables */
  nmatch = 0;

  /* Save the PDB id */
  full_pdbentry_ptr->pdb_id = first_pdbentry_ptr->pdb_id;

  /* Transfer all the bonds read in from the .drw file */
  full_pdbentry_ptr->first_bond_ptr = first_pdbentry_ptr->first_bond_ptr;
  full_pdbentry_ptr->nbonds = first_pdbentry_ptr->nbonds;

  /* Get the full entry's first residue */
  residue_ptr = full_pdbentry_ptr->first_residue_ptr;

  /* Loop over the residues to find corresponding records in entry read
     in from the ensemble.drw file */
  while (residue_ptr != NULL)
    {
      /* Find the corresponding residue from the .drw file */
      match_residue_ptr
        = find_residue(first_pdbentry_ptr->first_residue_ptr,
                       residue_ptr->res_name,residue_ptr->res_num,
                       residue_ptr->chain);

      /* If match found, then update the relevant details */
      if (match_residue_ptr != NULL)
        {
          /* Update current residue record */
          residue_ptr->type = match_residue_ptr->type;
          if (match_residue_ptr->label != &null)
            store_string(&(residue_ptr->label),match_residue_ptr->label);
          residue_ptr->label_x = match_residue_ptr->label_x;
          residue_ptr->label_y = match_residue_ptr->label_y;
          residue_ptr->aa_code = match_residue_ptr->aa_code;
          residue_ptr->interaction = match_residue_ptr->interaction;
/* v.2.0 --> */
          residue_ptr->chain_type = match_residue_ptr->chain_type;
/* <-- v.2.0 */
          if (residue_ptr->residue_id != &null)
            store_string(&(residue_ptr->residue_id),
                         match_residue_ptr->residue_id);

          /* Get the first of this residue's atom records */
          atom_ptr = residue_ptr->first_atom_ptr;
          natoms = residue_ptr->natoms;
          iatom = 0;
      
          /* Loop through all the atoms */
          while (atom_ptr != NULL && iatom < natoms)
            {
              /* Find the corresponding atom from the .drw file */
              match_atom_ptr
                = find_atom_name(match_residue_ptr,atom_ptr->atom_name);

              /* If found, then transfer appropriate details */
              if (match_atom_ptr != NULL)
                {
                  /* Transfer all the data from the .drw atom */
                  atom_ptr->atom_number = match_atom_ptr->atom_number;
                  strcpy(atom_ptr->element,match_atom_ptr->element);
                  strcpy(atom_ptr->fake_name,match_atom_ptr->fake_name);
                  atom_ptr->original_coords[0]
                    = match_atom_ptr->original_coords[0];
                  atom_ptr->original_coords[1]
                    = match_atom_ptr->original_coords[1];
                  atom_ptr->original_coords[2]
                    = match_atom_ptr->original_coords[2];
                  atom_ptr->coords[0] = match_atom_ptr->coords[0];
                  atom_ptr->coords[1] = match_atom_ptr->coords[1];
                  atom_ptr->coords[2] = match_atom_ptr->coords[2];
                  atom_ptr->label_coords[0]
                    = match_atom_ptr->label_coords[0];
                  atom_ptr->label_coords[1]
                    = match_atom_ptr->label_coords[1];
                  atom_ptr->bvalue = match_atom_ptr->bvalue;
                  atom_ptr->occupancy = match_atom_ptr->occupancy;
                  atom_ptr->interacting = match_atom_ptr->interacting;

                  /* Replace the full atom in all the bond records */
                  replace_bond_atom(full_pdbentry_ptr,match_atom_ptr,
                                    atom_ptr);
                }
                  
              /* Get next atom record */
              atom_ptr = atom_ptr->next_atom_ptr;
              iatom++;
            }

          /* Mark the current residue as being in the plot */
          residue_ptr->in_plot = TRUE;

          /* Increment count of matched residues */
          nmatch++;
        }

      /* Get pointer to next residue record */
      residue_ptr = residue_ptr->next_residue_ptr;
    }

  /* Show number of residue matched */
  printf("%d residues matched from .drw file to original PDB file\n",
         nmatch);
}
/* <-- v.1.3 */
/***********************************************************************
  
check_bond  -  Check whether the given two atoms have a bond joining
               them

***********************************************************************/

int check_bond(struct bond *first_bond_ptr,struct atom *at1_ptr,
               struct atom *at2_ptr)
{
  int have_bond;

  struct atom *atom_ptr, *atom1_ptr, *atom2_ptr;
  struct bond *bond_ptr;

  /* Initialise variables */
  have_bond = FALSE;

  /* Save atom pointers */
  atom1_ptr = at1_ptr;
  atom2_ptr = at2_ptr;

  /* Swap atoms if necessary such that first atom is the lower-numbered
     one */
  if (atom2_ptr->atom_number < atom1_ptr->atom_number)
    {
      atom_ptr = atom1_ptr;
      atom1_ptr = atom2_ptr;
      atom2_ptr = atom_ptr;
    }

  /* Get pointer to the first bond record */
  bond_ptr = first_bond_ptr;

  /* Loop over the bond records to see if we have one for these two
     atoms */
  while (bond_ptr != NULL && have_bond == FALSE)
    {
      /* Check if this is the required record */
      if (bond_ptr->atom1_ptr == atom1_ptr &&
          bond_ptr->atom2_ptr == atom2_ptr)
        have_bond = TRUE;

      /* Get pointer to next bond record */
      bond_ptr = bond_ptr->next_bond_ptr;
    }

  /* Return whether have a bond record for these two atoms */
  return(have_bond);
}
/* v.2.0 --> */
/***********************************************************************

assign_domains  -  Assign residues to their interface domains

***********************************************************************/

void assign_domains(struct pdbentry *pdbentry_ptr,
                    struct domain *first_domain_ptr)
{
  int chain_type;

  struct domain *domain_ptr;
  struct residue *residue_ptr;

  /* Get the first domain */
  domain_ptr = first_domain_ptr;

  /* Loop through the domains to assign them to the structure's residues */
  while (domain_ptr != NULL)
    {
      /* Get this domain's chain type */
      if (domain_ptr->domain_number == 1)
        chain_type = CHAIN1;
      else if (domain_ptr->domain_number == 2)
        chain_type = CHAIN2;
      else
        chain_type = NON_INTERFACE;

      /* Get the first residue */
      residue_ptr = pdbentry_ptr->first_residue_ptr;

      /* Loop over the residues to mark those in the current domain */
      while (residue_ptr != NULL)
        {
          /* If this residue belongs to the current domain, save its
             chain type */
          if (residue_ptr->chain == domain_ptr->chain &&
              strcmp(residue_ptr->res_num,domain_ptr->start_residue_number)
              >= 0 &&
              strcmp(residue_ptr->res_num,domain_ptr->end_residue_number)
              <= 0)
            residue_ptr->chain_type = chain_type;

          /* Get pointer to next residue record */
          residue_ptr = residue_ptr->next_residue_ptr;
        }

      /* Get the next domain */
      domain_ptr = domain_ptr->next_domain_ptr;
    }
}
/* <-- v.2.0 */
/***********************************************************************

update_max_min  -  Update the given residue's maximum and minimum x-, y-,
                   and z-coords

***********************************************************************/

void update_max_min(struct residue *residue_ptr)
{
  int iatom, loop, natoms;
  int first[2];

  float x, y, z;

  struct atom *atom_ptr;

  /* Initialise variables */
  first[0] = first[1] = TRUE;

  /* Get the first atom of this residue */
  atom_ptr = residue_ptr->first_atom_ptr;
  natoms = residue_ptr->natoms;
  iatom = 0;

  /* Loop over this residue's atoms to update coord limits */
  while (atom_ptr != NULL && iatom < natoms)
    {
      /* Loop over the two types of coords: original and mapped */
      for (loop = 0; loop < 2; loop++)
        {
          /* Get the appropriate coordinates */
          if (loop == ORIGINAL)
            {
              x = atom_ptr->original_coords[X];
              y = atom_ptr->original_coords[Y];
              z = atom_ptr->original_coords[Z];
            }
          else
            {
              x = atom_ptr->coords[X];
              y = atom_ptr->coords[Y];
              z = atom_ptr->coords[Z];
            }

          /* If this is the first atom encountered, store its
             coords as the present max and min values */
          if (first[loop] == TRUE)
            {
              residue_ptr->valmin[loop][X]
                = residue_ptr->valmax[loop][X] = x;
              residue_ptr->valmin[loop][Y]
                = residue_ptr->valmax[loop][Y] = y;
              residue_ptr->valmin[loop][Z]
                = residue_ptr->valmax[loop][Z] = z;

              /* Unset flag */
              first[loop] = FALSE;
            }

          /* Otherwise, update the limits if coordinates are outside them */
          else
            {
              if (x < residue_ptr->valmin[loop][X])
                residue_ptr->valmin[loop][X] = x;
              if (x > residue_ptr->valmax[loop][X])
                residue_ptr->valmax[loop][X] = x;
              if (y < residue_ptr->valmin[loop][Y])
                residue_ptr->valmin[loop][Y] = y;
              if (y > residue_ptr->valmax[loop][Y])
                residue_ptr->valmax[loop][Y] = y;
              if (z < residue_ptr->valmin[loop][Z])
                residue_ptr->valmin[loop][Z] = z;
              if (z > residue_ptr->valmax[loop][Z])
                residue_ptr->valmax[loop][Z] = z;
            }
        }

      /* Get the next atom record and increment atom count */
      atom_ptr = atom_ptr->next_atom_ptr;
      iatom++;
    }
}
/***********************************************************************
  
update_all_maxmin -  Update all the residue max and in values

***********************************************************************/

void update_all_maxmin(struct pdbentry *first_pdbentry_ptr)
{
  int is_metal;

  struct atom *atom_ptr;
  struct pdbentry *pdbentry_ptr;
  struct residue *residue_ptr;

  /* Get the first PDB entry */
  pdbentry_ptr = first_pdbentry_ptr;

  /* Loop through the PDB entries */
  while (pdbentry_ptr != NULL)
    {
      /* Get this entry's first residue */
      residue_ptr = pdbentry_ptr->first_residue_ptr;

      /* Loop over the residues to update max-min coords */
      while (residue_ptr != NULL)
        {
          /* Update the coordinate limits */
          update_max_min(residue_ptr);

          /* If this residue is not a water, but has only a single
             atom, check if it is a metal */
          if (residue_ptr->type != WATER && residue_ptr->natoms == 1)
            {
              /* Get the atom */
              atom_ptr = residue_ptr->first_atom_ptr;

              /* Check if the atom is a metal */
              is_metal
                = check_for_metal(atom_ptr->atom_name,atom_ptr->element);

              /* If metal, then mark the residue as such */
              if (is_metal == TRUE)
                residue_ptr->type = METAL_ION;
            }

          /* Get pointer to next residue record */
          residue_ptr = residue_ptr->next_residue_ptr;
        }

      /* Get the next entry */
      pdbentry_ptr = pdbentry_ptr->next_pdbentry_ptr;
    }
}
/***********************************************************************
  
calc_bonds  -  Calculate all covalent bonds for the residues read in

***********************************************************************/

void calc_bonds(struct pdbentry *pdbentry_ptr)
{
/* v.2.0 --> */
//  int iatom, jatom, iresid, jresid, natoms1, natoms2;
  int iatom, jatom, natoms1, natoms2;
/* <-- v.2.0 */
  int bond_type, nbonds;
  int have_conect, in_range;

  float x1, x2, y1, y2, z1, z2;
  float distance_x, distance_y, distance_z, dist_comp, dist_sqrd;
  float length;

  struct atom *atom1_ptr, *atom2_ptr;
  struct bond *bond_ptr, *first_bond_ptr, *last_bond_ptr;
  struct residue *residue1_ptr, *residue2_ptr;

  /* Initailise variables */
  bond_ptr = first_bond_ptr = last_bond_ptr = NULL;
  nbonds = 0;

  /* Get the first residue */
  residue1_ptr = pdbentry_ptr->first_residue_ptr;

  /* Loop over all the residues to calculate all covalent bonds within
     the residue and to other residues */
  while(residue1_ptr != NULL)
    {
      /* Start at the same residue */
      residue2_ptr = residue1_ptr;

      /* Loop over all other residues, starting with the current one,
         to calculate covalent bonds */
      while(residue2_ptr != NULL)
        {
          /* If this is the same residue, then continue */
          if (residue2_ptr == residue1_ptr)
            in_range = TRUE;

          /* Otherwise, check whether these two residues are
             close enough for the atoms to be covalently bonded */
          else
            {
              /* Check whether the residue limits overlap */
              in_range = TRUE;
              if (residue1_ptr->valmin[ORIGINAL][0]
                  > residue2_ptr->valmax[ORIGINAL][0] + LONGBOND_DIST ||
                  residue1_ptr->valmax[ORIGINAL][0]
                  < residue2_ptr->valmin[ORIGINAL][0] - LONGBOND_DIST ||
                  residue1_ptr->valmin[ORIGINAL][1]
                  > residue2_ptr->valmax[ORIGINAL][1] + LONGBOND_DIST ||
                  residue1_ptr->valmax[ORIGINAL][1]
                  < residue2_ptr->valmin[ORIGINAL][1] - LONGBOND_DIST ||
                  residue1_ptr->valmin[ORIGINAL][2]
                  > residue2_ptr->valmax[ORIGINAL][2] + LONGBOND_DIST ||
                  residue1_ptr->valmax[ORIGINAL][2]
                  < residue2_ptr->valmin[ORIGINAL][2] - LONGBOND_DIST)
                in_range = FALSE;
            }

          /* If residues in range, calculate all-vs-all atom
             distances */
          if (in_range == TRUE)
            {
              /* Get the first of the first residue's atoms */
              atom1_ptr = residue1_ptr->first_atom_ptr;
              iatom = 0;
              natoms1 = residue1_ptr->natoms;

              /* Loop over all the atoms in the first residue */
              while (atom1_ptr != NULL && iatom < natoms1)
                {
                  /* Get this atom's coordinates */
                  x1 = atom1_ptr->coords[X];
                  y1 = atom1_ptr->coords[Y];
                  z1 = atom1_ptr->coords[Z];

                  /* If second residue is the same as the first,
                     then start atom is just the next atom in the
                     residue */
/* v.2.0 --> */
//                  if (iresid == jresid)
                  if (residue2_ptr == residue1_ptr)
/* <-- v.2.0 */
                    {
                      /* Get the next atom */
                      atom2_ptr = atom1_ptr->next_atom_ptr;
                      jatom = iatom + 1;
                    }

                  /* Otherwise, start with first atom in the residue */
                  else
                    {
                      /* Get the first atom */
                      atom2_ptr = residue2_ptr->first_atom_ptr;
                      jatom = 0;
                    }

                  /* Get the number of atoms in the second residue */
                  natoms2 = residue2_ptr->natoms;

                  /* Loop over all the atoms in the second residue */
                  while (atom2_ptr != NULL && jatom < natoms2)
                    {
                      /* Get this atom's coordinates */
                      x2 = atom2_ptr->coords[X];
                      y2 = atom2_ptr->coords[Y];
                      z2 = atom2_ptr->coords[Z];

                      /* Check whether have a CONECT record joining
                         these two atoms */
                      have_conect
                        = check_bond(pdbentry_ptr->first_bond_ptr,
                                     atom1_ptr,atom2_ptr);

                      /* Calculate distance between atoms 1 & 2 */
                      distance_x = (x1 - x2) * (x1 - x2);
                      distance_y = (y1 - y2) * (y1 - y2);
                      distance_z = (z1 - z2) * (z1 - z2);
                      dist_sqrd
                        = distance_x + distance_y + distance_z;

                      /* Calculate cut-off distance for covalent
                         connectivity, based on atom types involved */
                      dist_comp = COVALENT_DIST_SQRD;

                      /* Allow extra distance for bond to sulphur
                         and iron */
                      if (!strcmp(atom1_ptr->element," S") ||
                          !strcmp(atom1_ptr->element,"FE") ||
                          !strcmp(atom2_ptr->element," S") ||
                          !strcmp(atom2_ptr->element,"FE"))
                        dist_comp = LONGBOND_DIST_SQRD;

                      /* Ditto for phosphorus */
                      if (!strcmp(atom1_ptr->element," P") ||
                          !strcmp(atom2_ptr->element," P"))
                        dist_comp = PHOSBOND_DIST_SQRD;

                      /* Initialise bond type */
                      bond_type = DELETED;

                      /* If distance satisfies the distance criteria
                         for covalent bonding, then atoms are bonded */
/* v.2.0 --> */
//                      if (dist_sqrd < dist_comp)
                      if (dist_sqrd < dist_comp && atom1_ptr != atom2_ptr)
/* <-- v.2.0 */
                        {
                          /* Calculate the bond length */
                          length = (float) sqrt((double) dist_sqrd);
                              
                          /* If this is a disulphide bond, store it as
                             an H-bond */
                          if (!strcmp(atom1_ptr->element," S") &&
                              !strcmp(atom2_ptr->element," S"))
                            bond_type = HBOND;

                          /* Otherwise, store as a covalent bond */
                          else
                            bond_type = COVALENT;
                        }

                      /* If have a CONECT record, then check whether it
                         corresponds to a reasonable distance */
                      else if (have_conect == TRUE)
                        {
                          /* Define bond type and source */
                          bond_type = COVALENT;

                          /* Calculate the bond length */
                          length = (float) sqrt((double) dist_sqrd);

                          /* If too long, show warning message */
                          if (length > LONG_BOND)
                            {
                              /* Show warning */
                              printf("*** Warning. Long bond from "
                                     "CONECT [%s %s %s %c (%d) -> "
                                     "%s %s %s %c (%d)] "
                                     "Length %7.2f",
                                     atom1_ptr->atom_name,
                                     atom1_ptr->residue_ptr->res_name,
                                     atom1_ptr->residue_ptr->res_num,
                                     atom1_ptr->residue_ptr->chain,
                                     atom1_ptr->atom_number,
                                     atom2_ptr->atom_name,
                                     atom2_ptr->residue_ptr->res_name,
                                     atom2_ptr->residue_ptr->res_num,
                                     atom2_ptr->residue_ptr->chain,
                                     atom2_ptr->atom_number,length);

                              /* If length ridiculously too long, then
                                 discard */
                              if (length > BOND_TOO_LONG)
                                {
                                  printf(" - discarded\n");
                                  bond_type = DELETED;
                                }
                              else
                                printf("\n");
                            }
                        }

                      /* If have a bond, then create a bond record between
                         these two atoms */
                      if (bond_type != DELETED)
                        {
                          /* Create a bond record */
                          bond_ptr
                            = create_bond_record(&first_bond_ptr,
                                                 &last_bond_ptr,bond_type,
                                                 length,SINGLE,
                                                 atom1_ptr,atom2_ptr);

                          /* Increment count of this pdbentry's bonds */
                          nbonds++;
                        }

                      /* Get the next atom record from second residue */
                      atom2_ptr = atom2_ptr->next_atom_ptr;
                      jatom++;
                    }

                  /* Get the next atom record from first residue */
                  atom1_ptr = atom1_ptr->next_atom_ptr;
                  iatom++;
                }
            }

          /* Get the next residue */
          residue2_ptr = residue2_ptr->next_residue_ptr;
        }

      /* Get the next residue */
      residue1_ptr = residue1_ptr->next_residue_ptr;
    }

  /* Store bond details in pdbentry record */
  pdbentry_ptr->first_bond_ptr = first_bond_ptr;
  pdbentry_ptr->nbonds = nbonds;
  
  /* Show number of bonds created */
  printf("Number of covalent bonds  %8d\n",nbonds);
}
/***********************************************************************

get_hbplus_info  -  Extract the atom details from the interaction just
                    read in from the .hhb or .nnb file

***********************************************************************/

int get_hbplus_info(char *line,char *atom_name1,char *atom_name2,
                    char *chain1,char *chain2,char *res_name1,
                    char *res_name2,char *res_num1,char *res_num2,
                    float *bond_length,int file_type,int ligplot_format)
{
  char bond_lenstr[5];

  int ipos, nonzero, valid;

  /* Initialise variables */
  valid = TRUE;

  /* Check for obviously unwanted records */
  if ((int) strlen(line) < 10 ||
      (ligplot_format == TRUE && !strncmp(line + 21,"   ",3)) ||
      (ligplot_format == FALSE && !strncmp(line + 10,"   ",3)))
    valid = FALSE;

  /* Otherwise, extract the data according to the file format */
  else
    {
      /* LIGPLOT-format of the .hhb and .nnb files */
      if (ligplot_format == TRUE)
        {
          /* Atom name */
          strncpy(atom_name1,line+12,4);
          atom_name1[4] = '\0';
          strncpy(atom_name2,line+33,4);
          atom_name2[4] = '\0';
                
          /* Chain-id */
          *chain1 = line[4];
          *chain2 = line[25];
          if ((*chain1) == '-')
            *chain1 = ' ';
          if ((*chain2) == '-')
            *chain2 = ' ';
          
          /* Residue details */
          strncpy(res_num1,line+6,5);
          res_num1[5] = '\0';
          strncpy(res_num2,line+27,5);
          res_num2[5] = '\0';
          strncpy(res_name1,line,3);
          res_name1[3] = '\0';
          strncpy(res_name2,line+21,3);
          res_name2[3] = '\0';

          /* Bond length */
          strncpy(bond_lenstr,line+41,4);
          bond_lenstr[4] = '\0';
          *bond_length = (float) atof(bond_lenstr);
        }

      /* HBPLUS-format hydrogen bonds file */
      else
        {
          /* Atom name */
          strncpy(atom_name1,line+19,4);
          atom_name1[4] = '\0';
          strncpy(atom_name2,line+39,4);
          atom_name2[4] = '\0';
                
          /* Chain-id */
          *chain1 = line[10];
          *chain2 = line[30];
          if ((*chain1) == '-')
            *chain1 = ' ';
          if ((*chain2) == '-')
            *chain2 = ' ';
          
          /* Residue details */
          strncpy(res_num1,line+11,5);
          res_num1[5] = '\0';
          strncpy(res_num2,line+31,5);
          res_num2[5] = '\0';
          strncpy(res_name1,line+16,3);
          res_name1[3] = '\0';
          strncpy(res_name2,line+36,3);
          res_name2[3] = '\0';

          /* Bond length */
          strncpy(bond_lenstr,line+45,4);
          bond_lenstr[4] = '\0';
          *bond_length = (float) atof(bond_lenstr);
        }

      /* For the residue numbers, may need to replace any leading
         zeros by spaces */

      /* First residue */
      nonzero = FALSE;
      for (ipos = 0; ipos < 3 && nonzero == FALSE; ipos ++)
        {
          if (res_num1[ipos] == '0')
            res_num1[ipos] = ' ';
          else
            nonzero = TRUE;
        }

      /* Replace hyphen in insertion code with a space */
      if (res_num1[4] == '-')
        res_num1[4] = ' ';
                
      /* Repeat for second residue */
      nonzero = FALSE;
      for (ipos = 0; ipos < 3 && nonzero == FALSE; ipos ++)
        {
          if (res_num2[ipos] == '0')
            res_num2[ipos] = ' ';
          else
            nonzero = TRUE;
        }

      /* Replace hyphen in insertion code with a space */
      if (res_num2[4] == '-')
        res_num2[4] = ' ';
    }

  /* Return status */
  return(valid);
}
/***********************************************************************

match_residue  -  Find the given residue record

***********************************************************************/

struct residue *match_residue(struct residue *first_residue_ptr,
                              char *res_name,char *res_num,char chain)
{
  struct residue *residue_ptr, *match_residue_ptr;

  /* Initialise variables */
  match_residue_ptr = NULL;

  /* Get pointer to the first residue */
  residue_ptr = first_residue_ptr;

  /* Loop over the residues to find match */
  while (residue_ptr != NULL && match_residue_ptr == NULL)
    {
      /* Check match on residue */
      if (!strcmp(residue_ptr->res_num,res_num) &&
          residue_ptr->chain == chain &&
          !strcmp(residue_ptr->res_name,res_name))
        match_residue_ptr = residue_ptr;

      /* Get pointer to next residue record */
      residue_ptr = residue_ptr->next_residue_ptr;
    }

  /* Return match */
  return(match_residue_ptr);
}
/***********************************************************************

match_atom  -  Find the given atom record

***********************************************************************/

struct atom *match_atom(struct atom *first_atom_ptr,int natoms,
                        char *atom_name)
{
  int iatom;

  struct atom *atom_ptr, *match_atom_ptr;

  /* Initialise variables */
  iatom = 0;
  match_atom_ptr = NULL;

  /* Get pointer to the first atom */
  atom_ptr = first_atom_ptr;

  /* Loop over the atoms to write out */
  while (atom_ptr != NULL && iatom < natoms && match_atom_ptr == NULL)
    {
      /* Check if atom name matches */
      if (!strcmp(atom_ptr->atom_name,atom_name))
        match_atom_ptr = atom_ptr;

      /* Get pointer to next atom record */
      atom_ptr = atom_ptr->next_atom_ptr;
      iatom++;
    }
      
  /* Return match */
  return(match_atom_ptr);
}
/***********************************************************************

read_hbplus_file  -  Read in the required H-bond and nonbonded contacts
                     data from the appropriate HBPLUS output file

***********************************************************************/

int read_hbplus_file(char *file_name,int bond_type,
                     struct pdbentry *pdbentry_ptr,
                     struct parameters params)
{
  char bond_desc[20], input_line[LINELEN + 1];
  char atom_name1[5], atom_name2[5], chain1, chain2, res_name1[4];
  char res_name2[4], res_num1[6], res_num2[6];
  char element[3];

  int have_file, valid, wanted;
  int nbonds, nline;

  float bond_length;

  struct atom *atom1_ptr, *atom2_ptr;
  struct bond *bond_ptr, *first_bond_ptr, *last_bond_ptr;
  struct residue *residue1_ptr, *residue2_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  have_file = TRUE;
  nbonds = 0;
  nline = 0;
  if (bond_type == HBOND)
    strcpy(bond_desc,"hydrogen bond");
  else
    strcpy(bond_desc,"non-bonded contact");

  /* Get the first of this sequence's bonds */
  bond_ptr = first_bond_ptr = pdbentry_ptr->first_bond_ptr;
  last_bond_ptr = NULL;

  /* Loop through the bonds to find the last bond */
  while (bond_ptr != NULL)
    {
      /* Store this bond as the current last one */
      last_bond_ptr = bond_ptr;

      /* Get the next bond */
      bond_ptr = bond_ptr->next_bond_ptr;
    }

  /* Open the HBPLUS results file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    {
      printf("*** Warning. No %ss returned by HBPLUS for %s\n",
             bond_desc,pdbentry_ptr->pdb_code); 
      return(0);
    }

  /* Read in the HBPLUS data */
  while (fgets(input_line,LINELEN,file_ptr) != NULL)
    {
      /* Initialise variables and increment line-count */
      atom1_ptr = atom2_ptr = NULL;
      nline++;
      wanted = FALSE;

      /* Extract the relevant information from the line just
         read in */
      valid = get_hbplus_info(input_line,atom_name1,atom_name2,&chain1,
                              &chain2,res_name1,res_name2,res_num1,
                              res_num2,&bond_length,bond_type,FALSE);

      /* If this is a valid interaction check whether it is wanted */
      if (valid == TRUE)
        {
          /* Identify both residues involved in this bond */
          residue1_ptr
            = match_residue(pdbentry_ptr->first_residue_ptr,res_name1,
                            res_num1,chain1);
          residue2_ptr
            = match_residue(pdbentry_ptr->first_residue_ptr,res_name2,
                            res_num2,chain2);

          /* Identify both atoms involved in this bond */
          if (residue1_ptr != NULL && residue2_ptr != NULL)
            {
              atom1_ptr
                = match_atom(residue1_ptr->first_atom_ptr,
                             residue1_ptr->natoms,atom_name1);
              atom2_ptr
                = match_atom(residue2_ptr->first_atom_ptr,
                             residue2_ptr->natoms,atom_name2);
            }
        }

      /* If both atoms identified, determine whether this bond is
         wanted */
      if (atom1_ptr != NULL && atom2_ptr != NULL)
        {
          /* Initialise flag which determines whether this interaction
             is required */
          wanted = TRUE;

          /* For non-bonded contacts, check whether this is a desired
             contact */
          if (bond_type == CONTACT)
            {
              /* If hydrophobics only (as in the old version) check
                 that both atoms are either a carbon or a sulphur */
              if (params.contact_type == HYDROPHOBIC_ONLY)
                {
                  if ((strcmp(atom1_ptr->element," C") &&
                       strcmp(atom1_ptr->element," S")) ||
                      (strcmp(atom2_ptr->element," C") &&
                       strcmp(atom2_ptr->element," S")))
                    wanted = FALSE;
                }

              /* If at least one of the atoms has to be a hydrophobic
                 then check for this */
              else if (params.contact_type == HYDROPHOBIC_ANY)
                {
                  if (!strcmp(atom1_ptr->element," C") ||
                      !strcmp(atom1_ptr->element," S") ||
                      !strcmp(atom2_ptr->element," C") ||
                      !strcmp(atom2_ptr->element," S"))
                    wanted = TRUE;
                  else
                    wanted = FALSE;
                }

              /* Otherwise, accept any contact */
              else
                wanted = TRUE;

              /* If either atom is a water, then don't want a
                 non-bonded contact */
              if (!strcmp(res_name1,"HOH") ||
                  !strcmp(res_name2,"HOH"))
                wanted = FALSE;

              /* If either is a metal, then again don't want this
                 contact */
              if (wanted == TRUE)
                {
                  if (check_for_metal(atom_name1,element) == TRUE)
                    wanted = FALSE;
                  else if (check_for_metal(atom_name2,element) == TRUE)
                    wanted = FALSE;
                }
            }

          /* Check if either residue is in the ligand */
          if (wanted == TRUE)
            {
              /* If both residues are ligand residues, then don't
                 want this interaction */
              if (residue1_ptr->type == LIGAND &&
                  residue2_ptr->type == LIGAND)
                wanted = FALSE;

/* v.2.0 --> */
              /* For DIMPLOT, if both residues belong to the same
                 chain, or to an unwanted chain, then don't want 
                 this bond */
              if (params.dimplot == TRUE)
                {
                  if ((residue1_ptr->chain_type == CHAIN1 &&
                       residue2_ptr->chain_type == CHAIN1) ||
                      (residue1_ptr->chain_type == CHAIN2 &&
                       residue2_ptr->chain_type == CHAIN2))
                    wanted = FALSE;
                  else if (residue1_ptr->chain_type == NON_INTERFACE ||
                           residue2_ptr->chain_type == NON_INTERFACE)
                    wanted = FALSE;
                }
/* <-- v.2.0 */
            }

          /* If this interaction is still wanted, then store it */
          if (wanted == TRUE)
            {
              /* Create a bond record */
              bond_ptr
                = create_bond_record(&first_bond_ptr,
                                     &last_bond_ptr,bond_type,
                                     bond_length,SINGLE,
                                     atom1_ptr,atom2_ptr);

              /* Increment count of this sequence's bonds */
              nbonds++;
            }
        }
    }

  /* Close the input file */
  fclose(file_ptr);

  /* Save this PDB entry's first bond pointer */
  pdbentry_ptr->first_bond_ptr = first_bond_ptr;
  pdbentry_ptr->nbonds = pdbentry_ptr->nbonds + nbonds;

  /* Show number of bonds read in from file */
  printf("Number of %ss read in from HBPLUS file: %8d\n",bond_desc,
         nbonds);

  /* Return the number of bonds stored */
  return(nbonds);
}
/***********************************************************************
  
check_interactions  -  Mark all residues that interact directly with the
                       ligand (or, for DIMPLOT, with the other chain)

***********************************************************************/

/* v.2.0 --> */
//void check_interactions(struct pdbentry *first_pdbentry_ptr)
void check_interactions(struct pdbentry *first_pdbentry_ptr,
                        int dimplot)
/* <-- v.2.0 */
{
  int interaction, tmp_int;

  struct atom *atom_ptr, *other_atom_ptr;
  struct bond *bond_ptr;
  struct pdbentry *pdbentry_ptr;
  struct residue *residue_ptr, *other_residue_ptr;

  /* Get the first PDB entry */
  pdbentry_ptr = first_pdbentry_ptr;

  /* Loop through the PDB entries */
  while (pdbentry_ptr != NULL)
    {
      /* Get this entry's first residue */
      residue_ptr = pdbentry_ptr->first_residue_ptr;

      /* Loop over the residues to mark the ligands and waters */
      while (residue_ptr != NULL)
        {
         /* Initialise this residue's interaction type */
          interaction = NO_INTERACTION;

          /* If residue is a ligand, mark its interaction type
             accordingly */
          if (residue_ptr->type == LIGAND)
            interaction = IS_LIGAND;

          /* Otherwise, check whether it interacts directly with the
             ligand */
          else
            {
              /* Get the first of this PDB entry's bonds */
              bond_ptr = pdbentry_ptr->first_bond_ptr;

              /* Loop over the bonds to find any involving the current
                 residue */
              while (bond_ptr != NULL)
                {
                  /* Initialise atom pointers */
                  atom_ptr = other_atom_ptr = NULL;

                  /* If either atom belongs to the current residue, then
                     store */
                  if (bond_ptr->atom1_ptr->residue_ptr == residue_ptr &&
                      bond_ptr->atom2_ptr->residue_ptr != residue_ptr)
                    {
                      /* Save the atoms */
                      atom_ptr = bond_ptr->atom1_ptr;
                      other_atom_ptr = bond_ptr->atom2_ptr;
                    }

                  else if (bond_ptr->atom2_ptr->residue_ptr == residue_ptr &&
                           bond_ptr->atom1_ptr->residue_ptr != residue_ptr)
                    {
                      /* Save the atoms */
                      atom_ptr = bond_ptr->atom2_ptr;
                      other_atom_ptr = bond_ptr->atom1_ptr;
                    }

                  /* If the other atom belongs to the ligand, store the
                     interaction type */
                  if (other_atom_ptr != NULL &&
                      other_atom_ptr->residue_ptr->type == LIGAND)
                    {
                      /* If have a water interaction, then store */
                      if (bond_ptr->type == HBOND &&
                          residue_ptr->type == WATER)
                        interaction = WATER_TO_LIGAND;

                      /* Ditto for a metal interaction */
                      else if (residue_ptr->type == METAL_ION)
                        interaction = METAL_TO_LIGAND;

                      /* Otherwise, for ordinary residue, store if
                         this is the highest-classed interaction type
                         so far */
                      else if (residue_ptr->type != WATER &&
                               residue_ptr->type != METAL_ION)
                        {
                          /* Store potential interaction type */
/* v.1.3 --> */
                          tmp_int = NO_INTERACTION;
/* <-- v.1.3 */
                          if (bond_ptr->type == COVALENT)
                            tmp_int = COVALENT_TO_LIGAND;
                          else if (bond_ptr->type == HBOND)
                            tmp_int = HBOND_TO_LIGAND;
                          else if (bond_ptr->type == CONTACT)
                            tmp_int = NONBOND_TO_LIGAND;

                          /* If this is "higher" than the interaction
                             we have, store it */
                          if (tmp_int > interaction)
                            interaction = tmp_int;
                        }
                    }

/* v.2.0 --> */
                  /* For DIMPLOT, check for interaction between the
                     two chains/domains */
                  else if (other_atom_ptr != NULL && dimplot == TRUE)
                    {
                      /* Get the other residue */
                      other_residue_ptr = other_atom_ptr->residue_ptr;

                      /* If have a water interaction, then store */
                      if (bond_ptr->type == HBOND &&
                          residue_ptr->type != other_residue_ptr->type &&
                          (residue_ptr->type == WATER ||
                           other_residue_ptr->type == WATER) &&
                          interaction != PROT_PROT_HBOND)
                        interaction = VIA_WATER;

                      /* Ditto for a metal interaction */
                      else if ((residue_ptr->type == METAL_ION ||
                               other_residue_ptr->type == METAL_ION) &&
                               interaction != PROT_PROT_HBOND)
                        interaction = VIA_METAL;

                      /* Otherwise, for ordinary residue, store if
                         this is the highest-classed interaction type
                         so far */
                      else if ((residue_ptr->chain_type == CHAIN1 &&
                                other_residue_ptr->chain_type == CHAIN2) ||
                               (residue_ptr->chain_type == CHAIN2 &&
                                other_residue_ptr->chain_type == CHAIN1))
                        {
                          /* Store potential interaction type */
                          tmp_int = NO_INTERACTION;
                          if (bond_ptr->type == HBOND)
                            tmp_int = PROT_PROT_HBOND;
                          else if (bond_ptr->type == CONTACT)
                            tmp_int = PROT_PROT_NONBOND;

                          /* If this is "higher" than the interaction
                             we have, store it */
                          if (tmp_int > interaction)
                            interaction = tmp_int;
                        }
                    }
/* <-- v.2.0 */

                  /* Get the next bond */
                  bond_ptr = bond_ptr->next_bond_ptr;
                }
            }

          /* Store the interaction type */
/* v.2.0 --> */
          if (interaction > residue_ptr->interaction)
/* <-- v.2.0 */
            residue_ptr->interaction = interaction;

          /* If residue makes an interaction, flag as being in the plot */
          if (interaction != NO_INTERACTION)
            {
              /* Set flag */
              residue_ptr->in_plot = TRUE;

              /* Determine the residue type */
              if (interaction == PROT_PROT_NONBOND)
                residue_ptr->type = HYDROPHOBIC;
            }

          /* Get pointer to next residue record */
          residue_ptr = residue_ptr->next_residue_ptr;
        }

      /* Get the next entry */
      pdbentry_ptr = pdbentry_ptr->next_pdbentry_ptr;
    }
}
/***********************************************************************
  
classify_residues  -  Classify all the residues by the type of
                      interaction they make

***********************************************************************/

void classify_residues(struct pdbentry *first_pdbentry_ptr)
{
  int found;

  struct atom *atom_ptr, *other_atom_ptr;
  struct bond *bond_ptr;
  struct pdbentry *pdbentry_ptr;
  struct residue *residue_ptr;

  /* Get the first PDB entry */
  pdbentry_ptr = first_pdbentry_ptr;

  /* Loop through the PDB entries */
  while (pdbentry_ptr != NULL)
    {
      /* Get this entry's first residue */
      residue_ptr = pdbentry_ptr->first_residue_ptr;

      /* Loop over the residues to look for residues that are still
         marked as non-interacting */
      while (residue_ptr != NULL)
        {
          /* If not interacting directly with the ligand, see if it
             interacts via a water*/
/* v.1.3--> */
//          if (residue_ptr->interaction == NO_INTERACTION &&
//              residue_ptr->type != WATER &&
//              residue_ptr->type != METAL_ION)
          if (residue_ptr->interaction == NO_INTERACTION ||
              residue_ptr->interaction == NONBOND_TO_LIGAND)
/* <--v.1.3 */
            {
              /* Get the first of this PDB entry's bonds */
              bond_ptr = pdbentry_ptr->first_bond_ptr;
              found = FALSE;

              /* Loop over the bonds to find any involving the current
                 residue */
              while (bond_ptr != NULL && found == FALSE)
                {
                  /* Initialise atom pointers */
                  atom_ptr = other_atom_ptr = NULL;

                  /* If either atom belongs to the current residue, then
                     store */
                  if (bond_ptr->atom1_ptr->residue_ptr == residue_ptr)
                    {
                      /* Save the atoms */
                      atom_ptr = bond_ptr->atom1_ptr;
                      other_atom_ptr = bond_ptr->atom2_ptr;
                    }

                  else if (bond_ptr->atom2_ptr->residue_ptr == residue_ptr)
                    {
                      /* Save the atoms */
                      atom_ptr = bond_ptr->atom2_ptr;
                      other_atom_ptr = bond_ptr->atom1_ptr;
                    }

                  /* If the other atom is a ligand-interacting water
                     then flag current residue as interacting via a water */
                  if (other_atom_ptr != NULL &&
                      other_atom_ptr->residue_ptr->interaction
                      == WATER_TO_LIGAND)
                    {
                      /* Flag as interacting via a water */
                      residue_ptr->interaction = VIA_WATER_TO_LIGAND;

                      /* Set flag that interaction found */
                      found = TRUE;
                    }

/* v.1.3 --> */
                  /* If the other atom is a ligand-interacting metal
                     then flag current residue as interacting via a metal */
                  if (other_atom_ptr != NULL &&
                      other_atom_ptr->residue_ptr->interaction
                      == METAL_TO_LIGAND)
                    {
                      /* Flag as interacting via a water */
                      residue_ptr->interaction = VIA_METAL_TO_LIGAND;

                      /* Set flag that interaction found */
                      found = TRUE;
                    }
/* <-- v.1.3 */

                  /* Get the next bond */
                  bond_ptr = bond_ptr->next_bond_ptr;
                }
            }

          /* Get pointer to next residue record */
          residue_ptr = residue_ptr->next_residue_ptr;
        }

      /* Get the next entry */
      pdbentry_ptr = pdbentry_ptr->next_pdbentry_ptr;
    }
}
/***********************************************************************
  
delete_unwanted_residues  -  Delete any residues that do not make any
                             interactions

***********************************************************************/

int delete_unwanted_residues(struct pdbentry *pdbentry_ptr)
{
  int nresid;

  struct residue *residue_ptr, *next_residue_ptr, *prev_residue_ptr;

  /* Initialise variables */
  nresid = 0;

  /* Get this entry's first residue */
  residue_ptr = pdbentry_ptr->first_residue_ptr;
  prev_residue_ptr = NULL;

  /* Loop through the residues until done */
  while (residue_ptr != NULL)
    {
      /* Get pointer to the next residue record */
      next_residue_ptr = residue_ptr->next_residue_ptr;

      /* If residue makes no interactions, then delete */
      if (residue_ptr->interaction == NO_INTERACTION)
        {
          /* If this is the first residue, then adjust pointer to
             point to next residue */
          if (prev_residue_ptr == NULL)
            pdbentry_ptr->first_residue_ptr = next_residue_ptr;

          /* Otherwise, get previous residue to point round it */
          else
            prev_residue_ptr->next_residue_ptr = next_residue_ptr;
        }

      /* Otherwise, save as previous, interacting residue */
      else
        {
          prev_residue_ptr = residue_ptr;
          nresid++;
        }

      /* Go to the next residue */
      residue_ptr = next_residue_ptr;
    }

  /* Return number of residues remaining */
  return(nresid);
}
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

get_aa_name  -  Get residue-name from single-letter amino acid code

***********************************************************************/

void get_aa_name(char aa_code,char res_name[4])
{
  int iamino, ichar;
  int namino;

  static char amino[] = "ACDEFGHIKLMNPQRSTVWYXBZXX";
  static char *amino_name[] = {
    "ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU",
    "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR",
    "PCA ","ASX ","GLX ", "UNK", "XXX"
  };

  /* Initialise variables */
  namino = strlen(amino) - 1;

  /* Get the corresponding three-letter code */
  iamino = namino;
  for (ichar = 0; ichar < namino; ichar++)
    if (aa_code == amino[ichar])
      iamino = ichar;
  strncpy(res_name,amino_name[iamino],3);
  res_name[3] = '\0';
}
/***********************************************************************

create_resseq_record  -  Create and initialise a new resseq record

***********************************************************************/

struct resseq *create_resseq_record(struct resseq **fst_resseq_ptr,
                                    struct resseq **lst_resseq_ptr,
                                    struct residue *residue_ptr)
{
  struct resseq *resseq_ptr;
  struct resseq *last_resseq_ptr, *first_resseq_ptr;

  /* Initialise resseq pointers */
  resseq_ptr = NULL;
  last_resseq_ptr = *lst_resseq_ptr;
  first_resseq_ptr = *fst_resseq_ptr;

  /* Allocate memory for structure to hold sequence info */
  resseq_ptr = (struct resseq *) malloc(sizeof(struct resseq));
  if (resseq_ptr == NULL)
    {
      printf("*** Can't allocate memory for struct");
      printf(" resseq\n");
      exit (1);
    }

  /* Store the details for this entry */
  resseq_ptr->residue_ptr = residue_ptr;

  /* Initialise remaining fields */
  resseq_ptr->number = 0;
  resseq_ptr->align_resseq_ptr = NULL;
  resseq_ptr->next_resseq_ptr = NULL;

  /* If this is the first entry stored, then update pointer to head
     of linked list */
  if (first_resseq_ptr == NULL)
    first_resseq_ptr = resseq_ptr;

  /* Otherwise, make previous entry point to the current
     one */
  else
    last_resseq_ptr->next_resseq_ptr = resseq_ptr;
                      
  /* Save the current resseq's pointer */
  last_resseq_ptr = resseq_ptr;

  /* Return current pointers */
  *fst_resseq_ptr = first_resseq_ptr;
  *lst_resseq_ptr = last_resseq_ptr;
  return(resseq_ptr);
}
/***********************************************************************

check_backbone  -  Check whether the given residue has all 4 backbone
                   atoms

***********************************************************************/

int check_backbone(struct residue *residue_ptr)
{
  char has_atom[5];

  int iatom, natoms;
  int has_backbone = FALSE;

  struct atom *atom_ptr;

  /* Initialise atom pattern */
  strcpy(has_atom,"    ");

  /* Get the first atom of this residue */
  atom_ptr = residue_ptr->first_atom_ptr;
  natoms = residue_ptr->natoms;
  iatom = 0;

  /* Loop over this residue's atoms to check for the backbone atoms */
  while (atom_ptr != NULL && iatom < natoms)
    {
      /* Check for backbone atom, setting appropriate character if found */
      if (strcmp(atom_ptr->atom_name," N  "))
        has_atom[0] = '*';
      else if (strcmp(atom_ptr->atom_name," CA "))
        has_atom[1] = '*';
      else if (strcmp(atom_ptr->atom_name," C  "))
        has_atom[2] = '*';
      else if (strcmp(atom_ptr->atom_name," O  "))
        has_atom[3] = '*';

      /* Get the next atom record and increment atom count */
      atom_ptr = atom_ptr->next_atom_ptr;
      iatom++;
    }

  /* If all backbone atoms present, set flag */
  if (!strcmp(has_atom,"****"))
    has_backbone = TRUE;

  /* Return whether all 4 backbone atoms are present */
  return(has_backbone);
}
/***********************************************************************

compare_interactions  -  Check whether the given atoms make equivalent
                         interactions in their respective structures

***********************************************************************/

int compare_interactions(struct atom *atom1_ptr,struct atom *atom2_ptr)
{
  int match = FALSE;

/* @@@ To do */
/* debug
  printf("In compare_interactions for  %s %s %s %c from %d and "
         "%s %s %s %c from %d\n",
         atom1_ptr->atom_name,
         atom1_ptr->residue_ptr->res_name,
         atom1_ptr->residue_ptr->res_num,
         atom1_ptr->residue_ptr->chain,
         atom1_ptr->residue_ptr->pdbentry_ptr->pdb_id,
         atom2_ptr->atom_name,
         atom2_ptr->residue_ptr->res_name,
         atom2_ptr->residue_ptr->res_num,
         atom2_ptr->residue_ptr->chain,
         atom2_ptr->residue_ptr->pdbentry_ptr->pdb_id);
 debug */

  /* Return whether we have a match */
  return(match);
}
/***********************************************************************

get_fit_coords  -  Retrieve the coords of the two residues in the
                   given alignment position that can be fitted

***********************************************************************/

int get_fit_coords(struct resseq *resseq_ptr[2],
                   float fit_coords[MAX_RES_FIT_COORDS][2][3])
{
  int iatom, jatom, natoms1, natoms2, nfit_coords;
  int compare;
  int aacid1, aacid2, found, wanted;

  struct atom *atom1_ptr, *atom2_ptr;
  struct residue *residue1_ptr, *residue2_ptr;

  /* Initialise flag */
  nfit_coords = 0;

  /* Get the two residues */
  residue1_ptr = resseq_ptr[REF_SEQUENCE]->residue_ptr;
  residue2_ptr = resseq_ptr[CURR_SEQUENCE]->residue_ptr;

  /* If residues are of identical types, then compare all atoms */
  if (!strcmp(residue1_ptr->res_name,residue2_ptr->res_name))
    compare = ALL_ATOMS;

  /* Otherwise, check if one or both are standard amino acids */
  else
    {
      /* See if either residue is a standard amino acid */
/* v.1.2 --> */
/*      if (strstr(residue1_ptr->res_name,AMINO_ACIDS) != NULL) */
      if (strstr(AMINO_ACIDS,residue1_ptr->res_name) != NULL)
/* <-- v.1.2 */
        aacid1 = TRUE;
      else
        aacid1 = check_backbone(residue1_ptr);
/* v.1.2 --> */
/*      if (strstr(residue2_ptr->res_name,AMINO_ACIDS) != NULL) */
      if (strstr(AMINO_ACIDS,residue2_ptr->res_name) != NULL)
/* <-- v.1.2 */
        aacid2 = TRUE;
      else
        aacid2 = check_backbone(residue2_ptr);

      /* If both have backbone atoms, then check on these */
      if (aacid1 == TRUE && aacid2 == TRUE)
        compare = BACKBONE_ONLY;

      /* Otherwise, match only interacting atoms */
      else
        compare = INTERACTING_ONLY;
    }

  /* Get the first atom of the first residue */
  atom1_ptr = residue1_ptr->first_atom_ptr;
  natoms1 = residue1_ptr->natoms;
  iatom = 0;

  /* Loop over this residue's atoms to write out coords of any matches to
     atoms in the other residue */
  while (atom1_ptr != NULL && iatom < natoms1)
    {
      /* Get the first atom of the other residue */
      atom2_ptr = residue2_ptr->first_atom_ptr;
      natoms2 = residue2_ptr->natoms;
      jatom = 0;
      found = FALSE;

      /* Loop over this residue's atoms to find a match */
      while (atom2_ptr != NULL && jatom < natoms2 && found == FALSE)
        {
          /* Initialise wanted flag */
          wanted = FALSE;

          /* If comparing the two atom's interactions, see if they are
             equivalent */
          if (compare == INTERACTING_ONLY)
            {
              if (atom1_ptr->interacting == TRUE &&
                  atom2_ptr->interacting == TRUE)
                {
                  /* Compare the interactions */
                  wanted = compare_interactions(atom1_ptr,atom2_ptr);
                }
              else
                wanted = FALSE;
            }

          /* Otherwise, see if atom names match */
          else if (!strcmp(atom1_ptr->atom_name,atom2_ptr->atom_name))
            {
              /* If only comparing on backbone atoms, check this is one */
              if (compare == BACKBONE_ONLY)
                {
                  /* If a backbone atom, then match is wanted */
                  if (!strcmp(atom1_ptr->atom_name," N  ") ||
                      !strcmp(atom1_ptr->atom_name," CA ") ||
                      !strcmp(atom1_ptr->atom_name," C  ") ||
                      !strcmp(atom1_ptr->atom_name," O  "))
                    wanted = TRUE;
                  else
                    wanted = FALSE;
                }

              /* If only comparing on C-alpha atoms, check this is one */
              else if (compare == CALPHA_ONLY)
                {
                  /* If a C-alpha atom, then match is wanted */
                  if (!strcmp(atom1_ptr->atom_name," CA "))
                    wanted = TRUE;
                  else
                    wanted = FALSE;
                }
              else
                wanted = TRUE;
            }

          /* If this pair of atoms match, then store coords */
          if (wanted == TRUE && nfit_coords < MAX_RES_FIT_COORDS)
            {
              /* Store this atom's coordinates */
              fit_coords[nfit_coords][REF_SEQUENCE][X]
                = atom1_ptr->original_coords[X];
              fit_coords[nfit_coords][REF_SEQUENCE][Y]
                = atom1_ptr->original_coords[Y];
              fit_coords[nfit_coords][REF_SEQUENCE][Z]
                = atom1_ptr->original_coords[Z];
              fit_coords[nfit_coords][CURR_SEQUENCE][X]
                = atom2_ptr->original_coords[X];
              fit_coords[nfit_coords][CURR_SEQUENCE][Y]
                = atom2_ptr->original_coords[Y];
              fit_coords[nfit_coords][CURR_SEQUENCE][Z]
                = atom2_ptr->original_coords[Z];

              /* Increment count of stored coords */
              nfit_coords++;

              /* Set flag that match found */
              found = TRUE;
            }

          /* Get the next atom record and increment atom count */
          atom2_ptr = atom2_ptr->next_atom_ptr;
          jatom++;
        }

      /* Get the next atom record and increment atom count */
      atom1_ptr = atom1_ptr->next_atom_ptr;
      iatom++;
    }

  /* Return number of pairs of matched coords */
  return(nfit_coords);
}
/***********************************************************************

centre_coords  -  Calculate the C of G of the given set of coords
                  and adjust

***********************************************************************/

void centre_coords(COOR *coords,int nstored,double *cofg)
{
  int istored;

  double mean_x, mean_y, mean_z;

  /* Initialise variables */
  mean_x = 0.0;
  mean_y = 0.0;
  mean_z = 0.0;

  /* Loop over the stored coordinates */
  for (istored = 0; istored < nstored; istored++)
    {
      /* Accumulate mean values */
      mean_x = mean_x + coords[istored].x;
      mean_y = mean_y + coords[istored].y;
      mean_z = mean_z + coords[istored].z;
    }

  /* Calculate mean values */
  mean_x = mean_x / nstored;
  mean_y = mean_y / nstored;
  mean_z = mean_z / nstored;

  /* Loop through the stored coordinates to centre them at the origin */
  for (istored = 0; istored < nstored; istored++)
    {
      /* Adjust both sets of coords */
      coords[istored].x = coords[istored].x - mean_x;
      coords[istored].y = coords[istored].y - mean_y;
      coords[istored].z = coords[istored].z - mean_z;
    }

  /* Store the initial CofG */
  cofg[0] = mean_x;
  cofg[1] = mean_y;
  cofg[2] = mean_z;
}
/***********************************************************************

Andrew Martin's matfit and qikfit routines for least-squares fitting of
one set of coordinates onto another (both sets centred around the origin)

***********************************************************************/
/*>static void qikfit(REAL umat[3][3], REAL rm[3][3], BOOL column)
   ---------------------------------------------------------------
   Input:   REAL  umat[3][3]     The U matrix
            BOOL  column         TRUE: Create a column-wise matrix
                                 (other way round from normal).
   Output:  REAL  rm[3][3]       The output rotation matrix
  
   Does the actual fitting for matfit().
   04.02.91 Original
   01.06.92 ANSIed & doc'd
   11.03.94 column changed to BOOL
*/
static void qikfit(REAL  umat[3][3],
                   REAL  rm[3][3],
                   BOOL  column)
{
   
   REAL  rot[3][3],
         turmat[3][3],
         c[3][3],
         coup[3],
         dir[3],
         step[3],
         v[3],
         rtsum,rtsump,
         stp,stcoup,
         ud,tr,ta,cs,sn,ac,
         delta,
         gfac,
         cle,clep;
   int   i,j,k,l,m,
         jmax,
         ncyc,
         nsteep,
         nrem;

   /* Rotate repeatedly to reduce couple about initial direction
      to zero.
      Clear the rotation matrix */
   for(l=0;l<3;l++)
   {
      for(m=0;m<3;m++)
         rot[l][m] = 0.0;
      rot[l][l] = 1.0;
   }

   /* Copy vmat[][] (sp) into umat[][] (dp) */
   jmax = 30;
   rtsum = umat[0][0] + umat[1][1] + umat[2][2];
   delta = 0.0;

   for(ncyc=0;ncyc<jmax;ncyc++)
   {
      /* Modified CG. For first and every NSTEEP cycles, set previous
         step as zero and do an SD step */
      nsteep = 3;
      nrem = ncyc-nsteep*(int)(ncyc/nsteep);

      if(!nrem)
      {
         for(i=0;i<3;i++) step[i]=0.0;
         clep = 1.0;
      }
      
      /* Couple */
      coup[0] = umat[1][2]-umat[2][1];
      coup[1] = umat[2][0]-umat[0][2];
      coup[2] = umat[0][1]-umat[1][0];
      cle = sqrt(coup[0]*coup[0] + coup[1]*coup[1] + coup[2]*coup[2]);

      /* Gradient vector is now -coup */
      gfac = (cle/clep)*(cle/clep);

      /* Value of rtsum from previous step */
      rtsump = rtsum;
      clep   = cle;
      if(cle < SMALL) break;

      /* Step vector conjugate to  previous */
      stp = 0.0;
      for(i=0;i<3;i++)
      {
         step[i]=coup[i]+step[i]*gfac;
         stp   += (step[i] * step[i]);
      }
      stp = 1.0/sqrt(stp);
         
      /* Normalised step */
      for(i=0;i<3;i++) dir[i] = stp*step[i];

      /* Couple resolved along step direction */
      stcoup = coup[0]*dir[0] + coup[1]*dir[1] + coup[2]*dir[2];

      /* Component of UMAT along direction */
      ud = 0.0;
      for(l=0;l<3;l++)
         for(m=0;m<3;m++)
            ud += umat[l][m]*dir[l]*dir[m];


      tr = umat[0][0]+umat[1][1]+umat[2][2]-ud;
      ta = sqrt(tr*tr + stcoup*stcoup);
      cs=tr/ta;
      sn=stcoup/ta;
         
      /* If cs<0 then posiiton is unstable, so don't stop */
      if((cs>0.0) && (ABS(sn)<SMALSN)) break;
            
      /* Turn matrix for correcting rotation:

         Symmetric part */
      ac = 1.0-cs;
      for(l=0;l<3;l++)
      {
         v[l] = ac*dir[l];
         for(m=0;m<3;m++)
            turmat[l][m] = v[l]*dir[m];
         turmat[l][l] += cs;
         v[l]=dir[l]*sn;
      }

      /* Asymmetric part */
      turmat[0][1] -= v[2];
      turmat[1][2] -= v[0];
      turmat[2][0] -= v[1];
      turmat[1][0] += v[2];
      turmat[2][1] += v[0];
      turmat[0][2] += v[1];

      /* Update total rotation matrix */
      for(l=0;l<3;l++)
      {
         for(m=0;m<3;m++)
         {
            c[l][m] = 0.0;
            for(k=0;k<3;k++)
               c[l][m] += turmat[l][k]*rot[k][m];
         }
      }

      for(l=0;l<3;l++)
         for(m=0;m<3;m++)
            rot[l][m] = c[l][m];

      /* Update umat tensor */
      for(l=0;l<3;l++)
         for(m=0;m<3;m++)
         {
            c[l][m] = 0.0;
            for(k=0;k<3;k++)
               c[l][m] += turmat[l][k]*umat[k][m];
         }

      for(l=0;l<3;l++)
         for(m=0;m<3;m++)
            umat[l][m] = c[l][m];

      rtsum = umat[0][0] + umat[1][1] + umat[2][2];
      delta = rtsum - rtsump;

      /* If no improvement in this cycle then stop */
      if(ABS(delta)<SMALL) break;

      /* Next cycle */
   }

   /* Copy rotation matrix for output */
   if(column)
   {
      for(i=0;i<3;i++)
         for(j=0;j<3;j++)
            rm[j][i] = rot[i][j];
   }
   else
   {
      for(i=0;i<3;i++)
         for(j=0;j<3;j++)
            rm[i][j] = rot[i][j];
   }
}
/**********************************************************************/
/**********************************************************************/
/*>BOOL matfit(COOR *x1, COOR *x2, REAL rm[3][3], int n,
               REAL *wt1, BOOL column)
   -----------------------------------------------------
   Input:   COOR  *x1         First array of coordinates
            COOR  *x2         Second array of coordinates
            int   n           Number of coordinates
            REAL  *wt1        Weight array or NULL
            int   column      TRUE: Output a column-wise matrix (as used
                                 by FRODO)
                              FALSE: Output a standard row-wise matrix.
   Output:  REAL  rm[3][3]    Returned rotation matrix
   Returns: rmsd              RMS distance between the two sets of
                              coordinates after fitting

   Fit coordinate array x1 to x2 both centred around the origin and of 
   length n. Optionally weighted with the wt1 array if iw flag is set.
   If flag is set, error messages will be issued internally. 
   If column is set the matrix will be returned column-wise rather
   than row-wise.

   04.02.91 Original
   01.06.92 ANSIed & doc'd
   17.06.93 various changes for release (including parameters)
   11.03.94 column changed to BOOL
   04.03.97 Calc of rmsd added (RAL)
*/
double matfit(COOR    *x1,        /* First coord array    */
             COOR    *x2,        /* Second coord array   */
             REAL    rm[3][3],   /* Rotation matrix      */
             int     n,          /* Number of points     */
             REAL    *wt1,       /* Weight array         */
             BOOL    column)     /* Column-wise output   */
{
   int  i, ipoint, j;
   double dist, rmsd, x, y, z;

   REAL umat[3][3];

   if(wt1)
   {
      for(i=0;i<3;i++)
      {
         for(j=0;j<3;j++) umat[i][j] = 0.0;

         for(j=0;j<n;j++)
         {
            switch(i)
            {
               case 0:
                  umat[i][0] += wt1[j] * x1[j].x * x2[j].x;
                  umat[i][1] += wt1[j] * x1[j].x * x2[j].y;
                  umat[i][2] += wt1[j] * x1[j].x * x2[j].z;
                  break;
               case 1:
                  umat[i][0] += wt1[j] * x1[j].y * x2[j].x;
                  umat[i][1] += wt1[j] * x1[j].y * x2[j].y;
                  umat[i][2] += wt1[j] * x1[j].y * x2[j].z;
                  break;
               case 2:
                  umat[i][0] += wt1[j] * x1[j].z * x2[j].x;
                  umat[i][1] += wt1[j] * x1[j].z * x2[j].y;
                  umat[i][2] += wt1[j] * x1[j].z * x2[j].z;
                  break;
            }
         }
      }
   } 
   else
   {
      for(i=0;i<3;i++)
      {
         for(j=0;j<3;j++) umat[i][j] = 0.0;

         for(j=0;j<n;j++)
         {
            switch(i)
            {
               case 0:
                  umat[i][0] += x1[j].x * x2[j].x;
                  umat[i][1] += x1[j].x * x2[j].y;
                  umat[i][2] += x1[j].x * x2[j].z;
                  break;
               case 1:
                  umat[i][0] += x1[j].y * x2[j].x;
                  umat[i][1] += x1[j].y * x2[j].y;
                  umat[i][2] += x1[j].y * x2[j].z;
                  break;
               case 2:
                  umat[i][0] += x1[j].z * x2[j].x;
                  umat[i][1] += x1[j].z * x2[j].y;
                  umat[i][2] += x1[j].z * x2[j].z;
                  break;
            }
         }
      }
   }
   qikfit(umat,rm,column);

   /* Calculate rmsd between the two fitted sets of coordinates */
   rmsd = 0.0;

   /* Loop over all the pairs of coordinates */
   for (ipoint = 0; ipoint < n; ipoint++)
     {
       /* Calculate the transformed coordinates of the current point */
       if (column)
         {
           x = rm[0][0] * x1[ipoint].x + rm[1][0] * x1[ipoint].y
             + rm[2][0] * x1[ipoint].z;
           y = rm[0][1] * x1[ipoint].x + rm[1][1] * x1[ipoint].y
             + rm[2][1] * x1[ipoint].z;
           z = rm[0][2] * x1[ipoint].x + rm[1][2] * x1[ipoint].y
             + rm[2][2] * x1[ipoint].z;
         }
       else
         {
           x = rm[0][0] * x1[ipoint].x + rm[0][1] * x1[ipoint].y
             + rm[0][2] * x1[ipoint].z;
           y = rm[1][0] * x1[ipoint].x + rm[1][1] * x1[ipoint].y
             + rm[1][2] * x1[ipoint].z;
           z = rm[2][0] * x1[ipoint].x + rm[2][1] * x1[ipoint].y
             + rm[2][2] * x1[ipoint].z;
         }

       /* Calculate squared distance between these two points */
       dist = (x - x2[ipoint].x) * (x - x2[ipoint].x)
         + (y - x2[ipoint].y) * (y - x2[ipoint].y)
         + (z - x2[ipoint].z) * (z - x2[ipoint].z);
       rmsd = rmsd + dist;
     }

   /* Calculate the rmsd */
   if (n > 0)
     {
       rmsd = rmsd / n;
       if (rmsd > 0.0)
         rmsd = sqrt(rmsd);
       else
         rmsd = 0.0;
     }

   /* Return the calculated rmsd */
   return(rmsd);
}
/***********************************************************************

transform_atom  -  Apply the given transformation to residue coords

***********************************************************************/

void transform_atom(double x,double y,double z,double matrix[3][3],
                    double cofg[2][3],double *newx,double *newy,
                    double *newz)
{
  /* Adjust coords being adjusted to CofG frame of reference */
  x = x - cofg[CURR_SEQUENCE][0];
  y = y - cofg[CURR_SEQUENCE][1];
  z = z - cofg[CURR_SEQUENCE][2];

  /* Apply the rotation */
  *newx = x * matrix[0][0] + y * matrix[0][1] + z * matrix[0][2];
  *newy = x * matrix[1][0] + y * matrix[1][1] + z * matrix[1][2];
  *newz = x * matrix[2][0] + y * matrix[2][1] + z * matrix[2][2];

  /* Move relative to C of G of frame we are fitting to */
  *newx = *newx + cofg[REF_SEQUENCE][0];
  *newy = *newy + cofg[REF_SEQUENCE][1];
  *newz = *newz + cofg[REF_SEQUENCE][2];
}
/***********************************************************************

calc_separation  -  Calculate the separation (squared) between the
                    C-alpha atoms of the given residue pair

***********************************************************************/

double calc_separation(float fit_coords[MAX_RES_FIT_COORDS][2][3],
                       int nfit_coords,
/* v.2.0.1 --> */
//                       double matrix[3][3],double cofg[2][3])
                       double matrix[3][3],double cofg[2][3],
                       struct parameters params)
/* <-- v.2.0.1 */
{
  int i;
/* v.2.0.1 --> */
//  int verbose = FALSE;
/* <-- v.2.0.1 */

  double ave_dist, dist, dist_sqrd;
  double x1, x2, xp, y1, y2, yp, z1, z2, zp;

  /* Initialise average distance */
  ave_dist = 0.0;

  /* Loop over the paired atoms */
  for (i = 0; i < nfit_coords; i++)
    {
      /* Get coords of the two residues */
      xp = fit_coords[i][CURR_SEQUENCE][X];
      yp = fit_coords[i][CURR_SEQUENCE][Y];
      zp = fit_coords[i][CURR_SEQUENCE][Z];
      x1 = fit_coords[i][REF_SEQUENCE][X];
      y1 = fit_coords[i][REF_SEQUENCE][Y];
      z1 = fit_coords[i][REF_SEQUENCE][Z];

      /* Apply transformation to other residue */
      transform_atom(xp,yp,zp,matrix,cofg,&x2,&y2,&z2);

      /* Calculate the distance between these two atoms */
      dist_sqrd
        = (x1 - x2) * (x1 - x2)
        + (y1 - y2) * (y1 - y2)
        + (z1 - z2) * (z1 - z2);
      dist = sqrt(dist_sqrd);
/* v.2.0.1 --> */
//      if (verbose == TRUE)
      if (params.verbose == TRUE)
/* <-- v.2.0.1 */
        {
          printf("Coords: %8.3f %8.3f %8.3f and %8.3f %8.3f %8.3f\n",
                 xp,yp,zp,x1,y1,z1);
          printf("After :                                "
                 "%8.3f %8.3f %8.3f (%8.3f)\n",x2,y2,z2,dist);
        }

      /* Accumulate */
      ave_dist = ave_dist + dist;
    }

  /* Calculate average distance */
  if (nfit_coords > 0)
    ave_dist = ave_dist / nfit_coords;

  /* Return the average distance separation */
  return(ave_dist);
}
/***********************************************************************

fit_coords  -  Fit the residues in 3D using the sequence alignment
               between the interacting residues in our PDB entry and
               those in the LIGPLOT from the ensemble

***********************************************************************/

int fit_coords(struct resseq **alignment,int align_length,
/* v.1.0.4 --> */
/*               REAL matrix[3][3],double cofg[2][3]) */
               REAL matrix[3][3],double cofg[2][3],int *nt,
/* v.2.0.1 --> */
//               double min_sep_cutoff)
               double min_sep_cutoff,struct parameters params)
/* <-- v.2.0.1 */
/* <-- v.1.0.4 */
{
  int i, ipos, j, nfit_coords, npaired, nstored;
  int new_int, old_int;
/* v.1.0.4 --> */
  int ntried;
/* <-- v.1.0.4 */
/* v.2.0.1 --> */
//  int verbose = FALSE;
/* <-- v.2.0.1 */

  float fit_coords[MAX_RES_FIT_COORDS][2][3];

  double ave_dist, rmsd, rms_cutoff;
/* v.1.0.4 --> */
  double percentage;
/* <-- v.1.0.4 */

  BOOL column;

  COOR coords1[MAX_FIT_RES], coords2[MAX_FIT_RES];

  struct residue *residue_ptr;
  struct resseq *resseq_ptr[2];

  /* Initialise variables */
  rmsd = 0.0;
  column = FALSE;
  npaired = 0;
  nstored = 0;
  rms_cutoff = RMSD_CUTOFF;
/* v.1.0.4 --> */
  ntried = 0;
/* <-- v.1.0.4 */

  /* Initialise transformation matrices */
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      if (i == j)
        matrix[i][j] = 1.0;
      else
        matrix[i][j] = 0.0;
  for (i = 0; i < 2; i++)
    for (j = 0; j < 3; j++)
      cofg[i][j] = 0.0;

  /* Loop through all the positions in the alignment to store
     the C-alpha coords of the residues */
  for (ipos = 0; ipos < align_length; ipos++)
    {
      /* Get the corresponding residue pointers at this point in
         the alignment */
      resseq_ptr[REF_SEQUENCE] = alignment[ipos];
      resseq_ptr[CURR_SEQUENCE] = alignment[ipos + align_length];

      /* If have two alignment residues, get their corresponding coords
         to add to the set for fitting */
      if (resseq_ptr[REF_SEQUENCE] != NULL &&
          resseq_ptr[CURR_SEQUENCE] != NULL &&
          resseq_ptr[REF_SEQUENCE]->residue_ptr->type != DUMMY &&
          resseq_ptr[CURR_SEQUENCE]->residue_ptr->type != DUMMY)
        {
          /* If in verbose mode, show the residues being fitted */
/* v.2.0.1 --> */
//          if (verbose == TRUE)
          if (params.verbose == TRUE)
/* <-- v.2.0.1 */
            printf("Fitting %s %s %c from %s to %s %s %c from %s",
                   resseq_ptr[CURR_SEQUENCE]->residue_ptr->res_name,
                   resseq_ptr[CURR_SEQUENCE]->residue_ptr->res_num,
                   resseq_ptr[CURR_SEQUENCE]->residue_ptr->chain,
                   resseq_ptr[CURR_SEQUENCE]
                   ->residue_ptr->pdbentry_ptr->pdb_code,
                   resseq_ptr[REF_SEQUENCE]->residue_ptr->res_name,
                   resseq_ptr[REF_SEQUENCE]->residue_ptr->res_num,
                   resseq_ptr[REF_SEQUENCE]->residue_ptr->chain,
                   resseq_ptr[REF_SEQUENCE]
                   ->residue_ptr->pdbentry_ptr->pdb_code);
/* v.1.0.4 --> */
          /* Initialise these residues' mappings */
          resseq_ptr[REF_SEQUENCE]->residue_ptr->map_residue_ptr = NULL;
          resseq_ptr[CURR_SEQUENCE]->residue_ptr->map_residue_ptr = NULL;
/* <-- v.1.0.4 */

          /* Get the two residues' coords for fitting */
          nfit_coords = get_fit_coords(resseq_ptr,fit_coords);
/* v.2.0.1 --> */
//          if (verbose == TRUE)
          if (params.verbose == TRUE)
/* <-- v.2.0.1 */
            printf(" nfit_coords = %d\n",nfit_coords);

          /* Loop over the pairs of matched coords to store */
          for (i = 0; i < nfit_coords && nstored < MAX_FIT_RES; i++)
            {
              /* Store the coords */
              coords1[nstored].x = fit_coords[i][CURR_SEQUENCE][X];
              coords1[nstored].y = fit_coords[i][CURR_SEQUENCE][Y];
              coords1[nstored].z = fit_coords[i][CURR_SEQUENCE][Z];
              coords2[nstored].x = fit_coords[i][REF_SEQUENCE][X];
              coords2[nstored].y = fit_coords[i][REF_SEQUENCE][Y];
              coords2[nstored].z = fit_coords[i][REF_SEQUENCE][Z];

              /* Increment count of stored coords */
              nstored++;
            }
        }
    }

  /* Centre coords of PDB entry */
  centre_coords(coords1,nstored,cofg[CURR_SEQUENCE]);

  /* Repeat for LIGPLOT entry */
  centre_coords(coords2,nstored,cofg[REF_SEQUENCE]);

  /* If have more than one pair of atoms stored, can perform the fit */
  if (nstored > 1)
    {
      /* Perform the least-squares fit of coords1 -> coords2 */
      rmsd = matfit(coords1,coords2,matrix,nstored,NULL,column);
      printf("  %d coords fitted with rmsd = %8.3f\n",nstored,rmsd);

      /* Loop through all the alignment positions again */
      for (ipos = 0; ipos < align_length; ipos++)
        {
          /* Get the two residues */
          resseq_ptr[REF_SEQUENCE] = alignment[ipos];
          resseq_ptr[CURR_SEQUENCE] = alignment[ipos + align_length];

          /* If have two alignment residues, get their appropriate coords
             to superpose and calc separation distance */
          if (resseq_ptr[REF_SEQUENCE] != NULL &&
              resseq_ptr[CURR_SEQUENCE] != NULL &&
              resseq_ptr[REF_SEQUENCE]->residue_ptr->type != DUMMY &&
              resseq_ptr[CURR_SEQUENCE]->residue_ptr->type != DUMMY)
            {
              /* Get the two residues' coords for fitting */
              nfit_coords = get_fit_coords(resseq_ptr,fit_coords);
            }
          else
            nfit_coords = 0;
          
          /* If have both C-alphas, check apply fitting transformation
             and check their separation */
          if (nfit_coords > 0)
            {
/* v.1.0.4 --> */
              /* Increment count of residue pairs tried */
              ntried++;
/* <-- v.1.0.4 */
/* v.2.0.1 --> */
//              if (verbose == TRUE)
              if (params.verbose == TRUE)
/* <-- v.2.0.1 */
                {
                  printf("Fitted  %s %s %c of %s  to  %s %s %c of %s  ",
                         resseq_ptr[CURR_SEQUENCE]->residue_ptr->res_name,
                         resseq_ptr[CURR_SEQUENCE]->residue_ptr->res_num,
                         resseq_ptr[CURR_SEQUENCE]->residue_ptr->chain,
                         resseq_ptr[CURR_SEQUENCE]
                         ->residue_ptr->pdbentry_ptr->pdb_code,
                         resseq_ptr[REF_SEQUENCE]->residue_ptr->res_name,
                         resseq_ptr[REF_SEQUENCE]->residue_ptr->res_num,
                         resseq_ptr[REF_SEQUENCE]->residue_ptr->chain,
                         resseq_ptr[REF_SEQUENCE]
                         ->residue_ptr->pdbentry_ptr->pdb_code);
                }

              /* Calculate the average separation between the two residues'
                 matched atoms */
              ave_dist
/* v.2.0.1 --> */
//                = calc_separation(fit_coords,nfit_coords,matrix,cofg);
//              if (verbose == TRUE)
                = calc_separation(fit_coords,nfit_coords,matrix,cofg,params);
              if (params.verbose == TRUE)
/* <-- v.2.0.1 */
                printf("Ave dist = %8.3f (%8.3f)\n",ave_dist,
                       min_sep_cutoff);

              /* If close enough to constitute an overlap, link the
                 two residues together */
/* v.1.0.4 --> */
/*              if (ave_dist < MIN_SEP_CUTOFF) */
              if (ave_dist < min_sep_cutoff)
/* <-- v.1.0.4 */
                {
                  /* Get our residue */
                  residue_ptr = resseq_ptr[CURR_SEQUENCE]->residue_ptr;

                  /* Initialise interaction type of LIGPLOT residue */
                  old_int
                    = new_int
                    = resseq_ptr[REF_SEQUENCE]->residue_ptr->interaction;

                  /* If our residue already mapped, see if this is a
                     better mapping */
                  if (residue_ptr->map_residue_ptr != NULL)
                    {
                      /* Get the residue type */
                      old_int
                        = residue_ptr->map_residue_ptr->interaction;
                      new_int
                        = resseq_ptr[REF_SEQUENCE]->residue_ptr->interaction;
                    }

                  /* If no current mapping, or this gives a better one,
                     perform the mapping */
                  if (residue_ptr->map_residue_ptr == NULL ||
                      new_int > old_int)
                    {
                      residue_ptr->map_residue_ptr
                        = resseq_ptr[REF_SEQUENCE]->residue_ptr;
                      resseq_ptr[REF_SEQUENCE]->residue_ptr->map_residue_ptr
                        = residue_ptr;

                      /* Increment count of mapped residues */
                      npaired++;
                    }
                }
/* v.1.0.4 --> */
              /* If this residue pair doesn't fit too well, then blank out
                 this equivalence in the alignment */
              else
                {
                  alignment[ipos] = NULL;
                  alignment[ipos + align_length] = NULL;
                }
/* <-- v.1.0.4 */
            }
        }
    }

/* v.1.0.4 --> */
  /* If number of residues fitted is fewer than the minimum, then set
     count to zero */
  if (npaired < MIN_ALIGNED)
    npaired = 0;

  /* Check that the percentage of fitted to aligned residues is above
     minimum */
  if (ntried > 0)
    percentage = 100.0 * (float) npaired / (float) ntried;
  else
    percentage = 0;
  if (percentage < MIN_PERCENTAGE)
    npaired = 0;

  /* Return number of pairs tried */
  *nt = ntried;
/* <-- v.1.0.4 */

  /* Return number of paired residues */
  return(npaired);
}
/***********************************************************************

find_atomap  -  Find the given atomap record

***********************************************************************/

struct atomap *find_atomap(struct atomap *first_atomap_ptr,
                           struct atom *atom_ptr)
{
  struct atomap *atomap_ptr, *match_atomap_ptr;

  /* Initialise variables */
  match_atomap_ptr = NULL;

  /* Get pointer to the first atomap */
  atomap_ptr = first_atomap_ptr;

  /* Loop over the atomaps to write out */
  while (atomap_ptr != NULL && match_atomap_ptr == NULL)
    {
      /* Check if this is the required record */
      if (atomap_ptr->map_atom_ptr == atom_ptr)
        {
          /* Save the atomap pointer */
          match_atomap_ptr = atomap_ptr;
        }

      /* Get pointer to next atomap record */
      atomap_ptr = atomap_ptr->next_atomap_ptr;
    }
      
  /* Return match */
  return(match_atomap_ptr);
}
/***********************************************************************

create_atomap_record  -  Create a new atomap record

***********************************************************************/

struct atomap *create_atomap_record(struct atomap **fst_atomap_ptr,
                                    struct atomap **lst_atomap_ptr,
                                    struct atom *atom_ptr,
                                    struct atom *map_atom_ptr,int score)
{
  struct atomap *atomap_ptr;
  struct atomap *first_atomap_ptr, *last_atomap_ptr;

  /* Initialise atomap pointers */
  atomap_ptr = NULL;
  last_atomap_ptr = *lst_atomap_ptr;
  first_atomap_ptr = *fst_atomap_ptr;

  /* Create atomap entry */
  atomap_ptr = (struct atomap*)malloc(sizeof(struct atomap));
  if (atomap_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for");
      printf(" struct atomap\n");
      printf("***        Program terminated with error.\n");
      exit (1);
    }

  /* Store the details for this entry */
  atomap_ptr->atom_ptr = atom_ptr;
  atomap_ptr->map_atom_ptr = map_atom_ptr;
  atomap_ptr->score = score;

  /* Initialise remaining fields */
  atomap_ptr->anchor_weight = 1.0;

  /* Initialise pointer to next atomap record */
  atomap_ptr->next_atomap_ptr = NULL;

  /* If this is the first entry stored, then update
     pointer to head of linked list */
  if (first_atomap_ptr == NULL)
    first_atomap_ptr = atomap_ptr;

  /* Otherwise, make previous entry point to the current
     one */
  else
    last_atomap_ptr->next_atomap_ptr = atomap_ptr;
                      
  /* Save the current atomap's pointer */
  last_atomap_ptr = atomap_ptr;

  /* Return current pointers */
  *fst_atomap_ptr = first_atomap_ptr;
  *lst_atomap_ptr = last_atomap_ptr;
  return(atomap_ptr);
}
/***********************************************************************

check_equiv_hbonds  -  Check whether the mapped atoms both have H-bonds
                       to equivalenced residues

***********************************************************************/

int check_equiv_hbonds(struct atom *atom_ptr,struct atom *match_atom_ptr)
{
  int equiv;

  struct atom *other_atom_ptr;
  struct bond *bond_ptr, *match_bond_ptr;
  struct residue *residue_ptr, *match_residue_ptr, *map_residue_ptr;

  /* Initialise variables */
  equiv = FALSE;

  /* Get the parent residues */
  residue_ptr = atom_ptr->residue_ptr;
  match_residue_ptr = match_atom_ptr->residue_ptr;

  /* Get pointer to the first bond record */
  bond_ptr = residue_ptr->pdbentry_ptr->first_bond_ptr;
  
  /* Loop over the bonds to find those involving the current atom */
  while (bond_ptr != NULL)
    {
      /* Initialise interacting atom */
      other_atom_ptr = NULL;

      /* If bond involves current atom, save it */
      if (bond_ptr->atom1_ptr == atom_ptr)
        other_atom_ptr = bond_ptr->atom2_ptr;
      else if (bond_ptr->atom2_ptr == atom_ptr)
        other_atom_ptr = bond_ptr->atom1_ptr;

      /* If have a hydrogen bond involving our atom, check for an 
         equivalent H-bond to the match atom */
      if (other_atom_ptr != NULL && bond_ptr->type == HBOND &&
          residue_ptr != other_atom_ptr->residue_ptr)
        {
          /* Save the mapping residue */
          map_residue_ptr = other_atom_ptr->residue_ptr->map_residue_ptr;

          /* Get pointer to the first bond record */
          match_bond_ptr = match_residue_ptr->pdbentry_ptr->first_bond_ptr;
  
          /* Loop over the bonds to find those involving the current atom */
          while (match_bond_ptr != NULL && map_residue_ptr != NULL)
            {
              /* Initialise interacting atom */
              other_atom_ptr = NULL;

              /* If bond involves current atom, save it */
              if (match_bond_ptr->atom1_ptr == match_atom_ptr)
                other_atom_ptr = match_bond_ptr->atom2_ptr;
              else if (match_bond_ptr->atom2_ptr == match_atom_ptr)
                other_atom_ptr = match_bond_ptr->atom1_ptr;

              /* If have a hydrogen bond involving our atom, check for an 
                 equivalent H-bond to the match atom */
              if (other_atom_ptr != NULL &&
                  match_bond_ptr->type == HBOND &&
                  match_residue_ptr != other_atom_ptr->residue_ptr)
                {
                  /* If the interacting residue is the mapped residue
                     from the alignment, then original atoms are
                     equivalent */
                  /*                  if (!strcmp(other_atom_ptr->residue_ptr->res_name,
                              residue_ptr->map_residue_ptr->res_name) &&
                      !strcmp(other_atom_ptr->residue_ptr->res_num,
                              residue_ptr->map_residue_ptr->res_num) &&
                      other_atom_ptr->residue_ptr->chain
                      == residue_ptr->map_residue_ptr->chain) */
                  if (other_atom_ptr->residue_ptr == map_residue_ptr)
                    equiv = TRUE;
                }

              /* Get pointer to next bond record */
              match_bond_ptr = match_bond_ptr->next_bond_ptr;
            }
        }

      /* Get pointer to next bond record */
      bond_ptr = bond_ptr->next_bond_ptr;
    }

  /* Return whether the atoms have equivalent interactions */
  return(equiv);
}
/***********************************************************************

open_output_file  -  Open the output file

***********************************************************************/

void open_output_file(FILE **out_ptr,char *out_name)
{
  FILE *outfile_ptr;

  /* Open the output file, or write to stdout */
  if ((outfile_ptr = fopen(out_name,"w")) ==  NULL)
    {
      printf("\n*** Unable to open output file [%s]\n",out_name);
      exit(1);
    }

  /* Return the pointer to the output file */
  *out_ptr = outfile_ptr;
}
/***********************************************************************

residue_overlap  -  Check for level of overlap between the given two
                    residues

***********************************************************************/

int residue_overlap(struct residue *residue1_ptr,
                    struct residue *residue2_ptr,REAL matrix[3][3],
                    double cofg[2][3],struct atomap **fst_atomap_ptr,
/* v.1.3.7 --> */
//                    struct atomap **lst_atomap_ptr,int *natm,int store)
                    struct atomap **lst_atomap_ptr,int *natm,int store,
                    FILE *file_ptr,struct parameters params)
/* <-- v.1.3.7 */
{
  int best_match, iatom, jatom, natoms, natoms1, natoms2, noverlap;
/* v.1.3 --> */
  int factor, min_best_match;
/* <-- v.1.3 */
  int in_range, overlap;

  double dist_sqrd, min_dist_sqrd, x, y, z;

  struct atom *atom1_ptr, *atom2_ptr, *best_atom_ptr;
  struct atomap *atomap_ptr, *first_atomap_ptr, *last_atomap_ptr;

  /* Initialise variables */
  overlap = FALSE;
  noverlap = 0;
  natoms = *natm;

  /* Initialise pointers */
  first_atomap_ptr = *fst_atomap_ptr;
  last_atomap_ptr = *lst_atomap_ptr;

  /* Get the first atom */
  atom1_ptr = residue1_ptr->first_atom_ptr;
  natoms1 = residue1_ptr->natoms;
  iatom = 0;

  /* Loop over the residue's atoms to compute their transformed coords
     and check for overlap with the atoms of the other residue */
  while (atom1_ptr != NULL && iatom < natoms1)
    {
      /* Transform this atom's coords */
      transform_atom(atom1_ptr->original_coords[X],
                     atom1_ptr->original_coords[Y],
                     atom1_ptr->original_coords[Z],matrix,cofg,&x,&y,&z);

       /* Check if the atom falls within the LIGPLOT residue's limits */
      in_range = TRUE;
      if (x > residue2_ptr->valmax[ORIGINAL][X] + 2 * MAX_RADIUS ||
          x < residue2_ptr->valmin[ORIGINAL][X] - 2 * MAX_RADIUS ||
          y > residue2_ptr->valmax[ORIGINAL][Y] + 2 * MAX_RADIUS ||
          y < residue2_ptr->valmin[ORIGINAL][Y] - 2 * MAX_RADIUS ||
          z > residue2_ptr->valmax[ORIGINAL][Z] + 2 * MAX_RADIUS ||
          z < residue2_ptr->valmin[ORIGINAL][Z] - 2 * MAX_RADIUS)
        in_range = FALSE;

      /* If in range, see if atom overlaps any of the other residue's
         atoms */
      if (in_range == TRUE)
        {
          /* Get the first atom of the LIGPLOT residue */
          atom2_ptr = residue2_ptr->first_atom_ptr;
          natoms2 = residue2_ptr->natoms;
          jatom = 0;
          overlap = FALSE;
          best_match = min_dist_sqrd = 0;
/* v.1.3 --> */
          min_best_match = 0;
/* <-- v.1.3 */
          best_atom_ptr = NULL;

          /* Loop over the residue's atoms */
          while (atom2_ptr != NULL && jatom < natoms2)
            {
              /* Calc distance between these two atoms */
              dist_sqrd
                = (x - atom2_ptr->original_coords[X])
                * (x - atom2_ptr->original_coords[X])
                + (y - atom2_ptr->original_coords[Y])
                * (y - atom2_ptr->original_coords[Y])
                + (z - atom2_ptr->original_coords[Z])
                * (z - atom2_ptr->original_coords[Z]);

              /* If atoms overlap, then set flag */
              if (dist_sqrd < MAX_RADIUS2)
                {
/* v.1.3 --> */
                  /* Factor in the similarity of the two atoms */
                  factor = 1;
                  if (strcmp(atom1_ptr->element,atom2_ptr->element))
                    factor = DISSIMILARITY_FACTOR;
                  best_match = (int) (factor * 1000 * dist_sqrd);

                  /* If this is the best match score so far, save the
                     atom */
                  if (overlap == FALSE ||
                      (overlap == TRUE && best_match < min_best_match))

//                  /* If this is the closest distance so far, save the
//                     atom */
//                  if (overlap == FALSE ||
//                      (overlap == TRUE && dist_sqrd < min_dist_sqrd))
/* <-- v.1.3 */
                    {
                      best_atom_ptr = atom2_ptr;
                      min_dist_sqrd = dist_sqrd;
/* v.1.3 --> */
//                      best_match = (int) (1000 * min_dist_sqrd);
                      min_best_match = best_match;
/* <-- v.1.3 */
                    }

                  /* Set flag */
                  overlap = TRUE;
                }

              /* Get the next atom record and increment atom count */
              atom2_ptr = atom2_ptr->next_atom_ptr;
              jatom++;
            }

          /* If we have an overlap, then increment count of overlapping
             atoms and save the mapping, if required */
          if (overlap == TRUE)
            {
              /* Increment count */
              noverlap++;

              /* If need to save the mapping, create an atomap record */
              if (store == TRUE)
                {
                  /* Check whether we already have a mapping to the mapped
                     atom */
                  atomap_ptr = find_atomap(first_atomap_ptr,best_atom_ptr);

                  /* If not found, create a new map record */
                  if (atomap_ptr == NULL)
                    {
                      /* Create the atom mapping */
                      atomap_ptr
                        = create_atomap_record(&first_atomap_ptr,
                                               &last_atomap_ptr,
                                               atom1_ptr,best_atom_ptr,
                                               best_match);

/* v.1.3.7 --> */
                      /* Write out the atom mapping */
                      fprintf(file_ptr,"MAPPED  %s %s %s %c (%s) "
                             " to %s %s %s %c (%s)\n",
                             atom1_ptr->atom_name,
                             atom1_ptr->residue_ptr->res_name,
                             atom1_ptr->residue_ptr->res_num,
                             atom1_ptr->residue_ptr->chain,
                             atom1_ptr->residue_ptr->pdbentry_ptr->pdb_code,
                             best_atom_ptr->atom_name,
                             best_atom_ptr->residue_ptr->res_name,
                             best_atom_ptr->residue_ptr->res_num,
                             best_atom_ptr->residue_ptr->chain,
                             best_atom_ptr->residue_ptr->
                             pdbentry_ptr->pdb_code);
/* <-- v.1.3.7 */

/* v.1.3 --> */
//                      /* Determine whether the two atoms both have H-bonds
//                         to equivalenced residues */
//                      equiv = check_equiv_hbonds(atom1_ptr,best_atom_ptr);
//                      
//                      /* If so, then increase anchor weight between them */
//                      if (equiv == TRUE)
//                        atomap_ptr->anchor_weight = ANCHOR_WEIGHT;
/* <-- v.1.3 */

                      /* Increment count of saved atoms */
                      natoms++;
                    }

                  /* If found, then update if current mapping has a higher
                     score */
                  else
                    {
                      /* If current mapping has a smaller distance, replace
                         old mapping */
                      if (best_match < atomap_ptr->score)
                        {
                          atomap_ptr->atom_ptr = atom1_ptr;
                          atomap_ptr->score = best_match;
                        }
                    }
                }
            }
        }

      /* Get the next atom record and increment atom count */
      atom1_ptr = atom1_ptr->next_atom_ptr;
      iatom++;
    }

  /* Return first of the atom matches */
  *fst_atomap_ptr = first_atomap_ptr;
  *lst_atomap_ptr = last_atomap_ptr;

  /* Return number of atom mappings stored */
  *natm = natoms;

  /* Return number of overlapping atoms */
  return(noverlap);
}
/***********************************************************************

get_file_name  -  Open give file in current or working directory

***********************************************************************/

void get_file_name(char *directory,char *name,char *file_name)
{
/* v.2.2.8 --> */
  int spaces;
/* <-- v.2.2.8 */


  /* Form the name of the output file */
  if (directory != &null)
    {
/* v.2.2.8 --> */
      /* See if directory contains any spaces */
      if (strchr(directory,' ') != NULL)
        spaces = TRUE;
      else
        spaces = FALSE;
      spaces = FALSE;
/* <-- v.2.2.8 */
      
      /* Form the file name in the working directory */
/* v.2.2.8 --> */
//      strcpy(file_name,directory);
      if (spaces == TRUE)
        {
          strcpy(file_name,"\"");
          strcat(file_name,directory);
        }
      else
        strcpy(file_name,directory);
/* <-- v.2.2.8 */
      strcat(file_name,DIR_SEP);
      strcat(file_name,name);
/* v.2.2.8 --> */
      if (spaces == TRUE)
        strcat(file_name,"\"");
/* <-- v.2.2.8 */
    }
  else
    strcat(file_name,name);
}
/***********************************************************************

open_file  -  Open given file in current or working directory

***********************************************************************/

FILE *open_file(char *name,struct parameters params)
{
  char file_name[FILENAME_LEN];

  FILE *file_ptr;

  /* Form the name of the output file */
  get_file_name(params.directory,name,file_name);

  /* Open the output file */
  open_output_file(&file_ptr,file_name);

  /* Return the file pointer */
  return file_ptr;
}
/***********************************************************************

map_extras  -  Map any unmapped residues if they significantly overlap
               any LIGPLOT residues

***********************************************************************/

int map_extras(struct residue *fit_first_residue_ptr,
               struct residue *first_residue_ptr,
               REAL matrix[3][3],double cofg[2][3],int npaired,
/* v.1.3.7 --> */
//               struct atomap **fst_atomap_ptr,int fit_type)
               struct atomap **fst_atomap_ptr,int fit_type,
               struct parameters params)
/* <-- v.1.3.7 */
{
  int natomap, noverlap, max_overlap;
  int store;
/* v.2.0 --> */
  int type1, type2;
/* <-- v.2.0 */

  struct atomap *first_atomap_ptr, *last_atomap_ptr, *start_atomap_ptr;
  struct residue *residue1_ptr, *residue2_ptr, *overlap_residue_ptr;

/* v.1.3.7 --> */
  FILE *file_ptr;
/* <-- v.1.3.7 */

  /* Initialise variables */
  start_atomap_ptr = *fst_atomap_ptr;
  first_atomap_ptr = last_atomap_ptr = NULL;
  natomap = 0;

/* v.1.3.7 --> */
  /* Open the output ligand.map file */
  file_ptr = open_file("ligand.map",params);
/* <-- v.1.3.7 */

  /* Get the first of the residues */
  residue1_ptr = fit_first_residue_ptr;

  /* Loop to find all unmapped residues */
  while (residue1_ptr != NULL)
    {
      /* If residue unmapped, compare its 3D location with all other
         unmapped residues in the LIGPLOT diagram */
      if (residue1_ptr->map_residue_ptr == NULL)
        {
          /* Initialise variables */
          max_overlap = 0;
          overlap_residue_ptr = NULL;

          /* Get the first of the LIGPLOT residues */
          residue2_ptr = first_residue_ptr;

          /* Loop to find any unmapped residues */
          while (residue2_ptr != NULL)
            {
/* v.2.0 --> */
              /* Get the two residue types */
              if (params.dimplot == TRUE &&
                  residue1_ptr->chain_type == CHAIN2)
                type1 = LIGAND;
              else
                type1 = residue1_ptr->type;
              if (params.dimplot == TRUE &&
                  residue2_ptr->chain_type == CHAIN2)
                type2 = LIGAND;
              else
                type2 = residue2_ptr->type;
/* <-- v.2.0 */

             /* If residue unmapped, see if it overlaps our residue */
              if (residue2_ptr->map_residue_ptr == NULL &&
/* v.2.0 --> */
//                  ((residue1_ptr->type != LIGAND && 
//                    residue2_ptr->type != LIGAND) || 
//                   (residue1_ptr->type == LIGAND && 
//                    residue2_ptr->type == LIGAND)))
                  ((type1 != LIGAND && type2 != LIGAND) || 
                   (type1 == LIGAND && type2 == LIGAND)))
/* <-- v.2.0 */
                {
                  /* Set flag so that atom matches aren't stored */
                  store = FALSE;

                  /* Check for level of overlap */
                  noverlap
                    = residue_overlap(residue1_ptr,residue2_ptr,matrix,
                                      cofg,&first_atomap_ptr,
                                      &last_atomap_ptr,&natomap,
/* v.1.3.7 --> */
//                                      store);
                                      store,file_ptr,params);
/* <-- v.1.3.7 */

                  /* If this is the best overlap so far, then store
                     overlapping residue */
                  if (noverlap > max_overlap)
                    {
                      /* Save the number of overlapping atoms */
                      max_overlap = noverlap;

                      /* Save the residue pointer */
                      overlap_residue_ptr = residue2_ptr;
                    }
                }

              /* Get the next residue */
              residue2_ptr = residue2_ptr->next_residue_ptr;
            }

          /* If have an overlapping residue, check whether the mapping
             is a sensible one */
          if (overlap_residue_ptr != NULL &&
              (max_overlap > 3 || (residue1_ptr->natoms == 1 &&
                                   overlap_residue_ptr->natoms == 1)))
            {
/* v.2.0 --> */
              /* Get the overlap residue's type */
              if (params.dimplot == TRUE &&
                  overlap_residue_ptr->chain_type == CHAIN2)
                type2 = LIGAND;
              else
                type2 = overlap_residue_ptr->type;
/* <-- v.2.0 */

              /* If one residue is the ligand, but other isn't, then
                 reject this mapping */
/* v.2.0 --> */
//              if ((residue1_ptr->type == LIGAND &&
//                   overlap_residue_ptr->type != LIGAND) ||
//                  (residue1_ptr->type != LIGAND &&
//                   overlap_residue_ptr->type == LIGAND))
              if ((type1 == LIGAND && type2 != LIGAND) ||
                  (type1 != LIGAND && type2 == LIGAND))
/* <-- v.2.0 */
                overlap_residue_ptr = NULL;

              /* Otherwise, if both are non-ligand residues, check for
                 equivalence of type */
/* v.2.0 --> */
//              else if (residue1_ptr->type != LIGAND &&
//                       overlap_residue_ptr->type != LIGAND)
              else if (type1 != LIGAND && type2 != LIGAND)
/* <-- v.2.0 */
                {
                  /* If one is a water, and the other is not, then not
                     a good match */
/* v.2.0 --> */
//                  if ((residue1_ptr->type == WATER &&
//                       overlap_residue_ptr->type != WATER) ||
//                      (residue1_ptr->type != WATER &&
//                       overlap_residue_ptr->type == WATER))
                  if ((type1 == WATER && type2 != WATER) ||
                      (type1 != WATER && type2 == WATER))
/* <-- v.2.0 */
                    overlap_residue_ptr = NULL;

                  /* Ditto for metal */
/* v.2.0 --> */
//                  else if ((residue1_ptr->type == METAL_ION &&
//                            overlap_residue_ptr->type != METAL_ION) ||
//                           (residue1_ptr->type != METAL_ION &&
//                            overlap_residue_ptr->type == METAL_ION))
                  else if ((type1 == METAL_ION && type2 != METAL_ION) ||
                           (type1 != METAL_ION && type2 == METAL_ION))
/* <-- v.2.0 */
                    overlap_residue_ptr = NULL;
                }

              /* If both are ligand residues, get the best atom-atom
                 mappings from the overlay */
/* v.2.0 --> */
//              else if (residue1_ptr->type == LIGAND &&
//                       overlap_residue_ptr->type == LIGAND)
              else if (type1 == LIGAND && type2 == LIGAND)
/* <-- v.2.0 */
                {
                  /* Set flag so that atom matches are stored */
                  store = TRUE;

                  /* Get the atom mappings */
                  noverlap
                    = residue_overlap(residue1_ptr,overlap_residue_ptr,
                                      matrix,cofg,&first_atomap_ptr,
                                      &last_atomap_ptr,&natomap,
/* v.1.3.7 --> */
//                                      store);
                                      store,file_ptr,params);
/* <-- v.1.3.7 */
                }              
            }
          else
            overlap_residue_ptr = NULL;

          /* If overlap looks good, store the mapping */
          if (overlap_residue_ptr != NULL)
            {
              /* Map the two residues to one another */
              residue1_ptr->map_residue_ptr = overlap_residue_ptr;
              overlap_residue_ptr->map_residue_ptr = residue1_ptr;

              /* Increment count of paired residues */
              npaired++;
            }
        }

      /* Get the next residue */
      residue1_ptr = residue1_ptr->next_residue_ptr;
    }

  /* Return first of the ligand atom matches (if didn't have one before) */
  if (start_atomap_ptr == NULL)
    start_atomap_ptr = first_atomap_ptr;
  *fst_atomap_ptr = start_atomap_ptr;

/* v.1.3.7 --> */
  /* Close the output file */
  fclose(file_ptr);
/* <-- v.1.3.7 */

  /* Return number of paired residues */
  return(npaired);
}
/* v.1.2 --> */
/* v.1.3 --> */
/***********************************************************************

store_alignment  -  Store the alignment as an array of paired resseq
                    records

***********************************************************************/

int store_alignment(struct resseq **first_resseq_ptr,
/* v.2.0.1 --> */
//                    struct resseq ***aln,int nfit)
                    struct resseq ***aln,int nfit,
                    struct parameters params)
/* <-- v.2.0.1 */
{
  int align_len, ipos;
/* v.2.0.1 --> */
//  int verbose;
/* <-- v.2.0.1 */

  struct resseq *resseq_ptr[2];

  struct resseq **alignment;

  /* Initialise variables */
/* v.2.0.1 --> */
//  verbose = FALSE;
/* <-- v.2.0.1 */

  /* Create the alignment array */
  align_len = nfit;
  alignment = (struct resseq **)
    malloc(align_len * 2 * sizeof(struct resseq *));
  if (alignment == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for alignment array ");
      exit(1);
    }

  /* Initialise alignment array with NULLs */
  for (ipos = 0; ipos < 2 * align_len; ipos++)
    alignment[ipos] = NULL;

  /* Get the first of the stored residues */
  resseq_ptr[0] = first_resseq_ptr[0];
  resseq_ptr[1] = first_resseq_ptr[1];
  ipos = 0;

  /* Loop over all the aligned residues and place into the alignment
     array */
  while (resseq_ptr[0] != NULL && resseq_ptr[1] != NULL && ipos < nfit)
    {
      /* Store in alignment */
      alignment[ipos] = resseq_ptr[0];
      alignment[align_len + ipos] = resseq_ptr[1];

      /* If in verbose mode, show alignment at this position */
/* v.2.0.1 --> */
//      if (verbose == TRUE)
      if (params.verbose == TRUE)
/* <-- v.2.0.1 */
        {
          printf("ALIGN:  ");
          if (resseq_ptr[0] != NULL &&
              resseq_ptr[0]->residue_ptr != NULL)
            printf(" %s %s %c from %s ",
                   resseq_ptr[0]->residue_ptr->res_name,
                   resseq_ptr[0]->residue_ptr->res_num,
                   resseq_ptr[0]->residue_ptr->chain,
                   resseq_ptr[0]->residue_ptr->pdbentry_ptr->pdb_code);
          else
            printf(" --- ----- - .... ---- ");
          if (resseq_ptr[1] != NULL &&
              resseq_ptr[1]->residue_ptr != NULL)
            printf(" %s %s %c from %s ",
                   resseq_ptr[1]->residue_ptr->res_name,
                   resseq_ptr[1]->residue_ptr->res_num,
                   resseq_ptr[1]->residue_ptr->chain,
                   resseq_ptr[1]->residue_ptr->pdbentry_ptr->pdb_code);
          else
            printf(" --- ----- - .... ---- ");
          printf("\n");
        }

      /* Get the next residues */
      resseq_ptr[0] = resseq_ptr[0]->next_resseq_ptr;
      resseq_ptr[1] = resseq_ptr[1]->next_resseq_ptr;
      ipos++;
    }

  /* Return the alignment */
  *aln = alignment;

  /* Return the alignment length */
  return(align_len);
}
/* <-- v.1.3 */
/***********************************************************************

get_cora_alignment  -  Read in the CORA-format alignment file

***********************************************************************/

int get_cora_alignment(char *file_name,
                       struct pdbentry *first_pdbentry_ptr,
/* v.2.0.1 --> */
//                       struct resseq ***aln)
                       struct resseq ***aln,struct parameters params)
/* <-- v.2.0.1 */
{
  char input_line[LINELEN + 1], pdb_code[5];
  char number_string[10], residue_details[12];
  char aa_code, chain, res_name[4], res_num[6], res_number[2][6];
  char token[MAX_TOKEN][TOKEN_LEN];

  int column, ipdb, ipos, iresid, iseq, line;
  int nchains, nfit, nmatch, npaired, nresid[2], nresidues, nseqs;
  int aln_pos[2], column_width, start_pos, version;
  int a, aa_num, z, i, itoken, ntokens;
  int fit;

  struct pdbentry *pdbentry_ptr[2];
  struct residue *residue_ptr, *match_residue_ptr[2];
  struct resseq *resseq_ptr[2], *first_resseq_ptr[2], *last_resseq_ptr[2];

  struct resseq **alignment;

  FILE *file_ptr;

  static char chain_list[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890";

  /* Initialise variables */
  alignment = NULL;
  a = 'a';
  z = 'z';
  aln_pos[0] = aln_pos[1] = -1;
  iresid = 0;
  line = 0;
  nchains = strlen(chain_list);
  nseqs = 0;
  npaired = 0;
  nresidues = nresid[0] = nresid[1] = 0;
  number_string[0] = '\0';
  nfit = 0;

  /* Initialise pointers */
  pdbentry_ptr[0] = pdbentry_ptr[1] = NULL;
  pdbentry_ptr[CURR_SEQUENCE] = first_pdbentry_ptr;
  if (first_pdbentry_ptr != NULL)
    pdbentry_ptr[REF_SEQUENCE] = first_pdbentry_ptr->next_pdbentry_ptr;
  if (pdbentry_ptr[0] == NULL || pdbentry_ptr[1] == NULL)
    return(nfit);
  first_resseq_ptr[0] = first_resseq_ptr[1] = NULL;
  last_resseq_ptr[0] = last_resseq_ptr[1] = NULL;

  /* Initialise CORA-format parameters */
  column_width = 10;
  start_pos = 13;
  version = OLD_CORA;

  /* Open the alignments file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    {
      printf("\n*** Warning. Unable to open alignments file [%s]\n",
             file_name);
      return(nfit);
    }

  /* Loop through the file, line by line */
  while (fgets(input_line,LINELEN,file_ptr) != NULL)
    {
      /* Truncate the string to strip off the newline */
      string_truncate(input_line,LINELEN);

      /* Check for new CORA version */
      if (!strncmp(input_line,"#FM CORA_FORMAT",15))
        {
          /* Set CORA-format parameters */
          column_width = 11;
          start_pos = 14;
          version = NEW_CORA;
        }

      /* Process only if this is not a comment line */
      if (input_line[0] != '#')
        {
          /* If this is the first record, read in number of sequences */
          if (line == 0)
            nseqs = atoi(input_line);

          /* If this is the second record, read in the sequence names and
             identify columns  corresponding to our PDB codes of interest */
          else if (line == 1)
            {
              /* Extract all the tokens on the input line */
              ntokens = get_tokens(input_line,token,' ');

              /* Loop over the tokens to get the PDB codes */
              for (itoken = 0; itoken < ntokens; itoken++)
                {
                  /* Extract the PDB code from this token */
                  if (strlen(token[itoken]) > 3)
                    {
                      /* Extract the wolf code containing PDB code,
                         chain ID and domain */
                      chain = ' ';
                      strncpy(pdb_code,token[itoken],4);
                      pdb_code[4] = '\0';
                      if (strlen(token[itoken]) > 4)
                        chain = token[itoken][4];

                      /* Convert to lower case */
                      to_lower(pdb_code,4);

                      /* Check to see if PDB code matches either of
                         our codes of interest */
                      for (ipdb = 0; ipdb < 2; ipdb++)
                        {
                          /* If this the right PDB code, then store 
                             sequence position */
                          if (aln_pos[ipdb] == -1 &&
                              !strcmp(pdb_code,pdbentry_ptr[ipdb]->pdb_code))
                            {
                              aln_pos[ipdb] = itoken;
                              pdbentry_ptr[ipdb]->chain = chain;
                            }
                        }
                    }
                }

              /* If either PDB code not found, then abandon */
              if (aln_pos[0] == -1 || aln_pos[1] == -1)
                return(nfit);
            }

          /* If this is the third record, read in the total number of 
             residues in the alignment */
          else if (line == 2)
            nresidues = atoi(input_line);

          /* Otherwise, get the alignment at this residue position */
          else if (line > 2 && iresid < nresidues)
            {
              /* Initialise variables for this alignment line */
              nmatch = 0;
              fit = FALSE;
              match_residue_ptr[0] = match_residue_ptr[1] = NULL;

              /* Loop over our two sequences to see if have equivalenced
                 residues at this position */
              for (iseq = 0; iseq < 2; iseq++)
                {
                  /* Retrieve the alignment column */
                  column = aln_pos[iseq];

                  /* Get the chain identifier */
                  chain = pdbentry_ptr[iseq]->chain;

                  /* Calculate start of this residue's details */
                  ipos = start_pos + column * column_width;

                  /* Extract the residue details for this residue */
                  strncpy(residue_details,input_line+ipos,column_width);
                  residue_details[column_width] = '\0';

                  /* Get the residue name, number and chain-id from the
                     residue details */
                  if (version == OLD_CORA)
                    aa_code = residue_details[6];
                  else
                    aa_code = residue_details[7];

                  /* Check for lower-case letters (which represent
                     cysteines involved in disulphide bonds) */
                  aa_num = aa_code;
                  if (aa_num >= a && aa_num <= z)
                    aa_code = 'C';

                  /* Get the corresponding three-letter code */
                  get_aa_name(aa_code,res_name);

                  /* Get the residue number and chain id */
                  if (version == OLD_CORA)
                    strncpy(res_num,residue_details,5);
                  else
                    strncpy(res_num,residue_details+1,5);
                  res_num[5] = '\0';
                  strcpy(res_number[iseq],res_num);

                  /* Identify this residue */
                  residue_ptr = pdbentry_ptr[iseq]->first_residue_ptr;

                  /* Loop over this PDB entry's residues to find current
                     one */
                  while (residue_ptr != NULL &&
                         match_residue_ptr[iseq] == NULL)
                    {
                      /* If residue matches, then save pointer */
                      if (!strcmp(residue_ptr->res_num,res_num) &&
                          !strcmp(residue_ptr->res_name,res_name) &&
                          (chain == ' ' || chain == residue_ptr->chain))
                        match_residue_ptr[iseq] = residue_ptr;

                      /* Get pointer to next residue record */
                      residue_ptr = residue_ptr->next_residue_ptr;
                    }
                }

              /* If both residues identified, then store the equivalence */
              if (match_residue_ptr[0] != NULL &&
                  match_residue_ptr[1] != NULL)
                {
                  /* Create a resseq record for both residues */
                  for (i = 0; i < 2; i++)
                    {
                      /* Create record */
                      resseq_ptr[i]
                        = create_resseq_record(&(first_resseq_ptr[i]),
                                               &(last_resseq_ptr[i]),
                                               match_residue_ptr[i]);

                      /* Increment count of store residues */
                      nresid[i]++;
                    }

                  /* Store pointers between the two resseq records */
                  resseq_ptr[0]->align_resseq_ptr = resseq_ptr[1];
                  resseq_ptr[1]->align_resseq_ptr = resseq_ptr[0];

                  /* Save the equivalence */
                  //                  match_residue_ptr[0]->map_residue_ptr
                  //  = match_residue_ptr[1];
                  //match_residue_ptr[1]->map_residue_ptr
                  //  = match_residue_ptr[0];

                  /* Increment count of equivalenced residue */
                  nfit++;
                }

              /* Increment residue count */
              iresid++;
            }

          /* Increment line count */
          line++;
        }
    }

  /* Close the alignments file */
  fclose(file_ptr);

  /* Show number of matched residue pairs that can be fitted */
  printf("Number of residues matched from alignment: %d\n",nfit);

  /* If have enough equivalenced residues, store them in an array for use
     in the fitting routines */
  if (nfit > 2)
/* v.2.0.1 --> */
//    store_alignment(first_resseq_ptr,&alignment,nfit);
    store_alignment(first_resseq_ptr,&alignment,nfit,params);
/* <-- v.2.0.1 */

  /* Return the alignment */
  *aln = alignment;

  /* Return number of residue equivalences found */
  return(nfit);
}
/* <-- v.1.2 */
/***********************************************************************

get_residue_names  -  Get the residue-similarity scores

***********************************************************************/

void get_residue_names(char residue_name[NAMINO+1][4])
{
  char residue_code;

  int ipos, ires;

  static char *aa_codes =
    "  X  A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V";

  /* Loop through all the residues to get residue names */
  ipos = 0;
  for (ires = 0; ires < NAMINO + 1; ires++)
    {
      /* Get the residue code */
      residue_code = aa_codes[ipos + 2];
      ipos = ipos + 3;

      /* Get the corresponding 3-letter code */
      get_aa_name(residue_code,residue_name[ires]);
    }
}
/* v.1.3 --> */
/***********************************************************************

get_sequence  -  Align the two protein sequences to get residue
                 equivalences

***********************************************************************/

int get_sequence(struct pdbentry *pdbentry_ptr,char **seq,
/* v.1.3.7 --> */
//                 struct residue ***residue_idx_ptr)
                 struct residue ***residue_idx_ptr,char wanted_chain)
/* <-- v.1.3.7 */
{
/* v.1.3.7 --> */
  // char chain, last_chain;
  char chain, last_chain, map_chain;
/* <-- v.1.3.7 */

/* v.1.3.6 --> */
  // char *sequence;
  char *sequence, *tmp_sequence;
/* <-- v.1.3.6 */

  int chain_end, chain_start, i, ipos, iresid, max_int, nint, nresid;
  int seq_start, seq_end;

  struct residue *residue_ptr;

  struct residue **residue_index_ptr;
  
  /* Initialise variables */
  chain_start = chain_end = seq_start = seq_end = -1;
  ipos = 0;
  iresid = 0;
  last_chain = chain = '\0';
/* v.1.3.7 --> */
  map_chain = '\0';
/* <-- v.1.3.7 */
  max_int = nint = 0;
  nresid = pdbentry_ptr->nresid;

  /* Create an array for storing the protein sequence */
  sequence = (char *) malloc(sizeof(char)*(nresid + 1));
/* v.1.3.6 --> */
  tmp_sequence = (char *) malloc(sizeof(char)*(nresid + 1));
/* <-- v.1.3.6 */
  
  /* Create an array for storing the residue pointers */
  residue_index_ptr
    = (struct residue **) malloc(nresid * sizeof(struct residue *));
  if (residue_index_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for residue index\n");
      exit(1);
    }

  /* Initialise index array */
  for (i = 0; i < nresid; i++)
    residue_index_ptr[i] = NULL;
  
  /* Get the full entry's first residue */
  residue_ptr = pdbentry_ptr->first_residue_ptr;

  /* Loop over the residues to find corresponding records in entry read
     in from the ensemble.drw file */
  while (residue_ptr != NULL)
    {
      /* Check that we're not going over the array limit */
      if (ipos > nresid)
        {
          printf("*** PROGRAM ERROR. Number of residues in %s exceeded: "
                 "%d (nresid = %d)\n",pdbentry_ptr->pdb_code,
                 iresid,nresid);
          exit(1);
        }

      /* If an amino acid, add to sequence */
      if (residue_ptr->aa_code != 'X' && residue_ptr->aa_code != '?' &&
          residue_ptr->interaction != IS_LIGAND)
        {
          /* Get the chain */
          chain = residue_ptr->chain;

          /* If new chain, then see if previous chain had the most
             interactions */
          if (chain != last_chain)
            {
              /* If last chain had the most interactions with ligand
                 so far, then store it as the representative chain
                 to use in the sequence alignments */
/* v.1.3.7 --> */
//              if (nint > max_int)
              if (last_chain == wanted_chain ||
                  (wanted_chain == '\0' && nint > max_int))
/* <-- v.1.3.7 */
                {
                  /* Save maximum number of interactions */
                  max_int = nint;

                  /* Save the last chain as the best to use for
                     alignment */
                  seq_start = chain_start;
                  seq_end = chain_end;

/* v.1.3.7 --> */
                  /* Save the chain to be used for the alignment */
                  map_chain = last_chain;
/* <-- v.1.3.7 */
                }

              /* Save current residue position as start of new chain */
              chain_start = ipos;

              /* Re-initialise count of interacting residues */
              nint = 0;

              /* Save the chain */
              last_chain = chain;
            }

          /* If residue is an interacting one, then increment count of
             interactions for current chain */
          if (residue_ptr->interaction != NO_INTERACTION)
            nint++;

          /* Save the residue code */
          sequence[ipos] = residue_ptr->aa_code;

          /* Store the corresponding residue pointer */
          residue_index_ptr[ipos] = residue_ptr;

          /* Save current chain end */
          chain_end = ipos;

          /* Increment sequence count */
          ipos++;
        }

      /* Get pointer to next residue record */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;
    }

  /* Terminate the sequence */
  sequence[ipos] = '\0';

  /* Store last chain if it had the most interactions */
/* v.1.3.7 --> */
//  if (nint > max_int)
  if (last_chain == wanted_chain ||
      (wanted_chain == '\0' && nint > max_int))
/* <-- v.1.3.7 */
    {
      /* Save the last chain as the best to use for alignment */
      seq_start = chain_start;
      seq_end = chain_end;

/* v.1.3.7 --> */
      /* Save the chain to be used for the alignment */
      map_chain = last_chain;
/* <-- v.1.3.7 */
    }

/* v.1.3.7 --> */
  /* Save the chain used for the mapping */
  pdbentry_ptr->map_chain = map_chain;
/* <-- v.1.3.7 */

  /* If have identified the best sequence to use, then return it */
  if (seq_start > -1)
    {
      /* Top and tail the sequence */
      sequence[seq_end + 1] = '\0';
      if (seq_start > 0)
/* v.1.3.6 --> */
        // strcpy(sequence,sequence + seq_start);
        {
          strcpy(tmp_sequence,sequence + seq_start);
          strcpy(sequence,tmp_sequence);
        }
/* <-- v.1.3.6 */

      /* Return the sequence */
      *seq = sequence;

      /* Return the array of residue pointers */
      *residue_idx_ptr = residue_index_ptr + seq_start;

      /* Get its length */
      nresid = seq_end - seq_start + 1;
    }
  else
    nresid = 0;

  /* Return the number of residue in the stored sequence */
  return(nresid);
}
/* v.2.0 --> */
/***********************************************************************

get_dimplot_sequence  -  Retrieve the sequence of the first chain/domain

***********************************************************************/

int get_dimplot_sequence(struct pdbentry *pdbentry_ptr,char **seq,
                         struct residue ***residue_idx_ptr,int in_plot)
{
  char chain, map_chain;

  char *sequence;

  int i, ipos, iresid, nresid;

  struct residue *residue_ptr;

  struct residue **residue_index_ptr;
  
  /* Initialise variables */
  ipos = 0;
  iresid = nresid = 0;
  map_chain = '\0';

  /* Get the PDB entry's first residue */
  residue_ptr = pdbentry_ptr->first_residue_ptr;

  /* Loop over the residues to count how many belong to the first
     chain/domain  */
  while (residue_ptr != NULL)
    {
      /* If this is a chain1 residue, increment count */
      if (residue_ptr->chain_type == CHAIN1 &&
          (in_plot == FALSE ||
           residue_ptr->interaction != NO_INTERACTION))
        nresid++;

      /* Get pointer to next residue record */
      residue_ptr = residue_ptr->next_residue_ptr;
    }

  /* If no residues found, then can't continue */
  if (nresid == 0)
    {
      printf("*** ERROR. No residues found from first chain/domain\n");
      exit(1);
    }

  /* Create an array for storing the protein sequence */
  sequence = (char *) malloc(sizeof(char)*(nresid + 1));
  
  /* Create an array for storing the residue pointers */
  residue_index_ptr
    = (struct residue **) malloc(nresid * sizeof(struct residue *));
  if (residue_index_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for residue index\n");
      exit(1);
    }

  /* Initialise index array */
  for (i = 0; i < nresid; i++)
    residue_index_ptr[i] = NULL;

  /* Get the full entry's first residue */
  residue_ptr = pdbentry_ptr->first_residue_ptr;

  /* Loop over the residues to find corresponding records in entry read
     in from the ensemble.drw file */
  while (residue_ptr != NULL)
    {
      /* Check that we're not going over the array limit */
      if (ipos > nresid)
        {
          printf("*** PROGRAM ERROR. Number of residues in %s exceeded: "
                 "%d (nresid = %d)\n",pdbentry_ptr->pdb_code,
                 iresid,nresid);
          exit(1);
        }

      /* If residue belongs to the first chain/domain, add to sequence */
      if (residue_ptr->chain_type == CHAIN1 &&
          (in_plot == FALSE ||
           residue_ptr->interaction != NO_INTERACTION))
        {
          /* Get the chain */
          chain = residue_ptr->chain;

          /* Save the chain to be used for the alignment */
          if (map_chain == '\0')
            map_chain = chain;

          /* Save the residue code */
          sequence[ipos] = residue_ptr->aa_code;

          /* Store the corresponding residue pointer */
          residue_index_ptr[ipos] = residue_ptr;

          /* Increment sequence count */
          ipos++;
        }

      /* Get pointer to next residue record */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;
    }

  /* Terminate the sequence */
  sequence[ipos] = '\0';

  /* Save the chain used for the mapping */
  pdbentry_ptr->map_chain = map_chain;

  /* Return the sequence */
  *seq = sequence;

  /* Return the array of residue pointers */
  *residue_idx_ptr = residue_index_ptr;

  /* Return the number of residue in the stored sequence */
  return(nresid);
}
/* <-- v.2.0 */
/***********************************************************************

create_alignment_arrays  -  Create the arrays for the sequence alignment

***********************************************************************/

int create_alignment_arrays(int nresid[2],int **arr,int **tmparr,
                            int **pth,int **maparr)
{
  int max_path, narray, nmap;

  int *array, *mapping, *path, *tmp_array;

  /* Initialise variables */
  narray = nresid[0] * nresid[1];
  nmap = nresid[1];
  max_path = nresid[0] + nresid[1];
 
  /* Create Needleman-Wunsch array for performing the alignment */
  array = (int *) malloc(narray * sizeof(int));
  if (array == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for NW array\n");
      exit(1);
    }

  /* Create an array for storing the path through the matrix that
     gives the best alignment score */
  path = (int *) malloc(max_path * sizeof(int));
  if (path == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for path\n");
      exit(1);
    }

  /* Create working array for use in Needleman-Wunsch process */
  tmp_array = (int *) malloc(narray * sizeof(int));
  if (tmp_array == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for temporary array\n");
      exit(1);
    }

  /* Create array mapping residues across the two sequences */
  mapping = (int *) malloc(nmap * sizeof(int));
  if (mapping == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for mapping array\n");
      exit(1);
    }

  /* Return the array pointers */
  *arr = array;
  *pth = path;
  *tmparr = tmp_array;
  *maparr = mapping;

  /* Return size of mapping array */
  return(nmap);
}
/***********************************************************************

fill_array  -  Initialise array with standard residue similarities

***********************************************************************/

void fill_array(char *sequence[2],int *nresid,int *array,int *tmp_array,
/* v.2.0 --> */
//                int *mapping)
                int *mapping,struct residue **residue_index_ptr[2])
/* <-- v.2.0 */
{
#define MIN_VALUE    0

  char aa_code1, aa_code2;

  char *string_ptr;

  int i, ipoint, ires, jres, max_res, res_type1, res_type2;

/* v.2.0 --> */
  struct residue *residue1_ptr, *residue2_ptr;
/* <-- v.2.0 */

  /* Initialise variables */
  max_res = nresid[0];
  if (nresid[1] > max_res)
    max_res = nresid[1];

  /* Initialise both sequence's mappings */
  for (i = 0; i < nresid[1]; i++)
    mapping[i] = -1;

  /* Loop over the residues in first sequence */
  for (ires = 0; ires < nresid[0]; ires++)
    {
      /* Get the corresponding residue */
      aa_code1 = sequence[0][ires];
      res_type1 = 0;
      string_ptr = strchr(AA_CODES,aa_code1);
      if (string_ptr != NULL)
        res_type1 = string_ptr - AA_CODES;
/* v.2.0 --> */
      residue1_ptr = residue_index_ptr[0][ires];
/* <-- v.2.0 */

      /* Loop over the residues in the second sequence */
      for (jres = 0; jres < nresid[1]; jres++)
        {
          /* Get the corresponding residue */
          aa_code2 = sequence[1][jres];
          res_type2 = 0;
          string_ptr = strchr(AA_CODES,aa_code2);
          if (string_ptr != NULL)
            res_type2 = string_ptr - AA_CODES;
/* v.2.0 --> */
          residue2_ptr = residue_index_ptr[1][jres];
/* <-- v.2.0 */

          /* Get the array location of this residue pair */
          ipoint = ires * nresid[1] + jres;

          /* Get the similarity between these two residues */
          array[ipoint]
            = SIMILARITY_MATRIX[res_type1][res_type2] - MIN_VALUE;

/* v.2.0 --> */
          /* If the residue types are identical, add extra to score
             if residue numbers are the same */
          if (res_type1 == res_type2 && residue1_ptr != NULL &&
              residue2_ptr != NULL &&
              residue1_ptr->aa_code == residue2_ptr->aa_code &&
              residue1_ptr->chain == residue2_ptr->chain &&
              !strcmp(residue1_ptr->res_num,residue2_ptr->res_num))
            {
              /* Add to the similarity score */
              array[ipoint] = array[ipoint] + IDENTITY_SCORE;
            }
/* <-- v.2.0 */

          /* Store in temporary array */
          tmp_array[ipoint] = array[ipoint];
        }
    }
}
/***********************************************************************

nw_sweep  -  Sweep through the alignment matrix, updating the cells
             according to the best path from any cell on the previous
             row or column - slow version

***********************************************************************/

void nw_sweep(int *array,int nresid[2])
{
  int ipoint, ires, jpoint, jres, kres, lres;
  int max_score, score;

  /* Loop backwards through the array, starting at the bottom right
     cell and working leftwards and upwards */
  for (ires = nresid[0] - 2; ires > -1; ires--)
    { 
      for (jres = nresid[1] - 2; jres > -1; jres--)
        {
          /* Get the array location of this residue pair */
          ipoint = ires * nresid[1] + jres;

          /* Initialise maximum score */
          max_score = -100;

          /* Loop over the cells to the right in the row below the
             current one */
          for (kres = ires + 1, lres = jres + 1;
               kres < nresid[0] && lres < nresid[1]; lres++)
            {
              /* Get this cell's array location */
              jpoint = kres * nresid[1] + lres;

              /* Calculate score to get from here to our current cell */
              score = array[jpoint];
              if (lres > jres + 1)
                score = score - GAP_PENALTY;

              /* If this is the largest score so far, then store */
              if (score > max_score)
                max_score = score;
            }

          /* Loop over the downward cells in the column to the right
             of the current one */
          for (kres = ires + 1, lres = jres + 1;
               kres < nresid[0] && lres < nresid[1]; kres++)
            {
              /* Get this cell's array location */
              jpoint = kres * nresid[1] + lres;

              /* Calculate score to get from here to our current cell */
              score = array[jpoint];
              if (kres > ires + 1)
                score = score - GAP_PENALTY;

              /* If this is the largest score so far, then store */
              if (score > max_score)
                max_score = score;
            }

          /* Add maximum score to cell's current score */
          array[ipoint] = array[ipoint] + max_score;
        }
    }
}
/***********************************************************************

find_best_path  -  Find the best path through the matrix, corresponding
                   to the best alignment of the two sequences

***********************************************************************/

int find_best_path(int *array,int nresid[2],int *path,int *alnlen)
{
  int ipoint, ires, jres, start_ires, start_jres;
  int align_len, last_pos[2], max_path, max_pos[2], path_len;
  int best_dist, dist, gap_penalty, max_score, score;
  int done, first;

  /* Initialise variables */
  align_len = 0;
  best_dist = 999;
  path_len = 0;
  start_ires = 0;
  start_jres = 0;
  done = FALSE;
  first = TRUE;
  last_pos[0] = -1;
  last_pos[1] = -1;
  max_path = nresid[0] + nresid[1];

  /* Loop until entire array has been traversed */
  while (done == FALSE)
    {
      /* Initialise maximum score and its position */
      gap_penalty = 0;
      max_score = -9999999;
      max_pos[0] = -1;
      max_pos[1] = -1;

      /* Loop through row to get maximum score and position */
      ires = start_ires;
      for (jres = start_jres; jres < nresid[1]; jres++)
        {
          /* Get the array location of this cell */
          ipoint = ires * nresid[1] + jres;
          dist = jres - start_jres;
          score = array[ipoint] - gap_penalty;

          /* If this is the highest score so far, store score and
             position */
          if (score > max_score)
            {
              max_score = score;
              max_pos[0] = ires;
              max_pos[1] = jres;
              best_dist = dist;
            }

          /* Set the gap penalty */
          if (first == FALSE)
            gap_penalty = GAP_PENALTY;
        }

      /* Re-initialise gap penalty */
      gap_penalty = 0;

      /* Loop column to get maximum score and position */
      jres = start_jres;
      for (ires = start_ires; ires < nresid[0]; ires++)
        {
          /* Get the array location of this cell */
          ipoint = ires * nresid[1] + jres;
          dist = ires - start_ires;
          score = array[ipoint] - gap_penalty;

          /* If this is the highest score so far, store score and
             position */
          if (score > max_score ||
              (score == max_score && dist < best_dist))
            {
              max_score = score;
              max_pos[0] = ires;
              max_pos[1] = jres;
              best_dist = dist;
            }

          /* Set the gap penalty */
          if (first == FALSE)
            gap_penalty = GAP_PENALTY;
        }

      /* Store best position in the path array */
      if (path_len < max_path - 1)
        {
          /* Get array position */
          ipoint = max_pos[0] * nresid[1] + max_pos[1];

          /* Check for gap in first sequence */
          if (max_pos[0] - last_pos[0] > 1)
            {
              /* Add insertion in first sequence to overall alignment
                 length */
              align_len = align_len + max_pos[0] - last_pos[0] - 1;
            }

          /* Check for gap in second sequence */
          else if (max_pos[1] - last_pos[1] > 1)
            {
              /* Add insertion in second sequence to overall alignment
                 length */
              align_len = align_len + max_pos[1] - last_pos[1] - 1;
            }

          /* Increment total alignment length */
          align_len++;

          /* Store current position */
          last_pos[0] = max_pos[0];
          last_pos[1] = max_pos[1];

          /* Store in path */
          path[path_len] = ipoint;
          path_len++;
        }

      /* If have exceeded maximum allowed path length, then stop
         here */
      else
        {
          printf("*** ERROR. Alignment path length exceeded!!\n");
          printf("***        Serious program error ...\n");
          exit(-1);
        }

      /* Set the new start-position as diagonally right and below
         the point we have just reached */
      start_ires = max_pos[0] + 1;
      start_jres = max_pos[1] + 1;

      /* Check if we have reached the end of the array */
      if (start_ires > nresid[0] - 1 || start_jres > nresid[1] - 1)
        {
          /* If haven't finished on the very last cell of the array,
             make path end one cell beyond it */
          if (start_ires != nresid[0] || start_jres != nresid[1])
            {
              /* Check for gap in first sequence */
              if (nresid[0] - last_pos[0] > 1)
                {
                  /* Add insertion in first sequence to overall alignment
                     length */
                  align_len = align_len + nresid[0] - last_pos[0] - 1;
                }

              /* Check for gap in second sequence */
              else if (nresid[1] - last_pos[1] > 1)
                {
                  /* Add insertion in second sequence to overall alignment
                     length */
                  align_len = align_len + nresid[1] - last_pos[1] - 1;
                }

              /* Store in path */
              path[path_len] = -1;
              path_len++;
            }

          /* Set flag that we're done */
          done = TRUE;
        }

      /* Unset flag */
      first = FALSE;
    }

  /* Return total length of the alignment */
  *alnlen = align_len;

  /* Return the alignment path length */
  return(path_len);
}
/***********************************************************************

get_mapping  -  Use the sequence alignment to get a mapping between
                the two sequences

***********************************************************************/

int get_mapping(char *sequence[2],int nresid[2],int *array,int *path,
                int path_len,int align_len,int *mapping,int nmap,
/* v.2.0.1 --> */
//                float *pc)
                float *pc,struct parameters params)
/* <-- v.2.0.1 */
{
  char aa_code1, aa_code2;

  int align_pos, array_pos, ipoint, ipos, nres[2];
  int ires, jres, kres, last_ires, last_jres;
  int naligned, nidentical, noverlap, nshort;
/* v.2.0.1 --> */
//  int ok, overlap, verbose;
  int ok, overlap;
/* <-- v.2.0.1 */

  int *alignment;

  float percentage;

  /* Initialise variables */
  align_pos = 0;
  last_ires = -1;
  last_jres = -1;
  ok = FALSE;
  overlap = FALSE;
  noverlap = 0;
  naligned = 0;
  nidentical = 0;
  nres[0] = nres[1] = 0;
/* v.2.0.1 --> */
//  verbose = FALSE;
/* <-- v.2.0.1 */

  /* Create array that will hold the alignment */
  alignment = (int *) malloc(2 * align_len * sizeof(int));
  if (alignment == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for alignment array\n");
      exit(1);
    }

  /* Prepare to list the alignment */
/* v.2.0.1 --> */
//  if (verbose == TRUE)
  if (params.verbose == TRUE)
/* <-- v.2.0.1 */
    {
      printf("\n");
      printf("Alignment\n");
      printf("---------\n");
    }

  /* Initialise alignment array with NULLs */
  for (ipos = 0; ipos < 2 * align_len; ipos++)
    alignment[ipos] = -1;

  /* Loop over the points along the path */
  for (ipos = 0; ipos < path_len && align_pos < align_len; ipos++)
    {
      /* Get this point in the path */
      ipoint = path[ipos];

      /* Extract row and column number of this cell */
      if (ipoint == -1)
        {
          ires = nresid[0];
          jres = nresid[1];
          overlap = FALSE;
        }
      else
        {
          ires = ipoint / nresid[1];
          jres = ipoint - ires * nresid[1];
        }

      /* Check for gap in second sequence */
      if (ires - last_ires > 1)
        {
          /* Store insertion from first sequence */
          for (kres = last_ires + 1; kres < ires; kres++)
            {
              /* Store alignment */
              array_pos = align_pos;
              alignment[array_pos] = kres;
              array_pos = align_len + align_pos;
              alignment[array_pos] = -1;

              /* Increment alignment position */
              align_pos++;

              /* Increment overlap length */
              if (overlap == TRUE)
                noverlap++;

              /* Get this residue from the first sequence */
              aa_code1 = sequence[REF_SEQUENCE][kres];

              /* Write out this alignment position */
/* v.2.0.1 --> */
//              if (verbose == TRUE)
              if (params.verbose == TRUE)
/* <-- v.2.0.1 */
                printf("%c  x\n",aa_code1);

              /* Increment count of residues in this sequence */
              if (aa_code1 != 'X')
                nres[REF_SEQUENCE]++;
            }
        }

      /* Check for gap in first sequence */
      else if (jres - last_jres > 1)
        {
          /* Write out insertion from second sequence */
          for (kres = last_jres + 1; kres < jres; kres++)
            {
              /* Store alignment */
              array_pos = align_pos;
              alignment[array_pos] = -1;
              array_pos = align_len + align_pos;
              alignment[array_pos] = kres;

              /* Increment alignment position */
              align_pos++;

              /* Increment overlap length */
              if (overlap == TRUE)
                noverlap++;

              /* Get this residue from the second sequence */
              aa_code2 = sequence[CURR_SEQUENCE][kres];

              /* Write out this alignment position */
/* v.2.0.1 --> */
//              if (verbose == TRUE)
              if (params.verbose == TRUE)
/* <-- v.2.0.1 */
                printf("x  %c\n",aa_code2);

              /* Increment count of residues in this sequence */
              if (aa_code2 != 'X')
                nres[CURR_SEQUENCE]++;
            }
        }

      /* If valid, store the current alignment position */
      if (ires < nresid[0] && jres < nresid[1])
        {
          /* Get the aligned residues */
          aa_code1 = sequence[REF_SEQUENCE][ires];
          aa_code2 = sequence[CURR_SEQUENCE][jres];

          /* Store alignment */
          array_pos = align_pos;
          alignment[array_pos] = ires;
          array_pos = align_len + align_pos;
          alignment[array_pos] = jres;

          /* Increment alignment position */
          align_pos++;
          
          /* If resseqs identical, increment count of sequence
             identities */
          if (aa_code1 != 'X' && aa_code2 != 'X')
            {
              if (aa_code1 == aa_code2)
                nidentical++;
              naligned++;
            }

          /* Increment count of residues in each sequence */
          if (aa_code1 != 'X')
            nres[REF_SEQUENCE]++;
          if (aa_code2 != 'X')
            nres[CURR_SEQUENCE]++;

          /* Store the mapping */
          mapping[jres] = ires;

          /* Increment overlap length */
          overlap = TRUE;
          noverlap++;

          /* If in verbose mode, write out the mapping */
/* v.2.0.1 --> */
//          if (verbose == TRUE)
          if (params.verbose == TRUE)
/* <-- v.2.0.1 */
            {
              /* Write out the first residue */
              printf("%c",aa_code1);

              /* Write out the similarity markers */
              if (aa_code1 == aa_code2)
                printf("  ===  ");
              else
                printf("  ---  ");

              /* Write out the second residue */
              printf("%c\n",aa_code2);
            }
        }

      /* Store current position */
      last_ires = ires;
      last_jres = jres;
    }

  /* Get length of shortest sequence */
  nshort = nres[0];
  if (nres[1] < nshort)
    nshort = nres[1];
  if (nshort == 0)
    nshort = 1;

  /* Calculate the percentage identity */
  if (nshort > 0)
    percentage = 100.0 * (float) nidentical / (float) nshort;
  else
    percentage = 0;

  /* Show stats */
  printf("Number of identical residues:     %8d\n",nidentical);
  printf("Residue overlap:                  %8d\n",noverlap);
  printf("Percentage identity:              %8.2f\n",percentage);

  /* Free up the memory used for the alignment */
  free(alignment);

  /* Determine whether the alignment is good enough to use */
/* v.2.0 --> */
//  if (percentage > MIN_PERCENTAGE && noverlap > MIN_OVERLAP)
  if (percentage > MIN_PERCENTAGE &&
      (noverlap > MIN_OVERLAP ||
       (noverlap > nshort / 3 && nshort > MIN_NSHORT)))
/* <-- v.2.0 */
    ok = TRUE;
  else
    percentage = 0;

  /* Return the percentage */
  *pc = percentage;

  /* Return number of residues aligned */
  return(naligned);
}
/***********************************************************************

perform_alignment  -  Perform a standard sequence alignment

***********************************************************************/

int perform_alignment(char *sequence[2],int nresid[2],int *array,
                      int *tmp_array,int *path,int *mapping,int nmap,
/* v.2.0 --> */
//                      float *pc)
/* v.2.0.1 --> */
//                      float *pc,struct residue **residue_index_ptr[2])
                      float *pc,struct residue **residue_index_ptr[2],
                      struct parameters params)
/* <-- v.2.0.1 */
/* <-- v.2.0 */
{
  int align_len, longest_fit, nfit, path_len;
  int naligned;

  float percentage;

  /* Initialise variables */
  align_len = 0;
  longest_fit = 0;
  nfit = 0;

  /* Initialise array with standard residue similarities */
/* v.2.0 --> */
//  fill_array(sequence,nresid,array,tmp_array,mapping);
  fill_array(sequence,nresid,array,tmp_array,mapping,
             residue_index_ptr);
/* <-- v.2.0 */

  /* Apply Needleman and Wunsch algorithm to adjust
     scores matrix */
  nw_sweep(array,nresid);

  /* Find best route through matrix to give best alignment
     of the two sequences */
  path_len = find_best_path(array,nresid,path,&align_len);

  /* Use the alignment to get a mapping between the two sequences */
  naligned = get_mapping(sequence,nresid,array,path,path_len,align_len,
/* v.2.0.1 --> */
//                         mapping,nmap,&percentage);
                         mapping,nmap,&percentage,params);
/* <-- v.2.0.1 */

  /* Return the percentage id */
  *pc = percentage;

  /* Return number of residues aligned */
  return(naligned);
}
/* v.1.3 --> */
/***********************************************************************

match_up_residues  -  Pair off the residues by matching up the alignment
                      positions

***********************************************************************/

int match_up_residues(struct residue **residue_index_ptr[2],int *seqlen,
/* v.2.0.1 --> */
//                      struct resseq ***aln)
                      struct resseq ***aln,struct parameters params)
/* <-- v.2.0.1 */
{
  int i, ires, jres, nfit, nresid[2];
  int found;

  struct residue *residue_ptr[2];
  struct resseq *resseq_ptr[2], *first_resseq_ptr[2], *last_resseq_ptr[2];

  struct resseq **alignment;

  /* Initialise variables */
  nfit = 0;
  nresid[0] = nresid[1] = 0;

  /* Initialise pointers */
  alignment = NULL;
  first_resseq_ptr[0] = first_resseq_ptr[1] = NULL;
  last_resseq_ptr[0] = last_resseq_ptr[1] = NULL;

  /* Loop over the residues in the first sequence */
  for (ires = 0; ires < seqlen[REF_SEQUENCE]; ires++)
    {
      /* Get this residue */
      residue_ptr[REF_SEQUENCE]
        = residue_index_ptr[REF_SEQUENCE][ires];

      /* If residue has an alignment position, look for the equivalent
         position on the other sequence */
      if (residue_ptr[REF_SEQUENCE]->aln_pos > -1)
        {
          /* Initialise flag */
          found = FALSE;

          /* Loop over the residues in the other sequence to find 
             residue at this alignment position */
          for (jres = 0; jres < seqlen[CURR_SEQUENCE] && found == FALSE;
               jres++)
            {
              /* Get this residue */
              residue_ptr[CURR_SEQUENCE]
                = residue_index_ptr[CURR_SEQUENCE][jres];

              /* If residue has an alignment position, look for the equivalent
                 position on the other sequence */
              if (residue_ptr[CURR_SEQUENCE]->aln_pos
                  == residue_ptr[REF_SEQUENCE]->aln_pos)
                {
                  /* Create a resseq record for both residues */
                  for (i = 0; i < 2; i++)
                    {
                      /* Create record */
                      resseq_ptr[i]
                        = create_resseq_record(&(first_resseq_ptr[i]),
                                               &(last_resseq_ptr[i]),
                                               residue_ptr[i]);

                      /* Increment count of store residues */
                      nresid[i]++;
                    }

                  /* Store pointers between the two resseq records */
                  resseq_ptr[REF_SEQUENCE]->align_resseq_ptr
                    = resseq_ptr[CURR_SEQUENCE];
                  resseq_ptr[CURR_SEQUENCE]->align_resseq_ptr
                    = resseq_ptr[REF_SEQUENCE];

                  /* Increment count of equivalenced residue */
                  nfit++;

                  /* Set flag that match found */
                  found = TRUE;
                }
            }
        }
    }

  /* If have enough equivalenced residues, store them in an array for
     use in the fitting routines */
  if (nfit > 2)
/* v.2.0.1 --> */
//    store_alignment(first_resseq_ptr,&alignment,nfit);
    store_alignment(first_resseq_ptr,&alignment,nfit,params);
/* <-- v.2.0.1 */

  /* Return the alignment */
  *aln = alignment;

  /* Return the number of fitted residues */
  return(nfit);
}
/***********************************************************************

match_alignment  -  Match each sequence from the FASTA-type alignment
                    file to the sequence in the corresponding PDB chain

***********************************************************************/

int match_alignment(struct pdbentry *pdbentry_ptr[2],
                    struct residue **residue_index_ptr[2],
                    char *aln_sequence[2],char *sequence[2],
/* v.1.3.7 --> */
//                    int aln_len[2],struct resseq ***aln)
                    int aln_len[2],char aln_chain[2],
/* v.2.0.1 --> */
//                    struct resseq ***aln)
                    struct resseq ***aln,struct parameters params)
/* <-- v.2.0.1 */
/* <-- v.1.3.7 */
{
  int i, imap, ipdb, ipos, naligned, nfit, nmap, seqlen[2];
  int ires, jres, nresid[2];

  int *array, *mapping, *orig_pos, *path, *tmp_array;

  float percentage;

  struct residue *residue_ptr;

  struct resseq **alignment;

  /* Initialise variables */
  nfit = 0;

  /* Loop over the sequences and align them to the residues in the
     corresponding PDB entries */
  for (ipdb = 0; ipdb < 2; ipdb++)
    {
      /* Get the relevant sequence from the corresponding PDB entry */
      seqlen[0]
        = nresid[ipdb]
        = get_sequence(pdbentry_ptr[ipdb],&(sequence[0]),
/* v.1.3.7 --> */
//                       &(residue_index_ptr[ipdb]));
                       &(residue_index_ptr[ipdb]),aln_chain[ipdb]);
/* <-- v.1.3.7 */

      /* Create a version of the sequence from the alignment file
         without gaps */
      sequence[1] = (char *) malloc(sizeof(char)*(aln_len[ipdb] + 1));
      orig_pos = (int *) malloc(sizeof(int)*(aln_len[ipdb] + 1));

      /* Loop over the sequence to copy across the non-gap residues */
      ipos = 0;
      for (i = 0; i < aln_len[ipdb]; i++)
        {
          /* If not a gap, copy across */
          if (aln_sequence[ipdb][i] != '-')
            {
              /* Copy this residue across */
              sequence[1][ipos] = aln_sequence[ipdb][i];

              /* Save its original position in the multiple alignment */
              orig_pos[ipos] = i;

              /* Increment position and null terminate */
              ipos++;
              sequence[1][ipos] = '\0';
            }
        }

      /* Get the length of this sequence */
      seqlen[1] = strlen(sequence[1]);

      /* Create the arrays for the sequence alignment */
      nmap
        = create_alignment_arrays(seqlen,&array,&tmp_array,&path,
                                  &mapping);

      /* Perform the alignment */
      naligned
        = perform_alignment(sequence,seqlen,array,tmp_array,path,
/* v.2.0 --> */
//                            mapping,nmap,&percentage);
/* v.2.0.1 --> */
//                            mapping,nmap,&percentage,residue_index_ptr);
                            mapping,nmap,&percentage,residue_index_ptr,
                            params);
/* <-- v.2.0.1 */
/* <-- v.2.0 */

      /* Free up memory */
      free(array);
      free(tmp_array);
      free(path);

      /* For each residue, store its alignment position */
      for (imap = 0; imap < nmap; imap++)
        {
          /* Get the residues mapped */
          jres = mapping[imap];
          ires = imap;

          /* If have a mapping, then store it */
          if (jres > -1)
            {
              /* Get the corresponding residue */
              residue_ptr = residue_index_ptr[ipdb][jres];

              /* Store this residue's original alignment position */
              residue_ptr->aln_pos = orig_pos[ires];
            }
        }
    }

  /* Pair off the residues by matching up the alignment positions */
/* v.2.0.1 --> */
//  nfit = match_up_residues(residue_index_ptr,nresid,&alignment);
  nfit = match_up_residues(residue_index_ptr,nresid,&alignment,params);
/* <-- v.2.0.1 */

  /* Return the alignment */
  *aln = alignment;

  /* Return the number of matched residues */
  return(nfit);
}
/***********************************************************************

get_caf_alignment  -  Read in the FASTA-format alignment file from CORA

***********************************************************************/

int get_caf_alignment(char *file_name,
                      struct pdbentry *first_pdbentry_ptr,
/* v.2.0.1 --> */
//                      struct resseq ***aln)
                      struct resseq ***aln,struct parameters params)
/* <-- v.2.0.1 */
{
/* v.1.3.7 --> */
//  char chain, input_line[LINELEN + 1], pdb_code[5];
  char aln_chain[2], chain, input_line[LINELEN + 1], pdb_code[5];
/* <-- v.1.3.7 */

  char *aln_sequence[2], *sequence[2], *tmp_sequence, *string_ptr;

  int i, ipdb, len, line, nfit, aln_len[2], start_pos;
  int skip;

  struct pdbentry *pdbentry_ptr[2];
  struct resseq *first_resseq_ptr[2], *last_resseq_ptr[2];

  struct residue **residue_index_ptr[2];
  struct resseq **alignment;

  FILE *file_ptr;

  /* Initialise variables */
  alignment = NULL;
  aln_len[0] = aln_len[1] = 0;
  line = 0;
  nfit = 0;
  aln_sequence[0] = aln_sequence[1] = &null;
/* v.1.3.7 --> */
  aln_chain[0] = aln_chain[1] = '\0';
/* <-- v.1.3.7 */
  skip = FALSE;

  /* Initialise pointers */
  pdbentry_ptr[0] = pdbentry_ptr[1] = NULL;
  pdbentry_ptr[CURR_SEQUENCE] = first_pdbentry_ptr;
  if (first_pdbentry_ptr != NULL)
    pdbentry_ptr[REF_SEQUENCE] = first_pdbentry_ptr->next_pdbentry_ptr;
  if (pdbentry_ptr[0] == NULL || pdbentry_ptr[1] == NULL)
    return(nfit);
  first_resseq_ptr[0] = first_resseq_ptr[1] = NULL;
  last_resseq_ptr[0] = last_resseq_ptr[1] = NULL;

  /* Open the alignments file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    {
      printf("\n*** Warning. Unable to open alignments file [%s]\n",
             file_name);
      return(nfit);
    }

  /* Loop through the file, line by line */
  while (fgets(input_line,LINELEN,file_ptr) != NULL)
    {
      /* Truncate the string to strip off the newline */
      string_truncate(input_line,LINELEN);

      /* Process if this is a sequence line */
      if (!strncmp(input_line,"pdb|",4))
        {
          /* Extract the PDB code and chain ID, if any */
          strncpy(pdb_code,input_line + 4,4);
          pdb_code[4] = '\0';
          chain = input_line[8];

          /* Convert to lower case */
          to_lower(pdb_code,4);

          /* Check to see if PDB code matches either of
             our codes of interest */
          ipdb = -1;
          for (i = 0; i < 2 && ipdb == -1; i++)
            {
              /* Check if this the right PDB code and chain id */
/* v.1.3.7 --> */
//              if (!strcmp(pdb_code,pdbentry_ptr[i]->pdb_code) &&
//                  (pdbentry_ptr[i]->chain == chain || chain == ' '))
//                ipdb = i;
              if (!strcmp(pdb_code,pdbentry_ptr[i]->pdb_code) &&
                  aln_chain[i] == '\0')
                {
                  ipdb = i;
                  aln_chain[i] = chain;
                }
/* <-- v.1.3.7 */
            }

          /* Initialise sequence start-position */
          start_pos = -1;

          /* If wanted, then process */
          if (ipdb > -1 && skip == FALSE)
            {
              /* Find where the sequence starts */
              string_ptr = strchr(input_line,' ');
              if (string_ptr != NULL)
                {
                  /* get length of line */
                  len = strlen(string_ptr);

                  /* Loop until we hit first non-blank character */
                  for (i = 0; i < len && start_pos == -1; i++)
                    if (string_ptr[i] != ' ')
                      start_pos = i;
                  
                }

              /* Set the skip flag so that the associated secondary
                 structure line is skipped */
              skip = TRUE;
            }

          /* Otherwise, unset the skip flag */
          else
            skip = FALSE;

          /* If have found start of sequence, append to what
             we already have */
          if (ipdb > -1 && start_pos > -1)
            {
              /* If this is the first chunk, then store */
              if (aln_sequence[ipdb] == &null)
                store_string(&(aln_sequence[ipdb]),
                             string_ptr + start_pos);

              /* Otherwise, append to end of current sequence */
              else
                {
                  /* Get length of sequence required */
                  len = strlen(aln_sequence[ipdb])
                    + strlen(string_ptr + start_pos);

                  /* Create enough memory to store new segment */
                  tmp_sequence
                    = (char *) malloc(sizeof(char)*(len + 1));

                  /* Copy the sequence fragments across */
                  strcpy(tmp_sequence,aln_sequence[ipdb]);
                  strcat(tmp_sequence,string_ptr + start_pos);

                  /* Free up memory of old fragment */
                  free(aln_sequence[ipdb]);

                  /* Copy new fragment across */
                  aln_sequence[ipdb] = tmp_sequence;

                  /* Get the length of the current sequence */
                  aln_len[ipdb] = strlen(aln_sequence[ipdb]);
                }
            }
        }

      /* Increment line count */
      line++;
    }

  /* Close the alignments file */
  fclose(file_ptr);

  /* If either sequence not found, the can't fir */
  if (aln_len[0] == 0 || aln_len[1] == 0)
    return(nfit);

  /* Match the read-in sequences to the residues read in from the
     corresponding PDB chain */
  nfit = match_alignment(pdbentry_ptr,residue_index_ptr,aln_sequence,
/* v.1.3.7 --> */
//                         sequence,aln_len,&alignment);
/* v.2.0.1 --> */
//                         sequence,aln_len,aln_chain,&alignment);
                         sequence,aln_len,aln_chain,&alignment,params);
/* <-- v.2.0.1 */
/* <-- v.1.3.7 */

  /* Return the alignment */
  *aln = alignment;

  /* Return number of residue equivalences found */
  return(nfit);
}
/***********************************************************************

get_fasta_alignment  -  Read in the FASTA-format multiple alignment

***********************************************************************/

int get_fasta_alignment(char *file_name,
                        struct pdbentry *first_pdbentry_ptr,
/* v.2.0.1 --> */
//                        struct resseq ***aln)
                        struct resseq ***aln,struct parameters params)
/* <-- v.2.0.1 */
{
/* v.1.3.7 --> */
//  char chain, input_line[LINELEN + 1], pdb_code[5];
  char aln_chain[2], chain, input_line[LINELEN + 1], pdb_code[5];
/* <-- v.1.3.7 */

  char *aln_sequence[2], *sequence[2], *tmp_sequence;

  int i, ipdb, ipos, len, line, nfit, aln_len[2];
  int read_sequence;

  struct pdbentry *pdbentry_ptr[2];
  struct resseq *first_resseq_ptr[2], *last_resseq_ptr[2];

  struct residue **residue_index_ptr[2];
  struct resseq **alignment;

  FILE *file_ptr;

  /* Initialise variables */
  alignment = NULL;
  aln_len[0] = aln_len[1] = 0;
  line = 0;
  nfit = 0;
/* v.1.3.7 --> */
  aln_chain[0] = aln_chain[1] = '\0';
/* <-- v.1.3.7 */
  aln_sequence[0] = aln_sequence[1] = &null;
  read_sequence = FALSE;

  /* Initialise pointers */
  pdbentry_ptr[0] = pdbentry_ptr[1] = NULL;
  pdbentry_ptr[CURR_SEQUENCE] = first_pdbentry_ptr;
  if (first_pdbentry_ptr != NULL)
    pdbentry_ptr[REF_SEQUENCE] = first_pdbentry_ptr->next_pdbentry_ptr;
  if (pdbentry_ptr[0] == NULL || pdbentry_ptr[1] == NULL)
    return(nfit);
  first_resseq_ptr[0] = first_resseq_ptr[1] = NULL;
  last_resseq_ptr[0] = last_resseq_ptr[1] = NULL;

  /* Open the alignments file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    {
      printf("\n*** Warning. Unable to open alignments file [%s]\n",
             file_name);
      return(nfit);
    }

  /* Loop through the file, line by line */
  while (fgets(input_line,LINELEN,file_ptr) != NULL)
    {
      /* Truncate the string to strip off the newline */
      len = string_truncate(input_line,LINELEN);

      /* If this is the start of a new sequence, see if it is wanted */
      if (input_line[0] == '>')
        {
          /* Initialise flag */
          read_sequence = FALSE;

          /* Check for 'pdb' prefix */
          ipos = 1;
          if (len > 4 && !strncmp(input_line,">pdb",4))
            ipos = 4;

          /* Check for divider */
          if (input_line[ipos] == '|')
            ipos++;

          /* If line long enough to hold PDB code, extract it from next
             4 characters */
          if (ipos + 4 <= len)
            {
              /* Extract the PDB code */
              strncpy(pdb_code,input_line + ipos,4);
              pdb_code[4] = '\0';

              /* Convert to lower case */
              to_lower(pdb_code,4);

              /* Shift position along */
              ipos = ipos + 4;

              // Check for possible chain identifier
              if (input_line[ipos] == ':' || input_line[ipos] == '|'
                  || input_line[ipos] == '(')
                ipos++;

              /* Extract the chain id, if there is one */
              chain = ' ';
              if (ipos < len)
                chain = input_line[ipos];

              /* Check to see if PDB code matches either of
                 our codes of interest */
              ipdb = -1;
              for (i = 0; i < 2 && ipdb == -1; i++)
                {
                  /* Check if this the right PDB code and chain id */
/* v.1.3.7 --> */
//                  if (!strcmp(pdb_code,pdbentry_ptr[i]->pdb_code) &&
//                      (pdbentry_ptr[i]->chain == chain || chain == ' '))
//                    ipdb = i;
                  if (!strcmp(pdb_code,pdbentry_ptr[i]->pdb_code) &&
                      aln_chain[i] == '\0')
                    {
                      ipdb = i;
                      aln_chain[i] = chain;
                    }
/* <-- v.1.3.7 */
                }

              /* If matched, then set flag */
              if (ipdb > -1)
                read_sequence = TRUE;
            }
        }

      /* Otherwise, if we're reading a sequence, append it to growing
         sequence */
      else if (read_sequence == TRUE)
        {
          /* If this is the first chunk, then store */
          if (aln_sequence[ipdb] == &null)
              store_string(&(aln_sequence[ipdb]),input_line);

          /* Otherwise, append to end of current sequence */
          else
            {
              /* Get length of sequence required */
              len = strlen(aln_sequence[ipdb]) + strlen(input_line);

              /* Create enough memory to store new segment */
              tmp_sequence = (char *) malloc(sizeof(char)*(len + 1));

              /* Copy the sequence fragments across */
              strcpy(tmp_sequence,aln_sequence[ipdb]);
              strcat(tmp_sequence,input_line);

              /* Free up memory of old fragment */
              free(aln_sequence[ipdb]);

              /* Copy new fragment across */
              aln_sequence[ipdb] = tmp_sequence;
            }

/* v.1.3.7 --> */
          /* Get the length of the current sequence */
          aln_len[ipdb] = strlen(aln_sequence[ipdb]);
/* <-- v.1.3.7 */
        }

      /* Increment line count */
      line++;
    }

  /* Close the alignments file */
  fclose(file_ptr);

  /* If either sequence not found, the can't fir */
  if (aln_len[0] == 0 || aln_len[1] == 0)
    return(nfit);

  /* Match the read-in sequences to the residues read in from the
     corresponding PDB chain */
  nfit = match_alignment(pdbentry_ptr,residue_index_ptr,aln_sequence,
/* v.1.3.7 --> */
//                         sequence,aln_len,&alignment);
/* v.2.0.1 --> */
//                         sequence,aln_len,aln_chain,&alignment);
                         sequence,aln_len,aln_chain,&alignment,params);
/* <-- v.2.0.1 */
/* <-- v.1.3.7 */

  /* Return the alignment */
  *aln = alignment;

  /* Return number of residue equivalences found */
  return(nfit);
}
/* <-- v.1.3 */
/***********************************************************************

get_residue_mappings  -  Extract the residue-residue mappings from the
                         alignment

***********************************************************************/

int get_residue_mappings(struct residue **residue_index_ptr[2],
/* v.2.0 --> */
//                         int *mapping,int nmap,struct resseq ***aln)
                         int *mapping,int nmap,struct resseq ***aln,
/* v.2.0.1 --> */
//                         int antibody)
                         struct parameters params)
/* <-- v.2.0.1 */
/* <-- v.2.0 */
{
  int i, imap, ires, jres, nfit, nresid[2];
/* v.2.0 --> */
  int wanted;
/* <-- v.2.0 */

  struct residue *residue_ptr[2];
  struct resseq *resseq_ptr[2], *first_resseq_ptr[2], *last_resseq_ptr[2];

  struct resseq **alignment;

  /* Initialise variables */
  nfit = 0;
  nresid[0] = nresid[1] = 0;

  /* Initialise pointers */
  alignment = NULL;
  first_resseq_ptr[0] = first_resseq_ptr[1] = NULL;
  last_resseq_ptr[0] = last_resseq_ptr[1] = NULL;
/* v.2.0.2 --> */
  resseq_ptr[0] = resseq_ptr[1] = NULL;
/* <-- v.2.0.2 */

  /* Loop over the residue mappings */
  for (imap = 0; imap < nmap; imap++)
    {
      /* Get the residues mapped */
      jres = mapping[imap];
      ires = imap;

      /* If have a mapping, then store it */
      if (jres > -1)
        {
          /* Get the two residue records */
          residue_ptr[REF_SEQUENCE]
            = residue_index_ptr[REF_SEQUENCE][jres];
          residue_ptr[CURR_SEQUENCE]
            = residue_index_ptr[CURR_SEQUENCE][ires];

          /* If both residues valid, then store mapping between them */
          if (residue_ptr[REF_SEQUENCE] != NULL &&
              residue_ptr[CURR_SEQUENCE] != NULL &&
              strcmp(residue_ptr[REF_SEQUENCE]->res_name,"DUM") &&
              strcmp(residue_ptr[CURR_SEQUENCE]->res_name,"DUM"))
            {
/* v.2.0 --> */
              /* If this is an antibody plot, check that the residue
                 numbers match */
              wanted = TRUE;
/* v.2.0.1 --> */
//              if (antibody == TRUE &&
              if (params.antibody == TRUE &&
/* <-- v.2.0.1 */
                  (strcmp(residue_ptr[REF_SEQUENCE]->res_name,
                          residue_ptr[CURR_SEQUENCE]->res_name) ||
                   strcmp(residue_ptr[REF_SEQUENCE]->res_num,
                          residue_ptr[CURR_SEQUENCE]->res_num)))
                wanted = FALSE;

              /* Create a resseq record for both residues */
//              for (i = 0; i < 2; i++)
              for (i = 0; i < 2 && wanted == TRUE; i++)
/* <-- v.2.0 */
                {
                  /* Create record */
                  resseq_ptr[i]
                    = create_resseq_record(&(first_resseq_ptr[i]),
                                           &(last_resseq_ptr[i]),
                                           residue_ptr[i]);

                  /* Increment count of store residues */
                  nresid[i]++;
                }

              /* Store pointers between the two resseq records */
/* v.2.0.2 --> */
              if (wanted == TRUE)
                {
/* <-- v.2.0.2 */
                  resseq_ptr[REF_SEQUENCE]->align_resseq_ptr
                    = resseq_ptr[CURR_SEQUENCE];
                  resseq_ptr[CURR_SEQUENCE]->align_resseq_ptr
                    = resseq_ptr[REF_SEQUENCE];

                  /* Increment count of equivalenced residue */
                  nfit++;
/* v.2.0.2 --> */
                }
/* <-- v.2.0.2 */
            }
        }
    }

  /* If have enough equivalenced residues, store them in an array for use
     in the fitting routines */
  if (nfit > 2)
/* v.2.0.1 --> */
//    store_alignment(first_resseq_ptr,&alignment,nfit);
    store_alignment(first_resseq_ptr,&alignment,nfit,params);
/* v.2.0.1 --> */

  /* Return the alignment */
  *aln = alignment;

  /* Return the number of fitted residues */
  return(nfit);
}
/***********************************************************************

get_seq_alignment  -  Align the two protein sequences to get residue
                      equivalences

***********************************************************************/

int get_seq_alignment(struct pdbentry *first_pdbentry_ptr,
/* v.2.0 --> */
//                      struct resseq ***aln)
                      struct resseq ***aln,struct parameters params)
/* <-- v.2.0 */
{
  char *sequence[2];

  int naligned, nfit, nmap, seqlen[2];
/* v.2.0 --> */
//  int ok;
  int in_plot;
/* <-- v.2.0 */

  int *array, *mapping, *path, *tmp_array;

  float percentage;

  struct pdbentry *pdbentry_ptr[2];

  struct residue **residue_index_ptr[2];

  struct resseq **alignment;

  /* Initialise variables */
/* v.2.0 --> */
  in_plot = FALSE;
/* <-- v.2.0 */
  nfit = 0;

  /* Initialise pointers */
  alignment = NULL;
  pdbentry_ptr[CURR_SEQUENCE] = pdbentry_ptr[REF_SEQUENCE] = NULL;
  pdbentry_ptr[CURR_SEQUENCE] = first_pdbentry_ptr;
  if (first_pdbentry_ptr != NULL)
    pdbentry_ptr[REF_SEQUENCE] = first_pdbentry_ptr->next_pdbentry_ptr;
  if (pdbentry_ptr[CURR_SEQUENCE] == NULL ||
      pdbentry_ptr[REF_SEQUENCE] == NULL)
    return(nfit);

  /* Get the relevant sequence from our PDB entry */
/* v.2.0 --> */
  if (params.dimplot == TRUE)
    seqlen[CURR_SEQUENCE]
      = get_dimplot_sequence(pdbentry_ptr[CURR_SEQUENCE],
                             &(sequence[CURR_SEQUENCE]),
                             &(residue_index_ptr[CURR_SEQUENCE]),in_plot);
  else
/* <-- v.2.0 */
    seqlen[CURR_SEQUENCE]
      = get_sequence(pdbentry_ptr[CURR_SEQUENCE],
                     &(sequence[CURR_SEQUENCE]),
                     &(residue_index_ptr[CURR_SEQUENCE]),'\0');

  /* If no interacting residues found, then return empty-handed */
  if (seqlen[CURR_SEQUENCE] < 2)
    {
      printf("*** ERROR. No interacting residues found in sequence for "
             "%s%c\n",pdbentry_ptr[CURR_SEQUENCE]->pdb_code,
             pdbentry_ptr[CURR_SEQUENCE]->chain);
      return(nfit);
    }

  /* Repeat for reference PDB entry */
/* v.2.0 --> */
  if (params.dimplot == TRUE)
    seqlen[REF_SEQUENCE]
      = get_dimplot_sequence(pdbentry_ptr[REF_SEQUENCE],
                             &(sequence[REF_SEQUENCE]),
                             &(residue_index_ptr[REF_SEQUENCE]),in_plot);
  else
/* <-- v.2.0 */
    seqlen[REF_SEQUENCE]
      = get_sequence(pdbentry_ptr[REF_SEQUENCE],
                     &(sequence[REF_SEQUENCE]),
                     &(residue_index_ptr[REF_SEQUENCE]),'\0');

  /* If no interacting residues found, then return empty-handed */
  if (seqlen[REF_SEQUENCE] < 2)
    {
      printf("*** ERROR. No interacting residues found in sequence for "
             "%s%c\n",pdbentry_ptr[REF_SEQUENCE]->pdb_code,
             pdbentry_ptr[REF_SEQUENCE]->chain);
      return(nfit);
    }

  /* Create the arrays for the sequence alignment */
/* v.1.3.7 --> */
  printf("Aligning chain %c of %s with chain %c of %s\n",
         pdbentry_ptr[REF_SEQUENCE]->map_chain,
         pdbentry_ptr[REF_SEQUENCE]->pdb_code,
         pdbentry_ptr[CURR_SEQUENCE]->map_chain,
         pdbentry_ptr[CURR_SEQUENCE]->pdb_code);
/* <-- v.1.3.7 */
  nmap
    = create_alignment_arrays(seqlen,&array,&tmp_array,&path,&mapping);

  /* Perform the alignment */
  naligned
    = perform_alignment(sequence,seqlen,array,tmp_array,path,mapping,
/* v.2.0 --> */
//                        nmap,&percentage);
/* v.2.0.1 --> */
//                        nmap,&percentage,residue_index_ptr);
                        nmap,&percentage,residue_index_ptr,params);
/* <-- v.2.0.1 */
/* <-- v.2.0 */

  /* Free up memory */
  free(array);
  free(tmp_array);
  free(path);

  /* If alignment is OK to use, determine the residue equivalences */
  if (naligned > MIN_ALIGNED)
    {
      /* Extract the residue mappings from the alignment */
      nfit = get_residue_mappings(residue_index_ptr,mapping,nmap,
/* v.2.0 --> */
//                                  &alignment);
/* v.2.0.1 --> */
//                                  &alignment,params.antibody);
                                  &alignment,params);
/* <-- v.2.0.1 */
/* <-- v.2.0 */
    }

  /* Free up memory */
  free(mapping);

  /* Return the alignment */
  *aln = alignment;

  /* Return number of residues equivalenced */
  return(nfit);
}
/***********************************************************************

get_bs_sequence  -  Align the two protein sequences to get residue
                    equivalences

***********************************************************************/

int get_bs_sequence(struct pdbentry *pdbentry_ptr,char **seq,
                    struct residue ***residue_idx_ptr,
                    struct residue *dummy_residue_ptr,int *nc)
{
  char chain, last_chain;
  char insertion, last_insertion, number_string[20];

  char *sequence;

  int i, ipos, iresid, nchains, nresid;
  int last_res_no, res_no;

  struct residue *residue_ptr;

  struct residue **residue_index_ptr;
  
  /* Initialise variables */
  ipos = 0;
  iresid = 0;
  last_chain = chain = '\0';
  last_insertion = '\0';
  nchains = 0;
  nresid = pdbentry_ptr->nresid;

  /* Create an array for storing the protein sequence */
  sequence = (char *) malloc(sizeof(char)*(nresid + 1));
  
  /* Create an array for storing the residue pointers */
  residue_index_ptr
    = (struct residue **) malloc(nresid * sizeof(struct residue *));
  if (residue_index_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for residue index\n");
      exit(1);
    }

  /* Initialise index array */
  for (i = 0; i < nresid; i++)
    residue_index_ptr[i] = NULL;

  /* Re-initialise residue count */
  nresid = 0;
  
  /* Get the full entry's first residue */
  residue_ptr = pdbentry_ptr->first_residue_ptr;

  /* Loop over the residues to find corresponding records in entry read
     in from the ensemble.drw file */
  while (residue_ptr != NULL)
    {
      /* Check that we're not going over the array limit */
      if (ipos > nresid)
        {
          printf("*** PROGRAM ERROR. Number of residues in %s exceeded: "
                 "%d (nresid = %d)\n",pdbentry_ptr->pdb_code,
                 iresid,nresid);
          exit(1);
        }

      /* If this is an interacting protein residue, then create an
         entry for it in the sequence */
      if ((residue_ptr->interaction == NONBOND_TO_LIGAND ||
          residue_ptr->interaction == HBOND_TO_LIGAND ||
          residue_ptr->interaction == COVALENT_TO_LIGAND) &&
          residue_ptr->type != DUMMY &&
          residue_ptr->aa_code != 'X')
        {
          /* Get the numeric residue number */
          strcpy(number_string,residue_ptr->res_num);
          number_string[4] = '\0';
          res_no = atoi(number_string);

          /* Get the insertion code */
          insertion = residue_ptr->res_num[4];

          /* If this is a new chain, increment count of chains */
          if (residue_ptr->chain != last_chain)
            nchains++;

          /* If this residue doesn't follow directly after the previous
             one, and is not the first, add in a dummy residue record */
          if ((residue_ptr->chain != last_chain ||
               (residue_ptr->chain == last_chain &&
                res_no > last_res_no + 1)) &&
              last_chain != '\0')
            {
              /* Save the residue code */
              sequence[ipos] = 'X';

              /* Store the dummy residue pointer */
              residue_index_ptr[ipos] = dummy_residue_ptr;

              /* Increment sequence count */
              ipos++;
                          
              /* Increment count of residues */
              nresid++;
            }

          /* Save the residue code */
          sequence[ipos] = residue_ptr->aa_code;

          /* Store the corresponding residue pointer */
          residue_index_ptr[ipos] = residue_ptr;

          /* Increment sequence count */
          ipos++;

          /* Increment count of residues */
          nresid++;

          /* Save residue number, chain and insertion code */
          last_res_no = res_no;
          last_insertion = insertion;
          last_chain = residue_ptr->chain;
        }

      /* Get pointer to next residue record */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;
    }

  /* Add a dummy residue on the end */
  sequence[ipos] = 'X';

  /* Store the corresponding residue pointer */
  residue_index_ptr[ipos] = dummy_residue_ptr;

  /* Increment sequence count */
  ipos++;

  /* Increment count of residues */
  nresid++;

 /* Terminate the sequence */
  sequence[ipos] = '\0';

  /* Return the sequence */
  *seq = sequence;

  /* Return the array of residue pointers */
  *residue_idx_ptr = residue_index_ptr;

  /* Return number of chains encountered */
  *nc = nchains;

  /* Return the number of residue in the stored sequence */
  return(nresid);
}
/***********************************************************************

copy_sequence  -  Save a copy of the given residue sequence

***********************************************************************/

void copy_sequence(struct residue **ref_residue_index_ptr,int nresid,
                   struct residue ***residue_idx_ptr)
{
  int i;

  struct residue **residue_index_ptr;
  
  /* Create the copy index array */
  residue_index_ptr
    = (struct residue **) malloc(nresid * sizeof(struct residue *));
  if (residue_index_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for residue index\n");
      exit(1);
    }

  /* Loop over the reference index to copy pointers */
  for (i = 0; i < nresid; i++)
    {
      /* Copy the pointer */
      residue_index_ptr[i] = ref_residue_index_ptr[i];
    }

  /* Return the copy */
  *residue_idx_ptr = residue_index_ptr;
}
/***********************************************************************

reverse_chains  -  Sort the residues into reverse chain order

***********************************************************************/

void reverse_chains(struct residue **ref_residue_index_ptr,
                    struct residue **residue_index_ptr,int nresid)
{
  char chain;

  int ipos, istart, jpos;

  struct residue *residue_ptr;

  /* Initialise variables */
  chain = '\0';
  istart = -1;
  jpos = 0;

  /* Loop over all the records to find where the first chain change
     occurs */
  for (ipos = 0; ipos < nresid; ipos++)
    {
      /* Get the corresponding residue */
      residue_ptr = ref_residue_index_ptr[ipos];

      /* If not a dummy record, check its chain */
      if (residue_ptr != NULL &&
          strcmp(residue_ptr->res_name,"DUM"))
        {
          /* If the is the first residue, store the chain id */
          if (chain == '\0')
            chain = residue_ptr->chain;

          /* Otherwise, if chain id has changed, then this is our
             cut-point */
          else if (residue_ptr->chain != chain &&
                   istart == -1)
            {
              /* Save the cut-point */
              istart = ipos;

              /* Save the new chain */
              chain = residue_ptr->chain;
            }
        }

      /* If we are copying across, then copy this residue to reverse
         index */
      if (istart != -1)
        {
          residue_index_ptr[jpos] = residue_ptr;
          jpos++;
        }
    }

  /* If we don't have a cut-point, use the whole sequence */
  if (istart == -1)
    istart = nresid;

  /* Loop over the residues in the first chain to add on end */
  for (ipos = 0; ipos < istart; ipos++)
    {
      residue_index_ptr[jpos] = ref_residue_index_ptr[ipos];
      jpos++;
    }
}
/***********************************************************************
  
get_bs_alignment  -  Align sequences using just the interacting residues
                     in each structure

***********************************************************************/

int get_bs_alignment(struct pdbentry *first_pdbentry_ptr,
/* v.2.0 --> */
//                     struct resseq ***aln)
                     struct resseq ***aln,struct parameters params)
/* <-- v.2.0 */
{
  char residue_name[NAMINO+1][4], residue_id[20];

  char *sequence[4];

  int best_fit, nfit, nmap, nresid[4];
  int loop, nloops, use;
  int nchains = 0;
/* v.2.0 --> */
//  int first;
  int first, in_plot;
/* <-- v.2.0 */

  int *array, *mapping, *path, *tmp_array;

  float best_percentage, percentage;

  struct pdbentry *pdbentry_ptr, *fit_pdbentry_ptr;
  struct residue *dummy_residue_ptr;
  struct residue *first_residue_ptr, *last_residue_ptr;

  struct residue **residue_index_ptr[4];
  struct resseq **alignment;

  /* Initialise variables */
  best_fit = nfit = 0;
  best_percentage = 0;
/* v.2.0 --> */
  in_plot = TRUE;
/* <-- v.2.0 */
  sprintf(residue_id,"[%s %s %c]","DUM","   1 ",'X');

  /* Read in the residue-similarity matrix */
  get_residue_names(residue_name);

  /* Create a dummy residue to represent the missing sequence segments
     between the residues we want aligned */
  first_residue_ptr = last_residue_ptr = NULL;
  dummy_residue_ptr
    = create_residue_record(&first_residue_ptr,&last_residue_ptr,"DUM",
                            "   1 ",'X',DUMMY,"N/A",0.0,0.0,residue_id,
                            first_pdbentry_ptr);

  /* Get the first PDB entry as the entry to be fitted */
  fit_pdbentry_ptr = first_pdbentry_ptr;

  /* Create a set of residue sequence records for the alignment */
/* v.2.0 --> */
  if (params.dimplot == TRUE)
    nresid[CURR_SEQUENCE]
      = get_dimplot_sequence(fit_pdbentry_ptr,&(sequence[CURR_SEQUENCE]),
                             &(residue_index_ptr[CURR_SEQUENCE]),in_plot);
  else
/* <-- v.2.0 */
    nresid[CURR_SEQUENCE]
      = get_bs_sequence(fit_pdbentry_ptr,&(sequence[CURR_SEQUENCE]),
                        &(residue_index_ptr[CURR_SEQUENCE]),
                        dummy_residue_ptr,&nchains);

  /* Get the first of the ensemble entries */
  pdbentry_ptr = fit_pdbentry_ptr->next_pdbentry_ptr;
  first = TRUE;

  /* Process the first of the other PDB entries */
  if (pdbentry_ptr != NULL)
    {
      /* Create a set of residue sequence records for the alignment */
/* v.2.0 --> */
      if (params.dimplot == TRUE)
        nresid[REF_SEQUENCE]
          = get_dimplot_sequence(pdbentry_ptr,&(sequence[REF_SEQUENCE]),
                                 &(residue_index_ptr[REF_SEQUENCE]),in_plot);
      else
/* <-- v.2.0 */
        nresid[REF_SEQUENCE]
          = get_bs_sequence(pdbentry_ptr,&(sequence[REF_SEQUENCE]),
                            &(residue_index_ptr[REF_SEQUENCE]),
                            dummy_residue_ptr,&nchains);

      /* If we have more than one chain, create a second copy
         of the sequence with the chains reversed to see which gives a
         better fit */
      if (nchains > 1)
        {
          /* Take a copy of the current sequence */
          copy_sequence(residue_index_ptr[CURR_SEQUENCE],
                        nresid[CURR_SEQUENCE],
                        &(residue_index_ptr[COPY_SEQUENCE]));
          nresid[COPY_SEQUENCE] = nresid[CURR_SEQUENCE];
          strcpy(sequence[COPY_SEQUENCE],sequence[CURR_SEQUENCE]);

          /* Create a copy of the reference sequence */
          copy_sequence(residue_index_ptr[REF_SEQUENCE],nresid[REF_SEQUENCE],
                        &(residue_index_ptr[REVERSE_SEQUENCE]));
          nresid[REVERSE_SEQUENCE] = nresid[REF_SEQUENCE];
          strcpy(sequence[REVERSE_SEQUENCE],sequence[REF_SEQUENCE]);

          /* Alter sequence of residues so that the chains are in reverse
             order */
          reverse_chains(residue_index_ptr[REF_SEQUENCE],
                         residue_index_ptr[REVERSE_SEQUENCE],
                         nresid[REF_SEQUENCE]);

          /* Need to loop three times */
          nloops = 3;
        }
      else
        nloops = 1;

      /* Generate trial alignments of binding site residues */
      for (loop = 0; loop < nloops; loop++)
        {
          /* Align reference sequence against current sequence */
          use = 0;

          /* If this is the second loop, use the reversed-chain
             sequence for the alignment*/
          if (loop == 1)
            use = 2;

          /* Create the arrays for the sequence alignment */
          nmap = create_alignment_arrays(nresid,&array,&tmp_array,
                                         &path,&mapping);

          /* Perform the alignment */
          nfit = perform_alignment(&sequence[use],nresid,array,tmp_array,
/* v.2.0 --> */
//                                   path,mapping,nmap,&percentage);
                                   path,mapping,nmap,&percentage,
/* v.2.0.1 --> */
//                                   residue_index_ptr);
                                   residue_index_ptr,params);
/* <-- v.2.0.1 */
/* <-- v.2.0 */

          /* If this is the best fit so far, store its details */
          if (nfit > MIN_ALIGNED && percentage > best_percentage)
            {
              /* Store the details */
              best_fit = nfit;
              best_percentage = percentage;

              /* If this is the second loop, then can use this alignment */
              if (loop == 1)
                nloops = 2;
            }
        }
    }

  /* If we have an alignment, then get the residue mappings */
  if (best_fit > MIN_ALIGNED)
    {
      /* Extract the residue mappings from the alignment */
      nfit = get_residue_mappings(residue_index_ptr,mapping,nmap,
/* v.2.0 --> */
//                                  &alignment);
/* v.2.0.1 --> */
//                                  &alignment,params.antibody);
                                  &alignment,params);
/* <-- v.2.0.1 */
/* <-- v.2.0 */
    }
  else
    nfit = 0;

  /* Free up memory */
  free(mapping);

  /* Return the alignment */
  *aln = alignment;

  /* Return the aligned residues in the best alignment */
  return(nfit);
}
/***********************************************************************

fit_structures  -  Use the equivalent residues to fit the structures

***********************************************************************/

int fit_structures(struct resseq **alignment,int align_len,
/* v.2.0.1 --> */
//                   double cofg[2][3],REAL matrix[3][3])
                   double cofg[2][3],REAL matrix[3][3],
                   struct parameters params)
/* <-- v.2.0.1 */
{
  int i, j, loop, nfit, npaired, ntried;
  int done;

  double min_sep_cutoff;

  /* Initialise variables */
  done = FALSE;
  loop = 0;
  nfit = 0;

  /* Initialise transformation matrices */
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      if (i == j)
        matrix[i][j] = 1.0;
      else
        matrix[i][j] = 0.0;
  for (i = 0; i < 2; i++)
    for (j = 0; j < 3; j++)
      cofg[i][j] = 0.0;

  /* Loop until reasonable fit achieved */
  while (loop < NCUTOFFS && done == FALSE)
    {
      /* Get the distance cutoff to be used in this loop */
      min_sep_cutoff = MIN_SEP_CUTOFF[loop];

      /* Use the alignments to fit the 3D coords */
      npaired
        = fit_coords(alignment,align_len,matrix,cofg,&ntried,
/* v.2.0.1 --> */
//                     min_sep_cutoff);
                     min_sep_cutoff,params);
/* <-- v.2.0.1 */

      /* If fit failed, or the number of potential pairs equals
         the number tried, then we are done */
      if (npaired == 0 || npaired == ntried)
        done = TRUE;

      /* Reduce the separation cut-off for next loop (if there is
         one) */
      loop++;
    }

  /* If fit looks OK, then save the fit type */
  nfit = npaired;

  /* Return number of successfully fitted residues */
  return(nfit);
}
/* <-- v.1.3 */
/***********************************************************************

create_sequence  -  Create a set of residue sequence records for the
                    alignment

***********************************************************************/

int create_sequence(struct pdbentry *pdbentry_ptr,
                    struct resseq **fst_resseq_ptr,
/* v.1.0.4 --> */
/*                    struct residue *dummy_residue_ptr) */
                    struct residue *dummy_residue_ptr,int *nc)
/* <-- v.1.0.4 */
{
  char insertion, last_chain, last_insertion, number_string[20];

  int last_res_no, nresid, res_no;
/* v.1.0.4 --> */
  int nchains;
/* <-- v.1.0.4 */

  struct resseq *resseq_ptr, *first_resseq_ptr, *last_resseq_ptr;
  struct residue *residue_ptr;

  /* Initialise */
  last_res_no = 0;
  last_chain = '\0';
  last_insertion = '\0';
/* v.1.0.4 --> */
  nchains = 0;
/* <-- v.1.0.4 */
  nresid = 0;

  /* Initialise pointers */
  resseq_ptr = first_resseq_ptr = last_resseq_ptr = NULL;
  *fst_resseq_ptr = NULL;

  /* Get the first of this PDB entry's residues */
  residue_ptr = pdbentry_ptr->first_residue_ptr;

  /* Loop over the residues to create entries in the sequence for all
     interacting protein residues */
  while (residue_ptr != NULL)
    {
      /* If this is an interacting protein residue, then create an
         entry for it in the sequence */
      if ((residue_ptr->interaction == NONBOND_TO_LIGAND ||
          residue_ptr->interaction == HBOND_TO_LIGAND ||
          residue_ptr->interaction == COVALENT_TO_LIGAND) &&
          residue_ptr->type != DUMMY &&
          residue_ptr->aa_code != 'X')
        {
          /* Get the numeric residue number */
          strcpy(number_string,residue_ptr->res_num);
          number_string[4] = '\0';
          res_no = atoi(number_string);

          /* Get the insertion code */
          insertion = residue_ptr->res_num[4];

/* v.1.0.4 --> */
          /* If this is a new chain, increment count of chains */
          if (residue_ptr->chain != last_chain)
            nchains++;
/* <-- v.1.0.4 */

          /* If this residue doesn't follow directly after the previous
             one, and is not the first, add in a dummy residue record */
          if ((residue_ptr->chain != last_chain ||
               (residue_ptr->chain == last_chain &&
                res_no > last_res_no + 1)) &&
              last_chain != '\0')
            {
              /* Add a dummy residue to this sequence */
              resseq_ptr
                = create_resseq_record(&first_resseq_ptr,&last_resseq_ptr,
                                       dummy_residue_ptr);
                          
              /* Increment count of residues read in */
              nresid++;
            }

          /* Add the current residue to this sequence */
          resseq_ptr
            = create_resseq_record(&first_resseq_ptr,&last_resseq_ptr,
                                   residue_ptr);
          
          /* Increment count of residues read in */
          nresid++;

          /* Save residue number, chain and insertion code */
          last_res_no = res_no;
          last_insertion = insertion;
          last_chain = residue_ptr->chain;
        }

      /* Get the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
    }

  /* Return pointer to first resseq record*/
  *fst_resseq_ptr = first_resseq_ptr;

/* v.1.0.4 --> */
  /* Return number of chains encountered */
  *nc = nchains;
/* <-- v.1.0.4 */

  /* Return number of residue records stored */
  return(nresid);
}
/***********************************************************************

get_residue_type  -  Get the residue type from its name

***********************************************************************/

int get_residue_type(char *res_nam,char residue_name[NAMINO+1][4])
{
#define NNONSTAND  5

  char res_name[4];

  int inon, ires, type;

  static char *non_standard[] = {
    "CSS", "MIL", "MSE", "PTR", "SEP"
  };
  static char *replacement[] = {
    "CYS", "MET", "MET", "TYR", "SER"
  };

  /* Initialise variables */
  type = 0;
  strcpy(res_name,res_nam);

  /* Loop through all the residues until find a match */
  for (ires = 1; ires < NAMINO + 1 && type == 0; ires++)
    if (!strcmp(res_name,residue_name[ires]))
      type = ires;

  /* If not a standard amino acid, check for any of the common
     modified forms */
  if (type == 0)
    {
      /* Set the non-standard identifier */
      inon = -1;

      /* Scan all the non-standard residue names */
      for (ires = 0; ires < NNONSTAND && inon == -1; ires++)
        if (!strcmp(res_name,non_standard[ires]))
          inon = ires;

      /* If have a non-standard name, replace by standard version
         and search the standard codes again */
      if (inon != -1)
        {
          /* Get standard a.a. code */
          strcpy(res_name,replacement[inon]);

          /* Search against the standard amino acid codes again */
          for (ires = 1; ires < NAMINO + 1 && type == 0; ires++)
            if (!strcmp(res_name,residue_name[ires]))
              type = ires;
        }
    }

  /* Return the residue type */
  return(type);
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

create_dummy_dic  -  Create a dummy Het Group Dictionary using the
                     fragment read in from the .sdf file

***********************************************************************/

void create_dummy_dic(struct residue *first_residue_ptr,
                      struct parameters params)
{
  char bond_order[5];
  char atom_name[10], stored_name[MAX_TYPES][5];

  char *res_name = "LIG";

  int atom_number, bond_number, i, iatom, natoms;
  int atom_type, count, name_count[MAX_TYPES], natom_types;

  struct atom *atom_ptr;
  struct bond *bond_ptr;
  struct residue *residue_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  atom_number = 0;
  bond_number = 0;
  for (i = 0; i < MAX_TYPES; i++)
    name_count[i] = 0;
  natom_types = 0;

  /* Open the dummy dictionary file */
  file_ptr = open_file("dummy.cif",params);

  /* Write the header records */
  fprintf(file_ptr,"data_%s\n",res_name);
  fprintf(file_ptr,"# \n");
  fprintf(file_ptr,"_chem_comp.id                                    "
         "%s \n",res_name);
  fprintf(file_ptr,"_chem_comp.name                                  "
         "Ligand\n");
  fprintf(file_ptr,"# \n");
  fprintf(file_ptr,"loop_\n");
  fprintf(file_ptr,"_chem_comp_atom.comp_id \n");
  fprintf(file_ptr,"_chem_comp_atom.atom_id \n");
  fprintf(file_ptr,"_chem_comp_atom.alt_atom_id \n");
  fprintf(file_ptr,"_chem_comp_atom.type_symbol \n");
  fprintf(file_ptr,"_chem_comp_atom.charge \n");
  fprintf(file_ptr,"_chem_comp_atom.pdbx_align \n");
  fprintf(file_ptr,"_chem_comp_atom.pdbx_aromatic_flag \n");
  fprintf(file_ptr,"_chem_comp_atom.pdbx_leaving_atom_flag \n");
  fprintf(file_ptr,"_chem_comp_atom.pdbx_stereo_config \n");
  fprintf(file_ptr,"_chem_comp_atom.model_Cartn_x \n");
  fprintf(file_ptr,"_chem_comp_atom.model_Cartn_y \n");
  fprintf(file_ptr,"_chem_comp_atom.model_Cartn_z \n");
  fprintf(file_ptr,"_chem_comp_atom.pdbx_model_Cartn_x_ideal \n");
  fprintf(file_ptr,"_chem_comp_atom.pdbx_model_Cartn_y_ideal \n");
  fprintf(file_ptr,"_chem_comp_atom.pdbx_model_Cartn_z_ideal \n");
  fprintf(file_ptr,"_chem_comp_atom.pdbx_ordinal \n");

  /* Get pointer to the first residue */
  residue_ptr = first_residue_ptr;

  /* Loop over the ligand residues to write out */
  while (residue_ptr != NULL)
    {
      /* Only write out if this is a ligand residue */
      if (residue_ptr->interaction == IS_LIGAND)
        {
          /* Get the first of this residue's atoms */
          atom_ptr = residue_ptr->first_atom_ptr;
          natoms = residue_ptr->natoms;
          iatom = 0;

          /* Loop over the residue's atoms to calculate centre of mass */
          while (atom_ptr != NULL && iatom < natoms)
            {
              /* Increment atom number */
              atom_number++;

              /* Generate a fake atom name based on the element number
                 (so that multiple-residue ligands are treated as a single
                 Het Group) */
              atom_type = get_attype(atom_ptr->element,stored_name,
                                     name_count,&natom_types);
              count = name_count[atom_type];

              /* Form this atom's name */
              if (strlen(atom_ptr->element) == 2)
                sprintf(atom_name,"%s%d  ",atom_ptr->element,count);
              else
                sprintf(atom_name," %s%d  ",atom_ptr->element,count);
              atom_name[4] = '\0';

              /* Store the fake name to identify the atom later */
              strcpy(atom_ptr->fake_name,atom_name);

              /* Write out this atom */
              if (count < 100)
                fprintf(file_ptr,"%s %s %s%s 0 1 N N N %8.3f %8.3f %8.3f "
                        "%8.3f %8.3f %8.3f %d\n",res_name,atom_name,
                        atom_name,atom_ptr->element,
                        atom_ptr->original_coords[0],
                        atom_ptr->original_coords[1],
                        atom_ptr->original_coords[2],
                        atom_ptr->original_coords[0],
                        atom_ptr->original_coords[1],
                        atom_ptr->original_coords[2],
                        atom_number);

              /* Get the next atom */
              atom_ptr = atom_ptr->next_atom_ptr;
              iatom++;
            }
        }

      /* Get the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
    }

  /* Start of the bond's section */
  fprintf(file_ptr,"# \n");
  fprintf(file_ptr,"loop_\n");
  fprintf(file_ptr,"_chem_comp_bond.comp_id \n");
  fprintf(file_ptr,"_chem_comp_bond.atom_id_1 \n");
  fprintf(file_ptr,"_chem_comp_bond.atom_id_2 \n");
  fprintf(file_ptr,"_chem_comp_bond.value_order \n");
  fprintf(file_ptr,"_chem_comp_bond.pdbx_aromatic_flag \n");
  fprintf(file_ptr,"_chem_comp_bond.pdbx_stereo_config \n");
  fprintf(file_ptr,"_chem_comp_bond.pdbx_ordinal \n");

  /* Get the current fragment's first atom */
  bond_ptr = first_residue_ptr->pdbentry_ptr->first_bond_ptr;

  /* Loop over the fragment's bonds to write out */
  while (bond_ptr != NULL)
    {
      /* Check that both this bond's atoms are from the ligand residue */
      if (bond_ptr->atom1_ptr->residue_ptr->interaction == IS_LIGAND &&
/* v.2.2.1 --> */
//          bond_ptr->atom2_ptr->residue_ptr->interaction == IS_LIGAND)
          bond_ptr->atom2_ptr->residue_ptr->interaction == IS_LIGAND &&
          bond_ptr->type == COVALENT)
/* <-- v.2.2.1 */
        {
          /* Determine bond order */
          strcpy(bond_order,"SING");
/* v.2.2.1 --> */
//          if (bond_ptr->type > 1 && bond_ptr->type <= NBOND_TYPES)
//            strcpy(bond_order,BOND_ORDER[bond_ptr->type]);
          if (bond_ptr->bond_order > 1 && bond_ptr->bond_order < NBOND_TYPES)
            strcpy(bond_order,BOND_ORDER[bond_ptr->bond_order]);
/* <-- v.2.2.1 */

          /* Increment bond counter */
          bond_number++;

          /* Write out this bond */
          fprintf(file_ptr,"%s %s %s %s N N %d\n",res_name,
                  bond_ptr->atom1_ptr->fake_name,
                  bond_ptr->atom2_ptr->fake_name,
                  bond_order,bond_number);
        }

      /* Get the next bond */
      bond_ptr = bond_ptr->next_bond_ptr;
    }

  /* Close the output file */
  fclose(file_ptr);
}
/***********************************************************************

write_dummy_pdb  -  Write out the current ligand to a PDB file for input
                    to hbadd

***********************************************************************/

void write_dummy_pdb(struct residue *first_ligand_residue_ptr,
                     struct parameters params)
{
  char atom_name[10], stored_name[MAX_TYPES][5];

  char *res_name = "LIG";

  int atom_number, i, iatom, natoms;
  int atom_type, count, name_count[MAX_TYPES], natom_types;

  struct atom *atom_ptr;
  struct residue *residue_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  atom_number = 0;
  for (i = 0; i < MAX_TYPES; i++)
    name_count[i] = 0;
  natom_types = 0;

  /* Open the ligand coordinate file, dummy.pdb */
  file_ptr = open_file("dummy.pdb",params);

  /* Get the first ligand residue */
  residue_ptr = first_ligand_residue_ptr;

  /* Loop over the ligand residues */
  while (residue_ptr != NULL && residue_ptr->interaction == IS_LIGAND)
    {
      /* Get the residue's first atom */
      atom_ptr = residue_ptr->first_atom_ptr;
      natoms = residue_ptr->natoms;
      iatom = 0;

      /* Loop over the residue's atoms to write out ligand coords */
      while (atom_ptr != NULL && iatom < natoms)
        {
          /* Increment atom number */
          atom_number++;

          /* Generate a fake atom name based on the element number
             (so that multiple-residue ligands are treated as a single
             Het Group) */
          atom_type = get_attype(atom_ptr->element,stored_name,
                                 name_count,&natom_types);
          count = name_count[atom_type];
          
          /* Form this atom's name */
          if (strlen(atom_ptr->element) == 2)
            sprintf(atom_name,"%s%d  ",atom_ptr->element,count);
          else
            sprintf(atom_name," %s%d  ",atom_ptr->element,count);
          atom_name[4] = '\0';

          /* Store the fake name to identify the atom later */
          strcpy(atom_ptr->fake_name,atom_name);

          /* Write out */
          if (count < 100)
            fprintf(file_ptr,"HETATM%5d %s %s %c%s   "
                    "%8.3f%8.3f%8.3f%6.2f%6.2f          %s\n",
                    atom_ptr->atom_number,atom_name,res_name,
                    residue_ptr->chain,
                    residue_ptr->res_num,
                    atom_ptr->original_coords[X],
                    atom_ptr->original_coords[Y],
                    atom_ptr->original_coords[Z],
                    atom_ptr->occupancy,
                    atom_ptr->bvalue,atom_ptr->element);

          /* Get the next atom record and increment atom count */
          atom_ptr = atom_ptr->next_atom_ptr;
          iatom++;
        }

      /* Get the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
    }
  
  /* Close the PDB file */
  fclose(file_ptr);
}
/* v.1.4.1 --> */
/***********************************************************************

copy_file  -  Copy the appropriate file

***********************************************************************/

void copy_file(struct parameters params,char *infile,char *outfile)
{
  char infile_name[FILENAME_LEN + 1], outfile_name[FILENAME_LEN + 1];
  char input_line[LINELEN + 1];

  FILE *infile_ptr, *outfile_ptr;

  /* Form name of the input file */
  get_file_name(params.directory,infile,infile_name);

  /* Form name of the output file */
  get_file_name(params.directory,outfile,outfile_name);

  /* Open output file */
  open_output_file(&outfile_ptr,outfile_name);

  /* Open the input file */
  if ((infile_ptr = fopen(infile_name,"r")) !=  NULL)
    {
      /* Read through the file, writing each line out */
      while (fgets(input_line,LINELEN,infile_ptr) != NULL)
        {
          /* Write to output file */
          fprintf(outfile_ptr,"%s",input_line);
        }
    }

  /* Close both files */
  fclose(infile_ptr);
  fclose(outfile_ptr);
}
/* <-- v.1.4.1 */
/***********************************************************************

run_hbadd  -  Run HBADD to map the two ligands

***********************************************************************/

void run_hbadd(struct parameters params)
{
  char command[FILENAME_LEN + 1], file_name[FILENAME_LEN + 1];
  char hbadd_exe[FILENAME_LEN + 1], cif_file[FILENAME_LEN + 1];

  //  char *res_name = "LIG";

  /* Form name of the hbadd executable */
  get_file_name(params.ligplot_exe_dir,"hbadd",hbadd_exe);

  /* Form name of the PDB file */
  get_file_name(params.directory,"dummy.pdb",file_name);

  /* Form name of the dummy.cif dictionary file */
  get_file_name(params.directory,"dummy.cif",cif_file);

  /* Form the parameters to run HBADD */
  sprintf(command,"\"%s\" \"%s\" \"%s\" -wkdir \"%s\"",hbadd_exe,
          file_name,cif_file,params.directory);

  /* Write out the HBADD command */
  printf("Ligfit's HBADD command: %s\n",command);
  fflush(stdout);

  /* Run HBADD */
  printf("HBADD command: %s\n",command);
  system(command);
}
/***********************************************************************

read_hbadd_map  -  Read in the mapping between the two sets of ligand
                   atoms from the hbadd.map file

***********************************************************************/

int read_hbadd_map(struct residue **first_ligand_residue_ptr,
                   struct atomap **fst_atomap_ptr,
                   struct parameters params)
{
  char file_name[FILENAME_LEN + 1], input_line[LINELEN + 1];
  char atom_name[2][5];

  char *res_name = "LIG";

  int iatom, len, loop, natoms, nfound, nmatched;
/* v.1.0.4 --> */
/*  int done, found, ok; */
  int done, found, ok;
/* <-- v.1.0.4 */

  struct atom *atom_ptr, *found_atom_ptr[2];
  struct atomap *atomap_ptr, *first_atomap_ptr, *last_atomap_ptr;
  struct residue *residue_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  done = FALSE;
  nfound = nmatched = 0;
  ok = TRUE;

  /* Initialise pointers */
  first_atomap_ptr = last_atomap_ptr = NULL;

  /* Form file name for the hbadd.map file */
  get_file_name(params.directory,"hbadd.map",file_name);

  /* Open the input file */
  printf("Opening %s ...\n",file_name);
  if ((file_ptr = fopen(file_name,"r")) !=  NULL)
    {
      /* Loop while reading in records from the file */
      while (fgets(input_line,LINELEN,file_ptr) != NULL && done == FALSE)
        {
          /* If line wanted, process it */
          if ((input_line[0] == '[' &&
               !strncmp(input_line + 1,res_name,3)) &&
              (strstr(input_line,"Mapped") != NULL ||
               strchr(input_line,'*') != NULL))
            {
              /* Truncate the string to strip off the newline */
              len = string_truncate(input_line,LINELEN);
              nmatched++;

              /* Get the two atom names */
              strncpy(atom_name[0],input_line+9,4);
              atom_name[0][4] = '\0';
              strncpy(atom_name[1],input_line+21,4);
              atom_name[1][4] = '\0';

              /* Loop through the two ligands to match the two atoms */
              for (loop = 0; loop < 2; loop++)
                {
                  /* Initialise flag */
                  found = FALSE;
                  found_atom_ptr[loop] = NULL;

                  /* Get this entry's first residue */
                  residue_ptr = first_ligand_residue_ptr[loop];

                  /* Loop over this ligand's residues */
                  while (residue_ptr != NULL &&
                         residue_ptr->interaction == IS_LIGAND &&
                         found == FALSE)
                    {
                      /* Get the residue's first atom */
                      atom_ptr = residue_ptr->first_atom_ptr;
                      natoms = residue_ptr->natoms;
                      iatom = 0;

                      /* Loop over the residue's atoms to find 
                         one read in from file */
                      while (atom_ptr != NULL && iatom < natoms &&
                             found == FALSE)
                        {
                          /* If this is the right atom, store */
                          if (!strcmp(atom_ptr->fake_name,atom_name[loop]))
                            {
                              found_atom_ptr[loop] = atom_ptr;
                              found = TRUE;
                            }

                          /* Get the next atom record and increment
                             atom count */
                          atom_ptr = atom_ptr->next_atom_ptr;
                          iatom++;
                        }

                      /* Get pointer to next residue record */
                      residue_ptr = residue_ptr->next_residue_ptr;
                    }
                }
            
              /* If both atoms found, create an atomap record to map
                 them to one another */
              if (found_atom_ptr[0] != NULL &&
                  found_atom_ptr[1] != NULL)
                {
                  /* Create the atom mapping */
                  atomap_ptr
                    = create_atomap_record(&first_atomap_ptr,
                                           &last_atomap_ptr,
                                           found_atom_ptr[0],
                                           found_atom_ptr[1],0);

                  /* Increment count of identified atoms */
                  nfound++;
                }
            }

          /* If end of residue reached, we're done */
          else if (input_line[0] == 'A' ||
                   ((input_line[0] == '[' &&
                     strncmp(input_line + 1,res_name,3)) &&
                    nmatched > 0))
            done = TRUE;
        }

      /* Close the file */
      fclose(file_ptr);
    }

  /* File not found, show error */
  else
    {
      printf("*** ERROR. FILE NOT FOUND: %s ...\n",file_name);
      ok = FALSE;
    }

  /* Show number of atoms matched */
  printf("Number of atoms matched = %d. Number found: %d\n",nmatched,
         nfound);

  /* If too few matched found, then don't use this mapping */
  if (nfound < MIN_ATOMS_MAPPED)
    ok = FALSE;

  /* Return the first atom-map pointer */
  *fst_atomap_ptr = first_atomap_ptr;

  /* Return if file read OK */
  return(ok);
}
/* v.1.3 --> */
/***********************************************************************

mark_ligands  -  Mark which ligand residues map to one another

***********************************************************************/

void mark_ligands(struct atomap *first_atomap_ptr)
{
  struct atom *atom1_ptr, *atom2_ptr;
  struct atomap *atomap_ptr;
  struct residue *residue1_ptr, *residue2_ptr;

  /* Get pointer to the first atomap */
  atomap_ptr = first_atomap_ptr;

  /* Loop over the atomap records */
  while (atomap_ptr != NULL)
    {
      /* Get the two atoms */
      atom1_ptr = atomap_ptr->atom_ptr;
      atom2_ptr = atomap_ptr->map_atom_ptr;

      /* Get their residues */
      residue1_ptr = atom1_ptr->residue_ptr;
      residue2_ptr = atom2_ptr->residue_ptr;

      /* If not already mapped, perform the mapping */
      if (residue1_ptr->map_residue_ptr == NULL &&
          residue2_ptr->map_residue_ptr == NULL)
        {
          residue1_ptr->map_residue_ptr = residue2_ptr;
          residue2_ptr->map_residue_ptr = residue1_ptr;
        }

      /* Get pointer to next atomap record */
      atomap_ptr = atomap_ptr->next_atomap_ptr;
    }
}
/* <-- v.1.3 */
/***********************************************************************

fit_ligands  -  Fit on the matched ligand atoms

***********************************************************************/

int fit_ligands(struct atomap *first_atomap_ptr,REAL matrix[3][3],
                double cofg[2][3])
{
  BOOL column;

  COOR coords1[MAX_FIT_RES], coords2[MAX_FIT_RES];

  int i, j, npaired, nstored;
  int ok;

  double rmsd, rms_cutoff;

  struct atom *atom1_ptr, *atom2_ptr;
  struct atomap *atomap_ptr;

  /* Initialise variables */
  ok = TRUE;
  column = FALSE;
  npaired = 0;
  nstored = 0;
  rmsd = 0.0;
  rms_cutoff = RMSD_CUTOFF;

  /* Initialise transformation matrices */
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      if (i == j)
        matrix[i][j] = 1.0;
      else
        matrix[i][j] = 0.0;
  for (i = 0; i < 2; i++)
    for (j = 0; j < 3; j++)
      cofg[i][j] = 0.0;

  /* Get pointer to the first atomap */
  atomap_ptr = first_atomap_ptr;

  /* Loop over the atomap records */
  while (atomap_ptr != NULL)
    {
      /* Get the two atoms */
      atom1_ptr = atomap_ptr->atom_ptr;
      atom2_ptr = atomap_ptr->map_atom_ptr;

      /* Store the coords */
      coords1[nstored].x = atom1_ptr->original_coords[X];
      coords1[nstored].y = atom1_ptr->original_coords[Y];
      coords1[nstored].z = atom1_ptr->original_coords[Z];
      coords2[nstored].x = atom2_ptr->original_coords[X];
      coords2[nstored].y = atom2_ptr->original_coords[Y];
      coords2[nstored].z = atom2_ptr->original_coords[Z];

      /* Increment count of stored coords */
      nstored++;

      /* Get pointer to next atomap record */
      atomap_ptr = atomap_ptr->next_atomap_ptr;
    }

  /* Centre coords of first ligand */
  centre_coords(coords1,nstored,cofg[CURR_SEQUENCE]);

  /* Repeat for second ligand */
  centre_coords(coords2,nstored,cofg[REF_SEQUENCE]);

  /* If have more than one pair of atoms stored, can perform the fit */
  if (nstored > 1)
    {
      /* Perform the least-squares fit */
      rmsd = matfit(coords1,coords2,matrix,nstored,NULL,column);

/* v.1.3 --> */
      /* Mark the ligand residues as mapped to one another */
      mark_ligands(first_atomap_ptr);
/* <-- v.1.3 */
    }

  /* Otherwise, unset flag */
  else
    ok = FALSE;

  /* Return whether OK */
  return(ok);
}
/***********************************************************************

perform_graph_match  -  Run hbadd to map the ligand atoms

***********************************************************************/

int perform_graph_match(struct residue **first_ligand_residue_ptr,
                        int *nresid,REAL matrix[3][3],double cofg[2][3],
                        struct atomap **fst_atomap_ptr,
                        struct parameters params)
{
  int ok;

  struct atomap *first_atomap_ptr;

  /* Initialise variables */
  *fst_atomap_ptr = first_atomap_ptr = NULL;

  /* Create a dummy Het Group Dictionary using the ligand from the
     earlier LIGPLOT run */
  create_dummy_dic(first_ligand_residue_ptr[1],params);

  /* Write out the current ligand to a PDB file for input to hbadd */
  write_dummy_pdb(first_ligand_residue_ptr[0],params);

/* v.1.4.1 --> */
  /* Save the hbadd.bonds file */
  copy_file(params,"hbadd.bonds","tmp.bonds");
/* <-- v.1.4.1 */

  /* Run HBADD to map the fragment to the current ligand */
  run_hbadd(params);

  /* Read in the mapping between the fragment and ligand atoms from
     the hbadd.map file */
  ok = read_hbadd_map(first_ligand_residue_ptr,&first_atomap_ptr,
                      params);

  /* If a good mapping found, fit the two structures based on the atom
     equivalences obtained */
  if (ok == TRUE)
    {
      /* Fit on the matched ligand atoms */
      ok = fit_ligands(first_atomap_ptr,matrix,cofg);
    }

/* v.1.4.1 --> */
  /* Restore the previous hbadd.bonds file */
  copy_file(params,"tmp.bonds","hbadd.bonds");
/* <-- v.1.4.1 */

  /* Return first of the atom matches */
  *fst_atomap_ptr = first_atomap_ptr;

  /* Return whether fit was accomplished */
  return(ok);
}
/***********************************************************************

match_ligands  -  Fit the ligands by calling hbadd 

***********************************************************************/

int match_ligands(struct pdbentry *first_pdbentry_ptr,
                  REAL matrix[3][3],double cofg[2][3],
                  struct atomap **fst_atomap_ptr,
                  struct parameters params)
{
  int loop, nresid[2];
  int ok;

  struct atomap *first_atomap_ptr;
  struct pdbentry *pdbentry_ptr;
  struct residue *residue_ptr;
  struct residue *first_ligand_residue_ptr[2];

  /* Initialise variables */
  printf("Matching ligands ...\n");
  ok = FALSE;
  first_ligand_residue_ptr[0] = first_ligand_residue_ptr[1] = NULL;
  nresid[0] = nresid[1] = 0;
  *fst_atomap_ptr = first_atomap_ptr = NULL;

  /* Get the first (fit) PDB entry */
  pdbentry_ptr = first_pdbentry_ptr;

  /* Loop to pick up the first ligand residue of each entry */
  for (loop = 0; loop < 2 && pdbentry_ptr != NULL; loop++)
    {
      /* Get the first of the PDB entry's residues */
      residue_ptr = pdbentry_ptr->first_residue_ptr;

      /* Loop over the residues to find the first ligand residue */
      while (residue_ptr != NULL)
        {
          /* If this is a ligand residue, then store pointer */
          if (residue_ptr->interaction == IS_LIGAND)
            {
              /* If first pointer null, then store */
              if (first_ligand_residue_ptr[loop] == NULL)
                first_ligand_residue_ptr[loop] = residue_ptr;

              /* Increment count of ligand residues */
              nresid[loop]++;
            }

          /* Get the next residue */
          residue_ptr = residue_ptr->next_residue_ptr;
        }

      /* Get the next PDB entry */
      pdbentry_ptr = pdbentry_ptr->next_pdbentry_ptr;
    }

  /* If either ligand-start is missing, abort */
  if (first_ligand_residue_ptr[0] == NULL ||
      first_ligand_residue_ptr[1] == NULL)
    return(ok);

  /* Graph match the ligands */
  ok = perform_graph_match(first_ligand_residue_ptr,nresid,matrix,cofg,
                           &first_atomap_ptr,params);

  /* Return first of the atom matches */
  *fst_atomap_ptr = first_atomap_ptr;

  /* Return whether match and fitting was successful */
  return(ok);
}
/* v.1.3 --> */
/***********************************************************************

weigh_anchors  -  Upweight equivalent atoms where they make H-bonds to
                  equivalent protein residues 

***********************************************************************/

void weigh_anchors(struct atomap *first_atomap_ptr)
{
  int equiv;

  struct atomap *atomap_ptr;

  /* Get the first of the atom mappings */
  atomap_ptr = first_atomap_ptr;

  /* Loop to look for equivalent interactions */
  while (atomap_ptr != NULL)
    {
      /* Determine whether the two atoms both have H-bonds
         to equivalenced residues */
      equiv = check_equiv_hbonds(atomap_ptr->atom_ptr,
                                 atomap_ptr->map_atom_ptr);
                      
      /* If so, then increase anchor weight between them */
      if (equiv == TRUE)
        {
          atomap_ptr->anchor_weight = ANCHOR_WEIGHT;
          atomap_ptr->atom_ptr->anchor_weight = ANCHOR_WEIGHT;
        }

      /* Get pointer to next record */
      atomap_ptr = atomap_ptr->next_atomap_ptr;
    }
}
/* <-- v.1.3 */
/***********************************************************************

write_coords  -  Write out the transformed coords to the ligplus.pdb
                 file

***********************************************************************/

void write_coords(struct pdbentry *pdbentry_ptr,
                  REAL matrix[3][3],double cofg[2][3],
                  struct parameters params)
{
/* v.2.0 --> */
  char last_chain_type;
/* <-- v.2.0 */

  int iatom, natoms, nout;
/* v.2.0.3 --> */
  int min_atno, max_atno;
/* <-- v.2.0.3 */

  double x, y, z;

  struct atom *atom_ptr;
/* v.2.0.3 --> */
  struct bond *bond_ptr;
/* <-- v.2.0.3 */
/* v.1.1 --> */
  struct hhbnnb *hhbnnb_ptr;
/* <-- v.1.1 */
  struct residue *residue_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  nout = 0;
/* v.2.0 --> */
  last_chain_type = '\0';
/* <-- v.2.0 */
/* v.2.0.3 --> */
  min_atno = max_atno = -1;
/* <-- v.2.0.3 */

  /* Open the output ligplus.pdb file */
  file_ptr = open_file("ligplus.pdb",params);

  /* Write out the header record */
  fprintf(file_ptr,"REMARK PDB CODE: %s\n",pdbentry_ptr->pdb_code);

/* v.1.1 --> */
  /* If have any HHB or NNB records, write these out */
  hhbnnb_ptr = pdbentry_ptr->first_hhbnnb_ptr;
  while (hhbnnb_ptr != NULL)
    {
      /* Write out this record */
      fprintf(file_ptr,"%s\n",hhbnnb_ptr->data);

      /* Get the next record */
      hhbnnb_ptr = hhbnnb_ptr->next_hhbnnb_ptr;
    }
/* <-- v.1.1 */

  /* Get the first of the PDB entry's residues */
  residue_ptr = pdbentry_ptr->first_residue_ptr;

  /* Loop over the residues to write out the transformed coords */
  while (residue_ptr != NULL)
    {
/* v.2.0 --> */
      /* If residue belongs to a new domain, write out a domain header */
      if (residue_ptr->chain_type != last_chain_type)
        {
          fprintf(file_ptr,"REMARK DOMAIN %d\n",residue_ptr->chain_type);
          last_chain_type = residue_ptr->chain_type;
        }
/* <-- v.2.0 */

      /* Get the first atom */
      atom_ptr = residue_ptr->first_atom_ptr;
      natoms = residue_ptr->natoms;
      iatom = 0;

      /* Loop over the residue's atoms to write out transformed coords */
      while (atom_ptr != NULL && iatom < natoms)
        {
          /* Transform the coords */
          transform_atom(atom_ptr->original_coords[X],
                         atom_ptr->original_coords[Y],
                         atom_ptr->original_coords[Z],matrix,cofg,
                         &x,&y,&z);

          /* Write out */
          fprintf(file_ptr,"ATOM  %5d %s %s %c%s   "
                  "%8.3f%8.3f%8.3f%6.2f%6.2f          %s\n",
                  atom_ptr->atom_number,atom_ptr->atom_name,
                  residue_ptr->res_name,
                  residue_ptr->chain,
                  residue_ptr->res_num,x,y,z,atom_ptr->occupancy,
                  atom_ptr->bvalue,atom_ptr->element);

/* v.2.0.3 --> */
          /* If this a ligand atom, update the maximum and minimum
             atom range */
          if (residue_ptr->interaction == IS_LIGAND)
            {
              /* Update ligand-atom range */
              if (min_atno == -1)
                min_atno = max_atno = atom_ptr->atom_number;
              else
                {
                  if (atom_ptr->atom_number < min_atno)
                    min_atno = atom_ptr->atom_number;
                  if (atom_ptr->atom_number > max_atno)
                    max_atno = atom_ptr->atom_number;
                }
            }
/* <-- v.2.0.3 */

          /* Increment count of atoms written out */
          nout++;

          /* Get the next atom record and increment atom count */
          atom_ptr = atom_ptr->next_atom_ptr;
          iatom++;
        }

      /* Get the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
    }

/* v.2.0.3 --> */
  /* Get the current fragment's first atom */
  bond_ptr = pdbentry_ptr->first_bond_ptr;

  /* Loop over the fragment's bonds to write out */
  while (bond_ptr != NULL)
    {
      /* Check that both this bond's atoms are from the ligand residue */
      if (bond_ptr->type == COVALENT &&
          bond_ptr->atom1_ptr->residue_ptr->interaction == IS_LIGAND &&
          bond_ptr->atom2_ptr->residue_ptr->interaction == IS_LIGAND)
        {
          /* Check that both atom numbers fall within the range of our
             ligand */
          if (bond_ptr->atom1_ptr->atom_number >= min_atno &&
              bond_ptr->atom1_ptr->atom_number <= max_atno &&
              bond_ptr->atom2_ptr->atom_number >= min_atno &&
              bond_ptr->atom2_ptr->atom_number <= max_atno)
            {
              /* Write out this bond as a CONECT record */
              fprintf(file_ptr,"CONECT%5d%5d\n",
                      bond_ptr->atom1_ptr->atom_number,
                      bond_ptr->atom2_ptr->atom_number);
            }
        }

      /* Get the next bond */
      bond_ptr = bond_ptr->next_bond_ptr;
    }
/* <-- v.2.0.3 */

  /* Close the output file */
  fclose(file_ptr);
}
/***********************************************************************

write_flat_coords  -  Write out the transformed coords to the
                      ligflat.pdb file

***********************************************************************/

void write_flat_coords(struct pdbentry *pdbentry_ptr,
                       struct parameters params)
{
  int iatom, natoms, nout;

  struct atom *atom_ptr;
  struct residue *residue_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  nout = 0;

  /* Open the output ligflat.pdb file */
  file_ptr = open_file("ligflat.pdb",params);

  /* Write out the header record */
  fprintf(file_ptr,"REMARK PDB CODE: %s\n",pdbentry_ptr->pdb_code);

  /* Get the first of the PDB entry's residues */
  residue_ptr = pdbentry_ptr->first_residue_ptr;

  /* Loop over the residues to write out the flattened coords */
  while (residue_ptr != NULL)
    {
      /* Get the first atom */
      atom_ptr = residue_ptr->first_atom_ptr;
      natoms = residue_ptr->natoms;
      iatom = 0;

      /* Loop over the residue's atoms to write out transformed coords */
      while (atom_ptr != NULL && iatom < natoms)
        {
          /* Write out */
          fprintf(file_ptr,"ATOM  %5d %s %s %c%s   "
                  "%8.3f%8.3f%8.3f%6.2f%6.2f          %s\n",
                  atom_ptr->atom_number,atom_ptr->atom_name,
                  residue_ptr->res_name,
                  residue_ptr->chain,
                  residue_ptr->res_num,
                  atom_ptr->coords[X],
                  atom_ptr->coords[Y],
                  atom_ptr->coords[Z],
                  atom_ptr->occupancy,
                  atom_ptr->bvalue,atom_ptr->element);

          /* Increment count of atoms written out */
          nout++;

          /* Get the next atom record and increment atom count */
          atom_ptr = atom_ptr->next_atom_ptr;
          iatom++;
        }

      /* Get the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
    }

  /* Close the output file */
  fclose(file_ptr);
}
/***********************************************************************

score_interactions  -  Score the interactions made by the end-atoms of
                       the current symmetrical side chain

***********************************************************************/

int score_interactions(struct atom *atom_swap_ptr[2][2],int flag)
{
  int iatom, best_score, score, total_score;

  struct atom *atom1_ptr, *atom2_ptr, *other_atom1_ptr, *other_atom2_ptr;
  struct bond *bond1_ptr, *bond2_ptr;
  struct residue *residue2_ptr;

  /* Initialise variables */
  total_score = 0;

  /* Loop over the two atoms in the first sidechain */
  for (iatom = 0; iatom < 2; iatom++)
    {
      /* Get this atom */
      atom1_ptr = atom_swap_ptr[0][iatom];

      /* Get the atom we are comparing against */
      if (flag == UNSWAPPED)
        atom2_ptr = atom_swap_ptr[1][iatom];
      else
        atom2_ptr = atom_swap_ptr[1][1 - iatom];

      /* Get pointer to the first bond record */
      bond1_ptr = atom1_ptr->residue_ptr->pdbentry_ptr->first_bond_ptr;

      /* Loop over the bonds to find those involving the current atom */
      while (bond1_ptr != NULL)
        {
          /* Initialise interacting atom */
          other_atom1_ptr = NULL;

          /* If bond involves current atom, save it */
          if (bond1_ptr->atom1_ptr == atom1_ptr)
            other_atom1_ptr = bond1_ptr->atom2_ptr;
          else if (bond1_ptr->atom2_ptr == atom1_ptr)
            other_atom1_ptr = bond1_ptr->atom1_ptr;

          /* If have an interacting atom, store it */
          if (other_atom1_ptr != NULL && (bond1_ptr->type == COVALENT ||
                                          bond1_ptr->type == HBOND ||
                                          bond1_ptr->type == CONTACT) &&
/* v.2.0.2 --> */
//              residue1_ptr != other_atom1_ptr->residue_ptr &&
              atom1_ptr->residue_ptr != other_atom1_ptr->residue_ptr &&
/* <-- v.2.0.2 */
              other_atom1_ptr->residue_ptr->interaction == IS_LIGAND)
            {
              /* Initialise best score */
              best_score = 0;

              /* Get pointer to the first bond record from the second
                 structure */
              bond2_ptr
                = atom2_ptr->residue_ptr->pdbentry_ptr->first_bond_ptr;

              /* Loop over the bonds to find any that are equivalent to the
                 current bond */
              while (bond2_ptr != NULL)
                {
                  /* Initialise interacting atom */
                  other_atom2_ptr = NULL;

                  /* If bond involves current atom, save it */
                  if (bond2_ptr->atom1_ptr == atom2_ptr)
                    other_atom2_ptr = bond2_ptr->atom2_ptr;
                  else if (bond2_ptr->atom2_ptr == atom2_ptr)
                    other_atom2_ptr = bond2_ptr->atom1_ptr;

                  /* If have an interacting atom, store it */
                  if (other_atom2_ptr != NULL &&
                      (bond2_ptr->type == COVALENT ||
                       bond2_ptr->type == HBOND ||
                       bond2_ptr->type == CONTACT) &&
                      residue2_ptr != other_atom2_ptr->residue_ptr &&
                      other_atom2_ptr->residue_ptr->interaction == IS_LIGAND)
                    {
                      /* Add contribution to score */
                      if (bond1_ptr->type == HBOND &&
                          bond2_ptr->type == HBOND)
                        score = 5;
                      else if (bond1_ptr->type == COVALENT &&
                          bond2_ptr->type == COVALENT)
                        score = 5;
                      else if (bond1_ptr->type == CONTACT &&
                          bond2_ptr->type == CONTACT)
                        score = 2;
                      else
                        score = 1;

                      /* If this is the best score so far, save it */
                      if (score > best_score)
                        best_score = score;
                    }

                  /* Get pointer to next bond record */
                  bond2_ptr = bond2_ptr->next_bond_ptr;
                }

              /* Add the best score to the total */
              total_score = total_score + best_score;
            }

          /* Get pointer to next bond record */
          bond1_ptr = bond1_ptr->next_bond_ptr;
        }
    }

  /* Return the score */
  return(total_score);
}
/***********************************************************************

check_swappable_atoms  -  Check whether the residues have symmetrical
                          side chains and, if so, whether the equiv atoms
                          should be swapped to give a better match to
                          the interactions

***********************************************************************/

int check_swappable_atoms(struct residue *residue1_ptr,
                          struct residue *residue2_ptr,
                          struct atom *atom_swap_ptr[2][2])
{
#define NRESID    3

  static const char SYMM_RESID[NRESID][3][5] = 
    {
      {"ARG", " NH1", " NH2"},
      {"ASP", " OD1", " OD2"},
      {"GLU", " OE1", " OE2"}
    };

  int i, iatom, iresid, j, loop, natoms, score[2], sresid;
  int swap;

  struct atom *atom_ptr;

  /* Initialise variables */
  for (i = 0; i < 2; i++)
    for (j = 0; j < 2; j++)
      atom_swap_ptr[i][j] = NULL;
  swap = FALSE;
  sresid = -1;

  /* If residues aren't the same, then return */
  if (strcmp(residue1_ptr->res_name,residue2_ptr->res_name))
    return(swap);

  /* Loop over the symmetrical side chains */
  for (iresid = 0; iresid < NRESID && sresid == -1; iresid++)
    {
      /* If this is the right residue, check atom interactions */
      if (!strcmp(residue1_ptr->res_name,SYMM_RESID[iresid][0]))
        {
          /* Initialise swapped and unswapped scores */
          score[0] = score[1] = 0;

          /* Loop twice, once for each residue */
          for (loop = 0; loop < 2; loop++)
            {
              /* If first loop, get the first atom of the first residue */
              if (loop == 0)
                {
                  atom_ptr = residue1_ptr->first_atom_ptr;
                  natoms = residue1_ptr->natoms;
                }
              else
                {
                  atom_ptr = residue2_ptr->first_atom_ptr;
                  natoms = residue2_ptr->natoms;
                }
              iatom = 0;

              /* Loop over the residue's atoms to find the symmetrical ones */
              while (atom_ptr != NULL && iatom < natoms)
                {
                  /* If this is one of the symmetrical atoms, save its
                     pointer */
                  if (!strcmp(atom_ptr->atom_name,SYMM_RESID[iresid][1]))
                    atom_swap_ptr[loop][0] = atom_ptr;
                  else if (!strcmp(atom_ptr->atom_name,SYMM_RESID[iresid][2]))
                    atom_swap_ptr[loop][1] = atom_ptr;

                  /* Get the next atom record and increment atom count */
                  atom_ptr = atom_ptr->next_atom_ptr;
                  iatom++;
                }
            }

          /* Save the residue position */
          sresid = iresid;
        }
    }

  /* Count the atoms for which we have valid pointers */
  natoms = 0;
  for (i = 0; i < 2; i++)
    for (j = 0; j < 2; j++)
      if (atom_swap_ptr[i][j] != NULL)
        natoms++;

  /* If have all the required atoms, then can check the interactions */
  if (natoms == 4)
    {
      /* Score the bond similarity in the unswapped arrangement */
      score[UNSWAPPED] = score_interactions(atom_swap_ptr,UNSWAPPED);
      score[SWAPPED] = score_interactions(atom_swap_ptr,SWAPPED);

      /* If the swapped score is better than the unswapped score, then
         set the flag to equiovalence the swapped orientation */
      if (score[SWAPPED] > score[UNSWAPPED])
        swap = TRUE;
    }

  /* Return whether swap to be made */
  return(swap);
}
/* v.1.3 --> */
/***********************************************************************

write_ligand_restraints  -  Write out the flattened ligand coords
                            as restraints for the full LIGPLOT run

***********************************************************************/

int write_ligand_restraints(FILE *file_ptr,struct residue *residue_ptr)
{
  int iatom, natoms, nout;

  struct atom *atom_ptr;

  /* Initialise variables */
  natoms = nout = 0;

  /* Get the first atom */
  atom_ptr = residue_ptr->first_atom_ptr;
  natoms = residue_ptr->natoms;
  iatom = 0;

  /* Loop over the residue's atoms to write out the flattened LIGPLOT
     coords */
  while (atom_ptr != NULL && iatom < natoms)
    {
      /* Write out */
      fprintf(file_ptr,"ATOM  %5d %s %s %c%s   %8.3f%8.3f%8.3f"
              "  1.00 %5.2f   1.000\n",
              atom_ptr->atom_number,
              atom_ptr->atom_name,
              atom_ptr->residue_ptr->res_name,
              atom_ptr->residue_ptr->chain,
              atom_ptr->residue_ptr->res_num,
              atom_ptr->coords[X],
              atom_ptr->coords[Y],
              atom_ptr->coords[Z],atom_ptr->anchor_weight);

      /* Increment count of atoms written out */
      nout++;

      /* Get the next atom record and increment atom count */
      atom_ptr = atom_ptr->next_atom_ptr;
      iatom++;
    }

  /* Return number of residue restraints written out */
  return(nout);
}
/* <-- v.1.3 */
/***********************************************************************

write_equiv_residues  -  Write out the residue's restraints to the
                         ligplus.rcm file based on equivalenced residue

***********************************************************************/

int write_equiv_residues(FILE *file_ptr,struct residue *residue1_ptr,
                         struct residue *residue2_ptr,int loop)
{
  char atom_name[5];

  int iatom, jatom, natoms1, natoms2, nout;
  int found, swap;

  struct atom *atom1_ptr, *atom2_ptr, *atom_swap_ptr[2][2];

  /* Initialise variables */
  nout = 0;

  /* Check whether the two residues have symmetrical side chains and, if
     so, a swap of the equivalnt end atoms gives a better match to the
     H-bonds */
  swap
    = check_swappable_atoms(residue1_ptr,residue2_ptr,atom_swap_ptr);

  /* Get the first atom */
  atom1_ptr = residue1_ptr->first_atom_ptr;
  natoms1 = residue1_ptr->natoms;
  iatom = 0;

  /* Loop over the residue's atoms to get the LIGPLOT coords from
     the corresponding atoms in the mapped residue */
  while (atom1_ptr != NULL && iatom < natoms1)
    {
      /* Get the atom name in the other residue that we are looking for */
      if (atom1_ptr == atom_swap_ptr[0][0] && swap == TRUE)
        strcpy(atom_name,atom_swap_ptr[1][1]->atom_name);
      else if (atom1_ptr == atom_swap_ptr[0][1] && swap == TRUE)
        strcpy(atom_name,atom_swap_ptr[1][0]->atom_name);
      else
        strcpy(atom_name,atom1_ptr->atom_name);

      /* Get the first atom of the LIGPLOT residue */
      atom2_ptr = residue2_ptr->first_atom_ptr;
      natoms2 = residue2_ptr->natoms;
      jatom = 0;
      found = FALSE;
 
      /* Loop over the atoms to find the equivalent to the one
         of interest */
      while (atom2_ptr != NULL && jatom < natoms2 && found == FALSE)
        {
          /* If the equivalent atom, write out the LIGPLOT coords */
          if (!strcmp(atom2_ptr->atom_name,atom_name))
            {
              /* Write out */
              if (loop == 1 || residue1_ptr->interaction == IS_LIGAND)
                fprintf(file_ptr,"ATOM  %5d %s %s %c%s   "
                        "%8.3f%8.3f%8.3f  1.00  1.00   1.000\n",
                        atom1_ptr->atom_number,atom1_ptr->atom_name,
                        residue1_ptr->res_name,
                        residue1_ptr->chain,
                        residue1_ptr->res_num,
                        atom2_ptr->coords[X],
                        atom2_ptr->coords[Y],
                        atom2_ptr->coords[Z]);
                      
              /* Set flag and increment count of atoms written out */
              found = TRUE;
              nout++;
            }

          /* Get the next atom record and increment atom count */
          atom2_ptr = atom2_ptr->next_atom_ptr;
          jatom++;
        }

      /* Get the next atom record and increment atom count */
      atom1_ptr = atom1_ptr->next_atom_ptr;
      iatom++;
    }

  /* Return number of residue restraints written out */
  return(nout);
}
/***********************************************************************

have_residue  -  See if we have the residue in the residue list

***********************************************************************/

int have_residue(struct resseq *first_resseq_ptr,
                 struct residue *residue_ptr)
{
  int found;

  struct resseq *resseq_ptr;

  /* Initialise variables */
  found = FALSE;

  /* Get pointer to the first resseq */
  resseq_ptr = first_resseq_ptr;

  /* Loop over the resseqs to write out */
  while (resseq_ptr != NULL && found == FALSE)
    {
      /* If this is the right residue, then set flag */
      if (resseq_ptr->residue_ptr == residue_ptr)
        found = TRUE;

      /* Get pointer to next resseq record */
      resseq_ptr = resseq_ptr->next_resseq_ptr;
    }

  /* Return whether residue found */
  return(found);
}
/***********************************************************************

free_resseq  -  Free up the memory taken up by the list of entries

***********************************************************************/

void free_resseq(struct resseq *first_resseq_ptr)
{
  struct resseq *resseq_ptr, *next_resseq_ptr;

  /* Get first resseq record */
  resseq_ptr = first_resseq_ptr;

  /* Loop through the resseq records until done */
  while (resseq_ptr != NULL)
    {
      /* Get pointer to the next resseq record */
      next_resseq_ptr = resseq_ptr->next_resseq_ptr;

      /* Free up memory taken by current resseq */
      free(resseq_ptr);

      /* Go to the next resseq */
      resseq_ptr = next_resseq_ptr;
    }
}
/***********************************************************************

free_bonds  -  Free up the memory taken up by the list of entries

***********************************************************************/

void free_bonds(struct bond *first_bond_ptr)
{
  struct bond *bond_ptr, *next_bond_ptr;

  /* Get first bond record */
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

find_bond  -  Find the given bond record

***********************************************************************/

struct bond *find_bond(struct bond *first_bond_ptr,struct atom *atom1_ptr,
                       struct atom *atom2_ptr)
{
  struct bond *bond_ptr, *match_bond_ptr;

  /* Initialise variables */
  match_bond_ptr = NULL;

  /* Get pointer to the first bond */
  bond_ptr = first_bond_ptr;

  /* Loop over the bonds to write out */
  while (bond_ptr != NULL && match_bond_ptr == NULL)
    {
      /* Check if this is the required record */
      if ((bond_ptr->atom1_ptr == atom1_ptr &&
           bond_ptr->atom2_ptr == atom2_ptr) ||
          (bond_ptr->atom1_ptr == atom2_ptr &&
           bond_ptr->atom2_ptr == atom1_ptr))
        {
          /* Save the bond pointer */
          match_bond_ptr = bond_ptr;
        }

      /* Get pointer to next bond record */
      bond_ptr = bond_ptr->next_bond_ptr;
    }
      
  /* Return match */
  return(match_bond_ptr);
}
/***********************************************************************

write_equiv_ligatoms  -  Write out the restraints for the ligand's atoms
                         by idententifying the "best" equivalent atom in
                         the reference ligand

***********************************************************************/

int write_equiv_ligatoms(FILE *file_ptr,struct residue *residue1_ptr,
                         struct residue *residue2_ptr)
{
  int iatom, jatom, natoms, natoms1, natoms2, nbonds, nout, nmatch, score;
  int best_match;

  struct atom *atom1_ptr, *atom2_ptr, *best_atom_ptr, *map_atom_ptr;
  struct atom *other_atom_ptr;
  struct atomap *atomap_ptr, *first_atomap_ptr, *last_atomap_ptr;
  struct bond *bond_ptr, *have_bond_ptr;
  struct bond *save_bond_ptr, *first_bond_ptr, *last_bond_ptr;
  struct residue *map_residue_ptr, *other_residue_ptr;

  /* Initialise variables */
  natoms = nout = 0;
  first_atomap_ptr = last_atomap_ptr = NULL;

  /* Get the first atom */
  atom1_ptr = residue1_ptr->first_atom_ptr;
  natoms1 = residue1_ptr->natoms;
  iatom = 0;

  /* Loop over the residue's atoms to get the LIGPLOT coords from
     the corresponding atoms in the mapped residue */
  while (atom1_ptr != NULL && iatom < natoms1)
    {
      /* Initialise list of interactions involving current atom */
      save_bond_ptr = first_bond_ptr = last_bond_ptr = NULL;
      nbonds = 0;

      /* Get pointer to the first bond record */
      bond_ptr = residue1_ptr->pdbentry_ptr->first_bond_ptr;

      /* Loop over the bonds to find those involving the current atom */
      while (bond_ptr != NULL)
        {
          /* Initialise interacting atom */
          other_atom_ptr = NULL;

          /* If bond involves current atom, save it */
          if (bond_ptr->atom1_ptr == atom1_ptr)
            other_atom_ptr = bond_ptr->atom2_ptr;
          else if (bond_ptr->atom2_ptr == atom1_ptr)
            other_atom_ptr = bond_ptr->atom1_ptr;

          /* If have an interacting atom, store it */
          if (other_atom_ptr != NULL && (bond_ptr->type == COVALENT ||
                                         bond_ptr->type == HBOND ||
                                         bond_ptr->type == CONTACT) &&
              residue1_ptr != other_atom_ptr->residue_ptr)
            {
              /* Check whether we already have this bond */
              have_bond_ptr
                = find_bond(first_bond_ptr,atom1_ptr,other_atom_ptr);

              /* Save the bond */
              if (have_bond_ptr == NULL)
                {
                  save_bond_ptr
                    = create_bond_record(&first_bond_ptr,
                                         &last_bond_ptr,bond_ptr->type,
                                         bond_ptr->length,
                                         bond_ptr->bond_order,
                                         atom1_ptr,other_atom_ptr);

                  /* Increment bond count */
                  nbonds++;
                }
            }

          /* Get pointer to next bond record */
          bond_ptr = bond_ptr->next_bond_ptr;
        }

      /* Get the first atom of the LIGPLOT residue */
      atom2_ptr = residue2_ptr->first_atom_ptr;
      natoms2 = residue2_ptr->natoms;
      jatom = 0;
      best_atom_ptr = NULL;
      best_match = 0;
      if (nbonds == 0)
        atom2_ptr = NULL;
 
      /* Loop over the atoms to find the equivalent to the one
         of interest */
      while (atom2_ptr != NULL && jatom < natoms2)
        {
          /* Initialise number of matching bonds */
          nmatch = 0;

          /* Get pointer to the first bond record */
          bond_ptr = residue2_ptr->pdbentry_ptr->first_bond_ptr;

          /* Loop over the bonds to find those involving the current atom
             to see how many are equivalent to the saved bonds */
          while (bond_ptr != NULL)
            {
              /* Initialise interacting atom */
              other_atom_ptr = NULL;

              /* If bond involves current atom, save it */
              if (bond_ptr->atom1_ptr == atom2_ptr)
                other_atom_ptr = bond_ptr->atom2_ptr;
              else if (bond_ptr->atom2_ptr == atom2_ptr)
                other_atom_ptr = bond_ptr->atom1_ptr;

              /* If have a bond involving current atom, check if the
                 bond is equivalent to any from the saved set */
              if (other_atom_ptr != NULL && (bond_ptr->type == COVALENT ||
                                             bond_ptr->type == HBOND ||
                                             bond_ptr->type == CONTACT) &&
                  residue2_ptr != other_atom_ptr->residue_ptr)
                {
                  /* Get the residue belonging to the interacting atom */
                  other_residue_ptr = other_atom_ptr->residue_ptr;

                  /* Get the residue this maps to */
                  map_residue_ptr = other_residue_ptr->map_residue_ptr;

                  /* Get the first of the saved bonds */
                  save_bond_ptr = first_bond_ptr;
                  if (map_residue_ptr == NULL)
                    save_bond_ptr = NULL;

                  /* Loop over the bonds to any involving the mapped
                     residue */
                  while (save_bond_ptr != NULL)
                    {
                      /* Initialise mapped atom */
                      map_atom_ptr = NULL;

                      /* If bond involves mapped residue, then check
                         if it matches the current bond */
                      if (save_bond_ptr->atom1_ptr->residue_ptr
                          == map_residue_ptr)
                        map_atom_ptr = save_bond_ptr->atom1_ptr;
                      else if (save_bond_ptr->atom2_ptr->residue_ptr
                          == map_residue_ptr)
                        map_atom_ptr = save_bond_ptr->atom2_ptr;

                      /* If have a mapped atom, increment match count */
                      if (map_atom_ptr != NULL)
                        {
                          /* Initialize match score */
                          score = 0;

                          /* If bond types are the same, increment
                             match score */
                          if (bond_ptr->type == save_bond_ptr->type)
                            score = score + 3;

                          /* Add extra for H-bond rather than non-bonded
                             contact */
                          if (bond_ptr->type == save_bond_ptr->type &&
                              bond_ptr->type == HBOND)
                            {
                              score = score + 10;

                              /* If both interacting atoms are the same
                                 element, increment match count */
                              if (!strcmp(other_atom_ptr->element,
                                          map_atom_ptr->element))
                                score++;
                            }

                          /* If atom and residue names match, increment
                             more */
                          if (score > 0)
                            {
                              if (!strcmp(other_atom_ptr->atom_name,
                                          map_atom_ptr->atom_name) &&
                                  !strcmp(other_atom_ptr->residue_ptr->res_name,
                                          map_atom_ptr->residue_ptr->res_name))
                                score = score + 2;
                            }

                          /* Add score to total match score */
                          nmatch = nmatch + score;
                        }

                      /* Get pointer to next saved bond record */
                      save_bond_ptr = save_bond_ptr->next_bond_ptr;
                    }
                }
              
              /* Get pointer to next bond record */
              bond_ptr = bond_ptr->next_bond_ptr;
            }

          /* If have any matches and both atoms are the same element,
             increment match count */
          if (nmatch > 0 && !strcmp(atom1_ptr->element,atom2_ptr->element))
            nmatch++;

          /* If best match so far, save the atom */
          if (nmatch > best_match)
            {
              /* Save atom */
              best_atom_ptr = atom2_ptr;

              /* Save best score */
              best_match = nmatch;
            }

          /* Get the next atom record and increment atom count */
          atom2_ptr = atom2_ptr->next_atom_ptr;
          jatom++;
        }

      /* If have a reasonable match, save the atom mapping */
      if (best_atom_ptr != NULL && best_match > MIN_BEST_MATCH)
        {
          /* Check whether we already have a mapping to the mapped atom */
          atomap_ptr = find_atomap(first_atomap_ptr,best_atom_ptr);

          /* If not found, create a new map record */
          if (atomap_ptr == NULL)
            {
              /* Create the atom mapping */
              atomap_ptr
                = create_atomap_record(&first_atomap_ptr,&last_atomap_ptr,
                                       atom1_ptr,best_atom_ptr,best_match);

              /* Increment count of saved atoms */
              natoms++;
            }

          /* If found, then update if current mapping has a higher score */
          else
            {
              /* If current mapping has a higher score, replace old
                 mapping */
              if (best_match > atomap_ptr->score)
                {
                  atomap_ptr->atom_ptr = atom1_ptr;
                  atomap_ptr->score = best_match;
                }
            }

          /* Increment count of atoms written out */
          nout++;
        }

      /* Free up memory taken by the bond records created here */
      free_bonds(first_bond_ptr);

      /* Get the next atom record and increment atom count */
      atom1_ptr = atom1_ptr->next_atom_ptr;
      iatom++;
    }

  /* Get the first of the atom mappings */
  atomap_ptr = first_atomap_ptr;

  /* Loop through to write out */
  while (atomap_ptr != NULL)
    {
      /* Write out */
      fprintf(file_ptr,"ATOM  %5d %s %s %c%s   %8.3f%8.3f%8.3f"
              "  1.00  1.00   1.000\n",
              atomap_ptr->atom_ptr->atom_number,
              atomap_ptr->atom_ptr->atom_name,
              atomap_ptr->atom_ptr->residue_ptr->res_name,
              atomap_ptr->atom_ptr->residue_ptr->chain,
              atomap_ptr->atom_ptr->residue_ptr->res_num,
              atomap_ptr->map_atom_ptr->coords[X],
              atomap_ptr->map_atom_ptr->coords[Y],
              atomap_ptr->map_atom_ptr->coords[Z]);

      /* Get pointer to next record */
      atomap_ptr = atomap_ptr->next_atomap_ptr;
    }

  /* Return number of residue restraints written out */
  return(nout);
}
/***********************************************************************

write_residue_cofm  -  Write out the residue's C of M restraint

***********************************************************************/

int write_residue_cofm(FILE *file_ptr,struct residue *residue1_ptr,
                       struct residue *residue2_ptr,int loop)
{
  int jatom, natoms2, nout;

  double cofg[2];

  struct atom *atom2_ptr;

  /* Initialise variables */
  nout = 0;

  /* Initialise centre of mass */
  cofg[0] = cofg[1] = 0.0;

/* v.2.0 --> */
  /* If this is a dummy residue, then use CofM from residue */
/* v.2.1 --> */
//  if (!strcmp(residue2_ptr->res_name,"DUM"))
  if (!strcmp(residue2_ptr->res_name,"DUM") || residue2_ptr->natoms == 0)
/* <-- v.2.1 */
    {
      /* Get the average x- and y-values */
      cofg[0] = residue2_ptr->valmin[MAPPED][X]
        + residue2_ptr->valmax[MAPPED][X];
      cofg[1] = residue2_ptr->valmin[MAPPED][Y]
        + residue2_ptr->valmax[MAPPED][Y];
      natoms2 = 2;
    }
  else
    {
/* <-- v.2.0 */
      /* Get the first of this residue's atoms */
      atom2_ptr = residue2_ptr->first_atom_ptr;
      natoms2 = residue2_ptr->natoms;
      jatom = 0;

      /* Loop over the residue's atoms to calculate centre of mass */
      while (atom2_ptr != NULL && jatom < natoms2)
        {
          /* Add this atom's 2D coords */
          cofg[0] = cofg[0] + atom2_ptr->coords[X];
          cofg[1] = cofg[1] + atom2_ptr->coords[Y];

          /* Get the next atom record and increment atom count */
          atom2_ptr = atom2_ptr->next_atom_ptr;
          jatom++;
        }
/* v.2.0 --> */
    }
/* <-- v.2.0 */
              
  /* Calculate centre of mass and write out */
  if (natoms2 > 0)
    {
      cofg[0] = cofg[0] / natoms2;
      cofg[1] = cofg[1] / natoms2;
              
      /* Write out */
      if (loop == 1 || residue1_ptr->interaction == IS_LIGAND)
        fprintf(file_ptr,"RESDUE      CofM %s %c%s   %8.3f%8.3f\n",
                residue1_ptr->res_name,
                residue1_ptr->chain,
                residue1_ptr->res_num,cofg[0],cofg[1]);
      nout++;
    }

  /* Return number of residue restraints written out */
  return(nout);
}
/***********************************************************************

write_residue  -  Write out the coords of the given residue

***********************************************************************/

int write_residue(FILE *file_ptr, struct residue *residue_ptr,
                  int atom_number)
{
  int iatom, natoms;

  struct atom *atom_ptr;

  /* Get the first atom */
  atom_ptr = residue_ptr->first_atom_ptr;
  natoms = residue_ptr->natoms;
  iatom = 0;

  /* Loop over the residue's atoms to write out transformed coords */
  while (atom_ptr != NULL && iatom < natoms)
    {
      /* Increment atom number */
      atom_number++;

      /* Write out */
      fprintf(file_ptr,"ATOM  %5d %s %s %c%s   "
              "%8.3f%8.3f%8.3f%6.2f%6.2f          %s\n",
              atom_number,
              atom_ptr->atom_name,
              residue_ptr->res_name,
              residue_ptr->chain,
              residue_ptr->res_num,
              atom_ptr->coords[X],
              atom_ptr->coords[Y],
              atom_ptr->coords[Z],
              atom_ptr->occupancy,
              atom_ptr->bvalue,atom_ptr->element);

      /* Get the next atom record and increment atom count */
      atom_ptr = atom_ptr->next_atom_ptr;
      iatom++;
    }

  /* Return the current atom number */
  return(atom_number);
}
/***********************************************************************

write_ligfit  -  Write out the equivalent ligand atoms obtained from
                 hbadd's substructure search

***********************************************************************/

void write_ligfit(FILE *file_ptr,struct atomap *first_atomap_ptr)
{
  float anchor_weight;

  struct atomap *atomap_ptr;

  /* Get the first of the atom mappings */
  atomap_ptr = first_atomap_ptr;

  /* Loop through to write out */
  while (atomap_ptr != NULL)
    {
      /* Determine what the weight of this anchor position is */
      anchor_weight = atomap_ptr->anchor_weight;

      /* Write out */
      fprintf(file_ptr,"ATOM  %5d %s %s %c%s   %8.3f%8.3f%8.3f"
              "  1.00 %5.2f   1.000\n",
              atomap_ptr->atom_ptr->atom_number,
              atomap_ptr->atom_ptr->atom_name,
              atomap_ptr->atom_ptr->residue_ptr->res_name,
              atomap_ptr->atom_ptr->residue_ptr->chain,
              atomap_ptr->atom_ptr->residue_ptr->res_num,
              atomap_ptr->map_atom_ptr->coords[X],
              atomap_ptr->map_atom_ptr->coords[Y],
              atomap_ptr->map_atom_ptr->coords[Z],anchor_weight);

      /* Get pointer to next record */
      atomap_ptr = atomap_ptr->next_atomap_ptr;
    }
}
/***********************************************************************

write_blocking_ligands  -  Write out previous ligand(s) as dummy records
                           to block this region of the plot

***********************************************************************/

void write_blocking_ligands(FILE *file_ptr,int atom_number,
                            struct pdbentry *first_pdbentry_ptr)
{
  char chain, res_name[4], res_num[6];

  struct pdbentry *pdbentry_ptr;
  struct residue *residue_ptr;

  /* Get the first PDB entry */
  pdbentry_ptr = first_pdbentry_ptr;

  /* Loop through the PDB entries */
  while (pdbentry_ptr != NULL)
    {
      /* Get this entry's first residue */
      residue_ptr = pdbentry_ptr->first_residue_ptr;

      /* Loop over the residues to write out ligand coords */
      while (residue_ptr != NULL)
        {
          /* If ligand, then write out */
          if (residue_ptr->interaction == IS_LIGAND)
            {
              /* Save the residue details */
              strcpy(res_name,residue_ptr->res_name);
              strcpy(res_num,residue_ptr->res_num);
              chain = residue_ptr->chain;

              /* Replace by fake details */
              strcpy(residue_ptr->res_name,"...");
              strcpy(residue_ptr->res_num,"   1 ");
              residue_ptr->chain = '.';

              /* Write out using fake residue name */
              atom_number
                = write_residue(file_ptr,residue_ptr,atom_number);

              /* Restore correct residue details */
              strcpy(residue_ptr->res_name,res_name);
              strcpy(residue_ptr->res_num,res_num);
              residue_ptr->chain = chain;
            }

          /* Get pointer to next residue record */
          residue_ptr = residue_ptr->next_residue_ptr;
        }

      /* Get the next entry */
      pdbentry_ptr = pdbentry_ptr->next_pdbentry_ptr;
    }
}
/***********************************************************************

get_equivalences  -  Use the fitting to determine the atom and residue
                     equivalences for use in restraining the second plot

***********************************************************************/

/* v.1.3.7 --> */
int get_equivalences(int loop,struct pdbentry *ref_pdbentry_ptr,
                     struct pdbentry *pdbentry_ptr,
                     struct residue **fst_ligand_residue_ptr,
                     int fit_type,struct atomap *first_atomap_ptr,
                     struct parameters params)
// int get_equivalences(int loop,struct pdbentry *pdbentry_ptr,
//                     struct residue **fst_ligand_residue_ptr,
//                     int fit_type,struct atomap *first_atomap_ptr,
//                     struct parameters params)
/* <-- v.1.3.7 */
{
  int atom_number, nmap;
/* v.1.3 --> */
//  int found, have_ligfit, overlay, wanted;
  int done, have_ligfit, overlay;
  int all_atom1, all_atom2;
/* <-- v.1.3 */

  struct pdbentry *pdbentry2_ptr;
  struct residue *residue1_ptr, *residue2_ptr;
  struct residue *first_ligand_residue_ptr;

  FILE *file_ptr, *lgffile_ptr;

  /* Initialise variables */
  atom_number = 0;
  *fst_ligand_residue_ptr = first_ligand_residue_ptr = NULL;
  overlay = FALSE;
  if (first_atomap_ptr != NULL)
    have_ligfit = TRUE;
  else
    have_ligfit = FALSE;

  /* Open the appropriate output files */
/* v.2.0 --> */
//if (loop == 0)
  if (loop == 0 || params.dimplot == TRUE)
/* <-- v.2.0 */
    {
/* v.2.0 --> */
      if (params.dimplot == TRUE)
        file_ptr = open_file("dimplot.rcm",params);
      else
/* <-- v.2.0 */
        file_ptr = open_file("ligand.rcm",params);
      lgffile_ptr = open_file("ligplus.lgf",params);
    }
  else
    {
      file_ptr = open_file("ligplus.rcm",params);
    }

  /* Write out the header records */
  fprintf(file_ptr,"                 Res.  Res       Flattened    \n");
  fprintf(file_ptr,"            Atom Name  Num      coordinates   \n");
  fprintf(file_ptr,"            ---- --- -----    ----------------\n");

  /* Write out how the fit was made to the .lgf file */
/* v.2.0 --> */
//if (loop == 0)
  if (loop == 0 || params.dimplot == TRUE)
/* <-- v.2.0 */
/* v.1.3.7 --> */
    // fprintf(lgffile_ptr,"FIT_TYPE=%s\n",FIT_TYPE[fit_type]);
    {
      fprintf(lgffile_ptr,"FIT_TYPE=%s\n",FIT_TYPE[fit_type]);
      fprintf(lgffile_ptr,"ALIGNED=Chain %c of %s and chain %c of %s\n",
              ref_pdbentry_ptr->map_chain,ref_pdbentry_ptr->pdb_code,
              ref_pdbentry_ptr->map_chain,pdbentry_ptr->pdb_code);
    }
/* <-- v.1.3.7 */

  /* Get the first of the PDB entry's residues */
  residue1_ptr = pdbentry_ptr->first_residue_ptr;

  /* Loop over the residues to map equivalent atoms or residues */
  while (residue1_ptr != NULL)
    {
/* v.1.3 --> */
      /* Initialise flag */
      done = FALSE;

      /* If this is the ligand to be plotted and we have its flattened
         coords from LIGPLOT, then write these out */
      if (loop == 1 && residue1_ptr->interaction == IS_LIGAND)
        {
          /* Get the residue mapping */
          residue2_ptr = residue1_ptr->map_residue_ptr;

          /* If have the flattened coords, then write these out */
          if (residue2_ptr != NULL &&
              !strcmp(residue2_ptr->res_name,residue1_ptr->res_name) &&
              !strcmp(residue2_ptr->res_num,residue1_ptr->res_num) &&
              residue2_ptr->chain == residue1_ptr->chain &&
              residue2_ptr->natoms == residue1_ptr->natoms)
            {
              /* Write out the flattened coordinates to the .rcm
                 restraints file */
              write_ligand_restraints(file_ptr,residue2_ptr);

              /* Set flag that this residue done */
              done = TRUE;
            }
        }
/* <-- v.1.3 */

      /* If this residue has a mapping, then write out */
/* v.1.3 --> */
//      if (residue1_ptr->map_residue_ptr != NULL)
      if ((residue2_ptr = residue1_ptr->map_residue_ptr) != NULL &&
          residue2_ptr->in_plot == TRUE && done == FALSE)
/* <-- v.1.3 */
        {
          /* Get the residue it maps to */
/* v.1.3 --> */
//          residue2_ptr = residue1_ptr->map_residue_ptr;
/* <-- v.1.3 */

          /* Initialise count of atoms written out */
          nmap = 0;

/* v.1.3 --> */
          /* Determine whether both residues are all-atom */
          all_atom1 = all_atom2 = FALSE;
          if (residue1_ptr->interaction == HBOND_TO_LIGAND ||
              residue1_ptr->interaction == VIA_WATER_TO_LIGAND ||
              residue1_ptr->interaction == VIA_METAL_TO_LIGAND ||
/* v.2.0 --> */
//              residue1_ptr->interaction == COVALENT_TO_LIGAND)
              residue1_ptr->interaction == COVALENT_TO_LIGAND ||
              residue1_ptr->interaction == PROT_PROT_HBOND ||
              residue1_ptr->interaction == VIA_WATER ||
              residue1_ptr->interaction == VIA_METAL)
/* <-- v.2.0 */
            all_atom1 = TRUE;
          if (residue2_ptr->interaction == HBOND_TO_LIGAND ||
              residue2_ptr->interaction == VIA_WATER_TO_LIGAND ||
              residue2_ptr->interaction == VIA_METAL_TO_LIGAND ||
/* v.2.0 --> */
//              residue2_ptr->interaction == COVALENT_TO_LIGAND)
              residue2_ptr->interaction == COVALENT_TO_LIGAND ||
              residue2_ptr->interaction == PROT_PROT_HBOND ||
              residue2_ptr->interaction == VIA_WATER ||
              residue2_ptr->interaction == VIA_METAL)
/* <-- v.2.0 */
            all_atom2 = TRUE;
/* <-- v.1.3 */

          /* Determine whether both residues are H-bond groups, or
             both are identical ligands */
/* v.1.3 --> */
//          if ((residue1_ptr->interaction == HBOND_TO_LIGAND &&
//               residue2_ptr->interaction == HBOND_TO_LIGAND) ||
          if ((all_atom1 == TRUE && all_atom2 == TRUE) ||
/* <-- v.1.3 */
              (residue1_ptr->interaction == IS_LIGAND &&
               residue2_ptr->interaction == IS_LIGAND &&
/* v.1.3 --> */
//               have_ligfit == FALSE))
               !strcmp(residue1_ptr->res_name,residue2_ptr->res_name) &&
               have_ligfit == FALSE))
/* <-- v.1.3 */
            {
              /* Write out the all-atom positional restraints for this
                 residue */
              nmap
                = write_equiv_residues(file_ptr,residue1_ptr,
                                       residue2_ptr,loop);
            }

          /* If both residues are ligand residues, but are not identical,
             then match equivalent atoms */
          else if (residue1_ptr->interaction == IS_LIGAND &&
                   residue2_ptr->interaction == IS_LIGAND &&
                   have_ligfit == FALSE)
            {
              /* Write out the all-atom positional restraints for
                 equivalent ligand atoms */
              nmap
                = write_equiv_ligatoms(file_ptr,residue1_ptr,residue2_ptr);

              /* If this is the first ligand residue, store its pointer */
              if (first_ligand_residue_ptr == NULL)
                first_ligand_residue_ptr = residue1_ptr;

              /* If number of matched atoms less then number in current
                 residue, set flag to overlay ligands */
              if (nmap < residue1_ptr->natoms)
                overlay = TRUE;
            }

          /* If have fitting for the two ligands, don't want to write
             the centre of mass out */
          else if (residue1_ptr->interaction == IS_LIGAND &&
                   residue2_ptr->interaction == IS_LIGAND &&
                   have_ligfit == TRUE)
            {
              /* If this is the first ligand residue, store its pointer */
              if (first_ligand_residue_ptr == NULL)
                first_ligand_residue_ptr = residue1_ptr;

              /* Set flags */
              nmap = -1;
              overlay = TRUE;
            }

          /* If no atoms written out, then write out centre of mass
             instead */
          if (nmap == 0)
            write_residue_cofm(file_ptr,residue1_ptr,residue2_ptr,loop);

          /* If neither residue is a ligand residue, write out to the
             residue equivalences file, ligplus.lgf */
/* v.2.0 --> */
//          if (loop == 0 && residue1_ptr->interaction != IS_LIGAND &&
//              residue2_ptr->interaction != IS_LIGAND)
          if ((loop == 0 && residue1_ptr->interaction != IS_LIGAND &&
              residue2_ptr->interaction != IS_LIGAND) ||
              params.dimplot == TRUE)
/* <-- v.2.0 */
            {
              /* Get other PDB entry */
              pdbentry2_ptr = residue2_ptr->pdbentry_ptr;

              /* If have PDB and residue identifiers, write out
                 the mapping */
              if (pdbentry2_ptr != NULL && pdbentry2_ptr->pdb_id > -1 &&
                  residue2_ptr->residue_id != &null)
                {
                  /* Write the mapping */
                  fprintf(lgffile_ptr,"%d\t%s\t[%s %s %c]\n",
                          pdbentry2_ptr->pdb_id,residue2_ptr->residue_id,
                          residue1_ptr->res_name,residue1_ptr->res_num,
                          residue1_ptr->chain);
                }
            }
        }

      /* Get the next residue */
      residue1_ptr = residue1_ptr->next_residue_ptr;
    }

  /* If have ligand-ligand fit, write out the equivalent atoms */
/* v.1.3--> */
//  if (have_ligfit == TRUE)
  if (have_ligfit == TRUE && loop == 0)
    write_ligfit(file_ptr,first_atomap_ptr);
/* <--v.1.3 */

  /* If writing ligplus.rcm, write out previous ligand(s) as dummy records
     to block this region of the plot */
  if (loop == 1)
    write_blocking_ligands(file_ptr,atom_number,
                           pdbentry_ptr->next_pdbentry_ptr);

  /* Close the output files */
  fclose(file_ptr);
/* v.2.0 --> */
//if (loop == 0)
  if (loop == 0 || params.dimplot == TRUE)
/* <-- v.2.0 */
    fclose(lgffile_ptr);

  /* Return the first ligand residue */
  *fst_ligand_residue_ptr = first_ligand_residue_ptr;

  /* Return flag indicating whether we need to overlay the ligands
     because they are not identical */
  return(overlay);
}
/* v.2.0 --> */
/***********************************************************************

create_layout_record  -  Create and initialise a new layout record

***********************************************************************/

struct layout *create_layout_record(struct layout **fst_layout_ptr,
                                    struct layout **lst_layout_ptr,
                                    int pos,int chain_no,int x,int y,
                                    struct residue *residue_ptr)
{
  struct layout *layout_ptr, *first_layout_ptr, *last_layout_ptr;

  /* Initialise layout pointers */
  layout_ptr = NULL;
  first_layout_ptr = *fst_layout_ptr;
  last_layout_ptr = *lst_layout_ptr;

  /* Allocate memory for structure to hold layout info */
  layout_ptr = (struct layout *) malloc(sizeof(struct layout));
  if (layout_ptr == NULL)
    {
      printf("*** Can't allocate memory for struct layout\n");
      exit (1);
    }

  /* If this is the very first, save its pointer */
  if (first_layout_ptr == NULL)
    first_layout_ptr = layout_ptr;

  /* Add link from previous layout to the current one */
  if (last_layout_ptr != NULL)
    last_layout_ptr->next_layout_ptr = layout_ptr;
  last_layout_ptr = layout_ptr;

  /* Initialise the elements of the structure */
  layout_ptr->pos = pos;
  layout_ptr->chain_no = chain_no;
  layout_ptr->x = x;
  layout_ptr->y = y;
  layout_ptr->residue_ptr = residue_ptr;

  /* Initialise pointer to the next layout */
  layout_ptr->next_layout_ptr = NULL;

  /* Return current pointers */
  *fst_layout_ptr = first_layout_ptr;
  *lst_layout_ptr = last_layout_ptr;

  /* Return new layout pointer */
  return(layout_ptr);
}
/***********************************************************************

read_rcm_file  -  Read in the layout as computed by dimhtml

***********************************************************************/

int read_rcm_file(struct residue *first_residue_ptr,
                  struct layout **fst_layout_ptr,
                  struct layout **lst_layout_ptr,
                  struct parameters params)
{
  char file_name[FILENAME_LEN];
  char input_line[LINELEN + 1], number_string[20];
  char chain, res_name[4], res_num[6];

  int chain_no, len, nlayouts, x, y;

  struct layout *layout_ptr, *first_layout_ptr, *last_layout_ptr;
  struct residue *residue_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  file_name[0] = '\0';
  nlayouts = 0;
  res_name[0] = '\0';
  res_num[0] = '\0';
  x = y = 0;

  /* Initialise pointers */
  layout_ptr = first_layout_ptr = last_layout_ptr = NULL;

  /* Form name of the .rcm file */
  if (params.directory[0] != '\0')
    {
      strcpy(file_name,params.directory);
      strcat(file_name,DIR_SEP);
    }
  strcat(file_name,"dimhtml.rcm");

  /* Open the .rcm file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    {
      printf("*** ERROR. Unable to open dimhtml.rcm file: %s\n",
             file_name);
      exit(-1);
    }

  /* Loop while reading in records from the file */
  while (fgets(input_line,LINELEN,file_ptr) != NULL)
    {
      /* Truncate the string to strip off the newline */
      len = string_truncate(input_line,LINELEN);

      /* If this is a layout line, extract and store the details */
      if (!strncmp(input_line,"RESDUE",6))
        {
          /* Get layout name and number */
          strncpy(res_name,input_line+17,3);
          res_name[3] = '\0';
          strncpy(res_num,input_line+22,5);
          res_num[5] = '\0';
          chain = input_line[21];
          if (chain == '-')
            chain = ' ';

          /* Identify this residue */
          residue_ptr
            = find_residue(first_residue_ptr,res_name,res_num,chain);

          /* Get the x- and y-coordinates */
          strcpy(number_string,input_line+30);
          number_string[8] = '\0';
          x = (int) atof(number_string);
          strcpy(number_string,input_line+38);
          number_string[8] = '\0';
          y = (int) atof(number_string);

          /* Determine which of the two chains this is */
          if (x > 0)
            chain_no = 0;
          else
            chain_no = 1;

          /* Create layout record */
          layout_ptr
            = create_layout_record(&first_layout_ptr,&last_layout_ptr,
                                   nlayouts,chain_no,x,y,residue_ptr);

/* v.2.1--> */
          /* Assign the residue to the appropriate chain */
          if (!strcmp(residue_ptr->res_name,"HOH"))
            residue_ptr->chain_type = INTERFACE_WATER;
          else if (chain_no == 0)
            residue_ptr->chain_type = CHAIN1;
          else if (chain_no == 1)
            residue_ptr->chain_type = CHAIN2;
          else
            residue_ptr->chain_type = NON_INTERFACE;
/* <-- v.2.1 */

          /* Increment count of layouts stored */
          nlayouts++;
        }
    }

  /* Close the file */
  fclose(file_ptr);

  /* Show number of layout read in */
  printf("Number of layouts stored:                  %8d\n",nlayouts);

  /* Return pointer to first layout */
  *fst_layout_ptr = first_layout_ptr;
  *lst_layout_ptr = last_layout_ptr;

  /* Return number of layouts read in */
  return(nlayouts);
}
/***********************************************************************

add_dummy_layouts  -  Read in the layout as computed by dimhtml

***********************************************************************/

int add_dummy_layouts(struct layout **fst_layout_ptr,
                      struct layout *last_layout_ptr,
                      struct residue **fst_residue_ptr,int nlayouts,
                      float ave_y[2])
{
  char residue_id[20];

  int i, ichain, ndummy, nvalues[2];
  int create_layout, create_residue;

  float y;

  struct layout *layout_ptr, *first_layout_ptr, *other_layout_ptr;
  struct layout *dummy_layout_ptr, *match_layout_ptr;
  struct residue *residue_ptr, *dummy_residue_ptr;
  struct residue *first_residue_ptr, *last_residue_ptr;

  /* Initialise variables */
  ave_y[0] = ave_y[1] = 0;
  ndummy = 0;
  nvalues[0] = nvalues[1] = 0;
  sprintf(residue_id,"[%s %s %c]","DUM","   1 ",'X');

  /* Initialise pointers */
  first_layout_ptr = *fst_layout_ptr;
  first_residue_ptr = *fst_residue_ptr;

  /* Get the first  residue */
  residue_ptr = first_residue_ptr;

  /* Loop to find the last residue pointer */
  while (residue_ptr != NULL)
    {
      /* Save as last */
      last_residue_ptr = residue_ptr;

      /* If residue is in the plot, get its y-coords */
      if (residue_ptr->in_plot == TRUE)
        {
          /* Identify which interface this residue belongs to */
          if (residue_ptr->chain_type == CHAIN1)
            ichain = 0;
          else if (residue_ptr->chain_type == CHAIN2)
            ichain = 1;
          else
            ichain = -1;

          /* Save the mean y-value */
          if (ichain != -1)
            {
              y = (residue_ptr->valmin[MAPPED][Y]
                   + residue_ptr->valmax[MAPPED][Y]) / 2;
              ave_y[ichain] = ave_y[ichain] + y;
              nvalues[ichain]++;
            }
        }

      /* Get the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
    }

  /* Calculate the average y-value */
  for (i = 0; i < 2; i++)
    {
      if (nvalues[i] > 0)
        ave_y[i] = ave_y[i] / nvalues[i];
    }        

  /* Get the pointer for the first layout */
  layout_ptr = first_layout_ptr;

  /* Loop through all the layouts to find the matching ones */
  while (layout_ptr != NULL)
    {
      /* Initialise flags */
      create_layout = create_residue = FALSE;

      /* Get the corresponding residue */
      residue_ptr = layout_ptr->residue_ptr;

      /* If this residue belongs to the first domain, create a mapping
         if it doesn't have one */
      if (layout_ptr->chain_no == 0)
        {
          /* If the residue has no mapping, then set flag to
             create one */
          if (residue_ptr->map_residue_ptr == NULL ||
              residue_ptr->map_residue_ptr->in_plot == FALSE)
            create_residue = TRUE;
        }

      /* Otherwise, for the second domain, check it has a partner */
      else if (layout_ptr->chain_no == 1)
        {
          /* If the residue has no mapping, then link it with its
             layout partner if possible */
          if (residue_ptr->map_residue_ptr == NULL ||
              residue_ptr->map_residue_ptr->in_plot == FALSE)
            {
              /* Set flag to create a dummy residue */
              create_residue = TRUE;

              /* Get the pointers to the first layout again */
              other_layout_ptr = first_layout_ptr;
              match_layout_ptr = NULL;

              /* Loop through all the layouts to find the matching ones */
              while (other_layout_ptr != NULL &&
                     match_layout_ptr == NULL)
                {
                  /* If this is the partner of our current layout record,
                     save it */
                  if (other_layout_ptr->chain_no == 0 &&
                      abs(layout_ptr->y) == abs(other_layout_ptr->y))
                    match_layout_ptr = other_layout_ptr;

                  /* Get the next layout record */
                  other_layout_ptr = other_layout_ptr->next_layout_ptr;
                }
            }
        }

      /* If dummy residue required, then create it */
      if (create_residue == TRUE)
        {
          /* Create a new, dummy residue */
          dummy_residue_ptr
            = create_residue_record(&first_residue_ptr,&last_residue_ptr,
                                    "DUM",residue_ptr->res_num,
                                    residue_ptr->chain,DUMMY,"N/A",
                                    0.0,0.0,residue_id,
                                    residue_ptr->pdbentry_ptr);
          ndummy++;

          /* Set this residue's in-plot flag */
          dummy_residue_ptr->in_plot = TRUE;
          dummy_residue_ptr->chain_type = residue_ptr->chain_type;

          /* Map the current residue to this dummy residue */
          residue_ptr->map_residue_ptr = dummy_residue_ptr;
        }

      /* If dummy layout required, then create it */
      if (create_layout == TRUE)
        {
          /* Create layout record */
          dummy_layout_ptr
            = create_layout_record(&first_layout_ptr,&last_layout_ptr,
                                   nlayouts,0,layout_ptr->x,
                                   -layout_ptr->y,dummy_residue_ptr);

          /* Increment count of layouts */
          nlayouts++;
        }

      /* Get the next layout record */
      layout_ptr = layout_ptr->next_layout_ptr;
    }

  /* Return pointer to first layout and first residue */
  *fst_layout_ptr = first_layout_ptr;
  *fst_residue_ptr = first_residue_ptr;

  /* Return number of dummy residues created */
  return(ndummy);
}
/***********************************************************************

find_layout  -  Loop over the layout records to find the one required

***********************************************************************/

struct layout *find_layout(struct layout *first_layout_ptr,char *res_name,
                           char *res_num,char chain,int search_type,
                           struct layout *other_layout_ptr,int *rp)
{
  int diff, min_diff, rel_pos, y_value;

  struct layout *closest_layout_ptr, *found_layout_ptr, *layout_ptr;

  /* Initialise variables */
  min_diff = 0;
  rel_pos = AFTER;
  y_value = 0;

  /* Initialise pointer */
  closest_layout_ptr = found_layout_ptr = NULL;

  /* If looking for layout partner, then get y-value */
  if (other_layout_ptr != NULL)
    {
      if (search_type == BY_PARTNER)
        y_value = abs(other_layout_ptr->y);
      else
        y_value = other_layout_ptr->y;
    }

  /* Get the pointers for these two layouts */
  layout_ptr = first_layout_ptr;

  /* Loop through all the layouts to find the matching ones */
  while (layout_ptr != NULL && found_layout_ptr == NULL)
    {
      /* See if this is the one we want */
      switch (search_type)
        {
        case BY_RESIDUE :

          /* Check residue name, number and chain */
          if (chain == layout_ptr->residue_ptr->chain &&
              !strcmp(res_num,layout_ptr->residue_ptr->res_num) &&
              !strcmp(res_name,layout_ptr->residue_ptr->res_name))
            found_layout_ptr = layout_ptr;
          break;

        case BY_PARTNER :

          /* If this residue has the same y-value as our residue, then
             must be its partner */
          if (layout_ptr != other_layout_ptr &&
              abs(layout_ptr->y) == y_value &&
              strcmp(layout_ptr->residue_ptr->res_name,"DUM"))
            found_layout_ptr = layout_ptr;
          break;

        case BY_NEAREST :

          /* Skip if this is our layout record or the wrong interface */
          if (layout_ptr != other_layout_ptr &&
              layout_ptr->chain_no == other_layout_ptr->chain_no &&
              layout_ptr->residue_ptr->map_residue_ptr != NULL &&
              strcmp(layout_ptr->residue_ptr->map_residue_ptr->res_name,"DUM"))
            {
              /* Calculate the y-distance between the two residues */
              diff = abs(layout_ptr->y - y_value);

              /* If this is the closest so far, then store */
              if (closest_layout_ptr == NULL || diff < min_diff)
                {
                  /* Save the closest distance */
                  min_diff = diff;
                  closest_layout_ptr = layout_ptr;

                  /* Determine whether it is before or after our
                     residue */
                  if (layout_ptr->y > y_value)
                    rel_pos = BEFORE;
                }
            }
          break;

        default :
          break;
        }

      /* Get the next layout record */
      layout_ptr = layout_ptr->next_layout_ptr;
    }

  /* If searching for the closest residue, save clostest on found */
  if (search_type == BY_NEAREST)
    found_layout_ptr = closest_layout_ptr;

  /* Return whether out residue is before or after */
  *rp = rel_pos;

  /* Return the layout found */
  return(found_layout_ptr);
}
/***********************************************************************

shift_residues  -  Apply given shift to all following residues

***********************************************************************/

int shift_residues(struct pdbentry *first_pdbentry_ptr,
                   struct pdbentry *current_pdbentry_ptr,
                   struct residue *residue1_ptr,
                   struct residue *partner_residue_ptr,float start_x,
                   float shift)
{
  int iatom, natoms;
  int shift_made;

  float ave_x;

  struct atom *atom_ptr;
  struct pdbentry *pdbentry_ptr;
  struct residue *residue_ptr;

  /* Initialise variables */
  shift_made = FALSE;

  /* Get the first PDB entry */
  pdbentry_ptr = first_pdbentry_ptr;

  /* Loop through the PDB entries */
  while (pdbentry_ptr != NULL)
    {
      /* Skip if this is the current PDB entry */
      if (pdbentry_ptr != current_pdbentry_ptr)
        {
          /* Get pointer to the first residue */
          residue_ptr = pdbentry_ptr->first_residue_ptr;

          /* Loop over the residues to shift any that need shifting */
          while (residue_ptr != NULL)
            {
              /* If this isn't one of the two residues just placed, then
                 see if it needs shifting */
              if (residue_ptr != residue1_ptr &&
                  residue_ptr != partner_residue_ptr &&
                  residue_ptr->in_plot == TRUE)
                {
                  /* Get the average x-value of this residue */
                  ave_x = (residue_ptr->valmin[MAPPED][X]
                           + residue_ptr->valmax[MAPPED][X]) / 2;

                  /* If this is a following residue, then apply shift
                     to all its atoms */
                  if (ave_x > start_x)
                    {
                      /* Shift the residue min and max values */
                      residue_ptr->valmin[MAPPED][X]
                        = residue_ptr->valmin[MAPPED][X] + shift;
                      residue_ptr->valmax[MAPPED][X]
                        = residue_ptr->valmax[MAPPED][X] + shift;

                      /* Get the residue's first atom */
                      atom_ptr = residue_ptr->first_atom_ptr;
                      natoms = residue_ptr->natoms;
                      iatom = 0;

                      /* Loop over the residue's atoms to apply shift */
                      while (atom_ptr != NULL && iatom < natoms)
                        {
                          /* Save the minimum x-value */
                          atom_ptr->coords[X] = atom_ptr->coords[X]
                            + shift;
 
                          /* Get the next atom record and increment
                             atom count */
                          atom_ptr = atom_ptr->next_atom_ptr;
                          iatom++;
                        }

                      /* Mark residue as shifted */
                      residue_ptr->shifted = TRUE;
                      shift_made = TRUE;
                    }
                }

              /* Get the next residue */
              residue_ptr = residue_ptr->next_residue_ptr;
            }
        }

      /* Get the next entry */
      pdbentry_ptr = pdbentry_ptr->next_pdbentry_ptr;
    }

  /* Return whether any residues were shifted */
  return(shift_made);
}
/* v.2.1 --> */
/***********************************************************************

get_no  -  Get the residue number, taking into account insertion code

***********************************************************************/

float get_no(char *res_num)
{
  char *string_ptr;

  float extra, res_no;

  /* Initialise residue number */
  res_no = 0;

  /* Get the numeric part */
  res_no = atof(res_num);

  /* If have an insertion code, then add a fraction to the residue
     number */
  if (res_num[4] != ' ')
    {
      /* Get which insertion code this is */
      string_ptr = strchr(INSERT_LIST,res_num[4]);
      if (string_ptr != NULL)
        {
          extra = (string_ptr - INSERT_LIST + 1) / 100.0;
          res_no = res_no + extra;
        }
    }

  /* Return the residue number */
  return(res_no);
}
/* <-- v.2.1 */
/***********************************************************************

place_residue  -  Find the best place to put this residue and return the
                  starting x-coordinate

***********************************************************************/

float place_residue(struct residue *first_residue_ptr,int nresid,
                    struct residue *current_residue_ptr,
                    struct residue *map_residue_ptr,int rel_pos,
                    int order)
{
  char chain;


/* v.2.1 --> */
//  int gap, gap_after, gap_before, iresid, other_res_no, res_no;
  int iresid;
/* <-- v.2.1 */

/* v.2.1 --> */
//  float start_x;
  float gap, gap_after, gap_before, other_res_no, res_no;
  float start_x, width;
/* <-- v.2.1 */

  struct residue *residue_ptr, *after_residue_ptr, *before_residue_ptr;

  /* If residues in order, and this is the first interface chain,
     then see if there is a better position for this residue */
/* v.2.1--> */
//  if (order == TRUE && current_residue_ptr->chain_type == CHAIN1)
  if (order == TRUE && current_residue_ptr->chain_type == CHAIN1)
/* <-- v.2.1 */
    {
      /* Get the residue number of the current residue */
/* v.2.1 --> */
//      res_no = atoi(current_residue_ptr->res_num);
      res_no = get_no(current_residue_ptr->res_num);
/* <-- v.2.1 */
      chain = current_residue_ptr->chain;

      /* Initialise possible residues that might be closest to ours */
      gap_after = -1;
      gap_before = -1;
      after_residue_ptr = before_residue_ptr = NULL;

        /* Get the first residue */
      residue_ptr = first_residue_ptr;
      iresid = 0;

      /* Loop over the residues to find the closest before and after
         ours */
      while (residue_ptr != NULL && iresid < nresid)
        {
          /* If this belongs to the same chain as our residue, then
             consider it */
          if (residue_ptr->chain_type == current_residue_ptr->chain_type &&
              residue_ptr->in_plot == TRUE &&
              (strcmp(residue_ptr->res_name,"DUM") ||
               residue_ptr->placed == TRUE))
            {
              /* Get the current residue number */
/* v.2.1 --> */
//              other_res_no = atoi(residue_ptr->res_num);
              other_res_no = get_no(residue_ptr->res_num);
/* <-- v.2.1 */

              /* Calculate the gap between this residue and the current
                 residue */
              if (residue_ptr->chain == chain)
                gap = other_res_no - res_no;

              /* If this is not the same chain, then compute the gap
                 in a different way */
              else
                {
                  /* If our chain comes first, use the residue number,
                   plus an offset */
                  if (chain < residue_ptr->chain)
                    gap = CHAIN_OFFSET + other_res_no;

                  /* Otherwise, invert the calculation */
                  else
                    gap = - (CHAIN_OFFSET - other_res_no);
                }

              /* If residue is before ours, see if the smallest
                 before gap */
              if (gap < 0)
                {
                  /* Change the sign of the gap */
                  gap = -gap;

                  /* If this is the first, then save */
                  if (gap_before == -1 || gap < gap_before)
                    {
                      gap_before = gap;
                      before_residue_ptr = residue_ptr;
                    }
                }

              /* Repeat for after residue */
              else if (gap > 0)
                {
                  /* If this is the first, then save */
                  if (gap_after == -1 || gap < gap_after)
                    {
                      gap_after = gap;
                      after_residue_ptr = residue_ptr;
                    }
                }
            }
    
          /* Get the next residue */
          residue_ptr = residue_ptr->next_residue_ptr;
          iresid++;
        }

      /* Determine which is closest, the before or after */
      if (gap_before != -1 && gap_after != -1)
        {
          /* Get the smaller gap */
          if (gap_before < gap_after)
            {
              /* Set the before residue */
              map_residue_ptr = before_residue_ptr;
              rel_pos = AFTER;
            }
          else
            {
              /* Set the after residue */
              map_residue_ptr = after_residue_ptr;
              rel_pos = BEFORE;
            }
        }

      /* If just have the before, then set that */
      else if (gap_before != -1)
        {
          /* Set the before residue */
          map_residue_ptr = before_residue_ptr;
          rel_pos = AFTER;
        }

      /* Ditto for after residue */
      else if (gap_after != -1)
        {
          /* Set the before residue */
          map_residue_ptr = after_residue_ptr;
          rel_pos = BEFORE;
        }
    }

/* v.2.1--> */
  /* If map residue has no width, use width of non-bonded contact */
  width = 0;
  if (fabsf(map_residue_ptr->valmax[MAPPED][X]
            - map_residue_ptr->valmin[MAPPED][X]) < 0.1)
    width = SIZE_NONBOND / 2;
/* <-- v.2.1 */

  /* Get the starting x-coordinate */
  if (rel_pos == BEFORE)
/* v.2.1--> */
//    start_x = map_residue_ptr->valmin[MAPPED][X] - IDEAL_GAP / 2;
    start_x = map_residue_ptr->valmin[MAPPED][X] - IDEAL_GAP / 2 - width;
/* <-- v.2.1 */
  else
/* v.2.1--> */
//    start_x = map_residue_ptr->valmax[MAPPED][X] + IDEAL_GAP / 2;
    start_x = map_residue_ptr->valmax[MAPPED][X] + IDEAL_GAP / 2 + width;
/* <-- v.2.1 */

  /* Return the starting x-coordinate */
  return(start_x);
}
/* v.2.1--> */
/***********************************************************************

adjust_atoms  -  Adjust all atom positions in the given residue to
                 correspond to the residue's shift

***********************************************************************/

void adjust_atoms(struct residue *residue_ptr)
{
  int iatom, natoms;

  float centre_x, centre_y;
  float max_x, max_y, min_x, min_y, shift_x, shift_y;

  struct atom *atom_ptr;

  /* Initialise variables */
  natoms = 0;
  max_x = max_y = min_x = min_y = 0;
 
  /* Work out the desired centre-point of the residue */
  centre_x = (residue_ptr->valmin[MAPPED][X]
              + residue_ptr->valmax[MAPPED][X]) / 2;
  centre_y = (residue_ptr->valmin[MAPPED][Y]
              + residue_ptr->valmax[MAPPED][Y]) / 2;

  /* Get the residue's first atom */
  atom_ptr = residue_ptr->first_atom_ptr;
  natoms = residue_ptr->natoms;
  iatom = 0;

  /* Loop over the residue's atoms to calculate
     the average atom position */
  while (atom_ptr != NULL && iatom < natoms)
    {
      /* Save the minimum and maximum values */
      if (iatom == 0)
        {
          min_x = max_x = atom_ptr->coords[X];
          min_y = max_y = atom_ptr->coords[Y];
        }

      /* Update maximum and minimum values */
      else
        {
          if (atom_ptr->coords[X] < min_x)
            min_x = atom_ptr->coords[X];
          if (atom_ptr->coords[X] > max_x)
            max_x = atom_ptr->coords[X];
          if (atom_ptr->coords[Y] < min_y)
            min_y = atom_ptr->coords[Y];
          if (atom_ptr->coords[Y] > max_y)
            max_y = atom_ptr->coords[Y];
        }

      /* Get the next atom record and increment
         atom count */
      atom_ptr = atom_ptr->next_atom_ptr;
      iatom++;
    }

  /* Calculate the shifts to be applied to each atom */
  shift_x = centre_x - (min_x + max_x) / 2;
  shift_y = centre_y - (min_y + max_y) / 2;

  /* Get the residue's first atom again */
  atom_ptr = residue_ptr->first_atom_ptr;
  natoms = residue_ptr->natoms;
  iatom = 0;

  /* Loop over the residue's atoms to calculate
     the average atom position */
  while (atom_ptr != NULL && iatom < natoms)
    {
      /* Apply the shift */
      atom_ptr->coords[X] = atom_ptr->coords[X] + shift_x;
      atom_ptr->coords[Y] = atom_ptr->coords[Y] + shift_y;

      /* Get the next atom record and increment atom count */
      atom_ptr = atom_ptr->next_atom_ptr;
      iatom++;
    }
}
/***********************************************************************
  
list_residues -  List the residues in their current layout positions

***********************************************************************/

void list_residues(struct pdbentry *first_pdbentry_ptr)
{
  int ires, ires1, ires2, nres;
  int done;

  float mid, min1, min2;

  struct pdbentry *pdbentry_ptr;
  struct residue *residue_ptr;

  struct residue **residue_index_ptr;
  
  /* Get the first PDB entry */
  pdbentry_ptr = first_pdbentry_ptr;

  /* Loop through the PDB entries */
  while (pdbentry_ptr != NULL)
    {
      /* Show the PDB entry being processed */
      printf("\n");
      printf("PDB entry - %s\n",pdbentry_ptr->pdb_code);

      /* Get the number of residues */
      nres = pdbentry_ptr->nresid;

      /* Create an index array to store the residues */
      residue_index_ptr
        = (struct residue **) malloc(nres * sizeof(struct residue *));
      if (residue_index_ptr == NULL)
        {
          printf("*** ERROR. Unable to allocate memory for residue index\n");
          exit(1);
        }

      /* Get this entry's first residue */
      residue_ptr = pdbentry_ptr->first_residue_ptr;
      ires = 0;

      /* Loop over the residues to save in the index */
      while (residue_ptr != NULL)
        {
          /* If not over the limit, then save */
          if (residue_ptr->in_plot == TRUE && ires < nres)
            {
              /* Save and increment count */
              residue_index_ptr[ires] = residue_ptr;
              ires++;
            }

          /* Get pointer to next residue record */
          residue_ptr = residue_ptr->next_residue_ptr;
        }

      /* Save the number of residues */
      nres = ires;

      /* Initialise flag */
      done = FALSE;
      
      /* Loop until all residues have been printed */
      while (done == FALSE)
        {
          /* Initialise first residues encountered */
          ires1 = ires2 = -1;
          min1 = min2 = 0;

          /* Loop to find the CHAIN1 and CHAIN2 residues with the minimum
             x-value */
          for (ires = 0; ires < nres; ires++)
            {
              /* Get this residue */
              residue_ptr = residue_index_ptr[ires];

              /* If CHAIN1 residue, then process */
              if (residue_ptr != NULL &&
                  residue_ptr->chain_type == CHAIN1)
                {
                  /* If this is the first, then save x-coord */
                  if (ires1 == -1)
                    {
                      /* Save the minimum valoues and this position */
                      min1 = (residue_ptr->valmin[MAPPED][X]
                              + residue_ptr->valmax[MAPPED][X]) / 2;
                      ires1 = ires;
                    }

                  /* Otherwise save is lowest so far */
                  else
                    {
                      /* Calculate the mid-point */
                      mid = (residue_ptr->valmin[MAPPED][X]
                             + residue_ptr->valmax[MAPPED][X]) / 2;
                      if (mid < min1)
                        {
                          /* Save the minimum valoues and this position */
                          min1 = mid;
                          ires1 = ires;
                        }
                    }
                }

              /* Repeat for CHAIN2 */
              else if (residue_ptr != NULL &&
                       residue_ptr->chain_type == CHAIN2)
                {
                  /* If this is the first, then save x-coord */
                  if (ires2 == -1)
                    {
                      /* Save the minimum valoues and this position */
                      min2 = (residue_ptr->valmin[MAPPED][X]
                              + residue_ptr->valmax[MAPPED][X]) / 2;
                      ires2 = ires;
                    }

                  /* Otherwise save is lowest so far */
                  else
                    {
                      /* Calculate the mid-point */
                      mid = (residue_ptr->valmin[MAPPED][X]
                             + residue_ptr->valmax[MAPPED][X]) / 2;
                      if (mid < min2)
                        {
                          /* Save the minimum valoues and this position */
                          min2 = mid;
                          ires2 = ires;
                        }
                    }
                }
            }

          /* If have both chain 1 and 2 residues, see if they are to be
             printed side by side */
          if (ires1 > -1 && ires2 > -1)
            {
              /* If minimum values are too far apart, just use whichever
                 is the smaller */
              if (fabsf(min2 - min1) > 0.5)
                {
                  /* Take the smaller */
                  if (min1 < min2)
                    ires2 = -1;
                  else
                    ires1 = -1;
                }
            }

          /* If first is available, print it */
          if (ires1 > -1)
            {
              /* Get this residue */
              residue_ptr = residue_index_ptr[ires1];

              /* Print this residue */
              printf("%8.3f %8.3f %s %s %c    ",
                     residue_ptr->valmin[MAPPED][X],
                     residue_ptr->valmax[MAPPED][X],
                     residue_ptr->res_name,
                     residue_ptr->res_num,
                     residue_ptr->chain);

              /* Blank out its pointer in the array */
              residue_index_ptr[ires1] = NULL;
            }
          else
            printf("                                 ");

          /* Repeat for second residue */
          if (ires2 > -1)
            {
              /* Get this residue */
              residue_ptr = residue_index_ptr[ires2];

              /* Print this residue */
              printf("%8.3f %8.3f %s %s %c    ",
                     residue_ptr->valmin[MAPPED][X],
                     residue_ptr->valmax[MAPPED][X],
                     residue_ptr->res_name,
                     residue_ptr->res_num,
                     residue_ptr->chain);

              /* Blank out its pointer in the array */
              residue_index_ptr[ires2] = NULL;
            }

          /* If either was printed, then close off line */
          if (ires1 > -1 || ires2 > -1)
            printf("\n");

          /* If neither found, then we are done */
          if (ires1 == -1 && ires2 == -1)
            done = TRUE;
        }

      /* Free the memory taken by the residue index */
      free(residue_index_ptr);

      /* Get the next entry */
      pdbentry_ptr = pdbentry_ptr->next_pdbentry_ptr;
    }

  /* Throw blank line */
  printf("\n");
}
/* <-- v.2.1 */
/***********************************************************************

arrange_dimplot  -  Use the fitting to determine the atom and residue
                    equivalences for use in restraining the second plot

***********************************************************************/

int arrange_dimplot(struct pdbentry *ref_pdbentry_ptr,
                    struct pdbentry *pdbentry_ptr,
                    struct layout *first_layout_ptr,float ave_y[2],
                    int order)
{
  char blank[2];

  int rel_pos, shifted, shift_made;
/* v.2.1--> */
  int verbose;
/* <-- v.2.1 */

/* v.2.1--> */
//  float shift, start_x;
  float shift, start_x, width;
/* <-- v.2.1 */

  struct layout *current_layout_ptr, *partner_layout_ptr;
  struct layout *closest_layout_ptr;
  struct residue *residue1_ptr, *residue2_ptr;
  struct residue *closest_residue_ptr, *partner_residue_ptr;
  struct residue *map_residue_ptr;

  /* Initialise variables */
  blank[0] = '\0';
  shifted = FALSE;
/* v.2.1--> */
  verbose = FALSE;
/* <-- v.2.1 */


  /* Get the first of the PDB entry's residues */
  residue1_ptr = pdbentry_ptr->first_residue_ptr;

  /* Loop over the residues to map equivalent atoms or residues */
  while (residue1_ptr != NULL)
    {
      /* Only worry about this residue if it is in the plot */
      if (residue1_ptr->in_plot == TRUE)
        {
          /* Get the residue this maps to in the original plot */
          residue2_ptr = residue1_ptr->map_residue_ptr;

          /* If residue has no mapping, show warning message */
          if (residue2_ptr == NULL)
            printf("*** Warning. Residue %s %s %c has no mapping to "
                   "reference residue\n",
                   residue1_ptr->res_name,
                   residue1_ptr->res_num,
                   residue1_ptr->chain);

          /* If this residue maps to a dummy residue, then find where
             to place it, shifting everything else to make way for it.
             Ditto if H-group maps to a hydrophobic group */
          else if (residue2_ptr->placed == FALSE &&
                   (!strcmp(residue2_ptr->res_name,"DUM") ||
                    (residue2_ptr->interaction == PROT_PROT_NONBOND &&
                     residue1_ptr->interaction == PROT_PROT_HBOND)))
            {
              /* Get the layout record for this residue */
              current_layout_ptr
                = find_layout(first_layout_ptr,residue1_ptr->res_name,
                              residue1_ptr->res_num,residue1_ptr->chain,
                              BY_RESIDUE,NULL,&rel_pos);

              /* If found, then process */
              if (current_layout_ptr != NULL)
                {
                  /* Find this residue's partner residue in the layout */
                  partner_layout_ptr
                    = find_layout(first_layout_ptr,blank,blank,' ',
                                  BY_PARTNER,current_layout_ptr,&rel_pos);

                  /* If we have a partner residue, get it */
                  if (partner_layout_ptr != NULL)
                    {
                      /* Get the partner residue */
                      partner_residue_ptr
                        = partner_layout_ptr->residue_ptr;
                    }
                  else
                    partner_residue_ptr = NULL;

/* v.2.1--> */
                  /* If mapping to a real residue, rather than a dummy,
                     take this as the closest layout position */
                  if (strcmp(residue2_ptr->res_name,"DUM"))
                    closest_layout_ptr = current_layout_ptr;

                  /* Find the closest mapped residue to it in the layout */
//                  if (residue1_ptr->chain_type == CHAIN1 ||
                  else if (residue1_ptr->chain_type == CHAIN1 ||
/* <-- v.2.1 */
                      partner_residue_ptr == NULL ||
                      partner_residue_ptr->map_residue_ptr != NULL)
                    closest_layout_ptr
                      = find_layout(first_layout_ptr,blank,blank,' ',
                                    BY_NEAREST,current_layout_ptr,&rel_pos);
                  else
                    closest_layout_ptr = NULL;

                  /* If closest found, then place just before or after it */
                  if (closest_layout_ptr != NULL)
                    {
                      /* Get the closest mapped residue */
                      closest_residue_ptr
                        = closest_layout_ptr->residue_ptr;

                      /* Get the residue it maps to */
                      map_residue_ptr
                        = closest_layout_ptr->residue_ptr->map_residue_ptr;

                      /* Get the starting x-coordinate */
/* v.2.1--> */
                       if (!strcmp(residue2_ptr->res_name,"DUM"))
/* <-- v.2.1 */
                         start_x
                           = place_residue(ref_pdbentry_ptr->first_residue_ptr,
                                           ref_pdbentry_ptr->nresid,
                                           residue1_ptr,map_residue_ptr,
                                           rel_pos,order);
/* v.2.1--> */
                       else
                         {
                           /* If residue has no width, use width of
                              non-bonded contact */
                           width = 0;
                           if (fabsf(map_residue_ptr->valmax[MAPPED][X]
                                     - map_residue_ptr->valmin[MAPPED][X])
                               < 0.1)
                             width = SIZE_NONBOND / 2;
                           start_x
                             = map_residue_ptr->valmin[MAPPED][X]
                             - IDEAL_GAP / 2 - width;
                         }
/* <-- v.2.1 */

                      /* Determine the size of the shift to be made */
                      if (!strcmp(residue2_ptr->res_name,"DUM"))
                        {
                          /* If residue, or its partner, is an H-bonded
                             group, allow more space for fitting */
                          if (residue1_ptr->interaction == PROT_PROT_HBOND ||
                              (partner_residue_ptr != NULL &&
                               partner_residue_ptr->interaction
                               == PROT_PROT_HBOND))
                            shift = SIZE_HBOND;
                          else
                            shift = SIZE_NONBOND;
                        }

                      /* Otherwise, we have an H-group replacing a
                         hydrophobic group, so allow a little extra space */
                      else
                        shift = SIZE_HBOND - SIZE_NONBOND;

                      /* Place the current residue at this position */
                      residue2_ptr->valmin[MAPPED][X]
/* v.2.1--> */
//                        = start_x;
                        = start_x + IDEAL_GAP / 2;
/* <-- v.2.1 */
                      residue2_ptr->valmax[MAPPED][X]
/* v.2.1--> */
//                        = start_x + shift;
                        = start_x + shift + IDEAL_GAP / 2;

                      /* Save the residue's name and set flags */
                      strcpy(residue2_ptr->res_name,residue1_ptr->res_name);
                      residue2_ptr->in_plot = TRUE;
                      residue2_ptr->placed = TRUE;
/* <-- v.2.1 */

                      /* Fix its y-position */
                      residue2_ptr->valmin[MAPPED][Y]
                        = residue2_ptr->valmax[MAPPED][Y]
                        = ave_y[current_layout_ptr->chain_no];
                      residue2_ptr->placed = TRUE;

/* v.2.1--> */
                      /* Adjust all atom positions to correspond to this
                         shift */
                      adjust_atoms(residue2_ptr);
/* <-- v.2.1 */

                      /* Reinitialise pointers */
                      map_residue_ptr = NULL;

                      /* Repeat for its partner, if it has one */
                      if (partner_layout_ptr != NULL)
                        {
                          /* If valid, and isn't mapped to a residue
                             in the reference plot, then adjust its
                             position */
                          if (partner_residue_ptr != NULL &&
                              partner_residue_ptr->map_residue_ptr != NULL &&
                              !strcmp(partner_residue_ptr->map_residue_ptr->
                                      res_name,"DUM"))
                            {
                              /* Get the dummy residue mapped to */
                              map_residue_ptr
                                = partner_residue_ptr->map_residue_ptr;
                            }
                        }


                      /* Move the layout partner residue, if there is one */
                      if (map_residue_ptr != NULL)
                        {
                          map_residue_ptr->valmin[MAPPED][X]
                            = start_x;
                          map_residue_ptr->valmax[MAPPED][X]
                            = start_x + shift;
                          map_residue_ptr->valmin[MAPPED][Y]
                            = map_residue_ptr->valmax[MAPPED][Y]
                            = ave_y[partner_layout_ptr->chain_no];
                          map_residue_ptr->placed = TRUE;
                        }

/* v.2.1--> */
                      /* Increas the shift for the following residues */
                      shift = shift + IDEAL_GAP;
/* <-- v.2.1 */

                      /* Shift all following residues along */
                      shift_made
                        = shift_residues(ref_pdbentry_ptr,pdbentry_ptr,
                                         residue2_ptr,map_residue_ptr,
                                         start_x,shift);

                      /* Update shift flag */
                      if (shift_made == TRUE)
                        shifted = TRUE;

                      /* Update all the maximum and minimum values */
                      update_all_maxmin(pdbentry_ptr);

/* v.2.1--> */
                      /* If debugging, show the current layout of the
                         residues */
                      if (verbose == TRUE)
                        list_residues(pdbentry_ptr);
/* <-- v.2.1 */
                    }
                }
            }
        }

      /* Get the next residue */
      residue1_ptr = residue1_ptr->next_residue_ptr;
    }

  /* Return whether any reference residues were shifted */
  return(shifted);
}
/***********************************************************************

dimplot_restraints  -  Write out the dimplot restraints

***********************************************************************/

void dimplot_restraints(struct pdbentry *ref_pdbentry_ptr,
                        struct pdbentry *pdbentry_ptr,
                        int fit_type,struct parameters params)
{
  int atom_number, loop, nmap;
  int done, overlay;
  int all_atom1, all_atom2;

  struct pdbentry *pdbentry2_ptr;
  struct residue *residue1_ptr, *residue2_ptr;

  FILE *file_ptr, *lgffile_ptr;

  /* Initialise variables */
  atom_number = 0;
  loop = 1;
  overlay = FALSE;

  /* Open the appropriate output file */
  file_ptr = open_file("dimplot.rcm",params);
  lgffile_ptr = open_file("ligplus.lgf",params);

  /* Write out the header records to the .rcm file */
  fprintf(file_ptr,"                 Res.  Res       Flattened    \n");
  fprintf(file_ptr,"            Atom Name  Num      coordinates   \n");
  fprintf(file_ptr,"            ---- --- -----    ----------------\n");

  /* Write out how the fit was made to the .lgf file */
  fprintf(lgffile_ptr,"FIT_TYPE=%s\n",FIT_TYPE[fit_type]);
  fprintf(lgffile_ptr,"ALIGNED=Chain %c of %s and chain %c of %s\n",
          ref_pdbentry_ptr->map_chain,ref_pdbentry_ptr->pdb_code,
          ref_pdbentry_ptr->map_chain,pdbentry_ptr->pdb_code);

  /* Get the first of the PDB entry's residues */
  residue1_ptr = pdbentry_ptr->first_residue_ptr;

  /* Loop over the residues to map equivalent atoms or residues */
  while (residue1_ptr != NULL)
    {
      /* Initialise flag */
      done = FALSE;

      /* If this residue has a mapping, then write out */
      if (residue1_ptr->in_plot == TRUE &&
          (residue2_ptr = residue1_ptr->map_residue_ptr) != NULL &&
          residue2_ptr->in_plot == TRUE && done == FALSE)
        {
          /* Initialise count of atoms written out */
          nmap = 0;

          /* Determine whether both residues are all-atom */
          all_atom1 = all_atom2 = FALSE;
          if (residue1_ptr->interaction == PROT_PROT_HBOND ||
              residue1_ptr->interaction == VIA_WATER ||
              residue1_ptr->interaction == VIA_METAL)
            all_atom1 = TRUE;
          if (residue2_ptr->interaction == PROT_PROT_HBOND ||
              residue2_ptr->interaction == VIA_WATER ||
              residue2_ptr->interaction == VIA_METAL)
            all_atom2 = TRUE;

          /* If both residues are H-bond groups, write out all atom
             coords */
          if (all_atom1 == TRUE && all_atom2 == TRUE)
            {
              /* Write out the all-atom positional restraints for this
                 residue */
              nmap = write_equiv_residues(file_ptr,residue1_ptr,
                                          residue2_ptr,loop);
            }

          /* If no atoms written out, then write out centre of mass
             instead */
          if (nmap == 0)
            write_residue_cofm(file_ptr,residue1_ptr,residue2_ptr,loop);

          /* Get other PDB entry */
          pdbentry2_ptr = residue2_ptr->pdbentry_ptr;

          /* If have PDB and residue identifiers, write out the mapping */
          if (pdbentry2_ptr != NULL && pdbentry2_ptr->pdb_id > -1 &&
              residue2_ptr->residue_id != &null)
            {
              /* Write the mapping */
              fprintf(lgffile_ptr,"%d\t%s\t[%s %s %c]\n",
                      pdbentry2_ptr->pdb_id,residue2_ptr->residue_id,
                      residue1_ptr->res_name,residue1_ptr->res_num,
                      residue1_ptr->chain);
            }
        }

      /* Get the next residue */
      residue1_ptr = residue1_ptr->next_residue_ptr;
    }

  /* Close the output files */
  fclose(file_ptr);
  fclose(lgffile_ptr);
}
/* <-- v.2.0 */
/***********************************************************************

write_ligand_files  -  Write out the ligand.pdb, .hhb and .nnb files

***********************************************************************/

/* v.2.0.3 --> */
//void write_ligand_files(struct residue *first_ligand_residue_ptr,
//                        struct parameters params)
void write_ligand_files(struct residue *first_ligand_residue_ptr,
                        struct bond *first_bond_ptr,
                        struct parameters params)
/* <-- v.2.0.3 */
{
  int iatom, natoms;
/* v.2.0.3 --> */
  int min_atno, max_atno;
/* <-- v.2.0.3 */
  int done, overlay;

  struct atom *atom_ptr;
/* v.2.0.3 --> */
  struct bond *bond_ptr;
/* <-- v.2.0.3 */
  struct residue *residue_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  done = overlay = FALSE;
/* v.2.0.3 --> */
  min_atno = max_atno = -1;
/* <-- v.2.0.3 */

  /* Open the ligand coordinate file */
  file_ptr = open_file("ligand.pdb",params);

  /* Get the first ligand residue */
  residue_ptr = first_ligand_residue_ptr;

  /* Loop over the ligand residues */
  while (residue_ptr != NULL && residue_ptr->interaction == IS_LIGAND)
    {
      /* Get the residue's first atom */
      atom_ptr = residue_ptr->first_atom_ptr;
      natoms = residue_ptr->natoms;
      iatom = 0;

      /* Loop over the residue's atoms to write out ligand coords */
      while (atom_ptr != NULL && iatom < natoms)
        {
          /* Write out */
          fprintf(file_ptr,"ATOM  %5d %s %s %c%s   "
                  "%8.3f%8.3f%8.3f%6.2f%6.2f          %s\n",
                  atom_ptr->atom_number,atom_ptr->atom_name,
                  residue_ptr->res_name,
                  residue_ptr->chain,
                  residue_ptr->res_num,
                  atom_ptr->original_coords[X],
                  atom_ptr->original_coords[Y],
                  atom_ptr->original_coords[Z],
                  atom_ptr->occupancy,
                  atom_ptr->bvalue,atom_ptr->element);

/* v.2.0.3 --> */
          /* Update ligand-atom range */
          if (min_atno == -1)
            min_atno = max_atno = atom_ptr->atom_number;
          else
            {
              if (atom_ptr->atom_number < min_atno)
                min_atno = atom_ptr->atom_number;
              if (atom_ptr->atom_number > max_atno)
                max_atno = atom_ptr->atom_number;
            }
/* <-- v.2.0.3 */

          /* Get the next atom record and increment atom count */
          atom_ptr = atom_ptr->next_atom_ptr;
          iatom++;
        }

      /* Get the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
    }

/* v.2.0.3 --> */
  /* Get the current fragment's first atom */
  bond_ptr = first_bond_ptr;

  /* Loop over the fragment's bonds to write out */
  while (bond_ptr != NULL)
    {
      /* Check that both this bond's atoms are from the ligand residue */
      if (bond_ptr->type == COVALENT &&
          bond_ptr->atom1_ptr->residue_ptr->interaction == IS_LIGAND &&
          bond_ptr->atom2_ptr->residue_ptr->interaction == IS_LIGAND)
        {
          /* Check that both atom numbers fall within the range of our
             ligand */
          if (bond_ptr->atom1_ptr->atom_number >= min_atno &&
              bond_ptr->atom1_ptr->atom_number <= max_atno &&
              bond_ptr->atom2_ptr->atom_number >= min_atno &&
              bond_ptr->atom2_ptr->atom_number <= max_atno)
            {
              /* Write out this bond as a CONECT record */
              fprintf(file_ptr,"CONECT%5d%5d\n",
                      bond_ptr->atom1_ptr->atom_number,
                      bond_ptr->atom2_ptr->atom_number);
            }
        }

      /* Get the next bond */
      bond_ptr = bond_ptr->next_bond_ptr;
    }
/* <-- v.2.0.3 */
  
  /* Close the PDB file */
  fclose(file_ptr);

  /* Create empty .hhb file */
  file_ptr = open_file("ligand.hhb",params);
  fprintf(file_ptr,"ligplot.hhb output:\n");
  fprintf(file_ptr,"\n");
  fprintf(file_ptr,"   Donor                  Acceptor     Distance\n");
  fclose(file_ptr);

  /* Repeat for .nnb file */
  file_ptr = open_file("ligand.nnb",params);
  fprintf(file_ptr,"ligplot.nnb output:\n");
  fprintf(file_ptr,"\n");
  fprintf(file_ptr,"   Atom 1               Atom 2        Distance\n");
  fclose(file_ptr);
}
/***********************************************************************

perform_splice  -  Splice out given string from line and replace with
                   new string

***********************************************************************/

int perform_splice(char *line,char *string,char *new_string)
{
  char *tmp;

  int end_pos, ipos, len, len_new, len_string;
  int done;

  /* Initialise variables */
  done = FALSE;
  end_pos = 0;

  /* Get the length of the input line */
  len = strlen(line);
  len_string = strlen(string);
  len_new = strlen(new_string);

  /* If any string empty, then return */
  if (len == 0 || len_string == 0)
    return(len);

  /* Loop through the line searching for the replace-string */
  for (ipos = 0; ipos < len - len_string + 1 && done == FALSE; ipos++)
    {
      /* Check the string at the current position */
      if (!strncmp(line+ipos,string,len_string))
        {
          /* Create a temporary string to hold the spliced string */
          tmp = (char *) malloc(sizeof(char)*(len + len_new + 1));

          /* Initialise the temporary string */
          tmp[0] = '\0';

          /* Transfer first part of line */
          strncpy(tmp,line,ipos);
          tmp[ipos] = '\0';

          /* Splice in the replacement path */
          strcat(tmp,new_string);

          /* Store end position of splice */
          end_pos = ipos + len_new;

          /* Add on the remainder of the line */
          strcat(tmp,line+ipos+len_string);

          /* Transfer back across to the input line */
          strcpy(line,tmp);

          /* Free up memory taken by the temporary string */
          free(tmp);

          /* Set flag that splice done */
          done = TRUE;
        }
    }

  /* Return splice end position */
  return(end_pos);
}
/***********************************************************************

repeat_splice  -  Splice out all occurrences of the given string from
                  line and replace with new string

***********************************************************************/

void repeat_splice(char *line,char *string,char *new_string)
{
  int end_pos;
  int finished;
  
  /* Initialise variables */
  end_pos = 0;
  finished = FALSE;
  if (line[0] == '\0')
    finished = TRUE;

  /* Loop until all replacements have been made */
  while (finished == FALSE)
    {
      /* Perform splice */
      end_pos = perform_splice(line + end_pos,string,new_string);

      /* If no replacements made, we're done */
      if (end_pos == 0)
        finished = TRUE;
    }
}
/***********************************************************************

run_ligplot  -  Run LIGPLOT on the ligand.pdb file

***********************************************************************/

void run_ligplot(struct parameters params)
{
  char command[FILENAME_LEN];
/* v.2.2.8 --> */
  char ch, exe_dir[FILENAME_LEN], tmp[FILENAME_LEN];

  int win;
  
  /* See if running on Windows */
#ifdef _WIN32
  win = TRUE;
#else
  win = FALSE;
#endif  
/* debug */
  printf("WINDOWS = %s\n",T_F[win]);
  fflush(stdout);
/* debug */

  /* Save the executable directory */
  strcpy(exe_dir,params.ligplot_exe_dir);
/* debug */
  printf("In  run_ligplot with exe: %s\n",exe_dir);
  fflush(stdout);
/* debug */
  
/* <-- v.2.2.8 */
/* debug */
  printf("Forming LIGPLOT command\n");
  fflush(stdout);
/* debug */

  /* Form the command to run LIGPLOT */
  if (win == TRUE)
    strcpy(command,"\"\"");
  else
    strcpy(command,"\"");
  strcat(command,exe_dir);
  //  strcat(command,params.ligplot_exe_dir);
  strcat(command,"/ligplot\" \"");
  strcat(command,params.directory);
  strcat(command,DIR_SEP);
  strcat(command,"ligand.pdb\" ");
  strcat(command,params.ligplot_params);
  if (params.directory != &null)
    {
      strcat(command," -wkdir \"");
      strcat(command,params.directory);
      strcat(command,"\"");
      if (params.prmfile != &null)
        {
          strcat(command," -prm \"");
          strcat(command,params.prmfile);
          strcat(command,"\"");
        }
    }
  strcat(command," -ligfit");
  if (win == TRUE)
    strcat(command,"\"");

  /* Write out the LIGPLOT command */
  printf("Ligfit's LIGPLOT command: %s\n",command);
  fflush(stdout);

  /* Run LIGPLOT */
  printf("LIGPLOT command: %s\n",command);
  system(command);
}
/***********************************************************************

check_overlap  -  Check whether the given residue from original plot
                  overlaps the current ligand

***********************************************************************/

int check_overlap(struct residue *residue1_ptr,
                  struct residue *residue2_ptr)
{
  int iatom, jatom, natoms1, natoms2, noverlap;
  int in_range, overlap;

  double dist_sqrd, x, y, z;

  struct atom *atom1_ptr, *atom2_ptr;

  /* Initialise variables */
  overlap = FALSE;
  noverlap = 0;

  /* Get the first atom */
  atom1_ptr = residue1_ptr->first_atom_ptr;
  natoms1 = residue1_ptr->natoms;
  iatom = 0;

  /* Loop over the residue's atoms to compute their transformed coords
     and check for overlap with the atoms of the other residue */
  while (atom1_ptr != NULL && iatom < natoms1)
    {
      /* Get this atom's coords */
      x = atom1_ptr->coords[X];
      y = atom1_ptr->coords[Y];
      z = atom1_ptr->coords[Z];

       /* Check if the atom falls within the LIGPLOT residue's limits */
      in_range = TRUE;
      if (x > residue2_ptr->valmax[MAPPED][X] + 2 * MAX_MAP_RADIUS ||
          x < residue2_ptr->valmin[MAPPED][X] - 2 * MAX_MAP_RADIUS ||
          y > residue2_ptr->valmax[MAPPED][Y] + 2 * MAX_MAP_RADIUS ||
          y < residue2_ptr->valmin[MAPPED][Y] - 2 * MAX_MAP_RADIUS ||
          z > residue2_ptr->valmax[MAPPED][Z] + 2 * MAX_MAP_RADIUS ||
          z < residue2_ptr->valmin[MAPPED][Z] - 2 * MAX_MAP_RADIUS)
        in_range = FALSE;

      /* If in range, see if atom overlaps any of the other residue's
         atoms */
      if (in_range == TRUE)
        {
          /* Get the first atom of the LIGPLOT residue */
          atom2_ptr = residue2_ptr->first_atom_ptr;
          natoms2 = residue2_ptr->natoms;
          jatom = 0;
          overlap = FALSE;

          /* Loop over the residue's atoms */
          while (atom2_ptr != NULL && jatom < natoms2 && overlap == FALSE)
            {
              /* Calc distance between these two atoms */
              dist_sqrd
                = (x - atom2_ptr->coords[X])
                * (x - atom2_ptr->coords[X])
                + (y - atom2_ptr->coords[Y])
                * (y - atom2_ptr->coords[Y])
                + (z - atom2_ptr->coords[Z])
                * (z - atom2_ptr->coords[Z]);

              /* If atoms overlap, then set flag */
              if (dist_sqrd < MAX_MAP_RADIUS2)
                overlap = TRUE;

              /* Get the next atom record and increment atom count */
              atom2_ptr = atom2_ptr->next_atom_ptr;
              jatom++;
            }

          /* If we have an overlap, then increment count of overlapping
             atoms */
          if (overlap == TRUE)
            noverlap++;
        }

      /* Get the next atom record and increment atom count */
      atom1_ptr = atom1_ptr->next_atom_ptr;
      iatom++;
    }

  /* Return number of overlapping atoms */
  return(noverlap);
}
/* v.1.2 --> */
/***********************************************************************

fit_2d_ligands  -  Fit the flattened ligands on their equivalenced
                   atoms

***********************************************************************/

void fit_2d_ligands(struct atomap *first_atomap_ptr,
                    struct pdbentry *ligand_pdbentry_ptr,
                    struct parameters params)
{
  BOOL column;

  COOR coords1[MAX_FIT_RES], coords2[MAX_FIT_RES];

  int i, iatom, j, natoms, nstored;
  int wanted;

  double cofg[2][3], jog_cofg[2][3], rmsd, x, y, z;
  double theta = 60.0 / RADDEG;

  REAL jog_matrix[3][3], matrix[3][3];
                
  struct atom *atom_ptr, *ligand_atom_ptr, *map_atom_ptr;
  struct atom *other_atom_ptr;
  struct atomap *atomap_ptr;
  struct residue *residue_ptr, *map_residue_ptr;
  struct residue *ligand_residue_ptr, *first_ligand_residue_ptr;

  /* Initialise variables */
  column = FALSE;
  nstored = 0;
  rmsd = 0;
  first_ligand_residue_ptr = ligand_pdbentry_ptr->first_residue_ptr;

  /* Initialise transformation matrices */
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      if (i == j)
        jog_matrix[i][j] = matrix[i][j] = 1.0;
      else
        jog_matrix[i][j] = matrix[i][j] = 0.0;
  for (i = 0; i < 2; i++)
    for (j = 0; j < 3; j++)
      jog_cofg[i][j] = cofg[i][j] = 0.0;

  /* Set up the jog matrix to flip ligand out of the x-y plane to make
     it easier to be fitted to the target ligand in the plane */
  jog_matrix[0][0] = cos(theta);
  jog_matrix[0][2] = sin(theta);
  jog_matrix[2][0] = - sin(theta);
  jog_matrix[2][2] = cos(theta);

  /* Get the first of the atom mappings */
  atomap_ptr = first_atomap_ptr;

  /* Loop to pick up ligand atom equivalences */
  while (atomap_ptr != NULL)
    {
      /* Get the two atoms */
      atom_ptr = atomap_ptr->atom_ptr;
      residue_ptr = atom_ptr->residue_ptr;
      map_atom_ptr = atomap_ptr->map_atom_ptr;
      map_residue_ptr = map_atom_ptr->residue_ptr;

      /* If this is a ligand residue, then process */
      if (residue_ptr->interaction == IS_LIGAND)
        {
          /* Get the first ligand residue */
          ligand_residue_ptr = first_ligand_residue_ptr;

          /* Loop over the ligand residues */
          while (ligand_residue_ptr != NULL)
            {
              /* If this residue matches the one we're interested in,
                 then set flag */
              wanted = FALSE;
              if (!strcmp(ligand_residue_ptr->res_num,
                          residue_ptr->res_num) &&
                  !strcmp(ligand_residue_ptr->res_name,
                          residue_ptr->res_name) &&
                  ligand_residue_ptr->chain
                  == residue_ptr->chain)
                wanted = TRUE;

              /* Get the residue's first atom */
              other_atom_ptr = ligand_residue_ptr->first_atom_ptr;
              ligand_atom_ptr = NULL;
              natoms = ligand_residue_ptr->natoms;
              iatom = 0;

              /* Loop over the ligand atoms to find the current one of
                 interest */
              while (wanted == TRUE && other_atom_ptr != NULL &&
                     iatom < natoms && ligand_atom_ptr == NULL)
                {
                  /* If this is the atom we're after, then save */
                  if (!strcmp(other_atom_ptr->atom_name,
                              atom_ptr->atom_name))
                    ligand_atom_ptr = other_atom_ptr;

                  /* Get the next atom record and increment atom count */
                  other_atom_ptr = other_atom_ptr->next_atom_ptr;
                  iatom++;
                }

              /* If match found, then save 2D coords of both atoms */
              if (ligand_atom_ptr != NULL)
                {
                  /* Transform the ligand out of the plane to make
                     it easier to fit back down onto it */
                  transform_atom(ligand_atom_ptr->coords[X],
                                 ligand_atom_ptr->coords[Y],
                                 ligand_atom_ptr->coords[Z],
                                 jog_matrix,jog_cofg,
                                 &x,&y,&z);

                  /* Store the coords */
                  coords1[nstored].x = x;
                  coords1[nstored].y = y;
                  coords1[nstored].z = z;

                  /* Store coords of atom to map to */
                  coords2[nstored].x = map_atom_ptr->coords[X];
                  coords2[nstored].y = map_atom_ptr->coords[Y];
                  coords2[nstored].z = map_atom_ptr->coords[Z];

                  /* Increment count of stored coords */
                  nstored++;
                }

              /* Get the next ligand residue */
              ligand_residue_ptr = ligand_residue_ptr->next_residue_ptr;
            }
        }

      /* Get pointer to next record */
      atomap_ptr = atomap_ptr->next_atomap_ptr;
    }

  /* If fewer than 2 coords stored, then can't fit */
  if (nstored < 2)
    return;

  /* Centre coords of first ligand */
  centre_coords(coords1,nstored,cofg[CURR_SEQUENCE]);

  /* Repeat for second ligand */
  centre_coords(coords2,nstored,cofg[REF_SEQUENCE]);

  /* Perform the least-squares fit */
  rmsd = matfit(coords1,coords2,matrix,nstored,NULL,column);

  /* Get the first ligand residue */
  ligand_residue_ptr = first_ligand_residue_ptr;

  /* Loop over the ligand residues */
  while (ligand_residue_ptr != NULL)
    {
      /* Get the residue's first atom */
      ligand_atom_ptr = ligand_residue_ptr->first_atom_ptr;
      natoms = ligand_residue_ptr->natoms;
      iatom = 0;

      /* Loop over the ligand atoms to transform their coordinates */
      while (ligand_atom_ptr != NULL && iatom < natoms)
        {
          /* Apply the jog transformation to this atom */
          transform_atom(ligand_atom_ptr->coords[X],
                         ligand_atom_ptr->coords[Y],
                         ligand_atom_ptr->coords[Z],
                         jog_matrix,jog_cofg,
                         &x,&y,&z);

          /* Store the coords */
          ligand_atom_ptr->coords[X] = x;
          ligand_atom_ptr->coords[Y] = y;
          ligand_atom_ptr->coords[Z] = z;

          /* Apply transformation to map onto previous ligand */
          transform_atom(ligand_atom_ptr->coords[X],
                         ligand_atom_ptr->coords[Y],
                         ligand_atom_ptr->coords[Z],matrix,cofg,&x,&y,&z);

          /* Save the transformed coords */
          ligand_atom_ptr->coords[X] = x;
          ligand_atom_ptr->coords[Y] = y;
          ligand_atom_ptr->coords[Z] = z;

          /* Get the next atom record and increment atom count */
          ligand_atom_ptr = ligand_atom_ptr->next_atom_ptr;
          iatom++;
        }

      /* Get the next ligand residue */
      ligand_residue_ptr = ligand_residue_ptr->next_residue_ptr;
    }
}
/* <-- v.1.2 */
/***********************************************************************

calc_cofm  -  Calculate the residue's (or ligand's) C of M

***********************************************************************/

void calc_cofm(struct residue *first_residue_ptr,int nresid,double *cofm)
{
  int iatom, iresid, natoms, ntot;

  struct atom *atom_ptr;
  struct residue *residue_ptr;

  /* Initialise variables */
  iresid = ntot = 0;

  /* Initialise centre of mass */
  cofm[X] = cofm[Y] = cofm[Z] = 0.0;

  /* Get the first residue */
  residue_ptr = first_residue_ptr;

  /* Loop over the residues to accumulate C of M */
  while (residue_ptr != NULL && iresid < nresid)
    {
      /* Get the first of this residue's atoms */
      atom_ptr = residue_ptr->first_atom_ptr;
      natoms = residue_ptr->natoms;
      iatom = 0;

      /* Loop over the residue's atoms to calculate centre of mass */
      while (atom_ptr != NULL && iatom < natoms)
        {
          /* Add this atom's coords */
          cofm[X] = cofm[X] + atom_ptr->coords[X];
          cofm[Y] = cofm[Y] + atom_ptr->coords[Y];
          cofm[Z] = cofm[Z] + atom_ptr->coords[Z];

          /* Get the next atom record and increment atom count */
          atom_ptr = atom_ptr->next_atom_ptr;
          iatom++;
          ntot++;
        }
    
      /* Get the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;
    }
          
  /* Calculate centre of mass and write out */
  if (ntot > 0)
    {
      cofm[X] = cofm[X] / ntot;
      cofm[Y] = cofm[Y] / ntot;
      cofm[Z] = cofm[Z] / ntot;
    }
}
/***********************************************************************

sign  -  Return the sign of the given value

***********************************************************************/

int sign(double value)
{
  int isign;

  /* Determine if +ve or -ve */
  isign = 2 * (value >= 0) - 1;

  return(isign);
}
/***********************************************************************

make_shift  -  Apply shift to given side chain away from ligand position

***********************************************************************/

void make_shift(struct residue *residue_ptr, double *cofm,
                double *ligand_cofm)
{
  int i, iatom, natoms;

  double dist, shift[2], shift_fraction;

  struct atom *atom_ptr;

  /* Initialise variables */
  dist = 0;
  shift[X] = shift[Y] = 0;

  /* Calculate the difference in the centres of mass */
  for (i = 0; i < 2; i++)
    dist = dist
      + (ligand_cofm[i] - cofm[i]) * (ligand_cofm[i] - cofm[i]);
  dist = sqrt(dist);

  /* Determine what fraction of the separation distance the shift is */
  if (dist > 0)
    shift_fraction = SHIFT_SIZE / dist;
  else
    shift_fraction = 1.0;

  /* Calculate shift to be applied */
  for (i = 0; i < 2; i++)
    {
      shift[i] = shift_fraction * (cofm[i] - ligand_cofm[i]);
      if (fabs(shift[i]) < MIN_SHIFT)
        shift[i] = MIN_SHIFT * sign(cofm[i] - ligand_cofm[i]);
    }

  /* Get the first of this residue's atoms */
  atom_ptr = residue_ptr->first_atom_ptr;
  natoms = residue_ptr->natoms;
  iatom = 0;

  /* Loop over the residue's atoms to apply shift */
  while (atom_ptr != NULL && iatom < natoms)
    {
      /* Apply the shift to the coordinates */
      for (i = 0; i < 2; i++)
        atom_ptr->coords[i] = atom_ptr->coords[i] + shift[i];

      /* Get the next atom record and increment atom count */
      atom_ptr = atom_ptr->next_atom_ptr;
      iatom++;
    }

  /* Flag residue as shifted */
  residue_ptr->shifted = TRUE;
}
/***********************************************************************

shift_sidechains  -  Shift any old side chains out of the way of the
                     ligand

***********************************************************************/

int shift_sidechains(struct pdbentry *first_pdbentry_ptr,
                     struct pdbentry *ligand_pdbentry_ptr,int nligres)
{
  int noverlap, tot_overlap;
  int done, shifted, shift_made;

  double cofm[3], ligand_cofm[3];

  struct pdbentry *pdbentry_ptr;
  struct residue *residue_ptr, *ligand_residue_ptr;
  struct residue *other_residue_ptr;

  /* Initialise variables */
  shift_made = FALSE;

  /* Update the maximum and minimum coordinates of the ligand */
  update_all_maxmin(ligand_pdbentry_ptr);

  /* Calculate the C of M of the ligand */
  calc_cofm(ligand_pdbentry_ptr->first_residue_ptr,nligres,ligand_cofm);

  /* Get the first of the ensemble entries */
  pdbentry_ptr = first_pdbentry_ptr->next_pdbentry_ptr;

  /* Loop through the PDB entries */
  while (pdbentry_ptr != NULL)
    {
      /* Get this entry's first residue */
      residue_ptr = pdbentry_ptr->first_residue_ptr;

      /* Loop over the non-ligand residues to see if any overlap our
         ligand, shifting them further and further away if they do */
      while (residue_ptr != NULL)
        {
          /* If not ligand, then check overlap */
          if (residue_ptr->interaction != IS_LIGAND)
            {
              /* Initialise flag */
              done = FALSE;
              shifted = FALSE;

              /* Loop while this residue overlaps the ligand */
              while (done == FALSE)
                {
                  /* Get the first ligand residue */
                  ligand_residue_ptr
                    = ligand_pdbentry_ptr->first_residue_ptr;

                  /* Initialise the overlap score */
                  tot_overlap = 0;

                  /* Loop over the ligand residues to check for overlap */
                  while (ligand_residue_ptr != NULL)
                    {
                      /* See if the residue overlaps the ligand */
                      noverlap
                        = check_overlap(residue_ptr,ligand_residue_ptr);

                      /* Add to total overlap score */
                      tot_overlap = tot_overlap + noverlap;

                      /* Get pointer to next residue record */
                      ligand_residue_ptr
                        = ligand_residue_ptr->next_residue_ptr;
                    }

                  /* If side chain doesn't overlap the ligand, but has
                     been subject to shifts, then check that it hasn't
                     bumped into any other side chains on its own plot */
                  if (tot_overlap == 0 && shifted == TRUE)
                    {
                      /* Get the first of its plot's residues */
                      other_residue_ptr
                        = pdbentry_ptr->first_residue_ptr;

                      /* Loop over the other side chains to check for
                         bumps */
                      while (other_residue_ptr != NULL)
                        {
                          /* If not ligand, and not our residue, then
                             check overlap */
                          if (other_residue_ptr != residue_ptr &&
                              other_residue_ptr->interaction != IS_LIGAND &&
                              tot_overlap == 0)
                            {
                              /* See if the residues overlap */
                              noverlap
                                = check_overlap(residue_ptr,
                                                other_residue_ptr);

                              /* Add to total overlap score */
                              tot_overlap = tot_overlap + noverlap;
                            }

                          /* Get pointer to next residue record */
                          other_residue_ptr
                            = other_residue_ptr->next_residue_ptr;
                        }
                    }

                  /* If have an overlap, need to shift the current side chain
                     away from the ligand */
                  if (tot_overlap > 0)
                    {
                      /* Calculate this residue's C of M */
                      calc_cofm(residue_ptr,1,cofm);

                      /* Shift the residue */
                      make_shift(residue_ptr,cofm,ligand_cofm);

                      /* Update this residue's maximum and minimum
                         coords */
                      update_max_min(residue_ptr);

                      /* Set flags that shift made */
                      shifted = shift_made = TRUE;
                    }
                  else
                    done = TRUE;
                }
            }

          /* Get pointer to next residue record */
          residue_ptr = residue_ptr->next_residue_ptr;
        }

      /* Get the next entry */
      pdbentry_ptr = pdbentry_ptr->next_pdbentry_ptr;
    }

  /* Return whether any shift was made */
  return(shift_made);
}
/***********************************************************************

map_ligands  -  Map the ligand residues just read in onto original
                ligand residues so that LIGPLOT coords just generated
                will be used for the final plot

***********************************************************************/

void map_ligands(struct residue *first_ligand_residue_ptr,
                 struct pdbentry *ligand_pdbentry_ptr)
{
/* v.1.3 --> */
  int iatom, jatom, natoms1, natoms2;
  int done;

  struct atom *atom1_ptr, *atom2_ptr;
/* <-- v.1.3 */
  struct residue *residue1_ptr, *residue2_ptr;

  /* Get the first ligand residue */
  residue1_ptr = first_ligand_residue_ptr;

  /* Loop over the ligand residues */
  while (residue1_ptr != NULL && residue1_ptr->interaction == IS_LIGAND)
    {
      /* Get the first ligand residue */
      residue2_ptr = ligand_pdbentry_ptr->first_residue_ptr;

      /* Loop over the ligand residues to check for overlap */
      while (residue2_ptr != NULL)
        {
          /* Mark residue as a ligand residue */
          residue2_ptr->interaction = IS_LIGAND;

          /* If this is the right residue, create mapping */
          if (!strcmp(residue1_ptr->res_num,residue2_ptr->res_num) &&
              !strcmp(residue1_ptr->res_name,residue2_ptr->res_name) &&
              residue1_ptr->chain == residue2_ptr->chain)
/* v.1.3 --> */
            {
/* <-- v.1.3 */
              residue1_ptr->map_residue_ptr = residue2_ptr;

/* v.1.3 --> */
              /* Get the first of this residue's atoms */
              atom1_ptr = residue1_ptr->first_atom_ptr;
              natoms1 = residue1_ptr->natoms;
              iatom = 0;

              /* Loop over the residue's atoms to transfer the weights */
              while (atom1_ptr != NULL && iatom < natoms1)
                {
                  /* Get the first of the other residue's atoms */
                  atom2_ptr = residue2_ptr->first_atom_ptr;
                  natoms2 = residue2_ptr->natoms;
                  jatom = 0;
                  done = FALSE;

                  /* Loop over the residue's atoms to transfer the weights */
                  while (atom2_ptr != NULL && jatom < natoms2 &&
                         done == FALSE)
                    {
                      /* If atoms match. transfer the weight */
                      if (!strcmp(atom1_ptr->atom_name,atom2_ptr->atom_name))
                        {
                          /* Copy across the weight */
                          atom2_ptr->anchor_weight
                            = atom1_ptr->anchor_weight;

                          /* Set flag */
                          done = TRUE;
                        }

                      /* Get the next atom record and increment atom count */
                      atom2_ptr = atom2_ptr->next_atom_ptr;
                      jatom++;
                    }

                  /* Get the next atom record and increment atom count */
                  atom1_ptr = atom1_ptr->next_atom_ptr;
                  iatom++;
                }
            }
/* <-- v.1.3 */

          /* Get pointer to next residue record */
          residue2_ptr = residue2_ptr->next_residue_ptr;
        }

      /* Get the next residue */
      residue1_ptr = residue1_ptr->next_residue_ptr;
    }
}
/***********************************************************************

write_shifts_pdb  -  Write out all the atom shifts to the shifts.pdb
                     file for reading in by LigPlus

***********************************************************************/

void write_shifts_pdb(struct pdbentry *first_pdbentry_ptr,
                      struct parameters params)
{
  char file_name[FILENAME_LEN];

  int atom_number;
  int shift_made;

  struct pdbentry *pdbentry_ptr;
  struct residue *residue_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  atom_number = 0;
  shift_made = FALSE;

  /* Get the first of the ensemble entries */
  pdbentry_ptr = first_pdbentry_ptr->next_pdbentry_ptr;

  /* Loop through the PDB entries */
  while (pdbentry_ptr != NULL)
    {
      /* Form name of this output file */
      sprintf(file_name,"shifts%d.pdb",pdbentry_ptr->pdb_id);

      /* Open the output file */
      file_ptr = open_file(file_name,params);

      /* Get this entry's first residue */
      residue_ptr = pdbentry_ptr->first_residue_ptr;

      /* Loop over the residues to write out */
      while (residue_ptr != NULL)
        {
          /* If residue was shifted, then write all its atoms out */
          if (residue_ptr->shifted == TRUE)
            atom_number
              = write_residue(file_ptr,residue_ptr,atom_number);

          /* Get pointer to next residue record */
          residue_ptr = residue_ptr->next_residue_ptr;
        }

      /* Close the current file */
      fclose(file_ptr);

      /* Get the next PDB entry */
      pdbentry_ptr = pdbentry_ptr->next_pdbentry_ptr;
    }
}

/***********************************************************************

                             M   A   I   N

***********************************************************************/

int main(int argc,char *argv[])
{
  char hhb_name[FILENAME_LEN], nnb_name[FILENAME_LEN];
/* v.1.3 --> */
//  char pdb_name[FILENAME_LEN];
  char ligall_pdb[FILENAME_LEN], ligplus_pdb[FILENAME_LEN];
/* <-- v.1.3 */

  int fit_type, loop, nfit, nhbonds, nlig, nligres, nnbonds, npdb, nresid;
/* v.2.0 --> */
  int domain1, domain2, loop_start, ndummy, nlayout;
/* <-- v.2.0 */
/* v.1.3 --> */
//  int have_cora, overlay, shifted;
  int have_cora, ok, overlay, shifted;
/* <-- v.1.3 */

/* v.2.0 --> */
  float ave_y[2];
/* <-- v.2.0 */

  double cofg[2][3];

  REAL matrix[3][3];

  struct parameters params;

  struct atomap *first_atomap_ptr;
/* v.2.0 --> */
  struct domain *first_domain_ptr;
  struct layout *first_layout_ptr, *last_layout_ptr;
/* <-- v.2.0 */
  struct pdbentry *first_pdbentry_ptr, *fit_pdbentry_ptr;
/* v.1.3 --> */
  struct pdbentry *full_pdbentry_ptr, *ref_pdbentry_ptr;
/* <-- v.1.3 */
  struct pdbentry *ligand_pdbentry_ptr;
  struct residue *first_ligand_residue_ptr;

/* v.1.3 --> */
  struct resseq **alignment;
/* <-- v.1.3 */

  /* Initialise variables */
  printf("Running ligfit ...\n");
  fflush(stdout);
/* v.2.0 --> */
  domain1 = 1;
  domain2 = 2;
  first_domain_ptr = NULL;
/* <-- v.2.0 */
  first_atomap_ptr = NULL;
/* v.1.3 --> */
  fit_type = NO_FIT;
/* <-- v.1.3 */
  have_cora = FALSE;
  ligand_pdbentry_ptr = NULL;
  nfit = 0;
  shifted = FALSE;

/* v.2.2.1 --> */
  /* Initialise pointers */
  first_ligand_residue_ptr = NULL;
/* <-- v.2.2.1 */

  /* Interpret the command-line parameters */
  get_command_arguments(argv,argc - 1,&params);

/* v.2.0 --> */
  /* For domain-based DIMPLOT, read in domain definitions */
  if (params.dimplot == TRUE && params.have_domains == TRUE)
    get_domains("ligplus.dom",domain1,domain2,&first_domain_ptr,
                params);
/* <-- v.2.0 */

  /* Read in the PDB ensemble from the ensemble.drw file */
  npdb = read_drw_file(&first_pdbentry_ptr,"ensemble.drw",&nresid,
                       params);
  if (npdb == 0)
    exit(1);

  /* Sort the residues in each PDB entry in the ensemble */
  sort_residues(first_pdbentry_ptr);

  /* Form the names of the input PDB and HBPLUS files */
/* v.1.3 --> */
  init_filenames(hhb_name,nnb_name,ligall_pdb,ligplus_pdb,params);
//  init_filenames(hhb_name,nnb_name,pdb_name,params);
/* <-- v.1.3 */

/* v.1.3 --> */
  /* Update all the residue max and min values */
  update_all_maxmin(first_pdbentry_ptr);

  /* Classify all the residues according to whether they interact directly
     with the ligand */
/* v.2.0 --> */
//  check_interactions(first_pdbentry_ptr);
  check_interactions(first_pdbentry_ptr,params.dimplot);
/* <-- v.2.0 */

  /* Classify all the residues by the type of interaction they make */
/* v.2.0 --> */
  if (params.dimplot == FALSE)
/* <-- v.2.0 */
    classify_residues(first_pdbentry_ptr);

  /* Read in the full coordinates of the reference structure */
/* v.2.2.1 --> */
//  full_pdbentry_ptr = read_pdb_file(ligall_pdb,params);
  full_pdbentry_ptr
    = read_pdb_file(ligall_pdb,LIGALL_FILE,&first_ligand_residue_ptr,
                    params);
/* <-- v.2.2.1 */
  if (full_pdbentry_ptr == NULL)
    full_pdbentry_ptr = first_pdbentry_ptr;

  /* Replace the residues by those from the plot */
  else
    replace_plot_residues(full_pdbentry_ptr,first_pdbentry_ptr);

/* v.2.0 --> */
  /* For DIMPLOT , assign the domains read in from the .dom file */
  if (params.dimplot == TRUE)
    assign_domains(full_pdbentry_ptr,first_domain_ptr);
/* <-- v.2.0 */

  /* Replace the PDB entry from the .drw file with the one from the
     full PDB file */
  full_pdbentry_ptr->next_pdbentry_ptr
    = first_pdbentry_ptr->next_pdbentry_ptr;
  ref_pdbentry_ptr = first_pdbentry_ptr = full_pdbentry_ptr;
/* <-- v.1.3 */

  /* Read in the coordinates of the ligand and any nearby residues
     from ligplus.pdb */
/* v.1.3 --> */
//  fit_pdbentry_ptr = read_pdb_file(pdb_name,params);
/* v.2.2.1 --> */
//  fit_pdbentry_ptr = read_pdb_file(ligplus_pdb,params);
  fit_pdbentry_ptr
    = read_pdb_file(ligplus_pdb,LIGPLUS_FILE,&first_ligand_residue_ptr,
                    params);
/* <-- v.2.2.1 */
/* <-- v.1.3 */
  if (fit_pdbentry_ptr == NULL)
    exit(1);

  /* Put the target PDB entry at front of linked list */
  fit_pdbentry_ptr->next_pdbentry_ptr = first_pdbentry_ptr;
  first_pdbentry_ptr = fit_pdbentry_ptr;

/* v.2.0 --> */
  /* Sort the residues in each PDB entry in the ensemble */
  sort_residues(first_pdbentry_ptr);
/* <-- v.2.0 */

  /* Update all the residue max and min values */
  update_all_maxmin(first_pdbentry_ptr);

  /* Calculate the covalent bonds */
  calc_bonds(fit_pdbentry_ptr);

  /* Read in all the H-bonds involving the ligand of interest */
  nhbonds = read_hbplus_file(hhb_name,HBOND,fit_pdbentry_ptr,params);

  /* Repeat for the non-bonded contacts */
  nnbonds = read_hbplus_file(nnb_name,CONTACT,fit_pdbentry_ptr,params);

  /* Classify all the residues according to whether they interact directly
     with the ligand */
/* v.2.0 --> */
//  check_interactions(first_pdbentry_ptr);
  check_interactions(first_pdbentry_ptr,params.dimplot);
/* <-- v.2.0 */

  /* Classify all the residues by the type of interaction they make */
/* v.2.0 --> */
  if (params.dimplot == FALSE)
/* <-- v.2.0 */
    classify_residues(first_pdbentry_ptr);

/* v.1.3 --> */
  /* If not fitting on ligand, perform fit based on proteins */
  if (params.fit_ligands != TRUE)
    {
/* <-- v.1.3 */
/* v.1.2 --> */
      /* If have a structural alignment, read it in from the CORA file */
/* v.1.3 --> */
//  if (params.alnfile != &null)
      if (fit_type == NO_FIT && params.alnfile != &null)
/* <-- v.1.3 */
        {
          /* Get the alignment in CORA format */
/* v.1.3 --> */
          if (strstr(params.alnfile,".cora") != NULL)
/* <-- v.1.3 */
            nfit = get_cora_alignment(params.alnfile,first_pdbentry_ptr,
/* v.2.0.1 --> */
//                                      &alignment);
                                      &alignment,params);
/* <-- v.2.0.1 */
/* v.1.3 --> */
          /* Get the alignment in CORA format */
          else if (strstr(params.alnfile,".caf") != NULL)
            nfit = get_caf_alignment(params.alnfile,first_pdbentry_ptr,
/* v.2.0.1 --> */
//                                     &alignment);
                                     &alignment,params);
/* <-- v.2.0.1 */

          /* Get the alignment in FASTA format */
          else if (strstr(params.alnfile,".fasta") != NULL)
            nfit = get_fasta_alignment(params.alnfile,first_pdbentry_ptr,
/* v.2.0.1 --> */
//                                       &alignment);
                                       &alignment,params);
/* <-- v.2.0.1 */
/* <-- v.1.3 */

          /* If have a minimum alignment, then have fit */
          if (nfit > MIN_ALIGNED)
            fit_type = FIT_CORA;
        }
/* v.1.3 --> */
//  else
//    {
///* <-- v.1.2 */
//      /* Delete any residues that do not make any interactions */
//      nresid = delete_unwanted_residues(fit_pdbentry_ptr);
//      if (nresid == 0)
//        {
//          printf("*** ERROR. No ligand-interacting residues found!\n");
//          exit(1);
//        }
///* v.1.2 --> */
//    }
///* <-- v.1.2 */

      /* If don't have a fit, try fitting based on whole-sequence
         alignment */
      if (fit_type == NO_FIT)
        {
          /* Get the sequence alignment */
/* v.2.0 --> */
//          nfit = get_seq_alignment(first_pdbentry_ptr,&alignment);
          nfit = get_seq_alignment(first_pdbentry_ptr,&alignment,params);
/* <-- v.2.0 */

          /* If have a minimum alignment, then have fit */
          if (nfit > MIN_ALIGNED)
            fit_type = FIT_ALIGNMENT;
        }

      /* If still don't have a fit, try sequence alignment based on binding
         site residues only */
      if (fit_type == NO_FIT)
        {
          /* Get the binding-site alignment */
/* v.2.0 --> */
//          nfit = get_bs_alignment(first_pdbentry_ptr,&alignment);
          nfit = get_bs_alignment(first_pdbentry_ptr,&alignment,params);
/* <-- v.2.0 */

          /* If have a minimum alignment, then have fit */
          if (nfit > MIN_ALIGNED)
            fit_type = FIT_ALIGNMENT;
        }

      /* If have a set of residue equivalences, use them to calculate
         transformation of current structure onto reference structure */
      if (fit_type != NO_FIT)
/* v.2.0.1 --> */
//        nfit = fit_structures(alignment,nfit,cofg,matrix);
        nfit = fit_structures(alignment,nfit,cofg,matrix,params);
/* <-- v.2.0.1 */

      /* If don't have enough fitted residues, mark as unfitted */
      if (nfit <= MIN_ALIGNED)
        fit_type = NO_FIT;
    }

  /* If fitting on ligands, or fit on proteins failed, map current ligand
     to the ligand in the reference structure */
/* v.2.0 --> */
//  if (params.fit_ligands == TRUE || fit_type == NO_FIT)
  if (params.dimplot == FALSE &&
      (params.fit_ligands == TRUE || fit_type == NO_FIT))
/* <-- v.2.0 */
    {
      /* Fit the ligands by called hbadd */
      ok = match_ligands(first_pdbentry_ptr,matrix,cofg,
                         &first_atomap_ptr,params);

      /* If fit was successful, set fit type */
      if (ok == TRUE)
        fit_type = FIT_LIGANDS;
      else
        {
          first_atomap_ptr = NULL;
          fit_type = NO_FIT;
        }
    }

  /* If have a set of residue equivalences, use them to calculate
     transformation of current structure onto reference structure */
  if (fit_type != NO_FIT)
    {
      /* See if any residues outside the alignment, plus ligands, can be
         mapped onto one another */
      nfit = map_extras(first_pdbentry_ptr->first_residue_ptr,
                        ref_pdbentry_ptr->first_residue_ptr,
/* v.1.3.7 --> */
//                        matrix,cofg,nfit,&first_atomap_ptr,fit_type);
                        matrix,cofg,nfit,&first_atomap_ptr,fit_type,params);
/* <-- v.1.3.7 */

      /* Identify equivalent interactions to give more weight to the
         corresponding atom anchors */
      weigh_anchors(first_atomap_ptr);

      /* Write the transformed coords to the ligplus.pdb file, overwriting
         the current version */
      write_coords(first_pdbentry_ptr,matrix,cofg,params);
    }
/* <-- v.1.3 */

  /* Align sequence of interacting residues to sequences from LIGPLOTs
     in the current ensemble */
//  if (nfit <= MIN_ALIGNED)
//    fit_type
//      = fit_to_ligplot(first_pdbentry_ptr,&first_atomap_ptr,params);
//  else
//    fit_type = FIT_CORA;


/* v.2.0 --> */
  if (params.dimplot == TRUE)
    {
      /* Read in the layout as computed by dimhtml */
      nlayout = read_rcm_file(fit_pdbentry_ptr->first_residue_ptr,
                              &first_layout_ptr,&last_layout_ptr,
                              params);

      /* Add dummy layouts where have empty slots in the primary chain */
      ndummy = add_dummy_layouts(&first_layout_ptr,last_layout_ptr,
                                 &(ref_pdbentry_ptr->first_residue_ptr),
                                 nlayout,ave_y);
      ref_pdbentry_ptr->nresid = ref_pdbentry_ptr->nresid + ndummy;

      /* Determine how to arrange the residues that haven't been mapped */
      shifted
        = arrange_dimplot(ref_pdbentry_ptr,fit_pdbentry_ptr,
                          first_layout_ptr,ave_y,params.order);

      /* Write out the dimplot restraints */
      dimplot_restraints(ref_pdbentry_ptr,fit_pdbentry_ptr,fit_type,
                         params);

      /* If any shift was made, need to write out all shifts
         to the shifts.pdb file to adjust existing LIGPLOTs */
      if (shifted == TRUE)
        write_shifts_pdb(first_pdbentry_ptr,params);

      /* Finish here */
      return(0);
    }
/* <-- v.2.0 */

  /* Determine number of fitting loops - two for LIGPLOT, and just one for
     DIMPLOT */
  loop_start = 0;
  if (params.dimplot == TRUE)
    loop_start = 1;

  /* Loop twice, first to write out ligand equivalences to ligand.rcm for
     fitting ligands (if required), and then write out final equivalences
     to ligplus.rcm */
/* v.2.0 --> */
//  for (loop = 0; loop < 2; loop++)
  for (loop = loop_start; loop < 2; loop++)
/* <-- v.2.0 */
    {
      /* Use the fitting to determine the atom and residue equivalences */
      overlay
/* v.1.3.7 --> */

        //  = get_equivalences(loop,fit_pdbentry_ptr,&first_ligand_residue_ptr,
        //                   fit_type,first_atomap_ptr,params);
        = get_equivalences(loop,ref_pdbentry_ptr,fit_pdbentry_ptr,
                           &first_ligand_residue_ptr,fit_type,
                           first_atomap_ptr,params);
/* <-- v.1.3.7 */

      /* If ligands to be overlaid, perform the overlay */
      if (loop == 0 && overlay == TRUE)
        {
          /* Write out the ligand.pdb, .hhb and .nnb files */
/* v.2.0.3 --> */
//          write_ligand_files(first_ligand_residue_ptr,params);
          write_ligand_files(first_ligand_residue_ptr,
                             fit_pdbentry_ptr->first_bond_ptr,params);
/* <-- v.2.0.3 */

          /* Run LIGPLOT */
          run_ligplot(params);

          /* Read in the resultant plot to get the flattened coords of
             the current ligand */
          nlig = read_drw_file(&ligand_pdbentry_ptr,"ligplot.drw",
                               &nligres,params);

          /* If read in, then see if any sidechains from the previous plots
             need to be shifted out of the way of ligand */
          if (nlig > 0 && ligand_pdbentry_ptr != NULL)
            {
/* v.1.2 --> */
              /* Fit the 2D coords of the ligand onto previous ligand
                 using the atom equivalences previously computed */
              fit_2d_ligands(first_atomap_ptr,ligand_pdbentry_ptr,params);
/* <-- v.1.2 */

              /* Shift any old side chains out of the way of the ligand */
              shifted
                = shift_sidechains(first_pdbentry_ptr,ligand_pdbentry_ptr,
                                   nligres);

              /* Map the ligand residues just read in onto original
                 ligand residues so that LIGPLOT coords just generated
                 will be used for the final plot */
              map_ligands(first_ligand_residue_ptr,ligand_pdbentry_ptr);

              /* If any shift was made, need to write out all shifts
                 to the shifts.pdb file to adjust existing LIGPLOTs */
              if (shifted == TRUE)
                write_shifts_pdb(first_pdbentry_ptr,params);
            }
        }

/* v.1.3 --> */
      /* Reinitialise first atom map pointer */
//      first_atomap_ptr = NULL;
/* <-- v.1.3 */
    }

  return(0);
}

