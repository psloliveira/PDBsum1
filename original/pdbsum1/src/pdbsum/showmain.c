/***********************************************************************

  showmain.c - Program to display new-style HTML pages for given PDBsum
               entry

************************************************************************

Date:-         8 Dec 2003.

Written by:-   Roman Laskowski


Testing statement in thumbnail_3Dmoljs.html:

    ~/src/pdbsum/showmain "1aaw" main.html -c 0 -f ~/public_html/showmain.html
    ~/src/pdbsum/showmain "1e11" main.html -c 0 -f ~/public_html/showmain.html
    ~/src/pdbsum/showmain "4wt1" main.html -c 0 -f ~/public_html/showmain.html


----------------------------------------------------------------------*/

/* Datafiles
   ---------

Actual name       Description
-----------       -----------
database.dat      Input file containing all relevant data for given
                  PDB file
interfaces.dat    As database.dat, but containing information about
                  any protein-protein interfaces
ligands.dat       As database.dat, but containing information about
                  any ligands
seqinfo.dat       As database.dat, but containing information from
                  MSD database
template_XXX.html Input HTML file acting as template for required
                  HTML page
showmain.html     Output HTML file
variants.dat      Input file containing any human variation data

*/

/* Compilation
   -----------

cc -c showmain.c
cc -c common.c
cc -c reapar.c
cc -o showmain showmain.o common.o reapar.o -lm -lz

Optimization flags:-

cc -O2 -funroll-loops -c showmain.c
cc -O2 -funroll-loops -c common.c
cc -O2 -funroll-loops -c reapar.c
cc -o showmain showmain.o common.o reapar.o -lm -lz


*/

/* Function calling-tree
   ---------------------

main
   -> create_field
   -> get_command_arguments
   -> set_loop_level
   -> get_all_paths
         -> get_next_param (in reapar.c)
         -> create_field_record
               -> get_subscripts
               -> store_string
               -> add_field <----------+
                     -> create_field   |
                     -> add_field -----+
   -> create_pdir
   -> get_pdb_filename
         -> find_pdb_file
         -> create_field_record
   -> create_vdir
   -> read_database_file
         -> get_field
               -> find_field <----------+
                     -> find_field -----+
         -> string_chop (in common.c)
         -> perform_splice (in common.c)
         -> create_field_record (as above)
   -> read_dbt_file
         -> string_chop (in common.c)
         -> get_tokens (in common.c)
         -> get_field (as above)
         -> store_motif_count
               -> create_field_record (as above)
               -> store_string
   -> define_specials
         -> create_field_record (as above)
         -> create_contents_flags
               -> create_field_record (as above)
         -> create_tab_flags
               -> create_field_record (as above)
         -> create_loop_flags
               -> create_field_record (as above)
         -> create_structure_flags
               -> create_field_record (as above)
         -> create_call_flags
               -> create_field_record (as above)
         -> create_colour_flags
               -> create_field_record (as above)
         -> create_counters
               -> create_field_record (as above)
   -> get_html_template
         -> get_field (as above)
         -> string_truncate
         -> annotate_line (see below)
         -> process_blocks
               -> store_string
         -> store_string
   -> write_html_file
         -> get_directories
         -> open_output_file (in common.c)
         -> get_image_names
         -> check_field
               -> check_for_brackets
                     -> perform_splice (in common.c)
               -> check_for_string
               -> get_field (as above)
         -> check_for_file
         -> check_in_file
         -> end_loop
               -> get_field (as above)
         -> process_gif_names
               -> perform_splice (in common.c)
               -> extract_image_name
                     -> annotate_line
                     -> find_image_name
                     -> create_iname_record
                           -> store_string (in common.c)
         -> process_hrefs
               -> perform_splice (in common.c)
         -> process_ebi_incs
               -> perform_splice (in common.c)
         -> insert_file
               -> perform_splice (in common.c)
         -> annotate_line
               -> apply_template
                     -> get_field
                     -> form_filename
               -> extract_tag
                     -> check_for_brackets
                           -> perform_splice (in common.c)
               -> get_field (as above)
               -> get_subvalue
               -> link_to_pdb
                     -> get_tokens (in common.c)
                     -> perform_splice (in common.c)
               -> link_to_ec
                     -> get_tokens (in common.c)
                     -> get_field (as above)
                     -> perform_splice (in common.c)
               -> punc_chain
               -> type_set (in common.c)
               -> fix_upper (in common.c)
               -> to_lower
               -> html_format_formula
               -> split_text
               -> store_counter
                     -> get_field (as above)
                     -> store_string
               -> store_text
                     -> get_field (as above)
                     -> store_string
               -> perform_splice (in common.c)
         -> res_conservation
         -> read_offset_index
               -> get_field (as above)
               -> create_field_record (as above)
               -> store_string
         -> read_database_file (see above)
         -> fix_plurals
               -> perform_splice (in common.c)
   -> write_image_names

*/

#include "common.h"
#include "showmain.h"

/* Prototypes
   ---------- */

/* In common.c */
void c_get_paths(struct c_file_data *file_name_ptr,int run_64);
void fix_upper(char *name);
char *form_filename(char *pdb_code,char *filename_template);
int get_tokens(char *line,char **token,int max_token,char token_sep);
FILE *open_output_file(char *out_name,int terminate);
int perform_splice(char *line,char *string,char *new_string);
void store_string(char **string_ptr,char *string);
int string_chop(char *string);
int string_chomp(char *string);
int string_truncate(char *string,int max_length);
void to_lower(char *search_string,int length);
void to_upper(char *search_string,int length);
void type_set(char *name,int cap_start);

/* In reapar.c */
int get_next_param(char *token_name,char *file_name,int filename_len,
                   int exit_on_error);
int get_param(char *token_name,char *file_name,int filename_len,
              int exit_on_error);

/***********************************************************************

get_command_arguments  -  Interpret the command line arguments

***********************************************************************/

void get_command_arguments(char *strptrs[],int num,
                           struct parameters *params)
{
  char ch, number_string[20], tmp_string[LINELEN], tmp[LINELEN];

  char *string_ptr;

  int i, itoken, len, level;
  int done;

  long value;

  /* Initialise parameters */
  params->version = 2;
  params->pdb_code[0] = '\0';
  params->chain = '\0';
  params->ec_code[0] = '\0';
  params->template_name[0] = '\0';
  params->template_dir[0] = '\0';
  params->template_real_dir[0] = '\0';
  params->database_dat[0] = '\0';
  for (i = 0; i < MAX_LEVEL; i++)
    {
      params->ec_number[i] = 0;
      params->instance[i] = 0;
    }
  for (i = 0; i < NPAR; i++)
    params->parameter[i][0] = '\0';
  params->output_file[0] = '\0';
  strcpy(params->option,"DEFAULT_OPTION");
  params->include_file[0] = '\0';
  params->fullname = FALSE;
  params->call_no = 0;
  params->start_no = 0;
  params->have_error = FALSE;
  params->nlevels = 0;
  params->enzymes = FALSE;
  params->drat = FALSE;
  params->drugport = FALSE;
  params->profunc = FALSE;
  params->ligsearch = FALSE;
  params->ligplus = FALSE;
  params->varsite = FALSE;
  params->sas = FALSE;
  params->highlights = FALSE;
  params->bottom_level = FALSE;
  params->pdb_type = STANDARD_PDB;
  params->relinks = FALSE;
  params->fake_entry = FALSE;
  params->user_id = &null;
  params->link_path = &null;
  params->score[0] = '\0';
  params->security_code[0] = '\0';
  params->dir_name[0] = '\0';
  params->subdir[0] = '\0';
  params->subdir_real[0] = '\0';
  params->subdir_dir[0] = '\0';
  params->work_dir[0] = '\0';
  params->from_ebi = FALSE;
  params->drug_id[0] = '\0';
  params->uniprot_acc[0] = '\0';
  params->ligand_ranges = &null;
  params->list_fields = FALSE;
  params->proprietary = FALSE;

  /* If nothing entered on the command line, show options available */
  if (num < 2 || (num == 1 && (!strcmp(strptrs[1],"-help") ||
                               !strcmp(strptrs[1],"-h"))))
    {
      printf("\n");
      printf("Correct usage is ...\n"
             "\n"
             "  showmain  pdb_code template  [-l ligand/chain/E.C. no]  "
             "[-f output_file]\n"
             "            [-option option]  [{-i | -I} include_file]  "
             "[-c call-no]\n"
             "            [-e]  [-enzymes]  [-profunc]  [-highlights]  "
             "[-u user-id]\n"
             "            [-s start-no]  [-code security-code]  "
             "[-dir dirname]\n"
             "            [UPLOAD | MCSG]  [-wkdir working_dir]  "
             "[-score score]\n"
             "            [-EBI]  [-par1 parameter1] [-par2 parameter2] "
             " ... \n"
             "            [-par10 parameter10] [-drat] [-drugport] [-sas] "
             "[-chain chain_id]\n"
             "            [-drug_id DBref] [-uniprot_acc UniProt]\n"
             "            [-db database_dat]  [-tdir template_dir]  "
             "[-trdir template_real_dir]\n"
             "            [-ligsearch]  [-ligand_ranges ranges]  [-prop]  "
             "[-ligplus]  [-disastr]  [-varsite]\n"
             "            [-pdbsum1]  [-relinks]  [-link_path path]\n "
             "where\n"
             "   * pdb_code        = PDB entry to be displayed\n"
             "   * template        = Name of template page to be "
             "displayed\n"
             "   * -l chain/ligand/E.C. no  = chain/ligand/E.C. number to "
             "be shown\n"
             "   * -f output_file  = write output to given file rather "
             "than to standard\n"
             "                       output\n"
             "   * -option option  = switch on appropriate OPTION_XXXX "
             "block in the template\n"
             "                       file\n"
             "   * -i include_file = file to be included at the #include "
             "line\n"
             "   * -I include_file = full path and filename of file to "
             "be included at\n"
             "                       the #include line\n"
             "   * -c call number, indicating which menu item has been "
             "called and should be\n"
             "                       highlighted\n"
             "   * -s start number, indicating which start number for "
             "given group of items\n"
             "   * -e              = have error in PDB code\n"
             "   * -enzymes        = pages are for Enzyme Structures "
             "Database\n"
             "   * -parN parameter N = additional parameters, numbered "
             "1-10\n"
             "   * -profunc        = pages are for ProFunc\n"
             "   * -highlights     = PDBsum highlights page\n"
             "   * -u user-id      = Identifier for temporary data dir\n"
             "   * -code security-code = Password for secure entries\n"
             "   * -dir dirname    = Name of temporary working directory\n"
             "   * UPLOAD          = Uploaded PDB file for PDBsum "
             "generation\n"
             "   * MCSG            = MCSG structure\n"
             "   * -wkdir working_dir = Name of working directory\n"
             "   * -score score    = Score required to be passed to some "
             "pages\n"
             "   * -EBI            = Flag indicating that called from EBI "
             "PDBsum page\n"
             "   * -drat           = Flag indicating that called from DraT DB\n"
             "   * -drugport       = Flag indicating that called from DrugPort\n"
             "   * -sas            = Flag indicating that called from SAS\n"
             "   * -chain chain_id = PDB chain identifier\n"
             "   * -drug_id        = DrugBank id for DrugPort drug page\n"
             "   * -uniprot_acc    = UniProt accession for DrugPort target "
             "page\n"
             "   * -db database_dat = Name and path of database.dat file to "
             "be read in at start\n"
             "   * -trdir template_dir = Real directory where HTML templates "
             "are to be found\n"
             "   * -tdir template_dir = Equivalent logical name\n"
             "   * -ligsearch      = pages are for LigSearch\n"
             "   * -ligand_ranges ranges = Ligands selected on ligand "
             "clusters page\n"
             "   * -prop           = proprietary installation of PDBsum\n"
             "   * -ligplus        = pages are for LigPlus licence system\n"
             "   * -disastr        = pages are for DisaStr\n"
             "   * -varsite        = pages are for VarSite\n"
             "   * -list           = list fields (for debugging)\n"
             "   * -pdbsum1        = Standalone version of PDBsum\n"
             "   * -relinks        = Make all image links relative\n"
             "   * -link_path      = relative path for common images, etc\n"
             "\n"
             "\n");
      exit(1);
    }

  /* Get the PDB code from the first parameter */
  strcpy(params->pdb_code,strptrs[1]);

  /* Check whether this is an Alpha Fold structure (first character is
     an upper-case letter) */
  ch = params->pdb_code[0];
  if (ch >= 'A' && ch <= 'Z')
    params->pdb_type = ALPHA_FOLD;

  /* Get the template name from the second parameter */
  strcpy(params->template_name,strptrs[2]);

  /* Loop through the remaining command arguments */
  for (itoken = 3; itoken < num + 1; itoken++)
    {
      /* Get the length of this token */
      len = strlen(strptrs[itoken]);

      /* Check for the -l option indicating required chain/ligand/E.C.
         number */
      if (!strcmp(strptrs[itoken],"-l"))
        {
          /* Get the ligand number from the next token, if there is one */
          if (itoken < num)
            {
              /* Prepare to unpack the level values (if more than one) */
              strcpy(tmp_string,strptrs[itoken + 1]);
              strcpy(params->ec_code,tmp_string);
              done = FALSE;
              level = 0;

              /* Loop until all numbers have been unpacked */
              while (done == FALSE)
                {
                  /* Get position of dot, if there is one */
                  string_ptr = strstr(tmp_string,".");
                  if (string_ptr != NULL)
                    {
                      /* Get length of number and extract */
                      len = string_ptr - tmp_string;
                      if (len > 0)
                        {
                          strncpy(number_string,tmp_string,len);
                          number_string[len] = '\0';
                          value = atol(number_string);
                        }
                      else
                        value = 0;

                      /* Shift string along on readiness for the next
                         number */
                      strcpy(tmp,string_ptr+1);
                      strcpy(tmp_string,tmp);
                    }

                  /* If no dot, then this is the last number */
                  else
                    {
                      /* Extract the value */
                      strcpy(number_string,tmp_string);
                      value = atol(number_string);

                      /* Set flag that we're done */
                      done = TRUE;
                    }

                  /* Store current value */
                  if (level < MAX_LEVEL)
                    {
                      params->instance[level] = value;
                      level++;
                      params->nlevels = level;
                    }
                }

              /* Shift token counter */
              itoken++;
            }
        }

      /* Check for the -f option giving the name of the file to
         which the HTML output is to be written */
      else if (!strcmp(strptrs[itoken],"-f"))
        {
          /* Get the file name from the next token, if there is one */
          if (itoken < num)
            {
              strcpy(params->output_file,strptrs[itoken + 1]);
              itoken++;
            }          
        }

      /* Check for the -chain option giving the chain identifier */
      else if (!strcmp(strptrs[itoken],"-chain"))
        {
          /* Get the chain id from the next token, if there is one */
          if (itoken < num)
            {
              params->chain = strptrs[itoken + 1][0];
              itoken++;
            }          
        }

      /* Check for the -db option giving the name of the database.dat
         file to be read in at the start */
      else if (!strcmp(strptrs[itoken],"-db"))
        {
          /* Get the file name from the next token, if there is one */
          if (itoken < num)
            {
              strcpy(params->database_dat,strptrs[itoken + 1]);
              itoken++;
            }          
        }

      /* Check for the -tdir option giving the HTML templates directory */
      else if (!strcmp(strptrs[itoken],"-trdir"))
        {
          /* Get the directory name from the next token, if there is one */
          if (itoken < num)
            {
              strcpy(params->template_real_dir,strptrs[itoken + 1]);
              itoken++;
            }          
        }

      /* Check for the -tdir option giving the HTML templates directory */
      else if (!strcmp(strptrs[itoken],"-tdir"))
        {
          /* Get the directory name from the next token, if there is one */
          if (itoken < num)
            {
              strcpy(params->template_dir,strptrs[itoken + 1]);
              itoken++;
            }          
        }

      /* Get template option */
      else if (!strcmp(strptrs[itoken],"-option"))
        {
          /* Get the option name from the next token, if there is one */
          if (itoken < num)
            {
              strcpy(params->option,strptrs[itoken + 1]);
              strcat(params->option,"_OPTION");
              itoken++;
            }          
        }

      /* Get the score passed from previous page */
      else if (!strcmp(strptrs[itoken],"-score"))
        {
          /* Get the score from the next token, if there is one */
          if (itoken < num)
            {
              strcpy(params->score,strptrs[itoken + 1]);
              itoken++;
            }          
        }

      /* Check for the -i option giving name of include file */
      else if (!strcmp(strptrs[itoken],"-i"))
        {
          /* Get the file name from the next token, if there is one */
          if (itoken < num)
            {
              strcpy(params->include_file,strptrs[itoken + 1]);
              itoken++;
            }          
        }

      /* Check for the -I option giving full name of include file */
      else if (!strcmp(strptrs[itoken],"-I"))
        {
          /* Get the file name from the next token, if there is one */
          if (itoken < num)
            {
              strcpy(params->include_file,strptrs[itoken + 1]);
              itoken++;

              /* Set flag that this is the full file name */
              params->fullname = TRUE;
            }          
        }

      /* Check for the -c option giving call number */
      else if (!strcmp(strptrs[itoken],"-c"))
        {
          /* Get the number from the next token, if there is one */
          if (itoken < num)
            {
              strcpy(number_string,strptrs[itoken + 1]);
              params->call_no = atoi(number_string);
              itoken++;
            }          
        }

      /* Check for the -s option giving start number */
      else if (!strcmp(strptrs[itoken],"-s"))
        {
          /* Get the number from the next token, if there is one */
          if (itoken < num)
            {
              strcpy(number_string,strptrs[itoken + 1]);
              params->start_no = atol(number_string);
              itoken++;
            }          
        }

      /* Check for the -enzymes option indicating that page is for
         Enzyme Structures Database */
      else if (!strcmp(strptrs[itoken],"-enzymes"))
        {
          /* Set flag */
          params->enzymes = TRUE;
        }

      /* Check for the -profunc option indicating that page is for
         ProFunc */
      else if (!strcmp(strptrs[itoken],"-profunc"))
        {
          /* Set flags */
          params->profunc = TRUE;
          params->pdb_type = PROFUNC;
        }

      /* Check for the -ligsearch option indicating that page is for
         LigSearch */
      else if (!strcmp(strptrs[itoken],"-ligsearch"))
        {
          /* Set flags */
          params->ligsearch = TRUE;
          params->pdb_type = LIGSEARCH;
        }

      /* Check for the -ligplus option indicating that page is for
         LigPlus licence entry system */
      else if (!strcmp(strptrs[itoken],"-ligplus"))
        {
          /* Set flags */
          params->ligplus = TRUE;
          params->pdb_type = LIGPLUS;
        }

      /* Check for the -disastr option indicating that page is for
         VarSite */
      else if (!strcmp(strptrs[itoken],"-disastr") ||
               !strcmp(strptrs[itoken],"-varsite"))
        {
          /* Set flags */
          params->varsite = TRUE;
          params->pdb_type = VARSITE;
        }

      /* Check for the -drat option indicating that page is for DraT */
      else if (!strcmp(strptrs[itoken],"-drat"))
        {
          /* Set flags */
          params->drat = TRUE;
          params->pdb_type = DRAT_DB;
        }

      /* Check for the -drugport option indicating that page is for
         DrugPort */
      else if (!strcmp(strptrs[itoken],"-drugport"))
        {
          /* Set flags */
          params->drugport = TRUE;
          params->pdb_type = DRUGPORT;
        }

      /* Check for the -sas option indicating that page is for SAS */
      else if (!strcmp(strptrs[itoken],"-sas"))
        {
          /* Set flags */
          params->sas = TRUE;
          params->pdb_type = SAS;
        }

      /* Check for the -highlights option indicating that PDBsum
         highlights page to be shown */
      else if (!strcmp(strptrs[itoken],"-highlights"))
        {
          /* Set flag */
          params->highlights = TRUE;
          strcpy(params->pdb_code,"n/a");
        }

      /* Check for the -e option indicating an error condition */
      else if (!strcmp(strptrs[itoken],"-e"))
        {
          /* Set flag */
          params->have_error = TRUE;
        }

      /* Check for the -u option giving user identifier */
      else if (!strcmp(strptrs[itoken],"-u"))
        {
          /* Get the number from the next token, if there is one */
          if (itoken < num)
            {
              store_string(&(params->user_id),strptrs[itoken + 1]);
              itoken++;
            }          
        }

      /* Check for the -parN option signifying parameter N */
      else if (!strncmp(strptrs[itoken],"-par",4))
        {
          /* Identify this parameter */
          strcpy(number_string,strptrs[itoken]+4);
          i = atoi(number_string) - 1;

          /* Get the number from the next token, if there is one */
          if (itoken < num && i > -1 && i < NPAR)
            {
              strcpy(params->parameter[i],strptrs[itoken + 1]);
              itoken++;
            }          
        }

      /* Check for the -code option giving security code */
      else if (!strcmp(strptrs[itoken],"-code"))
        {
          /* Get the code from the next token, if there is one */
          if (itoken < num)
            {
              strcpy(params->security_code,strptrs[itoken + 1]);
              itoken++;
            }          
        }

      /* Check for the -dir option giving the working directory */
      else if (!strcmp(strptrs[itoken],"-dir"))
        {
          /* Get the name from the next token, if there is one */
          if (itoken < num)
            {
              strcpy(params->dir_name,strptrs[itoken + 1]);
              itoken++;
            }          
        }

      /* Check for the -wkdir option giving special working directory */
      else if (!strcmp(strptrs[itoken],"-wkdir"))
        {
          /* Get the name from the next token, if there is one */
          if (itoken < num)
            {
              strcpy(params->work_dir,strptrs[itoken + 1]);
              itoken++;
            }          
        }

      /* Check for the -drug_id option giving DrugBank identifier */
      else if (!strcmp(strptrs[itoken],"-drug_id"))
        {
          /* Get the name from the next token, if there is one */
          if (itoken < num)
            {
              strcpy(params->drug_id,strptrs[itoken + 1]);
              itoken++;
            }          
        }

      /* Check for the -uniprot_acc option giving UniProt accession
         code */
      else if (!strcmp(strptrs[itoken],"-uniprot_acc"))
        {
          /* Get the name from the next token, if there is one */
          if (itoken < num)
            {
              strcpy(params->uniprot_acc,strptrs[itoken + 1]);
              itoken++;
            }          
        }

      /* Check for the -ligand_ranges option giving ligands selected from
         the ligand clusters page */
      else if (!strcmp(strptrs[itoken],"-ligand_ranges"))
        {
          /* Get the name from the next token, if there is one */
          if (itoken < num)
            {
              store_string(&(params->ligand_ranges),strptrs[itoken + 1]);
              itoken++;
            }          
        }

      /* Check for the -link_path parameter */
      else if (!strcmp(strptrs[itoken],"-link_path"))
        {
          /* Get the name from the next token, if there is one */
          if (itoken < num)
            {
              store_string(&(params->link_path),strptrs[itoken + 1]);
              itoken++;
            }          
        }

      /* Check for the UPLOAD option */
      else if (!strcmp(strptrs[itoken],"UPLOAD"))
        params->pdb_type = UPLOAD_PDB;

      /* Check for the MCSG option */
      else if (!strcmp(strptrs[itoken],"MCSG"))
        params->pdb_type = MCSG_PDB;

      /* Check for the -EBI option */
      else if (!strcmp(strptrs[itoken],"-EBI"))
        params->from_ebi = TRUE;

      /* Check for the -prop option, indicating this is a Prioprietary
         installation of PDBsum */
      else if (!strcmp(strptrs[itoken],"-prop"))
        params->proprietary = TRUE;

      /* Check for the -list option, to list all the fields for
         debugging */
      else if (!strcmp(strptrs[itoken],"-list"))
        params->list_fields = TRUE;

      /* Check for the -pdbsum1 option, signifying the standalone
         version is being run */
      else if (!strcmp(strptrs[itoken],"-pdbsum1"))
        params->pdb_type = PDBSUM1;

      /* Check for the -relinks option, signifying PDBsum1 links are
         to be relative rather than absolute */
      else if (!strcmp(strptrs[itoken],"-relinks"))
        params->relinks = TRUE;
    }

  /* For enzymes database, save EC number */
  if (params->enzymes == TRUE)
    {
      params->bottom_level = TRUE;
      for (i = 0; i < MAX_LEVEL && i < 4; i++)
        {
          params->ec_number[i] = params->instance[i];
          params->instance[i] = 0;
          if (params->ec_number[i] == 0)
            params->bottom_level = FALSE;
        }
      for (i = 0; i < params->nlevels + 1 && i < MAX_LEVEL; i++)
        params->instance[i] = 1;
      params->instance[0] = params->nlevels;
    }
}
/***********************************************************************

get_subscripts  -  Retrieve any subscripts associated with the current
                   field name

***********************************************************************/

int get_subscripts(char *string,int *instance,int have_pread,int *loop,
                   int plevel)
{
  char number_string[20], tmp_string[LINELEN];

  char *string_ptr, *string1_ptr;

  int i, level;
  int done;

  /* Initialise variables */
  done = FALSE;
  for (level = 0; level < MAX_LEVEL; level++)
    instance[level] = 0;
  level = 0;

  /* If field is one that has just been read in using the pread command,
     then want to prefix each new field with the current loop info */
  if (have_pread == TRUE && level > -1)
    {
      /* Transfer existing loop info to new field */
      for (i = 0; i < plevel + 1; i++)
        {
          /* Store the value */
          instance[level] = loop[i];

          /* Increment level */
          level++;
        }
    }

  /* Copy current line */
  strcpy(tmp_string,string);
  string_ptr = tmp_string;

  /* Loop while still have numbers to extract */
  while (done == FALSE && level < MAX_LEVEL)
    {
      /* Get position of next number */
      string_ptr = strstr(string_ptr,"[");

      /* If present, then extract */
      if (string_ptr != NULL)
        {
          /* Extract the value */
          strcpy(number_string,string_ptr + 1);
          string1_ptr = strstr(number_string,"]");
          if (string1_ptr != NULL)
            string1_ptr[0] = '\0';

          /* Store the value */
          instance[level] = atoi(number_string);

          /* Increment level */
          level++;

          /* Shift search start-point */
          string_ptr = strstr(string_ptr,"]");
          if (string_ptr == NULL)
            done = TRUE;
        }

      /* Otherwise, we're done */
      else
        done = TRUE;
    }

  /* Return number of levels */
  return(level);
}
/***********************************************************************

create_field  -  Create a new field record

***********************************************************************/

struct field *create_field(char *field_id)
{
  int level;

  struct field *field_ptr;

  /* Create the field */
  field_ptr = (struct field*) malloc(sizeof(struct field));
  if (field_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for");
      printf(" struct field\n");
      printf("***        Program terminated with error.\n");
      exit (1);
    }

  /* Initialise fields */
  store_string(&(field_ptr->field_id),field_id);
  field_ptr->read = FALSE;
  field_ptr->nlevel = 0;
  field_ptr->value = &null;
  field_ptr->left_field_ptr = NULL;
  field_ptr->right_field_ptr = NULL;

  /* Return the new field */
  return(field_ptr);
}
/***********************************************************************

add_field  -  Add a new field to the binary tree

***********************************************************************/

struct field *add_field(char *field_id,struct field *field_ptr)
{
  struct field *new_field_ptr;

  /* Initialise new field */
  new_field_ptr = field_ptr;

  /* If this is the root field, update it */
  if (field_ptr->field_id[0] == '\0')
    store_string(&(field_ptr->field_id),field_id);

  /* If already have this field, then return its pointer so that it can
     be updated with new data */
  else if (!strcmp(field_ptr->field_id,field_id))
    new_field_ptr = field_ptr;

  /* If field id lower than current number, check left branch */
  else if (strcmp(field_ptr->field_id,field_id) > 0)
    {
      /* If left-hand child is empty, place field there */
      if (field_ptr->left_field_ptr == NULL)
        {
          /* Create a new field */
          new_field_ptr = create_field(field_id);

          /* Make this the left-hand child of current field */
          field_ptr->left_field_ptr = new_field_ptr;
        }

      /* Otherwise call current routine to continue search */
      else
        new_field_ptr
          = add_field(field_id,field_ptr->left_field_ptr);
    }

  /* If field id higher than current number, check right branch */
  else if (strcmp(field_ptr->field_id,field_id) < 0)
    {
      /* If right-hand child is empty, place field there */
      if (field_ptr->right_field_ptr == NULL)
        {
          /* Create a new field */
          new_field_ptr = create_field(field_id);

          /* Make this the right-hand child of current field */
          field_ptr->right_field_ptr = new_field_ptr;
        }

      /* Otherwise call current routine to continue search */
      else
        new_field_ptr
          = add_field(field_id,field_ptr->right_field_ptr);
    }

  /* Return current pointer */
  return(new_field_ptr);
}
/***********************************************************************

list_fields  -  List the stored fields (for testing purposes only)

***********************************************************************/

void list_fields(FILE *file_ptr,struct field *field_ptr)
{
  /* If left-hand child is not empty, go down it to write out all the
     fields lower than the current one */
  if (field_ptr->left_field_ptr != NULL)
    list_fields(file_ptr,field_ptr->left_field_ptr);

  /* Write out the current field id and value */
  fprintf(file_ptr,"%s %s\n",field_ptr->field_id,field_ptr->value);

  /* If have a right branch, go down it to write out all fields higher
     than this one */
  if (field_ptr->right_field_ptr != NULL)
    list_fields(file_ptr,field_ptr->right_field_ptr);
}
/***********************************************************************

create_field_record  -  Create a new field record from the given input
                        string

***********************************************************************/

struct field *create_field_record(struct field *first_field_ptr,
                                  char *string,int read,int have_pread,
                                  int *loop,int plevel)
{
  char field_id[NAMLEN + 1], name[NAMLEN + 1], number_string[20];

  char *string_ptr;

  int ipos, instance[MAX_LEVEL], len, level, nlevel, value_start;
  int found;

  struct field *field_ptr;

  /* Initialise variables */
  value_start = 0;
  len = strlen(string);

  /* Determine start of field value */
  string_ptr = strchr(string,' ');
  if (string_ptr != NULL)
    {
      /* Initialise flag */
      found = FALSE;

      /* Set value-start position */
      value_start = string_ptr - string + 1;

      /* Loop until first non-blank character located */
      for (ipos = value_start; ipos < len && found == FALSE; ipos++)
        {
          /* If this is a non-blank character, set start-point of value
             here */
          if (string[ipos] != ' ' && string[ipos] != '\0')
            {
              /* Set value-start and flag */
              value_start = ipos;
              found = TRUE;
            }
        }

      /* If non-blank character not found, then have no value */
      if (found == FALSE)
        value_start = 0;
    }

  /* Get field name from first part of string */
  if (strlen(string) > NAMLEN)
    {
      strncpy(field_id,string,NAMLEN);
      field_id[NAMLEN] = '\0';
    }
  else
    strcpy(field_id,string);

  /* Remove field value from end of string */
  string_ptr = strchr(field_id,' ');
  if (string_ptr != NULL)
    string_ptr[0] = '\0';

  /* Retrieve the subscript values */
  nlevel = get_subscripts(field_id,instance,have_pread,loop,plevel);

  /* Lop the subscripts off from the field name */
  string_ptr = strchr(field_id,'[');
  if (string_ptr != NULL)
    string_ptr[0] = '\0';

  /* Save the field name */
  strcpy(name,field_id);

  /* Form the field identifier by adding the instance identifiers,
     if any, to the field name */
  for (level = 0; level < nlevel; level++)
    {
      /* Form number string storing this instance number */
      sprintf(number_string,".%6.6d",instance[level]);
  
      /* Add on to end of field id */
      strcat(field_id,number_string);
    }

  /* Add field to binary tree */
  field_ptr = add_field(field_id,first_field_ptr);

  /* Initialise field details */
  field_ptr->read = read;

  /* Initialise record */
  field_ptr->nlevel = nlevel;

  /* Store field name and field value */
  if (field_ptr->value != &null)
    free(field_ptr->value);
  if (value_start > 0)
    {
      if (!strcmp(string + value_start,"/null/"))
        field_ptr->value = &null;
      else
        store_string(&(field_ptr->value),string + value_start);
    }
  else
    store_string(&(field_ptr->value),"NONE");

  /* Return pointer to field record just created */
  return(field_ptr);
}
/***********************************************************************

set_loop_level  -  Set the loop level according to the entered
                   command-line parameters

***********************************************************************/

int set_loop_level(int *loop,struct parameters *params)
{
  int ilevel, level;

  /* Initialise variables */
  level = -1;

  /* Initialise loop values at each level */
  for (ilevel = 0; ilevel < MAX_LEVEL; ilevel++)
    loop[ilevel] = 0;

  /* If generating a page for a specific chain or ligand, then
     starting level is up one, and current level corresponds to
     the selected chain/ligand */
  for (ilevel = 0; ilevel < MAX_LEVEL; ilevel++)
    {
      if (params->instance[ilevel] > 0)
        {
          /* Set starting level */
          level = ilevel;

          /* Set starting value */
          loop[level] = params->instance[ilevel];
        }
    }

  /* Return starting level */
  return(level);
}
/***********************************************************************

get_all_paths  -  Read in all the paths and environments in the
                  CATHPARAM parameter file

***********************************************************************/

void get_all_paths(struct field *first_field_ptr)
{
  char parameter[FILENAME_LEN + 1], value[FILENAME_LEN + 1];
  char string[LINELEN + 1];

  int have_error;

  struct field *field_ptr;

  /* Initialise variables */
  have_error = FALSE;

  /* Loop until all the parameters have been picked up from the
     CATHPARAM file */
  while (have_error == FALSE)
    {
      /* Retrieve the next parameter */
      have_error = get_next_param(parameter,value,FILENAME_LEN,TRUE);

      /* Create a field record for this parameter */
      if (have_error == FALSE)
        {
          /* From the string containing keyword-value pair */
          sprintf(string,"%s %s",parameter,value);

          /* Create a field record for this parameter */
          field_ptr
            = create_field_record(first_field_ptr,string,FALSE,FALSE,
                                  LOOP,-1);
        }
    }
}
/***********************************************************************

find_field  -  Find the given field in the binary tree

***********************************************************************/

struct field *find_field(char *field_id,struct field *field_ptr)
{
  struct field *found_field_ptr;

  /* Initialise found field */
  found_field_ptr = field_ptr;

  /* If field-id lower than current number, check left branch */
  if (strcmp(field_id,field_ptr->field_id) < 0)
    {
      /* If left-hand child is empty, then number not in tree */
      if (field_ptr->left_field_ptr == NULL)
        found_field_ptr = NULL;

      /* Otherwise call current routine to continue search */
      else
        found_field_ptr
          = find_field(field_id,field_ptr->left_field_ptr);
    }

  /* If field_id higher than current number, check right branch */
  else if (strcmp(field_id,field_ptr->field_id) > 0)
    {
      /* If right-hand child is empty, then number not in tree */
      if (field_ptr->right_field_ptr == NULL)
        found_field_ptr = NULL;

      /* Otherwise call current routine to continue search */
      else
        found_field_ptr
          = find_field(field_id,field_ptr->right_field_ptr);
    }

  /* Return current pointers */
  return(found_field_ptr);
}
/***********************************************************************

get_field  -  Determine the value of the given field

***********************************************************************/

struct field *get_field(struct field *first_field_ptr,char *name,
                        int *loop,int *instance)
{
  char field_id[LINELEN + 1], number_string[20];

  int i, level, max_level, param_top, top_level, try;
  int done;

  struct field *field_ptr;

  /* Initialise variables */
  done = FALSE;
  field_ptr = NULL;

  /* First try searching on just the name */
  strcpy(field_id,name);
  field_ptr = find_field(field_id,first_field_ptr);

  /* If found, then we are done */
  if (field_ptr != NULL)
    return(field_ptr);

  /* Determine the top level from supplied parameters */
  param_top = -1;
  top_level = -1;
  
  /* Determine what the top levels are */
  for (level = 0; level < MAX_LEVEL; level++)
    {
      /* If current level is non-zero update top level */
      if (loop[level] > 0)
        top_level = level;

      /* Repeat for parameter-supplied level */
      if (instance[level] > 0)
        param_top = level;
    }

  /* Determine which is the higher of the two levels */
  max_level = param_top;
  if (top_level > max_level)
    max_level = top_level;

  /* Loop over attempts to match field, first to parameter and
     current loop levels, and then, if that fails try again with
     first level set to zero (for Enzyme reactions display in PDBsum) */
  for (try = 0; try < 2 && field_ptr == NULL; try++)
    {
      /* Loop over possible searches, starting from most specific */
      for (level = max_level; level > -1 && field_ptr == NULL; level--)
        {
          /* Re-initialise field-id */
          strcpy(field_id,name);

          /* Loop to add all subscripts up to this level */
          for (i = 0; i < level + 1; i++)
            {
              /* Store loop number or instance */
              if (try == 0 || (try == 1 && i > 0))
                {
                  /* If at a parameter level, use that */
                  if (loop[i] <= 0)
                    sprintf(number_string,".%6.6d",instance[i]);
                  else
                    sprintf(number_string,".%6.6d",loop[i]);
                }
              else
                sprintf(number_string,".%6.6d",0);
  
              /* Add on to end of field id */
              strcat(field_id,number_string);
            }

          /* Search for the given field */
          field_ptr = find_field(field_id,first_field_ptr);
        }
    }

  /* If field found, but its value is "FALSE" or "NONE", then treat
     as though the field does not exist */
  if (field_ptr != NULL && (!strcmp(field_ptr->value,"FALSE") ||
                            !strcmp(field_ptr->value,"NONE")))
    field_ptr = NULL;

  /* Return field */
  return(field_ptr);
}
/***********************************************************************

create_pdir  -  Create PDIR, the name of the PDBsum directory for this
                PDB code

***********************************************************************/

void create_pdir(struct field *first_field_ptr,char *string,char *tmp,
                 struct parameters *params)
{
  struct field *field_ptr;

  /* Initialise variables */
  field_ptr = NULL;
  if (params->pdb_code[0] == '\0')
    return;

  /* Get the appropriate directory template */
  if (params->pdb_type == UPLOAD_PDB || params->pdb_type == MCSG_PDB)
    field_ptr
        = get_field(first_field_ptr,"UPLOAD_PDBSUM_TEMPLATE",LOOP,
                    params->instance);
  else
    field_ptr
        = get_field(first_field_ptr,"PDBSUM_REAL_TEMPLATE",LOOP,
                    params->instance);

  /* If found, create a PDIR record */
  if (field_ptr != NULL)
    {
      /* Form the full directory name */
      tmp = form_filename(params->pdb_code,field_ptr->value);

      /* Form the string to be saved */
      sprintf(string,"PDIR %s",tmp);

      /* Create a field record */
      field_ptr
        = create_field_record(first_field_ptr,string,FALSE,FALSE,
                              LOOP,-1);
    }

  /* Repeat for the URL corresponding to the PDBsum directory */
  if (params->pdb_type == STANDARD_PDB ||
      params->pdb_type == ALPHA_FOLD || params->pdb_type == PDBSUM1)
    {
      /* Get the template */
      field_ptr
        = get_field(first_field_ptr,"PDBSUM_TEMPLATE",LOOP,
                    params->instance);

      /* If found, create a PDBSUM_DIR record */
      if (field_ptr != NULL)
        {
          /* Form the full directory name */
          tmp = form_filename(params->pdb_code,field_ptr->value);

          /* Form the string to be saved */
          sprintf(string,"PDBSUM_DIR %s",tmp);

          /* Create a field record */
          field_ptr
            = create_field_record(first_field_ptr,string,FALSE,FALSE,
                                  LOOP,-1);
        }
    }
}
/***********************************************************************

find_pdb_file  -  Locate the PDB file from its code, trying the
                  standard PDB directories, or MCSG and Uploads
                  directories, as appropriate

***********************************************************************/

int find_pdb_file(char *pdb_code,struct field *first_field_ptr,
                  char *file_name,int pdb_type,int *instance)
{
  char *pdb_template[NPDB_DIR];

  char *fname;

  int dir, idir, start;

  struct field *field_ptr;

  FILE *file_ptr;

  /* Initialise */
  file_name[0] = '\0';
  dir = NOT_FOUND;

  /* Loop through the fields from CATHPARAM to pick up and store all
     the PDB template names */
  for (idir = 0; idir < NPDB_DIR; idir++)
    {
      /* Initialise the template name */
      pdb_template[idir] = &null;

      /* Get this template */
      field_ptr
        = get_field(first_field_ptr,PDB_DIR_TEMPLATE[idir],LOOP,
                    instance);

      /* If template found, then save */
      if (field_ptr != NULL)
        store_string(&(pdb_template[idir]),field_ptr->value);
    }

  /* Determine the starting PDB directory to search */
  start = 0;
  if (pdb_type == UPLOAD_PDB)
    start = PDB_UPLOAD;
  else if (pdb_type == MCSG_PDB)
    start = PDB_MCSG;
  else if (pdb_type == ALPHA_FOLD)
    start = PDB_ALPHA_FOLD;

  /* Loop to locate the PDB file */
  for (idir = start; idir < NPDB_DIR && dir == NOT_FOUND; idir++)
    {
      /* Form the trial name */
      fname = form_filename(pdb_code,pdb_template[idir]);

      /* Open the file, if it exists */
      if ((file_ptr = fopen(fname,"r")) !=  NULL)
        {
          /* Set the appropriate flag */
          dir = idir;

          /* Save the file name */
          strcpy(file_name,fname);

          /* Close the file */
          fclose(file_ptr);
        }

      /* Free memory of file name */
      free(fname);
    }

  /* Return the file type */
  return(dir);
}
/***********************************************************************

get_pdb_filename  -  Get the name of the PDB file

***********************************************************************/

void get_pdb_filename(struct field *first_field_ptr,char *string,
                      char *file_name,struct parameters params)
{
  int found;

  struct field *field_ptr;

  /* Initialise variables */
  field_ptr = NULL;
  found = FALSE;
  if (params.pdb_code[0] == '\0')
    return;

  /* Find this PDB file */
  if (strcmp(params.pdb_code,"none") && strcmp(params.pdb_code,"n/a"))
      found = find_pdb_file(params.pdb_code,first_field_ptr,file_name,
                            params.pdb_type,params.instance);

  /* If found, create a PDIR record */
  if (found == NOT_FOUND)
    strcpy(file_name,"NONE");

  /* Form the string to be saved */
  sprintf(string,"PDB_FILENAME %s",file_name);

  /* Create a field record */
  field_ptr
    = create_field_record(first_field_ptr,string,FALSE,FALSE,LOOP,-1);
}
/***********************************************************************

create_vdir  -  Create VDIR, the name of the variants directory for this
                PDB code

***********************************************************************/

void create_vdir(struct field *first_field_ptr,char *string,char *tmp,
                 struct parameters *params)
{
  struct field *field_ptr;

  /* Initialise variables */
  field_ptr = NULL;
  if (params->pdb_code[0] == '\0')
    return;

  /* If an upload file, then won't have variants data */
  if (params->pdb_type == UPLOAD_PDB || params->pdb_type == MCSG_PDB)
    return;

  /* Get the variants template */
  field_ptr
    = get_field(first_field_ptr,"VARIANTS_TEMPLATE",LOOP,
                params->instance);

  /* If found, create a VDIR record */
  if (field_ptr != NULL)
    {
      /* Form the full directory name */
      tmp = form_filename(params->pdb_code,field_ptr->value);

      /* Form the string to be saved */
      sprintf(string,"VDIR %s",tmp);

      /* Create a field record */
      field_ptr
        = create_field_record(first_field_ptr,string,FALSE,FALSE,
                              LOOP,-1);
    }
}
/***********************************************************************

read_database_file  -  Read in all the data fields for this PDB entry
                       from the database.dat file in its PDBsum directory

***********************************************************************/

int read_database_file(struct field *first_field_ptr,
                       char *database_file,int full_name,int read,
                       int have_pread,char *key,long offset,
                       int start_record,int nrecords,int *loop,
                       int level,char *replace_from,char *replace_by,
                       struct parameters *params)
{
  char dir_name[FILENAME_LEN + 1], file_name[FILENAME_LEN + 1];
  char split_dir[3], input_line[LINELEN + 1], tmp_line[LINELEN + 1];
  char subdir[FILENAME_LEN + 1], string[FILENAME_LEN + 1];
  char number_string[20], subdir_dir[FILENAME_LEN + 1];

  char *string_ptr, *tmp;

  int end_record, len, ndata, replen, irecord;
  int done, error, have_fields, have_file, replace, selected, wanted;

  struct field *field_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  if (nrecords == -1)
    end_record = -1;
  else
    end_record = start_record + nrecords - 1;
  ndata = 0;
  replace = FALSE;
  replen = 0;
  if (replace_from[0] != '\0' && replace_by[0] != '\0')
    {
      replace = TRUE;
      replen = strlen(replace_from);
    }

  /* Initialise field pointers */
  field_ptr = NULL;
  file_name[0] = '\0';

  /* If have complete file name, use that */
  if (full_name == TRUE)
    strcpy(file_name,database_file);

  /* If generating page for Enzyme Structures database, form full
     filename */
  else if (params->enzymes == TRUE)
    {
      /* Get path to enzymes directory */
      field_ptr = get_field(first_field_ptr,"ENZYME_REAL_DIR",LOOP,
                            params->instance);

      /* Form name of the database.dat file */
      if (field_ptr != NULL)
        {
          strcpy(file_name,field_ptr->value);
          strcat(file_name,DIR_SEP);
        }

      /* Form appropriate directory name */
      if (params->nlevels == 1)
        sprintf(dir_name,"ec%1d",params->ec_number[0]);
      else if (params->nlevels == 2)
        sprintf(dir_name,"ec%1d%sec%2.2d",params->ec_number[0],DIR_SEP,
                params->ec_number[1]);
      else if (params->nlevels == 3)
        sprintf(dir_name,"ec%1d%sec%2.2d%sec%2.2d",params->ec_number[0],
                DIR_SEP,params->ec_number[1],DIR_SEP,params->ec_number[2]);
      else if (params->nlevels == 4)
        sprintf(dir_name,"ec%1d%sec%2.2d%sec%2.2d%sec%4.4d",
                params->ec_number[0],DIR_SEP,params->ec_number[1],DIR_SEP,
                params->ec_number[2],DIR_SEP,params->ec_number[3]);

      /* Append file names */
      strcat(file_name,dir_name);
      strcat(file_name,DIR_SEP);
      strcat(file_name,database_file);
    }

  /* If generating PDBsum highlights page, form full filename */
  else if (params->highlights == TRUE)
    {
      /* Determine which PDBsum data directory is the appropriate one */
      if (params->pdb_type == STANDARD_PDB)
        strcpy(dir_name,"PDBSUM_REAL_DATA_DIR");

      /* Get path to PDBsum data directory */
      field_ptr = get_field(first_field_ptr,dir_name,LOOP,
                            params->instance);

      /* Form name of the highlights.dat file */
      if (field_ptr != NULL)
        {
          strcpy(file_name,field_ptr->value);
          strcat(file_name,DIR_SEP);
          strcat(file_name,"highlights.dat");
        }
    }

  /* For ProFunc results database.dat will either be in the temporary
     directory or in the PDBsum upload directory, or in the PDBsum
     directory */
  else if (params->profunc == TRUE)
    {
      /* If all data in the temporary directory (ie old system), define
         names of both directories here */
      if (params->user_id != &null && !strcmp(params->pdb_code,"none"))
        {
          /* Get the name of the temporary directory */
          field_ptr = get_field(first_field_ptr,"TMP_DIR",LOOP,
                                params->instance);
          if (field_ptr != NULL)
            {
              /* Form full name of database.dat file */
              strcpy(file_name,field_ptr->value);
              strcat(file_name,DIR_SEP);
              strcat(file_name,"profunc");
              strcat(file_name,params->user_id);
              strcat(file_name,DIR_SEP);
              strcat(file_name,"database.dat");
            }
          else
            strcpy(file_name,"database.dat");

          /* Create a flag indicating that data in temporary dir */
          sprintf(string,"IN_TMP_DIR TRUE");
          field_ptr
            = create_field_record(first_field_ptr,string,FALSE,FALSE,
                                  LOOP,-1);
        }

      /* Otherwise, PDBsum directory is in the uploads area, or in
         the PDBsum directory */
      else
        {
          /* Get the PDBsum directory */
          field_ptr
            = get_field(first_field_ptr,"PDIR",LOOP,params->instance);

          /* Form name of the database.dat file */
          if (field_ptr != NULL)
            {
              strcpy(file_name,field_ptr->value);
              strcat(file_name,DIR_SEP);
              strcat(file_name,"database.dat");
            }
          else
            strcpy(file_name,"database.dat");
        }
    }

  /* For LigSearch results database.dat will either be in the temporary
     directory or in the LigSearch directory */
  else if (params->ligsearch == TRUE)
    {
      /* If all data in the temporary directory, define names of both
         directories here */
      if (params->user_id != &null &&
          (!strcmp(params->pdb_code,"none") ||
           !strcmp(params->pdb_code,"n/a")))
        {
          /* Get the name of the temporary directory */
          field_ptr = get_field(first_field_ptr,"TMP_DIR",LOOP,
                                params->instance);
          if (field_ptr != NULL)
            {
              /* Form full name of database.dat file */
              strcpy(file_name,field_ptr->value);
              strcat(file_name,DIR_SEP);
              strcat(file_name,"ligsearch");
              strcat(file_name,params->user_id);
              strcat(file_name,DIR_SEP);
              strcat(file_name,"database.dat");
            }
          else
            strcpy(file_name,"database.dat");
        }

      /* Otherwise, LigSearch directory is in the uploads area */
      else
        {
          /* Get the LigSearch directory */
          field_ptr
            = get_field(first_field_ptr,"LIGSEARCH_REAL_RESULTS",LOOP,
                        params->instance);

          /* Form name of the database.dat file */
          if (field_ptr != NULL)
            {
              /* Get the length of the code */
              len = strlen(params->pdb_code);

              /* If code is 4 characters long, take it to be a
                 PDB code */
              if (len == 4)
                {
                  /* Get middle 2 digits of the PDB code */
                  strncpy(split_dir,params->pdb_code+1,2);
                  split_dir[2] = '\0';
                }

              /* Otherwise, use the last two characters */
              else
                {
                  /* Get last 2 digits of the UniProt code */
                  strncpy(split_dir,params->pdb_code + len - 2,2);
                  split_dir[2] = '\0';
                }

              /* Form full name of database.dat file */
              strcpy(file_name,field_ptr->value);
              strcat(file_name,DIR_SEP);
              strcat(file_name,split_dir);
              strcat(file_name,DIR_SEP);
              strcat(file_name,params->pdb_code);
              strcat(file_name,DIR_SEP);
              strcat(file_name,"database.dat");

              /* Save the name of the split directory */
              sprintf(string,"SUBDIR %s",split_dir);
              field_ptr
                = create_field_record(first_field_ptr,string,FALSE,FALSE,
                                      LOOP,-1);
            }
          else
            strcpy(file_name,"database.dat");
        }
    }

  /* For VarSite, results database.dat will either be in the temporary
     directory or in the VarSite directory */
  else if (params->varsite == TRUE)
    {
      /* If all data in the temporary directory, define names of both
         directories here */
      if (params->user_id != &null &&
          params->uniprot_acc[0] == '\0')
        {
          /* Get the name of the temporary directory */
          field_ptr = get_field(first_field_ptr,"TMP_DIR",LOOP,
                                params->instance);
          if (field_ptr != NULL)
            {
              /* Form full name of database.dat file */
              strcpy(file_name,field_ptr->value);
              strcat(file_name,DIR_SEP);
              strcat(file_name,"varsite");
              strcat(file_name,params->user_id);
              strcat(file_name,DIR_SEP);
              strcat(file_name,"database.dat");
            }
          else
            strcpy(file_name,"database.dat");
        }

      /* Otherwise, VarSite directory is in the data area */
      else
        {
          /* Get the VarSite directory */
          field_ptr
            = get_field(first_field_ptr,"VARSITE_REAL_UNIPROT_DIR",LOOP,
                        params->instance);

          /* Form name of the database.dat file */
          if (field_ptr != NULL)
            {
              /* Get the length of the code */
              len = strlen(params->uniprot_acc);

              /* Get last 2 digits of the UniProt code */
              strncpy(split_dir,params->uniprot_acc + len - 2,2);
              split_dir[2] = '\0';

              /* Form full name of database.dat file */
              strcpy(file_name,field_ptr->value);
              strcat(file_name,DIR_SEP);
              strcat(file_name,split_dir);
              strcat(file_name,DIR_SEP);
              strcat(file_name,params->uniprot_acc);
              strcat(file_name,DIR_SEP);
              strcpy(subdir,file_name);
              strcat(file_name,"database.dat");

              /* Save the name of the split directory */
              sprintf(string,"SUBDIR %s",split_dir);
              field_ptr
                = create_field_record(first_field_ptr,string,FALSE,FALSE,
                                      LOOP,-1);

              /* Save the name of this entry's directory */
              sprintf(string,"DDIR %s",subdir);
              field_ptr
                = create_field_record(first_field_ptr,string,FALSE,FALSE,
                                      LOOP,-1);
            }
          else
            strcpy(file_name,"database.dat");
        }
    }

  /* For DrugPort, read either the drug's or the target's database.dat
     file */
  else if (params->drugport == TRUE)
    {
      /* Get the name of the DrugBank data directory */
      field_ptr = get_field(first_field_ptr,"DRUGPORT_REAL_DIR",LOOP,
                            params->instance);

      /* If DrugBank id, form name of appropriate subdirectory */
      if (field_ptr != NULL && params->drug_id[0] != '\0')
        {
          /* Get last 2 digits of the drug id */
          len = strlen(params->drug_id);
          if (len >= 2)
            strncpy(split_dir,params->drug_id+len-2,2);
          else if (len == 1)
            sprintf(split_dir,"0%c",params->drug_id[0]);
          else
            return(0);
          split_dir[2] = '\0';

          /* Form name of subdirectory */
          strcpy(subdir,DIR_SEP);
          strcpy(subdir_dir,"drugs");
          strcat(subdir_dir,DIR_SEP);
          strcat(subdir_dir,split_dir);
          strcat(subdir_dir,DIR_SEP);
          strcat(subdir_dir,params->drug_id);
          strcat(subdir,subdir_dir);
          strcat(subdir,DIR_SEP);

          /* Form name of subdirectory and full name of the
             database.new file */
          strcpy(file_name,field_ptr->value);
          strcat(file_name,subdir);
          strcpy(params->subdir_dir,subdir_dir);
          strcpy(params->subdir_real,file_name);
          strcat(file_name,"database.dat");

          /* Get the name of the DrugBank data directory */
          field_ptr = get_field(first_field_ptr,"DRUGPORT_DIR",LOOP,
                                params->instance);

          /* Store name of subdirectory for images and .dat files */
          if (field_ptr != NULL)
            {
              strcpy(params->subdir,field_ptr->value);
              strcat(params->subdir,subdir);
            }
        }

      /* If have target id, form name of appropriate subdirectory */
      else if (field_ptr != NULL && params->uniprot_acc[0] != '\0')
        {
          /* Get last 2 digits of the drug id */
          len = strlen(params->uniprot_acc);
          if (len >= 2)
            strncpy(split_dir,params->uniprot_acc+len-2,2);
          else if (len == 1)
            sprintf(split_dir,"0%c",params->uniprot_acc[0]);
          else
            return(0);
          split_dir[2] = '\0';

          /* Form name of subdirectory */
          strcpy(subdir,DIR_SEP);
          strcpy(subdir_dir,"targets");
          strcat(subdir_dir,DIR_SEP);
          strcat(subdir_dir,split_dir);
          strcat(subdir_dir,DIR_SEP);
          strcat(subdir_dir,params->uniprot_acc);
          strcat(subdir,subdir_dir);
          strcat(subdir,DIR_SEP);

          /* Form name of subdirectory and full name of the
             database.new file */
          strcpy(file_name,field_ptr->value);
          strcat(file_name,subdir);
          strcpy(params->subdir_dir,subdir_dir);
          strcpy(params->subdir_real,file_name);
          strcat(file_name,"database.dat");

          /* Get the name of the DrugBank data directory */
          field_ptr = get_field(first_field_ptr,"DRUGPORT_DIR",LOOP,
                                params->instance);

          /* Store name of subdirectory for images and .dat files */
          if (field_ptr != NULL)
            {
              strcpy(params->subdir,field_ptr->value);
              strcat(params->subdir,subdir);
            }
        }
    }

  /* If looking for the variants.dat file, find it in the appropriate
     directory */
  else if (!strcmp(database_file,"variants.dat"))
    {
      /* Get the template for the variants directory */
      field_ptr
        = get_field(first_field_ptr,"VARIANTS_REAL_TEMPLATE",LOOP,
                    params->instance);

      /* If found, form the full directory name */
      if (field_ptr != NULL)
        {
          /* Form the full directory name */
          tmp = form_filename(params->pdb_code,field_ptr->value);

          /* Append name of the file */
          strcpy(file_name,tmp);
          strcat(file_name,DIR_SEP);
          strcat(file_name,database_file);
        }
    }

  /* Otherwise, for PDBsum database */
  else
    {
      /* Get the PDBsum directory */
      field_ptr
        = get_field(first_field_ptr,"PDIR",LOOP,params->instance);

      /* Form name of the database.dat file */
      if (field_ptr != NULL)
        {
          strcpy(file_name,field_ptr->value);
          strcat(file_name,DIR_SEP);
          strcat(file_name,database_file);
        }

      /* If have an invalid PDB code, then no file to find */
      if (!strcmp(params->pdb_code,"...."))
        file_name[0] = '\0';
    }

  /* Open the database.dat file */
  if (file_name[0] != '\0' && (file_ptr = fopen(file_name,"r")) !=  NULL)
    {
      /* Initialise flags */
      done = FALSE;
      wanted = FALSE;
      have_fields = FALSE;

      /* If have file offset, then start read at that position */
      if (offset > -1)
        {
          /* Find start-point in data file using fseek */
          error = fseek(file_ptr,offset,SEEK_SET);

          /* Define a dummy key */
          strcpy(key,"KEY: N/A");
        }

      /* Loop while reading in records from the file */
      while (fgets(input_line,LINELEN,file_ptr) != NULL && done == FALSE)
        {
          /* Truncate the string to strip off the newline */
          string_chop(input_line);

          /* If have a KEY, check whether this is a key record */
          if (key[0] != '\0' && !strncmp(input_line,"KEY: ",5))
            {
              /* If we already have some data fields, then we're done */
              if (have_fields == TRUE)
                {
                  /* Set flags */
                  done = TRUE;
                  wanted = FALSE;
                }

              /* Otherwise, check whether this is our key */
              else if (!strcmp(input_line+5,key) || offset > -1)
                wanted = TRUE;

              /* Convert into a comment line */
              input_line[0] = '#';
            }

          /* If have a key and are not reading in its data, make line
             a comment line */
          else if (key[0] != '\0' && wanted == FALSE)
              input_line[0] = '#';

          /* Otherwise if this is a not-key field, mark that we have
             some data fields for our key */
          else if (key[0] != '\0' && wanted == TRUE)
            have_fields = TRUE;

          /* Process if not a comment line */
          if (input_line[0] != '#')
            {
              /* Initialise flag */
              selected = TRUE;

              /* If we are selecting records, check this record's
                 subscript */
              if (end_record > -1)
                {
                  /* Check for first subscript */
                  string_ptr = strstr(input_line,"[");

                  /* If present, then extract */
                  if (string_ptr != NULL)
                    {
                      /* Extract the value */
                      strcpy(tmp_line,string_ptr + 1);
                      string_ptr = strstr(tmp_line,"]");
                      if (string_ptr != NULL)
                        string_ptr[0] = '\0';
                      irecord = atoi(tmp_line);

                      /* Check whether record falls within required
                         range */
                      if (irecord < start_record || irecord > end_record)
                        selected = FALSE;

                      /* Otherwise, need to replace initial subscript
                         by sequential number commencing at start-recor */
                      else
                        {
                          /* Determine the new subscript number */
                          irecord = irecord - start_record + 1;
                          sprintf(number_string,"%d",irecord);

                          /* Get the field name */
                          strcpy(tmp_line,input_line);

                          /* Chop a start of subscript */
                          string_ptr = strstr(tmp_line,"[");
                          if (string_ptr != NULL)
                            string_ptr[1] = '\0';

                          /* Append new subscript */
                          strcat(tmp_line,number_string);

                          /* Get rest of input line */
                          string_ptr = strstr(input_line,"]");
                          if (string_ptr != NULL)
                            strcat(tmp_line,string_ptr);

                          /* Transfer back to input line */
                          strcpy(input_line,tmp_line);
                        }
                    }
                }

              /* If field is wanted, save it */
              if (selected == TRUE)
                {
                  /* Increment count of lines read in */
                  ndata++;

                  /* If replacing part of field_name check if have string to
                     be replaced */
                  if (replace == TRUE)
                    {
                      /* If we have string to be replaced, splice in its
                         replacement */
                      if (!strncmp(input_line,replace_from,replen))
                        {
                          /* Get string to be spliced out */
                          strncpy(tmp_line,input_line,replen);
                          tmp_line[replen] = '\0';

                          /* Perform the splice */
                          perform_splice(input_line,tmp_line,replace_by);
                        }
                    }

                  /* Create a record for this data item */
                  field_ptr
                    = create_field_record(first_field_ptr,input_line,
                                          read,have_pread,loop,level);
                }
            }
        }

      /* Close the file */
      fclose(file_ptr);
    }

  /* No database.dat file found, show warning */
  else
    {
      if (file_name[0] != '\0' && !strcmp(database_file,"database.dat") &&
          strcmp(params->pdb_code,"n/a"))
        {
          printf("*** Warning. database.dat file not found ...\n");
          printf("             [%s]\n",file_name);
        }
    }

  /* Check that PDB code hasn't been corrupted */
  field_ptr = get_field(first_field_ptr,"PDB_CODE",LOOP,params->instance);
  if (field_ptr != NULL && !strcmp(field_ptr->value,"none") &&
      strcmp(params->pdb_code,"n/a") && strcmp(params->pdb_code,"none")
      && strcmp(params->pdb_code,"...."))
    {
      /* Reinstate PDB code overwritten by invalid value in file */
      sprintf(string,"PDB_CODE %s",params->pdb_code);
      field_ptr
        = create_field_record(first_field_ptr,string,FALSE,FALSE,
                              LOOP,-1);
    }

  /* Return number of data items read in */
  return(ndata);
}
/***********************************************************************

store_motif_count  -  Create a recording storing the given motif count

***********************************************************************/

struct field *store_motif_count(struct field *first_field_ptr,
                                char *motif_name,int line,int count)
{
  char string[LINELEN + 1];

  struct field *field_ptr;

  /* Initialise field pointers */
  field_ptr = NULL;

  /* From the string containing keyword-value pair */
  if (line > 0)
    sprintf(string,"N%s[%d] %d",motif_name,line,count);
  else
    sprintf(string,"N%s %d",motif_name,count);

  /* Create a field record */
  field_ptr
    = create_field_record(first_field_ptr,string,FALSE,FALSE,LOOP,-1);

  /* Return pointer to created record */
  return(field_ptr);
}
/***********************************************************************

read_dbt_file  -  Read in all the double-colon-separated data fields
                  from given .dbt file

***********************************************************************/

int read_dbt_file(struct field *first_field_ptr,char *file_name,
                  int read,int promotif,long offset,int nrecords,
                  int *loop,struct parameters *params)
{
  char input_line[LINELEN + 1], field[LINELEN + 1], subfield[LINELEN + 1];
  char ch, chain, inchn, name[LINELEN + 1], number_string[20];
  char key[LINELEN + 1], tmp[LINELEN + 1];
  char motif_name[LINELEN + 1], subscripts[20], string[LINELEN + 1];

  char *field_name[MAX_TOKEN], *string_ptr;

  int field_len, i, ipos, item, jpos, len, line, nnames, subitem;
  int fst_pos, lpos, ndata, nstored;
  int all_done, done, finished, get_all, have_semicolons, started, wanted;
  int error, have_field_names;

  struct field *field_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  all_done = FALSE;
  get_all = FALSE;
  key[0] = '\0';
  have_field_names = FALSE;
  line = 0;
  motif_name[0] = '\0';
  ndata = 0;
  nnames = 0;
  nstored = 0;
  started = TRUE;
  wanted = TRUE;

  /* Initialise field names */
  for (i = 0; i < MAX_TOKEN; i++)
    field_name[i] = &null;

  /* Initialise field pointers */
  field_ptr = NULL;

  /* If reading in promotif.dat file, get the current chain */
  if (promotif == TRUE)
    {
      /* Set started flag and intialise chain */
      started = FALSE;
      chain = ' ';

      /* Get field value */
      field_ptr
        = get_field(first_field_ptr,"CHAIN",loop,params->instance);

      /* If valid, save the chain identifier */
      if (field_ptr != NULL && field_ptr->value[0] != '\0')
        chain = field_ptr->value[0];

      /* Get the key defining the current promotif option
         (eg BETA_TURNS) */
      strcpy(key,params->option);
      string_ptr = strstr(key,"_OPTION");
      if (string_ptr != NULL)
        string_ptr[0] = '\0';

      /* If reading in all the PROMOTIF data for this chain, set flag */
      if (!strcmp(key,"SUMMARY"))
        {
          get_all = TRUE;
          started = TRUE;
        }
    }
  else
    get_all = TRUE;

  /* Open the .dbt file */
  if ((file_ptr = fopen(file_name,"r")) !=  NULL)
    {
      /* Loop while reading in records from the file */
      while (fgets(input_line,LINELEN,file_ptr) != NULL &&
             all_done == FALSE)
        {
          /* Truncate the string to strip off the newline */
          string_chop(input_line);

          /* Assume line wanted */
          if (promotif == TRUE)
            wanted = FALSE;
          else
            wanted = TRUE;

          /* Check for end of data */
          if (!strncmp(input_line,"##END_DATA",10))
            {
              all_done = TRUE;
              wanted = FALSE;
            }

          /* Check for the START record */
          else if (!strncmp(input_line,"START:",6))
            {
              /* If have stored data from a previous motif, create
                 a recording storing the count */
              if (get_all == TRUE && nstored > 0)
                store_motif_count(first_field_ptr,motif_name,0,nstored);

              /* Get the name of these data */
              strcpy(motif_name,input_line+6);

              /* Initialise count of motif stored */
              nstored = 0;

              /* If not started, check for start */
              if (started == FALSE)
                {
                  /* Check whether this is the correct option */
                  if (promotif == TRUE && !strcmp(motif_name,key))
                    {
                      /* Set flag that data about to start */
                      started = TRUE;
                    }
                  else if (promotif == FALSE)
                    get_all = started = TRUE;
                }

              /* Check for next START */
              else if (get_all == FALSE)
                all_done = TRUE;

              /* Line not wanted for further processing */
              wanted = FALSE;
            }

          /* For PROMOTIF data file, determine whether this line belongs
             to the required chain */
          else if (promotif == TRUE && started == TRUE)
            {
              /* Check whether this is the field-names line */
              if (strncmp(input_line+1,"::",2))
                {
                  /* Extract all the field names */
                  nnames
                    = get_tokens(input_line,field_name,MAX_TOKEN,':');
                  have_field_names = TRUE;
                }

              /* Otherwise, get the chain identifier */
              else
                {
                  inchn = input_line[0];

                  /* If this is the right chain, line is wanted */
                  if (inchn == chain)
                    wanted = TRUE;
                }
            }

          /* If not Promotif, but have a motif-name, then get the field
             names */
          else if (promotif == FALSE && motif_name[0] != '\0' &&
                   nnames == 0)
            {
              /* Extract all the field names */
              nnames
                = get_tokens(input_line,field_name,MAX_TOKEN,':');
              have_field_names = TRUE;
              wanted = FALSE;
            }

          /* If line wanted, extract and store the fields */
          if (wanted == TRUE)
            {
              /* Increment count of lines read in */
              line++;
              nstored++;

              /* Initialise field counter and flag */
              finished = FALSE;
              item = 0;
              ipos = 0;

              /* Loop until all data fields have been extracted */
              while (finished == FALSE)
                {
                  /* Look for next double-colon */
                  string_ptr = strstr(input_line+ipos,"::");

                  /* If have double-colon, extract field */
                  if (string_ptr != NULL)
                    {
                      /* Transfer field and store */
                      len = string_ptr - (input_line+ipos);
                      strncpy(field,input_line+ipos,len);
                      field[len] = '\0';

                      /* Increment starting position */
                      ipos = ipos + len + 2;
                    }

                  /* Otherwise, transfer rest of line, and set flag that
                     we're done */
                  else
                    {
                      /* Transfer across */
                      strcpy(field,input_line+ipos);

                      /* Set flag that this is the last field on the
                         line */
                      finished = TRUE;
                    }

                  /* Increment item count */
                  item++;
                  jpos = 0;
                  subitem = 0;
                  have_semicolons = FALSE;
                  done = FALSE;

                  /* Remove any leading blanks */
                  fst_pos = -1;
                  field_len = strlen(field);
                  for (lpos = 0; lpos < field_len; lpos++)
                    {
                      /* Get this character */
                      ch = field[lpos];

                      /* If first non-blank, store its position */
                      if (fst_pos == -1 && ch != ' ')
                        fst_pos = lpos;
                    }

                  /* If have leading blanks, shift field */
                  if (fst_pos > 0)
                    {
                      strcpy(tmp,field+fst_pos);
                      strcpy(field,tmp);
                    }

                  /* Get the current field length */
                  field_len = string_truncate(field,strlen(field) + 1);

                  /* Loop until all sub-items have been processed */
                  while (done == FALSE)
                    {
                      /* Check for &nbsp; strings */
                      string_ptr = strstr(field+jpos,"&nbsp;");

                      /* Look for next semi-colon */
                      if (string_ptr == NULL)
                        string_ptr = strstr(field+jpos,";");
                      else
                        string_ptr = NULL;
                  
                      /* If have semi-colon, extract sub-item */
                      if (string_ptr != NULL)
                        {
                          /* Transfer field and store */
                          len = string_ptr - (field+jpos);
                          strncpy(subfield,field+jpos,len);
                          subfield[len] = '\0';

                          /* Set flag that we have sub-items */
                          have_semicolons = TRUE;

                          /* Increment starting position */
                          jpos = jpos + len + 1;

                          /* If this is beyond the end of the current
                             field, we're done */
                          if (jpos > field_len - 1)
                            done = TRUE;
                        }

                      /* Otherwise, transfer rest of field, and set flag
                         that we're done */
                      else
                        {
                          /* Transfer across */
                          strcpy(subfield,field+jpos);

                          /* Set flag that this is the last field on
                             the line */
                          done= TRUE;
                        }

                      /* Increment sub-item count */
                      subitem++;
                      ndata++;

                      /* Form subscripts to define this field */
                      if (have_semicolons == TRUE)
                        sprintf(subscripts,"[%d][%d]",nstored,subitem);
                      else
                        sprintf(subscripts,"[%d]",nstored);

                      /* Form field name */
                      if (item <= nnames)
                        sprintf(string,"%s%s %s",field_name[item - 1],
                                subscripts,subfield);
                      else
                        sprintf(string,"DBT_FIELD%d%s %s",item,
                                subscripts,subfield);

                      /* Create a field record */
                      field_ptr
                        = create_field_record(first_field_ptr,string,
                                              FALSE,FALSE,LOOP,-1);
                    }

                  /* If field contained any sub-items, then store count */
                  if (have_semicolons == TRUE)
                    {
                      /* Form field name */
                      sprintf(name,"SUB_FIELDS%d",item);

                      /* Create a field for number of sub-items */
                      field_ptr
                        = store_motif_count(first_field_ptr,name,line,
                                            subitem);
                    }
                }

              /* If have the required number of records, set flag that done */
              if (nrecords > -1 && nstored >= nrecords)
                all_done = TRUE;
            }

          /* If we have just picked up the field names, and have an offset,
             set file pointer to that to start reading in the data from the
             required file position */
          if (offset > -1 && have_field_names == TRUE)
            {
              /* If have file offset, then start read at that position */
              error = fseek(file_ptr,offset,SEEK_SET);

              /* Unset the offset */
              offset = -1;
            }
        }

      /* Close the file */
      fclose(file_ptr);

      /* If have stored data from a previous motif, create
         a record storing the count */
      if (motif_name[0] != '\0' && nstored > 0)
        store_motif_count(first_field_ptr,motif_name,0,nstored);

      /* If only read in a single set of data, store count */
      else if (line > 0)
        store_motif_count(first_field_ptr,"DBT_FIELDS",0,line);
    }

  /* Return number of data items stored */
  return(ndata);
}
/***********************************************************************

create_contents_flags  -  Create all the flags relating to the Contents
                          keywords

***********************************************************************/

void create_contents_flags(struct field *first_field_ptr,
                           char *string,struct parameters *params)
{
  int iheader, last_header;

  struct field *field_ptr;

  /* Initialise variables */
  last_header = -1;

  /* Loop through the Contents keywords and check for presence in
     the data */
  for (iheader = 0; iheader < NCONTENT_HEADERS; iheader++)
    {
      /* Check whether this item is present */
      field_ptr
        = get_field(first_field_ptr,CONTENTS_HEADER[iheader],LOOP,
                    params->instance);

      /* If present, set keyword */
      if (field_ptr != NULL)
        {
          /* Store position of last field present */
          last_header = iheader;

          /* If this is the PROTEIN header, then create a WIRPLOT field
             as these are no longer generated by pdbhtml */
          if (!strcmp(CONTENTS_HEADER[iheader],"PROTEIN"))
            {
              /* From the string containing keyword-value pair */
              strcpy(string,"WIRPLOT TRUE");

              /* Create a field record */
              field_ptr
                = create_field_record(first_field_ptr,string,FALSE,
                                      FALSE,LOOP,-1);
            }
        }
    }

  /* Loop through a second time to create all the flags */
  for (iheader = 0; iheader < NCONTENT_HEADERS; iheader++)
    {
      /* Copy current key */
      strcpy(string,CONTENTS_HEADER[iheader]);

      /* Append whether this is the last item */
      if (iheader < last_header)
        strcat(string,"_MID");
      else
        strcat(string,"_LAST");
      strcat(string," TRUE");

      /* Create a field record */
      field_ptr
        = create_field_record(first_field_ptr,string,FALSE,FALSE,LOOP,-1);
    }
}
/***********************************************************************

create_tab_flags  -  Create all the flags relating to the tabs at the
                     top of each html page

***********************************************************************/

void create_tab_flags(struct field *first_field_ptr,char *string,
                      struct parameters *params)
{
#define NFLAGS     3

  int iflag, iheader, ncols, ntabs;

  static int plus_cols[NFLAGS] = { 0, 1, 5 };

  struct field *field_ptr;

  /* Initialise variables */
  ntabs = 0;

  /* Loop through the header tab keywords and determine how many tab
     are needed */
  for (iheader = 0; iheader < NTAB_HEADERS; iheader++)
    {
      /* Check whether this item is present */
      field_ptr
        = get_field(first_field_ptr,TAB_HEADER[iheader],LOOP,
                    params->instance);

      /* If present, increment count */
      if (field_ptr != NULL)
        ntabs++;
    }

  /* Create flags containing required tab parameters (nos. of columns,
     etc) */
  ncols = 9 - NTAB_HEADERS + ntabs;

  /* Loop over the flags that need to be created */
  for (iflag = 0; iflag < NFLAGS; iflag++)
    {
      /* From the string containing keyword-value pair */
      if (plus_cols[iflag] > 0)
        sprintf(string,"NCOLS%d %d",plus_cols[iflag],
                ncols + plus_cols[iflag]);
      else
        sprintf(string,"NCOLS %d",ncols + plus_cols[iflag]);

      /* Create a field record */
      field_ptr
        = create_field_record(first_field_ptr,string,FALSE,FALSE,LOOP,-1);
    }
}
/***********************************************************************

create_loop_flags  -  Create all the endloop/midloop flags

***********************************************************************/

void create_loop_flags(struct field *first_field_ptr,char *string,
                       struct parameters *params)
{
  int level;

  struct field *field_ptr;

  /* Loop to create the mid-loop/end-loop fields */
  for (level = 0; level < MAX_LEVEL; level++)
    {
      /* Create mid-loop record for this loop, set to TRUE */
      sprintf(string,"MIDLOOP%d TRUE",level + 1);

      /* Create a field record */
      field_ptr
        = create_field_record(first_field_ptr,string,FALSE,FALSE,
                              LOOP,-1);

      /* Repeat for end-loop record, set to FALSE */
      sprintf(string,"ENDLOOP%d FALSE",level + 1);

      /* Create a field record */
      field_ptr
        = create_field_record(first_field_ptr,string,FALSE,FALSE,
                              LOOP,-1);
    }
}
/***********************************************************************

create_structure_flags  -  Create all the flags relating to the
                           Structure keywords

***********************************************************************/

void create_structure_flags(struct field *first_field_ptr,
                            char *string,struct parameters *params)
{
  int iheader, last_header;

  struct field *field_ptr;

  /* Initialise variables */
  last_header = -1;

  /* Loop through the Structure keywords and check for presence in
     the data */
  for (iheader = 0; iheader < NSTRUCTURE_HEADERS; iheader++)
    {
      /* Check whether this item is present */
      field_ptr
        = get_field(first_field_ptr,STRUCTURE_HEADER[iheader],LOOP,
                    params->instance);

      /* If present, set keyword */
      if (field_ptr != NULL)
        {
          /* Store position of last field present */
          last_header = iheader;
        }
    }

  /* Loop through a second time to create all the flags */
  for (iheader = 0; iheader < NSTRUCTURE_HEADERS; iheader++)
    {
      /* Copy current key */
      strcpy(string,STRUCTURE_HEADER[iheader]);

      /* Append whether this is the last item */
      if (iheader < last_header)
        strcat(string,"_MID");
      else
        strcat(string,"_LAST");
      strcat(string," TRUE");

      /* Create a field record */
      field_ptr
        = create_field_record(first_field_ptr,string,FALSE,FALSE,LOOP,-1);
    }
}
/***********************************************************************

create_call_flags  -  Create the flags defining the call options and
                      call values (ie program's input parameters)

***********************************************************************/

void create_call_flags(struct field *first_field_ptr,char *string,
                       struct parameters *params)
{
  char number_string[20];

  int i;

  struct field *field_ptr;

  /* Create a field record for the call option parameter */
  sprintf(string,"%s TRUE",params->option);

  /* Create a field record */
  field_ptr
    = create_field_record(first_field_ptr,string,FALSE,FALSE,LOOP,-1);

  /* From the string containing keyword-value pair */
  sprintf(string,"CALL_NO%d TRUE",params->call_no);

  /* Create a field record */
  field_ptr
    = create_field_record(first_field_ptr,string,FALSE,FALSE,LOOP,-1);

  /* For PDBsum, create a call value flag for the call value */
  if (params->enzymes == FALSE)
    {
      /* Get the call value */
      sprintf(string,"CALL_VAL%d",params->instance[0]);
      for (i = 1; i < MAX_LEVEL; i++)
        {
          if (params->instance[i] > 0)
            {
              sprintf(number_string,".%d",params->instance[i]);
              strcat(string,number_string);
            }
        }
      strcat(string," TRUE");

      /* Create a field record */
      field_ptr
        = create_field_record(first_field_ptr,string,FALSE,FALSE,
                              LOOP,-1);
    }

  /* For enzymes, create a call value flag for the EC number */
  else
    {
      sprintf(string,"CALL_VAL%d",params->ec_number[0]);
      for (i = 1; i < MAX_LEVEL; i++)
        {
          if (params->ec_number[i] > 0)
            sprintf(number_string,".%d",params->ec_number[i]);
          else
            sprintf(number_string,".-");
          strcat(string,number_string);
        }
      strcat(string," TRUE");

      /* Create a field record */
      field_ptr
        = create_field_record(first_field_ptr,string,FALSE,FALSE,
                              LOOP,-1);
    }

  /* Get the start number */
  sprintf(string,"START_NO%ld TRUE",params->start_no);

  /* Create a field record */
  field_ptr
    = create_field_record(first_field_ptr,string,FALSE,FALSE,LOOP,-1);
}
/***********************************************************************

create_colour_flags  -  Create the flags defining certain colours 
                        to be used in the HTML pages

***********************************************************************/

void create_colour_flags(struct field *first_field_ptr,char *string,
                         struct parameters *params)
{
#define NCOLFLAGS  2

  int icol;

  static char *colour[NCOLFLAGS] = { "#EFEFF7", "#FFFFFF" };
  struct field *field_ptr;

  /* Loop over all the colours to be created */
  for (icol = 0; icol < NCOLFLAGS; icol++)
    {
      /* From the string containing keyword-value pair */
      sprintf(string,"ALT_COLOUR%d %s",icol,colour[icol]);

      /* Create a field record */
      field_ptr
        = create_field_record(first_field_ptr,string,FALSE,FALSE,
                              LOOP,-1);
    }
}
/***********************************************************************

create_counters  -  Create the all-purpose counters and text fields that
                    can be used within the templates

***********************************************************************/

void create_counters(struct field *first_field_ptr,char *string,
                     struct parameters *params)
{
  int icount;

  struct field *field_ptr;

  /* Loop over all the available counters creating a record for each */
  for (icount = 0; icount < NCOUNTERS; icount++)
    {
      /* From the string containing keyword-value pair */
      sprintf(string,"COUNTER%d 0",icount);

      /* Create a field record */
      field_ptr
        = create_field_record(first_field_ptr,string,FALSE,FALSE,LOOP,-1);

      /* Repeat for the text string */
      sprintf(string,"TEXT%d /null/",icount);

      /* Create a field record */
      field_ptr
        = create_field_record(first_field_ptr,string,FALSE,FALSE,LOOP,-1);
    }
}
/***********************************************************************

define_specials  -  Define the special flags used in the templates to
                    draw the tramlines, plus others

***********************************************************************/

void define_specials(struct field *first_field_ptr,
                     struct parameters *params)
{
  char name[NAMLEN + 1];
  char varsite_results_dir[FILENAME_LEN];
  char ligsearch_results_dir[FILENAME_LEN];
  char pdbsum_results_dir[FILENAME_LEN];
  char profunc_results_dir[FILENAME_LEN], string[FILENAME_LEN + 1];
  char split_dir[3], tmp[FILENAME_LEN + 1];

  char *string_ptr;

  int i, len;

  struct field *field_ptr;

  /* Initialise variables */
  varsite_results_dir[0] = '\0';
  ligsearch_results_dir[0] = '\0';
  pdbsum_results_dir[0] = '\0';
  profunc_results_dir[0] = '\0';

  /* Check whether this is the fake entry */
  field_ptr = get_field(first_field_ptr,"FAKE_PDB_ENTRY",LOOP,
                        params->instance);
  if (field_ptr != NULL && !strcmp(field_ptr->value,params->pdb_code))
    params->fake_entry = TRUE;

  /* If this is a standard PDB file, create a LINKS field so that the
     links tab is displayed */
  if (params->pdb_type == STANDARD_PDB)
    {
      /* From the string containing keyword-value pair */
      strcpy(string,"STANDARD_PDB TRUE");

      /* Create a field record */
      field_ptr
        = create_field_record(first_field_ptr,string,FALSE,FALSE,LOOP,-1);

      /* From the string containing keyword-value pair */
      strcpy(string,"LINKS TRUE");

      /* Create a field record */
      if (params->proprietary == FALSE)
        field_ptr
          = create_field_record(first_field_ptr,string,FALSE,FALSE,LOOP,-1);

      /* If have a security code, add that */
      if (params->security_code[0] != '\0')
        {
          /* From the string containing keyword-value pair */
          strcpy(string,"LOCK &code=");
          strcat(string,params->security_code);

          /* Create a field record */
          field_ptr
            = create_field_record(first_field_ptr,string,FALSE,FALSE,LOOP,-1);
        }
    }

  /* If this is an uploaded or MCSG file, form the appropriate identifier
     and associated security code */
  else if (params->pdb_type != DRAT_DB && params->pdb_type != DRUGPORT &&
           params->pdb_type != SAS && params->pdb_type != ALPHA_FOLD)
    {
      /* From the string containing keyword-value pair */
      strcpy(string,"LOCK ");
      strcat(string,"&pdb_type=");
      if (params->pdb_type == UPLOAD_PDB)
        strcat(string,"UPLOAD");
      else if (params->pdb_type == PROFUNC)
        strcat(string,"PROFUNC");
      else if (params->pdb_type == LIGSEARCH)
        strcat(string,"LIGSEARCH");
      else if (params->pdb_type == VARSITE)
        strcat(string,"VARSITE");
      else if (params->pdb_type == LIGPLUS)
        strcat(string,"LIGPLUS");
      else
        strcat(string,"MCSG");
      strcat(string,"&code=");
      strcat(string,params->security_code);

      /* Create a field record */
      field_ptr
        = create_field_record(first_field_ptr,string,FALSE,FALSE,LOOP,-1);

      /* From the string containing keyword-value pair */
      strcpy(string,"PDB_TYPE ");
      if (params->pdb_type == UPLOAD_PDB)
        strcat(string,"UPLOAD");
      else if (params->pdb_type == PROFUNC)
        strcat(string,"PROFUNC");
      else if (params->pdb_type == LIGSEARCH)
        strcat(string,"LIGSEARCH");
      else if (params->pdb_type == VARSITE)
        strcat(string,"VARSITE");
      else if (params->pdb_type == LIGPLUS)
        strcat(string,"LIGPLUS");
      else
        strcat(string,"MCSG");

      /* Create a field record */
      field_ptr
        = create_field_record(first_field_ptr,string,FALSE,FALSE,LOOP,-1);
    }

  /* Create all the Contents flags */
  create_contents_flags(first_field_ptr,string,params);

  /* Create all the tab flags */
  create_tab_flags(first_field_ptr,string,params);

  /* Create all the loop flags */
  create_loop_flags(first_field_ptr,string,params);

  /* Create all the structure flags */
  create_structure_flags(first_field_ptr,string,params);

  /* Create all the call option and call value flags */
  create_call_flags(first_field_ptr,string,params);

  /* Create all the colour flags */
  create_colour_flags(first_field_ptr,string,params);

  /* Create all the all-purpose counters */
  create_counters(first_field_ptr,string,params);

  /* For enzymes, if this is the bottom level, create bottom-level
     flag */
  if (params->enzymes == TRUE && params->bottom_level == TRUE)
    {
      /* From the string containing keyword-value pair */
      sprintf(string,"BOTTOM_LEVEL TRUE");

      /* Create a field record */
      field_ptr
        = create_field_record(first_field_ptr,string,FALSE,FALSE,
                              LOOP,-1);
    }

  /* If have an error condition, create the PDB_CODE field */
  if (params->have_error == TRUE || params->highlights == TRUE ||
      params->ligsearch == TRUE || params->ligplus == TRUE)
    {
      /* If missing PDB code, then replace by blank */
      if (!strcmp(params->pdb_code,"....") ||
          !strcmp(params->pdb_code,"n/a") ||
          !strcmp(params->pdb_code,"none"))
        params->pdb_code[0] = '\0';

      /* From the string containing keyword-value pair */
      if (params->enzymes == TRUE)
        sprintf(string,"EC_CODE %s",params->ec_code);

      /* Otherwise, use PDB code */
      else
        sprintf(string,"PDB_CODE %s",params->pdb_code);

      /* Create a field record */
      if (params->pdb_code[0] != '\0')
        field_ptr
          = create_field_record(first_field_ptr,string,FALSE,FALSE,
                                LOOP,-1);
    }

  /* For VarSite, create a UNIPROT_ACC field and others, if present */
  if (params->varsite == TRUE)
    {
      /* Create a UniProt accession field */
      sprintf(string,"UNIPROT_ACC %s",params->uniprot_acc);
      field_ptr
        = create_field_record(first_field_ptr,string,FALSE,FALSE,
                              LOOP,-1);

      /* If have a PDB code, create that */
      if (params->pdb_code[0] != '\0')
        {
          sprintf(string,"PDB_CODE %s",params->pdb_code);
          field_ptr
            = create_field_record(first_field_ptr,string,FALSE,FALSE,
                                  LOOP,-1);
        }

      /* Repeat for the chain id */
      if (params->chain != '\0')
        {
          sprintf(string,"CHAIN %c",params->chain);
          field_ptr
            = create_field_record(first_field_ptr,string,FALSE,FALSE,
                                  LOOP,-1);
        }
    }

  /* Save the chain id */
  if (params->chain != '\0')
    {
      sprintf(string,"PCHAIN %c",params->chain);
      field_ptr
        = create_field_record(first_field_ptr,string,FALSE,FALSE,
                                  LOOP,-1);
    }

  /* If have a user identifier, store */
  if (params->user_id != &null)
    {
      /* From the string containing keyword-value pair */
      sprintf(string,"USER_ID %s",params->user_id);

      /* Create a field record */
      field_ptr
        = create_field_record(first_field_ptr,string,FALSE,FALSE,
                              LOOP,-1);
    }

  /* If have a security code, store */
  if (params->security_code[0] != '\0')
    {
      /* From the string containing keyword-value pair */
      sprintf(string,"CODE %s",params->security_code);

      /* Create a field record */
      field_ptr
        = create_field_record(first_field_ptr,string,FALSE,FALSE,
                              LOOP,-1);
    }

  /* If have a working directory, store */
  if (params->dir_name[0] != '\0')
    {
      /* From the string containing keyword-value pair */
      sprintf(string,"DIRECTORY_NAME %s",params->dir_name);

      /* Create a field record */
      field_ptr
        = create_field_record(first_field_ptr,string,FALSE,FALSE,
                              LOOP,-1);
    }

  /* For DrugPort, define directory name */
  else if (params->drugport == TRUE)
    {
      /* From the string containing keyword-value pair */
      sprintf(string,"DIRECTORY_NAME %s",params->subdir_real);

      /* Create a field record */
      field_ptr
        = create_field_record(first_field_ptr,string,FALSE,FALSE,
                              LOOP,-1);
    }

  /* Repeat for special working directory */
  if (params->work_dir[0] != '\0')
    {
      /* From the string containing keyword-value pair */
      sprintf(string,"WORK_DIR %s",params->work_dir);

      /* Create a field record */
      field_ptr
        = create_field_record(first_field_ptr,string,FALSE,FALSE,
                              LOOP,-1);
    }

  /* Repeat for ligand ranges */
  if (params->ligand_ranges != &null)
    {
      /* From the string containing keyword-value pair */
      sprintf(string,"LIGAND_RANGES %s",params->ligand_ranges);

      /* Create a field record */
      field_ptr
        = create_field_record(first_field_ptr,string,FALSE,FALSE,
                              LOOP,-1);
    }

  /* Loop over the extra parameters to store any non-blank ones */
  for (i = 0; i < NPAR; i++)
    {
      /* If have a parameter defined, store */
      if (params->parameter[i][0] != '\0')
        {
          /* From the string containing keyword-value pair */
          sprintf(string,"PARAM%d %s",i + 1,params->parameter[i]);

          /* Create a field record */
          field_ptr
            = create_field_record(first_field_ptr,string,FALSE,FALSE,
                                  LOOP,-1);
        }
    }

  /* If called from EBI PDBsum home page, create FROM_EBI flag */
  if (params->from_ebi == TRUE)
    {
      /* From the string containing keyword-value pair */
      sprintf(string,"FROM_EBI TRUE");

      /* Create a field record */
      field_ptr
        = create_field_record(first_field_ptr,string,FALSE,FALSE,
                              LOOP,-1);

      /* From the string containing keyword-value pair */
      sprintf(string,"EBI_PARAM &EBI=TRUE");

      /* Create a field record */
      field_ptr
        = create_field_record(first_field_ptr,string,FALSE,FALSE,
                              LOOP,-1);
    }

  /* For ProFunc, define names of directories where the PDBsum
     database.dat file and ProFunc results are */
  if (params->profunc == TRUE)
    {
      /* If all data in the temporary directory (ie old system), define
         names of both directories here */
      if (strlen(params->user_id) > 4)
        {
          /* Get the name of the temporary directory */
          field_ptr = get_field(first_field_ptr,"TMP_DIR",LOOP,
                                params->instance);
          if (field_ptr != NULL)
            {
              /* Form name of ProFunc results directory */
              strcpy(profunc_results_dir,field_ptr->value);
              strcat(profunc_results_dir,DIR_SEP);
              strcat(profunc_results_dir,"profunc");
              strcat(profunc_results_dir,params->user_id);
              strcat(profunc_results_dir,DIR_SEP);

              /* PDBsum results directory is the same */
              strcpy(pdbsum_results_dir,profunc_results_dir);
            }
        }

      /* Otherwise, PDBsum directory is in the uploads area, and ProFunc
         results are in the PDBsum-ProFunc directories */
      else
        {
          /* Get the name of the PDBSUM-PROFUNC directory */
          field_ptr = get_field(first_field_ptr,"PDBSUM_PROFUNC_DIR",LOOP,
                                params->instance);
          if (field_ptr != NULL)
            {
              /* Form name of ProFunc results directory */
              strcpy(profunc_results_dir,field_ptr->value);
              strcat(profunc_results_dir,DIR_SEP);
              strcat(profunc_results_dir,"profunc");
              strcat(profunc_results_dir,params->pdb_code);
              strcat(profunc_results_dir,DIR_SEP);
            }

          /* Get the name of the PDBsum uploads directory */
          field_ptr = get_field(first_field_ptr,"UPLOAD_PDBSUM_DIR",LOOP,
                                params->instance);
          if (field_ptr != NULL)
            {
              /* Form name of ProFunc results directory */
              strcpy(pdbsum_results_dir,field_ptr->value);
              strcat(pdbsum_results_dir,DIR_SEP);
              strcat(pdbsum_results_dir,params->pdb_code);
              strcat(pdbsum_results_dir,DIR_SEP);
            }
          }

      /* From the string containing keyword-value pair */
      sprintf(string,"PROFUNC_RESULTS_DIR %s",profunc_results_dir);

      /* Create a field record */
      field_ptr
        = create_field_record(first_field_ptr,string,FALSE,FALSE,
                              LOOP,-1);

      /* From the string containing keyword-value pair */
      sprintf(string,"PDBSUM_RESULTS_DIR %s",pdbsum_results_dir);

      /* Create a field record */
      field_ptr
        = create_field_record(first_field_ptr,string,FALSE,FALSE,
                              LOOP,-1);
    }

  /* For LigSearch, define names of the directory where the LigSearch
     results are */
  if (params->ligsearch == TRUE)
    {
      /* Get the name of the LigSearch directory */
      field_ptr = get_field(first_field_ptr,"LIGSEARCH_REAL_RESULTS",LOOP,
                            params->instance);
      if (field_ptr != NULL)
        {
          /* Get the length of the code */
          len = strlen(params->pdb_code);

          /* If code is 4 characters long, take it to be a
             PDB code */
          if (len == 4)
            {
              /* Get middle 2 digits of the PDB code */
              strncpy(split_dir,params->pdb_code+1,2);
              split_dir[2] = '\0';
            }

          /* Otherwise, use the last two characters */
          else
            {
              /* Get last 2 digits of the UniProt code */
              strncpy(split_dir,params->pdb_code + len - 2,2);
              split_dir[2] = '\0';
            }

          /* Form name of LigSearch results directory */
          strcpy(ligsearch_results_dir,field_ptr->value);
          strcat(ligsearch_results_dir,DIR_SEP);
          strcat(ligsearch_results_dir,split_dir);
          strcat(ligsearch_results_dir,DIR_SEP);
          strcat(ligsearch_results_dir,params->pdb_code);
          strcat(ligsearch_results_dir,DIR_SEP);
        }

      /* From the string containing keyword-value pair */
      sprintf(string,"LIGSEARCH_RESULTS_DIR %s",ligsearch_results_dir);

      /* Create a field record */
      field_ptr
        = create_field_record(first_field_ptr,string,FALSE,FALSE,
                              LOOP,-1);
    }

  /* For VarSite, define names of the directory where the VarSite
     results are */
  if (params->varsite == TRUE)
    {
      /* Get the name of the VarSite directory */
      field_ptr = get_field(first_field_ptr,"VARSITE_REAL_UNIPROT_DIR",LOOP,
                            params->instance);
      if (field_ptr != NULL)
        {
          /* Get the length of the code */
          len = strlen(params->uniprot_acc);

          /* Get last 2 digits of the UniProt code */
          strncpy(split_dir,params->uniprot_acc + len - 2,2);
          split_dir[2] = '\0';

          /* Form name of VarSite results directory */
          strcpy(varsite_results_dir,field_ptr->value);
          strcat(varsite_results_dir,DIR_SEP);
          strcat(varsite_results_dir,split_dir);
          strcat(varsite_results_dir,DIR_SEP);
          strcat(varsite_results_dir,params->uniprot_acc);
          strcat(varsite_results_dir,DIR_SEP);
        }

      /* From the string containing keyword-value pair */
      sprintf(string,"DDIR %s",varsite_results_dir);

      /* Create a field record */
      field_ptr
        = create_field_record(first_field_ptr,string,FALSE,FALSE,
                              LOOP,-1);
    }

  /* If this is a proprietary version of PDBsum, create PROPRIETARY
     flag and change URL of home page */
  if (params->proprietary == TRUE)
    {
      /* From the string containing keyword-value pair */
      sprintf(string,"PROPRIETARY TRUE");

      /* Create a field record */
      field_ptr
        = create_field_record(first_field_ptr,string,FALSE,FALSE,
                              LOOP,-1);

      /* Get URL of the PDBsum Proprietary home page */
      field_ptr = get_field(first_field_ptr,"PDBSUM_PROP_HOME",
                            LOOP,params->instance);

      /* If we have it, overwrite PDBsum home page URL */
      if (field_ptr != NULL)
        {
          /* Save the URL */
          strcpy(string,field_ptr->value);

          /* Find the PDBsum home page */
          field_ptr = get_field(first_field_ptr,"PDBSUM_HOME_PAGE",
                                LOOP,params->instance);

          /* If found, replace the value */
          if (field_ptr != NULL)
            {
              /* Clear the old value */
              if (field_ptr->value != &null)
                free(field_ptr->value);

              /* Store the new value */
              store_string(&(field_ptr->value),string);
            }
        }
    }

  /* Create a field defining the colour of the bar on the PDBsum page */
  if (params->fake_entry == TRUE)
    sprintf(string,"BAR_COLOUR %s",BAR_COLOUR[STANDARD_PDB]);
  else
    sprintf(string,"BAR_COLOUR %s",BAR_COLOUR[params->pdb_type]);

  /* Create a field record */
  field_ptr
    = create_field_record(first_field_ptr,string,FALSE,FALSE,LOOP,-1);

  /* Remove tail end and any underscores from option */
  strcpy(name,params->option);
  string_ptr = strstr(name,"_OPTION");
  if (string_ptr != NULL)
    string_ptr[0] = '\0';
  for (i = 0; i < strlen(name); i++)
    if (name[i] == '_')
      name[i] = ' ';

  /* Create field defining run-time option */
  sprintf(string,"CALL_OPTION %s",name);

  /* Create a field record */
  field_ptr
    = create_field_record(first_field_ptr,string,FALSE,FALSE,LOOP,-1);

  /* Store score value */
  strcpy(name,params->score);
  for (i = 0; i < strlen(name); i++)
    if (name[i] == '_')
      name[i] = ' ';

  /* Create field defining score passed from previous page */
  sprintf(string,"CALL_SCORE %s",name);

  /* Create a field record */
  field_ptr
    = create_field_record(first_field_ptr,string,FALSE,FALSE,LOOP,-1);

  /* Create field for storing read-status from database reads during
     interpretation of template */
  sprintf(string,"READ_OK TRUE");

  /* Create a field record */
  field_ptr
    = create_field_record(first_field_ptr,string,FALSE,FALSE,LOOP,-1);

  /* If program called with the -relinks flag in PDBsum1, create the
     RELINKS field */
  if (params->relinks == TRUE)
    {
      /* Create the relinks flag */
      sprintf(string,"RELINKS TRUE");
      field_ptr
        = create_field_record(first_field_ptr,string,FALSE,FALSE,
                              LOOP,-1);
    }
}
/***********************************************************************

check_for_brackets  -  Check whether field name contains the given type
                       of brackets. If it does, splice out the brackets
                       and return the value they contain

***********************************************************************/

int check_for_brackets(char *field_name,char *open_bracket,
                       char *close_bracket,char *ftype)
{
  int i, ipos, len, value;

  char number_string[20], replace_string[LINELEN + 1];

  char *string1_ptr, *string2_ptr;

  /* Initialise value */
  ftype[0] = '\0';
  value = 0;

  /* Check for presence of open and close brackets */
  string1_ptr = strstr(field_name,open_bracket);
  string2_ptr = strstr(field_name,close_bracket);

  /* If both present, then process */
  if (string1_ptr != NULL && string2_ptr != NULL &&
      string2_ptr > string1_ptr)
    {
      /* Have both brackets, so extract the number between them */
      ipos = string1_ptr - field_name + 1;
      len = string2_ptr - string1_ptr - 1;
      strncpy(number_string,field_name + ipos,len);
      number_string[len] = '\0';

      /* Check for special "pdb" format command */
      if (!strcmp(number_string,"pdb"))
        {
          strcpy(ftype,"pdb");
          value = 0;
        }

      /* Check for special "chn" format command */
      else if (!strcmp(number_string,"chn"))
        {
          strcpy(ftype,"chn");
          value = 0;
        }

      /* Check for special "bchn" format command */
      else if (!strcmp(number_string,"bchn"))
        {
          strcpy(ftype,"bchn");
          value = 0;
        }

      /* Check for special "empty" check command */
      else if (!strcmp(number_string,"empty"))
        {
          strcpy(ftype,"empty");
          value = 0;
        }

      /* Check for special "zchn" format command */
      else if (!strcmp(number_string,"zchn"))
        {
          strcpy(ftype,"zchn");
          value = 0;
        }

      /* Check for special "hx" format command */
      else if (!strcmp(number_string,"hx"))
        {
          strcpy(ftype,"hx");
          value = 0;
        }

      /* Check for special "_" format command */
      else if (!strcmp(number_string,"_"))
        {
          strcpy(ftype,"_");
          value = 0;
        }

      /* Check for special "-" format command to replace zeroes by
         hyphens */
      else if (!strcmp(number_string,"-"))
        {
          strcpy(ftype,"-");
          value = 0;
        }

      /* Check for special character-replace format */
      else if (!strncmp(number_string,"s/",2))
        {
          strcpy(ftype,number_string);
          value = 0;
        }

      /* Check for special string-truncate indicator */
      else if (!strncmp(number_string,"d/",2))
        {
          strcpy(ftype,number_string);
          value = 0;
        }

      /* Otherwise, interpret the format command */
      else
        {
          /* Check if first two characters are not numeric (eg a 'b', 'c',
             'e', 'i', 'f', 'm', 'p' or 't', or the operators eq, gt, etc) */
          for (i = 0; i < 2 && i < len; i++)
            {
              if (number_string[i] < '0' || number_string[i] > '9')
                {
                  /* Store the format type and blank out */
                  ftype[i] = number_string[i];
                  ftype[i + 1] = '\0';
                  number_string[i] = ' ';
                }
            }

          /* Get the value itself */
          value = atoi(number_string);
        }

      /* Splice out the bracketed number from field name */
      strncpy(replace_string,string1_ptr,len + 2);
      replace_string[len + 2] = '\0';

      /* Perform splice */
      perform_splice(field_name,replace_string,"");
    }

  /* Return the value, if any */
  return(value);
}
/***********************************************************************

extract_tag  -  Retrieve the field name enclosed within the tags
                as well as the tag-containing string that needs to
                be replaced

                Tag formats:

                <LNK%>ref_to_url</LNK%>
                <FLD%>field_name</FLD>
                <FLD%>field_name{b2}</FLD>  ... Beginning 2 characaters
                                                from field
                <FLD%>field_name{e2}</FLD>  ... End 2 characaters from
                                                field
                <FLD%>field_name(n)</FLD>   ... nth item from field
                <FLD%>field_name{f}</FLD>   ... format chemical formula
                <FLD%>field_name{i2}</FLD>  ... integer field of two
                                                digits (can also have
                                                i3, and i4)
                <FLD%>field_name{c12}</FLD> ... character field to be
                                                chopped into segments
                                                of maximum length 12
                                                (can be any length)
                <FLD%>field_name{t12}</FLD> ... as for c12, but text is
                                                type-set
                <FLD%>field_name{r2}</FLD> ...  format value with 2
                                                decimal places
                <FLD%>field_name{s12}</FLD> ... as for c12, but width
                                                is fixed
                <FLD%>field_name{hx}</FLD>  ... remove any HTML tags
                                                from field
                <FLD%>field_name{x}</FLD>   ... remove leading spaces
                                                and underscores
                <FLD%>field_name{-}</FLD>   ... replace zero by hyphen
                <FLD%>field_name{empty}</FLD> ... return TRUE if field
                                                is blank or '-', FALSE
                                                otherwise

***********************************************************************/

int extract_tag(char *tmp_line,char *field_name,char *replace_string,
                int *ndig,int *sub,char *stype,char *ftype)
{
  char *string1_ptr, *string2_ptr;

  int len, itag, ndigits, tag_len, subscript;
  int got_tag;

  /* Initialise variables */
  got_tag = FALSE;
  ndigits = 0;
  subscript = 0;
  ftype[0] = stype[0] = '\0';

  /* Loop over all the tag types */
  for (itag = 0; itag < NTAGS && got_tag == FALSE; itag++)
    {
      /* Search for first occurrence of start tag */
      string1_ptr = strstr(tmp_line,TAG[itag][START]);

      /* If valid, get end tag */
      if (string1_ptr != NULL)
        {
          /* Get end tag */
          string2_ptr = strstr(string1_ptr,TAG[itag][END]);

          /* If this also valid, get the field name and string
             to be replaced */
          if (string2_ptr != NULL)
            {
              /* Get the length of the start-tag */
              tag_len = strlen(TAG[itag][START]);

              /* Get the length of the field name */
              len = string2_ptr - (string1_ptr + tag_len);

              /* Extract the name of the field between the bounding
                 tags */
              strncpy(field_name,string1_ptr + tag_len,len);
              field_name[len] = '\0';

              /* Check for a number in brackets indicating that a
                 subscript is present */
              subscript = check_for_brackets(field_name,"(",")",stype);

              /* Check for a number in curly brackets defining the
                 formatting to be applied to this field */
              ndigits = check_for_brackets(field_name,"{","}",ftype);

              /* Get the length of the end-tag */
              tag_len = strlen(TAG[itag][END]);

              /* Get the length of the string to be replaced */
              len = string2_ptr + tag_len - string1_ptr;

              /* Extract the string to be replaced */
              strncpy(replace_string,string1_ptr,len);
              replace_string[len] = '\0';

              /* Set flag that we have a tag */
              got_tag = TRUE;
            }
        }
    }

  /* Return number of digits/characters and format type */
  *ndig = ndigits;

  /* Return the string subscript, if there is one */
  *sub = subscript;

  /* Return whether a tag has been correctly identified */
  return(got_tag);
}
/***********************************************************************

get_subvalue  -  Extract the appropriate substring from the field

***********************************************************************/

void get_subvalue(char *value,int subscript)
{
  char ch, tmp_value[LINELEN + 1];

  int field, field_start, ipos, len;
  int done;

  /* Initialise variables */
  done = FALSE;
  field = 1;
  field_start = 0;

  /* Copy full value into temporary string */
  strcpy(tmp_value,value);
  len = strlen(tmp_value);

  /* Loop through all the characters in the string, counting dot- or
     space-separated fields */
  for (ipos = 0; ipos < len + 1 && done == FALSE; ipos++)
    {
      /* Get this character */
      ch = tmp_value[ipos];

      /* If this is a field separator, increment field count */
      if (ch == ' ' || ch == '.' || ipos == len)
        {
          /* If current field is the one we're after, we can end here */
          if (field == subscript)
            {
              /* Terminate the string here */
              tmp_value[ipos] = '\0';

              /* Copy field from the field start */
              strcpy(value,tmp_value + field_start);

              /* Set flag that we're done */
              done = TRUE;
            }

          /* Otherwise, reset the field start */
          else
            field_start = ipos + 1;

          /* Increment field counter */
          field++;
        }
    } 
}
/***********************************************************************

link_to_pdb  -  Replace any PDB-code-like strings by full URLs

***********************************************************************/

void link_to_pdb(char *string)
{
  char ch, replace_string[6];
  char punct[2], value[LINELEN + 1];

  char *token[MAX_TOKEN];

  int end_pos, ipos, itoken, len, ntokens;
  int valid;

  /* Initialise variables */
  end_pos = 0;

  /* Add space to end of string */
  strcat(string," ");

  /* Retrieve all the tokens from the current line */
  ntokens = get_tokens(string,token,MAX_TOKEN,' ');

  /* Loop through the tokens to see if any resemble PDB codes */
  for (itoken = 0; itoken < ntokens; itoken++)
    {
      /* Get the length of this token */
      len = strlen(token[itoken]);
      strcpy(punct," ");

      /* If length is 5, check for a punctuation mark at the end */
      if (len == 5)
        {
          /* Get the fifth character */
          ch = token[itoken][4];

          /* Check whether this is a punctuation mark */
          if (ch == ',' || ch == '.' || ch == '?' || ch == '!' || 
              ch == ';' || ch == ':')
            {
              /* Shorten token, and store punctuation mark */
              len = 4;
              token[itoken][len] = '\0';
              punct[0] = ch;
              punct[1] = '\0';
            }
        }

      /* Get the first character */
      ch = token[itoken][0];

      /* If right length, check that first character is numeric */
      if (len == 4 && ch >= '0' && ch <= '9')
        {
          /* Initialise validation flag */
          valid = TRUE;

          /* Check the remaining characters */
          for (ipos = 1; ipos < 4 && valid == TRUE; ipos++)
            {
              /* Get this character */
              ch = token[itoken][ipos];

              /* Check that is a number or character */
              if ((ch >= '0' && ch <= '9') || (ch >= 'a' && ch <= 'z'))
                valid = TRUE;
              else
                valid = FALSE;
            }

          /* If a valid PDB code, then replace by full URL */
          if (valid == TRUE)
            {
              /* Look for occurrence of this code in the original
                 string */
              strcpy(replace_string,token[itoken]);
              strcat(replace_string,punct);

              /* Form the replacement string, including URL */
              strcpy(value,"<b><a href=\"<LNK%>GET_PDBSUM_ENTRY</LNK%>");
              strcat(value,token[itoken]);
              strcat(value,"\" class=menuClass>");
              strcat(value,token[itoken]);
              strcat(value,"</a></b>");
              strcat(value,punct);

              /* Perform the splice */
              end_pos
                = perform_splice(string+end_pos,replace_string,value);
            }
        }
    }
}
/***********************************************************************

link_to_ec  -  Replace any E.C.-code-like strings by full URLs

***********************************************************************/

void link_to_ec(char *text,struct field *first_field_ptr,
                struct parameters *params)
{
  char file_name[FILENAME_LEN + 1], ec_code[EC_LEN];
  char ch, string[TOKEN_LEN], tmp[TOKEN_LEN];
  char punct[2], number_string[20], value[LINELEN + 1];

  char *token[MAX_TOKEN];

  int ec_number[4], end_pos, ipos, itoken, len, ndigits, npos, ntokens;
  int have_file, ok;

  struct field *field_ptr;

  FILE *file_ptr;

  /* Add space to end of string */
  strcat(text," ");

  /* Retrieve all the tokens from the current line */
  ntokens = get_tokens(text,token,MAX_TOKEN,' ');

  /* Remove final space */
  len = strlen(text);
  text[len] = '\0';

  /* Initialise variables */
  end_pos = 0;

  /* Loop through the tokens to see if any resemble PDB codes */
  for (itoken = 0; itoken < ntokens; itoken++)
    {
      /* Get this token */
      strcpy(string,token[itoken]);

      /* Get its length */
      len = strlen(string);
      strcpy(punct," ");

      /* Check for last character being a punctuation mark */
      ch = string[len - 1];

      /* Check whether this is a punctuation mark */
      if (ch == ',' || ch == '.' || ch == '?' || ch == '!' || 
          ch == ';' || ch == ':')
        {
          /* Shorten token, and store punctuation mark */
          len = len - 1;
          string[len] = '\0';
          punct[0] = ch;
          punct[1] = '\0';
        }

      /* Check for EC or E.C. token at start */
      if (!strncmp(string,"EC",2))
        {
          strcpy(tmp,string+2);
          strcpy(string,tmp);
        }
      else if (!strncmp(string,"E.C.",4))
        {
          strcpy(tmp,string+4);
          strcpy(string,tmp);
        }
      len = strlen(string);

      /* Add dot to end of the string */
      string[len] = '.';
      string[len + 1] = '\0';
      len = len + 1;

      /* Initialise counts of dots and numbers */
      ok = TRUE;
      ndigits = 0;
      npos = 0;
      number_string[0] = '\0';

      /* Loop through the characters in this string to count dots and
         digits */
      for (ipos = 0; ipos < len && ok == TRUE; ipos++)
        {
          /* Get this character */
          ch = string[ipos];

          /* Check for dot ... */
          if (ch == '.')
            {
              /* If have a built-up digit, then interpret and store */
              if (npos > 0 && ndigits < 4)
                {
                  /* Terminate string and store */
                  number_string[npos] = '\0';
                  ec_number[ndigits] = atoi(number_string);

                  /* Increment count of digits in the E.C. number */
                  ndigits++;
                }

              /* If have more than 4 digits making up this number,
                 it can't be an E.C. number so skip it */
              else if (npos > 0 && ndigits >= 4)
                ok = FALSE;

              /* Reinitialise for next number */
              number_string[0] = '\0';
              npos = 0;
            }

          /* ... or for digit */
          else if (ch >= '0' && ch <= '9')
            {
              /* Store current digit */
              if (npos < 20)
                {
                  number_string[npos] = ch;
                  npos++;
                }
              else
                ok = FALSE;
            }

          /* Otherwise, this can't be an E.C. number */
          else
            ok = FALSE;
        }

      /* If valid E.C. number, form URL to link to its page in the
         Enzyme Structures Database */
      if (ok == TRUE && ndigits == 4)
        {
          /* Remove the final dot from the string */
          string[len - 1] = '\0';

          /* Initialised flag */
          have_file = FALSE;

          /* Get path to enzymes directory */
          field_ptr = get_field(first_field_ptr,"ENZYME_REAL_DIR",
                                LOOP,params->instance);

          /* Form name of this E.C. number's database.dat file and
             check whether it exists */
          if (field_ptr != NULL)
            {
              /* Form file name */
              sprintf(file_name,
                      "%s%sec%1d%sec%2.2d%sec%2.2d%sec%4.4d%sdatabase."
                      "dat",
                      field_ptr->value,DIR_SEP,ec_number[0],DIR_SEP,
                      ec_number[1],DIR_SEP,ec_number[2],DIR_SEP,
                      ec_number[3],DIR_SEP);

              /* Check whether the file exists */
              if ((file_ptr = fopen(file_name,"r")) !=  NULL)
                {
                  /* Set flag that have file */
                  have_file = TRUE;

                  /* Close the file */
                  fclose(file_ptr);
                }
            }

          /* If database.dat file present, form url and splice into
             original text string */
          if (have_file == TRUE)
            {
              /* Form standard form of E.C. number */
              sprintf(ec_code,"%d.%d.%d.%d",ec_number[0],ec_number[1],
                      ec_number[2],ec_number[3]);

              /* Get path to enzymes directory */
              field_ptr = get_field(first_field_ptr,"GET_ECPAGE_PL",
                                    LOOP,params->instance);

              /* Form name of this E.C. number's database.dat file and
                 check whether it exists */
              if (field_ptr != NULL)
                {
                  /* Form URL */
                  strcpy(value,"<b><A HREF=\"");
                  strcat(value,field_ptr->value);
                  strcat(value,ec_code);
                  strcat(value,"\" class=menuClass>");
                  strcat(value,string);
                  strcat(value,"</A></b>");

                  /* Perform the splice */
                  end_pos = perform_splice(text+end_pos,string,value);
                }
            }
        }
    }
}
/***********************************************************************

punc_chain  -  Check for specific punctuation marks in the chain
               identifier and convert for accessing files with the
               chain identifier in their name

***********************************************************************/

void punc_chain(char *string)
{
#define NPROBLEM_CHAINS       5

  int ipos;

  static char *PROBLEM_CHAIN = "#:<>|";
  static char *CONVERT_CHAIN[NPROBLEM_CHAINS] =
    {
      "_hash", "_colon", "_lt", "_gt", "_pipe"
    };

  /* Loop through the rogue chain identifiers */
  for (ipos = 0; ipos < NPROBLEM_CHAINS; ipos++)
    {
      /* Check for match */
      if (string[0] == PROBLEM_CHAIN[ipos])
        strcpy(string,CONVERT_CHAIN[ipos]);
    }
}
/***********************************************************************

html_format_formula  -  Format formula using HTML subscript and
                        superscript tags

***********************************************************************/

void html_format_formula(char *formula)
{
  char ch, tmp[TMP_STRING_LEN], tmp_string[TMP_STRING_LEN];

  char *string_ptr;

  int len, ipos, last_state, outpos, state;

  /* Initialise variables */
  outpos = 0;
  tmp_string[0] = '\0';
  last_state = BLANK;
  state = BLANK;

  /* Get length of input formula */
  len = strlen(formula);

  /* If formula is in the format N(formula), then remove the outer
     brackets */
  string_ptr = strstr(formula,"(");
  if (string_ptr != NULL)
    {
      /* If close-bracket is last character, then knock out */
      if (formula[len - 1] == ')')
        formula[len - 1] = '\0';

      /* Remove lead bracket */
      strcpy(formula,string_ptr + 1);

      /* Recalculate formula length */
      len = strlen(formula);
    }

  /* Loop through the characters in the formula */
  for (ipos = 0; ipos < len + 1; ipos++)
    {
      /* Get the next character position */
      ch = formula[ipos];

      /* Determine what we have */
      if ((ch >= 'A' && ch <= 'Z') || (ch >= 'a' && ch <= 'z'))
        state = ELEMENT;
      else if (ch >= '0' && ch <= '9' && state != CHARGE)
        state = NUMBER;
      else if (ch == '[')
        state = CHARGE;
      else if (ch == ']')
        state = BLANK;
      else if (ch == ' ' || ch == '\0' || state != CHARGE)
        state = BLANK;

      /* If have ended a number, close it off */
      if (last_state == NUMBER && state != NUMBER)
        {
          /* Close the subscript command */
          strcpy(tmp_string+outpos,"</SUB>");

          /* Adjust output position */
          outpos = outpos + 6;
        }

      /* If have ended a charge, close it off */
      else if (last_state == CHARGE && state != CHARGE)
        {
          /* Close the superscript command */
          strcpy(tmp_string+outpos,"</SUP>");

          /* Adjust output position */
          outpos = outpos + 6;
        }

      /* If have the continuation of an element, put into lower-case */
      else if (state == ELEMENT && last_state == ELEMENT)
        {
          /* Make character lower-case if need-be */
          if (ch >= 'A' && ch <= 'Z')
            ch = ch - 'A' + 'a';
        }

      /* If starting a number, then set subscript on */
      else if (state == NUMBER && last_state != NUMBER)
        {
          /* Start the subscript command */
          strcpy(tmp_string+outpos,"<SUB>");

          /* Adjust output position */
          outpos = outpos + 5;
        }

      /* If starting a charge, then set superscript on */
      else if (state == CHARGE && last_state != CHARGE)
        {
          /* Start the superscript command */
          strcpy(tmp_string+outpos,"<SUP>");

          /* Adjust output position */
          outpos = outpos + 5;
        }

      /* If character to be transferred, then do so */
      if (state == ELEMENT || state == NUMBER || state == CHARGE)
        {
          /* Transfer across to the output line */
          tmp_string[outpos] = ch;

          /* Increment output position */
          outpos++;
          tmp_string[outpos] = '\0';
        }

      /* Save the current state */
      last_state = state;

      /* Check whether output line getting too long to continue */
      if (outpos + 6 + len - ipos > TMP_STRING_LEN - 1)
        {
          /* Make sure loop ends here */
          ipos = len + 2;

          /* Use the formula just as it stands */
          strcpy(tmp_string,formula);
        }
    }

  /* Post-process the string to remove any singletons in the formula */
  string_ptr = strstr(tmp_string,"<SUB>1</SUB>");

  /* Loop until all the 1's have been removed */
  while (string_ptr != NULL)
    {
      /* Excise the offending string */
      strcpy(tmp,string_ptr+12);
      strcpy(string_ptr,tmp);

      /* Check for any more */
      string_ptr = strstr(tmp_string,"<SUB>1</SUB>");
    }

  /* Return the formatted line */
  strcpy(formula,tmp_string);
}
/***********************************************************************

split_text  -  Split up the given character string into smaller chunks

***********************************************************************/

int split_text(char *string,int segment_len,int fixed)
{
  char ch;
  char out_string[LINELEN + 1], tmp_string[LINELEN + 1];

  char *string_ptr;

  static char vowels[] = { "aeiouAEIOU" };
  
  int ipos, last_consonant, last_hyphen, last_pos;
  int nleft, len, line_len, nchunks;
  int end_pos, start_pos;
  int done, got_line;

  /* Initialise variables */
  done = FALSE;
  nchunks = 0;
  ipos = 0;
  len = strlen(string);
  out_string[0] = '\0';

  /* Loop until list has been split up */
  while (done == FALSE)
    {
      /* Determine how much of the line is left to do */
      nleft = len - ipos;

      /* If this is smaller than the segment size, transfer across
         and we are done */
      if (nleft < segment_len)
        {
          /* Check that some part of the string still left */
          if (nleft > 0)
            {
              /* If this is not the first segment, add a break after
                 previous one */
              if (out_string[0] != '\0')
                strcat(out_string,"<BR>");

              /* Append last part of line */
              strcat(out_string,string + ipos);

              /* Increment line count */
              nchunks++;
            }

          /* Set flag that we're done */
          done = TRUE;
        }

      /* Otherwise, need to cut out a segment-sized chunk */
      else
        {
          /* Initialise position of last hyphen and consonant */
          last_hyphen = -1;
          last_consonant = -1;
          last_pos = -1;
          got_line = FALSE;
          line_len = 0;
          start_pos = ipos;

          /* Scan string until required line length obtained */
          while (got_line == FALSE)
            {
              /* Get current character */
              ch = string[ipos];
              line_len = ipos - start_pos + 1;

              /* If segment length exceeded, then terminate at last
                 hyphen or consonant */
              if (line_len > segment_len || ipos == len)
                {
                  /* If this is the last character in the string,
                     then end here */
                  if (ipos == len)
                    end_pos = ipos;

                  /* If hyphen encountered, then end at hyphen */
                  else if (last_hyphen > 0)
                    end_pos = last_hyphen;

                  /* If consonant encountered, then end there */
                  else if (last_consonant > 0 && fixed == FALSE)
                    end_pos = last_consonant;

                  /* Otherwise, end here */
                  else
                    end_pos = ipos;
                  
                  /* If this is not the first segment, add a break
                     after previous one */
                  if (out_string[0] != '\0')
                    strcat(out_string,"<BR>");

                  /* Get segment length */
                  line_len = end_pos - start_pos + 1;

                  /* Append current segment */
                  strncpy(tmp_string,string + start_pos,line_len);
                  tmp_string[line_len] = '\0';
                  strcat(out_string,tmp_string);

                  /* Set flag that we have a line */
                  got_line = TRUE;

                  /* Increment line count */
                  nchunks++;

                  /* Reset starting position */
                  ipos = start_pos + line_len;
                  start_pos = ipos;
                }

              /* Otherwise, copy character across */
              else
                {
                  /* If this is a hyphen, store its position */
                  if (ch == '-')
                    last_hyphen = ipos;

                  /* If punctuation mark, then is a good break
                     point */
                  else if (ch == '.' || ch == ',' || ch == ';' ||
                           ch == ':' || ch == ' ')
                    last_consonant = ipos;

                  /* If this is a vowel, store position of last
                     consonant */
                  string_ptr = strchr(vowels,ch);
                  if (string_ptr != NULL)
                    last_consonant = last_pos;
                  else
                    last_pos = ipos;

                  /* Increment character position and line length */
                  ipos++;
                }
            }
        }
    }

  /* Return the chopped-up string */
  strcpy(string,out_string);

  /* Return number of segments */
  return(nchunks);
}
/***********************************************************************

store_counter  -  Store the value of the given counter in the
                  corresponding counter field

***********************************************************************/

void store_counter(struct field *first_field_ptr,int icount,
                   long ival,struct parameters *params)
{
  char field_name[NAMLEN + 1], value[10];

  struct field *field_ptr, *matched_field_ptr;

  /* Form name of field to search for */
  sprintf(field_name,"COUNTER%d",icount);

  /* Locate the counter field */
  field_ptr = get_field(first_field_ptr,field_name,LOOP,params->instance);

  /* If have a match, update the field's value */
  if (field_ptr != NULL)
    {
      /* Format value */
      sprintf(value,"%ld",ival);

      /* Update the field's value */
      if (field_ptr->value != &null)
        free(field_ptr->value);
      store_string(&(field_ptr->value),value);
    }
}
/***********************************************************************

store_text  -  Store the text value in the given text field

***********************************************************************/

void store_text(struct field *first_field_ptr,int icount,char *string,
                struct parameters *params)
{
  char field_name[NAMLEN + 1], value[10];

  struct field *field_ptr, *matched_field_ptr;

  /* Form name of field to search for */
  sprintf(field_name,"TEXT%d",icount);

  /* Locate the counter field */
  field_ptr = get_field(first_field_ptr,field_name,LOOP,params->instance);

  /* If have a match, update the field's value */
  if (field_ptr != NULL)
    {
      /* Update the field's value */
      if (field_ptr->value != &null)
        free(field_ptr->value);
      if (!strcmp(string,"/null/"))
        field_ptr->value = &null;
      else
        store_string(&(field_ptr->value),string);
    }
}
/***********************************************************************

apply_template  -  Apply the PDBsum template using given PDB code

***********************************************************************/

void apply_template(char *string,struct field *first_field_ptr,
                    char *value,struct parameters *params)
{
  char pdb_code[5];

  char *string_ptr, *tmp;
  
  struct field *field_ptr;

  /* Initialise the output value */
  value[0] = '\0';

  /* Extract the PDB code from the string */
  string_ptr = strchr(string,'{');
  if (string_ptr == NULL)
    string_ptr = strchr(string,'[');
  if (string_ptr == NULL)
    string_ptr = strchr(string,'(');
  if (string_ptr != NULL)
    {
      strncpy(pdb_code,string_ptr + 1,4);
      pdb_code[4] = '\0';

      /* Get the directory template */
      field_ptr
        = get_field(first_field_ptr,"PDBSUM_REAL_TEMPLATE",LOOP,
                    params->instance);

      /* If found, create a the path to the PDBsum directory */
      if (field_ptr != NULL)
        {
          /* Form the full directory name */
          tmp = form_filename(pdb_code,field_ptr->value);
          strcpy(value,tmp);
        }
    }
}
/***********************************************************************

annotate_line  -  Replace all tags by appropriate substitutions

***********************************************************************/

int annotate_line(char *tmp_line,struct field *first_field_ptr,
                  int *loop,int *cno,int *onf,long *counter,
                  struct parameters *params)
{
  char field_name[NAMLEN + 1];
  char ftype[60], stype[60];
  char new_string[LINELEN + 1], old_string[LINELEN + 1];
  char replace_string[LINELEN + 1];
  char value[LINELEN + 1], number_string[20], tmp[LINELEN + 1];

  char *string_ptr;

  int call_no, i, icount, ivalue, len, level, ndigits, nloop, number;
  int end_pos, fpos, on_off, subscript;
  int done, finished, fixed, got_tag, splice;

  double rvalue;

  struct field *field_ptr;

  /* Initialise variables */
  call_no = *cno;
  on_off = *onf;
  done = FALSE;
  nloop = 0;
  splice = FALSE;

  /* Loop until no more tags to be replaced */
  while (done == FALSE)
    {
      /* Assume we're done */
      done = TRUE;

      /* Check for any tag ends still in the line */
      if (strstr(tmp_line,TAG_END) != NULL && nloop < 100)
        {
          /* Increment loop count to ensure we don't hang in the case
             of an error in the template */
          nloop++;

          /* Extract the tag, field name and string to be
             replaced */
          got_tag = extract_tag(tmp_line,field_name,replace_string,
                                &ndigits,&subscript,stype,ftype);

          /* If have a valid tag, get the replacement string and
             perform the splice */
          if (got_tag == TRUE)
            {
              /* Initialise replace value */
              value[0] = '\0';

              /* Special case of PSUM directory */
              if (!strncmp(field_name,"PSUM[",5))
                {
                  apply_template(replace_string,first_field_ptr,
                                 value,params);
                  field_ptr = NULL;
                }

              /* Get field value */
              else
                field_ptr = get_field(first_field_ptr,field_name,loop,
                                      params->instance);

              /* If valid, perform the splice */
              if (field_ptr != NULL)
                {
                  /* Get this field value */
                  strcpy(value,field_ptr->value);

                  /* If looking for a specific character position, extract
                     just that character */
                  if (!strcmp(stype,"p"))
                    {
                      /* Get the length of the current field */
                      len = strlen(value);

                      /* If position within length of the field, get the
                         required character */
                      if (subscript - 1 < len)
                        {
                          value[0] = value[subscript - 1];
                          value[1] = '\0';
                        }

                      /* Otherwise, blank out */
                      else
                        value[0] = '\0';
                    }

                  /* If looking for substring of given length, extract that */
                  else if (!strcmp(stype,"s"))
                    {
                      /* Get the length of the current field */
                      len = strlen(value);

                      /* If length less than length of field, get the
                         substring */
                      if (subscript <= len)
                        value[subscript] = '\0';

                      /* Otherwise, blank out */
                      else
                        value[0] = '\0';
                    }

                  /* If looking for a subscript, extract the
                     appropriate substring from the field */
                  else if (subscript > 0)
                    get_subvalue(value,subscript);

                  /* Check for "pdb" format instruction which requires
                     replacement of any PDB-code-like strings by full
                     URLs */
                  if (!strcmp(ftype,"pdb"))
                    link_to_pdb(value);

                  /* Check for E.C. codes in the string and add URLs
                     to Enzyme Structures Database */
                  else if (!strcmp(ftype,"ec"))
                    link_to_ec(value,first_field_ptr,params);

                  /* If have a chain identifier, check for punctuation
                     markes that need converting */
                  else if (!strcmp(ftype,"chn"))
                    punc_chain(value);

                  /* If have a chain identifier, replace blank chain
                     by underscore */
                  else if (!strcmp(ftype,"bchn") &&
                           (value[0] == ' ' || value[0] == '\0'))
                    {
                      strcpy(value,"_");
                    }

                  /* If have a chain identifier, replace blank chain
                     by zero */
                  else if (!strcmp(ftype,"zchn") &&
                           (value[0] == ' ' || value[0] == '\0'))
                    {
                      strcpy(value,"0");
                    }

                  /* If have a field with spaces, replace by underscores */
                  else if (!strcmp(ftype,"_"))
                    {
                      for (i = 0; i < strlen(value); i++)
                        if (value[i] == ' ')
                          value[i] = '_';
                    }

                  /* Extract first ndigits characters from the string */
                  else if (!strcmp(ftype,"b"))
                    {
                      /* Get the length of the field */
                      len = strlen(value);

                      /* Extract the substring */
                      if (ndigits < len)
                        {
                          strncpy(tmp,value,ndigits);
                          value[ndigits] = '\0';
                        }
                    }

                  /* Extract last ndigits characters from the string */
                  else if (!strcmp(ftype,"e"))
                    {
                      /* Get the length of the field */
                      len = strlen(value);

                      /* Determine start-pos of required substring */
                      if (len > ndigits)
                        fpos = len - ndigits;
                      else
                        fpos = 0;

                      /* Extract the substring */
                      if (fpos > 0)
                        {
                          strcpy(tmp,value+fpos);
                          strcpy(value,tmp);
                        }
                    }

                  /* Extract midle two characters from the string (for
                     a PDB code's parent directory) */
                  else if (!strcmp(ftype,"m"))
                    {
                      /* Get the length of the field */
                      len = strlen(value);

                      /* Determine start-pos of required substring */
                      if (len == 4)
                        fpos = 1;
                      else
                        fpos = 0;

                      /* Extract the substring */
                      if (fpos > 0)
                        {
                          strncpy(tmp,value + fpos,2);
                          tmp[2] = '\0';
                          strcpy(value,tmp);
                        }
                    }

                  /* For partial E.C. number, strip off trailing dots
                     and hyphens */
                  else if (!strcmp(ftype,"E"))
                    {
                      /* Check for first ".-" */
                      string_ptr = strstr(value,".-");

                      /* If found, then strip off */
                      if (string_ptr != NULL)
                        string_ptr[0] = '\0';
                    }

                  /* Remove leading spaces or underscores */
                  else if (!strcmp(ftype,"x"))
                    {
                      /* Initialise starting position */
                      fpos = -1;

                      /* Find first non-blank character */
                      for (i = 0; i < strlen(value) && fpos == -1; i++)
                        if (value[i] != ' ' && value[i] != '_')
                          fpos = i;

                      /* If have leading spaces, shift string along */
                      if (fpos > 0)
                        {
                          strcpy(tmp,value+fpos);
                          strcpy(value,tmp);
                        }
                    }

                  /* Remove any HTML tags in the string */
                  else if (!strcmp(ftype,"hx"))
                    {
                      /* Initialise flag */
                      finished = FALSE;
                      new_string[0] = '\0';

                      /* Loop until all replacements have been made */
                      while (finished == FALSE)
                        {
                          /* Search for first open bracket */
                          string_ptr = strchr(value,'<');

                          /* If not found, then we are done */
                          if (string_ptr == NULL)
                            finished = TRUE;

                          /* Otherwise, find end of string to be removed */
                          else
                            {
                              /* Copy to replace string */
                              strcpy(old_string,string_ptr);

                              /* Find the close-bracket */
                              string_ptr = strchr(old_string,'>');

                              /* If found, perform the splice */
                              if (string_ptr != NULL)
                                {
                                  /* Terminate the string to be removed */
                                  string_ptr[1] = '\0';

                                  /* Perform splice */
                                  perform_splice(value,old_string,
                                                 new_string);
                                }

                              /* If not found, then leave as is */
                              else
                                finished = TRUE;
                            }
                        }
                    }

                  /* If have a substitution command, effect the change */
                  else if (!strncmp(ftype,"s/",2))
                    {
                      /* Initialise both strings */
                      old_string[0] = new_string[0] = '\0';

                      /* Determine the substring to be replaced */
                      strcpy(old_string,ftype + 2);
                      string_ptr = strchr(old_string,'/');
                      if (string_ptr != NULL)
                        string_ptr[0] = '\0';

                      /* Get the replacement string */
                      string_ptr = strchr(ftype + 2,'/');
                      if (string_ptr != NULL)
                        {
                          strcpy(new_string,string_ptr + 1);
                          string_ptr = strchr(new_string,'/');
                          if (string_ptr != NULL)
                            string_ptr[0] = '\0';
                        }

                      /* Initialise flag */
                      finished = FALSE;
                      if (old_string[0] == '\0')
                        finished = TRUE;

                      /* Loop until all replacements have been made */
                      while (finished == FALSE)
                        {
                          /* Perform splice */
                          end_pos
                            = perform_splice(value,old_string,new_string);

                          /* If no replacements made, then we're done */
                          if (end_pos == 0)
                            finished = TRUE;
                        }
                    }

                  /* If have a string-truncation command, chop end off
                     string */
                  else if (!strncmp(ftype,"d/",2) && strlen(ftype) > 2)
                    {
                      /* If chop-string is present, then truncate the
                         string */
                      if ((string_ptr = strstr(value,ftype + 2)) != NULL)
                        string_ptr[0] = '\0';
                    }

                  /* Check whether field is a text field that requires
                     type setting */
                  else if (!strcmp(ftype,"t"))
                    {
                      type_set(value,TRUE);

                      /* Check for certain strings which should be in
                         upper-case */
                      fix_upper(value);
                    }

                  /* Check for text field to be set in lower case */
                  else if (!strcmp(ftype,"l"))
                    to_lower(value,strlen(value));

                  /* Check for text field to be set in upper case */
                  else if (!strcmp(ftype,"u"))
                    to_upper(value,strlen(value));

                  /* Check for formula field that needs formatting */
                  else if (!strcmp(ftype,"f"))
                    html_format_formula(value);

                  /* If need to chop up a character string, then
                     chop into required lengths */
                  if (ndigits > 0 && (!strcmp(ftype,"c") ||
                                      !strcmp(ftype,"s") ||
                                      !strcmp(ftype,"t")))
                    {
                      /* Determine whether splitting fixed, or variable */
                      if (!strcmp(ftype,"s"))
                        fixed = TRUE;
                      else
                        fixed = FALSE;

                      /* Chop the character string up into seqments */
                      split_text(value,ndigits,fixed);
                    }

                  /* If a numeric value, format according to number
                     of required digits */
                  else if (ndigits == 0 && !strcmp(ftype,"r"))
                    {
                      rvalue = atof(value);
                      sprintf(value,"%.0f",rvalue);
                    }

                  /* If a replacing zeroes by hyphens, check value */
                  else if (!strcmp(ftype,"-"))
                    {
                      ivalue = atoi(value);
                      if (ivalue == 0)
                        strcpy(value,"-");
                    }

                  /* If a numeric value, format according to number
                     of required digits */
                  else if (ndigits == 2)
                    {
                      if (!strcmp(ftype,"r"))
                        {
                          rvalue = atof(value);
                          sprintf(value,"%.2f",rvalue);
                        }
                      else if (!strcmp(ftype,"i"))
                        {
                          ivalue = atoi(value);
                          sprintf(value,"%2.2d",ivalue);
                        }
                    }
                  else if (ndigits == 3)
                    {
                      ivalue = atoi(value);
                      sprintf(value,"%3.3d",ivalue);
                    }
                  else if (ndigits == 4)
                    {
                      ivalue = atoi(value);
                      sprintf(value,"%4.4d",ivalue);
                    }

                  /* Set flag that a valid splice made in current line */
                  splice = TRUE;
                }

              /* Check for special loop tag which is to be replaced by
                 a loop counter */
              else if (!strncmp(field_name,"_LOOP",5))
                {
                  /* Extract which level we need */
                  strcpy(number_string,field_name+5);
                  level = atoi(number_string);

                  /* If a valid level, store the value */
                  if (level > 0 && level < MAX_LEVEL + 1)
                    if (loop[level - 1] > 0)
                      number = loop[level - 1];
                    else
                      number = params->instance[level - 1];
                  else
                    number = 1;

                  /* If a numeric value, format according to number
                     of required digits */
                  if (ndigits == 2)
                    sprintf(value,"%2.2d",number);
                  else if (ndigits == 3)
                    sprintf(value,"%3.3d",number);
                  else if (ndigits == 4)
                    sprintf(value,"%4.4d",number);

                  /* Otherwise format in free form */
                  else
                    sprintf(value,"%d",number);
                }

              /* Check for special loop tag which is to be replaced by
                 call counter */
              else if (!strncmp(field_name,"_CALL_NO",8))
                {
                  /* Store the call count */
                  sprintf(value,"%d",call_no);

                  /* Increment call counter */
                  if (!strcmp(field_name+8,"++"))
                    call_no++;
                }

              /* Check for special tag which is to be replaced by
                 start counter */
              else if (!strncmp(field_name,"_START_NO",9))
                {
                  /* Store the start count */
                  sprintf(value,"%ld",params->start_no);
                }

              /* Check for special tag which is to be replaced by
                 current option */
              else if (!strcmp(field_name,"_OPTION"))
                {
                  /* Store the option */
                  sprintf(value,"%s",params->option);

                  /* Strip off the "_OPTION" suffix */
                  string_ptr = strstr(value,"_OPTION");
                  if (string_ptr != NULL)
                    string_ptr[0] = '\0';
                }

              /* Check for special counter tag which is to be replaced by
                 relevant count */
              else if (!strncmp(field_name,"_COUNTER",8))
                {
                  /* Determine which counter this is */
                  number_string[0] = field_name[8];
                  number_string[1] = '\0';
                  icount = atoi(number_string);

                  /* If valid, store the current counter value */
                  if (icount >= 0 && icount < NCOUNTERS)
                    {
                      /* Check whether any formatting of the value is
                         required */
                      if (!strcmp(ftype,"i") && ndigits == 2)
                        sprintf(value,"%2.2ld",counter[icount]);

                      /* Store the counter value */
                      else
                        sprintf(value,"%ld",counter[icount]);

                      /* Increment call counter if required */
                      if (!strncmp(field_name+9,"+",1) ||
                          !strncmp(field_name+9,"-",1) ||
                          !strncmp(field_name+9,"=",1) ||
                          !strncmp(field_name+9,"*",1) ||
                          !strncmp(field_name+9,"/",1) ||
                          !strncmp(field_name+9,"%",1))
                        {
                          /* Increment counter */
                          if (!strcmp(field_name+9,"++"))
                            counter[icount]++;
                          else if (!strcmp(field_name+9,"--"))
                            counter[icount]--;
                          else if (!strncmp(field_name+9,"+",1))
                            counter[icount]
                              = counter[icount] + atol(field_name+10);
                          else if (!strncmp(field_name+9,"-",1))
                            counter[icount]
                              = counter[icount] - atol(field_name+10);
                          else if (!strncmp(field_name+9,"*",1))
                            counter[icount]
                              = counter[icount] * atol(field_name+10);
                          else if (!strncmp(field_name+9,"/",1))
                            {
                              if (atol(field_name+10) != 0)
                                counter[icount]
                                  = counter[icount] / atol(field_name+10);
                            }
                          else if (!strncmp(field_name+9,"%",1))
                            counter[icount]
                              = counter[icount] % atol(field_name+10);
                          else
                            counter[icount] = atol(field_name+10);

                          /* Store in the corresponding COUNTER field */
                          store_counter(first_field_ptr,icount,
                                        counter[icount],params);

                          /* If merely incrementing, then don't want to
                             display value */
                          value[0] = '\0';
                        }
                    }
                }

              /* Check for special text tag which is to be replaced by
                 relevant string */
              else if (!strncmp(field_name,"_TEXT",5))
                {
                  /* Determine which counter this is */
                  number_string[0] = field_name[5];
                  number_string[1] = '\0';
                  icount = atoi(number_string);

                  /* If valid, store the current text string */
                  if (icount >= 0 && icount < NCOUNTERS)
                    {
                      /* If text value being set, then store the
                         string */
                      if (!strncmp(field_name+6,"=",1))
                        {
                          strcpy(tmp,field_name+7);

                          /* Store in the corresponding TEXT field */
                          store_text(first_field_ptr,icount,tmp,params);

                          /* If storing text value, then don't want to
                             display it here */
                          value[0] = '\0';
                        }

                      /* Otherwise, display the given string */
                      else
                        {
                          /* Form name of field to search for */
                          sprintf(field_name,"TEXT%d",icount);

                          /* Locate the text field */
                          field_ptr
                            = get_field(first_field_ptr,field_name,LOOP,
                                        params->instance);

                          /* If have a match, update the field's value */
                          if (field_ptr != NULL)
                            strcpy(value,field_ptr->value);
                        }
                    }
                }

              /* Check for special on/off tag to be replaced by 0 or 1 */
              else if (!strncmp(field_name,"_ON_OFF",7))
                {
                  /* Store the on-off value */
                  sprintf(value,"%d",on_off);

                  /* Change on-off value if required */
                  if (!strcmp(field_name+7,"++"))
                    {
                      on_off++;
                      if (on_off > 1)
                        on_off = 0;
                    }
                }


              /* Check for zchn format for chain-id, making chain zero */
              else if (!strcmp(ftype,"zchn"))
                strcpy(value,"0");

              /* Perform splice */
              perform_splice(tmp_line,replace_string,value);

              /* Unset flag as there may be more tags left to find */
              done = FALSE;
            }
        }
    }

  /* Return current call counter and on-off value */
  *cno = call_no;
  *onf = on_off;

  /* Return whether a valid splice made in current line */
  return(splice);
}
/***********************************************************************

process_blocks  -  Process the current line, checking for block starts
                   and ends

***********************************************************************/

int process_blocks(char *input_line,int skip,int *instance,int *loop,
                   struct field *first_field_ptr,char **html_line,
                   int *nlines,char *skip_end,struct parameters *params)
{
  char field_name[LINELEN + 1], value[TMP_STRING_LEN + 1];
  char tmp[LINELEN + 1];
  char key[2], string[LINELEN + 1];

  char *string_ptr;

  char replace_from[2] = { "" };
  char replace_by[2] = { "" };

  int call_no, icount, on_off;
  int not, valid, wanted;

  long counter[NCOUNTERS];

  static const int start_record = 0;
  static const int nrecords = -1;

  static const long offset = -1;

  struct field *field_ptr;

  /* Initialise variables */
  call_no = 0;
  key[0] = '\0';
  for (icount = 0; icount < NCOUNTERS; icount++)
    counter[icount] = 0;
  on_off = 0;
  wanted = TRUE;

  /* If not skipping, check whether this might be the start of a region
     to be skipped because it is not relevant to this particular entry */
  if (skip == FALSE && !strncmp(input_line,"<BLOCK ",7))
    {
      /* Extract the field name identifying this block */
      strcpy(field_name,input_line+7);
      string_ptr = strstr(field_name,">");
      if (string_ptr != NULL)
        string_ptr[0] = '\0';

      /* Check for "!" prefix, signifying a negation */
      not = FALSE;
      if (input_line[7] == '!')
        {
          /* Set flag and remove the prefix from the name */
          not = TRUE;
          strcpy(tmp,field_name+1);
          strcpy(field_name,tmp);
        }

      /* Locate this field name */
      field_ptr = get_field(first_field_ptr,field_name,loop,instance);

      /* If field found, but its value is "FALSE" or "NONE", then treat
         as though the field does not exist */
      if (field_ptr != NULL && (!strcmp(field_ptr->value,"FALSE") ||
                                !strcmp(field_ptr->value,"NONE")))
        field_ptr = NULL;

      /* If not found, or value is FALSE, then want to skip this block */
      if ((field_ptr == NULL && not == FALSE) ||
          (field_ptr != NULL && not == TRUE))
        skip = TRUE;

      /* If skipping, then define end key */
      if (skip == TRUE)
        {
          /* Define the end key that will define the end of the lock
             to be skipped */
          strcpy(skip_end,"</BLOCK ");
          if (not == TRUE)
            strcat(skip_end,"!");
          strcat(skip_end,field_name);
          strcat(skip_end,">");
        }

      /* Don't want to store the BLOCK line itself */
      wanted = FALSE;
    }

  /* Check for SET command defining field to be set */
  else if (skip == FALSE && !strncmp(input_line,"<SET ",5))
    {
      /* Extract the field name */
      strcpy(field_name,input_line+5);
      string_ptr = strstr(field_name,"=");
      if (string_ptr != NULL)
        string_ptr[0] = '\0';

      /* Extract the field value */
      string_ptr = strstr(input_line,"=\"");
      if (string_ptr != NULL)
        {
          /* Get field value and remove final bracket */
          strcpy(value,string_ptr+2);
          string_ptr = strstr(value,"\">");
          if (string_ptr != NULL)
            string_ptr[0] = '\0';

          /* Annotate line */
          valid
            = annotate_line(value,first_field_ptr,loop,&call_no,
                            &on_off,counter,params);
        }
      else
        strcpy(value,"TRUE");

      /* From the string containing keyword-value pair */
      sprintf(string,"%s %s",field_name,value);

      /* Create a field record */
      field_ptr
        = create_field_record(first_field_ptr,string,FALSE,FALSE,
                              LOOP,-1);

      /* Don't want to store the SET line itself */
      wanted = FALSE;
    }

  /* If file containing more parameters to be read in, pick
     up and store the additional data */
  else if (skip == FALSE && !strncmp(input_line,"#data:",6))
    {
      /* Process line to replace all tags by their appropriate
         substitutions */
      valid
        = annotate_line(input_line,first_field_ptr,loop,&call_no,
                        &on_off,counter,params);

      /* Read in the parameter set */
      read_database_file(first_field_ptr,input_line+6,TRUE,TRUE,FALSE,
                         key,offset,start_record,nrecords,loop,-1,
                         replace_from,replace_by,params);

      /* Set flag to exclude current line */
      wanted = FALSE;
    }

  /* If skipping lines, check for the skip end */
  if (skip == TRUE)
    {
      /* If have the skip-end tag, unset flag */
      if (!strcmp(input_line,skip_end))
        {
          /* Unset flag */
          skip = FALSE;

          /* Don't want to store the /BLOCK line itself */
          wanted = FALSE;
        }
    }

  /* Check for line indicating this is the newest set of templates */
  if (!strcmp(input_line,"<!-- PDBsum template v.3 -->"))
    {
      /* Set parameter flag */
      params->version = 3;

      /* Get the bar colour field */
      field_ptr
        = get_field(first_field_ptr,"BAR_COLOUR",LOOP,params->instance);

      /* If valid and is default colour, replace by new colour */
      if (field_ptr != NULL &&
          !strcmp(field_ptr->value,BAR_COLOUR[STANDARD_PDB]))
        {
          /* Free up memory taken by previous value */
          free(field_ptr->value);

          /* Store new value */
          store_string(&(field_ptr->value),BAR_COLOUR[EBI_COLOUR]);
        }

      /* Don't need to write this line out */
      wanted = FALSE;
    }

  /* If not skipping, then ignore any block ends */
  else if (!strncmp(input_line,"</BLOCK ",8))
    wanted = FALSE;

  /* If line is wanted, store it */
  if (skip == FALSE && wanted == TRUE)
    {
      /* Store this line */
      store_string(&html_line[*nlines],input_line);

      /* Increment line count */
      (*nlines)++;
    }

  /* Return skip flag */
  return(skip);
}
/***********************************************************************

get_html_template  -  Read in and store the HTML template file

***********************************************************************/

int get_html_template(char **html_line,struct field *first_field_ptr,
                      int *loop,struct parameters *params)
{
  char file_name[FILENAME_LEN + 1], input_line[LINELEN + 1];
  char templates_dir[FILENAME_LEN + 1];
  char key[2], skip_end[100];

  char replace_from[2] = { "" };
  char replace_by[2] = { "" };

  int call_no, icount, len, line, nlines, on_off;
  int skip, valid;

  long counter[NCOUNTERS];

  static const int start_record = 0;
  static const int nrecords = -1;

  static const long offset = -1;

  struct field *field_ptr;

  FILE *file_ptr, *insert_file_ptr;

  /* Initialise variables */
  call_no = 0;
  for (icount = 0; icount < NCOUNTERS; icount++)
    counter[icount] = 0;
  key[0] = '\0';
  file_name[0] = '\0';
  line = 0;
  nlines = 0;
  on_off = 0;
  skip = FALSE;
  skip_end[0] = '\0';

  /* Get path to relevant templates directory */
  if (params->enzymes == TRUE)
    field_ptr
      = get_field(first_field_ptr,"ENZYME_REAL_TEMPLATES_DIR",LOOP,
                  params->instance);
  else if (params->profunc == TRUE)
    field_ptr
      = get_field(first_field_ptr,"PROFUNC_REAL_TEMPLATES_DIR",LOOP,
                  params->instance);
  else if (params->ligsearch == TRUE)
    field_ptr
      = get_field(first_field_ptr,"LIGSEARCH_REAL_TEMPLATES_DIR",LOOP,
                  params->instance);
  else if (params->varsite == TRUE)
    field_ptr
      = get_field(first_field_ptr,"VARSITE_REAL_TEMPLATES_DIR",LOOP,
                  params->instance);
  else if (params->ligplus == TRUE)
    field_ptr
      = get_field(first_field_ptr,"LIGPLUS_REAL_TEMPLATES_DIR",LOOP,
                  params->instance);
  else if (params->drat == TRUE)
    field_ptr
      = get_field(first_field_ptr,"DRAT_REAL_TEMPLATES_DIR",LOOP,
                  params->instance);
  else if (params->drugport == TRUE)
    field_ptr
      = get_field(first_field_ptr,"DRUGPORT_REAL_TEMPLATES_DIR",LOOP,
                  params->instance);
  else if (params->sas == TRUE)
    field_ptr
      = get_field(first_field_ptr,"SAS_REAL_TEMPLATES_DIR",LOOP,
                  params->instance);
  else
    field_ptr
      = get_field(first_field_ptr,"PDBSUM_REAL_TEMPLATES_DIR",LOOP,
                  params->instance);

  /* Form name of the template file */
  if (params->template_real_dir[0] != '\0')
    {
      strcpy(file_name,params->template_real_dir);
      strcat(file_name,"/");
    }
  else if (field_ptr != NULL)
    {
      strcpy(file_name,field_ptr->value);
      strcat(file_name,DIR_SEP);
    }
  strcpy(templates_dir,file_name);
  strcat(file_name,params->template_name);

  /* Open the template file */
  if ((file_ptr = fopen(file_name,"r")) !=  NULL)
    {
      /* Read through the template file, storing each line */
      while (fgets(input_line,LINELEN,file_ptr)!= NULL &&
             nlines < MAX_HTML_LINES - 1)
        {
          /* Truncate the string to strip off the newline */
          len = string_truncate(input_line,LINELEN);
          line++;

          /* Check for the "param" keyword giving name of parameter
             file associated with current HTML page, or the full name
             of the database.dat file to be read in */
          if (!strncmp(input_line,"#param:",7) ||
              !strncmp(input_line,"#database:",10))
            {
              /* Annotate line */
              valid
                = annotate_line(input_line,first_field_ptr,loop,&call_no,
                                &on_off,counter,params);

              /* Retrieve name of parameter file */
              if (!strncmp(input_line,"#param:",7))
                {
                  strcpy(file_name,templates_dir);
                  strcat(file_name,input_line+7);
                }

              /* Alternatively, use given file name */
              else
                strcpy(file_name,input_line+10);

              /* Read in all the parameters */
              read_database_file(first_field_ptr,file_name,TRUE,
                                 TRUE,FALSE,key,offset,start_record,
                                 nrecords,loop,-1,replace_from,
                                 replace_by,params);
            }

          /* Check for the "insert" keyword */
          else if (!strncmp(input_line,"#insert:",8))
            {
              /* Annotate line */
              valid
                = annotate_line(input_line,first_field_ptr,loop,&call_no,
                                &on_off,counter,params);

              /* Retrieve name of file to be inserted */
              strcpy(file_name,templates_dir);
              strcat(file_name,input_line+8);

              /* Open file and process all its lines */
              if ((insert_file_ptr = fopen(file_name,"r")) !=  NULL)
                {
                  /* Read through the template file, storing each line */
                  while (fgets(input_line,LINELEN,insert_file_ptr)!= NULL &&
                         nlines < MAX_HTML_LINES - 1)
                    {
                      /* Truncate the string to strip off the newline */
                      len = string_truncate(input_line,LINELEN);

                      /* Process line to check for starts and ends of
                         blocks */
                      skip
                        = process_blocks(input_line,skip,params->instance,
                                         loop,first_field_ptr,
                                         html_line,&nlines,skip_end,
                                         params);
                    }

                  /* Close the template file */
                  fclose(insert_file_ptr);
                }
            }

          /* Otherwise, process current line */
          else
            skip = process_blocks(input_line,skip,params->instance,loop,
                                  first_field_ptr,html_line,&nlines,
                                  skip_end,params);
        }

      /* Close the template file */
      fclose(file_ptr);
    }

  /* If no template file found, create HTML commands to report the
     error */
  else
    {
      store_string(&html_line[0],
                   "<H1><FONT COLOR=RED>*** ERROR</FONT></H1>");
      store_string(&html_line[1],
                   "<H2><FONT COLOR=BLUE>");
      store_string(&html_line[2],"HTML template missing!!</FONT>");
      store_string(&html_line[3],file_name);
      nlines = 4;
    }

  /* Return the number of lines read in */
  return(nlines);
}
/***********************************************************************

get_directories  -  Get the directory names required for the
                    process_gif_names and process_hrefs routines

***********************************************************************/

void get_directories(struct field *first_field_ptr,
                     char *templates_gif_dir,char *pdbsum_dir,
                     char *pdbsum_real_dir,char *templates_dir,
                     char *getpage_pl,char *getimage_pl,
                     char *getprofunc_image_pl,char *getligsearch_image_pl,
                     char *getvarsite_image_pl,char *cgi_dir,
                     char *consurf_dir,char *profunc_results_dir,
                     char *ligsearch_results_dir,char *varsite_results_dir,
                     char *ebi_include_dir,struct parameters *params)
{
  char db_name[NAMLEN + 1], field_name[NAMLEN + 1];
  char dir_name[FILENAME_LEN + 1];

  struct field *field_ptr;

  /* Get name of database */
  if (params->enzymes == TRUE)
    strcpy(db_name,"ENZYME");
  else if (params->profunc == TRUE)
    strcpy(db_name,"PROFUNC");
  else if (params->ligsearch == TRUE)
    strcpy(db_name,"LIGSEARCH");
  else if (params->varsite == TRUE)
    strcpy(db_name,"VARSITE");
  else if (params->ligplus == TRUE)
    strcpy(db_name,"LIGPLUS");
  else if (params->drat == TRUE)
    strcpy(db_name,"DRAT");
  else if (params->drugport == TRUE)
    strcpy(db_name,"DRUGPORT");
  else if (params->sas == TRUE)
    strcpy(db_name,"SAS");
  else
    strcpy(db_name,"PDBSUM");
 
  /* Get path to PDBsum directory */
  field_ptr = get_field(first_field_ptr,"PDBSUM_DIR",LOOP,
                        params->instance);

  /* Form directory name */
  pdbsum_dir[0] = '\0';
  if (field_ptr != NULL)
    {
      strcpy(pdbsum_dir,field_ptr->value);
      strcat(pdbsum_dir,DIR_SEP);
    }

  /* Get the PDBsum directory */
  field_ptr
    = get_field(first_field_ptr,"PDIR",LOOP,params->instance);

  /* Form directory name */
  pdbsum_real_dir[0] = '\0';
  if (field_ptr != NULL)
    {
      strcpy(pdbsum_real_dir,field_ptr->value);
      strcat(pdbsum_real_dir,DIR_SEP);
    }

  /* Get the directory holding the template GIF files */
  strcpy(field_name,db_name);
  strcat(field_name,"_TEMPLATES_GIF_DIR");
  field_ptr
    = get_field(first_field_ptr,field_name,LOOP,params->instance);

  /* Form full name */
  templates_gif_dir[0] = '\0';
  if (field_ptr != NULL)
    {
      if (params->pdb_type == PDBSUM1)
        {
          if (params->relinks == TRUE)
            strcpy(templates_gif_dir,"../gif/");
          else
            {
              strcpy(templates_gif_dir,"file://");
              strcat(templates_gif_dir,field_ptr->value);
            }
        }
      else
        strcpy(templates_gif_dir,field_ptr->value);
    }

  /* Form name of the template file */
  if (params->template_dir[0] != '\0')
    {
      if (!strncmp(params->template_dir,"C:",2) ||
          !strncmp(params->template_dir,"c:",2))
        {
          strcpy(templates_gif_dir,"file:///");
          strcat(templates_gif_dir,params->template_dir);
        }
      else
        strcpy(templates_gif_dir,params->template_dir);
      strcat(templates_gif_dir,"/gif");
    }

  /* Get the directory holding the template files */
  strcpy(field_name,db_name);
  strcat(field_name,"_TEMPLATES_DIR");
  field_ptr
    = get_field(first_field_ptr,field_name,LOOP,params->instance);

  /* Form full name */
  templates_dir[0] = '\0';
  if (field_ptr != NULL)
    strcpy(templates_dir,field_ptr->value);

  /* Form name of the template file */
  if (params->template_dir[0] != '\0')
    {
      if (!strncmp(params->template_dir,"C:",2) ||
          !strncmp(params->template_dir,"c:",2))
        {
          strcpy(templates_dir,"file:///");
          strcat(templates_dir,params->template_dir);
        }
      else
        strcpy(templates_dir,params->template_dir);
    }

  /* Get the directory holding the EBI template include files */
  strcpy(field_name,"EBI_INCLUDE_DIR");
  field_ptr
    = get_field(first_field_ptr,field_name,LOOP,params->instance);

  /* Form full name */
  ebi_include_dir[0] = '\0';
  if (field_ptr != NULL)
    strcpy(ebi_include_dir,field_ptr->value);

  /* Get call to Getpage perl script */
  if (params->drugport == TRUE)
    field_ptr = get_field(first_field_ptr,"GET_DRUGPORTPAGE_PL",LOOP,
                          params->instance);
  else if (params->ligsearch == TRUE)
    field_ptr = get_field(first_field_ptr,"GET_LIGSEARCH_PAGE_PL",LOOP,
                          params->instance);
  else if (params->varsite == TRUE)
    field_ptr = get_field(first_field_ptr,"GET_VARSITE_PAGE_PL",LOOP,
                          params->instance);
  else if (params->ligplus == TRUE)
    field_ptr = get_field(first_field_ptr,"GET_LIGPLUS_PAGE_PL",LOOP,
                          params->instance);
  else
    field_ptr = get_field(first_field_ptr,"GETPAGE_PL",LOOP,
                          params->instance);

  /* Form full name */
  getpage_pl[0] = '\0';
  if (field_ptr != NULL)
    strcpy(getpage_pl,field_ptr->value);

  /* Get path to getimg perl script */
  field_ptr = get_field(first_field_ptr,"GET_IMAGE",LOOP,params->instance);

  /* Form full name */
  getimage_pl[0] = '\0';
  if (field_ptr != NULL)
    strcpy(getimage_pl,field_ptr->value);

  /* Get path to perl script that retrieves ProFunc images */
  field_ptr = get_field(first_field_ptr,"GET_PROFUNC_IMG",LOOP,
                        params->instance);

  /* Form full name */
  getprofunc_image_pl[0] = '\0';
  if (field_ptr != NULL)
    strcpy(getprofunc_image_pl,field_ptr->value);

  /* Get path to perl script that retrieves LigSearch images */
  field_ptr = get_field(first_field_ptr,"GET_LIGSEARCH_IMAGE",LOOP,
                        params->instance);

  /* Form full name */
  getligsearch_image_pl[0] = '\0';
  if (field_ptr != NULL)
    strcpy(getligsearch_image_pl,field_ptr->value);

  /* Get path to perl script that retrieves VarSite images */
  field_ptr = get_field(first_field_ptr,"GET_VARSITE_IMAGE",LOOP,
                        params->instance);

  /* Form full name */
  getvarsite_image_pl[0] = '\0';
  if (field_ptr != NULL)
    strcpy(getvarsite_image_pl,field_ptr->value);

  /* Get name of ConSurf data directory */
  field_ptr = get_field(first_field_ptr,"NEW_CONSURF_DIR",LOOP,
                        params->instance);

  /* Form full name */
  consurf_dir[0] = '\0';
  if (field_ptr != NULL)
    strcpy(consurf_dir,field_ptr->value);

  /* Get path to ProFunc results directory */
  field_ptr = get_field(first_field_ptr,"PROFUNC_RESULTS_DIR",LOOP,
                        params->instance);

  /* Form full name */
  profunc_results_dir[0] = '\0';
  if (field_ptr != NULL)
    strcpy(profunc_results_dir,field_ptr->value);

  /* Get path to LigSearch results directory */
  field_ptr = get_field(first_field_ptr,"LIGSEARCH_REAL_RESULTS",LOOP,
                        params->instance);

  /* Form full name */
  ligsearch_results_dir[0] = '\0';
  if (field_ptr != NULL)
    strcpy(ligsearch_results_dir,field_ptr->value);

  /* Get path to VarSite results directory */
  field_ptr = get_field(first_field_ptr,"VARSITE_REAL_UNIPROT_DIR",LOOP,
                        params->instance);

  /* Form full name */
  varsite_results_dir[0] = '\0';
  if (field_ptr != NULL)
    strcpy(varsite_results_dir,field_ptr->value);

  /* Get path for PDBsum cgi-bin directory */
  strcpy(field_name,db_name);
  strcat(field_name,"_CGI_DIR");
  field_ptr = get_field(first_field_ptr,field_name,LOOP,
                        params->instance);

  /* Form full name */
  cgi_dir[0] = '\0';
  if (field_ptr != NULL)
    strcpy(cgi_dir,field_ptr->value);
}
/***********************************************************************

check_for_string  -  Check for comparison string within curly brackets
                     in a SKIP command

***********************************************************************/

void check_for_string(char *field_name,char *ftype,char *string)
{
  int i, ipos, len, value;

  char replace_string[LINELEN + 1];

  char *string1_ptr, *string2_ptr;

  /* Initialise value */
  ftype[0] = '\0';
  string[0] = '\0';

  /* Check for presence of open and close brackets */
  string1_ptr = strstr(field_name,"{");
  string2_ptr = strstr(field_name,"}");

  /* If both present, then process */
  if (string1_ptr != NULL && string2_ptr != NULL &&
      string2_ptr > string1_ptr)
    {
      /* Have both brackets, so extract the string between them */
      ipos = string1_ptr - field_name + 1;
      len = string2_ptr - string1_ptr - 1;
      strncpy(string,field_name + ipos,len);
      string[len] = '\0';

      /* Splice out the bracketed string from field name */
      strncpy(replace_string,string1_ptr,len + 2);
      replace_string[len + 2] = '\0';

      /* Perform splice */
      perform_splice(field_name,replace_string,"");

      /* Mark as a string equality comparison */
      strcpy(ftype,"str");
    }
}
/***********************************************************************

check_field  -  Test for existence of a given field

***********************************************************************/

int check_field(struct field *first_field_ptr,char *keyword,char *name,
                int *loop,char *skip_end,int *ival,char *ftype,
                struct parameters *params)
{
  char field_name[NAMLEN + 1], tmp_string[LINELEN + 1];
  char string[LINELEN + 1], tmp[LINELEN + 1];

  char *string_ptr;

  int ivalue, value;
  int not, skip;

  struct field dummy_field;

  struct field *field_ptr;

  /* Initialise variables */
  ftype[0] = '\0';
  ivalue = 0;
  skip = FALSE;

  /* Get start of field name */
  string_ptr = strstr(name,keyword);

  /* If have keyword, look for space prior to field name */
  if (string_ptr != NULL)
    {
      /* Save for forming skip-end string, if required later */
      strcpy(tmp_string,string_ptr);

      /* Find location of space */
      string_ptr = strstr(string_ptr," ");

      /* If found, then have start of field name */
      if (string_ptr != NULL)
        {
          /* Extract field name */
          strcpy(field_name,string_ptr + 1);

          /* Lop off close bracket */
          string_ptr = strstr(field_name,">");

          /* If valid, terminate string here */
          if (string_ptr != NULL)
            string_ptr[0] = '\0';

          /* Check for "!" prefix, signifying a negation */
          not = FALSE;
          if (field_name[0] == '!')
            {
              /* Set flag and remove the prefix from the name */
              not = TRUE;
              strcpy(tmp,field_name+1);
              strcpy(field_name,tmp);
            }

          /* Check for an operator and number in brackets after the
             field name giving logical operation to be performed
             (eg eq1, gt1, etc) */
          ivalue = check_for_brackets(field_name,"[","]",ftype);

          /* Check for a string comparison in curly brackets */
          if (ftype[0] == '\0')
            check_for_string(field_name,ftype,string);

          /* Check for existence of this field */
          field_ptr = get_field(first_field_ptr,field_name,loop,
                                params->instance);

          /* If field found, but its value is "FALSE" or "NONE", then treat
             as though the field does not exist */
          if (field_ptr != NULL && (!strcmp(field_ptr->value,"FALSE") ||
                                    !strcmp(field_ptr->value,"NONE")))
            field_ptr = NULL;

          /* If checking for a blank field, set to NULL if empty */
          if (field_ptr != NULL && !strcmp(ftype,"empty"))
            {
              /* Check for possible ways field can be blank */
              if (!strcmp(field_ptr->value,"-") ||
                  strlen(field_ptr->value) == 0)
                field_ptr = NULL;
            }

          /* If field not found, then skip the block if field required */
          if (field_ptr == NULL)
            {
              /* Set skip flag accordingly */
              if (not == TRUE)
                skip = TRUE;
              else
                skip = FALSE;

              /* If field name is "TRUE", or checking for empty,
                 reverse flag */
              if (!strcmp(field_name,"TRUE") || !strcmp(ftype,"empty"))
                skip = 1 - skip;
            }

          /* Otherwise, see if skip depends on field value */
          else if (ftype[0] != '\0')
            {
              /* Get the field value */
              value = atoi(field_ptr->value);

              /* Check for string comparison */
              if (!strcmp(ftype,"str"))
                {
                  if (!strcmp(string,field_ptr->value))
                    skip = TRUE;
                  else
                    skip = FALSE;

                  /* If field negated, then reverse result */
                  if (not == TRUE)
                    skip = 1 - skip;
                }

              /* Check for non blank field */
              else if (!strcmp(ftype,"empty"))
                {
                  /* Set flag to skip */
                  skip = FALSE;

                  /* If field negated, then reverse result */
                  if (not == TRUE)
                    skip = 1 - skip;
                }

              /* Compare numerical field and required value depending on
                 operator */
              else if (!strcmp(ftype,"eq") && value == ivalue)
                skip = TRUE;
              else if (!strcmp(ftype,"ne") && value != ivalue)
                skip = TRUE;
              else if (!strcmp(ftype,"gt") && value > ivalue)
                skip = TRUE;
              else if (!strcmp(ftype,"ge") && value >= ivalue)
                skip = TRUE;
              else if (!strcmp(ftype,"lt") && value < ivalue)
                skip = TRUE;
              else if (!strcmp(ftype,"le") && value <= ivalue)
                skip = TRUE;
              else
                skip = FALSE;
            }

          /* Otherwise, skip unless "not" selected */
          else
            {
              if (not == FALSE)
                skip = TRUE;
              else
                skip = FALSE;
            }

          /* If skipping, then define end key */
          if (skip == TRUE)
            {
              /* Store name of end-skip field */
              strcpy(skip_end,"</");
              strcat(skip_end,tmp_string+1);
              string_ptr = strstr(skip_end,">");
              if (string_ptr != NULL)
                string_ptr[1] = '\0';
            }
        }
    }

  /* Return value associated with any field operator */
  *ival = ivalue;

  /* Return whether field found or not */
  return(skip);
}
/***********************************************************************

check_for_file  -  Test for existence of a given file

***********************************************************************/

int check_for_file(char *name,char *skip_end)
{
  char file_name[NAMLEN + 1];

  char *string_ptr;

  int done, have_file, not, skip;

  FILE *file_ptr;

  /* Initialise variables */
  skip = FALSE;
  if (name[6] == '!')
    not = TRUE;
  else
    not = FALSE;

  /* Extract the file name */
  strcpy(file_name,name);
  string_ptr = strstr(name,"\"");
  if (string_ptr != NULL)
    strcpy(file_name,string_ptr+1);

  /* Close the file name */
  string_ptr = strstr(file_name,"\"");
  if (string_ptr != NULL)
    string_ptr[0] = '\0';

  /* If not on a unix system, replace any slashes with directory
     separators */
  if (strcmp(DIR_SEP,"/"))
    {
      /* Initialise flag */
      done = FALSE;

      /* Loop until all directory separators have been replaced */
      while (done == FALSE)
        {
          /* Locate next slash */
          string_ptr = strstr(file_name,"/");
          
          /* If found, then replace */
          if (string_ptr != NULL)
            string_ptr[0] = DIR_SEP[0];

          /* Otherwise, we're done */
          else
            done = TRUE;
        }
    }

  /* Assume file not there */
  have_file = FALSE;

  /* Open file to check whether it exists */
  if ((file_ptr = fopen(file_name,"r")) !=  NULL)
    {
      /* Set flag that have file */
      have_file = TRUE;

      /* Close the file */
      fclose(file_ptr);
    }

  /* Determine whether skip section is to be skipped */
  if (not == TRUE)
    if (have_file == TRUE)
      skip = FALSE;
    else
      skip = TRUE;
  else
    if (have_file == TRUE)
      skip = TRUE;
    else
      skip = FALSE;

  /* If skipping, then define end key */
  if (skip == TRUE)
    {
      /* Store name of end-skip field */
      strcpy(skip_end,"</");
      strcat(skip_end,name+1);
      string_ptr = strstr(skip_end,">");
      if (string_ptr != NULL)
        string_ptr[1] = '\0';
    }

  /* Return whether file found or not */
  return(skip);
}
/***********************************************************************

check_in_file  -  Test for presence of a given string in a given file

***********************************************************************/

int check_in_file(char *name,char *skip_end,int *line_pos)
{
  char file_name[NAMLEN + 1], string[NAMLEN + 1], tmp[NAMLEN + 1];
  char input_line[LINELEN + 1];

  char *string_ptr;

  int at_line, line;
  int done, have_string, not, skip;

  FILE *file_ptr;

  /* Initialise variables */
  *line_pos = -1;
  skip = FALSE;
  if (name[6] == '!')
    not = TRUE;
  else
    not = FALSE;

  /* Extract the file name */
  strcpy(file_name,name);
  string_ptr = strstr(name,"\"");
  if (string_ptr != NULL)
    strcpy(file_name,string_ptr+1);

  /* Close the file name */
  string_ptr = strstr(file_name,"\"");
  if (string_ptr != NULL)
    string_ptr[0] = '\0';

  /* Extract the string to be searched */
  strcpy(string,name);
  string_ptr = strstr(name,"[");
  if (string_ptr != NULL)
    {
      strcpy(tmp,string_ptr+1);
      strcpy(string,tmp);
    }

  /* Close the string */
  string_ptr = strstr(string,"]");
  if (string_ptr != NULL)
    string_ptr[0] = '\0';

  /* If not on a unix system, replace any slashes with directory
     separators */
  if (strcmp(DIR_SEP,"/"))
    {
      /* Initialise flag */
      done = FALSE;

      /* Loop until all directory separators have been replaced */
      while (done == FALSE)
        {
          /* Locate next slash */
          string_ptr = strstr(file_name,"/");
          
          /* If found, then replace */
          if (string_ptr != NULL)
            string_ptr[0] = DIR_SEP[0];

          /* Otherwise, we're done */
          else
            done = TRUE;
        }
    }

  /* Assume file and string not there */
  have_string = FALSE;
  line = 0;
  at_line = -1;

  /* Open file to check whether it exists */
  if ((file_ptr = fopen(file_name,"r")) !=  NULL)
    {
      /* Loop while reading through the file */
      while (fgets(input_line,LINELEN,file_ptr) != NULL &&
             at_line == -1)
        {
          /* Increment line count */
          line++;

          /* Test for presence of string in this line */
          string_ptr = strstr(input_line,string);
          if (string_ptr != NULL)
            at_line = line;
        }

      /* Close the file */
      fclose(file_ptr);
    }

  /* If item found, then set flag */
  if (at_line > -1)
    have_string = TRUE;

  /* Determine whether skip section is to be skipped */
  if (not == TRUE)
    if (have_string == TRUE)
      skip = FALSE;
    else
      skip = TRUE;
  else
    if (have_string == TRUE)
      skip = TRUE;
    else
      skip = FALSE;

  /* If skipping, then define end key */
  if (skip == TRUE)
    {
      /* Store name of end-skip field */
      strcpy(skip_end,"</");
      strcat(skip_end,name+1);
      string_ptr = strstr(skip_end,">");
      if (string_ptr != NULL)
        string_ptr[1] = '\0';
    }

  /* Return line position at which string was located */
  *line_pos = at_line;

  /* Return whether string found or not */
  return(skip);
}
/***********************************************************************

end_loop  -  Determine whether we are currently in the middle of, or at
             the end of the current loop

***********************************************************************/

void end_loop(struct field *first_field_ptr,int level,int loop,
              int nloops,struct parameters *params)
{
  char on_name[NAMLEN + 1], off_name[NAMLEN + 1], string[NAMLEN + 1];

  struct field *field_ptr;

  /* Determine whether in mid-loop or at end */
  if (loop < nloops)
    {
      /* In mid-loop */
      strcpy(on_name,"MID");
      strcpy(off_name,"END");
    }
  else
    {
      /* At end of the loop */
      strcpy(off_name,"MID");
      strcpy(on_name,"END");
    }

  /* Switch on-field to TRUE */
  sprintf(string,"%sLOOP%d TRUE",on_name,level + 1);

  /* Create a field record */
  field_ptr
    = create_field_record(first_field_ptr,string,FALSE,FALSE,LOOP,-1);

  /* Switch off-field to TRUE */
  sprintf(string,"%sLOOP%d FALSE",off_name,level + 1);

  /* Create a field record */
  field_ptr
    = create_field_record(first_field_ptr,string,FALSE,FALSE,LOOP,-1);
}
/***********************************************************************

fix_plurals  -  Replace all <s% tags by "s" or blank depending on whether
                the field in question is one or more than one

***********************************************************************/

void fix_plurals(char *tmp_line,struct field *first_field_ptr,
                 int *loop,struct parameters *params)
{
  char field_name[NAMLEN + 1], replace_string[LINELEN + 1];
  char value[LINELEN + 1], number_string[20], s[10];

  char *string_ptr;

  int ivalue, nloop;
  int and, done, has_have, is_are, splice;

  struct field *field_ptr;

  /* Initialise variables */
  done = FALSE;
  nloop = 0;
  splice = FALSE;

  /* Loop until no more tags to be replaced */
  while (done == FALSE)
    {
      /* Assume we're done */
      done = TRUE;

      /* Check for any plurals still in the line */
      if ((strstr(tmp_line,"<s% ") != NULL ||
           strstr(tmp_line,"<S% ") != NULL ||
           strstr(tmp_line,"<&% ") != NULL ||
           strstr(tmp_line,"<h% ") != NULL ||
           strstr(tmp_line,"<i% ") != NULL) && nloop < 100)
        {
          /* Increment loop count to ensure we don't hang in the case
             of an error in the template */
          nloop++;
          and = FALSE;
          is_are = FALSE;
          has_have = FALSE;

          /* Get the location of the s to be replaced */
          string_ptr = strstr(tmp_line,"<s% ");
          if (string_ptr == NULL)
            string_ptr = strstr(tmp_line,"<S% ");
          if (string_ptr == NULL)
            {
              string_ptr = strstr(tmp_line,"<&% ");
              if (string_ptr != NULL)
                and = TRUE;
            }
          if (string_ptr == NULL)
            {
              string_ptr = strstr(tmp_line,"<h% ");
              if (string_ptr != NULL)
                has_have = TRUE;
            }
          if (string_ptr == NULL)
            {
              string_ptr = strstr(tmp_line,"<i% ");
              if (string_ptr != NULL)
                is_are = TRUE;
            }

          /* If found, then proceed */
          if (string_ptr != NULL)
            {
              /* Determine whether lower- or upper-case s to be used */
              s[0] = string_ptr[1];
              s[1] = '\0';

              /* Form the string to be replaced */
              strcpy(replace_string,string_ptr);
              string_ptr = strstr(replace_string,">");
              if (string_ptr != NULL)
                string_ptr[1] = '\0';

              /* Extract the field name to be checked */
              strcpy(value,replace_string + 4);
              string_ptr = strstr(value,">");
              if (string_ptr != NULL)
                string_ptr[0] = '\0';
              strcpy(field_name,value);

              /* Get field */
              field_ptr = get_field(first_field_ptr,field_name,loop,
                                    params->instance);

              /* If valid, perform the splice */
              if (field_ptr != NULL)
                {
                  /* Get field value */
                  strcpy(number_string,field_ptr->value);
                  ivalue = atoi(number_string);
                }
              else
                ivalue = 0;

              /* If greater than 1, want to add s to make plural */
              if ((and == TRUE && ivalue == 0) ||
                  (and == FALSE && ivalue < 2 && ivalue != 0))
                s[0] = '\0';

              /* Check for special case of is/are */
              if (is_are == TRUE)
                {
                  /* If a single value, then use "is" */
                  if (ivalue == 1)
                    strcpy(s,"is");
                  else
                    strcpy(s,"are");
                }

              /* Check for special case of has/have */
              else if (has_have == TRUE)
                {
                  /* If a single value, then use "has" */
                  if (ivalue == 1)
                    strcpy(s,"has");
                  else
                    strcpy(s,"have");
                }

             /* Perform splice */
              perform_splice(tmp_line,replace_string,s);

              /* Unset flag as there may be more tags left to find */
              done = FALSE;
            }
        }
    }
}
/***********************************************************************

insert_file  -  Pick up the required file and insert into the output
                file here

***********************************************************************/

void insert_file(FILE *outfile_ptr,char *file_name,char *pdbsum_real_dir,
                 char *pdbsum_dir,char *getpage_pl,char *getimage_pl,
                 char *templates_gif_dir,char *profunc_results_dir,
                 char *ligsearch_results_dir,char *varsite_results_dir,
                 struct field *first_field_ptr,int *loop,int *cno,
                 int *onf,long *counter,struct parameters *params)
{
  char input_line[LINELEN + 1], number_string[20];
  char full_link[LINELEN + 1], replace_string[LINELEN + 1];
  char html_file[FILENAME_LEN + 1];

  char *string_ptr;

  int call_no, on_off;
  int done, valid, wanted;

  FILE *file_ptr;

  /* Initialise variables */
  on_off = *onf;
  call_no = *cno;
  done = FALSE;

  /* If file is called "filename", use name given in command-line
     parameters */
  if (!strcmp(file_name,"filename"))
    {
      /* If have the full file name, use it as it stands */
      if (params->fullname == TRUE)
        strcpy(file_name,params->include_file);

      /* Otherwise, assume file is in this entry's PDBsum directory,
         so form the full path */
      else
        {
          strcpy(file_name,pdbsum_real_dir);
          strcat(file_name,DIR_SEP);
          strcat(file_name,params->include_file);
          strcat(file_name,".html");
        }
    }

  /* For ProFunc, get the file from the ProFunc results directory */
  else if (params->profunc == TRUE)
    {
      strcpy(full_link,profunc_results_dir);
      strcat(full_link,DIR_SEP);
      strcat(full_link,file_name);
      strcpy(file_name,full_link);
    }

  /* For LigSearch, get the file from the LigSearch results directory */
  else if (params->ligsearch == TRUE)
    {
      strcpy(full_link,ligsearch_results_dir);
      strcat(full_link,DIR_SEP);
      strcat(full_link,file_name);
      strcpy(file_name,full_link);
    }

  /* For VarSite, get the file from the VarSite results directory */
  else if (params->varsite == TRUE)
    {
      strcpy(full_link,varsite_results_dir);
      strcat(full_link,DIR_SEP);
      strcat(full_link,file_name);
      strcpy(file_name,full_link);
    }

  /* Open the summary file */
  if ((file_ptr = fopen(file_name,"r")) !=  NULL)
    {
      /* Read through the file, writing each line out */
      while (fgets(input_line,LINELEN,file_ptr) != NULL &&
             done == FALSE)
        {
          /* Assume line is wanted */
          wanted = TRUE;

          /* Check for one of the lines to be skipped */
          if (!strncmp(input_line,"<html",5) ||
              !strncmp(input_line,"<head",5) ||
              !strncmp(input_line,"<title",6) ||
              !strncmp(input_line,"</head",6) ||
              !strncmp(input_line,"<body",5))
            wanted = FALSE;

          /* Check for end */
          if (params->profunc == FALSE && params->sas == FALSE &&
              params->ligsearch == FALSE && params->varsite == FALSE &&
              strstr(input_line,"</PRE>"))
              done = TRUE;

          /* Check for link to PROMOTIF pictures */
          if ((string_ptr = strstr(input_line,"<a href =\"pdb")) != NULL)
            {
              /* Get current level */
              sprintf(number_string,"%d",params->instance[0]);

              /* Get the name of the html_file to be inserted */
              strcpy(html_file,string_ptr+10);
              string_ptr = strstr(html_file,".html\"");
              if (string_ptr != NULL)
                string_ptr[0] = '\0';

              /* Form full link to call page with new inserted file */
              strcpy(full_link,"href=\"");
              strcat(full_link,getpage_pl);
              string_ptr = strstr(full_link,"<FLD");
              if (string_ptr != NULL)
                string_ptr[0] = '\0';
              strcat(full_link,params->pdb_code);

              /* Add security details for non-standard files */
              if (params->pdb_type != STANDARD_PDB &&
                  params->pdb_type != ALPHA_FOLD)
                {
                  strcat(full_link,"&pdb_type=");
                  if (params->pdb_type == UPLOAD_PDB)
                    strcat(full_link,"UPLOAD");
                  else if (params->pdb_type == PROFUNC)
                    strcat(full_link,"PROFUNC");
                  else if (params->pdb_type == LIGSEARCH)
                    strcat(full_link,"LIGSEARCH");
                  else if (params->pdb_type == VARSITE)
                    strcat(full_link,"VARSITE");
                  else if (params->pdb_type == LIGPLUS)
                    strcat(full_link,"LIGPLUS");
                  else
                    strcat(full_link,"MCSG");
                }

              /* Add security code if there is one */
              if (params->security_code[0] != '\0')
                {
                  strcat(full_link,"&code=");
                  strcat(full_link,params->security_code);
                }

              /* Add PROMOTIF details */
              strcat(full_link,"&template=protein.html&o=PROMOTIF&l=");
              strcat(full_link,number_string);
              strcat(full_link,"&i=");
              strcat(full_link,html_file);
              strcat(full_link,"\"");

              /* Get the part of the string that is to be replaced */
              replace_string[0] = '\0';
              string_ptr = strstr(input_line,"href =\"pdb");
              if (string_ptr != NULL)
                strcpy(replace_string,string_ptr);
              string_ptr = strstr(replace_string,">");
              if (string_ptr != NULL)
                string_ptr[0] = '\0';

              /* Splice out current path and replace with correct path */
              perform_splice(input_line,replace_string,full_link);
            }

          /* Check for PROMOTIF GIF image */
          if ((string_ptr = strstr(input_line,"<img src=pdb")) != NULL)
            {
              /* Form the full path to the GIF file */
              strcpy(full_link,"<img src=");

              /* If image in standard PDBsum directory, just refrerence
                 it directly */
              if (params->pdb_type == STANDARD_PDB || 
                  params->pdb_type == ALPHA_FOLD)
                {
                  strcat(full_link,pdbsum_dir);
                  strcat(full_link,DIR_SEP);
                }

              /* Otherwise need to retrieve the image on the fly using the
                 getimg.pl script */
              else
                {
                  strcat(full_link,getimage_pl);
                  strcat(full_link,"&source=pdbsum&pdb_type=");
                  if (params->pdb_type == UPLOAD_PDB)
                    strcat(full_link,"UPLOAD");
                  else
                    strcat(full_link,"MCSG");
                  strcat(full_link,"&pdb_code=");
                  strcat(full_link,params->pdb_code);
                  strcat(full_link,"&code=");
                  strcat(full_link,params->security_code);
                  strcat(full_link,"&file=");
                }
              strcat(full_link,"pdb");

              /* Splice out current path and replace with correct path */
              perform_splice(input_line,"<img src=pdb",full_link);
            }

          /* Check for sites.html GIF image */
          if ((string_ptr = strstr(input_line,"<IMG SRC=\"../")) != NULL)
            {
              /* Form the full path to the GIF file */
              strcpy(full_link,"<img src=\"");
              strcat(full_link,templates_gif_dir);
              strcat(full_link,DIR_SEP);

              /* Splice out current path and replace with correct path */
              perform_splice(input_line,"<IMG SRC=\"../",full_link);
            }

          /* Check for other GIF images in PDBsum directory */
          if ((string_ptr = strstr(input_line,"<IMG SRC=\"./")) != NULL)
            {
              /* Form the full path to the GIF file */
              strcpy(full_link,"<img src=\"");

              /* Need to retrieve the image on the fly using the
                 getimg.pl script */
              strcat(full_link,getimage_pl);
              strcat(full_link,"&source=pdbsum");
              if (params->pdb_type == UPLOAD_PDB)
                strcat(full_link,"&pdb_type=UPLOAD");
              else if (params->pdb_type == MCSG_PDB)
                strcat(full_link,"&pdb_type=MCSG");
              strcat(full_link,"&pdb_code=");
              strcat(full_link,params->pdb_code);

              /* Add security code if there is one */
              if (params->security_code[0] != '\0')
                {
                  strcat(full_link,"&code=");
                  strcat(full_link,params->security_code);
                }
              strcat(full_link,"&file=");

              /* Splice out current path and replace with correct path */
              perform_splice(input_line,"<IMG SRC=\"./",full_link);
            }

          /* Process line to replace any tags by their appropriate
             substitutions */
          if (strstr(input_line,TAG_END) != NULL)
            valid
              = annotate_line(input_line,first_field_ptr,loop,&call_no,
                              &on_off,counter,params);

          /* Write to output file */
          if (wanted == TRUE)
            fprintf(outfile_ptr,"%s",input_line);
        }

      /* Close the file */
      fclose(file_ptr);
    }

  /* If file not there, write out a blank */
  else
    {
      fprintf(outfile_ptr,"&nbsp;\n");
      //printf("*** Warning. Unable to find include file: %s\n",
      //     file_name);
    }

  /* Return updated variables */
  *cno = call_no;
  *onf = on_off;
}
/***********************************************************************

create_iname_record  -  Create and initialise a new image name record

***********************************************************************/

struct iname *create_iname_record(struct iname **fst_iname_ptr,
                                  struct iname **lst_iname_ptr,
                                  char *image_name)
{
  int i;

  struct iname *iname_ptr, *first_iname_ptr, *last_iname_ptr;

  /* Initialise image pointers */
  iname_ptr = NULL;
  first_iname_ptr = *fst_iname_ptr;
  last_iname_ptr = *lst_iname_ptr;

  /* Allocate memory for structure to hold image info */
  iname_ptr = (struct iname *) malloc(sizeof(struct iname));
  if (iname_ptr == NULL)
    {
      printf("*** Can't allocate memory for struct iname\n");
      exit (1);
    }

  /* Store the image details */
  store_string(&(iname_ptr->image_name),image_name);

  /* If this is the very first image, then store pointer */
  if (first_iname_ptr == NULL)
    first_iname_ptr = iname_ptr;

  /* Add link from previous image to the current one */
  if (last_iname_ptr != NULL)
    last_iname_ptr->next_iname_ptr = iname_ptr;
  last_iname_ptr = iname_ptr;
  iname_ptr->next_iname_ptr = NULL;

  /* Return current pointers */
  *fst_iname_ptr = first_iname_ptr;
  *lst_iname_ptr = last_iname_ptr;

  /* Return new image pointer */
  return(iname_ptr);
}
/***********************************************************************

find_image_name  -  See if we already have this image name

***********************************************************************/

struct iname *find_image_name(struct iname *first_iname_ptr,
                              char *image_name)
{
  struct iname *iname_ptr, *found_iname_ptr;

  /* Initialise variables */
  found_iname_ptr = NULL;

  /* Get pointer to the first image_name */
  iname_ptr = first_iname_ptr;

  /* Loop over all the image_names to find the one of interest */
  while (iname_ptr != NULL && found_iname_ptr == NULL)
    {
      /* If have a match, store the pointer */
      if (!strcmp(iname_ptr->image_name,image_name))
        found_iname_ptr = iname_ptr;

      /* Get pointer to the next image_name */
      iname_ptr = iname_ptr->next_iname_ptr;
    }

  /* Return the matching image_name pointer */
  return(found_iname_ptr);
}
/***********************************************************************

extract_image_name  -  Extract the name of the image

***********************************************************************/

int extract_image_name(char *tmp_line,struct field *first_field_ptr,
                       int *loop,int *cno,int *onf,long *counter,
                       struct iname **fst_iname_ptr,
                       struct iname **lst_iname_ptr,int nimage_names,
                       struct parameters *params)
{
  char image_name[LINELEN + 1];

  char *string_ptr;

  int call_no, dot_pos;
  int on_off;

  struct iname *iname_ptr, *first_iname_ptr, *last_iname_ptr;

  /* Initialise variables */
  call_no = *cno;
  on_off = *onf;

  /* Initialise image pointers */
  iname_ptr = NULL;
  first_iname_ptr = *fst_iname_ptr;
  last_iname_ptr = *lst_iname_ptr;

  /* Extract the start of the image name */
  image_name[0] = '\0';
  string_ptr = strstr(tmp_line,"gif/");
  if (string_ptr != NULL)
    strcpy(image_name,string_ptr + 4);

  /* Find the dot before the extension */
  string_ptr = strchr(image_name,'.');
  dot_pos = 0;
  if (string_ptr != NULL)
    {
      /* Find end of image name */
      string_ptr[dot_pos + 4] = '\0';

      /* If name contains any annotations, then convert */
      if (strstr(image_name,"FLD") != NULL)
        annotate_line(image_name,first_field_ptr,loop,&call_no,
                      &on_off,counter,params);

      /* See if we already have this image name */
      iname_ptr = find_image_name(first_iname_ptr,image_name);

      /* If this is a new name, then save it */
      if (iname_ptr == NULL)
        {
          /* Create a new iname record */
          iname_ptr
            = create_iname_record(&first_iname_ptr,&last_iname_ptr,
                                  image_name);

          /* Increment count of unique image names */
          nimage_names++;
        }
    }

  /* Return current pointers */
  *fst_iname_ptr = first_iname_ptr;
  *lst_iname_ptr = last_iname_ptr;

  /* Return the number of unique image names */
  return(nimage_names);
}
/***********************************************************************

process_gif_names  -  Replace all tags by appropriate substitutions

***********************************************************************/

int process_gif_names(char *tmp_line,char *templates_gif_dir,
                      char *pdbsum_dir,char *getimage_pl,
                      char *getprofunc_image_pl,
                      char *getligsearch_image_pl,
                      char *getvarsite_image_pl,
                      struct field *first_field_ptr,int *loop,
                      int *cno,int *onf,long *counter,
                      struct iname **fst_iname_ptr,
                      struct iname **lst_iname_ptr,int nimage_names,
                      struct parameters *params)
{
  char full_link[LINELEN + 1];

  int call_no;
  int on_off, replaced;

  struct iname *iname_ptr, *first_iname_ptr, *last_iname_ptr;

  /* Set flag indicating whether a GIF URL has been replaced */
  call_no = *cno;
  on_off = *onf;
  replaced = TRUE;

  /* Initialise image pointers */
  iname_ptr = NULL;
  first_iname_ptr = *fst_iname_ptr;
  last_iname_ptr = *lst_iname_ptr;

  /* Loop until no more GIFs to be redirected */
  while (replaced == TRUE)
    {
      /* Unset flag */
      replaced = FALSE;
      full_link[0] = '\0';

      /* If line contains a standard PDBsum GIF image, then need to
         replace with the full path */
      if (strstr(tmp_line,"src=\"gif/") != NULL)
        {
          /* Form the full path to the GIF file */
          strcat(full_link,"src=\"");
          strcat(full_link,templates_gif_dir);
          strcat(full_link,DIR_SEP);

          /* For relative links, point to local directory */
          if (params->relinks == TRUE)
            {
              /* Set the link */
              strcpy(full_link,"src=\"../gif/");

	      /* If this is a common image, keep the path but change
	         src to SRC */
	      if (strstr(tmp_line," c=\"T\" ") != NULL)
		strcpy(full_link,"SRC=\"gif/");

              /* Extract the image name */
              nimage_names
                = extract_image_name(tmp_line,first_field_ptr,loop,
                                     &call_no,&on_off,counter,
                                     &first_iname_ptr,&last_iname_ptr,
                                     nimage_names,params);
            }

          /* Splice out current path and replace with correct path */
          perform_splice(tmp_line,"src=\"gif/",full_link);

          /* Set flag */
          replaced = TRUE;
        }

      /* If line contains a page "background image", add full path */
      else if (strstr(tmp_line,"background=\"gif/"))
        {
          /* Form the full path to the GIF file */
          strcat(full_link,"background=\"");
          strcat(full_link,templates_gif_dir);
          strcat(full_link,DIR_SEP);

          /* For relative links, point to local directory */
          if (params->relinks == TRUE)
            {
              /* Set the link */
              strcpy(full_link,"background=\"../gif/");

              /* Extract the image name */
              nimage_names
                = extract_image_name(tmp_line,first_field_ptr,loop,
                                     &call_no,&on_off,counter,
                                     &first_iname_ptr,&last_iname_ptr,
                                     nimage_names,params);
            }

          /* Splice out current path and replace with correct path */
          perform_splice(tmp_line,"background=\"gif/",full_link);

          /* Set flag */
          replaced = TRUE;
        }

      /* If line contains a "background image", add full path */
      else if (strstr(tmp_line,"url(gif/") != NULL)
        {
          /* Form the full path to the GIF file */
          strcat(full_link,"url(");
          strcat(full_link,templates_gif_dir);
          strcat(full_link,DIR_SEP);

          /* Splice out current path and replace with correct path */
          perform_splice(tmp_line,"url(gif/",full_link);

          /* Set flag */
          replaced = TRUE;
        }

      /* If GIF image is from the PDB entry's PDBsum directory, then
         set parameters to fetch it from there */
      else if (strstr(tmp_line,"src=\"tmpgif/") != NULL ||
               strstr(tmp_line,"src=tmpgif/") != NULL ||
               strstr(tmp_line,"href=\"tmpgif/") != NULL ||
               strstr(tmp_line,"background=\"tmpgif/") != NULL)
        {
          /* For ProFunc image, call GetImage.pl script */
          if (params->profunc == TRUE)
            {
              strcat(full_link,getprofunc_image_pl);
              if (!strcmp(params->pdb_code,"none") ||
                  !strcmp(params->pdb_code,"NONE"))
                {
                  strcat(full_link,"&u=");
                  strcat(full_link,params->user_id);
                }
              strcat(full_link,"&image=");
            }

          /* For LigSearch image, call GetImage.pl script */
          else if (params->ligsearch == TRUE)
            {
              strcat(full_link,getligsearch_image_pl);
              if (params->user_id != &null)
                {
                  strcat(full_link,"&temp=");
                  strcat(full_link,params->user_id);
                }
              strcat(full_link,"&source=ligsearch");
              strcat(full_link,"&file=");
            }

          /* For VarSite image, call GetImage.pl script */
          else if (params->varsite == TRUE)
            {
              strcat(full_link,getvarsite_image_pl);
              if (params->user_id != &null)
                {
                  strcat(full_link,"&temp=");
                  strcat(full_link,params->user_id);
                }
              strcat(full_link,"&source=varsite");
              strcat(full_link,"&file=");
            }

          /* For DrugPort image, link to appropriate subdirectory */
          else if (params->drugport == TRUE)
            {
              strcat(full_link,getimage_pl);
              strcat(full_link,"temp=drug");
              strcat(full_link,"&gif_dir=");
              strcat(full_link,params->subdir_dir);
              strcat(full_link,"&file=");
            }

          /* For PDBsum1, get the name of the directory */
          else if (params->pdb_type == PDBSUM1)
            {
              //              strcpy(full_link,"file://");
              //              strcat(full_link,pdbsum_dir);

              /* For relative links, point up a directory */
              //              if (params->relinks == TRUE)
                strcpy(full_link,"../");
            }

          /* If image in standard PDBsum directory, just refrerence
             it directly - no longer used. Now all images displayed
             via getimg.pl script */
          //else if (params->pdb_type == STANDARD_PDB ||
          //         params->pdb_type == ALPHA_FOLD)
          //{
          //strcat(full_link,pdbsum_dir);
          //  strcat(full_link,DIR_SEP);
          //}

          /* Otherwise need to retrieve the image on the fly using the
             getimg.pl script */
          else
            {
              strcat(full_link,getimage_pl);
              strcat(full_link,"source=pdbsum");
              if (params->pdb_type == UPLOAD_PDB)
                strcat(full_link,"&pdb_type=UPLOAD");
              else if (params->pdb_type == MCSG_PDB)
                strcat(full_link,"&pdb_type=MCSG");
              strcat(full_link,"&pdb_code=");
              strcat(full_link,params->pdb_code);
              if (params->security_code[0] != '\0')
                {
                  strcat(full_link,"&code=");
                  strcat(full_link,params->security_code);
                }
              strcat(full_link,"&file=");
            }

          /* Splice out current path and replace with correct path */
          perform_splice(tmp_line,"tmpgif/",full_link);

          /* Set flag */
          replaced = TRUE;
        }

      /* If have a link to a GIF in the PDBsum directory, then
         set full path */
      else if (strstr(tmp_line,"href=\"gif/") != NULL)
        {
          /* Form the full path to the GIF file */
          strcat(full_link,"href=\"");
          strcat(full_link,templates_gif_dir);
          strcat(full_link,DIR_SEP);

          /* For relative links, point to local directory */
          if (params->relinks == TRUE)
            {
              /* Set the link */
              strcpy(full_link,"href=\"../gif/");

              /* Extract the image name */
              nimage_names
                = extract_image_name(tmp_line,first_field_ptr,loop,
                                     &call_no,&on_off,counter,
                                     &first_iname_ptr,&last_iname_ptr,
                                     nimage_names,params);
            }

          /* Splice out current path and replace with correct path */
          perform_splice(tmp_line,"href=\"gif/",full_link);

          /* Set flag */
          replaced = TRUE;
        }
    }

  /* Return current pointers */
  *fst_iname_ptr = first_iname_ptr;
  *lst_iname_ptr = last_iname_ptr;

  /* Return the number of unique image names */
  return(nimage_names);
  
}
/***********************************************************************

process_hrefs  -  Replace all local HREF tags by appropriate CGI calls

***********************************************************************/

void process_hrefs(char *tmp_line,char *templates_dir,
                   char *getpage_pl,char *cgi_dir,int proprietary,
                   int pdb_type,int relinks,char *link_path)
{
  char full_link[LINELEN + 1];

  int replaced;

  /* Set flag indicating whether a tag has been replaced */
  replaced = TRUE;

  /* Loop until no more tags to be redirected */
  while (replaced == TRUE)
    {
      /* Unset flag */
      replaced = FALSE;
      full_link[0] = '\0';

      /* If line contains a HREF record pointing to the css directory,
         then redirect to the templates URL */
      if (strstr(tmp_line,"href=\"css/") != NULL)
        {
          /* Form the full path to the css directory */
          strcat(full_link,"href=\"");
          if (pdb_type == PDBSUM1)
            strcat(full_link,"file://");
          strcat(full_link,templates_dir);
          strcat(full_link,DIR_SEP);
          strcat(full_link,"css");
          strcat(full_link,DIR_SEP);

          /* For relative links, point to local css directory */
          if (relinks == TRUE)
	    {
	      /* Default links to the directory above */
	      strcpy(full_link,"href=\"../css/");

              /* If this is a common link,  keep the path but change
                 css to CSS */
              if (strstr(tmp_line," c=\"T\" ") != NULL)
		strcpy(full_link,"HREF=\"css/");
	    }
	  
	  /* Splice out current path and replace with correct path */
          perform_splice(tmp_line,"href=\"css/",full_link);

          /* Set flag */
          replaced = TRUE;
        }

      /* Check for a javascript directory */
      else if (strstr(tmp_line,"src=\"js/") != NULL)
        {
          /* Form the full path to the js directory */
          strcat(full_link,"src=\"");
          if (pdb_type == PDBSUM1)
            strcat(full_link,"file://");
          strcat(full_link,templates_dir);
          strcat(full_link,DIR_SEP);
          strcat(full_link,"js");
          strcat(full_link,DIR_SEP);

          /* For relative links, point to local js directory */
          if (relinks == TRUE)
            strcpy(full_link,"src=\"../js/");

          /* Splice out current path and replace with correct path */
          perform_splice(tmp_line,"src=\"js/",full_link);

          /* Set flag */
          replaced = TRUE;
        }

      /* Check for an iframe include file and replace by call to
         Getpage.pl CGI-script */
      else if (strstr(tmp_line,"iframe SRC=\"./") != NULL)
        {
          /* Form the full path to the GIF file */
          strcat(full_link,"src=\"");
          strcat(full_link,getpage_pl);

          /* Splice out current path and replace with correct path */
          perform_splice(tmp_line,"SRC=\"./",full_link);

          /* Set flag */
          replaced = TRUE;
        }

      /* Check for an iframe include file (old system) */
      else if (strstr(tmp_line,"src=\"./") != NULL)
        {
          /* Form the full path to the directory */
          strcat(full_link,"src=\"");
          strcat(full_link,templates_dir);
          strcat(full_link,DIR_SEP);

          /* Splice out current path and replace with correct path */
          perform_splice(tmp_line,"src=\"./",full_link);

          /* Set flag */
          replaced = TRUE;
        }

      /* If HREF calls one of the template pages, call it via the
         Getpage.pl CGI-script */
      else if (strstr(tmp_line,"href=\"./") != NULL)
        {
          /* Form the full path to the GIF file */
          strcat(full_link,"href=\"");
          strcat(full_link,getpage_pl);

          /* Splice out current path and replace with correct path */
          perform_splice(tmp_line,"href=\"./",full_link);

          /* If this is PDBsum Proprietary, splice in the flag */
          if (proprietary == TRUE && strstr(tmp_line,"?") != NULL)
            perform_splice(tmp_line,"?","?prop=TRUE&");

          /* Set flag */
          replaced = TRUE;
        }

      /* Check for call to script in cgi-bin directory */
      else if (strstr(tmp_line,"=\"/cgi-bin/") != NULL)
        {
          /* Form the full path to the GIF file */
          strcat(full_link,"=\"");
          strcat(full_link,cgi_dir);
          strcat(full_link,DIR_SEP);

          /* Splice out current path and replace with correct path */
          perform_splice(tmp_line,"=\"/cgi-bin/",full_link);

          /* If this is PDBsum Proprietary, splice in the flag */
          if (proprietary == TRUE && strstr(tmp_line,"?") != NULL)
            perform_splice(tmp_line,"?","?prop=TRUE&");

          /* Set flag */
          replaced = TRUE;
        }
    }
}
/***********************************************************************

process_ebi_incs  -  Replace EBI include paths with full URL

***********************************************************************/

void process_ebi_incs(char *tmp_line,char *ebi_include_dir)
{
  char full_link[LINELEN + 1];

  int replaced;

  /* Set flag indicating whether a tag has been replaced */
  replaced = TRUE;

  /* Loop until no more tags to be redirected */
  while (replaced == TRUE)
    {
      /* Unset flag */
      replaced = FALSE;
      full_link[0] = '\0';

      /* If line contains a HREF include dir, replace by full URL */
      if (strstr(tmp_line,"href=\"inc/") != NULL)
        {
          /* Form the full path to the EBI include directory */
          strcat(full_link,"href=\"");
          strcat(full_link,ebi_include_dir);
          strcat(full_link,DIR_SEP);

          /* Splice out current path and replace with correct path */
          perform_splice(tmp_line,"href=\"inc/",full_link);

          /* Set flag */
          replaced = TRUE;
        }

      /* Repeat for SRC refs */
      else if (strstr(tmp_line,"src=\"inc/") != NULL)
        {
          /* Form the full path to the EBI include directory */
          strcat(full_link,"src=\"");
          strcat(full_link,ebi_include_dir);
          strcat(full_link,DIR_SEP);

          /* Splice out current path and replace with correct path */
          perform_splice(tmp_line,"src=\"inc/",full_link);

          /* Set flag */
          replaced = TRUE;
        }
    }
}
/***********************************************************************

get_html_colour  -  Convert the given colour name to hexadecimal

***********************************************************************/

void get_html_colour(char *c_name,char *c_string)
{
#define NCOLS   15

  int colour, icol;

  static char *colour_name[] = { 
    "blue", "purple", "skyblue", "cyan", "green", "yellow",
    "orange", "pink", "magenta", "red", "khaki", "grey",
    "lime_green", "black", "white" };

  static char *colour_code[] = {
    "#0000FF", "#CC66FF", "#0099FF", "#33CCFF", "#00FF00", "#FFFF00",
    "#FF6600", "#FF69B4", "#FF00FF", "#FF0000", "#8B864E", "#555555",
    "#7DFF00", "#000000", "#FFFFFF" };

  /* Initialise colour to first colour */
  colour = 0;

  /* Search colours list for given colour name */
  for (icol = 0; icol < NCOLS; icol++)
    if (!strcmp(c_name,colour_name[icol]))
      colour = icol;

  /* Return the hexadecimal equivalent of the given colour */
  strcpy(c_string,colour_code[colour]);
}
/***********************************************************************

res_conservation  -  Read in the residue conservation data from the
                     appropriate  ConSurf file and output the
                     coloured sequence

***********************************************************************/

void res_conservation(FILE *outfile_ptr,char *code,char *consurf_dir)
{
  char file_name[FILENAME_LEN + 1], input_line[LINELEN + 1];
  char aa_code, c_name[20], c_string[8], middle[3];

  char *token[MAX_TOKEN];

  static char *colour_name[] = { 
    "black", "blue", "purple", "skyblue", "cyan", "green",
    "yellow", "orange", "pink", "red" };

  int conserv, last_conserv, nres, ntoken;
  int reading;

  FILE *file_ptr;

  /* Initialise variables */
  last_conserv = -1;
  nres = 0;
  reading = FALSE;

  /* Get middle part of PDB code */
  strncpy(middle,code + 1,2);
  middle[2] = '\0';

  /* Form name of the ConSurf data file */
  strcpy(file_name,consurf_dir);
  strcat(file_name,DIR_SEP);
  strcat(file_name,middle);
  strcat(file_name,DIR_SEP);
  strcat(file_name,code);
  strcat(file_name,".gradesPE");

  /* Open the ConSurf file */
  if ((file_ptr = fopen(file_name,"r")) !=  NULL)
    {
      /* Prepare to write out the annotated sequence */
      fprintf(outfile_ptr,"<PRE>");
      fprintf(outfile_ptr,"<B>");

      /* Loop while reading in records from the file */
      while (fgets(input_line,LINELEN,file_ptr) != NULL)
        {
          /* If reading in the residue data, process this line */
          if (reading == TRUE)
            {
              /* Retrieve all the tokens from the current line */
              ntoken = get_tokens(input_line,token,MAX_TOKEN,'\t');

              /* If have sufficient data, then store details */
              if (ntoken > 6)
                {
                  /* Get the residue's single-letter a.a. code */
                  aa_code = token[1][0];

                  /* Get the conservation score */
                  conserv = atoi(token[4]);

                  /* If conservation score has changed, write out new
                     colour */
                  if (conserv != last_conserv)
                    {
                      /* If have a previous colour, close off its font */
                      if (last_conserv > -1)
                        fprintf(outfile_ptr,"</FONT>");

                      /* Get associated colour */
                      strcpy(c_name,colour_name[conserv]);
                      get_html_colour(c_name,c_string);

                      /* Set the required colour */
                      fprintf(outfile_ptr,"<FONT COLOR=\"%s\">",c_string);

                      /* Save the current conservation score */
                      last_conserv = conserv;
                    }

                  /* Write out this residue */
                  fprintf(outfile_ptr,"%c",aa_code);

                  /* Increment residue count */
                  nres++;

                  /* If have a line-full, throw line */
                  if (nres % RES_PER_LINE == 0)
                    fprintf(outfile_ptr,"\n");
                }
            }

          /* Otherwise check for column titles line signifying the start
             of the data */
          else if (!strncmp(input_line," POS",4))
            reading = TRUE;
        }

      /* If have a previous colour, close off its font */
      if (last_conserv > -1)
        fprintf(outfile_ptr,"</FONT>");

      /* Close off the sequence */
      fprintf(outfile_ptr,"</B>");
      fprintf(outfile_ptr,"</PRE>\n");

      /* Close the input file */
      fclose(file_ptr);
    }
}
/***********************************************************************

read_offset_index  -  Read in the offset value for the given key field
                      in the given index file

***********************************************************************/

int read_offset_index(struct field *first_field_ptr,char *file_name,
                      char *field_name,char *key,
                      struct parameters *params)
{
  char number_string[NAMLEN + 1];
  char input_line[LINELEN + 1], string[LINELEN + 1], value[LINELEN + 1];

  int key_len;
  int found;

  struct field *field_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  found = FALSE;
  key_len = strlen(key);
  strcpy(value,"-1");

  /* Initialise field pointers */
  field_ptr = NULL;

  /* Open file to check whether it exists */
  if ((file_ptr = fopen(file_name,"r")) !=  NULL)
    {
      /* Loop while reading through the file */
      while (fgets(input_line,LINELEN,file_ptr) != NULL &&
             found == FALSE)
        {
          /* Test for key */
          if (!strncmp(input_line,key,key_len))
            {
              /* Truncate the string to strip off the newline */
              string_truncate(input_line,LINELEN);

              /* Store the field value */
              strcpy(value,input_line+key_len);

              /* Set flag that record found */
              found = TRUE;
            }
        }

      /* Close the file */
      fclose(file_ptr);
    }

  /* Create a field storing the offset index */
  sprintf(string,"%s %s",field_name,value);

  /* Create a field record */
  field_ptr
    = create_field_record(first_field_ptr,string,FALSE,FALSE,LOOP,-1);

  /* Return whether string was found */
  return(found);
}
/***********************************************************************

get_image_names  -  Read in any image names that may be in relinks.txt

***********************************************************************/

int get_image_names(struct iname **fst_iname_ptr,
                    struct iname **lst_iname_ptr)
{
  char file_name[FILENAME_LEN], input_line[LINELEN + 1];

  int len, nimage_names;

  struct iname *iname_ptr, *first_iname_ptr, *last_iname_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  nimage_names = 0;

  /* Initialise pointers */
  iname_ptr = first_iname_ptr = last_iname_ptr = NULL;

  /* Form the name of the input file */
  strcpy(file_name,"relinks.txt");

  /* Open the file */
  file_ptr = fopen(file_name,"r");

  /* If file doesn't exist, then return */
  if (file_ptr == NULL)
    return(nimage_names);

  /* Loop through the file, line by line */
  while (fgets(input_line,LINELEN,file_ptr) != NULL)
    {
      /* Truncate the input line */
      len = string_truncate(input_line,LINELEN);

      /* Save this name */
      iname_ptr
        = create_iname_record(&first_iname_ptr,&last_iname_ptr,
                              input_line);

      /* Increment count of unique image names */
      nimage_names++;
    }

  /* Return current pointers */
  *fst_iname_ptr = first_iname_ptr;
  *lst_iname_ptr = last_iname_ptr;

  /* Return number of images read in */
  return(nimage_names);
}
/***********************************************************************

write_html_file  -  Annotate and write out the HTML file

***********************************************************************/

void write_html_file(char **html_line,int nlines,
                     struct field *first_field_ptr,
                     struct field *last_field_ptr,
                     struct iname **fst_iname_ptr,
                     struct parameters *params)
{
  char file_name[FILENAME_LEN + 1], ftype[6], string[LINELEN + 1];
  char loop_end[NAMLEN + 1], end_skip[NAMLEN + 1];
  char field_name[NAMLEN + 1], fname[NAMLEN + 1], tmp_line[LINELEN + 1];
  char replace_from[LINELEN + 1], replace_by[LINELEN + 1];
  char consurf_dir[FILENAME_LEN + 1], key[LINELEN + 1];
  char pdbsum_dir[FILENAME_LEN + 1], pdbsum_real_dir[FILENAME_LEN + 1];
  char templates_gif_dir[FILENAME_LEN + 1];
  char getpage_pl[FILENAME_LEN + 1], getimage_pl[FILENAME_LEN + 1];
  char getprofunc_image_pl[FILENAME_LEN + 1];
  char profunc_results_dir[FILENAME_LEN + 1];
  char getligsearch_image_pl[FILENAME_LEN + 1];
  char ligsearch_results_dir[FILENAME_LEN + 1];
  char getvarsite_image_pl[FILENAME_LEN + 1];
  char varsite_results_dir[FILENAME_LEN + 1];
  char ebi_include_dir[FILENAME_LEN + 1];
  char templates_dir[FILENAME_LEN + 1];
  char cgi_dir[FILENAME_LEN + 1], number_string[20];

  char *string_ptr, *string1_ptr, *string2_ptr;

  int call_no, icount, ivalue, line, on_off;
  int current_level, ilevel, len, line_pos, level, spos;
  int loop[MAX_LEVEL], loop_start[MAX_LEVEL], nloops[MAX_LEVEL];
  int ndata, nimage_names,nrecords, start_record;
  int dbt_format, have_value, promotif, skip, skip_loop, valid, wanted;
  int gloop, have_oread, have_pread;

  long counter[NCOUNTERS];

  static long offset = -1;

  struct field *field_ptr;
  struct iname *first_iname_ptr, *last_iname_ptr;

  FILE *outfile_ptr;

  /* Initialise variables */
  call_no = 0;
  for (icount = 0; icount < NCOUNTERS; icount++)
    counter[icount] = 0;
  loop_end[0] = '\0';
  end_skip[0] = '\0';
  fname[0] = '\0';
  gloop = FALSE;
  line = 0;
  nimage_names = 0;
  nrecords = -1;
  start_record = 0;
  on_off = 0;
  replace_from[0] = replace_by[0] = '\0';
  current_level = -1;
  level = -1;
  for (ilevel = 0; ilevel < MAX_LEVEL; ilevel++)
    {
      loop[ilevel] = 0;
      nloops[ilevel] = 0;
      loop_start[ilevel] = -1;
    }
  skip = FALSE;
  skip_loop = FALSE;

  /* Initialise image pointers */
  first_iname_ptr = last_iname_ptr = NULL;

  /* Get the directory names required for the process_gif_names
     routine */
  get_directories(first_field_ptr,templates_gif_dir,pdbsum_dir,
                  pdbsum_real_dir,templates_dir,getpage_pl,getimage_pl,
                  getprofunc_image_pl,getligsearch_image_pl,
                  getvarsite_image_pl,cgi_dir,consurf_dir,
                  profunc_results_dir,ligsearch_results_dir,
                  varsite_results_dir,ebi_include_dir,params);

  /* Form name of the output file */
  if (params->output_file[0] != '\0')
    {
      strcpy(file_name,params->output_file);

      /* Open output HTML file */
      outfile_ptr = open_output_file(file_name,TRUE);
    }
  else
    outfile_ptr = stdout;

  /* If relative links are required, read in that may have been written
     out in an earlier run */
  if (params->relinks == TRUE)
    nimage_names = get_image_names(&first_iname_ptr,&last_iname_ptr);

  /* Loop through all the lines of the HTML file, annotating and
     writing out */
  while (line < nlines)
    {
      /* Get the current line */
      strcpy(tmp_line,html_line[line]);

      /* Assume it is wanted */
      wanted = TRUE;

      /* Check for start of block that might need to be skipped */
      if ((string_ptr = strstr(tmp_line,"<SKIP ")) != NULL &&
          skip == FALSE)
        {
          /* Process line to replace any tags by their appropriate
             substitutions */
          if (strstr(tmp_line,TAG_END) != NULL)
            valid
              = annotate_line(tmp_line,first_field_ptr,loop,&call_no,
                              &on_off,counter,params);

          /* If section to be skipped if given file not found, test
             whether file exists */
          if (!strncmp(string_ptr+6,"\"",1) ||
              !strncmp(string_ptr+6,"!\"",2))
            skip = check_for_file(string_ptr,end_skip);

          /* Check for presences of given string in file */
          else if (!strncmp(string_ptr+6,"[",1) ||
              !strncmp(string_ptr+6,"![",2))
            {
              /* Check if string is present in file and retrieve its
                 position */
              skip = check_in_file(string_ptr,end_skip,&line_pos);

              /* Create/update the FILE_INDEX_POS field */
              strcpy(field_name,"FILE_INDEX_POS");

              /* From the string containing keyword-value pair */
              sprintf(string,"%s %d",field_name,line_pos);

              /* Create a field record */
              field_ptr
                = create_field_record(first_field_ptr,string,FALSE,
                                      FALSE,LOOP,-1);
            }

          /* Otherwise, check for presence of the current field */
          else
            skip = check_field(first_field_ptr,"<SKIP ",tmp_line,
                               loop,end_skip,&ivalue,ftype,params);

          /* Set flag that this line itself is not wanted for output */
          wanted = FALSE;
        }

      /* Check for end of skipped block */
      else if ((string_ptr = strstr(tmp_line,"</SKIP ")) != NULL)
        {
          /* Process line to replace any tags by their appropriate
             substitutions */
          if (strstr(tmp_line,TAG_END) != NULL)
            valid
              = annotate_line(tmp_line,first_field_ptr,loop,&call_no,
                              &on_off,counter,params);

          /* If this is the end of the current skip, unset the skip
             flag */
          if (strstr(tmp_line,end_skip) != NULL)
            skip = FALSE;

          /* Set flag that this line itself is not wanted for output */
          wanted = FALSE;
        }

      /* If skipping, then line is not wanted */
      if (skip == TRUE || skip_loop == TRUE)
        wanted = FALSE;

      /* If not skipping, process current line */
      if (skip == FALSE)
        {
          /* Check for start of block of HTML that is to be repeated for
             each sequence */
          if (((string1_ptr = strstr(tmp_line,"<LOOP ")) != NULL ||
               (string2_ptr = strstr(tmp_line,"<GLOOP ")) != NULL) &&
              skip_loop == FALSE)
            {
              /* Check for "general" loop */
              if (strstr(tmp_line,"<GLOOP ") != NULL)
                {
                  gloop = TRUE;
                  string_ptr = string2_ptr;
                }
              else
                {
                  gloop = FALSE;
                  string_ptr = string1_ptr;
                }

              /* Process line to replace any tags by their appropriate
                 substitutions */
              if (strstr(tmp_line,TAG_END) != NULL)
                valid
                  = annotate_line(tmp_line,first_field_ptr,loop,&call_no,
                                  &on_off,counter,params);

              /* Set flags */
              wanted = FALSE;

              /* Get the field-name associated with the loop */
              if (gloop == TRUE)
                strcpy(field_name,string_ptr + 7);
              else
                strcpy(field_name,string_ptr + 6);
              string_ptr = strstr(field_name,">");
              if (string_ptr != NULL)
                {
                  /* Initialise flags */
                  have_value = FALSE;

                  /* Terminate field name */
                  string_ptr[0] = '\0';

                  /* Identify this field */
                  field_ptr = get_field(first_field_ptr,field_name,loop,
                                        params->instance);

                  /* If valid, get the number of loops required */
                  if (field_ptr != NULL)
                    {
                      /* Get the level at which this loop operates */
                      level = field_ptr->nlevel;

                      /* If field's level is zero, but we are higher
                         than this, set it one above the current level */
                      if (gloop == TRUE && current_level > -1)
                        level = current_level + 1;

                      /* Check that level is a valid one */
                      if (level < 0 || level > MAX_LEVEL - 1)
                        {
                          printf("*** TEMPLATE ERROR. Field %s has too "
                                 "many levels!\n",field_name);
                          printf("    Line %d of template %s\n",line,
                                 params->template_name);
                          exit(1);
                        }

                      /* Store the value */
                      current_level = level;
                      ivalue = atoi(field_ptr->value);

                      /* Set flag */
                      have_value = TRUE;
                    }

                  /* Check for special loop-counter field */
                  else if (!strncmp(field_name,"_LOOP",5))
                    {
                      /* Determine which loop counter this is */
                      strcpy(number_string,field_name+5);
                      ilevel = atoi(number_string);

                      /* If a valid level, store the value */
                      if (ilevel > 0 && ilevel < MAX_LEVEL + 1)
                        if (loop[ilevel - 1] > 0)
                          ivalue = loop[ilevel - 1];
                        else
                          ivalue = params->instance[ilevel - 1];
                      else
                        ivalue = 1;

                      /* Set flag */
                      have_value = TRUE;
                      level++;
                    }

                  /* If have a value, then set up this loop */
                  if (have_value == TRUE)
                    {
                      /* Set number of loop at this level */
                      nloops[level] = ivalue;

                      /* Initialise loop number */
                      loop[level] = 1;

                      /* Store block's start-point */
                      loop_start[level] = line;

                      /* If have only a single loop, then this first loop is
                         also the end loop, otherwise we're in midloop */
                      end_loop(first_field_ptr,level,loop[level],
                               nloops[level],params);
                    }

                  /* If field not found, then need to skip this loop
                     entirely */
                  else
                    {
                      /* Set flag */
                      skip_loop = TRUE;
                    }
                }

              /* If field-name not found, then want to skip the loop */
              else
                skip_loop = TRUE;

              /* If skipping this loop, form the string marking the
                 end of the loop */
              if (skip_loop == TRUE)
                {
                  if (gloop == TRUE)
                    {
                      strcpy(loop_end,"</GLOOP ");
                      strcat(loop_end,tmp_line+7);
                    }
                  else
                    {
                      strcpy(loop_end,"</LOOP ");
                      strcat(loop_end,tmp_line+6);
                    }
                  string_ptr = strstr(loop_end,">");
                  if (string_ptr != NULL)
                    string_ptr[1] = '\0';
                }
            }

          /* Check for end of block */
          else if (strstr(tmp_line,"</LOOP ") != NULL ||
                   strstr(tmp_line,"</GLOOP ") != NULL)
            {
              /* Process line to replace any tags by their appropriate
                 substitutions */
              if (strstr(tmp_line,TAG_END) != NULL)
                valid
                  = annotate_line(tmp_line,first_field_ptr,loop,&call_no,
                                  &on_off,counter,params);

              /* If skipping one or more loops, check whether this is
                 the end of our skip-loop */
              if (skip_loop == TRUE && strstr(tmp_line,loop_end) != NULL)
                {
                  /* Have the end of the loop we've been skipping,
                     so unset flag */
                  skip_loop = FALSE;
                  loop_end[0] = '\0';
                }

              /* Otherwise, must be the end of one of our loops, so
                 need to increment loop counter and go back to start
                 of loop */
              else if (skip_loop == FALSE)
                {
                  /* Check for error in template */
                  if (level < 0)
                    {
                      printf("*** TEMPLATE ERROR. Loop ends before its "
                             "start!\n");
                      printf("    Line %d of template %s\n",line,
                             params->template_name);
                      printf("    Line: %s\n",tmp_line);
                      //list_fields(first_field_ptr);
                      exit(1);
                    }

                  /* If all loops haven't yet been done, increment loop
                     counter and loop back */
                  loop[level]++;
                  if (loop[level] < nloops[level] + 1)
                    {
                      /* Return to block's start-point */
                      line = loop_start[level];

                      /* Check for loop end */
                      end_loop(first_field_ptr,level,loop[level],
                               nloops[level],params);
                    }

                  /* Otherwise, reset loop number and reduce level */
                  else
                    {
                      /* Reset */
                      loop[level] = 0;
                      skip_loop = FALSE;
                      level--;
                      current_level = level;
                    }
                }

              /* This line not wanted */
              wanted = FALSE;
            }

          /* Check for any of the tags that need to be replaced */
          else if (level < 0 || skip_loop == FALSE)
            {
              /* Check for any GIF image names and replace by CGI script
                 to retrieve the images */
              nimage_names
                = process_gif_names(tmp_line,templates_gif_dir,pdbsum_dir,
                                    getimage_pl,getprofunc_image_pl,
                                    getligsearch_image_pl,
                                    getvarsite_image_pl,first_field_ptr,
                                    loop,&call_no,&on_off,counter,
                                    &first_iname_ptr,&last_iname_ptr,
                                    nimage_names,params);

              /* Check for any HREF lines referring to template files */
              process_hrefs(tmp_line,templates_dir,getpage_pl,cgi_dir,
                            params->proprietary,params->pdb_type,
                            params->relinks,params->link_path);

              /* Check for any HREF lines referring to template files */
              process_ebi_incs(tmp_line,ebi_include_dir);

              /* If line refers to include file, perform the include */
              if (!strncmp(tmp_line,"#include:",9))
                {
                  /* Process line to replace all tags by their appropriate
                     substitutions */
                  valid
                    = annotate_line(tmp_line,first_field_ptr,loop,&call_no,
                                    &on_off,counter,params);

                  /* Pick up the include file and insert into the output
                     file here */
                  insert_file(outfile_ptr,tmp_line+9,pdbsum_real_dir,
                              pdbsum_dir,getpage_pl,getimage_pl,
                              templates_gif_dir,profunc_results_dir,
                              ligsearch_results_dir,
                              varsite_results_dir,first_field_ptr,
                              loop,&call_no,&on_off,counter,params);

                  /* Set flag to exclude current line */
                  wanted = FALSE;
                }

              /* Special residue conservation command calling for read-in
                 of residue conservation data from the appropriate 
                 ConSurf file and insertion of annotated sequence */
              else if (!strncmp(tmp_line,"#rescons:",9))
                {
                  /* Process line to replace all tags by their appropriate
                     substitutions */
                  valid
                    = annotate_line(tmp_line,first_field_ptr,loop,&call_no,
                                    &on_off,counter,params);

                  /* Read in the residue conservation data and output
                     the coloured sequence */
                  res_conservation(outfile_ptr,tmp_line+9,consurf_dir);

                  /* Set flag to exclude current line */
                  wanted = FALSE;
                }

              /* If file containing more parameters to be read in, pick
                 up and store the additional data */
              else if (!strncmp(tmp_line,"#read:",6) ||
                       !strncmp(tmp_line,"#oread:",7) ||
                       !strncmp(tmp_line,"#pread:",7) ||
                       !strncmp(tmp_line,"#kread:",7) ||
                       !strncmp(tmp_line,"#kpread:",8))
                {
                  /* Initialise flags */
                  have_oread = FALSE;
                  have_pread = FALSE;

                  /* Initialise replacement strings */
                  replace_from[0] = replace_by[0] = '\0';

                  /* Process line to replace all tags by their appropriate
                     substitutions */
                  valid
                    = annotate_line(tmp_line,first_field_ptr,loop,&call_no,
                                    &on_off,counter,params);
                  key[0] = '\0';
                  offset = -1;

                  /* Check for replacement command */
                  string_ptr = strstr(tmp_line,"::REPLACE:");
                  if (string_ptr != NULL)
                    {
                      /* Save the replacement string and terminate
                         filename */
                      strcpy(replace_from,string_ptr+10);
                      string_ptr[0] = '\0';

                      /* Split the string into the part to be replaced
                         and the part to replace by */
                      string_ptr = strchr(replace_from,':');
                      if (string_ptr != NULL)
                        {
                          /* Save the replacement string and terminate
                             string to be replaced */
                          strcpy(replace_by,string_ptr+1);
                          string_ptr[0] = '\0';
                        }
                    }

                  /* Determine whether this is a pread command requiring
                     that all fields read in be prefixed by the current
                     loop */
                  if (!strncmp(tmp_line,"#pread:",7))
                    {
                      have_pread = TRUE;
                      spos = 7;
                    }

                  /* If reading data for specified key value only,
                     separate the key field from the file name */
                  else if (!strncmp(tmp_line,"#kread:",7) ||
                           !strncmp(tmp_line,"#kpread:",8) ||
                           !strncmp(tmp_line,"#oread:",7))
                    {
                      /* Initialise variables */
                      nrecords = -1;
                      start_record = 0;

                      /* If oread command, set flag */
                      if (!strncmp(tmp_line,"#oread:",7))
                        have_oread = TRUE;

                      /* Initialise the field-name */
                      strcpy(field_name,"OFFSET");

                      /* Define start-position for file name */
                      have_pread = FALSE;
                      spos = 7;
                      if (!strncmp(tmp_line,"#kpread:",8))
                        {
                          spos = 8;
                          have_pread = TRUE;
                        }

                      /* Check for FIELD_NAME record associated with
                         the oread command */
                      string_ptr = strstr(tmp_line,"::FIELD_NAME:");
                      if (string_ptr != NULL)
                        {
                          /* Get the field_name */
                          strcpy(field_name,string_ptr+13);
                          string_ptr[0] = '\0';
                        }

                      /* Check for FOR record */
                      string_ptr = strstr(tmp_line,"::FOR:");
                      if (string_ptr != NULL)
                        {
                          /* Get the number of records to read in */
                          strcpy(key,string_ptr+6);
                          string_ptr[0] = '\0';
                          nrecords = atoi(key);
                          key[0] = '\0';
                        }

                      /* Check for START record */
                      string_ptr = strstr(tmp_line,"::START:");
                      if (string_ptr != NULL)
                        {
                          /* Get the starting record to read in */
                          strcpy(key,string_ptr+8);
                          string_ptr[0] = '\0';
                          start_record = atoi(key);
                          if (start_record == 0)
                            start_record++;
                          key[0] = '\0';
                        }

                      /* Look for key to search on */
                      string_ptr = strstr(tmp_line,"::KEY:");
                      if (string_ptr != NULL)
                        {
                          /* Save the key and terminate filename */
                          strcpy(key,string_ptr+6);
                          string_ptr[0] = '\0';

                          /* If oread, append a space to the key */
                          if (have_oread == TRUE)
                            strcat(key," ");
                        }

                      /* If no KEY, check whether we have an OFFSET */
                      else
                        {
                          /* Look for file offset */
                          string_ptr = strstr(tmp_line,"::OFFSET:");
                          if (string_ptr != NULL)
                            {
                              /* Save the offset and terminate filename */
                              strcpy(number_string,string_ptr+9);
                              offset = atol(number_string);
                              string_ptr[0] = '\0';
                            }
                        }
                    }

                  /* Otherwise reading standard .dat file */
                  else
                    {
                      have_pread = FALSE;
                      spos = 6;
                    }

                  /* Check for Victor-style double-colon format */
                  dbt_format = FALSE;
                  promotif = FALSE;
                  len = strlen(tmp_line);
                  if (len > 4 && !strcmp(tmp_line+len-4,".dbt"))
                    dbt_format = TRUE;

                  /* Check for promotif.dat file which has the double-colon
                     format, plus header records giving field names */
                  if (strstr(tmp_line,"promotif.dat") != NULL)
                    {
                      dbt_format = TRUE;
                      promotif = TRUE;
                    }

                  /* If reading in offset from an index file, get it
                     for the current key */
                  if (have_oread == TRUE)
                    ndata
                      = read_offset_index(first_field_ptr,tmp_line+spos,
                                          field_name,key,params);

                  /* If have standard .dat file, read in data */
                  else if (dbt_format == FALSE)
                    ndata
                      = read_database_file(first_field_ptr,
                                           tmp_line+spos,TRUE,TRUE,
                                           have_pread,key,offset,
                                           start_record,nrecords,loop,
                                           level,replace_from,replace_by,
                                           params);
                  else
                    ndata
                      = read_dbt_file(first_field_ptr,tmp_line+spos,TRUE,
                                      promotif,offset,nrecords,loop,params);

                  /* Update the special READ_OK field */
                  field_ptr
                    = get_field(first_field_ptr,"READ_OK",LOOP,
                                params->instance);

                  /* Update the field's value */
                  if (field_ptr != NULL)
                    {
                      /* If have an existing value, clear it out */
                      if (field_ptr->value != &null)
                        free(field_ptr->value);

                      /* Store read-status */
                      if (ndata > 0)
                        store_string(&(field_ptr->value),"TRUE");
                      else
                        store_string(&(field_ptr->value),"FALSE");
                    }

                  /* Set flag to exclude current line */
                  wanted = FALSE;
                }

              /* Process line to replace all tags by their appropriate
                 substitutions */
              if (strstr(tmp_line,TAG_END) != NULL)
                valid
                  = annotate_line(tmp_line,first_field_ptr,loop,&call_no,
                                  &on_off,counter,params);

              /* Fix any plural terms that need fixing */
              fix_plurals(tmp_line,first_field_ptr,loop,params);
            }
        }

      /* Check for blank lines or comment */
      if (wanted == TRUE)
        {
          /* If line is blank, don't need to write it out */
          if (tmp_line[0] == '\0')
            wanted = FALSE;

          /* Otherwise, check whether line is a comment line */
          else if (!strncmp(tmp_line,"<!--",4) &&
                   strstr(tmp_line,"-->") != NULL)
            wanted = FALSE;
        }

      /* If line is wanted, write out to file */
      if (wanted == TRUE)
        {
          /* Check for a "no line throw" marker at then end of the
             line */
          string_ptr = strstr(tmp_line,"<nlt>");
          if (string_ptr != NULL)
            {
              /* Remove the marker and print line without a line-throw */
              string_ptr[0] = '\0';
              fprintf(outfile_ptr,"%s",tmp_line);
            }

          /* Otherwise, write to file as normal */
          else
            fprintf(outfile_ptr,"%s\n",tmp_line);
        }

      /* Increment line count */
      line++;
    }

  /* Close the output file */
  fclose(outfile_ptr);

  /* Return the first of the image names */
  *fst_iname_ptr = first_iname_ptr;
}
/***********************************************************************

write_image_names  -  Write out the image names

***********************************************************************/

void write_image_names(struct iname *first_iname_ptr)
{
  struct iname *iname_ptr;

  FILE *file_ptr;

  /* Open the output file */
  file_ptr = open_output_file("relinks.txt",TRUE);

  /* Get pointer to the first image_name */
  iname_ptr = first_iname_ptr;

  /* Loop over all the image_names to write out */
  while (iname_ptr != NULL)
    {
      /* Write this name out */
      fprintf(file_ptr,"%s\n",iname_ptr->image_name);

      /* Get pointer to the next image_name */
      iname_ptr = iname_ptr->next_iname_ptr;
    }

  /* Close the output file */
  fclose(file_ptr);
}

/***********************************************************************

                             M   A   I   N

***********************************************************************/

int main(int argc,char *argv[])
{
  char string[FILENAME_LEN + 1], tmp[FILENAME_LEN + 1];
  char file_name[FILENAME_LEN + 1];

  char *html_line[MAX_HTML_LINES], key[2];

  char replace_from[2] = { "" };
  char replace_by[2] = { "" };

  int ifile, loop[MAX_LEVEL], level, nimage_names, nlines;

  static const int start_record = 0;
  static const int nrecords = -1;

  static const long offset = -1;

  struct parameters params;

  struct field *field_ptr, *first_field_ptr, *last_field_ptr;
  struct iname *first_iname_ptr, *last_iname_ptr;

#define NDAT_FILES              9

  static char *DAT_FILE[NDAT_FILES]
    = {
    "motifs.dat",
    "ligands.dat",
    "interfaces.dat",
    "seqinfo.dat",
    "molecules.dat",
    "variants.dat",
    "*_mole_params.dat",
    "*_mole2_params.dat",
    "*_mole2_pore_params.dat"
  };

  static char *DAT_DIR[NDAT_FILES]
    = {
    "promotif",
    "newsec",
    "newsec",
    ".",
    ".",
    ".",
    ".",
    ".",
    "."
  };

  FILE *file_ptr;

  /* Initialise variables */
  first_field_ptr = last_field_ptr = NULL;
  first_iname_ptr = last_iname_ptr = NULL;
  key[0] = '\0';
  nimage_names = 0;

  /* Initialise binary tree holding field records */
  first_field_ptr = create_field("");

  /* Read in the command-line parameters and store */
  get_command_arguments(argv,argc - 1,&params);

  /* Set the loop level according to command-line parameters */
  level = set_loop_level(loop,&params);

  /* Get filenames, paths and URLs from the CATHPARAM file */
  get_all_paths(first_field_ptr);

  /* Check for proprietary PDBsum */
  field_ptr = get_field(first_field_ptr,"EBI",LOOP,params.instance);
  if (field_ptr != NULL && !strcmp(field_ptr->value,"FALSE"))
    params.proprietary = TRUE;

  /* Create PDIR field, the name of the PDBsum directory for this PDB code */
  create_pdir(first_field_ptr,string,tmp,&params);

  /* Get the name of the corresponding PDB file and create a PDB_FILE
     record */
  get_pdb_filename(first_field_ptr,string,file_name,params);

  /* Repeat for VDIR, the name of the variants directory for this PDB code */
  create_vdir(first_field_ptr,string,tmp,&params);

  /* If not displaying an error page, pick up the data from the
     relevant database file(s) */
  if (params.have_error == FALSE)
    {
      /* If have a database file given as a command-line option,
         read it in */
      if (params.database_dat[0] != '\0')
        {
          /* Read in all the data fields from the given database.dat
             file */
          read_database_file(first_field_ptr,params.database_dat,TRUE,
                             FALSE,FALSE,key,offset,start_record,nrecords,
                             loop,-1,replace_from,replace_by,&params);
        }

      /* If showing pages for Enzyme Structures Database or the
         PDBsum highlights page, get the appropriate database.dat file */
      else if (params.enzymes == TRUE || params.highlights == TRUE)
        {
          /* Read in the database.dat file */
          read_database_file(first_field_ptr,"database.dat",FALSE,FALSE,
                             FALSE,key,offset,start_record,nrecords,loop,
                             -1,replace_from,replace_by,&params);

          /* Read in the paths.dat file */
          read_database_file(first_field_ptr,"paths.dat",FALSE,FALSE,
                             FALSE,key,offset,start_record,nrecords,loop,
                             -1,replace_from,replace_by,&params);
        }

      /* Otherwise, read in the data for the given PDB code */
      else if (params.drat == FALSE && params.sas == FALSE)
        {
          /* Read in all the data fields for this PDB entry from the 
             database.dat file in its PDBsum directory */
          read_database_file(first_field_ptr,"database.dat",FALSE,FALSE,
                             FALSE,key,offset,start_record,nrecords,loop,
                             -1,replace_from,replace_by,&params);


          /* For PDBsum, read in the remaining data files */
          if (params.profunc == FALSE && params.drugport == FALSE &&
              params.ligsearch == FALSE && params.varsite == FALSE &&
              params.ligplus == FALSE)
            {
              /* Loop over all the .dat files to be read in */
              for (ifile = 0; ifile < NDAT_FILES; ifile++)
                {
                  /* If PDBsum1, add relevant subdirectory */
                  file_name[0] = '\0';
                  if (params.pdb_type == PDBSUM1)
                    {
                      strcpy(file_name,DAT_DIR[ifile]);
                      strcat(file_name,"/");
                    }

                  /* If have an asterisk at the start, replace by the
                     PDB code */
                  if (DAT_FILE[ifile][0] == '*')
                    {
                      /* Prefix with the PDB code */
                      strcat(file_name,params.pdb_code);
                      strcat(file_name,DAT_FILE[ifile] + 1);
                    }
                  else
                    strcat(file_name,DAT_FILE[ifile]);

                  /* Read in the data for this .dat file */
                  read_database_file(first_field_ptr,file_name,
                                     FALSE,FALSE,FALSE,key,offset,
                                     start_record,nrecords,loop,-1,
                                     replace_from,replace_by,&params);
                }
            }
        }
    }

  /* Use the database data to define certain special flags used in the
     template */
  define_specials(first_field_ptr,&params);

  /* Read in the whole template file into memory */
  nlines = get_html_template(html_line,first_field_ptr,loop,&params);

  /* Write out the HTML file */
  write_html_file(html_line,nlines,first_field_ptr,last_field_ptr,
                  &first_iname_ptr,&params);

  /* Write out the image names, if required */
  if (params.relinks == TRUE)
    write_image_names(first_iname_ptr);

  /* For debugging, list the fields */
  if (params.list_fields == TRUE)
    {
      /* Open output file, showmain.list */
      file_ptr = open_output_file("showmain.list",TRUE);

      /* Write the fields out */
      list_fields(file_ptr,first_field_ptr);

      /* Close the file */
      fclose(file_ptr);
    }
}

