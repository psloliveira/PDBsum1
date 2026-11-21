/***********************************************************************

  pdbsum1.c - Program to generate the PDBsum1 pages for a single PDB
              file, or a list of files

************************************************************************

Date:-         18 Feb 2022

Written by:-   Roman Laskowski

Available PDB codes: 726,329

Version 1.1 - 30 Jan 2026



v.1.0    Original version                                  18 Feb 2022

v.1.0.1  Change order of commands for convert/magick as the latter needs
         them to go after the name of the input image.
                                                     30 Jan 2026 (RAL)

----------------------------------------------------------------------*/

/* Datafiles
   ---------

Actual name           Description
-----------           -----------
dataset.lst           Input list of PDB files to be processed
[filename]            Input PDB file
runs.dbt              Output list of all past PDBsum1 runs

Compilation
-----------

cc -c pdbsum1.c
cc -c common.c
cc -c reapar.c
cc -c readpdb.c
cc -o pdbsum1 pdbsum1.o common.o reapar.o readpdb.o -lm -lz

Optimization flags:-  

cc -O2 -funroll-loops -c pdbsum1.c
cc -O2 -funroll-loops -c common.c
cc -O2 -funroll-loops -c reapar.c
cc -O2 -funroll-loops -c readpdb.c
cc -o pdbsum1 pdbsum1.o common.o reapar.o readpdb.o -lm -lz

*/

/* Function calling-tree
   ---------------------

main
   -> get_command_arguments
         -> store_string (in common.c)
   -> c_get_paths (in common.c)
         -> get_param (in reapar.c)
         -> store_string (in common.c)
   -> initialise_run
         -> open_output_file (in common.c)
         -> create_directory (in common.c)
         -> check_file (in common.c)
   -> read_dataset_lst
         -> string_chomp (in common.c)
         -> get_tokens (in common.c)
         -> extract_pdb_code
               -> check_file (in common.c)
               -> gen_pdb_code
         -> create_pdb_entry
               -> store_string (in common.c)
         -> free_tokens (in common.c)
         -> process_file
               -> initialise
                     -> create_directory (in common.c)
               -> get_pdb_details
                     -> string_chomp (in common.c)
               -> create_directory (in common.c)
               -> get_date_and_time (in common.c)
               -> run_fixpdb
                     -> create_directory (in common.c)
               -> run_clean
                     -> create_directory (in common.c)
                     -> write_int_file
                           -> open_output_file (in common.c)
               -> run_procheck
                     -> list_plots
                           -> open_output_file (in common.c)
                           -> string_chomp (in common.c)
                           -> check_file (in common.c)
               -> run_hbplus
                     -> create_directory (in common.c)
               -> run_newsec
                     -> create_directory (in common.c)
               -> wait_for_file
                     -> replace_char (in common.c)
                     -> check_file (in common.c)
               -> get_interfaces
                     -> string_truncate (in common.c)
                     -> find_chain
                     -> create_chain_record
                     -> create_interface_record
               -> run_ppcont
                     -> create_directory (in common.c)
                     -> gen_pp_plot
                           -> move_file
                                 -> check_file (in common.c)
                                 -> call_system (in common.c)
                     -> gen_iface_plot
                           -> move_file (as above)
                     -> write_iface_txt
                           -> open_output_file (in common.c)
               -> r_get_dfn_data (in readpdb.c)
               -> read_chain_equivs
                     -> string_truncate (in common.c)
               -> count_chains
               -> run_promotif
                     -> create_directory (in common.c)
               -> run_psplot
                     -> get_motifs
                           -> string_chomp (in common.c)
                     -> create_directory (in common.c)
                     -> open_output_file (in common.c)
                     -> check_file (in common.c)
                     -> move_file (as above)
               -> run_resdata
                     -> create_directory (in common.c)
               -> run_hera
                     -> create_directory (in common.c)
                     -> open_output_file (in common.c)
                     -> genarate_topology
                           -> check_file (in common.c)
                           -> move_file (as above)
               -> run_speedfill
                     -> create_directory (in common.c)
                     -> gen_clefts_image
                           -> open_output_file (in common.c)
                           -> string_chomp (in common.c)
                           -> run_pymol
                                 -> call_system (in common.c)
               -> read_ligmet_run
                     -> create_ligmet_entry
                     -> create_resid_entry
               -> run_ligplots
                     -> create_directory (in common.c)
                     -> write_ligint_txt
                           -> open_output_file (in common.c)
                           -> remove_char (in common.c)
                           -> string_chomp (in common.c)
                     -> call_system (in common.c)
                     -> move_file (as above)
               -> run_nucplot
                     -> create_directory (in common.c)
               -> run_seqdata
                     -> create_directory (in common.c)
               -> run_uniplot
               -> run_wiring
                     -> create_directory (in common.c)
               -> run_pdbhtml
                     -> create_directory (in common.c)
               -> rasmol_scripts
                     -> create_directory (in common.c)
                     -> promotif_rasmol
                     -> main_rasmol
               -> pymol_scripts
                     -> create_directory (in common.c)
                     -> main_pymol
              -> run_showmain
                     -> create_directory (in common.c)
                     -> gen_protein_html
                     -> gen_procheck_html
                     -> gen_dna_html
                     -> gen_ligands_html
                           -> rewrite_js_file
                                 -> open_output_file (in common.c)
                                 -> string_chomp (in common.c)
                                 -> insert_file
                                       -> string_chomp (in common.c)
                     -> gen_interfaces_html
                     -> gen_wirpfam_html
                     -> gen_clefts_html
                     -> gen_index_html
                           -> rewrite_js_file (as above)
                     -> copy_images
   -> extract_pdb_code (as above)
   -> process_file (as above)
   -> write_index_html
   -> rewrite_runs_dbt
         -> string_chomp (in common.c)
         -> open_output_file (in common.c)
         -> remove_excess_blanks
         -> append_string (in common.c)
         -> replace_char (in common.c)

*/

#include "common.h"
#include "pdbsum1.h"
#include "readpdb.h"

/* Prototypes
   ---------- */

/* In common.c */
void append_string(char **string_ptr,char *string);
void call_system(char *command,int dos);
void c_get_paths(struct c_file_data *file_name_ptr,int run_64);
int check_file(char *file_name);
int create_directory(char *directory);
void free_tokens(char **token,int ntokens);
void get_date_and_time(char *date);
int get_tokens(char *line,char **token,int max_token,char token_sep);
FILE *open_output_file(char *out_name,int terminate);
void remove_char(char *string,char ch);
void replace_char(char *string,char replace,char by);
void store_string(char **string_ptr,char *string);
int string_chomp(char *string);
int string_truncate(char *string,int max_length);

/* In reapar.c */
int get_param(char *token_name,char *file_name,int filename_len,
              int exit_on_error);

/* In readpdb.c */
int r_get_dfn_data(char *pdb_code,int pdb_type,char *pdbsum_template,
                   struct r_rasdef **fst_r_rasdef_ptr,
                   struct r_ligand **fst_r_ligand_ptr,int *nlig,
                   int *nmet);

/***********************************************************************

get_command_arguments  -  Interpret the command line arguments

***********************************************************************/

void get_command_arguments(char *strptrs[],int num,
                           struct parameters *params)
{
  char file_name[FILENAME_LEN];

  char *string_ptr;

  int itoken, len;

/* v.1.0.1 --> */
  /* Save the program version */
  store_string(&(params->version),VERSION);
/* <-- v.1.0.1 */

  /* Initialise parameters */
  params->single_file = FALSE;
  params->max_cleft_chains = MAX_CLEFT_CHAINS;
  params->file = &null;
  params->pdb_code[0] = '\0';
  params->run_romlas = TRUE;
  params->relinks = TRUE;

  /* If nothing entered on the command line, show options available */
  if (num < 1 || (num == 1 && (!strcmp(strptrs[1],"-help") ||
                               !strcmp(strptrs[1],"-h"))))
    {
      printf("\n");
      printf("Usage is:\n");
      printf("\n");
      printf("pdbsum1 {file [pdb_code] | -l dataset.lst} [-nomaxchn]  "
             "[-noras]  [-help]\n");
      printf("\n");
      printf("where\n");
      printf("     * file           = Name of single PDB file to be "
             "processed\n");
      printf("     * pdb_code       = Optional 4-character identifier\n");
      printf("     * -l dataset.lst = Name of file listing PDB files to be "
             "processed\n");
      printf("     * -nomaxchn      = Clefts to be generated even if more "
             "than %d chains in structure\n",MAX_CLEFT_CHAINS);
      printf("     * -noras         = Don't run romlas to generate RasMol "
             "scripts\n");
      printf("\n");
      printf("\n");
      exit(1);
    }

  /* Loop through any additional command arguments */
  for (itoken = 1; itoken < num + 1; itoken++)
    {
      /* Check for the -l option */
      if (!strcmp(strptrs[itoken],"-l"))
        {
          /* Get the name of the dataet file from the next parameter */
          if (itoken < num)
            {
              /* Save the file name */
              store_string(&(params->file),strptrs[itoken + 1]);
            }
        }

      /* Check for the -nomaxchn option */
      else if (!strcmp(strptrs[itoken],"-nomaxchn"))
        params->max_cleft_chains = -1;

      /* Check for the -noras option, signifying that romlas doesn't
         need to be run to generate RasMol scripts */
      else if (!strcmp(strptrs[itoken],"-noras"))
        params->run_romlas = FALSE;

      /* Check for the -relinks option, signifying image links are
         to be relative rather than absolute, and that stock images
         are to be copied into gif directory */
      else if (!strcmp(strptrs[itoken],"-relinks"))
        params->relinks = TRUE;

      /* Otherwise, if this is the first parameter, take it to be the
         name of the PDB file to be processed */
      else if (itoken == 1)
        {
          /* Save the file name */
          store_string(&(params->file),strptrs[itoken]);

          /* Set flag that a single file to be processed */
          params->single_file = TRUE;
        }

      /* If second parameter, take it to be the assigned PDB code */
      else if (itoken == 2 && strlen(strptrs[itoken]) == 4)
        {
          /* Save the PDB code */
          strcpy(params->pdb_code,strptrs[itoken]);
        }
    }
}
/***********************************************************************

initialise_run  -  Initialise for current run, creating run.dat file

***********************************************************************/

FILE *initialise_run(int *run_no,char *date,
                     struct c_file_data file_name_ptr,
                     struct parameters params)
{
#define NIMAGE_FILES              3

  static char *IMAGE_FILE[NIMAGE_FILES]
    = { "spacer.gif", "PDBsum1_logo.png", "1x1.gif" };

  char command[LINELEN + 1], input_line[LINELEN + 1];
  char pdbsum1_dir[FILENAME_LEN], runs_dir[FILENAME_LEN];
  char file_name[FILENAME_LEN], gif_dir[FILENAME_LEN];
  char directory[FILENAME_LEN], image_file[FILENAME_LEN];

  int dos, err_no, image, loop, run_number;
  int have_file, ok;

  FILE *file_ptr, *runfile_ptr;

  /* Initialise variables */
  *run_no = run_number = 0;
  dos = file_name_ptr.dos;

  /* Initialise directory names */
  strcpy(pdbsum1_dir,file_name_ptr.other_dir[PDBSUM1_RESULTS_DIR]);

  /* Get the current date and time */
  get_date_and_time(date);

  /* Create the run data directory */
  sprintf(runs_dir,"%s/runs",pdbsum1_dir);

  /* Create the directory if necessary */
  create_directory(runs_dir);
  
  /* If using relative image links, create the gif directory for the
     common images */
  if (params.relinks == TRUE)
    {
      /* Loop over the two sets of image files required */
      for (loop = 0; loop < 2; loop++)
        {
          /* Form the name of the images directory */
          if (loop == 0)
            sprintf(gif_dir,"%s/gif",pdbsum1_dir);
          else
            sprintf(gif_dir,"%s/runs/gif",pdbsum1_dir);

          /* Create the directory if necessary */
          create_directory(gif_dir);

          /* Loop over the image files to be copied */
          for (image = 0; image < NIMAGE_FILES; image++)
            {
              /* Form the name of the image file */
              sprintf(image_file,"%s/%s",
                      file_name_ptr.other_dir[PDBSUM_TEMPLATES_GIF_DIR],
                      IMAGE_FILE[image]);

              /* Form destination name */
              sprintf(file_name,"%s/%s",gif_dir,IMAGE_FILE[image]);

              /* See if this file exists */
              have_file = check_file(file_name);

              /* Copy the required files across */
              if (have_file == FALSE)
                {
                  sprintf(command,"%s %s %s",CP_COMMAND,image_file,
                          gif_dir);
                  call_system(command,dos);
                }
            }
        }
      
      /* Create the css directory */
      sprintf(directory,"%s/css",pdbsum1_dir);
      create_directory(directory);

      /* Copy across all the css files */
      sprintf(command,"%s %s/css/* %s",CP_COMMAND,
              file_name_ptr.other_dir[PDBSUM_REAL_TEMPLATES_DIR],
              directory);
      call_system(command,dos);

      /* Repeat for the run subdirectory */
      sprintf(directory,"%s/runs/css",pdbsum1_dir);
      create_directory(directory);

      /* Copy across all the css files */
      sprintf(command,"%s %s/css/* %s",CP_COMMAND,
              file_name_ptr.other_dir[PDBSUM_REAL_TEMPLATES_DIR],
              directory);
      call_system(command,dos);
    }
  
  /* Form the name of the run counter, runno.txt */
  sprintf(file_name,"%s/runno.txt",runs_dir);

  /* Open the dataset.lst file */
  if ((file_ptr = fopen(file_name,"r")) !=  NULL)
    {
      /* Read in the first line containing the number */
      if (fgets(input_line,LINELEN,file_ptr) != NULL)
        run_number = atoi(input_line);

      /* Close the file */
      fclose(file_ptr);
    }

  /* Increment the run number */
  run_number++;

  /* Write the number out */
  file_ptr = open_output_file(file_name,TRUE);
  fprintf(file_ptr,"%d\n",run_number);

  /* Close the file */
  fclose(file_ptr);

  /* Form the name of the run.dat file */
  sprintf(file_name,"%s/run.dat",runs_dir);

  /* Open the file */
  runfile_ptr = open_output_file(file_name,TRUE);

  /* Write current run number and date of run */
  fprintf(runfile_ptr,"RUN_NUMBER %d\n",run_number);
  fprintf(runfile_ptr,"RUN_TIME %s\n",date);

  /* Return the current run number */
  *run_no = run_number;

  /* Return file pointer */
  return(runfile_ptr);
}
/***********************************************************************

create_pdb_entry  -  Create a pdb entry

***********************************************************************/

struct pdb *create_pdb_entry(struct pdb **fst_pdb_ptr,
                             struct pdb **lst_pdb_ptr,char *pdb_code)
{
  int isort;

  struct pdb *pdb_ptr, *first_pdb_ptr, *last_pdb_ptr;

  /* Initialise pdb pointers */
  pdb_ptr = NULL;
  first_pdb_ptr = *fst_pdb_ptr;
  last_pdb_ptr = *lst_pdb_ptr;

  /* Create pdb entry */
  pdb_ptr = (struct pdb *) malloc(sizeof(struct pdb));
  if (pdb_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for struct pdb\n");
      printf("***        Program analrefs terminated with error.\n");
      exit (1);
    }

  /* Store the PDB code */
  strcpy(pdb_ptr->pdb_code,pdb_code);

  /* Initialise remaining fields */
  pdb_ptr->file_name = &null;
  pdb_ptr->title = &null;

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

gen_pdb_code  -  Generate a PDB code

***********************************************************************/

int gen_pdb_code(char *pdb_code,struct c_file_data file_name_ptr)
{
  char command[LINELEN + 1], input_line[LINELEN + 1];
  char file_name[FILENAME_LEN], runs_dir[FILENAME_LEN];
  char parent_dir[FILENAME_LEN], pdir[FILENAME_LEN], sub_dir[3];
  char code_file[FILENAME_LEN];
  char ch, number_string[20];

  char *string_ptr;

  int i, ilet, inum, j, len;
  int done, have_file, ok, valid;

  FILE *file_ptr;

  /* Initialise variables */
  strcpy(pdb_code,"a000");
  ok = TRUE;
  valid = FALSE;

  /* Create the run data directory */
  sprintf(runs_dir,"%s/runs",
          file_name_ptr.other_dir[PDBSUM1_RESULTS_DIR]);
  create_directory(runs_dir);

  /* Form the name of the code.txt file holding the last-used code */
  sprintf(code_file,"%s/code.txt",runs_dir);

  /* Open the code.txt file */
  if ((file_ptr = fopen(code_file,"r")) !=  NULL)
    {
      /* Read in the first line containing the number */
      if (fgets(input_line,LINELEN,file_ptr) != NULL)
        {
          /* Truncate the input line */
          len = string_chomp(input_line);
          
          /* If looks OK, then save */
          if (len == 4)
            strcpy(pdb_code,input_line);
        }

      /* Close the file */
      fclose(file_ptr);
    }

  /* Loop until a valid code has been generated */
  while (valid == FALSE)
    {
      /* Initialise flag */
      done = FALSE;
  
      /* Check for end-codes */
      if (!strcmp(pdb_code,"z999"))
        strcpy(pdb_code,"aa00");
      else if (!strcmp(pdb_code,"zz99"))
        strcpy(pdb_code,"aaa0");
      else if (!strcmp(pdb_code,"zzz9"))
        strcpy(pdb_code,"aaaa");
      else if (!strcmp(pdb_code,"zzzz"))
        {
          /* End of the line reached - no more codes to generate */
          printf("*** ERROR. Cannot generate next PDB code from %s\n",
                 pdb_code);
          ok = FALSE;
          return(ok);;
        }

      /* Generate the next code */
      for (i = 3; i > -1 && done == FALSE; i--)
        {
          /* Get this character */
          ch = pdb_code[i];
          ilet = inum = 0;
      
          /* If this is a number, increment it */
          string_ptr = strchr(NUMBERS,ch);
          if (string_ptr != NULL)
            {
              /* Get the number */
              inum = string_ptr - NUMBERS;

              /* Increment the position */
              inum++;

              /* If not off end, then we're done */
              if (inum < strlen(NUMBERS))
                {
                  /* Save this number */
                  pdb_code[i] = NUMBERS[inum];
                  
                  /* Set flag that done */
                  done = TRUE;
                }
            }

          /* Do the same if position is a letter */
          else
            {
              /* Find letter */
              string_ptr = strchr(LETTERS,ch);

              /* Get the corresponding letter */
              ilet = string_ptr - LETTERS;

              /* Increment the position */
              ilet++;
              
              /* If not off end, then we're done */
              if (ilet < strlen(LETTERS))
                {
                  /* Save this number */
                  pdb_code[i] = LETTERS[ilet];
                  
                  /* Set flag that done */
                  done = TRUE;
                }
            }

          /* If not done, need to reset current position and all to the
             right */
          if (done == FALSE)
            {
              /* Loop for all characters to the right */
              for (j = i; j < 4; j++)
                {
                  /* Get this character */
                  ch = pdb_code[j];

                  /* If this is a number, reset to first */
                  string_ptr = strchr(NUMBERS,ch);
                  if (string_ptr != NULL)
                    pdb_code[j] = NUMBERS[0];

                  /* Otherwise, reset to first letter */
                  else
                    pdb_code[j] = LETTERS[0];
                }
            }
        }

      /* Form the name of the PDBsum directory and check that this
         code hasn't already been used */
      len = strlen(pdb_code);
      strncpy(sub_dir,pdb_code + 1,2);
      sub_dir[2] = '\0';

      /* Form the name of the parent directory */
      sprintf(parent_dir,"%s/%s",file_name_ptr.other_dir[PDBSUM1_RESULTS_DIR],
              sub_dir);

      /* Form the name of the PDBsum directory */
      sprintf(pdir,"%s/%s",parent_dir,pdb_code);

      /* Form the name of the corresponding database.dat file */
      sprintf(file_name,"%s/database.dat",pdir);
      
      /* See if this file exists */
      have_file = check_file(file_name);
      if (have_file == FALSE)
        valid = TRUE;
    }

  /* Write the new code out */
  file_ptr = open_output_file(code_file,TRUE);
  fprintf(file_ptr,"%s\n",pdb_code);

  /* Close the file */
  fclose(file_ptr);

  /* Return that OK */
  return(ok);
}
/***********************************************************************

extract_pdb_code  -  Extract the PDB code from the file name

***********************************************************************/

int extract_pdb_code(char *file_name,char *dir_name,char *pdb_code,
                     struct c_file_data file_name_ptr)
{
  char ch, tmp[FILENAME_LEN];

  char *string_ptr;

  int i, name_start, len;
  int done, have_file, ok;

  /* Initialise variables */
  have_file = ok = FALSE;
  strcpy(tmp,file_name);

  /* Check whether this file exists */
  have_file = check_file(file_name);
  if (have_file == FALSE)
    {
      /* Show error message */
      printf("*** ERROR. File not found: %s\n",file_name);
      
      /* Return flag */
      return(ok);
    }

  /* Get the length of the file name */
  len = strlen(tmp);
  ok = TRUE;

  /* Split name into directory name and file name */
  dir_name[0] = '\0';
  done = FALSE;
  name_start = 0;
  
  /* Find where file name starts */
  for (i = 0; i < len & done == FALSE; i++)
    {
      /* Get this character */
      ch = tmp[i];

      /* If have the dot, then end here */
      if (ch == '.')
        done = TRUE;

      /* If a directory separator, then reset name start */
      else if (ch == '/' || ch == '\\') 
        name_start = i + 1;
    }

  /* If directory separator located, then split the name */
  if (name_start > 0)
    {
      /* Separate the file name into directory name and file name */
      strcpy(dir_name,tmp);
      strcpy(tmp,dir_name + name_start);
      dir_name[name_start] = '\0';
    }

  /* Check the length of the filename is reasonable */
  if (len > 8 && len < FILENAME_LEN && !strncmp(tmp,"pdb",3))
    {
      /* Extract the PDB code from the name */
      string_ptr = strchr(tmp,'.');

      /* If not too long, then save */
      if (string_ptr != NULL)
        {
          /* Truncate name at dot */
          string_ptr[0] = '\0';

          /* Get the length of the code */
          len = strlen(tmp + 3);

          /* If code OK, then save the PDB code and create a
             new record */
          if (len > 0 && len <= MAX_PDB_CODE_LEN)
            {
              /* Save the code */
              strcpy(pdb_code,tmp + 3);
            }

          /* Cannot locate a valid PDB code */
          else
            ok = FALSE;
        }

      /* File name has no extension */
      else
        ok = FALSE;
    }

  /* Filename not of correct format to contain a PDB code */
  else
    ok = FALSE;

  /* If we don't have a PDB code from the file name, generate one */
  if (ok == FALSE)
    ok = gen_pdb_code(pdb_code,file_name_ptr);

  /* Return the flag */
  return(ok);
}
/***********************************************************************

move_file  -  Rename the given file, deleting the destination file if it
              already exists

***********************************************************************/

void move_file(char *directory,char *old_name,char *new_name,int dos)
{
  char command[LINELEN + 1], file_name[FILENAME_LEN];
  
  int have_file;
  
  /* Form the full name of the destination file */
  sprintf(file_name,"%s/%s",directory,new_name);

  /* See if the destination file exists */
  have_file = check_file(file_name);

  /* If file exists, then delete it */
  if (have_file == TRUE)
    {
      /* Form the delete command and delete file */
      sprintf(command,"%s %s",RM_COMMAND,file_name);
      call_system(command,dos);
    }

  /* Rename the file */
  sprintf(command,"cd %s && %s %s %s",directory,MV_COMMAND,old_name,
          new_name);
  call_system(command,dos);
}
/***********************************************************************

initialise  -  Create PDBsum directory and transfer file

***********************************************************************/

int initialise(char *file_name,char *dir_name,char *pdb_code,char *pdir,
               char *fname,char *orig_fname,
               struct c_file_data file_name_ptr)
{
  char command[LINELEN + 1], parent_dir[FILENAME_LEN], sub_dir[3];
  char ch, new_name[FILENAME_LEN], term[5];

  char *string_ptr;

  int i, ipos, len;
  int dos, gzipped, ok, zipped;

  /* Initialise variables */
  dos = file_name_ptr.dos;
  ok = TRUE;
  gzipped = zipped = FALSE;
  term[0] = '\0';

  /* Get the 2nd and 3rd digits of the PDB code */
  len = strlen(pdb_code);
  strncpy(sub_dir,pdb_code + 1,2);
  sub_dir[2] = '\0';

  /* Form the name of the parent directory */
  sprintf(parent_dir,"%s/%s",
          file_name_ptr.other_dir[PDBSUM1_RESULTS_DIR],sub_dir);

  /* Create this directory if it doesn't exist */
  create_directory(parent_dir);

  /* Form the name of the PDBsum directory */
  sprintf(pdir,"%s/%s",parent_dir,pdb_code);

  /* Create this directory if it doesn't exist */
  create_directory(pdir);

  /* Save the original file name */
  strcpy(orig_fname,file_name);

  /* Remove any path */
  len = strlen(orig_fname);
  ipos = 0;
  for (i = 0; i < len; i++)
    {
      /* If directory separator, save position */
      ch = orig_fname[i];
      if (ch == '\\' || ch == '/')
        ipos = i + 1;
    }

  /* If have directory separators, remove path */
  if (ipos > 0)
    {
      /* Transfer just the file name */
      strcpy(new_name,orig_fname + ipos);
      strcpy(orig_fname,new_name);
      new_name[0] = '\0';
    }

  /* See if file is zipped or gzipped */
  string_ptr = strstr(file_name,".gz");
  if (string_ptr != NULL)
    {
      gzipped = TRUE;
      strcpy(term,".gz");
    }
  string_ptr = strstr(file_name,".zip");
  if (string_ptr != NULL)
    {
      zipped = TRUE;
      strcpy(term,".zip");
    }

  /* Form the new file name */
  sprintf(new_name,"%s/%s.pdb",pdir,pdb_code);
  strcat(new_name,term);

  /* Copy file to the PDBsum directory */
  sprintf(command,"%s %s %s",CP_COMMAND,file_name,new_name);
  call_system(command,dos);

  /* Unzip file if required */
  if (gzipped == TRUE)
    {
      sprintf(command,"gunzip -f %s",new_name);
      call_system(command,dos);
    }
  else if (zipped == TRUE)
    {
      sprintf(command,"unzip -f %s",new_name);
      call_system(command,dos);
    }

  /* Form new name without the directory name or extension */
  sprintf(fname,"%s.pdb",pdb_code);

  /* Return if all OK */
  return(ok);
}
/***********************************************************************

run_fixpdb  -  Run fixpdb on the input PDB file

***********************************************************************/

void run_fixpdb(char *fno,int process,char *pdir,char *pdb_code,
                char *file_name,char *file_new,char *log_dir,
                struct c_file_data file_name_ptr,int dos)
{
  char command[LINELEN + 1], directory[FILENAME_LEN];
  char log_file[FILENAME_LEN];

  int len;
  int ok;

  /* Initialise variables */
  ok = TRUE;

  /* Show progress */
  printf("   %s%d. In run_fixpdb ...\n",fno,process);

  /* Save the original file */
  sprintf(command,"cd %s && %s %s original.pdb",pdir,CP_COMMAND,
          file_name);
  call_system(command,dos);

  /* Form name of the log file */
  sprintf(log_file,"%s/fixpdb.log",log_dir);

  /* Show progress */
  printf("          - Running fixpdb ...\n");

  /* Run the program */
  sprintf(command,"cd %s && %s/fixpdb %s > %s",pdir,
          file_name_ptr.other_dir[PDBSUM_EXE_DIR],file_name,
          log_file);
  if (dos == TRUE)
    replace_char(command,'/','\\');
  call_system(command,dos);
}
/***********************************************************************

write_in_file  -  Write the given commands to the in file

***********************************************************************/

void write_in_file(char *directory,char *commands)
{
  char file_name[FILENAME_LEN];

  FILE *file_ptr;

  /* Form the output file name */
  sprintf(file_name,"%s/in",directory);

  /* Open the file */
  file_ptr = open_output_file(file_name,TRUE);

  /* Write out the commands */
  fprintf(file_ptr,"%s\n",commands);

  /* Close the file */
  fclose(file_ptr);
}
/***********************************************************************

run_clean  -  Run CLEAN on the input PDB file

***********************************************************************/

void run_clean(char *fno,int process,char *pdir,char *pdb_code,
               char *file_name,char *file_new,char *log_dir,
               struct c_file_data file_name_ptr)
{
  char command[LINELEN + 1], directory[FILENAME_LEN];
  char log_file[FILENAME_LEN], param_dir[FILENAME_LEN];
  char fix_name[FILENAME_LEN], in_file[FILENAME_LEN];
  char from_name[FILENAME_LEN], copy_name[FILENAME_LEN];

  char *string_ptr;

  int len;
  int ok, dos;

  /* Initialise variables */
  dos = file_name_ptr.dos;
  ok = TRUE;
  strcpy(param_dir,file_name_ptr.other_dir[PDBSUM_PARAMS_DIR]);

  /* Show progress */
  printf("   %s%d. In run_clean ...\n",fno,process);

  /* Form name of the log file */
  sprintf(log_file,"%s/clean.log",log_dir);

  /* Create a directory called procheck */
  sprintf(directory,"%s/procheck",pdir);
  create_directory(directory);

  /* Copy across the required parameter files */
  sprintf(command,"%s %s/procheck.dat %s/prodata",CP_COMMAND,
          param_dir,directory);
  call_system(command,dos);
  sprintf(command,"%s %s/resdefs.dat %s",CP_COMMAND,param_dir,
          directory);
  call_system(command,dos);
  sprintf(command,"%s %s/runprocheck.prm %s/procheck.prm",CP_COMMAND,
          param_dir,directory);
  call_system(command,dos);

  /* Define input parameter file */
  sprintf(in_file,"%s/in",directory);

  /* Get the name of the fixed PDB file */
  strcpy(fix_name,file_name);
  string_ptr = strstr(fix_name,".pdb");
  if (string_ptr != NULL)
    {
      string_ptr[0] = '\0';
      strcat(fix_name,".fix");
    }

  /* Copy fixed file */
  if (strcmp(fix_name,file_name))
    {
      sprintf(command,"cd %s && %s %s %s",pdir,CP_COMMAND,fix_name,
              file_name);
      call_system(command,dos);
    }

  /* Write parameters for clean */
  sprintf(command,"../%s\n",fix_name);
  write_in_file(directory,command);

  /* Show progress */
  printf("          - Running clean ...\n");

  /* Run the program */
  sprintf(command,"cd %s && %s/clean < %s > %s",directory,
          file_name_ptr.other_dir[PDBSUM_EXE_DIR],in_file,log_file);
  call_system(command,dos);

  /* Form name of the .new file created */
  sprintf(file_new,"%s/%s.new",directory,pdb_code);

  /* Form name of the secstr log file */
  sprintf(log_file,"%s/secstr.log",log_dir);

  /* Write parameters for secstr */
  sprintf(command,"+%s\n",fix_name);
  write_in_file(directory,command);
  //  sprintf(command,"echo +%s > %s",fix_name,in_file);
  //  call_system(command,dos);

  /* Show progress */
  printf("          - Running secstr ...\n");

  /* Run the program */
  sprintf(command,"cd %s && %s/secstr < %s > %s",directory,
          file_name_ptr.other_dir[PDBSUM_EXE_DIR],in_file,log_file);
  call_system(command,dos);

  /* Copy the .new file to the directory above */
  sprintf(command,"%s %s %s",CP_COMMAND,file_new,pdir);
  call_system(command,dos);

  /* Rename the .new file */
  sprintf(file_new,"%s.new",pdb_code);
}
/***********************************************************************

list_plots  -  Create a database.dat file of the PostScript files
               produced and convert to PDF

***********************************************************************/

void list_plots(char *directory,char *pdb_code,char *ps2pdf,
                struct c_file_data file_name_ptr)
{
  char file_name[FILENAME_LEN], input_line[LINELEN + 1];
  char out_file[FILENAME_LEN], pdf_file[FILENAME_LEN];
  char command[FILENAME_LEN];

  char *string_ptr;
  
  int len, nfile, npdf;
  int done, finish, found, have_file, dos;

  FILE *file_ptr, *outfile_ptr;

  /* Initialise variables */
  done = finish = FALSE;
  nfile = 0;
  npdf = 0;
  dos = file_name_ptr.dos;

  /* Show progress */
  printf("          - Converting files ...\n");

  /* Form the name of the output database.dat file */
  strcpy(out_file,directory);
  strcat(out_file,"/database.dat");

  /* Open the output file */
  outfile_ptr = open_output_file(out_file,TRUE);

  /* Loop while going through the PostScript files */
  while (done == FALSE)
    {
      /* Initialise finish flag */
      finish = FALSE;

      /* Increment file counter */
      nfile++;
      found = FALSE;

      /* Form the name of the PostScript file */
      sprintf(file_name,"%s/%s_%2.2d.ps",directory,pdb_code,nfile);
  
      /* Open this file */
      if ((file_ptr = fopen(file_name,"r")) !=  NULL)
        {
          /* Read through the file to get the title */
          while (fgets(input_line,LINELEN,file_ptr) != NULL &&
                 done == FALSE && finish == FALSE)
            {
              /* Check line-length */
              len = string_chomp(input_line);

              /* If this is the title line, extract it */
              if (!strncmp(input_line,"%%Title: ",9))
                {
                  /* Remove everything from the comma onwards */
                  string_ptr = strchr(input_line,',');
                  if (string_ptr != NULL)
                    string_ptr[0] = '\0';

                  /* Write to the database.dat file */
                  fprintf(outfile_ptr,"PPLOT[%d] %s_%2.2d.ps\n",nfile,
                          pdb_code,nfile);
                  fprintf(outfile_ptr,"PPLOT_NAME[%d] %s\n",nfile,
                          input_line + 9);

                  /* Set flag that we're done with this file */
                  finish = TRUE;
                  found = TRUE;
                }
            }

          /* Close the file */
          fclose(file_ptr);

          /* If OK, convert Postcript file to PDF */
          if (found == TRUE && strcmp(ps2pdf,"NONE"))
            {
              /* Form the name of the PDF file */
              sprintf(pdf_file,"%s/%s_%2.2d.pdf",directory,pdb_code,
                      nfile);
  
              /* Convert the PostScript file */
              sprintf(command,"cd %s && %s %s %s",directory,ps2pdf,
                      file_name,pdf_file);
              call_system(command,dos);

              /* Increment count of PDF files */
              npdf++;
            }
        }

      else
        done = TRUE;
    }

  /* Write out the number of plots there are */
  fprintf(outfile_ptr,"NPPLOT %d\n",nfile);

  /* Add flag to database.dat file that have PDF file */
  if (npdf > 0)
    fprintf(outfile_ptr,"HAVE_PDF TRUE\n");

  /* Check whether the residue listing file exists */
  sprintf(out_file,"%s/%s.out",directory,pdb_code);
  have_file = check_file(out_file);
  if (have_file == TRUE)
    fprintf(outfile_ptr,"OUT_FILE TRUE\n");

  /* Repeat for the gfact.rasmol file */
  sprintf(out_file,"%s/%s_gfact.rasmol",directory,pdb_code);
  have_file = check_file(out_file);
  if (have_file == TRUE)
    fprintf(outfile_ptr,"GFAC2PDB TRUE\n");

  /* Finish off */
  fprintf(outfile_ptr,"PLOTS_DONE TRUE\n");

  /* Close output file */
  fclose(outfile_ptr);
}
/***********************************************************************

run_procheck  -  Run the remainder of the PROCHECK programs

***********************************************************************/

void run_procheck(char *fno,int process,char *pdir,char *pdb_code,
                  char *file_name,char *file_new,float resol,
                  char *log_dir,struct c_file_data file_name_ptr)
{
  char command[LINELEN + 1], directory[FILENAME_LEN];
  char log_file[FILENAME_LEN], prog_params[FILENAME_LEN];
  char in_file[FILENAME_LEN], out_file[FILENAME_LEN];

  int len;
  int ok, dos;

  /* Initialise variables */
  ok = TRUE;
  dos = file_name_ptr.dos;

  /* Show progress */
  printf("   %s%d. In run_procheck ...\n",fno,process);

  /* Define directory */
  sprintf(directory,"%s/procheck",pdir);

  /* == RASCRIPT == */

  /* Show progress */
  //printf("          - Running rascript ...\n");

  /* Form name of the rascript.rasmol output file */
  //sprintf(out_file,"%s/rascript.rasmol",directory);

  /* Define the input parameters for rascript */
  //sprintf(params,"%s -dir %s -arad -1.0",pdb_code,directory,out_file);

  /* Run the program */
  //sprintf(command,"cd %s && %s/rascript %s > %s",directory,
  //        file_name_ptr.other_dir[PDBSUM_EXE_DIR],prog_params,out_file);
  //call_system(command,dos);

  /* == NB == */

  /* Show progress */
  printf("          - Running nb ...\n");

  /* Define input parameter file */
  sprintf(in_file,"%s/in",directory);

  /* Form name of the log file */
  sprintf(log_file,"%s/nb.log",log_dir);

  /* Write out the input parameters */
  sprintf(command,"%s\n",file_name);
  write_in_file(directory,command);

  //  sprintf(command,"echo %s > %s",file_name,in_file);
  //  call_system(command,dos);

  /* Run the program */
  sprintf(command,"cd %s && %s/nb < %s > %s",directory,
          file_name_ptr.other_dir[PDBSUM_EXE_DIR],in_file,log_file);
  call_system(command,dos);

  /* == ANGLEN == */

  /* Show progress */
  printf("          - Running anglen ...\n");

  /* Form name of the log file */
  sprintf(log_file,"%s/anglen.log",log_dir);

  /* Run the program */
  sprintf(command,"cd %s && %s/anglen < %s > %s",directory,
          file_name_ptr.other_dir[PDBSUM_EXE_DIR],in_file,log_file);
  call_system(command,dos);

  /* == TPLOT == */

  /* Show progress */
  printf("          - Running tplot ...\n");

  /* Form name of the log file */
  sprintf(log_file,"%s/tplot.log",log_dir);

  /* Append resolution to input parameters */
  sprintf(command,"%s\n \n%8.3f\nN\n",file_name,resol);
  write_in_file(directory,command);

  //  sprintf(command,"echo  >> %s",in_file);
  //  call_system(command,dos);
  //  sprintf(command,"echo %8.3f >> %s",resol,in_file);
  //  call_system(command,dos);
  //  sprintf(command,"echo N >> %s",in_file);
  //  call_system(command,dos);

  /* Run the program */
  sprintf(command,"cd %s && %s/tplot < %s > %s",directory,
          file_name_ptr.other_dir[PDBSUM_EXE_DIR],in_file,log_file);
  call_system(command,dos);

  /* Define parameters to convert the Ramachandran plot to a small
     GIF image */
  sprintf(prog_params,
          "-antialias -crop 500x470+5+39 -transparent white +repage");

  /* Convert .ps file to GIF */
/* v.1.0.1 --> */
//  sprintf(command,"cd %s && %s %s ps:%s_01.ps gif:ramach.gif",directory,
//          file_name_ptr.exe_name[CONVERT_EXE],prog_params,pdb_code);
  sprintf(command,"cd %s && %s ps:%s_01.ps %s gif:ramach.gif",directory,
          file_name_ptr.exe_name[CONVERT_EXE],pdb_code,prog_params);
/* <-- v.1.0.1 */
  call_system(command,dos);

  /* == PPLOT == */

  /* Show progress */
  printf("          - Running pplot ...\n");

  /* Form name of the log file */
  sprintf(log_file,"%s/pplot.log",log_dir);

  /* Append resolution to input parameters */
  sprintf(command,"%s\n \n%8.3f\n",file_name,resol);
  write_in_file(directory,command);

  //  sprintf(command,"echo %s > %s",file_new,in_file);
  //  call_system(command,dos);
  //  sprintf(command,"echo  >> %s",resol,in_file);
  //  call_system(command,dos);
  //  sprintf(command,"echo %8.3f >> %s",resol,in_file);
  //  call_system(command,dos);

  /* Run the program */
  sprintf(command,"cd %s && %s/pplot < %s > %s",directory,
          file_name_ptr.other_dir[PDBSUM_EXE_DIR],in_file,log_file);
  call_system(command,dos);

  /* == BPLOT == */

  /* Show progress */
  printf("          - Running bplot ...\n");

  /* Form name of the log file */
  sprintf(log_file,"%s/bplot.log",log_dir);

  /* Run the program */
  sprintf(command,"cd %s && %s/bplot < %s > %s",directory,
          file_name_ptr.other_dir[PDBSUM_EXE_DIR],in_file,log_file);
  call_system(command,dos);

  /* == GFAC2PDB == */

  /* Show progress */
  printf("          - Running gfac2pdb ...\n");

  /* Form name of the log file */
  sprintf(log_file,"%s/gfac2pdb.log",log_dir);

  /* Run the program */
  sprintf(command,"cd %s && %s/gfac2pdb < %s > %s",directory,
          file_name_ptr.other_dir[PDBSUM_EXE_DIR],in_file,log_file);
  call_system(command,dos);

  /* List the plots in a database.dat file, and convert to PDF */
  list_plots(directory,pdb_code,file_name_ptr.exe_name[PS2PDF],
             file_name_ptr);
}
/***********************************************************************

run_hbplus  -  Run the remainder of the HBPLUS programs

***********************************************************************/

void run_hbplus(char *fno,int process,char *pdir,char *pdb_code,
                char *file_name,char *file_new,char *log_dir,
                struct c_file_data file_name_ptr,
                struct parameters params)
{
  char command[LINELEN + 1], directory[FILENAME_LEN];
  char log_file[FILENAME_LEN], prog_params[FILENAME_LEN];
  char add_param[FILENAME_LEN], out_file[FILENAME_LEN];

  int len;
  int have_file, ok, dos;

  float resol;

  /* Initialise variables */
  add_param[0] = '\0';
  ok = TRUE;
  dos = file_name_ptr.dos;

  /* Show progress */
  printf("   %s%d. In run_hbplus ...\n",fno,process);

  /* Create the HBPLUS directory */
  sprintf(directory,"%s/hbplus",pdir);
  create_directory(directory);

  /* == HBADD == */

  /* Show progress */
  printf("          - Running hbadd ...\n");

  /* Form name of the log file */
  sprintf(log_file,"%s/hbadd.log",log_dir);

  /* Define the command-line parameters */
  sprintf(prog_params,"%s/%s %s -pdbsum",pdir,file_new,
          file_name_ptr.other_dir[HETGROUP_DICTIONARY]);

  /* Run the program */
  sprintf(command,"cd %s && %s/hbadd %s > %s",directory,
          file_name_ptr.other_dir[PDBSUM_EXE_DIR],prog_params,log_file);
  call_system(command,dos);

  /* == HHB == */

  /* Show progress */
  printf("          - Running hhb ...\n");

  /* Form name of the log file */
  sprintf(log_file,"%s/hhb.log",log_dir);

  /* See if the hbplus.rc file exists */
  sprintf(out_file,"%s/hbplus.rc",directory);
  have_file = check_file(out_file);

  /* If have the file, add extra parameter to HBPLUS call */
  if (have_file == TRUE)
    sprintf(add_param,"-f %s",out_file);
 
  /* Define the command-line parameters */
  sprintf(prog_params,"-L %s -h 2.7 -d 3.35 %s/%s",add_param,pdir,file_new);

  /* Run the program */
  sprintf(command,"cd %s && %s/hbplus %s >> %s",directory,
          file_name_ptr.other_dir[PDBSUM_EXE_DIR],prog_params,log_file);
  call_system(command,dos);

  /* Repeat for the .hb2 output */
  sprintf(prog_params,"%s -h 2.7 -d 3.35 %s/%s",add_param,pdir,file_new);

  /* Run the program */
  sprintf(command,"cd %s && %s/hbplus %s >> %s",directory,
          file_name_ptr.other_dir[PDBSUM_EXE_DIR],prog_params,log_file);
  call_system(command,dos);

  /* == NNB == */

  /* Show progress */
  printf("          - Running nnb ...\n");

  /* Form name of the log file */
  sprintf(log_file,"%s/nnb.log",log_dir);

  /* Define the command-line parameters */
  sprintf(prog_params,"-L %s -h 2.9 -d 3.9 -N %s/%s",add_param,pdir,file_new);

  /* Run the program */
  sprintf(command,"cd %s && %s/hbplus %s > %s",directory,
          file_name_ptr.other_dir[PDBSUM_EXE_DIR],prog_params,log_file);
  call_system(command,dos);

  /* Repeat for the .nb2 output */
  sprintf(prog_params,"%s -h 2.9 -d 3.9 -N %s/%s",add_param,pdir,file_new);

  /* Run the program */
  sprintf(command,"cd %s && %s/hbplus %s > %s",directory,
          file_name_ptr.other_dir[PDBSUM_EXE_DIR],prog_params,log_file);
  call_system(command,dos);
}
/***********************************************************************

run_newsec  -  Run newsec to analyse contents of the PDB file

***********************************************************************/

void run_newsec(char *fno,int process,char *pdir,char *pdb_code,
                char *file_name,char *file_new,char *log_dir,
                char *param_dir,char *pdbsum_exe_dir,char *pymol_exe,
                struct c_file_data file_name_ptr)
{
  char command[LINELEN + 1], directory[FILENAME_LEN];
  char log_file[FILENAME_LEN], prog_params[FILENAME_LEN];
  char add_param[FILENAME_LEN], dos_exit[20], dos_start[20];

  int len;
  int have_file, ok, dos;

  float resol;

  /* Initialise variables */
  add_param[0] = '\0';
  dos = file_name_ptr.dos;
  dos_exit[0] = dos_start[0] = '\0';
  ok = TRUE;

  /* Show progress */
  printf("   %s%d. In run_newsec ...\n",fno,process);

  /* Create the NEWSEC directory */
  sprintf(directory,"%s/newsec",pdir);
  create_directory(directory);

  /* Copy files required by NACCESS into current directory */
  sprintf(command,"%s %s/vdw.radii %s",CP_COMMAND,param_dir,
          directory);
  call_system(command,dos);
  sprintf(command,"%s %s/standard.data %s",CP_COMMAND,param_dir,
          directory);
  call_system(command,dos);

  /* Form name of the log file */
  sprintf(log_file,"%s/newsec.log",log_dir);

  /* For DOS, need to use the start command to prevent PyMOL shutting
     current window */
  if (dos == TRUE)
    {
      strcpy(dos_start,"start \"\" ");
      strcpy(dos_exit," && exit");
    }

  /* If PyMOL executable not defined, set flag accordingly */
  if (!strcmp(pymol_exe,"NONE"))
    strcpy(add_param,"-nogifs");
  else
    strcpy(add_param,"-pymol");

  /* Define the command-line parameters */
  sprintf(prog_params,"%s %s/%s %s -pdbsum1",pdb_code,pdir,file_name,
          add_param);

  /* Run the program */
  sprintf(command,"cd %s && %s %s/newsec %s > %s %s",directory,dos_start,
          pdbsum_exe_dir,prog_params,log_file,dos_exit);
  call_system(command,dos);
}
/***********************************************************************

wait_for_file  -  In DOS need to wait until newsec has finished and the
                  rasmol.dfn file has been created

***********************************************************************/

void wait_for_file(char *pdir,int dos)
{
  char chains_dat[FILENAME_LEN], rasmol_dfn[FILENAME_LEN];
  char command[LINELEN + 1];

  int nfound, ntries;
  int have_file;

  /* Initialise variables */
  nfound = ntries = 0;

  /* Form the name of the rasmol.difn and chains.dat files */
  sprintf(rasmol_dfn,"%s/newsec/rasmol.dfn",pdir);
  replace_char(rasmol_dfn,'/','\\');
  sprintf(chains_dat,"%s/newsec/chains.dat",pdir);
  replace_char(chains_dat,'/','\\');

  /* Loop until both files have been created */
  while (nfound < 2 && ntries < 100)
    {
      /* Initialise count */
      nfound = 0;
      
      /* Check whether the rasmol.dfn file exists */
      have_file = check_file(rasmol_dfn);
      if (have_file == TRUE)
        nfound++;
      
      /* Repeat for the chains.dat file */
      have_file = check_file(chains_dat);
      if (have_file == TRUE)
        nfound++;

      /* Increment count of tries */
      ntries++;

      /* Pause before next try */
      sprintf(command,"%s",SLEEP_COMMAND);
      call_system(command,dos);
    }
}
/***********************************************************************

find_chain  -  Find the given chain

***********************************************************************/

struct chain *find_chain(struct chain *first_chain_ptr,char chain)
{
  struct chain *chain_ptr, *found_chain_ptr;

  /* Initialise variables */
  found_chain_ptr = NULL;

  /* Get pointer to the first chain */
  chain_ptr = first_chain_ptr;

  /* Loop over all the chains to find the one just read in */
  while (chain_ptr != NULL && found_chain_ptr == NULL)
    {
      /* If have a match, store the pointer */
      if (chain == chain_ptr->chain)
        found_chain_ptr = chain_ptr;

      /* Get pointer to the next chain */
      chain_ptr = chain_ptr->next_chain_ptr;
    }

  /* Return the matching chain pointer */
  return(found_chain_ptr);
}
/***********************************************************************

create_chain_record  -  Create and initialise a new chain record

***********************************************************************/

struct chain *create_chain_record(struct chain **fst_chain_ptr,
                                  struct chain **lst_chain_ptr,
                                  char chain,int number)
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
  chain_ptr->chain = chain;
  chain_ptr->number = number;
  chain_ptr->chain_no = number;
  chain_ptr->first_iface_ptr = NULL;
  chain_ptr->last_iface_ptr = NULL;
  chain_ptr->nifaces = 0;
  chain_ptr->next_chain_ptr = NULL;

  /* Return current pointers */
  *fst_chain_ptr = first_chain_ptr;
  *lst_chain_ptr = last_chain_ptr;

  /* Return new chain pointer */
  return(chain_ptr);
}
/***********************************************************************

create_iface_record  -  Create and initialise a new iface record

***********************************************************************/

struct iface *create_iface_record(struct iface **fst_iface_ptr,
                                          struct iface **lst_iface_ptr,
                                          struct chain *other_chain_ptr)
{
  int i;

  struct iface *iface_ptr, *first_iface_ptr, *last_iface_ptr;

  /* Initialise iface pointers */
  iface_ptr = NULL;
  first_iface_ptr = *fst_iface_ptr;
  last_iface_ptr = *lst_iface_ptr;

  /* Allocate memory for structure to hold iface info */
  iface_ptr = (struct iface *) malloc(sizeof(struct iface));
  if (iface_ptr == NULL)
    {
      printf("*** Can't allocate memory for struct iface\n");
      exit (1);
    }

  /* If this is the very first iface, then store pointer */
  if (first_iface_ptr == NULL)
    first_iface_ptr = iface_ptr;

  /* Add link from previous iface to the current one */
  if (last_iface_ptr != NULL)
    last_iface_ptr->next_iface_ptr = iface_ptr;
  last_iface_ptr = iface_ptr;

  /* Initialise interface details */
  iface_ptr->iface_no = 0;
  for (i = 0; i < 2; i++)
    {
      iface_ptr->nresid[i] = 0;
      iface_ptr->chain_area[i] = 0.0;
      iface_ptr->iface_area[i] = 0.0;
    }
  for (i = 0; i < NPP_TYPES; i++)
    iface_ptr->count[i] = 0;
  iface_ptr->ntypes = 0;
  iface_ptr->chain_ptr = other_chain_ptr;
  iface_ptr->next_iface_ptr = NULL;

  /* Return current pointers */
  *fst_iface_ptr = first_iface_ptr;
  *lst_iface_ptr = last_iface_ptr;

  /* Return new iface pointer */
  return(iface_ptr);
}
/***********************************************************************

get_interfaces  -  Read in the all the chains and their interface data
                   from the interfaces.dat file

***********************************************************************/

int get_interfaces(char *pdir,struct chain **fst_chain_ptr,
                   struct chain **lst_chain_ptr,int *nchn)
{
#define NKEY  8

  char input_line[LINELEN + 1], file_name[FILENAME_LEN];
  char tmp_line[LINELEN + 1];
  char value[LINELEN + 1];
  char chains[2], number_string[20];

  char *string_ptr;

  static char *keyword[] = {
    "PPIFACE_CHAIN", "PPIFACE_NRESID", "PPIFACE_AREA", "PPIFACE_BURIED",
    "PPIFACE_NHBONDS", "PPIFACE_NCONTACTS", "PPIFACE_SALT_BRIDGES",
    "PPIFACE_DISULPHIDES",
  };

  int ichain, iface_no, ikey, key, klen, len, nchains, niface;

  static int keylen[NKEY] = { 13, 14, 12, 14, 15, 17, 20, 19 };

  struct chain *chain_ptr[2], *first_chain_ptr, *last_chain_ptr;
  struct iface *iface_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  nchains = *nchn;
  niface = 0;

  /* Initialise chain pointers */
  chain_ptr[0] = NULL;
  chain_ptr[1] = NULL;
  first_chain_ptr = *fst_chain_ptr;
  last_chain_ptr = *lst_chain_ptr;
  
  /* Form name of the interfaces.dat file */
  sprintf(file_name,"%s/newsec/interfaces.dat",pdir);

  /* Open the file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    return(niface);

  /* Loop while reading in records from the file */
  while (fgets(input_line,LINELEN,file_ptr) != NULL)
    {
      /* Truncate the string to strip off the newline */
      len = string_truncate(input_line,LINELEN);

      /* Check for a "copy" interface */
      if (!strncmp(input_line,"COPY_",5))
        {
          /* Remove the "COPY_" part of the keyword */
          strcpy(tmp_line,input_line+5);
          strcpy(input_line,tmp_line);
          len = len -5;
        }

      /* Initisliae flag */
      key = -1;

      /* Check whether have one of the required key words */
      for (ikey = 0; ikey < NKEY && key == -1; ikey++)
        {
          /* Check whether this is the required keyword */
          if (!strncmp(input_line,keyword[ikey],keylen[ikey]))
            key = ikey;
        }

      /* If have one of the keywords applying to a specific chain, then
         extract which chain it is */
      if (key > -1 && key < 4)
        {
          /* Get length of current keyword */
          klen = strlen(keyword[key]);

          /* Get which chain this is */
          number_string[0] = input_line[klen];
          number_string[1] = '\0';
          ichain = atoi(number_string) - 1;
        }
      else
        ichain = -1;

      /* Extract the field value */
      string_ptr = strstr(input_line," ");
      if (string_ptr != NULL)
        strcpy(value,string_ptr+1);
      else
        value[0] = '\0';

      /* If this is a chain identifier, find/create chain record */
      if (key == 0)
        {
          /* Initialise interface pointer */
          iface_ptr = NULL;

          /* Store the chain */
          chains[ichain] = value[0];

          /* Check whether we already have this chain */
          chain_ptr[ichain] = find_chain(first_chain_ptr,chains[ichain]);

          /* If don't have chain, create chain record */
          if (chain_ptr[ichain] == NULL)
            {
              /* Create chain record */
              chain_ptr[ichain]
                = create_chain_record(&first_chain_ptr,&last_chain_ptr,
                                      chains[ichain],nchains);

             /* Increment count of new chains */
              nchains++;
            }

          /* If this is the second chain, create a link record representing
             the interface between the two chains */
          if (ichain == 1 && chain_ptr[0] != NULL && chain_ptr[1] != NULL)
            {
              /* Create the iface record for first chain */
              iface_ptr
                = create_iface_record(&(chain_ptr[0]->first_iface_ptr),
                                      &(chain_ptr[0]->last_iface_ptr),
                                      chain_ptr[1]);

              /* Increment each chain's count of ifaces */
              chain_ptr[0]->nifaces++;
              chain_ptr[1]->nifaces++;

              /* Increment count of interfaces */
              niface++;

              /* Save the interface number */
              iface_ptr->iface_no = niface;
            }
        }

      /* If have wanted keyword and valid interface record store details */
      else if (key > -1 && iface_ptr != NULL)
        {
          /* Number of interface residues */
          if (key == 1 && ichain > -1)
            iface_ptr->nresid[ichain] = atoi(value);

          /* Chain area */
          else if (key == 2 && ichain > -1 && chain_ptr[ichain] != NULL)
            {
              /* Store chain's area */
              iface_ptr->chain_area[ichain] = atof(value);
            }

          /* Interface area */
          else if (key == 3 && ichain > -1 && chain_ptr[ichain] != NULL)
            {
              /* Store interface area */
              iface_ptr->iface_area[ichain] = atof(value);
            }

          /* Number of hydrogen bonds */
          else if (key == 4)
            {
              /* Store number of hydrogen bonds */
              iface_ptr->count[PP_HBONDS] = atoi(value);
              if (iface_ptr->count[PP_HBONDS] > 0)
                iface_ptr->ntypes++;
            }

          /* Number of non-bonded contacts */
          else if (key == 5)
            {
              /* Store number of non-bonded contacts */
              iface_ptr->count[PP_CONTACTS] = atoi(value);
              if (iface_ptr->count[PP_CONTACTS] > 0)
                iface_ptr->ntypes++;
            }

          /* Number of salt bridges */
          else if (key == 6)
            {
              /* Store number of salt bridges */
              iface_ptr->count[PP_SALT_BRIDGES] = atoi(value);
              if (iface_ptr->count[PP_SALT_BRIDGES] > 0)
                iface_ptr->ntypes++;
            }

          /* Number of disulphide bonds */
          else if (key == 7)
            {
              /* Store number of disulphide bonds */
              iface_ptr->count[PP_DISULPHIDES] = atoi(value);
              if (iface_ptr->count[PP_DISULPHIDES] > 0)
                iface_ptr->ntypes++;
            }
        }
    }

  /* Close the file */
  fclose(file_ptr);

  /* Return pointer to first chain record */
  *fst_chain_ptr = first_chain_ptr;

  /* Return number of chains read in */
  *nchn = nchains;

  /* Return number of interfaces */
  return(niface);
}
/***********************************************************************

gen_pp_plot  -  Run ppcont on the given interface, or on all

***********************************************************************/

void gen_pp_plot(char *directory,char *pdb_code,char *chains,
                 char *log_file,struct c_file_data file_name_ptr,
                 struct parameters params)
{
  char command[LINELEN + 1], prog_params[FILENAME_LEN];
  char new_name[FILENAME_LEN], old_name[FILENAME_LEN];
  char fname[FILENAME_LEN];

  int dos;

  /* Initialise variables */
  dos = file_name_ptr.dos;
  
  /* Define the command-line parameters to generate all-chains plot */
  sprintf(prog_params,"%s%s -pdbsum1",pdb_code,chains);

  /* Show progress */
  printf("          ppcont for chains %s\n",chains);

  /* Run the program */
  sprintf(command,"cd %s && %s/ppcont %s > %s",directory,
          file_name_ptr.other_dir[PDBSUM_EXE_DIR],prog_params,log_file);
  call_system(command,dos);

  /* Form the output file name */
  if (strstr(chains,"all") != NULL)
    strcpy(fname,"iface");
  else
    sprintf(fname,"iface_%s",chains);

  /* Rename the PostScript file */
  sprintf(old_name,"ppcont.ps");
  sprintf(new_name,"%s.ps",fname);
  move_file(directory,old_name,new_name,dos);
  //  sprintf(command,"cd %s && %s ppcont.ps %s.ps",directory,MV_COMMAND,
  //          fname);
  //  call_system(command,dos);

  /* Show progress */
  printf("             Converting files: PDF ");

  /* Convert the PostScript file to PDF */
  sprintf(command,"cd %s && %s %s.ps %s.pdf",directory,
          file_name_ptr.exe_name[PS2PDF],fname,fname);
  call_system(command,dos);

  /* Show progress */
  printf("GIF ");
  fflush(stdout);

  /* Form the input parameters */
  sprintf(prog_params,"-trim +repage -geometry 50%% -density 200x200");

  /* Convert .ps file to GIF */
/* v.1.0.1 --> */
//  sprintf(command,"cd %s && %s %s %s.ps %s.gif",directory,
//          file_name_ptr.exe_name[CONVERT_EXE],prog_params,fname,fname);
  sprintf(command,"cd %s && %s %s.ps %s %s.gif",directory,
          file_name_ptr.exe_name[CONVERT_EXE],fname,prog_params,fname);
/* <-- v.1.0.1 */
  call_system(command,dos);

  /* Close off progress line */
  printf("\n");
  fflush(stdout);
}
/***********************************************************************

gen_iface_plot  -  Run dimhtml on the given interface

***********************************************************************/

void gen_iface_plot(char *directory,char *pdb_code,char *chains,
                    char *log_file,struct c_file_data file_name_ptr,
                    struct parameters params)
{
  char command[LINELEN + 1], prog_params[FILENAME_LEN];
  char new_name[FILENAME_LEN], old_name[FILENAME_LEN];
  char fname[FILENAME_LEN];

  int dos;

  /* Initialise variables */
  dos = file_name_ptr.dos;

  /* Show progress */
  printf("          dimhtml for chains %s\n",chains);

  /* Define the command-line parameters to generate all-chains plot */
  sprintf(prog_params,"%s%s -pdbsum1",pdb_code,chains);

  /* Run the program */
  sprintf(command,"cd %s && %s/dimhtml %s > %s",directory,
          file_name_ptr.other_dir[PDBSUM_EXE_DIR],prog_params,log_file);
  call_system(command,dos);

  /* Form the output file name */
  sprintf(fname,"dimplot_%s",chains);

  /* Rename the PostScript file */
  sprintf(old_name,"dimhtml.ps");
  sprintf(new_name,"%s.ps",fname);
  move_file(directory,old_name,new_name,dos);

  //  sprintf(command,"cd %s && %s dimhtml.ps %s.ps",directory,MV_COMMAND,
  //          fname);
  //  call_system(command,dos);

  /* Show progress */
  printf("             Converting files: PDF ");

  /* Convert the PostScript file to PDF */
  sprintf(command,"cd %s && %s %s.ps %s.pdf",directory,
          file_name_ptr.exe_name[PS2PDF],fname,fname);
  call_system(command,dos);

  /* Show progress */
  printf("GIF ");
  fflush(stdout);

  /* Form the input parameters */
  sprintf(prog_params,"-trim +repage -geometry 40%% -density 200x200");

  /* Convert .ps file to GIF */
/* v.1.0.1 --> */
//  sprintf(command,"cd %s && %s %s %s.ps %s.gif",directory,
//          file_name_ptr.exe_name[CONVERT_EXE],prog_params,fname,fname);
  sprintf(command,"cd %s && %s %s.ps %s %s.gif",directory,
          file_name_ptr.exe_name[CONVERT_EXE],fname,prog_params,fname);
/* <-- v.1.0.1 */
  call_system(command,dos);

  /* Close off progress line */
  printf("\n");
  fflush(stdout);
}
/***********************************************************************

write_iface_txt  -  Output the listing file giving the interactions

***********************************************************************/

void write_iface_txt(char *pdir,char *directory,char *pdb_code,
                     char *chains,char *log_file,
                     struct c_file_data file_name_ptr)
{
  char file_name[FILENAME_LEN], out_file[FILENAME_LEN];
  char input_line[LINELEN + 1];
  char chain1, chain2, inchain1, inchain2, last_type, type;
  char atom_name1[5], atom_name2[5], atom_no1[6], atom_no2[6];
  char res_name1[4], res_name2[4], res_num1[6], res_num2[6];
  char count[20], dist[6];

  int len, ncontacts, ndisulphides, nhbonds, nsalt_bridges;
  int swap;

  FILE *file_ptr, *outfile_ptr;

  /* Initialise variables */
  chain1 = chains[0];
  chain2 = chains[1];
  last_type = ' ';
  nhbonds = 0;
  ncontacts = 0;
  nsalt_bridges = 0;
  ndisulphides = 0;

  /* Show progress */
  printf("          listing for chains %s\n",chains);

  /* Form the output file name */
  sprintf(out_file,"%s/iface_%s.txt",directory,chains);

  /* Open the output file */
  outfile_ptr = open_output_file(out_file,TRUE);

  /* Write out the header records */
  fprintf(outfile_ptr,"List of atom-atom interactions across "
          "protein-protein interface\n");
  fprintf(outfile_ptr,"--------------------------------------"
          "-------------------------\n");
  fprintf(outfile_ptr,"\n");
  fprintf(outfile_ptr,"                 PDB code: %s   Chains %c }{ %c\n",
          pdb_code,chains[0],chains[1]);
  fprintf(outfile_ptr,"                 ------------------------------\n");
  fprintf(outfile_ptr,"\n");

  /* Form name of the igrow.out file */
  sprintf(file_name,"%s/newsec/igrow.out",pdir);
  
  /* Open the input file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    {
      printf("*** Warning. Unable to open file [%s]\n",file_name);
      return;
    }

  /* Loop through the file, line by line */
  while (fgets(input_line,LINELEN,file_ptr) != NULL)
    {
      /* Check line-length */
      len = string_chomp(input_line);

      /* Check that line is long enough */
      if (len > 60)
        {
          /* Extract the data from the line */
          type = input_line[0];
          inchain1 = input_line[6];
          strncpy(res_num1,input_line + 9, 5);
          res_num1[5] = '\0';
          strncpy(res_name1,input_line + 14, 3);
          res_name1[3] = '\0';
          strncpy(atom_no1,input_line + 18, 5);
          atom_no1[5] = '\0';
          strncpy(atom_name1,input_line + 24, 4);
          atom_name1[4] = '\0';
          strncpy(atom_name2,input_line + 33, 4);
          atom_name2[4] = '\0';
          strncpy(atom_no2,input_line + 38, 5);
          atom_no2[5] = '\0';
          strncpy(res_name2,input_line + 44, 3);
          res_name2[3] = '\0';
          strncpy(res_num2,input_line + 49, 5);
          res_num2[5] = '\0';
          inchain2 = input_line[54];
          strncpy(dist,input_line + 56, 5);
          dist[5] = '\0';

          /* If interaction is between required chains, process */
          if ((inchain1 == chain1 && inchain2 == chain2) ||
              (inchain2 == chain1 && inchain1 == chain2))
            {
              /* If type has changed, write heading */
              if (type != last_type)
                {
                  /* Salt bridges */
                  if (type == 'S')
                    {
                      fprintf(outfile_ptr,"\n");
                      fprintf(outfile_ptr,"Salt bridges\n");
                      fprintf(outfile_ptr,"------------\n");
                      fprintf(outfile_ptr,"\n");
                    }

                  /* Disulphide bonds */
                  else if (type == 'D')
                    {
                      fprintf(outfile_ptr,"\n");
                      fprintf(outfile_ptr,"Disulphide bonds\n");
                      fprintf(outfile_ptr,"----------------\n");
                      fprintf(outfile_ptr,"\n");
                    }

                  /* H-bonds */
                  else if (type == 'H')
                    {
                      fprintf(outfile_ptr,"\n");
                      fprintf(outfile_ptr,"Hydrogen bonds\n");
                      fprintf(outfile_ptr,"--------------\n");
                      fprintf(outfile_ptr,"\n");
                    }

                  /* Non-bonded contacts */
                  else
                    {
                      fprintf(outfile_ptr,"\n");
                      fprintf(outfile_ptr,"Non-bonded contacts\n");
                      fprintf(outfile_ptr,"-------------------\n");
                      fprintf(outfile_ptr,"\n");
                    }

                  /* Column headings */
                  fprintf(outfile_ptr,"       <----- A T O M   1 ----->"
                          "<----- A T O M   2 ----->\n");
                  fprintf(outfile_ptr,"\n");
                  fprintf(outfile_ptr,"       Atom Atom Res  Res            "
                          "  Atom Atom Res  Res\n");
                  fprintf(outfile_ptr,"        no. name name no.  Chain     "
                          "   no. name name no.  Chain  Distance\n");

                  /* Save interaction type */
                  last_type = type;
                }

              /* Determine whether atoms need to be swapped */
              if (inchain2 == chain1 && inchain1 == chain2)
                {
                  swap = TRUE;
                }
              else
                {
                  swap = FALSE;
                }

              /* Increment count of this interaction type */
              if (type == 'S')
                {
                  nsalt_bridges++;
                  sprintf(count,"%3d.",nsalt_bridges);
                  }
              else if (type == 'D')
                {
                    ndisulphides++;
                    sprintf(count,"%3d.",ndisulphides);
                }
              else if (type == 'H')
                {
                  nhbonds++;
                  sprintf(count,"%3d.",nhbonds);
                }
              else
                {
                  ncontacts++;
                  sprintf(count,"%3d.",ncontacts);
                  }

              /* Print the current line */
              if (swap == FALSE)
                {
                  fprintf(outfile_ptr,"%s  %s %s %s %s   %c   <--> ",
                          count,atom_no1,atom_name1,res_name1,res_num1,
                          inchain1);
                  fprintf(outfile_ptr,"%s %s %s %s   %c     %s\n",
                          atom_no2,atom_name2,res_name2,res_num2,
                          inchain2,dist);
                  }
                else
                  {
                    fprintf(outfile_ptr,"%s  %s %s %s %s   %c   <--> ",
                            count,atom_no2,atom_name2,res_name2,res_num2,
                            inchain2);
                    fprintf(outfile_ptr,"%s %s %s %s   %c     %s\n",
                            atom_no1,atom_name1,res_name1,res_num1,
                            inchain1,dist);
                  }
              }
          }
      }

  /* Write counts of interactions */
  if (nsalt_bridges > 0)
    fprintf(outfile_ptr,"Number of salt bridges:          %3d\n",
            nsalt_bridges);
  if (ndisulphides > 0)
    fprintf(outfile_ptr,"Number of disulphide bonds:      %3d\n",
            ndisulphides);
  if (nhbonds > 0)
    fprintf(outfile_ptr,"Number of hydrogen bonds:        %3d\n",
            nhbonds);
  if (ncontacts > 0)
    fprintf(outfile_ptr,"Number of non-bonded contacts:   %3d\n",
            ncontacts);

  /* Close both files */
  fclose(file_ptr);
  fclose(outfile_ptr);
}
/***********************************************************************

run_ppcont  -  Run ppcont and dimhtml to generate the interface plots

***********************************************************************/

void run_ppcont(char *fno,int process,char *pdir,char *pdb_code,
                char *file_name,char *file_new,char *log_dir,
                struct chain *first_chain_ptr,
                struct c_file_data file_name_ptr,struct parameters params)
{
  char command[LINELEN + 1], directory[FILENAME_LEN];
  char chains[6], log_file[FILENAME_LEN];

  int len;
  int have_file, dos;

  float resol;

  struct chain *chain_ptr;
  struct iface *iface_ptr;

  /* Initialise variables */
  dos = file_name_ptr.dos;

  /* Show progress */
  printf("   %s%d. In run_ppcont ...\n",fno,process);

  /* Create the protprot directory */
  sprintf(directory,"%s/protprot",pdir);
  create_directory(directory);

  /* Form name of the log file */
  sprintf(log_file,"%s/ppcont.log",log_dir);

  /* Run the all-chains ppcont plot */ 
  strcpy(chains," -all");
  gen_pp_plot(directory,pdb_code,chains,log_file,file_name_ptr,params);

  /* Get first record */
  chain_ptr = first_chain_ptr;

  /* Loop through all the records */
  while (chain_ptr != NULL)
    {
      /* Get first interface record */
      iface_ptr = chain_ptr->first_iface_ptr;

      /* Loop through all the records */
      while (iface_ptr != NULL)
        {
          /* Join the two chains */
          sprintf(chains,"%c%c",chain_ptr->chain,
                  iface_ptr->chain_ptr->chain);

          /* Run ppcont for this interface */
          gen_pp_plot(directory,pdb_code,chains,log_file,
                      file_name_ptr,params);

          /* Run dimhtml for this interface */
          gen_iface_plot(directory,pdb_code,chains,log_file,
                         file_name_ptr,params);

          /* Output the listing file giving the interactions */
          write_iface_txt(pdir,directory,pdb_code,chains,log_file,
                          file_name_ptr);

          /* Get pointer to the next record */
          iface_ptr = iface_ptr->next_iface_ptr;
        }  

      /* Get pointer to the next record */
      chain_ptr = chain_ptr->next_chain_ptr;
    }  
}
/***********************************************************************

read_chain_equivs  -  Read in the chain equivalences from the chains.dat
                      file

***********************************************************************/

int read_chain_equivs(char *directory,
                      struct r_rasdef *first_r_rasdef_ptr,
                      struct chain **fst_chain_ptr,
                      struct chain **lst_chain_ptr,int dos)
{
  char file_name[FILENAME_LEN + 1], input_line[LINELEN + 1];
  char chain_id, copy_of;

  int len, nchains, nlines, nunique;
  int found;

  struct chain *chain_ptr, *first_chain_ptr, *last_chain_ptr;
  struct r_rasdef *r_rasdef_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  nchains = 0;
  nlines = 0;
  nunique = 0;

  /* Initialise chain pointers */
  chain_ptr = first_chain_ptr = last_chain_ptr = NULL;

  /* Form the name of the chains.dat file */
  sprintf(file_name,"%s/chains.dat",directory);

  /* For DOS, change the directory separators */
  if (dos == TRUE)
    replace_char(file_name,'/','\\');

  /* Open the chains.dat input file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    {
      printf("\n*** Warning. Unable to open chains.dat file [%s]\n",
             file_name);
      return(nchains);
    }

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

      /* If this is a copy chain, mark it as such in the corresponding
         rasdef record */

      /* Get first rasdef record */
      r_rasdef_ptr = first_r_rasdef_ptr;
      found = FALSE;
      if (copy_of == '\0')
        found = TRUE;

      /* Loop through the rasdefs until done */
      while (r_rasdef_ptr != NULL && found == FALSE)
        {
          /* If this is a protein chain, then process */
          if (r_rasdef_ptr->type == R_PROTEIN)
            {
              /* If this is the copy chain, then mark it */
              if (chain_id == r_rasdef_ptr->chain)
                {
                  /* Set copy number to 1 */
                  r_rasdef_ptr->copy = 1;

                  /* Set flag that found */
                  found = TRUE;
                }
            }

          /* Get pointer to the next rasdef record */
          r_rasdef_ptr = r_rasdef_ptr->next_r_rasdef_ptr;
        }  
    }

  /* Get first rasdef record */
  r_rasdef_ptr = first_r_rasdef_ptr;

  /* Loop through the rasdefs to create a chain record form each 
     protein chain */
  while (r_rasdef_ptr != NULL)
    {
      /* If this is a protein chain, then create chain */
      if (r_rasdef_ptr->type == R_PROTEIN)
        {
          /* Create chain record */
          chain_ptr
            = create_chain_record(&first_chain_ptr,&last_chain_ptr,
                                  r_rasdef_ptr->chain,nchains);

          /* Increment count of new chains */
          nchains++;
        }

      /* Get pointer to the next rasdef record */
      r_rasdef_ptr = r_rasdef_ptr->next_r_rasdef_ptr;
    }  

  /* Close the input file */
  fclose(file_ptr);

  /* Return current pointers */
  *fst_chain_ptr = first_chain_ptr;
  *lst_chain_ptr = last_chain_ptr;

  /* Return the number of chains stored */
  return(nchains);
}
/***********************************************************************

count_chains  -  Count protein and DNA chains

***********************************************************************/

int count_chains(struct r_rasdef *first_r_rasdef_ptr,int *nn)
{
  int nnucl, nprot;

  struct r_rasdef *r_rasdef_ptr;

  /* Initialise variables */
  *nn = nnucl = nprot = 0;

  /* Get first rasdef record */
  r_rasdef_ptr = first_r_rasdef_ptr;

  /* Loop through the rasdefs to count chain types */
  while (r_rasdef_ptr != NULL)
    {
      /* If this is a protein chain, then increment count */
      if (r_rasdef_ptr->type == R_PROTEIN)
        nprot++;

      /* Ditto for DNA/RNA chain */
      else if (r_rasdef_ptr->type == R_NUCLEIC_ACID)
        nnucl++;

      /* Get pointer to the next rasdef record */
      r_rasdef_ptr = r_rasdef_ptr->next_r_rasdef_ptr;
    }  

  /* Return the number of DNA/RNA chains */
  *nn = nnucl;

  /* Return the number of protein chains */
  return(nprot);
}
/***********************************************************************

run_promotif  -  Run the remainder of the PROMOTIF programs

***********************************************************************/

void run_promotif(char *fno,int process,char *pdir,char *pdb_code,
                  char *file_name,char *file_new,char *log_dir,
                  struct c_file_data file_name_ptr,
                  struct parameters params)
{
  char command[LINELEN + 1], directory[FILENAME_LEN];
  char log_file[FILENAME_LEN], prog_params[FILENAME_LEN];
  char in_file[FILENAME_LEN], out_file[FILENAME_LEN];
  char param_dir[FILENAME_LEN];

  int len;
  int ok, dos;

  float resol;

  /* Initialise variables */
  ok = TRUE;
  strcpy(param_dir,file_name_ptr.other_dir[PDBSUM_PARAMS_DIR]);
  dos = file_name_ptr.dos;

  /* Show progress */
  printf("   %s%d. In run_promotif ...\n",fno,process);

  /* Define directory */
  sprintf(directory,"%s/promotif",pdir);
  create_directory(directory);

  /* Copy across the required data files */
  sprintf(command,"%s %s/phipsi.mat %s",CP_COMMAND,param_dir,
          directory);
  call_system(command,dos);
  sprintf(command,"%s %s/promotif2.prm %s",CP_COMMAND,param_dir,
          directory);
  call_system(command,dos);
  sprintf(command,"%s %s/hera_colours %s",CP_COMMAND,param_dir,
          directory);
  call_system(command,dos);

  /* == SSTRUC == */

  /* Show progress */
  printf("          - Running sstruc ...\n");

  /* Form name of the log file */
  sprintf(log_file,"%s/promotif.log",log_dir);

  /* Define input parameter file */
  sprintf(in_file,"%s/in",directory);

  /* Write out the input parameters */
  sprintf(command,"../%s\n",file_new);
  write_in_file(directory,command);

  //  sprintf(command,"echo ../%s > %s",file_new,in_file);
  //  call_system(command,dos);

  /* Run the program */
  sprintf(command,"cd %s && %s/p_sstruc3 < %s > %s",directory,
          file_name_ptr.other_dir[PDBSUM_EXE_DIR],in_file,log_file);
  call_system(command,dos);

  /* == TURN == */

  /* Show progress */
  printf("          - Running turn ...\n");

  /* Run the program */
  sprintf(command,"cd %s && %s/p_turn3 < %s >> %s",directory,
          file_name_ptr.other_dir[PDBSUM_EXE_DIR],in_file,log_file);
  call_system(command,dos);

  /* == HAIRPIN == */

  /* Show progress */
  printf("          - Running hairpin ...\n");

  /* Run the program */
  sprintf(command,"cd %s && %s/p_hairpin3 < %s >> %s",directory,
          file_name_ptr.other_dir[PDBSUM_EXE_DIR],in_file,log_file);
  call_system(command,dos);

  /* == BULGE == */

  /* Show progress */
  printf("          - Running bulge ...\n");

  /* Run the program */
  sprintf(command,"cd %s && %s/p_bulge3 < %s >> %s",directory,
          file_name_ptr.other_dir[PDBSUM_EXE_DIR],in_file,log_file);
  call_system(command,dos);

  /* == SHEET == */

  /* Show progress */
  printf("          - Running sheet ...\n");

  /* Run the program */
  sprintf(command,"cd %s && %s/p_sheet3 < %s >> %s",directory,
          file_name_ptr.other_dir[PDBSUM_EXE_DIR],in_file,log_file);
  call_system(command,dos);

  /* == HELIX == */

  /* Show progress */
  printf("          - Running helix ...\n");

  /* Run the program */
  sprintf(command,"cd %s && %s/p_helix3 < %s >> %s",directory,
          file_name_ptr.other_dir[PDBSUM_EXE_DIR],in_file,log_file);
  call_system(command,dos);

  /* == DISULPH == */

  /* Show progress */
  printf("          - Running disulph ...\n");

  /* Run the program */
  sprintf(command,"cd %s && %s/p_disulph3 < %s >> %s",directory,
          file_name_ptr.other_dir[PDBSUM_EXE_DIR],in_file,log_file);
  call_system(command,dos);

  /* == PROMFILT == */

  /* Show progress */
  printf("          - Running promfilt ...\n");

  /* Run the program */
  sprintf(command,"cd %s && %s/promfilt %s -pdbsum1 >> %s",directory,
          file_name_ptr.other_dir[PDBSUM_EXE_DIR],pdb_code,log_file);
  call_system(command,dos);
}
/***********************************************************************

create_plot_entry  -  Create a plot entry

***********************************************************************/

struct plot *create_plot_entry(struct plot **fst_plot_ptr,
                               struct plot **lst_plot_ptr,
                               char chain_no)
{
  int i;

  struct plot *plot_ptr;
  struct plot *first_plot_ptr, *last_plot_ptr;

  /* Initialise pointers */
  plot_ptr = NULL;
  first_plot_ptr = *fst_plot_ptr;
  last_plot_ptr = *lst_plot_ptr;

  /* Create entry */
  plot_ptr = (struct plot *) malloc(sizeof(struct plot));
  if (plot_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for struct plot\n");
      printf("***        Program terminated with error.\n");
      exit (1);
    }

  /* Store the known info */
  plot_ptr->chain_no = chain_no;

  /* Initialise remaining entries */
  plot_ptr->chain = ' ';
  for (i = 0; i < NPLOTS; i++)
    plot_ptr->count[i] = 0;

  /* Initialise pointer to next entry */
  plot_ptr->next_plot_ptr = NULL;

  /* If this is the first entry stored, then update
     pointer to head of linked list */
  if (first_plot_ptr == NULL)
    first_plot_ptr = plot_ptr;

  /* Otherwise, make previous entry point to the current
     one */
  else
    last_plot_ptr->next_plot_ptr = plot_ptr;
                      
  /* Save the current residue's pointer */
  last_plot_ptr = plot_ptr;

  /* Return current pointers */
  *fst_plot_ptr = first_plot_ptr;
  *lst_plot_ptr = last_plot_ptr;
  return(plot_ptr);
}
/***********************************************************************

get_motifs  -  Get the motifs defined for each chain

***********************************************************************/

int get_motifs(char *file_name,struct plot ***plot_idx_ptr)
{
  char input_line[LINELEN + 1], number_string[LINELEN + 1];

  char *string_ptr;
  
  int chain_no, count, ichain, len, nchains, plot_no;

  struct plot *plot_ptr, *first_plot_ptr, *last_plot_ptr;

  struct plot **plot_index_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  ichain = 0;
  nchains = 0;
  plot_no = 0;

  /* Initialise pointers */
  first_plot_ptr = last_plot_ptr = NULL;

  /* Create the index array of records per possible chain */
  plot_index_ptr
    = (struct plot **) malloc(MAX_CHAINS * sizeof(struct plot *));

  /* Check that sufficient memory was allocated */
  if (plot_index_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for index array\n");
      exit(1);
    }

  /* Loop to initialise all the pointers */
  for (ichain = 0; ichain < MAX_CHAINS; ichain++)
    plot_index_ptr[ichain] = NULL;

  /* Open the motifs.dat file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    {
      printf("\n*** Warning. Unable to open file [%s]\n",file_name);
      return(nchains);
    }

  /* Loop through the file, line by line */
  while (fgets(input_line,LINELEN,file_ptr) != NULL)
    {
      /* Initialise variables */
      chain_no = ichain = 0;
      plot_ptr = NULL;

      /* Truncate the input line */
      len = string_chomp(input_line);

      /* If plot type, then extract details */
      if (!strncmp(input_line,"P_TYPE[",7))
        {
          /* Get the chain number */
          strcpy(number_string,input_line + 7);

          /* Get the close bracket */
          string_ptr = strchr(number_string,']');
          if (string_ptr != NULL)
            {
              string_ptr[0] = '\0';
              chain_no = atoi(number_string);
              ichain = chain_no - 1;
            }

          /* Get the plot number */
          string_ptr = strchr(input_line,' ');
          if (string_ptr != NULL)
            plot_no = atoi(string_ptr + 1);

          /* If have chain number and plot number, then save */
          if (chain_no > 0 && plot_no > 0)
            {
              /* Get the corresponding plot record */
              plot_ptr = plot_index_ptr[ichain];

              /* If we haven't created a plot record for this chain,
                 create it now */
              if (plot_ptr == NULL)
                {
                  /* Create a new plot entry */
                  plot_ptr
                    = create_plot_entry(&first_plot_ptr,&last_plot_ptr,
                                        chain_no);

                  /* Save this pointer */
                  plot_index_ptr[ichain] = plot_ptr;

                  /* Increment count of chain records created */
                  nchains++;
                }
            }
        }

      /* If count, then store */
      else if (!strncmp(input_line,"P_COUNT[",8))
        {
          /* Get the chain number */
          strcpy(number_string,input_line + 8);

          /* Get the close bracket */
          string_ptr = strchr(number_string,']');
          if (string_ptr != NULL)
            {
              string_ptr[0] = '\0';
              chain_no = atoi(number_string);
              ichain = chain_no - 1;
            }

          /* Get the count */
          string_ptr = strchr(input_line,' ');
          if (string_ptr != NULL)
            count = atoi(string_ptr + 1);

          /* Get the corresponding plot record */
          plot_ptr = plot_index_ptr[ichain];

          /* If OK, then save the count */
          if (plot_ptr != NULL)
            {
              /* Save the count */
              plot_ptr->count[plot_no - 1] = count;
            }
        }
    }

  /* Close file */
  fclose(file_ptr);

  /* Return the index of plots */
  *plot_idx_ptr = plot_index_ptr;

  /* Return number of chains read in */
  return(nchains);
}
/***********************************************************************

run_psplot  -  Run psplot to generate each of the required plots

***********************************************************************/

int run_psplot(char *fno,int process,char *pdir,char *pdb_code,
               char *file_name,char *file_new,char *log_dir,
               struct r_rasdef *first_r_rasdef_ptr,
               struct plot ***plot_idx_ptr,
               struct c_file_data file_name_ptr,
               struct parameters params)
{
  char command[LINELEN + 1], directory[FILENAME_LEN];
  char log_file[FILENAME_LEN], prog_params[FILENAME_LEN];
  char in_file[FILENAME_LEN], gif_file[FILENAME_LEN];
  char chain, fname[FILENAME_LEN], param_dir[FILENAME_LEN];
  char new_name[FILENAME_LEN], old_name[FILENAME_LEN];

  int chain_no, count, ichain, iplot, len, nchains, plot_no;
  int have_file, ok, dos;

  struct plot *plot_ptr;
  struct r_rasdef *r_rasdef_ptr;

  struct plot **plot_index_ptr;

  FILE *file_ptr, *outfile_ptr;

  /* Initialise variables */
  chain_no = 0;
  dos = file_name_ptr.dos;
  ok = TRUE;
  strcpy(param_dir,file_name_ptr.other_dir[PDBSUM_PARAMS_DIR]);

  /* Show progress */
  printf("   %s%d. In run_psplot ...\n",fno,process);

  /* Form the name of the motif.dat file */
  sprintf(file_name,"%s/promotif/motifs.dat",pdir);

  /* Read in all the motifs that we have */
  nchains = get_motifs(file_name,&plot_index_ptr);

  /* Save the index of plots */
  *plot_idx_ptr = plot_index_ptr;

  /* If no chain data read in, then have no plots to generate */
  if (nchains == 0)
    return(nchains);

  /* Define plot directory */
  sprintf(directory,"%s/psplots",pdir);
  create_directory(directory);

  /* Copy across the required data files */
  sprintf(command,"%s %s/promotif2.prm %s",CP_COMMAND,param_dir,
          directory);
  call_system(command,dos);

  /* Copy CATHPARAM file into this directory - no longer required */
  //sprintf(command,"%s %s/CATHPARAM %s",CP_COMMAND,param_dir,
  //        directory);
  //call_system(command,dos);

  /* Form name of the log file */
  sprintf(log_file,"%s/psplot.log",log_dir);

  /* Define input parameter file */
  sprintf(in_file,"%s/in",directory);

  /* Get first rasdef record */
  r_rasdef_ptr = first_r_rasdef_ptr;

  /* Loop through the rasdefs until done */
  while (r_rasdef_ptr != NULL)
    {
      /* If this is a protein chain, then process */
      if (r_rasdef_ptr->type == R_PROTEIN && r_rasdef_ptr->copy == 0)
        {
          /* Get the chain and the chain number */
          chain = r_rasdef_ptr->chain;
          chain_no++;
          ichain = chain_no - 1;

          /* Show progress */
          printf("           Chain %c : ",chain);
          fflush(stdout);

          /* Get the corresponding plot record */
          plot_ptr = plot_index_ptr[ichain];

          /* Save the chain id */
          plot_ptr->chain = chain;

          /* Generate all the plots for this chain */
          for (iplot = 0; iplot < NPLOTS && plot_ptr != NULL; iplot++)
            {
              /* If have a non-zero count for this motif, and Promotif
                 produces PostScript plot for it, then generate */
              if (plot_ptr->count[iplot] > 0 && 
                  !strcmp(HAVE_PS[iplot],"Yes"))
                {
                  /* Show progress */
                  printf("%s ",PLOT_DESC[iplot]);
                  fflush(stdout);

                  /* Get the number of PostScript pages */
                  count = plot_ptr->count[iplot];

                  /* Open the in file */
                  outfile_ptr = open_output_file(in_file,TRUE);

                  /* Write out the input parameters */
                  fprintf(outfile_ptr,"%s\n",pdb_code);
                  fprintf(outfile_ptr,"%c\n",chain);
                  fprintf(outfile_ptr,"%s\n",PLOT_NAME[iplot]);
                  fprintf(outfile_ptr," \n");
                  fprintf(outfile_ptr,"1\n");
                  fprintf(outfile_ptr,"ONE\n");
                  fprintf(outfile_ptr,"PDBsum1\n");

                  /* Close the file */
                  fclose(outfile_ptr);

                  /* Run the program */
                  sprintf(command,"cd %s && %s/psplot3 < %s > %s",
                          directory,
                          file_name_ptr.other_dir[PDBSUM_EXE_DIR],
                          in_file,log_file);
                  call_system(command,dos);

                  /* Form the file name */
                  sprintf(fname,"%s_%c",PLOT_DESC[iplot],chain);

                  /* Check if file produced */
                  sprintf(file_name,"%s/promotif001.ps",directory);
                  have_file = check_file(file_name);

                  /* If file found, then rename and convert */
                  if (have_file == TRUE)
                    {
                      /* Form the name of the destination file and see
                         if it exists */
                      sprintf(new_name,"%s/%s.ps",directory,fname);
                      have_file = check_file(new_name);

                      /* If file exists, then delete it */
                      if (have_file == TRUE)
                        {
                          /* For the delete command and delete file */
                          sprintf(command,"%s %s",RM_COMMAND,new_name);
                          call_system(command,dos);
                        }

                      /* Rename the PostScript file */
                      sprintf(old_name,"promotif001.ps");
                      sprintf(new_name,"%s.ps",fname);
                      move_file(directory,old_name,new_name,dos);
                      
                      // sprintf(command,"cd %s && %s promotif001.ps %s.ps",
                      //              directory,MV_COMMAND,fname);
                      // call_system(command,dos);

                      /* Convert the PostScript file to PDF */
                      sprintf(command,"cd %s && %s %s.ps %s.pdf",directory,
                              file_name_ptr.exe_name[PS2PDF],fname,fname);
                      call_system(command,dos);
                    }
                }
            }

          /* Close off progress line */
          printf("\n");
          fflush(stdout);
        }

      /* Get pointer to the next rasdef record */
      r_rasdef_ptr = r_rasdef_ptr->next_r_rasdef_ptr;
    }

  /* Return number of chains with Promotif plots */
  return(nchains);
}
/***********************************************************************

run_resdata  -  Run the remainder of the RESDATA programs

***********************************************************************/

void run_resdata(char *fno,int process,char *pdir,char *pdb_code,
                 char *file_name,char *file_new,char *log_dir,
                 struct c_file_data file_name_ptr,
                 struct parameters params)
{
  char command[LINELEN + 1], directory[FILENAME_LEN];
  char log_file[FILENAME_LEN];

  int len;
  int ok, dos;

  float resol;

  /* Initialise variables */
  ok = TRUE;
  dos = file_name_ptr.dos;

  /* Show progress */
  printf("   %s%d. In run_resdata ...\n",fno,process);

  /* Define directory */
  sprintf(directory,"%s/newsec",pdir);
  create_directory(directory);

  /* Form name of the log file */
  sprintf(log_file,"%s/resdata.log",log_dir);

  /* Run the program */
  sprintf(command,"cd %s && %s/resdata %s -pdbsum1 > %s",directory,
          file_name_ptr.other_dir[PDBSUM_EXE_DIR],pdb_code,log_file);
  call_system(command,dos);
}
/***********************************************************************

generate_topology  -  Run Hera and generate the topology diagram

***********************************************************************/

void generate_topology(char *pdb_file,char *pdb_code,char chain,
                       int domain,char *directory,char *pdbsum_exe_dir,
                       char *log_file,char *ps2pdf,char *convert_exe,
                       struct c_file_data file_name_ptr)
{
#define NTERMS       11

  char command[LINELEN + 1], domain_file[FILENAME_LEN];
  char in_file[FILENAME_LEN], prog_params[FILENAME_LEN];
  char domain_id[FILENAME_LEN], file_name[FILENAME_LEN];
  char new_name[FILENAME_LEN], old_name[FILENAME_LEN];

  static char *TERM[NTERMS]
    = { "1", "1", "1", "2", "1", "1", "1", "0", "0", "0", "Y" };

  int i;
  int have_file, dos;

  FILE *outfile_ptr;

  /* Initialise variables */
  dos = file_name_ptr.dos;

  /* Form name of the domain file for this chain and domain */
  sprintf(domain_file,"../newsec/dom%c%2.2d.pdb",chain,domain);

  /* == SSTRUC == */

  /* Show progress */
  printf("          - Running sstruc ...\n");

  /* Define input parameter file */
  sprintf(in_file,"%s/in",directory);

  /* Write out the input parameters */
  sprintf(command,"%s\n",domain_file);
  write_in_file(directory,command);

  //  sprintf(command,"echo %s > %s",domain_file,in_file);
  //  call_system(command,dos);

  /* Run the program */
  sprintf(command,"cd %s && %s/p_sstruc3 < %s > %s",directory,
          pdbsum_exe_dir,in_file,log_file);
  call_system(command,dos);

  /* == HERA == */

  /* Show progress */
  printf("          - Running hera ...\n");

  /* Define input parameter file */
  sprintf(in_file,"%s/hera.in",directory);

  /* Open the output file */
  outfile_ptr = open_output_file(in_file,TRUE);

  /* Loop to write out the Hera input parameters */
  for (i = 0; i < NTERMS; i++)
    fprintf(outfile_ptr,"%s\n",TERM[i]);

  /* Add the chain is and the file name */
  fprintf(outfile_ptr,"%c\n",chain);
  fprintf(outfile_ptr,"%s\n",domain_file);

  /* Close the file */
  fclose(outfile_ptr);

  /* Run the program */
  sprintf(command,"cd %s && %s/p_hera3 < %s >> %s",directory,
          pdbsum_exe_dir,in_file,log_file);
  call_system(command,dos);

  /* Get the filename stub */
  sprintf(domain_id,"%c%2.2d",chain,domain);

  /* Rename the output .her file to .ps */
  sprintf(old_name,"dom%s.her",domain_id);
  sprintf(new_name,"hera%s.ps",domain_id);
  move_file(directory,old_name,new_name,dos);

  /* Convert .ps file to PDF */
  if (strcmp(ps2pdf,"NONE"))
    {
      sprintf(command,"cd %s && %s hera%s.ps hera%s.pdf",directory,ps2pdf,
              domain_id,domain_id);
      call_system(command,dos);
    }

  /* == TOPOLDIA == */

  /* Show progress */
  printf("          - Running topoldia ...\n");

  /* Form the input parameters */
  sprintf(prog_params,"%s -save -resno -chain %c -domain %d -pdbsum1",
          pdb_code,chain,domain);

  /* Run the program */
  sprintf(command,"cd %s && %s/topoldia %s >> %s",directory,pdbsum_exe_dir,
          prog_params,log_file);
  call_system(command,dos);

  /* Show progress */
  printf("               ... Converting files\n");
  fflush(stdout);

  /* Check whether the output file was created */
  sprintf(file_name,"%s/topoldia.ps",directory);
  have_file = check_file(file_name);

  /* If have file, then convert */
  if (have_file == TRUE)
    {
      /* Show progress */
      printf("                   ");
      fflush(stdout);

      /* Form the new name for the topoldia output file */
      sprintf(file_name,"topoldia_%c_%2.2d",chain,domain);

      /* Copy the PostScript file */
      sprintf(command,"%s %s/topoldia.ps %s/%s.ps",
              CP_COMMAND,directory,directory,file_name);
      call_system(command,dos);

      /* Convert .ps file to PDF */
      if (strcmp(ps2pdf,"NONE"))
        {
          sprintf(command,"cd %s && %s topoldia.ps %s.pdf",
                  directory,ps2pdf,file_name);
          call_system(command,dos);
        }

      /* Form the parameters for the GIF image */
      sprintf(prog_params,"-density 200x200 -trim +repage");

      /* Generate the GIF image */
/* v.1.0.1 --> */
//      sprintf(command,"cd %s && %s %s %s.ps %s.gif",directory,
//              convert_exe,prog_params,file_name,file_name);
      sprintf(command,"cd %s && %s %s.ps %s %s.gif",directory,
              convert_exe,file_name,prog_params,file_name);
/* <-- v.1.0.1 */
      call_system(command,dos);

      /* Form the parameters for the mini GIF image */
      sprintf(prog_params,
              "-trim +repage -geometry 150x150 -density 200x200");

      /* Generate the mini image */
/* v.1.0.1 --> */
//      sprintf(command,"cd %s && %s %s %s.ps %s_small.png",directory,
//              convert_exe,prog_params,file_name,file_name);
      sprintf(command,"cd %s && %s %s.ps %s %s_small.png",directory,
              convert_exe,file_name,prog_params,file_name);
/* <-- v.1.0.1 */
      call_system(command,dos);

      /* Generate the mini GIF image */
      //      sprintf(command,"cd %s && %s %s %s.gif %s_small.gif",directory,
      //        convert_exe,prog_params,file_name,file_name);
      //call_system(command,dos);

      /* Rename the output topoldia.ps file to domX00.ps */
      sprintf(old_name,"topoldia.ps");
      sprintf(new_name,"dom%s.ps",domain_id);
      move_file(directory,old_name,new_name,dos);

      /* Repeat for the .plt file */
      sprintf(old_name,"topoldia.plt");
      sprintf(new_name,"dom%s.plt",domain_id);
      move_file(directory,old_name,new_name,dos);

      /* Convert .ps file to PDF */
      if (strcmp(ps2pdf,"NONE"))
        {
          /* Show progress */
          printf("PDF ");
          fflush(stdout);

          sprintf(command,"cd %s && %s dom%s.ps dom%s.pdf",directory,
                  ps2pdf,domain_id,domain_id);
          call_system(command,dos);
        }

      /* Form the input parameters */
      sprintf(prog_params,"-density 200x200 -trim +repage");

      /* Show progress */
      printf("GIF1 ");
      fflush(stdout);

      /* Convert .ps file to GIF */
/* v.1.0.1 --> */
//      sprintf(command,"cd %s && %s %s dom%s.ps dom%s.gif",directory,
//              convert_exe,prog_params,domain_id,domain_id);
      sprintf(command,"cd %s && %s dom%s.ps %s dom%s.gif",directory,
              convert_exe,domain_id,prog_params,domain_id);
/* <-- v.1.0.1 */
      call_system(command,dos);

      /* Repeat to create a mini version for the web page */

      /* Form the input parameters */
      sprintf(prog_params,"-trim +repage -geometry 150x150 -density 200x200");

      /* Show progress */
      printf("GIF2 ");
      fflush(stdout);

      /* Convert .ps file to GIF */
/* v.1.0.1 --> */
//      sprintf(command,"cd %s && %s %s dom%s.ps dom%s_small.gif",directory,
//              convert_exe,prog_params,domain_id,domain_id);
      sprintf(command,"cd %s && %s dom%s.ps %s dom%s_small.gif",directory,
              convert_exe,domain_id,prog_params,domain_id);
/* <-- v.1.0.1 */
      call_system(command,dos);

      /* Close off progress line */
      printf("\n");
      fflush(stdout);
    }

  /* Otherwise, show warning message */
  else
    printf("                   *** File missing: topoldia.ps\n");
}
/***********************************************************************

run_hera  -  Generate the hera diagrams for all the protein chains

***********************************************************************/

void run_hera(char *fno,int process,char *pdir,char *pdb_code,
              char *file_name,char *file_new,char *log_dir,
              struct r_rasdef *first_r_rasdef_ptr,
              struct c_file_data file_name_ptr)
{
  char command[LINELEN + 1], directory[FILENAME_LEN];
  char log_file[FILENAME_LEN], param_dir[FILENAME_LEN];

  int domain, len;
  int ok, dos;

  struct r_rasdef *r_rasdef_ptr;

  /* Initialise variables */
  domain = 1;
  ok = TRUE;
  strcpy(param_dir,file_name_ptr.other_dir[PDBSUM_PARAMS_DIR]);
  dos = file_name_ptr.dos;

  /* Show progress */
  printf("   %s%d. In run_hera ...\n",fno,process);

  /* Define directory */
  sprintf(directory,"%s/hera",pdir);
  create_directory(directory);

  /* Copy across the required data files */
  sprintf(command,"%s %s/phipsi.mat %s",CP_COMMAND,param_dir,
          directory);
  call_system(command,dos);
  sprintf(command,"%s %s/promotif2.prm %s",CP_COMMAND,param_dir,
          directory);
  call_system(command,dos);
  sprintf(command,"%s %s/hera_colours %s",CP_COMMAND,param_dir,
          directory);
  call_system(command,dos);

  /* Form name of the log file */
  sprintf(log_file,"%s/hera.log",log_dir);

  /* Get first rasdef record */
  r_rasdef_ptr = first_r_rasdef_ptr;

  /* Loop through the rasdefs until done */
  while (r_rasdef_ptr != NULL)
    {
      /* If this is a protein chain, then process */
      if (r_rasdef_ptr->type == R_PROTEIN && r_rasdef_ptr->copy == 0)
        {
          /* Run hera and generate the topology diagram for this
             chain */
          generate_topology(file_new,pdb_code,r_rasdef_ptr->chain,
                            domain,directory,
                            file_name_ptr.other_dir[PDBSUM_EXE_DIR],
                            log_file,file_name_ptr.exe_name[PS2PDF],
                            file_name_ptr.exe_name[CONVERT_EXE],
                            file_name_ptr);
        }

      /* Get pointer to the next rasdef record */
      r_rasdef_ptr = r_rasdef_ptr->next_r_rasdef_ptr;
    }  
}
/***********************************************************************

run_wiring  -  Generate the wiring diagrams for all the protein chains

***********************************************************************/

void run_wiring(char *fno,int process,char *pdir,char *pdb_code,
                char *file_name,char *file_new,char *log_dir,
                struct r_rasdef *first_r_rasdef_ptr,
                char *pdbsum_exe_dir,struct c_file_data file_name_ptr)
{
  char command[LINELEN + 1], directory[FILENAME_LEN];
  char chain, log_file[FILENAME_LEN], prog_params[FILENAME_LEN];

  int len;
  int ok, dos;

  struct r_rasdef *r_rasdef_ptr;

  /* Initialise variables */
  ok = TRUE;
  dos = file_name_ptr.dos;

  /* Show progress */
  printf("   %s%d. In run_wiring ...\n",fno,process);

  /* Define directory */
  sprintf(directory,"%s/wiring",pdir);
  create_directory(directory);

  /* Form name of the log file */
  sprintf(log_file,"%s/wiring.log",log_dir);

  /* Get first rasdef record */
  r_rasdef_ptr = first_r_rasdef_ptr;

  /* Loop through the rasdefs until done */
  while (r_rasdef_ptr != NULL)
    {
      /* If this is a protein chain, then process */
      if (r_rasdef_ptr->type == R_PROTEIN && r_rasdef_ptr->copy == 0)
        {
          /* Get the chain id */
          chain = r_rasdef_ptr->chain;
          
          /* Run the program to generate the wiring diagram */
          sprintf(command,"cd %s && %s/wiring %s %c -pdbsum1 > %s",
                  directory,pdbsum_exe_dir,pdb_code,chain,log_file);
          call_system(command,dos);

          /* Repeat for the enlarged version */
          sprintf(command,"cd %s && %s/wiring %s %c -pdbsum1 -large > %s",
                  directory,pdbsum_exe_dir,pdb_code,chain,log_file);
          call_system(command,dos);

          /* Form the parameters for the mini domain wiring diagram */
          sprintf(prog_params,"%s %c -pfam -no_motifs -nresid 500 -u 0 -64 "
                  "-pdbsum1",pdb_code,chain);

          /* Generate the enlarged version of the mini domain wiring
             diagram */
          sprintf(command,"cd %s && %s/wiring %s > %s",
                  directory,pdbsum_exe_dir,prog_params,log_file);
          call_system(command,dos);
        }

      /* Get pointer to the next rasdef record */
      r_rasdef_ptr = r_rasdef_ptr->next_r_rasdef_ptr;
    }  
}
/***********************************************************************

run_pymol  -  Run PyMol to generate all the image files defined in the
              pymol.pml script file

***********************************************************************/

void run_pymol(char *directory,char *pml_file,char *pymol_exe,
               char *log_dir,struct c_file_data file_name_ptr)
{
  char command[FILENAME_LEN], log_file[FILENAME_LEN];
  char dos_exit[20], dos_start[20];

  int dos;

  /* Initialise variables */
  dos = file_name_ptr.dos;
  dos_exit[0] = dos_start[0] = '\0';

  /* Show progress */
  printf("          - Running PyMOL on clefts_img.pml\n");

  /* Form name of log file */
  sprintf(log_file,"%s/pymol_clefts.log",log_dir);
  
  /* For DOS, need to use the start command to prevent PyMOL shutting
     current window */
  if (dos == TRUE)
    {
      strcpy(dos_start,"start \"\" ");
      strcpy(dos_exit," && exit");
    }
  
  /* Form the parameters to run PyMol */
  sprintf(command,"cd %s && %s %s -c %s > %s %s",directory,dos_start,
          pymol_exe,pml_file,log_file,dos_exit);

  /* Run PyMol */
  call_system(command,dos);
}
/***********************************************************************

gen_clefts_image  -  Run PyMOL to enerate the image for the clefts page

***********************************************************************/

void gen_clefts_image(char *directory,char *pymol_exe,char *log_dir,
                      struct c_file_data file_name_ptr)
{
  char infile[FILENAME_LEN], input_line[LINELEN + 1], tmp[LINELEN + 1];
  char png_name[FILENAME_LEN], outfile[FILENAME_LEN];

  FILE *file_ptr, *outfile_ptr;

  /* Form the name of the output script file */
  sprintf(outfile,"%s/clefts_img.pml",directory);

  /* Open the file */
  outfile_ptr = open_output_file(outfile,TRUE);

  /* Form the name of the input script file */
  sprintf(infile,"%s/clefts.pml",directory);

  /* Form the name of the .png image file */
  sprintf(png_name,"%s/clefts.png",directory);

  /* Open the input file */
  if ((file_ptr = fopen(infile,"r")) !=  NULL)
    {
      /* Loop while reading in records from the file */
      while (fgets(input_line,LINELEN,file_ptr) != NULL)
        {
          /* Truncate the line at the last non-blank character */
          string_chomp(input_line);

          /* If this is the load line, replace by relative file name */
          if (!strncmp(input_line,"load ",5))
            strcpy(input_line,"load clefts.pdb");

          /* If vertices making a gap region, replace by surface */
          else if (!strncmp(input_line,"show line, ",11))
            {
              /* Make the change */
              strcpy(tmp,"show surface, ");
              strcat(tmp,input_line + 11);
              strcpy(input_line,tmp);
            }
          
          /* Write the line out */
          fprintf(outfile_ptr,"%s\n",input_line);
        }

      /* Close the file */
      fclose(file_ptr);

      /* Append the rendering commands */
      fprintf(outfile_ptr,"\n");
      fprintf(outfile_ptr,"# Render %s\n",png_name);

      /* Write out ray-trace command */
      fprintf(outfile_ptr,"ray 240, 240\n");
      fprintf(outfile_ptr,"viewport 240, 240\n");

      /* Write out command to generate png file */
      fprintf(outfile_ptr,"png %s\n",png_name);
    }

  /* Show error message */
  else
    {
      printf("*** Warning. Unable to open input file [%s]\n",
             infile);
    }

  /* Close the output file */
  fclose(outfile_ptr);

  /* Run PyMOL to generate the image */
  run_pymol(directory,outfile,pymol_exe,log_dir,file_name_ptr);
}
/***********************************************************************

run_speedfill  -  Run speedfill to generate the clefts

***********************************************************************/

void run_speedfill(char *fno,int process,char *pdir,char *pdb_code,
                   char *log_dir,int nchains,char *relinks,
                   struct c_file_data file_name_ptr)
{
  char bv[10], command[LINELEN + 1], directory[FILENAME_LEN];
  char log_file[FILENAME_LEN], prog_params[FILENAME_LEN];
  char pml_file[FILENAME_LEN], bv_file[FILENAME_LEN];

  int len, loop;
  int ok, dos;

  float resol;

  /* Initialise variables */
  strcpy(bv,"-bv");
  ok = TRUE;
  dos = file_name_ptr.dos;

  /* Show progress */
  printf("   %s%d. In run_speedfill ...\n",fno,process);

  /* Define directory */
  sprintf(directory,"%s/clefts",pdir);
  create_directory(directory);

  /* Form name of the log file */
  sprintf(log_file,"%s/speedfill.log",log_dir);

  /* Loop twice to calculate the surfaces and the clefts */
  for (loop = 0; loop < 2; loop++)
    {
      /* Form the input parameters */
      sprintf(prog_params,"-f ../%s.new -col 0 -ntop 10 -r3d -log -flip"
              " %s -pdbsum1 %s -pymol",pdb_code,bv,relinks);

      /* Run the program */
      sprintf(command,"cd %s && %s/speedfill %s >> %s",directory,
              file_name_ptr.other_dir[PDBSUM_EXE_DIR],prog_params,log_file);
      call_system(command,dos);

      /* Reset the binding site vertices */
      bv[0] = '\0';

      /* On the first loop, rename the clefts files */
      if (loop == 0)
        {
          /* Rename the clefts.pml file to clefts_bv.pml */
          sprintf(pml_file,"%s/clefts.pml",directory);
          sprintf(bv_file,"%s/clefts_bv.pml",directory);
          sprintf(command,"%s %s %s",CP_COMMAND,pml_file,bv_file);
          call_system(command,dos);
      
          /* Repeat for the clefts.rasmol file */
          sprintf(pml_file,"%s/clefts.rasmol",directory);
          sprintf(bv_file,"%s/clefts_bv.rasmol",directory);
          sprintf(command,"%s %s %s",CP_COMMAND,pml_file,bv_file);
          call_system(command,dos);
        }
    }
  
  /* Run PyMOL to enerate the image for the clefts page */
  gen_clefts_image(directory,file_name_ptr.exe_name[PYMOL_EXE],log_dir,
                   file_name_ptr);
}
/***********************************************************************

create_ligmet_entry  -  Create a ligmet entry

***********************************************************************/

struct ligmet *create_ligmet_entry(struct ligmet **fst_ligmet_ptr,
                                   struct ligmet **lst_ligmet_ptr,
                                   char *id)
{
  int i;

  struct ligmet *ligmet_ptr, *first_ligmet_ptr, *last_ligmet_ptr;

  /* Initialise pointers */
  ligmet_ptr = NULL;
  first_ligmet_ptr = *fst_ligmet_ptr;
  last_ligmet_ptr = *lst_ligmet_ptr;

  /* Create entry */
  ligmet_ptr = (struct ligmet *) malloc(sizeof(struct ligmet));
  if (ligmet_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for struct ligmet\n");
      printf("***        Program terminated with error.\n");
      exit (1);
    }

  /* Store the known info */
  store_string(&ligmet_ptr->id,id);

  /* Initialise remaining entries */
  ligmet_ptr->ligand_seq = &null;
  ligmet_ptr->het_name = &null;
  ligmet_ptr->resname1 = &null;
  ligmet_ptr->res1 = &null;
  ligmet_ptr->resname2 = &null;
  ligmet_ptr->res2 = &null;
  ligmet_ptr->chain = &null;
  ligmet_ptr->colour = &null;
  ligmet_ptr->from_dic = FALSE;
  ligmet_ptr->natoms = 0;
  ligmet_ptr->first_resid_ptr = NULL;

  /* Initialise pointer to next entry */
  ligmet_ptr->next_ligmet_ptr = NULL;

  /* If this is the first entry stored, then update
     pointer to head of linked list */
  if (first_ligmet_ptr == NULL)
    first_ligmet_ptr = ligmet_ptr;

  /* Otherwise, make previous entry point to the current
     one */
  else
    last_ligmet_ptr->next_ligmet_ptr = ligmet_ptr;
                      
  /* Save the current residue's pointer */
  last_ligmet_ptr = ligmet_ptr;

  /* Return current pointers */
  *fst_ligmet_ptr = first_ligmet_ptr;
  *lst_ligmet_ptr = last_ligmet_ptr;
  return(ligmet_ptr);
}
/***********************************************************************

create_resid_entry  -  Create a resid entry

***********************************************************************/

struct resid *create_resid_entry(struct resid **fst_resid_ptr,
                                 struct resid **lst_resid_ptr,
                                 char *id)
{
  int i;

  struct resid *resid_ptr, *first_resid_ptr, *last_resid_ptr;

  /* Initialise pointers */
  resid_ptr = NULL;
  first_resid_ptr = *fst_resid_ptr;
  last_resid_ptr = *lst_resid_ptr;

  /* Create entry */
  resid_ptr = (struct resid *) malloc(sizeof(struct resid));
  if (resid_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for struct resid\n");
      printf("***        Program terminated with error.\n");
      exit (1);
    }

  /* Store the known info */
  store_string(&resid_ptr->id,id);

  /* Initialise pointer to next entry */
  resid_ptr->next_resid_ptr = NULL;

  /* If this is the first entry stored, then update
     pointer to head of linked list */
  if (first_resid_ptr == NULL)
    first_resid_ptr = resid_ptr;

  /* Otherwise, make previous entry point to the current
     one */
  else
    last_resid_ptr->next_resid_ptr = resid_ptr;
                      
  /* Save the current residue's pointer */
  last_resid_ptr = resid_ptr;

  /* Return current pointers */
  *fst_resid_ptr = first_resid_ptr;
  *lst_resid_ptr = last_resid_ptr;
  return(resid_ptr);
}
/***********************************************************************

read_ligmet_run  -  Read in the ligands and metals to be plotted

***********************************************************************/

int read_ligmet_run(char *pdir,struct ligmet **fst_ligmet_ptr)
{
  char file_name[FILENAME_LEN], input_line[LINELEN + 1];
  char id[LINELEN], number_string[20], res_id[LINELEN];

  char *string_ptr;
  
  int len, nligmet;

  struct resid *resid_ptr, *first_resid_ptr, *last_resid_ptr;
  struct ligmet *ligmet_ptr, *first_ligmet_ptr, *last_ligmet_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  nligmet = 0;

  /* Initialise pointers */
  *fst_ligmet_ptr = first_ligmet_ptr = last_ligmet_ptr = NULL;
  last_resid_ptr = NULL;

  /* Form the name of the ligmet.run file */
  sprintf(file_name,"%s/newsec/ligmet.run",pdir);

  /* Open the file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    return(nligmet);

  /* Loop through the file, line by line */
  while (fgets(input_line,LINELEN,file_ptr) != NULL)
    {
      /* Truncate the input line */
      len = string_chomp(input_line);

      /* If the is a new entry, create a new record */
      if (!strncmp(input_line,"Name: ",6))
        {
          /* Retrieve the id */
          strcpy(id,input_line + 6);

          /* Create a new ligmet entry */
          ligmet_ptr
            = create_ligmet_entry(&first_ligmet_ptr,&last_ligmet_ptr,
                                  id);
          /* Increment count of etries read */
          nligmet++;

          /* Re-initialise pointer */
          last_resid_ptr = NULL;
        }

      /* Store given field */
      else if (!strncmp(input_line,"#residue ",9))
        {
          /* Extract the residue id */
          strcpy(res_id,input_line + 9);

          /* Create a resid record */
          resid_ptr
            = create_resid_entry(&(ligmet_ptr->first_resid_ptr),
                                 &last_resid_ptr,res_id);
        }

      /* Store given field */
      else if (!strncmp(input_line,"set ligand_seq = ",17))
        store_string(&(ligmet_ptr->ligand_seq),input_line + 17);

      /* Store given field */
      else if (!strncmp(input_line,"set het_name = ",15))
        store_string(&(ligmet_ptr->het_name),input_line + 15);

      /* Store given field */
      else if (!strncmp(input_line,"set resname1 = ",15))
        store_string(&(ligmet_ptr->resname1),input_line + 15);

      /* Store given field */
      else if (!strncmp(input_line,"set res1 = ",11))
        store_string(&(ligmet_ptr->res1),input_line + 11);

      /* Store given field */
      else if (!strncmp(input_line,"set resname2 = ",19))
        store_string(&(ligmet_ptr->resname2),input_line + 19);

      /* Store given field */
      else if (!strncmp(input_line,"set res2 = ",11))
        store_string(&(ligmet_ptr->res2),input_line + 11);

      /* Store given field */
      else if (!strncmp(input_line,"set chain = ",12))
        store_string(&(ligmet_ptr->chain),input_line + 12);

      /* Store given field */
      else if (!strncmp(input_line,"set colour = ",13))
        store_string(&(ligmet_ptr->colour),input_line + 13);

      /* Store given field */
      else if (!strcmp(input_line,"set from_dic = TRUE"))
        ligmet_ptr->from_dic = TRUE;

      /* Store given field */
      else if (!strncmp(input_line,"@ natoms = ",11))
        {
          strcpy(number_string,input_line + 11);
          ligmet_ptr->natoms = atoi(number_string);
        }
    }

  /* Close file */
  fclose(file_ptr);

  /* Return first ligmet record */
  *fst_ligmet_ptr = first_ligmet_ptr;

  /* Return the number of entries stored */
  return(nligmet);
}
/***********************************************************************

write_ligint_txt  -  Output the listing file giving the ligand-protein
                     or metal-protein interactions

***********************************************************************/

void write_ligint_txt(char *pdir,char *directory,char *pdb_code,
                      char *ligand_seq,int num1,int num2,int metal,
                      struct resid *first_resid_ptr)
{
  char file_name[FILENAME_LEN], out_file[FILENAME_LEN];
  char input_line[LINELEN + 1];
  char inchain1, inchain2;
  char atom_name1[5], atom_name2[5], atom_no1[6], atom_no2[6];
  char res_name1[4], res_name2[4], res_num1[6], res_num2[6];
  char count[20], dist[6], last_type, type;

  char *string_ptr;

  int i, ipos, len;
  int ncontacts, ndisulphides, nhbonds, nsalt_bridges;
  int swap, wanted;

  struct resid *resid_ptr;

  FILE *file_ptr, *outfile_ptr;

  /* Initialise variables */
  nhbonds = 0;
  ncontacts = 0;
  nsalt_bridges = 0;
  ndisulphides = 0;
  type = last_type = '\0';

  /* Form the output file name */
  if (metal == TRUE)
    sprintf(out_file,"%s/metint_%2.2d_%2.2d.txt",directory,num1,num2);
  else
    sprintf(out_file,"%s/ligint_%2.2d_%2.2d.txt",directory,num1,num2);

  /* Open the output file */
  outfile_ptr = open_output_file(out_file,TRUE);

  /* Write out the header records */
  if (metal == TRUE)
    {
      fprintf(outfile_ptr,"List of protein-metal interactions\n");
      fprintf(outfile_ptr,"----------------------------------\n");
    }
  else
    {
      fprintf(outfile_ptr,"List of protein-ligand interactions\n");
      fprintf(outfile_ptr,"-----------------------------------\n");
    }
  fprintf(outfile_ptr,"\n");
  fprintf(outfile_ptr,"                 PDB code: %s   ",pdb_code);
  if (metal == TRUE)
    fprintf(outfile_ptr,"Metal  ");
  else
    fprintf(outfile_ptr,"Ligand ");

  /* Remove quotes from the ligand sequence and get lenth */
  remove_char(ligand_seq,'"');
  len = strlen(ligand_seq);

  /* Write out the ligand sequence */
  fprintf(outfile_ptr,"%s\n",ligand_seq);

  /* Write out the dashes */
  fprintf(outfile_ptr,"                 ------------------------");
  for (i = 0; i < len; i++)
    fprintf(outfile_ptr,"-");
  fprintf(outfile_ptr,"\n");

  /* Form name of the grow.out file */
  sprintf(file_name,"%s/newsec/grow.out",pdir);
  
  /* Open the input file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    {
      printf("*** Warning. Unable to open file [%s]\n",file_name);
      return;
    }

  /* Loop through the file, line by line */
  while (fgets(input_line,LINELEN,file_ptr) != NULL)
    {
      /* Check line-length */
      len = string_chomp(input_line);
      swap = wanted = FALSE;

      /* Check that line is long enough */
      if (len > 60)
        {
          /* Get first residue record */
          resid_ptr = first_resid_ptr;

          /* Loop through all the records */
          while (resid_ptr != NULL && wanted == FALSE)
            {
              /* See if this residue is on the line */
              string_ptr = strstr(input_line,resid_ptr->id);

              /* If we have this residue in the line, then line
                 is wanted */
              if (string_ptr != NULL)
                {
                  /* Get ligand position */
                  ipos = string_ptr - input_line;

                  /* If it is in the first column, then want to swap */
                  if (ipos < 20)
                    swap = TRUE;

                  /* Set flag */
                  wanted = TRUE;
                }

              /* Get pointer to the next record */
              resid_ptr = resid_ptr->next_resid_ptr;
            }  
        }

      /* If line wanted, extract all the details */
      if (wanted == TRUE)
        {
          /* Extract the data from the line */
          type = input_line[0];
          inchain1 = input_line[6];
          strncpy(res_num1,input_line + 9, 5);
          res_num1[5] = '\0';
          strncpy(res_name1,input_line + 14, 3);
          res_name1[3] = '\0';
          strncpy(atom_no1,input_line + 18, 5);
          atom_no1[5] = '\0';
          strncpy(atom_name1,input_line + 24, 4);
          atom_name1[4] = '\0';
          strncpy(atom_name2,input_line + 33, 4);
          atom_name2[4] = '\0';
          strncpy(atom_no2,input_line + 38, 5);
          atom_no2[5] = '\0';
          strncpy(res_name2,input_line + 44, 3);
          res_name2[3] = '\0';
          strncpy(res_num2,input_line + 49, 5);
          res_num2[5] = '\0';
          inchain2 = input_line[54];
          strncpy(dist,input_line + 56, 5);
          dist[5] = '\0';

          /* If type has changed, write heading */
          if (type != last_type)
            {
              /* Salt bridges */
              if (type == 'S')
                {
                  fprintf(outfile_ptr,"\n");
                  fprintf(outfile_ptr,"Salt bridges\n");
                  fprintf(outfile_ptr,"------------\n");
                  fprintf(outfile_ptr,"\n");
                }

              /* Disulphide bonds */
              else if (type == 'D')
                {
                  fprintf(outfile_ptr,"\n");
                  fprintf(outfile_ptr,"Disulphide bonds\n");
                  fprintf(outfile_ptr,"----------------\n");
                  fprintf(outfile_ptr,"\n");
                }

              /* H-bonds */
              else if (type == 'H')
                {
                  fprintf(outfile_ptr,"\n");
                  fprintf(outfile_ptr,"Hydrogen bonds\n");
                  fprintf(outfile_ptr,"--------------\n");
                  fprintf(outfile_ptr,"\n");
                }

              /* Non-bonded contacts */
              else
                {
                  fprintf(outfile_ptr,"\n");
                  fprintf(outfile_ptr,"Non-bonded contacts\n");
                  fprintf(outfile_ptr,"-------------------\n");
                  fprintf(outfile_ptr,"\n");
                }

              /* Column headings */
              fprintf(outfile_ptr,"       <----- A T O M   1 ----->"
                      "<----- A T O M   2 ----->\n");
              fprintf(outfile_ptr,"\n");
              fprintf(outfile_ptr,"       Atom Atom Res  Res            "
                      "  Atom Atom Res  Res\n");
              fprintf(outfile_ptr,"        no. name name no.  Chain     "
                      "   no. name name no.  Chain  Distance\n");

              /* Save interaction type */
              last_type = type;
            }

          /* Increment count of this interaction type */
          if (type == 'S')
            {
              nsalt_bridges++;
              sprintf(count,"%3d.",nsalt_bridges);
            }
          else if (type == 'D')
            {
              ndisulphides++;
              sprintf(count,"%3d.",ndisulphides);
            }
          else if (type == 'H')
            {
              nhbonds++;
              sprintf(count,"%3d.",nhbonds);
            }
          else
            {
              ncontacts++;
              sprintf(count,"%3d.",ncontacts);
            }

          /* Print the current line */
          if (swap == FALSE)
            {
              fprintf(outfile_ptr,"%s  %s %s %s %s   %c   <--> ",
                      count,atom_no1,atom_name1,res_name1,res_num1,
                      inchain1);
              fprintf(outfile_ptr,"%s %s %s %s   %c     %s\n",
                      atom_no2,atom_name2,res_name2,res_num2,
                      inchain2,dist);
            }
          else
            {
              fprintf(outfile_ptr,"%s  %s %s %s %s   %c   <--> ",
                      count,atom_no2,atom_name2,res_name2,res_num2,
                      inchain2);
              fprintf(outfile_ptr,"%s %s %s %s   %c     %s\n",
                      atom_no1,atom_name1,res_name1,res_num1,
                      inchain1,dist);
            }
        }
      }

  /* Write counts of interactions */
  fprintf(outfile_ptr,"\n");
  if (nsalt_bridges > 0)
    fprintf(outfile_ptr,"Number of salt bridges:          %3d\n",
            nsalt_bridges);
  if (ndisulphides > 0)
    fprintf(outfile_ptr,"Number of disulphide bonds:      %3d\n",
            ndisulphides);
  if (nhbonds > 0)
    fprintf(outfile_ptr,"Number of hydrogen bonds:        %3d\n",
            nhbonds);
  if (ncontacts > 0)
    fprintf(outfile_ptr,"Number of non-bonded contacts:   %3d\n",
            ncontacts);

  /* Close both files */
  fclose(file_ptr);
  fclose(outfile_ptr);
}
/***********************************************************************

run_ligplots  -  Generate the ligplot diagrams for all the ligands and
                 metals

***********************************************************************/

void run_ligplots(char *fno,int process,char *pdir,char *pdb_code,
                  char *file_new,char *log_dir,char *pdbsum_exe_dir,
                  char *param_dir,char *ps2pdf,char *convert_exe,
                  struct ligmet *first_ligmet_ptr,int nligmet,
                  char *relinks,int run_romlas,
                  struct c_file_data file_name_ptr)
{
  char command[LINELEN + 1], directory[FILENAME_LEN];
  char log_file[FILENAME_LEN], prog_params[FILENAME_LEN];
  char file_name[FILENAME_LEN], id[20], number_string[20];
  char new_name[FILENAME_LEN], old_name[FILENAME_LEN];
  char out_file[FILENAME_LEN];

  char *string_ptr;

  int len, loop, iplot, nacross, nloops, num1, num2;
  int het_group, metal, ok, dos;

  struct ligmet *ligmet_ptr;

  /* Initialise variables */
  iplot = 0;
  het_group = metal = FALSE;
  nacross = 6;
  ok = TRUE;
  dos = file_name_ptr.dos;

  /* Show progress */
  printf("   %s%d. In run_ligplots ...\n",fno,process);

  /* Define directories */
  sprintf(directory,"%s/rasmol",pdir);
  create_directory(directory);
  sprintf(directory,"%s/ligplot",pdir);
  create_directory(directory);

  /* Copy across the required parameter file */
  sprintf(command,"%s %s/runhetchem.prm %s/ligplot.prm",CP_COMMAND,
          param_dir,directory);
  call_system(command,dos);

  /* Form name of the log file */
  sprintf(log_file,"%s/ligplot.log",log_dir);

  /* Show progress */
  printf("              ");
  fflush(stdout);

  /* Get first ligand record */
  ligmet_ptr = first_ligmet_ptr;

  /* Loop through the ligands until done */
  while (ligmet_ptr != NULL)
    {
      /* Initialise flags */
      het_group = metal = FALSE;

      /* Get this entry's id */
      strcpy(id,ligmet_ptr->id);
      string_ptr = strstr(id,".run");
      if (string_ptr != NULL)
        string_ptr[0] = '\0';

      /* Show progress */
      iplot++;
      if (iplot > nacross)
        {
          printf("\n");
          printf("              ");
          iplot = 1;
        }
      printf("%s ",id);
      fflush(stdout);

      /* If this is a metal or Het Group, then set flag */
      if (!strncmp(ligmet_ptr->id,"met",3))
        metal = TRUE;
      else if (!strncmp(ligmet_ptr->id,"het",3))
        het_group = TRUE;

      /* Determine how many loops are required */
      if (metal == TRUE)
        nloops = 1;
      else if (het_group == TRUE)
        nloops = 0;
      else
        nloops = 2;

      /* If not a Het Group, extract the two numbers, num1 and num2 */
      if (het_group != TRUE)
        {
          /* Extract the two numbers from the string */
          strcpy(number_string,id + 4);
          string_ptr = strchr(number_string,'_');
          if (string_ptr != NULL)
            {
              num2 = atoi(string_ptr + 1);
              string_ptr[0] = '\0';
            }
          num1 = atoi(number_string);

          /* If not the first copy of this ligand, then don't
             want the chemical ligplot */
          if (num2 > 1)
            nloops = 1;
        }

      /* Loop the required number of times */
      for (loop = 0; loop < nloops; loop++)
        {
          /* Get the parameter file for this loop */
          if (loop == 0)
            {
              /* Copy the standard parameter file */
              sprintf(command,"%s %s/runligplot.prm %s/ligplot.prm",
                      CP_COMMAND,param_dir,directory);
              call_system(command,dos);

              /* Get the file name */
              strcpy(file_name,file_new);
            }
          else
            {
              /* Copy the chemical ligplot parameter file */
              sprintf(command,"%s %s/runligchem.prm %s/ligplot.prm",
                      CP_COMMAND,param_dir,directory);
              call_system(command,dos);

              /* Change the id */
              sprintf(id,"clig_%2.2d",num1);

              /* Get the file name */
              sprintf(file_name,"newsec/lig_%2.2d.pdb",num1);
            }

          /* Form the command line parameters for ligplot */
          sprintf(prog_params,"../%s %s %s %s %s %s -r -pdbsum1",file_name,
                  ligmet_ptr->resname1,ligmet_ptr->res1,
                  ligmet_ptr->resname2,ligmet_ptr->res2,ligmet_ptr->chain);

          /* Add additional parameter, depending on whether metal or
             ligand */
          if (metal == TRUE)
            strcat(prog_params," -r -wat -m");
          else
            strcat(prog_params," -r -t");
          if (loop > 0)
            strcat(prog_params," -l -t");

          /* Run the program to generate the ligplot diagram */
          sprintf(command,"cd %s && %s/ligplot %s >> %s",directory,
                  pdbsum_exe_dir,prog_params,log_file);
          call_system(command,dos);

          /* Rename the PostScript file */
          sprintf(old_name,"ligplot.ps");
          sprintf(new_name,"%s.ps",id);
          move_file(directory,old_name,new_name,dos);
  //          sprintf(command,"cd %s && %s ligplot.ps %s.ps",directory,
  //                  MV_COMMAND,id);
  //          call_system(command,dos);

          /* Repeat for the .drw file */
          sprintf(old_name,"ligplot.drw");
          sprintf(new_name,"%s.drw",id);
          move_file(directory,old_name,new_name,dos);
  //          sprintf(command,"cd %s && %s ligplot.drw %s.drw",directory,
  //                  MV_COMMAND,id);
  //          call_system(command,dos);

          /* Repeat for the .pdb file */
          sprintf(old_name,"ligplot.pdb");
          sprintf(new_name,"%s.pdb",id);
          move_file(directory,old_name,new_name,dos);
  //          sprintf(command,"cd %s && %s ligplot.pdb %s.pdb",directory,
  //                  MV_COMMAND,id);
  //          call_system(command,dos);

          /* Repeat for the .res file */
          sprintf(old_name,"ligplot.res");
          sprintf(new_name,"%s.res",id);
          move_file(directory,old_name,new_name,dos);
  //          sprintf(command,"cd %s && %s ligplot.res %s.res",directory,
  //                  MV_COMMAND,id);
  //          call_system(command,dos);

          /* Convert the PostScript file to PDF */
          if (strcmp(ps2pdf,"NONE"))
            {
              /* Convert the PostScript file */
              sprintf(command,"cd %s && %s %s.ps %s.pdf",directory,ps2pdf,
                      id,id);
              call_system(command,dos);
            }

          /* Form the command line parameters to convert to gif */
          sprintf(prog_params," -density 200x200 -trim +repage -transparent "
                  "white ");

          /* Convert the PostScript file to GIF */
/* v.1.0.1 --> */
//          sprintf(command,"cd %s && %s %s %s.ps %s.gif",directory,
//                  convert_exe,prog_params,id,id);
          sprintf(command,"cd %s && %s %s.ps %s %s.gif",directory,
                  convert_exe,id,prog_params,id);
/* <-- v.1.0.1 */
          call_system(command,dos);

          /* If this is the first loop, write out the text file of
             interactions and the RasMol script */
          if (loop == 0)
            {
              /*  Form the name of the output script file */
              sprintf(out_file,"../rasmol/%s.rasmol",id);

              /* Define the parameters */
              sprintf(prog_params,"%s%s -lr -ll -lp -pdbsum1 %s -pdb %s.pdb",
                      pdb_code,id + 4,relinks,id);

              /* Run the program */
              if (run_romlas == TRUE)
                {
                  sprintf(command,"cd %s && %s/romlas %s > %s",directory,
                          pdbsum_exe_dir,prog_params,out_file);
                  call_system(command,dos);
                }
              
              /* Write out the text file of interactions */
              write_ligint_txt(pdir,directory,pdb_code,
                               ligmet_ptr->ligand_seq,num1,num2,metal,
                               ligmet_ptr->first_resid_ptr);

              /* Delete the corresponding .pdb and .res files */
              sprintf(command,"%s %s/%s.pdb",RM_COMMAND,directory,id);
              sprintf(command,"%s %s/%s.res",RM_COMMAND,directory,id);
            }
        }

      /* Get pointer to the next ligand record */
      ligmet_ptr = ligmet_ptr->next_ligmet_ptr;
    }

  /* Close off progress line */
  printf("\n");
  fflush(stdout);
}
/***********************************************************************

run_nucplot  -  Generate the nucplot diagrams for any protein-DNA
                interactions

***********************************************************************/

void run_nucplot(char *fno,int process,char *pdir,char *pdb_code,
                 char *file_new,char *log_dir,char *pdbsum_exe_dir,
                 char *param_dir,char *ps2pdf,
                 struct c_file_data file_name_ptr)
{
  char command[LINELEN + 1], directory[FILENAME_LEN];
  char log_file[FILENAME_LEN], prog_params[FILENAME_LEN];

  int ok, dos;

  /* Initialise variables */
  ok = TRUE;
  dos = file_name_ptr.dos;

  /* Show progress */
  printf("   %s%d. In run_nucplot ...\n",fno,process);

  /* Define directory */
  sprintf(directory,"%s/nucplot",pdir);
  create_directory(directory);

  /* Copy across the required parameter file */
  sprintf(command,"%s %s/runnucplot.par %s/nucplot.par",CP_COMMAND,
          param_dir,directory);
  call_system(command,dos);

  /* Form name of the log file */
  sprintf(log_file,"%s/nucplot.log",log_dir);

  /* Form the command line parameters for nucplot */
  sprintf(prog_params,"../%s -pdbsum1",file_new);

  /* Run the program to generate the nucplot diagram */
  sprintf(command,"cd %s && %s/nucplot %s >> %s",directory,
          pdbsum_exe_dir,prog_params,log_file);
  call_system(command,dos);

  /* Convert the PostScript file to PDF */
  if (strcmp(ps2pdf,"NONE"))
    {
      /* Convert the PostScript file */
      sprintf(command,"cd %s && %s nucplot.ps nucplot.pdf",directory,ps2pdf);
      call_system(command,dos);
    }
}
/***********************************************************************

run_seqdata  -  Run seqdata to extract and UniProt and E.C. data

***********************************************************************/

void run_seqdata(char *fno,int process,char *pdir,char *pdb_code,
                 char *file_new,char *log_dir,char *pdbsum_exe_dir,
                 struct c_file_data file_name_ptr)
{
  char command[LINELEN + 1], directory[FILENAME_LEN];
  char log_file[FILENAME_LEN], prog_params[FILENAME_LEN];

  int ok, dos;

  /* Initialise variables */
  ok = TRUE;
  dos = file_name_ptr.dos;

  /* Show progress */
  printf("   %s%d. In run_seqdata ...\n",fno,process);

  /* Define directory */
  sprintf(directory,"%s/seqdata",pdir);
  create_directory(directory);

  /* Form name of the log file */
  sprintf(log_file,"%s/seqdata.log",log_dir);

  /* Form the command line parameters for seqdata */
  sprintf(prog_params,"-pdbsum1 %s",pdb_code);

  /* Run the program to generate the seqdata diagram */
  sprintf(command,"cd %s && %s/seqdata %s >> %s",directory,
          pdbsum_exe_dir,prog_params,log_file);
  call_system(command,dos);
}
/***********************************************************************

run_uniplot  -  Run uniplot to extract and UniProt and E.C. data

***********************************************************************/

void run_uniplot(char *fno,int process,char *pdir,char *pdb_code,
                 char *file_new,char *log_dir,char *pdbsum_exe_dir,
                 struct c_file_data file_name_ptr)
{
  char command[LINELEN + 1], directory[FILENAME_LEN];
  char log_file[FILENAME_LEN], prog_params[FILENAME_LEN];

  int ok, dos;

  /* Initialise variables */
  ok = TRUE;
  dos = file_name_ptr.dos;

  /* Show progress */
  printf("   %s%d. In run_uniplot ...\n",fno,process);

  /* Define directory */
  sprintf(directory,"%s/seqdata",pdir);

  /* Form name of the log file */
  sprintf(log_file,"%s/uniplot.log",log_dir);

  /* Form the command line parameters for uniplot */
  sprintf(prog_params,"-pdbsum1 %s -wget",pdb_code);

  /* Run the program to generate the uniplot diagram */
  sprintf(command,"cd %s && %s/uniplot %s >> %s",directory,
          pdbsum_exe_dir,prog_params,log_file);
  call_system(command,dos);
}
/***********************************************************************

run_pdbhtml  -  Run pdbhtml to generate the database.dat file

***********************************************************************/

void run_pdbhtml(char *fno,int process,char *pdir,char *pdb_code,
                 char *log_dir,int run_number,
                 struct c_file_data file_name_ptr,
                 struct parameters params)
{
  char command[LINELEN + 1], directory[FILENAME_LEN];
  char log_file[FILENAME_LEN], param_dir[FILENAME_LEN];
  char in_file[FILENAME_LEN], string[LINELEN + 1];
  char file_name[FILENAME_LEN];

  int len;
  int have_file, ok, dos;

  float resol;

  FILE *file_ptr;

  /* Initialise variables */
  ok = TRUE;
  strcpy(param_dir,file_name_ptr.other_dir[PDBSUM_PARAMS_DIR]);
  dos = file_name_ptr.dos;

  /* Show progress */
  printf("   %s%d. In run_pdbhtml ...\n",fno,process);

  /* Create the PDBHTML directory */
  sprintf(directory,"%s/pdbhtml",pdir);
  create_directory(directory);

  /* Copy CATHPARAM file into this directory - no longer required */
  //sprintf(command,"%s %s/CATHPARAM %s",CP_COMMAND,param_dir,directory);
  //call_system(command,dos);

  /* Copy the data file across */
  sprintf(command,"%s %s/%s.new %s/file.pdb",CP_COMMAND,pdir,pdb_code,
          directory);
  call_system(command,dos);

  /* Form name of the log file */
  sprintf(log_file,"%s/pdbhtml.log",log_dir);

  /* Define input parameter file */
  sprintf(in_file,"%s/in",directory);

  /* Write parameters for pdbhtml */
  sprintf(command,"%s\n1\n",pdb_code);
  write_in_file(directory,command);

  //  sprintf(command,"echo %s > %s",pdb_code,in_file);
  //  call_system(command,dos);
  //  sprintf(command,"echo 1 >> %s",in_file);
  //  call_system(command,dos);

  /* Run the program */
  sprintf(command,"cd %s && %s/pdbhtml < %s > %s",directory,
          file_name_ptr.other_dir[PDBSUM_EXE_DIR],in_file,log_file);
  call_system(command,dos);

  /* Form the name of the database.dat file */
  sprintf(file_name,"%s/database.dat",directory);

  /* Open the file for append */
  file_ptr = fopen(file_name,"a");

  /* If OK, then append the run number */
  if (file_ptr != NULL)
    fprintf(file_ptr,"RUN_NUMBER %d\n",run_number);

  /* Close the file */
  fclose(file_ptr);

  /* Copy the database.dat file up */
  sprintf(command,"%s %s/database.dat %s",CP_COMMAND,directory,pdir);
  call_system(command,dos);
}
/***********************************************************************

promotif_rasmol  -  Generate the PROMOTIF RasMol scripts

***********************************************************************/

void promotif_rasmol(char *directory,char *pdb_code,
                     char *pdbsum_exe_dir,char *log_file,
                     struct plot **plot_index_ptr,int nchains,
                     char *relinks,struct c_file_data file_name_ptr)
{
  char chain, out_file[FILENAME_LEN];
  char command[FILENAME_LEN], prog_params[FILENAME_LEN];

  int ichain, imotif, iplot;
  int dos;

  struct plot *plot_ptr;

  FILE *logfile_ptr;

  /* Initialise variables */
  dos = file_name_ptr.dos;

  /* Show progress */
  printf("         - promotif scripts ...\n");

  /* Open the log file */
  logfile_ptr = open_output_file(log_file,TRUE);

  /* Loop over the chains to be plotted */
  for (ichain = 0; ichain < nchains; ichain++)
    {
      /* Get the plot record for this chain */
      plot_ptr = plot_index_ptr[ichain];

      /* Get the corresponding chain */
      chain = plot_ptr->chain;

      /* Loop over all the plots to generate a separate page for
         each */
      for (iplot = 0; iplot < NPLOTS; iplot++)
        {
          /* If generating a script for this motif type, generate one
             for each instance */
          if (MOTIF_NO[iplot] > -1)
            {
              /* Loop over all the motifs that require scripts */
              for (imotif = 0; imotif < plot_ptr->count[iplot]; imotif++)
                {
                  /* Form the name of the output script file */
                  sprintf(out_file,"%s/%s_%c_%d.rasmol",directory,
                          PLOT_DESC[iplot],plot_ptr->chain,imotif + 1);

                  /* Define the parameters */
                  sprintf(prog_params,
                          "%s%c -mc -chain -promotif %d %d -pdbsum1 %s",
                          pdb_code,chain,MOTIF_NO[iplot],imotif + 1,
                          relinks);

                  /* Run the program */
                  sprintf(command,"cd %s && %s/romlas %s > %s",directory,
                          pdbsum_exe_dir,prog_params,out_file);
                  call_system(command,dos);

                  /* Write the command to the log file */
                  fprintf(logfile_ptr,"%s\n",command);
                }
            }
        }
    }

  /* Close the log file */
  fclose(logfile_ptr);
}
/***********************************************************************

main_rasmol  -  Generate the main RasMol scripts

***********************************************************************/

void main_rasmol(char *directory,char *pdb_code,
                 char *pdbsum_exe_dir,struct chain *first_chain_ptr,
                 int nchains,char *relinks,
                 struct c_file_data file_name_ptr)
{
  char command[LINELEN + 1], out_file[FILENAME_LEN];
  char prog_params[FILENAME_LEN];
  char chain, chains[3];

  int len;
  int ok, dos;

  struct chain *chain_ptr;
  struct iface *iface_ptr;

  /* Initialise variables */
  dos = file_name_ptr.dos;

  /* Form the name of the output script file */
  sprintf(out_file,"%s/%s.rasmol",directory,pdb_code);

  /* Define the parameters for rascript */
  sprintf(prog_params,"%s -arad -1.0 -pdbsum1 %s",pdb_code,relinks);

  /* Run the rascript program */
  sprintf(command,"cd %s && %s/rascript %s > %s",directory,pdbsum_exe_dir,
          prog_params,out_file);
  call_system(command,dos);

  /* Get first record */
  chain_ptr = first_chain_ptr;

  /* Loop through all the records */
  while (chain_ptr != NULL)
    {
      /* Get the chain */
      chain = chain_ptr->chain;

      /* Generate the RasMol script for this chain */
      sprintf(out_file,"%s/chain_%c.rasmol",directory,chain);

      /* Form the rascript parameters */
      sprintf(prog_params,"%s -c%c -arad -1.0 -pdbsum1 %s",pdb_code,
              chain,relinks);

      /* Run the rascript program */
      sprintf(command,"cd %s && %s/rascript %s > %s",directory,
              pdbsum_exe_dir,prog_params,out_file);
      call_system(command,dos);

      /* Get first interface record */
      iface_ptr = chain_ptr->first_iface_ptr;

      /* Loop through all the records */
      while (iface_ptr != NULL)
        {
          /* Join the two chains */
          sprintf(chains,"%c%c",chain_ptr->chain,
                  iface_ptr->chain_ptr->chain);

          /* Generate the RasMol script for this chain pair */
          sprintf(out_file,"%s/chains_%s.rasmol",directory,chains);

          /* Form the rascript parameters */
          sprintf(prog_params,"%s -c%c -c%c -arad -1.0 -pdbsum1 %s",
                  pdb_code,chain,iface_ptr->chain_ptr->chain,relinks);

          /* Run the rascript program */
          sprintf(command,"cd %s && %s/rascript %s > %s",directory,
                  pdbsum_exe_dir,prog_params,out_file);
          call_system(command,dos);

          /* Get pointer to the next record */
          iface_ptr = iface_ptr->next_iface_ptr;
        }  

      /* Get pointer to the next record */
      chain_ptr = chain_ptr->next_chain_ptr;
    }  
}
/***********************************************************************

rasmol_scripts  -  Generate all the required RasMol scripts

***********************************************************************/

void rasmol_scripts(char *fno,int process,char *pdir,char *pdb_code,
                    char *file_name,char *file_new,char *log_dir,
                    struct r_rasdef *first_r_rasdef_ptr,
                    struct plot **plot_index_ptr,int nuniq,
                    struct chain *first_chain_ptr,int nchains,
                    char *relinks,struct c_file_data file_name_ptr)
{
  char command[LINELEN + 1], directory[FILENAME_LEN];
  char log_file[FILENAME_LEN], out_file[FILENAME_LEN];
  char prog_params[FILENAME_LEN], pdbsum_exe_dir[FILENAME_LEN];

  int len;
  int have_file, ok, dos;

  float resol;

  /* Initialise variables */
  ok = TRUE;
  strcpy(pdbsum_exe_dir,file_name_ptr.other_dir[PDBSUM_EXE_DIR]);
  dos = file_name_ptr.dos;

  /* Show progress */
  printf("   %s%d. In rasmol_scripts ...\n",fno,process);

  /* Create the rasmol directory */
  sprintf(directory,"%s/rasmol",pdir);
  create_directory(directory);

  /* Form name of the log file */
  sprintf(log_file,"%s/rasmol_scripts.log",log_dir);

  /* Generate the main rasmol scripts */
  main_rasmol(directory,pdb_code,pdbsum_exe_dir,first_chain_ptr,
              nchains,relinks,file_name_ptr);

  /* If have any Promotif plots, generate all their scripts */
  if (nuniq > 0)
    promotif_rasmol(directory,pdb_code,pdbsum_exe_dir,log_file,
                    plot_index_ptr,nuniq,relinks,file_name_ptr);
}
/***********************************************************************

main_pymol  -  Generate the main Pymol scripts

***********************************************************************/

void main_pymol(char *directory,char *pdb_code,
                char *pdbsum_exe_dir,struct chain *first_chain_ptr,
                int nchains,char *log_file,char *relinks,
                struct c_file_data file_name_ptr)
{
  char command[LINELEN + 1], out_file[FILENAME_LEN];
  char prog_params[FILENAME_LEN];
  char chain, chains[3];

  int len;
  int ok, dos;

  struct chain *chain_ptr;
  struct iface *iface_ptr;

  /* Initialise variables */
  dos = file_name_ptr.dos;

  /* Form the name of the output script file */
  sprintf(out_file,"%s/%s.pml",directory,pdb_code);

  /* Define the parameters for genpymol */
  sprintf(prog_params,"%s -pdbsum1 %s -pml %s",pdb_code,relinks,
          out_file);

  /* Run the genpymol program */
  sprintf(command,"cd %s && %s/genpymol %s >> %s",directory,
          pdbsum_exe_dir,prog_params,log_file);
  call_system(command,dos);

  /* Get first record */
  chain_ptr = first_chain_ptr;

  /* Loop through all the records */
  while (chain_ptr != NULL)
    {
      /* Get the chain */
      chain = chain_ptr->chain;

      /* Generate the Pymol script for this chain */
      sprintf(out_file,"%s/chain_%c.pml",directory,chain);

      /* Form the genpymol parameters */
      sprintf(prog_params,"%s %c -pdbsum1 %s -pml %s",pdb_code,chain,
              relinks,out_file);

      /* Run the genpymol program */
      sprintf(command,"cd %s && %s/genpymol %s >> %s",directory,
              pdbsum_exe_dir,prog_params,log_file);
      call_system(command,dos);

      /* Get first interface record */
      iface_ptr = chain_ptr->first_iface_ptr;

      /* Loop through all the records */
      while (iface_ptr != NULL)
        {
          /* Join the two chains */
          sprintf(chains,"%c%c",chain_ptr->chain,
                  iface_ptr->chain_ptr->chain);

          /* Generate the Pymol script for this chain pair */
          sprintf(out_file,"%s/chains_%s.pymol",directory,chains);

          /* Form the rascript parameters */
          sprintf(prog_params,"%s -c%c -c%c -arad -1.0 -pdbsum1 %s",
                  pdb_code,chain,iface_ptr->chain_ptr->chain,relinks);

          /* Run the rascript program */
          sprintf(command,"cd %s && %s/rascript %s > %s",directory,
                  pdbsum_exe_dir,prog_params,out_file);
          call_system(command,dos);

          /* Get pointer to the next record */
          iface_ptr = iface_ptr->next_iface_ptr;
        }  

      /* Get pointer to the next record */
      chain_ptr = chain_ptr->next_chain_ptr;
    }  
}
/***********************************************************************

pymol_scripts  -  Generate all the required Pymol scripts

***********************************************************************/

void pymol_scripts(char *fno,int process,char *pdir,char *pdb_code,
                   char *file_name,char *file_new,char *log_dir,
                   struct r_rasdef *first_r_rasdef_ptr,
                   struct plot **plot_index_ptr,int nuniq,
                   struct chain *first_chain_ptr,int nchains,
                   char *relinks,struct c_file_data file_name_ptr)
{
  char command[LINELEN + 1], directory[FILENAME_LEN];
  char log_file[FILENAME_LEN], out_file[FILENAME_LEN];
  char prog_params[FILENAME_LEN], pdbsum_exe_dir[FILENAME_LEN];

  int len;
  int have_file, ok;

  float resol;

  /* Initialise variables */
  ok = TRUE;
  strcpy(pdbsum_exe_dir,file_name_ptr.other_dir[PDBSUM_EXE_DIR]);

  /* Show progress */
  printf("   %s%d. In pymol_scripts ...\n",fno,process);

  /* Create the pymol directory */
  sprintf(directory,"%s/pymol",pdir);
  create_directory(directory);

  /* Form name of the log file */
  sprintf(log_file,"%s/pymol_scripts.log",log_dir);

  /* Generate the main pymol scripts */
  main_pymol(directory,pdb_code,pdbsum_exe_dir,first_chain_ptr,
             nchains,log_file,relinks,file_name_ptr);

  /* If have any Promotif plots, generate all their scripts */
  //  if (nuniq > 0)
  //    promotif_pymol(directory,pdb_code,pdbsum_exe_dir,log_file,
  //                    plot_index_ptr,nuniq);
}
/***********************************************************************

gen_protein_html  -  Generate the PROTEIN page

***********************************************************************/

void gen_protein_html(char *directory,char *pdb_code,
                      char *pdbsum_exe_dir,char *log_file,char *relinks,
                      struct r_rasdef *first_r_rasdef_ptr,
                      struct c_file_data file_name_ptr)
{
  char out_file[FILENAME_LEN], template[FILENAME_LEN];
  char command[FILENAME_LEN], prog_params[FILENAME_LEN];
  char chain;

  int level;
  int dos;

  struct r_rasdef *r_rasdef_ptr;

  /* Initialise variables */
  level = 0;
  dos = file_name_ptr.dos;

  /* Define the name of the template and output file */
  strcpy(template,"protein.html");

  /* Get first rasdef record */
  r_rasdef_ptr = first_r_rasdef_ptr;

  /* Loop through the rasdefs until done */
  while (r_rasdef_ptr != NULL)
    {
      /* If this is a protein chain, then process */
      if (r_rasdef_ptr->type == R_PROTEIN && r_rasdef_ptr->copy == 0)
        {
          /* Get the chain id */
          chain = r_rasdef_ptr->chain;

          /* Form the name of the HTML file */
          sprintf(out_file,"%s/protein_%c.html",directory,chain);

          /* Show progress */
          printf("         - %s ...\n",template);

          /* Define the parameters */
          level++;
          sprintf(prog_params,"%s %s -chain %c -l %d -c 0 %s -f %s "
                  "-pdbsum1",
                  pdb_code,template,chain,level,relinks,out_file);

          /* Run the program */
          sprintf(command,"cd %s && %s/showmain %s > %s",directory,
                  pdbsum_exe_dir,prog_params,log_file);
          call_system(command,dos);
        }

      /* Get pointer to the next rasdef record */
      r_rasdef_ptr = r_rasdef_ptr->next_r_rasdef_ptr;
    }  
}
/***********************************************************************

gen_promotif_html  -  Generate the PROMOTIF pages

***********************************************************************/

void gen_promotif_html(char *directory,char *pdb_code,
                       char *pdbsum_exe_dir,char *log_file,
                       struct plot **plot_index_ptr,int nuniq,
                       char *relinks,struct c_file_data file_name_ptr)
{
  char out_file[FILENAME_LEN], template[FILENAME_LEN];
  char command[FILENAME_LEN], prog_params[FILENAME_LEN];

  int call_no, ichain, iplot;
  int dos;

  struct plot *plot_ptr;

  /* Initialise variables */
  call_no = 0;
  dos = file_name_ptr.dos;

  /* Define the name of the template and output file */
  strcpy(template,"protein.html");

  /* Show progress */
  printf("         - promotif plots ...\n");

  /* Loop over the chains to be plotted */
  for (ichain = 0; ichain < nuniq; ichain++)
    {
      /* Get the plot record for this chain */
      plot_ptr = plot_index_ptr[ichain];

      /* Loop over all the plots to generate a separate page for
         each */
      for (iplot = 0; iplot < NPLOTS; iplot++)
        {
          /* If this plot has motifs, then generate page */
          if (plot_ptr->count[iplot] > 0)
            {
              /* Increment call number */
              call_no++;

              /* Form the names of the HTML file */
              sprintf(out_file,"%s/%s_%c.html",directory,
                      PLOT_DESC[iplot],plot_ptr->chain);

              /* Define the parameters */
              sprintf(prog_params,"%s %s -f %s -pdbsum1 -l %d -c %d "
                      "-option %s %s",
                      pdb_code,template,out_file,ichain + 1,
                      call_no,PLOT_NAME[iplot],relinks);

              /* Run the program */
              sprintf(command,"cd %s && %s/showmain %s > %s",directory,
                      pdbsum_exe_dir,prog_params,log_file);
              call_system(command,dos);
            }
        }
    }
}
/***********************************************************************

gen_procheck_html  -  Generate the PROCHECK page

***********************************************************************/

void gen_procheck_html(char *directory,char *pdb_code,
                       char *pdbsum_exe_dir,char *log_file,
                       char *relinks,struct c_file_data file_name_ptr)
{
  char out_file[FILENAME_LEN], template[FILENAME_LEN];
  char command[FILENAME_LEN], prog_params[FILENAME_LEN];

  int dos;

  /* Initialise variables */
  dos = file_name_ptr.dos;
  
  /* Define the name of the template and output file */
  strcpy(template,"procheck.html");
  sprintf(out_file,"%s/procheck.html",directory);

  /* Show progress */
  printf("         - %s ...\n",template);

  /* Define the parameters */
  sprintf(prog_params,"%s %s -f %s -pdbsum1 %s",pdb_code,template,
          out_file,relinks);

  /* Run the program */
  sprintf(command,"cd %s && %s/showmain %s > %s",directory,pdbsum_exe_dir,
          prog_params,log_file);
  call_system(command,dos);
}
/***********************************************************************

gen_dna_html  -  Generate the DNA/RNA page

***********************************************************************/

void gen_dna_html(char *directory,char *pdb_code,
                  char *pdbsum_exe_dir,char *log_file,char *relinks,
                  struct c_file_data file_name_ptr)
{
  char out_file[FILENAME_LEN], template[FILENAME_LEN];
  char command[FILENAME_LEN], prog_params[FILENAME_LEN];

  int dos;

  /* Initialise variables */
  dos = file_name_ptr.dos;

  /* Define the name of the template and output file */
  strcpy(template,"dna.html");
  sprintf(out_file,"%s/dna.html",directory);

  /* Show progress */
  printf("         - %s ...\n",template);

  /* Define the parameters */
  sprintf(prog_params,"%s %s -f %s -pdbsum1 %s",pdb_code,template,
          out_file,relinks);

  /* Run the program */
  sprintf(command,"cd %s && %s/showmain %s > %s",directory,pdbsum_exe_dir,
          prog_params,log_file);
  call_system(command,dos);
}
/***********************************************************************

insert_file  -  Insert the given PDB of .dfn file at this point

***********************************************************************/

int insert_file(FILE *outfile_ptr,char *file_name,int is_pdb)
{
  char input_line[LINELEN + 1];

  int len, nout;

  FILE *file_ptr;

  /* Initialise variables */
  nout = 0;

  /* Open the input file */
  file_ptr = fopen(file_name,"r");
  if (file_ptr == NULL)
    {
      printf("*** Warning. Unable to open file [%s] \"\n",file_name);
      return(nout);
    }

  /* Loop through the file, copying to output file */
  while (fgets(input_line,LINELEN,file_ptr) != NULL)
    {
      /* Truncate the string to strip off the newline */
      len = string_chomp(input_line);

      /* If a PDB file, on;y want ATOM, HETATM, and CONECT records */
      if (!strncmp(input_line,"ATOM  ",6) ||
          !strncmp(input_line,"HETATM",6) ||
          !strncmp(input_line,"CONECT",6) ||
          is_pdb == FALSE)
        {
          /* If line too long, then truncate */
          if (len > 66)
            input_line[66] = '\0';
      
          /* Write the line out */
          fprintf(outfile_ptr,"%c%s%c%c%c +\n",'"',input_line,'\\',
                  'n','"');

          /* Increment count of lines output */
          nout++;
        }
    }

  /* Close the input file */
  fclose(file_ptr);

  /* Return number of lines written out */
  return(nout);
}
/***********************************************************************

rewrite_js_file  -  Rewrite the given javascript file to add in the atom
                    coords required

***********************************************************************/

void rewrite_js_file(char *tmp_file,char *js_file,char *pdb_file,
                     char *dfn_file)
{
  char input_line[LINELEN + 1];

  int len;

  FILE *infile_ptr, *outfile_ptr;

  /* Initialise variables */

  /* Open the input javascipt file */
  infile_ptr = fopen(tmp_file,"r");
  if (infile_ptr ==  NULL)
    {
      printf("*** ERROR. Unable to open js file [%s] \"\n",tmp_file);
      return;
    }

  /* Open the output file */
  outfile_ptr = open_output_file(js_file,TRUE);

  /* Loop through the file, copying to output file */
  while (fgets(input_line,LINELEN,infile_ptr) != NULL)
    {
      /* Truncate the string to strip off the newline */
      string_chomp(input_line);

      /* Check for the insertion point for given file */
      if (strstr(input_line,"#insert PDB") != NULL)
        insert_file(outfile_ptr,pdb_file,TRUE);
      else if (strstr(input_line,"#insert dfn") != NULL)
        insert_file(outfile_ptr,dfn_file,FALSE);

      /* Otherwise, write the javascript line straight out */
      else
        fprintf(outfile_ptr,"%s\n",input_line);
    }

  /* Close both files */
  fclose(infile_ptr);
  fclose(outfile_ptr);
}
/***********************************************************************

gen_ligands_html  -  Generate the ligands and metals pages

***********************************************************************/

void gen_ligands_html(char *directory,char *pdb_code,
                      char *pdbsum_exe_dir,char *log_file,
                      struct ligmet *first_ligmet_ptr,int nligmets,
                      char *relinks,struct c_file_data file_name_ptr)
{
  char out_file[FILENAME_LEN], template[FILENAME_LEN], id[20];
  char ligmets[3], command[FILENAME_LEN], prog_params[FILENAME_LEN];
  char number_string[20], tmp_file[FILENAME_LEN];
  char pdb_file[FILENAME_LEN], dfn_file[FILENAME_LEN];
  char js_file[FILENAME_LEN];

  char *string_ptr;

  int iplot, nacross, num1, num2;
  int het_group, metal, dos;

  struct ligmet *ligmet_ptr;

  /* Initialise variables */
  dfn_file[0] = '\0';
  iplot = 0;
  nacross = 6;
  num1 = num2 = 0;
  dos = file_name_ptr.dos;

  /* Define the name of the ligand_plugin.js file */
  sprintf(js_file,"%s/js/ligand_plugin.js",
          file_name_ptr.other_dir[PDBSUM_REAL_TEMPLATES_DIR]);

  /* Define the name of the template */
  strcpy(template,"ligands.html");

  /* Show progress */
  printf("         - ligand pages ...\n");
  printf("              ");
  fflush(stdout);

  /* Get first record */
  ligmet_ptr = first_ligmet_ptr;

  /* Loop through all the records */
  while (ligmet_ptr != NULL)
    {
      /* Initialise flags */
      het_group = metal = FALSE;

      /* Get this entry's id */
      strcpy(id,ligmet_ptr->id);
      string_ptr = strstr(id,".run");
      if (string_ptr != NULL)
        string_ptr[0] = '\0';

      /* Show progress */
      iplot++;
      if (iplot > nacross)
        {
          printf("\n");
          printf("              ");
          iplot = 1;
        }
      printf("%s ",id);
      fflush(stdout);

      /* If this is a metal or Het Group, then set flag */
      if (!strncmp(ligmet_ptr->id,"met",3))
        metal = TRUE;
      else if (!strncmp(ligmet_ptr->id,"het",3))
        het_group = TRUE;

      /* Only generate page for ligand or metal */
      if (het_group == FALSE)
        {
          /* Extract the two numbers from the string */
          strcpy(number_string,id + 4);
          string_ptr = strchr(number_string,'_');
          if (string_ptr != NULL)
            {
              num2 = atoi(string_ptr + 1);
              string_ptr[0] = '\0';
            }
          num1 = atoi(number_string);

          /* Form the name of the PDB file with the ligand coordinates */
          sprintf(pdb_file,"%s/../ligplot/clig_%2.2d.pdb",directory,
                  num1);

          /* Form the name of the output html file */
          sprintf(out_file,"%s/%s.html",directory,id);

          /* Define the parameters */
          sprintf(prog_params,"%s %s -f %s -l %d.%d -c 0 -pdbsum1 %s",
                  pdb_code,template,out_file,num1,num2,relinks);
          if (metal == TRUE)
            strcat(prog_params," -option METAL");

          /* Run the program */
          sprintf(command,"cd %s && %s/showmain %s > %s",directory,
                  pdbsum_exe_dir,prog_params,log_file);
          call_system(command,dos);

          /* Rewrite the ligand_plugin.js javascript file to add in
             the ligand PDB coords where required */
          if (metal == FALSE)
            {
              /* Form the name of the output file */
              sprintf(out_file,"%s/ligand_%2.2d_plugin.js",
                      directory,num1);

              /* Rewrite the ligand_plugin.js file */
              rewrite_js_file(js_file,out_file,pdb_file,dfn_file);
            }
        }

      /* Get pointer to the next record */
      ligmet_ptr = ligmet_ptr->next_ligmet_ptr;
    }  

  /* Close progress output */
  printf("\n");
  fflush(stdout);
}
/***********************************************************************

gen_interfaces_html  -  Generate the INTERFACES pages

***********************************************************************/

void gen_interfaces_html(char *directory,char *pdb_code,
                         char *pdbsum_exe_dir,char *log_file,
                         struct chain *first_chain_ptr,int nchains,
                         char *relinks,struct c_file_data file_name_ptr)
{
  char out_file[FILENAME_LEN], template[FILENAME_LEN];
  char chains[3], command[FILENAME_LEN], prog_params[FILENAME_LEN];

  int iface;
  int dos;

  struct chain *chain_ptr;
  struct iface *iface_ptr;

  /* Initialise variables */
  iface = 0;
  dos = file_name_ptr.dos;

  /* Define the name of the template and output file */
  strcpy(template,"interfaces.html");
  sprintf(out_file,"%s/interfaces.html",directory);

  /* Show progress */
  printf("         - interface pages ...\n");

  /* Define the parameters */
  sprintf(prog_params,"%s %s -f %s -pdbsum1 %s -c 999",pdb_code,template,
          out_file,relinks);

  /* Run the program */
  sprintf(command,"cd %s && %s/showmain %s > %s",directory,pdbsum_exe_dir,
          prog_params,log_file);
  call_system(command,dos);

  /* Get first record */
  chain_ptr = first_chain_ptr;

  /* Loop through all the records */
  while (chain_ptr != NULL)
    {
      /* Get first interface record */
      iface_ptr = chain_ptr->first_iface_ptr;

      /* Loop through all the records */
      while (iface_ptr != NULL)
        {
          /* Increment interface number */
          iface++;

          /* Join the two chains */
          sprintf(chains,"%c%c",chain_ptr->chain,
                  iface_ptr->chain_ptr->chain);

          /* Form the name of the output html file */
          sprintf(out_file,"%s/interfaces_%s.html",directory,
                  chains);

          /* Define the parameters */
          sprintf(prog_params,
                  "%s %s -f %s -l %d -option RESIDUE -c 0 -pdbsum1 %s",
                  pdb_code,template,out_file,iface_ptr->iface_no,
                  relinks);

          /* Run the program */
          sprintf(command,"cd %s && %s/showmain %s > %s",directory,
                  pdbsum_exe_dir,prog_params,log_file);
          call_system(command,dos);

          /* Get pointer to the next record */
          iface_ptr = iface_ptr->next_iface_ptr;
        }  

      /* Get pointer to the next record */
      chain_ptr = chain_ptr->next_chain_ptr;
    }  
}
/***********************************************************************

gen_align_html  -  Generate the align pages

***********************************************************************/

void gen_align_html(char *directory,char *pdb_code,
                    char *pdbsum_exe_dir,char *log_file,
                    struct chain *first_chain_ptr,int nchains,
                    char *relinks,struct c_file_data file_name_ptr)
{
  char out_file[FILENAME_LEN], template[FILENAME_LEN];
  char chains[3], command[FILENAME_LEN], prog_params[FILENAME_LEN];

  int ichain;
  int dos;

  struct chain *chain_ptr;
  struct iface *iface_ptr;

  /* Initialise variables */
  ichain = 0;
  dos = file_name_ptr.dos;

  /* Define the name of the template and output file */
  strcpy(template,"align.html");

  /* Show progress */
  printf("         - alignment page ...\n");

  /* Define the parameters */
  //sprintf(params,"%s %s -f %s -pdbsum1 -c 999",pdb_code,template,
  //        out_file);

  /* Run the program */
  //sprintf(command,"cd %s && %s/showmain %s > %s",directory,pdbsum_exe_dir,
  //        prog_params,log_file);
  //call_system(command,dos);

  /* Get first record */
  chain_ptr = first_chain_ptr;

  /* Loop through all the records */
  while (chain_ptr != NULL)
    {
      /* Increment chain counter */
      ichain++;

      /* Form the name of the output file */
      sprintf(out_file,"%s/align_%c.html",directory,chain_ptr->chain);

      /* Define the parameters */
      sprintf(prog_params,
              "%s %s -f %s -l %d -chain %c -c 0 -pdbsum1 %s",pdb_code,
              template,out_file,ichain,chain_ptr->chain,relinks);

      /* Run the program */
      sprintf(command,"cd %s && %s/showmain %s >> %s",directory,
              pdbsum_exe_dir,prog_params,log_file);
      call_system(command,dos);

      /* Get pointer to the next record */
      chain_ptr = chain_ptr->next_chain_ptr;
    }  
}
/***********************************************************************

gen_clefts_html  -  Generate the clefts page

***********************************************************************/

void gen_clefts_html(char *directory,char *pdb_code,
                     char *pdbsum_exe_dir,char *log_file,char *relinks,
                     struct c_file_data file_name_ptr)
{
  char out_file[FILENAME_LEN], template[FILENAME_LEN];
  char command[FILENAME_LEN], prog_params[FILENAME_LEN];

  int dos;

  /* Initialise variables */
  dos = file_name_ptr.dos;

  /* Define the name of the template and output file */
  strcpy(template,"clefts.html");
  sprintf(out_file,"%s/clefts.html",directory);

  /* Show progress */
  printf("         - %s ...\n",template);

  /* Define the parameters */
  sprintf(prog_params,"%s %s -f %s -pdbsum1 %s",pdb_code,template,
          out_file,relinks);

  /* Run the program */
  sprintf(command,"cd %s && %s/showmain %s > %s",directory,pdbsum_exe_dir,
          prog_params,log_file);
  call_system(command,dos);
}
/***********************************************************************

gen_index_html  -  Generate the index page

***********************************************************************/

void gen_index_html(char *directory,char *pdb_code,
                    char *pdbsum_exe_dir,char *log_file,char *relinks,
                    struct c_file_data file_name_ptr,
                    struct parameters params)
{
  char out_file[FILENAME_LEN], template[FILENAME_LEN];
  char command[FILENAME_LEN], prog_params[FILENAME_LEN];
  char pdb_file[FILENAME_LEN], dfn_file[FILENAME_LEN];
  char js_file[FILENAME_LEN];

  int dos;

  /* Initialise variables */
  dos = file_name_ptr.dos;

  /* Define the name of the template and output file */
  strcpy(template,"main.html");
  sprintf(out_file,"%s/index.html",directory);

  /* Define the name of the pdbsum_plugin.js file */
  sprintf(js_file,"%s/js/pdbsum_plugin.js",
          file_name_ptr.other_dir[PDBSUM_REAL_TEMPLATES_DIR]);

  /* Show progress */
  printf("         - %s ...\n",template);

  /* Define the parameters */
  sprintf(prog_params,"%s %s -f %s -pdbsum1 %s",pdb_code,template,
          out_file,relinks);

  /* Run the program */
  sprintf(command,"cd %s && %s/showmain %s > %s",directory,
          pdbsum_exe_dir,prog_params,log_file);
  call_system(command,dos);

  /* Form the name of the PDB file */
  sprintf(pdb_file,"%s/../%s.pdb",directory,pdb_code);

  /* Form the name of the rasmol.dfn file */
  sprintf(dfn_file,"%s/../newsec/rasmol.dfn",directory);

  /* Form the name of the output .js file */
  sprintf(out_file,"%s/PDBsum1_plugin.js",directory);

  /* Rewrite the pdbsum_plugin.js file */
  rewrite_js_file(js_file,out_file,pdb_file,dfn_file);
}
/***********************************************************************

copy_images  -  Copy across all the image files that have been used


***********************************************************************/

void copy_images(char *pdir,struct c_file_data file_name_ptr)
{
  char command[FILENAME_LEN], file_name[FILENAME_LEN];
  char image_file[FILENAME_LEN], copy_file[FILENAME_LEN];
  char directory[FILENAME_LEN], templates_dir[FILENAME_LEN];
  char input_line[LINELEN + 1];

  int len;
  int dos;

  FILE *file_ptr;

  /* Initialise variables */
  dos = file_name_ptr.dos;
  strcpy(templates_dir,
         file_name_ptr.other_dir[PDBSUM_REAL_TEMPLATES_DIR]);

  /* Show progress */
  printf("         - copying images ...\n");

  /* Create the image directory */
  sprintf(directory,"%s/gif",pdir);
  create_directory(directory);

  /* Form the name of the relinks.txt file */
  sprintf(file_name,"%s/html/relinks.txt",pdir);

  /* Open the input file */
  file_ptr = fopen(file_name,"r");
  if (file_ptr == NULL)
    {
      printf("*** Warning. Unable to open file [%s] \"\n",file_name);
      return;
    }

  /* Loop through the file, copying to output file */
  while (fgets(input_line,LINELEN,file_ptr) != NULL)
    {
      /* Truncate the string to strip off the newline */
      len = string_chomp(input_line);

      /* Form the name of the image file */
      sprintf(image_file,"%s/gif/%s",templates_dir,input_line);

      /* Form the name of the copy file */
      sprintf(copy_file,"%s/%s",directory,input_line);

      /* Form command to perform the copy operation */
      sprintf(command,"%s %s %s",CP_COMMAND,image_file,copy_file);

      /* Perform the copy */
      call_system(command,dos);
    }

  /* Close the input file */
  fclose(file_ptr);

  /* Create the css directory */
  sprintf(directory,"%s/css",pdir);
  create_directory(directory);

  /* Copy across all the css files */
  sprintf(command,"%s %s/css/* %s",CP_COMMAND,templates_dir,directory);
  call_system(command,dos);

  /* Create the js directory */
  sprintf(directory,"%s/js",pdir);
  create_directory(directory);

  /* Copy across all the js files */
  sprintf(command,"%s %s/js/* %s",CP_COMMAND,templates_dir,directory);
  call_system(command,dos);
}
/***********************************************************************

run_showmain  -  Run showmain to generate all the HTML pages

***********************************************************************/

void run_showmain(char *fno,int process,char *pdir,char *pdb_code,
                  char *file_name,char *file_new,char *log_dir,int nprot,
                  int nnucl,struct r_rasdef *first_r_rasdef_ptr,
                  struct plot **plot_index_ptr,int nuniq,
                  struct chain *first_chain_ptr,int nchains,
                  struct ligmet *first_ligmet_ptr,int nligmet,
                  struct c_file_data file_name_ptr,
                  struct parameters params)
{
  char command[LINELEN + 1], directory[FILENAME_LEN];
  char log_file[FILENAME_LEN], relinks[20];
  char pdbsum_exe_dir[FILENAME_LEN];

  int len;
  int have_file, ok;

  float resol;

  /* Initialise variables */
  ok = TRUE;
  strcpy(pdbsum_exe_dir,file_name_ptr.other_dir[PDBSUM_EXE_DIR]);
  relinks[0] = '\0';

  /* Show progress */
  printf("   %s%d. In run_showmain ...\n",fno,process);

  /* Create the html directory */
  sprintf(directory,"%s/html",pdir);
  create_directory(directory);

  /* Form name of the log file */
  sprintf(log_file,"%s/showmain.log",log_dir);

  /* If image links to be relative, set the flag */
  if (params.relinks == TRUE)
    strcpy(relinks,"-relinks");

  /* == PROTEIN page(s) == */
  gen_protein_html(directory,pdb_code,pdbsum_exe_dir,log_file,
                   relinks,first_r_rasdef_ptr,file_name_ptr);

  /* == PROMOTIF pages == */
  if (nuniq > 0)
    gen_promotif_html(directory,pdb_code,pdbsum_exe_dir,log_file,
                      plot_index_ptr,nuniq,relinks,file_name_ptr);

  /* == PROCHECK page == */
  gen_procheck_html(directory,pdb_code,pdbsum_exe_dir,log_file,
                    relinks,file_name_ptr);

  /* == DNA/RNA page == */
  if (nprot > 0 && nnucl > 0)
    gen_dna_html(directory,pdb_code,pdbsum_exe_dir,log_file,relinks,
                 file_name_ptr);

  /* == Ligand and metal pages == */
  if (nligmet > 0)
    gen_ligands_html(directory,pdb_code,pdbsum_exe_dir,log_file,
                     first_ligmet_ptr,nligmet,relinks,file_name_ptr);

  /* == Prot-prot interface pages == */
  if (nchains > 1)
    gen_interfaces_html(directory,pdb_code,pdbsum_exe_dir,log_file,
                        first_chain_ptr,nchains,relinks,file_name_ptr);

  /* == Sequence alignment vs UniProt sequence pages == */
  gen_align_html(directory,pdb_code,pdbsum_exe_dir,log_file,
                 first_chain_ptr,nchains,relinks,file_name_ptr);

  /* == Clefts page == */
  gen_clefts_html(directory,pdb_code,pdbsum_exe_dir,log_file,relinks,
                  file_name_ptr);

  /* == Main index page == */
  gen_index_html(directory,pdb_code,pdbsum_exe_dir,log_file,relinks,
                 file_name_ptr,params);

  /* If have defined image ;inks as relative, copy across all the image
     files that have been used */
  if (params.relinks == TRUE)
    copy_images(pdir,file_name_ptr);
}
/***********************************************************************

get_pdb_details  -  Retrieve the structure's name and resolution, if
                    any, from the PDB file

***********************************************************************/

float get_pdb_details(char *fno,int file_no,char *pdir,
                      struct pdb *pdb_ptr,char *fname,
                      char *protein_name,int dos)
{
  char file_name[FILENAME_LEN], input_line[LINELEN + 1];

  char *string_ptr;
  
  int len;
  int done;

  float resol;

  FILE *file_ptr;

  /* Initialise variables */
  done = FALSE;
  resol = 2;

  /* Initialise protein name */
  strcpy(protein_name,"No name given");

  /* Form the name of the file */
  sprintf(file_name,"%s/%s",pdir,fname);
  if (file_no > 1)
    printf("\n\n");
  printf("%s Reading: %s\n",fno,file_name);
  fflush(stdout);

  /* Open the input file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    return(resol);

  /* Re-nitialise protein name */
  protein_name[0] = '\0';

  /* Loop through the file, line by line */
  while (fgets(input_line,LINELEN,file_ptr) != NULL && done == FALSE)
    {
      /* Truncate the input line */
      len = string_chomp(input_line);

      /* If this is the title, then save */
      if (!strncmp(input_line,"TITLE ",6))
        {
          /* Save or append the title */
          if (protein_name[0] == '\0')
            strcpy(protein_name,input_line + 10);
          else
            {
              strcat(protein_name," ");
              strcat(protein_name,input_line + 10);
            }
        }

      /* If this is the resolution, store the value */
      else if (!strncmp(input_line,"REMARK   2 RESOLUTION. ",23))
        {
          /* Remove the Angstroms at the end */
          string_ptr = strstr(input_line,"ANGSTROMS");
          if (string_ptr != NULL)
            {
              /* Terminate the line */
              string_ptr[0] = '\0';
                  
              /* Extract the resolution */
              resol = atof(input_line + 23);
            }
        }

      /* If we've reached the coordinate records, then we're done */
      else if (!strncmp(input_line,"ATOM  ",6) ||
               !strncmp(input_line,"HETATM",6))
        done = TRUE;
    }

  /* Close the file */
  fclose(file_ptr);

  /* If no name, then replace */
  if (protein_name[0] == '\0')
    strcpy(protein_name,"No name given");

  /* If have name, save it */
  else
    store_string(&(pdb_ptr->title),protein_name);

  /* Store the file name */
  store_string(&(pdb_ptr->file_name),fname);

  /* Return the resolution */
  return(resol);
}
/***********************************************************************

process_file  -  Process the given PDB file

***********************************************************************/

int process_file(int run_number,char *file_name,int file_no,
                 char *dir_name,char *pdb_code,struct pdb *pdb_ptr,
                 int npdb,FILE *runfile_ptr,
                 struct c_file_data file_name_ptr,
                 struct parameters params)
{
  char command[LINELEN + 1], fname[FILENAME_LEN], pdir[FILENAME_LEN];
  char file_new[FILENAME_LEN], log_dir[FILENAME_LEN];
  char directory[FILENAME_LEN], param_dir[FILENAME_LEN];
  char pdbsum_exe_dir[FILENAME_LEN], ps2pdf[FILENAME_LEN];
  char convert_exe[FILENAME_LEN], pymol_exe[FILENAME_LEN];
  char pdbsum1_dir[FILENAME_LEN], runs_dir[FILENAME_LEN];
  char orig_fname[FILENAME_LEN], protein_name[LINELEN + 1];
  char fno[20], relinks[20], sub_dir[3];

  int len, nchains, nligand, nligmet, nmetal, niface, nuniq, process;
  int nprot, nnucl;
  int ok, dos;

  float resol;

  struct chain *first_chain_ptr, *last_chain_ptr;
  struct ligmet *first_ligmet_ptr;
  struct r_ligand *first_r_ligand_ptr;
  struct r_rasdef *first_r_rasdef_ptr;

  struct plot **plot_index_ptr;

  /* Initialise variables */
  nchains = niface = nligmet = nuniq = 0;
  process = 0;
  dos = file_name_ptr.dos;
  ok = TRUE;
  strcpy(relinks," ");
  if (params.relinks == TRUE)
    strcpy(relinks,"-relinks");

  /* Initialise pointers */
  first_ligmet_ptr = NULL;
  first_r_ligand_ptr = NULL;
  first_r_rasdef_ptr = NULL;

  /* Initialise directory names */
  strcpy(param_dir,file_name_ptr.other_dir[PDBSUM_PARAMS_DIR]);
  strcpy(pdbsum_exe_dir,file_name_ptr.other_dir[PDBSUM_EXE_DIR]);
  strcpy(ps2pdf,file_name_ptr.exe_name[PS2PDF]);
  strcpy(pymol_exe,file_name_ptr.exe_name[PYMOL_EXE]);
  strcpy(convert_exe,file_name_ptr.exe_name[CONVERT_EXE]);
  strcpy(pdbsum1_dir,file_name_ptr.other_dir[PDBSUM1_RESULTS_DIR]);

  /* If more than file being processed, add file number to process
     lines */
  if (npdb == 1)
    fno[0] = '\0';
  else
    sprintf(fno,"%d.",file_no);

  /* Create PDBsum directory and transfer file */
  ok = initialise(file_name,dir_name,pdb_code,pdir,fname,
                  orig_fname,file_name_ptr);

  /* Get the 2nd and 3rd digits of the PDB code */
  len = strlen(pdb_code);
  strncpy(sub_dir,pdb_code + 1,2);
  sub_dir[2] = '\0';

  /* Write the file name and PDB code to the run file */
  fprintf(runfile_ptr,"ORIGINAL_FILE_NAME[%d] %s\n",file_no,orig_fname);
  fprintf(runfile_ptr,"PDB_CODE[%d] %s\n",file_no,pdb_code);
  fprintf(runfile_ptr,"PDB_DIR[%d] %s/%s\n",file_no,sub_dir,pdb_code);
  fprintf(runfile_ptr,"OK[%d] %s\n",file_no,T_F[ok]);

  /* If not OK, then return */
  if (ok == FALSE)
    return(ok);

  /* Get the structure name and resolution from the PDB file */
  resol = get_pdb_details(fno,file_no,pdir,pdb_ptr,fname,
                          protein_name,dos);

  /* Write title and file name to the run file */
  fprintf(runfile_ptr,"TITLE[%d] %s\n",file_no,protein_name);
  fprintf(runfile_ptr,"PDB_FILE_NAME[%d] %s\n",file_no,fname);

  /* Create a log directory */
  sprintf(log_dir,"%s/logs",pdir);
  create_directory(log_dir);

  /* Create the run data directory */
  sprintf(runs_dir,"%s/runs",pdbsum1_dir);
  create_directory(runs_dir);

  /* Run fixpdb on the input PDB file */
  process++;
  run_fixpdb(fno,process,pdir,pdb_code,fname,file_new,log_dir,
             file_name_ptr,dos);
  
  /* Run CLEAN on the input PDB file */
  process++;
  run_clean(fno,process,pdir,pdb_code,fname,file_new,log_dir,
            file_name_ptr);

  /* Run the remainder of the PROCHECK programs */
  process++;
  run_procheck(fno,process,pdir,pdb_code,fname,file_new,resol,log_dir,
               file_name_ptr);

  /* Run HBPLUS to compute all the H-bonds and non-bonded contacts */
  process++;
  run_hbplus(fno,process,pdir,pdb_code,fname,file_new,log_dir,
             file_name_ptr,params);

  /* Run newsec to analyse contents of the PDB file */
  process++;
  run_newsec(fno,process,pdir,pdb_code,fname,file_new,log_dir,param_dir,
             pdbsum_exe_dir,pymol_exe,file_name_ptr);

  /* Wait until newsec has created the rasmol.dfn file */
  if (dos == TRUE)
    wait_for_file(pdir,dos);
  
  /* Read in the definitions from the rasmol.dfn file */
  r_get_dfn_data(pdb_code,PDBSUM1,pdir,&first_r_rasdef_ptr,
                 &first_r_ligand_ptr,&nligand,&nmetal);

  /* Form the name of the newsec directory */
  sprintf(directory,"%s/newsec",pdir);

  /* Read in the chain equivalences from chains.dat */
  nchains = read_chain_equivs(directory,first_r_rasdef_ptr,
                              &first_chain_ptr,&last_chain_ptr,dos);

  /* Count protein and DNA chains */
  nprot = count_chains(first_r_rasdef_ptr,&nnucl);

  /* Read in the protein-protein interfaces, if any */
  niface = get_interfaces(pdir,&first_chain_ptr,&last_chain_ptr,
                          &nchains);

  /* If have any interfaces, generate the interface plots
     and PyMOL images */
  if (niface > 0)
    {
      /* Run ppcont to generate the schematic interaction plots */
      process++;
      run_ppcont(fno,process,pdir,pdb_code,fname,file_new,log_dir,
                 first_chain_ptr,file_name_ptr,params);
    }

  /* Run promotif to compute all the secondary structure motifs */
  process++;
  run_promotif(fno,process,pdir,pdb_code,fname,file_new,log_dir,
               file_name_ptr,params);
      
  /* Generate all the Promotif plots */
  process++;
  nuniq
    = run_psplot(fno,process,pdir,pdb_code,fname,file_new,log_dir,
                 first_r_rasdef_ptr,&plot_index_ptr,
                 file_name_ptr,params);

  /* Run resdata to compile the residue information */
  process++;
  run_resdata(fno,process,pdir,pdb_code,fname,file_new,log_dir,
              file_name_ptr,params);

  /* Run hera and topoldia to generate topology diagrams */
  process++;
  run_hera(fno,process,pdir,pdb_code,fname,file_new,log_dir,
           first_r_rasdef_ptr,file_name_ptr);

  /* Run speedfill to compute all the clefts */
  process++;
  if (params.max_cleft_chains == -1 ||
      nchains <= params.max_cleft_chains)
    run_speedfill(fno,process,pdir,pdb_code,log_dir,nchains,relinks,
                  file_name_ptr);

  /* Read in all the ligand and metal records from ligmet.run */
  nligmet = read_ligmet_run(pdir,&first_ligmet_ptr);

  /* Run ligplot for each ligand */
  if (nligmet > 0)
    {
      process++;
      run_ligplots(fno,process,pdir,pdb_code,file_new,log_dir,
                   pdbsum_exe_dir,param_dir,ps2pdf,convert_exe,
                   first_ligmet_ptr,nligmet,relinks,params.run_romlas,
                   file_name_ptr);
    }

  /* If file contains any DNA and protein sequences, run NUCPLOT */
  if (nprot > 0 && nnucl > 0)
    {
      process++;
      run_nucplot(fno,process,pdir,pdb_code,file_new,log_dir,
                  pdbsum_exe_dir,param_dir,ps2pdf,file_name_ptr);
    }

  /* Run seqdata to pull out UniProt refs and E.C. numbers  */
  process++;
  run_seqdata(fno,process,pdir,pdb_code,file_new,log_dir,pdbsum_exe_dir,
              file_name_ptr);

  /* Run uniplot to generate the mini domain diagram */
  process++;
  run_uniplot(fno,process,pdir,pdb_code,file_new,log_dir,pdbsum_exe_dir,
              file_name_ptr);

  /* Generate the wiring diagrams for all the protein chains */
  process++;
  run_wiring(fno,process,pdir,pdb_code,fname,file_new,log_dir,
             first_r_rasdef_ptr,pdbsum_exe_dir,file_name_ptr);

  /* Run pdbhtml to generate the database.dat file */
  process++;
  run_pdbhtml(fno,process,pdir,pdb_code,log_dir,run_number,
              file_name_ptr,params);

  /* Generate all the required RasMol scripts */
  if (params.run_romlas == TRUE)
    {
      process++;
      rasmol_scripts(fno,process,pdir,pdb_code,fname,file_new,log_dir,
                     first_r_rasdef_ptr,plot_index_ptr,nuniq,
                     first_chain_ptr,nchains,relinks,file_name_ptr);
    }

  /* Repeat for PyMOL */
  process++;
  pymol_scripts(fno,process,pdir,pdb_code,fname,file_new,log_dir,
                first_r_rasdef_ptr,plot_index_ptr,nuniq,
                first_chain_ptr,nchains,relinks,file_name_ptr);

  /* Run showmain to generate all the HTML pages */
  process++;
  run_showmain(fno,process,pdir,pdb_code,fname,file_new,log_dir,nprot,
               nnucl,first_r_rasdef_ptr,plot_index_ptr,nuniq,
               first_chain_ptr,nchains,first_ligmet_ptr,nligmet,
               file_name_ptr,params);

  return(ok);
}
/***********************************************************************

read_dataset_lst  -  Read in the files in dataset.lst

***********************************************************************/

int read_dataset_lst(int run_number,char *dataset,FILE *runfile_ptr,
                     struct pdb **fst_pdb_ptr,
                     struct c_file_data file_name_ptr,
                     struct parameters params)
{
  char dir_name[FILENAME_LEN], file_name[FILENAME_LEN];
  char input_line[LINELEN + 1];
  char pdb_code[MAX_PDB_CODE_LEN + 1];

  char *string_ptr, *token[MAX_TOKEN];

  int ipdb, len, nerrors, ntokens, npdb;
  int ok;

  struct pdb *pdb_ptr, *first_pdb_ptr, *last_pdb_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  ipdb = nerrors = npdb = 0;
  ok = TRUE;

  /* Initialise pdb pointers */
  pdb_ptr = first_pdb_ptr = last_pdb_ptr = NULL;

  /* Open the dataset.lst file */
  if ((file_ptr = fopen(dataset,"r")) ==  NULL)
    {
      printf("\n*** ERROR. Unable to open file [%s]\n",dataset);
      exit(1);
    }

  /* Read through the file to count number of entries */
  while (fgets(input_line,LINELEN,file_ptr) != NULL)
    {
      /* Check line-length */
      len = string_chomp(input_line);

      /* If line not blank, increment count */
      if (len > 0)
        npdb++;
    }

  /* Close file */
  fclose(file_ptr);

  /* Reopen it to read through properly */
  if ((file_ptr = fopen(dataset,"r")) ==  NULL)
    {
      printf("\n*** ERROR. Unable to open file [%s]\n",dataset);
      exit(1);
    }

  /* Loop through the file, line by line */
  while (fgets(input_line,LINELEN,file_ptr) != NULL)
    {
      /* Check line-length */
      len = string_chomp(input_line);

      /* Read in all the tokens of this line */
      ntokens = get_tokens(input_line,token,MAX_TOKEN,' ');

      /* If have at least one token, proceed */
      if (ntokens > 0)
        {
          /* Take the first token to be the file name */
          strcpy(file_name,token[0]);

          /* If have two tokens, take second to be the PDB code */
          if (ntokens == 2 && strlen(token[1]) == 4)
            strcpy(pdb_code,token[1]);

          /* Otherwise try to extract the PDB code from the file name,
           or else generate a new code */
          else
            ok = extract_pdb_code(file_name,dir_name,pdb_code,
                                  file_name_ptr);

          /* If file name OK, process this entry */
          if (ok == TRUE)
            {
              /* Increment count */
              ipdb++;

              /* Create a new pdb entry */
              pdb_ptr
                = create_pdb_entry(&first_pdb_ptr,&last_pdb_ptr,
                                   pdb_code);

              /* Process to file */
              ok = process_file(run_number,file_name,ipdb,dir_name,
                                pdb_code,pdb_ptr,npdb,runfile_ptr,
                                file_name_ptr,params);
            }
          else
            nerrors++;
        }

      /* Free up memory taken up by the tokens */
      free_tokens(token,ntokens);
    }

  /* Close file */
  fclose(file_ptr);

  /* Show number of PDB codes to be processed */
  printf("Number of PDB codes to be processed: %8d\n",npdb);

  /* If no codes found, abort */
  if (npdb == 0)
    {
      printf("*** ERROR. No PDB codes to be processed!\n");
      exit(1);
    }

  /* Show number of files processed */
  printf("Number of PDB files processed:    %8d\n",npdb);
  printf("Number of filename errors:        %8d\n",nerrors);

  /* Return the first of the PDB entries */
  *fst_pdb_ptr = first_pdb_ptr;

  /* Return the number of files processed in */
  return(npdb);
}
/***********************************************************************

write_index_html  -  Run showmain to generate the index.html and runX.html
                     pages

***********************************************************************/

void write_index_html(int run_number,struct c_file_data file_name_ptr,
                      struct parameters params)
{
  char command[LINELEN + 1], directory[FILENAME_LEN];
  char log_file[FILENAME_LEN], template[FILENAME_LEN];
  char pdbsum_exe_dir[FILENAME_LEN], pdbsum1_dir[FILENAME_LEN];
  char prog_params[FILENAME_LEN], out_file[FILENAME_LEN];

  int len;
  int have_file, ok, dos;

  float resol;

  /* Initialise variables */
  ok = TRUE;
  strcpy(pdbsum_exe_dir,file_name_ptr.other_dir[PDBSUM_EXE_DIR]);
  strcpy(pdbsum1_dir,file_name_ptr.other_dir[PDBSUM1_RESULTS_DIR]);
  dos = file_name_ptr.dos;

  /* Show progress */
  printf("In write_index_html ...\n");

  /* Define the run directory */
  sprintf(directory,"%s/runs",pdbsum1_dir);
  create_directory(directory);

  /* Form name of the log file */
  sprintf(log_file,"%s/run_index.log",directory);

  /* Define the name of the template and output file */
  strcpy(template,"run.html");
  sprintf(out_file,"%s/index.html",pdbsum1_dir);

  /* Show progress */
  printf("         - %s ...\n",template);

  /* Define the parameters */
  sprintf(prog_params,"n/a %s -f %s -pdbsum1 -doc -s %d",template,out_file,
          run_number);

  /* If using relative image links, define relative path to common
     images */
  if (params.relinks == TRUE)
    strcat(prog_params," -relinks");
  
  /* Run the program */
  sprintf(command,"cd %s && %s/showmain %s >> %s",directory,pdbsum_exe_dir,
          prog_params,log_file);
  call_system(command,dos);

  /* Form name of the copy file, runs/runXX.html, where XX is the
     run number */
  sprintf(out_file,"run%d.html",run_number);

  /* Redefine the parameters */
  sprintf(prog_params,"n/a %s -f %s -pdbsum1 -doc -s %d -option UPDIR",
          template,out_file,run_number);

  /* Run showmain again */
  sprintf(command,"cd %s && %s/showmain %s >> %s",directory,pdbsum_exe_dir,
          prog_params,log_file);
  call_system(command,dos);
}
/***********************************************************************

create_run_entry  -  Create a run entry

***********************************************************************/

struct run *create_run_entry(struct run **fst_run_ptr,
                             struct run **lst_run_ptr,char *input_line)
{
  int isort;

  struct run *run_ptr, *first_run_ptr, *last_run_ptr;

  /* Initialise run pointers */
  run_ptr = NULL;
  first_run_ptr = *fst_run_ptr;
  last_run_ptr = *lst_run_ptr;

  /* Create run entry */
  run_ptr = (struct run *) malloc(sizeof(struct run));
  if (run_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for struct run\n");
      printf("***        Program analrefs terminated with error.\n");
      exit (1);
    }

  /* Store the run details */
  store_string(&(run_ptr->input_line),input_line);

  /* Initialise pointer to next entry */
  run_ptr->next_run_ptr = NULL;

  /* If this is the first entry stored, then update
     pointer to head of linked list */
  if (first_run_ptr == NULL)
    first_run_ptr = run_ptr;

  /* Otherwise, make previous entry point to the current
     one */
  else
    last_run_ptr->next_run_ptr = run_ptr;
                      
  /* Save the current residue's pointer */
  last_run_ptr = run_ptr;

  /* Return current pointers */
  *fst_run_ptr = first_run_ptr;
  *lst_run_ptr = last_run_ptr;
  return(run_ptr);
}
/***********************************************************************

remove_excess_blanks  -  Replace any runs of blanks with a single one

***********************************************************************/

void remove_excess_blanks(char *string)
{
  char ch, last_ch, tmp[LINELEN];

  int i, len, next_pos;
  int done;

  /* Initialise variables */
  done = FALSE;
  i = next_pos = 0;
  last_ch = ' ';
  tmp[0] = '\0';

  /* Get the length of the string */
  len = strlen(string);

  /* Loop through the characters in the string */
  for (i = 0; i < len; i++)
    {
      /* Get this character */
      ch = string[i];

      /* If not a space, or last character was not a space, copy
         across */
      if (ch != ' ' || last_ch != ' ')
        {
          tmp[next_pos] = ch;
          next_pos++;
        }

      /* Save character  */
      last_ch = ch;
    }

  /* Close off the strig */
  tmp[next_pos] = '\0';

  /* Copy back */
  strcpy(string,tmp);
}
/***********************************************************************

rewrite_runs_dbt  -  Rewrite the data file holding details of all past
                     runs of PDBsum1

***********************************************************************/
    
void rewrite_runs_dbt(int run_number,char *date,struct pdb *first_pdb_ptr,
                      int npdb,struct c_file_data file_name_ptr)
{
  char command[LINELEN + 1], input_line[LONG_LINELEN + 1];
  char pdbsum1_dir[FILENAME_LEN], runs_dir[FILENAME_LEN], sub_dir[3];
  char file_name[FILENAME_LEN], string[LINELEN + 1], title[LINELEN + 1];

  char *directories, *string_ptr, *titles;

  int len, ncodes, nruns;
  int ok;

  struct pdb *pdb_ptr;
  struct run *run_ptr, *first_run_ptr, *last_run_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  nruns = 0;
  directories = titles = &null;

  /* Initialise pointers */
  run_ptr = first_run_ptr = last_run_ptr = NULL;
  
  /* Initialise directory name */
  strcpy(pdbsum1_dir,file_name_ptr.other_dir[PDBSUM1_RESULTS_DIR]);

  /* Define name of the run data directory */
  sprintf(runs_dir,"%s/runs",pdbsum1_dir);

  /* Form the name of the run list, runs.dbt */
  sprintf(file_name,"%s/runs.dbt",runs_dir);

  /* Open the runs.dbt file */
  if ((file_ptr = fopen(file_name,"r")) !=  NULL)
    {
      /* Loop while reading in the data */
      while (fgets(input_line,LONG_LINELEN,file_ptr) != NULL)
        {
          /* Truncate the input line */
          len = string_chomp(input_line);
          
          /* Create a record to store the run details */
          run_ptr
            = create_run_entry(&first_run_ptr,&last_run_ptr,
                               input_line);

          /* Increment count of run details */
          nruns++;
        }

      /* Close the file */
      fclose(file_ptr);
    }

  /* Open the output file */
  file_ptr = open_output_file(file_name,TRUE);

  /* If this is the first run, then write out the header records */
  if (nruns == 0)
    {
      /* Write out the header records */
      fprintf(file_ptr,"START:RESULTS\n");
      fprintf(file_ptr,
              "RUN_NUMBER::DATE::NPDB::PDB_CODES::PDB_DIR::TITLES\n");
    }

  /* Otherwise, write out all the run details just read in from
     run.dbt */
  else
    {
      /* Get pointer to the first run record */
      run_ptr = first_run_ptr;

      /* Loop over all the runs to find the one just read in */
      while (run_ptr != NULL)
        {
          /* Write out this line */
          fprintf(file_ptr,"%s\n",run_ptr->input_line);

          /* Get pointer to the next run */
          run_ptr = run_ptr->next_run_ptr;
        }
    }

  /* Write out the current run's details */
  fprintf(file_ptr,"%d::%s::%d::",run_number,date,npdb);

  /* Get the first of the PDB entries */
  pdb_ptr = first_pdb_ptr;
  ncodes = 0;

  /* Loop over all the PDBs to write out the details */
  while (pdb_ptr != NULL && ncodes < MAX_RUN_CODES)
    {
      /* Increment count of codes */
      ncodes++;

      /* Write out the PDB code */
      fprintf(file_ptr,"%s;",pdb_ptr->pdb_code);

      /* Get the 2nd and 3rd digits of the PDB code */
      len = strlen(pdb_ptr->pdb_code);
      strncpy(sub_dir,pdb_ptr->pdb_code + 1,2);
      sub_dir[2] = '\0';

      /* Save the directory */
      sprintf(string,"%s/%s;",sub_dir,pdb_ptr->pdb_code);
      append_string(&directories,string);

      /* If we don't have a title, use the original file name */
      if (pdb_ptr->title != &null)
        strcpy(title,pdb_ptr->title);
      else
        strcpy(title,pdb_ptr->file_name);

      /* Check for any semi-colons and remove if found */
      string_ptr = strchr(title,';');
      if (string_ptr != NULL)
        replace_char(title,';',',');

      /* Remove excess blanks */
      remove_excess_blanks(title);

      /* Save the title */
      sprintf(string,"%s;",title);
      append_string(&titles,string);

      /* Get pointer to the next pdb */
      pdb_ptr = pdb_ptr->next_pdb_ptr;
    }

  /* Write out the directories and titles */
  fprintf(file_ptr,"::%s::%s\n",directories,titles);

  /* Close the file */
  fclose(file_ptr);
}
/***********************************************************************

write_allruns_html  -  Run showmain to generate the list of all PDBsum1
                       runs, allruns.html

***********************************************************************/

void write_allruns_html(struct c_file_data file_name_ptr,
                        struct parameters params)
{
  char command[LINELEN + 1], directory[FILENAME_LEN];
  char log_file[FILENAME_LEN], template[FILENAME_LEN];
  char pdbsum_exe_dir[FILENAME_LEN], pdbsum1_dir[FILENAME_LEN];
  char prog_params[FILENAME_LEN], out_file[FILENAME_LEN];

  int len;
  int have_file, ok, dos;

  float resol;

  /* Initialise variables */
  ok = TRUE;
  strcpy(pdbsum_exe_dir,file_name_ptr.other_dir[PDBSUM_EXE_DIR]);
  strcpy(pdbsum1_dir,file_name_ptr.other_dir[PDBSUM1_RESULTS_DIR]);
  dos = file_name_ptr.dos;

  /* Show progress */
  printf("In write_allruns_html ...\n");

  /* Define the run directory */
  sprintf(directory,"%s/runs",pdbsum1_dir);
  create_directory(directory);

  /* Form name of the log file */
  sprintf(log_file,"%s/run_index.log",directory);

  /* Define the name of the template and output file */
  strcpy(template,"allruns.html");
  sprintf(out_file,"%s/runs/allruns.html",pdbsum1_dir);

  /* Show progress */
  printf("         - %s ...\n",template);

  /* Define the parameters */
  sprintf(prog_params,"n/a %s -f %s -pdbsum1 -doc",template,out_file);

  /* Run the program */
  sprintf(command,"cd %s && %s/showmain %s >> %s",directory,
          pdbsum_exe_dir,prog_params,log_file);
  call_system(command,dos);
}

/***********************************************************************

                             M   A   I   N

***********************************************************************/

int main(int argc,char *argv[])
{
  char dir_name[FILENAME_LEN], pdb_code[MAX_PDB_CODE_LEN + 1];
  char date[64], index_html[FILENAME_LEN];
/* v.1.0.1 --> */
  char heading[LONG_LINELEN];
/* <-- v.1.0.1 */

  char *string_ptr;

  char dquote = '"';

  int i, len, npdb, operating_system, os_type, run_number;
  int ok;

  struct c_file_data file_name_ptr;
  struct parameters params;

  struct pdb *pdb_ptr, *first_pdb_ptr, *last_pdb_ptr;

  FILE *runfile_ptr;

  /* Initialise */
  ok = TRUE;
  pdb_ptr = NULL;

  /* Initialise pointers */
  pdb_ptr = first_pdb_ptr = last_pdb_ptr = NULL;

  /* Interpret the command-line parameters */
  get_command_arguments(argv,argc - 1,&params);

  /* Read in all the directory paths from the parameter file, CATHPARAM */
  c_get_paths(&file_name_ptr,TRUE);

  /* Initialise operating system specific commands */
  operating_system = -1;
  for (i = 0; i < NOP_SYS; i++)
    if (!strcmp(file_name_ptr.other_dir[OPERATING_SYSTEM],OS[i]))
      operating_system = i;

  /* Define the os-specific commands */
  if (operating_system > -1)
    {
      /* Set the commands depending on whether DOS or non-DOS */
      os_type = NOT_DOS;
      if (file_name_ptr.dos == TRUE)
        os_type = DOS;

      /* Set the commands accordingly */
      strcpy(CP_COMMAND,CP_COPY[os_type]);
      strcpy(MV_COMMAND,MV_REN[os_type]);
      strcpy(RM_COMMAND,RM_DEL[os_type]);
      strcpy(SLEEP_COMMAND,SLEEP[os_type]);
    }
  else
    {
      printf("*** ERROR. Unrecognised Operating System in CATHPARAM "
             "parameter file: %s\n",
             file_name_ptr.other_dir[OPERATING_SYSTEM]);
      printf("           Must be one of: ");
      for (i = 0; i < NOP_SYS; i++)
        {
          if (i > 0)
            printf(", ");
          printf("%s",OS[i]);
        }
      printf("\n");
      exit(1);
    }

  /* Initialise for current run */
  runfile_ptr = initialise_run(&run_number,date,file_name_ptr,
                               params);

  /* If only a single file to process, call process routine */
  if (params.single_file == TRUE)
    {
      /* Show selection */
/* v.1.0.1 --> */
//      printf("Running PDBsum1 on file %s\n",params.file);
//      printf("------------------------");
//      len = strlen(params.file);
//      for (i = 0; i < len; i++)
//        printf("-");
//      printf("\n");
      sprintf(heading,"Running PDBsum1 (%s) on file %s",
	      params.version,params.file);
      printf("%s\n",heading);

      /* Draw underscores */
      len = strlen(heading);
      for (i = 0; i < len; i++)
	printf("-");
      printf("\n");
/* <-- v.1.0.1 */

      /* Show PDB code, if we have one */
      if (params.pdb_code[0] != '\0')
        printf("PDB code:  %s\n",params.pdb_code);
      printf("\n");

      /* If we don't already have the PDB code, try to extract it from
         the file name, or generate a new one */
      if (params.pdb_code[0] == '\0')
        ok = extract_pdb_code(params.file,dir_name,pdb_code,
                              file_name_ptr);
      else
        strcpy(pdb_code,params.pdb_code);

      /* Create a new pdb entry */
      pdb_ptr
        = create_pdb_entry(&first_pdb_ptr,&last_pdb_ptr,pdb_code);

      /* Set number of files processed to 1 */
      npdb = 1;

      /* If file name OK, process this entry */
      if (ok == TRUE)
        process_file(run_number,params.file,npdb,dir_name,pdb_code,
                     pdb_ptr,npdb,runfile_ptr,file_name_ptr,params);
    }

  /* Otherwise, process each of the files listed in dataset.lst in turn */
  else
    {
      /* Show selection */
      printf("Running PDBsum1 on dataset of files\n");
      printf("-----------------------------------\n");
      printf("\n");

      /* Read through dataset.lst and process each file */
      npdb = read_dataset_lst(run_number,params.file,runfile_ptr,
                              &first_pdb_ptr,file_name_ptr,params);
    }

  /* Write number of files processed to the run file */
  fprintf(runfile_ptr,"NPDB %d\n",npdb);

  /* Close the run file */
  fclose(runfile_ptr);

  /* Generate the html page with links to the results */
  write_index_html(run_number,file_name_ptr,params);

  /* Update the list of runs file, runs.dbt, to add in the details
     of the current run */
  rewrite_runs_dbt(run_number,date,first_pdb_ptr,npdb,file_name_ptr);

  /* Use the runs.dbt to generate the page listing all PDBsum1 runs */
  write_allruns_html(file_name_ptr,params);

  /* Form link to results page */
  strcpy(index_html,file_name_ptr.other_dir[PDBSUM_TEMPLATE]);
  string_ptr = strstr(index_html,"[bc]");
  if (string_ptr != NULL)
    string_ptr[0] = '\0';
  strcat(index_html,"index.html");

  /* Show the link */
  printf("\n");
  printf("Results at: file://%s\n",index_html);

  return(0);
}

