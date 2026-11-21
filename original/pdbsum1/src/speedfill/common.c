/***********************************************************************

  common.c - Common routines used by all programs

***********************************************************************

Date:-         22 Oct 2009

Written by:-   Roman Laskowski

----------------------------------------------------------------------*/

/* Compilation
   -----------

cc -c common.c

*/

/* Routines
   --------

   -> store_string
   -> append_string
   -> apply_pixel
         -> convert_to_rgb
   -> list_paths (for testing only)
   -> call_system
   -> c_get_paths
         -> get_param (in reapar.c)
   -> check_file
   -> colorise
   -> convert_rgb_to_hex
   -> convert_to_grey
   -> copy_file
         -> open_output_file
         -> string_truncate
   -> create_directory
         -> check_directory (no longer used)
         -> stat
   -> empty_check
   -> extract_res_name
         -> to_upper
         -> to_lower
         -> get_aa_code
         -> get_res_name
   -> fix_upper
   -> format_exp_no
   -> form_filename
         -> store_string
   -> file_diff
   -> find_pdbsum_dir
         -> form_filename
   -> file_months
   -> fix_token
   -> free_tokens
   -> get_aa_code
   -> get_chain_colour
         -> get_colour_rgb
   -> get_date
         -> to_upper
   -> get_date_and_time
         -> to_upper
   -> get_font
         -> read_ppm_file
               -> string_chomp
               -> tokenize
         -> get_tokens
         -> store_string
         -> create_letter_record
   -> get_pdbsum1_dir
   -> get_random_number
   -> get_res_name
   -> get_rgb
   -> get_irgb_colour
   -> get_pixel_shade
   -> get_tokens
   -> is_number
   -> is_numeric
   -> move_file
   -> open_output_file
   -> output_ppm_file
   -> parse_uniprot_name
   -> perform_splice
   -> read_cath_titles
         -> string_chomp
         -> create_c_cath_entry
   -> repeat_splice
         -> perform_splice
   -> read_ppm_file
         -> get_tokens
   -> remove_blanks
   -> remove_char
   -> remove_leading_blanks
   -> replace_char
   -> rgb_to_int
   -> run_wget
         -> call_system (in common.c)
   -> strccnt
   -> string_chomp
   -> string_chop
   -> string_truncate
   -> to_lower
   -> to_upper
   -> tokenize
   -> type_set
   -> write_ppm_file
   -> convert_to_gif
   -> extract_domain_info
         -> get_tokens
         -> format_residue_number
         -> free_tokens
   -> get_c_scores
         -> create_c_score_entry


*/


#include "common.h"
#include "dcommon.h"

/* Prototypes
   ---------- */

/* In reapar.c */
int get_param(char *token_name,char *file_name,int filename_len,
              int exit_on_error);

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

convert_to_rgb  -  Convert colour to RGB triplet

***********************************************************************/

void convert_to_rgb(int colour,int *rgb)
{
  /* Perform the conversion */
  rgb[0] = colour / (256 * 256);
  colour = colour % (256 * 256);
  rgb[1] = colour / 256;
  rgb[2] = colour % 256;
}
/***********************************************************************

apply_pixel  -  Apply given pixel colour over background colour

***********************************************************************/

int apply_pixel(int bg_colour,float grey_shade,int *text_rgb)
{
  int colour, i, rgb[3];

  /* Extract the components of the background colour */
  convert_to_rgb(bg_colour,rgb);

  /* Compute the new colour */
  for (i = 0; i < 3; i++)
    {
      /* Adjust according to which is larger */
      if (text_rgb[i] > rgb[i])
        rgb[i]
          = (int) (rgb[i] + (text_rgb[i] - rgb[i]) * grey_shade);
      else
        rgb[i]
          = (int) (rgb[i] - (rgb[i] - text_rgb[i]) * grey_shade);
    }

  /* Convert the new RGB value to colour value */
  colour = (int) (256 * 256 * rgb[0] + 256 * rgb[1] + rgb[2]);

  /* Return the pixel's new colour value */
  return(colour);
}
/***********************************************************************

append_string  -  Append given string to existing string, allocating
                  memory for the expanded string, storing pointer to it

***********************************************************************/

void append_string(char **string_ptr,char *string)
{
  char *new;

  int len;

  /* Get the total length of the two strings */
  len = strlen(*string_ptr) + strlen(string);

  /* Allocate memory for the new string */
  new = (char *) malloc(sizeof(char)*(len + 1));

  /* Append second string onto first */
  strcpy(new,*string_ptr);
  strcat(new,string);

  /* Free up memory taken by old string */
  if (strlen(*string_ptr) > 0)
    free(*string_ptr);

  /* Return the new string pointer */
  *string_ptr = new;
}
/***********************************************************************

list_paths  -  Print out the parameters read in from CATHPARAM

***********************************************************************/

void list_paths(struct c_file_data *file_name_ptr)
{
  char *string_ptr;

  int i;

  /* Loop over the PDB templates */
  printf("PDB templates\n");
  for (i = 0; i < NPDB_DIR; i++)
    {
      printf("%d. %s = %s\n",i,PDB_DIR_TEMPLATE[i],
             file_name_ptr->pdb_template[i]);
      fflush(stdout);
    }

  /* PDB exe and gen directories */
  printf("\n");
  printf("PDB exe: %s\n",file_name_ptr->pdbsum_exe_dir);
  fflush(stdout);
  printf("PDB gen: %s\n",file_name_ptr->pdbsum_gen_dir);
  fflush(stdout);

  /* PDBsum templates */
  printf("\n");
  printf("PDBsum template and data dirs:\n");
  for (i = 0; i < NPDBSUM_DIR; i++)
    {
      printf("%2d. %s = %s\n",i,PDBSUM_DIR_TEMPLATE[i],
             file_name_ptr->pdbsum_template[i]);
      printf("     %s = %s\n",PDBSUM_DATA_DIR[i],
             file_name_ptr->pdbsum_data_dir[i]);
      fflush(stdout);
    }
  
  /* Executable programs */
  printf("\n");
  printf("EXE dirs:\n");
  for (i = 0; i < NEXE_DIR; i++)
    {
      printf("%2d. %s = %s\n",i,EXE_NAME[i],
             file_name_ptr->exe_name[i]);
      fflush(stdout);
    }

  /* Other directories */
  printf("\n");
  printf("Other dirs:\n");
  for (i = 0; i < NOTHER_DIR; i++)
    {
      printf("%2d. %s = %s\n",i,OTHER_DIR_NAME[i],
             file_name_ptr->other_dir[i]);
      fflush(stdout);
    }
}
/***********************************************************************

c_get_paths  -  Read in the paths for the common PDBsum and PDB files

***********************************************************************/

void c_get_paths(struct c_file_data *file_name_ptr,int run_64)
{
  char file_name[C_FILENAME_LEN];

  int error_status, idir;
  static int exit_on_error = FALSE;

  /* Initialise DOS flag */
  file_name_ptr->dos = FALSE;
  
  /* PDB directory pathnames */
  for (idir = 0; idir < NPDB_DIR; idir++)
    {
      error_status = get_param(PDB_DIR_TEMPLATE[idir],file_name,
                               C_FILENAME_LEN,exit_on_error);
      store_string(&(file_name_ptr->pdb_template[idir]),file_name);
    }

  /* PDBsum directory pathnames */
  for (idir = 0; idir < NPDBSUM_DIR; idir++)
    {
      /* Get the PDBsum directory template */
      error_status = get_param(PDBSUM_DIR_TEMPLATE[idir],file_name,
                               C_FILENAME_LEN,exit_on_error);
      store_string(&(file_name_ptr->pdbsum_template[idir]),file_name);

      /* Get the name of the corresponding data directory */
      error_status = get_param(PDBSUM_DATA_DIR[idir],file_name,
                               C_FILENAME_LEN,exit_on_error);
      store_string(&(file_name_ptr->pdbsum_data_dir[idir]),file_name);
    }

  /* PDBsum generate executable directories */
  if (run_64 == TRUE)
    error_status = get_param("PDBSUM_EXE_64_DIR",file_name,
                             C_FILENAME_LEN,exit_on_error);
  else
    error_status = get_param("PDBSUM_EXE_DIR",file_name,
                             C_FILENAME_LEN,exit_on_error);
  store_string(&(file_name_ptr->pdbsum_exe_dir),file_name);

  /* PDBsum generate executable directories */
  if (run_64 == TRUE)
    error_status = get_param("PDBSUM_GEN_64_DIR",file_name,
                             C_FILENAME_LEN,exit_on_error);
  else
    error_status = get_param("PDBSUM_GEN_DIR",file_name,
                             C_FILENAME_LEN,exit_on_error);
  store_string(&(file_name_ptr->pdbsum_gen_dir),file_name);

  /* Executables */
  for (idir = 0; idir < NEXE_DIR; idir++)
    {
      error_status = get_param(EXE_NAME[idir],file_name,
                               C_FILENAME_LEN,exit_on_error);
      store_string(&(file_name_ptr->exe_name[idir]),file_name);
    }

  /* Other file and pathnames */
  for (idir = 0; idir < NOTHER_DIR; idir++)
    {
      error_status = get_param(OTHER_DIR_NAME[idir],file_name,
                               C_FILENAME_LEN,exit_on_error);
      store_string(&(file_name_ptr->other_dir[idir]),file_name);
    }

  /* Check if the path is a DOS-style one with backslashes */
  if (strstr(file_name_ptr->other_dir[PDBSUM_EXE_DIR],":/") != NULL)
    file_name_ptr->dos = TRUE;
  
  /* For checking */
  // list_paths(file_name_ptr);
}
/***********************************************************************

check_file  -  Check whether the given file exists

***********************************************************************/

int check_file(char *file_name)
{
  int have_file;

  FILE *file_ptr;

  /* Initialise */
  have_file = FALSE;

  /* Open the file to see if it exists */
  if ((file_ptr = fopen(file_name,"r")) !=  NULL)
    {
      /* Have the file */
      have_file = TRUE;

      /* Close it */
      fclose(file_ptr);
    }

  /* Return whether file exists */
  return(have_file);
}
/***********************************************************************

empty_check  -  Check whether the given file exists and isn't empty

***********************************************************************/

int empty_check(char *file_name)
{
  char input_line[C_LINELEN];

  int len, lines;

  int have_file;

  FILE *file_ptr;

  /* Initialise */
  have_file = FALSE;
  lines = 0;

  /* Open the file to see if it exists */
  if ((file_ptr = fopen(file_name,"r")) !=  NULL)
    {
      /* Loop through the file processing each record */
      while (fgets(input_line,C_LINELEN,file_ptr) != NULL &&
             lines == 0)
        {
          /* If line not empty, increment line count */
          len = strlen(input_line);
          if (len > 0)
            lines++;
        }

      /* Close it */
      fclose(file_ptr);
    }

  /* If file not empty, increment line count */
  if (lines > 0)
    have_file = TRUE;

  /* Return whether file exists */
  return(have_file);
}
/***********************************************************************

colorise  -  Convert given grey shade to colour

***********************************************************************/

int colorise(int grey_value,int *ref_rgb)
{
  int colour, i, rgb[3];

  float calc;

  /* Perform conversion */
  for (i = 0; i < 3; i++)
    {
      rgb[i] = ref_rgb[i]
        + (255 - ref_rgb[i]) * (float) grey_value / 255.0;
      if (rgb[i] > 255)
        rgb[i] = 255;
    }

  /* Convert to colour number */
  colour = 256 * 256 * rgb[0] + 256 * rgb[1] + rgb[2];

  /* Return colour */
  return(colour);
}
/***********************************************************************

convert_rgb_to_hex  -  Convert RGB triplet to the corresponding hex code

***********************************************************************/

void convert_rgb_to_hex(int *irgb,char *hex)
{
  static char conv_ch[] = "0123456789ABCDEF";

  char code[3];

  int i, number;

  /* Initialize */
  strcpy(hex,"#");

  /* Loop over the RGB components of the colour */
  for (i = 0; i < 3; i++)
    {
      /* Initialise code */
      strcpy(code,"FF");

      /* If value is valid, compute the code */
      if (irgb[i] > -1 && irgb[i] < 256)
        {
          /* Get first character */
          number = irgb[i] / 16;
          code[0] = conv_ch[number];

          /* Get the second character */
          number = irgb[i] % 16;
          code[1] = conv_ch[number];

          /* Terminate the code */
          code[2] = '\0';
        }

      /* Add to end of hex string */
      strcat(hex,code);
    }
}
/***********************************************************************

convert_to_grey  -  Convert given colour to grey shade

***********************************************************************/

int convert_to_grey(int colour,int *gv)
{
  int i, grey_value, rgb[3];

  static const int grey_component[3] = { 76, 150, 29 };

  /* Convert colour into an RGB triplet */
  rgb[0] = colour / (256 * 256);
  colour = colour % (256 * 256);
  rgb[1] = colour / 256;
  rgb[2] = colour % 256;

  /* Calculate new RGB value */
  grey_value = 0;
  for (i = 0; i < 3; i++)
    grey_value
      = grey_value + grey_component[i] * (float) rgb[i] / 255.0;
  if (grey_value > 255)
    grey_value = 255;

  /* Compute colour number of this grey value */
  colour = 256 * 256 * grey_value + 256 * grey_value + grey_value;

  /* Return grey value */
  *gv = grey_value;

  /* Return colour */
  return(colour);
}
/***********************************************************************

open_output_file  -  Open the output file

***********************************************************************/

FILE *open_output_file(char *out_name,int terminate)
{
  FILE *outfile_ptr;

  /* Open the output file, or write to stdout */
  if (out_name[0] == '\0')
    outfile_ptr = stdout;
  else
    outfile_ptr = fopen(out_name,"w");

  /* Warn or terminate */
  if (outfile_ptr == NULL)
    {
      if (terminate == TRUE)
        {
          printf("*** ERROR. Unable to open output file: %s\n",
                 out_name);
          printf("***        Program aborted\n");
          exit(-1);
        }
      else
        printf("*** Warning. Unable to open output file: %s\n",
               out_name);
    }

  /* Return the pointer to the output file */
  return(outfile_ptr);
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

copy_file  -  Copy the given file across

***********************************************************************/

void copy_file(char *infile,char *outfile)
{
  char input_line[C_LINELEN + 1];

  FILE *file_ptr, *outfile_ptr;

  /* Open the output file */
  outfile_ptr = open_output_file(outfile,TRUE);

  /* Open the input file */
  if ((file_ptr = fopen(infile,"r")) !=  NULL)
    {
      /* Loop while reading in records from the file */
      while (fgets(input_line,C_LINELEN,file_ptr) != NULL)
        {
          /* Truncate the line at the last non-blank character */
          string_truncate(input_line,C_LINELEN);

          /* Write the line out */
          fprintf(outfile_ptr,"%s\n",input_line);
        }

      /* Close the file */
      fclose(file_ptr);
    }

  /* Show error message */
  else
    {
      printf("*** Warning. Unable to open input file [%s]\n",
             infile);
    }

  /* Close the output file */
  fclose(outfile_ptr);
}
/***********************************************************************

check_directory  -  Check whether the given directory exists

***********************************************************************/

int check_directory(char *directory,int dos)
{
  char command[C_FILENAME_LEN];

  char *string_ptr;

  int error, exists;

  /* Initialise variables */
  error = FALSE;
  exists = FALSE;

  /* Check whether directory exists */
  if (dos == TRUE)
    strcpy(command,"dir /D ");
  else
    strcpy(command,"ls -1d ");
  strcat(command,directory);

  /* Append null redirect if output not already being redirected */
  if (strchr(command,'>') == NULL)
    {
      if (dos == TRUE)
	strcat(command," > NUL");
      else
	strcat(command," >& /dev/null");
    }
#ifdef NOSYST
  error = FALSE;
#else
  error = system(command);
#endif

  /* If directory doesn't exist, create it */
  if (error == FALSE)
    exists = TRUE;

  /* Return whether the directory exists */
  return(exists);
}
/***********************************************************************

replace_char  -  Replace all occurrences of one character by another

***********************************************************************/

void replace_char(char *string,char replace,char by)
{
  int i, len;

  /* Get the length of the string */
  len = strlen(string);

  /* Loop through the characters in the string to perform replacements */
  for (i = 0; i < len; i++)
    {
      /* If current character is to be replaced, perform replacement */
      if (string[i] == replace)
        string[i] = by;
    }
}
/***********************************************************************

call_system  -  Run the given system command

***********************************************************************/

int call_system(char *cmd,int dos)
{
  char *command;
  
  int error;
  int len;

  /* Initialise flag */
  error = FALSE;

  /* If Dos, change all the directory separators */
  if (dos == TRUE && strncmp(cmd,"timeout",7))
    replace_char(cmd,'/','\\');

  /* Get the length of the command */
  len = strlen(cmd);

  /* Allocate memory for the string */
  command = (char *) malloc(sizeof(char)*(len + 10));
  strcpy(command,cmd);

  /* For DOS, append the NUL redirect */
  if (dos == TRUE && strchr(command,'>') == NULL)
    strcat(command," > NUL");

/* debug 
  printf("COMMAND: %s\n",command);
  fflush(stdout);
  printf("  error = %d\n",error);
  fflush(stdout);
 debug */

  /* Run the command */
  error = system(command);

  return(error);
}
/***********************************************************************

create_directory  -  Create the given directory

***********************************************************************/

int create_directory(char *directory)
{
  char command[C_FILENAME_LEN], dir_name[C_FILENAME_LEN];

  char *string_ptr;

  int error, missing, dos;

  struct stat stats;
   
  /* Initialise variables */
  error = FALSE;
  dos = FALSE;

  /* See if this looks like a DOS directory */
  string_ptr = strstr(directory,":/");
  if (string_ptr != NULL)
    dos = TRUE;
  
  /* Check whether the directory exists */
  missing = stat(directory,&stats);

  /* If directory doesn't exist, create it */
  if (missing != FALSE)
    {
      /* Get directory name */
      strcpy(dir_name,directory);

      /* For DOS, change all the slashes to backslashes */
      if (dos == TRUE)
        replace_char(dir_name,'/','\\');

      /* Create the directory */
#ifdef WIN
      error = mkdir(dir_name);
#else
      error = mkdir(dir_name, S_IRWXU | S_IRWXG | S_IRWXO);
#endif

      /* Show warning if create failed */
      if (error != FALSE)
        printf("*** Warning. Unable to create directory: %s\n",
               directory);
    }

  /* Return error flag */
  return(error);
}
/***********************************************************************

create_directory_old  -  Create the given directory

***********************************************************************/

int create_directory_old(char *directory)
{
  char command[C_FILENAME_LEN];

  char *string_ptr;

  int error, missing, dos;

  struct stat stats;
  
  /* Initialise variables */
  error = FALSE;
  dos = FALSE;

  /* See if this looks like a DOS directory */
  string_ptr = strstr(directory,":/");
  if (string_ptr != NULL)
    dos = TRUE;
  
  /* Check whether the directory exists */
  missing = stat(directory,&stats);

  /* If directory doesn't exist, create it */
  if (missing != FALSE)
    {
      /* Create the directory */
      strcpy(command,"mkdir ");
      strcat(command,directory);
      if (dos == TRUE)
        {
          /* For DOS, change all the slashes to backslashes */
          replace_char(command,'/','\\');
          strcat(command," > NUL");
        }
      else
        strcat(command," >& /dev/null");
#ifdef NOSYST
      error = FALSE;
#else
      error = system(command);
#endif
      if (error != FALSE)
        printf("*** Warning. Unable to create directory: %s\n",
               directory);

      /* Change directory premission */
      if (dos == FALSE)
        {
          strcpy(command,"chmod 777 ");
          strcat(command,directory);
#ifdef TEST
          error = FALSE;
#else
          system(command);
#endif
        }
    }

  /* Return error flag */
  return(error);
}
/***********************************************************************

file_months  -  Get how old the given file is in months

***********************************************************************/

int file_months(char *file_name)
{
#define FACTOR     (60 * 60 * 24 * 30)

  char date_string[64];

  int current_time, diff, nmonths, file_time;

  struct stat st;

  struct tm *current_date;

  size_t i;

  time_t now = time(NULL);

  /* Get the time file was last modified */
  if (stat(file_name,&st) != 0)
    {
      /* If file not found, return -1 */
      return(-1);
    }

  /* Convert time to seconds */
//  sprintf(date_string,"%i\n", st.st_mtime);
//  file_time = atoi(date_string);
  sprintf(date_string,"%ld\n", st.st_mtime);
  file_time = atol(date_string);

  /* Get the current date */
  current_date = localtime(&now);
  current_time = now;

  /* Get the number of seconds since the file was modified */
  diff = current_time - file_time;

  /* Convert to months */
  nmonths = diff / FACTOR;

  /* Return number of whole months since the file was modified */
  return(nmonths);
}
/***********************************************************************

fix_upper  -  Search for certain strings that must be set in upper-case

***********************************************************************/

void fix_upper(char *name)
{
#define NCAPWORDS 30

  char ch;

  char *string_ptr;

  int ipos, jpos, len, nchars, word;
  int check_word[NCAPWORDS], done, is_word;

  static char *lowword[] = {
    "adp", "amp", "atp", "coa", "dna", "ec", "e.c", "fab", "fad", "hiv",
    "mrna", "nad", "nadp", "nadph", "nag", "pdb", "rna",  "trna",
    "i", "iii", "iii", "iv", "v", "vi", "vii", "viii", "ix", "x",
    "xi", "xii"
  };

  static char *capword[] = {
    "ADP", "AMP", "ATP", "CoA", "DNA", "EC", "E.C", "FAB", "FAD", "HIV",
    "mRNA", "NAD", "NADP", "NADPH", "NAG", "PDB", "RNA", "tRNA",
    "I", "III", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X",
    "XI", "XII"
  };

  /* Initialise variables */
  done = FALSE;
  nchars = strlen(name);
  for (word = 0; word < NCAPWORDS; word++)
    check_word[word] = TRUE;

  /* Loop until all cap-words have been replaced */
  while (done == FALSE)
    {
      /* Assume that we are done */
      done = TRUE;

      /* Loop through all the key words */
      for (word = 0; word < NCAPWORDS; word++)
        {
          /* If this word to be checked, then do so */
          if (check_word[word] == TRUE)
            {
              /* Search for presence of this word in the string */
              string_ptr = strstr(name,lowword[word]);
              if (string_ptr != NULL)
                {
                  /* Assume this is a complete word */
                  is_word = TRUE;

                  /* Get the string's start-position */
                  ipos = string_ptr - name;
                  len = strlen(lowword[word]);
                  jpos = len + ipos;

                  /* Check character to left of match to test for word
                     separator */
                  if (ipos > 0)
                    {
                      /* Get character to left */
                      ch = name[ipos - 1];

                      /* If character is a letter, then match does
                         not represent a word */
                      if ((ch >= 'a' && ch <= 'z') ||
                          (ch >= 'A' && ch <= 'Z'))
                        is_word = FALSE;
                    }

                  /* If still possibly a word, check the character
                     to the right of the word */
                  if (is_word == TRUE && jpos < nchars)
                    {
                      /* Get character to right */
                      ch = name[jpos];

                      /* If character is a letter, then match does
                         not represent a word */
                      if ((ch >= 'a' && ch <= 'z') ||
                          (ch >= 'A' && ch <= 'Z'))
                        is_word = FALSE;
                    }

                  /* If have a capital-word, then splice into string */
                  if (is_word == TRUE)
                    {
                      /* Splice into string */
                      strncpy(string_ptr,capword[word],len);

                      /* Set flag that we're not done yet */
                      done = FALSE;
                    }
                }

              /* If word not present, then don't need to check for
                 it again */
              else
                check_word[word] = FALSE;
            }
        }
    }
}
/***********************************************************************

form_filename  -  Use the filename template and entered PDB code to
                  determine full name of file

***********************************************************************/

char *form_filename(char *pdb_code,char *filename_template)
{
  char ch, copy_char[2], file_name[C_FILENAME_LEN];

  char *final_name = &null;
  char *string_ptr;

  int i, ichar, len, len_code;
  int incode;

  /* Initialise variables */
  file_name[0] = '\0';
  incode = FALSE;

  /* Get length of template and PDB code */
  len = strlen(filename_template);
  len_code = strlen(pdb_code);

  /* If not too long, form file name */
  if (len < C_FILENAME_LEN)
    {
      /* Find any strings between square brackets in this template */
      for (i = 0; i < len; i++)
        {
          /* Get this character */
          ch = filename_template[i];
          copy_char[0]= '\0';

          /* If this an open-bracket, set flag */
          if (ch == '[')
            incode = TRUE;

          /* If close-bracket, unset flag */
          else if (ch == ']')
            incode = FALSE;
          else
            {
              /* If in code, transfer appropriate part of PDB code to
                 file name */
              if (incode == TRUE)
                {
                  /* See if this character is a template character */
                  string_ptr = strchr(TEMPLATE_CHARS,ch);

                  /* If present, transfer corresponding character from
                     PDB code */
                  if (string_ptr != NULL)
                    {
                      /* Get this character's array position */
                      ichar = string_ptr - TEMPLATE_CHARS;

                      /* If valid, then save character as the copy
                         character */
                      if (ichar > -1 && ichar < len_code)
                        copy_char[0] = pdb_code[ichar];
                    }
                }

              /* Otherwise, just transfer character from template */
              else
                copy_char[0] = ch;

              /* If have a character to copy across, perform the transfer */
              if (copy_char[0] != '\0')
                {
                  /* Append null and tack onto end of file name */
                  copy_char[1] = '\0';
                  strcat(file_name,copy_char);
                }
            }
        }

      /* Store the final file name */
      store_string(&final_name,file_name);
    }

  /* Otherwise, show error message */
  else
    printf("*** ERROR. Filename template too long!\n"
           "           Template: %s\n",filename_template);

  /* Return name of file */
  return(final_name);
}
/***********************************************************************

find_pdbsum_dir  -  Locate the PDBsum directory for this PDB entry

***********************************************************************/

int find_pdbsum_dir(char *pdb_code,char *pdbsum_template[NPDBSUM_DIR],
                    int pdb_type,char *dir_name)
{
  char file_name[C_FILENAME_LEN];
  char *fname;

  int dir_type, end, len, start, type;
  int error;

  struct stat st;

  /* Initialise */
  strcpy(dir_name,"./");
  dir_type = NOT_FOUND;

  /* Determine the directory type */
  if (pdb_type == STANDARD_PDB || pdb_type == ALPHA_FOLD_PDB ||
      pdb_type == PDBSUM1)
    {
      start = PDBSUM_CURRENT;
      end = PDBSUM_CURRENT + 1;
    }
  else if (pdb_type == FIND_FIRST)
    {
      start = PDBSUM_CURRENT;
      end = NPDBSUM_DIR;
    }
  else
    {
      start = PDBSUM_UPLOAD;
      end = NPDBSUM_DIR;
    }

  /* Loop over the PDBsum directories to try */
  for (type = start; type < end && dir_type == NOT_FOUND; type++)
    {
      /* Form the name of the PDBsum directory */
      fname = form_filename(pdb_code,pdbsum_template[type]);
      strcpy(file_name,fname);

      /* If searching for PDBsum dir, see if this one exists */
      if (pdb_type == FIND_FIRST)
	error = stat(file_name,&st);
      else
	error = FALSE;

      /* If directory exists, save its name */
      if (error == FALSE)
        {
          /* Save the directory name */
          strcpy(dir_name,fname);

          /* Set flag */
          dir_type = type;
        }
    }

  /* Append directory separator if required */
  if (dir_type != NOT_FOUND)
    {
      len = strlen(dir_name);
      if (dir_name[len - 1] != '/')
        strcat(dir_name,"/");
    }

  /* Return the file type */
  return(dir_type);
}
/***********************************************************************

format_exp_no  -  Convert exponential-format number to HTML notation

***********************************************************************/

int format_exp_no(char *number_string,char *html_string)
{
  char tmp[50];

  char *string_ptr;

  int band_no, iband, exponent;

  double band_val, number;

  /* Get the number to determine which band this value lies in */
  number = atof(number_string);
  band_val = 1;
  band_no = -1;
  
  /* Loop until we find the corresponding band */
  for (iband = 0; iband < MIN_ALLELE_FREQ && band_no == -1; iband++)
    {
      /* If value above the current band value, then have band */
      if (number > band_val)
        band_no = iband;

      /* Reduce band value by a factor of 10 */
      band_val = band_val / 10;
    }

  /* If no band identified, then make zero */
  if (band_no == -1)
    band_no = 0;
  
  /* Get the exponent */
  exponent = 9;
  strcpy(html_string,number_string);
  string_ptr = strchr(html_string,'e');
  if (string_ptr == NULL)
    string_ptr = strchr(html_string,'E');
  if (string_ptr != NULL)
    {
      /* Read the exponent */
      exponent = atoi(string_ptr + 1);

      /* Terminate the string */
      string_ptr[0] = '\0';
    }

  /* Get the number part */
  number = atof(html_string);

  /* Form the HTML string */
  sprintf(html_string,"%3.1f",number);
  if (exponent != 9)
    {
      if (exponent != 0)
        {
          sprintf(tmp,"&times;10<SUP>%d</SUP>",exponent);
          strcat(html_string,tmp);
        }
    }
  else
    strcpy(html_string,number_string);

  /* Return the band number */
  return(band_no);
}
/***********************************************************************

free_tokens - Free up the memory used for storing token strings

**********************************************************************/  

void free_tokens(char **token,int ntokens)
{
  int itoken;

  /* Loop over the tokens to free up the memory */
  for (itoken = 0; itoken < ntokens; itoken++)
    if (token[itoken][0] != '\0')
      free(token[itoken]);
}
/***********************************************************************

get_aa_code  -  Get single-letter code from amino acid name

***********************************************************************/

char get_aa_code(char *res_name)
{
  char aa_code;
  int iamino, icode;
  int namino;

  static char amino[] = "ACDEFGHIKLMNPQRSTVWYXBZXMX";
  static char *amino_name[] = {
    "ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU",
    "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR",
    "PCA ","ASX ","GLX ", "UNK", "MSE", "XXX"
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

get_colour_rgb  -  Get the colour name for the given colour number and
                   the corresponding RGB values

***********************************************************************/

void get_colour_rgb(char *c_name,double rgb[3])
{
  int colour, icol;

  /* Initialise colour to first colour */
  colour = 0;

  /* Search colours list for given colour name */
  for (icol = 0; icol < NCOLOURS; icol++)
    if (!strcmp(c_name,COLOUR_LIST[icol]))
      colour = icol;

  /* Get the name and RGB triplet for this colour */
  rgb[0] = RGB_TRIPLET[colour][0];
  rgb[1] = RGB_TRIPLET[colour][1];
  rgb[2] = RGB_TRIPLET[colour][2];
}
/***********************************************************************

get_chain_colour  -  Get the given chain's colour

***********************************************************************/

int get_chain_colour(char chain,double *rgb,double *rgb_light)
{
  int colour, i, ichain, ipos, nchains;

  /* Determine number of available chain identifiers */
  nchains = strlen(CHAIN_LIST);

  /* Identify the chain-id */
  ichain = -1;
  for (ipos = 0; ipos < nchains && ichain == -1; ipos++)
    if (chain == CHAIN_LIST[ipos])
      ichain = ipos;

  /* If chain not found, use default */
  if (ichain == -1)
    ichain = 0;

  /* Get the sphere colour for this chain */
  colour = ichain % NCHAIN_COLOURS;

  /* Get rgb value for this colour */
  get_colour_rgb(CHAIN_COLOUR[colour],rgb);

  /* Calculate "light" version of this colour */
  for (i = 0; i < 3; i++)
    {
      /* Compute lighter colour */
      rgb_light[i] = rgb[i] + (1.0 - rgb[i]) / 2.0;
      if (rgb_light[i] > 1.0)
        rgb_light[i] = 1.0;
    }

  /* Return the colour number */
  return(colour);
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

get_date  -  Get today's date from the system clock

***********************************************************************/

void get_date(char *date)
{
  char date_string[64];

  struct tm *current_date;
 
  time_t now = time(NULL);

  /* Set the current date */
  current_date = localtime(&now);

  /* Format the date string */
  strftime(date_string,sizeof date_string,"%d-%b-%y",current_date);

  /* Convert to upper case for use in PDB header */
  to_upper(date_string,strlen(date_string));

  /* Save the date string */
  if (strlen(date_string) == 9)
    strcpy(date,date_string);
  else
    {
      /* If no leading zero on date, then add one */
      strcpy(date,"0");
      strcat(date,date_string);
    }
}
/***********************************************************************

get_date_and_time  -  Get today's date from the system clock

***********************************************************************/

void get_date_and_time(char *date)
{
  char date_string[64];

  struct tm *current_date;
 
  time_t now = time(NULL);

  /* Initialise */
  date[0] = date_string[0] = '\0';
  
  /* Set the current date */
  current_date = localtime(&now);

  /* Format the date string */
  strftime(date_string,sizeof date_string,"%H:%M:%S  %a %d %b %C%y",
           current_date);

  /* If format has failed, use sprintf instead */
  if (date_string[0] == '\0')
    sprintf(date_string,"%d:%d:%d  %d %s %d",
            (*current_date).tm_hour,
            (*current_date).tm_min,
            (*current_date).tm_sec,
            (*current_date).tm_mday,
            MONTH[(*current_date).tm_mon],
            1900 + (*current_date).tm_year);

  /* Save the date string */
  strcpy(date,date_string);
}
/***********************************************************************

get_res_name  -  Get the three-letter residue name from the single-letter
                 code

***********************************************************************/

int get_res_name(char aa_code,char *res_name)
{
  char *string_ptr;

  int ires;

  static char amino[] = "ACDEFGHIKLMNPQRSTVWYXBZXXM*#%-+";
  static char *amino_name[] = {
    "Ala", "Cys", "Asp", "Glu", "Phe", "Gly", "His", "Ile", "Lys", "Leu",
    "Met", "Asn", "Pro", "Gln", "Arg", "Ser", "Thr", "Val", "Trp", "Tyr",
    "Ter", "Asx", "Glx", "Unk", "Ter", "Mse", "Ter", "Seg", "Mis", "Del",
    "Ins"
  };

  /* Initialise residue number */
  ires = 30;
  string_ptr = NULL;

  /* Find this code in the list */
  if (aa_code != '\0')
    string_ptr = strchr(amino,aa_code);

  /* If found, return the corresponding residue name */
  if (string_ptr != NULL)
    {
      /* Get the corresponding residue position */
      ires = string_ptr - amino;

      /* Return the corresponding residue name */
      strcpy(res_name,amino_name[ires]);
    }

  /* Otherwise, return unknown residue name */
  else
    strcpy(res_name,"XXX");

  /* Return the residue number */
  return(ires);
}
/***********************************************************************

get_rgb  -  Get the colour'S RGB values

***********************************************************************/

int get_rgb(int colour,int *irgb)
{
  int i;
  int is_grey;

  double rgb[3], maxrgb, minrgb, span;

  /* Initialise grey flag */
  is_grey = FALSE;

  /* Get the RGB triplet for this colour */
  rgb[0] = colour / (256 * 256);
  colour = colour % (256 * 256);
  rgb[1] = colour / 256;
  rgb[2] = colour % 256;

  /* Get the maximum and minimum rgb values */
  minrgb = maxrgb = rgb[0];
  for (i = 1; i < 3; i++)
    {
      if (rgb[i] < minrgb)
        minrgb = rgb[i];
      if (rgb[i] > maxrgb)
        maxrgb = rgb[i];
    }

  /* Get the span of rgb values */
  span = maxrgb - minrgb;

  /* If difference is small, then colour is some shade of grey */
  if (span < SPAN_CUTOFF)
    is_grey = TRUE;

  /* Convert into int from double format */
  for (i = 0; i < 3; i++)
    irgb[i] = (int) rgb[i];

  /* Return whether this colour is a shade of grey */
  return(is_grey);
}
/***********************************************************************

get_irgb_colour  -  Get the integer componants for the given colour name

***********************************************************************/

void get_irgb_colour(char *c_name,int *irgb)
{
  int i;

  double rgb[3];

  /* Get the corresponding floating-point RGB components */
  get_colour_rgb(c_name,rgb);

  /* Convert to integer equivalents (ie 0-255) */
  for (i = 0; i < 3; i++)
    irgb[i] = (int) (255 * rgb[i] + 0.5);
}
/***********************************************************************

get_pdbsum1_dir  -  Get the PDBsum1 directory for the given entry

***********************************************************************/

void get_pdbsum1_dir(char *pdir,char *pdb_code,
                     struct c_file_data file_name_ptr)
{
  char parent_dir[C_FILENAME_LEN], sub_dir[3];
  char ch, new_name[C_FILENAME_LEN], term[5];

  char *string_ptr;

  int i, ipos, len;
  int gzipped, ok, zipped;

  /* Get the 2nd and 3rd digits of the PDB code */
  len = strlen(pdb_code);
  strncpy(sub_dir,pdb_code + 1,2);
  sub_dir[2] = '\0';

  /* Form the name of the parent directory */
  strcpy(parent_dir,file_name_ptr.other_dir[PDBSUM1_RESULTS_DIR]);
  strcat(parent_dir,"/");
  strcat(parent_dir,sub_dir);

  /* Form the name of the PDBsum directory */
  strcpy(pdir,parent_dir);
  strcat(pdir,"/");
  strcat(pdir,pdb_code);
}
/***********************************************************************

get_random_number  -  Routine for generating a random number from the
                      supplied seed. Returns a uniform random deviate
                      between 0.0 and 1.0. On the first call, the
                      seed must be set to any negative value to
                      initialise or reinitialise the sequence.
                      (Algorithm is that described for routine ran2
                      in Numerical Recipes).

***********************************************************************/

float get_random_number(int *iseed,int random_start)
{
  static int ir[97], iy, j;

  static float random_number, seed;

  static int first_pass = TRUE;
  static int m = 714025;
  static int ia = 1366;
  static int ic = 150889;
  static float rm = (float) (1.4005112E-6);

  FILE *file_ptr;

  /* Initialise the random number */
  random_number = 0.0;

  /* If this is the first pass, or sequence to be re-initialised, then
     pick up seed from file */
  if (*iseed < 0 || first_pass == TRUE)
    {
      /* Initialise the seed */
      seed = 0.0;

      /* If a completely random start is required, then read in the
         number written out to the file ranseed.dat last time this
         routine was run */
      if (random_start == TRUE)
        {
          /* Open the random seed file and read in seed */
          if ((file_ptr = fopen("ranseed.dat","r")) !=  NULL)
            {
              /* If file exists, then read in the seed */
              fscanf(file_ptr,"%f",&seed);

              /* Close the file */
              fclose(file_ptr);
            }

          /* Store integer version of the seed */
          *iseed = (int) seed;
        }

      /* Otherwise, use a standard starting point */
      else
        *iseed = 54813;

      /* Check that the seed is non-zero */
      if (*iseed == 0) 
        *iseed = 1;

      /* Generate random number from seed */
      *iseed = (ic - *iseed) % m;
      for (j = 0; j < 97; j++)
        {
          *iseed = (ia * (*iseed) + ic) % m;
          ir[j] = *iseed;
        }
      *iseed = (ia * (*iseed) + ic) % m;
      iy = *iseed;
    }

  /* Check for error */
  j = (97 * iy) / m;
  if (j >= 97 || j < 0)
    {
      printf("*** Random number error\n");
      exit(1);
    }

  /* Get the random number */
  iy = ir[j];
  random_number = iy * rm;
  *iseed = (ia * (*iseed) + ic) % m;
  ir[j] = *iseed;

  /* If this is the first call, use the generated random number to
     create a new seed and write to disk */
  if (first_pass == TRUE)
    {
      first_pass = FALSE;
      seed = iy * rm * 100000;

      /* Open the random seed file */
      if ((file_ptr = fopen("ranseed.dat","w")) !=  NULL)
        {
          /* If no file error, then write out the seed */
          fprintf(file_ptr,"%f\n",seed);

          /* Close the file */
          fclose(file_ptr);
        }
    }

  /* Return the random number generated */
  return(random_number);
}
/***********************************************************************

get_pixel_shade  -  Convert letter-pixel to appropriate shade

***********************************************************************/

float get_pixel_shade(char ch)
{
  int shade;

  float grey_shade;

  /* Initialise shade */
  grey_shade = 0.0;

  /* If this character is an 'X', set its grey-value to zero */
  if (ch == 'X')
    grey_shade = 1.0;

  /* If a shade of grey, adjust the colour shade */
  else if (ch != ' ')
    {
      /* Get the corresponding shade */
      shade = ch - 'A';
      if (shade < 0)
        shade = 0;
      else if (shade > NSHADES - 1)
        shade = NSHADES - 1;

      /* Convert shade to grey value */
      grey_shade = 1.0 - (float) shade / (float) NSHADES;
    }

  /* Return the grey shade */
  return(grey_shade);
}
/***********************************************************************

get_tokens - Get the space- or tab-delimited token from the given
             input line

**********************************************************************/  

int get_tokens(char *line,char **token,int max_token,char token_sep)
{
  char ch, last_ch, token_sep2;

  char *tmp, *tmp_cpy;

  int intok, ipos, len, ntoken, tokenlen;
  int done, keep_token;

  /* Initialise variables */
  done = FALSE;
  keep_token = FALSE;
  last_ch = '\0';
  ntoken = 0;
  intok = 0;
  ipos = 0;
  if (token_sep == '[')
    {
      token_sep2 = ']';
      keep_token = TRUE;
    }
  else if (token_sep == ':')
    token_sep2 = ':';
  else
    token_sep2 = '\0';

  /* Get the length of the line */
  len = strlen(line);
  if (len == 0)
    return(ntoken);

  /* Allocate space for temporary string */
  tmp = (char *) malloc(sizeof(char)*(len + 1));
  tmp[0] = '\0';
  tmp_cpy = (char *) malloc(sizeof(char)*(len + 1));
  tmp_cpy[0] = '\0';

  /* Loop until end of line reached */
  while (done == FALSE && ipos < len + 1 && ntoken < max_token)
    {
      /* Get the next character in the string */
      ch = line[ipos];

      /* Check for end of line */
      if (ch == '\0' || ch == '\n')
        {
          done = TRUE;
          ch = token_sep;
        }

      /* Check for double-colon separator */
      if (token_sep == ':' && ch == token_sep && last_ch == token_sep)
        ch = '\0';

      /* If have a token separator, then take to be the end of the
         token if have already got something in it */
      if (ch == token_sep || ch == token_sep2)
        {
          /* If have a token, increment token count */
          if (tmp[0] != '\0' || token_sep != ' ')
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
                  strcpy(tmp_cpy,tmp+1);
                  strcpy(tmp,tmp_cpy);
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
             start next token with it */
          else if (keep_token == TRUE && ch == token_sep)
            {
              tmp[intok] = ch;
              intok++;
            }
        }

      /* Otherwise, add character to token-string provided not
         too long */
      else if (ch != '\0')
        {
          if (ntoken < max_token && intok < len)
            {
              tmp[intok] = ch;
              intok++;
            }
        }

      /* Save the current character */
      last_ch = ch;

      /* Go to the next character in the line */
      ipos++;
    }

  /* Free up memory taken by temporary line */
  free(tmp);
  free(tmp_cpy);

  /* Return the number of tokens picked up */
  return(ntoken);
} 
/***********************************************************************

is_number  -  Determine whether the given string contains a number only

***********************************************************************/

int is_number(char *string)
{
  char ch;

  int ipos;
  int numeric;

  /* Initialise variables */
  numeric = TRUE;

  /* Loop through all the character positions, looking for any non-
     numeric characters */
  for (ipos = 0; ipos < strlen(string) && numeric == TRUE; ipos++)
    {
      /* Check the current character */
      ch = string[ipos];
      if ((ch >= '0' && ch <= '9') || ch == '.' || ch == '+' ||
          ch == '-')
        numeric = TRUE;
      else
        numeric = FALSE;
    }

  /* Return the numeric flag */
  return(numeric);
}
/***********************************************************************

is_numeric  -  Determine whether the given string contains any numerals

***********************************************************************/

int is_numeric(char *string)
{
  char ch;

  int ipos;
  int numeric;

  /* Initialise variables */
  numeric = FALSE;

  /* Loop through all the character positions, checking whether any
     numbers are present */
  for (ipos = 0; ipos < strlen(string) && numeric == FALSE; ipos++)
    {
      /* Check the current character */
      ch = string[ipos];
      if (ch >= '0' && ch <= '9')
        numeric = TRUE;
    }

  /* Return the numeric flag */
  return(numeric);
}
/***********************************************************************

old_move_file  -  Move the given file to new name (not used)

***********************************************************************/

void old_move_file(char *file_name1,char *file_name2)
{
  char *command;

  int len;
  int changed;

  /* Get the length of the command line */
  len = strlen(file_name1) + strlen(file_name2) + 10;

  /* Allocate memory for the command */
  command = (char *) malloc(sizeof(char)*(len + 1));

  /* Form the move command */
  sprintf(command,"mv %s %s > /dev/null",file_name1,file_name2);

  /* Run the command */
#ifdef TEST
#else
  system(command);
#endif

  /* Free up memory taken by the command */
  free(command);
}
/***********************************************************************

output_ppm_file  -  Write out the given image array as a .ppm file

***********************************************************************/

void output_ppm_file(int *array,int ncols,int nrows,char *file_name,
                     char *plot_name)
{
  int colour, ntrip, rgb[3];

  long iloc, size;

  FILE *file_ptr;

  /* Get the image size */
  size = ncols * nrows;

  /* Open the output file */
  file_ptr = open_output_file(file_name,FALSE);
  if (file_ptr == NULL)
    return;

  /* Write out header records */
  fprintf(file_ptr,"P3\n");
  fprintf(file_ptr,"# %s %s\n",plot_name,file_name);
  fprintf(file_ptr,"%d %d\n",ncols,nrows);
  fprintf(file_ptr,"255\n");

  /* Initialise count of RGB triplets on line */
  ntrip = 0;

  /* Loop through the pixels making up the image */
  for (iloc = 0; iloc < size; iloc++)
    {
      /* Get the colour at this pixel position */
      colour = array[iloc];

      /* Convert colour into an RGB triplet */
      rgb[0] = colour / (256 * 256);
      colour = colour % (256 * 256);
      rgb[1] = colour / 256;
      rgb[2] = colour % 256;

      /* Write out this pixel */
      fprintf(file_ptr,"%d %d %d ",rgb[0],rgb[1],rgb[2]);

      /* Increment count of triplets on this line */
      ntrip++;

      /* If have too many triplets on this line, throw
         new line */
      if (ntrip > 4)
        {
          fprintf(file_ptr,"\n");
          ntrip = 0;
        }
    }

  /* Close the .ppm image file */
  fclose(file_ptr);
}
/***********************************************************************

parse_uniprot_name  -  Parse the UniProt title to split up into its
                       components: protein name, species, etc

***********************************************************************/

void parse_uniprot_name(char *string,char *protein,char *species,
                        char *uniprot_acc,char *uniprot_id)
{
  char *tmp, *string_ptr, *start_string_ptr;

  int i, last_space, len;
  int done;

  /* Initialise variables */
  strcpy(protein,string);
  species[0] = '\0';
  uniprot_acc[0] = '\0';
  uniprot_id[0] = '\0';

  /* Get the length of the string */
  len = strlen(string);

  /* Allocate memory for temporary string */
  tmp = (char *) malloc(sizeof(char)*(len + 1));

  /* Identify the start of the codes */
  start_string_ptr = string_ptr = strstr(string,"sp|");
  if (string_ptr == NULL)
    start_string_ptr = string_ptr = strstr(string,"tr|");

  /* If found, then extract them */
  if (string_ptr != NULL)
    {
      /* Copy to temporary string */
      strcpy(tmp,string_ptr + 3);

      /* Get the position of the divider */
      string_ptr = strchr(tmp,'|');

      /* If found, we have the two codes, so can store */
      if (string_ptr != NULL)
        {
          /* Get the length of the accession code */
          len = string_ptr - tmp;

          /* Copy across */
          strncpy(uniprot_acc,tmp,len);
          uniprot_acc[len] = '\0';

          /* Get the UniProt id */
          string_ptr = strchr(tmp,' ');

          /* If end of id found, store the id */
          if (string_ptr != NULL)
            {
              /* Terminate the string */
              string_ptr[0] = '\0';

              /* Copy across */
              strcpy(uniprot_id,tmp + len + 1);
            }
        }

      /* Look for first blank after codes block */
      string_ptr = strchr(start_string_ptr,' ');

      /* If found, then transfer */
      if (string_ptr != NULL)
        {
          /* Save the start of the name string */
          start_string_ptr = string_ptr + 1;

          /* Copy protein name to temporary string */
          strcpy(tmp,start_string_ptr);

          /* Find the start of the species */
          string_ptr = strstr(tmp," OS=");

          /* If found, can save the protein name */
          if (string_ptr != NULL)
            string_ptr[0] = '\0';

          /* Save the protein name */
          strcpy(protein,tmp);
        }
    }

  /* Extract the species name */
  string_ptr = strstr(string," OS=");

  /* If found, then extract species name */
  if (string_ptr != NULL)
    {
      /* Copy to temporary area */
      strcpy(tmp,string_ptr + 4);

      /* Get the length of the current string */
      len = strlen(tmp);

      /* Loop through the characters to find the last space before the
         next item */
      last_space = 0;
      done = FALSE;
      for (i = 0; i < len && done == FALSE; i++)
        {
          /* Check if this is a space or a new field */
          if (tmp[i] == ' ')
            last_space = i;
          else if (tmp[i] == '=')
            done = TRUE;
        }

      /* Terminate the species string and copy */
      tmp[last_space] = '\0';
      strcpy(species,tmp);
    }

  /* Free up memory taken by the temporary string */
  free(tmp);
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

tokenize - Split the given string into tokens, replacing each instance
           of the token separator with a null and returning an array of
           string start-points for each token

           Note: this corrupts the original string

**********************************************************************/  

int tokenize(char *string,char **token,int max_tokens,char token_sep)
{
  char *string_ptr;

  int ntokens;
  int done;

  /* Initialise variables */
  done = FALSE;

  /* Set first token at rhe start of the string */
  token[0] = string;
  ntokens = 1;

  /* Loop until all tokens have been extracted */
  while (done == FALSE && ntokens < max_tokens)
    {
      /* Get location of the next token separator */
      string_ptr = strchr(string,token_sep);

      /* If have a token, null terminate and keep looking */
      if (string_ptr != NULL)
        {
          /* Shift start-pos of string and next token */
          token[ntokens] = string = string_ptr + 1;
          ntokens++;

          /* Null terminate current token */
          string_ptr[0] = '\0';
        }

      /* Otherwise, we're done */
      else
        done = TRUE;
    }

  /* Return number of tokens */
  return(ntokens);
}
/***********************************************************************

string_chomp  -  Replace line-feed by null

***********************************************************************/

int string_chomp(char *string)
{
  char cr = 13;

  char *string_ptr;

  /* Initialise variables */
  string_ptr = strchr(string,'\n');
  if (string_ptr != NULL)
    string_ptr[0] = '\0';
  else
    {
      /* Check for carriage return */
      string_ptr = strchr(string,cr);
      if (string_ptr != NULL)
        string_ptr[0] = '\0';
    }

  /* Return the string length */
  return(strlen(string));
}
/***********************************************************************

read_ppm_file  -  Read in the image in the given .ppm file

***********************************************************************/

int read_ppm_file(char *file_name,int **arr,int *nr,int *nc,
                  char *comment)
{
  char input_line[PPM_LINELEN + 1];

  char *token[PPM_MAX_TOKEN];

  int colour, ipos, irgb, itoken, len, line, ncols, nrows, ntokens;
  int rgb[3], size;
  int ok, skip;

  int *array;

  FILE *file_ptr;

  /* Initialise variables */
  comment[0] = '\0';
  ipos = 0;
  line = 0;
  ncols = nrows = size = 0;
  ok = FALSE;
  skip = FALSE;

  /* Open the file */
  if ((file_ptr = fopen(file_name,"r")) !=  NULL)
    {
      /* Loop through the file processing each record */
      while (fgets(input_line,PPM_LINELEN,file_ptr) != NULL)
        {
          /* Increment line-count */
          line++;

          /* Truncate the string to strip off the newline */
          len = string_chomp(input_line);

          /* If this is a comment line, store the comment */
          if (input_line[0] == '#')
            {
              /* Use if not blank */
              if (len > 1)
                {
                  /* If current comment is empty, then copy across */
                  if (comment[0] == '\0')
                    {
                      /* Skip initial blank if there is one */
                      if (input_line[1] == ' ')
                        strcpy(comment,input_line + 2);
                      else
                        strcpy(comment,input_line + 1);
                    }

                  /* Otherwise, add after tab if not too long */
                  else if (len + strlen(comment) < PPM_LINELEN - 1)
                    {
                      /* Add tab */
                      strcat(comment,TAB);

                      /* Skip initial blank if there is one */
                      if (input_line[1] == ' ')
                        strcat(comment,input_line + 2);
                      else
                        strcat(comment,input_line + 1);
                    }
                }
            }

          /* If line contains array size, read in parameters */
          else if (line > 1 && ncols == 0)
            {
              /* Get the array size */
              sscanf(input_line,"%d %d",&ncols,&nrows);
              size = ncols * nrows;

              /* Create array to hold image */
              array = (int *) malloc(size * sizeof(int));
              if (array == NULL)
                {
                  printf("*** ERROR. Unable to allocate memory for image "
                         "array: %d x %d in %s\n",ncols,nrows,file_name);
                  exit(1);
                }

              /* Initialise array */
              for (ipos = 0; ipos < size; ipos++)
                array[ipos] = 0;
              irgb = 0;
              ipos = 0;

              /* Set flag that OK */
              ok = TRUE;

              /* Set flag to skip next line */
              skip = TRUE;
            }

          /* If this is part of the image, retrieve the RGB numbers on
             this line  */
          else if (ncols > 0 && ipos < size)
            {
              /* Retrieve the numbers on this line */
//              ntokens = get_tokens(input_line,token,PPM_MAX_TOKEN,' ');
              ntokens = tokenize(input_line,token,PPM_MAX_TOKEN,' ');

              /* If skipping this line, set number of tokens to zero */
              if (skip == TRUE)
                ntokens = 0;
              skip = FALSE;

              /* Loop over the tokens to form RGB triplets */
              for (itoken = 0; itoken < ntokens; itoken++)
                {
                  /* Extract the colour number and store as part of
                     an RGB triplet */
                  if (token[itoken][0] != '\0')
                    {
                      rgb[irgb] = atoi(token[itoken]);
                      irgb++;
                    }

                  /* If have a triplet, then store the colour */
                  if (irgb == 3)
                    {
                      /* Retrieve the colour */
                      colour = 256 * 256 * rgb[0] + 256 * rgb[1] + rgb[2];

                      /* Store in the array */
                      array[ipos] = colour;
                      ipos++;

                      /* Reset the colour */
                      irgb = 0;
                    }
                }
            }
        }

      /* Close the file */
      fclose(file_ptr);
    }

  /* Show warning that couldn't find the given .ppm file */
  else
    {
      printf("*** Warning. Unable to find .ppm file: %s\n",file_name);
      ok = FALSE;
    }

  /* Return numbers of rows and columns */
  *nr = nrows;
  *nc = ncols;

  /* Return image array */
  *arr = array;

  /* Return whether file OK */
  return(ok);
}
/***********************************************************************

remove_blanks  -  Remove all blanks from the given string

***********************************************************************/

void remove_blanks(char *string)
{
  char *tmp;

  int i, len;
  int done;

  /* Initialise variables */
  done = FALSE;
  i = 0;

  /* Get the length of the string */
  len = strlen(string);

  /* Allocate memory for temporary string */
  tmp = (char *) malloc(sizeof(char)*(len + 1));

  /* Loop through the characters squeezing out the blanks */
  while (done == FALSE)
    {
      /* If current character is NULL, then we're done */
      if (string[i] == '\0')
        done = TRUE;

      /* If this is a blank, then shift string along */
      else if (string[i] == ' ')
        {
          strcpy(tmp,string + i + 1);
          strcpy(string + i,tmp);
        }

      /* Otherwise, increment position */
      else
        i++;
    }

  /* Free up memory taken by the temporary string */
  free(tmp);
}
/***********************************************************************

remove_char  -  Remove specified character from the given string

***********************************************************************/

void remove_char(char *string,char ch)
{
  char *tmp;

  int ipos, len;

  /* Get the length of the string */
  len = strlen(string);
  if (len == 0)
    return;

  /* Allocate memory for temporary string */
  tmp = (char *) malloc(sizeof(char)*(len + 1));

  /* Initialise variables */
  ipos = 0;

  /* Loop until string has been processed */
  while (ipos < len)
    {
      /* If have the given charaters at this position, remove */
      if (string[ipos] == ch)
        {
          /* Shift string along */
          strcpy(tmp,string + ipos + 1);
          strcpy(string + ipos,tmp);
          len--;
        }

      /* Otherwise, move on to next position */
      else
        ipos++;
    }

  /* Free up memory taken by the temporary string */
  free(tmp);
}
/***********************************************************************

remove_leading_blanks  -  Remove leading blank spaces from the given
                          string

***********************************************************************/

void remove_leading_blanks(char *string)
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

remove_marks  -  Remove all type-setting characters from the given
                 string

***********************************************************************/

void remove_marks(char *string)
{
  char *tmp;

  int ipos, len;

  /* Get the length of the string */
  len = strlen(string);
  if (len == 0)
    return;

  /* Allocate memory for temporary string */
  tmp = (char *) malloc(sizeof(char)*(len + 1));

  /* Initialise variables */
  ipos = 0;

  /* Loop until string has been processed */
  while (ipos < len)
    {
      /* If have a blank space at this position, remove */
      if (string[ipos] == '$' || string[ipos] == '/' ||
          string[ipos] == '\\')
        {
          /* Shift string along */
          strcpy(tmp,string+ipos+1);
          strcpy(string+ipos,tmp);
          len--;
        }

      /* Otherwise, move on to next position */
      else
        ipos++;
    }

  /* Free up memory taken by the temporary string */
  free(tmp);
}
/***********************************************************************

right_justify  -  Right justify the given string

***********************************************************************/

void right_justify(char *string,int length)
{
  int ipos, fpos, len;

  /* Get the length of the string */
  len = strlen(string);

  /* If string already filled, then return */
  if (len >= length - 1)
    return;

  /* Initialise final position */
  fpos = length - 2;

  /* Loop through the characters to shift to end of string */
  for (ipos = len - 1; ipos > -1; ipos--, fpos--)
    {
      /* Shift this character */
      string[fpos] = string[ipos];

      /* Blank out old character */
      string[ipos] = ' ';
    }

  /* Blank out any character at the start */
  for (ipos = 0; ipos < fpos + 1; ipos++)
    string[ipos] = ' ';

  /* Terminate the string */
  string[length - 1] = '\0';
}
/***********************************************************************

rgb_to_int  -  Convert RGB value to integers

***********************************************************************/

void rgb_to_int(double *rgb,int *irgb)
{
  int i;

  /* Loop over the three components */
  for (i = 0; i < 3; i++)
    irgb[i] = (int) (255 * rgb[i]);
}
/***********************************************************************

run_wget  -  Run wget to download appropriate page

***********************************************************************/

int run_wget(char *wget_exe,char *url,char *out_file,int dos)
{
  char command[C_FILENAME_LEN];

  int error, ok;

  /* Initialise variables */
  ok = TRUE;

  /* Form the lynx command to retrieve the page */
  strcpy(command,wget_exe);

  /* For DOS, change the slashes */
  if (dos == TRUE)
    replace_char(command,'/','\\');

  /* Append remainder of the parameters */
  strcat(command," -O ");
  strcat(command,out_file);
  strcat(command," -o wget.log --no-check-certificate ");
  strcat(command,url);

  /* Run wget to retrieve the HTML file */
  error = system(command);
  if (error != 0)
    ok = FALSE;

  /* Check if file downloaded */
  else
    ok = empty_check(out_file);

  /* Return whether download was successful */
  return(ok);
}
/***********************************************************************

strccnt  -  Count the number of occurrences of a given character in the
            given string

***********************************************************************/

int strccnt(const char *s,char c)
{
  int n = 0;

  while(*s != '\0') {
    if(*s++ == c)
      n++;
  }

  return n;
}
/***********************************************************************

string_chop  -  Truncate the given string at the last non-blank
                character

***********************************************************************/

int string_chop(char *string)
{
  int iend, ipos, len;

  /* Get length of string */
  len = strlen(string);

  /* If empty, then return */
  if (len == 0)
    return(len);

  /* Remove endline, if there is one */
  if (string[len - 1] == 13 || string[len - 1] == '\n')
    {
      string[len - 1] = '\0';
      len--;
    }

  /* Initialise end character */
  iend = 0;

  /* Loop back through the string until we hit the first non-blank
     character */
  for (ipos = len - 1; ipos > -1 && iend == 0; ipos--)
    {
      /* If not a space at this position, then have end of string */
      if (string[ipos] != ' ')
        iend = ipos + 1;
    }

  /* Truncate the string and set its length */
  string[iend] = '\0';
  len = iend;

  /* Return the string length */
  return(len);
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

type_set  -  Remove typographical characters and convert into upper and
             lower case.

***********************************************************************/

void type_set(char *name,int cap_start)
{
  char ch, last_ch;

  int ifree, ipos, nchars;
  int capital, dollar, skip, slash, upcase, up_one;

  /* Initialise variables */
  capital = TRUE;
  last_ch = ' ';
  upcase = FALSE;
  up_one = FALSE;
  ifree = 0;
  nchars = strlen(name);

  /* Check whether the string contains a dollar after a slash */
  dollar = FALSE;
  slash = FALSE;
  for (ipos = 0; ipos < nchars; ipos++)
    {
      ch = name[ipos];
      if (ch == '/')
        slash = TRUE;
      if (ch == '$' && slash == TRUE)
        dollar = TRUE;
    }

  /* Loop over all the characters */
  for (ipos = 0; ipos < nchars; ipos++)
    {
      /* Retrieve this character */
      ch = name[ipos];
      skip = FALSE;

      /* If in lower case, convert to upper */
      if (ch >= 'a' && ch <= 'z')
        ch = ch - 'a' + 'A';

      /* Check if this is a letter */
      if (ch >= 'A' && ch <= 'Z')
        {
          /* If first character then leave it in upper-case */
          if (capital == TRUE || up_one == TRUE)
            {
              capital = FALSE;
              up_one = FALSE;
            }

          /* Otherwise, set in lower-case, if required */
          else if (upcase == FALSE && up_one == FALSE)
            ch = ch - 'A' + 'a';
        }

      /* Check if this is a number */
      else if (ch >= '0' && ch <= '9')
        {
          /* If first character then set whole word in lower-case (as
             possibly a PDB code) */
          if (capital == TRUE)
            {
              capital = FALSE;
              up_one = FALSE;
            }
        }

      /* Check if this is one of the typesetting characters */

      /* Upper-case the next letter only */
      else if (ch == '*')
        {
          up_one = TRUE;
          skip = TRUE;
        }

      /* Upper-case until dollar */
      else if (ch == '/' && dollar == TRUE)
        {
          upcase = TRUE;
          skip = TRUE;
        }

      /* End upper-case */
      else if (ch == '$')
        {
          upcase = FALSE;
          skip = TRUE;
        }

      /* After a full-stop, need to capitalise next character */
      else if (ch == '.')
        capital = TRUE;

      /* If in author-list, capitalise next letter after a comma, dash
         or apostrophe too */
      else if (cap_start == TRUE && (ch == ',' || ch == '-' ||
                                     ch == '\''))
        capital = TRUE;

      /* Save the current character */
      if (skip == FALSE)
        {
          name[ifree] = ch;
          last_ch = ch;
          ifree = ifree + 1;
        }
    }

  /* Terminate string */
  name[ifree] = '\0';
}
/***********************************************************************

write_ppm_file  -  Write out the given image array as a .ppm file

***********************************************************************/

void write_ppm_file(int *array,int ncols,int nrows,char *file_name,
                    char *plot_name)
{
  int colour, icol, iloc, irow, ntrip, rgb[3];

  FILE *file_ptr;

  /* Open the output file */
  file_ptr = open_output_file(file_name,FALSE);
  if (file_ptr == NULL)
    return;

  /* Write out header records */
  fprintf(file_ptr,"P3\n");
  fprintf(file_ptr,"# %s %s\n",plot_name,file_name);
  fprintf(file_ptr,"%d %d\n",ncols,nrows);
  fprintf(file_ptr,"255\n");

  /* Initialise count of RGB triplets on line */
  ntrip = 0;

  /* Loop through the rows making up the image */
  for (irow = nrows - 1; irow > -1; irow--)
    {
      /* Loop over the columns to be plotted */
      for (icol = 0; icol < ncols; icol++)
        {
          /* Get position of this pixel in the image array */
          iloc = irow * ncols + icol;

          /* Get the colour at this pixel position */
          colour = array[iloc];

          /* Convert colour into an RGB triplet */
          rgb[0] = colour / (256 * 256);
          colour = colour % (256 * 256);
          rgb[1] = colour / 256;
          rgb[2] = colour % 256;

          /* Write out this pixel */
          fprintf(file_ptr,"%d %d %d ",rgb[0],rgb[1],rgb[2]);

          /* Increment count of triplets on this line */
          ntrip++;

          /* If have too many triplets on this line, throw
             new line */
          if (ntrip > 4)
            {
              fprintf(file_ptr,"\n");
              ntrip = 0;
            }
        }
    }

  /* Close the .ppm image file */
  fclose(file_ptr);
}
/***********************************************************************

convert_to_gif  -  Convert .ppm file to a GIF image

***********************************************************************/

void convert_to_gif(char *convert_exe,char *file_name,char *out_name,
                    int trim)
{
  char trim_option[10];
  
  char *command;

  int len;

  /* Select trim option of required */
  trim_option[0] = '\0';
  if (trim == TRUE)
    strcpy(trim_option,"-trim");

  /* Get the length of the command line */
  len = strlen(convert_exe) + strlen(file_name) + strlen(out_name) + 50;

  /* Show name of image file to be converted */
  printf("Converting .ppm image to %s ...\n",out_name);

  /* Allocate memory for the command */
  command = (char *) malloc(sizeof(char)*(len + 1));

  /* Form the parameters to run convert */
  sprintf(command,"%s %s %s %s  > convert.log",convert_exe,trim_option,
          file_name,out_name);

  /* Run convert to generate large GIF image */
  printf("   Running convert on %s ...\n",file_name);
  printf("     Command: %s\n",command);
 
  /* Run convert */
#ifdef TEST
#else
  system(command);
#endif

  /* Free up memory taken by the command */
  free(command);
}
/***********************************************************************

create_letter_record  -  Create a new letter record

***********************************************************************/

struct letter *create_letter_record(char character,int height,int width)
{
  int ipos, size;

  struct letter *letter_ptr;

  /* Initialise letter pointers */
  letter_ptr = NULL;

  /* Create letter entry */
  letter_ptr = (struct letter*) malloc(sizeof(struct letter));
  if (letter_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for");
      printf(" struct letter\n");
      printf("***        Program terminated with error.\n");
      exit (1);
    }

  /* Store the details for this entry */
  letter_ptr->character = character;
  letter_ptr->height = height;
  letter_ptr->width = width;

  /* Create array to hold this character's image definition */
  size = height * width;
  letter_ptr->array = (char *) malloc(size * sizeof(char));
  if (letter_ptr->array == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for letter array\n");
      exit(1);
    }

  /* Blank out array to initialise */
  for (ipos = 0; ipos < size; ipos++)
    letter_ptr->array[ipos] = ' ';

  /* Initialise remaining fields */
  letter_ptr->miny = letter_ptr->maxy = -1;

  /* Return pointer */
  return(letter_ptr);
}
/***********************************************************************

file_diff  -  Use the system diff command to determine whether files
              are identical or differ

***********************************************************************/

int file_diff(char *file_name1,char *file_name2)
{
  char *command;

  int len;
  int changed;

  /* Get the length of the command line */
  len = strlen(file_name1) + strlen(file_name2) + 10;

  /* Allocate memory for the command */
  command = (char *) malloc(sizeof(char)*(len + 1));

  /* Form the diff command */
  sprintf(command,"diff %s %s > /dev/null",file_name1,file_name2);

  /* Run diff */
#ifdef TEST
#else
  changed = system(command);
  if (changed != 0)
    changed = TRUE;
#endif

  /* Free up memory taken by the command */
  free(command);

  /* Return the result */
  return(changed);
}
/***********************************************************************

get_font  -  Pick up the character definitions from the given .ppm font
             file

***********************************************************************/

int get_font(struct font *font_ptr,int font,
             struct c_file_data file_name_ptr)
{
  char file_name[C_FILENAME_LEN], comment[C_LINELEN + 1];

  char *token[C_MAX_TOKEN];

  int icol, iend, iletter, iloc, ipos, irow, istart;
  int line, narray, ncols, nletters, nonblank, nrows, ntokens;
  int height, start_col, width, xpos, ypos;
  int colour, grey_value;
  int done, ok, skip_first;

  int *array;

  struct letter *letter_ptr, *first_letter_ptr, *last_letter_ptr;

  /* Initialise variables */
  line = 0;
  ncols = 0;
  iletter = -1;
  nletters = 0;

  /* Initialise pointers */
  letter_ptr = NULL;
  first_letter_ptr = NULL;
  last_letter_ptr = NULL;

  /* Form full name of wirletters.ppm file */
  strcpy(file_name,file_name_ptr.other_dir[PDBSUM_GENDATA_DIR]);
  strcat(file_name,DIR_SEP);
  strcat(file_name,FONT_FILE[font]);
  strcat(file_name,".ppm");

  /* Read in the image in the file */
  ok = read_ppm_file(file_name,&array,&nrows,&ncols,comment);
  narray = ncols * nrows;

  /* If not OK, then abort */
  if (ok == FALSE)
    exit(1);

  /* Extract the font name and character list from the comment lines */
  ntokens = get_tokens(comment,token,C_MAX_TOKEN,'\t');

  /* Store the font name */
  store_string(&(font_ptr->name),token[0]);

  /* Store the characters represented by the font definition */
  store_string(&(font_ptr->character),token[1]);

  /* Store the number of characters */
  font_ptr->ncharacters = strlen(font_ptr->character);

  /* Store the font height */
  font_ptr->height = nrows;

  /* Create index to store pointers to each defined character */
  font_ptr->letter_index_ptr
    = (struct letter **) malloc(font_ptr->ncharacters
                                * sizeof(struct letter *));

  /* Check that sufficient memory was allocated */
  if (font_ptr->letter_index_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for letter_index\n");
      exit(1);
    }

  /* Save the height of each character */
  height = nrows;

  /* Initialise flag and character count */
  done = FALSE;
  start_col = 0;

  /* Loop until all characters have been extracted */
  while (done == FALSE && nletters < font_ptr->ncharacters)
    {
      /* Initialise start and end positions of the next character */
      istart = iend = -1;
      skip_first = FALSE;
      if (font_ptr->character[nletters] == '"')
        skip_first = TRUE;

      /* Loop over the columns to find start and end of next character */
      for (icol = start_col; icol < ncols && (istart == -1 || iend == -1);
           icol++)
        {
          /* Initialise count of non-empty pixels */
          nonblank = 0;

          /* Loop over the image's rows to count the number of 
             non-empty pixels */
          for (irow = 0; irow < nrows; irow++)
            {
              /* Get the location in the character's array */
              iloc = irow * ncols + icol;
              
              /* If non-blank, increment count */
              if (array[iloc] != WHITE)
                nonblank++;
            }

          /* If line was blank and looking for end, then we have it */
          if (nonblank == 0 && istart != -1)
            {
              /* If looking for double-quotes, need to include first
                 blank */
              if (skip_first == TRUE)
                skip_first = FALSE;
              else
                iend = icol - 1;
            }

          /* If line non-blank and we're looking for start, then have
             it */
          else if (nonblank > 0 && istart == -1)
            istart = icol;
        }

      /* If have a character, then create a letter record to store its
         definition */
      if (istart != -1 && iend != -1 && iend >= istart)
        {
          /* Get the width of this character */
          width = iend - istart + 1;

          /* Create a letter record for this character */
          letter_ptr
            = create_letter_record(font_ptr->character[nletters],
                                   height,width);

          /* Store in the font record */
          font_ptr->letter_index_ptr[nletters] = letter_ptr;

          /* Initialise letter's x-position */
          xpos = 0;

          /* Loop over the relevant rows and columns of the image
             to extract the character definition */
          for (icol = istart; icol < iend + 1; icol++)
            {
              /* Re-initialise y-position */
              ypos = 0;

              /* Loop over the image's rows to count the number of 
                 non-empty pixels */
              for (irow = 0; irow < nrows; irow++)
                {
                  /* Get the location in the image array */
                  iloc = irow * ncols + icol;

                  /* Get the location in the character's array */
                  ipos = ypos * width + xpos;

                  /* If not white, then update the letter's y-limits */
                  if (array[iloc] != WHITE)
                    {
                      if (letter_ptr->miny == -1 || ypos < letter_ptr->miny)
                        letter_ptr->miny = ypos;
                      if (letter_ptr->maxy == -1 || ypos > letter_ptr->maxy)
                        letter_ptr->maxy = ypos;
                    }

                  /* If black, then store an X */
                  if (array[iloc] == BLACK)
                    letter_ptr->array[ipos] = 'X';

                  /* Otherwise, if grey, convert to letter A-T */
                  else if (array[iloc] != WHITE)
                    {
                      /* Get the grey level encoded in this number */
                      grey_value = array[iloc] / (256 * 256);
                      colour = (NSHADES * grey_value) / 255;
                      if (colour > NSHADES - 1)
                        colour = NSHADES - 1;

                      /* Save the corresponding letter */
                      letter_ptr->array[ipos] = LETTER[colour];
                    }

                  /* Increment y-position */
                  ypos++;
                }

              /* Increment x-position */
              xpos++;
            }

          /* Increment count of characters found */
          nletters++;
        }

      /* Otherwise, set flag that we're done */
      else
        done = TRUE;

      /* Set start-column after end of current character */
      start_col = iend + 1;
    }

  /* If number of characters read in doesn't agree with number in
     definition then have error */
  if (nletters != font_ptr->ncharacters)
    {
      printf("*** ERROR in font file %s\n",file_name);
      printf("***       Only %d of %d characters found\n",
             nletters,font_ptr->ncharacters);
      fflush(stdout);
      exit(1);
    }

  /* Return number of letters stored */
  return(nletters);
}
/***********************************************************************

extract_res_name  -  Extract the residue name and number

***********************************************************************/

int extract_res_name(char *residue,char *aac,char *aam,char **aa_ex,
                     char **mut_ex,int *typ,int *mut_len,char **add_var)
{
  char aa_code, aa_mut, ch, res_name[4], string[3][20];
  char residue_name[20];

  char *aa_extra, *added_variant, *mut_extra, *string_ptr;

  int i, ichar, ipos[3], istring, len, loop, mut_length, seq_no, type;
  int done, ok;

  /* Initialise variables */
  aa_code = aa_mut = '\0';
  aa_extra = mut_extra = &null;
  added_variant = &null;
  mut_length = 0;
  seq_no = 0;
  type = NSNP;
  istring = 0;
  for (i = 0; i < 3; i++)
    {
      ipos[i] = 0;
      string[i][0] = '\0';
    }
  ok = TRUE;
  
  /* Get the length of the string */
  len = strlen(residue);

  /* Loop over the characters in the string */
  for (ichar = 0; ichar < len; ichar++)
    {
      /* Get this character */
      ch = residue[ichar];

      /* If it is a number, add to residue sequence */
      if (ch >= '0' && ch <= '9')
        {
          /* Set the string we're adding to */
          istring = 1;
        }

      /* Otherwise, add to appropriate string */
      else if (istring != 0)
        istring = 2;

      /* Get the next position in the current string */
      i = ipos[istring];

      /* Save this character */
      string[istring][i] = ch;
      string[istring][i + 1] = '\0';

      /* Increment this string's position */
      ipos[istring]++;
    }

  /* Get the sequence number */
  if (string[1][0] != '\0')
    seq_no = atoi(string[1]);

  /* Otherwise, we have no sequence number, so abort */
  else
    {
      printf("*** ERROR. Invalid variant name: %s\n",residue);
      exit(1);
    }

  /* Initialise */
  istring = 0;
  residue_name[0] = '\0';

  /* Loop to process the original and variant amino acid codes */
  for (loop = 0; loop < 2; loop++)
    {
      /* Get the length of the string */
      len = strlen(string[istring]);
      done = FALSE;
      
      /* If have a 3-letter code, get the 1-letter version */
      if (len == 3)
        {
          /* If second and third characters are lower case, we have a
             residue name */
          if (string[istring][1] >= 'a' && string[istring][1] <= 'z' &&
              string[istring][2] >= 'a' && string[istring][2] <= 'z')
            {
              /* Save the three-character code */
              strcpy(res_name,string[istring]);

              /* Check for termination */
              if (loop == 1 && !strcmp(string[istring],"Ter"))
                {
                  /* Have a termination, so set type */
                  type = TERMINATION;
                  aa_mut = MUT_CODE[TERMINATION];
                }

              /* Otherwise, proces the residue name */
              else
                {
                  /* Convert to upper case */
                  to_upper(res_name,3);

                  /* Get the corresponding 1-character code */
                  if (loop == 0)
                    aa_code = get_aa_code(res_name);
                  else
                    aa_mut = get_aa_code(res_name);
                }

              /* Set flag */
              done = TRUE;
            }
        }

      /* If not processed, then take to be a single a.a. code or sequence
         of codes */
      if (done == FALSE)
        {
          /* Save the single-character code */
          if (loop == 0)
            aa_code = string[istring][0];
          else
            aa_mut = string[istring][0];

          /* Get the corresponding residue name */
          get_res_name(string[istring][0],res_name);

          /* If have more a.a. codes, then save as extension */
          if (len > 1)
            {
              if (loop == 0)
                store_string(&aa_extra,string[istring] + 1);
              else
                /* Copy the extra sequence */
                store_string(&mut_extra,string[istring] + 1);

              /* Save variant type */
              type = NS_SEQ;
            }
        }
        
      /* Format residue name */
      to_lower(res_name,3);
      to_upper(res_name,1);

      /* Add to name of variant */
      if (loop == 0)
        {
          strcpy(residue_name,res_name);
          strcat(residue_name,string[1]);
        }
      else
        strcat(residue_name,res_name);

      /* Increment with string is being processed */
      istring = istring + 2;
    }

  /* If we have variant name, save it */
  if (residue_name[0] != '\0')
    store_string(&added_variant,residue_name);

  /* If we have a sequence substitution, check lengths are the same */
  if (type == NS_SEQ)
    {
      /* Check that added lengths agree */
      if (strlen(aa_extra) != strlen(mut_extra))
        {
          /* Delete this one */
          seq_no = 0;
                      
          /* Show warning */
          printf("*** Warning. Mismatch in variant lengths: "
                 "[%s] (%d)  !=  [%s] (%d). Ignored\n",
                 aa_extra,(int) strlen(aa_extra),mut_extra,
                 (int) strlen(mut_extra));
        }
      else
        mut_length = strlen(aa_extra) + 1;
    }

  /* Return parameters */
  *aac = aa_code;
  *aam = aa_mut;
  *add_var = added_variant;
  *aa_ex = aa_extra;
  *mut_ex = mut_extra;
  *typ = type;
  *mut_len = mut_length;

  /* Return the sequence number */
  return(seq_no);
}
/***********************************************************************

create_c_cath_entry  -  Create and initialise a new c_cath record

***********************************************************************/

struct c_cath *create_c_cath_entry(struct c_cath **fst_c_cath_ptr,
                                   struct c_cath **lst_c_cath_ptr,
                                   char *code,char *name)
{
  struct c_cath *c_cath_ptr, *first_c_cath_ptr, *last_c_cath_ptr;

  /* Initialise c_cath pointers */
  c_cath_ptr = NULL;
  first_c_cath_ptr = *fst_c_cath_ptr;
  last_c_cath_ptr = *lst_c_cath_ptr;

  /* Allocate memory for structure to hold c_cath info */
  c_cath_ptr = (struct c_cath *) malloc(sizeof(struct c_cath));
  if (c_cath_ptr == NULL)
    {
      printf("*** Can't allocate memory for struct c_cath\n");
      exit (1);
    }

  /* If this is the very first, save its pointer */
  if (first_c_cath_ptr == NULL)
    first_c_cath_ptr = c_cath_ptr;

  /* Add link from previous c_cath to the current one */
  if (last_c_cath_ptr != NULL)
    last_c_cath_ptr->next_c_cath_ptr = c_cath_ptr;
  last_c_cath_ptr = c_cath_ptr;

  /* Store data */
  store_string(&(c_cath_ptr->code),code);
  store_string(&(c_cath_ptr->name),name);

  /* Initialise pointer to the next c_cath */
  c_cath_ptr->next_c_cath_ptr = NULL;

  /* Return current pointers */
  *fst_c_cath_ptr = first_c_cath_ptr;
  *lst_c_cath_ptr = last_c_cath_ptr;

  /* Return new c_cath pointer */
  return(c_cath_ptr);
}
/***********************************************************************

read_cath_titles  -  Read in the CATH titles from the titles file

***********************************************************************/

int read_cath_titles(char *file_name,struct c_cath **fst_c_cath_ptr)
{
  char input_line[C_LINELEN + 1];
  char code[C_STRING_LEN + 1], name[C_STRING_LEN + 1];

  char *string_ptr;

  int len, line, ncath;

  struct c_cath *c_cath_ptr, *first_c_cath_ptr, *last_c_cath_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  line = ncath = 0;

  /* Initialise pointers */
  c_cath_ptr = first_c_cath_ptr = last_c_cath_ptr = NULL;
  *fst_c_cath_ptr = NULL;

  /* Open the file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    {
      /* If can't find file, then return empty-handed */
      printf("*** Warning. Unable to find input file: %s\n",
             file_name);
      return(ncath);
    }

  /* Loop while reading in records from the file */
  while (fgets(input_line,C_LINELEN,file_ptr) != NULL)
    {
      /* Truncate the string to strip off the newline */
      len = string_chomp(input_line);

      /* Find first space */
      string_ptr = strchr(input_line,' ');

      /* If found, save this record */
      if (string_ptr != NULL)
        {
          /* Save the CATH name */
          strcpy(name,string_ptr + 1);
          string_ptr[0] = '\0';
          strcpy(code,input_line);

          /* Create a c_cath record */
          c_cath_ptr
            = create_c_cath_entry(&first_c_cath_ptr,&last_c_cath_ptr,
                                  code,name);

          /* Increment count of CATH titles */
          ncath++;
        }
    }

  /* Return current pointer */
  *fst_c_cath_ptr = first_c_cath_ptr;

  /* Close the input file */
  fclose(file_ptr);

  /* Return number of CATH titles read in */
  return(ncath);
}
/***********************************************************************

format_residue_number - Right-justify residue number, taking into
                        account whether it includes an insertion code

***********************************************************************/

void format_residue_number(char residue_number[6],char *token)
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

extract_domain_info - Extract the domain data from CATH domain details
                      line

**********************************************************************/

int extract_domain_info(char *input_line,char *wolf_code,char *cath_code,
                        char domain_details[C_MAXDOMCHAIN][C_DOMLEN],
                        int ndetails)
{
  char chain, number_string[20], first_res_num[6], last_res_num[6];
  char tmp[C_STRING_LEN];

  char *string_ptr, *token[C_MAX_TOKEN], *subtoken[C_MAX_TOKEN];

  int domain, itoken, len, nfrags, ntokens, nsubtokens;

  /* Initialise variables */
  cath_code[0] = '\0';
  chain = input_line[4];
  nfrags = 0;
  wolf_code[0] = '\0';

  /* Split line up into space-separated tokens */
  ntokens = get_tokens(input_line,token,C_MAX_TOKEN,' ');

  /* If have enough tokens, extract the data */
  if (ntokens > 3)
    {
      /* Save the Wolf code */
      strcpy(wolf_code,token[0]);

      /* Save the CATH code */
      strcpy(cath_code,token[2]);

      /* Pull out the domain number */
      len = strlen(wolf_code);
      if (len == 7)
        {
          strcpy(number_string,wolf_code + 5);
          number_string[2] = '\0';
          domain = atoi(number_string);
        }
      else
        domain = 0;

      /* Extract the domain ranges from the last token */
      nsubtokens = get_tokens(token[3],subtoken,C_MAX_TOKEN,',');

      /* Loop over the subtokens to process */
      for (itoken = 0; itoken < nsubtokens; itoken++)
        {
          /* Save the token */
          strcpy(tmp,subtoken[itoken]);

          /* Get the chain id */
          string_ptr = strchr(tmp,':');

          /* If valid, then strip out the chain id */
          if (string_ptr != NULL)
            {
              /* Save the chain */
              chain = string_ptr[1];
              string_ptr[0] = '\0';

              /* Get the residue range */
              string_ptr = strchr(tmp,'-');

              /* If OK, separate the residue numbers */
              if (string_ptr != NULL)
                {
                  /* Get the end residue */
                  strcpy(number_string,string_ptr + 1);
                  string_ptr[0] = '\0';

                  /* Format the residue number as a 5-character string */
                  format_residue_number(last_res_num,number_string);

                  /* Get the start residue */
                  strcpy(number_string,tmp);

                  /* Format the residue number as a 5-character string */
                  format_residue_number(first_res_num,number_string);

                  /* Form the domain details string */
                  sprintf(domain_details[ndetails + nfrags],"D%2d%c%s%c%s",
                          domain,chain,first_res_num,chain,last_res_num);

                  /* Increment count of domain fragments */
                  nfrags++;
                }
            }
        }

      /* Free the subtokens */
      free_tokens(subtoken,nsubtokens);
    }

  /* Free the tokens */
  free_tokens(token,ntokens);

  /* Return numbers of domain segments stored */
  return(nfrags);
}  
/***********************************************************************

create_c_score_entry  -  Create an Alpha Fold score record

***********************************************************************/

struct c_score *create_c_score_entry(struct c_score **fst_c_score_ptr,
                                     struct c_score **lst_c_score_ptr,
                                     char *res_name,char *res_num,
                                     char chain,double score,int range)
{
  struct c_score *c_score_ptr;
  struct c_score *first_c_score_ptr, *last_c_score_ptr;

  /* Initialise c_score pointers */
  c_score_ptr = NULL;
  first_c_score_ptr = *fst_c_score_ptr;
  last_c_score_ptr = *lst_c_score_ptr;

  /* Create c_score entry */
  c_score_ptr = (struct c_score*)malloc(sizeof(struct c_score));
  if (c_score_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for");
      printf(" struct c_score\n");
      printf("***        Program terminated with");
      printf(" error.\n");
      exit (1);
    }

  /* Store the UniProt accession code and nunmber of sequences */
  strcpy(c_score_ptr->res_name,res_name);
  strcpy(c_score_ptr->res_num,res_num);
  c_score_ptr->chain = chain;
  c_score_ptr->score = score;
  c_score_ptr->range = range;

  /* Initialise pointer to next entry */
  c_score_ptr->next_c_score_ptr = NULL;

  /* If this is the first entry stored, then update
     pointer to head of linked list */
  if (first_c_score_ptr == NULL)
    first_c_score_ptr = c_score_ptr;

  /* Otherwise, make previous entry point to the current
     one */
  else
    last_c_score_ptr->next_c_score_ptr = c_score_ptr;
                      
  /* Save the current residue's pointer */
  last_c_score_ptr = c_score_ptr;

  /* Return current pointers */
  *fst_c_score_ptr = first_c_score_ptr;
  *lst_c_score_ptr = last_c_score_ptr;
  return(c_score_ptr);
}
/***********************************************************************

get_c_scores  -  Read in the confidence scores from the AF_scores.dat
                 file

***********************************************************************/

int get_c_scores(char *pdb_code,char *pdbsum_dir,
                 struct c_score **fst_c_score_ptr)
{
  char file_name[C_FILENAME_LEN], input_line[C_LINELEN + 1];
  char chain, res_name[4], res_num[6], number_string[20];

  char *string_ptr;

  int i, len, nc_score, range;

  double score;

  struct c_score *c_score_ptr, *first_c_score_ptr, *last_c_score_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  nc_score = 0;

  /* Initialise c_score pointers */
  *fst_c_score_ptr = first_c_score_ptr = last_c_score_ptr = NULL;

  /* Form full name of the AF_scores.dat file */
  strcpy(file_name,pdbsum_dir);
  strcat(file_name,"AF_scores.dat");
  printf("Reading %s ...\n",file_name);
  fflush(stdout);

  /* Open the file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    return(nc_score);    

  /* Loop through the file, line by line */
  while (fgets(input_line,C_LINELEN,file_ptr) != NULL)
    {
      /* Truncate the input line */
      len = string_truncate(input_line,C_LINELEN);

      /* Extract the residue details */
      strncpy(res_name,input_line,3);
      res_name[3] = '\0';
      strncpy(res_num,input_line + 4,5);
      res_num[5] = '\0';
      chain = input_line[10];

      /* Get the score */
      strncpy(number_string,input_line + 13,7);
      number_string[7] = '\0';
      score = atof(number_string);

      /* Get the corresponding confidence range */
      number_string[0] = input_line[21];
      number_string[1] = '\0';
      range = atoi(number_string);

      /* Create new c_score record to store the details */
      c_score_ptr
        = create_c_score_entry(&first_c_score_ptr,
                               &last_c_score_ptr,res_name,
                               res_num,chain,score,range);
                  
      /* Increment count of records created */
      nc_score++;
    }

  /* Close the file */
  fclose(file_ptr);

  /* Return current pointer */
  *fst_c_score_ptr = first_c_score_ptr;

  /* Return number of c_score records created */
  printf("Number of AF scores retrieved       %8d\n",nc_score);
  return(nc_score);
}
