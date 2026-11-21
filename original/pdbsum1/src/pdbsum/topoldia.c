/***********************************************************************

  topoldia.c - Program to read in promotif.dat and hera.dat files and
               write out a PostScript plot showing a topology diagram

************************************************************************

Date:-         28 Oct 2005

Written by:-   Roman Laskowski
Dir update:-   19 Nov 2012

----------------------------------------------------------------------*/

/* Datafiles
   ---------

Actual name      Description
-----------      -----------
CATHPARAM        Input parameter file giving all the  directory paths for
                 the files used by the program
hera.dat         Input file created by p_hera3.f giving the layout of the
                 secondary structure elements as computed by HERA
promotif.dat     Input file containing all the PROMOTIF motifs
topoldia.ps      Output PostScript file of the topology diagram
topoldia.tdf     Output tdf file version of the topology diagram

*/

/* Function calling-tree
   ---------------------

main
   -> get_command_arguments
         -> store_string
   -> c_get_paths (in common.c)
         -> get_param (in reapar.c)
         -> store_string (in common.c)
   -> find_pdbsum_dir (in common.c)
   -> read_save_file
         -> string_truncate (in common.c)
         -> get_tokens (in common.c)
         -> create_sse_record
   -> get_promotif_data
         -> string_truncate (in common.c)
         -> create_helint_record
         -> get_strand_no
         -> find_promotif_record
         -> create_promotif_record
   -> transfer_labels
   -> get_hera_data
         -> string_truncate (in common.c)
         -> create_sse_record
         -> remove_blanks (in common.c)
         -> create_termini
               -> create_sse_record
   -> separate_sses
   -> layout_helices
   -> create_sort_index
   -> sort_sses
   -> get_max_min
   -> create_mask_array
   -> add_joins
         -> find_path1
               -> check_crossover
         -> find_path2
               -> check_crossover
         -> find_path3
               -> check_crossover
         -> try_paths
               -> check_crossover
         -> create_sse_record
         -> update_mask_array
   -> get_max_min
   -> initialise_psparams
   -> write_postscript
         -> init_pspage
               -> adjust_for_landscape
         -> start_new_page
               -> psendp_
               -> psclos_
               -> pspage_
                     -> psrot_
               -> psopen_
                     -> get_col_rgb
         -> plot_sses
               -> plot_strand
                     -> pscomm_
                     -> pslwid_
                     -> pshade_
                     -> psubox_
                     -> pscolb_
                     -> psline_
                     -> psutri_
                     -> psitxt_
               -> plot_helix
                     -> pscomm_
                     -> pslwid_
                     -> pshade_
                     -> psubox_
                     -> pscolb_
                     -> psline_
                     -> pselip_
                     -> psearc_
                     -> psitxt_
               -> plot_terminus
                     -> pscomm_
                     -> pslwid_
                     -> pshade_
                     -> psubox_
                     -> pscolb_
                     -> psline_
                     -> psctxt_
               -> plot_coil
                     -> pscomm_
                     -> pslwid_
                     -> pscolb_
                     -> psline_
         -> close_page
   -> write_save_file
         -> open_output_file (in common.c)


Compile as follows:-
------------------

cc -c topoldia.c
cc -c common.c
cc -c reapar.c
cc -o topoldia topoldia.o common.o reapar.o -lm -lz

Optimization flags:-  

cc -O2 -funroll-loops -c topoldia.c
cc -O2 -funroll-loops -c common.c
cc -O2 -funroll-loops -c reapar.c
cc -o topoldia topoldia.o common.o reapar.o -lm -lz

*/

#include "common.h"
#include "topoldia.h"

/* Prototypes
   ---------- */

/* In common.c */
void c_get_paths(struct c_file_data *file_name_ptr,int run_64);
int find_pdbsum_dir(char *pdb_code,char *pdbsum_template[NPDBSUM_DIR],
                    int pdb_type,char *dir_name);
int get_tokens(char *line,char **token,int max_token,char token_sep);
FILE *open_output_file(char *out_name,int terminate);
void remove_blanks(char *string);
void store_string(char **string_ptr,char *string);
int string_truncate(char *string,int max_length);

/* Function prototypes from reapar.c */
int get_param(char *token_name,char *file_name,int filename_len,
              int exit_on_error);


/***********************************************************************

get_command_arguments  -  Interpret the command line arguments

***********************************************************************/

void get_command_arguments(char *strptrs[],int num,
                           struct parameters *params)
{
  int itoken;

  /* Initialise variables */
  params->pdb_code[0] = '\0';
  params->chain = ' ';
  params->domain = 0;
  params->show_residue_nos = FALSE;
  params->show_sse_labels = FALSE;
  params->plt_file = &null;
  params->use_data = FALSE;
  params->save_data = FALSE;
  params->local = FALSE;
  params->new = FALSE;
  params->pdb_type = STANDARD_PDB;
  params->run_64 = FALSE;

  /* If nothing entered on the command line, show options available */
  if (num < 1)
    {
      printf("\n");
      printf("Correct usage is:\n");
      printf("\n");
      printf("    topoldia  pdb_code  [-save]  [-use filename]  [-resno]  "
             "[-labs]  [-chain X]  [-domain n]\n"
             "               [-new]  [UPLOAD | MCSG | AF] [-local]  [-64]  "
             "[-pdbsum1]\n");
      printf("\n");
      printf("where\n");
      printf("  * pdb_code     = 4-character PDB code\n");
      printf("  * -chain X     = Chain id\n");
      printf("  * -domain n    = Colour by given domain number\n");
      printf("  * -resno       = Show residue numbers on the SSEs\n");
      printf("  * -labs        = Show labels on the SSEs\n");
      printf("  * -save        = Save data in topoldia.plt file\n");
      printf("  * -use filename = Read data from given .plt file in PDBsum ");
      printf("directory\n");
      printf("  * UPLOAD       = Uploaded PDB file for PDBsum generation\n");
      printf("  * MCSG         = MCSG structure\n");
      printf("  * AF           = Alpha Fold structure\n");
      printf("  * -pdbsum1     = PDBsum1 structure\n");
      printf("  * -local       = PDBsum files are in current directory\n");
      printf("  * -new         = use new layout algorithm rather than Hera "
             "plot\n");
      printf("  * -64          = run on 64-bit processor\n");
      printf("\n");
      printf("\n");
      exit(1);
    }

  /* Get the PDB code */
  strcpy(params->pdb_code,strptrs[1]);

  /* Loop through any remaining command arguments */
  for (itoken = 2; itoken < num + 1; itoken++)
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

      /* Check for local files option */
      else if (!strcmp(strptrs[itoken],"-local"))
        params->local = TRUE;

      /* Check for the -resno option */
      else if (!strcmp(strptrs[itoken],"-resno"))
        params->show_residue_nos = TRUE;

      /* Check for the -labs option */
      else if (!strcmp(strptrs[itoken],"-labs"))
        params->show_sse_labels = TRUE;

      /* Check for chain id */
      else if (!strcmp(strptrs[itoken],"-chain"))
        {
          /* Get the id from the next token, if there is one */
          if (itoken < num)
            {
              params->chain = strptrs[itoken + 1][0];
              itoken++;
            }          
        }

      /* Check for domain number */
      else if (!strcmp(strptrs[itoken],"-domain"))
        {
          /* Get the number from the next token, if there is one */
          if (itoken < num)
            {
              params->domain = atoi(strptrs[itoken + 1]);
              itoken++;
            }          
        }

      /* Check for the -save option */
      else if (!strcmp(strptrs[itoken],"-save"))
        params->save_data = TRUE;

      /* Check for the -new option */
      else if (!strcmp(strptrs[itoken],"-new"))
        params->new = TRUE;

      /* Check for the -use option */
      else if (!strcmp(strptrs[itoken],"-use"))
        {
          /* Get the file name from the next token, if there is one */
          if (itoken < num)
            {
              store_string(&(params->plt_file),strptrs[itoken + 1]);
              itoken++;

              /* Set flag */
              params->use_data = TRUE;
            }          
        }

      /* Check for the 64-bit option */
      else if (!strcmp(strptrs[itoken],"-64"))
        params->run_64 = TRUE;
    }
}
/***********************************************************************

create_sse_record  -  Create and initialise a new sse record

***********************************************************************/

struct sse *create_sse_record(struct sse **fst_sse_ptr,
                              struct sse **lst_sse_ptr,int type,
                              int first_pos,int last_pos,int *xpos,
                              int *ypos,int endseg)
{
  struct sse *sse_ptr, *first_sse_ptr, *last_sse_ptr;

  /* Initialise sse pointers */
  sse_ptr = NULL;
  first_sse_ptr = *fst_sse_ptr;
  last_sse_ptr = *lst_sse_ptr;

  /* Allocate memory for structure to hold sse info */
  sse_ptr = (struct sse *) malloc(sizeof(struct sse));
  if (sse_ptr == NULL)
    {
      printf("*** Can't allocate memory for struct sse\n");
      exit (1);
    }

  /* If this is the very first sse, then store pointer */
  if (first_sse_ptr == NULL)
    first_sse_ptr = sse_ptr;

  /* Add link from previous sse to the current one */
  if (last_sse_ptr != NULL)
    last_sse_ptr->next_sse_ptr = sse_ptr;
  last_sse_ptr = sse_ptr;

  /* Save the original element coords */
  sse_ptr->xorig[0] = xpos[0];
  sse_ptr->xorig[1] = xpos[1];
  sse_ptr->yorig[0] = ypos[0];
  sse_ptr->yorig[1] = ypos[1];

  /* Determine the orientation of this SSE */
  if (xpos[0] == xpos[1])
    sse_ptr->orientation = VERTICAL;
  else
    sse_ptr->orientation = HORIZONTAL;

  /* Transfer the data fields to the new sse record */
  sse_ptr->type = type;
  sse_ptr->sort = 0;
  sse_ptr->endseg = endseg;
  if (sse_ptr->orientation == VERTICAL)
    {
      sse_ptr->xpos[0] = xpos[0];
      sse_ptr->xpos[1] = xpos[1];
      if (first_pos > last_pos)
        {
          sse_ptr->first_pos = last_pos;
          sse_ptr->last_pos = first_pos;
          sse_ptr->ypos[0] = ypos[1];
          sse_ptr->ypos[1] = ypos[0];
        }
      else
        {
          sse_ptr->first_pos = first_pos;
          sse_ptr->last_pos = last_pos;
          sse_ptr->ypos[0] = ypos[0];
          sse_ptr->ypos[1] = ypos[1];
        }
      if (sse_ptr->ypos[0] < sse_ptr->ypos[1])
        sse_ptr->direction = FORWARD;
      else
        sse_ptr->direction = BACKWARD;
    }
  else
    {
      sse_ptr->ypos[0] = ypos[0];
      sse_ptr->ypos[1] = ypos[1];
      if (first_pos > last_pos)
        {
          sse_ptr->first_pos = last_pos;
          sse_ptr->last_pos = first_pos;
          sse_ptr->xpos[0] = xpos[1];
          sse_ptr->xpos[1] = xpos[0];
        }
      else
        {
          sse_ptr->first_pos = first_pos;
          sse_ptr->last_pos = last_pos;
          sse_ptr->xpos[0] = xpos[0];
          sse_ptr->xpos[1] = xpos[1];
        }
      if (sse_ptr->xpos[0] < sse_ptr->xpos[1])
        sse_ptr->direction = FORWARD;
      else
        sse_ptr->direction = BACKWARD;
    }

  /* Determine extent of this SSE for use in cross-over checks */
  sse_ptr->xmin = sse_ptr->xpos[0];
  sse_ptr->xmax = sse_ptr->xpos[0];
  if (sse_ptr->xpos[1] < sse_ptr->xmin)
    sse_ptr->xmin = sse_ptr->xpos[1];
  if (sse_ptr->xpos[1] > sse_ptr->xmax)
    sse_ptr->xmax = sse_ptr->xpos[1];
  sse_ptr->ymin = sse_ptr->ypos[0];
  sse_ptr->ymax = sse_ptr->ypos[0];
  if (sse_ptr->ypos[1] < sse_ptr->ymin)
    sse_ptr->ymin = sse_ptr->ypos[1];
  if (sse_ptr->ypos[1] > sse_ptr->ymax)
    sse_ptr->ymax = sse_ptr->ypos[1];

  /* Add SSE size */
  sse_ptr->xmin = sse_ptr->xmin - SSE_WIDTH;
  sse_ptr->xmax = sse_ptr->xmax + SSE_WIDTH;

  /* Initialise remaining fields */
  sse_ptr->first_resnum[0] = '\0';
  sse_ptr->last_resnum[0] = '\0';
  sse_ptr->sst_id = &null;

  /* Initialise pointer to next record */
  sse_ptr->next_sse_ptr = NULL;

  /* Return current pointers */
  *fst_sse_ptr = first_sse_ptr;
  *lst_sse_ptr = last_sse_ptr;

  /* Return pointer to the new sse record */
  return(sse_ptr);
}
/***********************************************************************

read_save_file  -  Read in the details from the given .plt file

***********************************************************************/

int read_save_file(struct sse **fst_sse_ptr,struct parameters *params)
{
  char filename[FILENAME_LEN], input_line[LINELEN + 1], number_string[20];
  char first_resnum[6], last_resnum[6];

  char *string_ptr, *token[MAX_TOKEN];

  int first_pos, last_pos, xpos[2], ypos[2];
  int itoken, len, nsse, ntokens, type;
  int endseg;

  struct sse *last_sse_ptr, *first_sse_ptr, *sse_ptr;

  struct residue *residue_array_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  strcpy(filename,params->plt_file);
  nsse = 0;

  /* Initialise pointers */
  first_sse_ptr = NULL;
  last_sse_ptr = NULL;
  sse_ptr = NULL;

  /* Open the input file */
  if ((file_ptr = fopen(filename,"r")) ==  NULL)
    {
      /* If file not found, have error */
      printf("*** ERROR. Unable to open file: %s\n",filename);
      exit(-1);
    }

  /* Read in the very first line giving the chain and domain numbers */
  if (fgets(input_line,LINELEN,file_ptr) != NULL)
    {
      /* Truncate the string to strip off the newline */
      len = string_truncate(input_line,LINELEN);

      /* If line correct, pick up the required details */
      if (!strncmp(input_line,"Chain ",6))
        {
          /* Extract the chain and domain */
          params->chain = input_line[7];
          strcpy(number_string,input_line+17);
          params->domain = atoi(number_string);
        }
    }

  /* Loop while reading in records from the file */
  while (fgets(input_line,LINELEN,file_ptr) != NULL)
    {
      /* Truncate the string to strip off the newline */
      len = string_truncate(input_line,LINELEN);

      /* Read in all the tokens of this line */
      ntokens = get_tokens(input_line,token,MAX_TOKEN,' ');

      /* If have at least 7 tokens, then create an SSE record */
      if (ntokens >= 8)
        {
          /* Extract the key elements */
          type = atoi(token[0]);
          xpos[0] = atoi(token[1]);
          ypos[0] = atoi(token[2]);
          xpos[1] = atoi(token[3]);
          ypos[1] = atoi(token[4]);
          first_pos = atoi(token[5]);
          last_pos = atoi(token[6]);
          endseg = atoi(token[7]);

          /* Get the start and end residues, if present */
          if (ntokens > 8)
            strcpy(first_resnum,token[8]);
          else
            first_resnum[0] = '\0';
          if (ntokens > 9)
            strcpy(last_resnum,token[9]);
          else
            last_resnum[0] = '\0';

          /* Create a line SSE */
          sse_ptr
            = create_sse_record(&first_sse_ptr,&last_sse_ptr,type,
                                first_pos,last_pos,xpos,
                                ypos,endseg);

          /* Store the first and last residue numbers */
          strcpy(sse_ptr->first_resnum,first_resnum);
          strcpy(sse_ptr->last_resnum,last_resnum);
        }
    }

  /* Close the input file */
  fclose(file_ptr);

  /* Return pointers */
  *fst_sse_ptr = first_sse_ptr;

  /* Return the number of SSEs */
  return(nsse);
}
/***********************************************************************

create_helint_record  -  Create and initialise a new helint record

***********************************************************************/

struct helint *create_helint_record(struct helint **fst_helint_ptr,
                                    struct helint **lst_helint_ptr,
                                    float distance,float omega,
                                    char int_type1,char int_type2,
                                    int nint)
{
  struct helint *helint_ptr, *first_helint_ptr, *last_helint_ptr;

  /* Initialise helint pointers */
  helint_ptr = NULL;
  first_helint_ptr = *fst_helint_ptr;
  last_helint_ptr = *lst_helint_ptr;

  /* Allocate memory for structure to hold helint info */
  helint_ptr = (struct helint *) malloc(sizeof(struct helint));
  if (helint_ptr == NULL)
    {
      printf("*** Can't allocate memory for struct helint\n");
      exit (1);
    }

  /* If this is the very first helint, then store pointer */
  if (first_helint_ptr == NULL)
    first_helint_ptr = helint_ptr;

  /* Add link from previous helint to the current one */
  if (last_helint_ptr != NULL)
    last_helint_ptr->next_helint_ptr = helint_ptr;
  last_helint_ptr = helint_ptr;

  /* Transfer the data fields to the new helint record */
  helint_ptr->distance = distance;
  helint_ptr->omega = omega;
  helint_ptr->int_type1 = int_type1;
  helint_ptr->int_type2 = int_type2;
  helint_ptr->nint = nint;

  /* Initialise remaining fields */
  helint_ptr->promotif1_ptr = NULL;
  helint_ptr->promotif2_ptr = NULL;

  /* Initialise pointer to next record */
  helint_ptr->next_helint_ptr = NULL;

  /* Return current pointers */
  *fst_helint_ptr = first_helint_ptr;
  *lst_helint_ptr = last_helint_ptr;

  /* Return pointer to the new helint record */
  return(helint_ptr);
}
/***********************************************************************

find_promotif_record  -  Find the required promotif record

***********************************************************************/

struct promotif *find_promotif_record(struct promotif *first_promotif_ptr,
                                      int helix)
{
  struct promotif *promotif_ptr, *match_promotif_ptr;

  /* Initialise variables */
  match_promotif_ptr = NULL;

  /* Get pointer to first promotif record */
  promotif_ptr = first_promotif_ptr;

  /* Loop through promotif records to find out */
  while (promotif_ptr != NULL && match_promotif_ptr == NULL)
    {
      /* If this is the right record, save the pointer */
      if (promotif_ptr->sec_struc == 'H' &&
          promotif_ptr->identifier[0] == helix)
        match_promotif_ptr = promotif_ptr;

      /* Get next promotif record */
      promotif_ptr = promotif_ptr->next_promotif_ptr;
    }

  /* Return the matched promotif record */
  return(match_promotif_ptr);
}
/***********************************************************************

create_promotif_record  -  Create and initialise a new promotif record

***********************************************************************/

struct promotif *create_promotif_record(struct promotif **fst_promotif_ptr,
                                        struct promotif **lst_promotif_ptr,
                                        int motif,char sec_struc,
                                        char *res_name1,char *res_name2,
                                        char *res_num1,char *res_num2,
                                        char chain1,char chain2,
                                        int *id)
{
  struct promotif *promotif_ptr, *first_promotif_ptr, *last_promotif_ptr;

  /* Initialise promotif pointers */
  promotif_ptr = NULL;
  first_promotif_ptr = *fst_promotif_ptr;
  last_promotif_ptr = *lst_promotif_ptr;

  /* Allocate memory for structure to hold promotif info */
  promotif_ptr = (struct promotif *) malloc(sizeof(struct promotif));
  if (promotif_ptr == NULL)
    {
      printf("*** Can't allocate memory for struct promotif\n");
      exit (1);
    }

  /* If this is the very first promotif, then store pointer */
  if (first_promotif_ptr == NULL)
    first_promotif_ptr = promotif_ptr;

  /* Add link from previous promotif to the current one */
  if (last_promotif_ptr != NULL)
    last_promotif_ptr->next_promotif_ptr = promotif_ptr;
  last_promotif_ptr = promotif_ptr;

  /* Transfer the data fields to the new promotif record */
  promotif_ptr->motif = motif;
  promotif_ptr->sec_struc = sec_struc;
  strcpy(promotif_ptr->res_name1,res_name1);
  strcpy(promotif_ptr->res_name2,res_name2);
  strcpy(promotif_ptr->res_num1,res_num1);
  strcpy(promotif_ptr->res_num2,res_num2);
  promotif_ptr->chain1 = chain1;
  promotif_ptr->chain2 = chain2;
  promotif_ptr->identifier[0] = id[0];
  promotif_ptr->identifier[1] = id[1];

  /* Initialise pointer to next record */
  promotif_ptr->next_promotif_ptr = NULL;

  /* Return current pointers */
  *fst_promotif_ptr = first_promotif_ptr;
  *lst_promotif_ptr = last_promotif_ptr;

  /* Return pointer to the new promotif record */
  return(promotif_ptr);
}
/***********************************************************************

get_strand_no  -  Get the strand number within this sheet 

***********************************************************************/

int get_strand_no(struct promotif *first_promotif_ptr,char chain,
                  int sheet_id)
{
  int strand;

  struct promotif *promotif_ptr;

  /* Initialise the strand counter */
  strand = 0;

  /* Get starting promotif record */
  promotif_ptr = first_promotif_ptr;

  /* Loop through promotifs to process all beta strands */
  while (promotif_ptr != NULL)
    {
      /* If this is a strand, see if it is our one */
      if (promotif_ptr->motif == STRANDS)
        {
          /* If strand belongs to our chain and sheet, then save strand
             number */
          if (promotif_ptr->chain1 == chain &&
              promotif_ptr->identifier[0] == sheet_id)
            strand = promotif_ptr->identifier[1];
        }

      /* Get the next promotif */
      promotif_ptr = promotif_ptr->next_promotif_ptr;
    }

  /* Increment strand number */
  strand++;

  /* Return the strand number in the sheet */
  return(strand);
}
/***********************************************************************

get_promotif_data  -  Read in the PROMOTIF analyses from the promotif.dat
                      file

***********************************************************************/

int get_promotif_data(char *pdb_code,char *pdbsum_dir,
                      struct promotif **fst_promotif_ptr,
                      struct helint **fst_helint_ptr,int *top_sht,
                      struct parameters params)
{
  char file_name[FILENAME_LEN], input_line[LINELEN], number_string[20];
  char chain, chain1, chain2, res_num1[6], res_num2[6];
  char res_name1[4], res_name2[4];
  char int_type1, int_type2;

  int i, id[2], motif, npromotif, top_sheet;
  int helix1, helix2, nhelint, nint;

  float distance, omega;

  struct helint *last_helint_ptr, *first_helint_ptr, *helint_ptr;
  struct promotif *last_promotif_ptr, *first_promotif_ptr;
  struct promotif *promotif_ptr;

  FILE *file_ptr;

  static char sec_struc[]={"HiEbgD"};

  /* Initialise variables */
  nhelint = npromotif = 0;
  *top_sht = top_sheet = 0;

  /* Initialise pointers */
  helint_ptr = first_helint_ptr = last_helint_ptr = NULL;
  promotif_ptr = first_promotif_ptr = last_promotif_ptr = NULL;

  /* Initialise variables */
  id[0] = id[1] = 0;
  motif = -1;
  strcpy(res_name1,"XXX");
  strcpy(res_name2,"XXX");

  /* Form name of promotif.dat file */
  strcpy(file_name,pdbsum_dir);
  if (params.pdb_type == PDBSUM1)
    strcat(file_name,"promotif/");
  strcat(file_name,"promotif.dat");

  /* Open the PROMOTIF file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    printf("*** Warning. Unable to find PROMOTIF file: %s\n",file_name);

  /* If have the file, read in all its data */
  else
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
                if (!strcmp(input_line+6,MOTIF_NAME[i]))
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
              id[0] = id[1] = 0;

              /* Process according to the different record formats */
              switch (motif)
                {
                  /* a. Helix */
                case 0 :
                  strncpy(number_string,input_line+3,3);
                  number_string[3] = '\0';
                  id[0] = atoi(number_string);
                  strncpy(res_num1,input_line+8,5);
                  res_num1[5] = '\0';
                  strncpy(res_num2,input_line+20,5);
                  res_num2[5] = '\0';
                  break;

                  /* b. Helix-helix interactions */
                case 1 :
                  chain1 = input_line[3];
                  chain2 = input_line[14];
                  strncpy(number_string,input_line+25,5);
                  number_string[5] = '\0';
                  distance = atof(number_string);
                  strncpy(number_string,input_line+38,6);
                  number_string[6] = '\0';
                  omega = atof(number_string);
                  int_type1 = input_line[32];
                  int_type2 = input_line[35];
                  strncpy(number_string,input_line+46,3);
                  number_string[3] = '\0';
                  nint = atoi(number_string);

                  /* Create a helix-interaction record */
                  helint_ptr =
                    create_helint_record(&first_helint_ptr,&last_helint_ptr,
                                         distance,omega,int_type1,int_type2,
                                         nint);

                  /* Increment count of helix interaction records */
                  nhelint++;

                  /* Get the helix numbers of the two interacting helices */
                  strncpy(number_string,input_line+6,3);
                  number_string[3] = '\0';
                  helix1 = atoi(number_string);
                  strncpy(number_string,input_line+17,3);
                  number_string[3] = '\0';
                  helix2 = atoi(number_string);

                  /* Locate the promotif records for these two helices, and
                     store their pointers */
                  helint_ptr->promotif1_ptr
                    = find_promotif_record(first_promotif_ptr,helix1);
                  helint_ptr->promotif2_ptr
                    = find_promotif_record(first_promotif_ptr,helix2);
                  break;

                  /* c. Beta strand */
                case 2 :
                  strncpy(res_num1,input_line+8,5);
                  res_num1[5] = '\0';
                  strncpy(res_num2,input_line+20,5);
                  res_num2[5] = '\0';

                  /* Get the sheet this strand belongs to */
                  strncpy(number_string,input_line+32,3);
                  number_string[3] = '\0';
                  id[0] = atoi(number_string);
                  if (id[0] > top_sheet)
                    top_sheet = id[0];

                  /* Get the strand number within this sheet */
                  id[1] = get_strand_no(first_promotif_ptr,chain,id[0]);
                  break;

                  /* d. Beta turns */
                case 3 :
                  strncpy(res_num1,input_line+3,5);
                  res_num1[5] = '\0';
                  strncpy(res_num2,input_line+39,5);
                  res_num2[5] = '\0';
                  break;

                  /* e. Gamma turns */
                case 4 :
                  strncpy(res_num1,input_line+3,5);
                  res_num1[5] = '\0';
                  strncpy(res_num2,input_line+27,5);
                  res_num2[5] = '\0';
                  break;

                  /* f. Disulphide bridge */
                case 5 :
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

                  /* g. Beta hairpin */
                case 6 :
                  strncpy(res_num1,input_line+15,5);
                  res_num1[5] = '\0';
                  strncpy(res_num2,input_line+27,5);
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
                    create_promotif_record(&first_promotif_ptr,
                                           &last_promotif_ptr,
                                           motif,sec_struc[motif],
                                           res_name1,res_name2,
                                           res_num1,res_num2,
                                           chain1,chain2,id);

                  /* Increment count of records created */
                  npromotif++;
                }
            }
        }

      /* Close the file */
      fclose(file_ptr);
    }

  /* Show number of PROMOTIF records created */
  printf("Number of PROMOTIF records created: %d\n",npromotif);

  /* Return pointers to first PROMOTIF and helint entries */
  *fst_promotif_ptr = first_promotif_ptr;
  *fst_helint_ptr = first_helint_ptr;

  /* Return the top sheet */
  *top_sht = top_sheet;

  /* Return number of helix interactions */
  return(nhelint);
}
/***********************************************************************

transfer_labels  -  Transfer the helix/strand labels to the corresponding
                    secondary structure elements

***********************************************************************/

void transfer_labels(struct promotif *first_promotif_ptr,
                     struct sse *first_sse_ptr,char chain,int top_sheet)
{
  char number_string[20], sst_id[20];

  int len, number, res1, res2, sse_res1, sse_res2;
  int done, overlap;

  struct promotif *promotif_ptr;
  struct sse *sse_ptr;

  /* Get starting promotif record */
  promotif_ptr = first_promotif_ptr;

  /* Loop through promotifs to process all beta strands */
  while (promotif_ptr != NULL)
    {
      /* If this is a helix or strand record, find the corresponding SSE
         record and transfer the id */
      if (promotif_ptr->chain1 == chain &&
          (promotif_ptr->motif == HELICES ||
           promotif_ptr->motif == STRANDS))
        {
          /* Get the residue range for this motif */
          strcpy(number_string,promotif_ptr->res_num1);
          if (number_string[4] != ' ')
            number_string[4] = ' ';
          res1 = atoi(number_string);
          strcpy(number_string,promotif_ptr->res_num2);
          if (number_string[4] != ' ')
            number_string[4] = ' ';
          res2 = atoi(number_string);

          /* Form the id string */
          if (promotif_ptr->motif == HELICES)
            sprintf(sst_id,"H%d",promotif_ptr->identifier[0]);
          else
            {
              /* Initialise sheet identifier */
              sst_id[0] = '\0';
              
              /* Get the sheet identifier */
              number = promotif_ptr->identifier[0];

              /* If using letters, convert accordingly */
              if (top_sheet < 27)
                {
                  sst_id[0] = CHAIN_LIST[number - 1];
                  sst_id[1] = '\0';
                }

              /* Otherwise, use number */
              else
                sprintf(sst_id,"%d/",number);

              /* Add strand number to current sst identifier */
              sprintf(number_string,"%d",promotif_ptr->identifier[1]);

              /* Add to overall strand identifier */
              strcat(sst_id,number_string);
            }

          /* Initialise flag */
          done = FALSE;

          /* Get pointer to the first sse entry */
          sse_ptr = first_sse_ptr;

          /* Loop through all entries to find the one corresponding
             to the current motif */
          while (sse_ptr != NULL && done == FALSE)
            {
              /* Check that a helix or strand element */
              if ((sse_ptr->type == HELIX &&
                   promotif_ptr->motif == HELICES) ||
                  (sse_ptr->type == STRAND &&
                   promotif_ptr->motif == STRANDS))
                {
                  /* Get the residue range for this secondary structure
                     element */
                  strcpy(number_string,sse_ptr->first_resnum);
                  len = strlen(number_string);
                  if (number_string[len - 1] != ' ')
                    number_string[len - 1] = ' ';
                  sse_res1 = atoi(number_string);
                  strcpy(number_string,sse_ptr->last_resnum);
                  len = strlen(number_string);
                  if (number_string[len - 1] != ' ')
                    number_string[len - 1] = ' ';
                  sse_res2 = atoi(number_string);

                  /* See if the residue ranges overlap */
                  overlap = FALSE;
                  if ((sse_res1 >= res1 && sse_res1 <= res2) ||
                      (sse_res2 >= res1 && sse_res2 <= res2) ||
                      (res1 >= sse_res1 && res1 <= sse_res2) ||
                      (res2 >= sse_res1 && res2 <= sse_res2))
                    overlap = TRUE;
                  /* If have overlap, then transfer the id */
                  if (overlap == TRUE)
                    {
                      /* Have match, so transfer id */
                      store_string(&(sse_ptr->sst_id),sst_id);

                      /* Set flag that found */
                      done = TRUE;
                    }
                }

              /* Get pointer to the next sse in the linked list */
              sse_ptr = sse_ptr->next_sse_ptr;
            }
        }

      /* Get the next promotif */
      promotif_ptr = promotif_ptr->next_promotif_ptr;
    }
}
/***********************************************************************

create_termini  -  Add the two SSEs representing the N- and C-termini

***********************************************************************/

int create_termini(struct sse **fst_sse_ptr,struct sse **lst_sse_ptr,
                   struct sse *pren_sse_ptr,struct sse *postn_sse_ptr,
                   struct sse *prec_sse_ptr,struct sse *postc_sse_ptr,
                   int nsse,struct residue *residue_array_ptr,int nresid)
{
  int first_pos, last_pos;
  int mid_point, xpos[2], ypos[2];

  struct sse *dummy1_sse_ptr, *dummy2_sse_ptr;
  struct sse *last_sse_ptr, *first_sse_ptr, *sse_ptr;

  /* Initialise pointers */
  first_sse_ptr = *fst_sse_ptr;
  last_sse_ptr = *lst_sse_ptr;
  sse_ptr = NULL;

  // N-terminus

  /* Determine where to place this element */
  xpos[0] = xpos[1] = postn_sse_ptr->xpos[0];

  /* Determine mid-point of terminus box and "residue positions" */
  if (postn_sse_ptr->direction == FORWARD)
    {
      mid_point = postn_sse_ptr->ypos[0] - 2 * CONVERTY;
      // first_pos = postn_sse_ptr->first_pos - 2;
      // last_pos = postn_sse_ptr->first_pos - 1;
    }
  else
    {
      mid_point = postn_sse_ptr->ypos[0] + 2 * CONVERTY;
      // first_pos = postn_sse_ptr->first_pos - 1;
      // last_pos = postn_sse_ptr->first_pos - 2;
    }

  /* Set residue range to the first residue */
  first_pos = last_pos = 0;

  /* Get y-coords of start and end of box */
  ypos[0] = mid_point - TERMINUS_HEIGHT;
  ypos[1] = mid_point + TERMINUS_HEIGHT;

  /* Create a SSE record */
  dummy1_sse_ptr = NULL;
  dummy2_sse_ptr = NULL;
  sse_ptr = create_sse_record(&dummy1_sse_ptr,&dummy2_sse_ptr,NTERMINUS,
                              first_pos,last_pos,xpos,ypos,FALSE);
  nsse++;

  /* Store the residue number of this terminus */
  strcpy(sse_ptr->first_resnum,residue_array_ptr[0].res_num);
  strcpy(sse_ptr->last_resnum,residue_array_ptr[0].res_num);

  /* Store the sort position */
  sse_ptr->sort = -1;

  /* Insert SSE record into appropriate position */
  if (postn_sse_ptr->direction == FORWARD)
    {
      /* Insert prior to first SSE */
      if (pren_sse_ptr == NULL)
        {
          /* Nothing prior to current SSE, so make new SSE the first one */
          first_sse_ptr = sse_ptr;
        }
      else
        {
          /* Get previous SSE to point to current SSE */
          pren_sse_ptr->next_sse_ptr = sse_ptr;
        }

      /* Point to the next SSE */
      sse_ptr->next_sse_ptr = postn_sse_ptr;
    }

  /* If first SSE pointing downwards, need to insert current SSE above
     it on the plot */
  else
    {
      /* Point to the next SSE */
      sse_ptr->next_sse_ptr = postn_sse_ptr->next_sse_ptr;

      /* If there is no SSE beyond this one, then update last SSE pointer */
      if (sse_ptr->next_sse_ptr == NULL)
        last_sse_ptr = sse_ptr;

      /* Alternatively, if the next SSE corresponds to the SSE after the
         C-terminus, need to make this the pre-C-terminus SSE */
      else if (sse_ptr->next_sse_ptr == postc_sse_ptr)
        prec_sse_ptr = sse_ptr;

      /* Point first SSE to current SSE */
      postn_sse_ptr->next_sse_ptr = sse_ptr;
    }

  // C-terminus

  /* Determine where to place this element */
  xpos[0] = xpos[1] = postc_sse_ptr->xpos[0];

  /* Determine mid-point of terminus box and its "residue positions" */
  if (postc_sse_ptr->direction == FORWARD)
    {
      mid_point = postc_sse_ptr->ypos[1] + 2 * CONVERTY;
      first_pos = postc_sse_ptr->last_pos + 1;
      last_pos = postc_sse_ptr->last_pos + 2;
    }
  else
    {
      mid_point = postc_sse_ptr->ypos[1] - 2 * CONVERTY;
      first_pos = postc_sse_ptr->last_pos + 2;
      last_pos = postc_sse_ptr->last_pos + 1;
    }

  /* Get y-coords of start and end of box */
  ypos[0] = mid_point - TERMINUS_HEIGHT;
  ypos[1] = mid_point + TERMINUS_HEIGHT;

  /* Create a SSE record */
  dummy1_sse_ptr = NULL;
  dummy2_sse_ptr = NULL;
  sse_ptr =
    create_sse_record(&dummy1_sse_ptr,&dummy2_sse_ptr,CTERMINUS,
                        first_pos,last_pos,xpos,ypos,FALSE);
  nsse++;

  /* Store the residue number of this terminus */
  strcpy(sse_ptr->first_resnum,residue_array_ptr[nresid - 1].res_num);
  strcpy(sse_ptr->last_resnum,residue_array_ptr[nresid - 1].res_num);

  /* Store the sort position */
  sse_ptr->sort = nresid;

  /* Insert SSE record into appropriate position */
  if (postc_sse_ptr->direction == FORWARD)
    {
      /* Point to the next SSE */
      sse_ptr->next_sse_ptr = postc_sse_ptr->next_sse_ptr;

      /* If there is no SSE beyond this one, then update last SSE pointer */
      if (sse_ptr->next_sse_ptr == NULL)
        last_sse_ptr = sse_ptr;

      /* Point last SSE to current SSE */
      postc_sse_ptr->next_sse_ptr = sse_ptr;
    }

  /* If first SSE pointing downwards, need to insert current SSE below
     it on the plot */
  else
    {
      /* Insert prior to last SSE */
      if (prec_sse_ptr == NULL)
        {
          /* Nothing prior to current SSE, so make new SSE the first one */
          first_sse_ptr = sse_ptr;
        }
      else
        {
          /* Get previous SSE to point to current SSE */
          prec_sse_ptr->next_sse_ptr = sse_ptr;
        }

      /* Point to the next SSE */
      sse_ptr->next_sse_ptr = postc_sse_ptr;
    }

  /* Return pointers */
  *fst_sse_ptr = first_sse_ptr;
  *lst_sse_ptr = last_sse_ptr;

  /* Return number of SSEs */
  return(nsse);
}
/***********************************************************************

get_hera_data  -  Read in the residues and SSE data from the hera.dat
                  file

***********************************************************************/

int get_hera_data(struct residue **res_array_ptr,int *nres,
                  struct sse **fst_sse_ptr,struct sse **lst_sse_ptr,
                  struct promotif *first_promotif_ptr,
                  struct parameters *params)
{
  char file_name[FILENAME_LEN], input_line[LINELEN];
  char chain, inchn, number_string[20], res_name[4], res_num[6];

  int i, ires, len, line, min_res, max_res, nresid, nsse, rcount, type;
  int first_pos, last_pos;
  int ctermx[2], ctermy[2], ntermx[2], ntermy[2];
  int last_posc, first_posn;
  int xpos[2], ypos[2];

  struct sse *last_sse_ptr, *first_sse_ptr, *sse_ptr;
  struct sse *previous_sse_ptr;
  struct sse *pren_sse_ptr, *postn_sse_ptr, *prec_sse_ptr, *postc_sse_ptr;

  struct residue *residue_array_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  chain = ' ';
  ires = 0;
  line = 0;
  min_res = 0;
  max_res = 0;
  nresid = 0;
  nsse = 0;
  rcount = 0;
  last_posc = first_posn = -9;
  for (i = 0; i < 2; i++)
    ctermx[i] = ctermy[i] = ntermx[i] = ntermy[i] = 0;

  /* Initialise pointers */
  first_sse_ptr = NULL;
  last_sse_ptr = NULL;
  sse_ptr = NULL;
  pren_sse_ptr = NULL;
  prec_sse_ptr = NULL;
  postn_sse_ptr = NULL;
  postc_sse_ptr = NULL;
  previous_sse_ptr = NULL;

  /* Form name of hera.dat file */
  strcpy(file_name,"hera.dat");

  /* Open the file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    printf("*** Warning. Unable to find file: %s\n",file_name);

  /* If have the file, read in all its data */
  else
    {
      /* Loop while reading in records from the file */
      while (fgets(input_line,LINELEN,file_ptr) != NULL)
        {
          /* Truncate the string to strip off the newline */
          len = string_truncate(input_line,LINELEN);

          /* Increment line count */
          line++;

          /* If this is the first line, get the chain id, number of
             residues and residue range */
          if (line == 1)
            {
              /* Check that line sufficiently long */
              if (len < 20)
                {
                  printf("*** ERROR. Line 1 of hera.dat too short\n");
                  exit(-1);
                }

              /* Extract the chain */
              chain = input_line[0];

              /* Extract the start and end residues */
              strncpy(number_string,input_line+8,6);
              number_string[6] = '\0';
              min_res = atoi(number_string);
              strncpy(number_string,input_line+14,6);
              number_string[6] = '\0';
              max_res = atoi(number_string);

              /* Determine number of residues */
              nresid = max_res - min_res + 1;

              /* Create array for storing residue details */
              residue_array_ptr
                = (struct residue *) malloc(nresid * sizeof(struct residue));
              if (residue_array_ptr == NULL)
                {
                  printf("*** ERROR. Unable to allocate memory for "
                         "residue array\n");
                  exit(1);
                }

              /* Initialise all the residue records */
              for (ires = 0; ires < nresid; ires++)
                {
                  residue_array_ptr[ires].res_name[0] = '\0';
                  residue_array_ptr[ires].res_num[0] = '\0';
                  residue_array_ptr[ires].chain = '\0';
                  residue_array_ptr[ires].sec_struc = COIL;
                }
            }              

          /* Otherwise, check for HELIX or STRAND record */
          else if (!strncmp(input_line,"HELIX",5) ||
                   !strncmp(input_line,"STRAND",6))
            {
              /* Determine type of secondary structure element */
              if (!strncmp(input_line,"HELIX",5))
                type = HELIX;
              else
                type = STRAND;

              /* Extract the first and last sequence positions for this
                 secondary structure element (SSE) */
              strncpy(number_string,input_line+6,5);
              number_string[5] = '\0';
              first_pos = atoi(number_string) - min_res;
              strncpy(number_string,input_line+11,5);
              number_string[5] = '\0';
              last_pos = atoi(number_string) - min_res;

              /* Extract the column number (giving the x-position of
                 this SSE) */
              strncpy(number_string,input_line+16,4);
              number_string[4] = '\0';
              xpos[0] = xpos[1] = CONVERTX * atoi(number_string);

              /* Repeat for the to row numbers (giving the start and
                 end y-positions of this SSE) */
              strncpy(number_string,input_line+20,4);
              number_string[4] = '\0';
              ypos[0] = atoi(number_string);
              strncpy(number_string,input_line+24,4);
              number_string[4] = '\0';
              ypos[1] = atoi(number_string);

              /* Adjust the y-positions to expand this SSE to its true
                 length */
              ypos[0] = CONVERTY * ypos[0];
              ypos[1] = CONVERTY * ypos[1];
              if (ypos[0] < ypos[1])
                {
                  ypos[0] = ypos[0] - 1;
                  ypos[1] = ypos[1] + 1;
                }
              else
                {
                  ypos[0] = ypos[0] + 1;
                  ypos[1] = ypos[1] - 1;
                }

              /* Create a SSE record */
              sse_ptr =
                create_sse_record(&first_sse_ptr,&last_sse_ptr,type,
                                  first_pos,last_pos,xpos,ypos,FALSE);

              /* Store the start and end residue numbers */
              ires = sse_ptr->first_pos;
              strcpy(sse_ptr->first_resnum,residue_array_ptr[ires].res_num);
              remove_blanks(sse_ptr->first_resnum);
              ires = sse_ptr->last_pos;
              strcpy(sse_ptr->last_resnum,residue_array_ptr[ires].res_num);
              remove_blanks(sse_ptr->last_resnum);

              /* If this is the lowest-numbered SSE so far, place N-terminus
                 before it */
              if (first_posn == -9 || first_pos < first_posn)
                {
                  /* Store previous and current SSE */
                  pren_sse_ptr = previous_sse_ptr;
                  postn_sse_ptr = sse_ptr;

                  /* Store the first position */
                  first_posn = first_pos;
                }

              /* If this is the highest-numbered SSE so far, place
                 C-terminus after it */
              if (last_posc == -9 || last_pos > last_posc)
                {
                  /* Store previous and current SSE */
                  prec_sse_ptr = previous_sse_ptr;
                  postc_sse_ptr = sse_ptr;

                  /* Store the last position */
                  last_posc = last_pos;
                }

              /* Save the SSE pointer */
              previous_sse_ptr = sse_ptr;

              /* Increment count of records created */
              nsse++;
            }

          /* Otherwise, take to be a residue record */
          else if (rcount < nresid)
            {
              /* Extract residue details */
              inchn = input_line[0];
              strncpy(res_name,input_line+1,3);
              res_name[3] = '\0';
              strcpy(res_num,input_line+4);
              if (strlen(res_num) < 5)
                strcat(res_num," ");

              /* Store the details */
              strcpy(residue_array_ptr[rcount].res_name,res_name);
              strcpy(residue_array_ptr[rcount].res_num,res_num);
              residue_array_ptr[rcount].chain = chain;

              /* Increment residue count */
              rcount++;
            }
        }

      /* Close the file */
      fclose(file_ptr);
    }

  /* Show number of SSE records created */
  printf("Number of SSE records created: %d\n",nsse);

  /* Create SSE records corresponding to N- and C-termini */
  if (nsse > 0)
    nsse = create_termini(&first_sse_ptr,&last_sse_ptr,pren_sse_ptr,
                          postn_sse_ptr,prec_sse_ptr,postc_sse_ptr,
                          nsse,residue_array_ptr,nresid);

  /* Store the chain-id as one of the parameters */
  params->chain = chain;

  /* Return pointers to residue array and first SSE entry as well as
     number of residues */
  *res_array_ptr = residue_array_ptr;
  *fst_sse_ptr = first_sse_ptr;
  *lst_sse_ptr = last_sse_ptr;
  *nres = nresid;

  /* Return number of SSEs read in */
  return(nsse);
}
/***********************************************************************

separate_sses  -  Add extra space between SSEs in the same column

***********************************************************************/

void separate_sses(struct sse *first_sse_ptr)
{
  int nsse, strand_pos;
  int adjust, x;
  int endcol;

  struct sse *sse_ptr, *next_sse_ptr, *start_sse_ptr;

  /* Initialise pointer to first SSE record */
  next_sse_ptr = first_sse_ptr;

  /* Loop until we have processed all the columns */
  while (next_sse_ptr != NULL)
    {
      /* Get the next column's starting SSE */
      start_sse_ptr = sse_ptr = next_sse_ptr;

      /* Initialise variables */
      endcol = FALSE;
      nsse = 0;
      strand_pos = -1;

      /* Get the x-value of this column */
      x = sse_ptr->xpos[0];

      /* Loop through all the SSE records to find those belonging to
         the current column */
      while (sse_ptr != NULL && endcol == FALSE)
        {
          /* If this is the same x-position, then determine whether
             it is a strand */
          if (x == sse_ptr->xpos[0])
            {
              /* If no strand encountered yet and this is the first
                 one, store its position */
              if (sse_ptr->type == STRAND && strand_pos == -1)
                strand_pos = nsse;

              /* Increment count of secondary structure elements in this
                 column */
              nsse++;
            }

          /* Otherwise we are done */
          else
            {
              /* Set flag */
              endcol = TRUE;

              /* Store this as the next SSE to start at */
              next_sse_ptr = sse_ptr;
            }

          /* Get next SSE */
          sse_ptr = sse_ptr->next_sse_ptr;
          if (sse_ptr == NULL)
            next_sse_ptr = NULL;
        }

      /* If have more than one SSE in this column, adjust y-coords
         of all but the first strand */
      if (nsse > 1)
        {
          /* If no strand, fit to middle helix */
          if (strand_pos == -1)
            strand_pos = nsse / 2;

          /* Get the first of this column's starting SSE */
          sse_ptr = start_sse_ptr;
          nsse = 0;

          /* Loop through all the SSEs in the current column */
          while (sse_ptr != next_sse_ptr)
            {
              /* Compute adjustment to be made */
              adjust = 2 * (nsse - strand_pos);

              /* Adjust the y-coords of this SSE */
              sse_ptr->ypos[0] = sse_ptr->ypos[0] + adjust;
              sse_ptr->ypos[1] = sse_ptr->ypos[1] + adjust;

              /* Adjust the min and max y-coords of this SSE */
              sse_ptr->ymin = sse_ptr->ymin + adjust;
              sse_ptr->ymax = sse_ptr->ymax + adjust;

              /* Increment SSE count */
              nsse++;

              /* Get next SSE */
              sse_ptr = sse_ptr->next_sse_ptr;
            }
        }
    }
}
/***********************************************************************

layout_helices  -  Calculate the layout of the SSEs

***********************************************************************/

void layout_helices(struct sse *first_sse_ptr,int nsse,
                    struct helint *first_helint_ptr,int nhelint)
{
  struct helint *helint_ptr;
  struct sse *sse_ptr;

  /* Initialise variables */

  /* Get pointer to first helint record */
  helint_ptr = first_helint_ptr;

  /* Loop through helint records to write out */
  while (helint_ptr != NULL)
    {
/* debug */
      printf("Helix interaction:  ");
      if (helint_ptr->promotif1_ptr == NULL)
        printf("NULL ");
      else
        printf("Helix %d ",helint_ptr->promotif1_ptr->identifier[0]);
      if (helint_ptr->promotif2_ptr == NULL)
        printf("NULL ");
      else
        printf("Helix %d ",helint_ptr->promotif2_ptr->identifier[0]);
      printf("  Distance %8.3f\n",helint_ptr->distance);
      fflush(stdout);
/* debug */

      /* Get next helint record */
      helint_ptr = helint_ptr->next_helint_ptr;
    }

/* debug */
  printf("Number of interactions : %d\n",nhelint);
      fflush(stdout);

  exit(1);
/* debug */



}
/***********************************************************************

create_sort_index  -  Create the array holding the sort-indices for
                      the SSE records

***********************************************************************/

void create_sort_index(struct sse *first_sse_ptr,
                       struct sse ***sse_index_ptr,int nsse)
{
  int nentry;

  struct sse **sse_ind_ptr;
  struct sse *sse_ptr;

  /* Initialise variables */
  printf("Creating sse index ...\n");
  nentry = 0;

  /* Create the index array as one large array */
  sse_ind_ptr = (struct sse **) malloc(nsse * sizeof(struct sse *));

  /* Check that sufficient memory was allocated */
  if (sse_ind_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for sse-index "
             "array\n");
      exit(1);
    }

  /* Get pointer to the first sse entry */
  sse_ptr = first_sse_ptr;

  /* Loop through all entries to place each one in the sorted index
     array */
  while (sse_ptr != NULL)
    {
      /* Store current pointer at end of the index array */
      if (nentry < nsse)
        sse_ind_ptr[nentry] = sse_ptr;

      /* Increment count of stored entries */
      nentry++;

      /* Get pointer to the next sse in the linked list */
      sse_ptr = sse_ptr->next_sse_ptr;
    }

  /* Return the array pointer */
  *sse_index_ptr = sse_ind_ptr;
}
/***********************************************************************

sort_sses  -  Sort the SSEs according to their order in the sequence
              using a simple bubble-sort

***********************************************************************/

void sort_sses(struct sse **sse_index_ptr,int nsse)
{
  int isse;
  int swapped;

  struct sse *swap_sse_ptr;

  /* Initialise swap flag */
  swapped = TRUE;

  /* Keep looping while any elements are still being swapped */
  while (swapped == TRUE)
    {
      /* Unset the swap flag */
      swapped = FALSE;

      /* Loop through all the SSE codes */
      for (isse = 0; isse < nsse - 1; isse++)
        {
          /* If the two sses are the wrong way round, then need to
             swap them */
          if (sse_index_ptr[isse]->first_pos
              > sse_index_ptr[isse + 1]->first_pos ||
              (sse_index_ptr[isse]->first_pos ==
               sse_index_ptr[isse + 1]->first_pos &&
               sse_index_ptr[isse]->type > sse_index_ptr[isse + 1]->type))
            {
              /* Swap the corresponding pointers */
              swap_sse_ptr = sse_index_ptr[isse + 1];
              sse_index_ptr[isse + 1] = sse_index_ptr[isse];
              sse_index_ptr[isse] = swap_sse_ptr;

              /* Set flag to indicate that a swap took place in
                 this loop */
              swapped = TRUE;
            }
        }
    }
}
/***********************************************************************

create_mask_array  -  Create the array storing locations of all SSEs
                      for use in placing the joins

***********************************************************************/

void create_mask_array(int **arr,int *npts,int *origin,
                       struct limits limit,struct sse **sse_index_ptr,
                       int nsse)
{
  int i, iend, ipos, isse, istart, j, jend, jstart, size;
  int valid;

  int *array;

  struct sse *sse_ptr;

  /* Determine size of array */
  npts[0] = limit.max_x - limit.min_x + 1 + 2 * ARRAY_MARGIN;
  npts[1] = limit.max_y - limit.min_y + 1 + 2 * ARRAY_MARGIN;
  size = npts[0] * npts[1];

  /* Allocate memory for the array */
  array = (int *) malloc(size * sizeof(int));

  /* Initialise array */
  for (i = 0; i < size; i++)
    array[i] = 0;

  /* Determine coords of array's origin */
  origin[0] = limit.min_x - ARRAY_MARGIN;
  origin[1] = limit.min_y - ARRAY_MARGIN;

  /* Loop through all the SSE codes to update mask array */
  for (isse = 0; isse < nsse; isse++)
    {
      /* Get this SSE */
      sse_ptr = sse_index_ptr[isse];

      /* Determine start- and end-points of this SSE */
      istart = sse_ptr->xmin - origin[0];
      iend = sse_ptr->xmax - origin[0];
      jstart = sse_ptr->ymin - origin[1];
      jend = sse_ptr->ymax - origin[1];

      /* Check that limits are OK */
      valid = TRUE;
      if (istart < 0 || istart > npts[0] - 1 ||
          iend < 0 || iend > npts[0] - 1 ||
          jstart < 0 || jstart > npts[1] - 1 ||
          jend < 0 || jend > npts[1] - 1)
        {
          printf("*** ERROR ... SSE off end of array!\n");
          printf("  Drawing limits:  x = %6.1f -> %6.1f;  "
                 "y = %6.1f -> %6.1f\n",
                 limit.min_x,limit.max_x,limit.min_y,limit.max_y);
          printf("  SSE limits:      x = %d -> %d;  y = %d -> %d\n",
                 sse_ptr->xmin,sse_ptr->xmax,sse_ptr->ymin,sse_ptr->ymax);
          printf("  Net SSE limits:  x = %d -> %d;  y = %d -> %d\n",
                 istart,iend,jstart,jend);
          valid = FALSE;
        }

      /* Loop over the array points occupied by this SSE to update */
      if (valid == TRUE)
        {
          for (i = istart; i <= iend; i++)
            for (j = jstart; j <= jend; j++)
              {
                /* Compute position of this point in the array */
                ipos = i * npts[1] + j;

                /* Set array value to 1 */
                if (ipos > -1 && ipos < size)
                  array[ipos] = SSE_SCORE;
              }
        }
    }

  /* Return pointer to the start of the array */
  *arr = array;
}
/***********************************************************************

check_crossover  -  Check whether the given line crosses any SSEs

***********************************************************************/

int check_crossover(int x1,int y1,int x2,int y2,int *array,int *npts,
                    int *origin)
{
  int i, iend, ipos, istart, j, jend, jstart;
  int score;

  /* Initialise variables */
  score = 0;

  /* Determine start and end of current line */
  if (x1 < x2)
    {
      istart = x1 - origin[0];
      iend = x2 - origin[0];
    }
  else
    {
      istart = x2 - origin[0];
      iend = x1 - origin[0];
    }
  if (y1 < y2)
    {
      jstart = y1 - origin[1];
      jend = y2 - origin[1];
    }
  else
    {
      jstart = y2 - origin[1];
      jend = y1 - origin[1];
    }

  /* If line lies within the bounds of the array, check whether it
     crosses any SSEs or other lines */
  if (istart < 0 || istart > npts[0] - 1 ||
      iend < 0 || iend > npts[0] - 1 ||
      jstart < 0 || jstart > npts[1] - 1 ||
      jend < 0 || jend > npts[1] - 1)
    score = 10;

  /* If line lies within the array, check for cross-overs */
  if (score == 0)
    {
      /* Loop over all points along the given line */
      for (i = istart; i <= iend; i++)
        {
          for (j = jstart; j <= jend; j++)
            {
              /* Calculate position in array */
              ipos = i * npts[1] + j;

              /* Add score to current score */
              score = score + array[ipos];
            }
        }
    }

  /* Return cross-over score */
  return(score);
}
/***********************************************************************

find_path1  -  Check whether a single-segment path between the two
               endpoints is OK

***********************************************************************/

int find_path1(int xstart,int ystart,int xend,int yend,
               int x[MAX_SEG][2],int y[MAX_SEG][2],int *array,int *npts,
               int *origin)
{
  int iseg, nseg, score;

  /* Initialise flag */
  iseg = 0;
  nseg = 0;
  score = 50;

  /* Check whether end can be joined by a single horizontal or vertical
     line */
  if (xstart == xend || ystart == yend)
    {
      /* Define segment coordinates */
      x[iseg][0] = xstart;
      x[iseg][1] = xend;
      y[iseg][0] = ystart;
      y[iseg][1] = yend;

      /* Check whether this segment crosses any SSEs */
      score = check_crossover(x[iseg][0],y[iseg][0],x[iseg][1],
                              y[iseg][1],array,npts,origin);

      /* If all segments are clear, can create the SSE records */
      if (score == 0)
        nseg = 1;
    }

  /* Return number of segments created */
  return(nseg);
}
/***********************************************************************

find_path2  -  Check whether a two-segment path between the two
               endpoints is OK

***********************************************************************/

int find_path2(int xstart,int ystart,int xend,int yend,
               int x[MAX_SEG][2],int y[MAX_SEG][2],int *array,int *npts,
               int *origin)
{
  int iseg, loop, nseg, dx, dy, score;
  int done;

  /* Initialise flag */
  done = FALSE;
  nseg = 0;
  score = 50;

  /* Calculate step-size in each direction */
  dx = xend - xstart;
  dy = yend - ystart;

  /* Loop over the two possibilities (along and up/down, or up/down
     and along) */
  for (loop = 0; loop < 2 && done == FALSE; loop++)
    {
      /* Re-initialise score */
      score = 0;

      /* Loop over this path's 2 segments */
      for (iseg = 0; iseg < 2 && score == 0; iseg++)
        {
          /* Define segment coordinates */
          if (iseg == 0)
            {
              x[iseg][0] = xstart;
              y[iseg][0] = ystart;
              x[iseg][1] = xstart + loop * dx;
              y[iseg][1] = ystart + (1 - loop) * dy;
            }
          else
            {
              x[iseg][0] = xstart + loop * dx;
              y[iseg][0] = ystart + (1 - loop) * dy;
              x[iseg][1] = xend;
              y[iseg][1] = yend;
            }

          /* Check whether this segment crosses any SSEs */
          score = check_crossover(x[iseg][0],y[iseg][0],x[iseg][1],
                                  y[iseg][1],array,npts,origin);
        }

      /* If score zero, then we're done */
      if (score == 0)
        done = TRUE;
    }

  /* If all segments are clear, can create the SSE records */
  if (score == 0)
    nseg = 2;

  /* Return number of segments created */
  return(nseg);
}
/***********************************************************************

find_path3  -  Check whether a three-segment path between the two
               endpoint is OK

***********************************************************************/

int find_path3(int xstart,int ystart,int xend,int yend,
               int x[MAX_SEG][2],int y[MAX_SEG][2],int *array,int *npts,
               int *origin)
{
  int iseg, loop, nseg, score;
  int done;

  /* Initialise flag */
  done = FALSE;
  nseg = 0;
  score = 50;

  /* Loop over the two possibilities (along and up/down, or up/down
     and along) */
  for (loop = 0; loop < 2 && done == FALSE; loop++)
    {
      /* Re-initialise score */
      score = 0;

      /* Loop over the path's 3 segments */
      for (iseg = 0; iseg < 3 && score == 0; iseg++)
        {
          /* Define segment coordinates */
          if (loop == 0)
            {
              x[iseg][0] = xstart;
              if (iseg > 0)
                x[iseg][0] = (xstart + xend) / 2;
              x[iseg][1] = (xstart + xend) / 2;
              if (iseg > 1)
                x[iseg][1] = xend;
              y[iseg][0] = ystart;
              if (iseg > 1)
                y[iseg][0] = yend;
              y[iseg][1] = ystart;
              if (iseg > 0)
                y[iseg][1] = yend;
            }
          else
            {
              x[iseg][0] = xstart;
              if (iseg > 1)
                x[iseg][0] = xend;
              x[iseg][1] = xstart;
              if (iseg > 0)
                x[iseg][1] = xend;
              y[iseg][0] = ystart;
              if (iseg > 0)
                y[iseg][0] = (ystart + yend) / 2;
              y[iseg][1] = (ystart + yend) / 2;
              if (iseg > 1)
                y[iseg][1] = yend;
            }

          /* Check whether this segment crosses any SSEs */
          score = check_crossover(x[iseg][0],y[iseg][0],x[iseg][1],
                                  y[iseg][1],array,npts,origin);
        }

      /* If score still zero, then we are done */
      if (score == 0)
        done = TRUE;
    }

  /* If all segments are clear, can create the SSE records */
  if (score == 0)
    nseg = 3;

  /* Return number of segments created */
  return(nseg);
}
/***********************************************************************

try_paths  -  Try all possible paths between the two given end-points

***********************************************************************/

void try_paths(int xstart,int ystart,int xend,int yend,
               int nseg,int step[MAX_SEG][2],int max_step,
               int xseg[MAX_SEG][2],int yseg[MAX_SEG][2],int *best_path,
               int *best_score,int *best_nseg,int *array,int *npts,
               int *origin,int path_score)
{
  int dir, dir_end, i, j;
  int curr_path, path_len, start_size, end_size, score, size;
  int diff, x1, x2, xstep, y1, y2, ystep;
  int x[MAX_SEG][2], y[MAX_SEG][2];
  int at_end, valid;

  /* Initialise variables */
  at_end = FALSE;
  path_len = 0;
  x1 = xstart;
  y1 = ystart;

  /* Add all segments up to current point in path */
  for (i = 0; i < nseg; i++)
    {
      /* Store segment's start-point */
      x[i][0] = x1;
      y[i][0] = y1;

      /* Get current size and direction */
      xstep = step[i][STEP_SIZE] * DIRECTION[step[i][STEP_DIR]][0];
      ystep = step[i][STEP_SIZE] * DIRECTION[step[i][STEP_DIR]][1];
      x1 = x1 + xstep;
      y1 = y1 + ystep;

      /* Store segment's end-point */
      x[i][1] = x1;
      y[i][1] = y1;

      /* Increment the current path-length */
      path_len = path_len + abs(xstep) + abs(ystep);
    }

  /* Loop over the 4 directions for this segment */
  for (dir = 0; dir < 4 && at_end == FALSE; dir++)
    {
      /* If this direction is the same as that of the previous segment,
         then skip it */
      if (nseg == 0 || (step[nseg - 1][STEP_DIR] != dir &&
                        DIR_AXIS[step[nseg - 1][STEP_DIR]] != DIR_AXIS[dir]))
        {
          /* Store current direction */
          step[nseg][STEP_DIR] = dir;

          /* Set default number of trials steps for this segment */
          start_size = 1;
          end_size = max_step;

          /* If can reach endpoint in a single step, then set to single
             loop */
          x2 = x1 + abs(xend - x1) * DIRECTION[dir][0];
          y2 = y1 + abs(yend - y1) * DIRECTION[dir][1];
          if (x2 == xend && y2 == yend)
            {
              /* Set the step-size to get to the end-point in a single
                 step */
              start_size = abs(xend - x1) + abs(yend - y1);
              end_size = start_size;

              /* Set flag to say that this segment takes us to the
                 end */
              at_end = TRUE;
            }

          /* If this is the final permissible segment and we can't reach
             the end by stepping in this direction, then no point trying */
          else if (nseg >= MAX_SEG - 1)
            end_size = 0;

          /* Set flag */
          valid = TRUE;

          /* Loop over the trial steps in the current direction */
          for (size = start_size; size <= end_size && valid == TRUE;
               size++)
            {
              /* Store the step size */
              step[nseg][STEP_SIZE] = size;

              /* Calculate end-point of segment */
              x2 = x1 + size * DIRECTION[dir][0];
              y2 = y1 + size * DIRECTION[dir][1];

              /* Check whether this segment is a valid one */
              score = check_crossover(x1,y1,x2,y2,array,npts,origin);

              /* Add to current score */
              score = path_score + score;

              /* If valid, then continue */
              if (score < MAX_SCORE + 1)
                {
                  /* If at end of path, see if this is the shortest path
                     so far */
                  if (at_end == TRUE)
                    {
                      /* Calculate total path length */
                      curr_path = path_len + size;

                      /* Add number of segments to score to bias towards
                         paths having fewer segments */
                      score = score + nseg;

                      /* If this is the first, or shortest path, store
                         as best path so far */
                      if (*best_path == 0 || score < *best_score ||
                          (score == *best_score &&
                           curr_path < *best_path))
                        {
                          /* Store current path length and its number
                             of segments and score */
                          *best_path = curr_path;
                          *best_nseg = nseg + 1;
                          *best_score = score;

                          /* Store last segment's start-point */
                          x[nseg][0] = x1;
                          y[nseg][0] = y1;

                          /* Store its end-point */
                          x[nseg][1] = xend;
                          y[nseg][1] = yend;

                          /* Save the segment coords */
                          for (i = 0; i < nseg + 1; i++)
                            {
                              xseg[i][0] = x[i][0];
                              yseg[i][0] = y[i][0];
                              xseg[i][1] = x[i][1];
                              yseg[i][1] = y[i][1];
                            }
                         }
                    }

                  /* If not at the end of the path, add another segment */
                  else if (nseg < MAX_SEG - 1)
                    try_paths(xstart,ystart,xend,yend,nseg + 1,step,
                              max_step,xseg,yseg,best_path,best_score,
                              best_nseg,array,npts,origin,score);
                }

              /* If current path invalid, then set flag */
              else
                valid = FALSE;
            }
        }
    }

  return;
}
/***********************************************************************

update_mask_array  -  Update the array with the line segment just
                      created

***********************************************************************/

void update_mask_array(int *x,int *y,int *array,int *npts,int *origin)
{
  int i, iend, ipos, istart, j, jend, jstart;

  /* Determine start and end of current line */
  if (x[0] < x[1])
    {
      istart = x[0] - origin[0];
      iend = x[1] - origin[0];
    }
  else
    {
      istart = x[1] - origin[0];
      iend = x[0] - origin[0];
    }
  if (y[0] < y[1])
    {
      jstart = y[0] - origin[1];
      jend = y[1] - origin[1];
    }
  else
    {
      jstart = y[1] - origin[1];
      jend = y[0] - origin[1];
    }

  /* If line lies outside bounds, then return */
  if (istart < 0 || istart > npts[0] - 1 ||
      iend < 0 || iend > npts[0] - 1 ||
      jstart < 0 || jstart > npts[1] - 1 ||
      jend < 0 || jend > npts[1] - 1)
    return;

  /* Loop over all points along the given line to update array */
  for (i = istart; i <= iend; i++)
    {
      for (j = jstart; j <= jend; j++)
        {
          /* Calculate position in array */
          ipos = i * npts[1] + j;

          /* Update mask array */
          array[ipos] = LINE_SCORE;
        }
    }
}
/***********************************************************************

add_joins  -  Add straight-line segment between consecutive SSEs

***********************************************************************/

int add_joins(struct sse **sse_index_ptr,int nsse,
              struct sse **fst_sse_ptr,struct sse **lst_sse_ptr,
              struct residue *residue_array_ptr,int nres,
              struct limits limit,int *array,int *npts,int *origin)
{
  int ires, iseg, isse, jres, nseg, tot_sse;
  int first_pos, last_pos, x[2], y[2], xover[2], yover[2];
  int max_step, step[MAX_SEG][2], xseg[MAX_SEG][2], yseg[MAX_SEG][2];
  int best_path, best_score;
  int endseg;

  struct sse *sse_ptr, *new_sse_ptr, *next_sse_ptr;
  struct sse *first_sse_ptr, *last_sse_ptr;

  /* Initialise variables */
  tot_sse = nsse;

  /* Initialise sse pointers */
  first_sse_ptr = *fst_sse_ptr;
  last_sse_ptr = *lst_sse_ptr;

  /* Loop through all the SSE codes */
  for (isse = 0; isse < nsse - 1; isse++)
    {
      /* Initialise variables */
      best_path = 0;
      best_score = MAX_SCORE + 1;

      /* Get this SSE */
      sse_ptr = sse_index_ptr[isse];

      /* Get the next one */
      next_sse_ptr = sse_index_ptr[isse + 1];

      /* Get coordinates between which we need join */
      x[0] = sse_ptr->xpos[1];
      y[0] = sse_ptr->ypos[1];
      x[1] = next_sse_ptr->xpos[0];
      y[1] = next_sse_ptr->ypos[0];

      /* Determine largest step size to try */
      max_step = abs(x[1] - x[0]);
      if (abs(y[1] - y[0]) > max_step)
        max_step = abs(y[1] - y[0]);
      max_step = 2 * max_step;

      /* Store the coil's first and last residues */
      first_pos = sse_ptr->last_pos;
      last_pos = next_sse_ptr->first_pos;

      /* Adjust, depending on orientation and direction of SSEs */
      if (sse_ptr->orientation == HORIZONTAL)
        {
          if (sse_ptr->direction == FORWARD)
            x[0] = x[0] + COIL_OVERHANG;
          else
            x[0] = x[0] - COIL_OVERHANG;
        }
      else
        {
          if (sse_ptr->direction == FORWARD)
            y[0] = y[0] + COIL_OVERHANG;
          else
            y[0] = y[0] - COIL_OVERHANG;
        }
      if (next_sse_ptr->orientation == HORIZONTAL)
        {
          if (next_sse_ptr->direction == FORWARD)
            x[1] = x[1] - COIL_OVERHANG;
          else
            x[1] = x[1] + COIL_OVERHANG;
        }
      else
        {
          if (next_sse_ptr->direction == FORWARD)
            y[1] = y[1] - COIL_OVERHANG;
          else
            y[1] = y[1] + COIL_OVERHANG;
        }

      /* Save the overhang coords */
      xover[0] = sse_ptr->xpos[1];
      yover[0] = sse_ptr->ypos[1];
      xover[1] = x[0];
      yover[1] = y[0];

      /* Create the first overhang segment */
      endseg = FALSE;
      new_sse_ptr
        = create_sse_record(&first_sse_ptr,&last_sse_ptr,OVERHANG1,
                            first_pos + 1,last_pos - 1,xover,
                            yover,endseg);

      /* Store the sort position */
      new_sse_ptr->sort = first_pos;

      /* Get the start and end residue numbers for this coil region */
      ires = first_pos + 1;
      jres = last_pos - 1;

      /* If this coil region ends on the C-terminus, adjust the end
         residue */
      if (next_sse_ptr->type == CTERMINUS)
        jres = nres - 1;

      /* If this is the first segment, get the coil region's 
         start and end residues */
      if (ires > -1 && ires < nres && jres > -1 && jres < nres &&
          ires <= jres)
        {
          strcpy(new_sse_ptr->first_resnum,
                 residue_array_ptr[ires].res_num);
          strcpy(new_sse_ptr->last_resnum,
                 residue_array_ptr[jres].res_num);
        }
      else
        {
          strcpy(new_sse_ptr->first_resnum,"   0 ");
          strcpy(new_sse_ptr->last_resnum,"   0 ");
        }

      /* Remove blanks from the residue numbers */
      remove_blanks(new_sse_ptr->first_resnum);
      remove_blanks(new_sse_ptr->last_resnum);

      /* See if can join the end with a single line */
      nseg = find_path1(x[0],y[0],x[1],y[1],xseg,yseg,array,npts,
                        origin);

      /* Identify a possible 2-segment path */
      if (nseg == 0)
        nseg = find_path2(x[0],y[0],x[1],y[1],xseg,yseg,array,npts,
                          origin);

      /* If not found, identify a possible 3-segment path */
      if (nseg == 0)
        nseg = find_path3(x[0],y[0],x[1],y[1],xseg,yseg,array,npts,
                          origin);

      /* Finally, try a general search using multi-segment paths */
      if (nseg == 0)
        try_paths(x[0],y[0],x[1],y[1],0,step,max_step,xseg,yseg,&best_path,
                  &best_score,&nseg,array,npts,origin,0);

      /* If have a path, create the required number of line SSEs */
      if (nseg > 0)
        {
          /* Loop over the line segments */
          for (iseg = 0; iseg < nseg; iseg++)
            {
              /* Determine whether this is an end segment */
              if (nseg != 2 && (iseg == 0 || iseg == nseg - 1))
                endseg = TRUE;
              else
                endseg = FALSE;

              /* Create a line SSE */
              new_sse_ptr
                = create_sse_record(&first_sse_ptr,&last_sse_ptr,COIL,
                                    first_pos,last_pos,xseg[iseg],
                                    yseg[iseg],endseg);

              /* Store the sort position */
              new_sse_ptr->sort = first_pos;

              /* Update the mask array */
              update_mask_array(xseg[iseg],yseg[iseg],array,npts,
                                origin);

              /* Increment count of records created */
              tot_sse++;
            }
        }

      /* Save the end overhang coords */
      xover[0] = x[1];
      yover[0] = y[1];
      xover[1] = next_sse_ptr->xpos[0];
      yover[1] = next_sse_ptr->ypos[0];

      /* Create the last overhang segment */
      endseg = TRUE;
      new_sse_ptr
        = create_sse_record(&first_sse_ptr,&last_sse_ptr,OVERHANG2,
                            first_pos + 1,last_pos - 1,xover,
                            yover,endseg);

      /* Store the sort position */
      new_sse_ptr->sort = first_pos;
    }

  /* Return current pointers */
  *fst_sse_ptr = first_sse_ptr;
  *lst_sse_ptr = last_sse_ptr;

  /* Return number of SSEs */
  return(tot_sse);
}
/***********************************************************************

get_max_min  -  Get the maximum and minimum x-, and y-coords of objects
                in the picture

***********************************************************************/

void get_max_min(struct sse *first_sse_ptr,struct limits *limit)
{
  int i;

  float xmin, xmax, ymin, ymax;

  struct sse *sse_ptr;

  /* Initialise pointer to first SSE record */
  sse_ptr = first_sse_ptr;

  /* Loop through all the SSE records to update coordinate limits */
  while (sse_ptr != NULL)
    {
      /* Loop over both ends of this SSE */
      for (i = 0; i < 2; i++)
        {
          /* Get the coordinates */
          xmin = sse_ptr->xmin;
          ymin = sse_ptr->ymin;
          xmax = sse_ptr->xmax;
          ymax = sse_ptr->ymax;

          /* If this is the first SSE encountered, store its
             coords as the present max and min values */
          if (limit->first == TRUE)
            {
              limit->min_x = xmin;
              limit->max_x = xmax;
              limit->min_y = ymin;
              limit->max_y = ymax;
              limit->first = FALSE;
            }

          /* Otherwise, update the limits if coordinates are
             outside them */
          else
            {
              if (xmin < limit->min_x)
                limit->min_x = xmin;
              if (ymin < limit->min_y)
                limit->min_y = ymin;
              if (xmax > limit->max_x)
                limit->max_x = xmax;
              if (ymax > limit->max_y)
                limit->max_y = ymax;
            }
        }

      /* Get pointer to next sse in linked-list */
      sse_ptr = sse_ptr->next_sse_ptr;
    }
}
/***********************************************************************

adjust_for_landscape  -  Perform boundary adjustments for printing
                         in Landscape, as opposed to Portrait,
                         orientation

***********************************************************************/

void adjust_for_landscape(struct pspage *pspage_ptr)
{
  /* Make necessary adjustments */
  pspage_ptr->Page_Min_x = BBOXY1 + XMARGIN;
  pspage_ptr->Page_Max_x = pspage_ptr->Page_Min_x + YHEIGHT;
  pspage_ptr->Page_Min_y = BBOXX1 + YMARGIN;
  pspage_ptr->Page_Max_y = pspage_ptr->Page_Min_y + YHEIGHT_LAND;
}
/***********************************************************************

initialise_psparams  -  Initialise PostScript parameters

***********************************************************************/

void initialise_psparams(struct psparams *psparams_ptr,
                         struct pspage *pspage_ptr)
{
  /* Initialise PostScript parameters */
  psparams_ptr->landscape = FALSE;
  psparams_ptr->splitps = FALSE;
  psparams_ptr->in_colour = TRUE;
  strcpy(psparams_ptr->program_name,"topoldia");

  /* Initialise PostScript plot variables */
  pspage_ptr->ps_name[0] = '\0';
  pspage_ptr->ps_file = NULL;
  pspage_ptr->tdf_file = NULL;
  pspage_ptr->Page_Min_x = BBOXX1 + XMARGIN;
  pspage_ptr->Page_Max_x = pspage_ptr->Page_Min_x + XWIDTH;
  pspage_ptr->Page_Min_y = BBOXY1 + YMARGIN;
  pspage_ptr->Page_Max_y = pspage_ptr->Page_Min_y + YHEIGHT;
}
/***********************************************************************

init_pspage  -  Initialise PostScript page variables

***********************************************************************/

void init_pspage(struct psparams psparams_ptr,struct pspage *pspage_ptr)
{
  /* If PostScript plot to be in landscape mode, adjust relevant
     page parameters */
  if (psparams_ptr.landscape == TRUE)
    adjust_for_landscape(pspage_ptr);

  /* Initialise page variables */
  pspage_ptr->page = 0;

  /* Initialise starting y-position */
  pspage_ptr->y = pspage_ptr->Page_Max_y;

  /* Initialise page plot limits */
  pspage_ptr->lim[0] = XMARGIN + XWIDTH;
  pspage_ptr->lim[1] = 0.0;
  pspage_ptr->lim[2] = YMARGIN + YHEIGHT;
  pspage_ptr->lim[3] = 0.0;
}
/***********************************************************************

psbtxt_   Write out bold centred text to PostScript file

***********************************************************************/

void psbtxt_(FILE *ps_file,float x,float y,float size,char *text_string)
{
  fprintf(ps_file,"%7.2f %7.2f moveto\n", x, y);
  fprintf(ps_file,"(%s) %7.2f Center\n",text_string,size);
  fprintf(ps_file,"(%s) %7.2f Bprint\n",text_string,size);
}
/***********************************************************************

psccol_  -  Define colour for filled circle

***********************************************************************/

void psccol_(FILE *ps_file,char *text_string)
{
    fprintf(ps_file,"/Circol { %s } def\n",text_string);
}
/***********************************************************************

pscolb_  -  Write background colour level (for text, etc)

***********************************************************************/

void pscolb_(FILE *ps_file,char *text_string)
{
    fprintf(ps_file,"%s\n",text_string);
}
/***********************************************************************

pscomm_   Write out a comment to the PostScript file

***********************************************************************/

void pscomm_(FILE *ps_file,char *text_string)
{
  fprintf(ps_file,"\n");
  fprintf(ps_file,"%% %s\n",text_string);
  fprintf(ps_file,"\n");
}
/***********************************************************************

psctxt_   Write out centred text to PostScript file

***********************************************************************/

void psctxt_(FILE *ps_file,float x,float y,float size,char *text_string)
{
  fprintf(ps_file,"%7.2f %7.2f moveto\n", x, y);
  fprintf(ps_file,"(%s) %7.2f Center\n",text_string,size);
  fprintf(ps_file,"(%s) %7.2f Print\n",text_string,size);
}
/***********************************************************************

psearc_  -  Write out an elliptical arc

***********************************************************************/

void psearc_(FILE *ps_file,float x,float y,float radius,
             double angle_start,double angle_end,float elong_x,
             float elong_y)
{
    fprintf(ps_file,"0.0 0.0 %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f "
            "Ellarc\n",radius,angle_start,angle_end,elong_x,elong_y,x,y);
}
/***********************************************************************

pselip_  -  Write out a coloured ellipse to PostScript file

***********************************************************************/

void pselip_(FILE *ps_file,float x,float y,float radius,float elong_x,
             float elong_y)
{
    fprintf(ps_file,"0.0 0.0 %7.2f %7.2f %7.2f %7.2f %7.2f Ellipse\n",
            radius,elong_x,elong_y,x,y);
}
/***********************************************************************

pshade_  -  Write colour/shade out to PostScript file

***********************************************************************/

void pshade_(FILE *ps_file,char *colour)
{
    fprintf(ps_file,"G %s\n",colour);
}
/***********************************************************************

pshcol_  -  Write colour/shade as RGB values out to PostScript file

***********************************************************************/

void pshcol_(FILE *ps_file,float rgb[3])
{
    fprintf(ps_file,"G %8.3f %8.3f %8.3f setrgbcolor\n",rgb[0],rgb[1],
            rgb[2]);
}
/***********************************************************************

psitxt_   Write out centred, italic text to PostScript file

***********************************************************************/

void psitxt_(FILE *ps_file,float x,float y,float size,char *text_string)
{
  fprintf(ps_file,"%7.2f %7.2f moveto\n", x, y);
  fprintf(ps_file,"(%s) %7.2f Center\n",text_string,size);
  fprintf(ps_file,"(%s) %7.2f Iprint\n",text_string,size);
}
/***********************************************************************

psline_  -   Write line out to PostScript file

***********************************************************************/

void psline_(FILE *ps_file,float x1,float y1,float x2,float y2)
{
  fprintf(ps_file,"%7.2f %7.2f %7.2f %7.2f L\n",x1, y1, x2, y2);
}
/***********************************************************************

pslwid_  -  Write line-width out to PostScript file

***********************************************************************/

void pslwid_(FILE *ps_file,float width)
{
    fprintf(ps_file,"%8.3f W\n",width);
}
/***********************************************************************

psrot_  -  Rotate graphics page through 90 degrees

***********************************************************************/

void psrot_(FILE *ps_file,float new_origin_x, float new_origin_y)
{
  fprintf(ps_file," %7.2f %7.2f moveto Rot90\n",new_origin_x,
          new_origin_y);
}
/***********************************************************************

psubox_  -  Write out unbounded box to PostScript file

***********************************************************************/

void psubox_(FILE *ps_file,float x1,float y1,float x2,float y2,float x3,
             float y3,float x4,float y4)
{
  fprintf(ps_file,"%7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f Pl4\n",
          x1, y1, x2, y2, x3, y3, x4, y4);
}
/***********************************************************************

psutri_  -  Write out unbounded triangle to PostScript file

***********************************************************************/

void psutri_(FILE *ps_file,float x1,float y1,float x2,float y2,float x3,
             float y3)
{
  fprintf(ps_file,"%7.2f %7.2f %7.2f %7.2f %7.2f %7.2f Pl3\n",
          x1, y1, x2, y2, x3, y3);
}
/***********************************************************************

get_col_rgb  -  Get the colour name for the given colour number and
                   the corresponding RGB values

***********************************************************************/

void get_col_rgb(int colour,char *c_name,float *rgb)
{
#define NCOLS   46

  char number_string[20];

  char *string_ptr;

  int number, ipos, len;

  static char *colour_defn[NCOLS + 1][3] = {
    "Blue",            "#0000FF", "0 0 255",     //  0
    "Red",             "#FF0000", "255 0 0",     //  1
    "Gold",            "#FFD700", "255 215 0",   //  2
    "Green",           "#00FF00", "0 255 0",     //  3
    "Purple",          "#A020F0", "160 32 240",  //  4
    "Cyan",            "#00FFFF", "0 255 255",   //  5
    "Firebrick",       "#B22222", "178 34 34",   //  6
    "LimeGreen",       "#32CD32", "50 205 50",   //  7
    "DeepSkyBlue",     "#00BFFF", "0 191 255",   //  8
    "HotPink",         "#FF69B4", "255 105 180", //  9
    "LightGreen",      "#90EE90", "144 238 144", // 10
    "LightGoldenrod1", "#FFEC8B", "255 236 139", // 11
    "SkyBlue",         "#87CEEB", "135 206 235", // 12
    "OrangeRed",       "#FF4500", "255 69 0",    // 13
    "SeaGreen",        "#2E8B57", "46 139 87",   // 14
    "MidnightBlue",    "#191970", "25 25 112",   // 15
    "DarkRed",         "#8B0000", "139 0 0",     // 16
    "Peru",            "#CD853F", "205 133 63",  // 17
    "SpringGreen",     "#00FF7F", "0 255 127",   // 18
    "Magenta",         "#FF00FF", "255 0 255",   // 19
    "DimGrey",         "#696969", "105 105 105", // 20
    "DarkGoldenrod",   "#B8860B", "184 134 11",  // 21
    "Maroon",          "#B03060", "176 48 96",   // 22
    "Orange",          "#FFA500", "255 165 0",   // 23
    "OliveDrab",       "#6B8E23", "107 142 35",  // 24
    "CornflowerBlue",  "#6495ED", "100 149 237", // 25
    "DeepPink",        "#FF1493", "255 20 147",  // 26
    "DarkSlateBlue",   "#483D8B", "72 61 139",   // 27
    "SandyBrown",      "#F4A460", "244 164 96",  // 28
    "DarkKhaki",       "#BDB76B", "189 183 107", // 29
    "Tomato",          "#FF6347", "255 99 71",   // 30
    "SteelBlue",       "#4682B4", "70 130 180",  // 31
    "Chocolate",       "#D2691E", "210 105 30",  // 32
    "PaleGreen",       "#98FB98", "152 251 152", // 33
    "MediumOrchid",    "#BA55D3", "186 85 211",  // 34
    "DarkSlateGray",   "#2F4F4F", "47 79 79",    // 35
    "Salmon",          "#FA8072", "250 128 114", // 36
    "SaddleBrown",     "#8B4513", "139 69 19",   // 37
    "SlateGrey",       "#708090", "112 128 144", // 38
    "LawnGreen",       "#7CFC00", "124 252 0",   // 39
    "Burlywood",       "#DEB887", "222 184 135", // 40
    "DarkGreen",       "#006400", "0 100 0",     // 41
    "SlateBlue",       "#6A5ACD", "106 90 205",  // 42
    "DarkOliveGreen",  "#556B2F", "85 107 47",   // 43
    "Grey",            "#BEBEBE", "190 190 190", // 44
    "Black",           "#000000", "0 0 0",       // 45
    "White",           "#FFFFFF", "255 255 255"  // 46
  };

  /* If colour invalid, set to white */
  if (colour < 0 || colour > NCOLS)
    colour = NCOLS;

  /* Get the colour name and corresponding HTML string */
  strcpy(c_name,colour_defn[colour][0]);

  /* Initialise for RGB conversion */
  len = 0;

  /* Get the RGB triplet for this colour */
  for (ipos = 0; ipos < 3; ipos++)
    {
      /* Extract the next number */
      strcpy(number_string,colour_defn[colour][2]+len);
      string_ptr = strstr(number_string," ");
      if (string_ptr != NULL)
        string_ptr[0] = '\0';
      len = len + strlen(number_string) + 1;
      number = atoi(number_string);

      /* Convert to float between 0.0 and 1.0 */
      rgb[ipos] = number / 255.0;
    }
}
/***********************************************************************

psopen_  -  Open PostScript file and write out header records

***********************************************************************/

FILE *psopen_(char *filename,char *program_name)
{
  char  c_name[COLNAM_LEN], colno[3], exclamation, percent;
  float grey_shade, rgb[3];
  int   icolour, igrey;

  FILE *ps_file;

  /* Initialise variables */
  percent = '%';
  exclamation = '!';

  /* Open the PostScript file */
  if ((ps_file = fopen(filename,"w")) == NULL)
    {
      printf("\n*** Unable to open PostScript file: %s\n",filename);
      exit(1);
    }

  /* Write out headings to PostScript file*/
  fprintf(ps_file,"%c%cPS-Adobe-3.0\n",percent,exclamation);
  fprintf(ps_file,"%%%%Creator: %s\n",program_name);
  fprintf(ps_file,"%%%%DocumentNeededResources: font Times-Roman Symbol\n");
  fprintf(ps_file,"%%%%BoundingBox: (atend)\n");
  fprintf(ps_file,"%%%%EndComments\n");
  fprintf(ps_file,"%%%%BeginProlog\n");
  fprintf(ps_file,"/L { moveto lineto stroke } bind def\n");
  fprintf(ps_file,"/G { gsave } bind def\n");
  fprintf(ps_file,"/W { setlinewidth } bind def\n");
  fprintf(ps_file,"/D { setdash } bind def\n");
  fprintf(ps_file,"/Col { setrgbcolor } bind def\n");
  fprintf(ps_file,"/Zero_linewidth { 0.0 } def\n");
  fprintf(ps_file,"/Sphcol { 1 setgray } def\n");
  fprintf(ps_file,"/Sphere {");
  fprintf(ps_file,"  newpath 3 copy 0 360 arc gsave Sphcol fill 0\n");
  fprintf(ps_file,"  setgray 0.5 setlinewidth\n");
  fprintf(ps_file,"  3 copy 0.94 mul 260 350 arc stroke 3 copy 0.87");
  fprintf(ps_file," mul 275 335 arc stroke\n");
  fprintf(ps_file,"  3 copy 0.79 mul 295 315 arc stroke 3 copy 0.8");
  fprintf(ps_file," mul 115 135 arc\n");
  fprintf(ps_file,"  3 copy 0.6 mul 135 115 arcn closepath gsave 1");
  fprintf(ps_file," setgray fill grestore stroke\n");
  fprintf(ps_file,"  3 copy 0.7 mul 115 135 arc stroke 3 copy 0.6");
  fprintf(ps_file," mul 124.9 125 arc\n");
  fprintf(ps_file,"  0.8 mul 125 125.1 arc stroke grestore stroke");
  fprintf(ps_file," } bind def\n");

  /* Write out the different colours */
  for (icolour = 0; icolour < NCOLS + 1; icolour++)
    {
      /* Get the colour name and RGB triplet for this colour */
      get_col_rgb(icolour,c_name,rgb);

      /* Write out the colour definition */
      fprintf(ps_file,"/%s    { %5.3f %5.3f %5.3f  setrgbcolor } def\n",
              c_name,rgb[0],rgb[1],rgb[2]);
    }

  /* Loop to write out the different grey-shades */
  for (igrey = 0; igrey < 11; igrey++)
  {
      grey_shade = (float)igrey / (float) 10.0;
      sprintf(colno,"%2d",igrey);
      colno[2] = '\0';

      if (colno[0] == ' ')
          colno[0] = '0';
      fprintf(ps_file,"/Grey%s { %6.2f setgray } def\n",colno,grey_shade);
  }

  /* Print the remainder of the definitions */
  fprintf(ps_file,"/Poly3 { moveto lineto lineto fill grestore } bind ");
  fprintf(ps_file,"def\n");
  fprintf(ps_file,"/Pl3 { 6 copy Poly3 moveto moveto moveto closepath ");
  fprintf(ps_file,"stroke } bind def\n");
  fprintf(ps_file,"/Pline3 { 6 copy Poly3 moveto lineto lineto closepa");
  fprintf(ps_file,"th stroke } bind def\n");
  fprintf(ps_file,"/Poly4 { moveto lineto lineto lineto fill grestore } ");
  fprintf(ps_file,"bind def\n");
  fprintf(ps_file,"/Pl4 { 8 copy Poly4 moveto moveto moveto moveto ");
  fprintf(ps_file,"closepath stroke } bind def\n");
  fprintf(ps_file,"/Pline4 { 8 copy Poly4 moveto lineto lineto lineto");
  fprintf(ps_file," closepath stroke } bind def\n");
  fprintf(ps_file,"/Circol { 1 setgray } def\n");
  fprintf(ps_file,"/Circle { gsave newpath 0 360 arc gsave Circol fill ");
  fprintf(ps_file,"grestore stroke grestore } bind def\n");
  fprintf(ps_file,"/Ocircle { gsave newpath 0 360 arc stroke grestore } ");
  fprintf(ps_file," bind def\n");
  fprintf(ps_file,"/PC { moveto arc gsave Circol fill grestore } bind def\n");
  fprintf(ps_file,"/Pcircle { 7 copy PC moveto arc closepath stroke } bind ");
  fprintf(ps_file,"def\n");
  fprintf(ps_file,"/Ucircle { gsave newpath 0 360 arc gsave Circol fill ");
  fprintf(ps_file,"grestore grestore } bind def\n");
  fprintf(ps_file,"/Arc { newpath arc stroke newpath } bind def\n");
  fprintf(ps_file,"/Ellipse { gsave translate scale Circle grestore } ");
  fprintf(ps_file,"bind def\n");
  fprintf(ps_file,"/Ellarc { gsave translate scale Arc grestore } ");
  fprintf(ps_file,"bind def\n");
  fprintf(ps_file,"/Print { /Times-Roman findfont exch scalefont setfont");
  fprintf(ps_file," show } bind def\n");
  fprintf(ps_file,"/Bprint { /Times-Bold findfont exch scalefont setfont");
  fprintf(ps_file," show } bind def\n");
  fprintf(ps_file,"/Iprint { /Times-Italic findfont exch scalefont setfont");
  fprintf(ps_file," show } bind def\n");
  fprintf(ps_file,"/Aprint { /Arial-Bold findfont exch scalefont setfont");
  fprintf(ps_file," show } bind def\n");
  fprintf(ps_file,"/Gprint { /Symbol findfont exch scalefont setfont show");
  fprintf(ps_file," } bind def\n");
  fprintf(ps_file,"/Lalign {\n");
  fprintf(ps_file,"  dup /Times-Bold findfont exch scalefont setfont\n");
  fprintf(ps_file,"  exch stringwidth pop pop -3 div 0.0 exch rmoveto\n");
  fprintf(ps_file," } bind def\n");
  fprintf(ps_file,"/Ralign {\n");
  fprintf(ps_file,"  dup /Times-Bold findfont exch scalefont setfont\n");
  fprintf(ps_file,"  exch stringwidth pop -1 mul exch -3 div rmoveto\n");
  fprintf(ps_file," } bind def\n");
  fprintf(ps_file,"/Center {\n");
  fprintf(ps_file,"  dup /Times-Roman findfont exch scalefont setfont\n");
  fprintf(ps_file,"  exch stringwidth pop -2 div exch -3 div rmoveto\n");
  fprintf(ps_file," } bind def\n");
  fprintf(ps_file,"/Centerx {\n");
  fprintf(ps_file,"  dup /Times-Roman findfont exch scalefont setfont\n");
  fprintf(ps_file,"  exch stringwidth pop -2 div exch 0 mul rmoveto\n");
  fprintf(ps_file," } bind def\n");
  fprintf(ps_file,"/CenterRot90 {\n");
  fprintf(ps_file,"  dup /Times-Roman findfont exch scalefont setfont\n");
  fprintf(ps_file,"  exch stringwidth pop -2 div exch 3 div exch rmoveto\n");
  fprintf(ps_file," } bind def\n");
  fprintf(ps_file,"/UncenterRot90 {\n");
  fprintf(ps_file,"  dup /Times-Roman findfont exch scalefont setfont\n");
  fprintf(ps_file,"  exch stringwidth } bind def\n");
  fprintf(ps_file,"/Rot90 { gsave currentpoint translate 90 rotate }");
  fprintf(ps_file," bind def\n");
  fprintf(ps_file,"%%%%EndProlog\n");
  fprintf(ps_file,"%%%%BeginSetup\n");
  fprintf(ps_file,"1 setlinecap 1 setlinejoin 1 setlinewidth 0 setgray");
  fprintf(ps_file," [ ] 0 setdash newpath\n");
  fprintf(ps_file,"%%%%EndSetup\n");

  /* Return file pointer */
  return(ps_file);
}
/***********************************************************************

pspage_  -  Write out start of new PostScript page

***********************************************************************/

void pspage_(FILE *ps_file,int page,int in_colour,int landscape,
             char *program_name)
{
  float xx1, xx2, xy1, xy2;

  /* Define border round picture */
  xx1 = (float)BBOXX1;
  xx2 = (float)BBOXX2;
  xy1 = (float)BBOXY1;
  xy2 = (float)BBOXY2;

  /* No border round picture */
  xx1 = -1.0;
  xx2 = 650.0;
  xy1 = -1.0;
  xy2 = 951.0;

  /* Print the page-opening commands */
  fprintf(ps_file,"%%%%Page: %d %d \n",page,page);
  fprintf(ps_file,"/%sSave save def\n",program_name);
  fprintf(ps_file,"gsave\n");

  /* If picture to be plotted in landscape mode, then rotate page */
  if (landscape == TRUE)
    psrot_(ps_file,BBOXX2 + BBOXX1,0.0);
}
/***********************************************************************

psendp_  -  Close off the current PostScript page

***********************************************************************/

void psendp_(FILE *ps_file,float lim[4],char *program_name)
{
  /* Write out the clip command (only works for Portrait mode) */
  fprintf(ps_file,"%%convert -clip %dx%d+%d-%d\n",
          (int) (lim[1] - lim[0] + 10),(int) (lim[3] - lim[2] + 10),
          (int) lim[0],(int) lim[2]);

  fprintf(ps_file,"\n\n");
  fprintf(ps_file,"%sSave restore\n",program_name);
  fprintf(ps_file,"showpage\n");

  /* Re-initialise limits */
  lim[1] = 0.0;
  lim[2] = YMARGIN + YHEIGHT;
}
/***********************************************************************

psclos_   Write final lines to PostScript file and close

***********************************************************************/

void psclos_(FILE *ps_file)
{
  /* Write out closing lines to PostScript file*/
  fprintf(ps_file,"%%%%Trailer\n");
  fprintf(ps_file,"%%%%BoundingBox: %d %d %d %d\n",(int) BBOXX1 - 1,
          (int) BBOXY1 - 1, (int) BBOXX2 + 1, (int) BBOXY2 + 1);
  fprintf(ps_file,"%%%%EOF");
  fclose(ps_file);
}
/***********************************************************************

start_new_page  -  Start a new PostScript page

***********************************************************************/

void start_new_page(struct psparams psparams_ptr,char *first_resnum,
                    char *last_resnum,struct pspage *pspage_ptr,
                    struct parameters params)
{
  char filename[FILENAME_LEN];

    /* If have a previous page, then close off */
  if (pspage_ptr->page > 0)
    psendp_(pspage_ptr->ps_file,pspage_ptr->lim,
            psparams_ptr.program_name);

  /* If writing each page to a separate PostScript file, then
     close previous file and form name of next one */
  if (psparams_ptr.splitps == TRUE)
    {
      /* Close off previous file, if there is one */
      if (pspage_ptr->page > 0)
        psclos_(pspage_ptr->ps_file);

      /* Form name of next file */
      sprintf(filename,"%s%4.4d.ps",pspage_ptr->ps_name,
              pspage_ptr->page + 1);
    }
  else
    strcpy(filename,pspage_ptr->ps_name);

  /* If this is the first page, or printing every page to a new
     file, then open the PostScript file */
  if (pspage_ptr->page == 0 || psparams_ptr.splitps == TRUE)
    {
      /* Open the file */
      pspage_ptr->ps_file = psopen_(filename,psparams_ptr.program_name);
    }

  /* Close off previous tdf file, if there is one */
  if (pspage_ptr->page > 0)
    fclose(pspage_ptr->tdf_file);

  /* Form name of next tdf file */
  sprintf(filename,"topoldia%4.4d.tdf",pspage_ptr->page + 1);

  /* Open the tdf file */
  pspage_ptr->tdf_file = open_output_file(filename,TRUE);

  /* Write out the headings to the file */
  fprintf(pspage_ptr->tdf_file,"# Topology diagram for PDB code %s,"
          " chain %c, domain %3.3d\n",params.pdb_code,params.chain,
          params.domain);

  /* Write out the residue range */
  fprintf(pspage_ptr->tdf_file,"RESIDUERANGE %s 1 %c %s %s\n",
          params.pdb_code,params.chain,first_resnum,last_resnum);

  /* Increment page count */
  pspage_ptr->page++;

  /* Initialise the new PostScript page */
  pspage_(pspage_ptr->ps_file,pspage_ptr->page,psparams_ptr.in_colour,
          psparams_ptr.landscape,psparams_ptr.program_name);

  /* Write out the headings for the plot */
          
  /* Initialise all the variables for the current page */
  pspage_ptr->y = pspage_ptr->Page_Max_y;
}
/***********************************************************************

close_page  -  Close current PostScript file if have an open one

***********************************************************************/

void close_page(struct psparams psparams_ptr,struct pspage *pspage_ptr)
{
  /* Close off previous file, if there is one */
  if (pspage_ptr->page > 0)
    {
      /* Finish off the PostScript file and close */
      psendp_(pspage_ptr->ps_file,pspage_ptr->lim,
            psparams_ptr.program_name);
      psclos_(pspage_ptr->ps_file);

      /* Close the tdf file */
      fclose(pspage_ptr->tdf_file);
    }
}
/***********************************************************************

plot_strand  -  Plot the given gene arrow

***********************************************************************/

void plot_strand(float xstart,float ystart,float xend,float yend,
                 int orientation,int direction,float scale_factor,
                 int shadow,char *sst_id,char *first_resnum,
                 char *last_resnum,int sse_len,int strand_no,
                 int domain_col,struct pspage *pspage_ptr,
                 struct parameters params)
{
  char arrow_col[COLNAM_LEN], c_name[COLNAM_LEN];
  char name[LINELEN];

  static char *coil_col = "Blue";
  static char *line_col = "Black";
  static char *shade = "Grey";
  static char *text_col = "Black";

  double arrow_cut, arrow_height, arrow_width, arrow_size, line_width;
  double coil_width, overhang, shadow_shift, text_size;
  double diff, len1, len2;
  double x1, x2, x3, x4, y1, y2, y3, y4, xmid, ymid;
  double vx1, vx2, vx3, vx4, vx5, vx6, vx7;
  double vy1, vy2, vy3, vy4, vy5, vy6, vy7;
  double xtext1, xtext2, ytext1, ytext2;

  /* Initialise variables */
  vx1 = vx2 = vx3 = vx4 = vx5 = vx6 = vx7 = 0;
  vy1 = vy2 = vy3 = vy4 = vy5 = vy6 = vy7 = 0;

  /* Get the colour of the arrow */
  strcpy(arrow_col,SSE_COLOUR[domain_col][1]);

  /* Apply scale factor to size parameters */
  coil_width = COIL_WIDTH * scale_factor;
  line_width = LINE_WIDTH * scale_factor;
  arrow_size = ARROW_SIZE * scale_factor;
  arrow_height = 1.5 * arrow_size;
  arrow_width = arrow_height;
  arrow_cut = ARROW_CUT * scale_factor;
  overhang = COIL_OVERHANG * scale_factor;
  shadow_shift = SHADOW_WIDTH * scale_factor;
  text_size = RES_NUMBERS * scale_factor;
  if (text_size < TEXT_SIZE_MIN)
    text_size = TEXT_SIZE_MIN;

  /* Get lengths of start and end residue numbers */
  len1 = strlen(first_resnum);
  len2 = strlen(last_resnum);
  if (sse_len < 3 || len1 > 3 || len2 > 3)
    text_size = 3 * text_size / 4; 

  /* Form strand name for PostScript comment */
  if (params.show_sse_labels == TRUE)
    sprintf(name,"Strand %s: %s to %s",sst_id,first_resnum,last_resnum);
  else
    sprintf(name,"Strand %s to %s",first_resnum,last_resnum);
  if (shadow == TRUE)
    strcat(name," - shadow");

  /* Start this strand's entry in the tdf file */
  if (shadow == FALSE)
    fprintf(pspage_ptr->tdf_file,"STRAND %s %s ",first_resnum,last_resnum);

  /* Print comment in PostScript file */
  pscomm_(pspage_ptr->ps_file,name);

  /* Calculate coordinates of bar representing the bar of the arrow */
  if (orientation == HORIZONTAL)
    {
      x1 = xstart;
      x2 = xend;
      xmid = (x1 + x2) / 2;
      y1 = ystart - arrow_size;
      y2 = yend + arrow_size;
      y3 = y4 = ytext1 = ytext2 = ystart;
      ymid = (y1 + y2) / 2;

      /* If forward arrow, then lop right-hand end */
      if (direction == FORWARD)
        {
          x2 = x2 - arrow_width;
          x3 = xstart - overhang;
          x4 = xend + overhang;
          xtext1 = x1 + (1 + len1 / 2.0) * text_size * CHAR_ASPECT;
          xtext2 = x2 + arrow_width / 2.0
            - (len2 / 2.0) * text_size * CHAR_ASPECT;

          /* Store the four vertices of the arrow shaft */
          vx1 = x1;
          vy1 = y2;
          vx2 = x2;
          vy2 = y2;
          vx6 = x2;
          vy6 = y1;
          vx7 = x1;
          vy7 = y1;
        }

      /* Otherwise, lop off left-hand end */
      else
        {
          x2 = x2 + arrow_width;
          x3 = xend - overhang;
          x4 = xstart + overhang;
          xtext1 = x1 - (1 + len1 / 2.0) * text_size * CHAR_ASPECT;
          xtext2 = x2 - arrow_width / 2.0
            + (len2 / 2.0) * text_size * CHAR_ASPECT;

          /* Store the four vertices of the arrow shaft */
          vx1 = x1;
          vy1 = y1;
          vx2 = x2;
          vy2 = y1;
          vx6 = x2;
          vy6 = y2;
          vx7 = x1;
          vy7 = y2;
        }
    }
  else
    {
      x1 = xstart - arrow_size;
      x2 = xend + arrow_size;
      xmid = (x1 + x2) / 2;
      x3 = x4 = xtext1 = xtext2 = xstart;
      y1 = ystart;
      y2 = yend;
      ymid = (y1 + y2) / 2;

      /* If forward arrow, then lop top end */
      if (direction == FORWARD)
        {
          y2 = y2 - arrow_width;
          y3 = ystart - overhang;
          y4 = yend + overhang;
          ytext1 = y1 + 2 * text_size / 3;
          ytext2 = y2 + arrow_width / 4.0;

          /* Adjust for overlaps */
          ymid = (ytext1 + ytext2) / 2;
          diff = (ytext2 - ytext1) - text_size;
          if (params.show_sse_labels == TRUE)
            diff = diff - text_size;
          if (diff < 0)
            {
              ytext1 = ytext1 + diff / 2;
              ytext2 = ytext2 - diff / 2;
            }

          /* Store the four vertices of the arrow shaft */
          vx1 = x1;
          vy1 = y1;
          vx2 = x1;
          vy2 = y2;
          vx6 = x2;
          vy6 = y2;
          vx7 = x2;
          vy7 = y1;
        }

      /* Otherwise, lop off bottom end */
      else
        {
          y2 = y2 + arrow_width;
          y3 = yend - overhang;
          y4 = ystart + overhang;
          ytext1 = y1 - 2 * text_size / 3;
          ytext2 = y2 - arrow_width / 4.0;

          /* Adjust for overlaps */
          ymid = (ytext1 + ytext2) / 2;
          diff = (ytext1 - ytext2) - text_size;
          if (params.show_sse_labels == TRUE)
            diff = diff - text_size;
          if (diff < 0)
            {
              ytext1 = ytext1 - diff / 2;
              ytext2 = ytext2 + diff / 2;
            }

          /* Store the four vertices of the arrow shaft */
          vx1 = x2;
          vy1 = y1;
          vx2 = x2;
          vy2 = y2;
          vx6 = x1;
          vy6 = y2;
          vx7 = x1;
          vy7 = y1;
        }
    }

  /* Determine arrow colour */
  if (shadow == TRUE)
    {
      /* Set colour to grey */
      strcpy(c_name,shade);

      /* Shift all the coords */
      x1 = x1 + shadow_shift;
      x2 = x2 + shadow_shift;
      y1 = y1 - shadow_shift;
      y2 = y2 - shadow_shift;
    }
  else
    {
      /* Set strand colour */
      strcpy(c_name,arrow_col);

      /* Plot coil, which extends either end of the strand */
      pscolb_(pspage_ptr->ps_file,coil_col);
      pslwid_(pspage_ptr->ps_file,coil_width);
      psline_(pspage_ptr->ps_file,x3,y3,x4,y4);
    }

  /* Draw the arrow shaft */
  pslwid_(pspage_ptr->ps_file,0.0);
  pshade_(pspage_ptr->ps_file,c_name);
  psubox_(pspage_ptr->ps_file,x1,y1,x1,y2,x2,y2,x2,y1);

  /* If not plotting the shadow, draw boundary lines */
  if (shadow == FALSE)
    {
      pscolb_(pspage_ptr->ps_file,line_col);
      pslwid_(pspage_ptr->ps_file,line_width);

      /* Plot lines according to orientation */
      if (orientation == HORIZONTAL)
        {
          psline_(pspage_ptr->ps_file,x1,y1,x2,y1);
          psline_(pspage_ptr->ps_file,x1,y2,x2,y2);
      
          /* Close off arrow's end */
          psline_(pspage_ptr->ps_file,x1,y1,x1,y2);
        }
      else
        {
          psline_(pspage_ptr->ps_file,x1,y1,x1,y2);
          psline_(pspage_ptr->ps_file,x2,y1,x2,y2);
      
          /* Close off arrow's end */
          psline_(pspage_ptr->ps_file,x1,y1,x2,y1);
        }
    }

  /* Place arrow-head according to orientation */
  if (orientation == HORIZONTAL)
    {
      /* If arrow on the right */
      if (direction == FORWARD)
        {
          x1 = x2;
          x3 = x1 + arrow_width;
        }
      else
        {
          x1 = x2;
          x3 = x1 - arrow_width;
        }
      y3 = y1 + arrow_size;
      y1 = y3 - arrow_height;
      y2 = y3 + arrow_height;

      /* Save the three vertices */
      vx3 = vx5 = x1;
      if (direction == FORWARD)
        {
          vy3 = y2;
          vy5 = y1;
        }
      else
        {
          vy3 = y1;
          vy5 = y2;
        }
      vx4 = x3;
      vy4 = y3;
    }
  else
    {
      /* If arrow on the top */
      if (direction == FORWARD)
        {
          y1 = y2;
          y3 = y1 + arrow_width;
        }
      else
        {
          y1 = y2;
          y3 = y1 - arrow_width;
        }
      x3 = x1 + arrow_size;
      x1 = x3 - arrow_height;
      x2 = x3 + arrow_height;

      /* Save the three vertices */
      vy3 = vy5 = y1;
      if (direction == FORWARD)
        {
          vx3 = x1;
          vx5 = x2;
        }
      else
        {
          vx3 = x2;
          vx5 = x1;
        }
      vx4 = x3;
      vy4 = y3;
    }

  /* Define colour */
  pslwid_(pspage_ptr->ps_file,0.0);
  pshade_(pspage_ptr->ps_file,c_name);

  /* Plot the arrow-head */
  psutri_(pspage_ptr->ps_file,x1,y1,x2,y2,x3,y3);

  /* If not showing the shadow, draw border lines and show residue
     numbers */
  if (shadow == FALSE)
    {
      pscolb_(pspage_ptr->ps_file,line_col);
      pslwid_(pspage_ptr->ps_file,line_width);
      if (orientation == HORIZONTAL)
        {
          psline_(pspage_ptr->ps_file,x1,y1,x1,y1 + arrow_size / 2);
          psline_(pspage_ptr->ps_file,x1,y1,x3,y3);
          psline_(pspage_ptr->ps_file,x3,y3,x2,y2);
          psline_(pspage_ptr->ps_file,x2,y2,x2,y2 - arrow_size / 2);
        }
      else
        {
          psline_(pspage_ptr->ps_file,x1,y1,x1 + arrow_size / 2,y1);
          psline_(pspage_ptr->ps_file,x1,y1,x3,y3);
          psline_(pspage_ptr->ps_file,x3,y3,x2,y2);
          psline_(pspage_ptr->ps_file,x2,y2,x2 - arrow_size / 2,y2);
        }

      /* Plot strand id, if required */
      if (params.show_sse_labels == TRUE)
        {
          /* Get coords */
          pscolb_(pspage_ptr->ps_file,text_col);
          psbtxt_(pspage_ptr->ps_file,xmid,ymid,text_size,sst_id);
        }

      /* Plot the residue numbers */
      if (params.show_residue_nos == TRUE)
        {

          /* Write out the residue numbers */
          pscolb_(pspage_ptr->ps_file,text_col);
          psitxt_(pspage_ptr->ps_file,xtext1,ytext1,text_size,
                  first_resnum);
          psitxt_(pspage_ptr->ps_file,xtext2,ytext2,text_size,
                  last_resnum);
        }
    }

  /* Write out the vertices to the tdf file */
  if (orientation == HORIZONTAL)
    {
      psline_(pspage_ptr->ps_file,x1,y1,x1,y1 + arrow_size / 2);
      psline_(pspage_ptr->ps_file,x1,y1,x3,y3);
      psline_(pspage_ptr->ps_file,x3,y3,x2,y2);
      psline_(pspage_ptr->ps_file,x2,y2,x2,y2 - arrow_size / 2);
    }
  else
    {
      psline_(pspage_ptr->ps_file,x1,y1,x1 + arrow_size / 2,y1);
      psline_(pspage_ptr->ps_file,x1,y1,x3,y3);
      psline_(pspage_ptr->ps_file,x3,y3,x2,y2);
      psline_(pspage_ptr->ps_file,x2,y2,x2 - arrow_size / 2,y2);
    }

  /* Close this strand's entry in the tdf file */
  if (shadow == FALSE)
    fprintf(pspage_ptr->tdf_file,
            " %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f"
            " %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",vx1, vy1,
            vx2, vy2, vx3, vy3, vx4, vy4, vx5, vy5, vx6, vy6, vx7, vy7);
}
/***********************************************************************

plot_helix  -  Plot a cylinder to represent the given helix

***********************************************************************/

void plot_helix(float xstart,float ystart,float xend,float yend,
                int orientation,int direction,float scale_factor,
                int shadow,char *sst_id,char *first_resnum,
                char *last_resnum,int sse_len,int helix_no,
                int domain_col,struct pspage *pspage_ptr,
                struct parameters params)
{
  char helix_col[COLNAM_LEN], c_name[COLNAM_LEN];
  char name[LINELEN];

  static char *coil_col = "Blue";
  static char *line_col = "Black";
  static char *shade = "Grey";
  static char *text_col = "Black";

  double hellipse, helix_width, line_width;
  double coil_width, overhang, shadow_shift, text_size;
  double diff, len1, len2, xmid, ymid;
  double x1, x2, x3, x4, x5, x6, y1, y2, y3, y4, y5, y6;
  double vx1, vx2, vy1, vy2, emin, emax;
  double angle_start1, angle_start2, angle_end1, angle_end2;
  double elong_x, elong_y;
  double xtext1, xtext2, ytext1, ytext2;

  /* Initial;ise variables */
  vx1 = vx2 = vy1 = vy2 = emin = emax = 0;

  /* Get the colour of the helix */
  strcpy(helix_col,SSE_COLOUR[domain_col][0]);

  /* Apply scale factor to size parameters */
  coil_width = COIL_WIDTH * scale_factor;
  hellipse = HELIX_ELLIPSE * scale_factor;
  helix_width = HELIX_WIDTH * scale_factor;
  line_width = LINE_WIDTH * scale_factor;
  overhang = COIL_OVERHANG * scale_factor;
  shadow_shift = SHADOW_WIDTH * scale_factor;
  text_size = RES_NUMBERS * scale_factor;
  if (text_size < TEXT_SIZE_MIN)
    text_size = TEXT_SIZE_MIN;

  /* Get lengths of start and end residue numbers */
  len1 = strlen(first_resnum);
  len2 = strlen(last_resnum);
  if (sse_len < 5 || len1 > 3 || len2 > 3)
    text_size = 3 * text_size / 4; 

  /* Form helix name for PostScript comment */
  if (params.show_sse_labels == TRUE)
    sprintf(name,"Helix %s: %s to %s",sst_id,first_resnum,last_resnum);
  else
    sprintf(name,"Helix %s to %s",first_resnum,last_resnum);
  if (shadow == TRUE)
    strcat(name," - shadow");

  /* Start this helix's entry in the tdf file */
  if (shadow == FALSE)
    fprintf(pspage_ptr->tdf_file,"HELIX %s %s ",first_resnum,last_resnum);

  /* Print comment in PostScript file */
  pscomm_(pspage_ptr->ps_file,name);

  /* Calculate coordinates of bar representing the bar of the helix */
  if (orientation == HORIZONTAL)
    {
      elong_x = 0.5;
      elong_y = 1.0;
      y1 = ystart - helix_width;
      y2 = yend + helix_width;
      y3 = y4 = y5 = y6 = ytext1 = ytext2 = ystart;
      ymid = (y1 + y2) / 2;

      /* Bring helix ends in */
      if (direction == FORWARD)
        {
          x1 = xstart + hellipse;
          x2 = xend - hellipse;
          xmid = (x1 + x2) / 2;
          x3 = xstart - overhang;
          x4 = xend + overhang;
          xtext1 = xstart + (1 + len1 / 2.0) * text_size * CHAR_ASPECT;
          xtext2 = x2 - hellipse
            - (len2 / 2.0) * text_size * CHAR_ASPECT;
        }

      /* Otherwise, lop off left-hand end */
      else
        {
          x1 = xstart - hellipse;
          x2 = xend + hellipse;
          xmid = (x1 + x2) / 2;
          x3 = xend - overhang;
          x4 = xstart + overhang;
          xtext1 = xstart - (1 + len1 / 2.0) * text_size * CHAR_ASPECT;
          xtext2 = x2 + hellipse
            + (len2 / 2.0) * text_size * CHAR_ASPECT;
        }

      /* Get coords of centre of ellipse */
      x5 = x1;
      x6 = x2;
    }
  else
    {
      elong_x = 1.0;
      elong_y = 0.5;
      x1 = xstart - helix_width;
      x2 = xend + helix_width;
      xmid = (x1 + x2) / 2;
      x3 = x4 = x5 = x6 = xtext1 = xtext2 = xstart;

      /* If forward helix, then lop top end */
      if (direction == FORWARD)
        {
          y1 = y5 = ystart + hellipse;
          y2 = y6 = yend - hellipse;
          y3 = ystart - overhang;
          y4 = yend + overhang;
          ytext1 = ystart + text_size;
          ytext2 = y2 - hellipse - text_size / 2;

          /* Adjust for overlaps */
          ymid = (ytext1 + ytext2) / 2;
          diff = (ytext2 - ytext1) - text_size;
          if (params.show_sse_labels == TRUE)
            diff = diff - text_size;
          if (diff < 0)
            {
              ytext1 = ytext1 + diff / 2;
              ytext2 = ytext2 - diff / 2;
            }
        }

      /* Otherwise, lop off bottom end */
      else
        {
          y1 = y6 = ystart - hellipse;
          y2 = y5 = yend + hellipse;
          y3 = yend - overhang;
          y4 = ystart + overhang;
          ytext1 = y1 - hellipse - text_size / 2;
          ytext2 = yend + text_size;

          /* Adjust for overlaps */
          ymid = (ytext1 + ytext2) / 2;
          diff = (ytext1 - ytext2) - text_size;
          if (params.show_sse_labels == TRUE)
            diff = diff - text_size;
          if (diff < 0)
            {
              ytext1 = ytext1 - diff / 2;
              ytext2 = ytext2 + diff / 2;
            }
        }

      /* Get coords of centre of ellipse */
      angle_start1 = 180.0;
      angle_start2 = 0.0;
      angle_end1 = angle_end2 = 360.0;
    }

  /* Determine helix colour */
  if (shadow == TRUE)
    {
      /* Set colour to grey */
      strcpy(c_name,shade);

      /* Shift all the coords */
      x1 = x1 + shadow_shift;
      x2 = x2 + shadow_shift;
      x5 = x5 + shadow_shift;
      x6 = x6 + shadow_shift;
      y1 = y1 - shadow_shift;
      y2 = y2 - shadow_shift;
      y5 = y5 - shadow_shift;
      y6 = y6 - shadow_shift;
    }
  else
    {
      /* Set helix colour */
      strcpy(c_name,helix_col);

      /* Plot coil, which extends either end of the helix */
      pscolb_(pspage_ptr->ps_file,coil_col);
      pslwid_(pspage_ptr->ps_file,coil_width);
      psline_(pspage_ptr->ps_file,x3,y3,x4,y4);
    }

  /* Draw the helix shaft */
  pslwid_(pspage_ptr->ps_file,0.0);
  pshade_(pspage_ptr->ps_file,c_name);
  psubox_(pspage_ptr->ps_file,x1,y1,x1,y2,x2,y2,x2,y1);

  /* If not plotting the shadow, draw boundary lines */
  if (shadow == FALSE)
    {
      pscolb_(pspage_ptr->ps_file,line_col);
      pslwid_(pspage_ptr->ps_file,line_width);

      /* Plot lines according to orientation */
      if (orientation == HORIZONTAL)
        {
          psline_(pspage_ptr->ps_file,x1,y1,x2,y1);
          psline_(pspage_ptr->ps_file,x1,y2,x2,y2);
          emax = fabs (y2 - y1) / 2;
        }
      else
        {
          psline_(pspage_ptr->ps_file,x1,y1,x1,y2);
          psline_(pspage_ptr->ps_file,x2,y1,x2,y2);
          emax = fabs (x2 - x1) / 2;
        }
    }

  /* Save the tdf vertices */
  vx1 = x1;
  vy1 = y1;
  vx2 = x2;
  vy2 = y2;
  emin = emax / 2;

  /* Define colour */
  pslwid_(pspage_ptr->ps_file,0.0);
  pscolb_(pspage_ptr->ps_file,c_name);
  psccol_(pspage_ptr->ps_file,c_name);

  /* Plot the ellipse at the helix head */
  pselip_(pspage_ptr->ps_file,x5,y5,helix_width,elong_x,elong_y);
  pselip_(pspage_ptr->ps_file,x6,y6,helix_width,elong_x,elong_y);

  /* If not showing the shadow, draw elliptical borders */
  if (shadow == FALSE)
    {
      pscolb_(pspage_ptr->ps_file,line_col);
      pslwid_(pspage_ptr->ps_file,line_width);
      psearc_(pspage_ptr->ps_file,x5,y5,helix_width,angle_start1,
              angle_end1,elong_x,elong_y);
      psearc_(pspage_ptr->ps_file,x6,y6,helix_width,angle_start2,
              angle_end2,elong_x,elong_y);

      /* Plot coil extending from top ellipse */
      pscolb_(pspage_ptr->ps_file,coil_col);
      pslwid_(pspage_ptr->ps_file,coil_width);
      psline_(pspage_ptr->ps_file,x6,y6,x4,y4);

      /* Plot helix id, if required */
      if (params.show_sse_labels == TRUE)
        {
          /* Get coords */
          pscolb_(pspage_ptr->ps_file,text_col);
          psbtxt_(pspage_ptr->ps_file,xmid,ymid,text_size,sst_id);
        }

      /* Plot the residue numbers */
      if (params.show_residue_nos == TRUE)
        {
          pscolb_(pspage_ptr->ps_file,text_col);
          psitxt_(pspage_ptr->ps_file,xtext1,ytext1,text_size,
                  first_resnum);
          psitxt_(pspage_ptr->ps_file,xtext2,ytext2,text_size,
                  last_resnum);
        }
    }

  /* Close this helix's entry in the tdf file */
  if (shadow == FALSE)
    fprintf(pspage_ptr->tdf_file," %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",
            vx1,vy1,vx2,vy2,emax,emin);
}
/***********************************************************************

plot_terminus  -  Plot a square to represent the given terminus

***********************************************************************/

void plot_terminus(float xstart,float ystart,float xend,float yend,
                   int direction,float scale_factor,int shadow,int type,
                   char *res_num,struct pspage *pspage_ptr)
{
  char c_name[COLNAM_LEN];
  char label[2], name[LINELEN];

  static char *coil_col = "Blue";
  static char *line_col = "Black";
  static char *shade = "Grey";
  static char *terminus_col = "LightGoldenrod1";
  static char *text_col = "Black";

  double box_size, coil_width, line_width, overhang;
  double shadow_shift, text_size;
  double x1, x2, x3, x4, y1, y2, y3, y4;
  double xtext, ytext;

  /* Apply scale factor to size parameters */
  box_size = TERMINUS_WIDTH * scale_factor;
  coil_width = COIL_WIDTH * scale_factor;
  overhang = COIL_OVERHANG * scale_factor;
  line_width = LINE_WIDTH * scale_factor;
  shadow_shift = SHADOW_WIDTH * scale_factor;
  text_size = TERMINUS_TEXT * scale_factor;
  if (text_size < TEXT_SIZE_MIN)
    text_size = TEXT_SIZE_MIN;

  /* Form terminus name for PostScript comment */
  if (type == NTERMINUS)
    {
      strcpy(name,"N-terminus");
      strcpy(label,"N");
    }
  else
    {
      strcpy(name,"C-terminus");
      strcpy(label,"C");
    }
  if (shadow == TRUE)
    strcat(name," - shadow");

  /* Start this terminus' entry in the tdf file */
  if (shadow == FALSE)
    fprintf(pspage_ptr->tdf_file,"TERM %s %s ",label,res_num);

  /* Print comment in PostScript file */
  pscomm_(pspage_ptr->ps_file,name);

  /* Calculate coordinates of box */
  x1 = xstart - box_size;
  x2 = xstart + box_size;
  x3 = x4 = xstart;
  y1 = ystart;
  y2 = yend;

  /* Calculate coords of label */
  xtext = xstart;
  y3 = ytext = (ystart + yend) / 2;

  /* Calculate coords of end of coil overhang */
  if (type == NTERMINUS)
    if (direction == FORWARD)
      y4 = yend + overhang;
    else
      y4 = yend - overhang;
  else
    if (direction == FORWARD)
      y4 = ystart - overhang;
    else
      y4 = ystart + overhang;

  /* Determine terminus colour */
  if (shadow == TRUE)
    {
      /* Set colour to grey */
      strcpy(c_name,shade);

      /* Shift all the coords */
      x1 = x1 + shadow_shift;
      x2 = x2 + shadow_shift;
      y1 = y1 - shadow_shift;
      y2 = y2 - shadow_shift;
    }
  else
    {
      /* Set terminus colour */
      strcpy(c_name,terminus_col);

      /* Plot coil, extending from one end of the terminus box */
      pscolb_(pspage_ptr->ps_file,coil_col);
      pslwid_(pspage_ptr->ps_file,coil_width);
      psline_(pspage_ptr->ps_file,x3,y3,x4,y4);
    }

  /* Draw the terminus box */
  pslwid_(pspage_ptr->ps_file,0.0);
  pshade_(pspage_ptr->ps_file,c_name);
  psubox_(pspage_ptr->ps_file,x1,y1,x1,y2,x2,y2,x2,y1);

  /* If not plotting the shadow, draw boundary lines */
  if (shadow == FALSE)
    {
      pscolb_(pspage_ptr->ps_file,line_col);
      pslwid_(pspage_ptr->ps_file,line_width);

      /* Plot lines */
      psline_(pspage_ptr->ps_file,x1,y1,x2,y1);
      psline_(pspage_ptr->ps_file,x1,y2,x2,y2);
      psline_(pspage_ptr->ps_file,x1,y1,x1,y2);
      psline_(pspage_ptr->ps_file,x2,y1,x2,y2);
    }

  /* If not showing the shadow, label to box */
  if (shadow == FALSE)
    {
      /* Plot the appropriate label */
      pscolb_(pspage_ptr->ps_file,text_col);
      psctxt_(pspage_ptr->ps_file,xtext,ytext,text_size,label);
    }

  /* Close this terminus' entry in the tdf file */
  if (shadow == FALSE)
    fprintf(pspage_ptr->tdf_file," %8.3f %8.3f %8.3f %8.3f\n",
            x1,y1,x2,y2);
}
/***********************************************************************

plot_coil  -  Plot the given gene arrow

***********************************************************************/

void plot_coil(float xstart,float ystart,float xend,float yend,
               float *xtdf,float *ytdf,int orientation,int direction,
               float scale_factor,int shadow,int curve,
               struct pspage *pspage_ptr)
{
  char name[LINELEN];

  static char *coil_col = "Blue";

  double arrow_height, arrow_width, coil_width;
  double x1, x2, x3, x4, x5, y1, y2, y3, y4, y5;

  /* Skip shadow plot */
  if (shadow == TRUE)
    return;

  /* Apply scale factor to size parameters */
  arrow_height = COIL_ARROW_HEIGHT * scale_factor;
  arrow_width = COIL_ARROW_WIDTH * scale_factor;
  coil_width = COIL_WIDTH * scale_factor;

  /* Form PostScript comment */
  sprintf(name,"Coil");

  /* Print comment in PostScript file */
  pscomm_(pspage_ptr->ps_file,name);

  /* Calculate coil's end coordinates */
  x1 = xstart;
  x2 = xend;
  y1 = ystart;
  y2 = yend;

  /* If this is a reasonably short end segment, plot as semi circle */
  //if (curve == TRUE)
  //  {
  //  }

  /* Otherwise, plot coil as a straight-line segment */
  //else
  {
    pscolb_(pspage_ptr->ps_file,coil_col);
    pslwid_(pspage_ptr->ps_file,coil_width);
    psline_(pspage_ptr->ps_file,x1,y1,x2,y2);
  }

  /* Write end-vertex to tdf file */
  fprintf(pspage_ptr->tdf_file," %8.3f %8.3f ",xtdf[1],ytdf[1]);

  /* Calculate coordinates of direction arrow */
  if (orientation == HORIZONTAL)
    {
      y3 = ystart;
      y4 = ystart - arrow_width;
      y5 = ystart + arrow_width;
      x3 = (xstart + xend) / 2;
      if (direction == FORWARD)
        {
          x4 = x5 = x3 - arrow_height / 2;
          x3 = x3 + arrow_height / 2;
        }
      else
        {
          x4 = x5 = x3 + arrow_height / 2;
          x3 = x3 - arrow_height / 2;
        }
    }
  else
    {
      x3 = xstart;
      x4 = xstart - arrow_width;
      x5 = xstart + arrow_width;
      y3 = (ystart + yend) / 2;
      if (direction == FORWARD)
        {
          y4 = y5 = y3 - arrow_height / 2;
          y3 = y3 + arrow_height / 2;
        }
      else
        {
          y4 = y5 = y3 + arrow_height / 2;
          y3 = y3 - arrow_height / 2;
        }
    }

  /* Plot the small arrow giving the direction of the coil segment */
  psline_(pspage_ptr->ps_file,x3,y3,x4,y4);
  psline_(pspage_ptr->ps_file,x3,y3,x5,y5);
}
/* TDF output format

## tdf stands for topology diagram format. A tdf file contains data for only 1
topology diagram.
## Anything in a line after a hash character is a comment.
## Blank lines are ignored.
## Multiple consecutive spaces are collapsed to single space, all starting and
trailing spaces are ignored.
## 
## There are following records, a record must sit on a single line:
## 
##     RESIDUERANGE pdbId modelNumber chainId resnumStart resnumEnd
## Note: all residue ranges in this file are inclusive of mentioned residue numbers.
## 
##     HELIX resnumStart resnumEnd startLeftX startLeftY endRightX endRightY
ellipseMajorAxis, ellipseMinorAxis
##
##     STRAND resnumStart resnumEnd v1x v1y v2x v2y v3x v3y v4x v4y v5x v5y v6x v6y
v7x v7y
## Note: coordinates should be ordered such that 4th point is tip of arrow & 1st,7th
points are start of the arrow.
##
##     LOOP resnumStart resnumEnd v1x v1y v2x v2y {...and so on for as many vertices
there are...}
##
##     TERM {N or C} resnum cornerX cornerY otherCornerX otherCornerY
##
## The records need not be arranged in any order.

RESIDUERANGE 1aac 1 A
TERM N 3 5 5 10 10
TERM C 53 95 5 100 10
LOOP 4 7 10 15 10 25
STRAND 8 15 25 25 25 35 20 35 27 33 ........
LOOP 16 17 .......
STRAND 18 25 .......
LOOP 26 28 ......
HELIX 29 36 .......
LOOP 37 42 ......
HELIX 43 49 .......
LOOP 50 52 ......

 */
/***********************************************************************

plot_sses  -  Plot all the secondary structure elements in the plot

***********************************************************************/

void plot_sses(struct sse *first_sse_ptr,float scale_factor,
               float *origin,struct pspage *pspage_ptr,
               struct limits limit,struct parameters params)
{
  int domain_col, loop, height, helix_no, sse_len, strand_no, width;
  int curve;

  float xstart, xend, xtdf[2], ystart, yend, ytdf[2];

  struct sse *sse_ptr;

  static int shadow[2] = { TRUE, FALSE };

  /* Initialise variables */
  helix_no = 0;
  strand_no = 0;

  /* Determine what the domain colour to use */
  if (params.domain == 0)
    domain_col = 0;
  else
    domain_col = (params.domain - 1) % MAX_DOMCOL + 1;

  /* Initialise pointer to first SSE record */
  sse_ptr = first_sse_ptr;

  /* Loop through all the SSE records to update coordinate limits */
  while (sse_ptr != NULL)
    {
      /* Get the coordinates of this SSE */
      xstart = origin[0]
        + scale_factor * (sse_ptr->xpos[0] - limit.min_x);
      xend = origin[0]
        + scale_factor * (sse_ptr->xpos[1] - limit.min_x);
      ystart = origin[1]
        + scale_factor * (sse_ptr->ypos[0] - limit.min_y);
      yend = origin[1]
        + scale_factor * (sse_ptr->ypos[1] - limit.min_y);

      /* Calculate the tdf coords */
      xtdf[0] = origin[0]
        + scale_factor * (sse_ptr->xorig[0] - limit.min_x);
      xtdf[1] = origin[0]
        + scale_factor * (sse_ptr->xorig[1] - limit.min_x);
      ytdf[0] = origin[1]
        + scale_factor * (sse_ptr->yorig[0] - limit.min_y);
      ytdf[1] = origin[1]
        + scale_factor * (sse_ptr->yorig[1] - limit.min_y);

      /* Get the length of this secondary structure unit */
      sse_len = sse_ptr->last_pos - sse_ptr->first_pos + 1;

      /* Increment SSE counter */
      if (sse_ptr->type == STRAND)
        strand_no++;
      else if (sse_ptr->type == HELIX)
        helix_no++;

      /* Loop twice, first to plot shadow and then SSE */
      for (loop = 0; loop < 2; loop++)
        {
          /* If this is a strand, then plot arrow */
          if (sse_ptr->type == STRAND)
            plot_strand(xstart,ystart,xend,yend,sse_ptr->orientation,
                        sse_ptr->direction,scale_factor,shadow[loop],
                        sse_ptr->sst_id,sse_ptr->first_resnum,
                        sse_ptr->last_resnum,sse_len,strand_no,
                        domain_col,pspage_ptr,params);

          /* For helix plot as cylinder */
          else if (sse_ptr->type == HELIX)
            plot_helix(xstart,ystart,xend,yend,sse_ptr->orientation,
                       sse_ptr->direction,scale_factor,
                       shadow[loop],sse_ptr->sst_id,sse_ptr->first_resnum,
                       sse_ptr->last_resnum,sse_len,helix_no,domain_col,
                       pspage_ptr,params);

          /* For N- or C-terminus plot as square */
          else if (sse_ptr->type == NTERMINUS ||
                   sse_ptr->type == CTERMINUS)
            plot_terminus(xstart,ystart,xend,yend,sse_ptr->direction,
                          scale_factor,shadow[loop],sse_ptr->type,
                          sse_ptr->first_resnum,pspage_ptr);

          /* For starting coil overhang, plot coil start */
          else if (sse_ptr->type == OVERHANG1 && shadow[loop] == FALSE)
            fprintf(pspage_ptr->tdf_file,
                    "COIL %s %s  %8.3f %8.3f %8.3f %8.3f",
                    sse_ptr->first_resnum,sse_ptr->last_resnum,
                    xtdf[0],ytdf[0],xtdf[1],ytdf[1]);

          /* For coil plot as joining */
          else if (sse_ptr->type == COIL)
            {
              /* If this is an end segment, determine whether it is
                 sufficiently short to plot as a semi-circle */
              if (sse_ptr->endseg == TRUE)
                {
                  /* Determine height and width */
                  height = abs(sse_ptr->ypos[1] - sse_ptr->ypos[0]);
                  width = abs(sse_ptr->xpos[1] - sse_ptr->xpos[0]);

                  /* If short enough, set curve flag */
                  if (height < MIN_CURVE_DIST && width < MIN_CURVE_DIST)
                    curve = TRUE;
                  else
                    curve = FALSE;
                }
              else
                curve = FALSE;

              /* Plot the coil segment */
              plot_coil(xstart,ystart,xend,yend,xtdf,ytdf,
                        sse_ptr->orientation,sse_ptr->direction,
                        scale_factor,shadow[loop],curve,pspage_ptr);
            }

          /* For closing coil overhang, plot coil end */
          else if (sse_ptr->type == OVERHANG2 && shadow[loop] == FALSE)
            {
              /* Close this coil's entry in the tdf file */
              fprintf(pspage_ptr->tdf_file," %8.3f %8.3f\n",xtdf[1],
                      ytdf[1]);
            }
        }

      /* Get pointer to next sse in linked-list */
      sse_ptr = sse_ptr->next_sse_ptr;
    }
}
/***********************************************************************

write_postscript  -  Write out a PostScript file showing the topology
                     disgram

***********************************************************************/

void write_postscript(struct sse *first_sse_ptr,
                      struct residue *residue_array_ptr,int nres,
                      struct psparams psparams_ptr,
                      struct pspage *pspage_ptr,struct limits limit,
                      struct parameters params)
{
  char first_resnum[6], last_resnum[6];

  float default_height, height, scale_factor, width;
  float frame_aspect, origin[2], picture_aspect;

  /* Initialise variables */
  first_resnum[0] = last_resnum[0] = '\0';

  /* Determine height and width of current picture */
  height = limit.max_y - limit.min_y;
  width = limit.max_x - limit.min_x;

  /* Calculate default height */
  default_height = PICTURE_HEIGHT * DEFAULT_WIDTH / PICTURE_WIDTH;

  /* Calculate aspect ratio for picture and for frame within picture to
     be placed */
  picture_aspect = height / width;
  frame_aspect = PICTURE_HEIGHT / PICTURE_WIDTH;

  /* Determine optimal scaling factor for the picture */
  if (picture_aspect > frame_aspect)
    {
      /* Bounded by y: calculate the scaling factor */
      if (height < default_height)
        height = default_height;
      scale_factor = PICTURE_HEIGHT / height;

      /* Determine origin of central point */
      origin[0] = (PICTURE_WIDTH - width * scale_factor) /2;
      origin[1] = 0.0;
    }
  else
    {
      /* Bounded by x: calculate the scaling factor */
      if (width < DEFAULT_WIDTH)
        width = DEFAULT_WIDTH;
      scale_factor = PICTURE_WIDTH / width;

      /* Determine origin of central point */
      origin[0] = 0.0;
      origin[1] = (PICTURE_HEIGHT - height * scale_factor) / 2;
    }

  /* Determine the picture origin */
  origin[0] = origin[0] + XMARGIN;
  origin[1] = origin[1] + YMARGIN;

  /* Form name of output PostScript file */
  strcat(pspage_ptr->ps_name,"topoldia.ps");

  /* Initialise PostScript page variables */
  init_pspage(psparams_ptr,pspage_ptr);

  /* Get the start and end residues */
  if (params.use_data == FALSE)
    {
      strcpy(first_resnum,residue_array_ptr[0].res_num);
      strcpy(last_resnum,residue_array_ptr[nres - 1].res_num);
    }

  /* Start the page holding the plot */
  start_new_page(psparams_ptr,first_resnum,last_resnum,pspage_ptr,
                 params);

  /* Plot all the strands in the plot */
  plot_sses(first_sse_ptr,scale_factor,origin,pspage_ptr,limit,params);

  /* Close the current page */
  close_page(psparams_ptr,pspage_ptr);

}
/***********************************************************************

write_save_file  -  Write out plot data to .plt file

***********************************************************************/

void write_save_file(struct sse *first_sse_ptr,struct parameters params)
{
  char filename[FILENAME_LEN];

  struct sse *sse_ptr;

  FILE *file_ptr;

  /* Form the name of the save file, topoldia.plt */
  strcpy(filename,"topoldia.plt");

  /* Open the output file */
  file_ptr = open_output_file(filename,TRUE);

  /* Write out the header record giving chain and domain */
  fprintf(file_ptr,"Chain [%c] domain %d\n",params.chain,params.domain);

  /* Get pointer to first sse record */
  sse_ptr = first_sse_ptr;

  /* Loop through sse records to write out */
  while (sse_ptr != NULL)
    {
      /* Write only if required */
      if (sse_ptr->type != OVERHANG1 && sse_ptr->type != OVERHANG2)
        fprintf(file_ptr,"%d %d %d %d %d %d %d %d %s %s\n",sse_ptr->type,
                sse_ptr->xpos[0],sse_ptr->ypos[0],sse_ptr->xpos[1],
                sse_ptr->ypos[1],sse_ptr->first_pos,sse_ptr->last_pos,
                sse_ptr->endseg,sse_ptr->first_resnum,
                sse_ptr->last_resnum);

      /* Get next sse record */
      sse_ptr = sse_ptr->next_sse_ptr;
    }

  /* Close the output file */
  fclose(file_ptr);
}

/***********************************************************************

                             M   A   I   N

***********************************************************************/

int main(int argc,char *argv[])
{
  char pdbsum_dir[FILENAME_LEN];

  int nhelint, npts[2], nres, nsse, origin[2], top_sheet;
  int dir_type;

  int *array;

  struct c_file_data file_name_ptr;
  struct limits limit;
  struct parameters params;

  struct helint *first_helint_ptr;
  struct promotif *first_promotif_ptr;
  struct psparams psparams_ptr;
  struct pspage pspage_ptr;
  struct residue *residue_array_ptr;
  struct sse *first_sse_ptr, *last_sse_ptr;

  struct sse **sse_index_ptr;

  /* Initialise variables */
  limit.first = TRUE;
  nsse = 0;
  strcpy(pdbsum_dir,"./");

  /* Initialise pointers */
  first_promotif_ptr = NULL;

  /* Interpret the command-line parameters */
  get_command_arguments(argv,argc - 1,&params);

  /* Get the paths */
  c_get_paths(&file_name_ptr,params.run_64);

  /* Get the name of the appropriate PDBsum directory */
  if (params.local == FALSE)
    {
      dir_type = find_pdbsum_dir(params.pdb_code,
                                 file_name_ptr.pdbsum_template,
                                 params.pdb_type,pdbsum_dir);
      if (dir_type == NOT_FOUND)
        {
          printf("*** ERROR. Unable to locate PDBsum directory for %s\n",
                 params.pdb_code);
          exit(-1);
        }
    }

  /* Reading in data from .plt file if have file defined */
  if (params.use_data == TRUE)
    {
      /* Read in the list of SSEs and their locations */
      nsse = read_save_file(&first_sse_ptr,&params);

      /* If helix and strand labels are required, read in the Promotif
         data to pick up the appropriate data */
      if (params.show_sse_labels == TRUE)
        {
          /* Read in the PROMOTIF secondary structure assignments */
          get_promotif_data(params.pdb_code,pdbsum_dir,&first_promotif_ptr,
                            &first_helint_ptr,&top_sheet,params);

          /* Transfer the helix/strand labels to the corresponding
             secondary structure elements */
          transfer_labels(first_promotif_ptr,first_sse_ptr,
                          params.chain,top_sheet);
        }
    }

  /* Otherwise, get data from hera.dat file and generate plot from that */
  else
    {
      /* Read in the PROMOTIF secondary structure assignments */
      nhelint
        = get_promotif_data(params.pdb_code,pdbsum_dir,&first_promotif_ptr,
                            &first_helint_ptr,&top_sheet,params);

      /* Read in the layout of the SSEs from the HERA output */
      nsse = get_hera_data(&residue_array_ptr,&nres,&first_sse_ptr,
                           &last_sse_ptr,first_promotif_ptr,&params);

      /* If no data to process, then end here */
      if (nsse == 0)
        {
          printf("*** Warning. Ne hera data to process\n");
          exit(1);
        }

      /* Add extra space between SSEs in the same column */
      separate_sses(first_sse_ptr);

      /* If using the new algorithm, lay out any interacting helices */
      if (params.new == TRUE)
        {
          /* Calculate the layout of the SSEs */
          layout_helices(first_sse_ptr,nsse,first_helint_ptr,nhelint);
        }

      /* Create a sort index for the SSEs */
      create_sort_index(first_sse_ptr,&sse_index_ptr,nsse);

      /* Sort the SSEs according to order in sequence */
      sort_sses(sse_index_ptr,nsse);

      /* Determine size of picture taken up by the SSEs */
      limit.first = TRUE;
      get_max_min(first_sse_ptr,&limit);

      /* Create an array for storing locations of all SSEs for use in
         adding the joins */
      create_mask_array(&array,npts,origin,limit,sse_index_ptr,nsse);

      /* Join SSE ends by straight-line loops */
      nsse = add_joins(sse_index_ptr,nsse,&first_sse_ptr,&last_sse_ptr,
                       residue_array_ptr,nres,limit,array,npts,origin);
    }

  /* Determine maximum and minimum extent of the picture */
  limit.first = TRUE;
  get_max_min(first_sse_ptr,&limit);

  /* Initialise PostScript parameters */
  initialise_psparams(&psparams_ptr,&pspage_ptr);

  /* Plot the topology diagram */
  write_postscript(first_sse_ptr,residue_array_ptr,nres,psparams_ptr,
                   &pspage_ptr,limit,params);

  /* Write out plot data to .plt file, if required */
  if (params.save_data == TRUE)
    write_save_file(first_sse_ptr,params);

  return(0);
}
