/***********************************************************************

  dimhtml.c - Program to read in the output from dimplot/ligplot and
              generate an HTML table showing the interactions across the
              given dimer interface

************************************************************************

Date:-         28 Oct 2004
Dir update:-   19 Oct 2012

v.5.4.4        26 Jan 2022

Written by:-   Roman Laskowski

*** To do:

  1. If interface is massive (like 2hydAB), limit number of residues
     to ~50.

Test run as:

    dimhtml 2hydAB -dir . -use -res GLY481:A -disastr $DISASTR_REAL_UNIPROT_DIR/97/P33897/

***************

Testing DisaStr plot in $TMP_DIR/dim:

     ~/src/ligplot/dimhtml 1aiyBH -dir . -use -res Ala14:B -disastr
         /ebi/research/thornton/data/DisaStr/08/P01308

v.5.4          Inclusion of salt-bridges.
                                                      28 Jan 2020 (RAL)
v.5.4.4        Reduce numbers of minimization loops for very large
               interfaces to speed up.
               Adjustment to line widths and text sizes for DIMPLOT output.
                                                      25 Jan 2022 (RAL)
v.5.4.5        Addition of PDBsum1 options.
                                                      15 Mar 2022 (RAL)
v.5.4.6        Bug-fix in read_save_file.
                                                      16 Jan 2023 (RAL)

----------------------------------------------------------------------*/

/* Datafiles
   ---------

Actual name      Description
-----------      -----------
CATHPARAM        Input parameter file
ligplot.rcm      Input file listing the CofM positions of the residues
                 in both chains
ligplot.hhb      Input file listing the H-bonds across the interface
ligplot.nnb      Input file listing the non-bonded contacts across the
                 interface
dimhtmlXY.dat    Input/output file listing the H-bonds and non-bonded
                 interactions across the interface
dimplot.rcm      Output file listing centre-of-mass restraints for each
                 residue, for use by LIGPLOT when generating DIMPLOT
dimhtml.rcm      As dimplot.rcm, for use by the ligfit program
XXXXc.gradesPE   Input file giving ConSurf residue conservation scores
dimplot.pdb      Input file created by dimer.c when run as part of the
                 DIMPLOT script

Compilation
-----------

cc -c dimhtml.c
cc -c dimps.c
cc -c common.c
cc -c reapar.c
cc -o dimhtml dimhtml.o dimps.o common.o reapar.o -lm -lz

Optimization flags:-  

cc -O2 -funroll-loops -c dimhtml.c
cc -O2 -funroll-loops -c dimps.c
cc -O2 -funroll-loops -c common.c
cc -O2 -funroll-loops -c reapar.c
cc -o dimhtml dimhtml.o dimps.o common.o reapar.o -lm -lz

For DIMPLOT script run as

   dimhtml none -dimp
ligplot dimplot.pdb -prm C:\Roman\java\LigPlus\lib\params\dimplot.prm -ctype 1
*/

/* Function calling-tree
   ---------------------

main
   -> get_command_arguments
         -> store_string(in common.c)
   -> c_get_paths (in common.c)
         -> get_param (in reapar.c)
         -> store_string(in common.c)
   -> find_pdbsum_dir (in common.c)
   -> initialise_psparams
   -> read_save_file
         -> string_truncate (in common.c)
         -> find_residue_byno
         -> create_residue_record
               -> get_residue_colour
         -> create_link_record
   -> get_random_number (in common.c)
   -> scan_rcm_file
         -> string_truncate (in common.c)
   -> read_rcm_file
         -> string_truncate (in common.c)
         -> create_residue_record
               -> get_residue_colour
   -> read_dimplot_pdb
         -> find_residue
         -> string_truncate (in common.c)
   -> sort_residues
         -> shell_sort
   -> read_hhb_file
         -> find_residue
         -> find_link
         -> create_link_record
   -> add_water_bridges
         -> find_link
         -> create_link_record
   -> delete_residues
   -> get_consurf_scores
         -> get_tokens (in common.c)
         -> format_res_num
         -> get_colour_conserv
   -> place_residues
         -> create_alignment_arrays
         -> init_arrays
         -> store_ordered_alignment
         -> nw_sweep
         -> find_best_path
         -> store_alignment
         -> improve_alignment
               -> fill_scores_array
               -> make_random_adjustment
               -> calc_alignment_score
               -> judge_move
                     -> get_random_number (in common.c)
               -> undo_adjustment
         -> close_up_gaps
         -> close_up_more
   -> plot_diagram
         -> init_pspage
               -> adjust_for_landscape
         -> start_new_page
               -> psendp_
               -> psclos_
               -> pspage_
                     -> psrot_
               -> psopen_
         -> get_chain_col
               -> get_col_rgb
         -> get_domain_colour
         -> pscolb_
         -> psrtxt_
         -> psltxt_
         -> plot_interactions
               -> draw_interactions
                     -> get_angle
                     -> psline_
         -> store_variants
               -> extract_residue_name
                     -> to_upper (in common.c)
                     -> format_res_number
               -> find_pdbsum_dir (in common.c)
               -> get_other_chain
                     -> string_truncate (in common.c)
               -> create_highres_entry
               -> get_variants
                     -> string_truncate
                     -> format_res_number
                     -> find_highres_entry
                     -> create_highres_entry
         -> trim_plot
         -> plot_ellipses
               -> get_residue_colour
               -> psccol_
               -> pselip_
               -> convert_residue_name
               -> psctxt_
               -> psrtxt_
               -> psltxt_
         -> close_page
   -> write_save_file
         -> open_output_file (in common.c)
   -> write_dimplot_rcm
         -> write_rcm_details

*/

#include "common.h"
#include "dimhtml.h"
#include "dimps.h"

/* Prototypes
   ---------- */

/* In common.c */
void c_get_paths(struct c_file_data *file_name_ptr,int run_64);
int get_tokens(char *line,char **token,int max_token,char token_sep);
int find_pdbsum_dir(char *pdb_code,char *pdbsum_template[NPDBSUM_DIR],
                    int pdb_type,char *dir_name);
float get_random_number(int *iseed,int random_start);
FILE *open_output_file(char *out_name,int terminate);
void store_string(char **string_ptr,char *string);
int string_truncate(char *string,int max_length);
void to_upper(char *search_string,int length);

/* Function prototypes from dimps.c */
void initialise_psparams(struct psparams *psparams_ptr,
                         struct pspage *pspage_ptr);
float plot_diagram(struct residue **alignment,int align_len,
                   int total_len,int *max_count,
                   struct psparams psparams_ptr,
                   struct pspage *pspage_ptr,char *pdbsum_dir,
                   struct parameters params);

/* Function prototypes from reapar.c */
int get_param(char *token_name,char *file_name,int filename_len,
              int exit_on_error);

/***********************************************************************

get_command_arguments  -  Interpret the command line arguments

***********************************************************************/

void get_command_arguments(char *strptrs[],int ntoken,
                           struct parameters *params)
{
  int itoken, number;

/* v.5.3.6--> */
  float size;
/* <--v.5.3.6 */

  /* Initialise parameters */
  params->pdb_code[0] = '\0';
  params->chain[0] = '\0';
  params->chain[1] = '\0';
  params->domain[0] = 0;
  params->domain[1] = 0;
  params->chains = TRUE;
  params->dimplot = FALSE;
  params->pdb_type = STANDARD_PDB;
  params->colour_scheme = RESIDUE_TYPE;
  params->use_data = FALSE;
  params->save_data = FALSE;
  params->work_dir[0] = '\0';
  params->labels_outside = FALSE;
  params->show_chains = FALSE;
  params->order = FALSE;
  params->flip_seqs = FALSE;
  params->include_waters = FALSE;
  params->contact_type = HYDROPHOBIC_ONLY;
/* v.5.3.6--> */
  params->size_hydrophobic = SIZE_HYDROPHOBIC;
  params->portrait = FALSE;
/* <--v.5.3.6 */
/* v.5.3.8--> */
  params->highlight_res[0] = '\0';
  params->disastr_dir = &null;
/* <--v.5.3.8 */
  params->run_64 = FALSE;

  /* If -help entered on the command line, show options available */
  if (ntoken < 2 || (ntoken == 1 && (!strcmp(strptrs[1],"-help") ||
                                     !strcmp(strptrs[1],"-h"))))
    {
      printf("\n");
      printf("Correct usage is ...\n");
      printf("\n");
      printf("    dimhtml  pdb_codeXY|none  [UPLOAD | MCSG]  [-cons]");
      printf("  [-dimp]  [-dir work_dir]\n");
      printf("             [-use | -save]  [-order]  [-flip]  [-ctype n]"
/* v.5.3.6--> */
//             "  [-waters]\n");
             "  [-waters]  [-hs size]  [-portrait]\n");
/* <--v.5.3.6 */
/* v.5.3.8--> */
/* v.5.4.5 --> */
//      printf("             [-res res_id]  [-disastr dir]\n");
      printf("             [-res res_id]  [-disastr dir]  [-pdbsum1]\n");
/* <-- v.5.4.5 */
/* <--v.5.3.8 */
      printf("where\n");
      printf("     * pdb_code = PDB code for file to be processed\n");
      printf("     * XY       = Chains defining interface (for DIMPLOT "
             "scripts pdb_codeXY has\n");
      printf("                  to be set to \"none\")\n");
      printf("     * UPLOAD   = Uploaded PDB file for PDBsum generation\n");
      printf("     * MCSG     = MCSG structure\n");
      printf("     * -cons    = Colour residues by residue conservation\n");
      printf("     * -dimp    = Read in DIMPLOT files rather than files ");
      printf("output by newsec\n");
      printf("     * -dir work_dir = Directory containing input and ");
      printf("output files\n");
      printf("     * -save    = Save data in dimhtmlXY.dat file in ");
      printf("PDBsum directory\n");
      printf("     * -use     = Read data from dimhtmlXY.dat file in ");
      printf("PDBsum directory\n");
      printf("     * -labout  = Force plot of labels outside residue ");
      printf("ellipses\n");
/* v.5.3.8--> */
      printf("     * -res res_id  = Highlight residue for DisaStr\n");
      printf("     * -disastr dir = DisaStr directory containing variants "
             "data\n");
/* <--v.5.3.8 */
      printf("     * -order   = Order the residues on the first interface "
             "by residue number\n");
      printf("     * -flip    = Flip sequences round - for use with "
             "LigPlus\n");
      printf("     * -ctype   = Type of hydrophobic contacts to use: 0="
             "Hphobic only 1=Hphobic-any\n");
      printf("                  2=Any-any\n");
      printf("     * -waters  = Include waters in calculating alignment\n");
/* v.5.3.6--> */
      printf("     * -hs size = Size of hydrophobic residue representations\n");
      printf("     * -portrait = Reverse order of residues for portrait "
             "plot\n");
/* <--v.5.3.6 */
/* v.5.4.5 --> */
      printf("     * -pdbsum1  = PDBsum1 structure\n");
/* <-- v.5.4.5 */
      printf("     * -64      = run on 64-bit processor\n");
      printf("\n");
      exit(1);
    }

  /* Take first parameter to be the PDB code */
  strncpy(params->pdb_code,strptrs[1],4);
  params->pdb_code[4] = '\0';
  if (strlen(strptrs[1]) > 5)
    {
      params->chain[0] = strptrs[1][4];
      params->chain[1] = strptrs[1][5];
    }

  /* Loop through any additional command arguments */
  for (itoken = 2; itoken < ntoken + 1; itoken++)
    {
      /* Check for the -dir option giving the name of the working
         directory */
      if (!strcmp(strptrs[itoken],"-dir"))
        {
          /* Get the directory name from the next token, if there is one */
          if (itoken < ntoken)
            {
              strcpy(params->work_dir,strptrs[itoken + 1]);
              itoken++;
            }          
        }

      /* Check for the UPLOAD option */
      else if (!strcmp(strptrs[itoken],"UPLOAD"))
        params->pdb_type = UPLOADED_PDB;

      /* Check for the MCSG option */
      else if (!strcmp(strptrs[itoken],"MCSG"))
        params->pdb_type = MCSG_PDB;

      /* Check for the PDBsum1 option */
      else if (!strcmp(strptrs[itoken],"-pdbsum1"))
        params->pdb_type = PDBSUM1;

      /* Check for the -cons option */
      else if (!strcmp(strptrs[itoken],"-cons"))
        params->colour_scheme = CONSERVATION;

      /* Check for the -dimp option */
      else if (!strcmp(strptrs[itoken],"-dimp"))
        params->dimplot = TRUE;

      /* Check for the -save option */
      else if (!strcmp(strptrs[itoken],"-save"))
        params->save_data = TRUE;

      /* Check for the -use option */
      else if (!strcmp(strptrs[itoken],"-use"))
        params->use_data = TRUE;

      /* Check for the -labout option */
      else if (!strcmp(strptrs[itoken],"-labout"))
        params->labels_outside = TRUE;

      /* Check for the -order option */
      else if (!strcmp(strptrs[itoken],"-order"))
        params->order = TRUE;

      /* Check for the -flip option */
      else if (!strcmp(strptrs[itoken],"-flip"))
        params->flip_seqs = TRUE;

      /* Check for the -waters option */
      else if (!strcmp(strptrs[itoken],"-waters"))
        params->include_waters = TRUE;

      /* Check for the -ctype option defining the non-bonded contacts to
         be included */
      else if (!strcmp(strptrs[itoken],"-ctype"))
        {
          /* Get the contact type from the next token, if there is one */
          if (itoken < ntoken)
            {
              number = atoi(strptrs[itoken + 1]);
              if (number > -1 && number < 3)
                params->contact_type = number;
              itoken++;
            }          
        }

/* v.5.3.8--> */
      /* Check for the -res option identifying residue to be
         highlighted */
      else if (!strcmp(strptrs[itoken],"-res"))
        {
          /* Retrieve the residue name from the next token */
          if (itoken < ntoken)
            {
              /* Increment token-number and get next token */
              strcpy(params->highlight_res,strptrs[itoken + 1]);
              itoken++;
            }
        }

      /* Check for the -disastr option giving relevant DisaStr 
         directory */
      else if (!strcmp(strptrs[itoken],"-disastr"))
        {
          /* Retrieve the directory name from the next token */
          if (itoken < ntoken)
            {
              /* Increment token-number and get next token */
              store_string(&(params->disastr_dir),strptrs[itoken + 1]);
              itoken++;
            }
        }
/* <--v.5.3.8 */

/* v.5.3.6--> */
      /* Check for the -hs option defining hydrophobic size */
      else if (!strcmp(strptrs[itoken],"-hs"))
        {
          /* Get the size from the next token, if there is one */
          if (itoken < ntoken)
            {
              size = atof(strptrs[itoken + 1]);
              if (size > 0)
                params->size_hydrophobic = size;
              itoken++;
            }          
        }

      /* Check for the -portrait option */
      else if (!strcmp(strptrs[itoken],"-portrait"))
        params->portrait = TRUE;
/* <--v.5.3.6 */

      /* Check for the 64-bit option */
      else if (!strcmp(strptrs[itoken],"-64"))
        params->run_64 = TRUE;
    }
}
/***********************************************************************

scan_rcm_file  -  Scan the ligplot.rcm file to determine the two chains
                  involved

***********************************************************************/

int scan_rcm_file(char *chains,struct parameters params)
{
  char input_line[LINELEN + 1];
  char file_name[FILENAME_LEN];
  char chain;

  int len, nchains;

  FILE *file_ptr;

  /* Initialise variables */
  chain = ' ';
  chains[0] = '\0';
  chains[1] = '\0';
  file_name[0] = '\0';
  nchains = 0;

  /* Form name of the .rcm file */
  if (params.work_dir[0] != '\0')
    {
      strcpy(file_name,params.work_dir);
      strcat(file_name,DIR_SEP);
    }
  strcat(file_name,"ligplot.rcm");

  /* Open the .rcm file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    {
      printf("*** ERROR. Unable to open ligplot.rcm file: %s\n",file_name);
      exit(-1);
    }

  /* Loop while reading in records from the file */
  while (fgets(input_line,LINELEN,file_ptr) != NULL && nchains < 2)
    {
      /* Truncate the string to strip off the newline */
      len = string_truncate(input_line,LINELEN);

      /* If this is a residue line, extract the chain identifier */
      if (!strncmp(input_line,"RESDUE",6))
        {
          /* Get chain */
          chain = input_line[21];

          /* If no chains so far, store this as the first */
          if (nchains == 0)
            {
              /* Store chain and increment count */
              chains[0] = chain;
              nchains++;
            }

          /* Otherwise, check whether we already have this chain */
          else if (chains[0] != chain)
            {
              /* Store chain and increment count */
              chains[1] = chain;
              nchains++;
            }
        }
    }

  /* Close the file */
  fclose(file_ptr);

  /* If file doesn't contain 2 chains, then abort */
  if (nchains != 2)
    {
      printf("*** ERROR. File ligplot.rcm does not contain 2 chains\n");
      exit(-1);
    }

  /* Check chain order, swapping if necessary */
  if (chains[0] > chains[1])
    {
      /* Swap chains */
      chain = chains[0];
      chains[0] = chains[1];
      chains[1] = chain;
    }

  /* Return number of chains */
  return(nchains);
}
/***********************************************************************

get_residue_colour  -  Get the colour of the given residue, according to
                       its residue type

***********************************************************************/

void get_residue_colour(char *res_name,char *colour)
{
#define NRES_TYPE    25

  int iamino, ires;

  static char *colour_name[] = { 
    "grey", "yellow", "red", "red", "purple", "orange", "sky_blue", "grey",
    "sky_blue", "grey", "grey", "green", "orange", "green", "sky_blue",
    "green", "green", "grey", "purple", "purple", "cyan", "red", "red",
    "cyan", "cyan" };
  static char *amino_name[] = {
    "ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU",
    "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR",
    "PCA ","ASX ","GLX ", "UNK", "XXX"
  };

  /* Find the position of this residue in the list */
  iamino = NRES_TYPE - 1;

  /* Find the residue name */
  for (ires = 0; ires < NRES_TYPE; ires++)
    if (!strcmp(res_name,amino_name[ires]))
      iamino = ires;
  strcpy(colour,colour_name[iamino]);
}
/***********************************************************************

create_residue_record  -  Create and initialise a new residue record

***********************************************************************/

struct residue *create_residue_record(struct residue **fst_residue_ptr,
                                      struct residue **lst_residue_ptr,
                                      char chain,char *res_name,
                                      char *res_num,int chain_no,
                                      double ycoord,int number)
{
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

  /* Initialise the elements of the structure */
  strcpy(residue_ptr->res_name,res_name);
  strcpy(residue_ptr->res_num,res_num);
  residue_ptr->chain = chain;
  residue_ptr->chain_no = chain_no;
  residue_ptr->ycoord = ycoord;
  residue_ptr->number = number;
  residue_ptr->conservation = 0;
  residue_ptr->ref_pos = -1;
  residue_ptr->nhbonds = 0;
/* v.5.3.8--> */
  residue_ptr->deleted = FALSE;
/* <--v.5.3.8 */
  residue_ptr->align_residue_ptr = NULL;
  residue_ptr->first_link_ptr = NULL;
  residue_ptr->last_link_ptr = NULL;

  /* Get this residue's colour */
  get_residue_colour(res_name,residue_ptr->colour);

  /* Initialise pointer to the next residue */
  residue_ptr->next_residue_ptr = NULL;

  /* Return current pointers */
  *fst_residue_ptr = first_residue_ptr;
  *lst_residue_ptr = last_residue_ptr;

  /* Return new residue pointer */
  return(residue_ptr);
}
/***********************************************************************

read_rcm_file  -  Read in the residues and their y-coords form the
                  ligplot.rcm file

***********************************************************************/

int read_rcm_file(struct residue **fst_residue_ptr,int *nresid,
                  struct parameters params)
{
  char input_line[LINELEN + 1], number_string[20];
  char file_name[FILENAME_LEN];
  char chain, res_name[4], res_num[6];

  int chain_no, len, nchains, nresidues;

  double ycoord;

  struct residue *residue_ptr, *first_residue_ptr, *last_residue_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  chain = ' ';
  file_name[0] = '\0';
  nchains = 0;
  nresidues = 0;
  nresid[0] = 0;
  nresid[1] = 0;
  nresid[2] = 0;
  res_name[0] = '\0';
  res_num[0] = '\0';

  /* Initialise pointers */
  residue_ptr = NULL;
  first_residue_ptr = NULL;
  last_residue_ptr = NULL;

  /* Form name of the .rcm file */
  if (params.work_dir[0] != '\0')
    {
      strcpy(file_name,params.work_dir);
      strcat(file_name,DIR_SEP);
    }
  strcat(file_name,"ligplot.rcm");

  /* Open the .rcm file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    {
      printf("*** ERROR. Unable to open ligplot.rcm file: %s\n",file_name);
      exit(-1);
    }

  /* Loop while reading in records from the file */
  while (fgets(input_line,LINELEN,file_ptr) != NULL)
    {
      /* Truncate the string to strip off the newline */
      len = string_truncate(input_line,LINELEN);

      /* If this is a residue line, extract and store the details */
      if (!strncmp(input_line,"RESDUE",6))
        {
          /* Get residue name and number */
          strncpy(res_name,input_line+17,3);
          res_name[3] = '\0';
          strncpy(res_num,input_line+22,5);
          res_num[5] = '\0';
          chain = input_line[21];
          if (chain == '-')
            chain = ' ';

          /* Determine which of the two chains this is */
          if (chain == params.chain[0])
            chain_no = 0;
          else
            chain_no = 1;

          /* Get the y-coordinate of this residue on the DIMPLOT */
          strcpy(number_string,input_line+38);
          ycoord = atof(number_string);

          /* Create residue record */
          residue_ptr
            = create_residue_record(&first_residue_ptr,&last_residue_ptr,
                                    chain,res_name,res_num,chain_no,
                                    ycoord,nresidues);

          /* Increment count of residues stored */
          nresidues++;
          nresid[chain_no]++;
        }
    }

  /* Close the file */
  fclose(file_ptr);

  /* Show number of residue read in */
  printf("Number of residues stored:                  %8d\n",nresidues);

  /* Return pointer to first residue */
  *fst_residue_ptr = first_residue_ptr;

  /* Return number of residues read in */
  return(nresidues);
}
/***********************************************************************

find_residue  -  Find the given residue

***********************************************************************/

struct residue *find_residue(struct residue *first_residue_ptr,
                             char *res_name,char *res_num,char chain)
{
  struct residue *residue_ptr, *found_residue_ptr;

  /* Initialise variables */
  found_residue_ptr = NULL;

  /* Get pointer to the first residue */
  residue_ptr = first_residue_ptr;

  /* Loop over all the residues to find the one just read in */
  while (residue_ptr != NULL && found_residue_ptr == NULL)
    {
      /* If have a match, store the pointer */
      if (chain == residue_ptr->chain && 
          !strcmp(res_num,residue_ptr->res_num) &&
          !strcmp(res_name,residue_ptr->res_name))
        found_residue_ptr = residue_ptr;

      /* Get pointer to the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
    }

  /* Return the matching residue pointer */
  return(found_residue_ptr);
}
/***********************************************************************

read_dimplot_pdb  -  Read in the dimplot.pdb file to determine which
                     residues belong to which interface

***********************************************************************/

int read_dimplot_pdb(struct residue **fst_residue_ptr,int *nresid,
                     struct parameters *params)
{
  char input_line[LINELEN + 1], number_string[20];
  char file_name[FILENAME_LEN];
  char chain, res_name[4], res_num[6];
  char last_chain, last_res_name[4], last_res_num[6];

  char *string_ptr;

  int len, nresidues, surface;
  int done;

  float ycoord;

  struct residue *residue_ptr, *first_residue_ptr, *last_residue_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  done = FALSE;
  file_name[0] = '\0';
  last_chain = '\0';
  last_res_num[0] = '\0';
  last_res_name[0] = '\0';
  surface = -1;
  nresid[0] = 0;
  nresid[1] = 0;
  nresid[2] = 0;
  nresidues = 0;
  ycoord = 0;

  /* Initialise pointers */
  residue_ptr = first_residue_ptr = last_residue_ptr = NULL;

  /* Form name of the dimplot.pdb file */
  if (params->work_dir[0] != '\0')
    {
      strcpy(file_name,params->work_dir);
      strcat(file_name,DIR_SEP);
    }
  strcat(file_name,"dimplot.pdb");

  /* Open the file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    {
      printf("*** ERROR. Unable to open file: %s\n",file_name);
      exit(-1);
    }

  /* Loop while reading in records from the file */
  while (fgets(input_line,LINELEN,file_ptr) != NULL && done == FALSE)
    {
      /* Truncate the string to strip off the newline */
      len = string_truncate(input_line,LINELEN);

      /* If this is the first line, determine whether the two surfaces
         are chains or domains */
      if (!strncmp(input_line,"Title:",6))
        {
          /* Determine whether plotting chains */
          string_ptr = strstr(input_line,"chains");
          if (string_ptr != NULL)
            {
              /* Get the two chains */
              params->chain[0] = string_ptr[7];
              params->chain[1] = string_ptr[13];
              params->chains = TRUE;
            }

          /* Otherwise, must be plotting domains */
          else
            {
              /* Set flag */
              params->chains = FALSE;

              /* Get the two domains */
              string_ptr = strstr(input_line,"domains");
              if (string_ptr != NULL)
                {
                  /* Get the first domain number */
                  strcpy(number_string,string_ptr+8);
                  string_ptr = strstr(number_string," and");
                  if (string_ptr != NULL)
                    string_ptr[0] = '\0';
                  params->domain[0] = atoi(number_string);

                  /* Get the second domain number */
                  string_ptr = strstr(input_line," and");
                  if (string_ptr != NULL)
                    strcpy(number_string,string_ptr+5);
                  params->domain[1] = atoi(number_string);
                }
            }
        }

      /* If this is a surface line, change which surface this is */
      else if (!strncmp(input_line,"SURFACE",7))
        surface++;

      /* If this is an ATOM record, determine which residue it belongs to */
      else if (!strncmp(input_line,"ATOM  ",6) ||
               !strncmp(input_line,"HETATM",6))
        {
          /* Extract the residue name, number and chain-id */
          strncpy(res_name,input_line+17,3);
          res_name[3] = '\0';
          strncpy(res_num,input_line+22,5);
          res_num[5] = '\0';
          chain = input_line[21];

          /* Check whether this is a new residue */
          if (chain != last_chain || strcmp(res_num,last_res_num) ||
              strcmp(res_name,last_res_name))
            {
              /* Find this residue */
              residue_ptr
                = find_residue(first_residue_ptr,res_name,res_num,chain);

              /* If not found, create a new residue */
              if (residue_ptr == NULL)
                {
                  /* Create residue record */
                  residue_ptr
                    = create_residue_record(&first_residue_ptr,
                                            &last_residue_ptr,
                                            chain,res_name,res_num,
                                            surface,ycoord,nresidues);

                  /* Increment residue count */
                  nresidues++;
                }

              /* Update surface assignment */
              residue_ptr->chain_no = surface;

              /* Increment count of residues belonging to this
                 surface */
              nresid[surface]++;

              /* Store the current residue number, name and chain ID */
              strcpy(last_res_num,res_num);
              strcpy(last_res_name,res_name);
              last_chain = chain;
            }
        }

      /* If this is the start of the waters, see if waters are wanted */
      else if (!strncmp(input_line,"WATERS",6))
        {
          /* If waters wanted, then save to third 'surface' */
          if (params->include_waters == TRUE)
            surface++;
          else
            done = TRUE;
        }
    }

  /* Close the file */
  fclose(file_ptr);

  /* Return the first and last residue pointers */
  *fst_residue_ptr = first_residue_ptr;

  /* Return the number of residues */
  return(nresidues);
}
/***********************************************************************

shell_sort  -  Shell-sort which sorts the float values into ascending
               order

***********************************************************************/

void shell_sort(int *index,struct residue **residue_index_ptr,int nindex,
                int order)
{
  int endloop, i, indi, indl, iswap, j, k, l, lognb2, m, n, nn;
  int swap;
  float aln2i = (float) 1.4426950, tiny = (float) 1.e-5;
  float calc;
  double dble, dble1;

  struct residue *residue1_ptr, *residue2_ptr;

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

              /* Initialise flag */
              swap = FALSE;

              /* Get the two residues */
              residue1_ptr = residue_index_ptr[index[indi]];
              residue2_ptr = residue_index_ptr[index[indl]];

              /* Get the two residue records */
              if (order == TRUE)
                {
                  if (residue1_ptr->chain < residue2_ptr->chain)
                    swap = TRUE;
                  else if (residue1_ptr->chain == residue2_ptr->chain &&
                           strcmp(residue1_ptr->res_num,residue2_ptr->res_num)
                           < 0)
                    swap = TRUE;
                }
              else if (residue1_ptr->ycoord > residue2_ptr->ycoord)
                swap = TRUE;

              /* If in wrong order, the swap sort-pointers */
              if (swap == TRUE)
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

sort_residues  -  Sort the residues into descending order of y-coord

***********************************************************************/

void sort_residues(struct residue **fst_residue_ptr,int nresid,int order)
{
  int iorder, ipos;

  int *index;

  struct residue *first_residue_ptr, *last_residue_ptr;
  struct residue *residue_ptr;
  struct residue **residue_index_ptr;

  /* Initialise pointers */
  first_residue_ptr = *fst_residue_ptr;
  last_residue_ptr = NULL;

  /* Grab memory for temporary sort arrays */
  index = (int *) malloc(nresid * sizeof(int));
  residue_index_ptr
    = (struct residue **) malloc(nresid * sizeof(struct residue *));
  if (residue_index_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for residue ");
      printf("index array ");
      exit(1);
    }

  /* Get pointer to first residue record */
  residue_ptr = first_residue_ptr;
  ipos = 0;

  /* Loop through residue records to initialise the sort index */
  while (residue_ptr != NULL)
    {
      /* Initialise sort index */
      index[ipos] = ipos;
      residue_index_ptr[ipos] = residue_ptr;

      /* Increment residue counter */
      ipos++;

      /* Get next residue record */
      residue_ptr = residue_ptr->next_residue_ptr;
    }

  /* Initialise residue counter */
  nresid = ipos;

  /* Sort the residues */
  if (nresid > 1)
    shell_sort(index,residue_index_ptr,nresid,order);

  /* Having sorted the residues, rearrange order of the residue-record
     pointers in descending order of residue */
  for (iorder = nresid - 1; iorder > -1; iorder--)
    {
      /* Get the array position of this residue record */
      ipos = index[iorder];

      /* Get the current residue record */
      residue_ptr = residue_index_ptr[ipos];

      /* If this is the first, then make the first record */
      if (iorder == nresid - 1)
        first_residue_ptr = residue_ptr;

      /* Otherwise, get the previous residue record to point to this one */
      else
        last_residue_ptr->next_residue_ptr = residue_ptr;

      /* Blank out current record's next residue pointer */
      residue_ptr->next_residue_ptr = NULL;

      /* Store current record as last encountered */
      last_residue_ptr = residue_ptr;
    }

  /* Free up memory used for temporary sort arrays */
  free(index);
  free(residue_index_ptr);

  /* Return current pointer */
  *fst_residue_ptr = first_residue_ptr;
}
/***********************************************************************

find_residue_byno  -  Find the given residue by reference number

***********************************************************************/

struct residue *find_residue_byno(struct residue *first_residue_ptr,
                                  int ires)
{
  struct residue *residue_ptr, *found_residue_ptr;

  /* Initialise variables */
  found_residue_ptr = NULL;

  /* Get pointer to the first residue */
  residue_ptr = first_residue_ptr;

  /* Loop over all the residues to find the one just read in */
  while (residue_ptr != NULL && found_residue_ptr == NULL && ires > -1)
    {
      /* If have a match, store the pointer */
      if (ires == residue_ptr->number)
        found_residue_ptr = residue_ptr;

      /* Get pointer to the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
    }

  /* Return the matching residue pointer */
  return(found_residue_ptr);
}
/***********************************************************************

find_link  -  Find the given link

***********************************************************************/

struct link *find_link(struct link *first_link_ptr,
                       struct residue **residue_ptr)
{
  struct link *link_ptr, *found_link_ptr;

  /* Initialise variables */
  found_link_ptr = NULL;

  /* Get pointer to the first link */
  link_ptr = first_link_ptr;

  /* Loop over all the links to find the one just read in */
  while (link_ptr != NULL && found_link_ptr == NULL)
    {
      /* If have a match, store the pointer */
      if (residue_ptr[0] == link_ptr->residue_ptr[0] && 
          residue_ptr[1] == link_ptr->residue_ptr[1])
        found_link_ptr = link_ptr;

      /* Get pointer to the next link */
      link_ptr = link_ptr->next_link_ptr;
    }

  /* Return the matching link pointer */
  return(found_link_ptr);
}
/***********************************************************************

create_link_record  -  Create and initialise a new link record

***********************************************************************/

struct link *create_link_record(struct link **fst_link_ptr,
                                struct link **lst_link_ptr,
                                struct residue **residue_ptr)
{
  int ires;

  struct link *link_ptr, *first_link_ptr, *last_link_ptr;

  /* Initialise link pointers */
  link_ptr = NULL;
  first_link_ptr = *fst_link_ptr;
  last_link_ptr = *lst_link_ptr;

  /* Allocate memory for structure to hold link info */
  link_ptr = (struct link *) malloc(sizeof(struct link));
  if (link_ptr == NULL)
    {
      printf("*** Can't allocate memory for struct link\n");
      exit (1);
    }

  /* If this is the very first, save its pointer */
  if (first_link_ptr == NULL)
    first_link_ptr = link_ptr;

  /* Add link from previous link to the current one */
  if (last_link_ptr != NULL)
    last_link_ptr->next_link_ptr = link_ptr;
  last_link_ptr = link_ptr;

  /* Initialise the elements of the structure */
  link_ptr->residue_ptr[0] = residue_ptr[0];
  link_ptr->residue_ptr[1] = residue_ptr[1];
  link_ptr->nextres_link_ptr[0] = NULL;
  link_ptr->nextres_link_ptr[1] = NULL;
  link_ptr->count[0] = 0;
  link_ptr->count[1] = 0;
  link_ptr->count[2] = 0;
  link_ptr->disulphide = FALSE;
  link_ptr->salt_bridge = FALSE;
  link_ptr->water_bridge = FALSE;

  /* Update each residue's link pointers */
  for (ires = 0; ires < 2; ires++)
    {
      /* If this is this residue's first link, store as first */
      if (residue_ptr[ires]->first_link_ptr == NULL)
        residue_ptr[ires]->first_link_ptr = link_ptr;

      /* Otherwise, get previous link to point to this one */
      else
        residue_ptr[ires]->last_link_ptr->nextres_link_ptr[ires]
          = link_ptr;
      residue_ptr[ires]->last_link_ptr = link_ptr;
    }

  /* Initialise pointer to the next link */
  link_ptr->next_link_ptr = NULL;

  /* Return current pointers */
  *fst_link_ptr = first_link_ptr;
  *lst_link_ptr = last_link_ptr;

  /* Return new link pointer */
  return(link_ptr);
}
/***********************************************************************

read_hhb_file  -  Read in the H-bonds from the ligplot.hhb file or the
                  non-bonded interactions from the ligplot.nnb file

***********************************************************************/

void read_hhb_file(struct residue *first_residue_ptr,
                   struct link **fst_link_ptr,struct link **lst_link_ptr,
                   int type,struct parameters params)
{
  char input_line[LINELEN + 1];
  char file_name[FILENAME_LEN];
  char atom_type1[5], atom_type2[5], chain, res_name[4], res_num[6];

/* v.5.4--> */
//  int ires, len, line, nlinks, offset;
  int bond_type, ires, len, line, nlinks, offset;
/* <--v.5.4 */
  int wanted;

  struct link *link_ptr, *first_link_ptr, *last_link_ptr;
  struct residue *residue_ptr[2], *swap_residue_ptr;

  FILE *file_ptr;

  /* Initialise variables */
/* v.5.4--> */
  bond_type = type;
/* <--v.5.4 */
  file_name[0] = '\0';
  line = 0;
  nlinks = 0;

  /* Initialise pointers */
  link_ptr = NULL;
  first_link_ptr = *fst_link_ptr;
  last_link_ptr = *lst_link_ptr;

  /* Form name of the .rcm file */
  if (params.work_dir[0] != '\0')
    {
      strcpy(file_name,params.work_dir);
      strcat(file_name,DIR_SEP);
    }
  if (strcmp(params.pdb_code,"none"))
    strcat(file_name,"ligplot.");
  else
    strcat(file_name,"dimplot.");
  if (type == HBONDS)
    strcat(file_name,"hhb");
  else
    strcat(file_name,"nnb");

  /* Open the file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    {
      printf("*** ERROR. Unable to open file: %s\n",file_name);
      exit(-1);
    }

  /* Loop while reading in records from the file */
  while (fgets(input_line,LINELEN,file_ptr) != NULL)
    {
      /* Truncate the string to strip off the newline */
      len = string_truncate(input_line,LINELEN);
      line++;

      /* Skip the first 3 lines */
      if (line > 3)
        {
          /* Initialise flag */
          wanted = TRUE;

          /* If reading in H-bonds, then want all */
          if (type == HBONDS)
            wanted = TRUE;

          /* Otherwise, check whether this is the right type of 
             non-bonded contact */
          else if (params.contact_type != ANY_ANY)
            {
              /* Get the two atoms */
              strncpy(atom_type1,input_line+12,4);
              atom_type1[4] = '\0';
              strncpy(atom_type2,input_line+33,4);
              atom_type2[4] = '\0';

              /* If only accepting hydrophobic-hydrophobic contacts,
                 check both atoms are either carbon or sulphur */
              if (params.contact_type == HYDROPHOBIC_ONLY)
                {
                  if ((atom_type1[1] != 'C' && atom_type1[1] != 'S') ||
                      (atom_type2[1] != 'C' && atom_type2[1] != 'S'))
                    wanted = FALSE;
                }

              /* If at least one of the atoms has to be a hydrophobic
                 then check for this */
              else if (params.contact_type == HYDROPHOBIC_ANY)
                {
                  if (atom_type1[1] == 'C' || atom_type1[1] == 'S' ||
                      atom_type2[1] == 'C' || atom_type2[1] == 'S')
                    wanted = TRUE;
                  else
                    wanted = FALSE;
                }
            }

          /* If this interaction is wanted, find the two residues
             involved, and store */
          if (wanted == TRUE)
            {
              /* Initialise pointers */
              residue_ptr[0] = NULL;
              residue_ptr[1] = NULL;

              /* Read in the two residues */
              offset = 0;
              for (ires = 0; ires < 2; ires++)
                {
                  /* Get the residue name and number */
                  strncpy(res_name,input_line+offset,3);
                  res_name[3] = '\0';
                  strncpy(res_num,input_line+offset+6,5);
                  res_num[5] = '\0';
                  chain = input_line[offset + 4];

                  /* Find this residue */
                  residue_ptr[ires]
                    = find_residue(first_residue_ptr,res_name,res_num,chain);

                  /* Adjust the line offset for the next residue */
                  offset = 21;
                }

/* v.5.4--> */
              /* See if this is a salt bridge or a disulphide bond */
              if (len > 49 && !strncmp(input_line + 47,"+/-",3))
                bond_type = SALT_BRIDGES;
              else if (len > 49 && !strncmp(input_line + 47,"S-S",3))
                bond_type = DISULPHIDES;
              else
                bond_type = type;
/* <--v.5.4 */

              /* If have found both residues, form a link between them */
              if (residue_ptr[0] != NULL && residue_ptr[1] != NULL)
                {
                  /* Check that the residues are the right way round */
                  if (residue_ptr[0]->chain_no > residue_ptr[1]->chain_no)
                    {
                      /* Residues wrong way round, so swap */
                      swap_residue_ptr = residue_ptr[0];
                      residue_ptr[0] = residue_ptr[1];
                      residue_ptr[1] = swap_residue_ptr;
                    }

                  /* Check whether we already have this link */
                  link_ptr = find_link(first_link_ptr,residue_ptr);

                  /* Create a link record linking these two residues */
                  if (link_ptr == NULL)
                    {
                      /* Create link */
                      link_ptr
                        = create_link_record(&first_link_ptr,&last_link_ptr,
                                             residue_ptr);

                      /* Increment count of link records */
                      nlinks++;
                    }

/* v.5.4--> */
                  /* If this is a salt bridge, set flag */
                  if (bond_type == SALT_BRIDGES)
                    link_ptr->salt_bridge = TRUE;
                  else if (bond_type == DISULPHIDES)
                    link_ptr->disulphide = TRUE;
/* <--v.5.4 */

                  /* Increment count of occurrences of this link */
/* v.5.4--> */
//                  link_ptr->count[type]++;
                  link_ptr->count[bond_type]++;
/* <--v.5.4 */

                  /* If this is an H-bond, increment count for both
                     residues */
/* v.5.4--> */
//                  if (type == HBONDS)
                  if (bond_type == HBONDS)
/* <--v.5.4 */
                    {
                      residue_ptr[0]->nhbonds++;
                      residue_ptr[1]->nhbonds++;
                    }
                }
            }
        }
    }

  /* Return first and last link pointers */
  *fst_link_ptr = first_link_ptr;
  *lst_link_ptr = last_link_ptr;

  /* Close the file */
  fclose(file_ptr);
}
/***********************************************************************

add_water_bridges  -  Increase scores between residues that are linked
                      via a bridging water

***********************************************************************/

int add_water_bridges(struct residue *first_residue_ptr,
                      struct link **fst_link_ptr,
                      struct link **lst_link_ptr)
{
  int nbridges;

  struct link *link_ptr, *new_link_ptr;
  struct link *first_link_ptr, *last_link_ptr, *other_link_ptr;
  struct residue *residue_ptr[2], *swap_residue_ptr;

  /* Initialise variables */
  nbridges = 0;

   /* Initialise pointers */
  link_ptr = NULL;
  first_link_ptr = *fst_link_ptr;
  last_link_ptr = *lst_link_ptr;

  /* Get the first interaction */
  link_ptr = first_link_ptr;

  /* Loop over all interaction, to process the water H-bonds */
  while (link_ptr != NULL)
    {
      /* If this is a water interaction, then process */
      if (!strcmp(link_ptr->residue_ptr[1]->res_name,"HOH"))
        {
          /* Get the residue this water links to */
          residue_ptr[0] = link_ptr->residue_ptr[0];

          /* Get the next link */
          other_link_ptr = link_ptr->next_link_ptr;

          /* Loop over the remaining links to identify residues on
             either side of the interface that are bridged by this water */
          while (other_link_ptr != NULL)
            {
              /* If this is also a water link, then get its interacting
                 residue */
              if (!strcmp(other_link_ptr->residue_ptr[1]->res_name,"HOH") &&
                  !strcmp(other_link_ptr->residue_ptr[1]->res_num,
                          link_ptr->residue_ptr[1]->res_num) &&
                  other_link_ptr->residue_ptr[1]->chain ==
                  link_ptr->residue_ptr[1]->chain)
                {
                  /* Get the residue this water links to */
                  residue_ptr[1] = other_link_ptr->residue_ptr[0];

                  /* If the two residues are on either side of the
                     interface, then have bridge */
                  if (residue_ptr[0]->chain_no != residue_ptr[1]->chain_no)
                    {
                      /* Check that the residues are the right way round */
                      if (residue_ptr[0]->chain_no >
                          residue_ptr[1]->chain_no)
                        {
                          /* Residues wrong way round, so swap */
                          swap_residue_ptr = residue_ptr[0];
                          residue_ptr[0] = residue_ptr[1];
                          residue_ptr[1] = swap_residue_ptr;
                        }

                      /* Check whether we already have this link */
                      new_link_ptr
                        = find_link(first_link_ptr,residue_ptr);

                      /* Create a link record linking these two residues */
                      if (new_link_ptr == NULL)
                        {
                          /* Create link */
                          new_link_ptr
                            = create_link_record(&first_link_ptr,
                                                 &last_link_ptr,
                                                 residue_ptr);

                          /* Set flag */
                          new_link_ptr->water_bridge = TRUE;

                          /* Increment number of water bridges */
                          nbridges++;
                        }
                    }
                }

              /* Get the next link record */
              other_link_ptr = other_link_ptr->next_link_ptr;
            }
        }

      /* Get the next link record */
      link_ptr = link_ptr->next_link_ptr;
    }

   /* Return pointers */
   *fst_link_ptr = first_link_ptr;
   *lst_link_ptr = last_link_ptr;

  /* Return number of water bridges */
  return(nbridges);
}
/***********************************************************************

delete_residues  -  Delete any residues that have been marked for
                    deletion

***********************************************************************/

int delete_residues(struct residue **fst_residue_ptr,int *nresid,
                    int *nchn)
{
#define MAX_CHAINS     30

  char chain[MAX_CHAINS + 1];

  int i, nchains, ndeleted, nresidues;
  int done;

  struct residue *residue_ptr, *first_residue_ptr, *next_residue_ptr;
  struct residue *prev_residue_ptr;

  /* Initialise variables */
  chain[0] = '\0';
  nchains = 0;
  ndeleted = 0;
  nresid[0] = nresid[1] = nresid[2] = 0;
  nresidues = 0;

  /* Initialise pointers */
  first_residue_ptr = *fst_residue_ptr;

  /* Get first residue */
  residue_ptr = first_residue_ptr;
  prev_residue_ptr = NULL;

  /* Loop through the residue records until done */
  while (residue_ptr != NULL)
    {
      /* Get pointer to next residue record */
      next_residue_ptr = residue_ptr->next_residue_ptr;

      /* If residue to be deleted, remove it from the linked list */
      if (residue_ptr->first_link_ptr == NULL)
        {
          /* If have a previous residue, then get it to point
             around this one */
          if (prev_residue_ptr != NULL)
            prev_residue_ptr->next_residue_ptr = next_residue_ptr;

          /* Otherwise, this must be the first residue, so make
             the next one the new first */
          else
            first_residue_ptr = next_residue_ptr;

          /* Free up the memory taken by the deleted residue */
          free(residue_ptr);

          /* Increment cou7nt of deleted residues */
          ndeleted++;
        }

      /* Otherwise, save as the new previous record */
      else
        {
          /* Save current pointer */
          prev_residue_ptr = residue_ptr;

          /* Increment count of undeleted entries */
          nresidues++;

          /* See if we have this chain already */
          done = FALSE;
          for (i = 0; i < nchains && done == FALSE; i++)
            if (chain[i] == residue_ptr->chain)
              done = TRUE;

          /* If chain not already stored, add it to the list and
             increment count */
          if (done == FALSE && nchains < MAX_CHAINS)
            {
              chain[nchains] = residue_ptr->chain;
              nchains++;
              chain[nchains] = '\0';
            }

          /* Increm,ent number of residues on each side of the interface */
          nresid[residue_ptr->chain_no]++;
        }

      /* Get pointer to next residue record */
      residue_ptr = next_residue_ptr;
    }

  /* Return the first residue pointer */
  *fst_residue_ptr = first_residue_ptr;

  /* Return the number of chains encountered */
  *nchn = nchains;

  /* Return number of residues remaining */
  return(nresidues);
}
/***********************************************************************

read_save_file  -  Read in the details from the dimhtmlXY.dat file

***********************************************************************/

int read_save_file(char *pdbsum_dir,struct residue ***algnmnt,
                   int *aln_len,char *pdb_code,char *chains,
                   int *max_count,struct residue **fst_residue_ptr,
                   int *nr,int *nresid,struct link **fst_link_ptr,
/* v.5.4.5 --> */
//                   char *extension,char *output_dir)
                   char *extension,char *output_dir,int pdb_type)
/* <-- v.5.4.5 */
{
  char filename[FILENAME_LEN], input_line[LINELEN + 1], number_string[20];
/* v.5.4.6 --> */
//  char directory[FILENAME_LEN];
  char directory[FILENAME_LEN], line_end[LINELEN + 1];
/* <-- v.5.4.6 */
  char chain, res_name[4], res_num[6];

  char *string_ptr;

  int align_len, chain_no, len, nresidues, number, total_len;
  int align_pos, i, ipos, ires, jres, ncontacts, nhbonds, nlines, nlinks;
  int disulphide, done, reading, salt_bridge, started;

  double ycoord;

  struct link *link_ptr, *first_link_ptr, *last_link_ptr;
  struct residue *residue_ptr, *first_residue_ptr, *last_residue_ptr;
  struct residue *residue_pair_ptr[2];

  struct residue **alignment;

  FILE *file_ptr;

  /* Initialise variables */
  align_len = 0;
  align_pos = 0;
  alignment = NULL;
  done = FALSE;
  for (i = 0; i < 2; i++)
    max_count[i] = 0;
  nlinks = 0;
  nresidues = 0;
  nresid[0] = 0;
  nresid[1] = 0;
  nresid[2] = 0;
  reading = NOTHING;
  started = FALSE;
  total_len = 0;
  ycoord = 0.0;

  /* Get output directory */
/* v.5.4.5 --> */
//  if (output_dir[0] != '\0')
  if (pdb_type == PDBSUM1)
    strcpy(directory,"../newsec/");
  else if (output_dir[0] != '\0')
/* <-- v.5.4.5 */
    {
      strcpy(directory,output_dir);
      strcat(directory,DIR_SEP);
    }
  else
    directory[0] = '\0';

  /* Initialise pointers */
  first_link_ptr = NULL;
  last_link_ptr = NULL;
  first_residue_ptr = NULL;
  last_residue_ptr = NULL;

  /* If reading .out file, form full name of dimhtml.out in the PDBsum
     directory */
  if (!strcmp(extension,"out"))
    sprintf(filename,"%sdimhtml.%s",pdbsum_dir,extension);

  /* Otherwise, take that we are reading the dimhtmlXY.dat file from
     the working directory */
  else
    {
      /* Define file name */
      sprintf(filename,"%sdimhtml%c%c.%s",directory,chains[0],chains[1],
              extension);

      /* Set flags to start reading data straight away */
      reading = RESIDUE_DATA;
      started = TRUE;
    }

  /* Open the input file */
  if ((file_ptr = fopen(filename,"r")) ==  NULL)
    {
      /* If dimhtml.out file not found, try the old dimhtmlXY.out file */
      if (!strcmp(extension,"out"))
        {
          /* Define the file name */
          sprintf(filename,"%sdimhtml%c%c.%s",pdbsum_dir,chains[0],
                  chains[1],extension);

          /* Set flags to start reading data straight away */
          reading = RESIDUE_DATA;
          started = TRUE;
        }

      /* Otherwise have error */
      else
        {
          printf("*** ERROR. Unable to open file: %s\n",filename);
          exit(-1);
        }

      /* Open the dimhtmlXY.out file */
      if ((file_ptr = fopen(filename,"r")) ==  NULL)
        {
          printf("*** ERROR. Unable to open file: %s\n",filename);
          exit(-1);
        }
    }

  /* Loop while reading in records from the file */
  while (fgets(input_line,LINELEN,file_ptr) != NULL && done == FALSE)
    {
      /* Truncate the string to strip off the newline */
      len = string_truncate(input_line,LINELEN);

      /* If not yet started, check whether this might be the start */
      if (started == FALSE)
        {
          /* If this is the START line, check if it is the right
             interface */
          if (!strncmp(input_line,"START: ",7))
            {
              /* Check for the right chain identifiers */
              if (input_line[14] == chains[0] &&
                  input_line[15] == chains[1])
                {
                  /* Set flags */
                  started = TRUE;
                  reading = RESIDUE_DATA;
                }
            }
        }

      /* Check for end of current interface */
      else if (!strncmp(input_line,"----------",10))
        done = TRUE;

      /* Check for interface key word prior to alignment */
      else if (!strncmp(input_line,"INTERFACE",9))
        {
          /* Set flag */
          reading = ALIGNMENT_DATA;

          /* Get the alignment length */
          strcpy(number_string,input_line+10);
          total_len = atoi(number_string);
          align_len = total_len;
          *aln_len = total_len;

          /* Allocate memory for storing the alignment */
          alignment
            = (struct residue **)
            malloc(align_len * 2 * sizeof(struct residue *));
          if (alignment == NULL)
            {
              printf("*** ERROR. Unable to allocate memory for ");
              printf("alignment array ");
              exit(1);
            }

          /* Initialise the alignment array */
          for (ipos = 0; ipos < 2 * align_len; ipos++)
            alignment[ipos] = NULL;
        }

      /* Check for interactions key word */
      else if (!strcmp(input_line,"INTERACTIONS"))
        reading = INTERACTION_DATA;

      /* If reading in residue data, get details and create residue
         record */
      else if (reading == RESIDUE_DATA)
        {
          /* Get residue details */
          strncpy(res_name,input_line,3);
          res_name[3] = '\0';
          strncpy(res_num,input_line+4,5);
          res_num[5] = '\0';
          chain = input_line[10];
          if (chain == chains[0])
            chain_no = 0;
          else
            chain_no = 1;

          /* Increment count of residues belonging to each chain */
          nresid[chain_no]++;

          /* Get the residue reference number */
          strcpy(number_string,input_line+12);
          string_ptr = strstr(number_string," ");
/* v.5.4.6 --> */
//          if (string_ptr != NULL)
//            string_ptr[0] = '\0';
          line_end[0] = '\0';
          if (string_ptr != NULL)
            {
              strcpy(line_end,string_ptr + 1);
              string_ptr[0] = '\0';
            }
/* <-- v.5.4.6 */
          number = atoi(number_string);

          /* If have a y-coordinate at the end of the line, read that
             in */
          if (line_end[0] != '\0')
            {
              /* Retrieve the y-coordinate for this residue */
/* v.5.4.6 --> */
//              strcpy(number_string,string_ptr+1);
              strcpy(number_string,line_end);
/* <-- v.5.4.6 */
              ycoord = atof(number_string);
            }
          else
            ycoord = 0.0;

          /* Create a residue record */
          residue_ptr
            = create_residue_record(&first_residue_ptr,&last_residue_ptr,
                                    chain,res_name,res_num,chain_no,
                                    ycoord,number);

          /* Increment count of residues stored */
          nresidues++;
        }

      /* If reading in interface data, store residue alignment */
      else if (reading == ALIGNMENT_DATA && align_pos < align_len)
        {
          /* Get the two residue numbers */
          sscanf(input_line,"%d %d",&ires,&jres);


          /* Get the first residue number */
          /*          strcpy(number_string,input_line);
          string_ptr = strstr(number_string," ");
          if (string_ptr != NULL)
            string_ptr[0] = '\0';
            ires = atoi(number_string);*/

          /* Locate the first residue */
          residue_ptr = find_residue_byno(first_residue_ptr,ires);

          /* Store at this position in the alignment */
          alignment[align_pos] = residue_ptr;

          /* Get the second residue number */
          /* string_ptr = strstr(input_line," ");
          if (string_ptr != NULL)
            {
              strcpy(number_string,string_ptr+1);
              ires = atoi(number_string); */

              /* Locate the residue */
              residue_ptr = find_residue_byno(first_residue_ptr,jres);

              /* Store at this position in the alignment */
              alignment[align_pos + align_len] = residue_ptr;
              // }

          /* Increment alignment position */
          align_pos++;
        }

      /* If reading in residue interactionsdata, get details and create residue
         record */
      else if (reading == INTERACTION_DATA)
        {
          /* Get the two residue numbers and their contact numbers */
          sscanf(input_line,"%d %d %d %d",&ires,&jres,&nhbonds,
                 &ncontacts);

          /* Check for markers indicating a salt-bridge or disulphide
             bond */
          string_ptr = strstr(input_line,"s");
          if (string_ptr != NULL)
            salt_bridge = TRUE;
          else
            salt_bridge = FALSE;
          string_ptr = strstr(input_line,"d");
          if (string_ptr != NULL)
            disulphide = TRUE;
          else
            disulphide = FALSE;
          
          /* Locate the two residues */
          residue_pair_ptr[0] = find_residue_byno(first_residue_ptr,ires);
          residue_pair_ptr[1] = find_residue_byno(first_residue_ptr,jres);

          /* If both valid, create a link record between them */
          if (residue_pair_ptr[0] != NULL && residue_pair_ptr[1] != NULL)
            {
              /* Create link */
              link_ptr
                = create_link_record(&first_link_ptr,&last_link_ptr,
                                     residue_pair_ptr);

              /* Store the interaction counts */
              link_ptr->count[HBONDS] = nhbonds;
              link_ptr->count[CONTACTS] = ncontacts;
              link_ptr->disulphide = disulphide;
              link_ptr->salt_bridge = salt_bridge;

              /* Get the number of lines to be drawn between these two
                 residues */
              nlines = nhbonds;
              if (disulphide == TRUE)
                nlines++;
              if (salt_bridge == TRUE)
                nlines++;

              /* Update maximum counts */
              if (nlines > max_count[HBONDS])
                max_count[HBONDS] = nlines;
              if (ncontacts > max_count[CONTACTS])
                max_count[CONTACTS] = ncontacts;

              /* Increment count of link records */
              nlinks++;
            }
        }
    }

  /* Close the input file */
  fclose(file_ptr);

  /* Return pointers */
  *algnmnt = alignment;
  *aln_len = align_len;
  *nr = nresidues;
  *fst_link_ptr = first_link_ptr;
  *fst_residue_ptr = first_residue_ptr;

  /* If maximum number of residues exceeded, then abort */
  if (nresidues > MAX_RESIDUES)
    {
      printf("*** ERROR. Maximum number of residues (%d) exceeded (%d)\n",
             MAX_RESIDUES,nresidues);
      exit(1);
    }

  /* Return the total alignment length */
  return(total_len);
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

get_colour_conserv  -  Get vertex colour from associated residue's
                       conservation score

***********************************************************************/

void get_colour_conserv(double conservation,char *colour)
{
  int iband, ipos;

  static char *colour_name[] = { 
    "white", "blue", "purple", "sky_blue", "cyan", "grey", "green",
    "yellow", "orange", "pink", "red" };

  /* Initialise variables */
  iband = 0;
  strcpy(colour,"white");

  /* Loop until find which band the conservation score falls into */
  for (ipos = 1; ipos < 11 && iband == 0; ipos++)
    {
      if (conservation <= RANGE[ipos])
        iband = ipos;
    }

  /* Otherwise, use colour corresponding to given band */
  strcpy(colour,colour_name[iband]);
}
/***********************************************************************

get_consurf_scores  -  Pick up all the residue conservation scores
                       from this PDB file's ConSurf data file

***********************************************************************/

void get_consurf_scores(struct c_file_data file_name_ptr,char *pdb_code,
                        char chain,struct residue *first_residue_ptr)
{
  char chain2[2];
  char filename[FILENAME_LEN], input_line[LINELEN + 1];
  char res_name[4], res_num[6];
  char last_chain, last_read_chain, last_res_name[4], last_res_num[6];
  char number_string[20];

  char *string_ptr, *token[MAX_TOKEN];

  int len, nlines, nmatched, ntoken;
  int found, get_next_residue, have_file, next_one, reading;

  double conserv;

  struct residue *residue_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  printf("Reading ConSurf data file %s%c ...\n",pdb_code,chain);
  have_file = TRUE;
  get_next_residue = FALSE;
  last_chain = '\0';
  last_read_chain = '\0';
  last_res_name[0] = '\0';
  last_res_num[0] = '\0';
  nlines = 0;
  next_one = FALSE;
  nmatched = 0;
  reading = FALSE;

  /* Form name of input file */
  strcpy(filename,file_name_ptr.other_dir[NEW_CONSURF_DIR]);
  strcat(filename,DIR_SEP);
  strcat(filename,pdb_code);
  if (chain == ' ')
    strcat(filename,"_");
  else
    {
      chain2[0] = chain;
      chain2[1] = '\0';
      strcat(filename,chain2);
    }
  strcat(filename,".gradesPE");

  /* Open input residue conservation scores file */
  if ((file_ptr = fopen(filename,"r")) ==  NULL)
    {
      printf("*** Warning. Unable to open file [%s]\n",filename);
      have_file = FALSE;
    }

  /* Loop while reading in records from the PDB file */
  while (have_file == TRUE && fgets(input_line,LINELEN,file_ptr)!= NULL)
    {
      /* Increment count of lines read in */
      nlines++;

      /* If reading in the residue data, process this line */
      if (reading == TRUE)
        {
          /* Set flag */
          get_next_residue = FALSE;

          /* Unpack the line of tab-separated fields */
          ntoken = get_tokens(input_line,token,MAX_TOKEN,'\t');

          /* If have sufficient data, then store details */
          if (ntoken > 6)
            {
              /* If residue name missing (ie just a hyphen) assume
                 it maps to insertion due to bug in ConSurf */
              if (!strcmp(token[2],"-"))
                {
                  /* Need to match residue following the one we've just
                     processed */
                  get_next_residue = TRUE;
                  res_name[0] = '\0';
                  res_num[0] = '\0';
                }

              /* Otherwise, get the residue details */
              else
                {
                  strncpy(res_name,token[2],3);
                  res_name[3] = '\0';
                  len = strlen(token[2]);
                  strcpy(number_string,token[2]+3);
                  string_ptr = strstr(number_string,":");
                  if (string_ptr != NULL)
                    string_ptr[0] = '\0';
                  strcpy(res_num,number_string);

                  /* Format the residue number, adding the relevant numbers
                     of blanks either side of the number */
                  format_res_num(res_num);
                }

              /* Get the conservation score */
              conserv = atof(token[4]) / 10.0 + 0.01;

              /* Initialise flags */
              found = FALSE;
              next_one = FALSE;

              /* Get pointer to the first of the stored residues */
              residue_ptr = first_residue_ptr;

              /* Loop through all the residues until current one is found */
              while (residue_ptr != NULL)
                {
                  /* Initialise flags */
                  found = FALSE;

                  /* If residue details match, then set flag */
                  if (!strcmp(residue_ptr->res_name,res_name) &&
                      !strcmp(residue_ptr->res_num,res_num) &&
                      residue_ptr->chain == chain)
                    found = TRUE;

                  /* If residue number and chain match, check whether we have
                     a special residue */
                  else if (!strcmp(residue_ptr->res_num,res_num) &&
                           residue_ptr->chain == chain)
                    {
                      /* Check for MSE */
                      if (!strcmp(residue_ptr->res_name,"MSE") &&
                          !strcmp(res_name,"MET"))
                        found = TRUE;
                    }

                  /* If looking for unnamed insertion code residue,
                     check whether this matches the residue prior to
                     it or is the one we're after */
                  if (get_next_residue == TRUE)
                    {
                      /* If this is the next one, then we're done */
                      if (next_one == TRUE)
                        {
                          /* Set flags */
                          found = TRUE;
                          next_one = FALSE;
                        }

                      /* Otherwise, check for match to previous residue */
                      else
                        {
                          /* If residue details match, then set flag */
                          if (!strcmp(residue_ptr->res_name,last_res_name) &&
                              !strcmp(residue_ptr->res_num,last_res_num) &&
                              residue_ptr->chain == chain)
                            next_one = TRUE;

                          /* If residue number and chain match, check
                             whether we have a special residue */
                          else if (!strcmp(residue_ptr->res_num,
                                           last_res_num) &&
                                   residue_ptr->chain == chain)
                            {
                              /* Check for MSE */
                              if (!strcmp(residue_ptr->res_name,"MSE") &&
                                  !strcmp(last_res_name,"MET"))
                                next_one = TRUE;
                            }
                        }
                    }

                  /* If residue found, then store conservation score */
                  if (found == TRUE)
                    {
                      /* Store conversation score */
                      residue_ptr->conservation = conserv;

                      /* Get the corresponding colour */
                      get_colour_conserv(conserv,residue_ptr->colour);

                      /* Store matched residue name and number */
                      strcpy(last_res_name,residue_ptr->res_name);
                      strcpy(last_res_num,residue_ptr->res_num);

                      /* Increment count of matched residues */
                      nmatched++;
                    }

                  /* Get pointer to the next residue */
                  residue_ptr = residue_ptr->next_residue_ptr;
                }
            }
        }

      /* Otherwise check for column titles line signifying the start
         of the data */
      else if (!strncmp(input_line," POS",4))
        reading = TRUE;
    }

  /* Show number of matched residues */
  printf("Number of matched residues = %d\n",nmatched);
}
/***********************************************************************

create_alignment_arrays  -  Create the arrays for the sequence alignment

***********************************************************************/

void create_alignment_arrays(int *nresid,int **arr,int **pth,
                             struct residue *first_residue_ptr,
                             struct residue **residue_index_ptr[3])
{
  int iresid, iseq, max_path, narray;

  int *array, *path;

  struct residue *residue_ptr;

  /* Initialise variables */
  narray = nresid[0] * nresid[1];
  max_path = nresid[0] + nresid[1];
  
  /* Loop over the three sequences */
  for (iseq = 0; iseq < 3; iseq++)
    {
      /* If sequence has any residues, create an index for it */
      if (nresid[iseq] > 0)
        {
          /* Grab memory for array of residue pointers */
          residue_index_ptr[iseq]
            = (struct residue **)
            malloc(nresid[iseq] * sizeof(struct residue *));
          if (residue_index_ptr[iseq] == NULL)
            {
              printf("*** ERROR. Unable to allocate memory for residue"
                     " array\n");
              exit(1);
            }

          /* Get pointer to the first of the residue records */
          residue_ptr = first_residue_ptr;
          iresid = 0;
  
          /* Loop through all the residuees to initialise index array */
          while (residue_ptr != NULL && iresid < nresid[iseq])
            {
              /* If this residues belongs to this sequence, add to array */
              if (residue_ptr->chain_no == iseq)
                {
                  /* Place in array */
                  residue_index_ptr[iseq][iresid] = residue_ptr;

                  /* Increment residue count */
                  iresid++;
                }

              /* Get the next residue record in the list */
              residue_ptr = residue_ptr->next_residue_ptr;
            }
        }
      else
        residue_index_ptr[iseq] = NULL;
    }

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

  /* Return the array pointers */
  *arr = array;
  *pth = path;
}
/***********************************************************************

init_arrays  -  Initialise alignment arrays with interaction data

***********************************************************************/

void init_arrays(int *nresid,struct residue **residue_index_ptr[3],
                int *array,int *max_count)
{
#define MIN_VALUE    0

  int i, ipoint, ires, jres, nlines;
  int ichain, ncontacts, nhbonds;

  int disulphide, found, salt_bridge, water_bridge;

  struct link *link_ptr, *first_link_ptr;
  struct residue *residue_ptr, *residue1_ptr, *residue2_ptr;

  /* Initialise variables */
  for (i = 0; i < 2; i++)
    max_count[i] = 0;

  /* Initialise both chains' alignment and fitting variables */
  for (ichain = 0; ichain < 2; ichain++)
    {
      /* Loop over the residues in this chain */
      for (ires = 0; ires < nresid[ichain]; ires++)
        {
          /* Get the corresponding residue */
          residue_ptr = residue_index_ptr[ichain][ires];

          /* Initialise pointer to aligned residue in other chain */
          residue_ptr->align_residue_ptr = NULL;
        }
    }

  /* Loop over the residues in first chain */
  for (ires = 0; ires < nresid[0]; ires++)
    {
      /* Get the corresponding residue */
      residue1_ptr = residue_index_ptr[0][ires];

      /* Get the first of this residue's interactions */
      first_link_ptr = residue1_ptr->first_link_ptr;

      /* Loop over the residues in the second chain */
      for (jres = 0; jres < nresid[1]; jres++)
        {
          /* Get the corresponding residue */
          residue2_ptr = residue_index_ptr[1][jres];

          /* Get the first link invloving the first residue */
          link_ptr = first_link_ptr;

          /* Initialise count of H-bonds and non-bonded interactions */
          nhbonds = 0;
          ncontacts = 0;
          disulphide = FALSE;
          salt_bridge = FALSE;
          water_bridge = FALSE;
          found = FALSE;

          /* Loop through all the links to see if first residue
             interacts with the second */
          while (link_ptr != NULL && found == FALSE)
            {
              /* If this is the right residue, store numbers of H-bonds
                 and non-bonded contacts */
              if (link_ptr->residue_ptr[1] == residue2_ptr)
                {
                  /* Store numbers of interactions */
                  nhbonds = link_ptr->count[HBONDS];
                  ncontacts = link_ptr->count[CONTACTS];
                  if (link_ptr->disulphide == TRUE)
                    disulphide = TRUE;
                  if (link_ptr->salt_bridge == TRUE)
                    salt_bridge = TRUE;
                  if (link_ptr->water_bridge == TRUE)
                    water_bridge = TRUE;

                  /* Get the number of lines to be drawn between these two
                     residues */
                  nlines = nhbonds;
                  if (disulphide == TRUE)
                    nlines++;
                  if (salt_bridge == TRUE)
                    nlines++;
                  if (water_bridge == TRUE)
                    {
                      nhbonds++;
                      nlines++;
                    }

                  /* Update maximum counts */
                  if (nlines > max_count[HBONDS])
                    max_count[HBONDS] = nlines;
                  if (ncontacts > max_count[CONTACTS])
                    max_count[CONTACTS] = ncontacts;

                  /* Set flag */
                  found = TRUE;
                }

              /* Get the next link */
              link_ptr = link_ptr->nextres_link_ptr[0];
            }

          /* Get the array location of this residue pair */
          ipoint = ires * nresid[1] + jres;

          /* Store the similarity between these two residues */
          array[ipoint] = HBOND_WEIGHT * nhbonds + ncontacts;
        }
    }
}
/***********************************************************************

show_array  -  Print out the scoring matrix (for debugging purposes)

***********************************************************************/

void show_array(int *array,struct residue **residue_index_ptr[3],
                int nresid[2])
{
  int score;
  int ipoint, ires, jres;

  struct residue *residue1_ptr, *residue2_ptr;

  printf("                ");
  for (jres = 0; jres < nresid[1]; jres++)
    {
      /* Get this residue */
      residue2_ptr = residue_index_ptr[1][jres];
      printf("%s%s",
             residue2_ptr->res_name,
             residue2_ptr->res_num);
    }
  printf("\n");
  for (ires = 0; ires < nresid[0]; ires++)
    {
      residue1_ptr = residue_index_ptr[0][ires];
      printf("%2d. %s %s %c",ires,
             residue1_ptr->res_name,
             residue1_ptr->res_num,
             residue1_ptr->chain);

      for (jres = 0; jres < nresid[1]; jres++)
        {
          ipoint = ires * nresid[1] + jres;
          score = array[ipoint];

          printf(" %7d",score);
        }
      printf("\n");
      fflush(stdout);
    }

}
/***********************************************************************

store_ordered_alignment  -  Store the sequence alignment as an array of
                            residue pointers in the case where the order
                            of the first sequence is to be retained

***********************************************************************/

int store_ordered_alignment(int *array,int nresid[2],
                            struct residue **residue_index_ptr[3],
                            struct residue ***align)
{
  int align_len, array_pos, dist, longest, max_score, score;
  int best_point, ipoint, ipos, ires, jpoint, jres, kres, npoints;
  int done;

  struct residue *residue1_ptr, *residue2_ptr;

  struct residue **alignment;

  /* Initialise variables */
  dist = 0;
  longest = nresid[0];
  if (nresid[1] > longest)
    longest = nresid[1];
  align_len = ALIGN_GAP * longest;
  npoints = nresid[0] * nresid[1];

  /* Create array that will hold the alignment */
  alignment = (struct residue **)
    malloc(align_len * 2 * sizeof(struct residue *));
  if (alignment == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for alignment array ");
      exit(1);
    }

  /* Initialise alignment array with NULLs */
  for (ipos = 0; ipos < 2 * align_len; ipos++)
    alignment[ipos] = NULL;

  /* Initialise the array position */
  array_pos = ALIGN_GAP - 1;

  /* Loop through the residues in the first sequence to place them,
     equally spaced out */
  for (ires = 0; ires < nresid[0]; ires++)
    {
      /* Get this residue */
      residue1_ptr = residue_index_ptr[0][ires];

      /* Place it in the current slot */
      alignment[array_pos] = residue1_ptr;

      /* Store this location */
      residue1_ptr->ref_pos = array_pos;

      /* Increment location */
      array_pos = array_pos + ALIGN_GAP;
    }

  /* Initialise flag */
  done = FALSE;

  /* Loop until all residues in other chain interacting with the 
     current one have been placed in order of interaction */
  while (done == FALSE)
    {
      /* Initialise the maximum interaction score */
      max_score = 0;
      best_point = -1;

      /* Loop through the scores array to find the highest-scoring
         interaction */
      for (ipoint = 0; ipoint < npoints; ipoint++)
        {
          /* If this is the highest score so far, then store */
          if (array[ipoint] > max_score)
            {
              max_score = array[ipoint];
              best_point = ipoint;
            }
        }

      /* If we have found a best point, then place the corresponding
         second residue in its optimal location */
      if (best_point > -1)
        {
          /* Find the two residues corresponding to this point */
          ires = best_point / nresid[1];
          jres = best_point % nresid[1];
          residue1_ptr = residue_index_ptr[0][ires];
          residue2_ptr = residue_index_ptr[1][jres];
          
          /* Initialise this residue's score */
          max_score = -9999;
          best_point = -1;

          /* Try the residue at various locations to see which gives
             the highest score */
          for (ipos = 0; ipos < align_len; ipos++)
            {
              /* If this position is available, then calculate a score
                 for placing the current residue there */
              if (alignment[ipos + align_len] == NULL)
                {
                  /* Initialise the score at this location */
                  score = 0;

                  /* Use this residue's interactions to calculate a score
                     for it */
                  kres = ires;
                    {
                      /* Get this point in the score array */
                      jpoint = kres * nresid[1] + jres;

                      /* If score non-zero, add it to the total, including
                         a gap penalty */
                      if (array[jpoint] > 0)
                        {
                          /* Calculate how far this point is from the 
                             residue in question */
                          dist = abs(ipos - (kres + 1) * ALIGN_GAP + 1);

                          /* Add to overall score */
                          score = score + ALIGN_GAP * array[jpoint] - dist;
                        }
                    }
                  /* If the score at this position is the best so far, then
                     store */
                  if (score > max_score)
                    {
                      /* Store details */
                      best_point = ipos;
                      max_score = score;
                    }
                }
            }

          /* Determine which was the best location */
          if (best_point > -1)
            {
              /* Place the residue at this position */
              alignment[best_point + align_len] = residue2_ptr;
              residue2_ptr->ref_pos = best_point;

              /* Loop through all this residue's scores to blank
                 them out */
              for (kres = 0; kres < nresid[0]; kres++)
                {
                  /* Get this residue's score */
                  jpoint = kres * nresid[1] + jres;
                  array[jpoint] = 0;
                }
            }
        }

      /* Otherwies, we are done */
      else
        done = TRUE;
    }

  /* Return the alignment array pointer */
  *align = alignment;

  /* Return the expanded alignment length */
  return(align_len);
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

store_alignment  -  Store the sequence alignment as an array of residue
                    pointers

***********************************************************************/

int store_alignment(int *array,int *path,int path_len,int nresid[2],
                    int align_len,
                    struct residue **residue_index_ptr[3],
                    struct residue ***align)
{
  int align_pos, array_pos, ipoint, ipos, old_align_len;
  int ires, jres, kres, last_ires, last_jres;
  int nextra, ngaps, noverlap, nshort;
  int overlap;

  struct residue *residue1_ptr, *residue2_ptr;

  struct residue **alignment;

  /* Initialise variables */
  align_pos = 0;
  last_ires = -1;
  last_jres = -1;
  ngaps = 0;
  overlap = FALSE;
  noverlap = 0;
  nshort = nresid[0];
  if (nresid[1] < nshort)
    nshort = nresid[1];
  if (nshort == 0)
    nshort = 1;

  /* Determine how many extra gaps to include in the alignment to help
     the optimization process */
  nextra = align_len / ALIGN_FRAG_LEN + 1;

  /* Increase the alignment length to allow improvement routines to have
     more 'room' to work with */
  old_align_len = align_len;
  align_len = align_len + 2 * ALIGN_MARGIN + nextra;

  /* Create array that will hold the alignment */
  alignment
    = (struct residue **)
    malloc(align_len * 2 * sizeof(struct residue *));
  if (alignment == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for alignment array ");
      exit(1);
    }

  /* Initialise alignment array with NULLs */
  for (ipos = 0; ipos < 2 * align_len; ipos++)
    alignment[ipos] = NULL;

  /* Loop over the points along the path */
  for (ipos = 0; ipos < path_len && align_pos < old_align_len; ipos++)
    {
      /* Determine how many extra gaps to be included */
      ngaps = ipos / ALIGN_FRAG_LEN;

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
              /* Retrieve the relevant residue */
              residue1_ptr = residue_index_ptr[0][kres];

              /* Store in the alignment array */
              array_pos = align_pos + ALIGN_MARGIN + ngaps;
              alignment[array_pos] = residue1_ptr;
              array_pos = align_len + align_pos + ALIGN_MARGIN + ngaps;
              alignment[array_pos] = NULL;

              /* Increment alignment position */
              align_pos++;

              /* Increment overlap length */
              if (overlap == TRUE)
                noverlap++;
            }
        }

      /* Check for gap in first sequence */
      else if (jres - last_jres > 1)
        {
          /* Write out insertion from second sequence */
          for (kres = last_jres + 1; kres < jres; kres++)
            {
              /* Retrieve the relevant residue */
              residue2_ptr = residue_index_ptr[1][kres];

              /* Store in the alignment array */
              array_pos = align_pos + ALIGN_MARGIN + ngaps;
              alignment[array_pos] = NULL;
              array_pos = align_len + align_pos + ALIGN_MARGIN + ngaps;
              alignment[array_pos] = residue2_ptr;

              /* Increment alignment position */
              align_pos++;

              /* Increment overlap length */
              if (overlap == TRUE)
                noverlap++;
            }
        }

      /* If valid, store the current alignment position */
      if (ires < nresid[0] && jres < nresid[1])
        {
          /* Retrieve the aligned residue from the first sequence */
          residue1_ptr = residue_index_ptr[0][ires];

          /* Retrieve the aligned residue from the second sequence */
          residue2_ptr = residue_index_ptr[1][jres];

          /* Store in the alignment array */
          array_pos = align_pos + ALIGN_MARGIN + ngaps;
          alignment[array_pos] = residue1_ptr;
          array_pos = align_len + align_pos + ALIGN_MARGIN + ngaps;
          alignment[array_pos] = residue2_ptr;

          /* Increment alignment position */
          align_pos++;

          /* Store pointers between the two residues */
          residue1_ptr->align_residue_ptr = residue2_ptr;
          residue2_ptr->align_residue_ptr = residue1_ptr;

          /* Increment overlap length */
          overlap = TRUE;
          noverlap++;
        }

      /* Store current position */
      last_ires = ires;
      last_jres = jres;
    }

  /* Return the alignment array pointer */
  *align = alignment;

  /* Return the expanded alignment length */
  return(align_len);
}
/***********************************************************************

fill_scores_array  -  Fill the scores array with the calculated contact
                      scores

***********************************************************************/

void fill_scores_array(struct residue **alignment,int *score_array,
                       int align_len)
{
  int array_pos, contact_score, i, ipos, j, score;

  struct link *link_ptr;
  struct residue *residue_ptr;
  struct residue *other_residue_ptr;

  /* Initialise alignment score */
  score = 0;

  /* Initialise the interact matrix */
  for (ipos = 0; ipos < align_len * align_len; ipos++)
    score_array[ipos] = 0;

  /* Loop through the alignment positions to assign reference numbers
     to the residues on both sides of the alignment */
  for (ipos = 0; ipos < align_len; ipos++)
    {
      /* Get the two residues */
      residue_ptr = alignment[ipos];
      other_residue_ptr = alignment[ipos + align_len];

      /* Save the alignment position as a reference for both */
      if (residue_ptr != NULL)
        residue_ptr->ref_pos = ipos;
      if (other_residue_ptr != NULL)
        other_residue_ptr->ref_pos = ipos;
    }

  /* Loop through all the positions in the alignment to get the
     interacting residues */
  for (ipos = 0; ipos < align_len; ipos++)
    {
      /* Get the left-hand residue pointer at this point in
         the alignment */
      residue_ptr = alignment[ipos];

      /* If have residue, get all the residues it interacts with */
      if (residue_ptr != NULL)
        {
          /* Get the first of this residue's interactions */
          link_ptr = residue_ptr->first_link_ptr;

          /* Loop through all the links to find interactions */
          while (link_ptr != NULL)
            {
              /* Get other residue */
              other_residue_ptr = link_ptr->residue_ptr[1];

              /* If not a water, then can include */
              if (strcmp(other_residue_ptr->res_name,"HOH"))
                {
                  /* Get the contact score based on the numbers of H-bonds
                     and non-bonded interactions */
                  contact_score = (HBOND_WEIGHT * link_ptr->count[HBONDS]
                                   + link_ptr->count[CONTACTS]);

                  /* Add contribution from salt-bridge/disulphide bond */
                  if (link_ptr->disulphide == TRUE)
                    contact_score = contact_score + DISULPHIDE_WEIGHT;
                  if (link_ptr->salt_bridge == TRUE)
                    contact_score = contact_score + SALT_BRIDGE_WEIGHT;
                  if (link_ptr->water_bridge == TRUE)
                    contact_score = contact_score + WATER_BRIDGE_WEIGHT;

                  /* Store the contact score at this alignment position */
                  i = residue_ptr->ref_pos;
                  j = other_residue_ptr->ref_pos;
                  if (i < 0 || j < 0)
                    {
                      printf("*** PROGRAM ERROR. Unknown array ref for");
                      if (i < 0)
                        printf(" %s %s %c (%d)",
                               residue_ptr->res_name,
                               residue_ptr->res_num,
                               residue_ptr->chain,
                               residue_ptr->ref_pos);
                      if (j < 0)
                        printf(" %s %s %c (%d)",
                               other_residue_ptr->res_name,
                               other_residue_ptr->res_num,
                               other_residue_ptr->chain,
                               other_residue_ptr->ref_pos);
                      printf("\n");
                    }
                  array_pos = i * align_len + j;
                  score_array[array_pos] = contact_score;
                }

              /* Get the next link */
              link_ptr = link_ptr->nextres_link_ptr[0];
            }
        }
    }
}
/***********************************************************************

make_random_adjustment  -  Make a random adjustment to the alignment to
                           see how it affects the overall score

***********************************************************************/

void make_random_adjustment(struct residue **alignment,int align_len,
                            int swap_pos[2][2],int *rseed,int order,
                            struct parameters params)
{
  int ipos, iseq, iseed, iswap, jpos;
  int done, swap;

  float random_no;

  struct residue *residue_ptr;

  /* Get the current random number generator seed */
  iseed = *rseed;

  /* Initialise swap positions */
  for (iseq = 0; iseq < 2; iseq++)
    for (iswap = 0; iswap < 2; iswap++)
      swap_pos[iseq][iswap] = -1;

  /* If keeping the residue order of the first sequence, only accept
     residue-swaps within the second sequence */
  if (order == TRUE)
    swap = RIGHT_SEQ;

  /* Otherwise, generate random number to determine whether to swap
     residues from one of the two interfaces, or from both */
  else
    {
      random_no = get_random_number(&iseed,RANDOM_START);
      if (random_no < LEFT_SEQ_CUTOFF)
        swap = LEFT_SEQ;
      else if (random_no < RIGHT_SEQ_CUTOFF)
        swap = RIGHT_SEQ;
      else
        swap = BOTH_SEQS;
    }

  /* Loop over the two sides of the interaction */
  for (iseq = 0; iseq < 2; iseq++)
    {
      /* If this sequence to be swapped, perform the swap */
      if (iseq == swap || swap == BOTH_SEQS)
        {
          /* Initialise flag */
          done = FALSE;

          /* Generate two random numbers to determine which positions
             to swap */
          while (done == FALSE)
            {
              /* Get first position */
              random_no = get_random_number(&iseed,RANDOM_START);
              ipos = (int) (align_len * random_no);
              if (ipos < 0)
                ipos = 0;
              else if (ipos > align_len - 1)
                ipos = align_len - 1;

              /* Get second position */
              random_no = get_random_number(&iseed,RANDOM_START);
              jpos = (int) (align_len * random_no);
              if (jpos < 0)
                jpos = 0;
              else if (jpos > align_len - 1)
                jpos = align_len - 1;

              /* If this is not the same position, and both positions aren't
                 gaps, we can use these selected positions */
              if (ipos != jpos &&
                  (alignment[ipos + iseq * align_len] != NULL ||
                   alignment[jpos + iseq * align_len] != NULL))
                done = TRUE;
            }

          /* Store the swapped positions */
          swap_pos[iseq][0] = ipos;
          swap_pos[iseq][1] = jpos;

          /* Get the two residue pointers at these points in the
             alignment and swap */
          residue_ptr = alignment[ipos + iseq * align_len];
          alignment[ipos + iseq * align_len]
            = alignment[jpos + iseq * align_len];
          alignment[jpos + iseq * align_len] = residue_ptr;
        }
    }

  /* Return current random number generator seed */
  *rseed = iseed;
}
/***********************************************************************

calc_alignment_score  -  Calculate given alignment based on numbers of
                         H-bonds and non-bonded contacts between aligned
                         residues and lengths of other interactions

***********************************************************************/

int calc_alignment_score(struct residue **alignment,int *score_array,
                         int *tmp,int align_len)
{
  int align_score, array_pos, contact_score, dist, i, ipos, j, jpos;
  int ii, jj, ncross;

  struct residue *residue_ptr;
  struct residue *other_residue_ptr;

  /* Initialise alignment score */
  align_score = 0;

  /* Initialise the temporary scores matrix */
  for (ipos = 0; ipos < align_len * align_len; ipos++)
    tmp[ipos] = 0;

  /* Loop through all the positions in the alignment to transfer
     appropriate scores to the temporary array */
  for (ipos = 0; ipos < align_len; ipos++)
    {
      /* Get the corresponding residue pointer */
      residue_ptr = alignment[ipos];

      /* Loop over all other residues */
       for (jpos = 0; jpos < align_len && residue_ptr != NULL; jpos++)
         {
           /* Get the two residues */
           other_residue_ptr = alignment[jpos + align_len];

           /* If both valid, get their interaction score */
           if (other_residue_ptr != NULL)
             {
               /* Get their positions in the scores array */
               i = residue_ptr->ref_pos;
               j = other_residue_ptr->ref_pos;

               /* Get the array position and contact score */
               array_pos = i * align_len + j;
               contact_score = score_array[array_pos];

               /* Get the relevant location in the temporary array */
               array_pos = ipos * align_len + jpos;
               tmp[array_pos] = contact_score;
             }
         }
    }

  /* Loop through the temporary array to calculate the total score */
  for (ipos = 0; ipos < align_len; ipos++)
    {
      /* Loop over all other residues */
       for (jpos = 0; jpos < align_len; jpos++)
         {
           /* Get the array position and contact score */
           array_pos = ipos * align_len + jpos;
           contact_score = tmp[array_pos];

           /* If not zero, proceed */
           if (contact_score > 0)
             {
               /* If residues not adjacent in the alignment, apply
                  adjustment */
               if (ipos == jpos)
                 align_score = align_score - contact_score;

              /* Otherwise, need to adjust for alignment separation */
              else
                {
                  /* Get the separation in the alignment */
                  dist = abs(ipos - jpos);

                  /* Add numbers of iteractions to overall score */
                  align_score = align_score + dist * contact_score;
                }

               /* Calculate number of other interactions that intersect
                  with this one and add penalty score */
               ncross = 0;

               /* Loop over all the other residues on the left-hand
                  side of the alignment to see how many have interactions
                  that cross the current one */
               for (ii = 0; ii < ipos; ii++)
                 {
                   /* Loop over all residues on the right-hand side
                      of the alignment */
                   for (jj = jpos + 1; jj < align_len; jj++)
                     {
                       /* Get the array position and contact score */
                       array_pos = ii * align_len + jj;
                       contact_score = tmp[array_pos];

                       /* If contact score non-zero then have a cross
                          over */
                       if (contact_score > 0)
                         ncross++;
                     }
                 }

               /* Repeat the other way round */
               for (ii = ipos + 1; ii < align_len; ii++)
                 {
                   /* Loop over all residues on the right-hand side
                      of the alignment */
                   for (jj = 0; jj < jpos; jj++)
                     {
                       /* Get the array position and contact score */
                       array_pos = ii * align_len + jj;
                       contact_score = tmp[array_pos];

                       /* If contact score non-zero then have a cross
                          over */
                       if (contact_score > 0)
                         ncross++;
                     }
                 }

               /* Add penalty for number of cross-overs we have */
               align_score = align_score + ncross * CROSS_PENALTY;
             }
         }
    }

  /* Return the alignment score */
  return(align_score);
}
/***********************************************************************

judge_move  -  Determine whether move should be accepted or rejected

***********************************************************************/

int judge_move(int current_score,int score, double temperature,
               int *rseed,struct parameters params)
{
  int iseed;
  int accept_move;

  float random_no;

  double energy_diff, prob;

  /* Initialise variables */
  accept_move = FALSE;

  /* Get the current random number generator seed */
  iseed = *rseed;

  /* Calculate change in energy of making the given move */
  energy_diff = score - current_score;

  /* If score has decreased, accept it */
  if (energy_diff < 0)
    accept_move = TRUE;

  /* Otherwise, calcuate the Boltzmann factor */
  else
    {
      /* Calculate probability of energy change */
      prob = exp(- energy_diff / temperature);

      /* Generate a random number */
      random_no = get_random_number(&iseed,RANDOM_START);

      /* If random number lower than probability value, then accept
         move */
      if (random_no < prob)
        accept_move = TRUE;
    }

  /* Return current random number generator seed */
  *rseed = iseed;

  /* Return whether move has been accepted */
  return(accept_move);
}
/***********************************************************************

undo_adjustment  -  Undo the adjustment just made to the alignment

***********************************************************************/

void undo_adjustment(struct residue **alignment,int align_len,
                     int swap_pos[2][2])
{
  int ipos, iseq, jpos;

  struct residue *residue_ptr;

  /* Loop over the two sides of the interaction */
  for (iseq = 0; iseq < 2; iseq++)
    {
      /* If this sequence to be swapped, perform the swap */
      if (swap_pos[iseq][0] > -1 && swap_pos[iseq][1] > -1)
        {
          /* Get swapped positions */
          ipos = swap_pos[iseq][0];
          jpos = swap_pos[iseq][1];

          /* Get the two residue pointers at these points in the
             alignment and swap */
          residue_ptr = alignment[ipos + iseq * align_len];
          alignment[ipos + iseq * align_len]
            = alignment[jpos + iseq * align_len];
          alignment[jpos + iseq * align_len] = residue_ptr;
        }
    }
}
/***********************************************************************

improve_alignment  -  Use simulated annealing to improve the alignment
                      to get the interacting residues as close to one
                      as possible

***********************************************************************/

void improve_alignment(struct residue **alignment,int align_len,
                       int *rseed,int order,struct parameters params)
{
  int current_score, score, starting_score;
  int n, nmoves, nsteps, nsuccesses, swap_pos[2][2];
  int max_steps, moves_at_const_temp, successes_at_const_temp;
  int accept_move, done;

  int *score_array, *tmp;

/* v.5.4.4 --> */
  float factor = 1.0;
/* <-- v.5.4.4 */

  double temperature, temperature_reduction_factor;

  /* Initialise variables */
  accept_move = FALSE;
  done = FALSE;
  max_steps = MAX_STEPS;
  nsteps = 0;
  n = align_len;
/* debug */
  printf("align_len = %d (max = %d)\n",align_len,MAX_ALIGN_LEN);
  fflush(stdout);
/* debug */

  /* If alignment too long, then minimization just takes forever, so
     peg back to a "maximum alignment length" to shorten the process */
  if (align_len > MAX_ALIGN_LEN)
    {
      n = MAX_ALIGN_LEN;
/* v.5.4.4 --> */
//      max_steps = 10;
      factor = (float) MAX_ALIGN_LEN / (float) align_len;
      max_steps = MIN_STEPS + factor * (MAX_STEPS - MIN_STEPS);
/* <-- v.5.4.4 */
    }
/* debug */
  printf("In improve_alignment with max_steps = %d\n",max_steps);
  fflush(stdout);
  printf("  Factor = %9.6f\n",factor);
  fflush(stdout);
/* debug */

  /* Annealing schedule parameters */
/* v.5.4.4 --> */
//  moves_at_const_temp = 100 * n;
//  successes_at_const_temp = 10 * n;
  moves_at_const_temp = factor * 100 * n;
  successes_at_const_temp = factor * 10 * n;
/* <-- v.5.4.4 */
  temperature_reduction_factor = 0.9;
  temperature = STARTING_TEMP;

  /* Create matrix to store all the residue-residue interaction scores */
  score_array = (int *) malloc(align_len * align_len * sizeof(int));
  if (score_array == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for score_array\n");
      exit(1);
    }

  /* Fill the scores array with the interaction scores between all
     residue pairs */
  fill_scores_array(alignment,score_array,align_len);

  /* Create matrix to store all the residue-residue interaction scores */
  tmp = (int *) malloc(align_len * align_len * sizeof(int));
  if (tmp == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for tmp array\n");
      exit(1);
    }

  /* Calculate initial alignment score */
  current_score
    = calc_alignment_score(alignment,score_array,tmp,align_len);
  starting_score = current_score;
/* debug */
  printf("Starting score = %d\n",starting_score);
  printf("moves_at_const_temp = %d, successes_at_const_temp = %d\n",
         moves_at_const_temp,successes_at_const_temp);
  fflush(stdout);
/* debug */

  /* Loop until reasonable score attained */
  while (done == FALSE && nsteps < max_steps)
    {
      /* Initialise counts for moves made at this temperature */
      nmoves = 0;
      nsuccesses = 0;

      /* Loop for required number of moves at this temperature */
      while (done == FALSE && nmoves < moves_at_const_temp &&
             nsuccesses < successes_at_const_temp)
        {
          /* Make a random change to the alignment */
          make_random_adjustment(alignment,align_len,swap_pos,rseed,
                                 order,params);

          /* Score alignment to assigne an energy value to it */
          score
            = calc_alignment_score(alignment,score_array,tmp,align_len);

          /* Determine whether move should be accepted or rejected */
          accept_move
            = judge_move(current_score,score,temperature,rseed,params);

          /* If move accepted, store current energy and increment 
             count of successes */
          if (accept_move == TRUE)
            {
              /* If energy changed, increment count of successes */
              if (current_score != score)
                nsuccesses++;

              /* Store current score */
              current_score = score;
            }

          /* If move rejected, undo move */
          else
            undo_adjustment(alignment,align_len,swap_pos);

          /* Increment count of moves made at this temperature */
          nmoves++;
        }

      /* If there were no successful moves accepted at the current
         temperature, then we're done */
      if (nsuccesses == 0)
        done = TRUE;

      /* Adjust temperature */
      temperature = temperature * temperature_reduction_factor;

      /* Increment count of temperature steps taken */
      nsteps++;

      /* Show the number of step taken and the current score */
      printf("Loop: %4d. Current score %d\n",nsteps,current_score);
      fflush(stdout);
    }

  /* Calculate final alignment score */
  score = calc_alignment_score(alignment,score_array,tmp,align_len);

  /* Show difference between start and end scores */
  printf("Alignment score at start = %d and at end = %d\n",
         starting_score,score);
}
/***********************************************************************

close_up_gaps  -  Close up any gaps in the alignment that might have
                  opened as a result of the minimization

***********************************************************************/

int close_up_gaps(struct residue **alignment,int align_len)
{
  int ipos, next_pos;

  struct residue *residue1_ptr, *residue2_ptr;

  /* Initialise variables */
  next_pos = 0;

  /* Loop through all the positions in the alignment to fill any blanks */
  for (ipos = 0; ipos < align_len; ipos++)
    {
      /* Get the corresponding residue pointers at this point in
         the alignment */
      residue1_ptr = alignment[ipos];
      residue2_ptr = alignment[ipos + align_len];

      /* If either residue pointers not NULL then have a valid slot
         in the alignment */
      if (residue1_ptr != NULL || residue2_ptr != NULL)
        {
          /* If have a free slot before the current one, move both
             pointers to close gap */
          if (next_pos < ipos)
            {
              /* Shift the pointers */
              alignment[next_pos] = residue1_ptr;
              alignment[next_pos + align_len] = residue2_ptr;
              alignment[ipos] = alignment[ipos + align_len] = NULL;
            }

          /* Increment location of next free position */
          next_pos++;
        }
    }

  /* Return the new alignment length */
  return(next_pos);
}
/***********************************************************************

close_up_more  -  Close up gaps by blocks, if can do so

***********************************************************************/

int close_up_more(struct residue **alignment,int align_len)
{
  int icol, ipos, jcol, kpos, loop, next_pos;
  int ires, jres, last_jres, next_ires, next_jres, nother, nup;
  int done;

  /* Initialise variables */
  next_pos = 0;

  /* Loop twice, first down the left column and then down the right */
  for (loop = 0; loop < 2; loop++)
    {
      /* Determine which column we're going down */
      icol = loop * align_len;
      jcol = (1 - loop) * align_len;

      /* Loop through all the positions in this column to find a terminal
         residue */
      for (ipos = 0; ipos < align_len - 1; ipos++)
        {
          /* Get this position */
          ires = ipos + icol;
          jres = ipos + jcol;
          next_ires = next_jres = -1;

          /* If there is a blank below it, and next to it, then we might
             be able to close up */
          if (alignment[ires] != NULL &&
              (alignment[ires + 1] == NULL && alignment[jres] == NULL))
            {

              /* Find the first non-blank position below it */
              for (kpos = ires + 1; kpos < align_len + icol &&
                     next_ires == -1; kpos++)
                {
                  /* If residue not blank, then store position */
                  if (alignment[kpos] != NULL)
                    next_ires = kpos;
                }

              /* Do the same for the other column */
              for (kpos = jres + 1; kpos < align_len + jcol &&
                     next_jres == -1; kpos++)
                {
                  /* If residue not blank, then store position */
                  if (alignment[kpos] != NULL)
                    next_jres = kpos;
                }

              /* If have found starting residues that can be moved,
                 perform the move */
              if (next_ires > -1 && next_jres > -1)
                {
                  /* Work out the number of positions we can shift current
                     columns */
                  nup = next_ires - ires - 1;
                  last_jres = jres;
                  done = FALSE;
                  
                  /* Work out how much the other column can be shifted by */
                  for (kpos = jres; kpos > -1 && done == FALSE; kpos--)
                    {
                      /* If residue not blank, then store position */
                      if (alignment[kpos] == NULL)
                        last_jres = kpos;
                      else
                        done = TRUE;
                    }

                  /* Work out how far this column can be moved up */
                  nother = next_jres - last_jres;


                  /* If other shift is smaller, have to use that */
                  if (nother < nup)
                    nup = nother;

                  /* Perform the shift of the current column */
                  for (kpos = next_ires; kpos < align_len + icol; kpos++)
                    {
                      /* Perform the shift */
                      alignment[kpos - nup] = alignment[kpos];
                      alignment[kpos] = NULL;
                    }

                  /* Do the same for the other column */
                  for (kpos = next_jres; kpos < align_len + jcol; kpos++)
                    {
                      /* Perform the shift */
                      alignment[kpos - nup] = alignment[kpos];
                      alignment[kpos] = NULL;
                    }
                }
            }
        }
    }

  /* Determine the new length of the alignment */
  for (ipos = 0; ipos < align_len; ipos++)
    {
      /* If either left- or right-hand residue exists, then this is
         the last position in the new alignment */
      if (alignment[ipos] != NULL || alignment[ipos + align_len])
        next_pos = ipos + 1;
    }

  /* Return the new alignment length */
  return(next_pos);
}
/***********************************************************************

place_residues  -  Pair up the residues on either side of the interface
                   to maximize pairing of most interacting residues

***********************************************************************/

int place_residues(struct residue *first_residue_ptr,int *nresid,
                   struct residue ***align,int *aln_len,int *max_count,
                   int *rseed,struct psparams *psparams_ptr,
                   struct pspage *pspage_ptr,struct parameters params)
{
  int align_len, nbridges, path_len, total_len;

  int *array, *path;

  struct residue **alignment, **residue_index_ptr[3];

  /* Initialise variables */
  nbridges = 0;

  /* Create the arrays that will be used in computing the alignment
     of the interface residues */
  create_alignment_arrays(nresid,&array,&path,first_residue_ptr,
                          residue_index_ptr);

  /* Initialise alignment arrays */
  init_arrays(nresid,residue_index_ptr,array,max_count);

  /* If first sequence to retain its order, align the residues of
     the second sequence to maximuze the neighbouring interaction
     scores */
  if (params.order == TRUE)
    align_len
      = store_ordered_alignment(array,nresid,residue_index_ptr,&alignment);

  /* Otherwise, perform a sequence alignment of the two sequences,
     weighted by the interactions between them */
  else
    {
      /* Apply Needleman and Wunsch algorithm to adjust scores matrix */
      nw_sweep(array,nresid);

      /* Find best route through matrix to give best alignment of the two
         chains */
      path_len = find_best_path(array,nresid,path,&align_len);

      /* Store the alignment as an array of residue pointers */
      align_len = store_alignment(array,path,path_len,nresid,align_len,
                                  residue_index_ptr,&alignment);

      /* Improve alignment by minimising its "energy" */
      if (align_len > 1)
        improve_alignment(alignment,align_len,rseed,params.order,params);
    }

  /* Close up any gaps in the alignment that might have opened as a
     result of the minimization*/
  total_len = close_up_gaps(alignment,align_len);

  /* Close them up even more */
  total_len = close_up_more(alignment,align_len);

  /* Return the final alignment */
  *align = alignment;
  *aln_len = align_len;

  /* Return the total length of the alignment */
  return(total_len);
}
/***********************************************************************

write_save_file  -  Write out the plot data to the dimhtmlXY.dat file

***********************************************************************/

void write_save_file(struct residue **alignment,int align_len,
                     int total_len,char *pdb_code,char *chains,
                     struct residue *first_residue_ptr,
                     struct link *first_link_ptr,
                     struct parameters params)
{
  char filename[FILENAME_LEN];

  int ipos;

  struct link *link_ptr;
  struct residue *residue_ptr, *residue1_ptr, *residue2_ptr;

  FILE *file_ptr;

  /* Form the name of the save file, dimhtmlXY.dat */
  if (params.work_dir[0] != '\0')
    sprintf(filename,"%s%sdimhtml%c%c.out",params.work_dir,DIR_SEP,
            chains[0],chains[1]);
  else
    sprintf(filename,"dimhtml%c%c.out",chains[0],chains[1]);

  /* Open the output file */
  file_ptr = open_output_file(filename,TRUE);

  /* Get pointer to first residue record */
  residue_ptr = first_residue_ptr;

  /* Loop through residue records to write out */
  while (residue_ptr != NULL)
    {
      /* Write out the residue name */
      fprintf(file_ptr,"%s %s %c %d\n",residue_ptr->res_name,
              residue_ptr->res_num,residue_ptr->chain,residue_ptr->number);

      /* Get next residue record */
      residue_ptr = residue_ptr->next_residue_ptr;
    }

  /* Write heading to signify alignment positions */
  fprintf(file_ptr,"INTERFACE %d\n",total_len);

  /* Loop through all the positions in the alignment to write out */
  for (ipos = 0; ipos < total_len; ipos++)
    {
      /* Get the first residue in the alignment and write out */
      residue_ptr = alignment[ipos];
      if (residue_ptr != NULL)
        fprintf(file_ptr,"%d ",residue_ptr->number);
      else
        fprintf(file_ptr,"-1 ");

      /* Repeat for second residue */
      residue_ptr = alignment[ipos + align_len];
      if (residue_ptr != NULL)
        fprintf(file_ptr,"%d\n",residue_ptr->number);
      else
        fprintf(file_ptr,"-1\n");
    }

  /* Write heading to signify interaction details to follow */
  fprintf(file_ptr,"INTERACTIONS\n");

  /* Get the first interaction link */
  link_ptr = first_link_ptr;

  /* Loop through all the interactions to write out */
  while (link_ptr != NULL)
    {
      /* Get the two residues */
      residue1_ptr = link_ptr->residue_ptr[0];
      residue2_ptr = link_ptr->residue_ptr[1];

      /* Write out the interactions between these two residues */
      fprintf(file_ptr,"%d %d %d %d",residue1_ptr->number,
              residue2_ptr->number,link_ptr->count[HBONDS],
              link_ptr->count[CONTACTS]);
      if (link_ptr->salt_bridge == TRUE)
        fprintf(file_ptr," s");
      if (link_ptr->disulphide == TRUE)
        fprintf(file_ptr," d");
      fprintf(file_ptr,"\n");

      /* Get the next link */
      link_ptr = link_ptr->next_link_ptr;
    }

  /* Close the output file */
  fclose(file_ptr);
}
/***********************************************************************

old_write_dimplot_rcm  -  Write out the residue centres of mass to the
                      dimplot.rcm file

***********************************************************************/

/* v.5.3.4--> */
//void write_dimplot_rcm(struct residue *first_residue_ptr,float ygap,
//                       int nbridges,struct parameters params)
void old_write_dimplot_rcm(struct residue *first_residue_ptr,float ygap,
                           int nbridges,struct parameters params)
/* <--v.5.3.4 */
{
  char file_name[FILENAME_LEN];

  char chain;

  float scale_factor_x, scale_factor_y, x, x_sep, y;

  struct residue *residue_ptr;

  FILE *file_ptr;

  /* Calculate the x and y scale factors to write out the coordinates
     in (approximately) Angstroms */
  if (ygap > 0)
    scale_factor_y = RCM_SEP_Y / ygap;
  else
    scale_factor_y = 1.0;
  x_sep = scale_factor_x = RCM_SEP_X;
  if (nbridges > 0)
    x_sep = scale_factor_x = scale_factor_x + RCM_WATER_SEP_X;
  if (params.flip_seqs == TRUE)
    scale_factor_x = - scale_factor_x;

  /* Form name of the .rcm file */
  file_name[0] = '\0';
  if (params.work_dir[0] != '\0')
    {
      strcpy(file_name,params.work_dir);
      strcat(file_name,DIR_SEP);
    }
  strcat(file_name,"dimplot.rcm");

  /* Open the output file */
  file_ptr = open_output_file(file_name,TRUE);

  /* Write out the header records */
  fprintf(file_ptr,"                 Res.  Res       Flattened    \n");
  fprintf(file_ptr,"            Atom Name  Num      coordinates   \n");
  fprintf(file_ptr,"            ---- --- -----    ----------------\n");

  /* Set the pointer to the first residue */
  residue_ptr = first_residue_ptr;

  /* Loop over all this object's residues */
  while (residue_ptr != NULL)
    {
      /* Write if not water */
      if (residue_ptr->chain_no < 2)
        {
          /* If chain is blank, write out as a hyphen */
          chain = residue_ptr->chain;
          if (chain == ' ')
            chain = '-';

          /* Determine the coordinates of this residue */
          x = scale_factor_x * (residue_ptr->chain_no - x_sep / 2);
          y = scale_factor_y * residue_ptr->ycoord;

          /* Write out the residue's centre-of-mass coords */
          fprintf(file_ptr,"RESDUE      CofM %s %c%s   %8.3f%8.3f\n",
                  residue_ptr->res_name,chain,residue_ptr->res_num,x,y);
        }

      /* Get pointer to the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
    }

  /* Close the output file */
  fclose(file_ptr);
}
/* v.5.3.4--> */
/***********************************************************************

write_rcm_details  -  Write out the details for the last line of
                      residues

***********************************************************************/

float write_rcm_details(FILE *file_ptr,struct residue **store_residue_ptr,
/* v.5.3.6--> */
//                        int npos,float y,int sign)
                        int npos,float y,int sign,int signy,
                        struct parameters params)
/* <--v.5.3.6 */
{
  char chain;

  int ipos;

  float size, x;

  struct residue *residue_ptr;

  /* Initialise variables */
/* v.5.3.6--> */
//  size = SIZE_HYDROPHOBIC;
  size = params.size_hydrophobic;
/* <--v.5.3.6 */

  /* Determine largest residue (ie hydrophobic or H-bonded) in this
     line */
  for (ipos = 0; ipos < npos; ipos++)
    if (store_residue_ptr[ipos]->nhbonds > 0)
      size = SIZE_HBOND;

  /* Calculate the y-coordinate */
/* v.5.3.6--> */
//  y = y + size + GAP_Y;
  y = y + signy * (size + GAP_Y);
/* <--v.5.3.6 */

  /* Loop over the residues to write out their CofM position */
  for (ipos = 0; ipos < npos; ipos++)
    {
      /* Get this residue */
      residue_ptr = store_residue_ptr[ipos];

      /* Get this residue's size */
      if (residue_ptr->nhbonds > 0)
        size = SIZE_HBOND / 2;
      else
/* v.5.3.6--> */
//        size = SIZE_HYDROPHOBIC / 2;
        size = params.size_hydrophobic / 2;
/* <--v.5.3.6 */

      /* Calculate this residue's x-coord */
      if (residue_ptr->chain_no < 2)
        x = sign * (2 * residue_ptr->chain_no - 1)
          * (size + IDEAL_GAP / 2);
      else
        x = 0;

      /* If chain is blank, write out as a hyphen */
      chain = residue_ptr->chain;
      if (chain == ' ')
        chain = '-';

      /* Write out the residue's centre-of-mass coords */
      fprintf(file_ptr,"RESDUE      CofM %s %c%s   %8.3f%8.3f\n",
              residue_ptr->res_name,chain,residue_ptr->res_num,x,y);
    }

  /* Return the current y-value */
  return(y);
}
/***********************************************************************

write_dimplot_rcm  -  Write out the residue centres of mass to the
                      dimplot.rcm file

***********************************************************************/

void write_dimplot_rcm(char *fname,struct residue *first_residue_ptr,
                       int nresidues,struct parameters params)
{
  char file_name[FILENAME_LEN];
  char char_y[20], last_char_y[20];

/* v.5.3.6--> */
//  int ipos, npos, sign;
  int ipos, npos, sign, signy;
/* <--v.5.3.6 */

  float y;

  struct residue *residue_ptr, *store_residue_ptr[2];

  FILE *file_ptr;

  /* Initialise variables */
  char_y[0] = last_char_y[0] = '\0';
  for (ipos = 0; ipos < 2; ipos++)
    store_residue_ptr[ipos] = NULL;
  npos = 0;
  sign = 1;
  y = - (SIZE_HBOND + GAP_Y) * nresidues / 2;
  if (params.flip_seqs == TRUE)
    sign = -1;
/* v.5.3.6--> */
  signy = 1;
  if (params.portrait == TRUE)
    signy = -1;
/* <--v.5.3.6 */

  /* Form name of the .rcm file */
  file_name[0] = '\0';
  if (params.work_dir[0] != '\0')
    {
      strcpy(file_name,params.work_dir);
      strcat(file_name,DIR_SEP);
    }
  strcat(file_name,fname);

  /* Open the output file */
  file_ptr = open_output_file(file_name,TRUE);

  /* Write out the header records */
  fprintf(file_ptr,"                 Res.  Res       Flattened    \n");
  fprintf(file_ptr,"            Atom Name  Num      coordinates   \n");
  fprintf(file_ptr,"            ---- --- -----    ----------------\n");

  /* Set the pointer to the first residue */
  residue_ptr = first_residue_ptr;

  /* Loop over all this object's residues */
  while (residue_ptr != NULL)
    {
/* v.5.3.8--> */
      /* Only process if the residue hasn't been deleted */
      if (residue_ptr->deleted == FALSE)
        {
/* <--v.5.3.8 */
          /* Format the y-coord of this residue */
          sprintf(char_y,"%8.1f",residue_ptr->ycoord);

          /* If it has changed, then process last line of residues,
             if there is one */
          if (strcmp(char_y,last_char_y))
            {
              /* If have residues, then process them */
              if (npos > 0)
                y = write_rcm_details(file_ptr,store_residue_ptr,npos,
/* v.5.3.6--> */
//                                  y,sign);
                                      y,sign,signy,params);
/* <--v.5.3.6 */

              /* Reinitialise */
              for (ipos = 0; ipos < 2; ipos++)
                store_residue_ptr[ipos] = NULL;
              npos = 0;

              /* Save the y-coord */
              strcpy(last_char_y,char_y);
            }

          /* Store residue in next available slot */
          if (npos < 2)
            {
              /* If an interface residue, then store */
              if (residue_ptr->chain_no < 2)
                {
                  /* Store this residue pointer */
                  store_residue_ptr[npos] = residue_ptr;

                  /* Increment slot position */
                  npos++;
                }
            }

          /* Otherwise, show warning message */
          else
            {
              printf("*** Warning. No slot for residue %s %s %c\n",
                     residue_ptr->res_name,
                     residue_ptr->res_num,
                     residue_ptr->chain);
            }
/* v.5.3.8--> */
        }
/* <--v.5.3.8 */

      /* Get pointer to the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
    }

  /* If have last line of residues, then process them */
  if (npos > 0)
/* v.5.3.6--> */
//    write_rcm_details(file_ptr,store_residue_ptr,npos,y,sign);
    write_rcm_details(file_ptr,store_residue_ptr,npos,y,sign,signy,
                      params);
/* <--v.5.3.6 */

  /* Close the output file */
  fclose(file_ptr);
}
/* <--v.5.3.4 */

/***********************************************************************

                             M   A   I   N

***********************************************************************/

int main(int argc,char *argv[])
{
  char pdbsum_dir[FILENAME_LEN + 1];

  int align_len, nbridges, nchains, nresid[3], nresidues, total_len;
  int dir_type, max_count[3];
  int iseed;

  float ygap;

  struct c_file_data file_name_ptr;

  struct parameters params;
  struct psparams psparams_ptr;
  struct pspage pspage_ptr;

  struct residue *first_residue_ptr;
  struct link *first_link_ptr, *last_link_ptr;

  struct residue **alignment;

  /* Initialise variables */
  first_link_ptr = NULL;
  last_link_ptr = NULL;
  nchains = 0;
  max_count[0] = 0;
  max_count[1] = 0;

  /* Read in the command-line arguments */
  get_command_arguments(argv,argc - 1,&params);

  /* Get required paths and file names */
/* v.5.3.4--> */
  if (strcmp(params.pdb_code,"none"))
    {
/* <--v.5.3.4 */
      c_get_paths(&file_name_ptr,params.run_64);

      /* Form the name of the PDBsum directory */
      dir_type = find_pdbsum_dir(params.pdb_code,
                                 file_name_ptr.pdbsum_template,
                                 params.pdb_type,pdbsum_dir);
/* v.5.3.4--> */
    }
/* <--v.5.3.4 */

  /* Initialise PostScript parameters */
  initialise_psparams(&psparams_ptr,&pspage_ptr);

  /* Generate plot from LIGPLOT data */
  if (params.use_data == FALSE)
    {
      /* Initialise random number generator */
      iseed = -1;
      get_random_number(&iseed,RANDOM_START);

      /* If reading in files created by DIMPLOT, read in the ligplot.rcm
         file and the ligplot.hhb and .nnb files */
      if (params.dimplot == TRUE)
        {
          /* Unless running as part of the DIMPLOT script, scan the
             ligplot.rcm file to pick up the two chains involved */
          if (strcmp(params.pdb_code,"none"))
            {
              nchains = scan_rcm_file(params.chain,params);

              /* Read in the list of interacting residues from the
                 ligplot.rcm file */
              nresidues = read_rcm_file(&first_residue_ptr,nresid,params);
            }

          /* If running as part of the DIMPLOT script, read through the
             dimplot.pdb file to assign residues to appropriate interface */
          else
            nresidues
              = read_dimplot_pdb(&first_residue_ptr,nresid,&params);

          /* Sort the residues by residue number of dimplot y-coord */
          sort_residues(&first_residue_ptr,nresidues,params.order);

          /* Read in the H-bonds from the ligplot.hhb file */
          read_hhb_file(first_residue_ptr,&first_link_ptr,&last_link_ptr,
                    HBONDS,params);

          /* Read in the non-bonded contacts from the ligplot.nnb file */
          read_hhb_file(first_residue_ptr,&first_link_ptr,&last_link_ptr,
                        CONTACTS,params);

          /* Add links representing water bridges */
          nbridges = add_water_bridges(first_residue_ptr,&first_link_ptr,
                                       &last_link_ptr);

          /* Delete any residues not involved in interactions across the
             interface */
          nresidues
            = delete_residues(&first_residue_ptr,nresid,&nchains);

          /* If plotting domains, and we have more than one chain, then
             set flag to include chain-id in final plot */
          if (params.domain[0] != 0 && nchains > 2)
            params.show_chains = TRUE;
        }

      /* Otherwise, read in the dimhtmlXY.dat file created by newsec */
      else
        {
          /* Read in the list of residues and their interactions */
          total_len
            = read_save_file(pdbsum_dir,&alignment,&align_len,
                             params.pdb_code,params.chain,max_count,
                             &first_residue_ptr,&nresidues,nresid,
/* v.5.4.5 --> */
//                             &first_link_ptr,"dat",params.work_dir);
                             &first_link_ptr,"dat",params.work_dir,
                             params.pdb_type);
/* <-- v.5.4.5 */

          /* Sort the residues by their y-coord */
          sort_residues(&first_residue_ptr,nresidues,FALSE);
        }

      /* Place the interface residues down either side, lining up interacting
         residues first */
      total_len
        = place_residues(first_residue_ptr,nresid,&alignment,&align_len,
                         max_count,&iseed,&psparams_ptr,
                         &pspage_ptr,params);
    }

  /* If prog has been run before, read in the plot data from the
     dimhtmlXY.dat file */
  else
    total_len
      = read_save_file(pdbsum_dir,&alignment,&align_len,
                       params.pdb_code,params.chain,max_count,
                       &first_residue_ptr,&nresidues,nresid,
/* v.5.4.5 --> */
//                       &first_link_ptr,"out",params.work_dir);
                       &first_link_ptr,"out",params.work_dir,
                       params.pdb_type);
/* <-- v.5.4.5 */

  /* If colouring residues by residue conservation, pick up the ConSurf
     conservation scores */
  if (params.colour_scheme == CONSERVATION)
    {
      /* Read in the conservation data for the first chain */
      get_consurf_scores(file_name_ptr,params.pdb_code,params.chain[0],
                         first_residue_ptr);

      /* Repeat for second chain */
      get_consurf_scores(file_name_ptr,params.pdb_code,params.chain[1],
                         first_residue_ptr);
    }

  /* Plot the alignment, if we have one */
  if (total_len > 0)
    ygap
      = plot_diagram(alignment,align_len,total_len,max_count,psparams_ptr,
                     &pspage_ptr,pdbsum_dir,params);

  /* Write out data file holding details of the interface */
  if (params.save_data == TRUE)
    write_save_file(alignment,align_len,total_len,params.pdb_code,
                    params.chain,first_residue_ptr,first_link_ptr,
                    params);

  /* Sort the residues by their y-coord */
  sort_residues(&first_residue_ptr,nresidues,FALSE);

  /* Write out the dimplot.rcm file for use with DIMPLOT */
/* v.5.3.4--> */
//  write_dimplot_rcm(first_residue_ptr,ygap,nbridges,params);
  write_dimplot_rcm("dimplot.rcm",first_residue_ptr,nresidues,params);

  /* Repeat to write out as "dimhtml.rcm" */
  write_dimplot_rcm("dimhtml.rcm",first_residue_ptr,nresidues,params);
/* <--v.5.3.4 */

  return(0);
}
