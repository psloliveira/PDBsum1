/***********************************************************************

  rascript.c - Program to generate a RasMol script, and appended PDB file,
               in PDBsum. Options allow the protein surface and binding
               site info to be included. New features include generation
               of MolScript script

************************************************************************

Date:-         17 Feb 1998
Dir update:-    2 Oct 2012

Written by:-   Roman Laskowski

----------------------------------------------------------------------*/

/* Datafiles
   ---------

Actual name      File name        Description
-----------      ---------        -----------
CATHPARAM        file_name        Input parameter file giving all the
                                  directory paths for the files used by
                                  the program
rasmol.dfn       dfn_name         Input file holding the RasMol script
                                  definitions from which the script
                                  commands are to be generated
<filename>.pdb   pdb_name         Input PDB file holding coords of
                                  protein to be processed
cath.domains     cath_domains_file CATH domains file giving residue-ranges
                                  for each domain in each PDB file
mask.vtx         vtx_name         Input file holding the vertex coords
                                  for generating the binding-site mask
                                  (as computed by the SURFNET and SURFACE
                                  programs)
rascript.out     out_name         Output file (created when -f option
                                  given - otherwise output goes to stdout)
sst.dat          sst_dat          Input file giving dominant secondary
                                  structure for each chain and governing
                                  how chains are to be coloured

*/

/* Function calling-tree
   ---------------------

main
   -> get_command_arguments
   -> c_get_paths
         -> get_param (in reapar.c)
         -> store_string (in common.c)
   -> form_filenames
         -> r_find_pdb_file (in readpdb.c)
         -> find_pdbsum_dir (in common.c)
   -> open_output_file (in common.c)
   -> get_dfn_data
         -> string_truncate (in common.c)
         -> create_rasdef_entry
         -> get_ligand_residues
               -> get_token
               -> format_rasdef_number
               -> create_rasdef_entry
         -> store_chain
   -> get_sst_dat
         -> string_truncate (in common.c)
   -> get_promotif_data
         -> read_promotif_dat
               -> create_promotif_record
         -> read_in_promotif_data
               -> create_promotif_record
   -> read_in_grow_data
         -> create_grow_record
   -> write_script_headers
         -> get_rasmol_colour
   -> write_script_domains
         -> read_in_domain_data
               -> check_chain
               -> extract_domain_info (in common.c)
                     -> get_tokens (in common.c)
                     -> free_tokens (in common.c)
               -> extract_domain_data -- No longer used
                     -> get_atoken
                     -> format_residue_number (in common.c)
         -> get_domain_colour
         -> add_residue
   -> generate_script
         -> check_chain
         -> write_script_site
               -> terminate_line
         -> write_script_protein
               -> check_chain
               -> get_chain_col
               -> get_split_colours
                     -> get_colour_rgb (in common.c)
               -> add_to_list
         -> write_script_ligand
               -> terminate_line
               -> get_rasmol_colour
         -> terminate_line
         -> get_rasmol_colour
         -> write_script_metal
               -> add_to_list
         -> write_script_nucleic
               -> get_chain_col
               -> add_to_list
         -> write_script_dots
               -> get_chain_col
         -> write_script_surface
               -> get_chain_col
               -> add_to_list
         -> write_script_mask
         -> write_script_sets
   -> write_script_tail
   -> get_atom_records
         -> assign_to_object
         -> check_chain
         -> check_interaction
         -> get_atom_details
         -> create_residue_record
         -> check_for_duplicate
         -> create_conect_record
         -> create_atom_record
         -> get_max_min_coords
         -> get_sec_struc
         -> r_check_for_metal (in readpdb.c)
         -> get_element (in readpdb.c)
               -> is_hydrogen_atom (in readpdb.c)
                      -> r_check_for_metal (in readpdb.c)
         -> process_conect_record
   -> get_interacting_ligands
   -> write_secstr
   -> write_atom_records
   -> write_conect_records
   -> write_molscript
         -> get_rasmol_colour
         -> write_moline
               -> get_chain_col
         -> write_molscript_metals
               -> get_colour_rgb (in common.c)
               -> string_truncate (in common.c)
         -> write_molscript_ligands
               -> get_colour_rgb (in common.c)
   -> get_mask_vertices


Compile as follows:-
------------------

cc -c rascript.c
cc -c common.c
cc -c reapar.c
cc -c readpdb.c
cc -o rascript rascript.o common.o reapar.o readpdb.o -lm -lz

Optimization flags:-  

cc -O2 -funroll-loops -c rascript.c
cc -O2 -funroll-loops -c common.c
cc -O2 -funroll-loops -c reapar.c
cc -O2 -funroll-loops -c readpdb.c
cc -o rascript rascript.o common.o reapar.o readpdb.o -lm -lz

*/

#include "common.h"
#include "readpdb.h"
#include "rascript.h"

/* Prototypes
   ---------- */

/* In common.c */
void c_get_paths(struct c_file_data *file_name_ptr,int run_64);
int extract_domain_info(char *input_line,char *wolf_code,char *cath_code,
                        char domain_details[C_MAXDOMCHAIN][C_DOMLEN],
                        int ndetails);
int find_pdbsum_dir(char *pdb_code,char *pdbsum_template[NPDBSUM_DIR],
                    int pdb_type,char *dir_name);
void format_residue_number(char residue_number[6],char *token);
char *form_filename(char *pdb_code,char *filename_template);
void get_colour_rgb(char *c_name,double rgb[3]);
FILE *open_output_file(char *out_name,int terminate);
void store_string(char **string_ptr,char *string);
int string_truncate(char *string,int max_length);

/* Function prototypes from reapar.c */
int get_param(char *token_name,char *file_name,int filename_len,
              int exit_on_error);

/* Function prototypes from readpdb.c */
int get_element(char *atom_name,char *element);
int r_check_for_metal(char *atom_name,char *element);
int r_find_pdb_file(char *pdb_code,char *pdb_template[NPDB_DIR],
                    char *file_name,int pdb_type);
FILE *r_open_pdb_file(char *file_name,int *gz);


/***********************************************************************

get_command_arguments  -  Interpret the command line arguments

***********************************************************************/

void get_command_arguments(char *strptrs[],int num,char pdb_code[5],
                           char ch[MAX_CHAINS],int *nchain,
                           char *pdb_name,char *work_dir,
                           struct parameters *params)
{
  char chain, number_string[20];
  char temp_string[LINELEN];
 
  int i, itoken, len;

  /* Initialise variables */
  pdb_code[0] = '\0';
  pdb_name[0] = '\0';
  work_dir[0] = '\0';
  chain = ' ';
  *nchain = 0;
  for (i = 0; i < MAX_CHAINS; i++)
    ch[i] = ' ';

  /* Switch colouring-by domain option off */
  params->colour_by_domain = FALSE;
  params->all_chains = TRUE;

  /* Set default values for which objects are to be defined */
  params->define_dots = FALSE;
  params->define_ligands = TRUE;
  params->define_mask = TRUE;
  params->define_metals = TRUE;
  params->define_nucleic = TRUE;
  params->define_protein = TRUE;
  params->define_sites = TRUE;
  params->define_surfaces = FALSE;

  /* Set default values for which objects are to be displayed */
  params->show_disulphides = TRUE;
  params->show_dots = FALSE;
  params->show_ligands = TRUE;
  params->show_mask = TRUE;
  params->show_metals = TRUE;
  params->show_nucleic = TRUE;
  params->show_protein = TRUE;
  params->show_sites = TRUE;
  params->show_surfaces = FALSE;

  /* Initialise other parameters */
  params->display_option = 0;
  params->write_to_file = FALSE;
  params->write_gif_only = FALSE;
  params->coords_only = FALSE;
  params->change_MSE = FALSE;
  params->revrot = FALSE;
  params->pdb_type = STANDARD_PDB;

  /* Options for figure-generating routines */
  params->background_colour[0] = '\0';
  strcpy(params->helix_colour,"magenta");
  strcpy(params->strand_colour,"yellow");
  strcpy(params->coil_colour,"white");
  strcpy(params->dhelix_colour,"blue");
  params->colour_option = COLOUR_BY_CHAIN;
  params->output_option = OUT_RASMOL;
  for (i = 0; i < 3; i++)
    params->rotation_angle[i] = (float) 0.0;
  params->have_transformations = FALSE;
  params->metatom_radius = (float) 1.5;
  params->ligatom_radius = (float) 0.6;
  params->ligbond_radius = (float) 0.2;
  strcpy(params->ligand_colour,"CPK");
  params->run_64 = FALSE;
  params->jmol = FALSE;

  /* If nothing entered on the command line, show options available */
  if (num < 1)
    {
      printf("\n");
      printf("*** ERROR. Correct usage is:\n");
      printf("\n");
      printf("    rascript  pdb_code  [-cX]  [-dX]  [-f]  [-s]  [-oN]");
      printf("  [-NS]  [-gifonly]\n");
      printf("              [-bg colour]  [-colsst]  [-mols | -rast]");
      printf("  [-rotx angle]  [-roty angle]  [-rotz angle]\n");
      printf("              [-revrot]  [-mrad metal_radius]  ");
      printf("[-arad atom_radius]  [-brad bond_radius]\n");
      printf("              [-lcol colour]  [-dir directory]");
      printf(" [UPLOAD | MCSG]  [-64]\n");
      printf("              [-jmol]  [-pdbsum1]\n");
      printf("\n");
      printf("where\n");
      printf("  * pdb_code     = 4-character PDB code\n");
      printf("  * [-cX]        = Show single chain (X - or space).\n");
      printf("  * [-dX]        = Show single chain (X - or space) ");
      printf("coloured by domain.\n");
      printf("                   For whole structure, enter -d*\n");
      printf("  * [-f]         = Optional flag indicating that ");
      printf("output to be written\n");
      printf("                   to rascript.out file\n");
      printf("  * [-s]         = Optional flag indicating that ");
      printf("surfaces to be included\n");
      printf("  * [-oN]        = Display option, where N is the option:-\n");
      printf("                       0 = default, dot surfaces, ligands, a");
      printf("nd binding-site mask\n");
      printf("                       1 = ligand, binding-site mask, and a");
      printf("ctive-site residue\n");
      printf("  * [-NS]        = No script (ie coords only)\n");
      printf("  * [-gifonly]   = Write script that generates a GIF file\n");
      printf("  * [-bg colour] = Background colour (black, white, cream)\n");
      printf("  * [-colsst]    = Colour by secondary structure\n");
      printf("  * [-mse]       = Change MSE to MET\n");
      printf("  * [-mols]      = Output MolScript script\n");
      printf("  * [-rast]      =   ''       ''      ''   for Raster3D\n");
      printf("  * [-rotx angle]= Rotation angle about x-axis (anti-clock");
      printf("wise)\n");
      printf("  * [-roty angle]=     ''     ''    ''  y-axis   ''     ''\n");
      printf("  * [-rotz angle]=     ''     ''    ''  z-axis   ''     ''\n");
      printf("  * [-revrot]    = Reverse rotation order\n");
      printf("  * [-brad rad]  = Bond radius for ligands\n");
      printf("  * [-arad rad]  = Atom   ''    ''    ''  \n");
      printf("  * [-mrad rad]  = Atom   ''    '' metals\n");
      printf("  * [-lcol string] = Ligand colour\n");
      printf("  * [-dir directory] = Directory where PDB file (called");
      printf(" file.new is located)\n");
      printf("  * UPLOAD       = Uploaded PDB file for PDBsum generation\n");
      printf("  * MCSG         = MCSG structure\n");
      printf("  * AF           = Alpha Fold structure\n");
      printf("  * -jmol        = Output is for JSmol\n");
      printf("  * -pdbsum1     = PDBsum1 structure\n");
      printf("  * -64          = Run on 64-bit machine\n");
      printf("\n");
      printf("\n");
      printf("Example:-\n");
      printf("\n");
      printf("   rascript  1eca  -f  -s\n");
      printf("\n");
      exit(1);
    }

  /* Get the PDB code */
  strcpy(pdb_code,strptrs[1]);

  /* Loop through any remaining command arguments */
  for (itoken = 2; itoken < num + 1; itoken++)
    {
      /* Check for background colour */
      if (!strncmp(strptrs[itoken],"-bg",3))
        {
          /* Get the colour from the next command-line parameter */
          if (itoken + 1 < num + 1)
            {
              /* Check that colour name not too long */
              len = strlen(strptrs[itoken + 1]);
              if (len > COLNAM_LEN - 1)
                {
                  strncpy(temp_string,strptrs[itoken + 1],COLNAM_LEN - 1);
                  temp_string[COLNAM_LEN - 1] = '\0';
                }
              else
                strcpy(temp_string,strptrs[itoken + 1]);
              strcpy(params->background_colour,temp_string);
            }
          itoken++;
        }
      
      /* Check for option to colour protein by secondary structure */
      else if (!strncmp(strptrs[itoken],"-colsst",7))
        params->colour_option = COLOUR_BY_SST;

      /* Check whether MSE residues to be renamed MET in output file */
      else if (!strncmp(strptrs[itoken],"-mse",4))
        params->change_MSE = TRUE;

      /* Check whether MolScript output required */
      else if (!strncmp(strptrs[itoken],"-mols",5))
        params->output_option = OUT_MOLSCRIPT;

      /* Check whether MolScript output for Raster3D required */
      else if (!strncmp(strptrs[itoken],"-rast",5))
        params->output_option = OUT_RASTER3D;

      /* Check rotation order */
      else if (!strncmp(strptrs[itoken],"-revrot",7))
        params->revrot = TRUE;

      /* Check whether a rotation angle has been entered */
      else if (!strncmp(strptrs[itoken],"-rot",4))
        {
          /* Determine which angle it is */
          if (strptrs[itoken][4] == 'x')
            i = 0;
          else if (strptrs[itoken][4] == 'y')
            i = 1;
          else if (strptrs[itoken][4] == 'z')
            i = 2;
          else
            i = -1;

          /* If a valid keyword, then retrieve the angle from then
             next token */
          if (i > -1)
            {
              /* Get the value string and convert to float */
              if (itoken + 1 < num + 1)
                {
                  /* Check that string not too long */
                  len = strlen(strptrs[itoken + 1]);
                  if (len > 19)
                    {
                      strncpy(temp_string,strptrs[itoken + 1],19);
                      temp_string[19] = '\0';
                    }
                  else
                    strcpy(temp_string,strptrs[itoken + 1]);
                  strcpy(number_string,temp_string);

                  /* Store the value */
                  params->rotation_angle[i] = (float) atof(number_string);

                  /* Set flag that we have transformations */
                  params->have_transformations = TRUE;
                }
              itoken++;
            }
        }

      /* Check for ligand colour */
      else if (!strncmp(strptrs[itoken],"-lcol",5))
        {
          /* Get the string from the next command-line parameter */
          if (itoken + 1 < num + 1)
            {
              /* Store string provided not too long */
              len = strlen(strptrs[itoken + 1]);
              if (len > COLNAM_LEN)
                {
                  strncpy(temp_string,strptrs[itoken + 1],COLNAM_LEN);
                  temp_string[COLNAM_LEN] = '\0';
                }
              else
                strcpy(temp_string,strptrs[itoken + 1]);
              strcpy(params->ligand_colour,temp_string);

              /* Increment token count */
              itoken++;
            }
        }

      /* Check whether a ligand bond radius has been entered */
      else if (!strncmp(strptrs[itoken],"-brad",5))
        {
          /* Get the value string and convert to float */
          if (itoken + 1 < num + 1)
            {
              /* Check that string not too long */
              len = strlen(strptrs[itoken + 1]);
              if (len > 19)
                {
                  strncpy(temp_string,strptrs[itoken + 1],19);
                  temp_string[19] = '\0';
                }
              else
                strcpy(temp_string,strptrs[itoken + 1]);
              strcpy(number_string,temp_string);

              /* Store the value */
              params->ligbond_radius = (float) atof(number_string);
            }
          itoken++;
        }

      /* Check whether a ligand atom radius has been entered */
      else if (!strncmp(strptrs[itoken],"-arad",5))
        {
          /* Get the value string and convert to float */
          if (itoken + 1 < num + 1)
            {
              /* Check that string not too long */
              len = strlen(strptrs[itoken + 1]);
              if (len > 19)
                {
                  strncpy(temp_string,strptrs[itoken + 1],19);
                  temp_string[19] = '\0';
                }
              else
                strcpy(temp_string,strptrs[itoken + 1]);
              strcpy(number_string,temp_string);

              /* Store the value */
              params->ligatom_radius = (float) atof(number_string);
            }
          itoken++;
        }

      /* Check for directory holding PDB file */
      else if (!strncmp(strptrs[itoken],"-dir",4))
        {
          /* Get the directory name */
          if (itoken + 1 < num + 1)
            {
              strcpy(work_dir,strptrs[itoken + 1]);

              /* Form name of PDB file */
              strcpy(pdb_name,work_dir);
              strcat(pdb_name,"/file.new");
            }
          itoken++;
        }

      /* Check whether a metal atom radius has been entered */
      else if (!strncmp(strptrs[itoken],"-mrad",5))
        {
          /* Get the value string and convert to float */
          if (itoken + 1 < num + 1)
            {
              /* Check that string not too long */
              len = strlen(strptrs[itoken + 1]);
              if (len > 19)
                {
                  strncpy(temp_string,strptrs[itoken + 1],19);
                  temp_string[19] = '\0';
                }
              else
                strcpy(temp_string,strptrs[itoken + 1]);
              strcpy(number_string,temp_string);

              /* Store the value */
              params->metatom_radius = (float) atof(number_string);
            }
          itoken++;
        }

      /* Check for the -d option indicating colouring by domain required */
      else if (!strncmp(strptrs[itoken],"-d",2))
        {
          params->colour_by_domain = TRUE;

          /* Check whether have a chain-id, or whole structure is
             required */
          if (strptrs[itoken][2] == '*')
            {
              params->all_chains = TRUE;
              chain = '\0';
            }
          else if (strptrs[itoken][2] != '\0')
            {
              /* Get chain */
              chain = strptrs[itoken][2];
            }
          else
            {
              /* Take chain to be blank */
              chain = ' ';
            }

          /* If valid chain, then add to list of chains */
          if (chain != '\0' && *nchain < MAX_CHAINS - 1)
            {
              ch[*nchain] = chain;
              (*nchain)++;

              /* Set flag to indicate that specific chains have
                 been selected */
              params->all_chains = FALSE;
            }
        }

      /* Check for -c option preceding the chain-identifier */
      else if (!strncmp(strptrs[itoken],"-c",2))
        {
          /* Retrieve the chain identifier */
          if (strptrs[itoken][2] != '\0')
            chain = strptrs[itoken][2];
          else
            chain = ' ';

          /* If can, then store in list of chains */
          if (*nchain < MAX_CHAINS - 1)
            {
              ch[*nchain] = chain;
              (*nchain)++;

              /* Set flag to indicate that specific chains have
                 been selected */
              params->all_chains = FALSE;
            }
        }

      /* Check for the -f option */
      else if (!strncmp(strptrs[itoken],"-f",2))
        params->write_to_file = TRUE;

      /* Check for the -s option */
      else if (!strncmp(strptrs[itoken],"-s",2))
        {
          params->define_dots = TRUE;
          params->define_mask = TRUE;
          params->define_surfaces = TRUE;
          params->show_dots = TRUE;
          params->show_mask = TRUE;
          params->show_surfaces = TRUE;
        }

      /* Check for the -o flag giving the display option */
      else if (!strncmp(strptrs[itoken],"-o",2))
        {
          /* Get the option number */
          number_string[0] = strptrs[itoken][2];
          number_string[1] = '\0';
          params->display_option = atoi(number_string);
        }

      /* Check for the -gifonly option used to generate gif image of the
         start-up view for inclusion in Web pages */
      else if (!strncmp(strptrs[itoken],"-gifonly",8))
        params->write_gif_only = TRUE;

      /* Check for the -NS no-script option (ie output coordinates only,
         with no script commands */
      else if (!strncmp(strptrs[itoken],"-NS",3))
        params->coords_only = TRUE;

      /* Check for the UPLOAD option */
      else if (!strcmp(strptrs[itoken],"UPLOAD"))
        params->pdb_type = UPLOADED_PDB;

      /* Check for the MCSG option */
      else if (!strcmp(strptrs[itoken],"MCSG"))
        params->pdb_type = MCSG_PDB;

      /* Check for the Alpha Fold option */
      else if (!strcmp(strptrs[itoken],"AF"))
        params->pdb_type = ALPHA_FOLD_PDB;

      /* Check for the PDBsum1 option */
      else if (!strcmp(strptrs[itoken],"-pdbsum1"))
        params->pdb_type = PDBSUM1;

      /* Check for the 64-bit option */
      if (!strcmp(strptrs[itoken],"-64"))
        params->run_64 = TRUE;

      /* Check for the JSmol option */
      if (!strcmp(strptrs[itoken],"-jmol"))
        params->jmol = TRUE;
    }

  /* If no chains identified, or only chain blank, then assume all
     chains are required */
  if (*nchain == 0 || (*nchain == 1 && ch[0] == ' '))
    {
      params->all_chains = TRUE;
      *nchain = 0;
    }
  
  /* If showing surfaces, then switch off all but the ligands and
     active sites */
  if (params->show_surfaces == TRUE)
    {
      /* Switch off unwanted objects */
      params->show_nucleic = FALSE;
      params->show_protein = FALSE;
      params->show_disulphides = FALSE;

      /* Switch other display objects depending on display option */
      if (params->display_option == 1)
        {
          params->show_dots = FALSE;
        }
    }

  /* If colouring by domain, switch all things off */
  if (params->colour_by_domain == TRUE)
    {
      params->define_nucleic = FALSE;
      params->define_protein = FALSE;
      params->define_ligands = FALSE;
      params->define_metals = FALSE;
      params->define_sites = FALSE;
      params->show_nucleic = FALSE;
      params->show_protein = FALSE;
      params->show_disulphides = FALSE;
      params->show_ligands = FALSE;
      params->show_metals = FALSE;
      params->show_sites = FALSE;
    }

  /* Hard-coded fix for a bug in RasMol on 2ckb */
  if (!strncmp(pdb_code,"2ckb",4))
    params->show_disulphides = FALSE;
}
/***********************************************************************

form_filenames - Form the various file names required by the program

***********************************************************************/

void form_filenames(char *pdb_code,char *pdbsum_dir,char *pdb_name,
                    char *vtx_name,char *promotif_dir,
                    char *cath_domains_file,
                    struct c_file_data file_name_ptr,
                    struct parameters params)
{
  char dir[20], letter[3], prefix[20], suffix[20];
  char filename[FILENAME_LEN];
  char first_pdbname[FILENAME_LEN];

  char *fname, *string_ptr;

  int error_status;
  int dir_type, ipdb;

  static int exit_on_error = FALSE;

  FILE *file_ptr;

  /* Initialise flag */
  first_pdbname[0] = '\0';
  strcpy(pdbsum_dir,"./");

  /* For PDBsum1, form name of the .new file */
  if (params.pdb_type == PDBSUM1)
    dir_type = PDBSUM1;

  /* If already have the required PDB filename, use that */
  else if (pdb_name[0] != '\0')
    strcpy(filename,pdb_name);

  /* Otherwise, use given PDB code */
  else
    dir_type = r_find_pdb_file(pdb_code,file_name_ptr.pdb_template,pdb_name,
                               params.pdb_type);

  /* If file not found, then show error message */
  if (dir_type == NOT_FOUND)
    {
      printf("\n*** Unable to find PDB file for code [%s]\n",pdb_code);
      exit(1);
    }

  if (params.pdb_type == PDBSUM1)
    {
      /* Get the name of the PDBsum directory */
      pdbsum_dir
        = form_filename(pdb_code,
                        file_name_ptr.other_dir[PDBSUM_TEMPLATE]);

      /* Form the name of the Promofit directory and PDB file name */
      strcat(promotif_dir,"promotif/");
      sprintf(pdb_name,"%s/%s.new",pdbsum_dir,pdb_code);
    }

  else
    {
      /* Form name of PDBsum path */
      dir_type = find_pdbsum_dir(pdb_code,file_name_ptr.pdbsum_template,
                                 params.pdb_type,pdbsum_dir);

      /* Form name of PROMOTIF paths */
      strcpy(promotif_dir,pdbsum_dir);
    }

  /* CATH domains file */
  strcpy(cath_domains_file,file_name_ptr.other_dir[CATH_DOMAINS]);

  /* Form name of input binding-site mask coordinates file */
  strcpy(vtx_name,pdbsum_dir);
  strcat(vtx_name,"mask.vtx");
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
  struct rasdef *rasdef_ptr;
  struct rasdef *first_rasdef_ptr, *last_rasdef_ptr;

  /* Initialise rasdef pointers */
  rasdef_ptr = NULL;
  first_rasdef_ptr = *fst_rasdef_ptr;
  last_rasdef_ptr = *lst_rasdef_ptr;

  /* Create rasdef entry */
  rasdef_ptr =
    (struct rasdef*)malloc(sizeof(struct rasdef));
  if (rasdef_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for");
      printf(" struct rasdef\n");
      printf("***        Program rasurf terminated with");
      printf(" error.\n");
      exit (1);
    }

  /* Store the details */
  rasdef_ptr->deleted = FALSE;
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
  rasdef_ptr->atom_limits.first = TRUE;
  rasdef_ptr->first_atom_ptr = NULL;
  rasdef_ptr->last_atom_ptr = NULL;
  strcpy(rasdef_ptr->ligand_from,ligand_from);
  strcpy(rasdef_ptr->ligand_to,ligand_to);

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

format_rasdef_number  -  Interpret the residue-number from the rasdef line

***********************************************************************/

void format_rasdef_number(char *tmpstr,char *res_num)
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

void get_ligand_residues(char *input_line,struct rasdef **fst_rasdef_ptr,
                         struct rasdef **lst_rasdef_ptr,int number)
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
      token_len = get_token(input_line,LINELEN,&lpos,token_string,
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
          format_rasdef_number(number_string,res_num);

          /* Copy chain into chain-string */
          chain_str[0] = chain;
          chain_str[1] = '\0';

          /* If this is the first residue of this ligand, store as first */
          if (first_residue[0] == '\0')
            {
              /* Copy across and add chain-id */
              strcpy(first_residue,chain_str);
              strcat(first_residue,number_string);
            }

          /* Copy residue number and chain-id across as last residue of
             current ligand */
          strcpy(last_residue,chain_str);
          strcat(last_residue,number_string);

          /* Increment count of residues belonging to this ligand */
          count++;

          /* Blank out unwanted fields */
          ligand_from[0] = '\0';
          ligand_to[0] = '\0';

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
                              input_line,first_residue,last_residue,
                              count);
    }

  /* Return current pointers */
  *fst_rasdef_ptr = first_rasdef_ptr;
  *lst_rasdef_ptr = last_rasdef_ptr;
}
/***********************************************************************

store_chain  -  Store current chain in list of this protein's chains

***********************************************************************/

void store_chain(char chain,char *last_chain,char *chains,int *nchains)
{
  int ipos;
  int got_it;

  /* If the chain is a new one, add it to the list
     of chains */
  if (chain != *last_chain)
    {
      /* Check whether we already have this chain */
      got_it = FALSE;
      for (ipos = 0; ipos < (*nchains); ipos++)
        if (chains[ipos] == chain)
          got_it = TRUE;

      /* Store this chain */
      if (got_it == FALSE && (*nchains) < MAX_CHAINS)
        {
          chains[*nchains] = chain;
          (*nchains)++;
        }
    }
}
/***********************************************************************

get_dfn_data  -  Read in and store the RasMol definitions from the
                 rasmol.dfn file

***********************************************************************/

int get_dfn_data(char *pdb_code,char *pdbsum_dir,char *dfn_name,
                 struct rasdef **fst_rasdef_ptr,char *chains,int *nchn,
                 char *work_dir,struct parameters params)
{
  char file_name[FILENAME_LEN], input_line[LINELEN + 1];
  char chain, last_chain, res_name[4], res_num[6];
  char colour[COLNAM_LEN];
  char ligand_from[7], ligand_to[7];

  int count, line, ndef, nchains, nligand, nmetal, type;
  int first, have_dfn;

  struct rasdef *rasdef_ptr;
  struct rasdef *first_rasdef_ptr, *last_rasdef_ptr;

  FILE *file_ptr;

  /* Initialise */
  first = TRUE;
  count = 0;
  have_dfn = FALSE;
  last_chain = '\0';
  line = 0;
  ligand_from[0] = '\0';
  ligand_to[0] = '\0';
  ndef = 0;
  nchains = *nchn;
  nligand = 0;
  nmetal = 0;

  /* Initialise rasdef pointers */
  rasdef_ptr = NULL;
  first_rasdef_ptr = NULL;
  last_rasdef_ptr = NULL;

  /* Form filename of the rasmol.dfn file */
  strcpy(file_name,pdbsum_dir);
  if (params.pdb_type == PDBSUM1)
    strcpy(file_name,"../newsec/");
  strcat(file_name,dfn_name);
  if (work_dir[0] != '\0')
    {
      strcpy(file_name,work_dir);
      strcat(file_name,DIR_SEP);
      strcat(file_name,dfn_name);
    }

  /* Open the rasmol.dfn input file */
  if ((file_ptr = fopen(file_name,"r")) !=  NULL)
    {
      /* Loop through the file processing each record */
      while (fgets(input_line,LINELEN,file_ptr)!= NULL)
        {
          /* Increment line-count */
          line++;

          /* Truncate the string to strip off the newline */
          string_truncate(input_line,LINELEN);

          /* Check line-type according to 1st character: C, L, M,
             N or P */
          /* If protein CA-only chain, then flag as such */
          if (input_line[0] == 'C')
            {
              /* Define type */
              type = CA_ONLY;

              /* Define residue details */
              res_name[0] = '\0';
              res_num[0] = '\0';
              chain = input_line[1];
              if (chain == '\0' || chain == '\n')
                chain = ' ';
              colour[0] = '\0';

              /* If reading in the entire structure, store current 
                 chain in chains list */
              if (params.all_chains == TRUE)
                store_chain(chain,&last_chain,chains,&nchains);

              /* Create a rasdef entry */
              rasdef_ptr
                = create_rasdef_entry(&first_rasdef_ptr,
                                      &last_rasdef_ptr,type,ndef,
                                      res_name,res_num,chain,colour,
                                      input_line,ligand_from,ligand_to,
                                      count);

              /* Increment count */
              ndef++;
            }

          /* Check for ligand identifier */
          else if (input_line[0] == 'L')
            {
              /* Increment ligand-count */
              nligand++;

              /* Extract the ligand residues */
              get_ligand_residues(input_line+4,&first_rasdef_ptr,
                                  &last_rasdef_ptr,nligand);
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
              if (input_line[4] == ' ')
                strcpy(colour,input_line+5);
              else
                strcpy(colour,input_line+4);

              /* Create a rasdef entry */
              rasdef_ptr
                = create_rasdef_entry(&first_rasdef_ptr,
                                      &last_rasdef_ptr,type,ndef,
                                      res_name,res_num,chain,colour,
                                      input_line,ligand_from,ligand_to,
                                      count);

              /* Increment count */
              ndef++;
            }
          
          /* Check for nucleic acid chain */
          else if (input_line[0] == 'N')
            {
              /* Define type */
              type = NUCLEIC_ACID;

              /* Define residue details */
              res_name[0] = '\0';
              res_num[0] = '\0';
              chain = input_line[1];
              if (chain == '\0' || chain == '\n')
                chain = ' ';
              colour[0] = '\0';

              /* If reading in the entire structure, store current 
                 chain in chains list */
              if (params.all_chains == TRUE)
                store_chain(chain,&last_chain,chains,&nchains);

              /* Create a rasdef entry */
              rasdef_ptr
                = create_rasdef_entry(&first_rasdef_ptr,
                                      &last_rasdef_ptr,type,ndef,
                                      res_name,res_num,chain,colour,
                                      input_line,ligand_from,ligand_to,
                                      count);

              /* Increment count */
              ndef++;
            }

          /* If protein chain, then write the script records to define
             and display the chain */
          else if (input_line[0] == 'P')
            {
              /* Define type */
              type = PROTEIN;

              /* Define residue details */
              res_name[0] = '\0';
              res_num[0] = '\0';
              chain = input_line[1];
              if (chain == '\0' || chain == '\n')
                chain = ' ';
              colour[0] = '\0';

              /* If reading in the entire structure, store current 
                 chain in chains list */
              if (params.all_chains == TRUE)
                store_chain(chain,&last_chain,chains,&nchains);

              /* Create a rasdef entry */
              rasdef_ptr
                = create_rasdef_entry(&first_rasdef_ptr,
                                      &last_rasdef_ptr,type,ndef,
                                      res_name,res_num,chain,colour,
                                      input_line,ligand_from,ligand_to,
                                      count);

              /* Increment count */
              ndef++;
            }
        }

      /* Close the file */
      fclose(file_ptr);

      /* Determine whether .dfn data read in */
      if (line > 0)
        have_dfn = TRUE;
    }

  /* Return pointer to start of RasMol definitions */
  *fst_rasdef_ptr = first_rasdef_ptr;

  /* Return number of protein chains for which there will be PROMOTIF
     data */
  *nchn = nchains;

  /* Return whether definitions read in */
  return(have_dfn);
}
/***********************************************************************

get_sst_dat  -  Read in the doiminant secondary structures in each chain
                from the sst.dat file

***********************************************************************/

int get_sst_dat(char *pdb_code,char *pdbsum_dir,char *chains,
                char *dom_sst,int nchains,struct parameters params)
{
  char file_name[FILENAME_LEN], input_line[LINELEN + 1];
  char chain, dominant_sst;

  int ichain, len, line;
  int found, have_sstdat;

  FILE *file_ptr;

  /* Initialise */
  have_sstdat = FALSE;
  line = 0;

  /* Form filename of the sst.dat file */
  strcpy(file_name,pdbsum_dir);
  if (params.pdb_type == PDBSUM1)
    strcat(file_name,"newsec/");
  strcat(file_name,"sst.dat");

  /* Open the input file */
  if ((file_ptr = fopen(file_name,"r")) !=  NULL)
    {
      /* Loop through the file processing each record */
      while (fgets(input_line,LINELEN,file_ptr)!= NULL)
        {
          /* Increment line-count */
          line++;

          /* Truncate the string to strip off the newline */
          len = string_truncate(input_line,LINELEN);

          /* Get the chain identifier */
          chain = input_line[0];

          /* Get the dominant secondary structure for this chain */
          dominant_sst = input_line[len - 1];

          /* Find the right chain record */
          found = FALSE;
          for (ichain = 0; ichain < nchains && found == FALSE; ichain++)
            {
              /* If this is the right chain, save dominant secondary 
                 structure */
              if (chain == chains[ichain])
                {
                  dom_sst[ichain] = dominant_sst;
                  found = TRUE;
                }
            }
        }

      /* Close the file */
      fclose(file_ptr);

      /* Set flag that have sst.dat file */
      have_sstdat = TRUE;
    }

  /* Return whether data read in */
  return(have_sstdat);
}
/***********************************************************************

create_promotif_record  -  Create and initialise a new promotif record

***********************************************************************/

struct promotif *create_promotif_record(struct promotif **fst_promotif_ptr,
                                        struct promotif **lst_promotif_ptr,
                                        int data_type,char sec_struc,
                                        char *res_name1,char *res_name2,
                                        char *res_num1,char *res_num2,
                                        char chain1,char chain2,int class)
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
  promotif_ptr->data_type = data_type;
  promotif_ptr->sec_struc = sec_struc;
  strcpy(promotif_ptr->res_name1,res_name1);
  strcpy(promotif_ptr->res_name2,res_name2);
  strcpy(promotif_ptr->res_num1,res_num1);
  strcpy(promotif_ptr->res_num2,res_num2);
  promotif_ptr->chain1 = chain1;
  promotif_ptr->chain2 = chain2;
  promotif_ptr->class = class;

  /* Initialise pointer to next record */
  promotif_ptr->next_promotif_ptr = NULL;

  /* Return current pointers */
  *fst_promotif_ptr = first_promotif_ptr;
  *lst_promotif_ptr = last_promotif_ptr;

  /* Return pointer to the new promotif record */
  return(promotif_ptr);
}
/***********************************************************************

read_promotif_dat  -  Read in the PROMOTIF analyses from the promotif.dat
                      file, if present

***********************************************************************/

int read_promotif_dat(char *pdb_code,
                      struct promotif **fst_promotif_ptr,
                      char *promotif_dir,char *chains,int nchains)
{
#define NMOTIFS 5

  char file_name[FILENAME_LEN], input_line[LINELEN];
  char chain, chain1, chain2, res_num1[6], res_num2[6];
  char res_name1[4], res_name2[4];

  int class, data_end, data_start, i, ichain, ifile, line, motif;
  int have_file;

  struct promotif *last_promotif_ptr, *first_promotif_ptr;
  struct promotif *promotif_ptr;

  FILE *file_ptr;

  static char *motif_name[NMOTIFS] = {
    "HELICES", "STRANDS", "BETA_TURNS", "GAMMA_TURNS", "DISULPHIDES"
  };
  static char sec_struc[]={"HEbgD"};

  /* Initialise pointers */
  first_promotif_ptr = NULL;
  last_promotif_ptr = NULL;

  /* Initialise variables */
  class = 1;
  have_file = FALSE;
  motif = -1;
  strcpy(res_name1,"XXX");
  strcpy(res_name2,"XXX");

  /* Form name of promotif.dat file */
  strcpy(file_name,promotif_dir);
  strcat(file_name,DIR_SEP);
  strcat(file_name,pdb_code);
  strcat(file_name,DIR_SEP);
  strcat(file_name,"promotif.dat");

  /* Open the PROMOTIF file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    have_file = FALSE;

  /* If have the file, read in all its data */
  if (have_file == TRUE)
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
                if (!strcmp(input_line+6,motif_name[i]))
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

              /* Process according to the different record formats */
              switch (motif)
                {
                  /* a. Helix */
                case 0 :
                  strncpy(res_num1,input_line+8,5);
                  res_num1[5] = '\0';
                  strncpy(res_num2,input_line+20,5);
                  res_num2[5] = '\0';
                  break;

                  /* b. Beta strand */
                case 1 :
                  strncpy(res_num1,input_line+8,5);
                  res_num1[5] = '\0';
                  strncpy(res_num2,input_line+20,5);
                  res_num2[5] = '\0';
                  break;

                  /* c. Beta turns */
                case 2 :
                  strncpy(res_num1,input_line+3,5);
                  res_num1[5] = '\0';
                  strncpy(res_num2,input_line+39,5);
                  res_num2[5] = '\0';
                  break;

                  /* d. Gamma turns */
                case 3 :
                  strncpy(res_num1,input_line+3,5);
                  res_num1[5] = '\0';
                  strncpy(res_num2,input_line+27,5);
                  res_num2[5] = '\0';
                  break;

                  /* e. Disulphide bridge */
                case 4 :
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

                  /* f. Beta hairpin */
                case 5 :
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
                                           ifile,sec_struc[motif],
                                           res_name1,res_name2,
                                           res_num1,res_num2,
                                           chain1,chain2,class);
                }
            }
        }

      /* Close the file */
      fclose(file_ptr);
    }

  /* Return pointer to first PROMOTIF entry */
  *fst_promotif_ptr = first_promotif_ptr;

  /* Return flag indicating whether promotif.dat file found */
  return(have_file);
}
/***********************************************************************

read_in_promotif_data  -  Read in the PROMOTIF analyses for the given
                          structure, if present

***********************************************************************/

void read_in_promotif_data(char *pdb_code,
                           struct promotif **fst_promotif_ptr,
                           char *promotif_dir,char *chains,int nchains)
{
  char chain_str[2], current_chain, file_name[FILENAME_LEN];
  char input_line[LINELEN];
  char chain1, chain2, res_num1[6], res_num2[6];
  char res_name1[4], res_name2[4];

  int class, data_end, data_start, have_file, ichain, ifile, line;
  int npromotif_files;

  struct promotif *last_promotif_ptr, *first_promotif_ptr;
  struct promotif *promotif_ptr;

  FILE *file_ptr;

  static char *file_extension[] = {
    "hlx", "str", "btn", "gtn", "dsf" };
  static char sec_struc[]={"HEbgD"};

  /* Initialise pointers */
  first_promotif_ptr = NULL;
  last_promotif_ptr = NULL;

  /* Initialise variables */
  class = 1;
  npromotif_files = strlen(sec_struc);
  strcpy(res_name1,"XXX");
  strcpy(res_name2,"XXX");

  /* Loop over all the chains to be processed */
  for (ichain = 0; ichain < nchains; ichain++)
    {
      /* Pick up the current chain */
      current_chain = chains[ichain];
      chain_str[0] = current_chain;
      chain_str[1] = '\0';

      /* Loop over the different type of PROMOTIF data that may be present */
      for (ifile = 0; ifile < npromotif_files; ifile++)
        {
          /* Initialise flag indicating whether PROMOTIF file present */
          have_file = TRUE;
          line = 0;

          /* Form the current filename */
          strcpy(file_name,promotif_dir);
          strcat(file_name,DIR_SEP);
          strcat(file_name,pdb_code);
          strcat(file_name,DIR_SEP);
          strcat(file_name,pdb_code);
          if (current_chain != ' ')
            strcat(file_name,chain_str);
          strcat(file_name,"_");
          strcat(file_name,file_extension[ifile]);
          strcat(file_name,".html");

          /* Open the PROMOTIF file */
          if ((file_ptr = fopen(file_name,"r")) ==  NULL)
            have_file = FALSE;

          /* If have the file, read in all its data */
          if (have_file == TRUE)
            {
              /* Initialise variables */
              data_start = FALSE;
              data_end = FALSE;

              /* Loop while reading in records from the file */
              while (fgets(input_line,LINELEN,file_ptr) != NULL &&
                     data_end == FALSE)
                {
                  /* If waiting for the start of the data, check for dashes */
                  if (data_start == FALSE)
                    {
                      if (!strncmp(input_line,"----------",10))
                        data_start = TRUE;
                    }

                  /* Otherwise, get this line of data */
                  else
                    {
                      /* Check for end of data */
                      if (!strncmp(input_line,"</PRE>",6) ||
                          !strncmp(input_line," <p>",4))
                        data_end = TRUE;
                      
                      /* Read in the data according to the file type */
                      else if ((int) strlen(input_line) > 10)
                        {
                          /* Initialise data */
                          chain1 = current_chain;
                          chain2 = current_chain;
                          res_num1[0] = '\0';
                          res_num2[0] = '\0';

                          /* Process according to the different file formats */
                          switch (ifile)
                            {
                              /* a. Helix */
                            case 0 :
                              strncpy(res_num1,input_line+10,5);
                              res_num1[5] = '\0';
                              strncpy(res_num2,input_line+16,5);
                              res_num2[5] = '\0';
                              break;

                              /* b. Beta strand */
                            case 1 :
                              strncpy(res_num1,input_line+11,5);
                              res_num1[5] = '\0';
                              strncpy(res_num2,input_line+17,5);
                              res_num2[5] = '\0';
                              break;

                              /* c. Beta or gamma turns */
                            case 2 :
                            case 3 :
                              strncpy(res_num1,input_line,5);
                              res_num1[5] = '\0';
                              strncpy(res_num2,input_line+6,5);
                              res_num2[5] = '\0';
                              break;

                              /* d. Disulphide bridge */
                            case 4 :
                              chain1 = input_line[5];
                              strncpy(res_num1,input_line+6,5);
                              res_num1[5] = '\0';
                              chain2 = input_line[18];
                              strncpy(res_num2,input_line+19,5);
                              res_num2[5] = '\0';

                              /* Check which of these belong to the
                                 same chain */
                              if (current_chain != ' ')
                                {
                                  if (chain1 != ' ' && chain1 != current_chain)
                                    res_num1[0] = '\0';
                                  if (chain2 != ' ' && chain2 != current_chain)
                                    res_num2[0] = '\0';
                                }
                              break;

                              /* e. Beta hairpin */
                            case 5 :
                              strncpy(res_num1,input_line+6,5);
                              res_num1[5] = '\0';
                              strncpy(res_num2,input_line+21,5);
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
                                                       ifile,sec_struc[ifile],
                                                       res_name1,res_name2,
                                                       res_num1,res_num2,
                                                       chain1,chain2,class);
                            }
                        }
                    }
                }

              /* Close the file */
              fclose(file_ptr);
            }
        }
    }

  /* Return pointer to first PROMOTIF entry */
  *fst_promotif_ptr = first_promotif_ptr;
}
/***********************************************************************

get_promotif_data  -  Read in the PROMOTIF analyses for the given
                      structure, if present

***********************************************************************/

void get_promotif_data(char *pdb_code,
                       struct promotif **fst_promotif_ptr,
                       char *promotif_dir,char *chains,int nchains)
{
  int have_file;

  struct promotif *first_promotif_ptr;

  /* Initialise pointer */
  first_promotif_ptr = *fst_promotif_ptr;

  /* Read in from new promotof.dat file, if present */
  have_file
    = read_promotif_dat(pdb_code,&first_promotif_ptr,promotif_dir,
                        chains,nchains);

  /* If promotif.dat not present, pick up the data from the old-style
     PROMOTIF .html files */
  if (have_file == FALSE)
    read_in_promotif_data(pdb_code,&first_promotif_ptr,promotif_dir,
                          chains,nchains);

  /* Return pointer to first PROMOTIF entry */
  *fst_promotif_ptr = first_promotif_ptr;
}
/***********************************************************************

read_rin_file  -  Read in the secondary structure definitions from the
                  file.rin file

***********************************************************************/

void read_rin_file(struct promotif **fst_promotif_ptr,char *work_dir)
{
  char file_name[FILENAME_LEN];
  char input_line[LINELEN];
  char chain, last_chain, res_num[6], res_num1[6], res_num2[6];
  char res_name[4], res_name1[4], res_name2[4];
  char last_secstr, secstr;

  int class;
  int have_file, length;

  struct promotif *last_promotif_ptr, *first_promotif_ptr;
  struct promotif *promotif_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  class = 1;
  last_secstr = ' ';
  last_chain = ' ';
  length = 0;
  res_num1[0] = '\0';
  res_num2[0] = '\0';
  strcpy(res_name1,"XXX");
  strcpy(res_name2,"XXX");

  /* Initialise pointers */
  first_promotif_ptr = NULL;
  last_promotif_ptr = NULL;

  /* Form the filename */
  strcpy(file_name,work_dir);
  strcat(file_name,DIR_SEP);
  strcat(file_name,"file.rin");

  /* Open the .rin file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    have_file = FALSE;
  else
    have_file = TRUE;

  /* If have the file, read in all its data */
  if (have_file == TRUE)
    {
      /* Loop while reading in records from the file */
      while (fgets(input_line,LINELEN,file_ptr) != NULL)
        {
          /* Get residue name */
          strncpy(res_name,input_line+4,3);
          res_name[3] = '\0';

          /* Get residue number */
          strncpy(res_num,input_line+9,5);
          res_num[5] = '\0';

          /* Get the secondary structure and chain-id */
          secstr = input_line[14];
          chain = input_line[8];

          /* Change secondary structure where necessary */
          if (secstr == 'g' || secstr == 'G' || secstr == 'i' ||
              secstr == 'I')
            secstr = 'H';
          if (secstr != 'H' && secstr != 'E')
            secstr = ' ';

          /* If either secondary structure or chain have changed,
             check for start or end */
          if (secstr != last_secstr || chain != last_chain)
            {
              /* If previous secondary structure helix or strand,
                 then store start and end residues */
              if (last_secstr != ' ' && length > 2)
                {
                  /* Create a PROMOTIF record */
                  promotif_ptr =
                    create_promotif_record(&first_promotif_ptr,
                                           &last_promotif_ptr,
                                           0,last_secstr,
                                           res_name1,res_name2,
                                           res_num1,res_num2,
                                           last_chain,last_chain,class);
                }

              /* If have secondary structure, start a new entry */
              if (secstr != ' ')
                {
                  /* Store the start residue name and number */
                  strcpy(res_name1,res_name);
                  strcpy(res_num1,res_num);

                  /* Initialise length of secondary structure element */
                  length = 0;
                }

              /* Store current secondary structure and chain-id */
              last_secstr = secstr;
              last_chain = chain;
            }

          /* If have secondary structure details, then store current
             residue as the end-residue */
          if (secstr != ' ')
            {
              /* Store the end residue name and number */
              strcpy(res_name2,res_name);
              strcpy(res_num2,res_num);

              /* Increment length of secondary structure element */
              length++;
            }
        }

      /* Close the file */
      fclose(file_ptr);

      /* If have a final entry to stroe, do so */
      if (last_secstr != ' ' && length > 2)
        {
          /* Create a PROMOTIF record */
          promotif_ptr =
            create_promotif_record(&first_promotif_ptr,
                                   &last_promotif_ptr,
                                   0,last_secstr,
                                   res_name1,res_name2,
                                   res_num1,res_num2,
                                   last_chain,last_chain,class);
        }
    }

  /* Return pointer to first PROMOTIF entry */
  *fst_promotif_ptr = first_promotif_ptr;
}
/***********************************************************************

create_grow_record  -  Create and initialise a new grow record

***********************************************************************/

struct grow *create_grow_record(struct grow **fst_grow_ptr,
                struct grow **lst_grow_ptr,
                int data_type,
                char *res_name1,char *res_num1,char chain1,
                char *res_name2,char *res_num2,char chain2)
{
  struct grow *grow_ptr, *first_grow_ptr, *last_grow_ptr;

  /* Initialise grow pointers */
  grow_ptr = NULL;
  first_grow_ptr = *fst_grow_ptr;
  last_grow_ptr = *lst_grow_ptr;

  /* Allocate memory for structure to hold grow info */
  grow_ptr = (struct grow *) malloc(sizeof(struct grow));
  if (grow_ptr == NULL)
    {
      printf("*** Can't allocate memory for struct grow\n");
      exit (1);
    }

  /* If this is the very first grow, then store pointer */
  if (first_grow_ptr == NULL)
    first_grow_ptr = grow_ptr;

  /* Add link from previous grow to the current one */
  if (last_grow_ptr != NULL)
    last_grow_ptr->next_grow_ptr = grow_ptr;
  last_grow_ptr = grow_ptr;

  /* Transfer the data fields to the new grow record */
  grow_ptr->data_type = data_type;
  strcpy(grow_ptr->res_name1,res_name1);
  strcpy(grow_ptr->res_name2,res_name2);
  strcpy(grow_ptr->res_num1,res_num1);
  strcpy(grow_ptr->res_num2,res_num2);
  grow_ptr->chain1 = chain1;
  grow_ptr->chain2 = chain2;

  /* Initialise pointer to next record */
  grow_ptr->next_grow_ptr = NULL;

  /* Return current pointers */
  *fst_grow_ptr = first_grow_ptr;
  *lst_grow_ptr = last_grow_ptr;

  /* Return pointer to the new grow record */
  return(grow_ptr);
}
/***********************************************************************

read_in_grow_data  -  Read in the GROW data on contacts with ligands
                      and metals for the given structure

***********************************************************************/

int read_in_grow_data(char *pdb_code,char *pdbsum_dir,char *work_dir,
                      char *grow_name,struct grow **fst_grow_ptr)
{
  char file_name[FILENAME_LEN];
  char input_line[LINELEN];
  char chain1, chain2;
  char res_name1[4], res_name2[4], res_num1[6], res_num2[6];

  int data_type, line, match1, match2;
  int have_grow, match;

  struct grow *last_grow_ptr, *first_grow_ptr;
  struct grow *grow_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  have_grow = FALSE;
  line = 0;

  /* Initialise pointers */
  first_grow_ptr = NULL;
  last_grow_ptr = NULL;

  /* Form the current filename */
  if (work_dir[0] != '\0')
    strcpy(file_name,work_dir);
  else
    strcpy(file_name,pdbsum_dir);
  strcat(file_name,DIR_SEP);
  strcat(file_name,"grow.out");

  /* Open the GROW file */
  if ((file_ptr = fopen(file_name,"r")) !=  NULL)
    {
      /* Loop while reading in records from the file */
      while (fgets(input_line,LINELEN,file_ptr) != NULL)
        {
          /* Check this is not a comment line */
          if (input_line[0] != '#')
            {
              /* Increment line count */
              line++;

              /* Extract the interaction type */
              data_type = input_line[0];

              /* Extract the first residue name */
              strncpy(res_name1,input_line+14,3);
              res_name1[3] = '\0';
              strncpy(res_num1,input_line+9,5);
              res_num1[5] = '\0';
              chain1 = input_line[6];

              /* Extract the second residue name */
              strncpy(res_name2,input_line+44,3);
              res_name2[3] = '\0';
              strncpy(res_num2,input_line+49,5);
              res_num2[5] = '\0';
              chain2 = input_line[54];

              /* Get pointer to the first grow record */
              grow_ptr = first_grow_ptr;
              match = FALSE;

              /* Loop through all the grow records to see if we
                 already have this residue pair */
              while (grow_ptr != NULL && match == FALSE)
                {
                  /* Initialise match flags */
                  match1 = 0;
                  match2 = 0;

                  /* Check for first residue matching either residue
                     from this grow record */
                  if (!strcmp(res_name1,grow_ptr->res_name1) &&
                      !strcmp(res_num1,grow_ptr->res_num1) &&
                      chain1 == grow_ptr->chain1)
                    match1 = 1;
                  else if (!strcmp(res_name1,grow_ptr->res_name2) &&
                           !strcmp(res_num1,grow_ptr->res_num2) &&
                           chain1 == grow_ptr->chain2)
                    match1 = 2;

                  /* Check for second residue matching either residue
                     from this grow record */
                  if (!strcmp(res_name2,grow_ptr->res_name1) &&
                      !strcmp(res_num2,grow_ptr->res_num1) &&
                      chain2 == grow_ptr->chain1)
                    match2 = 1;
                  else if (!strcmp(res_name2,grow_ptr->res_name2) &&
                           !strcmp(res_num2,grow_ptr->res_num2) &&
                           chain2 == grow_ptr->chain2)
                    match2 = 2;

                  /* If both match, then already have the interaction
                     between these two residues */
                  if (match1 > 0 && match2 > 0 && match1 != match2)
                    match = TRUE;

                  /* Get pointer to next grow record */
                  grow_ptr = grow_ptr->next_grow_ptr;
                }

              /* If didn't find a match, create a new grow record */
              if (match == FALSE)
                grow_ptr
                  = create_grow_record(&first_grow_ptr,&last_grow_ptr,
                                       data_type,res_name1,res_num1,
                                       chain1,res_name2,res_num2,chain2);
            }
        }

      /* Close the file */
      fclose(file_ptr);

      /* Set flag if grow data read in */
      if (line > 0)
        have_grow = TRUE;
    }

  /* Return pointer to first grow entry */
  *fst_grow_ptr = first_grow_ptr;

  /* Return whether we have read in grow data */
  return(have_grow);
}
/***********************************************************************

get_chain_col  -  Determine the colour for the given chain

***********************************************************************/

void get_chain_col(char chain,int *icol,char colour_name[COLNAM_LEN],
                      char colour_numbers[COLNO_LEN])
{
  int colour, ichain, nchains;

  static char chain_list[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ123456 7890";
  static char *col_name[] = {
    "purple    ", "red       ", "brown     ", "pink      ", "blue      ",
    "green     ", "cyan      ", "yellow    "
  };
  static char *col_numbers[] = {
    "[178,51,255]", "[255,0,0]", "[128,0,0]", "[255,128,255]", "[0,0,255]",
    "[0,255,0]", "[128,255,255]", "[255,255,0]"
  };

  /* Initialise variable */
  colour = 0;
  colour_name[0] = '\0';
  colour_numbers[0] = '\0';

  /* Get number of chain types */
  nchains = strlen(chain_list);

  /* Get the chain's position in the list */
  for (ichain = 0 ; ichain < nchains && colour == 0; ichain++)
    {
      if (chain == chain_list[ichain])
        {
          /* Store the colour number */
          colour = ichain;

          /* Get the corresponding colour name */
          colour = colour % MAX_OBJ_COL;
          strcpy(colour_name,col_name[colour]);
          strcpy(colour_numbers,col_numbers[colour]);
        }
    }

  /* Return the colour number */
  *icol = colour;
}
/***********************************************************************

get_split_colours  -  Determine the 3 colours for the given chain

***********************************************************************/

int get_split_colours(char chain,char domsst,
                      char colour_name[3][COLNAM_LEN],
                      char colour_numbers[3][COLNO_LEN])
{
  char suffix[3][10];
  int colour, i, ichain, j, nchains;

  double rgb[3][3], swap;

  static char chain_list[]
    = "ABCDEFGHIJKLMNOPQRSTUVWXYZ123456 7890!#_-=:<>|"
    "abcdefghijklmnopqrstuvwxyz";
  static char *col_name[] = {
    "purple", "red", "brown", "pink", "blue", "green", "cyan", "yellow"
  };
  static char *SUFFIX[] = { "_coil", "_helix", "_strand" };

  /* Initialise variables */
  colour = -1;
  for (i = 0; i < 3; i++)
    {
      colour_name[i][0] = '\0';
      colour_numbers[i][0] = '\0';
      strcpy(suffix[i],SUFFIX[i]);
    }

  /* Get number of chain types */
  nchains = strlen(chain_list);

  /* Get the chain's position in the list */
  for (ichain = 0 ; ichain < nchains && colour == -1; ichain++)
    {
      if (chain == chain_list[ichain])
        {
          /* Store the colour number */
          colour = ichain;

          /* Get the corresponding colour name */
          colour = colour % MAX_OBJ_COL;

          /* Get rgb value for this colour */
          get_colour_rgb(col_name[colour],rgb[0]);

          /* Calculate helix and strand colours */
          for (j = 0; j < 3; j++)
            {
              /* Save chain RGB values as helix colours */
              rgb[HELIX][j] = rgb[0][j];

              /* Compute coil colour to be darker */
              rgb[COIL][j] = 0.6 * rgb[HELIX][j];

              /* Compute strand colour to be lighter */
              rgb[STRAND][j] = 0.4 + rgb[HELIX][j];
              if (rgb[STRAND][j] > 1.0)
                rgb[STRAND][j] = 1.0;

              /* If dominant secondary structure is strand, then swap
                 the helix and strand colours */
              if (domsst == 'S')
                {
                  /* Swap helix and strand RGB values */
                  swap = rgb[HELIX][j];
                  rgb[HELIX][j] = rgb[STRAND][j];
                  rgb[STRAND][j] = swap;
                }
            }

          /* Form the RasMol colour numbers for each of these colours */
          for (i = 0; i < 3; i++)
            {
              /* Write colour number string */
              sprintf(colour_numbers[i],"[%d,%d,%d]",
                      (int) (255 * rgb[i][0]),
                      (int) (255 * rgb[i][1]),
                      (int) (255 * rgb[i][2]));

              /* Write out the colour name */
              sprintf(colour_name[i],"%s%s",col_name[colour],suffix[i]);
            }
        }
    }

  /* Return the colour number */
  return(colour);
}
/***********************************************************************

get_domain_colour  -  Get the colour of the given residue, according to
                      its CATH domain assignment

***********************************************************************/

void get_domain_colour(int domain,char colour_name[COLNAM_LEN],
                       char colour_numbers[COLNO_LEN])
{
  int c_numb;

  static char *col_name[] = { 
    "black", "red", "blue", "green", "purple",
    "orange", "cyan", "pink" };
  static char *col_numbers[] = {
    "[0,0,0]", "[255,0,0]", "[0,0,255]", "[0,255,0]", "[179,51,255]",
    "[204,128,0]", "[128,255,255]", "[255,128,255]"
  };

  /* Initialise colour */
  c_numb = 0;

  /* Find the colour corresponding to the given domain */
  if (domain < 0)
    c_numb = 0;
  else if (domain > UPPER_DOMAIN_COL)
    c_numb = domain % UPPER_DOMAIN_COL + 1;
  else
    c_numb = domain;

  /* Get the corresponding colour name */
  strcpy(colour_name,col_name[c_numb]);
  strcpy(colour_numbers,col_numbers[c_numb]);
}
/***********************************************************************

add_to_list  -  Add the current object to list of objects

***********************************************************************/

int add_to_list(char object_list[LINELEN + 1],int list_len,
                char object[OBJ_NAME_LEN])
{

  /* If not first entry, then add comma after last entry prior to
     adding current entry*/
  if (list_len > 0)
    {
      /* Add comma  */
      strcat(object_list,", ");
      list_len = list_len + 2;
    }

  /* Add object's name to list of chains defined */
  if (list_len + 10 < LINELEN)
    {
      /* Add current object and increase string length */
      strcat(object_list,object);
      list_len = list_len + strlen(object);
    }

  /* Return current list length */
  return(list_len);
}
/***********************************************************************

add_residue  -  Add the current residue to the list of residues

***********************************************************************/

void add_residue(char residue_list[LINELEN + 1],char res_num[6])
{
  char residue[6];

  int ipos, first_char_pos;

  /* Initialise variables */
  strncpy(residue,res_num,5);
  residue[5] = '\0';

  /* Strip off the insertion code as RasMol doesn't recognize it at
     present */
  residue[4] = '\0';

  /* Find for non-blank character */
  first_char_pos = -1;
  for (ipos = 0; ipos < 4; ipos++)
    {
      /* If first non-blank then store start-position */
      if (first_char_pos == -1 && residue[ipos] != ' ')
        first_char_pos = ipos;
    }

  /* Add residue to list */
  strcat(residue_list,residue+first_char_pos);
}
/***********************************************************************

terminate_line  -  Terminate line after last non-blank character

***********************************************************************/

void terminate_line(char line[LINELEN + 1],int line_len)
{
  char ch;
  int ipos, last_char_pos, nspaces;
  int done;

  /* Initialise variables */
  done = FALSE;
  last_char_pos = 0;
  nspaces = 0;

  /* Loop over the line */
  for (ipos = 0; ipos < line_len - 1 && done == FALSE; ipos++)
    {
      /* Get the current character */
      ch = line[ipos];

      /* If a blank space, then increment count of blanks */
      if (ch == ' ')
        {
          nspaces++;

          /* If have hit 2 blanks in a row, then can stop */
          if (nspaces > 1)
            done = TRUE;
        }

      /* Check for non-blank character */
      else if (ch != '\n' && ch != '\0')
        {
          /* Save current position */
          last_char_pos = ipos;
          nspaces = 0;
        }

      /* Otherwise have hit the end of the line */
      else
        done = TRUE;
    }

  /* Terminate line after last character position */
  line[last_char_pos + 1] = '\0';
}
/***********************************************************************

get_rasmol_colour  -  Retrieve the corresponding RasMol colour triplet
                      from the supplied colour name

***********************************************************************/

void get_rasmol_colour(char *colour_name,char *colour_numbers)
{

  int icol;
  int match;

  static char *col_name[] = {
    "black", "white", "cream", "paleblue",
    "purple", "yellow", "blue", "red", "green",
    "brown", "pink", "greenblue", "magenta",
    "cyan", "redorange", "violet", "orange",
    "peach", "sky_blue", "grey", "dark_grey"
  };

  static char *col_numbers[] = {
    "[0,0,0]", "[255,255,255]", "[255,255,179]", "[173,216,230]",
    "[178,51,255]", "[255,255,0]", "[0,0,255]", "[255,0,0]", "[0,255,0]",
    "[128,0,0]", "[255,128,255]", "[51,153,102]", "[255,0,255]",
    "[128,255,255]", "[255,102,0]", "[238,130,238]", "[204,128,0]",
    "[255,204,153]", "[0,153,255]", "[128,128,128]", "[51,51,51]"
  };

  /* Initialise variables */
  match = FALSE;
  strcpy(colour_numbers,col_numbers[0]);

  /* Loop until find colour */
  for (icol = 0; icol < NCOLOURS && match == FALSE; icol++)
    {
      /* If have match, get colour number */
      if (!strcmp(colour_name,col_name[icol]))
        {
          /* Transfer colour numbers across */
          strcpy(colour_numbers,col_numbers[icol]);
          match = TRUE;
        }
    }
}
/***********************************************************************

write_script_headers  -  Write all the header records to the output
                         RasMol script file

***********************************************************************/

void write_script_headers(FILE *outfile_ptr,char pdb_code[5],
                          char *pdb_name,char *chains,int nchains,
                          int have_sstdat,char *pdbsum_dir,
                          struct parameters params)
{
  int ichain;

  char colour_name[COLNAM_LEN], colour_numbers[COLNO_LEN];

  /* Write out the starting lines */
  fprintf(outfile_ptr,"#!rasmol -script\n");
  fprintf(outfile_ptr,"# Generated by PDBsum - R A Laskowski, Feb 1998\n");
  fprintf(outfile_ptr,"\n");

  /* Show PDB code and PDB file name */
  fprintf(outfile_ptr,"# PDB code: %s\n",pdb_code);
  fprintf(outfile_ptr,"# PDB file: %s\n",pdb_name);
  fprintf(outfile_ptr,"\n");

  /* Initilise and indicate that ATOM records to follow */
  fprintf(outfile_ptr,"zap\n");
  if (params.pdb_type == PDBSUM1)
    fprintf(outfile_ptr,"load %s\n",pdb_name);
  else
    fprintf(outfile_ptr,"load inline\n");

  /* Set background colour according to where prog has been called from */
  if (params.background_colour[0] != '\0')
    strcpy(colour_name,params.background_colour);
  else
    {
      if (params.show_surfaces == TRUE)
    strcpy(colour_name,"paleblue");
      else if (params.colour_by_domain == TRUE)
    strcpy(colour_name,"peach");
      else
    strcpy(colour_name,"cream");
    }

  /* Get the corresponding RasMol colour triplet */
  get_rasmol_colour(colour_name,colour_numbers);

  /* Write out the script line defining the background colour */
  fprintf(outfile_ptr,"background %s\n",colour_numbers);

  /* Other rendering parameters */
  //fprintf(outfile_ptr,"set ambient 60\n");
  //fprintf(outfile_ptr,"set specular off\n");
  //fprintf(outfile_ptr,"slab off\n");

  /* Apply any supplied rotations */
  if (params.have_transformations == TRUE)
    {
      /* Check for reverse-order rotations */
      if (params.revrot == TRUE)
        {
          /* Write out the rotations */
          fprintf(outfile_ptr,"rotate x %8.2f\n",params.rotation_angle[0]);
          fprintf(outfile_ptr,"rotate y %8.2f\n",params.rotation_angle[1]);
          fprintf(outfile_ptr,"rotate z %8.2f\n",-params.rotation_angle[2]);
        }

      /* Otherwise, have standard order */
      else
        {
          /* Write out the rotations */
          fprintf(outfile_ptr,"rotate z %8.2f\n",-params.rotation_angle[2]);
          fprintf(outfile_ptr,"rotate y %8.2f\n",params.rotation_angle[1]);
          fprintf(outfile_ptr,"rotate x %8.2f\n",params.rotation_angle[0]);
        }
    }

  /* Otherwise, rotate by 180 about x-axis to be in same orientation as
     in PDB file */
  else if (have_sstdat == FALSE)
    fprintf(outfile_ptr,"rotate x 180\n");
  fprintf(outfile_ptr,"\n");

  /* Switch all objects off */
  //fprintf(outfile_ptr,"set bonds off\n");
  //fprintf(outfile_ptr,"set axes off\n");
  //fprintf(outfile_ptr,"set boundingbox off\n");
  //fprintf(outfile_ptr,"set unitcell off\n");
  //fprintf(outfile_ptr,"set bondmode and\n");
  //fprintf(outfile_ptr,"dots off\n");
  fprintf(outfile_ptr,"\n");

  /* Initialise all object colours */
  fprintf(outfile_ptr,"# Initialise colours\n");
  fprintf(outfile_ptr,"select all\n");
  //fprintf(outfile_ptr,"colour bonds none\n");
  //fprintf(outfile_ptr,"colour backbone none\n");
  //fprintf(outfile_ptr,"colour hbonds none\n");
  //fprintf(outfile_ptr,"colour ssbonds none\n");
  //fprintf(outfile_ptr,"colour ribbons none\n");

  /* Switch everything off */
  fprintf(outfile_ptr,"wireframe off\n");
  fprintf(outfile_ptr,"spacefill off\n");
  fprintf(outfile_ptr,"\n");
  fprintf(outfile_ptr,"\n");

  /* Initialise RasMol's x-term output commands */
  fprintf(outfile_ptr,"echo \" \"\n");
  if (params.all_chains == FALSE)
    {
      /* If only one chain, show that */
      if (nchains == 1)
        {
          fprintf(outfile_ptr,"echo \"Structure of %s. Chain [%c]\"\n",
                  pdb_code,chains[0]);
        }

      /* Otherwise, list them all */
      else
        {
          fprintf(outfile_ptr,"echo \"Structure of %s. Chains: ",pdb_code);
          for (ichain = 0; ichain < nchains; ichain++)
            {
              fprintf(outfile_ptr,"%c",chains[ichain]);
              if (ichain < nchains - 1)
                fprintf(outfile_ptr,", ");
            }
          fprintf(outfile_ptr,".\"\n");
        }
    }
  else
    {
      fprintf(outfile_ptr,"echo \"Structure of %s.\"\n",pdb_code);
    }

  /* If colouring by domain, add note */
  if (params.colour_by_domain == TRUE)
    fprintf(outfile_ptr,"echo \"   Coloured by domain\"\n");

  /* Extra line-throws */
  fprintf(outfile_ptr,"echo \" \"\n");
  fprintf(outfile_ptr,"\n");
  fprintf(outfile_ptr,"\n");
}
/***********************************************************************

write_script_site  -  Write script records to display active site
                      residues

***********************************************************************/

void write_script_site(FILE *outfile_ptr,char *site_name,
                       char *input_line,char *object,char *display_mode,
                       struct parameters params)
{

  /* Chop the input line after the last character */
  terminate_line(input_line,LINELEN - 5);

  /* Write comment line */
  fprintf(outfile_ptr,"# Active site: %s\n",site_name);

  /* Define the site object */
  strcpy(object,"site");
  strcat(object,site_name);
  display_mode[0] = '\0';

  /* Write out the site definition */
  fprintf(outfile_ptr,"define %s %s\n",object,input_line);

  /* If displaying active sites, then issue script commands */
  if (params.show_sites)
    {
      /* Show just the side-chain and C-alpha atoms as thick sticks */
      fprintf(outfile_ptr,"select %s & (sidechain | alpha) ",object);
      fprintf(outfile_ptr,"& !hydrogen\n");
      strcpy(display_mode,"wireframe");
      fprintf(outfile_ptr,"%s 0.3\n",display_mode);

      /* Colour just the sidechain atoms as CPK */
      fprintf(outfile_ptr,"select %s & sidechain\n",object);
      fprintf(outfile_ptr,"colour atoms cpk\n");
    }

  fprintf(outfile_ptr,"\n");
}
/***********************************************************************

write_script_protein  -  Write script records defining current protein
                         chain

***********************************************************************/

void write_script_protein(FILE *outfile_ptr,char chain,int ca_only,
                          char object[OBJ_NAME_LEN],
                          char display_mode[DISPLAY_MODE_LEN],
                          char chains_list[LINELEN + 1],int *list_len,
                          char domsst,struct parameters params)
{
  char add_on[11];
  char colour_name[3][COLNAM_LEN], colour_numbers[3][COLNO_LEN];
  char copy_line[LINELEN];

  int colour, factor, i, ichain, lower_b, upper_b;
  int wanted;

  /* Initialise strings */
  copy_line[0] = '\0';
  display_mode[0] = '\0';
  object[0] = '\0';
  add_on[0] = '\0';
  if (ca_only == TRUE)
    strcpy(add_on," (CA-only)");
  if (params.jmol == TRUE)
    factor = 1;
  else
    factor = 100;

  /* If chain-id is blank, then use whole protein */
  if (chain == ' ')
    {
      /* Write comment line and define protein object */
      fprintf(outfile_ptr,"# Protein%s\n",add_on);
      strcpy(object,"protein");

      /* Select the protein */
      if (params.show_protein == TRUE)
        fprintf(outfile_ptr,"select protein\n");
      sprintf(copy_line,"select protein");
    }

  /* Otherwise, consider the given protein chain */
  else
    {
      /* Write comment line */
      fprintf(outfile_ptr,"# Chain %c - protein%s\n",chain,add_on);
      sprintf(object,"chn%c",chain);

      /* Define this chain and then select it */
      fprintf(outfile_ptr,"define %s *:%c\n",object,chain);
      if (params.show_protein == TRUE)
        fprintf(outfile_ptr,"select %s & protein\n",object);
      sprintf(copy_line,"select %s & protein",object);
    }

  /* If protein to be shown, then display */
  if (params.show_protein == TRUE)
    {
      /* Get this chain's colour */
      get_chain_col(chain,&colour,colour_name[0],colour_numbers[0]);

      /* Set chain colour */
      fprintf(outfile_ptr,"wireframe off\n");
      if (params.colour_option == COLOUR_BY_SST)
        fprintf(outfile_ptr,"colour structure\n");

      /* If this is Alpha Fold model, colour by confidence score */
      else if (params.pdb_type == ALPHA_FOLD_PDB)
        {
          /* Loop over the confidence ranges */
          for (i = 0; i < NAF_RANGES; i++)
            {
              /* Get the B-factor range */
              if (i == 0)
                lower_b = 0;
              else
                lower_b = factor * AF_SCORE[i - 1];
              upper_b = factor * AF_SCORE[i];

              /* Colour this range */
              fprintf(outfile_ptr,"select temperature > %d and "
                      "temperature <= %d\n",lower_b,upper_b);
              fprintf(outfile_ptr,"colour [%d,%d,%d]\n",
                      AF_COLOUR_IRGB[i][0],
                      AF_COLOUR_IRGB[i][1],
                      AF_COLOUR_IRGB[i][2]);
            }

          /* Select protein again */
          fprintf(outfile_ptr,"%s\n",copy_line);
        }

      /* If have a dominant secondary structure, adjust colours by
         secondary structure type */
      else if (domsst != ' ')
        {
          /* Get the split colours */
          get_split_colours(chain,domsst,colour_name,colour_numbers);

          /* Colour by secondary structure type */
          fprintf(outfile_ptr,"colour %s  # (%s)\n",
                  colour_numbers[COIL],colour_name[COIL]);
          fprintf(outfile_ptr,"select %s & protein & helix\n",object);
          fprintf(outfile_ptr,"colour %s  # (%s)\n",
                  colour_numbers[HELIX],colour_name[HELIX]);
          fprintf(outfile_ptr,"select %s & protein & sheet\n",object);
          fprintf(outfile_ptr,"colour %s  # (%s)\n",
                  colour_numbers[STRAND],colour_name[STRAND]);
          fprintf(outfile_ptr,"select %s & protein\n",object);
        }

      /* Otherwise, colour chain by single colour */
      else
        fprintf(outfile_ptr,"colour %s  # (%s)\n",colour_numbers[0],
                colour_name[0]);

      /* For full protein display as cartoon */
      if (ca_only == FALSE)
        {
          fprintf(outfile_ptr,"cartoon on\n");
          strcpy(display_mode,"cartoon");
        }

      /* Otherwise, for CA-only protein chain, display as backbone */
      else
        {
          fprintf(outfile_ptr,"backbone on\n");
          strcpy(display_mode,"backbone");
        }
    }
  fprintf(outfile_ptr,"\n");

  /* Add current chain to the list of chains defined */
  *list_len = add_to_list(chains_list,*list_len,object);
}
/***********************************************************************

write_script_ligand  -  Write script records to display ligands

***********************************************************************/

void write_script_ligand(FILE *outfile_ptr,char ligand_num[3],
                         char input_line[LINELEN + 1],
                         char object[OBJ_NAME_LEN],
                         char display_mode[DISPLAY_MODE_LEN],
                         struct parameters params)
{
  char colour_numbers[COLNO_LEN];

  /* Chop the input line after the last character */
  terminate_line(input_line,LINELEN - 4);

  /* Write comment line */
  fprintf(outfile_ptr,"# Ligand %s\n",ligand_num);

  /* Define the ligand object */
  strcpy(object,"lig");
  strcat(object,ligand_num);
  display_mode[0] = '\0';

  /* Write out the ligand definition */
  fprintf(outfile_ptr,"define %s %s\n",object,input_line);

  /* If ligands to be shown, then display */
  if (params.show_ligands == TRUE)
    {
      /* Show ligand atoms (excluding hydrogens) */
      fprintf(outfile_ptr,"select %s & !hydrogen\n",object);

      /* Show bonds as sticks */
      sprintf(display_mode,"wireframe %5.2f",params.ligbond_radius);

      /* Write out the wireframe command */
      fprintf(outfile_ptr,"%s\n",display_mode);

      /* Set atom radius */
      if (params.ligatom_radius > 0)
        fprintf(outfile_ptr,"spacefill %5.2f\n",params.ligatom_radius);

      /* If radius is negative, use van der Waals atom radii */
      else if (params.ligatom_radius < 0)
        fprintf(outfile_ptr,"spacefill on\n");

      /* Colour as CPK and show in spacefill mode */
      if (!strcmp(params.ligand_colour,"CPK"))
        {
          fprintf(outfile_ptr,"colour atoms cpk\n");

          /* Colour the carbon atoms a dark grey */
          fprintf(outfile_ptr,"select %s & carbon\n",object);
          fprintf(outfile_ptr,"colour [56,56,70]\n");
        }
      else
        {
          /* Get RGB values for ligand colour */
          get_rasmol_colour(params.ligand_colour,colour_numbers);

          /* Write out the colour */
          fprintf(outfile_ptr,"colour atoms %s\n",colour_numbers);
        }
    }
  fprintf(outfile_ptr,"\n");
}
/***********************************************************************

write_script_metal  -  Write script records to define and display metal
                       ions

***********************************************************************/

void write_script_metal(FILE *outfile_ptr,char *metal_name,
                        char *metal_colour,char *object,
                        char *display_mode,char *metals_list,
                        int *list_len,struct parameters params)
{
  char full_name[8], tmp_name[5];
  int ipos;
  float radius;

  /* Write comment line */
  fprintf(outfile_ptr,"# Metal: %s\n",metal_name);

  /* If metal name starts with a space, shift it along */
  if (metal_name[0] == ' ')
    {
      strcpy(tmp_name,metal_name + 1);
      strcpy(metal_name,tmp_name);
    }

  /* Truncate metal's name where it ends */
  for (ipos = 0; ipos < 4; ipos++)
    if (metal_name[ipos] == ' ')
      metal_name[ipos] = '\0';

  /* Truncate metal's colour where it ends */
  for (ipos = 0; ipos < COLNAM_LEN - 1; ipos++)
    if (metal_colour[ipos] == ' ' || metal_colour[ipos] == '\n')
      metal_colour[ipos] = '\0';

  /* Form the metal name */
  strcpy(full_name,"*.");
  strcat(full_name,metal_name);

  /* Deal with special cases: Ca and Cd */
  if (!strncmp(metal_name,"CA",2))
    strcpy(full_name,"CA.CA");
  if (!strncmp(metal_name,"CD",2))
    strcpy(full_name,"CD");

  /* Special case of single-letter metals such as K, U, etc */
  if (metal_name[0] == '_')
    {
      full_name[0] = metal_name[1];
      full_name[1] = '\0';
    }

  /* Set radius */
  radius = params.metatom_radius;

  /* Form object name for this metal type */
  strcpy(object,"metal");
  strcat(object,metal_name);
  display_mode[0] = '\0';

  /* Define metal */
  fprintf(outfile_ptr,"define %s %s\n",object,full_name);

  /* If metals to be shown, then display */
  if (params.show_metals == TRUE)
    {
      /* Select and colour the metal */
      fprintf(outfile_ptr,"select %s & !hydrogen\n",object);
      fprintf(outfile_ptr,"colour atoms %s\n",metal_colour);
      strcpy(display_mode,"spacefill");
      fprintf(outfile_ptr,"%s %8.3f\n",display_mode,radius);
    }
  fprintf(outfile_ptr,"\n");

  /* Add current chain to the list of chains */
  *list_len = add_to_list(metals_list,*list_len,object);
}
/***********************************************************************

write_script_nucleic  -  Write script records defining current nucleic
                         acid chain

***********************************************************************/

void write_script_nucleic(FILE *outfile_ptr,char chain,
                          char object[OBJ_NAME_LEN],
                          char display_mode[DISPLAY_MODE_LEN],
                          char chains_list[LINELEN + 1],int *list_len,
                          struct parameters params)
{
  char colour_name[COLNAM_LEN], colour_numbers[COLNO_LEN];
  int colour;

  /* Initialise strings */
  display_mode[0] = '\0';
  object[0] = '\0';

  /* If chain-id is blank, then use whole object */
  if (chain == ' ')
    {
      /* Write comment line and define nucleic acid object */
      fprintf(outfile_ptr,"# Nucleic acid chain\n");
      strcpy(object,"nucleic");
    }

  /* Otherwise, consider the given protein chain */
  else
    {
      /* Write comment line */
      fprintf(outfile_ptr,"# Chain %c - nucleic acid\n",chain);
      sprintf(object,"chn%c",chain);

      /* Define this chain */
      fprintf(outfile_ptr,"define %s *:%c\n",object,chain);
    }


  /* If nucleic acid chain to be shown, then display */
  if (params.show_nucleic == TRUE)
    {
      /* Select the given object */
      fprintf(outfile_ptr,"select %s & !hydrogen\n",object);

      /* Colour by atom-colour and show as thin sticks */
      fprintf(outfile_ptr,"colour atoms cpk\n");
      fprintf(outfile_ptr,"wireframe 0.2\n");
      strcpy(display_mode,"wireframe");

      /* Get chain-colour for colouring DNA backbone */
      get_chain_col(chain,&colour,colour_name,colour_numbers);

      /* Select backbone and colour accordingly */
      fprintf(outfile_ptr,"select %s & *.P\n",object);
      fprintf(outfile_ptr,"colour atoms %s  # (%s)\n",colour_numbers,
              colour_name);
      fprintf(outfile_ptr,"backbone on\n");
    }
  fprintf(outfile_ptr,"\n");

  /* Add current chain to the list of chains */
  *list_len = add_to_list(chains_list,*list_len,object);
}
/***********************************************************************

write_script_dots  -  Define colours for this chain's dot surface

***********************************************************************/

void write_script_dots(FILE *outfile_ptr,char chain,int atom_from,
                       int atom_to,
                       char dot_surfaces_list[LINELEN + 1])
{
  char colour_name[COLNAM_LEN], colour_numbers[COLNO_LEN];
  char object[OBJ_NAME_LEN];
  int colour, line_len;

  /* Write comment line */
  fprintf(outfile_ptr,"# Dot surface for chain %c\n",chain);

  /* Define this dot surface */
  sprintf(object,"surf%c",chain);

  fprintf(outfile_ptr,"define %s atomno >= %d & atomno <= %d\n",
          object,atom_from,atom_to);

  /* Get the colour for this surface according to its chain */
  get_chain_col(chain,&colour,colour_name,colour_numbers);

  /* Define the chain colour */
  fprintf(outfile_ptr,"select %s\n",object);
  fprintf(outfile_ptr,"colour atoms %s  # (%s)\n",colour_numbers,
          colour_name);
  fprintf(outfile_ptr,"\n");

  /* Add current dot surface to the list */

  /* If list empty, then put current dot surface at start */
  if (dot_surfaces_list[0] == '\0')
    strcpy(dot_surfaces_list,object);

  /* Otherwise, check that list not too long */
  else
    {
      /* Get length of list */
      line_len = strlen(dot_surfaces_list);

      /* If not too long, then add current dot surface to it */
      if (line_len < LINELEN - 9)
        {
          strcat(dot_surfaces_list," | ");
          strcat(dot_surfaces_list,object);
        }
    }
}
/***********************************************************************

write_script_surface  -  Write script records to define and display
                         protein surface

***********************************************************************/

void write_script_surface(FILE *outfile_ptr,char chain,int surface,
                          char surfaces_list[LINELEN + 1],int *list_len,
                          struct parameters params)
{
  char colour_name[COLNAM_LEN], colour_numbers[COLNO_LEN];
  char object[OBJ_NAME_LEN];
  int colour;

  /* Write comment line */
  fprintf(outfile_ptr,"# Surface %d - chain %c\n",surface,chain);

  /* Define this surface */
  sprintf(object,"surf%c",chain);

  /* If surfaces to be shown, then display */
  if (params.show_surfaces == TRUE)
    {
      fprintf(outfile_ptr,"define %s MAP%d:%d\n",object,surface,surface);

      /* Get the colour for this surface according to its chain */
      get_chain_col(chain,&colour,colour_name,colour_numbers);

      /* Define the chain colour */
      fprintf(outfile_ptr,"select %s\n",object);
      fprintf(outfile_ptr,"wireframe on\n");
      fprintf(outfile_ptr,"colour bonds %s  # (%s)\n",colour_numbers,
              colour_name);
    }
  fprintf(outfile_ptr,"\n");

  /* Add current surfaces to the list of surfaces */
  *list_len = add_to_list(surfaces_list,*list_len,object);
}
/***********************************************************************

write_script_mask  -  Write script records to define and display
                      the binding-site mask

***********************************************************************/

void write_script_mask(FILE *outfile_ptr,int surface,
                       struct parameters params)
{
  char colour_name[COLNAM_LEN];
  char object[OBJ_NAME_LEN];

  /* Write comment line */
  fprintf(outfile_ptr,"# Surface %d - mask\n",surface);

  /* Define this surface */
  sprintf(object,"mask");
  fprintf(outfile_ptr,"define %s MAP%d:%d\n",object,surface,surface);

  /* If surfaces to be shown, then display */
  if (params.show_mask == TRUE)
    {
      /* Set mask colour */
      strcpy(colour_name,"green");

      /* Define the chain colour */
      fprintf(outfile_ptr,"select %s\n",object);
      fprintf(outfile_ptr,"wireframe on\n");
      fprintf(outfile_ptr,"colour bonds %s\n",colour_name);
    }
  fprintf(outfile_ptr,"\n");
}
/***********************************************************************

write_script_sets  -  Write out the script records that will display the
                      user-defined sets in the file

***********************************************************************/

void write_script_sets(FILE *outfile_ptr,
                       char ligands_list[MAXLIGANDS][MAX_RANGE_LEN],
                       int nligands,
                       char mchains_list[LINELEN + 1],int mlist_len,
                       char nchains_list[LINELEN + 1],int nlist_len,
                       char pchains_list[LINELEN + 1],int plist_len,
                       char sites_list[MAXSITE][4],int nsites,
                       char surfaces_list[LINELEN + 1],int slist_len,
                       char dot_surfaces_list[LINELEN + 1],
                       int have_backbone,int have_cartoon,int have_dots,
                       int have_mask,struct parameters params)
{
  char ligand_num[3], display_mode[DISPLAY_MODE_LEN];
  int ligand, site;


  /* If have dot surfaces, then switch them on */
  if (have_dots == TRUE)
    {
      fprintf(outfile_ptr,"# Switch dot surfaces on\n");
      fprintf(outfile_ptr,"define dot_surfaces %s\n",dot_surfaces_list);
      fprintf(outfile_ptr,"select dot_surfaces\n");
      if (params.show_dots == TRUE)
        fprintf(outfile_ptr,"dots on\n");
      fprintf(outfile_ptr,"\n");
    }

  /* Output start of list of user-defined objects in this RasMol script */
  fprintf(outfile_ptr,"echo \" \"\n");
  fprintf(outfile_ptr,"echo \"* User-defined sets present in structure:-\"\n");
  fprintf(outfile_ptr,"echo \" \"\n");

  /* Write out the list of protein chains, if any */
  if (plist_len > 0)
    {
      /* Determine how protein chains are shown */
      if (have_backbone == TRUE)
        strcpy(display_mode,"backbone");
      if (have_cartoon == TRUE)
        strcpy(display_mode,"cartoon");
      if (have_backbone == TRUE && have_cartoon == TRUE)
        strcpy(display_mode,"cartoon & backbone");

      /* Say how protein chains are shown */
      fprintf(outfile_ptr,"echo \"    protein = All protein atoms   ");
      fprintf(outfile_ptr,"[shown as: %s]\"\n",display_mode);

      /* If list of chains is long, then truncate it */
      if (plist_len > 15)
        {
          strncpy(pchains_list+15," ...",4);
          pchains_list[19] = '\0';
        }

      /* Show the list of protein chains */
      if (!strncmp(pchains_list,"chn",3))
        {
          fprintf(outfile_ptr,"echo \"    %s = separate protein chains",
                  pchains_list);
          fprintf(outfile_ptr,"   [shown as: %s]\"\n",display_mode);
        }
    }

  /* Write out the list of nucleic acid chains, if any */
  if (nlist_len > 0)
    {
      /* Say how nucleic acid chains are shown */
      strcpy(display_mode,"wireframe 0.2");
      fprintf(outfile_ptr,"echo \"    nucleic = All nucleic acid atoms   ");
      fprintf(outfile_ptr,"[shown as: %s]\"\n",display_mode);

      /* If list of chains is long, then truncate it */
      if (nlist_len > 15)
        {
          strncpy(nchains_list+15," ...",4);
          nchains_list[19] = '\0';
        }

      /* Show the list of nucleic acid chains */
      if (!strncmp(nchains_list,"chain",5))
        {
          fprintf(outfile_ptr,"echo \"    %s = separate nucleic acid chains",
                  nchains_list);
          fprintf(outfile_ptr,"   [shown as: %s]\"\n",display_mode);
        }
    }

  /* Write out the list of metals, if any */
  if (mlist_len > 0)
    {
      /* Say how metals are shown */
      strcpy(display_mode,"spacefill");

      /* Show the list of metals */
      fprintf(outfile_ptr,"echo \"    %s = metal ion(s)",
              mchains_list);
      fprintf(outfile_ptr,"   [shown as: %s]\"\n",display_mode);
    }

  /* Write out the list of ligand ranges, if any */
  for (ligand = 0; ligand < nligands; ligand++)
    {
      /* Say how this ligand is shown */
      strcpy(display_mode,"spacefill");

      /* Form the ligand number */
      sprintf(ligand_num,"%2d",ligand + 1);
      if (ligand_num[0] == ' ')
        ligand_num[0] = '0';

      /* Show the ligand definition */
      fprintf(outfile_ptr,"echo \"    lig%s = ligand residue",ligand_num);

      /* If more than one residue, print "s" */
      if (strstr(ligands_list[ligand]," to ") != NULL)
        fprintf(outfile_ptr,"s");

      /* Show residue range */
      fprintf(outfile_ptr," %s",ligands_list[ligand]);
      fprintf(outfile_ptr," [shown as: %s]\"\n",display_mode);
    }

  /* Write out the list of active sites, if any */
  for (site = 0; site < nsites; site++)
    {
      /* Say how this active site is shown */
      strcpy(display_mode,"wireframe 0.3");

      /* Show the active site */
      fprintf(outfile_ptr,"echo \"    site%s = active site %s ",
              sites_list[site],sites_list[site]);
      fprintf(outfile_ptr,"(from PDB file)   [shown as: %s]\"\n",
              display_mode);
    }

  /* Write out the list of surfaces, if any */
  if (slist_len > 0)
    {
      /* If list of surfaces is long, then truncate it */
      if (slist_len > 12)
        {
          strncpy(surfaces_list+12," ...",4);
          surfaces_list[16] = '\0';
        }

      /* Show the list of surfaces */
      fprintf(outfile_ptr,"echo \"    %s = van der Waals surfaces",
              surfaces_list);
      fprintf(outfile_ptr,"   [shown as: wireframe]\"\n");
    }

  /* If have a mask surface, then show that */
  if (have_mask == TRUE)
    {
      /* Show mask definition */
      fprintf(outfile_ptr,"echo \"    mask = wire-cage outline of binding");
      fprintf(outfile_ptr," site   [shown as: green wireframe]\"\n");
    }

  /* If have dot surfaces, then show how to switch on and off */
  if (have_dots == TRUE)
    {
      /* Show mask definition */
      fprintf(outfile_ptr,"echo \"    Dot surfaces can be switched on and");
      fprintf(outfile_ptr," off with: 'dots on' and 'dots off'\"\n");
    }

  fprintf(outfile_ptr,"\n");
}
/***********************************************************************

write_script_tail  -  Write all the final records to the output RasMol
                      script file

***********************************************************************/

void write_script_tail(FILE *outfile_ptr,char object[OBJ_NAME_LEN],
                       char display_mode[DISPLAY_MODE_LEN],
                       struct parameters params)
{

  /* Define all disulphides */
  if (params.show_disulphides == TRUE)
    {
      fprintf(outfile_ptr,"# Disulphide bonds\n");
      fprintf(outfile_ptr,"select cystine & (sidechain | alpha) &");
      fprintf(outfile_ptr,"!hydrogen\n");
      fprintf(outfile_ptr,"wireframe 0.3\n");
      fprintf(outfile_ptr,"select cystine & sidechain\n");
      fprintf(outfile_ptr,"colour atoms cpk\n");
      fprintf(outfile_ptr,"colour ssbonds [128,0,0]\n");
      fprintf(outfile_ptr,"ssbonds on\n");
      fprintf(outfile_ptr,"\n");
    }

  /* Switch all waters off */
  fprintf(outfile_ptr,"# Waters off\n");
  fprintf(outfile_ptr,"select HOH\n");
  fprintf(outfile_ptr,"wireframe off\n");
  fprintf(outfile_ptr,"\n");

  /* Show some examples of how the user-definable sets might be used */
  fprintf(outfile_ptr,"# Example of use\n");
  fprintf(outfile_ptr,"echo \" \"\n");
  fprintf(outfile_ptr,"echo \"* Examples of use:-\"\n");
  fprintf(outfile_ptr,"echo \" \"\n");
  fprintf(outfile_ptr,"echo \"    To switch %s off:-\"\n",object);
  fprintf(outfile_ptr,"echo \"          > select %s\"\n",object);
  fprintf(outfile_ptr,"echo \"          > %s off\"\n",display_mode);

  /* If display mode is wireframe, then show alternative as backbone */
  if (!strncmp(display_mode,"wireframe",9))
    {
      fprintf(outfile_ptr,"echo \"    Then, to display as backbone:-\"\n");
      fprintf(outfile_ptr,"echo \"          > backbone on\"\n");
      fprintf(outfile_ptr,"echo \" \"\n");
    }

  /* If display mode is backbone, then show alternative as spacefill */
  if (!strncmp(display_mode,"backbone",8))
    {
      fprintf(outfile_ptr,"echo \"    Then, to display as spacefill:-\"\n");
      fprintf(outfile_ptr,"echo \"          > spacefill on\"\n");
      fprintf(outfile_ptr,"echo \"      or, say\"\n");
      fprintf(outfile_ptr,"echo \"          > spacefill %4.2f\"\n",
              params.metatom_radius);
      fprintf(outfile_ptr,"echo \" \"\n");
    }

  /* Otherwise, show both wireframe options as alternatives */
  else
    {
      fprintf(outfile_ptr,"echo \"    Then, to display as wireframe:-\"\n");
      fprintf(outfile_ptr,"echo \"          > wireframe on\"\n");
      fprintf(outfile_ptr,"echo \"      or, say\"\n");
      fprintf(outfile_ptr,"echo \"          > wireframe 0.3\"\n");
      fprintf(outfile_ptr,"echo \" \"\n");
    }

  /* If generating a GIF image, then create it here */
  if (params.write_gif_only == TRUE)
    fprintf(outfile_ptr,"write gif rascript.gif\n");

  /* Close off script */
  fprintf(outfile_ptr,"\n");
  fprintf(outfile_ptr,"select all\n");
  fprintf(outfile_ptr,"exit\n");
  fprintf(outfile_ptr,"\n");
}
/***********************************************************************

get_atoken - Get the next space- or tab-delimited token from the given
             input line, starting at the given position

**********************************************************************/  

int get_atoken(char line[LINELEN +1],int *ipos,char token[TOKEN_LEN + 1],
           int *done)
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
  while (*done == FALSE && have_token == FALSE && *ipos < LINELEN)
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
  token_len = i;
  if (i > 0)
    *done = FALSE;

  /* Return the token */
  return(token_len);
} 
/***********************************************************************

extract_domain_data - Extract the domain data from the free-format line
                      of the CATH domains file and store in a fixed-
                      format array

**********************************************************************/  

int extract_domain_data(char input_line[LINELEN + 1],
                        char domain_details[C_MAXDOMCHAIN][C_DOMLEN],
                        int *ndomains,int *nfrags,
                        char error_message[LINELEN + 1])
{
  char chain, number_string[20], res_num[6], token_string[TOKEN_LEN + 1];

  int idomain, ifrag, isegment, lpos, ndetails, nsegments;
  int len, start_pos, token_len;
  int done, have_hyphen, next_is_frag, state;

  /* Initialise variables */
  chain = input_line[4];
  error_message[0] = '\0';
  idomain = 0;
  ifrag = 0;
  ndetails = 0;
  *ndomains = 0;
  *nfrags = 0;
  lpos = 6;
  state = NDOMAINS;
  token_len = 0;
  token_string[0] = '\0';
  done = FALSE;

  /* Loop through line, character by character, until end reached */
  while (done == FALSE)
    {
      /* Get next token string */
      start_pos = lpos;
      token_len = get_atoken(input_line,&lpos,token_string,&done);

      /* Process the token-string according to what it is */

      /* If token is empty, then skip */
      if (token_len == 0 || (token_len > 1 && token_string[0] == '('))
        continue;

      /* Number of domains */
      else if (state == NDOMAINS)
        {
          /* Get the number of domains */
          if (token_string[0] == 'D')
            {
              /* Extract the number */
              strcpy(number_string,token_string+1);
              *ndomains = atoi(number_string);
            }

          /* If have wrong token, then generate error message */
          else
            {
              sprintf(error_message,"Missing Dxx token. Instead have [%s]",
                      token_string);
              done = TRUE;
            }

          /* Switch state */
          state = NFRAGS;
        }

      /* Number of fragments */
      else if (state == NFRAGS)
        {
          /* Get the number of inter-domain fragments */
          if (token_string[0] == 'F')
            {
              /* Extract the number */
              strcpy(number_string,token_string+1);
              *nfrags = atoi(number_string);
            }

          /* If have wrong token, then generate error message */
          else
            {
              sprintf(error_message,"Missing Fxx token. Instead have [%s]",
                      token_string);
              done = TRUE;
            }

          /* Switch state */
          state = NSEGMENTS;
        }

      /* Number of segments */
      else if (state == NSEGMENTS)
        {
          /* Get the number of segments */
          strcpy(number_string,token_string);
          nsegments = atoi(number_string);

          /* Initialise the details line */
          sprintf(domain_details[ndetails],"D%2d",idomain);
          isegment = 0;

          /* Switch state */
          state = CHAIN1;
        }

      /* Chain identifier */
      else if (state == CHAIN1 || state == CHAIN2)
        {
          /* If too long to be a chain-id, then show error message */
          if (token_len > 1)
            {
              sprintf(error_message,"Expected chain-id is too long [%s]",
                      token_string);
              done = TRUE;
            }

          /* Otherwise, store the chain-id */
          else 
            strcat(domain_details[ndetails],token_string);

          /* Check whether the chain agrees with that at the start of
             the record */
          if (token_string[0] != chain)
            {
              sprintf(error_message,"Chain mismatch [%c] != [%c]",
                      token_string[0],chain);
              done = TRUE;
            }

          /* Switch state */
          if (state == CHAIN1)
            state = RES1;
          else
            state = RES2;
        }

      /* Residue number */
      else if (state == RES1 || state == RES2)
        {
          /* Check whether have a hyphen on the end of the residue number */
          if (token_string[token_len - 1] == '-')
            {
              /* Have hyphen, so strip it off */
              have_hyphen = TRUE;
              token_len = token_len - 1;
              token_string[token_len] = '\0';
            }
          else
            have_hyphen = FALSE;

          /* Format the residue number as a 5-character string */
          format_residue_number(res_num,token_string);

          /* Store the residue number */
          strcat(domain_details[ndetails],res_num);

          /* Switch state */
          if (state == RES1)
            {
              if (have_hyphen == TRUE)
                state = CHAIN2;
              else
                state = HYPHEN1;
            }
          else
            {
              if (have_hyphen == TRUE)
                {
                  strcpy(token_string,"-");
                  token_len = 1;
                  state = DONE_SEG;
                }
              else
                state = HYPHEN2;
            }
        }

      /* Hyphen/insertion-code between residues */
      else if (state == HYPHEN1)
        {
          /* Check that is really a hyphen, then switch state */
          if (token_string[0] == '-')
            state = CHAIN2;

          /* Otherwise, take this to be the insertion code */
          else if (res_num[4] == ' ')
            {
              domain_details[ndetails][8] = token_string[0];
              state = CHAIN2;
            }

          /* Otherwise, if already have an insertion code, then there's
             something wrong */
          else
            {
              sprintf(error_message,
                      "Duplicate insertion code on 1st residue [%s] - [%s]",
                      res_num,token_string);
              done = TRUE;
            }
        }

      /* Hyphen/insertion-code after second residue */
      else if (state == HYPHEN2)
        {
          /* Check that is really a hyphen, then switch state */
          if (token_string[0] == '-')
            state = DONE_SEG;

          /* Otherwise, take this to be the insertion code */
          else if (res_num[4] == ' ')
            {
              domain_details[ndetails][14] = token_string[0];
              state = DONE_SEG;
            }

          /* Otherwise, if already have an insertion code, then there's
             something wrong */
          else
            {
              sprintf(error_message,
                      "Duplicate insertion code on 2nd residue [%s] - [%s]",
                      res_num,token_string);
              done = TRUE;
            }
        }

      /* If details complete, prepare for the next */
      if (state == DONE_SEG)
        {
          /* Check that length exactly 15 characters */
          len = strlen(domain_details[ndetails]);
          if (len == 15)
            {
              /* Determine what's coming next */

              /* Increment segment count */
              isegment++;

              /* If have finished segments for this domain/fragment
                 then expect number of segments next */
              if (isegment >= nsegments)
                {
                  /* Initialise flag */
                  next_is_frag = FALSE;

                  /* If domain just completed, the increment count */
                  if (domain_details[ndetails][0] == 'D')
                    {
                      /* Increment domains count */
                      idomain++;

                      /* If all domains just done, next info will be
                         about the fragments */
                      if (idomain >= *ndomains)
                        next_is_frag = TRUE;

                      /* Otherwise, expecting number of segments making
                         up the domain next */
                      else
                        state = NSEGMENTS;
                    }

                  /* If domain just completed was a fragment, increment
                     count */
                  else if (domain_details[ndetails][0] == 'F')
                    {
                      /* Increment fragment count and set flag */
                      ifrag++;
                      next_is_frag = TRUE;

                      /* If all fragments done, then can stop here */
                      if (ifrag >= *nfrags)
                        {
                          done = TRUE;
                          next_is_frag = FALSE;
                        }
                    }

                  /* For a fragment, initialise the details */
                  if (next_is_frag == TRUE)
                    {
                      /* Start the next details line */
                      if (ndetails < C_MAXDOMCHAIN - 1)
                        sprintf(domain_details[ndetails + 1],"F%2d",ifrag);

                      /* Set state to look for a chain-id */
                      state = CHAIN1;
                      isegment = 0;
                    }
                }

              /* Otherwise, expect chain-id of first residue next */
              else
                {
                  /* Copy the domain details for next residue range */
                  if (ndetails < C_MAXDOMCHAIN - 1)
                    {
                      strncpy(domain_details[ndetails + 1],
                              domain_details[ndetails],3);
                      domain_details[ndetails + 1][3] = '\0';
                    }

                  /* Set state to look for chain-id */
                  state = CHAIN1;
                }

              /* Increment count of details stored */
              if (ndetails < C_MAXDOMCHAIN - 1)
                ndetails++;
              else
                {
                  sprintf(error_message,
                          "Maximum number of domain details exceeded %d",
                          ndetails);
                  done = TRUE;
                }
            }

          /* Otherwise, give an error message */
          else
            {
              sprintf(error_message,
                      "Domain details wrong length [%s] from pstn %d",
                      token_string,start_pos);
              done = TRUE;
            }
        }
    }

  /* Return numbers of domain segments stored */
  return(ndetails);
}  
/***********************************************************************

check_chain  -  Check whether given chain is one of those wanted

***********************************************************************/

int check_chain(char chain,char *chains,int nchains,int *chn_no)
{
  int ichain;
  int wanted;

  /* Initialise flag */
  *chn_no = -1;
  wanted = FALSE;

  /* Loop through wanted chains to see if this is one of them */
  for (ichain = 0; ichain < nchains && wanted == FALSE; ichain++)
    {
      /* If have a match, then set flag */
      if (chain == chains[ichain])
        wanted = TRUE;

      /* Allow for CATH assignment of blank chains to chain '0' */
      else if (chain == '0' && (chains[ichain] == ' ' ||
                                chains[ichain] == 'A'))
        wanted = TRUE;

      /* Store chain's position */
      *chn_no = ichain;
    }

  /* Return whether this chain is wanted */
  return(wanted);
}
/***********************************************************************

read_in_domain_data  -  Read in the domain data on contacts with ligands
                        for the given structure, if present

***********************************************************************/

int read_in_domain_data(char pdb_code[5],char *chains,int nchains,
                        char cath_domains_file[FILENAME_LEN],
                        char domain_details[C_MAXDOMCHAIN][C_DOMLEN],
                        char error_message[LINELEN + 1])
{
  char chain_id;
  char cath_code[C_CATH_LEN], wolf_code[7];
  char input_line[LINELEN + 1];

  int ichain, line, line_length;
  int ndetails, ndomains, nfrags;
  int done, got_rec, wanted;

  FILE *file_ptr;

  /* Initialise variables */
  done = FALSE;
  error_message[0] = '\0';
  got_rec = FALSE;
  line = 0;
  ndetails = 0;

  /* Open the domains file */
  if ((file_ptr = fopen(cath_domains_file,"r")) ==  NULL)
    {
      /* Generate error message about domains file being missing */
      strcpy(error_message,"CATH domains file not found.");
      return(ndetails);
    }

  /* Loop while reading in records from the file until get the
     one we're after */
  while (fgets(input_line,LINELEN,file_ptr) != NULL && done == FALSE)
    {
      /* Increment line-coint */
      line++;

      /* Get the line length */
      line_length = strlen(input_line);

      /* Get the chain-id */
      if (line_length > 14)
        {
          chain_id = input_line[4];
          wanted = FALSE;

          /* Check whether this record applies to the PDB code and chain
             that we want */
          if (!strncmp(input_line,pdb_code,4))
            {
              /* Determine whether chain is one of those we want */
              wanted = check_chain(chain_id,chains,nchains,&ichain);
            }

          /* If record is wanted, then process */
          if (wanted == TRUE)
            {
              /* Set flag to indicate we have found the record we're
                 after */
              got_rec = TRUE;

              /* Extract the domain data from the free-format line and
                 store in a fixed-format array for ease of use */
              /* Old version for old format CATH files
              ndetails = extract_domain_data(input_line,domain_details,
                                             &ndomains,&nfrags,
                                             error_message);
              */
              nfrags
                = extract_domain_info(input_line,wolf_code,cath_code,
                                      domain_details,ndetails);

              /* Increment count of details */
              ndetails = ndetails + nfrags;

              /* Extract the Wolf code from the front of the line */
              //strncpy(wolf_code,input_line,6);
              //wolf_code[6] = '\0';
            }

          /* If not wanted, see if we are done */
          else if (got_rec == TRUE)
            done = TRUE;
        }
    }

  /* Close the file */
  fclose(file_ptr);

  /* Return number of detail items found */
  return(ndetails);
}
/***********************************************************************

write_script_domains  -  Read in the domain definitions from the CATH
                        domains file and generate the RasMol script

***********************************************************************/

void write_script_domains(FILE *outfile_ptr,char pdb_code[5],
                          char cath_domains_file[FILENAME_LEN],
                          char first_object[OBJ_NAME_LEN],
                          char first_display_mode[DISPLAY_MODE_LEN],
                          char *chains,int nchains,
                          struct parameters params)
{
  char colour_name[COLNAM_LEN], colour_numbers[COLNO_LEN];
  char number_string[20];
  char residue_list[LINELEN + 1], object[OBJ_NAME_LEN];
  char domain_details[C_MAXDOMCHAIN][C_DOMLEN], error_message[LINELEN + 1];
  char domains_list[LINELEN + 1], fragments_list[LINELEN + 1];
  char last_detail[4], heading[20];

  int domain, domains_list_len, fragment, fragments_list_len;
  int idetail, ndetails;

  /* Initialise variables */
  domains_list[0] = '\0';
  fragments_list[0] = '\0';
  heading[0] = '\0';
  last_detail[0] = '\0';
  domains_list_len = 0;
  fragment = 0;
  fragments_list_len = 0;
  object[0] = '\0';
  residue_list[0] = '\0';
  strcpy(colour_name,"white");
  strcpy(colour_numbers,"white");

  /* Read in the domain definitions for the given chain */
  ndetails = read_in_domain_data(pdb_code,chains,nchains,cath_domains_file,
                                 domain_details,error_message);

  /* Loop through all the domain details, colouring each domain */
  for (idetail = 0; idetail < ndetails; idetail++)
    {
      /* If the details are for a new domain/fragment, then write out
         the last one */
      if (strncmp(domain_details[idetail],last_detail,3))
        {
          /* Check that there is something to write out */
          if (object[0] != '\0')
            {
              /* Write out heading */
              fprintf(outfile_ptr,"\n");
              fprintf(outfile_ptr,"%s\n",heading);

              /* Write out the object definition */
              fprintf(outfile_ptr,"define %s %s\n",object,residue_list);

              /* Select the range of residues defining this domain/fragment */
              fprintf(outfile_ptr,"select %s\n",object);
              fprintf(outfile_ptr,"cartoon on\n");
              fprintf(outfile_ptr,"colour atoms %s  # (%s)\n",
                      colour_numbers,colour_name);

              /* If this is the first object, then save its details */
              if (first_object[0] == '\0')
                {
                  strcpy(first_object,object);
                  strcpy(first_display_mode,"cartoon");
                }

              /* Add object to appropriate list */
              if (object[0] == 'd')
                domains_list_len = add_to_list(domains_list,
                                               domains_list_len,object);
              else
                fragments_list_len = add_to_list(fragments_list,
                                                 fragments_list_len,object);

              /* Reinitialise object details */
              object[0] = '\0';
              residue_list[0] = '\0';
            }

          /* Save the domain/fragment identifier */
          strncpy(last_detail,domain_details[idetail],3);
          last_detail[3] = '\0';
        }

      /* Extract the domain number */
      strncpy(number_string,domain_details[idetail]+1,2);
      number_string[2] = '\0';
      domain = atoi(number_string);

      /* If this is a domain definition, get its colour */
      if (domain_details[idetail][0] == 'D')
        {
          /* Get the corresponding domain colour */
          get_domain_colour(domain,colour_name,colour_numbers);

          /* Store heading */
          sprintf(heading,"# Domain %d",domain);

          /* Form the object name */
          sprintf(number_string,"%2d",domain);
          if (number_string[0] == ' ')
            number_string[0] = '0';
          strcpy(object,"domain");
          strcat(object,number_string);
        }

      /* Otherwise, for a fragment, set the colour to white */
      else
        {
          /* Get the fragment number */
          fragment++;

          /* Set colour to white */
          strcpy(colour_name,"darkgrey");
          strcpy(colour_numbers,"[128,128,144]");

          /* Store heading */
          sprintf(heading,"# Fragment %d",fragment);

          /* Form the object name */
          sprintf(number_string,"%2d",fragment);
          if (number_string[0] == ' ')
            number_string[0] = '0';
          strcpy(object,"frag");
          strcat(object,number_string);
        }

      /* Extract residue-range and add to residue-list */
      if (residue_list[0] != '\0')
        strcat(residue_list," or ");
      add_residue(residue_list,domain_details[idetail]+4);
      strcat(residue_list,"-");
      add_residue(residue_list,domain_details[idetail]+10);
    }

  /* Write out last domain/fragment definition, if there is one */
  if (object[0] != '\0')
    {
      /* Write out heading */
      fprintf(outfile_ptr,"\n");
      fprintf(outfile_ptr,"%s\n",heading);

      /* Write out the object definition */
      fprintf(outfile_ptr,"define %s %s\n",object,residue_list);

      /* Select the range of residues defining this domain/fragment */
      fprintf(outfile_ptr,"select %s\n",object);
      fprintf(outfile_ptr,"cartoon on\n");
      fprintf(outfile_ptr,"colour atoms %s  # (%s)\n",colour_numbers,
              colour_name);

      /* Add object to appropriate list */
      if (object[0] == 'd')
        domains_list_len = add_to_list(domains_list,
                                       domains_list_len,object);
      else
        fragments_list_len = add_to_list(fragments_list,
                                         fragments_list_len,object);

      /* Reinitialise object details */
      /* If this is the first object, then save its details */
      if (first_object[0] == '\0')
        {
          strcpy(first_object,object);
          strcpy(first_display_mode,"cartoon");
        }
    }

  /* Show the domain and fragment lists defined */

  /* If list of domains is long, then truncate it */
  if (domains_list_len > 18)
    {
      strncpy(domains_list+18," ...",4);
      domains_list[22] = '\0';
    }

  /* Show the list of domains */
  if (!strncmp(domains_list,"domain",6))
    {
      fprintf(outfile_ptr,"echo \"    %s     = separate domains",
              domains_list);
      fprintf(outfile_ptr,"   [shown as: cartoon]\"\n");
    }

  /* If list of fragments is long, then truncate it */
  if (fragments_list_len > 22)
    {
      strncpy(fragments_list+22," ...",4);
      fragments_list[26] = '\0';
    }

  /* Show the list of fragments */
  if (!strncmp(fragments_list,"frag",4))
    {
      fprintf(outfile_ptr,"echo \"    %s = inter-domain fragments",
              fragments_list);
      fprintf(outfile_ptr,"   [shown in grey]\"\n");
    }

  /* If no domain details read in, then show error meassage */
  if (ndetails == 0 || error_message[0] != '\0')
    {
      /* If have an error message from reading the CATH domains file,
         then show that */
      if (error_message[0] != '\0')
        {
          fprintf(outfile_ptr,"\n");
          fprintf(outfile_ptr,"echo \" \"\n");
          fprintf(outfile_ptr,"echo \"*** ERROR reading CATH domains file:-");
          fprintf(outfile_ptr,"\"\n");
          fprintf(outfile_ptr,"echo \"          %s\n",error_message);
          fprintf(outfile_ptr,"echo \" \"\n");
          fprintf(outfile_ptr,"\n");
        }

      /* Otherwise, say that no domain info found */
      else
        {
          fprintf(outfile_ptr,"\n");
          fprintf(outfile_ptr,"echo \" \"\n");
          fprintf(outfile_ptr,"echo \"*** ERROR. No domains data found ");
          fprintf(outfile_ptr,"in CATH domains file for this chain\"\n");
          fprintf(outfile_ptr,"echo \" \"\n");
          fprintf(outfile_ptr,"\n");
        }

      /* Select whole structure for display */
      fprintf(outfile_ptr,"\n");
      fprintf(outfile_ptr,"# All\n");
      fprintf(outfile_ptr,"select all\n");
      fprintf(outfile_ptr,"wireframe on\n");
      fprintf(outfile_ptr,"\n");

      /* Note objects selected */
      strcpy(first_object,"all");
      strcpy(first_display_mode,"wireframe");
    }
}
/***********************************************************************

generate_script  -  Read in the RasMol definitions from the rasmol.dfn
                    file and generate the RasMol script

***********************************************************************/

void generate_script(FILE *outfile_ptr,char *pdbsum_dir,char *dfn_name,
                     char *work_dir,char *first_object,
                     char *first_display_mode,char *chains,int nchains,
                     char *domsst,struct parameters params)
{
  char ch, dominant_sst;
  char file_name[FILENAME_LEN], input_line[LINELEN + 1];
  char tmp_line[LINELEN + 1];
  char chain, ligand_num[4], metal_name[5], site_name[4];
  char metal_colour[COLNAM_LEN], colour_numbers[COLNO_LEN];
  char number_string[10];
  char object[OBJ_NAME_LEN], display_mode[DISPLAY_MODE_LEN];
  char nchains_list[LINELEN + 1], pchains_list[LINELEN + 1];
  char mchains_list[LINELEN + 1], surfaces_list[LINELEN + 1];
  char ligands_list[MAXLIGANDS][MAX_RANGE_LEN], sites_list[MAXSITE][4];
  char dot_surfaces_list[LINELEN + 1];

  char *string_ptr;

  int atom_from, atom_to, i, ichain, line;
  int mlist_len, nlist_len, nnuc_chain, plist_len, slist_len;
  int nligands, nprot_chain, nsites, nsurfaces;
  int first, have_backbone, have_cartoon, have_dots, have_mask, wanted;

  FILE *file_ptr;

  /* Initialise */
  mchains_list[0] = '\0';
  nchains_list[0] = '\0';
  pchains_list[0] = '\0';
  dot_surfaces_list[0] = '\0';
  surfaces_list[0] = '\0';
  mlist_len = 0;
  nlist_len = 0;
  plist_len = 0;
  slist_len = 0;
  nligands = 0;
  nsites = 0;
  nsurfaces = 0;
  display_mode[0] = '\0';
  first = TRUE;
  have_backbone = FALSE;
  have_cartoon = FALSE;
  have_dots = FALSE;
  have_mask = FALSE;
  line = 0;
  nprot_chain = 0;
  nnuc_chain = 0;
  object[0] = '\0';

  /* Form filename of the rasmol.dfn file */
  if (work_dir[0] == '\0')
    strcpy(file_name,pdbsum_dir);
  else
    {
      strcpy(file_name,work_dir);
      strcat(file_name,DIR_SEP);
    }
  if (params.pdb_type == PDBSUM1)
    strcpy(file_name,"../newsec/");
  strcat(file_name,dfn_name);

  /* Open the rasmol.dfn input file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    {
      /* Unable to open RasMol definitions file, so display an error
         message, unless colouring by domain */
      if (params.colour_by_domain == FALSE)
        {
          /* Show error message */
          fprintf(outfile_ptr,"echo \"*** Unable to open RasMol definitions");
          fprintf(outfile_ptr," file:-\"\n");
          fprintf(outfile_ptr,"echo \"       [%s]\"\n",file_name);

          /* Select whole structure for display */
          fprintf(outfile_ptr,"\n");
          fprintf(outfile_ptr,"# All\n");
          fprintf(outfile_ptr,"select all\n");
          fprintf(outfile_ptr,"wireframe on\n");
          fprintf(outfile_ptr,"\n");

          /* Note objects selected */
          strcpy(first_object,"all");
          strcpy(first_display_mode,"wireframe");
        }
      return;
    }

  /* Loop through the file processing each record */
  while (fgets(input_line,LINELEN,file_ptr)!= NULL)
    {
      /* Increment line-count */
      line++;

      /* For single-chain display, only process record if it belongs
         to one of the selected chains */
      if (params.all_chains == FALSE)
        {
          /* If protein or DNA, check if this is the right chain */
          if (input_line[0] == 'P' || input_line[0] == 'C' ||
              input_line[0] == 'N')
            {
              /* Get chain identifier from input line */
              ch = input_line[1];

              /* Check whether this chain is wanted */
              wanted = check_chain(ch,chains,nchains,&ichain);

              /* If don't have a match on chain, knock this entry out */
              if (wanted == FALSE)
                input_line[0] = ' ';
            }

          /* If ligand or metal, then don't include */
          else if (input_line[0] == 'l' || input_line[0] == 'L' ||
                   input_line[0] == 'M')
            input_line[0] = ' ';

          /* Otherwise, check for active-site residues */
          else if (input_line[0] == 'A')
            {
              /* Check that have a chain defined */
              if (params.all_chains == FALSE)
                {
                  /* Check for colons */
                  string_ptr = strstr(input_line,":");
                  if (string_ptr != NULL)
                    {
                      /* Get next character */
                      ch = string_ptr[1];

                      /* Check whether this chain is wanted */
                      wanted = check_chain(ch,chains,nchains,&ichain);

                      /* If doesn't match chain, then drop */
                      if (wanted == FALSE)
                        input_line[0] = ' ';
                    }
                }
            }
        }

      /* Check line-type according to 1st character: A, C, L, l, M, N, P
         or S */

      /* If active site, then write the script records to define
         and display the relevant residues */
      if (input_line[0] == 'A' && params.define_sites == TRUE)
        {
          /* Extract the 3-character site name */
          strncpy(site_name,input_line+1,3);
          site_name[3] = '\0';

          /* Write out the active-site definition and display commands
             to the script */
          write_script_site(outfile_ptr,site_name,input_line+5,object,
                            display_mode,params);

          /* Store this active site's name */
          if (nsites < MAXSITE - 1)
            {
              strcpy(sites_list[nsites],site_name);
              nsites++;
            }
        }

      /* If protein CA-only chain, then write corresponding script
         records */
      else if (input_line[0] == 'C' && params.define_protein == TRUE)
        {
          dominant_sst = ' ';
          write_script_protein(outfile_ptr,input_line[1],TRUE,object,
                               display_mode,pchains_list,&plist_len,
                               dominant_sst,params);
          have_backbone = TRUE;
          nprot_chain++;
        }

      /* Check for ligand identifier */
      else if (input_line[0] == 'L' && params.define_ligands == TRUE)
        {
          /* Extract the ligand number */
          strcpy(tmp_line,input_line+1);
          string_ptr = strstr(tmp_line," ");
          if (string_ptr != NULL)
            string_ptr[0] = '\0';
          else
            string_ptr[2] = '\0';
          strcpy(ligand_num,tmp_line);

          /* Extract the ligand name */
          string_ptr = strstr(input_line," ");
          if (string_ptr != NULL)
            strcpy(tmp_line,string_ptr+1);
          else
            strcpy(tmp_line,"error");

          /* Write out the ligand definition and display commands to the
             script */
          write_script_ligand(outfile_ptr,ligand_num,tmp_line,object,
                              display_mode,params);
        }

      /* Check for record holding the ligand residue range */
      else if (input_line[0] == 'l' && params.define_ligands == TRUE)
        {
          /* Check that have space to store this range */
          if (nligands < MAXLIGANDS - 1)
            {
              /* Terminate the line after last non-blank character */
              terminate_line(input_line,MAX_RANGE_LEN + 4);

              /* Extract the ligand residue range */
              string_ptr = strstr(input_line," ");
              if (string_ptr != NULL)
                strcpy(tmp_line,string_ptr+1);
              else
                strcpy(tmp_line,"error");

              /* Store the ligand residue range */
              strcpy(ligands_list[nligands],tmp_line);
              nligands++;
            }
        }

      /* Check for metal */
      else if (input_line[0] == 'M' && params.define_metals == TRUE)
        {
          /* Extract the metal name */
          strncpy(metal_name,input_line+1,2);
          metal_name[2] = '\0';
          strcat(metal_name,"  ");

          /* Get the metal colour */
          string_ptr = strstr(input_line + 2," ");
          if (string_ptr != NULL)
            strcpy(tmp_line,string_ptr+1);
          else
            strcpy(tmp_line,"white");

          /* Terminate the colour after last non-blank character */
          terminate_line(tmp_line,COLNAM_LEN - 1);
          strcpy(metal_colour,tmp_line);

          /* Get the corresponding RasMol colour triplet */
          get_rasmol_colour(metal_colour,colour_numbers);

          /* Write out the metal definition and display commands to the
             script */
          write_script_metal(outfile_ptr,metal_name,colour_numbers,object,
                             display_mode,mchains_list,&mlist_len,params);
        }

      /* Check for nucleic acid chain */
      else if (input_line[0] == 'N' && params.define_nucleic == TRUE)
        {
          write_script_nucleic(outfile_ptr,input_line[1],object,
                               display_mode,nchains_list,&nlist_len,params);
          nnuc_chain++;
        }

      /* If protein chain, then write the script records to define
         and display the chain */
      else if (input_line[0] == 'P' && params.show_protein == TRUE)
        {
          /* Determine whether this chain has a dominant secondary
             structure */
          dominant_sst = ' ';
          for (i = 0; i < nchains && dominant_sst == ' '; i++)
            {
              /* If this is the right chain, save dominant secondary
                 structure */
              if (input_line[1] == chains[i])
                dominant_sst = domsst[i];
            }

          /* Write out to script file */
          write_script_protein(outfile_ptr,input_line[1],FALSE,object,
                               display_mode,pchains_list,&plist_len,
                               dominant_sst,params);
          have_cartoon = TRUE;
          nprot_chain++;
        }

      /* Check for dot surfaces */
      else if (input_line[0] == 'S' && params.define_dots == TRUE)
        {
          /* Extract the chain-id of this surface */
          chain = input_line[1];

          /* Extract the range of atom numbers defining this surface */
          strncpy(number_string,input_line+2,8);
          number_string[8] = '\0';
          atom_from = atoi(number_string);
          strncpy(number_string,input_line+10,8);
          number_string[8] = '\0';
          atom_to = atoi(number_string);

          /* Write out the script commands for this surface */
          write_script_dots(outfile_ptr,chain,atom_from,atom_to,
                            dot_surfaces_list);
          have_dots = TRUE;
        }

      /* Check for surfaces (not currently used) */
      else if (input_line[0] == 'X' && params.define_surfaces == TRUE)
        {
          /* Extract the chain-id of this surface */
          chain = input_line[1];

          /* Increment count of surfaces */
          nsurfaces++;

          /* Write out the script commands for this surface */
          write_script_surface(outfile_ptr,chain,nsurfaces,surfaces_list,
                               &slist_len,params);
        }

      /* If this is the first object, store its name and display-mode */
      if (first == TRUE && object[0] != '\0' && display_mode[0] != '\0')
        {
          /* Save object and display */
          strcpy(first_object,object);
          strcpy(first_display_mode,display_mode);

          /* Set flag to indicate that first object done */
          first = FALSE;
        }
    }

  /* If have encountered any surfaces, write out the mask surface */
  if (have_dots == TRUE)
    {
      nsurfaces++;
      write_script_mask(outfile_ptr,nsurfaces,params);
      have_mask = TRUE;
    }

  /* Write out all the user-defined sets we have identified */
  write_script_sets(outfile_ptr,ligands_list,nligands,mchains_list,
                    mlist_len,nchains_list,nlist_len,pchains_list,
                    plist_len,sites_list,nsites,surfaces_list,slist_len,
                    dot_surfaces_list,have_backbone,have_cartoon,
                    have_dots,have_mask,params);
}
/***********************************************************************

process_conect_record  -  Check whether the given CONECT record applies
                          to one of the ATOM or HETATM records already
                          written out and hence is required

***********************************************************************/

int process_conect_record(char input_line[LINELEN + 1],
                          struct atom *first_atom_ptr)
{
  char atmnum[10];

  int atom_number;
  int wanted, write_out;

  struct atom *atom_ptr;

  /* Initialise variables */
  write_out = FALSE;

  /* Get the atom-number of the first number in the CONECT record */
  strncpy(atmnum,input_line+6,5);
  atmnum[5] = '\0';
  atom_number = atoi(atmnum);

  /* Check whether this atom number is one of those already written out */
  wanted = FALSE;

  /* Initialise pointer to the first record */
  atom_ptr = first_atom_ptr;

  /* Loop through all the atom records */
  while (atom_ptr != NULL && wanted == FALSE)
    {
      /* If have a match, then set flag */
      if (atom_number == atom_ptr->atom_number)
          wanted = TRUE;

      /* Get pointer to the next atom record */
      atom_ptr = atom_ptr->next_atom_ptr;
    }

  /* If record is wanted, set flag to write it out */
  if (wanted == TRUE)
    write_out = TRUE;

  /* Return whether this CONECT record is wanted or not */
  return(write_out);
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

get_atom_details  -  Extract the atom details from the given line

***********************************************************************/

void get_atom_details(char *input_line,int *atom_number,char *atom_name,
                      char *chain,char *res_name,char *res_num,
                      double *xval,double *yval,double *zval,
                      double *bval,double *occup)
{
  char number_string[20];

  int len;

  double bvalue, occupancy, x, y, z;

  /* Get the length of the input line */
  len = strlen(input_line);

  /* Get the residue details */
  strncpy(res_name,input_line+17,3);
  res_name[3] = '\0';
  strncpy(res_num,input_line+22,5);
  res_num[5] = '\0';
  *chain = input_line[21];

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

  /* Store the retrieved values */
  *xval = x;
  *yval = y;
  *zval = z;
  *bval = bvalue;
  *occup = occupancy;
}
/***********************************************************************

create_residue_record  -  Create a new residue record

***********************************************************************/

struct residue *create_residue_record(struct residue **fst_residue_ptr,
                                      struct residue **lst_residue_ptr,
                                      char res_name[4],char res_num[6],
                                      char chain,int type,
                                      struct rasdef *rasdef_ptr)
{
  struct residue *residue_ptr;
  struct residue *last_residue_ptr, *first_residue_ptr;


  /* Initialise residue pointers */
  residue_ptr = NULL;
  last_residue_ptr = *lst_residue_ptr;
  first_residue_ptr = *fst_residue_ptr;

  /* Create residue entry */
  residue_ptr = (struct residue*) malloc(sizeof(struct residue));
  if (residue_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for");
      printf(" struct residue\n");
      printf("***        Program rascript terminated with error.\n");
      exit (1);
    }

  /* Store the details for this entry */
  strcpy(residue_ptr->res_name,res_name);
  strcpy(residue_ptr->res_num,res_num);
  residue_ptr->chain = chain;
  residue_ptr->type = type;

  /* Initialise the remaining residue details */
  residue_ptr->sec_struc = COIL;
  residue_ptr->natoms = 0;
  residue_ptr->rasdef_ptr = rasdef_ptr;
  residue_ptr->first_atom_ptr = NULL;
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

check_for_duplicate  -  Check for an earlier copy of the current atom

***********************************************************************/

int check_for_duplicate(struct residue *residue_ptr,char *atom_name,
                        char *res_name,char *res_num,char chain)
{
  int iatom, natoms;
  int have_atom;

  struct atom *atom_ptr;

  /* Initialise variables */
  have_atom = FALSE;

  /* If have a current residue, and its name matches, check its atoms
     for an existing copy of the current atom */
  if (residue_ptr != NULL &&
      !strcmp(res_name,residue_ptr->res_name) &&
      !strcmp(res_num,residue_ptr->res_num) &&
      chain == residue_ptr->chain)
    {
      /* Get number of atoms */
      natoms = residue_ptr->natoms;

      /* Initialise pointer to first atom of this residue */
      atom_ptr = residue_ptr->first_atom_ptr;
      iatom = 0;
  
      /* Loop through all the residue's atoms to look for match */
      while (iatom < natoms && atom_ptr != NULL && have_atom == FALSE)
        {
          /* Check for duplication */
          if (!strcmp(atom_name,atom_ptr->atom_name))
            have_atom = TRUE;

          /* Get pointer to next atom in linked-list */
          atom_ptr = atom_ptr->next_atom_ptr;
          iatom++;
        }
    }

  /* Return flag */
  return(have_atom);
}
/***********************************************************************

create_atom_record  -  Create a new atom record

***********************************************************************/

struct atom *create_atom_record(struct atom **fst_atom_ptr,
                                struct atom **lst_atom_ptr,
                                struct residue *residue_ptr,
                                int atom_number,char *atom_name,
                                double x,double y,double z,
                                double bvalue,double occupancy,
                                char *element)
{
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
      printf("***        Program contour terminated with error.\n");
      exit (1);
    }

  /* Store the details for this entry */
  atom_ptr->atom_number = atom_number;
  strcpy(atom_ptr->atom_name,atom_name);
  strcpy(atom_ptr->element,element);
  atom_ptr->x = x;
  atom_ptr->y = y;
  atom_ptr->z = z;
  atom_ptr->bvalue = bvalue;
  atom_ptr->occupancy = occupancy;
  atom_ptr->residue_ptr = residue_ptr;

  /* Increment count of this residue's atoms */
  residue_ptr->natoms++;

  /* If residue's first atom, store pointer to it */
  if (residue_ptr->first_atom_ptr == NULL)
    residue_ptr->first_atom_ptr = atom_ptr;

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

get_max_min_coords  -  Get the maximum and minimum x-, y-, and z-coords
                       of all atoms

***********************************************************************/

void get_max_min_coords(struct limits *limit,double x,double y,double z)
{
  /* If this is the first atom encountered, store its
     coords as the present max and min values */
  if (limit->first == TRUE)
    {
      limit->valmin[0] = x;
      limit->valmax[0] = x;
      limit->valmin[1] = y;
      limit->valmax[1] = y;
      limit->valmin[2] = z;
      limit->valmax[2] = z;
      limit->first = FALSE;
      limit->npoints = 1;
    }

  /* Otherwise, update the limits if coordinates are outside them */
  else
    {
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

check_interaction  -  Check whether the current residue interacts with
                      any of the protein chains under consideration

***********************************************************************/

int check_interaction(char *res_name,char *res_num,char chain,
                      struct grow *first_grow_ptr,char *chains,
                      int nchains)
{
  char other_chain;

  int wanted;
  int ichain, match;

  struct grow *grow_ptr;

  /* Assume not wanted */
  wanted = FALSE;

  /* Get pointer to the first grow record */
  grow_ptr = first_grow_ptr;

  /* Loop through all the grow records */
  while (grow_ptr != NULL && wanted == FALSE)
    {
      /* Initialise flag */
      match = 0;

      /* Check whether current residue is present in this grow record */
      if (!strcmp(res_name,grow_ptr->res_name1) &&
          !strcmp(res_num,grow_ptr->res_num1) &&
          chain == grow_ptr->chain1)
        match = 1;
      else if (!strcmp(res_name,grow_ptr->res_name2) &&
               !strcmp(res_num,grow_ptr->res_num2) &&
               chain == grow_ptr->chain2)
        match = 2;

      /* If it is present, see if the residue it interacts with belongs
         to one of the wanted chains */
      if (match > 0)
        {
          /* Get the chain of the other residue */
          if (match == 1)
            other_chain = grow_ptr->chain2;
          else
            other_chain = grow_ptr->chain1;
        }

      /* Loop through all the wanted chains to see if
         this is among them */
      for (ichain = 0; ichain < nchains && wanted == FALSE; ichain++)
        {
          /* Check this chain */
          if (other_chain == chains[ichain])
            wanted = TRUE;
        }

      /* Get pointer to next grow record */
      grow_ptr = grow_ptr->next_grow_ptr;
    } 

  /* Return whether residue is wanted */
  return(wanted);
}
/***********************************************************************

create_conect_record  -  Create a new conect record

***********************************************************************/

struct conect *create_conect_record(struct conect **fst_conect_ptr,
                                    struct conect **lst_conect_ptr,
                                    int atom_number[MAXCONECT],
                                    int nconect)
{
  int iconect;

  struct conect *conect_ptr;
  struct conect *first_conect_ptr, *last_conect_ptr;

  /* Initialise conect pointers */
  conect_ptr = NULL;
  last_conect_ptr = *lst_conect_ptr;
  first_conect_ptr = *fst_conect_ptr;

  /* Create conect entry */
  conect_ptr = (struct conect*)malloc(sizeof(struct conect));
  if (conect_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for");
      printf(" struct conect\n");
      printf("***        Program rasurf terminated with error.\n");
      exit (1);
    }

  /* Store the details for this entry */
  conect_ptr->nconect = nconect;
  for (iconect = 0; iconect < nconect; iconect++)
    conect_ptr->atom_number[iconect] = atom_number[iconect];

  /* Initialise pointer to next conect record */
  conect_ptr->next_conect_ptr = NULL;

  /* If this is the first entry stored, then update
     pointer to head of linked list */
  if (first_conect_ptr == NULL)
    first_conect_ptr = conect_ptr;

  /* Otherwise, make previous entry point to the current
     one */
  else
    last_conect_ptr->next_conect_ptr = conect_ptr;
                      
  /* Save the current conect's pointer */
  last_conect_ptr = conect_ptr;

  /* Return current pointers */
  *fst_conect_ptr = first_conect_ptr;
  *lst_conect_ptr = last_conect_ptr;
  return(conect_ptr);
}
/***********************************************************************

get_sec_struc  -  Use the data read in from the PROMOTIF file to determine
                  given residue's secondary structure

***********************************************************************/

int get_sec_struc(struct promotif *first_promotif_ptr,char *res_num,
                  char chain,int *last_sec_struc)
{
  int sec_struc;

  struct promotif *promotif_ptr;

  /* Initialise variables */
  sec_struc = *last_sec_struc;

  /* Get pointer to the first of the PROMOTIF definitions */
  promotif_ptr = first_promotif_ptr;

  /* Loop through the PROMOTIF definitions */
  while (promotif_ptr != NULL)
    {
      /* Only consider secondary structure assignments of interest */
      if (promotif_ptr->sec_struc == 'H' ||
          promotif_ptr->sec_struc == 'E')
        {
          /* If our residue matches the start of a PROMOTIF definition,
             then set its secondary structure accordingly */
          if (!strcmp(res_num,promotif_ptr->res_num1) &&
              chain == promotif_ptr->chain1)
            {
              /* Get current secondary structure type */
              if (promotif_ptr->sec_struc == 'H')
                sec_struc = HELIX;
              else if (promotif_ptr->sec_struc == 'E')
                sec_struc = STRAND;
              else
                sec_struc = COIL;

              /* Set last secondary structure to match this */
              *last_sec_struc = sec_struc;
            }

          /* Otherwise, check whether it matches the last residue in the
             PROMOTIF definition */
          else if (!strcmp(res_num,promotif_ptr->res_num2) &&
                   chain == promotif_ptr->chain2)
            {
              /* Set next residue to revert to COIL */
              *last_sec_struc = COIL;
            }
        }

      /* Get pointer to the next PROMOTIF data item */
      promotif_ptr = promotif_ptr->next_promotif_ptr;
    }

  /* Return the residue's secondary structure */
  return(sec_struc);
}
/***********************************************************************

get_atom_records  -  Read through the PDB file, writing out all the
                     atom and CONECT records for the residue-list read
                     in from ligplot.res

***********************************************************************/

void get_atom_records(char *pdb_name,struct rasdef *first_rasdef_ptr,
                      struct promotif *first_promotif_ptr,
                      struct grow *first_grow_ptr,
                      struct atom **fst_atom_ptr,
                      struct residue **fst_residue_ptr,char *chains,
                      int nchains,int have_dfn,int have_grow,
                      struct conect **fst_conect_ptr,
                      struct limits *atom_limits,
                      struct parameters params)
{
  char ch, input_line[LINELEN + 1], number_string[20];
  char atom_name[5], element[3];
  char res_name[4], res_num[6], last_res_name[4], last_res_num[6];
  char chain, last_chain;

  int atom_counter, atom_number, ichain, nresid, type;
  int iconect, ipos, len, nconect;
  int last_sec_struc;
  int keep_reading, is_metal;
  int done, eof, have_atom, wanted_residue, write_out;
  int atom_record, have_sst;
  int store_number[MAXCONECT];
  int gzipped;

  double bvalue, occupancy, x, y, z;

  struct atom *atom_ptr, *first_atom_ptr, *last_atom_ptr;
  struct conect *conect_ptr, *first_conect_ptr, *last_conect_ptr;
  struct rasdef *rasdef_ptr;
  struct residue *first_residue_ptr, *last_residue_ptr;
  struct residue *residue_ptr;
    
#ifdef TEST
  FILE *file_ptr;
#else
  FILE *file_ptr;
  gzFile gzfile_ptr;
#endif

  /* Initialise variables */
  eof = FALSE;
  last_sec_struc = COIL;

  /* Determine whether the PDB file is a gzipped one */
  gzipped = FALSE;
  len = strlen(pdb_name);
  if (!strncmp(pdb_name + len - 3,".gz",3))
    gzipped = TRUE;
#ifdef TEST
  gzipped = FALSE;
#endif

  /* Open PDB file */
#ifdef TEST
  file_ptr = fopen(pdb_name,"r");
  if (file_ptr ==  NULL)
    eof = TRUE;
#else
  if (gzipped == TRUE)
    {
      gzfile_ptr = gzopen(pdb_name,"r");
      if (gzfile_ptr ==  NULL)
        eof = TRUE;
    }
  else
    {
      file_ptr = fopen(pdb_name,"r");
      if (file_ptr ==  NULL)
        eof = TRUE;
    }
#endif

  /* If file not opened, then abort */
  if (eof ==  TRUE)
    {
      printf("\n*** Unable to open your PDB file [%s]\n",pdb_name);
      exit(1);
    }

  /* Initialise variables */
  atom_counter = 0;
  atom_limits->first = TRUE;
  keep_reading = TRUE;
  first_atom_ptr = NULL;
  first_conect_ptr = NULL;
  last_conect_ptr = NULL;
  have_sst = FALSE;
  last_atom_ptr = NULL;
  first_residue_ptr = NULL;
  last_residue_ptr = NULL;
  last_res_name[0] = '\0';
  last_res_num[0] = '\0';
  last_chain = '\0';
  nresid = 0;
  type = UNKNOWN;
  wanted_residue = FALSE;

  /* Loop through the PDB file, picking up all the atoms that are
     required */
  while (eof == FALSE && keep_reading == TRUE)
    {
      /* Read in the next record */
      if (gzipped == TRUE)
        {
#ifdef TEST
#else
          if (gzgets(gzfile_ptr,input_line,LINELEN) == NULL)
            eof = TRUE;
#endif
        }
      else
        {
          if (fgets(input_line,LINELEN,file_ptr) == NULL)
            eof = TRUE;
        }

      /* If record read in OK, then process it */
      if (eof == FALSE)
        {
          /* Initialise flags */
          atom_record = FALSE;
          write_out = FALSE;

          /* Check for the end of an NMR structure model */
          if (!strncmp(input_line,"ENDMDL",6))
            keep_reading = FALSE;

          /* If one of the records giving the secondary structure
             assignments, then write out */
          else if (!strncmp("SHEET ",input_line,6) ||
                   !strncmp("HELIX ",input_line,6))
            {
              /* Skip if writing Molscript output */
              if (params.output_option == OUT_MOLSCRIPT ||
                  params.output_option == OUT_RASTER3D)
                write_out = FALSE;
          
              else
                {
                  have_sst = TRUE;
                  write_out = TRUE;
                }
            }

          /* If this is a CONECT record, then see if it is wanted */
          else if (!strncmp(input_line,"CONECT",6))
            {
              /* Initialise variables */
              done = FALSE;
              ipos = 6;
              nconect = 0;

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
                      atom_number = atoi(number_string);

                      /* Store, if can */
                      if (nconect < MAXCONECT)
                        {
                          /* Store number */
                          store_number[nconect] = atom_number;
                          nconect++;
                        }
                    }

                  /* Increment position in line */
                  ipos = ipos + 5;
                }

              /* If have at least two atom numbers stored, create a
                 conect record */
              if (nconect > 1)
                conect_ptr
                  = create_conect_record(&first_conect_ptr,
                                         &last_conect_ptr,store_number,
                                         nconect);
            }

          /* If this is an ATOM or HETATM record, process it */
          else if (keep_reading == TRUE &&
                   (!strncmp(input_line,"ATOM",4) ||
                    !strncmp(input_line,"HETATM",6)))
            {
              /* Set flag that this is an atom record */
              atom_record = TRUE;

              /* Get the atom details from the current line */
              get_atom_details(input_line,&atom_number,atom_name,
                               &chain,res_name,res_num,&x,&y,&z,
                               &bvalue,&occupancy);

              /* Check whether atom belongs to a new residue */
              if (chain != last_chain || strncmp(res_num,last_res_num,5) ||
                  strncmp(res_name,last_res_name,3))
                {
                  /* Initialise pointer and residue type */
                  rasdef_ptr = NULL;
                  type = UNKNOWN;
                  wanted_residue = TRUE;

                  /* If residue is a water, then don't want it */
                  if (!strncmp(res_name,"HOH",3))
                    wanted_residue = FALSE;

                  /* Determine which of the objects, if any, in the
                     rasmol.dfn file this residue belongs to */
                  else
                    rasdef_ptr
                      = assign_to_object(first_rasdef_ptr,atom_name,
                                         res_name,res_num,chain);

                  /* If have an assignment, see what it is */
                  if (rasdef_ptr != NULL)
                    {
                      /* Store residue's object type */
                      type = rasdef_ptr->type;

                      /* If a protein chain, then check whether it
                         belongs to a chain that is wanted */
                      if (type == PROTEIN || type == NUCLEIC_ACID ||
                          type == CA_ONLY)
                        {
                          /* If all chains required, then want this residue
                             whatever it is */
                          if (params.all_chains == TRUE)
                            wanted_residue = TRUE;

                          /* Otherwise, check chain */
                          else
                            wanted_residue
                              = check_chain(chain,chains,nchains,&ichain);
                        }

                      /* If not a protein chain, then need to check the
                         grow records to see whether residue interacts
                         with one of the desired chains */
                      else if (type == LIGAND || type == METAL)
                        {
                          /* If have grow data, check whether this metal or
                             ligand interacts with one of the chosen chains */
                          if (have_grow == TRUE)
                            wanted_residue
                              = check_interaction(res_name,res_num,chain,
                                                  first_grow_ptr,chains,
                                                  nchains);

                          /* Otherwises, take as wanted for the time being */
                          else
                            wanted_residue = TRUE;
                        }
                    }
                  else if (wanted_residue == TRUE)
                    {
                      /* If all chains required, then want this residue
                         whatever it is */
                      if (params.all_chains == TRUE)
                        wanted_residue = TRUE;

                      /* Otherwise, check chain */
                      else
                        wanted_residue
                          = check_chain(chain,chains,nchains,&ichain);
                    }

                  /* If residue is wanted, then create a residue record */
                  if (wanted_residue == TRUE)
                    {
                      residue_ptr
                        = create_residue_record(&first_residue_ptr,
                                                &last_residue_ptr,res_name,
                                                res_num,chain,type,
                                                rasdef_ptr);

                      /* Increment count of residues read in */
                      nresid++;

                      /* If residue is a protein residue, use the data
                         read in from the PROMOTIF file to determine its
                         secondary structure */
                      if (type == PROTEIN)
                        residue_ptr->sec_struc
                          = get_sec_struc(first_promotif_ptr,res_num,chain,
                                          &last_sec_struc);
                      else if (type == NUCLEIC_ACID)
                        residue_ptr->sec_struc = DHELIX;
                    }

                  /* Save the current residue number and chain */
                  strncpy(last_res_name,res_name,3);
                  last_res_name[3] = '\0';
                  strncpy(last_res_num,res_num,5);
                  last_res_num[5] = '\0';
                  last_chain = chain;
                }

              /* Set flag to write this record out */
              write_out = TRUE;

              /* If residue not wanted, then atom isn't, either */
              if (wanted_residue == FALSE)
                write_out = FALSE;

              /* Otherwise check atom */
              else
                {
                  /* Check whether this is an alternate occupancy atom that
                     we already have */
                  if (input_line[16] != ' ' || occupancy < 1.0)
                    {
                      /* Check for earlier copy of this atom */
                      have_atom
                        = check_for_duplicate(residue_ptr,atom_name,
                                              res_name,res_num,chain);
          
                      /* If have this atom already, then don't want this
                         second copy */
                      if (have_atom == TRUE)
                        write_out = FALSE;
                    }
                }

              /* If writing for RasMol script, exclude any unwanted atoms */
              if (write_out == TRUE)
                {
                  /* Determine the atom's element from it name */
                  get_element(atom_name,element);

                  /* Check whether this atom is a metal ion */
                  is_metal = FALSE;
                  if (atom_name[0] != ' ')
                    is_metal = r_check_for_metal(atom_name,element);

                  /* If this is a metal ion, then change its chain ID */
                  if (is_metal == TRUE && strcmp(res_name,"MSE"))
                    {
                      input_line[21] = 'M';
                      strncpy(input_line,"ATOM  ",6);
                    }

                  /* Otherwise, if hydrogen atom, then don't want to
                     write it out */
                  else if (!strcmp(element," H"))
                    write_out = FALSE;

                  /* Check for unusual atom-naming of some residues
                     which cause Rasmol to go awry */
                  if (input_line[12] == 'A' || input_line[12] == 'N' ||
                      input_line[12] == 'C' || input_line[12] == 'G')
                    {
                      /* If this is C, N, O, or P, then blank out the initial
                         A or N */
                      ch = input_line[13];
                      if (ch == 'C' || ch == 'N' || ch == 'O' || ch == 'P')
                        input_line[12] = ' ';
                    }

                  /* If modifying MSE residues, then change to Met */
                  if (params.change_MSE == TRUE && !strcmp(res_name,"MSE"))
                    {
                      /* Change residue name */
                      strncpy(input_line+17,"MET",3);

                      /* If this is the SE record, drop altogether */
                      if (!strcmp(atom_name,"SE  "))
                        write_out = FALSE;
                    }
                }

              /* If atom wanted, then store its details */
              if (write_out == TRUE)
                {
                  /* Create an atom record */
                  atom_ptr
                    = create_atom_record(&first_atom_ptr,&last_atom_ptr,
                                         residue_ptr,atom_number,atom_name,
                                         x,y,z,bvalue,occupancy,element);

                  /* If have an associated rasdef definition, add atom
                     to it */
                  if (rasdef_ptr != NULL)
                    {
                      /* If first atom, save */
                      if (rasdef_ptr->first_atom_ptr == NULL)
                        rasdef_ptr->first_atom_ptr = atom_ptr;
                      rasdef_ptr->last_atom_ptr = atom_ptr;
                    }

                  /* Increment count of atoms read in */
                  atom_counter++;

                  /* If ligand or metal, update its atom limits */
                  if (type == LIGAND || type == METAL)
                    get_max_min_coords(&(rasdef_ptr->atom_limits),x,y,z);
          
                  /* Otherwise, update the maximum and minimum atom
                     coordinates for protein */
                  else
                    get_max_min_coords(atom_limits,x,y,z);
                }
            }
        }
    }

  /* Close PDB file */
  if (gzipped == TRUE)
    {
#ifdef TEST
#else
      gzclose(gzfile_ptr);
#endif
    }
  else
    fclose(file_ptr);

  /* Return pointers to the first of the atoms and residues */
  *fst_atom_ptr = first_atom_ptr;
  *fst_residue_ptr = first_residue_ptr;
  *fst_conect_ptr = first_conect_ptr;
}
/***********************************************************************

get_interacting_ligands  -  Delete all ligands and metals that aren't
                            interacting with the selected chain(s)

***********************************************************************/

void get_interacting_ligands(struct rasdef *first_rasdef_ptr,
                             struct limits atom_limits)
{
  int include;

  struct rasdef *rasdef_ptr;

  /* Get pointer to the first RasMol definition record */
  rasdef_ptr = first_rasdef_ptr;

  /* Loop through all the definition records */
  while (rasdef_ptr != NULL)
    {
      /* Check for a match to ligand or metal */
      if (rasdef_ptr->type == LIGAND || rasdef_ptr->type == METAL)
        {
          /* Check whether the atom limits overlap */
          include = TRUE;
          if (rasdef_ptr->atom_limits.valmin[0] > atom_limits.valmax[0] ||
              rasdef_ptr->atom_limits.valmax[0] < atom_limits.valmin[0] ||
              rasdef_ptr->atom_limits.valmin[1] > atom_limits.valmax[1] ||
              rasdef_ptr->atom_limits.valmax[1] < atom_limits.valmin[1] ||
              rasdef_ptr->atom_limits.valmin[2] > atom_limits.valmax[2] ||
              rasdef_ptr->atom_limits.valmax[2] < atom_limits.valmin[2])
            include = FALSE;

          /* If not within the limits, then exclude */
          if (include == FALSE)
            rasdef_ptr->deleted = TRUE;
        }

      /* Get next definition record */
      rasdef_ptr = rasdef_ptr->next_rasdef_ptr;
    }
}
/***********************************************************************

write_secstr  -  Write out the PROMOTIF data as HELIX and SHEET records

***********************************************************************/

void write_secstr(FILE *file_ptr,struct promotif *first_promotif_ptr)
{
  char href[4], sref[4];

  int nhelix, nstrand;

  struct promotif *promotif_ptr;

  /* Initialise variables */
  nhelix = 0;
  nstrand = 0;

  /* Get pointer to the first of the PROMOTIF definitions */
  promotif_ptr = first_promotif_ptr;

  /* Loop through the PROMOTIF definitions to write out the helices */
  while (promotif_ptr != NULL)
    {
      /* If a helix, write out the HELIX line */
      if (promotif_ptr->sec_struc == 'H')
        {
          /* Increment helix count */
          nhelix++;

          /* Form the helix identifier */
          sprintf(href,"%d",nhelix);
          if (nhelix < 10)
            strcat(href,"  ");
          else if (nhelix < 100)
            strcat(href," ");

          /* Write out the helix start and end residues */
          fprintf(file_ptr,"HELIX  %3d %s ",nhelix,href);
          fprintf(file_ptr,"%s %c %s ",promotif_ptr->res_name1,
                  promotif_ptr->chain1,promotif_ptr->res_num1);
          fprintf(file_ptr,"%s %c %s",promotif_ptr->res_name2,
                  promotif_ptr->chain2,promotif_ptr->res_num2);
          fprintf(file_ptr,"%2d\n",promotif_ptr->class);
        }

      /* Get pointer to the next PROMOTIF data item */
      promotif_ptr = promotif_ptr->next_promotif_ptr;
    }

  /* Get pointer to the first of the PROMOTIF definitions again */
  promotif_ptr = first_promotif_ptr;

  /* Loop through the PROMOTIF definitions to write out the strands */
  while (promotif_ptr != NULL)
    {
      /* If a strand, write out the SHEET line */
      if (promotif_ptr->sec_struc == 'E')
        {
          /* Increment strand count */
          nstrand++;

          /* Form the strand identifier */
          sprintf(sref,"%d",nstrand);
          if (nstrand < 10)
            strcat(sref,"  ");
          else if (nstrand < 100)
            strcat(sref," ");

          /* Write out the strand start and end residues */
          fprintf(file_ptr,"SHEET %3d %s 1 ",nstrand,sref);
          fprintf(file_ptr,"%s %c%s ",promotif_ptr->res_name1,
                  promotif_ptr->chain1,promotif_ptr->res_num1);
          fprintf(file_ptr,"%s %c%s ",promotif_ptr->res_name2,
                  promotif_ptr->chain2,promotif_ptr->res_num2);
          fprintf(file_ptr,"%2d\n",0);
        }

      /* Get pointer to the next PROMOTIF data item */
      promotif_ptr = promotif_ptr->next_promotif_ptr;
    }
}
/***********************************************************************

write_atom_records  -  Write out all the stored ATOM records

***********************************************************************/

int write_atom_records(FILE *file_ptr,struct atom *first_atom_ptr,
                       int *have_atom,struct parameters params)
{
  char atom_name[5], chain;

  int atom_number, type;
  int wanted;

  struct atom *atom_ptr;
  struct rasdef *rasdef_ptr, *last_rasdef_ptr;
  struct residue *residue_ptr;

  /* Initialise variables */
  atom_number = 0;
  last_rasdef_ptr = NULL;

  /* Initialise pointer to first atom */
  atom_ptr = first_atom_ptr;
          
  /* Loop through all the atoms to write out */
  while (atom_ptr != NULL)
    {
      /* Get this atom's residue record */
      residue_ptr = atom_ptr->residue_ptr;

      /* Assume atom is wanted */
      wanted = TRUE;
      type = UNKNOWN;

      /* Check for deleted ligands and metals */
      rasdef_ptr = residue_ptr->rasdef_ptr;
      if (rasdef_ptr != NULL)
        {
          if (rasdef_ptr->deleted == TRUE)
            wanted = FALSE;
          else
            {
              /* Get the object type */
              type = rasdef_ptr->type;
              
              /* If have a change of object, throw a TER record */
              if (rasdef_ptr != last_rasdef_ptr)
                {
                  /* Write out a TER record */
                  if (last_rasdef_ptr != NULL)
                    fprintf(file_ptr,"TER\n");

                  /* Save the object */
                  last_rasdef_ptr = rasdef_ptr;
                }
            }
        }

      /* If atom wanted, write out */
      if (wanted == TRUE)
        {
          /* Get this atom type and chain-id */
          strcpy(atom_name,atom_ptr->atom_name);
          chain = residue_ptr->chain;

          /* If writing for MolScript and have a calcium, need to change
             its name to AM to get round MolScript bug */
          if ((params.output_option == OUT_MOLSCRIPT ||
              params.output_option == OUT_RASTER3D) &&
              !strncmp(atom_name,"CA",2))
            {
              atom_name[0] = 'A';
              atom_name[1] = 'M';
            }

          /* Write whether this is an ATOM or HETATM record */
          if (type == LIGAND || type == METAL)
            fprintf(file_ptr,"HETATM");
          else
            fprintf(file_ptr,"ATOM  ");

          /* Get the atom number */
          atom_number = atom_ptr->atom_number;

          /* Write out the details */
          fprintf(file_ptr,"%5d %s %s %c%s   %8.3f%8.3f%8.3f%6.2f%6.2f\n",
                  atom_number,atom_name,residue_ptr->res_name,chain,
                  residue_ptr->res_num,atom_ptr->x,atom_ptr->y,atom_ptr->z,
              atom_ptr->occupancy,atom_ptr->bvalue);

          /* Set flag indicating atom has been written out */
          if (atom_ptr->atom_number <= MAXATOMS)
            have_atom[atom_ptr->atom_number - 1] = TRUE;
        }

      /* Get pointer to next atom in linked-list */
      atom_ptr = atom_ptr->next_atom_ptr;
    }

  /* Throw a TER separator at the end */
  fprintf(file_ptr,"TER\n");

  /* Return the last atom number */
  return(atom_number);
}
/***********************************************************************

write_conect_records  -  Write out all the stored CONECT records

***********************************************************************/

void write_conect_records(FILE *outfile_ptr,
                          struct conect *first_conect_ptr,
                          int *have_atom)
{
  int iatom, iconect, nconect, nstore, store_number[MAXCONECT];

  struct conect *conect_ptr;

  /* Initialise pointer to first CONECT record */
  conect_ptr = first_conect_ptr;

  /* Loop through all the CONECT records, writing each one out */
  while (conect_ptr != NULL)
    {
      /* Initialise for this CONECT record */
      nconect = conect_ptr->nconect;
      nstore = 0;

      /* Get first atom number */
      iatom = conect_ptr->atom_number[0];

      /* Check that it is valid */
      if (have_atom[iatom] == TRUE)
        {
          /* Store the first atom number */
          store_number[nstore] = iatom;
          nstore++;

          /* Check that at least one of the atoms the current atom
             is connected to also exist */
          for (iconect = 1; iconect < nconect; iconect++)
            {
              /* Get this atom number */
              iatom = conect_ptr->atom_number[iconect];

              /* If it has been written out, then store its number */
              if (have_atom[iatom] == TRUE)
                {
                  /* Store this atom number */
                  store_number[nstore] = iatom;
                  nstore++;
                }
            }

          /* If have at least one atom connected to the current one,
             then write the CONECT record */
          if (nstore > 1)
            {
              /* Start this record */
              fprintf(outfile_ptr,"CONECT");

              /* Loop through all the atom numbers in this record */
              for (iconect = 0; iconect < nstore; iconect++)
                {
                  /* Write out */
                  fprintf(outfile_ptr,"%5d",store_number[iconect]);
                }

              /* Finish this record */
              fprintf(outfile_ptr,"\n");
            }
        }

      /* Get pointer to next conect in linked-list */
      conect_ptr = conect_ptr->next_conect_ptr;
    }
}
/***********************************************************************

write_moline  -  Write out all the current MolScript secondary structure
                 definition

***********************************************************************/

void write_moline(FILE *outfile_ptr,int *last_sec_struc,
                  char *first_chain,char *first_res_num,
                  char *last_chain,char *last_res_num,char *prev_chain,
                  struct parameters params)
{
  char from_res[7], to_res[7], sec_struc_name[15];
  char colour_name[COLNAM_LEN], colour_numbers[COLNO_LEN];

  int i, icol, ipos;

  /* If colouring by chain, see if colour needs changing */
  if (params.colour_option == COLOUR_BY_CHAIN &&
      *last_chain != *prev_chain)
    {
      /* Chain has changed, so get colour of new one */
      get_chain_col(*last_chain,&icol,colour_name,colour_numbers);

      /* Define colour */
      fprintf(outfile_ptr,"  set planecolour %s;\n",colour_name);

      /* Save current chain */
      *prev_chain = *last_chain;
    }

  /* Determine last secondary structure and define colour */

  /* Coil */
  if (*last_sec_struc == COIL)
    {
      /* Store sec struc name */
      strcpy(sec_struc_name,"coil");
      if (params.colour_option == COLOUR_BY_SST)
        fprintf(outfile_ptr,"  set planecolour %s;\n",
                params.coil_colour);
    }
  
  /* Helix */
  else if (*last_sec_struc == HELIX)
    {
      /* Store sec struc name */
      strcpy(sec_struc_name,"helix");
      if (params.colour_option == COLOUR_BY_SST)
        fprintf(outfile_ptr,"  set planecolour %s;\n",
                params.helix_colour);
    }
  
  /* Strand */
  else if (*last_sec_struc == STRAND)
    {
      /* Store sec struc name */
      strcpy(sec_struc_name,"strand");
      if (params.colour_option == COLOUR_BY_SST)
    fprintf(outfile_ptr,"  set planecolour %s;\n",
        params.strand_colour);
    }

  /* Double helix */
  else if (*last_sec_struc == DHELIX)
    {
      /* Store sec struc name */
      strcpy(sec_struc_name,"double-helix");
      if (params.colour_option == COLOUR_BY_SST)
        fprintf(outfile_ptr,"  set planecolour %s;\n",
                params.dhelix_colour);
    }
  
  /* Form the output residue name for start of range */
  from_res[0] = *first_chain;
  from_res[1] = '\0';
  ipos = 1;
  for (i = 0; i < (int) strlen(first_res_num); i++)
    {
      if (first_res_num[i] != ' ')
        {
          from_res[ipos] = first_res_num[i];
          ipos++;
        }
    }
  from_res[ipos] = '\0';

  /* Repeat for end of range */
  to_res[0] = *last_chain;
  to_res[1] = '\0';
  ipos = 1;
  for (i = 0; i < (int) strlen(last_res_num); i++)
    {
      if (last_res_num[i] != ' ')
        {
          to_res[ipos] = last_res_num[i];
          ipos++;
        }
    }
  to_res[ipos] = '\0';

  /* Write to the MolScript input file unless start and end residues
     are identical */
  if (strcmp(from_res,to_res))
    {
      /* For double-helix, show bonds */
      if (!strcmp(sec_struc_name,"double-helix"))
        {
          fprintf(outfile_ptr,"  set atomradius in from %s to %s 0.2;\n",
                  from_res,to_res);
          fprintf(outfile_ptr,"  ball-and-stick in from %s to %s;\n",
                  from_res,to_res);
        }

      /* Print out the secondary structure type */
      fprintf(outfile_ptr,"  %s from %s to %s;\n",sec_struc_name,from_res,
              to_res);
    }

  /* Re-initialise variables */
  *last_sec_struc = COIL;
  *first_chain = ' ';
  first_res_num[0] = '\0';
}
/***********************************************************************

write_molscript_metals  -  Write out all the metal definitions to the
                           MolScript file

***********************************************************************/

void write_molscript_metals(FILE *outfile_ptr,
                            struct rasdef *first_rasdef_ptr,
                            struct parameters params)
{
  char element[5];

  double rgb[3];

  struct rasdef *rasdef_ptr;

  /* Get pointer to the first RasMol definition record */
  rasdef_ptr = first_rasdef_ptr;

  /* Loop through all the definition records */
  while (rasdef_ptr != NULL)
    {
      /* Check for a match to a metal (ignoring Se) */
      if (rasdef_ptr->type == METAL &&
          strncmp(rasdef_ptr->res_num,"SE",2))
        {
          /* Get metal's element name and truncate */
          strcpy(element,rasdef_ptr->res_num);
          string_truncate(element,4);

          /* Get this metal's colour */
          get_colour_rgb(rasdef_ptr->colour,rgb);

          /* If have a calcium, then need to get round MolScript bug
             by renaming to "AM" */
          if (!strncmp(element,"CA",2))
            strcpy(element,"AM  ");

          /* Write out atom radius for metals */
          fprintf(outfile_ptr,"  set atomradius element %s %4.2f;\n",
                  element,params.metatom_radius);

          /* Write out the colour for this metal */
          fprintf(outfile_ptr,"  set atomcolour element %s ",element);
          fprintf(outfile_ptr,"rgb %8.4f %8.4f %8.4f;\n",rgb[0],
                  rgb[1],rgb[2]);

          /* Write out as sphere */
          fprintf(outfile_ptr,"  cpk element %s;\n",element);
        }

      /* Get next definition record */
      rasdef_ptr = rasdef_ptr->next_rasdef_ptr;
    }
}
/***********************************************************************

write_molscript_ligands  -  Write out all the ligand definitions to the
                            MolScript file

***********************************************************************/

void write_molscript_ligands(FILE *outfile_ptr,
                             struct rasdef *first_rasdef_ptr,
                             struct parameters params)
{
  /* Standard atom types and radii */
#define NATOM_TYPES           6
  static char *ELEMENT[NATOM_TYPES] = {
    "XX", "N", "C", "O", "S", "P"
  };

  int iatom;

  double rgb[3];

  static double VDW_RADIUS[NATOM_TYPES] = {
    1.80,  /*  0. Unknown */
    1.65,  /*  1. Nitrogen */
    1.87,  /*  2. Carbon */
    1.40,  /*  3. Oxygen */
    1.85,  /*  4. Sulphur */
    1.90,  /*  5. Phosphorus */
  };

  struct rasdef *rasdef_ptr;

  /* Define bond radius */
  fprintf(outfile_ptr,"  set stickradius %4.2f;\n",params.ligbond_radius);

  /* Define bond colour */
  if (!strcmp(params.ligand_colour,"CPK"))
    {
      /* Colour by atom colour */
      fprintf(outfile_ptr,"  set colourparts on;\n");
    }

  /* Otherwise, use user-defined colour */
  else
    {
      /* Get RGB values for ligand colour */
      get_colour_rgb(params.ligand_colour,rgb);

      /* Write out the colour */
      fprintf(outfile_ptr,"  set planecolour rgb %8.4f %8.4f %8.4f;\n",
              rgb[0],rgb[1],rgb[2]);
    }

  /* Get pointer to the first RasMol definition record */
  rasdef_ptr = first_rasdef_ptr;

  /* Loop through all the definition records */
  while (rasdef_ptr != NULL)
    {
      /* Check for a ligand range */
      if (rasdef_ptr->type == LIGAND_RANGE)
        {
          /* Define command for single residue */
          if (!strcmp(rasdef_ptr->ligand_from,rasdef_ptr->ligand_to))
            {
              /* If have a fixed atom size, then define it for this residue */
              if (params.ligatom_radius >= 0.0)
                fprintf(outfile_ptr,
                        "  set atomradius in residue %s  %4.2f;\n",
                        rasdef_ptr->ligand_from,
                        4 * params.ligatom_radius);

              /* Otherwise, define van der Waals radii */
              else
                {
                  /* van der Walls radii (times 4 as MolScript reduces
                     ball-and-stick atom radii by a factor of 4) */
                  for (iatom = 1; iatom < NATOM_TYPES; iatom++)
                    fprintf(outfile_ptr,
                            "  set atomradius res-atom %s %s* %4.2f;\n",
                            rasdef_ptr->ligand_from,ELEMENT[iatom],
                            4 * VDW_RADIUS[iatom]);
                }

              /* If atom colours not CPK, then define ligand colour */
              if (strcmp(params.ligand_colour,"CPK"))
                {
                  fprintf(outfile_ptr,
                          "  set atomcolour in residue %s "        
                          "rgb %8.4f %8.4f %8.4f;\n",
                          rasdef_ptr->ligand_from,rgb[0],rgb[1],rgb[2]);
                }

              /* Define ligand as ball-and-stick */
              fprintf(outfile_ptr,
                      "  ball-and-stick in residue %s;\n",
                      rasdef_ptr->ligand_from);
            }

          /* Commands for range of residues */
          else
            {
              /* If have a fixed atom size, then define it for this
                 residue range */
              if (params.ligatom_radius >= 0.0)
                fprintf(outfile_ptr,
                        "  set atomradius in from %s to %s  %4.2f;\n",
                        rasdef_ptr->ligand_from,rasdef_ptr->ligand_to,
                        4 * params.ligatom_radius);

              /* Otherwise, define van der Waals radii */
              else
                {
                  /* van der Walls radii (times 4 as MolScript reduces
                     ball-and-stick atom radii by a factor of 4) */
                  for (iatom = 1; iatom < NATOM_TYPES; iatom++)
                    fprintf(outfile_ptr,
                            "  set atomradius require in from %s to %s "
                            "and atom %s* %4.2f;\n",
                            rasdef_ptr->ligand_from,rasdef_ptr->ligand_to,
                            ELEMENT[iatom],4 * VDW_RADIUS[iatom]);
                }

              /* If atom colours not CPK, then define ligand colour */
              if (strcmp(params.ligand_colour,"CPK"))
                {
                  fprintf(outfile_ptr,
                          "  set atomcolour in from %s to  %s "
                          "rgb %8.4f %8.4f %8.4f;\n",
                          rasdef_ptr->ligand_from,rasdef_ptr->ligand_to,
                          rgb[0],rgb[1],rgb[2]);
                }

              /* Define ligand as ball-and-stick */
              fprintf(outfile_ptr,
                      "  ball-and-stick in from %s to %s;\n",
                      rasdef_ptr->ligand_from,rasdef_ptr->ligand_to);
            }
        }

      /* Get next definition record */
      rasdef_ptr = rasdef_ptr->next_rasdef_ptr;
    }
}
/***********************************************************************

write_molscript  -  Write out all the MolScript commands to the
                    molscript.in file

***********************************************************************/

void write_molscript(FILE *outfile_ptr,char *input_file,
                     struct residue *first_residue_ptr,
                     struct rasdef *first_rasdef_ptr,
                     struct parameters params)
{
  char first_chain, last_chain, prev_chain;
  char first_res_num[6], last_res_num[6];
  char colour_numbers[COLNO_LEN];

  int i, icol[3];
  int last_sec_struc;

  float rgb[3];

  struct residue *residue_ptr;

  /* Initialise variables */
  for (i = 0; i < 3; i++)
    rgb[i] = (float) 0.0;

  /* Write out the header records to the molscript.in file */
  fprintf(outfile_ptr,"! File: molscript.in\n");
  fprintf(outfile_ptr,"! Creator: rascript - Roman Laskowski Jan 2002\n");
  fprintf(outfile_ptr,"! Version: MolScript v1.3\n");
  fprintf(outfile_ptr,"\n");
  fprintf(outfile_ptr,"plot\n");

  /* Get the corresponding RasMol colour triplet */
  get_rasmol_colour(params.background_colour,colour_numbers);

  /* Get the rgb equivalent of the background colour */
  sscanf(colour_numbers,"[%d,%d,%d]",&icol[0],&icol[1],&icol[2]);
  for (i = 0; i < 3; i++)
    rgb[i] = (float) icol[i] / (float) 255.0;

  /* Define the background colour */
  fprintf(outfile_ptr,"  background rgb %8.3f %8.3f %8.3f;\n",rgb[0],
          rgb[1],rgb[2]);
  fprintf(outfile_ptr,"\n");

  /* Define the input file */
  fprintf(outfile_ptr,"  read mol \"%s\";\n",input_file);

  /* Write out the rotations */
  fprintf(outfile_ptr,"  transform atom *\n");
  fprintf(outfile_ptr,"    by centre position atom *\n");
  fprintf(outfile_ptr,"    by rotation x 180.0\n");

  /* Check for reverse order rotation */
  if (params.revrot == TRUE)
    {
      fprintf(outfile_ptr,"    by rotation x %8.3f\n",
              params.rotation_angle[0]);
      fprintf(outfile_ptr,"    by rotation y %8.3f\n",
              params.rotation_angle[1]);
      fprintf(outfile_ptr,"    by rotation z %8.3f;\n",
              params.rotation_angle[2]);
      fprintf(outfile_ptr,"\n");
    }

  /* Otherwise, use standard order */
  else
    {
      fprintf(outfile_ptr,"    by rotation z %8.3f\n",
              params.rotation_angle[2]);
      fprintf(outfile_ptr,"    by rotation y %8.3f\n",
              params.rotation_angle[1]);
      fprintf(outfile_ptr,"    by rotation x %8.3f;\n",
              params.rotation_angle[0]);
      fprintf(outfile_ptr,"\n");
    }

  /* Define plane colours */
  fprintf(outfile_ptr,"  set specularcolour white;\n");
  fprintf(outfile_ptr,"  set stickradius 0.2;\n");
  fprintf(outfile_ptr,"\n");
  fprintf(outfile_ptr,"  set plane2colour grey 0.8;\n");
  fprintf(outfile_ptr,"  set planecolour white;\n");
  fprintf(outfile_ptr,"\n");

  /* Initialise secondary structure */
  last_sec_struc = COIL;
  first_chain = ' ';
  first_res_num[0] = '\0';
  last_chain = ' ';
  last_res_num[0] = '\0';
  prev_chain = '\0';

  /* Get pointer to the first residue */
  residue_ptr = first_residue_ptr;

  /* Loop over all the residues to write out secondary structure
     assignments */
  while (residue_ptr != NULL)
    {
      /* If residue belongs to protein, get its secondary structure */
      if (residue_ptr->type == PROTEIN ||
          residue_ptr->type == NUCLEIC_ACID)
        {
          /* If secondary structure has changed then want to write
             previous definition out */
          if ((residue_ptr->sec_struc != last_sec_struc ||
              residue_ptr->chain != last_chain) &&
              first_res_num[0] != '\0')
            {
              /* If previous secondary structure was of the same chain
                 want it to extend to the current secondary structure */
              if (residue_ptr->chain == last_chain)
                {
                  strcpy(last_res_num,residue_ptr->res_num);
                  last_chain = residue_ptr->chain;
                }

              /* Write out the previous definition */
              write_moline(outfile_ptr,&last_sec_struc,&first_chain,
                           first_res_num,&last_chain,last_res_num,
                           &prev_chain,params);
            }

          /* If first residue of range is empty, then fill it */
          if (first_res_num[0] == '\0')
            {
              /* Fill in all the details */
              strcpy(first_res_num,residue_ptr->res_num);
              first_chain = residue_ptr->chain;
              last_sec_struc = residue_ptr->sec_struc;

              /* If current secondary structure is coil, make it
                 start at last valid residue */
              if (residue_ptr->sec_struc == COIL &&
                  last_res_num[0] != '\0' &&
                  residue_ptr->chain == last_chain)
                {
                  strcpy(first_res_num,last_res_num);
                  first_chain = last_chain;
                }
            }

          /* Fill in details for last residue of range */
          strcpy(last_res_num,residue_ptr->res_num);
          last_chain = residue_ptr->chain;
        }

      /* Otherwise, write out any previous definition */
      else if (first_res_num[0] != '\0')
        {
          /* Write out the definition */
          write_moline(outfile_ptr,&last_sec_struc,&first_chain,
                       first_res_num,&last_chain,last_res_num,
                       &prev_chain,params);

          /* Reset first residue of range */
          first_res_num[0] = '\0';
        }

      /* Get pointer to the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
    }

  /* If have any remaining defintions, then write out */
  if (first_res_num[0] != '\0')
    write_moline(outfile_ptr,&last_sec_struc,&first_chain,
                 first_res_num,&last_chain,last_res_num,&prev_chain,
                 params);

  /* Write out any metals as spheres */
  write_molscript_metals(outfile_ptr,first_rasdef_ptr,params);

  /* Write out ligands */
  write_molscript_ligands(outfile_ptr,first_rasdef_ptr,params);

  /* Close of the script */
  fprintf(outfile_ptr,"end_plot\n");

  /* Close the file */
  fclose(outfile_ptr);
}
/***********************************************************************

get_mask_vertices  -  Read through the .vtx file, created by SURFNET and
                      SURFACE, giving the vertex coordinates of the
                      binding-site mask. Write each vertex out as a 
                      HETATM record in PDB format for display by RasMol

***********************************************************************/

void get_mask_vertices(FILE *outfile_ptr,char vtx_name[FILENAME_LEN],
                       int atom_number)
{
  char input_line[LINELEN + 1];

  float x, y, z;

  FILE *file_ptr;

  /* Open the .vtx file */
  if ((file_ptr = fopen(vtx_name,"r")) ==  NULL)
    {
      printf("\n*** Unable to open the .vtx vertices file [%s]\n",vtx_name);
      exit(1);
    }

  /* Loop through the .vtx file, picking up all the vertex coords  */
  while (fgets(input_line,LINELEN,file_ptr)!= NULL)
    {
      /* Get the x-, y- and z-coordinates */
      sscanf(input_line,"%f %f %f",&x,&y,&z);

      /* Increment atom number */
      atom_number++;

      /* Write out in PDB format */
      if (atom_number < 100000)
        {
          fprintf(outfile_ptr,"HETATM%5d  C   MAP 1   1    ",atom_number);
          fprintf(outfile_ptr,"%8.3f%8.3f%8.3f  1.00 10.00\n",x,y,z);
        }
    }

  /* Close mask vertices file */
  fclose(file_ptr);
}
/***********************************************************************

                             M   A   I   N

***********************************************************************/

int main(int argc,char *argv[])
{
  char chains[MAX_CHAINS], domsst[MAX_CHAINS];
  char out_name[FILENAME_LEN], vtx_name[FILENAME_LEN];
  char grow_name[FILENAME_LEN], mol_name[FILENAME_LEN];
  char pdb_code[5], pdb_name[FILENAME_LEN], pdbsum_dir[FILENAME_LEN];
  char cath_domains_file[FILENAME_LEN], dfn_name[FILENAME_LEN];
  char promotif_dir[FILENAME_LEN], work_dir[FILENAME_LEN];
  char first_object[OBJ_NAME_LEN], first_display_mode[DISPLAY_MODE_LEN];

  int atom_number, i, iatom, nchains;
  int have_dfn, have_grow, have_sstdat;
  int have_atom[MAXATOMS];

  struct c_file_data file_name_ptr;
  struct parameters params;

  struct atom *first_atom_ptr;
  struct conect *first_conect_ptr;
  struct grow *first_grow_ptr;
  struct limits atom_limits;
  struct promotif *first_promotif_ptr;
  struct rasdef *first_rasdef_ptr;
  struct residue *first_residue_ptr;

  FILE *molfile_ptr, *outfile_ptr;

  /* Initialise variables */
  atom_number = 0;
  for (iatom = 0; iatom < MAXATOMS; iatom++)
    have_atom[iatom] = FALSE;
  for (i = 0; i < MAX_CHAINS; i++)
    domsst[i] = ' ';
  first_display_mode[0] = '\0';
  first_object[0] = '\0';
  have_dfn = FALSE;
  have_grow = FALSE;
  have_sstdat = FALSE;
  out_name[0] = '\0';
  strcpy(dfn_name,"rasmol.dfn");
  strcpy(mol_name,"molscript.in");
  strcpy(grow_name,"grow.out");

  /* Initialise pointers */
  first_grow_ptr = NULL;
  first_promotif_ptr = NULL;
  first_rasdef_ptr = NULL;
  first_residue_ptr = NULL;

  /* Interpret the command-line parameters */
  get_command_arguments(argv,argc - 1,pdb_code,chains,&nchains,
                        pdb_name,work_dir,&params);

  /* If writing to file, then set file name */
  if (params.write_to_file == TRUE)
    strcpy(out_name,"rascript.out");

  /* Read in all the directory paths from the parameter file, CATHPARAM */
  c_get_paths(&file_name_ptr,params.run_64);

  /* Get directories holding PDBsum and PDB data from CATHPARAM file */
  form_filenames(pdb_code,pdbsum_dir,pdb_name,vtx_name,
                 promotif_dir,cath_domains_file,file_name_ptr,params);

  /* Open the output file */
  outfile_ptr = open_output_file(out_name,TRUE);

  /* Read in the RasMol definitions and generate the RasMol script */
  if (params.coords_only == FALSE)
    {
      /* Read in the RasMol definitions */
      have_dfn
        = get_dfn_data(pdb_code,pdbsum_dir,dfn_name,&first_rasdef_ptr,
                       chains,&nchains,work_dir,params);

      /* Read in the secondary structure counts from the sst.dat file */
      have_sstdat
        = get_sst_dat(pdb_code,pdbsum_dir,chains,domsst,nchains,params);

      /* If called from Web, read in secondary structures from
         local file.rin file */
      if (work_dir[0] != '\0')
        read_rin_file(&first_promotif_ptr,work_dir);

      /* Otherwise, read in the PROMOTIF secondary structure
         assignments */
      else
        get_promotif_data(pdb_code,&first_promotif_ptr,promotif_dir,
                          chains,nchains);

      /* Get the ligand- and metal-interaction data from the Grow files */
      have_grow
        = read_in_grow_data(pdb_code,pdbsum_dir,work_dir,grow_name,
                            &first_grow_ptr);

      /* For RasMol script output write out all the script commands */
      if (params.output_option == OUT_RASMOL)
        {
          /* Write out the script header records */
          write_script_headers(outfile_ptr,pdb_code,pdb_name,chains,
                               nchains,have_sstdat,pdbsum_dir,params);

          /* If colouring by domain, get the domain info from the CATH domains
             file and colour accordingly */
          if (params.colour_by_domain == TRUE)
            write_script_domains(outfile_ptr,pdb_code,cath_domains_file,
                                 first_object,first_display_mode,
                                 chains,nchains,params);

          /* Generate the script commands from the data in the .dfn RasMol
             definitions file */
          generate_script(outfile_ptr,pdbsum_dir,dfn_name,work_dir,
                          first_object,first_display_mode,chains,nchains,
                          domsst,params);

          /* Write closing lines to script */
          write_script_tail(outfile_ptr,first_object,first_display_mode,
                            params);
        }
    }

  /* Read through the PDB file writing out all the required ATOM and HETATM
     records and associated CONECT records */
  if (params.pdb_type != PDBSUM1)
    {
      get_atom_records(pdb_name,first_rasdef_ptr,first_promotif_ptr,
                       first_grow_ptr,&first_atom_ptr,&first_residue_ptr,
                       chains,nchains,have_dfn,have_grow,&first_conect_ptr,
                       &atom_limits,params);

      /* If no grow data read in, and not dealing with all chains,
         check which ligands and metals are in the proximity of the
         selected chain(s) */
      if (have_grow == FALSE && params.all_chains == FALSE)
        get_interacting_ligands(first_rasdef_ptr,atom_limits);
      
      /* If have PROMOTIF details, write out as HELIX and SHEET records */
      // write_secstr(outfile_ptr,first_promotif_ptr);

      /* Write out the atom coordinate records */
      atom_number
        = write_atom_records(outfile_ptr,first_atom_ptr,have_atom,params);

      /* Write out any relevant CONECT records */
      write_conect_records(outfile_ptr,first_conect_ptr,have_atom);
    }

  /* If MolScript script required, then generate it */
  if (params.output_option == OUT_MOLSCRIPT ||
      params.output_option == OUT_RASTER3D)
    {
      /* Open the MolScript output file */
      molfile_ptr = open_output_file(mol_name,TRUE);

      /* Write out all the script commands */
      write_molscript(molfile_ptr,out_name,first_residue_ptr,
                      first_rasdef_ptr,params);
    }

  /* If adding binding-site mask vertices, as generated by SURFNET, then
     read in their coordinates and write out as PDB HETATM records */
  if (params.show_surfaces == TRUE)
    get_mask_vertices(outfile_ptr,vtx_name,atom_number);

  return(0);
}
