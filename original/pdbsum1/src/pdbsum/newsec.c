/***********************************************************************

  newsec.c - Replacement for secsurf.f

             Program determines the contents of the given PDB file,
             writing the same output files as secsurf, as well as
             identifying unique and copy protein chains (as per chuniq)
             and determining which are the unique ligands in the
             structure (as per liguniq)

Date:-         15 Jan 2004
Dir update:-   16 Oct 2012

Written by:-   Roman Laskowski

----------------------------------------------------------------------*/

/* Compilation
   -----------

cc -c newsec.c
cc -c common.c
cc -c reapar.c
cc -o newsec newsec.o common.o reapar.o -lm -lz

Optimization flags:-

cc -O2 -funroll-loops -c newsec.c
cc -O2 -funroll-loops -c common.c
cc -O2 -funroll-loops -c reapar.c
cc -o newsec newsec.o common.o reapar.o -lm -lz



Datafiles
---------

Actual name                Description
-----------                -----------
AF_scores.dat              Output Alpha Fold confidence scores for
                           each residue
AF_ranges.dat              Output Alpha Fold ranges for colouring
                           structure
CATHPARAM                  Input parameter file
chains.dat                 Output chains file listing unique and
                           duplicate protein chains
dimhtmlXY                  Output file of residue interactions
                           across a P-P or P-DNA interface, used by
                           program dimhtml to generate in interface
                           plot
grow.out                   Output file of H-bonds and non-bonded
                           interactions between protein and ligand,
                           metal and DNA. Is a summary of the
                           key data in the .hhb and .nnb files
igrow.out                  Output file of H-bonds and non-bonded
                           interactions between different protein
                           chains
hbadd.sum                  Input file defining all the Het Groups,
                           as defined by the hbadd program
hetXXX.en                  Input file giving current "energy" score
                           of JPEG associated with Het Group XXX
hetXXX.jpg                 Output image of Het Group XXX as rendered
                           by Raster3D
interfaces.dat             Output file for use by showmain.c
                           to define all unique interfaces and
                           their duplicates
lig_MM_NN.pdb              Output PDB file of each unique
                           ligand-protein instance
ligands.dat                Output file for use by showmain.c
                           to define all unique ligands and
                           their duplicates
newsec.out                 Output file giving all the molecule
                           definitions
pdbXXXX.hhb                Input .hhb file, as created by HBPLUS
                           listing all polar H-bonds in the
                           structure
pdbXXXX.nnb                Input .nnb file, as created by HBPLUS
                           listing all non-bonded contacts in
                           the structure
<pdb_code>_<chain_id>.csq  Output file giving a compressed
                           representation of the residue numbers
                           for use by SASHTML
rasmol.dfn                 Output definition file of PDB contents
secsurf.out                Output file mimicking output produced
                           by secsurf.f and required by programs
                           such as pdbhtml.f and grow
seqXXX.fasta               Output FASTA-format file - one for
                           each unique sequence
seqXXX.html                Output HTML file - one for each
                           unique sequence
seqXXX.trn                 Output file of the residue numbers
                           corresponding to the sequence written
                           to seqXXX.fasta
seqXXX.chain               Output file giving sequence's chain-id
seqXXX.pdb                 Output file giving sequence's coords
surfp.par                  Output file giving surfnet-style info
                           about atom ranges corresponding to
                           each chain
sst.dat                    Counts of residues in each secondary
                           structure type for each chain 
BIOLIST                    PQS file listing all PDB entries and
                           their biological unit
ASALIST                    PQS file listing assessments of PQS calcs
authmol.dat                File containing author-defined biomolecules

*/

/* Function calling-tree
   ---------------------

main
   -> get_command_arguments
   -> c_get_paths (in common.c)
         -> get_param (in reapar.c)
         -> store_string (in common.c)
   -> read_pdb
         -> string_truncate (in common.c)
         -> get_atom_details
         -> is_hydrogen_atom
               -> check_for_hmetal
         -> check_for_duplicate
         -> create_chain_record
         -> create_residue_record
         -> scan_sites
               -> store_string (in common.c)
         -> create_atom_record
               -> get_atom_type
         -> process_conect_record
               -> create_connect_record
                     -> create_bond_record
         -> create_site_record
         -> create_hetdesc_record
         -> store_string (in common.c)
         -> get_qstruc
   -> classify_residues
         -> get_residue_type
         -> get_het_desc
         -> excise_metal
              -> create_chain_record
              -> update_residues
   -> shift_metals
   -> adjust_residues
   -> find_pdbsum_dir (in common.c)
   -> read_rin_file
   -> get_disulphides
         -> string_truncate (in common.c)
         -> get_sulphur
   -> classify_chains
   -> join_chains
   -> reclassify_chains
   -> split_ligands
         -> check_for_bond
               -> find_joins
                     -> create_bond_record
         -> identify_molecule
               -> create_list_record
               -> free_list_memory
         -> create_chain_record
         -> update_residues
   -> single_residue_ligands
         -> check_for_metal
   -> calc_chain_limits
   -> write_defs
   -> check_biomol
         -> check_pqs
               -> string_truncate (in common.c)
               -> get_qstruc
         -> check_asa
               -> string_truncate (in common.c)
         -> check_authmol
               -> string_truncate (in common.c)
               -> get_qstruc
   -> write_biomol_dat
   -> read_hbadd_sum
         -> string_truncate (in common.c)
         -> create_list_record
   -> write_hbadd_files
         -> check_het_group
         -> write_en_file
               -> open_output_file (in common.c)
         -> render_hetgroup
               -> calc_covalent_bonds
               -> open_output_file (in common.c)
               -> write_raster3d_headers
               -> align_molecule
                     -> centre_of_mass
                     -> adjust_to_origin
                     -> principal_components
                           -> calc_eigen_values
                           -> eigen_sort
                           -> orien3
                                 -> vecprd
                                 -> dot3
                     -> rotate_molecule
               -> calc_r3d_scale
               -> write_raster3d_residue
                     -> write_raster3d_atom
                           -> get_atom_colour
                                 -> get_colour_rgb (in common.c)
                     -> write_raster3d_bond
                           -> get_atom_colour
                                 -> get_colour_rgb (in common.c)
         -> run_raster3d
         -> write_coords
               -> open_output_file (in common.c)
   -> extract_sequences
         -> store_string (in common.c)
   -> compare_sequences
         -> perform_comparison
   -> write_chain_equivs
         -> write_seq_fasta
         -> write_seq_html
         -> write_seq_chain
         -> write_seq_pdb
   -> form_ligands
         -> create_ligand_record
         -> add_res_name
         -> format_residue_no
         -> store_string (in common.c)
   -> get_ligand_neighbours
         -> create_list_record
   -> identify_unique_ligands
         -> count_equivs
   -> sort_ligands
         -> shell_sort
   -> write_ligmet_run
   -> write_af_scores
         -> open_output_file (in common.c)
         -> add_range
   -> run_hbplus
   -> create_grow_out
         -> residue_sort
         -> open_output_file (in common.c)
         -> process_hbplus_file
               -> open_output_file (in common.c)
               -> to_upper (in common.c)
               -> reformat_hbplus_line
               -> fix_res_number
               -> binary_residue_search
               -> find_atom
               -> write_grow_line
               -> process_interface
                     -> find_interface
                     -> create_interface_record
                     -> find_respair
                     -> create_respair_record
         -> store_iface_ssbonds
               -> process_interface (as above)
         -> get_salt_bridges
               -> check_residue_pairs
                     -> check_atom_pairs
                     -> process_interface (as above)
   -> save_interface_data
         -> sort_interfaces
               -> interface_sort
         -> unique_interfaces
               -> equiv_chain
               -> count_pairs
         -> write_interfaces_dat
               -> calc_surface_area
                     -> run_naccess
               -> write_details
         -> write_dimhtml_file
               -> calc_iface_transformation
                     -> create_atom_record
                           -> get_atom_type
                     -> align_molecule (as above)
                     -> free_atom_memory
   -> write_pymol_pdb
         -> open_output_file (in common.c)
         -> write_pdb_headers
               -> apply_transformation
               -> get_date (in common.c)
                     -> to_upper (in common.c)
         -> write_pdb_secstruc
         -> write_pdb_coords
   -> write_pymol_pml
         -> open_output_file
         -> remove_residue_blanks
         -> write_pymol_chains
               -> create_pobject_record
         -> write_pymol_ligands
               -> create_pobject_record
         -> write_pymol_metals
               -> create_pobject_record
         -> write_pymol_disulph
               -> create_pobject_record
   -> write_sst_dat
         -> open_output_file (in common.c)
   -> render_structure
         -> open_output_file (in common.c)
         -> write_raster3d_headers
         -> centre_of_mass
         -> adjust_to_origin
         -> rotate_structure
         -> get_max_min
         -> calc_r3d_scale
         -> render_protein_chain
               -> material_properties
               -> get_chain_col
               -> find_sec_struc
                     -> render_sec_struc
                           -> create_atom_record
                                 -> get_atom_type
                           -> generate_helix
                                 -> reverse_transform
                           -> generate_strand
                                 -> store_vertex
                                 -> write_polygon
                                       -> reverse_transform
                           -> modify_calphas
                           -> free_atom_memory
               -> render_backbone
                     -> write_raster3d_bond
                           -> get_atom_colour
                                 -> get_colour_rgb (in common.c)
               -> label_nterminus
                     -> store_vertex
                     -> write_line
         -> render_ca_only_chain
               -> material_properties
               -> get_chain_col
               -> render_backbone
                     -> write_raster3d_bond
                           -> get_atom_colour
                                 -> get_colour_rgb (in common.c)
               -> label_nterminus
                     -> get_colour_rgb (in common.c)
                     -> store_vertex
                     -> write_line
         -> render_dna_chain
               -> material_properties
               -> get_chain_col
               -> render_backbone
                     -> write_raster3d_bond
                           -> get_atom_colour
                                 -> get_colour_rgb (in common.c)
               -> calc_covalent_bonds
               -> write_raster3d_allbonds
                     -> write_raster3d_bond
                           -> get_atom_colour
                                 -> get_colour_rgb (in common.c)
         -> render_ligand_metal
               -> material_properties
               -> write_raster3d_atom
                     -> get_atom_colour
                           -> get_colour_rgb (in common.c)
         -> render_ssbonds
               -> material_properties
         -> render_sites
               -> material_properties
               -> calc_covalent_bonds
               -> write_raster3d_allbonds
                     -> write_raster3d_bond
                           -> get_atom_colour
                                 -> get_colour_rgb (in common.c)
         -> run_raster3d
         -> pymol_render
         -> create_pymol_record
   -> render_all_dna
         -> transparent_image
               -> open_output_file (in common.c)
               -> material_properties
               -> get_tokens (in common.c)
               -> run_raster3d
               -> select_pymol
               -> pymol_render
               -> create_pymol_record
               -> free_tokens (in common.c)
         -> open_output_file (in common.c)
         -> write_raster3d_headers
         -> centre_of_mass
         -> adjust_to_origin
         -> get_max_min
         -> calc_r3d_scale
         -> render_dna_chain (as above)
   -> render_indiv_chains
         -> open_output_file (in common.c)
         -> write_raster3d_headers
         -> centre_of_mass
         -> adjust_to_origin
         -> get_max_min
         -> calc_r3d_scale
         -> render_protein_chain (as above)
         -> render_ca_only_chain (as above)
         -> render_ssbonds
         -> render_sites (as above)
         -> run_raster3d
         -> transparent_image (as above)
   -> run_pymol
         -> call_system (in common.c)
   -> convert_png_images
   -> write_dat_file
         -> open_output_file (in common.c)
         -> split_residue_sequence
         -> write_het_groups
   -> render_ligands
         -> calc_covalent_bonds
         -> open_output_file (in common.c)
         -> write_raster3d_headers
         -> align_molecule (as above)
         -> calc_r3d_scale
         -> write_raster3d_ligand
               -> write_raster3d_atom
                     -> get_atom_colour
                           -> get_colour_rgb (in common.c)
               -> write_raster3d_bond
                     -> get_atom_colour
                           -> get_colour_rgb (in common.c)
         -> render_pymol_ligand
               -> open_output_file (in common.c)
               -> run_pymol
         -> run_raster3d
         -> transparent_image (as above)
         -> write_run_scr_params
               -> store_string (in common.c)
               -> append_string (in common.c)
               -> to_upper
   -> render_interfaces
         -> render_chain_pair
               -> open_output_file (in common.c)
               -> write_raster3d_headers
               -> centre_of_mass2
               -> adjust_to_origin
               -> get_max_min2
               -> calc_r3d_scale
               -> render_protein_chain
               -> run_raster3d
         -> transparent_image (as above)
   -> write_secsurf_out
         -> open_output_file (in common.c)
         -> write_secsurf_chain
               -> get_chain_colname
               -> write_secsurf_sequence
                     -> create_hetgroup_entry
                     -> free_hetgroup_memory
         -> write_surf_defs
               -> get_chain_colname
   -> write_surfp_par
         -> write_surf_defs
               -> get_chain_colname
   -> write_rasmol_dfn
         -> open_output_file (in common.c)
   -> get_biomol_pdb
         -> open_output_file (in common.c)
         -> string_truncate (in common.c)
         -> get_qstruc
         -> get_atom_details
         -> is_hydrogen_atom
               -> check_for_hmetal
         -> get_residue_type
         -> find_chain
         -> create_chain_record
         -> find_residue
         -> create_residue_record
         -> find_atom
         -> create_atom_record
   -> mirror_molecules
         -> create_chain_record
   -> form_ligands
   -> copy_disulphides
   -> align_molecule (as above)
   -> write_pymol_pdb(as above)
   -> write_biomol_dat
         -> open_output_file

*/

#include "common.h"
#include "readpdb.h"
#include "newsec.h"


/* Prototypes
   ---------- */

/* In common.c */
void append_string(char **string_ptr,char *string);
void call_system(char *command,int dos);
void c_get_paths(struct c_file_data *file_name_ptr,int run_64);
int find_pdbsum_dir(char *pdb_code,char *pdbsum_template[NPDBSUM_DIR],
                    int pdb_type,char *dir_name);
char *form_filename(char *pdb_code,char *filename_template);
void free_tokens(char **token,int ntokens);
void get_colour_rgb(char *c_name,double rgb[3]);
void get_date(char *date);
int get_tokens(char *line,char **token,int max_token,char token_sep);
FILE *open_output_file(char *out_name,int terminate);
void store_string(char **string_ptr,char *string);
int string_truncate(char *string,int max_length);
void to_upper(char *search_string,int length);

/* Function prototypes from reapar.c */
int get_param(char *token_name,char *file_name,int filename_len,
              int exit_on_error);


/***********************************************************************

get_command_arguments  -  Interpret the command line arguments

***********************************************************************/

void get_command_arguments(char *strptrs[],int ntoken,
                           struct parameters *params)
{
  int itoken;

  /* Initialise parameters */
  params->file_name[0] = '\0';
  params->pdb_code[0] = '\0';
  params->no_gifs = FALSE;
  params->uploaded_file = FALSE;
  params->pdb_type = STANDARD_PDB;
  params->calc_interfaces = TRUE;
  params->molecule_align = TRUE;
  params->gen_pymol = FALSE;
  params->biomol_only = FALSE;
  params->ligand_length = LIGAND_LENGTH;
  params->run_64 = FALSE;
  params->run_raster3d = RUN_RASTER3D;

  /* If -help entered on the command line, show options available */
  if (ntoken < 2 || (ntoken == 1 && (!strcmp(strptrs[1],"-help") ||
                                     !strcmp(strptrs[1],"-h"))))
    {
      printf("\n");
      printf("Correct usage is ...\n");
      printf("\n");
      printf("    newsec  pdb_code  filename  [-nogifs]  [-uplo]\n");
      printf("  [UPLOAD | MCSG | AF]\n");
      printf("            [-noint]  [-biomol]  [-64]  [-pdbsum1]\n");
      printf("where\n");
      printf("     * pdb_code = PDB code for file to be processed\n");
      printf("     * filename = Name of PDB file to be processed\n");
      printf("     * -nogifs  = Images not required\n");
      printf("     * -uplo    = File is an uploaded PDB file\n");
      printf("     * UPLOAD   = Uploaded PDB file for PDBsum generation\n");
      printf("     * MCSG     = MCSG structure\n");
      printf(      " AF       = Alpha Fold structure\n");
      printf("     * -noalign = skip molecule alignment routines\n");
      printf("     * -noint   = skip interface calculation routines\n");
      printf("     * -pymol   = use PyMol to generate thumbnail images\n");
      printf("     * -biomol  = generate biomolecule images only\n");
      printf("     * -pdbsum1 = PDBsum1 structure\n");
      printf("     * -llen    = Default ligand length\n");
      printf("     * -64      = run on 64-bit processor\n");
      printf("\n");
      exit(1);
    }

  /* Get the PDB code from the first token */
  strcpy(params->pdb_code,strptrs[1]);
  params->pdb_code[4] = '\0';

  /* Get name of PDB file from the second token */
  strcpy(params->file_name,strptrs[2]);

  /* Loop through any additional command arguments */
  for (itoken = 3; itoken < ntoken + 1; itoken++)
    {
      /* Check for the -nogifs option */
      if (!strcmp(strptrs[itoken],"-nogifs"))
        params->no_gifs = TRUE;

      /* Check for the -uplo option */
      else if (!strcmp(strptrs[itoken],"-uplo"))
        params->uploaded_file = TRUE;

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
        {
          params->pdb_type = PDBSUM1;
          params->run_raster3d = FALSE;
        }

      /* Check for the -noint option */
      else if (!strcmp(strptrs[itoken],"-noint"))
        params->calc_interfaces = FALSE;

      /* Check for the -noalign option */
      else if (!strcmp(strptrs[itoken],"-noalign"))
        params->molecule_align = FALSE;

      /* Check for the PyMol option */
      else if (!strcmp(strptrs[itoken],"-pymol"))
        params->gen_pymol = TRUE;

      /* Check for the biomolecules option */
      else if (!strcmp(strptrs[itoken],"-biomol"))
        params->biomol_only = TRUE;

      /* Check for the -llen option giving the default ligand length */
      else if (!strcmp(strptrs[itoken],"-llen"))
        {
          /* Get the length from the next token, if there is one */
          if (itoken < ntoken)
            {
              params->ligand_length = atoi(strptrs[itoken + 1]);
              itoken++;
            }          
        }

      /* Check for the 64-bit option */
      else if (!strcmp(strptrs[itoken],"-64"))
        params->run_64 = TRUE;
    }

  /* If no images required, then don't run PyMol */
  if (params->no_gifs == TRUE)
    params->gen_pymol = FALSE;
}
/***********************************************************************

get_atom_details  -  Extract the atom details from the given line

***********************************************************************/

void get_atom_details(char input_line[LINELEN + 1],int *atom_number,
                      char *atom_name,
                      char *chain,char *res_name,char *res_num,
                      double *xval,double *yval,double *zval,
                      double *bval,double *occup,char *alt,char *element)
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

  /* Get the altyernate occupancy identifier */
  *alt = input_line[16];

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

  /* Retrieve the element, if present */
  if (len > 77)
    {
      strncpy(element,input_line+77,2);
      element[2] = '\0';
    }
  else
    element[0] = '\0';

  /* Store the retrieved values */
  *xval = x;
  *yval = y;
  *zval = z;
  *bval = bvalue;
  *occup = occupancy;
}
/***********************************************************************

check_for_hmetal  -  Check whether the given atom is one of the metal
                     ions whose name begins with H

***********************************************************************/

int check_for_hmetal(char atom_name[5])
{
  int imetal;
  int metal;

  static int nmetals = 4;

  static char *metal_name[] = { "HE  ", "HO  ", "HF  ", "HG  " };

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
        metal = TRUE;
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
  char atom_name[5];

  int hydrogen, metal;

  /* Initialise */
  hydrogen = FALSE;
  strncpy(atom_name,atname,4);
  atom_name[4] = '\0';

  /* Check the atom name for any indicators to suggest this is either
     a hydrogen or a pseudo atom */
  if (atom_name[0] == 'H' || atom_name[1] == 'H' || atom_name[1] == 'Q' ||
      (atom_name[1] == 'D' && (atom_name[0] == ' ' || atom_name[0] == '1' ||
                               atom_name[0] == '2' || atom_name[0] == '3')))
    hydrogen = TRUE;

  /* Check for possible strange hydrogens */
  if (atom_name[0] >= '0' && atom_name[0] <= '9' && atom_name[2] == 'H')
    hydrogen = TRUE;

  /* If this looks likely to be a hydrogen, check that it isn't in
     fact a metal */
  if (hydrogen == TRUE)
    {
      /* Check whether this is a metal */
      metal = check_for_hmetal(atom_name);

      /* If a metal, then isn't a hydrogen */
      if (metal == TRUE)
        hydrogen = FALSE;
    }

  return hydrogen;
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

create_chain_record  -  Create and initialise a new chain record

***********************************************************************/

struct chain *create_chain_record(struct chain **fst_chain_ptr,
                                  struct chain **lst_chain_ptr,
                                  char chain_id)
{
  int i;

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
  chain_ptr->chain_id = chain_id;
  chain_ptr->offset = 0;
  sprintf(chain_ptr->name,"Chain %c",chain_id);
  chain_ptr->sequence = &null;
  chain_ptr->type = UNKNOWN;
  chain_ptr->nresid = 0;
  chain_ptr->nstandard = 0;
  chain_ptr->natoms = 0;
  chain_ptr->first_atom_number = 0;
  chain_ptr->last_atom_number = 0;
  chain_ptr->nhelices = 0;
  chain_ptr->nstrands = 0;
  chain_ptr->nstand = 0;
  chain_ptr->nbase = 0;
  chain_ptr->copy_of = '\0';
  chain_ptr->ncopies = 0;
  chain_ptr->got_bonds = FALSE;
  chain_ptr->got_sites = FALSE;
  chain_ptr->same_length = FALSE;
  chain_ptr->natom_recs = 0;
  chain_ptr->nhetatm_recs = 0;
  chain_ptr->asa = 0.0;
  chain_ptr->dominant_sst = HELIX;
  for (i = 0; i < 3; i++)
    chain_ptr->nres_sst[i] = 0;
  for (i = 0; i < 3; i++)
    {
      chain_ptr->coord_min[i] = 0.0;
      chain_ptr->coord_max[i] = 0.0;
    }
  chain_ptr->symm_chain = FALSE;
  chain_ptr->used = FALSE;
  chain_ptr->link_chain_ptr = NULL;
  chain_ptr->first_residue_ptr = NULL;
  chain_ptr->next_chain_ptr = NULL;

  /* Return current pointers */
  *fst_chain_ptr = first_chain_ptr;
  *lst_chain_ptr = last_chain_ptr;

  /* Return new chain pointer */
  return(chain_ptr);
}
/***********************************************************************

create_residue_record  -  Create and initialise a new residue record

***********************************************************************/

struct residue *create_residue_record(struct residue **fst_residue_ptr,
                                      struct residue **lst_residue_ptr,
                                      struct chain *chain_ptr,
                                      char *res_name,char *res_num,
                                      char chain,int number)
{
  int i;

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
  residue_ptr->number = number;
  residue_ptr->type = UNKNOWN;
  residue_ptr->aa_code = 'X';
  residue_ptr->ref = -1;
  strcpy(residue_ptr->res_name,res_name);
  strcpy(residue_ptr->res_num,res_num);
  residue_ptr->chain = chain;
  residue_ptr->sec_struc = COIL;
  residue_ptr->natoms = 0;
  residue_ptr->got_bonds = FALSE;
  residue_ptr->connected = FALSE;
  residue_ptr->nmain_chain = 0;
  residue_ptr->ncentroid_atoms = 0;
  residue_ptr->niface_atoms = 0;
  residue_ptr->have_calpha = FALSE;
  residue_ptr->range = 0;
  residue_ptr->first_atom_ptr = NULL;
  residue_ptr->first_rbond_ptr = NULL;
  residue_ptr->hetdesc_ptr = NULL;
  residue_ptr->ligand_ptr = NULL;
  residue_ptr->link_residue_ptr = NULL;
  residue_ptr->attach_residue_ptr = NULL;
  residue_ptr->chain_ptr = chain_ptr;
  for (i = 0; i < 3; i++)
    {
      residue_ptr->coord_min[i] = 0.0;
      residue_ptr->coord_max[i] = 0.0;
      residue_ptr->centroid[i] = 0.0;
      residue_ptr->ave_iface_coords[i] = 0.0;
    }
  residue_ptr->site = &null;
  residue_ptr->disulphide = FALSE;

  /* If this is the first residue for this chain, store its pointer */
  if (chain_ptr->first_residue_ptr == NULL)
    chain_ptr->first_residue_ptr = residue_ptr;

  /* Increment count of this chain's residues */
  chain_ptr->nresid++;

  /* Initialise pointer to the next residue */
  residue_ptr->next_residue_ptr = NULL;

  /* Return current pointers */
  *fst_residue_ptr = first_residue_ptr;
  *lst_residue_ptr = last_residue_ptr;

  /* Return new residue pointer */
  return(residue_ptr);
}
/***********************************************************************

scan_sites  -  Determine whether given residue is part of any active
               site

***********************************************************************/

int scan_sites(struct residue *residue_ptr,struct site *first_site_ptr)
{
  int nsite_res;
  int done;

  struct site *site_ptr;

  /* Initialise variables */
  done = FALSE;
  nsite_res = 0;

  /* Get pointer to the first site */
  site_ptr = first_site_ptr;

  /* Loop over all the sites */
  while (site_ptr != NULL && done == FALSE)
    {
      /* Check residue name, number and chain id */
      if (!strcmp(site_ptr->res_num,residue_ptr->res_num) &&
          !strcmp(site_ptr->res_name,residue_ptr->res_name) &&
          site_ptr->chain == residue_ptr->chain)
        {
          /* Store the site name */
          store_string(&(residue_ptr->site),site_ptr->name);

          /* Increment count of site residues */
          nsite_res++;

          /* Set flag that we're done */
          done = TRUE;
        }

      /* Get the next site */
      site_ptr = site_ptr->next_site_ptr;
    }

  /* Return number of identified site residues */
  return(nsite_res);
}
/***********************************************************************

check_for_metal  -  Check whether the given atom is a metal ion

***********************************************************************/

int check_for_metal(char *atom_name,char *element)
{
  int imetal;
  int metal;

  /* Initialise variables */
  metal = FALSE;

  /* Loop through all the metals */
  for (imetal = 0; imetal < NMETALS && metal == FALSE; imetal++)
    {
      /* If name matches in first two characters, and third character
         is either a space or a number, then have a metal */
      if (!strncmp(atom_name,METAL_NAME[imetal],2) &&
          (atom_name[2] == ' ' ||
           (atom_name[2] >= '0' && atom_name[2] <= '9' &&
            atom_name[3] == ' ')))
        {
          /* Set flag */
          metal = TRUE;

          /* Return the metal name matched */
          strncpy(element,METAL_NAME[imetal],2);
          element[2] = '\0';
        }
    }

  /* Return the answer */
  return(metal);
}
/***********************************************************************

get_atom_type  -  Get atom's van der Waals radius

***********************************************************************/

int get_atom_type(char *atom_name,char stored_metal[NMETALS][3],
                  int *metal_count,int *nmetyp,int *nmet,char *element)
{
  int atom_type, itype, nmetal_types, nmetals;
  int found, metal;

  /* Initialise variables */
  atom_type = 0;
  element[0] = '\0';
  found = FALSE;
  nmetals = *nmet;
  nmetal_types = *nmetyp;

  /* Loop over the defined atom names */
  for (itype = 1; itype < NATOM_TYPES && found == FALSE; itype++)
    {
      /* Check for match - either exact, or for the non-standard
         cases as in FAD and NAD molecules where C can be AC, NC or
         PC, etc */
      if (!strncmp(atom_name,ELEMENT[itype],2) ||
          (atom_name[1] == ELEMENT[itype][1] &&
           (atom_name[0] == 'A' || atom_name[0] == 'N' ||
            atom_name[0] == 'P')))
        {
          /* Store type and element and set flag */
          atom_type = itype;
          strcpy(element,ELEMENT[itype]);
          found = TRUE;
        }
    }

  /* If no match found, check whether this might be a metal */
  if (found == FALSE)
    {
      /* Check whether this is a metal */
      metal = check_for_metal(atom_name,element);

      /* If metal, then see if we already have this one and increment
         counts */
      if (metal == TRUE)
        {
          /* Loop over the stored metal names */
          for (itype = 0; itype < nmetal_types && found == FALSE;
               itype++)
            {
              /* If have a match, update the count */
              if (!strcmp(element,stored_metal[itype]))
                {
                  /* Increment count */
                  metal_count[itype]++;

                  /* Set the atom type */
                  atom_type = NATOM_TYPES + itype;

                  /* Set flag */
                  found = TRUE;
                }
            }

          /* If don't already have this one, create a new slot for it */
          if (found == FALSE)
            {
              /* Store metal name and increment count */
              strcpy(stored_metal[nmetal_types],element);
                     
              /* Increment count */
              metal_count[nmetal_types]++;

              /* Set the atom type */
              atom_type = NATOM_TYPES + nmetal_types;

              /* Increment count of metal types */
              nmetal_types++;

              /* Set flag */
              found = TRUE;
            }

          /* Increment count of metals */
          nmetals++;
        }

      /* Otherwise, can't say what this atom is */
      else
        {
          /* Show warning message */
          printf("* Warning. Unrecognised atom type [%s]\n",
                 atom_name);

          /* Set atom type as unknown */
          atom_type = 0;
        }
    }

  /* Return current number of metal types */
  *nmet = nmetals;
  *nmetyp = nmetal_types;

  /* Return the atom type */
  return(atom_type);
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
                                char stored_metal[NMETALS][3],
                                int *metal_count,int *nmetyp,int *nmet,
                                int rec_type)
{
  int i, nmetal_types, nmetals;

  struct atom *atom_ptr;
  struct atom *first_atom_ptr, *last_atom_ptr;

  /* Initialise variables */
  nmetals = *nmet;
  nmetal_types = *nmetyp;

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
      printf("***        Program newsec terminated with error.\n");
      exit (1);
    }

  /* Store the details for this entry */
  atom_ptr->atom_number = atom_number;
  strcpy(atom_ptr->atom_name,atom_name);
  atom_ptr->original_coords[0] = x;
  atom_ptr->original_coords[1] = y;
  atom_ptr->original_coords[2] = z;
  for (i = 0; i < 3; i++)
    atom_ptr->coords[i]
      = atom_ptr->view_coords[i] = atom_ptr->original_coords[i];
  atom_ptr->bvalue = bvalue;
  atom_ptr->occupancy = occupancy;
  atom_ptr->residue_ptr = residue_ptr;
  atom_ptr->rec_type = rec_type;
  atom_ptr->alt = ' ';

  /* Increment count of this residue's atoms */
  residue_ptr->natoms++;

  /* If residue's first atom, store pointer to it */
  if (residue_ptr->first_atom_ptr == NULL)
    {
      /* Store atom pointer */
      residue_ptr->first_atom_ptr = atom_ptr;

      /* Store atom's coords as the max- and min-coords for the residue */
      for (i = 0; i < 3; i++)
        {
          residue_ptr->coord_min[i] = atom_ptr->coords[i];
          residue_ptr->coord_max[i] = atom_ptr->coords[i];
        }
    }

  /* Otherwise, update residue's min and max coordinates */
  else
    {
      for (i = 0; i < 3; i++)
        {
          if (atom_ptr->coords[i] < residue_ptr->coord_min[i])
            residue_ptr->coord_min[i] = atom_ptr->coords[i];
          if (atom_ptr->coords[i] > residue_ptr->coord_max[i])
            residue_ptr->coord_max[i] = atom_ptr->coords[i];
        }
    }

  /* If atom is one of the 3 mainchain atoms (CA, N or C), increment
     residue's count of these */
  if (!strcmp(atom_name," CA ") || !strcmp(atom_name," N  ") ||
      !strcmp(atom_name," C  "))
    residue_ptr->nmain_chain++;

  /* Check for C-alpha */
  if (!strcmp(atom_name," CA "))
    residue_ptr->have_calpha = TRUE;

  /* If this is a sidechain atom, store coords for calculation of
     side-chain's centroid */
  if (strcmp(atom_name," CA ") && strcmp(atom_name," N  ") &&
      strcmp(atom_name," C  ") && strcmp(atom_name," O  "))
    {
      /* Accumulate atom coords */
      for (i = 0; i < 3; i++)
        residue_ptr->centroid[i]
          = residue_ptr->centroid[i] + atom_ptr->coords[i];

      /* Increment count of centroid atoms */
      residue_ptr->ncentroid_atoms++;
    }

  /* Get the atom radius */
  atom_ptr->type
    = get_atom_type(atom_name,stored_metal,metal_count,&nmetal_types,
                    &nmetals,atom_ptr->element);

  /* Initialise pointers */
  atom_ptr->first_conect_ptr = NULL;
  atom_ptr->last_conect_ptr = NULL;
  atom_ptr->link_atom_ptr = NULL;

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

  /* Return updated count of metal types */
  *nmet = nmetals;
  *nmetyp = nmetal_types;

  /* Return current pointers */
  *fst_atom_ptr = first_atom_ptr;
  *lst_atom_ptr = last_atom_ptr;
  return(atom_ptr);
}
/***********************************************************************

create_rbond_record  -  Create a connection from the first residue to
                          the second

***********************************************************************/

void create_rbond_record(struct residue *residue1_ptr,
                         struct residue *residue2_ptr)
{
  struct rbond *rbond_ptr;

  /* Create rbond entry */
  rbond_ptr = (struct rbond*)malloc(sizeof(struct rbond));
  if (rbond_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for");
      printf(" struct rbond\n");
      printf("***        Program newsec terminated with error.\n");
      exit (1);
    }

  /* Store the details for this entry */
  rbond_ptr->to_residue_ptr = residue2_ptr;
  rbond_ptr->next_rbond_ptr = NULL;

  /* If this is the first entry stored, then update pointer to
     head of linked list */
  if (residue1_ptr->first_rbond_ptr == NULL)
    residue1_ptr->first_rbond_ptr = rbond_ptr;

  /* Otherwise, make current entry point to the first one and make it
     the new first */
  else
    {
      rbond_ptr->next_rbond_ptr = residue1_ptr->first_rbond_ptr;
      residue1_ptr->first_rbond_ptr = rbond_ptr;
    }                 
}
/***********************************************************************

create_conect_record  -  Create a connection from the first atom to
                         the second

***********************************************************************/

void create_conect_record(struct atom *atom1_ptr,struct atom *atom2_ptr,
                          int from_conect)
{
  struct conect *conect_ptr;

  /* Create conect entry */
  conect_ptr = (struct conect*)malloc(sizeof(struct conect));
  if (conect_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for");
      printf(" struct conect\n");
      printf("***        Program newsec terminated with error.\n");
      exit (1);
    }

  /* Store the details for this entry */
  conect_ptr->to_atom_ptr = atom2_ptr;
  conect_ptr->from_conect = from_conect;
  conect_ptr->next_conect_ptr = NULL;

  /* If this is the first entry stored, then update pointer to
     head of linked list */
  if (atom1_ptr->first_conect_ptr == NULL)
    atom1_ptr->first_conect_ptr = conect_ptr;

  /* Otherwise, make previous entry point to the current one */
  else
    atom1_ptr->last_conect_ptr->next_conect_ptr = conect_ptr;

  /* If atoms belong to different residues, create a connect record
     between the two residues */
  if (atom1_ptr->residue_ptr != atom2_ptr->residue_ptr)
    create_rbond_record(atom1_ptr->residue_ptr,atom2_ptr->residue_ptr);
                      
  /* Save the current conect's pointer */
  atom1_ptr->last_conect_ptr = conect_ptr;
}
/***********************************************************************

process_conect_record  -  Check whether the given CONECT record applies
                          to one of the ATOM or HETATM records stored
                          written out and hence is required

***********************************************************************/

void process_conect_record(char *input_line,struct atom **atom_index_ptr,
                           int natoms)
{
  char atmnum[10];

  int atom_number, iatom, iconect, ipos;
  int done;

  double cutoff_sqrd, dist_sqrd, x1, x2, y1, y2, z1, z2;

  struct atom *atom_ptr, *to_atom_ptr;

  /* Initialise variables */
  atom_ptr = NULL;
  to_atom_ptr = NULL;
  iatom = 0;
  done = FALSE;
  ipos = 6;

  /* Loop through all the possible atom-number positions */
  for (iconect = 0; iconect < MAXCONECT && done == FALSE; iconect++)
    {
      /* If this position is blank, then have reached end of atom-numbers
         on this line */
      if (!strncmp(input_line+ipos,"     ",5))
        done = TRUE;

      /* Otherwise, check for this atom number */
      else
        {
          /* Get the atom-number at this position */
          strncpy(atmnum,input_line+ipos,5);
          atmnum[5] = '\0';
          atom_number = atoi(atmnum);

          /* Check whether this atom number is valid */
          if (atom_number > 0 && atom_number <= natoms)
            {
              /* If the first atom number, store as first pointer */
              if (iatom == 0)
                {
                  /* Get the corresponding atom pointer */
                  atom_ptr = atom_index_ptr[atom_number - 1];

                  /* If this is null, then can't make use of the
                     current CONECT line */
                  if (atom_ptr == NULL)
                    done = TRUE;

                  /* Otherwise, retrieve atom's coordinates */
                  else
                    {
                      x1 = atom_ptr->coords[0];
                      y1 = atom_ptr->coords[1];
                      z1 = atom_ptr->coords[2];
                    }
                }

              /* Otherwise, get atom first is bonded to and check out
                 this bond */
              else
                {
                  /* If first position is null, then no point looking
                     at this, or any further, positions */
                  if (atom_ptr == NULL)
                    done = TRUE;

                  /* Otherwise, get the bonded atom */
                  else
                    {
                      /* Get the corresponding atom pointer */
                      to_atom_ptr = atom_index_ptr[atom_number - 1];

                      /* If this is not null, then check out the
                         connection */
                      if (to_atom_ptr != NULL)
                        {
                          /* Retrieve atom's coordinates */
                          x2 = to_atom_ptr->coords[0];
                          y2 = to_atom_ptr->coords[1];
                          z2 = to_atom_ptr->coords[2];

                          /* Calculate squared distance between
                             these two atoms */
                          dist_sqrd
                            = (x2 - x1) * (x2 - x1)
                            + (y2 - y1) * (y2 - y1)
                            + (z2 - z1) * (z2 - z1);

                          /* Check whether either atom is a metal */
                          if (atom_ptr->type > NATOM_TYPES - 1 ||
                              to_atom_ptr->type > NATOM_TYPES - 1)
                            cutoff_sqrd = MAX_METAL_LEN_SQRD;
                          else
                            cutoff_sqrd = MAX_CONECT_LEN_SQRD;

                          /* If distance is a reasonable one, store
                             a connection between these two atoms */
                          if (dist_sqrd < cutoff_sqrd)
                            {
                              /* Create a connection from the first
                                 atom to the second */
                              create_conect_record(atom_ptr,to_atom_ptr,
                                                   TRUE);

                              /* Create a connection the other way
                                 round */
                              create_conect_record(to_atom_ptr,atom_ptr,
                                                   TRUE);
                            }

                          /* If distance too great, show warning message */
                          else
                            {
                              printf("* Warning. Long bond from CONECT: "
                                     "%5d. %s %s %s %c -> %5d. %s %s %s %c"
                                     "%8.3f\n",
                                     atom_ptr->atom_number,
                                     atom_ptr->atom_name,
                                     atom_ptr->residue_ptr->res_name,
                                     atom_ptr->residue_ptr->res_num,
                                     atom_ptr->residue_ptr->chain,
                                     to_atom_ptr->atom_number,
                                     to_atom_ptr->atom_name,
                                     to_atom_ptr->residue_ptr->res_name,
                                     to_atom_ptr->residue_ptr->res_num,
                                     to_atom_ptr->residue_ptr->chain,
                                     sqrt(dist_sqrd));
                            }
                        }
                    }
                }

              /* Increment count of atom position */
              iatom++;
            }
        }

      /* Shift position along input line */
      ipos = ipos + 5;
    }
}
/***********************************************************************

create_site_record  -  Create and initialise a new site record

***********************************************************************/

struct site *create_site_record(struct site **fst_site_ptr,
                                struct site **lst_site_ptr,char *name,
                                int number,char *res_name,char *res_num,
                                char chain)
{
  struct site *site_ptr, *first_site_ptr, *last_site_ptr;

  /* Initialise site pointers */
  site_ptr = NULL;
  first_site_ptr = *fst_site_ptr;
  last_site_ptr = *lst_site_ptr;

  /* Allocate memory for structure to hold site info */
  site_ptr = (struct site *) malloc(sizeof(struct site));
  if (site_ptr == NULL)
    {
      printf("*** Can't allocate memory for struct site\n");
      exit (1);
    }

  /* If this is the very first site, then store pointer */
  if (first_site_ptr == NULL)
    first_site_ptr = site_ptr;

  /* Add link from previous site to the current one */
  if (last_site_ptr != NULL)
    last_site_ptr->next_site_ptr = site_ptr;
  last_site_ptr = site_ptr;

  /* Store the site details */
  store_string(&(site_ptr->name),name);
  site_ptr->number = number;
  strcpy(site_ptr->res_name,res_name);
  strcpy(site_ptr->res_num,res_num);
  site_ptr->chain = chain;
  site_ptr->next_site_ptr = NULL;

  /* Return current pointers */
  *fst_site_ptr = first_site_ptr;
  *lst_site_ptr = last_site_ptr;

  /* Return new site pointer */
  return(site_ptr);
}
/***********************************************************************

create_hetdesc_record  -  Create a new hetdesc record or update
                           existing one if we already have it

***********************************************************************/

struct hetdesc *create_hetdesc_record(struct hetdesc **fst_hetdesc_ptr,
                                      struct hetdesc **lst_hetdesc_ptr,
                                      char *res_name,int *new)
{
  int found;

  struct hetdesc *hetdesc_ptr, *match_hetdesc_ptr;
  struct hetdesc *first_hetdesc_ptr, *last_hetdesc_ptr;

  /* Initialise variables */
  found = FALSE;
  *new = FALSE;

  /* Initialise hetdesc pointers */
  last_hetdesc_ptr = *lst_hetdesc_ptr;
  first_hetdesc_ptr = *fst_hetdesc_ptr;
  match_hetdesc_ptr = NULL;

  /* Get pointer to first entry */
  hetdesc_ptr = first_hetdesc_ptr;

  /* Loop through all these to determine whether we already have this
     entry */
  while (hetdesc_ptr != NULL && found == FALSE)
    {
      /* If have match on name, then already have this entry */
      if (!strcmp(hetdesc_ptr->res_name,res_name))
        {
          /* Store pointer to located entry */
          match_hetdesc_ptr = hetdesc_ptr;

          /* Increment count of matches */
          hetdesc_ptr->count++;

          /* Set flag that found */
          found = TRUE;
        }

      /* Get next hetdesc entry */
      hetdesc_ptr = hetdesc_ptr->next_hetdesc_ptr;
    }

  /* If match found, then retrieve pointer */
  if (found == TRUE)
    hetdesc_ptr = match_hetdesc_ptr;

  /* Otherwise, create a new hetdesc entry */
  else
    {
      hetdesc_ptr = (struct hetdesc*)malloc(sizeof(struct hetdesc));
      if (hetdesc_ptr == NULL)
        {
          printf("*** ERROR. Unable to allocate memory for");
          printf(" struct hetdesc\n");
          printf("***        Program newsec terminated with error.\n");
          exit (1);
        }

      /* Set flag to indicate that this is a new entry */
      *new = TRUE;

      /* Store the details for this entry */
      strcpy(hetdesc_ptr->res_name,res_name);
      hetdesc_ptr->details[DESCRIPTION] = &null;
      hetdesc_ptr->details[SYNONYM] = &null;
      hetdesc_ptr->details[FORMULA] = &null;
      hetdesc_ptr->count = 0;

      /* Initialise pointer to next hetdesc record */
      hetdesc_ptr->next_hetdesc_ptr = NULL;

      /* If this is the first entry stored, then update pointer to
         head of linked list */
      if (first_hetdesc_ptr == NULL)
        first_hetdesc_ptr = hetdesc_ptr;

      /* Otherwise, make previous entry point to the current
         one */
      else
        last_hetdesc_ptr->next_hetdesc_ptr = hetdesc_ptr;
                      
      /* Save the current hetdesc's pointer */
      last_hetdesc_ptr = hetdesc_ptr;
    }

  /* Return current pointers */
  *fst_hetdesc_ptr = first_hetdesc_ptr;
  *lst_hetdesc_ptr = last_hetdesc_ptr;
  return(hetdesc_ptr);
}
/***********************************************************************

get_qstruc  -  Determine the quaternary strcutrue from the description

***********************************************************************/

int get_qstruc(char *unit_type,char *prefix,char *qstruc_desc)
{
  char ch, number_string[20];

  char *string_ptr;

  int i, itype, len, qstruc;
  int done, found;

  /* Initialise variables */
  found = FALSE;
  qstruc = NOT_GIVEN;

  /* Remove common prefixes */
  if (!strncmp(unit_type,"HOMO",4))
    {
      strcpy(prefix,"HOMO-");
      if (unit_type[4] == '-')
        strcpy(unit_type,unit_type+5);
      else
        strcpy(unit_type,unit_type+4);
    }
  else if (!strncmp(unit_type,"HETERO",6))
    {
      strcpy(prefix,"HETERO-");
      if (unit_type[6] == '-')
        strcpy(unit_type,unit_type+7);
      else
        strcpy(unit_type,unit_type+6);
    }

  /* Remove common suffixes */
  string_ptr = strstr(unit_type,"MERIC");
  if (string_ptr != NULL)
    string_ptr[3] = '\0';

  /* Remove semi-colon, where present */
  string_ptr = strchr(unit_type,';');
  if (string_ptr != NULL)
    string_ptr[0] = '\0';

  /* Remove space, where present */
  string_ptr = strchr(unit_type,' ');
  if (string_ptr != NULL)
    string_ptr[0] = '\0';

  /* Prepare to identify unit type */
  done = FALSE;
  number_string[0] = '\0';
  len = strlen(unit_type);

  /* Check for number at start */
  for (i = 0; i < len && done == FALSE; i++)
    {
      /* Get this character */
      ch = unit_type[i];

      /* Check whether this character is numeric */
      if (ch >= '0' && ch <= '9')
        {
          number_string[i] = ch;
          number_string[i + 1] = '\0';
        }
      else
        done = TRUE;
    }

  /* If have a number, store equivalent quaternary structure */
  if (number_string[0] != '\0')
    {
      /* Get number */
      qstruc = atoi(number_string);

      /* Set flag that found */
      found = TRUE;
    }

  /* Identify the unit type */
  for (itype = 1; itype < NUNIT_TYPES + 1 && found == FALSE; itype++)
    {
      /* If matches, then store number and set flag */
      if (!strcmp(unit_type,UNIT_TYPE[itype]))
        {
          qstruc = itype;
          found = TRUE;
        }
    }

  /* Check the extra unit type names */
  for (itype = 0; itype < NEXTRA_TYPES && found == FALSE; itype++)
    {
      /* If matches, then store number */
      if (!strcmp(unit_type,EXTRA_UNIT_TYPE[itype]))
        {
          qstruc = EXTRA_NUMBER[itype];
          found = TRUE;
        }
    }

  /* If not found, show warning */
  if (found == FALSE)
    printf("* Warning. Unidentified quaternary structure type: %s\n",
           unit_type);

  /* Otherwise, form name of quaternary structure */
  else
    {
      /* Get type */
      itype = qstruc;

      /* If from 1 to 12, give full name */
      if (itype > 0 && itype < 13)
        strcpy(qstruc_desc,UNIT_TYPE[itype]);

      /* Otherwise, use number version */
      else
        sprintf(qstruc_desc,"%dMER",itype);
    }

  /* Return the quaternary structure */
  return(qstruc);
}
/***********************************************************************

read_pdb  -  Read in the given PDB file

***********************************************************************/

int read_pdb(char *pdb_name,struct chain **fst_chain_ptr,
             struct chain **lst_chain_ptr,int *nchn,
             struct atom **fst_atom_ptr,
             struct residue **fst_residue_ptr,int *nres,int *nsres,
             char stored_metal[NMETALS][3],int *metal_count,
             int *nmetyp,int *nmet,struct data *data_item,
             struct hetdesc **fst_hetdesc_ptr,int *mcsg,
             struct biomol *biomol_ptr,struct parameters params)
{
  char alt, ch, input_line[LINELEN + 1], number_string[20];
  char tmp_line[LINELEN + 1];
  char description[MAX_DESC_LEN + 1], tmp_string[100];
  char atom_name[5], chain, element[3], res_name[4], res_num[6], occup[4];
  char last_chain, last_res_name[4], last_res_num[6], site_name[4];
  char last_het[4], last_num[6], unit_type[LINELEN + 1];
  char chain_list[LINELEN + 1];

  char *string_ptr;

  int atom_number, natoms, nlines, nread, nresid, top_atom_no;
  int count, htype, i, iatom, ipos, ires, itype, j, k, len, nchains;
  int nmatrix, nmetal_types, nmetals, nmodels, nsite_res, ntransform;
  int nidentity, number;
  int chain_end, have_atom, have_chain, keep_reading;
  int rec_type;
  int unit;
  int done, ensemble, eof, gzipped, first_conect, found, hydrogen, new;
  int new_transform, nmr, refresh, wanted, warned;

  float matrix[3][3], translation[3];

  double bvalue, occupancy, inx, iny, inz, x, y, z;

  struct atom *atom_ptr, *first_atom_ptr, *last_atom_ptr;
  struct chain *chain_ptr, *first_chain_ptr, *last_chain_ptr;
  struct hetdesc *hetdesc_ptr, *first_hetdesc_ptr, *last_hetdesc_ptr;
  struct residue *residue_ptr, *first_residue_ptr, *last_residue_ptr;
  struct site *site_ptr, *first_site_ptr, *last_site_ptr;

  struct atom **atom_index_ptr;

#ifdef TEST
  FILE *file_ptr;
#else
  FILE *file_ptr;
  gzFile gzfile_ptr;
#endif
  /* Initialise counts */
  printf("Reading PDB file (%s) ...\n",pdb_name);
  chain_end = FALSE;
  count = 0;
  ensemble = FALSE;
  first_conect = TRUE;
  have_chain = FALSE;
  keep_reading = TRUE;
  last_res_name[0] = '\0';
  last_res_num[0] = '\0';
  last_het[0] = '\0';
  last_num[0] = '\0';
  last_chain = '\0';
  *mcsg = FALSE;
  natoms = 0;
  nchains = 0;
  new_transform = TRUE;
  nidentity = 0;
  nmetal_types = 0;
  nread = 0;
  nlines = 0;
  nmetals = 0;
  nmodels = 0;
  nmr = FALSE;
  nresid = 0;
  nsite_res = 0;
  number = 0;
  top_atom_no = 0;
  for (i = 0; i < NMETALS; i++)
    {
      if (i < NCOMMON_METALS)
        strcpy(stored_metal[i],COMMON_METAL[i]);
      else
        stored_metal[i][0] = '\0';
      metal_count[i] = 0;
    }
  nmetal_types = NCOMMON_METALS;
  warned = FALSE;

  /* Initialise pointers */
  atom_ptr = NULL;
  first_atom_ptr = NULL;
  last_atom_ptr = NULL;
  chain_ptr = NULL;
  first_chain_ptr = NULL;
  last_chain_ptr = NULL;
  hetdesc_ptr = NULL;
  first_hetdesc_ptr = NULL;
  last_hetdesc_ptr = NULL;
  residue_ptr = NULL;
  first_residue_ptr = NULL;
  last_residue_ptr = NULL;
  site_ptr = NULL;
  first_site_ptr = NULL;
  last_site_ptr = NULL;

  /* Initialise all the data items */
  strcpy(data_item->pdb_code,params.pdb_code);
  data_item->resolution = 0.0;
  data_item->rfactor = 0.0;
  data_item->rfree = 0.0;
  data_item->nmr = FALSE;
  data_item->theoretical_model = FALSE;
  data_item->ensemble = FALSE;
  data_item->nmodels = 0;
  data_item->nwaters = 0;
  data_item->nhydrogens = 0;
  data_item->nchains = 0;
  data_item->nprotein = 0;
  data_item->nca_only = 0;
  data_item->ndna = 0;
  data_item->nligands = 0;
  data_item->nmetals = 0;
  data_item->nattach = 0;
  data_item->nunknown = 0;
  data_item->nresidues = 0;
  data_item->natom_records = 0;
  data_item->nhetatm_records = 0;

  /* Initialise biomolecule data */
  biomol_ptr->ntransform = 0;
  biomol_ptr->nsymmetry = 0;
  biomol_ptr->nmolecules = 0;
  biomol_ptr->nbiomolecule = 0;
  biomol_ptr->is_biomolecule = FALSE;
  biomol_ptr->nmatrix = 0;
  biomol_ptr->nidentity = 0;
  biomol_ptr->new_style = FALSE;
  biomol_ptr->split_asu = FALSE;
  biomol_ptr->skip = FALSE;
  biomol_ptr->prefix[0] = '\0';
  biomol_ptr->quaternary_structure = NOT_GIVEN;
  biomol_ptr->qstruc_desc[0] = '\0';
  biomol_ptr->nsymop = 0;
  biomol_ptr->in_pqs = FALSE;
  biomol_ptr->use_pqs = FALSE;
  biomol_ptr->pqs_split_asu = FALSE;
  biomol_ptr->pqs_skip = FALSE;
  biomol_ptr->pqs_prefix[0] = '\0';
  biomol_ptr->pqs_quaternary_structure = NOT_GIVEN;
  biomol_ptr->pqs_qstruc_desc[0] = '\0';
  biomol_ptr->have_author_info = FALSE;
  biomol_ptr->author_skip = FALSE;
  biomol_ptr->author_prefix[0] = '\0';
  biomol_ptr->author_quaternary_structure = NOT_GIVEN;
  biomol_ptr->author_qstruc_desc[0] = '\0';
  for (i = 0; i < MAX_BIOMT; i++)
    {
      biomol_ptr->first_matrix[i] = 0;
      biomol_ptr->type[i] = NO_MATRIX;
      biomol_ptr->tchains[i] = &null;
      biomol_ptr->prefix[0] = '\0';
      for (j = 0; j < 3; j++)
        for (k = 0; k < 3; k++)
          {
            if (j == k)
              biomol_ptr->biomol_transform[i][j][k] = 1.0;
            else
              biomol_ptr->biomol_transform[i][j][k] = 0.0;
          }
    }
  biomol_ptr->have_biomol_image = FALSE;

  /* Initialise matrices */
  for (j = 0; j < 3; j++)
    {
      for (k = 0; k < 3; k++)
        {
          if (j == k)
            matrix[j][k] = 1.0;
          else
            matrix[j][k] = 0.0;
        }
      translation[j] = 0.0;
    }
  chain_list[0] = '\0';

  /* If first digit of the PDB code is non-numeric, take it to be
     an MCSG structure */
  if (params.pdb_code[0] < '0' || params.pdb_code[0] > '9')
    *mcsg = TRUE;

  /* Determine whether PDB file is a gzipped one */
  eof = gzipped = FALSE;
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
  if (eof == TRUE)
     {
      printf("echo \"*** ERROR. Unable to open PDB file ");
      printf("[%s] \"\n",pdb_name);
      printf("echo \" \"\n");
      exit(-1);
    }

  /* Loop through the PDB file, picking up all the atoms that are
     required */
  while (eof == FALSE)
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
          /* Truncate the string to strip off the newline */
          len = string_truncate(input_line,LINELEN);

          /* Increment count of lines read in */
          nlines++;

          /* Check for ATOM or HETATM record */
          if (keep_reading == TRUE &&
              (!strncmp(input_line,"ATOM",4) ||
               !strncmp(input_line,"HETATM",6) ||
               !strncmp(input_line,"ATZERO",6) ||
               !strncmp(input_line,"HEZERO",6)))
            {
              /* Increment count of coord records read in */
              nread++;

              /* Increment count record types read in */
              if (!strncmp(input_line,"ATOM",4) ||
                  !strncmp(input_line,"ATZERO",6))
                {
                  rec_type = ATOM_REC;
                  data_item->natom_records++;
                }
              if (!strncmp(input_line,"HETATM",6) ||
                  !strncmp(input_line,"HEZERO",6))
                {
                  rec_type = HETATM_REC;
                  data_item->nhetatm_records++;
                }

              /* Get all the atom details */
              get_atom_details(input_line,&atom_number,atom_name,
                               &chain,res_name,res_num,&inx,&iny,&inz,
                               &bvalue,&occupancy,&alt,element);

              /* Rotate through 180 degrees about the x-axis */
              x = inx;
              y = - iny;
              z = - inz;

              /* Determine whether this is a wanted ATOM or HETATM record */
              wanted = TRUE;

              /* NMR structures in the pdb have hydrogens and pseudo-atoms, 
                 so we need to remove these */
              if (!strcmp(element,"H"))
                hydrogen = TRUE;
              else
                hydrogen = is_hydrogen_atom(atom_name);
              if (hydrogen == TRUE)
                {
                  data_item->nhydrogens++;
                  wanted = FALSE;
                }

              /* If this is a water, then don't want it */
              if (!strcmp(res_name,"HOH") || !strcmp(res_name,"WAT"))
                {
                  data_item->nwaters++;
                  wanted = FALSE;
                }

              /* If this is an alternate-occupancy atom, check whether
                 we already have it */
              strncpy(occup,input_line+56,3);
              occup[3] = '\0';
              if (wanted == TRUE &&
                  (alt != ' ' || strncmp(occup,"1.0",3)))
                {
                  /* Check for earlier copy of this atom */
                  have_atom
                    = check_for_duplicate(residue_ptr,atom_name,
                                          res_name,res_num,chain);
          
                  /* If have this atom already, then don't want this
                     second copy */
                  if (have_atom == TRUE)
                    wanted = FALSE;
                }

              /* If atom still wanted, then prepare to store */
              if (wanted == TRUE)
                {
                  /* Check whether atom belongs to a new chain */
                  if (chain != last_chain)
                    {
                      /* Create a chain record */
                      chain_ptr
                        = create_chain_record(&first_chain_ptr,
                                              &last_chain_ptr,chain);

                      /* Increment count of chains */
                      nchains++;
                    }

                  /* Check whether atom belongs to a new residue */
                  if (chain != last_chain ||
                      strncmp(res_num,last_res_num,5) ||
                      strncmp(res_name,last_res_name,3))
                    {
                      /* Create a residue record */
                      residue_ptr
                        = create_residue_record(&first_residue_ptr,
                                                &last_residue_ptr,
                                                chain_ptr,res_name,
                                                res_num,chain,nresid);

                      /* Increment count of residues read in */
                      nresid++;

                      /* Determine whether this residue belongs to an
                         active site (as defined by the previously
                         read in SITE records) */
                      scan_sites(residue_ptr,first_site_ptr);

                      /* Save the current residue details */
                      strcpy(last_res_name,res_name);
                      strcpy(last_res_num,res_num);
                      last_chain = chain;
                    }

                  /* If have an object to which this atom's residue has
                     been assigned, then can store the atom */
                  if (residue_ptr != NULL)
                    {
                      /* Create atom record */
                      atom_ptr
                        = create_atom_record(&first_atom_ptr,&last_atom_ptr,
                                             residue_ptr,atom_number,
                                             atom_name,x,y,z,bvalue,occupancy,
                                             stored_metal,metal_count,
                                             &nmetal_types,&nmetals,rec_type);

                      /* Increment count of atoms retained */
                      if (atom_ptr != NULL)
                        natoms++;

                      /* If this is the highest atom number so far, store */
                      if (atom_number > top_atom_no)
                        top_atom_no = atom_number;

                      /* Update counts of this chain's ATOM/HETATM records */
                      if (rec_type == ATOM_REC)
                        atom_ptr->residue_ptr->chain_ptr->natom_recs++;
                      else
                        atom_ptr->residue_ptr->chain_ptr->nhetatm_recs++;

                      /* Store the alternate occupancy identifier */
                      atom_ptr->alt = alt;
                    }
                }
            }

          /* If this is a CONECT record, then store if it is wanted */
          else if (!strncmp(input_line,"CONECT",6))
            {
              /* If this is the first of the CONECT records set up array
                 linking atom numbers to atom pointers */
              if (first_conect == TRUE && top_atom_no > 0)
                {
                  /* Create an atom pointers array */
                  atom_index_ptr
                    = (struct atom **)
                    malloc(top_atom_no * sizeof(struct atom *));

                  /* Loop through the array to initialise */
                  for (iatom = 0; iatom < top_atom_no; iatom++)
                    atom_index_ptr[iatom] = NULL;

                  /* Get pointer to first atom record */
                  atom_ptr = first_atom_ptr;

                  /* Loop through all the atoms to update the pointers
                     array */
                  while (atom_ptr != NULL)
                    {
                      /* Get the atom number */
                      iatom = atom_ptr->atom_number;

                      /* If valid, stoe the atom's pointer in the array */
                      if (iatom > 0 && iatom <= top_atom_no)
                        {
                          /* Store atom pointer */
                          atom_index_ptr[iatom - 1] = atom_ptr;
                        }

                      /* Get the next atom record */
                      atom_ptr = atom_ptr->next_atom_ptr;
                    }

                  /* Unset keep-reading flag to ensure no more atoms are
                     read in */
                  keep_reading = FALSE;

                  /* Unset the flag */
                  first_conect = FALSE;
                }

              /* Process the current record */
              process_conect_record(input_line,atom_index_ptr,top_atom_no);
            }

          /* Check for the end of an NMR structure model */
          else if (!strncmp("ENDMDL",input_line,6) && natoms > 0)
            keep_reading = FALSE;

          /* If this is a model record, increment count of models */
          else if (!strncmp("MODEL ",input_line,6))
            {
              /* Set flag */
              ensemble = TRUE;

              /* Increment count of models */
              nmodels++;

              /* Check for notorious cases where first model actually
                 contains no atoms(!) */
              if (nmodels > 1 && natoms == 0)
                nmodels--;
            }

          /* If this is a TER record, have an end to the current chain */
          else if (!strncmp("TER",input_line,3))
            {
              /* Set last chain to NULL */
              last_chain = '\0';
            }

          /* Check the experimental method */
          else if (!strncmp("EXPDTA",input_line,6))
            {
              /* Check whether this structure solved by NMR */
              string_ptr = strstr(input_line,"NMR");
              if (string_ptr != NULL)
                nmr = TRUE;

              /* Check for theoretical model */
              string_ptr = strstr(input_line,"THEORETICAL MODEL");
              if (string_ptr != NULL)
                data_item->theoretical_model = TRUE;
            }

          /* Pick up the resolution */
          else if (!strncmp("REMARK   2 RESOLUTION",input_line,21))
            {
              /* Check for case where resolution not given */
              if (!strncmp(input_line+23,"NOT",3))
                data_item->resolution = -9.9;

              /* Otherwise, extract value */
              else
                {
                  /* Find the first space after keyword */
                  string_ptr = strchr(input_line+20,' ');

                  /* If found, then extract value */
                  if (string_ptr != NULL)
                    {
                      /* Copy to temporary line */
                      strcpy(tmp_line,string_ptr);

                      /* Remove ANGSTROM keyword */
                      string_ptr = strstr(tmp_line," ANGSTROM");

                      /* If valid, extract number and store */
                      if (string_ptr != NULL)
                        {
                          /* Terminate number string and convert to float */
                          string_ptr[0] = '\0';
                          data_item->resolution = atof(tmp_line);
                        }
                    }
                }
            }

          /* Pick up the R-factor(s) */
          else if (!strncmp("REMARK   3",input_line,10))
            {
              /* See if R-value given in this line */
              string_ptr = strstr(input_line,"R VALUE");

              /* If present, and we don't already have it, check for
                 different various formats */
              if (string_ptr != NULL && data_item->rfactor <= 0.0)
                {
                  /* Get start position */
                  ipos = string_ptr - input_line;

                  /* Check X-PLOR format */
                  if (ipos == 13)
                    {
                      /* Check for new format */
                      string_ptr = strstr(input_line,":");
                      if (string_ptr != NULL)
                        ipos = string_ptr - input_line;
                      else
                        ipos = 0;
                      if (ipos > 19)
                        ipos = ipos + 2;
                      else
                        ipos = 40;

                      /* Read in the R-factor value */
                      strncpy(number_string,input_line+ipos,5);
                      number_string[5] = '\0';
                      data_item->rfactor = atof(number_string);
                    }

                  /* Check text version */
                  else if (!strncmp(input_line+ipos+7," IS ",4) &&
                           ipos + 15 < 72)
                    {
                      /* Read in the R-factor value */
                      strncpy(number_string,input_line+ipos+11,5);
                      number_string[5] = '\0';
                      data_item->rfactor = atof(number_string);
                    }
                }

              /* See if R-free given in this line */
              string_ptr = strstr(input_line,"FREE R VALUE   ");

              /* If present, and we don't already have it, check for
                 different various formats */
              if (string_ptr != NULL && data_item->rfree <= 0.0)
                {
                  /* Get start position */
                  ipos = string_ptr - input_line;

                  /* Check X-PLOR format */
                  if (ipos == 13)
                    {
                      /* Check for new format */
                      string_ptr = strstr(input_line,":");
                      if (string_ptr != NULL)
                        ipos = string_ptr - input_line;
                      else
                        ipos = 0;
                      if (ipos > 19)
                        ipos = ipos + 2;
                      else
                        ipos = 40;

                      /* Read in the R-free value */
                      strncpy(number_string,input_line+ipos,5);
                      number_string[5] = '\0';
                      data_item->rfree = atof(number_string);
                    }
                }
            }

          /* If this is a SITE record, store the site details */
          else if (!strncmp("SITE  ",input_line,6))
            {
              /* Add blank to end of line to ensure that truncated lines
                 are correctly interpreted */
              strcat(input_line,"                                      ");

              /* Get the site name */
              strncpy(site_name,input_line+11,3);
              site_name[3] = '\0';

              /* Initialise the record position for reading in the residues
                 in this record */
              ipos = 18;
              done = FALSE;

              /* Process each of the potential site residues in turn */
              for (ires = 0; ires < 4 && done == FALSE; ires++)
                {
                  /* Get the chain-ID, residue-name and number */
                  strncpy(res_name,input_line+ipos,3);
                  res_name[3] = '\0';
                  strncpy(res_num,input_line+ipos+5,5);
                  res_num[5] = '\0';
                  chain = input_line[ipos+4];

                  /* If all is blank, then we are done */
                  if (!strcmp(res_name,"   ") || !strcmp(res_num,"     "))
                    done = TRUE;

                  /* Otherwise, create a site record to store these details */
                  else
                    {
                      /* Increment count of site residues */
                      nsite_res++;

                      /* Create the record */
                      site_ptr
                        = create_site_record(&first_site_ptr,
                                             &last_site_ptr,site_name,
                                             nsite_res,res_name,res_num,
                                             chain);
                    }

                  /* Increment position for next residue */
                  ipos = ipos + 11;
                }
            }

          /* If this is a record describing a Het Group get the details */
          else if (!strncmp("HET   ",input_line,6) ||
                   !strncmp("HETNAM",input_line,6) ||
                   !strncmp("HETSYN",input_line,6) ||
                   !strncmp("FORMUL",input_line,6))
            {
              /* Get the residue name from\ the appropriate position in the
                 line */
              tmp_string[0] = '\0';
              refresh = FALSE;
              if (!strncmp("HET   ",input_line,6) && len > 30)
                {
                  strncpy(res_name,input_line+7,3);
                  res_name[3] = '\0';
                  strncpy(res_num,input_line+13,5);
                  res_num[5] = '\0';
                  htype = DESCRIPTION;
                  strcpy(tmp_string,input_line+30);
                  tmp_string[40] = '\0';

                  /* If this is the same Het Group as we've just read, but
                     for a different residue, then don't want this
                     description */
                  if (!strcmp(res_name,last_het) && strcmp(res_num,last_num))
                    tmp_string[0] = '\0';

                  /* Store current Het Group name and residue number */
                  strcpy(last_het,res_name);
                  strcpy(last_num,res_num);
                }
              else if (!strncmp("HETNAM",input_line,6) && len > 15)
                {
                  strncpy(res_name,input_line+11,3);
                  res_name[3] = '\0';
                  htype = DESCRIPTION;
                  strcpy(tmp_string,input_line+15);
                  tmp_string[55] = '\0';

                  /* If this is the very first HETNAM record for this
                     residue type, then overwrite any description that
                     may have come from the HET record */
                  if (input_line[9] == ' ')
                    refresh = TRUE;
                }
              else if (!strncmp("HETSYN",input_line,6) && len > 15)
                {
                  strncpy(res_name,input_line+11,3);
                  res_name[3] = '\0';
                  htype = SYNONYM;
                  strcpy(tmp_string,input_line+15);
                  tmp_string[55] = '\0';
                }
              else if (!strncmp("FORMUL",input_line,6) && len > 19)
                {
                  strncpy(res_name,input_line+12,3);
                  res_name[3] = '\0';
                  htype = FORMULA;
                  strcpy(tmp_string,input_line+19);
                  tmp_string[51] = '\0';
                }

              /* Truncate the description */
              string_truncate(tmp_string,60);

              /* If have a description, store it */
              if (tmp_string[0] != '\0')
                {
                  /* Create a hetdesc record, if we don't already have it */
                  hetdesc_ptr
                    = create_hetdesc_record(&first_hetdesc_ptr,
                                            &last_hetdesc_ptr,res_name,&new);

                  /* Retrieve the appropriate field */
                  if (refresh == FALSE)
                    strcpy(description,hetdesc_ptr->details[htype]);
                  else
                    description[0] = '\0';

                  /* If have an existing description, then free memory
                     associated with it */
                  if (hetdesc_ptr->details[htype][0] != '\0')
                    free(hetdesc_ptr->details[htype]);

                  /* Append current description */
                  strcat(description," ");
                  strcat(description,tmp_string);

                  /* Store the new string */
                  store_string(&(hetdesc_ptr->details[htype]),description);
                }
            }

          /* Check for quaternary structure */
          else if (!strncmp(input_line,"REMARK 300 QUATERNARY STRUCTURE",31) ||
                   (!strncmp(input_line,"COMPND",6) &&
                    !strncmp(input_line+11,"BIOLOGICAL_UNIT:",16)) ||
                   !strncmp(input_line,"REMARK 300 BIOLOGICAL UNIT:",27))
            {
              /* Get the location of the colon */
              string_ptr = strchr(input_line,':');

              /* If present, get description of quaternary structure */
              if (string_ptr != NULL)
                {
                  /* Extract description */
                  strcpy(unit_type,string_ptr + 2);

                  /* Extract quaternary structure type and any prefixes */
                  biomol_ptr->quaternary_structure
                    = get_qstruc(unit_type,biomol_ptr->prefix,
                                 biomol_ptr->qstruc_desc);
                }
            }

          /* Check for start of biomolecule transformation */
          else if (!strncmp(input_line,
                            "REMARK 350 APPLY THE FOLLOWING TO CHAIN",39))
            {
              /* Initialise matrices */
              for (j = 0; j < 3; j++)
                {
                  for (k = 0; k < 3; k++)
                    {
                      if (j == k)
                        matrix[j][k] = 1.0;
                      else
                        matrix[j][k] = 0.0;
                    }
                  translation[j] = 0.0;
                }
              chain_list[0] = '\0';

              /* Get the location of the colon */
              string_ptr = strchr(input_line,':');

              /* If present, get list of chains to which following
                 transformation is to apply */
              if (string_ptr != NULL)
                {
                  /* Extract chain list */
                  strcpy(chain_list,string_ptr + 2);
                }

              /* Set flag to indicate that we have the start of a new
                 set of transformations */
              new_transform = TRUE;
            }

          /* Check for start of biomolecule transformation */
          else if (!strncmp(input_line,"REMARK 350   BIOMT",18) ||
                   !strncmp(input_line,"REMARK 290   SMTRY",18))
            {
              /* Get the matrix row */
              strncpy(number_string,input_line+18,1);
              number_string[1] = '\0';
              j = atoi(number_string) - 1;

              /* If valid, extract the transformation values */
              if (j > -1 && j < 3)
                {
                  /* Extract the current row for the rotation matrix */
                  strncpy(number_string,input_line+24,9);
                  number_string[9] = '\0';
                  matrix[j][0] = atof(number_string);
                  strncpy(number_string,input_line+34,9);
                  number_string[9] = '\0';
                  matrix[j][1] = atof(number_string);
                  strncpy(number_string,input_line+44,9);
                  number_string[9] = '\0';
                  matrix[j][2] = atof(number_string);

                  /* Extract the translation */
                  strncpy(number_string,input_line+54,14);
                  number_string[14] = '\0';
                  translation[j] = atof(number_string);

                  /* If this is the last row, can store this transformation */
                  if (j == 2)
                    {
                      /* Check whether we simply have the unit matrix */
                      unit = TRUE;

                      /* Loop through matrix elements to check */
                      for (j = 0; j < 3; j++)
                        {
                          for (k = 0; k < 3; k++)
                            {
                              if ((j == k && matrix[j][k] != 1.0) ||
                                  (j != k && matrix[j][k] != 0.0))
                                unit = FALSE;
                            }
                          if (translation[j] != 0.0)
                            unit = FALSE;
                        }

                      /* Set flag that we have new-style transformation
                         data */
                      biomol_ptr->new_style = TRUE;

                      /* If transformation is not the identity transform,
                         store it */
                      if (unit == FALSE)
                        {
                          /* Get which biomolecule this is */
                          nmatrix = biomol_ptr->nmatrix;

                          /* If within allowed range, store matrix */
                          if (nmatrix < MAX_BIOMT)
                            {
                              /* If this is a new transformation, store its
                                 starting position */
                              if (new_transform == TRUE)
                                {
                                  /* Store array position of starting
                                     transformation */
                                  ntransform = biomol_ptr->ntransform;
                                  biomol_ptr->first_matrix[ntransform]
                                    = nmatrix;

                                  /* Store list of chains corresponding to
                                     this set of transformations */
                                  store_string(&(biomol_ptr->tchains[ntransform]),
                                               chain_list);

                                  /* Increment transform number */
                                  biomol_ptr->ntransform++;
                                }

                              /* Store matrix type */
                              if (!strncmp(input_line+13,"SMTRY",5))
                                {
                                  biomol_ptr->type[nmatrix] = SYMMETRY;
                                  biomol_ptr->nsymmetry++;
                                }
                              else
                                {
                                  biomol_ptr->type[nmatrix] = BIOMOLECULE;
                                  biomol_ptr->nbiomolecule++;
                                }

                              /* Loop through matrix elements to store */
                              for (j = 0; j < 3; j++)
                                for (k = 0; k < 3; k++)
                                  biomol_ptr->biomol_transform[nmatrix][j][k]
                                    = matrix[j][k];

                              /* Reset flag */
                              new_transform = FALSE;

                              /* Increment biomolecule number */
                              biomol_ptr->nmatrix++;
                            }

                          /* Otherwise, show warning */
                          else if (warned == FALSE)
                            {
                              printf("* Warning. Maximum number of biomol "
                                     "transformations exceeded. Some data "
                                     "lost\n");
                              warned = TRUE;
                            }
                        }

                      /* If have identity transformation, increment count */
                      else
                        biomol_ptr->nidentity++;
                    }
                }

              /* If not valid, show warning */
              else
                printf("* Warning. Invalid matrix row: %s\n",input_line);
            }

          /* Check for old-style transformation */
          else if (!strncmp(input_line,"REMARK   5 MTRIX",16))
            {
              /* Get the matrix row */
              strncpy(number_string,input_line+16,1);
              number_string[1] = '\0';
              j = atoi(number_string) - 1;

              /* If valid, extract the transformation values */
              if (j > -1 && j < 3)
                {
                  /* Extract the current row for the rotation matrix */
                  strncpy(number_string,input_line+22,9);
                  number_string[9] = '\0';
                  matrix[j][0] = atof(number_string);
                  strncpy(number_string,input_line+32,9);
                  number_string[9] = '\0';
                  matrix[j][1] = atof(number_string);
                  strncpy(number_string,input_line+42,9);
                  number_string[9] = '\0';
                  matrix[j][2] = atof(number_string);

                  /* Extract the translation */
                  strncpy(number_string,input_line+52,14);
                  number_string[14] = '\0';
                  translation[j] = atof(number_string);

                  /* If this is the last row, can store this transformation */
                  if (j == 2)
                    {
                      /* Check whether we simply have the unit matrix */
                      unit = TRUE;

                      /* Loop through matrix elements to check */
                      for (j = 0; j < 3; j++)
                        {
                          for (k = 0; k < 3; k++)
                            {
                              if ((j == k && matrix[j][k] != 1.0) ||
                                  (j != k && matrix[j][k] != 0.0))
                                unit = FALSE;
                            }
                          if (translation[j] != 0.0)
                            unit = FALSE;
                        }

                      /* If transformation is not the identity transform,
                         store it */
                      if (unit == FALSE)
                        {
                          /* Get which biomolecule this is */
                          nmatrix = biomol_ptr->nmatrix;

                          /* If within allowed range, store matrix */
                          if (nmatrix < MAX_BIOMT)
                            {
                              /* If this is a new transformation, store its
                                 starting position */
                              if (new_transform == TRUE)
                                {
                                  /* Store array position of starting
                                     transformation */
                                  ntransform = biomol_ptr->ntransform;
                                  biomol_ptr->first_matrix[ntransform]
                                    = nmatrix;

                                  /* Store list of chains corresponding to
                                     this set of transformations */
                                  store_string(&(biomol_ptr->tchains[ntransform]),
                                               chain_list);

                                  /* Increment transform number */
                                  biomol_ptr->ntransform++;
                                }

                              /* Loop through matrix elements to store */
                              for (j = 0; j < 3; j++)
                                for (k = 0; k < 3; k++)
                                  biomol_ptr->biomol_transform[nmatrix][j][k]
                                    = matrix[j][k];

                              /* Reset flag */
                              new_transform = FALSE;

                              /* Increment biomolecule number */
                              biomol_ptr->nmatrix++;
                            }

                          /* Otherwise, show warning */
                          else if (warned == FALSE)
                            {
                              printf("* Warning. Maximum number of biomol "
                                     "transformations exceeded. Some data "
                                     "lost\n");
                              warned = TRUE;
                            }
                        }

                      /* If have identity transformation, increment count */
                      else
                        biomol_ptr->nidentity++;
                    }
                }

              /* If not valid, show warning */
              else
                printf("* Warning. Invalid matrix row: %s\n",input_line);
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

  /* Determine whether this is an NMR structure or a crystallographic
     ensemble */
  if (ensemble == TRUE && nmr == FALSE &&
      params.pdb_type != ALPHA_FOLD_PDB)
    {
      /* If don't have a resolution, R-factor or R-free, can assume
         this to be an ensemble of NMR structures */
      if (data_item->resolution == 0.0 && data_item->rfactor == 0.0 &&
          data_item->rfree == 0.0)
        nmr = TRUE;
    }

  /* If have an NMR structure with zero models, make it 1 model */
  if (nmr == TRUE && nmodels == 0)
    nmodels = 1;

  /* Set flags for NMR or crystallographic ensemble */
  data_item->nmr = nmr;
  if (nmr == FALSE && ensemble == TRUE)
    data_item->ensemble = ensemble;
  data_item->nmodels = nmodels;

  /* Show number of atoms read in */
  printf("  Number of atoms read in:           %8d\n",nread);
  printf("  Number stored:                     %8d\n",natoms);
  printf("  Number of residues stored:         %8d\n",nresid);
  printf("  Number of 'chains':                %8d\n",nchains);
  printf("\n");

  /* Store total number of residues */
  data_item->nresidues = nresid;

  /* If had any CONECT records, then free memory taken by temporary
     atom pointers array */
  if (first_conect == FALSE)
    free(atom_index_ptr);

  /* Store the first atom- and residue-record pointers */
  *fst_atom_ptr = first_atom_ptr;
  *fst_chain_ptr = first_chain_ptr;
  *lst_chain_ptr = last_chain_ptr;
  *fst_hetdesc_ptr = first_hetdesc_ptr;
  *fst_residue_ptr = first_residue_ptr;

  /* Return number of residues stored */
  *nchn = nchains;
  *nmetyp = nmetal_types;
  *nmet = nmetals;
  *nres = nresid;
  *nsres = nsite_res;

  /* Return the number of atoms retained */
  return(natoms);
}
/***********************************************************************

get_residue_type  -  Get the residue type from its name

***********************************************************************/

void get_residue_type(struct residue *residue_ptr)
{
  int iamino, ires, namino, nbase, nonamino;

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
  namino = strlen(amino);
  nonamino = strlen(non_amino);
  nbase = strlen(base);

  /* Loop through all the residues until find a match */
  for (ires = 0; ires < namino && iamino == -1; ires++)
    if (!strcmp(residue_ptr->res_name,amino_name[ires]))
      {
        /* Found residue name, so retrieve single letter code and
           set as standard amino acid */
        iamino = ires;
        residue_ptr->aa_code = amino[iamino];
        residue_ptr->type = STANDARD;
      }

  /* If still not identified, repeat for the non-standard amino acids */
  for (ires = 0; ires < nonamino && iamino == -1; ires++)
    if (!strcmp(residue_ptr->res_name,non_amino_name[ires]))
      {
        /* Found residue name, so retrieve single letter code */
        iamino = ires;
        residue_ptr->aa_code = non_amino[iamino];
        residue_ptr->type = NON_STANDARD;
      }

  /* Loop through all the residues until find a match */
  for (ires = 0; ires < nbase && iamino == -1; ires++)
    if (!strcmp(residue_ptr->res_name,base_name[ires]))
      {
        /* Found residue name, so retrieve single letter code */
        iamino = ires;
        residue_ptr->aa_code = base[iamino];
        residue_ptr->type = NUCLEIC_ACID;
      }

  /* If still not found, check whether this might be a non-standard
     amino acid (ie one containing all three mainchain atoms N, CA
     and C) */
  if (iamino == -1 && (residue_ptr->nmain_chain == 3 ||
                       residue_ptr->have_calpha == TRUE))
    {
      residue_ptr->type = NON_STANDARD;
      iamino = 0;
    }

  /* Check whether this might be a metal ion */
  if (iamino == -1 && residue_ptr->natoms == 1)
    {
      if (residue_ptr->first_atom_ptr->type > NATOM_TYPES - 1)
        {
          /* Set flags to mark residue and chain as a single
             metal ion*/
          residue_ptr->type = METAL_ION;
        }
    }
}
/***********************************************************************

get_het_desc  -  Get Het Group description for given residue

***********************************************************************/

void get_het_desc(struct residue *residue_ptr,
                  struct hetdesc *first_hetdesc_ptr)
{
  int found;

  struct hetdesc *hetdesc_ptr;

  /* Initialise variables */
  found = FALSE;

  /* Get pointer to the first hetdesc record */
  hetdesc_ptr = first_hetdesc_ptr;

  /* Loop over all the hetdesc records */
  while (hetdesc_ptr != NULL && found == FALSE)
    {
      /* Check residue name for match */
      if (!strcmp(hetdesc_ptr->res_name,residue_ptr->res_name))
        {
          /* Store the pounter to this hetdesc record */
          residue_ptr->hetdesc_ptr = hetdesc_ptr;

          /* Set flag that we're found */
          found = TRUE;
        }

      /* Get the next hetdesc record */
      hetdesc_ptr = hetdesc_ptr->next_hetdesc_ptr;
    }
}
/***********************************************************************

update_residues -  Update chain pointers for all residues in the current
                   chain

***********************************************************************/

void update_residues(struct chain *chain_ptr)
{
  int iresid, nresid;

  struct residue *residue_ptr;

  /* Get the first of this chain's residues */
  residue_ptr = chain_ptr->first_residue_ptr;
  nresid = chain_ptr->nresid;
  iresid = 0;

  /* Loop through the chain's residues, determining to which
     other residues in this chain they are covalently bonded to */
  while (residue_ptr != NULL && iresid < nresid)
    {
      /* Set chain pointer */
      residue_ptr->chain_ptr = chain_ptr;

      /* Get the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;
    }
}
/***********************************************************************

excise_metal  -  Place metal ion as a separate chain

***********************************************************************/

void excise_metal(struct residue *residue_ptr,
                  struct chain **lst_chain_ptr)
{
  int iresid, nresid, nbefore;

  struct chain *chain_ptr, *new_chain_ptr;
  struct chain *last_chain_ptr;
  struct chain *first_dummy_chain_ptr, *last_dummy_chain_ptr;
  struct residue *other_residue_ptr;

  /* Initialise end chain pointer */
  last_chain_ptr = *lst_chain_ptr;

  /* Get this residue's chain */
  chain_ptr = residue_ptr->chain_ptr;

  /* If it is not the first residue of the chain start a new chain
     for it */
  if (residue_ptr != chain_ptr->first_residue_ptr)
    {
      /* Initialise count of preceding residues */
      other_residue_ptr = chain_ptr->first_residue_ptr;
      nresid = chain_ptr->nresid;
      iresid = 0;
      nbefore = 0;

      /* Loop through the residues preceding it in the current chain
         to determine how many there are */
      while (other_residue_ptr != NULL &&
             other_residue_ptr != residue_ptr && iresid < nresid)
        {
          /* Increment count */
          nbefore++;

          /* Get the next residue */
          other_residue_ptr = other_residue_ptr->next_residue_ptr;
          iresid++;
        }

      /* Initialise dummy pointers */
      first_dummy_chain_ptr = NULL;
      last_dummy_chain_ptr = NULL;

      /* Create a new chain record */
      new_chain_ptr = create_chain_record(&first_dummy_chain_ptr,
                                          &last_dummy_chain_ptr,
                                          chain_ptr->chain_id);

      /* Store current residue as the new chain's first */
      new_chain_ptr->first_residue_ptr = residue_ptr;

      /* Create links to the new chain */
      new_chain_ptr->next_chain_ptr = chain_ptr->next_chain_ptr;
      chain_ptr->next_chain_ptr = new_chain_ptr;

      /* If the current chain is the last in the linked list, make
         the new chain the new last */
      if (last_chain_ptr == chain_ptr)
        last_chain_ptr = new_chain_ptr;

      /* Update first residues and residue counts */
      new_chain_ptr->nresid = chain_ptr->nresid - nbefore;
      chain_ptr->nresid = nbefore;

      /* Copy across other parameters to the new chain */
      new_chain_ptr->type = chain_ptr->type;

      /* Update all this chain's residues so that they point to the
         new chain */
      update_residues(new_chain_ptr);

      /* Store new chain pointer */
      chain_ptr = new_chain_ptr;
    }

  /* Split off any residues after the current one, creating a new
     chain */
  if (chain_ptr->nresid > 1)
    {
      /* Initialise dummy pointers */
      first_dummy_chain_ptr = NULL;
      last_dummy_chain_ptr = NULL;

      /* Create a new chain record */
      new_chain_ptr = create_chain_record(&first_dummy_chain_ptr,
                                          &last_dummy_chain_ptr,
                                          chain_ptr->chain_id);

      /* Store residue after current residue as the new chain's first */
      new_chain_ptr->first_residue_ptr = residue_ptr->next_residue_ptr;

      /* Create links to the new chain */
      new_chain_ptr->next_chain_ptr = chain_ptr->next_chain_ptr;
      chain_ptr->next_chain_ptr = new_chain_ptr;

      /* If the current chain is the last in the linked list, make
         the new chain the new last */
      if (last_chain_ptr == chain_ptr)
        last_chain_ptr = new_chain_ptr;

      /* Update first residues and residue counts */
      new_chain_ptr->nresid = chain_ptr->nresid - 1;
      chain_ptr->nresid = 1;

      /* Copy across other parameters to the new chain */
      new_chain_ptr->type = chain_ptr->type;
      chain_ptr->type = METAL;

      /* Update all this chain's residues so that they point to the
         new chain */
      update_residues(new_chain_ptr);
    }

  /* Return last chain pointers */
  *lst_chain_ptr = last_chain_ptr;
}
/***********************************************************************

classify_residues  -  Classify all the residues according to type

***********************************************************************/

void classify_residues(struct residue *first_residue_ptr,
                       struct hetdesc *first_hetdesc_ptr,
                       struct chain **lst_chain_ptr)
{
  int i;

  struct chain *last_chain_ptr;
  struct residue *residue_ptr;

  /* Initialise end chain pointer */
  last_chain_ptr = *lst_chain_ptr;

  /* Get pointer to the first residue */
  residue_ptr = first_residue_ptr;

  /* Loop through all the residues to classify */
  while (residue_ptr != NULL)
    {
      /* Get the residue type and amino acid code from its name */
      get_residue_type(residue_ptr);

      /* If this a standard amino acid, increment count for this chain */
      if (residue_ptr->type == STANDARD)
        residue_ptr->chain_ptr->nstandard++;

      /* If residue non-standard, get its description */
      if (residue_ptr->type == NON_STANDARD ||
          residue_ptr->type == METAL_ION ||
          residue_ptr->type == UNKNOWN)
        get_het_desc(residue_ptr,first_hetdesc_ptr);

      /* If this residue is a metal ion, then want to place it in a
         separate chain of its own */
      if (residue_ptr->type == METAL_ION)
        excise_metal(residue_ptr,&last_chain_ptr);

      /* If have any sidechain atoms, compute the side-chain's centroid */
      if (residue_ptr->ncentroid_atoms > 0)
        {
          for (i = 0; i < 3; i++)
            residue_ptr->centroid[i]
              = residue_ptr->centroid[i] / residue_ptr->ncentroid_atoms;
        }

      /* Get the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
    }

  /* Return last chain pointers */
  *lst_chain_ptr = last_chain_ptr;
}
/***********************************************************************

shift_metals  -  Read in the given PDB file

***********************************************************************/

int shift_metals(struct chain **fst_chain_ptr,
                 struct chain **lst_chain_ptr)
{
  int nshift;

  struct chain *chain_ptr, *first_chain_ptr, *last_chain_ptr;
  struct chain *done_chain_ptr, *next_chain_ptr, *prev_chain_ptr;

  /* Initialise variables */
  done_chain_ptr = NULL;
  next_chain_ptr = NULL;
  prev_chain_ptr = NULL;
  nshift = 0;

  /* Initialise start and end chain pointers */
  first_chain_ptr = *fst_chain_ptr;
  last_chain_ptr = *lst_chain_ptr;

  /* Get pointer to the first chain */
  chain_ptr = first_chain_ptr;

  /* Loop over all the chains */
  while (chain_ptr != done_chain_ptr)
    {
      /* Get the next chain */
      next_chain_ptr = chain_ptr->next_chain_ptr;

      /* If this is a metal chain, shift to end of linked list */
      if (chain_ptr->type == METAL)
        {
          /* If have a previous chain, make its next pointer point
             around the metal */
          if (prev_chain_ptr != NULL)
            prev_chain_ptr->next_chain_ptr = next_chain_ptr;

          /* Otherwise, this must be the first chain, so adjust
             pointer to first chain */
          else
            first_chain_ptr = next_chain_ptr;

          /* Place current metal chain at end of linked list */
          if (last_chain_ptr != NULL)
            {
              /* Change pointers */
              last_chain_ptr->next_chain_ptr = chain_ptr;
              chain_ptr->next_chain_ptr = NULL;

              /* Make current chain the new last chain */
              last_chain_ptr = chain_ptr;
            }

          /* If this is the first shifted metal, save its pointer */
          if (done_chain_ptr == NULL)
            done_chain_ptr = chain_ptr;

          /* Increment count of chains shifted */
          nshift++;
        }

      /* If not metal, save as previous chain */
      else
        prev_chain_ptr = chain_ptr;

      /* Get the next chain */
      chain_ptr = next_chain_ptr;
    }

  /* Return current chain pointers */
  *fst_chain_ptr = first_chain_ptr;
  *lst_chain_ptr = last_chain_ptr;

  /* Return number of chains shifted */
  return(nshift);
}
/***********************************************************************

adjust_residues  -  Adjust the residues so that their order corresponds
                    to the current chain order

***********************************************************************/

int adjust_residues(struct residue **fst_residue_ptr,
                    struct chain *first_chain_ptr)
{
  int iresid, nres, nresid;

  struct chain *chain_ptr;
  struct residue *residue_ptr, *first_residue_ptr, *last_residue_ptr;

  /* Initialise variables */
  first_residue_ptr = NULL;
  last_residue_ptr = NULL;
  nres = 0;

  /* Get pointer to the first chain */
  chain_ptr = first_chain_ptr;

  /* Loop over the chains */
  while (chain_ptr != NULL)
    {
      /* Get this chain's first residue */
      residue_ptr = chain_ptr->first_residue_ptr;
      nresid = chain_ptr->nresid;
      iresid = 0;
  
      /* If this is the very first residue, store it */
      if (first_residue_ptr == NULL)
        first_residue_ptr = residue_ptr;

      /* If have a last residue, get it to point to this one */
      if (last_residue_ptr != NULL)
        last_residue_ptr->next_residue_ptr = residue_ptr;

      /* Initialise pointer to this chain's last residue */
      last_residue_ptr = NULL;

      /* Loop through the chain's residues to find the last */
      while (residue_ptr != NULL && iresid < nresid)
        {
          /* Update residue pointer */
          last_residue_ptr = residue_ptr;

          /* Get the next residue */
          residue_ptr = residue_ptr->next_residue_ptr;
          iresid++;
          nres++;
        }

      /* Get the next chain */
      chain_ptr = chain_ptr->next_chain_ptr;
    }

  /* Make last residue terminate the linked list */
  if (last_residue_ptr != NULL)
    last_residue_ptr->next_residue_ptr = NULL;

  /* Return current residue pointer */
  *fst_residue_ptr = first_residue_ptr;

  /* Return number of residues */
  return(nres);
}
/***********************************************************************

read_rin_file  -  Read in the secondary structure definitions from the
                  file.rin file

***********************************************************************/

void read_rin_file(char *pdb_code,char *pdbsum_dir,
                   struct residue *first_residue_ptr,int pdb_type)
{
  char input_line[LINELEN], file_name[FILENAME_LEN];
  char chain, last_chain, res_num[6], res_num1[6], res_num2[6];
  char res_name[4], res_name1[4], res_name2[4], secstr;

  int length, nread;
  int last_sec_struc;
  int found, have_file;

  struct residue *residue_ptr, *blank_residue_ptr, *start_residue_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  last_sec_struc = COIL;
  last_chain = '\0';
  length = 0;
  nread = 0;
  res_num1[0] = '\0';
  res_num2[0] = '\0';
  strcpy(res_name1,"XXX");
  strcpy(res_name2,"XXX");
  start_residue_ptr = NULL;

  /* Get pointer to the first residue */
  residue_ptr = first_residue_ptr;

  /* Form name of input .rin file */
  strcpy(file_name,"file.rin");
  if (pdb_type == PDBSUM1)
    sprintf(file_name,"../procheck/%s.rin",pdb_code);

  /* Open the .rin file */
  printf("Trying .rin file (%s) ...\n",file_name);
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    {
      /* File not found in current directory, so look into this
         PDB entry's PDBsum directory */
      strcpy(file_name,pdbsum_dir);
      if (pdb_type == PDBSUM1)
        strcat(file_name,"procheck/");
      strcat(file_name,"pdb");
      strcat(file_name,pdb_code);
      strcat(file_name,".rin");

      /* Open this file */
      printf("Trying .rin file (%s) ...\n",file_name);
      if ((file_ptr = fopen(file_name,"r")) ==  NULL)
        {
          printf("* Warning. Can't locate .rin file\n");
          have_file = FALSE;
        }
      else
        have_file = TRUE;
    }
  else
    have_file = TRUE;

  /* If have the file, read in all its data */
  if (have_file == TRUE)
    {
      /* Loop while reading in records from the file */
      while (fgets(input_line,LINELEN,file_ptr) != NULL &&
             residue_ptr != NULL)
        {
          /* Increment count of records read in */
          nread++;

          /* Get residue name */
          strncpy(res_name,input_line+4,3);
          res_name[3] = '\0';

          /* Get residue number */
          strncpy(res_num,input_line+9,5);
          res_num[5] = '\0';

          /* Get the chain-id */
          chain = input_line[8];

          /* Locate the corresponding residue */
          found = FALSE;
          while (found == FALSE && residue_ptr != NULL)
            {
              /* Check for match */
              if (!strcmp(res_num,residue_ptr->res_num) && 
                  !strcmp(res_name,residue_ptr->res_name) &&
                  chain == residue_ptr->chain)
                found = TRUE;

              /* If this isn't the right residue, move on to the
                 next residue */
              else
                residue_ptr = residue_ptr->next_residue_ptr;
            }

          /* If failed to find residue, show warning message */
          if (found == FALSE)
            printf("* Warning. Failed to locate residue %s %s %c "
                   "from .rin file\n",res_name,res_num,chain);

          /* Get the secondary structure and chain-id */
          secstr = input_line[14];

          /* Change secondary structure where necessary */
          if (secstr == 'g' || secstr == 'G' || secstr == 'i' ||
              secstr == 'I')
            secstr = 'H';
          if (secstr != 'H' && secstr != 'E')
            secstr = ' ';

          /* Store the secondary structure */
          if (residue_ptr != NULL && secstr != ' ')
            {
              if (secstr == 'H')
                residue_ptr->sec_struc = HELIX;
              else if (secstr == 'E')
                residue_ptr->sec_struc = STRAND;
              else
                residue_ptr->sec_struc = COIL;
            }
        }

      /* Close the file */
      fclose(file_ptr);

      /* Show number of residues read in */
      printf("  Number of residue records read:    %8d\n",nread);

      /* Get pointer to the first residue again */
      residue_ptr = first_residue_ptr;

      /* Loop through all the residues to process runs of secondary
         structure assignments */
      while (residue_ptr != NULL)
        {
          /* If either secondary structure or chain have changed,
             check for end of old run and/or start of new run */
          if (residue_ptr->sec_struc != last_sec_struc || chain != last_chain)
            {
              /* If previous secondary structure was a helix or strand
                 of fewer than 3 residues in length, need to rub it out */
              if (last_sec_struc != COIL && length < 3 &&
                  start_residue_ptr != NULL)
                {
                  /* Get pointer to start of the offending segment */
                  blank_residue_ptr = start_residue_ptr;

                  /* Loop over the previous residues to rub out secondary
                     structure */
                  while (blank_residue_ptr != residue_ptr &&
                         blank_residue_ptr != NULL)
                    {
                      /* Blank out this secondary structure */
                      blank_residue_ptr->sec_struc = COIL;

                      /* Get the next residue */
                      blank_residue_ptr
                        = blank_residue_ptr->next_residue_ptr;
                    }
                }

              /* If have secondary structure, start a new entry */
              if (residue_ptr->sec_struc != COIL)
                {
                  /* Store the start residue */
                  start_residue_ptr = residue_ptr;

                  /* Initialise length of secondary structure element */
                  length = 0;
                }

              /* Store current secondary structure and chain-id */
              last_sec_struc = residue_ptr->sec_struc;
              last_chain = chain;
            }

          /* If have secondary structure details, then increment length
             of current run */
          if (residue_ptr->sec_struc != COIL)
            length++;

          /* Get the next residue */
          residue_ptr = residue_ptr->next_residue_ptr;
        }
    }
}
/***********************************************************************

get_sulphur  -  Locate the sulphur atom in the given residue

***********************************************************************/

struct atom *get_sulphur(struct residue *residue_ptr)
{
  int iatom, natoms;

  struct atom *atom_ptr, *sulphur_atom_ptr;

  /* Initialise variables */
  sulphur_atom_ptr = NULL;

  /* Get number of atoms in this residue */
  natoms = residue_ptr->natoms;

  /* Initialise pointer to first atom of this residue */
  atom_ptr = residue_ptr->first_atom_ptr;
  iatom = 0;
  
  /* Loop through all the residue's atoms to look for sulphur */
  while (iatom < natoms && atom_ptr != NULL && sulphur_atom_ptr == NULL)
    {
      /* If this is the sulphur, then store */
      if (!strcmp(atom_ptr->atom_name," SG "))
        sulphur_atom_ptr = atom_ptr;

      /* Get pointer to next atom in linked-list */
      atom_ptr = atom_ptr->next_atom_ptr;
      iatom++;
    }

  /* Return the sulphur atom */
  return(sulphur_atom_ptr);
}
/***********************************************************************

create_ssbond_record  -  Create and initialise a new ssbond record

***********************************************************************/

struct ssbond *create_ssbond_record(struct ssbond **fst_ssbond_ptr,
                                    struct ssbond **lst_ssbond_ptr,
                                    struct residue *residue1_ptr,
                                    struct residue *residue2_ptr,
                                    struct atom *atom1_ptr,
                                    struct atom *atom2_ptr)
{
  struct ssbond *ssbond_ptr, *first_ssbond_ptr, *last_ssbond_ptr;

  /* Initialise ssbond pointers */
  ssbond_ptr = NULL;
  first_ssbond_ptr = *fst_ssbond_ptr;
  last_ssbond_ptr = *lst_ssbond_ptr;

  /* Allocate memory for structure to hold ssbond info */
  ssbond_ptr = (struct ssbond *) malloc(sizeof(struct ssbond));
  if (ssbond_ptr == NULL)
    {
      printf("*** Can't allocate memory for struct ssbond\n");
      exit (1);
    }

  /* If this is the very first ssbond, then store pointer */
  if (first_ssbond_ptr == NULL)
    first_ssbond_ptr = ssbond_ptr;

  /* Add link from previous ssbond to the current one */
  if (last_ssbond_ptr != NULL)
    last_ssbond_ptr->next_ssbond_ptr = ssbond_ptr;
  last_ssbond_ptr = ssbond_ptr;

  /* Store the ssbond details */
  ssbond_ptr->residue1_ptr = residue1_ptr;
  ssbond_ptr->residue2_ptr = residue2_ptr;
  ssbond_ptr->atom1_ptr = atom1_ptr;
  ssbond_ptr->atom2_ptr = atom2_ptr;
  ssbond_ptr->next_ssbond_ptr = NULL;

  /* Tag residues as belonging to disulphide bonds */
  residue1_ptr->disulphide = TRUE;
  residue2_ptr->disulphide = TRUE;

  /* Return current pointers */
  *fst_ssbond_ptr = first_ssbond_ptr;
  *lst_ssbond_ptr = last_ssbond_ptr;

  /* Return new ssbond pointer */
  return(ssbond_ptr);
}
/***********************************************************************

get_disulphides  -  Get all the disulphide bonds in the structure

***********************************************************************/

int get_disulphides(char *pdb_code,char *pdbsum_dir,
                    struct ssbond **fst_ssbond_ptr,
                    struct residue *first_residue_ptr,int pdb_type)
{
  char file_name[FILENAME_LEN + 1], input_line[LINELEN + 1];
  char res_name1[4], res_num1[6], chain1;
  char res_name2[4], res_num2[6], chain2;

  int len, nfound, nssbonds;

  struct atom *atom1_ptr, *atom2_ptr;
  struct residue *residue_ptr, *residue1_ptr, *residue2_ptr;
  struct ssbond *ssbond_ptr, *first_ssbond_ptr, *last_ssbond_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  nssbonds = 0;

  /* Initialise ssbond pointers */
  ssbond_ptr = NULL;
  first_ssbond_ptr = NULL;
  last_ssbond_ptr = NULL;

  /* Form the file name of the .dis disulphide bonds file */
  strcpy(file_name,pdbsum_dir);
  strcat(file_name,"pdb");
  strcat(file_name,pdb_code);
  strcat(file_name,".dis");
  if (pdb_type == PDBSUM1)
    sprintf(file_name,"../procheck/%s.dis",pdb_code);
  printf("Reading .dis file (%s) ...\n",file_name);

  /* Open disulphides file */
  if ((file_ptr = fopen(file_name,"r")) != NULL)
    {
      /* If present, loop while reading in records */
      while (fgets(input_line,LINELEN,file_ptr)!= NULL)
        {
          /* Truncate the string to strip off the newline */
          len = string_truncate(input_line,LINELEN);

          /* Read in the first residue's details */
          strncpy(res_name1,input_line+7,3);
          res_name1[3] = '\0';
          strncpy(res_num1,input_line+1,5);
          res_num1[5] = '\0';
          chain1 = input_line[0];

          /* Repeat for the second residue */
          strncpy(res_name2,input_line+21,3);
          res_name2[3] = '\0';
          strncpy(res_num2,input_line+15,5);
          res_num2[5] = '\0';
          chain2 = input_line[14];

          /* Get pointer to the first residue */
          residue_ptr = first_residue_ptr;
          residue1_ptr = NULL;
          residue2_ptr = NULL;
          nfound = 0;

          /* Loop through all the residues, to locate these two */
          while (residue_ptr != NULL && nfound < 2)
            {
              /* Check for match to first residue */
              if (residue1_ptr == NULL &&
                  !strcmp(res_name1,residue_ptr->res_name) &&
                  !strcmp(res_num1,residue_ptr->res_num) &&
                  chain1 == residue_ptr->chain)
                {
                  /* Store residue pointer */
                  residue1_ptr = residue_ptr;

                  /* Increment count */
                  nfound++;
                }

              /* Check for match to first residue */
              if (residue2_ptr == NULL &&
                  !strcmp(res_name2,residue_ptr->res_name) &&
                  !strcmp(res_num2,residue_ptr->res_num) &&
                  chain2 == residue_ptr->chain)
                {
                  /* Store residue pointer */
                  residue2_ptr = residue_ptr;

                  /* Increment count */
                  nfound++;
                }

              /* Get the next residue */
              residue_ptr = residue_ptr->next_residue_ptr;
            }

          /* If found both residues, then create a ssbond record */
          if (nfound == 2 && residue1_ptr != NULL &&
              residue2_ptr != NULL)
            {
              /* Identify the two sulphur atoms involved */
              atom1_ptr = get_sulphur(residue1_ptr);
              atom2_ptr = get_sulphur(residue2_ptr);

              /* Create ssbond record */
              ssbond_ptr
                = create_ssbond_record(&first_ssbond_ptr,&last_ssbond_ptr,
                                       residue1_ptr,residue2_ptr,
                                       atom1_ptr,atom2_ptr);

              /* Increment count of disulphide bonds stored */
              nssbonds++;
            }
        }

      /* Close the input file */
      fclose(file_ptr);
    }

  /* Show number of disulphides read in */
  printf("  Number of disulphide bonds:        %8d\n",nssbonds);

  /* Return pointer to first bond */
  *fst_ssbond_ptr = first_ssbond_ptr;

  /* Return number of disulphides */
  return(nssbonds);
}
/***********************************************************************

classify_chains  -  Initial analysis of each chain to determine what it
                    might be (protein/DNA/ligand)

***********************************************************************/

int classify_chains(struct chain *first_chain_ptr)
{
  int iatom, iresid, last_sec_struc, niatoms, natoms, nchains, nresid;
  int have_non_blank;

  struct atom *atom_ptr;
  struct chain *chain_ptr;
  struct residue *residue_ptr;

  /* Initialise variables */
  have_non_blank = FALSE;
  last_sec_struc = COIL;
  nchains = 0;

  /* Get pointer to the first chain */
  chain_ptr = first_chain_ptr;

  /* Loop over all the currently defined chains */
  while (chain_ptr != NULL)
    {
      /* If this is a non-blank chain identifier, the set flag */
      if (chain_ptr->chain_id != ' ')
        have_non_blank = TRUE;

      /* Initialise counts */
      chain_ptr->nhelices = 0;
      chain_ptr->nstrands = 0;
      chain_ptr->nstand = 0;
      chain_ptr->nbase = 0;
      natoms = 0;

      /* Get this chain's first residue */
      residue_ptr = chain_ptr->first_residue_ptr;
      nresid = chain_ptr->nresid;
      iresid = 0;

      /* Initialise chain counts */
      chain_ptr->natom_recs = 0;
      chain_ptr->nhetatm_recs = 0;

      /* Loop through this chain's residues, updating all the counts */
      while (residue_ptr != NULL && iresid < nresid)
        {
          /* Increment total atom count for this chain */
          natoms = natoms + residue_ptr->natoms;

          /* Get the first of the first residue's atoms */
          atom_ptr = residue_ptr->first_atom_ptr;
          iatom = 0;
          niatoms = residue_ptr->natoms;
  
          /* Loop over the residue's atoms to update the chain counts */
          while (atom_ptr != NULL && iatom < niatoms)
            {
              /* Update ATOM or HETATM count */
              if (atom_ptr->rec_type == ATOM_REC)
                chain_ptr->natom_recs++;
              else
                chain_ptr->nhetatm_recs++;

              /* Get pointer to next atom in linked-list */
              atom_ptr = atom_ptr->next_atom_ptr;
              iatom++;
            }

          /* Update chain counts depending on this residue's type */
          if (residue_ptr->type == STANDARD ||
              residue_ptr->type == NON_STANDARD)
            chain_ptr->nstand++;
          else if (residue_ptr->type == NUCLEIC_ACID)
            chain_ptr->nbase++;

          /* Update secondary structure counts */
          if (residue_ptr->sec_struc != last_sec_struc)
            {
              /* Update chain's helix and strand counts */
              if (residue_ptr->sec_struc == HELIX)
                chain_ptr->nhelices++;
              else if (residue_ptr->sec_struc == STRAND)
                chain_ptr->nstrands++;

              /* Save the current secondary structure */
              last_sec_struc = residue_ptr->sec_struc;
            }

          /* Get the next residue */
          residue_ptr = residue_ptr->next_residue_ptr;
          iresid++;
        }

      /* Consider non-metal chains */
      if (chain_ptr->type != METAL)
        {
          /* Check for obvious ligands: no ATOM records, only HETATM ones */
          if (chain_ptr->natom_recs == 0)
            chain_ptr->type = LIGAND;

          /* If chain has a blank chain-id, and we already have see non-blank
             chain identifiers, then this is not a protein */
          else if (chain_ptr->chain_id == ' ' && have_non_blank == TRUE &&
                   chain_ptr->nhelices + chain_ptr->nstrands == 0)
            chain_ptr->type = LIGAND;

          /* If chain contains nucleic acids, then mark down as a nucleic
             acid chain */
          else if (chain_ptr->nbase > 0)
            chain_ptr->type = DNA;

          /* Otherwise, if contains standard, or unusual, amino acids, then
             mark down as a polypeptide (or all-CA) chain */
          else if (chain_ptr->nstand > 0)
            {
              /* Determine whether this looks like an all-CA chain */
              if (natoms < 2 * chain_ptr->nresid)
                chain_ptr->type = CALPHA_ONLY;

              /* Otherwise, mark this chain down as a polypeptide chain */
              else
                chain_ptr->type = PROTEIN;
            }

          /* If no standard amino acids, and no nucleic acids, then this
             must be a ligand chain */
          else
            chain_ptr->type = LIGAND;
        }

      /* Increment count of chains */
      nchains++;

      /* Get the next chain */
      chain_ptr = chain_ptr->next_chain_ptr;
    }

  /* Return number of chains */
  return(nchains);
}
/***********************************************************************

join_chains  -  Join any adjacent chains which seem to belong to the
                same polypeptide

***********************************************************************/

void join_chains(struct chain *first_chain_ptr,struct parameters params)
{
  int iresid, nresid;

  struct chain *chain_ptr, *next_chain_ptr;
  struct residue *residue_ptr, *last_residue_ptr;

  /* Get pointer to the first chain */
  chain_ptr = first_chain_ptr;
  if (chain_ptr != NULL)
    next_chain_ptr = chain_ptr->next_chain_ptr;

  /* Loop over all the currently defined chains */
  while (chain_ptr != NULL && next_chain_ptr != NULL)
    {
      /* If both chains have the same chain-id and are of the same
         type (ie protein, DNA or ligand) then join together. The
         second chain needs to be longer than standard ligand length */
      if (chain_ptr->chain_id == next_chain_ptr->chain_id &&
          ((chain_ptr->type == next_chain_ptr->type &&
            next_chain_ptr->nresid > params.ligand_length) ||
           (chain_ptr->type == PROTEIN &&
            next_chain_ptr->type == CALPHA_ONLY)))
        {
          /* Get the first of the second chain's residues */
          residue_ptr = next_chain_ptr->first_residue_ptr;
          nresid = next_chain_ptr->nresid;
          iresid = 0;

          /* Loop through the chain's residues, marking them as now
             belonging to the first chain */
          while (residue_ptr != NULL && iresid < nresid)
            {
              /* Update chain pointer */
              residue_ptr->chain_ptr = chain_ptr;

              /* Get the next residue */
              residue_ptr = residue_ptr->next_residue_ptr;
              iresid++;
            }

          /* Join the second chain onto the first */
          chain_ptr->nresid
            = chain_ptr->nresid + next_chain_ptr->nresid;
          chain_ptr->nstandard
            = chain_ptr->nstandard + next_chain_ptr->nstandard;
          chain_ptr->nhelices
            = chain_ptr->nhelices + next_chain_ptr->nhelices;
          chain_ptr->nstrands
            = chain_ptr->nstrands + next_chain_ptr->nstrands;
          chain_ptr->nstand
            = chain_ptr->nstand + next_chain_ptr->nstand;
          chain_ptr->nbase
            = chain_ptr->nbase + next_chain_ptr->nbase;

          /* Get elongated chain to point past the removed second
             fragment */
          chain_ptr->next_chain_ptr = next_chain_ptr->next_chain_ptr;

          /* Free up memory from the discarded chain */
          free(next_chain_ptr);

          /* Get the new next chain */
          next_chain_ptr = chain_ptr->next_chain_ptr;
        }

      /* Otherwise, go on to the next chain */
      else
        {
          chain_ptr = next_chain_ptr;
          if (chain_ptr != NULL)
            next_chain_ptr = chain_ptr->next_chain_ptr;
        }
    }
}
/***********************************************************************

reclassify_chains  -  Reclassify very short polypeptide chains or
                      fragments of DNA as ligands

***********************************************************************/

void reclassify_chains(struct chain *first_chain_ptr,
                       struct parameters params)
{
  struct chain *chain_ptr;

  /* Get pointer to the first chain */
  chain_ptr = first_chain_ptr;

  /* Loop over all the currently defined chains */
  while (chain_ptr != NULL)
    {
      /* If short polypeptide chain, class as a ligand - unless it is
         the very first chain */
      if (chain_ptr->type == PROTEIN)
        {
          /* If chain shorter than ligand-length cutoff, reclassify as
             ligand */
          if (chain_ptr->nresid <= params.ligand_length)
            chain_ptr->type = LIGAND;

          /* Otherwise, for longer sequences, check for how many are
             non-standard residues */
          // else
          //   {
          //     /* If more than half the residues are non-standard,
          //        classify as ligand */
          //     if (chain_ptr->nstandard < chain_ptr->nresid / 2)
          //       chain_ptr->type = LIGAND;
          //   }
        }

      /* Similarly, a fragment of DNA of fewer than 3 bases can be
         classed a ligand */
      else if (chain_ptr->type == DNA && chain_ptr->nresid < 3)
        chain_ptr->type = LIGAND;

      /* Get the next chain */
      chain_ptr = chain_ptr->next_chain_ptr;
    }
}
/***********************************************************************

check_for_bond  -  Determine whether the given two residues are
                   covalently bonded

***********************************************************************/

int check_for_bond(struct residue *residue_ptr,
                   struct residue *other_residue_ptr)
{
  int i, iatom, niatoms, jatom, njatoms;
  int bonded, toofar;

  double dist_sqrd;

  struct atom *atom_ptr, *other_atom_ptr;
  struct rbond *rbond_ptr;

  /* Initialise variables */
  bonded = FALSE;

  /* Get the first of the residue's CONECT records */
  rbond_ptr = residue_ptr->first_rbond_ptr;

  /* Check for a bond defined by the CONECT records */
  while (rbond_ptr != NULL && bonded == FALSE)
    {
      /* If this bond joins the current two residues, then we're done */
      if (rbond_ptr->to_residue_ptr == other_residue_ptr)
        bonded = TRUE;

      /* Get the next residue CONECT record */
      rbond_ptr = rbond_ptr->next_rbond_ptr;
    }

  /* If no existing bond found, then check residue limits to see if
     residues might be within covalent bonding distance */
  if (bonded == FALSE)
    {
      /* Check whether these residue are within reach */
      toofar = FALSE;
      for (i = 0; i < 3 && toofar == FALSE; i++)
        {
          if (residue_ptr->coord_min[i]
              - other_residue_ptr->coord_max[i] > TOO_FAR_DIST)
            toofar = TRUE;
          if (other_residue_ptr->coord_min[i]
              - residue_ptr->coord_max[i] > TOO_FAR_DIST)
            toofar = TRUE;
        }
      
      /* If residues within range calculate all atom-atom distances */
      if (toofar == FALSE)
        {
          /* Get the first of the first residue's atoms */
          atom_ptr = residue_ptr->first_atom_ptr;
          iatom = 0;
          niatoms = residue_ptr->natoms;
  
          /* Loop through all the residue's atoms */
          while (atom_ptr != NULL && iatom < niatoms && bonded == FALSE)
            {
              /* Check whether this atom within range of the other
                 residue */
              toofar = FALSE;
              for (i = 0; i < 3 && toofar == FALSE; i++)
                {
                  if (atom_ptr->coords[i]
                      - other_residue_ptr->coord_max[i] > TOO_FAR_DIST ||
                      other_residue_ptr->coord_min[i]
                      - atom_ptr->coords[i] > TOO_FAR_DIST)
                    toofar = TRUE;
                }

              /* If atom not too far from the other residue, check its
                 distance from all that residue's atoms */
              if (toofar == FALSE)
                {
                  /* Get the first of the other residue's atoms */
                  other_atom_ptr = other_residue_ptr->first_atom_ptr;
                  jatom = 0;
                  njatoms = other_residue_ptr->natoms;
  
                  /* Loop through all the residue's atoms */
                  while (jatom < njatoms && other_atom_ptr != NULL &&
                         bonded == FALSE &&
                         atom_ptr->alt == other_atom_ptr->alt)
                    {
                      /* Calculate the squared distance between these
                         two atoms */
                      dist_sqrd = 0.0;
                      for (i = 0; i < 3; i++)
                        dist_sqrd
                          = dist_sqrd
                          + (atom_ptr->coords[i]
                             - other_atom_ptr->coords[i])
                          * (atom_ptr->coords[i]
                             - other_atom_ptr->coords[i]);

                      /* If within covalent bonding distance, have a
                         bond between the two residues */
                      if (dist_sqrd < COVALENT_DIST_SQRD)
                        {
                          /* Create a bond between these two residues */
                          create_rbond_record(residue_ptr,
                                              other_residue_ptr);

                          /* Set flag that residues bonded */
                          bonded = TRUE;
                        }

                      /* Get pointer to next atom in linked-list */
                      other_atom_ptr = other_atom_ptr->next_atom_ptr;
                      jatom++;
                    }
                }

              /* Get pointer to next atom in linked-list */
              atom_ptr = atom_ptr->next_atom_ptr;
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

void find_joins(struct chain *chain_ptr)
{
  int iresid, jresid, nresid;
  int bonded;

  struct residue *residue_ptr, *other_residue_ptr;

  /* Get the first of this chain's residues */
  residue_ptr = chain_ptr->first_residue_ptr;
  nresid = chain_ptr->nresid;
  iresid = 0;

  /* Loop through the chain's residues, determining to which
     other residues in this chain they are covalently bonded to */
  while (residue_ptr != NULL && iresid < nresid)
    {
      /* Get pointer to the next residue */
      other_residue_ptr = residue_ptr->next_residue_ptr;
      jresid = iresid + 1;

      /* Loop through each of the chain's other residues
         to calculate residue connectivities */
      while (other_residue_ptr != NULL && jresid < nresid)
        {
          /* Determine whether these two residues are covalently
             bonded */
          bonded = check_for_bond(residue_ptr,other_residue_ptr);

          /* Get the next residue */
          other_residue_ptr = other_residue_ptr->next_residue_ptr;
          jresid++;
        }

      /* Get the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;
    }
}
/***********************************************************************

create_list_record  -  Create and initialise a new list record

***********************************************************************/

struct list *create_list_record(struct list **fst_list_ptr,
                                struct list **lst_list_ptr,
                                struct residue *residue_ptr)
{
  struct list *list_ptr, *first_list_ptr, *last_list_ptr;

  /* Initialise list pointers */
  list_ptr = NULL;
  first_list_ptr = *fst_list_ptr;
  last_list_ptr = *lst_list_ptr;

  /* Allocate memory for structure to hold list info */
  list_ptr = (struct list *) malloc(sizeof(struct list));
  if (list_ptr == NULL)
    {
      printf("*** Can't allocate memory for struct list\n");
      exit (1);
    }

  /* If this is the very first list, then store pointer */
  if (first_list_ptr == NULL)
    first_list_ptr = list_ptr;

  /* Add link from previous list to the current one */
  if (last_list_ptr != NULL)
    last_list_ptr->next_list_ptr = list_ptr;
  last_list_ptr = list_ptr;

  /* Store the list details */
  list_ptr->residue_ptr = residue_ptr;
  list_ptr->next_list_ptr = NULL;

  /* Return current pointers */
  *fst_list_ptr = first_list_ptr;
  *lst_list_ptr = last_list_ptr;

  /* Return new list pointer */
  return(list_ptr);
}
/***********************************************************************

free_list_memory  -  Free up memory used by residue list

***********************************************************************/

void free_list_memory(struct list *first_list_ptr)
{
  struct list *list_ptr, *next_list_ptr;

  /* Get first list */
  list_ptr = first_list_ptr;

  /* Loop through the list records until done */
  while (list_ptr != NULL)
    {
      /* Get pointer to the next list record */
      next_list_ptr = list_ptr->next_list_ptr;

      /* Free up memory taken by current list */
      free(list_ptr);

      /* Go to the next list */
      list_ptr = next_list_ptr;
    }
}
/***********************************************************************

identify_molecule  -  Identify all residues in this chain that are in
                      some way connected to the first residue

***********************************************************************/

void identify_molecule(struct chain *chain_ptr)
{
  int iresid, nresid;

  struct list *list_ptr;
  struct list *first_list_ptr, *last_list_ptr;
  struct rbond *rbond_ptr;
  struct residue *residue_ptr, *to_residue_ptr;

  /* Initialise pointers */
  list_ptr = NULL;
  first_list_ptr = NULL;
  last_list_ptr = NULL;

  /* Get the first of this chain's residues */
  residue_ptr = chain_ptr->first_residue_ptr;
  nresid = chain_ptr->nresid;
  iresid = 0;

  /* Loop through the chain's residues to initialise connected flag */
  while (residue_ptr != NULL && iresid < nresid)
    {
      /* Initialise flag */
      residue_ptr->connected = FALSE;

      /* Get the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;
    }

  /* Get the first residue again */
  residue_ptr = chain_ptr->first_residue_ptr;
  iresid = 0;
  
  /* Create a list record for the first residue */
  list_ptr
    = create_list_record(&first_list_ptr,&last_list_ptr,residue_ptr);

  /* Loop while there are list records on the stack */
  while (list_ptr != NULL)
    {
      /* Get the residue stored */
      residue_ptr = list_ptr->residue_ptr;

      /* Process only if not already marked */
      if (residue_ptr->connected == FALSE)
        {
          /* Mark this residue as connected */
          residue_ptr->connected = TRUE;

          /* Get its first residue bond */
          rbond_ptr = residue_ptr->first_rbond_ptr;

          /* Loop through the list to create a list record for each */
          while (rbond_ptr != NULL)
            {
              /* Get the connected residue */
              to_residue_ptr = rbond_ptr->to_residue_ptr;

              /* If this residue not already connected, add to list */
              if (to_residue_ptr->connected == FALSE &&
                  to_residue_ptr->chain_ptr == chain_ptr)
                {
                  /* Create a list record for the residue to which the current
                     one is bonded */
                  create_list_record(&first_list_ptr,&last_list_ptr,
                                     rbond_ptr->to_residue_ptr);
                }

              /* Get the next residue bond */
              rbond_ptr = rbond_ptr->next_rbond_ptr;
            }
        }

      /* Get the next list record */
      list_ptr = list_ptr->next_list_ptr;
    }

  /* Free up memory taken by the list pointers */
  free_list_memory(first_list_ptr);
}
/***********************************************************************

split_ligands  -  Calculate the connectivities between all residues in
                  each ligand chain to determine which are phyically
                  distinct molecules

***********************************************************************/

void split_ligands(struct chain *first_chain_ptr)
{
  int iresid, nresid, start_nresid;
  int done;

  struct chain *chain_ptr, *new_chain_ptr;
  struct chain *fst_chain_ptr, *lst_chain_ptr;
  struct residue *residue_ptr, *start_residue_ptr;

  /* Get pointer to the first chain */
  chain_ptr = first_chain_ptr;

  /* Loop over all the currently defined chains */
  while (chain_ptr != NULL)
    {
      /* Only process if this is a ligand chain */
      if (chain_ptr->type == LIGAND)
        {
          /* Determine which residues in the chain are connected by
             covalent bonds */
          find_joins(chain_ptr);

          /* Initialise flag */
          done = FALSE;

          /* Loop while splitting the chain up */
          while (done == FALSE)
            {
              /* Mark all residues that can be reached from the first
                 residue */
              identify_molecule(chain_ptr);

              /* Get the first of this chain's residues */
              residue_ptr = chain_ptr->first_residue_ptr;
              nresid = chain_ptr->nresid;
              iresid = 0;
              start_residue_ptr = NULL;

              /* Loop through the chain's residues to retain all those
                 flagged as connected */
              while (residue_ptr != NULL && iresid < nresid &&
                     start_residue_ptr == NULL)
                {
                  /* If not connected, need to terminate here */
                  if (residue_ptr->connected == FALSE)
                    {
                      /* Save this residue as being the first of
                         a new chain */
                      start_residue_ptr = residue_ptr;

                      /* Calculate how many residues the new chain
                         will have */
                      start_nresid = nresid - iresid;

                      /* Save the number of residues the truncated
                         chain has */
                      chain_ptr->nresid = iresid;
                    }

                  /* Get the next residue */
                  residue_ptr = residue_ptr->next_residue_ptr;
                  iresid++;
                }

              /* If any residues left over, create a new chain for them
                 and transfer residues to it */
              if (start_residue_ptr != NULL)
                {
                  /* Initialise pointers */
                  fst_chain_ptr = NULL;
                  lst_chain_ptr = NULL;

                  /* Create a new chain record */
                  new_chain_ptr
                    = create_chain_record(&fst_chain_ptr,
                                          &lst_chain_ptr,
                                          chain_ptr->chain_id);

                  /* Store the chain's first residue and number of
                     residues */
                  new_chain_ptr->first_residue_ptr = start_residue_ptr;
                  new_chain_ptr->nresid = start_nresid;
                  new_chain_ptr->type = chain_ptr->type;

                  /* Update all this chain's residues so that they
                     point to the new chain */
                  update_residues(new_chain_ptr);

                  /* Adjust chain pointers to insert the new chain
                     into the linked list */
                  new_chain_ptr->next_chain_ptr
                    = chain_ptr->next_chain_ptr;
                  chain_ptr->next_chain_ptr = new_chain_ptr;
                  chain_ptr = new_chain_ptr;
                }

              /* Otherwise, if none left, then we are done */
              else
                done = TRUE;
            }
        }

      /* Get the next chain */
      chain_ptr = chain_ptr->next_chain_ptr;
    }
}
/***********************************************************************

single_residue_ligands  -  Identify which single-residue ligands are
                           metals and which are actually residue
                           attachments

***********************************************************************/

void single_residue_ligands(struct chain *first_chain_ptr)
{
  int iresid, nresid;
  int bonded, found, is_metal;

  struct atom *atom_ptr;
  struct chain *chain_ptr, *other_chain_ptr;
  struct residue *residue_ptr, *other_residue_ptr;

  /* Get pointer to the first chain */
  chain_ptr = first_chain_ptr;

  /* Loop over all the currently defined chains */
  while (chain_ptr != NULL)
    {
      /* Only process if this is a single-residue ligand chain */
      if (chain_ptr->type == LIGAND && chain_ptr->nresid == 1)
        {
          /* Get the residue */
          residue_ptr = chain_ptr->first_residue_ptr;
          is_metal = FALSE;

          /* If only a single atom, check whether its is a metal */
          if (residue_ptr->natoms == 1)
            {
              /* Get the atom */
              atom_ptr = residue_ptr->first_atom_ptr;

              /* Check whether this atom is a metal */
              if (atom_ptr->type > NATOM_TYPES - 1)
                {
                  /* Set flags to mark residue and chain as a single
                     metal ion*/
                  residue_ptr->type = METAL_ION;
                  chain_ptr->type = METAL;
                }
            }

          /* Get first chain */
          other_chain_ptr = first_chain_ptr;

          /* Loop over all the protein and DNA chains to see if this
             might be a residue attachment */
          while (other_chain_ptr != NULL)
            {
              /* Only process if this is a protein, ligand or DNA chain
                 with the same chain-id as our current residue */
              if (other_chain_ptr != chain_ptr &&
                  (other_chain_ptr->type == PROTEIN ||
                   other_chain_ptr->type == DNA ||
                   other_chain_ptr->type == LIGAND) &&
                  other_chain_ptr->chain_id == residue_ptr->chain)
                {
                  /* Get the first residue */
                  other_residue_ptr = other_chain_ptr->first_residue_ptr;
                  nresid = other_chain_ptr->nresid;
                  iresid = 0;
                  found = FALSE;

                  /* Loop through the chain's residues to find a
                     match on the residue number */
                  while (other_residue_ptr != NULL && iresid < nresid &&
                         found == FALSE)
                    {
                      /* Check for match on residue number */
                      if (!strcmp(other_residue_ptr->res_num,
                                  residue_ptr->res_num))
                        {
                          /* Check for a covalent bond between these
                             two residues */
                          bonded
                            = check_for_bond(residue_ptr,other_residue_ptr);

                          /* If bonded, have a residue attachment, so set
                             the residue and chain types accordingly */
                          if (bonded == TRUE)
                            {
                              /* Set chain type */
                              chain_ptr->type = ATTACHMENT;

                              /* Store pointers between the two
                                 residue */
                              residue_ptr->attach_residue_ptr
                                = other_residue_ptr;
                              other_residue_ptr->attach_residue_ptr
                                = residue_ptr;
                            }

                          /* Set flag that found a match */
                          found = TRUE;
                        }

                      /* Get this chain's next residue */
                      other_residue_ptr
                        = other_residue_ptr->next_residue_ptr;
                      iresid++;
                    }
                }

              /* Get the next chain */
              other_chain_ptr = other_chain_ptr->next_chain_ptr;
            }
        }

      /* Get the next chain */
      chain_ptr = chain_ptr->next_chain_ptr;
    }
}
/***********************************************************************

calc_chain_limits  -  Calculate coordinate limits of each chain

***********************************************************************/

int calc_chain_limits(struct chain *first_chain_ptr,double *coord_min,
                      double *coord_max,int *havdna,int *cando_nucplot)
{
  int i, iatom, iresid, max_dna_len, natoms, nresid, nviews;
  int have_dna, have_protein;

  double xwidth, ywidth, zwidth;

  struct atom *atom_ptr;
  struct chain *chain_ptr;
  struct residue *residue_ptr;

  /* Initialise variables */
  *cando_nucplot = FALSE;
  have_dna = FALSE;
  have_protein = FALSE;
  max_dna_len = 0;
  nviews = 1;

  /* Get pointer to the first chain */
  chain_ptr = first_chain_ptr;

  /* Loop over all the currently defined chains */
  while (chain_ptr != NULL)
    {
      /* If this is a protein or DNA/RNA chain, set flag */
      if (chain_ptr->type == PROTEIN)
        have_protein = TRUE;
      else if (chain_ptr->type == DNA)
        {
          /* Set flag */
          have_dna = TRUE;

          /* Store maximum DNA length */
          if (chain_ptr->nresid > max_dna_len)
            max_dna_len = chain_ptr->nresid;
        }

      /* Get the first residue */
      residue_ptr = chain_ptr->first_residue_ptr;
      nresid = chain_ptr->nresid;
      iresid = 0;
      chain_ptr->natoms = 0;

      /* Store chain's first atom number */
      if (residue_ptr != NULL)
        chain_ptr->first_atom_number
          = residue_ptr->first_atom_ptr->atom_number;

      /* Store first residue's limits */
      for (i = 0; i < 3; i++)
        {
          chain_ptr->coord_min[i] = residue_ptr->coord_min[i];
          chain_ptr->coord_max[i] = residue_ptr->coord_max[i];
        }

      /* Loop through the chain's residues to adjust coordinate limits */
      while (residue_ptr != NULL && iresid < nresid)
        {
          /* Update the chain's coordinate limits */
          for (i = 0; i < 3; i++)
            {
              if (residue_ptr->coord_min[i] < chain_ptr->coord_min[i])
                chain_ptr->coord_min[i] = residue_ptr->coord_min[i];
              if (residue_ptr->coord_max[i] > chain_ptr->coord_max[i])
                chain_ptr->coord_max[i] = residue_ptr->coord_max[i];
            }

          /* Add residue's number of atoms to chain's total */
          chain_ptr->natoms = chain_ptr->natoms + residue_ptr->natoms;

          /* If last residue in this chain, get the atom number of
             its last atom */
          if (iresid == nresid - 1)
            {
              /* Get this residue's first atom */
              atom_ptr = residue_ptr->first_atom_ptr;
              iatom = 0;
              natoms = residue_ptr->natoms;
  
              /* Loop through all the atoms to get last one */
              while (iatom < natoms && atom_ptr != NULL)
                {
                  /* Store the atom number */
                  chain_ptr->last_atom_number = atom_ptr->atom_number;

                  /* Get pointer to next atom in linked-list */
                  atom_ptr = atom_ptr->next_atom_ptr;
                  iatom++;
                }
            }

          /* Get this chain's next residue */
          residue_ptr = residue_ptr->next_residue_ptr;
          iresid++;
        }

      /* If first chain, store its limits as overall coordinate limits */
      if (chain_ptr == first_chain_ptr)
        {
          /* Store chain's limits */
          for (i = 0; i < 3; i++)
            {
              coord_min[i] = chain_ptr->coord_min[i];
              coord_max[i] = chain_ptr->coord_max[i];
            }
        }

      /* Otherwise, update overall coordinate limits */
      else
        {
          /* Update the chain's coordinate limits */
          for (i = 0; i < 3; i++)
            {
              if (chain_ptr->coord_min[i] < coord_min[i])
                coord_min[i] = chain_ptr->coord_min[i];
              if (chain_ptr->coord_max[i] > coord_max[i])
                coord_max[i] = chain_ptr->coord_max[i];
            }
        }

      /* Get the next chain */
      chain_ptr = chain_ptr->next_chain_ptr;
    }

  /* Show overall coordinate limits */
  printf("Overall coordinate limits:\n");
  printf("    X:  %8.3f   %8.3f\n",coord_min[0],coord_max[0]);
  printf("    Y:  %8.3f   %8.3f\n",coord_min[1],coord_max[1]);
  printf("    Z:  %8.3f   %8.3f\n",coord_min[2],coord_max[2]);

  /* Determine number of views (1 or 3) of whole structure to show */
  if (have_dna == TRUE)
    nviews = 3;

  /* Determine relative extent of the structure in each coordinate
     direction */
  else
    {
      /* Calculate extent in each coordinate direction */
      xwidth = fabs(coord_max[0] - coord_min[0]);
      ywidth = fabs(coord_max[1] - coord_min[1]);
      zwidth = fabs(coord_max[2] - coord_min[2]);

      /* If z-width (which in the single image is the direction into
         the page) is significantly greater than either of the two
         visible directions, then want to show all three orthogonal
         views */
      if (zwidth> 1.25 * xwidth || zwidth > 1.25 * ywidth)
        nviews = 3;
    }

  /* If no protein and no DNA, only need 1 view */
  //if (have_protein == FALSE && have_dna == FALSE)
  //  nviews = 1;

  /* Return whether we have one of more DNA chains */
  *havdna = have_dna;

  /* Return whether can generate a NUCPLOT diagram */
  if (have_protein == TRUE && have_dna == TRUE &&
      max_dna_len < MAX_NUCPLOT_DNA_LEN)
    *cando_nucplot = TRUE;

  /* Return number of views to show */
  return(nviews);
}
/***********************************************************************

write_defs  -  Write out all the molecule definitions

***********************************************************************/

void write_defs(char *file_name,struct chain *first_chain_ptr,
                struct data *data_item)
{
  int iresid, nresid;

  struct chain *chain_ptr;
  struct residue *residue_ptr, *first_residue_ptr, *last_residue_ptr;

  FILE *file_ptr;

  /* Open the output file */
  file_ptr = open_output_file(file_name,TRUE);

  /* Get pointer to the first chain */
  chain_ptr = first_chain_ptr;

  /* Loop over all the currently defined chains */
  while (chain_ptr != NULL)
    {
      /* Get the first of this chain's residues */
      first_residue_ptr = chain_ptr->first_residue_ptr;
      residue_ptr = first_residue_ptr;
      nresid = chain_ptr->nresid;
      iresid = 0;

      /* Loop through all the residues to identify the last residue
         in this chain */
      while (residue_ptr != NULL && iresid < nresid)
        {
          /* Store current residue as the last so far */
          last_residue_ptr = residue_ptr;

          /* Get the next residue */
          residue_ptr = residue_ptr->next_residue_ptr;
          iresid++;
        }

      /* Write out this chain's type, number of residues and start and
         end residues */
      if (chain_ptr->type == LIGAND)
        {
          fprintf(file_ptr,"LIGAND  ");
          data_item->nligands++;
        }
      else if (chain_ptr->type == PROTEIN)
        {
          fprintf(file_ptr,"PROTEIN ");
          data_item->nprotein++;
        }
      else if (chain_ptr->type == DNA)
        {
          fprintf(file_ptr,"DNA     ");
          data_item->ndna++;
        }
      else if (chain_ptr->type == CALPHA_ONLY)
        {
          fprintf(file_ptr,"CALPHA  ");
          data_item->nca_only++;
        }
      else if (chain_ptr->type == METAL)
        {
          fprintf(file_ptr,"METAL   ");
          data_item->nmetals++;
        }
      else if (chain_ptr->type == ATTACHMENT)
        {
          fprintf(file_ptr,"ATTACH  ");
          data_item->nattach++;
        }
      else
        {
          fprintf(file_ptr,"??????  ");
          data_item->nunknown++;
        }
      fprintf(file_ptr,"%6d resid  ",chain_ptr->nresid);
      fprintf(file_ptr," %s %s %c",
              first_residue_ptr->res_name,
              first_residue_ptr->res_num,
              first_residue_ptr->chain);
      if (chain_ptr->nresid > 1)
        fprintf(file_ptr," to %s %s %c",
                last_residue_ptr->res_name,
                last_residue_ptr->res_num,
                last_residue_ptr->chain);
      fprintf(file_ptr,"\n");

      /* Increment chain count */
      data_item->nchains++;

      /* Get the next chain */
      chain_ptr = chain_ptr->next_chain_ptr;
    }

  /* Close the output file */
  fclose(file_ptr);
}
/***********************************************************************

check_pqs  -  Check whether we have a definition for this structure's
              biological unit in the PQS

***********************************************************************/

int check_pqs(char *pdb_code,struct biomol *biomol_ptr,
              struct c_file_data file_name_ptr)
{
  char input_line[LINELEN], file_name[FILENAME_LEN];

  char *string_ptr;

  int i, ipos, last_space, len;
  int have_biomol, found;

  FILE *file_ptr;

  /* Initialise variables */
  have_biomol = FALSE;
  found = FALSE;

  /* Form the file name of the PQS BIOLIST file */
  strcpy(file_name,file_name_ptr.other_dir[PQS_DIR]);
  strcat(file_name,DIR_SEP);
  strcat(file_name,"BIOLIST");

  /* Open file */
  if ((file_ptr = fopen(file_name,"r")) != NULL)
    {
      /* If present, loop while reading in records */
      while (fgets(input_line,LINELEN,file_ptr)!= NULL && found == FALSE)
        {
          /* Truncate the string to strip off the newline */
          len = string_truncate(input_line,LINELEN);

          /* Check whether this is our PDB code */
          if (!strncmp(input_line,pdb_code,4))
            {
              /* Set flag that entry in PQS */
              biomol_ptr->in_pqs = TRUE;

              /* Check whether this entry to be skipped */
              string_ptr = strstr(input_line,"SKIP");

              /* If not being skipped, determine whether this is a
                 split ASU, or a generated bilogical unit */
              if (string_ptr == NULL)
                {
                  /* Set flag that we have PQS data */
                  have_biomol = TRUE;

                  /* Check whether ASU to be split */
                  string_ptr = strstr(input_line,"SPLIT");

                  /* If to be split, then set flag */
                  if (string_ptr != NULL)
                    biomol_ptr->pqs_split_asu = TRUE;

                  /* Test to see if we have polymer type */
                  string_ptr = strstr(input_line,"MERIC");
                  if (string_ptr != NULL)
                    {
                      /* Find current position */
                      ipos = string_ptr - input_line;

                      /* Loop through line to store position of last
                         space prior to above string */
                      last_space = 0;
                      for (i = 0; i < ipos; i++)
                        if (input_line[i] == ' ')
                          last_space = i;

                      /* Extract quaternary structure type */
                      biomol_ptr->pqs_quaternary_structure
                        = get_qstruc(input_line+last_space+1,
                                     biomol_ptr->pqs_prefix,
                                     biomol_ptr->pqs_qstruc_desc);
                    }

                  /* Otherwise, if a complex, store description */
                  else
                    {
                      /* Extract the oligomerization state */
                      string_ptr = strstr(input_line,"COMPLEX");
                      if (string_ptr != NULL)
                        {
                          /* Check for virus complex */
                          if (strstr(input_line,"VIRUS-COMPLEX") != NULL)
                            {
                              strcpy(biomol_ptr->pqs_qstruc_desc,
                                     "VIRUS-COMPLEX");
                              biomol_ptr->pqs_quaternary_structure
                                = VIRUS_COMPLEX;
                            }

                          /* Otherwise, store as complex */
                          else
                            {
                              strcpy(biomol_ptr->pqs_qstruc_desc,"COMPLEX");
                              biomol_ptr->pqs_quaternary_structure
                                = COMPLEX;
                            }
                        }
                    }
                }

              /* Mark as entry to be skipped */
              else
                biomol_ptr->pqs_skip = TRUE;

              /* Set flag that we've found our code */
              found = TRUE;
            }
        }

      /* Close the input file */
      fclose(file_ptr);
    }

  /* If file missing, show warning message */
  else
    printf("* Warning. PQS BIOLIST file is missing\n");

  /* Return whether we have a biological unit */
  return(have_biomol);
}
/***********************************************************************

check_asa  -  Check reliability of PQS assignment, as given in the
              ASALIST file

***********************************************************************/

void check_asa(char *pdb_code,struct biomol *biomol_ptr,
               struct c_file_data file_name_ptr)
{
  char input_line[LINELEN], file_name[FILENAME_LEN];

  char *string_ptr;

  int i, ipos, last_space, len;
  int have_biomol, found;

  FILE *file_ptr;

  /* Initialise variables */
  found = FALSE;

  /* Form the file name of the PQS ASALIST file */
  strcpy(file_name,file_name_ptr.other_dir[PQS_DIR]);
  strcat(file_name,DIR_SEP);
  strcat(file_name,"ASALIST");

  /* Open file */
  if ((file_ptr = fopen(file_name,"r")) != NULL)
    {
      /* If present, loop while reading in records */
      while (fgets(input_line,LINELEN,file_ptr)!= NULL && found == FALSE)
        {
          /* Truncate the string to strip off the newline */
          len = string_truncate(input_line,LINELEN);

          /* Check whether this is our PDB code */
          if (!strncmp(input_line,pdb_code,4))
            {
              /* Check whether this entry is marked OK */
              string_ptr = strstr(input_line," OK ");

              /* If not OK, mark entry as one to skip */
              if (string_ptr == NULL)
                {
                  /* Set flags */
                  biomol_ptr->pqs_skip = TRUE;
                  biomol_ptr->pqs_quaternary_structure = NONE;
                }

              /* Set flag that we've found our code */
              found = TRUE;
            }
        }

      /* Close the input file */
      fclose(file_ptr);
    }

  /* If file missing, show warning message */
  else
    printf("* Warning. PQS ASALIST file is missing\n");
}
/***********************************************************************

check_authmol  -  Check whether this PDB entry has an author-defined
                  biomolecule in the authmol.dat file

***********************************************************************/

void check_authmol(char *pdb_code,struct biomol *biomol_ptr,
                   struct c_file_data file_name_ptr)
{
  char input_line[LINELEN], file_name[FILENAME_LEN];

  char *string_ptr;

  int i, ipos, last_space, len;
  int have_biomol, found;

  FILE *file_ptr;

  /* Initialise variables */
  found = FALSE;

  /* Form the file name of the authmol.dat file */
  strcpy(file_name,file_name_ptr.other_dir[PDBSUM_REAL_DATA_DIR]);
  strcat(file_name,DIR_SEP);
  strcat(file_name,"authmol.dat");

  /* Open file */
  if ((file_ptr = fopen(file_name,"r")) != NULL)
    {
      /* If present, loop while reading in records */
      while (fgets(input_line,LINELEN,file_ptr)!= NULL && found == FALSE)
        {
          /* Truncate the string to strip off the newline */
          len = string_truncate(input_line,LINELEN);

          /* If this is our PDB code retrieve the biomolecule definition */
          if (!strncmp(input_line,pdb_code,4))
            {
              /* Test to see if we have polymer type */
              string_ptr = strstr(input_line,"MERIC");
              if (string_ptr != NULL)
                {
                  /* Find current position */
                  ipos = string_ptr - input_line;

                  /* Loop through line to store position of last
                     space prior to above string */
                  last_space = 0;
                  for (i = 0; i < ipos; i++)
                    if (input_line[i] == ' ')
                      last_space = i;

                  /* Extract quaternary structure type */
                  biomol_ptr->author_quaternary_structure
                    = get_qstruc(input_line+last_space+1,
                                 biomol_ptr->author_prefix,
                                 biomol_ptr->author_qstruc_desc);
                }

              /* Otherwise, if a complex, store description */
              else
                {
                  /* Extract the oligomerization state */
                  string_ptr = strstr(input_line,"COMPLEX");
                  if (string_ptr != NULL)
                    {
                      /* Check for virus complex */
                      if (strstr(input_line,"VIRUS-COMPLEX") != NULL)
                        {
                          strcpy(biomol_ptr->author_qstruc_desc,
                                 "VIRUS-COMPLEX");
                          biomol_ptr->author_quaternary_structure
                            = VIRUS_COMPLEX;
                        }

                      /* Otherwise, store as complex */
                      else
                        {
                          strcpy(biomol_ptr->author_qstruc_desc,
                                 "COMPLEX");
                          biomol_ptr->author_quaternary_structure
                            = COMPLEX;
                        }
                    }

                  /* Otherwise, check for SKIP */
                  else
                    {
                      /* Check whether this entry to be skipped */
                      string_ptr = strstr(input_line,"SKIP");

                      /* If to be skipped, set flag */
                      if (string_ptr != NULL)
                        biomol_ptr->author_skip = TRUE;
                    }
                }

              /* Set flag that we've found our code */
              found = TRUE;
              biomol_ptr->have_author_info = TRUE;
            }
        }

      /* Close the input file */
      fclose(file_ptr);
    }
}
/***********************************************************************

check_biomol  -  Check whether we have a definition for this structure's
                 biological unit, whether from the PDB file or from PQS

***********************************************************************/

int check_biomol(struct biomol *biomol_ptr,struct data data_item,
                 struct c_file_data file_name_ptr,
                 struct parameters params)
{
  char qstruc[30], style[4];

  int itrans, itype, nmatrix, nmolecules;
  int have_biomol;

  /* Initialise variables */
  have_biomol = FALSE;

  /* Get number of protein/DNA molecules in the PDB file */
  nmolecules
    = data_item.nprotein + data_item.nca_only + data_item.ndna;
  biomol_ptr->nmolecules = nmolecules;

  /* Determine style of biomolecule transformation matrices */
  if (biomol_ptr->new_style == TRUE)
    strcpy(style,"NEW");
  else
    strcpy(style,"OLD");
  itype = biomol_ptr->quaternary_structure;
  if (itype > 0 && itype < 13)
    strcpy(qstruc,UNIT_TYPE[itype]);
  else
    sprintf(qstruc,"%dMER",itype);
  if (biomol_ptr->prefix[0] != '\0')
    {
      strcat(qstruc," prefix: ");
      strcat(qstruc,biomol_ptr->prefix);
    }

  /* Determine whether structure type means that we can have no biomolecule
     information */
  if (data_item.nmr == TRUE ||
      data_item.nprotein + data_item.nca_only == 0 ||
      data_item.theoretical_model == TRUE)
    biomol_ptr->skip = TRUE;

  /* Check the PQS entry for this PDB code, if present */
  if (params.pdb_type == STANDARD_PDB)
    {
      /* Check whether this PDB code has en entry in the author-defined
         biomolecules file, authmol.dat */
      check_authmol(params.pdb_code,biomol_ptr,file_name_ptr);

      /* If don't have author assignment, check PQS */
      if (biomol_ptr->have_author_info == FALSE)
        {
          /* Get PQS quaternary structure assignments from BIOLIST file */
          check_pqs(params.pdb_code,biomol_ptr,file_name_ptr);

          /* If have info, check the ASALIST file to determine whether
             it is reliable */
          if (biomol_ptr->in_pqs == TRUE && biomol_ptr->pqs_skip == FALSE
              && biomol_ptr->pqs_quaternary_structure != VIRUS_COMPLEX)
            check_asa(params.pdb_code,biomol_ptr,file_name_ptr);
        }

      /* Otherwise set flag that we have info on the biomolecule */
      else
        have_biomol = TRUE;
    }

  /* If have any transformations from the PDB file, show what they are */
  if (biomol_ptr->ntransform > 0)
    {
      /* Print the numbers of transforms and matrices */
      printf("PDB %s: Transforms %2d (%s) Matrices %2d  %s\n",
             params.pdb_code,biomol_ptr->ntransform,style,
             biomol_ptr->nmatrix,qstruc);

      /* Break down number of matrics by each transform */
      for (itrans = 0; itrans < biomol_ptr->ntransform; itrans++)
        {
          if (itrans == biomol_ptr->ntransform - 1)
            nmatrix
              = biomol_ptr->nmatrix - biomol_ptr->first_matrix[itrans];
          else
            nmatrix
              = biomol_ptr->first_matrix[itrans + 1]
              - biomol_ptr->first_matrix[itrans];
          printf("PDB %s: Transform %2d Matrices %2d  Chains: %s\n",
                 params.pdb_code,itrans + 1,nmatrix,
                 biomol_ptr->tchains[itrans]);
        }

      /* Set flag that have biological unit transformations */
      have_biomol = TRUE;
    }

  /* If no transforms given in the PDB file, print message */
  else
    {
      /* Only have identity operators, so if more than one, need to
         generate image of biological molecule */
      if (biomol_ptr->nidentity > 0)
        {
          /* Show number */
          printf("PDB %s: %2d identity transforms only (%s) %s\n",
                 params.pdb_code,biomol_ptr->nidentity,style,qstruc);

          /* If more than one, then have split asymmetric unit, so
             want to generate image for biomolecule */
          if (biomol_ptr->nidentity > 1)
            {
              /* Set flags */
              biomol_ptr->split_asu = TRUE;
              have_biomol = TRUE;
            }
        }

      /* If no transformations given */
      else
        {
          /* Show message */
          printf("PDB %s: No transforms given %s\n",params.pdb_code,
                 qstruc);
        }
    }

  /* Show quaternary structure definition from PDB file, if any */
  printf("  PDB %s: ",params.pdb_code);
  printf("P=%d C=%d D=%d ",data_item.nprotein,data_item.nca_only,
         data_item.ndna);

  /* If have biomolecule info, show what it is */
  if (biomol_ptr->skip == FALSE)
    {
      /* Show the quaternary structure, as reported in the PDB file */
      if (biomol_ptr->quaternary_structure != NOT_GIVEN)
        {
          printf("%3d ",biomol_ptr->quaternary_structure);
          if (biomol_ptr->prefix[0] != '\0')
            printf("%s",biomol_ptr->prefix);
          if (biomol_ptr->qstruc_desc[0] != '\0')
            printf("%s",biomol_ptr->qstruc_desc);
        }
      else
        printf("NOT GIVEN");

      /* If have any transforms, show numbers of matrices */
      if (biomol_ptr->ntransform > 0)
        {
          if (biomol_ptr->split_asu == TRUE)
            printf(" split ASU");
          printf("  Transforms %3d  Matrices %3d",
                 biomol_ptr->ntransform,biomol_ptr->nmatrix);
          if (biomol_ptr->nsymmetry > 0)
            printf("  Symmetry %3d",biomol_ptr->nsymmetry);
          else
            printf("              ");
          if (biomol_ptr->nbiomolecule > 0)
            printf("  Biomol %3d",biomol_ptr->nbiomolecule);
        }
    }
  else
    printf("SKIP");

  /* Close off this line */
  printf("\n");

  /* If have PQS data, show that */
  if (biomol_ptr->in_pqs == TRUE)
    {
      /* Print PDB code */
      printf("  PQS %s: ",params.pdb_code);

      /* If have biomolecule info, show what it is */
      if (biomol_ptr->pqs_skip == FALSE)
        {
          /* Show the quaternary structure, as reported in PQS */
          if (biomol_ptr->pqs_quaternary_structure != NOT_GIVEN)
            {
              printf("%3d ",biomol_ptr->pqs_quaternary_structure);
              if (biomol_ptr->pqs_prefix[0] != '\0')
                printf("%s",biomol_ptr->pqs_prefix);
              if (biomol_ptr->pqs_qstruc_desc[0] != '\0')
                printf("%s",biomol_ptr->pqs_qstruc_desc);
              printf("\n");

              /* Compute number of symmetry operations this corresponds
                 to */
              if (nmolecules > 0)
                biomol_ptr->nsymop
                  = biomol_ptr->pqs_quaternary_structure / nmolecules;

              /* Set flag that have the biomolecule */
              have_biomol = TRUE;
            }
          else
            printf("NOT GIVEN");

          /* Show whether this is a split ASU */
          if (biomol_ptr->pqs_split_asu == TRUE)
            printf(" split ASU");
        }
      else
        printf("SKIP");

      /* Close off this line */
      printf("\n");
    }

  /* If ASU corresponds to the biomolecule, then already have the
     biological molecule */
  if (nmolecules == biomol_ptr->quaternary_structure ||
      nmolecules == biomol_ptr->pqs_quaternary_structure ||
      nmolecules == biomol_ptr->author_quaternary_structure)
    biomol_ptr->is_biomolecule = TRUE;

  /* Otherwise determine whether we need to use PQS coords */
  else if (biomol_ptr->in_pqs == TRUE && biomol_ptr->pqs_skip == FALSE)
    biomol_ptr->use_pqs = TRUE;

  /* Return whether we have a biological unit */
  return(have_biomol);
}
/***********************************************************************

create_hbadd_record  -  Create and initialise a new hbadd record

***********************************************************************/

struct hbadd *create_hbadd_record(struct hbadd **fst_hbadd_ptr,
                                  struct hbadd **lst_hbadd_ptr,
                                  struct residue *residue_ptr,int nmissed)
{
  struct hbadd *hbadd_ptr, *first_hbadd_ptr, *last_hbadd_ptr;

  /* Initialise hbadd pointers */
  hbadd_ptr = NULL;
  first_hbadd_ptr = *fst_hbadd_ptr;
  last_hbadd_ptr = *lst_hbadd_ptr;

  /* Allocate memory for structure to hold hbadd info */
  hbadd_ptr = (struct hbadd *) malloc(sizeof(struct hbadd));
  if (hbadd_ptr == NULL)
    {
      printf("*** Can't allocate memory for struct hbadd\n");
      exit (1);
    }

  /* If this is the very first hbadd, then store pointer */
  if (first_hbadd_ptr == NULL)
    first_hbadd_ptr = hbadd_ptr;

  /* Add link from previous hbadd to the current one */
  if (last_hbadd_ptr != NULL)
    last_hbadd_ptr->next_hbadd_ptr = hbadd_ptr;
  last_hbadd_ptr = hbadd_ptr;

  /* Store the hbadd details */
  hbadd_ptr->residue_ptr = residue_ptr;
  hbadd_ptr->nmissed = nmissed;
  hbadd_ptr->next_hbadd_ptr = NULL;

  /* Return current pointers */
  *fst_hbadd_ptr = first_hbadd_ptr;
  *lst_hbadd_ptr = last_hbadd_ptr;

  /* Return new hbadd pointer */
  return(hbadd_ptr);
}
/***********************************************************************

read_hbadd_sum  -  Read in all the Het Groups identified by the hbadd
                   program

***********************************************************************/

void read_hbadd_sum(char *pdb_code,char *pdbsum_dir,
                    struct chain *first_chain_ptr,
                    struct hbadd **fst_hbadd_ptr,int pdb_type)
{
  char input_line[LINELEN], file_name[FILENAME_LEN];
  char chain, res_num[6], res_name[4];

  int iresid, nresid;
  int nmiss, nmissed, nmissing;
  int found;

  struct chain *chain_ptr;
  struct hbadd *hbadd_ptr;
  struct hbadd *first_hbadd_ptr, *last_hbadd_ptr;
  struct residue *residue_ptr;

  FILE *file_ptr;

  /* Initialise hbadd pointers */
  hbadd_ptr = NULL;
  last_hbadd_ptr = NULL;
  first_hbadd_ptr = NULL;

  /* Form name of input hbadd.sum file */
  strcpy(file_name,pdbsum_dir);
  if (pdb_type == PDBSUM1)
    strcat(file_name,"hbplus/");
  strcat(file_name,"hbadd.sum");

  /* Open this file */
  printf("Opening hbadd.sum file (%s) ...\n",file_name);
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    printf("* Warning. Can't locate input hbadd.sum file\n");

  /* If have the file, read in all its data */
  else
    {
      /* Loop while reading in records from the file */
      while (fgets(input_line,LINELEN,file_ptr))
        {
          /* Get residue name */
          strncpy(res_name,input_line,3);
          res_name[3] = '\0';

          /* Get the chain-id */
          chain = input_line[4];

          /* Get residue number */
          strncpy(res_num,input_line+5,5);
          res_num[5] = '\0';

          /* Get number of missed and missing atoms */
          sscanf(input_line+38," %d %d",&nmiss,&nmissing);
          nmissed = nmiss + nmissing;

          /* Get pointer to the first chain */
          chain_ptr = first_chain_ptr;
          found = FALSE;

          /* Loop over all the chains */
          while (chain_ptr != NULL && found == FALSE)
            {
              /* Only consider if this is the right chain identifier */
              if (chain == chain_ptr->chain_id)
                {
                  /* Get this chain's first residue */
                  residue_ptr = chain_ptr->first_residue_ptr;
                  nresid = chain_ptr->nresid;
                  iresid = 0;

                  /* Loop through this chain's residues to find match */
                  while (residue_ptr != NULL && iresid < nresid &&
                         found == FALSE)
                    {
                      /* Check for match */
                      if (!strcmp(res_num,residue_ptr->res_num) && 
                          !strcmp(res_name,residue_ptr->res_name) &&
                          chain == residue_ptr->chain)
                        {
                          /* Create a new hbadd record to store this 
                             residue */
                          hbadd_ptr
                            = create_hbadd_record(&first_hbadd_ptr,
                                                  &last_hbadd_ptr,
                                                  residue_ptr,nmissed);

                          /* Set flag that residue found */
                          found = TRUE;
                        }

                      /* Get the next residue */
                      residue_ptr = residue_ptr->next_residue_ptr;
                      iresid++;
                    }
                }

              /* Get the next chain */
              chain_ptr = chain_ptr->next_chain_ptr;
            }
        }

      /* Close the file */
      fclose(file_ptr);
    }

  /* Return pointer to the first of the hbadd records */
  *fst_hbadd_ptr = first_hbadd_ptr;
}
/***********************************************************************

check_het_group  -  Check for any existing representative of the
                    current Het Group

***********************************************************************/

int check_het_group(char *res_name,int nmissed,char *hetgif_dir)
{
  char file_name[FILENAME_LEN], input_line[LINELEN];
  char first_char[2], number_string[20];

  int len, nmissing;
  int best;

  FILE *file_ptr;

  /* Initialise variables */
  best = FALSE;
  nmissing = 999;

  /* Get the first character */
  first_char[0] = res_name[0];
  first_char[1] = '\0';

  /* Form the file name for this Het Group's .en file */
  strcpy(file_name,hetgif_dir);
  strcat(file_name,DIR_SEP);
  strcat(file_name,first_char);
  strcat(file_name,DIR_SEP);
  strcat(file_name,res_name);
  strcat(file_name,DIR_SEP);
  strcat(file_name,"het");
  strcat(file_name,res_name);
  strcat(file_name,".en");

  /* Open input .en file */
  if ((file_ptr = fopen(file_name,"r")) !=  NULL)
    {
      /* Loop while reading in records from the input file */
      while (fgets(input_line,LINELEN,file_ptr)!= NULL)
        {
          /* Truncate the string to strip off the newline */
          len = string_truncate(input_line,LINELEN);

          /* Read in the score associated with this Het Group */
          strcpy(number_string,input_line+17);
          nmissing = atoi(number_string);
        }

      /* If current Het Group has a better match to the dictionary
         than the representation on disk, then it is the best so
         far */
      if (nmissed < nmissing)
        best = TRUE;

      /* Close the file */
      fclose(file_ptr);
    }

  /* If no file, then our current Het Group is the best (and only)
     do far, so need to genherate the files */
  else
    best = TRUE;

  /* Return whether this is the best example so far */
  return(best);
}
/***********************************************************************

write_en_file  -  Write out the .en file for the current Het Group

***********************************************************************/

void write_en_file(char *res_name,char *pdb_code,
                   struct residue *residue_ptr,int nmissed,
                   char *output_dir)
{
  char file_name[FILENAME_LEN];

  FILE *file_ptr;

  /* Form the file name for this Het Group's .en file */
  strcpy(file_name,output_dir);
  strcat(file_name,"het");
  strcat(file_name,res_name);
  strcat(file_name,".en");

  /* Open the output file */
  file_ptr = open_output_file(file_name,TRUE);

  /* Write out the details of this Het Group */
  fprintf(file_ptr,"%s %s %s %c %6d\n",pdb_code,residue_ptr->res_name,
          residue_ptr->res_num,residue_ptr->chain,nmissed);

  /* Close the file */
  fclose(file_ptr);
}
/***********************************************************************

centre_of_mass  -  Find the centre of mass of the supplied coordinates

***********************************************************************/

void centre_of_mass(struct atom *first_atom_ptr,int natoms,
                    double *centre_x,double *centre_y,double *centre_z)
{
  int iatom;

  double centre_mass_x, centre_mass_y, centre_mass_z;
  double mass_x = 0.0, mass_y = 0.0, mass_z = 0.0;

  struct atom *atom_ptr;

  /* Get pointer to the first atom */
  atom_ptr = first_atom_ptr;
  iatom = 0;

  /* Loop through all the atoms */
  while (atom_ptr != NULL && iatom < natoms)
    {
      /* Add coordinates to centre of mass */
      mass_x = mass_x + atom_ptr->coords[0];
      mass_y = mass_y + atom_ptr->coords[1];
      mass_z = mass_z + atom_ptr->coords[2];

      /* Get pointer to the next atom */
      atom_ptr = atom_ptr->next_atom_ptr;
      iatom++;
    }

  /* Calculate centre of mass */
  centre_mass_x = mass_x / natoms;
  centre_mass_y = mass_y / natoms;
  centre_mass_z = mass_z / natoms;
  *centre_x = centre_mass_x;
  *centre_y = centre_mass_y;
  *centre_z = centre_mass_z;
}
/***********************************************************************

centre_of_mass2  -  Find the centre of mass of the two sets of supplied
                    coordinates

***********************************************************************/

void centre_of_mass2(struct atom *first_atom1_ptr,int natoms1,
                     struct atom *first_atom2_ptr,int natoms2,
                     double *centre_x,double *centre_y,double *centre_z)
{
  int iatom, loop, natoms, ntotal;

  double centre_mass_x, centre_mass_y, centre_mass_z;
  double mass_x = 0.0, mass_y = 0.0, mass_z = 0.0;

  struct atom *atom_ptr;

  /* Initialise overall atom count */
  ntotal = natoms1 + natoms2;

  /* Loop over the two sets of atoms */
  for (loop = 0; loop < 2; loop++)
    {
      /* Get pointer to the first atom */
      if (loop == 0)
        {
          atom_ptr = first_atom1_ptr;
          natoms = natoms1;
        }
      else
        {
          atom_ptr = first_atom2_ptr;
          natoms = natoms2;
        }
      iatom = 0;

      /* Loop through all the atoms */
      while (atom_ptr != NULL && iatom < natoms)
        {
          /* Add coordinates to centre of mass */
          mass_x = mass_x + atom_ptr->coords[0];
          mass_y = mass_y + atom_ptr->coords[1];
          mass_z = mass_z + atom_ptr->coords[2];

          /* Get pointer to the next atom */
          atom_ptr = atom_ptr->next_atom_ptr;
          iatom++;
        }
    }

  /* Calculate centre of mass */
  centre_mass_x = mass_x / ntotal;
  centre_mass_y = mass_y / ntotal;
  centre_mass_z = mass_z / ntotal;
  *centre_x = centre_mass_x;
  *centre_y = centre_mass_y;
  *centre_z = centre_mass_z;
}
/***********************************************************************

adjust_to_origin  -  Translate all the atoms by the given amount

***********************************************************************/

void adjust_to_origin(struct atom *first_atom_ptr,int natoms,
                      double centre_x,double centre_y,double centre_z,
                      int store)
{
  int iatom;

  struct atom *atom_ptr;

  /* Get pointer to the first atom */
  atom_ptr = first_atom_ptr;
  iatom = 0;

  /* Loop through all the atoms */
  while (atom_ptr != NULL && iatom < natoms)
    {
      /* Apply the transformation */
      atom_ptr->view_coords[0] = atom_ptr->coords[0] - centre_x;
      atom_ptr->view_coords[1] = atom_ptr->coords[1] - centre_y;
      atom_ptr->view_coords[2] = atom_ptr->coords[2] - centre_z;

      /* If change to be permanent, store the new coords */
      if (store == TRUE)
        {
          atom_ptr->coords[0] = atom_ptr->view_coords[0];
          atom_ptr->coords[1] = atom_ptr->view_coords[1];
          atom_ptr->coords[2] = atom_ptr->view_coords[2];
        }

      /* Get pointer to the next atom */
      atom_ptr = atom_ptr->next_atom_ptr;
      iatom++;
    }
}
/***********************************************************************
  
calc_eigen_values  -  Subroutine to compute eigenvalues & eigenvectors
                      of a real symmetric matrix, from IBM SSP manual
                      (see p165).

                      Ian Tickle April 1992
 
                      (modified by David Moss February 1993)

Eigenvalues & vectors of real symmetric matrix stored in triangular form.
 
Arguments:
 
mv = 0 to compute eigenvalues only.
mv = 1 to compute eigenvalues & eigenvectors.
 
n - supplied dimension of n*n matrix.
 
a - supplied array of size n*(n+1)/2 containing n*n matrix in the order:
 
        1      2      3    ...
     1    a[0]
     2    a[1]   a[2]
     3    a[3]   a[4]   a[5]   ...
     ...

     NOTE a is used as working space and is overwritten on return.
     Eigenvalues are written into diagonal elements of a
     i.e.  a[0]  a[2]  a[5]  for a 3*3 matrix.
 
r - Resultant matrix of eigenvectors stored columnwise in the same
    order as eigenvalues.

***********************************************************************/

void calc_eigen_values(int mv, int n, double* a, double* r)
{
  int i, il, ilq, ilr, im, imq, imr, ind, iq, j, l, ll, lm, lq, m, mm, mq;
  double anorm, anrmx, cosx, cosx2, sincs, sinx, sinx2, thr, x, y;

  if (mv)
    { for (i=1; i<n*n; i++) r[i]=0.;
    for (i=0; i<n*n; i+=n+1) r[i]=1.;
    }

  /* Initial and final norms (anorm & anrmx). */
  anorm=0.;
  iq=0;
  for (i=0; i<n; i++) for (j=0; j<=i; j++)
    { if (j!=i) anorm+=a[iq]*a[iq];
    ++iq;
    }

  if (anorm>0.)
    { anorm=sqrt(2.0*anorm);
    anrmx=(1.0e-6*anorm/n);

    /* Compute threshold and initialise flag. */
    thr=anorm;
    do
      { thr/=n;
      do
        { ind=0;
        l=0;
        do
          { lq=l*(l+1)/2;
          ll=l+lq;
          m=l+1;
          ilq=n*l;
          do

            /* Compute sin & cos. */
            { mq=m*(m+1)/2;
            lm=l+mq;
            if (fabs(a[lm])>=thr)
              { ind=1;
              mm=m+mq;
              x=0.5*(a[ll]-a[mm]);
              y=-a[lm]/sqrt(a[lm]*a[lm]+x*x);
              if (x<0.0) y=-y;
              sinx=y/sqrt(2.0*(1.0+(sqrt(1.0-y*y))));
              sinx2=sinx*sinx;
              cosx=sqrt(1.0-sinx2);
              cosx2=cosx*cosx;
              sincs=sinx*cosx;

              /* Rotate l & m columns. */
              imq=n*m;
              for (i=0; i<n; i++)
                { iq=i*(i+1)/2;
                if (i!=l && i!=m)
                  { if (i<m) im=i+mq;
                  else im=m+iq;
                  if (i<l) il=i+lq;
                  else il=l+iq;
                  x=a[il]*cosx-a[im]*sinx;
                  a[im]=a[il]*sinx+a[im]*cosx;
                  a[il]=x;
                  }
                if (mv)
                  { ilr=ilq+i;
                  imr=imq+i;
                  x=r[ilr]*cosx-r[imr]*sinx;
                  r[imr]=r[ilr]*sinx+r[imr]*cosx;
                  r[ilr]=x;
                  }
                }

              x=2.*a[lm]*sincs;
              y=a[ll]*cosx2+a[mm]*sinx2-x;
              x=a[ll]*sinx2+a[mm]*cosx2+x;
              a[lm]=(a[ll]-a[mm])*sincs+a[lm]*(cosx2-sinx2);
              a[ll]=y;
              a[mm]=x;
              }

            /* Tests for completion.
               Test for m = last column. */
            } while (++m!=n);

          /* Test for l =penultimate column. */
          } while (++l!=n-1);
        } while (ind);
      
      /* Compare threshold with final norm. */
      } while (thr>anrmx);
    }
}

/***********************************************************************

eigen_sort  -  Sort eigenvalues & eigenvectors (if mv=1) in order of
               descending eigenvalue.
               Arguments are as for eigen, which must have been called.

***********************************************************************/

void eigen_sort(int mv, int n, double* a, double* r)
{
  int i, im, j, k, km, l, m;
  double am;

  k=0;
  for (i=0; i<n-1; i++)
    { im=i;
    km=k;
    am=a[k];
    l=0;
    for (j=0; j<n; j++)
      { if (j>i && a[l]>am)
        { im=j;
        km=l;
        am=a[l];
        }
      l+=j+2;
      }
    if (im!=i)
      { a[km]=a[k];
      a[k]=am;
      if (mv)
        { l=n*i;
        m=n*im;
        for (j=0; j<n; j++)
          { am=r[l];
          r[l++]=r[m];
          r[m++]=am;
          }
        }
      }
    k+=i+2;
    }
}
/***********************************************************************

vecprd  -  Calculates the vector product (vx,vy,vz) of two vectors 
           (px,py,pz) and (qx,qy,qz)

***********************************************************************/

void vecprd(double px, double py, double pz, double qx, double qy, double qz,
            double *vx, double *vy, double *vz)
{
  *vx = py * qz - pz * qy;
  *vy = pz * qx - px * qz;
  *vz = px * qy - py * qx;
}
/***********************************************************************

dot3  -  Returns the scalar product of the two vectors (x1,y1,z1) 
         and (x2,y2,z2)

***********************************************************************/

double dot3(double x1, double y1, double z1, double x2, double y2, double z2)
{
  double dot_prod;

  dot_prod = x1 * x2 + y1 * y2 + z1 * z2;
  return(dot_prod);
}
/***********************************************************************

orien3  -  Returns the orientation of the polygon with the consecutive
           vertices (x1,y1,z1), (x2,y2,z2), and (x3,y3,z3) as viewed
           from (ex,ey,ez).
                      -1  :  clockwise orientation
                      +1  :  anti-clockwise orientation (desired order)
                       0  :  degenerate (line or point)

 [Adapted from Angell I O & Griffith G (1987). High-resolution Computer
  Graphics using FORTRAN 77. Macmillan Education Ltd, Basingstoke, UK.]

***********************************************************************/

int orien3(double x1, double y1, double z1, double x2, double y2, double z2,
           double x3, double y3, double z3, double ex, double ey, double ez)
{
  double a, b, c, dx1, dy1, dz1, dx2, dy2, dz2, fe;
  int orientation;

  /* Initialise variables */
  dx1 = x2 - x1;
  dy1 = y2 - y1;
  dz1 = z2 - z1;
  dx2 = x3 - x2;
  dy2 = y3 - y2;
  dz2 = z3 - z2;

  /* Calculate triple scalar product */
  vecprd(dx1,dy1,dz1,dx2,dy2,dz2,&a,&b,&c);
  fe = dot3(a,b,c,ex-x1,ey-y1,ez-z1);
  if (fe < 0.0)
    orientation = -1;
  else
    orientation = 1;
  return(orientation);
}
/***********************************************************************

principal_components  -  Calculate the transformation to align the
                         supplied atom coords with their principal axes
                         along the x-, y-, and z-directions

***********************************************************************/

void principal_components(struct atom *first_atom_ptr,int natoms,
                          double transformation[3][3],
                          double eigen_value[3])

{
  int iatom, nzeroy, nzeroz;
  int array_size, i, iorder, ipos, j, npoints, nsize;

  double amatrix[6], matrix[3][3], rmatrix[9], x, y, z;
  double xxsum, xysum, xzsum, yysum, yzsum, zzsum;

  struct atom *atom_ptr;

  /* Initialise counts */
  iatom = 0;
  npoints = 0;
  nzeroy = 0;
  nzeroz = 0;
  xxsum = 0.0;
  xysum = 0.0;
  xzsum = 0.0;
  yysum = 0.0;
  yzsum = 0.0;
  zzsum = 0.0;

  /* Initialise transformation matrix */
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      if (i == j)
        transformation[i][j] = 1.0;
      else
        transformation[i][j] = 0.0;

  /* Initialise the eigen values */
  for (i = 0; i < 3; i++)
    eigen_value[i] = 1.0;

  /* If only have a single atom, then return */
  if (natoms < 2)
    return;

  /* Get pointer to the first atom */
  atom_ptr = first_atom_ptr;

  /* Loop through all vertices */
  while (atom_ptr != NULL && iatom < natoms)
    {
      /* Retrieve the coordinates for this atom */
      x = atom_ptr->view_coords[0];
      y = atom_ptr->view_coords[1];
      z = atom_ptr->view_coords[2];

      /* Count instances where the y- and z-coords are zero */
      if (fabs(y) < 0.000001)
        nzeroy++;
      if (fabs(z) < 0.000001)
        nzeroz++;

      /* Compute the terms of the scalar products matrix */
      xxsum = xxsum+x*x;
      xysum = xysum+x*y;
      xzsum = xzsum+x*z;
      yysum = yysum+y*y;
      yzsum = yzsum+y*z;
      zzsum = zzsum+z*z;

      /* Increment count of vertices used */
      npoints++;

      /* Get pointer to the next atom */
      atom_ptr = atom_ptr->next_atom_ptr;
      iatom++;
    }
    
  /* If only have a single point, then return */
  if (npoints < 2)
    return;

  /* Form the scalar products matrix */
  matrix[0][0] = xxsum;
  matrix[0][1] = xysum;
  matrix[0][2] = xzsum;
  matrix[1][0] = xysum;
  matrix[1][1] = yysum;
  matrix[1][2] = yzsum;
  matrix[2][0] = xzsum;
  matrix[2][1] = yzsum;
  matrix[2][2] = zzsum;

  /* If there are only two points, and both their z-coords are zero, or
     both their y- and z-coords are zero, adjust the order of the scalar
     products matrix whose eigenvalues are required */
  nsize = 3;
  if ((npoints == 2) && (nzeroz == npoints))
    {
      if (nzeroy == npoints)
        nsize = 1;
      else
        nsize = 2;
    }

  /* Calculate the eigen vectors and eigen values providing order of the
     matrix is not 1 */
  if (nsize > 1)
    {
      /* Convert matrix into form required by the eigen routine */
      ipos = 0;
      for (i = 0; i < nsize; i++)
        for (j = 0; j <= i; j++)
          {
            amatrix[ipos] = matrix[i][j];
            ipos++;
          }

      /* Calculate the eigen values */
      array_size = nsize * (nsize + 1) / 2;
      calc_eigen_values(1,nsize,amatrix,rmatrix);

      /* Sort the eigen vectors according to the eigen values */
      eigen_sort(1,nsize,amatrix,rmatrix);

      /* Retrieve the eigen values */
      i = j = 0;
      for (ipos = 0; ipos < array_size; ipos++)
        {
          /* If this is a diagonal term, then is one of the eigen values */
          if (i == j)
            eigen_value[i] = amatrix[ipos];

          /* Increment array indices */
          j++;
          if (j > i)
            {
              i++;
              j = 0;
            }
        }

      /* Retrieve the eigen vectors and hence transformation matrix */
      i = j = 0;
      for (ipos = 0; ipos < nsize * nsize; ipos++)
        {
          /* Store the transformation matrix */
          transformation[j][i] = rmatrix[ipos];

          /* Increment array indices */
          j++;
          if (j > nsize - 1)
            {
              i++;
              j = 0;
            }
        }

      /* If the matrix was of order 2, fix up row 3 and column 3 of the
         output transformation matrix */
      if (nsize == 2)
        {
          for (i = 0; i < 3; i++)
            {
              for (j = 0; j < 3; j++)
                {
                  if ((i == 2) && (j == 2))
                    transformation[i][j] = 1.0;
                  else if ((i == 2) || (j == 2))
                    transformation[i][j] = 0.0;
                }
            }
        }

      /* If the matrix was of order 2, fix up row 3 and column 3 of the
         output transformation matrix */
      if (nsize == 2)
        {
          for (i = 0; i < 3; i++)
            {
              for (j = 0; j < 3; j++)
                {
                  if ((i == 3) && (j == 3))
                    transformation[i][j] = 1.0;
                  else if ((i == 3) || (j == 3))
                    transformation[i][j] = 0.0;
                }
            }
        }
    }

  /* If the matrix was of order 1, make the transformation matrix the
     identity matrix */
  else
    {
      for (i = 0; i < 3; i++)
        {
          for (j = 0; j < 3; j++)
            {
              if (i == j)
                transformation[i][j] = 1.0;
              else
                transformation[i][j] = 0.0;
            }
        }
    }
    
  /* Check whether the sorted eigenvectors make a right-handed set */
  iorder = orien3(transformation[0][0],transformation[0][1],
                  transformation[0][2],0.0,0.0,0.0,
                  transformation[1][0],transformation[1][1],
                  transformation[1][2],transformation[2][0],
                  transformation[2][1],transformation[2][2]);
  
  /* If the sorted vectors don't make a right-handed set, then reverse the
     signs on the third vector */
  if (iorder == 1)
    {
      transformation[0][2] = - transformation[0][2];
      transformation[1][2] = - transformation[1][2];
      transformation[2][2] = - transformation[2][2];
    }
}
/***********************************************************************

rotate_molecule  -  Apply the given rotation to the molecule

***********************************************************************/

void rotate_molecule(struct atom *first_atom_ptr,int natoms,
                     double transformation[3][3],double *coord_min,
                     double *coord_max)
{
  int iatom, icoord;

  double x, y, z;

  struct atom *atom_ptr;

  /* Get pointer to the first atom */
  atom_ptr = first_atom_ptr;
  iatom = 0;

  /* Loop through all the atoms */
  while (atom_ptr != NULL && iatom < natoms)
    {
      /* Retrieve the coordinates for this atom */
      x = atom_ptr->view_coords[0];
      y = atom_ptr->view_coords[1];
      z = atom_ptr->view_coords[2];

      /* Apply the transformation */
      for (icoord = 0; icoord < 3; icoord++)
        {
          /* Calculate transformed coordinates */
          atom_ptr->view_coords[icoord]
            = ((x * transformation[0][icoord])
               + (y * transformation[1][icoord])
               + (z * transformation[2][icoord]));
        }

      /* If first atom, store new coords as the coord limits */
      if (atom_ptr == first_atom_ptr)
        {
          for (icoord = 0; icoord < 3; icoord++)
            {
              coord_min[icoord] = atom_ptr->view_coords[icoord];
              coord_max[icoord] = atom_ptr->view_coords[icoord];
            }
        }

      /* Otherwise, update limits */
      else
        {
          for (icoord = 0; icoord < 3; icoord++)
            {
              if (atom_ptr->view_coords[icoord] < coord_min[icoord])
                coord_min[icoord] = atom_ptr->view_coords[icoord];
              if (atom_ptr->view_coords[icoord] > coord_max[icoord])
                coord_max[icoord] = atom_ptr->view_coords[icoord];
            }
        }

      /* Get pointer to the next atom */
      atom_ptr = atom_ptr->next_atom_ptr;
      iatom++;
    }
}
/***********************************************************************

align_molecule  -  Use principal components to align molecule's
                   principal axes with x-, y- and z-axes

***********************************************************************/

void align_molecule(struct atom *first_atom_ptr,int natoms,
                    double *coord_min,double *coord_max,
                    double transformation[3][3],double cofg[3])
{
  double mass_x, mass_y, mass_z, eigen_value[3];

  /* Calculate the centre of mass of molecule's atoms */
  centre_of_mass(first_atom_ptr,natoms,&mass_x,&mass_y,&mass_z);

  /* Place the centre of mass at the origin */
  adjust_to_origin(first_atom_ptr,natoms,mass_x,mass_y,mass_z,FALSE);

  /* Store the centre of mass */
  cofg[0] = mass_x;
  cofg[1] = mass_y;
  cofg[2] = mass_z;

  /* Calculate transformation to orient the principal axes along
     the x-, y-, z-directions */
  principal_components(first_atom_ptr,natoms,transformation,eigen_value);

  /* Apply this transformation to the entire molecule */
  rotate_molecule(first_atom_ptr,natoms,transformation,coord_min,
                  coord_max);
}
/***********************************************************************

get_atom_colour  -  Get the given atom's colour

***********************************************************************/

void get_atom_colour(int type,double *rgb)
{
  char colour_name[COLNAM_LEN];

  int colour;

  /* Get the colour corresponding to the current atom type */
  if (type < NATOM_TYPES)
    strcpy(colour_name,ATOM_COLOUR[type]);

  /* If not a standard atom, must be a metal, so get metal colour */
  else
    {
      /* Get the metal type */
      type = type - NATOM_TYPES;

      /* Get the corresponding colour */
      colour = type % NMETAL_COLOURS;
      strcpy(colour_name,METAL_COLOUR[colour]);
    }

  /* Get rgb value for this colour */
  get_colour_rgb(colour_name,rgb);
}
/***********************************************************************

write_raster3d_atom  -  Write out the given atom to the Raster3D
                        output file

***********************************************************************/

double write_raster3d_atom(FILE *file_ptr,struct atom *atom_ptr,
                           double scale_factor,double percen,
                           double min_z)
{
  double radius, rgb[3], x, y, z;

  /* Transform the current atom's coordinates into the Raster3D
     reference frame */
  x = scale_factor * atom_ptr->view_coords[0];
  y = scale_factor * atom_ptr->view_coords[1];
  z = scale_factor * atom_ptr->view_coords[2];

  /* Get the atom colour */
  get_atom_colour(atom_ptr->type,rgb);

  /* Get the atom's radius */
  if (atom_ptr->type < NATOM_TYPES)
    radius = VDW_RADIUS[atom_ptr->type];
  else
    radius = METAL_RADIUS;

  /* Get the sphere radius */
  radius = scale_factor * radius * percen;

  /* Adjust minimum z-vaolue, if necessary */
  if (z - radius < min_z)
    min_z = z - radius;

  /* Write out a sphere object */
  fprintf(file_ptr,"2\n");
  fprintf(file_ptr,"%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n",
          x,y,z,radius,rgb[0],rgb[1],rgb[2]);

  /* Return minimum z-value */
  return(min_z);
}
/***********************************************************************

write_raster3d_bond  -  Write out the given bond to the Raster3D
                        output file

***********************************************************************/

void write_raster3d_bond(FILE *file_ptr,struct atom *atom1_ptr,
                         struct atom *atom2_ptr,double scale_factor,
                         double radius,int split_bond,double *rgb,
                         int flat_ends,int grey_carbons)
{
  double x1, y1, z1, x2, y2, z2, xm, ym, zm;

  /* Transform the two atoms' coordinates into the Raster3D
     reference frame */
  x1 = scale_factor * atom1_ptr->view_coords[0];
  y1 = scale_factor * atom1_ptr->view_coords[1];
  z1 = scale_factor * atom1_ptr->view_coords[2];
  x2 = scale_factor * atom2_ptr->view_coords[0];
  y2 = scale_factor * atom2_ptr->view_coords[1];
  z2 = scale_factor * atom2_ptr->view_coords[2];

  /* Get the bond radius */
  radius = scale_factor * radius;

  /* If showing split-bonds, calculate the midpoint */
  if (split_bond == TRUE)
    {
      xm = (x1 + x2) / 2;
      ym = (y1 + y2) / 2;
      zm = (z1 + z2) / 2;

      /* Get the first atom colour */
      get_atom_colour(atom1_ptr->type,rgb);
      if (grey_carbons == TRUE && atom1_ptr->type < NATOM_TYPES &&
          !strcmp(ELEMENT[atom1_ptr->type]," C"))
        {
          rgb[0] = LIGHT_GREY_RED;
          rgb[1] = LIGHT_GREY_GREEN;
          rgb[2] = LIGHT_GREY_BLUE;
        }

      /* Write out a flat-ended cylinder representing the first
         half-bond */
      fprintf(file_ptr,"5\n");
      fprintf(file_ptr,"%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n",
              x1,y1,z1,radius,xm,ym,zm);
      fprintf(file_ptr,"0.0 %6.3f %6.3f %6.3f\n",rgb[0],rgb[1],rgb[2]);

      /* Get the second atom colour */
      get_atom_colour(atom2_ptr->type,rgb);

      /* Write out a flat-ended cylinder representing the second
         half-bond */
      if (flat_ends == TRUE)
        fprintf(file_ptr,"5\n");
      else
        fprintf(file_ptr,"3\n");
      fprintf(file_ptr,"%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n",
              x2,y2,z2,radius,xm,ym,zm);
      fprintf(file_ptr,"0.0 %6.3f %6.3f %6.3f\n",rgb[0],rgb[1],rgb[2]);
    }

  /* Otherwise, use the given colour and write out round-ended
     cylinders for the bonds */
  else
    {
      /* Write out a round-ended cylinder representing the whole
         bond */
      fprintf(file_ptr,"3\n");
      fprintf(file_ptr,"%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n",
              x1,y1,z1,radius,x2,y2,z2);
      fprintf(file_ptr,"0.0 %6.3f %6.3f %6.3f\n",rgb[0],rgb[1],rgb[2]);
    }
}
/***********************************************************************

write_raster3d_residue  -  Write out the given residue to the Raster3D
                           output file

***********************************************************************/

double write_raster3d_residue(FILE *file_ptr,struct residue *residue_ptr,
                              double scale_factor)
{
  int atom_number, iatom, natoms;

  double min_z;
  static double rgb[3] = { 0.0, 0.0, 0.0 };

  struct atom *atom_ptr;
  struct conect *conect_ptr;
  struct residue *other_residue_ptr;

  /* Initialise variables */
  min_z = 0.0;

  /* Get pointer to the first atom record */
  atom_ptr = residue_ptr->first_atom_ptr;
  natoms = residue_ptr->natoms;
  iatom = 0;

  /* Loop through all the residue's atoms */
  while (atom_ptr != NULL && iatom < natoms)
    {
      /* Get the current atom number */
      atom_number = atom_ptr->atom_number;

      /* Write out atom as a sphere in Raster3D format */
      min_z = write_raster3d_atom(file_ptr,atom_ptr,scale_factor,
                                  LIGAND_ATOM_PERCEN,min_z);

      /* Get the first of this atom's covalent bonds */
      conect_ptr = atom_ptr->first_conect_ptr;

      /* Loop over this atom's covalent bonds and write out */
      while (conect_ptr != NULL)
        {
          /* If this connects to a higher-numbered atom, then
             want this bond */
          if (conect_ptr->to_atom_ptr->atom_number > atom_number)
            {
              /* Get this atom's residue */
              other_residue_ptr = conect_ptr->to_atom_ptr->residue_ptr;

              /* Check that the atom belongs to the same residue and is not
                 from a covalent bond to protein/DNA/metal */
              if (other_residue_ptr == residue_ptr)
                {
                  /* Write out this bond */
                  write_raster3d_bond(file_ptr,atom_ptr,
                                      conect_ptr->to_atom_ptr,scale_factor,
                                      LIGAND_BOND_RADIUS,TRUE,rgb,TRUE,
                                      FALSE);
                }
            }

          /* Get the next conect */
          conect_ptr = conect_ptr->next_conect_ptr;
        }
  
      /* Get pointer to next atom in linked-list */
      atom_ptr = atom_ptr->next_atom_ptr;
      iatom++;
    }

  /* Return minimum z-value for placing of back wall */
  return(min_z);
}
/***********************************************************************

check_conect  -  Check whether we already have a connection between
                 the given two atoms

***********************************************************************/

int check_conect(struct atom *atom1_ptr,struct atom *atom2_ptr)
{
  int have_bond;

  struct conect *conect_ptr;

  /* Initialise variables */
  have_bond = FALSE;

  /* Get pointer to the first atom's first conect record */
  conect_ptr = atom1_ptr->first_conect_ptr;

  /* Loop through all the conect records to see if we have one
     linking to the second atom*/
  while (conect_ptr != NULL && have_bond == FALSE)
    {
      /* If this connects to the second atom, then have a bond */
      if (conect_ptr->to_atom_ptr == atom2_ptr)
        have_bond = TRUE;

      /* Get the next conect */
      conect_ptr = conect_ptr->next_conect_ptr;
    }

  /* Return whether bond for these two atoms already exists */
  return(have_bond);
}
/***********************************************************************

calc_covalent_bonds  -  Calculate all the covalent bonds in the given
                      ligand

***********************************************************************/

void calc_covalent_bonds(struct residue *first_residue_ptr,int nresid)
{
  int i, iatom, jatom, iresid, jresid, niatoms, njatoms;
  int have_bond, toofar;

  double dist_sqrd;

  struct atom *atom_ptr, *other_atom_ptr;
  struct residue *residue_ptr, *other_residue_ptr;

  /* Get the first residue */
  residue_ptr = first_residue_ptr;
  iresid = 0;

  /* Loop through the residues, considering all within-range residues
     which may contain covalent bonds */
  while (residue_ptr != NULL && iresid < nresid)
    {
      /* Initially have other residue same as current residue */
      other_residue_ptr = residue_ptr;
      jresid = iresid;

      /* Loop through all remaining residues to find in-range residue */
      while (other_residue_ptr != NULL && jresid < nresid)
        {
          /* Initialise flag */
          toofar = FALSE;

          /* Determine whether these two residues are within range */
          if (other_residue_ptr != residue_ptr)
            {
              /* Check whether these residue are within reach */
              for (i = 0; i < 3 && toofar == FALSE; i++)
                {
                  if (residue_ptr->coord_min[i]
                      - other_residue_ptr->coord_max[i] > TOO_FAR_DIST)
                    toofar = TRUE;
                  if (other_residue_ptr->coord_min[i]
                      - residue_ptr->coord_max[i] > TOO_FAR_DIST)
                    toofar = TRUE;
                }
            }
      
          /* If residues within range calculate all atom-atom distances */
          if (toofar == FALSE)
            {
              /* Get the first of the first residue's atoms */
              atom_ptr = residue_ptr->first_atom_ptr;
              iatom = 0;
              niatoms = residue_ptr->natoms;
  
              /* Loop through all the residue's atoms */
              while (atom_ptr != NULL && iatom < niatoms)
                {
                  /* Initialise flag */
                  toofar = FALSE;

                  /* Determine whether these two residues are within range */
                  if (other_residue_ptr != residue_ptr)
                    {
                      /* Check whether this atom within range of the other
                         residue */
                      for (i = 0; i < 3 && toofar == FALSE; i++)
                        {
                          if (atom_ptr->coords[i]
                              - other_residue_ptr->coord_max[i]
                              > TOO_FAR_DIST ||
                              other_residue_ptr->coord_min[i]
                              - atom_ptr->coords[i] > TOO_FAR_DIST)
                            toofar = TRUE;
                        }
                    }
                  /* If atom not too far from the other residue, check its
                     distance from all that residue's atoms */
                  if (toofar == FALSE)
                    {
                      /* If looking same residue, then get next atom */
                      if (other_residue_ptr == residue_ptr)
                        {
                          /* Get this residue's next atom */
                          other_atom_ptr = atom_ptr->next_atom_ptr;
                          jatom = iatom + 1;
                          njatoms = niatoms;
                        }

                      /* Otherwise, et the first of the other residue's
                         atoms */
                      else
                        {
                          other_atom_ptr
                            = other_residue_ptr->first_atom_ptr;
                          jatom = 0;
                          njatoms = other_residue_ptr->natoms;
                        }
  
                      /* Loop through all the second residue's atoms */
                      while (jatom < njatoms && other_atom_ptr != NULL)
                        {
                          /* Calculate the squared distance between these
                             two atoms */
                          dist_sqrd = 0.0;
                          for (i = 0; i < 3; i++)
                            dist_sqrd
                              = dist_sqrd
                              + (atom_ptr->coords[i]
                                 - other_atom_ptr->coords[i])
                              * (atom_ptr->coords[i]
                                 - other_atom_ptr->coords[i]);

                          /* If within covalent bonding distance, have a
                             bond between the two residues */
                          if (dist_sqrd < COVALENT_DIST_SQRD)
                            {
                              /* Check whether this bond is already
                                 covered by a previous CONECT record */
                              have_bond
                                = check_conect(atom_ptr,other_atom_ptr);

                              /* If don't already have this bond, create it */
                              if (have_bond == FALSE)
                                {
                                  /* Create a connection from the first
                                     atom to the second */
                                  create_conect_record(atom_ptr,
                                                       other_atom_ptr,FALSE);

                                  /* Create a connection the other way round */
                                  create_conect_record(other_atom_ptr,
                                                       atom_ptr,FALSE);
                                }
                            }

                          /* Get pointer to next atom in linked-list */
                          other_atom_ptr = other_atom_ptr->next_atom_ptr;
                          jatom++;
                        }
                    }

                  /* Get pointer to next atom in linked-list */
                  atom_ptr = atom_ptr->next_atom_ptr;
                  iatom++;
                }
            }

          /* Get the next residue */
          other_residue_ptr = other_residue_ptr->next_residue_ptr;
          jresid++;
        }

      /* Get the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;
    }
}
/***********************************************************************

write_raster3d_headers  -  Write out the header records to the output
                           Raster3D file

***********************************************************************/

void write_raster3d_headers(FILE *file_ptr)
{
  /* Initialise variables */
  fprintf(file_ptr,"Generated by newsec. R A Laskowski (2004).\n");
  fprintf(file_ptr,"15 15     tiles in x,y\n");
  fprintf(file_ptr,"24 24     pixels (x,y) per tile\n");
  fprintf(file_ptr,"3         3x3 virtual pixels -> 2x2 pixels\n");
  fprintf(file_ptr,"1.0 1.0 1.0   Background colour\n");
  fprintf(file_ptr,"T         shadows cast\n");
  fprintf(file_ptr,"25        Phong power\n");
  fprintf(file_ptr,"0.15      secondary light contribution\n");
  fprintf(file_ptr,"0.05      ambient light contribution\n");
  fprintf(file_ptr,"0.25      specular reflection component\n");
  fprintf(file_ptr,"40.0      eye position\n");
  fprintf(file_ptr,"0.5 0 1   main light source position\n");
  fprintf(file_ptr,"1 0 0 0   input coordinate, radius ");
  fprintf(file_ptr,"transformation\n");
  fprintf(file_ptr,"0 1 0 0\n");
  fprintf(file_ptr,"0 0 1 0\n");
  fprintf(file_ptr,"0 0 0 1\n");
  fprintf(file_ptr,"3         mixed objects\n");
  fprintf(file_ptr,"*\n");
  fprintf(file_ptr,"*\n");
  fprintf(file_ptr,"*\n");
}
/***********************************************************************

calc_r3d_scale  -  Calculate Raster3D scaling factor

***********************************************************************/

double calc_r3d_scale(double *coord_min,double *coord_max)
{
  double scale_factor;
  double max_width, tmp, x_width, y_width;

  /* Initialise variables */
  scale_factor = 1.0;

  /* Calculate the extent of the molecule in each direction */
  x_width = 2 * fabs(coord_max[0]);
  tmp = 2 * fabs(coord_min[0]);
  if (tmp > x_width)
    x_width = tmp;
  y_width = 2 * fabs(coord_max[1]);
  tmp = 2 * fabs(coord_min[1]);
  if (tmp > y_width)
    y_width = tmp;

  /* Compute the limiting dimension */
  max_width = x_width;
  if (y_width > max_width)
    max_width = y_width;

  /* Calculate the Raster3D scaling factor */
  if (max_width < 0.001)
    max_width = 1.0;
  scale_factor = (RASTER3D_MAX - RASTER3D_MIN) / max_width;

  /* Return the scaling factor */
  return(scale_factor);
}
/***********************************************************************

run_raster3d  -  Run Raster3d for current image and convert to GIF

***********************************************************************/

void run_raster3d(char *r3d_name,char *gif_name,char *render_exe,
                  char *convert_exe,char *ppmquant_exe,
                  char *ppmtogif_exe,char *image,int crop,int tiny)
{
  char command[FILENAME_LEN + 1], image_size[15], quality[15];
  char out_name[FILENAME_LEN];

  /* Form the parameters to run the render program */
  printf("Running Raster3D on %s ...\n",image);

  /* If image to be cropped, want to generate it at 100% quality */
  if (crop == TRUE)
    {
      strcpy(quality,"-quality 100");
      strcpy(out_name,"lx.jpg");
    }
  else
    {
      quality[0] = '\0';
      strcpy(out_name,gif_name);
    }

  /* If tiny image required (for MCSG gallery) define size */
  if (tiny == TRUE)
    strcpy(image_size," ");
  else
    image_size[0] = '\0';

  /* If using JPEG files, create directly */
  if (IMAGE_FORMAT == JPEG)
    {
/* debug */
//      sprintf(command,"%s -jpeg %s %s %s < %s >& /dev/null",render_exe,
//              out_name,quality,image_size,r3d_name);
      sprintf(command,"%s -jpeg %s %s %s < %s >& concert_a.log",
              render_exe,out_name,quality,image_size,r3d_name);
/* debug */
      system(command);

      /* If cropping this image, use convert to do so */
      if (crop == TRUE)
        {
          sprintf(command,"%s -crop 0x0 %s %s",convert_exe,
                  out_name,gif_name);
          system(command);

          /* Tidy up by removing the intermediate files */
          sprintf(command,"%s -f lx.jpg",RM_COMMAND);
          system(command);
        }
    }

  /* Otherwise, generate GIF file via TIFF */
  else
    {
      /* Run the Raster3D render program */
      sprintf(command,"%s -tiff s.tiff < %s >& /dev/null",render_exe,
              r3d_name);
      system(command);

      /* Form the parameters to convert the TIFF file just created to PPM
         format */
      printf("   Running convert on hit ...\n");
      sprintf(command,"%s tiff:s.tiff ppm:s.pbm",convert_exe);

      /* Run convert to generate intermediate PBM file */
      system(command);

      /* Limit number of colours to 256 so that ppmtogif will work */
      sprintf(command,"%s 256 s.pbm > s.ppm",ppmquant_exe);
      system(command);

      /* Finally, use ppmtogif to create the required GIF file */
      sprintf(command,"%s s.ppm > %s",ppmtogif_exe,gif_name);
      system(command);

      /* Tidy up by removing the intermediate files */
      sprintf(command,"%s %s s.jpg s.pbm s.ppm",RM_COMMAND,r3d_name);
      system(command);
    }
}
/***********************************************************************

render_hetgroup  -  Render the given Het Group

***********************************************************************/

void render_hetgroup(char *res_name,struct residue *residue_ptr,
                     char *render_exe,char *convert_exe,
                     char *ppmquant_exe,char *ppmtogif_exe,
                     char *output_dir,int molecule_align,
                     struct parameters params)
{
  char gif_name[FILENAME_LEN + 1], file_r3d[FILENAME_LEN + 1];
  char image[20];

  int i;

  double coord_min[3], coord_max[3], min_z, scale_factor;
  double transformation[3][3], cofg[3];

  FILE *file_ptr;

  /* Form file names */
  strcpy(file_r3d,"h.r3d");
  if (IMAGE_FORMAT == GIF)
    sprintf(gif_name,"het%s.gif",res_name);
  else
    sprintf(gif_name,"het%s.jpg",res_name);
  sprintf(image,"Het Group %s",res_name);

  /* Calculate all the covalent bonds in the residue */
  if (residue_ptr->got_bonds == FALSE)
    {
      calc_covalent_bonds(residue_ptr,1);
      residue_ptr->got_bonds = TRUE;
    }

  /* Open the l.r3d output file */
  file_ptr = open_output_file(file_r3d,TRUE);

  /* Write out the Raster3D header records */
  write_raster3d_headers(file_ptr);

  /* Use principal components to align residue along its
     principal axes */
  if (molecule_align == TRUE)
    align_molecule(residue_ptr->first_atom_ptr,residue_ptr->natoms,
                   coord_min,coord_max,transformation,cofg);
  
  /* Adjust max and min coords for largest atom radius */
  for (i = 0; i < 3; i++)
    {
      coord_min[i] = coord_min[i] - 2.0;
      coord_max[i] = coord_max[i] + 2.0;
    }          

  /* Calculate Raster3D scaling factor */
  scale_factor = calc_r3d_scale(coord_min,coord_max);

  /* Write out Het Group to Raster3D file */
  min_z = write_raster3d_residue(file_ptr,residue_ptr,scale_factor);

  /* Place the back wall just behind the molecule */
  min_z = min_z - scale_factor * BACK_WALL_GAP;
  // if (residue_ptr->metal == FALSE)
  //  {
  // fprintf(file_ptr,"6\n");
  // fprintf(file_ptr,"0.0 0.0 %6.3f 0.2 0.0 %6.3f 0.2 0.2 %6.3f\n",
  //         min_z,min_z,min_z);
  // fprintf(file_ptr,"1.0 1.0 1.0\n");
  //  }
  
  /* Close the l.r3d output file */
  fclose(file_ptr);

  /* Run Raster3D to render the image */
  if (params.run_raster3d == TRUE)
    run_raster3d(file_r3d,gif_name,render_exe,convert_exe,
                 ppmquant_exe,ppmtogif_exe,image,FALSE,FALSE);
}
/***********************************************************************

write_coords  -  Write the atom coords in PDB format

***********************************************************************/

void write_coords(char *file_name,struct residue *first_residue_ptr,
                  int nresid,int plus_conect)
{
  int iatom, iresid, natoms, ntotal_atoms;

  struct atom *atom_ptr, *first_atom_ptr;
  struct conect *conect_ptr;
  struct residue *residue_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  ntotal_atoms = 0;

  /* Open the output PDB file */
  file_ptr = open_output_file(file_name,TRUE);

  /* Get pointer to the first residue record */
  residue_ptr = first_residue_ptr;
  iresid = 0;
  if (residue_ptr != NULL)
    first_atom_ptr = residue_ptr->first_atom_ptr;
  else
    first_atom_ptr = NULL;

  /* Loop over the residues to write out */
  while (residue_ptr != NULL && iresid < nresid)
    {
      /* Get the first of this residue's atoms */
      atom_ptr = residue_ptr->first_atom_ptr;
      natoms = residue_ptr->natoms;
      iatom = 0;

      /* Loop over this residue's atoms to write out */
      while (atom_ptr != NULL && iatom < natoms)
        {
          /* Start the ATOM record */
          fprintf(file_ptr,"ATOM  ");

          /* Write out the details */
          fprintf(file_ptr,"%5d %s %s %c%s   %8.3f%8.3f%8.3f",
                  atom_ptr->atom_number,atom_ptr->atom_name,
                  residue_ptr->res_name,
                  residue_ptr->chain,residue_ptr->res_num,
                  atom_ptr->view_coords[0],atom_ptr->view_coords[1],
                  atom_ptr->view_coords[2]);
          fprintf(file_ptr,"%6.2f%6.2f\n",atom_ptr->occupancy,
                  atom_ptr->bvalue);

          /* Get pointer to next atom record and increment
             atom counter */
          atom_ptr = atom_ptr->next_atom_ptr;
          iatom++;
          ntotal_atoms++;
        }

      /* Get the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;
    }

  /* If need to write out the CONECT records, loop through all the
     atoms again */
  if (plus_conect == TRUE)
    {
      /* Write out a TER record */
      fprintf(file_ptr,"TER\n");

      /* Get the first atom */
      atom_ptr = first_atom_ptr;
      natoms = ntotal_atoms;
      iatom = 0;

      /* Loop over all the atoms to get any CONECT information */
      while (atom_ptr != NULL && iatom < natoms)
        {
          /* Get pointer to the first atom's first conect record */
          conect_ptr = atom_ptr->first_conect_ptr;

          /* Loop through all the conect records to find those that
             came from the PDB file */
          while (conect_ptr != NULL)
            {
              /* If this came from the original PDB file, write out */
              if (conect_ptr->from_conect == TRUE)
                {
                  /* Write out a CONECT record */
                  fprintf(file_ptr,"CONECT%5d%5d\n",atom_ptr->atom_number,
                          conect_ptr->to_atom_ptr->atom_number);
                }

              /* Get the next conect */
              conect_ptr = conect_ptr->next_conect_ptr;
            }


          /* Get pointer to next atom record and increment
             atom counter */
          atom_ptr = atom_ptr->next_atom_ptr;
          iatom++;
        }

    }
      
  /* Close the output file */
  fclose(file_ptr);
}
/***********************************************************************

write_run_scr_params  -  Write out all the ##run.scr lines to the .run
                         file

***********************************************************************/

void write_run_scr_params(char *file_name,
                          struct residue *first_residue_ptr,int nresid,
                          int metal,int colour,FILE *ligmet_file_ptr)
{
  char colour_name[20], het_name[4], run_scr[12];
  char chain, res_name1[4], res_name2[4], res_num1[6], res_num2[6];

  char *ligand_seq;

  int i, iresid, ntotal_atoms;

  struct residue *residue_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  ligand_seq = &null;
  ntotal_atoms = 0;
  chain = ' ';
  res_name1[0] = '\0';
  res_name2[0] = '\0';
  res_num1[0] = '\0';
  res_num2[0] = '\0';
  strcpy(run_scr,"  ##run.scr");

  /* Open the output .run file */
  if (ligmet_file_ptr == NULL)
    file_ptr = open_output_file(file_name,TRUE);
  else
    {
      /* Write to ligmet.run file */
      file_ptr = ligmet_file_ptr;

      /* Write out the file name */
      fprintf(file_ptr,"Name: %s\n",file_name);
      run_scr[0] = '\0';
    }

  /* Get pointer to the first residue record */
  residue_ptr = first_residue_ptr;
  iresid = 0;

  /* Loop over the residues to write out */
  while (residue_ptr != NULL && iresid < nresid)
    {
      /* If this is the first residue, store its name and number */
      if (iresid == 0)
        {
          /* Store name and number */
          strcpy(res_name1,residue_ptr->res_name);
          strcpy(res_num1,residue_ptr->res_num);

          /* Add to ligand_sequence string */
          store_string(&ligand_seq,residue_ptr->res_name);
        }

      /* Otherwise, append to ligand sequence string */
      else
        {
          append_string(&ligand_seq,"-");
          append_string(&ligand_seq,residue_ptr->res_name);
        }

      /* Store name and number as last residue so far */
      strcpy(res_name2,residue_ptr->res_name);
      strcpy(res_num2,residue_ptr->res_num);
      chain = residue_ptr->chain;

      /* Write out thi residue */
      fprintf(file_ptr,"#residue %s  %s%c\n",residue_ptr->res_name,
              residue_ptr->res_num,residue_ptr->chain);

      /* Increment atom count */
      ntotal_atoms = ntotal_atoms + residue_ptr->natoms;

      /* Get the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;
    }

  /* Write out the opening lines to the ##run.scr file */
  fprintf(file_ptr,"set from_dic = TRUE%s\n",run_scr);
  fprintf(file_ptr,"@ natoms = %d%s\n",ntotal_atoms,run_scr);

  /* First residue details */
  strcpy(het_name,res_name1);
  for (i = 0; i < 3; i++)
    if (het_name[i] == ' ')
      het_name[i] = '_';
  fprintf(file_ptr,"set ligand_seq = \"%s\"%s\n",ligand_seq,run_scr);
  fprintf(file_ptr,"set het_name = \"%s\"%s\n",het_name,run_scr);
  fprintf(file_ptr,"set resname1 = \"-n%s\"%s\n",res_name1,run_scr);
  fprintf(file_ptr,"set res1 = \"%s\"%s\n",res_num1,run_scr);

  /* Second residue details */
  fprintf(file_ptr,"set resname2 = \"-n%s\"%s\n",res_name2,run_scr);
  fprintf(file_ptr,"set res2 = \"%s\"%s\n",res_num2,run_scr);

  /* Write out the chain-id */
  fprintf(file_ptr,"set chain = \"%c\"%s\n",chain,run_scr);

  /* If metal, write out its colour */
  if (metal == TRUE)
    {
      /* Get the colour */
      strcpy(colour_name,METAL_COLOUR[colour]);

      /* Convert to upper case */
      to_upper(colour_name,strlen(colour_name));

      /* Write out */
      fprintf(file_ptr,"set colour = \"%s\"%s\n",colour_name,run_scr);
    }
      
  /* Close the output file */
  if (ligmet_file_ptr == NULL)
    fclose(file_ptr);
}
/***********************************************************************

write_hbadd_files  -  For each Het Group for which we have no
                      representative, or for which this is the best
                      representative, write out the .pdb and .rcm files
                      and generate a PyMOL or Raster3D image

***********************************************************************/

void write_hbadd_files(struct hbadd *first_hbadd_ptr,char *output_dir,
                       char *hetgif_dir,char *pdb_code,char *render_exe,
                       char *convert_exe,char *ppmquant_exe,
                       char *ppmtogif_exe,int molecule_align,
                       FILE *ligmet_file_ptr,struct parameters params)
{
  char file_name[FILENAME_LEN];
  char res_name[4];

  int i, het_count;
  int best;

  struct hbadd *hbadd_ptr;
  struct residue *residue_ptr;

  /* Initialise variables */
  het_count = 0;

  /* Get first Het Group */
  hbadd_ptr = first_hbadd_ptr;

  /* Loop through all the Het Groups, determining which are to be
     written out */
  while (hbadd_ptr != NULL)
    {
      /* Get the corresponding residue record */
      residue_ptr = hbadd_ptr->residue_ptr;

      /* Get Het Group name */
      strcpy(res_name,residue_ptr->res_name);

      /* Replace all blank spaces by underscores */
      for (i = 0; i < 3; i++)
        if (res_name[i] == ' ')
          res_name[i] = '_';

      /* Check whether we already have a representative of this Het
         Group in the hetgroups directory */
      best = check_het_group(res_name,hbadd_ptr->nmissed,hetgif_dir);

      /* If this is the best (or only) representative of this Het
         Group that we have, generate all the necessary files */
      if (best == TRUE && residue_ptr->type != METAL_ION)
        {
          /* Increment Het Group count */
          het_count++;

          /* Write out the .en file */
          write_en_file(res_name,pdb_code,residue_ptr,hbadd_ptr->nmissed,
                        output_dir);

          /* Write out the .r3d file */
          render_hetgroup(res_name,residue_ptr,render_exe,convert_exe,
                          ppmquant_exe,ppmtogif_exe,output_dir,
                          molecule_align,params);

          /* Form the file name for this Het Group */
          strcpy(file_name,output_dir);
          strcat(file_name,"het");
          strcat(file_name,res_name);
          strcat(file_name,".pdb");

          /* Write out the stored coordinates as a PDB file */
          write_coords(file_name,residue_ptr,1,TRUE);

          /* Repeat for the .rcm file */
          strcpy(file_name,output_dir);
          strcat(file_name,"het");
          strcat(file_name,res_name);
          strcat(file_name,".rcm");

          /* Write out the .rcm file */
          write_coords(file_name,residue_ptr,1,FALSE);

          /* Form name of .run file */
          sprintf(file_name,"%shet_%2.2d.run",output_dir,het_count);

          /* Write the .run file */
          write_run_scr_params(file_name,residue_ptr,1,FALSE,0,
                               ligmet_file_ptr);
        }

      /* Get the next Het Group */
      hbadd_ptr = hbadd_ptr->next_hbadd_ptr;
    }
}
/***********************************************************************

extract_sequences  -  Extract all the protein sequences

***********************************************************************/

void extract_sequences(struct chain *first_chain_ptr)
{
  char sequence[MAX_SEQLEN + 1];

  int ipos, iresid, nresid;

  struct chain *chain_ptr;
  struct residue *residue_ptr;

  /* Get pointer to the first chain */
  chain_ptr = first_chain_ptr;

  /* Loop over all the chains */
  while (chain_ptr != NULL)
    {
      /* Only process if this is a protein sequence */
      if (chain_ptr->type == PROTEIN || chain_ptr->type == CALPHA_ONLY)
        {
          /* Get the first residue */
          residue_ptr = chain_ptr->first_residue_ptr;
          nresid = chain_ptr->nresid;
          iresid = 0;

          /* Check whether sequence too long to store in its entirety */
          if (chain_ptr->nresid > MAX_SEQLEN)
            printf("* Warning. Chain [%c] sequence too long: %d\n",
                   chain_ptr->chain_id,chain_ptr->nresid);

          /* Initialise sequence */
          ipos = 0;
          sequence[ipos] = '\0';

          /* Loop through the chain's residues to extract the full
             sequence of amino acids */
          while (residue_ptr != NULL && iresid < nresid)
            {
              /* Store this residue's a.a. code */
              if (ipos < MAX_SEQLEN)
                {
                  sequence[ipos] = residue_ptr->aa_code;
                  ipos++;
                }

              /* Get the next residue */
              residue_ptr = residue_ptr->next_residue_ptr;
              iresid++;
            }

          /* Close off the sequence */
          sequence[ipos] = '\0';

          /* Store the sequence */
          store_string(&(chain_ptr->sequence),sequence);
        }

      /* Get the next chain */
      chain_ptr = chain_ptr->next_chain_ptr;
    }
}
/***********************************************************************

perform_comparison  -  Compare the two given sequences

***********************************************************************/

int perform_comparison(struct chain *chain_ptr,
                       struct chain *other_chain_ptr)
{
  char  block[TUPLEN + 1];

  char *string_ptr;

  int diff, ipos, match_count, minlen, nresid, other_nresid;
  int similar;

  float percen;

  /* Initialise variables */
  similar = FALSE;
  nresid = chain_ptr->nresid;
  other_nresid = other_chain_ptr->nresid;

  /* Get the length of the shortest chain */
  diff = (int) abs(nresid - other_nresid);
  minlen = nresid;
  if (other_nresid < nresid)
    minlen = other_nresid;
  percen = (float) (100.0 * (float) diff / (float) minlen);

  /* If sequence lengths are identical, check whether sequences are also
     identical */
  if (diff == 0)
    {
      /* Check for identity of sequences */
      if (!strcmp(other_chain_ptr->sequence,chain_ptr->sequence))
        similar = TRUE;
    }

  /* If don't have a perfect match, and chains are close in length, then
     perform quick comparison */
  if (similar == FALSE && percen < PERCEN_SEQLEN_DIFF && minlen > TUPLEN)
    {
      /* Initialise count of matching blocks */
      match_count = 0;

      /* Loop through the first sequence in blocks */
      for (ipos = 0; ipos < nresid - TUPLEN; ipos++)
        {
          /* Extract a sequence block from the first sequence */
          strncpy(block,other_chain_ptr->sequence + ipos,TUPLEN);
          block[TUPLEN] = '\0';

          /* Test whether this block is in the second sequence */
          string_ptr = strstr(chain_ptr->sequence,block);
          if (string_ptr != NULL)
            match_count++;
        }

      /* If have enough matches to consider the two sequences
         equal, flag a match */
      if (match_count >= (minlen - 2 * TUPLEN * minlen / 100))
        similar = TRUE;
    }

  /* Return whether the two sequences are similar */
  return(similar);
}
/***********************************************************************

compare_sequences  -  Compare the sequences to see whether they are
                      sufficiently similar to be considered duplicates

***********************************************************************/

int compare_sequences(struct chain *first_chain_ptr)
{
  char  add_on[4];

  int nunique;
  int match;

  struct chain *chain_ptr, *other_chain_ptr;

  /* Initialise variables */
  nunique = 0;

  /* Get pointer to first chain */
  chain_ptr = first_chain_ptr;

  /* Loop over all chains */
  while (chain_ptr != NULL)
    {
      /* Only process if this is a protein sequence */
      if (chain_ptr->type == PROTEIN || chain_ptr->type == CALPHA_ONLY)
        {
          /* Check that chain hasn't already been found to be a copy
             of another chain */
          if (chain_ptr->copy_of == '\0')
            {
              /* Get pointer to next chain */
              other_chain_ptr = chain_ptr->next_chain_ptr;

              /* Loop over all following chains */
              while (other_chain_ptr != NULL)
                {
                  /* Only process if this is a protein sequence */
                  if (other_chain_ptr->type == PROTEIN ||
                      other_chain_ptr->type == CALPHA_ONLY)
                    {
                      /* Perform the comparison */
                      match = perform_comparison(chain_ptr,other_chain_ptr);

                      /* If chains are sufficiently similar to count as copies
                         of the same chain, then mark second version as a copy
                         of the first */
                      if (match == TRUE)
                        {
                          other_chain_ptr->copy_of = chain_ptr->chain_id;
                          chain_ptr->ncopies++;

                          /* Update name of current chain */
                          if (!strncmp(chain_ptr->name,"Chain ",6))
                            sprintf(chain_ptr->name,"Chains %c",
                                    chain_ptr->chain_id);

                          /* Add the copy chain's id onto current chain's
                             name */
                          sprintf(add_on,", %c",other_chain_ptr->chain_id);
                          strcat(chain_ptr->name,add_on);

                          /* If chain lengths differ, set flag */
                          if (chain_ptr->nresid != other_chain_ptr->nresid)
                            chain_ptr->same_length = FALSE;
                        }
                    }

                  /* Get pointer to next chain */
                  other_chain_ptr = other_chain_ptr->next_chain_ptr;
                }

              /* If chain is a unique one, then increment count */
              if (chain_ptr->copy_of == '\0')
                nunique++;
            }
        }

      /* Get pointer to next chain */
      chain_ptr = chain_ptr->next_chain_ptr;
    }

  /* Write out the number of unique chains */
  printf("Total number of unique protein chains:  %8d\n",nunique);

  /* Return total number of unique chains */
  return(nunique);
}
/***********************************************************************

old_compare_sequences  -  Compare the sequences to see which are duplicates

***********************************************************************/

int old_compare_sequences(struct chain *first_chain_ptr)
{
  char add_on[4];

  int diff, iresid, jresid, niresid, njresid, minlen, nunique;
  int start_jresid;
  int nmatch, nmismatch;
  int found, match;

  float percen;

  struct chain *chain_ptr, *other_chain_ptr;
  struct residue *residue_ptr, *other_residue_ptr;
  struct residue *other_start_residue_ptr;

  /* Initialise variables */
  nunique = 0;

  /* Get pointer to first chain */
  chain_ptr = first_chain_ptr;

  /* Loop over all chains */
  while (chain_ptr != NULL)
    {
      /* Only process if this is a protein sequence */
      if (chain_ptr->type == PROTEIN || chain_ptr->type == CALPHA_ONLY)
        {
          /* Check that chain hasn't already been found to be a copy
             of another chain */
          if (chain_ptr->copy_of == '\0')
            {
              /* Get pointer to next chain */
              other_chain_ptr = chain_ptr->next_chain_ptr;

              /* Loop over all following chains */
              while (other_chain_ptr != NULL)
                {
                  /* Only process if this is a protein sequence */
                  if (other_chain_ptr->type == PROTEIN ||
                      other_chain_ptr->type == CALPHA_ONLY)
                    {
                      /* Initialise flag */
                      match = FALSE;

                      /* Get the length of the shortest chain */
                      diff = (int) abs(chain_ptr->nresid
                                       - other_chain_ptr->nresid);
                      minlen = chain_ptr->nresid;
                      if (other_chain_ptr->nresid < chain_ptr->nresid)
                        minlen = other_chain_ptr->nresid;
                      percen
                        = (float) (100.0 * (float) diff / (float) minlen);

                      /* Get the first of the first chain's residues */
                      residue_ptr = chain_ptr->first_residue_ptr;
                      niresid = chain_ptr->nresid;
                      iresid = 0;

                      /* Get the first of the second chain's residues */
                      other_start_residue_ptr
                        = other_chain_ptr->first_residue_ptr;
                      njresid = other_chain_ptr->nresid;
                      start_jresid = 0;
                      
                      /* Initialise count of matches and mismatches */
                      nmatch = 0;
                      nmismatch = 0;

                      /* Loop through all the residues to count the number
                         of residue number/name equivalences in the two
                         sequences */
                      while (residue_ptr != NULL && iresid < niresid)
                        {
                          /* Get the other chain's residues */
                          other_residue_ptr = other_start_residue_ptr;
                          jresid = start_jresid;

                          /* Initialise flag */
                          found = FALSE;

                          /* Loop through the other chain's residues until
                             hit a match on residue number */
                          while (other_residue_ptr != NULL &&
                                 jresid < njresid &&
                                 found == FALSE)
                            {
                              /* If residue numbers match, check whether
                                 residue names also match */
                              if (!strcmp(residue_ptr->res_num,
                                          other_residue_ptr->res_num))
                                {
                                  /* If names match, increment count of
                                     matches */
                                  if (!strcmp(residue_ptr->res_name,
                                              other_residue_ptr->res_name))
                                    nmatch++;

                                  /* Or mismatches */
                                  else
                                    nmismatch++;

                                  /* Move second sequence's start position
                                     along one */
                                  other_start_residue_ptr
                                    = other_residue_ptr->next_residue_ptr;
                                  start_jresid = jresid;
                                  
                                  /* Set flag to indicate we found what we
                                     were looking for */
                                  found = TRUE;
                                }
                              
                              /* Get the next residue */
                              other_residue_ptr
                                = other_residue_ptr->next_residue_ptr;
                              jresid++;
                            }
                          
                          /* Get the next residue */
                          residue_ptr = residue_ptr->next_residue_ptr;
                          iresid++;
                        }

                      /* Determine whether chains are similar */
                      if (nmismatch < 5 && nmatch > 0.8 * (float) minlen &&
                          percen < 20.0)
                        match = TRUE;

                      /* If chains are sufficiently similar to count as copies
                         of the same chain, then mark second version as a copy
                         of the first */
                      if (match == TRUE)
                        {
                          other_chain_ptr->copy_of = chain_ptr->chain_id;
                          chain_ptr->ncopies++;

                          /* Update name of current chain */
                          if (!strncmp(chain_ptr->name,"Chain ",6))
                            sprintf(chain_ptr->name,"Chains %c",
                                    chain_ptr->chain_id);

                          /* Add the copy chain's id onto current chain's
                             name */
                          sprintf(add_on,", %c",other_chain_ptr->chain_id);
                          strcat(chain_ptr->name,add_on);

                          /* If chain lengths differ, set flag */
                          if (chain_ptr->nresid != other_chain_ptr->nresid)
                            chain_ptr->same_length = FALSE;
                        }
                    }

                  /* Get pointer to next chain */
                  other_chain_ptr = other_chain_ptr->next_chain_ptr;
                }

              /* If chain is a unique one, then increment count */
              if (chain_ptr->copy_of == '\0')
                nunique++;
            }
        }

      /* Get pointer to next chain */
      chain_ptr = chain_ptr->next_chain_ptr;
    }

  /* Write out the number of unique chains */
  printf("Total number of unique protein chains:  %8d\n",nunique);

  /* Return total number of unique chains */
  return(nunique);
}
/***********************************************************************

write_seq_fasta  -  Write out the seqXXX.fasta, seqXXX.trn and .csq
                    files

***********************************************************************/

void write_seq_fasta(char *pdb_code,int ichain,struct chain *chain_ptr,
                     char *output_dir,int uploaded_file)
{
  char filename[FILENAME_LEN + 1];
  char number_string[5], res_number[6];

  int icode, last_icode;
  int count, last_number, residue_number;
  int iresid, nresid;

  struct residue *residue_ptr;

  FILE *file_ptr, *csq_file_ptr, *trn_file_ptr;

  /* Initialise variables */
  count = 0;
  last_icode = '\0';
  last_number = 0;
  res_number[0] = '\0';

  /* Form name of sequence file */
  sprintf(filename,"%sseq%3.3d.fasta",output_dir,ichain);

  /* If this is an uploaded PDB file and this is its first chain,
     write out as seq.fasta */
  if (uploaded_file == TRUE && ichain == 1)
    sprintf(filename,"%sseq.fasta",output_dir);

  /* Open the sequence file */
  file_ptr = open_output_file(filename,TRUE);

  /* Write out the sequence name */
  if (uploaded_file == TRUE)
    fprintf(file_ptr,">Uploaded file [%c]\n",chain_ptr->chain_id);
  else
    fprintf(file_ptr,">%s:%c\n",pdb_code,chain_ptr->chain_id);

  /* Form name of .trn file */
  sprintf(filename,"%sseq%3.3d.trn",output_dir,ichain);

  /* Open the .trn file */
  trn_file_ptr = open_output_file(filename,TRUE);

  /* Form name of .csq file */
  if (chain_ptr->chain_id != ' ')
    sprintf(filename,"%s%s_%c.csq",output_dir,pdb_code,
            chain_ptr->chain_id);
  else
    sprintf(filename,"%s%s_.csq",output_dir,pdb_code);

  /* Open the .csq file */
  csq_file_ptr = open_output_file(filename,TRUE);

  /* Write out the sequence name */
  fprintf(trn_file_ptr,">%s\n",chain_ptr->name);

  /* Get the first of the chain's residues */
  residue_ptr = chain_ptr->first_residue_ptr;
  nresid = chain_ptr->nresid;
  iresid = 0;

  /* Loop through all the residues writing out the sequence to the
     .fasta file and the corresponding residue numbers to the .trn file */
  while (residue_ptr != NULL && iresid < nresid)
    {
      /* Skip any "ACE" residues */
      if (strcmp(residue_ptr->res_name,"ACE"))
        {
          /* Write this residue out to the FASTA format file */
          fprintf(file_ptr,"%c",residue_ptr->aa_code);

          /* Write the residue number to the .trn file */
          fprintf(trn_file_ptr,"%s",residue_ptr->res_num);

          /* Get the numerical part of this residue number */
          strncpy(number_string,residue_ptr->res_num,4);
          number_string[4] = '\0';
          residue_number = atoi(number_string);

          /* Get the insertion code */
          icode = residue_ptr->res_num[4];

          /* If this is not the next sequential residue number, then we
             have a jump */
          if (icode != last_icode || residue_number != (last_number + 1))
            {
              /* If have a residue to write out, then do so */
              if (count > 0)
                fprintf(csq_file_ptr,"%s%d\n",res_number,count);

              /* Save current residue number and zero the count */
              strcpy(res_number,residue_ptr->res_num);
              count = 0;
            }

          /* Increment the count */
          count++;

          /* Save the current residue number and insertion code */
          last_icode = icode;
          last_number = residue_number;
        }

      /* Get the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;
    }

  /* End the lines in .fasta and .trn files */
  fprintf(file_ptr,"\n");
  fprintf(trn_file_ptr,"\n");

  /* If have a residue to write out to the .csq file, then do so */
  if (count > 0)
    fprintf(csq_file_ptr,"%s%d\n",res_number,count);

  /* Close the files */
  fclose(file_ptr);
  fclose(csq_file_ptr);
  fclose(trn_file_ptr);
}
/***********************************************************************

write_seq_html  -  Write out the seqXXX.html file

***********************************************************************/

void write_seq_html(int ichain,struct chain *chain_ptr,char *output_dir)
{
  char filename[FILENAME_LEN + 1];

  int i, len;

  FILE *file_ptr;

  /* Form name of sequence file */
  sprintf(filename,"%sseq%3.3d.html",output_dir,ichain);

  /* Open the file */
  file_ptr = open_output_file(filename,TRUE);

  /* Write out the sequence name */
  // fprintf(file_ptr,">%s\n",chain_ptr->name);

  /* Get the sequence length */
  len = strlen(chain_ptr->sequence);

  /* Loop while writing the sequence out */
  for (i = 0; i < len; i++)
    {
      /* Write out this character */
      fprintf(file_ptr,"%c",chain_ptr->sequence[i]);

      /* If at end of the line, thne throw line */
      if ((i + 1) % MAX_SEQ_LINE == 0)
        fprintf(file_ptr,"\n");
    }

  /* Throw final new-line, if necessary */
  if (len % MAX_SEQ_HTML_LINE != 0)
    fprintf(file_ptr,"\n");

  /* Close the sequence file */
  fclose(file_ptr);
}
/***********************************************************************

write_seq_chain  -  Write out the seqXXX.chain file

***********************************************************************/

void write_seq_chain(int ichain,char chain,char *output_dir)
{
  char filename[FILENAME_LEN + 1];

  FILE *file_ptr;

  /* Form name of file */
  sprintf(filename,"%sseq%3.3d.chain",output_dir,ichain);

  /* Open the file */
  file_ptr = open_output_file(filename,TRUE);

  /* Write out this sequence's chain id */
  fprintf(file_ptr,"set chain = \"%c\"\n",chain);

  /* Close the sequence file */
  fclose(file_ptr);
}
/***********************************************************************

write_seq_pdb  -  Write out the seqXXX.pdb file

***********************************************************************/

void write_seq_pdb(int ichain,struct chain *ref_chain_ptr,
                   struct chain *first_chain_ptr,char *output_dir,
                   int uploaded_file)
{
  char filename[FILENAME_LEN + 1];

  int i, iatom, iresid, natoms, nresid;
  int toofar, wanted;

  struct atom *atom_ptr;
  struct chain *chain_ptr;
  struct residue *residue_ptr;

  FILE *file_ptr;

  /* Form name of file */
  sprintf(filename,"%sseq%3.3d.pdb",output_dir,ichain);

  /* Open the file */
  file_ptr = open_output_file(filename,TRUE);

  /* If have an uploaded file get pointer to the first chain */
  if (uploaded_file == TRUE)
    chain_ptr = first_chain_ptr;

  /* Otherwise, get the required chain */
  else
    chain_ptr = ref_chain_ptr;

  /* If chain wanted, write out its atom coords */
  while (chain_ptr != NULL)
    {
      /* Determine whether this chain is wanted */
      wanted = TRUE;
      if (uploaded_file == TRUE)
        {
          /* If this is the reference chain, then by definition it is
             wanted */
          if (chain_ptr == ref_chain_ptr)
            wanted = TRUE;

          /* If this is any other protein chain it is not wanted */
          else if (chain_ptr->type == PROTEIN ||
                   chain_ptr->type == CALPHA_ONLY)
            wanted = FALSE;

          /* If this is a residue attachment, keep only if it belongs
             to the right chain */
          else if (chain_ptr->type == ATTACHMENT)
            {
              if (chain_ptr->chain_id == ref_chain_ptr->chain_id)
                wanted = TRUE;
              else
                wanted = FALSE;
            }

          /* Otherwise, for DNA, metal or ligand, check whether chain
             is within range of the reference chain */
          else
            {
              /* Check whether these chains are within reach */
              toofar = FALSE;
              for (i = 0; i < 3 && toofar == FALSE; i++)
                {
                  if (chain_ptr->coord_min[i]
                      - ref_chain_ptr->coord_max[i] > TOO_FAR_CONTACT)
                    toofar = TRUE;
                  if (ref_chain_ptr->coord_min[i]
                      - chain_ptr->coord_max[i] > TOO_FAR_CONTACT)
                    toofar = TRUE;
                }

              /* If chains too far apart, then don't want current one */
              if (toofar == TRUE)
                wanted = FALSE;
            }
        }

      /* If chain is wanted, write out */
      if (wanted == TRUE)
        {
          /* Get pointer to the first residue */
          residue_ptr = chain_ptr->first_residue_ptr;
          nresid = chain_ptr->nresid;
          iresid = 0;

          /* Loop through all the residues to write out */
          while (residue_ptr != NULL && iresid < nresid)
            {
              /* Get number of atoms */
              natoms = residue_ptr->natoms;

              /* Initialise pointer to first atom of this residue */
              atom_ptr = residue_ptr->first_atom_ptr;
              iatom = 0;
  
              /* Loop through all the residue's atoms to look for match */
              while (iatom < natoms && atom_ptr != NULL)
                {
                  /* Write out this atom */
                  fprintf(file_ptr,"ATOM  %5d %s %s %c%s   %8.3f%8.3f%8.3f",
                          atom_ptr->atom_number,atom_ptr->atom_name,
                          atom_ptr->residue_ptr->res_name,
                          atom_ptr->residue_ptr->chain,
                          atom_ptr->residue_ptr->res_num,atom_ptr->coords[0],
                          atom_ptr->coords[1],atom_ptr->coords[2]);
                  fprintf(file_ptr,"%6.2f%6.2f\n",atom_ptr->occupancy,
                          atom_ptr->bvalue);

                  /* Get pointer to next atom in linked-list */
                  atom_ptr = atom_ptr->next_atom_ptr;
                  iatom++;
                }

              /* Get the next residue */
              residue_ptr = residue_ptr->next_residue_ptr;
              iresid++;
            }
        }
      /* Get pointer to the next chain */
      chain_ptr = chain_ptr->next_chain_ptr;

      /* If this is not an uploaded file, we are done */
      if (uploaded_file == FALSE)
        chain_ptr = NULL;
    }

  /* Close the sequence file */
  fclose(file_ptr);
}
/***********************************************************************

write_chain_equivs  -  Write out the chain equivalences to the chains.dat
                       file

***********************************************************************/

void write_chain_equivs(struct chain *first_chain_ptr,
                        char *chains_dat,char *output_dir,
                        struct parameters params)
{
  char last_chain;

  int nchains;

  struct chain *chain_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  last_chain = '\0';
  nchains = 0;

  /* Open the output chains.dat file */
  file_ptr = open_output_file(chains_dat,TRUE);

  /* Get pointer to first chain */
  chain_ptr = first_chain_ptr;

  /* Loop over all chains to write out their groupings */
  while (chain_ptr != NULL)
    {
      /* Only process if this is a protein sequence */
      if ((chain_ptr->type == PROTEIN || chain_ptr->type == CALPHA_ONLY) &&
          chain_ptr->chain_id != last_chain)
        {
          /* Write out the chain and which chain it is a copy of */
          fprintf(file_ptr,"%c",chain_ptr->chain_id);
          if (chain_ptr->copy_of != '\0')
            fprintf(file_ptr," %c",chain_ptr->copy_of);
          fprintf(file_ptr,"\n");

          /* Check that chain isn't a copy of another chain */
          if (chain_ptr->copy_of == '\0')
            {
              /* Increment chain count */
              nchains++;

              /* Write out the seqXXX.fasta file */
              write_seq_fasta(params.pdb_code,nchains,chain_ptr,
                              output_dir,params.uploaded_file);

              /* Write out the seqXXX.html file */
              write_seq_html(nchains,chain_ptr,output_dir);

              /* Write out the seqXXX.chain file */
              write_seq_chain(nchains,chain_ptr->chain_id,output_dir);

              /* Write out the seqXXX.pdb file */
              write_seq_pdb(nchains,chain_ptr,first_chain_ptr,output_dir,
                            params.uploaded_file);
            }

          /* Save the current chain */
          last_chain = chain_ptr->chain_id;
        }

      /* Get pointer to next chain */
      chain_ptr = chain_ptr->next_chain_ptr;
    }

  /* Close the output file */
  fclose(file_ptr);
}
/***********************************************************************

create_ligand_record  -  Create and initialise a new ligand record

***********************************************************************/

struct ligand *create_ligand_record(struct ligand **fst_ligand_ptr,
                                    struct ligand **lst_ligand_ptr,
                                    struct chain *chain_ptr,int number,
                                    int metal)
{
  int i;

  struct ligand *ligand_ptr, *first_ligand_ptr, *last_ligand_ptr;

  /* Initialise ligand pointers */
  ligand_ptr = NULL;
  first_ligand_ptr = *fst_ligand_ptr;
  last_ligand_ptr = *lst_ligand_ptr;

  /* Allocate memory for structure to hold ligand info */
  ligand_ptr = (struct ligand *) malloc(sizeof(struct ligand));
  if (ligand_ptr == NULL)
    {
      printf("*** Can't allocate memory for struct ligand\n");
      exit (1);
    }

  /* If this is the very first ligand, then store pointer */
  if (first_ligand_ptr == NULL)
    first_ligand_ptr = ligand_ptr;

  /* Add link from previous ligand to the current one */
  if (last_ligand_ptr != NULL)
    last_ligand_ptr->next_ligand_ptr = ligand_ptr;
  last_ligand_ptr = ligand_ptr;

  /* Store the ligand details */
  ligand_ptr->number = number;
  ligand_ptr->metal = metal;
  ligand_ptr->metal_colour = 0;
  ligand_ptr->metal_name[0] = '\0';
  ligand_ptr->type = 0;
  ligand_ptr->count = 1;
  ligand_ptr->unique = TRUE;
  ligand_ptr->unique_count = 0;
  ligand_ptr->ncopy = 0;
  ligand_ptr->duplicate_count = 0;
  ligand_ptr->nduplicates = 0;
  ligand_ptr->residue_list = &null;
  ligand_ptr->residue_sequence = &null;
  ligand_ptr->residue_range = &null;
  ligand_ptr->nresid = chain_ptr->nresid;
  ligand_ptr->natoms = chain_ptr->natoms;
  ligand_ptr->sort_score = 0;
  ligand_ptr->got_bonds = FALSE;
  ligand_ptr->first_residue_ptr = chain_ptr->first_residue_ptr;
  ligand_ptr->first_list_ptr = NULL;
  ligand_ptr->duplicate_ligand_ptr = NULL;
  ligand_ptr->next_ligand_ptr = NULL;

  /* Copy chain's coordinate limits */
  for (i = 0; i < 3; i++)
    {
      ligand_ptr->coord_min[i] = chain_ptr->coord_min[i];
      ligand_ptr->coord_max[i] = chain_ptr->coord_max[i];
    }

  /* Return current pointers */
  *fst_ligand_ptr = first_ligand_ptr;
  *lst_ligand_ptr = last_ligand_ptr;

  /* Return new ligand pointer */
  return(ligand_ptr);
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
  if (len + strlen(residue_sequence) > LINELEN - 1)
    {
      printf("*** ERROR. Ligand residue sequence too long! Program "
             "aborted.\n");
      exit(-1);
    }

  /* Add onto name of current residue sequence */
  if (residue_sequence[0] != '\0')
    strcat(residue_sequence,"-");
  strcat(residue_sequence,residue_name);

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

  /* If residue list too long, then abort this entry */
  if (strlen(residue_name) + strlen(residue_list) > LINELEN - 3)
    {
      printf("*** ERROR. Ligand residue list too long! Program "
             "aborted.\n");
      printf("***        Residue list %s + residue name %s\n",
             residue_list,residue_name);
      exit(-1);
    }

  /* Add onto end of current residue list */
  if (residue_list[0] != '\0')
    strcat(residue_list,", ");
  strcat(residue_list,residue_name);
}
/***********************************************************************

format_residue_no  -  Format the ligand's residue range

***********************************************************************/

void format_residue_no(char *residue_number,
                           struct residue *residue_ptr)
{
  char chain_id[4], tmp_residue_number[15];
  char *string_ptr;

  int ipos, len;
  int done;

  /* Initialise variables */
  done = FALSE;
  residue_number[0] = '\0';

  /* Add residue name and number */
  strcat(residue_number,residue_ptr->res_name);
  strcat(residue_number," ");
  strcat(residue_number,residue_ptr->res_num);

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
    {
      strcpy(tmp_residue_number,residue_number+1);
      strcpy(residue_number,tmp_residue_number);
    }

  /* Ditto for final space */
  len = strlen(residue_number);
  if (residue_number[len - 1] == ' ')
    residue_number[len - 1] = '\0';

  /* If chain-id not blank, then add */
  if (residue_ptr->chain != ' ')
    {
      /* Chain id enclosed in brackets */
      sprintf(chain_id,"(%c)",residue_ptr->chain);
      strcat(residue_number,chain_id);
    }
}
/***********************************************************************

form_ligands  -  Create all the ligand records for each of the ligand
                 chains

***********************************************************************/

int form_ligands(struct chain *first_chain_ptr,
                 struct ligand **fst_ligand_ptr,int *nlig,int *nmet,
                 int biomol)
{
  char residue_list[LINELEN + 1], residue_range[LINELEN + 1];
  char residue_sequence[LINELEN + 1];
  char start_residue[15], end_residue[15];

  int iresid, nligands, nligmet, nmetals, nresid;
  int metal;

  struct atom *atom_ptr;
  struct chain *chain_ptr;
  struct residue *residue_ptr, *end_residue_ptr, *start_residue_ptr;
  struct ligand *ligand_ptr, *first_ligand_ptr, *last_ligand_ptr;

  /* Initialise variables */
  nligands = 0;
  nligmet = 0;
  nmetals = 0;

  /* Initialise ligand pointers */
  ligand_ptr = NULL;
  first_ligand_ptr = NULL;
  last_ligand_ptr = NULL;

  /* Get pointer to the first chain */
  chain_ptr = first_chain_ptr;

  /* Loop over all the currently defined chains */
  while (chain_ptr != NULL)
    {
      /* Only process if this is a ligand or metal chain */
      if (chain_ptr->type == LIGAND || chain_ptr->type == METAL)
        {
          /* Determine which type this is */
          if (chain_ptr->type == LIGAND)
            metal = FALSE;
          else
            metal = TRUE;

          /* Initialise variables */
          start_residue_ptr = NULL;
          end_residue_ptr = NULL;

          /* Create a ligand record */
          ligand_ptr
            = create_ligand_record(&first_ligand_ptr,
                                   &last_ligand_ptr,chain_ptr,nligands,
                                   metal);

          /* Increment count of ligand residues read in */
          if (metal == TRUE)
            nmetals++;
          else
            nligands++;
          nligmet++;

          /* If this is a biomolecule ligand, assign number to type */
          if (biomol == TRUE && metal == FALSE)
            ligand_ptr->type = nligands;
          else if (biomol == TRUE && metal == TRUE)
            ligand_ptr->type = nmetals;

          /* Initialise residue list and sequence */
          residue_list[0] = '\0';
          residue_sequence[0] = '\0';

          /* Get the first of this chain's residues */
          residue_ptr = chain_ptr->first_residue_ptr;
          nresid = chain_ptr->nresid;
          iresid = 0;
          
          /* If a metal, get its name */
          if (ligand_ptr->metal == TRUE && residue_ptr != NULL )
            {
              /* Get the atom corresponding to this metal ion */
              atom_ptr = ligand_ptr->first_residue_ptr->first_atom_ptr;

              /* Store the metal name */
              strcpy(ligand_ptr->metal_name,atom_ptr->element);
            }

          /* Loop through the chain's residues to add to current ligand */
          while (residue_ptr != NULL && iresid < nresid)
            {
              /* Add current residue name to the ligand residue list */
              add_res_name(residue_sequence,residue_list,
                           residue_ptr->res_name,residue_ptr->res_num,
                           residue_ptr->chain);

              /* If this is the first or last residue of the ligand,
                 store */
              if (iresid == 0)
                start_residue_ptr = residue_ptr;
              if (iresid == nresid - 1)
                end_residue_ptr = residue_ptr;

              /* Store ligand pointer */
              residue_ptr->ligand_ptr = ligand_ptr;

              /* Get the next residue */
              residue_ptr = residue_ptr->next_residue_ptr;
              iresid++;
            }

          /* Store the ligand's residue list and sequence */
          store_string(&(ligand_ptr->residue_list),residue_list);
          store_string(&(ligand_ptr->residue_sequence),residue_sequence);

          /* Format the ligand's residue range */
          format_residue_no(start_residue,start_residue_ptr);
          format_residue_no(end_residue,end_residue_ptr);

          /* Form the range */
          strcpy(residue_range,start_residue);
          if (nresid > 1)
            {
              strcat(residue_range," to ");
              strcat(residue_range,end_residue);
            }

          /* Store the ligand's residue range */
          store_string(&(ligand_ptr->residue_range),residue_range);
        }

      /* Get the next chain */
      chain_ptr = chain_ptr->next_chain_ptr;
    }

  /* Write out the number of ligands */
  printf("Total number of ligand molecules:       %8d\n",nligands);
  printf("Total number of metals:                 %8d\n",nmetals);

  /* Return pointer to the first ligand */
  *fst_ligand_ptr = first_ligand_ptr;

  /* Return total numbers of ligand and metals */
  *nlig = nligands;
  *nmet = nmetals;
  return(nligmet);
}
/***********************************************************************

get_ligand_neighbours  -  Identify the protein residues that are within
                          interaction range of each ligand

***********************************************************************/

void get_ligand_neighbours(struct ligand *first_ligand_ptr,
                           struct chain *first_chain_ptr)
{
  int i, in_range, iresid, lresid, nlresid, nresid;
  int toofar, wanted;

  struct chain *chain_ptr;
  struct ligand *ligand_ptr;
  struct list *list_ptr, *last_list_ptr;
  struct residue *residue_ptr, *lig_residue_ptr;

  /* Initialise variables */

  /* Get pointer to the first ligand record */
  ligand_ptr = first_ligand_ptr;

  /* Loop through all the ligands */
  while (ligand_ptr != NULL)
    {
      /* Initialise variables */
      last_list_ptr = NULL;
      wanted = TRUE;

      /* Get pointer to the first chain */
      chain_ptr = first_chain_ptr;

      /* Loop over all the chains */
      while (chain_ptr != NULL)
        {
          /* Check whether ligand is too far from this chain */
          if (chain_ptr->type == PROTEIN)
            {
              /* Check whether ligand is within bounds of
                 this chain */
              toofar = FALSE;
              for (i = 0; i < 3 && toofar == FALSE; i++)
                {
                  if (ligand_ptr->coord_min[i] - chain_ptr->coord_max[i]
                      > TOO_FAR_DIST)
                    toofar = TRUE;
                  if (chain_ptr->coord_min[i] - ligand_ptr->coord_max[i]
                      > TOO_FAR_DIST)
                    toofar = TRUE;
                }
            }

          /* Only consider if this is a protein chain */
          if (chain_ptr->type == PROTEIN && toofar == FALSE)
            {
              /* Get pointer to first of this protein chain's residues */
              residue_ptr = chain_ptr->first_residue_ptr;
              nresid = chain_ptr->nresid;
              iresid = 0;

              /* Loop through all the protein residues to identify all that
                 lie close to the ligand */
              while (residue_ptr != NULL && iresid < nresid)
                {
                  /* Check whether this residue is within bounds of
                     the ligand */
                  toofar = FALSE;
                  for (i = 0; i < 3 && toofar == FALSE; i++)
                    {
                      if (ligand_ptr->coord_min[i] - residue_ptr->coord_max[i]
                          > TOO_FAR_DIST)
                        toofar = TRUE;
                      if (residue_ptr->coord_min[i] - ligand_ptr->coord_max[i]
                          > TOO_FAR_DIST)
                        toofar = TRUE;
                    }

                  /* If residue within range of ligand bounds, check
                     whether it is also within range of each of the
                     ligand's residue bounds */
                  if (toofar == FALSE && ligand_ptr->nresid > 1)
                    {
                      /* Get pointer to first of the ligand residues */
                      lig_residue_ptr = ligand_ptr->first_residue_ptr;
                      nlresid = ligand_ptr->nresid;
                      lresid = 0;
                      toofar = TRUE;

                      /* Loop through all the ligand residues to determine
                         whether current residue is within range of one or
                         more */
                      while (lig_residue_ptr != NULL && lresid < nlresid &&
                             toofar == TRUE)
                        {
                          /* Initialise flag */
                          in_range = 0;

                          /* Check whether the residue bounds are within
                             range */
                          for (i = 0; i < 3; i++)
                            {
                              if ((lig_residue_ptr->coord_min[i] >
                                   residue_ptr->coord_min[i] - TOO_FAR_DIST &&
                                   lig_residue_ptr->coord_min[i] <
                                   residue_ptr->coord_max[i] + TOO_FAR_DIST) ||
                                  (lig_residue_ptr->coord_max[i] >
                                   residue_ptr->coord_min[i] - TOO_FAR_DIST &&
                                   lig_residue_ptr->coord_max[i] <
                                   residue_ptr->coord_max[i] + TOO_FAR_DIST))
                                in_range++;
                            }

                          /* If all coords within range of this residue, then
                             residue is a neighbour of the ligand */
                          if (in_range == 3)
                            toofar = FALSE;

                          /* Get the next ligand residue */
                          lig_residue_ptr = lig_residue_ptr->next_residue_ptr;
                          lresid++;
                        }
                    }

                  /* If residue within range, store pointer to it */
                  if (toofar == FALSE)
                    {
                      /* Create a list record containing a pointer to the
                         stored residue */
                      list_ptr
                        = create_list_record(&(ligand_ptr->first_list_ptr),
                                             &last_list_ptr,residue_ptr);
                    }

                  /* Get the next residue */
                  residue_ptr = residue_ptr->next_residue_ptr;
                  iresid++;
                }
            }

          /* Get the next chain */
          chain_ptr = chain_ptr->next_chain_ptr;
        }

      /* Get next ligand record */
      ligand_ptr = ligand_ptr->next_ligand_ptr;
    }
}
/***********************************************************************

count_equivs  -  Count the number of equivalent residues that the two
                 ligands intract with

***********************************************************************/

int count_equivs(struct ligand *ligand_ptr,
                 struct ligand *other_ligand_ptr,
                 struct chain *first_chain_ptr,int *nint)
{
  char chain1, chain2;

  int nequiv, ninteract1, ninteract2;
  int chain_match, matched;

  struct chain *chain_ptr;
  struct list *list_ptr, *other_list_ptr;
  struct residue *residue_ptr, *other_residue_ptr;

  /* Initialise variables */
  nequiv = 0;
  ninteract1 = 0;

  /* Get first ligand's first residue */
  list_ptr = ligand_ptr->first_list_ptr;

  /* Loop over the residues in the list to compare */
  while (list_ptr != NULL)
    {
      /* Get this residue */
      residue_ptr = list_ptr->residue_ptr;
      if (residue_ptr->ligand_ptr != ligand_ptr)
        ninteract1++;

      /* Get other ligand's first residue */
      other_list_ptr = other_ligand_ptr->first_list_ptr;
      matched = FALSE;
      ninteract2 = 0;

      /* Loop until a match is found */
      while (other_list_ptr != NULL && matched == FALSE)
        {
          /* Get this residue */
          other_residue_ptr = other_list_ptr->residue_ptr;
          if (other_residue_ptr->ligand_ptr != other_ligand_ptr)
            ninteract2++;

          /* Check for a match on the residue name and number */
          if (!strcmp(residue_ptr->res_num,other_residue_ptr->res_num) &&
              !strcmp(residue_ptr->res_name,other_residue_ptr->res_name))
            {
              /* Check whether the chains are equivalent */
              chain1 = residue_ptr->chain;
              chain2 = other_residue_ptr->chain;

              /* Assume we have no match */
              chain_match = FALSE;

              /* Get first chain pointer */
              chain_ptr = first_chain_ptr;

              /* Loop through the chains to see if we find an
                 equivalence */
              while (chain_ptr != NULL && chain_match == FALSE)
                {
                  /* Check for the equivalence */
                  if ((chain1 == chain_ptr->chain_id &&
                       chain2 == chain_ptr->copy_of) ||
                      (chain1 == chain_ptr->copy_of &&
                       chain2 == chain_ptr->chain_id))
                    chain_match = TRUE;

                  /* Get the next chain pointer */
                  chain_ptr = chain_ptr->next_chain_ptr;
                }

              /* If have a chain match, then increment count of
                 equivalent residues */
              if (chain_match == TRUE)
                {
                  /* Increment count */
                  nequiv++;

                  /* Set flag that match found */
                  matched = TRUE;
                }
            }

          /* Get the next list item */
          other_list_ptr = other_list_ptr->next_list_ptr;
        }

      /* Get the next list item */
      list_ptr = list_ptr->next_list_ptr;
    }

  /* Store the lower count of interacting residues */
  *nint = ninteract1;
  if (ninteract2 < ninteract1)
    *nint = ninteract2;

  /* Return number of equivalent residues */
  return(nequiv);
}
/***********************************************************************

identify_unique_ligands  -  Determine which ligands are unique

***********************************************************************/

int identify_unique_ligands(struct ligand *first_ligand_ptr,
                            struct chain *first_chain_ptr,int *nmetyp)
{
  char element[3];

  int nequiv, max_equiv, ninter, ninteract, nmetypes, ntypes, nunique;
  int new_type;

  struct atom *atom_ptr;
  struct ligand *ligand_ptr;
  struct ligand *copy_ligand_ptr, *equiv_ligand_ptr, *other_ligand_ptr;

  /* Initialise variables */
  nmetypes = 0;
  ntypes = 0;
  nunique = 0;

  /* Get pointer to the first ligand record */
  ligand_ptr = first_ligand_ptr;

  if (ligand_ptr != NULL)
    {
      /* The first ligand is, by definition, unique */
      ligand_ptr->unique = TRUE;
      if (ligand_ptr->metal == TRUE)
        {
          nmetypes++;
          ligand_ptr->type = nmetypes;
        }
      else
        {
          ntypes++;
          ligand_ptr->type = ntypes;
        }
      nunique++;

      /* Get the next ligand */
      ligand_ptr = ligand_ptr->next_ligand_ptr;
    }

  /* Loop through all the ligands to determine which are unique
     and which are symmetry-related copies */
  while (ligand_ptr != NULL)
    {
      /* Initialise variables */
      max_equiv = 0;
      ninteract = 0;
      copy_ligand_ptr = NULL;
      equiv_ligand_ptr = NULL;
      new_type = TRUE;

      /* If a metal, get its name */
      if (ligand_ptr->metal == TRUE)
        {
          /* Get the atom corresponding to this metal ion */
          atom_ptr = ligand_ptr->first_residue_ptr->first_atom_ptr;

          /* Store the metal name */
          strcpy(ligand_ptr->metal_name,atom_ptr->element);
        }

      /* Get pointer to the first ligand */
      other_ligand_ptr = first_ligand_ptr;

      /* Loop through all the previous ligands to see if this is a
         copy of any one of them */
      while (other_ligand_ptr != NULL &&
             other_ligand_ptr != ligand_ptr)
        {
          /* Only consider this ligand if it is a unique one */
          if (other_ligand_ptr->unique == TRUE)
            {
              /* Check whether these two ligands have the same string
                 of residues */
              if (!strcmp(ligand_ptr->residue_sequence,
                          other_ligand_ptr->residue_sequence))
                {
                  /* Assign same ligand type */
                  ligand_ptr->type = other_ligand_ptr->type;
                  new_type = FALSE;

                  /* Save ligand that this is a copy of */
                  if (copy_ligand_ptr == NULL)
                    copy_ligand_ptr = other_ligand_ptr;

                  /* Count the number of equivalent residues that the
                     two ligands intract with */
                  nequiv = count_equivs(ligand_ptr,other_ligand_ptr,
                                        first_chain_ptr,&ninter);

                  /* If this is the largest number of equivalences so
                     far, save the ligand matched */
                  if (nequiv > max_equiv)
                    {
                      /* Save equiv ligand */
                      equiv_ligand_ptr = other_ligand_ptr;
                      max_equiv = nequiv;
                      ninteract = ninter;
                    }
                }
            }

          /* Get next ligand record */
          other_ligand_ptr = other_ligand_ptr->next_ligand_ptr;
        }

      /* If this is a new ligand type, assign it the next type
         number */
      if (new_type == TRUE)
        {
          /* Increment type count and assign to current ligand */
          if (ligand_ptr->metal == TRUE)
            {
              nmetypes++;
              ligand_ptr->type = nmetypes;
            }
          else
            {
              ntypes++;
              ligand_ptr->type = ntypes;
            }
        }

      /* If have enough matching residues, then current ligand is
         a duplicate */
      if (max_equiv > ninteract / 2)
        {
          /* Mark as non-unique */
          ligand_ptr->unique = FALSE;

          /* Save ligand of which it is a duplicate */
          equiv_ligand_ptr->nduplicates++;

          /* Store this as the duplicate number of this ligand */
          ligand_ptr->unique_count = equiv_ligand_ptr->unique_count;
          ligand_ptr->duplicate_count = equiv_ligand_ptr->nduplicates;
          ligand_ptr->duplicate_ligand_ptr = equiv_ligand_ptr;
        }

      /* Otherwise, deem it to be unique */
      else
        {
          /* Mark as unique */
          ligand_ptr->unique = TRUE;

          /* Get its copy number */
          if (copy_ligand_ptr != NULL)
            {
              /* Increment count of copies of original ligand type */
              copy_ligand_ptr->ncopy++;

              /* Store this as the copy number of this ligand */
              ligand_ptr->unique_count = copy_ligand_ptr->ncopy;
            }
        }

      /* Increment the total count of this ligand type */
      if (copy_ligand_ptr != NULL)
        copy_ligand_ptr->count++;

      /* Get next ligand record */
      ligand_ptr = ligand_ptr->next_ligand_ptr;
    }

  /* Show the number of ligand types */
  printf("Number of different ligand types:       %8d\n",ntypes);
  printf("Number of different metal types:        %8d\n",nmetypes);

  /* Return the numbers of ligand and metal types */
  *nmetyp = nmetypes;
  return(ntypes);
}
/***********************************************************************

shell_sort  -  Shell-sort which sorts the ligand types into ascending
               order ot type, copy and duplicate

***********************************************************************/

void shell_sort(int *index,struct ligand **ligand_index_ptr,int nindex)
{
  int endloop, i, indi, indl, iswap, j, k, l, lognb2, m, n, nn;
  float aln2i = (float) 1.4426950, tiny = (float) 1.e-5;
  float calc;
  double dble, dble1;

  struct ligand *ligand1_ptr, *ligand2_ptr;

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
              ligand1_ptr = ligand_index_ptr[index[indi]];
              ligand2_ptr = ligand_index_ptr[index[indl]];

              /* If in wrong order, the swap sort-pointers */
              if (ligand1_ptr->sort_score > ligand2_ptr->sort_score)
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

sort_ligands  -  Sort the ligands by their score and delete any of
                 the lower-scoring ones that overlap the higher
                 scoring ones significantly

***********************************************************************/

int sort_ligands(struct ligand **fst_ligand_ptr,int nligands)
{
  int iligand, iorder, ipos;
  int *index;

  struct ligand *first_ligand_ptr, *last_ligand_ptr;
  struct ligand *ligand_ptr;
  struct ligand **ligand_index_ptr;

  /* Initialise pointers */
  first_ligand_ptr = *fst_ligand_ptr;
  last_ligand_ptr = NULL;

  /* Grab memory for temporary sort arrays */
  index = (int *) malloc(nligands * sizeof(int));
  ligand_index_ptr
    = (struct ligand **) malloc(nligands * sizeof(struct ligand *));

  /* Get pointer to first ligand record */
  ligand_ptr = first_ligand_ptr;
  ipos = 0;

  /* Loop through ligand records to initialise the sort index */
  while (ligand_ptr != NULL)
    {
      /* Initialise sort index */
      index[ipos] = ipos;
      ligand_index_ptr[ipos] = ligand_ptr;

      /* Calculate the score on which to sort, involving the type,
         instance and duplicate copy */
      ligand_ptr->sort_score
        = 1000 * 1000 * ligand_ptr->type
        + 1000 * ligand_ptr->unique_count
        + ligand_ptr->duplicate_count;

      /* Want ligands to come ahead of metals */
      if (ligand_ptr->metal == TRUE)
        ligand_ptr->sort_score = ligand_ptr->sort_score + 10000000;

      /* Increment ligand counter */
      ipos++;

      /* Get next ligand record */
      ligand_ptr = ligand_ptr->next_ligand_ptr;
    }

  /* Initialise ligand counter */
  iligand = 0;
  nligands = ipos;

  /* Sort the ligands */
  shell_sort(index,ligand_index_ptr,nligands);

  /* Having sorted the ligands, rearrange order of the ligand-record
     pointers in descending order of ligand */
  for (iorder = 0; iorder < nligands; iorder++)
    {
      /* Get the array position of this ligand record */
      ipos = index[iorder];

      /* Get the current ligand record */
      ligand_ptr = ligand_index_ptr[ipos];

      /* If this is the first, then make the first record */
      if (iorder == 0)
        first_ligand_ptr = ligand_ptr;

      /* Otherwise, get the previous ligand record to point to this one */
      else
        last_ligand_ptr->next_ligand_ptr = ligand_ptr;

      /* Blank out current record's next ligand pointer */
      ligand_ptr->next_ligand_ptr = NULL;

      /* Increment ligand counter and store */
      iligand++;
      ligand_ptr->number = iligand;

      /* Store current record as last encountered */
      last_ligand_ptr = ligand_ptr;
    }

  /* Free up memory used for temporary sort arrays */
  free(index);
  free(ligand_index_ptr);

  /* Return current pointer */
  *fst_ligand_ptr = first_ligand_ptr;

  /* Return number of sorted ligands */
  return(nligands);
}
/***********************************************************************

write_ligmet_run  -  Write out the .run files for ligands and metals

***********************************************************************/

void write_ligmet_run(struct ligand *first_ligand_ptr,
                      FILE *ligmet_file_ptr)
{
  char file_name[FILENAME_LEN + 1], output_dir[FILENAME_LEN + 1];
  char ltype[6];

  struct ligand *ligand_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  output_dir[0] = '\0';

  /* Get pointer to the first ligand record */
  ligand_ptr = first_ligand_ptr;

  /* Loop over the ligands */
  while (ligand_ptr != NULL)
    {
      /* Determine whether metal or ligand */
      if (ligand_ptr->metal == TRUE)
        strcpy(ltype,"met");
      else
        strcpy(ltype,"lig");

      /* If this is a unique ligand, write out the run.scr parameters
         to the .run file */
      if (ligand_ptr->unique == TRUE &&
          ligand_ptr->unique_count < MAX_COPIES)
        {
          /* Form name of .run file */
          sprintf(file_name,"%s%s_%2.2d_%2.2d.run",output_dir,ltype,
                  ligand_ptr->type,ligand_ptr->unique_count + 1);

          /* Write the .run file */
          write_run_scr_params(file_name,ligand_ptr->first_residue_ptr,
                               ligand_ptr->nresid,ligand_ptr->metal,
                               ligand_ptr->metal_colour,ligmet_file_ptr);
        }

      /* Get the next ligand */
      ligand_ptr = ligand_ptr->next_ligand_ptr;
    }
}
/***********************************************************************

add_range  -  Add the range to string, writing out if too long

***********************************************************************/

void add_range(FILE *range_file_ptr,char *residue_ranges,int range,
               int first_res,int last_res)
{
  char current_range[MAX_OUTLEN];

  int len, lenstring;

  /* Form the string to add */
  sprintf(current_range,"\"%d-%d\"",first_res,last_res);

  /* Get the string lengths */
  len = strlen(current_range);
  lenstring = strlen(residue_ranges);

  /* If current string already long, then write out */
  if (lenstring + len + 5 > MAX_OUTLEN)
    {
      /* Write out the range string for the last
         range */
      fprintf(range_file_ptr,"%d:%s\n",range,residue_ranges);

      /* Copy current range to new string */
      strcpy(residue_ranges,current_range);
    }

  /* Otherwise, append current string */
  else
    {
      if (residue_ranges[0] != '\0')
        strcat(residue_ranges,", ");
      strcat(residue_ranges,current_range);
    }
}
/***********************************************************************

write_af_scores  -  Write out the AF_scores.dat file giving the
                    confidence scores for each residue

***********************************************************************/

void write_af_scores(struct residue *first_residue_ptr,char *pdbsum_dir)
{
  char file_name[FILENAME_LEN + 1];
  char residue_ranges[NAF_RANGES][MAX_OUTLEN];

  int first_res, i, iresid, range, last_range, last_res;

  double bvalue;

  struct atom *atom_ptr;
  struct residue *residue_ptr;

  FILE *file_ptr, *range_file_ptr;

  /* Initialise variables */
  last_range = -1;
  first_res = last_res = -1;

  /* Initialise ranges */
  for (i = 0; i < NAF_RANGES; i++)
    residue_ranges[i][0] = '\0';

  /* Form output file name of AF_scores.dat */
  strcpy(file_name,pdbsum_dir);
  strcat(file_name,"AF_scores.dat");

  /* Open the rasmol.dfn output file */
  file_ptr = open_output_file(file_name,TRUE);

  /* Repeat for AF_ranges.dat */
  strcpy(file_name,pdbsum_dir);
  strcat(file_name,"AF_ranges.dat");

  /* Open the rasmol.dfn output file */
  range_file_ptr = open_output_file(file_name,TRUE);

  /* Get pointer to the first residue */
  residue_ptr = first_residue_ptr;

  /* Loop over all the residues */
  while (residue_ptr != NULL)
    {
      /* Get the first atom */
      atom_ptr = residue_ptr->first_atom_ptr;

      /* If valid, get the B-value */
      if (atom_ptr != NULL)
        {
          /* Get the B-value */
          bvalue = atom_ptr->bvalue;

          /* Determine which range this score falls into */
          range = -1;
          for (i = 0; i < NAF_RANGES && range == -1; i++)
            {
              if (bvalue < (double) AF_SCORE[i])
                range = i;
            }
          if (range == -1)
            range = NAF_RANGES - 1;

          /* Write to the output file */
          fprintf(file_ptr,"%s %s %c %8.3f %d\n",residue_ptr->res_name,
                  residue_ptr->res_num,residue_ptr->chain,bvalue,range);

          /* Save the range */
          residue_ptr->range = range;

          /* Get the residue number */
          iresid = atoi(residue_ptr->res_num);

          /* If this is a new range, write out last one */
          if (range != last_range)
            {
              /* If we have a previous range, then append to list */
              if (last_range > -1)
                add_range(range_file_ptr,residue_ranges[last_range],
                          last_range,first_res,last_res);

              /* Save the range */
              last_range = range;

              /* Save current residue number as start of new range */
              first_res = iresid;
            }

          /* Save current residue number as last of range */
          last_res = iresid;
        }

      /* Get the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
    }

  /* If have a last range, then process it */
  if (last_range > -1)
    add_range(range_file_ptr,residue_ranges[last_range],last_range,
              first_res,last_res);

  /* Write out the last parts of all ranges */
  for (i = 0; i < NAF_RANGES; i++)
    {
      if (residue_ranges[i][0] != '\0')
        fprintf(range_file_ptr,"%d:%s\n",i,residue_ranges[i]);
    }

  /* Close the output files */
  fclose(file_ptr);
  fclose(range_file_ptr);
}
/***********************************************************************

split_residue_sequence  -  Split up the residue list into smaller chunks

***********************************************************************/

int split_residue_sequence(char residue_sequence[MAX_LINES][MAX_LEN + 1],
                           char *full_list)
{
  char ch;

  int ipos, last_hyphen, len, line_len, nlines, start_pos;
  int got_line;

  /* Initialise variables */
  nlines = 0;
  ipos = 0;
  len = strlen(full_list);

  /* Loop until list has been split up */
  while (ipos < len - 1 && nlines < MAX_LINES)
    {
      /* Initialise position of last hyphen */
      last_hyphen = -1;
      line_len = 0;
      got_line = FALSE;
      start_pos = ipos;

      /* Scan list until required line length obtained */
      while (got_line == FALSE)
        {
          /* Get current character */
          ch = full_list[ipos];

          /* If maximum line length exceeded, then terminate at last
             hyphen */
          if (line_len > MAX_LEN - 1 || ipos == len)
            {
              /* If no hyphen encountered, then end line here */
              if (ipos == len || last_hyphen < 1)
                {
                  residue_sequence[nlines][line_len] = '\0';
                }

              /* Otherwise, end after the last hyphen position */
              else
                {
                  /* Get lst hyphen position */
                  residue_sequence[nlines][last_hyphen + 1] = '\0';

                  /* Save the resultant line length */
                  line_len = last_hyphen + 1;
                }

              /* Increment line count */
              nlines++;

              /* Set flag that we have a line */
              got_line = TRUE;

              /* Reset starting position */
              ipos = start_pos + line_len;
              start_pos = ipos;
            }

          /* Otherwise, copy character across */
          else
            {
              /* Store current character */
              residue_sequence[nlines][line_len] = ch;

              /* If this is a hyphen, store its position */
              if (ch == '-')
                last_hyphen = line_len;

              /* Increment character position and line length */
              ipos++;
              line_len++;
            }
        }
    }

  /* Return number of lines */
  return(nlines);
}
/***********************************************************************

write_het_groups  -  Write out the Het Groups and their descriptions
                     to the given file

***********************************************************************/

void write_het_groups(FILE *file_ptr,struct residue *first_residue_ptr,
                      int nresid,struct hetdesc *first_hetdesc_ptr,
                      char *lighet,char *item1)
{
  char item2[10];

  int ihet, iresid, nhets;

  struct hetdesc *hetdesc_ptr;
  struct residue *residue_ptr;

  /* Initialise variables */
  nhets = 0;

  /* Get pointer to the first hetdesc record */
  hetdesc_ptr = first_hetdesc_ptr;

  /* Loop over all the hetdesc records to initialise counts */
  while (hetdesc_ptr != NULL)
    {
      /* Zero count */
      hetdesc_ptr->count = 0;

      /* Get the next hetdesc record */
      hetdesc_ptr = hetdesc_ptr->next_hetdesc_ptr;
    }

  /* Get pointer to the first residue record */
  residue_ptr = first_residue_ptr;
  iresid = 0;

  /* Loop over all the residue records to find Het Groups */
  while (residue_ptr != NULL && iresid < nresid)
    {
      /* If this is a Het Group, increment count */
      if (residue_ptr->hetdesc_ptr != NULL)
        {
          /* Get the appropriate Het Group */
          hetdesc_ptr = residue_ptr->hetdesc_ptr;

          /* If this is the first occurrence of this Het Group,
             increment count of type */
          if (hetdesc_ptr->count == 0)
            nhets++;

          /* Increment count */
          hetdesc_ptr->count++;
        }

      /* Get the next residue record */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;
    }

  /* If have any Het Groups, write out the details */
  if (nhets > 0)
    {
      /* Write out the count of Het Groups */
      fprintf(file_ptr,"N%sGROUPS%s %d\n",lighet,item1,nhets);

      /* Initialise Het Group count */
      ihet = 0;

      /* Get pointer to the first hetdesc record */
      hetdesc_ptr = first_hetdesc_ptr;

      /* Loop over all the hetdesc records to initialise counts */
      while (hetdesc_ptr != NULL)
        {
          /* If count for this group is non-zero, write out */
          if (hetdesc_ptr->count > 0)
            {
              /* Increment count of Het Groups */
              ihet++;

              /* Form the item number for this line */
              sprintf(item2,"[%d]",ihet);

              /* Write out the Het Group name */
              fprintf(file_ptr,"%s_NAME%s%s %s\n",lighet,item1,item2,
                      hetdesc_ptr->res_name);

              /* Write out description, if there is one */
              if (hetdesc_ptr->details[DESCRIPTION][0] != '\0')
                fprintf(file_ptr,"%s_DESC%s%s %s\n",lighet,item1,item2,
                        hetdesc_ptr->details[DESCRIPTION]);

              /* Write out synonym, if there is one */
              if (hetdesc_ptr->details[SYNONYM][0] != '\0')
                fprintf(file_ptr,"%s_SYNONYM%s%s %s\n",lighet,item1,item2,
                        hetdesc_ptr->details[SYNONYM]);

              /* Write out formula, if there is one */
              if (hetdesc_ptr->details[FORMULA][0] != '\0')
                fprintf(file_ptr,"%s_FORMULA%s%s %s\n",lighet,item1,item2,
                        hetdesc_ptr->details[FORMULA]);
            }

          /* Get the next hetdesc record */
          hetdesc_ptr = hetdesc_ptr->next_hetdesc_ptr;
        }
    }
}
/***********************************************************************

write_raster3d_ligand  -  Write out the given ligand to the Raster3D
                          output file

***********************************************************************/

double write_raster3d_ligand(FILE *file_ptr,struct ligand *ligand_ptr,
                             double scale_factor)
{
  int atom_number, iatom, natoms;

  double min_z;
  static double rgb[3] = { 0.0, 0.0, 0.0 };

  struct atom *atom_ptr;
  struct conect *conect_ptr;
  struct residue *residue_ptr;

  /* Initialise variables */
  min_z = 0.0;

  /* Get pointer to the first atom record */
  atom_ptr = ligand_ptr->first_residue_ptr->first_atom_ptr;
  natoms = ligand_ptr->natoms;
  iatom = 0;

  /* Loop through all the residue's atoms */
  while (atom_ptr != NULL && iatom < natoms)
    {
      /* Get the current atom number */
      atom_number = atom_ptr->atom_number;

      /* Write out atom as a sphere in Raster3D format */
      min_z = write_raster3d_atom(file_ptr,atom_ptr,scale_factor,
                                  LIGAND_ATOM_PERCEN,min_z);

      /* Get the first of this atom's covalent bonds */
      conect_ptr = atom_ptr->first_conect_ptr;

      /* Loop over this atom's covalent bonds and write out */
      while (conect_ptr != NULL)
        {
          /* If this connects to a higher-numbered atom, then
             want this bond */
          if (conect_ptr->to_atom_ptr->atom_number > atom_number)
            {
              /* Get this atom's residue */
              residue_ptr = conect_ptr->to_atom_ptr->residue_ptr;

              /* Check that the atom belongs to the ligand and is not
                 from a covalent bond to protein/DNA/metal */
              if (residue_ptr->ligand_ptr == ligand_ptr)
                {
                  /* Write out this bond */
                  write_raster3d_bond(file_ptr,atom_ptr,
                                      conect_ptr->to_atom_ptr,scale_factor,
                                      LIGAND_BOND_RADIUS,TRUE,rgb,TRUE,
                                      FALSE);
                }
            }

          /* Get the next conect */
          conect_ptr = conect_ptr->next_conect_ptr;
        }
  
      /* Get pointer to next atom in linked-list */
      atom_ptr = atom_ptr->next_atom_ptr;
      iatom++;
    }

  /* Return minimum z-value for placing of back wall */
  return(min_z);
}
/***********************************************************************

run_pymol  -  Run PyMol to generate all the image files defined in the
              pymol.pml script file

***********************************************************************/

void run_pymol(char *file,char *pymol_exe,int dos)
{
  char command[FILENAME_LEN + 1], tmp[FILENAME_LEN + 1];

  /* Form the parameters to run PyMol */
  printf("Running PyMol on %s.pml file ...\n",file);
  sprintf(command,"%s -c %s.pml > pymol.log",pymol_exe,file);

  /* Under DOS, start a new command window to prevent PyMOL closing the
     current one */
  if (dos == TRUE)
    {
      sprintf(tmp,"start \"%s\" ",file);
      strcat(tmp,command);
      //      strcpy(command,tmp);
    }
  
  /* Run PyMol */
  call_system(command,dos);

  /* For DOS, pause to allow PyMOL to finish */
  if (dos == TRUE)
    {
      strcpy(command,"timeout /t 10");
      //      system(command);
    }
}
/***********************************************************************

remove_residue_blanks  -  Remove all blanks from the given residue range

***********************************************************************/

int remove_residue_blanks(char *string)
{
  char *tmp;

  int i, len, nhyphens;
  int done, valid;

  /* Initialise variables */
  done = FALSE;
  i = 0;
  nhyphens = 0;
  valid = TRUE;

  /* Get the length of the string */
  len = strlen(string);
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
        {
          /* If this is a hyphen, increment hyphen count */
          if (string[i] == '-')
            nhyphens++;

          /* Otherwise, if not a number flag as invalid */
          else if (string[i] < '0' || string[i] > '9')
            valid = FALSE;

          /* Increment position */
          i++;
        }
    }

  /* If don't have one hyphen, then range is not valid */
  if (nhyphens != 1)
    valid = FALSE;

  /* Free up memory taken up by temporary string */
  free(tmp);

  /* Return whether range is valid */
  return(valid);
}
/***********************************************************************

render_pymol_ligand  -  Write PyMol script to render the given ligand as
                        sticks

***********************************************************************/

void render_pymol_ligand(char *pdb_file,struct ligand *ligand_ptr,
                         char *pymol_exe,int dos)
{
  char chain, file_name[FILENAME_LEN], png_name[FILENAME_LEN];
  char name[3], res_num1[6], res_num2[6], residue_range[30], tmp[30];

  int colour, iresid, nresid, rnum1, rnum2;
  int valid;

  struct residue *residue_ptr;

  FILE *file_ptr;

  /* Initialise variables */

  /* Form the name of the pymol script file */
  strcpy(file_name,"ligand.pml");
  sprintf(png_name,"ligand_%2.2d.png",ligand_ptr->type);

  /* Open output file */
  file_ptr = open_output_file(file_name,TRUE);

  /* Write out the header records */
  fprintf(file_ptr,"#\n");
  fprintf(file_ptr,"# PDBsum PyMol script for ligand %d\n",
          ligand_ptr->type);
  fprintf(file_ptr,"#\n");
  fprintf(file_ptr,"load %s\n",pdb_file);
  fprintf(file_ptr,"bg_color white\n");
  fprintf(file_ptr,"set depth_cue=0\n");
  fprintf(file_ptr,"set ray_trace_fog=0\n");
  fprintf(file_ptr,"set cartoon_ring_mode,1\n");
  fprintf(file_ptr,"hide all\n");
  fprintf(file_ptr,"set auto_color, 0\n");

  /* Define all the RGB colours */
  for (colour = 0; colour < NCOLOURS; colour++)
    {
      if (!strcmp(COLOUR_LIST[colour],"grey"))
        fprintf(file_ptr,"set_color %s, [ %6.3f, %6.3f, %6.3f ]\n",
                COLOUR_LIST[colour],RGB_TRIPLET[colour][0],
                RGB_TRIPLET[colour][1],RGB_TRIPLET[colour][2]);
    }

  /* Define the ligand */
  fprintf(file_ptr,"\n");
  fprintf(file_ptr,"# Ligand %d. %s\n",ligand_ptr->number,
          ligand_ptr->residue_sequence);

  /* Get this ligand's first residue */
  residue_ptr = ligand_ptr->first_residue_ptr;
  nresid = ligand_ptr->nresid;
  iresid = 0;
  res_num1[0] = res_num2[0] = '\0';

  /* Loop through this ligand's residues, to find start and
     end residues */
  while (residue_ptr != NULL && iresid < nresid)
    {
      /* If this is the first residue, save number */
      if (res_num1[0] == '\0')
        strcpy(res_num1,residue_ptr->res_num);

      /* Save number as last residue */
      strcpy(res_num2,residue_ptr->res_num);

      /* Store chain id */
      chain = residue_ptr->chain;

      /* Get the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;
    }

  /* If have a residue range, write out definition */
  if (res_num1[0] != '\0' && res_num2[0] != '\0')
    {
      /* Remove insertion codes as PyMol doesn't understand them */
      res_num1[4] = '\0';
      res_num2[4] = '\0';

      /* Define the residue range */
      sprintf(residue_range,"%s-%s",res_num1,res_num2);

      /* Check that residue numbers are in ascending order */
      res_num1[4] = '\0';
      res_num2[4] = '\0';
      rnum1 = atoi(res_num1);
      rnum2 = atoi(res_num2);
      if (rnum1 > rnum2)
        sprintf(residue_range,"%s-%s",res_num2,res_num1);

      /* Remove blanks from the range */
      valid = remove_residue_blanks(residue_range);

      /* If chain not blank, then add chain */
      if (chain != ' ')
        {
          sprintf(tmp,"%s and chain %c",residue_range,
                  chain);
          strcpy(residue_range,tmp);
        }

      /* Write out this ligand definition */
      if (valid == TRUE)
        fprintf(file_ptr,"create ligand, (resi %s)\n",
                residue_range);
    }

  /* Show ligand as sticks */
  fprintf(file_ptr,"show sticks, ligand\n");

  /* Colour by atom type */
  fprintf(file_ptr,"util.cbag ligand\n");
  fprintf(file_ptr,"color grey, ligand and elem C\n");

  /* Zoom in on final picture */
  fprintf(file_ptr,"zoom all, complete=1, buffer=2\n");

  /* Render the ligand */
  fprintf(file_ptr,"ray 240, 240\n");
  fprintf(file_ptr,"viewport 240, 240\n");

  /* Write out command to generate png file */
  fprintf(file_ptr,"png %s\n",png_name);

  /* Close the file */
  fclose(file_ptr);

  /* Run PyMOL to render the ligand */
  run_pymol("ligand",pymol_exe,dos);
}
/***********************************************************************

rotate_structure  -  Apply the given rotation to whole structure

***********************************************************************/

void rotate_structure(struct atom *first_atom_ptr,int natoms,
                      double matrix[3][3])
{
  int iatom, icoord;

  double x, y, z;

  struct atom *atom_ptr;

  /* Get pointer to the first atom */
  atom_ptr = first_atom_ptr;
  iatom = 0;

  /* Loop through all the atoms */
  while (atom_ptr != NULL && iatom < natoms)
    {
      /* Retrieve the coordinates for this atom */
      x = atom_ptr->coords[0];
      y = atom_ptr->coords[1];
      z = atom_ptr->coords[2];

      /* Apply the transformation */
      for (icoord = 0; icoord < 3; icoord++)
        {
          /* Calculate transformed coordinates */
          atom_ptr->coords[icoord]
            = ((x * matrix[0][icoord])
               + (y * matrix[1][icoord])
               + (z * matrix[2][icoord]));

          /* Store as new view coords */
          atom_ptr->view_coords[icoord] = atom_ptr->coords[icoord];
        }

      /* Get pointer to the next atom */
      atom_ptr = atom_ptr->next_atom_ptr;
      iatom++;
    }
}
/***********************************************************************

get_max_min  -  Get the maximum and minimum coords of the current
                molecule

***********************************************************************/

void get_max_min(struct atom *first_atom_ptr,int natoms,
                 double *coord_min,double *coord_max)
{
  int iatom, icoord;

  struct atom *atom_ptr;

  /* Get pointer to the first atom */
  atom_ptr = first_atom_ptr;
  iatom = 0;

  /* If have valid atom, store its coords as the current maximum and
     minimum */
  if (atom_ptr != NULL)
    {
      for (icoord = 0; icoord < 3; icoord++)
        {
          coord_min[icoord] = atom_ptr->view_coords[icoord];
          coord_max[icoord] = atom_ptr->view_coords[icoord];
        }
    }

  /* Otherwise, set all to zero */
  else
    {
      for (icoord = 0; icoord < 3; icoord++)
        {
          coord_min[icoord] = 0.0;
          coord_max[icoord] = 0.0;
        }
    }

  /* Loop through all the atoms */
  while (atom_ptr != NULL && iatom < natoms)
    {
      /* Update coordinate limits */
      for (icoord = 0; icoord < 3; icoord++)
        {
          if (atom_ptr->view_coords[icoord] < coord_min[icoord])
            coord_min[icoord] = atom_ptr->view_coords[icoord];
          if (atom_ptr->view_coords[icoord] > coord_max[icoord])
            coord_max[icoord] = atom_ptr->view_coords[icoord];
        }

      /* Get pointer to the next atom */
      atom_ptr = atom_ptr->next_atom_ptr;
      iatom++;
    }
}
/***********************************************************************

get_max_min2  -  Get the maximum and minimum coords of the two sets of
                 atoms

***********************************************************************/

void get_max_min2(struct atom *first_atom1_ptr,int natoms1,
                  struct atom *first_atom2_ptr,int natoms2,
                  double *coord_min,double *coord_max)
{
  int iatom, icoord, loop, natoms;

  struct atom *atom_ptr;

  /* Loop over the two sets of atoms */
  for (loop = 0; loop < 2; loop++)
    {
      /* Get pointer to the first atom */
      if (loop == 0)
        {
          atom_ptr = first_atom1_ptr;
          natoms = natoms1;
        }
      else
        {
          atom_ptr = first_atom2_ptr;
          natoms = natoms2;
        }
      iatom = 0;

      /* If have valid atom, store its coords as the current maximum and
         minimum */
      if (atom_ptr != NULL)
        {
          for (icoord = 0; icoord < 3; icoord++)
            {
              coord_min[icoord] = atom_ptr->view_coords[icoord];
              coord_max[icoord] = atom_ptr->view_coords[icoord];
            }
        }

      /* Otherwise, set all to zero */
      else
        {
          for (icoord = 0; icoord < 3; icoord++)
            {
              coord_min[icoord] = 0.0;
              coord_max[icoord] = 0.0;
            }
        }

      /* Loop through all the atoms */
      while (atom_ptr != NULL && iatom < natoms)
        {
          /* Update coordinate limits */
          for (icoord = 0; icoord < 3; icoord++)
            {
              if (atom_ptr->view_coords[icoord] < coord_min[icoord])
                coord_min[icoord] = atom_ptr->view_coords[icoord];
              if (atom_ptr->view_coords[icoord] > coord_max[icoord])
                coord_max[icoord] = atom_ptr->view_coords[icoord];
            }

          /* Get pointer to the next atom */
          atom_ptr = atom_ptr->next_atom_ptr;
          iatom++;
        }
    }
}
/***********************************************************************

get_chain_col  -  Get the given chain's colour

***********************************************************************/

int get_chain_col(char chain,double *rgb,double *rgb_hel,
                  double *rgb_str,int split)
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

  /* If splitting by secondary structure, compute the corresponding
     helix and strand colours */
  if (split == TRUE)
    {
      /* Calculate helix and strand colours */
      for (i = 0; i < 3; i++)
        {
          /* Save chain RGB values as helix colours */
          rgb_hel[i] = rgb[i];

          /* Compute coil colour to be darker */
          rgb[i] = 0.6 * rgb_hel[i];

          /* Compute strand colour to be lighter */
          rgb_str[i] = 0.4 + rgb_hel[i];
          if (rgb_str[i] > 1.0)
            rgb_str[i] = 1.0;
        }
    }

  /* Otherwise return same RGB values to strand and helix as for chain */
  else
    {
      /* Copy the RGB values across */
      for (i = 0; i < 3; i++)
        rgb_hel[i] = rgb_str[i] = rgb[i];
    }

  /* Return the colour number */
  return(colour);
}
/***********************************************************************

get_chain_colname  -  Get the given chain's colour name

***********************************************************************/

void get_chain_colname(char chain,char *cname,char *ccoil,char *cmajor,
                       char *cminor)
{
  int colour, ichain, ipos, nchains;

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

  /* Get the colour name */
  strcpy(cname,CHAIN_COLOUR[colour]);
  strcpy(ccoil,CHAIN_COLOUR[colour]);
  strcat(ccoil,"_coil");
  strcpy(cmajor,CHAIN_COLOUR[colour]);
  strcat(cmajor,"_major");
  strcpy(cminor,CHAIN_COLOUR[colour]);
  strcat(cminor,"_minor");
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

reverse_transform  -  Transform given point back into the correct
                      position in the structure

***********************************************************************/

void reverse_transform(double x,double y,double z,double *newx,
                       double *newy,double *newz,
                       double transformation[3][3],double cofg[3])
{
  /* Apply the reverse transformation to the current point */
  *newx = x * transformation[0][0] + y * transformation[0][1]
    + z * transformation[0][2] + cofg[0];
  *newy = x * transformation[1][0] + y * transformation[1][1]
    + z * transformation[1][2] + cofg[1];
  *newz = x * transformation[2][0] + y * transformation[2][1]
    + z * transformation[2][2] + cofg[2];
}
/***********************************************************************

generate_helix  -  Render helix as a single, solid flat-ended
                   cylinder

***********************************************************************/

void generate_helix(FILE *file_ptr,double length1,double length2,
                    double transformation[3][3],double cofg[3],
                    double scale_factor,double *rgb)
{
  double radius, x1, x2, y1, y2, z1, z2;
  double newx1, newx2, newy1, newy2, newz1, newz2;

  /* Initialise variables */
  radius = scale_factor * HELIX_RADIUS;

  /* Define the helix end-points on the x-axis */
  x1 = - length1;
  y1 = 0.0;
  z1 = 0.0;
  x2 = length2;
  y2 = 0.0;
  z2 = 0.0;

  /* Transform them back into the correct position in the structure */
  reverse_transform(x1,y1,z1,&newx1,&newy1,&newz1,transformation,cofg);
  reverse_transform(x2,y2,z2,&newx2,&newy2,&newz2,transformation,cofg);

  /* Apply scale factor to the helix end-points */
  newx1 = scale_factor * newx1;
  newy1 = scale_factor * newy1;
  newz1 = scale_factor * newz1;
  newx2 = scale_factor * newx2;
  newy2 = scale_factor * newy2;
  newz2 = scale_factor * newz2;

  /* Write out to Raster3D file as a flat-ended cylinder */
  fprintf(file_ptr,"5\n");
  fprintf(file_ptr,"%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n",
          newx1,newy1,newz1,radius,newx2,newy2,newz2);
  fprintf(file_ptr,"0.0 %6.3f %6.3f %6.3f\n",rgb[0],rgb[1],rgb[2]);
}
/***********************************************************************

store_vertex  -  Store the given coords as the nhext strand vertex

***********************************************************************/

int store_vertex(double vertex[NSTRAND_VERTICES][3],int nvertices,
                 double x,double y,double z)
{
  /* Store current coords */
  vertex[nvertices][0] = x;
  vertex[nvertices][1] = y;
  vertex[nvertices][2] = z;

  /* Increment vertex count */
  nvertices++;

  /* Return current number of vertices */
  return(nvertices);
}
/***********************************************************************

write_polygon  -  Write out given polygon to Raster3D file

***********************************************************************/

void write_polygon(FILE *file_ptr,double vertex1[3],double vertex2[3],
                   double vertex3[3],double transformation[3][3],
                   double cofg[3],double scale_factor,double *rgb)
{
  double newx1, newx2, newx3, newy1, newy2, newy3, newz1, newz2, newz3;

  /* Apply reverse transformation to each of the polygon's points */
  reverse_transform(vertex1[0],vertex1[1],vertex1[2],&newx1,&newy1,
                    &newz1,transformation,cofg);
  reverse_transform(vertex2[0],vertex2[1],vertex2[2],&newx2,&newy2,
                    &newz2,transformation,cofg);
  reverse_transform(vertex3[0],vertex3[1],vertex3[2],&newx3,&newy3,
                    &newz3,transformation,cofg);

  /* Apply scale factor to each point */
  newx1 = scale_factor * newx1;
  newy1 = scale_factor * newy1;
  newz1 = scale_factor * newz1;
  newx2 = scale_factor * newx2;
  newy2 = scale_factor * newy2;
  newz2 = scale_factor * newz2;
  newx3 = scale_factor * newx3;
  newy3 = scale_factor * newy3;
  newz3 = scale_factor * newz3;

  /* Write out to Raster3D file as a triangle */
  fprintf(file_ptr,"1\n");
  fprintf(file_ptr,
          "%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n",
          newx1,newy1,newz1,newx2,newy2,newz2,newx3,newy3,newz3);
  fprintf(file_ptr,"%6.3f %6.3f %6.3f\n",rgb[0],rgb[1],rgb[2]);
}
/***********************************************************************

generate_strand  -  Render strand as a large arrow

***********************************************************************/

void generate_strand(FILE *file_ptr,double length1,double length2,
                     double transformation[3][3],double cofg[3],
                     double scale_factor,double *rgb,int reverse)
{
  int nvertices;

  double x1, x2, x3, y1, y2, y3, y4, y5, z1, z2;
  double vertex[NSTRAND_VERTICES][3];

  /* Initialise variables */
  nvertices = 0;

  /* Calculate where the starts and ends of the strand and its head are */
  x1 = - length1;
  x3 =   length2;
  x2 = x3 - STRAND_HEAD_LEN;
  if (x2 < 0.0)
    x2 = 0.0;
  if (reverse == TRUE)
    {
      x1 =   length2;
      x3 = - length1;
      x2 = x3 + STRAND_HEAD_LEN;
      if (x2 > 0.0)
        x2 = 0.0;
    }

  /* Calculate the strand width and thickness */
  y1 = STRAND_HEAD_WIDTH / 2.0;
  y2 = STRAND_WIDTH / 2.0;
  y3 = - y2;
  y4 = - y1;
  y5 = 0.0;
  z1 = STRAND_THICKNESS / 2.0;
  z2 = - z1;
  if (reverse == TRUE)
    {
      z1 = - z1;
      z2 = - z2;
    }

  /* Store the vertices making up the end face of the strand body */
  nvertices = store_vertex(vertex,nvertices,x1,y2,z2);
  nvertices = store_vertex(vertex,nvertices,x1,y3,z2);
  nvertices = store_vertex(vertex,nvertices,x1,y3,z1);
  nvertices = store_vertex(vertex,nvertices,x1,y2,z1);

  /* Store the vertices for the body-end */
  nvertices = store_vertex(vertex,nvertices,x2,y2,z2);
  nvertices = store_vertex(vertex,nvertices,x2,y3,z2);
  nvertices = store_vertex(vertex,nvertices,x2,y3,z1);
  nvertices = store_vertex(vertex,nvertices,x2,y2,z1);

  /* Polygon 1: Form the end-face polygon */
  write_polygon(file_ptr,vertex[0],vertex[1],vertex[2],transformation,
                cofg,scale_factor,rgb);
  write_polygon(file_ptr,vertex[0],vertex[2],vertex[3],transformation,
                cofg,scale_factor,rgb);

  /* Polygon 2: Side face 1 */
  write_polygon(file_ptr,vertex[2],vertex[1],vertex[5],transformation,
                cofg,scale_factor,rgb);
  write_polygon(file_ptr,vertex[2],vertex[5],vertex[6],transformation,
                cofg,scale_factor,rgb);

  /* Polygon 3: Top face */
  write_polygon(file_ptr,vertex[3],vertex[2],vertex[6],transformation,
                cofg,scale_factor,rgb);
  write_polygon(file_ptr,vertex[3],vertex[6],vertex[7],transformation,
                cofg,scale_factor,rgb);

  /* Polygon 4: Side face 2 */
  write_polygon(file_ptr,vertex[0],vertex[3],vertex[7],transformation,
                cofg,scale_factor,rgb);
  write_polygon(file_ptr,vertex[0],vertex[7],vertex[4],transformation,
                cofg,scale_factor,rgb);

  /* Polygon 5: Bottom face */
  write_polygon(file_ptr,vertex[1],vertex[0],vertex[4],transformation,
                cofg,scale_factor,rgb);
  write_polygon(file_ptr,vertex[1],vertex[4],vertex[5],transformation,
                cofg,scale_factor,rgb);

  /* Store the vertices making up the strand head */
  nvertices = store_vertex(vertex,nvertices,x2,y1,z2);
  nvertices = store_vertex(vertex,nvertices,x2,y4,z2);
  nvertices = store_vertex(vertex,nvertices,x2,y4,z1);
  nvertices = store_vertex(vertex,nvertices,x2,y1,z1);
  nvertices = store_vertex(vertex,nvertices,x3,y5,z2);
  nvertices = store_vertex(vertex,nvertices,x3,y5,z1);

  /* Form all the polygons making up the strand head */

  /* Polygon 6: Top triangle */
  write_polygon(file_ptr,vertex[10],vertex[13],vertex[11],transformation,
                cofg,scale_factor,rgb);

  /* Polygon 7: Wing face 1 */
  write_polygon(file_ptr,vertex[5],vertex[9],vertex[10],transformation,
                cofg,scale_factor,rgb);
  write_polygon(file_ptr,vertex[5],vertex[10],vertex[6],transformation,
                cofg,scale_factor,rgb);

  /* Polygon 8: Wing face 2 */
  write_polygon(file_ptr,vertex[8],vertex[4],vertex[7],transformation,
                cofg,scale_factor,rgb);
  write_polygon(file_ptr,vertex[8],vertex[7],vertex[11],transformation,
                cofg,scale_factor,rgb);

  /* Polygon 9: Bottom triangle */
  write_polygon(file_ptr,vertex[8],vertex[12],vertex[9],transformation,
                cofg,scale_factor,rgb);

  /* Polygon 10: Triangle side 1 */
  write_polygon(file_ptr,vertex[10],vertex[9],vertex[12],transformation,
                cofg,scale_factor,rgb);
  write_polygon(file_ptr,vertex[10],vertex[12],vertex[13],transformation,
                cofg,scale_factor,rgb);

  /* Polygon 11: Triangle side 2 */
  write_polygon(file_ptr,vertex[11],vertex[13],vertex[12],transformation,
                cofg,scale_factor,rgb);
  write_polygon(file_ptr,vertex[11],vertex[12],vertex[8],transformation,
                cofg,scale_factor,rgb);
}
/***********************************************************************

modify_calphas  -  Modify all the C-alpha positions belonging to the
                   current helix/strand so that they lie within the
                   helix/strand

***********************************************************************/

void modify_calphas(struct residue *first_residue_ptr,int nresid,
                    double length1,double length2,int reverse,
                    double transformation[3][3],double cofg[3])
{

  int iatom, iresid, natoms;
  int found;

  double end1, end2, fraction;
  double x, x1, x2, y, y1, y2, z, z1, z2;

  struct atom *atom_ptr;
  struct residue *residue_ptr;

  /* Initialise variables */

  /* Get the two ends of the helix/strand, depending on its direction */
  if (reverse == TRUE)
    {
      end1 = length2;
      end2 = - length1;
    }
  else
    {
      end1 = - length1;
      end2 = length2;
    }

  /* Form the coords of the untransformed end-points of this helix/strand
     and apply the reverse transformation to them */
  x = end1;
  y = 0.0;
  z = 0.0;
  reverse_transform(x,y,z,&x1,&y1,&z1,transformation,cofg);
  x = end2;
  y = 0.0;
  z = 0.0;
  reverse_transform(x,y,z,&x2,&y2,&z2,transformation,cofg);

  /* Get the first residue of this secondary structure element */
  residue_ptr = first_residue_ptr;
  iresid = 0;

  /* Loop through the residues to locate all the C-alpha atoms */
  while (residue_ptr != NULL && iresid < nresid)
    {
      /* Get the first of this residue's atoms */
      atom_ptr = residue_ptr->first_atom_ptr;
      iatom = 0;
      natoms = residue_ptr->natoms;
      found = FALSE;
  
      /* Loop through all the residue's atoms to find the C-alphas */
      while (atom_ptr != NULL && iatom < natoms && found == FALSE)
        {
          /* If this is a C-alpha atom, modify it coordinates */
          if (!strcmp(atom_ptr->atom_name," CA "))
            {
              /* Calculate the modified coords of this atom */
              fraction = (float) iresid / (float) (nresid - 1);
              x = x1 + fraction * (x2 - x1);
              y = y1 + fraction * (y2 - y1);
              z = z1 + fraction * (z2 - z1);

              /* Store these coords */
              atom_ptr->view_coords[0] = x;
              atom_ptr->view_coords[1] = y;
              atom_ptr->view_coords[2] = z;

              /* Update flag */
              found = TRUE;
            }

          /* Get pointer to next atom in linked-list */
          atom_ptr = atom_ptr->next_atom_ptr;
          iatom++;
        }

      /* Get the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;
    }
}
/***********************************************************************

render_sec_struc  -  Render given secondary structure

***********************************************************************/

void render_sec_struc(FILE *file_ptr,struct residue *first_residue_ptr,
                      int nresid,double scale_factor,double *rgb,
                      int sec_struc,int nsse)
{
  char dummy_name[NMETALS][3];

  int dummy_count[NMETALS], ndummy_types;
  int i, iatom, iresid, natoms, ncreated, ndummy;
  int reverse;

  double coord_min[3], coord_max[3], length1, length2;
  double transformation[3][3], cofg[3];

  struct residue dummy_residue_ptr;

  struct atom *atom_ptr, *new_atom_ptr;
  struct atom *first_ca_atom_ptr, *last_ca_atom_ptr;
  struct atom *first_atom_ptr, *last_atom_ptr;
  struct residue *residue_ptr;

  /* Initialise variables */
  ndummy_types = 0;
  ncreated = 0;
  reverse = FALSE;
  for (i = 0; i < NMETALS; i++)
    {
      dummy_name[i][0] = '\0';
      dummy_count[i] = 0;
    }

  /* Write comment to Raster3D file */
  if (sec_struc == HELIX)
    fprintf(file_ptr,"#   Helix %d\n",nsse);
  else
    fprintf(file_ptr,"#   Strand %d\n",nsse);

  /* Initialise atom pointers */
  atom_ptr = NULL;
  first_atom_ptr = NULL;
  last_atom_ptr = NULL;
  first_ca_atom_ptr = NULL;
  last_ca_atom_ptr = NULL;

  /* Initialise dummy residue record */
  dummy_residue_ptr.natoms = 0;
  dummy_residue_ptr.first_atom_ptr = NULL;
  dummy_residue_ptr.nmain_chain = 0;

  /* Get the first residue of this secondary structure element */
  residue_ptr = first_residue_ptr;
  iresid = 0;

  /* Loop through the residues to extract all the main-chain atoms and
     create copy atoms of each */
  while (residue_ptr != NULL && iresid < nresid)
    {
      /* Get the first of this residue's atoms */
      atom_ptr = residue_ptr->first_atom_ptr;
      iatom = 0;
      natoms = residue_ptr->natoms;
  
      /* Loop through all the residue's atoms to find the mainchain
         ones */
      while (atom_ptr != NULL && iatom < natoms)
        {
          /* If this is a mainchain atom, create a copy */
          if (!strcmp(atom_ptr->atom_name," N  ") ||
              !strcmp(atom_ptr->atom_name," CA ") ||
              !strcmp(atom_ptr->atom_name," C  "))
            {
              /* Create a new atom */
              new_atom_ptr
                = create_atom_record(&first_atom_ptr,&last_atom_ptr,
                                     &dummy_residue_ptr,
                                     atom_ptr->atom_number,
                                     atom_ptr->atom_name,
                                     atom_ptr->coords[0],
                                     atom_ptr->coords[1],
                                     atom_ptr->coords[2],
                                     atom_ptr->bvalue,
                                     atom_ptr->occupancy,dummy_name,
                                     dummy_count,&ndummy_types,&ndummy,
                                     ATOM_REC);

              /* Increment count of atoms created */
              ncreated++;

              /* If this is a C-alpha atom, store first and last */
              if (!strcmp(new_atom_ptr->atom_name," CA "))
                {
                  /* If first, store as first */
                  if (first_ca_atom_ptr == NULL)
                    first_ca_atom_ptr = new_atom_ptr;

                  /* Store as last */
                  last_ca_atom_ptr = new_atom_ptr;
                }
            }

          /* Get pointer to next atom in linked-list */
          atom_ptr = atom_ptr->next_atom_ptr;
          iatom++;
        }

      /* Get the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;
    }

  /* If have sufficient info to generate the secondary structure
     element, then do so */
  if (ncreated > 3 && first_ca_atom_ptr != NULL &&
      last_ca_atom_ptr != NULL && first_ca_atom_ptr != last_ca_atom_ptr)
    {
      /* Align helix or strand parallel to the x-axis */
      align_molecule(first_atom_ptr,ncreated,coord_min,coord_max,
                     transformation,cofg);

      /* Determine in which direction the polypeptide chain lies in
         the transformed helix/strand */
      length1 = fabs(first_ca_atom_ptr->view_coords[0]);
      length2 = fabs(last_ca_atom_ptr->view_coords[0]);
      if (last_ca_atom_ptr->view_coords[0]
          < first_ca_atom_ptr->view_coords[0])
        reverse = TRUE;
      
      /* Generate the standard helix/strand using the given length */
      if (sec_struc == HELIX)
        generate_helix(file_ptr,length1,length2,transformation,cofg,
                       scale_factor,rgb);
      else
        generate_strand(file_ptr,length1,length2,transformation,cofg,
                        scale_factor,rgb,reverse);

      /* Modify the coordinates of all C-alpha atoms in this secondary
         structure element such that they lie within the helix/strand */
      modify_calphas(first_residue_ptr,nresid,length1,length2,reverse,
                     transformation,cofg);
    }

  /* Free up memory taken by atom records created here */
  free_atom_memory(first_atom_ptr);
}
/***********************************************************************

find_sec_struc  -  Locate each instance of the given secondary structure

***********************************************************************/

void find_sec_struc(FILE *file_ptr,struct chain *chain_ptr,
                    double scale_factor,double *rgb,int sec_struc)
{
  int iresid, nresid, nsse, nstresid;

  struct residue *residue_ptr, *start_residue_ptr;

  /* Initialise variables */
  start_residue_ptr = NULL;
  nsse = 0;
  nstresid = 0;

  /* Write comment to Raster3D file */
  if (sec_struc == HELIX)
    fprintf(file_ptr,"#  Helices\n");
  else
    fprintf(file_ptr,"#  Strands\n");

  /* Get this chain's first residue */
  residue_ptr = chain_ptr->first_residue_ptr;
  nresid = chain_ptr->nresid;
  iresid = 0;

  /* Loop through this chain's residues, updating all the counts */
  while (residue_ptr != NULL && iresid < nresid)
    {
      /* Check if this residue belongs to the required secondary
         structure */
      if (residue_ptr->sec_struc == sec_struc)
        {
          /* If this is the first residue, store its pointer */
          if (nstresid == 0)
            start_residue_ptr = residue_ptr;

          /* Increment residue count */
          nstresid++;
        }

      /* If not the right secondary structure, or this is the last
         residue, see if need to process previous secondary structure */
      if (residue_ptr->sec_struc != sec_struc || iresid == nresid - 1)
        {
          /* If have a secondary structure element to process, then
             do so */
          if (nstresid > 0)
            {
              /* Increment count of secondary structure elements done */
              nsse++;

              /* Process current secondary structure */
              render_sec_struc(file_ptr,start_residue_ptr,nstresid,
                               scale_factor,rgb,sec_struc,nsse);

              /* Re-initialise for next one */
              start_residue_ptr = NULL;
              nstresid = 0;
            }
        }

      /* Get the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;
    }
}
/***********************************************************************

render_backbone  -  Render the C-alpha backbone

***********************************************************************/

void render_backbone(FILE *file_ptr,struct chain *chain_ptr,
                     double scale_factor,double *rgb)
{
  char backbone_atom[5];

  int iatom, iresid, natoms, nresid, type;
  int last_sec_struc;
  int found, wanted;

  double cutoff_sqrd, dist_sqrd, x1, x2, y1, y2, z1, z2;

  struct atom *atom_ptr, *backbone_atom_ptr, *prev_backbone_atom_ptr;
  struct residue *residue_ptr, *next_residue_ptr;

  /* Initialise variables */
  last_sec_struc = OTHER;
  prev_backbone_atom_ptr = NULL;
  type = chain_ptr->type;

  /* Write comment to Raster3D file */
  if (type == PROTEIN)
    fprintf(file_ptr,"#  Coil\n");
  else if (type == CALPHA_ONLY)
    fprintf(file_ptr,"#  C-alpha backbone\n");
  else if (type == DNA)
    fprintf(file_ptr,"#  DNA phosphate backbone\n");

  /* Get name of backbone atom that we're looking for and the cutoff
     distance */
  if (type == DNA)
    {
      strcpy(backbone_atom," P  ");
      cutoff_sqrd = P_P_CUTOFF_SQRD;
    }
  else
    {
      strcpy(backbone_atom," CA ");
      cutoff_sqrd = CA_CA_CUTOFF_SQRD;
    }

  /* Get this chain's first residue */
  residue_ptr = chain_ptr->first_residue_ptr;
  nresid = chain_ptr->nresid;
  iresid = 0;

  /* Loop through this chain's residues to find the coil regions */
  while (residue_ptr != NULL && iresid < nresid)
    {
      /* Assume residue is not wanted */
      wanted = FALSE;

      /* If residue is in a coil region then it is wanted */
      if (residue_ptr->sec_struc == COIL)
        wanted = TRUE;

      /* Otherwise, determine whether it is a secondary structure start
         or end */
      else
        {
          /* If last secondary structure was a coil and this is a
             helix or strand residue, then it is wanted */
          if (last_sec_struc == COIL)
            wanted = TRUE;

          /* Otherwise, determine whether next residue is a coil
             residue */
          else
            {
              /* Get the next residue */
              next_residue_ptr = residue_ptr->next_residue_ptr;

              /* If it exists, check its secondary structure */
              if (next_residue_ptr != NULL &&
                  next_residue_ptr->sec_struc == COIL)
                wanted = TRUE;
            }
        }

      /* Save the current secondary structure */
      last_sec_struc = residue_ptr->sec_struc;

      /* If residue wanted, draw backbone to previous residue, if 
         that exists */
      if (wanted == TRUE)
        {
          /* Initialise pointer to this residue's C-alpha atom */
          backbone_atom_ptr = NULL;
          found = FALSE;

          /* Get the first of this residue's atoms */
          atom_ptr = residue_ptr->first_atom_ptr;
          iatom = 0;
          natoms = residue_ptr->natoms;
  
          /* Save pointer to the first atom, just in case we can't find
             the c-alpha itself */
          backbone_atom_ptr = atom_ptr;

          /* Loop through all the residue's atoms to find the C-alpha
             or phosphate atom */
          while (atom_ptr != NULL && iatom < natoms && found == FALSE)
            {
              /* If this is C-alpha atom, store its pointer */
              if (!strcmp(atom_ptr->atom_name,backbone_atom))
                {
                  /* Store pointer and set flag */
                  backbone_atom_ptr = atom_ptr;
                  found = TRUE;
                }

              /* Get pointer to next atom in linked-list */
              atom_ptr = atom_ptr->next_atom_ptr;
              iatom++;
            }

          /* If have previous and current residue's C-alphas, draw a
             backbone bond between them */
          if (backbone_atom_ptr != NULL && prev_backbone_atom_ptr != NULL)
            {
              /* Check that these two atoms aren't too far apart (eg due
                 to a chain break bat this point) */

              /* Get both atoms' coordinates */
              x1 = prev_backbone_atom_ptr->coords[0];
              y1 = prev_backbone_atom_ptr->coords[1];
              z1 = prev_backbone_atom_ptr->coords[2];
              x2 = backbone_atom_ptr->coords[0];
              y2 = backbone_atom_ptr->coords[1];
              z2 = backbone_atom_ptr->coords[2];

              /* Calculate squared distance between these two atoms */
              dist_sqrd
                = (x2 - x1) * (x2 - x1)
                + (y2 - y1) * (y2 - y1)
                + (z2 - z1) * (z2 - z1);
              
              /* If distance is below the cutoff, then can join the
                 atoms */
              if (dist_sqrd <= cutoff_sqrd)
                write_raster3d_bond(file_ptr,backbone_atom_ptr,
                                    prev_backbone_atom_ptr,scale_factor,
                                    BACKBONE_RADIUS,FALSE,rgb,FALSE,
                                    FALSE);
            }

          /* Store current C-alpha atom and its coordinates */
          prev_backbone_atom_ptr = backbone_atom_ptr;
        }

      /* Otherwise, reinitialise C-alpha pointer */
      else
        {
          prev_backbone_atom_ptr = NULL;
        }

      /* Get the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;
    }
}
/***********************************************************************

write_line  -  Write out given line as a bond in the Raster3D file

***********************************************************************/

void write_line(FILE *file_ptr,double vertex1[3],double vertex2[3],
                double scale_factor,double *rgb,double radius)
{
  double x1, x2, y1, y2, z1, z2;

  /* Apply scale factor to each point */
  x1 = scale_factor * vertex1[0];
  y1 = scale_factor * vertex1[1];
  z1 = scale_factor * vertex1[2];
  x2 = scale_factor * vertex2[0];
  y2 = scale_factor * vertex2[1];
  z2 = scale_factor * vertex2[2];

  /* Get the line radius */
  radius = scale_factor * radius;

  /* Write out to Raster3D file as a line */
  fprintf(file_ptr,"3\n");
  fprintf(file_ptr,"%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n",
          x1,y1,z1,radius,x2,y2,z2);
  fprintf(file_ptr,"0.0 %6.3f %6.3f %6.3f\n",rgb[0],rgb[1],rgb[2]);
}
/***********************************************************************

label_nterminus  -  Label the N-terminus of the given protein chain

***********************************************************************/

void label_nterminus(FILE *file_ptr,struct chain *chain_ptr,
                     double scale_factor)
{
  int iatom, natoms, nvertices;
  int found;

  double centre[3], rgb[3], x, y, z;
  double vertex[NSTRAND_VERTICES][3];

  struct atom *atom_ptr, *ca_atom_ptr;
  struct residue *residue_ptr;

  /* Initialise variables */
  ca_atom_ptr = NULL;
  found = FALSE;
  nvertices = 0;

  /* Get rgb value for colour of the "N" */
  get_colour_rgb("black",rgb);

  /* Get this chain's first residue */
  residue_ptr = chain_ptr->first_residue_ptr;

  /* Get the first of this residue's atoms */
  atom_ptr = residue_ptr->first_atom_ptr;
  iatom = 0;
  natoms = residue_ptr->natoms;
  
  /* Save pointer to the first atom, just in case we can't find the
     c-alpha itself */
  ca_atom_ptr = atom_ptr;

  /* Loop through all the residue's atoms to find the C-alpha atom */
  while (atom_ptr != NULL && iatom < natoms && found == FALSE)
    {
      /* If this is C-alpha atom, store its pointer */
      if (!strcmp(atom_ptr->atom_name," CA "))
        {
          /* Store pointer and set flag */
          ca_atom_ptr = atom_ptr;
          found = TRUE;
        }

      /* Get pointer to next atom in linked-list */
      atom_ptr = atom_ptr->next_atom_ptr;
      iatom++;
    }

  /* If have the C-alpha atom, generate the N to mark this end of the
     protein chain */
  if (ca_atom_ptr != NULL)
    {
      /* Write comment to Raster3D file */
      fprintf(file_ptr,"# N-terminus of chain %c\n",chain_ptr->chain_id);

      /* Get this atom's coordinates */
      x = ca_atom_ptr->coords[0];
      y = ca_atom_ptr->coords[1];
      z = ca_atom_ptr->coords[2];

      /* Form the coordinates for the letter N */
      centre[0] = x + 0.6;
      centre[1] = y + 0.6;
      centre[2] = z + 0.6;

      /* Create and store the 4 vertices defining the letter N */

      /* Vertex 1 */
      x = centre[0] - N_WIDTH / 2.0;
      y = centre[1] - N_HEIGHT / 2.0;
      z = centre[2];
      nvertices = store_vertex(vertex,nvertices,x,y,z);

      /* Vertex 2 */
      y = y + N_HEIGHT;
      nvertices = store_vertex(vertex,nvertices,x,y,z);

      /* Vertex 3 */
      x = x + N_WIDTH;
      nvertices = store_vertex(vertex,nvertices,x,y,z);

      /* Vertex 4 */
      y = y - N_HEIGHT;
      nvertices = store_vertex(vertex,nvertices,x,y,z);

      /* Write out the 3 "bonds" making up the "N" */

      /* Line 1 */
      write_line(file_ptr,vertex[0],vertex[1],scale_factor,rgb,N_RADIUS);

      /* Line 2 */
      write_line(file_ptr,vertex[1],vertex[3],scale_factor,rgb,N_RADIUS);

      /* Line 3 */
      write_line(file_ptr,vertex[2],vertex[3],scale_factor,rgb,N_RADIUS);
    }
}
/***********************************************************************

material_properties  -  Write out the default material properties

***********************************************************************/

void material_properties(FILE *file_ptr,int transparency)
{
  static char *material_dfn[3] =
    { 
      "25 0.25 1.0 1.0 1.0 0.0 0 0 0 0", 
      "25 0.25 1.0 1.0 1.0 0.6 0 0 0 0", 
      "25 0.25 1.0 1.0 1.0 0.9 0 0 0 0"
    };

  if (transparency == NONE)
    {
      fprintf(file_ptr,"#8\n");
      fprintf(file_ptr,"#-\n");
    }
  else
    {
      fprintf(file_ptr,"8\n");
      fprintf(file_ptr,"%s\n",material_dfn[transparency]);
    }
}
/***********************************************************************

render_protein_chain  -  Render given protein chain

***********************************************************************/

void render_protein_chain(FILE *file_ptr,struct chain *chain_ptr,
                          double scale_factor,int colour_as_orig,
                          int transparent_copies,int biomol)
{
  char chain;

  int transparency;
  int split;

  double rgb[3], rgb_hel[3], rgb_str[3];

  /* Initialise variables */
  split = FALSE;
  transparency = NONE;
  if (biomol == TRUE && chain_ptr->symm_chain == TRUE &&
      transparent_copies == TRUE)
    transparency = PARTIAL;

  /* Write comment to Raster3D file */
  fprintf(file_ptr,"# Start of protein chain %c.\n",chain_ptr->chain_id);

  /* Write out the default material properties */
  material_properties(file_ptr,transparency);

  /* Get the chain */
  chain = chain_ptr->chain_id;
  if (biomol == TRUE && colour_as_orig == TRUE)
    chain = chain_ptr->link_chain_ptr->chain_id;

  /* Get this chain's colour */
  get_chain_col(chain,rgb,rgb_hel,rgb_str,split);

  /* Locate and render all the alpha-helices */
  find_sec_struc(file_ptr,chain_ptr,scale_factor,rgb,HELIX);

  /* Repeat for all the beta strands */
  find_sec_struc(file_ptr,chain_ptr,scale_factor,rgb,STRAND);

  /* Render the C-alpha backbone */
  render_backbone(file_ptr,chain_ptr,scale_factor,rgb);

  /* Label the N-terminal end of this protein chain */
  label_nterminus(file_ptr,chain_ptr,scale_factor);

  /* Reset default material */
  fprintf(file_ptr,"9\n");
}
/***********************************************************************

render_ca_only_chain  -  Render given C-alpha-only protein chain

***********************************************************************/

void render_ca_only_chain(FILE *file_ptr,struct chain *chain_ptr,
                          double scale_factor)
{
  int split;

  double rgb[3], rgb_hel[3], rgb_str[3];

  /* Initialise variables */
  split = FALSE;

  /* Write comment to Raster3D file */
  fprintf(file_ptr,"# Start of protein chain %c (C-alpha only)\n",
          chain_ptr->chain_id);

  /* Write out the default material properties */
  material_properties(file_ptr,NONE);

  /* Get this chain's colour */
  get_chain_col(chain_ptr->chain_id,rgb,rgb_hel,rgb_str,split);

  /* Render the C-alpha backbone */
  render_backbone(file_ptr,chain_ptr,scale_factor,rgb);

  /* Label the N-terminal end of this protein chain */
  label_nterminus(file_ptr,chain_ptr,scale_factor);

  /* Reset default material */
  fprintf(file_ptr,"9\n");
}
/***********************************************************************

write_raster3d_allbonds  -  Write out the given DNA chain to the Raster3D
                            output file as CPK bonds

***********************************************************************/

double write_raster3d_allbonds(FILE *file_ptr,struct atom *first_atom_ptr,
                               int natoms,double scale_factor,
                               int no_mchain_bonds,int grey_carbons)
{
  int atom_number, iatom;
  int mainchain1, mainchain2;

  double min_z;
  static double rgb[3] = { 0.0, 0.0, 0.0 };

  struct atom *atom_ptr, *other_atom_ptr;
  struct conect *conect_ptr;

  /* Write comment to Raster3D file */
  fprintf(file_ptr,"#  Bonds\n");

  /* Initialise variables */
  mainchain1 = FALSE;
  mainchain2 = FALSE;
  min_z = 0.0;

  /* Get pointer to the first atom record */
  atom_ptr = first_atom_ptr;
  iatom = 0;

  /* Loop through all the residue's atoms */
  while (atom_ptr != NULL && iatom < natoms)
    {
      /* Determine whether this is a mainchain atom */
      if (!strcmp(atom_ptr->atom_name," CA ") ||
          !strcmp(atom_ptr->atom_name," N  ") ||
          !strcmp(atom_ptr->atom_name," C  ") ||
          !strcmp(atom_ptr->atom_name," O  "))
        mainchain1 = TRUE;
      else
        mainchain1 = FALSE;

      /* Get the current atom number */
      atom_number = atom_ptr->atom_number;

      /* Get the first of this atom's covalent bonds */
      conect_ptr = atom_ptr->first_conect_ptr;

      /* Loop over this atom's covalent bonds and write out */
      while (conect_ptr != NULL)
        {
          /* Get the other atom */
          other_atom_ptr = conect_ptr->to_atom_ptr;

          /* If this connects to a higher-numbered atom, then
             want this bond */
          if (other_atom_ptr->atom_number > atom_number)
            {
              /* Determine whether this is a mainchain atom */
              if (!strcmp(other_atom_ptr->atom_name," CA ") ||
                  !strcmp(other_atom_ptr->atom_name," N  ") ||
                  !strcmp(other_atom_ptr->atom_name," C  ") ||
                  !strcmp(other_atom_ptr->atom_name," O  "))
                mainchain2 = TRUE;
              else
                mainchain2 = FALSE;

              /* Write out this bond */
              if (no_mchain_bonds == FALSE || mainchain1 == FALSE ||
                  mainchain2 == FALSE)
                write_raster3d_bond(file_ptr,atom_ptr,other_atom_ptr,
                                    scale_factor,BOND_RADIUS,TRUE,rgb,
                                    FALSE,grey_carbons);
            }

          /* Get the next conect */
          conect_ptr = conect_ptr->next_conect_ptr;
        }
  
      /* Get pointer to next atom in linked-list */
      atom_ptr = atom_ptr->next_atom_ptr;
      iatom++;
    }

  /* Return minimum z-value for placing of back wall */
  return(min_z);
}
/***********************************************************************

render_dna_chain  -  Render given DNA chain

***********************************************************************/

void render_dna_chain(FILE *file_ptr,struct chain *chain_ptr,
                      double scale_factor)
{
  int natoms, nresid;
  int split;

  double rgb[3], rgb_hel[3], rgb_str[3];

  struct atom *first_atom_ptr;
  struct residue *first_residue_ptr;

  /* Initialise variables */
  split = FALSE;

  /* Write comment to Raster3D file */
  fprintf(file_ptr,"# Start of DNA chain %c.\n",chain_ptr->chain_id);

  /* Write out the default material properties */
  material_properties(file_ptr,NONE);

  /* Get this chain's colour */
  get_chain_col(chain_ptr->chain_id,rgb,rgb_hel,rgb_str,split);

  /* Calculate all the covalent bonds in the DNA chain */
  first_residue_ptr = chain_ptr->first_residue_ptr;
  nresid = chain_ptr->nresid;
  if (chain_ptr->got_bonds == FALSE)
    {
      calc_covalent_bonds(first_residue_ptr,nresid);
      chain_ptr->got_bonds = TRUE;
    }

  /* Render as all bonds between atoms */
  first_atom_ptr = first_residue_ptr->first_atom_ptr;
  natoms = chain_ptr->natoms;
  write_raster3d_allbonds(file_ptr,first_atom_ptr,natoms,scale_factor,
                          FALSE,FALSE);

  /* Render the phosphate backbone */
  render_backbone(file_ptr,chain_ptr,scale_factor,rgb);

  /* Reset default material */
  fprintf(file_ptr,"9\n");
}
/***********************************************************************

render_ligand_metal  -  Render given ligand or metal

***********************************************************************/

void render_ligand_metal(FILE *file_ptr,struct chain *chain_ptr,
                         double scale_factor)
{
  int iatom, natoms, type;

  static int last_type = -1;
  static int last_metal = -1;

  double min_z;

  struct atom *atom_ptr;
  struct residue *residue_ptr;

  /* Initialise variables */
  min_z = 0.0;

  /* Get the ligand or metal type */
  type = 0;
  residue_ptr = chain_ptr->first_residue_ptr;
  if (residue_ptr->ligand_ptr != NULL)
    type = residue_ptr->ligand_ptr->type;

  /* Write comment to Raster3D file */
  if (chain_ptr->type == LIGAND)
    fprintf(file_ptr,"# Ligand ");
  else
    fprintf(file_ptr,"# Metal ");
  fprintf(file_ptr,"type %3d. %s %s %c.\n",type,residue_ptr->res_name,
          residue_ptr->res_num,residue_ptr->chain);

  /* If new ligand or metal type, write out the default material
     properties */
  if (type != last_type || chain_ptr->type != last_metal)
    {
      /* If previous chain was ligand or metal close off its material
         properties */
      if (last_type != -1)
        fprintf(file_ptr,"9\n");

      /* Set new properties */
      material_properties(file_ptr,NONE);

      /* Save ligand or metal type */
      last_metal = chain_ptr->type;
      last_type = type;
    }

  /* Get the first atom */
  atom_ptr = chain_ptr->first_residue_ptr->first_atom_ptr;
  natoms = chain_ptr->natoms;
  iatom = 0;

  /* Loop through all the residue's atoms */
  while (atom_ptr != NULL && iatom < natoms)
    {
      /* Render as all atoms, shown in spacefill mode */
      min_z = write_raster3d_atom(file_ptr,atom_ptr,scale_factor,
                                  SPACEFILL_PERCEN,min_z);
  
      /* Get pointer to next atom in linked-list */
      atom_ptr = atom_ptr->next_atom_ptr;
      iatom++;
    }
}
/***********************************************************************

render_ssbonds  -  Render any disulphide bonds

***********************************************************************/

void render_ssbonds(FILE *file_ptr,struct ssbond *first_ssbond_ptr,
                    char chain,double scale_factor)
{
  int iatom, loop, natoms;
  int found;

  static double rgb[3] = { 1.0, 1.0, 0.0 }; 

  struct atom *atom_ptr, *ca_atom_ptr, *ca_atom1_ptr, *ca_atom2_ptr;
  struct residue *residue_ptr, *residue1_ptr, *residue2_ptr;
  struct ssbond *ssbond_ptr;

  /* Get the first disulphide bond */
  ssbond_ptr = first_ssbond_ptr;

  /* Loop through the bonds, writing each out to the Raster3D file */
  while (ssbond_ptr != NULL)
    {
      /* Initialise pointers to their two C-alpha atoms */
      ca_atom1_ptr = NULL;
      ca_atom2_ptr = NULL;

      /* Loop over the two residues to locate their C-alpha atoms */
      for (loop = 0; loop < 2; loop++)
        {
          /* Get appropriate residue */
          if (loop == 0)
            residue_ptr = ssbond_ptr->residue1_ptr;
          else
            residue_ptr = ssbond_ptr->residue2_ptr;

          /* Get the first of this residue's atoms */
          atom_ptr = residue_ptr->first_atom_ptr;
          iatom = 0;
          natoms = residue_ptr->natoms;
  
          /* Save pointer to the first atom, just in case we can't find
             the c-alpha itself */
          ca_atom_ptr = atom_ptr;
          found = FALSE;

          /* Loop through all the residue's atoms to find the C-alpha
             or phosphate atom */
          while (atom_ptr != NULL && iatom < natoms && found == FALSE)
            {
              /* If this is C-alpha atom, store its pointer */
              if (!strcmp(atom_ptr->atom_name," CA "))
                {
                  /* Store pointer and set flag */
                  ca_atom_ptr = atom_ptr;
                  found = TRUE;
                }

              /* Get pointer to next atom in linked-list */
              atom_ptr = atom_ptr->next_atom_ptr;
              iatom++;
            }

          /* If found, store in appropriate variable */
          if (ca_atom_ptr != NULL && loop == 0)
            ca_atom1_ptr = ca_atom_ptr;
          else
            ca_atom2_ptr = ca_atom_ptr;
        }

      /* If have both atoms, draw a line between them */
      if (ca_atom1_ptr != NULL && ca_atom2_ptr != NULL)
        {
          /* Get the two residue */
          residue1_ptr = ca_atom1_ptr->residue_ptr;
          residue2_ptr = ca_atom2_ptr->residue_ptr;

          /* If plotting bonds from a single chain only, make sure
             this disulphide bond belongs to that chain */
          if (chain == '\0' || (chain == residue1_ptr->chain &&
                                chain == residue2_ptr->chain))
            {
              /* Write comment to Raster3D file */
              fprintf(file_ptr,
                      "#  Disulphide bond: %s %s %c to %s %s %c.\n",
                      residue1_ptr->res_name,residue1_ptr->res_num,
                      residue1_ptr->chain,residue2_ptr->res_name,
                      residue2_ptr->res_num,residue2_ptr->chain);

              /* Write out the default material properties */
              material_properties(file_ptr,NONE);

              /* Write out the bond */
              write_raster3d_bond(file_ptr,ca_atom1_ptr,ca_atom2_ptr,
                                  scale_factor,SSBOND_RADIUS,FALSE,
                                  rgb,FALSE,FALSE);

              /* Reset default material */
              fprintf(file_ptr,"9\n");
            }
        }

      /* Get the next ssbond */
      ssbond_ptr = ssbond_ptr->next_ssbond_ptr;
    }
}
/***********************************************************************

render_sites  -  Render active site residues as all-atoms

***********************************************************************/

void render_sites(FILE *file_ptr,struct chain *first_chain_ptr,
                  char chain,double scale_factor)
{
  int iresid, nresid;

  struct chain *chain_ptr;
  struct residue *residue_ptr;

  /* Get pointer to the first chain */
  chain_ptr = first_chain_ptr;

  /* Loop over all chains */
  while (chain_ptr != NULL)
    {
      /* Get this chain's first residue */
      residue_ptr = chain_ptr->first_residue_ptr;
      nresid = chain_ptr->nresid;
      iresid = 0;

      /* Loop through this chain's residues to locate any active site
         residues */
      while (residue_ptr != NULL && iresid < nresid)
        {
          /* If this is an active site residue, then process */
          if (residue_ptr->site[0] != '\0' &&
              (chain == '\0' || chain == residue_ptr->chain))
            {
              /* Calculate all the covalent bonds in this residue */
              if (residue_ptr->got_bonds == FALSE)
                {
                  calc_covalent_bonds(residue_ptr,1);
                  residue_ptr->got_bonds = TRUE;
                }

              /* Write comment to Raster3D file */
              fprintf(file_ptr,"#  Active site: Chain %c Residue %s %s.\n",
                      chain_ptr->chain_id,residue_ptr->res_name,
                      residue_ptr->res_num);

              /* Write out the default material properties */
              material_properties(file_ptr,NONE);

              /* Render as all bonds between atoms */
              write_raster3d_allbonds(file_ptr,
                                      residue_ptr->first_atom_ptr,
                                      residue_ptr->natoms,
                                      scale_factor,TRUE,TRUE);

              /* Reset default material */
              fprintf(file_ptr,"9\n");
            }

          /* Get the next residue */
          residue_ptr = residue_ptr->next_residue_ptr;
          iresid++;
        }

      /* Get the next chain */
      chain_ptr = chain_ptr->next_chain_ptr;
    }
}
/***********************************************************************

pymol_render  -  Write ray-trace command to PyMol script

***********************************************************************/

void pymol_render(FILE *pymol_script_ptr,char *png_name,int tiny)
{
  /* Write comment */
  fprintf(pymol_script_ptr,"# Render %s\n",png_name);

  /* Write out ray-trace command */
  if (tiny == TRUE)
    {
      fprintf(pymol_script_ptr,"ray 120, 120\n");
      fprintf(pymol_script_ptr,"viewport 120, 120\n");
    }
  else
    {
      fprintf(pymol_script_ptr,"ray 240, 240\n");
      fprintf(pymol_script_ptr,"viewport 240, 240\n");
    }

  /* Write out command to generate png file */
  fprintf(pymol_script_ptr,"png %s\n",png_name);
}
/***********************************************************************

create_pymol_record  -  Create and initialise a new pymol record

***********************************************************************/

struct pymol *create_pymol_record(struct pymol **lst_pymol_ptr,
                                  char *gif_name,char *png_name,int crop)
{
  struct pymol *pymol_ptr, *last_pymol_ptr;

  /* Initialise pymol pointers */
  pymol_ptr = NULL;
  last_pymol_ptr = *lst_pymol_ptr;

  /* Allocate memory for structure to hold pymol info */
  pymol_ptr = (struct pymol *) malloc(sizeof(struct pymol));
  if (pymol_ptr == NULL)
    {
      printf("*** Can't allocate memory for struct pymol\n");
      exit (1);
    }

  /* Store pointer to previous record */
  pymol_ptr->next_pymol_ptr = last_pymol_ptr;
  last_pymol_ptr = pymol_ptr;

  /* Store the pymol details */
  store_string(&(pymol_ptr->gif_name),gif_name);
  store_string(&(pymol_ptr->png_name),png_name);
  pymol_ptr->crop = crop;

  /* Return current pointer */
  *lst_pymol_ptr = last_pymol_ptr;

  /* Return new pymol pointer */
  return(pymol_ptr);
}
/***********************************************************************

render_structure  -  Render whole structure in one or three orthogonal
                     views

***********************************************************************/

int render_structure(struct chain *first_chain_ptr,
                     struct atom *first_atom_ptr,int natoms,int nviews,
                     struct ssbond *first_ssbond_ptr,int nssbonds,
                     int nsite_res,char *render_exe,char *convert_exe,
                     char *ppmquant_exe,char *ppmtogif_exe,int mcsg,
                     double *grid_sep,int nligtypes,int nmetypes,
                     FILE *pymol_script_ptr,struct pymol **lst_pymol_ptr,
                     int colour_as_orig,int transparent_copies,
                     int biomol,struct parameters params)
{
  char file_r3d[FILENAME_LEN + 1], gif_name[FILENAME_LEN + 1];
  char image[20], png_name[FILENAME_LEN + 1];

  int iview, ntype[2];
  int crop, first_dna, last_type, tiny, transparent;

  double coord_min[3], coord_max[3], mass_x, mass_y, mass_z;
  double grid_separation, scale_factor, volume;

  static double matrix[2][3][3]
    = {  0,  0,  1,
         0,  1,  0,
         -1,  0,  0,
         1,  0,  0,
         0,  0, -1,
         0,  1,  0
    };

  static double rmatrix[2][3][3]
    = {  0,  0, -1,
         0,  1,  0,
         1,  0,  0,
         1,  0,  0,
         0,  0,  1,
         0, -1,  0
    };

  struct chain *chain_ptr;
  struct pymol *last_pymol_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  crop = FALSE;
  grid_separation = 0.8;

  /* Initialise pointers */
  last_pymol_ptr = *lst_pymol_ptr;

  /* Name of temporary Raster3D input file */
  strcpy(file_r3d,"s.r3d");

  /* Loop over the required number of views */
  for (iview = 0; iview < nviews; iview++)
    {
      /* Initialise count of object types (chains or ligand/metals) */
      ntype[CHAIN] = 0;
      ntype[LIGMET] = nligtypes + nmetypes;
      first_dna = TRUE;
      last_type = -1;

      /* Form name of the GIF file to be generated */
      strcpy(gif_name,"traces");
      if (iview == 1)
        strcat(gif_name,"_r");
      else if (iview == 2)
        strcat(gif_name,"_b");
      if (biomol == TRUE)
        strcpy(gif_name,"biomol");
      strcpy(png_name,gif_name);
      strcat(png_name,".png");
      if (IMAGE_FORMAT == GIF)
        strcat(gif_name,".gif");
      else
        strcat(gif_name,".jpg");

      /* Form name of this view */
      sprintf(image,"View %d",iview+1);

      /* Get centre of mass of all the atoms */
      centre_of_mass(first_atom_ptr,natoms,&mass_x,&mass_y,&mass_z);

      /* Place the centre of mass at the origin */
      adjust_to_origin(first_atom_ptr,natoms,mass_x,mass_y,mass_z,
                       TRUE);

      /* If generating images using PyMol, write appropriate turn
         command */
      if (pymol_script_ptr != NULL)
        {
          /* Turn to get right view */
          if (iview == 1)
            fprintf(pymol_script_ptr,"turn y, -90\n");

          /* Otherwise, reverse previous turn before going for bottom
             view */
          else if (iview == 2)
            {
              /* Apply the two turns */
              fprintf(pymol_script_ptr,"turn y, 90\n");
              fprintf(pymol_script_ptr,"turn x, -90\n");
            }
        }

      /* For old method, open Raster3D file */
      else
        {
          /* Open the .r3d output file */
          file_ptr = open_output_file(file_r3d,TRUE);

          /* Write out the Raster3D header records */
          write_raster3d_headers(file_ptr);

          /* For 2nd or 3rd view, apply appropriate orthogonal rotation
             to all the atoms */
          if (iview > 0)
            {
              /* Apply appropriate rotation matrix */
              rotate_structure(first_atom_ptr,natoms,matrix[iview - 1]);
            }
        }

      /* Get maximum and minimum atom coords */
      get_max_min(first_atom_ptr,natoms,coord_min,coord_max);

      /* Calculate Reaster3D scaling factor */
      scale_factor = calc_r3d_scale(coord_min,coord_max);

      /* On first loop, calculate the most appropriate grid-spacing
         for the SURFNET run on the current structure */
      if (iview == 0)
        {
          volume = (coord_max[0] - coord_min[0])
            * (coord_max[1] - coord_min[1])
            * (coord_max[2] - coord_min[2]);
          grid_separation = pow((volume / 175000),(1.0/3.0));
          if (grid_separation > 1.0)
            grid_separation = pow(grid_separation,1.35);
          if (grid_separation < 0.6)
            grid_separation = 0.6;
        }

      /* If writing to Raster3D file, write drawing commands */
      if (pymol_script_ptr == NULL)
        {
          /* Get pointer to the first chain record */
          chain_ptr = first_chain_ptr;

          /* Loop over the chains */
          while (chain_ptr != NULL)
            {
              /* If have reached end of ligands/metals, close off
                 material properties */
              if (chain_ptr->type != LIGAND &&
                  chain_ptr->type != METAL &&
                  (last_type == LIGAND || last_type == METAL))
                fprintf(file_ptr,"9\n");

              /* If this is a protein chain, show secondary structure
                 schematically */
              if (chain_ptr->type == PROTEIN)
                {
                  /* Render this protein chain */
                  render_protein_chain(file_ptr,chain_ptr,scale_factor,
                                       colour_as_orig,transparent_copies,
                                       biomol);
                  ntype[CHAIN]++;
                }

              /* If this is a C-alpha only protein chain, show C-alpha
                 backbone */
              else if (chain_ptr->type == CALPHA_ONLY)
                {
                  /* Render this protein chain as a C-alpha backbone
                     trace */
                  render_ca_only_chain(file_ptr,chain_ptr,scale_factor);
                  ntype[CHAIN]++;
                }

              /* If this is a DNA/RNA chain, draw as all atoms, plus
                 backbone */
              else if (chain_ptr->type == DNA)
                {
                  /* Render this DNA chain as all atoms plus backbone
                     trace */
                  render_dna_chain(file_ptr,chain_ptr,scale_factor);
                  if (first_dna == TRUE)
                    ntype[CHAIN]++;
                  first_dna = FALSE;
                }

              /* If this is a ligand or metal, show as all atoms */
              else if (chain_ptr->type == LIGAND ||
                       chain_ptr->type == METAL)
                {
                  /* Render this ligand or metal as spacefill */
                  render_ligand_metal(file_ptr,chain_ptr,scale_factor);
                }

              /* Save the chain type */
              last_type = chain_ptr->type;

              /* Get the next chain */
              chain_ptr = chain_ptr->next_chain_ptr;
            }

          /* If the last type was metal or ligand, close off material
             properties */
          if (last_type == LIGAND || last_type == METAL)
            fprintf(file_ptr,"9\n");

          /* If structure contains any disulphide bonds, show those */
          if (nssbonds > 0)
            render_ssbonds(file_ptr,first_ssbond_ptr,'\0',scale_factor);

          /* If structure contains any active site residues, show those */
          if (nsite_res > 0)
            render_sites(file_ptr,first_chain_ptr,'\0',scale_factor);

          /* Close the .r3d output file */
          fclose(file_ptr);
        }

      /* Otherwise, loop through the chains to count up numbers */
      else
        {
          /* Get pointer to the first chain record */
          chain_ptr = first_chain_ptr;

          /* Loop over the chains */
          while (chain_ptr != NULL)
            {
              if (chain_ptr->type == PROTEIN ||
                  chain_ptr->type == CALPHA_ONLY ||
                  chain_ptr->type == DNA)
                ntype[CHAIN]++;

              /* Get the next chain */
              chain_ptr = chain_ptr->next_chain_ptr;
            }
        }

      /* Set flag */
      tiny = FALSE;

      /* If generating PyMol image, add render commands to pymol.pml */
      if (pymol_script_ptr != NULL)
        {
          /* Add render command */
          pymol_render(pymol_script_ptr,png_name,tiny);

          /* Save image names */
          create_pymol_record(&last_pymol_ptr,gif_name,png_name,crop);
        }

      /* Otherwise, use old Raster3D method to render structure */
      else if (params.run_raster3d == TRUE)
        run_raster3d(file_r3d,gif_name,render_exe,convert_exe,
                     ppmquant_exe,ppmtogif_exe,image,crop,tiny);

      /* If this is an MCSG structure, generate a tiny version of
         the main view image */
      if (params.run_raster3d == TRUE && mcsg == TRUE && iview == 0)
        {
          /* Form name for this image file */
          if (IMAGE_FORMAT == GIF)
            strcpy(gif_name,"mcsg_trace.gif");
          else
            strcpy(gif_name,"mcsg_trace.jpg");
          strcpy(png_name,"mcsg_trace.png");

          /* Set flag */
          tiny = TRUE;

          /* If generating PyMol image, add render commands to pymol.pml */
          if (pymol_script_ptr != NULL)
            {
              /* Add render command */
              pymol_render(pymol_script_ptr,png_name,tiny);

              /* Save image names */
              create_pymol_record(&last_pymol_ptr,gif_name,png_name,crop);
            }

          /* Otherwise, use old Raster3D method to render structure */
          else
            run_raster3d(file_r3d,gif_name,render_exe,convert_exe,
                         ppmquant_exe,ppmtogif_exe,image,crop,tiny);
        }

      /* If have rotated the structure, rotate it back */
      if (pymol_script_ptr == NULL && iview > 0)
        rotate_structure(first_atom_ptr,natoms,rmatrix[iview - 1]);

      /* Rename the temporary Raster3D input file */
      strcpy(file_r3d,"l.r3d");
      sprintf(file_r3d,"s%d.r3d",iview);
    }

  /* Based on the number of chains and the number of ligands/metals
     determine whether need to generate transparent images */
  transparent = FALSE;
  if (ntype[CHAIN] > 1 || (ntype[CHAIN] > 0 && ntype[LIGMET] > 0))
    transparent = TRUE;

  /* If generating PyMol image, reset to original orientation */
  if (pymol_script_ptr != NULL)
    fprintf(pymol_script_ptr,"reset\n");
  else
    transparent = FALSE;

  /* Return the calculated grid separation for use by SURFNET */
  *grid_sep = grid_separation;

  /* Return the last PyMol image record */
  *lst_pymol_ptr = last_pymol_ptr;

  /* Return whether transparent plots to be generated */
  return(transparent);
}
/***********************************************************************

select_pymol  -  Select the objects to be hidden and made visible, and
                 add transparency where required

***********************************************************************/

void select_pymol(struct pobject *first_pobject_ptr,int type,
                  char *chains,char *copy_chain,int ncopies,
                  int ligmet_type,FILE *pymol_script_ptr,
                  struct pymol **lst_pymol_ptr,char *gif_name,
                  int transparent)
{
  char display_type[TOKEN_LEN + 1], inchain, object[TOKEN_LEN + 1];
  char png_name[FILENAME_LEN + 1];

  int ichain, icopy, len, object_type, transparency;
  int crop, first, found, show, tiny;

  struct pobject *pobject_ptr;
  struct pymol *last_pymol_ptr;

  /* Initialise variables */
  crop = FALSE;
  first = TRUE;
  tiny = FALSE;

  /* Initialise pointers */
  last_pymol_ptr = *lst_pymol_ptr;

  /* Form name of the .png file for this image */
  strcpy(png_name,gif_name);
  len = strlen(png_name);
  png_name[len - 3] = '\0';
  strcat(png_name,"png");

  /* Initialise picture by hiding everything */
  fprintf(pymol_script_ptr,"hide all\n");

  /* Get pointer to the first pobject record */
  pobject_ptr = first_pobject_ptr;
  
  /* Loop over the pobjects */
  while (pobject_ptr != NULL)
    {
      /* Get this object's chain-id and type */
      inchain = pobject_ptr->chain_id;
      object_type = pobject_ptr->type;

      /* Initialise flags */
      show = TRUE;
      transparency = NONE;

      /* If showing protein or protein-protein interface, determine
         whether this object is to be shown */
      if (type == PROTEIN || type == PP_INTERFACE)
        {
          /* If object is not a protein, or part of a protein
             hide */
          if (object_type != PROTEIN && object_type != CALPHA_ONLY &&
              object_type != DISULPH && object_type != ACTSITE)
            show = FALSE;

          /* If this is the correct protein chain then show */
          else if (inchain == chains[0] || inchain == chains[1])
            show = TRUE;

          /* Otherwise, check whether this is a copy chain */
          else if (transparent == TRUE)
            {
              /* Initialise found flag */
              found = FALSE;

              /* Loop over the copy chains */
              for (icopy = 0; icopy < ncopies; icopy++)
                {
                  /* If this is the right chain, mark as found */
                  if (inchain == copy_chain[icopy])
                    found = TRUE;
                }

              /* If chain identified as a copy chain, show partially
                 transparent */
              if (found == TRUE)
                transparency = PARTIAL;

              /* Otherwise, make fully transparent */
              else
                transparency = FULL;
            }

          /* Otherwise, don't show this object */
          else
            show = FALSE;
        }

      /* If showing DNA, determine whether to show object or not */
      else if (type == DNA)
        {
          /* If object is DNA, then show */
          if (object_type == DNA)
            show = TRUE;

          /* Otherwise, don't show */
          else
            show = FALSE;
        }


      /* If showing ligand or metal, determine whether this is the
         right one */
      else if (type == LIGAND || type == METAL)
        {
          /* If object is ligand, check ligand number */
          if (object_type == type &&
              pobject_ptr->ligmet_type == ligmet_type)
            show = TRUE;

          /* Otherwise, don't show */
          else
            show = FALSE;
        }

      /* If showing transparent image, define level of transparency */
      if (transparent == TRUE)
        {
          /* If transparency hasn't yet been set determine what it
             should be */
          if (transparency == NONE)
            {
              /* If object not to be shown set transparency as full */
              if (show == FALSE)
                {
                  transparency = FULL;
                  show = TRUE;
                }
            }
        }

      /* Get the object's name and display type */
      strcpy(object,pobject_ptr->name);
      strcpy(display_type,pobject_ptr->display_type);

      /* If object to be shown, display it */
      if (show == TRUE)
        {
          /* Display object */
          fprintf(pymol_script_ptr,"show %s, %s\n",display_type,object);

          /* If display mode ends in an 's', remove it for the transparency
             command */
          len = strlen(display_type);
          if (display_type[len - 1] == 's')
            display_type[len - 1] = '\0';

          /* Define object's transparency */
          fprintf(pymol_script_ptr,"set %s_transparency, %3.1f, %s\n",
                  display_type,PYMOL_TRANSPARENCY[transparency],object);
          
          /* Add to selection for zoom */
          if (first == TRUE)
            {
              if (transparent == FALSE)
                fprintf(pymol_script_ptr,"select vismol, %s\n",object);
              first = FALSE;
            }
          else if (transparent == FALSE)
            fprintf(pymol_script_ptr,"select vismol, (vismol or %s)\n",
                    object);
        }

      /* Get the next pobject */
      pobject_ptr = pobject_ptr->next_pobject_ptr;
    }

  /* If have a view object, zoom in on it */
  if (first == FALSE)
    {
      /* Zoom on appropriate object */
      if (transparent == TRUE)
        fprintf(pymol_script_ptr,"zoom all, ");
      else
        fprintf(pymol_script_ptr,"zoom vismol, ");

      /* Add additional zoom parameters to ensure whole object is in
         field of view */
      fprintf(pymol_script_ptr,"complete=1, buffer=2\n");
      
    }

  /* Add render command */
  pymol_render(pymol_script_ptr,png_name,tiny);

  /* Save image names */
  create_pymol_record(&last_pymol_ptr,gif_name,png_name,crop);

  /* Return the last PyMol image record */
  *lst_pymol_ptr = last_pymol_ptr;
}
/***********************************************************************

transparent_image  -  Edit the s.r3d Raster3D input file containing the
                      entire structure to make all objects except the
                      ones of interest transparent

***********************************************************************/

int transparent_image(struct chain *first_chain_ptr,
                      struct ligand *first_ligand_ptr,int type,
                      char *chains,char *copy_chain,int ncopies,
                      int ligmet_type,char *gif,char *render_exe,
                      char *convert_exe,char *ppmquant_exe,
                      char *ppmtogif_exe,FILE *pymol_script_ptr,
                      struct pobject *first_pobject_ptr,
                      struct pymol **lst_pymol_ptr,
                      struct parameters params)
{
  char s_r3d[FILENAME_LEN + 1], l_r3d[FILENAME_LEN + 1];
  char input_line[LINELEN], gif_name[FILENAME_LEN + 1];
  char image[60], inchain, number_string[20], tmp_string[60];

  char *token[MAX_TOKEN];

  int icopy, itoken, len, ltype, ntokens, object_type;
  int transparency, transparent;
  int found, have_file, have_sphere, skip_next, translucent, wanted;

  double radius;

  struct pymol *last_pymol_ptr;

  FILE *infile_ptr, *outfile_ptr;

  /* Initialise variables */
  have_file = TRUE;
  have_sphere = FALSE;
  skip_next = FALSE;
  object_type = -1;
  transparency = NONE;
  transparent = TRUE;
  wanted = TRUE;

  /* Initialise pointers */
  last_pymol_ptr = *lst_pymol_ptr;

  /* Form comment describing image type */
  strcpy(image,"Transparent plot for ");
  if (type == PROTEIN)
    strcat(image,"protein");
  else if (type == CALPHA_ONLY)
    strcat(image,"C-alpha only");
  else if (type == DNA)
    strcat(image,"DNA");
  else if (type == LIGAND)
    {
      strcat(image,"LIGAND type ");
      sprintf(number_string,"%d",ligmet_type);
      strcat(image,number_string);
    }
  else if (type == METAL)
    {
      strcat(image,"METAL type ");
      sprintf(number_string,"%d",ligmet_type);
      strcat(image,number_string);
    }
  else if (type == PP_INTERFACE)
    strcat(image,"Protein-protein interface ");
  else
    strcat(image,"UNKNOWN");
  if (type == PP_INTERFACE)
    {
      sprintf(tmp_string,"%c-%c",chains[0],chains[1]);
      strcat(image,tmp_string);
    }
  else if (type == PROTEIN || type == CALPHA_ONLY)
    {
      strcat(image," chain [");
      strcat(image,chains);
      strcat(image,"]");
    }

  /* If generating PyMol image, write comment and render transparent 
     image */
  if (pymol_script_ptr != NULL)
    {
/* debug */
  printf("  ... Writing PyMOL script\n");
  fflush(stdout);
/* debug */
      /* If not a ligand or metal, generate solid version of this
         molecule */
      if (type != LIGAND && type != METAL)
        {
          /* Form name of output image file */
          strcpy(gif_name,gif);

          /* For protein-protein interface, name is different */
          if (type == PP_INTERFACE)
            {
              strcpy(gif_name,"chains");
              strcat(gif_name,gif+7);
            }

          /* Write comment for solid plot */
          fprintf(pymol_script_ptr,"# Solid %s\n",image+12);

          /* Set transparency flag */
          translucent = FALSE;

          /* Generate PyMol image for select part of structure only */
          select_pymol(first_pobject_ptr,type,chains,copy_chain,ncopies,
                       ligmet_type,pymol_script_ptr,&last_pymol_ptr,
                       gif_name,translucent);
        }

      /* Form name of output image file for transparent plot */
      strcpy(gif_name,"t");
      strcat(gif_name,gif);

      /* Write comment for transparent plot */
      fprintf(pymol_script_ptr,"# %s\n",image);

      /* Set transparency flag */
      translucent = TRUE;

      /* Generate PyMol image with corresponding parts transparent */
      select_pymol(first_pobject_ptr,type,chains,copy_chain,ncopies,
                   ligmet_type,pymol_script_ptr,&last_pymol_ptr,
                   gif_name,translucent);
    }

  /* Otherwise, use old system with Raster3D */
  else
    {
      /* Form name of output GIF (or JPEG) file */
      strcpy(gif_name,"t");
      strcat(gif_name,gif);

      /* Names of input and output Raster3D files */
      strcpy(s_r3d,"s.r3d");
      strcpy(l_r3d,"l.r3d");

      /* Open input s.r3d file */
      if ((infile_ptr = fopen(s_r3d,"r")) !=  NULL)
        {
          /* Open the l.r3d output file */
          outfile_ptr = open_output_file(l_r3d,TRUE);

          /* Loop while reading in records from the input file */
          while (fgets(input_line,LINELEN,infile_ptr)!= NULL)
            {
              /* Truncate the string to strip off the newline */
              len = string_truncate(input_line,LINELEN);

              /* Assume current line is wanted */
              wanted = TRUE;

              /* If skipping current line set flags */
              if (skip_next == TRUE)
                {
                  wanted = FALSE;
                  skip_next = FALSE;
                }

              /* Check for comment line */
              if (input_line[0] == '#')
                {
                  /* Check whether this comment defines the start of a new
                     protein chain, active site residue or disulphide
                     bond */
                  if (!strncmp(input_line,"# Start of protein chain",24) ||
                      !strncmp(input_line,"#  Active site: Chain",21) ||
                      !strncmp(input_line,"#  Disulphide bond:",19))
                    {
                      /* Get chain identifier */
                      if (!strncmp(input_line,
                                   "# Start of protein chain",24))
                        inchain = input_line[25];
                      else if (!strncmp(input_line,
                                        "#  Active site: Chain",21))
                        inchain = input_line[22];
                      else if (input_line[30] == input_line[45])
                        inchain = input_line[30];
                      else
                        inchain = '\0';

                      /* Set object type */
                      object_type = PROTEIN;

                      /* If highlighting something other than a protein
                         chain, make this chain fully transparent */
                      if (type != PROTEIN && type != PP_INTERFACE)
                        transparency = FULL;

                      /* If this is the protein chain to be highlighted,
                         make non-transparent */
                      else if (inchain == chains[0] || inchain == chains[1])
                        transparency = NONE;

                      /* Otherwise, check whether this is a copy chain
                         to be made partially transparent */
                      else
                        {
                          /* Initialise found flag */
                          found = FALSE;

                          /* Loop over the copy chains */
                          for (icopy = 0; icopy < ncopies; icopy++)
                            {
                              /* If this is the right chain, mark as
                                 found */
                              if (inchain == copy_chain[icopy])
                                found = TRUE;
                            }

                          /* If chain identified as a copy chain, show
                             partially transparent */
                          if (found == TRUE)
                            transparency = PARTIAL;

                          /* Otherwise, make fully transparent */
                          else
                            transparency = FULL;
                        }
                    }

                  /* Start of DNA chain */
                  else if (!strncmp(input_line,"# Start of DNA chain",20))
                    {
                      /* Set object type */
                      object_type = DNA;

                      /* If highlighting DNA chains, make non-transparent */
                      if (type == DNA)
                        transparency = NONE;

                      /* Otherwise, make fully transparent */
                      else
                        transparency = FULL;
                    }

                  /* Start of ligand */
                  else if (!strncmp(input_line,"# Ligand",8))
                    {
                      /* Set object type */
                      object_type = LIGAND;

                      /* If anything other than a ligand, make transparent */
                      if (type != LIGAND)
                        transparency = FULL;

                      /* Otherwise, get the ligand type and check if it is
                         the right one */
                      else
                        {
                          /* Get the ligand type */
                          strncpy(number_string,input_line+14,3);
                          number_string[3] = '\0';
                          ltype = atoi(number_string);

                          /* If correct ligand type, make visible */
                          if (ltype == ligmet_type)
                            transparency = NONE;

                          /* Otherwise, make fully transparent */
                          else
                            transparency = FULL;
                        }
                    }

                  /* Start of metal */
                  else if (!strncmp(input_line,"# Metal",7))
                    {
                      /* Set object type */
                      object_type = METAL;

                      /* If anything other than a metal, make transparent */
                      if (type != METAL)
                        transparency = FULL;

                      /* Otherwise, get the metal type and check if it is
                         the right one */
                      else
                        {
                          /* Get the metal type */
                          strncpy(number_string,input_line+13,3);
                          number_string[3] = '\0';
                          ltype = atoi(number_string);

                          /* If correct metal type, make visible */
                          if (ltype == ligmet_type)
                            transparency = NONE;

                          /* Otherwise, make fully transparent */
                          else
                            transparency = FULL;
                        }
                    }

                  /* Check for the start of the material properties
                     definition */
                  else if (!strcmp(input_line,"#8"))
                    {
                      /* Write out the material properties for the current
                         object */
                      material_properties(outfile_ptr,transparency);

                      /* Current line and next line not to be written
                         out to file */
                      wanted = FALSE;
                      skip_next = TRUE;
                    }
                }

              /* If have a sphere for a ligand or metal atom, need
                 to check it's not below the minimum visible size */
              if (have_sphere == TRUE &&
                  (object_type == METAL || object_type == LIGAND))
                {
                  /* Read in all the tokens of this line */
                  ntokens = get_tokens(input_line,token,MAX_TOKEN,' ');

                  /* If have exactly 7 tokens, then can check atom
                     radius */
                  if (ntokens == 7)
                    {
                      /* Get the atom radius (in Raster3D coords) */
                      radius = atof(token[3]);

                      /* If less than the minimum visible, replace by
                         minimum value */
                      if (radius < MIN_ATOM_SIZE)
                        {
                          sprintf(token[3],"%6.3f",MIN_ATOM_SIZE);

                          /* Recreate the input line */
                          input_line[0] = '\0';
                          for (itoken = 0; itoken < ntokens; itoken++)
                            {
                              strcat(input_line,token[itoken]);
                              strcat(input_line," ");
                            }
                        }
                    }

                  /* Free up memory taken up by the tokens */
                  free_tokens(token,ntokens);
                }

              /* Check for start of sphere definition */
              if (len == 1 && input_line[0] == '2')
                have_sphere = TRUE;
              else
                have_sphere = FALSE;

              /* If wanted, write line out to output file */
              if (wanted == TRUE)
                fprintf(outfile_ptr,"%s\n",input_line);
            }

          /* Close both files */
          fclose(infile_ptr);
          fclose(outfile_ptr);
        }

      /* If file has gone missing, show error */
      else
        {
          printf("* Warning. s.r3d file has gone missing!! Can't "
                 "generate transparent images\n");
          transparent = FALSE;
          have_file = FALSE;
        }

      /* Run Raster3D to render the image */
      if (have_file == TRUE && params.run_raster3d == TRUE)
        run_raster3d(l_r3d,gif_name,render_exe,convert_exe,ppmquant_exe,
                     ppmtogif_exe,image,FALSE,FALSE);
    }

  /* Return the last PyMol image record */
  *lst_pymol_ptr = last_pymol_ptr;

  /* Return whether to continue with transparency */
  return(transparent);
}
/***********************************************************************

render_all_dna  -  Render whole structure in one or three orthogonal
                     views

***********************************************************************/

int render_all_dna(struct chain *first_chain_ptr,
                   struct ligand *first_ligand_ptr,
                   struct atom *first_atom_ptr,int natoms,
                   char *render_exe,char *convert_exe,
                   char *ppmquant_exe,char *ppmtogif_exe,
                   FILE *pymol_script_ptr,
                   struct pobject *first_pobject_ptr,
                   struct pymol **lst_pymol_ptr,int transparent,
                   struct parameters params)
{
  char file_r3d[FILENAME_LEN + 1], gif_name[FILENAME_LEN + 1];
  char chains[3], image[20];

  int i, ndna_atoms;
  int first;

  double cmin[3], cmax[3], cofg[3], coord_min[3], coord_max[3];
  double mass_x, mass_y, mass_z;
  double scale_factor;

  struct chain *chain_ptr;
  struct pymol *last_pymol_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  chains[0] = '\0';
  first = TRUE;
  for (i = 0; i < 3; i++)
    cofg[i] = 0.0;
  ndna_atoms = 0;
  strcpy(image,"All DNA chains");

  /* Initialise pointers */
  last_pymol_ptr = *lst_pymol_ptr;

  /* Name of temporary Raster3D input file */
  strcpy(file_r3d,"l.r3d");

  /* Form name of the GIF file to be generated */
  if (IMAGE_FORMAT == GIF)
    strcpy(gif_name,"dna.gif");
  else
    strcpy(gif_name,"dna.jpg");

  /* If need to generate transparent image, do so now */
  if (transparent == TRUE || pymol_script_ptr != NULL)
    transparent
      = transparent_image(first_chain_ptr,first_ligand_ptr,DNA,chains,
                          "0",0,0,gif_name,render_exe,convert_exe,
                          ppmquant_exe,ppmtogif_exe,pymol_script_ptr,
                          first_pobject_ptr,&last_pymol_ptr,params);

  /* For old method, generate images using Raster3D */
  if (pymol_script_ptr == NULL)
    {
      /* Open the l.r3d output file */
      file_ptr = open_output_file(file_r3d,TRUE);

      /* Write out the Raster3D header records */
      write_raster3d_headers(file_ptr);

      /* Get pointer to the first chain record */
      chain_ptr = first_chain_ptr;

      /* Loop over the chains */
      while (chain_ptr != NULL)
        {
          /* Include only if this is a DNA chain */
          if (chain_ptr->type == DNA)
            {
              /* Get centre of mass of all atoms in this chain */
              centre_of_mass(chain_ptr->first_residue_ptr->first_atom_ptr,
                             chain_ptr->natoms,&mass_x,&mass_y,&mass_z);

              /* Update overall centre of mass position */
              cofg[0] = cofg[0] + chain_ptr->natoms * mass_x;
              cofg[1] = cofg[1] + chain_ptr->natoms * mass_y;
              cofg[2] = cofg[2] + chain_ptr->natoms * mass_z;
              ndna_atoms = ndna_atoms + chain_ptr->natoms;
            }

          /* Get the next chain */
          chain_ptr = chain_ptr->next_chain_ptr;
        }

      /* Calculate overall centre of mass */
      if (ndna_atoms > 0)
        {
          for (i = 0; i < 3; i++)
            cofg[i] = cofg[i] / ndna_atoms;
        }

      /* Adjust whole structure to centre of mass of just the DNA
         chains */
      adjust_to_origin(first_atom_ptr,natoms,cofg[0],cofg[1],cofg[2],
                       TRUE);

      /* Get pointer to the first chain record again */
      chain_ptr = first_chain_ptr;

      /* Loop over the chains to calculate maximum and minimum coords of
         just the DNA chains */
      while (chain_ptr != NULL)
        {
          /* Include only if this is a DNA chain */
          if (chain_ptr->type == DNA)
            {
              /* Get maximum and minimum atom coords */
              get_max_min(chain_ptr->first_residue_ptr->first_atom_ptr,
                          chain_ptr->natoms,cmin,cmax);

              /* If first time, save the max and min values */
              if (first == TRUE)
                {
                  for (i = 0; i < 3; i++)
                    {
                      coord_min[i] = cmin[i];
                      coord_max[i] = cmax[i];
                    }
                }

              /* Otherwise, update max and min values */
              else
                {
                  for (i = 0; i < 3; i++)
                    {
                      if (cmin[i] < coord_min[i])
                        coord_min[i] = cmin[i];
                      if (cmax[i] > coord_max[i])
                        coord_max[i] = cmax[i];
                    }
                }

              /* Unset flag */
              first = FALSE;
            }

          /* Get the next chain */
          chain_ptr = chain_ptr->next_chain_ptr;
        }

      /* Calculate Reaster3D scaling factor */
      scale_factor = calc_r3d_scale(coord_min,coord_max);

      /* Get pointer to the first chain record yet again */
      chain_ptr = first_chain_ptr;

      /* Loop over the chains to render all the DNA chains */
      while (chain_ptr != NULL)
        {
          /* If this is a DNA/RNA chain, draw as all atoms, plus
             backbone */
          if (chain_ptr->type == DNA)
            {
              /* Render this DNA chain as all atoms plus backbone
                 trace */
              render_dna_chain(file_ptr,chain_ptr,scale_factor);
            }

          /* Get the next chain */
          chain_ptr = chain_ptr->next_chain_ptr;
        }

      /* Close the l.r3d output file */
      fclose(file_ptr);

      /* Run Raster3D to render the image */
      if (params.run_raster3d == TRUE)
        run_raster3d(file_r3d,gif_name,render_exe,convert_exe,
                     ppmquant_exe,ppmtogif_exe,image,FALSE,FALSE);
    }

  /* Return the last PyMol image record */
  *lst_pymol_ptr = last_pymol_ptr;

  /* Return whether still have transparent image */
  return(transparent);
}
/***********************************************************************

render_indiv_chains  -  Render each protein chain separately

***********************************************************************/

int render_indiv_chains(struct chain *first_chain_ptr,
                        struct ligand *first_ligand_ptr,
                        struct atom *first_atom_ptr,int natoms,
                        struct ssbond *first_ssbond_ptr,
                        char *render_exe,char *convert_exe,
                        char *ppmquant_exe,char *ppmtogif_exe,
                        FILE *pymol_script_ptr,
                        struct pobject *first_pobject_ptr,
                        struct pymol **lst_pymol_ptr,int transparent,
                        struct parameters params)
{
  char file_r3d[FILENAME_LEN + 1], gif_name[FILENAME_LEN + 1];
  char chains[3], copy_chain[MAX_CHAINS];
  char image[20];

  int biomol, colour_as_orig, transparent_copies, len, ncopies;

  double coord_min[3], coord_max[3], mass_x, mass_y, mass_z;
  double scale_factor;

  struct chain *chain_ptr, *other_chain_ptr;
  struct pymol *last_pymol_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  biomol = FALSE;
  colour_as_orig = FALSE;
  transparent_copies = FALSE;

  /* Initialise pointers */
  last_pymol_ptr = *lst_pymol_ptr;

  /* Name of temporary Raster3D input file */
  strcpy(file_r3d,"l.r3d");

  /* Get pointer to the first chain record */
  chain_ptr = first_chain_ptr;

  /* Loop over the chains */
  while (chain_ptr != NULL)
    {
      /* Include only if this is a unique protein chain */
      if ((chain_ptr->type == PROTEIN || chain_ptr->type == CALPHA_ONLY)
          && chain_ptr->copy_of == '\0')
        {
          /* Form image name */
          sprintf(image,"Protein chain [%c]",chain_ptr->chain_id);

          /* Form name of the GIF file to be generated */
          if (IMAGE_FORMAT == GIF)
            {
              if (chain_ptr->chain_id == ' ')
                strcpy(gif_name,"chain.gif");
              else
                sprintf(gif_name,"chain%c.gif",chain_ptr->chain_id);
            }
          else
            {
              if (chain_ptr->chain_id == ' ')
                strcpy(gif_name,"chain.jpg");
              else
                sprintf(gif_name,"chain%c.jpg",chain_ptr->chain_id);
            }

          /* For old method, generate images using Raster3D */
          if (pymol_script_ptr == NULL)
            {
              /* Open the l.r3d output file */
              file_ptr = open_output_file(file_r3d,TRUE);

              /* Write out the Raster3D header records */
              write_raster3d_headers(file_ptr);

              /* Get centre of mass of all atoms in this chain */
              centre_of_mass(chain_ptr->first_residue_ptr->first_atom_ptr,
                             chain_ptr->natoms,&mass_x,&mass_y,&mass_z);

              /* Adjust whole structure to centre of mass of this chain */
              adjust_to_origin(first_atom_ptr,natoms,mass_x,mass_y,mass_z,
                               TRUE);

              /* Get maximum and minimum atom coords */
              get_max_min(chain_ptr->first_residue_ptr->first_atom_ptr,
                          chain_ptr->natoms,coord_min,coord_max);

              /* Calculate Reaster3D scaling factor */
              scale_factor = calc_r3d_scale(coord_min,coord_max);

              /* Render according to whether have full protein or only
                 C-alpha coords */
              if (chain_ptr->type == PROTEIN)
                render_protein_chain(file_ptr,chain_ptr,scale_factor,
                                     colour_as_orig,transparent_copies,
                                     biomol);
              else
                render_ca_only_chain(file_ptr,chain_ptr,scale_factor);

              /* Render any disulphide bonds belonging only to the
                 current chain */
              render_ssbonds(file_ptr,first_ssbond_ptr,
                             chain_ptr->chain_id,scale_factor);

              /* Render any of this chain's active site residues */
              render_sites(file_ptr,first_chain_ptr,chain_ptr->chain_id,
                           scale_factor);

              /* Close the l.r3d output file */
              fclose(file_ptr);

              /* Run Raster3D to render the image */
              if (params.run_raster3d == TRUE)
                run_raster3d(file_r3d,gif_name,render_exe,convert_exe,
                             ppmquant_exe,ppmtogif_exe,image,FALSE,
                             FALSE);
            }

          /* If generating a transparent image, determine which chains
             are copies of this one, and generate the transparent image */
          if (transparent == TRUE || pymol_script_ptr != NULL)
            {
              /* Initialise variables */
              ncopies = 0;
              copy_chain[ncopies] = '\0';

              /* Get pointer to the first chain record */
              other_chain_ptr = first_chain_ptr;

              /* Loop over the chains to find copies of the current one */
              while (other_chain_ptr != NULL)
                {
                  /* If a copy, add chain identifier to list */
                  if (other_chain_ptr != chain_ptr &&
                      other_chain_ptr->copy_of == chain_ptr->chain_id)
                    {
                      /* If have space, store chain identifier */
                      if (ncopies < MAX_CHAINS - 2)
                        {
                          /* Store chain identifier */
                          copy_chain[ncopies] = other_chain_ptr->chain_id;

                          /* Increment copy count */
                          ncopies++;
                          copy_chain[ncopies] = '\0';
                        }
                    }

                  /* Get the next chain */
                  other_chain_ptr = other_chain_ptr->next_chain_ptr;
                }

              /* Store chain identifier */
              chains[0] = chain_ptr->chain_id;
              chains[1] = '\0';

              /* Generate the transparent image */
              transparent
                = transparent_image(first_chain_ptr,first_ligand_ptr,
                                    PROTEIN,chains,copy_chain,ncopies,0,
                                    gif_name,render_exe,convert_exe,
                                    ppmquant_exe,ppmtogif_exe,
                                    pymol_script_ptr,first_pobject_ptr,
                                    &last_pymol_ptr,params);
            }
        }

      /* Get the next chain */
      chain_ptr = chain_ptr->next_chain_ptr;
    }

  /* Return the last PyMol image record */
  *lst_pymol_ptr = last_pymol_ptr;

  /* Return whether still have transparent image */
  return(transparent);
}
/***********************************************************************

create_hetgroup_record  -  Create a new hetgroup record or update
                           existing one if we already have it

***********************************************************************/

struct hetgroup *create_hetgroup_record(struct hetgroup **fst_hetgroup_ptr,
                                        struct hetgroup **lst_hetgroup_ptr,
                                        int type,char *name,
                                        int *het_count)
{
  int found;

  struct hetgroup *hetgroup_ptr, *match_hetgroup_ptr;
  struct hetgroup *first_hetgroup_ptr, *last_hetgroup_ptr;

  /* Initialise variables */
  found = FALSE;

  /* Initialise hetgroup pointers */
  last_hetgroup_ptr = *lst_hetgroup_ptr;
  first_hetgroup_ptr = *fst_hetgroup_ptr;
  match_hetgroup_ptr = NULL;

  /* Get pointer to first entry */
  hetgroup_ptr = first_hetgroup_ptr;

  /* Loop through all these to determine whether we already have this
     entry */
  while (hetgroup_ptr != NULL && found == FALSE)
    {
      /* If have match on name and type, then already have this
         entry */
      if (hetgroup_ptr->type == type && !strcmp(hetgroup_ptr->name,name))
        {
          /* Store pointer to located entry */
          match_hetgroup_ptr = hetgroup_ptr;

          /* Increment count of matches */
          hetgroup_ptr->count++;

          /* Set flag that found */
          found = TRUE;
        }

      /* Get next hetgroup entry */
      hetgroup_ptr = hetgroup_ptr->next_hetgroup_ptr;
    }

  /* If match found, then retrieve pointer */
  if (found == TRUE)
    hetgroup_ptr = match_hetgroup_ptr;

  /* Otherwise, create a new hetgroup entry */
  else
    {
      hetgroup_ptr = (struct hetgroup*)malloc(sizeof(struct hetgroup));
      if (hetgroup_ptr == NULL)
        {
          printf("*** ERROR. Unable to allocate memory for");
          printf(" struct hetgroup\n");
          printf("***        Program newsec terminated with error.\n");
          exit (1);
        }

      /* Store the details for this entry */
      hetgroup_ptr->type = type;
      strcpy(hetgroup_ptr->name,name);
      hetgroup_ptr->count = 1;

      /* Initialise pointer to next hetgroup record */
      hetgroup_ptr->next_hetgroup_ptr = NULL;

      /* If this is the first entry stored, then update pointer to
         head of linked list */
      if (first_hetgroup_ptr == NULL)
        first_hetgroup_ptr = hetgroup_ptr;

      /* Otherwise, make previous entry point to the current
         one */
      else
        last_hetgroup_ptr->next_hetgroup_ptr = hetgroup_ptr;
                      
      /* Save the current hetgroup's pointer */
      last_hetgroup_ptr = hetgroup_ptr;

      /* Increment count of this type of hetgroup */
      het_count[type]++;
    }

  /* Return current pointers */
  *fst_hetgroup_ptr = first_hetgroup_ptr;
  *lst_hetgroup_ptr = last_hetgroup_ptr;
  return(hetgroup_ptr);
}
/***********************************************************************

free_hetgroup_memory  -  Free up memory used by residue hetgroup

***********************************************************************/

void free_hetgroup_memory(struct hetgroup *first_hetgroup_ptr)
{
  struct hetgroup *hetgroup_ptr, *next_hetgroup_ptr;

  /* Get first hetgroup */
  hetgroup_ptr = first_hetgroup_ptr;

  /* Loop through the hetgroup records until done */
  while (hetgroup_ptr != NULL)
    {
      /* Get pointer to the next hetgroup record */
      next_hetgroup_ptr = hetgroup_ptr->next_hetgroup_ptr;

      /* Free up memory taken by current hetgroup */
      free(hetgroup_ptr);

      /* Go to the next hetgroup */
      hetgroup_ptr = next_hetgroup_ptr;
    }
}
/***********************************************************************

render_ligands  -  For each ligand type, write out a Raster3D format
                   file and run Raster3D and convert to GIF

***********************************************************************/

void render_ligands(struct chain *first_chain_ptr,
                    struct ligand *first_ligand_ptr,char *render_exe,
                    char *convert_exe,char *ppmquant_exe,
                    char *ppmtogif_exe,char *output_dir,
                    FILE *pymol_script_ptr,
                    struct pobject *first_pobject_ptr,
                    struct pymol **lst_pymol_ptr,int transparent,
                    int images,struct c_file_data file_name_ptr,
                    struct parameters params)
{
  char file_name[FILENAME_LEN + 1], pdb_name[FILENAME_LEN + 1];
  char gif_name[FILENAME_LEN + 1], file_r3d[FILENAME_LEN + 1];
  char pymol_exe[FILENAME_LEN + 1];
  char chains[3], image[20], ltype[6], rastype[10];

  int i, last_metal, last_type, len, ligmet_type, nresid, ntypes, type;
  int dos;

  double coord_min[3], coord_max[3], min_z, scale_factor;
  double transformation[3][3], cofg[3];

  struct ligand *ligand_ptr;
  struct pymol *last_pymol_ptr;
  struct residue *first_residue_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  chains[0] = '\0';
  dos = file_name_ptr.dos;
  last_type = -1;
  last_metal = FALSE;
  ntypes = 0;
  strcpy(pymol_exe,file_name_ptr.exe_name[PYMOL_EXE]);

  /* Initialise pointers */
  last_pymol_ptr = *lst_pymol_ptr;

  /* Get pointer to the first ligand record */
  ligand_ptr = first_ligand_ptr;

  /* Loop over the ligands */
  while (ligand_ptr != NULL)
    {
      /* Determine whether metal or ligand */
      if (ligand_ptr->metal == TRUE)
        {
          /* Form name and description of image */
          strcpy(rastype,"rastermet");
          sprintf(image,"Metal type %d",ligand_ptr->type);
          strcpy(ltype,"met");
        }
      else
        {
          /* Form name and description of image */
          strcpy(rastype,"raster");
          sprintf(image,"Ligand type %d",ligand_ptr->type);
          strcpy(ltype,"lig");
        }

      /* If ligand is of a new type, render this ligand */
      if ((ligand_ptr->type != last_type ||
           ligand_ptr->metal != last_metal) && ligand_ptr->natoms > 0)
        {
          /* Align ligand along principal axes */
          if (ligand_ptr->metal == FALSE)
            {
              /* Get pointer to the first residue record */
              first_residue_ptr = ligand_ptr->first_residue_ptr;
              nresid = ligand_ptr->nresid;

              /* Calculate all the covalent bonds in the ligand */
              if (ligand_ptr->got_bonds == FALSE)
                {
                  calc_covalent_bonds(first_residue_ptr,nresid);
                  ligand_ptr->got_bonds = TRUE;
                }

              /* Use principal components to align ligand along its
                 principal axes */
              align_molecule(ligand_ptr->first_residue_ptr->first_atom_ptr,
                             ligand_ptr->natoms,coord_min,coord_max,
                             transformation,cofg);

              /* Adjust max and min coords for largest atom radius */
              for (i = 0; i < 3; i++)
                {
                  coord_min[i] = coord_min[i] - 2.0;
                  coord_max[i] = coord_max[i] + 2.0;
                }          
            }

          /* Form the name of the ligand's PDB file */
          sprintf(pdb_name,"%s_%2.2d.pdb",ltype,ligand_ptr->type);

          /* Write out the stored coordinates as a PDB file */
          sprintf(file_name,"%s%s",output_dir,pdb_name);
          write_coords(file_name,ligand_ptr->first_residue_ptr,
                       ligand_ptr->nresid,TRUE);

          /* Write out as a .rcm file for use in driving the LIGPLOT
             diagram */
          sprintf(file_name,"%s%s_%2.2d.rcm",output_dir,ltype,
                  ligand_ptr->type);
          write_coords(file_name,ligand_ptr->first_residue_ptr,
                       ligand_ptr->nresid,FALSE);

          /* Form file names */
          if (ligand_ptr->metal == FALSE)
            strcpy(file_r3d,"lig.r3d");
          else
            strcpy(file_r3d,"m.r3d");
          if (IMAGE_FORMAT == GIF)
            sprintf(gif_name,"%s%2.2d.gif",rastype,ligand_ptr->type);
          else
            sprintf(gif_name,"%s%2.2d.jpg",rastype,ligand_ptr->type);

          /* If generating PyMOL version */
          if (params.gen_pymol == TRUE && ligand_ptr->metal == FALSE)
            render_pymol_ligand(pdb_name,ligand_ptr,pymol_exe,dos);

          /* If not metal, render ligand */
          if (ligand_ptr->metal == FALSE && images == GIFS)
            {
              /* Open the .r3d output file */
              file_ptr = open_output_file(file_r3d,TRUE);

              /* Write out the Raster3D header records */
              write_raster3d_headers(file_ptr);

              /* Calculate Raster3D scaling factor */
              scale_factor = calc_r3d_scale(coord_min,coord_max);

              /* Write out ligand to Raster3D file */
              if (ligand_ptr->metal == FALSE)
                min_z = write_raster3d_ligand(file_ptr,ligand_ptr,
                                              scale_factor);

              /* Place the back wall just behind the molecule */
              min_z = min_z - scale_factor * BACK_WALL_GAP;
              // if (ligand_ptr->metal == FALSE)
              //  {
              // fprintf(file_ptr,"6\n");
              // fprintf(file_ptr,"0.0 0.0 %6.3f 0.2 0.0 %6.3f 0.2 0.2 %6.3f\n",
              //         min_z,min_z,min_z);
              // fprintf(file_ptr,"1.0 1.0 1.0\n");
              //  }

              /* Close the .r3d output file */
              fclose(file_ptr);

              /* Run Raster3D to render the image */
              if (params.run_raster3d == TRUE &&
                  ligand_ptr->metal == FALSE)
                run_raster3d(file_r3d,gif_name,render_exe,convert_exe,
                             ppmquant_exe,ppmtogif_exe,image,transparent,
                             FALSE);
            }

          /* If need to generate transparent image, do so now */
          if ((transparent == TRUE || pymol_script_ptr != NULL) &&
              images == GIFS)
            {
              /* Determine object type */
              if (ligand_ptr->metal == TRUE)
                type = METAL;
              else
                type = LIGAND;
              ligmet_type = ligand_ptr->type;

              /* Generate the transparent image */
              transparent
                = transparent_image(first_chain_ptr,first_ligand_ptr,
                                    type,chains,"0",0,ligmet_type,
                                    gif_name,render_exe,convert_exe,
                                    ppmquant_exe,ppmtogif_exe,
                                    pymol_script_ptr,first_pobject_ptr,
                                    &last_pymol_ptr,params);
            }

          /* Save the current ligand type and whether metal or not */
          last_type = ligand_ptr->type;
          last_metal = ligand_ptr->metal;
        }

      /* Get the next ligand */
      ligand_ptr = ligand_ptr->next_ligand_ptr;
    }

  /* Return the last PyMol image record */
  *lst_pymol_ptr = last_pymol_ptr;
}
/***********************************************************************

render_chain_pair  -  Render the given protein chain-pair

***********************************************************************/

void render_chain_pair(struct chain *first_chain_ptr,char *chains,
                       struct atom *first_atom_ptr,int natoms,
                       char *gif_name,char *render_exe,char *convert_exe,
                       char *ppmquant_exe,char *ppmtogif_exe,
                       FILE *pymol_script_ptr,
                       struct pobject *first_pobject_ptr,
                       struct pymol **lst_pymol_ptr,
                       struct parameters params)
{
  char file_r3d[FILENAME_LEN + 1];
  char image[50];

  int biomol, colour_as_orig, transparent_copies;

  double coord_min[3], coord_max[3], mass_x, mass_y, mass_z;
  double scale_factor;

  struct chain *chain_ptr, *chains_ptr[2];
  struct pymol *last_pymol_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  biomol = FALSE;
  chains_ptr[0] = NULL;
  chains_ptr[1] = NULL;
  colour_as_orig = FALSE;
  transparent_copies = FALSE;

  /* Initialise pointers */
  last_pymol_ptr = *lst_pymol_ptr;

  /* Name of temporary Raster3D input file */
  strcpy(file_r3d,"i.r3d");

  /* Get pointer to the first chain record */
  chain_ptr = first_chain_ptr;

  /* Loop to identify the two chains in question */
  while (chain_ptr != NULL)
    {
      /* Check for first chain */
      if (chain_ptr->type == PROTEIN && chain_ptr->chain_id == chains[0])
        {
          /* Have the first chain */
          chains_ptr[0] = chain_ptr;
        }

      /* Check for second chain */
      if (chain_ptr->type == PROTEIN && chain_ptr->chain_id == chains[1])
        {
          /* Have the second chain */
          chains_ptr[1] = chain_ptr;
        }

      /* Get the next chain */
      chain_ptr = chain_ptr->next_chain_ptr;
    }

  /* If have both chains, then continue */
  if (chains_ptr[0] != NULL && chains_ptr[1] != NULL)
    {
      /* Open the i.r3d output file */
      file_ptr = open_output_file(file_r3d,TRUE);

      /* Write out the Raster3D header records */
      write_raster3d_headers(file_ptr);

      /* Form image name */
      sprintf(image,"Protein chains %c and %c",chains[0],chains[1]);

      /* Get centre of mass of all atoms in this chain */
      centre_of_mass2(chains_ptr[0]->first_residue_ptr->first_atom_ptr,
                      chains_ptr[0]->natoms,
                      chains_ptr[1]->first_residue_ptr->first_atom_ptr,
                      chains_ptr[1]->natoms,&mass_x,&mass_y,&mass_z);

      /* Adjust whole structure to centre of mass of this chain */
      adjust_to_origin(first_atom_ptr,natoms,mass_x,mass_y,mass_z,
                       TRUE);

      /* Get maximum and minimum atom coords */
      get_max_min2(chains_ptr[0]->first_residue_ptr->first_atom_ptr,
                   chains_ptr[0]->natoms,
                   chains_ptr[1]->first_residue_ptr->first_atom_ptr,
                   chains_ptr[1]->natoms,coord_min,coord_max);

      /* Calculate Reaster3D scaling factor */
      scale_factor = calc_r3d_scale(coord_min,coord_max);

      /* Render the two protein chains */
      render_protein_chain(file_ptr,chains_ptr[0],scale_factor,
                           colour_as_orig,transparent_copies,biomol);
      render_protein_chain(file_ptr,chains_ptr[1],scale_factor,
                           colour_as_orig,transparent_copies,biomol);

      /* Close the i.r3d output file */
      fclose(file_ptr);

      /* Run Raster3D to render the image */
      if (params.run_raster3d == TRUE)
        run_raster3d(file_r3d,gif_name,render_exe,convert_exe,
                     ppmquant_exe,ppmtogif_exe,image,FALSE,
                     FALSE);
    }

  /* Return the last PyMol image record */
  *lst_pymol_ptr = last_pymol_ptr;
}
/***********************************************************************

render_interfaces  -  Render each protein chain separately

***********************************************************************/

void render_interfaces(struct interface *first_interface_ptr,
                       struct chain *first_chain_ptr,
                       struct ligand *first_ligand_ptr,
                       struct atom *first_atom_ptr,int natoms,
                       char *render_exe,char *convert_exe,
                       char *ppmquant_exe,char *ppmtogif_exe,
                       FILE *pymol_script_ptr,
                       struct pobject *first_pobject_ptr,
                       struct pymol **lst_pymol_ptr,
                       struct parameters params)
{
  char gif_name[FILENAME_LEN + 1], tgif_name[FILENAME_LEN + 1];
  char copy_interface[2 * MAX_INTERFACES];

  int len, ncopies;

  struct interface *interface_ptr, *other_interface_ptr;
  struct pymol *last_pymol_ptr;

  /* Initialise pointers */
  last_pymol_ptr = *lst_pymol_ptr;

  /* Get pointer to the first interface record */
  interface_ptr = first_interface_ptr;

  /* Loop over the interfaces */
  while (interface_ptr != NULL)
    {
      /* Include only if this is a unique protein-protein interface
         (protein-DNA interfaces aren't currently plotted) */
      if (interface_ptr->type == PROTEIN_PROTEIN)
        {
/* debug */
          printf("INTERFACE: Chains %s\n",interface_ptr->chains);
          fflush(stdout);
/* debug */
          /* Form name of output image files */
          if (IMAGE_FORMAT == GIF)
            {
              /* File name for chains image and transparent image */
              sprintf(gif_name,"chains%c%c.gif",interface_ptr->chains[0],
                      interface_ptr->chains[1]);
              sprintf(tgif_name,"ppiface%c%c.gif",interface_ptr->chains[0],
                      interface_ptr->chains[1]);
            }
          else
            {
              /* File name for chains image and transparent image */
              sprintf(gif_name,"chains%c%c.jpg",interface_ptr->chains[0],
                      interface_ptr->chains[1]);
              sprintf(tgif_name,"ppiface%c%c.jpg",interface_ptr->chains[0],
                      interface_ptr->chains[1]);
            }

          /* Render the pair of chains */
          if (pymol_script_ptr == NULL)
            render_chain_pair(first_chain_ptr,interface_ptr->chains,
                              first_atom_ptr,natoms,gif_name,
                              render_exe,convert_exe,ppmquant_exe,
                              ppmtogif_exe,pymol_script_ptr,
                              first_pobject_ptr,&last_pymol_ptr,params);

          /* Initialise variables */
          ncopies = 0;
          copy_interface[ncopies] = '\0';

          /* Get pointer to the next interface record */
          other_interface_ptr = interface_ptr->next_interface_ptr;

          /* Loop over the interfaces to find copies of the current one */
          while (other_interface_ptr != NULL)
            {
              /* If a copy, add interface identifier to list */
              if (other_interface_ptr->copyof_interface_ptr
                  == interface_ptr)
                {
                  /* If have space, store interface identifier */
                  if (ncopies < 2 * MAX_INTERFACES - 2)
                    {
                      /* Store the two chain identifiers */
                      copy_interface[ncopies]
                        = other_interface_ptr->chains[0];
                      ncopies++;
                      copy_interface[ncopies]
                        = other_interface_ptr->chains[1];
                      ncopies++;

                      /* Terminate copy chains string */
                      copy_interface[ncopies] = '\0';
                    }
                }

              /* Get the next interface */
              other_interface_ptr
                = other_interface_ptr->next_interface_ptr;
            }

          /* Generate the transparent image */
          transparent_image(first_chain_ptr,first_ligand_ptr,PP_INTERFACE,
                            interface_ptr->chains,copy_interface,ncopies,
                            0,tgif_name,render_exe,convert_exe,
                            ppmquant_exe,ppmtogif_exe,pymol_script_ptr,
                            first_pobject_ptr,&last_pymol_ptr,params);
        }

      /* Get the next interface */
      interface_ptr = interface_ptr->next_interface_ptr;
    }
/* debug */
  printf("All interfaces done\n");
  fflush(stdout);
/* debug */

  /* Return the last PyMol image record */
  *lst_pymol_ptr = last_pymol_ptr;
}
/***********************************************************************

write_dat_file  -  Write out the ligands.dat file detailing all the
                   ligands and their duplicates

***********************************************************************/

void write_dat_file(char *file_name,struct ligand *first_ligand_ptr,
                    int ntypes,int nmetypes,
                    struct hetdesc *first_hetdesc_ptr,
                    char stored_metal[NMETALS][3],int *metal_count,
                    int nmetyp,int transparent)
{
  char item1[10], item2[10], item3[10];
  char element[3], lighet[5], ligtype[8], ltype[6];

  int colour, i, last_metal, last_type, ncopy, nlines, metal_type, nmetals;

  struct atom *atom_ptr;
  struct ligand *ligand_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  last_type = -1;
  last_metal = FALSE;
  nlines = 0;
  nmetals = 0;

  /* Re-initialise metal counts */
  for (i = 0; i < NMETALS; i++)
    metal_count[i] = 0;

  /* Open the output ligands.dat file */
  file_ptr = open_output_file(file_name,TRUE);

  /* Write out transparency flag, if relevant */
  if (transparent == TRUE)
    fprintf(file_ptr,"TRANSPARENT TRUE\n");

  /* Write out the number of ligand types */
  if (ntypes > 0)
    {
      fprintf(file_ptr,"LIGAND TRUE\n");
      fprintf(file_ptr,"NLIGAND_TYPES %d\n",ntypes);
    }
  if (nmetypes > 0)
    {
      fprintf(file_ptr,"METAL TRUE\n");
      fprintf(file_ptr,"NMETAL_TYPES %d\n",nmetypes);
    }

  /* Get pointer to the first ligand record */
  ligand_ptr = first_ligand_ptr;

  /* Loop over the ligands */
  while (ligand_ptr != NULL)
    {
      /* If ligand is of a new type, write the residue list */
      if (ligand_ptr->type != last_type ||
          ligand_ptr->metal != last_metal)
        {
          /* Determine whether metal or ligand */
          if (ligand_ptr->metal == TRUE)
            {
              strcpy(ligtype,"METAL");
              strcpy(ltype,"MET");
              strcpy(lighet,"MET");
            }
          else
            {
              strcpy(ligtype,"LIGAND");
              strcpy(ltype,"LIG");
              strcpy(lighet,"HET");
            }

          /* Form the subscript for this type */
          sprintf(item1,"[%d]",ligand_ptr->type);

          /* Write out the residue type */
          fprintf(file_ptr,"%s_TYPE%s %2.2d\n",ligtype,item1,
                  ligand_ptr->type);

          /* Write out the full ligand sequence */
          fprintf(file_ptr,"%s_SEQ%s %s\n",ligtype,item1,
                  ligand_ptr->residue_sequence);

          /* Write out the number of residues */
          fprintf(file_ptr,"%s_NRESID%s %d\n",ligtype,item1,
                  ligand_ptr->nresid);

          /* If a metal, write out the colour */
          if (ligand_ptr->metal == TRUE)
            {
              /* Get the atom corresponding to this metal ion */
              atom_ptr
                = ligand_ptr->first_residue_ptr->first_atom_ptr;

              /* Get the element */
              check_for_metal(atom_ptr->atom_name,element);

              /* Store the metal name */
              strcpy(ligand_ptr->metal_name,element);

              /* Get the corresponding colour */
              colour = -1;
              for (metal_type = 0; metal_type < nmetyp && colour == -1;
                   metal_type++)
                {
                  /* If element name matches, store colour */
                  if (!strcmp(stored_metal[metal_type],element))
                    {
                      /* Retrieve the metal's colour */
                      colour = metal_type % NMETAL_COLOURS;

                      /* Update count of metal ions in structure */
                      metal_count[metal_type]
                        = metal_count[metal_type] + ligand_ptr->count;
                    }
                }

              /* Write out the colour */
              if (colour < 0)
                colour = 0;
              fprintf(file_ptr,"METAL_COLOUR%s %s\n",item1,
                      METAL_COLOUR[colour]);

              /* Store the metal's colour */
              ligand_ptr->metal_colour = colour;
            }

          /* Split the ligand's residue list into smaller chunks */
          // nlines
          //   = split_residue_sequence(residue_sequence,
          //                            ligand_ptr->residue_sequence);

          /* Write out the number of lines making up the ligand name */
          // fprintf(file_ptr,"NLIGLINE%s %d\n",item1,nlines);

          /* Loop over the chunks to write out */
          // for (iline = 0; iline < nlines; iline++)
          //   {
          //     /* Form the item number for this line */
          //     sprintf(item2,"[%d]",iline + 1);

          //     /* Write out the line */
          //     fprintf(file_ptr,"LIGAND_NAME%s%s %s\n",item1,item2,
          //             residue_sequence[iline]);
          //   }

          /* Write out the number of occurrences of this ligand */
          fprintf(file_ptr,"N%s_TIMES%s %d\n",ligtype,item1,
                  ligand_ptr->count);

          /* Write out the count of unique copies of this ligand */
          ncopy = ligand_ptr->ncopy + 1;
          if (ncopy > MAX_COPIES)
            {
              /* If more than the maximum allowed, show the maximum
                 and relegate the rest to "... plus others" */
              fprintf(file_ptr,"N%sUNIQ%s %d\n",ltype,item1,MAX_COPIES);
              fprintf(file_ptr,"N%sMORE%s %d\n",ltype,item1,
                      ncopy - MAX_COPIES);
            }
          else
            fprintf(file_ptr,"N%sUNIQ%s %d\n",ltype,item1,ncopy);

          /* Write out any Het Groups and their descriptions */
          write_het_groups(file_ptr,ligand_ptr->first_residue_ptr,
                           ligand_ptr->nresid,first_hetdesc_ptr,lighet,
                           item1);

          /* Save the current ligand type and whether metal or not */
          last_type = ligand_ptr->type;
          last_metal = ligand_ptr->metal;
        }

      /* If this is a unique ligand, write out residue range */
      if (ligand_ptr->unique == TRUE)
        {
          /* Form the subscript for this copy */
          sprintf(item2,"[%d]",ligand_ptr->unique_count+1);

          /* Write out the residue range */
          fprintf(file_ptr,"%s_RANGE%s%s %s\n",ltype,item1,item2,
                  ligand_ptr->residue_range);

          /* Write out the number of duplicates, if any */
          if (ligand_ptr->nduplicates > 0)
            fprintf(file_ptr,"N%s_DUPLIC%s%s %d\n",ltype,item1,item2,
                    ligand_ptr->nduplicates);
        }

      /* Otherwise, write out as a copy */
      else
        {
          /* Form the subscript for this duplicate */
          sprintf(item3,"[%d]",ligand_ptr->duplicate_count);

          /* Write out the residue range */
          fprintf(file_ptr,"%s_DUPLIC%s%s%s %s\n",ltype,item1,item2,item3,
                  ligand_ptr->residue_range);

        }

      /* Get the next ligand */
      ligand_ptr = ligand_ptr->next_ligand_ptr;
    }

  /* Loop over the metals, writing out names and colours */
  // for (metal_type = 0; metal_type < nmetyp; metal_type++)
  //   {
  //     /* If have at least one instance of this metal, write out
  //        the metal record */
  //     if (metal_count[metal_type] > 0)
  //       {
  //         /* Increment count of metal types written out */
  //         nmetals++;

  //         /* Get the metal colour */
  //         colour = metal_type % NMETAL_COLOURS;

  //         /* Write out this metal type */
  //         fprintf(file_ptr,"METAL_ID[%d] %s\n",nmetals,
  //                 stored_metal[metal_type]);
  //         fprintf(file_ptr,"METAL_COL[%d] %s\n",nmetals,
  //                 METAL_COLOUR[colour]);
  //         fprintf(file_ptr,"METAL_NUM[%d] %d\n",nmetals,
  //                 metal_count[metal_type]);
  //       }
  //   }

  /* If had any metals, write out the number of types */
  // if (nmetals > 0)
  //   {
  //     fprintf(file_ptr,"METAL TRUE\n");
  //     fprintf(file_ptr,"NMETAL %d\n",nmetals);
  //   }

  /* Close the output file */
  fclose(file_ptr);
}
/***********************************************************************

write_secsurf_sequence  -  Write out the sequence to the secsurf.out
                           file

***********************************************************************/

void write_secsurf_sequence(FILE *file_ptr,
                            struct residue *first_residue_ptr,int nresid,
                            int type)
{
  char res_code[15];
 
  int het_count[NHETYPE], het_type, i, ipos, ires, iresid, len;

  struct hetgroup *hetgroup_ptr, *first_hetgroup_ptr, *last_hetgroup_ptr;
  struct residue *residue_ptr, *attach_residue_ptr;

  /* Initialise variables */
  ipos = 0;
  ires = 0;
  for (i = 0; i < NHETYPE; i++)
    het_count[i] = 0;

  /* Initialise pointers */
  hetgroup_ptr = NULL;
  first_hetgroup_ptr = NULL;
  last_hetgroup_ptr = NULL;

  /* Get the first residues */
  residue_ptr = first_residue_ptr;
  iresid = 0;

  /* Loop through all the residues forming the output string to
     write out */
  while (residue_ptr != NULL && iresid < nresid)
    {
      /* If this is a standard amino acid code, or DNA nucleotide,
         then use single letter code */
      if (residue_ptr->type == STANDARD ||
          residue_ptr->type == NON_STANDARD ||
          residue_ptr->type == NUCLEIC_ACID)
        {
          /* Get the single-letter code unless in ligand */
          if (type != LIGAND)
            {
              res_code[0] = residue_ptr->aa_code;
              res_code[1] = '\0';
            }
          else
            strcpy(res_code,residue_ptr->res_name);
          het_type = -1;

          /* If residue is a ligand residue, then consider it special
             that have a standard amino acid code here */
          if (type == LIGAND && residue_ptr->type == STANDARD)
            het_type = AMINO_ACID;

          /* If non-standard residue in protein chain, then mark
             as unusual */
          if (residue_ptr->type == NON_STANDARD)
            het_type = NON_STANDARD_AA;

          /* Create a hetgroup record, if required */
          if (het_type > -1)
            hetgroup_ptr
              = create_hetgroup_record(&first_hetgroup_ptr,
                                       &last_hetgroup_ptr,
                                       het_type,residue_ptr->res_name,
                                       het_count);
        }

      /* Otherwise use the three letter name */
      else
        {
          /* Copy name */
          strcpy(res_code,residue_ptr->res_name);

          /* Create a hetgroup record */
          hetgroup_ptr
            = create_hetgroup_record(&first_hetgroup_ptr,
                                     &last_hetgroup_ptr,
                                     HET_GROUP,residue_ptr->res_name,
                                     het_count);
        }

      /* Check whether this residue has an attachment */
      if (residue_ptr->attach_residue_ptr != NULL)
        {
          /* Get pointer to the attaching residue */
          attach_residue_ptr = residue_ptr->attach_residue_ptr;

          /* Create a hetgroup record for this attachment */
          hetgroup_ptr
            = create_hetgroup_record(&first_hetgroup_ptr,
                                     &last_hetgroup_ptr,
                                     RESIDUE_ATTACHMENT,
                                     attach_residue_ptr->res_name,
                                     het_count);

          /* Add residue name to current residue name as (+XXX) */
          strcat(res_code,"(+");
          strcat(res_code,attach_residue_ptr->res_name);
          strcat(res_code,")");
        }

      /* Replace any blanks with underscores */
      len = strlen(res_code);
      for (i = 0; i < len; i++)
        if (res_code[i] == ' ')
          res_code[i] = '_';

      /* Check whether too close to end */
      if (ipos + len + 1 > MAX_SEQ_LINE - 1)
        {
          /* Close off current line */
          fprintf(file_ptr,"\n");

          /* Re-initialise variables */
          ipos = 0;
          ires = 0;
        }

      /* If at start of the line, write out an asterisk */
      if (ipos == 0)
        {
          fprintf(file_ptr,"*");
          ipos++;
        }

      /* If not first residue on line, precede by hyphen */
      if (ires > 0)
        {
          fprintf(file_ptr,"-");
          ipos++;
        }

      /* Write out the current residue code */
      fprintf(file_ptr,"%s",res_code);
      ires++;
      ipos = ipos + len;

      /* Get the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;
    }

  /* If have a line to close off, then do so */
  if (ires > 0)
    fprintf(file_ptr,"\n");

  /* Gt pointer to the first of the hetgroups in this chain, if any */
  hetgroup_ptr = first_hetgroup_ptr;

  /* Write out any unusual residues: hetgroups in protein or DNA chains,
     amino acids in ligands, and residue attachments */
  while (hetgroup_ptr != NULL)
    {
      /* Write out type */
      if (hetgroup_ptr->type == HET_GROUP)
        fprintf(file_ptr,"HET group:-            %s%6d\n",
                hetgroup_ptr->name,hetgroup_ptr->count);
      else if (hetgroup_ptr->type == AMINO_ACID)
        fprintf(file_ptr,"Amino acid:-           %s%6d\n",
                hetgroup_ptr->name,hetgroup_ptr->count);
      else if (hetgroup_ptr->type == NON_STANDARD_AA)
        fprintf(file_ptr,"Non-standard residue:- %s%6d\n",
                hetgroup_ptr->name,hetgroup_ptr->count);
      else if (hetgroup_ptr->type == RESIDUE_ATTACHMENT)
        fprintf(file_ptr,"Attached HET group:-   %s%6d\n",
                hetgroup_ptr->name,hetgroup_ptr->count);

      /* Get the next het group */
      hetgroup_ptr = hetgroup_ptr->next_hetgroup_ptr;
    }

  /* Free up memory taken by het group pointers */
  free_hetgroup_memory(first_hetgroup_ptr);
}
/***********************************************************************

write_secsurf_chain  -  Write out the detailed information on each
                          protein chain

***********************************************************************/

void write_secsurf_chain(FILE *file_ptr,struct chain *first_chain_ptr)
{
#define NCTYPES    3

  char ca_only[8], colour_name[COLNAM_LEN];
  char colour_helix[COLNAM_LEN + 4],  colour_strand[COLNAM_LEN + 4];
  char colour_coil[COLNAM_LEN + 5];

  static char *chain_desc[NCTYPES]
    = { "CHAIN", "NUCLE", "LIGAN" };

  int len, loop, nchain;
  int wanted;

  static int chain_type[NCTYPES] = { PROTEIN, DNA, LIGAND };
  static int altch_type[NCTYPES] = { CALPHA_ONLY, DNA, LIGAND };

  struct chain *chain_ptr;

  /* Initialise variables */
  nchain = 0;

  /* Loop over the different chain types to be written out */
  for (loop = 0; loop < NCTYPES; loop++)
    {
      /* Re-initialise chain counter */
      nchain = 0;

      /* Get pointer to the first chain */
      chain_ptr = first_chain_ptr;

      /* Loop over all the chains */
      while (chain_ptr != NULL)
        {
          /* Check if this is the right chain type */
          wanted = FALSE;
          if (chain_ptr->type == chain_type[loop] ||
              chain_ptr->type == altch_type[loop])
            wanted = TRUE;

          /* For ligand, only take the originasl copy */
          if (chain_ptr->type == LIGAND && chain_ptr->copy_of != '\0')
            wanted = FALSE;

          /* If appropriate chain type, then write out details */
          if (wanted == TRUE)
            {
              /* Increment chain count */
              nchain++;

              /* Get this chain's colour */
              get_chain_colname(chain_ptr->chain_id,colour_name,
                                colour_coil,colour_helix,colour_strand);

              /* Right-pad with blanks */
              len = strlen(colour_name);
              if (len < 11)
                strcat(colour_name,"         ");
              colour_name[10] = '\0';

              /* If C-alpha only, set marker */
              if (chain_ptr->type == CALPHA_ONLY)
                strcpy(ca_only,"CA only");
              else
                strcpy(ca_only,"       ");

              /* Write out chain name and summary details */
              fprintf(file_ptr,"%s%3d [%c] %8d%8d   ",
                      chain_desc[loop],nchain,chain_ptr->chain_id,
                      chain_ptr->natoms,chain_ptr->nresid);

              /* For protein/DNA write chain colour */
              if (chain_ptr->type == PROTEIN ||
                  chain_ptr->type == CALPHA_ONLY ||
                  chain_ptr->type == DNA)
                fprintf(file_ptr,"%s %s ",ca_only,colour_name);
              else
                fprintf(file_ptr,"Copies:%4d",chain_ptr->ncopies);

              /* For protein, write out the number of helices and
                 strands */
              if (chain_ptr->type == PROTEIN)
                fprintf(file_ptr,"Helices:%4d  Strands:%4d",
                        chain_ptr->nhelices,chain_ptr->nstrands);
              fprintf(file_ptr,"\n");

              /* Write out the sequence, including any residue attachments */
              write_secsurf_sequence(file_ptr,chain_ptr->first_residue_ptr,
                                     chain_ptr->nresid,chain_ptr->type);
            }

          /* Get the next chain */
          chain_ptr = chain_ptr->next_chain_ptr;
        }
    }
}
/***********************************************************************

write_surf_defs  -  Write out the surfnet.par definitions for each
                    protein or DNA chain

***********************************************************************/

void write_surf_defs(FILE *file_ptr,struct chain *first_chain_ptr)
{
  char colour_name[COLNAM_LEN], stype;
  char colour_helix[COLNAM_LEN + 4],  colour_strand[COLNAM_LEN + 4];
  char colour_coil[COLNAM_LEN + 5];

  int icol, surfnet_col;

  struct chain *chain_ptr;

  /* Get pointer to the first chain */
  chain_ptr = first_chain_ptr;

  /* Loop over all the chains */
  while (chain_ptr != NULL)
    {
      /* If chain is protein or DNA then write out the atom range */
      if (chain_ptr->type == PROTEIN || chain_ptr->type == CALPHA_ONLY ||
          chain_ptr->type == DNA)
        {
          /* For C-alpha only chains, write out as backbone */
          if (chain_ptr->type == CALPHA_ONLY)
            stype = 'J';
          else
            stype = 'S';

          /* Get the chain colour */
          get_chain_colname(chain_ptr->chain_id,colour_name,
                            colour_coil,colour_helix,colour_strand);

          /* Determine what the corresponding SURFNET colour number
             is */
          surfnet_col = 0;
          for (icol = 0; icol < NSURFNET_COLS && surfnet_col == 0; icol++)
            {
              if (!strcmp(colour_name,SURFNET_COLOUR[icol]))
                surfnet_col = icol + 1;
            }

          /* Write out the surfnet output line */
          fprintf(file_ptr,"surface%c      %c %8d%8d  100.0",
                  chain_ptr->chain_id,stype,
                  chain_ptr->first_atom_number,
                  chain_ptr->last_atom_number);
          fprintf(file_ptr,"      %2d        10        0.1  SOLID\n",
                  surfnet_col);
        }

      /* Get the next chain */
      chain_ptr = chain_ptr->next_chain_ptr;
    }
}
/***********************************************************************

write_secsurf_out  -  Write out the secsurf.out file

***********************************************************************/

void write_secsurf_out(struct chain *first_chain_ptr,
                       struct ligand *first_ligand_ptr,int nresid,
                       int natoms,char stored_metal[NMETALS][3],
                       int *metal_count,int nmetyp,struct data data_item,
                       char *output_dir)
{
  char file_name[FILENAME_LEN + 1], metal_name[5];

  int i, itype, metal_type;
  int found;

  struct ligand *ligand_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  for (i = 0; i < NMETALS; i++)
    metal_count[i] = 0;
  
  /* Form output file names */
  strcpy(file_name,output_dir);
  strcat(file_name,"secsurf.out");

  /* Open the secsurf.out output file */
  file_ptr = open_output_file(file_name,TRUE);

  /* Write out the header details */
  fprintf(file_ptr,"PDB code: %s\n",data_item.pdb_code);
  fprintf(file_ptr,"Resolution:  %6.3f\n",data_item.resolution);
  fprintf(file_ptr,"R-factor:    %6.3f\n",data_item.rfactor);
  fprintf(file_ptr,"R-free:      %6.3f\n",data_item.rfree);
  if (data_item.nmr == TRUE)
    {
      fprintf(file_ptr,"NMR structure: %6d model",data_item.nmodels);
      if (data_item.nmodels > 1)
        fprintf(file_ptr,"s\n");
      else
        fprintf(file_ptr,"\n");
    }
  if (data_item.ensemble == TRUE)
    {
      fprintf(file_ptr,"Ensemble of %6d model",data_item.nmodels);
      if (data_item.nmodels > 1)
        fprintf(file_ptr,"s\n");
      else
        fprintf(file_ptr,"\n");
    }
  fprintf(file_ptr,"\n");

  /* Print counts of the different chain types and atoms read in */
  fprintf(file_ptr,
          "Number of chains:             %8d\n",data_item.nchains);
  fprintf(file_ptr,
          "Number of protein chains      %8d\n",data_item.nprotein);
  fprintf(file_ptr,
          "Number of nucleic acid chains %8d\n",data_item.ndna);
  fprintf(file_ptr,
          "Number of C-alpha only chains %8d\n",data_item.nca_only);
  fprintf(file_ptr,
          "Number of residues:           %8d\n",data_item.nresidues);
  fprintf(file_ptr,
          "Number of ligands:            %8d\n",data_item.nligands);
  fprintf(file_ptr,
          "Number of metal ions:         %8d\n",data_item.nmetals);
  fprintf(file_ptr,
          "Number of attachments:        %8d\n",data_item.nattach);
  fprintf(file_ptr,
          "Number of unknowns:           %8d\n",data_item.nunknown);
  fprintf(file_ptr,
          "Number of waters:             %8d\n",data_item.nwaters);
  fprintf(file_ptr,
          "Number of ATOM records:       %8d\n",data_item.natom_records);
  fprintf(file_ptr,
          "Number of HETATM records:     %8d\n",data_item.nhetatm_records);
  fprintf(file_ptr,
          "Number of hydrogens:          %8d\n",data_item.nhydrogens);
  fprintf(file_ptr,"\n");

  /* Write info on each chain */
  write_secsurf_chain(file_ptr,first_chain_ptr);

  /* Get pointer to the first ligand record again */
  ligand_ptr = first_ligand_ptr;

  /* Loop over the ligand records to identify and count */
  while (ligand_ptr != NULL)
    {
      /* If metal, then process */
      if (ligand_ptr->metal == TRUE)
        {
          /* Get this metal */
          strcpy(metal_name,ligand_ptr->metal_name);

          /* Initialise flag */
          found = FALSE;

          /* Loop over the stored metal names */
          for (itype = 0; itype < nmetyp && found == FALSE; itype++)
            {
              /* If have a match, update the count */
              if (!strcmp(metal_name,stored_metal[itype]))
                {
                  /* Increment count */
                  metal_count[itype]++;

                  /* Set flag */
                  found = TRUE;
                }
            }
        }

      /* Get the next ligand */
      ligand_ptr = ligand_ptr->next_ligand_ptr;
    }

  /* Write out the counts of all the metals in the structure */
  for (metal_type = 0; metal_type < nmetyp; metal_type++)
    {
      /* If have at least one instance of this metal, write out
         the metal record */
      if (metal_count[metal_type] > 0)
        {
          /* Write out this metal type */
          fprintf(file_ptr,"METAL %s  %6d\n",stored_metal[metal_type],
                  metal_count[metal_type]);
        }
    }

  /* Write out the surfnet.par-type definitions for the surfaces */
  write_surf_defs(file_ptr,first_chain_ptr);

  /* Close the secsurf.out file */
  fclose(file_ptr);
}
/***********************************************************************

write_surfp_par  -  Write out the surfp.par file

***********************************************************************/

void write_surfp_par(struct chain *first_chain_ptr,char *output_dir,
                     char *pdb_file,double grid_separation)
{
  char file_name[FILENAME_LEN + 1];

  FILE *file_ptr;

  /* Form output file name */
  strcpy(file_name,output_dir);
  strcat(file_name,"surfp.par");

  /* Open the surfp.par output file */
  file_ptr = open_output_file(file_name,TRUE);

  /* Write out the header records */
  fprintf(file_ptr,"Parameters for SURFNET              (surfnet.par)\n");
  fprintf(file_ptr,"----------------------\n");
  fprintf(file_ptr,"TITLE\n");
  fprintf(file_ptr,"\n");
  fprintf(file_ptr,"\n");
  fprintf(file_ptr,"INPUT FILE\n");
  fprintf(file_ptr,"%s\n",pdb_file);
  fprintf(file_ptr,"\n");
  fprintf(file_ptr,"OUTPUT FILES\n");

  /* Write out the surfnet.par-type definitions for the surfaces */
  write_surf_defs(file_ptr,first_chain_ptr);

  /* Write out the first few SURFNET parameters to get the grid-spacing
     in */
  fprintf(file_ptr,"\n");
  fprintf(file_ptr,"\n");
  fprintf(file_ptr,"Program parameters\n");
  fprintf(file_ptr,"------------------\n");
  fprintf(file_ptr,"\n");
  fprintf(file_ptr,"SURFNET\n");
  fprintf(file_ptr,"  0.0   0.0   0.0\n");
  fprintf(file_ptr,"-1000 500000\n");
  fprintf(file_ptr,"2.0 2.0 2.0\n");
  fprintf(file_ptr,"%6.2f\n",grid_separation);

  /* Close the secsurf.out file */
  fclose(file_ptr);
}
/***********************************************************************

write_rasmol_dfn  -  Write out the rasmol.dfn file

***********************************************************************/

void write_rasmol_dfn(struct chain *first_chain_ptr,
                      struct ligand *first_ligand_ptr,int nresid,
                      int natoms,char stored_metal[NMETALS][3],
                      int *metal_count,int nmetyp,struct data data_item,
                      char *output_dir)
{
  char file_name[FILENAME_LEN + 1], last_sequence[10];
  char ctype;

  static char ctypes[3] = { 'C', 'P', 'N' };

  int colour, loop, metal_type, nligands;

  struct chain *chain_ptr;
  struct ligand *ligand_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  nligands = 0;
  
  /* Form output file names */
  strcpy(file_name,output_dir);
  strcat(file_name,"rasmol.dfn");

  /* Open the rasmol.dfn output file */
  file_ptr = open_output_file(file_name,TRUE);

  /* Loop over protein and DNA chains */
  for (loop = 0; loop < 3; loop++)
    {
      /* Get chain type marker */
      ctype = ctypes[loop];

      /* Get pointer to the first chain */
      chain_ptr = first_chain_ptr;

      /* Loop over all the chains */
      while (chain_ptr != NULL)
        {
          /* If appropriate chain type, then write out details */
          if ((loop == 0 && chain_ptr->type == CALPHA_ONLY) ||
              (loop == 1 && chain_ptr->type == PROTEIN) ||
              (loop == 2 && chain_ptr->type == DNA))
            {
              /* Write out this chain type and chain id */
              fprintf(file_ptr,"%c%c\n",ctype,chain_ptr->chain_id);
            }

          /* Get the next chain */
          chain_ptr = chain_ptr->next_chain_ptr;
        }
    }

  /* Loop over the metals, writing out names and colours */
  // for (metal_type = 0; metal_type < nmetyp; metal_type++)
  //   {
  //     /* If have at least one instance of this metal, write out
  //        the metal record */
  //     if (metal_count[metal_type] > 0)
  //       {
  //         /* Get the metal colour */
  //         colour = metal_type % NMETAL_COLOURS;

  //         /* Write out this metal type */
  //         fprintf(file_ptr,"M%s  %s\n",stored_metal[metal_type],
  //                 METAL_COLOUR[colour]);
  //       }
  //   }

  /* Get pointer to the first ligand record */
  ligand_ptr = first_ligand_ptr;

  /* Loop over the ligands */
  while (ligand_ptr != NULL)
    {
      /* Write out the details of ligand */
      if (ligand_ptr->metal == FALSE && ligand_ptr->natoms > 0)
        {
          /* Increment ligand counter */
          nligands++;

          /* Write out the residue list */
          fprintf(file_ptr,"L%3.3d %s\n",nligands,
                  ligand_ptr->residue_list);

          /* Write out the residue range */
          fprintf(file_ptr,"l%3.3d %s\n",nligands,
                  ligand_ptr->residue_range);
        }

      /* Get the next ligand */
      ligand_ptr = ligand_ptr->next_ligand_ptr;
    }

  /* Get pointer to the first ligand record again */
  ligand_ptr = first_ligand_ptr;
  last_sequence[0] = '\0';
  metal_type = 0;

  /* Loop over the ligand records to write out the metals */
  while (ligand_ptr != NULL)
    {
      /* Write out the details of the metal */
      if (ligand_ptr->metal == TRUE &&
          strcmp(ligand_ptr->residue_sequence,last_sequence))
        {
          /* Increment metal counter */
          metal_type++;

          /* Get the metal colour */
          colour = ligand_ptr->metal_colour;

          /* Write out the residue list */
          fprintf(file_ptr,"M%s %s\n",ligand_ptr->metal_name,
                  METAL_COLOUR[colour]);

          /* Store metal name */
          strcpy(last_sequence,ligand_ptr->residue_sequence);
        }

      /* Get the next ligand */
      ligand_ptr = ligand_ptr->next_ligand_ptr;
    }

  /* @@@ Write out the active site definitions */

  /* Close the rasmol.dfn file */
  fclose(file_ptr);
}
/***********************************************************************

reformat_hbplus_line  -  Write out the given HBPLUS line in .hb2 or
                         .nb2 format

***********************************************************************/

void reformat_hbplus_line(FILE *file_ptr,char *input_line)
{
  char output_line[LINELEN], string[LINELEN];

  /* Initialise output line */
  output_line[0] = '\0';

  /* Splice the relevant pieces of the old line into the new one */
  strncpy(string,input_line+10,14);
  string[14] = '\0';
  strcat(output_line,string);
  strncpy(string,input_line+30,14);
  string[14] = '\0';
  strcat(output_line,string);
  strncpy(string,input_line+45,4);
  string[4] = '\0';
  strcat(output_line,string);
  strcat(output_line," ");
  strncpy(string,input_line+49,6);
  string[6] = '\0';
  strcat(output_line,string);
  strcat(output_line," ");
  strncpy(string,input_line+55,5);
  string[5] = '\0';
  strcat(output_line,string);
  strcat(output_line," ");
  strncpy(string,input_line+60,5);
  string[5] = '\0';
  strcat(output_line,string);
  strcat(output_line," ");
  strncpy(string,input_line+65,5);
  string[5] = '\0';
  strcat(output_line,string);
  strcat(output_line," ");
  strncpy(string,input_line+70,5);
  string[5] = '\0';
  strcat(output_line,string);
  strncpy(string,input_line+75,12);
  string[12] = '\0';
  strcat(output_line,string);

  /* Write the new line out */
  fprintf(file_ptr,"%s\n",output_line);
}
/***********************************************************************

residue_sort  -  Shell-sort which sorts the residue numbers and names
                 into ascending order

***********************************************************************/

void residue_sort(struct residue **residue_index_ptr,int nindex)
{
  char residue_name1[10], residue_name2[10];

  int endloop, i, indi, indl, j, k, l, lognb2, m, n, nn;
  float aln2i = (float) 1.4426950, tiny = (float) 1.e-5;
  float calc;
  double dble, dble1;

  struct residue *residue1_ptr, *residue2_ptr, *swap_residue_ptr;

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

              /* Get the two residue records */
              residue1_ptr = residue_index_ptr[indi];
              residue2_ptr = residue_index_ptr[indl];

              /* Form full residue names on which sort is based */
              residue_name1[0] = residue1_ptr->chain;
              residue_name1[1] = '\0';
              strcat(residue_name1,residue1_ptr->res_num);
              strcat(residue_name1,residue1_ptr->res_name);
              residue_name2[0] = residue2_ptr->chain;
              residue_name2[1] = '\0';
              strcat(residue_name2,residue2_ptr->res_num);
              strcat(residue_name2,residue2_ptr->res_name);

              /* If in wrong order, the swap sort-pointers */
              if (strcmp(residue_name1,residue_name2) > 0)
                {
                  swap_residue_ptr = residue_index_ptr[indi];
                  residue_index_ptr[indi] = residue_index_ptr[indl];
                  residue_index_ptr[indl] = swap_residue_ptr;

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

binary_residue_search  -  Search for the given residue in the sorted
                          list using a binary search

***********************************************************************/

struct residue *binary_residue_search(struct residue **residue_index_ptr,
                                      int nresid,char *res_name,
                                      char *res_num,char chain)
{
  char residue_name[10], match_residue_name[10];

  int bottom_point, icomp, ipoint, new_point, top_point;
  int done, found;

  struct residue *match_residue_ptr;

  /* Form full residue name to be found */
  residue_name[0] = chain;
  residue_name[1] = '\0';
  strcat(residue_name,res_num);
  strcat(residue_name,res_name);

  /* Make starting point half-way into the index */
  ipoint = nresid / 2;
  bottom_point = 0;
  top_point = nresid - 1;

  /* Initialise variables */
  done = FALSE;
  found = FALSE;

  /* Loop until required array position found */
  while (done == FALSE)
    {
      /* Get the associated pointer */
      match_residue_ptr = residue_index_ptr[ipoint];

      /* Form the full name of this residue */
      match_residue_name[0] = match_residue_ptr->chain;
      match_residue_name[1] = '\0';
      strcat(match_residue_name,match_residue_ptr->res_num);
      strcat(match_residue_name,match_residue_ptr->res_name);

      /* Compare the two residue names */
      icomp = strcmp(residue_name,match_residue_name);

      /* Check if array value is equal to the value we want */
      if (icomp == 0)
        {
          /* Set flag to indicate that match found */
          done = TRUE;
          found = TRUE;
        }

      /* If array value is higher than the one we want, then check
         whether value of entry before is lower */
      else if (icomp < 0)
        {
          /* If this is the first entry in the sorted list, then
             can't go any lower */
          if (ipoint == 0)
            {
              /* Set flag that we are done */
              done = TRUE;
            }

          /* Otherwise, check whether the previous entry's name
             is below the one we're searching for */
          else
            {
              /* Get the pointer to the previous entry */
              match_residue_ptr = residue_index_ptr[ipoint - 1];

              /* Form the full name of this residue */
              match_residue_name[0] = match_residue_ptr->chain;
              match_residue_name[1] = '\0';
              strcat(match_residue_name,match_residue_ptr->res_num);
              strcat(match_residue_name,match_residue_ptr->res_name);

              /* If value is below the search value, we've found the
                 position in the index we're after, but don't have an exact
                 match*/
              if (strcmp(residue_name,match_residue_name) > 0)
                done = TRUE;

              /* Otherwise, we need to search earlier in the list */
              else
                {
                  /* Save the current position as the uppermost position
                     to search on */
                  top_point = ipoint;

                  /* Set pointer to halfway between the last search
                     low-search and here */
                  ipoint = (ipoint + bottom_point) / 2;
                }
            }

        }

      /* If current entry's value is below the search value, have to
         look further up the list */
      else
        {
          /* Save the current position as the lowest position to search
             on */
          bottom_point = ipoint;

          /* If have reached the top of the list, then our value is
             higher than any value in the list */
          if (ipoint >= nresid - 1)
            {
              /* Set flag as done */
              done = TRUE;
            }

          /* Otherwise, can search further up the list */
          else
            {
              /* Set pointer to halfway between the last search high-search
                 and here */
              new_point = (ipoint + top_point) / 2;

              /* Adjust if going back to the same point */
              if (new_point == ipoint)
                new_point++;
              ipoint = new_point;
            }
        }
    }

  /* Return the required entry in the index array */
  if (found == FALSE)
    match_residue_ptr = NULL;
  return(match_residue_ptr);
}
/***********************************************************************

fix_res_number  -  Change the HBPLUS-format residue number into the
                   correct representation

***********************************************************************/

void fix_res_number(char *res_num)
{
  int ipos;
  int minus, nonzero;

  /* Initialise variables */
  minus = FALSE;
  nonzero = FALSE;

  /* Loop to replace any leading zeros by spaces */
  for (ipos = 0; ipos < 4 && nonzero == FALSE; ipos ++)
    {
      if (res_num[ipos] == '0' && ipos < 3)
        res_num[ipos] = ' ';
      else if (res_num[ipos] == '-')
        {
          minus = TRUE;
          res_num[ipos] = ' ';
        }
      else
        {
          nonzero = TRUE;
          if (minus == TRUE && ipos > 0)
            res_num[ipos - 1] = '-';
        }
    }

  /* Replace hyphen in insertion code with a space */
  if (res_num[4] == '-')
    res_num[4] = ' ';
}
/***********************************************************************

find_atom  -  Find given atom in given residue

***********************************************************************/

struct atom *find_atom(struct residue *residue_ptr,char *atom_name)
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

write_grow_line  -  Write out the given interaction to the grow.out
                    file

***********************************************************************/

void write_grow_line(FILE *file_ptr,char intype,int *type,
                     struct residue *residue1_ptr,struct atom *atom1_ptr,
                     struct residue *residue2_ptr,struct atom *atom2_ptr,
                     double bond_length)
{
  char arrows[4], ctype1, ctype2, cswap;

  static char mol_type[6] = { 'P', 'N', 'L', 'P', 'M', 'X' };

  int swap;

  struct atom *swap_atom_ptr;
  struct residue *swap_residue_ptr;

  /* Get the molecule types */
  ctype1 = mol_type[type[0]];
  ctype2 = mol_type[type[1]];

  /* See if need to swap the two atoms round */
  swap = FALSE;
  if (ctype1 != 'P')
    swap = TRUE;

  /* If required, perform the swap */
  if (swap == TRUE)
    {
      /* Swap pointers */
      swap_atom_ptr = atom1_ptr;
      atom1_ptr = atom2_ptr;
      atom2_ptr = swap_atom_ptr;
      swap_residue_ptr = residue1_ptr;
      residue1_ptr = residue2_ptr;
      residue2_ptr = swap_residue_ptr;
      cswap = ctype1;
      ctype1 = ctype2;
      ctype2 = cswap;
    }

  /* Check if one of the atoms is a metal */
  if (intype == 'N' && (ctype1 == 'M' || ctype2 == 'M'))
    {
      /* Interaction with metal */
      intype = 'M';
      strcpy(arrows,":-:");
    }

  /* Non-bonded interaction */
  else if (intype == 'N')
    strcpy(arrows,"---");

  /* Disulphide bond */
  else if (intype == 'D')
    strcpy(arrows,"S-S");

  /* Salt bridge */
  else if (intype == 'S')
    {
      /* If first residue is negative, take other to be positive */
      if (!strcmp(residue1_ptr->res_name,"ASP") ||
          !strcmp(residue1_ptr->res_name,"GLU"))
        strcpy(arrows,"-.+");

      /* Otherwise, take the charges to be the other way round */
      else
        strcpy(arrows,"+.-");
    }

  /* For H-bond, determine its "direction" */
  else if (intype == 'H')
    {
      if (swap == FALSE)
        strcpy(arrows,"-->");
      else
        strcpy(arrows,"<--");
    }

  /* Write out start of line and first atom */
  fprintf(file_ptr,"%c %c %c %c  %s%s%6d %s ",intype,ctype1,ctype2,
          residue1_ptr->chain,residue1_ptr->res_num,
          residue1_ptr->res_name,atom1_ptr->atom_number,
          atom1_ptr->atom_name);

  /* Write out "direction" of interaction */
  fprintf(file_ptr,"%s",arrows);

  /* Write out the second atom and the bond distance */
  fprintf(file_ptr," %s%6d %s  %s%c%6.2f    -\n",atom2_ptr->atom_name,
          atom2_ptr->atom_number,residue2_ptr->res_name,
          residue2_ptr->res_num,residue2_ptr->chain,bond_length);
}
/***********************************************************************

find_interface  -  Find the given interface

***********************************************************************/

struct interface *find_interface(struct interface *first_interface_ptr,
                                 int type,char *chains)
{
  struct interface *interface_ptr, *found_interface_ptr;

  /* Initialise variables */
  found_interface_ptr = NULL;

  /* Get pointer to the first interface */
  interface_ptr = first_interface_ptr;

  /* Loop over all the interfaces to find the one just read in */
  while (interface_ptr != NULL && found_interface_ptr == NULL)
    {
      /* If have a match, store the pointer */
      if (type == interface_ptr->type &&
          chains[0] == interface_ptr->chains[0] && 
          chains[1] == interface_ptr->chains[1])
        found_interface_ptr = interface_ptr;

      /* Get pointer to the next interface */
      interface_ptr = interface_ptr->next_interface_ptr;
    }

  /* Return the matching interface pointer */
  return(found_interface_ptr);
}
/***********************************************************************

create_interface_record  -  Create and initialise a new interface record

***********************************************************************/

struct interface *create_interface_record(struct interface **fst_interface_ptr,
                                          struct interface **lst_interface_ptr,
                                          int type,char *chains)
{
  int i;

  struct interface *interface_ptr, *first_interface_ptr, *last_interface_ptr;

  /* Initialise interface pointers */
  interface_ptr = NULL;
  first_interface_ptr = *fst_interface_ptr;
  last_interface_ptr = *lst_interface_ptr;

  /* Allocate memory for structure to hold interface info */
  interface_ptr = (struct interface *) malloc(sizeof(struct interface));
  if (interface_ptr == NULL)
    {
      printf("*** Can't allocate memory for struct interface\n");
      exit (1);
    }

  /* If this is the very first interface, then store pointer */
  if (first_interface_ptr == NULL)
    first_interface_ptr = interface_ptr;

  /* Add link from previous interface to the current one */
  if (last_interface_ptr != NULL)
    last_interface_ptr->next_interface_ptr = interface_ptr;
  last_interface_ptr = interface_ptr;

  /* Store the interface details */
  interface_ptr->type = type;
  interface_ptr->chains[0] = chains[0];
  interface_ptr->chains[1] = chains[1];
  interface_ptr->chains[2] = '\0';
  for (i = 0; i < NTYPES; i++)
    interface_ptr->count[i] = 0;
  interface_ptr->nresid[0] = 0;
  interface_ptr->nresid[1] = 0;
  interface_ptr->nlinks = 0;
  interface_ptr->unique = TRUE;
  interface_ptr->ncopies = 0;
  interface_ptr->chain_area[0] = 0;
  interface_ptr->chain_area[1] = 0;
  interface_ptr->interface_area[0] = 0;
  interface_ptr->interface_area[1] = 0;
  interface_ptr->copyof_interface_ptr = NULL;
  interface_ptr->first_respair_ptr = NULL;
  interface_ptr->last_respair_ptr = NULL;

  /* Initialise pointer to next record */
  interface_ptr->next_interface_ptr = NULL;

  /* Return current pointers */
  *fst_interface_ptr = first_interface_ptr;
  *lst_interface_ptr = last_interface_ptr;

  /* Return new interface pointer */
  return(interface_ptr);
}
/***********************************************************************

find_respair  -  Find the given link

***********************************************************************/

struct respair *find_respair(struct respair *first_respair_ptr,
                             struct residue **residue_ptr,int *new)
{
  struct respair *respair_ptr, *found_respair_ptr;

  /* Initialise variables */
  found_respair_ptr = NULL;
  new[0] = TRUE;
  new[1] = TRUE;

  /* Get pointer to the first link */
  respair_ptr = first_respair_ptr;

  /* Loop over all the links to find the one just read in */
  while (respair_ptr != NULL)
    {
      /* If have a match, store the pointer */
      if (residue_ptr[0] == respair_ptr->residue_ptr[0] && 
          residue_ptr[1] == respair_ptr->residue_ptr[1])
        found_respair_ptr = respair_ptr;

      /* Check if residues are already in interface records */
      if (residue_ptr[0] == respair_ptr->residue_ptr[0])
        new[0] = FALSE;
      if (residue_ptr[1] == respair_ptr->residue_ptr[1])
        new[1] = FALSE;

      /* Get pointer to the next link */
      respair_ptr = respair_ptr->next_respair_ptr;
    }

  /* Return the matching link pointer */
  return(found_respair_ptr);
}
/***********************************************************************

create_respair_record  -  Create and initialise a new link record

***********************************************************************/

struct respair *create_respair_record(struct respair **fst_respair_ptr,
                                      struct respair **lst_respair_ptr,
                                      struct residue **residue_ptr)
{
  int i, j;

  struct respair *respair_ptr, *first_respair_ptr, *last_respair_ptr;

  /* Initialise link pointers */
  respair_ptr = NULL;
  first_respair_ptr = *fst_respair_ptr;
  last_respair_ptr = *lst_respair_ptr;

  /* Allocate memory for structure to hold link info */
  respair_ptr = (struct respair *) malloc(sizeof(struct respair));
  if (respair_ptr == NULL)
    {
      printf("*** Can't allocate memory for struct respair\n");
      exit (1);
    }

  /* If this is the very first link, then store pointer */
  if (first_respair_ptr == NULL)
    first_respair_ptr = respair_ptr;

  /* Add link from previous link to the current one */
  if (last_respair_ptr != NULL)
    last_respair_ptr->next_respair_ptr = respair_ptr;
  last_respair_ptr = respair_ptr;

  /* Store the link details */
  respair_ptr->residue_ptr[0] = residue_ptr[0];
  respair_ptr->residue_ptr[1] = residue_ptr[1];
  for (i = 0; i < NTYPES; i++)
    respair_ptr->count[i] = 0;

  /* Initialise coordinates */
  for (i = 0; i < 2; i++)
    {
      for (j = 0; j < 3; j++)
        respair_ptr->ave_coords[i][j] = 0.0;
      respair_ptr->natoms[i] = 0;
    }

  /* Initialise pointer to next record */
  respair_ptr->next_respair_ptr = NULL;

  /* Return current pointers */
  *fst_respair_ptr = first_respair_ptr;
  *lst_respair_ptr = last_respair_ptr;

  /* Return new link pointer */
  return(respair_ptr);
}
/***********************************************************************

process_interface  -  Store the interface information just read in from
                      the HBPLUS file

***********************************************************************/

void process_interface(struct interface **fst_interface_ptr,
                       struct interface **lst_interface_ptr,int *niface,
                       struct atom *atom1_ptr,struct atom *atom2_ptr,
                       int interface_type,int type)
{
  char chain, chains[2];

  int i, j, new[2], ninterfaces;

  struct atom *atom_ptr[2], *swap_atom_ptr;
  struct interface *interface_ptr, *first_interface_ptr, *last_interface_ptr;
  struct residue *residue_ptr[2], *swap_residue_ptr;
  struct respair *respair_ptr;

  /* Initialise interface data */
  first_interface_ptr = *fst_interface_ptr;
  last_interface_ptr = *lst_interface_ptr;
  ninterfaces = *niface;

  /* Get the two atoms and residues */
  atom_ptr[0] = atom1_ptr;
  atom_ptr[1] = atom2_ptr;
  residue_ptr[0] = atom1_ptr->residue_ptr;
  residue_ptr[1] = atom2_ptr->residue_ptr;

  /* Get the two chains */
  chains[0] = residue_ptr[0]->chain;
  chains[1] = residue_ptr[1]->chain;

  /* If chains are different, then check whether we
     already have this interface */
  if (chains[0] != chains[1])
    {
      /* Swap the chains and residues if necessary */
      if (chains[0] > chains[1])
        {
          /* Perform the swap on the chains */
          chain = chains[0];
          chains[0] = chains[1];
          chains[1] = chain;

          /* Swap the atoms */
          swap_atom_ptr = atom_ptr[0];
          atom_ptr[0] = atom_ptr[1];
          atom_ptr[1] = swap_atom_ptr;

          /* Swap the residues */
          swap_residue_ptr = residue_ptr[0];
          residue_ptr[0] = residue_ptr[1];
          residue_ptr[1] = swap_residue_ptr;
        }

      /* Check to see if we already have this interface */
      interface_ptr
        = find_interface(first_interface_ptr,interface_type,chains);

      /* If not found, create a new interface record */
      if (interface_ptr == NULL)
        {
          /* Create interface record */
          interface_ptr
            = create_interface_record(&first_interface_ptr,
                                      &last_interface_ptr,
                                      interface_type,chains);

          /* Increment interface count */
          ninterfaces++;
        }

      /* Check for the link between these two residues */
      respair_ptr = find_respair(interface_ptr->first_respair_ptr,
                                 residue_ptr,new);

      /* If either residue new, increment count of residues on each side
         of the interface */
      if (new[0] == TRUE)
        interface_ptr->nresid[0]++;
      if (new[1] == TRUE)
        interface_ptr->nresid[1]++;

      /* If not found, create a new respair record */
      if (respair_ptr == NULL)
        {
          /* Create a respair record */
          respair_ptr
            = create_respair_record(&(interface_ptr->first_respair_ptr),
                                    &(interface_ptr->last_respair_ptr),
                                    residue_ptr);

          /* Update count of links for this interface */
          interface_ptr->nlinks++;
        }

      /* Updated the current respair record with the current interaction
         type */
      respair_ptr->count[type]++;

      /* Repeat for overall count of interactions across this interface */
      interface_ptr->count[type]++;

      /* Accumulate coordinates for obtaining each residue's average
         interaction position (used by dimhtml for plotting) */
      for (i = 0; i < 2; i++)
        {
          /* Add coordinates */
          for (j = 0; j < 3; j++)
            respair_ptr->ave_coords[i][j]
              = respair_ptr->ave_coords[i][j] + atom_ptr[i]->coords[j];

          /* Update atom count */
          respair_ptr->natoms[i]++;
        }
    }

  /* Return interface data */
  *fst_interface_ptr = first_interface_ptr;
  *lst_interface_ptr = last_interface_ptr;
  *niface = ninterfaces;
}
/***********************************************************************

process_hbplus_file  -  Read in the H-bond and non-bonded interactions
                        from the .hhb or .nnb file and write out all
                        protein-ligand, protein-DNA and protein-metal
                        interactions to the grow.out file, plus all the
                        interactions across protein-protein interfaces

***********************************************************************/

int process_hbplus_file(FILE *outfile_ptr,FILE *outfile_i_ptr,
                        char *pdb_code,char *pdb_file,
                        struct residue **residue_index_ptr,int nresid,
                        struct interface **fst_interface_ptr,
                        struct interface **lst_interface_ptr,int *niface,
                        char *extension,int cando_nucplot,
                        char *output_dir,int *ng,int *nig,
                        struct parameters params)
{
  char file_name[FILENAME_LEN + 1], input_line[LINELEN], intype;
  char extension2[5], upp_code[5];
  char chain, res_name[4], res_num[6];
  char atom_name[5], number_string[20];

  int i, line, interaction_type, interface_type, ninterfaces, type[2];
  int ngrow, nigrow;
  int have_file, have_protein, have_ligand, have_dna, have_metal;
  int have_calpha, have_other;

  double bond_length;

  struct atom *atom1_ptr, *atom2_ptr;
  struct interface *first_interface_ptr, *last_interface_ptr;
  struct residue *residue1_ptr, *residue2_ptr;

  FILE *file_ptr, *newfile_ptr;

  /* Initialise variables */
  have_file = FALSE;
  line = 0;
  ninterfaces = *niface;
  ngrow = *ng;
  nigrow = *nig;

  /* Initialise interface pointers */
  first_interface_ptr = *fst_interface_ptr;
  last_interface_ptr = *lst_interface_ptr;

  /* Determine types of interactions we are dealing with */
  if (!strcmp(extension,".hhb"))
    {
      intype = 'H';
      interaction_type = HBONDS;
      strcpy(extension2,".hb2");
    }
  else
    {
      intype = 'N';
      interaction_type = CONTACTS;
      strcpy(extension2,".nb2");
    }

  /* Form name of input HBPLUS file */
  if (params.pdb_type == PDBSUM1)
    sprintf(file_name,"../hbplus/%s%s",pdb_code,extension);
  else
    {
      strcpy(file_name,output_dir);
      strcat(file_name,"file");
      strcat(file_name,extension);
    }

  /* If have an uploaded PDB file, input file will be the seq001.pdb
     file created by this program */
  if (params.uploaded_file == TRUE)
    {
      strcpy(file_name,output_dir);
      strcat(file_name,"seq001");
      strcat(file_name,extension);
    }

  /* Open this file */
  printf("Opening HBPLUS file (%s) ...\n",file_name);
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    {
      printf("* Warning. Unable to locate input HBPLUS file ... %s\n",
             file_name);
      return(have_file);
    }

  /* Open the output .hb2 or .nb2 file, being the new-format version
     of the .hhb and .nnb files. (Written out here to save rerunning
     HBPLUS for NUCPLOT, if necessary) */
  if (cando_nucplot == TRUE)
    {
      /* Form name of the output .hb2 or .nb2 file */
      strcpy(file_name,output_dir);
      strcat(file_name,"file");
      strcat(file_name,extension2);

      /* Open this file */
      newfile_ptr = open_output_file(file_name,TRUE);

      /* Convert PDB code to upper case */
      strcpy(upp_code,pdb_code);
      to_upper(upp_code,4);

      /* Write out the header records */
      fprintf(newfile_ptr,
              "HBPLUS Hydrogen Bond Calculator v 3.15            Ja"
              "n 29 16:08:01 GMT 2004\n"
              "(c) I McDonald, D Naylor, D Jones and J Thornton 199"
              "3 All Rights Reserved.\n"
              "Citing HBPLUS in publications that use these results"
              "is condition of use.\n"
              "%s <- Brookhaven Code \"%s\" <- PDB file\n"
              "<---DONOR---> <-ACCEPTOR-->    atom                 "
              "^               \n"
              "c    i                          cat <-CA-CA->   ^   "
              "H-A-AA   ^      H- \n"
              "h    n   atom  resd res      DA  || num        DHA  "
              "H-A  angle D-A-AA Bond\n"
              "n    s   type  num  typ     dist DA aas  dist angle "
              "dist       angle   num\n",upp_code,pdb_file);
    }

  /* Loop while reading in records from the file */
  while (fgets(input_line,LINELEN,file_ptr))
    {
      /* Increment line counter */
      line++;

      /* If writing out to the .hb2 or .nb2 files, reformat this
         line and write out */
      if (cando_nucplot == TRUE)
        reformat_hbplus_line(newfile_ptr,input_line);

      /* Get the residue name of the first residue on the line */
      strncpy(res_name,input_line+16,3);
      res_name[3] = '\0';
      strncpy(res_num,input_line+11,5);
      res_num[5] = '\0';
      chain = input_line[10];
      if (chain == '-')
        chain = ' ';

      /* Fix the residue number's format */
      fix_res_number(res_num);

      /* Find this residue if not water */
      if (!strcmp(res_name,"HOH"))
        residue1_ptr = NULL;
      else
        residue1_ptr
          = binary_residue_search(residue_index_ptr,nresid,res_name,
                                  res_num,chain);

      /* If first residue found, look for second residue */
      if (residue1_ptr != NULL)
        {
          /* Get the residue name of the second residue on the line */
          strncpy(res_name,input_line+36,3);
          res_name[3] = '\0';
          strncpy(res_num,input_line+31,5);
          res_num[5] = '\0';
          chain = input_line[30];
          if (chain == '-')
            chain = ' ';

          /* Fix the residue number's format */
          fix_res_number(res_num);

          /* Find this residue if not water */
          if (!strcmp(res_name,"HOH"))
            residue2_ptr = NULL;
          else
            residue2_ptr
              = binary_residue_search(residue_index_ptr,nresid,res_name,
                                      res_num,chain);

          /* If residue found, get the molecule type that each
             belongs to */
          if (residue2_ptr != NULL)
            {
              /* Initialise flags */
              have_calpha = FALSE;
              have_protein = FALSE;
              have_ligand = FALSE;
              have_dna = FALSE;
              have_metal = FALSE;
              have_other = FALSE;

              /* Get the two chain types */
              type[0] = residue1_ptr->chain_ptr->type;
              type[1] = residue2_ptr->chain_ptr->type;

              /* Determine molecule to which each residue belongs */
              for (i = 0; i < 2; i++)
                {
                  if (type[i] == PROTEIN || type[i] == CALPHA_ONLY)
                    have_protein = TRUE;
                  else if (type[i] == DNA)
                    have_dna = TRUE;
                  else if (type[i] == LIGAND)
                    have_ligand = TRUE;
                  else if (type[i] == METAL)
                    have_metal = TRUE;
                  else
                    have_other = TRUE;
                  if (type[i] == CALPHA_ONLY)
                    have_calpha = TRUE;
                }

              /* Get the first atom name */
              strncpy(atom_name,input_line+19,4);
              atom_name[4] = '\0';

              /* Identify the first atom */
              atom1_ptr = find_atom(residue1_ptr,atom_name);

              /* Get the second atom name */
              strncpy(atom_name,input_line+39,4);
              atom_name[4] = '\0';

              /* Identify the second atom */
              atom2_ptr = find_atom(residue2_ptr,atom_name);

              /* Bond length */
              strncpy(number_string,input_line+45,4);
              number_string[4] = '\0';
              bond_length = atof(number_string);

              /* If have a protein and ligand/metal/dna then write
                 this interaction out */
              if ((have_protein == TRUE && have_dna == TRUE) ||
                  ((have_protein == TRUE || have_dna == TRUE) &&
                   (have_ligand == TRUE || have_metal == TRUE)))
                {
                  /* Write out the line to the grow.out file */
                  if (atom1_ptr != NULL && atom2_ptr != NULL)
                    {
                      write_grow_line(outfile_ptr,intype,type,
                                      residue1_ptr,atom1_ptr,
                                      residue2_ptr,atom2_ptr,bond_length);
                      ngrow++;
                    }
                }

              /* If have a protein-protein or protein-DNA interface,
                 save the chains involved */
              /* Currently only interested in protein-protein
                 interfaces */
              if (residue1_ptr->chain != residue2_ptr->chain &&
                  // ((have_protein == TRUE && have_dna == TRUE) ||
                  (have_protein == TRUE && have_dna == FALSE &&
                   have_ligand == FALSE && have_metal == FALSE &&
                   have_other == FALSE && have_calpha == FALSE))
                {
                  /* Set interface type */
                  if (have_protein == TRUE && have_dna == TRUE)
                    interface_type = PROTEIN_DNA;
                  else
                    interface_type = PROTEIN_PROTEIN;

                  /* Have interface, so store details */
                  if (atom1_ptr != NULL && atom2_ptr != NULL)
                    process_interface(&first_interface_ptr,
                                      &last_interface_ptr,&ninterfaces,
                                      atom1_ptr,atom2_ptr,
                                      interface_type,interaction_type);

                  /* Write out the interaction to the igrow.out file */
                  if (atom1_ptr != NULL && atom2_ptr != NULL)
                    {
                      write_grow_line(outfile_i_ptr,intype,type,
                                      residue1_ptr,atom1_ptr,
                                      residue2_ptr,atom2_ptr,bond_length);
                      nigrow++;
                    }
                }
            }
        }
    }

  /* Close the file */
  fclose(file_ptr);

  /* Set flag saying whether we have file */
  if (line > 0)
    have_file = TRUE;

  /* Close the output .hb2 or .nb2 file, if we have one */
  if (cando_nucplot == TRUE)
    fclose(newfile_ptr);

  /* Return interface data */
  *fst_interface_ptr = first_interface_ptr;
  *lst_interface_ptr = last_interface_ptr;
  *niface = ninterfaces;

  /* Return counts of lines written out */
  *ng = ngrow;
  *nig = nigrow;

  /* Return whether have the HBPLUS file */
  return(have_file);
}
/***********************************************************************

run_hbplus  -  Run HBPLUS on the seq001.pdb file to compute all the
               H-bonds and non-bonded contacts between the first chain
               and any ligands or DNA

***********************************************************************/

void run_hbplus(char *hbplus_exe,char *output_dir)
{
  char command[FILENAME_LEN + 1], filename[FILENAME_LEN + 1];

  /* Form name of seq001.pdb file */
  sprintf(filename,"%sseq001.pdb",output_dir);
 
  /* Form the parameters to run the HBPLUS program to compute all
     non-bonded contacts */
  sprintf(command,"%s -L -h 2.9 -d 3.9 -N %s > nnb.log",hbplus_exe,
          filename);

  /* Run the HBPLUS program */
  system(command);

  /* Run HBPLUS a second time to compute all H-bonds */
  sprintf(command,"%s -L -h 2.7 -d 3.35 %s > nnb.log",hbplus_exe,
          filename);

  /* Run the HBPLUS program */
  system(command);
}
/***********************************************************************

store_iface_ssbonds  -  Identify any disulphide bonds across protein-
                        protein interfaces and store

***********************************************************************/

void store_iface_ssbonds(FILE *outfile_ptr,
                         struct interface **fst_interface_ptr,
                         struct interface **lst_interface_ptr,
                         int *niface,struct ssbond *first_ssbond_ptr)
{
  char intype;

  int ninterfaces;
  int interaction_type, interface_type, type[2];

  double bond_length, dist_sqrd, x1, x2, y1, y2, z1, z2;

  struct atom *atom1_ptr, *atom2_ptr;
  struct interface *first_interface_ptr, *last_interface_ptr;
  struct interface *interface_ptr;
  struct residue *residue1_ptr, *residue2_ptr;
  struct ssbond *ssbond_ptr;

  /* Initialise variables */
  ninterfaces = *niface;
  interface_type = PROTEIN_PROTEIN;
  interaction_type = DISULPHIDES;
  intype = 'D';

  /* Initialise interface pointers */
  interface_ptr = NULL;
  first_interface_ptr = *fst_interface_ptr;
  last_interface_ptr = *lst_interface_ptr;

  /* Get the first disulphide bond */
  ssbond_ptr = first_ssbond_ptr;

  /* Loop through the disulphide bonds */
  while (ssbond_ptr != NULL)
    {
      /* Get the two residues */
      residue1_ptr = ssbond_ptr->residue1_ptr;
      residue2_ptr = ssbond_ptr->residue2_ptr;

      /* Get the two atoms */
      atom1_ptr = ssbond_ptr->atom1_ptr;
      atom2_ptr = ssbond_ptr->atom2_ptr;

      /* Get the two chain types */
      type[0] = residue1_ptr->chain_ptr->type;
      type[1] = residue2_ptr->chain_ptr->type;

      /* If residues belong to separate chains, and both chains are
         proteins, store bond across the interface */
      if (residue1_ptr->chain != residue2_ptr->chain &&
          type[0] == PROTEIN && type[1] == PROTEIN && atom1_ptr != NULL &&
          atom2_ptr != NULL)
        {
          /* Have interface, so store details */
          process_interface(&first_interface_ptr,
                            &last_interface_ptr,&ninterfaces,
                            atom1_ptr,atom2_ptr,
                            interface_type,interaction_type);


          /* Get the coords of the two atoms */
          x1 = atom1_ptr->coords[0];
          y1 = atom1_ptr->coords[1];
          z1 = atom1_ptr->coords[2];
          x2 = atom2_ptr->coords[0];
          y2 = atom2_ptr->coords[1];
          z2 = atom2_ptr->coords[2];

          /* Calculate the length of this disulphide bond */
          dist_sqrd
            = (x2 - x1) * (x2 - x1)
            + (y2 - y1) * (y2 - y1)
            + (z2 - z1) * (z2 - z1);
          bond_length = sqrt(dist_sqrd);

          /* Write out the interaction to the igrow.out file */
          write_grow_line(outfile_ptr,intype,type,
                          residue1_ptr,atom1_ptr,
                          residue2_ptr,atom2_ptr,bond_length);
        }

      /* Get the next ssbond */
      ssbond_ptr = ssbond_ptr->next_ssbond_ptr;
    }

  /* Return number of interfaces */
  *niface = ninterfaces;
  
  /* Return current pointers */
  *fst_interface_ptr = first_interface_ptr;
  *lst_interface_ptr = last_interface_ptr;
}
/***********************************************************************

check_atom_pairs  -  Check whether carboxyl oxygens and side-chain
                     nitrogens of the given two residues are within the
                     cut-off distance for a salt bridge

***********************************************************************/

int check_atom_pairs(struct residue *residue1_ptr,
                     struct residue *residue2_ptr,
                     struct atom **sbridge_atom1_ptr,
                     struct atom **sbridge_atom2_ptr)
{
  int iatom, jatom, niatoms, njatoms;
  int salt_bridge;

  double dist_sqrd, min_dist_sqrd, x1, x2, y1, y2, z1, z2;

  struct atom *atom1_ptr, *atom2_ptr;

  /* Initialise variables */
  salt_bridge = FALSE;

  /* Initialise atom pointers */
  *sbridge_atom1_ptr = NULL;
  *sbridge_atom2_ptr = NULL;
  min_dist_sqrd = SBRIDGE_DIST_SQRD;

  /* Get the first of the first residue's atoms */
  atom1_ptr = residue1_ptr->first_atom_ptr;
  iatom = 0;
  niatoms = residue1_ptr->natoms;
  
  /* Loop through all the residue's atoms */
  while (atom1_ptr != NULL && iatom < niatoms)
    {
      /* Check whether this is the right kind of atom */
      if (!strcmp(atom1_ptr->atom_name," OD1") ||
          !strcmp(atom1_ptr->atom_name," OD2") ||
          !strcmp(atom1_ptr->atom_name," OE1") ||
          !strcmp(atom1_ptr->atom_name," OE2") ||
          !strcmp(atom1_ptr->atom_name," NE ") ||
          !strcmp(atom1_ptr->atom_name," NH1") ||
          !strcmp(atom1_ptr->atom_name," NH2") ||
          !strcmp(atom1_ptr->atom_name," NZ ") ||
          !strcmp(atom1_ptr->atom_name," ND1") ||
          !strcmp(atom1_ptr->atom_name," NE2"))
        {
          /* Store this atom's coords */
          x1 = atom1_ptr->coords[0];
          y1 = atom1_ptr->coords[1];
          z1 = atom1_ptr->coords[2];

          /* Get the first of the other residue's atoms */
          atom2_ptr = residue2_ptr->first_atom_ptr;
          jatom = 0;
          njatoms = residue2_ptr->natoms;
  
          /* Loop through all the residue's atoms */
          while (atom2_ptr != NULL && jatom < njatoms)
            {
              /* Check whether this is the right kind of atom */
              if (!strcmp(atom2_ptr->atom_name," OD1") ||
                  !strcmp(atom2_ptr->atom_name," OD2") ||
                  !strcmp(atom2_ptr->atom_name," OE1") ||
                  !strcmp(atom2_ptr->atom_name," OE2") ||
                  !strcmp(atom2_ptr->atom_name," NE ") ||
                  !strcmp(atom2_ptr->atom_name," NH1") ||
                  !strcmp(atom2_ptr->atom_name," NH2") ||
                  !strcmp(atom2_ptr->atom_name," NZ ") ||
                  !strcmp(atom2_ptr->atom_name," ND1") ||
                  !strcmp(atom2_ptr->atom_name," NE2"))
                {
                  /* Store this atom's coords */
                  x2 = atom2_ptr->coords[0];
                  y2 = atom2_ptr->coords[1];
                  z2 = atom2_ptr->coords[2];

                  /* Calculate the squared distance between these
                     two atoms */
                  dist_sqrd
                    = (x2 - x1) * (x2 - x1)
                    + (y2 - y1) * (y2 - y1)
                    + (z2 - z1) * (z2 - z1);

                  /* If within salt-bridge distance, we have a
                     salt-bridge */
                  if (dist_sqrd < SBRIDGE_DIST_SQRD)
                    {
                      /* Store the atoms */
                      if (dist_sqrd < min_dist_sqrd)
                        {
                          *sbridge_atom1_ptr = atom1_ptr;
                          *sbridge_atom2_ptr = atom2_ptr;
                          min_dist_sqrd = dist_sqrd;
                        }

                      /* Set flag */
                      salt_bridge = TRUE;
                    }
                }

              /* Get pointer to next atom in linked-list */
              atom2_ptr = atom2_ptr->next_atom_ptr;
              jatom++;
            }
        }
          
      /* Get pointer to next atom in linked-list */
      atom1_ptr = atom1_ptr->next_atom_ptr;
      iatom++;
    }

  /* Return whether atoms are sufficiently close for salt-bridge */
  return(salt_bridge);
}
/***********************************************************************

check_residue_pairs  -  Find any salt bridges across the protein-protein
                     interfaces

***********************************************************************/

void check_residue_pairs(FILE *outfile_ptr,struct chain *chain_ptr,
                         struct chain *other_chain_ptr,
                         struct interface **fst_interface_ptr,
                         struct interface **lst_interface_ptr,
                         int *niface)
{
  char intype;

  int ninterfaces;
  int interaction_type, interface_type, type[2];
  int i, iresid, jresid, niresid, njresid;
  int res_type1, res_type2;
  int toofar, wanted;

  double dist, dist_sqrd, x1, x2, y1, y2, z1, z2;

  struct atom *atom1_ptr, *atom2_ptr;
  struct interface *first_interface_ptr, *last_interface_ptr;
  struct interface *interface_ptr;
  struct residue *residue1_ptr, *residue2_ptr;

  /* Initialise variables */
  ninterfaces = *niface;
  interface_type = PROTEIN_PROTEIN;
  interaction_type = SALT_BRIDGES;
  intype = 'S';

  /* Initialise interface pointers */
  interface_ptr = NULL;
  first_interface_ptr = *fst_interface_ptr;
  last_interface_ptr = *lst_interface_ptr;

  /* Get the first of the first chain's residues */
  residue1_ptr = chain_ptr->first_residue_ptr;
  niresid = chain_ptr->nresid;
  iresid = 0;

  /* Loop through the residues */
  while (residue1_ptr != NULL && iresid < niresid)
    {
      /* Check if residue wanted */
      res_type1 = NOT_WANTED;
      if (!strcmp(residue1_ptr->res_name,"ASP") ||
          !strcmp(residue1_ptr->res_name,"GLU"))
        res_type1 = NEGATIVE;
      else if (!strcmp(residue1_ptr->res_name,"ARG") ||
               !strcmp(residue1_ptr->res_name,"LYS") ||
               !strcmp(residue1_ptr->res_name,"HIS"))
        res_type1 = POSITIVE;

      /* If have a charged amino acid, look for oppositely charged
         residues on the other chain */
      if (res_type1 != NOT_WANTED)
        {
          /* Store the coords of this residue's sidechain centroid */
          x1 = residue1_ptr->centroid[0];
          y1 = residue1_ptr->centroid[1];
          z1 = residue1_ptr->centroid[2];

          /* Get the first of the first chain's residues */
          residue2_ptr = other_chain_ptr->first_residue_ptr;
          njresid = other_chain_ptr->nresid;
          jresid = 0;

          /* Loop through the residues to check for possible salt
             bridge partners */
          while (residue2_ptr != NULL && jresid < njresid)
            {
              /* Check if residue wanted */
              res_type2 = NOT_WANTED;
              if (!strcmp(residue2_ptr->res_name,"ASP") ||
                  !strcmp(residue2_ptr->res_name,"GLU"))
                res_type2 = NEGATIVE;
              else if (!strcmp(residue2_ptr->res_name,"ARG") ||
                       !strcmp(residue2_ptr->res_name,"LYS") ||
                       !strcmp(residue2_ptr->res_name,"HIS"))
                res_type2 = POSITIVE;

              /* If have the right combination of residue types, check
                 their proximity */
              if ((res_type1 == POSITIVE && res_type2 == NEGATIVE) ||
                  (res_type2 == POSITIVE && res_type1 == NEGATIVE))
                {
                  /* Check whether these residue are within reach */
                  toofar = FALSE;
                  for (i = 0; i < 3 && toofar == FALSE; i++)
                    {
                      if (residue1_ptr->coord_min[i]
                          - residue2_ptr->coord_max[i] > SALT_BRIDGE_DIST)
                        toofar = TRUE;
                      if (residue2_ptr->coord_min[i]
                          - residue1_ptr->coord_max[i] > SALT_BRIDGE_DIST)
                        toofar = TRUE;
                    }
      
                  /* If residues within range check actual separation */
                  if (toofar == FALSE)
                    {
                      /* Get the coords of this residue's sidechain
                         centroid */
                      x2 = residue2_ptr->centroid[0];
                      y2 = residue2_ptr->centroid[1];
                      z2 = residue2_ptr->centroid[2];

                      /* Calculate distance between the centroids of these
                         two residue sidechains */
                      dist_sqrd
                        = (x2 - x1) * (x2 - x1)
                        + (y2 - y1) * (y2 - y1)
                        + (z2 - z1) * (z2 - z1);

                      /* If centroid is within the target distance, need
                         to check atoms as well */
                      // OUT if (dist_sqrd < SBRIDGE_DIST_SQRD)
                        {
                          /* Check the atoms of these two residues to see
                             if charged atoms also within the cut-off for
                             a salt-bridge */
                          wanted
                            = check_atom_pairs(residue1_ptr,residue2_ptr,
                                               &atom1_ptr,&atom2_ptr);

                          /* If have salt-bridge store details */
                          if (wanted == TRUE)
                            {
                              /* Store this interface */
                              process_interface(&first_interface_ptr,
                                                &last_interface_ptr,
                                                &ninterfaces,
                                                atom1_ptr,atom2_ptr,
                                                interface_type,
                                                interaction_type);

                              /* Get the coords of the two atoms */
                              x1 = atom1_ptr->coords[0];
                              y1 = atom1_ptr->coords[1];
                              z1 = atom1_ptr->coords[2];
                              x2 = atom2_ptr->coords[0];
                              y2 = atom2_ptr->coords[1];
                              z2 = atom2_ptr->coords[2];

                              /* Calculate distance between the two
                                 closest charged atoms */
                              dist_sqrd
                                = (x2 - x1) * (x2 - x1)
                                + (y2 - y1) * (y2 - y1)
                                + (z2 - z1) * (z2 - z1);
                              dist = sqrt(dist_sqrd);

                              /* Get the two chain types */
                              type[0] = residue1_ptr->chain_ptr->type;
                              type[1] = residue2_ptr->chain_ptr->type;

                              /* Write out the interaction to the
                                 igrow.out file */
                              write_grow_line(outfile_ptr,intype,type,
                                              residue1_ptr,atom1_ptr,
                                              residue2_ptr,atom2_ptr,
                                              dist);
                            }
                        }
                    }
                }

              /* Get the next residue */
              residue2_ptr = residue2_ptr->next_residue_ptr;
              jresid++;
            }
        }

      /* Get the next residue */
      residue1_ptr = residue1_ptr->next_residue_ptr;
      iresid++;
    }

  /* Return number of interfaces */
  *niface = ninterfaces;
  
  /* Return current pointers */
  *fst_interface_ptr = first_interface_ptr;
  *lst_interface_ptr = last_interface_ptr;
}
/***********************************************************************

get_salt_bridges  -  Find any salt bridges across the protein-protein
                     interfaces

                     Salt bridge definition taken from
                     Kumar & Nussinov (1999). J.Mol.Biol., 293,1241-1255

***********************************************************************/

void get_salt_bridges(FILE *outfile_ptr,
                      struct interface **fst_interface_ptr,
                      struct interface **lst_interface_ptr,
                      int *niface,struct chain *first_chain_ptr)
{
  int ninterfaces;

  struct chain *chain_ptr, *other_chain_ptr;
  struct interface *first_interface_ptr, *last_interface_ptr;
  struct interface *interface_ptr;

  /* Initialise variables */
  ninterfaces = *niface;

  /* Initialise interface pointers */
  interface_ptr = NULL;
  first_interface_ptr = *fst_interface_ptr;
  last_interface_ptr = *lst_interface_ptr;

  /* Get pointer to the first chain */
  chain_ptr = first_chain_ptr;

  /* Loop over all the chains */
  while (chain_ptr != NULL)
    {
      /* If this is a protein chain, and is one of the interface chains,
         then process */
      if (chain_ptr->type == PROTEIN)
        {
          /* Get next chain */
          other_chain_ptr = chain_ptr->next_chain_ptr;

          /* Loop over all the other chains */
          while (other_chain_ptr != NULL)
            {
              /* If this is a protein chain */
              if (other_chain_ptr->type == PROTEIN)
                {
                  /* Compare the residues of each chain to check which
                     pairs might be forming a salt bridge across the
                     protein-protein interface */
                  check_residue_pairs(outfile_ptr,chain_ptr,other_chain_ptr,
                                      &first_interface_ptr,
                                      &last_interface_ptr,
                                      &ninterfaces);
                }

              /* Get the next chain */
              other_chain_ptr = other_chain_ptr->next_chain_ptr;
            }
        }
      /* Get the next chain */
      chain_ptr = chain_ptr->next_chain_ptr;
    }

  /* Return number of interfaces */
  *niface = ninterfaces;
  
  /* Return current pointers */
  *fst_interface_ptr = first_interface_ptr;
  *lst_interface_ptr = last_interface_ptr;
}
/***********************************************************************

create_grow_out  -  Summarise the H-bond and non-bonded interactions in
                    the .hhb and .nnb files to create a grow.out file

***********************************************************************/

void create_grow_out(char *pdb_code,char *pdb_file,
                     struct residue *first_residue_ptr,int nresid,
                     struct interface **fst_interface_ptr,int *niface,
                     struct ssbond *first_ssbond_ptr,
                     struct chain *first_chain_ptr,int cando_nucplot,
                     char *output_dir,struct parameters params)
{
  char file_name[FILENAME_LEN + 1], ifile_name[FILENAME_LEN + 1];
  char command[FILENAME_LEN + 1];

  int ires, ninterfaces, ngrow, nigrow;
  int have_hhb, have_nnb;

  struct interface *first_interface_ptr, *last_interface_ptr;
  struct residue *residue_ptr;

  struct residue **residue_index_ptr;
  
  FILE *file_ptr, *ifile_ptr;

  /* Initialise interface pointers */
  first_interface_ptr = NULL;
  last_interface_ptr = NULL;
  ninterfaces = 0;
  ngrow = nigrow = 0;

  /* Create arrays of residue pointers and sort index */
  residue_index_ptr
    = (struct residue **) malloc(nresid * sizeof(struct residue *));
  if (residue_index_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for residue-sort\n");
      exit(1);
    }
  
  /* Get the first residue */
  residue_ptr = first_residue_ptr;
  ires = 0;

  /* Loop through the residues to place their pointers in the array */
  while (residue_ptr != NULL && ires < nresid)
    {
      /* Place residue pointer on the array */
      residue_index_ptr[ires] = residue_ptr;

      /* Get the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
      ires++;
    }

  /* Sort the residue pointers into order of residue number and name */
  if (nresid > 1)
    residue_sort(residue_index_ptr,nresid);

  /* Form name of output grow.out */
  strcpy(file_name,output_dir);
  strcat(file_name,"grow.out");

  /* Form name of output igrow.out */
  strcpy(ifile_name,output_dir);
  strcat(ifile_name,"igrow.out");

  /* Open the grow.out and igrow.out output files */
  file_ptr = open_output_file(file_name,TRUE);
  ifile_ptr = open_output_file(ifile_name,TRUE);

  /* Write out header line */
  fprintf(file_ptr,"# Grow v2.00  output from newsec\n");

  /* Process the .hhb file, writing all the H-bonds between protein
     and ligands, DNA or metals out to the grow.out file */
  have_hhb = process_hbplus_file(file_ptr,ifile_ptr,pdb_code,pdb_file,
                                 residue_index_ptr,nresid,
                                 &first_interface_ptr,
                                 &last_interface_ptr,&ninterfaces,".hhb",
                                 cando_nucplot,output_dir,&ngrow,&nigrow,
                                 params);

  /* Repeat for the non-bonded interactions in the .nnb file */
  have_nnb = process_hbplus_file(file_ptr,ifile_ptr,pdb_code,pdb_file,
                                 residue_index_ptr,nresid,
                                 &first_interface_ptr,
                                 &last_interface_ptr,&ninterfaces,".nnb",
                                 cando_nucplot,output_dir,&ngrow,&nigrow,
                                 params);

  /* Check for any disulphide bonds across protein-protein interfaces */
  store_iface_ssbonds(ifile_ptr,&first_interface_ptr,&last_interface_ptr,
                      &ninterfaces,first_ssbond_ptr);

  /* Check for any salt bridge across protein-protein interfaces */
  get_salt_bridges(ifile_ptr,&first_interface_ptr,&last_interface_ptr,
                   &ninterfaces,first_chain_ptr);

  /* Close the grow.out and igrow.out files */
  fclose(file_ptr);
  fclose(ifile_ptr);

  /* If neither file present, or both are empty, delete the grow.out
     files just created */
  if (have_hhb == FALSE && have_nnb == FALSE)
    {
      /* Create the system commands to delete the files */
      sprintf(command,"%s -f grow.out",RM_COMMAND);
      system(command);
      sprintf(command,"%s -f igrow.out",RM_COMMAND);
      system(command);
    }
  else if (ngrow == 0 || nigrow == 0)
    {
      if (ngrow == 0)
        {
          sprintf(command,"%s -f grow.out",RM_COMMAND);
          system(command);
        }
      if (nigrow == 0)
        {
          sprintf(command,"%s -f igrow.out",RM_COMMAND);
          system(command);
        }
    }

  /* Free up the memory grabbed at the start of the function */
  free(residue_index_ptr);

  /* Return interface data */
  *fst_interface_ptr = first_interface_ptr;
  *niface = ninterfaces;
}
/***********************************************************************

interface_sort  -  Shell-sort to sorts the interfaces alphabetically

***********************************************************************/

void interface_sort(int *index,struct interface **interface_index_ptr,
                    int nindex)
{
  int endloop, i, indi, indl, iswap, j, k, l, lognb2, m, n, nn;
  float aln2i = (float) 1.4426950, tiny = (float) 1.e-5;
  float calc;
  double dble, dble1;

  struct interface *interface1_ptr, *interface2_ptr;

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

              /* Get the two interface records */
              interface1_ptr = interface_index_ptr[index[indi]];
              interface2_ptr = interface_index_ptr[index[indl]];

              /* If in wrong order, the swap sort-pointers */
              if (strcmp(interface1_ptr->chains,interface2_ptr->chains) > 0)
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

sort_interfaces  -  Sort the interfaces alphabetically

***********************************************************************/

void sort_interfaces(struct interface **fst_interface_ptr,int ninterfaces)
{
  int iorder, ipos;

  int *index;

  struct interface *first_interface_ptr, *last_interface_ptr;
  struct interface *interface_ptr;
  struct interface **interface_index_ptr;

  /* Initialise pointers */
  first_interface_ptr = *fst_interface_ptr;
  last_interface_ptr = NULL;

  /* Grab memory for temporary sort arrays */
  index = (int *) malloc(ninterfaces * sizeof(int));
  interface_index_ptr
    = (struct interface **) malloc(ninterfaces * sizeof(struct interface *));

  /* Get pointer to first interface record */
  interface_ptr = first_interface_ptr;
  ipos = 0;

  /* Loop through interface records to initialise the sort index */
  while (interface_ptr != NULL)
    {
      /* Initialise sort index */
      index[ipos] = ipos;
      interface_index_ptr[ipos] = interface_ptr;

      /* Increment interface counter */
      ipos++;

      /* Get next interface record */
      interface_ptr = interface_ptr->next_interface_ptr;
    }

  /* Initialise interface counter */
  ninterfaces = ipos;

  /* Sort the interfaces */
  interface_sort(index,interface_index_ptr,ninterfaces);

  /* Having sorted the interfaces, rearrange order of the interface-record
     pointers */
  for (iorder = 0; iorder < ninterfaces; iorder++)
    {
      /* Get the array position of this interface record */
      ipos = index[iorder];

      /* Get the current interface record */
      interface_ptr = interface_index_ptr[ipos];

      /* If this is the first, then make the first record */
      if (iorder == 0)
        first_interface_ptr = interface_ptr;

      /* Otherwise, get the previous interface record to point to this one */
      else
        last_interface_ptr->next_interface_ptr = interface_ptr;

      /* Blank out current record's next interface pointer */
      interface_ptr->next_interface_ptr = NULL;

      /* Store current record as last encountered */
      last_interface_ptr = interface_ptr;
    }

  /* Free up memory used for temporary sort arrays */
  free(index);
  free(interface_index_ptr);

  /* Return current pointer */
  *fst_interface_ptr = first_interface_ptr;
}
/***********************************************************************

equiv_chain  -  If this is a copy chain, return name of original chain

***********************************************************************/

char equiv_chain(char chain,struct chain *first_chain_ptr)
{
  char orig_chain;

  int found;

  struct chain *chain_ptr;

  /* Initialise variables */
  orig_chain = chain;
  found = FALSE;

  /* Get pointer to the first chain */
  chain_ptr = first_chain_ptr;

  /* Loop over all the chains */
  while (chain_ptr != NULL && found == FALSE)
    {
      /* If this is the right chain, retrieve original */
      if (chain_ptr->chain_id == chain)
        {
          /* If this is a copy chain, retrieve the original */
          if (chain_ptr->copy_of != '\0')
            orig_chain = chain_ptr->copy_of;

          /* Set flag that we're done */
          found = TRUE;
        }

      /* Get the next chain */
      chain_ptr = chain_ptr->next_chain_ptr;
    }

  /* Return the original chain */
  return(orig_chain);
}
/***********************************************************************

count_pairs  -  Check number of common residue pairings the two given
                interfaces have to determine whether they can be taken
                as equivalent

***********************************************************************/

int count_pairs(struct interface *interface_ptr,
                struct interface *other_interface_ptr,int equivalence)
{
  int nequiv, nlinks;
  int equivalent, found;

  struct residue *residue_ptr[2], *other_residue_ptr[2];
  struct respair *respair_ptr, *other_respair_ptr;

  /* Initialise variables */
  equivalent = FALSE;
  nequiv = 0;
  nlinks = interface_ptr->nlinks;
  if (other_interface_ptr->nlinks > nlinks)
    nlinks = other_interface_ptr->nlinks;

  /* Get pointer to the first link */
  respair_ptr = interface_ptr->first_respair_ptr;

  /* Loop through all the interactions of the first interface */
  while (respair_ptr != NULL)
    {
      /* Get the two residues */
      residue_ptr[0] = respair_ptr->residue_ptr[0];
      residue_ptr[1] = respair_ptr->residue_ptr[1];

      /* Get pointer to the first link of second interface */
      other_respair_ptr = other_interface_ptr->first_respair_ptr;
      found = FALSE;

      /* Loop through all the interactions of the first interface */
      while (other_respair_ptr != NULL && found == FALSE)
        {
          /* Get the two residues */
          other_residue_ptr[0] = other_respair_ptr->residue_ptr[0];
          other_residue_ptr[1] = other_respair_ptr->residue_ptr[1];

          /* If residue name and number match, then have an equivalence */
          if ((equivalence == DIRECT &&
               !strcmp(residue_ptr[0]->res_num,
                       other_residue_ptr[0]->res_num) &&
               !strcmp(residue_ptr[0]->res_name,
                       other_residue_ptr[0]->res_name) &&
               !strcmp(residue_ptr[1]->res_num,
                       other_residue_ptr[1]->res_num) &&
               !strcmp(residue_ptr[1]->res_name,
                       other_residue_ptr[1]->res_name)) ||
              (equivalence == ACROSS &&
               !strcmp(residue_ptr[0]->res_num,
                       other_residue_ptr[1]->res_num) &&
               !strcmp(residue_ptr[0]->res_name,
                       other_residue_ptr[1]->res_name) &&
               !strcmp(residue_ptr[1]->res_num,
                       other_residue_ptr[1]->res_num) &&
               !strcmp(residue_ptr[0]->res_name,
                       other_residue_ptr[0]->res_name)))
            {
              /* Increment count of equivalences */
              nequiv++;

              /* Set flag to end search here */
              found = TRUE;
            }

          /* Get the next link */
          other_respair_ptr = other_respair_ptr->next_respair_ptr;
        }

      /* Get the next link */
      respair_ptr = respair_ptr->next_respair_ptr;
    }

  /* Determine whether the number of ewquivalences constitutes equivalence */
  if (nequiv > nlinks / 2)
    equivalent = TRUE;

  /* Return whether the two given interfaces are equivalent */
  return(equivalent);
}
/***********************************************************************

unique_interfaces  -  Determine which are the unique interfaces and
                      which are copies

***********************************************************************/

void unique_interfaces(struct interface *first_interface_ptr,
                       struct chain *first_chain_ptr,int *nunique,
                       int *niface)
{
  char chain1_1, chain1_2, chain2_1, chain2_2;

  int equiv, equivalence;

  struct interface *interface_ptr, *other_interface_ptr;

  /* Initialise variables */
  nunique[PROTEIN_PROTEIN] = 0;
  nunique[PROTEIN_DNA] = 0;
  niface[PROTEIN_PROTEIN] = 0;
  niface[PROTEIN_DNA] = 0;

  /* Get pointer to the first interface */
  interface_ptr = first_interface_ptr;

  /* Loop through all the interfaces */
  while (interface_ptr != NULL)
    {
      /* Get this interface's two chains */
      chain1_1 = equiv_chain(interface_ptr->chains[0],first_chain_ptr);
      chain1_2 = equiv_chain(interface_ptr->chains[1],first_chain_ptr);

      /* Get pointer to the first interface */
      other_interface_ptr = first_interface_ptr;

      /* Loop through all the previous interfaces to see if this is a
         copy of any one of them */
      while (interface_ptr->unique == TRUE &&
             other_interface_ptr != NULL &&
             other_interface_ptr != interface_ptr)
        {
          /* Only consider this interface if it is a unique one */
          if (other_interface_ptr->unique == TRUE)
            {
              /* Get this interface's chains */
              chain2_1
                = equiv_chain(other_interface_ptr->chains[0],first_chain_ptr);
              chain2_2
                = equiv_chain(other_interface_ptr->chains[1],first_chain_ptr);

              /* Check for equivalence */
              equivalence = NO;
              if (chain1_1 == chain2_1 && chain1_2 == chain2_2)
                equivalence = DIRECT;
              else if (chain1_1 == chain2_2 && chain1_2 == chain2_1)
                equivalence = ACROSS;

              /* If have an equivalence between the chains need
                 to check whether they involve the same groups of
                 residues */
              if (equivalence != NO)
                {
                  equiv = count_pairs(interface_ptr,other_interface_ptr,
                                      equivalence);

                  /* If have sufficient residue-pair equivalences, then
                     can consider interfaces to be equivalent */
                  if (equiv == TRUE)
                    {
                      /* Mark first interface as non-unique and as a copy
                         of the current interface */
                      interface_ptr->unique = FALSE;
                      interface_ptr->copyof_interface_ptr
                        = other_interface_ptr;

                      /* Increment count of copies the unique one has */
                      other_interface_ptr->ncopies++;
                    }
                }
            }

          /* Get next interface record */
          other_interface_ptr = other_interface_ptr->next_interface_ptr;
        }

      /* If interface is a unique one, increment count */
      if (interface_ptr->unique == TRUE)
        nunique[interface_ptr->type]++;

      /* Increment count of all interfaces */
      niface[interface_ptr->type]++;

      /* Get the next interface */
      interface_ptr = interface_ptr->next_interface_ptr;
    }
}
/***********************************************************************

run_naccess  -  Run NACCESS to calculate surface area of given protein
                chain(s)

***********************************************************************/

double run_naccess(struct c_file_data file_name_ptr,char *file_name,
                   char *output_dir,double *area2,char *chains,
                   int pdb_type)
{
  char command[FILENAME_LEN + 1], out_file[FILENAME_LEN + 1];
  char log_file[FILENAME_LEN + 1];
  char chain, input_line[LINELEN + 1], number_string[20];

  char *string_ptr;

  int ichain, len;

  double area;

  FILE *file_ptr;

  /* Initialise variables */
  area = 0.0;
  area2[0] = 0.0;
  area2[1] = 0.0;
  ichain = 0;

  /* Form name of NACCESS log file */
  strcpy(log_file,output_dir);
  strcat(log_file,"naccess.log");

  /* Form the parameters to run the NACCESS program */
  sprintf(command,"%s %s > %s",file_name_ptr.exe_name[NACCESS_EXE],
          file_name,log_file);

  /* For PDBsum1, need to first write the input parameters to an input
     file */
  if (pdb_type == PDBSUM1)
    {
      /* Write out the input parameters */
      sprintf(command,"echo pdbf %s > in",file_name);
      system(command);
      sprintf(command,"echo vdwf vdw.radii >> in");
      system(command);
      sprintf(command,"echo stdf standard.data >> in");
      system(command);
      sprintf(command,"echo prob 1.40 >> in");
      system(command);
      sprintf(command,"echo zsli 0.05 >> in");
      system(command);

      /* Define the command to run naccess */
      sprintf(command,"%s < in >> ../logs/%s",
              file_name_ptr.exe_name[NACCESS_EXE],log_file);
    }

  /* Run the NACCESS program */
  system(command);

  /* Form name of NACCESS output file */
  strcpy(out_file,file_name);
  string_ptr = strstr(out_file,".pdb");
  if (string_ptr != NULL)
    {
      /* Replace .pdb extension with .rsa extension */
      string_ptr[0] = '\0';
      strcat(out_file,".rsa");
    }

  /* Open this file */
  printf("Opening NACCESS file (%s) ...\n",out_file);
  if ((file_ptr = fopen(out_file,"r")) !=  NULL)
    {
      /* Loop while reading in records from the input file */
      while (fgets(input_line,LINELEN,file_ptr)!= NULL)
        {
          /* Truncate the string to strip off the newline */
          len = string_truncate(input_line,LINELEN);

          /* If this is the total area line, retrieve the total surface
             area */
          if (!strncmp(input_line,"TOTAL",5))
            {
              /* Retrieve the total surface area */
              strncpy(number_string,input_line+10,11);
              number_string[11] = '\0';
              area = atof(number_string);
            }

          /* If this is the area of one of the chains, store */
          else if (!strncmp(input_line,"CHAIN",5))
            {
              /* Identify the chain */
              chain = input_line[9];
              if (chain == chains[0])
                ichain = 0;
              else
                ichain = 1;

              /* Retrieve the chain's surface area */
              strncpy(number_string,input_line+10,11);
              number_string[11] = '\0';
              area2[ichain] = atof(number_string);
            }
        }

      /* Close the input file */
      fclose(file_ptr);
    }

  /* Show warning message */
  else
    printf("Unable to open NACCESS output file: %s ...\n",out_file);

  /* Return the surface area */
  return(area);
}
/***********************************************************************

calc_surface_area  -  Use NACCESS to calculate surface area of interface

***********************************************************************/

void calc_surface_area(struct c_file_data file_name_ptr,
                       struct interface *interface_ptr,
                       struct chain *first_chain_ptr,char *output_dir,
                       int pdb_type)
{
  char file_name[FILENAME_LEN], file_name2[FILENAME_LEN];

  int iatom, ichain, iresid, natoms, nout[2], nresid, ntotal;

  double area[2], area2[2], total_area;

  struct atom *atom_ptr;
  struct chain *chain_ptr;
  struct residue *residue_ptr;

  FILE *file_ptr, *file2_ptr;

  /* Initialise variables */
  nout[0] = 0;
  nout[1] = 0;
  ntotal = 0;

  /* Form name of joint chains file */
  strcpy(file_name2,output_dir);
  strcat(file_name2,"ch2.pdb");

  /* Open joint chains file */
  file2_ptr = open_output_file(file_name2,TRUE);

  /* Get pointer to the first chain */
  chain_ptr = first_chain_ptr;

  /* Loop over all the chains */
  while (chain_ptr != NULL)
    {
      /* If this is a protein chain, and is one of the interface chains,
         then process */
      if (chain_ptr->type == PROTEIN &&
          (chain_ptr->chain_id == interface_ptr->chains[0] ||
           chain_ptr->chain_id == interface_ptr->chains[1]))
        {
          /* Determine which one it is */
          if (chain_ptr->chain_id == interface_ptr->chains[0])
            ichain = 0;
          else
            ichain = 1;

          /* Form name of output file */
          strcpy(file_name,output_dir);
          strcat(file_name,"ch.pdb");

          /* Open single chain file */
          file_ptr = open_output_file(file_name,TRUE);

          /* Get the first of this chain's residues */
          residue_ptr = chain_ptr->first_residue_ptr;
          nresid = chain_ptr->nresid;
          iresid = 0;

          /* Loop through the residues to write out */
          while (residue_ptr != NULL && iresid < nresid)
            {
              /* Get number of atoms */
              natoms = residue_ptr->natoms;

              /* Initialise pointer to first atom of this residue */
              atom_ptr = residue_ptr->first_atom_ptr;
              iatom = 0;
  
              /* Loop through all the residue's atoms to write out */
              while (iatom < natoms && atom_ptr != NULL)
                {
                  /* Increment count of atoms written out */
                  nout[ichain]++;
                  ntotal++;

                  /* Write this atom out to the first file */
                  fprintf(file_ptr,"ATOM  %5d %s %s %c%s   %8.3f%8.3f%8.3f",
                          nout[ichain],atom_ptr->atom_name,
                          residue_ptr->res_name,residue_ptr->chain,
                          residue_ptr->res_num,atom_ptr->coords[0],
                          atom_ptr->coords[1],atom_ptr->coords[2]);
                  fprintf(file_ptr,"%6.2f%6.2f\n",atom_ptr->occupancy,
                          atom_ptr->bvalue);

                  /* Write to second file */
                  fprintf(file2_ptr,"ATOM  %5d %s %s %c%s   %8.3f%8.3f%8.3f",
                          ntotal,atom_ptr->atom_name,
                          residue_ptr->res_name,residue_ptr->chain,
                          residue_ptr->res_num,atom_ptr->coords[0],
                          atom_ptr->coords[1],atom_ptr->coords[2]);
                  fprintf(file2_ptr,"%6.2f%6.2f\n",atom_ptr->occupancy,
                          atom_ptr->bvalue);

                  /* Get pointer to next atom in linked-list */
                  atom_ptr = atom_ptr->next_atom_ptr;
                  iatom++;
                }

              /* Get the next residue */
              residue_ptr = residue_ptr->next_residue_ptr;
              iresid++;
            }

          /* Close the single-chain file */
          fclose(file_ptr);

          /* If have already calculated surface area of this chain,
             store area */
          if (chain_ptr->asa > 0.0)
            area[ichain] = chain_ptr->asa;

          /* Otherwise, run NACCESS to compute solvent accessible surface
             area of current chain */
          else
            {
              /* Run NACCESS */
              area[ichain]
                = run_naccess(file_name_ptr,file_name,output_dir,area2,
                              interface_ptr->chains,pdb_type);

              /* Store the computed surface area */
              chain_ptr->asa = area[ichain];
            }
        }

      /* Get the next chain */
      chain_ptr = chain_ptr->next_chain_ptr;
    }

  /* Close the joint chains file */
  fclose(file2_ptr);

  /* If have file containing both chains, then can run NACCESS on the
     complex */
  if (nout[0] > 0 && nout[1] > 0)
    {
      /* Run NACCESS to compute solvent accessible surface area of
         complex of both chains */
      total_area
        = run_naccess(file_name_ptr,file_name2,output_dir,area2,
                      interface_ptr->chains,pdb_type);

      /* Compute the total surface area of each chain */
      interface_ptr->chain_area[0] = area[0];
      interface_ptr->chain_area[1] = area[1];

      /* Calculate the surface area buried by the other chain */
      interface_ptr->interface_area[0] = area[0] - area2[0];
      interface_ptr->interface_area[1] = area[1] - area2[1];
    }
}
/***********************************************************************

write_details  -  Write the current interface's details to the
                  interfaces.dat file

***********************************************************************/

void write_details(FILE *file_ptr,struct interface *interface_ptr,
                   int ncount,char *face_name,int copy,int last_copy)
{
  char item1[10];

  /* Form the subscript for this interface */
  sprintf(item1,"[%d]",ncount);

  /* Write out chains */
  fprintf(file_ptr,"%sIFACE_CHAIN1%s %c\n",face_name,item1,
          interface_ptr->chains[0]);
  fprintf(file_ptr,"%sIFACE_CHAIN2%s %c\n",face_name,item1,
          interface_ptr->chains[1]);

  /* Write out numbers of residues */
  fprintf(file_ptr,"%sIFACE_NRESID1%s %d\n",face_name,item1,
          interface_ptr->nresid[0]);
  fprintf(file_ptr,"%sIFACE_NRESID2%s %d\n",face_name,item1,
          interface_ptr->nresid[1]);

  /* Write out numbers of interactions */
  fprintf(file_ptr,"%sIFACE_NHBONDS%s %d\n",face_name,item1,
          interface_ptr->count[HBONDS]);
  fprintf(file_ptr,"%sIFACE_NCONTACTS%s %d\n",face_name,item1,
          interface_ptr->count[CONTACTS]);
  fprintf(file_ptr,"%sIFACE_SALT_BRIDGES%s %d\n",face_name,item1,
          interface_ptr->count[SALT_BRIDGES]);
  fprintf(file_ptr,"%sIFACE_DISULPHIDES%s %d\n",face_name,item1,
          interface_ptr->count[DISULPHIDES]);

  /* Write out total surface area of each chain */
  fprintf(file_ptr,"%sIFACE_AREA1%s %8.1f\n",face_name,item1,
          interface_ptr->chain_area[0]);
  fprintf(file_ptr,"%sIFACE_AREA2%s %8.1f\n",face_name,item1,
          interface_ptr->chain_area[1]);

  /* Write out surface area of interaction */
  fprintf(file_ptr,"%sIFACE_BURIED1%s %8.1f\n",face_name,item1,
          interface_ptr->interface_area[0]);
  fprintf(file_ptr,"%sIFACE_BURIED2%s %8.1f\n",face_name,item1,
          interface_ptr->interface_area[1]);

  /* In interface has any copies, write out count */
  if (interface_ptr->ncopies > 0)
    fprintf(file_ptr,"%sIFACE_NCOPY%s %d\n",face_name,item1,
            interface_ptr->ncopies);

  /* If this is a copy chain, write out flag */
  if (copy == TRUE)
    fprintf(file_ptr,"%sIFACE_COPY%s TRUE\n",face_name,item1);

  /* If this is the last copy chain, write out flag */
  if (last_copy == TRUE)
    fprintf(file_ptr,"%sIFACE_ENDCOPY%s TRUE\n",face_name,item1);
}
/***********************************************************************

write_interfaces_dat  -  Write the interfaces.dat file giving database
                         info on this structure's interfaces

***********************************************************************/

int write_interfaces_dat(char *file_name,struct c_file_data file_name_ptr,
                         struct interface *first_interface_ptr,
                         struct interface ***interface_idx_ptr,
                         int *nunique,int *niface,
                         struct chain *first_chain_ptr,char *output_dir,
                         int pdb_type)
{
  static char *interface_name[] = { "PPINTERFACE", "PDINTERFACE" };
  static char *face_name[] = { "PP", "PD" };

  int interface_type, iunique, ncount, ntotal;
  int i, ncopy, nindex;
  int last_copy;

  struct interface *interface_ptr, *other_interface_ptr;

  struct interface **interface_index_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  nindex = 0;

  /* Get the total number of unique interfaces */
  // ntotal = nunique[PROTEIN_PROTEIN] + nunique[PROTEIN_DNA];
  ntotal = niface[PROTEIN_PROTEIN] + niface[PROTEIN_DNA];

  /* Grab memory for interface index array and initialise */
  interface_index_ptr
    = (struct interface **) malloc(ntotal * sizeof(struct interface *));
  for (i = 0; i < ntotal; i++)
    interface_index_ptr[i] = NULL;

  /* Open the  output file */
  file_ptr = open_output_file(file_name,TRUE);

  /* Write out primary flag */
  if (nunique[PROTEIN_PROTEIN] > 0)
    fprintf(file_ptr,"INTERFACES TRUE\n");

  /* Write out total count of protein-protein interfaces */
  if (niface[PROTEIN_PROTEIN] > 0)
    fprintf(file_ptr,"NINTERFACES %d\n",niface[PROTEIN_PROTEIN]);

  /* Loop over the two interface types */
  for (interface_type = 0; interface_type < 2; interface_type++)
    {
      /* Write out the number of unique interfaces */
      if (nunique[interface_type] > 0)
        fprintf(file_ptr,"NUNIQUE_%sS %d\n",interface_name[interface_type],
                nunique[interface_type]);

      /* Get pointer to first interface record */
      interface_ptr = first_interface_ptr;

      /* Loop through interfaces to calculate their surface accessible
         areas */
      while (interface_ptr != NULL)
        {
          /* If this is an interface of the right type, run NACCESS */
          if (interface_ptr->type == interface_type)
            {
              /* Use NACCESS to calculate surface area of interface */
              calc_surface_area(file_name_ptr,interface_ptr,
                                first_chain_ptr,output_dir,pdb_type);
            }

          /* Get next interface record */
          interface_ptr = interface_ptr->next_interface_ptr;
        }

      /* Get pointer to first interface record */
      interface_ptr = first_interface_ptr;
      iunique = 0;
      ncount = 0;

      /* Loop through interface records to write out all unique interfaces
         and their copies */
      while (interface_ptr != NULL)
        {
          /* If this is a unique interface of the right type, write out */
          if (interface_ptr->unique == TRUE &&
              interface_ptr->type == interface_type)
            {
              /* Increment count of unique interfaces */
              iunique++;

              /* Save pointer in index */
              if (nindex < ntotal)
                {
                  interface_index_ptr[nindex] = interface_ptr;
                  nindex++;
                }

              /* Increment count of interfaces written out */
              ncount++;

              /* Write the details to the interfaces.dat file */
              write_details(file_ptr,interface_ptr,ncount,
                            face_name[interface_type],FALSE,FALSE);

              /* In interface has any copies, write out */
              if (interface_ptr->ncopies > 0)
                {
                  /* Get pointer to next interface record */
                  other_interface_ptr = interface_ptr->next_interface_ptr;
                  ncopy = 0;

                  /* Loop through all the other interfaces to write out
                     the copies */
                  while (other_interface_ptr != NULL)
                    {
                      /* If this is a copy of the current interface,
                         write out */
                      if (other_interface_ptr->copyof_interface_ptr
                          == interface_ptr)
                        {
                          /* Increment copy counter */
                          ncopy++;

                          /* If this is the last copy, set appropriate
                             flag */
                          if (ncopy == interface_ptr->ncopies)
                            last_copy = TRUE;
                          else
                            last_copy = FALSE;

                          /* Save pointer in index */
                          if (nindex < ntotal)
                            {
                              interface_index_ptr[nindex]
                                = other_interface_ptr;
                              nindex++;
                            }

                          /* Increment count of interfaces written out */
                          ncount++;

                          /* Write the details to the interfaces.dat file */
                          write_details(file_ptr,other_interface_ptr,ncount,
                                        face_name[interface_type],TRUE,
                                        last_copy);

                        }

                      /* Get next interface record */
                      other_interface_ptr
                        = other_interface_ptr->next_interface_ptr;
                    }
                }
            }

          /* Get next interface record */
          interface_ptr = interface_ptr->next_interface_ptr;
        }
    }

  /* Close the output file */
  fclose(file_ptr);

  /* Return index of unique interfaces */
  *interface_idx_ptr = interface_index_ptr;

  /* Return number of unique interfaces */
  return(nindex);
}
/***********************************************************************

calc_iface_transformation  -  Compute transformation to align interface
                              interaction points

***********************************************************************/

void calc_iface_transformation(struct respair *first_respair_ptr,
                               double transformation[3][3],double cofg[3])
{
  char stored_metal[NMETALS][3];

  static char atom_name[] = " CA ";

  int i, j, natoms, nresid;
  int metal_count[NMETALS], nmetal_types, nmetals;

  static int rec_type = ATOM_REC;

  double coord_min[3], coord_max[3];

  static double bvalue = 10.0;
  static double occupancy = 1.0;

  struct atom *atom_ptr, *first_atom_ptr, *last_atom_ptr;
  struct chain *dummy_chain_ptr, *first_chain_ptr, *last_chain_ptr;
  struct residue *residue_ptr, *first_residue_ptr, *last_residue_ptr;
  struct residue *dummy_residue_ptr;
  struct respair *respair_ptr;

  /* Initialise variables */
  natoms = nresid = 0;

  /* Initialise atom pointers */
  first_atom_ptr = NULL;
  last_atom_ptr = NULL;

  /* Create a dummy chain record */
  first_chain_ptr = last_chain_ptr = NULL;
  dummy_chain_ptr
    = create_chain_record(&first_chain_ptr,&last_chain_ptr,'X');
  
  /* Create a dummy residue record */
  first_residue_ptr = last_residue_ptr = NULL;
  nresid++;
  dummy_residue_ptr
    = create_residue_record(&first_residue_ptr,&last_residue_ptr,
                            dummy_chain_ptr,"DUM","   1 ",'X',nresid);

  /* Get pointer to the first link */
  respair_ptr = first_respair_ptr;

  /* Loop over all the links to create dujmmy atoms for each averaged
     residue interaction point */
  while (respair_ptr != NULL)
    {
      /* Loop over the two residues making up this interaction */
      for (i = 0; i< 2; i++)
        {
          /* Get this residue */
          residue_ptr = respair_ptr->residue_ptr[i];

          /* Update this residue's average interaction position */
          for (j = 0; j < 3; j++)
            residue_ptr->ave_iface_coords[j]
              = residue_ptr->ave_iface_coords[j]
              + respair_ptr->ave_coords[i][j];
          residue_ptr->niface_atoms
            = residue_ptr->niface_atoms + respair_ptr->natoms[i];

          /* Calculate average position for this residue's interactions */
          if (respair_ptr->natoms[i] > 0)
            for (j = 0; j < 3; j++)
              respair_ptr->ave_coords[i][j]
                = respair_ptr->ave_coords[i][j] / respair_ptr->natoms[i];

          /* Create a dummy record to represent this point */
          atom_ptr
            = create_atom_record(&first_atom_ptr,&last_atom_ptr,
                                 dummy_residue_ptr,natoms,atom_name,
                                 respair_ptr->ave_coords[i][0],
                                 respair_ptr->ave_coords[i][1],
                                 respair_ptr->ave_coords[i][2],
                                 bvalue,occupancy,stored_metal,
                                 metal_count,&nmetal_types,&nmetals,
                                 rec_type);

          /* Increment count of dummy atoms */
          natoms++;
        }

      /* Get pointer to the next link */
      respair_ptr = respair_ptr->next_respair_ptr;
    }

  /* Use principal components to align residue along its
     principal axes */
  align_molecule(first_atom_ptr,natoms,coord_min,coord_max,
                 transformation,cofg);

  /* Free up the memory taken up by the dumm atoms */
  free_atom_memory(first_atom_ptr);
}
/***********************************************************************

write_dimhtml_file  -  Write out the dimhtmlXY.dat file for plotting of
                       the interface residues and their interactions by
                       the dinhtml program

***********************************************************************/

void write_dimhtml_file(char *output_dir,struct interface *interface_ptr,
                        struct residue *first_residue_ptr)
{
  char file_name[FILENAME_LEN + 1];

  int i, icoord, nresid;

  double cofg[3], coords[3], transformation[3][3], x, y, z;

  struct residue *residue_ptr;
  struct respair *respair_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  nresid = 0;

  /* Get the first residues */
  residue_ptr = first_residue_ptr;

  /* Loop through all the residues to initialise ref numbers */
  while (residue_ptr != NULL)
    {
      /* Initialise reference number */
      residue_ptr->ref = -1;

      /* Initialise average interacting coords */
      for (i = 0; i < 3; i++)
        residue_ptr->ave_iface_coords[i] = 0.0;
      residue_ptr->niface_atoms = 0;

      /* Get the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
    }

  /* Compute principal component transformation on all the average
     interaction positions for this interface */
  calc_iface_transformation(interface_ptr->first_respair_ptr,
                            transformation,cofg);

  /* Form name of output dimhtmlXY.dat file */
  strcpy(file_name,output_dir);
  strcat(file_name,"dimhtml");
  strcat(file_name,interface_ptr->chains);
  strcat(file_name,".dat");

  /* Open the output file */
  file_ptr = open_output_file(file_name,TRUE);

  /* Get pointer to the first link */
  respair_ptr = interface_ptr->first_respair_ptr;

  /* Loop over all the links to write out the residus involved */
  while (respair_ptr != NULL)
    {
      /* Loop over the two residues making up this interaction */
      for (i = 0; i < 2; i++)
        {
          /* Get this residue */
          residue_ptr = respair_ptr->residue_ptr[i];

          /* If residue hasn't already been written out, assign a
             reference number and write out */
          if (residue_ptr->ref == -1)
            {
              /* Assign reference number */
              residue_ptr->ref = nresid;

              /* Transform this residue's average interacting coords */
              if (residue_ptr->niface_atoms > 0)
                {
                  /* Average the coords */
                  x = residue_ptr->ave_iface_coords[0]
                    / residue_ptr->niface_atoms;
                  y = residue_ptr->ave_iface_coords[1]
                    / residue_ptr->niface_atoms;
                  z = residue_ptr->ave_iface_coords[2]
                    / residue_ptr->niface_atoms;
                }
              else
                x = y = z = 0.0;

              /* Perform CofG transformation */
              x = x - cofg[0];
              y = y - cofg[1];
              z = z - cofg[2];

              /* Apply interface transformation */
              for (icoord = 0; icoord < 3; icoord++)
                {
                  /* Calculate transformed coordinates */
                  coords[icoord]
                    = ((x * transformation[0][icoord])
                       + (y * transformation[1][icoord])
                       + (z * transformation[2][icoord]));
                }

              /* Write this residue out */
              fprintf(file_ptr,"%s %s %c %d %8.3f\n",residue_ptr->res_name,
                      residue_ptr->res_num,residue_ptr->chain,
                      residue_ptr->ref,coords[0]);

              /* Increment count of used residues */
              nresid++;
            }
        }

      /* Get pointer to the next link */
      respair_ptr = respair_ptr->next_respair_ptr;
    }

  /* Write heading to signify interaction details to follow */
  fprintf(file_ptr,"INTERACTIONS\n");

  /* Get pointer to the first link again */
  respair_ptr = interface_ptr->first_respair_ptr;

  /* Loop through all the interactions to write out */
  while (respair_ptr != NULL)
    {
      /* Get the first residue and write out */
      residue_ptr = respair_ptr->residue_ptr[0];
      fprintf(file_ptr,"%d",residue_ptr->ref);

      /* Repeat for second residue */
      residue_ptr = respair_ptr->residue_ptr[1];
      fprintf(file_ptr," %d",residue_ptr->ref);

      /* Write out the counts of interactions between these two
         residues */
      for (i = 0; i < 2; i++)
        fprintf(file_ptr," %d",respair_ptr->count[i]);
      if (respair_ptr->count[SALT_BRIDGES] > 0)
        fprintf(file_ptr," s");
      if (respair_ptr->count[DISULPHIDES] > 0)
        fprintf(file_ptr," d");
      fprintf(file_ptr,"\n");

      /* Get the next link */
      respair_ptr = respair_ptr->next_respair_ptr;
    }

  /* Close the output file */
  fclose(file_ptr);
}
/***********************************************************************

save_interface_data  -  Write out the interface datafiles

***********************************************************************/

void save_interface_data(char *interfaces_dat,
                         struct c_file_data file_name_ptr,
                         struct interface **fst_interface_ptr,
                         int ninterfaces,
                         struct chain *first_chain_ptr,
                         struct residue *first_residue_ptr,
                         char *output_dir,int pdb_type)
{
  int i, niface[2], nindex, nunique[2];

  struct interface *first_interface_ptr;

  struct interface **interface_index_ptr;

  /* Get the first interface pointer */
  first_interface_ptr = *fst_interface_ptr;

  /* Sort the interfaces alphabetically */
  sort_interfaces(&first_interface_ptr,ninterfaces);

  /* Determine which are the unique interfaces and which are the copies */
  unique_interfaces(first_interface_ptr,first_chain_ptr,nunique,niface);

  /* Write the interfaces.dat file */
  nindex
    = write_interfaces_dat(interfaces_dat,file_name_ptr,first_interface_ptr,
                           &interface_index_ptr,nunique,niface,
                           first_chain_ptr,output_dir,pdb_type);

  /* Loop over all interfaces to write out the residue files for each */
  for (i = 0; i < nindex; i++)
    write_dimhtml_file(output_dir,interface_index_ptr[i],
                       first_residue_ptr);

  /* Return first interface pointer */
  *fst_interface_ptr = first_interface_ptr;
}
/***********************************************************************

write_pdb_headers  -  Write out the header records for the pymol.pdb
                      file

***********************************************************************/

void write_pdb_headers(FILE *file_ptr,char *pdb_code,
                       struct chain *first_chain_ptr,int biomol)
{
  char chain_name[5], code[5], compnd_no[5], date[64];

  int ncompnd, nmol, type;

  struct chain *chain_ptr;

  /* Initialise variables */
  ncompnd = 0;
  nmol = 0;

  /* Get the current date */
  get_date(date);

  /* Convert PDB code to upper case */
  strcpy(code,pdb_code);
  to_upper(code,4);

  /* Write out the PDB header records */
  fprintf(file_ptr,
          "HEADER    "
          "PDB STRUCTURE                           %s   %s\n"
          "TITLE     ",date,code);
  if (biomol == TRUE)
    fprintf(file_ptr,"BIOMOLECULE TAKEN FROM PQS DIRECTORY\n");
  else
    fprintf(file_ptr,
            "AUTO-GENERATED PDB FILE FOR GENERATING PYMOL IMAGES     \n");

  /* Get pointer to the first chain */
  chain_ptr = first_chain_ptr;

  /* Loop over all the chains to write out the molecule-id records */
  while (chain_ptr != NULL)
    {
      /* Get the chain type */
      type = chain_ptr->type;

      /* If protein or DNA, write out the MOL-ID records */
      if (type == CALPHA_ONLY || type == PROTEIN || type ==  DNA)
        {
          /* Get the chain identifier for this record */
          if (chain_ptr->chain_id == ' ')
            strcpy(chain_name,"NULL");
          else
            {
              chain_name[0] = chain_ptr->chain_id;
              chain_name[1] = '\0';
            }

          /* Increment molecule counter */
          nmol++;

          /* Write the records */
          ncompnd++;
          sprintf(compnd_no,"%4d",ncompnd);
          if (ncompnd == 1)
            strcpy(compnd_no,"   ");
          fprintf(file_ptr,"COMPND%s MOL_ID: %d;\n",compnd_no,nmol);
          ncompnd++;
          fprintf(file_ptr,"COMPND%4d MOLECULE: %s;\n",ncompnd,
                  CHAIN_TYPE[type]);
          ncompnd++;
          fprintf(file_ptr,"COMPND%4d CHAIN: %s;\n",ncompnd,chain_name);
        }

      /* Get the next chain */
      chain_ptr = chain_ptr->next_chain_ptr;
    }
}
/***********************************************************************

write_pdb_secstruc  -  Write out the secondary structure assignments to
                       the pymol.pdb file

***********************************************************************/

void write_pdb_secstruc(FILE *file_ptr,struct chain *first_chain_ptr,
                        struct ssbond *first_ssbond_ptr,int nssbonds)
{
  int i, iresid, loop, nresid, nsse, nstresid, sec_struc, sst, type;

  float percent;

  struct chain *chain_ptr;
  struct residue *residue_ptr, *last_residue_ptr, *start_residue_ptr;
  struct residue *residue1_ptr, *residue2_ptr;
  struct ssbond *ssbond_ptr;

  /* Get pointer to the first chain again */
  chain_ptr = first_chain_ptr;

  /* Loop over any protein chains to write out the secondary structure
     assignments */
  while (chain_ptr != NULL)
    {
      /* Get the chain type */
      type = chain_ptr->type;

      /* If protein chain, write out secondary structure */
      if (type == PROTEIN)
        {
          /* Loop twice to print out the helix and strand assignments */
          for (loop = 0; loop < 2; loop++)
            {
              /* Initialise variables */
              start_residue_ptr = NULL;
              last_residue_ptr = NULL;
              nsse = 0;
              nstresid = 0;
              sec_struc = loop + 1;

              /* Get this chain's first residue */
              residue_ptr = chain_ptr->first_residue_ptr;
              nresid = chain_ptr->nresid;
              iresid = 0;

              /* Loop through this chain's residues, to find start and
                 end of the current secondary structure type */
              while (residue_ptr != NULL && iresid < nresid)
                {
                  /* Get this residue's secondary structure assignment */
                  if (residue_ptr->link_residue_ptr == NULL)
                    sst = residue_ptr->sec_struc;
                  else
                    sst = residue_ptr->link_residue_ptr->sec_struc;

                  /* Check if this residue belongs to the required secondary
                     structure */
                  if (sst == sec_struc)
                    {
                      /* If this is the first residue, store its pointer */
                      if (nstresid == 0)
                        start_residue_ptr = residue_ptr;

                      /* Increment residue count */
                      nstresid++;

                      /* Increment number of residues belonging to this
                         secondary structure */
                      chain_ptr->nres_sst[sec_struc]++;
                    }

                  /* Update coil count */
                  else if (loop == 0)
                    chain_ptr->nres_sst[COIL]++;

                  /* If not the right secondary structure, or this is
                     the last residue, see if need to write out previous
                     secondary structure */
                  if (sst != sec_struc || iresid == nresid - 1)
                    {
                      /* If this is the last residue and belongs to the
                         current secondary, then save its pointer */
                      if (iresid == nresid - 1 && sst != sec_struc)
                        last_residue_ptr = residue_ptr;

                      /* If have a secondary structure element to
                         write out, then do so */
                      if (nstresid > 0)
                        {
                          /* Increment count of secondary structure
                             elements done */
                          nsse++;

                          /* Write out */
                          if (sec_struc == HELIX)
                            {
                              fprintf(file_ptr,"HELIX %4d%4d ",nsse,nsse);
                              fprintf(file_ptr,"%s %c %s ",
                                      start_residue_ptr->res_name,
                                      start_residue_ptr->chain,
                                      start_residue_ptr->res_num);
                              fprintf(file_ptr,"%s %c %s",
                                      last_residue_ptr->res_name,
                                      last_residue_ptr->chain,
                                      last_residue_ptr->res_num);
                              fprintf(file_ptr,"%2d\n",1);
                            }
                          else
                            {
                              /* Write out the strand start and end
                                 residues */
                              fprintf(file_ptr,"SHEET %4d%4d 1 ",nsse,nsse);
                              fprintf(file_ptr,"%s %c%s ",
                                      start_residue_ptr->res_name,
                                      start_residue_ptr->chain,
                                      start_residue_ptr->res_num);
                              fprintf(file_ptr,"%s %c%s ",
                                      last_residue_ptr->res_name,
                                      last_residue_ptr->chain,
                                      last_residue_ptr->res_num);
                              fprintf(file_ptr,"%2d\n",0);
                            }

                          /* Re-initialise for next one */
                          start_residue_ptr = NULL;
                          nstresid = 0;
                        }
                    }

                  /* Save the current residue pointer */
                  last_residue_ptr = residue_ptr;

                  /* Get the next residue */
                  residue_ptr = residue_ptr->next_residue_ptr;
                  iresid++;
                }
            }

          /* Determine which secondary structure type is the dominant
             one in this chain */
          if (chain_ptr->nres_sst[HELIX] + chain_ptr->nres_sst[STRAND] > 0)
            {
              /* Calculate the percentage of helix */
              percent
                = 100.0 * (float) chain_ptr->nres_sst[HELIX]
                / (float) (chain_ptr->nres_sst[HELIX]
                           + chain_ptr->nres_sst[STRAND]);

              /* If strands are "dominant", then set flag in the chain
                 record */
              if (percent < MIN_HELIX &&
                  chain_ptr->nres_sst[STRAND] > MIN_STRAND)
                chain_ptr->dominant_sst = STRAND;
            }
        }

      /* Get the next chain */
      chain_ptr = chain_ptr->next_chain_ptr;
    }

  /* Get the first disulphide bond */
  ssbond_ptr = first_ssbond_ptr;
  nsse = 0;

  /* Loop through the disulphide bonds to write the SSBOND records */
  while (ssbond_ptr != NULL)
    {
      /* Get the two residues */
      residue1_ptr = ssbond_ptr->residue1_ptr;
      residue2_ptr = ssbond_ptr->residue2_ptr;

      /* Increment count */
      nsse++;

      /* Write out this disulphide bond */
      fprintf(file_ptr,"SSBOND%4d %s %c %s   %s %c %s\n",nsse,
              residue1_ptr->res_name,residue1_ptr->chain,
              residue1_ptr->res_num,
              residue2_ptr->res_name,residue2_ptr->chain,
              residue2_ptr->res_num);

      /* Flag residues as belonging to a disulphide bond */
      residue1_ptr->disulphide = TRUE;
      residue2_ptr->disulphide = TRUE;

      /* Get the next ssbond */
      ssbond_ptr = ssbond_ptr->next_ssbond_ptr;
    }
}
/***********************************************************************

apply_transformation  -  Apply the given transformation to all the
                         atoms of the current PDB file

***********************************************************************/

void apply_transformation(float *new_coords,double *coords,
                          float matrix[4][3])
{
  /* Apply the rotation */
  new_coords[0]
    = coords[0] * matrix[0][0]
    + coords[1] * matrix[0][1]
    + coords[2] * matrix[0][2];
  new_coords[1]
    = coords[0] * matrix[1][0]
    + coords[1] * matrix[1][1]
    + coords[2] * matrix[1][2];
  new_coords[2]
    = coords[0] * matrix[2][0]
    + coords[1] * matrix[2][1]
    + coords[2] * matrix[2][2];

  /* Apply translation */
  new_coords[0] = new_coords[0] + matrix[3][0];
  new_coords[1] = new_coords[1] + matrix[3][1];
  new_coords[2] = new_coords[2] + matrix[3][2];
}
/***********************************************************************

write_pdb_coords  -  Write the atom coords to the pymol.pdb file

***********************************************************************/

void write_pdb_coords(FILE *file_ptr,struct chain *first_chain_ptr,
                      int biomol)
{
  char atom_type[7], number_string[20], res_num[6];

  int i, iatom, iresid, natoms, nresid, type;
  int atom_number, residue_number;
  int wanted;

  float bvalue, coords[3], occupancy;

  struct atom *atom_ptr;
  struct chain *chain_ptr;
  struct residue *residue_ptr;

  /* Initialise variables */
  atom_number = 0;
  bvalue = 10.0;
  occupancy = 1.0;

  /* Get pointer to the first chain again */
  chain_ptr = first_chain_ptr;

  /* Loop over all chains to write out atomm coords */
  while (chain_ptr != NULL)
    {
      /* Get the chain type */
      type = chain_ptr->type;
      wanted = FALSE;

      /* If protein chain, write out secondary structure */
      if (type == PROTEIN || type == DNA || type == CALPHA_ONLY ||
          type == ATTACHMENT)
        {
          strcpy(atom_type,"ATOM  ");
          wanted = TRUE;
        }
      else if (type == LIGAND || type == METAL)
        {
          strcpy(atom_type,"HETATM");
          wanted = TRUE;
        }

      /* If writing out biomolecule, want all atoms */
      if (biomol == TRUE)
        wanted = TRUE;

      /* If chain is wanted, loop over required number of transformations
         to write out its atoms */
      if (wanted == TRUE)
        {
          /* Get pointer to the first residue record */
          residue_ptr = chain_ptr->first_residue_ptr;
          nresid = chain_ptr->nresid;
          iresid = 0;

          /* Loop through this chain's residues, to write out */
          while (residue_ptr != NULL && iresid < nresid)
            {
              /* Get the first of this residue's atoms */
              atom_ptr = residue_ptr->first_atom_ptr;
              natoms = residue_ptr->natoms;
              iatom = 0;

              /* Get residue number, adjusting if this is a
                 symmetry-related copy of the original */
              strcpy(res_num,residue_ptr->res_num);

              /* Loop over this residue's atoms to write out */
              while (atom_ptr != NULL && iatom < natoms)
                {
                  /* Determine whether this atom is wanted */
                  wanted = TRUE;
                  if (type == PROTEIN &&
                      strcmp(atom_ptr->atom_name," N  ") &&
                      strcmp(atom_ptr->atom_name," CA ") &&
                      strcmp(atom_ptr->atom_name," C  ") &&
                      strcmp(atom_ptr->atom_name," O  "))
                    wanted = FALSE;

                  /* If residue belongs to a disulphide bond, or is one of
                     the SITE residues, we want all its atoms */
                  if (residue_ptr->disulphide == TRUE ||
                      residue_ptr->site[0] != '\0')
                    wanted = TRUE;

                  /* If this is a biomolecule, all atoms are wanted */
                  if (biomol == TRUE)
                    {
                      /* Set flag that wanted */
                      wanted = TRUE;

                      /* Get residue type */
                      type = residue_ptr->chain_ptr->type;
                      if (type == LIGAND || type == METAL)
                        strcpy(atom_type,"HETATM");
                      else
                        strcpy(atom_type,"ATOM  ");
                    }

                  /* If atom is wanted, write out */
                  if (wanted == TRUE)
                    {
                      /* Increment atom count */
                      atom_number++;
                      if (atom_number > 99999)
                        atom_number = 1;

                      /* Get this atoms coords */
                      for (i = 0; i < 3; i++)
                        coords[i] = atom_ptr->coords[i];

                      /* Write out the details */
                      fprintf(file_ptr,"%s%5d %s %s %c%s   %8.3f%8.3f%8.3f",
                              atom_type,atom_number,atom_ptr->atom_name,
                              residue_ptr->res_name,
                              residue_ptr->chain,res_num,coords[0],
                              coords[1],coords[2]);
                      fprintf(file_ptr,"%6.2f%6.2f",occupancy,bvalue);
                      fprintf(file_ptr,"          %s\n",atom_ptr->element);
                    }

                  /* Get pointer to next atom record and increment atom
                     counter */
                  atom_ptr = atom_ptr->next_atom_ptr;
                  iatom++;
                }

              /* Get the next residue */
              residue_ptr = residue_ptr->next_residue_ptr;
              iresid++;
            }

          /* Write out a TER record */
          fprintf(file_ptr,"TER\n");
        }

      /* Get the next chain */
      chain_ptr = chain_ptr->next_chain_ptr;
    }
}
/***********************************************************************

write_pymol_pdb  -  Create a fresh version of the PDB file, called
                    pymol.pdb, incorporating all the secondary structure
                    definitions

***********************************************************************/

void write_pymol_pdb(char *pdb_code,struct chain *first_chain_ptr,
                     struct ssbond *first_ssbond_ptr,int nssbonds,
                     char *file,int biomol)
{
  char file_name[FILENAME_LEN + 1];

  FILE *file_ptr;

  /* Form name of pymol.pdb file */
  strcpy(file_name,file);

  /* Open the output file */
  file_ptr = open_output_file(file_name,TRUE);

  /* Write out the PDB header records */
  write_pdb_headers(file_ptr,pdb_code,first_chain_ptr,biomol);

  /* Write out the secondary structure assignments */
  write_pdb_secstruc(file_ptr,first_chain_ptr,first_ssbond_ptr,
                     nssbonds);

  /* Write out all the atom coordinates */
  write_pdb_coords(file_ptr,first_chain_ptr,biomol);

  /* Close the output pymol.pdb file */
  fclose(file_ptr);
}
/***********************************************************************

write_sst_dat  -  Write out the sst.dat files giving residue counts for
                  each secondary structure type

***********************************************************************/

void write_sst_dat(struct chain *first_chain_ptr,char *sst_dat)
{
  struct chain *chain_ptr;

  FILE *file_ptr;

  /* Open the output sst.dat file */
  file_ptr = open_output_file(sst_dat,TRUE);

  /* Get pointer to first chain */
  chain_ptr = first_chain_ptr;

  /* Loop over all chains to write out */
  while (chain_ptr != NULL)
    {
      /* Only process if this is a protein sequence */
      if (chain_ptr->type == PROTEIN || chain_ptr->type == CALPHA_ONLY)
        {
          /* Write out the chain and its secondary structure counts */
          fprintf(file_ptr,"%c %d %d %d ",chain_ptr->chain_id,
                  chain_ptr->nres_sst[COIL],chain_ptr->nres_sst[HELIX],
                  chain_ptr->nres_sst[STRAND]);
          if (chain_ptr->dominant_sst == HELIX)
            fprintf(file_ptr,"H");
          else
            fprintf(file_ptr,"S");
          fprintf(file_ptr,"\n");
        }

      /* Get pointer to next chain */
      chain_ptr = chain_ptr->next_chain_ptr;
    }

  /* Close the output file */
  fclose(file_ptr);
}
/***********************************************************************

create_pobject_record  -  Create and initialise a new link record

***********************************************************************/

struct pobject *create_pobject_record(struct pobject **fst_pobject_ptr,
                                      struct pobject **lst_pobject_ptr,
                                      char *name,int type,char chain_id)
{
  struct pobject *pobject_ptr, *first_pobject_ptr, *last_pobject_ptr;

  /* Initialise pointers */
  pobject_ptr = NULL;
  first_pobject_ptr = *fst_pobject_ptr;
  last_pobject_ptr = *lst_pobject_ptr;

  /* Allocate memory for structure to hold link info */
  pobject_ptr = (struct pobject *) malloc(sizeof(struct pobject));
  if (pobject_ptr == NULL)
    {
      printf("*** Can't allocate memory for struct pobject\n");
      exit (1);
    }

  /* If this is the very first link, then store pointer */
  if (first_pobject_ptr == NULL)
    first_pobject_ptr = pobject_ptr;

  /* Add link from previous link to the current one */
  if (last_pobject_ptr != NULL)
    last_pobject_ptr->next_pobject_ptr = pobject_ptr;
  last_pobject_ptr = pobject_ptr;

  /* Store the known details */
  store_string(&(pobject_ptr->name),name);
  pobject_ptr->type = type;
  pobject_ptr->chain_id = chain_id;

  /* Initialise remaining fields */
  pobject_ptr->display_type = &null;
  pobject_ptr->ligmet_type = 0;
  pobject_ptr->got_sites = FALSE;
  pobject_ptr->got_disulph = FALSE;

  /* Initialise pointer to next record */
  pobject_ptr->next_pobject_ptr = NULL;

  /* Return current pointers */
  *fst_pobject_ptr = first_pobject_ptr;
  *lst_pobject_ptr = last_pobject_ptr;

  /* Return new link pointer */
  return(pobject_ptr);
}
/***********************************************************************

write_pymol_chains  -  Create PyMol objects for each protein and DNA
                       chain in the structure

***********************************************************************/

void write_pymol_chains(struct chain *first_chain_ptr,
                        FILE *pymol_script_ptr,
                        struct pobject **fst_pobject_ptr,
                        struct pobject **lst_pobject_ptr,
                        int colour_as_orig,int transparent_copies,
                        int biomol,int pdb_type)
{
  char chain, object[TOKEN_LEN + 1], chain_name[TOKEN_LEN + 1];
  char colour_name[COLNAM_LEN], display_type[20];
  char colour_major[COLNAM_LEN + 6],  colour_minor[COLNAM_LEN + 6];
  char colour_coil[COLNAM_LEN + 5], number_string[20];
  char name[3], res_num1[6], res_num2[6], residue_range[13];

  int colour, iresid, metal_type, nresid, nsite, number, type, valid;
  int first_res, last_res, last_range;
  int transparency;
  int split;

  double rgb[3], rgb_hel[3], rgb_str[3];

  struct chain *chain_ptr;
  struct pobject *pobject_ptr, *first_pobject_ptr, *last_pobject_ptr;
  struct residue *residue_ptr;

  /* Initialise variables */
  split = TRUE;
  transparency = PARTIAL;

  /* Initialise pointers */
  first_pobject_ptr = *fst_pobject_ptr;
  last_pobject_ptr = *lst_pobject_ptr;

  /* Get pointer to the first chain */
  chain_ptr = first_chain_ptr;

  /* Loop through all the chains to define the various objects */
  while (chain_ptr != NULL)
    {
      /* Get chain type */
      type = chain_ptr->type;

      /* If protein or DNA, create a name for this chain */
      if (type == PROTEIN || type == DNA || type == CALPHA_ONLY)
        {
          /* Show comment */
          fprintf(pymol_script_ptr,"# Definition for chain [%c] - ",
                  chain_ptr->chain_id);
          if (type == PROTEIN)
            fprintf(pymol_script_ptr,"PROTEIN\n");
          else if (type == DNA)
            fprintf(pymol_script_ptr,"DNA\n");
          else
            fprintf(pymol_script_ptr,"CALPHA_ONLY\n");

          /* Initialise the first residue and last residues */
          res_num1[0] = '\0';
          res_num2[0] = '\0';

          /* Get this chain's first residue */
          residue_ptr = chain_ptr->first_residue_ptr;
          nresid = chain_ptr->nresid;
          iresid = 0;
          valid = FALSE;

          /* Loop through this chain's residues, to find start and
             end residues */
          while (residue_ptr != NULL && iresid < nresid)
            {
              /* If this is the first residue, save number */
              if (res_num1[0] == '\0')
                strcpy(res_num1,residue_ptr->res_num);

              /* Save number as last residue */
              strcpy(res_num2,residue_ptr->res_num);

              /* Get the next residue */
              residue_ptr = residue_ptr->next_residue_ptr;
              iresid++;
            }

          /* If have a residue range, form name of range */
          if (res_num1[0] != '\0' && res_num2[0] != '\0')
            {
              /* Define the residue range */
              sprintf(residue_range,"%s-%s",res_num1,res_num2);

              /* Remove blanks from the range */
              valid = remove_residue_blanks(residue_range);
            }

          /* If this is chain blank, define it in terms of its residue
             range */
          if (chain_ptr->chain_id == ' ')
            {
              /* Define the object name as chain_ */
              strcpy(object,"chain_");

              /* Print out chain definition */
              if (valid == TRUE)
                fprintf(pymol_script_ptr,
                        "create %s, (resi %s and not het)\n",object,
                        residue_range);
              else
                fprintf(pymol_script_ptr,"create %s, not het\n",object);
            }

          /* Otherwise, for non-blank chain-id, create name chain_X */
          else
            {
              /* Define chain name */
              sprintf(object,"chain_%c",chain_ptr->chain_id);

              /* If this is a duplicate biomol chain, include
                 residue definition in its object name */
              if (biomol == TRUE && chain_ptr->symm_chain == TRUE &&
                  valid == TRUE)
                fprintf(pymol_script_ptr,
                        "create %s, (resi %s and chain %c and not het) \n",
                        object,residue_range,chain_ptr->chain_id);

              /* Otherwise write out the definition using chain name
                 only */
              else
                fprintf(pymol_script_ptr,
                        "create %s, (chain %c and not het) \n",
                        object,chain_ptr->chain_id);
            }

          /* Create a pymol object record */
          pobject_ptr
            = create_pobject_record(&first_pobject_ptr,&last_pobject_ptr,
                                    object,type,chain_ptr->chain_id);

          /* Get the chain */
          chain = chain_ptr->chain_id;
          if (biomol == TRUE && colour_as_orig == TRUE)
            chain = chain_ptr->link_chain_ptr->chain_id;

          /* Define the chain colour */
          get_chain_colname(chain,colour_name,colour_coil,colour_major,
                            colour_minor);
          if (type == DNA)
            fprintf(pymol_script_ptr,"color %s, %s\n",colour_name,object);
          else
            fprintf(pymol_script_ptr,"color %s, %s\n",colour_coil,object);

          /* If not Alpha Fold, colour according to chain and secondary
             structure */
          if (pdb_type != ALPHA_FOLD_PDB)
            {
              /* Define secondary structure colours according to which
                 secondary structure type is the "dominant" one */
              if (type != DNA && chain_ptr->dominant_sst == HELIX)
                {
                  fprintf(pymol_script_ptr,"color %s, (%s and ss h)\n",
                          colour_major,object);
                  fprintf(pymol_script_ptr,"color %s, (%s and ss s)\n",
                          colour_minor,object);
                }
              else if (type != DNA)
                {
                  fprintf(pymol_script_ptr,"color %s, (%s and ss h)\n",
                          colour_minor,object);
                  fprintf(pymol_script_ptr,"color %s, (%s and ss s)\n",
                          colour_major,object);
                }
            }

          /* Define display type for protein */
          if (type == PROTEIN)
            strcpy(display_type,"cartoon");

          /* Cartoon for CA */
          else if (type == CALPHA_ONLY)
            {
              fprintf(pymol_script_ptr,"set cartoon_trace, 1, %s\n",
                      object);
              strcpy(display_type,"cartoon");
            }

          /* Sticks for DNA/RNA */
          else if (type == DNA)
            strcpy(display_type,"cartoon");

          /* Determine whether chain to be shown transparent */
          if (biomol == TRUE && chain_ptr->symm_chain == TRUE &&
              transparent_copies == TRUE)
            {
              /* Show as partially transparent */
              fprintf(pymol_script_ptr,
                      "set %s_transparency, %3.1f, %s\n",
                      display_type,PYMOL_TRANSPARENCY[transparency],
                      object);
            }

          /* Print out the display type */
          fprintf(pymol_script_ptr,"show %s, %s\n",display_type,object);

          /* Save the display type */
          store_string(&(pobject_ptr->display_type),display_type);
        }

      /* For Alpha fold protein, colour according to AF score */
      if (pdb_type == ALPHA_FOLD_PDB && type == PROTEIN)
        {
          /* Initialise variables */
          last_range = -1;

          /* Get pointer to the first of this chain's residues */
          residue_ptr = chain_ptr->first_residue_ptr;
          nresid = chain_ptr->nresid;
          iresid = 0;

          /* Initialise start and end residues */
          first_res = last_res = 0;

          /* Loop through this chain's residues to colour according to
             AF score */
          while (residue_ptr != NULL && iresid < nresid)
            {
              /* If the range has changed, write out the previous range's
                 colour */
              if (residue_ptr->range != last_range && first_res > 0)
                {
                  /* Write out the last range */
                  fprintf(pymol_script_ptr,"colour %s, (resi %d-%d "
                          "and not het)\n",AF_COLOUR_NAME[last_range],
                          first_res,last_res);

                  /* Initialise first and last residues of new range */
                  first_res = last_res = 0;
                }

              /* Save the range */
              last_range = residue_ptr->range;

              /* Get the residue number */
              strcpy(number_string,residue_ptr->res_num);
              number_string[4] = '\0';
              number = atoi(number_string);

              /* If valid, then store */
              if (number > 0)
                {
                  /* If we don't have a first residue, save it */
                  if (first_res == 0)
                    first_res = number;

                  /* Save as the last residue */
                  last_res = number;
                }

              /* Get the next residue */
              residue_ptr = residue_ptr->next_residue_ptr;
              iresid++;
            }

          /* If we have a final range, then write out */
          if (first_res > 0)
            fprintf(pymol_script_ptr,"colour %s, (resi %d-%d "
                    "and not het)\n",AF_COLOUR_NAME[last_range],
                    first_res,last_res);
        }

      /* If this is a protein chain, check for any residues belonging
         to active sites */
      if (type == PROTEIN)
        {
          /* Save the name of the chain */
          strcpy(chain_name,object);

          /* Define name of site */
          if (chain_ptr->chain_id == ' ')
            strcpy(object,"sites_");
          else
            sprintf(object,"sites_%c",chain_ptr->chain_id);

          /* Initialiase residue counter */
          nsite = 0;

          /* Get pointer to the first of this chain's residues */
          residue_ptr = chain_ptr->first_residue_ptr;
          nresid = chain_ptr->nresid;
          iresid = 0;

          /* Loop through this chain's residues to find any active site
             residues */
          while (residue_ptr != NULL && iresid < nresid)
            {
              /* If this is an active site residue, add to definition */
              if (residue_ptr->site != &null)
                {
                  /* Get the residue number and check that valid */
                  strcpy(number_string,residue_ptr->res_num);
                  number_string[4] = '\0';
                  number = atoi(number_string);

                  /* If a PyMol-friendly residue number (ie not negative)
                     add to site definition */
                  if (number >= 0)
                    {
                      /* If first site residue, start definition */
                      if (nsite == 0)
                        {
                          fprintf(pymol_script_ptr,
                                  "# Active sites for chain [%c]\n",
                                  chain_ptr->chain_id);
                          fprintf(pymol_script_ptr,
                                  "create %s, "
                                  "(%s and resi %s and not het)\n",
                                  object,chain_name,residue_ptr->res_num);
                        }

                      /* Otherwise, add residue to current definition */
                      else
                        {
                          fprintf(pymol_script_ptr,
                                  "create %s, (%s or "
                                  "(%s and resi %s and not het))\n",
                                  object,object,chain_name,
                                  residue_ptr->res_num);
                        }

                      /* Increment count of site residues in this chain */
                      nsite++;
                    }
                }

              /* Get the next residue */
              residue_ptr = residue_ptr->next_residue_ptr;
              iresid++;
            }

          /* If any active site residues found for this chain, show
             side chains */
          if (nsite > 0)
            {
              /* Show sidechains of active site residues as sticks */
              fprintf(pymol_script_ptr,
                      "show sticks, (%s and not name n+c+o)\n",
                      object);

              /* Colour by atom type */
              fprintf(pymol_script_ptr,
                      "util.cbag (%s and not name n+ca+c+o)\n",
                      object);
              fprintf(pymol_script_ptr,
                      "color grey, %s and elem C and not name n+ca+c+o\n",
                      object);

              /* Create a pymol object record */
              pobject_ptr
                = create_pobject_record(&first_pobject_ptr,
                                        &last_pobject_ptr,
                                        object,ACTSITE,
                                        chain_ptr->chain_id);
              store_string(&(pobject_ptr->display_type),"sticks");


              /* Flag chain as having active sites */
              chain_ptr->got_sites = TRUE;
              pobject_ptr->got_sites = TRUE;
            }
        }

      /* Get the next chain */
      chain_ptr = chain_ptr->next_chain_ptr;
    }

  /* Return current pobject pointers */
  *fst_pobject_ptr = first_pobject_ptr;
  *lst_pobject_ptr = last_pobject_ptr;
}
/***********************************************************************

write_pymol_ligands  -  Create PyMol object for each ligand

***********************************************************************/

void write_pymol_ligands(struct ligand *first_ligand_ptr,int ntypes,
                         FILE *pymol_script_ptr,
                         struct pobject **fst_pobject_ptr,
                         struct pobject **lst_pobject_ptr,
                         int transparent_copies,int biomol)
{
  char chain, lobject[TOKEN_LEN + 1], object[TOKEN_LEN + 1];
  char name[3], res_num1[6], res_num2[6], residue_range[30], tmp[30];

  int iresid, nrange, nresid, rnum1, rnum2, transparency, type;
  int first, valid;

  struct chain *chain_ptr;
  struct ligand *ligand_ptr;
  struct pobject *pobject_ptr, *first_pobject_ptr, *last_pobject_ptr;
  struct residue *residue_ptr;

  /* Initialise variables */
  nrange = 0;

  /* Initialise pointers */
  first_pobject_ptr = *fst_pobject_ptr;
  last_pobject_ptr = *lst_pobject_ptr;

  /* Loop over all the ligand types */
  if (ntypes > 0)
    fprintf(pymol_script_ptr,"# Ligands\n");
  for (type = 1; type < ntypes + 1; type++)
    {
      /* Define a ligand object for this ligand type */
      sprintf(lobject,"ligand_%d",type);

      /* Get pointer to first ligand record */
      chain_ptr = NULL;
      ligand_ptr = first_ligand_ptr;
      first = TRUE;

      /* Loop through ligand records to find the one matching this type */
      while (ligand_ptr != NULL)
        {
          /* If not a metal, and this is the right type of ligand
             then add to definition */
          if (ligand_ptr->metal == FALSE && ligand_ptr->type == type)
            {
              /* Increment count of residue ranges */
              nrange++;

              /* Define the name of this ligand's residue range */
              sprintf(object,"lrange_%d",nrange);

              /* Initialise the first residue and last residues */
              res_num1[0] = '\0';
              res_num2[0] = '\0';

              /* Get this chain's first residue */
              residue_ptr = ligand_ptr->first_residue_ptr;
              nresid = ligand_ptr->nresid;
              iresid = 0;
              valid = FALSE;

              /* Loop through this chain's residues, to find start and
                 end residues */
              while (residue_ptr != NULL && iresid < nresid)
                {
                  /* If this is the first residue, save number */
                  if (res_num1[0] == '\0')
                    strcpy(res_num1,residue_ptr->res_num);

                  /* Save number as last residue */
                  strcpy(res_num2,residue_ptr->res_num);

                  /* Store chain id */
                  chain = residue_ptr->chain;
                  chain_ptr = residue_ptr->chain_ptr;

                  /* Get the next residue */
                  residue_ptr = residue_ptr->next_residue_ptr;
                  iresid++;
                }

              /* If have a residue range, write out definition */
              if (res_num1[0] != '\0' && res_num2[0] != '\0')
                {
                  /* Remove insertion codes as PyMol doesn't understand
                     them */
                  res_num1[4] = '\0';
                  res_num2[4] = '\0';

                  /* Define the residue range */
                  sprintf(residue_range,"%s-%s",res_num1,res_num2);

                  /* Check that residue numbers are in ascending order */
                  res_num1[4] = '\0';
                  res_num2[4] = '\0';
                  rnum1 = atoi(res_num1);
                  rnum2 = atoi(res_num2);
                  if (rnum1 > rnum2)
                    sprintf(residue_range,"%s-%s",res_num2,res_num1);

                  /* Remove blanks from the range */
                  valid = remove_residue_blanks(residue_range);

                  /* If chain not blank, then add chain */
                  if (chain != ' ')
                    {
                      sprintf(tmp,"%s and chain %c",residue_range,
                              chain);
                      strcpy(residue_range,tmp);
                    }

                  /* Write out this ligand definition */
                  if (valid == TRUE)
                    fprintf(pymol_script_ptr,
                            "create %s, (resi %s and het)\n",
                            object,residue_range);
                }

              /* If have a valid residue range then add to ligand
                 definition */
              if (valid == TRUE)
                {
                  /* If this is the first of this type, start definition */
                  if (first == TRUE)
                    {
                      fprintf(pymol_script_ptr,"create %s, %s\n",
                              lobject,object);

                      /* Unset flag */
                      first = FALSE;
                    }

                  /* Otherwise, add current range to existing
                     definition */
                  else
                    fprintf(pymol_script_ptr,"create %s, (%s or %s)\n",
                            lobject,lobject,object);
                }

              /* If not valid, then reset residue range count */
              else
                nrange--;
            }
          
          /* Get next ligand record */
          ligand_ptr = ligand_ptr->next_ligand_ptr;
        }

      /* If have a ligand defined, finish the definition and create
         a pobject record */
      if (first == FALSE)
        {
          /* Create a pymol object record */
          pobject_ptr
            = create_pobject_record(&first_pobject_ptr,&last_pobject_ptr,
                                    lobject,LIGAND,chain);
          store_string(&(pobject_ptr->display_type),"spheres");

          /* Store the ligand-type number */
          pobject_ptr->ligmet_type = type;

          /* Define how ligand is to be displayed */
          fprintf(pymol_script_ptr,"show spheres, %s\n",lobject);
          fprintf(pymol_script_ptr,"util.cbag %s\n",lobject);
          fprintf(pymol_script_ptr,"color black, (%s and elem C)\n",
                  lobject);

          /* If showing biomoloecule, determine whether this ligand to
             be shown transparent */
          if (biomol == TRUE && transparent_copies == TRUE)
            {
              /* Assume no transparency */
              transparency = NONE;

              /* Get corresponding chain */
              if (chain_ptr != NULL && chain_ptr->symm_chain == TRUE)
                transparency = PARTIAL;

              /* Show as transparent, if required */
              if (transparency != NONE)
                fprintf(pymol_script_ptr,"set sphere_transparency, %3.1f,"
                        " %s\n",PYMOL_TRANSPARENCY[transparency],lobject);
            }
        }
    }

  /* Return current pobject pointers */
  *fst_pobject_ptr = first_pobject_ptr;
  *lst_pobject_ptr = last_pobject_ptr;
}
/***********************************************************************

write_pymol_metals  -  Create PyMol object for each metal type

***********************************************************************/

void write_pymol_metals(char stored_metal[NMETALS][3],int *metal_count,
                        int nmetyp,struct ligand *first_ligand_ptr,
                        FILE *pymol_script_ptr,
                        struct pobject **fst_pobject_ptr,
                        struct pobject **lst_pobject_ptr,
                        int transparent_copies,int biomol)
{
  char definition[50], number_string[20], object[TOKEN_LEN + 1];
  char name[3], res_num1[6], res_num2[6], residue_range[13];
  char res_name[30];

  int colour, iresid, metal_type, nresid, nsite, number, type, valid;
  int ligmet_type, nmetals, transparency;
  int done, split;

  double rgb[3], rgb_hel[3], rgb_str[3];

  struct chain *chain_ptr;
  struct ligand *ligand_ptr;
  struct pobject *pobject_ptr, *first_pobject_ptr, *last_pobject_ptr;
  struct residue *residue_ptr;

  /* Initialise variables */
  ligmet_type = 0;

  /* Initialise pointers */
  first_pobject_ptr = *fst_pobject_ptr;
  last_pobject_ptr = *lst_pobject_ptr;

  /* Process all the metals, giving each its PDBsum colour */
  if (nmetyp > 0)
    fprintf(pymol_script_ptr,"# Metals\n");
  for (metal_type = 0; metal_type < nmetyp; metal_type++)
    {
      /* Get pointer to first ligand record */
      ligand_ptr = first_ligand_ptr;
      done = FALSE;
      nmetals = 0;

      /* Loop through ligand records to find the one matching
         this metal */
      while (ligand_ptr != NULL && done == FALSE)
        {
          /* If this is the right metal, create PyMol commands to
             display it as a sphere */
          if (ligand_ptr->metal == TRUE &&
              !strcmp(ligand_ptr->metal_name,stored_metal[metal_type]))
            {
              /* Store type */
              ligmet_type = ligand_ptr->type;

              /* Assume no transparency */
              transparency = NONE;

              /* Get the metal name, and replace leading blank if
                 necessary */
              strcpy(name,stored_metal[metal_type]);
              if (name[0] == ' ')
                name[0] = '_';

              /* Form name of metal object */
              strcpy(object,"metal");
              strcat(object,name);

              /* Define metal using element name */
              sprintf(definition,"(elem %s)",stored_metal[metal_type]);

              /* For biomolecule with transparency, need to show each
                 metal separately */
              if (biomol == TRUE && transparent_copies == TRUE)
                {
                  /* Create unique number to identify this metal */
                  nmetals++;
                  sprintf(number_string,"_%d",nmetals);

                  /* Add to object name */
                  strcat(object,number_string);

                  /* Get corresponding residue and chain */
                  residue_ptr = ligand_ptr->first_residue_ptr;
                  chain_ptr = residue_ptr->chain_ptr;

                  /* Define residue name */
                  if (chain_ptr->chain_id == ' ')
                    sprintf(res_name,"resi %s",residue_ptr->res_num);
                  else
                    sprintf(res_name,"chain %c and resi %s",
                            chain_ptr->chain_id,residue_ptr->res_num);

                  /* Include residue number and chain in metal's
                     definition */
                  sprintf(definition,"(%s and elem %s)",res_name,
                          stored_metal[metal_type]);

                  /* If chain is a transparent one, set flag */
                  if (chain_ptr->symm_chain == TRUE)
                    transparency = PARTIAL;
                }

              /* Create metal object */
              fprintf(pymol_script_ptr,"create %s, %s\n",object,
                      definition);

              /* Get the metal colour */
              colour = metal_type % NMETAL_COLOURS;

              /* Write out this metal's colour */
              fprintf(pymol_script_ptr,"colour %s, %s\n",
                      METAL_COLOUR[colour],object);

              /* Define its size */
              fprintf(pymol_script_ptr,"alter %s, vdw=%3.1f\n",
                      object,METAL_RADIUS);
              fprintf(pymol_script_ptr,"rebuild\n");

              /* Show as sphere */
              fprintf(pymol_script_ptr,"show spheres, %s\n",object);

              /* Show as transparent, if required */
              if (transparency != NONE)
                fprintf(pymol_script_ptr,"set sphere_transparency, %3.1f,"
                        " %s\n",PYMOL_TRANSPARENCY[transparency],object);

              /* Create a pymol object record */
              pobject_ptr
                = create_pobject_record(&first_pobject_ptr,
                                        &last_pobject_ptr,object,
                                        METAL,' ');
              store_string(&(pobject_ptr->display_type),"spheres");

              /* Store the ligand-type number */
              pobject_ptr->ligmet_type = ligmet_type;

              /* If not processing biomolecule, we're done for this
                 metal */
              if (biomol == FALSE || transparent_copies == FALSE)
                done = TRUE;
            }
          
          /* Get next ligand record */
          ligand_ptr = ligand_ptr->next_ligand_ptr;
        }
    }

  /* Return current pobject pointers */
  *fst_pobject_ptr = first_pobject_ptr;
  *lst_pobject_ptr = last_pobject_ptr;
}
/***********************************************************************

write_pymol_disulph  -  Create PyMol objects for cysteines involved in
                        disulphide bonds

***********************************************************************/

void write_pymol_disulph(struct chain *first_chain_ptr,
                         struct ssbond *first_ssbond_ptr,int nssbonds,
                         FILE *pymol_script_ptr,
                         struct pobject **fst_pobject_ptr,
                         struct pobject **lst_pobject_ptr)
{
  char object[TOKEN_LEN + 1], res_name[50];

  int ires, nres, type;

  struct chain *chain_ptr;
  struct pobject *pobject_ptr, *first_pobject_ptr, *last_pobject_ptr;
  struct residue *residue_ptr[2];
  struct ssbond *ssbond_ptr;

  /* Initialise pointers */
  first_pobject_ptr = *fst_pobject_ptr;
  last_pobject_ptr = *lst_pobject_ptr;

  /* Get pointer to the first chain */
  fprintf(pymol_script_ptr,"# Disulphide bonds\n");
  chain_ptr = first_chain_ptr;

  /* Loop through all the chains to define the various objects */
  while (chain_ptr != NULL)
    {
      /* Get chain type */
      type = chain_ptr->type;

      /* If protein chain then process */
      if (type == PROTEIN)
        {
          /* Define name for this chain's disulphides */
          if (chain_ptr->chain_id == ' ')
            strcpy(object,"disulph_");
          else
            sprintf(object,"disulph_%c",chain_ptr->chain_id);

          /* Get the first disulphide bond */
          ssbond_ptr = first_ssbond_ptr;
          nres = 0;

          /* Loop through the disulphide bonds to write the SSBOND records */
          while (ssbond_ptr != NULL)
            {
              /* Get the two residues */
              residue_ptr[0] = ssbond_ptr->residue1_ptr;
              residue_ptr[1] = ssbond_ptr->residue2_ptr;

              /* Loop over the two residues */
              for (ires = 0; ires < 2; ires++)
                {
                  /* If residue belongs to current chain, then
                     write out */
                  if (residue_ptr[ires]->chain == chain_ptr->chain_id)
                    {
                      /* Define residue name */
                      if (chain_ptr->chain_id == ' ')
                        sprintf(res_name,"resi %s",
                                residue_ptr[ires]->res_num);
                      else
                        sprintf(res_name,"chain %c and resi %s",
                                chain_ptr->chain_id,
                                residue_ptr[ires]->res_num);

                      /* If this is the first residue, then start
                         definition */
                      if (nres == 0)
                        fprintf(pymol_script_ptr,"create %s, %s\n",
                                object,res_name);

                      /* Otherwise, add to existing definition */
                      else
                        fprintf(pymol_script_ptr,"create %s, (%s or %s)\n",
                                object,object,res_name);

                      /* Increment count of residues */
                      nres++;
                    }
                }

              /* Get the next ssbond */
              ssbond_ptr = ssbond_ptr->next_ssbond_ptr;
            }

          /* Create a pymol object record */
          if (nres > 0)
            {
              /* Create a pymol object record */
              pobject_ptr
                = create_pobject_record(&first_pobject_ptr,
                                        &last_pobject_ptr,object,
                                        DISULPH,chain_ptr->chain_id);
              store_string(&(pobject_ptr->display_type),"sticks");

              /* Write out the disulphide definition for this chain */
              fprintf(pymol_script_ptr,"show sticks, %s\n",object);
              fprintf(pymol_script_ptr,"color yellow, (%s and name sg)\n",
                      object);
              fprintf(pymol_script_ptr,"color grey, (%s and name cb)\n",
                      object);
            }
        }

      /* Get the next chain */
      chain_ptr = chain_ptr->next_chain_ptr;
    }

  /* Return current pobject pointers */
  *fst_pobject_ptr = first_pobject_ptr;
  *lst_pobject_ptr = last_pobject_ptr;
}
/***********************************************************************

write_pymol_pml  -  Write command to PyMol script to define all objects
                    and their representation

***********************************************************************/

void write_pymol_pml(char *file,char *pdb_code,
                     struct chain *first_chain_ptr,
                     struct ligand *first_ligand_ptr,
                     struct ssbond *first_ssbond_ptr,int nssbonds,
                     char stored_metal[NMETALS][3],int *metal_count,
                     int nmetyp,int ntypes,FILE **pymol_script_ptr,
                     struct pobject **fst_pobject_ptr,
                     int colour_as_orig,int transparent_copies,
                     int biomol,int pdb_type)
{
  char file_name[FILENAME_LEN + 1], number_string[20];
  char object[TOKEN_LEN + 1];
  char colour_name[COLNAM_LEN], display_type[20];
  char colour_major[COLNAM_LEN + 6], colour_minor[COLNAM_LEN + 6];
  char colour_coil[COLNAM_LEN + 5];

  int colour, iresid, metal_type, nresid, nsite, number, type, valid;
  int split;

  double rgb[3], rgb_hel[3], rgb_str[3];

  struct pobject *pobject_ptr, *first_pobject_ptr, *last_pobject_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  split = TRUE;

  /* Initialise pointers */
  pobject_ptr = NULL;
  first_pobject_ptr = NULL;
  last_pobject_ptr = NULL;

  /* Form name of pymol script file, pymol.pml */
  strcpy(file_name,file);
  strcat(file_name,".pml");

  /* Open the output file */
  file_ptr = open_output_file(file_name,TRUE);

  /* Write out the header records */
  fprintf(file_ptr,"#\n");
  fprintf(file_ptr,"# PDBsum PyMol script for PDB entry %s\n",pdb_code);
  fprintf(file_ptr,"#\n");
  fprintf(file_ptr,"load %s.pdb\n",file);
  fprintf(file_ptr,"bg_color white\n");
  fprintf(file_ptr,"set depth_cue=0\n");
  fprintf(file_ptr,"set ray_trace_fog=0\n");
  fprintf(file_ptr,"set cartoon_ring_mode,1\n");
  fprintf(file_ptr,"hide all\n");
  fprintf(file_ptr,"set auto_color, 0\n");

  /* Define all the RGB colours */
  for (colour = 0; colour < NCOLOURS; colour++)
    fprintf(file_ptr,"set_color %s, [ %6.3f, %6.3f, %6.3f ]\n",
            COLOUR_LIST[colour],RGB_TRIPLET[colour][0],
            RGB_TRIPLET[colour][1],RGB_TRIPLET[colour][2]);

  /* Define all the major and minor secondary structure colours */
  for (colour = 0; colour < NCHAIN_COLOURS; colour++)
    {
      /* Get this chain's colour */
      get_chain_col(CHAIN_LIST[colour],rgb,rgb_hel,rgb_str,split);

      /* Define the coil colour */
      fprintf(file_ptr,"set_color %s_coil, [ %6.3f, %6.3f, %6.3f ]\n",
              CHAIN_COLOUR[colour],rgb[0],rgb[1],rgb[2]);

      /* Define the helix colour */
      fprintf(file_ptr,"set_color %s_major, [ %6.3f, %6.3f, %6.3f ]\n",
              CHAIN_COLOUR[colour],rgb_hel[0],rgb_hel[1],rgb_hel[2]);

      /* Define the strand colour */
      fprintf(file_ptr,"set_color %s_minor, [ %6.3f, %6.3f, %6.3f ]\n",
              CHAIN_COLOUR[colour],rgb_str[0],rgb_str[1],rgb_str[2]);
    }

  /* If this is an Alpha Fold structure, define the confidence colours */
  if (pdb_type == ALPHA_FOLD_PDB)
    {
      for (colour = 0; colour < NAF_RANGES; colour++)
        {
          /* Define the coil colour */
          fprintf(file_ptr,"set_color %s, [ %s ]\n",
                  AF_COLOUR_NAME[colour],AF_COLOUR_RGB[colour]);
        }
    }

  /* Write out the definitions for all the protein and DNA chains */
  write_pymol_chains(first_chain_ptr,file_ptr,&first_pobject_ptr,
                     &last_pobject_ptr,colour_as_orig,transparent_copies,
                     biomol,pdb_type);

  /* Write out the definitions for all the ligands */
  write_pymol_ligands(first_ligand_ptr,ntypes,file_ptr,
                      &first_pobject_ptr,&last_pobject_ptr,
                      transparent_copies,biomol);

  /* Write out the definitions for all the metals */
  write_pymol_metals(stored_metal,metal_count,nmetyp,first_ligand_ptr,
                     file_ptr,&first_pobject_ptr,&last_pobject_ptr,
                     transparent_copies,biomol);

  /* Add disulphide bonds */
  if (nssbonds > 0)
    write_pymol_disulph(first_chain_ptr,first_ssbond_ptr,nssbonds,
                        file_ptr,&first_pobject_ptr,&last_pobject_ptr);

  /* fprintf(file_ptr,"create disulph, (cys/ca+cb+sg) and byres "
     "(cys/sg and bound_to cys/sg)\n");
     fprintf(file_ptr,"show sticks, disulph\n");
     fprintf(file_ptr,"color yellow, (disulph and name sg)\n");
     fprintf(file_ptr,"color grey, (disulph and name cb)\n"); */

  /* Zoom in on final picture */
  fprintf(file_ptr,"zoom all, complete=1, buffer=2\n");

  /* Return first pobject pointer */
  *fst_pobject_ptr = first_pobject_ptr;

  /* Return pointer to script file */
  *pymol_script_ptr = file_ptr;
}
/***********************************************************************

convert_png_images  -  Convert all .png files generated by PyMol to
                       .gif or .jpg, as required

***********************************************************************/

void convert_png_images(struct pymol *last_pymol_ptr,
                        char *convert_exe)
{
  char command[FILENAME_LEN + 1];

  struct pymol *pymol_ptr;

  /* Get the first image to be converted */
  printf("Converting PyMol .png images ...\n");
  sprintf(command,"echo 'Converting PyMol .png images ...' > convert.log\n");
  system(command);
  pymol_ptr = last_pymol_ptr;

  /* Loop through all images to be converted */
  while (pymol_ptr != NULL)
    {
      /* Write convert command to log file */
      sprintf(command,"echo '%s %s %s' >> convert.log\n",convert_exe,
              pymol_ptr->png_name,pymol_ptr->gif_name);
      system(command);

      /* Form the parameters to run convert */
      sprintf(command,"%s %s %s  >> convert.log",convert_exe,
              pymol_ptr->png_name,pymol_ptr->gif_name);

      /* Run convert */
      system(command);

      /* Get the next image file */
      pymol_ptr = pymol_ptr->next_pymol_ptr;
    }
}
/***********************************************************************

write_biomol_dat  -  Write out the biomol.dat file

***********************************************************************/

void write_biomol_dat(int have_biomol,struct biomol *biomol_ptr,
                      char *output_dir)
{
  char file_name[FILENAME_LEN + 1];

  FILE *file_ptr;

  /* Initialise variables */

  /* Form output file name */
  strcpy(file_name,output_dir);
  strcat(file_name,"biomol.dat");

  /* Open the biomol.dat output file */
  file_ptr = open_output_file(file_name,TRUE);

  /* If we have biomolecule data, write out */
  if (have_biomol == TRUE)
    {
      /* Write flag indicating that we have biomolecule data */
      fprintf(file_ptr,"BIO_HAVE_BIOMOL TRUE\n");

      /* Write out number of biomolecule transforms given in the PDB
         file */
      if (biomol_ptr->nbiomolecule > 0)
        fprintf(file_ptr,"BIO_NBIOMAT %d\n",biomol_ptr->nbiomolecule);

      /* Indicate whether data come from the PQS database */
      if (biomol_ptr->in_pqs == TRUE)
        fprintf(file_ptr,"BIO_IN_PQS TRUE\n");
      else
        fprintf(file_ptr,"BIO_IN_PQS FALSE\n");

      /* Indicate whether data come from author information */
      if (biomol_ptr->have_author_info == TRUE)
        fprintf(file_ptr,"BIO_AUTHOR TRUE\n");
      else
        fprintf(file_ptr,"BIO_AUTHOR FALSE\n");

      /* If ASU equivalent to biological unit, set flag */
      if (biomol_ptr->is_biomolecule == TRUE)
        fprintf(file_ptr,"BIO_ASU_EQ_BIO TRUE\n");

      /* Show quaternary structure, if defined */
      if (biomol_ptr->quaternary_structure == NOT_GIVEN)
        fprintf(file_ptr,"BIO_QSTRUC NONE\n");
      else
        {
          /* Show quaternary structure */
          fprintf(file_ptr,"BIO_QSTRUC ");
          if (biomol_ptr->prefix[0] != '\0')
            fprintf(file_ptr,"%s",biomol_ptr->prefix);
          if (biomol_ptr->qstruc_desc[0] != '\0')
            fprintf(file_ptr,"%s",biomol_ptr->qstruc_desc);
          if (biomol_ptr->prefix[0] == '\0' &&
              biomol_ptr->qstruc_desc[0] == '\0')
            fprintf(file_ptr,"NOT GIVEN\n");
          else
            fprintf(file_ptr,"\n");

          /* Show number of units in quaternary structure */
          fprintf(file_ptr,"BIO_NQUATERN %3d\n",
                  biomol_ptr->quaternary_structure);

          /* Show whether have split asymmetric unit */
          if (biomol_ptr->split_asu == TRUE)
            fprintf(file_ptr,"BIO_SPLIT_ASU TRUE\n");
        }

      /* Repeat for PQS definitions */
      if (biomol_ptr->pqs_quaternary_structure == NOT_GIVEN)
        fprintf(file_ptr,"BIO_PQS_QSTRUC NONE\n");
      else
        {
          /* Show quaternary structure */
          fprintf(file_ptr,"BIO_PQS_QSTRUC ");
          if (biomol_ptr->pqs_prefix[0] != '\0')
            fprintf(file_ptr,"%s",biomol_ptr->pqs_prefix);
          if (biomol_ptr->pqs_qstruc_desc[0] != '\0')
            fprintf(file_ptr,"%s",biomol_ptr->pqs_qstruc_desc);
          if (biomol_ptr->pqs_prefix[0] == '\0' &&
              biomol_ptr->pqs_qstruc_desc[0] == '\0')
            fprintf(file_ptr,"NOT GIVEN\n");
          else
            fprintf(file_ptr,"\n");

          /* Show number of units in quaternary structure */
          fprintf(file_ptr,"BIO_PQS_NQUATERN %3d\n",
                  biomol_ptr->pqs_quaternary_structure);

          /* Show whether have split asymmetric unit */
          if (biomol_ptr->pqs_split_asu == TRUE)
            fprintf(file_ptr,"BIO_PQS_SPLIT_ASU TRUE\n");
        }

      /* If have generated a thumbnail image of the biomolecule,
         display it */
      if (biomol_ptr->have_biomol_image == TRUE)
        fprintf(file_ptr,"BIO_HAVE_IMAGE TRUE\n");

      /* Repeat for author definitions */
      if (biomol_ptr->author_quaternary_structure == NOT_GIVEN)
        fprintf(file_ptr,"BIO_AUTHOR_QSTRUC NONE\n");
      else
        {
          /* Show quaternary structure */
          fprintf(file_ptr,"BIO_AUTHOR_QSTRUC ");
          if (biomol_ptr->author_prefix[0] != '\0')
            fprintf(file_ptr,"%s",biomol_ptr->author_prefix);
          if (biomol_ptr->author_qstruc_desc[0] != '\0')
            fprintf(file_ptr,"%s",biomol_ptr->author_qstruc_desc);
          if (biomol_ptr->author_prefix[0] == '\0' &&
              biomol_ptr->author_qstruc_desc[0] == '\0')
            fprintf(file_ptr,"NOT GIVEN\n");
          else
            fprintf(file_ptr,"\n");

          /* Show number of units in quaternary structure */
          fprintf(file_ptr,"BIO_AUTHOR_NQUATERN %3d\n",
                  biomol_ptr->author_quaternary_structure);
        }
    }

  /* Otherwise, show that don't have biomol data */
  else
    fprintf(file_ptr,"BIO_HAVE_BIOMOL FALSE\n");

  /* Close the biomol.dat file */
  fclose(file_ptr);
}
/***********************************************************************

get_transforms  -  Compile the transformations to be applied to the
                   atomic coords

***********************************************************************/

int get_transforms(struct biomol *biomol_ptr,
                   float matrix[MAX_BIOMT][4][3])
{
  int i, j, ntrans;

  /* Initialise variables */
  ntrans = 0;
  for (i = 0; i < 4; i++)
    {
      for (j = 0; j < 3; j++)
        {
          if (i == j)
            matrix[ntrans][i][j] = 1.0;
          else
            matrix[ntrans][i][j] = 0.0;
        }
    }
  ntrans++;

  /* Return number of transformations stored */
  return(ntrans);
}
/***********************************************************************

find_chain  -  Find the given chain

***********************************************************************/

struct chain *find_chain(struct chain *first_chain_ptr,char chain,
                         int protein)
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
      if (chain == chain_ptr->chain_id)
        found_chain_ptr = chain_ptr;

      /* If chain has to be a protein chain, check that it is */
      if (protein == TRUE && (chain_ptr->type != PROTEIN &&
                              chain_ptr->type != CALPHA_ONLY &&
                              chain_ptr->type != DNA))
        found_chain_ptr = NULL;

      /* Get pointer to the next chain */
      chain_ptr = chain_ptr->next_chain_ptr;
    }

  /* Return the matching chain pointer */
  return(found_chain_ptr);
}
/***********************************************************************

find_residue  -  Find the given residue

***********************************************************************/

struct residue *find_residue(struct residue *first_residue_ptr,
                             char *res_name,char *res_num,char chain,
                             int offset,int top_offset)
{
  char adj_num[6], chain_id, number_string[20];

  int iresid;
  int adjust_num;

  struct residue *residue_ptr, *found_residue_ptr;

  /* Initialise variables */
  adjust_num = FALSE;
  if (offset == 0 && top_offset != 0)
    adjust_num = TRUE;
  found_residue_ptr = NULL;

  /* If have a residue offset, adjust the residue name to reflect
     this */
  if (offset != 0)
    {
      /* Get the numerical part of the residue number */
      strcpy(number_string,res_num);
      number_string[4] = '\0';
      iresid = atoi(number_string);

      /* Apply adjustment */
      iresid = iresid - offset;

      /* Generate new residue number */
      sprintf(adj_num,"%4d",iresid);
      strcpy(adj_num+4,res_num+4);
    }
  else
    strcpy(adj_num,res_num);

  /* Get pointer to the first residue */
  residue_ptr = first_residue_ptr;

  /* Loop over all the residues to find the one just read in */
  while (residue_ptr != NULL && found_residue_ptr == NULL)
    {
      /* Get this residue's chain */
      chain_id = residue_ptr->chain_ptr->chain_id;

      /* If have a match, store the pointer */
      if (chain == chain_id &&
          !strcmp(residue_ptr->res_num,adj_num) &&
          !strcmp(residue_ptr->res_name,res_name))
        {
          /* Store residue pointer */
          found_residue_ptr = residue_ptr;

          /* If need to adjust residue number, add offset to
             number to prevent duplication of residue numbers */
          if (adjust_num == TRUE)
            {
              strcpy(number_string,res_num);
              number_string[4] = '\0';
              iresid = atoi(number_string);

              /* Apply adjustment */
              iresid = iresid + top_offset;

              /* Generate new residue number */
              sprintf(adj_num,"%4d",iresid);
              strcpy(adj_num+4,res_num+4);
              strcpy(res_num,adj_num);
            }
        }

      /* Get pointer to the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
    }

  /* Return the matching residue pointer */
  return(found_residue_ptr);
}
/***********************************************************************

get_biomol_pdb  -  Copy across the parts of the PQS .mmol file that
                   are needed to generate PyMol image of biological
                   macromolecule

***********************************************************************/

int get_biomol_pdb(struct biomol *biomol_ptr,char *pdb_code,
                   struct chain *first_chain_ptr,
                   struct residue *first_residue_ptr,
                   char stored_metal[NMETALS][3],int *metal_count,
                   int *nmetyp,int *nmet,
                   struct chain **fst_biochain_ptr,
                   struct residue **fst_bioresidue_ptr,
                   struct atom **fst_bioatom_ptr,
                   float matrix[MAX_BIOMT][4][3],int *nmat,int *nats,
                   struct c_file_data file_name_ptr)
{
  char alt, atom_name[5], chain, chn[2], element[3], symm_chain;
  char res_name[4], res_num[6], residue_num[6], occup[4];
  char last_atom_name[5], last_chain, last_res_name[4], last_res_num[6];
  char input_line[LINELEN + 1], file_name[FILENAME_LEN + 1];
  char tmp_line[LINELEN + 1], number_string[20];

  char *string_ptr;

  int atom_no, atom_number, iresid, len, natoms, nbioresid, nchains;
  int nmatrix, nread, nresid, nstored, offset, rec_type, top_offset;
  int nkept, nmetal_types, nmetals, type;
  int first_A, found, have_atom, hydrogen, is_copy, ok, wanted;
  int have_map, protein;

  double bvalue, occupancy, inx, iny, inz, rmsd, x, y, z;

  struct atom *atom_ptr, *have_atom_ptr;
  struct atom *bioatom_ptr, *first_bioatom_ptr, *last_bioatom_ptr;
  struct chain *chain_ptr, *start_chain_ptr;
  struct chain *biochain_ptr, *first_biochain_ptr, *last_biochain_ptr;
  struct residue *residue_ptr, *have_residue_ptr;
  struct residue *bioresidue_ptr;
  struct residue *first_bioresidue_ptr, *last_bioresidue_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  atom_number = 0;
  chain_ptr = NULL;
  first_A = TRUE;
  have_map = FALSE;
  last_chain = '\0';
  last_atom_name[0] = '\0';
  last_res_name[0] = '\0';
  last_res_num[0] = '\0';
  nbioresid = 0;
  natoms = 0;
  nmatrix = 0;
  nread = 0;
  nstored = 0;
  offset = 0;
  protein = TRUE;
  rmsd = 0.0;
  top_offset = 0;
  residue_ptr = NULL;
  nkept = 0;
  nmetal_types = *nmetyp;
  nmetals = *nmet;
  ok = TRUE;

  /* Initialise pointers */
  first_bioatom_ptr = NULL;
  first_bioresidue_ptr = NULL;
  first_biochain_ptr = NULL;
  last_bioatom_ptr = NULL;
  last_bioresidue_ptr = NULL;
  last_biochain_ptr = NULL;

  /* Form name of appropriate PQS file */
  strcpy(file_name,file_name_ptr.other_dir[PQS_DATA_DIR]);
  strcat(file_name,DIR_SEP);
  strcat(file_name,pdb_code);
  if (biomol_ptr->pqs_split_asu == TRUE)
    strcat(file_name,"_1");
  strcat(file_name,".mmol");

  /* Open the input file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    {
      /* Show error message */
      printf("* Warning. Unable to open PQS file [%s]\n",
             file_name);
      return(FALSE);
    }

  /* Show message */
  printf("Getting PQS biomolecule: %s\n",file_name);

  /* Loop while reading in records from the file */
  while (fgets(input_line,LINELEN,file_ptr) != NULL)
    {
      /* Truncate the line at the last non-blank character */
      len = string_truncate(input_line,LINELEN);

      /* Assume line is wanted */
      wanted = TRUE;

      /* Check for REMARK 300 record */
      if (!strncmp(input_line,"REMARK 300 ",11))
        {
          /* Initialise chains */
          chain = '\0';
          symm_chain = '\0';
          chain_ptr = NULL;
          offset = 0;

          /* Check for original chain */
          if (!strncmp(input_line+11,"KEEPING",7))
            {
              /* Get retained chain */
              chain = input_line[25];
              if (chain == '-')
                chain = ' ';
              symm_chain = chain;
              nkept++;
            }

          /* Check for symmetry-related copy of chain */
          else if (!strncmp(input_line+11,"RENAMING",8))
            {
              /* Get retained chain */
              chain = input_line[20];
              if (chain == '-')
                chain = ' ';
              symm_chain = input_line[25];

              /* Check for any residue-number offset */
              if ((string_ptr = strstr(input_line,"ADDED TO RESNUM"))
                  != NULL)
                {
                  /* Retrieve the offset */
                  strncpy(number_string,string_ptr-7,7);
                  number_string[7] = '\0';
                  offset = atoi(number_string);
                }
            }

          /* If have a chain mapping, create record and link */
          if (symm_chain != '\0')
            {
              /* Create new chain record */
              biochain_ptr
                = create_chain_record(&first_biochain_ptr,
                                      &last_biochain_ptr,symm_chain);
              biochain_ptr->type = PROTEIN;

              /* Store the offset, in any and whether this is a
                 symmetry-related copy */
              biochain_ptr->offset = offset;

              /* Set flag that have initial chain map */
              have_map = TRUE;
              protein = TRUE;

              /* Find the corresponding chain read in from the PDB file */
              chain_ptr = find_chain(first_chain_ptr,chain,protein);

              /* If chain not found and we were looking for blank chain,
                 try again with search for chain A */
              if (chain_ptr == NULL && chain == ' ')
                chain_ptr = find_chain(first_chain_ptr,'A',protein);

              /* If found, then create a link record between the two */
              if (chain_ptr != NULL)
                {
                  /* Store the link */
                  biochain_ptr->link_chain_ptr = chain_ptr;

                  /* If chain not yet been 'copied' then this is the
                     same as the original PDB chain */
                  if (chain_ptr->used == FALSE)
                    {
                      /* Define bio chain as the original */
                      biochain_ptr->symm_chain = FALSE;

                      /* Flag chain as now having been used */
                      chain_ptr->used = TRUE;
                    }

                  /* Otherwise, treat this as a symmetry related copy
                     of the original chain */
                  else
                    biochain_ptr->symm_chain = TRUE;

                  /* Copy the chain type */
                  biochain_ptr->type = chain_ptr->type;
                }
            }

          /* Re-initialise chain pointers */
          chain_ptr = NULL;
          biochain_ptr = NULL;
        }

      /* Check for REMARK 350 record */
      if (!strncmp(input_line,"REMARK 350 ",11))
        {
          /* Check for oligomerization line */
          if (!strncmp(input_line+11,"AN OLIGOMER OF TYPE :",21))
            {
              strcpy(tmp_line,input_line+32);
              string_ptr = strstr(tmp_line," : ");
              string_ptr[0] = '\0';

              /* Extract quaternary structure type */
              if (string_ptr != NULL)
                biomol_ptr->pqs_quaternary_structure
                  = get_qstruc(tmp_line,biomol_ptr->pqs_prefix,
                               biomol_ptr->pqs_qstruc_desc);

              /* Re-compute number of symmetry operations this corresponds
                 to */
              if (biomol_ptr->nmolecules > 0)
                biomol_ptr->nsymop
                  = biomol_ptr->pqs_quaternary_structure
                  / biomol_ptr->nmolecules;
            }
        }

      /* Check for ATOM or HETATM record */
      else if (!strncmp(input_line,"ATOM  ",6) ||
               !strncmp(input_line,"HETATM",6))
        {
          /* Assume atom is wanted */
          wanted = TRUE;
          nread++;

          /* Get record type */
          if (!strncmp(input_line,"ATOM",4))
            rec_type = ATOM_REC;
          else
            rec_type = HETATM_REC;

          /* Get all the atom details */
          get_atom_details(input_line,&atom_no,atom_name,
                           &chain,res_name,res_num,&inx,&iny,&inz,
                           &bvalue,&occupancy,&alt,element);

          /* Rotate through 180 degrees about the x-axis */
          x = inx;
          y = - iny;
          z = - inz;

          /* Check whether this is a hydrogen atom */
          hydrogen = is_hydrogen_atom(atom_name);
          if (hydrogen == TRUE)
            wanted = FALSE;

          /* If this is a water, then don't want it */
          else if (!strcmp(res_name,"HOH") || !strcmp(res_name,"WAT"))
            wanted = FALSE;

          /* Otherwise, check whether this is a new chain */
          else if (chain != last_chain)
            {
              /* Initialise variables */
              found = FALSE;
              biochain_ptr = NULL;
              bioresidue_ptr = NULL;
              protein = TRUE;
              start_chain_ptr = first_biochain_ptr;

              /* Loop until current chain identified */
              while (found == FALSE && start_chain_ptr != NULL)
                {
                  /* Look for this chain */
                  if (have_map == TRUE)
                    biochain_ptr
                      = find_chain(start_chain_ptr,chain,protein);

                  /* Otherwise, need to create record for each new
                     chain */
                  else
                    {
                      /* If this is chain 'A' check if it is the
                         first */
                      if (chain == 'A')
                        {
                          /* If not first, then increment top offset */
                          if (first_A == FALSE)
                            top_offset = top_offset + 1000;
                          else
                            first_A = FALSE;
                        }

                      /* Create new chain record */
                      biochain_ptr
                        = create_chain_record(&first_biochain_ptr,
                                              &last_biochain_ptr,chain);
                      biochain_ptr->type = PROTEIN;

                      /* Find the corresponding chain read in from the
                         PDB file */
                      protein = TRUE;
                      if (first_A == TRUE)
                        chain_ptr
                          = find_chain(first_chain_ptr,chain,protein);
                      else
                        chain_ptr = NULL;

                      /* If found, then create a link record between
                         the two */
                      if (chain_ptr != NULL)
                        {
                          /* Store the link */
                          biochain_ptr->link_chain_ptr = chain_ptr;
                        }

                      /* Otherwise, need to identify which chain this is
                         a symmetry-related copy of */
                      else
                        {
                          /* @@@ */
                        }
                    }

                  /* If found, check whether chain has already
                     been processed and have a new version (ie
                     large multimers where we get multiple chain As,
                     each with a different residue offset) */
                  if (biochain_ptr != NULL)
                    {
                      /* If used, then restart search at next
                         record */
                      if (biochain_ptr->used == TRUE)
                        start_chain_ptr
                          = biochain_ptr->next_chain_ptr;

                      /* Otherwise, have the chain record we're
                         after */
                      else
                        {
                          /* Mark record as used */
                          biochain_ptr->used = TRUE;

                          /* Set flag that found */
                          found = TRUE;

                          /* If chain found, get its offset */
                          offset = biochain_ptr->offset;
                          if (offset != 0)
                            top_offset = offset;
                      
                          /* Set offset to current top value */
                          biochain_ptr->offset = top_offset;
                        }
                    }

                  /* If record not found, then need to skip
                     this entire chain as unable to handle it */
                  else
                    {
                      /* Set pointer so that loop ends */
                      start_chain_ptr = NULL;

                      /* Show warning message */
                      printf("* Warning. Unable to find record "
                             "for chain [%c] from line:\n",chain);
                      printf("*** %s\n",input_line);
                    }
                }

              /* Store the current chain */
              last_chain = chain;
            }

          /* If have an invalid chain, then don't want this atom */
          if (biochain_ptr == NULL ||
              biochain_ptr->link_chain_ptr == NULL)
            wanted = FALSE;

          /* If atom potentially still wanted, check residue type */
          if (wanted == TRUE)
            {
              /* If first residue or different from previous,
                 identify this residue */
              if (bioresidue_ptr == NULL ||
                  strcmp(res_name,last_res_name) ||
                  strcmp(res_num,last_res_num))
                {
                  /* Create a residue record */
                  bioresidue_ptr
                    = create_residue_record(&first_bioresidue_ptr,
                                            &last_bioresidue_ptr,
                                            biochain_ptr,res_name,
                                            res_num,
                                            biochain_ptr->chain_id,
                                            nbioresid);

                  /* Increment count of residues */
                  nbioresid++;

                  /* Get the corresponding chain in the original PDB
                     file */
                  chain_ptr = biochain_ptr->link_chain_ptr;

                  /* Find the corresponding residue in the PDB file */
                  strcpy(residue_num,res_num);
                  residue_ptr
                    = find_residue(first_residue_ptr,res_name,residue_num,
                                   chain_ptr->chain_id,offset,top_offset);

                  /* If not found, and we are searching for chain A, 
                     try blank chain */
                  if (residue_ptr == NULL && chain_ptr->chain_id == 'A')
                    residue_ptr
                    = find_residue(first_residue_ptr,res_name,residue_num,
                                   ' ',offset,top_offset);

                  /* Store residue number in case it was changed */
                  strcpy(bioresidue_ptr->res_num,residue_num);

                  /* Add link between the two residues */
                  bioresidue_ptr->link_residue_ptr = residue_ptr;

                  /* If not found, show warning */
                  if (residue_ptr == NULL)
                    {
                      printf("* Warning. Residue %s %s %c not found"
                             " in original\n",res_name,residue_num,
                             chain_ptr->chain_id);
                      wanted = FALSE;
                    }

                  /* Otherwise, copy across information from the
                     original residue record */
                  else
                    {
                      bioresidue_ptr->disulphide
                        = residue_ptr->disulphide;
                      if (residue_ptr->site[0] != '\0')
                        store_string(&(bioresidue_ptr->site),
                                     residue_ptr->site);
                      bioresidue_ptr->sec_struc
                        = residue_ptr->sec_struc;
                      bioresidue_ptr->ligand_ptr
                        = residue_ptr->ligand_ptr;
                    }

                  /* Save the current residue details */
                  strcpy(last_res_name,res_name);
                  strcpy(last_res_num,res_num);
                }

            }

          /* If have an invalid residue, then don't want this
             atom */
          if (bioresidue_ptr == NULL ||
              bioresidue_ptr->link_residue_ptr == NULL)
            wanted = FALSE;

          /* If atom still wanted, check whether to keep it */
          if (wanted == TRUE)
            {
              /* Get corresponding atom and residue from the original
                 PDB file */
              chain_ptr = biochain_ptr->link_chain_ptr;
              residue_ptr = bioresidue_ptr->link_residue_ptr;

              /* Get the chain type and */
              type = residue_ptr->chain_ptr->type;

              /* For protein atom, discard if not a mainchain atom */
              if (type == PROTEIN &&
                  strcmp(atom_name," N  ") &&
                  strcmp(atom_name," CA ") &&
                  strcmp(atom_name," C  ") &&
                  strcmp(atom_name," O  "))
                wanted = FALSE;

              /* If residue belongs to a disulphide bond, or is one of
                 the SITE residues, we want all its atoms */
              if (residue_ptr->disulphide == TRUE ||
                  residue_ptr->site[0] != '\0')
                wanted = TRUE;
            }

          /* If atom still wanted, then check for duplicates */
          if (wanted == TRUE)
            {
              /* Extract the occupancy */
              strncpy(occup,input_line+56,3);
              occup[3] = '\0';

              /* If this is an alternate-occupancy atom, check whether
                 we already have it */
              if (input_line[16] != ' ' || strncmp(occup,"1.0",3))
                {
                  /* Check for earlier copy of this atom */
                  have_atom
                    = check_for_duplicate(bioresidue_ptr,atom_name,
                                          bioresidue_ptr->res_name,
                                          bioresidue_ptr->res_num,
                                          chain);

                  /* If have this atom already, then don't want this
                     second copy */
                  if (have_atom == TRUE)
                    wanted = FALSE;
                }
            }

          /* If atom still wanted, then prepare to store */
          if (wanted == TRUE)
            {
              /* Increment atom number */
              atom_number++;
              if (atom_number > 99999)
                atom_number = 1;
              nstored++;

              /* Create an atom record */
              bioatom_ptr
                = create_atom_record(&first_bioatom_ptr,&last_bioatom_ptr,
                                     bioresidue_ptr,atom_number,
                                     atom_name,x,y,z,bvalue,occupancy,
                                     stored_metal,metal_count,
                                     &nmetal_types,&nmetals,rec_type);

              /* Identify the corresponding atom in the original PDB
                 file */
              atom_ptr = find_atom(residue_ptr,atom_name);

              /* Store link between the two */
              bioatom_ptr->link_atom_ptr = atom_ptr;

              /* If atoms are the same, compute rmsd between them */
              if (atom_ptr != NULL &&
                  biochain_ptr->symm_chain == FALSE)
                {
                  /* Get the atom coords */
                  inx = atom_ptr->original_coords[0];
                  iny = atom_ptr->original_coords[1];
                  inz = atom_ptr->original_coords[2];

                  /* Compute rmsd */
                  rmsd = rmsd
                    + (x - inx) * (x - inx)
                    + (y - iny) * (y - iny)
                    + (z - inz) * (z - inz);

                  /* Increment count of atoms compared */ 
                  natoms++;
                }
            }
        }
    }

  /* Close the file */
  fclose(file_ptr);

  /* Show stats */
  printf("Number of atoms read in from .mmol file:   %8d\n",nread);
  printf("Number of atoms stored:                    %8d\n",nstored);

  /* Compute the rmsd of all initial atoms */
  if (natoms > 0)
    rmsd = sqrt(rmsd / natoms);
  printf("\n");
  printf("Rmsd of %6d retained atoms:             %8.1f\n",natoms,
         rmsd);

  /* If rmsd too high, then abandon here */
  if (rmsd > 2.0)
    {
      printf("* Warning. Rmsd too high to rely on PQS file. Abandoning "
             "PDB entry %s\n",pdb_code);
      ok = FALSE;
    }

  /* Return numbers of atoms and matrices stored */
  *nats = nstored;
  *nmat = nmatrix;

  /* Return metal counts */
  *nmetyp = nmetal_types;
  *nmet = nmetals;

  /* Return pointers */
  *fst_bioatom_ptr = first_bioatom_ptr;
  *fst_bioresidue_ptr = first_bioresidue_ptr;
  *fst_biochain_ptr = first_biochain_ptr;

  /* Return whether all OK */
  return(ok);
}
/***********************************************************************

mirror_molecules  -  Match chain-structure in biomolecule to that in
                     original PDB file

***********************************************************************/

void mirror_molecules(struct chain *first_biochain_ptr,
                      struct residue *first_bioresidue_ptr,
                      struct chain *first_chain_ptr)
{
  int iresid, nresid;
  int found, symm_chain;

  struct chain *chain_ptr, *biochain_ptr, *last_biochain_ptr;
  struct chain *prev_chain_ptr;
  struct residue *residue_ptr, *bioresidue_ptr;

  /* Initialise last biomolecule chain */
  last_biochain_ptr = NULL;

  /* Get pointer to the first chain */
  biochain_ptr = first_biochain_ptr;

  /* Loop through all the biochains to get pointer to last one */
  while (biochain_ptr != NULL)
    {
      /* Save this record as last one */
      last_biochain_ptr = biochain_ptr;

      /* Get pointer to the next chain */
      biochain_ptr = biochain_ptr->next_chain_ptr;
    }

  /* Get pointer to the first of the original PDB chains */
  chain_ptr = first_chain_ptr;

  /* Loop through all the chains */
  while (chain_ptr != NULL)
    {
      /* Get pointer to the first of the biomolecule chains */
      biochain_ptr = first_biochain_ptr;
      found = FALSE;

      /* Loop to find the one corresponding to current PDB chain  */
      while (biochain_ptr != NULL)
        {
          /* If this is the correct copy, then adjust residue
             count, if necessary */
          if (biochain_ptr->link_chain_ptr == chain_ptr)
            {
              /* Check that no more residues in original */
              if (chain_ptr->nresid < biochain_ptr->nresid)
                biochain_ptr->nresid = chain_ptr->nresid;

              /* Set flag that match found */
              found = TRUE;
            }

          /* Get pointer to the next chain */
          biochain_ptr = biochain_ptr->next_chain_ptr;
        }

      /* If no matches found, then need to create appropriate number
         of mirror chains */
      if (found == FALSE)
        {
          /* Get this chain's first residue */
          biochain_ptr = NULL;
          residue_ptr = chain_ptr->first_residue_ptr;
          nresid = chain_ptr->nresid;
          iresid = 0;
          symm_chain = FALSE;
          prev_chain_ptr = NULL;

          /* Get pointer to the first bioresidue */
          bioresidue_ptr = first_bioresidue_ptr;

          /* Loop through all the residues to classify */
          while (bioresidue_ptr != NULL)
            {
              /* If this residue points to the original PDB residue,
                 then is a copy, so create new biochain to hold it */
              if (bioresidue_ptr->link_residue_ptr == residue_ptr)
                {
                  /* Create a new biomolecule chain */
                  biochain_ptr
                    = create_chain_record(&first_biochain_ptr,
                                          &last_biochain_ptr,
                                          bioresidue_ptr->chain);

                  /* Store the offset */
                  biochain_ptr->offset
                    = bioresidue_ptr->chain_ptr->offset;

                  /* Copy across other details */
                  biochain_ptr->type = chain_ptr->type;

                  /* Store pointer to first residue */
                  biochain_ptr->first_residue_ptr = bioresidue_ptr;
                  biochain_ptr->nresid++;

                  /* Save the chain it was previously pointing to */
                  prev_chain_ptr = bioresidue_ptr->chain_ptr;

                  /* Get residue to point to new chain */
                  bioresidue_ptr->chain_ptr = biochain_ptr;

                  /* Store flag indicating whether this is an original
                     or symmetry-related copy */
                  biochain_ptr->symm_chain = symm_chain;

                  /* Set flag to indicate that all subsequent copies
                     are symmetry-generated copies */
                  symm_chain = TRUE;
                }

              /* Otherwise, if residue points to same chain as
                 previous one did, add it to the new biochain */
              else if (bioresidue_ptr->chain_ptr == prev_chain_ptr &&
                       biochain_ptr->nresid < nresid)
                {
                  /* Get residue to point to new chain */
                  bioresidue_ptr->chain_ptr = biochain_ptr;
                  biochain_ptr->nresid++;
                }

              /* Otherwise, re-initialise previous chain */
              else
                prev_chain_ptr = NULL;


              /* Get the next residue */
              bioresidue_ptr = bioresidue_ptr->next_residue_ptr;
            }
        }

      /* Get pointer to the next chain */
      chain_ptr = chain_ptr->next_chain_ptr;
    }
}
/***********************************************************************

copy_disulphides  -  Copy the view coordinate over original coords

***********************************************************************/

int copy_disulphides(struct ssbond *first_ssbond_ptr,
                     struct ssbond **fst_biossbond_ptr,
                     struct residue *first_bioresidue_ptr)
{
  int loop, nmatch, nssbond, try;

  double distsq, x1, x2, y1, y2, z1, z2;

  struct atom *atom1_ptr, *atom2_ptr;
  struct ssbond *ssbond_ptr;
  struct ssbond *biossbond_ptr, *first_biossbond_ptr, *last_biossbond_ptr;
  struct residue *residue_ptr[2], *link_residue_ptr;
  struct residue *bioresidue_ptr, *match_bioresidue_ptr[2];
  struct residue *start_bioresidue_ptr;

  /* Initialise variables */
  nssbond = 0;

  /* Get pointer to the first bio residue */
  bioresidue_ptr = first_bioresidue_ptr;

  /* Loop through all the residues to re-initialise disulphide flag */
  while (bioresidue_ptr != NULL)
    {
      /* Set flag */
      bioresidue_ptr->disulphide = FALSE;
  
      /* Get pointer to the next residue */
      bioresidue_ptr = bioresidue_ptr->next_residue_ptr;
    }

  /* Initialise pointers */
  biossbond_ptr = NULL;
  first_biossbond_ptr = NULL;
  last_biossbond_ptr = NULL;

  /* Get pointer to the first ssbond */
  ssbond_ptr = first_ssbond_ptr;

  /* Loop through all the ssbonds */
  while (ssbond_ptr != NULL)
    {
      /* Get the two disulphides making up this residue */
      residue_ptr[0] = ssbond_ptr->residue1_ptr;
      residue_ptr[1] = ssbond_ptr->residue2_ptr;

      /* Initialise matches */
      match_bioresidue_ptr[0] = NULL;
      match_bioresidue_ptr[1] = NULL;
      start_bioresidue_ptr = NULL;
      try = 0;
      nmatch = 0;

      /* Get pointer to the first residue */
      bioresidue_ptr = first_bioresidue_ptr;

      /* Loop through all the bio residues to find those equivalent
         to either of these two */
      while (bioresidue_ptr != NULL)
        {
          /* If residue is a Cys and is not already marked as belonging
             to a disulphide bond, check it out */
          if (!strcmp(bioresidue_ptr->res_name,"CYS") &&
              bioresidue_ptr->disulphide == FALSE)
            {
              /* Get the residue this links to */
              link_residue_ptr = bioresidue_ptr->link_residue_ptr;

              /* Loop to check for match to either */
              for (loop = 0; loop < 2; loop++)
                {
                  /* If have match, then store */
                  if (link_residue_ptr == residue_ptr[loop] &&
                      match_bioresidue_ptr[loop] == NULL)
                    {
                      /* Store this residue */
                      match_bioresidue_ptr[loop] = bioresidue_ptr;

                      /* Increment match count */
                      nmatch++;

                      /* If have both residues, check that distance
                         between them is a sensible one for a disulphide
                         bond */
                      if (nmatch == 2)
                        {
                          /* Identify the two sulphur atoms involved */
                          atom1_ptr
                            = get_sulphur(match_bioresidue_ptr[0]);
                          atom2_ptr
                            = get_sulphur(match_bioresidue_ptr[1]);

                          /* If both atoms present, calculate the
                             distance between them */
                          if (atom1_ptr != NULL && atom2_ptr != NULL)
                            {
                              /* Get coords of the two atoms */
                              x1 = atom1_ptr->coords[0];
                              y1 = atom1_ptr->coords[1];
                              z1 = atom1_ptr->coords[2];
                              x2 = atom2_ptr->coords[0];
                              y2 = atom2_ptr->coords[1];
                              z2 = atom2_ptr->coords[2];

                              /* Calculate the squared distance between
                                 them */
                              distsq
                                = (x1 - x2) * (x1 - x2)
                                + (y1 - y2) * (y1 - y2)
                                + (z1 - z2) * (z1 - z2);
                            }

                          /* Otherwise, set distance too high */
                          else
                            distsq = 1.0 + MAX_SSDIST2;

                          /* If distance too long, then reject the
                             second residue */
                          if (distsq > MAX_SSDIST2)
                            {
                              /* Blank out current residue */
                              match_bioresidue_ptr[loop] = NULL;
                              nmatch = 1;
                            }
                        }

                      /* If this is the first residue matched,
                         save the current point in the bioresidue
                         list so that we can return to it later */
                      else if (nmatch == 1)
                        start_bioresidue_ptr = bioresidue_ptr;
                    }
                }

              /* If have links to both residues, then create a
                 disulphide bond record */
              if (nmatch == 2)
                {
                  /* Identify the two sulphur atoms involved */
                  atom1_ptr = get_sulphur(match_bioresidue_ptr[0]);
                  atom2_ptr = get_sulphur(match_bioresidue_ptr[1]);

                  /* Create ssbond record */
                  biossbond_ptr
                    = create_ssbond_record(&first_biossbond_ptr,
                                           &last_biossbond_ptr,
                                           match_bioresidue_ptr[0],
                                           match_bioresidue_ptr[1],
                                           atom1_ptr,atom2_ptr);

                  /* Increment count of disulphide bonds stored */
                  nssbond++;

                  /* Reset position of residue scan to first disulphide
                     residue of the current pair */
                  if (start_bioresidue_ptr != NULL)
                    bioresidue_ptr = start_bioresidue_ptr;

                  /* Re-initialise */
                  match_bioresidue_ptr[0] = NULL;
                  match_bioresidue_ptr[1] = NULL;
                  start_bioresidue_ptr = NULL;
                  nmatch = 0;
                }
            }

          /* If about to hit end with only one disulphide residue
             found, then reset to look for more */
          if (bioresidue_ptr->next_residue_ptr == NULL &&
              start_bioresidue_ptr != NULL && try < 100)
            {
              bioresidue_ptr = start_bioresidue_ptr;
              try++;
            }

          /* Get pointer to the next residue */
          bioresidue_ptr = bioresidue_ptr->next_residue_ptr;
        }

      /* Get pointer to the next ssbond */
      ssbond_ptr = ssbond_ptr->next_ssbond_ptr;
    }

  /* Return pointer to first disulphide */
  *fst_biossbond_ptr = first_biossbond_ptr;

  /* Return number of records created */
  return(nssbond);
}
/***********************************************************************

copy_coords  -  Copy the view coordinate over original coords

***********************************************************************/

void copy_coords(struct atom *first_atom_ptr,int natoms)
{
  int iatom;

  struct atom *atom_ptr;

  /* Get pointer to the first atom */
  atom_ptr = first_atom_ptr;
  iatom = 0;

  /* Loop through all the atoms */
  while (atom_ptr != NULL && iatom < natoms)
    {
      /* Copy the coordinates */
      atom_ptr->coords[0] = atom_ptr->view_coords[0];
      atom_ptr->coords[1] = atom_ptr->view_coords[1];
      atom_ptr->coords[2] = atom_ptr->view_coords[2];

      /* Get pointer to the next atom */
      atom_ptr = atom_ptr->next_atom_ptr;
      iatom++;
    }
}
/***********************************************************************

wait_for_file  -  Keep checking that file has been fully written before
                  calling next stage

***********************************************************************/

void wait_for_file(char *file_name,int dos)
{
  char command[FILENAME_LEN + 1], input_line[LINELEN + 1];

  int len, line, ntries;
  int file_done, have_file;

  FILE *file_ptr;

  /* Initialise variables */
  file_done = FALSE;
  ntries = 0;

  /* Loop until file has been written */
  while (file_done == FALSE && ntries < 100)
    {
      /* Initialise line count */
      line = 0;
          
      /* Open the file */
      file_ptr = fopen(file_name,"r");

      /* If have file, then read through to the end */
      if (file_ptr != NULL)
        {
          /* Loop while reading in records from the file */
          while (fgets(input_line,LINELEN,file_ptr) != NULL)
            {
              /* Truncate the line at the last non-blank character */
              len = string_truncate(input_line,LINELEN);

              /* Increment line counter */
              line++;

              /* See if we have hit the end record */
              if (!strcmp(input_line,"## END"))
                file_done = TRUE;
            }

          /* Close the file */
          fclose(file_ptr);
        }
      
      /* Increment count of tries */
      ntries++;

/* debug */
      printf("Try %d.  lines %d  file done %s\n",ntries,line,T_F[file_done]);
      fflush(stdout);
/* debug */

      /* Pause before next try */
      sprintf(command,"%s",SLEEP_COMMAND);
      call_system(command,dos);
    }
/* debug */
  printf("Found after %d tries  file done %s\n",ntries,T_F[file_done]);
  fflush(stdout);
/* debug */
  
}

/***********************************************************************

                             M   A   I   N

***********************************************************************/

int main(int argc,char *argv[])
{
  char chains_dat[FILENAME_LEN + 1];
  char interfaces_dat[FILENAME_LEN + 1];
  char ligands_dat[FILENAME_LEN + 1];
  char ligmet_run[FILENAME_LEN + 1];
  char newsec_out[FILENAME_LEN];
  char output_dir[FILENAME_LEN + 1];
  char pdbsum_dir[FILENAME_LEN + 1];
  char sst_dat[FILENAME_LEN];
  char stored_metal[NMETALS][3];

  char *fname;

  int natoms, nchains, nligands, nligmet, nmetals, nresid, nsite_res;
  int ninterfaces, nssbonds, nmetypes, nshift, ntypes, nunique, nviews;
  int nbioatoms, nbiossbonds, nbioligands, nbiometals;
  int metal_count[NMETALS], nmatrix, nmetal_types;
  int cando_nucplot, have_biomol, have_dna, mcsg, ok, transparent;
  int colour_as_orig, transparent_copies, os_type;
  int dir_type;
  int dos;

  float matrix[MAX_BIOMT][4][3];

  double cofg[3], coord_min[3], coord_max[3], grid_sep;
  double transformation[3][3];

  struct biomol biomol_ptr;
  struct data data_item;
  struct c_file_data file_name_ptr;
  struct parameters params;

  struct atom *first_atom_ptr;
  struct atom *first_bioatom_ptr;
  struct chain *first_chain_ptr, *last_chain_ptr;
  struct chain *first_biochain_ptr;
  struct hbadd *first_hbadd_ptr;
  struct hetdesc *first_hetdesc_ptr;
  struct interface *first_interface_ptr;
  struct ligand *first_ligand_ptr, *first_bioligand_ptr;
  struct pobject *first_pobject_ptr;
  struct pymol *last_pymol_ptr;
  struct residue *first_residue_ptr;
  struct residue *first_bioresidue_ptr;
  struct ssbond *first_ssbond_ptr, *first_biossbond_ptr;

  FILE *ligmet_file_ptr, *pymol_script_ptr;

  /* Initialise variables */
  colour_as_orig = FALSE;
  transparent_copies = FALSE;
  first_ligand_ptr = NULL;
  output_dir[0] = '\0';
  ligmet_file_ptr = NULL;
  pymol_script_ptr = NULL;
  last_pymol_ptr = NULL;
  nmatrix = 0;

  /* Form output file names */
  strcpy(chains_dat,output_dir);
  strcpy(interfaces_dat,output_dir);
  strcpy(ligands_dat,output_dir);
  strcpy(ligmet_run,output_dir);
  strcpy(newsec_out,output_dir);
  strcpy(sst_dat,output_dir);

  strcat(chains_dat,"chains.dat");
  strcat(interfaces_dat,"interfaces.dat");
  strcat(ligands_dat,"ligands.dat");
  strcat(ligmet_run,"ligmet.run");
  strcat(newsec_out,"newsec.out");
  strcat(sst_dat,"sst.dat");
  nviews = 1;
  ntypes = 0;
  nmetypes = 0;
  ninterfaces = 0;

  /* Interpret the command-line parameters */
  get_command_arguments(argv,argc - 1,&params);

  /* Get required paths and file names */
  c_get_paths(&file_name_ptr,params.run_64);

  /* Get whether running under DOS */
  dos = file_name_ptr.dos;
  
  /* Set the commands depending on whether DOS or non-DOS */
  os_type = NOT_DOS;
  if (file_name_ptr.dos == TRUE)
    os_type = DOS;

  /* Set the commands accordingly */
  strcpy(CP_COMMAND,CP_COPY[os_type]);
  strcpy(MV_COMMAND,MV_REN[os_type]);
  strcpy(RM_COMMAND,RM_DEL[os_type]);
  strcpy(SLEEP_COMMAND,SLEEP[os_type]);

  /* Read in the structure from the file.new file */
  natoms = read_pdb(params.file_name,&first_chain_ptr,
                    &last_chain_ptr,&nchains,
                    &first_atom_ptr,&first_residue_ptr,
                    &nresid,&nsite_res,stored_metal,metal_count,
                    &nmetal_types,&nmetals,&data_item,
                    &first_hetdesc_ptr,&mcsg,&biomol_ptr,params);

  /* If structure excessively large, generate thumbnails using Raster3D
     rather than PyMol */
  if (natoms > MAX_PYMOL_ATOMS)
    params.gen_pymol = FALSE;

  /* Classify the residues according to type */
  classify_residues(first_residue_ptr,first_hetdesc_ptr,&last_chain_ptr);

  /* Shift all chains containing single metal ions to the end of
     the chains linked list */
  nshift = shift_metals(&first_chain_ptr,&last_chain_ptr);

  /* Adjust the residues so that their order corresponds to the current
     chain order */
  if (nshift > 0)
    nresid = adjust_residues(&first_residue_ptr,first_chain_ptr);

  /* Form the name of the PDBsum directory */
  dir_type = find_pdbsum_dir(params.pdb_code,file_name_ptr.pdbsum_template,
                             params.pdb_type,pdbsum_dir);
  if (dir_type == NOT_FOUND)
    {
      printf("*** ERROR. Unable to locate PDBsum directory for %s\n",
             params.pdb_code);
      exit(-1);
    }

  /* Read in the secondary structure definitions from the file.rin
     file */
  read_rin_file(params.pdb_code,pdbsum_dir,first_residue_ptr,
                params.pdb_type);

  /* Get all the disulphides in the structure */
  nssbonds = get_disulphides(params.pdb_code,pdbsum_dir,
                             &first_ssbond_ptr,first_residue_ptr,
                             params.pdb_type);

  /* Initial classification of the "chains" */
  nchains = classify_chains(first_chain_ptr);

  /* If have more than 1 chain, check for consecutive chains with the
     same chain-ID which can be combine together into a single chain
     (currently split into separate chains by spurious TER records in
     the PDB file identifying breaks in the chain) */
  if (nchains > 1)
    join_chains(first_chain_ptr,params);

  /* Loop through all the chains to reclassify all very short chains
     as ligands */
  reclassify_chains(first_chain_ptr,params);

  /* Split the ligand chains into physically distinct molecules */
  split_ligands(first_chain_ptr);

  /* Identify which single-residue ligands are metals and which are
     actually residue attachments */
  single_residue_ligands(first_chain_ptr);

  /* Calculate the coordinate bounds of each chain */
  nviews
    = calc_chain_limits(first_chain_ptr,coord_min,coord_max,&have_dna,
                        &cando_nucplot);

  /* Write out all the molecule definitions */
  write_defs(newsec_out,first_chain_ptr,&data_item);

  /* Check whether we have a definition for this structure's biological
     unit */
  if (params.no_gifs == FALSE)
    have_biomol = check_biomol(&biomol_ptr,data_item,file_name_ptr,params);
  else
    have_biomol = FALSE;

  /* Read in the Het Groups defined in hbadd.sum */
  read_hbadd_sum(params.pdb_code,pdbsum_dir,first_chain_ptr,
                 &first_hbadd_ptr,params.pdb_type);

  /* For PDBsum1, open the ligmet.run file */
  if (params.pdb_type == PDBSUM1)
    ligmet_file_ptr = open_output_file(ligmet_run,TRUE);

  /* If hbadd.sum file found, write out the .pdb and .rcm files for
     each Het Group */
  write_hbadd_files(first_hbadd_ptr,output_dir,
                    file_name_ptr.other_dir[HETGIF_REAL_DIR],
                    params.pdb_code,file_name_ptr.exe_name[RASTER3D],
                    file_name_ptr.exe_name[CONVERT_EXE],
                    file_name_ptr.exe_name[PPMQUANT_EXE],
                    file_name_ptr.exe_name[PPMTOGIF_EXE],
                    params.molecule_align,ligmet_file_ptr,params);

  /* For protein chains, extract the amino acid sequences */
  extract_sequences(first_chain_ptr);

  /* Compare the sequences */
  nunique = compare_sequences(first_chain_ptr);

  /* Write out the chain details to the chains.dat file */
  write_chain_equivs(first_chain_ptr,chains_dat,output_dir,params);

  /* Create all the ligand records for each of the ligand chains */
  nligmet = form_ligands(first_chain_ptr,&first_ligand_ptr,&nligands,
                         &nmetals,FALSE);

  /* Loop through all the ligands to pick up all neighbouring residues */
  get_ligand_neighbours(first_ligand_ptr,first_chain_ptr);

  /* If have any ligands, process them to determine which are unique,
     and write out the data to the ligands.dat file */
  if (nligmet > 0)
    {
      /* Loop through all the ligands to determine which are unique */
      ntypes
        = identify_unique_ligands(first_ligand_ptr,first_chain_ptr,
                                  &nmetypes);

      /* Sort ligands by type, uniqueness and duplicate number */
      if (nligmet > 1)
        sort_ligands(&first_ligand_ptr,nligmet);

      /* Write out the .run files for the ligands and metals */
      write_ligmet_run(first_ligand_ptr,ligmet_file_ptr);
    }

  /* If this is an Alpha Fold structure, write out the AF_scores.dat file
     giving the confidence scores for each residue */
  if (params.pdb_type == ALPHA_FOLD_PDB)
    write_af_scores(first_residue_ptr,pdbsum_dir);

  /* If generating all images, run all routines */
  if (params.biomol_only == FALSE)
    {
      /* If have an uploaded PDB file, run HBPLUS on the first chain */
      if (params.uploaded_file == TRUE)
        run_hbplus(file_name_ptr.exe_name[HBPLUS_EXE],output_dir);

      /* Read in the H-bonds and non-bonded contacts from the
         .hhb and .nnb files and write out a grow.out file summarising the
         interactions */
      create_grow_out(params.pdb_code,params.file_name,first_residue_ptr,
                      nresid,&first_interface_ptr,&ninterfaces,
                      first_ssbond_ptr,first_chain_ptr,cando_nucplot,
                      output_dir,params);

      /* Write out all the interface data */
      if (params.calc_interfaces == TRUE && ninterfaces > 0)
        save_interface_data(interfaces_dat,file_name_ptr,
                            &first_interface_ptr,ninterfaces,
                            first_chain_ptr,first_residue_ptr,output_dir,
                            params.pdb_type);

      /* If generating images, create all files needed */
      if (params.no_gifs == FALSE)
        {
          /* If generating PyMol images, prepare a fresh version of the
             PDB file incorporating all the secondary structure definitions */
          if (params.gen_pymol == TRUE)
            {
              /* Write out pymol.pdb */
              write_pymol_pdb(params.pdb_code,first_chain_ptr,
                              first_ssbond_ptr,nssbonds,"pymol.pdb",
                              FALSE);

              /* Write out the secondary structure counts to the sst.dat
                 file */
              write_sst_dat(first_chain_ptr,sst_dat);

              /* Create a pymol.pml script and write out the header records
                 defining all the objects */
              write_pymol_pml("pymol",params.pdb_code,first_chain_ptr,
                              first_ligand_ptr,first_ssbond_ptr,nssbonds,
                              stored_metal,metal_count,nmetal_types,ntypes,
                              &pymol_script_ptr,&first_pobject_ptr,
                              colour_as_orig,transparent_copies,FALSE,
                              params.pdb_type);
            }

          /* Generate the required number of Raster3D views of the structure
             as a whole */
          transparent
            = render_structure(first_chain_ptr,first_atom_ptr,natoms,nviews,
                               first_ssbond_ptr,nssbonds,nsite_res,
                               file_name_ptr.exe_name[RASTER3D],
                               file_name_ptr.exe_name[CONVERT_EXE],
                               file_name_ptr.exe_name[PPMQUANT_EXE],
                               file_name_ptr.exe_name[PPMTOGIF_EXE],mcsg,
                               &grid_sep,ntypes,nmetypes,pymol_script_ptr,
                               &last_pymol_ptr,colour_as_orig,
                               transparent_copies,FALSE,params);

          /* Render only the DNA chains */
          if (have_dna == TRUE)
            transparent
              = render_all_dna(first_chain_ptr,first_ligand_ptr,
                               first_atom_ptr,natoms,
                               file_name_ptr.exe_name[RASTER3D],
                               file_name_ptr.exe_name[CONVERT_EXE],
                               file_name_ptr.exe_name[PPMQUANT_EXE],
                               file_name_ptr.exe_name[PPMTOGIF_EXE],
                               pymol_script_ptr,first_pobject_ptr,
                               &last_pymol_ptr,transparent,params);
  
          /* Render each protein chain separately */
          transparent
            = render_indiv_chains(first_chain_ptr,first_ligand_ptr,
                                  first_atom_ptr,natoms,
                                  first_ssbond_ptr,
                                  file_name_ptr.exe_name[RASTER3D],
                                  file_name_ptr.exe_name[CONVERT_EXE],
                                  file_name_ptr.exe_name[PPMQUANT_EXE],
                                  file_name_ptr.exe_name[PPMTOGIF_EXE],
                                  pymol_script_ptr,
                                  first_pobject_ptr,&last_pymol_ptr,
                                  transparent,params);
  
          /* If have any ligands, render the unique ones */
          if (nligmet > 0)
            {
              /* Generate a Raster3D image for each ligand type */
              render_ligands(first_chain_ptr,first_ligand_ptr,
                             file_name_ptr.exe_name[RASTER3D],
                             file_name_ptr.exe_name[CONVERT_EXE],
                             file_name_ptr.exe_name[PPMQUANT_EXE],
                             file_name_ptr.exe_name[PPMTOGIF_EXE],
                             output_dir,pymol_script_ptr,
                             first_pobject_ptr,&last_pymol_ptr,
                             transparent,GIFS,file_name_ptr,params);
            }

          /* If have any protein-protein or protein-DNA interfaces, render
             the molecules involved */
          if (params.calc_interfaces == TRUE && ninterfaces > 0 &&
              transparent == TRUE)
            {
              /* Generate a Raster3D image for each unique interface */
              render_interfaces(first_interface_ptr,first_chain_ptr,
                                first_ligand_ptr,first_atom_ptr,natoms,
                                file_name_ptr.exe_name[RASTER3D],
                                file_name_ptr.exe_name[CONVERT_EXE],
                                file_name_ptr.exe_name[PPMQUANT_EXE],
                                file_name_ptr.exe_name[PPMTOGIF_EXE],
                                pymol_script_ptr,first_pobject_ptr,
                                &last_pymol_ptr,params);
            }

          /* If using PyMol, run script to generate all the required
             images */
          if (params.gen_pymol == TRUE)
            {
              /* Close theWrite end record to PyMol script file */
              fprintf(pymol_script_ptr,"## END\n");
              
              /* Close the PyMol script file */
              fclose(pymol_script_ptr);

              /* Check that file has been completely written */
              wait_for_file("pymol.pml",dos);

/* debug */
  printf("RUNNING PyMOL ...\n");
  fflush(stdout);
/* debug */

              /* Run PyMol to generate all the ray-traced .png files */
              run_pymol("pymol",file_name_ptr.exe_name[PYMOL_EXE],dos);
/* debug */
  printf("CONVERTING images\n");
  fflush(stdout);
/* debug */

              /* Convert all the generated files to .gif or .jpg,
                 as required */
              if (params.pdb_type != PDBSUM1)
                convert_png_images(last_pymol_ptr,
                                   file_name_ptr.exe_name[CONVERT_EXE]);
            }
        }

      /* If don't need GIFs, run render ligand routine so that all ligand
         coords are written out for use by LIGPLOT */
      else
        {
          /* Determine whether have transparent images */
          transparent = FALSE;
          if (nunique > 1 || (nunique > 0 && nligmet > 0))
            transparent = TRUE;
          if (params.gen_pymol == FALSE)
            transparent = FALSE;

          /* Run render ligands routine */
          render_ligands(first_chain_ptr,first_ligand_ptr,
                         file_name_ptr.exe_name[RASTER3D],
                         file_name_ptr.exe_name[CONVERT_EXE],
                         file_name_ptr.exe_name[PPMQUANT_EXE],
                         file_name_ptr.exe_name[PPMTOGIF_EXE],output_dir,
                         pymol_script_ptr,first_pobject_ptr,
                         &last_pymol_ptr,transparent,NO_GIFS,
                         file_name_ptr,params);
        }

      /* Write out the ligands.dat file detailing all the ligands and their
         duplicates */
      if (nligmet > 0 || transparent == TRUE)
        write_dat_file(ligands_dat,first_ligand_ptr,ntypes,nmetypes,
                       first_hetdesc_ptr,stored_metal,metal_count,
                       nmetal_types,transparent);

      /* Write out the secsurf.out file */
      write_secsurf_out(first_chain_ptr,first_ligand_ptr,nresid,natoms,
                        stored_metal,metal_count,nmetal_types,data_item,
                        output_dir);

      /* Write out the surfp.par file, if reasonable to do so */
      if ((data_item.nprotein + data_item.ndna > 0) &&
          (data_item.nprotein + data_item.ndna < MAX_SURF_CHAINS) &&
          data_item.natom_records < MAX_SURF_ATOMS &&
          data_item.nca_only == 0)
        write_surfp_par(first_chain_ptr,output_dir,params.file_name,
                        grid_sep);

      /* Write out the rasmol.dfn file */
      write_rasmol_dfn(first_chain_ptr,first_ligand_ptr,nresid,natoms,
                       stored_metal,metal_count,nmetal_types,data_item,
                       output_dir);
    }

  /* If don't already have the biological unit in the PDB file, generate
     PyMol image of the biological molecule as given in PQS */
  if (biomol_ptr.use_pqs == TRUE)
    {
      /* Get the cut-down coordinates from the PQS macmol directory */
      ok = get_biomol_pdb(&biomol_ptr,params.pdb_code,first_chain_ptr,
                          first_residue_ptr,stored_metal,metal_count,
                          &nmetal_types,&nmetals,&first_biochain_ptr,
                          &first_bioresidue_ptr,&first_bioatom_ptr,
                          matrix,&nmatrix,&nbioatoms,file_name_ptr);

      /* Write biomol.pdb file and generate PyMol scripts to render
         molecule */
      if (ok == TRUE)
        {
          /* Set flag that have image of biomolecule */
          biomol_ptr.have_biomol_image = TRUE;

          /* Initialise variables */
          nbiossbonds = 0;
          nbioligands = 0;
          nbiometals = 0;
          first_bioligand_ptr = NULL;
          first_biossbond_ptr = NULL;
          last_pymol_ptr = NULL;
          pymol_script_ptr = NULL;

          /* Match chain-structure in biomolecule to that in
             original PDB file */
          if (first_biochain_ptr != NULL)
            mirror_molecules(first_biochain_ptr,first_bioresidue_ptr,
                             first_chain_ptr);

          /* Create all the ligand records for each of the ligand chains */
          form_ligands(first_biochain_ptr,&first_bioligand_ptr,&nbioligands,
                       &nbiometals,TRUE);

          /* Create corresponding disulphide bonds in all the symmetry
             related copies */
          nbiossbonds
            = copy_disulphides(first_ssbond_ptr,&first_biossbond_ptr,
                               first_bioresidue_ptr);

          /* Align molecule's principal axes to get better view */
          if (params.molecule_align == TRUE)
            align_molecule(first_bioatom_ptr,nbioatoms,coord_min,
                           coord_max,transformation,cofg);

          /* Copy view coordinate over original coords */
          copy_coords(first_bioatom_ptr,nbioatoms);

          /* If small number of symmetry-related units, render copies
             as transparent */
          // if (biomol_ptr.nsymop <= MAX_SYMOP)
          //   transparent_copies = TRUE;

          /* If only a single molecule, colour symmetry-related copies using
             colour of original chain */
          // if (biomol_ptr.nmolecules > 2 ||
          //     (biomol_ptr.nmolecules == 1 && transparent_copies == TRUE))
          //   colour_as_orig = TRUE;
          transparent_copies = FALSE;
          if (biomol_ptr.nmolecules > NCHAIN_COLOURS ||
              biomol_ptr.nsymop > MAX_SYMOP)
            colour_as_orig = TRUE;
          else
            colour_as_orig = FALSE;

          /* If structure not too large, generate PyMol script and
             write out the header records defining all the objects */
          if (nbioatoms <= MAX_PYMOL_ATOMS)
            {
              /* Write out biomol.pdb */
              write_pymol_pdb(params.pdb_code,first_biochain_ptr,
                              first_biossbond_ptr,nbiossbonds,
                              "biomol.pdb",TRUE);

              /* Create a pymol.pml script and write out the header records
                 defining all the objects */
              write_pymol_pml("biomol",params.pdb_code,first_biochain_ptr,
                              first_bioligand_ptr,first_biossbond_ptr,
                              nbiossbonds,stored_metal,metal_count,
                              nmetal_types,nbioligands,&pymol_script_ptr,
                              &first_pobject_ptr,colour_as_orig,
                              transparent_copies,TRUE,params.pdb_type);

              /* Set flag */
              params.gen_pymol = TRUE;
            }
          else
            params.gen_pymol = FALSE;

          /* Initialise parameters */
          nviews = 1;

          /* Generate the image of the whole biomolecule */
          transparent
            = render_structure(first_biochain_ptr,first_bioatom_ptr,
                               nbioatoms,nviews,first_biossbond_ptr,
                               nbiossbonds,nsite_res,
                               file_name_ptr.exe_name[RASTER3D],
                               file_name_ptr.exe_name[CONVERT_EXE],
                               file_name_ptr.exe_name[PPMQUANT_EXE],
                               file_name_ptr.exe_name[PPMTOGIF_EXE],mcsg,
                               &grid_sep,ntypes,nmetypes,pymol_script_ptr,
                               &last_pymol_ptr,colour_as_orig,
                               transparent_copies,TRUE,params);

          /* If using Pymol, run it to generate images */
          if (params.gen_pymol == TRUE && params.pdb_type != PDBSUM1)
            {
              /* Close the PyMol script file */
              fclose(pymol_script_ptr);

              /* Run PyMol to generate all the ray-teced .png files */
              run_pymol("biomol",file_name_ptr.exe_name[PYMOL_EXE],dos);

              /* Convert all the generated files to .gif or .jpg, as
                 required */
              convert_png_images(last_pymol_ptr,
                                 file_name_ptr.exe_name[CONVERT_EXE]);
            }
        }
    }
/* debug */
  printf("Done ...\n");
  fflush(stdout);
/* debug */

  /* Write out the biomol.dat file defining the biological unit, where
     known */
  write_biomol_dat(have_biomol,&biomol_ptr,output_dir);
/* debug */
  printf("END\n");
  fflush(stdout);
/* debug */


  return(0);
}
